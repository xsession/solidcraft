//! Design-of-experiments campaigns: factorial and Latin-hypercube case generation, template
//! rendering, case folders, a TSV case registry, solver command plans (native, Docker,
//! Apptainer, Slurm), a runtime doctor, output discovery, local process launching, log parsing
//! and finite-difference gradient descent. Ported from upstream `cfd/workflow/campaign`
//! (MIT, cfd_solvers contributors).
use super::text::{
    create_dirs, error_text, format_g17, getline_split, io_err, lower, parse_size, path_text,
    stod_full, trim_ascii, upper, write_file, ws_tokens,
};
use crate::{require, Result, SolverError};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum CampaignParameterKind {
    #[default]
    Continuous,
    Discrete,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CampaignParameter {
    pub name: String,
    pub kind: CampaignParameterKind,
    pub minimum: f64,
    pub maximum: f64,
    pub levels: Vec<String>,
}

impl CampaignParameter {
    pub fn continuous(name: impl Into<String>, minimum: f64, maximum: f64) -> Self {
        Self {
            name: name.into(),
            kind: CampaignParameterKind::Continuous,
            minimum,
            maximum,
            levels: Vec::new(),
        }
    }
    pub fn discrete<S: Into<String>>(
        name: impl Into<String>,
        levels: impl IntoIterator<Item = S>,
    ) -> Self {
        Self {
            name: name.into(),
            kind: CampaignParameterKind::Discrete,
            minimum: 0.0,
            maximum: 0.0,
            levels: levels.into_iter().map(Into::into).collect(),
        }
    }
}

/// One design row. Values are kept sorted by parameter name.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CampaignCase {
    pub case_id: String,
    pub values: BTreeMap<String, String>,
}

/// 64-bit LCG used for Latin-hypercube sampling (upstream default seed).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CampaignRandom {
    pub state: u64,
}

impl Default for CampaignRandom {
    fn default() -> Self {
        Self {
            state: 0x243f_6a88_85a3_08d3,
        }
    }
}

impl CampaignRandom {
    /// Uniform sample in the open interval (0, 1).
    pub fn uniform_open01(&mut self) -> f64 {
        self.state = self
            .state
            .wrapping_mul(2_862_933_555_777_941_757)
            .wrapping_add(3_037_000_493);
        let mantissa = (self.state >> 11) | 1;
        let v = mantissa as f64 * (1.0 / 9_007_199_254_740_992.0);
        v.clamp(f64::MIN_POSITIVE, 1.0 - f64::EPSILON / 2.0)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SolverCapability {
    Residuals,
    Probes,
    Restart,
    Compare,
    Performance,
    Control,
    Gui,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SolverAdapterDescriptor {
    pub name: String,
    pub native_binary: String,
    pub container_binary: String,
    pub capabilities: Vec<SolverCapability>,
    pub control_actions: Vec<String>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum CampaignCaseStatus {
    #[default]
    Pending,
    Running,
    Done,
    Failed,
    Stopped,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CampaignRegistryEntry {
    pub case_id: String,
    pub status: CampaignCaseStatus,
    pub objective: f64,
    pub iterations: usize,
    pub message: String,
}

impl CampaignRegistryEntry {
    pub fn new(
        case_id: impl Into<String>,
        status: CampaignCaseStatus,
        objective: f64,
        iterations: usize,
        message: impl Into<String>,
    ) -> Self {
        Self {
            case_id: case_id.into(),
            status,
            objective,
            iterations,
            message: message.into(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CampaignSummary {
    pub pending: usize,
    pub running: usize,
    pub done: usize,
    pub failed: usize,
    pub stopped: usize,
    pub best_objective: f64,
    pub best_case_id: String,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDifferenceGradientConfig {
    pub absolute_step: f64,
    pub relative_step: f64,
    pub central: bool,
}

impl Default for FiniteDifferenceGradientConfig {
    fn default() -> Self {
        Self {
            absolute_step: 1.0e-6,
            relative_step: 1.0e-4,
            central: true,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CampaignTemplateFile {
    pub relative_path: PathBuf,
    pub contents: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CampaignWriteOptions {
    pub overwrite: bool,
    pub registry_filename: String,
    pub manifest_filename: String,
    pub doe_row_filename: String,
}

impl Default for CampaignWriteOptions {
    fn default() -> Self {
        Self {
            overwrite: false,
            registry_filename: "registry.tsv".into(),
            manifest_filename: "template_manifest.tsv".into(),
            doe_row_filename: "doe_row.csv".into(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CampaignPersistenceResult {
    pub root: PathBuf,
    pub cases_written: usize,
    pub files_written: usize,
    pub registry_entries: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SolverRuntime {
    #[default]
    Native,
    Docker,
    Singularity,
    Slurm,
}

impl SolverRuntime {
    /// Upstream enum ordinal (written into TSV tables).
    pub fn ordinal(self) -> i32 {
        match self {
            Self::Native => 0,
            Self::Docker => 1,
            Self::Singularity => 2,
            Self::Slurm => 3,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SolverRunRequest {
    pub adapter: SolverAdapterDescriptor,
    pub case_directory: PathBuf,
    pub runtime: SolverRuntime,
    pub mpi_ranks: usize,
    pub threads: usize,
    pub extra_args: Vec<String>,
    pub docker_image: String,
    pub singularity_image: String,
    pub slurm_partition: String,
}

impl Default for SolverRunRequest {
    fn default() -> Self {
        Self {
            adapter: SolverAdapterDescriptor::default(),
            case_directory: PathBuf::new(),
            runtime: SolverRuntime::Native,
            mpi_ranks: 1,
            threads: 1,
            extra_args: Vec::new(),
            docker_image: String::new(),
            singularity_image: String::new(),
            slurm_partition: String::new(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SolverCommandPlan {
    pub runtime: SolverRuntime,
    pub argv: Vec<String>,
    pub environment: Vec<(String, String)>,
    pub display_command: String,
    pub requires_scheduler: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResidualSample {
    pub iteration: usize,
    pub equation: String,
    pub residual: f64,
}

impl Default for ResidualSample {
    fn default() -> Self {
        Self {
            iteration: 0,
            equation: "residual".into(),
            residual: 0.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PerformanceMetric {
    pub key: String,
    pub value: f64,
    pub units: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RuntimeDoctorCheck {
    pub name: String,
    pub ok: bool,
    pub message: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RuntimeDoctorReport {
    pub checks: Vec<RuntimeDoctorCheck>,
    pub ok: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CampaignOutputDiscovery {
    pub log_files: Vec<PathBuf>,
    pub residual_files: Vec<PathBuf>,
    pub performance_files: Vec<PathBuf>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalSolverRunOptions {
    pub write_launcher_logs: bool,
    pub stdout_filename: String,
    pub stderr_filename: String,
    pub max_capture_bytes: usize,
}

impl Default for LocalSolverRunOptions {
    fn default() -> Self {
        Self {
            write_launcher_logs: true,
            stdout_filename: "solver.stdout.log".into(),
            stderr_filename: "solver.stderr.log".into(),
            max_capture_bytes: 1 << 20,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct LocalSolverRunResult {
    pub case_id: String,
    pub status: CampaignCaseStatus,
    pub exit_code: i32,
    pub stdout_text: String,
    pub stderr_text: String,
    pub residuals: Vec<ResidualSample>,
    pub performance: Vec<PerformanceMetric>,
    pub discovered_outputs: CampaignOutputDiscovery,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalCampaignRunOptions {
    pub adapter: SolverAdapterDescriptor,
    pub runtime: SolverRuntime,
    pub mpi_ranks: usize,
    pub threads: usize,
    pub extra_args: Vec<String>,
    pub pending_only: bool,
    pub solver_run: LocalSolverRunOptions,
}

impl Default for LocalCampaignRunOptions {
    fn default() -> Self {
        Self {
            adapter: SolverAdapterDescriptor::default(),
            runtime: SolverRuntime::Native,
            mpi_ranks: 1,
            threads: 1,
            extra_args: Vec::new(),
            pending_only: true,
            solver_run: LocalSolverRunOptions::default(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct LocalCampaignRunResult {
    pub case_results: Vec<LocalSolverRunResult>,
    pub summary: CampaignSummary,
    pub launched: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CampaignOptimizationConfig {
    pub max_iterations: usize,
    pub step_size: f64,
    pub tolerance: f64,
    pub finite_difference: FiniteDifferenceGradientConfig,
}

impl Default for CampaignOptimizationConfig {
    fn default() -> Self {
        Self {
            max_iterations: 32,
            step_size: 0.1,
            tolerance: 1.0e-6,
            finite_difference: FiniteDifferenceGradientConfig::default(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CampaignOptimizationResult {
    pub parameters: Vec<f64>,
    pub objective: f64,
    pub iterations: usize,
    pub converged: bool,
    pub last_gradient: Vec<f64>,
}

fn case_id(index: usize) -> String {
    format!("case{:04}", index + 1)
}

fn truthy(v: &str) -> bool {
    !(v.is_empty() || v == "false" || v == "False" || v == "FALSE" || v == "0")
}

fn validate_parameter(p: &CampaignParameter) -> Result<()> {
    require(!p.name.is_empty(), "campaign parameter requires a name")?;
    require(p.name != "case_id", "case_id is reserved")?;
    match p.kind {
        CampaignParameterKind::Continuous => require(
            p.maximum >= p.minimum && p.minimum.is_finite() && p.maximum.is_finite(),
            "invalid continuous campaign parameter",
        ),
        CampaignParameterKind::Discrete => require(
            !p.levels.is_empty(),
            "discrete campaign parameter requires levels",
        ),
    }
}

/// Full factorial over discrete parameters; the first parameter varies slowest.
pub fn generate_factorial_campaign(parameters: &[CampaignParameter]) -> Result<Vec<CampaignCase>> {
    require(
        !parameters.is_empty(),
        "factorial campaign requires parameters",
    )?;
    for p in parameters {
        validate_parameter(p)?;
        require(
            p.kind == CampaignParameterKind::Discrete,
            "factorial campaign only supports discrete parameters",
        )?;
    }
    let mut cases = vec![CampaignCase::default()];
    for p in parameters {
        let mut next = Vec::with_capacity(cases.len() * p.levels.len());
        for base in &cases {
            for level in &p.levels {
                let mut row = base.clone();
                row.values.insert(p.name.clone(), level.clone());
                next.push(row);
            }
        }
        cases = next;
    }
    for (i, c) in cases.iter_mut().enumerate() {
        c.case_id = case_id(i);
    }
    Ok(cases)
}

/// Latin hypercube: one jittered stratum per sample, Fisher–Yates shuffled per continuous
/// parameter (values written with 12 significant digits); discrete parameters cycle levels.
pub fn generate_latin_hypercube_campaign(
    parameters: &[CampaignParameter],
    samples: usize,
    mut rng: CampaignRandom,
) -> Result<Vec<CampaignCase>> {
    require(
        !parameters.is_empty() && samples > 0,
        "latin hypercube campaign requires parameters and samples",
    )?;
    for p in parameters {
        validate_parameter(p)?;
    }
    let mut cases: Vec<CampaignCase> = (0..samples)
        .map(|i| CampaignCase {
            case_id: case_id(i),
            values: BTreeMap::new(),
        })
        .collect();
    for p in parameters {
        match p.kind {
            CampaignParameterKind::Continuous => {
                let mut strata: Vec<f64> = (0..samples)
                    .map(|i| (i as f64 + rng.uniform_open01()) / samples as f64)
                    .collect();
                for i in (2..=samples).rev() {
                    let j = ((rng.uniform_open01() * i as f64) as usize).min(i - 1);
                    strata.swap(i - 1, j);
                }
                for (c, s) in cases.iter_mut().zip(&strata) {
                    let v = p.minimum + (p.maximum - p.minimum) * s;
                    c.values
                        .insert(p.name.clone(), crate::lbm::format_general(v, 12));
                }
            }
            CampaignParameterKind::Discrete => {
                for (i, c) in cases.iter_mut().enumerate() {
                    c.values
                        .insert(p.name.clone(), p.levels[i % p.levels.len()].clone());
                }
            }
        }
    }
    Ok(cases)
}

fn replace_placeholders(line: &str, row: &CampaignCase, strict: bool) -> Result<String> {
    let mut line = line.to_string();
    let mut pos = 0;
    while let Some(off) = line.get(pos..).and_then(|s| s.find('{')) {
        let start = pos + off;
        let Some(end_off) = line.get(start + 1..).and_then(|s| s.find('}')) else {
            break;
        };
        let end = start + 1 + end_off;
        let key = line[start + 1..end].to_string();
        match row.values.get(&key) {
            Some(value) => {
                line.replace_range(start..=end, value);
                pos = start + value.len();
            }
            None => {
                require(
                    !strict,
                    format!("missing campaign template placeholder: {key}"),
                )?;
                pos = end + 1;
            }
        }
    }
    Ok(line)
}

/// Render `{name}` placeholders and `<!-- IF key -->`, `<!-- IF key=value -->`,
/// `<!-- IF key!=value -->` … `<!-- ENDIF -->` line blocks. Every emitted line ends with `\n`.
pub fn render_campaign_template(templ: &str, row: &CampaignCase, strict: bool) -> Result<String> {
    let mut out = String::new();
    let mut stack: Vec<bool> = Vec::new();
    for line in getline_split(templ, b'\n') {
        if let Some(if_pos) = line.find("<!-- IF ") {
            let Some(close_off) = line[if_pos..].find("-->") else {
                return Err(SolverError::InvalidInput(
                    "unterminated campaign IF block".into(),
                ));
            };
            let close = if_pos + close_off;
            let expr = trim_ascii(line.get(if_pos + 8..close).unwrap_or(""));
            let include = if let Some(neq) = expr.find("!=") {
                let key = trim_ascii(&expr[..neq]);
                let expected = trim_ascii(&expr[neq + 2..]);
                row.values.get(key).is_some_and(|v| v != expected)
            } else if let Some(eq) = expr.find('=') {
                let key = trim_ascii(&expr[..eq]);
                let expected = trim_ascii(&expr[eq + 1..]);
                row.values.get(key).is_some_and(|v| v == expected)
            } else {
                row.values.get(expr).is_some_and(|v| truthy(v))
            };
            stack.push(include);
            continue;
        }
        if line.contains("<!-- ENDIF -->") {
            require(stack.pop().is_some(), "campaign ENDIF without IF")?;
            continue;
        }
        if stack.iter().all(|&v| v) {
            out.push_str(&replace_placeholders(line, row, strict)?);
            out.push('\n');
        }
    }
    require(stack.is_empty(), "unclosed campaign IF block")?;
    Ok(out)
}

pub fn validate_solver_adapter_descriptor(adapter: &SolverAdapterDescriptor) -> Result<()> {
    require(
        !adapter.name.is_empty()
            && !adapter.native_binary.is_empty()
            && !adapter.container_binary.is_empty(),
        "solver adapter descriptor requires identity and binaries",
    )?;
    require(
        adapter.control_actions.iter().all(|a| !a.is_empty()),
        "solver control action must not be empty",
    )
}

pub fn solver_has_capability(
    adapter: &SolverAdapterDescriptor,
    capability: SolverCapability,
) -> bool {
    adapter.capabilities.contains(&capability)
}

/// Insert or replace the entry with the same case id.
pub fn update_case_status(
    registry: &mut Vec<CampaignRegistryEntry>,
    entry: CampaignRegistryEntry,
) -> Result<()> {
    require(
        !entry.case_id.is_empty(),
        "campaign registry entry needs case_id",
    )?;
    match registry.iter_mut().find(|e| e.case_id == entry.case_id) {
        Some(slot) => *slot = entry,
        None => registry.push(entry),
    }
    Ok(())
}

/// Status counts plus the done case with the lowest objective (first wins ties).
pub fn summarize_campaign_registry(registry: &[CampaignRegistryEntry]) -> CampaignSummary {
    let mut out = CampaignSummary::default();
    let mut have_best = false;
    for item in registry {
        match item.status {
            CampaignCaseStatus::Pending => out.pending += 1,
            CampaignCaseStatus::Running => out.running += 1,
            CampaignCaseStatus::Done => out.done += 1,
            CampaignCaseStatus::Failed => out.failed += 1,
            CampaignCaseStatus::Stopped => out.stopped += 1,
        }
        if item.status == CampaignCaseStatus::Done
            && (!have_best || item.objective < out.best_objective)
        {
            have_best = true;
            out.best_objective = item.objective;
            out.best_case_id = item.case_id.clone();
        }
    }
    out
}

fn escape_tsv(v: &str) -> String {
    let mut o = String::with_capacity(v.len());
    for c in v.chars() {
        match c {
            '\\' => o.push_str("\\\\"),
            '\t' => o.push_str("\\t"),
            '\n' => o.push_str("\\n"),
            _ => o.push(c),
        }
    }
    o
}

fn unescape_tsv(v: &str) -> String {
    let mut o = String::with_capacity(v.len());
    let mut it = v.chars().peekable();
    while let Some(c) = it.next() {
        if c == '\\' {
            if let Some(n) = it.next() {
                o.push(match n {
                    't' => '\t',
                    'n' => '\n',
                    other => other,
                });
                continue;
            }
        }
        o.push(c);
    }
    o
}

fn ensure_safe_relative_path(relative: &Path) -> Result<()> {
    require(
        !relative.as_os_str().is_empty() && !relative.is_absolute() && !relative.has_root(),
        "campaign template path must be relative",
    )?;
    require(
        !relative
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir)),
        "campaign template path must not escape case directory",
    )
}

/// Campaign-style command display: quote an argument with single quotes only when it contains
/// whitespace or a quote.
fn join_display(argv: &[String]) -> String {
    let mut out = String::new();
    for (i, a) in argv.iter().enumerate() {
        if i > 0 {
            out.push(' ');
        }
        if a.contains([' ', '\t', '\n', '"', '\'']) {
            out.push('\'');
            for c in a.chars() {
                if c == '\'' {
                    out.push_str("'\\''");
                } else {
                    out.push(c);
                }
            }
            out.push('\'');
        } else {
            out.push_str(a);
        }
    }
    out
}

/// Finite `std::stod` of a whole token.
fn parse_double(token: &str) -> Option<f64> {
    stod_full(token).filter(|v| v.is_finite())
}

fn normalized_tokens(line: &str) -> Vec<String> {
    let mapped: String = line
        .chars()
        .map(|c| {
            if matches!(c, '=' | ',' | ';' | ':') {
                ' '
            } else {
                c
            }
        })
        .collect();
    ws_tokens(&mapped)
}

fn is_executable_file(path: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::metadata(path).is_ok_and(|m| m.permissions().mode() & 0o111 != 0)
    }
    #[cfg(not(unix))]
    {
        path.is_file()
    }
}

pub(crate) fn find_executable_on_path(name: &str) -> Option<PathBuf> {
    if name.is_empty() {
        return None;
    }
    if name.contains('/') || name.contains('\\') {
        let p = PathBuf::from(name);
        return is_executable_file(&p).then_some(p);
    }
    let env = std::env::var("PATH").ok()?;
    getline_split(&env, b':').into_iter().find_map(|dir| {
        let dir = if dir.is_empty() { "." } else { dir };
        let candidate = Path::new(dir).join(name);
        is_executable_file(&candidate).then_some(candidate)
    })
}

fn read_text_limited(path: &Path, max_bytes: usize) -> String {
    use std::io::Read;
    let Ok(file) = std::fs::File::open(path) else {
        return String::new();
    };
    let mut buf = Vec::new();
    let limit = u64::try_from(max_bytes).unwrap_or(u64::MAX);
    if file.take(limit).read_to_end(&mut buf).is_err() {
        return String::new();
    }
    String::from_utf8_lossy(&buf).into_owned()
}

pub(crate) fn read_text_full(path: &Path) -> String {
    std::fs::read(path)
        .map(|b| String::from_utf8_lossy(&b).into_owned())
        .unwrap_or_default()
}

fn objective_from_residuals(residuals: &[ResidualSample], exit_code: i32) -> f64 {
    match residuals.last() {
        Some(r) => r.residual,
        None if exit_code == 0 => 0.0,
        None => f64::from(exit_code.max(1)),
    }
}

pub(crate) fn max_iteration(residuals: &[ResidualSample]) -> usize {
    residuals.iter().map(|s| s.iteration).max().unwrap_or(0)
}

pub(crate) fn sort_residuals(residuals: &mut [ResidualSample]) {
    residuals.sort_by(|a, b| {
        a.iteration
            .cmp(&b.iteration)
            .then_with(|| a.equation.cmp(&b.equation))
            .then_with(|| a.residual.total_cmp(&b.residual))
    });
}

pub fn campaign_status_name(status: CampaignCaseStatus) -> &'static str {
    match status {
        CampaignCaseStatus::Pending => "pending",
        CampaignCaseStatus::Running => "running",
        CampaignCaseStatus::Done => "done",
        CampaignCaseStatus::Failed => "failed",
        CampaignCaseStatus::Stopped => "stopped",
    }
}

/// Case-insensitive; accepts the aliases SUCCESS/COMPLETED, ERROR and CANCELLED.
pub fn campaign_status_from_name(name: &str) -> Result<CampaignCaseStatus> {
    Ok(match upper(name).as_str() {
        "PENDING" => CampaignCaseStatus::Pending,
        "RUNNING" => CampaignCaseStatus::Running,
        "DONE" | "SUCCESS" | "COMPLETED" => CampaignCaseStatus::Done,
        "FAILED" | "ERROR" => CampaignCaseStatus::Failed,
        "STOPPED" | "CANCELLED" => CampaignCaseStatus::Stopped,
        _ => {
            return Err(SolverError::InvalidInput(
                "unknown campaign case status".into(),
            ))
        }
    })
}

/// Text of the TSV registry (`case_id status objective iterations message`, objective `%.17g`).
pub fn campaign_registry_tsv(registry: &[CampaignRegistryEntry]) -> String {
    let mut out = String::from("case_id\tstatus\tobjective\titerations\tmessage\n");
    for e in registry {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            escape_tsv(&e.case_id),
            campaign_status_name(e.status),
            format_g17(e.objective),
            e.iterations,
            escape_tsv(&e.message)
        ));
    }
    out
}

pub fn save_campaign_registry(path: &Path, registry: &[CampaignRegistryEntry]) -> Result<()> {
    require(
        !path.as_os_str().is_empty(),
        "campaign registry path is empty",
    )?;
    if let Some(parent) = path.parent() {
        create_dirs(parent)?;
    }
    write_file(
        path,
        &campaign_registry_tsv(registry),
        "failed to write campaign registry",
    )
}

/// Parse registry text. Unlike upstream, an infinite objective (written by the local campaign
/// runner for cases that could not be launched) is accepted so such a registry can be reloaded.
pub fn parse_campaign_registry(text: &str) -> Result<Vec<CampaignRegistryEntry>> {
    let mut registry = Vec::new();
    let mut header = true;
    for line in getline_split(text, b'\n') {
        if line.is_empty() {
            continue;
        }
        if header {
            header = false;
            if line.starts_with("case_id\t") {
                continue;
            }
        }
        let fields: Vec<&str> = line.split('\t').collect();
        require(fields.len() >= 5, "malformed campaign registry row")?;
        let status = campaign_status_from_name(fields[1])?;
        let objective = stod_full(fields[2]).filter(|v| !v.is_nan());
        let iterations = parse_size(fields[3]);
        let (Some(objective), Some(iterations)) = (objective, iterations) else {
            return Err(SolverError::InvalidInput(
                "malformed numeric campaign registry field".into(),
            ));
        };
        registry.push(CampaignRegistryEntry {
            case_id: unescape_tsv(fields[0]),
            status,
            objective,
            iterations,
            message: unescape_tsv(fields[4]),
        });
    }
    Ok(registry)
}

pub fn load_campaign_registry(path: &Path) -> Result<Vec<CampaignRegistryEntry>> {
    let text = std::fs::read(path).map_err(io_err("failed to read campaign registry"))?;
    parse_campaign_registry(&String::from_utf8_lossy(&text))
}

/// Write one folder per case with every rendered template, a `doe_row.csv`, the registry
/// (all cases `pending`, message `prepared`) and the template manifest.
pub fn write_campaign_case_folders(
    root: &Path,
    cases: &[CampaignCase],
    templates: &[CampaignTemplateFile],
    options: &CampaignWriteOptions,
) -> Result<CampaignPersistenceResult> {
    require(
        !root.as_os_str().is_empty() && !cases.is_empty() && !templates.is_empty(),
        "campaign folder writer requires root, cases and templates",
    )?;
    create_dirs(root)?;
    for t in templates {
        ensure_safe_relative_path(&t.relative_path)?;
    }
    let mut registry = Vec::with_capacity(cases.len());
    let mut result = CampaignPersistenceResult {
        root: root.to_path_buf(),
        ..Default::default()
    };
    for row in cases {
        require(!row.case_id.is_empty(), "campaign case requires case_id")?;
        let case_dir = root.join(&row.case_id);
        if case_dir.exists() && !options.overwrite {
            return Err(SolverError::Io(format!(
                "campaign case directory already exists: {}",
                case_dir.display()
            )));
        }
        create_dirs(&case_dir)?;
        for t in templates {
            let target = case_dir.join(&t.relative_path);
            if let Some(parent) = target.parent() {
                create_dirs(parent)?;
            }
            let text = render_campaign_template(&t.contents, row, true)?;
            write_file(&target, &text, "failed to write campaign template file")?;
            result.files_written += 1;
        }
        let mut doe = String::from("case_id");
        for k in row.values.keys() {
            doe.push(',');
            doe.push_str(k);
        }
        doe.push('\n');
        doe.push_str(&row.case_id);
        for v in row.values.values() {
            doe.push(',');
            doe.push_str(v);
        }
        doe.push('\n');
        write_file(
            &case_dir.join(&options.doe_row_filename),
            &doe,
            "failed to write campaign doe row",
        )?;
        result.files_written += 1;
        registry.push(CampaignRegistryEntry::new(
            row.case_id.clone(),
            CampaignCaseStatus::Pending,
            0.0,
            0,
            "prepared",
        ));
        result.cases_written += 1;
    }
    save_campaign_registry(&root.join(&options.registry_filename), &registry)?;
    result.registry_entries = registry.len();
    let mut manifest = String::from("relative_path\n");
    for t in templates {
        manifest.push_str(&path_text(&t.relative_path));
        manifest.push('\n');
    }
    write_file(
        &root.join(&options.manifest_filename),
        &manifest,
        "failed to write campaign manifest",
    )?;
    Ok(result)
}

/// argv, environment and display string for a native, Docker, Apptainer or Slurm launch.
pub fn build_solver_command_plan(request: &SolverRunRequest) -> Result<SolverCommandPlan> {
    validate_solver_adapter_descriptor(&request.adapter)?;
    require(
        !request.case_directory.as_os_str().is_empty()
            && request.mpi_ranks > 0
            && request.threads > 0,
        "invalid solver run request",
    )?;
    let mut plan = SolverCommandPlan {
        runtime: request.runtime,
        environment: vec![
            ("CFD_SOLVER_ADAPTER".into(), request.adapter.name.clone()),
            ("OMP_NUM_THREADS".into(), request.threads.to_string()),
        ],
        ..Default::default()
    };
    let case_path = request.case_directory.to_string_lossy().into_owned();
    let append_common = |argv: &mut Vec<String>, case_arg: &str| {
        argv.extend([
            "--case".to_string(),
            case_arg.to_string(),
            "--nprocs".to_string(),
            request.mpi_ranks.to_string(),
            "--threads".to_string(),
            request.threads.to_string(),
        ]);
        argv.extend(request.extra_args.iter().cloned());
    };
    let s = |v: &str| v.to_string();
    match request.runtime {
        SolverRuntime::Native => {
            plan.argv.push(request.adapter.native_binary.clone());
            append_common(&mut plan.argv, &case_path);
        }
        SolverRuntime::Docker => {
            require(
                !request.docker_image.is_empty(),
                "docker solver run requires image",
            )?;
            plan.argv = vec![
                s("docker"),
                s("run"),
                s("--rm"),
                s("-v"),
                format!("{case_path}:/case"),
                s("-w"),
                s("/case"),
                request.docker_image.clone(),
                request.adapter.container_binary.clone(),
            ];
            append_common(&mut plan.argv, "/case");
        }
        SolverRuntime::Singularity => {
            require(
                !request.singularity_image.is_empty(),
                "singularity solver run requires image",
            )?;
            plan.argv = vec![
                s("apptainer"),
                s("exec"),
                s("--bind"),
                format!("{case_path}:/case"),
                request.singularity_image.clone(),
                request.adapter.container_binary.clone(),
            ];
            append_common(&mut plan.argv, "/case");
        }
        SolverRuntime::Slurm => {
            plan.requires_scheduler = true;
            let mut inner = vec![request.adapter.native_binary.clone()];
            append_common(&mut inner, &case_path);
            plan.argv = vec![
                s("sbatch"),
                s("--ntasks"),
                request.mpi_ranks.to_string(),
                s("--cpus-per-task"),
                request.threads.to_string(),
            ];
            if !request.slurm_partition.is_empty() {
                plan.argv.push(s("--partition"));
                plan.argv.push(request.slurm_partition.clone());
            }
            plan.argv.push(s("--wrap"));
            plan.argv.push(join_display(&inner));
        }
    }
    plan.display_command = join_display(&plan.argv);
    Ok(plan)
}

/// Check adapter, case directory, launcher on `PATH` and container/scheduler configuration.
pub fn doctor_solver_runtime(request: &SolverRunRequest) -> RuntimeDoctorReport {
    let mut report = RuntimeDoctorReport::default();
    let mut add = |name: &str, ok: bool, message: String| {
        report.checks.push(RuntimeDoctorCheck {
            name: name.into(),
            ok,
            message,
        })
    };
    if let Err(e) = validate_solver_adapter_descriptor(&request.adapter) {
        add("adapter", false, error_text(&e));
        report.ok = false;
        return report;
    }
    add("adapter", true, "solver adapter descriptor is valid".into());
    let case_ok = !request.case_directory.as_os_str().is_empty() && request.case_directory.is_dir();
    add(
        "case_directory",
        case_ok,
        if case_ok {
            "case directory exists".into()
        } else {
            "case directory is missing or not a directory".into()
        },
    );
    match build_solver_command_plan(request) {
        Ok(plan) => {
            let launcher = plan.argv.first().cloned().unwrap_or_default();
            match find_executable_on_path(&launcher) {
                Some(found) => add(
                    "launcher",
                    true,
                    format!("found executable: {}", found.display()),
                ),
                None => add(
                    "launcher",
                    false,
                    format!("executable not found: {launcher}"),
                ),
            }
            match request.runtime {
                SolverRuntime::Docker => {
                    let ok = !request.docker_image.is_empty();
                    add(
                        "docker_image",
                        ok,
                        if ok {
                            "docker image configured"
                        } else {
                            "docker image not configured"
                        }
                        .into(),
                    );
                }
                SolverRuntime::Singularity => {
                    let ok = !request.singularity_image.is_empty();
                    add(
                        "singularity_image",
                        ok,
                        if ok {
                            "singularity image configured"
                        } else {
                            "singularity image not configured"
                        }
                        .into(),
                    );
                }
                SolverRuntime::Slurm => add(
                    "scheduler",
                    plan.requires_scheduler,
                    "slurm command plan requires scheduler submission".into(),
                ),
                SolverRuntime::Native => {}
            }
        }
        Err(e) => add("command_plan", false, error_text(&e)),
    }
    report.ok = report.checks.iter().all(|c| c.ok);
    report
}

fn walk_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        match entry.file_type() {
            Ok(ft) if ft.is_dir() => walk_files(&path, out),
            Ok(_) => {
                if path.is_file() {
                    out.push(path);
                }
            }
            Err(_) => {}
        }
    }
}

/// Recursively classify `.log/.out/.txt/.dat/.csv` files into logs (`.log`, `.out`, or names
/// containing stdout/stderr), residual histories (`residual`, `history`) and performance files
/// (`perf`, `timing`, `time`). Each list is sorted.
pub fn discover_campaign_outputs(case_directory: &Path) -> CampaignOutputDiscovery {
    let mut out = CampaignOutputDiscovery::default();
    if case_directory.as_os_str().is_empty() || !case_directory.exists() {
        return out;
    }
    let mut files = Vec::new();
    walk_files(case_directory, &mut files);
    for path in files {
        let name = lower(
            &path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
        );
        let ext = path
            .extension()
            .map(|e| format!(".{}", lower(&e.to_string_lossy())))
            .unwrap_or_default();
        if !matches!(ext.as_str(), ".log" | ".out" | ".txt" | ".dat" | ".csv") {
            continue;
        }
        if ext == ".log" || ext == ".out" || name.contains("stdout") || name.contains("stderr") {
            push_unique(&mut out.log_files, &path);
        }
        if name.contains("residual") || name.contains("history") {
            push_unique(&mut out.residual_files, &path);
        }
        if name.contains("perf") || name.contains("timing") || name.contains("time") {
            push_unique(&mut out.performance_files, &path);
        }
    }
    out.log_files.sort();
    out.residual_files.sort();
    out.performance_files.sort();
    out
}

fn push_unique(paths: &mut Vec<PathBuf>, path: &Path) {
    if !paths.iter().any(|p| p == path) {
        paths.push(path.to_path_buf());
    }
}

/// Launch one case locally (native, Docker or Apptainer; Slurm plans are rejected) after a
/// passing runtime doctor, capture stdout/stderr into the case folder, then parse residuals and
/// performance metrics from the logs and discovered files.
pub fn run_local_solver_case(
    request: &SolverRunRequest,
    options: &LocalSolverRunOptions,
) -> Result<LocalSolverRunResult> {
    let plan = build_solver_command_plan(request)?;
    require(
        !plan.requires_scheduler,
        "local campaign runner does not submit scheduler jobs",
    )?;
    let doctor = doctor_solver_runtime(request);
    if !doctor.ok {
        let mut msg = String::from("runtime doctor failed");
        for c in doctor.checks.iter().filter(|c| !c.ok) {
            msg.push_str(&format!("; {}: {}", c.name, c.message));
        }
        return Err(SolverError::NumericalFailure(msg));
    }
    let mut result = LocalSolverRunResult {
        case_id: request
            .case_directory
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default(),
        ..Default::default()
    };
    let stdout_path = request.case_directory.join(&options.stdout_filename);
    let stderr_path = request.case_directory.join(&options.stderr_filename);
    let out_file =
        std::fs::File::create(&stdout_path).map_err(io_err("failed to open solver stdout log"))?;
    let err_file =
        std::fs::File::create(&stderr_path).map_err(io_err("failed to open solver stderr log"))?;
    let (program, args) = plan
        .argv
        .split_first()
        .ok_or_else(|| SolverError::InvalidInput("empty solver command".into()))?;
    let spawned = std::process::Command::new(program)
        .args(args)
        .current_dir(&request.case_directory)
        .envs(
            plan.environment
                .iter()
                .map(|(k, v)| (k.as_str(), v.as_str())),
        )
        .stdin(std::process::Stdio::null())
        .stdout(out_file)
        .stderr(err_file)
        .status();
    result.exit_code = match spawned {
        Ok(status) => exit_code_of(status),
        Err(_) => {
            write_file(
                &stderr_path,
                &format!("execvp failed for {program}\n"),
                "failed to write solver stderr log",
            )?;
            127
        }
    };
    result.stdout_text = read_text_limited(&stdout_path, options.max_capture_bytes);
    result.stderr_text = read_text_limited(&stderr_path, options.max_capture_bytes);
    result.discovered_outputs = discover_campaign_outputs(&request.case_directory);
    let combined = format!("{}\n{}", result.stdout_text, result.stderr_text);
    result.residuals = parse_residual_history(&combined);
    result.performance = parse_performance_metrics(&combined);
    for f in &result.discovered_outputs.residual_files {
        result
            .residuals
            .extend(parse_residual_history(&read_text_full(f)));
    }
    for f in &result.discovered_outputs.performance_files {
        result
            .performance
            .extend(parse_performance_metrics(&read_text_full(f)));
    }
    sort_residuals(&mut result.residuals);
    let outcome = detect_solver_outcome_from_log(&combined);
    result.status = if result.exit_code == 0 {
        if outcome == CampaignCaseStatus::Running {
            CampaignCaseStatus::Done
        } else {
            outcome
        }
    } else if outcome == CampaignCaseStatus::Stopped {
        CampaignCaseStatus::Stopped
    } else {
        CampaignCaseStatus::Failed
    };
    if !options.write_launcher_logs {
        std::fs::remove_file(&stdout_path).map_err(io_err("cannot remove solver stdout log"))?;
        std::fs::remove_file(&stderr_path).map_err(io_err("cannot remove solver stderr log"))?;
    }
    Ok(result)
}

fn exit_code_of(status: std::process::ExitStatus) -> i32 {
    if let Some(code) = status.code() {
        return code;
    }
    #[cfg(unix)]
    {
        use std::os::unix::process::ExitStatusExt;
        if let Some(sig) = status.signal() {
            return 128 + sig;
        }
    }
    127
}

/// Run every (pending) registry case locally, saving the registry before and after each launch.
pub fn run_local_campaign(
    root: &Path,
    options: &LocalCampaignRunOptions,
    registry_filename: &str,
) -> Result<LocalCampaignRunResult> {
    require(
        !root.as_os_str().is_empty() && !registry_filename.is_empty(),
        "local campaign runner requires a root and registry filename",
    )?;
    validate_solver_adapter_descriptor(&options.adapter)?;
    let registry_path = root.join(registry_filename);
    let mut registry = load_campaign_registry(&registry_path)?;
    let mut result = LocalCampaignRunResult::default();
    for i in 0..registry.len() {
        if options.pending_only && registry[i].status != CampaignCaseStatus::Pending {
            continue;
        }
        let request = SolverRunRequest {
            adapter: options.adapter.clone(),
            runtime: options.runtime,
            mpi_ranks: options.mpi_ranks,
            threads: options.threads,
            extra_args: options.extra_args.clone(),
            case_directory: root.join(&registry[i].case_id),
            ..Default::default()
        };
        registry[i].status = CampaignCaseStatus::Running;
        registry[i].message = "launched".into();
        save_campaign_registry(&registry_path, &registry)?;
        match run_local_solver_case(&request, &options.solver_run) {
            Ok(run) => {
                result.launched += 1;
                let e = &mut registry[i];
                e.status = run.status;
                e.objective = objective_from_residuals(&run.residuals, run.exit_code);
                e.iterations = max_iteration(&run.residuals);
                e.message = format!(
                    "exit={} residuals={} perf={}",
                    run.exit_code,
                    run.residuals.len(),
                    run.performance.len()
                );
                result.case_results.push(run);
            }
            Err(err) => {
                let e = &mut registry[i];
                e.status = CampaignCaseStatus::Failed;
                e.objective = f64::INFINITY;
                e.message = error_text(&err);
            }
        }
        save_campaign_registry(&registry_path, &registry)?;
    }
    result.summary = summarize_campaign_registry(&registry);
    Ok(result)
}

/// FAILED/FATAL/ABORT/ERROR → failed; STOPPED/CANCELLED → stopped; DONE/SUCCESS/COMPLETED →
/// done; otherwise running (case-insensitive substring search, in that priority order).
pub fn detect_solver_outcome_from_log(log: &str) -> CampaignCaseStatus {
    let t = upper(log);
    let any = |keys: &[&str]| keys.iter().any(|k| t.contains(k));
    if any(&["FAILED", "FATAL", "ABORT", "ERROR"]) {
        CampaignCaseStatus::Failed
    } else if any(&["STOPPED", "CANCELLED"]) {
        CampaignCaseStatus::Stopped
    } else if any(&["DONE", "SUCCESS", "COMPLETED"]) {
        CampaignCaseStatus::Done
    } else {
        CampaignCaseStatus::Running
    }
}

/// Residual lines are either `iteration residual …` or contain `iter[ation] N` and
/// `res[idual] R` pairs (optionally `equation|field NAME`); `= , ; :` act as separators.
pub fn parse_residual_history(log: &str) -> Vec<ResidualSample> {
    let mut out = Vec::new();
    for line in getline_split(log, b'\n') {
        let tokens = normalized_tokens(line);
        if tokens.is_empty() {
            continue;
        }
        if tokens.len() >= 2 {
            if let (Some(it), Some(r)) = (parse_size(&tokens[0]), parse_double(&tokens[1])) {
                out.push(ResidualSample {
                    iteration: it,
                    equation: "residual".into(),
                    residual: r,
                });
                continue;
            }
        }
        let mut sample = ResidualSample::default();
        let (mut have_iter, mut have_res) = (false, false);
        for i in 0..tokens.len() - 1 {
            let key = upper(&tokens[i]);
            let next = &tokens[i + 1];
            if key == "ITER" || key == "ITERATION" {
                if let Some(v) = parse_size(next) {
                    sample.iteration = v;
                    have_iter = true;
                    continue;
                }
            }
            if key == "RES" || key == "RESIDUAL" {
                if let Some(v) = parse_double(next) {
                    sample.residual = v;
                    have_res = true;
                }
            } else if key == "EQUATION" || key == "FIELD" {
                sample.equation = next.clone();
            }
        }
        if have_iter && have_res {
            out.push(sample);
        }
    }
    out
}

/// `perf.KEY value [units]` (or `PERF.`), and the bare keys `wall_time`, `cpu_time`,
/// `cells_per_s`, `iterations_per_s` (case-insensitive).
pub fn parse_performance_metrics(log: &str) -> Vec<PerformanceMetric> {
    let mut out = Vec::new();
    for line in getline_split(log, b'\n') {
        let tokens = normalized_tokens(line);
        if tokens.len() < 2 {
            continue;
        }
        let key = if let Some(k) = tokens[0]
            .strip_prefix("perf.")
            .or_else(|| tokens[0].strip_prefix("PERF."))
        {
            k.to_string()
        } else if matches!(
            upper(&tokens[0]).as_str(),
            "WALL_TIME" | "CPU_TIME" | "CELLS_PER_S" | "ITERATIONS_PER_S"
        ) {
            tokens[0].clone()
        } else {
            continue;
        };
        let Some(value) = parse_double(&tokens[1]) else {
            continue;
        };
        out.push(PerformanceMetric {
            key,
            value,
            units: tokens.get(2).cloned().unwrap_or_default(),
        });
    }
    out
}

/// Steepest descent with finite-difference gradients and step halving (12 backtracks).
pub fn run_gradient_descent_campaign(
    objective: &dyn Fn(&[f64]) -> f64,
    initial: Vec<f64>,
    config: &CampaignOptimizationConfig,
) -> Result<CampaignOptimizationResult> {
    require(
        !initial.is_empty()
            && config.max_iterations > 0
            && config.step_size > 0.0
            && config.tolerance > 0.0,
        "invalid campaign optimization request",
    )?;
    let mut result = CampaignOptimizationResult {
        objective: objective(&initial),
        parameters: initial,
        ..Default::default()
    };
    let mut step = config.step_size;
    for iter in 0..config.max_iterations {
        result.last_gradient =
            finite_difference_gradient(objective, &result.parameters, &config.finite_difference)?;
        let norm = result
            .last_gradient
            .iter()
            .map(|g| g * g)
            .sum::<f64>()
            .sqrt();
        result.iterations = iter + 1;
        if norm < config.tolerance {
            result.converged = true;
            break;
        }
        let mut accepted = false;
        for _ in 0..12 {
            let trial = gradient_descent_update(&result.parameters, &result.last_gradient, step)?;
            let value = objective(&trial);
            if value < result.objective {
                result.parameters = trial;
                result.objective = value;
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        if !accepted || step < config.tolerance {
            result.converged = true;
            break;
        }
    }
    Ok(result)
}

/// Central (default) or forward differences with `h = max(absolute, relative·|x_i|)`.
pub fn finite_difference_gradient(
    objective: &dyn Fn(&[f64]) -> f64,
    parameters: &[f64],
    config: &FiniteDifferenceGradientConfig,
) -> Result<Vec<f64>> {
    require(
        !parameters.is_empty() && config.absolute_step > 0.0 && config.relative_step >= 0.0,
        "invalid finite difference gradient request",
    )?;
    let x = parameters.to_vec();
    let base = objective(&x);
    let mut gradient = vec![0.0; x.len()];
    for i in 0..x.len() {
        let rel = config.relative_step * x[i].abs();
        let h = if config.absolute_step < rel {
            rel
        } else {
            config.absolute_step
        };
        let mut xp = x.clone();
        xp[i] += h;
        gradient[i] = if config.central {
            let mut xm = x.clone();
            xm[i] -= h;
            (objective(&xp) - objective(&xm)) / (2.0 * h)
        } else {
            (objective(&xp) - base) / h
        };
    }
    Ok(gradient)
}

pub fn gradient_descent_update(
    parameters: &[f64],
    gradient: &[f64],
    step_size: f64,
) -> Result<Vec<f64>> {
    require(
        parameters.len() == gradient.len() && !parameters.is_empty() && step_size >= 0.0,
        "invalid gradient descent update",
    )?;
    Ok(parameters
        .iter()
        .zip(gradient)
        .map(|(p, g)| p - step_size * g)
        .collect())
}
