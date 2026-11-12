//! Solver-independent Multiphysics authoring and portable run evidence.
//!
//! This crate does not compute geometry or invoke a numerical solver. Saved incomplete
//! cases are bounded authoring records. The companion solver validates its full schema
//! and physics when execution is requested. See `README.md` for the admission boundary.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

use std::collections::BTreeSet;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};

mod templates;

pub use templates::templates;

/// Pinned producer whose case and run schemas this initial adapter understands.
pub const SOLVER_VERSION: &str = "0.2.0";
pub const SOLVER_REVISION: &str = "7cceaaaddafaf2515e01a86406e13d6037fdbf9a";
pub const MAX_STUDIES: usize = 64;
pub const MAX_RUNS: usize = 12;
pub const MAX_CASE_BYTES: usize = 2 * 1024 * 1024;
pub const MAX_ARTIFACT_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_RUN_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_ARTIFACTS: usize = 128;
const MAX_PROJECT_BYTES: usize = 64 * 1024 * 1024;
const MAX_DEPTH: usize = 32;
const MAX_NODES: usize = 200_000;
const MAX_CELLS: u64 = 32_768;
const MAX_STEPS: u64 = 20_000;
const MAX_WORK: u64 = 50_000_000;

#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum Error {
    #[error("invalid Multiphysics authoring data: {0}")]
    Invalid(String),
    #[error("Multiphysics resource limit: {0}")]
    Limit(String),
    #[error("Multiphysics study not found: {0}")]
    NotFound(String),
    #[error("the study changed; reload the draft before applying")]
    Stale,
    #[error("invalid Multiphysics JSON: {0}")]
    Json(String),
}

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Project {
    pub studies: Vec<Study>,
    pub next_id: u64,
}

impl Default for Project {
    fn default() -> Self {
        Self { studies: vec![], next_id: 1 }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Study {
    pub id: String,
    pub name: String,
    #[serde(default)]
    pub question: String,
    pub profile: String,
    pub revision: u64,
    pub case: Value,
    #[serde(default)]
    pub runs: Vec<RunEvidence>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunEvidence {
    pub id: String,
    pub study_revision: u64,
    pub case_hash: String,
    pub input_case: Value,
    pub manifest: Value,
    pub artifacts: Vec<Artifact>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Artifact {
    pub name: String,
    pub text: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Template {
    pub id: String,
    pub name: String,
    pub profile: String,
    pub case: Value,
}

/// Bounded model checks only. An empty errors list does not qualify numerical execution.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Validation {
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
}

impl Project {
    pub fn is_empty(&self) -> bool {
        self.studies.is_empty()
    }
    pub fn study(&self, id: &str) -> Option<&Study> {
        self.studies.iter().find(|study| study.id == id)
    }

    pub fn create(&mut self, template: &str, name: &str) -> Result<Study> {
        self.create_with_budget(template, name, MAX_PROJECT_BYTES)
    }

    fn create_with_budget(&mut self, template: &str, name: &str, budget: usize) -> Result<Study> {
        self.validate_with_budget(budget)?;
        if self.studies.len() >= MAX_STUDIES {
            return Err(Error::Limit("too many studies".into()));
        }
        let template =
            templates()?.into_iter().find(|item| item.id == template).ok_or_else(|| Error::Invalid(format!("unknown case template {template}")))?;
        let next_id = self.next_id.checked_add(1).ok_or_else(|| Error::Limit("study identifiers exhausted".into()))?;
        let study = Study {
            id: format!("mp-study-{}", self.next_id),
            name: name.trim().into(),
            question: String::new(),
            profile: template.profile,
            revision: 1,
            case: template.case,
            runs: vec![],
        };
        validate_study(&study)?;
        let mut candidate = self.clone();
        candidate.studies.push(study.clone());
        candidate.next_id = next_id;
        candidate.validate_with_budget(budget)?;
        *self = candidate;
        Ok(study)
    }

    /// Apply one authored revision atomically. Existing immutable runs become historical.
    pub fn update(&mut self, id: &str, expected_revision: u64, name: &str, question: &str, case: Value) -> Result<Study> {
        self.update_with_budget(id, expected_revision, name, question, case, MAX_PROJECT_BYTES)
    }

    fn update_with_budget(&mut self, id: &str, expected_revision: u64, name: &str, question: &str, case: Value, budget: usize) -> Result<Study> {
        self.validate_with_budget(budget)?;
        let current = self.study(id).ok_or_else(|| Error::NotFound(id.into()))?;
        if current.revision != expected_revision {
            return Err(Error::Stale);
        }
        if current.name == name.trim() && current.question == question && current.case == case {
            return Ok(current.clone());
        }
        let mut proposed = current.clone();
        proposed.name = name.trim().into();
        proposed.question = question.into();
        proposed.case = case;
        proposed.revision = proposed.revision.checked_add(1).ok_or_else(|| Error::Limit("study revisions exhausted".into()))?;
        validate_study(&proposed)?;
        let mut candidate = self.clone();
        let target = candidate.studies.iter_mut().find(|study| study.id == id).ok_or_else(|| Error::NotFound(id.into()))?;
        *target = proposed.clone();
        candidate.validate_with_budget(budget)?;
        *self = candidate;
        Ok(proposed)
    }

    pub fn remove(&mut self, id: &str) -> Result<Study> {
        self.validate()?;
        let index = self.studies.iter().position(|study| study.id == id).ok_or_else(|| Error::NotFound(id.into()))?;
        let mut candidate = self.clone();
        let removed = candidate.studies.remove(index);
        candidate.validate()?;
        *self = candidate;
        Ok(removed)
    }

    /// Admit one completed, current-revision run after the host validates its origin.
    pub fn attach_evidence(&mut self, study_id: &str, evidence: RunEvidence) -> Result<Study> {
        self.attach_with_budget(study_id, evidence, MAX_PROJECT_BYTES)
    }

    fn attach_with_budget(&mut self, study_id: &str, evidence: RunEvidence, budget: usize) -> Result<Study> {
        self.validate_with_budget(budget)?;
        let current = self.study(study_id).ok_or_else(|| Error::NotFound(study_id.into()))?;
        if current.runs.len() >= MAX_RUNS {
            return Err(Error::Limit("run history is full; remove an old run explicitly".into()));
        }
        if current.runs.iter().any(|run| run.id == evidence.id) {
            return Err(Error::Invalid("duplicate run identifier".into()));
        }
        validate_evidence(current, &evidence)?;
        let mut candidate = self.clone();
        let study = candidate.studies.iter_mut().find(|study| study.id == study_id).ok_or_else(|| Error::NotFound(study_id.into()))?;
        study.runs.push(evidence);
        let result = study.clone();
        candidate.validate_with_budget(budget)?;
        *self = candidate;
        Ok(result)
    }

    pub fn validate(&self) -> Result<()> {
        self.validate_with_budget(MAX_PROJECT_BYTES)
    }

    fn validate_with_budget(&self, budget: usize) -> Result<()> {
        if self.studies.len() > MAX_STUDIES {
            return Err(Error::Limit("too many studies".into()));
        }
        if self.next_id == 0 {
            return Err(Error::Invalid("next study identifier must be positive".into()));
        }
        let mut ids = BTreeSet::new();
        for study in &self.studies {
            if !ids.insert(&study.id) {
                return Err(Error::Invalid("duplicate study identifier".into()));
            }
            let ordinal = study
                .id
                .strip_prefix("mp-study-")
                .and_then(|id| id.parse::<u64>().ok())
                .ok_or_else(|| Error::Invalid("malformed study identifier".into()))?;
            if ordinal == 0 || ordinal >= self.next_id {
                return Err(Error::Invalid("study identifier is outside the project sequence".into()));
            }
            validate_study(study)?;
        }
        // Count the complete serialized payload, including both copies of run case
        // data, manifest metadata and JSON escaping. No extra payload buffer is made.
        let mut counter = ByteBudget { used: 0, limit: budget };
        serde_json::to_writer(&mut counter, self).map_err(|_| Error::Limit(format!("project payload exceeds {budget} bytes")))?;
        Ok(())
    }
}

struct ByteBudget {
    used: usize,
    limit: usize,
}
impl std::io::Write for ByteBudget {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let next = self
            .used
            .checked_add(bytes.len())
            .filter(|next| *next <= self.limit)
            .ok_or_else(|| std::io::Error::other("project byte budget exceeded"))?;
        self.used = next;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Original MIT design inventory, not installed-worker discovery or a solver qualification.
pub fn catalog() -> Result<Value> {
    parse_json(include_str!("../data/solver-catalog.json"))
}

pub fn case_hash(case: &Value) -> Result<String> {
    validate_json(case, MAX_CASE_BYTES)?;
    let bytes = serde_json::to_vec(case).map_err(|error| Error::Json(error.to_string()))?;
    Ok(format!("{:x}", Sha256::digest(bytes)))
}

/// Only bounded object structure and resource requests are checked here.
/// Missing physical settings and unknown solver-specific keys remain savable drafts.
pub fn validate_case(case: &Value) -> Result<()> {
    validate_json(case, MAX_CASE_BYTES)?;
    if !case.is_object() {
        return Err(Error::Invalid("case must be a JSON object".into()));
    }
    let mut cell_bound = 1u64;
    let mut step_bound = 1u64;
    let mut stack = vec![case];
    while let Some(value) = stack.pop() {
        match value {
            Value::Object(object) => {
                let mut local_cells = 1u64;
                let mut has_dimensions = false;
                for key in ["nx", "ny", "nz", "elements", "cells", "nodes"] {
                    if let Some(number) = object.get(key).and_then(Value::as_u64) {
                        has_dimensions = true;
                        local_cells = local_cells.checked_mul(number).ok_or_else(|| Error::Limit("mesh size overflow".into()))?;
                    }
                }
                if let Some(dimensions) = object.get("dimensions").and_then(Value::as_array) {
                    if dimensions.len() > 3 {
                        return Err(Error::Limit("mesh dimensions exceed three axes".into()));
                    }
                    has_dimensions = true;
                    for dimension in dimensions {
                        if let Some(number) = dimension.as_u64() {
                            local_cells = local_cells.checked_mul(number).ok_or_else(|| Error::Limit("mesh size overflow".into()))?;
                        }
                    }
                }
                if has_dimensions {
                    cell_bound = cell_bound.max(local_cells);
                }
                for (key, value) in object {
                    if matches!(key.as_str(), "steps" | "iterations" | "max_iterations" | "max_steps")
                        && let Some(number) = value.as_u64()
                    {
                        if number > MAX_STEPS {
                            return Err(Error::Limit(format!("{key} exceeds {MAX_STEPS}")));
                        }
                        // Linear iterations are bounded independently, not multiplied by physical steps.
                        if matches!(key.as_str(), "steps" | "iterations") {
                            step_bound = step_bound.max(number);
                        }
                    }
                    stack.push(value);
                }
            }
            Value::Array(array) => stack.extend(array),
            _ => {}
        }
    }
    if cell_bound > MAX_CELLS {
        return Err(Error::Limit(format!("generated mesh exceeds {MAX_CELLS} cells/nodes")));
    }
    if cell_bound.checked_mul(step_bound).is_none_or(|work| work > MAX_WORK) {
        return Err(Error::Limit(format!("mesh times step count exceeds {MAX_WORK}")));
    }
    Ok(())
}

pub fn preflight(study: &Study) -> Validation {
    let mut validation = Validation::default();
    if let Err(error) = validate_authoring(study) {
        validation.errors.push(error.to_string());
    }
    match study.case.get("solver").and_then(Value::as_str) {
        Some(solver) => match profile_for_solver(solver) {
            Some(profile) if profile == study.profile => {}
            Some(_) => validation.errors.push("case solver does not belong to the selected profile".into()),
            None => validation.errors.push("this solver has no native case adapter".into()),
        },
        None => validation.errors.push("choose a solver in the case object".into()),
    }
    validation.warnings.push(
        "Only bounded authoring structure is checked. The companion validates the full numerical case and physical prerequisites at execution."
            .into(),
    );
    validation
        .warnings
        .push("Cases use their declared SI or lattice units. CAD-to-volume meshing and arbitrary CAD study execution are unavailable.".into());
    validation.warnings.push("Successful execution is separate from linear convergence, steady convergence and physical validation.".into());
    validation
}

/// Validate current-result association and portable artifact structure, not physical truth.
pub fn validate_evidence(study: &Study, evidence: &RunEvidence) -> Result<()> {
    validate_record(evidence)?;
    if evidence.study_revision != study.revision || evidence.case_hash != case_hash(&study.case)? || evidence.input_case != study.case {
        return Err(Error::Stale);
    }
    if study.case.get("solver").and_then(Value::as_str).and_then(profile_for_solver) != Some(study.profile.as_str()) {
        return Err(Error::Invalid("run does not belong to the selected profile".into()));
    }
    Ok(())
}

/// Validate retained current or historical evidence without attaching it to a new revision.
pub fn validate_record(evidence: &RunEvidence) -> Result<()> {
    validate_run_structure(evidence)
}

pub fn profile_for_solver(solver: &str) -> Option<&'static str> {
    Some(match solver {
        "poly_scalar_transport3d" | "collocated_flow3d" | "rans_transport3d" | "diffusion2d" => "fvm-incompressible",
        "solid_heat3d" | "thermal_transport3d" => "thermal-fvm",
        "poisson1d" | "poisson2d" | "poisson3d" | "scalar_diffusion2d" | "heat2d" => "fem-heat2d",
        "lbm_d2q9" => "lbm-fluid",
        "maxwell1d" => "fdtd-em",
        "euler1d" => "fvm-compressible",
        "rf_tem" => "si-pi-emc",
        "optics_sequential" => "optical-rays",
        "chemistry_batch" => "chemical-reactors",
        "acoustic_kspace2d" => "acoustic-kspace2d",
        _ => return None,
    })
}

fn validate_study(study: &Study) -> Result<()> {
    validate_authoring(study)?;
    if study.runs.len() > MAX_RUNS {
        return Err(Error::Limit("too many retained runs".into()));
    }
    let mut ids = BTreeSet::new();
    for run in &study.runs {
        if !ids.insert(&run.id) {
            return Err(Error::Invalid("duplicate run identifier".into()));
        }
        if run.study_revision > study.revision {
            return Err(Error::Invalid("run belongs to a future study revision".into()));
        }
        if run.study_revision == study.revision {
            validate_evidence(study, run)?;
        } else {
            validate_record(run)?;
        }
    }
    Ok(())
}

fn validate_authoring(study: &Study) -> Result<()> {
    check_text(&study.id, "study identifier", 80, false)?;
    check_text(&study.name, "study name", 160, false)?;
    check_text(&study.question, "engineering question", 4000, true)?;
    check_text(&study.profile, "profile", 80, false)?;
    if study.revision == 0 {
        return Err(Error::Invalid("study revision must be positive".into()));
    }
    validate_case(&study.case)
}

fn validate_run_structure(run: &RunEvidence) -> Result<()> {
    check_text(&run.id, "run identifier", 100, false)?;
    if run.study_revision == 0
        || run.case_hash.len() != 64
        || !run.case_hash.bytes().all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        return Err(Error::Invalid("run has an invalid input identity".into()));
    }
    validate_case(&run.input_case)?;
    if run.case_hash != case_hash(&run.input_case)? {
        return Err(Error::Invalid("retained input case digest differs".into()));
    }
    validate_json(&run.manifest, MAX_CASE_BYTES + 256 * 1024)?;
    let manifest = &run.manifest;
    for (key, expected) in [
        ("schema", "cfd-solvers-rs/run/1"),
        ("rust_port_version", SOLVER_VERSION),
        ("upstream_revision", SOLVER_REVISION),
        ("backend", "native_rust_serial_cpu"),
    ] {
        if manifest.get(key).and_then(Value::as_str) != Some(expected) {
            return Err(Error::Invalid(format!("unsupported run {key}")));
        }
    }
    if manifest.get("execution_completed").and_then(Value::as_bool) != Some(true)
        || manifest.get("complete_clone").and_then(Value::as_bool) != Some(false)
    {
        return Err(Error::Invalid("run must be completed and declare its partial-port scope".into()));
    }
    for key in ["linear_solver_converged", "steady_state_converged"] {
        if let Some(value) = manifest.get(key)
            && !(value.is_null() || value.is_boolean())
        {
            return Err(Error::Invalid(format!("invalid {key} evidence")));
        }
    }
    let case = manifest.get("case").ok_or_else(|| Error::Invalid("run has no case".into()))?;
    // Producer defaults describe its schema, not work requested by the author.
    // The original input above retains all execution-resource admission checks.
    validate_json(case, MAX_CASE_BYTES)?;
    if !case.is_object() {
        return Err(Error::Invalid("run case must be a JSON object".into()));
    }
    if !normalized_authored_matches(&run.input_case, case) {
        return Err(Error::Invalid("run case differs from explicitly authored inputs".into()));
    }
    let solver = manifest.get("solver").and_then(Value::as_str).ok_or_else(|| Error::Invalid("run has no solver".into()))?;
    if profile_for_solver(solver).is_none()
        || case.get("solver").and_then(Value::as_str) != Some(solver)
        || run.input_case.get("solver").and_then(Value::as_str) != Some(solver)
    {
        return Err(Error::Invalid("run solver/case mismatch".into()));
    }
    let metrics = manifest.get("metrics").and_then(Value::as_object).ok_or_else(|| Error::Invalid("run has no metric table".into()))?;
    if metrics.len() > 256 || metrics.values().any(|value| value.as_f64().is_none_or(|number| !number.is_finite())) {
        return Err(Error::Invalid("run metrics must be finite scalars".into()));
    }
    let files = manifest.get("files").and_then(Value::as_array).ok_or_else(|| Error::Invalid("run has no artifact list".into()))?;
    if files.is_empty() || files.len() > MAX_ARTIFACTS || run.artifacts.len() != files.len() {
        return Err(Error::Limit("incomplete or excessive run artifacts".into()));
    }
    let mut expected = BTreeSet::new();
    for file in files {
        let file = file.as_str().ok_or_else(|| Error::Invalid("artifact name must be text".into()))?;
        validate_artifact_name(file)?;
        if !expected.insert(file) {
            return Err(Error::Invalid("duplicate manifest artifact".into()));
        }
    }
    if !expected.contains("manifest.json") {
        return Err(Error::Invalid("run must retain its manifest artifact".into()));
    }
    let mut seen = BTreeSet::new();
    let mut total = 0usize;
    for artifact in &run.artifacts {
        validate_artifact_name(&artifact.name)?;
        if !expected.contains(artifact.name.as_str()) || !seen.insert(artifact.name.as_str()) {
            return Err(Error::Invalid("unexpected or duplicate artifact".into()));
        }
        if artifact.text.len() > MAX_ARTIFACT_BYTES {
            return Err(Error::Limit("artifact exceeds 8 MiB".into()));
        }
        total = total.checked_add(artifact.text.len()).ok_or_else(|| Error::Limit("artifact byte count overflow".into()))?;
        if total > MAX_RUN_BYTES {
            return Err(Error::Limit("portable run exceeds 16 MiB".into()));
        }
        if artifact.name == "manifest.json" && parse_json(&artifact.text)? != *manifest {
            return Err(Error::Invalid("manifest artifact differs from retained metadata".into()));
        }
    }
    Ok(())
}

/// Portable artifacts are flat filenames. The host rechecks directory and symlink ownership.
pub fn validate_artifact_name(name: &str) -> Result<()> {
    if name.is_empty()
        || name.len() > 160
        || name.starts_with('.')
        || name.ends_with('.')
        || !name.bytes().all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
    {
        return Err(Error::Invalid("artifact must have a safe ordinary filename".into()));
    }
    let stem = name.split('.').next().unwrap_or_default().to_ascii_uppercase();
    if matches!(stem.as_str(), "CON" | "PRN" | "AUX" | "NUL")
        || stem
            .strip_prefix("COM")
            .or_else(|| stem.strip_prefix("LPT"))
            .is_some_and(|number| matches!(number, "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"))
    {
        return Err(Error::Invalid("artifact name is a reserved device filename".into()));
    }
    Ok(())
}

fn check_text(value: &str, label: &str, limit: usize, empty: bool) -> Result<()> {
    if value.len() > limit {
        return Err(Error::Limit(format!("{label} is too long")));
    }
    if (!empty && value.trim().is_empty()) || value.chars().any(|character| character.is_control() && !matches!(character, '\n' | '\t')) {
        return Err(Error::Invalid(format!("invalid {label}")));
    }
    Ok(())
}

fn validate_json(value: &Value, max_bytes: usize) -> Result<()> {
    let mut stack = vec![(value, 0usize)];
    let mut nodes = 0usize;
    let mut text_bytes = 0usize;
    while let Some((value, depth)) = stack.pop() {
        nodes += 1;
        if depth > MAX_DEPTH || nodes > MAX_NODES {
            return Err(Error::Limit("JSON nesting or node limit exceeded".into()));
        }
        match value {
            Value::Object(object) => {
                for (key, child) in object {
                    text_bytes = text_bytes.saturating_add(key.len());
                    stack.push((child, depth + 1));
                }
            }
            Value::Array(array) => {
                for child in array {
                    stack.push((child, depth + 1));
                }
            }
            Value::String(text) => {
                text_bytes = text_bytes.saturating_add(text.len());
            }
            Value::Number(number) if number.as_f64().is_none_or(|number| !number.is_finite()) => {
                return Err(Error::Invalid("nonfinite number".into()));
            }
            _ => {}
        }
        if text_bytes > max_bytes {
            return Err(Error::Limit("JSON text limit exceeded".into()));
        }
    }
    let bytes = serde_json::to_vec(value).map_err(|error| Error::Json(error.to_string()))?;
    if bytes.len() > max_bytes {
        return Err(Error::Limit("JSON byte limit exceeded".into()));
    }
    Ok(())
}

fn authored_subset(authored: &Value, produced: &Value) -> bool {
    match (authored, produced) {
        (Value::Object(a), Value::Object(b)) => a.iter().all(|(key, value)| b.get(key).is_some_and(|other| authored_subset(value, other))),
        (Value::Array(a), Value::Array(b)) => a.len() == b.len() && a.iter().zip(b).all(|(a, b)| authored_subset(a, b)),
        (Value::Number(a), Value::Number(b)) => {
            if a == b {
                return true;
            }
            // Preserve large exact integers; only compare mixed encodings in the exact f64 range.
            const EXACT_INTEGER_LIMIT: u64 = 1 << 53;
            if a.as_u64().is_some_and(|number| number > EXACT_INTEGER_LIMIT)
                || b.as_u64().is_some_and(|number| number > EXACT_INTEGER_LIMIT)
                || a.as_i64().is_some_and(|number| number < -(EXACT_INTEGER_LIMIT as i64))
                || b.as_i64().is_some_and(|number| number < -(EXACT_INTEGER_LIMIT as i64))
            {
                return false;
            }
            a.as_f64().zip(b.as_f64()).is_some_and(|(a, b)| a == b)
        }
        _ => authored == produced,
    }
}

fn normalized_authored_matches(authored: &Value, produced: &Value) -> bool {
    if authored.get("solver").and_then(Value::as_str) != Some("lbm_d2q9") {
        return authored_subset(authored, produced);
    }
    // The pinned D2Q9 schema declares precisely these scalars as f32. Match its
    // actual rounding, without a tolerance on unrelated settings or identities.
    let mut normalized = authored.clone();
    let mut normalized_produced = produced.clone();
    for path in [
        "/config/tau",
        "/config/acceleration_x",
        "/config/acceleration_y",
        "/initial/rho",
        "/initial/ux",
        "/initial/uy",
        "/initial/amplitude",
        "/top_wall_speed",
    ] {
        if let Some(value) = normalized.pointer_mut(path) {
            let Some(number) = value.as_f64() else {
                return false;
            };
            let Some(number) = serde_json::Number::from_f64(number as f32 as f64) else {
                return false;
            };
            *value = Value::Number(number);
            let Some(produced_value) = normalized_produced.pointer_mut(path) else {
                return false;
            };
            let Some(produced_number) = produced_value.as_f64() else {
                return false;
            };
            let Some(produced_number) = serde_json::Number::from_f64(produced_number as f32 as f64) else {
                return false;
            };
            *produced_value = Value::Number(produced_number);
        }
    }
    authored_subset(&normalized, &normalized_produced)
}

fn parse_json(text: &str) -> Result<Value> {
    serde_json::from_str(text).map_err(|error| Error::Json(error.to_string()))
}

#[cfg(test)]
mod tests;
