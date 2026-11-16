//! Kernel-tuning database (TSV), NUMA discovery and thread-placement planning, and a CPU
//! reference of the SYCL CSR linear-algebra facade, ported from upstream `cfd/core/{autotune,
//! numa, sycl_sparse}` (MIT, cfd_solvers contributors).
use super::{
    bicgstab, preconditioned_conjugate_gradient, CsrMatrix, IdentityPreconditioner,
    IterativeSolverResult, LinearConfig,
};
use crate::{lbm::format_general, require, Result, SolverError};
use std::path::Path;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DeviceTuningKey {
    pub vendor: String,
    pub name: String,
    pub driver: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KernelTuningCandidate {
    pub work_group_size: usize,
    pub vector_width: usize,
    pub fused_stages: usize,
}

impl Default for KernelTuningCandidate {
    fn default() -> Self {
        Self {
            work_group_size: 1,
            vector_width: 1,
            fused_stages: 1,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct KernelTuningRecord {
    pub device: DeviceTuningKey,
    pub kernel_group: String,
    pub candidate: KernelTuningCandidate,
    pub score_seconds: f64,
}

fn clean_field(v: &str) -> String {
    v.replace(['\t', '\n', '\r'], " ")
}

pub fn device_tuning_key_string(key: &DeviceTuningKey) -> String {
    format!(
        "{}|{}|{}",
        clean_field(&key.vendor),
        clean_field(&key.name),
        clean_field(&key.driver)
    )
}

fn parse_size(text: &str, what: &str) -> Result<usize> {
    match text.parse::<usize>() {
        Ok(v) if v > 0 => Ok(v),
        _ => Err(SolverError::Io(format!("invalid tuning {what}"))),
    }
}

/// Best-known launch parameters per (device, kernel group).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct KernelTuningDatabase {
    records: Vec<KernelTuningRecord>,
}

impl KernelTuningDatabase {
    /// Parses the TSV text written by [`Self::to_tsv`] (`#` lines are comments).
    pub fn from_tsv(text: &str) -> Result<Self> {
        let mut db = Self::default();
        for (n, line) in text.split('\n').enumerate() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            require(
                f.len() == 8,
                format!("invalid tuning database row at line {}", n + 1),
            )?;
            let candidate = KernelTuningCandidate {
                work_group_size: parse_size(f[4], "work-group size")?,
                vector_width: parse_size(f[5], "vector width")?,
                fused_stages: parse_size(f[6], "fused stage count")?,
            };
            let score = f[7]
                .parse::<f64>()
                .ok()
                .filter(|s| *s > 0.0 && s.is_finite())
                .ok_or_else(|| {
                    SolverError::Io(format!("invalid tuning score at line {}", n + 1))
                })?;
            db.upsert(KernelTuningRecord {
                device: DeviceTuningKey {
                    vendor: f[0].into(),
                    name: f[1].into(),
                    driver: f[2].into(),
                },
                kernel_group: f[3].into(),
                candidate,
                score_seconds: score,
            })?;
        }
        Ok(db)
    }
    /// A missing file gives an empty database, as upstream.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        match std::fs::read_to_string(path) {
            Ok(text) => Self::from_tsv(&text),
            Err(_) => Ok(Self::default()),
        }
    }
    pub fn to_tsv(&self) -> String {
        let mut o = String::from("# vendor\tdevice\tdriver\tkernel_group\twork_group\tvector_width\tfused_stages\tscore_seconds\n");
        for r in &self.records {
            o += &format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                clean_field(&r.device.vendor),
                clean_field(&r.device.name),
                clean_field(&r.device.driver),
                clean_field(&r.kernel_group),
                r.candidate.work_group_size,
                r.candidate.vector_width,
                r.candidate.fused_stages,
                format_general(r.score_seconds, 17)
            );
        }
        o
    }
    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
            std::fs::create_dir_all(parent).map_err(|e| {
                SolverError::Io(format!("cannot create tuning database directory: {e}"))
            })?;
        }
        std::fs::write(path, self.to_tsv()).map_err(|e| {
            SolverError::Io(format!(
                "cannot create tuning database: {} ({e})",
                path.display()
            ))
        })
    }
    /// Keeps the faster record per (device, kernel group).
    pub fn upsert(&mut self, record: KernelTuningRecord) -> Result<()> {
        let c = record.candidate;
        require(
            !record.kernel_group.is_empty()
                && c.work_group_size > 0
                && c.vector_width > 0
                && c.fused_stages > 0
                && record.score_seconds > 0.0
                && record.score_seconds.is_finite(),
            "invalid kernel tuning record",
        )?;
        match self
            .records
            .iter_mut()
            .find(|r| r.device == record.device && r.kernel_group == record.kernel_group)
        {
            None => self.records.push(record),
            Some(existing) if record.score_seconds <= existing.score_seconds => *existing = record,
            Some(_) => {}
        }
        Ok(())
    }
    pub fn lookup(
        &self,
        device: &DeviceTuningKey,
        kernel_group: &str,
    ) -> Option<KernelTuningRecord> {
        self.records
            .iter()
            .find(|r| r.device == *device && r.kernel_group == kernel_group)
            .cloned()
    }
    pub fn records(&self) -> &[KernelTuningRecord] {
        &self.records
    }
    /// Benchmarks every candidate (seconds, lower is better) and stores the fastest.
    pub fn autotune(
        &mut self,
        device: &DeviceTuningKey,
        kernel_group: &str,
        candidates: &[KernelTuningCandidate],
        benchmark: &mut dyn FnMut(&KernelTuningCandidate) -> f64,
    ) -> Result<KernelTuningRecord> {
        require(
            !kernel_group.is_empty() && !candidates.is_empty(),
            "autotune requires a kernel group and at least one candidate",
        )?;
        let mut best: Option<KernelTuningRecord> = None;
        for c in candidates {
            require(
                c.work_group_size > 0 && c.vector_width > 0 && c.fused_stages > 0,
                "autotune candidate values must be positive",
            )?;
            let score = benchmark(c);
            require(
                score > 0.0 && score.is_finite(),
                "autotune benchmark returned an invalid score",
            )?;
            if best.as_ref().is_none_or(|b| score < b.score_seconds) {
                best = Some(KernelTuningRecord {
                    device: device.clone(),
                    kernel_group: kernel_group.into(),
                    candidate: *c,
                    score_seconds: score,
                });
            }
        }
        let best =
            best.ok_or_else(|| SolverError::InvalidInput("autotune produced no record".into()))?;
        self.upsert(best.clone())?;
        Ok(best)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NumaPlacementPolicy {
    #[default]
    Inherit,
    Compact,
    Spread,
    ExplicitNode,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NumaNodeInfo {
    pub id: i32,
    pub cpus: Vec<u32>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NumaPlacementPlan {
    pub policy: NumaPlacementPolicy,
    pub thread_cpus: Vec<u32>,
    pub thread_nodes: Vec<i32>,
}

impl NumaPlacementPlan {
    pub fn binds_threads(&self) -> bool {
        !self.thread_cpus.is_empty()
    }
    pub fn thread_count(&self) -> usize {
        self.thread_cpus.len()
    }
}

fn parse_uint(token: &str) -> Result<u32> {
    require(!token.is_empty(), "empty CPU index")?;
    require(
        token.bytes().all(|b| b.is_ascii_digit()),
        "invalid CPU list token",
    )?;
    token
        .parse()
        .map_err(|_| SolverError::InvalidInput("invalid CPU list token".into()))
}

/// Linux `cpulist` syntax (`0-3,8,10-11`), sorted and deduplicated.
pub fn parse_cpu_list(text: &str) -> Result<Vec<u32>> {
    let ws = [' ', '\t', '\n'];
    let mut cpus = Vec::new();
    let mut rest = text;
    loop {
        rest = rest.trim_start_matches(ws);
        if rest.is_empty() {
            break;
        }
        let (token, next) = match rest.find(',') {
            Some(c) => (&rest[..c], Some(&rest[c + 1..])),
            None => (rest, None),
        };
        let token = token.trim_end_matches(ws);
        match token.split_once('-') {
            None => cpus.push(parse_uint(token)?),
            Some((a, b)) => {
                let (first, last) = (parse_uint(a)?, parse_uint(b)?);
                require(last >= first, "descending CPU range")?;
                require(u64::from(last - first) < 1 << 22, "CPU range too large")?;
                cpus.extend(first..=last);
            }
        }
        match next {
            Some(n) => rest = n,
            None => break,
        }
    }
    cpus.sort_unstable();
    cpus.dedup();
    Ok(cpus)
}

fn allowed_cpus() -> Vec<u32> {
    let count = std::thread::available_parallelism().map_or(1, |n| n.get()) as u32;
    (0..count).collect()
}

/// `/sys/devices/system/node/node*/cpulist` on Linux, else one node with all CPUs. The allowed
/// CPU set is `0..available_parallelism` (affinity masks need `unsafe` FFI).
pub fn discover_numa_nodes() -> Vec<NumaNodeInfo> {
    let allowed = allowed_cpus();
    let mut nodes = Vec::new();
    if let Ok(entries) = std::fs::read_dir("/sys/devices/system/node") {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            let Some(id) = name
                .strip_prefix("node")
                .filter(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()))
                .and_then(|s| s.parse::<i32>().ok())
            else {
                continue;
            };
            let Ok(text) = std::fs::read_to_string(entry.path().join("cpulist")) else {
                continue;
            };
            let Ok(mut cpus) = parse_cpu_list(text.lines().next().unwrap_or("")) else {
                continue;
            };
            cpus.retain(|c| allowed.contains(c));
            if !cpus.is_empty() {
                nodes.push(NumaNodeInfo { id, cpus });
            }
        }
    }
    if nodes.is_empty() {
        nodes.push(NumaNodeInfo {
            id: 0,
            cpus: allowed,
        });
    }
    nodes.sort_by_key(|n| n.id);
    nodes
}

/// Thread-to-CPU plan over the given nodes (`thread_count = 0` uses every CPU).
pub fn plan_numa_thread_placement_on(
    nodes: &[NumaNodeInfo],
    thread_count: usize,
    policy: NumaPlacementPolicy,
    node_id: Option<i32>,
) -> Result<NumaPlacementPlan> {
    let mut plan = NumaPlacementPlan {
        policy,
        ..Default::default()
    };
    if policy == NumaPlacementPolicy::Inherit {
        return Ok(plan);
    }
    let available: usize = nodes.iter().map(|n| n.cpus.len()).sum();
    require(available > 0, "no CPUs available for NUMA placement")?;
    let threads = if thread_count == 0 {
        available
    } else {
        thread_count
    };
    match policy {
        NumaPlacementPolicy::ExplicitNode => {
            let id = node_id.ok_or_else(|| {
                SolverError::InvalidInput("explicit NUMA placement requires a node id".into())
            })?;
            let node = nodes
                .iter()
                .find(|n| n.id == id && !n.cpus.is_empty())
                .ok_or_else(|| {
                    SolverError::InvalidInput("requested NUMA node is not available".into())
                })?;
            for i in 0..threads {
                plan.thread_cpus.push(node.cpus[i % node.cpus.len()]);
                plan.thread_nodes.push(node.id);
            }
        }
        NumaPlacementPolicy::Compact => {
            let flat: Vec<(u32, i32)> = nodes
                .iter()
                .flat_map(|n| n.cpus.iter().map(move |c| (*c, n.id)))
                .collect();
            for i in 0..threads {
                let (c, n) = flat[i % flat.len()];
                plan.thread_cpus.push(c);
                plan.thread_nodes.push(n);
            }
        }
        NumaPlacementPolicy::Spread => {
            let mut next = vec![0usize; nodes.len()];
            for i in 0..threads {
                let ni = i % nodes.len();
                let node = &nodes[ni];
                require(!node.cpus.is_empty(), "empty NUMA node")?;
                let ci = next[ni] % node.cpus.len();
                next[ni] += 1;
                plan.thread_cpus.push(node.cpus[ci]);
                plan.thread_nodes.push(node.id);
            }
        }
        NumaPlacementPolicy::Inherit => {}
    }
    Ok(plan)
}

pub fn plan_numa_thread_placement(
    thread_count: usize,
    policy: NumaPlacementPolicy,
    node_id: Option<i32>,
) -> Result<NumaPlacementPlan> {
    plan_numa_thread_placement_on(&discover_numa_nodes(), thread_count, policy, node_id)
}

/// Thread affinity needs OS FFI, which `forbid(unsafe_code)` rules out: always `false`.
pub fn bind_current_thread_to_cpu(_cpu: u32) -> bool {
    false
}

/// Succeeds only for plans that do not bind threads.
pub fn apply_numa_thread_placement(plan: &NumaPlacementPlan) -> bool {
    !plan.binds_threads()
}

/// First-touch initialisation (the OS places pages on the touching thread's node).
pub fn numa_first_touch<T: Clone + Send + Sync>(values: &mut [T], value: &T) {
    values.iter_mut().for_each(|v| *v = value.clone());
}

/// CPU reference of upstream's SYCL CSR facade with the same calls; "device" buffers are host
/// slices.
#[derive(Clone, Debug)]
pub struct SyclCsrLinearAlgebra {
    matrix: CsrMatrix,
}

impl SyclCsrLinearAlgebra {
    pub fn new(matrix: &CsrMatrix) -> Self {
        Self {
            matrix: matrix.clone(),
        }
    }
    pub fn rows(&self) -> usize {
        self.matrix.rows()
    }
    pub fn cols(&self) -> usize {
        self.matrix.cols()
    }
    pub fn nonzeros(&self) -> usize {
        self.matrix.nonzeros()
    }
    pub fn device_name(&self) -> String {
        "host CPU (SYCL not linked)".into()
    }
    pub fn multiply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
        self.matrix.multiply(x, y)
    }
    pub fn update_values(&mut self, values: &[f64]) -> Result<()> {
        require(
            values.len() == self.matrix.nonzeros(),
            "SYCL CSR value update size mismatch",
        )?;
        self.matrix = CsrMatrix::new(
            self.matrix.rows(),
            self.matrix.cols(),
            self.matrix.row_offsets().to_vec(),
            self.matrix.columns().to_vec(),
            values.to_vec(),
        )?;
        Ok(())
    }
    pub fn conjugate_gradient(
        &self,
        rhs: &[f64],
        x: &mut [f64],
        max_iterations: usize,
        relative_tolerance: f64,
    ) -> Result<IterativeSolverResult> {
        let c = LinearConfig {
            max_iterations,
            relative_tolerance,
            ..Default::default()
        };
        preconditioned_conjugate_gradient(
            rhs,
            x,
            |v, o| self.matrix.multiply(v, o),
            &IdentityPreconditioner,
            c,
        )
    }
    pub fn bicgstab(
        &self,
        rhs: &[f64],
        x: &mut [f64],
        max_iterations: usize,
        relative_tolerance: f64,
    ) -> Result<IterativeSolverResult> {
        let c = LinearConfig {
            max_iterations,
            relative_tolerance,
            ..Default::default()
        };
        bicgstab(
            rhs,
            x,
            |v, o| self.matrix.multiply(v, o),
            &IdentityPreconditioner,
            c,
        )
    }
}
