//! Runtime infrastructure: backend selection, CPU parallel helpers, slab decomposition, device
//! transfer accounting, distributed CSR partitions and owner-routed assembly, f32-inner
//! iterative refinement, a thread-safe profiler, scratch arenas, SIMD packs, static SoA fields,
//! dependency task graphs, Tri3 mesh partitioning, distributed reverse-mode AD and an AMG
//! preconditioner adapter, ported from upstream `cfd/core/{backend, parallel, decomposition,
//! device_residency, distributed_csr, distributed_sparse_assembly, mixed_precision, profiler,
//! scratch_arena, simd, soa_field, task_graph, mesh_partition, reverse_ad_distributed,
//! amg_adapter, aligned_allocator}` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::{CsrMatrix, Preconditioner, ReverseScalar, ReverseTape};
use crate::{require, Result, SolverError};
use std::{
    collections::HashMap,
    ops::{Add, AddAssign, Index, IndexMut, Mul, Sub},
    sync::{Mutex, OnceLock},
    time::Instant,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RuntimeBackend {
    Serial,
    /// CPU threads: rayon when built with the `parallel` feature.
    OpenMp,
    Mpi,
    Sycl,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BackendInfo {
    pub backend: RuntimeBackend,
    pub name: &'static str,
    pub available: bool,
}

pub fn backend_name(b: RuntimeBackend) -> &'static str {
    match b {
        RuntimeBackend::Serial => "serial",
        RuntimeBackend::OpenMp => "openmp",
        RuntimeBackend::Mpi => "mpi",
        RuntimeBackend::Sycl => "sycl",
    }
}

/// Serial always; CPU threads with the `parallel` feature; MPI and SYCL are never linked here.
pub fn backend_available(b: RuntimeBackend) -> bool {
    match b {
        RuntimeBackend::Serial => true,
        RuntimeBackend::OpenMp => cfg!(feature = "parallel"),
        RuntimeBackend::Mpi | RuntimeBackend::Sycl => false,
    }
}

pub fn compiled_backends() -> Vec<BackendInfo> {
    [
        RuntimeBackend::Serial,
        RuntimeBackend::OpenMp,
        RuntimeBackend::Mpi,
        RuntimeBackend::Sycl,
    ]
    .into_iter()
    .map(|b| BackendInfo {
        backend: b,
        name: backend_name(b),
        available: backend_available(b),
    })
    .collect()
}

pub fn parse_backend(name: &str) -> Result<RuntimeBackend> {
    match name {
        "serial" | "cpu-serial" => Ok(RuntimeBackend::Serial),
        "openmp" | "omp" | "cpu" => Ok(RuntimeBackend::OpenMp),
        "mpi" => Ok(RuntimeBackend::Mpi),
        "sycl" | "gpu" => Ok(RuntimeBackend::Sycl),
        _ => Err(SolverError::InvalidInput("unknown runtime backend".into())),
    }
}

/// `"auto"` picks the most capable available backend.
pub fn select_backend(request: &str) -> Result<RuntimeBackend> {
    if request == "auto" {
        for b in [
            RuntimeBackend::Sycl,
            RuntimeBackend::Mpi,
            RuntimeBackend::OpenMp,
            RuntimeBackend::Serial,
        ] {
            if backend_available(b) {
                return Ok(b);
            }
        }
    }
    let b = parse_backend(request)?;
    if !backend_available(b) {
        return Err(SolverError::Unsupported(
            "requested runtime backend was not compiled in".into(),
        ));
    }
    Ok(b)
}

pub fn cpu_thread_capacity() -> usize {
    #[cfg(feature = "parallel")]
    {
        rayon::current_num_threads().max(1)
    }
    #[cfg(not(feature = "parallel"))]
    {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    }
}

/// Rayon's global pool is fixed after first use, so this is advisory (as upstream without OpenMP).
pub fn set_cpu_threads(_count: usize) {}

/// At least 512 iterations per thread and more than one thread (only with `parallel`).
pub fn cpu_parallel_worthwhile(count: usize) -> bool {
    if !cfg!(feature = "parallel") {
        return false;
    }
    let threads = cpu_thread_capacity();
    threads > 1 && count / threads >= 512
}

pub fn parallel_for(count: usize, f: impl Fn(usize) + Sync + Send) {
    #[cfg(feature = "parallel")]
    if cpu_parallel_worthwhile(count) {
        use rayon::prelude::*;
        (0..count).into_par_iter().for_each(f);
        return;
    }
    (0..count).for_each(f);
}

/// Serial sums are bit-identical to upstream's serial build; parallel sums reorder additions.
pub fn parallel_sum(count: usize, f: impl Fn(usize) -> f64 + Sync + Send) -> f64 {
    #[cfg(feature = "parallel")]
    if cpu_parallel_worthwhile(count) {
        use rayon::prelude::*;
        return (0..count).into_par_iter().map(f).sum();
    }
    let mut v = 0.0;
    for i in 0..count {
        v += f(i);
    }
    v
}

/// Maximum of `f(i)` and zero.
pub fn parallel_max(count: usize, f: impl Fn(usize) -> f64 + Sync + Send) -> f64 {
    #[cfg(feature = "parallel")]
    if cpu_parallel_worthwhile(count) {
        use rayon::prelude::*;
        return (0..count).into_par_iter().map(f).reduce(|| 0.0, f64::max);
    }
    (0..count).fold(0.0, |m, i| f64::max(m, f(i)))
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SlabPartition {
    pub begin: usize,
    pub end: usize,
    pub halo: usize,
}

impl SlabPartition {
    pub fn size(&self) -> usize {
        self.end - self.begin
    }
}

/// Balanced contiguous block for `rank` of `ranks` (upstream default halo 1).
pub fn partition_1d(
    global_size: usize,
    rank: usize,
    ranks: usize,
    halo: usize,
) -> Result<SlabPartition> {
    require(ranks > 0 && rank < ranks, "invalid partition rank/count")?;
    let (base, extra) = (global_size / ranks, global_size % ranks);
    let begin = rank * base + rank.min(extra);
    let local = base + usize::from(rank < extra);
    Ok(SlabPartition {
        begin,
        end: begin + local,
        halo,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DeviceTransferStats {
    pub host_to_device_bytes: u64,
    pub device_to_host_bytes: u64,
    pub device_to_device_bytes: u64,
    pub synchronization_points: u64,
}

impl DeviceTransferStats {
    pub fn explicit_transfer_bytes(&self) -> u64 {
        self.host_to_device_bytes + self.device_to_host_bytes + self.device_to_device_bytes
    }
    pub fn host_transfer_bytes(&self) -> u64 {
        self.host_to_device_bytes + self.device_to_host_bytes
    }
    pub fn record_host_to_device(&mut self, bytes: usize) {
        self.host_to_device_bytes += bytes as u64;
    }
    pub fn record_device_to_host(&mut self, bytes: usize) {
        self.device_to_host_bytes += bytes as u64;
    }
    pub fn record_device_to_device(&mut self, bytes: usize) {
        self.device_to_device_bytes += bytes as u64;
    }
    pub fn record_synchronization(&mut self, count: u64) {
        self.synchronization_points += count;
    }
    pub fn reset(&mut self) {
        *self = Self::default();
    }
}

/// Owned row block of a square CSR matrix with sorted unique off-block (halo) columns.
#[derive(Clone, Debug)]
pub struct DistributedCsrPartition {
    matrix: CsrMatrix,
    row_begin: usize,
    row_end: usize,
    halo_columns: Vec<usize>,
}

impl DistributedCsrPartition {
    pub fn new(global: CsrMatrix, row_begin: usize, row_end: usize) -> Result<Self> {
        require(
            row_begin < row_end && row_end <= global.rows() && global.rows() == global.cols(),
            "invalid distributed CSR row partition",
        )?;
        let mut halo = Vec::new();
        for r in row_begin..row_end {
            for k in global.row_offsets()[r]..global.row_offsets()[r + 1] {
                let c = global.columns()[k];
                if c < row_begin || c >= row_end {
                    halo.push(c);
                }
            }
        }
        halo.sort_unstable();
        halo.dedup();
        Ok(Self {
            matrix: global,
            row_begin,
            row_end,
            halo_columns: halo,
        })
    }
    pub fn local_rows(&self) -> usize {
        self.row_end - self.row_begin
    }
    pub fn row_begin(&self) -> usize {
        self.row_begin
    }
    pub fn halo_columns(&self) -> &[usize] {
        &self.halo_columns
    }
    pub fn multiply(&self, local_x: &[f64], halo_x: &[f64], y: &mut [f64]) -> Result<()> {
        require(
            local_x.len() == self.local_rows()
                && halo_x.len() == self.halo_columns.len()
                && y.len() == self.local_rows(),
            "distributed CSR buffer mismatch",
        )?;
        let m = &self.matrix;
        for lr in 0..self.local_rows() {
            let r = self.row_begin + lr;
            let mut sum = 0.0;
            for k in m.row_offsets()[r]..m.row_offsets()[r + 1] {
                let c = m.columns()[k];
                let x = if c >= self.row_begin && c < self.row_end {
                    local_x[c - self.row_begin]
                } else {
                    let slot = self.halo_columns.binary_search(&c).map_err(|_| {
                        SolverError::NumericalFailure("distributed CSR missing halo column".into())
                    })?;
                    halo_x[slot]
                };
                sum += m.values()[k] * x;
            }
            y[lr] = sum;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DistributedTriplet {
    pub row: usize,
    pub column: usize,
    pub value: f64,
    pub source_rank: usize,
    pub sequence: usize,
}

#[derive(Clone, Debug)]
pub struct DistributedCsrBlock {
    pub row_begin: usize,
    pub row_end: usize,
    pub matrix: CsrMatrix,
}

/// Triplet assembly routed to row owners; duplicates are summed in (row, column, source rank,
/// insertion) order so every rank count gives bit-identical owned blocks.
#[derive(Clone, Debug)]
pub struct DistributedCsrAssembler {
    rows: usize,
    columns: usize,
    row_boundaries: Vec<usize>,
    entries: Vec<DistributedTriplet>,
}

impl DistributedCsrAssembler {
    pub fn new(rows: usize, columns: usize, row_boundaries: Vec<usize>) -> Result<Self> {
        require(
            row_boundaries.len() >= 2
                && row_boundaries[0] == 0
                && row_boundaries.last() == Some(&rows),
            "distributed assembly boundaries must cover all rows",
        )?;
        require(
            row_boundaries.windows(2).all(|w| w[1] > w[0]),
            "distributed assembly partitions must be non-empty",
        )?;
        Ok(Self {
            rows,
            columns,
            row_boundaries,
            entries: Vec::new(),
        })
    }
    /// Upstream default `source_rank = 0`.
    pub fn add(&mut self, row: usize, column: usize, value: f64, source_rank: usize) -> Result<()> {
        require(
            row < self.rows && column < self.columns,
            "distributed assembly index out of range",
        )?;
        require(
            value.is_finite(),
            "distributed assembly value must be finite",
        )?;
        let sequence = self.entries.len();
        self.entries.push(DistributedTriplet {
            row,
            column,
            value,
            source_rank,
            sequence,
        });
        Ok(())
    }
    pub fn parts(&self) -> usize {
        self.row_boundaries.len() - 1
    }
    pub fn owner_of_row(&self, row: usize) -> Result<usize> {
        require(row < self.rows, "distributed assembly row out of range")?;
        Ok(self.row_boundaries.partition_point(|b| *b <= row) - 1)
    }
    fn sort(entries: &mut [DistributedTriplet]) {
        entries.sort_by(|l, r| {
            (l.row, l.column, l.source_rank, l.sequence).cmp(&(
                r.row,
                r.column,
                r.source_rank,
                r.sequence,
            ))
        });
    }
    pub fn routed_entries(&self) -> Result<Vec<Vec<DistributedTriplet>>> {
        let mut routed = vec![Vec::new(); self.parts()];
        for e in &self.entries {
            routed[self.owner_of_row(e.row)?].push(*e);
        }
        routed.iter_mut().for_each(|b| Self::sort(b));
        Ok(routed)
    }
    fn build(
        rows: usize,
        columns: usize,
        entries: &[DistributedTriplet],
        begin: usize,
    ) -> Result<CsrMatrix> {
        let mut offsets = vec![0; rows + 1];
        let (mut cols, mut vals) = (
            Vec::with_capacity(entries.len()),
            Vec::with_capacity(entries.len()),
        );
        let mut cursor = 0;
        for row in 0..rows {
            let global = begin + row;
            while cursor < entries.len() && entries[cursor].row < global {
                cursor += 1;
            }
            while cursor < entries.len() && entries[cursor].row == global {
                let column = entries[cursor].column;
                let mut sum = 0.0;
                while cursor < entries.len()
                    && entries[cursor].row == global
                    && entries[cursor].column == column
                {
                    sum += entries[cursor].value;
                    cursor += 1;
                }
                cols.push(column);
                vals.push(sum);
            }
            offsets[row + 1] = vals.len();
        }
        require(
            cursor == entries.len(),
            "distributed assembly entry lies outside owned rows",
        )?;
        CsrMatrix::new(rows, columns, offsets, cols, vals)
    }
    pub fn build_owned(&self, rank: usize) -> Result<DistributedCsrBlock> {
        require(
            rank < self.parts(),
            "distributed assembly rank out of range",
        )?;
        let routed = self.routed_entries()?;
        let (begin, end) = (self.row_boundaries[rank], self.row_boundaries[rank + 1]);
        Ok(DistributedCsrBlock {
            row_begin: begin,
            row_end: end,
            matrix: Self::build(end - begin, self.columns, &routed[rank], begin)?,
        })
    }
    pub fn build_global_reference(&self) -> Result<CsrMatrix> {
        let mut e = self.entries.clone();
        Self::sort(&mut e);
        Self::build(self.rows, self.columns, &e, 0)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RefinementResult {
    pub converged: bool,
    pub outer_iterations: usize,
    pub relative_residual: f64,
}

/// Double-precision residual with f32 Jacobi corrections (upstream defaults 12, 100, 1e-11).
pub fn mixed_precision_iterative_refinement(
    a: &CsrMatrix,
    b: &[f64],
    x: &mut [f64],
    outer: usize,
    inner: usize,
    tolerance: f64,
) -> Result<RefinementResult> {
    require(
        a.rows() == a.cols()
            && b.len() == a.rows()
            && x.len() == a.cols()
            && outer > 0
            && inner > 0
            && tolerance > 0.0,
        "invalid iterative-refinement inputs",
    )?;
    let diag = a.diagonal();
    require(
        diag.iter().all(|d| *d != 0.0),
        "zero diagonal in iterative refinement",
    )?;
    let d: Vec<f32> = diag.iter().map(|v| *v as f32).collect();
    let n = b.len();
    let (mut ax, mut r) = (vec![0.0; n], vec![0.0; n]);
    let (mut rf, mut delta, mut next) = (vec![0.0f32; n], vec![0.0f32; n], vec![0.0f32; n]);
    let mut bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    if bnorm == 0.0 {
        bnorm = 1.0;
    }
    for out in 0..outer {
        a.multiply(x, &mut ax)?;
        let mut rn = 0.0;
        for i in 0..n {
            r[i] = b[i] - ax[i];
            rf[i] = r[i] as f32;
            rn += r[i] * r[i];
            delta[i] = 0.0;
        }
        let rn = rn.sqrt() / bnorm;
        if rn <= tolerance {
            return Ok(RefinementResult {
                converged: true,
                outer_iterations: out,
                relative_residual: rn,
            });
        }
        for _ in 0..inner {
            for row in 0..a.rows() {
                let mut sum = rf[row];
                for k in a.row_offsets()[row]..a.row_offsets()[row + 1] {
                    let col = a.columns()[k];
                    if col != row {
                        sum -= a.values()[k] as f32 * delta[col];
                    }
                }
                next[row] = sum / d[row];
            }
            std::mem::swap(&mut delta, &mut next);
        }
        for i in 0..n {
            x[i] += f64::from(delta[i]);
        }
    }
    a.multiply(x, &mut ax)?;
    let rn = (0..n)
        .map(|i| (b[i] - ax[i]) * (b[i] - ax[i]))
        .sum::<f64>()
        .sqrt()
        / bnorm;
    Ok(RefinementResult {
        converged: rn <= tolerance,
        outer_iterations: outer,
        relative_residual: rn,
    })
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProfileCounter {
    pub name: String,
    pub calls: usize,
    pub total_seconds: f64,
    pub minimum_seconds: f64,
    pub maximum_seconds: f64,
}

/// Thread-safe named timing counters.
#[derive(Debug)]
pub struct Profiler {
    counters: Mutex<HashMap<String, ProfileCounter>>,
    enabled: std::sync::atomic::AtomicBool,
}

impl Default for Profiler {
    fn default() -> Self {
        Self {
            counters: Mutex::new(HashMap::new()),
            enabled: std::sync::atomic::AtomicBool::new(true),
        }
    }
}

impl Profiler {
    pub fn global() -> &'static Profiler {
        static P: OnceLock<Profiler> = OnceLock::new();
        P.get_or_init(Profiler::default)
    }
    pub fn set_enabled(&self, enabled: bool) {
        self.enabled
            .store(enabled, std::sync::atomic::Ordering::Relaxed);
    }
    pub fn enabled(&self) -> bool {
        self.enabled.load(std::sync::atomic::Ordering::Relaxed)
    }
    pub fn record(&self, name: &str, seconds: f64) {
        if !self.enabled() {
            return;
        }
        let mut map = self.counters.lock().unwrap_or_else(|p| p.into_inner());
        let c = map.entry(name.to_string()).or_default();
        if c.calls == 0 {
            c.name = name.to_string();
            c.minimum_seconds = seconds;
            c.maximum_seconds = seconds;
        } else {
            c.minimum_seconds = c.minimum_seconds.min(seconds);
            c.maximum_seconds = c.maximum_seconds.max(seconds);
        }
        c.calls += 1;
        c.total_seconds += seconds;
    }
    /// Counters sorted by name.
    pub fn snapshot(&self) -> Vec<ProfileCounter> {
        let map = self.counters.lock().unwrap_or_else(|p| p.into_inner());
        let mut out: Vec<ProfileCounter> = map.values().cloned().collect();
        out.sort_by(|a, b| a.name.cmp(&b.name));
        out
    }
    pub fn reset(&self) {
        self.counters
            .lock()
            .unwrap_or_else(|p| p.into_inner())
            .clear();
    }
}

/// Records its lifetime into a profiler on drop.
pub struct ScopedProfile<'a> {
    name: String,
    profiler: &'a Profiler,
    start: Instant,
}

impl<'a> ScopedProfile<'a> {
    pub fn new(name: &str, profiler: &'a Profiler) -> Self {
        Self {
            name: name.to_string(),
            profiler,
            start: Instant::now(),
        }
    }
}

impl Drop for ScopedProfile<'_> {
    fn drop(&mut self) {
        self.profiler
            .record(&self.name, self.start.elapsed().as_secs_f64());
    }
}

/// Bump allocator of default-initialised `T` slots. Allocations return index ranges (handles)
/// rather than borrowed spans so several can be live at once without `unsafe`.
#[derive(Clone, Debug, Default)]
pub struct ScratchArena<T: Clone + Default> {
    storage: Vec<T>,
    offset: usize,
}

impl<T: Clone + Default> ScratchArena<T> {
    pub fn new(capacity: usize) -> Self {
        Self {
            storage: vec![T::default(); capacity],
            offset: 0,
        }
    }
    pub fn reserve(&mut self, capacity: usize) -> Result<()> {
        require(self.offset == 0, "cannot reserve active scratch arena")?;
        self.storage.resize(capacity, T::default());
        Ok(())
    }
    pub fn reset(&mut self) {
        self.offset = 0;
    }
    pub fn capacity(&self) -> usize {
        self.storage.len()
    }
    pub fn used(&self) -> usize {
        self.offset
    }
    pub fn allocate(&mut self, count: usize) -> Result<std::ops::Range<usize>> {
        let end = self
            .offset
            .checked_add(count)
            .filter(|e| *e <= self.storage.len())
            .ok_or_else(|| SolverError::InvalidInput("scratch arena exhausted".into()))?;
        let range = self.offset..end;
        self.storage[range.clone()].fill(T::default());
        self.offset = end;
        Ok(range)
    }
    pub fn slice(&self, range: std::ops::Range<usize>) -> &[T] {
        &self.storage[range]
    }
    pub fn slice_mut(&mut self, range: std::ops::Range<usize>) -> &mut [T] {
        &mut self.storage[range]
    }
}

/// Fixed-width lane pack with element-wise arithmetic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimdPack<T: Copy + Default, const N: usize> {
    pub lane: [T; N],
}

impl<T: Copy + Default, const N: usize> Default for SimdPack<T, N> {
    fn default() -> Self {
        Self {
            lane: [T::default(); N],
        }
    }
}

impl<T: Copy + Default + AddAssign, const N: usize> SimdPack<T, N> {
    pub fn splat(x: T) -> Self {
        Self { lane: [x; N] }
    }
    pub fn sum(&self) -> T {
        let mut r = T::default();
        for x in self.lane {
            r += x;
        }
        r
    }
}

impl<T: Copy + Default + Add<Output = T>, const N: usize> Add for SimdPack<T, N> {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            lane: std::array::from_fn(|i| self.lane[i] + b.lane[i]),
        }
    }
}
impl<T: Copy + Default + Sub<Output = T>, const N: usize> Sub for SimdPack<T, N> {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            lane: std::array::from_fn(|i| self.lane[i] - b.lane[i]),
        }
    }
}
impl<T: Copy + Default + Mul<Output = T>, const N: usize> Mul for SimdPack<T, N> {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            lane: std::array::from_fn(|i| self.lane[i] * b.lane[i]),
        }
    }
}
impl<T: Copy + Default + Mul<Output = T>, const N: usize> Mul<T> for SimdPack<T, N> {
    type Output = Self;
    fn mul(self, b: T) -> Self {
        Self {
            lane: std::array::from_fn(|i| self.lane[i] * b),
        }
    }
}

/// Structure-of-arrays storage with `C` contiguous components.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StaticSoA<T: Clone + Default, const C: usize> {
    size: usize,
    data: Vec<T>,
}

impl<T: Clone + Default, const C: usize> StaticSoA<T, C> {
    pub fn new(size: usize) -> Self {
        Self {
            size,
            data: vec![T::default(); size * C],
        }
    }
    pub fn resize(&mut self, size: usize) {
        self.size = size;
        self.data.resize(size * C, T::default());
    }
    pub fn size(&self) -> usize {
        self.size
    }
    pub const fn components() -> usize {
        C
    }
    pub fn component(&self, c: usize) -> Result<&[T]> {
        require(c < C, "StaticSoA component")?;
        Ok(&self.data[c * self.size..(c + 1) * self.size])
    }
    pub fn component_mut(&mut self, c: usize) -> Result<&mut [T]> {
        require(c < C, "StaticSoA component")?;
        Ok(&mut self.data[c * self.size..(c + 1) * self.size])
    }
    pub fn raw(&self) -> &[T] {
        &self.data
    }
}

impl<T: Clone + Default, const C: usize> Index<(usize, usize)> for StaticSoA<T, C> {
    type Output = T;
    fn index(&self, (c, i): (usize, usize)) -> &T {
        &self.data[c * self.size + i]
    }
}
impl<T: Clone + Default, const C: usize> IndexMut<(usize, usize)> for StaticSoA<T, C> {
    fn index_mut(&mut self, (c, i): (usize, usize)) -> &mut T {
        &mut self.data[c * self.size + i]
    }
}

/// `Vec<T>`; Rust's safe allocator API cannot request 64-byte alignment without `unsafe`.
pub type AlignedVector<T> = Vec<T>;

pub type Task = Box<dyn Fn() -> Result<()> + Send + Sync>;

/// Dependency graph of tasks run in topological order, optionally level-by-level on threads.
#[derive(Default)]
pub struct TaskGraph {
    nodes: Vec<(String, Task, Vec<usize>)>,
}

impl TaskGraph {
    pub fn add(&mut self, name: &str, task: Task, dependencies: Vec<usize>) -> Result<usize> {
        require(
            dependencies.iter().all(|d| *d < self.nodes.len()),
            "task dependency does not exist",
        )?;
        self.nodes.push((name.to_string(), task, dependencies));
        Ok(self.nodes.len() - 1)
    }
    pub fn size(&self) -> usize {
        self.nodes.len()
    }
    /// Depth-first post-order, as upstream. Cycles cannot be formed because dependencies must
    /// already exist, but the check is kept.
    pub fn topological_order(&self) -> Result<Vec<usize>> {
        let mut state = vec![0u8; self.nodes.len()];
        let mut out = Vec::with_capacity(self.nodes.len());
        for root in 0..self.nodes.len() {
            let mut stack = vec![(root, 0usize)];
            while let Some((n, k)) = stack.pop() {
                if k == 0 {
                    if state[n] == 2 {
                        continue;
                    }
                    require(state[n] != 1, "task graph cycle")?;
                    state[n] = 1;
                }
                if let Some(&d) = self.nodes[n].2.get(k) {
                    stack.push((n, k + 1));
                    if state[d] != 2 {
                        require(state[d] != 1, "task graph cycle")?;
                        stack.push((d, 0));
                    }
                } else {
                    state[n] = 2;
                    out.push(n);
                }
            }
        }
        Ok(out)
    }
    /// Runs every task once after its dependencies; `asynchronous` runs independent tasks of the
    /// same dependency level on scoped threads. The first task error is returned.
    pub fn run(&self, asynchronous: bool) -> Result<()> {
        let order = self.topological_order()?;
        if !asynchronous {
            for id in order {
                (self.nodes[id].1)()?;
            }
            return Ok(());
        }
        let mut level = vec![0usize; self.nodes.len()];
        for &id in &order {
            level[id] = self.nodes[id]
                .2
                .iter()
                .map(|d| level[*d] + 1)
                .max()
                .unwrap_or(0);
        }
        let depth = level.iter().copied().max().map_or(0, |m| m + 1);
        for l in 0..depth {
            let ids: Vec<usize> = order.iter().copied().filter(|id| level[*id] == l).collect();
            let results: Vec<Result<()>> = std::thread::scope(|s| {
                let handles: Vec<_> = ids
                    .iter()
                    .map(|id| s.spawn(|| (self.nodes[*id].1)()))
                    .collect();
                handles
                    .into_iter()
                    .map(|h| {
                        h.join().unwrap_or_else(|_| {
                            Err(SolverError::NumericalFailure("task panicked".into()))
                        })
                    })
                    .collect()
            });
            for r in results {
                r?;
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeshPartitionRange {
    pub element_begin: usize,
    pub element_end: usize,
    pub owned_nodes: Vec<usize>,
    pub halo_nodes: Vec<usize>,
    pub shared_owned_nodes: Vec<usize>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeshPartitionPlan {
    pub global_node_count: usize,
    pub element_boundaries: Vec<usize>,
    /// `None` for nodes not used by any element.
    pub node_owner: Vec<Option<usize>>,
    pub partitions: Vec<MeshPartitionRange>,
}

impl MeshPartitionPlan {
    pub fn parts(&self) -> usize {
        self.partitions.len()
    }
    pub fn owner_of_node(&self, node: usize) -> Result<usize> {
        self.node_owner
            .get(node)
            .copied()
            .flatten()
            .ok_or_else(|| SolverError::InvalidInput("mesh node has no deterministic owner".into()))
    }
}

/// Weighted contiguous element partition; nodes are owned by the first part that uses them.
pub fn partition_tri3_mesh(
    elements: &[[usize; 3]],
    node_count: usize,
    parts: usize,
    element_weights: &[f64],
) -> Result<MeshPartitionPlan> {
    require(
        !elements.is_empty() && node_count > 0 && parts > 0 && elements.len() >= parts,
        "invalid Tri3 mesh partition dimensions",
    )?;
    require(
        element_weights.is_empty() || element_weights.len() == elements.len(),
        "mesh partition weight count mismatch",
    )?;
    require(
        elements.iter().flatten().all(|n| *n < node_count),
        "mesh partition node index out of range",
    )?;
    let weights = if element_weights.is_empty() {
        vec![1.0; elements.len()]
    } else {
        element_weights.to_vec()
    };
    let weighted = crate::distributed::plan_contiguous_repartition(&weights, parts, 1)?;
    let mut result = MeshPartitionPlan {
        global_node_count: node_count,
        element_boundaries: weighted.boundaries,
        node_owner: vec![None; node_count],
        partitions: vec![MeshPartitionRange::default(); parts],
    };
    let mut used_nodes = vec![Vec::new(); parts];
    for part in 0..parts {
        let (b, e) = (
            result.element_boundaries[part],
            result.element_boundaries[part + 1],
        );
        let range = &mut result.partitions[part];
        range.element_begin = b;
        range.element_end = e;
        let used = &mut used_nodes[part];
        for element in &elements[b..e] {
            for &node in element {
                if result.node_owner[node].is_none() {
                    result.node_owner[node] = Some(part);
                }
                used.push(node);
            }
        }
        used.sort_unstable();
        used.dedup();
        for &node in used.iter() {
            if result.node_owner[node] == Some(part) {
                range.owned_nodes.push(node);
            } else {
                range.halo_nodes.push(node);
            }
        }
    }
    for part in 0..parts {
        let shared: Vec<usize> = result.partitions[part]
            .owned_nodes
            .iter()
            .copied()
            .filter(|n| {
                used_nodes
                    .iter()
                    .enumerate()
                    .any(|(p, peer)| p != part && peer.binary_search(n).is_ok())
            })
            .collect();
        result.partitions[part].shared_owned_nodes = shared;
    }
    Ok(result)
}

pub type ReverseObjective = Box<dyn Fn(&ReverseTape, &[ReverseScalar]) -> Result<ReverseScalar>>;

pub struct ReverseAdPartition {
    pub global_indices: Vec<usize>,
    pub objective: ReverseObjective,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReverseAdDistributedResult {
    pub objective: f64,
    pub partitions: usize,
}

/// Sum of per-partition reverse-mode objectives with gradients scattered to global indices.
pub struct DistributedReverseOperator {
    global_size: usize,
    partitions: Vec<ReverseAdPartition>,
}

impl DistributedReverseOperator {
    pub fn new(global_size: usize) -> Result<Self> {
        require(global_size > 0, "reverse AD global size must be non-zero")?;
        Ok(Self {
            global_size,
            partitions: Vec::new(),
        })
    }
    pub fn add_partition(&mut self, partition: ReverseAdPartition) -> Result<()> {
        require(
            !partition.global_indices.is_empty(),
            "invalid reverse AD partition",
        )?;
        require(
            partition
                .global_indices
                .iter()
                .all(|i| *i < self.global_size),
            "reverse AD global index",
        )?;
        self.partitions.push(partition);
        Ok(())
    }
    pub fn evaluate(
        &self,
        global_state: &[f64],
        global_gradient: &mut [f64],
    ) -> Result<ReverseAdDistributedResult> {
        require(
            global_state.len() == self.global_size
                && global_gradient.len() == self.global_size
                && !self.partitions.is_empty(),
            "reverse AD distributed size mismatch",
        )?;
        global_gradient.fill(0.0);
        let mut objective = 0.0;
        for p in &self.partitions {
            let tape = ReverseTape::new();
            let vars: Vec<ReverseScalar> = p
                .global_indices
                .iter()
                .map(|i| tape.variable(global_state[*i]))
                .collect();
            let out = (p.objective)(&tape, &vars)?;
            let local = tape.backward(&out, &vars)?;
            objective += out.value();
            for (g, i) in local.iter().zip(&p.global_indices) {
                global_gradient[*i] += g;
            }
        }
        require(
            objective.is_finite() && global_gradient.iter().all(|v| v.is_finite()),
            "non-finite distributed reverse AD result",
        )?;
        Ok(ReverseAdDistributedResult {
            objective,
            partitions: self.partitions.len(),
        })
    }
}

pub type AmgApply = Box<dyn Fn(&[f64], &mut [f64]) -> Result<()>>;

/// Wraps an external AMG V-cycle as a [`Preconditioner`].
pub struct AmgPreconditionerAdapter {
    apply: AmgApply,
}

impl AmgPreconditionerAdapter {
    pub fn new(apply: AmgApply) -> Self {
        Self { apply }
    }
}

impl Preconditioner for AmgPreconditionerAdapter {
    fn apply(&self, rhs: &[f64], result: &mut [f64]) -> Result<()> {
        require(rhs.len() == result.len(), "AMG adapter size mismatch")?;
        (self.apply)(rhs, result)
    }
}
