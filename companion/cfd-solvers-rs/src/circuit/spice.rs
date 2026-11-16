//! Modified-nodal-analysis circuit kernel: linear elements, controlled sources, diode, MOS
//! level 1, Ebers-Moll BJT, JFET, hysteretic switch, frequency-domain N-ports and transmission
//! lines, static / charge-based dynamic / electrothermal compact devices, behavioral sources,
//! DC Newton with source and gmin homotopy, AC, sweeps, fixed-step BE/trapezoidal/BDF2 and
//! step-doubling adaptive transient, and small-signal output noise. Ported from upstream
//! `cfd/circuit/spice` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::analysis::limit_pn_junction_voltage;
use crate::core::{restarted_gmres, CsrBuilder, Ilu0Preconditioner, LinearConfig};
use crate::rf::{self, ComplexMatrix, NPortPoint, TemTransmissionLine, TransmissionLine};
pub use crate::tcad::{StaticDeviceEvaluation, StaticDeviceEvaluator};
use crate::{require, Result, SolverError};
use num_complex::Complex64;
use std::{
    collections::HashMap,
    f64::consts::PI,
    ops::{AddAssign, Mul, Neg, SubAssign},
    sync::Arc,
};

pub type Complex = Complex64;

pub(crate) const BOLTZMANN: f64 = 1.380649e-23;
pub(crate) const ELECTRON_CHARGE: f64 = 1.602176634e-19;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiodeModel {
    pub saturation_current: f64,
    pub emission_coefficient: f64,
    pub temperature_k: f64,
}
impl Default for DiodeModel {
    fn default() -> Self {
        Self {
            saturation_current: 1.0e-14,
            emission_coefficient: 1.0,
            temperature_k: 300.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MosLevel1Model {
    pub pmos: bool,
    pub threshold_voltage: f64,
    /// A/V².
    pub transconductance_parameter: f64,
    /// 1/V.
    pub channel_length_modulation: f64,
}
impl Default for MosLevel1Model {
    fn default() -> Self {
        Self {
            pmos: false,
            threshold_voltage: 0.7,
            transconductance_parameter: 1.0e-3,
            channel_length_modulation: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BjtModel {
    pub pnp: bool,
    pub saturation_current: f64,
    pub forward_beta: f64,
    pub reverse_beta: f64,
    pub temperature_k: f64,
}
impl Default for BjtModel {
    fn default() -> Self {
        Self {
            pnp: false,
            saturation_current: 1.0e-15,
            forward_beta: 100.0,
            reverse_beta: 1.0,
            temperature_k: 300.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JfetModel {
    pub pchannel: bool,
    /// V; the n-channel convention is negative.
    pub pinch_off_voltage: f64,
    /// A/V².
    pub beta: f64,
    pub channel_length_modulation: f64,
}
impl Default for JfetModel {
    fn default() -> Self {
        Self {
            pchannel: false,
            pinch_off_voltage: -2.0,
            beta: 1.0e-3,
            channel_length_modulation: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SwitchModel {
    pub on_resistance: f64,
    pub off_resistance: f64,
    pub threshold_voltage: f64,
    pub hysteresis_voltage: f64,
}
impl Default for SwitchModel {
    fn default() -> Self {
        Self {
            on_resistance: 1.0,
            off_resistance: 1.0e12,
            threshold_voltage: 0.0,
            hysteresis_voltage: 0.0,
        }
    }
}

/// Newton and linear-solver controls. Real DC/transient systems with at least
/// `sparse_mna_threshold` unknowns use CSR + ILU(0)-GMRES (0 forces the dense path).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NewtonConfig {
    pub max_iterations: usize,
    pub voltage_tolerance: f64,
    pub residual_tolerance: f64,
    pub gmin: f64,
    pub sparse_mna_threshold: usize,
    pub sparse_max_iterations: usize,
    pub sparse_gmres_restart: usize,
    pub sparse_relative_tolerance: f64,
    pub sparse_dense_fallback: bool,
}
impl Default for NewtonConfig {
    fn default() -> Self {
        Self {
            max_iterations: 80,
            voltage_tolerance: 1.0e-10,
            residual_tolerance: 1.0e-10,
            gmin: 1.0e-12,
            sparse_mna_threshold: 64,
            sparse_max_iterations: 600,
            sparse_gmres_restart: 40,
            sparse_relative_tolerance: 1.0e-11,
            sparse_dense_fallback: true,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OperatingPoint {
    /// Includes node 0 (ground).
    pub node_voltage: Vec<f64>,
    /// Branch order: independent V, behavioral V, VCVS, CCVS, then inductors.
    pub branch_current: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AcSolution {
    pub frequency_hz: f64,
    pub node_voltage: Vec<Complex>,
    pub branch_current: Vec<Complex>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransientPoint {
    pub time: f64,
    pub node_voltage: Vec<f64>,
    pub branch_current: Vec<f64>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum TransientMethod {
    #[default]
    BackwardEuler,
    Trapezoidal,
    Bdf2,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveTransientConfig {
    pub initial_step: f64,
    pub minimum_step: f64,
    pub maximum_step: f64,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
    pub safety_factor: f64,
    pub minimum_scale: f64,
    pub maximum_scale: f64,
    pub max_step_attempts: usize,
    /// Backward Euler or BDF2 (trapezoidal is rejected, as upstream).
    pub method: TransientMethod,
}
impl Default for AdaptiveTransientConfig {
    fn default() -> Self {
        Self {
            initial_step: 1.0e-6,
            minimum_step: 1.0e-12,
            maximum_step: 1.0,
            relative_tolerance: 1.0e-4,
            absolute_tolerance: 1.0e-7,
            safety_factor: 0.9,
            minimum_scale: 0.2,
            maximum_scale: 5.0,
            max_step_attempts: 100_000,
            method: TransientMethod::BackwardEuler,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AdaptiveTransientResult {
    pub points: Vec<TransientPoint>,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub minimum_accepted_step: f64,
    pub maximum_accepted_step: f64,
}

/// Time-dependent source value.
pub type SourceWaveform = Arc<dyn Fn(f64) -> Result<f64> + Send + Sync>;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PulseWaveform {
    pub initial: f64,
    pub pulsed: f64,
    pub delay: f64,
    pub rise: f64,
    pub fall: f64,
    pub width: f64,
    pub period: f64,
}

impl PulseWaveform {
    pub fn value(&self, time: f64) -> Result<f64> {
        require(
            time.is_finite()
                && self.initial.is_finite()
                && self.pulsed.is_finite()
                && self.delay >= 0.0
                && self.rise >= 0.0
                && self.fall >= 0.0
                && self.width >= 0.0
                && self.period >= 0.0,
            "invalid PULSE waveform",
        )?;
        if time < self.delay {
            return Ok(self.initial);
        }
        let mut local = time - self.delay;
        if self.period > 0.0 {
            local %= self.period;
        }
        if self.rise > 0.0 && local < self.rise {
            return Ok(self.initial + (self.pulsed - self.initial) * (local / self.rise));
        }
        local -= self.rise;
        if local < self.width {
            return Ok(self.pulsed);
        }
        local -= self.width;
        if self.fall > 0.0 && local < self.fall {
            return Ok(self.pulsed + (self.initial - self.pulsed) * (local / self.fall));
        }
        Ok(self.initial)
    }
    pub fn waveform(self) -> SourceWaveform {
        Arc::new(move |t| self.value(t))
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SineWaveform {
    pub offset: f64,
    pub amplitude: f64,
    pub frequency_hz: f64,
    pub delay: f64,
    pub damping_per_s: f64,
    pub phase_deg: f64,
}

impl SineWaveform {
    pub fn value(&self, time: f64) -> Result<f64> {
        require(
            time.is_finite()
                && self.offset.is_finite()
                && self.amplitude.is_finite()
                && self.frequency_hz >= 0.0
                && self.delay >= 0.0
                && self.damping_per_s >= 0.0
                && self.phase_deg.is_finite(),
            "invalid SIN waveform",
        )?;
        if time < self.delay {
            return Ok(self.offset);
        }
        let local = time - self.delay;
        let phase = self.phase_deg * PI / 180.0;
        Ok(self.offset
            + self.amplitude
                * (-self.damping_per_s * local).exp()
                * (2.0 * PI * self.frequency_hz * local + phase).sin())
    }
    pub fn waveform(self) -> SourceWaveform {
        Arc::new(move |t| self.value(t))
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PwlWaveform {
    pub points: Vec<(f64, f64)>,
}

impl PwlWaveform {
    pub fn value(&self, time: f64) -> Result<f64> {
        let (Some(first), Some(last)) = (self.points.first(), self.points.last()) else {
            return Err(SolverError::InvalidInput("invalid PWL waveform".into()));
        };
        require(time.is_finite(), "invalid PWL waveform")?;
        if time <= first.0 {
            return Ok(first.1);
        }
        for w in self.points.windows(2) {
            let ((t0, v0), (t1, v1)) = (w[0], w[1]);
            require(
                t1 > t0 && t0.is_finite() && t1.is_finite() && v0.is_finite() && v1.is_finite(),
                "PWL times must be finite and strictly increasing",
            )?;
            if time <= t1 {
                return Ok(v0 + (v1 - v0) * ((time - t0) / (t1 - t0)));
            }
        }
        Ok(last.1)
    }
    pub fn waveform(self) -> SourceWaveform {
        Arc::new(move |t| self.value(t))
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DcSweepPoint {
    pub source_value: f64,
    pub operating_point: OperatingPoint,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TemperatureSweepPoint {
    pub temperature_k: f64,
    pub operating_point: OperatingPoint,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct NoiseContribution {
    pub source: String,
    pub output_noise_density_v2_per_hz: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct NoiseResult {
    pub frequency_hz: f64,
    pub output_noise_v_per_sqrt_hz: f64,
    pub output_noise_density_v2_per_hz: f64,
    pub contributions: Vec<NoiseContribution>,
}

/// DAE form per terminal `I(V) + dQ(V)/dt = 0`; currents and charges are positive out of each
/// terminal; Jacobians are row-major.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DynamicDeviceEvaluation {
    pub terminal_current: Vec<f64>,
    pub current_jacobian: Vec<f64>,
    pub terminal_charge: Vec<f64>,
    pub charge_jacobian: Vec<f64>,
}
pub type DynamicDeviceEvaluator =
    Arc<dyn Fn(&[f64]) -> Result<DynamicDeviceEvaluation> + Send + Sync>;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ElectroThermalDeviceEvaluation {
    pub electrical_current: Vec<f64>,
    /// dI(row)/dV(column).
    pub electrical_jacobian: Vec<f64>,
    /// dI/dT.
    pub electrical_temperature_derivative: Vec<f64>,
    pub dissipated_power_w: f64,
    /// dP/dV(column).
    pub power_voltage_derivative: Vec<f64>,
    pub power_temperature_derivative_w_per_k: f64,
}
/// Evaluated with the electrical terminal voltages and the device temperature in K.
pub type ElectroThermalDeviceEvaluator =
    Arc<dyn Fn(&[f64], f64) -> Result<ElectroThermalDeviceEvaluation> + Send + Sync>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectroThermalDeviceConfig {
    pub ambient_temperature_k: f64,
    pub thermal_resistance_k_per_w: f64,
    pub thermal_capacitance_j_per_k: f64,
}
impl Default for ElectroThermalDeviceConfig {
    fn default() -> Self {
        Self {
            ambient_temperature_k: 300.0,
            thermal_resistance_k_per_w: 1.0,
            thermal_capacitance_j_per_k: 0.0,
        }
    }
}

/// SPICE B-source: value from all node voltages (index 0 = ground) and time.
pub type BehavioralEvaluator = Arc<dyn Fn(&[f64], f64) -> Result<f64> + Send + Sync>;

#[derive(Clone, Debug)]
pub(crate) struct TwoTerminal {
    pub name: String,
    pub p: usize,
    pub n: usize,
    pub value: f64,
}
#[derive(Clone, Debug)]
pub(crate) struct MutualInductance {
    pub first: String,
    pub second: String,
    pub coupling: f64,
}
#[derive(Clone)]
pub(crate) struct IndependentSource<T> {
    pub name: String,
    pub p: usize,
    pub n: usize,
    pub dc: f64,
    pub ac: T,
    pub transient: Option<SourceWaveform>,
}
#[derive(Clone, Debug)]
pub(crate) struct VoltageControlled {
    pub name: String,
    pub p: usize,
    pub n: usize,
    pub cp: usize,
    pub cn: usize,
    pub gain: f64,
}
#[derive(Clone, Debug)]
pub(crate) struct CurrentControlled {
    pub name: String,
    pub p: usize,
    pub n: usize,
    pub control: String,
    pub gain: f64,
}
#[derive(Clone, Debug)]
pub(crate) struct ModelDevice {
    pub name: String,
    pub nodes: [usize; 4],
    pub model: String,
}
#[derive(Clone, Debug)]
pub(crate) struct NPortDevice {
    pub ports: Vec<(usize, usize)>,
    pub samples: Vec<NPortPoint>,
}
#[derive(Clone)]
pub(crate) struct StaticDevice {
    pub terminals: Vec<usize>,
    pub evaluator: StaticDeviceEvaluator,
}
#[derive(Clone)]
pub(crate) struct DynamicDevice {
    pub terminals: Vec<usize>,
    pub evaluator: DynamicDeviceEvaluator,
}
#[derive(Clone)]
pub(crate) struct BehavioralSource {
    pub name: String,
    pub p: usize,
    pub n: usize,
    pub evaluator: BehavioralEvaluator,
}

/// A netlist. Nodes are numbered from 0 (ground, also reachable as `GND`); node names are
/// case-insensitive.
#[derive(Clone)]
pub struct Circuit {
    pub(crate) node_names: Vec<String>,
    pub(crate) node_lookup: HashMap<String, usize>,
    pub(crate) resistors: Vec<TwoTerminal>,
    pub(crate) capacitors: Vec<TwoTerminal>,
    pub(crate) inductors: Vec<TwoTerminal>,
    pub(crate) mutual_inductances: Vec<MutualInductance>,
    pub(crate) voltage_sources: Vec<IndependentSource<Complex>>,
    pub(crate) current_sources: Vec<IndependentSource<Complex>>,
    pub(crate) vccs: Vec<VoltageControlled>,
    pub(crate) vcvs: Vec<VoltageControlled>,
    pub(crate) cccs: Vec<CurrentControlled>,
    pub(crate) ccvs: Vec<CurrentControlled>,
    pub(crate) diodes: Vec<ModelDevice>,
    pub(crate) mosfets: Vec<ModelDevice>,
    pub(crate) bjts: Vec<ModelDevice>,
    pub(crate) jfets: Vec<ModelDevice>,
    pub(crate) switches: Vec<ModelDevice>,
    pub(crate) nports: Vec<NPortDevice>,
    pub(crate) static_devices: Vec<StaticDevice>,
    pub(crate) dynamic_devices: Vec<DynamicDevice>,
    pub(crate) behavioral_current_sources: Vec<BehavioralSource>,
    pub(crate) behavioral_voltage_sources: Vec<BehavioralSource>,
    pub(crate) diode_models: HashMap<String, DiodeModel>,
    pub(crate) mos_models: HashMap<String, MosLevel1Model>,
    pub(crate) bjt_models: HashMap<String, BjtModel>,
    pub(crate) jfet_models: HashMap<String, JfetModel>,
    pub(crate) switch_models: HashMap<String, SwitchModel>,
    pub(crate) initial_voltage: HashMap<usize, f64>,
    pub(crate) nodeset_voltage: HashMap<usize, f64>,
}

impl std::fmt::Debug for Circuit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Circuit")
            .field("nodes", &self.node_names)
            .field("resistors", &self.resistors.len())
            .field("capacitors", &self.capacitors.len())
            .field("inductors", &self.inductors.len())
            .field("voltage_sources", &self.voltage_sources.len())
            .field("current_sources", &self.current_sources.len())
            .finish_non_exhaustive()
    }
}

impl Default for Circuit {
    fn default() -> Self {
        Self::new()
    }
}

pub(crate) fn upper(s: &str) -> String {
    s.to_ascii_uppercase()
}

fn err_runtime(msg: impl Into<String>) -> SolverError {
    SolverError::NumericalFailure(msg.into())
}

// ----------------------------------------------------------------------------------------------
// Matrix stamping helpers
// ----------------------------------------------------------------------------------------------

pub(crate) trait MnaValue:
    Copy + Default + AddAssign + SubAssign + Neg<Output = Self> + Mul<Output = Self> + From<f64>
{
    fn magnitude(self) -> f64;
    fn divide(self, other: Self) -> Self;
}
impl MnaValue for f64 {
    fn magnitude(self) -> f64 {
        self.abs()
    }
    fn divide(self, other: Self) -> Self {
        self / other
    }
}
impl MnaValue for Complex {
    fn magnitude(self) -> f64 {
        self.norm()
    }
    fn divide(self, other: Self) -> Self {
        self / other
    }
}

pub(crate) trait MnaMatrix {
    type T: MnaValue;
    fn at(&mut self, index: usize) -> &mut Self::T;
}

impl<T: MnaValue> MnaMatrix for Vec<T> {
    type T = T;
    fn at(&mut self, index: usize) -> &mut T {
        &mut self[index]
    }
}

pub(crate) struct RealMna {
    dimension: usize,
    dense: Vec<f64>,
    sparse: Option<HashMap<usize, f64>>,
}

impl RealMna {
    fn new(dimension: usize, sparse: bool) -> Result<Self> {
        require(dimension > 0, "MNA matrix dimension must be positive")?;
        Ok(if sparse {
            Self {
                dimension,
                dense: Vec::new(),
                sparse: Some(HashMap::new()),
            }
        } else {
            Self {
                dimension,
                dense: vec![0.0; dimension * dimension],
                sparse: None,
            }
        })
    }
    fn dense_copy(&self) -> Vec<f64> {
        match &self.sparse {
            None => self.dense.clone(),
            Some(map) => {
                let mut out = vec![0.0; self.dimension * self.dimension];
                for (&i, &v) in map {
                    out[i] = v;
                }
                out
            }
        }
    }
}

impl MnaMatrix for RealMna {
    type T = f64;
    fn at(&mut self, index: usize) -> &mut f64 {
        match &mut self.sparse {
            None => &mut self.dense[index],
            Some(map) => map.entry(index).or_insert(0.0),
        }
    }
}

/// Gauss-Jordan elimination with partial pivoting (upstream `solve_dense`).
pub(crate) fn solve_dense<T: MnaValue>(
    mut matrix: Vec<T>,
    mut rhs: Vec<T>,
    n: usize,
    what: &str,
) -> Result<Vec<T>> {
    require(
        matrix.len() == n * n && rhs.len() == n,
        "dense circuit matrix size mismatch",
    )?;
    for column in 0..n {
        let mut pivot = column;
        let mut best = matrix[column * n + column].magnitude();
        for row in column + 1..n {
            let candidate = matrix[row * n + column].magnitude();
            if candidate > best {
                best = candidate;
                pivot = row;
            }
        }
        if best < 1.0e-24 || best.is_nan() {
            return Err(SolverError::Singular(what.into()));
        }
        if pivot != column {
            for j in column..n {
                matrix.swap(column * n + j, pivot * n + j);
            }
            rhs.swap(column, pivot);
        }
        let diagonal = matrix[column * n + column];
        for j in column..n {
            matrix[column * n + j] = matrix[column * n + j].divide(diagonal);
        }
        rhs[column] = rhs[column].divide(diagonal);
        for row in 0..n {
            if row == column {
                continue;
            }
            let factor = matrix[row * n + column];
            if factor.magnitude() == 0.0 {
                continue;
            }
            for j in column..n {
                let v = factor * matrix[column * n + j];
                matrix[row * n + j] -= v;
            }
            let v = factor * rhs[column];
            rhs[row] -= v;
        }
    }
    Ok(rhs)
}

fn solve_real_mna(matrix: &RealMna, rhs: Vec<f64>, config: &NewtonConfig) -> Result<Vec<f64>> {
    let n = matrix.dimension;
    require(rhs.len() == n, "MNA right-hand-side size mismatch")?;
    let Some(entries) = &matrix.sparse else {
        return solve_dense(matrix.dense_copy(), rhs, n, "singular circuit matrix");
    };
    require(
        config.sparse_max_iterations > 0
            && config.sparse_gmres_restart > 0
            && config.sparse_relative_tolerance > 0.0,
        "invalid sparse MNA controls",
    )?;
    let attempt = || -> Result<Option<Vec<f64>>> {
        let mut builder = CsrBuilder::new(n, n);
        // Keep a structural diagonal for ILU(0) even for ideal voltage-source branch equations.
        for row in 0..n {
            builder.add(row, row, f64::from_bits(1))?;
        }
        for (&index, &value) in entries {
            builder.add(index / n, index % n, value)?;
        }
        let csr = builder.build(0.0)?;
        let ilu = Ilu0Preconditioner::new(&csr, 1.0e-30)?;
        let mut solution = vec![0.0; n];
        let result = restarted_gmres(
            &rhs,
            &mut solution,
            |x, y| csr.multiply(x, y),
            &ilu,
            LinearConfig {
                max_iterations: config.sparse_max_iterations,
                relative_tolerance: config.sparse_relative_tolerance,
                restart: config.sparse_gmres_restart.min(n),
            },
        )?;
        Ok(result.converged.then_some(solution))
    };
    match attempt() {
        Ok(Some(solution)) => return Ok(solution),
        Ok(None) if !config.sparse_dense_fallback => {
            return Err(err_runtime("sparse MNA GMRES did not converge"));
        }
        Err(e) if !config.sparse_dense_fallback => return Err(e),
        _ => {}
    }
    solve_dense(matrix.dense_copy(), rhs, n, "singular circuit matrix")
}

fn stamp_admittance<M: MnaMatrix>(m: &mut M, dim: usize, p: usize, n: usize, y: M::T) {
    if p > 0 {
        *m.at((p - 1) * dim + (p - 1)) += y;
    }
    if n > 0 {
        *m.at((n - 1) * dim + (n - 1)) += y;
    }
    if p > 0 && n > 0 {
        *m.at((p - 1) * dim + (n - 1)) -= y;
        *m.at((n - 1) * dim + (p - 1)) -= y;
    }
}

fn stamp_current<T: MnaValue>(rhs: &mut [T], p: usize, n: usize, current: T) {
    if p > 0 {
        rhs[p - 1] -= current;
    }
    if n > 0 {
        rhs[n - 1] += current;
    }
}

#[allow(clippy::too_many_arguments)]
fn stamp_voltage_branch<M: MnaMatrix>(
    m: &mut M,
    rhs: &mut [M::T],
    dim: usize,
    node_count: usize,
    branch: usize,
    p: usize,
    n: usize,
    voltage: M::T,
    branch_diagonal: M::T,
) {
    let eq = (node_count - 1) + branch;
    let one = M::T::from(1.0);
    if p > 0 {
        *m.at((p - 1) * dim + eq) += one;
        *m.at(eq * dim + (p - 1)) += one;
    }
    if n > 0 {
        *m.at((n - 1) * dim + eq) -= one;
        *m.at(eq * dim + (n - 1)) -= one;
    }
    *m.at(eq * dim + eq) += branch_diagonal;
    rhs[eq] += voltage;
}

fn stamp_vccs<M: MnaMatrix>(m: &mut M, dim: usize, s: &VoltageControlled, gain: M::T) {
    let mut add = |row: usize, col: usize, v: M::T| {
        if row > 0 && col > 0 {
            *m.at((row - 1) * dim + (col - 1)) += v;
        }
    };
    add(s.p, s.cp, gain);
    add(s.p, s.cn, -gain);
    add(s.n, s.cp, -gain);
    add(s.n, s.cn, gain);
}

fn stamp_cccs<M: MnaMatrix>(
    m: &mut M,
    dim: usize,
    node_count: usize,
    p: usize,
    n: usize,
    control: usize,
    gain: M::T,
) {
    let col = (node_count - 1) + control;
    if p > 0 {
        *m.at((p - 1) * dim + col) += gain;
    }
    if n > 0 {
        *m.at((n - 1) * dim + col) -= gain;
    }
}

fn stamp_vcvs_control<M: MnaMatrix>(
    m: &mut M,
    dim: usize,
    node_count: usize,
    branch: usize,
    cp: usize,
    cn: usize,
    gain: M::T,
) {
    let row = (node_count - 1) + branch;
    if cp > 0 {
        *m.at(row * dim + (cp - 1)) -= gain;
    }
    if cn > 0 {
        *m.at(row * dim + (cn - 1)) += gain;
    }
}

fn stamp_ccvs_control<M: MnaMatrix>(
    m: &mut M,
    dim: usize,
    node_count: usize,
    branch: usize,
    control: usize,
    r: M::T,
) {
    let row = (node_count - 1) + branch;
    let col = (node_count - 1) + control;
    *m.at(row * dim + col) -= r;
}

fn nv(state: &[f64], node: usize) -> f64 {
    if node > 0 {
        state[node - 1]
    } else {
        0.0
    }
}

fn local_values(terminals: &[usize], nodes: &[f64], what: &str) -> Result<Vec<f64>> {
    terminals
        .iter()
        .map(|&t| {
            nodes.get(t).copied().ok_or_else(|| {
                SolverError::InvalidInput(format!("{what} terminal node out of range"))
            })
        })
        .collect()
}

fn stamp_static_device<M: MnaMatrix>(
    m: &mut M,
    rhs: Option<&mut [M::T]>,
    dim: usize,
    device: &StaticDevice,
    nodes: &[f64],
) -> Result<()> {
    let t = &device.terminals;
    let local = local_values(t, nodes, "static-device")?;
    let ev = (device.evaluator)(&local)?;
    if ev.terminal_current.len() != t.len() || ev.jacobian.len() != t.len() * t.len() {
        return Err(err_runtime(
            "static-device evaluator returned invalid dimensions",
        ));
    }
    let mut rhs = rhs;
    for r in 0..t.len() {
        let row = t[r];
        if row == 0 {
            continue;
        }
        let mut ieq = ev.terminal_current[r];
        for c in 0..t.len() {
            let g = ev.jacobian[r * t.len() + c];
            if !g.is_finite() {
                return Err(err_runtime("static-device Jacobian is non-finite"));
            }
            ieq -= g * local[c];
            if t[c] > 0 {
                *m.at((row - 1) * dim + (t[c] - 1)) += M::T::from(g);
            }
        }
        if !ieq.is_finite() {
            return Err(err_runtime(
                "static-device equivalent current is non-finite",
            ));
        }
        if let Some(rhs) = rhs.as_deref_mut() {
            rhs[row - 1] -= M::T::from(ieq);
        }
    }
    Ok(())
}

fn evaluate_dynamic_device(
    device: &DynamicDevice,
    nodes: &[f64],
) -> Result<DynamicDeviceEvaluation> {
    let local = local_values(&device.terminals, nodes, "dynamic-device")?;
    let ev = (device.evaluator)(&local)?;
    let n = device.terminals.len();
    if ev.terminal_current.len() != n
        || ev.terminal_charge.len() != n
        || ev.current_jacobian.len() != n * n
        || ev.charge_jacobian.len() != n * n
    {
        return Err(err_runtime(
            "dynamic-device evaluator returned invalid dimensions",
        ));
    }
    let finite = |v: &[f64]| v.iter().all(|x| x.is_finite());
    if !finite(&ev.terminal_current)
        || !finite(&ev.terminal_charge)
        || !finite(&ev.current_jacobian)
        || !finite(&ev.charge_jacobian)
    {
        return Err(err_runtime(
            "dynamic-device evaluator returned non-finite values",
        ));
    }
    Ok(ev)
}

fn stamp_dynamic_dc(
    m: &mut RealMna,
    rhs: &mut [f64],
    dim: usize,
    device: &DynamicDevice,
    nodes: &[f64],
) -> Result<()> {
    let ev = evaluate_dynamic_device(device, nodes)?;
    let t = &device.terminals;
    for r in 0..t.len() {
        let row = t[r];
        if row == 0 {
            continue;
        }
        let mut ieq = ev.terminal_current[r];
        for c in 0..t.len() {
            let g = ev.current_jacobian[r * t.len() + c];
            ieq -= g * nodes[t[c]];
            if t[c] > 0 {
                *m.at((row - 1) * dim + (t[c] - 1)) += g;
            }
        }
        rhs[row - 1] -= ieq;
    }
    Ok(())
}

fn stamp_dynamic_ac(
    m: &mut [Complex],
    dim: usize,
    device: &DynamicDevice,
    nodes: &[f64],
    omega: f64,
) -> Result<()> {
    let ev = evaluate_dynamic_device(device, nodes)?;
    let t = &device.terminals;
    let jw = Complex::new(0.0, omega);
    for r in 0..t.len() {
        if t[r] == 0 {
            continue;
        }
        for c in 0..t.len() {
            if t[c] == 0 {
                continue;
            }
            let k = r * t.len() + c;
            m[(t[r] - 1) * dim + (t[c] - 1)] +=
                Complex::new(ev.current_jacobian[k], 0.0) + jw * ev.charge_jacobian[k];
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn stamp_dynamic_transient(
    m: &mut RealMna,
    rhs: &mut [f64],
    dim: usize,
    device: &DynamicDevice,
    current_nodes: &[f64],
    old_nodes: &[f64],
    older_nodes: Option<&[f64]>,
    a: [f64; 3],
) -> Result<()> {
    let current = evaluate_dynamic_device(device, current_nodes)?;
    let old = evaluate_dynamic_device(device, old_nodes)?;
    let older = match older_nodes {
        Some(nodes) => Some(evaluate_dynamic_device(device, nodes)?),
        None => None,
    };
    let t = &device.terminals;
    for r in 0..t.len() {
        let row = t[r];
        if row == 0 {
            continue;
        }
        let mut residual = current.terminal_current[r]
            + a[0] * current.terminal_charge[r]
            + a[1] * old.terminal_charge[r];
        if let Some(older) = &older {
            residual += a[2] * older.terminal_charge[r];
        }
        let mut ieq = residual;
        for c in 0..t.len() {
            let k = r * t.len() + c;
            let jac = current.current_jacobian[k] + a[0] * current.charge_jacobian[k];
            ieq -= jac * current_nodes[t[c]];
            if t[c] > 0 {
                *m.at((row - 1) * dim + (t[c] - 1)) += jac;
            }
        }
        rhs[row - 1] -= ieq;
    }
    Ok(())
}

pub(crate) struct BehavioralLinearization {
    value: f64,
    derivative: Vec<f64>,
}

fn behavioral_linearize(
    evaluator: &BehavioralEvaluator,
    nodes: &[f64],
    time: f64,
) -> Result<BehavioralLinearization> {
    let value = evaluator(nodes, time)?;
    if !value.is_finite() {
        return Err(err_runtime("behavioral source returned non-finite value"));
    }
    let mut derivative = vec![0.0; nodes.len()];
    let mut perturbed = nodes.to_vec();
    for node in 1..nodes.len() {
        let h = 1.0e-7 * nodes[node].abs().max(1.0);
        perturbed[node] = nodes[node] + h;
        let plus = evaluator(&perturbed, time)?;
        perturbed[node] = nodes[node] - h;
        let minus = evaluator(&perturbed, time)?;
        perturbed[node] = nodes[node];
        if !plus.is_finite() || !minus.is_finite() {
            return Err(err_runtime("behavioral-source Jacobian is non-finite"));
        }
        derivative[node] = (plus - minus) / (2.0 * h);
    }
    Ok(BehavioralLinearization { value, derivative })
}

#[allow(clippy::too_many_arguments)]
fn stamp_behavioral_current<M: MnaMatrix>(
    m: &mut M,
    rhs: Option<&mut [M::T]>,
    dim: usize,
    p: usize,
    n: usize,
    lin: &BehavioralLinearization,
    nodes: &[f64],
    scale: f64,
) {
    let mut ieq = scale * lin.value;
    for node in 1..lin.derivative.len() {
        let g = scale * lin.derivative[node];
        ieq -= g * nodes[node];
        if p > 0 {
            *m.at((p - 1) * dim + (node - 1)) += M::T::from(g);
        }
        if n > 0 {
            *m.at((n - 1) * dim + (node - 1)) -= M::T::from(g);
        }
    }
    if let Some(rhs) = rhs {
        stamp_current(rhs, p, n, M::T::from(ieq));
    }
}

#[allow(clippy::too_many_arguments)]
fn stamp_behavioral_voltage<M: MnaMatrix>(
    m: &mut M,
    rhs: Option<&mut [M::T]>,
    dim: usize,
    node_count: usize,
    branch: usize,
    p: usize,
    n: usize,
    lin: &BehavioralLinearization,
    nodes: &[f64],
    scale: f64,
) {
    let mut dummy = vec![M::T::default(); dim];
    let has_rhs = rhs.is_some();
    let target: &mut [M::T] = match rhs {
        Some(r) => r,
        None => &mut dummy,
    };
    stamp_voltage_branch(
        m,
        target,
        dim,
        node_count,
        branch,
        p,
        n,
        M::T::default(),
        M::T::default(),
    );
    let row = (node_count - 1) + branch;
    let mut ieq = scale * lin.value;
    for node in 1..lin.derivative.len() {
        let g = scale * lin.derivative[node];
        ieq -= g * nodes[node];
        *m.at(row * dim + (node - 1)) -= M::T::from(g);
    }
    if has_rhs {
        target[row] += M::T::from(ieq);
    }
}

fn node_values_from_state(state: &[f64], node_count: usize) -> Vec<f64> {
    let mut values = vec![0.0; node_count];
    values[1..node_count].copy_from_slice(&state[..node_count - 1]);
    values
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct MosLinearization {
    pub gm: f64,
    pub gds: f64,
    pub ieq: f64,
}

fn square_law(kp: f64, lambda: f64, overdrive: f64, vds: f64) -> (f64, f64, f64) {
    if vds < overdrive {
        let base = kp * (overdrive * vds - 0.5 * vds * vds);
        (
            base * (1.0 + lambda * vds),
            kp * vds * (1.0 + lambda * vds),
            kp * ((overdrive - vds) * (1.0 + lambda * vds)
                + (overdrive * vds - 0.5 * vds * vds) * lambda),
        )
    } else {
        (
            0.5 * kp * overdrive * overdrive * (1.0 + lambda * vds),
            kp * overdrive * (1.0 + lambda * vds),
            0.5 * kp * overdrive * overdrive * lambda,
        )
    }
}

pub(crate) fn mos_linearize(model: &MosLevel1Model, vd: f64, vg: f64, vs: f64) -> MosLinearization {
    let polarity = if model.pmos { -1.0 } else { 1.0 };
    let vgs = polarity * (vg - vs);
    let vds = polarity * (vd - vs);
    let threshold = model.threshold_voltage.abs();
    if vgs <= threshold || vds <= 0.0 {
        return MosLinearization::default();
    }
    let (eff, gm, gds) = square_law(
        model.transconductance_parameter,
        model.channel_length_modulation,
        vgs - threshold,
        vds,
    );
    let id = polarity * eff;
    MosLinearization {
        gm,
        gds,
        ieq: id - gm * (vg - vs) - gds * (vd - vs),
    }
}

pub(crate) fn jfet_linearize(model: &JfetModel, vd: f64, vg: f64, vs: f64) -> MosLinearization {
    let polarity = if model.pchannel { -1.0 } else { 1.0 };
    let vgs = polarity * (vg - vs);
    let vds = polarity * (vd - vs);
    let pinch = if model.pchannel {
        -model.pinch_off_voltage.abs()
    } else {
        model.pinch_off_voltage
    };
    if vgs <= pinch || vds <= 0.0 {
        return MosLinearization::default();
    }
    let (eff, gm, gds) = square_law(
        model.beta,
        model.channel_length_modulation,
        vgs - pinch,
        vds,
    );
    let id = polarity * eff;
    MosLinearization {
        gm,
        gds,
        ieq: id - gm * (vg - vs) - gds * (vd - vs),
    }
}

fn stamp_transistor<M: MnaMatrix>(
    m: &mut M,
    dim: usize,
    lin: &MosLinearization,
    d: usize,
    g: usize,
    s: usize,
) {
    let mut add = |row: usize, col: usize, v: f64| {
        if row > 0 && col > 0 {
            *m.at((row - 1) * dim + (col - 1)) += M::T::from(v);
        }
    };
    add(d, d, lin.gds);
    add(d, g, lin.gm);
    add(d, s, -lin.gds - lin.gm);
    add(s, d, -lin.gds);
    add(s, g, -lin.gm);
    add(s, s, lin.gds + lin.gm);
}

pub(crate) fn bjt_currents(model: &BjtModel, vc: f64, vb: f64, ve: f64) -> [f64; 3] {
    let polarity = if model.pnp { -1.0 } else { 1.0 };
    let vt = BOLTZMANN * model.temperature_k / ELECTRON_CHARGE;
    let vbe = polarity * (vb - ve);
    let vbc = polarity * (vb - vc);
    let forward = model.saturation_current * ((vbe / vt).clamp(-40.0, 40.0).exp() - 1.0);
    let reverse = model.saturation_current * ((vbc / vt).clamp(-40.0, 40.0).exp() - 1.0);
    let ic = forward - reverse * (1.0 + 1.0 / model.reverse_beta);
    let ib = forward / model.forward_beta + reverse / model.reverse_beta;
    let ie = -ic - ib;
    [polarity * ic, polarity * ib, polarity * ie]
}

struct BjtLinearization {
    ieq: [f64; 3],
    jacobian: [[f64; 3]; 3],
}

fn bjt_linearize(model: &BjtModel, vc: f64, vb: f64, ve: f64) -> BjtLinearization {
    let voltage = [vc, vb, ve];
    let current = bjt_currents(model, vc, vb, ve);
    let mut jacobian = [[0.0; 3]; 3];
    let h = 1.0e-6;
    for col in 0..3 {
        let (mut plus, mut minus) = (voltage, voltage);
        plus[col] += h;
        minus[col] -= h;
        let ip = bjt_currents(model, plus[0], plus[1], plus[2]);
        let im = bjt_currents(model, minus[0], minus[1], minus[2]);
        for row in 0..3 {
            jacobian[row][col] = (ip[row] - im[row]) / (2.0 * h);
        }
    }
    let mut ieq = [0.0; 3];
    for row in 0..3 {
        ieq[row] = current[row];
        for col in 0..3 {
            ieq[row] -= jacobian[row][col] * voltage[col];
        }
    }
    BjtLinearization { ieq, jacobian }
}

fn stamp_bjt_jacobian<M: MnaMatrix>(
    m: &mut M,
    dim: usize,
    lin: &BjtLinearization,
    nodes: [usize; 3],
) {
    for row in 0..3 {
        if nodes[row] == 0 {
            continue;
        }
        for col in 0..3 {
            if nodes[col] == 0 {
                continue;
            }
            *m.at((nodes[row] - 1) * dim + (nodes[col] - 1)) += M::T::from(lin.jacobian[row][col]);
        }
    }
}

fn stamp_bjt_dc(
    m: &mut RealMna,
    rhs: &mut [f64],
    dim: usize,
    lin: &BjtLinearization,
    nodes: [usize; 3],
) {
    stamp_bjt_jacobian(m, dim, lin, nodes);
    for t in 0..3 {
        if nodes[t] > 0 {
            rhs[nodes[t] - 1] -= lin.ieq[t];
        }
    }
}

fn switch_conductance(model: &SwitchModel, control: f64) -> f64 {
    let up = model.threshold_voltage + model.hysteresis_voltage.abs();
    let low = model.threshold_voltage - model.hysteresis_voltage.abs();
    let on = control >= up || (control > low && control >= model.threshold_voltage);
    1.0 / if on {
        model.on_resistance
    } else {
        model.off_resistance
    }
}

fn interpolate_nport(samples: &[NPortPoint], frequency: f64) -> Result<NPortPoint> {
    let (Some(first), Some(last)) = (samples.first(), samples.last()) else {
        return Err(SolverError::InvalidInput(
            "N-port requires at least one frequency sample".into(),
        ));
    };
    if samples.len() == 1 || frequency <= first.frequency_hz {
        return Ok(first.clone());
    }
    if frequency >= last.frequency_hz {
        return Ok(last.clone());
    }
    let upper_index = samples.partition_point(|p| p.frequency_hz <= frequency);
    let (lo, hi) = (&samples[upper_index - 1], &samples[upper_index]);
    if lo.s.size() != hi.s.size() || lo.reference_impedance != hi.reference_impedance {
        return Err(err_runtime(
            "N-port interpolation topology/reference mismatch",
        ));
    }
    let fraction = (frequency - lo.frequency_hz) / (hi.frequency_hz - lo.frequency_hz);
    let n = lo.s.size();
    let mut s = ComplexMatrix::new(n);
    for r in 0..n {
        for c in 0..n {
            s[(r, c)] = lo.s[(r, c)] * (1.0 - fraction) + hi.s[(r, c)] * fraction;
        }
    }
    Ok(NPortPoint {
        frequency_hz: frequency,
        s,
        reference_impedance: lo.reference_impedance.clone(),
    })
}

fn stamp_nport(m: &mut [Complex], dim: usize, device: &NPortDevice, frequency: f64) -> Result<()> {
    let point = interpolate_nport(&device.samples, frequency)?;
    if point.s.size() != device.ports.len() {
        return Err(err_runtime("N-port device port count mismatch"));
    }
    let z = rf::nport_s_to_z(&point.s, &point.reference_impedance)?;
    let y = rf::inverse(&z)?;
    let mut add = |row: usize, col: usize, v: Complex| {
        if row > 0 && col > 0 {
            m[(row - 1) * dim + (col - 1)] += v;
        }
    };
    for (i, &(pi, ni)) in device.ports.iter().enumerate() {
        for (j, &(pj, nj)) in device.ports.iter().enumerate() {
            let v = y[(i, j)];
            add(pi, pj, v);
            add(pi, nj, -v);
            add(ni, pj, -v);
            add(ni, nj, v);
        }
    }
    Ok(())
}

fn diode_vt(model: &DiodeModel) -> f64 {
    model.emission_coefficient * BOLTZMANN * model.temperature_k / ELECTRON_CHARGE
}

fn model_of<'a, T>(map: &'a HashMap<String, T>, name: &str, kind: &str) -> Result<&'a T> {
    map.get(name)
        .ok_or_else(|| err_runtime(format!("unknown {kind} model: {name}")))
}

/// Integration companion-model selection for one Newton assembly.
#[derive(Clone, Copy)]
enum Companion<'a> {
    Dc {
        scale: f64,
    },
    Fixed {
        method: TransientMethod,
        dt: f64,
        time: f64,
        capacitor_current: &'a [f64],
        inductor_voltage: &'a [f64],
    },
    Adaptive {
        time: f64,
        a: [f64; 3],
    },
}

// ----------------------------------------------------------------------------------------------
// Circuit construction
// ----------------------------------------------------------------------------------------------

impl Circuit {
    pub fn new() -> Self {
        let mut node_lookup = HashMap::new();
        node_lookup.insert("0".to_string(), 0);
        node_lookup.insert("GND".to_string(), 0);
        Self {
            node_names: vec!["0".into()],
            node_lookup,
            resistors: Vec::new(),
            capacitors: Vec::new(),
            inductors: Vec::new(),
            mutual_inductances: Vec::new(),
            voltage_sources: Vec::new(),
            current_sources: Vec::new(),
            vccs: Vec::new(),
            vcvs: Vec::new(),
            cccs: Vec::new(),
            ccvs: Vec::new(),
            diodes: Vec::new(),
            mosfets: Vec::new(),
            bjts: Vec::new(),
            jfets: Vec::new(),
            switches: Vec::new(),
            nports: Vec::new(),
            static_devices: Vec::new(),
            dynamic_devices: Vec::new(),
            behavioral_current_sources: Vec::new(),
            behavioral_voltage_sources: Vec::new(),
            diode_models: HashMap::new(),
            mos_models: HashMap::new(),
            bjt_models: HashMap::new(),
            jfet_models: HashMap::new(),
            switch_models: HashMap::new(),
            initial_voltage: HashMap::new(),
            nodeset_voltage: HashMap::new(),
        }
    }

    /// Node id for `name`, creating it if needed.
    pub fn node(&mut self, name: &str) -> usize {
        let mut key = upper(name);
        if key == "GND" {
            key = "0".into();
        }
        if let Some(&id) = self.node_lookup.get(&key) {
            return id;
        }
        let id = self.node_names.len();
        self.node_names.push(name.to_string());
        self.node_lookup.insert(key, id);
        id
    }
    pub fn find_node(&self, name: &str) -> Result<usize> {
        let mut key = upper(name);
        if key == "GND" {
            key = "0".into();
        }
        self.node_lookup
            .get(&key)
            .copied()
            .ok_or_else(|| SolverError::InvalidInput("unknown circuit node".into()))
    }
    pub fn node_count(&self) -> usize {
        self.node_names.len()
    }
    pub fn node_names(&self) -> &[String] {
        &self.node_names
    }
    pub fn initial_voltages(&self) -> &HashMap<usize, f64> {
        &self.initial_voltage
    }
    pub fn nodeset_voltages(&self) -> &HashMap<usize, f64> {
        &self.nodeset_voltage
    }

    fn check_nodes(&self, nodes: &[usize]) -> Result<()> {
        require(
            nodes.iter().all(|&n| n < self.node_count()),
            "circuit element node out of range",
        )
    }

    pub fn add_resistor(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        value: f64,
    ) -> Result<()> {
        require(value > 0.0 && value.is_finite(), "invalid resistance")?;
        self.check_nodes(&[p, n])?;
        self.resistors.push(TwoTerminal {
            name: name.into(),
            p,
            n,
            value,
        });
        Ok(())
    }
    pub fn add_capacitor(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        value: f64,
    ) -> Result<()> {
        require(value >= 0.0 && value.is_finite(), "invalid capacitance")?;
        self.check_nodes(&[p, n])?;
        self.capacitors.push(TwoTerminal {
            name: name.into(),
            p,
            n,
            value,
        });
        Ok(())
    }
    pub fn add_inductor(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        value: f64,
    ) -> Result<()> {
        require(value > 0.0 && value.is_finite(), "invalid inductance")?;
        self.check_nodes(&[p, n])?;
        self.inductors.push(TwoTerminal {
            name: name.into(),
            p,
            n,
            value,
        });
        Ok(())
    }
    pub fn add_mutual_inductance(
        &mut self,
        _name: impl Into<String>,
        first: &str,
        second: &str,
        coupling: f64,
    ) -> Result<()> {
        require(
            coupling.is_finite() && coupling.abs() <= 1.0,
            "mutual-inductor coupling must be in [-1,1]",
        )?;
        self.mutual_inductances.push(MutualInductance {
            first: upper(first),
            second: upper(second),
            coupling,
        });
        Ok(())
    }
    pub fn add_voltage_source(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        dc: f64,
        ac: Complex,
        transient: Option<SourceWaveform>,
    ) -> Result<()> {
        require(
            dc.is_finite() && ac.re.is_finite() && ac.im.is_finite(),
            "invalid voltage source",
        )?;
        self.check_nodes(&[p, n])?;
        self.voltage_sources.push(IndependentSource {
            name: name.into(),
            p,
            n,
            dc,
            ac,
            transient,
        });
        Ok(())
    }
    pub fn add_current_source(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        dc: f64,
        ac: Complex,
        transient: Option<SourceWaveform>,
    ) -> Result<()> {
        require(
            dc.is_finite() && ac.re.is_finite() && ac.im.is_finite(),
            "invalid current source",
        )?;
        self.check_nodes(&[p, n])?;
        self.current_sources.push(IndependentSource {
            name: name.into(),
            p,
            n,
            dc,
            ac,
            transient,
        });
        Ok(())
    }
    pub fn add_vccs(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        cp: usize,
        cn: usize,
        gain: f64,
    ) -> Result<()> {
        require(gain.is_finite(), "invalid VCCS gain")?;
        self.check_nodes(&[p, n, cp, cn])?;
        self.vccs.push(VoltageControlled {
            name: name.into(),
            p,
            n,
            cp,
            cn,
            gain,
        });
        Ok(())
    }
    pub fn add_vcvs(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        cp: usize,
        cn: usize,
        gain: f64,
    ) -> Result<()> {
        require(gain.is_finite(), "invalid VCVS gain")?;
        self.check_nodes(&[p, n, cp, cn])?;
        self.vcvs.push(VoltageControlled {
            name: name.into(),
            p,
            n,
            cp,
            cn,
            gain,
        });
        Ok(())
    }
    pub fn add_cccs(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        control: &str,
        gain: f64,
    ) -> Result<()> {
        require(gain.is_finite(), "invalid CCCS gain")?;
        self.check_nodes(&[p, n])?;
        self.cccs.push(CurrentControlled {
            name: name.into(),
            p,
            n,
            control: upper(control),
            gain,
        });
        Ok(())
    }
    pub fn add_ccvs(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        control: &str,
        gain: f64,
    ) -> Result<()> {
        require(gain.is_finite(), "invalid CCVS transresistance")?;
        self.check_nodes(&[p, n])?;
        self.ccvs.push(CurrentControlled {
            name: name.into(),
            p,
            n,
            control: upper(control),
            gain,
        });
        Ok(())
    }
    /// Ideal lossless transformer (`turns_ratio = Vprimary/Vsecondary`) built from a VCVS and a
    /// power-conserving CCCS named `<name>$V` and `<name>$I`.
    pub fn add_ideal_transformer(
        &mut self,
        name: &str,
        pp: usize,
        pn: usize,
        sp: usize,
        sn: usize,
        turns_ratio: f64,
    ) -> Result<()> {
        require(
            turns_ratio > 0.0 && turns_ratio.is_finite(),
            "invalid ideal-transformer turns ratio",
        )?;
        let v = format!("{name}$V");
        self.add_vcvs(v.clone(), pp, pn, sp, sn, turns_ratio)?;
        self.add_cccs(format!("{name}$I"), sp, sn, &v, -turns_ratio)
    }
    pub fn add_diode_model(&mut self, name: &str, model: DiodeModel) -> Result<()> {
        require(
            model.saturation_current > 0.0
                && model.emission_coefficient > 0.0
                && model.temperature_k > 0.0,
            "invalid diode model",
        )?;
        self.diode_models.insert(upper(name), model);
        Ok(())
    }
    pub fn add_diode(
        &mut self,
        name: impl Into<String>,
        anode: usize,
        cathode: usize,
        model: &str,
    ) -> Result<()> {
        self.check_nodes(&[anode, cathode])?;
        self.diodes.push(ModelDevice {
            name: name.into(),
            nodes: [anode, cathode, 0, 0],
            model: upper(model),
        });
        Ok(())
    }
    pub fn add_mos_model(&mut self, name: &str, model: MosLevel1Model) -> Result<()> {
        require(
            model.transconductance_parameter > 0.0
                && model.channel_length_modulation >= 0.0
                && model.threshold_voltage.is_finite(),
            "invalid MOS level-1 model",
        )?;
        self.mos_models.insert(upper(name), model);
        Ok(())
    }
    pub fn add_mosfet(
        &mut self,
        name: impl Into<String>,
        d: usize,
        g: usize,
        s: usize,
        b: usize,
        model: &str,
    ) -> Result<()> {
        self.check_nodes(&[d, g, s, b])?;
        self.mosfets.push(ModelDevice {
            name: name.into(),
            nodes: [d, g, s, b],
            model: upper(model),
        });
        Ok(())
    }
    pub fn add_bjt_model(&mut self, name: &str, model: BjtModel) -> Result<()> {
        require(
            model.saturation_current > 0.0
                && model.forward_beta > 0.0
                && model.reverse_beta > 0.0
                && model.temperature_k > 0.0,
            "invalid BJT model",
        )?;
        self.bjt_models.insert(upper(name), model);
        Ok(())
    }
    pub fn add_bjt(
        &mut self,
        name: impl Into<String>,
        c: usize,
        b: usize,
        e: usize,
        model: &str,
    ) -> Result<()> {
        self.check_nodes(&[c, b, e])?;
        self.bjts.push(ModelDevice {
            name: name.into(),
            nodes: [c, b, e, 0],
            model: upper(model),
        });
        Ok(())
    }
    pub fn add_jfet_model(&mut self, name: &str, model: JfetModel) -> Result<()> {
        require(
            model.beta > 0.0
                && model.channel_length_modulation >= 0.0
                && model.pinch_off_voltage.is_finite(),
            "invalid JFET model",
        )?;
        self.jfet_models.insert(upper(name), model);
        Ok(())
    }
    pub fn add_jfet(
        &mut self,
        name: impl Into<String>,
        d: usize,
        g: usize,
        s: usize,
        model: &str,
    ) -> Result<()> {
        self.check_nodes(&[d, g, s])?;
        self.jfets.push(ModelDevice {
            name: name.into(),
            nodes: [d, g, s, 0],
            model: upper(model),
        });
        Ok(())
    }
    pub fn add_switch_model(&mut self, name: &str, model: SwitchModel) -> Result<()> {
        require(
            model.on_resistance > 0.0
                && model.off_resistance > model.on_resistance
                && model.threshold_voltage.is_finite()
                && model.hysteresis_voltage.is_finite(),
            "invalid voltage-controlled switch model",
        )?;
        self.switch_models.insert(upper(name), model);
        Ok(())
    }
    pub fn add_switch(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        cp: usize,
        cn: usize,
        model: &str,
    ) -> Result<()> {
        self.check_nodes(&[p, n, cp, cn])?;
        self.switches.push(ModelDevice {
            name: name.into(),
            nodes: [p, n, cp, cn],
            model: upper(model),
        });
        Ok(())
    }
    /// Frequency-domain N-port between each (positive, negative) node pair; sampled S-parameters
    /// are linearly interpolated in frequency (clamped outside the sampled band).
    pub fn add_nport(
        &mut self,
        _name: impl Into<String>,
        ports: Vec<(usize, usize)>,
        mut samples: Vec<NPortPoint>,
    ) -> Result<()> {
        require(
            !ports.is_empty() && !samples.is_empty(),
            "N-port device requires ports and samples",
        )?;
        for &(p, n) in &ports {
            self.check_nodes(&[p, n])?;
        }
        samples.sort_by(|a, b| a.frequency_hz.total_cmp(&b.frequency_hz));
        for s in &samples {
            require(
                s.frequency_hz > 0.0
                    && s.s.size() == ports.len()
                    && s.reference_impedance.len() == ports.len(),
                "invalid N-port sample",
            )?;
        }
        self.nports.push(NPortDevice { ports, samples });
        Ok(())
    }
    pub fn add_tem_transmission_line(
        &mut self,
        name: impl Into<String>,
        first_port: (usize, usize),
        second_port: (usize, usize),
        line: &TemTransmissionLine,
        frequencies: &[f64],
        reference_impedance: f64,
    ) -> Result<()> {
        require(
            !frequencies.is_empty() && reference_impedance > 0.0 && reference_impedance.is_finite(),
            "transmission-line samples/reference impedance are invalid",
        )?;
        let mut samples = Vec::with_capacity(frequencies.len());
        for &f in frequencies {
            let two = line.s_parameters(f, reference_impedance)?;
            samples.push(two_port_point(f, reference_impedance, two));
        }
        self.add_nport(name, vec![first_port, second_port], samples)
    }
    #[allow(clippy::too_many_arguments)]
    pub fn add_microstrip_line(
        &mut self,
        name: impl Into<String>,
        first_port: (usize, usize),
        second_port: (usize, usize),
        width: f64,
        height: f64,
        er: f64,
        length: f64,
        frequencies: &[f64],
        reference_impedance: f64,
    ) -> Result<()> {
        let f0 = first_sample(frequencies, length, "microstrip")?;
        let q = rf::microstrip_quasi_static(width, height, er, f0)?;
        let line = TemTransmissionLine {
            characteristic_impedance: q.characteristic_impedance,
            phase_velocity: q.phase_velocity,
            attenuation_nepers_per_m: 0.0,
            length,
        };
        self.add_tem_transmission_line(
            name,
            first_port,
            second_port,
            &line,
            frequencies,
            reference_impedance,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn add_stripline_line(
        &mut self,
        name: impl Into<String>,
        first_port: (usize, usize),
        second_port: (usize, usize),
        width: f64,
        spacing: f64,
        er: f64,
        length: f64,
        frequencies: &[f64],
        reference_impedance: f64,
    ) -> Result<()> {
        let f0 = first_sample(frequencies, length, "stripline")?;
        let q = rf::stripline_quasi_static(width, spacing, er, f0)?;
        let line = TemTransmissionLine {
            characteristic_impedance: q.characteristic_impedance,
            phase_velocity: q.phase_velocity,
            attenuation_nepers_per_m: 0.0,
            length,
        };
        self.add_tem_transmission_line(
            name,
            first_port,
            second_port,
            &line,
            frequencies,
            reference_impedance,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn add_coplanar_waveguide_line(
        &mut self,
        name: impl Into<String>,
        first_port: (usize, usize),
        second_port: (usize, usize),
        width: f64,
        gap: f64,
        er: f64,
        length: f64,
        frequencies: &[f64],
        reference_impedance: f64,
    ) -> Result<()> {
        let f0 = first_sample(frequencies, length, "coplanar-waveguide")?;
        let q = rf::coplanar_waveguide_quasi_static(width, gap, er, f0)?;
        let line = TemTransmissionLine {
            characteristic_impedance: q.characteristic_impedance,
            phase_velocity: q.phase_velocity,
            attenuation_nepers_per_m: 0.0,
            length,
        };
        self.add_tem_transmission_line(
            name,
            first_port,
            second_port,
            &line,
            frequencies,
            reference_impedance,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn add_rectangular_waveguide_te10_line(
        &mut self,
        name: impl Into<String>,
        first_port: (usize, usize),
        second_port: (usize, usize),
        broad: f64,
        narrow: f64,
        er: f64,
        mur: f64,
        length: f64,
        frequencies: &[f64],
        reference_impedance: f64,
    ) -> Result<()> {
        require(
            !frequencies.is_empty()
                && length >= 0.0
                && length.is_finite()
                && reference_impedance > 0.0,
            "waveguide circuit element requires samples and non-negative length",
        )?;
        let mut samples = Vec::with_capacity(frequencies.len());
        for &f in frequencies {
            let mode = rf::rectangular_waveguide_te10(broad, narrow, er, mur, f)?;
            let line = TransmissionLine {
                characteristic_impedance: mode.wave_impedance,
                propagation_constant: Complex::new(0.0, mode.propagation_constant_rad_per_m),
                length,
            };
            samples.push(two_port_point(
                f,
                reference_impedance,
                line.s_parameters(reference_impedance)?,
            ));
        }
        self.add_nport(name, vec![first_port, second_port], samples)
    }
    /// Memoryless nonlinear multi-terminal compact device (currents and Jacobian).
    pub fn add_static_device(
        &mut self,
        _name: impl Into<String>,
        terminals: Vec<usize>,
        evaluator: StaticDeviceEvaluator,
    ) -> Result<()> {
        require(
            !terminals.is_empty(),
            "static device requires terminals and evaluator",
        )?;
        self.check_nodes(&terminals)?;
        self.static_devices.push(StaticDevice {
            terminals,
            evaluator,
        });
        Ok(())
    }
    /// Charge-based DAE device `I(V) + dQ(V)/dt = 0`. Trapezoidal transient is rejected for
    /// these devices (as upstream).
    pub fn add_dynamic_device(
        &mut self,
        _name: impl Into<String>,
        terminals: Vec<usize>,
        evaluator: DynamicDeviceEvaluator,
    ) -> Result<()> {
        require(
            !terminals.is_empty(),
            "dynamic device requires terminals and evaluator",
        )?;
        self.check_nodes(&terminals)?;
        self.dynamic_devices.push(DynamicDevice {
            terminals,
            evaluator,
        });
        Ok(())
    }
    /// Electrothermal wrapper: the thermal node voltage is the temperature rise above ambient;
    /// Rth/Cth are stamped as thermal conductance/charge and dissipated power drives the node.
    pub fn add_electrothermal_device(
        &mut self,
        name: &str,
        electrical_terminals: Vec<usize>,
        thermal_node: usize,
        config: ElectroThermalDeviceConfig,
        evaluator: ElectroThermalDeviceEvaluator,
    ) -> Result<()> {
        require(
            !name.is_empty()
                && !electrical_terminals.is_empty()
                && thermal_node != 0
                && thermal_node < self.node_count()
                && config.ambient_temperature_k > 0.0
                && config.ambient_temperature_k.is_finite()
                && config.thermal_resistance_k_per_w > 0.0
                && config.thermal_resistance_k_per_w.is_finite()
                && config.thermal_capacitance_j_per_k >= 0.0
                && config.thermal_capacitance_j_per_k.is_finite(),
            "invalid electrothermal device configuration",
        )?;
        self.check_nodes(&electrical_terminals)?;
        let ne = electrical_terminals.len();
        let mut terminals = electrical_terminals;
        terminals.push(thermal_node);
        let dae: DynamicDeviceEvaluator = Arc::new(move |local: &[f64]| {
            electrothermal_evaluation(local, ne, &config, &evaluator)
        });
        self.add_dynamic_device(name, terminals, dae)
    }
    pub fn add_behavioral_current_source(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        evaluator: BehavioralEvaluator,
    ) -> Result<()> {
        let name = name.into();
        require(
            !name.is_empty() && p < self.node_count() && n < self.node_count(),
            "invalid behavioral current source",
        )?;
        self.behavioral_current_sources.push(BehavioralSource {
            name,
            p,
            n,
            evaluator,
        });
        Ok(())
    }
    pub fn add_behavioral_voltage_source(
        &mut self,
        name: impl Into<String>,
        p: usize,
        n: usize,
        evaluator: BehavioralEvaluator,
    ) -> Result<()> {
        let name = name.into();
        require(
            !name.is_empty() && p < self.node_count() && n < self.node_count(),
            "invalid behavioral voltage source",
        )?;
        self.behavioral_voltage_sources.push(BehavioralSource {
            name,
            p,
            n,
            evaluator,
        });
        Ok(())
    }

    pub fn set_voltage_source_dc(&mut self, name: &str, voltage: f64) -> Result<()> {
        require(voltage.is_finite(), "invalid voltage-source DC value")?;
        let target = upper(name);
        let source = self
            .voltage_sources
            .iter_mut()
            .find(|s| upper(&s.name) == target)
            .ok_or_else(|| SolverError::InvalidInput("unknown voltage source".into()))?;
        source.dc = voltage;
        Ok(())
    }
    pub fn voltage_source_dc(&self, name: &str) -> Result<f64> {
        let target = upper(name);
        self.voltage_sources
            .iter()
            .find(|s| upper(&s.name) == target)
            .map(|s| s.dc)
            .ok_or_else(|| SolverError::InvalidInput("unknown voltage source".into()))
    }
    pub fn set_resistance(&mut self, name: &str, value: f64) -> Result<()> {
        require(value > 0.0 && value.is_finite(), "invalid resistance")?;
        let target = upper(name);
        let r = self
            .resistors
            .iter_mut()
            .find(|r| upper(&r.name) == target)
            .ok_or_else(|| SolverError::InvalidInput("unknown resistor".into()))?;
        r.value = value;
        Ok(())
    }
    pub fn resistance(&self, name: &str) -> Result<f64> {
        let target = upper(name);
        self.resistors
            .iter()
            .find(|r| upper(&r.name) == target)
            .map(|r| r.value)
            .ok_or_else(|| SolverError::InvalidInput("unknown resistor".into()))
    }
    /// Update temperature-dependent compact models (diode and BJT).
    pub fn set_device_temperature(&mut self, temperature: f64) -> Result<()> {
        require(
            temperature > 0.0 && temperature.is_finite(),
            "device temperature must be finite and positive",
        )?;
        for m in self.diode_models.values_mut() {
            m.temperature_k = temperature;
        }
        for m in self.bjt_models.values_mut() {
            m.temperature_k = temperature;
        }
        Ok(())
    }
    pub fn set_initial_voltage(&mut self, node: usize, voltage: f64) -> Result<()> {
        require(
            node < self.node_count() && voltage.is_finite(),
            "invalid circuit initial condition",
        )?;
        if node == 0 {
            return require(voltage == 0.0, "ground initial condition must be zero");
        }
        self.initial_voltage.insert(node, voltage);
        Ok(())
    }
    pub fn set_nodeset_voltage(&mut self, node: usize, voltage: f64) -> Result<()> {
        require(
            node < self.node_count() && voltage.is_finite(),
            "invalid circuit nodeset",
        )?;
        if node == 0 {
            return require(voltage == 0.0, "ground nodeset must be zero");
        }
        self.nodeset_voltage.insert(node, voltage);
        Ok(())
    }

    pub(crate) fn branch_count(&self) -> usize {
        self.voltage_sources.len()
            + self.behavioral_voltage_sources.len()
            + self.vcvs.len()
            + self.ccvs.len()
            + self.inductors.len()
    }
    fn inductor_branch(&self, index: usize) -> usize {
        self.voltage_sources.len()
            + self.behavioral_voltage_sources.len()
            + self.vcvs.len()
            + self.ccvs.len()
            + index
    }
    pub(crate) fn branch_index(&self, name: &str) -> Result<usize> {
        let target = upper(name);
        let names = self
            .voltage_sources
            .iter()
            .map(|s| &s.name)
            .chain(self.behavioral_voltage_sources.iter().map(|s| &s.name))
            .chain(self.vcvs.iter().map(|s| &s.name))
            .chain(self.ccvs.iter().map(|s| &s.name))
            .chain(self.inductors.iter().map(|s| &s.name));
        for (i, n) in names.enumerate() {
            if upper(n) == target {
                return Ok(i);
            }
        }
        Err(SolverError::InvalidInput(
            "unknown branch-current control element".into(),
        ))
    }
    fn nonlinear(&self) -> bool {
        !self.diodes.is_empty()
            || !self.mosfets.is_empty()
            || !self.jfets.is_empty()
            || !self.bjts.is_empty()
            || !self.switches.is_empty()
            || !self.static_devices.is_empty()
            || !self.dynamic_devices.is_empty()
            || !self.behavioral_current_sources.is_empty()
            || !self.behavioral_voltage_sources.is_empty()
    }
    fn mutual_pair(&self, m: &MutualInductance) -> Result<(usize, usize, f64)> {
        let (mut first, mut second) = (None, None);
        for (i, l) in self.inductors.iter().enumerate() {
            let name = upper(&l.name);
            if name == m.first {
                first = Some(i);
            }
            if name == m.second {
                second = Some(i);
            }
        }
        let (Some(a), Some(b)) = (first, second) else {
            return Err(err_runtime("mutual inductance references unknown inductor"));
        };
        Ok((
            a,
            b,
            m.coupling * (self.inductors[a].value * self.inductors[b].value).sqrt(),
        ))
    }
    fn state_to_point(&self, state: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let nu = self.node_count() - 1;
        (
            node_values_from_state(state, self.node_count()),
            state[nu..].to_vec(),
        )
    }

    // ------------------------------------------------------------------------------------------
    // Real (DC / transient) assembly
    // ------------------------------------------------------------------------------------------

    /// One Newton assembly. `previous` is the Newton iterate, `old`/`older` the accepted history.
    fn assemble_real(
        &self,
        cfg: &NewtonConfig,
        mode: Companion<'_>,
        previous: &[f64],
        old: &[f64],
        older: Option<&[f64]>,
    ) -> Result<(RealMna, Vec<f64>)> {
        let nc = self.node_count();
        let nu = nc - 1;
        let dim = nu + self.branch_count();
        let sparse = cfg.sparse_mna_threshold > 0 && dim >= cfg.sparse_mna_threshold;
        let mut m = RealMna::new(dim, sparse)?;
        let mut rhs = vec![0.0; dim];
        let (time, scale) = match mode {
            Companion::Dc { scale } => (None, scale),
            Companion::Fixed { time, .. } | Companion::Adaptive { time, .. } => (Some(time), 1.0),
        };
        let source_value = |dc: f64, w: &Option<SourceWaveform>| -> Result<f64> {
            match (time, w) {
                (None, _) => Ok(scale * dc),
                (Some(t), Some(w)) => w(t),
                (Some(_), None) => Ok(dc),
            }
        };
        for r in &self.resistors {
            stamp_admittance(&mut m, dim, r.p, r.n, 1.0 / r.value);
        }
        for node in 1..nc {
            *m.at((node - 1) * dim + (node - 1)) += cfg.gmin;
        }
        for s in &self.current_sources {
            stamp_current(&mut rhs, s.p, s.n, source_value(s.dc, &s.transient)?);
        }
        for s in &self.vccs {
            stamp_vccs(&mut m, dim, s, s.gain);
        }
        let bnodes = node_values_from_state(previous, nc);
        let btime = time.unwrap_or(0.0);
        for s in &self.behavioral_current_sources {
            let lin = behavioral_linearize(&s.evaluator, &bnodes, btime)?;
            stamp_behavioral_current(&mut m, Some(&mut rhs), dim, s.p, s.n, &lin, &bnodes, scale);
        }
        if !matches!(mode, Companion::Dc { .. }) {
            for (k, c) in self.capacitors.iter().enumerate() {
                let v1 = nv(old, c.p) - nv(old, c.n);
                let v2 = older.map_or(0.0, |o| nv(o, c.p) - nv(o, c.n));
                let (g, history) = match mode {
                    Companion::Fixed {
                        method: TransientMethod::BackwardEuler,
                        dt,
                        ..
                    } => {
                        let g = c.value / dt;
                        (g, -g * v1)
                    }
                    Companion::Fixed {
                        method: TransientMethod::Trapezoidal,
                        dt,
                        capacitor_current,
                        ..
                    } => {
                        let g = 2.0 * c.value / dt;
                        (g, -g * v1 - capacitor_current[k])
                    }
                    Companion::Fixed {
                        method: TransientMethod::Bdf2,
                        dt,
                        ..
                    } => (1.5 * c.value / dt, c.value / (2.0 * dt) * (-4.0 * v1 + v2)),
                    Companion::Adaptive { a, .. } => {
                        (c.value * a[0], c.value * (a[1] * v1 + a[2] * v2))
                    }
                    Companion::Dc { .. } => (0.0, 0.0),
                };
                stamp_admittance(&mut m, dim, c.p, c.n, g);
                stamp_current(&mut rhs, c.p, c.n, history);
            }
        }
        let mut branch = 0;
        for s in &self.voltage_sources {
            let v = source_value(s.dc, &s.transient)?;
            stamp_voltage_branch(&mut m, &mut rhs, dim, nc, branch, s.p, s.n, v, 0.0);
            branch += 1;
        }
        for s in &self.behavioral_voltage_sources {
            let lin = behavioral_linearize(&s.evaluator, &bnodes, btime)?;
            stamp_behavioral_voltage(
                &mut m,
                Some(&mut rhs),
                dim,
                nc,
                branch,
                s.p,
                s.n,
                &lin,
                &bnodes,
                scale,
            );
            branch += 1;
        }
        for s in &self.vcvs {
            stamp_voltage_branch(&mut m, &mut rhs, dim, nc, branch, s.p, s.n, 0.0, 0.0);
            stamp_vcvs_control(&mut m, dim, nc, branch, s.cp, s.cn, s.gain);
            branch += 1;
        }
        for s in &self.ccvs {
            stamp_voltage_branch(&mut m, &mut rhs, dim, nc, branch, s.p, s.n, 0.0, 0.0);
            let control = self.branch_index(&s.control)?;
            stamp_ccvs_control(&mut m, dim, nc, branch, control, s.gain);
            branch += 1;
        }
        for (k, l) in self.inductors.iter().enumerate() {
            let ib = self.inductor_branch(k);
            let i1 = old[nu + ib];
            let i2 = older.map_or(0.0, |o| o[nu + ib]);
            let (coefficient, history) = match mode {
                Companion::Dc { .. } => (0.0, 0.0),
                Companion::Fixed {
                    method: TransientMethod::BackwardEuler,
                    dt,
                    ..
                } => {
                    let c = l.value / dt;
                    (c, -c * i1)
                }
                Companion::Fixed {
                    method: TransientMethod::Trapezoidal,
                    dt,
                    inductor_voltage,
                    ..
                } => {
                    let c = 2.0 * l.value / dt;
                    (c, -c * i1 - inductor_voltage[k])
                }
                Companion::Fixed {
                    method: TransientMethod::Bdf2,
                    dt,
                    ..
                } => (1.5 * l.value / dt, l.value / (2.0 * dt) * (-4.0 * i1 + i2)),
                Companion::Adaptive { a, .. } => {
                    (l.value * a[0], l.value * (a[1] * i1 + a[2] * i2))
                }
            };
            stamp_voltage_branch(
                &mut m,
                &mut rhs,
                dim,
                nc,
                ib,
                l.p,
                l.n,
                history,
                -coefficient,
            );
        }
        if !matches!(mode, Companion::Dc { .. }) {
            for mi in &self.mutual_inductances {
                let (first, second, mutual) = self.mutual_pair(mi)?;
                let (b1, b2) = (self.inductor_branch(first), self.inductor_branch(second));
                let (r1, r2) = (nu + b1, nu + b2);
                let (coefficient, h1, h2) = match mode {
                    Companion::Fixed {
                        method: TransientMethod::BackwardEuler,
                        dt,
                        ..
                    } => {
                        let c = mutual / dt;
                        (c, -c * old[r2], -c * old[r1])
                    }
                    Companion::Fixed {
                        method: TransientMethod::Trapezoidal,
                        dt,
                        ..
                    } => {
                        let c = 2.0 * mutual / dt;
                        (c, -c * old[r2], -c * old[r1])
                    }
                    Companion::Fixed {
                        method: TransientMethod::Bdf2,
                        dt,
                        ..
                    } => {
                        let o = older.unwrap_or(old);
                        (
                            1.5 * mutual / dt,
                            mutual / (2.0 * dt) * (-4.0 * old[r2] + o[r2]),
                            mutual / (2.0 * dt) * (-4.0 * old[r1] + o[r1]),
                        )
                    }
                    Companion::Adaptive { a, .. } => {
                        let (os, of) = (older.map_or(0.0, |o| o[r2]), older.map_or(0.0, |o| o[r1]));
                        (
                            mutual * a[0],
                            mutual * (a[1] * old[r2] + a[2] * os),
                            mutual * (a[1] * old[r1] + a[2] * of),
                        )
                    }
                    Companion::Dc { .. } => (0.0, 0.0, 0.0),
                };
                *m.at(r1 * dim + r2) -= coefficient;
                *m.at(r2 * dim + r1) -= coefficient;
                rhs[r1] += h1;
                rhs[r2] += h2;
            }
        }
        for s in &self.cccs {
            let control = self.branch_index(&s.control)?;
            stamp_cccs(&mut m, dim, nc, s.p, s.n, control, s.gain);
        }
        for d in &self.diodes {
            let model = model_of(&self.diode_models, &d.model, "diode")?;
            let vt = diode_vt(model);
            let [a, k, ..] = d.nodes;
            let raw = nv(previous, a) - nv(previous, k);
            let vd = if matches!(mode, Companion::Dc { .. }) {
                raw
            } else {
                limit_pn_junction_voltage(
                    raw,
                    nv(old, a) - nv(old, k),
                    vt,
                    model.saturation_current,
                )?
            };
            let ex = (vd / vt).clamp(-40.0, 40.0).exp();
            let current = model.saturation_current * (ex - 1.0);
            let g = model.saturation_current * ex / vt;
            stamp_admittance(&mut m, dim, a, k, g);
            stamp_current(&mut rhs, a, k, current - g * vd);
        }
        for dev in &self.mosfets {
            let model = model_of(&self.mos_models, &dev.model, "MOS")?;
            let [d, g, s, _] = dev.nodes;
            let lin = mos_linearize(model, nv(previous, d), nv(previous, g), nv(previous, s));
            stamp_transistor(&mut m, dim, &lin, d, g, s);
            stamp_current(&mut rhs, d, s, lin.ieq);
        }
        for dev in &self.jfets {
            let model = model_of(&self.jfet_models, &dev.model, "JFET")?;
            let [d, g, s, _] = dev.nodes;
            let lin = jfet_linearize(model, nv(previous, d), nv(previous, g), nv(previous, s));
            stamp_transistor(&mut m, dim, &lin, d, g, s);
            stamp_current(&mut rhs, d, s, lin.ieq);
        }
        for dev in &self.bjts {
            let model = model_of(&self.bjt_models, &dev.model, "BJT")?;
            let [c, b, e, _] = dev.nodes;
            let lin = bjt_linearize(model, nv(previous, c), nv(previous, b), nv(previous, e));
            stamp_bjt_dc(&mut m, &mut rhs, dim, &lin, [c, b, e]);
        }
        for dev in &self.switches {
            let model = model_of(&self.switch_models, &dev.model, "switch")?;
            let [p, n, cp, cn] = dev.nodes;
            let control = nv(previous, cp) - nv(previous, cn);
            stamp_admittance(&mut m, dim, p, n, switch_conductance(model, control));
        }
        if !self.static_devices.is_empty() {
            let nodes = node_values_from_state(previous, nc);
            for dev in &self.static_devices {
                stamp_static_device(&mut m, Some(&mut rhs), dim, dev, &nodes)?;
            }
        }
        if !self.dynamic_devices.is_empty() {
            let current_nodes = node_values_from_state(previous, nc);
            match mode {
                Companion::Dc { .. } => {
                    for dev in &self.dynamic_devices {
                        stamp_dynamic_dc(&mut m, &mut rhs, dim, dev, &current_nodes)?;
                    }
                }
                Companion::Fixed { method, dt, .. } => {
                    let old_nodes = node_values_from_state(old, nc);
                    let older_nodes = older.map(|o| node_values_from_state(o, nc));
                    let bdf2 = method == TransientMethod::Bdf2;
                    let a = if bdf2 {
                        [1.5 / dt, -2.0 / dt, 0.5 / dt]
                    } else {
                        [1.0 / dt, -1.0 / dt, 0.0]
                    };
                    for dev in &self.dynamic_devices {
                        let history = if bdf2 { older_nodes.as_deref() } else { None };
                        stamp_dynamic_transient(
                            &mut m,
                            &mut rhs,
                            dim,
                            dev,
                            &current_nodes,
                            &old_nodes,
                            history,
                            a,
                        )?;
                    }
                }
                Companion::Adaptive { a, .. } => {
                    let old_nodes = node_values_from_state(old, nc);
                    let older_nodes = older.map(|o| node_values_from_state(o, nc));
                    for dev in &self.dynamic_devices {
                        stamp_dynamic_transient(
                            &mut m,
                            &mut rhs,
                            dim,
                            dev,
                            &current_nodes,
                            &old_nodes,
                            older_nodes.as_deref(),
                            a,
                        )?;
                    }
                }
            }
        }
        Ok((m, rhs))
    }

    fn dc_operating_point_scaled(
        &self,
        cfg: &NewtonConfig,
        scale: f64,
        initial: Option<&OperatingPoint>,
    ) -> Result<OperatingPoint> {
        require(
            cfg.max_iterations > 0
                && cfg.voltage_tolerance > 0.0
                && cfg.residual_tolerance > 0.0
                && cfg.gmin >= 0.0
                && scale.is_finite(),
            "invalid circuit Newton controls",
        )?;
        let nc = self.node_count();
        let nu = nc - 1;
        let dim = nu + self.branch_count();
        let mut state = vec![0.0; dim];
        match initial {
            Some(op)
                if op.node_voltage.len() == nc
                    && op.branch_current.len() == self.branch_count() =>
            {
                state[..nc - 1].copy_from_slice(&op.node_voltage[1..nc]);
                state[nu..].copy_from_slice(&op.branch_current);
            }
            _ => {
                for (&node, &v) in &self.nodeset_voltage {
                    if node > 0 && node < nc {
                        state[node - 1] = v;
                    }
                }
            }
        }
        let mut previous = state.clone();
        let nonlinear = self.nonlinear();
        for iteration in 0..cfg.max_iterations {
            let (m, rhs) =
                self.assemble_real(cfg, Companion::Dc { scale }, &previous, &previous, None)?;
            let mut candidate = solve_real_mna(&m, rhs, cfg)?;
            for d in &self.diodes {
                let model = model_of(&self.diode_models, &d.model, "diode")?;
                let vt = diode_vt(model);
                let [a, k, ..] = d.nodes;
                let old_v = nv(&previous, a) - nv(&previous, k);
                let proposed = nv(&candidate, a) - nv(&candidate, k);
                let limited =
                    limit_pn_junction_voltage(proposed, old_v, vt, model.saturation_current)?;
                let correction = limited - proposed;
                if a > 0 && k > 0 {
                    candidate[a - 1] += 0.5 * correction;
                    candidate[k - 1] -= 0.5 * correction;
                } else if a > 0 {
                    candidate[a - 1] += correction;
                } else if k > 0 {
                    candidate[k - 1] -= correction;
                }
            }
            let mut delta = 0.0f64;
            for i in 0..nu {
                delta = delta.max((candidate[i] - previous[i]).abs());
            }
            let damping = if nonlinear && delta > 0.5 {
                0.5 / delta
            } else {
                1.0
            };
            for i in 0..dim {
                state[i] = previous[i] + damping * (candidate[i] - previous[i]);
            }
            previous.clone_from(&state);
            if delta * damping < cfg.voltage_tolerance {
                let (node_voltage, branch_current) = self.state_to_point(&state);
                return Ok(OperatingPoint {
                    node_voltage,
                    branch_current,
                    iterations: iteration + 1,
                    converged: true,
                });
            }
        }
        let (node_voltage, branch_current) = self.state_to_point(&previous);
        Ok(OperatingPoint {
            node_voltage,
            branch_current,
            iterations: cfg.max_iterations,
            converged: false,
        })
    }

    /// Plain damped Newton from the nodeset (or zero) guess; check `converged`.
    pub fn dc_operating_point(&self, cfg: &NewtonConfig) -> Result<OperatingPoint> {
        self.dc_operating_point_scaled(cfg, 1.0, None)
    }

    /// Source stepping with gmin ≥ 1e-3, then geometric gmin stepping down to `max(gmin,1e-15)`,
    /// then a final solve with the caller's controls. Returns early (unconverged) on failure.
    pub fn dc_operating_point_homotopy(
        &self,
        cfg: &NewtonConfig,
        source_steps: usize,
        gmin_steps: usize,
    ) -> Result<OperatingPoint> {
        require(
            source_steps > 0 && gmin_steps > 0,
            "homotopy step counts must be positive",
        )?;
        let mut staged = *cfg;
        staged.gmin = cfg.gmin.max(1.0e-3);
        let mut current = OperatingPoint::default();
        for step in 0..source_steps {
            let scale = (step + 1) as f64 / source_steps as f64;
            let init = if current.node_voltage.is_empty() {
                None
            } else {
                Some(&current)
            };
            current = self.dc_operating_point_scaled(&staged, scale, init)?;
            if !current.converged {
                return Ok(current);
            }
        }
        let target = cfg.gmin.max(1.0e-15);
        for step in 0..gmin_steps {
            let f = (step + 1) as f64 / gmin_steps as f64;
            staged.gmin = (1.0e-3f64.ln() * (1.0 - f) + target.ln() * f).exp();
            current = self.dc_operating_point_scaled(&staged, 1.0, Some(&current))?;
            if !current.converged {
                return Ok(current);
            }
        }
        self.dc_operating_point_scaled(cfg, 1.0, Some(&current))
    }

    /// Upstream default homotopy (8 source steps, 6 gmin steps).
    pub fn dc_operating_point_default_homotopy(
        &self,
        cfg: &NewtonConfig,
    ) -> Result<OperatingPoint> {
        self.dc_operating_point_homotopy(cfg, 8, 6)
    }

    // ------------------------------------------------------------------------------------------
    // Complex (AC / noise) assembly
    // ------------------------------------------------------------------------------------------

    fn assemble_ac(
        &self,
        frequency: f64,
        op: &OperatingPoint,
        noise: bool,
    ) -> Result<(Vec<Complex>, Vec<Complex>)> {
        let nc = self.node_count();
        let nu = nc - 1;
        let dim = nu + self.branch_count();
        let mut m = vec![Complex::default(); dim * dim];
        let mut rhs = vec![Complex::default(); dim];
        let omega = 2.0 * PI * frequency;
        let jw = Complex::new(0.0, omega);
        let nodes = &op.node_voltage;
        require(nodes.len() == nc, "operating point node size mismatch")?;
        for r in &self.resistors {
            stamp_admittance(&mut m, dim, r.p, r.n, Complex::new(1.0 / r.value, 0.0));
        }
        for c in &self.capacitors {
            stamp_admittance(&mut m, dim, c.p, c.n, jw * c.value);
        }
        if !noise {
            for s in &self.current_sources {
                stamp_current(&mut rhs, s.p, s.n, s.ac);
            }
        }
        for s in &self.vccs {
            stamp_vccs(&mut m, dim, s, Complex::new(s.gain, 0.0));
        }
        for s in &self.behavioral_current_sources {
            let lin = behavioral_linearize(&s.evaluator, nodes, 0.0)?;
            stamp_behavioral_current(&mut m, None, dim, s.p, s.n, &lin, nodes, 1.0);
        }
        let zero = Complex::default();
        let mut branch = 0;
        for s in &self.voltage_sources {
            let v = if noise { zero } else { s.ac };
            stamp_voltage_branch(&mut m, &mut rhs, dim, nc, branch, s.p, s.n, v, zero);
            branch += 1;
        }
        for s in &self.behavioral_voltage_sources {
            let lin = behavioral_linearize(&s.evaluator, nodes, 0.0)?;
            stamp_behavioral_voltage(&mut m, None, dim, nc, branch, s.p, s.n, &lin, nodes, 1.0);
            branch += 1;
        }
        for s in &self.vcvs {
            stamp_voltage_branch(&mut m, &mut rhs, dim, nc, branch, s.p, s.n, zero, zero);
            stamp_vcvs_control(
                &mut m,
                dim,
                nc,
                branch,
                s.cp,
                s.cn,
                Complex::new(s.gain, 0.0),
            );
            branch += 1;
        }
        for s in &self.ccvs {
            stamp_voltage_branch(&mut m, &mut rhs, dim, nc, branch, s.p, s.n, zero, zero);
            let control = self.branch_index(&s.control)?;
            stamp_ccvs_control(&mut m, dim, nc, branch, control, Complex::new(s.gain, 0.0));
            branch += 1;
        }
        for l in &self.inductors {
            stamp_voltage_branch(
                &mut m,
                &mut rhs,
                dim,
                nc,
                branch,
                l.p,
                l.n,
                zero,
                -jw * l.value,
            );
            branch += 1;
        }
        if !noise {
            for mi in &self.mutual_inductances {
                let (first, second, mutual) = self.mutual_pair(mi)?;
                let (r1, r2) = (
                    nu + self.inductor_branch(first),
                    nu + self.inductor_branch(second),
                );
                m[r1 * dim + r2] += -jw * mutual;
                m[r2 * dim + r1] += -jw * mutual;
            }
        }
        for s in &self.cccs {
            let control = self.branch_index(&s.control)?;
            stamp_cccs(
                &mut m,
                dim,
                nc,
                s.p,
                s.n,
                control,
                Complex::new(s.gain, 0.0),
            );
        }
        for d in &self.diodes {
            let model = model_of(&self.diode_models, &d.model, "diode")?;
            let vt = diode_vt(model);
            let [a, k, ..] = d.nodes;
            let vd = nodes[a] - nodes[k];
            let g = model.saturation_current * (vd / vt).clamp(-40.0, 40.0).exp() / vt;
            stamp_admittance(&mut m, dim, a, k, Complex::new(g, 0.0));
        }
        for dev in &self.mosfets {
            let [d, g, s, _] = dev.nodes;
            let lin = mos_linearize(
                model_of(&self.mos_models, &dev.model, "MOS")?,
                nodes[d],
                nodes[g],
                nodes[s],
            );
            stamp_transistor(&mut m, dim, &lin, d, g, s);
        }
        for dev in &self.jfets {
            let [d, g, s, _] = dev.nodes;
            let lin = jfet_linearize(
                model_of(&self.jfet_models, &dev.model, "JFET")?,
                nodes[d],
                nodes[g],
                nodes[s],
            );
            stamp_transistor(&mut m, dim, &lin, d, g, s);
        }
        for dev in &self.bjts {
            let [c, b, e, _] = dev.nodes;
            let lin = bjt_linearize(
                model_of(&self.bjt_models, &dev.model, "BJT")?,
                nodes[c],
                nodes[b],
                nodes[e],
            );
            stamp_bjt_jacobian(&mut m, dim, &lin, [c, b, e]);
        }
        for dev in &self.switches {
            let [p, n, cp, cn] = dev.nodes;
            let g = switch_conductance(
                model_of(&self.switch_models, &dev.model, "switch")?,
                nodes[cp] - nodes[cn],
            );
            stamp_admittance(&mut m, dim, p, n, Complex::new(g, 0.0));
        }
        if !noise {
            for dev in &self.static_devices {
                stamp_static_device(&mut m, None, dim, dev, nodes)?;
            }
        }
        for dev in &self.dynamic_devices {
            stamp_dynamic_ac(&mut m, dim, dev, nodes, omega)?;
        }
        for dev in &self.nports {
            stamp_nport(&mut m, dim, dev, frequency)?;
        }
        Ok((m, rhs))
    }

    /// Small-signal AC solve at `frequency`. Without an operating point, nonlinear circuits are
    /// biased with the default homotopy; linear circuits use a zero bias.
    pub fn ac(
        &self,
        frequency: f64,
        operating_point: Option<&OperatingPoint>,
        cfg: &NewtonConfig,
    ) -> Result<AcSolution> {
        require(
            frequency > 0.0 && frequency.is_finite(),
            "AC frequency must be positive",
        )?;
        let owned;
        let op = match operating_point {
            Some(op) => op,
            None => {
                owned = if self.nonlinear() {
                    let op = self.dc_operating_point_default_homotopy(cfg)?;
                    if !op.converged {
                        return Err(err_runtime("AC operating point did not converge"));
                    }
                    op
                } else {
                    OperatingPoint {
                        node_voltage: vec![0.0; self.node_count()],
                        branch_current: vec![0.0; self.branch_count()],
                        iterations: 0,
                        converged: true,
                    }
                };
                &owned
            }
        };
        let (m, rhs) = self.assemble_ac(frequency, op, false)?;
        let dim = rhs.len();
        let solution = solve_dense(m, rhs, dim, "singular circuit matrix")?;
        let nu = self.node_count() - 1;
        let mut node_voltage = vec![Complex::default(); self.node_count()];
        node_voltage[1..].copy_from_slice(&solution[..nu]);
        Ok(AcSolution {
            frequency_hz: frequency,
            node_voltage,
            branch_current: solution[nu..].to_vec(),
        })
    }

    /// Logarithmic AC sweep around one homotopy operating point (last point exactly `stop`).
    pub fn ac_log_sweep(
        &self,
        start: f64,
        stop: f64,
        points: usize,
        cfg: &NewtonConfig,
    ) -> Result<Vec<AcSolution>> {
        require(
            start > 0.0 && stop >= start && points > 0,
            "invalid AC sweep range",
        )?;
        let op = self.dc_operating_point_default_homotopy(cfg)?;
        if !op.converged {
            return Err(err_runtime("AC sweep operating point did not converge"));
        }
        if points == 1 {
            return Ok(vec![self.ac(start, Some(&op), cfg)?]);
        }
        let ratio = (stop / start).powf(1.0 / (points - 1) as f64);
        let mut f = start;
        let mut out = Vec::with_capacity(points);
        for i in 0..points {
            out.push(self.ac(if i + 1 == points { stop } else { f }, Some(&op), cfg)?);
            f *= ratio;
        }
        Ok(out)
    }

    /// Sweep one voltage source's DC value (restored afterwards) with homotopy at each point.
    pub fn dc_sweep_voltage_source(
        &mut self,
        source: &str,
        start: f64,
        stop: f64,
        points: usize,
        cfg: &NewtonConfig,
    ) -> Result<Vec<DcSweepPoint>> {
        require(
            points > 0 && start.is_finite() && stop.is_finite(),
            "invalid DC sweep range",
        )?;
        let target = upper(source);
        let index = self
            .voltage_sources
            .iter()
            .position(|s| upper(&s.name) == target)
            .ok_or_else(|| {
                SolverError::InvalidInput("unknown voltage source for DC sweep".into())
            })?;
        let original = self.voltage_sources[index].dc;
        let mut out = Vec::with_capacity(points);
        let mut outcome = Ok(());
        for i in 0..points {
            let t = if points == 1 {
                0.0
            } else {
                i as f64 / (points - 1) as f64
            };
            let value = start + (stop - start) * t;
            self.voltage_sources[index].dc = value;
            match self.dc_operating_point_default_homotopy(cfg) {
                Ok(op) => out.push(DcSweepPoint {
                    source_value: value,
                    operating_point: op,
                }),
                Err(e) => {
                    outcome = Err(e);
                    break;
                }
            }
        }
        self.voltage_sources[index].dc = original;
        outcome.map(|_| out)
    }

    /// Re-solve with every diode/BJT model at each temperature (on a copy).
    pub fn dc_temperature_sweep(
        &self,
        start: f64,
        stop: f64,
        points: usize,
        cfg: &NewtonConfig,
    ) -> Result<Vec<TemperatureSweepPoint>> {
        require(
            points > 0 && start > 0.0 && stop > 0.0 && start.is_finite() && stop.is_finite(),
            "invalid temperature sweep range",
        )?;
        let mut out = Vec::with_capacity(points);
        for i in 0..points {
            let f = if points == 1 {
                0.0
            } else {
                i as f64 / (points - 1) as f64
            };
            let temperature = start + (stop - start) * f;
            let mut copy = self.clone();
            copy.set_device_temperature(temperature)?;
            out.push(TemperatureSweepPoint {
                temperature_k: temperature,
                operating_point: copy.dc_operating_point_default_homotopy(cfg)?,
            });
        }
        Ok(out)
    }

    fn initial_state(&self) -> Vec<f64> {
        let nc = self.node_count();
        let mut state = vec![0.0; nc - 1 + self.branch_count()];
        for (&node, &v) in &self.initial_voltage {
            if node > 0 && node < nc {
                state[node - 1] = v;
            }
        }
        state
    }

    fn transient_point(&self, time: f64, state: &[f64]) -> TransientPoint {
        let (node_voltage, branch_current) = self.state_to_point(state);
        TransientPoint {
            time,
            node_voltage,
            branch_current,
        }
    }

    /// Fixed-step transient from the `.IC` state (no DC solve). BDF2 bootstraps with one
    /// backward-Euler step.
    pub fn transient(
        &self,
        dt: f64,
        steps: usize,
        method: TransientMethod,
        cfg: &NewtonConfig,
    ) -> Result<Vec<TransientPoint>> {
        require(
            dt > 0.0 && dt.is_finite(),
            "transient time step must be positive",
        )?;
        require(
            !(method == TransientMethod::Trapezoidal && !self.dynamic_devices.is_empty()),
            "charge-based dynamic devices currently support backward Euler or BDF2 transient integration",
        )?;
        let nu = self.node_count() - 1;
        let dim = nu + self.branch_count();
        let mut state = self.initial_state();
        let mut older = state.clone();
        let mut capacitor_current = vec![0.0; self.capacitors.len()];
        let mut inductor_voltage = vec![0.0; self.inductors.len()];
        let mut out = Vec::with_capacity(steps + 1);
        out.push(self.transient_point(0.0, &state));
        for step in 0..steps {
            let next_time = (step + 1) as f64 * dt;
            let old = state.clone();
            let mut previous = state.clone();
            let active = if method == TransientMethod::Bdf2 && step == 0 {
                TransientMethod::BackwardEuler
            } else {
                method
            };
            let mut converged = false;
            for _ in 0..cfg.max_iterations {
                let mode = Companion::Fixed {
                    method: active,
                    dt,
                    time: next_time,
                    capacitor_current: &capacitor_current,
                    inductor_voltage: &inductor_voltage,
                };
                let (m, rhs) = self.assemble_real(cfg, mode, &previous, &old, Some(&older))?;
                let candidate = solve_real_mna(&m, rhs, cfg)?;
                let mut delta = 0.0f64;
                for i in 0..dim {
                    delta = delta.max((candidate[i] - previous[i]).abs());
                }
                let damping = if delta > 0.5 { 0.5 / delta } else { 1.0 };
                for i in 0..dim {
                    state[i] = previous[i] + damping * (candidate[i] - previous[i]);
                }
                previous.clone_from(&state);
                if delta * damping < cfg.voltage_tolerance {
                    converged = true;
                    break;
                }
            }
            if !converged {
                return Err(err_runtime(
                    "transient circuit Newton solve did not converge",
                ));
            }
            if active == TransientMethod::Trapezoidal {
                for (k, c) in self.capacitors.iter().enumerate() {
                    let v0 = nv(&old, c.p) - nv(&old, c.n);
                    let v1 = nv(&state, c.p) - nv(&state, c.n);
                    capacitor_current[k] = 2.0 * c.value / dt * (v1 - v0) - capacitor_current[k];
                }
                for (k, l) in self.inductors.iter().enumerate() {
                    inductor_voltage[k] = nv(&state, l.p) - nv(&state, l.n);
                }
            }
            older = old;
            out.push(self.transient_point(next_time, &state));
        }
        Ok(out)
    }

    #[allow(clippy::too_many_arguments)]
    fn adaptive_advance(
        &self,
        cfg: &NewtonConfig,
        old: &[f64],
        older: Option<&[f64]>,
        previous_dt: f64,
        next_time: f64,
        dt: f64,
        active: TransientMethod,
    ) -> Result<Vec<f64>> {
        let mut a = [1.0 / dt, -1.0 / dt, 0.0];
        if active == TransientMethod::Bdf2 {
            require(
                older.is_some() && previous_dt > 0.0 && previous_dt.is_finite(),
                "BDF2 adaptive step requires valid history",
            )?;
            let sum = dt + previous_dt;
            a = [
                (2.0 * dt + previous_dt) / (dt * sum),
                -sum / (dt * previous_dt),
                dt / (previous_dt * sum),
            ];
        }
        let nu = self.node_count() - 1;
        let dim = old.len();
        let nonlinear = self.nonlinear();
        let mut current = old.to_vec();
        let mut previous = old.to_vec();
        for _ in 0..cfg.max_iterations {
            let mode = Companion::Adaptive { time: next_time, a };
            let (m, rhs) = self.assemble_real(cfg, mode, &previous, old, older)?;
            let candidate = solve_real_mna(&m, rhs, cfg)?;
            let mut delta = 0.0f64;
            for i in 0..nu {
                delta = delta.max((candidate[i] - previous[i]).abs());
            }
            let damping = if nonlinear && delta > 0.5 {
                0.5 / delta
            } else {
                1.0
            };
            for i in 0..dim {
                current[i] = previous[i] + damping * (candidate[i] - previous[i]);
            }
            previous.clone_from(&current);
            if delta * damping < cfg.voltage_tolerance {
                return Ok(current);
            }
        }
        Err(err_runtime(
            "adaptive transient Newton solve did not converge",
        ))
    }

    /// Variable-step backward Euler or BDF2 with step-doubling (Richardson) LTE control.
    pub fn transient_adaptive(
        &self,
        stop_time: f64,
        controls: &AdaptiveTransientConfig,
        cfg: &NewtonConfig,
    ) -> Result<AdaptiveTransientResult> {
        let c = controls;
        require(
            stop_time > 0.0
                && stop_time.is_finite()
                && c.initial_step > 0.0
                && c.minimum_step > 0.0
                && c.maximum_step >= c.minimum_step
                && c.relative_tolerance > 0.0
                && c.absolute_tolerance > 0.0
                && c.safety_factor > 0.0
                && c.minimum_scale > 0.0
                && c.maximum_scale >= c.minimum_scale
                && c.max_step_attempts > 0,
            "invalid adaptive transient controls",
        )?;
        require(
            c.method != TransientMethod::Trapezoidal,
            "adaptive transient currently supports backward Euler or BDF2",
        )?;
        let nu = self.node_count() - 1;
        let dim = nu + self.branch_count();
        let mut state = self.initial_state();
        let mut older_state = state.clone();
        let mut result = AdaptiveTransientResult {
            minimum_accepted_step: f64::INFINITY,
            ..Default::default()
        };
        result.points.push(self.transient_point(0.0, &state));
        let (mut time, mut previous_dt, mut has_history) = (0.0, 0.0, false);
        let mut dt = c.initial_step.max(c.minimum_step).min(c.maximum_step);
        let mut attempts = 0;
        while time < stop_time {
            attempts += 1;
            if attempts > c.max_step_attempts {
                return Err(err_runtime(
                    "adaptive transient step-attempt budget exhausted",
                ));
            }
            let remaining = stop_time - time;
            dt = dt.min(remaining);
            let use_bdf2 = c.method == TransientMethod::Bdf2 && has_history;
            let active = if use_bdf2 {
                TransientMethod::Bdf2
            } else {
                TransientMethod::BackwardEuler
            };
            let history = if use_bdf2 {
                Some(older_state.as_slice())
            } else {
                None
            };
            let history_dt = if use_bdf2 { previous_dt } else { 0.0 };
            let full =
                self.adaptive_advance(cfg, &state, history, history_dt, time + dt, dt, active)?;
            let half = self.adaptive_advance(
                cfg,
                &state,
                history,
                history_dt,
                time + 0.5 * dt,
                0.5 * dt,
                active,
            )?;
            let refined = self.adaptive_advance(
                cfg,
                &half,
                if use_bdf2 {
                    Some(state.as_slice())
                } else {
                    None
                },
                if use_bdf2 { 0.5 * dt } else { 0.0 },
                time + dt,
                0.5 * dt,
                active,
            )?;
            let richardson = if use_bdf2 { 3.0 } else { 1.0 };
            let unknowns = if nu > 0 { nu } else { dim };
            let mut error = 0.0f64;
            for i in 0..unknowns {
                let scale = c.absolute_tolerance
                    + c.relative_tolerance * full[i].abs().max(refined[i].abs());
                error = error.max((refined[i] - full[i]).abs() / (richardson * scale));
            }
            if !error.is_finite() {
                return Err(err_runtime(
                    "adaptive transient produced non-finite LTE estimate",
                ));
            }
            if error <= 1.0 {
                older_state = std::mem::replace(&mut state, refined);
                previous_dt = dt;
                has_history = true;
                time += dt;
                result.accepted_steps += 1;
                result.minimum_accepted_step = result.minimum_accepted_step.min(dt);
                result.maximum_accepted_step = result.maximum_accepted_step.max(dt);
                result.points.push(self.transient_point(time, &state));
            } else {
                result.rejected_steps += 1;
                if dt <= c.minimum_step * (1.0 + 16.0 * f64::EPSILON) && remaining > c.minimum_step
                {
                    return Err(err_runtime(
                        "adaptive transient LTE tolerance cannot be met at minimum step",
                    ));
                }
            }
            let exponent = if use_bdf2 { 1.0 / 3.0 } else { 0.5 };
            let raw = if error > 0.0 {
                c.safety_factor * (1.0 / error).powf(exponent)
            } else {
                c.maximum_scale
            };
            let factor = raw.max(c.minimum_scale).min(c.maximum_scale);
            dt = (dt * factor).max(c.minimum_step).min(c.maximum_step);
        }
        if result.accepted_steps == 0 {
            result.minimum_accepted_step = 0.0;
        }
        Ok(result)
    }

    /// Output noise PSD at `output_node` from resistor thermal, diode/BJT shot and MOS/JFET
    /// channel noise, each propagated by a unit-current transfer solve.
    pub fn output_noise(
        &self,
        frequency: f64,
        output_node: &str,
        temperature: f64,
        cfg: &NewtonConfig,
    ) -> Result<NoiseResult> {
        require(
            frequency > 0.0 && temperature > 0.0,
            "invalid noise analysis controls",
        )?;
        let op = self.dc_operating_point_default_homotopy(cfg)?;
        if !op.converged {
            return Err(err_runtime("noise operating point did not converge"));
        }
        let output = self.find_node(output_node)?;
        let (m, _) = self.assemble_ac(frequency, &op, true)?;
        let dim = self.node_count() - 1 + self.branch_count();
        let transfer = |p: usize, n: usize| -> Result<Complex> {
            let mut rhs = vec![Complex::default(); dim];
            stamp_current(&mut rhs, p, n, Complex::new(1.0, 0.0));
            let sol = solve_dense(m.clone(), rhs, dim, "singular circuit matrix")?;
            Ok(if output > 0 {
                sol[output - 1]
            } else {
                Complex::default()
            })
        };
        let mut result = NoiseResult {
            frequency_hz: frequency,
            ..Default::default()
        };
        let mut add = |name: String, p: usize, n: usize, density: f64| -> Result<()> {
            if !(density >= 0.0 && density.is_finite()) {
                return Err(err_runtime("non-finite device noise density"));
            }
            let contribution = transfer(p, n)?.norm_sqr() * density;
            result.output_noise_density_v2_per_hz += contribution;
            result.contributions.push(NoiseContribution {
                source: name,
                output_noise_density_v2_per_hz: contribution,
            });
            Ok(())
        };
        let nodes = &op.node_voltage;
        for r in &self.resistors {
            add(
                r.name.clone(),
                r.p,
                r.n,
                4.0 * BOLTZMANN * temperature / r.value,
            )?;
        }
        for d in &self.diodes {
            let model = model_of(&self.diode_models, &d.model, "diode")?;
            let vt = diode_vt(model);
            let [a, k, ..] = d.nodes;
            let vd = nodes[a] - nodes[k];
            let current = model.saturation_current * ((vd / vt).clamp(-40.0, 40.0).exp() - 1.0);
            add(d.name.clone(), a, k, 2.0 * ELECTRON_CHARGE * current.abs())?;
        }
        for dev in &self.mosfets {
            let [d, g, s, _] = dev.nodes;
            let lin = mos_linearize(
                model_of(&self.mos_models, &dev.model, "MOS")?,
                nodes[d],
                nodes[g],
                nodes[s],
            );
            add(
                dev.name.clone(),
                d,
                s,
                4.0 * BOLTZMANN * temperature * (2.0 / 3.0) * lin.gm.max(0.0),
            )?;
        }
        for dev in &self.jfets {
            let [d, g, s, _] = dev.nodes;
            let lin = jfet_linearize(
                model_of(&self.jfet_models, &dev.model, "JFET")?,
                nodes[d],
                nodes[g],
                nodes[s],
            );
            add(
                dev.name.clone(),
                d,
                s,
                4.0 * BOLTZMANN * temperature * (2.0 / 3.0) * lin.gm.max(0.0),
            )?;
        }
        for dev in &self.bjts {
            let [c, b, e, _] = dev.nodes;
            let currents = bjt_currents(
                model_of(&self.bjt_models, &dev.model, "BJT")?,
                nodes[c],
                nodes[b],
                nodes[e],
            );
            add(
                format!("{}:collector", dev.name),
                c,
                e,
                2.0 * ELECTRON_CHARGE * currents[0].abs(),
            )?;
            add(
                format!("{}:base", dev.name),
                b,
                e,
                2.0 * ELECTRON_CHARGE * currents[1].abs(),
            )?;
        }
        result.output_noise_v_per_sqrt_hz = result.output_noise_density_v2_per_hz.max(0.0).sqrt();
        Ok(result)
    }

    pub fn voltage(&self, point: &OperatingPoint, node: &str) -> Result<f64> {
        let id = self.find_node(node)?;
        require(
            point.node_voltage.len() == self.node_count(),
            "operating point node size mismatch",
        )?;
        Ok(point.node_voltage[id])
    }
    pub fn ac_voltage(&self, point: &AcSolution, node: &str) -> Result<Complex> {
        let id = self.find_node(node)?;
        require(
            point.node_voltage.len() == self.node_count(),
            "AC point node size mismatch",
        )?;
        Ok(point.node_voltage[id])
    }
    pub fn branch_current(&self, point: &AcSolution, branch: &str) -> Result<Complex> {
        require(
            point.branch_current.len() == self.branch_count(),
            "AC point branch-current size mismatch",
        )?;
        Ok(point.branch_current[self.branch_index(branch)?])
    }
}

fn first_sample(frequencies: &[f64], length: f64, what: &str) -> Result<f64> {
    match frequencies.first() {
        Some(&f) if length >= 0.0 && length.is_finite() => Ok(f),
        _ => Err(SolverError::InvalidInput(format!(
            "{what} circuit element requires samples and non-negative length"
        ))),
    }
}

fn two_port_point(frequency: f64, z0: f64, two: rf::Matrix2C) -> NPortPoint {
    let mut s = ComplexMatrix::new(2);
    s[(0, 0)] = two.a11;
    s[(0, 1)] = two.a12;
    s[(1, 0)] = two.a21;
    s[(1, 1)] = two.a22;
    NPortPoint {
        frequency_hz: frequency,
        s,
        reference_impedance: vec![z0, z0],
    }
}

fn electrothermal_evaluation(
    local: &[f64],
    ne: usize,
    config: &ElectroThermalDeviceConfig,
    evaluator: &ElectroThermalDeviceEvaluator,
) -> Result<DynamicDeviceEvaluation> {
    if local.len() != ne + 1 {
        return Err(err_runtime("electrothermal local state size mismatch"));
    }
    let theta = local[ne];
    let temperature = config.ambient_temperature_k + theta;
    if !(temperature > 0.0 && temperature.is_finite()) {
        return Err(err_runtime(
            "electrothermal device reached non-physical temperature",
        ));
    }
    let e = evaluator(&local[..ne], temperature)?;
    if e.electrical_current.len() != ne
        || e.electrical_jacobian.len() != ne * ne
        || e.electrical_temperature_derivative.len() != ne
        || e.power_voltage_derivative.len() != ne
        || !e.dissipated_power_w.is_finite()
        || !e.power_temperature_derivative_w_per_k.is_finite()
    {
        return Err(err_runtime(
            "electrothermal evaluator returned invalid dimensions/values",
        ));
    }
    let total = ne + 1;
    let mut r = DynamicDeviceEvaluation {
        terminal_current: vec![0.0; total],
        current_jacobian: vec![0.0; total * total],
        terminal_charge: vec![0.0; total],
        charge_jacobian: vec![0.0; total * total],
    };
    for row in 0..ne {
        if !e.electrical_current[row].is_finite()
            || !e.electrical_temperature_derivative[row].is_finite()
        {
            return Err(err_runtime(
                "electrothermal electrical current is non-finite",
            ));
        }
        r.terminal_current[row] = e.electrical_current[row];
        for col in 0..ne {
            let g = e.electrical_jacobian[row * ne + col];
            if !g.is_finite() {
                return Err(err_runtime(
                    "electrothermal electrical Jacobian is non-finite",
                ));
            }
            r.current_jacobian[row * total + col] = g;
        }
        r.current_jacobian[row * total + ne] = e.electrical_temperature_derivative[row];
    }
    let th = ne;
    r.terminal_current[th] = theta / config.thermal_resistance_k_per_w - e.dissipated_power_w;
    for col in 0..ne {
        let d = e.power_voltage_derivative[col];
        if !d.is_finite() {
            return Err(err_runtime("electrothermal power Jacobian is non-finite"));
        }
        r.current_jacobian[th * total + col] = -d;
    }
    r.current_jacobian[th * total + th] =
        1.0 / config.thermal_resistance_k_per_w - e.power_temperature_derivative_w_per_k;
    r.terminal_charge[th] = config.thermal_capacitance_j_per_k * theta;
    r.charge_jacobian[th * total + th] = config.thermal_capacitance_j_per_k;
    Ok(r)
}
