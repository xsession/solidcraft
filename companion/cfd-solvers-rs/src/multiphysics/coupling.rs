//! Partitioned fixed-point coupling with Aitken relaxation and EM/circuit port iteration.
use crate::{require, Complex64, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PartitionedCouplerConfig {
    pub max_iterations: usize,
    pub absolute_tolerance: f64,
    pub initial_relaxation: f64,
    pub aitken: bool,
    pub minimum_relaxation: f64,
    pub maximum_relaxation: f64,
}

impl Default for PartitionedCouplerConfig {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            absolute_tolerance: 1e-10,
            initial_relaxation: 1.0,
            aitken: true,
            minimum_relaxation: 0.05,
            maximum_relaxation: 1.5,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PartitionedCouplerResult {
    pub iterations: usize,
    pub residual_rms: f64,
    pub relaxation: f64,
    pub converged: bool,
}

/// Solve `state = update(state)` by relaxed fixed-point iteration. `update` writes its
/// candidate into the second slice. Convergence is on the RMS of `candidate - state`.
pub fn solve_partitioned_fixed_point(
    state: &mut [f64],
    mut update: impl FnMut(&[f64], &mut [f64]) -> Result<()>,
    config: PartitionedCouplerConfig,
) -> Result<PartitionedCouplerResult> {
    require(!state.is_empty(), "partitioned coupling requires state")?;
    require(
        config.max_iterations > 0
            && config.absolute_tolerance > 0.0
            && config.minimum_relaxation > 0.0
            && config.maximum_relaxation >= config.minimum_relaxation
            && config.initial_relaxation >= config.minimum_relaxation
            && config.initial_relaxation <= config.maximum_relaxation,
        "invalid partitioned coupling controls",
    )?;
    let n = state.len();
    let mut candidate = vec![0.0; n];
    let mut residual = vec![0.0; n];
    let mut previous = vec![0.0; n];
    let mut omega = config.initial_relaxation;
    let mut rms = 0.0;
    for iteration in 0..config.max_iterations {
        update(state, &mut candidate)?;
        let mut rr = 0.0;
        for i in 0..n {
            residual[i] = candidate[i] - state[i];
            rr += residual[i] * residual[i];
        }
        rms = (rr / n as f64).sqrt();
        if !rms.is_finite() {
            return Ok(PartitionedCouplerResult {
                iterations: iteration,
                residual_rms: rms,
                relaxation: omega,
                converged: false,
            });
        }
        if rms <= config.absolute_tolerance {
            return Ok(PartitionedCouplerResult {
                iterations: iteration,
                residual_rms: rms,
                relaxation: omega,
                converged: true,
            });
        }
        if config.aitken && iteration > 0 {
            let (mut num, mut den) = (0.0, 0.0);
            for i in 0..n {
                let delta = residual[i] - previous[i];
                num += previous[i] * delta;
                den += delta * delta;
            }
            if den > 1e-30 {
                omega = (-omega * num / den)
                    .clamp(config.minimum_relaxation, config.maximum_relaxation);
            }
        }
        for (x, r) in state.iter_mut().zip(&residual) {
            *x += omega * r;
        }
        previous.copy_from_slice(&residual);
    }
    Ok(PartitionedCouplerResult {
        iterations: config.max_iterations,
        residual_rms: rms,
        relaxation: omega,
        converged: false,
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EmSpiceCouplingConfig {
    pub max_iterations: usize,
    pub relaxation: f64,
    pub tolerance: f64,
}

impl Default for EmSpiceCouplingConfig {
    fn default() -> Self {
        Self {
            max_iterations: 50,
            relaxation: 0.7,
            tolerance: 1e-10,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EmSpiceCouplingResult {
    pub port_voltage: Complex64,
    pub port_current: Complex64,
    pub iterations: usize,
    pub converged: bool,
}

/// Iterate a field solve (current -> port voltage) against a circuit solve
/// (voltage -> port current) with under-relaxation of the current.
pub fn iterate_em_spice_port(
    mut current: Complex64,
    field_solve: impl Fn(Complex64) -> Complex64,
    circuit_solve: impl Fn(Complex64) -> Complex64,
    config: EmSpiceCouplingConfig,
) -> Result<EmSpiceCouplingResult> {
    require(
        config.max_iterations > 0
            && config.relaxation > 0.0
            && config.relaxation <= 1.0
            && config.tolerance.is_finite()
            && config.tolerance > 0.0,
        "invalid EM/SPICE coupling controls",
    )?;
    let mut result = EmSpiceCouplingResult {
        port_voltage: Complex64::new(0.0, 0.0),
        port_current: current,
        iterations: 0,
        converged: false,
    };
    for iteration in 0..config.max_iterations {
        result.port_voltage = field_solve(current);
        let proposed = circuit_solve(result.port_voltage);
        let updated = current * (1.0 - config.relaxation) + proposed * config.relaxation;
        result.iterations = iteration + 1;
        if (updated - current).norm() <= config.tolerance * updated.norm().max(1.0) {
            result.port_current = updated;
            result.converged = true;
            return Ok(result);
        }
        current = updated;
    }
    result.port_current = current;
    Ok(result)
}
