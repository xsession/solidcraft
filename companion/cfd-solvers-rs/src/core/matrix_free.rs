use super::{
    l2_norm, restarted_gmres, rms, IdentityPreconditioner, IterativeSolverResult, LinearConfig,
};
use crate::{finite_slice, positive, require, Result};
use serde::{Deserialize, Serialize};

pub type ResidualAction<'a> = dyn Fn(&[f64], &mut [f64]) -> Result<()> + 'a;
pub type JacobianAction<'a> = dyn Fn(&[f64], &[f64], &mut [f64]) -> Result<()> + 'a;
pub struct MatrixFreeNonlinearSystem<'a> {
    pub size: usize,
    pub residual: &'a ResidualAction<'a>,
    pub jacobian_action: Option<&'a JacobianAction<'a>>,
    pub transpose_jacobian_action: Option<&'a JacobianAction<'a>>,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct MatrixFreeNewtonConfig {
    pub max_iterations: usize,
    pub linear_max_iterations: usize,
    pub linear_restart: usize,
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub linear_relative_tolerance: f64,
    pub line_search_steps: usize,
    pub minimum_step: f64,
    pub finite_difference_step: f64,
}
impl Default for MatrixFreeNewtonConfig {
    fn default() -> Self {
        Self {
            max_iterations: 30,
            linear_max_iterations: 200,
            linear_restart: 30,
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-8,
            linear_relative_tolerance: 1e-8,
            line_search_steps: 12,
            minimum_step: 1e-7,
            finite_difference_step: 1e-7,
        }
    }
}
impl MatrixFreeNewtonConfig {
    pub fn validate(&self) -> Result<()> {
        require(
            self.max_iterations > 0
                && self.linear_max_iterations > 0
                && self.linear_restart > 0
                && self.line_search_steps > 0,
            "matrix-free iteration controls",
        )?;
        for (x, n) in [
            (self.absolute_tolerance, "absolute tolerance"),
            (self.relative_tolerance, "relative tolerance"),
            (self.linear_relative_tolerance, "linear tolerance"),
            (self.minimum_step, "minimum step"),
            (self.finite_difference_step, "finite difference step"),
        ] {
            positive(x, n)?;
        }
        require(self.minimum_step <= 1.0, "minimum step exceeds one")
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct MatrixFreeNewtonResult {
    pub iterations: usize,
    pub linear_iterations: usize,
    pub initial_residual_rms: f64,
    pub residual_rms: f64,
    pub accepted_step: f64,
    pub converged: bool,
}
pub fn finite_difference_jacobian_action(
    residual: &ResidualAction<'_>,
    x: &[f64],
    base: &[f64],
    direction: &[f64],
    output: &mut [f64],
    relative_step: f64,
) -> Result<()> {
    let n = x.len();
    require(
        n == base.len() && n == direction.len() && n == output.len(),
        "Jacobian action dimensions",
    )?;
    positive(relative_step, "Jacobian step")?;
    let norm = l2_norm(direction);
    if norm <= 0.0 {
        output.fill(0.0);
        return Ok(());
    }
    let step = relative_step * l2_norm(x).max(1.0) / norm;
    let trial: Vec<_> = x.iter().zip(direction).map(|(x, d)| x + step * d).collect();
    residual(&trial, output)?;
    for i in 0..n {
        output[i] = (output[i] - base[i]) / step;
    }
    Ok(())
}
pub fn finite_difference_transpose_jacobian_action(
    residual: &ResidualAction<'_>,
    x: &[f64],
    base: &[f64],
    direction: &[f64],
    output: &mut [f64],
    relative_step: f64,
) -> Result<()> {
    let n = x.len();
    require(
        n == base.len() && n == direction.len() && n == output.len(),
        "transpose Jacobian action dimensions",
    )?;
    positive(relative_step, "transpose Jacobian step")?;
    let mut trial = x.to_vec();
    let mut result = vec![0.0; n];
    for col in 0..n {
        let h = relative_step * x[col].abs().max(1.0);
        trial[col] += h;
        residual(&trial, &mut result)?;
        output[col] = (0..n)
            .map(|i| direction[i] * (result[i] - base[i]) / h)
            .sum();
        trial[col] = x[col];
    }
    Ok(())
}
pub fn solve_matrix_free_newton(
    x: &mut [f64],
    system: &MatrixFreeNonlinearSystem<'_>,
    c: MatrixFreeNewtonConfig,
) -> Result<MatrixFreeNewtonResult> {
    c.validate()?;
    require(
        !x.is_empty() && system.size == x.len(),
        "matrix-free system size",
    )?;
    finite_slice(x, "matrix-free initial guess")?;
    let n = x.len();
    let mut r = vec![0.0; n];
    (system.residual)(x, &mut r)?;
    let mut norm = rms(&r);
    let initial = norm;
    let target = c.absolute_tolerance.max(c.relative_tolerance * initial);
    let mut result = MatrixFreeNewtonResult {
        iterations: 0,
        linear_iterations: 0,
        initial_residual_rms: initial,
        residual_rms: norm,
        accepted_step: 0.0,
        converged: false,
    };
    if !norm.is_finite() {
        return Ok(result);
    }
    if norm <= target {
        result.converged = true;
        result.accepted_step = 1.0;
        return Ok(result);
    }
    for iteration in 0..c.max_iterations {
        let rhs: Vec<_> = r.iter().map(|v| -v).collect();
        let mut delta = vec![0.0; n];
        let apply = |v: &[f64], out: &mut [f64]| {
            if let Some(j) = system.jacobian_action {
                j(x, v, out)
            } else {
                finite_difference_jacobian_action(
                    system.residual,
                    x,
                    &r,
                    v,
                    out,
                    c.finite_difference_step,
                )
            }
        };
        let linear = restarted_gmres(
            &rhs,
            &mut delta,
            apply,
            &IdentityPreconditioner,
            LinearConfig {
                max_iterations: c.linear_max_iterations,
                relative_tolerance: c.linear_relative_tolerance,
                restart: c.linear_restart.min(n),
            },
        )?;
        result.linear_iterations += linear.iterations;
        if !linear.converged {
            result.iterations = iteration;
            return Ok(result);
        }
        let mut step = 1.0;
        let mut trial = x.to_vec();
        let mut trial_r = vec![0.0; n];
        let mut accepted = false;
        for _ in 0..c.line_search_steps {
            for i in 0..n {
                trial[i] = x[i] + step * delta[i];
            }
            (system.residual)(&trial, &mut trial_r)?;
            let next = rms(&trial_r);
            if next.is_finite() && next < norm {
                norm = next;
                accepted = true;
                break;
            }
            step *= 0.5;
            if step < c.minimum_step {
                break;
            }
        }
        result.iterations = iteration + 1;
        result.accepted_step = step;
        if !accepted {
            return Ok(result);
        }
        x.copy_from_slice(&trial);
        r = trial_r;
        result.residual_rms = norm;
        if norm <= target {
            result.converged = true;
            return Ok(result);
        }
    }
    Ok(result)
}
pub fn solve_matrix_free_adjoint(
    x: &[f64],
    system: &MatrixFreeNonlinearSystem<'_>,
    rhs: &[f64],
    adjoint: &mut [f64],
    c: MatrixFreeNewtonConfig,
) -> Result<IterativeSolverResult> {
    c.validate()?;
    require(
        !x.is_empty() && system.size == x.len() && rhs.len() == x.len() && adjoint.len() == x.len(),
        "adjoint system size",
    )?;
    let mut base = vec![0.0; x.len()];
    (system.residual)(x, &mut base)?;
    finite_slice(&base, "adjoint base residual")?;
    adjoint.fill(0.0);
    restarted_gmres(
        rhs,
        adjoint,
        |v, out| {
            if let Some(j) = system.transpose_jacobian_action {
                j(x, v, out)
            } else {
                finite_difference_transpose_jacobian_action(
                    system.residual,
                    x,
                    &base,
                    v,
                    out,
                    c.finite_difference_step,
                )
            }
        },
        &IdentityPreconditioner,
        LinearConfig {
            max_iterations: c.linear_max_iterations,
            relative_tolerance: c.linear_relative_tolerance,
            restart: c.linear_restart,
        },
    )
}

pub type DaeResidual<'a> = dyn Fn(f64, &[f64], &[f64], &mut [f64]) -> Result<()> + 'a;
pub struct DaeSystem<'a> {
    pub size: usize,
    pub residual: &'a DaeResidual<'a>,
}
pub struct DaeEvent<'a> {
    pub name: String,
    pub value: &'a dyn Fn(f64, &[f64]) -> f64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DaeTrajectory {
    pub times: Vec<f64>,
    pub states: Vec<Vec<f64>>,
    pub events: Vec<String>,
    pub converged: bool,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct DaeStepResult {
    pub time: f64,
    pub nonlinear_iterations: usize,
    pub residual_rms: f64,
    pub converged: bool,
}
pub fn implicit_euler_step(
    system: &DaeSystem<'_>,
    time: f64,
    step: f64,
    previous: &[f64],
    state: &mut [f64],
    c: MatrixFreeNewtonConfig,
) -> Result<DaeStepResult> {
    require(
        system.size > 0
            && previous.len() == system.size
            && state.len() == system.size
            && time.is_finite(),
        "DAE system/time",
    )?;
    positive(step, "DAE time step")?;
    let next = time + step;
    require(next.is_finite(), "DAE time overflow")?;
    let residual = |x: &[f64], out: &mut [f64]| {
        let derivative: Vec<_> = x
            .iter()
            .zip(previous)
            .map(|(x, old)| (x - old) / step)
            .collect();
        (system.residual)(next, x, &derivative, out)
    };
    let nonlinear = MatrixFreeNonlinearSystem {
        size: system.size,
        residual: &residual,
        jacobian_action: None,
        transpose_jacobian_action: None,
    };
    let r = solve_matrix_free_newton(state, &nonlinear, c)?;
    Ok(DaeStepResult {
        time: next,
        nonlinear_iterations: r.iterations,
        residual_rms: r.residual_rms,
        converged: r.converged,
    })
}
pub fn integrate_implicit_euler(
    system: &DaeSystem<'_>,
    time: f64,
    dt: f64,
    initial: &[f64],
    steps: usize,
    events: &[DaeEvent<'_>],
    c: MatrixFreeNewtonConfig,
) -> Result<DaeTrajectory> {
    require(
        steps > 0 && initial.len() == system.size,
        "DAE trajectory size",
    )?;
    positive(dt, "DAE dt")?;
    finite_slice(initial, "DAE state")?;
    let mut result = DaeTrajectory {
        times: vec![time],
        states: vec![initial.to_vec()],
        events: vec![],
        converged: false,
    };
    let mut state = initial.to_vec();
    let mut current = time;
    let mut old_events = Vec::new();
    for e in events {
        require(!e.name.is_empty(), "DAE event needs name")?;
        let value = (e.value)(time, initial);
        require(value.is_finite(), "nonfinite DAE event")?;
        old_events.push(value);
    }
    for _ in 0..steps {
        let old = state.clone();
        let r = implicit_euler_step(system, current, dt, &old, &mut state, c)?;
        if !r.converged {
            return Ok(result);
        }
        current = r.time;
        for (i, e) in events.iter().enumerate() {
            let value = (e.value)(current, &state);
            require(value.is_finite(), "nonfinite DAE event")?;
            if (old_events[i] < 0.0 && value >= 0.0) || (old_events[i] > 0.0 && value <= 0.0) {
                result.events.push(e.name.clone());
            }
            old_events[i] = value;
        }
        result.times.push(current);
        result.states.push(state.clone());
    }
    result.converged = true;
    Ok(result)
}
