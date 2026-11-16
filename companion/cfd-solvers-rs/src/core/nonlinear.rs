use super::{restarted_gmres, rms, CsrMatrix, Ilu0Preconditioner, LinearConfig};
use crate::{finite_slice, positive, require, Result};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct NewtonConfig {
    pub max_iterations: usize,
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub linear_max_iterations: usize,
    pub gmres_restart: usize,
    pub linear_relative_tolerance: f64,
    pub line_search_steps: usize,
    pub minimum_step: f64,
}
impl Default for NewtonConfig {
    fn default() -> Self {
        Self {
            max_iterations: 30,
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-9,
            linear_max_iterations: 1000,
            gmres_restart: 40,
            linear_relative_tolerance: 1e-11,
            line_search_steps: 12,
            minimum_step: 1e-6,
        }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct NewtonResult {
    pub iterations: usize,
    pub residual_rms: f64,
    pub accepted_step: f64,
    pub converged: bool,
}
pub fn newton_solve<R, J>(
    x: &mut [f64],
    residual: R,
    jacobian: J,
    c: NewtonConfig,
) -> Result<NewtonResult>
where
    R: Fn(&[f64], &mut [f64]) -> Result<()>,
    J: Fn(&[f64]) -> Result<CsrMatrix>,
{
    require(
        !x.is_empty()
            && c.max_iterations > 0
            && c.linear_max_iterations > 0
            && c.gmres_restart > 0
            && c.line_search_steps > 0,
        "Newton controls/system incomplete",
    )?;
    positive(c.absolute_tolerance, "Newton absolute tolerance")?;
    positive(c.relative_tolerance, "Newton relative tolerance")?;
    positive(c.linear_relative_tolerance, "Newton linear tolerance")?;
    positive(c.minimum_step, "Newton minimum step")?;
    require(c.minimum_step <= 1.0, "Newton minimum step exceeds one")?;
    finite_slice(x, "Newton initial guess")?;
    let n = x.len();
    let (mut r, mut trial_r, mut trial_x) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    residual(x, &mut r)?;
    let mut error = rms(&r);
    let target = c
        .absolute_tolerance
        .max(c.relative_tolerance * error.max(f64::MIN_POSITIVE));
    if !error.is_finite() {
        return Ok(NewtonResult {
            iterations: 0,
            residual_rms: error,
            accepted_step: 0.0,
            converged: false,
        });
    }
    if error <= target {
        return Ok(NewtonResult {
            iterations: 0,
            residual_rms: error,
            accepted_step: 1.0,
            converged: true,
        });
    }
    let mut accepted_step = 1.0;
    for iteration in 0..c.max_iterations {
        let j = jacobian(x)?;
        require(j.rows() == n && j.cols() == n, "Newton Jacobian size")?;
        let ilu = Ilu0Preconditioner::new(&j, 1e-12)?;
        let rhs: Vec<f64> = r.iter().map(|v| -v).collect();
        let mut delta = vec![0.0; n];
        let linear = restarted_gmres(
            &rhs,
            &mut delta,
            |v, y| j.multiply(v, y),
            &ilu,
            LinearConfig {
                max_iterations: c.linear_max_iterations,
                relative_tolerance: c.linear_relative_tolerance,
                restart: c.gmres_restart.min(n),
            },
        )?;
        if !linear.converged {
            return Ok(NewtonResult {
                iterations: iteration,
                residual_rms: error,
                accepted_step: 0.0,
                converged: false,
            });
        }
        let mut step = 1.0;
        let mut accepted = false;
        let mut trial_error = error;
        for _ in 0..c.line_search_steps {
            for i in 0..n {
                trial_x[i] = x[i] + step * delta[i];
            }
            residual(&trial_x, &mut trial_r)?;
            trial_error = rms(&trial_r);
            if trial_error.is_finite() && trial_error < error {
                accepted = true;
                break;
            }
            step *= 0.5;
            if step < c.minimum_step {
                break;
            }
        }
        if !accepted {
            return Ok(NewtonResult {
                iterations: iteration + 1,
                residual_rms: error,
                accepted_step: step,
                converged: false,
            });
        }
        x.copy_from_slice(&trial_x);
        std::mem::swap(&mut r, &mut trial_r);
        error = trial_error;
        accepted_step = step;
        if error <= target {
            return Ok(NewtonResult {
                iterations: iteration + 1,
                residual_rms: error,
                accepted_step,
                converged: true,
            });
        }
    }
    Ok(NewtonResult {
        iterations: c.max_iterations,
        residual_rms: error,
        accepted_step,
        converged: false,
    })
}
pub fn finite_difference_jacobian<R>(
    x: &[f64],
    residual: &R,
    relative_step: f64,
) -> Result<CsrMatrix>
where
    R: Fn(&[f64], &mut [f64]) -> Result<()>,
{
    positive(relative_step, "Jacobian step")?;
    let n = x.len();
    let mut base = vec![0.0; n];
    residual(x, &mut base)?;
    let mut a = vec![0.0; n * n];
    let mut trial = x.to_vec();
    let mut output = vec![0.0; n];
    for j in 0..n {
        let h = relative_step * x[j].abs().max(1.0);
        trial[j] += h;
        residual(&trial, &mut output)?;
        for i in 0..n {
            a[i * n + j] = (output[i] - base[i]) / h;
        }
        trial[j] = x[j];
    }
    CsrMatrix::from_dense(&a, n, n, 0.0)
}

pub fn finite_difference_newton<R>(
    x: &mut [f64],
    residual: R,
    c: NewtonConfig,
) -> Result<NewtonResult>
where
    R: Fn(&[f64], &mut [f64]) -> Result<()>,
{
    // Explicit finite-difference Jacobian variant. JFNK is tracked separately;
    // callers never receive a false claim that a matrix-free solve occurred.
    newton_solve(
        x,
        |v, r| residual(v, r),
        |v| finite_difference_jacobian(v, &residual, 1e-7),
        c,
    )
}
