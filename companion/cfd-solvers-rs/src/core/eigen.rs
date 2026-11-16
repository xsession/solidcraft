use super::{
    dot, preconditioned_conjugate_gradient, rms, CsrMatrix, JacobiPreconditioner, LinearConfig,
};
use crate::{positive, require, Result};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct GeneralizedEigenConfig {
    pub mode_count: usize,
    pub max_outer_iterations: usize,
    pub relative_tolerance: f64,
    pub linear_max_iterations: usize,
    pub linear_relative_tolerance: f64,
}
impl Default for GeneralizedEigenConfig {
    fn default() -> Self {
        Self {
            mode_count: 1,
            max_outer_iterations: 200,
            relative_tolerance: 1e-9,
            linear_max_iterations: 4000,
            linear_relative_tolerance: 1e-12,
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeneralizedEigenpair {
    pub eigenvalue: f64,
    pub vector: Vec<f64>,
    pub iterations: usize,
    pub relative_residual: f64,
    pub converged: bool,
}
pub fn csr_bilinear(a: &CsrMatrix, x: &[f64], y: &[f64]) -> Result<f64> {
    require(
        x.len() == y.len() && a.rows() == a.cols() && a.rows() == x.len(),
        "bilinear dimensions",
    )?;
    Ok(dot(x, &a.apply(y)?))
}
pub fn lowest_generalized_eigenpairs(
    k: &CsrMatrix,
    m: &CsrMatrix,
    c: GeneralizedEigenConfig,
) -> Result<Vec<GeneralizedEigenpair>> {
    let n = k.rows();
    require(
        n > 0
            && k.cols() == n
            && m.rows() == n
            && m.cols() == n
            && c.mode_count > 0
            && c.mode_count <= n
            && c.max_outer_iterations > 0
            && c.linear_max_iterations > 0,
        "eigen dimensions/iterations",
    )?;
    positive(c.relative_tolerance, "eigen tolerance")?;
    positive(c.linear_relative_tolerance, "eigen linear tolerance")?;
    let jacobi = JacobiPreconditioner::new(&k.diagonal(), 1e-30)?;
    let mut modes: Vec<GeneralizedEigenpair> = vec![];
    for mode in 0..c.mode_count {
        let orthogonalize = |v: &mut [f64]| -> Result<()> {
            for previous in &modes {
                let alpha = csr_bilinear(m, &previous.vector, v)?;
                for i in 0..n {
                    v[i] -= alpha * previous.vector[i];
                }
            }
            let norm2 = csr_bilinear(m, v, v)?;
            positive(norm2, "eigen M norm")?;
            for x in v {
                *x /= norm2.sqrt();
            }
            Ok(())
        };
        let mut x: Vec<_> = (0..n)
            .map(|i| {
                let p = ((i + 1) * (mode + 1)) as f64;
                (0.731 * p).sin() + 0.25 * (0.371 * p).cos()
            })
            .collect();
        orthogonalize(&mut x)?;
        let mut pair = GeneralizedEigenpair {
            eigenvalue: 0.0,
            vector: x.clone(),
            iterations: 0,
            relative_residual: f64::INFINITY,
            converged: false,
        };
        for iteration in 0..c.max_outer_iterations {
            let rhs = m.apply(&x)?;
            let mut y = vec![0.0; n];
            let linear = preconditioned_conjugate_gradient(
                &rhs,
                &mut y,
                |v, out| k.multiply(v, out),
                &jacobi,
                LinearConfig {
                    max_iterations: c.linear_max_iterations,
                    relative_tolerance: c.linear_relative_tolerance,
                    restart: 40,
                },
            )?;
            if !linear.converged {
                break;
            }
            orthogonalize(&mut y)?;
            x = y;
            let kx = k.apply(&x)?;
            let mx = m.apply(&x)?;
            let lambda = dot(&x, &kx) / dot(&x, &mx);
            let residual: Vec<_> = kx.iter().zip(&mx).map(|(a, b)| a - lambda * b).collect();
            let relative = rms(&residual) / (lambda.abs() * rms(&mx)).max(f64::MIN_POSITIVE);
            pair = GeneralizedEigenpair {
                eigenvalue: lambda,
                vector: x.clone(),
                iterations: iteration + 1,
                relative_residual: relative,
                converged: relative <= c.relative_tolerance,
            };
            if pair.converged {
                break;
            }
        }
        if !pair.converged {
            return Ok(modes);
        }
        modes.push(pair);
    }
    Ok(modes)
}
