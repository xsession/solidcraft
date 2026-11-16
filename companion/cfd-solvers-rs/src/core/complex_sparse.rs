use super::{
    restarted_gmres, CsrBuilder, CsrMatrix, IdentityPreconditioner, Ilu0Preconditioner,
    IterativeSolverResult, LinearConfig,
};
use crate::{positive, require, Complex64 as C, Result};
#[derive(Clone, Debug)]
pub struct ComplexCsrBuilder {
    dimension: usize,
    entries: Vec<Vec<(usize, C)>>,
}
impl ComplexCsrBuilder {
    pub fn new(dimension: usize) -> Result<Self> {
        require(
            dimension > 0 && dimension.checked_mul(2).is_some(),
            "invalid complex dimension",
        )?;
        Ok(Self {
            dimension,
            entries: vec![vec![]; dimension],
        })
    }
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    pub fn add(&mut self, row: usize, column: usize, value: C) -> Result<()> {
        require(
            row < self.dimension
                && column < self.dimension
                && value.re.is_finite()
                && value.im.is_finite(),
            "invalid complex sparse entry",
        )?;
        self.entries[row].push((column, value));
        Ok(())
    }
    pub fn build_real_block(&self, drop_tolerance: f64) -> Result<CsrMatrix> {
        require(
            drop_tolerance.is_finite() && drop_tolerance >= 0.0,
            "complex drop tolerance",
        )?;
        let n = self.dimension;
        let mut b = CsrBuilder::new(2 * n, 2 * n);
        for row in 0..n {
            for &(col, v) in &self.entries[row] {
                if v.re.abs() > drop_tolerance || row == col {
                    let re = if row == col && v.re.abs() <= drop_tolerance {
                        f64::MIN_POSITIVE
                    } else {
                        v.re
                    };
                    b.add(row, col, re)?;
                    b.add(n + row, n + col, re)?;
                }
                if v.im.abs() > drop_tolerance {
                    b.add(row, n + col, -v.im)?;
                    b.add(n + row, col, v.im)?;
                }
            }
        }
        b.build(0.0)
    }
    pub fn apply(&self, x: &[C]) -> Result<Vec<C>> {
        require(x.len() == self.dimension, "complex vector size")?;
        let mut y = vec![C::new(0.0, 0.0); self.dimension];
        for (r, row) in self.entries.iter().enumerate() {
            for &(c, v) in row {
                y[r] += v * x[c];
            }
        }
        Ok(y)
    }
}
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ComplexSparseSolveConfig {
    pub max_iterations: usize,
    pub gmres_restart: usize,
    pub relative_tolerance: f64,
    pub ilu_diagonal_floor: f64,
}
impl Default for ComplexSparseSolveConfig {
    fn default() -> Self {
        Self {
            max_iterations: 4000,
            gmres_restart: 60,
            relative_tolerance: 1e-10,
            ilu_diagonal_floor: 1e-24,
        }
    }
}
#[derive(Clone, Debug)]
pub struct ComplexSparseSolveResult {
    pub solution: Vec<C>,
    pub linear_result: IterativeSolverResult,
}
pub fn solve_complex_sparse(
    matrix: &ComplexCsrBuilder,
    rhs: &[C],
    c: ComplexSparseSolveConfig,
) -> Result<ComplexSparseSolveResult> {
    let n = matrix.dimension;
    require(
        rhs.len() == n && c.max_iterations > 0 && c.gmres_restart > 0,
        "complex solve dimensions/controls",
    )?;
    positive(c.relative_tolerance, "complex relative tolerance")?;
    positive(c.ilu_diagonal_floor, "complex ILU floor")?;
    require(
        rhs.iter().all(|v| v.re.is_finite() && v.im.is_finite()),
        "nonfinite complex RHS",
    )?;
    let a = matrix.build_real_block(0.0)?;
    let ilu = Ilu0Preconditioner::new(&a, c.ilu_diagonal_floor)?;
    let mut b = vec![0.0; 2 * n];
    for i in 0..n {
        b[i] = rhs[i].re;
        b[n + i] = rhs[i].im;
    }
    let mut x = vec![0.0; 2 * n];
    let mut result = restarted_gmres(
        &b,
        &mut x,
        |v, y| a.multiply(v, y),
        &ilu,
        LinearConfig {
            max_iterations: c.max_iterations,
            restart: c.gmres_restart.min(2 * n),
            relative_tolerance: c.relative_tolerance,
        },
    )?;
    if !result.converged {
        x.fill(0.0);
        result = restarted_gmres(
            &b,
            &mut x,
            |v, y| a.multiply(v, y),
            &IdentityPreconditioner,
            LinearConfig {
                max_iterations: c.max_iterations,
                restart: c.gmres_restart.max(2 * n).min(c.max_iterations),
                relative_tolerance: c.relative_tolerance,
            },
        )?;
    }
    Ok(ComplexSparseSolveResult {
        solution: (0..n).map(|i| C::new(x[i], x[n + i])).collect(),
        linear_result: result,
    })
}
