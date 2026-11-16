use super::{dot, l2_norm, rms};
use crate::{finite_slice, positive, require, Result, SolverError};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(try_from = "RawCsrMatrix")]
pub struct CsrMatrix {
    rows: usize,
    cols: usize,
    row_offsets: Vec<usize>,
    column_indices: Vec<usize>,
    values: Vec<f64>,
}
#[derive(Deserialize)]
struct RawCsrMatrix {
    rows: usize,
    cols: usize,
    row_offsets: Vec<usize>,
    column_indices: Vec<usize>,
    values: Vec<f64>,
}
impl TryFrom<RawCsrMatrix> for CsrMatrix {
    type Error = SolverError;
    fn try_from(r: RawCsrMatrix) -> Result<Self> {
        Self::new(r.rows, r.cols, r.row_offsets, r.column_indices, r.values)
    }
}
impl CsrMatrix {
    pub fn new(
        rows: usize,
        cols: usize,
        row_offsets: Vec<usize>,
        column_indices: Vec<usize>,
        values: Vec<f64>,
    ) -> Result<Self> {
        require(
            rows.checked_add(1) == Some(row_offsets.len()) && column_indices.len() == values.len(),
            "CSR shape mismatch",
        )?;
        require(
            row_offsets.first() == Some(&0) && row_offsets.last() == Some(&values.len()),
            "CSR offsets invalid",
        )?;
        finite_slice(&values, "CSR values")?;
        require(
            row_offsets.windows(2).all(|w| w[0] <= w[1])
                && row_offsets.iter().all(|&v| v <= values.len()),
            "CSR offsets out of range",
        )?;
        for row in 0..rows {
            require(
                row_offsets[row] <= row_offsets[row + 1],
                "CSR offsets not monotone",
            )?;
            let slice = &column_indices[row_offsets[row]..row_offsets[row + 1]];
            require(
                slice.iter().all(|&c| c < cols) && slice.windows(2).all(|w| w[0] < w[1]),
                "CSR columns must be in range and strictly increasing",
            )?;
        }
        Ok(Self {
            rows,
            cols,
            row_offsets,
            column_indices,
            values,
        })
    }
    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }
    pub fn nonzeros(&self) -> usize {
        self.values.len()
    }
    pub fn row_offsets(&self) -> &[usize] {
        &self.row_offsets
    }
    pub fn columns(&self) -> &[usize] {
        &self.column_indices
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    pub fn multiply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
        require(
            x.len() == self.cols && y.len() == self.rows,
            "CSR multiply dimensions",
        )?;
        for (i, out) in y.iter_mut().enumerate() {
            *out = (self.row_offsets[i]..self.row_offsets[i + 1])
                .map(|k| self.values[k] * x[self.column_indices[k]])
                .sum();
        }
        Ok(())
    }
    pub fn apply(&self, x: &[f64]) -> Result<Vec<f64>> {
        let mut y = vec![0.0; self.rows];
        self.multiply(x, &mut y)?;
        Ok(y)
    }
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.rows)
            .map(|r| self.position(r, r).map_or(0.0, |k| self.values[k]))
            .collect()
    }
    pub fn position(&self, row: usize, col: usize) -> Option<usize> {
        if row >= self.rows {
            return None;
        }
        let begin = self.row_offsets[row];
        self.column_indices[begin..self.row_offsets[row + 1]]
            .binary_search(&col)
            .ok()
            .map(|i| begin + i)
    }
    pub fn transpose(&self) -> Result<Self> {
        let mut b = CsrBuilder::new(self.cols, self.rows);
        for r in 0..self.rows {
            for k in self.row_offsets[r]..self.row_offsets[r + 1] {
                b.add(self.column_indices[k], r, self.values[k])?;
            }
        }
        b.build(0.0)
    }
    pub fn to_dense(&self) -> Vec<f64> {
        let mut a = vec![0.0; self.rows * self.cols];
        for r in 0..self.rows {
            for k in self.row_offsets[r]..self.row_offsets[r + 1] {
                a[r * self.cols + self.column_indices[k]] = self.values[k];
            }
        }
        a
    }
    pub fn from_dense(a: &[f64], rows: usize, cols: usize, drop: f64) -> Result<Self> {
        require(
            rows.checked_mul(cols) == Some(a.len()),
            "dense-to-CSR dimensions",
        )?;
        let mut b = CsrBuilder::new(rows, cols);
        for i in 0..rows {
            for j in 0..cols {
                b.add(i, j, a[i * cols + j])?;
            }
        }
        b.build(drop)
    }
}
#[derive(Clone, Debug)]
pub struct CsrBuilder {
    rows: usize,
    cols: usize,
    entries: Vec<BTreeMap<usize, f64>>,
}
impl CsrBuilder {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            entries: vec![BTreeMap::new(); rows],
        }
    }
    pub fn add(&mut self, row: usize, col: usize, value: f64) -> Result<()> {
        require(
            row < self.rows && col < self.cols && value.is_finite(),
            "CSR add index/value invalid",
        )?;
        *self.entries[row].entry(col).or_default() += value;
        Ok(())
    }
    pub fn set(&mut self, row: usize, col: usize, value: f64) -> Result<()> {
        require(
            row < self.rows && col < self.cols && value.is_finite(),
            "CSR set index/value invalid",
        )?;
        self.entries[row].insert(col, value);
        Ok(())
    }
    pub fn dirichlet(&mut self, rhs: &mut [f64], node: usize, value: f64) -> Result<()> {
        require(
            self.rows == self.cols
                && rhs.len() == self.rows
                && node < self.rows
                && value.is_finite(),
            "Dirichlet dimensions/value",
        )?;
        for i in 0..self.rows {
            if i != node {
                if let Some(coefficient) = self.entries[i].remove(&node) {
                    rhs[i] -= coefficient * value;
                }
            }
        }
        self.entries[node].clear();
        self.entries[node].insert(node, 1.0);
        rhs[node] = value;
        Ok(())
    }
    pub fn build(&self, drop: f64) -> Result<CsrMatrix> {
        require(drop.is_finite() && drop >= 0.0, "CSR drop tolerance")?;
        let mut offsets = Vec::with_capacity(self.rows + 1);
        let (mut columns, mut values) = (vec![], vec![]);
        offsets.push(0);
        for row in &self.entries {
            for (&col, &v) in row {
                if v.abs() > drop {
                    columns.push(col);
                    values.push(v);
                }
            }
            offsets.push(values.len());
        }
        CsrMatrix::new(self.rows, self.cols, offsets, columns, values)
    }
}

pub trait Preconditioner {
    fn apply(&self, rhs: &[f64], result: &mut [f64]) -> Result<()>;
}
#[derive(Clone, Debug, Default)]
pub struct IdentityPreconditioner;
impl Preconditioner for IdentityPreconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<()> {
        require(r.len() == z.len(), "identity preconditioner size")?;
        z.copy_from_slice(r);
        Ok(())
    }
}
#[derive(Clone, Debug)]
pub struct JacobiPreconditioner {
    inverse: Vec<f64>,
}
impl JacobiPreconditioner {
    pub fn new(diagonal: &[f64], floor: f64) -> Result<Self> {
        positive(floor, "Jacobi floor")?;
        finite_slice(diagonal, "Jacobi diagonal")?;
        Ok(Self {
            inverse: diagonal
                .iter()
                .map(|&x| if x.abs() >= floor { 1.0 / x } else { 1.0 })
                .collect(),
        })
    }
}
impl Preconditioner for JacobiPreconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<()> {
        require(
            r.len() == self.inverse.len() && z.len() == r.len(),
            "Jacobi size",
        )?;
        for i in 0..r.len() {
            z[i] = r[i] * self.inverse[i];
        }
        Ok(())
    }
}
#[derive(Clone, Debug)]
pub struct Ilu0Preconditioner {
    matrix: CsrMatrix,
    diagonal: Vec<usize>,
}
impl Ilu0Preconditioner {
    pub fn new(matrix: &CsrMatrix, floor: f64) -> Result<Self> {
        require(matrix.rows == matrix.cols, "ILU0 needs a square matrix")?;
        positive(floor, "ILU0 floor")?;
        let mut a = matrix.clone();
        let mut diag = vec![0; a.rows];
        for (i, d) in diag.iter_mut().enumerate() {
            *d = a.position(i, i).ok_or_else(|| {
                SolverError::InvalidInput("ILU0 requires an explicit diagonal".into())
            })?;
        }
        for i in 0..a.rows {
            for k in a.row_offsets[i]..diag[i] {
                let j = a.column_indices[k];
                let p = diag[j];
                if a.values[p].abs() < floor {
                    a.values[p] = if a.values[p] < 0.0 { -floor } else { floor };
                }
                let value = a.values[k] / a.values[p];
                a.values[k] = value;
                for m in p + 1..a.row_offsets[j + 1] {
                    if let Some(target) = a.position(i, a.column_indices[m]) {
                        a.values[target] -= value * a.values[m];
                    }
                }
            }
            let d = diag[i];
            if a.values[d].abs() < floor {
                a.values[d] = if a.values[d] < 0.0 { -floor } else { floor };
            }
        }
        finite_slice(&a.values, "ILU0 factors")?;
        Ok(Self {
            matrix: a,
            diagonal: diag,
        })
    }
}
impl Preconditioner for Ilu0Preconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<()> {
        let a = &self.matrix;
        require(r.len() == a.rows && z.len() == a.rows, "ILU0 size")?;
        for i in 0..a.rows {
            let mut value = r[i];
            for k in a.row_offsets[i]..self.diagonal[i] {
                value -= a.values[k] * z[a.column_indices[k]];
            }
            z[i] = value;
        }
        for i in (0..a.rows).rev() {
            let mut value = z[i];
            for k in self.diagonal[i] + 1..a.row_offsets[i + 1] {
                value -= a.values[k] * z[a.column_indices[k]];
            }
            z[i] = value / a.values[self.diagonal[i]];
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LinearConfig {
    pub max_iterations: usize,
    pub relative_tolerance: f64,
    pub restart: usize,
}
impl Default for LinearConfig {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            relative_tolerance: 1e-10,
            restart: 40,
        }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct IterativeSolverResult {
    pub iterations: usize,
    pub residual_rms: f64,
    pub converged: bool,
}
fn validate_linear(rhs: &[f64], x: &[f64], c: LinearConfig) -> Result<()> {
    require(
        !rhs.is_empty() && x.len() == rhs.len() && c.max_iterations > 0 && c.restart > 0,
        "linear solve dimensions/iterations",
    )?;
    positive(c.relative_tolerance, "linear tolerance")?;
    finite_slice(rhs, "rhs")?;
    finite_slice(x, "initial guess")
}

pub fn preconditioned_conjugate_gradient<F, P>(
    rhs: &[f64],
    x: &mut [f64],
    apply: F,
    precondition: &P,
    c: LinearConfig,
) -> Result<IterativeSolverResult>
where
    F: Fn(&[f64], &mut [f64]) -> Result<()>,
    P: Preconditioner + ?Sized,
{
    validate_linear(rhs, x, c)?;
    let n = rhs.len();
    let (mut q, mut z) = (vec![0.0; n], vec![0.0; n]);
    apply(x, &mut q)?;
    let mut r: Vec<f64> = rhs.iter().zip(&q).map(|(b, a)| b - a).collect();
    let target = c.relative_tolerance * rms(rhs).max(1.0);
    let mut residual = rms(&r);
    if residual <= target {
        return Ok(IterativeSolverResult {
            iterations: 0,
            residual_rms: residual,
            converged: true,
        });
    }
    precondition.apply(&r, &mut z)?;
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    for k in 0..c.max_iterations {
        apply(&p, &mut q)?;
        let pq = dot(&p, &q);
        if pq <= 0.0 || !pq.is_finite() || !rz.is_finite() {
            return Ok(IterativeSolverResult {
                iterations: k,
                residual_rms: residual,
                converged: false,
            });
        }
        let alpha = rz / pq;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }
        residual = rms(&r);
        if residual <= target {
            return Ok(IterativeSolverResult {
                iterations: k + 1,
                residual_rms: residual,
                converged: true,
            });
        }
        precondition.apply(&r, &mut z)?;
        let new_rz = dot(&r, &z);
        if rz == 0.0 || !new_rz.is_finite() {
            return Ok(IterativeSolverResult {
                iterations: k + 1,
                residual_rms: residual,
                converged: false,
            });
        }
        let beta = new_rz / rz;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz = new_rz;
    }
    Ok(IterativeSolverResult {
        iterations: c.max_iterations,
        residual_rms: residual,
        converged: false,
    })
}
pub fn conjugate_gradient(
    a: &CsrMatrix,
    rhs: &[f64],
    x: &mut [f64],
    c: LinearConfig,
) -> Result<IterativeSolverResult> {
    require(a.rows == a.cols && a.rows == rhs.len(), "CG matrix size")?;
    let jacobi = JacobiPreconditioner::new(&a.diagonal(), 1e-30)?;
    preconditioned_conjugate_gradient(rhs, x, |v, y| a.multiply(v, y), &jacobi, c)
}

pub fn bicgstab<F, P>(
    rhs: &[f64],
    x: &mut [f64],
    apply: F,
    precondition: &P,
    c: LinearConfig,
) -> Result<IterativeSolverResult>
where
    F: Fn(&[f64], &mut [f64]) -> Result<()>,
    P: Preconditioner + ?Sized,
{
    validate_linear(rhs, x, c)?;
    let n = rhs.len();
    let (mut v, mut p, mut s, mut t, mut z, mut q) = (
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
    );
    apply(x, &mut q)?;
    let mut r: Vec<f64> = rhs.iter().zip(&q).map(|(b, a)| b - a).collect();
    let r0 = r.clone();
    let target = c.relative_tolerance * rms(rhs).max(1.0);
    let mut residual = rms(&r);
    if residual <= target {
        return Ok(IterativeSolverResult {
            iterations: 0,
            residual_rms: residual,
            converged: true,
        });
    }
    let (mut rho_old, mut alpha, mut omega) = (1.0, 1.0, 1.0);
    let tiny = 64.0 * f64::EPSILON;
    for k in 0..c.max_iterations {
        let rho = dot(&r0, &r);
        if rho.abs() <= tiny || !rho.is_finite() {
            return Ok(IterativeSolverResult {
                iterations: k,
                residual_rms: residual,
                converged: false,
            });
        }
        let beta = (rho / rho_old) * (alpha / omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        precondition.apply(&p, &mut z)?;
        apply(&z, &mut v)?;
        let denom = dot(&r0, &v);
        if denom.abs() <= tiny || !denom.is_finite() {
            return Ok(IterativeSolverResult {
                iterations: k,
                residual_rms: residual,
                converged: false,
            });
        }
        alpha = rho / denom;
        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }
        if rms(&s) <= target {
            for i in 0..n {
                x[i] += alpha * z[i];
            }
            return Ok(IterativeSolverResult {
                iterations: k + 1,
                residual_rms: rms(&s),
                converged: true,
            });
        }
        precondition.apply(&s, &mut q)?;
        apply(&q, &mut t)?;
        let tt = dot(&t, &t);
        if tt <= tiny || !tt.is_finite() {
            return Ok(IterativeSolverResult {
                iterations: k,
                residual_rms: residual,
                converged: false,
            });
        }
        omega = dot(&t, &s) / tt;
        if omega.abs() <= tiny || !omega.is_finite() {
            return Ok(IterativeSolverResult {
                iterations: k,
                residual_rms: residual,
                converged: false,
            });
        }
        for i in 0..n {
            x[i] += alpha * z[i] + omega * q[i];
            r[i] = s[i] - omega * t[i];
        }
        residual = rms(&r);
        if residual <= target {
            return Ok(IterativeSolverResult {
                iterations: k + 1,
                residual_rms: residual,
                converged: true,
            });
        }
        rho_old = rho;
    }
    Ok(IterativeSolverResult {
        iterations: c.max_iterations,
        residual_rms: residual,
        converged: false,
    })
}

pub fn restarted_gmres<F, P>(
    rhs: &[f64],
    x: &mut [f64],
    apply: F,
    precondition: &P,
    c: LinearConfig,
) -> Result<IterativeSolverResult>
where
    F: Fn(&[f64], &mut [f64]) -> Result<()>,
    P: Preconditioner + ?Sized,
{
    validate_linear(rhs, x, c)?;
    let n = rhs.len();
    let restart = c.restart.min(c.max_iterations);
    let (mut ax, mut z, mut w) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let mut r = vec![0.0; n];
    let target = c.relative_tolerance * rms(rhs).max(1.0);
    let mut completed = 0;
    apply(x, &mut ax)?;
    for i in 0..n {
        r[i] = rhs[i] - ax[i];
    }
    let mut residual = rms(&r);
    if residual <= target {
        return Ok(IterativeSolverResult {
            iterations: 0,
            residual_rms: residual,
            converged: true,
        });
    }
    while completed < c.max_iterations {
        precondition.apply(&r, &mut z)?;
        let beta = l2_norm(&z);
        if beta <= 0.0 || !beta.is_finite() {
            return Ok(IterativeSolverResult {
                iterations: completed,
                residual_rms: residual,
                converged: false,
            });
        }
        let mut v = vec![vec![0.0; n]; restart + 1];
        for i in 0..n {
            v[0][i] = z[i] / beta;
        }
        let mut h = vec![0.0; (restart + 1) * restart];
        let (mut cs, mut sn, mut g) = (
            vec![0.0; restart],
            vec![0.0; restart],
            vec![0.0; restart + 1],
        );
        g[0] = beta;
        let limit = restart.min(c.max_iterations - completed);
        let mut used = 0;
        for j in 0..limit {
            apply(&v[j], &mut ax)?;
            precondition.apply(&ax, &mut w)?;
            for i in 0..=j {
                let hij = dot(&w, &v[i]);
                h[i * restart + j] = hij;
                for k in 0..n {
                    w[k] -= hij * v[i][k];
                }
            }
            let hn = l2_norm(&w);
            h[(j + 1) * restart + j] = hn;
            if hn > f64::EPSILON {
                for k in 0..n {
                    v[j + 1][k] = w[k] / hn;
                }
            }
            for i in 0..j {
                let a = h[i * restart + j];
                let b = h[(i + 1) * restart + j];
                h[i * restart + j] = cs[i] * a + sn[i] * b;
                h[(i + 1) * restart + j] = -sn[i] * a + cs[i] * b;
            }
            let a = h[j * restart + j];
            let b = h[(j + 1) * restart + j];
            let d = a.hypot(b);
            cs[j] = if d > 0.0 { a / d } else { 1.0 };
            sn[j] = if d > 0.0 { b / d } else { 0.0 };
            h[j * restart + j] = cs[j] * a + sn[j] * b;
            h[(j + 1) * restart + j] = 0.0;
            let gj = g[j];
            g[j] = cs[j] * gj;
            g[j + 1] = -sn[j] * gj;
            used = j + 1;
            completed += 1;
            if g[j + 1].abs() / (n as f64).sqrt() <= target
                || hn <= f64::EPSILON
                || completed == c.max_iterations
            {
                break;
            }
        }
        let mut y = vec![0.0; used];
        for i in (0..used).rev() {
            let mut sum = g[i];
            for j in i + 1..used {
                sum -= h[i * restart + j] * y[j];
            }
            let d = h[i * restart + i];
            if d.abs() <= f64::EPSILON {
                return Ok(IterativeSolverResult {
                    iterations: completed,
                    residual_rms: residual,
                    converged: false,
                });
            }
            y[i] = sum / d;
        }
        for j in 0..used {
            for i in 0..n {
                x[i] += y[j] * v[j][i];
            }
        }
        apply(x, &mut ax)?;
        for i in 0..n {
            r[i] = rhs[i] - ax[i];
        }
        residual = rms(&r);
        if residual <= target {
            return Ok(IterativeSolverResult {
                iterations: completed,
                residual_rms: residual,
                converged: true,
            });
        }
    }
    Ok(IterativeSolverResult {
        iterations: completed,
        residual_rms: residual,
        converged: false,
    })
}

/// Outcome of [`matrix_free_conjugate_gradient`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ConjugateGradientResult {
    pub iterations: usize,
    pub residual_rms: f64,
    pub converged: bool,
}

/// Unpreconditioned matrix-free CG for symmetric positive (semi-)definite operators, ported from
/// upstream `cfd/core/conjugate_gradient` (MIT, cfd_solvers contributors). Convergence is
/// `rms(r) <= tol * max(rms(rhs), 1)`. For a semidefinite operator the caller supplies a
/// right-hand side and initial guess in the compatible subspace. As upstream, the reported
/// iteration counter advances twice per CG update, so `max_iterations` bounds that counter.
pub fn matrix_free_conjugate_gradient(
    rhs: &[f64],
    x: &mut [f64],
    mut apply: impl FnMut(&[f64], &mut [f64]),
    max_iterations: usize,
    relative_tolerance: f64,
) -> Result<ConjugateGradientResult> {
    require(
        !rhs.is_empty() && x.len() == rhs.len(),
        "CG vector size mismatch",
    )?;
    require(
        max_iterations > 0 && relative_tolerance > 0.0,
        "invalid CG controls",
    )?;
    let n = rhs.len();
    let nf = n as f64;
    let mut q = vec![0.0; n];
    apply(x, &mut q);
    let mut r: Vec<f64> = rhs.iter().zip(&q).map(|(b, a)| b - a).collect();
    let mut d = r.clone();
    let rhs2: f64 = rhs.iter().map(|v| v * v).sum();
    let mut rr: f64 = r.iter().map(|v| v * v).sum();
    let target = relative_tolerance * (rhs2 / nf).sqrt().max(1.0);
    let mut rms = (rr / nf).sqrt();
    if rms <= target {
        return Ok(ConjugateGradientResult {
            iterations: 0,
            residual_rms: rms,
            converged: true,
        });
    }
    let mut iteration = 0;
    while iteration < max_iterations {
        apply(&d, &mut q);
        let dq: f64 = d.iter().zip(&q).map(|(a, b)| a * b).sum();
        if dq.is_nan() || dq <= 0.0 || !dq.is_finite() {
            break;
        }
        let alpha = rr / dq;
        for i in 0..n {
            x[i] += alpha * d[i];
            r[i] -= alpha * q[i];
        }
        let rr_new: f64 = r.iter().map(|v| v * v).sum();
        rms = (rr_new / nf).sqrt();
        iteration += 1;
        if rms <= target {
            return Ok(ConjugateGradientResult {
                iterations: iteration,
                residual_rms: rms,
                converged: true,
            });
        }
        let beta = rr_new / rr;
        for (di, ri) in d.iter_mut().zip(&r) {
            *di = ri + beta * *di;
        }
        rr = rr_new;
        iteration += 1;
    }
    Ok(ConjugateGradientResult {
        iterations: iteration,
        residual_rms: rms,
        converged: false,
    })
}
