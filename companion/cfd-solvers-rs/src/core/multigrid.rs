//! The upstream red/black smoother, full-weighting restriction and cross-stencil correction V-cycle.
use crate::{finite_slice, positive, require, Result};
#[derive(Clone, Copy, Debug)]
pub struct MultigridResult {
    pub converged: bool,
    pub cycles: usize,
    pub relative_residual: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct PoissonMultigrid2D {
    n: usize,
}
impl PoissonMultigrid2D {
    pub fn new(n: usize) -> Result<Self> {
        require(
            n >= 3
                && n.checked_add(1).is_some_and(|v| v.is_power_of_two())
                && n.checked_mul(n).is_some(),
            "multigrid interior size must be 2^k-1",
        )?;
        Ok(Self { n })
    }
    pub fn size(self) -> usize {
        self.n
    }
    fn neighbors(n: usize, u: &[f64], i: usize, j: usize) -> f64 {
        let q = j * n + i;
        let mut v = 0.0;
        if i > 0 {
            v += u[q - 1];
        }
        if i + 1 < n {
            v += u[q + 1];
        }
        if j > 0 {
            v += u[q - n];
        }
        if j + 1 < n {
            v += u[q + n];
        }
        v
    }
    fn smooth(n: usize, f: &[f64], u: &mut [f64], sweeps: usize) {
        let h = 1.0 / (n + 1) as f64;
        for _ in 0..sweeps {
            for color in 0..2 {
                for j in 0..n {
                    for i in 0..n {
                        if (i + j) % 2 == color {
                            u[j * n + i] =
                                0.25 * (h * h * f[j * n + i] + Self::neighbors(n, u, i, j));
                        }
                    }
                }
            }
        }
    }
    fn residual(n: usize, f: &[f64], u: &[f64]) -> Vec<f64> {
        let invh2 = ((n + 1) as f64).powi(2);
        let mut r = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                let q = j * n + i;
                r[q] = f[q] - (4.0 * u[q] - Self::neighbors(n, u, i, j)) * invh2;
            }
        }
        r
    }
    fn vcycle(n: usize, f: &[f64], u: &mut [f64]) {
        if n <= 3 {
            Self::smooth(n, f, u, 20);
            return;
        }
        Self::smooth(n, f, u, 3);
        let r = Self::residual(n, f, u);
        let nc = (n - 1) / 2;
        let mut fc = vec![0.0; nc * nc];
        let mut ec = vec![0.0; nc * nc];
        for j in 0..nc {
            for i in 0..nc {
                let q = (2 * j + 1) * n + 2 * i + 1;
                fc[j * nc + i] = (4.0 * r[q]
                    + 2.0 * (r[q - 1] + r[q + 1] + r[q - n] + r[q + n])
                    + r[q - n - 1]
                    + r[q - n + 1]
                    + r[q + n - 1]
                    + r[q + n + 1])
                    / 16.0;
            }
        }
        Self::vcycle(nc, &fc, &mut ec);
        for j in 0..nc {
            for i in 0..nc {
                let x = 2 * i + 1;
                let y = 2 * j + 1;
                let q = y * n + x;
                let e = ec[j * nc + i];
                u[q] += e;
                if x > 0 {
                    u[q - 1] += 0.5 * e;
                }
                if x + 1 < n {
                    u[q + 1] += 0.5 * e;
                }
                if y > 0 {
                    u[q - n] += 0.5 * e;
                }
                if y + 1 < n {
                    u[q + n] += 0.5 * e;
                }
            }
        }
        Self::smooth(n, f, u, 3);
    }
    pub fn solve(
        self,
        rhs: &[f64],
        solution: &mut [f64],
        max_cycles: usize,
        tolerance: f64,
    ) -> Result<MultigridResult> {
        require(
            rhs.len() == self.n * self.n && solution.len() == rhs.len() && max_cycles > 0,
            "multigrid dimensions/cycle count",
        )?;
        positive(tolerance, "multigrid tolerance")?;
        finite_slice(rhs, "multigrid RHS")?;
        finite_slice(solution, "multigrid initial solution")?;
        let mut u = solution.to_vec();
        let norm = rhs.iter().map(|v| v * v).sum::<f64>().sqrt();
        let norm = if norm == 0.0 { 1.0 } else { norm };
        let mut relative = f64::INFINITY;
        for c in 0..max_cycles {
            Self::vcycle(self.n, rhs, &mut u);
            let r = Self::residual(self.n, rhs, &u);
            relative = r.iter().map(|v| v * v).sum::<f64>().sqrt() / norm;
            require(relative.is_finite(), "multigrid residual became nonfinite")?;
            if relative < tolerance {
                solution.copy_from_slice(&u);
                return Ok(MultigridResult {
                    converged: true,
                    cycles: c + 1,
                    relative_residual: relative,
                });
            }
        }
        solution.copy_from_slice(&u);
        Ok(MultigridResult {
            converged: false,
            cycles: max_cycles,
            relative_residual: relative,
        })
    }
}
