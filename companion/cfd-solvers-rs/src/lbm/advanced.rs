//! Periodic D2Q9 regularized, cumulant and MRT collision baselines, ported from upstream
//! `cfd/solvers/lbm/advanced_d2q9` (MIT, cfd_solvers contributors).
use super::{equilibrium, CX, CY, WEIGHT};
use crate::{require, Result, SolverError};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdvancedD2Q9Config {
    pub nx: usize,
    pub ny: usize,
    pub tau: f32,
}

impl Default for AdvancedD2Q9Config {
    fn default() -> Self {
        Self {
            nx: 64,
            ny: 64,
            tau: 0.7,
        }
    }
}

/// Shared periodic lattice storage for the three collision models.
#[derive(Clone, Debug)]
struct Lattice {
    c: AdvancedD2Q9Config,
    f: Vec<[f32; 9]>,
    next: Vec<[f32; 9]>,
}

fn macroscopic(f: &[f32; 9]) -> [f32; 3] {
    let (mut rho, mut ux, mut uy) = (0.0f32, 0.0f32, 0.0f32);
    for d in 0..9 {
        rho += f[d];
        ux += f[d] * CX[d] as f32;
        uy += f[d] * CY[d] as f32;
    }
    if rho > 0.0 {
        ux /= rho;
        uy /= rho;
    }
    [rho, ux, uy]
}

impl Lattice {
    fn new(c: AdvancedD2Q9Config, what: &str) -> Result<Self> {
        require(
            c.tau > 0.5 && c.tau.is_finite(),
            format!("invalid {what} tau"),
        )?;
        let cells =
            c.nx.checked_mul(c.ny)
                .filter(|n| *n > 0)
                .ok_or_else(|| SolverError::InvalidInput("invalid lattice size".into()))?;
        let mut l = Self {
            c,
            f: vec![[0.0; 9]; cells],
            next: vec![[0.0; 9]; cells],
        };
        l.initialize_uniform(1.0, 0.0, 0.0);
        Ok(l)
    }
    fn fill(&mut self, state: impl Fn(usize, usize) -> [f32; 3]) {
        let nx = self.c.nx;
        for (n, f) in self.f.iter_mut().enumerate() {
            let s = state(n % nx, n / nx);
            for (d, v) in f.iter_mut().enumerate() {
                *v = equilibrium(d, s[0], s[1], s[2]);
            }
        }
    }
    fn initialize_uniform(&mut self, rho: f32, ux: f32, uy: f32) {
        self.fill(|_, _| [rho, ux, uy]);
    }
    fn initialize_taylor_green(&mut self, amplitude: f32) {
        let two_pi = 2.0 * std::f32::consts::PI;
        let (nx, ny) = (self.c.nx as f32, self.c.ny as f32);
        self.fill(|x, y| {
            let xf = (x as f32 + 0.5) / nx;
            let yf = (y as f32 + 0.5) / ny;
            [
                1.0,
                amplitude * (two_pi * xf).sin() * (two_pi * yf).cos(),
                -amplitude * (two_pi * xf).cos() * (two_pi * yf).sin(),
            ]
        });
    }
    /// Pulls the post-collision populations of the upwind neighbours (periodic).
    fn gather(&self, n: usize) -> [f32; 9] {
        let (nx, ny) = (self.c.nx as i64, self.c.ny as i64);
        let (x, y) = ((n as i64) % nx, (n as i64) / nx);
        let mut fin = [0.0; 9];
        for (d, v) in fin.iter_mut().enumerate() {
            let sx = (x - i64::from(CX[d]) + nx) % nx;
            let sy = (y - i64::from(CY[d]) + ny) % ny;
            *v = self.f[(sy * nx + sx) as usize][d];
        }
        fin
    }
    fn step_with(&mut self, collide: impl Fn([f32; 9]) -> [f32; 9]) {
        for n in 0..self.f.len() {
            self.next[n] = collide(self.gather(n));
        }
        std::mem::swap(&mut self.f, &mut self.next);
    }
    fn mass(&self) -> f64 {
        self.f
            .iter()
            .map(|f| f.iter().map(|v| f64::from(*v)).sum::<f64>())
            .sum()
    }
    fn kinetic_energy(&self) -> f64 {
        self.f
            .iter()
            .map(|f| {
                let s = macroscopic(f);
                let (ux, uy) = (f64::from(s[1]), f64::from(s[2]));
                0.5 * f64::from(s[0]) * (ux * ux + uy * uy)
            })
            .sum()
    }
}

const MOMENT_MATRIX: [[f64; 9]; 9] = [
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [-4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0],
    [4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0],
    [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],
    [0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0],
    [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0],
    [0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0],
    [0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0],
];

/// Gauss-Jordan inverse with partial pivoting, mirroring upstream's elimination order.
fn invert9(m: &[[f64; 9]; 9], what: &str) -> Result<[[f64; 9]; 9]> {
    let mut a = [[0.0f64; 18]; 9];
    for i in 0..9 {
        a[i][..9].copy_from_slice(&m[i]);
        a[i][9 + i] = 1.0;
    }
    for col in 0..9 {
        let mut pivot = col;
        for row in col + 1..9 {
            if a[row][col].abs() > a[pivot][col].abs() {
                pivot = row;
            }
        }
        if a[pivot][col].abs() < 1e-14 {
            return Err(SolverError::Singular(what.into()));
        }
        a.swap(pivot, col);
        let diagonal = a[col][col];
        for v in &mut a[col] {
            *v /= diagonal;
        }
        let pivot_row = a[col];
        for (row, r) in a.iter_mut().enumerate() {
            if row == col {
                continue;
            }
            let factor = r[col];
            for (v, p) in r.iter_mut().zip(&pivot_row) {
                *v -= factor * p;
            }
        }
    }
    let mut inv = [[0.0; 9]; 9];
    for i in 0..9 {
        inv[i].copy_from_slice(&a[i][9..]);
    }
    Ok(inv)
}

fn raw_moment_matrix() -> [[f64; 9]; 9] {
    const POWERS: [[i32; 2]; 9] = [
        [0, 0],
        [1, 0],
        [0, 1],
        [2, 0],
        [0, 2],
        [1, 1],
        [2, 1],
        [1, 2],
        [2, 2],
    ];
    let ipow = |v: f64, p: i32| match p {
        0 => 1.0,
        1 => v,
        _ => v * v,
    };
    let mut a = [[0.0; 9]; 9];
    for (row, pw) in POWERS.iter().enumerate() {
        for d in 0..9 {
            a[row][d] = ipow(f64::from(CX[d]), pw[0]) * ipow(f64::from(CY[d]), pw[1]);
        }
    }
    a
}

macro_rules! common_api {
    () => {
        pub fn initialize_uniform(&mut self, rho: f32, ux: f32, uy: f32) {
            self.lattice.initialize_uniform(rho, ux, uy);
        }
        pub fn initialize_taylor_green(&mut self, amplitude: f32) {
            self.lattice.initialize_taylor_green(amplitude);
        }
        pub fn step(&mut self, count: usize) {
            for _ in 0..count {
                self.step_once();
            }
        }
        pub fn mass(&self) -> f64 {
            self.lattice.mass()
        }
        pub fn kinetic_energy(&self) -> f64 {
            self.lattice.kinetic_energy()
        }
        pub fn config(&self) -> AdvancedD2Q9Config {
            self.lattice.c
        }
        /// Populations per cell, `[cell][direction]` in the standard D2Q9 ordering.
        pub fn populations(&self) -> &[[f32; 9]] {
            &self.lattice.f
        }
    };
}

/// Regularized BGK: non-equilibrium reconstructed from its second-order moment.
#[derive(Clone, Debug)]
pub struct RegularizedD2Q9Solver {
    lattice: Lattice,
}

impl RegularizedD2Q9Solver {
    pub fn new(config: AdvancedD2Q9Config) -> Result<Self> {
        Ok(Self {
            lattice: Lattice::new(config, "regularized")?,
        })
    }
    common_api!();
    fn step_once(&mut self) {
        let omega = 1.0f32 / self.lattice.c.tau;
        self.lattice.step_with(|fin| {
            let s = macroscopic(&fin);
            let mut feq = [0.0f32; 9];
            let (mut pxx, mut pxy, mut pyy) = (0.0f64, 0.0f64, 0.0f64);
            for d in 0..9 {
                feq[d] = equilibrium(d, s[0], s[1], s[2]);
                let neq = f64::from(fin[d] - feq[d]);
                let (cx, cy) = (f64::from(CX[d]), f64::from(CY[d]));
                pxx += neq * cx * cx;
                pxy += neq * cx * cy;
                pyy += neq * cy * cy;
            }
            let mut out = [0.0f32; 9];
            for d in 0..9 {
                let (cx, cy) = (f64::from(CX[d]), f64::from(CY[d]));
                let regularized = 4.5
                    * f64::from(WEIGHT[d])
                    * ((cx * cx - 1.0 / 3.0) * pxx
                        + 2.0 * cx * cy * pxy
                        + (cy * cy - 1.0 / 3.0) * pyy);
                out[d] = (f64::from(feq[d]) + (1.0 - f64::from(omega)) * regularized) as f32;
            }
            out
        });
    }
}

/// Central-moment/cumulant collision: second-order cumulants relax at `1/tau`, the
/// non-hydrodynamic higher cumulants relax fully to the factorized equilibrium.
#[derive(Clone, Debug)]
pub struct CumulantD2Q9Solver {
    lattice: Lattice,
    raw_inverse: [[f64; 9]; 9],
}

impl CumulantD2Q9Solver {
    pub fn new(config: AdvancedD2Q9Config) -> Result<Self> {
        let raw_inverse = invert9(&raw_moment_matrix(), "singular D2Q9 raw moment transform")?;
        Ok(Self {
            lattice: Lattice::new(config, "cumulant")?,
            raw_inverse,
        })
    }
    common_api!();
    fn step_once(&mut self) {
        let omega = 1.0 / f64::from(self.lattice.c.tau);
        let cs2 = 1.0 / 3.0;
        let inv = self.raw_inverse;
        self.lattice.step_with(|fin| {
            let s = macroscopic(&fin);
            let (rho, ux, uy) = (f64::from(s[0]), f64::from(s[1]), f64::from(s[2]));
            if rho.is_nan() || rho <= 0.0 {
                return fin;
            }
            let (mut mu20, mut mu02, mut mu11, mut mu21, mut mu12, mut mu22) =
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            for d in 0..9 {
                let dx = f64::from(CX[d]) - ux;
                let dy = f64::from(CY[d]) - uy;
                let v = f64::from(fin[d]);
                let (dx2, dy2) = (dx * dx, dy * dy);
                mu20 += v * dx2;
                mu02 += v * dy2;
                mu11 += v * dx * dy;
                mu21 += v * dx2 * dy;
                mu12 += v * dx * dy2;
                mu22 += v * dx2 * dy2;
            }
            let (mut k20, mut k02, mut k11) = (mu20 / rho, mu02 / rho, mu11 / rho);
            // k21, k12 and k22 = mu22/rho - k20 k02 - 2 k11^2 are reset to zero below.
            let _ = (mu21, mu12, mu22);
            k20 -= omega * (k20 - cs2);
            k02 -= omega * (k02 - cs2);
            k11 -= omega * k11;
            let (k21, k12, k22) = (0.0, 0.0, 0.0);
            mu20 = rho * k20;
            mu02 = rho * k02;
            mu11 = rho * k11;
            mu21 = rho * k21;
            mu12 = rho * k12;
            mu22 = rho * (k22 + k20 * k02 + 2.0 * k11 * k11);
            let raw = [
                rho,
                rho * ux,
                rho * uy,
                mu20 + rho * ux * ux,
                mu02 + rho * uy * uy,
                mu11 + rho * ux * uy,
                mu21 + 2.0 * ux * mu11 + uy * mu20 + rho * ux * ux * uy,
                mu12 + 2.0 * uy * mu11 + ux * mu02 + rho * ux * uy * uy,
                mu22 + 2.0 * uy * mu21
                    + uy * uy * mu20
                    + 2.0 * ux * mu12
                    + 4.0 * ux * uy * mu11
                    + ux * ux * mu02
                    + rho * ux * ux * uy * uy,
            ];
            let mut out = [0.0f32; 9];
            for d in 0..9 {
                let mut value = 0.0;
                for row in 0..9 {
                    value += inv[d][row] * raw[row];
                }
                out[d] = value as f32;
            }
            out
        });
    }
}

/// Lallemand-Luo MRT with fixed non-hydrodynamic rates (1.1, 1.0, 1.2, 1.2) and `1/tau` shear.
#[derive(Clone, Debug)]
pub struct MrtD2Q9Solver {
    lattice: Lattice,
    inverse: [[f64; 9]; 9],
}

impl MrtD2Q9Solver {
    pub fn new(config: AdvancedD2Q9Config) -> Result<Self> {
        let inverse = invert9(&MOMENT_MATRIX, "singular MRT transform")?;
        Ok(Self {
            lattice: Lattice::new(config, "MRT")?,
            inverse,
        })
    }
    common_api!();
    fn step_once(&mut self) {
        let omega = 1.0 / f64::from(self.lattice.c.tau);
        let relaxation = [0.0, 1.1, 1.0, 0.0, 1.2, 0.0, 1.2, omega, omega];
        let inv = self.inverse;
        self.lattice.step_with(|fin| {
            let s = macroscopic(&fin);
            let feq: [f32; 9] = std::array::from_fn(|d| equilibrium(d, s[0], s[1], s[2]));
            let mut moments = [0.0f64; 9];
            for row in 0..9 {
                let mut eq = 0.0;
                for d in 0..9 {
                    moments[row] += MOMENT_MATRIX[row][d] * f64::from(fin[d]);
                    eq += MOMENT_MATRIX[row][d] * f64::from(feq[d]);
                }
                moments[row] -= relaxation[row] * (moments[row] - eq);
            }
            let mut out = [0.0f32; 9];
            for d in 0..9 {
                let mut value = 0.0;
                for row in 0..9 {
                    value += inv[d][row] * moments[row];
                }
                out[d] = value as f32;
            }
            out
        });
    }
}
