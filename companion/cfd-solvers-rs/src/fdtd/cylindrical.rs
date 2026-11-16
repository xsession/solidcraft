//! Axisymmetric TM(r, z) FDTD and cylindrical multigrid level construction, ported from
//! upstream `cfd/solvers/fdtd/{cylindrical_tm, cylindrical_multigrid}` (MIT, cfd_solvers
//! contributors).
use super::{EPSILON0, MU0};
use crate::{require, Result, SolverError};
use std::f64::consts::PI;

const C0: f64 = 299792458.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylindricalTMConfig {
    pub nr: usize,
    pub nz: usize,
    pub dr: f64,
    pub dz: f64,
    pub courant: f64,
    pub epsilon_r: f64,
    pub mu_r: f64,
}

impl Default for CylindricalTMConfig {
    fn default() -> Self {
        Self {
            nr: 64,
            nz: 128,
            dr: 1e-3,
            dz: 1e-3,
            courant: 0.45,
            epsilon_r: 1.0,
            mu_r: 1.0,
        }
    }
}

/// m = 0 TM mode (Er, Ez, Hphi) with PEC outer-radius and axial walls; the axis Ez update uses
/// the analytic `r -> 0` limit `4 Hphi / dr`.
#[derive(Clone, Debug)]
pub struct CylindricalTM {
    c: CylindricalTMConfig,
    dt: f64,
    eps: f64,
    mu: f64,
    er: Vec<f64>,
    ez: Vec<f64>,
    hphi: Vec<f64>,
}

impl CylindricalTM {
    pub fn new(c: CylindricalTMConfig) -> Result<Self> {
        require(
            c.nr >= 3
                && c.nz >= 3
                && c.dr > 0.0
                && c.dz > 0.0
                && c.courant > 0.0
                && c.courant <= 0.5
                && c.epsilon_r > 0.0
                && c.mu_r > 0.0,
            "invalid cylindrical FDTD config",
        )?;
        require(
            [c.dr, c.dz, c.epsilon_r, c.mu_r]
                .iter()
                .all(|v| v.is_finite()),
            "invalid cylindrical FDTD config",
        )?;
        let dt = c.courant
            / (C0 / (c.epsilon_r * c.mu_r).sqrt()
                * (1.0 / (c.dr * c.dr) + 1.0 / (c.dz * c.dz)).sqrt());
        Ok(Self {
            c,
            dt,
            eps: EPSILON0 * c.epsilon_r,
            mu: MU0 * c.mu_r,
            er: vec![0.0; (c.nr - 1) * c.nz],
            ez: vec![0.0; c.nr * (c.nz - 1)],
            hphi: vec![0.0; (c.nr - 1) * (c.nz - 1)],
        })
    }
    fn er_idx(&self, i: usize, j: usize) -> usize {
        j * (self.c.nr - 1) + i
    }
    fn ez_idx(&self, i: usize, j: usize) -> usize {
        j * self.c.nr + i
    }
    fn h_idx(&self, i: usize, j: usize) -> usize {
        j * (self.c.nr - 1) + i
    }
    /// Gaussian Ez at `(r0, z0)`; clears Er and Hphi (upstream default amplitude 1).
    pub fn initialize_gaussian_ez(
        &mut self,
        r0: f64,
        z0: f64,
        width: f64,
        amplitude: f64,
    ) -> Result<()> {
        require(width > 0.0, "Gaussian width must be positive")?;
        for j in 0..self.c.nz - 1 {
            for i in 0..self.c.nr {
                let (r, z) = (i as f64 * self.c.dr, (j as f64 + 0.5) * self.c.dz);
                let d2 = (r - r0) * (r - r0) + (z - z0) * (z - z0);
                let n = self.ez_idx(i, j);
                self.ez[n] = amplitude * (-d2 / (2.0 * width * width)).exp();
            }
        }
        self.er.fill(0.0);
        self.hphi.fill(0.0);
        Ok(())
    }
    pub fn add_soft_ez_source(&mut self, ir: usize, iz: usize, value: f64) -> Result<()> {
        if ir >= self.c.nr || iz >= self.c.nz - 1 {
            return Err(SolverError::InvalidInput("cylindrical source index".into()));
        }
        let n = self.ez_idx(ir, iz);
        self.ez[n] += value;
        Ok(())
    }
    pub fn step(&mut self, count: usize) {
        for _ in 0..count {
            self.step_once();
        }
    }
    fn step_once(&mut self) {
        let (nr, nz, dr, dz, dt) = (self.c.nr, self.c.nz, self.c.dr, self.c.dz, self.dt);
        for n in 0..(nr - 1) * (nz - 1) {
            let (i, j) = (n % (nr - 1), n / (nr - 1));
            let dez = (self.ez[self.ez_idx(i + 1, j)] - self.ez[self.ez_idx(i, j)]) / dr;
            let der = (self.er[self.er_idx(i, j + 1)] - self.er[self.er_idx(i, j)]) / dz;
            self.hphi[n] += dt / self.mu * (dez - der);
        }
        for n in 0..(nr - 1) * (nz - 2) {
            let (i, j) = (n % (nr - 1), 1 + n / (nr - 1));
            let e = self.er_idx(i, j);
            self.er[e] -= dt / self.eps
                * (self.hphi[self.h_idx(i, j)] - self.hphi[self.h_idx(i, j - 1)])
                / dz;
        }
        for j in 0..nz - 1 {
            let n0 = self.ez_idx(0, j);
            self.ez[n0] += dt / self.eps * 4.0 * self.hphi[self.h_idx(0, j)] / dr;
            for i in 1..nr - 1 {
                let (rp, rm, r) = ((i as f64 + 0.5) * dr, (i as f64 - 0.5) * dr, i as f64 * dr);
                let n = self.ez_idx(i, j);
                self.ez[n] += dt / self.eps
                    * (rp * self.hphi[self.h_idx(i, j)] - rm * self.hphi[self.h_idx(i - 1, j)])
                    / (r * dr);
            }
            let last = self.ez_idx(nr - 1, j);
            self.ez[last] = 0.0;
        }
        for i in 0..nr - 1 {
            let (a, b) = (self.er_idx(i, 0), self.er_idx(i, nz - 1));
            self.er[a] = 0.0;
            self.er[b] = 0.0;
        }
    }
    /// Volume-weighted field energy with the axis cell radius floored at `dr / 2`.
    pub fn energy(&self) -> f64 {
        let (nr, nz, dr, dz) = (self.c.nr, self.c.nz, self.c.dr, self.c.dz);
        let mut e = 0.0;
        for j in 0..nz - 1 {
            for i in 0..nr {
                let (r, v) = ((i as f64 * dr).max(0.5 * dr), self.ez[self.ez_idx(i, j)]);
                e += 0.5 * self.eps * v * v * 2.0 * PI * r * dr * dz;
            }
        }
        for j in 0..nz {
            for i in 0..nr - 1 {
                let (r, v) = ((i as f64 + 0.5) * dr, self.er[self.er_idx(i, j)]);
                e += 0.5 * self.eps * v * v * 2.0 * PI * r * dr * dz;
            }
        }
        for j in 0..nz - 1 {
            for i in 0..nr - 1 {
                let (r, v) = ((i as f64 + 0.5) * dr, self.hphi[self.h_idx(i, j)]);
                e += 0.5 * self.mu * v * v * 2.0 * PI * r * dr * dz;
            }
        }
        e
    }
    /// Largest `|Ez|` or `|Er|`.
    pub fn max_field(&self) -> f64 {
        self.ez
            .iter()
            .chain(&self.er)
            .fold(0.0_f64, |m, v| m.max(v.abs()))
    }
    pub fn dt(&self) -> f64 {
        self.dt
    }
    pub fn ez(&self) -> &[f64] {
        &self.ez
    }
    pub fn er(&self) -> &[f64] {
        &self.er
    }
    pub fn hphi(&self) -> &[f64] {
        &self.hphi
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CylindricalGridLevel {
    pub radial_cells: usize,
    pub axial_cells: usize,
    pub radial_spacing: f64,
    pub axial_spacing: f64,
}

/// Successive halvings (never below two cells) until `levels` or the coarsest grid is reached.
pub fn make_cylindrical_multigrid_levels(
    radial_cells: usize,
    axial_cells: usize,
    radius_m: f64,
    height_m: f64,
    levels: usize,
) -> Result<Vec<CylindricalGridLevel>> {
    require(
        radial_cells >= 2 && axial_cells >= 2 && radius_m > 0.0 && height_m > 0.0 && levels > 0,
        "invalid cylindrical multigrid controls",
    )?;
    let (mut radial, mut axial) = (radial_cells, axial_cells);
    let mut out = Vec::with_capacity(levels.min(64));
    for _ in 0..levels {
        out.push(CylindricalGridLevel {
            radial_cells: radial,
            axial_cells: axial,
            radial_spacing: radius_m / radial as f64,
            axial_spacing: height_m / axial as f64,
        });
        if radial <= 2 || axial <= 2 {
            break;
        }
        radial = (radial / 2).max(2);
        axial = (axial / 2).max(2);
    }
    Ok(out)
}
