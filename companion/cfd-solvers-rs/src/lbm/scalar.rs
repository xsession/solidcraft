//! Periodic D2Q5 BGK advection-diffusion of a passive scalar (or temperature without buoyancy),
//! ported from upstream `cfd/solvers/lbm/passive_scalar` (MIT, cfd_solvers contributors).
// Index loops mirror the upstream lattice stencils direction by direction.
#![allow(clippy::needless_range_loop)]
use crate::{require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PassiveScalarConfig {
    pub nx: usize,
    pub ny: usize,
    pub diffusivity: f64,
    pub ux: f64,
    pub uy: f64,
}

impl Default for PassiveScalarConfig {
    fn default() -> Self {
        Self {
            nx: 64,
            ny: 64,
            diffusivity: 0.01,
            ux: 0.0,
            uy: 0.0,
        }
    }
}

const CX5: [i64; 5] = [0, 1, 0, -1, 0];
const CY5: [i64; 5] = [0, 0, 1, 0, -1];
const W5: [f64; 5] = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];

/// D2Q5 lattice with `tau = 0.5 + 3 D` and a uniform advecting velocity.
#[derive(Clone, Debug)]
pub struct PassiveScalarD2Q5 {
    config: PassiveScalarConfig,
    cells: usize,
    tau: f64,
    g: Vec<f64>,
    next: Vec<f64>,
}

/// Temperature is the same conserved scalar equation; no buoyancy feedback is implied.
pub type ThermalD2Q5 = PassiveScalarD2Q5;

impl PassiveScalarD2Q5 {
    pub fn new(config: PassiveScalarConfig) -> Result<Self> {
        require(
            config.nx >= 2 && config.ny >= 2 && config.diffusivity > 0.0,
            "invalid passive scalar lattice",
        )?;
        let cells = config.nx * config.ny;
        let mut s = Self {
            config,
            cells,
            tau: 0.5 + 3.0 * config.diffusivity,
            g: vec![0.0; 5 * cells],
            next: vec![0.0; 5 * cells],
        };
        s.initialize(|_, _| 0.0);
        Ok(s)
    }
    fn equilibrium(&self, q: usize, phi: f64) -> f64 {
        W5[q]
            * phi
            * (1.0 + 3.0 * (CX5[q] as f64 * self.config.ux + CY5[q] as f64 * self.config.uy))
    }
    /// Sets the equilibrium for `value(x, y)` at normalised cell centres `((i+0.5)/nx, (j+0.5)/ny)`.
    pub fn initialize(&mut self, value: impl Fn(f64, f64) -> f64) {
        let (nx, ny) = (self.config.nx, self.config.ny);
        for n in 0..self.cells {
            let (x, y) = (n % nx, n / nx);
            let phi = value((x as f64 + 0.5) / nx as f64, (y as f64 + 0.5) / ny as f64);
            for q in 0..5 {
                self.g[q * self.cells + n] = self.equilibrium(q, phi);
            }
        }
    }
    pub fn step(&mut self, count: usize) {
        let (nx, ny) = (self.config.nx as i64, self.config.ny as i64);
        for _ in 0..count {
            let omega = 1.0 / self.tau;
            for n in 0..self.cells {
                let (x, y) = (n as i64 % nx, n as i64 / nx);
                let mut fin = [0.0; 5];
                let mut phi = 0.0;
                for q in 0..5 {
                    let xs = (x - CX5[q] + nx) % nx;
                    let ys = (y - CY5[q] + ny) % ny;
                    fin[q] = self.g[q * self.cells + (ys * nx + xs) as usize];
                    phi += fin[q];
                }
                for q in 0..5 {
                    self.next[q * self.cells + n] =
                        fin[q] - omega * (fin[q] - self.equilibrium(q, phi));
                }
            }
            std::mem::swap(&mut self.g, &mut self.next);
        }
    }
    pub fn scalar(&self) -> Vec<f64> {
        (0..self.cells)
            .map(|n| (0..5).map(|q| self.g[q * self.cells + n]).sum())
            .collect()
    }
    pub fn total(&self) -> f64 {
        self.scalar().iter().sum()
    }
    pub fn config(&self) -> &PassiveScalarConfig {
        &self.config
    }
}
