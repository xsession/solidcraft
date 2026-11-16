//! Native single-precision D2Q9 pull-streaming BGK, Guo forcing and wall/inlet conditions.
mod advanced;
mod descriptors;
mod diagnostics;
mod distributed_pull;
mod free_surface;
mod ib_particles;
mod pull;
mod scalar;
mod visualization;
mod voxelize;
use crate::{require, Result};
pub use advanced::*;
pub use descriptors::*;
pub use diagnostics::*;
pub use distributed_pull::*;
pub use free_surface::*;
pub use ib_particles::*;
pub use pull::*;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
pub use scalar::*;
use serde::{Deserialize, Serialize};
pub use visualization::*;
pub use voxelize::*;
pub const CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
pub const CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
pub const OPPOSITE: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
pub const WEIGHT: [f32; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
pub fn equilibrium(q: usize, rho: f32, ux: f32, uy: f32) -> f32 {
    let cu = 3.0 * (CX[q] as f32 * ux + CY[q] as f32 * uy);
    let uu = 1.5 * (ux * ux + uy * uy);
    WEIGHT[q] * rho * (1.0 + cu + 0.5 * cu * cu - uu)
}
pub fn bouzidi_interpolated_bounce_back(
    fraction: f32,
    opposite_here: f32,
    opposite_away: f32,
    incoming_here: f32,
) -> Result<f32> {
    require(
        fraction > 0.0 && fraction <= 1.0,
        "Bouzidi fraction must be in (0,1]",
    )?;
    Ok(if fraction <= 0.5 {
        2.0 * fraction * opposite_here + (1.0 - 2.0 * fraction) * opposite_away
    } else {
        let inv = 0.5 / fraction;
        inv * opposite_here + (1.0 - inv) * incoming_here
    })
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct D2Q9Config {
    pub nx: usize,
    pub ny: usize,
    pub tau: f32,
    pub acceleration_x: f32,
    pub acceleration_y: f32,
}
impl Default for D2Q9Config {
    fn default() -> Self {
        Self {
            nx: 128,
            ny: 128,
            tau: 0.6,
            acceleration_x: 0.0,
            acceleration_y: 0.0,
        }
    }
}
#[derive(Clone, Debug)]
pub struct D2Q9Solver {
    config: D2Q9Config,
    f: Vec<[f32; 9]>,
    rho: Vec<f32>,
    ux: Vec<f32>,
    uy: Vec<f32>,
    solid: Vec<bool>,
    wall_velocity: Vec<[f32; 2]>,
    wall_fraction: Vec<[f32; 9]>,
    inlet: Option<[f32; 2]>,
    outlet: Option<[f32; 2]>,
}
impl D2Q9Solver {
    pub fn new(config: D2Q9Config) -> Result<Self> {
        require(
            config.nx > 0
                && config.ny > 0
                && config.nx <= i32::MAX as usize
                && config.ny <= i32::MAX as usize,
            "invalid LBM dimensions",
        )?;
        require(
            config.tau.is_finite()
                && config.tau > 0.5
                && config.acceleration_x.is_finite()
                && config.acceleration_y.is_finite(),
            "invalid LBM relaxation/forcing",
        )?;
        let n = config
            .nx
            .checked_mul(config.ny)
            .ok_or_else(|| crate::SolverError::InvalidInput("LBM size overflow".into()))?;
        let mut s = Self {
            config,
            f: vec![[0.0; 9]; n],
            rho: vec![1.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            solid: vec![false; n],
            wall_velocity: vec![[0.0; 2]; n],
            wall_fraction: vec![[-1.0; 9]; n],
            inlet: None,
            outlet: None,
        };
        s.initialize_uniform(1.0, 0.0, 0.0)?;
        Ok(s)
    }
    pub fn config(&self) -> D2Q9Config {
        self.config
    }
    pub fn populations(&self) -> &[[f32; 9]] {
        &self.f
    }
    pub fn density(&self) -> &[f32] {
        &self.rho
    }
    pub fn velocity_x(&self) -> &[f32] {
        &self.ux
    }
    pub fn velocity_y(&self) -> &[f32] {
        &self.uy
    }
    pub fn solid_mask(&self) -> &[bool] {
        &self.solid
    }
    fn index(&self, x: usize, y: usize) -> Result<usize> {
        require(
            x < self.config.nx && y < self.config.ny,
            "LBM coordinate outside grid",
        )?;
        Ok(y * self.config.nx + x)
    }
    pub fn initialize_uniform(&mut self, rho: f32, ux: f32, uy: f32) -> Result<()> {
        require(
            rho.is_finite() && rho > 0.0 && ux.is_finite() && uy.is_finite(),
            "invalid initial LBM field",
        )?;
        for i in 0..self.f.len() {
            let u = if self.solid[i] {
                self.wall_velocity[i]
            } else {
                [ux, uy]
            };
            self.rho[i] = rho;
            self.ux[i] = u[0];
            self.uy[i] = u[1];
            for q in 0..9 {
                self.f[i][q] = equilibrium(q, rho, u[0], u[1]);
            }
        }
        Ok(())
    }
    pub fn initialize_taylor_green(&mut self, amplitude: f32) -> Result<()> {
        require(amplitude.is_finite(), "nonfinite Taylor-Green amplitude")?;
        let two_pi = 2.0 * std::f32::consts::PI;
        for i in 0..self.f.len() {
            let x = (i % self.config.nx) as f32 + 0.5;
            let y = (i / self.config.nx) as f32 + 0.5;
            let x = two_pi * x / self.config.nx as f32;
            let y = two_pi * y / self.config.ny as f32;
            let ux = amplitude * x.sin() * y.cos();
            let uy = -amplitude * x.cos() * y.sin();
            self.rho[i] = 1.0;
            self.ux[i] = ux;
            self.uy[i] = uy;
            for q in 0..9 {
                self.f[i][q] = equilibrium(q, 1.0, ux, uy);
            }
        }
        Ok(())
    }
    pub fn set_body_acceleration(&mut self, ax: f32, ay: f32) -> Result<()> {
        require(ax.is_finite() && ay.is_finite(), "nonfinite acceleration")?;
        self.config.acceleration_x = ax;
        self.config.acceleration_y = ay;
        Ok(())
    }
    pub fn set_solid(&mut self, x: usize, y: usize, solid: bool) -> Result<()> {
        let i = self.index(x, y)?;
        self.solid[i] = solid;
        self.wall_velocity[i] = [0.0; 2];
        Ok(())
    }
    pub fn set_wall_velocity(&mut self, x: usize, y: usize, ux: f32, uy: f32) -> Result<()> {
        require(ux.is_finite() && uy.is_finite(), "nonfinite wall velocity")?;
        let i = self.index(x, y)?;
        self.solid[i] = true;
        self.wall_velocity[i] = [ux, uy];
        Ok(())
    }
    pub fn set_interpolated_wall_link(
        &mut self,
        x: usize,
        y: usize,
        q: usize,
        fraction: f32,
    ) -> Result<()> {
        let i = self.index(x, y)?;
        require(
            q > 0 && q < 9 && fraction > 0.0 && fraction <= 1.0,
            "invalid interpolated wall link",
        )?;
        self.wall_fraction[i][q] = fraction;
        Ok(())
    }
    pub fn clear_interpolated_wall_links(&mut self) {
        self.wall_fraction.fill([-1.0; 9]);
    }
    pub fn set_velocity_inlet_left(&mut self, ux: f32, uy: f32) -> Result<()> {
        require(
            ux.is_finite() && ux < 1.0 && uy.is_finite(),
            "invalid velocity inlet",
        )?;
        self.inlet = Some([ux, uy]);
        Ok(())
    }
    pub fn set_pressure_outlet_right(&mut self, rho: f32, uy: f32) -> Result<()> {
        require(
            rho.is_finite() && rho > 0.0 && uy.is_finite(),
            "invalid pressure outlet",
        )?;
        self.outlet = Some([rho, uy]);
        Ok(())
    }
    pub fn clear_x_boundaries(&mut self) {
        self.inlet = None;
        self.outlet = None;
    }
    fn update_cell(&self, i: usize) -> ([f32; 9], f32, f32, f32) {
        if self.solid[i] {
            return (self.f[i], self.rho[i], self.ux[i], self.uy[i]);
        }
        let (nx, ny) = (self.config.nx as i64, self.config.ny as i64);
        let (x, y) = ((i % self.config.nx) as i64, (i / self.config.nx) as i64);
        let mut f = [0.0; 9];
        for q in 0..9 {
            let (sx, sy) = (x - CX[q] as i64, y - CY[q] as i64);
            if (self.inlet.is_some() && x == 0 && sx < 0)
                || (self.outlet.is_some() && x + 1 == nx && sx >= nx)
            {
                f[q] = self.f[i][q];
                continue;
            }
            let source = (sy.rem_euclid(ny) * nx + sx.rem_euclid(nx)) as usize;
            if self.solid[source] {
                let opposite = OPPOSITE[q];
                let fraction = self.wall_fraction[i][q];
                let mut value = self.f[i][opposite];
                if fraction > 0.0 {
                    let away = ((y + CY[q] as i64).rem_euclid(ny) * nx
                        + (x + CX[q] as i64).rem_euclid(nx))
                        as usize;
                    let oa = if self.solid[away] {
                        value
                    } else {
                        self.f[away][opposite]
                    };
                    value = if fraction <= 0.5 {
                        2.0 * fraction * value + (1.0 - 2.0 * fraction) * oa
                    } else {
                        let inv = 0.5 / fraction;
                        inv * value + (1.0 - inv) * self.f[i][q]
                    };
                }
                let wall_dot = CX[q] as f32 * self.wall_velocity[source][0]
                    + CY[q] as f32 * self.wall_velocity[source][1];
                f[q] = value + 6.0 * WEIGHT[q] * self.rho[i].max(1e-8) * wall_dot;
            } else {
                f[q] = self.f[source][q];
            }
        }
        if x == 0 {
            if let Some([ux, uy]) = self.inlet {
                let r = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / (1.0 - ux);
                f[1] = f[3] + (2.0 / 3.0) * r * ux;
                f[5] = f[7] + 0.5 * (f[4] - f[2]) + (1.0 / 6.0) * r * ux + 0.5 * r * uy;
                f[8] = f[6] + 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * r * ux - 0.5 * r * uy;
            }
        }
        if x + 1 == nx {
            if let Some([r, uy]) = self.outlet {
                let ux = -1.0 + (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / r;
                f[3] = f[1] - (2.0 / 3.0) * r * ux;
                f[6] = f[8] + 0.5 * (f[4] - f[2]) - (1.0 / 6.0) * r * ux + 0.5 * r * uy;
                f[7] = f[5] + 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * r * ux - 0.5 * r * uy;
            }
        }
        let (mut rho, mut mx, mut my) = (0.0, 0.0, 0.0);
        for q in 0..9 {
            rho += f[q];
            mx += f[q] * CX[q] as f32;
            my += f[q] * CY[q] as f32;
        }
        let (ax, ay) = (self.config.acceleration_x, self.config.acceleration_y);
        let (mut ux, mut uy) = if rho > 0.0 {
            (mx / rho + 0.5 * ax, my / rho + 0.5 * ay)
        } else {
            (0.0, 0.0)
        };
        if x == 0 {
            if let Some(u) = self.inlet {
                ux = u[0];
                uy = u[1];
            }
        }
        if x + 1 == nx {
            if let Some(u) = self.outlet {
                uy = u[1];
            }
        }
        let omega = 1.0 / self.config.tau;
        for q in 0..9 {
            let (cx, cy) = (CX[q] as f32, CY[q] as f32);
            let cu = cx * ux + cy * uy;
            let cf = cx * rho * ax + cy * rho * ay;
            let uf = ux * rho * ax + uy * rho * ay;
            let source = if ax == 0.0 && ay == 0.0 {
                0.0
            } else {
                WEIGHT[q] * (1.0 - 0.5 * omega) * (3.0 * (cf - uf) + 9.0 * cu * cf)
            };
            f[q] = f[q] - omega * (f[q] - equilibrium(q, rho, ux, uy)) + source;
        }
        (f, rho, ux, uy)
    }
    pub fn step(&mut self, count: usize) -> Result<()> {
        self.advance(count, true)
    }
    pub fn step_serial(&mut self, count: usize) -> Result<()> {
        self.advance(count, false)
    }
    fn advance(&mut self, count: usize, _allow_parallel: bool) -> Result<()> {
        for _ in 0..count {
            #[cfg(feature = "parallel")]
            let values: Vec<_> = if _allow_parallel && self.f.len() >= 8192 {
                (0..self.f.len())
                    .into_par_iter()
                    .map(|i| self.update_cell(i))
                    .collect()
            } else {
                (0..self.f.len()).map(|i| self.update_cell(i)).collect()
            };
            #[cfg(not(feature = "parallel"))]
            let values: Vec<_> = (0..self.f.len()).map(|i| self.update_cell(i)).collect();
            require(
                values.iter().all(|(f, r, x, y)| {
                    r.is_finite()
                        && *r > 0.0
                        && x.is_finite()
                        && y.is_finite()
                        && f.iter().all(|v| v.is_finite())
                }),
                "LBM became nonfinite or nonpositive",
            )?;
            for (i, (f, r, x, y)) in values.into_iter().enumerate() {
                self.f[i] = f;
                self.rho[i] = r;
                self.ux[i] = x;
                self.uy[i] = y;
            }
        }
        self.compute_macroscopic();
        Ok(())
    }
    pub fn compute_macroscopic(&mut self) {
        for i in 0..self.f.len() {
            if self.solid[i] {
                self.rho[i] = 1.0;
                self.ux[i] = self.wall_velocity[i][0];
                self.uy[i] = self.wall_velocity[i][1];
                continue;
            }
            let (mut r, mut mx, mut my) = (0.0, 0.0, 0.0);
            for q in 0..9 {
                let v = self.f[i][q];
                r += v;
                mx += v * CX[q] as f32;
                my += v * CY[q] as f32;
            }
            self.rho[i] = r;
            self.ux[i] = if r > 0.0 {
                mx / r + 0.5 * self.config.acceleration_x
            } else {
                0.0
            };
            self.uy[i] = if r > 0.0 {
                my / r + 0.5 * self.config.acceleration_y
            } else {
                0.0
            };
            let x = i % self.config.nx;
            if x == 0 {
                if let Some(u) = self.inlet {
                    self.ux[i] = u[0];
                    self.uy[i] = u[1];
                }
            }
            if x + 1 == self.config.nx {
                if let Some(u) = self.outlet {
                    self.uy[i] = u[1];
                }
            }
        }
    }
    pub fn mass(&self) -> f64 {
        (0..self.f.len())
            .filter(|&i| !self.solid[i])
            .map(|i| self.rho[i] as f64)
            .sum()
    }
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * (0..self.f.len())
            .filter(|&i| !self.solid[i])
            .map(|i| {
                self.rho[i] as f64 * ((self.ux[i] as f64).powi(2) + (self.uy[i] as f64).powi(2))
            })
            .sum::<f64>()
    }
    pub fn max_speed(&self) -> f32 {
        (0..self.f.len())
            .filter(|&i| !self.solid[i])
            .map(|i| self.ux[i].hypot(self.uy[i]))
            .fold(0.0, f32::max)
    }
    pub fn average_ux(&self, xb: usize, xe: usize, yb: usize, ye: usize) -> Result<f64> {
        require(
            xb < xe && yb < ye && xe <= self.config.nx && ye <= self.config.ny,
            "LBM averaging bounds",
        )?;
        let (mut sum, mut n) = (0.0, 0);
        for y in yb..ye {
            for x in xb..xe {
                let i = y * self.config.nx + x;
                if !self.solid[i] {
                    sum += self.ux[i] as f64;
                    n += 1;
                }
            }
        }
        Ok(if n == 0 { 0.0 } else { sum / n as f64 })
    }
}
