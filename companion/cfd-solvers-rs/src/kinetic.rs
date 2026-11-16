//! Positivity-preserving 1x1v Vlasov and periodic Vlasov-Poisson baselines, ported from
//! upstream `cfd/kinetic/{vlasov1d1v, vlasov_poisson1d1v}` (MIT, cfd_solvers contributors).
use crate::{require, Complex64, Result};

#[derive(Clone, Copy, Debug)]
struct PhaseGrid {
    nx: usize,
    nv: usize,
    length: f64,
    vmax: f64,
}

impl PhaseGrid {
    fn dx(&self) -> f64 {
        self.length / self.nx as f64
    }
    fn dv(&self) -> f64 {
        2.0 * self.vmax / self.nv as f64
    }
    fn velocity(&self, j: usize) -> f64 {
        -self.vmax + (j as f64 + 0.5) * self.dv()
    }
    fn fill(&self, f: &mut [f64], init: impl Fn(f64, f64) -> f64) -> Result<()> {
        for i in 0..self.nx {
            let x = (i as f64 + 0.5) * self.dx();
            for j in 0..self.nv {
                let v = init(x, self.velocity(j));
                require(v.is_finite(), "non-finite Vlasov initial value")?;
                f[i * self.nv + j] = v.max(0.0);
            }
        }
        Ok(())
    }
    /// First-order upwind update; `accel(i)` is the velocity-space acceleration in cell column i.
    fn upwind(&self, f: &[f64], next: &mut [f64], dt: f64, accel: impl Fn(usize) -> f64) {
        let (dx, dv, nv) = (self.dx(), self.dv(), self.nv);
        for i in 0..self.nx {
            let left = if i == 0 { self.nx - 1 } else { i - 1 };
            let right = if i + 1 == self.nx { 0 } else { i + 1 };
            let a = accel(i);
            for j in 0..nv {
                let v = self.velocity(j);
                let here = f[i * nv + j];
                let xl = if v >= 0.0 {
                    v * f[left * nv + j]
                } else {
                    v * here
                };
                let xr = if v >= 0.0 {
                    v * here
                } else {
                    v * f[right * nv + j]
                };
                let vl = if j == 0 {
                    0.0
                } else if a >= 0.0 {
                    a * f[i * nv + j - 1]
                } else {
                    a * here
                };
                let vr = if j + 1 == nv {
                    0.0
                } else if a >= 0.0 {
                    a * here
                } else {
                    a * f[i * nv + j + 1]
                };
                next[i * nv + j] = (here - dt * ((xr - xl) / dx + (vr - vl) / dv)).max(0.0);
            }
        }
    }
    fn mass(&self, f: &[f64]) -> f64 {
        let volume = self.dx() * self.dv();
        f.iter().map(|v| v * volume).sum()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vlasov1D1VConfig {
    pub spatial_cells: usize,
    pub velocity_cells: usize,
    pub spatial_length: f64,
    pub maximum_velocity: f64,
}

impl Default for Vlasov1D1VConfig {
    fn default() -> Self {
        Self {
            spatial_cells: 64,
            velocity_cells: 64,
            spatial_length: 1.0,
            maximum_velocity: 4.0,
        }
    }
}

/// `f_t + v f_x + a f_v = 0`, periodic in x and closed in v; distribution `[x * nv + v]`.
#[derive(Clone, Debug)]
pub struct Vlasov1D1V {
    grid: PhaseGrid,
    f: Vec<f64>,
    next: Vec<f64>,
}

impl Vlasov1D1V {
    pub fn new(c: Vlasov1D1VConfig) -> Result<Self> {
        require(
            c.spatial_cells >= 3
                && c.velocity_cells >= 3
                && c.spatial_length > 0.0
                && c.maximum_velocity > 0.0,
            "invalid Vlasov 1x1v configuration",
        )?;
        let n = c.spatial_cells * c.velocity_cells;
        Ok(Self {
            grid: PhaseGrid {
                nx: c.spatial_cells,
                nv: c.velocity_cells,
                length: c.spatial_length,
                vmax: c.maximum_velocity,
            },
            f: vec![0.0; n],
            next: vec![0.0; n],
        })
    }
    /// Non-positive initial values are clipped to zero (non-finite values are rejected).
    pub fn initialize(&mut self, f: impl Fn(f64, f64) -> f64) -> Result<()> {
        self.grid.fill(&mut self.f, f)
    }
    pub fn step(&mut self, dt: f64, acceleration: f64) -> Result<()> {
        require(
            dt > 0.0 && dt.is_finite() && acceleration.is_finite(),
            "invalid Vlasov timestep",
        )?;
        self.grid
            .upwind(&self.f, &mut self.next, dt, |_| acceleration);
        std::mem::swap(&mut self.f, &mut self.next);
        Ok(())
    }
    pub fn distribution(&self) -> &[f64] {
        &self.f
    }
    pub fn mass(&self) -> f64 {
        self.grid.mass(&self.f)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VlasovPoisson1D1VConfig {
    pub spatial_cells: usize,
    pub velocity_cells: usize,
    pub spatial_length: f64,
    pub maximum_velocity: f64,
    pub neutralizing_density: f64,
}

impl Default for VlasovPoisson1D1VConfig {
    fn default() -> Self {
        Self {
            spatial_cells: 32,
            velocity_cells: 32,
            spatial_length: 1.0,
            maximum_velocity: 4.0,
            neutralizing_density: 0.0,
        }
    }
}

/// Electrostatic Vlasov-Poisson with a direct-DFT periodic field solve (`E_hat = -i rho_hat / k`).
#[derive(Clone, Debug)]
pub struct VlasovPoisson1D1V {
    grid: PhaseGrid,
    neutralizing_density: f64,
    f: Vec<f64>,
    next: Vec<f64>,
    field: Vec<f64>,
    charge: Vec<f64>,
}

impl VlasovPoisson1D1V {
    pub fn new(c: VlasovPoisson1D1VConfig) -> Result<Self> {
        require(
            c.spatial_cells >= 4
                && c.velocity_cells >= 4
                && c.spatial_length > 0.0
                && c.maximum_velocity > 0.0
                && c.neutralizing_density >= 0.0
                && c.neutralizing_density.is_finite(),
            "invalid Vlasov-Poisson configuration",
        )?;
        let n = c.spatial_cells * c.velocity_cells;
        Ok(Self {
            grid: PhaseGrid {
                nx: c.spatial_cells,
                nv: c.velocity_cells,
                length: c.spatial_length,
                vmax: c.maximum_velocity,
            },
            neutralizing_density: c.neutralizing_density,
            f: vec![0.0; n],
            next: vec![0.0; n],
            field: vec![0.0; c.spatial_cells],
            charge: vec![0.0; c.spatial_cells],
        })
    }
    pub fn initialize(&mut self, f: impl Fn(f64, f64) -> f64) -> Result<()> {
        self.grid.fill(&mut self.f, f)?;
        self.compute_field();
        Ok(())
    }
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(
            dt > 0.0 && dt.is_finite(),
            "invalid Vlasov-Poisson timestep",
        )?;
        self.compute_field();
        let field = &self.field;
        self.grid.upwind(&self.f, &mut self.next, dt, |i| field[i]);
        std::mem::swap(&mut self.f, &mut self.next);
        self.compute_field();
        Ok(())
    }
    fn compute_field(&mut self) {
        let (nx, nv, dv) = (self.grid.nx, self.grid.nv, self.grid.dv());
        for i in 0..nx {
            let mut density = 0.0;
            for j in 0..nv {
                density += self.f[i * nv + j] * dv;
            }
            self.charge[i] = density - self.neutralizing_density;
        }
        let pi = std::f64::consts::PI;
        let mut rho_hat = vec![Complex64::new(0.0, 0.0); nx];
        for (mode, r) in rho_hat.iter_mut().enumerate().skip(1) {
            for cell in 0..nx {
                let angle = -2.0 * pi * (mode * cell) as f64 / nx as f64;
                *r += self.charge[cell] * Complex64::new(angle.cos(), angle.sin());
            }
        }
        let mut field_hat = vec![Complex64::new(0.0, 0.0); nx];
        for mode in 1..nx {
            let signed = if mode <= nx / 2 {
                mode as f64
            } else {
                mode as f64 - nx as f64
            };
            let k = 2.0 * pi * signed / self.grid.length;
            field_hat[mode] = Complex64::new(0.0, -k) * rho_hat[mode] / (k * k);
        }
        for cell in 0..nx {
            let mut value = Complex64::new(0.0, 0.0);
            for (mode, fh) in field_hat.iter().enumerate().skip(1) {
                let angle = 2.0 * pi * (mode * cell) as f64 / nx as f64;
                value += fh * Complex64::new(angle.cos(), angle.sin());
            }
            self.field[cell] = value.re / nx as f64;
        }
    }
    pub fn distribution(&self) -> &[f64] {
        &self.f
    }
    pub fn electric_field(&self) -> &[f64] {
        &self.field
    }
    pub fn charge_density(&self) -> &[f64] {
        &self.charge
    }
    pub fn mass(&self) -> f64 {
        self.grid.mass(&self.f)
    }
    pub fn field_energy(&self) -> f64 {
        let dx = self.grid.dx();
        self.field.iter().map(|e| 0.5 * e * e * dx).sum()
    }
}
