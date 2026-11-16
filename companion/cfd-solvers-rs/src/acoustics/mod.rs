//! 2D FFT k-space acoustics, split-field damping, power-law loss and sensor time reversal.
mod coupling;
mod postprocess;
use crate::{
    core::{fft2d, wavenumbers},
    finite_slice, positive, require, Complex64 as C, Result,
};
pub use coupling::*;
pub use postprocess::*;
use serde::{Deserialize, Serialize};
use std::{f64::consts::PI, rc::Rc};
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct KSpace2DConfig {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    pub dt: f64,
    pub cfl: f64,
    pub reference_sound_speed: f64,
    pub pml_cells: usize,
    pub pml_order: f64,
    pub pml_target_reflection: f64,
}
impl Default for KSpace2DConfig {
    fn default() -> Self {
        Self {
            nx: 128,
            ny: 128,
            dx: 1e-4,
            dy: 1e-4,
            dt: 0.0,
            cfl: 0.2,
            reference_sound_speed: 0.0,
            pml_cells: 12,
            pml_order: 3.0,
            pml_target_reflection: 1e-8,
        }
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PowerLawAbsorption2D {
    pub alpha0_np_per_m: f64,
    pub exponent: f64,
    pub reference_frequency_hz: f64,
}
impl Default for PowerLawAbsorption2D {
    fn default() -> Self {
        Self {
            alpha0_np_per_m: 0.0,
            exponent: 1.5,
            reference_frequency_hz: 1e6,
        }
    }
}
impl PowerLawAbsorption2D {
    fn validate(self) -> Result<()> {
        require(
            self.alpha0_np_per_m.is_finite() && self.alpha0_np_per_m >= 0.0,
            "invalid absorption coefficient",
        )?;
        positive(self.exponent, "absorption exponent")?;
        positive(self.reference_frequency_hz, "reference frequency")
    }
}
#[derive(Clone, Debug)]
pub struct AcousticSensor2D {
    pub cell: usize,
    pub pressure: Vec<f64>,
}
fn spectrum(
    field: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
) -> Result<(Vec<C>, Vec<f64>, Vec<f64>)> {
    require(
        nx.checked_mul(ny) == Some(field.len()) && nx.is_power_of_two() && ny.is_power_of_two(),
        "spectral grid dimensions",
    )?;
    finite_slice(field, "spectral input")?;
    let kx = wavenumbers(nx, dx)?;
    let ky = wavenumbers(ny, dy)?;
    let mut f: Vec<_> = field.iter().map(|&v| C::new(v, 0.0)).collect();
    fft2d(&mut f, nx, ny, false)?;
    Ok((f, kx, ky))
}
pub fn fractional_laplacian_periodic_2d(
    field: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    order: f64,
) -> Result<Vec<f64>> {
    require(
        order.is_finite() && order >= 0.0,
        "invalid fractional order",
    )?;
    let (mut f, kx, ky) = spectrum(field, nx, ny, dx, dy)?;
    for j in 0..ny {
        for i in 0..nx {
            let k2 = kx[i] * kx[i] + ky[j] * ky[j];
            let factor = if k2 == 0.0 {
                if order == 0.0 {
                    1.0
                } else {
                    0.0
                }
            } else {
                k2.powf(order)
            };
            f[j * nx + i] *= factor;
        }
    }
    fft2d(&mut f, nx, ny, true)?;
    Ok(f.iter().map(|v| v.re).collect())
}
pub fn power_law_attenuate_periodic_2d(
    field: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    c_ref: f64,
    dt: f64,
    a: PowerLawAbsorption2D,
) -> Result<Vec<f64>> {
    positive(c_ref, "reference sound speed")?;
    require(dt.is_finite() && dt >= 0.0, "invalid attenuation timestep")?;
    a.validate()?;
    let (mut f, kx, ky) = spectrum(field, nx, ny, dx, dy)?;
    for j in 0..ny {
        for i in 0..nx {
            let frequency = c_ref * kx[i].hypot(ky[j]) / (2.0 * PI);
            let alpha = a.alpha0_np_per_m * (frequency / a.reference_frequency_hz).powf(a.exponent);
            f[j * nx + i] *= (-alpha * c_ref * dt).exp();
        }
    }
    fft2d(&mut f, nx, ny, true)?;
    Ok(f.iter().map(|v| v.re).collect())
}
type PressureSource<'a> = (usize, Rc<dyn Fn(f64) -> f64 + 'a>);
pub struct KSpaceAcoustic2D<'a> {
    config: KSpace2DConfig,
    dt: f64,
    time: f64,
    c0: Vec<f64>,
    rho0: Vec<f64>,
    bon_a: Vec<f64>,
    ux: Vec<f64>,
    uy: Vec<f64>,
    rho_x: Vec<f64>,
    rho_y: Vec<f64>,
    rho: Vec<f64>,
    pressure: Vec<f64>,
    damp_x: Vec<f64>,
    damp_y: Vec<f64>,
    kx: Vec<f64>,
    ky: Vec<f64>,
    kappa: Vec<f64>,
    sources: Vec<PressureSource<'a>>,
    sensors: Vec<AcousticSensor2D>,
    absorption: Option<PowerLawAbsorption2D>,
}
impl<'a> KSpaceAcoustic2D<'a> {
    pub fn new(config: KSpace2DConfig) -> Result<Self> {
        require(
            config.nx.is_power_of_two()
                && config.ny.is_power_of_two()
                && config.pml_cells < config.nx.min(config.ny).div_ceil(2),
            "invalid acoustic dimensions/layer",
        )?;
        positive(config.dx, "dx")?;
        positive(config.dy, "dy")?;
        positive(config.cfl, "CFL")?;
        positive(config.pml_order, "PML order")?;
        require(
            config.dt.is_finite()
                && config.dt >= 0.0
                && config.reference_sound_speed.is_finite()
                && config.reference_sound_speed >= 0.0
                && config.pml_target_reflection > 0.0
                && config.pml_target_reflection < 1.0,
            "invalid acoustic timestep/PML",
        )?;
        let n = config
            .nx
            .checked_mul(config.ny)
            .ok_or_else(|| crate::SolverError::InvalidInput("acoustic size overflow".into()))?;
        let mut s = Self {
            config,
            dt: 0.0,
            time: 0.0,
            c0: vec![1500.0; n],
            rho0: vec![1000.0; n],
            bon_a: vec![0.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            rho_x: vec![0.0; n],
            rho_y: vec![0.0; n],
            rho: vec![0.0; n],
            pressure: vec![0.0; n],
            damp_x: vec![1.0; config.nx],
            damp_y: vec![1.0; config.ny],
            kx: wavenumbers(config.nx, config.dx)?,
            ky: wavenumbers(config.ny, config.dy)?,
            kappa: vec![0.0; n],
            sources: vec![],
            sensors: vec![],
            absorption: None,
        };
        s.rebuild()?;
        Ok(s)
    }
    pub fn config(&self) -> KSpace2DConfig {
        self.config
    }
    pub fn dt(&self) -> f64 {
        self.dt
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn cells(&self) -> usize {
        self.pressure.len()
    }
    pub fn pressure(&self) -> &[f64] {
        &self.pressure
    }
    pub fn velocity_x(&self) -> &[f64] {
        &self.ux
    }
    pub fn velocity_y(&self) -> &[f64] {
        &self.uy
    }
    pub fn density_perturbation(&self) -> &[f64] {
        &self.rho
    }
    pub fn sound_speed(&self) -> &[f64] {
        &self.c0
    }
    pub fn medium_density(&self) -> &[f64] {
        &self.rho0
    }
    pub fn sensors(&self) -> &[AcousticSensor2D] {
        &self.sensors
    }
    fn rebuild(&mut self) -> Result<()> {
        let cmax = self.c0.iter().copied().fold(0.0, f64::max);
        let cref = if self.config.reference_sound_speed > 0.0 {
            self.config.reference_sound_speed
        } else {
            cmax
        };
        let spacing = self.config.dx.min(self.config.dy);
        let dt = if self.config.dt > 0.0 {
            self.config.dt
        } else {
            self.config.cfl * spacing / cmax
        };
        require(
            dt.is_finite() && dt > 0.0 && dt * cmax / spacing <= 1.0,
            "unstable acoustic timestep",
        )?;
        self.dt = dt;
        for j in 0..self.config.ny {
            for i in 0..self.config.nx {
                let x = 0.5 * cref * self.kx[i].hypot(self.ky[j]) * dt;
                self.kappa[j * self.config.nx + i] =
                    if x.abs() < 1e-12 { 1.0 } else { x.sin() / x };
            }
        }
        self.damp_x.fill(1.0);
        self.damp_y.fill(1.0);
        let layer = self.config.pml_cells;
        if layer > 0 {
            let sigma =
                -(self.config.pml_order + 1.0) * cmax * self.config.pml_target_reflection.ln()
                    / (2.0 * layer as f64 * spacing);
            for d in [&mut self.damp_x, &mut self.damp_y] {
                let n = d.len();
                for i in 0..n {
                    let depth = i.min(n - 1 - i);
                    if depth < layer {
                        let s = (layer - depth) as f64 / layer as f64;
                        d[i] = (-sigma * s.powf(self.config.pml_order) * dt).exp();
                    }
                }
            }
        }
        Ok(())
    }
    pub fn set_uniform_medium(&mut self, c: f64, rho: f64) -> Result<()> {
        positive(c, "sound speed")?;
        positive(rho, "density")?;
        self.c0.fill(c);
        self.rho0.fill(rho);
        self.rebuild()
    }
    pub fn set_medium(&mut self, c: &[f64], rho: &[f64]) -> Result<()> {
        require(
            c.len() == self.cells()
                && rho.len() == self.cells()
                && c.iter().chain(rho).all(|v| v.is_finite() && *v > 0.0),
            "invalid acoustic medium",
        )?;
        self.c0.copy_from_slice(c);
        self.rho0.copy_from_slice(rho);
        self.rebuild()
    }
    pub fn set_nonlinearity(&mut self, b: &[f64]) -> Result<()> {
        require(b.len() == self.cells(), "BonA size")?;
        finite_slice(b, "BonA")?;
        self.bon_a.copy_from_slice(b);
        Ok(())
    }
    pub fn set_uniform_nonlinearity(&mut self, b: f64) -> Result<()> {
        require(b.is_finite(), "nonfinite BonA")?;
        self.bon_a.fill(b);
        Ok(())
    }
    pub fn set_power_law_absorption(&mut self, a: PowerLawAbsorption2D) -> Result<()> {
        a.validate()?;
        self.absorption = if a.alpha0_np_per_m > 0.0 {
            Some(a)
        } else {
            None
        };
        Ok(())
    }
    pub fn clear_absorption(&mut self) {
        self.absorption = None;
    }
    pub fn initialize_pressure(&mut self, p: &[f64]) -> Result<()> {
        require(p.len() == self.cells(), "initial pressure size")?;
        finite_slice(p, "initial pressure")?;
        self.pressure.copy_from_slice(p);
        for q in 0..self.cells() {
            self.rho[q] = p[q] / self.c0[q].powi(2);
            self.rho_x[q] = 0.5 * self.rho[q];
            self.rho_y[q] = 0.5 * self.rho[q];
        }
        self.ux.fill(0.0);
        self.uy.fill(0.0);
        self.time = 0.0;
        self.record_sensors();
        Ok(())
    }
    pub fn initialize_gaussian_pressure(
        &mut self,
        cx: f64,
        cy: f64,
        sigma: f64,
        amplitude: f64,
    ) -> Result<()> {
        positive(sigma, "Gaussian sigma")?;
        require(
            [cx, cy, amplitude].iter().all(|v| v.is_finite()),
            "invalid Gaussian parameters",
        )?;
        let x0 = cx * self.config.nx as f64 * self.config.dx;
        let y0 = cy * self.config.ny as f64 * self.config.dy;
        let mut p = vec![0.0; self.cells()];
        for j in 0..self.config.ny {
            for i in 0..self.config.nx {
                let x = (i as f64 + 0.5) * self.config.dx;
                let y = (j as f64 + 0.5) * self.config.dy;
                p[j * self.config.nx + i] = amplitude
                    * (-0.5 * ((x - x0).powi(2) + (y - y0).powi(2)) / (sigma * sigma)).exp();
            }
        }
        self.initialize_pressure(&p)
    }
    pub fn add_pressure_source(
        &mut self,
        cell: usize,
        value: impl Fn(f64) -> f64 + 'a,
    ) -> Result<()> {
        require(cell < self.cells(), "acoustic source index")?;
        self.sources.push((cell, Rc::new(value)));
        Ok(())
    }
    pub fn add_plane_source_x(&mut self, i: usize, value: impl Fn(f64) -> f64 + 'a) -> Result<()> {
        require(i < self.config.nx, "acoustic plane source index")?;
        let f: Rc<dyn Fn(f64) -> f64> = Rc::new(value);
        for j in 0..self.config.ny {
            self.sources.push((j * self.config.nx + i, f.clone()));
        }
        Ok(())
    }
    pub fn add_delayed_source_array(
        &mut self,
        cells: &[usize],
        delays: &[f64],
        waveform: impl Fn(f64) -> f64 + 'a,
    ) -> Result<()> {
        require(
            !cells.is_empty()
                && cells.len() == delays.len()
                && cells.iter().all(|&q| q < self.cells())
                && delays.iter().all(|v| v.is_finite() && *v >= 0.0),
            "invalid delayed sources",
        )?;
        let f = Rc::new(waveform);
        for (&q, &delay) in cells.iter().zip(delays) {
            let waveform = f.clone();
            self.sources.push((
                q,
                Rc::new(move |t| if t >= delay { waveform(t - delay) } else { 0.0 }),
            ));
        }
        Ok(())
    }
    pub fn add_sensor(&mut self, cell: usize) -> Result<usize> {
        require(cell < self.cells(), "sensor index")?;
        self.sensors.push(AcousticSensor2D {
            cell,
            pressure: vec![],
        });
        Ok(self.sensors.len() - 1)
    }
    pub fn add_sensor_array(&mut self, cells: &[usize]) -> Result<Vec<usize>> {
        require(cells.iter().all(|&q| q < self.cells()), "sensor index")?;
        cells.iter().map(|&q| self.add_sensor(q)).collect()
    }
    pub fn clear_sensors(&mut self) {
        self.sensors.clear();
    }
    fn record_sensors(&mut self) {
        for s in &mut self.sensors {
            s.pressure.push(self.pressure[s.cell]);
        }
    }
    fn gradient(&self, input: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let mut f: Vec<_> = input.iter().map(|&v| C::new(v, 0.0)).collect();
        fft2d(&mut f, nx, ny, false)?;
        let mut fx = f.clone();
        let mut fy = f;
        for j in 0..ny {
            for i in 0..nx {
                let q = j * nx + i;
                let kx = if i == nx / 2 { 0.0 } else { self.kx[i] };
                let ky = if j == ny / 2 { 0.0 } else { self.ky[j] };
                fx[q] *= C::new(0.0, kx * self.kappa[q]);
                fy[q] *= C::new(0.0, ky * self.kappa[q]);
            }
        }
        fft2d(&mut fx, nx, ny, true)?;
        fft2d(&mut fy, nx, ny, true)?;
        Ok((
            fx.iter().map(|v| v.re).collect(),
            fy.iter().map(|v| v.re).collect(),
        ))
    }
    fn update_pressure(&mut self) {
        for q in 0..self.cells() {
            self.rho[q] = self.rho_x[q] + self.rho_y[q];
            self.pressure[q] = self.c0[q].powi(2)
                * self.rho[q]
                * (1.0 + 0.5 * self.bon_a[q] * self.rho[q] / self.rho0[q]);
        }
    }
    fn advance_fields(&mut self) -> Result<()> {
        let (gx, gy) = self.gradient(&self.pressure)?;
        let nx = self.config.nx;
        for j in 0..self.config.ny {
            for i in 0..nx {
                let q = j * nx + i;
                self.ux[q] = (self.ux[q] - self.dt * gx[q] / self.rho0[q]) * self.damp_x[i];
                self.uy[q] = (self.uy[q] - self.dt * gy[q] / self.rho0[q]) * self.damp_y[j];
            }
        }
        let (dux, _) = self.gradient(&self.ux)?;
        let (_, duy) = self.gradient(&self.uy)?;
        for j in 0..self.config.ny {
            for i in 0..nx {
                let q = j * nx + i;
                self.rho_x[q] = (self.rho_x[q] - self.dt * self.rho0[q] * dux[q]) * self.damp_x[i];
                self.rho_y[q] = (self.rho_y[q] - self.dt * self.rho0[q] * duy[q]) * self.damp_y[j];
            }
        }
        self.update_pressure();
        Ok(())
    }
    pub fn step(&mut self, count: usize) -> Result<()> {
        for _ in 0..count {
            self.advance_fields()?;
            if let Some(a) = self.absorption {
                let cref = if self.config.reference_sound_speed > 0.0 {
                    self.config.reference_sound_speed
                } else {
                    self.c0.iter().copied().fold(0.0, f64::max)
                };
                self.pressure = power_law_attenuate_periodic_2d(
                    &self.pressure,
                    self.config.nx,
                    self.config.ny,
                    self.config.dx,
                    self.config.dy,
                    cref,
                    self.dt,
                    a,
                )?;
                for q in 0..self.cells() {
                    self.rho[q] = self.pressure[q] / self.c0[q].powi(2);
                    self.rho_x[q] = 0.5 * self.rho[q];
                    self.rho_y[q] = 0.5 * self.rho[q];
                }
            }
            let time = self.time + self.dt;
            for (q, f) in &self.sources {
                let dp = f(time);
                require(dp.is_finite(), "nonfinite acoustic source")?;
                self.pressure[*q] += dp;
                let dr = dp / self.c0[*q].powi(2);
                self.rho[*q] += dr;
                self.rho_x[*q] += 0.5 * dr;
                self.rho_y[*q] += 0.5 * dr;
            }
            finite_slice(&self.pressure, "acoustic pressure")?;
            self.time = time;
            self.record_sensors();
        }
        Ok(())
    }
    pub fn time_reversal_reconstruct(
        &mut self,
        cells: &[usize],
        traces: &[Vec<f64>],
        positivity: bool,
    ) -> Result<Vec<f64>> {
        require(
            !cells.is_empty()
                && cells.len() == traces.len()
                && cells.iter().all(|&q| q < self.cells()),
            "time reversal sensors",
        )?;
        let samples = traces[0].len();
        require(
            samples >= 2
                && traces
                    .iter()
                    .all(|t| t.len() == samples && t.iter().all(|v| v.is_finite())),
            "time reversal trace size/values",
        )?;
        for f in [
            &mut self.ux,
            &mut self.uy,
            &mut self.rho_x,
            &mut self.rho_y,
            &mut self.rho,
            &mut self.pressure,
        ] {
            f.fill(0.0);
        }
        self.time = 0.0;
        for n in 0..samples {
            self.advance_fields()?;
            for (k, &q) in cells.iter().enumerate() {
                let p = traces[k][samples - 1 - n];
                self.pressure[q] = p;
                let dr = p / self.c0[q].powi(2);
                self.rho[q] = dr;
                self.rho_x[q] = 0.5 * dr;
                self.rho_y[q] = 0.5 * dr;
            }
            self.time += self.dt;
        }
        if positivity {
            for p in &mut self.pressure {
                *p = p.max(0.0);
            }
        }
        finite_slice(&self.pressure, "time reversal pressure")?;
        Ok(self.pressure.clone())
    }
    pub fn total_energy(&self) -> f64 {
        let dv = self.config.dx * self.config.dy;
        (0..self.cells())
            .map(|q| {
                0.5 * (self.rho0[q] * (self.ux[q] * self.ux[q] + self.uy[q] * self.uy[q])
                    + self.pressure[q] * self.pressure[q]
                        / (self.rho0[q] * self.c0[q] * self.c0[q]))
                    * dv
            })
            .sum()
    }
}
