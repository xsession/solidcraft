use super::{EPSILON0, MU0};
use crate::{finite_slice, positive, require, Result};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Boundary1D {
    Pec,
    Mur1,
    Pml,
    Cpml,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DispersionModel1D {
    None,
    Debye,
    Drude,
    Lorentz,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Maxwell1DConfig {
    pub cells: usize,
    pub dx: f64,
    pub courant: f64,
    pub epsilon_r: f64,
    pub mu_r: f64,
    pub boundary: Boundary1D,
    pub pml_cells: usize,
    pub pml_order: f64,
    pub pml_target_reflection: f64,
    pub cpml_kappa_max: f64,
    pub cpml_alpha_fraction: f64,
}
impl Default for Maxwell1DConfig {
    fn default() -> Self {
        Self {
            cells: 1024,
            dx: 1e-3,
            courant: 0.99,
            epsilon_r: 1.0,
            mu_r: 1.0,
            boundary: Boundary1D::Pec,
            pml_cells: 16,
            pml_order: 3.0,
            pml_target_reflection: 1e-8,
            cpml_kappa_max: 5.0,
            cpml_alpha_fraction: 0.0,
        }
    }
}
#[derive(Clone, Debug)]
struct Cpml {
    k: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    psi: Vec<f64>,
}
impl Cpml {
    fn new(n: usize) -> Self {
        Self {
            k: vec![1.0; n],
            b: vec![1.0; n],
            c: vec![0.0; n],
            psi: vec![0.0; n],
        }
    }
    fn apply(&mut self, i: usize, derivative: f64) -> f64 {
        self.psi[i] = self.b[i] * self.psi[i] + self.c[i] * derivative;
        derivative / self.k[i] + self.psi[i]
    }
}
struct Tfsf<'a> {
    cell: usize,
    incident: Box<dyn Fn(f64) -> f64 + 'a>,
    impedance: f64,
}
pub struct Maxwell1D<'a> {
    config: Maxwell1DConfig,
    dt: f64,
    mu: f64,
    time: f64,
    epsilon: Vec<f64>,
    conductivity: Vec<f64>,
    e_decay: Vec<f64>,
    e_curl: Vec<f64>,
    pml_sigma: Vec<f64>,
    h_decay: Vec<f64>,
    h_curl: Vec<f64>,
    ez: Vec<f64>,
    hy: Vec<f64>,
    dispersion: Vec<DispersionModel1D>,
    polarization: Vec<f64>,
    previous_p: Vec<f64>,
    ade_a: Vec<f64>,
    ade_b: Vec<f64>,
    ade_c: Vec<f64>,
    cpml_e: Cpml,
    cpml_h: Cpml,
    lumped_eps: Vec<f64>,
    lumped_sigma: Vec<f64>,
    inductor_drive: Vec<f64>,
    inductor_current: Vec<f64>,
    tfsf: Option<Tfsf<'a>>,
}
impl<'a> Maxwell1D<'a> {
    pub fn new(c: Maxwell1DConfig) -> Result<Self> {
        require(c.cells >= 4, "Maxwell1D needs at least 4 cells")?;
        positive(c.dx, "dx")?;
        positive(c.epsilon_r, "epsilon_r")?;
        positive(c.mu_r, "mu_r")?;
        require(
            c.courant > 0.0 && c.courant < 1.0,
            "Courant must be in (0,1)",
        )?;
        if matches!(c.boundary, Boundary1D::Pml | Boundary1D::Cpml) {
            require(
                c.pml_cells >= 2
                    && c.pml_cells <= ((c.cells - 4) / 2)
                    && c.pml_order.is_finite()
                    && c.pml_order >= 1.0
                    && c.pml_target_reflection > 0.0
                    && c.pml_target_reflection < 1.0,
                "invalid absorbing layer",
            )?;
        }
        if c.boundary == Boundary1D::Cpml {
            require(
                c.cpml_kappa_max.is_finite()
                    && c.cpml_kappa_max >= 1.0
                    && c.cpml_alpha_fraction.is_finite()
                    && c.cpml_alpha_fraction >= 0.0,
                "invalid CPML controls",
            )?;
        }
        let n = c.cells;
        let mu = MU0 * c.mu_r;
        let dt = c.courant * c.dx * (EPSILON0 * c.epsilon_r * mu).sqrt();
        let mut s = Self {
            config: c,
            dt,
            mu,
            time: 0.0,
            epsilon: vec![EPSILON0 * c.epsilon_r; n],
            conductivity: vec![0.0; n],
            e_decay: vec![1.0; n],
            e_curl: vec![0.0; n],
            pml_sigma: vec![0.0; n],
            h_decay: vec![1.0; n - 1],
            h_curl: vec![0.0; n - 1],
            ez: vec![0.0; n],
            hy: vec![0.0; n - 1],
            dispersion: vec![DispersionModel1D::None; n],
            polarization: vec![0.0; n],
            previous_p: vec![0.0; n],
            ade_a: vec![0.0; n],
            ade_b: vec![0.0; n],
            ade_c: vec![0.0; n],
            cpml_e: Cpml::new(n),
            cpml_h: Cpml::new(n - 1),
            lumped_eps: vec![0.0; n],
            lumped_sigma: vec![0.0; n],
            inductor_drive: vec![0.0; n],
            inductor_current: vec![0.0; n],
            tfsf: None,
        };
        s.initialize_profiles();
        s.preprocess();
        Ok(s)
    }
    pub fn config(&self) -> Maxwell1DConfig {
        self.config
    }
    pub fn dt(&self) -> f64 {
        self.dt
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn electric(&self) -> &[f64] {
        &self.ez
    }
    pub fn magnetic(&self) -> &[f64] {
        &self.hy
    }
    pub fn polarization(&self) -> &[f64] {
        &self.polarization
    }
    fn depth(&self, i: usize, electric: bool) -> f64 {
        let (cells, layer) = (self.config.cells as f64, self.config.pml_cells as f64);
        let x = i as f64 + 0.5;
        let v = if electric {
            if i < self.config.pml_cells {
                (layer - i as f64 - 0.5) / layer
            } else if i + self.config.pml_cells >= self.config.cells {
                (i as f64 - (cells - layer) + 0.5) / layer
            } else {
                0.0
            }
        } else if x < layer {
            (layer - x) / layer
        } else if x > cells - layer {
            (x - (cells - layer)) / layer
        } else {
            0.0
        };
        v.clamp(0.0, 1.0)
    }
    fn initialize_profiles(&mut self) {
        self.pml_sigma.fill(0.0);
        self.cpml_e = Cpml::new(self.config.cells);
        self.cpml_h = Cpml::new(self.config.cells - 1);
        if !matches!(self.config.boundary, Boundary1D::Pml | Boundary1D::Cpml) {
            return;
        }
        let eta = (self.mu / self.epsilon[0]).sqrt();
        let sigma_max = -(self.config.pml_order + 1.0) * self.config.pml_target_reflection.ln()
            / (2.0 * eta * self.config.pml_cells as f64 * self.config.dx);
        if self.config.boundary == Boundary1D::Pml {
            for i in 0..self.config.cells {
                self.pml_sigma[i] = sigma_max * self.depth(i, true).powf(self.config.pml_order);
            }
            return;
        }
        for e in [true, false] {
            let count = if e {
                self.config.cells
            } else {
                self.config.cells - 1
            };
            for i in 0..count {
                let depth = self.depth(i, e);
                if depth <= 0.0 {
                    continue;
                }
                let grade = depth.powf(self.config.pml_order);
                let k = 1.0 + (self.config.cpml_kappa_max - 1.0) * grade;
                let factor = if e {
                    1.0
                } else {
                    self.mu / (0.5 * (self.epsilon[i] + self.epsilon[i + 1]))
                };
                let sigma = sigma_max * grade * factor;
                let alpha = self.config.cpml_alpha_fraction * sigma_max * (1.0 - depth) * factor;
                let constitutive = if e { self.epsilon[i] } else { self.mu };
                let b = (-(sigma / k + alpha) * self.dt / constitutive).exp();
                let denom = sigma * k + k * k * alpha;
                let c = if denom > 0.0 {
                    sigma * (b - 1.0) / denom
                } else {
                    0.0
                };
                let p = if e {
                    &mut self.cpml_e
                } else {
                    &mut self.cpml_h
                };
                p.k[i] = k;
                p.b[i] = b;
                p.c[i] = c;
            }
        }
    }
    fn preprocess(&mut self) {
        for i in 0..self.config.cells {
            let eps = self.epsilon[i] + self.lumped_eps[i];
            let sigma = self.conductivity[i] + self.lumped_sigma[i] + self.pml_sigma[i];
            let loss = sigma * self.dt / (2.0 * eps);
            self.e_decay[i] = (1.0 - loss) / (1.0 + loss);
            self.e_curl[i] = self.dt / (eps * self.config.dx) / (1.0 + loss);
        }
        for i in 0..self.hy.len() {
            let eps = 0.5 * (self.epsilon[i] + self.epsilon[i + 1]);
            let sigma = 0.5 * (self.pml_sigma[i] + self.pml_sigma[i + 1]) * self.mu / eps;
            let loss = sigma * self.dt / (2.0 * self.mu);
            self.h_decay[i] = (1.0 - loss) / (1.0 + loss);
            self.h_curl[i] = self.dt / (self.mu * self.config.dx) / (1.0 + loss);
        }
    }
    fn range(&self, begin: usize, end: usize, eps: f64, sigma: f64) -> Result<()> {
        require(
            begin < end
                && end <= self.config.cells
                && eps.is_finite()
                && eps > 0.0
                && sigma.is_finite()
                && sigma >= 0.0,
            "invalid material range/parameters",
        )
    }
    fn clear_dispersion(&mut self, begin: usize, end: usize) {
        for i in begin..end {
            self.dispersion[i] = DispersionModel1D::None;
            self.polarization[i] = 0.0;
            self.previous_p[i] = 0.0;
            self.ade_a[i] = 0.0;
            self.ade_b[i] = 0.0;
            self.ade_c[i] = 0.0;
        }
    }
    pub fn set_material(
        &mut self,
        begin: usize,
        end: usize,
        epsilon_r: f64,
        sigma: f64,
    ) -> Result<()> {
        self.range(begin, end, epsilon_r, sigma)?;
        self.clear_dispersion(begin, end);
        for i in begin..end {
            self.epsilon[i] = EPSILON0 * epsilon_r;
            self.conductivity[i] = sigma;
        }
        self.initialize_profiles();
        self.preprocess();
        Ok(())
    }
    fn dispersive_material(
        &mut self,
        range: std::ops::Range<usize>,
        epsilon_r: f64,
        sigma: f64,
        model: DispersionModel1D,
        a: f64,
        b: f64,
        c: f64,
    ) -> Result<()> {
        self.range(range.start, range.end, epsilon_r, sigma)?;
        for i in range.clone() {
            require(
                self.lumped_eps[i] == 0.0
                    && self.lumped_sigma[i] == 0.0
                    && self.inductor_drive[i] == 0.0,
                "dispersive and lumped models cannot overlap",
            )?;
        }
        self.clear_dispersion(range.start, range.end);
        for i in range {
            self.epsilon[i] = EPSILON0 * epsilon_r;
            self.conductivity[i] = sigma;
            self.dispersion[i] = model;
            self.ade_a[i] = a;
            self.ade_b[i] = b;
            self.ade_c[i] = c;
        }
        self.initialize_profiles();
        self.preprocess();
        Ok(())
    }
    pub fn set_debye_material(
        &mut self,
        range: std::ops::Range<usize>,
        epsilon_infinity_r: f64,
        delta_epsilon_r: f64,
        relaxation_time: f64,
        sigma: f64,
    ) -> Result<()> {
        require(
            delta_epsilon_r.is_finite() && delta_epsilon_r >= 0.0,
            "invalid Debye increment",
        )?;
        positive(relaxation_time, "Debye relaxation time")?;
        let q = self.dt / (2.0 * relaxation_time);
        self.dispersive_material(
            range,
            epsilon_infinity_r,
            sigma,
            DispersionModel1D::Debye,
            (1.0 - q) / (1.0 + q),
            EPSILON0 * delta_epsilon_r * q / (1.0 + q),
            0.0,
        )
    }
    pub fn set_drude_material(
        &mut self,
        range: std::ops::Range<usize>,
        epsilon_infinity_r: f64,
        plasma_frequency: f64,
        collision_frequency: f64,
        sigma: f64,
    ) -> Result<()> {
        require(
            plasma_frequency.is_finite()
                && plasma_frequency >= 0.0
                && collision_frequency.is_finite()
                && collision_frequency >= 0.0,
            "invalid Drude frequencies",
        )?;
        let d = 1.0 + 0.5 * collision_frequency * self.dt;
        self.dispersive_material(
            range,
            epsilon_infinity_r,
            sigma,
            DispersionModel1D::Drude,
            2.0 / d,
            -(1.0 - 0.5 * collision_frequency * self.dt) / d,
            EPSILON0 * plasma_frequency.powi(2) * self.dt.powi(2) / d,
        )
    }
    pub fn set_lorentz_material(
        &mut self,
        range: std::ops::Range<usize>,
        epsilon_infinity_r: f64,
        delta_epsilon_r: f64,
        resonance_frequency: f64,
        damping: f64,
        sigma: f64,
    ) -> Result<()> {
        require(
            delta_epsilon_r.is_finite()
                && delta_epsilon_r >= 0.0
                && damping.is_finite()
                && damping >= 0.0,
            "invalid Lorentz parameters",
        )?;
        positive(resonance_frequency, "Lorentz resonance")?;
        let od = resonance_frequency * self.dt;
        let d = 1.0 + 0.5 * damping * self.dt;
        self.dispersive_material(
            range,
            epsilon_infinity_r,
            sigma,
            DispersionModel1D::Lorentz,
            (2.0 - od * od) / d,
            -(1.0 - 0.5 * damping * self.dt) / d,
            EPSILON0 * delta_epsilon_r * od * od / d,
        )
    }
    pub fn dispersion_model(&self, cell: usize) -> Result<DispersionModel1D> {
        require(cell < self.config.cells, "cell out of range")?;
        Ok(self.dispersion[cell])
    }
    pub fn set_tfsf_source(
        &mut self,
        cell: usize,
        incident: impl Fn(f64) -> f64 + 'a,
        impedance: f64,
    ) -> Result<()> {
        require(
            cell >= 1 && cell < self.config.cells - 2,
            "invalid TFSF interface",
        )?;
        if matches!(self.config.boundary, Boundary1D::Pml | Boundary1D::Cpml) {
            require(
                cell > self.config.pml_cells && cell + self.config.pml_cells < self.config.cells,
                "TFSF interface in absorbing layer",
            )?;
        }
        require(
            self.dispersion[cell] == DispersionModel1D::None
                && self.lumped_eps[cell] == 0.0
                && self.lumped_sigma[cell] == 0.0
                && self.inductor_drive[cell] == 0.0,
            "TFSF interface in dispersive/lumped material",
        )?;
        let impedance = if impedance > 0.0 {
            impedance
        } else {
            (self.mu / self.epsilon[cell]).sqrt()
        };
        positive(impedance, "TFSF impedance")?;
        self.tfsf = Some(Tfsf {
            cell,
            incident: Box::new(incident),
            impedance,
        });
        Ok(())
    }
    pub fn clear_tfsf_source(&mut self) {
        self.tfsf = None;
    }
    pub fn set_parallel_lumped_rlc(
        &mut self,
        cell: usize,
        length: f64,
        area: f64,
        resistance: f64,
        inductance: f64,
        capacitance: f64,
    ) -> Result<()> {
        require(cell > 0 && cell < self.config.cells - 1, "invalid RLC cell")?;
        positive(length, "RLC length")?;
        positive(area, "RLC area")?;
        require(
            resistance > 0.0 && inductance > 0.0 && capacitance.is_finite() && capacitance >= 0.0,
            "invalid RLC values",
        )?;
        require(
            self.dispersion[cell] == DispersionModel1D::None,
            "RLC and dispersive models cannot overlap",
        )?;
        self.lumped_eps[cell] = capacitance * length / area;
        self.lumped_sigma[cell] = if resistance.is_finite() {
            length / (resistance * area)
        } else {
            0.0
        };
        self.inductor_drive[cell] = if inductance.is_finite() {
            length / (inductance * area)
        } else {
            0.0
        };
        self.inductor_current[cell] = 0.0;
        self.preprocess();
        Ok(())
    }
    pub fn clear_lumped_rlc(&mut self, cell: usize) -> Result<()> {
        require(cell < self.config.cells, "RLC cell out of range")?;
        self.lumped_eps[cell] = 0.0;
        self.lumped_sigma[cell] = 0.0;
        self.inductor_drive[cell] = 0.0;
        self.inductor_current[cell] = 0.0;
        self.preprocess();
        Ok(())
    }
    pub fn lumped_inductor_current_density(&self, cell: usize) -> Result<f64> {
        require(cell < self.config.cells, "cell out of range")?;
        Ok(self.inductor_current[cell])
    }
    pub fn initialize_gaussian(&mut self, center_fraction: f64, width_fraction: f64) -> Result<()> {
        positive(width_fraction, "Gaussian width")?;
        require(center_fraction.is_finite(), "Gaussian center")?;
        let center = center_fraction * (self.config.cells - 1) as f64;
        let width = width_fraction * self.config.cells as f64;
        for (i, e) in self.ez.iter_mut().enumerate() {
            let z = (i as f64 - center) / width;
            *e = (-0.5 * z * z).exp();
        }
        if self.config.boundary == Boundary1D::Pec {
            self.ez[0] = 0.0;
            self.ez[self.config.cells - 1] = 0.0;
        }
        Ok(())
    }
    pub fn add_soft_source(&mut self, cell: usize, value: f64) -> Result<()> {
        require(
            cell < self.config.cells && value.is_finite(),
            "invalid soft source",
        )?;
        self.ez[cell] += value;
        Ok(())
    }
    pub fn set_hard_source(&mut self, cell: usize, value: f64) -> Result<()> {
        require(
            cell < self.config.cells && value.is_finite(),
            "invalid hard source",
        )?;
        self.ez[cell] = value;
        Ok(())
    }
    pub fn local_wave_speed(&self, cell: usize) -> Result<f64> {
        require(cell < self.config.cells, "cell out of range")?;
        Ok(1.0 / (self.epsilon[cell] * self.mu).sqrt())
    }
    pub fn step(&mut self, count: usize) -> Result<()> {
        let n = self.config.cells;
        for _ in 0..count {
            let old = [self.ez[0], self.ez[1], self.ez[n - 1], self.ez[n - 2]];
            for i in 0..self.hy.len() {
                let mut derivative = (self.ez[i] - self.ez[i + 1]) / self.config.dx;
                if self.config.boundary == Boundary1D::Cpml {
                    derivative = self.cpml_h.apply(i, derivative);
                }
                self.hy[i] =
                    self.h_decay[i] * self.hy[i] + self.h_curl[i] * self.config.dx * derivative;
            }
            if let Some(t) = &self.tfsf {
                let inc = (t.incident)(self.time);
                require(inc.is_finite(), "nonfinite incident field")?;
                self.hy[t.cell - 1] -= self.dt / (self.mu * self.config.dx) * inc;
            }
            for i in 1..n - 1 {
                let e = self.ez[i];
                self.inductor_current[i] += self.dt * self.inductor_drive[i] * e;
                let mut curl = (self.hy[i - 1] - self.hy[i]) / self.config.dx;
                if self.config.boundary == Boundary1D::Cpml {
                    curl = self.cpml_e.apply(i, curl);
                }
                match self.dispersion[i] {
                    DispersionModel1D::None => {
                        self.ez[i] = self.e_decay[i] * e
                            + self.e_curl[i] * self.config.dx * (curl - self.inductor_current[i]);
                    }
                    DispersionModel1D::Debye => {
                        let p = self.polarization[i];
                        let (a, b) = (self.ade_a[i], self.ade_b[i]);
                        let lhs =
                            self.epsilon[i] / self.dt + 0.5 * self.conductivity[i] + b / self.dt;
                        let rhs = (self.epsilon[i] / self.dt - 0.5 * self.conductivity[i]) * e
                            + curl
                            - ((a - 1.0) * p + b * e) / self.dt;
                        let new_e = rhs / lhs;
                        self.previous_p[i] = p;
                        self.polarization[i] = a * p + b * (new_e + e);
                        self.ez[i] = new_e;
                    }
                    DispersionModel1D::Drude | DispersionModel1D::Lorentz => {
                        let p = self.polarization[i];
                        let new_p = self.ade_a[i] * p
                            + self.ade_b[i] * self.previous_p[i]
                            + self.ade_c[i] * e;
                        self.previous_p[i] = p;
                        self.polarization[i] = new_p;
                        self.ez[i] = ((self.epsilon[i] - 0.5 * self.conductivity[i] * self.dt) * e
                            + self.dt * curl
                            - (new_p - p))
                            / (self.epsilon[i] + 0.5 * self.conductivity[i] * self.dt);
                    }
                }
            }
            if let Some(t) = &self.tfsf {
                let delay = 0.5 * self.config.dx * (self.epsilon[t.cell] * self.mu).sqrt();
                let inc = (t.incident)(self.time + 0.5 * self.dt + delay);
                require(inc.is_finite(), "nonfinite incident field")?;
                self.ez[t.cell] += self.e_curl[t.cell] * (-inc / t.impedance);
            }
            if self.config.boundary == Boundary1D::Mur1 {
                let cl = 1.0 / (self.epsilon[0] * self.mu).sqrt();
                let cr = 1.0 / (self.epsilon[n - 1] * self.mu).sqrt();
                let kl = (cl * self.dt - self.config.dx) / (cl * self.dt + self.config.dx);
                let kr = (cr * self.dt - self.config.dx) / (cr * self.dt + self.config.dx);
                self.ez[0] = old[1] + kl * (self.ez[1] - old[0]);
                self.ez[n - 1] = old[3] + kr * (self.ez[n - 2] - old[2]);
            } else {
                self.ez[0] = 0.0;
                self.ez[n - 1] = 0.0;
            }
            self.time += self.dt;
            finite_slice(&self.ez, "electric field")?;
            finite_slice(&self.hy, "magnetic field")?;
        }
        Ok(())
    }
    pub fn energy(&self) -> f64 {
        self.ez
            .iter()
            .enumerate()
            .map(|(i, e)| 0.5 * (self.epsilon[i] + self.lumped_eps[i]) * e * e * self.config.dx)
            .sum::<f64>()
            + self
                .hy
                .iter()
                .map(|h| 0.5 * self.mu * h * h * self.config.dx)
                .sum::<f64>()
    }
}
