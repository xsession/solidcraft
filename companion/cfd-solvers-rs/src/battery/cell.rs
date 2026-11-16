//! Particles, OCV curves, SPM/SPMe, thermal slab, drive cycles and EIS.
use crate::{require, Complex64 as C, Result, SolverError};
use std::{f64::consts::PI, sync::Arc};

pub const FARADAY_CONSTANT: f64 = 96485.33212;
pub const GAS_CONSTANT: f64 = 8.31446261815324;

/// Open-circuit voltage as a function of surface stoichiometry.
pub type OcvFunction = Arc<dyn Fn(f64) -> f64 + Send + Sync>;

pub(crate) fn clamp01(x: f64) -> f64 {
    x.clamp(1e-6, 1.0 - 1e-6)
}

/// Thomas algorithm; `rhs` is overwritten with the solution.
pub(crate) fn solve_tridiagonal(
    lower: &[f64],
    diagonal: &[f64],
    upper: &[f64],
    rhs: &mut [f64],
) -> Result<()> {
    let n = rhs.len();
    require(
        n > 0 && diagonal.len() == n && lower.len() + 1 == n && upper.len() + 1 == n,
        "tridiagonal system size mismatch",
    )?;
    let mut cprime = vec![0.0; n - 1];
    let mut pivot = diagonal[0];
    require(pivot.abs() >= 1e-30, "singular tridiagonal system")?;
    if n > 1 {
        cprime[0] = upper[0] / pivot;
    }
    rhs[0] /= pivot;
    for i in 1..n {
        pivot = diagonal[i] - lower[i - 1] * cprime[i - 1];
        require(pivot.abs() >= 1e-30, "singular tridiagonal system")?;
        if i + 1 < n {
            cprime[i] = upper[i] / pivot;
        }
        rhs[i] = (rhs[i] - lower[i - 1] * rhs[i - 1]) / pivot;
    }
    for i in (1..n).rev() {
        rhs[i - 1] -= cprime[i - 1] * rhs[i];
    }
    Ok(())
}

pub(crate) fn effective_transport(bulk: f64, porosity: f64, exponent: f64) -> f64 {
    bulk * porosity.clamp(1e-6, 1.0).powf(exponent)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleConfig {
    pub radius_m: f64,
    pub diffusivity_m2_per_s: f64,
    pub maximum_concentration_mol_per_m3: f64,
    pub initial_stoichiometry: f64,
    pub radial_cells: usize,
}

impl Default for ParticleConfig {
    fn default() -> Self {
        Self {
            radius_m: 5e-6,
            diffusivity_m2_per_s: 1e-14,
            maximum_concentration_mol_per_m3: 3e4,
            initial_stoichiometry: 0.5,
            radial_cells: 24,
        }
    }
}

/// Finite-volume spherical diffusion with implicit Euler and a surface flux boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct SphericalDiffusionParticle {
    config: ParticleConfig,
    concentration: Vec<f64>,
}

impl SphericalDiffusionParticle {
    pub fn new(config: ParticleConfig) -> Result<Self> {
        require(
            config.radius_m > 0.0
                && config.diffusivity_m2_per_s > 0.0
                && config.maximum_concentration_mol_per_m3 > 0.0
                && config.radial_cells >= 3
                && config.radial_cells <= 1 << 20,
            "invalid spherical-particle configuration",
        )?;
        let mut p = Self {
            config,
            concentration: Vec::new(),
        };
        p.reset(config.initial_stoichiometry);
        Ok(p)
    }
    pub fn reset(&mut self, stoichiometry: f64) {
        self.concentration = vec![
            clamp01(stoichiometry)
                * self.config.maximum_concentration_mol_per_m3;
            self.config.radial_cells
        ];
    }
    /// `outward_flux` in mol/(m^2 s) leaving the particle surface.
    pub fn step(&mut self, outward_flux: f64, dt: f64) -> Result<()> {
        require(
            dt > 0.0 && outward_flux.is_finite(),
            "invalid particle step",
        )?;
        let n = self.concentration.len();
        let dr = self.config.radius_m / n as f64;
        let d = self.config.diffusivity_m2_per_s;
        let (mut lower, mut diag, mut upper) = (vec![0.0; n - 1], vec![1.0; n], vec![0.0; n - 1]);
        let mut rhs = self.concentration.clone();
        for i in 0..n - 1 {
            let rf = (i + 1) as f64 * dr;
            let area = 4.0 * PI * rf * rf;
            let (rlo_i, rhi_i) = (i as f64 * dr, (i + 1) as f64 * dr);
            let vi = (4.0 / 3.0) * PI * (rhi_i * rhi_i * rhi_i - rlo_i * rlo_i * rlo_i);
            let (rlo_j, rhi_j) = (rhi_i, (i + 2) as f64 * dr);
            let vj = (4.0 / 3.0) * PI * (rhi_j * rhi_j * rhi_j - rlo_j * rlo_j * rlo_j);
            let (ki, kj) = (d * area / (dr * vi), d * area / (dr * vj));
            diag[i] += dt * ki;
            upper[i] -= dt * ki;
            lower[i] -= dt * kj;
            diag[i + 1] += dt * kj;
        }
        let r = self.config.radius_m;
        let surface_area = 4.0 * PI * r * r;
        let rlo = r - dr;
        let v_last = (4.0 / 3.0) * PI * (r * r * r - rlo * rlo * rlo);
        rhs[n - 1] -= dt * surface_area * outward_flux / v_last;
        solve_tridiagonal(&lower, &diag, &upper, &mut rhs)?;
        let m = self.config.maximum_concentration_mol_per_m3;
        rhs.iter_mut()
            .for_each(|c| *c = c.clamp(1e-12 * m, (1.0 - 1e-12) * m));
        self.concentration = rhs;
        Ok(())
    }
    pub fn average_concentration(&self) -> f64 {
        let n = self.concentration.len();
        let dr = self.config.radius_m / n as f64;
        let (mut amount, mut volume) = (0.0, 0.0);
        for (i, c) in self.concentration.iter().enumerate() {
            let (r0, r1) = (i as f64 * dr, (i + 1) as f64 * dr);
            let v = r1 * r1 * r1 - r0 * r0 * r0;
            amount += c * v;
            volume += v;
        }
        amount / volume
    }
    /// Half-cell extrapolation of the outer cell value with the given surface flux.
    pub fn surface_concentration(&self, outward_flux: f64) -> f64 {
        let dr = self.config.radius_m / self.concentration.len() as f64;
        let m = self.config.maximum_concentration_mol_per_m3;
        let last = self.concentration.last().copied().unwrap_or(0.0);
        (last - outward_flux * 0.5 * dr / self.config.diffusivity_m2_per_s)
            .clamp(1e-12 * m, (1.0 - 1e-12) * m)
    }
    pub fn average_stoichiometry(&self) -> f64 {
        self.average_concentration() / self.config.maximum_concentration_mol_per_m3
    }
    pub fn surface_stoichiometry(&self, outward_flux: f64) -> f64 {
        self.surface_concentration(outward_flux) / self.config.maximum_concentration_mol_per_m3
    }
    pub fn concentrations(&self) -> &[f64] {
        &self.concentration
    }
    pub fn config(&self) -> &ParticleConfig {
        &self.config
    }
}

/// Graphite open-circuit potential fit.
pub fn graphite_ocv_v(x: f64) -> f64 {
    let x = clamp01(x);
    0.08 + 0.75 * (-12.0 * x).exp()
        + 0.12 * (-60.0 * (x - 0.12) * (x - 0.12)).exp()
        + 0.04 * (1.0 - x)
}

/// NMC open-circuit potential fit.
pub fn nmc_ocv_v(x: f64) -> f64 {
    let x = clamp01(x);
    4.45 - 1.15 * x - 0.10 * ((x - 0.55) / 0.08).tanh() + 0.03 * (6.0 * x).sin()
}

/// `2 R T / F asinh(j / (2 i0))`.
pub fn symmetric_butler_volmer_overpotential_v(j: f64, i0: f64, temperature: f64) -> Result<f64> {
    require(
        i0 > 0.0 && temperature > 0.0,
        "invalid Butler-Volmer parameters",
    )?;
    Ok(2.0 * GAS_CONSTANT * temperature / FARADAY_CONSTANT * (j / (2.0 * i0)).asinh())
}

#[derive(Clone)]
pub struct ElectrodeConfig {
    pub particle: ParticleConfig,
    pub thickness_m: f64,
    pub active_volume_fraction: f64,
    pub solid_conductivity_s_per_m: f64,
    pub exchange_current_density_a_per_m2: f64,
    /// `None` uses graphite (negative) or NMC (positive).
    pub open_circuit_voltage_v: Option<OcvFunction>,
}

impl Default for ElectrodeConfig {
    fn default() -> Self {
        Self {
            particle: ParticleConfig::default(),
            thickness_m: 80e-6,
            active_volume_fraction: 0.60,
            solid_conductivity_s_per_m: 100.0,
            exchange_current_density_a_per_m2: 3.0,
            open_circuit_voltage_v: None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrolyteConfig {
    pub negative_thickness_m: f64,
    pub separator_thickness_m: f64,
    pub positive_thickness_m: f64,
    pub initial_concentration_mol_per_m3: f64,
    pub diffusivity_m2_per_s: f64,
    pub conductivity_s_per_m: f64,
    pub transference_number: f64,
    pub negative_porosity: f64,
    pub separator_porosity: f64,
    pub positive_porosity: f64,
    pub bruggeman_exponent: f64,
    pub cells: usize,
}

impl Default for ElectrolyteConfig {
    fn default() -> Self {
        Self {
            negative_thickness_m: 80e-6,
            separator_thickness_m: 25e-6,
            positive_thickness_m: 75e-6,
            initial_concentration_mol_per_m3: 1000.0,
            diffusivity_m2_per_s: 2e-10,
            conductivity_s_per_m: 1.0,
            transference_number: 0.38,
            negative_porosity: 0.30,
            separator_porosity: 0.50,
            positive_porosity: 0.30,
            bruggeman_exponent: 1.5,
            cells: 45,
        }
    }
}

pub(crate) fn region_porosity(e: &ElectrolyteConfig, x: f64) -> f64 {
    if x < e.negative_thickness_m {
        e.negative_porosity
    } else if x < e.negative_thickness_m + e.separator_thickness_m {
        e.separator_porosity
    } else {
        e.positive_porosity
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedThermalConfig {
    pub mass_kg: f64,
    pub heat_capacity_j_per_kg_k: f64,
    pub heat_transfer_coefficient_w_per_m2_k: f64,
    pub cooling_area_m2: f64,
    pub ambient_temperature_k: f64,
}

impl Default for LumpedThermalConfig {
    fn default() -> Self {
        Self {
            mass_kg: 0.045,
            heat_capacity_j_per_kg_k: 900.0,
            heat_transfer_coefficient_w_per_m2_k: 8.0,
            cooling_area_m2: 0.012,
            ambient_temperature_k: 298.15,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalSlabConfig {
    pub thickness_m: f64,
    pub cells: usize,
    pub density_kg_per_m3: f64,
    pub heat_capacity_j_per_kg_k: f64,
    pub conductivity_w_per_m_k: f64,
    pub convection_w_per_m2_k: f64,
    pub ambient_temperature_k: f64,
}

impl Default for ThermalSlabConfig {
    fn default() -> Self {
        Self {
            thickness_m: 180e-6,
            cells: 32,
            density_kg_per_m3: 2200.0,
            heat_capacity_j_per_kg_k: 1000.0,
            conductivity_w_per_m_k: 1.0,
            convection_w_per_m2_k: 8.0,
            ambient_temperature_k: 298.15,
        }
    }
}

/// Implicit 1-D conduction through the cell stack with convection on both faces.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalSlab1D {
    config: ThermalSlabConfig,
    temperature: Vec<f64>,
}

impl ThermalSlab1D {
    pub fn new(config: ThermalSlabConfig) -> Result<Self> {
        require(
            config.thickness_m > 0.0
                && config.cells >= 3
                && config.cells <= 1 << 22
                && config.density_kg_per_m3 > 0.0
                && config.heat_capacity_j_per_kg_k > 0.0
                && config.conductivity_w_per_m_k >= 0.0,
            "invalid thermal-slab configuration",
        )?;
        Ok(Self {
            config,
            temperature: vec![config.ambient_temperature_k; config.cells],
        })
    }
    pub fn reset(&mut self, temperature_k: f64) {
        self.temperature = vec![temperature_k; self.config.cells];
    }
    pub fn step(&mut self, dt: f64, volumetric_heat_w_per_m3: &[f64]) -> Result<()> {
        let n = self.temperature.len();
        require(
            dt > 0.0 && volumetric_heat_w_per_m3.len() == n,
            "invalid thermal-slab step",
        )?;
        let c = self.config;
        let dx = c.thickness_m / n as f64;
        let rho_cp = c.density_kg_per_m3 * c.heat_capacity_j_per_kg_k;
        let a = c.conductivity_w_per_m_k / (rho_cp * dx * dx);
        let (mut lower, mut diag, mut upper) = (
            vec![-dt * a; n - 1],
            vec![1.0 + 2.0 * dt * a; n],
            vec![-dt * a; n - 1],
        );
        let mut rhs: Vec<f64> = (0..n)
            .map(|i| self.temperature[i] + dt * volumetric_heat_w_per_m3[i] / rho_cp)
            .collect();
        let beta = c.convection_w_per_m2_k / (rho_cp * dx);
        diag[0] = 1.0 + dt * (a + beta);
        upper[0] = -dt * a;
        rhs[0] += dt * beta * c.ambient_temperature_k;
        diag[n - 1] = 1.0 + dt * (a + beta);
        lower[n - 2] = -dt * a;
        rhs[n - 1] += dt * beta * c.ambient_temperature_k;
        solve_tridiagonal(&lower, &diag, &upper, &mut rhs)?;
        self.temperature = rhs;
        Ok(())
    }
    pub fn temperature(&self) -> &[f64] {
        &self.temperature
    }
    pub fn average_temperature(&self) -> f64 {
        self.temperature.iter().sum::<f64>() / self.temperature.len() as f64
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DegradationConfig {
    pub enable_sei: bool,
    pub enable_plating: bool,
    pub enable_active_material_loss: bool,
    pub enable_particle_cracking: bool,
    pub initial_sei_thickness_m: f64,
    pub sei_rate_m2_per_s: f64,
    pub plating_efficiency: f64,
    pub active_material_loss_per_coulomb: f64,
    pub cracking_rate_per_gradient_second: f64,
}

impl Default for DegradationConfig {
    fn default() -> Self {
        Self {
            enable_sei: true,
            enable_plating: true,
            enable_active_material_loss: true,
            enable_particle_cracking: true,
            initial_sei_thickness_m: 5e-9,
            sei_rate_m2_per_s: 2e-20,
            plating_efficiency: 0.02,
            active_material_loss_per_coulomb: 2e-10,
            cracking_rate_per_gradient_second: 5e-10,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DegradationState {
    pub sei_thickness_m: f64,
    pub plated_lithium_mol: f64,
    pub active_material_fraction: f64,
    pub crack_damage: f64,
}

impl Default for DegradationState {
    fn default() -> Self {
        Self {
            sei_thickness_m: 0.0,
            plated_lithium_mol: 0.0,
            active_material_fraction: 1.0,
            crack_damage: 0.0,
        }
    }
}

#[derive(Clone)]
pub struct LithiumIonCellConfig {
    pub area_m2: f64,
    pub negative: ElectrodeConfig,
    pub positive: ElectrodeConfig,
    pub electrolyte: ElectrolyteConfig,
    pub thermal: LumpedThermalConfig,
    pub degradation: DegradationConfig,
    pub contact_resistance_ohm: f64,
    pub initial_temperature_k: f64,
}

impl Default for LithiumIonCellConfig {
    fn default() -> Self {
        Self {
            area_m2: 0.010,
            negative: ElectrodeConfig::default(),
            positive: ElectrodeConfig::default(),
            electrolyte: ElectrolyteConfig::default(),
            thermal: LumpedThermalConfig::default(),
            degradation: DegradationConfig::default(),
            contact_resistance_ohm: 0.010,
            initial_temperature_k: 298.15,
        }
    }
}

/// Graphite / NMC parameter set.
pub fn default_graphite_nmc_config() -> LithiumIonCellConfig {
    let mut c = LithiumIonCellConfig::default();
    c.negative.particle = ParticleConfig {
        radius_m: 5.5e-6,
        diffusivity_m2_per_s: 3e-14,
        maximum_concentration_mol_per_m3: 3.1e4,
        initial_stoichiometry: 0.82,
        radial_cells: 24,
    };
    c.negative.thickness_m = 80e-6;
    c.negative.active_volume_fraction = 0.62;
    c.negative.solid_conductivity_s_per_m = 100.0;
    c.negative.exchange_current_density_a_per_m2 = 4.0;
    c.negative.open_circuit_voltage_v = Some(Arc::new(graphite_ocv_v));
    c.positive.particle = ParticleConfig {
        radius_m: 5e-6,
        diffusivity_m2_per_s: 1e-13,
        maximum_concentration_mol_per_m3: 5.1e4,
        initial_stoichiometry: 0.42,
        radial_cells: 24,
    };
    c.positive.thickness_m = 75e-6;
    c.positive.active_volume_fraction = 0.58;
    c.positive.solid_conductivity_s_per_m = 10.0;
    c.positive.exchange_current_density_a_per_m2 = 4.0;
    c.positive.open_circuit_voltage_v = Some(Arc::new(nmc_ocv_v));
    c.electrolyte.negative_thickness_m = c.negative.thickness_m;
    c.electrolyte.positive_thickness_m = c.positive.thickness_m;
    c.contact_resistance_ohm = 0.012;
    c
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BatteryStepResult {
    pub time_s: f64,
    /// Positive = discharge.
    pub current_a: f64,
    pub terminal_voltage_v: f64,
    pub open_circuit_voltage_v: f64,
    pub reaction_overpotential_v: f64,
    pub electrolyte_overpotential_v: f64,
    pub ohmic_drop_v: f64,
    pub negative_surface_stoichiometry: f64,
    pub positive_surface_stoichiometry: f64,
    pub electrolyte_min_concentration_mol_per_m3: f64,
    pub electrolyte_max_concentration_mol_per_m3: f64,
    pub temperature_k: f64,
    pub heat_generation_w: f64,
}

/// Models that a drive cycle or pack can step.
pub trait BatteryModel {
    fn step(&mut self, current_a: f64, dt_s: f64) -> Result<BatteryStepResult>;
}

/// Single-particle model with degradation and lumped thermal.
#[derive(Clone)]
pub struct SingleParticleModel {
    config: LithiumIonCellConfig,
    negative_ocv: OcvFunction,
    positive_ocv: OcvFunction,
    negative: SphericalDiffusionParticle,
    positive: SphericalDiffusionParticle,
    degradation: DegradationState,
    temperature_k: f64,
    time_s: f64,
    external_thermal: bool,
}

impl SingleParticleModel {
    pub fn new(config: LithiumIonCellConfig) -> Result<Self> {
        require(config.area_m2 > 0.0, "invalid battery area")?;
        let negative = SphericalDiffusionParticle::new(config.negative.particle)?;
        let positive = SphericalDiffusionParticle::new(config.positive.particle)?;
        let negative_ocv = config
            .negative
            .open_circuit_voltage_v
            .clone()
            .unwrap_or_else(|| Arc::new(graphite_ocv_v));
        let positive_ocv = config
            .positive
            .open_circuit_voltage_v
            .clone()
            .unwrap_or_else(|| Arc::new(nmc_ocv_v));
        let mut m = Self {
            config,
            negative_ocv,
            positive_ocv,
            negative,
            positive,
            degradation: DegradationState::default(),
            temperature_k: 0.0,
            time_s: 0.0,
            external_thermal: false,
        };
        m.reset();
        Ok(m)
    }
    pub fn reset(&mut self) {
        self.external_thermal = false;
        self.negative
            .reset(self.config.negative.particle.initial_stoichiometry);
        self.positive
            .reset(self.config.positive.particle.initial_stoichiometry);
        self.degradation = DegradationState {
            sei_thickness_m: self.config.degradation.initial_sei_thickness_m,
            ..Default::default()
        };
        self.temperature_k = self.config.initial_temperature_k;
        self.time_s = 0.0;
    }
    pub fn set_external_temperature(&mut self, t: f64) -> Result<()> {
        require(
            t.is_finite() && t > 0.0,
            "invalid external cell temperature",
        )?;
        self.temperature_k = t;
        self.external_thermal = true;
        Ok(())
    }
    fn negative_flux(&self, current: f64) -> f64 {
        let c = &self.config;
        let a = 3.0 * c.negative.active_volume_fraction / c.negative.particle.radius_m;
        current / (FARADAY_CONSTANT * c.area_m2 * c.negative.thickness_m * a)
    }
    fn positive_flux(&self, current: f64) -> f64 {
        let c = &self.config;
        let a = 3.0 * c.positive.active_volume_fraction / c.positive.particle.radius_m;
        -current / (FARADAY_CONSTANT * c.area_m2 * c.positive.thickness_m * a)
    }
    fn evaluate(
        &self,
        current: f64,
        electrolyte_eta: f64,
        cmin: f64,
        cmax: f64,
    ) -> Result<BatteryStepResult> {
        let c = &self.config;
        let xn = clamp01(
            self.negative
                .surface_stoichiometry(self.negative_flux(current)),
        );
        let xp = clamp01(
            self.positive
                .surface_stoichiometry(self.positive_flux(current)),
        );
        let ocv = (self.positive_ocv)(xp) - (self.negative_ocv)(xn);
        let an = 3.0 * c.negative.active_volume_fraction / c.negative.particle.radius_m;
        let ap = 3.0 * c.positive.active_volume_fraction / c.positive.particle.radius_m;
        let jn = current / (c.area_m2 * c.negative.thickness_m * an);
        let jp = current / (c.area_m2 * c.positive.thickness_m * ap);
        let i0n = c.negative.exchange_current_density_a_per_m2 * 2.0 * (xn * (1.0 - xn)).sqrt();
        let i0p = c.positive.exchange_current_density_a_per_m2 * 2.0 * (xp * (1.0 - xp)).sqrt();
        let reaction =
            symmetric_butler_volmer_overpotential_v(jn, i0n.max(1e-12), self.temperature_k)?.abs()
                + symmetric_butler_volmer_overpotential_v(jp, i0p.max(1e-12), self.temperature_k)?
                    .abs();
        let solid_r = c.negative.thickness_m / (c.negative.solid_conductivity_s_per_m * c.area_m2)
            + c.positive.thickness_m / (c.positive.solid_conductivity_s_per_m * c.area_m2);
        let ohmic = current * (c.contact_resistance_ohm + solid_r);
        let voltage =
            ocv - (reaction + electrolyte_eta.abs()).copysign(if current == 0.0 {
                1.0
            } else {
                current
            }) - ohmic;
        Ok(BatteryStepResult {
            time_s: self.time_s,
            current_a: current,
            terminal_voltage_v: voltage,
            open_circuit_voltage_v: ocv,
            reaction_overpotential_v: reaction,
            electrolyte_overpotential_v: electrolyte_eta,
            ohmic_drop_v: ohmic,
            negative_surface_stoichiometry: xn,
            positive_surface_stoichiometry: xp,
            electrolyte_min_concentration_mol_per_m3: cmin,
            electrolyte_max_concentration_mol_per_m3: cmax,
            temperature_k: self.temperature_k,
            heat_generation_w: current.abs() * (ocv - voltage).abs().max(0.0),
        })
    }
    fn advance_degradation(&mut self, current: f64, dt: f64, gradient: f64) {
        let d = self.config.degradation;
        let s = &mut self.degradation;
        if d.enable_sei {
            s.sei_thickness_m += dt * d.sei_rate_m2_per_s / s.sei_thickness_m.max(1e-10);
        }
        if d.enable_plating && current < 0.0 {
            s.plated_lithium_mol += d.plating_efficiency * (-current * dt) / FARADAY_CONSTANT;
        }
        if d.enable_active_material_loss {
            s.active_material_fraction = 0.5f64.max(
                s.active_material_fraction
                    - d.active_material_loss_per_coulomb * current.abs() * dt,
            );
        }
        if d.enable_particle_cracking {
            s.crack_damage = (s.crack_damage
                + d.cracking_rate_per_gradient_second * gradient.abs() * dt)
                .clamp(0.0, 1.0);
        }
    }
    fn advance_thermal(&mut self, current: f64, voltage: f64, ocv: f64, dt: f64) {
        if self.external_thermal {
            return;
        }
        let t = self.config.thermal;
        let q = current.abs() * (ocv - voltage).abs();
        let cooling = t.heat_transfer_coefficient_w_per_m2_k
            * t.cooling_area_m2
            * (self.temperature_k - t.ambient_temperature_k);
        self.temperature_k += dt * (q - cooling) / (t.mass_kg * t.heat_capacity_j_per_kg_k);
    }
    fn finish_step(
        &mut self,
        mut result: BatteryStepResult,
        current: f64,
        dt: f64,
        jn: f64,
    ) -> BatteryStepResult {
        let gradient =
            self.negative.surface_concentration(jn) - self.negative.average_concentration();
        self.advance_degradation(
            current,
            dt,
            gradient
                / self
                    .config
                    .negative
                    .particle
                    .maximum_concentration_mol_per_m3,
        );
        self.advance_thermal(
            current,
            result.terminal_voltage_v,
            result.open_circuit_voltage_v,
            dt,
        );
        result.temperature_k = self.temperature_k;
        result
    }
    /// Transactional: on error the model is unchanged.
    pub fn step(&mut self, current: f64, dt: f64) -> Result<BatteryStepResult> {
        require(dt > 0.0 && current.is_finite(), "invalid SPM step")?;
        let mut trial = self.clone();
        let (jn, jp) = (trial.negative_flux(current), trial.positive_flux(current));
        trial.negative.step(jn, dt)?;
        trial.positive.step(jp, dt)?;
        trial.time_s += dt;
        let c0 = trial.config.electrolyte.initial_concentration_mol_per_m3;
        let result = trial.evaluate(current, 0.0, c0, c0)?;
        let result = trial.finish_step(result, current, dt, jn);
        *self = trial;
        Ok(result)
    }
    pub fn state(&self, current: f64) -> Result<BatteryStepResult> {
        let c0 = self.config.electrolyte.initial_concentration_mol_per_m3;
        self.evaluate(current, 0.0, c0, c0)
    }
    pub fn negative_particle(&self) -> &SphericalDiffusionParticle {
        &self.negative
    }
    pub fn positive_particle(&self) -> &SphericalDiffusionParticle {
        &self.positive
    }
    pub fn degradation_state(&self) -> &DegradationState {
        &self.degradation
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub fn config(&self) -> &LithiumIonCellConfig {
        &self.config
    }
}

impl BatteryModel for SingleParticleModel {
    fn step(&mut self, current_a: f64, dt_s: f64) -> Result<BatteryStepResult> {
        SingleParticleModel::step(self, current_a, dt_s)
    }
}

/// SPM plus 1-D electrolyte diffusion and potential drop (SPMe).
#[derive(Clone)]
pub struct SingleParticleElectrolyteModel {
    base: SingleParticleModel,
    electrolyte_concentration: Vec<f64>,
    electrolyte_potential: Vec<f64>,
    solid_potential: Vec<f64>,
}

impl SingleParticleElectrolyteModel {
    pub fn new(config: LithiumIonCellConfig) -> Result<Self> {
        require(
            config.electrolyte.cells >= 3 && config.electrolyte.cells <= 1 << 22,
            "SPMe needs at least three electrolyte cells",
        )?;
        let mut m = Self {
            base: SingleParticleModel::new(config)?,
            electrolyte_concentration: Vec::new(),
            electrolyte_potential: Vec::new(),
            solid_potential: Vec::new(),
        };
        m.reset();
        Ok(m)
    }
    pub fn reset(&mut self) {
        self.base.reset();
        let e = self.base.config.electrolyte;
        self.electrolyte_concentration = vec![e.initial_concentration_mol_per_m3; e.cells];
        self.electrolyte_potential = vec![0.0; e.cells];
        self.solid_potential = vec![0.0; e.cells];
    }
    pub fn set_external_temperature(&mut self, t: f64) -> Result<()> {
        self.base.set_external_temperature(t)
    }
    fn length(&self) -> f64 {
        let e = &self.base.config.electrolyte;
        e.negative_thickness_m + e.separator_thickness_m + e.positive_thickness_m
    }
    fn minmax(&self) -> (f64, f64) {
        let lo = self
            .electrolyte_concentration
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let hi = self
            .electrolyte_concentration
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    }
    pub fn state(&self, current: f64) -> Result<BatteryStepResult> {
        require(current.is_finite(), "invalid SPMe current")?;
        let (lo, hi) = self.minmax();
        self.base
            .evaluate(current, self.electrolyte_overpotential(current), lo, hi)
    }
    fn step_electrolyte(&mut self, current: f64, dt: f64) -> Result<()> {
        let e = self.base.config.electrolyte;
        let n = self.electrolyte_concentration.len();
        let dx = self.length() / n as f64;
        let (mut lower, mut diag, mut upper) = (vec![0.0; n - 1], vec![1.0; n], vec![0.0; n - 1]);
        let mut rhs = self.electrolyte_concentration.clone();
        let porosity: Vec<f64> = (0..n)
            .map(|i| region_porosity(&e, (i as f64 + 0.5) * dx))
            .collect();
        let deff: Vec<f64> = porosity
            .iter()
            .map(|p| effective_transport(e.diffusivity_m2_per_s, *p, e.bruggeman_exponent))
            .collect();
        for i in 0..n - 1 {
            let df = 2.0 * deff[i] * deff[i + 1] / (deff[i] + deff[i + 1]).max(1e-30);
            let (ki, kj) = (
                df / (porosity[i] * dx * dx),
                df / (porosity[i + 1] * dx * dx),
            );
            diag[i] += dt * ki;
            upper[i] -= dt * ki;
            lower[i] -= dt * kj;
            diag[i + 1] += dt * kj;
        }
        let area = self.base.config.area_m2;
        let source_n = (1.0 - e.transference_number) * current
            / (FARADAY_CONSTANT * area * e.negative_thickness_m);
        let source_p = -(1.0 - e.transference_number) * current
            / (FARADAY_CONSTANT * area * e.positive_thickness_m);
        for i in 0..n {
            let x = (i as f64 + 0.5) * dx;
            let source = if x < e.negative_thickness_m {
                source_n
            } else if x > e.negative_thickness_m + e.separator_thickness_m {
                source_p
            } else {
                0.0
            };
            rhs[i] += dt * source / porosity[i];
        }
        solve_tridiagonal(&lower, &diag, &upper, &mut rhs)?;
        rhs.iter_mut().for_each(|c| *c = c.max(1e-6));
        self.electrolyte_concentration = rhs;
        Ok(())
    }
    fn electrolyte_overpotential(&self, current: f64) -> f64 {
        let e = &self.base.config.electrolyte;
        let n = self.electrolyte_concentration.len();
        let dx = self.length() / n as f64;
        let mut r = 0.0;
        for i in 0..n {
            let kappa = effective_transport(
                e.conductivity_s_per_m,
                region_porosity(e, (i as f64 + 0.5) * dx),
                e.bruggeman_exponent,
            );
            r += dx / (kappa * self.base.config.area_m2);
        }
        let q = (n / 8).max(1);
        let cn = self.electrolyte_concentration[..q].iter().sum::<f64>() / q as f64;
        let cp = self.electrolyte_concentration[n - q..].iter().sum::<f64>() / q as f64;
        let concentration =
            2.0 * GAS_CONSTANT * self.base.temperature_k * (1.0 - e.transference_number)
                / FARADAY_CONSTANT
                * (cp.max(1e-9) / cn.max(1e-9)).ln();
        current * r - concentration
    }
    fn update_potentials(&mut self, current: f64) {
        let c = &self.base.config;
        let e = &c.electrolyte;
        let n = self.electrolyte_concentration.len();
        let dx = self.length() / n as f64;
        self.electrolyte_potential[0] = 0.0;
        self.solid_potential[0] = (self.base.negative_ocv)(
            self.base
                .negative
                .surface_stoichiometry(self.base.negative_flux(current)),
        );
        for i in 1..n {
            let kappa = effective_transport(
                e.conductivity_s_per_m,
                region_porosity(e, (i as f64 - 0.5) * dx),
                e.bruggeman_exponent,
            );
            self.electrolyte_potential[i] =
                self.electrolyte_potential[i - 1] - current * dx / (kappa * c.area_m2);
            let xpos = (i as f64 + 0.5) * dx;
            let sigma = if xpos < e.negative_thickness_m {
                c.negative.solid_conductivity_s_per_m
            } else if xpos > e.negative_thickness_m + e.separator_thickness_m {
                c.positive.solid_conductivity_s_per_m
            } else {
                1e30
            };
            self.solid_potential[i] = self.solid_potential[i - 1]
                - if sigma > 1e20 {
                    0.0
                } else {
                    current * dx / (sigma * c.area_m2)
                };
        }
    }
    /// Transactional: on error the model is unchanged.
    pub fn step(&mut self, current: f64, dt: f64) -> Result<BatteryStepResult> {
        require(dt > 0.0 && current.is_finite(), "invalid SPMe step")?;
        let mut t = self.clone();
        let (jn, jp) = (t.base.negative_flux(current), t.base.positive_flux(current));
        t.base.negative.step(jn, dt)?;
        t.base.positive.step(jp, dt)?;
        t.step_electrolyte(current, dt)?;
        t.base.time_s += dt;
        t.update_potentials(current);
        let (lo, hi) = t.minmax();
        let eta = t.electrolyte_overpotential(current);
        let result = t.base.evaluate(current, eta, lo, hi)?;
        let result = t.base.finish_step(result, current, dt, jn);
        *self = t;
        Ok(result)
    }
    pub fn electrolyte_concentration(&self) -> &[f64] {
        &self.electrolyte_concentration
    }
    pub fn electrolyte_potential(&self) -> &[f64] {
        &self.electrolyte_potential
    }
    pub fn solid_potential(&self) -> &[f64] {
        &self.solid_potential
    }
    pub fn base(&self) -> &SingleParticleModel {
        &self.base
    }
}

impl BatteryModel for SingleParticleElectrolyteModel {
    fn step(&mut self, current_a: f64, dt_s: f64) -> Result<BatteryStepResult> {
        SingleParticleElectrolyteModel::step(self, current_a, dt_s)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DriveCyclePoint {
    pub duration_s: f64,
    pub current_a: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DriveCycleResult {
    pub samples: Vec<BatteryStepResult>,
    pub discharged_capacity_ah: f64,
    pub electrical_energy_wh: f64,
}

/// Steps each segment in sub-steps of at most `maximum_step_s` (upstream default 1 s).
pub fn simulate_drive_cycle(
    model: &mut dyn BatteryModel,
    segments: &[DriveCyclePoint],
    maximum_step_s: f64,
) -> Result<DriveCycleResult> {
    require(
        maximum_step_s > 0.0,
        "drive-cycle maximum step must be positive",
    )?;
    let mut out = DriveCycleResult::default();
    for seg in segments {
        require(
            seg.duration_s >= 0.0 && seg.current_a.is_finite(),
            "invalid drive-cycle segment",
        )?;
        require(
            seg.duration_s / maximum_step_s < 1e8,
            "drive-cycle segment needs too many steps",
        )?;
        let mut remaining = seg.duration_s;
        while remaining > 1e-12 {
            let dt = maximum_step_s.min(remaining);
            let s = model.step(seg.current_a, dt)?;
            out.samples.push(s);
            if seg.current_a > 0.0 {
                out.discharged_capacity_ah += seg.current_a * dt / 3600.0;
            }
            out.electrical_energy_wh += s.terminal_voltage_v * seg.current_a * dt / 3600.0;
            remaining -= dt;
        }
    }
    Ok(out)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatteryEisConfig {
    pub series_resistance_ohm: f64,
    pub negative_charge_transfer_resistance_ohm: f64,
    pub positive_charge_transfer_resistance_ohm: f64,
    pub negative_double_layer_capacitance_f: f64,
    pub positive_double_layer_capacitance_f: f64,
    pub diffusion_warburg_ohm_sqrt_s: f64,
    pub diffusion_time_s: f64,
}

impl Default for BatteryEisConfig {
    fn default() -> Self {
        Self {
            series_resistance_ohm: 0.010,
            negative_charge_transfer_resistance_ohm: 0.020,
            positive_charge_transfer_resistance_ohm: 0.025,
            negative_double_layer_capacitance_f: 2.0,
            positive_double_layer_capacitance_f: 2.0,
            diffusion_warburg_ohm_sqrt_s: 0.015,
            diffusion_time_s: 200.0,
        }
    }
}

/// `R0 + two parallel RC arcs + finite-length Warburg`.
/// Overflow-free complex tanh: `(sinh 2x + i sin 2y) / (cosh 2x + cos 2y)` with the large-|x| limit.
fn stable_tanh(z: C) -> C {
    if z.re.abs() > 20.0 {
        let e = (-2.0 * z.re.abs()).exp();
        return C::new(z.re.signum(), 2.0 * (2.0 * z.im).sin() * e);
    }
    z.tanh()
}

pub fn battery_impedance(frequency_hz: f64, c: &BatteryEisConfig) -> Result<C> {
    require(frequency_hz > 0.0, "battery EIS frequency must be positive")?;
    let omega = 2.0 * PI * frequency_hz;
    let j = C::new(0.0, 1.0);
    let one = C::new(1.0, 0.0);
    let parallel_rc = |r: f64, cap: f64| C::new(r, 0.0) / (one + j * omega * r * cap);
    let root = (j * omega * c.diffusion_time_s).sqrt();
    let warburg = C::new(c.diffusion_warburg_ohm_sqrt_s, 0.0) / (j * omega).sqrt()
        * if root.norm() > 1e-12 {
            stable_tanh(root) / root
        } else {
            one
        };
    Ok(C::new(c.series_resistance_ohm, 0.0)
        + parallel_rc(
            c.negative_charge_transfer_resistance_ohm,
            c.negative_double_layer_capacitance_f,
        )
        + parallel_rc(
            c.positive_charge_transfer_resistance_ohm,
            c.positive_double_layer_capacitance_f,
        )
        + warburg)
}

pub fn battery_impedance_spectrum(frequency_hz: &[f64], c: &BatteryEisConfig) -> Result<Vec<C>> {
    frequency_hz
        .iter()
        .map(|f| battery_impedance(*f, c))
        .collect()
}

pub(crate) fn invalid(msg: &str) -> SolverError {
    SolverError::InvalidInput(msg.into())
}
