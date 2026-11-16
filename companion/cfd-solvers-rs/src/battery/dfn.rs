//! Pseudo-two-dimensional Doyle-Fuller-Newman reference model.
#![allow(clippy::needless_range_loop)]
use super::{
    clamp01, default_graphite_nmc_config, effective_transport, graphite_ocv_v, nmc_ocv_v,
    solve_tridiagonal, symmetric_butler_volmer_overpotential_v, OcvFunction, ParticleConfig,
    SphericalDiffusionParticle, FARADAY_CONSTANT, GAS_CONSTANT,
};
use crate::{require, Result};
use std::sync::Arc;

#[derive(Clone)]
pub struct DfnConfig {
    pub area_m2: f64,
    pub negative_thickness_m: f64,
    pub separator_thickness_m: f64,
    pub positive_thickness_m: f64,
    pub negative_particle: ParticleConfig,
    pub positive_particle: ParticleConfig,
    pub negative_solid_conductivity_s_per_m: f64,
    pub positive_solid_conductivity_s_per_m: f64,
    pub negative_exchange_current_density_a_per_m2: f64,
    pub positive_exchange_current_density_a_per_m2: f64,
    pub negative_porosity: f64,
    pub positive_porosity: f64,
    pub separator_porosity: f64,
    pub bruggeman_exponent: f64,
    pub negative_active_volume_fraction: f64,
    pub positive_active_volume_fraction: f64,
    pub electrolyte_initial_concentration_mol_per_m3: f64,
    pub electrolyte_diffusivity_m2_per_s: f64,
    pub electrolyte_conductivity_s_per_m: f64,
    pub transference_number: f64,
    pub negative_ocv_v: Option<OcvFunction>,
    pub positive_ocv_v: Option<OcvFunction>,
    pub negative_nodes: usize,
    pub separator_nodes: usize,
    pub positive_nodes: usize,
    pub contact_resistance_ohm: f64,
    pub initial_temperature_k: f64,
}

impl Default for DfnConfig {
    fn default() -> Self {
        Self {
            area_m2: 0.010,
            negative_thickness_m: 80e-6,
            separator_thickness_m: 25e-6,
            positive_thickness_m: 75e-6,
            negative_particle: ParticleConfig::default(),
            positive_particle: ParticleConfig::default(),
            negative_solid_conductivity_s_per_m: 100.0,
            positive_solid_conductivity_s_per_m: 10.0,
            negative_exchange_current_density_a_per_m2: 4.0,
            positive_exchange_current_density_a_per_m2: 4.0,
            negative_porosity: 0.30,
            positive_porosity: 0.30,
            separator_porosity: 0.50,
            bruggeman_exponent: 1.5,
            negative_active_volume_fraction: 0.62,
            positive_active_volume_fraction: 0.58,
            electrolyte_initial_concentration_mol_per_m3: 1000.0,
            electrolyte_diffusivity_m2_per_s: 2e-10,
            electrolyte_conductivity_s_per_m: 1.0,
            transference_number: 0.38,
            negative_ocv_v: None,
            positive_ocv_v: None,
            negative_nodes: 24,
            separator_nodes: 9,
            positive_nodes: 22,
            contact_resistance_ohm: 0.012,
            initial_temperature_k: 298.15,
        }
    }
}

pub fn default_graphite_nmc_dfn_config() -> DfnConfig {
    let b = default_graphite_nmc_config();
    DfnConfig {
        area_m2: b.area_m2,
        negative_thickness_m: b.negative.thickness_m,
        positive_thickness_m: b.positive.thickness_m,
        negative_particle: b.negative.particle,
        positive_particle: b.positive.particle,
        negative_solid_conductivity_s_per_m: b.negative.solid_conductivity_s_per_m,
        positive_solid_conductivity_s_per_m: b.positive.solid_conductivity_s_per_m,
        negative_exchange_current_density_a_per_m2: b.negative.exchange_current_density_a_per_m2,
        positive_exchange_current_density_a_per_m2: b.positive.exchange_current_density_a_per_m2,
        negative_porosity: b.electrolyte.negative_porosity,
        positive_porosity: b.electrolyte.positive_porosity,
        separator_porosity: b.electrolyte.separator_porosity,
        bruggeman_exponent: b.electrolyte.bruggeman_exponent,
        negative_active_volume_fraction: b.negative.active_volume_fraction,
        positive_active_volume_fraction: b.positive.active_volume_fraction,
        electrolyte_initial_concentration_mol_per_m3: b
            .electrolyte
            .initial_concentration_mol_per_m3,
        electrolyte_diffusivity_m2_per_s: b.electrolyte.diffusivity_m2_per_s,
        electrolyte_conductivity_s_per_m: b.electrolyte.conductivity_s_per_m,
        transference_number: b.electrolyte.transference_number,
        negative_ocv_v: Some(Arc::new(graphite_ocv_v)),
        positive_ocv_v: Some(Arc::new(nmc_ocv_v)),
        negative_nodes: 20,
        separator_nodes: 9,
        positive_nodes: 18,
        contact_resistance_ohm: b.contact_resistance_ohm,
        initial_temperature_k: b.initial_temperature_k,
        ..Default::default()
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DfnStepResult {
    pub time_s: f64,
    pub current_a: f64,
    pub terminal_voltage_v: f64,
    pub open_circuit_voltage_v: f64,
    pub reaction_overpotential_v: f64,
    /// Solid + electrolyte + contact.
    pub ohmic_drop_v: f64,
    pub negative_surface_stoichiometry: f64,
    pub positive_surface_stoichiometry: f64,
    pub electrolyte_min_concentration_mol_per_m3: f64,
    pub electrolyte_max_concentration_mol_per_m3: f64,
    pub temperature_k: f64,
    pub heat_generation_w: f64,
}

/// Particles at every electrode node, a 1-D electrolyte and solid/electrolyte potentials; the
/// reaction current is partitioned by local exchange-current weights.
#[derive(Clone)]
pub struct DoyleFullerNewmanModel {
    config: DfnConfig,
    negative_ocv: OcvFunction,
    positive_ocv: OcvFunction,
    negative_particles: Vec<SphericalDiffusionParticle>,
    positive_particles: Vec<SphericalDiffusionParticle>,
    electrolyte_concentration: Vec<f64>,
    electrolyte_potential: Vec<f64>,
    solid_potential: Vec<f64>,
    local_reaction_current: Vec<f64>,
    temperature_k: f64,
    time_s: f64,
}

impl DoyleFullerNewmanModel {
    pub fn new(config: DfnConfig) -> Result<Self> {
        require(config.area_m2 > 0.0, "invalid DFN area")?;
        require(
            config.negative_thickness_m > 0.0
                && config.separator_thickness_m > 0.0
                && config.positive_thickness_m > 0.0,
            "invalid DFN thickness",
        )?;
        require(
            config.negative_nodes >= 2 && config.separator_nodes >= 1 && config.positive_nodes >= 2,
            "DFN needs at least two electrode nodes per side",
        )?;
        require(
            config.negative_nodes + config.separator_nodes + config.positive_nodes <= 1 << 16,
            "DFN node count too large",
        )?;
        SphericalDiffusionParticle::new(config.negative_particle)?;
        SphericalDiffusionParticle::new(config.positive_particle)?;
        let negative_ocv = config
            .negative_ocv_v
            .clone()
            .unwrap_or_else(|| Arc::new(graphite_ocv_v));
        let positive_ocv = config
            .positive_ocv_v
            .clone()
            .unwrap_or_else(|| Arc::new(nmc_ocv_v));
        let mut m = Self {
            config,
            negative_ocv,
            positive_ocv,
            negative_particles: vec![],
            positive_particles: vec![],
            electrolyte_concentration: vec![],
            electrolyte_potential: vec![],
            solid_potential: vec![],
            local_reaction_current: vec![],
            temperature_k: 0.0,
            time_s: 0.0,
        };
        m.reset()?;
        Ok(m)
    }
    pub fn reset(&mut self) -> Result<()> {
        let c = &self.config;
        let n = self.total_nodes();
        self.negative_particles =
            vec![SphericalDiffusionParticle::new(c.negative_particle)?; c.negative_nodes];
        self.positive_particles =
            vec![SphericalDiffusionParticle::new(c.positive_particle)?; c.positive_nodes];
        self.electrolyte_concentration = vec![c.electrolyte_initial_concentration_mol_per_m3; n];
        self.electrolyte_potential = vec![0.0; n];
        self.solid_potential = vec![0.0; n];
        self.local_reaction_current = vec![0.0; n];
        self.temperature_k = c.initial_temperature_k;
        self.time_s = 0.0;
        Ok(())
    }
    fn total_nodes(&self) -> usize {
        self.config.negative_nodes + self.config.separator_nodes + self.config.positive_nodes
    }
    fn node_spacing(&self) -> f64 {
        let c = &self.config;
        (c.negative_thickness_m + c.separator_thickness_m + c.positive_thickness_m)
            / self.total_nodes() as f64
    }
    fn region_porosity(&self, x: f64) -> f64 {
        let c = &self.config;
        if x < c.negative_thickness_m {
            c.negative_porosity
        } else if x < c.negative_thickness_m + c.separator_thickness_m {
            c.separator_porosity
        } else {
            c.positive_porosity
        }
    }
    fn solid_conductivity(&self, x: f64) -> f64 {
        let c = &self.config;
        if x < c.negative_thickness_m {
            c.negative_solid_conductivity_s_per_m
        } else if x < c.negative_thickness_m + c.separator_thickness_m {
            0.0
        } else {
            c.positive_solid_conductivity_s_per_m
        }
    }
    fn update_reaction_partition(&mut self, current: f64) {
        let (nn, ns) = (self.config.negative_nodes, self.config.separator_nodes);
        let i_areal = current / self.config.area_m2;
        self.local_reaction_current.fill(0.0);
        let partition =
            |particles: &[SphericalDiffusionParticle], i0: f64, sign: f64, out: &mut [f64]| {
                let w: Vec<f64> = particles
                    .iter()
                    .map(|p| {
                        let x = clamp01(p.surface_stoichiometry(0.0));
                        i0 * 2.0 * (x * (1.0 - x)).sqrt()
                    })
                    .collect();
                let mut sum = w.iter().fold(0.0, |a, b| a + b);
                if sum <= 0.0 {
                    sum = 1.0;
                }
                for (o, wi) in out.iter_mut().zip(&w) {
                    *o = sign * i_areal * wi / sum;
                }
            };
        partition(
            &self.negative_particles,
            self.config.negative_exchange_current_density_a_per_m2,
            1.0,
            &mut self.local_reaction_current[..nn],
        );
        partition(
            &self.positive_particles,
            self.config.positive_exchange_current_density_a_per_m2,
            -1.0,
            &mut self.local_reaction_current[nn + ns..],
        );
    }
    fn advance_particles(&mut self, dt: f64) -> Result<()> {
        let c = &self.config;
        let dx = self.node_spacing();
        let a_neg = 3.0 * c.negative_active_volume_fraction / c.negative_particle.radius_m;
        let a_pos = 3.0 * c.positive_active_volume_fraction / c.positive_particle.radius_m;
        let offset = c.negative_nodes + c.separator_nodes;
        for (i, p) in self.negative_particles.iter_mut().enumerate() {
            p.step(
                self.local_reaction_current[i] / (a_neg * dx * FARADAY_CONSTANT),
                dt,
            )?;
        }
        for (i, p) in self.positive_particles.iter_mut().enumerate() {
            p.step(
                self.local_reaction_current[offset + i] / (a_pos * dx * FARADAY_CONSTANT),
                dt,
            )?;
        }
        Ok(())
    }
    fn advance_electrolyte(&mut self, dt: f64) -> Result<()> {
        let n = self.total_nodes();
        let dx = self.node_spacing();
        let c = &self.config;
        let porosity: Vec<f64> = (0..n)
            .map(|i| self.region_porosity((i as f64 + 0.5) * dx))
            .collect();
        let deff: Vec<f64> = porosity
            .iter()
            .map(|p| {
                effective_transport(c.electrolyte_diffusivity_m2_per_s, *p, c.bruggeman_exponent)
            })
            .collect();
        let mut face = vec![0.0; n + 1];
        for i in 1..n {
            let df = 2.0 * deff[i - 1] * deff[i] / (deff[i - 1] + deff[i]).max(1e-30);
            face[i] = df / (dx * dx);
        }
        let (mut lower, mut diag, mut upper, mut rhs) = (
            vec![0.0; n - 1],
            vec![1.0; n],
            vec![0.0; n - 1],
            vec![0.0; n],
        );
        for i in 0..n {
            let (left, right) = (face[i], face[i + 1]);
            diag[i] = 1.0 + dt / porosity[i] * (left + right);
            if i > 0 {
                lower[i - 1] = -dt / porosity[i] * left;
            }
            if i + 1 < n {
                upper[i] = -dt / porosity[i] * right;
            }
            let source = (1.0 - c.transference_number) * self.local_reaction_current[i]
                / (dx * FARADAY_CONSTANT);
            rhs[i] = self.electrolyte_concentration[i] + dt * source / porosity[i];
        }
        solve_tridiagonal(&lower, &diag, &upper, &mut rhs)?;
        rhs.iter_mut().for_each(|v| *v = v.max(1e-6));
        self.electrolyte_concentration = rhs;
        Ok(())
    }
    fn update_potentials(&mut self, current: f64) {
        let n = self.total_nodes();
        let dx = self.node_spacing();
        let c = &self.config;
        let mut i_e = vec![0.0; n];
        let mut running = 0.0;
        for i in 0..n {
            i_e[i] = running;
            running += self.local_reaction_current[i];
        }
        self.electrolyte_potential[0] = 0.0;
        for i in 1..n {
            let kappa = effective_transport(
                c.electrolyte_conductivity_s_per_m,
                self.region_porosity((i as f64 - 0.5) * dx),
                c.bruggeman_exponent,
            );
            let i_face = 0.5 * (i_e[i - 1] + i_e[i]);
            self.electrolyte_potential[i] =
                self.electrolyte_potential[i - 1] - i_face * dx / (kappa * c.area_m2);
        }
        self.solid_potential[0] =
            (self.negative_ocv)(self.negative_particles[0].surface_stoichiometry(0.0));
        let mut i_s = current / c.area_m2;
        for i in 1..n {
            let sigma = self.solid_conductivity((i as f64 - 0.5) * dx);
            if sigma > 0.0 {
                self.solid_potential[i] =
                    self.solid_potential[i - 1] - i_s * dx / (sigma * c.area_m2);
                i_s -= self.local_reaction_current[i - 1];
            } else {
                self.solid_potential[i] = self.solid_potential[i - 1];
            }
        }
    }
    fn evaluate(&self, current: f64) -> Result<DfnStepResult> {
        let c = &self.config;
        let (nn, ns, np) = (c.negative_nodes, c.separator_nodes, c.positive_nodes);
        let dx = self.node_spacing();
        let a_neg = 3.0 * c.negative_active_volume_fraction / c.negative_particle.radius_m;
        let a_pos = 3.0 * c.positive_active_volume_fraction / c.positive_particle.radius_m;
        let (mut ocv_neg, mut ocv_pos, mut reaction) = (0.0, 0.0, 0.0);
        for i in 0..nn {
            let x = clamp01(self.negative_particles[i].surface_stoichiometry(0.0));
            ocv_neg += (self.negative_ocv)(x);
            let j = self.local_reaction_current[i] / (a_neg * dx);
            let i0 = c.negative_exchange_current_density_a_per_m2 * 2.0 * (x * (1.0 - x)).sqrt();
            reaction +=
                symmetric_butler_volmer_overpotential_v(j, i0.max(1e-12), self.temperature_k)?
                    .abs();
        }
        for i in 0..np {
            let x = clamp01(self.positive_particles[i].surface_stoichiometry(0.0));
            ocv_pos += (self.positive_ocv)(x);
            let j = self.local_reaction_current[nn + ns + i] / (a_pos * dx);
            let i0 = c.positive_exchange_current_density_a_per_m2 * 2.0 * (x * (1.0 - x)).sqrt();
            reaction +=
                symmetric_butler_volmer_overpotential_v(j, i0.max(1e-12), self.temperature_k)?
                    .abs();
        }
        ocv_neg /= nn as f64;
        ocv_pos /= np as f64;
        reaction /= (nn + np) as f64;
        let ocv = ocv_pos - ocv_neg;
        let total = self.total_nodes();
        let mut r_e = 0.0;
        for i in 0..total {
            let kappa = effective_transport(
                c.electrolyte_conductivity_s_per_m,
                self.region_porosity((i as f64 + 0.5) * dx),
                c.bruggeman_exponent,
            );
            r_e += dx / (kappa * c.area_m2);
        }
        let r_s = c.negative_thickness_m / (c.negative_solid_conductivity_s_per_m * c.area_m2)
            + c.positive_thickness_m / (c.positive_solid_conductivity_s_per_m * c.area_m2);
        let q = (total / 8).max(1);
        let cn = self.electrolyte_concentration[..q].iter().sum::<f64>() / q as f64;
        let cp = self.electrolyte_concentration[total - q..]
            .iter()
            .sum::<f64>()
            / q as f64;
        let conc = 2.0 * GAS_CONSTANT * self.temperature_k * (1.0 - c.transference_number)
            / FARADAY_CONSTANT
            * (cp.max(1e-9) / cn.max(1e-9)).ln();
        let eta_e = current * r_e - conc;
        let ohmic = current * (c.contact_resistance_ohm + r_s);
        let voltage = ocv
            - (reaction + eta_e.abs()).copysign(if current == 0.0 { 1.0 } else { current })
            - ohmic;
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
        Ok(DfnStepResult {
            time_s: self.time_s,
            current_a: current,
            terminal_voltage_v: voltage,
            open_circuit_voltage_v: ocv,
            reaction_overpotential_v: reaction,
            ohmic_drop_v: ohmic + eta_e,
            negative_surface_stoichiometry: self.negative_particles[0].surface_stoichiometry(0.0),
            positive_surface_stoichiometry: self.positive_particles[np - 1]
                .surface_stoichiometry(0.0),
            electrolyte_min_concentration_mol_per_m3: lo,
            electrolyte_max_concentration_mol_per_m3: hi,
            temperature_k: self.temperature_k,
            heat_generation_w: current.abs() * (ocv - voltage).abs().max(0.0),
        })
    }
    /// Transactional: on error the model is unchanged.
    pub fn step(&mut self, current: f64, dt: f64) -> Result<DfnStepResult> {
        require(dt > 0.0 && current.is_finite(), "invalid DFN step")?;
        let mut t = self.clone();
        t.update_reaction_partition(current);
        t.advance_particles(dt)?;
        t.advance_electrolyte(dt)?;
        t.time_s += dt;
        t.update_potentials(current);
        let r = t.evaluate(current)?;
        *self = t;
        Ok(r)
    }
    pub fn state(&self, current: f64) -> Result<DfnStepResult> {
        require(current.is_finite(), "invalid DFN current")?;
        self.evaluate(current)
    }
    pub fn total_lithium_mol(&self) -> f64 {
        let c = &self.config;
        let dx = self.node_spacing();
        let a = c.area_m2;
        let mut total = 0.0;
        for p in &self.negative_particles {
            total += p.average_concentration() * c.negative_active_volume_fraction * a * dx;
        }
        for p in &self.positive_particles {
            total += p.average_concentration() * c.positive_active_volume_fraction * a * dx;
        }
        for (i, ce) in self.electrolyte_concentration.iter().enumerate() {
            total += ce * self.region_porosity((i as f64 + 0.5) * dx) * a * dx;
        }
        total
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub fn config(&self) -> &DfnConfig {
        &self.config
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
    pub fn local_reaction_current(&self) -> &[f64] {
        &self.local_reaction_current
    }
    pub fn negative_particles(&self) -> &[SphericalDiffusionParticle] {
        &self.negative_particles
    }
    pub fn positive_particles(&self) -> &[SphericalDiffusionParticle] {
        &self.positive_particles
    }
}
