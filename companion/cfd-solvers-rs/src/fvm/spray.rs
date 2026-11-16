//! Coupled spray parcels (Stokes drag, D^2-law evaporation, latent heat), ported from upstream
//! `cfd/solvers/fvm/spray` (MIT, cfd_solvers contributors).
use super::{
    evaporate_d2_law,
    particles::{finite_vec, nearest_cell},
    PolyMesh,
};
use crate::{core::Vec3, require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SprayParcel {
    pub position: Vec3,
    pub velocity: Vec3,
    pub diameter: f64,
    pub mass: f64,
    pub temperature: f64,
}

impl Default for SprayParcel {
    fn default() -> Self {
        Self {
            position: Vec3::default(),
            velocity: Vec3::default(),
            diameter: 1e-4,
            mass: 1e-9,
            temperature: 300.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SprayInjection {
    pub position: Vec3,
    pub velocity: Vec3,
    pub diameter: f64,
    pub temperature: f64,
    pub mass_flow_rate: f64,
    pub start_time: f64,
    pub duration: f64,
}

impl Default for SprayInjection {
    fn default() -> Self {
        Self {
            position: Vec3::default(),
            velocity: Vec3::default(),
            diameter: 1e-4,
            temperature: 300.0,
            mass_flow_rate: 1e-6,
            start_time: 0.0,
            duration: 1e-2,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SprayConfig {
    pub dt: f64,
    pub liquid_density: f64,
    pub carrier_dynamic_viscosity: f64,
    pub evaporation_constant: f64,
    pub latent_heat: f64,
    pub thermal_relaxation_time: f64,
    pub minimum_diameter: f64,
    pub maximum_parcels: usize,
}

impl Default for SprayConfig {
    fn default() -> Self {
        Self {
            dt: 1e-4,
            liquid_density: 1000.0,
            carrier_dynamic_viscosity: 1.8e-5,
            evaporation_constant: 1e-8,
            latent_heat: 2.2e6,
            thermal_relaxation_time: 1e-3,
            minimum_diameter: 1e-8,
            maximum_parcels: 100_000,
        }
    }
}

/// Parcel spray on a PolyMesh; the carrier receives vapour mass, drag reaction and latent-heat sinks.
#[derive(Clone, Debug)]
pub struct SprayInjectionEvaporation {
    mesh: PolyMesh,
    config: SprayConfig,
    injection: Option<SprayInjection>,
    parcels: Vec<SprayParcel>,
    vapor_mass_source: Vec<f64>,
    carrier_energy_source: Vec<f64>,
    carrier_momentum_source: Vec<Vec3>,
    injected_mass: f64,
    evaporated_mass: f64,
    time: f64,
    steps: usize,
}

fn positive_finite(v: f64) -> bool {
    v > 0.0 && v.is_finite()
}
fn nonnegative_finite(v: f64) -> bool {
    v >= 0.0 && v.is_finite()
}

impl SprayInjectionEvaporation {
    pub fn new(mesh: PolyMesh, config: SprayConfig) -> Result<Self> {
        let c = config;
        require(
            positive_finite(c.dt)
                && positive_finite(c.liquid_density)
                && positive_finite(c.carrier_dynamic_viscosity)
                && nonnegative_finite(c.evaporation_constant)
                && nonnegative_finite(c.latent_heat)
                && positive_finite(c.thermal_relaxation_time)
                && positive_finite(c.minimum_diameter)
                && c.maximum_parcels > 0,
            "invalid spray controls",
        )?;
        let n = mesh.cell_count();
        Ok(Self {
            mesh,
            config,
            injection: None,
            parcels: Vec::new(),
            vapor_mass_source: vec![0.0; n],
            carrier_energy_source: vec![0.0; n],
            carrier_momentum_source: vec![Vec3::default(); n],
            injected_mass: 0.0,
            evaporated_mass: 0.0,
            time: 0.0,
            steps: 0,
        })
    }

    pub fn set_injection(&mut self, i: SprayInjection) -> Result<()> {
        require(
            finite_vec(i.position) && finite_vec(i.velocity),
            "spray injection position and velocity must be finite",
        )?;
        require(
            i.diameter >= self.config.minimum_diameter
                && i.diameter.is_finite()
                && nonnegative_finite(i.temperature)
                && nonnegative_finite(i.mass_flow_rate)
                && nonnegative_finite(i.start_time)
                && positive_finite(i.duration),
            "invalid spray injection",
        )?;
        self.injection = Some(i);
        Ok(())
    }
    pub fn clear_injection(&mut self) {
        self.injection = None;
    }
    pub fn add_parcel(&mut self, p: SprayParcel) -> Result<()> {
        require(
            self.parcels.len() < self.config.maximum_parcels,
            "spray parcel capacity exceeded",
        )?;
        require(
            finite_vec(p.position) && finite_vec(p.velocity),
            "spray parcel position and velocity must be finite",
        )?;
        require(
            p.diameter >= self.config.minimum_diameter
                && p.diameter.is_finite()
                && positive_finite(p.mass)
                && p.temperature.is_finite(),
            "invalid spray parcel",
        )?;
        self.parcels.push(p);
        Ok(())
    }
    pub fn clear_parcels(&mut self) {
        self.parcels.clear();
    }

    /// Advance with a uniform carrier temperature.
    pub fn step_uniform(
        &mut self,
        carrier_velocity: &[Vec3],
        carrier_temperature: f64,
    ) -> Result<()> {
        require(
            nonnegative_finite(carrier_temperature),
            "spray carrier temperature must be non-negative and finite",
        )?;
        let t = vec![carrier_temperature; self.mesh.cell_count()];
        self.step(carrier_velocity, &t)
    }

    pub fn step(&mut self, carrier_velocity: &[Vec3], carrier_temperature: &[f64]) -> Result<()> {
        let n = self.mesh.cell_count();
        require(
            carrier_velocity.len() == n && carrier_temperature.len() == n,
            "spray carrier field size mismatch",
        )?;
        require(
            carrier_velocity.iter().all(|v| finite_vec(*v)),
            "spray carrier velocity must be finite",
        )?;
        require(
            carrier_temperature.iter().all(|t| nonnegative_finite(*t)),
            "spray carrier temperature must be non-negative and finite",
        )?;
        self.vapor_mass_source.fill(0.0);
        self.carrier_energy_source.fill(0.0);
        self.carrier_momentum_source.fill(Vec3::default());
        let c = self.config;
        if let Some(inj) = self.injection {
            if self.time >= inj.start_time
                && self.time < inj.start_time + inj.duration
                && inj.mass_flow_rate > 0.0
            {
                self.add_parcel(SprayParcel {
                    position: inj.position,
                    velocity: inj.velocity,
                    diameter: inj.diameter,
                    mass: inj.mass_flow_rate * c.dt,
                    temperature: inj.temperature,
                })?;
                self.injected_mass += inj.mass_flow_rate * c.dt;
            }
        }
        let mut active = Vec::with_capacity(self.parcels.len());
        for mut p in self.parcels.iter().copied() {
            let cell = nearest_cell(&self.mesh, p.position);
            let old_mass = p.mass;
            let tau =
                c.liquid_density * p.diameter * p.diameter / (18.0 * c.carrier_dynamic_viscosity);
            let decay = (-c.dt / tau).exp();
            let old_velocity = p.velocity;
            p.velocity = carrier_velocity[cell] + (p.velocity - carrier_velocity[cell]) * decay;
            let drag_on_parcel = (p.velocity - old_velocity) * (old_mass / c.dt);
            let volume = self.mesh.cells()[cell].volume;
            self.carrier_momentum_source[cell] -= drag_on_parcel / volume;
            let thermal_decay = (-c.dt / c.thermal_relaxation_time).exp();
            p.temperature = carrier_temperature[cell]
                + (p.temperature - carrier_temperature[cell]) * thermal_decay;
            let new_d = evaporate_d2_law(p.diameter, c.evaporation_constant, c.dt)?;
            let ratio = new_d / p.diameter;
            p.mass = old_mass * ratio * ratio * ratio;
            let vaporized = old_mass - p.mass;
            if vaporized > 0.0 {
                self.vapor_mass_source[cell] += vaporized / (c.dt * volume);
                self.carrier_energy_source[cell] -= c.latent_heat * vaporized / (c.dt * volume);
                self.evaporated_mass += vaporized;
            }
            p.diameter = new_d;
            p.position += p.velocity * c.dt;
            if p.diameter > c.minimum_diameter && p.mass > 0.0 {
                active.push(p);
            }
        }
        self.parcels = active;
        self.time += c.dt;
        self.steps += 1;
        Ok(())
    }

    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn config(&self) -> &SprayConfig {
        &self.config
    }
    pub fn parcels(&self) -> &[SprayParcel] {
        &self.parcels
    }
    pub fn vapor_mass_source(&self) -> &[f64] {
        &self.vapor_mass_source
    }
    pub fn carrier_momentum_source(&self) -> &[Vec3] {
        &self.carrier_momentum_source
    }
    pub fn carrier_energy_source(&self) -> &[f64] {
        &self.carrier_energy_source
    }
    pub fn total_liquid_mass(&self) -> f64 {
        self.parcels.iter().map(|p| p.mass).sum()
    }
    pub fn injected_mass(&self) -> f64 {
        self.injected_mass
    }
    pub fn evaporated_mass(&self) -> f64 {
        self.evaporated_mass
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn steps(&self) -> usize {
        self.steps
    }
}
