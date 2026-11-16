//! Reacting-flow building blocks ported from upstream `cfd/solvers/fvm`
//! (`species_transport`, `reacting`, `premixed_flame1d`, `flame_models`, `soot_radiation`;
//! MIT, cfd_solvers contributors).
use super::{optically_thin_radiation_source, PolyMesh, ScalarTransport, ScalarTransportConfig};
use crate::{chemistry::ReactionNetwork, core::Vec3, require, Result, SolverError};

/// Independent implicit scalar transport per species on one mesh and one face flux.
pub struct MultiSpeciesTransport {
    species: Vec<ScalarTransport>,
}

impl MultiSpeciesTransport {
    pub fn new(mesh: PolyMesh, diffusivity: &[f64], dt: f64) -> Result<Self> {
        require(
            !diffusivity.is_empty(),
            "species transport requires species",
        )?;
        let mut species = Vec::with_capacity(diffusivity.len());
        for d in diffusivity {
            require(d.is_finite() && *d >= 0.0, "invalid species diffusivity")?;
            species.push(ScalarTransport::new(
                mesh.clone(),
                ScalarTransportConfig {
                    dt,
                    diffusivity: *d,
                    ..Default::default()
                },
            )?);
        }
        Ok(Self { species })
    }
    fn get(&self, s: usize) -> Result<&ScalarTransport> {
        self.species
            .get(s)
            .ok_or_else(|| SolverError::InvalidInput("species index".into()))
    }
    pub fn initialize(&mut self, s: usize, value: impl Fn(Vec3) -> f64) -> Result<()> {
        self.species
            .get_mut(s)
            .ok_or_else(|| SolverError::InvalidInput("species index".into()))?
            .initialize_profile(value)
    }
    pub fn set_face_flux(&mut self, flux: &[f64]) -> Result<()> {
        for s in &mut self.species {
            s.set_face_flux(flux.to_vec())?;
        }
        Ok(())
    }
    pub fn step(&mut self) -> Result<()> {
        for s in &mut self.species {
            s.step()?;
        }
        Ok(())
    }
    pub fn run(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }
    pub fn species_count(&self) -> usize {
        self.species.len()
    }
    pub fn values(&self, s: usize) -> Result<&[f64]> {
        Ok(self.get(s)?.values())
    }
    pub fn inventory(&self, s: usize) -> Result<f64> {
        Ok(self.get(s)?.volume_integral())
    }
}

/// Chemical source terms of every cell from its concentrations and temperature.
pub fn reacting_cell_sources(
    network: &ReactionNetwork,
    cell_concentrations: &[Vec<f64>],
    temperature: &[f64],
) -> Result<Vec<Vec<f64>>> {
    require(
        cell_concentrations.len() == temperature.len(),
        "reacting CFD cell-state size mismatch",
    )?;
    cell_concentrations
        .iter()
        .zip(temperature)
        .map(|(c, t)| network.source_terms(c, *t))
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PremixedFlame1DConfig {
    pub cells: usize,
    pub length: f64,
    pub diffusivity: f64,
    pub reaction_rate: f64,
    pub dt: f64,
    pub unburned_temperature: f64,
    pub burned_temperature: f64,
}

impl Default for PremixedFlame1DConfig {
    fn default() -> Self {
        Self {
            cells: 400,
            length: 1.0,
            diffusivity: 1e-3,
            reaction_rate: 1.0,
            dt: 1e-4,
            unburned_temperature: 300.0,
            burned_temperature: 1800.0,
        }
    }
}

/// Fisher-KPP progress-variable laminar flame reference (c = 1 burned, 0 unburned),
/// explicit in time with the diffusion stability bound enforced at construction.
pub struct PremixedFlame1D {
    config: PremixedFlame1DConfig,
    c: Vec<f64>,
    next: Vec<f64>,
}

impl PremixedFlame1D {
    pub fn new(config: PremixedFlame1DConfig) -> Result<Self> {
        require(config.cells >= 4, "invalid flame controls")?;
        let dx = config.length / (config.cells - 1) as f64;
        require(
            config.length > 0.0
                && config.diffusivity > 0.0
                && config.reaction_rate > 0.0
                && config.dt > 0.0
                && config.dt <= 0.45 * dx * dx / config.diffusivity,
            "invalid flame controls",
        )?;
        Ok(Self {
            config,
            c: vec![0.0; config.cells],
            next: vec![0.0; config.cells],
        })
    }
    fn dx(&self) -> f64 {
        self.config.length / (self.config.cells - 1) as f64
    }
    pub fn initialize_front(&mut self, location: f64, thickness: f64) -> Result<()> {
        require(thickness > 0.0, "invalid flame thickness")?;
        let dx = self.dx();
        for (i, c) in self.c.iter_mut().enumerate() {
            *c = 1.0 / (1.0 + ((i as f64 * dx - location) / thickness).exp());
        }
        Ok(())
    }
    pub fn step(&mut self, count: usize) {
        let dx = self.dx();
        let a = self.config.diffusivity / (dx * dx);
        let n = self.config.cells;
        for _ in 0..count {
            self.next[0] = 1.0;
            self.next[n - 1] = 0.0;
            for i in 1..n - 1 {
                let diffusion = a * (self.c[i - 1] - 2.0 * self.c[i] + self.c[i + 1]);
                let reaction = self.config.reaction_rate * self.c[i] * (1.0 - self.c[i]);
                self.next[i] =
                    (self.c[i] + self.config.dt * (diffusion + reaction)).clamp(0.0, 1.0);
            }
            std::mem::swap(&mut self.c, &mut self.next);
        }
    }
    /// Position where the progress crosses 0.5, linearly interpolated.
    pub fn front_location(&self) -> f64 {
        let dx = self.dx();
        for i in 1..self.c.len() {
            if self.c[i] <= 0.5 && self.c[i - 1] > 0.5 {
                let f = (self.c[i - 1] - 0.5) / (self.c[i - 1] - self.c[i]);
                return ((i - 1) as f64 + f) * dx;
            }
        }
        if self.c.last().is_some_and(|c| *c > 0.5) {
            self.config.length
        } else {
            0.0
        }
    }
    /// Fisher-KPP minimum wave speed `2 sqrt(D k)`.
    pub fn theoretical_speed(&self) -> f64 {
        2.0 * (self.config.diffusivity * self.config.reaction_rate).sqrt()
    }
    pub fn progress(&self) -> &[f64] {
        &self.c
    }
    pub fn temperature(&self, i: usize) -> Result<f64> {
        let c = self
            .c
            .get(i)
            .ok_or_else(|| SolverError::InvalidInput("flame cell".into()))?;
        Ok(self.config.unburned_temperature
            + (self.config.burned_temperature - self.config.unburned_temperature) * c)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FlameModelState {
    pub progress: f64,
    pub mixture_fraction: f64,
    pub temperature: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FlameModelResult {
    pub progress_source: f64,
    pub mixture_fraction_source: f64,
    pub heat_release: f64,
}

pub trait FlameModel {
    fn evaluate(&self, state: FlameModelState) -> Result<FlameModelResult>;
}

fn validate_flame_state(s: FlameModelState) -> Result<()> {
    require(
        (0.0..=1.0).contains(&s.progress)
            && (0.0..=1.0).contains(&s.mixture_fraction)
            && s.temperature.is_finite()
            && s.temperature >= 0.0,
        "invalid flame model state",
    )
}

fn nonnegative(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite() && *v >= 0.0)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PremixedFlameModelConfig {
    pub reaction_rate: f64,
    pub heat_release: f64,
    pub ignition_temperature: f64,
}

impl Default for PremixedFlameModelConfig {
    fn default() -> Self {
        Self {
            reaction_rate: 1.0,
            heat_release: 1e6,
            ignition_temperature: 300.0,
        }
    }
}

pub struct PremixedFlameModel(PremixedFlameModelConfig);

impl PremixedFlameModel {
    pub fn new(c: PremixedFlameModelConfig) -> Result<Self> {
        require(
            nonnegative(&[c.reaction_rate, c.heat_release, c.ignition_temperature]),
            "invalid premixed flame controls",
        )?;
        Ok(Self(c))
    }
}

impl FlameModel for PremixedFlameModel {
    fn evaluate(&self, s: FlameModelState) -> Result<FlameModelResult> {
        validate_flame_state(s)?;
        if s.temperature < self.0.ignition_temperature {
            return Ok(FlameModelResult::default());
        }
        let source = self.0.reaction_rate * s.progress * (1.0 - s.progress);
        Ok(FlameModelResult {
            progress_source: source,
            mixture_fraction_source: 0.0,
            heat_release: self.0.heat_release * source,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonPremixedFlameModelConfig {
    pub stoichiometric_mixture_fraction: f64,
    pub mixing_width: f64,
    pub mixing_rate: f64,
    pub heat_release: f64,
    pub ignition_temperature: f64,
}

impl Default for NonPremixedFlameModelConfig {
    fn default() -> Self {
        Self {
            stoichiometric_mixture_fraction: 0.5,
            mixing_width: 0.08,
            mixing_rate: 10.0,
            heat_release: 1e6,
            ignition_temperature: 300.0,
        }
    }
}

/// Mixture-fraction/progress interface for diffusion flames with a Gaussian
/// stoichiometric window, replaceable later by a flamelet table.
pub struct NonPremixedFlameModel(NonPremixedFlameModelConfig);

impl NonPremixedFlameModel {
    pub fn new(c: NonPremixedFlameModelConfig) -> Result<Self> {
        require(
            c.stoichiometric_mixture_fraction > 0.0
                && c.stoichiometric_mixture_fraction < 1.0
                && c.mixing_width.is_finite()
                && c.mixing_width > 0.0
                && nonnegative(&[c.mixing_rate, c.heat_release, c.ignition_temperature]),
            "invalid non-premixed flame controls",
        )?;
        Ok(Self(c))
    }
}

impl FlameModel for NonPremixedFlameModel {
    fn evaluate(&self, s: FlameModelState) -> Result<FlameModelResult> {
        validate_flame_state(s)?;
        let c = self.0;
        if s.temperature < c.ignition_temperature {
            return Ok(FlameModelResult::default());
        }
        let x = (s.mixture_fraction - c.stoichiometric_mixture_fraction) / c.mixing_width;
        let progress_source = c.mixing_rate * (-x * x).exp() * (1.0 - s.progress);
        Ok(FlameModelResult {
            progress_source,
            mixture_fraction_source: -c.mixing_rate
                * (s.mixture_fraction - c.stoichiometric_mixture_fraction),
            heat_release: c.heat_release * progress_source,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SootRadiationConfig {
    pub soot_yield: f64,
    pub oxidation_rate: f64,
    pub absorption_per_mass: f64,
    pub emissivity_floor: f64,
}

impl Default for SootRadiationConfig {
    fn default() -> Self {
        Self {
            soot_yield: 0.02,
            oxidation_rate: 1.0,
            absorption_per_mass: 2e4,
            emissivity_floor: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SootRadiationState {
    pub soot_mass_fraction: f64,
    pub temperature: f64,
    pub heat_release: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SootRadiationResult {
    pub soot_source: f64,
    pub absorption_coefficient: f64,
    pub radiative_source: f64,
}

/// Compact soot production/oxidation and optically thin radiation hook.
pub struct SootRadiationModel(SootRadiationConfig);

impl SootRadiationModel {
    pub fn new(c: SootRadiationConfig) -> Result<Self> {
        require(
            nonnegative(&[
                c.soot_yield,
                c.oxidation_rate,
                c.absorption_per_mass,
                c.emissivity_floor,
            ]),
            "invalid soot-radiation controls",
        )?;
        Ok(Self(c))
    }
    pub fn evaluate(
        &self,
        s: SootRadiationState,
        environment_temperature: f64,
    ) -> Result<SootRadiationResult> {
        require(
            (0.0..=1.0).contains(&s.soot_mass_fraction)
                && s.temperature.is_finite()
                && s.temperature >= 0.0
                && s.heat_release.is_finite()
                && environment_temperature.is_finite()
                && environment_temperature >= 0.0,
            "invalid soot-radiation state",
        )?;
        let source = self.0.soot_yield * s.heat_release.max(0.0)
            - self.0.oxidation_rate * s.soot_mass_fraction;
        let absorption = self.0.absorption_per_mass * s.soot_mass_fraction;
        let radiation =
            optically_thin_radiation_source(s.temperature, environment_temperature, absorption)?;
        Ok(SootRadiationResult {
            soot_source: source,
            absorption_coefficient: absorption,
            radiative_source: radiation,
        })
    }
}
