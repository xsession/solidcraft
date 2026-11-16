//! Tabulated cross sections, a coupled neutron/photon depletion kernel and fixed-source 1-D slab
//! Monte Carlo transport. Ported from upstream `cfd/particle/{cross_sections, monte_carlo}`
//! (MIT, cfd_solvers contributors).
use crate::core::Mt19937_64;
use crate::{require, Result};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ParticleKind {
    #[default]
    Neutron,
    Photon,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CrossSectionPoint {
    pub energy: f64,
    pub area: f64,
}

/// Piecewise-linear cross section in energy, clamped outside the table.
#[derive(Clone, Debug, PartialEq)]
pub struct CrossSectionTable {
    points: Vec<CrossSectionPoint>,
}

impl CrossSectionTable {
    pub fn new(points: Vec<CrossSectionPoint>) -> Result<Self> {
        require(!points.is_empty(), "cross-section table is empty")?;
        for (i, p) in points.iter().enumerate() {
            require(
                p.energy > 0.0
                    && p.area >= 0.0
                    && p.energy.is_finite()
                    && p.area.is_finite()
                    && (i == 0 || p.energy > points[i - 1].energy),
                "invalid cross-section point",
            )?;
        }
        Ok(Self { points })
    }
    pub fn value(&self, energy: f64) -> Result<f64> {
        require(
            energy > 0.0 && energy.is_finite(),
            "invalid particle energy",
        )?;
        let p = &self.points;
        let (first, last) = (p[0], p[p.len() - 1]);
        if p.len() == 1 || energy <= first.energy {
            return Ok(first.area);
        }
        if energy >= last.energy {
            return Ok(last.area);
        }
        for w in p.windows(2) {
            if energy <= w[1].energy {
                let weight = (energy - w[0].energy) / (w[1].energy - w[0].energy);
                return Ok(w[0].area + weight * (w[1].area - w[0].area));
            }
        }
        Ok(last.area)
    }
    pub fn points(&self) -> &[CrossSectionPoint] {
        &self.points
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DepletionReaction {
    pub parent: usize,
    pub daughter: usize,
    pub kind: ParticleKind,
    pub cross_section: CrossSectionTable,
    pub branching_fraction: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DepletionDiagnostics {
    pub total_inventory: f64,
    pub absorbed_rate: f64,
    pub photon_rate: f64,
    pub neutron_rate: f64,
}

/// Explicit depletion: every reaction uses the old inventory and all changes apply together.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiSpeciesDepletion {
    inventory: Vec<f64>,
    reactions: Vec<DepletionReaction>,
}

impl MultiSpeciesDepletion {
    pub fn new(inventory: Vec<f64>) -> Result<Self> {
        require(!inventory.is_empty(), "depletion inventory is empty")?;
        require(
            inventory.iter().all(|v| *v >= 0.0 && v.is_finite()),
            "invalid depletion inventory",
        )?;
        Ok(Self {
            inventory,
            reactions: Vec::new(),
        })
    }
    pub fn add_reaction(&mut self, reaction: DepletionReaction) -> Result<()> {
        require(
            reaction.parent < self.inventory.len()
                && reaction.daughter < self.inventory.len()
                && reaction.branching_fraction >= 0.0
                && reaction.branching_fraction <= 1.0,
            "invalid depletion reaction",
        )?;
        self.reactions.push(reaction);
        Ok(())
    }
    pub fn advance(
        &mut self,
        dt: f64,
        neutron_flux: f64,
        photon_flux: f64,
        energy: f64,
    ) -> Result<DepletionDiagnostics> {
        require(
            dt > 0.0
                && neutron_flux >= 0.0
                && photon_flux >= 0.0
                && dt.is_finite()
                && neutron_flux.is_finite()
                && photon_flux.is_finite(),
            "invalid depletion step",
        )?;
        let mut delta = vec![0.0; self.inventory.len()];
        let mut d = DepletionDiagnostics::default();
        for r in &self.reactions {
            let flux = if r.kind == ParticleKind::Neutron {
                neutron_flux
            } else {
                photon_flux
            };
            let rate = self.inventory[r.parent] * r.cross_section.value(energy)? * flux;
            let amount = self.inventory[r.parent].min(dt * rate);
            delta[r.parent] -= amount;
            delta[r.daughter] += amount * r.branching_fraction;
            d.absorbed_rate += rate;
            if r.kind == ParticleKind::Photon {
                d.photon_rate += rate;
            } else {
                d.neutron_rate += rate;
            }
        }
        for (v, dv) in self.inventory.iter_mut().zip(&delta) {
            *v = (*v + dv).max(0.0);
        }
        d.total_inventory = self.inventory.iter().sum();
        Ok(d)
    }
    pub fn inventory(&self) -> &[f64] {
        &self.inventory
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabTransportConfig {
    pub width: f64,
    pub total_cross_section: f64,
    pub absorption_cross_section: f64,
}

impl Default for SlabTransportConfig {
    fn default() -> Self {
        Self {
            width: 1.0,
            total_cross_section: 1.0,
            absorption_cross_section: 0.2,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SlabTransportTally {
    pub histories: usize,
    pub absorbed: usize,
    pub transmitted: usize,
    pub reflected: usize,
    pub total_path: f64,
}

impl SlabTransportTally {
    pub fn absorption_probability(&self) -> f64 {
        if self.histories == 0 {
            0.0
        } else {
            self.absorbed as f64 / self.histories as f64
        }
    }
}

/// Rod-model (±x) slab transport with isotropic-in-direction scattering, using the exact
/// `std::mt19937_64` + `std::uniform_real_distribution` sequence of upstream.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabTransport1D {
    config: SlabTransportConfig,
}

impl SlabTransport1D {
    pub fn new(config: SlabTransportConfig) -> Result<Self> {
        require(
            config.width > 0.0
                && config.total_cross_section > 0.0
                && config.absorption_cross_section >= 0.0
                && config.absorption_cross_section <= config.total_cross_section,
            "invalid slab transport configuration",
        )?;
        Ok(Self { config })
    }
    pub fn run(&self, histories: usize, seed: u64) -> Result<SlabTransportTally> {
        require(histories > 0, "Monte Carlo requires histories")?;
        let c = &self.config;
        let mut rng = Mt19937_64::new(seed);
        let mut tally = SlabTransportTally {
            histories,
            ..Default::default()
        };
        let pa = c.absorption_cross_section / c.total_cross_section;
        for _ in 0..histories {
            let (mut x, mut dir) = (0.0f64, 1.0f64);
            loop {
                let sample = rng.uniform(0.0, 1.0).max(1.0e-15);
                let free_path = -sample.ln() / c.total_cross_section;
                let boundary = if dir > 0.0 { c.width - x } else { x };
                if free_path >= boundary {
                    tally.total_path += boundary;
                    if dir > 0.0 {
                        tally.transmitted += 1;
                    } else {
                        tally.reflected += 1;
                    }
                    break;
                }
                x += dir * free_path;
                tally.total_path += free_path;
                if rng.uniform(0.0, 1.0) < pa {
                    tally.absorbed += 1;
                    break;
                }
                dir = if rng.uniform(0.0, 1.0) < 0.5 {
                    -1.0
                } else {
                    1.0
                };
            }
        }
        Ok(tally)
    }
}
