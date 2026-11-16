//! Periodic 1-D molecular dynamics: harmonic chain (velocity Verlet), neighbor list,
//! Lennard-Jones forces, velocity-rescaling thermostat and isotropic barostat, ported from
//! upstream `cfd/atomistics/{md1d, neighbor_list}` (MIT, cfd_solvers contributors).
use crate::{require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicChainConfig {
    pub particles: usize,
    pub length: f64,
    pub mass: f64,
    pub spring_constant: f64,
}

impl Default for HarmonicChainConfig {
    fn default() -> Self {
        Self {
            particles: 32,
            length: 1.0,
            mass: 1.0,
            spring_constant: 100.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicChainEnergy {
    pub kinetic: f64,
    pub potential: f64,
}

impl HarmonicChainEnergy {
    pub fn total(&self) -> f64 {
        self.kinetic + self.potential
    }
}

/// Periodic nearest-neighbour spring chain on `[0, L)` with equilibrium spacing `L/N`.
#[derive(Clone, Debug)]
pub struct HarmonicChain1D {
    config: HarmonicChainConfig,
    positions: Vec<f64>,
    velocities: Vec<f64>,
    forces: Vec<f64>,
}

impl HarmonicChain1D {
    pub fn new(config: HarmonicChainConfig) -> Result<Self> {
        require(
            config.particles >= 3
                && config.length > 0.0
                && config.mass > 0.0
                && config.spring_constant > 0.0,
            "invalid atomistic chain configuration",
        )?;
        let spacing = config.length / config.particles as f64;
        let mut s = Self {
            config,
            positions: (0..config.particles).map(|i| spacing * i as f64).collect(),
            velocities: vec![0.0; config.particles],
            forces: vec![0.0; config.particles],
        };
        s.compute_forces();
        Ok(s)
    }
    pub fn set_velocity(&mut self, particle: usize, value: f64) -> Result<()> {
        require(
            particle < self.velocities.len() && value.is_finite(),
            "atomistic particle index",
        )?;
        self.velocities[particle] = value;
        Ok(())
    }
    pub fn positions(&self) -> &[f64] {
        &self.positions
    }
    pub fn velocities(&self) -> &[f64] {
        &self.velocities
    }
    fn distance(&self, i: usize) -> (usize, f64) {
        let next = if i + 1 == self.positions.len() {
            0
        } else {
            i + 1
        };
        let mut d = self.positions[next] - self.positions[i];
        if d < 0.0 {
            d += self.config.length;
        }
        (next, d)
    }
    fn compute_forces(&mut self) {
        self.forces.fill(0.0);
        let spacing = self.config.length / self.config.particles as f64;
        for i in 0..self.positions.len() {
            let (next, d) = self.distance(i);
            let f = self.config.spring_constant * (d - spacing);
            self.forces[i] += f;
            self.forces[next] -= f;
        }
    }
    /// Velocity-Verlet step with periodic wrapping.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(
            dt > 0.0 && dt.is_finite(),
            "invalid molecular dynamics timestep",
        )?;
        let inverse_mass = 1.0 / self.config.mass;
        let length = self.config.length;
        for i in 0..self.positions.len() {
            self.velocities[i] += 0.5 * dt * self.forces[i] * inverse_mass;
            let mut x = (self.positions[i] + dt * self.velocities[i]) % length;
            if x < 0.0 {
                x += length;
            }
            self.positions[i] = x;
        }
        self.compute_forces();
        for i in 0..self.positions.len() {
            self.velocities[i] += 0.5 * dt * self.forces[i] * inverse_mass;
        }
        Ok(())
    }
    pub fn energy(&self) -> HarmonicChainEnergy {
        let spacing = self.config.length / self.config.particles as f64;
        let kinetic = self
            .velocities
            .iter()
            .map(|v| 0.5 * self.config.mass * v * v)
            .sum();
        let potential = (0..self.positions.len())
            .map(|i| {
                let e = self.distance(i).1 - spacing;
                0.5 * self.config.spring_constant * e * e
            })
            .sum();
        HarmonicChainEnergy { kinetic, potential }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NeighborList1DConfig {
    pub periodic_length: f64,
    pub cutoff: f64,
    pub skin: f64,
}

impl Default for NeighborList1DConfig {
    fn default() -> Self {
        Self {
            periodic_length: 1.0,
            cutoff: 1.0,
            skin: 0.1,
        }
    }
}

/// All-pairs Verlet list with skin under the minimum-image convention.
#[derive(Clone, Debug)]
pub struct PeriodicNeighborList1D {
    config: NeighborList1DConfig,
    pairs: Vec<(usize, usize)>,
}

impl PeriodicNeighborList1D {
    pub fn new(config: NeighborList1DConfig) -> Result<Self> {
        let c = config;
        require(
            c.periodic_length > 0.0
                && c.cutoff > 0.0
                && c.skin >= 0.0
                && c.periodic_length.is_finite()
                && c.cutoff.is_finite()
                && c.skin.is_finite()
                && c.cutoff + c.skin < 0.5 * c.periodic_length,
            "invalid periodic neighbor-list configuration",
        )?;
        Ok(Self {
            config,
            pairs: Vec::new(),
        })
    }
    pub fn rebuild(&mut self, positions: &[f64]) -> Result<()> {
        require(
            positions.iter().all(|p| p.is_finite()),
            "non-finite atomistic position",
        )?;
        self.pairs.clear();
        let reach = self.config.cutoff + self.config.skin;
        for a in 0..positions.len() {
            for b in a + 1..positions.len() {
                if self.minimum_image(positions[b] - positions[a]).abs() <= reach {
                    self.pairs.push((a, b));
                }
            }
        }
        Ok(())
    }
    pub fn pairs(&self) -> &[(usize, usize)] {
        &self.pairs
    }
    pub fn minimum_image(&self, mut distance: f64) -> f64 {
        let half = 0.5 * self.config.periodic_length;
        while distance > half {
            distance -= self.config.periodic_length;
        }
        while distance < -half {
            distance += self.config.periodic_length;
        }
        distance
    }
    pub fn cutoff(&self) -> f64 {
        self.config.cutoff
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LennardJonesParameters {
    pub epsilon: f64,
    pub sigma: f64,
}

impl Default for LennardJonesParameters {
    fn default() -> Self {
        Self {
            epsilon: 1.0,
            sigma: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ForceFieldDiagnostics {
    pub potential_energy: f64,
    pub virial: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct LennardJones1DForceField {
    p: LennardJonesParameters,
}

impl LennardJones1DForceField {
    pub fn new(parameters: LennardJonesParameters) -> Result<Self> {
        let p = parameters;
        require(
            p.epsilon > 0.0 && p.sigma > 0.0 && p.epsilon.is_finite() && p.sigma.is_finite(),
            "invalid Lennard-Jones parameters",
        )?;
        Ok(Self { p })
    }
    /// Truncated (unshifted) 12-6 forces over the listed pairs within the cutoff.
    pub fn compute(
        &self,
        positions: &[f64],
        neighbors: &PeriodicNeighborList1D,
        forces: &mut [f64],
    ) -> Result<ForceFieldDiagnostics> {
        require(forces.len() == positions.len(), "force-field size mismatch")?;
        require(
            neighbors
                .pairs()
                .iter()
                .all(|&(a, b)| a < positions.len() && b < positions.len()),
            "neighbor list does not match positions",
        )?;
        forces.fill(0.0);
        let mut d = ForceFieldDiagnostics::default();
        for &(a, b) in neighbors.pairs() {
            let distance = neighbors.minimum_image(positions[b] - positions[a]);
            let r = distance.abs();
            if r.is_nan() || r <= 0.0 || r > neighbors.cutoff() {
                continue;
            }
            let inv = self.p.sigma / r;
            let inv2 = inv * inv;
            let inv6 = inv2 * inv2 * inv2;
            let inv12 = inv6 * inv6;
            d.potential_energy += 4.0 * self.p.epsilon * (inv12 - inv6);
            let radial = 24.0 * self.p.epsilon * (2.0 * inv12 - inv6) / r;
            let sign = if distance >= 0.0 { 1.0 } else { -1.0 };
            forces[a] -= sign * radial;
            forces[b] += sign * radial;
            d.virial += radial * r;
        }
        Ok(d)
    }
}

/// Berendsen-type velocity rescaling toward `k_B T = 2 KE / N` (1-D, unit Boltzmann constant).
pub fn velocity_rescaling_thermostat(
    velocities: &mut [f64],
    mass: f64,
    target_temperature: f64,
    dt: f64,
    relaxation_time: f64,
) -> Result<()> {
    require(
        mass > 0.0
            && target_temperature >= 0.0
            && dt > 0.0
            && relaxation_time > 0.0
            && mass.is_finite()
            && target_temperature.is_finite()
            && dt.is_finite()
            && relaxation_time.is_finite(),
        "invalid thermostat controls",
    )?;
    if velocities.is_empty() {
        return Ok(());
    }
    let kinetic: f64 = velocities.iter().map(|v| 0.5 * mass * v * v).sum();
    let temperature = 2.0 * kinetic / velocities.len() as f64;
    if temperature.is_nan() || temperature <= 0.0 {
        return Ok(());
    }
    let factor = (1.0 + dt / relaxation_time * (target_temperature / temperature - 1.0))
        .max(0.0)
        .sqrt();
    for v in velocities {
        *v *= factor;
    }
    Ok(())
}

/// Isotropic box/coordinate scaling by `clamp(1 - beta dP dt, 0.5, 1.5)`; returns the factor.
pub fn isotropic_barostat_scale(
    positions: &mut [f64],
    periodic_length: &mut f64,
    pressure_error: f64,
    compressibility: f64,
    dt: f64,
) -> Result<f64> {
    require(
        *periodic_length > 0.0
            && compressibility >= 0.0
            && dt > 0.0
            && periodic_length.is_finite()
            && pressure_error.is_finite()
            && compressibility.is_finite()
            && dt.is_finite(),
        "invalid barostat controls",
    )?;
    let factor = (1.0 - compressibility * pressure_error * dt).clamp(0.5, 1.5);
    for p in positions {
        *p *= factor;
    }
    *periodic_length *= factor;
    Ok(factor)
}
