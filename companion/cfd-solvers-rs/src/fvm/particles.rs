//! Lagrangian particles, pairwise interactions and dense-phase rheology, ported from upstream
//! `cfd/solvers/fvm/{particles, particle_interactions, dense_particle_rheology}`
//! (MIT, cfd_solvers contributors).
use super::PolyMesh;
use crate::{core::Vec3, require, Result};
use std::{collections::HashMap, f64::consts::PI};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Particle {
    pub position: Vec3,
    pub velocity: Vec3,
    pub diameter: f64,
    pub density: f64,
}

impl Default for Particle {
    fn default() -> Self {
        Self {
            position: Vec3::default(),
            velocity: Vec3::default(),
            diameter: 1e-4,
            density: 1000.0,
        }
    }
}

impl Particle {
    pub fn volume(&self) -> f64 {
        PI * self.diameter * self.diameter * self.diameter / 6.0
    }
    pub fn mass(&self) -> f64 {
        self.density * self.volume()
    }
    fn finite(&self) -> bool {
        finite_vec(self.position) && finite_vec(self.velocity)
    }
}

pub(crate) fn finite_vec(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

/// Index of the cell whose centre is nearest (first one on ties), as upstream.
pub(crate) fn nearest_cell(mesh: &PolyMesh, p: Vec3) -> usize {
    let mut best = 0;
    let mut d2 = f64::MAX;
    for (c, cell) in mesh.cells().iter().enumerate() {
        let d = cell.center - p;
        let v = d.dot(d);
        if v < d2 {
            d2 = v;
            best = c;
        }
    }
    best
}

/// One-way or two-way coupled Stokes-drag particle cloud with exact exponential relaxation.
#[derive(Clone, Debug, Default)]
pub struct ParticleCloud {
    particles: Vec<Particle>,
    momentum_source: Vec<Vec3>,
}

impl ParticleCloud {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn add(&mut self, p: Particle) -> Result<()> {
        require(p.diameter > 0.0 && p.density > 0.0, "invalid particle")?;
        self.particles.push(p);
        Ok(())
    }
    pub fn advance(
        &mut self,
        mesh: &PolyMesh,
        fluid_velocity: &[Vec3],
        dynamic_viscosity: f64,
        dt: f64,
        two_way: bool,
    ) -> Result<()> {
        require(
            fluid_velocity.len() == mesh.cell_count() && dynamic_viscosity > 0.0 && dt > 0.0,
            "invalid particle advance",
        )?;
        self.momentum_source = vec![Vec3::default(); mesh.cell_count()];
        for p in &mut self.particles {
            let c = nearest_cell(mesh, p.position);
            let old = p.velocity;
            let tau = p.density * p.diameter * p.diameter / (18.0 * dynamic_viscosity);
            let decay = (-dt / tau).exp();
            p.velocity = fluid_velocity[c] + (p.velocity - fluid_velocity[c]) * decay;
            p.position += p.velocity * dt;
            if two_way {
                let mass = p.density * p.volume();
                self.momentum_source[c] -=
                    (p.velocity - old) * (mass / (dt * mesh.cells()[c].volume));
            }
        }
        Ok(())
    }
    pub fn particles(&self) -> &[Particle] {
        &self.particles
    }
    pub fn fluid_momentum_source(&self) -> &[Vec3] {
        &self.momentum_source
    }
}

/// Saffman shear lift `C rho d^2 sqrt(nu |omega|) (u_rel x omega)/|omega|`; upstream default C = 1.615.
pub fn saffman_lift_force(
    relative_velocity: Vec3,
    vorticity: Vec3,
    fluid_density: f64,
    kinematic_viscosity: f64,
    particle_diameter: f64,
    coefficient: f64,
) -> Result<Vec3> {
    require(
        fluid_density > 0.0
            && kinematic_viscosity > 0.0
            && particle_diameter > 0.0
            && coefficient >= 0.0,
        "invalid Saffman lift state",
    )?;
    let om = vorticity.norm();
    if om == 0.0 {
        return Ok(Vec3::default());
    }
    let scale = coefficient
        * fluid_density
        * particle_diameter
        * particle_diameter
        * (kinematic_viscosity * om).sqrt();
    Ok(relative_velocity.cross(vorticity) * (scale / om))
}

/// Added-mass force on the displaced fluid volume; upstream default coefficient 0.5.
pub fn virtual_mass_force(
    fluid_acceleration: Vec3,
    particle_acceleration: Vec3,
    fluid_density: f64,
    particle_diameter: f64,
    coefficient: f64,
) -> Result<Vec3> {
    require(
        fluid_density > 0.0 && particle_diameter > 0.0 && coefficient >= 0.0,
        "invalid virtual-mass state",
    )?;
    let displaced =
        fluid_density * PI * particle_diameter * particle_diameter * particle_diameter / 6.0;
    Ok((fluid_acceleration - particle_acceleration) * (coefficient * displaced))
}

/// Impulse-based collision along the line of centres for approaching pairs.
pub fn collide_particles_elastic(
    a: &mut Particle,
    b: &mut Particle,
    restitution: f64,
) -> Result<()> {
    require((0.0..=1.0).contains(&restitution), "invalid restitution")?;
    let mut n = b.position - a.position;
    let dist = n.norm();
    if dist == 0.0 {
        return Ok(());
    }
    n = n / dist;
    let rel = (b.velocity - a.velocity).dot(n);
    if rel >= 0.0 {
        return Ok(());
    }
    let (ma, mb) = (a.mass(), b.mass());
    let j = -(1.0 + restitution) * rel / (1.0 / ma + 1.0 / mb);
    a.velocity -= n * (j / ma);
    b.velocity += n * (j / mb);
    Ok(())
}

/// Mass- and volume-conserving merge with mass-weighted position and velocity.
pub fn coalesce_particles(a: &Particle, b: &Particle) -> Particle {
    let (va, vb) = (a.volume(), b.volume());
    let (ma, mb) = (a.density * va, b.density * vb);
    let (m, v) = (ma + mb, va + vb);
    Particle {
        diameter: (6.0 * v / PI).cbrt(),
        density: m / v,
        position: (a.position * ma + b.position * mb) / m,
        velocity: (a.velocity * ma + b.velocity * mb) / m,
    }
}

/// `We = rho u^2 d / sigma > critical`; upstream default critical Weber number 12.
pub fn breakup_by_weber(
    fluid_density: f64,
    relative_speed: f64,
    particle_diameter: f64,
    surface_tension: f64,
    critical_weber: f64,
) -> Result<bool> {
    require(
        fluid_density > 0.0
            && particle_diameter > 0.0
            && surface_tension > 0.0
            && critical_weber > 0.0,
        "invalid breakup state",
    )?;
    Ok(
        fluid_density * relative_speed * relative_speed * particle_diameter / surface_tension
            > critical_weber,
    )
}

/// D-squared law: `d^2(t+dt) = max(0, d^2 - K dt)`.
pub fn evaporate_d2_law(diameter: f64, evaporation_constant: f64, dt: f64) -> Result<f64> {
    require(
        diameter >= 0.0 && evaporation_constant >= 0.0 && dt >= 0.0,
        "invalid evaporation state",
    )?;
    Ok((diameter * diameter - evaporation_constant * dt)
        .max(0.0)
        .sqrt())
}

/// Nearest-cell particle volume fraction, clamped to [0, 1].
pub fn particle_volume_fraction(mesh: &PolyMesh, particles: &[Particle]) -> Vec<f64> {
    let mut a = vec![0.0; mesh.cell_count()];
    for p in particles {
        let c = nearest_cell(mesh, p.position);
        a[c] += p.volume() / mesh.cells()[c].volume;
    }
    for v in &mut a {
        *v = v.clamp(0.0, 1.0);
    }
    a
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleInteractionConfig {
    pub collision_distance_factor: f64,
    pub restitution: f64,
    pub fluid_density: f64,
    pub surface_tension: f64,
    pub coalescence_weber: f64,
    pub breakup_weber: f64,
    /// Retained for upstream configuration compatibility; upstream does not use it either.
    pub breakup_velocity_fraction: f64,
    pub enable_coalescence: bool,
    pub enable_breakup: bool,
    pub maximum_particles: usize,
}

impl Default for ParticleInteractionConfig {
    fn default() -> Self {
        Self {
            collision_distance_factor: 1.0,
            restitution: 0.85,
            fluid_density: 1.2,
            surface_tension: 0.072,
            coalescence_weber: 6.0,
            breakup_weber: 12.0,
            breakup_velocity_fraction: 0.15,
            enable_coalescence: true,
            enable_breakup: true,
            maximum_particles: 1_000_000,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ParticleInteractionStats {
    pub candidate_pairs: usize,
    pub collisions: usize,
    pub coalescences: usize,
    pub breakups: usize,
    pub removed_particles: usize,
    pub momentum_before: f64,
    pub momentum_after: f64,
    pub kinetic_energy_before: f64,
    pub kinetic_energy_after: f64,
}

/// Deterministic cell-linked collision, coalescence and Weber breakup.
#[derive(Clone, Debug)]
pub struct ParticleInteractionSystem {
    config: ParticleInteractionConfig,
    stats: ParticleInteractionStats,
}

fn momentum_magnitude(p: &[Particle]) -> f64 {
    let mut m = Vec3::default();
    for x in p {
        m += x.velocity * x.mass();
    }
    m.norm()
}

fn kinetic_energy(p: &[Particle]) -> f64 {
    p.iter()
        .map(|x| {
            let s = x.velocity.norm();
            0.5 * x.mass() * s * s
        })
        .sum()
}

fn unit_or_x(v: Vec3) -> Vec3 {
    let l = v.norm();
    if l > 1e-14 {
        v / l
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    }
}

impl ParticleInteractionSystem {
    pub fn new(config: ParticleInteractionConfig) -> Result<Self> {
        let c = config;
        let finite = [
            c.collision_distance_factor,
            c.restitution,
            c.fluid_density,
            c.surface_tension,
            c.coalescence_weber,
            c.breakup_weber,
            c.breakup_velocity_fraction,
        ]
        .iter()
        .all(|v| v.is_finite());
        require(
            finite
                && c.collision_distance_factor > 0.0
                && (0.0..=1.0).contains(&c.restitution)
                && c.fluid_density > 0.0
                && c.surface_tension > 0.0
                && c.coalescence_weber >= 0.0
                && c.breakup_weber > c.coalescence_weber
                && c.breakup_velocity_fraction >= 0.0
                && c.maximum_particles > 0,
            "invalid particle interaction controls",
        )?;
        Ok(Self {
            config,
            stats: ParticleInteractionStats::default(),
        })
    }
    pub fn config(&self) -> &ParticleInteractionConfig {
        &self.config
    }
    pub fn last_stats(&self) -> &ParticleInteractionStats {
        &self.stats
    }

    /// Resolves each overlapping, approaching pair once, in sorted pair order.
    pub fn resolve(&mut self, particles: &mut Vec<Particle>) -> Result<()> {
        let cfg = self.config;
        require(
            particles.len() <= cfg.maximum_particles,
            "particle population exceeds configured maximum",
        )?;
        require(
            particles.iter().all(|p| {
                p.finite()
                    && p.diameter > 0.0
                    && p.diameter.is_finite()
                    && p.density > 0.0
                    && p.density.is_finite()
            }),
            "invalid particle state",
        )?;
        let mut stats = ParticleInteractionStats {
            momentum_before: momentum_magnitude(particles),
            kinetic_energy_before: kinetic_energy(particles),
            ..Default::default()
        };
        if particles.len() < 2 {
            stats.momentum_after = stats.momentum_before;
            stats.kinetic_energy_after = stats.kinetic_energy_before;
            self.stats = stats;
            return Ok(());
        }
        let largest = particles.iter().fold(0.0_f64, |m, p| m.max(p.diameter));
        let cell = (largest * cfg.collision_distance_factor).max(1e-14);
        let key = |p: Vec3| {
            (
                (p.x / cell).floor() as i64,
                (p.y / cell).floor() as i64,
                (p.z / cell).floor() as i64,
            )
        };
        let mut buckets: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        for (i, p) in particles.iter().enumerate() {
            buckets.entry(key(p.position)).or_default().push(i);
        }
        let mut pairs = Vec::new();
        for (&(kx, ky, kz), members) in &buckets {
            for dz in -1..=1 {
                for dy in -1..=1 {
                    for dx in -1..=1 {
                        let Some(found) = buckets.get(&(
                            kx.saturating_add(dx),
                            ky.saturating_add(dy),
                            kz.saturating_add(dz),
                        )) else {
                            continue;
                        };
                        for &i in members {
                            for &j in found {
                                if i < j {
                                    pairs.push((i, j));
                                }
                            }
                        }
                    }
                }
            }
        }
        pairs.sort_unstable();
        pairs.dedup();

        let n = particles.len();
        let (mut claimed, mut removed, mut breakup) =
            (vec![false; n], vec![false; n], vec![false; n]);
        for &(i, j) in &pairs {
            stats.candidate_pairs += 1;
            if claimed[i] || claimed[j] || removed[i] || removed[j] {
                continue;
            }
            let delta = particles[j].position - particles[i].position;
            let distance = delta.norm();
            let collision_distance = 0.5
                * (particles[i].diameter + particles[j].diameter)
                * cfg.collision_distance_factor;
            if distance > collision_distance {
                continue;
            }
            let normal = unit_or_x(delta);
            let approach = (particles[j].velocity - particles[i].velocity).dot(normal);
            if approach >= 0.0 || approach.is_nan() {
                continue;
            }
            stats.collisions += 1;
            let mean_d = 0.5 * (particles[i].diameter + particles[j].diameter);
            let weber = cfg.fluid_density * approach * approach * mean_d / cfg.surface_tension;
            if cfg.enable_coalescence && weber <= cfg.coalescence_weber {
                particles[i] = coalesce_particles(&particles[i], &particles[j]);
                claimed[i] = true;
                claimed[j] = true;
                removed[j] = true;
                stats.coalescences += 1;
                continue;
            }
            let (head, tail) = particles.split_at_mut(j);
            collide_particles_elastic(&mut head[i], &mut tail[0], cfg.restitution)?;
            claimed[i] = true;
            claimed[j] = true;
            if cfg.enable_breakup && weber >= cfg.breakup_weber {
                let target = if particles[i].diameter >= particles[j].diameter {
                    i
                } else {
                    j
                };
                breakup[target] = true;
            }
        }

        let mut result = Vec::with_capacity(n);
        for (i, parent) in particles.iter().enumerate() {
            if removed[i] {
                stats.removed_particles += 1;
                continue;
            }
            if !breakup[i] {
                result.push(*parent);
                continue;
            }
            require(
                result.len() + 2 <= cfg.maximum_particles,
                "particle breakup exceeds configured maximum",
            )?;
            let direction = unit_or_x(parent.velocity);
            let mut child = *parent;
            child.diameter = parent.diameter * 0.5_f64.cbrt();
            let offset = 0.26 * parent.diameter;
            child.position = parent.position - direction * offset;
            result.push(child);
            child.position = parent.position + direction * offset;
            result.push(child);
            stats.breakups += 1;
        }
        *particles = result;
        stats.momentum_after = momentum_magnitude(particles);
        stats.kinetic_energy_after = kinetic_energy(particles);
        self.stats = stats;
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DenseRheologyConfig {
    pub fluid_dynamic_viscosity: f64,
    pub particle_density: f64,
    pub grain_diameter: f64,
    pub maximum_solid_fraction: f64,
    pub static_friction: f64,
    pub dynamic_friction: f64,
    pub inertial_scale: f64,
    pub regularization: f64,
}

impl Default for DenseRheologyConfig {
    fn default() -> Self {
        Self {
            fluid_dynamic_viscosity: 1.8e-5,
            particle_density: 1000.0,
            grain_diameter: 1e-4,
            maximum_solid_fraction: 0.64,
            static_friction: 0.38,
            dynamic_friction: 0.64,
            inertial_scale: 0.30,
            regularization: 1e-9,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DenseRheologyResult {
    pub solid_fraction: Vec<f64>,
    pub granular_pressure: Vec<f64>,
    pub effective_viscosity: Vec<f64>,
    pub carrier_momentum_source: Vec<Vec3>,
    pub particle_momentum_source: Vec<Vec3>,
}

/// Hindered (Richardson-Zaki/Ergun-like) drag with a regularised mu(I) granular stress;
/// carrier and particle sources are equal and opposite.
#[derive(Clone, Debug)]
pub struct DenseParticleRheology {
    mesh: PolyMesh,
    config: DenseRheologyConfig,
}

impl DenseParticleRheology {
    pub fn new(mesh: PolyMesh, config: DenseRheologyConfig) -> Result<Self> {
        let c = config;
        let finite = [
            c.fluid_dynamic_viscosity,
            c.particle_density,
            c.grain_diameter,
            c.maximum_solid_fraction,
            c.static_friction,
            c.dynamic_friction,
            c.inertial_scale,
            c.regularization,
        ]
        .iter()
        .all(|v| v.is_finite());
        require(
            finite
                && c.fluid_dynamic_viscosity > 0.0
                && c.particle_density > 0.0
                && c.grain_diameter > 0.0
                && c.maximum_solid_fraction > 0.0
                && c.maximum_solid_fraction < 1.0
                && c.static_friction >= 0.0
                && c.dynamic_friction >= c.static_friction
                && c.inertial_scale > 0.0
                && c.regularization > 0.0,
            "invalid dense rheology controls",
        )?;
        Ok(Self { mesh, config })
    }
    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn config(&self) -> &DenseRheologyConfig {
        &self.config
    }
    pub fn evaluate(
        &self,
        particles: &[Particle],
        carrier_velocity: &[Vec3],
    ) -> Result<DenseRheologyResult> {
        let cells = self.mesh.cell_count();
        require(
            carrier_velocity.len() == cells,
            "dense rheology carrier field size mismatch",
        )?;
        let c = self.config;
        let mut r = DenseRheologyResult {
            solid_fraction: vec![0.0; cells],
            granular_pressure: vec![0.0; cells],
            effective_viscosity: vec![0.0; cells],
            carrier_momentum_source: vec![Vec3::default(); cells],
            particle_momentum_source: vec![Vec3::default(); cells],
        };
        let mut mean_velocity = vec![Vec3::default(); cells];
        let mut weighted_volume = vec![0.0; cells];
        for p in particles {
            require(
                p.diameter > 0.0 && p.density > 0.0 && p.finite(),
                "invalid dense particle state",
            )?;
            let cell = nearest_cell(&self.mesh, p.position);
            let w = p.volume() / self.mesh.cells()[cell].volume;
            r.solid_fraction[cell] += w;
            mean_velocity[cell] += p.velocity * w;
            weighted_volume[cell] += w;
        }
        for cell in 0..cells {
            let raw = r.solid_fraction[cell];
            if raw <= 0.0 {
                continue;
            }
            let phi = raw.min(c.maximum_solid_fraction * (1.0 - 1e-10));
            r.solid_fraction[cell] = phi;
            let mean = mean_velocity[cell] / weighted_volume[cell].max(c.regularization);
            let slip = mean - carrier_velocity[cell];
            let shear_rate = slip.norm() / c.grain_diameter.max(c.regularization);
            let gap = (c.maximum_solid_fraction - phi).max(c.regularization);
            let pressure = c.particle_density
                * c.grain_diameter
                * c.grain_diameter
                * shear_rate
                * shear_rate
                * phi
                / gap;
            let inertial = c.grain_diameter * shear_rate
                / (pressure / c.particle_density).max(c.regularization).sqrt();
            let friction = c.static_friction
                + (c.dynamic_friction - c.static_friction) * c.inertial_scale
                    / (c.inertial_scale + inertial);
            let viscosity = pressure * friction / shear_rate.max(c.regularization);
            let hindered = (1.0 - phi / c.maximum_solid_fraction).max(1e-6).powf(-2.0);
            let beta = 18.0 * c.fluid_dynamic_viscosity * phi * (1.0 + 2.65 * phi) * hindered
                / (c.grain_diameter * c.grain_diameter);
            let force = slip * beta;
            r.granular_pressure[cell] = pressure;
            r.effective_viscosity[cell] = viscosity;
            r.particle_momentum_source[cell] = force;
            r.carrier_momentum_source[cell] = force * -1.0;
        }
        Ok(r)
    }
}
