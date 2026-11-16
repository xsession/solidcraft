//! Charged-particle pushers (Boris, relativistic Boris, Vay), null-collision Monte Carlo
//! collisions, secondary-electron yield and axis-aligned wall interaction. Ported from upstream
//! `cfd/particle/electromagnetic` (MIT, cfd_solvers contributors).
use crate::core::Mt19937_64;
use crate::{require, Result, SolverError};
use std::{f64::consts::PI, sync::Arc};

pub(crate) const EPSILON0: f64 = 8.8541878128e-12;
pub(crate) const MU0: f64 = 1.25663706212e-6;
pub(crate) const SPEED_OF_LIGHT: f64 = 299_792_458.0;
pub(crate) const ELEMENTARY_CHARGE: f64 = 1.602176634e-19;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

pub(crate) fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}
pub(crate) fn scale(a: Vec3, s: f64) -> Vec3 {
    Vec3::new(a.x * s, a.y * s, a.z * s)
}
pub(crate) fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}
pub(crate) fn norm2(a: Vec3) -> f64 {
    a.x * a.x + a.y * a.y + a.z * a.z
}
pub(crate) fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ElectromagneticField {
    pub electric_v_per_m: Vec3,
    pub magnetic_t: Vec3,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChargedParticle {
    pub position_m: Vec3,
    pub velocity_m_per_s: Vec3,
    pub charge_c: f64,
    pub mass_kg: f64,
    pub weight: f64,
}

impl Default for ChargedParticle {
    fn default() -> Self {
        Self {
            position_m: Vec3::default(),
            velocity_m_per_s: Vec3::default(),
            charge_c: 1.0,
            mass_kg: 1.0,
            weight: 1.0,
        }
    }
}

/// Field at (position, time).
pub type FieldSampler = Arc<dyn Fn(Vec3, f64) -> Result<ElectromagneticField> + Send + Sync>;

fn check_push(p: &ChargedParticle, dt: f64, what: &str) -> Result<()> {
    require(
        p.mass_kg > 0.0
            && p.mass_kg.is_finite()
            && p.charge_c.is_finite()
            && dt > 0.0
            && dt.is_finite(),
        what,
    )
}

/// Non-relativistic Boris push (speed-preserving in a pure magnetic field).
pub fn boris_push(p: &mut ChargedParticle, field: &ElectromagneticField, dt: f64) -> Result<()> {
    check_push(p, dt, "invalid Boris particle parameters")?;
    let h = 0.5 * p.charge_c * dt / p.mass_kg;
    let v_minus = add(p.velocity_m_per_s, scale(field.electric_v_per_m, h));
    let t = scale(field.magnetic_t, h);
    let s = scale(t, 2.0 / (1.0 + norm2(t)));
    let v_prime = add(v_minus, cross(v_minus, t));
    let v_plus = add(v_minus, cross(v_prime, s));
    p.velocity_m_per_s = add(v_plus, scale(field.electric_v_per_m, h));
    p.position_m = add(p.position_m, scale(p.velocity_m_per_s, dt));
    Ok(())
}

/// Boris steps with the field sampled at mid-step time.
pub fn track_particle(
    p: &mut ChargedParticle,
    field: &dyn Fn(Vec3, f64) -> Result<ElectromagneticField>,
    start_time: f64,
    dt: f64,
    steps: usize,
) -> Result<()> {
    let mut time = start_time;
    for _ in 0..steps {
        let f = field(p.position_m, time + 0.5 * dt)?;
        boris_push(p, &f, dt)?;
        time += dt;
    }
    Ok(())
}

pub fn lorentz_gamma(p: &ChargedParticle) -> Result<f64> {
    let v2 = norm2(p.velocity_m_per_s);
    let c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    require(v2 < c2, "particle speed must be below c")?;
    Ok(1.0 / (1.0 - v2 / c2).sqrt())
}

pub fn kinetic_energy_j(p: &ChargedParticle) -> Result<f64> {
    require(p.mass_kg > 0.0, "particle mass must be positive")?;
    Ok((lorentz_gamma(p)? - 1.0) * p.mass_kg * SPEED_OF_LIGHT * SPEED_OF_LIGHT)
}

/// Relativistic Boris push on proper velocity `u = γv`.
pub fn relativistic_boris_push(
    p: &mut ChargedParticle,
    field: &ElectromagneticField,
    dt: f64,
) -> Result<()> {
    check_push(p, dt, "invalid relativistic Boris parameters")?;
    let c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    let gamma0 = lorentz_gamma(p)?;
    let u = scale(p.velocity_m_per_s, gamma0);
    let h = 0.5 * p.charge_c * dt / p.mass_kg;
    let u_minus = add(u, scale(field.electric_v_per_m, h));
    let gamma_minus = (1.0 + norm2(u_minus) / c2).sqrt();
    let t = scale(field.magnetic_t, h / gamma_minus);
    let s = scale(t, 2.0 / (1.0 + norm2(t)));
    let u_prime = add(u_minus, cross(u_minus, t));
    let u_plus = add(u_minus, cross(u_prime, s));
    let u = add(u_plus, scale(field.electric_v_per_m, h));
    let gamma_new = (1.0 + norm2(u) / c2).sqrt();
    p.velocity_m_per_s = scale(u, 1.0 / gamma_new);
    p.position_m = add(p.position_m, scale(p.velocity_m_per_s, dt));
    Ok(())
}

/// Vay (2008) relativistic push, free of the spurious E×B drift at high γ.
pub fn vay_push(p: &mut ChargedParticle, field: &ElectromagneticField, dt: f64) -> Result<()> {
    check_push(p, dt, "invalid Vay particle parameters")?;
    let c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    let gamma0 = lorentz_gamma(p)?;
    let u0 = scale(p.velocity_m_per_s, gamma0);
    let h = 0.5 * p.charge_c * dt / p.mass_kg;
    let tau = scale(field.magnetic_t, h);
    let u_minus = add(u0, scale(field.electric_v_per_m, h));
    let v0 = p.velocity_m_per_s;
    let u_prime = add(
        add(u_minus, scale(field.electric_v_per_m, h)),
        scale(cross(v0, field.magnetic_t), h),
    );
    let gamma_prime = (1.0 + norm2(u_prime) / c2).sqrt();
    let tau2 = norm2(tau);
    let sigma = gamma_prime * gamma_prime - tau2;
    let u_star = dot(u_prime, tau) / SPEED_OF_LIGHT;
    let gamma_new =
        (0.5 * (sigma + (sigma * sigma + 4.0 * (tau2 + u_star * u_star)).sqrt())).sqrt();
    let t = scale(tau, 1.0 / gamma_new);
    let t2 = norm2(t);
    let numerator = add(add(u_prime, cross(u_prime, t)), scale(t, dot(u_prime, t)));
    let u_new = scale(numerator, 1.0 / (1.0 + t2));
    let gamma_final = (1.0 + norm2(u_new) / c2).sqrt();
    p.velocity_m_per_s = scale(u_new, 1.0 / gamma_final);
    p.position_m = add(p.position_m, scale(p.velocity_m_per_s, dt));
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NeutralCollisionModel {
    pub neutral_density_per_m3: f64,
    pub elastic_cross_section_m2: f64,
    pub ionization_cross_section_m2: f64,
    pub ionization_energy_ev: f64,
    pub random_seed: u64,
}

impl Default for NeutralCollisionModel {
    fn default() -> Self {
        Self {
            neutral_density_per_m3: 0.0,
            elastic_cross_section_m2: 0.0,
            ionization_cross_section_m2: 0.0,
            ionization_energy_ev: 15.0,
            random_seed: 0xC0FFEE,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CollisionStatistics {
    pub elastic_events: usize,
    pub ionization_events: usize,
    pub energy_loss_j: f64,
}

impl CollisionStatistics {
    pub(crate) fn accumulate(&mut self, s: &CollisionStatistics) {
        self.elastic_events += s.elastic_events;
        self.ionization_events += s.ionization_events;
        self.energy_loss_j += s.energy_loss_j;
    }
}

/// Generic velocity view used by the collision kernel (3-D velocity, mass, weight).
pub(crate) trait CollisionParticle: Clone {
    fn velocity(&self) -> Vec3;
    fn set_velocity(&mut self, v: Vec3);
    fn mass(&self) -> f64;
    fn weight(&self) -> f64;
    fn as_charged(&self) -> ChargedParticle;
}

impl CollisionParticle for ChargedParticle {
    fn velocity(&self) -> Vec3 {
        self.velocity_m_per_s
    }
    fn set_velocity(&mut self, v: Vec3) {
        self.velocity_m_per_s = v;
    }
    fn mass(&self) -> f64 {
        self.mass_kg
    }
    fn weight(&self) -> f64 {
        self.weight
    }
    fn as_charged(&self) -> ChargedParticle {
        *self
    }
}

pub(crate) fn validate_collision_model(dt: f64, m: &NeutralCollisionModel) -> Result<()> {
    require(
        dt > 0.0
            && dt.is_finite()
            && m.neutral_density_per_m3 >= 0.0
            && m.elastic_cross_section_m2 >= 0.0
            && m.ionization_cross_section_m2 >= 0.0
            && m.ionization_energy_ev >= 0.0,
        "invalid Monte-Carlo collision model",
    )
}

/// Upstream collision kernel shared by the 1-D/2-D/3-D PIC classes.
pub(crate) fn collide<P: CollisionParticle>(
    particles: &mut Vec<P>,
    dt: f64,
    m: &NeutralCollisionModel,
    what: &str,
) -> Result<CollisionStatistics> {
    validate_collision_model(dt, m)?;
    let mut rng = Mt19937_64::new(m.random_seed);
    let mut stats = CollisionStatistics::default();
    let initial = particles.len();
    for i in 0..initial {
        let p = &mut particles[i];
        if p.mass().is_nan() || p.mass() <= 0.0 {
            return Err(SolverError::InvalidInput(what.into()));
        }
        let speed = norm2(p.velocity()).sqrt();
        let sigma = m.elastic_cross_section_m2 + m.ionization_cross_section_m2;
        if speed == 0.0 || sigma == 0.0 || m.neutral_density_per_m3 == 0.0 {
            continue;
        }
        let probability = 1.0 - (-m.neutral_density_per_m3 * sigma * speed * dt).exp();
        if rng.uniform(0.0, 1.0) > probability {
            continue;
        }
        let ionize = rng.uniform(0.0, 1.0) < m.ionization_cross_section_m2 / sigma;
        if ionize {
            let energy = kinetic_energy_j(&p.as_charged())?;
            let loss = m.ionization_energy_ev * ELEMENTARY_CHARGE;
            if energy > loss {
                let gamma_new =
                    1.0 + (energy - loss) / (p.mass() * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
                let beta2 = (1.0 - 1.0 / (gamma_new * gamma_new)).max(0.0);
                let factor = SPEED_OF_LIGHT * beta2.sqrt() / speed;
                p.set_velocity(scale(p.velocity(), factor));
                let mut secondary = p.clone();
                secondary.set_velocity(scale(p.velocity(), -0.05));
                stats.energy_loss_j += loss * p.weight();
                particles.push(secondary);
                stats.ionization_events += 1;
                continue;
            }
        }
        let mu = 2.0 * rng.uniform(0.0, 1.0) - 1.0;
        let phi = 2.0 * PI * rng.uniform(0.0, 1.0);
        let sint = (1.0 - mu * mu).max(0.0).sqrt();
        particles[i].set_velocity(Vec3::new(
            speed * sint * phi.cos(),
            speed * sint * phi.sin(),
            speed * mu,
        ));
        stats.elastic_events += 1;
    }
    Ok(stats)
}

/// Null-collision Monte Carlo: elastic events scatter isotropically at fixed speed; ionization
/// removes the threshold energy and appends a slow secondary macro-particle.
pub fn apply_monte_carlo_collisions(
    particles: &mut Vec<ChargedParticle>,
    dt: f64,
    model: &NeutralCollisionModel,
) -> Result<CollisionStatistics> {
    collide(particles, dt, model, "invalid collision particle mass")
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisAlignedParticleBox {
    pub minimum_m: Vec3,
    pub maximum_m: Vec3,
}

impl Default for AxisAlignedParticleBox {
    fn default() -> Self {
        Self {
            minimum_m: Vec3::default(),
            maximum_m: Vec3::new(1.0, 1.0, 1.0),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ParticleWallMode {
    #[default]
    Absorb,
    SpecularReflect,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SecondaryEmissionModel {
    pub maximum_yield: f64,
    pub energy_at_maximum_ev: f64,
    pub threshold_energy_ev: f64,
}

impl Default for SecondaryEmissionModel {
    fn default() -> Self {
        Self {
            maximum_yield: 1.5,
            energy_at_maximum_ev: 300.0,
            threshold_energy_ev: 20.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleWallInteraction {
    pub alive: bool,
    pub impacted: bool,
    pub impact_normal: Vec3,
    pub secondary_macro_weight: f64,
}

impl Default for ParticleWallInteraction {
    fn default() -> Self {
        Self {
            alive: true,
            impacted: false,
            impact_normal: Vec3::default(),
            secondary_macro_weight: 0.0,
        }
    }
}

/// `δ(E) = δmax·x·e^{1−x}` with `x = (E − E_th)/(E_max − E_th)`, zero below threshold.
pub fn secondary_electron_yield(energy_ev: f64, m: &SecondaryEmissionModel) -> Result<f64> {
    require(
        m.maximum_yield >= 0.0
            && m.energy_at_maximum_ev > m.threshold_energy_ev
            && m.threshold_energy_ev >= 0.0
            && energy_ev.is_finite(),
        "invalid secondary-emission model/energy",
    )?;
    if energy_ev <= m.threshold_energy_ev {
        return Ok(0.0);
    }
    let x = (energy_ev - m.threshold_energy_ev) / (m.energy_at_maximum_ev - m.threshold_energy_ev);
    Ok(m.maximum_yield * x * (1.0 - x).exp())
}

/// Mirror `position`/`velocity` into `[lo, hi]` repeatedly (specular reflection).
pub(crate) fn reflect_into(position: &mut f64, velocity: &mut f64, lo: f64, hi: f64) {
    while *position < lo || *position > hi {
        if *position < lo {
            *position = 2.0 * lo - *position;
            *velocity = -*velocity;
        } else if *position > hi {
            *position = 2.0 * hi - *position;
            *velocity = -*velocity;
        }
    }
}

/// Apply an axis-aligned wall after a push: absorbing walls report an expected secondary
/// macro-weight; reflecting walls mirror position and the normal velocity.
pub fn apply_particle_box_boundary(
    p: &mut ChargedParticle,
    b: &AxisAlignedParticleBox,
    mode: ParticleWallMode,
    secondary: Option<&SecondaryEmissionModel>,
) -> Result<ParticleWallInteraction> {
    require(
        b.maximum_m.x > b.minimum_m.x
            && b.maximum_m.y > b.minimum_m.y
            && b.maximum_m.z > b.minimum_m.z,
        "invalid particle boundary box",
    )?;
    let mut r = ParticleWallInteraction::default();
    let incident_ev = kinetic_energy_j(p)? / ELEMENTARY_CHARGE;
    let outside = |v: f64, lo: f64, hi: f64| v < lo || v > hi;
    if !outside(p.position_m.x, b.minimum_m.x, b.maximum_m.x)
        && !outside(p.position_m.y, b.minimum_m.y, b.maximum_m.y)
        && !outside(p.position_m.z, b.minimum_m.z, b.maximum_m.z)
    {
        return Ok(r);
    }
    r.impacted = true;
    let axis = |position: &mut f64,
                velocity: &mut f64,
                lo: f64,
                hi: f64,
                low: Vec3,
                high: Vec3,
                r: &mut ParticleWallInteraction| {
        if *position >= lo && *position <= hi {
            return;
        }
        r.impact_normal = if *position < lo { low } else { high };
        if mode == ParticleWallMode::Absorb {
            r.alive = false;
            return;
        }
        reflect_into(position, velocity, lo, hi);
    };
    axis(
        &mut p.position_m.x,
        &mut p.velocity_m_per_s.x,
        b.minimum_m.x,
        b.maximum_m.x,
        Vec3::new(-1.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        &mut r,
    );
    if r.alive {
        axis(
            &mut p.position_m.y,
            &mut p.velocity_m_per_s.y,
            b.minimum_m.y,
            b.maximum_m.y,
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            &mut r,
        );
    }
    if r.alive {
        axis(
            &mut p.position_m.z,
            &mut p.velocity_m_per_s.z,
            b.minimum_m.z,
            b.maximum_m.z,
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, 0.0, 1.0),
            &mut r,
        );
    }
    if !r.alive {
        if let Some(s) = secondary {
            r.secondary_macro_weight = p.weight * secondary_electron_yield(incident_ev, s)?;
        }
    }
    Ok(r)
}

pub(crate) fn wrap_periodic(value: f64, length: f64) -> f64 {
    let v = value % length;
    if v < 0.0 {
        v + length
    } else {
        v
    }
}
