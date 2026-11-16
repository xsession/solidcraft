//! Transport-agnostic distributed DEM: slab ownership/ghost planning, lower-ID pair ownership
//! for contacts and bonds with explicit remote reactions, force accumulation and integration.
use super::{
    cross, div, dot, integrate_orientation, inverse_rotate, norm, norm_squared, rotate,
    BondedParticleModel, HertzMindlinContactModel, RigidBody, RigidBodyState, SphereShape, Vec3,
};
use crate::{require, Result, SolverError};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistributedDemParticle {
    pub global_id: u64,
    pub state: RigidBodyState,
    pub mass: f64,
    pub inertia_diagonal: Vec3,
    pub radius: f64,
    pub owner_rank: i32,
}

impl Default for DistributedDemParticle {
    fn default() -> Self {
        Self {
            global_id: 0,
            state: RigidBodyState::default(),
            mass: 1.0,
            inertia_diagonal: Vec3::new(1.0, 1.0, 1.0),
            radius: 0.5,
            owner_rank: 0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DemSlabDecomposition {
    pub x_min: f64,
    pub x_max: f64,
    pub ranks: i32,
}

impl Default for DemSlabDecomposition {
    fn default() -> Self {
        Self {
            x_min: 0.0,
            x_max: 1.0,
            ranks: 1,
        }
    }
}

impl DemSlabDecomposition {
    pub fn validate(&self) -> Result<()> {
        require(
            self.x_max > self.x_min && self.x_min.is_finite() && self.x_max.is_finite(),
            "DEM slab bounds must be finite and ordered",
        )?;
        require(
            self.ranks > 0,
            "DEM slab decomposition requires at least one rank",
        )
    }
    pub fn slab_width(&self) -> Result<f64> {
        self.validate()?;
        Ok((self.x_max - self.x_min) / self.ranks as f64)
    }
    pub fn owner_rank(&self, x: f64) -> Result<i32> {
        self.validate()?;
        require(x.is_finite(), "DEM particle coordinate must be finite")?;
        if x <= self.x_min {
            return Ok(0);
        }
        if x >= self.x_max {
            return Ok(self.ranks - 1);
        }
        let rel = (x - self.x_min) / (self.x_max - self.x_min);
        Ok(((rel * self.ranks as f64).floor() as i32).clamp(0, self.ranks - 1))
    }
    pub fn slab_min(&self, rank: i32) -> Result<f64> {
        self.validate()?;
        require(rank >= 0 && rank < self.ranks, "DEM slab rank out of range")?;
        Ok(self.x_min + rank as f64 * self.slab_width()?)
    }
    pub fn slab_max(&self, rank: i32) -> Result<f64> {
        self.validate()?;
        require(rank >= 0 && rank < self.ranks, "DEM slab rank out of range")?;
        Ok(if rank + 1 == self.ranks {
            self.x_max
        } else {
            self.x_min + (rank + 1) as f64 * self.slab_width()?
        })
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DemExchangePlan {
    pub migrate_indices_by_rank: Vec<Vec<usize>>,
    pub ghost_indices_by_rank: Vec<Vec<usize>>,
}

/// Canonical `particle_a < particle_b`; the owner of `particle_a` owns the state.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DistributedDemContactHistory {
    pub particle_a: u64,
    pub particle_b: u64,
    pub tangential_displacement: Vec3,
    pub age: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DistributedDemBond {
    pub particle_a: u64,
    pub particle_b: u64,
    pub rest_length: f64,
    pub model: BondedParticleModel,
    pub tangential_displacement: Vec3,
    pub damage: f64,
    pub broken: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DistributedDemForceContribution {
    pub global_id: u64,
    pub owner_rank: i32,
    pub force: Vec3,
    pub torque: Vec3,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DistributedDemInteractionResult {
    pub local: Vec<DistributedDemForceContribution>,
    pub remote: Vec<DistributedDemForceContribution>,
    pub contacts: usize,
    pub active_bonds: usize,
    pub broken_bonds: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DistributedDemStepStats {
    pub owned: usize,
    pub ghosts: usize,
    pub contacts: usize,
    pub active_bonds: usize,
    pub broken_bonds: usize,
    pub remote_force_messages: usize,
}

fn finite(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

fn canonical(a: u64, b: u64) -> Result<(u64, u64)> {
    require(a != b, "distributed DEM pair requires distinct global IDs")?;
    Ok(if a < b { (a, b) } else { (b, a) })
}

fn validate_particle(p: &DistributedDemParticle) -> Result<()> {
    let i = p.inertia_diagonal;
    require(
        p.mass > 0.0
            && p.mass.is_finite()
            && p.radius > 0.0
            && p.radius.is_finite()
            && i.x > 0.0
            && i.y > 0.0
            && i.z > 0.0
            && finite(p.state.position)
            && finite(p.state.linear_velocity)
            && finite(p.state.angular_velocity),
        "invalid distributed DEM particle state",
    )
}

fn validate_contact_model(m: &HertzMindlinContactModel) -> Result<()> {
    let nn = |v: f64| v >= 0.0 && v.is_finite();
    require(
        nn(m.normal_stiffness)
            && nn(m.normal_damping)
            && nn(m.tangential_stiffness)
            && nn(m.tangential_damping)
            && nn(m.friction)
            && nn(m.rolling_resistance)
            && nn(m.cohesion_force),
        "invalid distributed DEM contact coefficients",
    )
}

fn validate_bond_model(m: &BondedParticleModel) -> Result<()> {
    let nn = |v: f64| v >= 0.0 && v.is_finite();
    require(
        nn(m.normal_stiffness)
            && nn(m.shear_stiffness)
            && nn(m.normal_damping)
            && nn(m.shear_damping)
            && m.tensile_failure_force > 0.0
            && m.tensile_failure_force.is_finite()
            && m.shear_failure_force > 0.0
            && m.shear_failure_force.is_finite()
            && m.damage_onset_ratio >= 0.0
            && m.damage_onset_ratio < 1.0,
        "invalid distributed DEM bond constitutive parameters",
    )
}

fn add_contribution(
    dest: &mut Vec<DistributedDemForceContribution>,
    id: u64,
    owner: i32,
    force: Vec3,
    torque: Vec3,
) -> Result<()> {
    if let Some(e) = dest.iter_mut().find(|e| e.global_id == id) {
        if e.owner_rank != owner {
            return Err(SolverError::InvalidInput(
                "inconsistent distributed DEM force owner".into(),
            ));
        }
        e.force += force;
        e.torque += torque;
        return Ok(());
    }
    dest.push(DistributedDemForceContribution {
        global_id: id,
        owner_rank: owner,
        force,
        torque,
    });
    Ok(())
}

fn point_velocity(p: &DistributedDemParticle, x: Vec3) -> Vec3 {
    p.state.linear_velocity + cross(p.state.angular_velocity, x - p.state.position)
}

pub fn make_distributed_dem_particle(
    global_id: u64,
    body: &RigidBody,
    sphere: &SphereShape,
    owner_rank: i32,
) -> Result<DistributedDemParticle> {
    require(
        sphere.radius > 0.0 && sphere.radius.is_finite(),
        "distributed DEM sphere radius must be finite and positive",
    )?;
    body.validate()?;
    Ok(DistributedDemParticle {
        global_id,
        state: body.state,
        mass: body.mass,
        inertia_diagonal: body.inertia_diagonal,
        radius: sphere.radius,
        owner_rank,
    })
}

/// `rest_length = 0` uses the current distance.
pub fn make_distributed_dem_bond(
    a: &DistributedDemParticle,
    b: &DistributedDemParticle,
    model: BondedParticleModel,
    rest_length: f64,
) -> Result<DistributedDemBond> {
    validate_particle(a)?;
    validate_particle(b)?;
    validate_bond_model(&model)?;
    let (ia, ib) = canonical(a.global_id, b.global_id)?;
    let d = norm(b.state.position - a.state.position);
    let rest = if rest_length == 0.0 { d } else { rest_length };
    require(
        rest > 0.0 && rest.is_finite(),
        "distributed DEM bond rest length must be finite and positive",
    )?;
    Ok(DistributedDemBond {
        particle_a: ia,
        particle_b: ib,
        rest_length: rest,
        model,
        ..Default::default()
    })
}

/// Ownership by particle center; ghosts expand the sphere by `ghost_width` into neighbor slabs.
pub fn plan_dem_slab_exchange(
    particles: &[DistributedDemParticle],
    decomposition: &DemSlabDecomposition,
    local_rank: i32,
    ghost_width: f64,
) -> Result<DemExchangePlan> {
    decomposition.validate()?;
    require(
        local_rank >= 0 && local_rank < decomposition.ranks,
        "local DEM rank out of range",
    )?;
    require(
        ghost_width >= 0.0
            && ghost_width.is_finite()
            && ghost_width <= decomposition.slab_width()?,
        "DEM ghost width must be finite and no larger than one slab",
    )?;
    require(decomposition.ranks <= 1 << 20, "too many DEM ranks")?;
    let r = decomposition.ranks as usize;
    let mut plan = DemExchangePlan {
        migrate_indices_by_rank: vec![vec![]; r],
        ghost_indices_by_rank: vec![vec![]; r],
    };
    let (lower, upper) = (
        decomposition.slab_min(local_rank)?,
        decomposition.slab_max(local_rank)?,
    );
    for (i, p) in particles.iter().enumerate() {
        validate_particle(p)?;
        let owner = decomposition.owner_rank(p.state.position.x)?;
        if owner != local_rank {
            plan.migrate_indices_by_rank[owner as usize].push(i);
            continue;
        }
        let reach = p.radius + ghost_width;
        if local_rank > 0 && p.state.position.x - reach < lower {
            plan.ghost_indices_by_rank[(local_rank - 1) as usize].push(i);
        }
        if local_rank + 1 < decomposition.ranks && p.state.position.x + reach > upper {
            plan.ghost_indices_by_rank[(local_rank + 1) as usize].push(i);
        }
    }
    Ok(plan)
}

#[allow(clippy::too_many_arguments)]
fn evaluate_pair(
    a: &DistributedDemParticle,
    b: &DistributedDemParticle,
    b_remote: bool,
    rank: i32,
    dt: f64,
    m: &HertzMindlinContactModel,
    old: &[DistributedDemContactHistory],
    next: &mut Vec<DistributedDemContactHistory>,
    result: &mut DistributedDemInteractionResult,
) -> Result<()> {
    let delta = a.state.position - b.state.position;
    let d2 = norm_squared(delta);
    let rs = a.radius + b.radius;
    if d2.is_nan() || d2 >= rs * rs {
        return Ok(());
    }
    if !(d2 > 1e-30 && d2.is_finite()) {
        return Err(SolverError::NumericalFailure(
            "distributed DEM overlapping particle centers are coincident/non-finite".into(),
        ));
    }
    let distance = d2.sqrt();
    let normal = div(delta, distance);
    let pen = rs - distance;
    let root = pen.max(0.0).sqrt();
    let cp = a.state.position - normal * (a.radius - 0.5 * pen);
    let (ra, rb) = (cp - a.state.position, cp - b.state.position);
    let rel = point_velocity(a, cp) - point_velocity(b, cp);
    let vn = dot(rel, normal);
    let vt = rel - normal * vn;
    let (ia, ib) = canonical(a.global_id, b.global_id)?;
    let mut state = old
        .iter()
        .find(|s| s.particle_a == ia && s.particle_b == ib)
        .copied()
        .unwrap_or(DistributedDemContactHistory {
            particle_a: ia,
            particle_b: ib,
            ..Default::default()
        });
    state.tangential_displacement += vt * dt;
    let td = state.tangential_displacement;
    state.tangential_displacement -= normal * dot(td, normal);
    state.age += 1;
    let fn_ = (m.normal_stiffness * pen * root - m.normal_damping * root * vn - m.cohesion_force)
        .max(-m.cohesion_force);
    let (kt, ct) = (m.tangential_stiffness * root, m.tangential_damping * root);
    let mut ft = -(state.tangential_displacement * kt + vt * ct);
    let ftm = norm(ft);
    let limit = m.friction * fn_.abs();
    if ftm > limit && ftm > 1e-30 {
        ft = ft * (limit / ftm);
        if kt > 1e-30 {
            state.tangential_displacement = div(-(ft + vt * ct), kt);
        }
    }
    next.push(state);
    let force = normal * fn_ + ft;
    let mut ta = cross(ra, force);
    let mut tb = cross(rb, -force);
    let rw = a.state.angular_velocity - b.state.angular_velocity;
    let wm = norm(rw);
    if wm > 1e-30 && m.rolling_resistance > 0.0 {
        let reff = a.radius * b.radius / (a.radius + b.radius);
        let rolling = rw * (-m.rolling_resistance * fn_.abs() * reff / wm);
        ta += rolling;
        tb -= rolling;
    }
    add_contribution(&mut result.local, a.global_id, rank, force, ta)?;
    if b_remote {
        add_contribution(&mut result.remote, b.global_id, b.owner_rank, -force, tb)?;
    } else {
        add_contribution(&mut result.local, b.global_id, rank, -force, tb)?;
    }
    result.contacts += 1;
    Ok(())
}

/// Local-local and local-ghost Hertz-Mindlin contacts; a cross-rank pair is evaluated only by the
/// owner of the lower global id, and remote reactions are returned for reverse exchange.
pub fn evaluate_distributed_dem_contacts(
    owned: &[DistributedDemParticle],
    ghosts: &[DistributedDemParticle],
    history: &mut Vec<DistributedDemContactHistory>,
    local_rank: i32,
    dt: f64,
    model: &HertzMindlinContactModel,
) -> Result<DistributedDemInteractionResult> {
    require(
        dt > 0.0 && dt.is_finite(),
        "distributed DEM contact dt must be finite and positive",
    )?;
    validate_contact_model(model)?;
    for p in owned {
        validate_particle(p)?;
        require(
            p.owner_rank == local_rank,
            "owned DEM particle has wrong owner rank",
        )?;
    }
    for p in ghosts {
        validate_particle(p)?;
    }
    let mut result = DistributedDemInteractionResult::default();
    let mut next = Vec::with_capacity(history.len() + owned.len());
    for i in 0..owned.len() {
        for j in i + 1..owned.len() {
            let (mut a, mut b) = (&owned[i], &owned[j]);
            if a.global_id > b.global_id {
                std::mem::swap(&mut a, &mut b);
            }
            evaluate_pair(
                a,
                b,
                false,
                local_rank,
                dt,
                model,
                history,
                &mut next,
                &mut result,
            )?;
        }
    }
    for l in owned {
        for g in ghosts {
            if l.global_id >= g.global_id {
                continue;
            }
            evaluate_pair(
                l,
                g,
                true,
                local_rank,
                dt,
                model,
                history,
                &mut next,
                &mut result,
            )?;
        }
    }
    *history = next;
    Ok(result)
}

/// Bonded-particle forces with the same lower-id ownership rule (partners may be ghosts).
pub fn evaluate_distributed_dem_bonds(
    owned: &[DistributedDemParticle],
    ghosts: &[DistributedDemParticle],
    bonds: &mut [DistributedDemBond],
    local_rank: i32,
    dt: f64,
) -> Result<DistributedDemInteractionResult> {
    require(
        dt > 0.0 && dt.is_finite(),
        "distributed DEM bond dt must be finite and positive",
    )?;
    let mut result = DistributedDemInteractionResult::default();
    for bond in bonds.iter_mut() {
        validate_bond_model(&bond.model)?;
        require(
            bond.particle_a < bond.particle_b,
            "distributed DEM bond IDs must be canonical",
        )?;
        if bond.broken {
            result.broken_bonds += 1;
            continue;
        }
        let Some(a) = owned.iter().find(|p| p.global_id == bond.particle_a) else {
            continue;
        };
        if a.owner_rank != local_rank {
            return Err(SolverError::InvalidInput(
                "distributed DEM bond owner mismatch".into(),
            ));
        }
        let (b, remote) = match owned.iter().find(|p| p.global_id == bond.particle_b) {
            Some(b) => (b, false),
            None => match ghosts.iter().find(|p| p.global_id == bond.particle_b) {
                Some(b) => (b, true),
                None => continue,
            },
        };
        let delta = b.state.position - a.state.position;
        let distance = norm(delta);
        if !(distance > 1e-15 && distance.is_finite()) {
            return Err(SolverError::NumericalFailure(
                "distributed bonded particle centers are coincident/non-finite".into(),
            ));
        }
        let n = div(delta, distance);
        let rv = b.state.linear_velocity - a.state.linear_velocity;
        let vn = dot(rv, n);
        let vt = rv - n * vn;
        bond.tangential_displacement += vt * dt;
        let td = bond.tangential_displacement;
        bond.tangential_displacement -= n * dot(td, n);
        let m = bond.model;
        let ne = m.normal_stiffness * (distance - bond.rest_length);
        let se = bond.tangential_displacement * m.shear_stiffness;
        let ratio = (ne.max(0.0) / m.tensile_failure_force).max(norm(se) / m.shear_failure_force);
        if ratio > m.damage_onset_ratio {
            bond.damage = bond.damage.max(
                ((ratio - m.damage_onset_ratio) / (1.0 - m.damage_onset_ratio)).clamp(0.0, 1.0),
            );
        }
        if ratio >= 1.0 || bond.damage >= 1.0 {
            bond.damage = 1.0;
            bond.broken = true;
            result.broken_bonds += 1;
            continue;
        }
        result.active_bonds += 1;
        let intact = 1.0 - bond.damage;
        let f = n * (intact * ne + m.normal_damping * vn) + (se * intact + vt * m.shear_damping);
        add_contribution(&mut result.local, a.global_id, local_rank, f, Vec3::ZERO)?;
        if remote {
            add_contribution(
                &mut result.remote,
                b.global_id,
                b.owner_rank,
                -f,
                Vec3::ZERO,
            )?;
        } else {
            add_contribution(&mut result.local, b.global_id, local_rank, -f, Vec3::ZERO)?;
        }
    }
    Ok(result)
}

pub fn accumulate_distributed_dem_forces(
    destination: &mut Vec<DistributedDemForceContribution>,
    source: &[DistributedDemForceContribution],
) -> Result<()> {
    for c in source {
        add_contribution(destination, c.global_id, c.owner_rank, c.force, c.torque)?;
    }
    Ok(())
}

/// Semi-implicit Euler for owned spheres with body-frame diagonal inertia; on error the
/// particles are unchanged.
pub fn integrate_distributed_dem_particles(
    owned: &mut [DistributedDemParticle],
    contributions: &[DistributedDemForceContribution],
    gravity: Vec3,
    dt: f64,
) -> Result<()> {
    require(
        dt > 0.0 && dt.is_finite() && finite(gravity),
        "invalid distributed DEM integration parameters",
    )?;
    let mut out = owned.to_vec();
    for p in out.iter_mut() {
        validate_particle(p)?;
        let mut force = gravity * p.mass;
        let mut torque = Vec3::ZERO;
        for c in contributions.iter().filter(|c| c.global_id == p.global_id) {
            force += c.force;
            torque += c.torque;
        }
        p.state.linear_velocity += force * (dt / p.mass);
        let tl = inverse_rotate(p.state.orientation, torque)?;
        let mut wl = inverse_rotate(p.state.orientation, p.state.angular_velocity)?;
        wl.x += dt * tl.x / p.inertia_diagonal.x;
        wl.y += dt * tl.y / p.inertia_diagonal.y;
        wl.z += dt * tl.z / p.inertia_diagonal.z;
        p.state.angular_velocity = rotate(p.state.orientation, wl)?;
        p.state.position += p.state.linear_velocity * dt;
        integrate_orientation(&mut p.state.orientation, p.state.angular_velocity, dt)?;
    }
    owned.copy_from_slice(&out);
    Ok(())
}
