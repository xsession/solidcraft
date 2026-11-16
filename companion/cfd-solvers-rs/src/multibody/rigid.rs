//! Rigid bodies, joints and motors (projected Gauss-Seidel rows), broad and narrow phase contact
//! (spheres, planes, triangle meshes, GJK/EPA convex hulls), penalty/Hertz-Mindlin/impulse
//! contact laws, bonded particles and the rigid-body and explicit DEM systems.
#![allow(clippy::needless_range_loop)]
use super::{
    any_perpendicular, component_mul, cross, div, dot, integrate_orientation, inverse_rotate, norm,
    norm_squared, normalized_or_zero, normalized_quaternion, rotate, Quaternion, Vec3,
};
use crate::{require, Result, SolverError};

fn nrm(v: Vec3) -> Vec3 {
    normalized_or_zero(v, 1e-15)
}

fn bad(msg: &str) -> SolverError {
    SolverError::InvalidInput(msg.into())
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RigidBodyState {
    pub position: Vec3,
    pub orientation: Quaternion,
    pub linear_velocity: Vec3,
    pub angular_velocity: Vec3,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidBody {
    pub state: RigidBodyState,
    pub mass: f64,
    pub inertia_diagonal: Vec3,
    pub fixed: bool,
    pub force: Vec3,
    pub torque: Vec3,
}

impl Default for RigidBody {
    fn default() -> Self {
        Self {
            state: RigidBodyState::default(),
            mass: 1.0,
            inertia_diagonal: Vec3::new(1.0, 1.0, 1.0),
            fixed: false,
            force: Vec3::ZERO,
            torque: Vec3::ZERO,
        }
    }
}

impl RigidBody {
    pub fn inverse_mass(&self) -> f64 {
        if self.fixed {
            0.0
        } else {
            1.0 / self.mass
        }
    }
    pub fn world_inverse_inertia_mul(&self, v: Vec3) -> Result<Vec3> {
        if self.fixed {
            return Ok(Vec3::ZERO);
        }
        let l = inverse_rotate(self.state.orientation, v)?;
        let i = self.inertia_diagonal;
        rotate(
            self.state.orientation,
            Vec3::new(l.x / i.x, l.y / i.y, l.z / i.z),
        )
    }
    pub fn point_velocity(&self, p: Vec3) -> Vec3 {
        self.state.linear_velocity + cross(self.state.angular_velocity, p - self.state.position)
    }
    pub fn add_force(&mut self, f: Vec3) {
        self.force += f;
    }
    pub fn add_torque(&mut self, t: Vec3) {
        self.torque += t;
    }
    pub fn add_force_at_point(&mut self, f: Vec3, p: Vec3) {
        self.force += f;
        self.torque += cross(p - self.state.position, f);
    }
    pub fn apply_impulse(&mut self, impulse: Vec3, p: Vec3) -> Result<()> {
        if self.fixed {
            return Ok(());
        }
        self.state.linear_velocity += impulse * self.inverse_mass();
        let dw = self.world_inverse_inertia_mul(cross(p - self.state.position, impulse))?;
        self.state.angular_velocity += dw;
        Ok(())
    }
    pub fn clear_accumulators(&mut self) {
        self.force = Vec3::ZERO;
        self.torque = Vec3::ZERO;
    }
    pub fn validate(&self) -> Result<()> {
        require(
            self.mass > 0.0 && self.mass.is_finite(),
            "rigid body mass must be finite and positive",
        )?;
        let i = self.inertia_diagonal;
        require(
            i.x > 0.0 && i.y > 0.0 && i.z > 0.0,
            "rigid body inertia must be positive",
        )?;
        require(
            norm_squared(self.state.position).is_finite()
                && norm_squared(self.state.linear_velocity).is_finite()
                && norm_squared(self.state.angular_velocity).is_finite(),
            "rigid body state must be finite",
        )?;
        normalized_quaternion(self.state.orientation).map(|_| ())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RigidBodyIntegrator {
    #[default]
    SemiImplicitEuler,
    VelocityVerlet,
}

/// Includes the gyroscopic term; fixed bodies only clear their accumulators.
pub fn integrate_rigid_body(
    body: &mut RigidBody,
    dt: f64,
    integrator: RigidBodyIntegrator,
) -> Result<()> {
    if body.fixed {
        body.clear_accumulators();
        return Ok(());
    }
    body.validate()?;
    require(
        dt > 0.0 && dt.is_finite(),
        "rigid body dt must be finite and positive",
    )?;
    let s = body.state;
    let a = body.force * body.inverse_mass();
    let gyro = cross(
        s.angular_velocity,
        rotate(
            s.orientation,
            component_mul(
                body.inertia_diagonal,
                inverse_rotate(s.orientation, s.angular_velocity)?,
            ),
        )?,
    );
    let alpha = body.world_inverse_inertia_mul(body.torque - gyro)?;
    let st = &mut body.state;
    match integrator {
        RigidBodyIntegrator::VelocityVerlet => {
            st.position += st.linear_velocity * dt + a * (0.5 * dt * dt);
            st.linear_velocity += a * dt;
            st.angular_velocity += alpha * dt;
            integrate_orientation(
                &mut st.orientation,
                st.angular_velocity - alpha * (0.5 * dt),
                dt,
            )?;
        }
        RigidBodyIntegrator::SemiImplicitEuler => {
            st.linear_velocity += a * dt;
            st.angular_velocity += alpha * dt;
            st.position += st.linear_velocity * dt;
            integrate_orientation(&mut st.orientation, st.angular_velocity, dt)?;
        }
    }
    body.clear_accumulators();
    Ok(())
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SphereShape {
    pub body: usize,
    pub radius: f64,
}

/// `dot(normal, x) = offset`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneShape {
    pub normal: Vec3,
    pub offset: f64,
}

impl Default for PlaneShape {
    fn default() -> Self {
        Self {
            normal: Vec3::new(0.0, 1.0, 0.0),
            offset: 0.0,
        }
    }
}

/// Jacobian row; `None` bodies are the static world.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConstraintRow {
    pub body_a: Option<usize>,
    pub body_b: Option<usize>,
    pub linear_a: Vec3,
    pub angular_a: Vec3,
    pub linear_b: Vec3,
    pub angular_b: Vec3,
    pub bias: f64,
    pub lower: f64,
    pub upper: f64,
    pub impulse: f64,
}

impl Default for ConstraintRow {
    fn default() -> Self {
        Self {
            body_a: None,
            body_b: None,
            linear_a: Vec3::ZERO,
            angular_a: Vec3::ZERO,
            linear_b: Vec3::ZERO,
            angular_b: Vec3::ZERO,
            bias: 0.0,
            lower: f64::NEG_INFINITY,
            upper: f64::INFINITY,
            impulse: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DistanceConstraint {
    pub body_a: usize,
    pub body_b: usize,
    pub local_anchor_a: Vec3,
    pub local_anchor_b: Vec3,
    pub distance: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SphericalJoint {
    pub body_a: usize,
    pub body_b: usize,
    pub local_anchor_a: Vec3,
    pub local_anchor_b: Vec3,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RevoluteJoint {
    pub body_a: usize,
    pub body_b: usize,
    pub local_anchor_a: Vec3,
    pub local_anchor_b: Vec3,
    pub local_axis_a: Vec3,
    pub local_axis_b: Vec3,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrismaticJoint {
    pub body_a: usize,
    pub body_b: usize,
    pub local_anchor_a: Vec3,
    pub local_anchor_b: Vec3,
    pub local_axis_a: Vec3,
    pub reference_relative_orientation: Quaternion,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FixedJoint {
    pub body_a: usize,
    pub body_b: usize,
    pub local_anchor_a: Vec3,
    pub local_anchor_b: Vec3,
    pub reference_relative_orientation: Quaternion,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GearConstraint {
    pub body_a: usize,
    pub body_b: usize,
    pub local_axis_a: Vec3,
    pub local_axis_b: Vec3,
    pub ratio: f64,
}

/// Velocity motor as a bounded row; the force limit becomes a per-step impulse limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearMotor {
    pub body_a: usize,
    pub body_b: usize,
    pub local_axis_a: Vec3,
    pub target_speed: f64,
    pub max_force: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularMotor {
    pub body_a: usize,
    pub body_b: usize,
    pub local_axis_a: Vec3,
    pub target_speed: f64,
    pub max_torque: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Joint {
    Distance(DistanceConstraint),
    Spherical(SphericalJoint),
    Revolute(RevoluteJoint),
    Prismatic(PrismaticJoint),
    Fixed(FixedJoint),
    Gear(GearConstraint),
    LinearMotor(LinearMotor),
    AngularMotor(AngularMotor),
}

fn body(bodies: &[RigidBody], i: usize) -> Result<&RigidBody> {
    bodies
        .get(i)
        .ok_or_else(|| bad("rigid body index out of range"))
}

fn body_point(b: &RigidBody, local: Vec3) -> Result<Vec3> {
    Ok(b.state.position + rotate(b.state.orientation, local)?)
}

fn local_point(b: &RigidBody, world: Vec3) -> Result<Vec3> {
    inverse_rotate(b.state.orientation, world - b.state.position)
}

fn world_axis(b: &RigidBody, local: Vec3) -> Result<Vec3> {
    Ok(nrm(rotate(b.state.orientation, local)?))
}

fn relative_orientation(a: &RigidBody, b: &RigidBody) -> Result<Quaternion> {
    normalized_quaternion(
        normalized_quaternion(a.state.orientation)?.conjugate()
            * normalized_quaternion(b.state.orientation)?,
    )
}

fn orientation_error(current: Quaternion, reference: Quaternion) -> Result<Vec3> {
    let mut e = normalized_quaternion(current * normalized_quaternion(reference)?.conjugate())?;
    if e.w < 0.0 {
        e = Quaternion {
            w: -e.w,
            x: -e.x,
            y: -e.y,
            z: -e.z,
        };
    }
    Ok(Vec3::new(2.0 * e.x, 2.0 * e.y, 2.0 * e.z))
}

pub fn make_distance_constraint(
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    world_anchor_a: Vec3,
    world_anchor_b: Vec3,
    distance: Option<f64>,
) -> Result<DistanceConstraint> {
    let d = match distance {
        Some(d) if d >= 0.0 => d,
        _ => norm(world_anchor_b - world_anchor_a),
    };
    Ok(DistanceConstraint {
        body_a: a,
        body_b: b,
        local_anchor_a: local_point(body(bodies, a)?, world_anchor_a)?,
        local_anchor_b: local_point(body(bodies, b)?, world_anchor_b)?,
        distance: d,
    })
}

pub fn make_spherical_joint(
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    world_anchor: Vec3,
) -> Result<SphericalJoint> {
    Ok(SphericalJoint {
        body_a: a,
        body_b: b,
        local_anchor_a: local_point(body(bodies, a)?, world_anchor)?,
        local_anchor_b: local_point(body(bodies, b)?, world_anchor)?,
    })
}

pub fn make_revolute_joint(
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    world_anchor: Vec3,
    world_axis: Vec3,
) -> Result<RevoluteJoint> {
    let axis = nrm(world_axis);
    let (ba, bb) = (body(bodies, a)?, body(bodies, b)?);
    Ok(RevoluteJoint {
        body_a: a,
        body_b: b,
        local_anchor_a: local_point(ba, world_anchor)?,
        local_anchor_b: local_point(bb, world_anchor)?,
        local_axis_a: inverse_rotate(ba.state.orientation, axis)?,
        local_axis_b: inverse_rotate(bb.state.orientation, axis)?,
    })
}

pub fn make_prismatic_joint(
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    world_anchor: Vec3,
    world_axis: Vec3,
) -> Result<PrismaticJoint> {
    let axis = nrm(world_axis);
    let (ba, bb) = (body(bodies, a)?, body(bodies, b)?);
    Ok(PrismaticJoint {
        body_a: a,
        body_b: b,
        local_anchor_a: local_point(ba, world_anchor)?,
        local_anchor_b: local_point(bb, world_anchor)?,
        local_axis_a: inverse_rotate(ba.state.orientation, axis)?,
        reference_relative_orientation: relative_orientation(ba, bb)?,
    })
}

pub fn make_fixed_joint(
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    world_anchor: Vec3,
) -> Result<FixedJoint> {
    let (ba, bb) = (body(bodies, a)?, body(bodies, b)?);
    Ok(FixedJoint {
        body_a: a,
        body_b: b,
        local_anchor_a: local_point(ba, world_anchor)?,
        local_anchor_b: local_point(bb, world_anchor)?,
        reference_relative_orientation: relative_orientation(ba, bb)?,
    })
}

pub fn make_gear_constraint(
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    world_axis_a: Vec3,
    world_axis_b: Vec3,
    ratio: f64,
) -> Result<GearConstraint> {
    require(
        ratio.is_finite() && ratio != 0.0,
        "gear ratio must be finite and non-zero",
    )?;
    Ok(GearConstraint {
        body_a: a,
        body_b: b,
        local_axis_a: inverse_rotate(body(bodies, a)?.state.orientation, nrm(world_axis_a))?,
        local_axis_b: inverse_rotate(body(bodies, b)?.state.orientation, nrm(world_axis_b))?,
        ratio,
    })
}

/// `max_force` may be infinite (upstream default).
pub fn make_linear_motor(
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    world_axis: Vec3,
    target_speed: f64,
    max_force: f64,
) -> Result<LinearMotor> {
    require(
        target_speed.is_finite() && max_force >= 0.0,
        "invalid linear motor parameters",
    )?;
    body(bodies, b)?;
    Ok(LinearMotor {
        body_a: a,
        body_b: b,
        local_axis_a: inverse_rotate(body(bodies, a)?.state.orientation, nrm(world_axis))?,
        target_speed,
        max_force,
    })
}

pub fn make_angular_motor(
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    world_axis: Vec3,
    target_speed: f64,
    max_torque: f64,
) -> Result<AngularMotor> {
    require(
        target_speed.is_finite() && max_torque >= 0.0,
        "invalid angular motor parameters",
    )?;
    body(bodies, b)?;
    Ok(AngularMotor {
        body_a: a,
        body_b: b,
        local_axis_a: inverse_rotate(body(bodies, a)?.state.orientation, nrm(world_axis))?,
        target_speed,
        max_torque,
    })
}

#[allow(clippy::too_many_arguments)]
fn add_point_rows(
    rows: &mut Vec<ConstraintRow>,
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    pa: Vec3,
    pb: Vec3,
    dt: f64,
    beta: f64,
) -> Result<()> {
    let ra = pa - body(bodies, a)?.state.position;
    let rb = pb - body(bodies, b)?.state.position;
    let err = pb - pa;
    for axis in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ] {
        let n = nrm(axis);
        rows.push(ConstraintRow {
            body_a: Some(a),
            body_b: Some(b),
            linear_a: -n,
            angular_a: -cross(ra, n),
            linear_b: n,
            angular_b: cross(rb, n),
            bias: (beta / dt) * dot(err, n),
            ..Default::default()
        });
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn add_angular_lock_rows(
    rows: &mut Vec<ConstraintRow>,
    bodies: &[RigidBody],
    a: usize,
    b: usize,
    reference: Quaternion,
    dt: f64,
    beta: f64,
) -> Result<()> {
    let (ba, bb) = (body(bodies, a)?, body(bodies, b)?);
    let local = orientation_error(relative_orientation(ba, bb)?, reference)?;
    let world = rotate(ba.state.orientation, local)?;
    for axis in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ] {
        let n = nrm(axis);
        rows.push(ConstraintRow {
            body_a: Some(a),
            body_b: Some(b),
            angular_a: -n,
            angular_b: n,
            bias: (beta / dt) * dot(world, n),
            ..Default::default()
        });
    }
    Ok(())
}

/// Upstream default Baumgarte factor 0.2.
pub fn build_constraint_rows(
    bodies: &[RigidBody],
    joint: &Joint,
    dt: f64,
    baumgarte: f64,
) -> Result<Vec<ConstraintRow>> {
    require(
        dt > 0.0 && baumgarte >= 0.0,
        "invalid constraint timestep/Baumgarte factor",
    )?;
    let beta = baumgarte;
    let mut rows = vec![];
    let ids = |a: usize, b: usize| -> Result<(&RigidBody, &RigidBody)> {
        Ok((body(bodies, a)?, body(bodies, b)?))
    };
    match *joint {
        Joint::Distance(j) => {
            let (a, b) = ids(j.body_a, j.body_b)?;
            let (pa, pb) = (
                body_point(a, j.local_anchor_a)?,
                body_point(b, j.local_anchor_b)?,
            );
            let d = pb - pa;
            let l = norm(d);
            let n = if l > 1e-12 {
                div(d, l)
            } else {
                Vec3::new(1.0, 0.0, 0.0)
            };
            let (ra, rb) = (pa - a.state.position, pb - b.state.position);
            rows.push(ConstraintRow {
                body_a: Some(j.body_a),
                body_b: Some(j.body_b),
                linear_a: -n,
                angular_a: -cross(ra, n),
                linear_b: n,
                angular_b: cross(rb, n),
                bias: (beta / dt) * (l - j.distance),
                ..Default::default()
            });
        }
        Joint::Spherical(j) => {
            let (a, b) = ids(j.body_a, j.body_b)?;
            add_point_rows(
                &mut rows,
                bodies,
                j.body_a,
                j.body_b,
                body_point(a, j.local_anchor_a)?,
                body_point(b, j.local_anchor_b)?,
                dt,
                beta,
            )?;
        }
        Joint::Revolute(j) => {
            let (a, b) = ids(j.body_a, j.body_b)?;
            add_point_rows(
                &mut rows,
                bodies,
                j.body_a,
                j.body_b,
                body_point(a, j.local_anchor_a)?,
                body_point(b, j.local_anchor_b)?,
                dt,
                beta,
            )?;
            let (aa, ab) = (
                world_axis(a, j.local_axis_a)?,
                world_axis(b, j.local_axis_b)?,
            );
            let t1 = any_perpendicular(aa);
            let t2 = nrm(cross(aa, t1));
            let err = cross(aa, ab);
            for t in [t1, t2] {
                rows.push(ConstraintRow {
                    body_a: Some(j.body_a),
                    body_b: Some(j.body_b),
                    angular_a: -t,
                    angular_b: t,
                    bias: (beta / dt) * dot(err, t),
                    ..Default::default()
                });
            }
        }
        Joint::Prismatic(j) => {
            let (a, b) = ids(j.body_a, j.body_b)?;
            let (pa, pb) = (
                body_point(a, j.local_anchor_a)?,
                body_point(b, j.local_anchor_b)?,
            );
            let axis = world_axis(a, j.local_axis_a)?;
            let t1 = any_perpendicular(axis);
            let t2 = nrm(cross(axis, t1));
            let (ra, rb, err) = (pa - a.state.position, pb - b.state.position, pb - pa);
            for t in [t1, t2] {
                rows.push(ConstraintRow {
                    body_a: Some(j.body_a),
                    body_b: Some(j.body_b),
                    linear_a: -t,
                    angular_a: -cross(ra, t),
                    linear_b: t,
                    angular_b: cross(rb, t),
                    bias: (beta / dt) * dot(err, t),
                    ..Default::default()
                });
            }
            add_angular_lock_rows(
                &mut rows,
                bodies,
                j.body_a,
                j.body_b,
                j.reference_relative_orientation,
                dt,
                beta,
            )?;
        }
        Joint::Fixed(j) => {
            let (a, b) = ids(j.body_a, j.body_b)?;
            add_point_rows(
                &mut rows,
                bodies,
                j.body_a,
                j.body_b,
                body_point(a, j.local_anchor_a)?,
                body_point(b, j.local_anchor_b)?,
                dt,
                beta,
            )?;
            add_angular_lock_rows(
                &mut rows,
                bodies,
                j.body_a,
                j.body_b,
                j.reference_relative_orientation,
                dt,
                beta,
            )?;
        }
        Joint::Gear(j) => {
            let (a, b) = ids(j.body_a, j.body_b)?;
            rows.push(ConstraintRow {
                body_a: Some(j.body_a),
                body_b: Some(j.body_b),
                angular_a: world_axis(a, j.local_axis_a)?,
                angular_b: world_axis(b, j.local_axis_b)? * j.ratio,
                ..Default::default()
            });
        }
        Joint::LinearMotor(j) => {
            let (a, _) = ids(j.body_a, j.body_b)?;
            let axis = world_axis(a, j.local_axis_a)?;
            let limit = if j.max_force.is_finite() {
                j.max_force * dt
            } else {
                f64::INFINITY
            };
            rows.push(ConstraintRow {
                body_a: Some(j.body_a),
                body_b: Some(j.body_b),
                linear_a: -axis,
                linear_b: axis,
                bias: -j.target_speed,
                lower: -limit,
                upper: limit,
                ..Default::default()
            });
        }
        Joint::AngularMotor(j) => {
            let (a, _) = ids(j.body_a, j.body_b)?;
            let axis = world_axis(a, j.local_axis_a)?;
            let limit = if j.max_torque.is_finite() {
                j.max_torque * dt
            } else {
                f64::INFINITY
            };
            rows.push(ConstraintRow {
                body_a: Some(j.body_a),
                body_b: Some(j.body_b),
                angular_a: -axis,
                angular_b: axis,
                bias: -j.target_speed,
                lower: -limit,
                upper: limit,
                ..Default::default()
            });
        }
    }
    Ok(rows)
}

fn row_velocity(bodies: &[RigidBody], r: &ConstraintRow) -> Result<f64> {
    let mut v = 0.0;
    if let Some(a) = r.body_a {
        let b = body(bodies, a)?;
        v += dot(r.linear_a, b.state.linear_velocity) + dot(r.angular_a, b.state.angular_velocity);
    }
    if let Some(i) = r.body_b {
        let b = body(bodies, i)?;
        v += dot(r.linear_b, b.state.linear_velocity) + dot(r.angular_b, b.state.angular_velocity);
    }
    Ok(v)
}

fn row_diag(bodies: &[RigidBody], r: &ConstraintRow) -> Result<f64> {
    let mut d = 0.0;
    if let Some(a) = r.body_a {
        let b = body(bodies, a)?;
        d += b.inverse_mass() * dot(r.linear_a, r.linear_a)
            + dot(r.angular_a, b.world_inverse_inertia_mul(r.angular_a)?);
    }
    if let Some(i) = r.body_b {
        let b = body(bodies, i)?;
        d += b.inverse_mass() * dot(r.linear_b, r.linear_b)
            + dot(r.angular_b, b.world_inverse_inertia_mul(r.angular_b)?);
    }
    Ok(d)
}

fn apply_row(bodies: &mut [RigidBody], r: &ConstraintRow, dl: f64) -> Result<()> {
    for (id, lin, ang) in [
        (r.body_a, r.linear_a, r.angular_a),
        (r.body_b, r.linear_b, r.angular_b),
    ] {
        if let Some(i) = id {
            let b = bodies
                .get_mut(i)
                .ok_or_else(|| bad("rigid body index out of range"))?;
            if !b.fixed {
                b.state.linear_velocity += lin * (dl * b.inverse_mass());
                let dw = b.world_inverse_inertia_mul(ang * dl)?;
                b.state.angular_velocity += dw;
            }
        }
    }
    Ok(())
}

/// Projected Gauss-Seidel over the rows (upstream default 12 iterations).
pub fn solve_constraint_rows(
    bodies: &mut [RigidBody],
    rows: &mut [ConstraintRow],
    iterations: usize,
) -> Result<()> {
    for _ in 0..iterations {
        for r in rows.iter_mut() {
            let d = row_diag(bodies, r)?;
            if d <= 1e-18 {
                continue;
            }
            let old = r.impulse;
            let dl = -(row_velocity(bodies, r)? + r.bias) / d;
            r.impulse = r.lower.max(r.upper.min(old + dl));
            apply_row(bodies, r, r.impulse - old)?;
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Aabb {
    pub minimum: Vec3,
    pub maximum: Vec3,
}

/// `normal` is the outward force direction on body A; `None` bodies are static world geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactPoint {
    pub body_a: Option<usize>,
    pub body_b: Option<usize>,
    pub point: Vec3,
    pub normal: Vec3,
    pub penetration: f64,
    pub feature_id: u64,
}

impl Default for ContactPoint {
    fn default() -> Self {
        Self {
            body_a: None,
            body_b: None,
            point: Vec3::ZERO,
            normal: Vec3::new(1.0, 0.0, 0.0),
            penetration: 0.0,
            feature_id: 0,
        }
    }
}

/// Vertices in body-local coordinates.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConvexHullShape {
    pub body: usize,
    pub vertices: Vec<Vec3>,
}

/// `body: None` is a static world-space mesh.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TriangleMeshShape {
    pub body: Option<usize>,
    pub vertices: Vec<Vec3>,
    pub triangles: Vec<[usize; 3]>,
    pub feature_namespace: u32,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GjkEpaResult {
    pub intersect: bool,
    pub contact: ContactPoint,
    pub gjk_iterations: usize,
    pub epa_iterations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PenaltyContactModel {
    pub normal_stiffness: f64,
    pub normal_damping: f64,
    pub tangential_damping: f64,
    pub friction: f64,
}

impl Default for PenaltyContactModel {
    fn default() -> Self {
        Self {
            normal_stiffness: 1e5,
            normal_damping: 100.0,
            tangential_damping: 50.0,
            friction: 0.5,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzMindlinContactModel {
    pub normal_stiffness: f64,
    pub normal_damping: f64,
    pub tangential_stiffness: f64,
    pub tangential_damping: f64,
    pub friction: f64,
    pub rolling_resistance: f64,
    pub cohesion_force: f64,
}

impl Default for HertzMindlinContactModel {
    fn default() -> Self {
        Self {
            normal_stiffness: 1e5,
            normal_damping: 20.0,
            tangential_stiffness: 5e4,
            tangential_damping: 20.0,
            friction: 0.5,
            rolling_resistance: 0.01,
            cohesion_force: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpulseContactModel {
    pub restitution: f64,
    pub friction: f64,
    pub baumgarte: f64,
    pub iterations: usize,
}

impl Default for ImpulseContactModel {
    fn default() -> Self {
        Self {
            restitution: 0.0,
            friction: 0.5,
            baumgarte: 0.15,
            iterations: 12,
        }
    }
}

/// Force-threshold bonded-particle law with progressive damage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondedParticleModel {
    pub normal_stiffness: f64,
    pub shear_stiffness: f64,
    pub normal_damping: f64,
    pub shear_damping: f64,
    pub tensile_failure_force: f64,
    pub shear_failure_force: f64,
    pub damage_onset_ratio: f64,
}

impl Default for BondedParticleModel {
    fn default() -> Self {
        Self {
            normal_stiffness: 1e5,
            shear_stiffness: 5e4,
            normal_damping: 50.0,
            shear_damping: 25.0,
            tensile_failure_force: 1e3,
            shear_failure_force: 1e3,
            damage_onset_ratio: 0.7,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ParticleBond {
    pub body_a: usize,
    pub body_b: usize,
    pub rest_length: f64,
    pub model: BondedParticleModel,
    pub tangential_displacement: Vec3,
    pub damage: f64,
    pub broken: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BondUpdateStats {
    pub active: usize,
    pub damaged: usize,
    pub broken_total: usize,
    pub broken_this_step: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PersistentContactState {
    pub contact: ContactPoint,
    pub tangential_displacement: Vec3,
    pub age: usize,
}

/// Persistent pair/feature contact state with frame-safe tangential history.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ContactManifoldCache {
    states: Vec<PersistentContactState>,
}

fn same_key(a: &ContactPoint, b: &ContactPoint) -> bool {
    a.body_a == b.body_a && a.body_b == b.body_b && a.feature_id == b.feature_id
}

impl ContactManifoldCache {
    pub fn update(&mut self, contacts: &[ContactPoint]) {
        let next = contacts
            .iter()
            .map(|c| {
                let found = self.states.iter().find(|s| same_key(&s.contact, c));
                let mut s = found.copied().unwrap_or_default();
                s.contact = *c;
                s.tangential_displacement -= c.normal * dot(s.tangential_displacement, c.normal);
                s.age = found.map_or(1, |f| f.age + 1);
                s
            })
            .collect();
        self.states = next;
    }
    pub fn clear(&mut self) {
        self.states.clear();
    }
    pub fn states(&self) -> &[PersistentContactState] {
        &self.states
    }
    pub fn states_mut(&mut self) -> &mut Vec<PersistentContactState> {
        &mut self.states
    }
}

pub fn sphere_aabb(b: &RigidBody, s: &SphereShape) -> Result<Aabb> {
    require(s.radius > 0.0, "sphere radius must be positive")?;
    let r = Vec3::new(s.radius, s.radius, s.radius);
    Ok(Aabb {
        minimum: b.state.position - r,
        maximum: b.state.position + r,
    })
}

fn vmin(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z))
}

fn vmax(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z))
}

pub fn convex_aabb(b: &RigidBody, shape: &ConvexHullShape) -> Result<Aabb> {
    require(
        !shape.vertices.is_empty(),
        "convex hull must contain vertices",
    )?;
    let first = b.state.position + rotate(b.state.orientation, shape.vertices[0])?;
    let mut bx = Aabb {
        minimum: first,
        maximum: first,
    };
    for v in &shape.vertices {
        let p = b.state.position + rotate(b.state.orientation, *v)?;
        bx.minimum = vmin(bx.minimum, p);
        bx.maximum = vmax(bx.maximum, p);
    }
    Ok(bx)
}

fn sphere_boxes(bodies: &[RigidBody], spheres: &[SphereShape]) -> Result<Vec<Aabb>> {
    spheres
        .iter()
        .map(|s| sphere_aabb(body(bodies, s.body)?, s))
        .collect()
}

pub fn broad_phase_sweep_and_prune(
    bodies: &[RigidBody],
    spheres: &[SphereShape],
) -> Result<Vec<(usize, usize)>> {
    let boxes = sphere_boxes(bodies, spheres)?;
    let mut e: Vec<(usize, Aabb)> = boxes.into_iter().enumerate().collect();
    e.sort_by(|a, b| a.1.minimum.x.total_cmp(&b.1.minimum.x));
    let mut out = vec![];
    for i in 0..e.len() {
        let mut j = i + 1;
        while j < e.len() && e[j].1.minimum.x <= e[i].1.maximum.x {
            let (a, b) = (e[i].1, e[j].1);
            let separated = a.maximum.y < b.minimum.y
                || b.maximum.y < a.minimum.y
                || a.maximum.z < b.minimum.z
                || b.maximum.z < a.minimum.z;
            if !separated && spheres[e[i].0].body != spheres[e[j].0].body {
                out.push((e[i].0, e[j].0));
            }
            j += 1;
        }
    }
    Ok(out)
}

fn overlaps(a: &Aabb, b: &Aabb) -> bool {
    !(a.maximum.x < b.minimum.x
        || b.maximum.x < a.minimum.x
        || a.maximum.y < b.minimum.y
        || b.maximum.y < a.minimum.y
        || a.maximum.z < b.minimum.z
        || b.maximum.z < a.minimum.z)
}

struct BvhNode {
    bx: Aabb,
    children: Option<(usize, usize)>,
    begin: usize,
    end: usize,
}

pub fn broad_phase_bvh(
    bodies: &[RigidBody],
    spheres: &[SphereShape],
) -> Result<Vec<(usize, usize)>> {
    if spheres.len() < 2 {
        return Ok(vec![]);
    }
    let boxes = sphere_boxes(bodies, spheres)?;
    let mut order: Vec<usize> = (0..spheres.len()).collect();
    let mut nodes: Vec<BvhNode> = Vec::with_capacity(2 * spheres.len());
    // Iterative build: (node id, begin, end).
    let make = |order: &[usize], begin: usize, end: usize| {
        let mut bx = boxes[order[begin]];
        for i in begin + 1..end {
            bx = Aabb {
                minimum: vmin(bx.minimum, boxes[order[i]].minimum),
                maximum: vmax(bx.maximum, boxes[order[i]].maximum),
            };
        }
        bx
    };
    nodes.push(BvhNode {
        bx: make(&order, 0, order.len()),
        children: None,
        begin: 0,
        end: order.len(),
    });
    let mut stack = vec![0usize];
    while let Some(id) = stack.pop() {
        let (begin, end, bx) = (nodes[id].begin, nodes[id].end, nodes[id].bx);
        if end - begin <= 2 {
            continue;
        }
        let ext = bx.maximum - bx.minimum;
        let axis = if ext.y > ext.x && ext.y >= ext.z {
            1
        } else if ext.z > ext.x && ext.z > ext.y {
            2
        } else {
            0
        };
        let center = |si: usize| {
            let b = &boxes[si];
            match axis {
                0 => 0.5 * (b.minimum.x + b.maximum.x),
                1 => 0.5 * (b.minimum.y + b.maximum.y),
                _ => 0.5 * (b.minimum.z + b.maximum.z),
            }
        };
        let mid = begin + (end - begin) / 2;
        order[begin..end]
            .select_nth_unstable_by(mid - begin, |a, b| center(*a).total_cmp(&center(*b)));
        let l = nodes.len();
        nodes.push(BvhNode {
            bx: make(&order, begin, mid),
            children: None,
            begin,
            end: mid,
        });
        nodes.push(BvhNode {
            bx: make(&order, mid, end),
            children: None,
            begin: mid,
            end,
        });
        nodes[id].children = Some((l, l + 1));
        stack.push(l + 1);
        stack.push(l);
    }
    let mut out = vec![];
    let mut work = vec![(0usize, 0usize)];
    while let Some((na, nb)) = work.pop() {
        let (a, b) = (&nodes[na], &nodes[nb]);
        if !overlaps(&a.bx, &b.bx) {
            continue;
        }
        match (a.children, b.children) {
            (None, None) => {
                for ia in a.begin..a.end {
                    for ib in b.begin..b.end {
                        let (sa, sb) = (order[ia], order[ib]);
                        // Upstream skips `sa >= sb`, which drops pairs whose lower index sits in
                        // the right subtree; emit each distinct pair in canonical order instead.
                        if sa == sb || (na == nb && sa > sb) || spheres[sa].body == spheres[sb].body
                        {
                            continue;
                        }
                        if overlaps(&boxes[sa], &boxes[sb]) {
                            out.push((sa.min(sb), sa.max(sb)));
                        }
                    }
                }
            }
            (Some((l, r)), _) if na == nb => {
                work.push((l, l));
                work.push((l, r));
                work.push((r, r));
            }
            (ac, bc) => {
                let split_a =
                    bc.is_none() || (ac.is_some() && (a.end - a.begin) >= (b.end - b.begin));
                match (split_a, ac, bc) {
                    (true, Some((l, r)), _) => {
                        work.push((l, nb));
                        work.push((r, nb));
                    }
                    (_, _, Some((l, r))) => {
                        work.push((na, l));
                        work.push((na, r));
                    }
                    _ => {}
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

/// Uniform cell list; every sphere diameter must fit in one cell.
pub fn broad_phase_cell_linked(
    bodies: &[RigidBody],
    spheres: &[SphereShape],
    cell_size: f64,
) -> Result<Vec<(usize, usize)>> {
    require(
        cell_size > 0.0 && cell_size.is_finite(),
        "cell-linked broad-phase cell size must be finite and positive",
    )?;
    if spheres.is_empty() {
        return Ok(vec![]);
    }
    let mut lo = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for s in spheres {
        require(
            s.body < bodies.len() && s.radius > 0.0 && 2.0 * s.radius <= cell_size,
            "cell-linked broad phase requires valid spheres with diameter <= cell size",
        )?;
        let p = bodies[s.body].state.position;
        require(p.is_finite(), "rigid body state must be finite")?;
        lo = vmin(lo, p);
        hi = vmax(hi, p);
    }
    let pad = Vec3::new(cell_size, cell_size, cell_size);
    lo -= pad;
    hi += pad;
    let extent = |a: f64, b: f64| -> Result<usize> {
        let c = ((b - a) / cell_size).ceil();
        require(
            c <= (1u64 << 24) as f64,
            "cell-linked broad-phase grid too large",
        )?;
        Ok((c as usize).max(1))
    };
    let (nx, ny, nz) = (
        extent(lo.x, hi.x)?,
        extent(lo.y, hi.y)?,
        extent(lo.z, hi.z)?,
    );
    let cells = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(nz))
        .filter(|c| *c <= 1 << 24)
        .ok_or_else(|| bad("cell-linked broad-phase grid too large"))?;
    let mut buckets: Vec<Vec<usize>> = vec![vec![]; cells];
    let axis =
        |v: f64, o: f64, n: usize| (((v - o) / cell_size) as i64).clamp(0, n as i64 - 1) as usize;
    for (i, s) in spheres.iter().enumerate() {
        let p = bodies[s.body].state.position;
        buckets[(axis(p.z, lo.z, nz) * ny + axis(p.y, lo.y, ny)) * nx + axis(p.x, lo.x, nx)]
            .push(i);
    }
    let mut pairs = vec![];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                for &i in &buckets[(z * ny + y) * nx + x] {
                    for dz in -1i64..=1 {
                        let zz = z as i64 + dz;
                        if zz < 0 || zz >= nz as i64 {
                            continue;
                        }
                        for dy in -1i64..=1 {
                            let yy = y as i64 + dy;
                            if yy < 0 || yy >= ny as i64 {
                                continue;
                            }
                            for dx in -1i64..=1 {
                                let xx = x as i64 + dx;
                                if xx < 0 || xx >= nx as i64 {
                                    continue;
                                }
                                for &j in
                                    &buckets[(zz as usize * ny + yy as usize) * nx + xx as usize]
                                {
                                    if j > i {
                                        let ai =
                                            sphere_aabb(&bodies[spheres[i].body], &spheres[i])?;
                                        let aj =
                                            sphere_aabb(&bodies[spheres[j].body], &spheres[j])?;
                                        if ai.minimum.x <= aj.maximum.x
                                            && ai.maximum.x >= aj.minimum.x
                                            && ai.minimum.y <= aj.maximum.y
                                            && ai.maximum.y >= aj.minimum.y
                                            && ai.minimum.z <= aj.maximum.z
                                            && ai.maximum.z >= aj.minimum.z
                                        {
                                            pairs.push((i, j));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    Ok(pairs)
}

pub fn sphere_contacts(
    bodies: &[RigidBody],
    spheres: &[SphereShape],
    pairs: &[(usize, usize)],
) -> Result<Vec<ContactPoint>> {
    let mut out = vec![];
    for &(ia, ib) in pairs {
        let (sa, sb) = (
            spheres
                .get(ia)
                .ok_or_else(|| bad("sphere index out of range"))?,
            spheres
                .get(ib)
                .ok_or_else(|| bad("sphere index out of range"))?,
        );
        let (a, b) = (body(bodies, sa.body)?, body(bodies, sb.body)?);
        let d = a.state.position - b.state.position;
        let l = norm(d);
        let pen = sa.radius + sb.radius - l;
        if pen.is_nan() || pen <= 0.0 {
            continue;
        }
        let n = if l > 1e-12 {
            div(d, l)
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        out.push(ContactPoint {
            body_a: Some(sa.body),
            body_b: Some(sb.body),
            point: a.state.position - n * (sa.radius - 0.5 * pen),
            normal: n,
            penetration: pen,
            feature_id: 0,
        });
    }
    Ok(out)
}

pub fn sphere_plane_contacts(
    bodies: &[RigidBody],
    spheres: &[SphereShape],
    planes: &[PlaneShape],
) -> Result<Vec<ContactPoint>> {
    let mut out = vec![];
    for s in spheres {
        let b = body(bodies, s.body)?;
        for (pi, pr) in planes.iter().enumerate() {
            let n = nrm(pr.normal);
            let gap = dot(n, b.state.position) - pr.offset - s.radius;
            if gap < 0.0 {
                out.push(ContactPoint {
                    body_a: Some(s.body),
                    body_b: None,
                    point: b.state.position - n * s.radius,
                    normal: n,
                    penetration: -gap,
                    feature_id: pi as u64 + 1,
                });
            }
        }
    }
    Ok(out)
}

fn closest_point_on_triangle(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = 1.0 / (va + vb + vc);
    a + ab * (vb * denom) + ac * (vc * denom)
}

fn mesh_vertex(bodies: &[RigidBody], mesh: &TriangleMeshShape, local: Vec3) -> Result<Vec3> {
    match mesh.body {
        None => Ok(local),
        Some(i) => body_point(body(bodies, i)?, local),
    }
}

pub fn sphere_triangle_mesh_contacts(
    bodies: &[RigidBody],
    spheres: &[SphereShape],
    mesh: &TriangleMeshShape,
) -> Result<Vec<ContactPoint>> {
    require(
        mesh.body.is_none_or(|b| b < bodies.len()),
        "triangle mesh references invalid body",
    )?;
    if mesh.vertices.is_empty() || mesh.triangles.is_empty() {
        return Ok(vec![]);
    }
    let mut out = vec![];
    for s in spheres {
        require(s.body < bodies.len(), "sphere references invalid body")?;
        if mesh.body == Some(s.body) {
            continue;
        }
        let center = bodies[s.body].state.position;
        for (ti, tri) in mesh.triangles.iter().enumerate() {
            require(
                tri.iter().all(|v| *v < mesh.vertices.len()),
                "triangle mesh index out of range",
            )?;
            let a = mesh_vertex(bodies, mesh, mesh.vertices[tri[0]])?;
            let b = mesh_vertex(bodies, mesh, mesh.vertices[tri[1]])?;
            let c = mesh_vertex(bodies, mesh, mesh.vertices[tri[2]])?;
            let cp = closest_point_on_triangle(center, a, b, c);
            let delta = center - cp;
            let distance = norm(delta);
            if distance.is_nan() || distance >= s.radius {
                continue;
            }
            let normal = if distance > 1e-12 {
                div(delta, distance)
            } else {
                let n = nrm(cross(b - a, c - a));
                let centroid = div(a + b + c, 3.0);
                if dot(n, center - centroid) < 0.0 {
                    -n
                } else {
                    n
                }
            };
            out.push(ContactPoint {
                body_a: Some(s.body),
                body_b: mesh.body,
                point: cp,
                normal,
                penetration: s.radius - distance,
                feature_id: ((mesh.feature_namespace as u64) << 32) | (ti as u64 + 1),
            });
        }
    }
    Ok(out)
}

fn velocity_at(bodies: &[RigidBody], id: Option<usize>, p: Vec3) -> Result<Vec3> {
    match id {
        None => Ok(Vec3::ZERO),
        Some(i) => Ok(body(bodies, i)?.point_velocity(p)),
    }
}

fn relative_velocity(bodies: &[RigidBody], c: &ContactPoint) -> Result<Vec3> {
    Ok(velocity_at(bodies, c.body_a, c.point)? - velocity_at(bodies, c.body_b, c.point)?)
}

fn apply_pair_force(bodies: &mut [RigidBody], c: &ContactPoint, f: Vec3) -> Result<()> {
    if let Some(a) = c.body_a {
        bodies
            .get_mut(a)
            .ok_or_else(|| bad("rigid body index out of range"))?
            .add_force_at_point(f, c.point);
    }
    if let Some(b) = c.body_b {
        bodies
            .get_mut(b)
            .ok_or_else(|| bad("rigid body index out of range"))?
            .add_force_at_point(-f, c.point);
    }
    Ok(())
}

fn rolling_torque(
    bodies: &mut [RigidBody],
    c: &ContactPoint,
    spheres: &[SphereShape],
    m: &HertzMindlinContactModel,
    fn_: f64,
) -> Result<()> {
    let wa = c.body_a.map_or(Ok(Vec3::ZERO), |i| {
        body(bodies, i).map(|b| b.state.angular_velocity)
    })?;
    let wb = c.body_b.map_or(Ok(Vec3::ZERO), |i| {
        body(bodies, i).map(|b| b.state.angular_velocity)
    })?;
    let wr = wa - wb;
    let wmag = norm(wr);
    if wmag > 1e-14 && m.rolling_resistance > 0.0 {
        let radius = |id: Option<usize>| {
            id.and_then(|b| spheres.iter().find(|s| s.body == b))
                .map_or(0.0, |s| s.radius)
        };
        let (ra, rb) = (radius(c.body_a), radius(c.body_b));
        let reff = if rb > 0.0 { ra * rb / (ra + rb) } else { ra };
        let tr = wr * (-m.rolling_resistance * fn_.abs() * reff / wmag);
        if let Some(a) = c.body_a {
            bodies[a].add_torque(tr);
        }
        if let Some(b) = c.body_b {
            bodies[b].add_torque(-tr);
        }
    }
    Ok(())
}

pub fn apply_penalty_contacts(
    bodies: &mut [RigidBody],
    contacts: &[ContactPoint],
    m: &PenaltyContactModel,
) -> Result<()> {
    for c in contacts {
        let rel = relative_velocity(bodies, c)?;
        let vn = dot(rel, c.normal);
        let fn_ = (m.normal_stiffness * c.penetration - m.normal_damping * vn).max(0.0);
        let vt = rel - c.normal * vn;
        let vtm = norm(vt);
        let ft = if vtm > 1e-14 {
            vt * (-(m.friction * fn_).min(m.tangential_damping * vtm) / vtm)
        } else {
            Vec3::ZERO
        };
        apply_pair_force(bodies, c, c.normal * fn_ + ft)?;
    }
    Ok(())
}

pub fn apply_hertz_mindlin_contacts(
    bodies: &mut [RigidBody],
    contacts: &[ContactPoint],
    spheres: &[SphereShape],
    m: &HertzMindlinContactModel,
) -> Result<()> {
    for c in contacts {
        let rel = relative_velocity(bodies, c)?;
        let vn = dot(rel, c.normal);
        let root = c.penetration.max(0.0).sqrt();
        let fn_ = (m.normal_stiffness * c.penetration * root
            - m.normal_damping * root * vn
            - m.cohesion_force)
            .max(-m.cohesion_force);
        let vt = rel - c.normal * vn;
        let vtm = norm(vt);
        let ft = if vtm > 1e-14 {
            vt * (-(m.friction * fn_.abs()).min(m.tangential_damping * root * vtm) / vtm)
        } else {
            Vec3::ZERO
        };
        apply_pair_force(bodies, c, c.normal * fn_ + ft)?;
        rolling_torque(bodies, c, spheres, m, fn_)?;
    }
    Ok(())
}

pub fn apply_history_dependent_mindlin_contacts(
    bodies: &mut [RigidBody],
    contacts: &[ContactPoint],
    spheres: &[SphereShape],
    cache: &mut ContactManifoldCache,
    dt: f64,
    m: &HertzMindlinContactModel,
) -> Result<()> {
    require(
        dt > 0.0 && dt.is_finite(),
        "Mindlin contact dt must be finite and positive",
    )?;
    require(
        m.tangential_stiffness >= 0.0 && m.tangential_stiffness.is_finite(),
        "Mindlin tangential stiffness must be finite and non-negative",
    )?;
    cache.update(contacts);
    for state in cache.states.iter_mut() {
        let c = state.contact;
        let rel = relative_velocity(bodies, &c)?;
        let vn = dot(rel, c.normal);
        let root = c.penetration.max(0.0).sqrt();
        let fn_ = (m.normal_stiffness * c.penetration * root
            - m.normal_damping * root * vn
            - m.cohesion_force)
            .max(-m.cohesion_force);
        let vt = rel - c.normal * vn;
        state.tangential_displacement += vt * dt;
        let td = state.tangential_displacement;
        state.tangential_displacement -= c.normal * dot(td, c.normal);
        let kt = m.tangential_stiffness * root;
        let ct = m.tangential_damping * root;
        let mut ft = state.tangential_displacement * (-kt) + vt * (-ct);
        let max_ft = m.friction * fn_.abs();
        let ftm = norm(ft);
        if ftm > max_ft && ftm > 1e-14 {
            ft = ft * (max_ft / ftm);
            if kt > 1e-14 {
                state.tangential_displacement = -div(ft + vt * ct, kt);
            }
        }
        apply_pair_force(bodies, &c, c.normal * fn_ + ft)?;
        rolling_torque(bodies, &c, spheres, m, fn_)?;
    }
    Ok(())
}

/// `rest_length = 0` uses the current center distance (upstream convention).
pub fn make_particle_bond(
    bodies: &[RigidBody],
    body_a: usize,
    body_b: usize,
    model: BondedParticleModel,
    rest_length: f64,
) -> Result<ParticleBond> {
    require(
        body_a < bodies.len() && body_b < bodies.len() && body_a != body_b,
        "invalid bonded-particle body pair",
    )?;
    let nn = |x: f64| x >= 0.0 && x.is_finite();
    let m = &model;
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
        "invalid bonded-particle constitutive parameters",
    )?;
    let distance = norm(bodies[body_b].state.position - bodies[body_a].state.position);
    let rest = if rest_length == 0.0 {
        distance
    } else {
        rest_length
    };
    require(
        rest > 0.0 && rest.is_finite(),
        "bond rest length must be finite and positive",
    )?;
    Ok(ParticleBond {
        body_a,
        body_b,
        rest_length: rest,
        model,
        ..Default::default()
    })
}

pub fn apply_particle_bonds(
    bodies: &mut [RigidBody],
    bonds: &mut [ParticleBond],
    dt: f64,
) -> Result<BondUpdateStats> {
    require(
        dt > 0.0 && dt.is_finite(),
        "bond update dt must be finite and positive",
    )?;
    let mut stats = BondUpdateStats::default();
    for bond in bonds.iter_mut() {
        require(
            bond.body_a < bodies.len() && bond.body_b < bodies.len() && bond.body_a != bond.body_b,
            "bond references invalid body",
        )?;
        if bond.broken {
            stats.broken_total += 1;
            continue;
        }
        stats.active += 1;
        let (a, b) = (bodies[bond.body_a], bodies[bond.body_b]);
        let delta = b.state.position - a.state.position;
        let distance = norm(delta);
        if !(distance > 1e-15 && distance.is_finite()) {
            return Err(SolverError::NumericalFailure(
                "bonded particle centers are coincident/non-finite".into(),
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
        let normal_elastic = m.normal_stiffness * (distance - bond.rest_length);
        let shear_elastic = bond.tangential_displacement * m.shear_stiffness;
        let ratio = (normal_elastic.max(0.0) / m.tensile_failure_force)
            .max(norm(shear_elastic) / m.shear_failure_force);
        if ratio > m.damage_onset_ratio {
            let target =
                ((ratio - m.damage_onset_ratio) / (1.0 - m.damage_onset_ratio)).clamp(0.0, 1.0);
            bond.damage = bond.damage.max(target);
        }
        if ratio >= 1.0 || bond.damage >= 1.0 {
            bond.damage = 1.0;
            bond.broken = true;
            stats.active -= 1;
            stats.broken_total += 1;
            stats.broken_this_step += 1;
            continue;
        }
        if bond.damage > 0.0 {
            stats.damaged += 1;
        }
        let intact = 1.0 - bond.damage;
        let f = n * (intact * normal_elastic + m.normal_damping * vn)
            + (shear_elastic * intact + vt * m.shear_damping);
        if !bodies[bond.body_a].fixed {
            bodies[bond.body_a].add_force(f);
        }
        if !bodies[bond.body_b].fixed {
            bodies[bond.body_b].add_force(-f);
        }
    }
    Ok(stats)
}

fn effective_mass(bodies: &[RigidBody], c: &ContactPoint, dir: Vec3) -> Result<f64> {
    let mut k = 0.0;
    for (id, sign) in [(c.body_a, 1.0), (c.body_b, -1.0)] {
        if let Some(i) = id {
            let b = body(bodies, i)?;
            if b.fixed {
                continue;
            }
            let ang = cross(c.point - b.state.position, dir * sign);
            k += b.inverse_mass() * dot(dir, dir) + dot(ang, b.world_inverse_inertia_mul(ang)?);
        }
    }
    Ok(k)
}

fn apply_contact_impulse(bodies: &mut [RigidBody], c: &ContactPoint, j: Vec3) -> Result<()> {
    if let Some(a) = c.body_a {
        bodies
            .get_mut(a)
            .ok_or_else(|| bad("rigid body index out of range"))?
            .apply_impulse(j, c.point)?;
    }
    if let Some(b) = c.body_b {
        bodies
            .get_mut(b)
            .ok_or_else(|| bad("rigid body index out of range"))?
            .apply_impulse(-j, c.point)?;
    }
    Ok(())
}

/// Sequential impulses with restitution, Baumgarte bias and Coulomb friction.
pub fn resolve_impulse_contacts(
    bodies: &mut [RigidBody],
    contacts: &[ContactPoint],
    dt: f64,
    m: &ImpulseContactModel,
) -> Result<()> {
    require(dt > 0.0, "contact dt must be positive")?;
    let mut acc = vec![0.0; contacts.len()];
    let mut target = vec![0.0; contacts.len()];
    for (ci, c) in contacts.iter().enumerate() {
        let vn = dot(relative_velocity(bodies, c)?, c.normal);
        target[ci] = ((m.baumgarte / dt) * c.penetration).max(-m.restitution * vn.min(0.0));
    }
    for _ in 0..m.iterations {
        for (ci, c) in contacts.iter().enumerate() {
            let vn = dot(relative_velocity(bodies, c)?, c.normal);
            let k = effective_mass(bodies, c, c.normal)?;
            if k <= 1e-18 {
                continue;
            }
            let old = acc[ci];
            let now = (old + (target[ci] - vn) / k).max(0.0);
            acc[ci] = now;
            apply_contact_impulse(bodies, c, c.normal * (now - old))?;
            let rel2 = relative_velocity(bodies, c)?;
            let vt = rel2 - c.normal * dot(rel2, c.normal);
            let vm = norm(vt);
            if vm > 1e-14 {
                let t = div(vt, vm);
                let kt = effective_mass(bodies, c, t)?;
                if kt > 1e-18 {
                    let jt = (-m.friction * now).max((m.friction * now).min(-vm / kt));
                    apply_contact_impulse(bodies, c, t * jt)?;
                }
            }
        }
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct SupportPoint {
    point: Vec3,
    point_a: Vec3,
    point_b: Vec3,
}

fn convex_support_world(b: &RigidBody, shape: &ConvexHullShape, direction: Vec3) -> Result<Vec3> {
    require(
        !shape.vertices.is_empty(),
        "convex hull must contain vertices",
    )?;
    let ld = inverse_rotate(b.state.orientation, direction)?;
    let mut best = shape.vertices[0];
    let mut best_dot = dot(best, ld);
    for v in &shape.vertices {
        let d = dot(*v, ld);
        if d > best_dot {
            best = *v;
            best_dot = d;
        }
    }
    Ok(b.state.position + rotate(b.state.orientation, best)?)
}

fn minkowski_support(
    bodies: &[RigidBody],
    a: &ConvexHullShape,
    b: &ConvexHullShape,
    mut d: Vec3,
) -> Result<SupportPoint> {
    if norm_squared(d) < 1e-28 {
        d = Vec3::new(1.0, 0.0, 0.0);
    }
    let pa = convex_support_world(body(bodies, a.body)?, a, d)?;
    let pb = convex_support_world(body(bodies, b.body)?, b, -d)?;
    Ok(SupportPoint {
        point: pa - pb,
        point_a: pa,
        point_b: pb,
    })
}

fn triple_cross(a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    cross(cross(a, b), c)
}

fn same_dir(a: Vec3, b: Vec3) -> bool {
    dot(a, b) > 0.0
}

fn update_simplex(s: &mut Vec<SupportPoint>, d: &mut Vec3) -> bool {
    let a = s[s.len() - 1];
    let ao = -a.point;
    match s.len() {
        2 => {
            let b = s[0];
            let ab = b.point - a.point;
            if same_dir(ab, ao) {
                *d = triple_cross(ab, ao, ab);
                if norm_squared(*d) < 1e-28 {
                    *d = any_perpendicular(ab);
                }
            } else {
                *s = vec![a];
                *d = ao;
            }
            false
        }
        3 => {
            let (b, c) = (s[1], s[0]);
            let (ab, ac) = (b.point - a.point, c.point - a.point);
            let abc = cross(ab, ac);
            if same_dir(cross(abc, ac), ao) {
                if same_dir(ac, ao) {
                    *s = vec![c, a];
                    *d = triple_cross(ac, ao, ac);
                } else if same_dir(ab, ao) {
                    *s = vec![b, a];
                    *d = triple_cross(ab, ao, ab);
                } else {
                    *s = vec![a];
                    *d = ao;
                }
                return false;
            }
            if same_dir(cross(ab, abc), ao) {
                if same_dir(ab, ao) {
                    *s = vec![b, a];
                    *d = triple_cross(ab, ao, ab);
                } else {
                    *s = vec![a];
                    *d = ao;
                }
                return false;
            }
            if same_dir(abc, ao) {
                *d = abc;
            } else {
                *s = vec![b, c, a];
                *d = -abc;
            }
            false
        }
        4 => {
            let (b, c, dd) = (s[2], s[1], s[0]);
            let faces = [
                (b, c, dd, [c, b, a]),
                (c, dd, b, [dd, c, a]),
                (dd, b, c, [b, dd, a]),
            ];
            for (p1, p2, opposite, kept) in faces {
                let mut n = cross(p1.point - a.point, p2.point - a.point);
                if dot(n, opposite.point - a.point) > 0.0 {
                    n = -n;
                }
                if dot(n, ao) > 0.0 {
                    *s = kept.to_vec();
                    *d = n;
                    return false;
                }
            }
            true
        }
        _ => {
            *d = ao;
            false
        }
    }
}

#[derive(Clone, Copy)]
struct EpaFace {
    a: usize,
    b: usize,
    c: usize,
    normal: Vec3,
    distance: f64,
}

fn make_face(v: &[SupportPoint], ia: usize, ib: usize, ic: usize) -> EpaFace {
    let mut f = EpaFace {
        a: ia,
        b: ib,
        c: ic,
        normal: Vec3::ZERO,
        distance: 0.0,
    };
    let n = cross(v[ib].point - v[ia].point, v[ic].point - v[ia].point);
    let l = norm(n);
    if l < 1e-14 {
        return EpaFace {
            a: ia,
            b: ib,
            c: ic,
            normal: Vec3::new(1.0, 0.0, 0.0),
            distance: f64::INFINITY,
        };
    }
    let mut n = div(n, l);
    let mut d = dot(n, v[ia].point);
    if d < 0.0 {
        std::mem::swap(&mut f.b, &mut f.c);
        n = -n;
        d = -d;
    }
    f.normal = n;
    f.distance = d;
    f
}

fn barycentric(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> [f64; 3] {
    let (v0, v1, v2) = (b - a, c - a, p - a);
    let (d00, d01, d11, d20, d21) = (
        dot(v0, v0),
        dot(v0, v1),
        dot(v1, v1),
        dot(v2, v0),
        dot(v2, v1),
    );
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-24 {
        return [1.0, 0.0, 0.0];
    }
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    [1.0 - v - w, v, w]
}

fn closest_face(faces: &[EpaFace]) -> Option<EpaFace> {
    let mut best: Option<EpaFace> = None;
    for f in faces {
        if best.is_none_or(|b| f.distance < b.distance) {
            best = Some(*f);
        }
    }
    best
}

/// GJK intersection test followed by EPA penetration depth for two convex hulls (upstream
/// defaults 32/64 iterations and tolerance 1e-9).
pub fn convex_gjk_epa_contact(
    bodies: &[RigidBody],
    a: &ConvexHullShape,
    b: &ConvexHullShape,
    max_gjk_iterations: usize,
    max_epa_iterations: usize,
    tolerance: f64,
) -> Result<GjkEpaResult> {
    require(
        a.body < bodies.len() && b.body < bodies.len() && a.body != b.body,
        "invalid convex hull body attachment",
    )?;
    require(
        a.vertices.len() >= 4 && b.vertices.len() >= 4,
        "GJK/EPA convex hull requires at least four vertices per shape",
    )?;
    require(
        max_gjk_iterations > 0 && max_epa_iterations > 0 && tolerance > 0.0,
        "invalid GJK/EPA iteration/tolerance settings",
    )?;
    let mut d = bodies[a.body].state.position - bodies[b.body].state.position;
    if norm_squared(d) < 1e-24 {
        d = Vec3::new(1.0, 0.0, 0.0);
    }
    let mut simplex = vec![minkowski_support(bodies, a, b, d)?];
    d = -simplex[0].point;
    let mut gjk = 0;
    let mut intersect = false;
    while gjk < max_gjk_iterations {
        if norm_squared(d) < tolerance * tolerance {
            intersect = true;
            break;
        }
        let p = minkowski_support(bodies, a, b, d)?;
        if dot(p.point, d) < 0.0 {
            return Ok(GjkEpaResult {
                gjk_iterations: gjk + 1,
                ..Default::default()
            });
        }
        simplex.push(p);
        if update_simplex(&mut simplex, &mut d) {
            intersect = true;
            gjk += 1;
            break;
        }
        gjk += 1;
    }
    if !intersect || simplex.len() < 4 {
        return Ok(GjkEpaResult {
            gjk_iterations: gjk,
            ..Default::default()
        });
    }
    let mut v = simplex;
    let mut faces = vec![
        make_face(&v, 0, 1, 2),
        make_face(&v, 0, 3, 1),
        make_face(&v, 0, 2, 3),
        make_face(&v, 1, 3, 2),
    ];
    let finish = |v: &[SupportPoint], f: &EpaFace, epa: usize| {
        let projected = f.normal * f.distance;
        let mut w = barycentric(projected, v[f.a].point, v[f.b].point, v[f.c].point);
        for wi in w.iter_mut() {
            *wi = wi.clamp(0.0, 1.0);
        }
        let ws = w[0] + w[1] + w[2];
        if ws > 1e-15 {
            w = [w[0] / ws, w[1] / ws, w[2] / ws];
        }
        let pa = v[f.a].point_a * w[0] + v[f.b].point_a * w[1] + v[f.c].point_a * w[2];
        let pb = v[f.a].point_b * w[0] + v[f.b].point_b * w[1] + v[f.c].point_b * w[2];
        GjkEpaResult {
            intersect: true,
            contact: ContactPoint {
                body_a: Some(a.body),
                body_b: Some(b.body),
                point: (pa + pb) * 0.5,
                normal: -f.normal,
                penetration: f.distance.max(0.0),
                feature_id: 0,
            },
            gjk_iterations: gjk,
            epa_iterations: epa,
        }
    };
    for iter in 0..max_epa_iterations {
        let Some(best) = closest_face(&faces) else {
            break;
        };
        if !best.distance.is_finite() {
            break;
        }
        let p = minkowski_support(bodies, a, b, best.normal)?;
        if dot(p.point, best.normal) - best.distance <= tolerance {
            return Ok(finish(&v, &best, iter + 1));
        }
        let new_index = v.len();
        v.push(p);
        let mut horizon: Vec<(usize, usize)> = vec![];
        let mut add_edge = |x: usize, y: usize| {
            if let Some(pos) = horizon.iter().position(|e| *e == (y, x)) {
                horizon.remove(pos);
            } else {
                horizon.push((x, y));
            }
        };
        let mut kept = Vec::with_capacity(faces.len() + 8);
        for f in &faces {
            if dot(f.normal, p.point - v[f.a].point) > tolerance {
                add_edge(f.a, f.b);
                add_edge(f.b, f.c);
                add_edge(f.c, f.a);
            } else {
                kept.push(*f);
            }
        }
        if horizon.is_empty() {
            return Ok(finish(&v, &best, iter + 1));
        }
        for (x, y) in horizon {
            let f = make_face(&v, x, y, new_index);
            if f.distance.is_finite() {
                kept.push(f);
            }
        }
        faces = kept;
    }
    if let Some(best) = closest_face(&faces) {
        if best.distance.is_finite() {
            return Ok(finish(&v, &best, max_epa_iterations));
        }
    }
    Ok(GjkEpaResult {
        gjk_iterations: gjk,
        epa_iterations: max_epa_iterations,
        ..Default::default()
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ContactMethod {
    #[default]
    SmoothPenalty,
    NonsmoothImpulse,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultibodyConfig {
    pub dt: f64,
    pub gravity: Vec3,
    pub integrator: RigidBodyIntegrator,
    pub contact_method: ContactMethod,
    pub constraint_iterations: usize,
    pub constraint_baumgarte: f64,
    pub penalty: PenaltyContactModel,
    pub impulse: ImpulseContactModel,
}

impl Default for MultibodyConfig {
    fn default() -> Self {
        Self {
            dt: 1e-3,
            gravity: Vec3::new(0.0, -9.81, 0.0),
            integrator: RigidBodyIntegrator::SemiImplicitEuler,
            contact_method: ContactMethod::SmoothPenalty,
            constraint_iterations: 16,
            constraint_baumgarte: 0.2,
            penalty: PenaltyContactModel::default(),
            impulse: ImpulseContactModel::default(),
        }
    }
}

fn gather_contacts(
    bodies: &[RigidBody],
    spheres: &[SphereShape],
    planes: &[PlaneShape],
    meshes: &[TriangleMeshShape],
) -> Result<Vec<ContactPoint>> {
    let pairs = broad_phase_sweep_and_prune(bodies, spheres)?;
    let mut contacts = sphere_contacts(bodies, spheres, &pairs)?;
    contacts.extend(sphere_plane_contacts(bodies, spheres, planes)?);
    for mesh in meshes {
        contacts.extend(sphere_triangle_mesh_contacts(bodies, spheres, mesh)?);
    }
    Ok(contacts)
}

/// Rigid bodies with spheres, planes, meshes and joints; each step is transactional.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidBodySystem {
    config: MultibodyConfig,
    bodies: Vec<RigidBody>,
    spheres: Vec<SphereShape>,
    planes: Vec<PlaneShape>,
    meshes: Vec<TriangleMeshShape>,
    joints: Vec<Joint>,
    time: f64,
}

impl RigidBodySystem {
    pub fn new(config: MultibodyConfig) -> Result<Self> {
        require(config.dt > 0.0, "multibody dt must be positive")?;
        Ok(Self {
            config,
            bodies: vec![],
            spheres: vec![],
            planes: vec![],
            meshes: vec![],
            joints: vec![],
            time: 0.0,
        })
    }
    pub fn add_body(&mut self, mut b: RigidBody) -> Result<usize> {
        b.validate()?;
        b.state.orientation = normalized_quaternion(b.state.orientation)?;
        self.bodies.push(b);
        Ok(self.bodies.len() - 1)
    }
    pub fn add_sphere(&mut self, body: usize, radius: f64) -> Result<()> {
        require(
            body < self.bodies.len() && radius > 0.0,
            "invalid sphere attachment",
        )?;
        self.spheres.push(SphereShape { body, radius });
        Ok(())
    }
    pub fn add_plane(&mut self, mut p: PlaneShape) {
        p.normal = nrm(p.normal);
        self.planes.push(p);
    }
    pub fn add_triangle_mesh(&mut self, mut mesh: TriangleMeshShape) -> Result<()> {
        require(
            mesh.body.is_none_or(|b| b < self.bodies.len()),
            "invalid triangle mesh body attachment",
        )?;
        if mesh.feature_namespace == 0 {
            mesh.feature_namespace = self.meshes.len() as u32 + 1;
        }
        self.meshes.push(mesh);
        Ok(())
    }
    pub fn add_joint(&mut self, j: Joint) {
        self.joints.push(j);
    }
    pub fn step(&mut self) -> Result<()> {
        let mut bodies = self.bodies.clone();
        let c = self.config;
        for b in bodies.iter_mut().filter(|b| !b.fixed) {
            b.add_force(c.gravity * b.mass);
        }
        let contacts = gather_contacts(&bodies, &self.spheres, &self.planes, &self.meshes)?;
        if c.contact_method == ContactMethod::SmoothPenalty {
            apply_penalty_contacts(&mut bodies, &contacts, &c.penalty)?;
        }
        for b in bodies.iter_mut() {
            integrate_rigid_body(b, c.dt, c.integrator)?;
        }
        let mut rows = vec![];
        for j in &self.joints {
            rows.extend(build_constraint_rows(
                &bodies,
                j,
                c.dt,
                c.constraint_baumgarte,
            )?);
        }
        solve_constraint_rows(&mut bodies, &mut rows, c.constraint_iterations)?;
        if c.contact_method == ContactMethod::NonsmoothImpulse {
            resolve_impulse_contacts(&mut bodies, &contacts, c.dt, &c.impulse)?;
        }
        self.bodies = bodies;
        self.time += c.dt;
        Ok(())
    }
    pub fn run(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }
    pub fn bodies(&self) -> &[RigidBody] {
        &self.bodies
    }
    pub fn bodies_mut(&mut self) -> &mut [RigidBody] {
        &mut self.bodies
    }
    pub fn spheres(&self) -> &[SphereShape] {
        &self.spheres
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn total_kinetic_energy(&self) -> Result<f64> {
        let mut e = 0.0;
        for b in self.bodies.iter().filter(|b| !b.fixed) {
            e += 0.5 * b.mass * norm_squared(b.state.linear_velocity);
            let w = inverse_rotate(b.state.orientation, b.state.angular_velocity)?;
            let i = b.inertia_diagonal;
            e += 0.5 * (i.x * w.x * w.x + i.y * w.y * w.y + i.z * w.z * w.z);
        }
        Ok(e)
    }
    pub fn linear_momentum(&self) -> Vec3 {
        let mut p = Vec3::ZERO;
        for b in self.bodies.iter().filter(|b| !b.fixed) {
            p += b.state.linear_velocity * b.mass;
        }
        p
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DemConfig {
    pub dt: f64,
    pub gravity: Vec3,
    pub contact: HertzMindlinContactModel,
    pub history_dependent_tangential: bool,
}

impl Default for DemConfig {
    fn default() -> Self {
        Self {
            dt: 1e-4,
            gravity: Vec3::new(0.0, -9.81, 0.0),
            contact: HertzMindlinContactModel::default(),
            history_dependent_tangential: true,
        }
    }
}

/// Explicit spherical DEM with Hertz-Mindlin contact, optional tangential history and bonds;
/// each step is transactional.
#[derive(Clone, Debug, PartialEq)]
pub struct ExplicitDemSystem {
    config: DemConfig,
    bodies: Vec<RigidBody>,
    spheres: Vec<SphereShape>,
    planes: Vec<PlaneShape>,
    meshes: Vec<TriangleMeshShape>,
    history: ContactManifoldCache,
    bonds: Vec<ParticleBond>,
    last_bond_stats: BondUpdateStats,
    time: f64,
}

impl ExplicitDemSystem {
    pub fn new(config: DemConfig) -> Result<Self> {
        require(config.dt > 0.0, "DEM dt must be positive")?;
        Ok(Self {
            config,
            bodies: vec![],
            spheres: vec![],
            planes: vec![],
            meshes: vec![],
            history: ContactManifoldCache::default(),
            bonds: vec![],
            last_bond_stats: BondUpdateStats::default(),
            time: 0.0,
        })
    }
    pub fn add_particle(
        &mut self,
        position: Vec3,
        radius: f64,
        density: f64,
        velocity: Vec3,
    ) -> Result<usize> {
        require(
            radius > 0.0 && density > 0.0,
            "DEM radius/density must be positive",
        )?;
        let m = density * (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
        let inertia = 0.4 * m * radius * radius;
        let id = self.bodies.len();
        self.bodies.push(RigidBody {
            state: RigidBodyState {
                position,
                linear_velocity: velocity,
                ..Default::default()
            },
            mass: m,
            inertia_diagonal: Vec3::new(inertia, inertia, inertia),
            ..Default::default()
        });
        self.spheres.push(SphereShape { body: id, radius });
        Ok(id)
    }
    pub fn add_plane(&mut self, mut p: PlaneShape) {
        p.normal = nrm(p.normal);
        self.planes.push(p);
    }
    pub fn add_triangle_mesh(&mut self, mut mesh: TriangleMeshShape) -> Result<()> {
        require(
            mesh.body.is_none_or(|b| b < self.bodies.len()),
            "invalid DEM triangle mesh body attachment",
        )?;
        if mesh.feature_namespace == 0 {
            mesh.feature_namespace = self.meshes.len() as u32 + 1;
        }
        self.meshes.push(mesh);
        Ok(())
    }
    pub fn add_bond(
        &mut self,
        body_a: usize,
        body_b: usize,
        model: BondedParticleModel,
        rest_length: f64,
    ) -> Result<usize> {
        self.bonds.push(make_particle_bond(
            &self.bodies,
            body_a,
            body_b,
            model,
            rest_length,
        )?);
        Ok(self.bonds.len() - 1)
    }
    pub fn step(&mut self) -> Result<()> {
        let c = self.config;
        let mut bodies = self.bodies.clone();
        let mut bonds = self.bonds.clone();
        let mut history = self.history.clone();
        for b in bodies.iter_mut() {
            b.add_force(c.gravity * b.mass);
        }
        let stats = apply_particle_bonds(&mut bodies, &mut bonds, c.dt)?;
        let contacts = gather_contacts(&bodies, &self.spheres, &self.planes, &self.meshes)?;
        if c.history_dependent_tangential {
            apply_history_dependent_mindlin_contacts(
                &mut bodies,
                &contacts,
                &self.spheres,
                &mut history,
                c.dt,
                &c.contact,
            )?;
        } else {
            history.clear();
            apply_hertz_mindlin_contacts(&mut bodies, &contacts, &self.spheres, &c.contact)?;
        }
        for b in bodies.iter_mut() {
            integrate_rigid_body(b, c.dt, RigidBodyIntegrator::SemiImplicitEuler)?;
        }
        self.bodies = bodies;
        self.bonds = bonds;
        self.history = history;
        self.last_bond_stats = stats;
        self.time += c.dt;
        Ok(())
    }
    pub fn run(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }
    pub fn bodies(&self) -> &[RigidBody] {
        &self.bodies
    }
    pub fn bodies_mut(&mut self) -> &mut [RigidBody] {
        &mut self.bodies
    }
    pub fn spheres(&self) -> &[SphereShape] {
        &self.spheres
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn contact_history(&self) -> &ContactManifoldCache {
        &self.history
    }
    pub fn bonds(&self) -> &[ParticleBond] {
        &self.bonds
    }
    pub fn last_bond_stats(&self) -> &BondUpdateStats {
        &self.last_bond_stats
    }
}

/// Stokes drag on a sphere; its negative is the reaction the fluid must receive.
pub fn stokes_drag_force(
    body: &RigidBody,
    radius: f64,
    fluid_velocity: Vec3,
    dynamic_viscosity: f64,
) -> Result<Vec3> {
    require(
        radius > 0.0 && dynamic_viscosity >= 0.0,
        "invalid Stokes drag parameters",
    )?;
    Ok((fluid_velocity - body.state.linear_velocity)
        * (6.0 * std::f64::consts::PI * dynamic_viscosity * radius))
}
