//! Vector/quaternion helpers with the exact upstream `cfd::multibody` arithmetic (the core
//! [`Quaternion`] type is reused as storage).
pub use crate::core::{Quaternion, Vec3};
use crate::{require, Result};

pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

pub fn norm_squared(a: Vec3) -> f64 {
    dot(a, a)
}

pub fn norm(a: Vec3) -> f64 {
    norm_squared(a).sqrt()
}

/// Unit vector, or zero when the norm is not above `eps` (upstream default 1e-15).
pub fn normalized_or_zero(a: Vec3, eps: f64) -> Vec3 {
    let n = norm(a);
    if n > eps {
        Vec3::new(a.x / n, a.y / n, a.z / n)
    } else {
        Vec3::ZERO
    }
}

pub fn component_mul(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x * b.x, a.y * b.y, a.z * b.z)
}

pub fn quaternion_norm_squared(q: Quaternion) -> f64 {
    q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z
}

pub fn normalized_quaternion(q: Quaternion) -> Result<Quaternion> {
    let n = quaternion_norm_squared(q).sqrt();
    require(
        n > 0.0 && n.is_finite(),
        "invalid zero/non-finite quaternion",
    )?;
    Ok(Quaternion {
        w: q.w / n,
        x: q.x / n,
        y: q.y / n,
        z: q.z / n,
    })
}

pub fn rotate(q: Quaternion, v: Vec3) -> Result<Vec3> {
    let q = normalized_quaternion(q)?;
    let r =
        q * Quaternion {
            w: 0.0,
            x: v.x,
            y: v.y,
            z: v.z,
        } * q.conjugate();
    Ok(Vec3::new(r.x, r.y, r.z))
}

pub fn inverse_rotate(q: Quaternion, v: Vec3) -> Result<Vec3> {
    rotate(normalized_quaternion(q)?.conjugate(), v)
}

pub fn quaternion_from_axis_angle(axis: Vec3, angle: f64) -> Result<Quaternion> {
    let n = normalized_or_zero(axis, 1e-15);
    let h = 0.5 * angle;
    let s = h.sin();
    normalized_quaternion(Quaternion {
        w: h.cos(),
        x: n.x * s,
        y: n.y * s,
        z: n.z * s,
    })
}

/// First-order world-frame update `q += dt/2 * (omega q)`, renormalized.
pub fn integrate_orientation(
    q: &mut Quaternion,
    angular_velocity_world: Vec3,
    dt: f64,
) -> Result<()> {
    let w = angular_velocity_world;
    let dq = Quaternion {
        w: 0.0,
        x: w.x,
        y: w.y,
        z: w.z,
    } * *q;
    *q = normalized_quaternion(Quaternion {
        w: q.w + 0.5 * dt * dq.w,
        x: q.x + 0.5 * dt * dq.x,
        y: q.y + 0.5 * dt * dq.y,
        z: q.z + 0.5 * dt * dq.z,
    })?;
    Ok(())
}

pub fn any_perpendicular(axis: Vec3) -> Vec3 {
    let n = normalized_or_zero(axis, 1e-15);
    let seed = if n.x.abs() < 0.7 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    normalized_or_zero(cross(n, seed), 1e-15)
}

/// Component-wise division (upstream `operator/`; the core `Vec3 / f64` multiplies by the
/// reciprocal instead).
pub fn div(a: Vec3, s: f64) -> Vec3 {
    Vec3::new(a.x / s, a.y / s, a.z / s)
}
