//! Immersed-boundary point particles on a periodic unit-spacing lattice, ported from upstream
//! `cfd/solvers/lbm/particles` (MIT, cfd_solvers contributors).
use crate::{require, Result};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ParticleVector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImmersedBoundaryParticle {
    pub position: ParticleVector3,
    pub velocity: ParticleVector3,
    pub mass: f64,
    pub response_time: f64,
}

impl Default for ImmersedBoundaryParticle {
    fn default() -> Self {
        Self {
            position: ParticleVector3::default(),
            velocity: ParticleVector3::default(),
            mass: 1.0,
            response_time: 1.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParticleCouplingResult {
    /// Empty unless two-way coupling was requested.
    pub acceleration_x: Vec<f32>,
    pub acceleration_y: Vec<f32>,
    pub acceleration_z: Vec<f32>,
    pub total_reaction_force: ParticleVector3,
}

fn validate(nx: usize, ny: usize, nz: usize, ux: &[f32], uy: &[f32], uz: &[f32]) -> Result<usize> {
    require(
        nx >= 2 && ny >= 2 && nz >= 2,
        "particle coupling grid must be at least 2x2x2",
    )?;
    let cells = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(nz))
        .unwrap_or(usize::MAX);
    require(
        ux.len() == cells && uy.len() == cells && uz.len() == cells,
        "particle coupling velocity field size mismatch",
    )?;
    Ok(cells)
}

fn fold(value: f64, extent: f64) -> f64 {
    let v = value % extent;
    if v < 0.0 {
        v + extent
    } else {
        v
    }
}

fn periodic_wrap(p: &mut ParticleVector3, nx: usize, ny: usize, nz: usize) {
    p.x = fold(p.x, nx as f64);
    p.y = fold(p.y, ny as f64);
    p.z = fold(p.z, nz as f64);
}

/// Eight `(cell index, weight)` pairs of the periodic trilinear kernel.
fn trilinear_kernel(nx: usize, ny: usize, nz: usize, p: ParticleVector3) -> [(usize, f64); 8] {
    let (fx0, fy0, fz0) = (p.x.floor(), p.y.floor(), p.z.floor());
    let (x0, y0, z0) = (fx0 as i64, fy0 as i64, fz0 as i64);
    let (tx, ty, tz) = (p.x - fx0, p.y - fy0, p.z - fz0);
    let mut points = [(0usize, 0.0f64); 8];
    let mut n = 0;
    for dz in 0..=1i64 {
        for dy in 0..=1i64 {
            for dx in 0..=1i64 {
                let wx = if dx == 1 { tx } else { 1.0 - tx };
                let wy = if dy == 1 { ty } else { 1.0 - ty };
                let wz = if dz == 1 { tz } else { 1.0 - tz };
                let x = (x0 + dx).rem_euclid(nx as i64) as usize;
                let y = (y0 + dy).rem_euclid(ny as i64) as usize;
                let z = (z0 + dz).rem_euclid(nz as i64) as usize;
                points[n] = ((z * ny + y) * nx + x, wx * wy * wz);
                n += 1;
            }
        }
    }
    points
}

fn sample(kernel: &[(usize, f64); 8], ux: &[f32], uy: &[f32], uz: &[f32]) -> ParticleVector3 {
    let mut r = ParticleVector3::default();
    for &(i, w) in kernel {
        r.x += w * f64::from(ux[i]);
        r.y += w * f64::from(uy[i]);
        r.z += w * f64::from(uz[i]);
    }
    r
}

/// Periodic trilinear interpolation of a cell-centred lattice velocity.
pub fn interpolate_particle_velocity(
    nx: usize,
    ny: usize,
    nz: usize,
    ux: &[f32],
    uy: &[f32],
    uz: &[f32],
    mut position: ParticleVector3,
) -> Result<ParticleVector3> {
    validate(nx, ny, nz, ux, uy, uz)?;
    periodic_wrap(&mut position, nx, ny, nz);
    Ok(sample(&trilinear_kernel(nx, ny, nz, position), ux, uy, uz))
}

/// Explicit linear-drag particles; with two-way coupling the equal-and-opposite drag is spread
/// back with the same kernel as lattice accelerations (unit cell volume and density).
#[allow(clippy::too_many_arguments)]
pub fn advance_immersed_boundary_particles(
    nx: usize,
    ny: usize,
    nz: usize,
    ux: &[f32],
    uy: &[f32],
    uz: &[f32],
    particles: &mut [ImmersedBoundaryParticle],
    dt: f64,
    two_way_coupling: bool,
    gravity: ParticleVector3,
) -> Result<ParticleCouplingResult> {
    let cells = validate(nx, ny, nz, ux, uy, uz)?;
    require(dt > 0.0, "particle coupling dt must be positive")?;
    require(
        particles
            .iter()
            .all(|p| p.mass > 0.0 && p.response_time > 0.0),
        "particle mass and response time must be positive",
    )?;
    let mut result = ParticleCouplingResult::default();
    if two_way_coupling {
        result.acceleration_x = vec![0.0; cells];
        result.acceleration_y = vec![0.0; cells];
        result.acceleration_z = vec![0.0; cells];
    }
    for p in particles.iter_mut() {
        periodic_wrap(&mut p.position, nx, ny, nz);
        let kernel = trilinear_kernel(nx, ny, nz, p.position);
        let fluid = sample(&kernel, ux, uy, uz);
        let drag = ParticleVector3 {
            x: p.mass * (fluid.x - p.velocity.x) / p.response_time,
            y: p.mass * (fluid.y - p.velocity.y) / p.response_time,
            z: p.mass * (fluid.z - p.velocity.z) / p.response_time,
        };
        p.velocity.x += dt * (drag.x / p.mass + gravity.x);
        p.velocity.y += dt * (drag.y / p.mass + gravity.y);
        p.velocity.z += dt * (drag.z / p.mass + gravity.z);
        p.position.x += dt * p.velocity.x;
        p.position.y += dt * p.velocity.y;
        p.position.z += dt * p.velocity.z;
        periodic_wrap(&mut p.position, nx, ny, nz);
        if two_way_coupling {
            result.total_reaction_force.x -= drag.x;
            result.total_reaction_force.y -= drag.y;
            result.total_reaction_force.z -= drag.z;
            for &(i, w) in &kernel {
                result.acceleration_x[i] += (-drag.x * w) as f32;
                result.acceleration_y[i] += (-drag.y * w) as f32;
                result.acceleration_z[i] += (-drag.z * w) as f32;
            }
        }
    }
    Ok(result)
}
