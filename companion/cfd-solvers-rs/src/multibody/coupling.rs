//! Resolved (surface-traction) and unresolved (volume-averaged drag, pressure-gradient and
//! Saffman lift) CFD/DEM coupling primitives.
use super::{cross, div, dot, RigidBody, SphereShape, Vec3};
use crate::{fvm::PolyMesh, require, Result};
use std::f64::consts::PI;

/// Symmetric viscous stress; resolved traction is `-p n + tau n`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SymmetricStress3 {
    pub xx: f64,
    pub yy: f64,
    pub zz: f64,
    pub xy: f64,
    pub xz: f64,
    pub yz: f64,
}

impl SymmetricStress3 {
    pub fn traction(&self, n: Vec3) -> Vec3 {
        Vec3::new(
            self.xx * n.x + self.xy * n.y + self.xz * n.z,
            self.xy * n.x + self.yy * n.y + self.yz * n.z,
            self.xz * n.x + self.yz * n.y + self.zz * n.z,
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidSurfaceSample {
    pub position: Vec3,
    pub normal: Vec3,
    pub area: f64,
}

impl Default for RigidSurfaceSample {
    fn default() -> Self {
        Self {
            position: Vec3::ZERO,
            normal: Vec3::new(1.0, 0.0, 0.0),
            area: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ResolvedSurfaceLoad {
    pub force: Vec3,
    pub torque: Vec3,
}

fn mag(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

fn nearest_cell(mesh: &PolyMesh, p: Vec3) -> Result<usize> {
    require(
        mesh.cell_count() > 0,
        "CFD/DEM coupling requires non-empty mesh",
    )?;
    let (mut best, mut d2) = (0, f64::INFINITY);
    for (c, cell) in mesh.cells().iter().enumerate() {
        let q = cell.center;
        let (dx, dy, dz) = (q.x - p.x, q.y - p.y, q.z - p.z);
        let r2 = dx * dx + dy * dy + dz * dz;
        if r2 < d2 {
            d2 = r2;
            best = c;
        }
    }
    Ok(best)
}

fn sphere_volume(r: f64) -> f64 {
    (4.0 / 3.0) * PI * r * r * r
}

/// Octahedral rule for 6 samples, otherwise antipodal Fibonacci pairs (upstream default 6).
pub fn sphere_surface_quadrature(
    body: &RigidBody,
    radius: f64,
    sample_count: usize,
) -> Result<Vec<RigidSurfaceSample>> {
    require(
        radius > 0.0 && radius.is_finite(),
        "sphere quadrature radius must be positive",
    )?;
    require(
        sample_count >= 6,
        "sphere quadrature requires at least six samples",
    )?;
    require(
        sample_count <= 1 << 24,
        "sphere quadrature sample count too large",
    )?;
    let mut normals = vec![];
    if sample_count == 6 {
        normals = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, -1.0),
        ];
    } else {
        let count = sample_count + sample_count % 2;
        let half = count / 2;
        let golden = PI * (3.0 - 5.0f64.sqrt());
        for i in 0..half {
            let z = (i as f64 + 0.5) / half as f64;
            let rr = (1.0 - z * z).max(0.0).sqrt();
            let phi = golden * i as f64;
            let n = Vec3::new(rr * phi.cos(), rr * phi.sin(), z);
            normals.push(n);
            normals.push(n * -1.0);
        }
    }
    let area = 4.0 * PI * radius * radius / normals.len() as f64;
    Ok(normals
        .into_iter()
        .map(|n| RigidSurfaceSample {
            position: body.state.position + n * radius,
            normal: n,
            area,
        })
        .collect())
}

/// Nearest-cell sampled pressure and (optional, may be empty) viscous stress traction.
pub fn integrate_resolved_surface_load(
    mesh: &PolyMesh,
    cell_pressure: &[f64],
    cell_viscous_stress: &[SymmetricStress3],
    surface: &[RigidSurfaceSample],
    moment_center: Vec3,
) -> Result<ResolvedSurfaceLoad> {
    require(
        cell_pressure.len() == mesh.cell_count(),
        "resolved CFD/body pressure size mismatch",
    )?;
    require(
        cell_viscous_stress.is_empty() || cell_viscous_stress.len() == mesh.cell_count(),
        "resolved CFD/body stress size mismatch",
    )?;
    let mut r = ResolvedSurfaceLoad::default();
    for s in surface {
        require(
            s.area >= 0.0 && s.area.is_finite(),
            "invalid resolved surface area",
        )?;
        let c = nearest_cell(mesh, s.position)?;
        let mut t = s.normal * (-cell_pressure[c]);
        if !cell_viscous_stress.is_empty() {
            t += cell_viscous_stress[c].traction(s.normal);
        }
        let df = t * s.area;
        r.force += df;
        r.torque += cross(s.position - moment_center, df);
    }
    Ok(r)
}

pub fn apply_resolved_surface_load(body: &mut RigidBody, load: &ResolvedSurfaceLoad) {
    body.add_force(load.force);
    body.add_torque(load.torque);
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnresolvedCfdDemConfig {
    pub fluid_density: f64,
    pub dynamic_viscosity: f64,
    pub minimum_void_fraction: f64,
    pub void_drag_exponent: f64,
    pub saffman_coefficient: f64,
    pub include_pressure_gradient: bool,
    pub include_saffman_lift: bool,
}

impl Default for UnresolvedCfdDemConfig {
    fn default() -> Self {
        Self {
            fluid_density: 1.0,
            dynamic_viscosity: 1e-3,
            minimum_void_fraction: 0.05,
            void_drag_exponent: 2.65,
            saffman_coefficient: 1.615,
            include_pressure_gradient: true,
            include_saffman_lift: true,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct UnresolvedCfdDemResult {
    pub solid_volume_fraction: Vec<f64>,
    pub void_fraction: Vec<f64>,
    pub particle_force: Vec<Vec3>,
    /// Reaction force density on the fluid (force / volume).
    pub fluid_momentum_source: Vec<Vec3>,
}

fn schiller_naumann(rel: Vec3, d: f64, rho: f64, mu: f64, eps: f64, exponent: f64) -> Vec3 {
    let speed = mag(rel);
    if speed == 0.0 {
        return Vec3::ZERO;
    }
    let re = rho * speed * d / mu;
    let cd = if re < 1000.0 {
        (24.0 / re.max(1e-30)) * (1.0 + 0.15 * re.powf(0.687))
    } else {
        0.44
    };
    let area = 0.25 * PI * d * d;
    let hindrance = eps.max(1e-12).powf(-exponent);
    rel * (0.5 * cd * rho * area * speed * hindrance)
}

fn saffman(rel: Vec3, vort: Vec3, rho: f64, nu: f64, d: f64, coefficient: f64) -> Vec3 {
    let om = mag(vort);
    if om == 0.0 || mag(rel) == 0.0 {
        return Vec3::ZERO;
    }
    cross(rel, vort) * (coefficient * rho * d * d * (nu * om).sqrt() / om)
}

/// Nearest-cell volume deposition with drag/pressure/lift and the equal-and-opposite fluid
/// source. `pressure_gradient` and `vorticity` may be empty.
#[allow(clippy::too_many_arguments)]
pub fn unresolved_cfd_dem_coupling(
    mesh: &PolyMesh,
    bodies: &[RigidBody],
    spheres: &[SphereShape],
    fluid_velocity: &[Vec3],
    pressure_gradient: &[Vec3],
    vorticity: &[Vec3],
    config: &UnresolvedCfdDemConfig,
) -> Result<UnresolvedCfdDemResult> {
    let n = mesh.cell_count();
    require(
        fluid_velocity.len() == n,
        "CFD/DEM velocity field size mismatch",
    )?;
    require(
        pressure_gradient.is_empty() || pressure_gradient.len() == n,
        "CFD/DEM pressure-gradient size mismatch",
    )?;
    require(
        vorticity.is_empty() || vorticity.len() == n,
        "CFD/DEM vorticity size mismatch",
    )?;
    let c = config;
    require(
        c.fluid_density > 0.0
            && c.dynamic_viscosity > 0.0
            && c.minimum_void_fraction > 0.0
            && c.minimum_void_fraction <= 1.0
            && c.void_drag_exponent >= 0.0,
        "invalid unresolved CFD/DEM configuration",
    )?;
    let mut out = UnresolvedCfdDemResult {
        solid_volume_fraction: vec![0.0; n],
        void_fraction: vec![1.0; n],
        particle_force: vec![Vec3::ZERO; spheres.len()],
        fluid_momentum_source: vec![Vec3::ZERO; n],
    };
    let mut cell_of = Vec::with_capacity(spheres.len());
    for s in spheres {
        require(
            s.body < bodies.len() && s.radius > 0.0,
            "invalid unresolved CFD/DEM sphere",
        )?;
        let cell = nearest_cell(mesh, bodies[s.body].state.position)?;
        cell_of.push(cell);
        out.solid_volume_fraction[cell] += sphere_volume(s.radius) / mesh.cells()[cell].volume;
    }
    for cell in 0..n {
        let alpha = out.solid_volume_fraction[cell].clamp(0.0, 1.0 - c.minimum_void_fraction);
        out.solid_volume_fraction[cell] = alpha;
        out.void_fraction[cell] = c.minimum_void_fraction.max(1.0 - alpha);
    }
    let nu = c.dynamic_viscosity / c.fluid_density;
    for (i, s) in spheres.iter().enumerate() {
        let b = &bodies[s.body];
        let cell = cell_of[i];
        let rel = fluid_velocity[cell] - b.state.linear_velocity;
        let (d, vp) = (2.0 * s.radius, sphere_volume(s.radius));
        let mut force = schiller_naumann(
            rel,
            d,
            c.fluid_density,
            c.dynamic_viscosity,
            out.void_fraction[cell],
            c.void_drag_exponent,
        );
        if c.include_pressure_gradient && !pressure_gradient.is_empty() {
            force -= pressure_gradient[cell] * vp;
        }
        if c.include_saffman_lift && !vorticity.is_empty() {
            force += saffman(
                rel,
                vorticity[cell],
                c.fluid_density,
                nu,
                d,
                c.saffman_coefficient,
            );
        }
        out.particle_force[i] = force;
        out.fluid_momentum_source[cell] -= div(force, mesh.cells()[cell].volume);
    }
    Ok(out)
}

pub fn apply_unresolved_particle_forces(
    bodies: &mut [RigidBody],
    spheres: &[SphereShape],
    sphere_forces: &[Vec3],
) -> Result<()> {
    require(
        sphere_forces.len() == spheres.len(),
        "unresolved particle-force size mismatch",
    )?;
    require(
        spheres.iter().all(|s| s.body < bodies.len()),
        "unresolved sphere body index out of range",
    )?;
    for (s, f) in spheres.iter().zip(sphere_forces) {
        bodies[s.body].add_force(*f);
    }
    Ok(())
}
