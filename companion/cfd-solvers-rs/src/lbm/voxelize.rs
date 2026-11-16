//! Ray-parity voxelization of closed triangle surfaces and rigid surface transforms, ported from
//! upstream `cfd/solvers/lbm/voxelize` (MIT, cfd_solvers contributors).
use crate::{
    core::Vec3,
    io::{Triangle3, TriangleSurface},
    require, Result,
};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct VoxelGrid3D {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub minimum: Vec3,
    pub maximum: Vec3,
}

/// Moller-Trumbore hit of the +x ray from `o`.
fn hit_x(o: Vec3, t: &Triangle3) -> bool {
    let d = Vec3::new(1.0, 0.0, 0.0);
    let e1 = t.b - t.a;
    let e2 = t.c - t.a;
    let h = d.cross(e2);
    let det = e1.dot(h);
    if det.abs() < 1e-12 {
        return false;
    }
    let inv = 1.0 / det;
    let s = o - t.a;
    let u = inv * s.dot(h);
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let q = s.cross(e1);
    let v = inv * d.dot(q);
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    inv * e2.dot(q) > 1e-10
}

/// Inside/outside mask (1 = inside) by +x ray crossing parity from each cell centre, with tiny
/// y/z offsets to avoid grazing edges. Index is `(k * ny + j) * nx + i`.
pub fn voxelize_surface(surface: &TriangleSurface, g: &VoxelGrid3D) -> Result<Vec<u8>> {
    require(
        !surface.triangles.is_empty()
            && g.nx > 0
            && g.ny > 0
            && g.nz > 0
            && g.maximum.x > g.minimum.x
            && g.maximum.y > g.minimum.y
            && g.maximum.z > g.minimum.z,
        "invalid voxelization inputs",
    )?;
    let cells = g.nx.checked_mul(g.ny).and_then(|v| v.checked_mul(g.nz));
    let Some(cells) = cells else {
        return Err(crate::SolverError::InvalidInput(
            "voxel grid too large".into(),
        ));
    };
    let mut mask = vec![0u8; cells];
    let dx = (g.maximum.x - g.minimum.x) / g.nx as f64;
    let dy = (g.maximum.y - g.minimum.y) / g.ny as f64;
    let dz = (g.maximum.z - g.minimum.z) / g.nz as f64;
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let p = Vec3::new(
                    g.minimum.x + (i as f64 + 0.5) * dx,
                    g.minimum.y + (j as f64 + 0.5) * dy + 1e-10 * dx,
                    g.minimum.z + (k as f64 + 0.5) * dz + 2e-10 * dx,
                );
                let hits = surface.triangles.iter().filter(|t| hit_x(p, t)).count();
                mask[(k * g.ny + j) * g.nx + i] = (hits & 1) as u8;
            }
        }
    }
    Ok(mask)
}

/// Rodrigues rotation about `axis` (through the origin; a near-zero axis leaves points), then translation.
pub fn transform_surface(
    surface: &TriangleSurface,
    translation: Vec3,
    rotation_axis: Vec3,
    angle_radians: f64,
) -> TriangleSurface {
    let rotate = |p: Vec3| {
        let n = rotation_axis.dot(rotation_axis).sqrt();
        if n < 1e-30 {
            return p;
        }
        let axis = rotation_axis * (1.0 / n);
        let (c, s) = (angle_radians.cos(), angle_radians.sin());
        p * c + axis.cross(p) * s + axis * (axis.dot(p) * (1.0 - c))
    };
    TriangleSurface {
        triangles: surface
            .triangles
            .iter()
            .map(|t| Triangle3 {
                a: rotate(t.a) + translation,
                b: rotate(t.b) + translation,
                c: rotate(t.c) + translation,
            })
            .collect(),
    }
}
