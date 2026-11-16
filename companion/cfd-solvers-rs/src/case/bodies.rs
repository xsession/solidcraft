//! Arbitrary solid bodies from closed, consistently oriented triangle surfaces (`surface_bodies_v1`).
//!
//! A body is rasterized onto the background grid: a cell is solid when its centre lies inside
//! the surface, decided by a winding number along exact grid-line rays with a consistent
//! tie-breaking rule, so shared edges count once. The result is a staircase representation;
//! no cut cells or curved-wall fidelity is implied. The represented and true volumes are
//! reported so the user can judge whether the grid resolves the body.
use super::PolyMeshSpec;
use crate::{require, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MAX_TRIANGLES: usize = 200_000;

/// A closed triangle surface in metres, triangles wound counter-clockwise seen from outside.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SurfaceBody {
    pub name: String,
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[usize; 3]>,
}

impl SurfaceBody {
    /// Watertightness, orientation and finiteness checks. Returns the enclosed volume (m^3).
    pub fn validate(&self) -> Result<f64> {
        require(
            !self.triangles.is_empty()
                && self.triangles.len() <= MAX_TRIANGLES
                && self.vertices.len() <= 3 * MAX_TRIANGLES,
            format!(
                "surface body {} needs 1..{MAX_TRIANGLES} triangles",
                self.name
            ),
        )?;
        require(
            self.vertices
                .iter()
                .flatten()
                .all(|v| v.is_finite() && v.abs() < 1e6),
            "surface vertices must be finite metres",
        )?;
        let mut edges: BTreeMap<(usize, usize), i32> = BTreeMap::new();
        let mut volume = 0.0;
        for t in &self.triangles {
            require(
                t.iter().all(|i| *i < self.vertices.len())
                    && t[0] != t[1]
                    && t[1] != t[2]
                    && t[0] != t[2],
                "surface triangles need three distinct valid vertex indices",
            )?;
            let [a, b, c] = t.map(|i| self.vertices[i]);
            volume += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
            for (u, v) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                *edges.entry((u, v)).or_default() += 1;
            }
        }
        // Closed and consistently oriented: every directed edge once, with its reverse once.
        require(
            edges
                .iter()
                .all(|((u, v), n)| *n == 1 && edges.get(&(*v, *u)) == Some(&1)),
            format!(
                "surface body {} must be closed, manifold and consistently oriented",
                self.name
            ),
        )?;
        volume /= 6.0;
        require(
            volume > 0.0,
            format!(
                "surface body {} encloses no volume or is wound inward",
                self.name
            ),
        )?;
        Ok(volume)
    }

    /// Solid flags per background cell (x fastest), from cell-centre winding numbers.
    pub fn voxelize(&self, spec: PolyMeshSpec) -> Result<Vec<bool>> {
        self.validate()?;
        require(spec.xy_shear == 0.0, "surface bodies require zero xy_shear")?;
        let [nx, ny, nz] = spec.dimensions;
        let (dx, dy, dz) = (
            spec.lengths.x / nx as f64,
            spec.lengths.y / ny as f64,
            spec.lengths.z / nz as f64,
        );
        let mut inside = vec![false; nx * ny * nz];
        // Each projected triangle contributes +-1 at its x crossing to every (y, z) column centre
        // it covers. Edges use a top-left rule so a centre on a shared edge counts once.
        let mut columns: BTreeMap<(usize, usize), Vec<(f64, i32)>> = BTreeMap::new();
        for t in &self.triangles {
            let [a, b, c] = t.map(|i| self.vertices[i]);
            let area = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]);
            if area == 0.0 {
                continue;
            }
            let (ymin, ymax) = (a[1].min(b[1]).min(c[1]), a[1].max(b[1]).max(c[1]));
            let (zmin, zmax) = (a[2].min(b[2]).min(c[2]), a[2].max(b[2]).max(c[2]));
            let j0 = ((ymin / dy - 0.5).ceil().max(0.0)) as usize;
            let k0 = ((zmin / dz - 0.5).ceil().max(0.0)) as usize;
            let j1 = ((ymax / dy - 0.5).floor()).min(ny as f64 - 1.0);
            let k1 = ((zmax / dz - 0.5).floor()).min(nz as f64 - 1.0);
            if j1 < 0.0 || k1 < 0.0 {
                continue;
            }
            let sign = if area > 0.0 { 1.0 } else { -1.0 };
            for k in k0..=(k1 as usize) {
                for j in j0..=(j1 as usize) {
                    let (py, pz) = ((j as f64 + 0.5) * dy, (k as f64 + 0.5) * dz);
                    let edge = |p: [f64; 3], q: [f64; 3]| {
                        sign * ((q[1] - p[1]) * (pz - p[2]) - (q[2] - p[2]) * (py - p[1]))
                    };
                    // Top-left rule in the oriented projection: an edge owns its boundary points
                    // when it points "up" (positive z) or is horizontal pointing "left" (negative y).
                    let owns = |p: [f64; 3], q: [f64; 3]| {
                        let (ey, ez) = (sign * (q[1] - p[1]), sign * (q[2] - p[2]));
                        ez > 0.0 || (ez == 0.0 && ey < 0.0)
                    };
                    let covered = [(a, b), (b, c), (c, a)].iter().all(|(p, q)| {
                        let w = edge(*p, *q);
                        w > 0.0 || (w == 0.0 && owns(*p, *q))
                    });
                    if !covered {
                        continue;
                    }
                    // Barycentric x of the crossing.
                    let w0 = (b[1] - py) * (c[2] - pz) - (b[2] - pz) * (c[1] - py);
                    let w1 = (c[1] - py) * (a[2] - pz) - (c[2] - pz) * (a[1] - py);
                    let x = (w0 * a[0] + w1 * b[0] + (area - w0 - w1) * c[0]) / area;
                    // Outward normal x component positive means the ray leaves the body here.
                    columns
                        .entry((j, k))
                        .or_default()
                        .push((x, if area > 0.0 { -1 } else { 1 }));
                }
            }
        }
        for ((j, k), mut crossings) in columns {
            crossings.sort_by(|p, q| p.0.total_cmp(&q.0));
            let mut winding = 0;
            let mut next = crossings.iter().peekable();
            for i in 0..nx {
                let x = (i as f64 + 0.5) * dx;
                while let Some((cx, w)) = next.peek() {
                    if *cx <= x {
                        winding += w;
                        next.next();
                    } else {
                        break;
                    }
                }
                if winding != 0 {
                    if let Some(cell) = inside.get_mut((k * ny + j) * nx + i) {
                        *cell = true;
                    }
                }
            }
        }
        Ok(inside)
    }
}
