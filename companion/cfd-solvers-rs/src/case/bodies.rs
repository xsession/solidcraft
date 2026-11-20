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
    pub fn voxelize(&self, spec: &PolyMeshSpec) -> Result<Vec<bool>> {
        self.validate()?;
        require(spec.xy_shear == 0.0, "surface bodies require zero xy_shear")?;
        let [nx, ny, nz] = spec.dimensions;
        let (dx, dy, dz) = (
            spec.lengths.x / nx as f64,
            spec.lengths.y / ny as f64,
            spec.lengths.z / nz as f64,
        );
        // Cell centres along each axis: evenly spaced, or the mid-points of explicit nodes.
        let centres = |axis: usize, n: usize, d: f64| -> Vec<f64> {
            match &spec.coordinates {
                Some(_) => spec.nodes(axis).windows(2).map(|w| 0.5 * (w[0] + w[1])).collect(),
                None => (0..n).map(|i| (i as f64 + 0.5) * d).collect(),
            }
        };
        let (cx, cy, cz) = (centres(0, nx, dx), centres(1, ny, dy), centres(2, nz, dz));
        // First and last centre index within [lo, hi] (None if no centre lies there).
        let span = |c: &[f64], lo: f64, hi: f64| -> Option<(usize, usize)> {
            let first = c.partition_point(|v| *v < lo);
            let end = c.partition_point(|v| *v <= hi);
            (first < end).then(|| (first, end - 1))
        };
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
            let (Some((j0, j1)), Some((k0, k1))) = (span(&cy, ymin, ymax), span(&cz, zmin, zmax)) else {
                continue;
            };
            let sign = if area > 0.0 { 1.0 } else { -1.0 };
            for k in k0..=k1 {
                for j in j0..=j1 {
                    let (py, pz) = (cy[j], cz[k]);
                    // Evaluated with the endpoints in a fixed order, so the two triangles
                    // sharing an edge get exactly opposite values and a centre lying on the
                    // edge is counted once (in either order the rounding differs, and such a
                    // centre was counted by both triangles or by neither).
                    let edge = |p: [f64; 3], q: [f64; 3]| {
                        let swap = (q[1], q[2]) < (p[1], p[2]);
                        let (p, q) = if swap { (q, p) } else { (p, q) };
                        let w = (q[1] - p[1]) * (pz - p[2]) - (q[2] - p[2]) * (py - p[1]);
                        sign * if swap { -w } else { w }
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
            for (i, &x) in cx.iter().enumerate() {
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

/// A boundary condition on part of a surface body (`body_face_conditions_v1`): of the body's
/// solid-to-fluid grid faces, those whose nearest body triangle is one of `triangles` form the
/// boundary patch `name`, which velocity and pressure conditions then name like any other patch.
/// The other faces of the body stay on its own no-slip wall patch. On the staircase surface a
/// face is assigned by the exact distance to the triangles, so a model face keeps its grid faces
/// however the grid cuts it.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct BodyFaceCondition {
    pub name: String,
    pub body: String,
    pub triangles: Vec<usize>,
}

/// The point of triangle `abc` nearest `p` (Ericson, Real-Time Collision Detection, 5.1.5).
fn closest_on_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let sub = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
    let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let at = |o: [f64; 3], d: [f64; 3], t: f64| [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t];
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return at(a, ab, d1 / (d1 - d3));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return at(a, ac, d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return at(b, sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = 1.0 / (va + vb + vc);
    let (v, w) = (vb * denom, vc * denom);
    [a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w, a[2] + ab[2] * v + ac[2] * w]
}

impl SurfaceBody {
    /// For each point, the index of the body triangle nearest it (exact point-triangle
    /// distance; the lower index on a tie). Triangles are binned on a uniform grid over the
    /// body's box and searched in growing shells around each point, so the search stops as
    /// soon as no unvisited bin can hold a nearer triangle.
    pub fn nearest_triangles(&self, points: &[[f64; 3]]) -> Result<Vec<usize>> {
        self.validate()?;
        let mut low = [f64::INFINITY; 3];
        let mut high = [f64::NEG_INFINITY; 3];
        for v in &self.vertices {
            for a in 0..3 {
                low[a] = low[a].min(v[a]);
                high[a] = high[a].max(v[a]);
            }
        }
        let bins = ((self.triangles.len() as f64).cbrt().ceil() as usize).clamp(1, 64);
        let size: [f64; 3] = std::array::from_fn(|a| ((high[a] - low[a]) / bins as f64).max(1e-12));
        let bin_of = |x: f64, a: usize| (((x - low[a]) / size[a]).floor().max(0.0) as usize).min(bins - 1);
        let mut binned: Vec<Vec<usize>> = vec![Vec::new(); bins * bins * bins];
        for (ti, t) in self.triangles.iter().enumerate() {
            let corners = t.map(|i| self.vertices[i]);
            let lo: [usize; 3] = std::array::from_fn(|a| bin_of(corners.iter().map(|c| c[a]).fold(f64::INFINITY, f64::min), a));
            let hi: [usize; 3] = std::array::from_fn(|a| bin_of(corners.iter().map(|c| c[a]).fold(f64::NEG_INFINITY, f64::max), a));
            for k in lo[2]..=hi[2] {
                for j in lo[1]..=hi[1] {
                    for i in lo[0]..=hi[0] {
                        binned[(k * bins + j) * bins + i].push(ti);
                    }
                }
            }
        }
        let smallest = size.iter().copied().fold(f64::INFINITY, f64::min);
        let mut nearest = Vec::with_capacity(points.len());
        for p in points {
            require(p.iter().all(|v| v.is_finite()), "nearest-triangle points must be finite")?;
            let home: [usize; 3] = std::array::from_fn(|a| bin_of(p[a], a));
            let mut best = (f64::INFINITY, usize::MAX);
            for r in 0..bins {
                for k in home[2].saturating_sub(r)..=(home[2] + r).min(bins - 1) {
                    for j in home[1].saturating_sub(r)..=(home[1] + r).min(bins - 1) {
                        for i in home[0].saturating_sub(r)..=(home[0] + r).min(bins - 1) {
                            let shell = [i.abs_diff(home[0]), j.abs_diff(home[1]), k.abs_diff(home[2])].into_iter().max() == Some(r);
                            if !shell {
                                continue;
                            }
                            for &ti in &binned[(k * bins + j) * bins + i] {
                                let [a, b, c] = self.triangles[ti].map(|v| self.vertices[v]);
                                let q = closest_on_triangle(*p, a, b, c);
                                let d = ((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2)).sqrt();
                                if d < best.0 || (d == best.0 && ti < best.1) {
                                    best = (d, ti);
                                }
                            }
                        }
                    }
                }
                // Every bin beyond shell r is at least r bin widths from the point (a point
                // outside the box is clamped to its nearest bin, which only adds distance).
                if best.1 != usize::MAX && best.0 <= r as f64 * smallest {
                    break;
                }
            }
            require(best.1 != usize::MAX, "surface body has no triangles")?;
            nearest.push(best.1);
        }
        Ok(nearest)
    }
}
