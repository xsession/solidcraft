//! Automatic cell size from the bodies' own features: the narrowest gap the fluid has to pass
//! and the thinnest wall a part has, measured on their triangle surfaces along the grid's axes.
//!
//! The grid is Cartesian, so what decides how many cells span a gap or a wall is its width
//! along an axis. From each sampled face that looks along an axis (its normal within about 25°
//! of it), one ray goes out into the fluid and one into the part, both along that axis:
//!
//! - out: the first surface it meets is the far side of a **gap** (to the same part or another);
//! - in: the first surface it meets on the same part is the part's far side, its **wall**.
//!
//! The widths are robust minima: the smallest width below which at least 1 % of the sampled
//! face area lies, so a sliver face (a chamfer, a thread's end) does not set the cell size of
//! the whole study. A gap wants [`GAP_CELLS`] cells across and a wall [`WALL_CELLS`].
//!
//! Portable: triangle surfaces in, widths out.
use crate::case_scene::MeshObject;

/// Cells a gap the fluid passes should have across, at least.
pub const GAP_CELLS: f64 = 2.0;
/// Cells a part should be thick, at least.
pub const WALL_CELLS: f64 = 2.0;
/// Faces sampled at most, spread over all the bodies.
pub const MAX_SAMPLES: usize = 6_000;
/// The share of sampled face area that sets a robust minimum.
const AREA_SHARE: f64 = 0.01;
/// A face looks along an axis when its normal's component on it is at least this.
const ALONG: f64 = 0.9;

/// One measured width: how wide (mm), along which axis, and between which bodies.
#[derive(Clone, Debug, PartialEq)]
pub struct Width {
    pub mm: f64,
    pub axis: usize,
    pub from: String,
    pub to: String,
}

/// The bodies' governing gap and wall, where they have one.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Features {
    pub gap: Option<Width>,
    pub wall: Option<Width>,
}

/// A triangle's projection on the plane across `axis`, binned for axis-aligned rays.
struct Bins {
    axis: usize,
    /// The two other axes.
    across: [usize; 2],
    low: [f64; 2],
    step: [f64; 2],
    n: [usize; 2],
    cells: Vec<Vec<(usize, usize)>>,
}

impl Bins {
    fn new(objects: &[MeshObject], axis: usize, bounds: ([f64; 3], [f64; 3]), triangles: usize) -> Bins {
        let across = [(axis + 1) % 3, (axis + 2) % 3];
        let side = ((triangles as f64).sqrt() / 2.0).clamp(1.0, 256.0) as usize;
        let low = across.map(|a| bounds.0[a]);
        let step = across.map(|a| ((bounds.1[a] - bounds.0[a]) / side as f64).max(1e-9));
        let mut bins = Bins { axis, across, low, step, n: [side; 2], cells: vec![Vec::new(); side * side] };
        for (o, object) in objects.iter().enumerate() {
            for (t, tri) in object.triangles.iter().enumerate() {
                let Some(p) = corners(object, tri) else { continue };
                let (mut lo, mut hi) = ([usize::MAX; 2], [0usize; 2]);
                for corner in &p {
                    for k in 0..2 {
                        let at = bins.at(k, corner[across[k]]);
                        lo[k] = lo[k].min(at);
                        hi[k] = hi[k].max(at);
                    }
                }
                for i in lo[0]..=hi[0] {
                    for j in lo[1]..=hi[1] {
                        if let Some(cell) = bins.cells.get_mut(j * side + i) {
                            cell.push((o, t));
                        }
                    }
                }
            }
        }
        bins
    }

    fn at(&self, k: usize, v: f64) -> usize {
        (((v - self.low[k]) / self.step[k]).floor().max(0.0) as usize).min(self.n[k] - 1)
    }

    /// The first surface a ray from `p` along `self.axis` in direction `sign` meets, other
    /// than `skip`: its distance (about 0 for a surface touching `p`) and object.
    fn first(&self, objects: &[MeshObject], p: [f64; 3], sign: f64, skip: (usize, usize), eps: f64) -> Option<(f64, usize)> {
        let [b, c] = self.across;
        let cell = self.at(1, p[c]) * self.n[0] + self.at(0, p[b]);
        let mut best: Option<(f64, usize)> = None;
        for &(o, t) in self.cells.get(cell)? {
            if (o, t) == skip {
                continue;
            }
            let Some(q) = objects.get(o).and_then(|object| object.triangles.get(t).and_then(|tri| corners(object, tri))) else { continue };
            // Barycentric coordinates of p in the triangle's projection across the axis.
            let (x1, y1, x2, y2, x3, y3) = (q[0][b], q[0][c], q[1][b], q[1][c], q[2][b], q[2][c]);
            let det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
            if det.abs() < 1e-18 {
                continue;
            }
            let l1 = ((y2 - y3) * (p[b] - x3) + (x3 - x2) * (p[c] - y3)) / det;
            let l2 = ((y3 - y1) * (p[b] - x3) + (x1 - x3) * (p[c] - y3)) / det;
            let l3 = 1.0 - l1 - l2;
            let slack = -1e-9;
            if l1 < slack || l2 < slack || l3 < slack {
                continue;
            }
            let hit = l1 * q[0][self.axis] + l2 * q[1][self.axis] + l3 * q[2][self.axis];
            let t = (hit - p[self.axis]) * sign;
            // A surface the face touches (t about 0) counts too: the caller tells contact apart.
            if t > -eps && best.is_none_or(|(d, _)| t < d) {
                best = Some((t, o));
            }
        }
        best
    }
}

fn corners(object: &MeshObject, tri: &[usize; 3]) -> Option<[[f64; 3]; 3]> {
    Some([*object.vertices.get(tri[0])?, *object.vertices.get(tri[1])?, *object.vertices.get(tri[2])?])
}

/// The smallest width below which at least [`AREA_SHARE`] of `total`, the area of every
/// sampled face, lies: a feature on a sliver of the bodies' faces does not count.
fn robust(mut samples: Vec<(f64, f64, Width)>, total: f64) -> Option<Width> {
    samples.retain(|s| s.0.is_finite() && s.0 > 0.0 && s.1 > 0.0);
    samples.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut area = 0.0;
    for (_, a, width) in samples {
        area += a;
        if area >= AREA_SHARE * total {
            return Some(width);
        }
    }
    None
}

/// The governing gap and wall of `objects` (mm).
pub fn features(objects: &[MeshObject]) -> Features {
    let all = objects.iter().flat_map(|o| o.vertices.iter());
    let (mut low, mut high) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for v in all {
        for a in 0..3 {
            low[a] = low[a].min(v[a]);
            high[a] = high[a].max(v[a]);
        }
    }
    if (0..3).any(|a| !(high[a] > low[a] && (high[a] - low[a]).is_finite())) {
        return Features::default();
    }
    let size = (0..3).map(|a| high[a] - low[a]).fold(0.0, f64::max);
    let eps = size * 1e-7;
    let triangles: usize = objects.iter().map(|o| o.triangles.len()).sum();
    // Faces that look along an axis: (object, triangle, axis, sign, centroid, area).
    let mut faces = Vec::new();
    for (o, object) in objects.iter().enumerate() {
        for (t, tri) in object.triangles.iter().enumerate() {
            let Some([a, b, c]) = corners(object, tri) else { continue };
            let (u, v) = ([0, 1, 2].map(|k| b[k] - a[k]), [0, 1, 2].map(|k| c[k] - a[k]));
            let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            let length = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            if !(length.is_finite() && length > 0.0) {
                continue;
            }
            let Some(axis) = (0..3).find(|k| (n[*k] / length).abs() >= ALONG) else { continue };
            let centroid = [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0);
            faces.push((o, t, axis, n[axis].signum(), centroid, length / 2.0));
        }
    }
    if faces.is_empty() {
        return Features::default();
    }
    let stride = faces.len().div_ceil(MAX_SAMPLES).max(1);
    let bins: Vec<Bins> = (0..3).map(|axis| Bins::new(objects, axis, (low, high), triangles)).collect();
    let name = |o: usize| objects.get(o).map_or_else(String::new, |object| object.name.clone());
    let (mut gaps, mut walls, mut total) = (Vec::new(), Vec::new(), 0.0);
    for &(o, t, axis, sign, centroid, area) in faces.iter().step_by(stride) {
        let Some(bins) = bins.get(axis) else { continue };
        total += area;
        // Out into the fluid: a face resting on another part (the hit at about 0) has no gap.
        if let Some((d, other)) = bins.first(objects, centroid, sign, (o, t), eps)
            && d > eps
        {
            gaps.push((d, area, Width { mm: d, axis, from: name(o), to: name(other) }));
        }
        if let Some((d, other)) = bins.first(objects, centroid, -sign, (o, t), eps)
            && d > eps
            && other == o
        {
            walls.push((d, area, Width { mm: d, axis, from: name(o), to: name(o) }));
        }
    }
    Features { gap: robust(gaps, total), wall: robust(walls, total) }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An axis-aligned box as an outward-facing closed triangle surface.
    fn cuboid(name: &str, low: [f64; 3], high: [f64; 3]) -> MeshObject {
        let v = |i: usize| {
            [if i & 1 == 0 { low[0] } else { high[0] }, if i & 2 == 0 { low[1] } else { high[1] }, if i & 4 == 0 { low[2] } else { high[2] }]
        };
        let vertices = (0..8).map(v).collect();
        // Outward: -x, +x, -y, +y, -z, +z.
        let quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]];
        let triangles = quads.iter().flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]).collect();
        MeshObject { name: name.into(), vertices, triangles }
    }

    #[test]
    fn a_narrow_gap_between_fins_and_their_thickness_are_measured() {
        // Two 2 mm fins, 40 × 30 mm, standing 3 mm apart along x.
        let objects = [cuboid("Fin1", [0.0, 0.0, 0.0], [2.0, 40.0, 30.0]), cuboid("Fin2", [5.0, 0.0, 0.0], [7.0, 40.0, 30.0])];
        let f = features(&objects);
        let gap = f.gap.expect("a gap");
        assert!((gap.mm - 3.0).abs() < 1e-9, "{gap:?}");
        assert_eq!(gap.axis, 0);
        assert!(matches!((gap.from.as_str(), gap.to.as_str()), ("Fin1", "Fin2") | ("Fin2", "Fin1")));
        let wall = f.wall.expect("a wall");
        assert!((wall.mm - 2.0).abs() < 1e-9 && wall.axis == 0, "{wall:?}");
    }

    #[test]
    fn a_chip_on_a_board_in_a_duct() {
        // A duct with 4 mm walls (one hollow body as an outer and an inner surface), a 1.6 mm
        // board 6 mm off the floor, and a chip on it 10 mm under the lid.
        let mut duct = cuboid("Duct", [0.0, 0.0, 0.0], [100.0, 60.0, 40.0]);
        let inner = cuboid("Duct", [4.0, 4.0, 4.0], [96.0, 56.0, 36.0]);
        let offset = duct.vertices.len();
        duct.vertices.extend(inner.vertices);
        // The inner surface faces into the duct's fluid: its triangles reversed.
        duct.triangles.extend(inner.triangles.iter().map(|t| [t[0] + offset, t[2] + offset, t[1] + offset]));
        let objects = [duct, cuboid("PCB1", [10.0, 10.0, 10.0], [90.0, 50.0, 11.6]), cuboid("U1", [40.0, 20.0, 11.6], [60.0, 40.0, 26.0])];
        let f = features(&objects);
        // Narrowest gaps: under the board 6 mm, over the chip 10 mm, beside the board 6 mm.
        let gap = f.gap.expect("a gap");
        assert!((gap.mm - 6.0).abs() < 1e-9, "{gap:?}");
        // Thinnest wall: the board, 1.6 mm, as robust as its area allows.
        let wall = f.wall.expect("a wall");
        assert!((wall.mm - 1.6).abs() < 1e-9 && wall.from == "PCB1", "{wall:?}");
    }

    #[test]
    fn a_sliver_face_does_not_set_the_size_and_nothing_is_none() {
        // A 50 mm cube and a 0.1 mm speck 0.2 mm from it: a tiny share of the face area.
        let objects = [cuboid("Big", [0.0, 0.0, 0.0], [50.0, 50.0, 50.0]), cuboid("Tiny", [50.2, 10.0, 10.0], [50.3, 10.1, 10.1])];
        let f = features(&objects);
        // The tiny part's faces are a sliver of the area: the cube's own width governs.
        assert!((f.wall.unwrap().mm - 50.0).abs() < 1e-9);
        assert!(f.gap.is_none(), "{:?}", f.gap);
        assert!(features(&[]).gap.is_none());
    }
}
