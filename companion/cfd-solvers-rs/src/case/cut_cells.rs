//! Cut cells for surface bodies (`cut_cells_v1`).
//!
//! Without cut cells a body is a staircase of whole cells (a cell is solid when its centre is
//! inside). With them, a cell the surface crosses keeps its fluid: each connected piece of fluid
//! in the cell (a *fragment*; a wall thinner than the cell splits it in two) becomes a solver cell
//! with its own volume and centroid, the grid faces keep the open area between fragments that
//! touch through them, and each fragment gets a wall face where the surface cuts it.
//!
//! The geometry is measured on a finer lattice: every cell is split `sub` times along each axis
//! and each sub-cell is classified by the same exact winding rule the staircase uses
//! (`SurfaceBody::voxelize`). Fragments are the 6-connected fluid sub-cells of a cell. Volumes and
//! open areas therefore converge like `1/sub` (one sub-cell layer); they are not exact polygon
//! clipping. A fragment's wall face is closed exactly: its area vector is minus the sum of the
//! fragment's open face area vectors, so every control volume is closed whatever the sampling.
//!
//! Fragments with less fluid than [`MERGE_BELOW`] of their cell's box are merged into the
//! neighbour they share the largest open face with, so no tiny cell makes the time step or the
//! pressure solve stiff. Two fragments of one cell never share a face, so a thin wall is never
//! merged through.
use super::bodies::SurfaceBody;
use super::PolyMeshSpec;
use crate::{core::Vec3, fvm::Face, require, Result};

/// Sub-cells per cell along each axis at most.
pub const MAX_SUB: usize = 8;
/// Most sub-cells of the measuring lattice in all.
const MAX_SUB_CELLS: usize = 16_000_000;
/// A fragment with less fluid than this fraction of its cell's box is merged into a neighbour.
pub const MERGE_BELOW: f64 = 0.3;
/// A sub-cell in no fragment (in a body or an obstacle cell).
const SOLID: u32 = u32::MAX;

/// One connected piece of a cell's fluid.
pub(super) struct Fragment {
    pub cell: usize,
    /// Its share of the cell's box.
    pub fraction: f64,
    pub centroid: Vec3,
}

/// Open area through one grid face from one fragment to the fragment across it (or to the
/// exterior or an obstacle when `across` is `None`).
pub(super) struct Link {
    pub from: usize,
    pub across: Option<usize>,
    /// The share of the grid face that is open.
    pub fraction: f64,
    pub centre: Vec3,
}

/// The cut geometry of a grid.
pub(super) struct Cut {
    sub: usize,
    dims: [usize; 3],
    nodes: [Vec<f64>; 3],
    /// Per sub-cell: 0 for fluid, else 1 + the body that holds it.
    owner: Vec<u8>,
    /// Per cell: grid obstacle (always whole cells, as without cut cells).
    blocked: Vec<bool>,
    /// Per sub-cell: its fragment, or [`SOLID`].
    fragment: Vec<u32>,
    fragments: Vec<Fragment>,
}

fn sub_nodes(coarse: &[f64], sub: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity((coarse.len() - 1) * sub + 1);
    for w in coarse.windows(2) {
        for s in 0..sub {
            out.push(w[0] + (w[1] - w[0]) * s as f64 / sub as f64);
        }
    }
    out.push(coarse[coarse.len() - 1]);
    out
}

fn mid(nodes: &[f64], i: usize) -> f64 {
    0.5 * (nodes[i] + nodes[i + 1])
}

fn set(p: &mut Vec3, axis: usize, value: f64) {
    match axis {
        0 => p.x = value,
        1 => p.y = value,
        _ => p.z = value,
    }
}

impl Cut {
    /// Measure the bodies on a lattice `sub` times finer than `spec` (at most [`MAX_SUB`], fewer
    /// when the lattice would be too large). `blocked` marks the grid-obstacle cells.
    pub(super) fn new(spec: &PolyMeshSpec, bodies: &[SurfaceBody], blocked: Vec<bool>) -> Result<Self> {
        let dims = spec.dimensions;
        let cells: usize = dims.iter().product();
        require(blocked.len() == cells, "cut cells: obstacle table size mismatch")?;
        let mut sub = MAX_SUB;
        while sub > 2 && cells.saturating_mul(sub * sub * sub) > MAX_SUB_CELLS {
            sub -= 1;
        }
        let coarse: [Vec<f64>; 3] = std::array::from_fn(|a| spec.nodes(a));
        let nodes: [Vec<f64>; 3] = std::array::from_fn(|a| sub_nodes(&coarse[a], sub));
        let fine = PolyMeshSpec {
            dimensions: std::array::from_fn(|a| dims[a] * sub),
            lengths: spec.lengths,
            xy_shear: 0.0,
            coordinates: Some(super::MeshCoordinates { x: nodes[0].clone(), y: nodes[1].clone(), z: nodes[2].clone() }),
        };
        let mut owner = vec![0u8; cells * sub * sub * sub];
        for (bi, body) in bodies.iter().enumerate() {
            require(bi < 255, "cut cells take at most 254 surface bodies")?;
            let inside = body.voxelize(&fine)?;
            for (slot, inside) in owner.iter_mut().zip(inside) {
                if inside {
                    require(*slot == 0, format!("surface body {} overlaps another body", body.name))?;
                    *slot = bi as u8 + 1;
                }
            }
        }
        let mut cut = Self { sub, dims, nodes, owner, blocked, fragment: Vec::new(), fragments: Vec::new() };
        cut.find_fragments();
        Ok(cut)
    }

    fn fine_dims(&self) -> [usize; 3] {
        std::array::from_fn(|a| self.dims[a] * self.sub)
    }

    fn fine_index(&self, at: [usize; 3]) -> usize {
        let d = self.fine_dims();
        (at[2] * d[1] + at[1]) * d[0] + at[0]
    }

    fn coarse_of(&self, at: [usize; 3]) -> usize {
        let c: [usize; 3] = std::array::from_fn(|a| at[a] / self.sub);
        (c[2] * self.dims[1] + c[1]) * self.dims[0] + c[0]
    }

    fn coarse_at(&self, cell: usize) -> [usize; 3] {
        [cell % self.dims[0], (cell / self.dims[0]) % self.dims[1], cell / (self.dims[0] * self.dims[1])]
    }

    fn sub_centre(&self, at: [usize; 3]) -> Vec3 {
        Vec3::new(mid(&self.nodes[0], at[0]), mid(&self.nodes[1], at[1]), mid(&self.nodes[2], at[2]))
    }

    /// Label each cell's 6-connected fluid sub-cells: every component is a fragment.
    fn find_fragments(&mut self) {
        let s = self.sub;
        let cells: usize = self.dims.iter().product();
        self.fragment = vec![SOLID; self.owner.len()];
        let mut stack = Vec::new();
        for cell in 0..cells {
            if self.blocked[cell] {
                continue;
            }
            let c = self.coarse_at(cell);
            let lo: [usize; 3] = std::array::from_fn(|a| c[a] * s);
            for k in lo[2]..lo[2] + s {
                for j in lo[1]..lo[1] + s {
                    for i in lo[0]..lo[0] + s {
                        let start = self.fine_index([i, j, k]);
                        if self.owner[start] != 0 || self.fragment[start] != SOLID {
                            continue;
                        }
                        let id = self.fragments.len() as u32;
                        let mut count = 0usize;
                        let mut sum = Vec3::ZERO;
                        self.fragment[start] = id;
                        stack.push([i, j, k]);
                        while let Some(at) = stack.pop() {
                            count += 1;
                            sum += self.sub_centre(at);
                            for axis in 0..3 {
                                for up in [false, true] {
                                    let mut n = at;
                                    if up {
                                        if n[axis] + 1 >= lo[axis] + s {
                                            continue;
                                        }
                                        n[axis] += 1;
                                    } else {
                                        if n[axis] == lo[axis] {
                                            continue;
                                        }
                                        n[axis] -= 1;
                                    }
                                    let ni = self.fine_index(n);
                                    if self.owner[ni] == 0 && self.fragment[ni] == SOLID {
                                        self.fragment[ni] = id;
                                        stack.push(n);
                                    }
                                }
                            }
                        }
                        self.fragments.push(Fragment {
                            cell,
                            fraction: count as f64 / (s * s * s) as f64,
                            centroid: sum * (1.0 / count as f64),
                        });
                    }
                }
            }
        }
    }

    pub(super) fn fragments(&self) -> &[Fragment] {
        &self.fragments
    }

    /// The open links through the grid face of `cell` on its `axis` side `upper`: from each of
    /// its fragments to the fragment across (when `across` names a fluid-capable cell), else to
    /// the exterior or an obstacle. A sample is open when the sub-cells on both sides are fluid
    /// (one side for the exterior and obstacles).
    pub(super) fn links(&self, cell: usize, axis: usize, upper: bool, across: bool) -> Vec<Link> {
        let c = self.coarse_at(cell);
        let s = self.sub;
        let own_layer = if upper { c[axis] * s + s - 1 } else { c[axis] * s };
        let plane = self.nodes[axis][if upper { (c[axis] + 1) * s } else { c[axis] * s }];
        let other_layer = if upper { own_layer + 1 } else { own_layer.wrapping_sub(1) };
        let [u, v] = match axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let mut found: Vec<(usize, Option<usize>, usize, Vec3)> = Vec::new();
        for b in c[v] * s..(c[v] + 1) * s {
            for a in c[u] * s..(c[u] + 1) * s {
                let mut at = [0; 3];
                at[u] = a;
                at[v] = b;
                at[axis] = own_layer;
                let mine = self.fragment[self.fine_index(at)];
                if mine == SOLID {
                    continue;
                }
                let mut p = self.sub_centre(at);
                set(&mut p, axis, plane);
                let other = if across {
                    at[axis] = other_layer;
                    let theirs = self.fragment[self.fine_index(at)];
                    if theirs == SOLID {
                        continue;
                    }
                    Some(theirs as usize)
                } else {
                    None
                };
                match found.iter_mut().find(|f| f.0 == mine as usize && f.1 == other) {
                    Some(f) => {
                        f.2 += 1;
                        f.3 += p;
                    }
                    None => found.push((mine as usize, other, 1, p)),
                }
            }
        }
        found
            .into_iter()
            .map(|(from, across, n, sum)| Link { from, across, fraction: n as f64 / (s * s) as f64, centre: sum * (1.0 / n as f64) })
            .collect()
    }

    /// For each fragment, where its fluid meets the bodies (the mean of the sub-faces between its
    /// sub-cells and body sub-cells) and the body most of them belong to.
    pub(super) fn walls(&self) -> Vec<Option<(Vec3, usize)>> {
        let d = self.fine_dims();
        let mut sum = vec![Vec3::ZERO; self.fragments.len()];
        let mut count = vec![0usize; self.fragments.len()];
        let mut votes: Vec<Vec<usize>> = vec![Vec::new(); self.fragments.len()];
        for k in 0..d[2] {
            for j in 0..d[1] {
                for i in 0..d[0] {
                    let at = [i, j, k];
                    let f = self.fragment[self.fine_index(at)];
                    if f == SOLID {
                        continue;
                    }
                    let f = f as usize;
                    for axis in 0..3 {
                        for up in [false, true] {
                            let mut n = at;
                            if up {
                                if n[axis] + 1 >= d[axis] {
                                    continue;
                                }
                                n[axis] += 1;
                            } else {
                                if n[axis] == 0 {
                                    continue;
                                }
                                n[axis] -= 1;
                            }
                            let body = self.owner[self.fine_index(n)];
                            if body == 0 || self.blocked[self.coarse_of(n)] {
                                continue;
                            }
                            count[f] += 1;
                            sum[f] += (self.sub_centre(at) + self.sub_centre(n)) * 0.5;
                            let b = usize::from(body - 1);
                            if votes[f].len() <= b {
                                votes[f].resize(b + 1, 0);
                            }
                            votes[f][b] += 1;
                        }
                    }
                }
            }
        }
        (0..self.fragments.len())
            .map(|f| {
                let body = votes[f].iter().enumerate().max_by_key(|(_, n)| **n).map(|(b, _)| b)?;
                (count[f] > 0).then(|| (sum[f] * (1.0 / count[f] as f64), body))
            })
            .collect()
    }
}

/// Merge fragments with less fluid than [`MERGE_BELOW`] into the neighbour they share the most
/// open area with: returns, per fragment, the fragment it now belongs to (itself if kept), or
/// `None` for a small fragment with no open neighbour (a sealed sliver, dropped).
pub(super) fn merge_small(fractions: &[f64], faces: &[Face]) -> Vec<Option<usize>> {
    let n = fractions.len();
    let mut target: Vec<Option<usize>> = (0..n).map(Some).collect();
    let mut order: Vec<usize> = (0..n).filter(|&c| fractions[c] > 0.0 && fractions[c] < MERGE_BELOW).collect();
    order.sort_by(|a, b| fractions[*a].total_cmp(&fractions[*b]));
    let find = |target: &[Option<usize>], mut c: usize| -> Option<usize> {
        for _ in 0..n {
            match target[c] {
                Some(t) if t == c => return Some(c),
                Some(t) => c = t,
                None => return None,
            }
        }
        None
    };
    for small in order {
        let Some(me) = find(&target, small) else { continue };
        let mut best: Option<(f64, usize)> = None;
        for f in faces {
            let Some(other) = f.neighbour else { continue };
            let partner = if f.owner == small {
                other
            } else if other == small {
                f.owner
            } else {
                continue;
            };
            let Some(root) = find(&target, partner) else { continue };
            if root == me {
                continue;
            }
            let area = f.area.norm();
            if best.is_none_or(|(e, _)| area > e) {
                best = Some((area, root));
            }
        }
        target[me] = best.map(|(_, root)| root);
    }
    (0..n).map(|c| find(&target, c)).collect()
}
