//! Axis-aligned hexahedral octree AMR: marked refinement and sibling coarsening with 2:1
//! balancing, exact overlap maps, conservative remapping, face-jump indicators, Dorfler marking
//! and a combined adapt step, ported from upstream `cfd/fvm/adaptive_mesh` (MIT, cfd_solvers
//! contributors).
use super::{BoundaryPatch, Face, PolyCell, PolyMesh};
use crate::{core::Vec3, require, Result, SolverError};
use std::collections::BTreeMap;

const EPS: f64 = 1e-12;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AdaptiveHexCell {
    pub lower: Vec3,
    pub upper: Vec3,
    pub level: usize,
    pub root: usize,
    /// Three bits per refinement level identifying the child octant path.
    pub lineage: u64,
}

fn finite(p: Vec3) -> bool {
    p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
}
fn cell_volume(c: &AdaptiveHexCell) -> f64 {
    (c.upper.x - c.lower.x) * (c.upper.y - c.lower.y) * (c.upper.z - c.lower.z)
}
fn cell_center(c: &AdaptiveHexCell) -> Vec3 {
    Vec3::new(
        (c.lower.x + c.upper.x) * 0.5,
        (c.lower.y + c.upper.y) * 0.5,
        (c.lower.z + c.upper.z) * 0.5,
    )
}
fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= EPS * 1f64.max(a.abs()).max(b.abs())
}
fn overlap_1d(a0: f64, a1: f64, b0: f64, b1: f64) -> f64 {
    0f64.max(a1.min(b1) - a0.max(b0))
}
fn overlap_volume(a: &AdaptiveHexCell, b: &AdaptiveHexCell) -> f64 {
    overlap_1d(a.lower.x, a.upper.x, b.lower.x, b.upper.x)
        * overlap_1d(a.lower.y, a.upper.y, b.lower.y, b.upper.y)
        * overlap_1d(a.lower.z, a.upper.z, b.lower.z, b.upper.z)
}

struct Contact {
    axis: usize,
    plane: f64,
    lo0: f64,
    hi0: f64,
    lo1: f64,
    hi1: f64,
    direction: f64,
}

fn comp(v: Vec3, i: usize) -> f64 {
    [v.x, v.y, v.z][i]
}

fn face_contact(a: &AdaptiveHexCell, b: &AdaptiveHexCell) -> Option<Contact> {
    for (axis, (t0, t1)) in [(0, (1, 2)), (1, (0, 2)), (2, (0, 1))] {
        let (au, al, bu, bl) = (
            comp(a.upper, axis),
            comp(a.lower, axis),
            comp(b.upper, axis),
            comp(b.lower, axis),
        );
        if !(near(au, bl) || near(bu, al)) {
            continue;
        }
        let o0 = overlap_1d(
            comp(a.lower, t0),
            comp(a.upper, t0),
            comp(b.lower, t0),
            comp(b.upper, t0),
        );
        let o1 = overlap_1d(
            comp(a.lower, t1),
            comp(a.upper, t1),
            comp(b.lower, t1),
            comp(b.upper, t1),
        );
        if o0 > EPS && o1 > EPS {
            let direction = if near(au, bl) { 1.0 } else { -1.0 };
            return Some(Contact {
                axis,
                plane: if direction > 0.0 { au } else { al },
                lo0: comp(a.lower, t0).max(comp(b.lower, t0)),
                hi0: comp(a.upper, t0).min(comp(b.upper, t0)),
                lo1: comp(a.lower, t1).max(comp(b.lower, t1)),
                hi1: comp(a.upper, t1).min(comp(b.upper, t1)),
                direction,
            });
        }
    }
    None
}

fn split_cell(c: &AdaptiveHexCell) -> [AdaptiveHexCell; 8] {
    let m = cell_center(c);
    std::array::from_fn(|child| {
        let (hx, hy, hz) = (child & 1 != 0, child & 2 != 0, child & 4 != 0);
        AdaptiveHexCell {
            lower: Vec3::new(
                if hx { m.x } else { c.lower.x },
                if hy { m.y } else { c.lower.y },
                if hz { m.z } else { c.lower.z },
            ),
            upper: Vec3::new(
                if hx { c.upper.x } else { m.x },
                if hy { c.upper.y } else { m.y },
                if hz { c.upper.z } else { m.z },
            ),
            level: c.level + 1,
            root: c.root,
            lineage: (c.lineage << 3) | child as u64,
        }
    })
}

fn merge_siblings(s: &[AdaptiveHexCell]) -> AdaptiveHexCell {
    let mut p = s[0];
    for c in s {
        p.lower = Vec3::new(
            p.lower.x.min(c.lower.x),
            p.lower.y.min(c.lower.y),
            p.lower.z.min(c.lower.z),
        );
        p.upper = Vec3::new(
            p.upper.x.max(c.upper.x),
            p.upper.y.max(c.upper.y),
            p.upper.z.max(c.upper.z),
        );
    }
    p.level = s[0].level - 1;
    p.lineage = s[0].lineage >> 3;
    p
}

/// Non-overlapping axis-aligned cells exactly tiling a box domain.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveHexMesh {
    cells: Vec<AdaptiveHexCell>,
    domain_lower: Vec3,
    domain_upper: Vec3,
}

impl AdaptiveHexMesh {
    pub fn new(
        cells: Vec<AdaptiveHexCell>,
        domain_lower: Vec3,
        domain_upper: Vec3,
    ) -> Result<Self> {
        let m = Self {
            cells,
            domain_lower,
            domain_upper,
        };
        m.validate()?;
        Ok(m)
    }
    /// Upstream default lengths 1.
    pub fn cartesian(nx: usize, ny: usize, nz: usize, lx: f64, ly: f64, lz: f64) -> Result<Self> {
        require(
            nx > 0 && ny > 0 && nz > 0 && lx > 0.0 && ly > 0.0 && lz > 0.0,
            "adaptive Cartesian mesh requires positive dimensions and lengths",
        )?;
        require(
            nx.checked_mul(ny)
                .and_then(|v| v.checked_mul(nz))
                .is_some_and(|n| n <= 1 << 22),
            "adaptive Cartesian mesh too large",
        )?;
        let (dx, dy, dz) = (lx / nx as f64, ly / ny as f64, lz / nz as f64);
        let mut cells = Vec::with_capacity(nx * ny * nz);
        let mut root = 0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let (x0, y0, z0) = (i as f64 * dx, j as f64 * dy, k as f64 * dz);
                    cells.push(AdaptiveHexCell {
                        lower: Vec3::new(x0, y0, z0),
                        upper: Vec3::new(x0 + dx, y0 + dy, z0 + dz),
                        level: 0,
                        root,
                        lineage: 0,
                    });
                    root += 1;
                }
            }
        }
        Self::new(cells, Vec3::new(0.0, 0.0, 0.0), Vec3::new(lx, ly, lz))
    }
    fn validate(&self) -> Result<()> {
        require(
            !self.cells.is_empty(),
            "AdaptiveHexMesh requires at least one cell",
        )?;
        let (lo, hi) = (self.domain_lower, self.domain_upper);
        require(
            finite(lo) && finite(hi) && hi.x > lo.x && hi.y > lo.y && hi.z > lo.z,
            "invalid adaptive mesh domain",
        )?;
        let mut volume = 0.0;
        for c in &self.cells {
            require(
                finite(c.lower)
                    && finite(c.upper)
                    && c.upper.x > c.lower.x
                    && c.upper.y > c.lower.y
                    && c.upper.z > c.lower.z,
                "invalid adaptive cell bounds",
            )?;
            require(
                !(c.lower.x < lo.x - EPS
                    || c.lower.y < lo.y - EPS
                    || c.lower.z < lo.z - EPS
                    || c.upper.x > hi.x + EPS
                    || c.upper.y > hi.y + EPS
                    || c.upper.z > hi.z + EPS),
                "adaptive cell lies outside domain",
            )?;
            require(
                c.level <= 20,
                "adaptive refinement level exceeds lineage capacity",
            )?;
            volume += cell_volume(c);
        }
        let domain = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
        require(
            (volume - domain).abs() <= 1e-9 * domain.max(1.0),
            "adaptive cells do not conserve domain volume",
        )?;
        for i in 0..self.cells.len() {
            for j in i + 1..self.cells.len() {
                require(
                    overlap_volume(&self.cells[i], &self.cells[j]) <= 1e-13 * domain,
                    "adaptive cells overlap",
                )?;
            }
        }
        Ok(())
    }
    pub fn cells(&self) -> &[AdaptiveHexCell] {
        &self.cells
    }
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
    pub fn domain_lower(&self) -> Vec3 {
        self.domain_lower
    }
    pub fn domain_upper(&self) -> Vec3 {
        self.domain_upper
    }
    pub fn total_volume(&self) -> f64 {
        self.cells.iter().map(cell_volume).fold(0.0, |a, b| a + b)
    }
    /// Polyhedral mesh with one face per face contact (hanging interfaces split exactly) and
    /// patches left/right/bottom/top/front/back.
    pub fn poly_mesh(&self) -> Result<PolyMesh> {
        let cells: Vec<PolyCell> = self
            .cells
            .iter()
            .map(|c| PolyCell {
                center: cell_center(c),
                volume: cell_volume(c),
            })
            .collect();
        let patches = ["left", "right", "bottom", "top", "front", "back"]
            .iter()
            .map(|n| BoundaryPatch { name: (*n).into() })
            .collect();
        let mut faces = Vec::new();
        for i in 0..self.cells.len() {
            for j in i + 1..self.cells.len() {
                let Some(c) = face_contact(&self.cells[i], &self.cells[j]) else {
                    continue;
                };
                let a = (c.hi0 - c.lo0) * (c.hi1 - c.lo1);
                let (m0, m1) = ((c.lo0 + c.hi0) * 0.5, (c.lo1 + c.hi1) * 0.5);
                let (center, area) = match c.axis {
                    0 => (
                        Vec3::new(c.plane, m0, m1),
                        Vec3::new(c.direction * a, 0.0, 0.0),
                    ),
                    1 => (
                        Vec3::new(m0, c.plane, m1),
                        Vec3::new(0.0, c.direction * a, 0.0),
                    ),
                    _ => (
                        Vec3::new(m0, m1, c.plane),
                        Vec3::new(0.0, 0.0, c.direction * a),
                    ),
                };
                faces.push(Face::internal(i, j, center, area));
            }
        }
        let (lo, hi) = (self.domain_lower, self.domain_upper);
        for (i, c) in self.cells.iter().enumerate() {
            let (dx, dy, dz) = (
                c.upper.x - c.lower.x,
                c.upper.y - c.lower.y,
                c.upper.z - c.lower.z,
            );
            let m = cell_center(c);
            if near(c.lower.x, lo.x) {
                faces.push(Face::boundary_face(
                    i,
                    Vec3::new(c.lower.x, m.y, m.z),
                    Vec3::new(-dy * dz, 0.0, 0.0),
                    0,
                ));
            }
            if near(c.upper.x, hi.x) {
                faces.push(Face::boundary_face(
                    i,
                    Vec3::new(c.upper.x, m.y, m.z),
                    Vec3::new(dy * dz, 0.0, 0.0),
                    1,
                ));
            }
            if near(c.lower.y, lo.y) {
                faces.push(Face::boundary_face(
                    i,
                    Vec3::new(m.x, c.lower.y, m.z),
                    Vec3::new(0.0, -dx * dz, 0.0),
                    2,
                ));
            }
            if near(c.upper.y, hi.y) {
                faces.push(Face::boundary_face(
                    i,
                    Vec3::new(m.x, c.upper.y, m.z),
                    Vec3::new(0.0, dx * dz, 0.0),
                    3,
                ));
            }
            if near(c.lower.z, lo.z) {
                faces.push(Face::boundary_face(
                    i,
                    Vec3::new(m.x, m.y, c.lower.z),
                    Vec3::new(0.0, 0.0, -dx * dy),
                    4,
                ));
            }
            if near(c.upper.z, hi.z) {
                faces.push(Face::boundary_face(
                    i,
                    Vec3::new(m.x, m.y, c.upper.z),
                    Vec3::new(0.0, 0.0, dx * dy),
                    5,
                ));
            }
        }
        PolyMesh::new(cells, faces, patches)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CellOverlapWeight {
    pub source_cell: usize,
    pub volume: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TopologyChangeMap {
    pub old_cell_count: usize,
    pub new_cell_count: usize,
    pub new_from_old: Vec<Vec<CellOverlapWeight>>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveTopologyResult {
    pub mesh: AdaptiveHexMesh,
    pub mapping: TopologyChangeMap,
}

pub trait MeshTopologyOperation {
    fn apply(&self, mesh: &AdaptiveHexMesh) -> Result<AdaptiveTopologyResult>;
}

fn refine_once(
    mesh: &AdaptiveHexMesh,
    marked: &[bool],
    max_level: usize,
) -> Result<AdaptiveHexMesh> {
    require(
        marked.len() == mesh.cell_count(),
        "adaptive refinement mark size mismatch",
    )?;
    let mut out = Vec::with_capacity(mesh.cell_count() * 2);
    for (c, m) in mesh.cells.iter().zip(marked) {
        if *m && c.level < max_level {
            out.extend(split_cell(c));
        } else {
            out.push(*c);
        }
    }
    AdaptiveHexMesh::new(out, mesh.domain_lower, mesh.domain_upper)
}

fn enforce_balance(mut mesh: AdaptiveHexMesh, max_level: usize) -> Result<AdaptiveHexMesh> {
    for _ in 0..64 {
        let mut mark = vec![false; mesh.cell_count()];
        let mut any = false;
        for i in 0..mesh.cell_count() {
            for j in i + 1..mesh.cell_count() {
                if face_contact(&mesh.cells[i], &mesh.cells[j]).is_none() {
                    continue;
                }
                let (li, lj) = (mesh.cells[i].level, mesh.cells[j].level);
                if li + 1 < lj {
                    require(
                        li < max_level,
                        "2:1 balance requires refinement above maximum level",
                    )?;
                    mark[i] = true;
                    any = true;
                } else if lj + 1 < li {
                    require(
                        lj < max_level,
                        "2:1 balance requires refinement above maximum level",
                    )?;
                    mark[j] = true;
                    any = true;
                }
            }
        }
        if !any {
            return Ok(mesh);
        }
        mesh = refine_once(&mesh, &mark, max_level)?;
    }
    Err(SolverError::NumericalFailure(
        "2:1 balancing did not converge".into(),
    ))
}

fn coarsen_once(
    mesh: &AdaptiveHexMesh,
    marked: &[bool],
    min_level: usize,
) -> Result<AdaptiveHexMesh> {
    require(
        marked.len() == mesh.cell_count(),
        "adaptive coarsening mark size mismatch",
    )?;
    let mut groups: BTreeMap<(usize, usize, u64), Vec<usize>> = BTreeMap::new();
    for (i, c) in mesh.cells.iter().enumerate() {
        if c.level > min_level {
            groups
                .entry((c.root, c.level, c.lineage >> 3))
                .or_default()
                .push(i);
        }
    }
    let mut consumed = vec![false; mesh.cell_count()];
    let mut parents = Vec::new();
    for indices in groups.values() {
        if indices.len() != 8 {
            continue;
        }
        let mut seen = [false; 8];
        let mut all = true;
        for &i in indices {
            all &= marked[i];
            seen[(mesh.cells[i].lineage & 7) as usize] = true;
        }
        if !all || !seen.iter().all(|s| *s) {
            continue;
        }
        let siblings: Vec<AdaptiveHexCell> = indices.iter().map(|i| mesh.cells[*i]).collect();
        indices.iter().for_each(|i| consumed[*i] = true);
        parents.push(merge_siblings(&siblings));
    }
    let mut out: Vec<AdaptiveHexCell> = mesh
        .cells
        .iter()
        .zip(&consumed)
        .filter(|(_, c)| !**c)
        .map(|(c, _)| *c)
        .collect();
    out.extend(parents);
    AdaptiveHexMesh::new(out, mesh.domain_lower, mesh.domain_upper)
}

/// Exact overlap volumes of every new cell with the old cells.
pub fn build_topology_change_map(
    old_mesh: &AdaptiveHexMesh,
    new_mesh: &AdaptiveHexMesh,
) -> Result<TopologyChangeMap> {
    let mut map = TopologyChangeMap {
        old_cell_count: old_mesh.cell_count(),
        new_cell_count: new_mesh.cell_count(),
        new_from_old: vec![Vec::new(); new_mesh.cell_count()],
    };
    for (j, nc) in new_mesh.cells.iter().enumerate() {
        let mut covered = 0.0;
        for (i, oc) in old_mesh.cells.iter().enumerate() {
            let v = overlap_volume(oc, nc);
            if v > EPS {
                map.new_from_old[j].push(CellOverlapWeight {
                    source_cell: i,
                    volume: v,
                });
                covered += v;
            }
        }
        let target = cell_volume(nc);
        require(
            (covered - target).abs() <= 1e-9 * target.max(1.0),
            "topology change map does not fully cover target cell",
        )?;
    }
    Ok(map)
}

/// Volume-weighted averaging over exact overlaps (conserves the integral).
pub fn conservative_adaptive_remap(
    old_mesh: &AdaptiveHexMesh,
    old_field: &[f64],
    new_mesh: &AdaptiveHexMesh,
) -> Result<Vec<f64>> {
    require(
        old_field.len() == old_mesh.cell_count(),
        "adaptive remap field size mismatch",
    )?;
    require(
        old_field.iter().all(|v| v.is_finite()),
        "adaptive remap requires finite field values",
    )?;
    let map = build_topology_change_map(old_mesh, new_mesh)?;
    Ok(new_mesh
        .cells
        .iter()
        .zip(&map.new_from_old)
        .map(|(c, ws)| {
            let volume = cell_volume(c);
            ws.iter().fold(0.0, |acc, w| {
                acc + old_field[w.source_cell] * w.volume / volume
            })
        })
        .collect())
}

/// Refine marked cells below `maximum_level`, then (by default) restore 2:1 face balance.
#[derive(Clone, Debug)]
pub struct MarkedHexRefinement {
    marked: Vec<bool>,
    maximum_level: usize,
    enforce_two_to_one: bool,
}

impl MarkedHexRefinement {
    pub fn new(marked: Vec<bool>, maximum_level: usize, enforce_two_to_one: bool) -> Result<Self> {
        require(
            maximum_level <= 20,
            "maximum adaptive level exceeds lineage capacity",
        )?;
        Ok(Self {
            marked,
            maximum_level,
            enforce_two_to_one,
        })
    }
}

impl MeshTopologyOperation for MarkedHexRefinement {
    fn apply(&self, mesh: &AdaptiveHexMesh) -> Result<AdaptiveTopologyResult> {
        let mut out = refine_once(mesh, &self.marked, self.maximum_level)?;
        if self.enforce_two_to_one {
            out = enforce_balance(out, self.maximum_level)?;
        }
        let mapping = build_topology_change_map(mesh, &out)?;
        Ok(AdaptiveTopologyResult { mesh: out, mapping })
    }
}

/// Merge complete octets of marked siblings above `minimum_level` (parents appended last).
#[derive(Clone, Debug)]
pub struct MarkedHexCoarsening {
    marked: Vec<bool>,
    minimum_level: usize,
    enforce_two_to_one: bool,
}

impl MarkedHexCoarsening {
    /// Upstream defaults `minimum_level = 0`, `enforce_two_to_one = true`.
    pub fn new(marked: Vec<bool>, minimum_level: usize, enforce_two_to_one: bool) -> Self {
        Self {
            marked,
            minimum_level,
            enforce_two_to_one,
        }
    }
}

impl MeshTopologyOperation for MarkedHexCoarsening {
    fn apply(&self, mesh: &AdaptiveHexMesh) -> Result<AdaptiveTopologyResult> {
        let mut out = coarsen_once(mesh, &self.marked, self.minimum_level)?;
        if self.enforce_two_to_one {
            let max_level = mesh.cells.iter().map(|c| c.level).max().unwrap_or(0);
            out = enforce_balance(out, max_level)?;
        }
        let mapping = build_topology_change_map(mesh, &out)?;
        Ok(AdaptiveTopologyResult { mesh: out, mapping })
    }
}

pub fn satisfies_two_to_one(mesh: &AdaptiveHexMesh) -> bool {
    for i in 0..mesh.cell_count() {
        for j in i + 1..mesh.cell_count() {
            if face_contact(&mesh.cells[i], &mesh.cells[j]).is_none() {
                continue;
            }
            let (a, b) = (mesh.cells[i].level, mesh.cells[j].level);
            if a > b + 1 || b > a + 1 {
                return false;
            }
        }
    }
    true
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// `sqrt(sum_faces |S| (jump)^2 / d)` per cell over internal faces.
pub fn scalar_jump_error_indicator(mesh: &AdaptiveHexMesh, field: &[f64]) -> Result<Vec<f64>> {
    require(
        field.len() == mesh.cell_count(),
        "AMR indicator field size mismatch",
    )?;
    require(
        field.iter().all(|v| v.is_finite()),
        "AMR indicator requires finite field values",
    )?;
    let mut eta = vec![0.0; mesh.cell_count()];
    let poly = mesh.poly_mesh()?;
    for f in poly.faces() {
        let Some(nb) = f.neighbour else { continue };
        let (a, b) = (poly.cells()[f.owner].center, poly.cells()[nb].center);
        let d = Vec3::new(b.x - a.x, b.y - a.y, b.z - a.z);
        let distance = dot(d, d).sqrt();
        if distance.is_nan() || distance <= 0.0 {
            continue;
        }
        let jump = field[nb] - field[f.owner];
        let contribution = dot(f.area, f.area).sqrt() * jump * jump / distance;
        eta[f.owner] += contribution;
        eta[nb] += contribution;
    }
    Ok(eta.iter().map(|v| v.max(0.0).sqrt()).collect())
}

/// Stable descending Dorfler bulk marking; an all-zero indicator marks nothing.
pub fn mark_dorfler(indicators: &[f64], theta: f64) -> Result<Vec<bool>> {
    require(
        theta > 0.0 && theta <= 1.0,
        "Dorfler theta must lie in (0,1]",
    )?;
    require(
        indicators.iter().all(|v| *v >= 0.0 && v.is_finite()),
        "Dorfler indicators must be finite non-negative",
    )?;
    let total: f64 = indicators.iter().map(|v| v * v).fold(0.0, |a, b| a + b);
    let mut mark = vec![false; indicators.len()];
    if total <= 0.0 {
        return Ok(mark);
    }
    let mut order: Vec<usize> = (0..indicators.len()).collect();
    order.sort_by(|a, b| indicators[*b].total_cmp(&indicators[*a]));
    let mut selected = 0.0;
    for i in order {
        mark[i] = true;
        selected += indicators[i] * indicators[i];
        if selected >= theta * total {
            break;
        }
    }
    Ok(mark)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveAmrConfig {
    pub minimum_level: usize,
    pub maximum_level: usize,
    pub refine_theta: f64,
    pub coarsen_relative_threshold: f64,
    pub enforce_two_to_one: bool,
}

impl Default for AdaptiveAmrConfig {
    fn default() -> Self {
        Self {
            minimum_level: 0,
            maximum_level: 4,
            refine_theta: 0.5,
            coarsen_relative_threshold: 0.05,
            enforce_two_to_one: true,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveAmrResult {
    pub mesh: AdaptiveHexMesh,
    pub field: Vec<f64>,
    pub indicator: Vec<f64>,
    pub refined_cell_count: usize,
    pub coarsened_cell_count: usize,
}

/// Coarsen low-error octets, then Dorfler-refine, remapping conservatively after each step.
pub fn adapt_scalar_field(
    mesh: &AdaptiveHexMesh,
    field: &[f64],
    cfg: &AdaptiveAmrConfig,
) -> Result<AdaptiveAmrResult> {
    require(field.len() == mesh.cell_count(), "AMR field size mismatch")?;
    require(
        cfg.minimum_level <= cfg.maximum_level && cfg.maximum_level <= 20,
        "invalid AMR level bounds",
    )?;
    require(
        cfg.refine_theta > 0.0
            && cfg.refine_theta <= 1.0
            && (0.0..=1.0).contains(&cfg.coarsen_relative_threshold),
        "invalid AMR marking controls",
    )?;
    let mut current = mesh.clone();
    let mut values = field.to_vec();
    let indicator = scalar_jump_error_indicator(&current, &values)?;
    let before_coarsen = current.cell_count();
    let peak = indicator.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let threshold = cfg.coarsen_relative_threshold * peak;
    let coarsen_mark: Vec<bool> = current
        .cells
        .iter()
        .zip(&indicator)
        .map(|(c, v)| c.level > cfg.minimum_level && *v <= threshold)
        .collect();
    let coarsened =
        MarkedHexCoarsening::new(coarsen_mark, cfg.minimum_level, cfg.enforce_two_to_one)
            .apply(&current)?;
    if coarsened.mesh.cell_count() != current.cell_count() {
        values = conservative_adaptive_remap(&current, &values, &coarsened.mesh)?;
    }
    current = coarsened.mesh;
    let after_coarsen = current.cell_count();
    let indicator = scalar_jump_error_indicator(&current, &values)?;
    let mut refine_mark = mark_dorfler(&indicator, cfg.refine_theta)?;
    for (m, c) in refine_mark.iter_mut().zip(&current.cells) {
        if c.level >= cfg.maximum_level {
            *m = false;
        }
    }
    let before_refine = current.cell_count();
    let refined = MarkedHexRefinement::new(refine_mark, cfg.maximum_level, cfg.enforce_two_to_one)?
        .apply(&current)?;
    if refined.mesh.cell_count() != current.cell_count() {
        values = conservative_adaptive_remap(&current, &values, &refined.mesh)?;
    }
    current = refined.mesh;
    let indicator = scalar_jump_error_indicator(&current, &values)?;
    let refined_cell_count = current.cell_count().saturating_sub(before_refine);
    Ok(AdaptiveAmrResult {
        mesh: current,
        field: values,
        indicator,
        refined_cell_count,
        coarsened_cell_count: before_coarsen.saturating_sub(after_coarsen),
    })
}
