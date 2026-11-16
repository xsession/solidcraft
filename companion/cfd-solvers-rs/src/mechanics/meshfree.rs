//! Deterministic uniform-cell meshfree neighbor graph (CSR) and MPI-free halo classification.
use crate::{core::Vec3, require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeshfreeNeighborGraphConfig {
    pub support_radius: f64,
    pub include_self: bool,
    pub maximum_cells: usize,
}

impl Default for MeshfreeNeighborGraphConfig {
    fn default() -> Self {
        Self {
            support_radius: 0.1,
            include_self: false,
            maximum_cells: 1 << 20,
        }
    }
}

/// Neighbors of each point are sorted by global point id.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeshfreeNeighborGraph {
    pub cell_counts: [usize; 3],
    pub cell_size: f64,
    pub origin: Vec3,
    pub offsets: Vec<usize>,
    pub indices: Vec<usize>,
}

impl MeshfreeNeighborGraph {
    pub fn point_count(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }
    pub fn neighbors(&self, point: usize) -> Result<&[usize]> {
        require(point < self.point_count(), "meshfree point index")?;
        self.indices
            .get(self.offsets[point]..self.offsets[point + 1])
            .ok_or_else(|| crate::SolverError::InvalidInput("corrupt meshfree graph".into()))
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeshfreeHaloDescriptor {
    pub owned: Vec<usize>,
    pub halo: Vec<usize>,
}

fn point_cell(p: &Vec3, origin: &Vec3, h: f64, counts: &[usize; 3]) -> [i64; 3] {
    let f = |c: f64, o: f64, n: usize| (((c - o) / h).floor() as i64).clamp(0, n as i64 - 1);
    [
        f(p.x, origin.x, counts[0]),
        f(p.y, origin.y, counts[1]),
        f(p.z, origin.z, counts[2]),
    ]
}

fn flatten(c: [usize; 3], n: &[usize; 3]) -> usize {
    c[0] + n[0] * (c[1] + n[1] * c[2])
}

pub fn build_meshfree_neighbor_graph(
    points: &[Vec3],
    config: &MeshfreeNeighborGraphConfig,
) -> Result<MeshfreeNeighborGraph> {
    require(
        config.support_radius > 0.0
            && config.support_radius.is_finite()
            && config.maximum_cells > 0,
        "invalid meshfree neighbor configuration",
    )?;
    let mut g = MeshfreeNeighborGraph {
        offsets: vec![0; points.len() + 1],
        ..Default::default()
    };
    if points.is_empty() {
        return Ok(g);
    }
    let mut lo = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in points {
        require(
            p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
            "meshfree points must be finite",
        )?;
        lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    g.origin = lo;
    g.cell_size = config.support_radius;
    let count = |extent: f64| -> Result<usize> {
        let c = (extent / config.support_radius).ceil() + 1.0;
        require(
            c >= 1.0 && c <= config.maximum_cells as f64,
            "meshfree grid exceeds maximum cell count",
        )?;
        Ok(c as usize)
    };
    g.cell_counts = [
        count(hi.x - lo.x)?,
        count(hi.y - lo.y)?,
        count(hi.z - lo.z)?,
    ];
    let total = g.cell_counts[0]
        .checked_mul(g.cell_counts[1])
        .and_then(|v| v.checked_mul(g.cell_counts[2]))
        .ok_or_else(|| {
            crate::SolverError::InvalidInput("meshfree grid cell count overflows".into())
        })?;
    require(
        total <= config.maximum_cells,
        "meshfree grid exceeds maximum cell count",
    )?;
    let cells: Vec<[i64; 3]> = points
        .iter()
        .map(|p| point_cell(p, &g.origin, g.cell_size, &g.cell_counts))
        .collect();
    let flat: Vec<usize> = cells
        .iter()
        .map(|c| {
            flatten(
                [c[0] as usize, c[1] as usize, c[2] as usize],
                &g.cell_counts,
            )
        })
        .collect();
    let mut cell_offsets = vec![0usize; total + 1];
    for f in &flat {
        cell_offsets[f + 1] += 1;
    }
    for c in 0..total {
        cell_offsets[c + 1] += cell_offsets[c];
    }
    let mut cursor = cell_offsets.clone();
    let mut cell_points = vec![0usize; points.len()];
    for (p, f) in flat.iter().enumerate() {
        cell_points[cursor[*f]] = p;
        cursor[*f] += 1;
    }
    let r2 = config.support_radius * config.support_radius;
    for (p, cell) in cells.iter().enumerate() {
        let mut row = vec![];
        for dz in -1..=1i64 {
            let z = cell[2] + dz;
            if z < 0 || z >= g.cell_counts[2] as i64 {
                continue;
            }
            for dy in -1..=1i64 {
                let y = cell[1] + dy;
                if y < 0 || y >= g.cell_counts[1] as i64 {
                    continue;
                }
                for dx in -1..=1i64 {
                    let x = cell[0] + dx;
                    if x < 0 || x >= g.cell_counts[0] as i64 {
                        continue;
                    }
                    let f = flatten([x as usize, y as usize, z as usize], &g.cell_counts);
                    for &cand in &cell_points[cell_offsets[f]..cell_offsets[f + 1]] {
                        if !config.include_self && cand == p {
                            continue;
                        }
                        let (a, b) = (points[cand], points[p]);
                        let (ddx, ddy, ddz) = (a.x - b.x, a.y - b.y, a.z - b.z);
                        if ddx * ddx + ddy * ddy + ddz * ddz <= r2 {
                            row.push(cand);
                        }
                    }
                }
            }
        }
        row.sort_unstable();
        g.offsets[p] = g.indices.len();
        g.indices.extend(row);
    }
    g.offsets[points.len()] = g.indices.len();
    Ok(g)
}

/// Rank-local owned points and the sorted halo they reference.
pub fn classify_meshfree_halo(
    graph: &MeshfreeNeighborGraph,
    owner_rank: &[i32],
    local_rank: i32,
) -> Result<MeshfreeHaloDescriptor> {
    require(
        owner_rank.len() == graph.point_count() && local_rank >= 0,
        "meshfree ownership size or rank mismatch",
    )?;
    require(
        owner_rank.iter().all(|r| *r >= 0),
        "meshfree owner rank must be nonnegative",
    )?;
    let mut d = MeshfreeHaloDescriptor {
        owned: (0..graph.point_count())
            .filter(|p| owner_rank[*p] == local_rank)
            .collect(),
        halo: vec![],
    };
    let mut marked = vec![false; graph.point_count()];
    for &p in &d.owned {
        for &n in graph.neighbors(p)? {
            require(n < marked.len(), "corrupt meshfree graph")?;
            if owner_rank[n] != local_rank && !marked[n] {
                marked[n] = true;
                d.halo.push(n);
            }
        }
    }
    d.halo.sort_unstable();
    Ok(d)
}
