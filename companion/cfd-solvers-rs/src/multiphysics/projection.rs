//! Field projection between Tri3 FEM nodes, FVM cell centres and structured grids.
use crate::{
    fem::{Mesh2D, Node2},
    fvm::PolyMesh,
    require, Result, SolverError,
};

fn barycentric(a: Node2, b: Node2, c: Node2, x: f64, y: f64) -> Option<[f64; 3]> {
    let det = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
    if det.abs() < 1e-30 {
        return None;
    }
    let w0 = ((b.y - c.y) * (x - c.x) + (c.x - b.x) * (y - c.y)) / det;
    let w1 = ((c.y - a.y) * (x - c.x) + (a.x - c.x) * (y - c.y)) / det;
    let w2 = 1.0 - w0 - w1;
    let eps = -1e-12;
    (w0 >= eps && w1 >= eps && w2 >= eps).then_some([w0, w1, w2])
}

/// Piecewise-linear interpolation of FEM nodal values to FVM cell centres.
pub fn fem_nodes_to_fvm_cells(
    source: &Mesh2D,
    nodal: &[f64],
    target: &PolyMesh,
) -> Result<Vec<f64>> {
    source.validate()?;
    require(
        nodal.len() == source.node_count(),
        "FEM projection size mismatch",
    )?;
    target
        .cells()
        .iter()
        .map(|cell| {
            let p = cell.center;
            source
                .triangles
                .iter()
                .find_map(|t| {
                    let w = barycentric(
                        source.nodes[t.node[0]],
                        source.nodes[t.node[1]],
                        source.nodes[t.node[2]],
                        p.x,
                        p.y,
                    )?;
                    Some(
                        w[0] * nodal[t.node[0]] + w[1] * nodal[t.node[1]] + w[2] * nodal[t.node[2]],
                    )
                })
                .ok_or_else(|| {
                    SolverError::InvalidInput(
                        "FVM cell centre lies outside the FEM projection domain".into(),
                    )
                })
        })
        .collect()
}

/// Inverse-distance reconstruction of FVM cell values at FEM nodes from the nearest
/// `neighbours` cells; exact for constants, and takes a coincident cell's value directly.
pub fn fvm_cells_to_fem_nodes(
    source: &PolyMesh,
    cell_values: &[f64],
    target: &Mesh2D,
    neighbours: usize,
) -> Result<Vec<f64>> {
    require(
        cell_values.len() == source.cell_count() && neighbours > 0,
        "FVM projection input mismatch",
    )?;
    let mut out = Vec::with_capacity(target.node_count());
    for n in &target.nodes {
        let mut distance: Vec<(f64, usize)> = source
            .cells()
            .iter()
            .enumerate()
            .map(|(c, cell)| {
                let d = cell.center;
                ((d.x - n.x).powi(2) + (d.y - n.y).powi(2) + d.z * d.z, c)
            })
            .collect();
        distance.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let Some(nearest) = distance.first() else {
            return Err(SolverError::InvalidInput("empty FVM source".into()));
        };
        if nearest.0 < 1e-28 {
            out.push(cell_values[nearest.1]);
            continue;
        }
        let (mut sum, mut weight) = (0.0, 0.0);
        for (d2, c) in distance.iter().take(neighbours) {
            let w = 1.0 / d2.max(1e-28);
            sum += w * cell_values[*c];
            weight += w;
        }
        out.push(sum / weight);
    }
    Ok(out)
}

/// Node values of a uniform structured grid over [0,lx]x[0,ly]x[0,lz] (x fastest).
#[derive(Clone, Debug, PartialEq)]
pub struct StructuredScalarGrid3D {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub lx: f64,
    pub ly: f64,
    pub lz: f64,
    pub values: Vec<f64>,
}

impl StructuredScalarGrid3D {
    pub fn validate(&self) -> Result<()> {
        require(
            self.nx >= 2
                && self.ny >= 2
                && self.nz >= 2
                && [self.lx, self.ly, self.lz]
                    .iter()
                    .all(|l| l.is_finite() && *l > 0.0)
                && Some(self.values.len())
                    == self
                        .nx
                        .checked_mul(self.ny)
                        .and_then(|v| v.checked_mul(self.nz)),
            "invalid structured scalar grid",
        )
    }
    /// Trilinear sample, clamped to the grid box.
    pub fn sample(&self, x: f64, y: f64, z: f64) -> Result<f64> {
        self.validate()?;
        let axis = |v: f64, l: f64, n: usize| {
            let g = v.clamp(0.0, l) / l * (n - 1) as f64;
            let i = (g as usize).min(n - 2);
            (i, g - i as f64)
        };
        let ((i, fx), (j, fy), (k, fz)) = (
            axis(x, self.lx, self.nx),
            axis(y, self.ly, self.ny),
            axis(z, self.lz, self.nz),
        );
        let mut sum = 0.0;
        for dk in 0..2 {
            for dj in 0..2 {
                for di in 0..2 {
                    let w = (if di == 1 { fx } else { 1.0 - fx })
                        * (if dj == 1 { fy } else { 1.0 - fy })
                        * (if dk == 1 { fz } else { 1.0 - fz });
                    sum += w * self.values[((k + dk) * self.ny + j + dj) * self.nx + i + di];
                }
            }
        }
        Ok(sum)
    }
}

pub fn structured_to_fvm_cells(
    source: &StructuredScalarGrid3D,
    target: &PolyMesh,
) -> Result<Vec<f64>> {
    target
        .cells()
        .iter()
        .map(|c| source.sample(c.center.x, c.center.y, c.center.z))
        .collect()
}

pub fn structured_to_fem_nodes(
    source: &StructuredScalarGrid3D,
    target: &Mesh2D,
    z: f64,
) -> Result<Vec<f64>> {
    target
        .nodes
        .iter()
        .map(|n| source.sample(n.x, n.y, z))
        .collect()
}
