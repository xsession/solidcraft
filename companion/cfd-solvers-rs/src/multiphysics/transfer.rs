//! Conservative 1-D cell-average remapping and face-flux to cell-rate conversion.
use crate::{fvm::PolyMesh, require, Result};

#[derive(Clone, Debug, PartialEq)]
pub struct CellGrid1D {
    pub edge: Vec<f64>,
}

impl CellGrid1D {
    pub fn validate(&self) -> Result<()> {
        require(
            self.edge.len() >= 2,
            "1-D transfer grid requires at least one cell",
        )?;
        require(
            self.edge.iter().all(|e| e.is_finite()) && self.edge.windows(2).all(|w| w[1] > w[0]),
            "1-D transfer edges must be strictly increasing",
        )
    }
    pub fn cell_count(&self) -> usize {
        self.edge.len().saturating_sub(1)
    }
}

/// Remap cell averages by exact interval overlap; both grids must span the same interval.
pub fn conservative_cell_average_transfer(
    source: &CellGrid1D,
    source_average: &[f64],
    target: &CellGrid1D,
) -> Result<Vec<f64>> {
    source.validate()?;
    target.validate()?;
    require(
        source_average.len() == source.cell_count(),
        "source field/grid size mismatch",
    )?;
    let (s0, s1) = (source.edge[0], source.edge[source.edge.len() - 1]);
    let (t0, t1) = (target.edge[0], target.edge[target.edge.len() - 1]);
    let scale = 1f64.max(s0.abs()).max(s1.abs());
    require(
        (s0 - t0).abs() <= 1e-12 * scale && (s1 - t1).abs() <= 1e-12 * scale,
        "conservative 1-D transfer requires identical physical extent",
    )?;
    let mut out = vec![0.0; target.cell_count()];
    let mut first = 0usize;
    for (t, value) in out.iter_mut().enumerate() {
        let (left, right) = (target.edge[t], target.edge[t + 1]);
        while first + 1 < source.edge.len() && source.edge[first + 1] <= left {
            first += 1;
        }
        let mut integral = 0.0;
        let mut s = first;
        while s < source.cell_count() && source.edge[s] < right {
            let overlap = right.min(source.edge[s + 1]) - left.max(source.edge[s]);
            if overlap > 0.0 {
                integral += source_average[s] * overlap;
            }
            if source.edge[s + 1] >= right {
                break;
            }
            s += 1;
        }
        *value = integral / (right - left);
    }
    Ok(out)
}

/// Oriented, face-integrated fluxes to a per-cell volumetric rate. Interior faces cancel
/// exactly in the volume-weighted global sum.
pub fn conservative_face_flux_to_cell_rate(
    mesh: &PolyMesh,
    integrated_face_flux: &[f64],
) -> Result<Vec<f64>> {
    require(
        integrated_face_flux.len() == mesh.face_count(),
        "face flux/mesh size mismatch",
    )?;
    let mut rate = vec![0.0; mesh.cell_count()];
    for (face, flux) in mesh.faces().iter().zip(integrated_face_flux) {
        rate[face.owner] += flux / mesh.cells()[face.owner].volume;
        if let Some(n) = face.neighbour {
            rate[n] -= flux / mesh.cells()[n].volume;
        }
    }
    Ok(rate)
}
