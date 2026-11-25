//! Partial cells for conjugate heat transfer (`partial_cells_v1`).
//!
//! Without them a cell belongs wholly to the solid whose surface encloses its centre (a
//! staircase), so a part's power goes evenly into its cells whatever their size, and a part
//! thinner than a cell has no cells at all. With `thermal.partial_cells`, every cell is measured
//! on a lattice up to [`MAX_SUB`] times finer per axis, by the same exact winding rule the
//! staircase uses, and the energy equation takes each cell for what it holds:
//!
//! - **power and heat capacity** go by solid volume: a part's power is spread over its sub-cells
//!   wherever they lie, so the total is exact and its density does not depend on the cells'
//!   sizes; a cell's heat capacity is the sum of its sub-cells';
//! - **conductivity** of a cell holding more than one material is, along each axis, the mean of
//!   the two classic bounds over its sub-cells: columns along the axis in series, in parallel
//!   with each other, and layers across the axis in parallel, in series with each other. A
//!   layer of board in a cell conducts across the layer like a board, along it like a board in
//!   parallel with air;
//! - a part with no cell of its own (a `lumped` body, thinner than a cell) is read out over the
//!   cells its sub-cells lie in, or over the cell holding its box's centre when it is thinner
//!   than a sub-cell: its power arrives there and its temperature is theirs.
//!
//! Cells keep their staircase region: the flow, contact resistances and the component table's
//! regions are unchanged.
use super::bodies::SurfaceBody;
use super::PolyMeshSpec;
use crate::{require, Result};

/// Sub-cells per cell along each axis at most.
pub const MAX_SUB: usize = 8;
/// Most sub-cells of the measuring lattice in all.
const MAX_SUB_CELLS: usize = 16_000_000;
/// A sub-cell of no body.
const FLUID: u8 = 0;

/// What each cell of the grid holds.
pub(super) struct Partial {
    /// Per cell, per solid: its share of the cell's volume, from 0 to 1.
    pub shares: Vec<Vec<(usize, f64)>>,
    /// Per cell holding more than one material: its conductivity along x, y and z.
    pub conductivity: Vec<Option<[f64; 3]>>,
}

fn sub_nodes(coarse: &[f64], sub: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity((coarse.len() - 1) * sub + 1);
    for w in coarse.windows(2) {
        for s in 0..sub {
            out.push(w[0] + (w[1] - w[0]) * s as f64 / sub as f64);
        }
    }
    if let Some(last) = coarse.last() {
        out.push(*last);
    }
    out
}

/// Measure the grid `spec`. `cell_solid` is the staircase layout (grid boxes are whole cells and
/// keep it); `boxes` is how many of its solids are grid boxes, the surface `bodies` follow them.
/// `fluid_k` and `solid_k` (per solid) are the conductivities along each axis.
pub(super) fn measure(
    spec: &PolyMeshSpec,
    cell_solid: &[Option<usize>],
    boxes: usize,
    bodies: &[SurfaceBody],
    fluid_k: [f64; 3],
    solid_k: &[[f64; 3]],
) -> Result<Partial> {
    let dims = spec.dimensions;
    let cells: usize = dims.iter().product();
    require(cell_solid.len() == cells, "partial cells: layout size mismatch")?;
    require(boxes + bodies.len() == solid_k.len(), "partial cells: one conductivity per solid")?;
    require(boxes + bodies.len() < 255, "partial cells take at most 254 solids")?;
    let mut sub = MAX_SUB;
    while sub > 2 && cells.saturating_mul(sub * sub * sub) > MAX_SUB_CELLS {
        sub -= 1;
    }
    let coarse: [Vec<f64>; 3] = std::array::from_fn(|a| spec.nodes(a));
    let nodes: [Vec<f64>; 3] = std::array::from_fn(|a| sub_nodes(&coarse[a], sub));
    let fine_dims: [usize; 3] = std::array::from_fn(|a| dims[a] * sub);
    let fine = PolyMeshSpec {
        dimensions: fine_dims,
        lengths: spec.lengths,
        xy_shear: 0.0,
        coordinates: Some(super::MeshCoordinates { x: nodes[0].clone(), y: nodes[1].clone(), z: nodes[2].clone() }),
    };
    let fine_index = |at: [usize; 3]| (at[2] * fine_dims[1] + at[1]) * fine_dims[0] + at[0];
    let coarse_index = |at: [usize; 3]| {
        let c: [usize; 3] = std::array::from_fn(|a| at[a] / sub);
        (c[2] * dims[1] + c[1]) * dims[0] + c[0]
    };
    // Per sub-cell: FLUID, or 1 + its solid. Grid boxes fill their cells; bodies their sub-cells.
    let mut owner = vec![FLUID; cells * sub * sub * sub];
    for k in 0..fine_dims[2] {
        for j in 0..fine_dims[1] {
            for i in 0..fine_dims[0] {
                if let Some(Some(s)) = cell_solid.get(coarse_index([i, j, k])) {
                    if *s < boxes {
                        owner[fine_index([i, j, k])] = *s as u8 + 1;
                    }
                }
            }
        }
    }
    for (bi, body) in bodies.iter().enumerate() {
        let inside = body.voxelize(&fine)?;
        let id = (boxes + bi) as u8 + 1;
        for (slot, inside) in owner.iter_mut().zip(inside) {
            // A grid box keeps its cells; the first body keeps a shared sub-cell.
            if inside && *slot == FLUID {
                *slot = id;
            }
        }
    }
    let per_cell = (sub * sub * sub) as f64;
    let k_of = |o: u8, axis: usize| -> f64 {
        if o == FLUID {
            fluid_k[axis]
        } else {
            solid_k.get(o as usize - 1).map_or(fluid_k[axis], |k| k[axis])
        }
    };
    let mut shares = vec![Vec::new(); cells];
    let mut conductivity = vec![None; cells];
    let mut counts = vec![0u32; boxes + bodies.len() + 1];
    for (cell, (slot_shares, slot_k)) in shares.iter_mut().zip(conductivity.iter_mut()).enumerate() {
        let c = [cell % dims[0], (cell / dims[0]) % dims[1], cell / (dims[0] * dims[1])];
        let lo: [usize; 3] = std::array::from_fn(|a| c[a] * sub);
        let at = |d: [usize; 3]| owner[fine_index([lo[0] + d[0], lo[1] + d[1], lo[2] + d[2]])];
        counts.iter_mut().for_each(|n| *n = 0);
        for k in 0..sub {
            for j in 0..sub {
                for i in 0..sub {
                    counts[at([i, j, k]) as usize] += 1;
                }
            }
        }
        for (o, n) in counts.iter().enumerate().skip(1) {
            if *n > 0 {
                slot_shares.push((o - 1, f64::from(*n) / per_cell));
            }
        }
        if counts.iter().filter(|n| **n > 0).count() < 2 {
            continue;
        }
        // Mixed: the mean of columns in series (in parallel) and layers in parallel (in series).
        let mut k = [0.0; 3];
        for (axis, slot) in k.iter_mut().enumerate() {
            let [b, c2] = [(axis + 1) % 3, (axis + 2) % 3];
            let point = |along: usize, u: usize, v: usize| {
                let mut d = [0usize; 3];
                d[axis] = along;
                d[b] = u;
                d[c2] = v;
                d
            };
            let mut columns = 0.0;
            for u in 0..sub {
                for v in 0..sub {
                    let resistance: f64 = (0..sub).map(|w| 1.0 / k_of(at(point(w, u, v)), axis)).sum();
                    columns += sub as f64 / resistance;
                }
            }
            let columns = columns / (sub * sub) as f64;
            let mut resistance = 0.0;
            for w in 0..sub {
                let mean: f64 = (0..sub).flat_map(|u| (0..sub).map(move |v| (u, v))).map(|(u, v)| k_of(at(point(w, u, v)), axis)).sum::<f64>()
                    / (sub * sub) as f64;
                resistance += 1.0 / mean;
            }
            let layers = sub as f64 / resistance;
            let value = 0.5 * (columns + layers);
            *slot = if value.is_finite() { value } else { 0.0 };
        }
        *slot_k = Some(k);
    }
    Ok(Partial { shares, conductivity })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Vec3;

    fn block(name: &str, lo: [f64; 3], hi: [f64; 3]) -> SurfaceBody {
        let v = |i: usize| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }];
        let quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]];
        SurfaceBody {
            name: name.into(),
            vertices: (0..8).map(v).collect(),
            triangles: quads.iter().flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]).collect(),
            lumped: true,
        }
    }

    fn grid(n: [usize; 3]) -> PolyMeshSpec {
        PolyMeshSpec {
            dimensions: n,
            lengths: Vec3::new(n[0] as f64 * 0.01, n[1] as f64 * 0.01, n[2] as f64 * 0.01),
            xy_shear: 0.0,
            coordinates: None,
        }
    }

    /// A slab filling the lower half of a 10 mm cell along x: its share is a half, and across the
    /// slab (x) it conducts like the two halves in series, along it (y, z) like them in parallel.
    #[test]
    fn a_half_filled_cell_has_half_the_volume_and_series_and_parallel_conductivity() {
        let spec = grid([1, 1, 1]);
        let slab = block("slab", [0.0, 0.0, 0.0], [0.005, 0.01, 0.01]);
        let partial = measure(&spec, &[None], 0, &[slab], [1.0; 3], &[[9.0; 3]]).unwrap();
        assert_eq!(partial.shares[0].len(), 1);
        assert!((partial.shares[0][0].1 - 0.5).abs() < 1e-12, "{:?}", partial.shares[0]);
        let k = partial.conductivity[0].expect("a mixed cell");
        // Layers across x: both bounds give the series value 2 / (1/1 + 1/9) = 1.8.
        assert!((k[0] - 1.8).abs() < 1e-12, "{k:?}");
        // Along the slab: both bounds give the parallel value (1 + 9) / 2 = 5.
        assert!((k[1] - 5.0).abs() < 1e-12 && (k[2] - 5.0).abs() < 1e-12, "{k:?}");
    }

    /// A 2.5 mm plate inside one cell of a 3-cell row: a quarter of that cell, none of the others.
    /// Grid boxes keep their whole cells.
    #[test]
    fn a_thin_plate_and_a_grid_box_are_measured() {
        let spec = grid([3, 1, 1]);
        let plate = block("plate", [0.0125, 0.0, 0.0], [0.015, 0.01, 0.01]);
        let layout = [Some(0), None, None];
        let partial = measure(&spec, &layout, 1, &[plate], [0.026; 3], &[[200.0; 3], [0.3; 3]]).unwrap();
        assert_eq!(partial.shares[0], vec![(0, 1.0)]);
        assert!(partial.conductivity[0].is_none(), "a whole box cell keeps its own conductivity");
        // Sub-cells are 1.25 mm: centres at 13.125 and 14.375 mm lie in the plate, two of eight.
        assert_eq!(partial.shares[1].len(), 1);
        assert!((partial.shares[1][0].1 - 0.25).abs() < 1e-12, "{:?}", partial.shares[1]);
        assert!(partial.shares[2].is_empty());
        let k = partial.conductivity[1].expect("mixed");
        assert!(k[0] < 0.026 * 1.5 && k[1] > 0.07, "air in series across the plate, in parallel along it: {k:?}");
    }
}
