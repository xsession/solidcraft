//! Solution-adaptive refinement: where a finished run's field changes most from cell to cell,
//! as boxes to refine before the next run.
//!
//! The indicator of a cell is the largest jump of the field (a vector's magnitude) to one of its
//! six face neighbours on the run's grid. The cells with the largest indicators, `fraction` of
//! them, are flagged; flagged cells that touch (also along edges and corners) form a cluster, and
//! each of the largest `max_regions` clusters gives its box. Cells of the same grid position (the
//! pieces of a cut cell) count as one, with their mean value.
use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::vtk::{BoxGrid, parse_vtu};
use crate::{Error, Result};

/// Most boxes one adaptation gives.
pub const MAX_REGIONS: usize = 16;

/// What an adaptation found.
#[derive(Clone, Debug, PartialEq)]
pub struct Adapted {
    /// The field it read.
    pub field: String,
    /// Grid positions with a value, and how many of them were flagged.
    pub cells: usize,
    pub flagged: usize,
    /// The clusters' boxes in the run's coordinates (case metres), largest cluster first.
    pub regions_m: Vec<([f64; 3], [f64; 3])>,
}

fn bad(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

type Key = [usize; 3];

/// The boxes to refine from the run field file `vtu`: `field` by name (default the velocity,
/// else the first field), the `fraction` (0 to 0.5] of cells with the largest jumps.
pub fn regions(vtu: &str, field: Option<&str>, fraction: f64, max_regions: usize) -> Result<Adapted> {
    if !(fraction.is_finite() && fraction > 0.0 && fraction <= 0.5) {
        return Err(bad("refine a fraction above 0 and at most 0.5 of the cells"));
    }
    if !(1..=MAX_REGIONS).contains(&max_regions) {
        return Err(bad(format!("give 1 to {MAX_REGIONS} regions")));
    }
    let grid = parse_vtu(vtu).and_then(BoxGrid::from_grid).map_err(|e| bad(format!("the run's field file: {}", e.0)))?;
    let index = match field {
        Some(name) => grid.fields.iter().position(|f| f.name == name).ok_or_else(|| bad(format!("the run has no field {name}")))?,
        None => grid.fields.iter().position(|f| f.name == "velocity_m_per_s").unwrap_or(0),
    };
    let name = grid.fields.get(index).map(|f| f.name.clone()).unwrap_or_default();
    // Each cell's grid position: its centre among the sorted centres along each axis.
    let mut sum: BTreeMap<Key, (f64, usize, [[f64; 2]; 3])> = BTreeMap::new();
    for (c, cell) in grid.cells.iter().enumerate() {
        let Some(value) = grid.value(index, c) else { continue };
        let mut key = [0usize; 3];
        for a in 0..3 {
            let centre = 0.5 * (cell[a][0] + cell[a][1]);
            key[a] = grid.planes[a].partition_point(|p| *p < centre);
        }
        let entry = sum.entry(key).or_insert((0.0, 0, *cell));
        entry.0 += value;
        entry.1 += 1;
    }
    let value: BTreeMap<Key, f64> = sum.iter().map(|(k, (s, n, _))| (*k, s / *n as f64)).collect();
    let neighbour = |k: Key, a: usize, step: isize| -> Option<Key> {
        let mut n = k;
        n[a] = k[a].checked_add_signed(step)?;
        value.contains_key(&n).then_some(n)
    };
    let jump: BTreeMap<Key, f64> = value
        .iter()
        .map(|(k, v)| {
            let largest = (0..3)
                .flat_map(|a| [neighbour(*k, a, -1), neighbour(*k, a, 1)])
                .flatten()
                .filter_map(|n| value.get(&n))
                .map(|w| (v - w).abs())
                .fold(0.0, f64::max);
            (*k, largest)
        })
        .collect();
    let mut sorted: Vec<f64> = jump.values().copied().collect();
    sorted.sort_by(|a, b| b.total_cmp(a));
    let wanted = ((fraction * sorted.len() as f64).ceil() as usize).clamp(1, sorted.len().max(1));
    // Fewer cells than the fraction may change at all: then every cell that changes.
    let threshold = sorted
        .get(wanted - 1)
        .copied()
        .filter(|t| *t > 0.0)
        .or_else(|| sorted.iter().copied().rfind(|t| *t > 0.0))
        .ok_or_else(|| bad(format!("{name} is the same in every cell; there is nothing to refine")))?;
    let flagged: BTreeSet<Key> = jump.iter().filter(|(_, j)| **j >= threshold).map(|(k, _)| *k).collect();
    // Clusters of flagged cells touching along faces, edges or corners.
    let mut seen = BTreeSet::new();
    let mut clusters: Vec<(usize, [f64; 3], [f64; 3])> = Vec::new();
    for start in &flagged {
        if !seen.insert(*start) {
            continue;
        }
        let mut queue = VecDeque::from([*start]);
        let (mut count, mut low, mut high) = (0usize, [f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        while let Some(k) = queue.pop_front() {
            count += 1;
            if let Some((_, _, cell)) = sum.get(&k) {
                for a in 0..3 {
                    low[a] = low[a].min(cell[a][0]);
                    high[a] = high[a].max(cell[a][1]);
                }
            }
            for dx in -1isize..=1 {
                for dy in -1isize..=1 {
                    for dz in -1isize..=1 {
                        let (Some(x), Some(y), Some(z)) = (k[0].checked_add_signed(dx), k[1].checked_add_signed(dy), k[2].checked_add_signed(dz))
                        else {
                            continue;
                        };
                        let n = [x, y, z];
                        if flagged.contains(&n) && seen.insert(n) {
                            queue.push_back(n);
                        }
                    }
                }
            }
        }
        clusters.push((count, low, high));
    }
    clusters.sort_by_key(|c| std::cmp::Reverse(c.0));
    Ok(Adapted {
        field: name,
        cells: value.len(),
        flagged: flagged.len(),
        regions_m: clusters.into_iter().take(max_regions).map(|(_, low, high)| (low, high)).collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An `nx × ny × nz` grid of 0.1 m cells with one scalar `f` from `value(i, j, k)`.
    fn vtu(n: [usize; 3], value: impl Fn(usize, usize, usize) -> f64) -> String {
        let (mut points, mut connectivity, mut offsets, mut values) = (String::new(), String::new(), String::new(), String::new());
        let mut count = 0;
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    for (dx, dy, dz) in [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)] {
                        points.push_str(&format!("{} {} {} ", (i + dx) as f64 * 0.1, (j + dy) as f64 * 0.1, (k + dz) as f64 * 0.1));
                    }
                    for v in 0..8 {
                        connectivity.push_str(&format!("{} ", count * 8 + v));
                    }
                    count += 1;
                    offsets.push_str(&format!("{} ", count * 8));
                    values.push_str(&format!("{} ", value(i, j, k)));
                }
            }
        }
        let types = "12 ".repeat(count);
        format!(
            r#"<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian"><UnstructuredGrid><Piece NumberOfPoints="{}" NumberOfCells="{count}"><CellData><DataArray type="Float64" Name="f" NumberOfComponents="1" format="ascii">{values}</DataArray></CellData><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{points}</DataArray></Points><Cells><DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray><DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray><DataArray type="UInt8" Name="types" format="ascii">{types}</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>"#,
            count * 8
        )
    }

    #[test]
    fn a_step_in_the_field_is_refined_where_it_jumps() {
        // A step between x cells 4 and 5 of a 10 × 2 × 2 grid: the cells on both sides of it
        // jump most, and they form one box from x = 0.4 to 0.6 m through the whole y-z section.
        let step = vtu([10, 2, 2], |i, _, _| if i < 5 { 0.0 } else { 1.0 });
        let adapted = regions(&step, Some("f"), 0.2, 4).expect("adapt");
        assert_eq!((adapted.cells, adapted.flagged), (40, 8));
        assert_eq!(adapted.regions_m.len(), 1);
        let (low, high) = adapted.regions_m[0];
        for (got, want) in low.iter().zip([0.4, 0.0, 0.0]).chain(high.iter().zip([0.6, 0.2, 0.2])) {
            assert!((got - want).abs() < 1e-12, "{low:?} {high:?}");
        }
        // Fewer changing cells than the fraction: just those.
        let narrow = regions(&step, Some("f"), 0.5, 4).expect("adapt");
        assert_eq!(narrow.flagged, 8);
        // Two steps far apart are two clusters.
        let two = vtu([20, 1, 1], |i, _, _| {
            if i < 3 {
                0.0
            } else if i < 15 {
                1.0
            } else {
                3.0
            }
        });
        let adapted = regions(&two, None, 0.2, 4).expect("adapt");
        assert_eq!(adapted.field, "f");
        assert_eq!(adapted.regions_m.len(), 2);
        let mut starts: Vec<f64> = adapted.regions_m.iter().map(|r| r.0[0]).collect();
        starts.sort_by(f64::total_cmp);
        assert!((starts[0] - 0.2).abs() < 1e-12 && (starts[1] - 1.4).abs() < 1e-12, "{:?}", adapted.regions_m);
    }

    #[test]
    fn uniform_fields_and_bad_settings_are_refused() {
        let flat = vtu([4, 1, 1], |_, _, _| 2.0);
        assert!(regions(&flat, None, 0.1, 4).unwrap_err().to_string().contains("nothing to refine"));
        assert!(regions(&flat, None, 0.9, 4).is_err());
        assert!(regions(&flat, None, 0.1, 0).is_err());
        assert!(regions(&flat, Some("speed"), 0.1, 4).unwrap_err().to_string().contains("no field speed"));
        assert!(regions("not xml", None, 0.1, 4).is_err());
    }
}
