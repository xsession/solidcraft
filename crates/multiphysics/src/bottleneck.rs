//! The heat bottleneck field: where conduction struggles. Heat flows down the temperature
//! gradient, `q = −k ∇T`; the product of the flux and the gradient it runs down,
//!
//! ```text
//! B = −q · ∇T = kx (∂T/∂x)² + ky (∂T/∂y)² + kz (∂T/∂z)²     (W·K/m³)
//! ```
//!
//! is large where much heat has to cross a short distance of poor conductor: a board under a
//! hot chip, a thin thermal interface, a plastic housing. It is how fast conduction degrades the
//! temperature there (`B / T²` is the local entropy generation), so lowering the peaks of `B`
//! (a better material, a wider path, a shorter one) is what cools the parts above them.
//!
//! Computed by SolidCraft from a run's `thermal.vtu` (cell temperatures and region indices) and
//! the run's own input case (each region's conductivity), on the result's rectilinear cells:
//! central differences between neighbouring cell centres, one-sided at the domain's edge.
//! Portable: no files, no UI.
use serde_json::Value;

use crate::vtk::{Association, BoxGrid, Field};

/// The derived field's name, as the plots list it.
pub const FIELD: &str = "heat_bottleneck_w_k_per_m3";

/// Each region's conductivity along x, y and z (W/(m·K)): the fluid's for region 0, then the
/// case's thermal solids in order, as the companion numbers its `region_index` (solid `i` is
/// region `i + 1`). A heat pipe's axial conductivity is the companion's own working, so a heat
/// pipe's region takes its material's.
fn conductivities(case: &Value) -> Option<Vec<[f64; 3]>> {
    let thermal = case.get("thermal")?;
    let fluid = thermal.get("fluid")?.get("conductivity_w_per_m_k")?.as_f64()?;
    let mut out = vec![[fluid; 3]];
    for solid in thermal.get("solids")?.as_array()? {
        let k = solid.get("conductivity_w_per_m_k")?;
        let axes = match k {
            Value::Array(values) => {
                let v: Vec<f64> = values.iter().filter_map(Value::as_f64).collect();
                [*v.first()?, *v.get(1)?, *v.get(2)?]
            }
            value => [value.as_f64()?; 3],
        };
        out.push(axes);
    }
    out.iter().flatten().all(|k| k.is_finite() && *k >= 0.0).then_some(out)
}

/// Add the bottleneck field to `grid`, a run's thermal result, with conductivities from the
/// run's input `case`. False, and the grid as it was, when the grid holds no temperatures or
/// regions or the case no thermal block that numbers them.
pub fn add(grid: &mut BoxGrid, case: &Value) -> bool {
    if grid.fields.iter().any(|f| f.name == FIELD) {
        return true;
    }
    let Some(k) = conductivities(case) else { return false };
    let scalar = |name: &str| grid.fields.iter().position(|f| f.name == name && f.components == 1);
    let (Some(temperature), Some(region)) = (scalar("temperature_k"), scalar("region_index")) else { return false };
    let n = grid.cells.len();
    // Each cell's place on the axes' lines of centres, and the cell at each place.
    let centre = |cell: &[[f64; 2]; 3], axis: usize| (cell[axis][0] + cell[axis][1]) * 0.5;
    let mut place = Vec::with_capacity(n);
    for cell in &grid.cells {
        let mut at = [0usize; 3];
        for (axis, slot) in at.iter_mut().enumerate() {
            let Ok(i) = grid.planes[axis].binary_search_by(|p| p.total_cmp(&centre(cell, axis))) else { return false };
            *slot = i;
        }
        place.push(at);
    }
    let index: std::collections::HashMap<[usize; 3], usize> = place.iter().enumerate().map(|(c, at)| (*at, c)).collect();
    let values = |field: usize| grid.fields.get(field).map(|f| f.values.as_slice()).unwrap_or(&[]);
    let (t, regions) = (values(temperature), values(region));
    let mut out = Vec::with_capacity(n);
    for (c, at) in place.iter().enumerate() {
        let Some(&r) = regions.get(c) else { return false };
        let Some(cell_k) = (r >= 0.0 && r.fract() == 0.0).then_some(r as usize).and_then(|r| k.get(r)) else { return false };
        let mut b = 0.0;
        for axis in 0..3 {
            let neighbour = |step: isize| {
                let mut next = *at;
                next[axis] = at[axis].checked_add_signed(step)?;
                index.get(&next).copied()
            };
            let point = |cell: usize| Some((centre(grid.cells.get(cell)?, axis), *t.get(cell)?));
            let here = point(c);
            let gradient = match (neighbour(-1).and_then(point), here, neighbour(1).and_then(point)) {
                (Some((x0, t0)), _, Some((x1, t1))) | (Some((x0, t0)), Some((x1, t1)), None) | (None, Some((x0, t0)), Some((x1, t1))) => {
                    (t1 - t0) / (x1 - x0)
                }
                _ => 0.0,
            };
            b += cell_k[axis] * gradient * gradient;
        }
        out.push(if b.is_finite() { b } else { 0.0 });
    }
    grid.fields.push(Field { name: FIELD.into(), association: Association::Cell, data_type: "Float64".into(), components: 1, values: out });
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// `nx` cells along x (1 mm each), one across y and z, with temperatures `t` and regions.
    fn row(t: &[f64], regions: &[f64]) -> BoxGrid {
        let cells: Vec<[[f64; 2]; 3]> = (0..t.len()).map(|i| [[i as f64 * 1e-3, (i + 1) as f64 * 1e-3], [0.0, 1e-3], [0.0, 1e-3]]).collect();
        let planes = [(0..t.len()).map(|i| (i as f64 + 0.5) * 1e-3).collect(), vec![0.5e-3], vec![0.5e-3]];
        let field = |name: &str, values: &[f64]| Field {
            name: name.into(),
            association: Association::Cell,
            data_type: "Float64".into(),
            components: 1,
            values: values.to_vec(),
        };
        BoxGrid {
            bounds: [[0.0, t.len() as f64 * 1e-3], [0.0, 1e-3], [0.0, 1e-3]],
            cells,
            planes,
            fields: vec![field("temperature_k", t), field("region_index", regions)],
        }
    }

    fn case() -> Value {
        json!({"thermal": {"fluid": {"conductivity_w_per_m_k": 0.025}, "solids": [
            {"obstacle": "Board", "conductivity_w_per_m_k": [0.3, 0.3, 0.3]},
            {"obstacle": "Chip", "conductivity_w_per_m_k": 150.0}
        ]}})
    }

    #[test]
    fn heat_through_a_poor_conductor_is_the_bottleneck() {
        // One watt through a 1 mm² bar, 1 mm per cell: a copper-like chip (150), a board (0.3)
        // and the chip again. The flux is the same everywhere; the gradient, and B = q²/k, is
        // 500 times larger in the board.
        let q = 1.0 / 1e-6; // W/m²
        let (chip, board) = (q / 150.0 * 1e-3, q / 0.3 * 1e-3); // K per cell
        let t = [400.0, 400.0 - chip, 400.0 - chip - board, 400.0 - 2.0 * chip - board];
        let mut grid = row(&t, &[2.0, 2.0, 1.0, 2.0]);
        assert!(add(&mut grid, &case()));
        let b = &grid.fields.last().unwrap();
        assert_eq!((b.name.as_str(), b.components), (FIELD, 1));
        // The board's cell: central difference over two cells, half chip and half board.
        let central = ((chip + board) / 2e-3).powi(2) * 0.3;
        assert!((b.values[2] - central).abs() < 1e-6 * central, "{} vs {central}", b.values[2]);
        // The chip's first cell, one-sided at the domain's edge: exactly q²/k.
        assert!((b.values[0] - q * q / 150.0).abs() < 1e-9 * q * q / 150.0);
        assert!(b.values[2] > 100.0 * b.values[0], "the board is the bottleneck: {:?}", b.values);
        assert_eq!(crate::results::unit_hint(FIELD), Some("W·K/m³"));
        // Added once.
        assert!(add(&mut grid, &case()));
        assert_eq!(grid.fields.iter().filter(|f| f.name == FIELD).count(), 1);
    }

    #[test]
    fn a_flat_field_has_no_bottleneck_and_anisotropy_counts() {
        let mut grid = row(&[300.0; 3], &[0.0, 1.0, 2.0]);
        assert!(add(&mut grid, &case()));
        assert!(grid.fields.last().unwrap().values.iter().all(|v| *v == 0.0));
        // Conduction only along x in a board that conducts 0 across: x is all that counts.
        let mut grid = row(&[310.0, 300.0], &[1.0, 1.0]);
        let along = json!({"thermal": {"fluid": {"conductivity_w_per_m_k": 0.025}, "solids": [
            {"obstacle": "Board", "conductivity_w_per_m_k": [2.0, 0.0, 0.0]}
        ]}});
        assert!(add(&mut grid, &along));
        let gradient: f64 = 10.0 / 1e-3;
        assert!((grid.fields.last().unwrap().values[0] - 2.0 * gradient * gradient).abs() < 1e-6);
    }

    #[test]
    fn a_result_or_case_without_what_it_needs_is_left_alone() {
        let mut grid = row(&[310.0, 300.0], &[1.0, 1.0]);
        grid.fields.retain(|f| f.name != "region_index");
        assert!(!add(&mut grid, &case()));
        let mut grid = row(&[310.0, 300.0], &[1.0, 7.0]);
        assert!(!add(&mut grid, &case()), "a region the case does not number");
        let before = grid.fields.len();
        assert!(!add(&mut grid, &json!({"thermal": null})));
        assert_eq!(grid.fields.len(), before);
    }
}
