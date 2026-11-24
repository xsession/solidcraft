//! Mesh study: is the answer the grid's or the design's? The same study on three grids, each
//! finer than the last, and per quantity the grid convergence index of Celik et al. (2008),
//! "Procedure for estimation and reporting of uncertainty due to discretization in CFD
//! applications", J. Fluids Eng. 130(7): Richardson extrapolation with the observed order.
//!
//! For grids 1 (finest), 2 and 3 (coarsest) of N₁ > N₂ > N₃ cells in a volume V, the
//! representative sizes are hᵢ = (V/Nᵢ)^⅓, and r₂₁ = h₂/h₁, r₃₂ = h₃/h₂ should be above 1.3.
//! With ε₃₂ = φ₃ − φ₂, ε₂₁ = φ₂ − φ₁ and s = sign(ε₃₂/ε₂₁), the observed order p solves
//!
//! ```text
//! p = | ln|ε₃₂/ε₂₁| + q(p) | / ln r₂₁,   q(p) = ln((r₂₁ᵖ − s)/(r₃₂ᵖ − s))
//! ```
//!
//! and then φ_ext = (r₂₁ᵖ φ₁ − φ₂)/(r₂₁ᵖ − 1), the relative error e₂₁ = |(φ₁ − φ₂)/φ₁|, and the
//! fine-grid index GCI = 1.25 e₂₁/(r₂₁ᵖ − 1): the band the grid-independent answer is expected
//! in, as a share of the fine grid's value. When the changes alternate in sign (s < 0) the
//! convergence is oscillatory and the index less reliable; it is reported as such.
//!
//! Temperatures are compared as rises above the inflow: a share of a rise means something, a
//! share of an absolute temperature in kelvin does not.
//!
//! Portable: numbers in, numbers out.
use schemars::JsonSchema;
use serde::Serialize;

use crate::component_table::ComponentTable;

/// The least refinement ratio between neighbouring grids the procedure asks for.
pub const MIN_RATIO: f64 = 1.3;
/// Safety factor of the fine-grid index for three grids.
const SAFETY: f64 = 1.25;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Behaviour {
    /// The changes shrink with the same sign: Richardson extrapolation applies.
    Monotonic,
    /// They alternate in sign: the index is a rougher guide.
    Oscillatory,
    /// The fine and medium grids give the same value: nothing left to resolve.
    Converged,
    /// The change grows from the coarse pair to the fine one, or the order cannot be found: the
    /// grids are not yet in the range where the answer settles.
    Diverging,
}

/// One quantity on three grids.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Estimate {
    /// Fine, medium and coarse values.
    pub values: [f64; 3],
    pub behaviour: Behaviour,
    /// Observed order of accuracy.
    pub order: Option<f64>,
    /// The extrapolated, grid-independent value.
    pub extrapolated: Option<f64>,
    /// Fine-grid convergence index, % of the fine value.
    pub gci_percent: Option<f64>,
}

/// The ratios of the grids' representative sizes, fine to medium and medium to coarse, from their
/// cell counts (same domain).
pub fn ratios(cells: [f64; 3]) -> Option<[f64; 2]> {
    let ok = cells.iter().all(|n| n.is_finite() && *n > 0.0);
    let (r21, r32) = ((cells[0] / cells[1]).cbrt(), (cells[1] / cells[2]).cbrt());
    (ok && r21.is_finite() && r32.is_finite()).then_some([r21, r32])
}

/// The estimate for `values` (fine, medium, coarse) on grids of `cells` cells.
pub fn estimate(cells: [f64; 3], values: [f64; 3]) -> Option<Estimate> {
    let [r21, r32] = ratios(cells)?;
    if !(r21 > 1.0 && r32 > 1.0) || values.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let [f1, f2, f3] = values;
    let (e21, e32) = (f2 - f1, f3 - f2);
    let scale = f1.abs().max(f2.abs()).max(f3.abs()).max(1e-300);
    let none = |behaviour| Estimate { values, behaviour, order: None, extrapolated: None, gci_percent: None };
    if e21.abs() <= 1e-12 * scale {
        return Some(Estimate { values, behaviour: Behaviour::Converged, order: None, extrapolated: Some(f1), gci_percent: Some(0.0) });
    }
    if e32.abs() <= 1e-12 * scale {
        // The coarse pair agrees but the fine pair does not: no order to find.
        return Some(none(Behaviour::Diverging));
    }
    let s = (e32 / e21).signum();
    let ratio = (e32 / e21).abs();
    // Fixed-point iteration for p (Celik et al., eqs. 3a–3c).
    let mut p = ratio.ln() / r21.ln();
    let mut found = false;
    for _ in 0..200 {
        let q = ((r21.powf(p) - s) / (r32.powf(p) - s)).ln();
        let next = (ratio.ln() + q).abs() / r21.ln();
        if !next.is_finite() {
            break;
        }
        if (next - p).abs() < 1e-10 {
            p = next;
            found = true;
            break;
        }
        p = next;
    }
    if !found || !(p > 0.0 && p.is_finite()) || (s > 0.0 && ratio < 1.0) {
        return Some(none(Behaviour::Diverging));
    }
    let rp = r21.powf(p);
    let extrapolated = (rp * f1 - f2) / (rp - 1.0);
    let gci = (f1.abs() > 0.0).then(|| SAFETY * ((f1 - f2) / f1).abs() / (rp - 1.0) * 100.0);
    Some(Estimate {
        values,
        behaviour: if s > 0.0 { Behaviour::Monotonic } else { Behaviour::Oscillatory },
        order: Some(p),
        extrapolated: extrapolated.is_finite().then_some(extrapolated),
        gci_percent: gci.filter(|g| g.is_finite()),
    })
}

/// One part's highest-temperature rise above the inflow on the three grids.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct PartEstimate {
    pub body: String,
    /// The fine grid's highest temperature, °C.
    pub fine_c: f64,
    /// The rise above the inflow, K, on each grid and extrapolated.
    pub estimate: Estimate,
}

/// The mesh study of three runs' component tables (fine, medium, coarse) on grids of `cells`
/// cells, with the inflow at `inflow_k`: one estimate per part all three report.
pub fn parts(tables: [&ComponentTable; 3], cells: [f64; 3], inflow_k: f64) -> Vec<PartEstimate> {
    let rise = |table: &ComponentTable, body: &str| {
        table.rows.iter().find(|r| r.in_study && r.body == body).and_then(|r| r.result.as_ref()).map(|x| x.max_temperature_k - inflow_k)
    };
    tables[0]
        .rows
        .iter()
        .filter(|r| r.in_study)
        .filter_map(|row| {
            let values = [rise(tables[0], &row.body)?, rise(tables[1], &row.body)?, rise(tables[2], &row.body)?];
            Some(PartEstimate { body: row.body.clone(), fine_c: values[0] + inflow_k - 273.15, estimate: estimate(cells, values)? })
        })
        .collect()
}

/// A built case's base cell size (mm): its largest cell along any axis, which local
/// refinements leave as it was; and how many cells it has.
pub fn grid_of(case: &serde_json::Value) -> Option<(f64, f64)> {
    let nodes = crate::case_scene::nodes(case).ok()?;
    let base = nodes.iter().flat_map(|axis| axis.windows(2).map(|w| w[1] - w[0])).fold(0.0, f64::max) * 1000.0;
    let cells = nodes.iter().map(|axis| axis.len().saturating_sub(1) as f64).product::<f64>();
    (base.is_finite() && base > 0.0 && cells > 0.0).then_some((base, cells))
}

/// How much coarser each grid of a mesh study is than the one before, along each axis: √2,
/// so each has about 2.8 times fewer cells, past the procedure's 1.3 with room for whole cells.
pub const COARSENING: f64 = std::f64::consts::SQRT_2;

#[cfg(test)]
mod tests {
    use super::*;

    /// A quantity that converges at second order: φ(h) = 10 + 4 h².
    #[test]
    fn a_second_order_quantity_extrapolates_to_its_limit() {
        let h = [1.0f64, 2f64.sqrt(), 2.0];
        let cells = h.map(|h| 1000.0 / (h * h * h));
        let values = h.map(|h| 10.0 + 4.0 * h * h);
        let e = estimate(cells, values).expect("estimate");
        assert_eq!(e.behaviour, Behaviour::Monotonic);
        assert!((e.order.unwrap() - 2.0).abs() < 1e-6, "{e:?}");
        assert!((e.extrapolated.unwrap() - 10.0).abs() < 1e-6, "{e:?}");
        // GCI = 1.25 · |(14 − 18)/14| / (2 − 1) = 35.7 %.
        assert!((e.gci_percent.unwrap() - 1.25 * 4.0 / 14.0 * 100.0).abs() < 1e-6, "{e:?}");
        assert!(ratios(cells).unwrap().iter().all(|r| (r - 2f64.sqrt()).abs() < 1e-9));
    }

    #[test]
    fn unequal_ratios_still_find_the_order() {
        // φ = 5 − 3 h^1.5 on h = 1, 1.5, 2.6.
        let h = [1.0f64, 1.5, 2.6];
        let cells = h.map(|h| 1.0 / (h * h * h));
        let values = h.map(|h| 5.0 - 3.0 * h.powf(1.5));
        let e = estimate(cells, values).expect("estimate");
        assert!((e.order.unwrap() - 1.5).abs() < 1e-6, "{e:?}");
        assert!((e.extrapolated.unwrap() - 5.0).abs() < 1e-6, "{e:?}");
    }

    #[test]
    fn oscillating_converged_and_diverging_grids_are_named() {
        let cells = [8000.0, 2828.4, 1000.0];
        let wobble = estimate(cells, [10.0, 10.4, 9.0]).expect("wobble");
        assert_eq!(wobble.behaviour, Behaviour::Oscillatory);
        assert_eq!(estimate(cells, [7.0, 7.0, 6.0]).unwrap().behaviour, Behaviour::Converged);
        // The change grows on the finer pair: not settling yet.
        assert_eq!(estimate(cells, [10.0, 12.0, 12.5]).unwrap().behaviour, Behaviour::Diverging);
        assert!(estimate([1000.0, 1000.0, 500.0], [1.0, 2.0, 3.0]).is_none(), "no refinement between the fine grids");
    }
}
