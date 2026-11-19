//! Engineering device models compiled into plain case fields, so they run on any producer
//! that has the underlying capability:
//!
//! - **PCB**: a layer stack becomes an anisotropic solid conductivity (in-plane layers in
//!   parallel, through-plane in series).
//! - **LED**: the heat is electrical power not emitted as light; an optional junction-to-case
//!   resistance makes it a two-resistor component.
//! - **Battery**: Joule heat I²R plus the reversible (entropic) heat I·T·dU/dT.
//! - **Perforated plate**: a one-cell porous layer on a domain face whose inertial resistance
//!   reproduces the loss coefficient of a thin sharp-edged plate with the given free-area
//!   ratio (Idelchik).
//!
//! Each returns a [`crate::case_scene::Edit`] for the draft: the derived numbers are written
//! where the solver reads them, and the record is pointed at for review.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::case_scene::{Edit, PATCHES, face_axes, finish, grid};
use crate::{Error, MAX_CASE_BYTES, Result, validate_json};

fn bad(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// Copper's conductivity, W/m·K.
pub const COPPER_W_PER_M_K: f64 = 390.0;
/// FR-4's conductivity, W/m·K.
pub const FR4_W_PER_M_K: f64 = 0.3;

/// One copper layer of a board: its thickness and the fraction of the board it covers.
#[derive(Clone, Copy, Debug, PartialEq, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CopperLayer {
    pub thickness_um: f64,
    /// 0..1 (a plane is about 0.9, a signal layer 0.2–0.4).
    pub coverage: f64,
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PcbStack {
    pub board_thickness_mm: f64,
    pub layers: Vec<CopperLayer>,
    #[serde(default = "fr4")]
    pub dielectric_w_per_m_k: f64,
    /// The board's normal axis; by default the part's thinnest direction.
    #[serde(default)]
    pub normal_axis: Option<usize>,
}

fn fr4() -> f64 {
    FR4_W_PER_M_K
}

/// In-plane and through-plane conductivity of a board.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, JsonSchema)]
pub struct PcbConductivity {
    pub in_plane_w_per_m_k: f64,
    pub through_plane_w_per_m_k: f64,
    /// Copper volume fraction.
    pub copper_fraction: f64,
}

/// Effective conductivity of a laminate: each copper layer is a mix of copper (its coverage)
/// and dielectric; the dielectric fills the rest of the thickness. In-plane, the layers act in
/// parallel; through the board, in series.
pub fn pcb_conductivity(stack: &PcbStack) -> Result<PcbConductivity> {
    let total = stack.board_thickness_mm * 1e-3;
    if !positive(total) || !positive(stack.dielectric_w_per_m_k) || stack.layers.len() > 64 {
        return Err(bad("a board needs a positive thickness and dielectric conductivity, and at most 64 copper layers"));
    }
    let kd = stack.dielectric_w_per_m_k;
    let mut copper_thickness = 0.0;
    let (mut parallel, mut series, mut copper) = (0.0, 0.0, 0.0);
    for layer in &stack.layers {
        let t = layer.thickness_um * 1e-6;
        if !positive(t) || !(0.0..=1.0).contains(&layer.coverage) {
            return Err(bad("copper layers need a positive thickness and a coverage from 0 to 1"));
        }
        let k = layer.coverage * COPPER_W_PER_M_K + (1.0 - layer.coverage) * kd;
        copper_thickness += t;
        parallel += k * t;
        series += t / k;
        copper += layer.coverage * t;
    }
    let dielectric = total - copper_thickness;
    if dielectric < 0.0 {
        return Err(bad("the copper layers are thicker than the board"));
    }
    parallel += kd * dielectric;
    series += dielectric / kd;
    Ok(PcbConductivity { in_plane_w_per_m_k: parallel / total, through_plane_w_per_m_k: total / series, copper_fraction: copper / total })
}

/// The thermal-solid record of `part`, created with starter values if it has none.
fn solid<'a>(case: &'a mut Value, part: &str) -> Result<(&'a mut serde_json::Map<String, Value>, String)> {
    let known = ["obstacles", "bodies"]
        .iter()
        .any(|key| case.get(*key).and_then(Value::as_array).is_some_and(|rows| rows.iter().any(|row| row["name"].as_str() == Some(part))));
    if !known {
        return Err(bad(format!("{part} is not an obstacle or surface body")));
    }
    let solids = case
        .pointer_mut("/thermal/solids")
        .and_then(Value::as_array_mut)
        .ok_or_else(|| bad("this device needs a thermal block (conjugate heat transfer)"))?;
    let index = match solids.iter().position(|row| row["obstacle"].as_str() == Some(part)) {
        Some(index) => index,
        None => {
            solids.push(json!({"obstacle": part, "conductivity_w_per_m_k": [1.0, 1.0, 1.0], "power_w": 0.0}));
            solids.len() - 1
        }
    };
    let row = solids.get_mut(index).and_then(Value::as_object_mut).ok_or_else(|| bad("thermal solids must be objects"))?;
    Ok((row, format!("/thermal/solids/{index}")))
}

/// Cell extent of grid part `part` (obstacle), or the vertex extent of a surface body.
fn part_extent(case: &Value, part: &str) -> Option<[f64; 3]> {
    let grid_nodes = crate::case_scene::nodes(case).ok()?;
    if let Some(row) = case.get("obstacles").and_then(Value::as_array)?.iter().find(|row| row["name"].as_str() == Some(part)) {
        let (min, max) = (crate::case_scene::index3(&row["min"])?, crate::case_scene::index3(&row["max"])?);
        let span = |a: usize| Some(grid_nodes[a].get(max[a])? - grid_nodes[a].get(min[a])?);
        return Some([span(0)?, span(1)?, span(2)?]);
    }
    let row = case.get("bodies").and_then(Value::as_array)?.iter().find(|row| row["name"].as_str() == Some(part))?;
    let vertices: Vec<[f64; 3]> = row["vertices"].as_array()?.iter().filter_map(|v| serde_json::from_value(v.clone()).ok()).collect();
    let span = |a: usize| {
        let (lo, hi) = vertices.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), v| (l.min(v[a]), h.max(v[a])));
        hi - lo
    };
    (!vertices.is_empty()).then(|| [span(0), span(1), span(2)])
}

/// Give `part` the conductivity of a board stack (power and other fields are kept).
pub fn apply_pcb(case: &Value, part: &str, stack: &PcbStack) -> Result<(Edit, PcbConductivity)> {
    validate_json(case, MAX_CASE_BYTES)?;
    let k = pcb_conductivity(stack)?;
    let normal = match stack.normal_axis {
        Some(axis) if axis < 3 => axis,
        Some(_) => return Err(bad("the board normal must be axis 0, 1 or 2")),
        None => {
            let extent = part_extent(case, part).ok_or_else(|| bad(format!("{part} is not a part of a 3-D grid case")))?;
            (0..3).min_by(|&a, &b| extent[a].total_cmp(&extent[b])).unwrap_or(2)
        }
    };
    let mut proposed = case.clone();
    let (row, pointer) = solid(&mut proposed, part)?;
    let mut conductivity = [k.in_plane_w_per_m_k; 3];
    conductivity[normal] = k.through_plane_w_per_m_k;
    row.insert("conductivity_w_per_m_k".into(), json!(conductivity));
    Ok((finish(proposed, Some(pointer))?, k))
}

#[derive(Clone, Copy, Debug, PartialEq, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Led {
    pub electrical_power_w: f64,
    /// Fraction of the electrical power leaving as light (0..1; white LEDs about 0.3–0.45).
    pub optical_efficiency: f64,
    /// Junction-to-case resistance from the data sheet; gives a junction temperature.
    #[serde(default)]
    pub junction_to_case_k_per_w: Option<f64>,
    #[serde(default)]
    pub max_junction_temperature_k: Option<f64>,
}

/// Heat an LED leaves in the part: the electrical power not emitted as light.
pub fn led_heat(led: &Led) -> Result<f64> {
    if !positive(led.electrical_power_w) || !(0.0..1.0).contains(&led.optical_efficiency) {
        return Err(bad("an LED needs positive electrical power and an optical efficiency from 0 to below 1"));
    }
    Ok(led.electrical_power_w * (1.0 - led.optical_efficiency))
}

pub fn apply_led(case: &Value, part: &str, led: &Led) -> Result<(Edit, f64)> {
    validate_json(case, MAX_CASE_BYTES)?;
    let heat = led_heat(led)?;
    if led.junction_to_case_k_per_w.is_some_and(|r| !(r.is_finite() && r >= 0.0)) || led.max_junction_temperature_k.is_some_and(|t| !positive(t)) {
        return Err(bad("the junction resistance and temperature limit must be positive"));
    }
    let mut proposed = case.clone();
    let (row, pointer) = solid(&mut proposed, part)?;
    row.insert("power_w".into(), json!(heat));
    if let Some(r) = led.junction_to_case_k_per_w {
        row.insert("junction_to_case_k_per_w".into(), json!(r));
    }
    if let Some(limit) = led.max_junction_temperature_k {
        row.insert("max_temperature_k".into(), json!(limit));
    }
    Ok((finish(proposed, Some(pointer))?, heat))
}

#[derive(Clone, Copy, Debug, PartialEq, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Battery {
    /// Positive while discharging.
    pub current_a: f64,
    pub internal_resistance_ohm: f64,
    /// dU/dT of the open-circuit voltage; 0 leaves out the reversible heat.
    #[serde(default)]
    pub entropic_coefficient_v_per_k: f64,
    /// Cell temperature for the reversible term (default: the case's inflow temperature).
    #[serde(default)]
    pub temperature_k: Option<f64>,
    /// Number of identical cells in the part.
    #[serde(default = "one")]
    pub cells: u32,
}

fn one() -> u32 {
    1
}

/// Battery heat: cells × (I²R − I·T·dU/dT). Discharging (I > 0) with dU/dT > 0 absorbs heat.
pub fn battery_heat(battery: &Battery, temperature_k: f64) -> Result<f64> {
    let (i, r, s) = (battery.current_a, battery.internal_resistance_ohm, battery.entropic_coefficient_v_per_k);
    if !(i.is_finite() && r.is_finite() && r >= 0.0 && s.is_finite() && positive(temperature_k)) || battery.cells == 0 || battery.cells > 100_000 {
        return Err(bad("a battery needs a finite current, a non-negative resistance and 1..100000 cells"));
    }
    Ok(f64::from(battery.cells) * (i * i * r - i * temperature_k * s))
}

pub fn apply_battery(case: &Value, part: &str, battery: &Battery) -> Result<(Edit, f64)> {
    validate_json(case, MAX_CASE_BYTES)?;
    let temperature = battery.temperature_k.or_else(|| case.pointer("/thermal/inflow_temperature_k").and_then(Value::as_f64)).unwrap_or(293.15);
    let heat = battery_heat(battery, temperature)?;
    if heat < 0.0 {
        return Err(bad(format!(
            "the battery absorbs {:.3} W at this operating point; the solver takes only heat sources (lower the entropic term)",
            -heat
        )));
    }
    let mut proposed = case.clone();
    let (row, pointer) = solid(&mut proposed, part)?;
    row.insert("power_w".into(), json!(heat));
    Ok((finish(proposed, Some(pointer))?, heat))
}

/// Loss coefficient K (Δp = K ½ρu², u the approach velocity) of a thin sharp-edged perforated
/// plate with free-area ratio φ (Idelchik): K = (0.707 (1−φ)^0.375 + 1 − φ)² / φ².
pub fn perforated_loss(free_area: f64) -> Result<f64> {
    if !(free_area > 0.0 && free_area < 1.0) {
        return Err(bad("the free-area ratio must be between 0 and 1"));
    }
    let open = 1.0 - free_area;
    Ok((0.707 * open.powf(0.375) + open).powi(2) / (free_area * free_area))
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Plate {
    pub loss_coefficient: f64,
    pub inertial_resistance_per_m: f64,
}

/// A perforated plate on domain face `patch` (a [`PATCHES`] index): a one-cell porous layer
/// against that face, over the whole face or `[region_min, region_max)` cells along
/// [`face_axes`]. Its inertial resistance C = K / h gives the plate's loss across the layer.
pub fn apply_perforated_plate(
    case: &Value,
    name: &str,
    patch: usize,
    free_area: f64,
    region: Option<([usize; 2], [usize; 2])>,
) -> Result<(Edit, Plate)> {
    validate_json(case, MAX_CASE_BYTES)?;
    let (dimensions, lengths) = grid(case)?;
    PATCHES.get(patch).ok_or_else(|| bad("unknown domain face"))?;
    if name.is_empty() || name.len() > 64 || !name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-') || PATCHES.contains(&name) {
        return Err(bad("the plate needs a new ASCII name"));
    }
    let taken = ["obstacles", "porous_blocks", "bodies"]
        .iter()
        .any(|key| case.get(*key).and_then(Value::as_array).is_some_and(|rows| rows.iter().any(|row| row["name"].as_str() == Some(name))));
    if taken {
        return Err(bad(format!("{name} is already used")));
    }
    let loss_coefficient = perforated_loss(free_area)?;
    let axis = patch / 2;
    let _ = lengths;
    let grid_nodes = crate::case_scene::nodes(case)?;
    let edge = &grid_nodes[axis];
    let h = if patch.is_multiple_of(2) { edge[1] - edge[0] } else { edge[edge.len() - 1] - edge[edge.len() - 2] };
    let plate = Plate { loss_coefficient, inertial_resistance_per_m: loss_coefficient / h };
    let [u, v] = face_axes(axis);
    let (lo, hi) = region.unwrap_or(([0, 0], [dimensions[u], dimensions[v]]));
    if lo[0] >= hi[0] || lo[1] >= hi[1] || hi[0] > dimensions[u] || hi[1] > dimensions[v] {
        return Err(bad("the plate region must be a non-empty rectangle inside the face"));
    }
    let layer = if patch.is_multiple_of(2) { (0, 1) } else { (dimensions[axis] - 1, dimensions[axis]) };
    let (mut min, mut max) = ([0usize; 3], [0usize; 3]);
    (min[axis], max[axis]) = layer;
    (min[u], max[u], min[v], max[v]) = (lo[0], hi[0], lo[1], hi[1]);
    let mut proposed = case.clone();
    let object = proposed.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
    let rows = object.entry("porous_blocks").or_insert_with(|| json!([])).as_array_mut().ok_or_else(|| bad("porous_blocks must be an array"))?;
    rows.push(
        json!({"name": name, "min": min, "max": max, "viscous_resistance_per_m2": 0.0, "inertial_resistance_per_m": plate.inertial_resistance_per_m}),
    );
    let pointer = format!("/porous_blocks/{}", rows.len() - 1);
    Ok((finish(proposed, Some(pointer))?, plate))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enclosure() -> Value {
        serde_json::from_str(include_str!("../data/templates/flow-electronics-enclosure-cht.json")).unwrap_or(Value::Null)
    }

    #[test]
    fn board_stacks_mix_in_parallel_along_and_in_series_across() {
        // 1.6 mm board, two full 35 µm copper planes.
        let stack = PcbStack {
            board_thickness_mm: 1.6,
            layers: vec![CopperLayer { thickness_um: 35.0, coverage: 1.0 }; 2],
            dielectric_w_per_m_k: FR4_W_PER_M_K,
            normal_axis: None,
        };
        let k = pcb_conductivity(&stack).expect("stack");
        let (tc, td) = (70e-6, 1.6e-3 - 70e-6);
        assert!((k.in_plane_w_per_m_k - (390.0 * tc + 0.3 * td) / 1.6e-3).abs() < 1e-9);
        assert!((k.through_plane_w_per_m_k - 1.6e-3 / (tc / 390.0 + td / 0.3)).abs() < 1e-9);
        assert!(k.in_plane_w_per_m_k > 17.0 && k.through_plane_w_per_m_k < 0.32);
        assert!((k.copper_fraction - tc / 1.6e-3).abs() < 1e-12);
        // No copper: plain FR-4 both ways.
        let bare = pcb_conductivity(&PcbStack { layers: vec![], ..stack.clone() }).expect("bare");
        assert!((bare.in_plane_w_per_m_k - 0.3).abs() < 1e-12 && (bare.through_plane_w_per_m_k - 0.3).abs() < 1e-12);
        assert!(pcb_conductivity(&PcbStack { layers: vec![CopperLayer { thickness_um: 2000.0, coverage: 1.0 }], ..stack.clone() }).is_err());
        assert!(pcb_conductivity(&PcbStack { layers: vec![CopperLayer { thickness_um: 35.0, coverage: 1.5 }], ..stack.clone() }).is_err());
        // On the template's board (16 × 10 × 1 cells): the thin axis is z.
        let (edit, _) = apply_pcb(&enclosure(), "pcb_block", &stack).expect("apply");
        assert!(edit.issues.is_empty(), "{:?}", edit.issues);
        let row = edit.case.pointer(edit.pointer.as_deref().unwrap_or("")).cloned().unwrap_or(Value::Null);
        let c: [f64; 3] = serde_json::from_value(row["conductivity_w_per_m_k"].clone()).unwrap_or_default();
        assert!((c[0] - k.in_plane_w_per_m_k).abs() < 1e-12 && (c[2] - k.through_plane_w_per_m_k).abs() < 1e-12);
        assert_eq!(row["power_w"], 0.5, "power is kept");
        assert!(apply_pcb(&enclosure(), "nothing", &stack).is_err());
    }

    #[test]
    fn leds_and_batteries_become_heat_sources() {
        let led =
            Led { electrical_power_w: 3.0, optical_efficiency: 0.35, junction_to_case_k_per_w: Some(4.5), max_junction_temperature_k: Some(408.15) };
        assert!((led_heat(&led).expect("heat") - 1.95).abs() < 1e-12);
        let (edit, heat) = apply_led(&enclosure(), "processor_block", &led).expect("led");
        assert!(edit.issues.is_empty(), "{:?}", edit.issues);
        let row = edit.case.pointer(edit.pointer.as_deref().unwrap_or("")).cloned().unwrap_or(Value::Null);
        assert_eq!((row["power_w"].as_f64(), row["junction_to_case_k_per_w"].as_f64()), (Some(heat), Some(4.5)));
        assert!(led_heat(&Led { optical_efficiency: 1.0, ..led }).is_err());
        // 10 A through 20 mΩ: 2 W per cell; 4 cells, no entropic term.
        let battery = Battery { current_a: 10.0, internal_resistance_ohm: 0.02, entropic_coefficient_v_per_k: 0.0, temperature_k: None, cells: 4 };
        assert!((battery_heat(&battery, 300.0).expect("heat") - 8.0).abs() < 1e-12);
        // dU/dT = −0.2 mV/K at 300 K adds I·T·0.2e-3 = 0.6 W per cell.
        let entropic = Battery { entropic_coefficient_v_per_k: -0.2e-3, temperature_k: Some(300.0), ..battery };
        assert!((battery_heat(&entropic, 300.0).expect("heat") - 4.0 * 2.6).abs() < 1e-12);
        let (edit, heat) = apply_battery(&enclosure(), "power_block", &entropic).expect("battery");
        assert!(edit.issues.is_empty() && (heat - 10.4).abs() < 1e-12);
        // A strongly endothermic point is refused rather than written as a negative source.
        let cooling = Battery { entropic_coefficient_v_per_k: 1e-3, internal_resistance_ohm: 0.0, ..entropic };
        assert!(apply_battery(&enclosure(), "power_block", &cooling).is_err());
        // Without a thermal block there is nowhere to put the heat.
        let mut cold = enclosure();
        cold.as_object_mut().map(|o| o.remove("thermal"));
        assert!(apply_led(&cold, "processor_block", &led).is_err());
    }

    #[test]
    fn perforated_plates_reproduce_their_loss_coefficient() {
        // Idelchik's thin-plate values: about 4 at 50 % open, about 51 at 20 %.
        assert!((perforated_loss(0.5).expect("k") - (0.707 * 0.5f64.powf(0.375) + 0.5).powi(2) / 0.25).abs() < 1e-12);
        assert!((4.0..4.7).contains(&perforated_loss(0.5).expect("k")));
        assert!((48.0..56.0).contains(&perforated_loss(0.2).expect("k")));
        assert!(perforated_loss(0.0).is_err() && perforated_loss(1.0).is_err());
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        let (edit, plate) = apply_perforated_plate(&channel, "outlet_plate", 1, 0.5, None).expect("plate");
        assert!(edit.issues.is_empty(), "{:?}", edit.issues);
        let h = channel["mesh"]["lengths"]["x"].as_f64().unwrap_or(0.0) / 8.0;
        assert!((plate.inertial_resistance_per_m - plate.loss_coefficient / h).abs() < 1e-9);
        // Δp/ρ across the layer = ½ C u² h = K ½ u².
        let u = 0.3;
        assert!((0.5 * plate.inertial_resistance_per_m * u * u * h - plate.loss_coefficient * 0.5 * u * u).abs() < 1e-12);
        let block = &edit.case["porous_blocks"][0];
        assert_eq!((block["min"].clone(), block["max"].clone()), (json!([7, 0, 0]), json!([8, 6, 2])));
        let (part, _) = apply_perforated_plate(&channel, "inlet_plate", 0, 0.3, Some(([1, 0], [5, 2]))).expect("region");
        assert_eq!(part.case["porous_blocks"][0]["min"], json!([0, 1, 0]));
        assert!(apply_perforated_plate(&edit.case, "outlet_plate", 1, 0.5, None).is_err(), "names are unique");
        assert!(apply_perforated_plate(&channel, "p", 1, 0.5, Some(([0, 0], [9, 2]))).is_err());
        assert!(apply_perforated_plate(&channel, "left", 1, 0.5, None).is_err());
    }
}
