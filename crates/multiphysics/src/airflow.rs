//! Admission checks for the external fan curve and porous resistance zones of the
//! primitive airflow domain, plus the producer capabilities each authored feature needs.
//! These are structural bounds. They do not establish a fan's operating point, a filter's
//! characterization or any thermal effect; the companion computes and reports those.
use crate::{Error, Result};
use serde::Deserialize;
use serde_json::Value;
use std::collections::BTreeSet;

const EXTERIOR: [&str; 6] = ["left", "right", "bottom", "top", "front", "back"];
const MAX_CURVE_POINTS: usize = 64;
const MAX_POROUS_BLOCKS: usize = 32;

/// A producer extension an authored case depends on, and the producer versions that provide it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Requirement {
    pub capability: &'static str,
    pub versions: &'static [&'static str],
}

pub const GRID_OBSTACLES: Requirement = Requirement {
    capability: "grid_obstacles_v1",
    versions: &[
        "0.3.0", "0.4.0", "0.5.0", "0.6.0", "0.7.0", "0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0",
        "0.17.0", "0.18.0",
    ],
};
pub const FAN_CURVE: Requirement = Requirement {
    capability: "fan_curve_v1",
    versions: &[
        "0.4.0", "0.5.0", "0.6.0", "0.7.0", "0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0",
        "0.18.0",
    ],
};
pub const POROUS_BLOCKS: Requirement = Requirement {
    capability: "porous_blocks_v1",
    versions: &[
        "0.4.0", "0.5.0", "0.6.0", "0.7.0", "0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0",
        "0.18.0",
    ],
};
pub const CONJUGATE_HEAT: Requirement = Requirement {
    capability: "conjugate_heat_v1",
    versions: &[
        "0.5.0", "0.6.0", "0.7.0", "0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0",
    ],
};
pub const SURFACE_BODIES: Requirement = Requirement {
    capability: "surface_bodies_v1",
    versions: &["0.6.0", "0.7.0", "0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};
pub const GOALS: Requirement = Requirement {
    capability: "goals_v1",
    versions: &["0.6.0", "0.7.0", "0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};
pub const NATURAL_CONVECTION: Requirement = Requirement {
    capability: "natural_convection_v1",
    versions: &["0.6.0", "0.7.0", "0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};
/// Electrical, mechanics, particle, kinetic and coupled domain cases added in 0.7.0.
pub const DOMAIN_CASES: Requirement = Requirement {
    capability: "domain_cases_v1",
    versions: &["0.7.0", "0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};
/// Contact resistance, two-resistor junctions, heat pipes and thermoelectric coolers (0.8.0).
pub const THERMAL_DEVICES: Requirement = Requirement {
    capability: "thermal_devices_v1",
    versions: &["0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};
/// Global, point and equation goals (0.8.0).
/// Explicit rectilinear mesh coordinates (local mesh refinement, 0.9.0).
pub const LOCAL_MESH: Requirement = Requirement {
    capability: "local_mesh_v1",
    versions: &["0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};
/// Per-face transferred boundary profiles (zoom-in from a parent run, 0.10.0).
pub const TRANSFERRED: Requirement = Requirement {
    capability: "transferred_bc_v1",
    versions: &["0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};
/// Rotating regions, the multiple reference frame model (0.11.0).
pub const ROTATING: Requirement =
    Requirement { capability: "rotating_regions_v1", versions: &["0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"] };
/// PISO correctors that rebuild H(u)/a_P (stable at larger time steps, 0.11.0).
pub const PISO_UPDATE: Requirement = Requirement {
    capability: "piso_neighbour_update_v1",
    versions: &["0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};
/// Internal fans: actuator blocks following a pressure-flow curve (0.12.0).
pub const INTERNAL_FANS: Requirement =
    Requirement { capability: "internal_fans_v1", versions: &["0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"] };
/// Balanced-force pressure and body-force coupling (0.12.0).
pub const BALANCED_FORCES: Requirement =
    Requirement { capability: "balanced_body_force_v1", versions: &["0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"] };
/// Boundary patches on parts of surface bodies: conditions on model faces (0.13.0).
pub const BODY_FACE_CONDITIONS: Requirement =
    Requirement { capability: "body_face_conditions_v1", versions: &["0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"] };
/// Cut cells for surface bodies (0.14.0).
pub const CUT_CELLS: Requirement = Requirement { capability: "cut_cells_v1", versions: &["0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"] };
/// Power on parts of surface bodies, and heat transfer with closed body face conditions (0.15.0).
pub const FACE_HEAT_SOURCES: Requirement = Requirement { capability: "face_heat_sources_v1", versions: &["0.15.0", "0.16.0", "0.17.0", "0.18.0"] };
/// Global flow goals: the fluid's largest and mean speed and pressure (0.16.0).
pub const GOALS_V3: Requirement = Requirement { capability: "goals_v3", versions: &["0.16.0", "0.17.0", "0.18.0"] };
/// Partial cells for heat: power, heat capacity and conductivity by what each cell holds, and
/// lumped bodies with no cell of their own (0.17.0).
pub const PARTIAL_CELLS: Requirement = Requirement { capability: "partial_cells_v1", versions: &["0.17.0", "0.18.0"] };
/// Momentum solved to a tolerance relative to its own size, for small and slow flows.
pub const SCALED_MOMENTUM: Requirement = Requirement { capability: "scaled_momentum_v1", versions: &["0.18.0"] };
pub const GOALS_V2: Requirement = Requirement {
    capability: "goals_v2",
    versions: &["0.8.0", "0.9.0", "0.10.0", "0.11.0", "0.12.0", "0.13.0", "0.14.0", "0.15.0", "0.16.0", "0.17.0", "0.18.0"],
};

/// Goal kinds added with `goals_v3`.
pub const GOAL_KINDS_V3: &[&str] = &["global_max_speed", "global_mean_speed", "global_max_pressure", "global_min_pressure", "global_mean_pressure"];

/// Goal kinds added with `goals_v2`.
pub const GOAL_KINDS_V2: &[&str] = &[
    "global_max_fluid_temperature",
    "global_mean_fluid_temperature",
    "global_max_solid_temperature",
    "point_temperature",
    "point_pressure",
    "equation",
];

/// Thermal-solid fields added with `thermal_devices_v1`.
const SOLID_DEVICE_FIELDS: &[&str] = &["contact_resistance_m2_k_per_w", "junction_to_case_k_per_w", "heat_pipe"];

/// Solvers whose case adapters first appear in producer 0.7.0.
pub const DOMAIN_SOLVERS: &[&str] = &[
    "bar_modes1d",
    "battery_cell",
    "cahn_hilliard1d",
    "cosserat_cantilever",
    "dem_granular",
    "dg_advection1d",
    "edge_cavity3d",
    "elasticity2d",
    "galvanic_corrosion",
    "j2_uniaxial",
    "joule_heating2d",
    "landau_damping1d",
    "lennard_jones_md1d",
    "maxwell3d",
    "multiphase_box",
    "particle_transport_slab",
    "pdn_ir_drop",
    "pennes_bioheat2d",
    "phase_field_fracture1d",
    "photonic_waveguide",
    "pic_two_stream1d",
    "power_network",
    "spice_circuit",
    "spray_injection",
    "tcad_diode1d",
    "thmc_column1d",
    "wire_antenna",
];

fn nonempty_array(case: &Value, key: &str) -> bool {
    case.get(key).and_then(Value::as_array).is_some_and(|array| !array.is_empty())
}

/// Producer extensions required by the authored case, in a stable order.
pub fn requirements(case: &Value) -> Vec<Requirement> {
    let mut out = vec![];
    if case.get("solver").and_then(Value::as_str).is_some_and(|solver| DOMAIN_SOLVERS.contains(&solver)) {
        out.push(DOMAIN_CASES);
    }
    if nonempty_array(case, "obstacles") {
        out.push(GRID_OBSTACLES);
    }
    if case.get("fan").is_some_and(|fan| !fan.is_null()) {
        out.push(FAN_CURVE);
    }
    if nonempty_array(case, "porous_blocks") {
        out.push(POROUS_BLOCKS);
    }
    if case.get("thermal").is_some_and(|thermal| !thermal.is_null()) {
        out.push(CONJUGATE_HEAT);
    }
    if nonempty_array(case, "bodies") {
        out.push(SURFACE_BODIES);
    }
    if nonempty_array(case, "goals") {
        out.push(GOALS);
    }
    if case.get("thermal").and_then(|thermal| thermal.get("buoyancy")).is_some_and(|b| !b.is_null()) {
        out.push(NATURAL_CONVECTION);
    }
    let thermal = case.get("thermal").filter(|t| !t.is_null());
    let devices = thermal.is_some_and(|t| {
        nonempty_array(t, "thermoelectric_coolers")
            || t.get("solids")
                .and_then(Value::as_array)
                .is_some_and(|solids| solids.iter().any(|s| SOLID_DEVICE_FIELDS.iter().any(|f| s.get(*f).is_some_and(|v| !v.is_null()))))
    });
    if devices {
        out.push(THERMAL_DEVICES);
    }
    let goals_v2 = case
        .get("goals")
        .and_then(Value::as_array)
        .is_some_and(|goals| goals.iter().any(|g| g.get("kind").and_then(Value::as_str).is_some_and(|k| GOAL_KINDS_V2.contains(&k))));
    if goals_v2 {
        out.push(GOALS_V2);
    }
    let goals_v3 = case
        .get("goals")
        .and_then(Value::as_array)
        .is_some_and(|goals| goals.iter().any(|g| g.get("kind").and_then(Value::as_str).is_some_and(|k| GOAL_KINDS_V3.contains(&k))));
    if goals_v3 {
        out.push(GOALS_V3);
    }
    if case.pointer("/mesh/coordinates").is_some_and(|c| !c.is_null()) {
        out.push(LOCAL_MESH);
    }
    if case.get("transferred").and_then(Value::as_array).is_some_and(|rows| !rows.is_empty()) {
        out.push(TRANSFERRED);
    }
    if nonempty_array(case, "rotating_regions") {
        out.push(ROTATING);
    }
    if case.pointer("/config/piso_neighbour_update").and_then(Value::as_bool) == Some(true) {
        out.push(PISO_UPDATE);
    }
    if case.pointer("/config/scaled_momentum").and_then(Value::as_bool) == Some(true) {
        out.push(SCALED_MOMENTUM);
    }
    if nonempty_array(case, "internal_fans") {
        out.push(INTERNAL_FANS);
    }
    if case.pointer("/config/balanced_body_force").and_then(Value::as_bool) == Some(true) {
        out.push(BALANCED_FORCES);
    }
    if nonempty_array(case, "face_conditions") {
        out.push(BODY_FACE_CONDITIONS);
    }
    if case.get("cut_cells").and_then(Value::as_bool) == Some(true) {
        out.push(CUT_CELLS);
    }
    let heat = case.get("thermal").filter(|t| !t.is_null());
    let lumped = case
        .get("bodies")
        .and_then(Value::as_array)
        .is_some_and(|bodies| bodies.iter().any(|b| b.get("lumped").and_then(Value::as_bool) == Some(true)));
    if lumped || heat.and_then(|t| t.get("partial_cells")).and_then(Value::as_bool) == Some(true) {
        out.push(PARTIAL_CELLS);
    }
    if heat.is_some_and(|t| nonempty_array(t, "face_sources")) || (heat.is_some() && nonempty_array(case, "face_conditions")) {
        out.push(FACE_HEAT_SOURCES);
    }
    out
}

/// The blowing directions of internal fans.
pub const FAN_DIRECTIONS: [&str; 6] = ["+x", "-x", "+y", "-y", "+z", "-z"];
const MAX_INTERNAL_FANS: usize = 16;

/// Internal fans: named fluid cell blocks inside the grid, apart from solids, porous blocks
/// and each other, blowing along a grid direction with a falling pressure-flow curve.
fn validate_internal_fans(case: &Value, rows: &Value) -> Result<()> {
    flow_solver(case, "internal fans")?;
    let rows = rows.as_array().ok_or_else(|| Error::Invalid("internal_fans must be an array".into()))?;
    if rows.len() > MAX_INTERNAL_FANS {
        return Err(Error::Invalid(format!("at most {MAX_INTERNAL_FANS} internal fans")));
    }
    let dimensions: [usize; 3] = serde_json::from_value(case.pointer("/mesh/dimensions").cloned().unwrap_or(Value::Null))
        .map_err(|_| Error::Invalid("internal fans need three mesh dimensions".into()))?;
    let block = |row: &Value| -> Option<([usize; 3], [usize; 3])> {
        Some((serde_json::from_value(row.get("min")?.clone()).ok()?, serde_json::from_value(row.get("max")?.clone()).ok()?))
    };
    let overlaps = |a: ([usize; 3], [usize; 3]), b: ([usize; 3], [usize; 3])| (0..3).all(|i| a.0[i] < b.1[i] && b.0[i] < a.1[i]);
    let others: Vec<([usize; 3], [usize; 3])> = ["obstacles", "porous_blocks"]
        .iter()
        .flat_map(|key| case.get(*key).and_then(Value::as_array).into_iter().flatten())
        .filter_map(block)
        .collect();
    let mut names = BTreeSet::new();
    let mut fans: Vec<([usize; 3], [usize; 3])> = vec![];
    for row in rows {
        let object = row.as_object().ok_or_else(|| Error::Invalid("an internal fan must be an object".into()))?;
        if object.keys().any(|k| !matches!(k.as_str(), "name" | "min" | "max" | "direction" | "curve" | "relaxation")) {
            return Err(Error::Invalid("internal fans take name, min, max, direction, curve and relaxation only".into()));
        }
        let name = row["name"].as_str().unwrap_or("");
        let bytes = name.as_bytes();
        if bytes.is_empty() || bytes.len() > 60 || !bytes.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_') || !names.insert(name) {
            return Err(Error::Invalid("internal fan names must be unique identifiers of letters, digits and underscores".into()));
        }
        let Some(cells) = block(row) else {
            return Err(Error::Invalid(format!("internal fan {name} needs min and max cell indices")));
        };
        if (0..3).any(|a| cells.0[a] >= cells.1[a] || cells.1[a] > dimensions[a]) {
            return Err(Error::Invalid(format!("internal fan {name} must be a non-empty block inside the grid")));
        }
        if others.iter().any(|b| overlaps(cells, *b)) {
            return Err(Error::Invalid(format!("internal fan {name} overlaps a solid or porous block; its cells must be open fluid")));
        }
        if fans.iter().any(|b| overlaps(cells, *b)) {
            return Err(Error::Invalid(format!("internal fan {name} overlaps another fan")));
        }
        fans.push(cells);
        if !row["direction"].as_str().is_some_and(|d| FAN_DIRECTIONS.contains(&d)) {
            return Err(Error::Invalid(format!("internal fan {name} needs a direction: +x, -x, +y, -y, +z or -z")));
        }
        let curve: Vec<(f64, f64)> = row["curve"]
            .as_array()
            .into_iter()
            .flatten()
            .map(|p| (p["volume_flow_m3_per_s"].as_f64().unwrap_or(f64::NAN), p["pressure_rise_pa"].as_f64().unwrap_or(f64::NAN)))
            .collect();
        let points_ok = row["curve"].as_array().is_some_and(|list| {
            list.iter()
                .all(|p| p.as_object().is_some_and(|o| o.len() == 2 && o.contains_key("volume_flow_m3_per_s") && o.contains_key("pressure_rise_pa")))
        });
        if !points_ok
            || !(2..=MAX_CURVE_POINTS).contains(&curve.len())
            || !curve.iter().all(|(q, p)| q.is_finite() && p.is_finite() && *q >= 0.0 && *p >= 0.0)
            || !curve.windows(2).all(|w| w[1].0 > w[0].0 && w[1].1 < w[0].1)
        {
            return Err(Error::Invalid(format!(
                "internal fan {name} needs a curve of 2..{MAX_CURVE_POINTS} points with rising flow and falling pressure rise"
            )));
        }
        if let Some(r) = row.get("relaxation")
            && !r.as_f64().is_some_and(|r| r.is_finite() && r > 0.0 && r <= 1.0)
        {
            return Err(Error::Invalid(format!("internal fan {name} relaxation must be in (0, 1]")));
        }
    }
    Ok(())
}

/// The largest |Ω|·Δt the companion accepts for a rotating region.
pub const MAX_ROTATION_PER_STEP: f64 = 0.2;
const MAX_ROTATING_REGIONS: usize = 16;

/// Rotating regions: named cell blocks inside the grid, apart from each other, with a
/// finite non-zero axis and a speed the time step can follow (|Ω|·Δt ≤ 0.2).
fn validate_rotating(case: &Value, rows: &Value) -> Result<()> {
    flow_solver(case, "rotating regions")?;
    let rows = rows.as_array().ok_or_else(|| Error::Invalid("rotating_regions must be an array".into()))?;
    if rows.len() > MAX_ROTATING_REGIONS {
        return Err(Error::Invalid(format!("at most {MAX_ROTATING_REGIONS} rotating regions")));
    }
    let dimensions: [usize; 3] = serde_json::from_value(case.pointer("/mesh/dimensions").cloned().unwrap_or(Value::Null))
        .map_err(|_| Error::Invalid("rotating regions need three mesh dimensions".into()))?;
    let dt = case.pointer("/config/dt").and_then(Value::as_f64).unwrap_or(1e-2);
    let vector = |v: &Value| -> Option<[f64; 3]> {
        let out = ["x", "y", "z"].map(|a| v.get(a).and_then(Value::as_f64).unwrap_or(f64::NAN));
        (v.as_object().is_some_and(|o| o.len() == 3) && out.iter().all(|c| c.is_finite())).then_some(out)
    };
    let mut names = BTreeSet::new();
    let mut boxes: Vec<([usize; 3], [usize; 3])> = vec![];
    for row in rows {
        let object = row.as_object().ok_or_else(|| Error::Invalid("a rotating region must be an object".into()))?;
        if object.keys().any(|k| !matches!(k.as_str(), "name" | "min" | "max" | "axis_point" | "axis" | "rpm")) {
            return Err(Error::Invalid("rotating regions take name, min, max, axis_point, axis and rpm only".into()));
        }
        let name = row["name"].as_str().unwrap_or("");
        let bytes = name.as_bytes();
        if bytes.is_empty() || bytes.len() > 60 || !bytes.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_') || !names.insert(name) {
            return Err(Error::Invalid("rotating region names must be unique identifiers of letters, digits and underscores".into()));
        }
        let min: Option<[usize; 3]> = serde_json::from_value(row["min"].clone()).ok();
        let max: Option<[usize; 3]> = serde_json::from_value(row["max"].clone()).ok();
        let (Some(min), Some(max)) = (min, max) else {
            return Err(Error::Invalid(format!("rotating region {name} needs min and max cell indices")));
        };
        if (0..3).any(|a| min[a] >= max[a] || max[a] > dimensions[a]) {
            return Err(Error::Invalid(format!("rotating region {name} must be a non-empty block inside the grid")));
        }
        if boxes.iter().any(|(lo, hi)| (0..3).all(|a| min[a] < hi[a] && lo[a] < max[a])) {
            return Err(Error::Invalid(format!("rotating region {name} overlaps another region")));
        }
        boxes.push((min, max));
        let axis = vector(&row["axis"]).filter(|a| a.iter().map(|c| c * c).sum::<f64>() > 0.0);
        if axis.is_none() || vector(&row["axis_point"]).is_none() {
            return Err(Error::Invalid(format!("rotating region {name} needs a finite axis point and a non-zero axis")));
        }
        let rpm =
            row["rpm"].as_f64().filter(|v| v.is_finite()).ok_or_else(|| Error::Invalid(format!("rotating region {name} needs a finite rpm")))?;
        let omega = rpm.abs() * std::f64::consts::TAU / 60.0;
        if omega * dt > MAX_ROTATION_PER_STEP {
            return Err(Error::Invalid(format!(
                "rotating region {name} turns {:.3} rad per time step; at most {MAX_ROTATION_PER_STEP} (shorten the time step or lower the speed)",
                omega * dt
            )));
        }
    }
    Ok(())
}

/// Transferred profiles: an exterior patch each (once, not the fan's), one value per
/// background face of the patch, on a fixed-value velocity or pressure condition.
fn validate_transferred(case: &Value, rows: &Value) -> Result<()> {
    flow_solver(case, "transferred boundary conditions")?;
    let rows = rows.as_array().ok_or_else(|| Error::Invalid("transferred must be an array".into()))?;
    let dimensions: [usize; 3] = serde_json::from_value(case.pointer("/mesh/dimensions").cloned().unwrap_or(Value::Null))
        .map_err(|_| Error::Invalid("transferred profiles need three mesh dimensions".into()))?;
    let fixed = |key: &str, patch: &str| {
        case.get(key)
            .and_then(Value::as_array)
            .is_some_and(|list| list.iter().any(|row| row["patch"].as_str() == Some(patch) && row["kind"].as_str() == Some("fixed_value")))
    };
    let mut seen = std::collections::BTreeSet::new();
    for row in rows {
        let object = row.as_object().ok_or_else(|| Error::Invalid("a transferred entry must be an object".into()))?;
        if object.keys().any(|key| !matches!(key.as_str(), "patch" | "velocity" | "pressure")) {
            return Err(Error::Invalid("transferred entries take patch, velocity and pressure only".into()));
        }
        let patch = row["patch"].as_str().unwrap_or("");
        let in_plane = match patch {
            "left" | "right" => [1, 2],
            "bottom" | "top" => [0, 2],
            "front" | "back" => [0, 1],
            _ => return Err(Error::Invalid(format!("transferred patch {patch:?} is not an exterior patch"))),
        };
        if !seen.insert(patch) {
            return Err(Error::Invalid(format!("patch {patch} is transferred twice")));
        }
        if case.pointer("/fan/patch").and_then(Value::as_str) == Some(patch) {
            return Err(Error::Invalid(format!("the fan owns patch {patch}")));
        }
        let expected = dimensions[in_plane[0]] * dimensions[in_plane[1]];
        let (velocity, pressure) = (row.get("velocity").filter(|v| !v.is_null()), row.get("pressure").filter(|v| !v.is_null()));
        if velocity.is_none() && pressure.is_none() {
            return Err(Error::Invalid(format!("transferred patch {patch} has no values")));
        }
        if let Some(values) = velocity {
            let list = values
                .as_array()
                .filter(|list| list.len() == expected)
                .ok_or_else(|| Error::Invalid(format!("transferred velocity on {patch} needs {expected} values")))?;
            let finite = list.iter().all(|v| ["x", "y", "z"].iter().all(|a| v.get(*a).and_then(Value::as_f64).is_some_and(f64::is_finite)));
            if !finite || !fixed("velocity_boundaries", patch) {
                return Err(Error::Invalid(format!(
                    "transferred velocity on {patch} needs finite x, y, z values and a fixed-value velocity condition"
                )));
            }
        }
        if let Some(values) = pressure {
            let list = values
                .as_array()
                .filter(|list| list.len() == expected)
                .ok_or_else(|| Error::Invalid(format!("transferred pressure on {patch} needs {expected} values")))?;
            if !list.iter().all(|v| v.as_f64().is_some_and(f64::is_finite)) || !fixed("pressure_boundaries", patch) {
                return Err(Error::Invalid(format!("transferred pressure on {patch} needs finite values and a fixed-value pressure condition")));
            }
        }
    }
    Ok(())
}

/// Most nodes along one axis of explicit mesh coordinates (as the producer allows).
pub const MAX_AXIS_NODES: usize = 4097;

/// Explicit mesh coordinates: three lists from 0 to each length, strictly increasing, one more
/// entry than the axis's cells, on an unsheared `collocated_flow3d` grid.
fn validate_coordinates(case: &Value, coordinates: &Value) -> Result<()> {
    flow_solver(case, "mesh coordinates")?;
    let object = coordinates.as_object().ok_or_else(|| Error::Invalid("mesh coordinates must be an object".into()))?;
    if object.keys().any(|key| !matches!(key.as_str(), "x" | "y" | "z")) || object.len() != 3 {
        return Err(Error::Invalid("mesh coordinates need exactly x, y and z".into()));
    }
    if case.pointer("/mesh/xy_shear").is_some_and(|v| v.as_f64() != Some(0.0)) {
        return Err(Error::Invalid("mesh coordinates require zero xy_shear".into()));
    }
    let dimensions: [usize; 3] = serde_json::from_value(case.pointer("/mesh/dimensions").cloned().unwrap_or(Value::Null))
        .map_err(|_| Error::Invalid("mesh coordinates need three mesh dimensions".into()))?;
    for (axis, name) in ["x", "y", "z"].into_iter().enumerate() {
        let nodes: Vec<f64> = serde_json::from_value(object.get(name).cloned().unwrap_or(Value::Null))
            .map_err(|_| Error::Invalid(format!("mesh coordinates {name} must be a list of numbers")))?;
        let length = case.pointer(&format!("/mesh/lengths/{name}")).and_then(Value::as_f64).unwrap_or(f64::NAN);
        let last = nodes.last().copied().unwrap_or(f64::NAN);
        if nodes.len() > MAX_AXIS_NODES
            || nodes.len() != dimensions[axis] + 1
            || nodes.first() != Some(&0.0)
            || !nodes.iter().all(|v| v.is_finite())
            || !nodes.windows(2).all(|w| w[1] > w[0])
            || !(last - length).abs().le(&(1e-9 * length.abs()))
        {
            return Err(Error::Invalid(format!(
                "mesh coordinates {name} must run strictly upward from 0 to the length with one entry more than the cells"
            )));
        }
    }
    Ok(())
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CurvePoint {
    volume_flow_m3_per_s: f64,
    pressure_rise_pa: f64,
}

#[derive(Deserialize)]
#[serde(rename_all = "snake_case")]
enum Convention {
    Static,
    Total,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Fan {
    patch: String,
    reference_patch: String,
    pressure_convention: Convention,
    curve: Vec<CurvePoint>,
    #[serde(default)]
    initial_volume_flow_m3_per_s: Option<f64>,
    #[serde(default)]
    relaxation: Option<f64>,
    #[serde(default)]
    update_interval: Option<u64>,
    #[serde(default)]
    settle_tolerance: Option<f64>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PorousBlock {
    name: String,
    min: [u64; 3],
    max: [u64; 3],
    viscous_resistance_per_m2: f64,
    inertial_resistance_per_m: f64,
}

#[derive(Deserialize)]
struct Bounds {
    name: String,
    min: [u64; 3],
    max: [u64; 3],
}

fn flow_solver(case: &Value, feature: &str) -> Result<()> {
    if case.get("solver").and_then(Value::as_str) != Some("collocated_flow3d") {
        return Err(Error::Invalid(format!("{feature} require collocated_flow3d")));
    }
    Ok(())
}

fn patch_entries<'a>(case: &'a Value, key: &str) -> impl Iterator<Item = (&'a str, &'a str)> {
    case.get(key)
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|entry| Some((entry.get("patch").and_then(Value::as_str)?, entry.get("kind").and_then(Value::as_str)?)))
}

fn validate_fan(case: &Value, value: &Value) -> Result<()> {
    flow_solver(case, "fan curves")?;
    let fan: Fan = serde_json::from_value(value.clone()).map_err(|error| Error::Invalid(format!("invalid fan: {error}")))?;
    if !EXTERIOR.contains(&fan.patch.as_str()) || !EXTERIOR.contains(&fan.reference_patch.as_str()) {
        return Err(Error::Invalid("fan and reference patches must be exterior box patches".into()));
    }
    if fan.patch == fan.reference_patch {
        return Err(Error::Invalid("fan and reference patches must differ".into()));
    }
    if !(2..=MAX_CURVE_POINTS).contains(&fan.curve.len()) {
        return Err(Error::Invalid(format!("a fan curve needs 2..{MAX_CURVE_POINTS} points")));
    }
    if fan.curve.iter().any(|point| {
        !point.volume_flow_m3_per_s.is_finite()
            || !point.pressure_rise_pa.is_finite()
            || point.volume_flow_m3_per_s < 0.0
            || point.pressure_rise_pa < 0.0
    }) {
        return Err(Error::Invalid("fan curve points must be finite and non-negative".into()));
    }
    if fan.curve.windows(2).any(|pair| match pair {
        [a, b] => b.volume_flow_m3_per_s <= a.volume_flow_m3_per_s || b.pressure_rise_pa >= a.pressure_rise_pa,
        _ => true,
    }) {
        return Err(Error::Invalid("fan curve flow must strictly increase while pressure rise strictly decreases".into()));
    }
    if let Some(flow) = fan.initial_volume_flow_m3_per_s {
        let low = fan.curve.first().map_or(0.0, |point| point.volume_flow_m3_per_s);
        let high = fan.curve.last().map_or(0.0, |point| point.volume_flow_m3_per_s);
        if !flow.is_finite() || flow < low || flow > high {
            return Err(Error::Invalid("initial fan flow must lie within the curve's flow range".into()));
        }
    }
    if fan.relaxation.is_some_and(|value| !(value.is_finite() && value > 0.0 && value <= 1.0)) {
        return Err(Error::Invalid("fan relaxation must be in (0, 1]".into()));
    }
    if fan.update_interval.is_some_and(|value| !(1..=10_000).contains(&value)) {
        return Err(Error::Invalid("fan update interval must be 1..10000 steps".into()));
    }
    if fan.settle_tolerance.is_some_and(|value| !(value.is_finite() && value > 0.0 && value <= 1.0)) {
        return Err(Error::Invalid("fan settle tolerance must be in (0, 1]".into()));
    }
    // One owner per boundary: an authored velocity on the fan patch would silently compete with it.
    if patch_entries(case, "velocity_boundaries").any(|(patch, _)| patch == fan.patch) {
        return Err(Error::Invalid("the fan owns its patch velocity; remove that patch from velocity_boundaries".into()));
    }
    if patch_entries(case, "pressure_boundaries").any(|(patch, kind)| patch == fan.patch && kind == "fixed_value") {
        return Err(Error::Invalid("a fan patch cannot also have a fixed pressure".into()));
    }
    if !patch_entries(case, "pressure_boundaries").any(|(patch, kind)| patch == fan.reference_patch && kind == "fixed_value") {
        return Err(Error::Invalid("the fan reference patch needs a fixed-value pressure boundary".into()));
    }
    // The convention must be explicit; the companion adds fan-face dynamic pressure for Total.
    match fan.pressure_convention {
        Convention::Static | Convention::Total => Ok(()),
    }
}

fn validate_porous(case: &Value, value: &Value) -> Result<()> {
    let array = value.as_array().ok_or_else(|| Error::Invalid("porous_blocks must be an array".into()))?;
    if array.is_empty() {
        return Ok(());
    }
    flow_solver(case, "porous blocks")?;
    if array.len() > MAX_POROUS_BLOCKS {
        return Err(Error::Limit(format!("porous blocks exceed {MAX_POROUS_BLOCKS}")));
    }
    let mesh = case.get("mesh").ok_or_else(|| Error::Invalid("porous blocks require an explicit mesh".into()))?;
    let dimensions: [u64; 3] = serde_json::from_value(mesh.get("dimensions").cloned().unwrap_or(Value::Null))
        .map_err(|_| Error::Invalid("porous mesh dimensions must be three positive integers".into()))?;
    if mesh.get("xy_shear").is_some_and(|value| value.as_f64() != Some(0.0)) {
        return Err(Error::Invalid("porous blocks require zero xy_shear".into()));
    }
    let solids: Vec<Bounds> = match case.get("obstacles") {
        Some(obstacles) => serde_json::from_value(obstacles.clone()).map_err(|_| Error::Invalid("invalid grid obstacles".into()))?,
        None => vec![],
    };
    let mut names: BTreeSet<String> = EXTERIOR.into_iter().map(String::from).collect();
    names.extend(solids.iter().map(|solid| solid.name.clone()));
    let overlaps = |a: ([u64; 3], [u64; 3]), b: ([u64; 3], [u64; 3])| (0..3).all(|axis| a.0[axis] < b.1[axis] && b.0[axis] < a.1[axis]);
    let mut accepted: Vec<([u64; 3], [u64; 3])> = vec![];
    for value in array {
        let block: PorousBlock = serde_json::from_value(value.clone()).map_err(|error| Error::Invalid(format!("invalid porous block: {error}")))?;
        let name = block.name.as_bytes();
        if name.is_empty()
            || name.len() > 80
            || !name.first().is_some_and(u8::is_ascii_alphanumeric)
            || !name.iter().all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-'))
            || !names.insert(block.name.clone())
        {
            return Err(Error::Invalid("porous block names must be unique ASCII identifiers distinct from patches and obstacles".into()));
        }
        if block.min.iter().zip(block.max).zip(dimensions).any(|((min, max), size)| *min >= max || max > size) {
            return Err(Error::Invalid("porous block bounds must satisfy 0 <= min < max <= dimensions".into()));
        }
        let (viscous, inertial) = (block.viscous_resistance_per_m2, block.inertial_resistance_per_m);
        if !(viscous.is_finite() && inertial.is_finite() && viscous >= 0.0 && inertial >= 0.0 && (viscous > 0.0 || inertial > 0.0)) {
            return Err(Error::Invalid("porous coefficients must be finite, non-negative and not both zero".into()));
        }
        let bounds = (block.min, block.max);
        if solids.iter().any(|solid| overlaps(bounds, (solid.min, solid.max))) {
            return Err(Error::Invalid("porous blocks must not overlap solid obstacles".into()));
        }
        if accepted.iter().any(|other| overlaps(bounds, *other)) {
            return Err(Error::Invalid("porous blocks must not overlap each other".into()));
        }
        accepted.push(bounds);
    }
    Ok(())
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ThermalFluid {
    specific_heat_j_per_kg_k: f64,
    conductivity_w_per_m_k: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ThermalSolid {
    obstacle: String,
    conductivity_w_per_m_k: [f64; 3],
    power_w: f64,
    #[serde(default)]
    max_temperature_k: Option<f64>,
    #[serde(default)]
    density_kg_per_m3: Option<f64>,
    #[serde(default)]
    specific_heat_j_per_kg_k: Option<f64>,
    #[serde(default)]
    contact_resistance_m2_k_per_w: Option<f64>,
    #[serde(default)]
    junction_to_case_k_per_w: Option<f64>,
    #[serde(default)]
    heat_pipe: Option<HeatPipeInput>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct HeatPipeInput {
    axis: u64,
    effective_resistance_k_per_w: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CoolerInput {
    name: String,
    cold_solid: String,
    hot_solid: String,
    heat_pumped_w: f64,
    electrical_power_w: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BuoyancyInput {
    gravity_m_per_s2: [f64; 3],
    expansion_per_k: f64,
    reference_temperature_k: f64,
}

#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum ThermalCondition {
    Adiabatic,
    Temperature { temperature_k: f64 },
    HeatFlux { w_per_m2: f64 },
    Convection { h_w_per_m2_k: f64, ambient_k: f64 },
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ThermalPatch {
    patch: String,
    condition: ThermalCondition,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Thermal {
    fluid: ThermalFluid,
    inflow_temperature_k: f64,
    solids: Vec<ThermalSolid>,
    #[serde(default)]
    boundaries: Vec<ThermalPatch>,
    #[serde(default)]
    buoyancy: Option<BuoyancyInput>,
    #[serde(default)]
    thermoelectric_coolers: Vec<CoolerInput>,
    #[serde(default)]
    face_sources: Vec<FaceSourceInput>,
    #[serde(default)]
    #[allow(dead_code)]
    partial_cells: bool,
}

/// Power on part of a surface body (`face_heat_sources_v1`).
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FaceSourceInput {
    name: String,
    body: String,
    triangles: Vec<usize>,
    power_w: f64,
}

/// Most face heat sources a case holds.
const MAX_FACE_SOURCES: usize = 64;

fn kelvin(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn validate_thermal(case: &Value, value: &Value) -> Result<()> {
    flow_solver(case, "conjugate heat transfer")?;
    let thermal: Thermal = serde_json::from_value(value.clone()).map_err(|error| Error::Invalid(format!("invalid thermal block: {error}")))?;
    if !(kelvin(thermal.fluid.specific_heat_j_per_kg_k) && kelvin(thermal.fluid.conductivity_w_per_m_k)) {
        return Err(Error::Invalid("fluid specific heat and conductivity must be finite and positive".into()));
    }
    if !kelvin(thermal.inflow_temperature_k) {
        return Err(Error::Invalid("inflow temperature is absolute, in kelvin, and must be positive".into()));
    }
    if case.get("mesh").and_then(|mesh| mesh.get("xy_shear")).is_some_and(|value| value.as_f64() != Some(0.0)) {
        return Err(Error::Invalid("conjugate heat transfer requires zero xy_shear".into()));
    }
    let obstacles: Vec<Bounds> = match case.get("obstacles") {
        Some(obstacles) => serde_json::from_value(obstacles.clone()).map_err(|_| Error::Invalid("invalid grid obstacles".into()))?,
        None => vec![],
    };
    let mut solid_names: Vec<String> = obstacles.iter().map(|o| o.name.clone()).collect();
    solid_names.extend(body_names(case));
    if thermal.solids.len() > 128 {
        return Err(Error::Limit("thermal solids exceed 128".into()));
    }
    let mut named = BTreeSet::new();
    for solid in &thermal.solids {
        if !solid_names.contains(&solid.obstacle) {
            return Err(Error::Invalid(format!("thermal solid {} is not an obstacle or surface body", solid.obstacle)));
        }
        if !named.insert(solid.obstacle.as_str()) {
            return Err(Error::Invalid("each obstacle needs exactly one solid entry".into()));
        }
        if !solid.conductivity_w_per_m_k.iter().all(|k| kelvin(*k)) || !solid.power_w.is_finite() || solid.power_w < 0.0 {
            return Err(Error::Invalid("solid conductivities must be positive and power non-negative".into()));
        }
        if solid.max_temperature_k.is_some_and(|limit| !kelvin(limit)) {
            return Err(Error::Invalid("temperature limits are absolute, in kelvin".into()));
        }
        if solid.density_kg_per_m3.is_some_and(|v| !kelvin(v)) || solid.specific_heat_j_per_kg_k.is_some_and(|v| !kelvin(v)) {
            return Err(Error::Invalid("solid density and specific heat must be positive".into()));
        }
        let resistance = |r: Option<f64>| r.is_none_or(|r| r.is_finite() && r >= 0.0);
        if !resistance(solid.contact_resistance_m2_k_per_w) || !resistance(solid.junction_to_case_k_per_w) {
            return Err(Error::Invalid(format!("{}: contact and junction resistances must be finite and non-negative", solid.obstacle)));
        }
        if solid.heat_pipe.as_ref().is_some_and(|pipe| pipe.axis > 2 || !kelvin(pipe.effective_resistance_k_per_w)) {
            return Err(Error::Invalid(format!("{}: a heat pipe needs axis 0, 1 or 2 and a positive resistance", solid.obstacle)));
        }
    }
    if thermal.thermoelectric_coolers.len() > 32 {
        return Err(Error::Limit("at most 32 thermoelectric coolers".into()));
    }
    let mut coolers = BTreeSet::new();
    for tec in &thermal.thermoelectric_coolers {
        if tec.name.is_empty() || tec.name.len() > 60 || !coolers.insert(tec.name.as_str()) {
            return Err(Error::Invalid("thermoelectric coolers need unique names of at most 60 characters".into()));
        }
        if tec.cold_solid == tec.hot_solid || !named.contains(tec.cold_solid.as_str()) || !named.contains(tec.hot_solid.as_str()) {
            return Err(Error::Invalid(format!("{}: cold and hot sides must be two different thermal solids", tec.name)));
        }
        let watts = |w: f64| w.is_finite() && w >= 0.0;
        if !watts(tec.heat_pumped_w) || !watts(tec.electrical_power_w) {
            return Err(Error::Invalid(format!("{}: pumped heat and electrical power must be finite and non-negative", tec.name)));
        }
    }
    // No solid material is implied: every obstacle must be described.
    if named.len() != solid_names.len() {
        return Err(Error::Invalid("every obstacle and surface body needs a thermal solid entry".into()));
    }
    if let Some(b) = &thermal.buoyancy
        && !(b.gravity_m_per_s2.iter().all(|g| g.is_finite()) && kelvin(b.expansion_per_k) && kelvin(b.reference_temperature_k))
    {
        return Err(Error::Invalid("buoyancy needs finite gravity, positive expansion and a kelvin reference".into()));
    }
    if thermal.face_sources.len() > MAX_FACE_SOURCES {
        return Err(Error::Limit(format!("at most {MAX_FACE_SOURCES} face heat sources")));
    }
    let mut sources = BTreeSet::new();
    for source in &thermal.face_sources {
        if source.name.is_empty() || source.name.len() > 80 || !sources.insert(source.name.as_str()) {
            return Err(Error::Invalid("face heat sources need unique names of at most 80 characters".into()));
        }
        let triangles = case
            .get("bodies")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .find(|b| b["name"].as_str() == Some(source.body.as_str()))
            .and_then(|b| b["triangles"].as_array())
            .map(Vec::len)
            .ok_or_else(|| Error::Invalid(format!("face heat source {} names no surface body {}", source.name, source.body)))?;
        if source.triangles.is_empty() || source.triangles.iter().any(|t| *t >= triangles) {
            return Err(Error::Invalid(format!("face heat source {}: list triangles of {}", source.name, source.body)));
        }
        if !source.power_w.is_finite() {
            return Err(Error::Invalid(format!("face heat source {}: power must be finite", source.name)));
        }
    }
    // The heat solve runs on the grid, where a body face lies between fluid and solid: flow
    // through a body face would carry heat out of the solid, so only closed ones are taken.
    for condition in case.get("face_conditions").and_then(Value::as_array).into_iter().flatten() {
        let Some(name) = condition["name"].as_str() else { continue };
        let rows = |key: &str| case.get(key).and_then(Value::as_array).into_iter().flatten().filter(move |row| row["patch"].as_str() == Some(name));
        let moving = |row: &Value| ["x", "y", "z"].iter().any(|k| row["value"][*k].as_f64().is_some_and(|v| v != 0.0));
        let open = rows("pressure_boundaries").any(|row| row["kind"] == "fixed_value")
            || rows("velocity_boundaries").any(|row| row["kind"] == "zero_gradient" || (row["kind"] == "fixed_value" && moving(row)));
        if open {
            return Err(Error::Invalid(format!(
                "{name}: flow through a body face does not work with heat transfer yet; put openings on faces on the domain boundary"
            )));
        }
    }
    let mut patches = BTreeSet::new();
    for boundary in &thermal.boundaries {
        if !EXTERIOR.contains(&boundary.patch.as_str()) || !patches.insert(boundary.patch.as_str()) {
            return Err(Error::Invalid("thermal boundaries name each exterior patch at most once".into()));
        }
        let ok = match boundary.condition {
            ThermalCondition::Adiabatic => true,
            ThermalCondition::Temperature { temperature_k } => kelvin(temperature_k),
            ThermalCondition::HeatFlux { w_per_m2 } => w_per_m2.is_finite(),
            ThermalCondition::Convection { h_w_per_m2_k, ambient_k } => kelvin(h_w_per_m2_k) && kelvin(ambient_k),
        };
        if !ok {
            return Err(Error::Invalid(format!("invalid thermal condition on {}", boundary.patch)));
        }
    }
    Ok(())
}

fn body_names(case: &Value) -> Vec<String> {
    case.get("bodies")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|b| b.get("name").and_then(Value::as_str).map(String::from))
        .collect()
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SurfaceBody {
    name: String,
    vertices: Vec<[f64; 3]>,
    triangles: Vec<[usize; 3]>,
    /// May enclose no cell centre (`partial_cells_v1`).
    #[serde(default)]
    #[allow(dead_code)]
    lumped: bool,
}

/// Closed, consistently oriented, outward-wound triangle surfaces, as the companion requires.
fn validate_bodies(case: &Value, value: &Value) -> Result<()> {
    let array = value.as_array().ok_or_else(|| Error::Invalid("bodies must be an array".into()))?;
    if array.is_empty() {
        return Ok(());
    }
    flow_solver(case, "surface bodies")?;
    let obstacles = case.get("obstacles").and_then(Value::as_array).map_or(0, Vec::len);
    if array.len() + obstacles > 128 {
        return Err(Error::Limit("obstacles and surface bodies exceed 128".into()));
    }
    if case.get("mesh").and_then(|mesh| mesh.get("xy_shear")).is_some_and(|value| value.as_f64() != Some(0.0)) {
        return Err(Error::Invalid("surface bodies require zero xy_shear".into()));
    }
    let mut names: BTreeSet<String> = EXTERIOR.into_iter().map(String::from).collect();
    for o in case.get("obstacles").and_then(Value::as_array).into_iter().flatten() {
        if let Some(name) = o.get("name").and_then(Value::as_str) {
            names.insert(name.into());
        }
    }
    for value in array {
        let body: SurfaceBody = serde_json::from_value(value.clone()).map_err(|error| Error::Invalid(format!("invalid surface body: {error}")))?;
        let name = body.name.as_bytes();
        if name.is_empty()
            || name.len() > 80
            || !name.first().is_some_and(u8::is_ascii_alphanumeric)
            || !name.iter().all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-'))
            || !names.insert(body.name.clone())
        {
            return Err(Error::Invalid("surface body names must be unique identifiers distinct from patches and obstacles".into()));
        }
        if body.triangles.is_empty() || body.triangles.len() > 200_000 || !body.vertices.iter().flatten().all(|v| v.is_finite() && v.abs() < 1e6) {
            return Err(Error::Invalid(format!("surface body {} needs 1..200000 triangles and finite metre coordinates", body.name)));
        }
        let mut edges = std::collections::BTreeMap::new();
        let mut volume = 0.0;
        for t in &body.triangles {
            if t.iter().any(|i| *i >= body.vertices.len()) || t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
                return Err(Error::Invalid(format!("surface body {} has an invalid triangle", body.name)));
            }
            let [a, b, c] = t.map(|i| body.vertices[i]);
            volume += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
            for edge in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                *edges.entry(edge).or_insert(0u32) += 1;
            }
        }
        if !edges.iter().all(|((u, v), n)| *n == 1 && edges.get(&(*v, *u)) == Some(&1)) {
            return Err(Error::Invalid(format!("surface body {} must be closed, manifold and consistently oriented", body.name)));
        }
        if volume <= 0.0 {
            return Err(Error::Invalid(format!("surface body {} encloses no volume or is wound inward", body.name)));
        }
    }
    Ok(())
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GoalInput {
    name: String,
    kind: GoalKindInput,
    #[serde(default)]
    target: String,
    #[serde(default)]
    point: Option<[f64; 3]>,
    #[serde(default)]
    expression: Option<String>,
}

#[derive(Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
enum GoalKindInput {
    PatchVolumeFlow,
    PatchMeanPressure,
    PatchHeatOut,
    SolidMaxTemperature,
    SolidMeanTemperature,
    GlobalMaxFluidTemperature,
    GlobalMeanFluidTemperature,
    GlobalMaxSolidTemperature,
    PointTemperature,
    PointPressure,
    GlobalMaxSpeed,
    GlobalMeanSpeed,
    GlobalMaxPressure,
    GlobalMinPressure,
    GlobalMeanPressure,
    Equation,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GoalControlInput {
    window_steps: u64,
    relative_tolerance: f64,
    #[serde(default)]
    stop_when_converged: bool,
    #[serde(default)]
    min_steps: u64,
}

fn validate_goals(case: &Value) -> Result<()> {
    let goals: Vec<GoalInput> = match case.get("goals") {
        Some(value) => serde_json::from_value(value.clone()).map_err(|error| Error::Invalid(format!("invalid goals: {error}")))?,
        None => vec![],
    };
    if let Some(control) = case.get("goal_control").filter(|v| !v.is_null()) {
        let control: GoalControlInput =
            serde_json::from_value(control.clone()).map_err(|error| Error::Invalid(format!("invalid goal control: {error}")))?;
        if !(1..=10_000).contains(&control.window_steps)
            || !(control.relative_tolerance.is_finite() && control.relative_tolerance > 0.0 && control.relative_tolerance < 1.0)
            || control.min_steps > 20_000
        {
            return Err(Error::Invalid("goal window must be 1..10000 steps and tolerance in (0, 1)".into()));
        }
        // An unused flag is harmless, but keep it typed.
        let _ = control.stop_when_converged;
    }
    if goals.is_empty() {
        return Ok(());
    }
    flow_solver(case, "goals")?;
    if goals.len() > 32 {
        return Err(Error::Limit("at most 32 goals".into()));
    }
    let thermal = case.get("thermal").is_some_and(|t| !t.is_null());
    let mut solids: Vec<String> = case
        .get("obstacles")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|o| o.get("name").and_then(Value::as_str).map(String::from))
        .collect();
    solids.extend(body_names(case));
    // Patches on body faces (`body_face_conditions_v1`) are boundary patches like the exterior ones.
    let face_patches: Vec<String> = case
        .get("face_conditions")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|c| c.get("name").and_then(Value::as_str).map(String::from))
        .collect();
    let lengths: [Option<f64>; 3] = ["x", "y", "z"]
        .map(|axis| case.get("mesh").and_then(|mesh| mesh.get("lengths")).and_then(|lengths| lengths.get(axis)).and_then(Value::as_f64));
    let mut names = BTreeSet::new();
    for goal in &goals {
        let n = goal.name.as_bytes();
        if n.is_empty() || n.len() > 60 || !n.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_') || !names.insert(goal.name.as_str()) {
            return Err(Error::Invalid("goal names must be unique identifiers of letters, digits and underscores".into()));
        }
        let ok = match goal.kind {
            GoalKindInput::SolidMaxTemperature | GoalKindInput::SolidMeanTemperature => thermal && solids.contains(&goal.target),
            GoalKindInput::PatchHeatOut => thermal && EXTERIOR.contains(&goal.target.as_str()),
            GoalKindInput::PatchVolumeFlow | GoalKindInput::PatchMeanPressure => {
                EXTERIOR.contains(&goal.target.as_str()) || solids.contains(&goal.target) || face_patches.contains(&goal.target)
            }
            GoalKindInput::GlobalMaxFluidTemperature | GoalKindInput::GlobalMeanFluidTemperature | GoalKindInput::GlobalMaxSolidTemperature => {
                thermal
            }
            GoalKindInput::PointTemperature | GoalKindInput::PointPressure => {
                (goal.kind == GoalKindInput::PointPressure || thermal)
                    && goal.point.is_some_and(|p| p.iter().zip(lengths).all(|(v, l)| v.is_finite() && *v >= 0.0 && l.is_none_or(|l| *v <= l)))
            }
            GoalKindInput::GlobalMaxSpeed
            | GoalKindInput::GlobalMeanSpeed
            | GoalKindInput::GlobalMaxPressure
            | GoalKindInput::GlobalMinPressure
            | GoalKindInput::GlobalMeanPressure => true,
            GoalKindInput::Equation => equation_refers_to_goals(goal.expression.as_deref(), &goals),
        };
        if !ok {
            return Err(Error::Invalid(format!("goal {} has an unknown target or needs a thermal block", goal.name)));
        }
    }
    Ok(())
}

/// An equation goal's text uses only `+ - * / ( )`, numbers and the names of non-equation goals
/// (the solver parses it fully; this keeps unknown names out of the authored case).
fn equation_refers_to_goals(expression: Option<&str>, goals: &[GoalInput]) -> bool {
    let Some(text) = expression.filter(|t| !t.trim().is_empty() && t.len() <= 400) else { return false };
    if !text.chars().all(|c| c.is_ascii_alphanumeric() || "_+-*/(). ".contains(c)) {
        return false;
    }
    text.split(|c: char| !(c.is_ascii_alphanumeric() || c == '_' || c == '.'))
        .filter(|word| word.starts_with(|c: char| c.is_ascii_alphabetic() || c == '_'))
        .all(|word| goals.iter().any(|g| g.name == word && g.kind != GoalKindInput::Equation))
}

pub(crate) fn validate(case: &Value) -> Result<()> {
    if let Some(update) = case.pointer("/config/piso_neighbour_update")
        && !update.is_boolean()
    {
        return Err(Error::Invalid("config.piso_neighbour_update must be true or false".into()));
    }
    if let Some(scaled) = case.pointer("/config/scaled_momentum")
        && !scaled.is_boolean()
    {
        return Err(Error::Invalid("config.scaled_momentum must be true or false".into()));
    }
    if let Some(rows) = case.get("rotating_regions").filter(|rows| !rows.is_null()) {
        validate_rotating(case, rows)?;
    }
    if let Some(flag) = case.pointer("/config/balanced_body_force")
        && !flag.is_boolean()
    {
        return Err(Error::Invalid("config.balanced_body_force must be true or false".into()));
    }
    if let Some(rows) = case.get("internal_fans").filter(|rows| !rows.is_null()) {
        validate_internal_fans(case, rows)?;
    }
    if let Some(rows) = case.get("transferred").filter(|rows| !rows.is_null()) {
        validate_transferred(case, rows)?;
    }
    if let Some(coordinates) = case.pointer("/mesh/coordinates").filter(|c| !c.is_null()) {
        validate_coordinates(case, coordinates)?;
    }
    if let Some(bodies) = case.get("bodies") {
        validate_bodies(case, bodies)?;
    }
    validate_goals(case)?;
    if let Some(thermal) = case.get("thermal").filter(|thermal| !thermal.is_null()) {
        validate_thermal(case, thermal)?;
    }
    if let Some(fan) = case.get("fan").filter(|fan| !fan.is_null()) {
        validate_fan(case, fan)?;
    }
    if let Some(blocks) = case.get("porous_blocks") {
        validate_porous(case, blocks)?;
    }
    Ok(())
}
