//! Read-only preparation of original, reviewable primitive-flow feature drafts.
use crate::{Error, MAX_CASE_BYTES, Result, validate_case, validate_json};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

const PATCHES: [&str; 6] = ["left", "right", "bottom", "top", "front", "back"];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum FeatureKind {
    Fan,
    Porous,
    SolidMaterial,
    Goal,
    /// Contact resistance on a thermal solid (`name` is the solid).
    ContactResistance,
    /// Two-resistor junction on a thermal solid (`name` is the solid).
    TwoResistor,
    /// Heat pipe along x on a thermal solid (`name` is the solid).
    HeatPipe,
    /// Thermoelectric cooler between the first two thermal solids (`name` is the cooler).
    ThermoelectricCooler,
    /// Global goal: maximum fluid temperature (thermal cases).
    GlobalGoal,
    /// Point goal at the domain centre: temperature in thermal cases, pressure otherwise.
    PointGoal,
    /// Volume goal: maximum temperature of the first thermal solid.
    VolumeGoal,
    /// Equation goal over the first existing goal.
    EquationGoal,
}

impl FeatureKind {
    fn goal(self) -> bool {
        matches!(self, Self::Goal | Self::GlobalGoal | Self::PointGoal | Self::VolumeGoal | Self::EquationGoal)
    }
}

#[derive(Clone, Debug, Serialize, JsonSchema)]
pub struct PreparedFeature {
    pub case: Value,
    /// RFC 6901 pointer to the proposed record.
    pub pointer: String,
    pub issues: Vec<String>,
    pub defaults_notice: String,
}

fn bad(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

fn named(case: &Value, key: &str, field: &str, name: &str) -> bool {
    case.get(key).and_then(Value::as_array).is_some_and(|items| items.iter().any(|item| item[field].as_str() == Some(name)))
}

fn push(case: &mut Value, key: &str, item: Value, maximum: usize) -> Result<String> {
    let object = case.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
    let rows = object.entry(key).or_insert_with(|| json!([])).as_array_mut().ok_or_else(|| bad(format!("{key} must be an array")))?;
    if rows.len() >= maximum {
        return Err(bad(format!("{key} exceeds {maximum} records")));
    }
    let pointer = format!("/{key}/{}", rows.len());
    rows.push(item);
    Ok(pointer)
}

/// The pointer of the thermal-solid record named `solid`.
fn solid_record(case: &Value, solid: &str) -> Result<String> {
    case.pointer("/thermal/solids")
        .and_then(Value::as_array)
        .and_then(|rows| rows.iter().position(|row| row["obstacle"].as_str() == Some(solid)))
        .map(|index| format!("/thermal/solids/{index}"))
        .ok_or_else(|| bad("this device needs a thermal solid of that name; insert its solid material first"))
}

/// Set `field` on the thermal solid `solid` unless it is already set.
fn solid_field(proposed: &mut Value, solid: &str, field: &str, value: Value) -> Result<String> {
    let record = solid_record(proposed, solid)?;
    let row = proposed.pointer_mut(&record).and_then(Value::as_object_mut).ok_or_else(|| bad("invalid thermal solid record"))?;
    if row.get(field).is_some_and(|v| !v.is_null()) {
        return Err(bad(format!("{solid} already has {field}")));
    }
    row.insert(field.into(), value);
    Ok(format!("{record}/{field}"))
}

fn thermal(case: &Value) -> bool {
    case.get("thermal").is_some_and(|t| !t.is_null())
}

fn starter_cell(case: &Value) -> Result<[u64; 3]> {
    if case["bodies"].as_array().is_some_and(|bodies| !bodies.is_empty()) {
        return Err(bad("surface-body porous zones require explicitly reviewed bounds"));
    }
    let dimensions: [u64; 3] = serde_json::from_value(case.pointer("/mesh/dimensions").cloned().unwrap_or(Value::Null))
        .map_err(|_| bad("porous creation requires three explicit mesh dimensions"))?;
    let count = dimensions
        .iter()
        .try_fold(1u64, |a, b| a.checked_mul(*b))
        .filter(|n| *n > 0 && *n <= 32768)
        .ok_or_else(|| bad("porous creation requires 1..32768 grid cells"))?;
    let mut zones = Vec::new();
    for key in ["obstacles", "porous_blocks"] {
        if let Some(value) = case.get(key) {
            let rows = value.as_array().filter(|rows| rows.len() <= 128).ok_or_else(|| bad("invalid existing grid-zone collection"))?;
            for row in rows {
                let min: [u64; 3] = serde_json::from_value(row["min"].clone()).map_err(|_| bad("invalid existing zone minimum"))?;
                let max: [u64; 3] = serde_json::from_value(row["max"].clone()).map_err(|_| bad("invalid existing zone maximum"))?;
                if (0..3).any(|axis| min[axis] >= max[axis] || max[axis] > dimensions[axis]) {
                    return Err(bad("invalid existing zone bounds"));
                }
                zones.push((min, max));
            }
        }
    }
    for index in 0..count {
        let cell = [index % dimensions[0], (index / dimensions[0]) % dimensions[1], index / (dimensions[0] * dimensions[1])];
        if !zones.iter().any(|(min, max)| (0..3).all(|axis| cell[axis] >= min[axis] && cell[axis] < max[axis])) {
            return Ok(cell);
        }
    }
    Err(bad("no unoccupied grid cell is available for a porous starter"))
}

pub fn prepare_feature(case: &Value, kind: FeatureKind, name: &str) -> Result<PreparedFeature> {
    prepare_feature_with_options(case, kind, name, false)
}

pub fn prepare_feature_with_options(case: &Value, kind: FeatureKind, name: &str, replace_fan_patch_conditions: bool) -> Result<PreparedFeature> {
    validate_json(case, MAX_CASE_BYTES)?;
    if !case.is_object() || case["solver"].as_str() != Some("collocated_flow3d") {
        return Err(bad("feature creation requires collocated_flow3d"));
    }
    let name = name.trim();
    if name.is_empty()
        || name.len() > 60
        || !name.as_bytes().first().is_some_and(u8::is_ascii_alphanumeric)
        || !name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || (!kind.goal() && b == b'-'))
    {
        return Err(bad("feature name requires 1..60 ASCII identifier characters; goals allow letters, digits and underscores"));
    }
    if replace_fan_patch_conditions && kind != FeatureKind::Fan {
        return Err(bad("boundary replacement is a fan-only option"));
    }
    let mut proposed = case.clone();
    let pointer = match kind {
        FeatureKind::Fan => {
            if case.get("fan").is_some_and(|fan| !fan.is_null()) {
                return Err(bad("the case already contains its external fan"));
            }
            let patch = case["velocity_boundaries"]
                .as_array()
                .and_then(|rows| rows.iter().find_map(|row| row["patch"].as_str().filter(|patch| PATCHES.contains(patch))))
                .unwrap_or("left")
                .to_string();
            let reference = case["pressure_boundaries"]
                .as_array()
                .and_then(|rows| {
                    rows.iter().find_map(|row| {
                        row["patch"].as_str().filter(|p| *p != patch && PATCHES.contains(p) && row["kind"].as_str() == Some("fixed_value"))
                    })
                })
                .unwrap_or(if patch == "right" { "left" } else { "right" })
                .to_string();
            if replace_fan_patch_conditions {
                for key in ["velocity_boundaries", "pressure_boundaries"] {
                    if let Some(rows) = proposed.get_mut(key).and_then(Value::as_array_mut) {
                        rows.retain(|row| {
                            row["patch"].as_str() != Some(patch.as_str())
                                || (key == "pressure_boundaries" && row["kind"].as_str() != Some("fixed_value"))
                        });
                    }
                }
            }
            proposed["fan"] = json!({"patch":patch,"reference_patch":reference,"pressure_convention":"static",
                "curve":[{"volume_flow_m3_per_s":0.0,"pressure_rise_pa":10.0},{"volume_flow_m3_per_s":0.001,"pressure_rise_pa":0.0}],
                "relaxation":0.7,"update_interval":25,"settle_tolerance":0.002});
            "/fan".to_string()
        }
        FeatureKind::Porous => {
            if PATCHES.contains(&name) || ["obstacles", "bodies", "porous_blocks"].iter().any(|key| named(case, key, "name", name)) {
                return Err(bad("porous name conflicts with an existing region"));
            }
            let min = starter_cell(case)?;
            let max = min.map(|n| n + 1);
            push(
                &mut proposed,
                "porous_blocks",
                json!({"name":name,"min":min,"max":max,"viscous_resistance_per_m2":1000000.0,"inertial_resistance_per_m":0.0}),
                32,
            )?
        }
        FeatureKind::SolidMaterial => {
            if !["obstacles", "bodies"].iter().any(|key| named(case, key, "name", name)) {
                return Err(bad("solid material must reference an existing numerical body"));
            }
            if case
                .pointer("/thermal/solids")
                .and_then(Value::as_array)
                .is_some_and(|rows| rows.iter().any(|row| row["obstacle"].as_str() == Some(name)))
            {
                return Err(bad("the numerical body already has a solid thermal record"));
            }
            if proposed.get("thermal").is_none_or(Value::is_null) {
                proposed["thermal"] =
                    json!({"fluid":{"specific_heat_j_per_kg_k":1000.0,"conductivity_w_per_m_k":0.1},"inflow_temperature_k":293.15,"solids":[]});
            }
            let thermal = proposed.get_mut("thermal").ok_or_else(|| bad("thermal block missing"))?;
            let pointer = push(
                thermal,
                "solids",
                json!({"obstacle":name,"conductivity_w_per_m_k":[1.0,1.0,1.0],"power_w":0.0,"density_kg_per_m3":1000.0,"specific_heat_j_per_kg_k":1000.0}),
                128,
            )?;
            format!("/thermal{pointer}")
        }
        FeatureKind::Goal => {
            if named(case, "goals", "name", name) {
                return Err(bad("goal name already exists"));
            }
            let target = case["pressure_boundaries"]
                .as_array()
                .and_then(|rows| rows.iter().find_map(|row| row["patch"].as_str().filter(|p| PATCHES.contains(p))))
                .unwrap_or("right");
            push(&mut proposed, "goals", json!({"name":name,"kind":"patch_volume_flow","target":target}), 32)?
        }
        FeatureKind::ContactResistance => solid_field(&mut proposed, name, "contact_resistance_m2_k_per_w", json!(0.0001))?,
        FeatureKind::TwoResistor => solid_field(&mut proposed, name, "junction_to_case_k_per_w", json!(1.0))?,
        FeatureKind::HeatPipe => solid_field(&mut proposed, name, "heat_pipe", json!({"axis":0,"effective_resistance_k_per_w":0.5}))?,
        FeatureKind::ThermoelectricCooler => {
            let solids: Vec<String> = case
                .pointer("/thermal/solids")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .filter_map(|row| row["obstacle"].as_str().map(str::to_string))
                .collect();
            let [cold, hot, ..] = solids.as_slice() else {
                return Err(bad("a thermoelectric cooler needs two thermal solids (cold and hot side)"));
            };
            if case
                .pointer("/thermal/thermoelectric_coolers")
                .and_then(Value::as_array)
                .is_some_and(|rows| rows.iter().any(|row| row["name"].as_str() == Some(name)))
            {
                return Err(bad("cooler name already exists"));
            }
            let thermal = proposed.get_mut("thermal").ok_or_else(|| bad("thermal block missing"))?;
            let pointer = push(
                thermal,
                "thermoelectric_coolers",
                json!({"name":name,"cold_solid":cold,"hot_solid":hot,"heat_pumped_w":1.0,"electrical_power_w":0.5}),
                32,
            )?;
            format!("/thermal{pointer}")
        }
        FeatureKind::GlobalGoal | FeatureKind::PointGoal | FeatureKind::VolumeGoal | FeatureKind::EquationGoal => {
            if named(case, "goals", "name", name) {
                return Err(bad("goal name already exists"));
            }
            let goal = match kind {
                FeatureKind::GlobalGoal => {
                    if !thermal(case) {
                        return Err(bad("global temperature goals need a thermal block"));
                    }
                    json!({"name":name,"kind":"global_max_fluid_temperature"})
                }
                FeatureKind::PointGoal => {
                    let half = |axis: &str| case.pointer(&format!("/mesh/lengths/{axis}")).and_then(Value::as_f64).map(|l| l / 2.0);
                    let (Some(x), Some(y), Some(z)) = (half("x"), half("y"), half("z")) else {
                        return Err(bad("point goals need explicit mesh lengths"));
                    };
                    let kind = if thermal(case) { "point_temperature" } else { "point_pressure" };
                    json!({"name":name,"kind":kind,"point":[x,y,z]})
                }
                FeatureKind::VolumeGoal => {
                    let solid =
                        case.pointer("/thermal/solids/0/obstacle").and_then(Value::as_str).ok_or_else(|| bad("volume goals need a thermal solid"))?;
                    json!({"name":name,"kind":"solid_max_temperature","target":solid})
                }
                _ => {
                    let first = case
                        .get("goals")
                        .and_then(Value::as_array)
                        .and_then(|rows| rows.iter().find(|row| row["kind"].as_str() != Some("equation")))
                        .and_then(|row| row["name"].as_str())
                        .ok_or_else(|| bad("an equation goal needs another goal to refer to"))?;
                    json!({"name":name,"kind":"equation","expression":format!("{first} * 1")})
                }
            };
            push(&mut proposed, "goals", goal, 32)?
        }
    };
    validate_json(&proposed, MAX_CASE_BYTES)?;
    let issues = validate_case(&proposed).err().map(|e| vec![e.to_string()]).unwrap_or_default();
    Ok(PreparedFeature {
        case: proposed,
        pointer,
        issues,
        defaults_notice: format!(
            "Original illustrative starter values require review. No physical characterization or solver run is implied. Fan patch replacement explicitly accepted: {replace_fan_patch_conditions}."
        ),
    })
}
