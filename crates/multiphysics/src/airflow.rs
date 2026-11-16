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

pub const GRID_OBSTACLES: Requirement = Requirement { capability: "grid_obstacles_v1", versions: &["0.3.0", "0.4.0", "0.5.0", "0.6.0", "0.7.0"] };
pub const FAN_CURVE: Requirement = Requirement { capability: "fan_curve_v1", versions: &["0.4.0", "0.5.0", "0.6.0", "0.7.0"] };
pub const POROUS_BLOCKS: Requirement = Requirement { capability: "porous_blocks_v1", versions: &["0.4.0", "0.5.0", "0.6.0", "0.7.0"] };
pub const CONJUGATE_HEAT: Requirement = Requirement { capability: "conjugate_heat_v1", versions: &["0.5.0", "0.6.0", "0.7.0"] };
pub const SURFACE_BODIES: Requirement = Requirement { capability: "surface_bodies_v1", versions: &["0.6.0", "0.7.0"] };
pub const GOALS: Requirement = Requirement { capability: "goals_v1", versions: &["0.6.0", "0.7.0"] };
pub const NATURAL_CONVECTION: Requirement = Requirement { capability: "natural_convection_v1", versions: &["0.6.0", "0.7.0"] };
/// Electrical, mechanics, particle, kinetic and coupled domain cases added in 0.7.0.
pub const DOMAIN_CASES: Requirement = Requirement { capability: "domain_cases_v1", versions: &["0.7.0"] };

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
    out
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
}

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
    target: String,
}

#[derive(Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
enum GoalKindInput {
    PatchVolumeFlow,
    PatchMeanPressure,
    PatchHeatOut,
    SolidMaxTemperature,
    SolidMeanTemperature,
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
                EXTERIOR.contains(&goal.target.as_str()) || solids.contains(&goal.target)
            }
        };
        if !ok {
            return Err(Error::Invalid(format!("goal {} has an unknown target or needs a thermal block", goal.name)));
        }
    }
    Ok(())
}

pub(crate) fn validate(case: &Value) -> Result<()> {
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
