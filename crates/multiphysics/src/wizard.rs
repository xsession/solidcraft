//! Electronics enclosure wizard: a short engineering description in millimetres, degrees
//! Celsius and watts becomes a complete, admissible `collocated_flow3d` case with grid,
//! time step, boundaries, materials, goals and early stop chosen automatically.
//!
//! Every automatic choice is reported, including how far each component moved when it was
//! snapped to the grid, so the user can see what the solver will actually represent.
//! The output is an ordinary case: it can be edited, saved and run like any other.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::{Error, Result, validate_case};

/// Enclosure face by outward axis direction. Gravity points along `-z` unless stated.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum Face {
    #[serde(rename = "-x")]
    MinusX,
    #[serde(rename = "+x")]
    PlusX,
    #[serde(rename = "-y")]
    MinusY,
    #[serde(rename = "+y")]
    PlusY,
    #[serde(rename = "-z")]
    MinusZ,
    #[serde(rename = "+z")]
    PlusZ,
}

impl Face {
    fn patch(self) -> &'static str {
        match self {
            Self::MinusX => "left",
            Self::PlusX => "right",
            Self::MinusY => "bottom",
            Self::PlusY => "top",
            Self::MinusZ => "front",
            Self::PlusZ => "back",
        }
    }
    fn axis(self) -> usize {
        match self {
            Self::MinusX | Self::PlusX => 0,
            Self::MinusY | Self::PlusY => 1,
            Self::MinusZ | Self::PlusZ => 2,
        }
    }
    fn high(self) -> bool {
        matches!(self, Self::PlusX | Self::PlusY | Self::PlusZ)
    }
}

const FACES: [Face; 6] = [Face::MinusX, Face::PlusX, Face::MinusY, Face::PlusY, Face::MinusZ, Face::PlusZ];

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CurvePoint {
    pub litres_per_s: f64,
    pub pressure_pa: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Cooling {
    /// A fan blowing inward through a rectangular opening on one face. Air leaves through the
    /// whole exhaust face, optionally through a characterized filter layer.
    Fan {
        inlet: Face,
        /// Opening centre and size on the inlet face, in the face's two in-plane axes (mm).
        opening_center_mm: [f64; 2],
        opening_size_mm: [f64; 2],
        curve: Vec<CurvePoint>,
        static_pressure_curve: bool,
        exhaust: Face,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        exhaust_filter: Option<Filter>,
    },
    /// A sealed box cooled by natural convection inside and convection to ambient outside.
    Natural {
        #[serde(default = "minus_z")]
        gravity_toward: Face,
    },
}

fn minus_z() -> Face {
    Face::MinusZ
}

/// Characterized filter or grille layer (one cell thick across the exhaust face).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Filter {
    pub viscous_resistance_per_m2: f64,
    pub inertial_resistance_per_m: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Component {
    pub name: String,
    pub min_mm: [f64; 3],
    pub max_mm: [f64; 3],
    pub power_w: f64,
    /// Per-axis conductivity; use three equal values for an isotropic material.
    pub conductivity_w_per_m_k: [f64; 3],
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_temperature_c: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EnclosureSpec {
    pub size_mm: [f64; 3],
    pub cooling: Cooling,
    pub ambient_c: f64,
    /// External convection coefficient on walls that are not openings.
    pub wall_h_w_per_m2_k: f64,
    pub components: Vec<Component>,
    /// 1 (coarse, quick) to 5 (finest the companion admits).
    pub mesh_level: u8,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Snapped {
    pub name: String,
    pub requested_mm: [[f64; 3]; 2],
    pub represented_mm: [[f64; 3]; 2],
    pub max_shift_mm: f64,
    pub min_cells_across: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct WizardReport {
    pub dimensions: [u64; 3],
    pub cell_size_mm: [f64; 3],
    pub time_step_s: f64,
    pub steps_limit: u64,
    pub components: Vec<Snapped>,
    pub warnings: Vec<String>,
}

const AIR_DENSITY: f64 = 1.2;
const AIR_CP: f64 = 1006.0;
const AIR_K: f64 = 0.026;
const AIR_NU: f64 = 1.5e-5;
const CELL_BUDGET: [f64; 5] = [2_000.0, 4_000.0, 8_000.0, 16_000.0, 30_000.0];

fn finite_positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}

fn name_ok(name: &str) -> bool {
    let b = name.as_bytes();
    !b.is_empty()
        && b.len() <= 60
        && b.first().is_some_and(u8::is_ascii_alphanumeric)
        && b.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_')
        && !["left", "right", "bottom", "top", "front", "back"].contains(&name)
        && !name.starts_with("shroud_")
}

/// The two in-plane axes of a face, in increasing order.
fn plane(face: Face) -> [usize; 2] {
    match face.axis() {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

/// Build the case and its report. The returned case passes [`validate_case`].
pub fn enclosure_case(spec: &EnclosureSpec) -> Result<(Value, WizardReport)> {
    let invalid = |m: &str| Error::Invalid(format!("enclosure wizard: {m}"));
    if !spec.size_mm.iter().all(|v| finite_positive(*v) && *v <= 10_000.0) {
        return Err(invalid("enclosure size must be 0..10000 mm on each axis"));
    }
    if !(1..=5).contains(&spec.mesh_level) {
        return Err(invalid("mesh level must be 1..5"));
    }
    if !(spec.ambient_c.is_finite() && spec.ambient_c > -100.0 && spec.ambient_c < 500.0) || !finite_positive(spec.wall_h_w_per_m2_k) {
        return Err(invalid("ambient must be -100..500 °C and the wall coefficient positive"));
    }
    if spec.components.is_empty() || spec.components.len() > 100 {
        return Err(invalid("describe 1..100 components"));
    }
    let size: Vec<f64> = spec.size_mm.iter().map(|v| v / 1000.0).collect();
    // Near-cubic cells within the level's cell budget.
    let budget = CELL_BUDGET.get(usize::from(spec.mesh_level) - 1).copied().unwrap_or(2_000.0);
    let volume: f64 = size.iter().product();
    let h = (volume / budget).cbrt();
    let mut dims = [0u64; 3];
    for (d, len) in dims.iter_mut().zip(&size) {
        *d = ((len / h).round() as u64).clamp(4, 200);
    }
    while dims.iter().product::<u64>() > 32_768 {
        for d in &mut dims {
            *d = (*d * 9 / 10).max(4);
        }
    }
    let cell: Vec<f64> = size.iter().zip(dims).map(|(len, n)| len / n as f64).collect();
    let mut warnings = vec![];
    let mut obstacles = vec![];
    let mut solids = vec![];
    let mut snapped = vec![];
    let mut names = std::collections::BTreeSet::new();
    let ambient_k = spec.ambient_c + 273.15;
    let index = |axis: usize, mm: f64| -> u64 {
        let n = dims[axis];
        ((mm / 1000.0 / cell[axis]).round().max(0.0) as u64).min(n)
    };
    for c in &spec.components {
        if !name_ok(&c.name) || !names.insert(c.name.clone()) {
            return Err(invalid("component names must be unique identifiers (letters, digits, underscores)"));
        }
        if !(0..3).all(|a| {
            c.min_mm[a].is_finite() && c.max_mm[a].is_finite() && c.min_mm[a] >= 0.0 && c.min_mm[a] < c.max_mm[a] && c.max_mm[a] <= spec.size_mm[a]
        }) {
            return Err(invalid(&format!("component {} must lie inside the enclosure with min < max", c.name)));
        }
        if !(c.power_w.is_finite() && c.power_w >= 0.0) || !c.conductivity_w_per_m_k.iter().all(|k| finite_positive(*k)) {
            return Err(invalid(&format!("component {} needs non-negative power and positive conductivity", c.name)));
        }
        if c.max_temperature_c.is_some_and(|t| !(t.is_finite() && t > spec.ambient_c)) {
            return Err(invalid(&format!("component {} limit must exceed ambient", c.name)));
        }
        let mut lo = [0u64; 3];
        let mut hi = [0u64; 3];
        for a in 0..3 {
            lo[a] = index(a, c.min_mm[a]);
            hi[a] = index(a, c.max_mm[a]).max(lo[a] + 1).min(dims[a]);
            if hi[a] == lo[a] {
                lo[a] = hi[a] - 1;
            }
        }
        let represented = [[0, 1, 2].map(|a| lo[a] as f64 * cell[a] * 1000.0), [0, 1, 2].map(|a| hi[a] as f64 * cell[a] * 1000.0)];
        let shift = (0..3).flat_map(|a| [(represented[0][a] - c.min_mm[a]).abs(), (represented[1][a] - c.max_mm[a]).abs()]).fold(0.0, f64::max);
        let across = (0..3).map(|a| hi[a] - lo[a]).min().unwrap_or(0);
        if across < 2 {
            warnings.push(format!("{} is one cell thick on some axis; raise the mesh level to resolve its temperature gradient", c.name));
        }
        if shift > 0.25 * cell.iter().cloned().fold(f64::MAX, f64::min) * 1000.0 {
            warnings.push(format!("{} moved up to {shift:.2} mm when snapped to the grid", c.name));
        }
        obstacles.push(json!({"name": c.name, "min": lo, "max": hi}));
        let mut solid = json!({"obstacle": c.name, "conductivity_w_per_m_k": c.conductivity_w_per_m_k, "power_w": c.power_w});
        if let Some(limit) = c.max_temperature_c {
            solid["max_temperature_k"] = json!(limit + 273.15);
        }
        solids.push(solid);
        snapped.push(Snapped {
            name: c.name.clone(),
            requested_mm: [c.min_mm, c.max_mm],
            represented_mm: represented,
            max_shift_mm: shift,
            min_cells_across: across,
        });
    }
    let mut case = json!({
        "solver": "collocated_flow3d",
        "method": "piso",
        "mesh": {"dimensions": dims, "lengths": {"x": size[0], "y": size[1], "z": size[2]}},
        "config": {"density": AIR_DENSITY, "kinematic_viscosity": AIR_NU, "pressure_preconditioner": "jacobi",
            "pressure_iterations": 4000, "pressure_tolerance": 1e-13},
    });
    let h_min = cell.iter().cloned().fold(f64::MAX, f64::min);
    let power: f64 = spec.components.iter().map(|c| c.power_w).sum();
    let wall = |face: Face| json!({"patch": face.patch(), "condition": {"kind": "convection", "h_w_per_m2_k": spec.wall_h_w_per_m2_k, "ambient_k": ambient_k}});
    let speed;
    let mut goals = vec![];
    match &spec.cooling {
        Cooling::Fan { inlet, opening_center_mm, opening_size_mm, curve, static_pressure_curve, exhaust, exhaust_filter } => {
            if inlet == exhaust {
                return Err(invalid("fan inlet and exhaust faces must differ"));
            }
            if curve.len() < 2
                || curve.len() > 64
                || !curve.iter().all(|p| p.litres_per_s.is_finite() && p.pressure_pa.is_finite() && p.litres_per_s >= 0.0 && p.pressure_pa >= 0.0)
            {
                return Err(invalid("a fan curve needs 2..64 non-negative points"));
            }
            let [u, v] = plane(*inlet);
            let inside = |axis: usize, centre: f64, width: f64| {
                centre.is_finite() && width.is_finite() && width > 0.0 && centre - width / 2.0 >= 0.0 && centre + width / 2.0 <= spec.size_mm[axis]
            };
            if !inside(u, opening_center_mm[0], opening_size_mm[0]) || !inside(v, opening_center_mm[1], opening_size_mm[1]) {
                return Err(invalid("the fan opening must lie on the inlet face"));
            }
            let open_lo = [index(u, opening_center_mm[0] - opening_size_mm[0] / 2.0), index(v, opening_center_mm[1] - opening_size_mm[1] / 2.0)];
            let open_hi = [
                index(u, opening_center_mm[0] + opening_size_mm[0] / 2.0).max(open_lo[0] + 1),
                index(v, opening_center_mm[1] + opening_size_mm[1] / 2.0).max(open_lo[1] + 1),
            ];
            if !opening_size_mm.iter().all(|s| finite_positive(*s)) || open_hi[0] > dims[u] || open_hi[1] > dims[v] {
                return Err(invalid("the fan opening must lie on the inlet face"));
            }
            // A one-cell shroud layer leaves only the opening exposed on the inlet face.
            let layer = if inlet.high() { dims[inlet.axis()] - 1 } else { 0 };
            let mut add_shroud = |name: &str, umin: u64, umax: u64, vmin: u64, vmax: u64| {
                if umin >= umax || vmin >= vmax {
                    return;
                }
                let mut lo = [0u64; 3];
                let mut hi = [0u64; 3];
                lo[inlet.axis()] = layer;
                hi[inlet.axis()] = layer + 1;
                lo[u] = umin;
                hi[u] = umax;
                lo[v] = vmin;
                hi[v] = vmax;
                obstacles.push(json!({"name": name, "min": lo, "max": hi}));
                solids.push(json!({"obstacle": name, "conductivity_w_per_m_k": [0.2, 0.2, 0.2], "power_w": 0.0}));
            };
            add_shroud("shroud_low_u", 0, open_lo[0], 0, dims[v]);
            add_shroud("shroud_high_u", open_hi[0], dims[u], 0, dims[v]);
            add_shroud("shroud_low_v", open_lo[0], open_hi[0], 0, open_lo[1]);
            add_shroud("shroud_high_v", open_lo[0], open_hi[0], open_hi[1], dims[v]);
            for o in &obstacles {
                let axis = inlet.axis();
                let overlaps_shroud = o["name"].as_str().is_some_and(|n| !n.starts_with("shroud_"))
                    && o["min"][axis].as_u64().is_some_and(|min| min <= layer)
                    && o["max"][axis].as_u64().is_some_and(|max| layer < max);
                if overlaps_shroud {
                    return Err(invalid("components must not touch the fan's inlet face"));
                }
            }
            let mut points: Vec<Value> =
                curve.iter().map(|p| json!({"volume_flow_m3_per_s": p.litres_per_s / 1000.0, "pressure_rise_pa": p.pressure_pa})).collect();
            points.sort_by(|a, b| a["volume_flow_m3_per_s"].as_f64().unwrap_or(0.0).total_cmp(&b["volume_flow_m3_per_s"].as_f64().unwrap_or(0.0)));
            let q_max = curve.iter().map(|p| p.litres_per_s / 1000.0).fold(0.0, f64::max);
            let opening_area = (open_hi[0] - open_lo[0]) as f64 * cell[u] * (open_hi[1] - open_lo[1]) as f64 * cell[v];
            speed = q_max / opening_area;
            case["fan"] = json!({"patch": inlet.patch(), "reference_patch": exhaust.patch(),
                "pressure_convention": if *static_pressure_curve { "static" } else { "total" }, "curve": points});
            let walls: Vec<Value> = FACES
                .iter()
                .filter(|f| *f != inlet)
                .map(|f| {
                    if f == exhaust {
                        json!({"patch": f.patch(), "kind": "zero_gradient"})
                    } else {
                        json!({"patch": f.patch(), "kind": "fixed_value"})
                    }
                })
                .collect();
            case["velocity_boundaries"] = json!(walls);
            case["pressure_boundaries"] = json!([{"patch": exhaust.patch(), "kind": "fixed_value", "value": 0.0}]);
            if let Some(filter) = exhaust_filter {
                if !(filter.viscous_resistance_per_m2.is_finite()
                    && filter.inertial_resistance_per_m.is_finite()
                    && filter.viscous_resistance_per_m2 >= 0.0
                    && filter.inertial_resistance_per_m >= 0.0
                    && (filter.viscous_resistance_per_m2 > 0.0 || filter.inertial_resistance_per_m > 0.0))
                {
                    return Err(invalid("filter coefficients must be non-negative and not both zero"));
                }
                let a = exhaust.axis();
                let layer = if exhaust.high() { dims[a] - 1 } else { 0 };
                let mut lo = [0u64; 3];
                let mut hi = dims;
                lo[a] = layer;
                hi[a] = layer + 1;
                if obstacles
                    .iter()
                    .any(|o| (0..3).all(|ax| o["min"][ax].as_u64().unwrap_or(0) < hi[ax] && lo[ax] < o["max"][ax].as_u64().unwrap_or(0)))
                {
                    return Err(invalid("components must not touch the filtered exhaust face"));
                }
                case["porous_blocks"] = json!([{"name": "exhaust_filter", "min": lo, "max": hi,
                    "viscous_resistance_per_m2": filter.viscous_resistance_per_m2, "inertial_resistance_per_m": filter.inertial_resistance_per_m}]);
            }
            case["thermal"] = json!({"fluid": {"specific_heat_j_per_kg_k": AIR_CP, "conductivity_w_per_m_k": AIR_K},
                "inflow_temperature_k": ambient_k, "solids": solids,
                "boundaries": FACES.iter().filter(|f| *f != inlet && *f != exhaust).map(|f| wall(*f)).collect::<Vec<_>>()});
            goals.push(json!({"name": "fan_flow", "kind": "patch_volume_flow", "target": exhaust.patch()}));
        }
        Cooling::Natural { gravity_toward } => {
            let mut g = [0.0; 3];
            g[gravity_toward.axis()] = if gravity_toward.high() { 9.81 } else { -9.81 };
            // Rough natural-convection velocity scale for the time step: sqrt(g beta dT L).
            let rise_guess = (power / (spec.wall_h_w_per_m2_k * 2.0 * (size[0] * size[1] + size[1] * size[2] + size[0] * size[2]))).max(1.0);
            let beta = 1.0 / ambient_k;
            speed = (9.81 * beta * rise_guess * size.iter().cloned().fold(0.0, f64::max)).sqrt().max(0.01);
            case["velocity_boundaries"] = json!(FACES.iter().map(|f| json!({"patch": f.patch(), "kind": "fixed_value"})).collect::<Vec<_>>());
            case["pressure_boundaries"] = json!([]);
            case["thermal"] = json!({"fluid": {"specific_heat_j_per_kg_k": AIR_CP, "conductivity_w_per_m_k": AIR_K},
                "inflow_temperature_k": ambient_k, "solids": solids,
                "boundaries": FACES.iter().map(|f| wall(*f)).collect::<Vec<_>>(),
                "buoyancy": {"gravity_m_per_s2": g, "expansion_per_k": beta, "reference_temperature_k": ambient_k}});
            goals.push(json!({"name": "wall_heat", "kind": "patch_heat_out", "target": Face::PlusZ.patch()}));
        }
    }
    // Courant number <= 0.5 at the expected peak speed and viscous diffusion number <= 0.05.
    let dt_cfl = 0.5 * h_min / speed.max(1e-6);
    let dt_diff = 0.05 * h_min * h_min / AIR_NU;
    let dt = dt_cfl.min(dt_diff).min(0.1);
    let cells: u64 = dims.iter().product();
    let steps = (20_000u64).min(50_000_000 / cells.max(1)).max(400);
    case["config"]["dt"] = json!(dt);
    case["output"] = json!({"steps": steps, "output_interval": (steps / 4).max(1)});
    for c in &spec.components {
        goals.push(json!({"name": format!("{}_max", c.name), "kind": "solid_max_temperature", "target": c.name}));
    }
    case["obstacles"] = json!(obstacles);
    case["goals"] = json!(goals);
    case["goal_control"] = json!({"window_steps": 50, "relative_tolerance": 1e-3, "stop_when_converged": true, "min_steps": 300});
    validate_case(&case)?;
    Ok((
        case,
        WizardReport {
            dimensions: dims,
            cell_size_mm: [cell[0] * 1000.0, cell[1] * 1000.0, cell[2] * 1000.0],
            time_step_s: dt,
            steps_limit: steps,
            components: snapped,
            warnings,
        },
    ))
}
