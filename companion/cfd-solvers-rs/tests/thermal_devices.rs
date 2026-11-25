//! Thermal devices (`thermal_devices_v1`) and goal kinds (`goals_v2`) of the collocated flow
//! case: contact resistance, two-resistor junctions, heat pipes, thermoelectric coolers, and
//! global, point and equation goals. Each test checks a quantity against an independent
//! expectation (series resistance, exact formula, energy balance or a sibling report).
use cfd_solvers::case::conjugate::*;
use cfd_solvers::case::goals::*;
use cfd_solvers::{case::*, core::Vec3};
use std::sync::atomic::{AtomicUsize, Ordering};
static NEXT: AtomicUsize = AtomicUsize::new(0);

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "cfd-thermal-devices-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn fan_channel(curve: &[(f64, f64)]) -> CollocatedFlowCase {
    use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec {
            dimensions: [12, 6, 1],
            lengths: Vec3::new(0.12, 0.02, 0.01),
            xy_shear: 0.0,
            coordinates: None,
        },
        method: FlowMethod::Piso,
        output: PolyOutput {
            steps: 300,
            output_interval: 300,
        },
        ..Default::default()
    };
    case.config.dt = 0.01;
    case.config.density = 1.2;
    case.config.kinematic_viscosity = 1.5e-5;
    case.config.pressure_tolerance = 1e-12;
    case.config.pressure_iterations = 4000;
    for (name, kind) in [
        ("right", VelocityBoundaryType::ZeroGradient),
        ("bottom", VelocityBoundaryType::FixedValue),
        ("top", VelocityBoundaryType::FixedValue),
        ("front", VelocityBoundaryType::Slip),
        ("back", VelocityBoundaryType::Slip),
    ] {
        case.velocity_boundaries.push(VelocityPatch {
            patch: name.into(),
            kind,
            value: Vec3::ZERO,
        });
    }
    case.pressure_boundaries.push(PressurePatch {
        patch: "right".into(),
        kind: PressureBoundaryType::FixedValue,
        value: 0.0,
    });
    case.fan = Some(FanInlet {
        patch: "left".into(),
        reference_patch: "right".into(),
        pressure_convention: FanPressureConvention::Static,
        curve: curve
            .iter()
            .map(|&(q, p)| FanCurvePoint {
                volume_flow_m3_per_s: q,
                pressure_rise_pa: p,
            })
            .collect(),
        initial_volume_flow_m3_per_s: None,
        relaxation: 0.7,
        update_interval: 25,
        settle_tolerance: 2e-3,
    });
    case
}


fn solid(name: &str, k: f64, power: f64) -> ThermalSolid {
    ThermalSolid {
        obstacle: name.into(),
        conductivity_w_per_m_k: [k; 3],
        power_w: power,
        max_temperature_k: None,
        density_kg_per_m3: None,
        specific_heat_j_per_kg_k: None,
        contact_resistance_m2_k_per_w: None,
        junction_to_case_k_per_w: None,
        heat_pipe: None,
    }
}

/// A chip on a spreader in the fan channel: two touching solids, the chip powered.
fn chip_on_spreader() -> CollocatedFlowCase {
    let mut case = fan_channel(&[(0.0, 0.02), (2e-5, 0.0)]);
    case.obstacles.push(GridObstacle { name: "spreader".into(), min: [4, 0, 0], max: [8, 1, 1] });
    case.obstacles.push(GridObstacle { name: "chip".into(), min: [5, 1, 0], max: [7, 2, 1] });
    case.thermal = Some(ConjugateHeat {
        fluid: ThermalFluid { specific_heat_j_per_kg_k: 1006.0, conductivity_w_per_m_k: 0.026 },
        inflow_temperature_k: 298.15,
        solids: vec![solid("spreader", 200.0, 0.0), solid("chip", 100.0, 0.01)],
        boundaries: vec![ThermalPatch {
            patch: "bottom".into(),
            condition: ThermalCondition::Temperature { temperature_k: 298.15 },
        }],
        buoyancy: None,
        thermoelectric_coolers: vec![],
        face_sources: vec![],
        partial_cells: false,
    });
    case
}

fn run(case: &CollocatedFlowCase, scratch: &Scratch, name: &str) -> (RunManifest, serde_json::Value) {
    let dir = scratch.0.join(name);
    let manifest = Case::CollocatedFlow3d(case.clone()).run(&dir).unwrap();
    let report = serde_json::from_slice(&std::fs::read(dir.join("thermal_report.json")).unwrap()).unwrap();
    (manifest, report)
}

fn component<'a>(report: &'a serde_json::Value, name: &str) -> &'a serde_json::Value {
    report["components"].as_array().unwrap().iter().find(|c| c["name"] == name).unwrap()
}

fn heat_mut(case: &mut CollocatedFlowCase) -> &mut ConjugateHeat {
    case.thermal.as_mut().unwrap()
}

#[test]
fn contact_resistance_adds_the_series_temperature_step_between_solids() {
    let scratch = Scratch::new();
    let base = chip_on_spreader();
    let (manifest, plain) = run(&base, &scratch, "plain");
    assert!(manifest.capabilities.iter().any(|c| c == "thermal_devices_v1"));
    let mut resisted = base.clone();
    // 2 cm^2 K/W over the 2 x 1 cell contact of 0.01 x 0.01 m cells: 1e-4 m^2 total area.
    let r = 2e-4;
    heat_mut(&mut resisted).solids[1].contact_resistance_m2_k_per_w = Some(r);
    let (_, with) = run(&resisted, &scratch, "resisted");
    let chip = |report: &serde_json::Value| component(report, "chip")["mean_temperature_k"].as_f64().unwrap();
    let spreader = |report: &serde_json::Value| component(report, "spreader")["mean_temperature_k"].as_f64().unwrap();
    // Nearly all of the chip's power crosses the contact into the cooled spreader, so the
    // chip warms by about P R / A while the spreader barely changes.
    let expected = 0.01 * r / 2e-4;
    let rise = chip(&with) - chip(&plain);
    assert!(rise > 0.8 * expected && rise < 1.05 * expected, "rise {rise}, series estimate {expected}");
    assert!((spreader(&with) - spreader(&plain)).abs() < 0.1 * expected);
    // Energy is still conserved.
    assert!(with["relative_energy_imbalance"].as_f64().unwrap() < 1e-6, "{with:#}");
    // A contact resistance with fluid faces only (a lone solid) changes nothing.
    let mut lone = fan_channel(&[(0.0, 0.02), (2e-5, 0.0)]);
    lone.obstacles.push(GridObstacle { name: "chip".into(), min: [5, 0, 0], max: [7, 2, 1] });
    let mut a = solid("chip", 100.0, 0.01);
    lone.thermal = Some(ConjugateHeat {
        fluid: ThermalFluid { specific_heat_j_per_kg_k: 1006.0, conductivity_w_per_m_k: 0.026 },
        inflow_temperature_k: 298.15,
        solids: vec![a.clone()],
        boundaries: vec![],
        buoyancy: None,
        thermoelectric_coolers: vec![],
        face_sources: vec![],
        partial_cells: false,
    });
    let (_, before) = run(&lone, &scratch, "lone");
    a.contact_resistance_m2_k_per_w = Some(1.0);
    heat_mut(&mut lone).solids[0] = a;
    let (_, after) = run(&lone, &scratch, "lone_resisted");
    assert!((chip(&before) - chip(&after)).abs() < 1e-9);
    // Negative resistance is refused.
    let mut bad = base;
    heat_mut(&mut bad).solids[1].contact_resistance_m2_k_per_w = Some(-1.0);
    assert!(Case::CollocatedFlow3d(bad).run(&scratch.0.join("bad")).is_err());
}

#[test]
fn a_two_resistor_junction_is_reported_and_judged_above_its_case() {
    let scratch = Scratch::new();
    let mut case = chip_on_spreader();
    let (_, plain) = run(&case, &scratch, "plain");
    let case_max = component(&plain, "chip")["max_temperature_k"].as_f64().unwrap();
    // Limit between the case maximum and the junction: the case passes, the junction fails.
    let rjc = 100.0;
    let junction = case_max + 0.01 * rjc;
    heat_mut(&mut case).solids[1].junction_to_case_k_per_w = Some(rjc);
    heat_mut(&mut case).solids[1].max_temperature_k = Some(case_max + 0.5);
    let (manifest, report) = run(&case, &scratch, "junction");
    let chip = component(&report, "chip");
    assert!((chip["junction_temperature_k"].as_f64().unwrap() - junction).abs() < 1e-9);
    assert!((chip["max_temperature_k"].as_f64().unwrap() - case_max).abs() < 1e-9, "the junction resistance does not change the field");
    assert_eq!(chip["status"], "fail");
    assert_eq!(manifest.metrics["thermal_components_failing"], 1.0);
}

#[test]
fn a_heat_pipe_has_the_conductivity_of_its_effective_resistance() {
    let scratch = Scratch::new();
    let mut case = fan_channel(&[(0.0, 0.02), (2e-5, 0.0)]);
    // A 6-cell bar along x (0.06 m) over one 0.01 x 0.00333 m cell section.
    case.obstacles.push(GridObstacle { name: "pipe".into(), min: [3, 0, 0], max: [9, 1, 1] });
    case.thermal = Some(ConjugateHeat {
        fluid: ThermalFluid { specific_heat_j_per_kg_k: 1006.0, conductivity_w_per_m_k: 0.026 },
        inflow_temperature_k: 298.15,
        solids: vec![solid("pipe", 1.0, 0.01)],
        boundaries: vec![],
        buoyancy: None,
        thermoelectric_coolers: vec![],
        face_sources: vec![],
        partial_cells: false,
    });
    let (_, slow) = run(&case, &scratch, "bar");
    let r = 0.5;
    heat_mut(&mut case).solids[0].heat_pipe = Some(HeatPipe { axis: 0, effective_resistance_k_per_w: r });
    let (_, pipe) = run(&case, &scratch, "pipe");
    let length = 0.06;
    let area = 0.01 * (0.02 / 6.0);
    let k = component(&pipe, "pipe")["heat_pipe"]["axial_conductivity_w_per_m_k"].as_f64().unwrap();
    assert!((k - length / (r * area)).abs() / k < 1e-12, "k {k}");
    // Spreading along the pipe evens out the bar: a smaller max-min spread.
    let spread = |report: &serde_json::Value| {
        let c = component(report, "pipe");
        c["max_temperature_k"].as_f64().unwrap() - c["min_temperature_k"].as_f64().unwrap()
    };
    assert!(spread(&pipe) < 0.2 * spread(&slow), "{} vs {}", spread(&pipe), spread(&slow));
    heat_mut(&mut case).solids[0].heat_pipe = Some(HeatPipe { axis: 3, effective_resistance_k_per_w: r });
    assert!(Case::CollocatedFlow3d(case).run(&scratch.0.join("bad")).is_err());
}

#[test]
fn a_thermoelectric_cooler_pumps_heat_and_keeps_the_energy_balance() {
    let scratch = Scratch::new();
    let mut case = chip_on_spreader();
    // The chip is the cold side, the spreader the hot side.
    let (_, plain) = run(&case, &scratch, "plain");
    heat_mut(&mut case).thermoelectric_coolers.push(ThermoelectricCooler {
        name: "tec".into(),
        cold_solid: "chip".into(),
        hot_solid: "spreader".into(),
        heat_pumped_w: 0.008,
        electrical_power_w: 0.004,
    });
    let (manifest, report) = run(&case, &scratch, "tec");
    // Total heat into the domain is the chip power plus the electrical power.
    assert!((report["total_power_w"].as_f64().unwrap() - 0.014).abs() < 1e-12);
    assert!(report["relative_energy_imbalance"].as_f64().unwrap() < 1e-6);
    assert!((manifest.metrics["thermal_total_power_w"] - 0.014).abs() < 1e-12);
    let tec = &report["thermoelectric_coolers"][0];
    assert!((tec["cop"].as_f64().unwrap() - 2.0).abs() < 1e-12);
    assert!((tec["heat_released_w"].as_f64().unwrap() - 0.012).abs() < 1e-12);
    let chip = |r: &serde_json::Value| component(r, "chip")["mean_temperature_k"].as_f64().unwrap();
    assert!(chip(&report) < chip(&plain), "the cold side is cooled");
    // The same solid on both sides is refused.
    heat_mut(&mut case).thermoelectric_coolers[0].hot_solid = "chip".into();
    assert!(Case::CollocatedFlow3d(case).run(&scratch.0.join("bad")).is_err());
}

#[test]
fn global_point_and_equation_goals_agree_with_the_reports() {
    let scratch = Scratch::new();
    let mut case = chip_on_spreader();
    case.output = PolyOutput { steps: 600, output_interval: 600 };
    let goal = |name: &str, kind| Goal { name: name.into(), kind, target: String::new(), point: None, expression: None };
    let mut chip_max = goal("chip_max", GoalKind::SolidMaxTemperature);
    chip_max.target = "chip".into();
    let mut at_chip = goal("at_chip", GoalKind::PointTemperature);
    // The centre of chip cell (5, 1, 0).
    at_chip.point = Some([0.055, 0.005, 0.005]);
    let mut probe_p = goal("probe_p", GoalKind::PointPressure);
    probe_p.point = Some([0.115, 0.01, 0.005]);
    let mut rise = goal("chip_rise", GoalKind::Equation);
    rise.expression = Some("chip_max - 298.15".into());
    let mut ratio = goal("ratio", GoalKind::Equation);
    ratio.expression = Some("(fluid_max - 298.15) / -(-(chip_max - 298.15))".into());
    case.goals = vec![
        chip_max,
        goal("fluid_max", GoalKind::GlobalMaxFluidTemperature),
        goal("fluid_mean", GoalKind::GlobalMeanFluidTemperature),
        goal("solid_max", GoalKind::GlobalMaxSolidTemperature),
        at_chip,
        probe_p,
        rise,
        ratio,
    ];
    case.goal_control = Some(GoalControl { window_steps: 50, relative_tolerance: 1e-3, stop_when_converged: false, min_steps: 0 });
    let dir = scratch.0.join("goals");
    let manifest = Case::CollocatedFlow3d(case.clone()).run(&dir).unwrap();
    assert!(manifest.capabilities.iter().any(|c| c == "goals_v2"));
    let report: serde_json::Value = serde_json::from_slice(&std::fs::read(dir.join("thermal_report.json")).unwrap()).unwrap();
    let m = &manifest.metrics;
    let chip = component(&report, "chip");
    let spreader = component(&report, "spreader");
    let case_max = chip["max_temperature_k"].as_f64().unwrap().max(spreader["max_temperature_k"].as_f64().unwrap());
    assert!((m["goal_solid_max"] - case_max).abs() < 1e-9);
    assert!((m["goal_fluid_max"] - report["max_fluid_temperature_k"].as_f64().unwrap()).abs() < 1e-9);
    assert!(m["goal_fluid_mean"] > 298.15 && m["goal_fluid_mean"] <= m["goal_fluid_max"]);
    let (low, high) = (chip["min_temperature_k"].as_f64().unwrap(), chip["max_temperature_k"].as_f64().unwrap());
    assert!(m["goal_at_chip"] >= low - 1e-9 && m["goal_at_chip"] <= high + 1e-9);
    assert!(m["goal_probe_p"].is_finite());
    assert!((m["goal_chip_rise"] - (m["goal_chip_max"] - 298.15)).abs() < 1e-12);
    assert!((m["goal_ratio"] - (m["goal_fluid_max"] - 298.15) / (m["goal_chip_max"] - 298.15)).abs() < 1e-12);
    let goals: serde_json::Value = serde_json::from_slice(&std::fs::read(dir.join("goals.json")).unwrap()).unwrap();
    let units: Vec<&str> = goals["goals"].as_array().unwrap().iter().map(|g| g["unit"].as_str().unwrap()).collect();
    assert_eq!(units, ["K", "K", "K", "K", "K", "Pa", "", ""]);
    // Equations may only use non-equation goals that exist; points must lie in the domain.
    let mut bad = case.clone();
    bad.goals[6].expression = Some("ratio + 1".into());
    assert!(Case::CollocatedFlow3d(bad).run(&scratch.0.join("bad1")).is_err());
    let mut bad = case.clone();
    bad.goals[6].expression = Some("nope * 2".into());
    assert!(Case::CollocatedFlow3d(bad).run(&scratch.0.join("bad2")).is_err());
    let mut bad = case;
    bad.goals[4].point = Some([1.0, 0.0, 0.0]);
    assert!(Case::CollocatedFlow3d(bad).run(&scratch.0.join("bad3")).is_err());
}

#[test]
fn equations_parse_with_precedence_and_refuse_malformed_text() {
    let value = |name: &str| match name {
        "a" => 2.0,
        "b" => 3.0,
        _ => f64::NAN,
    };
    for (text, expected) in [
        ("a + b * 2", 8.0),
        ("(a + b) * 2", 10.0),
        ("-a - -b", 1.0),
        ("a / b * 3", 2.0),
        ("1.5e1 - a", 13.0),
        ("2e-1 * 10", 2.0),
    ] {
        let got = Expr::parse(text).unwrap().eval(&value);
        assert!((got - expected).abs() < 1e-12, "{text} = {got}");
    }
    assert!(Expr::parse("x").unwrap().eval(&value).is_nan());
    for bad in ["", "a +", "(a", "a b", "a $ b", "1..2"] {
        assert!(Expr::parse(bad).is_err(), "{bad} should be refused");
    }
    let deep = format!("{}a{}", "(".repeat(40), ")".repeat(40));
    assert!(Expr::parse(&deep).is_err());
}
