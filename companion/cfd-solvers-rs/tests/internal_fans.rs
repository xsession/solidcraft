//! Internal fans (`internal_fans_v1`): an actuator block that follows its pressure-flow curve.
//! A slip-walled channel carries an internal fan, with ambient pressure at the inlet end and
//! a back-pressure at the outlet. The walls add no loss and the speed is the same at both
//! ends, so the fan must deliver exactly the back-pressure, at the flow its curve gives.

use cfd_solvers::case::internal_fan::{curve_pressure, FanDirection, InternalFan};
use cfd_solvers::case::*;
use cfd_solvers::core::Vec3;
use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("internal-fans-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        Self(dir)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

const RHO: f64 = 1.2;
const NU: f64 = 1.5e-5;
/// Channel 0.2 × 0.04 × 0.01 m on 20 × 4 × 1 cells; cross-section 4e-4 m².
const AREA: f64 = 0.04 * 0.01;

fn curve() -> Vec<FanCurvePoint> {
    [(0.0, 30.0), (2e-4, 22.0), (4e-4, 10.0), (5e-4, 0.0)]
        .into_iter()
        .map(|(q, p)| FanCurvePoint {
            volume_flow_m3_per_s: q,
            pressure_rise_pa: p,
        })
        .collect()
}

/// The channel with `back` pascals at the downstream end (right for +x, left for −x).
fn channel(direction: FanDirection) -> CollocatedFlowCase {
    channel_against(direction, 10.0)
}

fn channel_against(direction: FanDirection, back: f64) -> CollocatedFlowCase {
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec {
            dimensions: [20, 4, 1],
            lengths: Vec3::new(0.2, 0.04, 0.01),
            xy_shear: 0.0,
            coordinates: None,
        },
        method: FlowMethod::Piso,
        output: PolyOutput {
            steps: 1500,
            output_interval: 1500,
        },
        ..Default::default()
    };
    case.config.dt = 2e-3;
    case.config.density = RHO;
    case.config.kinematic_viscosity = NU;
    case.config.pressure_tolerance = 1e-10;
    case.config.pressure_iterations = 4000;
    case.config.piso_neighbour_update = true;
    case.config.balanced_body_force = true;
    for (patch, kind) in [
        ("left", VelocityBoundaryType::ZeroGradient),
        ("right", VelocityBoundaryType::ZeroGradient),
        ("bottom", VelocityBoundaryType::Slip),
        ("top", VelocityBoundaryType::Slip),
        ("front", VelocityBoundaryType::Slip),
        ("back", VelocityBoundaryType::Slip),
    ] {
        case.velocity_boundaries.push(VelocityPatch {
            patch: patch.into(),
            kind,
            value: Vec3::ZERO,
        });
    }
    let downstream = if direction == FanDirection::PlusX {
        "right"
    } else {
        "left"
    };
    for patch in ["left", "right"] {
        // Kinematic pressure: pascals over density.
        let value = if patch == downstream { back / RHO } else { 0.0 };
        case.pressure_boundaries.push(PressurePatch {
            patch: patch.into(),
            kind: PressureBoundaryType::FixedValue,
            value,
        });
    }
    case.internal_fans.push(InternalFan {
        name: "blower".into(),
        min: [4, 0, 0],
        max: [6, 4, 1],
        direction,
        curve: curve(),
        relaxation: 0.2,
    });
    case
}

/// Streamwise cell velocities of the middle row of the last frame.
fn row_speeds(dir: &std::path::Path) -> Vec<f64> {
    let mut frames: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().to_string()))
        .filter(|n| n.ends_with(".vtu"))
        .collect();
    frames.sort();
    let text = std::fs::read_to_string(dir.join(frames.last().unwrap())).unwrap();
    let start = text.find("Name=\"velocity_m_per_s\"").unwrap();
    let body = &text[start..];
    let body = &body[body.find('>').unwrap() + 1..body.find("</DataArray>").unwrap()];
    let v: Vec<f64> = body
        .split_whitespace()
        .map(|x| x.parse().unwrap())
        .collect();
    (0..20).map(|i| v[3 * (20 + i)]).collect()
}

fn report(dir: &std::path::Path) -> serde_json::Value {
    serde_json::from_slice(&std::fs::read(dir.join("internal_fans.json")).unwrap()).unwrap()
}

#[test]
fn an_internal_fan_delivers_the_flow_its_curve_gives_at_the_back_pressure() {
    for (back, expected) in [(10.0, 4e-4), (22.0, 2e-4), (26.0, 1e-4)] {
        let scratch = Scratch::new(&format!("back{back}"));
        let manifest = Case::CollocatedFlow3d(channel_against(FanDirection::PlusX, back))
            .run(&scratch.0)
            .unwrap();
        assert!(manifest
            .capabilities
            .iter()
            .any(|c| c == "internal_fans_v1"));
        let fan = &report(&scratch.0)["fans"][0];
        let q = fan["volume_flow_m3_per_s"].as_f64().unwrap();
        let rise = fan["pressure_rise_pa"].as_f64().unwrap();
        eprintln!("back {back} Pa: Q {q:.5e} (curve {expected:.5e}), rise {rise:.4} Pa");
        assert!((q - expected).abs() < 0.01 * expected, "{q} vs {expected}");
        assert!((rise - back).abs() < 0.01 * back, "{rise} vs {back}");
        assert!(
            (curve_pressure(&curve(), q).0 - rise).abs() < 0.01 * back,
            "settled on the curve"
        );
        assert_eq!(fan["within_curve_range"], true);
        assert_eq!(fan["cells"], 8);
        assert!((fan["area_m2"].as_f64().unwrap() - AREA).abs() < 1e-15);
        assert!((manifest.metrics["internal_fan_blower_volume_flow_m3_per_s"] - q).abs() < 1e-15);
    }
}

#[test]
fn balanced_forces_keep_the_cell_velocities_free_of_checkerboarding_at_the_fan_edges() {
    let (balanced, cell_gradient) = (Scratch::new("balanced"), Scratch::new("cellgrad"));
    Case::CollocatedFlow3d(channel_against(FanDirection::PlusX, 22.0))
        .run(&balanced.0)
        .unwrap();
    let mut plain = channel_against(FanDirection::PlusX, 22.0);
    plain.config.balanced_body_force = false;
    Case::CollocatedFlow3d(plain).run(&cell_gradient.0).unwrap();
    let spread = |v: &[f64]| {
        v.iter().cloned().fold(f64::MIN, f64::max) - v.iter().cloned().fold(f64::MAX, f64::min)
    };
    let (a, b) = (row_speeds(&balanced.0), row_speeds(&cell_gradient.0));
    eprintln!(
        "balanced spread {:.2e}, cell-gradient spread {:.2e}",
        spread(&a),
        spread(&b)
    );
    // Plug flow at 0.5 m/s everywhere with balanced forces.
    assert!(a.iter().all(|u| (u - 0.5).abs() < 1e-6), "{a:?}");
    // The cell-gradient form puts the force jump's error into the cells beside the block
    // (about d·f/4 = 0.002 s × 917 m/s² / 4), and the fan reads the wrong flow from it.
    assert!(spread(&b) > 0.5, "{b:?}");
    let q = report(&cell_gradient.0)["fans"][0]["volume_flow_m3_per_s"]
        .as_f64()
        .unwrap();
    assert!((q - 2e-4).abs() > 0.02 * 2e-4, "{q}");
}

#[test]
fn a_fan_against_more_than_its_shut_off_pressure_is_reported_off_the_curve() {
    let scratch = Scratch::new("stall");
    Case::CollocatedFlow3d(channel_against(FanDirection::PlusX, 40.0))
        .run(&scratch.0)
        .unwrap();
    let fan = &report(&scratch.0)["fans"][0];
    // The fan holds its 30 Pa shut-off rise and the 40 Pa back-pressure drives air backwards.
    assert!(fan["volume_flow_m3_per_s"].as_f64().unwrap() < 0.0);
    assert!((fan["pressure_rise_pa"].as_f64().unwrap() - 30.0).abs() < 0.3);
    assert_eq!(fan["within_curve_range"], false);
}

#[test]
fn turning_the_fan_round_reverses_the_flow() {
    let (a, b) = (Scratch::new("plus"), Scratch::new("minus"));
    Case::CollocatedFlow3d(channel(FanDirection::PlusX))
        .run(&a.0)
        .unwrap();
    Case::CollocatedFlow3d(channel(FanDirection::MinusX))
        .run(&b.0)
        .unwrap();
    // The flow is measured along each fan's own direction, so both report the same.
    let qa = report(&a.0)["fans"][0]["volume_flow_m3_per_s"]
        .as_f64()
        .unwrap();
    let qb = report(&b.0)["fans"][0]["volume_flow_m3_per_s"]
        .as_f64()
        .unwrap();
    assert!(qa > 0.0 && (qa - qb).abs() < 1e-3 * qa, "{qa} vs {qb}");
}

#[test]
fn bad_internal_fans_are_refused() {
    let run = |edit: &dyn Fn(&mut CollocatedFlowCase), tag: &str| {
        let mut case = channel(FanDirection::PlusX);
        case.output = PolyOutput {
            steps: 2,
            output_interval: 2,
        };
        edit(&mut case);
        let scratch = Scratch::new(tag);
        Case::CollocatedFlow3d(case).run(&scratch.0).is_err()
    };
    assert!(run(&|c| c.internal_fans[0].max = [21, 4, 1], "outside"));
    assert!(run(&|c| c.internal_fans[0].min = [6, 0, 0], "empty"));
    assert!(run(&|c| c.internal_fans[0].curve.reverse(), "rising"));
    assert!(run(&|c| c.internal_fans[0].curve.truncate(1), "short"));
    assert!(run(&|c| c.internal_fans[0].relaxation = 0.0, "relax"));
    assert!(run(
        &|c| c.porous_blocks.push(PorousBlock {
            name: "filter".into(),
            min: [5, 0, 0],
            max: [7, 4, 1],
            viscous_resistance_per_m2: 1e6,
            inertial_resistance_per_m: 0.0,
        }),
        "porous"
    ));
    assert!(run(
        &|c| c.obstacles.push(GridObstacle {
            name: "box".into(),
            min: [5, 1, 0],
            max: [6, 2, 1]
        }),
        "solid"
    ));
    assert!(run(
        &|c| {
            let mut twin = c.internal_fans[0].clone();
            twin.min = [5, 0, 0];
            twin.max = [7, 4, 1];
            twin.name = "twin".into();
            c.internal_fans.push(twin);
        },
        "overlap"
    ));
    // The JSON spells directions as signed axes and rejects unknown members.
    let json = serde_json::to_value(channel(FanDirection::MinusY)).unwrap();
    assert_eq!(json["internal_fans"][0]["direction"], "-y");
    let mut bad = json.clone();
    bad["internal_fans"][0]["swirl"] = serde_json::json!(1.0);
    assert!(serde_json::from_value::<CollocatedFlowCase>(bad).is_err());
}
