//! Conditions on body faces (`body_face_conditions_v1`): a block body inside a channel blows
//! fluid out of one of its faces. The block, 0.04 × 0.02 × 0.02 m, lies on whole cells of a
//! 0.01 m grid, so its +x face is exactly 4e-4 m² of grid faces and its other five faces
//! 3.6e-3 m². The face's inflow, 0.05 m/s × 4e-4 m² = 2e-5 m³/s, must leave through the outlet.

use cfd_solvers::case::*;
use cfd_solvers::case::conjugate::ConjugateHeat;
use cfd_solvers::core::Vec3;
use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};
use serde_json::json;

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("body-face-conditions-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        Self(dir)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

const SPEED: f64 = 0.05;
const FACE_AREA: f64 = 0.02 * 0.02;

/// A box from `lo` to `hi` as a closed outward triangle surface; triangles 2 and 3 are its +x face.
fn block(lo: [f64; 3], hi: [f64; 3]) -> SurfaceBody {
    let v = |i: usize| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }];
    let vertices: Vec<[f64; 3]> = (0..8).map(v).collect();
    // Faces -x, +x, -y, +y, -z, +z, each as two triangles wound outward.
    let quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]];
    let triangles = quads.iter().flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]).collect();
    SurfaceBody { name: "block".into(), vertices, triangles }
}

fn channel(conditions: Vec<BodyFaceCondition>) -> CollocatedFlowCase {
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec { dimensions: [20, 4, 4], lengths: Vec3::new(0.2, 0.04, 0.04), xy_shear: 0.0, coordinates: None },
        method: FlowMethod::Piso,
        output: PolyOutput { steps: 40, output_interval: 40 },
        bodies: vec![block([0.04, 0.01, 0.01], [0.08, 0.03, 0.03])],
        face_conditions: conditions,
        ..Default::default()
    };
    case.config.dt = 1e-2;
    case.config.kinematic_viscosity = 1.5e-5;
    case.config.pressure_tolerance = 1e-12;
    case.config.pressure_iterations = 4000;
    for patch in ["left", "bottom", "top", "front", "back"] {
        case.velocity_boundaries.push(VelocityPatch { patch: patch.into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::ZERO });
    }
    case.velocity_boundaries.push(VelocityPatch { patch: "right".into(), kind: VelocityBoundaryType::ZeroGradient, value: Vec3::ZERO });
    case.pressure_boundaries.push(PressurePatch { patch: "right".into(), kind: PressureBoundaryType::FixedValue, value: 0.0 });
    case
}

fn jet() -> BodyFaceCondition {
    BodyFaceCondition { name: "jet".into(), body: "block".into(), triangles: vec![2, 3] }
}

fn patch_report(dir: &std::path::Path, name: &str) -> serde_json::Value {
    let text = std::fs::read_to_string(dir.join("mesh_report.json")).unwrap();
    let report: serde_json::Value = serde_json::from_str(&text).unwrap();
    report["patches"].as_array().unwrap().iter().find(|p| p["name"] == name).cloned().unwrap_or_else(|| panic!("no patch {name}: {report}"))
}

#[test]
fn a_body_face_blows_its_flow_out_of_the_channel() {
    let mut case = channel(vec![jet()]);
    case.velocity_boundaries.push(VelocityPatch { patch: "jet".into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::new(SPEED, 0.0, 0.0) });
    case.pressure_boundaries.push(PressurePatch { patch: "jet".into(), kind: PressureBoundaryType::ZeroGradient, value: 0.0 });
    let scratch = Scratch::new("jet");
    let manifest = Case::CollocatedFlow3d(case).run(&scratch.0).unwrap();
    assert!(manifest.capabilities.iter().any(|c| c == "body_face_conditions_v1"));
    // The condition holds exactly the +x face; the block keeps its other five faces as wall.
    let jet = patch_report(&scratch.0, "jet");
    assert!((jet["area_m2"].as_f64().unwrap() - FACE_AREA).abs() < 1e-15, "{jet}");
    assert!((jet["outward_volume_flux_m3_per_s"].as_f64().unwrap() + SPEED * FACE_AREA).abs() < 1e-15, "{jet}");
    let wall = patch_report(&scratch.0, "block");
    assert!((wall["area_m2"].as_f64().unwrap() - 3.6e-3).abs() < 1e-15, "{wall}");
    assert!(wall["outward_volume_flux_m3_per_s"].as_f64().unwrap().abs() < 1e-15, "{wall}");
    // What the face blows in leaves through the outlet, to the pressure solve's continuity
    // residual after 40 steps (about 2e-7 of the flow here).
    let right = patch_report(&scratch.0, "right");
    let out = right["outward_volume_flux_m3_per_s"].as_f64().unwrap();
    assert!((out - SPEED * FACE_AREA).abs() < 1e-5 * SPEED * FACE_AREA, "{out}");
}

#[test]
fn body_face_conditions_are_checked() {
    let scratch = Scratch::new("checked");
    let run = |conditions: Vec<BodyFaceCondition>| Case::CollocatedFlow3d(channel(conditions)).run(&scratch.0).unwrap_err().to_string();
    let named = |name: &str, body: &str, triangles: Vec<usize>| BodyFaceCondition { name: name.into(), body: body.into(), triangles };
    assert!(run(vec![named("jet", "nobody", vec![2])]).contains("names no surface body"));
    assert!(run(vec![named("jet", "block", vec![99])]).contains("does not have"));
    assert!(run(vec![named("jet", "block", vec![])]).contains("lists no triangles"));
    assert!(run(vec![named("right", "block", vec![2])]).contains("unique ASCII identifiers"));
    assert!(run(vec![named("a", "block", vec![2]), named("b", "block", vec![2, 3])]).contains("in two body face conditions"));
    // With conjugate heat a body face that lets flow through is refused rather than left
    // uncoupled; a closed one is accepted (see the face heat source test).
    let mut thermal = channel(vec![jet()]);
    thermal.velocity_boundaries.push(VelocityPatch { patch: "jet".into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::new(SPEED, 0.0, 0.0) });
    thermal.thermal = Some(heat(json!([])));
    let why = Case::CollocatedFlow3d(thermal).run(&scratch.0).unwrap_err().to_string();
    assert!(why.contains("not supported with conjugate heat"), "{why}");
}

fn heat(face_sources: serde_json::Value) -> ConjugateHeat {
    serde_json::from_value(json!({
        "fluid": {"specific_heat_j_per_kg_k": 1005.0, "conductivity_w_per_m_k": 0.026},
        "inflow_temperature_k": 300.0, "solids": [{"obstacle": "block", "conductivity_w_per_m_k": [200.0, 200.0, 200.0], "power_w": 1.0}],
        "face_sources": face_sources
    }))
    .unwrap()
}

/// Heat on a body face (`face_heat_sources_v1`): 2 W on the block's +x face add to its 1 W of
/// volume power, the air carries all 3 W out, and a closed condition on that face (a wall to
/// monitor) runs with heat transfer.
#[test]
fn a_face_heat_source_heats_its_body_and_a_closed_face_condition_runs_with_heat() {
    let mut case = channel(vec![BodyFaceCondition { name: "monitor".into(), body: "block".into(), triangles: vec![2, 3] }]);
    case.velocity_boundaries.retain(|b| b.patch != "left");
    case.velocity_boundaries.push(VelocityPatch { patch: "left".into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::new(SPEED, 0.0, 0.0) });
    case.velocity_boundaries.push(VelocityPatch { patch: "monitor".into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::ZERO });
    case.thermal = Some(heat(json!([{"name": "chip", "body": "block", "triangles": [2, 3], "power_w": 2.0}])));
    let scratch = Scratch::new("face-heat");
    let manifest = Case::CollocatedFlow3d(case.clone()).run(&scratch.0).unwrap();
    assert!(manifest.capabilities.iter().any(|c| c == "face_heat_sources_v1"));
    let report: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(scratch.0.join("thermal_report.json")).unwrap()).unwrap();
    assert!((report["total_power_w"].as_f64().unwrap() - 3.0).abs() < 1e-9, "{report}");
    assert_eq!(report["face_sources"][0]["name"], "chip");
    let relative = report["relative_energy_imbalance"].as_f64().unwrap();
    assert!(relative.abs() < 1e-3, "{report}");
    // A source on no triangle, an unknown body, or a face the grid does not reach is refused.
    for (i, (sources, why)) in [
        (json!([{"name": "chip", "body": "block", "triangles": [], "power_w": 1.0}]), "lists no triangles"),
        (json!([{"name": "chip", "body": "nobody", "triangles": [2], "power_w": 1.0}]), "names no surface body"),
        (json!([{"name": "chip", "body": "block", "triangles": [99], "power_w": 1.0}]), "does not have"),
        (json!([{"name": "a", "body": "block", "triangles": [2], "power_w": 1.0}, {"name": "a", "body": "block", "triangles": [3], "power_w": 1.0}]), "unique names"),
    ]
    .into_iter()
    .enumerate()
    {
        let mut bad = case.clone();
        bad.thermal = Some(heat(sources));
        let error = Case::CollocatedFlow3d(bad).run(&scratch.0.join(i.to_string())).unwrap_err().to_string();
        assert!(error.contains(why), "{error}");
    }
}

#[test]
fn the_nearest_triangle_is_found_inside_and_outside_the_body() {
    let body = block([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    // Just off the +x face, just off the -z face, and far out beyond a corner of the -x face.
    let nearest = body.nearest_triangles(&[[1.1, 0.5, 0.2], [0.3, 0.6, -0.05], [-5.0, 0.1, 0.9]]).unwrap();
    assert!(matches!(nearest[0], 2 | 3), "{nearest:?}");
    assert!(matches!(nearest[1], 8 | 9), "{nearest:?}");
    assert!(matches!(nearest[2], 0 | 1), "{nearest:?}");
}
