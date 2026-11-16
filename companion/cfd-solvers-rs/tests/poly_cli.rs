use cfd_solvers::{case::*, core::Vec3, io::CellType};
use std::sync::atomic::{AtomicUsize, Ordering};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "cfd-poly-cli-{}-{}",
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
#[test]
fn sheared_hex_exports_match_actual_control_volumes() {
    let spec = PolyMeshSpec {
        dimensions: [3, 4, 2],
        lengths: Vec3::new(1.2, 0.8, 0.3),
        xy_shear: 0.35,
    };
    let (mesh, export) = spec.build().unwrap();
    assert_eq!(mesh.cell_count(), export.cells.len());
    for (c, hex) in mesh.cells().iter().zip(&export.cells) {
        assert_eq!(hex.cell_type, CellType::Hexahedron);
        let center = hex
            .nodes
            .iter()
            .map(|&n| export.points[n])
            .fold(Vec3::ZERO, |a, b| a + b)
            / 8.0;
        assert!((center - c.center).norm() < 1e-14);
        let a = export.points[hex.nodes[0]];
        let ex = export.points[hex.nodes[1]] - a;
        let ey = export.points[hex.nodes[3]] - a;
        let ez = export.points[hex.nodes[4]] - a;
        assert!((ex.dot(ey.cross(ez)) - c.volume).abs() < 1e-14);
    }
    assert!(PolyMeshSpec {
        dimensions: [0, 2, 2],
        ..spec
    }
    .build()
    .is_err());
}
#[test]
fn invalid_new_cases_and_linear_failure_leave_no_success_manifest() {
    let f = Scratch::new();
    assert!(serde_json::from_str::<Case>(
        r#"{"solver":"poly_scalar_transport3d","diffusvity":1.0}"#
    )
    .is_err());
    assert!(serde_json::from_str::<Case>(
        r#"{"solver":"solid_heat3d","mesh":{"xy_shear":0.1,"typo":true}}"#
    )
    .is_err());
    assert!(serde_json::from_str::<Case>(r#"{"solver":"rans_transport3d","turbulence":{"model":"k_epsilon","k":0.1,"epsilon":0.01,"unused":true}}"#)
        .is_err());
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("../examples/flow-cavity-pimple.json")).unwrap();
    value["config"]["momentum_iterations"] = 1.into();
    value["config"]["momentum_tolerance"] = 1e-30.into();
    let case: Case = serde_json::from_value(value).unwrap();
    let out = f.0.join("failure");
    assert!(case.run(&out).is_err());
    assert!(!out.join("manifest.json").exists());
    let case = Case::PolyScalarTransport3d(PolyScalarCase {
        output: PolyOutput {
            steps: 1,
            output_interval: 0,
        },
        ..Default::default()
    });
    assert!(case.run(&f.0.join("interval")).is_err());
}
#[test]
fn simple_iteration_completion_is_distinct_from_steady_convergence() {
    let f = Scratch::new();
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("../examples/flow-cavity-simple.json")).unwrap();
    value["output"]["steps"] = 2.into();
    value["steady_velocity_tolerance"] = 1e-30.into();
    let case: Case = serde_json::from_value(value).unwrap();
    let m = case.run(&f.0).unwrap();
    assert!(m.execution_completed);
    assert_eq!(m.linear_solver_converged, Some(true));
    assert_eq!(m.steady_state_converged, Some(false));
    assert_eq!(m.metrics["time_s"], 0.0);
    let pvd = std::fs::read_to_string(f.0.join("steady_iterations.pvd")).unwrap();
    assert!(pvd.contains("timestep=\"2.00000000000000000e0\""));
}

fn blocked_flow() -> CollocatedFlowCase {
    CollocatedFlowCase {
        mesh: PolyMeshSpec {
            dimensions: [8, 6, 4],
            lengths: Vec3::new(0.08, 0.06, 0.04),
            xy_shear: 0.0,
        },
        obstacles: vec![GridObstacle {
            name: "component".into(),
            min: [3, 2, 1],
            max: [5, 4, 3],
        }],
        ..Default::default()
    }
}

#[test]
fn obstacle_mesh_excludes_solids_and_shares_exact_export_order() {
    let case = blocked_flow();
    let (mesh, export) = case.build_mesh().unwrap();
    assert_eq!(mesh.cell_count(), 8 * 6 * 4 - 8);
    assert_eq!(mesh.cell_count(), export.cells.len());
    let patch = mesh.patch_index("component").unwrap();
    assert_eq!(
        mesh.faces()
            .iter()
            .filter(|f| f.patch == Some(patch))
            .count(),
        24
    );
    for (ci, (cell, hex)) in mesh.cells().iter().zip(&export.cells).enumerate() {
        let center = hex
            .nodes
            .iter()
            .map(|&n| export.points[n])
            .fold(Vec3::ZERO, |a, b| a + b)
            / 8.0;
        assert!((center - cell.center).norm() < 1e-14);
        let p = cell.center;
        assert!(
            !(p.x > 0.03 && p.x < 0.05 && p.y > 0.02 && p.y < 0.04 && p.z > 0.01 && p.z < 0.03)
        );
        let closure = mesh.cell_faces()[ci]
            .iter()
            .map(|&fi| mesh.faces()[fi].outward_area(ci))
            .fold(Vec3::ZERO, |a, b| a + b);
        assert!(closure.norm() < 1e-16);
    }
    for f in mesh.faces() {
        if let Some(n) = f.neighbour {
            assert!(
                f.area
                    .dot(mesh.cells()[n].center - mesh.cells()[f.owner].center)
                    > 0.0
            );
        }
    }
    let legacy: Case = serde_json::from_str(include_str!("../examples/flow-channel.json")).unwrap();
    let Case::CollocatedFlow3d(legacy) = legacy else {
        panic!("wrong case")
    };
    assert!(legacy.obstacles.is_empty());
    let (old, _) = legacy.mesh.build().unwrap();
    let (compatible, _) = legacy.build_mesh().unwrap();
    assert_eq!(old.cell_count(), compatible.cell_count());
    assert_eq!(old.face_count(), compatible.face_count());
    assert!(serde_json::to_value(legacy)
        .unwrap()
        .get("obstacles")
        .is_none());
}

#[test]
fn invalid_obstacle_domains_are_rejected_before_allocation_or_solving() {
    let mut case = blocked_flow();
    case.obstacles.push(case.obstacles[0].clone());
    assert!(case.build_mesh().is_err());
    case = blocked_flow();
    case.obstacles[0].name = "right".into();
    assert!(case.build_mesh().is_err());
    case.obstacles[0].name = "bad/name".into();
    assert!(case.build_mesh().is_err());
    case = blocked_flow();
    case.obstacles[0].max = [9, 4, 3];
    assert!(case.build_mesh().is_err());
    case = blocked_flow();
    case.obstacles[0].min = [3, 0, 0];
    case.obstacles[0].max = [4, 6, 4];
    assert!(case
        .build_mesh()
        .unwrap_err()
        .to_string()
        .contains("disconnected"));
    case.obstacles[0].min = [0, 0, 0];
    case.obstacles[0].max = [8, 6, 4];
    assert!(case.build_mesh().is_err());
    case = blocked_flow();
    case.mesh.dimensions = [usize::MAX, usize::MAX, 2];
    assert!(case.build_mesh().is_err());
    case.mesh.dimensions = [64, 64, 64];
    assert!(case.build_mesh().is_err());
    case = blocked_flow();
    case.mesh.xy_shear = 0.1;
    assert!(case.build_mesh().is_err());
    case = blocked_flow();
    case.obstacles[0].min = [7, 0, 0];
    case.obstacles[0].max = [8, 6, 4];
    case.pressure_boundaries.push(PressurePatch {
        patch: "right".into(),
        kind: cfd_solvers::fvm::PressureBoundaryType::FixedValue,
        value: 0.0,
    });
    let scratch = Scratch::new();
    assert!(Case::CollocatedFlow3d(case)
        .run(&scratch.0.join("blocked_outlet"))
        .unwrap_err()
        .to_string()
        .contains("every face"));
}

#[test]
fn obstacle_flow_records_real_conservation_and_stationary_wall_flux() {
    use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};
    let scratch = Scratch::new();
    let mut case = blocked_flow();
    case.output = PolyOutput {
        steps: 4,
        output_interval: 2,
    };
    case.config.dt = 0.002;
    case.config.kinematic_viscosity = 1.5e-5;
    case.config.pressure_tolerance = 1e-13;
    case.initial_velocity = Vec3::new(0.01, 0.0, 0.0);
    for name in ["left", "right", "bottom", "top", "front", "back"] {
        case.velocity_boundaries.push(VelocityPatch {
            patch: name.into(),
            kind: if name == "right" {
                VelocityBoundaryType::ZeroGradient
            } else {
                VelocityBoundaryType::FixedValue
            },
            value: if name == "left" {
                case.initial_velocity
            } else {
                Vec3::ZERO
            },
        });
    }
    case.pressure_boundaries.push(PressurePatch {
        patch: "right".into(),
        kind: PressureBoundaryType::FixedValue,
        value: 0.0,
    });
    let manifest = Case::CollocatedFlow3d(case.clone())
        .run(&scratch.0.join("flow"))
        .unwrap();
    assert!(manifest.execution_completed);
    assert_eq!(manifest.linear_solver_converged, Some(true));
    assert_eq!(manifest.metrics["fluid_cell_count"], 184.0);
    assert_eq!(manifest.metrics["excluded_solid_cell_count"], 8.0);
    assert!(manifest.metrics["boundary_inflow_m3_per_s"] > 0.0);
    assert!(
        manifest.metrics["relative_boundary_flux_imbalance"] < 1e-6,
        "{:?}",
        manifest.metrics
    );
    let report: serde_json::Value =
        serde_json::from_slice(&std::fs::read(scratch.0.join("flow/mesh_report.json")).unwrap())
            .unwrap();
    let wall = report["patches"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "component")
        .unwrap();
    assert_eq!(wall["face_count"], 24);
    assert_eq!(wall["outward_volume_flux_m3_per_s"], 0.0);
    case.velocity_boundaries.push(VelocityPatch {
        patch: "component".into(),
        kind: VelocityBoundaryType::Slip,
        value: Vec3::ZERO,
    });
    assert!(Case::CollocatedFlow3d(case)
        .run(&scratch.0.join("invalid_wall"))
        .is_err());
    assert!(!scratch.0.join("invalid_wall/manifest.json").exists());
}

fn fan_channel(curve: &[(f64, f64)]) -> CollocatedFlowCase {
    use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec {
            dimensions: [12, 6, 1],
            lengths: Vec3::new(0.12, 0.02, 0.01),
            xy_shear: 0.0,
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

fn fan_point(
    case: &CollocatedFlowCase,
    dir: &std::path::Path,
) -> (f64, f64, serde_json::Value, RunManifest) {
    let manifest = Case::CollocatedFlow3d(case.clone()).run(dir).unwrap();
    let report: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("fan_operating_point.json")).unwrap())
            .unwrap();
    (
        report["volume_flow_m3_per_s"].as_f64().unwrap(),
        report["system_pressure_rise_pa"].as_f64().unwrap(),
        report,
        manifest,
    )
}

#[test]
fn fan_curve_lookup_never_extrapolates() {
    let fan = fan_channel(&[(0.0, 0.02), (1e-5, 0.015), (2e-5, 0.0)])
        .fan
        .unwrap();
    assert!(fan.validate().is_ok());
    assert_eq!(fan.pressure_at(5e-6), 0.0175);
    let (q, within) = fan.flow_at(0.0175);
    assert!(within && (q - 5e-6).abs() < 1e-18);
    assert_eq!(fan.flow_at(0.05), (0.0, false));
    assert_eq!(fan.flow_at(-1.0), (2e-5, false));
    // 0.015 - k q = 0 crossing on the second segment for k = 1000: q = 1.2e-5.
    assert!((fan.intersect_resistance(1000.0) - 1.2e-5).abs() < 1e-18);
    assert_eq!(fan.intersect_resistance(f64::INFINITY), 0.0);
    assert_eq!(fan.intersect_resistance(-1.0), 2e-5);
    assert!(fan.intersect_resistance(1e12) < 1e-13);
    for curve in [
        vec![(0.0, 0.02)],
        vec![(0.0, 0.02), (1e-5, 0.03)],
        vec![(1e-5, 0.02), (1e-5, 0.01)],
        vec![(0.0, 0.02), (f64::NAN, 0.0)],
        vec![(-1e-5, 0.02), (1e-5, 0.0)],
    ] {
        let points: Vec<(f64, f64)> = curve;
        assert!(
            fan_channel(&points).fan.unwrap().validate().is_err(),
            "{points:?}"
        );
    }
    let mut bad = fan;
    bad.relaxation = 0.0;
    assert!(bad.validate().is_err());
    bad.relaxation = 0.5;
    bad.initial_volume_flow_m3_per_s = Some(1.0);
    assert!(bad.validate().is_err());
}

#[test]
fn ambiguous_fan_boundaries_are_rejected_without_a_manifest() {
    use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};
    let scratch = Scratch::new();
    let curve = [(0.0, 0.02), (2e-5, 0.0)];
    let mut owned = fan_channel(&curve);
    owned.velocity_boundaries.push(VelocityPatch {
        patch: "left".into(),
        kind: VelocityBoundaryType::FixedValue,
        value: Vec3::ZERO,
    });
    let mut no_reference = fan_channel(&curve);
    no_reference.pressure_boundaries.clear();
    let mut fixed_fan_pressure = fan_channel(&curve);
    fixed_fan_pressure.pressure_boundaries.push(PressurePatch {
        patch: "left".into(),
        kind: PressureBoundaryType::FixedValue,
        value: 0.0,
    });
    let mut same_patch = fan_channel(&curve);
    if let Some(fan) = same_patch.fan.as_mut() {
        fan.reference_patch = "left".into();
    }
    for (name, case) in [
        ("owned", owned),
        ("reference", no_reference),
        ("fixed", fixed_fan_pressure),
        ("same", same_patch),
    ] {
        let dir = scratch.0.join(name);
        assert!(Case::CollocatedFlow3d(case).run(&dir).is_err(), "{name}");
        assert!(!dir.join("manifest.json").exists());
    }
    let json = r#"{"solver":"collocated_flow3d","fan":{"patch":"left","reference_patch":"right","pressure_convention":"static","curve":[],"typo":1}}"#;
    assert!(serde_json::from_str::<Case>(json).is_err());
    let json = r#"{"solver":"collocated_flow3d","fan":{"patch":"left","reference_patch":"right","curve":[]}}"#;
    assert!(
        serde_json::from_str::<Case>(json).is_err(),
        "the pressure convention must be explicit"
    );
}

#[test]
fn fan_operating_point_follows_the_system_resistance() {
    let scratch = Scratch::new();
    let weak = [(0.0, 0.02), (2e-5, 0.0)];
    let strong = [(0.0, 0.04), (1e-5, 0.03), (3e-5, 0.0)];
    let (q_weak, p_weak, report, manifest) =
        fan_point(&fan_channel(&weak), &scratch.0.join("weak"));
    assert!(manifest.capabilities.iter().any(|c| c == "fan_curve_v1"));
    assert!(manifest
        .files
        .iter()
        .any(|f| f == "fan_operating_point.csv"));
    assert_eq!(manifest.metrics["fan_operating_point_within_curve"], 1.0);
    assert!(
        report["relative_operating_point_mismatch"]
            .as_f64()
            .unwrap()
            < 0.01,
        "{report}"
    );
    assert!(report["corrections"].as_u64().unwrap() >= 3);
    let (q_strong, p_strong, report, _) =
        fan_point(&fan_channel(&strong), &scratch.0.join("strong"));
    assert!(
        report["relative_operating_point_mismatch"]
            .as_f64()
            .unwrap()
            < 0.01,
        "{report}"
    );
    assert!(q_strong > 1.2 * q_weak);
    // Each operating point must lie on the system curve measured independently: hold the
    // same flow fixed (no corrections within the run) and compare the settled pressure rise.
    for (name, curve, q, p) in [
        ("weak_fixed", &weak[..], q_weak, p_weak),
        ("strong_fixed", &strong[..], q_strong, p_strong),
    ] {
        let mut fixed = fan_channel(curve);
        if let Some(fan) = fixed.fan.as_mut() {
            fan.initial_volume_flow_m3_per_s = Some(q);
            fan.update_interval = 10_000;
        }
        let (q_fixed, p_fixed, report, _) = fan_point(&fixed, &scratch.0.join(name));
        assert_eq!(report["corrections"], 0);
        assert_eq!(q_fixed, q);
        assert!((p_fixed - p).abs() < 0.01 * p, "{name}: {p_fixed} {p}");
    }
    // Adding resistance moves the same fan back along its curve.
    let mut blocked = fan_channel(&weak);
    blocked.obstacles.push(GridObstacle {
        name: "baffle".into(),
        min: [6, 0, 0],
        max: [7, 4, 1],
    });
    let (q_blocked, p_blocked, _, _) = fan_point(&blocked, &scratch.0.join("blocked"));
    assert!(
        q_blocked < 0.9 * q_weak && p_blocked > p_weak,
        "{q_blocked} {q_weak}"
    );
    // A total-pressure curve must also overcome fan-face dynamic pressure.
    let mut total = fan_channel(&weak);
    if let Some(fan) = total.fan.as_mut() {
        fan.pressure_convention = FanPressureConvention::Total;
    }
    let (q_total, _, _, _) = fan_point(&total, &scratch.0.join("total"));
    assert!(q_total < q_weak);
}

fn porous_plug(viscous: f64, inertial: f64) -> CollocatedFlowCase {
    use cfd_solvers::fvm::VelocityBoundaryType;
    let mut case = fan_channel(&[(0.0, 1.0), (1e-3, 0.0)]);
    case.mesh = PolyMeshSpec {
        dimensions: [20, 2, 2],
        lengths: Vec3::new(0.2, 0.02, 0.02),
        xy_shear: 0.0,
    };
    case.output = PolyOutput {
        steps: 200,
        output_interval: 200,
    };
    case.config.dt = 0.02;
    for b in &mut case.velocity_boundaries {
        if b.patch != "right" {
            b.kind = VelocityBoundaryType::Slip;
        }
    }
    if let Some(fan) = case.fan.as_mut() {
        // A fixed flow: 0.05 m/s through 4e-4 m^2. The fan only measures the pressure rise.
        fan.initial_volume_flow_m3_per_s = Some(2e-5);
        fan.update_interval = 10_000;
    }
    case.porous_blocks.push(PorousBlock {
        name: "vent".into(),
        min: [8, 0, 0],
        max: [12, 2, 2],
        viscous_resistance_per_m2: viscous,
        inertial_resistance_per_m: inertial,
    });
    case
}

#[test]
fn porous_plug_matches_darcy_forchheimer_pressure_drop_and_power() {
    let scratch = Scratch::new();
    for (name, d, c) in [
        ("darcy", 1e6, 0.0),
        ("forchheimer", 0.0, 400.0),
        ("both", 1e6, 100.0),
    ] {
        // The collocated interpolation at an abrupt resistance jump is first order: every
        // error must shrink under refinement and be small on the finer grid.
        let mut errors = vec![];
        for n in [20, 40] {
            let mut case = porous_plug(d, c);
            case.mesh.dimensions[0] = n;
            case.porous_blocks[0].min[0] = 2 * n / 5;
            case.porous_blocks[0].max[0] = 3 * n / 5;
            let dir = scratch.0.join(format!("{name}{n}"));
            let (q, rise, _, manifest) = fan_point(&case, &dir);
            assert!(manifest
                .capabilities
                .iter()
                .any(|c| c == "porous_blocks_v1"));
            assert!(manifest.files.iter().any(|f| f == "porous_zones.json"));
            let u = q / 4e-4;
            let beta = 1.5e-5 * d + 0.5 * c * u;
            let expected = 1.2 * beta * u * 0.04;
            let report: serde_json::Value =
                serde_json::from_slice(&std::fs::read(dir.join("porous_zones.json")).unwrap())
                    .unwrap();
            let zone = &report["zones"][0];
            assert_eq!(zone["cells"], 4 * n / 5);
            let speed = zone["mean_superficial_speed_m_per_s"].as_f64().unwrap();
            // Mechanical energy balance for uniform superficial flow: dissipated power = Q * dp.
            let power = report["total_dissipated_power_W"].as_f64().unwrap();
            errors.push([
                (rise - expected).abs() / expected,
                (speed - u).abs() / u,
                (power - q * expected).abs() / (q * expected),
            ]);
        }
        for k in 0..3 {
            assert!(
                errors[0][k] < 0.15 && errors[1][k] < 0.06 && errors[1][k] < 0.6 * errors[0][k],
                "{name}: {errors:?}"
            );
        }
    }
}

#[test]
fn invalid_porous_zones_are_rejected_and_a_vent_lowers_the_fan_flow() {
    let scratch = Scratch::new();
    let mut bad = vec![];
    let mut case = porous_plug(1e6, 0.0);
    case.porous_blocks[0].viscous_resistance_per_m2 = 0.0;
    bad.push(case);
    let mut case = porous_plug(-1.0, 0.0);
    case.porous_blocks[0].inertial_resistance_per_m = 1.0;
    bad.push(case);
    let mut case = porous_plug(1e6, 0.0);
    case.porous_blocks.push(case.porous_blocks[0].clone());
    case.porous_blocks[1].name = "second".into();
    bad.push(case);
    let mut case = porous_plug(1e6, 0.0);
    case.obstacles.push(GridObstacle {
        name: "block".into(),
        min: [9, 0, 0],
        max: [10, 1, 1],
    });
    bad.push(case);
    let mut case = porous_plug(1e6, 0.0);
    case.porous_blocks[0].name = "right".into();
    bad.push(case);
    let mut case = porous_plug(1e6, 0.0);
    case.porous_blocks[0].max = [21, 2, 2];
    bad.push(case);
    for (i, case) in bad.into_iter().enumerate() {
        let dir = scratch.0.join(format!("bad{i}"));
        assert!(Case::CollocatedFlow3d(case).run(&dir).is_err(), "case {i}");
        assert!(!dir.join("manifest.json").exists());
    }
    let curve = [(0.0, 0.02), (2e-5, 0.0)];
    let (open, _, _, _) = fan_point(&fan_channel(&curve), &scratch.0.join("open"));
    let mut vented = fan_channel(&curve);
    vented.porous_blocks.push(PorousBlock {
        name: "filter".into(),
        min: [10, 0, 0],
        max: [11, 6, 1],
        viscous_resistance_per_m2: 2e5,
        inertial_resistance_per_m: 0.0,
    });
    let (q, p, report, _) = fan_point(&vented, &scratch.0.join("vented"));
    assert!(
        report["relative_operating_point_mismatch"]
            .as_f64()
            .unwrap()
            < 0.01,
        "{report}"
    );
    assert!(q < 0.9 * open && p > 0.0);
}

mod conjugate_heat {
    use super::*;
    use cfd_solvers::case::conjugate::*;

    fn wall(patch: &str, temperature_k: f64) -> ThermalPatch {
        ThermalPatch {
            patch: patch.into(),
            condition: ThermalCondition::Temperature { temperature_k },
        }
    }

    #[test]
    fn energy_operator_is_second_order_for_a_manufactured_solution() {
        let mut errors = vec![];
        for n in [8usize, 16, 32] {
            let spec = PolyMeshSpec {
                dimensions: [n, n, 1],
                lengths: Vec3::new(1.0, 1.0, 0.1),
                xy_shear: 0.0,
            };
            let (mesh, _) = spec.build().unwrap();
            let pi = std::f64::consts::PI;
            let exact: Vec<f64> = mesh
                .cells()
                .iter()
                .map(|c| (pi * c.center.x).sin() * (pi * c.center.y).sin())
                .collect();
            let source: Vec<f64> = mesh
                .cells()
                .iter()
                .zip(&exact)
                .map(|(c, e)| 2.0 * pi * pi * e * c.volume)
                .collect();
            let patches: Vec<_> = ["left", "right", "bottom", "top"]
                .iter()
                .map(|p| wall(p, 300.0))
                .collect();
            let s = solve_energy(
                spec,
                &vec![[1.0; 3]; n * n],
                &vec![0.0; mesh.face_count()],
                &source,
                &patches,
                300.0,
                1.0,
            )
            .unwrap();
            let volume: f64 = mesh.cells().iter().map(|c| c.volume).sum();
            let l2 = (mesh
                .cells()
                .iter()
                .zip(&s.temperature_k)
                .zip(&exact)
                .map(|((c, t), e)| (t - 300.0 - e).powi(2) * c.volume)
                .sum::<f64>()
                / volume)
                .sqrt();
            errors.push(l2);
        }
        let orders: Vec<f64> = errors.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
        assert!(
            orders.iter().all(|o| (o - 2.0).abs() < 0.1),
            "{errors:?} {orders:?}"
        );
    }

    #[test]
    fn layered_conduction_is_exact_with_harmonic_interfaces() {
        let spec = PolyMeshSpec {
            dimensions: [6, 1, 1],
            lengths: Vec3::new(0.06, 0.01, 0.01),
            xy_shear: 0.0,
        };
        let (mesh, _) = spec.build().unwrap();
        let k = [1.0, 1.0, 0.1, 0.1, 10.0, 10.0];
        let conductivity: Vec<[f64; 3]> = k.iter().map(|k| [*k, 5.0, 5.0]).collect();
        let s = solve_energy(
            spec,
            &conductivity,
            &vec![0.0; mesh.face_count()],
            &[0.0; 6],
            &[wall("left", 400.0), wall("right", 300.0)],
            300.0,
            1.0,
        )
        .unwrap();
        let resistance: f64 = k.iter().map(|k| 0.01 / k).sum();
        let q = 100.0 / resistance * 1e-4;
        let out = |name: &str| {
            s.conduction_out_w
                .iter()
                .find(|(p, _)| p == name)
                .unwrap()
                .1
        };
        assert!(
            (out("right") - q).abs() < 1e-10 * q && (out("left") + q).abs() < 1e-10 * q,
            "{:?}",
            s.conduction_out_w
        );
        // Exact piecewise-linear profile at every cell centre.
        let mut x_t = 400.0;
        for (i, t) in s.temperature_k.iter().enumerate() {
            let centre = x_t - q / 1e-4 * 0.005 / k[i];
            assert!((t - centre).abs() < 1e-9, "cell {i}: {t} vs {centre}");
            x_t -= q / 1e-4 * 0.01 / k[i];
        }
    }

    fn heated_channel(power: f64, limit: Option<f64>) -> CollocatedFlowCase {
        let mut case = fan_channel(&[(0.0, 0.02), (2e-5, 0.0)]);
        case.obstacles.push(GridObstacle {
            name: "chip".into(),
            min: [5, 0, 0],
            max: [7, 2, 1],
        });
        case.thermal = Some(ConjugateHeat {
            fluid: ThermalFluid {
                specific_heat_j_per_kg_k: 1006.0,
                conductivity_w_per_m_k: 0.026,
            },
            inflow_temperature_k: 298.15,
            solids: vec![ThermalSolid {
                obstacle: "chip".into(),
                conductivity_w_per_m_k: [150.0; 3],
                power_w: power,
                max_temperature_k: limit,
                density_kg_per_m3: None,
                specific_heat_j_per_kg_k: None,
            }],
            boundaries: vec![],
            buoyancy: None,
        });
        case
    }

    fn report(dir: &std::path::Path) -> serde_json::Value {
        serde_json::from_slice(&std::fs::read(dir.join("thermal_report.json")).unwrap()).unwrap()
    }

    #[test]
    fn adiabatic_enclosure_advects_all_power_and_balances_energy() {
        let scratch = Scratch::new();
        let manifest = Case::CollocatedFlow3d(heated_channel(0.01, Some(400.0)))
            .run(&scratch.0.join("a"))
            .unwrap();
        assert!(manifest
            .capabilities
            .iter()
            .any(|c| c == "conjugate_heat_v1"));
        assert!(manifest.files.iter().any(|f| f == "thermal.vtu"));
        let r = report(&scratch.0.join("a"));
        assert!(
            r["relative_energy_imbalance"].as_f64().unwrap() < 1e-9,
            "{r}"
        );
        let advected = r["advected_heat_out_w"].as_f64().unwrap();
        assert!(
            (advected - 0.01).abs() < 1e-10,
            "all heat leaves with the flow: {advected}"
        );
        // Bulk rise = P / (rho cp Q) for adiabatic walls.
        let q = r["inflow_volume_m3_per_s"].as_f64().unwrap();
        let rise = r["bulk_outflow_rise_k"].as_f64().unwrap();
        assert!((rise - 0.01 / (1.2 * 1006.0 * q)).abs() < 1e-9 * rise);
        let chip = &r["components"][0];
        assert_eq!(chip["status"], "pass");
        let hot = chip["max_temperature_k"].as_f64().unwrap();
        assert!(hot > r["max_fluid_temperature_k"].as_f64().unwrap() && hot > 298.15);
        // The steady problem is linear in power, and the flow does not depend on temperature.
        let r2 = {
            Case::CollocatedFlow3d(heated_channel(
                0.02,
                Some(298.15 + 2.0 * (hot - 298.15) * 1.01),
            ))
            .run(&scratch.0.join("b"))
            .unwrap();
            report(&scratch.0.join("b"))
        };
        let hot2 = r2["components"][0]["max_temperature_k"].as_f64().unwrap();
        assert!(
            ((hot2 - 298.15) - 2.0 * (hot - 298.15)).abs() < 1e-8 * (hot - 298.15),
            "{hot} {hot2}"
        );
        assert_eq!(r2["components"][0]["status"], "warn");
        Case::CollocatedFlow3d(heated_channel(0.02, Some(hot)))
            .run(&scratch.0.join("c"))
            .unwrap();
        assert_eq!(
            report(&scratch.0.join("c"))["components"][0]["status"],
            "fail"
        );
    }

    #[test]
    fn cooled_walls_take_their_share_and_invalid_thermal_cases_are_refused() {
        let scratch = Scratch::new();
        let mut case = heated_channel(0.01, None);
        if let Some(t) = case.thermal.as_mut() {
            t.boundaries.push(ThermalPatch {
                patch: "bottom".into(),
                condition: ThermalCondition::Convection {
                    h_w_per_m2_k: 50.0,
                    ambient_k: 298.15,
                },
            });
        }
        Case::CollocatedFlow3d(case)
            .run(&scratch.0.join("cooled"))
            .unwrap();
        let r = report(&scratch.0.join("cooled"));
        let conducted: f64 = r["conducted_heat_out_w"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| p["w"].as_f64().unwrap())
            .sum();
        assert!(conducted > 0.0 && r["advected_heat_out_w"].as_f64().unwrap() < 0.01);
        assert!(r["relative_energy_imbalance"].as_f64().unwrap() < 1e-9);
        assert_eq!(r["components"][0]["status"], "no_limit");
        let edits: Vec<Box<dyn Fn(&mut ConjugateHeat)>> = vec![
            Box::new(|t| t.solids.clear()),
            Box::new(|t| t.solids[0].obstacle = "ghost".into()),
            Box::new(|t| t.solids[0].conductivity_w_per_m_k[2] = 0.0),
            Box::new(|t| t.solids[0].power_w = -1.0),
            Box::new(|t| t.solids.push(t.solids[0].clone())),
            Box::new(|t| t.inflow_temperature_k = -1.0),
            Box::new(|t| t.fluid.conductivity_w_per_m_k = 0.0),
            Box::new(|t| {
                t.boundaries.push(ThermalPatch {
                    patch: "chip".into(),
                    condition: ThermalCondition::Adiabatic,
                })
            }),
            Box::new(|t| {
                t.boundaries.push(ThermalPatch {
                    patch: "top".into(),
                    condition: ThermalCondition::Adiabatic,
                });
                t.boundaries.push(ThermalPatch {
                    patch: "top".into(),
                    condition: ThermalCondition::Adiabatic,
                });
            }),
            Box::new(|t| {
                t.boundaries.push(ThermalPatch {
                    patch: "top".into(),
                    condition: ThermalCondition::Convection {
                        h_w_per_m2_k: 0.0,
                        ambient_k: 300.0,
                    },
                })
            }),
        ];
        for (i, edit) in edits.iter().enumerate() {
            let mut case = heated_channel(0.01, None);
            if let Some(t) = case.thermal.as_mut() {
                edit(t);
            }
            let dir = scratch.0.join(format!("bad{i}"));
            assert!(Case::CollocatedFlow3d(case).run(&dir).is_err(), "edit {i}");
            assert!(!dir.join("manifest.json").exists());
        }
    }
}

#[test]
fn offline_verification_suite_passes_and_refuses_an_existing_folder() {
    let scratch = Scratch::new();
    let report = cfd_solvers::case::verification::run_suite(&scratch.0.join("verify")).unwrap();
    assert_eq!(report["passed"], true, "{report:#}");
    assert_eq!(report["checks"].as_array().unwrap().len(), 7);
    assert!(scratch.0.join("verify/verification.json").exists());
    assert!(cfd_solvers::case::verification::run_suite(&scratch.0.join("verify")).is_err());
}

#[test]
fn buoyancy_inputs_are_validated_and_drive_flow_in_a_heated_cavity() {
    use cfd_solvers::case::conjugate::*;
    use cfd_solvers::fvm::VelocityBoundaryType;
    let scratch = Scratch::new();
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec {
            dimensions: [10, 10, 1],
            lengths: Vec3::new(1.0, 1.0, 0.1),
            xy_shear: 0.0,
        },
        method: FlowMethod::Piso,
        output: PolyOutput {
            steps: 200,
            output_interval: 200,
        },
        ..Default::default()
    };
    case.config.dt = 0.05;
    case.config.density = 1.0;
    case.config.kinematic_viscosity = 0.0071;
    for (patch, kind) in [
        ("left", VelocityBoundaryType::FixedValue),
        ("right", VelocityBoundaryType::FixedValue),
        ("bottom", VelocityBoundaryType::FixedValue),
        ("top", VelocityBoundaryType::FixedValue),
        ("front", VelocityBoundaryType::Slip),
        ("back", VelocityBoundaryType::Slip),
    ] {
        case.velocity_boundaries.push(VelocityPatch {
            patch: patch.into(),
            kind,
            value: Vec3::ZERO,
        });
    }
    case.thermal = Some(ConjugateHeat {
        fluid: ThermalFluid {
            specific_heat_j_per_kg_k: 1.0,
            conductivity_w_per_m_k: 0.01,
        },
        inflow_temperature_k: 300.0,
        solids: vec![],
        boundaries: vec![
            ThermalPatch {
                patch: "left".into(),
                condition: ThermalCondition::Temperature {
                    temperature_k: 300.5,
                },
            },
            ThermalPatch {
                patch: "right".into(),
                condition: ThermalCondition::Temperature {
                    temperature_k: 299.5,
                },
            },
        ],
        buoyancy: Some(Buoyancy {
            gravity_m_per_s2: [0.0, -1.0, 0.0],
            expansion_per_k: 0.071,
            reference_temperature_k: 300.0,
        }),
    });
    let manifest = Case::CollocatedFlow3d(case.clone())
        .run(&scratch.0.join("cavity"))
        .unwrap();
    assert!(
        manifest.metrics["kinetic_energy_J"] > 1e-8,
        "buoyancy must set the fluid moving"
    );
    assert!(manifest.metrics["max_viscous_diffusion_number"] < 0.1);
    let report: serde_json::Value = serde_json::from_slice(
        &std::fs::read(scratch.0.join("cavity/thermal_report.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(report["model"], "transient_two_way_boussinesq_cht_upwind");
    assert!(
        report["relative_energy_imbalance"].as_f64().unwrap() < 1e-8,
        "{report}"
    );
    let frame = std::fs::read_to_string(scratch.0.join("cavity/flow_000200.vtu")).unwrap();
    assert!(
        frame.contains("Name=\"temperature_k\""),
        "coupled frames carry temperature"
    );
    // Without buoyancy the same cavity stays still: pure conduction.
    let mut still = case.clone();
    if let Some(t) = still.thermal.as_mut() {
        t.buoyancy = None;
    }
    let manifest = Case::CollocatedFlow3d(still)
        .run(&scratch.0.join("still"))
        .unwrap();
    assert!(manifest.metrics["kinetic_energy_J"] < 1e-20);
    for (i, edit) in [
        |b: &mut Buoyancy| b.expansion_per_k = 0.0,
        |b: &mut Buoyancy| b.reference_temperature_k = -1.0,
        |b: &mut Buoyancy| b.gravity_m_per_s2[1] = f64::NAN,
    ]
    .iter()
    .enumerate()
    {
        let mut bad = case.clone();
        if let Some(b) = bad.thermal.as_mut().and_then(|t| t.buoyancy.as_mut()) {
            edit(b);
        }
        let dir = scratch.0.join(format!("bad{i}"));
        assert!(Case::CollocatedFlow3d(bad).run(&dir).is_err());
        assert!(!dir.join("manifest.json").exists());
    }
}

mod surface_bodies {
    use super::*;

    /// Axis-aligned box as 12 outward-wound triangles.
    pub fn box_body(name: &str, lo: [f64; 3], hi: [f64; 3]) -> SurfaceBody {
        let v = |i: usize| {
            [
                if i & 1 == 0 { lo[0] } else { hi[0] },
                if i & 2 == 0 { lo[1] } else { hi[1] },
                if i & 4 == 0 { lo[2] } else { hi[2] },
            ]
        };
        let vertices = (0..8).map(v).collect();
        let quads = [
            [0, 2, 3, 1],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 4, 6, 2],
            [1, 3, 7, 5],
        ];
        let triangles = quads
            .iter()
            .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
            .collect();
        SurfaceBody {
            name: name.into(),
            vertices,
            triangles,
        }
    }

    /// Subdivided icosahedron projected to a sphere.
    pub fn sphere(name: &str, centre: [f64; 3], radius: f64, levels: usize) -> SurfaceBody {
        let t = (1.0 + 5f64.sqrt()) / 2.0;
        let mut vertices: Vec<[f64; 3]> = [
            [-1.0, t, 0.0],
            [1.0, t, 0.0],
            [-1.0, -t, 0.0],
            [1.0, -t, 0.0],
            [0.0, -1.0, t],
            [0.0, 1.0, t],
            [0.0, -1.0, -t],
            [0.0, 1.0, -t],
            [t, 0.0, -1.0],
            [t, 0.0, 1.0],
            [-t, 0.0, -1.0],
            [-t, 0.0, 1.0],
        ]
        .to_vec();
        let mut triangles: Vec<[usize; 3]> = vec![
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ];
        for _ in 0..levels {
            let mut cache = std::collections::HashMap::new();
            let mut mid = |a: usize, b: usize, vertices: &mut Vec<[f64; 3]>| {
                *cache.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    let (p, q) = (vertices[a], vertices[b]);
                    vertices.push([
                        (p[0] + q[0]) / 2.0,
                        (p[1] + q[1]) / 2.0,
                        (p[2] + q[2]) / 2.0,
                    ]);
                    vertices.len() - 1
                })
            };
            let mut next = vec![];
            for [a, b, c] in triangles {
                let (ab, bc, ca) = (
                    mid(a, b, &mut vertices),
                    mid(b, c, &mut vertices),
                    mid(c, a, &mut vertices),
                );
                next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
            }
            triangles = next;
        }
        let vertices = vertices
            .iter()
            .map(|p| {
                let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                [
                    centre[0] + radius * p[0] / n,
                    centre[1] + radius * p[1] / n,
                    centre[2] + radius * p[2] / n,
                ]
            })
            .collect();
        SurfaceBody {
            name: name.into(),
            vertices,
            triangles,
        }
    }

    #[test]
    fn closed_surfaces_are_checked_and_voxelized_by_cell_centres() {
        let spec = PolyMeshSpec {
            dimensions: [10, 10, 10],
            lengths: Vec3::new(1.0, 1.0, 1.0),
            xy_shear: 0.0,
        };
        // Box faces between cell centres: exactly 3 x 4 x 2 cells.
        let b = box_body("block", [0.21, 0.11, 0.42], [0.49, 0.48, 0.61]);
        assert!((b.validate().unwrap() - 0.28 * 0.37 * 0.19).abs() < 1e-12);
        let inside = b.voxelize(spec).unwrap();
        assert_eq!(inside.iter().filter(|v| **v).count(), 3 * 4 * 2);
        // A face passing exactly through cell centres is counted once, by the tie rule.
        let edge = box_body("edge", [0.25, 0.25, 0.25], [0.75, 0.75, 0.75])
            .voxelize(spec)
            .unwrap();
        assert_eq!(edge.iter().filter(|v| **v).count(), 5 * 5 * 5);
        let mut flipped = b.clone();
        for t in &mut flipped.triangles {
            t.swap(1, 2);
        }
        assert!(flipped.validate().is_err(), "inward winding is refused");
        let mut open = b.clone();
        open.triangles.pop();
        assert!(open.validate().is_err(), "an open surface is refused");
        let mut bad = b.clone();
        bad.triangles[0][0] = 99;
        assert!(bad.validate().is_err());
        // The represented sphere volume converges to the true volume under refinement.
        let s = sphere("ball", [0.5, 0.5, 0.5], 0.3, 3);
        let mut errors = vec![];
        for n in [10usize, 20, 30] {
            let spec = PolyMeshSpec {
                dimensions: [n, n, n],
                lengths: Vec3::new(1.0, 1.0, 1.0),
                xy_shear: 0.0,
            };
            let cells = s.voxelize(spec).unwrap().iter().filter(|v| **v).count();
            errors.push(
                (cells as f64 / (n * n * n) as f64 - s.validate().unwrap()).abs()
                    / s.validate().unwrap(),
            );
        }
        assert!(errors[2] < 0.03 && errors[2] < errors[0], "{errors:?}");
    }

    #[test]
    fn a_surface_body_is_a_conducting_powered_solid_in_the_flow() {
        use cfd_solvers::case::conjugate::*;
        let scratch = Scratch::new();
        let mut case = fan_channel(&[(0.0, 0.02), (2e-5, 0.0)]);
        case.bodies
            .push(sphere("ball", [0.06, 0.01, 0.005], 0.0075, 2));
        case.thermal = Some(ConjugateHeat {
            fluid: ThermalFluid {
                specific_heat_j_per_kg_k: 1006.0,
                conductivity_w_per_m_k: 0.026,
            },
            inflow_temperature_k: 298.15,
            solids: vec![ThermalSolid {
                obstacle: "ball".into(),
                conductivity_w_per_m_k: [200.0; 3],
                power_w: 0.005,
                max_temperature_k: None,
                density_kg_per_m3: None,
                specific_heat_j_per_kg_k: None,
            }],
            boundaries: vec![],
            buoyancy: None,
        });
        let manifest = Case::CollocatedFlow3d(case.clone())
            .run(&scratch.0.join("ball"))
            .unwrap();
        assert!(manifest
            .capabilities
            .iter()
            .any(|c| c == "surface_bodies_v1"));
        let report: serde_json::Value = serde_json::from_slice(
            &std::fs::read(scratch.0.join("ball/mesh_report.json")).unwrap(),
        )
        .unwrap();
        let body = &report["surface_bodies"][0];
        assert_eq!(
            body["cells"], 8,
            "x of 55 and 65 mm by y of 5, 8.3, 11.7 and 15 mm lie inside"
        );
        let thermal: serde_json::Value = serde_json::from_slice(
            &std::fs::read(scratch.0.join("ball/thermal_report.json")).unwrap(),
        )
        .unwrap();
        assert!((thermal["advected_heat_out_w"].as_f64().unwrap() - 0.005).abs() < 1e-10);
        // Overlap with a box obstacle is refused.
        let mut overlap = case.clone();
        overlap.obstacles.push(GridObstacle {
            name: "box".into(),
            min: [5, 2, 0],
            max: [7, 4, 1],
        });
        assert!(overlap.build_mesh().is_err());
        // A body that covers no cell centre must be resolved by a finer grid.
        let mut tiny = case;
        tiny.bodies[0] = sphere("ball", [0.0605, 0.0105, 0.0055], 0.0002, 1);
        assert!(tiny
            .build_mesh()
            .unwrap_err()
            .to_string()
            .contains("refine"));
    }
}

#[test]
fn goals_converge_stop_the_run_and_match_the_final_reports() {
    use cfd_solvers::case::conjugate::*;
    use cfd_solvers::case::goals::*;
    let scratch = Scratch::new();
    let mut case = fan_channel(&[(0.0, 0.02), (2e-5, 0.0)]);
    case.output = PolyOutput {
        steps: 2000,
        output_interval: 500,
    };
    case.obstacles.push(GridObstacle {
        name: "chip".into(),
        min: [5, 0, 0],
        max: [7, 2, 1],
    });
    case.thermal = Some(ConjugateHeat {
        fluid: ThermalFluid {
            specific_heat_j_per_kg_k: 1006.0,
            conductivity_w_per_m_k: 0.026,
        },
        inflow_temperature_k: 298.15,
        solids: vec![ThermalSolid {
            obstacle: "chip".into(),
            conductivity_w_per_m_k: [150.0; 3],
            power_w: 0.01,
            max_temperature_k: None,
            density_kg_per_m3: None,
            specific_heat_j_per_kg_k: None,
        }],
        boundaries: vec![ThermalPatch {
            patch: "top".into(),
            condition: ThermalCondition::Convection {
                h_w_per_m2_k: 10.0,
                ambient_k: 298.15,
            },
        }],
        buoyancy: None,
    });
    let goal = |name: &str, kind, target: &str| Goal {
        name: name.into(),
        kind,
        target: target.into(),
    };
    case.goals = vec![
        goal("chip_max", GoalKind::SolidMaxTemperature, "chip"),
        goal("chip_mean", GoalKind::SolidMeanTemperature, "chip"),
        goal("outflow", GoalKind::PatchVolumeFlow, "right"),
        goal("inlet_pressure", GoalKind::PatchMeanPressure, "left"),
        goal("top_heat", GoalKind::PatchHeatOut, "top"),
    ];
    case.goal_control = Some(GoalControl {
        window_steps: 25,
        relative_tolerance: 1e-3,
        stop_when_converged: true,
        min_steps: 100,
    });
    let dir = scratch.0.join("goals");
    let manifest = Case::CollocatedFlow3d(case.clone()).run(&dir).unwrap();
    let goals: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("goals.json")).unwrap()).unwrap();
    assert_eq!(goals["all_converged"], true, "{goals:#}");
    let stopped = goals["stopped_early_at_step"].as_u64().unwrap();
    assert!(stopped >= 100 && stopped < 2000, "stopped at {stopped}");
    assert_eq!(manifest.metrics["goals_converged"], 1.0);
    assert!(
        dir.join(format!("flow_{stopped:06}.vtu")).exists(),
        "the stopping step is the last snapshot"
    );
    // Goals agree with the final reports of the same run.
    let thermal: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("thermal_report.json")).unwrap()).unwrap();
    let chip_max = thermal["components"][0]["max_temperature_k"]
        .as_f64()
        .unwrap();
    assert!((manifest.metrics["goal_chip_max"] - chip_max).abs() < 1e-9);
    let mesh: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("mesh_report.json")).unwrap()).unwrap();
    let right = mesh["patches"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "right")
        .unwrap();
    assert!(
        (manifest.metrics["goal_outflow"]
            - right["outward_volume_flux_m3_per_s"].as_f64().unwrap())
        .abs()
            < 1e-15
    );
    let fan: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("fan_operating_point.json")).unwrap())
            .unwrap();
    assert!(
        (manifest.metrics["goal_inlet_pressure"]
            - fan["system_pressure_rise_pa"].as_f64().unwrap())
        .abs()
            < 1e-6 * fan["system_pressure_rise_pa"].as_f64().unwrap().abs()
    );
    assert!(manifest.metrics["goal_top_heat"] > 0.0);
    // Invalid goals are refused before solving.
    let edits: Vec<Box<dyn Fn(&mut CollocatedFlowCase)>> = vec![
        Box::new(|c| c.goals[0].target = "ghost".into()),
        Box::new(|c| c.goals[4].target = "chip".into()),
        Box::new(|c| c.goals[2].target = "nowhere".into()),
        Box::new(|c| c.goals[1].name = c.goals[0].name.clone()),
        Box::new(|c| c.goals[1].name = "bad name".into()),
        Box::new(|c| c.thermal = None),
        Box::new(|c| {
            c.goal_control = Some(GoalControl {
                window_steps: 0,
                ..GoalControl::default()
            })
        }),
        Box::new(|c| {
            c.goal_control = Some(GoalControl {
                relative_tolerance: 1.5,
                ..GoalControl::default()
            })
        }),
    ];
    for (i, edit) in edits.iter().enumerate() {
        let mut bad = case.clone();
        edit(&mut bad);
        let dir = scratch.0.join(format!("bad{i}"));
        assert!(Case::CollocatedFlow3d(bad).run(&dir).is_err(), "edit {i}");
        assert!(!dir.join("manifest.json").exists());
    }
}
