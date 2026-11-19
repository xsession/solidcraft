//! Transferred boundary conditions (`transferred_bc_v1`): a child domain cut out of a parent
//! run takes the parent's velocity on its inlet face by face, and reproduces the parent's
//! flow far better than a uniform inlet does.
use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};
use cfd_solvers::{case::*, core::Vec3};
use std::sync::atomic::{AtomicUsize, Ordering};
static NEXT: AtomicUsize = AtomicUsize::new(0);

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!("cfd-transferred-{}-{}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

const H: f64 = 0.01;
const NY: usize = 6;

/// A 2-D channel of `nx` cells along x (0.01 m cells), walls top and bottom, a block on the
/// floor at cells `block`, outlet at the right. The inlet is a uniform 0.02 m/s unless
/// `transferred` gives a per-face profile.
fn channel(nx: usize, block: (usize, usize), transferred: Option<Vec<Vec3>>) -> CollocatedFlowCase {
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec { dimensions: [nx, NY, 1], lengths: Vec3::new(H * nx as f64, H * NY as f64, H), xy_shear: 0.0, coordinates: None },
        method: FlowMethod::Piso,
        output: PolyOutput { steps: 400, output_interval: 400 },
        ..Default::default()
    };
    case.config.dt = 0.02;
    case.config.density = 1.2;
    case.config.kinematic_viscosity = 1e-4;
    case.config.pressure_tolerance = 1e-12;
    case.config.pressure_iterations = 4000;
    for (name, kind, value) in [
        ("left", VelocityBoundaryType::FixedValue, Vec3::new(0.02, 0.0, 0.0)),
        ("right", VelocityBoundaryType::ZeroGradient, Vec3::ZERO),
        ("bottom", VelocityBoundaryType::FixedValue, Vec3::ZERO),
        ("top", VelocityBoundaryType::FixedValue, Vec3::ZERO),
        ("front", VelocityBoundaryType::Slip, Vec3::ZERO),
        ("back", VelocityBoundaryType::Slip, Vec3::ZERO),
    ] {
        case.velocity_boundaries.push(VelocityPatch { patch: name.into(), kind, value });
    }
    case.pressure_boundaries.push(PressurePatch { patch: "right".into(), kind: PressureBoundaryType::FixedValue, value: 0.0 });
    case.obstacles.push(GridObstacle { name: "block".into(), min: [block.0, 0, 0], max: [block.1, 3, 1] });
    if let Some(profile) = transferred {
        case.transferred.push(TransferredPatch { patch: "left".into(), velocity: Some(profile), pressure: None });
    }
    case
}

/// Cell velocities of the last VTU frame, by background cell (None inside the block).
fn velocities(dir: &std::path::Path, nx: usize, block: (usize, usize)) -> Vec<Option<Vec3>> {
    let mut frames: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().to_string()))
        .filter(|n| n.ends_with(".vtu") && !n.starts_with("thermal"))
        .collect();
    frames.sort();
    let text = std::fs::read_to_string(dir.join(frames.last().unwrap())).unwrap();
    let start = text.find("Name=\"velocity_m_per_s\"").unwrap();
    let body = &text[start..];
    let body = &body[body.find('>').unwrap() + 1..body.find("</DataArray>").unwrap()];
    let values: Vec<f64> = body.split_whitespace().map(|v| v.parse().unwrap()).collect();
    let mut cells = values.chunks(3).map(|c| Vec3::new(c[0], c[1], c[2]));
    (0..nx * NY)
        .map(|c| {
            let (i, j) = (c % nx, c / nx);
            let solid = (block.0..block.1).contains(&i) && j < 3;
            if solid { None } else { cells.next() }
        })
        .collect()
}

#[test]
fn a_child_domain_fed_by_its_parent_reproduces_the_parent_flow() {
    let scratch = Scratch::new();
    // Parent: 16 cells, block at x cells 9..11. Child: the parent's cells 4..16.
    let (parent_nx, parent_block) = (16, (9, 11));
    let (child_nx, child_block, offset) = (12, (5, 7), 4);
    let dir = scratch.0.join("parent");
    let manifest = Case::CollocatedFlow3d(channel(parent_nx, parent_block, None)).run(&dir).unwrap();
    assert!(manifest.capabilities.iter().any(|c| c == "transferred_bc_v1"));
    let parent = velocities(&dir, parent_nx, parent_block);
    // The child's inlet face lies between parent cells 3 and 4: linear interpolation.
    let at = |i: usize, j: usize| parent[j * parent_nx + i].unwrap();
    let profile: Vec<Vec3> = (0..NY).map(|j| (at(offset - 1, j) + at(offset, j)) * 0.5).collect();
    let inflow: f64 = profile.iter().map(|v| v.x * H * H).sum();
    assert!(profile[0].x < profile[2].x, "the parent's boundary layer is in the profile");

    let run_child = |transferred: Option<Vec<Vec3>>, name: &str| {
        let dir = scratch.0.join(name);
        Case::CollocatedFlow3d(channel(child_nx, child_block, transferred)).run(&dir).unwrap();
        let report: serde_json::Value = serde_json::from_slice(&std::fs::read(dir.join("mesh_report.json")).unwrap()).unwrap();
        (velocities(&dir, child_nx, child_block), report)
    };
    let (fed, report) = run_child(Some(profile.clone()), "fed");
    let (uniform, _) = run_child(None, "uniform");
    // The inlet carries exactly the transferred flow, and it leaves at the outlet.
    let left = report["patches"].as_array().unwrap().iter().find(|p| p["name"] == "left").unwrap();
    assert!((left["outward_volume_flux_m3_per_s"].as_f64().unwrap() + inflow).abs() < 1e-15, "{left}");
    eprintln!("imbalance {}", report["relative_boundary_flux_imbalance"]);
    assert!(report["relative_boundary_flux_imbalance"].as_f64().unwrap() < 1e-6);
    // Away from the inlet, the fed child matches the parent; the uniform inlet does not.
    let error = |child: &[Option<Vec3>]| {
        let mut worst: f64 = 0.0;
        for j in 0..NY {
            for i in 0..child_nx {
                if let (Some(c), Some(p)) = (child[j * child_nx + i], parent[j * parent_nx + i + offset]) {
                    worst = worst.max((c - p).norm());
                }
            }
        }
        worst / 0.02
    };
    let (fed_error, uniform_error) = (error(&fed), error(&uniform));
    eprintln!("max velocity error vs parent: transferred {fed_error:.4}, uniform {uniform_error:.4} (of the inflow speed)");
    assert!(fed_error < 0.05, "{fed_error}");
    assert!(fed_error < 0.25 * uniform_error, "{fed_error} vs {uniform_error}");
}

#[test]
fn transferred_profiles_are_checked() {
    let scratch = Scratch::new();
    let good = vec![Vec3::new(0.02, 0.0, 0.0); NY];
    let run = |case: CollocatedFlowCase, name: &str| Case::CollocatedFlow3d(case).run(&scratch.0.join(name));
    // Wrong count, non-finite values, a patch without a fixed-value condition, an unknown
    // patch, a duplicate patch and an empty entry are all refused.
    assert!(run(channel(8, (5, 6), Some(vec![Vec3::ZERO; NY - 1])), "count").is_err());
    assert!(run(channel(8, (5, 6), Some(vec![Vec3::new(f64::NAN, 0.0, 0.0); NY])), "nan").is_err());
    let mut outlet = channel(8, (5, 6), None);
    outlet.transferred.push(TransferredPatch { patch: "right".into(), velocity: Some(good.clone()), pressure: None });
    assert!(run(outlet, "outlet").is_err());
    let mut unknown = channel(8, (5, 6), None);
    unknown.transferred.push(TransferredPatch { patch: "block".into(), velocity: Some(good.clone()), pressure: None });
    assert!(run(unknown, "unknown").is_err());
    let mut twice = channel(8, (5, 6), Some(good.clone()));
    twice.transferred.push(TransferredPatch { patch: "left".into(), velocity: Some(good.clone()), pressure: None });
    assert!(run(twice, "twice").is_err());
    let mut empty = channel(8, (5, 6), None);
    empty.transferred.push(TransferredPatch { patch: "left".into(), velocity: None, pressure: None });
    assert!(run(empty, "empty").is_err());
    // A transferred outlet pressure is accepted on a fixed-value pressure patch.
    let mut pressure = channel(8, (5, 6), None);
    pressure.transferred.push(TransferredPatch { patch: "right".into(), velocity: None, pressure: Some(vec![0.0; NY]) });
    assert!(run(pressure, "pressure").is_ok());
    // JSON round trip, and unknown keys are refused.
    let json = serde_json::to_value(channel(8, (5, 6), Some(good))).unwrap();
    assert_eq!(json["transferred"][0]["velocity"].as_array().unwrap().len(), NY);
    let mut extra = json.clone();
    extra["transferred"][0]["temperature"] = serde_json::json!([]);
    assert!(serde_json::from_value::<CollocatedFlowCase>(extra).is_err());
}
