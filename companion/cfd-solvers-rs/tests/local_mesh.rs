//! Rectilinear (locally refined) grids for the collocated flow case (`local_mesh_v1`): the
//! mesh builder agrees with the uniform one, and grid parts, porous zones, surface bodies and
//! conjugate heat follow the graded cells while conserving mass and energy.
use cfd_solvers::case::conjugate::*;
use cfd_solvers::fvm::{make_cartesian_hexa_mesh, make_rectilinear_hexa_mesh, PressureBoundaryType, VelocityBoundaryType};
use cfd_solvers::{case::*, core::Vec3};
use std::sync::atomic::{AtomicUsize, Ordering};
static NEXT: AtomicUsize = AtomicUsize::new(0);

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!("cfd-local-mesh-{}-{}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn even(n: usize, length: f64) -> Vec<f64> {
    (0..=n).map(|i| length * i as f64 / n as f64).collect()
}

/// x nodes of a 0.12 m channel of 0.01 m cells, with the cells over 0.04..0.08 halved.
fn graded_x() -> Vec<f64> {
    let mut x: Vec<f64> = (0..=4).map(|i| 0.01 * i as f64).collect();
    x.extend((1..=8).map(|i| 0.04 + 0.005 * i as f64));
    x.extend((1..=4).map(|i| 0.08 + 0.01 * i as f64));
    x
}

/// Old uniform node index (0..=12) to the graded grid's node index.
fn remap(i: usize) -> usize {
    match i {
        0..=4 => i,
        5..=8 => 4 + 2 * (i - 4),
        _ => i + 4,
    }
}

fn channel(graded: bool) -> CollocatedFlowCase {
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec {
            dimensions: [if graded { 16 } else { 12 }, 6, 1],
            lengths: Vec3::new(0.12, 0.02, 0.01),
            xy_shear: 0.0,
            coordinates: graded.then(|| MeshCoordinates { x: graded_x(), y: even(6, 0.02), z: even(1, 0.01) }),
        },
        method: FlowMethod::Piso,
        output: PolyOutput { steps: 300, output_interval: 300 },
        ..Default::default()
    };
    case.config.dt = 0.01;
    case.config.density = 1.2;
    case.config.kinematic_viscosity = 1.5e-5;
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
    // A chip on a spreader, given in uniform cells and remapped onto the graded grid.
    let x = |lo: usize, hi: usize| if graded { (remap(lo), remap(hi)) } else { (lo, hi) };
    let (s0, s1) = x(4, 8);
    let (c0, c1) = x(5, 7);
    case.obstacles.push(GridObstacle { name: "spreader".into(), min: [s0, 0, 0], max: [s1, 1, 1] });
    case.obstacles.push(GridObstacle { name: "chip".into(), min: [c0, 1, 0], max: [c1, 2, 1] });
    let solid = |name: &str, k: f64, power: f64| ThermalSolid {
        obstacle: name.into(),
        conductivity_w_per_m_k: [k; 3],
        power_w: power,
        max_temperature_k: None,
        density_kg_per_m3: None,
        specific_heat_j_per_kg_k: None,
        contact_resistance_m2_k_per_w: None,
        junction_to_case_k_per_w: None,
        heat_pipe: None,
    };
    case.thermal = Some(ConjugateHeat {
        fluid: ThermalFluid { specific_heat_j_per_kg_k: 1006.0, conductivity_w_per_m_k: 0.026 },
        inflow_temperature_k: 298.15,
        solids: vec![solid("spreader", 200.0, 0.0), solid("chip", 100.0, 0.01)],
        boundaries: vec![],
        buoyancy: None,
        thermoelectric_coolers: vec![],
        face_sources: vec![],
        partial_cells: false,
    });
    case
}

fn read(dir: &std::path::Path, name: &str) -> serde_json::Value {
    serde_json::from_slice(&std::fs::read(dir.join(name)).unwrap()).unwrap()
}

#[test]
fn rectilinear_meshes_reproduce_the_uniform_mesh_and_grade_cells() {
    let uniform = make_cartesian_hexa_mesh(5, 4, 3, 1.0, 0.8, 0.3).unwrap();
    let nodes = make_rectilinear_hexa_mesh(&even(5, 1.0), &even(4, 0.8), &even(3, 0.3)).unwrap();
    assert_eq!(uniform.cells().len(), nodes.cells().len());
    assert_eq!(uniform.faces().len(), nodes.faces().len());
    for (a, b) in uniform.cells().iter().zip(nodes.cells()) {
        assert!((a.center - b.center).norm() < 1e-14 && (a.volume - b.volume).abs() < 1e-15);
    }
    for (a, b) in uniform.faces().iter().zip(nodes.faces()) {
        assert_eq!((a.owner, a.neighbour, a.patch), (b.owner, b.neighbour, b.patch));
        assert!((a.center - b.center).norm() < 1e-14 && (a.area - b.area).norm() < 1e-15);
    }
    // Graded: volumes add up to the box, and each cell is closed (its face areas sum to zero).
    let graded = make_rectilinear_hexa_mesh(&[0.0, 0.1, 0.15, 0.175, 1.0], &[0.0, 0.5, 1.0], &[0.0, 0.2]).unwrap();
    let total: f64 = graded.cells().iter().map(|c| c.volume).sum();
    assert!((total - 0.2).abs() < 1e-15);
    let mut closure = vec![Vec3::ZERO; graded.cells().len()];
    for f in graded.faces() {
        closure[f.owner] = closure[f.owner] + f.area;
        if let Some(n) = f.neighbour {
            closure[n] = closure[n] - f.area;
        }
    }
    assert!(closure.iter().all(|v| v.norm() < 1e-15));
    // Bad node lists are refused.
    assert!(make_rectilinear_hexa_mesh(&[0.0, 0.2, 0.1], &[0.0, 1.0], &[0.0, 1.0]).is_err());
    assert!(make_rectilinear_hexa_mesh(&[0.1, 0.2], &[0.0, 1.0], &[0.0, 1.0]).is_err());
    assert!(make_rectilinear_hexa_mesh(&[0.0], &[0.0, 1.0], &[0.0, 1.0]).is_err());
}

#[test]
fn a_locally_refined_channel_conserves_mass_and_energy_and_agrees_with_the_uniform_grid() {
    let scratch = Scratch::new();
    let mut results = vec![];
    for graded in [false, true] {
        let dir = scratch.0.join(if graded { "graded" } else { "uniform" });
        let manifest = Case::CollocatedFlow3d(channel(graded)).run(&dir).unwrap();
        assert!(manifest.capabilities.iter().any(|c| c == "local_mesh_v1"));
        let report = read(&dir, "thermal_report.json");
        assert!(report["relative_energy_imbalance"].as_f64().unwrap() < 1e-6, "{report:#}");
        assert!(manifest.metrics["relative_boundary_flux_imbalance"] < 1e-6, "{:?}", manifest.metrics);
        let chip = report["components"].as_array().unwrap().iter().find(|c| c["name"] == "chip").unwrap()["mean_temperature_k"].as_f64().unwrap();
        // The VTU carries the graded nodes.
        let vtu = std::fs::read_to_string(dir.join("thermal.vtu")).unwrap();
        results.push((chip, vtu));
    }
    let (uniform, graded) = (&results[0], &results[1]);
    let rise = |t: f64| t - 298.15;
    eprintln!("chip rise uniform {:.6} K graded {:.6} K", rise(uniform.0), rise(graded.0));
    assert!(rise(uniform.0) > 0.0 && rise(graded.0) > 0.0);
    assert!((rise(graded.0) - rise(uniform.0)).abs() < 0.05 * rise(uniform.0), "uniform {} graded {}", uniform.0, graded.0);
    assert!(graded.1.contains("NumberOfPoints=\"238\""), "17 × 7 × 2 graded nodes");
}

#[test]
fn porous_zones_and_surface_bodies_follow_graded_cells() {
    let scratch = Scratch::new();
    let mut case = channel(true);
    case.thermal = None;
    case.obstacles.clear();
    // A porous zone over the refined cells 4..12 (x 0.04..0.08): 8 × 6 cells.
    case.porous_blocks.push(PorousBlock {
        name: "filter".into(),
        min: [4, 0, 0],
        max: [12, 6, 1],
        viscous_resistance_per_m2: 1e6,
        inertial_resistance_per_m: 0.0,
    });
    let dir = scratch.0.join("porous");
    Case::CollocatedFlow3d(case.clone()).run(&dir).unwrap();
    let zones = read(&dir, "porous_zones.json");
    assert_eq!(zones["zones"][0]["cells"], 48);
    assert!((zones["zones"][0]["volume_m3"].as_f64().unwrap() - 0.04 * 0.02 * 0.01).abs() < 1e-15);
    // A surface box over x 0.042..0.058 and y 0..0.0067 in the refined region.
    let body = SurfaceBody {
        name: "block".into(),
        vertices: (0..8)
            .map(|i: usize| {
                [
                    if i & 1 == 0 { 0.042 } else { 0.058 },
                    if i & 2 == 0 { 0.0 } else { 0.0067 },
                    if i & 4 == 0 { 0.0 } else { 0.01 },
                ]
            })
            .collect(),
        triangles: [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]]
            .iter()
            .flat_map(|q: &[usize; 4]| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
            .collect(),
        lumped: false,
    };
    let inside = body.voxelize(&case.mesh).unwrap();
    // Fine x centres in range: 0.0425, 0.0475, 0.0525, 0.0575 → 4; y centres below 0.0067: 2.
    assert_eq!(inside.iter().filter(|v| **v).count(), 4 * 2);
    case.porous_blocks.clear();
    case.bodies.push(body);
    let dir = scratch.0.join("body");
    let manifest = Case::CollocatedFlow3d(case).run(&dir).unwrap();
    assert!(manifest.metrics["relative_boundary_flux_imbalance"] < 1e-6);
}

#[test]
fn coordinates_must_agree_with_the_grid() {
    let scratch = Scratch::new();
    let mut wrong_count = channel(true);
    wrong_count.mesh.dimensions[0] = 12;
    assert!(Case::CollocatedFlow3d(wrong_count).run(&scratch.0.join("a")).is_err());
    let mut wrong_length = channel(true);
    wrong_length.mesh.lengths = Vec3::new(0.13, 0.02, 0.01);
    assert!(Case::CollocatedFlow3d(wrong_length).run(&scratch.0.join("b")).is_err());
    let mut sheared = channel(true);
    sheared.mesh.xy_shear = 0.1;
    assert!(Case::CollocatedFlow3d(sheared).run(&scratch.0.join("c")).is_err());
    let mut decreasing = channel(true);
    if let Some(c) = decreasing.mesh.coordinates.as_mut() {
        c.x.swap(2, 3);
    }
    assert!(Case::CollocatedFlow3d(decreasing).run(&scratch.0.join("d")).is_err());
    // The JSON form round-trips, and unknown keys are refused.
    let json = serde_json::to_value(&channel(true).mesh).unwrap();
    assert_eq!(json["coordinates"]["x"].as_array().unwrap().len(), 17);
    let mut extra = json.clone();
    extra["coordinates"]["w"] = serde_json::json!([0.0, 1.0]);
    assert!(serde_json::from_value::<PolyMeshSpec>(extra).is_err());
    assert!(serde_json::to_value(&channel(false).mesh).unwrap().get("coordinates").is_none());
}
