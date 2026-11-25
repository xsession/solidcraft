//! Cut cells (`cut_cells_v1`): cells the surface crosses keep their fluid part.
//!
//! A sphere in a box: the fluid the cells hold matches the exact fluid volume (the box minus the
//! triangle surface's own enclosed volume) far better than the staircase of whole cells. A block
//! whose face lies 0.7 of a cell into the grid leaves the exact gap open at the constriction, where
//! the staircase is off by up to half a cell. Every cell stays closed, small cells are merged,
//! and a flow through the constriction runs and balances.

use cfd_solvers::case::*;
use cfd_solvers::core::Vec3;
use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("cut-cells-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        Self(dir)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A closed outward triangle sphere: an octahedron split `levels` times, pushed to the radius.
fn sphere(centre: [f64; 3], radius: f64, levels: usize) -> SurfaceBody {
    let mut vertices: Vec<[f64; 3]> = vec![[1., 0., 0.], [-1., 0., 0.], [0., 1., 0.], [0., -1., 0.], [0., 0., 1.], [0., 0., -1.]];
    let mut triangles: Vec<[usize; 3]> =
        vec![[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4], [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]];
    for _ in 0..levels {
        let mut midpoints = std::collections::HashMap::new();
        let mut next = Vec::new();
        let mut mid = |a: usize, b: usize, vertices: &mut Vec<[f64; 3]>| -> usize {
            let key = (a.min(b), a.max(b));
            *midpoints.entry(key).or_insert_with(|| {
                let p: [f64; 3] = std::array::from_fn(|k| 0.5 * (vertices[a][k] + vertices[b][k]));
                let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                vertices.push(p.map(|v| v / n));
                vertices.len() - 1
            })
        };
        for [a, b, c] in triangles {
            let (ab, bc, ca) = (mid(a, b, &mut vertices), mid(b, c, &mut vertices), mid(c, a, &mut vertices));
            next.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
        }
        triangles = next;
    }
    let vertices = vertices.into_iter().map(|p| std::array::from_fn(|k| centre[k] + radius * p[k])).collect();
    SurfaceBody { name: "ball".into(), vertices, triangles, lumped: false }
}

/// A box from `lo` to `hi` as a closed outward triangle surface.
fn block(name: &str, lo: [f64; 3], hi: [f64; 3]) -> SurfaceBody {
    let v = |i: usize| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }];
    let quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]];
    SurfaceBody {
        name: name.into(),
        vertices: (0..8).map(v).collect(),
        triangles: quads.iter().flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]).collect(),
        lumped: false,
    }
}

fn case(dimensions: [usize; 3], lengths: Vec3, bodies: Vec<SurfaceBody>, cut: bool) -> CollocatedFlowCase {
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec { dimensions, lengths, xy_shear: 0.0, coordinates: None },
        method: FlowMethod::Piso,
        output: PolyOutput { steps: 40, output_interval: 40 },
        bodies,
        cut_cells: cut,
        ..Default::default()
    };
    case.config.dt = 2e-3;
    case.config.kinematic_viscosity = 1.5e-5;
    case.config.pressure_tolerance = 1e-12;
    case.config.pressure_iterations = 4000;
    case
}

fn fluid_volume(case: &CollocatedFlowCase) -> f64 {
    case.build_mesh().unwrap().0.cells().iter().map(|c| c.volume).sum()
}

#[test]
fn a_sphere_keeps_its_volume_far_better_than_the_staircase() {
    let ball = sphere([0.05, 0.05, 0.05], 0.031, 4);
    let exact = 0.1f64.powi(3) - ball.validate().unwrap();
    let lengths = Vec3::new(0.1, 0.1, 0.1);
    let staircase = fluid_volume(&case([10, 10, 10], lengths, vec![ball.clone()], false));
    let cut = fluid_volume(&case([10, 10, 10], lengths, vec![ball], true));
    let (stair_error, cut_error) = ((staircase - exact).abs() / exact, (cut - exact).abs() / exact);
    // On 10 mm cells the staircase is off by percents; the cut cells by the sub-lattice's
    // 1.25 mm layer at most, which mostly cancels round the surface.
    assert!(cut_error < 2e-3, "cut cells: {cut_error}");
    assert!(cut_error * 5.0 < stair_error, "cut {cut_error} against staircase {stair_error}");
}

#[test]
fn every_cell_is_closed_and_no_small_cell_is_left() {
    let ball = sphere([0.05, 0.05, 0.05], 0.031, 4);
    let lengths = Vec3::new(0.1, 0.1, 0.1);
    let (mesh, export) = case([10, 10, 10], lengths, vec![ball], true).build_mesh().unwrap();
    assert_eq!(export.cells.len(), mesh.cell_count());
    let hex = 1e-6; // one 10 mm cell's box
    for (ci, cell) in mesh.cells().iter().enumerate() {
        // Small cells were merged into a neighbour: what is left holds at least the merge
        // threshold of one cell's box.
        assert!(cell.volume >= 0.3 * hex * (1.0 - 1e-9), "cell {ci} holds only {} m³", cell.volume);
        let mut closure = Vec3::ZERO;
        let mut total = 0.0;
        for &f in &mesh.cell_faces()[ci] {
            let a = mesh.faces()[f].outward_area(ci);
            closure += a;
            total += a.norm();
        }
        assert!(closure.norm() <= 1e-12 * total, "cell {ci} is not closed: {closure:?}");
    }
}

#[test]
fn a_constriction_keeps_its_exact_gap_and_the_flow_balances() {
    // Channel 0.2 × 0.05 × 0.01 m on 10 mm cells; a block from y = 0 to y = 0.027 (2.7 cells)
    // between x = 0.08 and 0.12 leaves a 0.023 m gap. The staircase opens 0.02 m (3 cells solid).
    let lengths = Vec3::new(0.2, 0.05, 0.01);
    // The block ends on the domain's faces, so its whole volume is in the domain.
    let wall = block("wall", [0.08, 0.0, 0.0], [0.12, 0.027, 0.01]);
    let mut flow = case([20, 5, 1], lengths, vec![wall], true);
    let speed = 0.02;
    for patch in ["bottom", "top", "front", "back"] {
        flow.velocity_boundaries.push(VelocityPatch { patch: patch.into(), kind: VelocityBoundaryType::Slip, value: Vec3::ZERO });
    }
    flow.velocity_boundaries.push(VelocityPatch { patch: "left".into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::new(speed, 0.0, 0.0) });
    flow.velocity_boundaries.push(VelocityPatch { patch: "right".into(), kind: VelocityBoundaryType::ZeroGradient, value: Vec3::ZERO });
    flow.pressure_boundaries.push(PressurePatch { patch: "right".into(), kind: PressureBoundaryType::FixedValue, value: 0.0 });
    let (mesh, _) = flow.build_mesh().unwrap();
    // The open area of the grid faces at x = 0.1 m, in the middle of the constriction.
    let open: f64 = mesh.faces().iter().filter(|f| (f.center.x - 0.1).abs() < 1e-12 && f.area.x.abs() > 0.0).map(|f| f.area.x.abs()).sum();
    let gap = (0.05 - 0.027) * 0.01;
    assert!((open - gap).abs() <= 0.01 / 8.0 * 0.01, "open {open} m², exact gap {gap} m²");
    let scratch = Scratch::new("constriction");
    let manifest = Case::CollocatedFlow3d(flow).run(&scratch.0).unwrap();
    assert!(manifest.capabilities.iter().any(|c| c == "cut_cells_v1"));
    let inflow = speed * 0.05 * 0.01;
    let out = manifest.metrics["boundary_outflow_m3_per_s"];
    assert!((out - inflow).abs() < 1e-5 * inflow, "outflow {out} against inflow {inflow}");
    let report: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(scratch.0.join("mesh_report.json")).unwrap()).unwrap();
    // The sampling moves a cut face by at most half a sub-cell (10 mm / 8 / 2) over its area: the
    // block's 40 × 10 mm face at y = 0.027 m against 8.92e-5 m³ of fluid, 0.28 %.
    let error = report["cut_cells"]["relative_fluid_volume_error"].as_f64().unwrap();
    let bound = 0.5 * 0.01 / 8.0 * 0.04 * 0.01 / (0.2 * 0.05 * 0.01 - 0.04 * 0.027 * 0.01);
    assert!(error.abs() <= bound, "{error} against {bound}");
}

#[test]
fn cut_cells_refuse_what_they_do_not_carry_yet() {
    let lengths = Vec3::new(0.1, 0.1, 0.1);
    let mut porous = case([10, 10, 10], lengths, vec![sphere([0.05, 0.05, 0.05], 0.03, 2)], true);
    porous.porous_blocks = vec![serde_json::from_value(serde_json::json!({
        "name": "filter", "min": [0, 0, 0], "max": [1, 1, 1], "viscous_resistance_per_m2": 1e6, "inertial_resistance_per_m": 10.0
    }))
    .unwrap()];
    let scratch = Scratch::new("refused");
    let why = Case::CollocatedFlow3d(porous).run(&scratch.0).unwrap_err().to_string();
    assert!(why.contains("cut cells are not supported with porous blocks"), "{why}");
}

#[test]
fn a_condition_on_a_cut_face_blows_its_true_area() {
    // The block's +x face (triangles 2 and 3) lies at x = 0.077 m, 0.7 of a 10 mm cell: cut.
    let lengths = Vec3::new(0.2, 0.04, 0.04);
    let mut flow = case([20, 4, 4], lengths, vec![block("block", [0.04, 0.01, 0.01], [0.077, 0.03, 0.03])], true);
    flow.face_conditions = vec![BodyFaceCondition { name: "jet".into(), body: "block".into(), triangles: vec![2, 3] }];
    let speed = 0.05;
    for patch in ["left", "bottom", "top", "front", "back"] {
        flow.velocity_boundaries.push(VelocityPatch { patch: patch.into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::ZERO });
    }
    flow.velocity_boundaries.push(VelocityPatch { patch: "jet".into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::new(speed, 0.0, 0.0) });
    flow.pressure_boundaries.push(PressurePatch { patch: "jet".into(), kind: PressureBoundaryType::ZeroGradient, value: 0.0 });
    flow.velocity_boundaries.push(VelocityPatch { patch: "right".into(), kind: VelocityBoundaryType::ZeroGradient, value: Vec3::ZERO });
    flow.pressure_boundaries.push(PressurePatch { patch: "right".into(), kind: PressureBoundaryType::FixedValue, value: 0.0 });
    let scratch = Scratch::new("jet");
    let manifest = Case::CollocatedFlow3d(flow).run(&scratch.0).unwrap();
    let report: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(scratch.0.join("mesh_report.json")).unwrap()).unwrap();
    let jet = report["patches"].as_array().unwrap().iter().find(|p| p["name"] == "jet").cloned().unwrap();
    // The jet blows through the face's whole 0.02 × 0.02 m: its inflow is the speed times the
    // face's area projected on the flow, to the sampling (cells along the face's edges have one
    // wall face that also holds a bit of the block's sides, so |area| alone is larger).
    let inflow = -jet["outward_volume_flux_m3_per_s"].as_f64().unwrap();
    assert!((inflow - speed * 4e-4).abs() < 0.02 * speed * 4e-4, "jet inflow {inflow}");
    let out = manifest.metrics["boundary_outflow_m3_per_s"];
    assert!((out - inflow).abs() < 1e-5 * inflow, "outflow {out} against the jet's {inflow}");
}

#[test]
fn a_wall_thinner_than_a_cell_stays_a_wall() {
    // A 4 mm plate across a 10 mm-cell channel covers no cell centre: the staircase has no wall
    // at all. Cut cells keep it, and the flow must go round its 30 mm height.
    let lengths = Vec3::new(0.2, 0.05, 0.01);
    // From x = 90.5 to 94.5 mm: between the cell centres at 85 and 95 mm.
    let plate = block("plate", [0.0905, 0.0, 0.0], [0.0945, 0.03, 0.01]);
    let staircase = case([20, 5, 1], lengths, vec![plate.clone()], false);
    assert!(staircase.build_mesh().unwrap_err().to_string().contains("covers no cell centre"));
    let mut flow = case([20, 5, 1], lengths, vec![plate], true);
    let speed = 0.02;
    for patch in ["bottom", "top", "front", "back"] {
        flow.velocity_boundaries.push(VelocityPatch { patch: patch.into(), kind: VelocityBoundaryType::Slip, value: Vec3::ZERO });
    }
    flow.velocity_boundaries.push(VelocityPatch { patch: "left".into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::new(speed, 0.0, 0.0) });
    flow.velocity_boundaries.push(VelocityPatch { patch: "right".into(), kind: VelocityBoundaryType::ZeroGradient, value: Vec3::ZERO });
    flow.pressure_boundaries.push(PressurePatch { patch: "right".into(), kind: PressureBoundaryType::FixedValue, value: 0.0 });
    let (mesh, _) = flow.build_mesh().unwrap();
    // The plate's wall faces carry its two 30 × 10 mm sides (and its 4 × 10 mm top), as area
    // pointing along x into the plate from both sides: |x| components sum to 2 × 3e-4 m².
    let plate_patch = mesh.patches().iter().position(|p| p.name == "plate").unwrap();
    let sides: f64 = mesh.faces().iter().filter(|f| f.patch == Some(plate_patch)).map(|f| f.area.x.abs()).sum();
    assert!((sides - 6e-4).abs() <= 0.05 * 6e-4, "plate sides {sides} m²");
    let scratch = Scratch::new("thin");
    let manifest = Case::CollocatedFlow3d(flow).run(&scratch.0).unwrap();
    let inflow = speed * 0.05 * 0.01;
    let out = manifest.metrics["boundary_outflow_m3_per_s"];
    assert!((out - inflow).abs() < 1e-5 * inflow, "outflow {out} against inflow {inflow}");
}
