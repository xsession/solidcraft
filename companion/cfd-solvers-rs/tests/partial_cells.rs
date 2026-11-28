//! Partial cells for conjugate heat (`partial_cells_v1`): a 1 W chip, 4 mm tall, sits on a board
//! in a channel of 10 mm cells. No cell centre lies inside the chip, so as a staircase it has no
//! cell at all: refused unless it is `lumped`, and then only with `thermal.partial_cells`, which
//! places its power in the cells its volume lies in and reads its temperature there.

use cfd_solvers::case::conjugate::ConjugateHeat;
use cfd_solvers::case::*;
use cfd_solvers::core::Vec3;
use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};
use serde_json::json;

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("partial-cells-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        Self(dir)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn block(name: &str, lo: [f64; 3], hi: [f64; 3], lumped: bool) -> SurfaceBody {
    let v = |i: usize| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }];
    let quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]];
    SurfaceBody {
        name: name.into(),
        vertices: (0..8).map(v).collect(),
        triangles: quads.iter().flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]).collect(),
        lumped,
    }
}

fn channel(lumped: bool, partial: bool) -> CollocatedFlowCase {
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec { dimensions: [20, 5, 4], lengths: Vec3::new(0.2, 0.05, 0.04), xy_shear: 0.0, coordinates: None },
        method: FlowMethod::Piso,
        output: PolyOutput { steps: 40, output_interval: 40 },
        bodies: vec![
            block("board", [0.04, 0.01, 0.01], [0.12, 0.03, 0.03], false),
            // On the board's top (y = 30 mm) up to 34 mm: below the next centres at 35 mm.
            block("chip", [0.075, 0.03, 0.015], [0.085, 0.034, 0.025], lumped),
        ],
        ..Default::default()
    };
    case.config.dt = 1e-2;
    case.config.kinematic_viscosity = 1.5e-5;
    case.config.pressure_tolerance = 1e-12;
    case.config.pressure_iterations = 4000;
    case.velocity_boundaries.push(VelocityPatch { patch: "left".into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::new(0.2, 0.0, 0.0) });
    for patch in ["bottom", "top", "front", "back"] {
        case.velocity_boundaries.push(VelocityPatch { patch: patch.into(), kind: VelocityBoundaryType::FixedValue, value: Vec3::ZERO });
    }
    case.velocity_boundaries.push(VelocityPatch { patch: "right".into(), kind: VelocityBoundaryType::ZeroGradient, value: Vec3::ZERO });
    case.pressure_boundaries.push(PressurePatch { patch: "right".into(), kind: PressureBoundaryType::FixedValue, value: 0.0 });
    let heat: ConjugateHeat = serde_json::from_value(json!({
        "fluid": {"specific_heat_j_per_kg_k": 1005.0, "conductivity_w_per_m_k": 0.026},
        "inflow_temperature_k": 300.0,
        "solids": [
            {"obstacle": "board", "conductivity_w_per_m_k": [10.0, 10.0, 10.0], "power_w": 0.0},
            {"obstacle": "chip", "conductivity_w_per_m_k": [150.0, 150.0, 150.0], "power_w": 1.0, "max_temperature_k": 400.0}
        ],
        "partial_cells": partial
    }))
    .unwrap();
    case.thermal = Some(heat);
    case
}

#[test]
fn a_chip_thinner_than_a_cell_is_lumped_with_partial_cells() {
    let scratch = Scratch::new("refused");
    // As a staircase the chip has no cell: refused, unless lumped, and then heat needs partial cells.
    let plain = Case::CollocatedFlow3d(channel(false, false)).run(&scratch.0.join("a")).unwrap_err().to_string();
    assert!(plain.contains("covers no cell centre"), "{plain}");
    let lumped = Case::CollocatedFlow3d(channel(true, false)).run(&scratch.0.join("b")).unwrap_err().to_string();
    assert!(lumped.contains("needs thermal.partial_cells"), "{lumped}");

    let scratch = Scratch::new("partial");
    // A goal on the lumped chip reads the cells it lies in, as its row in the component table
    // does (it read no cells and gave NaN, which failed the run when the goals were written).
    let mut case = channel(true, true);
    case.goals = serde_json::from_value(json!([{"name": "chip_max", "kind": "solid_max_temperature", "target": "chip"},
                                               {"name": "chip_mean", "kind": "solid_mean_temperature", "target": "chip"}])).unwrap();
    case.goal_control = serde_json::from_value(json!({"window_steps": 10, "relative_tolerance": 1e-9, "stop_when_converged": false, "min_steps": 20})).unwrap();
    let manifest = Case::CollocatedFlow3d(case).run(&scratch.0).unwrap();
    assert!(manifest.capabilities.iter().any(|c| c == "partial_cells_v1"));
    // The chip's watt is all there, and the air carries it out.
    assert!((manifest.metrics["thermal_total_power_w"] - 1.0).abs() < 1e-9, "{:?}", manifest.metrics);
    assert!(manifest.metrics["thermal_relative_energy_imbalance"] < 1e-6, "{:?}", manifest.metrics);
    let report: serde_json::Value = serde_json::from_slice(&std::fs::read(scratch.0.join("thermal_report.json")).unwrap()).unwrap();
    let part = |name: &str| report["components"].as_array().unwrap().iter().find(|c| c["name"] == name).cloned().unwrap();
    let (chip, board) = (part("chip"), part("board"));
    let chip_max = chip["max_temperature_k"].as_f64().unwrap();
    // It reads the cells it lies in, which it heats: warmer than the inflow and than the board under it.
    assert!(chip_max > 300.0 && chip_max < 400.0, "{chip}");
    assert!(chip_max >= board["max_temperature_k"].as_f64().unwrap() - 1e-9, "{chip} {board}");
    assert_eq!(chip["status"], "pass");
    let goals: serde_json::Value = serde_json::from_slice(&std::fs::read(scratch.0.join("goals.json")).unwrap()).unwrap();
    let goal = |name: &str| goals["goals"].as_array().unwrap().iter().find(|g| g["name"] == name).unwrap()["final_value"].as_f64().unwrap();
    assert!((goal("chip_max") - chip_max).abs() < 1e-6, "{goals}");
    assert!(goal("chip_mean").is_finite() && goal("chip_mean") > 300.0 && goal("chip_mean") <= chip_max + 1e-9, "{goals}");
}

fn patch_report(dir: &std::path::Path, name: &str) -> serde_json::Value {
    let text = std::fs::read_to_string(dir.join("mesh_report.json")).unwrap();
    let report: serde_json::Value = serde_json::from_str(&text).unwrap();
    report["patches"].as_array().unwrap().iter().find(|p| p["name"] == name).cloned().unwrap_or_else(|| panic!("no patch {name}: {report}"))
}

/// A 2 mm baffle at x = 100 mm across the lower 30 mm of the channel, between the cell centres
/// at 95 and 105 mm: lumped, it holds no cell but closes the 3 × 4 grid faces it crosses, a
/// wall on both sides; the flow goes over it and all of it still leaves.
#[test]
fn a_lumped_plate_is_a_wall_on_the_grid_faces_it_crosses() {
    let mut case = channel(true, false);
    case.thermal = None;
    case.bodies = vec![block("plate", [0.1, 0.0, 0.0], [0.102, 0.03, 0.04], true)];
    let scratch = Scratch::new("plate");
    Case::CollocatedFlow3d(case).run(&scratch.0).unwrap();
    let plate = patch_report(&scratch.0, "plate");
    assert!((plate["area_m2"].as_f64().unwrap() - 2.0 * 12.0 * 1e-4).abs() < 1e-12, "{plate}");
    assert!(plate["outward_volume_flux_m3_per_s"].as_f64().unwrap().abs() < 1e-15, "{plate}");
    let inflow = 0.2 * 0.05 * 0.04;
    let out = patch_report(&scratch.0, "right")["outward_volume_flux_m3_per_s"].as_f64().unwrap();
    assert!((out - inflow).abs() < 1e-5 * inflow, "{out} {inflow}");
}
