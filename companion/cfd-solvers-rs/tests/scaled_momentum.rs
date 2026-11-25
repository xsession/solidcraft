//! Scaled momentum solves (`scaled_momentum_v1`): laminar air in an 8 × 8 mm square duct, 0.8 mm
//! cells, 0.05 m/s. The momentum right-hand side is about cell volume × u / dt ≈ 5e-10 × 50, far
//! below the linear solvers' absolute floor of 1: without the flag the momentum is not solved and
//! the flow stays at its start; with it (and the pressure's own scaling,
//! `piso_neighbour_update`) the developed pressure gradient approaches the exact laminar value
//! dp/dx = (f·Re)·ν·U / (2·a²), f·Re = 56.91 for a square duct (Shah & London 1978).

use cfd_solvers::case::*;
use serde_json::json;

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("scaled-momentum-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        Self(dir)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

const NU: f64 = 1.5e-5;
const U: f64 = 0.05;
const A: f64 = 0.008;

/// The kinematic pressure gradient between x = 24 and 40 mm of a 48 mm duct.
fn gradient(scaled: bool) -> f64 {
    let case = json!({
        "mesh": {"dimensions": [60, 10, 10], "lengths": {"x": 0.048, "y": A, "z": A}},
        "method": "piso",
        "config": {"dt": 1e-3, "kinematic_viscosity": NU, "density": 1.2, "pressure_preconditioner": "jacobi",
                   "pressure_tolerance": 1e-10, "piso_neighbour_update": scaled, "scaled_momentum": scaled},
        "output": {"steps": 1500, "output_interval": 1500},
        "initial_velocity": {"x": 0.0, "y": 0.0, "z": 0.0},
        "velocity_boundaries": [
            {"patch": "left", "kind": "fixed_value", "value": {"x": U, "y": 0.0, "z": 0.0}},
            {"patch": "right", "kind": "zero_gradient"},
            {"patch": "bottom", "kind": "fixed_value"}, {"patch": "top", "kind": "fixed_value"},
            {"patch": "front", "kind": "fixed_value"}, {"patch": "back", "kind": "fixed_value"}
        ],
        "pressure_boundaries": [{"patch": "right", "kind": "fixed_value", "value": 0.0}],
        "goals": [
            {"name": "p0", "kind": "point_pressure", "point": [0.0244, A / 2.0 + 0.0004, A / 2.0 + 0.0004]},
            {"name": "p1", "kind": "point_pressure", "point": [0.0404, A / 2.0 + 0.0004, A / 2.0 + 0.0004]}
        ]
    });
    let case: CollocatedFlowCase = serde_json::from_value(case).unwrap();
    let scratch = Scratch::new(if scaled { "on" } else { "off" });
    let manifest = Case::CollocatedFlow3d(case).run(&scratch.0).unwrap();
    if scaled {
        assert!(manifest.capabilities.iter().any(|c| c == "scaled_momentum_v1"));
    }
    // Point pressure goals are in pascals.
    (manifest.metrics["goal_p0"] - manifest.metrics["goal_p1"]) / 1.2 / 0.016
}

#[test]
fn small_slow_flows_are_solved_with_scaled_momentum() {
    let exact = 56.91 * NU * U / (2.0 * A * A);
    let scaled = gradient(true);
    // 1.5 s is about seven viscous time constants a²/(2π²ν): developed. Ten cells across a
    // collocated grid leave a few per cent.
    assert!((scaled / exact - 1.0).abs() < 0.1, "scaled {scaled} exact {exact}");
    // Unscaled the momentum is never solved: nowhere near.
    let plain = gradient(false);
    assert!((plain / exact - 1.0).abs() > 0.2, "plain {plain} exact {exact}");
}
