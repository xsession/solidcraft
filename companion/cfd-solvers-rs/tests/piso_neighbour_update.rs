//! PISO correctors that rebuild H(u)/a_P from the corrected velocity
//! (`piso_neighbour_update`): stable well past the diffusion-number limit of the
//! C++-parity correctors, and the default stays the C++ behaviour.

use cfd_solvers::fvm::{CollocatedIncompressible, CollocatedIncompressibleConfig, VelocityBoundaryType};
use cfd_solvers::{case::*, core::Vec3};

/// Lid-driven cavity, n × n, side 0.1 m, lid 0.05 m/s; returns the largest speed
/// after `steps` PISO steps, or the error.
fn cavity(n: usize, dt: f64, nu: f64, steps: usize, update: bool) -> Result<f64, String> {
    let case = CollocatedFlowCase {
        mesh: PolyMeshSpec { dimensions: [n, n, 1], lengths: Vec3::new(0.1, 0.1, 0.1 / n as f64), xy_shear: 0.0, coordinates: None },
        ..Default::default()
    };
    let (mesh, _) = case.build_mesh().map_err(|e| e.to_string())?;
    let config = CollocatedIncompressibleConfig {
        dt,
        kinematic_viscosity: nu,
        pressure_tolerance: 1e-10,
        pressure_iterations: 4000,
        piso_neighbour_update: update,
        ..Default::default()
    };
    let mut s = CollocatedIncompressible::new(mesh, config).map_err(|e| e.to_string())?;
    for (name, kind) in [
        ("left", VelocityBoundaryType::FixedValue),
        ("right", VelocityBoundaryType::FixedValue),
        ("bottom", VelocityBoundaryType::FixedValue),
        ("top", VelocityBoundaryType::FixedValue),
        ("front", VelocityBoundaryType::Slip),
        ("back", VelocityBoundaryType::Slip),
    ] {
        let value = if name == "top" { Vec3::new(0.05, 0.0, 0.0) } else { Vec3::ZERO };
        s.set_velocity_boundary(name, kind, value).map_err(|e| e.to_string())?;
    }
    for _ in 0..steps {
        s.step_piso().map_err(|e| e.to_string())?;
    }
    Ok(s.velocity().iter().map(|v| v.norm()).fold(0.0, f64::max))
}

#[test]
fn neighbour_updates_keep_a_diffusion_dominated_cavity_stable() {
    // ν Δt / h² = 1e-2 · 1e-3 / 0.005² = 0.4.
    let fixed = cavity(20, 1e-3, 1e-2, 200, true).unwrap();
    assert!(fixed > 0.03 && fixed < 0.05, "settled flow stays below the lid speed: {fixed}");
    // The C++-parity correctors run away on the same case: either the speed grows
    // past the lid speed or the momentum solve gives up.
    match cavity(20, 1e-3, 1e-2, 200, false) {
        Ok(speed) => assert!(speed > 0.05, "{speed}"),
        Err(e) => assert!(e.contains("did not converge"), "{e}"),
    }
}

#[test]
fn neighbour_updates_agree_with_the_parity_correctors_where_both_are_stable() {
    // ν Δt / h² = 1e-3 · 1e-3 / 0.01² = 0.01: both settle to the same cavity flow.
    let a = cavity(10, 1e-3, 1e-3, 300, true).unwrap();
    let b = cavity(10, 1e-3, 1e-3, 300, false).unwrap();
    assert!((a - b).abs() < 1e-3 * b, "{a} vs {b}");
}

#[test]
fn the_flag_defaults_off_and_old_cases_still_load() {
    assert!(!CollocatedIncompressibleConfig::default().piso_neighbour_update);
    let mut json = serde_json::to_value(CollocatedIncompressibleConfig::default()).unwrap();
    json.as_object_mut().unwrap().remove("piso_neighbour_update");
    let back: CollocatedIncompressibleConfig = serde_json::from_value(json).unwrap();
    assert!(!back.piso_neighbour_update);
}
