//! Value-by-value comparison with the C++ upstream multiphysics units.
//! Goldens come from `parity/oracle_multiphysics.cpp` built against the pinned upstream revision.
use cfd_solvers::{
    core::LinearConfig,
    fem::*,
    fvm::make_sheared_cartesian_hexa_mesh,
    multiphysics::*,
    optics::{SequentialOpticalSystem, SequentialSurface, SurfaceType},
    Complex64,
};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-multiphysics.json")).unwrap()
}

#[track_caller]
fn same(name: &str, rust: &[f64], tol: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= tol * c.abs().max(1.0),
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn scalar(name: &str) -> f64 {
    goldens()[name].as_f64().unwrap()
}

#[test]
fn transfers_and_projections_match_upstream() {
    let src = CellGrid1D {
        edge: vec![0.0, 0.3, 0.5, 1.0],
    };
    let dst = CellGrid1D {
        edge: vec![0.0, 0.25, 0.6, 0.8, 1.0],
    };
    same(
        "transfer_1d",
        &conservative_cell_average_transfer(&src, &[1.0, 3.0, 2.0], &dst).unwrap(),
        1e-15,
    );
    let mesh = make_sheared_cartesian_hexa_mesh(3, 2, 1, 0.3, 0.2, 0.1, 0.0).unwrap();
    let flux: Vec<f64> = (0..mesh.face_count())
        .map(|i| 0.1 * i as f64 - 0.7)
        .collect();
    same(
        "face_flux_rate",
        &conservative_face_flux_to_cell_rate(&mesh, &flux).unwrap(),
        1e-14,
    );
    let grid = StructuredScalarGrid3D {
        nx: 3,
        ny: 3,
        nz: 2,
        lx: 1.0,
        ly: 1.0,
        lz: 1.0,
        values: (0..18).map(|i| (0.7 * i as f64).sin()).collect(),
    };
    same(
        "structured_samples",
        &[
            grid.sample(0.3, 0.6, 0.25).unwrap(),
            grid.sample(1.2, -0.1, 0.9).unwrap(),
        ],
        1e-15,
    );
    let cells = make_sheared_cartesian_hexa_mesh(4, 3, 1, 1.0, 1.0, 0.1, 0.0).unwrap();
    same(
        "structured_to_cells",
        &structured_to_fvm_cells(&grid, &cells).unwrap(),
        1e-14,
    );
    let tri = make_rectangle_tri_mesh(4, 4, 1.0, 1.0).unwrap();
    let nodal: Vec<f64> = tri
        .nodes
        .iter()
        .map(|n| 1.0 + 2.0 * n.x - 0.5 * n.y)
        .collect();
    same(
        "fem_to_cells",
        &fem_nodes_to_fvm_cells(&tri, &nodal, &cells).unwrap(),
        1e-14,
    );
    let cellv: Vec<f64> = cells
        .cells()
        .iter()
        .map(|c| c.center.x * c.center.x)
        .collect();
    same(
        "cells_to_fem",
        &fvm_cells_to_fem_nodes(&cells, &cellv, &tri, 4).unwrap(),
        1e-13,
    );
}

#[test]
fn coupling_algorithms_match_upstream() {
    let mut state = vec![0.5, 1.0, -0.2];
    let r = solve_partitioned_fixed_point(
        &mut state,
        |x, y| {
            y[0] = x[0].cos();
            y[1] = 0.5 * x[2] + 0.3;
            y[2] = -0.4 * x[1] + 0.1 * x[0].sin();
            Ok(())
        },
        PartitionedCouplerConfig::default(),
    )
    .unwrap();
    same("partitioned_state", &state, 1e-13);
    assert_eq!(r.iterations as f64, scalar("partitioned_iterations"));
    assert!((r.relaxation - scalar("partitioned_relaxation")).abs() < 1e-12);
    let mut b = BlockCoupledLinearSystem::new(vec![2, 1]).unwrap();
    for (rb, row, cb, col, v) in [
        (0, 0, 0, 0, 4.0),
        (0, 0, 0, 1, -1.0),
        (0, 1, 0, 0, -1.0),
        (0, 1, 0, 1, 4.0),
        (0, 1, 1, 0, -1.0),
        (1, 0, 0, 1, -1.0),
        (1, 0, 1, 0, 3.0),
    ] {
        b.add(rb, row, cb, col, v).unwrap();
    }
    for (blk, row, v) in [(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0)] {
        b.add_rhs(blk, row, v).unwrap();
    }
    let mut x = vec![0.0; 3];
    b.solve(&mut x, 1000, 30, 1e-10).unwrap();
    same("block_solution", &x, 1e-9);
    let z = Complex64::new(50.0, 10.0);
    let em = iterate_em_spice_port(
        Complex64::new(0.0, 0.0),
        |i| z * i,
        |v| (Complex64::new(5.0, 0.0) - v) / 75.0,
        EmSpiceCouplingConfig::default(),
    )
    .unwrap();
    same(
        "em_spice",
        &[
            em.port_voltage.re,
            em.port_voltage.im,
            em.port_current.re,
            em.port_current.im,
            em.iterations as f64,
        ],
        1e-13,
    );
}

#[test]
fn interface_laws_match_upstream() {
    let h = conjugate_interface_heat_flux(350.0, 200.0, 0.002, 300.0, 0.6, 0.001).unwrap();
    same("cht", &[h.heat_flux, h.interface_temperature], 1e-15);
    same(
        "electrochemical_heat",
        &electrochemical_heat_source(&[10.0, -5.0], &[0.1, 0.2], &[1.0, 2.0]).unwrap(),
        1e-15,
    );
    let mut s = PorousElectroThermalFlowState {
        species_concentration_mol_m3: 2.0,
        temperature_k: 300.0,
        pressure_pa: 0.0,
        darcy_velocity_m_per_s: 0.0,
    };
    for _ in 0..3 {
        s = advance_porous_electrothermal_flow(
            s,
            &PorousElectroThermalFlowConfig::default(),
            0.1,
            -1000.0,
        )
        .unwrap();
    }
    same(
        "porous_electrothermal",
        &[
            s.species_concentration_mol_m3,
            s.temperature_k,
            s.pressure_pa,
            s.darcy_velocity_m_per_s,
        ],
        1e-15,
    );
    let mut sys = SequentialOpticalSystem::new(1.0).unwrap();
    sys.add_surface(SequentialSurface {
        surface_type: SurfaceType::Sphere,
        vertex_z: 0.0,
        radius: 0.05,
        refractive_index_after: 1.5,
        ..Default::default()
    })
    .unwrap();
    sys.add_surface(SequentialSurface {
        vertex_z: 0.01,
        ..Default::default()
    })
    .unwrap();
    let moved = deform_optical_system_axially(&sys, |z| 1e-4 + 0.01 * z).unwrap();
    same(
        "optics_vertices",
        &[moved.surfaces()[0].vertex_z, moved.surfaces()[1].vertex_z],
        1e-15,
    );
}

#[test]
fn joule_heating_matches_upstream() {
    let mesh = make_rectangle_tri_mesh(4, 3, 0.4, 0.3).unwrap();
    // Upstream defaults: scalar diffusion 4000 iterations at 1e-11; Heat2D linear defaults as ported.
    let mut joule = JouleHeatingCoupler2D::new(
        mesh,
        Heat2DConfig {
            conductivity: 1.5,
            volumetric_heat_capacity: 2.0,
            dt: 0.05,
            ..Default::default()
        },
        LinearConfig {
            max_iterations: 4000,
            relative_tolerance: 1e-11,
            restart: 30,
        },
    )
    .unwrap();
    joule.set_electrical_boundary(0, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 1.0)));
    joule.set_electrical_boundary(1, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 0.0)));
    joule.initialize_temperature(|_| 300.0).unwrap();
    joule.set_thermal_dirichlet(|_, _| 300.0);
    joule.solve_electrical(&|p| 2.0 + p.x, None).unwrap();
    joule.thermal_run(4).unwrap();
    same("joule_density", joule.joule_heating_density(), 1e-8);
    same("joule_temperature", &joule.thermal().temperature, 1e-10);
}

#[test]
fn bioheat_matches_upstream() {
    assert!(
        (sar_from_rms_electric_field(0.8, 1050.0, 40.0).unwrap() - scalar("sar")).abs() < 1e-15
    );
    let f = temperature_feedback(&TemperatureFeedbackConfig::default(), 312.0).unwrap();
    same(
        "feedback",
        &[f.conductivity_s_per_m, f.perfusion_per_s],
        1e-15,
    );
    let grid = VoxelTissueGrid3D {
        nx: 3,
        ny: 2,
        nz: 2,
        dx_m: 0.01,
        dy_m: 0.01,
        dz_m: 0.01,
        materials: vec![
            VoxelTissueMaterial::default(),
            VoxelTissueMaterial {
                density_kg_per_m3: 1100.0,
                specific_heat_j_per_kg_k: 3500.0,
                thermal_conductivity_w_per_m_k: 0.4,
                electrical_conductivity_s_per_m: 0.3,
                ..Default::default()
            },
        ],
        material_index: (0..12).map(|i| i % 2).collect(),
        electric_rms_v_per_m: (0..12).map(|i| 10.0 + 3.0 * i as f64).collect(),
    };
    let sar = voxel_sar_w_per_kg(&grid).unwrap();
    same("voxel_sar", &sar, 1e-15);
    assert!(
        (max_mass_averaged_sar_w_per_kg(&grid, &sar, 0.002).unwrap() - scalar("peak_sar_2g")).abs()
            < 1e-12 * scalar("peak_sar_2g")
    );
    same(
        "sar_projection",
        &project_voxel_sar_to_pennes2d(&grid, &sar, 2, 2).unwrap(),
        1e-14,
    );
    for (boundary, name) in [
        (BioheatBoundary::Periodic, "pennes_periodic"),
        (BioheatBoundary::FixedTemperature, "pennes_fixed"),
    ] {
        let mut p = PennesBioheat2D::new(PennesBioheat2DConfig {
            nx: 6,
            ny: 5,
            blood_perfusion_per_s: 0.005,
            metabolic_heat_w_per_m3: 400.0,
            dt_s: 2.0,
            boundary,
            ..Default::default()
        })
        .unwrap();
        let map: Vec<f64> = (0..30).map(|i| 0.5 * (i % 7) as f64).collect();
        p.set_sar(&map).unwrap();
        p.set_uniform_volumetric_heating(1000.0).unwrap();
        p.step(3).unwrap();
        same(name, p.temperature_k(), 1e-12);
    }
}
