//! Analytic and property checks for the ported multiphysics coupling group.
use cfd_solvers::{
    core::LinearConfig,
    fem::*,
    fvm::make_sheared_cartesian_hexa_mesh,
    multiphysics::*,
    optics::{SequentialOpticalSystem, SequentialSurface, SurfaceType},
    Complex64,
};

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs().max(1.0)
}

#[test]
fn field_registry_keeps_units_sizes_and_unique_names() {
    let mut r = FieldRegistry::default();
    let mut temperature = FieldMetadata::new("temperature", 1, 4);
    temperature.units = UnitSignature::kelvin();
    r.add(temperature.clone(), 300.0).unwrap();
    let mut heat = FieldMetadata::new("heat", 1, 4);
    heat.units = UnitSignature::watt_per_cubic_metre();
    r.add(heat, 0.0).unwrap();
    assert!(r.add(temperature, 0.0).is_err(), "duplicate");
    assert!(r.add(FieldMetadata::new("", 1, 1), 0.0).is_err());
    assert!(r.add(FieldMetadata::new("zero", 0, 1), 0.0).is_err());
    assert_eq!(r.names(), ["heat", "temperature"]);
    r.data_mut("temperature").unwrap()[2] = 310.0;
    assert_eq!(r.data("temperature").unwrap(), [300.0, 300.0, 310.0, 300.0]);
    assert!(r
        .metadata("temperature")
        .unwrap()
        .units
        .dimensionally_compatible(&UnitSignature::kelvin()));
    assert!(!UnitSignature::kelvin().dimensionally_compatible(&UnitSignature::metre()));
    assert!(r.data("missing").is_err());
}

#[test]
fn transfers_conserve_integrals_exactly() {
    let src = CellGrid1D {
        edge: vec![0.0, 0.3, 0.5, 1.0],
    };
    let dst = CellGrid1D {
        edge: vec![0.0, 0.25, 0.6, 0.8, 1.0],
    };
    let out = conservative_cell_average_transfer(&src, &[1.0, 3.0, 2.0], &dst).unwrap();
    let expected = [1.0, (0.05 + 0.6 + 0.2) / 0.35, 2.0, 2.0];
    for (o, e) in out.iter().zip(expected) {
        assert!(close(*o, e, 1e-14), "{out:?}");
    }
    let integral = |g: &CellGrid1D, v: &[f64]| {
        g.edge
            .windows(2)
            .zip(v)
            .map(|(w, v)| (w[1] - w[0]) * v)
            .sum::<f64>()
    };
    assert!(close(
        integral(&dst, &out),
        integral(&src, &[1.0, 3.0, 2.0]),
        1e-14
    ));
    assert!(conservative_cell_average_transfer(
        &src,
        &[1.0, 3.0, 2.0],
        &CellGrid1D {
            edge: vec![0.0, 0.9]
        }
    )
    .is_err());
    assert!(CellGrid1D {
        edge: vec![0.0, 0.0]
    }
    .validate()
    .is_err());
    let mesh = make_sheared_cartesian_hexa_mesh(3, 2, 1, 0.3, 0.2, 0.1, 0.0).unwrap();
    let flux: Vec<f64> = (0..mesh.face_count())
        .map(|i| 0.1 * i as f64 - 0.7)
        .collect();
    let rate = conservative_face_flux_to_cell_rate(&mesh, &flux).unwrap();
    let total: f64 = rate
        .iter()
        .zip(mesh.cells())
        .map(|(r, c)| r * c.volume)
        .sum();
    let boundary: f64 = mesh
        .faces()
        .iter()
        .zip(&flux)
        .filter(|(f, _)| f.neighbour.is_none())
        .map(|(_, q)| q)
        .sum();
    assert!(close(total, boundary, 1e-12), "interior fluxes cancel");
}

#[test]
fn partitioned_coupling_finds_the_fixed_point_with_aitken() {
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
    assert!(r.converged, "{r:?}");
    assert!(
        close(state[0], 0.739_085_133_215_160_6, 1e-9),
        "Dottie number"
    );
    let s = 0.1 * state[0].sin();
    let x1 = (0.3 + 0.5 * s) / 1.2;
    assert!(close(state[1], x1, 1e-9) && close(state[2], -0.4 * x1 + s, 1e-9));
    let plain = solve_partitioned_fixed_point(
        &mut [0.5, 1.0, -0.2],
        |x, y| {
            y[0] = x[0].cos();
            y[1] = 0.5 * x[2] + 0.3;
            y[2] = -0.4 * x[1] + 0.1 * x[0].sin();
            Ok(())
        },
        PartitionedCouplerConfig {
            aitken: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        r.iterations < plain.iterations,
        "Aitken accelerates: {} vs {}",
        r.iterations,
        plain.iterations
    );
    assert!(solve_partitioned_fixed_point(
        &mut [],
        |_, _| Ok(()),
        PartitionedCouplerConfig::default()
    )
    .is_err());
    assert!(solve_partitioned_fixed_point(
        &mut [1.0],
        |_, _| Ok(()),
        PartitionedCouplerConfig {
            initial_relaxation: 9.0,
            ..Default::default()
        }
    )
    .is_err());
}

#[test]
fn block_system_and_port_iteration_solve_their_coupled_equations() {
    let mut b = BlockCoupledLinearSystem::new(vec![2, 1]).unwrap();
    for (rb, r, cb, c, v) in [
        (0, 0, 0, 0, 4.0),
        (0, 0, 0, 1, -1.0),
        (0, 1, 0, 0, -1.0),
        (0, 1, 0, 1, 4.0),
        (0, 1, 1, 0, -1.0),
        (1, 0, 0, 1, -1.0),
        (1, 0, 1, 0, 3.0),
    ] {
        b.add(rb, r, cb, c, v).unwrap();
    }
    for (blk, r, v) in [(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0)] {
        b.add_rhs(blk, r, v).unwrap();
    }
    let mut x = vec![0.0; 3];
    assert!(b.solve(&mut x, 100, 10, 1e-12).unwrap().converged);
    let a = b.matrix().unwrap();
    let mut ax = vec![0.0; 3];
    a.multiply(&x, &mut ax).unwrap();
    for (l, r) in ax.iter().zip(b.rhs()) {
        assert!(close(*l, *r, 1e-10));
    }
    assert!(b.add(2, 0, 0, 0, 1.0).is_err() && b.add(1, 1, 0, 0, 1.0).is_err());
    assert!(BlockCoupledLinearSystem::new(vec![2, 0]).is_err());
    let z = Complex64::new(50.0, 10.0);
    let r = iterate_em_spice_port(
        Complex64::new(0.0, 0.0),
        |i| z * i,
        |v| (Complex64::new(5.0, 0.0) - v) / 75.0,
        EmSpiceCouplingConfig::default(),
    )
    .unwrap();
    let exact = Complex64::new(5.0, 0.0) / (Complex64::new(75.0, 0.0) + z);
    assert!(
        r.converged && (r.port_current - exact).norm() < 1e-9,
        "{r:?}"
    );
}

#[test]
fn interface_laws_match_closed_forms() {
    let q = (350.0 - 300.0) / (0.002 / 200.0 + 0.001 / 0.6);
    let h = conjugate_interface_heat_flux(350.0, 200.0, 0.002, 300.0, 0.6, 0.001).unwrap();
    assert!(
        close(h.heat_flux, q, 1e-14) && close(h.interface_temperature, 350.0 - q * 1e-5, 1e-14)
    );
    assert!(conjugate_interface_heat_flux(1.0, 0.0, 1.0, 1.0, 1.0, 1.0).is_err());
    assert_eq!(
        electrochemical_heat_source(&[10.0, -5.0], &[0.1, 0.2], &[1.0, 2.0]).unwrap(),
        [2.0, 1.0]
    );
    assert!(electrochemical_heat_source(&[1.0], &[1.0, 2.0], &[]).is_err());
    let c = PorousElectroThermalFlowConfig::default();
    let mut s = PorousElectroThermalFlowState {
        species_concentration_mol_m3: 2.0,
        ..Default::default()
    };
    let mut expect = (2.0f64, 300.0f64, 0.0f64);
    for _ in 0..3 {
        s = advance_porous_electrothermal_flow(s, &c, 0.1, -1000.0).unwrap();
        let rate = expect.0;
        expect.0 = (expect.0 - 0.1 * rate).max(0.0);
        expect.1 += 0.1 * rate * 1e4 / 1e6;
        expect.2 += 1000.0 * (1e-12 / 1e-3 * 1000.0) * 0.1;
    }
    assert!(
        close(s.species_concentration_mol_m3, expect.0, 1e-14)
            && close(s.temperature_k, expect.1, 1e-14)
            && close(s.pressure_pa, expect.2, 1e-12)
    );
    assert!(advance_porous_electrothermal_flow(s, &c, -1.0, 0.0).is_err());
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
    assert!(
        close(moved.surfaces()[0].vertex_z, 1e-4, 1e-14)
            && close(moved.surfaces()[1].vertex_z, 0.01 + 1e-4 + 1e-4, 1e-14)
    );
    assert!(deform_optical_system_axially(&sys, |_| f64::NAN).is_err());
}

#[test]
fn projections_are_exact_for_linear_and_constant_fields() {
    let lin = |x: f64, y: f64, z: f64| 1.0 + 2.0 * x - 0.5 * y + 0.25 * z;
    let mut grid = StructuredScalarGrid3D {
        nx: 3,
        ny: 3,
        nz: 2,
        lx: 1.0,
        ly: 1.0,
        lz: 1.0,
        values: vec![],
    };
    for k in 0..2 {
        for j in 0..3 {
            for i in 0..3 {
                grid.values
                    .push(lin(i as f64 / 2.0, j as f64 / 2.0, k as f64));
            }
        }
    }
    assert!(close(
        grid.sample(0.3, 0.6, 0.25).unwrap(),
        lin(0.3, 0.6, 0.25),
        1e-14
    ));
    assert!(
        close(
            grid.sample(1.2, -0.1, 0.9).unwrap(),
            lin(1.0, 0.0, 0.9),
            1e-14
        ),
        "clamped"
    );
    let cells = make_sheared_cartesian_hexa_mesh(4, 3, 1, 1.0, 1.0, 0.1, 0.0).unwrap();
    for (v, c) in structured_to_fvm_cells(&grid, &cells)
        .unwrap()
        .iter()
        .zip(cells.cells())
    {
        assert!(close(*v, lin(c.center.x, c.center.y, c.center.z), 1e-13));
    }
    let tri = make_rectangle_tri_mesh(4, 4, 1.0, 1.0).unwrap();
    let nodal: Vec<f64> = tri.nodes.iter().map(|n| lin(n.x, n.y, 0.0)).collect();
    for (v, c) in fem_nodes_to_fvm_cells(&tri, &nodal, &cells)
        .unwrap()
        .iter()
        .zip(cells.cells())
    {
        assert!(
            close(*v, lin(c.center.x, c.center.y, 0.0), 1e-13),
            "P1 is exact for linear data"
        );
    }
    let constant = vec![4.25; cells.cell_count()];
    assert!(fvm_cells_to_fem_nodes(&cells, &constant, &tri, 4)
        .unwrap()
        .iter()
        .all(|v| close(*v, 4.25, 1e-14)));
    let outside = make_sheared_cartesian_hexa_mesh(2, 2, 1, 2.0, 1.0, 0.1, 0.0).unwrap();
    assert!(fem_nodes_to_fvm_cells(&tri, &nodal, &outside).is_err());
    assert!(StructuredScalarGrid3D {
        nx: 1,
        ..grid.clone()
    }
    .validate()
    .is_err());
}

#[test]
fn joule_heating_matches_uniform_current_and_heats_the_interior() {
    let mesh = make_rectangle_tri_mesh(4, 3, 0.4, 0.3).unwrap();
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
            relative_tolerance: 1e-13,
            restart: 30,
        },
    )
    .unwrap();
    joule.set_electrical_boundary(0, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 1.0)));
    joule.set_electrical_boundary(1, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 0.0)));
    assert!(
        joule.thermal_step().is_err(),
        "the electrical field must be solved first"
    );
    joule.initialize_temperature(|_| 300.0).unwrap();
    joule.set_thermal_dirichlet(|_, _| 300.0);
    joule.solve_electrical(&|_| 2.0, None).unwrap();
    // Uniform field 1 V / 0.4 m: J = 5 A/m^2 and q = J^2 / sigma = 12.5 W/m^3 in every element.
    assert!(
        joule
            .joule_heating_density()
            .iter()
            .all(|q| close(*q, 12.5, 1e-9)),
        "{:?}",
        joule.joule_heating_density()
    );
    let j = joule
        .electrical()
        .element_current_density(&|_| 2.0)
        .unwrap();
    assert!(j.iter().all(|v| close(v.x, 5.0, 1e-9) && v.y.abs() < 1e-9));
    joule.thermal_run(4).unwrap();
    let t = &joule.thermal().temperature;
    let mesh = &joule.thermal().mesh;
    for (i, temp) in t.iter().enumerate() {
        if mesh.boundary_node[i] {
            assert!(close(*temp, 300.0, 1e-12));
        } else {
            assert!(
                *temp > 300.0 && *temp < 300.0 + 4.0 * 0.05 * 12.5 / 2.0 + 1e-12,
                "bounded by adiabatic heating: {temp}"
            );
        }
    }
}

#[test]
fn electrostatics_gives_the_uniform_capacitor_field() {
    let mut es = Electrostatics2D::new(
        make_rectangle_tri_mesh(3, 3, 1.0, 1.0).unwrap(),
        LinearConfig {
            max_iterations: 1000,
            relative_tolerance: 1e-13,
            restart: 30,
        },
    )
    .unwrap();
    es.set_boundary(2, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 0.0)));
    es.set_boundary(3, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 10.0)));
    es.solve(&|_| 4.0, &|_| 0.0).unwrap();
    assert!(es
        .element_electric_field()
        .unwrap()
        .iter()
        .all(|e| e.x.abs() < 1e-9 && close(e.y, -10.0, 1e-9)));
}

#[test]
fn bioheat_sar_feedback_averaging_and_pennes_balance() {
    assert!(close(
        sar_from_rms_electric_field(0.8, 1050.0, 40.0).unwrap(),
        0.8 * 1600.0 / 1050.0,
        1e-15
    ));
    assert!(sar_from_rms_electric_field(-1.0, 1.0, 1.0).is_err());
    let f = temperature_feedback(&TemperatureFeedbackConfig::default(), 312.0).unwrap();
    let d = 312.0 - 310.15;
    assert!(
        close(f.conductivity_s_per_m, 1.0 + 0.01 * d, 1e-15)
            && close(f.perfusion_per_s, 0.01 * (1.0 + 0.02 * d), 1e-15)
    );
    let materials = vec![
        VoxelTissueMaterial::default(),
        VoxelTissueMaterial {
            density_kg_per_m3: 1100.0,
            electrical_conductivity_s_per_m: 0.3,
            ..Default::default()
        },
    ];
    let grid = VoxelTissueGrid3D {
        nx: 3,
        ny: 2,
        nz: 2,
        dx_m: 0.01,
        dy_m: 0.01,
        dz_m: 0.01,
        materials,
        material_index: (0..12).map(|i| i % 2).collect(),
        electric_rms_v_per_m: (0..12).map(|i| 10.0 + 3.0 * i as f64).collect(),
    };
    let sar = voxel_sar_w_per_kg(&grid).unwrap();
    // Averaging mass below one voxel's mass returns the largest voxel SAR.
    let peak = max_mass_averaged_sar_w_per_kg(&grid, &sar, 1e-4).unwrap();
    assert!(close(peak, sar.iter().cloned().fold(0.0, f64::max), 1e-14));
    // Averaging the whole grid's mass returns the global mass-weighted mean for every centre.
    let total_mass: f64 = grid
        .material_index
        .iter()
        .map(|m| grid.materials[*m].density_kg_per_m3 * 1e-6)
        .sum();
    let mean = grid
        .material_index
        .iter()
        .zip(&sar)
        .map(|(m, s)| grid.materials[*m].density_kg_per_m3 * 1e-6 * s)
        .sum::<f64>()
        / total_mass;
    assert!(close(
        max_mass_averaged_sar_w_per_kg(&grid, &sar, total_mass).unwrap(),
        mean,
        1e-12
    ));
    let one = project_voxel_sar_to_pennes2d(&grid, &sar, 1, 1).unwrap();
    assert!(
        close(one[0], mean, 1e-12),
        "a 1x1 projection is the mass-weighted mean"
    );
    // Periodic Pennes with uniform heating and no perfusion: implicit Euler adds Q dt/(rho c) per step.
    let config = PennesBioheat2DConfig {
        nx: 6,
        ny: 5,
        dt_s: 2.0,
        boundary: BioheatBoundary::Periodic,
        ..Default::default()
    };
    let mut p = PennesBioheat2D::new(config).unwrap();
    p.set_uniform_volumetric_heating(1000.0).unwrap();
    p.step(3).unwrap();
    let rise = 3.0 * 2.0 * 1000.0 / (1000.0 * 3600.0);
    assert!(p
        .temperature_k()
        .iter()
        .all(|t| close(*t, 310.15 + rise, 1e-10)));
    // With perfusion the uniform implicit update is T' = (C T + w Tb + Q) / (C + w).
    let mut p = PennesBioheat2D::new(PennesBioheat2DConfig {
        blood_perfusion_per_s: 0.005,
        ..config
    })
    .unwrap();
    p.initialize(305.0).unwrap();
    p.set_uniform_volumetric_heating(1000.0).unwrap();
    p.step(1).unwrap();
    let (cap, w) = (1000.0 * 3600.0 / 2.0, 1060.0 * 3770.0 * 0.005);
    assert!(p.temperature_k().iter().all(|t| close(
        *t,
        (cap * 305.0 + w * 310.15 + 1000.0) / (cap + w),
        1e-10
    )));
    // Fixed-temperature boundary is held and the interior heats.
    let mut p = PennesBioheat2D::new(PennesBioheat2DConfig {
        boundary: BioheatBoundary::FixedTemperature,
        ..config
    })
    .unwrap();
    p.set_uniform_volumetric_heating(1000.0).unwrap();
    p.step(3).unwrap();
    for y in 0..5 {
        for x in 0..6 {
            let t = p.temperature_k()[y * 6 + x];
            if x == 0 || y == 0 || x == 5 || y == 4 {
                assert!(close(t, 310.15, 1e-12));
            } else {
                assert!(t > 310.15 && t < 310.15 + rise + 1e-12);
            }
        }
    }
    assert!(PennesBioheat2D::new(PennesBioheat2DConfig { nx: 2, ..config }).is_err());
}
