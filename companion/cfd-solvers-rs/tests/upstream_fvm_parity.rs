//! Numerical comparisons to original C++ source functions, not a second Rust implementation.
use cfd_solvers::{core::*, fvm::*, Complex64 as C};
use serde_json::Value;
fn fixture() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-fvm.json")).unwrap()
}
fn check(g: &Value, name: &str, v: &[f64], atol: f64, rtol: f64) {
    let e = g[name].as_array().unwrap();
    assert_eq!(e.len(), v.len(), "{name} shape");
    for (i, (x, y)) in v.iter().zip(e).enumerate() {
        let y = y.as_f64().unwrap();
        let limit = atol + rtol * y.abs();
        assert!(
            x.is_finite() && (x - y).abs() <= limit,
            "{name}[{i}] Rust={x:.17e}, C++={y:.17e}, error={:.4e}, bound={limit:.4e}",
            (x - y).abs()
        );
    }
}
fn mesh(shear: f64) -> PolyMesh {
    make_sheared_cartesian_hexa_mesh(4, 3, 2, 1.2, 0.8, 0.5, shear).unwrap()
}
fn scalar(p: Vec3) -> f64 {
    0.3 + 1.2 * p.x - 0.7 * p.y + 0.4 * p.z + 0.08 * p.x * p.x
}
fn vector(p: Vec3) -> Vec3 {
    Vec3::new(0.2 + p.x * p.x, 0.3 - 0.1 * p.y, 0.5 * p.z * p.z)
}
fn vectors(out: &mut Vec<f64>, v: &[Vec3]) {
    out.extend(v.iter().flat_map(|p| p.as_array()));
}
#[test]
fn block_complex_and_multigrid_match_cpp() {
    let g = fixture();
    let mut b = BlockCsrBuilder::<2>::new(3, 3).unwrap();
    for i in 0..3 {
        b.add(i, i, [[4.0, 0.7], [-0.2, 3.0]]).unwrap();
        if i + 1 < 3 {
            b.add(i, i + 1, [[-1.0, 0.1], [0.2, -0.5]]).unwrap();
            b.add(i + 1, i, [[-0.8, 0.0], [0.1, -0.6]]).unwrap();
        }
    }
    b.add(0, 0, [[1.0, 0.0], [0.0, 1.0]]).unwrap();
    let a = b.build().unwrap();
    let rhs = a.apply(&[1.0, -2.0, 0.3, 0.7, -1.0, 2.0]).unwrap();
    let mut x = vec![0.0; 6];
    let result = restarted_block_gmres(
        &rhs,
        &mut x,
        |v, y| a.multiply(v, y),
        &a.diagonal_blocks(),
        LinearConfig {
            max_iterations: 100,
            restart: 8,
            relative_tolerance: 1e-12,
        },
    )
    .unwrap();
    assert!(result.converged);
    let mut out = rhs;
    out.extend(x);
    check(&g, "block_sparse", &out, 1e-11, 1e-11);
    let mut a = ComplexCsrBuilder::new(3).unwrap();
    for (i, j, re, im) in [
        (0, 0, 3.0, 1.0),
        (0, 1, 0.2, -0.3),
        (1, 0, -0.4, 0.1),
        (1, 1, 2.0, 0.4),
        (1, 2, 0.0, 0.2),
        (2, 1, 0.1, 0.0),
        (2, 2, 4.0, -0.5),
    ] {
        a.add(i, j, C::new(re, im)).unwrap();
    }
    let rhs = a
        .apply(&[C::new(1.0, 0.2), C::new(-0.3, 0.8), C::new(0.7, -0.1)])
        .unwrap();
    let r = solve_complex_sparse(&a, &rhs, Default::default()).unwrap();
    assert!(r.linear_result.converged);
    let out: Vec<_> = r.solution.iter().flat_map(|z| [z.re, z.im]).collect();
    check(&g, "complex_sparse", &out, 1e-9, 1e-9);
    let n = 15;
    let h = 1.0 / (n + 1) as f64;
    let lambda = 8.0 * (0.5 * std::f64::consts::PI * h).sin().powi(2) / (h * h);
    let rhs: Vec<_> = (0..n)
        .flat_map(|j| {
            (0..n).map(move |i| {
                lambda
                    * (std::f64::consts::PI * (i + 1) as f64 * h).sin()
                    * (std::f64::consts::PI * (j + 1) as f64 * h).sin()
            })
        })
        .collect();
    let mut u = vec![0.0; n * n];
    assert!(
        PoissonMultigrid2D::new(n)
            .unwrap()
            .solve(&rhs, &mut u, 200, 1e-10)
            .unwrap()
            .converged
    );
    check(&g, "multigrid", &u, 1e-10, 1e-10);
}
#[test]
fn mesh_geometry_and_corrected_operators_match_cpp() {
    let g = fixture();
    for scenario in 0..2 {
        let m = mesh(if scenario == 0 { 0.0 } else { 0.35 });
        let mut geometry = vec![];
        for c in m.cells() {
            geometry.extend([c.center.x, c.center.y, c.center.z, c.volume]);
        }
        for f in m.faces() {
            geometry.extend([
                f.center.x,
                f.center.y,
                f.center.z,
                f.area.x,
                f.area.y,
                f.area.z,
                f.owner as f64,
                f.neighbour.map_or(-1.0, |i| i as f64),
                f.patch.map_or(-1.0, |i| i as f64),
            ]);
        }
        check(&g, &format!("mesh_{scenario}"), &geometry, 1e-14, 1e-14);
        let p: Vec<_> = m.cells().iter().map(|c| scalar(c.center)).collect();
        let b: Vec<_> = m.faces().iter().map(|f| scalar(f.center)).collect();
        let v: Vec<_> = m.cells().iter().map(|c| vector(c.center)).collect();
        let vb: Vec<_> = m.faces().iter().map(|f| vector(f.center)).collect();
        let mut out = vec![];
        vectors(&mut out, &gauss_gradient_scalar(&m, &p, &b).unwrap());
        vectors(
            &mut out,
            &least_squares_gradient_scalar(&m, &p, &b).unwrap(),
        );
        out.extend(orthogonal_laplacian_scalar(&m, &p, 0.3, &b).unwrap());
        out.extend(corrected_laplacian_scalar(&m, &p, 0.3, &b).unwrap());
        out.extend(corrected_face_normal_gradient_scalar(&m, &p, &b).unwrap());
        out.extend(gauss_divergence_vector(&m, &v, &vb).unwrap());
        check(&g, &format!("operators_{scenario}"), &out, 2e-10, 2e-10);
    }
}
#[test]
fn all_interpolation_limiters_and_rhie_chow_match_cpp() {
    let g = fixture();
    let m = mesh(0.2);
    let signal = |p: Vec3| (5.0 * p.x).sin() + 0.2 * p.y - 0.1 * p.z;
    let flow = |p: Vec3| Vec3::new(-p.y + 0.2, p.x - 0.1, 0.05 * p.z);
    let p: Vec<_> = m.cells().iter().map(|c| signal(c.center)).collect();
    let b: Vec<_> = m.faces().iter().map(|f| signal(f.center)).collect();
    let u: Vec<_> = m.cells().iter().map(|c| flow(c.center)).collect();
    let ub: Vec<_> = m.faces().iter().map(|f| flow(f.center)).collect();
    let flux = face_flux_from_velocity(&m, &u, &ub).unwrap();
    for (index, s) in [
        FaceInterpolationScheme::Linear,
        FaceInterpolationScheme::Upwind,
        FaceInterpolationScheme::BoundedLinear,
        FaceInterpolationScheme::Minmod,
        FaceInterpolationScheme::VanLeer,
    ]
    .into_iter()
    .enumerate()
    {
        let mut out = interpolate_scalar_to_faces(&m, &p, s, &flux, &b).unwrap();
        out.extend(convective_divergence_scalar(&m, &p, &flux, s, &b).unwrap());
        vectors(
            &mut out,
            &convective_divergence_vector(&m, &u, &flux, s, &ub).unwrap(),
        );
        out.push(courant_limited_timestep(&m, &flux, 0.4, 1.0).unwrap());
        check(&g, &format!("scheme_{index}"), &out, 2e-10, 2e-10);
    }
    let p: Vec<_> = m.cells().iter().map(|c| scalar(c.center)).collect();
    let mobility: Vec<_> = m.cells().iter().map(|c| 0.03 + 0.01 * c.center.x).collect();
    let h: Vec<_> = m
        .cells()
        .iter()
        .map(|c| Vec3::new(0.2 + 0.1 * c.center.y, -0.1 * c.center.x, 0.01 * c.center.z))
        .collect();
    let mut ub = vec![VelocityBoundaryCondition::default(); 6];
    ub[0] = VelocityBoundaryCondition {
        boundary_type: VelocityBoundaryType::FixedValue,
        value: Vec3::new(0.2, 0.0, 0.0),
    };
    ub[2] = VelocityBoundaryCondition {
        boundary_type: VelocityBoundaryType::Slip,
        value: Vec3::ZERO,
    };
    let mut pb = vec![PressureBoundaryCondition::default(); 6];
    pb[1] = PressureBoundaryCondition {
        boundary_type: PressureBoundaryType::FixedValue,
        value: 1.3,
    };
    for corrected in 0..2 {
        let mut out = vec![];
        vectors(&mut out, &boundary_velocity_values(&m, &h, &ub).unwrap());
        out.extend(boundary_pressure_values(&m, &p, &pb).unwrap());
        out.extend(rhie_chow_face_flux(&m, &h, &mobility, &p, &ub, &pb, corrected == 1).unwrap());
        check(&g, &format!("rhie_{corrected}"), &out, 2e-11, 2e-11);
    }
}
#[test]
fn scalar_transport_history_and_operator_cache_match_cpp() {
    let g = fixture();
    for (mode, scheme) in [
        TemporalScheme::Euler,
        TemporalScheme::BackwardBdf2,
        TemporalScheme::CrankNicolson,
    ]
    .into_iter()
    .enumerate()
    {
        let m = mesh(0.2);
        let mut s = ScalarTransport::new(
            m.clone(),
            ScalarTransportConfig {
                dt: 0.003,
                diffusivity: 0.02,
                temporal_scheme: scheme,
                linear: LinearConfig {
                    max_iterations: 400,
                    restart: 30,
                    relative_tolerance: 1e-12,
                },
                ..Default::default()
            },
        )
        .unwrap();
        s.initialize_profile(|p| 1.0 + 0.1 * p.x + 0.05 * (5.0 * p.y).cos())
            .unwrap();
        s.set_boundary("left", ScalarBoundaryType::FixedValue, 1.0)
            .unwrap();
        s.set_boundary("right", ScalarBoundaryType::FixedValue, 0.3)
            .unwrap();
        s.set_source_profile(|p| 0.1 + 0.02 * p.z).unwrap();
        s.set_face_flux(
            m.faces()
                .iter()
                .map(|f| Vec3::new(0.1, 0.03, 0.0).dot(f.area))
                .collect(),
        )
        .unwrap();
        s.run(3).unwrap();
        s.set_boundary("right", ScalarBoundaryType::FixedValue, 1.0)
            .unwrap();
        s.run(2).unwrap();
        let mut out = s.values().to_vec();
        out.extend([
            s.time(),
            s.volume_integral(),
            s.minimum(),
            s.maximum(),
            s.operator_assemblies() as f64,
        ]);
        check(&g, &format!("scalar_{mode}"), &out, 2e-11, 2e-11);
    }
}
#[test]
fn simple_piso_pimple_and_temporal_flow_paths_match_cpp() {
    let g = fixture();
    for mode in 0..5 {
        let m = make_sheared_cartesian_hexa_mesh(5, 4, 1, 1.0, 1.0, 0.1, 0.15).unwrap();
        let mut s = CollocatedIncompressible::new(
            m.clone(),
            CollocatedIncompressibleConfig {
                dt: 0.003,
                density: 2.0,
                kinematic_viscosity: 0.02,
                momentum_tolerance: 1e-11,
                pressure_tolerance: 1e-12,
                pressure_preconditioner: PressurePreconditionerKind::Jacobi,
                temporal_scheme: if mode == 2 {
                    TemporalScheme::BackwardBdf2
                } else if mode == 3 {
                    TemporalScheme::CrankNicolson
                } else {
                    TemporalScheme::Euler
                },
                use_krylov_momentum: mode != 4,
                ..Default::default()
            },
        )
        .unwrap();
        for name in ["left", "right", "bottom"] {
            s.set_velocity_boundary(name, VelocityBoundaryType::FixedValue, Vec3::ZERO)
                .unwrap();
        }
        s.set_velocity_boundary(
            "top",
            VelocityBoundaryType::FixedValue,
            Vec3::new(0.1, 0.0, 0.0),
        )
        .unwrap();
        for name in ["front", "back"] {
            s.set_velocity_boundary(name, VelocityBoundaryType::Slip, Vec3::ZERO)
                .unwrap();
        }
        let a: Vec<_> = m
            .cells()
            .iter()
            .map(|c| Vec3::new(0.02 * c.center.y, -0.01 * c.center.x, 0.0))
            .collect();
        s.set_body_acceleration(&a).unwrap();
        for _ in 0..3 {
            let info = if mode == 0 || mode == 4 {
                s.iterate_simple().unwrap()
            } else if mode == 2 {
                s.step_pimple().unwrap()
            } else {
                s.step_piso().unwrap()
            };
            assert!(info.pressure.converged);
        }
        let mut out = vec![];
        vectors(&mut out, s.velocity());
        out.extend_from_slice(s.pressure());
        out.extend_from_slice(s.face_flux());
        out.extend([
            s.time(),
            s.kinetic_energy(),
            s.continuity_l2().unwrap(),
            s.continuity_max().unwrap(),
        ]);
        check(&g, &format!("flow_{mode}"), &out, 1e-8, 1e-8);
    }
}
#[test]
fn variable_property_reactive_radiative_heat_and_solid_conduction_match_cpp() {
    use std::sync::Arc;
    let g = fixture();
    let m = mesh(0.2);
    let mut s = ThermalTransport::new(
        m.clone(),
        ThermalTransportConfig {
            dt: 1e-5,
            density: 1.3,
            default_heat_capacity: 2.0,
            default_conductivity: 0.05,
        },
    )
    .unwrap();
    s.initialize_profile(|p| 300.0 + 0.3 * p.x + 0.2 * p.y)
        .unwrap();
    s.set_boundary(
        "left",
        ThermalBoundary {
            fixed_temperature: true,
            temperature: 310.0,
        },
    )
    .unwrap();
    s.set_boundary(
        "right",
        ThermalBoundary {
            fixed_temperature: true,
            temperature: 295.0,
        },
    )
    .unwrap();
    s.set_material(|_, t| 2.0 + 0.001 * (t - 300.0), |_, t| 0.05 + 0.0001 * t);
    s.set_face_mass_flux(
        m.faces()
            .iter()
            .map(|f| 1.3 * Vec3::new(0.03, 0.01, 0.0).dot(f.area))
            .collect(),
    )
    .unwrap();
    s.set_reactive_radiative_source(
        Some(Box::new(|p, time, t| {
            Ok(0.02 + 0.01 * p.x + 0.001 * time + 0.0001 * (t - 300.0))
        })),
        Arc::new(OpticallyThinRadiation::new(0.001).unwrap()),
        293.0,
    )
    .unwrap();
    s.run(7).unwrap();
    let mut out = s.temperature().to_vec();
    out.extend([s.time(), s.total_sensible_enthalpy().unwrap()]);
    vectors(
        &mut out,
        &boussinesq_acceleration(
            &m,
            s.temperature(),
            300.0,
            0.003,
            Vec3::new(0.0, -9.81, 0.0),
        )
        .unwrap(),
    );
    check(&g, "thermal", &out, 2e-10, 2e-12);
    let mut s = SolidHeatConduction::new(
        m,
        SolidHeatConfig {
            dt: 0.004,
            density: 2.0,
            heat_capacity: 1000.0,
            conductivity: 4.0,
            ..Default::default()
        },
    )
    .unwrap();
    s.initialize_profile(|p| 300.0 + 0.1 * p.x).unwrap();
    s.set_boundary("left", SolidThermalBoundaryType::FixedTemperature, 300.0)
        .unwrap();
    s.set_boundary("right", SolidThermalBoundaryType::FixedTemperature, 320.0)
        .unwrap();
    s.set_volumetric_heat_source_constant(200.0).unwrap();
    s.run(6).unwrap();
    let mut out = s.temperature().to_vec();
    out.extend([s.time(), s.total_sensible_energy()]);
    check(&g, "solid_heat", &out, 2e-9, 2e-11);
}
#[test]
fn turbulence_postprocessing_closures_and_radiosity_match_cpp() {
    let g = fixture();
    let m = mesh(0.2);
    let u: Vec<_> = m
        .cells()
        .iter()
        .map(|c| {
            Vec3::new(
                c.center.x * c.center.x + c.center.y,
                -2.0 * c.center.x + c.center.y * c.center.z,
                0.3 * c.center.z - c.center.y,
            )
        })
        .collect();
    for (mode, model) in [
        TurbulenceModel::Laminar,
        TurbulenceModel::MixingLength,
        TurbulenceModel::Smagorinsky,
        TurbulenceModel::Wale,
    ]
    .into_iter()
    .enumerate()
    {
        let out = eddy_viscosity(
            &m,
            &u,
            TurbulenceConfig {
                model,
                ..Default::default()
            },
        )
        .unwrap();
        check(&g, &format!("turbulence_{mode}"), &out, 2e-12, 2e-12);
    }
    let d = velocity_derived_fields(&m, &u).unwrap();
    let mut out = vec![];
    vectors(&mut out, &d.vorticity);
    out.extend(d.q_criterion);
    out.extend(d.lambda2);
    out.extend(d.enstrophy);
    out.push(wall_y_plus(0.001, 2.0, 1.0, 1e-5).unwrap());
    let force = pressure_force_on_patch(&m, &vec![2.0; m.cell_count()], 0).unwrap();
    out.extend(force.as_array());
    check(&g, "derived_fields", &out, 2e-10, 2e-10);
    let mut out = vec![
        spalart_allmaras_eddy_viscosity(1e-5, 2e-5, 7.1).unwrap(),
        k_epsilon_eddy_viscosity(2.0, 0.5, 0.1, 0.09).unwrap(),
        k_omega_sst_eddy_viscosity(1.2, 0.1, 3.0, 2.0, 1.0, 0.31).unwrap(),
        des_length_scale(0.1, 0.02, 0.65).unwrap(),
        drift_flux_slip_velocity(1000.0, 2000.0, 1e-4, 1e-3, 9.81).unwrap(),
        euler_euler_drag_source(1.0, 2.0, 3.0, 4.0).unwrap(),
        lubrication_thin_film_flux(0.001, 20.0, 0.01, 0.1).unwrap(),
    ];
    let p = linear_mushy_phase_change(305.0, 300.0, 310.0, 2.0, 4.0, 100.0).unwrap();
    out.extend([p.liquid_fraction, p.effective_cp, p.enthalpy]);
    let r = boussinesq_reynolds_stress(1.2, 0.5, 0.03, 0.1, -0.2, 0.4).unwrap();
    out.extend([r.xx, r.yy, r.xy]);
    check(&g, "advanced_models", &out, 2e-12, 2e-12);
    let mut out = vec![
        gray_surface_exchange(350.0, 0.7, 290.0, 0.8, 0.5).unwrap(),
        optically_thin_radiation_source(330.0, 295.0, 0.1).unwrap(),
    ];
    let r = solve_gray_radiosity(
        &[350.0, 300.0, 280.0],
        &[0.7, 0.8, 1.0],
        &[1.0; 3],
        &[0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0],
    )
    .unwrap();
    out.extend(r.radiosity);
    out.extend(r.net_heat);
    check(&g, "radiation", &out, 2e-9, 2e-11);
}
