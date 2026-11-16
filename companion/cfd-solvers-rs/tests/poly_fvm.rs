use cfd_solvers::{core::Vec3, fvm::*};

#[test]
fn limited_schemes_are_bounded_conservative_and_reuse_workspace() {
    let m = make_cartesian_hexa_mesh(20, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let v: Vec<_> = m
        .cells()
        .iter()
        .map(|c| if c.center.x < 0.5 { 0.0 } else { 1.0 })
        .collect();
    let velocity = vec![Vec3::new(0.2, 0.0, 0.0); m.cell_count()];
    let bvec = vec![Vec3::new(0.2, 0.0, 0.0); m.face_count()];
    let flux = face_flux_from_velocity(&m, &velocity, &bvec).unwrap();
    let boundary: Vec<_> = m
        .faces()
        .iter()
        .map(|f| if f.center.x < 0.5 { 0.0 } else { 1.0 })
        .collect();
    for scheme in [
        FaceInterpolationScheme::Upwind,
        FaceInterpolationScheme::BoundedLinear,
        FaceInterpolationScheme::Minmod,
        FaceInterpolationScheme::VanLeer,
    ] {
        let faces = interpolate_scalar_to_faces(&m, &v, scheme, &flux, &boundary).unwrap();
        assert!(faces.iter().all(|v| *v >= 0.0 && *v <= 1.0));
        let mut workspace = ScalarSchemeWorkspace::default();
        let mut divergence = vec![0.0; m.cell_count()];
        convective_divergence_scalar_into(
            &m,
            &v,
            &flux,
            &mut divergence,
            &mut workspace,
            scheme,
            &boundary,
        )
        .unwrap();
        let bytes = workspace.allocated_bytes();
        convective_divergence_scalar_into(
            &m,
            &v,
            &flux,
            &mut divergence,
            &mut workspace,
            scheme,
            &boundary,
        )
        .unwrap();
        assert_eq!(bytes, workspace.allocated_bytes());
        let integrated: f64 = divergence
            .iter()
            .zip(m.cells())
            .map(|(v, c)| v * c.volume)
            .sum();
        assert!((integrated - 0.2).abs() < 1e-13);
    }
    let dt = courant_limited_timestep(&m, &flux, 0.5, 10.0).unwrap();
    assert!((dt - 0.0625).abs() < 1e-14);
}

#[test]
fn scalar_temporal_schemes_match_the_discrete_diffusion_mode() {
    let n = 16;
    let dx = 1.0 / n as f64;
    let diffusivity = 0.03;
    let dt = 0.02;
    let lambda =
        4.0 * diffusivity * (std::f64::consts::PI / (2.0 * n as f64)).sin().powi(2) / (dx * dx);
    for scheme in [
        TemporalScheme::Euler,
        TemporalScheme::BackwardBdf2,
        TemporalScheme::CrankNicolson,
    ] {
        let m = make_cartesian_hexa_mesh(n, 1, 1, 1.0, 1.0, 1.0).unwrap();
        let mut s = ScalarTransport::new(
            m,
            ScalarTransportConfig {
                dt,
                diffusivity,
                temporal_scheme: scheme,
                ..Default::default()
            },
        )
        .unwrap();
        s.initialize_profile(|p| (std::f64::consts::PI * p.x).cos())
            .unwrap();
        let mut amplitude = 1.0;
        let mut older = 1.0;
        for step in 0..8 {
            let next = match scheme {
                TemporalScheme::Euler => amplitude / (1.0 + dt * lambda),
                TemporalScheme::CrankNicolson => {
                    amplitude * (1.0 - 0.5 * dt * lambda) / (1.0 + 0.5 * dt * lambda)
                }
                TemporalScheme::BackwardBdf2 => {
                    if step == 0 {
                        amplitude / (1.0 + dt * lambda)
                    } else {
                        (2.0 * amplitude - 0.5 * older) / (1.5 + dt * lambda)
                    }
                }
            };
            older = amplitude;
            amplitude = next;
            s.step().unwrap();
        }
        for (c, v) in s.mesh().cells().iter().zip(s.values()) {
            assert!((v - amplitude * (std::f64::consts::PI * c.center.x).cos()).abs() < 1e-8);
        }
        assert!(s.volume_integral().abs() < 1e-12);
        assert_eq!(s.operator_assemblies(), 1);
        assert!((s.time() - 8.0 * dt).abs() < 1e-14);
        s.set_boundary("left", ScalarBoundaryType::FixedValue, 1.0)
            .unwrap();
        s.step().unwrap();
        assert_eq!(s.operator_assemblies(), 2);
    }
}

#[test]
fn scalar_sources_and_fixed_inflow_preserve_constant_solutions() {
    for scheme in [
        TemporalScheme::Euler,
        TemporalScheme::BackwardBdf2,
        TemporalScheme::CrankNicolson,
    ] {
        let m = make_cartesian_hexa_mesh(8, 2, 1, 1.0, 0.5, 0.1).unwrap();
        let mut s = ScalarTransport::new(
            m,
            ScalarTransportConfig {
                dt: 0.01,
                temporal_scheme: scheme,
                ..Default::default()
            },
        )
        .unwrap();
        s.initialize_constant(2.0).unwrap();
        s.set_source_constant(4.0).unwrap();
        s.run(5).unwrap();
        assert!(s.values().iter().all(|v| (v - 2.2).abs() < 1e-10));
    }
}

#[test]
fn collocated_channel_piso_conserves_volume_flux() {
    let m = make_cartesian_hexa_mesh(8, 4, 1, 1.0, 0.5, 0.1).unwrap();
    let mut s = CollocatedIncompressible::new(
        m,
        CollocatedIncompressibleConfig {
            include_convection: false,
            pressure_tolerance: 1e-12,
            pressure_preconditioner: PressurePreconditionerKind::Jacobi,
            ..Default::default()
        },
    )
    .unwrap();
    s.set_velocity_boundary(
        "left",
        VelocityBoundaryType::FixedValue,
        Vec3::new(0.2, 0.0, 0.0),
    )
    .unwrap();
    s.set_pressure_boundary("right", PressureBoundaryType::FixedValue, 0.0)
        .unwrap();
    for name in ["bottom", "top", "front", "back"] {
        s.set_velocity_boundary(name, VelocityBoundaryType::Slip, Vec3::ZERO)
            .unwrap();
    }
    s.initialize_uniform(Vec3::new(0.2, 0.0, 0.0), 0.0).unwrap();
    for _ in 0..3 {
        let info = s.step_piso().unwrap();
        assert!(info.pressure.converged);
        assert!(info.continuity_l2 < 1e-8, "{}", info.continuity_l2);
    }
    assert!(s
        .velocity()
        .iter()
        .all(|u| (u.x - 0.2).abs() < 1e-8 && u.y.abs() < 1e-8));
    assert_eq!(s.steps(), 3);
    assert!((s.time() - 0.03).abs() < 1e-14);
}

#[test]
fn collocated_lid_flow_has_a_pressure_gauge_and_real_correctors() {
    let m = make_sheared_cartesian_hexa_mesh(5, 5, 1, 1.0, 1.0, 0.1, 0.2).unwrap();
    let mut s = CollocatedIncompressible::new(
        m,
        CollocatedIncompressibleConfig {
            include_convection: false,
            nonorthogonal_correctors: 3,
            pressure_preconditioner: PressurePreconditionerKind::Jacobi,
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
    for _ in 0..4 {
        let info = s.step_pimple().unwrap();
        assert!(info.pressure.converged);
        assert!(info.pressure_solves >= s.config().pressure_correctors);
        assert!(info.continuity_l2 < 1e-3, "{}", info.continuity_l2);
    }
    assert!(s.kinetic_energy() > 0.0);
    assert!(s.pressure().iter().sum::<f64>().abs() < 1e-10);
    assert!(s.velocity().iter().all(|u| u.is_finite()));
}

#[test]
fn external_pressure_preconditioner_must_be_supplied_and_unsupported_knobs_fail() {
    let m = make_cartesian_hexa_mesh(3, 2, 1, 1.0, 1.0, 1.0).unwrap();
    assert!(CollocatedIncompressible::new(
        m.clone(),
        CollocatedIncompressibleConfig {
            convection_scheme: FaceInterpolationScheme::VanLeer,
            ..Default::default()
        }
    )
    .is_err());
    let mut s = CollocatedIncompressible::new(
        m,
        CollocatedIncompressibleConfig {
            pressure_preconditioner: PressurePreconditionerKind::ExternalAmg,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(s.step_piso().is_err());
    s.set_pressure_amg_cycle(|r, z| {
        z.copy_from_slice(r);
        Ok(())
    });
    let result = s.step_piso().unwrap();
    assert!(result.pressure.converged);
}
fn linear(p: Vec3) -> f64 {
    0.3 + 1.2 * p.x - 0.7 * p.y + 0.4 * p.z
}
#[test]
fn sheared_poly_mesh_has_closed_faces_positive_metrics_and_round_trips() {
    let m = make_sheared_cartesian_hexa_mesh(4, 3, 2, 1.2, 0.8, 0.5, 0.4).unwrap();
    assert_eq!(m.cell_count(), 24);
    assert_eq!(m.face_count(), (5 * 3 * 2) + (4 * 4 * 2) + (4 * 3 * 3));
    assert_eq!(m.patch_index("back").unwrap(), 5);
    for c in 0..m.cell_count() {
        let sum = m.cell_faces()[c]
            .iter()
            .map(|&fi| m.faces()[fi].outward_area(c))
            .fold(Vec3::ZERO, |a, b| a + b);
        assert!(sum.norm() < 1e-14);
    }
    for fi in 0..m.face_count() {
        let d = decompose_face_area(&m, fi).unwrap();
        assert!(d.orthogonal_metric > 0.0);
        assert!((d.orthogonal_area + d.nonorthogonal_area - m.faces()[fi].area).norm() < 1e-14);
    }
    let text = serde_json::to_string(&m).unwrap();
    let decoded: PolyMesh = serde_json::from_str(&text).unwrap();
    assert_eq!(decoded.cell_faces(), m.cell_faces());
    assert!(serde_json::from_str::<PolyMesh>(
        r#"{"cells":[{"center":{"x":0,"y":0,"z":0},"volume":-1}],"faces":[],"patches":[]}"#
    )
    .is_err());
}
#[test]
fn poly_operators_reproduce_affine_fields_and_conserve_flux() {
    let m = make_sheared_cartesian_hexa_mesh(4, 3, 2, 1.2, 0.8, 0.5, 0.4).unwrap();
    let p: Vec<_> = m.cells().iter().map(|c| linear(c.center)).collect();
    let b: Vec<_> = m.faces().iter().map(|f| linear(f.center)).collect();
    for g in [
        gauss_gradient_scalar(&m, &p, &b).unwrap(),
        least_squares_gradient_scalar(&m, &p, &b).unwrap(),
    ] {
        for v in g {
            assert!((v - Vec3::new(1.2, -0.7, 0.4)).norm() < 2e-12);
        }
    }
    let lap = corrected_laplacian_scalar(&m, &p, 0.3, &b).unwrap();
    assert!(lap.iter().all(|v| v.abs() < 1e-10));
    let u = |p: Vec3| Vec3::new(p.x, 2.0 * p.y, -0.5 * p.z);
    let cells: Vec<_> = m.cells().iter().map(|c| u(c.center)).collect();
    let boundary: Vec<_> = m.faces().iter().map(|f| u(f.center)).collect();
    let div = gauss_divergence_vector(&m, &cells, &boundary).unwrap();
    assert!(div.iter().all(|v| (v - 2.5).abs() < 1e-12));
}
#[test]
fn slip_boundaries_and_rhie_chow_are_pressure_sensitive() {
    let m = make_cartesian_hexa_mesh(4, 3, 2, 1.0, 1.0, 1.0).unwrap();
    let n = m.cell_count();
    let u = vec![Vec3::new(0.2, 0.3, 0.1); n];
    let mut ub = vec![VelocityBoundaryCondition::default(); 6];
    ub[0] = VelocityBoundaryCondition {
        boundary_type: VelocityBoundaryType::Slip,
        value: Vec3::new(0.0, 0.1, 0.0),
    };
    let v = boundary_velocity_values(&m, &u, &ub).unwrap();
    for (i, f) in m.faces().iter().enumerate() {
        if f.patch == Some(0) {
            assert!(v[i].x.abs() < 1e-14);
            assert!((v[i].y - 0.3).abs() < 1e-14);
        }
    }
    let pressure = vec![PressureBoundaryCondition::default(); 6];
    let base =
        rhie_chow_face_flux(&m, &u, &vec![0.1; n], &vec![0.0; n], &ub, &pressure, true).unwrap();
    let p: Vec<_> = m.cells().iter().map(|c| c.center.x).collect();
    let changed = rhie_chow_face_flux(&m, &u, &vec![0.1; n], &p, &ub, &pressure, true).unwrap();
    assert!(changed.iter().zip(base).any(|(a, b)| (a - b).abs() > 1e-3));
    let flux: Vec<_> = m
        .faces()
        .iter()
        .map(|f| Vec3::new(0.2, 0.3, 0.1).dot(f.area))
        .collect();
    assert!(flux_divergence_l2(&m, &flux).unwrap() < 1e-13);
}
