use cfd_solvers::{core::*, fvm::*};
use serde_json::Value;
fn mesh() -> PolyMesh {
    make_sheared_cartesian_hexa_mesh(4, 3, 2, 1.2, 0.8, 0.5, 0.2).unwrap()
}
fn check(name: &str, v: &[f64]) {
    let g: Value = serde_json::from_str(include_str!("../parity/goldens-fvm.json")).unwrap();
    let e = g[name].as_array().unwrap();
    assert_eq!(v.len(), e.len(), "{name}");
    for (i, (x, y)) in v.iter().zip(e).enumerate() {
        let y = y.as_f64().unwrap();
        assert!(
            x.is_finite() && (x - y).abs() <= 2e-11 + 2e-9 * y.abs(),
            "{name}[{i}]: {x} vs {y}"
        );
    }
}
fn controls(mode: usize) -> RansTransportControls {
    RansTransportControls {
        dt: 0.001,
        linear: LinearConfig {
            max_iterations: 400,
            restart: 30,
            relative_tolerance: 1e-12,
        },
        temporal_scheme: [
            TemporalScheme::Euler,
            TemporalScheme::BackwardBdf2,
            TemporalScheme::CrankNicolson,
        ][mode],
        ..Default::default()
    }
}
fn fields(
    m: &PolyMesh,
) -> (
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<Vec3>,
) {
    let flux = m
        .faces()
        .iter()
        .map(|f| Vec3::new(0.03, 0.01, 0.0).dot(f.area))
        .collect();
    let nt = m.cells().iter().map(|c| 2e-5 + 1e-6 * c.center.x).collect();
    let k = m
        .cells()
        .iter()
        .map(|c| 0.02 + 0.001 * c.center.y)
        .collect();
    let e = m
        .cells()
        .iter()
        .map(|c| 0.01 + 0.001 * c.center.z)
        .collect();
    let w = m.cells().iter().map(|c| 0.6 + 0.03 * c.center.x).collect();
    let wall = m.cells().iter().map(|c| 0.05 + 0.01 * c.center.y).collect();
    let u = m
        .cells()
        .iter()
        .map(|c| Vec3::new(0.1 * c.center.y, -0.03 * c.center.x, 0.02 * c.center.z))
        .collect();
    (flux, nt, k, e, w, wall, u)
}
#[test]
fn sa_and_kepsilon_all_time_schemes_match_cpp() {
    for mode in 0..3 {
        let m = mesh();
        let (flux, nt, k, e, _, wall, u) = fields(&m);
        let mut sa = SpalartAllmarasTransport::new(
            m.clone(),
            SpalartAllmarasConfig {
                transport: controls(mode),
                ..Default::default()
            },
        )
        .unwrap();
        sa.initialize_profile(&nt).unwrap();
        sa.set_wall_distance_profile(&wall).unwrap();
        sa.set_velocity(&u).unwrap();
        sa.set_face_flux(flux.clone()).unwrap();
        sa.set_boundary("left", ScalarBoundaryType::FixedValue, 2e-5)
            .unwrap();
        assert!(sa.run(4).unwrap().unwrap().converged);
        let mut out = sa.nu_tilde().to_vec();
        out.extend(sa.kinematic_eddy_viscosity().unwrap());
        out.extend_from_slice(sa.strain_rate());
        out.push(sa.time());
        check(&format!("rans_sa_{mode}"), &out);
        let mut s = KEpsilonTransport::new(
            m,
            KEpsilonConfig {
                transport: controls(mode),
                density: 1.2,
                ..Default::default()
            },
        )
        .unwrap();
        s.initialize_profile(&k, &e).unwrap();
        s.set_velocity(&u).unwrap();
        s.set_face_flux(flux).unwrap();
        s.set_k_boundary("left", ScalarBoundaryType::FixedValue, 0.02)
            .unwrap();
        s.set_epsilon_boundary("left", ScalarBoundaryType::FixedValue, 0.01)
            .unwrap();
        let r = s.run(4).unwrap().unwrap();
        assert!(r.first.converged && r.second.converged);
        let mut out = s.k().to_vec();
        out.extend_from_slice(s.epsilon());
        out.extend(s.kinematic_eddy_viscosity().unwrap());
        out.extend(s.dynamic_eddy_viscosity().unwrap());
        out.push(s.time());
        check(&format!("rans_ke_{mode}"), &out);
    }
}
#[test]
fn sst_time_schemes_and_three_des_filters_match_cpp() {
    for mode in 0..3 {
        let m = mesh();
        let (flux, _, k, _, w, wall, u) = fields(&m);
        let mut s = KOmegaSSTTransport::new(
            m,
            KOmegaSSTConfig {
                transport: controls(mode),
                density: 1.2,
                ..Default::default()
            },
        )
        .unwrap();
        s.initialize_profile(&k, &w).unwrap();
        s.set_wall_distance_profile(&wall).unwrap();
        s.set_velocity(&u).unwrap();
        s.set_face_flux(flux).unwrap();
        s.set_k_boundary("left", ScalarBoundaryType::FixedValue, 0.02)
            .unwrap();
        s.set_omega_boundary("left", ScalarBoundaryType::FixedValue, 0.6)
            .unwrap();
        s.run(4).unwrap();
        let mut out = s.k().to_vec();
        out.extend_from_slice(s.omega());
        out.extend(s.kinematic_eddy_viscosity().unwrap());
        out.extend_from_slice(s.blending_f1());
        out.extend_from_slice(s.blending_f2());
        out.extend_from_slice(s.hybrid_dissipation_factor());
        out.push(s.time());
        check(&format!("rans_sst_{mode}"), &out);
    }
    for (i, filter) in [
        SSTDESZonalFilter::None,
        SSTDESZonalFilter::F1,
        SSTDESZonalFilter::F2,
    ]
    .into_iter()
    .enumerate()
    {
        let mut s = KOmegaSSTTransport::new(
            mesh(),
            KOmegaSSTConfig {
                transport: controls(0),
                des_enabled: true,
                des_zonal_filter: filter,
                ..Default::default()
            },
        )
        .unwrap();
        s.initialize(0.05, 0.2).unwrap();
        s.set_wall_distance(10.0).unwrap();
        s.set_grid_scale(0.01).unwrap();
        s.set_strain_rate(0.3).unwrap();
        s.run(3).unwrap();
        let mut out = s.k().to_vec();
        out.extend_from_slice(s.omega());
        out.extend_from_slice(s.blending_f1());
        out.extend_from_slice(s.blending_f2());
        out.extend_from_slice(s.hybrid_dissipation_factor());
        assert!(s.hybrid_dissipation_factor().iter().all(|v| *v > 1.0));
        check(&format!("rans_des_{i}"), &out);
    }
}
#[test]
fn reynolds_stress_transport_all_time_schemes_match_cpp() {
    for mode in 0..3 {
        let m = mesh();
        let (flux, _, k, e, _, _, u) = fields(&m);
        let mut s = ReynoldsStressTransport::new(
            m,
            ReynoldsStressConfig {
                transport: controls(mode),
                ..Default::default()
            },
        )
        .unwrap();
        let r: Vec<_> = k
            .iter()
            .map(|&k| SymmetricTensor3 {
                xx: 0.8 * k,
                yy: 0.7 * k,
                zz: 0.5 * k,
                xy: 0.02 * k,
                xz: -0.01 * k,
                yz: 0.03 * k,
            })
            .collect();
        s.initialize_profile(&r, &e).unwrap();
        s.set_velocity(&u).unwrap();
        s.set_face_flux(flux).unwrap();
        s.set_additional_source(Some(Box::new(|p, t| {
            Ok(SymmetricTensor3 {
                xx: 0.0001 * p.x,
                xy: 0.00001 + t * 0.001,
                ..Default::default()
            })
        })));
        let r = s.run(4).unwrap().unwrap();
        assert!(r.stress.iter().all(|r| r.converged) && r.epsilon.converged);
        let stresses = s.reynolds_stress();
        assert!(stresses
            .iter()
            .all(|r| reynolds_stress_is_realizable(*r, 1e-12)));
        let mut out: Vec<_> = stresses.iter().flat_map(|r| r.components()).collect();
        out.extend_from_slice(s.epsilon());
        out.extend(s.turbulent_kinetic_energy());
        out.extend(s.kinematic_eddy_viscosity().unwrap());
        out.push(s.time());
        check(&format!("rans_rsm_{mode}"), &out);
    }
}
#[test]
fn homogeneous_decay_is_positive_and_matches_the_discrete_equations() {
    let m = make_cartesian_hexa_mesh(1, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let c = KEpsilonConfig {
        transport: RansTransportControls {
            dt: 0.01,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut s = KEpsilonTransport::new(m.clone(), c).unwrap();
    s.initialize(0.1, 0.02).unwrap();
    s.step().unwrap();
    assert!((s.k()[0] - 0.1 / (1.0 + 0.01 * 0.02 / 0.1)).abs() < 1e-13);
    assert!((s.epsilon()[0] - 0.02 / (1.0 + 0.01 * 1.92 * 0.02 / 0.1)).abs() < 1e-13);
    s.run(200).unwrap();
    assert!(s.k()[0] > 0.0 && s.k()[0] < 0.1 && s.epsilon()[0] < 0.02);
    let mut r = ReynoldsStressTransport::new(m, ReynoldsStressConfig::default()).unwrap();
    r.initialize_isotropic(0.1, 0.02).unwrap();
    r.run(3).unwrap();
    let v = r.reynolds_stress()[0];
    assert!((v.xx - v.yy).abs() < 1e-14 && (v.yy - v.zz).abs() < 1e-14 && v.xy.abs() < 1e-14);
    assert!(turbulent_kinetic_energy(v) < 0.1);
}
#[test]
fn invalid_inputs_and_failed_multiequation_steps_do_not_commit() {
    let m = mesh();
    assert!(SpalartAllmarasTransport::new(
        m.clone(),
        SpalartAllmarasConfig {
            molecular_viscosity: f64::INFINITY,
            ..Default::default()
        }
    )
    .is_err());
    let mut s = KOmegaSSTTransport::new(m.clone(), Default::default()).unwrap();
    assert!(s.set_grid_scale(0.0).is_err());
    assert!(s.set_wall_distance(-1.0).is_err());
    assert!(s.set_face_flux(vec![0.0]).is_err());
    let mut r = ReynoldsStressTransport::new(m.clone(), Default::default()).unwrap();
    assert!(r
        .initialize_profile(
            &vec![
                SymmetricTensor3 {
                    xx: 1.0,
                    yy: 1.0,
                    zz: 1.0,
                    xy: 2.0,
                    ..Default::default()
                };
                m.cell_count()
            ],
            &vec![1.0; m.cell_count()]
        )
        .is_err());
    r.initialize_isotropic(0.1, 0.02).unwrap();
    let before = r.reynolds_stress()[0].components();
    r.set_additional_source(Some(Box::new(|_, _| {
        Ok(SymmetricTensor3 {
            yz: f64::NAN,
            ..Default::default()
        })
    })));
    assert!(r.step().is_err());
    assert_eq!(r.time(), 0.0);
    assert_eq!(r.reynolds_stress()[0].components(), before);
    assert!(r.linear_result().is_none());
}
