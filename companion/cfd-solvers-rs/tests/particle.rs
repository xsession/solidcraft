//! Analytic and invariant checks for the particle module.
#![allow(clippy::needless_range_loop)]
use cfd_solvers::particle::*;
use std::f64::consts::PI;

#[test]
fn boris_gyration_conserves_speed_and_closes_orbit() {
    let (q, m, bz) = (1.0, 2.0, 0.5);
    let omega = q * bz / m;
    let period = 2.0 * PI / omega;
    let steps = 2000;
    let dt = period / steps as f64;
    let mut p = ChargedParticle {
        position_m: Vec3::new(0.0, 0.0, 0.0),
        velocity_m_per_s: Vec3::new(1.0, 0.0, 0.0),
        charge_c: q,
        mass_kg: m,
        weight: 1.0,
    };
    let f = ElectromagneticField {
        electric_v_per_m: Vec3::default(),
        magnetic_t: Vec3::new(0.0, 0.0, bz),
    };
    for _ in 0..steps {
        boris_push(&mut p, &f, dt).unwrap();
    }
    let v = p.velocity_m_per_s;
    assert!(((v.x * v.x + v.y * v.y).sqrt() - 1.0).abs() < 1e-12);
    // Boris rotates by 2·atan(ωdt/2) per step, slightly less than ωdt.
    assert!((v.x - 1.0).abs() < 1e-4 && v.y.abs() < 1e-2, "{v:?}");
    assert!(p.position_m.x.abs() < 1e-2 && p.position_m.y.abs() < 1e-2);
}

#[test]
fn spectral_poisson_satisfies_gauss_law() {
    let (n, l, eps_r) = (64, 2.0, 1.0);
    let k = 2.0 * PI * 3.0 / l;
    let rho: Vec<f64> = (0..n)
        .map(|i| (k * l * i as f64 / n as f64).sin())
        .collect();
    let e = periodic_electric_field_from_charge_density(&rho, l, eps_r, true).unwrap();
    let eps = 8.8541878128e-12 * eps_r;
    for i in 0..n {
        let x = l * i as f64 / n as f64;
        let exact = -(k * x).cos() / (k * eps);
        assert!(
            (e[i] - exact).abs() < 1e-9 * exact.abs().max(1.0 / (k * eps)),
            "{i}: {} vs {exact}",
            e[i]
        );
    }
}

#[test]
fn charge_conserving_current_satisfies_continuity() {
    let old: Vec<PicParticle3D> = (0..9)
        .map(|i| PicParticle3D {
            position_m: Vec3::new(0.1 * i as f64, 0.05 + 0.1 * i as f64, 0.9 - 0.08 * i as f64),
            charge_c: 1.0 + i as f64,
            ..Default::default()
        })
        .collect();
    let mut new = old.clone();
    for p in &mut new {
        p.position_m.x += 0.04;
        p.position_m.z -= 0.03;
    }
    for local in [false, true] {
        let cc = if local {
            deposit_charge_conserving_current_3d_local(&old, &new, 6, 5, 4, 1.0, 1.0, 1.0, 1.0e-3)
        } else {
            deposit_charge_conserving_current_3d(&old, &new, 6, 5, 4, 1.0, 1.0, 1.0, 1.0e-3)
        }
        .unwrap();
        let scale = cc
            .charge_density_new_c_per_m3
            .iter()
            .fold(0.0f64, |a, b| a.max(b.abs()))
            / 1.0e-3;
        assert!(
            cc.continuity_linf_residual <= 1e-9 * scale,
            "local={local}: {}",
            cc.continuity_linf_residual
        );
    }
}

#[test]
fn staggered_vacuum_mode_keeps_energy() {
    let c = StaggeredElectromagneticPic3DConfig {
        nx: 8,
        ny: 8,
        nz: 8,
        dt_s: 1.5e-10,
        ..Default::default()
    };
    let mut st = StaggeredElectromagneticPic3D::new(c).unwrap();
    st.initialize_z_polarized_mode(10.0, 1, 1, 0).unwrap();
    let e0 = st.diagnostics().unwrap().field_energy_j;
    st.step(40, None).unwrap();
    let e1 = st.diagnostics().unwrap().field_energy_j;
    assert!((e1 - e0).abs() < 0.1 * e0, "{e0} -> {e1}");
    // a 0.98 CFL violation is an error, not a panic
    let bad = StaggeredElectromagneticPic3DConfig {
        nx: 12,
        ny: 12,
        nz: 12,
        dt_s: 1.0e-8,
        ..Default::default()
    };
    assert!(StaggeredElectromagneticPic3D::new(bad).is_err());
    assert!(st.initialize_z_polarized_mode(1.0, 0, 0, 0).is_err());
}

#[test]
fn sponge_damps_only_near_faces() {
    let n = 7 * 7 * 7;
    let mut f: Vec<Vec<f64>> = (0..6).map(|_| vec![1.0; n]).collect();
    let [ex, ey, ez, bx, by, bz] = &mut f[..] else {
        unreachable!()
    };
    let bc = GridBoundary3DConfig {
        mode: GridBoundaryMode3D::AbsorbingSponge,
        ..Default::default()
    };
    apply_electromagnetic_field_boundary_3d(ex, ey, ez, bx, by, bz, 7, 7, 7, &bc).unwrap();
    let centre = (3 * 7 + 3) * 7 + 3;
    assert_eq!(ex[centre], 1.0);
    assert!((ex[0] - (-3.0f64).exp()).abs() < 1e-15);
    assert!(apply_electromagnetic_field_boundary_3d(
        ex,
        ey,
        ez,
        bx,
        by,
        &mut [0.0; 3],
        7,
        7,
        7,
        &bc
    )
    .is_err());
}

#[test]
fn slab_transmission_matches_exponential_attenuation() {
    let slab = SlabTransport1D::new(SlabTransportConfig {
        width: 1.0,
        total_cross_section: 1.2,
        absorption_cross_section: 1.2,
    })
    .unwrap();
    let t = slab.run(40_000, 5).unwrap();
    let transmitted = t.transmitted as f64 / t.histories as f64;
    assert!(
        (transmitted - (-1.2f64).exp()).abs() < 0.01,
        "{transmitted}"
    );
    assert_eq!(t.reflected, 0);
}

fn slab_world() -> TransportWorld {
    let bx = |a: f64, b: f64| AxisAlignedParticleBox {
        minimum_m: Vec3::new(a, 0.0, 0.0),
        maximum_m: Vec3::new(b, 1.0, 1.0),
    };
    TransportWorld {
        particles: vec![TransportParticleDefinition::default()],
        materials: vec![
            TransportMaterial::default(),
            TransportMaterial {
                name: "dense".into(),
                density_kg_per_m3: 5000.0,
                production_cut_energy_ev: 0.0,
            },
        ],
        regions: vec![
            TransportRegion {
                name: "a".into(),
                r#box: bx(0.0, 0.2),
                material_index: 0,
                sensitive_detector_index: Some(0),
            },
            TransportRegion {
                name: "b".into(),
                r#box: bx(0.2, 0.5),
                material_index: 1,
                sensitive_detector_index: Some(0),
            },
        ],
        processes: vec![
            TransportProcess {
                name: "dedx".into(),
                stopping_power_ev_per_m: 1.0e5,
                ..Default::default()
            },
            TransportProcess {
                name: "hard".into(),
                kind: TransportProcessKind::DiscreteInteraction,
                material_index: Some(1),
                physical_interaction_length_m: 0.07,
                energy_loss_fraction: 0.05,
                secondary_energy_fraction: 0.1,
                ..Default::default()
            },
        ],
        sensitive_detectors: vec![TransportSensitiveDetector::default()],
    }
}

#[test]
fn transport_conserves_energy_and_scores_dose() {
    let world = slab_world();
    let t = TransportTrack {
        position_m: Vec3::new(0.0, 0.5, 0.5),
        direction: Vec3::new(1.0, 0.0, 0.0),
        kinetic_energy_ev: 1.0e5,
        ..Default::default()
    };
    let r = transport_track(
        &world,
        t,
        &TransportConfig {
            max_step_m: 0.09,
            ..Default::default()
        },
        500,
    )
    .unwrap();
    assert_eq!(
        r.primary.status,
        TrackStatus::Escaped,
        "{:?} {:?}",
        r.primary,
        r.steps
            .iter()
            .map(|s| (s.step_length_m, s.limiter))
            .collect::<Vec<_>>()
    );
    assert!(!r.secondaries.is_empty());
    let secondary: f64 = r
        .secondaries
        .iter()
        .map(|s| s.track.kinetic_energy_ev)
        .sum();
    let total = r.scoring.total_energy_deposit_ev + r.primary.kinetic_energy_ev + secondary;
    assert!((total - 1.0e5).abs() < 1e-6, "{total}");
    assert!(
        (r.scoring.total_track_length_m - 0.5).abs() < 1e-9,
        "{}",
        r.scoring.total_track_length_m
    );

    let hits = collect_transport_hits(&world, &r.steps).unwrap();
    assert!((hits.total_energy_deposit_ev - r.scoring.total_energy_deposit_ev).abs() < 1e-6);
    let b = AxisAlignedParticleBox {
        minimum_m: Vec3::default(),
        maximum_m: Vec3::new(0.5, 1.0, 1.0),
    };
    let mut grid = make_transport_dose_grid(5, 1, 1, b, 1000.0).unwrap();
    score_hits_to_dose_grid(&hits, &mut grid).unwrap();
    let voxel_mass = 1000.0 * 0.1;
    let energy_j: f64 = grid.dose_gy.iter().map(|d| d * voxel_mass).sum();
    assert!((energy_j - hits.total_energy_deposit_ev * 1.602176634e-19).abs() < 1e-24);
    let sar = project_dose_grid_to_pennes2d_sar(&grid, 2.0, 1, 1).unwrap();
    let mean: f64 = grid.dose_gy.iter().sum::<f64>() / 5.0 / 2.0;
    assert!((sar[0] - mean).abs() < 1e-30 + 1e-12 * mean);
}

#[test]
fn bvh_agrees_with_linear_search_and_errors_are_values() {
    let mut world = slab_world();
    world.regions.clear();
    for i in 0..5 {
        for j in 0..4 {
            for k in 0..3 {
                let (x, y, z) = (i as f64, j as f64, k as f64);
                world.regions.push(TransportRegion {
                    name: format!("r{i}{j}{k}"),
                    r#box: AxisAlignedParticleBox {
                        minimum_m: Vec3::new(x, y, z),
                        maximum_m: Vec3::new(x + 1.0, y + 1.0, z + 1.0),
                    },
                    material_index: 0,
                    sensitive_detector_index: None,
                });
            }
        }
    }
    let bvh = build_region_bvh(&world, 3).unwrap();
    let mut rng = TransportRandom::default();
    for _ in 0..500 {
        let p = Vec3::new(
            5.5 * rng.uniform_open01() - 0.2,
            4.4 * rng.uniform_open01(),
            3.3 * rng.uniform_open01(),
        );
        let a = locate_region(&world, p, 0.0);
        let b = locate_region_bvh(&world, &bvh, p, 0.0);
        // Regions only share faces, so both must agree on containment.
        assert_eq!(a.is_some(), b.is_some());
        if let Some(b) = b {
            assert!(contains(&world.regions[b].r#box, p, 0.0));
        }
    }
    assert!(build_region_bvh(&world, 0).is_err());
    assert!(normalize_direction(Vec3::default()).is_err());
    world.regions[0].material_index = 9;
    let t = TransportTrack {
        direction: Vec3::new(1.0, 0.0, 0.0),
        ..Default::default()
    };
    assert!(transport_track(&world, t, &TransportConfig::default(), 5).is_err());
}
