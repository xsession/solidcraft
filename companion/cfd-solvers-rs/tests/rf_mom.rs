//! Analytic checks for thin-wire/surface MoM, antennas and asymptotic methods.
use cfd_solvers::{core::Vec3, optics::Ray, rf::*, Complex64 as C};
use std::f64::consts::PI;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

#[test]
fn half_wave_dipole_impedance_and_reciprocity() {
    let c0 = 299792458.0;
    let f = 300e6;
    let r = solve_center_fed_thin_wire(&ThinWireMomConfig {
        length_m: 0.48 * c0 / f,
        radius_m: 1e-3,
        frequency_hz: f,
        segments: 31,
        ..Default::default()
    })
    .unwrap();
    // Upstream's pulse-basis/point-matched Pocklington reference is only claimed to give a finite,
    // passive, mirror-symmetric feed solution (its reactance is far from the textbook value).
    assert!(
        r.feed_impedance_ohm.re > 0.0
            && r.feed_impedance_ohm.im.is_finite()
            && r.feed_current_a > 0.0
    );
    // Current symmetric about the feed.
    let n = r.current_a.len();
    for i in 0..n / 2 {
        assert!((r.current_a[i] - r.current_a[n - 1 - i]).norm() < 1e-9 * r.feed_current_a);
    }
    let a = OrientedThinWire {
        start_m: WirePoint3::new(0.0, 0.0, -0.24),
        end_m: WirePoint3::new(0.0, 0.0, 0.24),
        radius_m: 2e-3,
        segments: 11,
        ..Default::default()
    };
    // Point matching is not Galerkin, so reciprocity holds only for congruent ports.
    let b = OrientedThinWire {
        start_m: WirePoint3::new(0.3, 0.1, -0.24),
        end_m: WirePoint3::new(0.3, 0.1, 0.24),
        ..a
    };
    let cfg = OrientedWireMomConfig {
        frequency_hz: 300e6,
        wires: vec![a, b],
        feeds: vec![
            ThinWireDeltaGapFeed {
                wire: 0,
                ..Default::default()
            },
            ThinWireDeltaGapFeed {
                wire: 1,
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    let z = oriented_wire_port_impedance_matrix(&cfg).unwrap();
    assert!((z[(0, 1)] - z[(1, 0)]).norm() < 1e-6 * z[(0, 1)].norm());
}

#[test]
fn lossless_wire_radiates_everything_and_losses_balance() {
    let wire = OrientedThinWire {
        start_m: WirePoint3::new(0.0, 0.0, -0.25),
        end_m: WirePoint3::new(0.0, 0.0, 0.25),
        radius_m: 1e-3,
        segments: 15,
        ..Default::default()
    };
    let mut cfg = OrientedWireMomConfig {
        frequency_hz: 290e6,
        wires: vec![wire],
        feeds: vec![ThinWireDeltaGapFeed::default()],
        ..Default::default()
    };
    let r = solve_oriented_thin_wires(&cfg).unwrap();
    assert!((r.radiation_efficiency - 1.0).abs() < 1e-12 && r.radiated_power_w > 0.0);
    cfg.wires[0].conductivity_s_per_m = 1e5;
    cfg.loads.push(ThinWireSeriesLoad {
        wire: 0,
        segment: 2,
        impedance_ohm: C::new(10.0, 0.0),
    });
    let l = solve_oriented_thin_wires(&cfg).unwrap();
    let total = l.conductor_loss_w + l.lumped_load_loss_w + l.radiated_power_w;
    assert!((total - l.accepted_power_w).abs() < 1e-12 * l.accepted_power_w);
    assert!(l.radiation_efficiency < 1.0 && l.radiation_efficiency > 0.0);
    // Crossing wires are rejected.
    let cross = OrientedThinWire {
        start_m: WirePoint3::new(-0.1, 0.0, 0.0),
        end_m: WirePoint3::new(0.1, 0.0, 0.0),
        ..wire
    };
    cfg.wires = vec![wire, cross];
    cfg.loads.clear();
    assert!(solve_oriented_thin_wires(&cfg).is_err());
}

#[test]
fn sinusoidal_and_hertzian_dipoles() {
    let half = SinusoidalDipoleSolver::new(0.5, 299792458.0, 1.0).unwrap();
    assert!((half.radiation_resistance(8192).unwrap() - 73.08).abs() < 0.05);
    assert!((half.directivity(8192).unwrap() - 1.6409).abs() < 1e-3);
    let p = half.normalized_pattern(181).unwrap();
    assert!((p[90].normalized_power - 1.0).abs() < 1e-12 && p[0].normalized_power == 0.0);
    let (l, i, f, r) = (0.01, 2.0, 1e9, 10.0);
    let h = hertzian_dipole_far_field(l, i, f, PI / 2.0, r).unwrap();
    let k = 2.0 * PI * f / 299792458.0;
    assert!((h.e_theta.norm() - 376.730313668 * k * i * l / (4.0 * PI * r)).abs() < 1e-12);
    assert!((shielding_effectiveness_db(1.0, 0.01).unwrap() - 40.0).abs() < 1e-12);
}

#[test]
fn polarization_matching_array_and_optimizer() {
    let rhc = polarization_metrics(C::new(1.0, 0.0), C::new(0.0, 1.0)).unwrap();
    assert!(
        (rhc.axial_ratio - 1.0).abs() < 1e-12
            && (rhc.ellipticity_angle_rad.abs() - PI / 4.0).abs() < 1e-12
    );
    let m = antenna_match_metrics(C::new(50.0, 0.0), 50.0).unwrap();
    assert!(m.return_loss_db.is_infinite() && m.vswr == 1.0);
    let m3 = antenna_match_metrics(C::new(150.0, 0.0), 50.0).unwrap();
    assert!((m3.vswr - 3.0).abs() < 1e-12);
    let el: Vec<ArrayElement> = (0..4)
        .map(|n| ArrayElement {
            x_m: 0.5 * n as f64,
            ..Default::default()
        })
        .collect();
    assert!((array_factor(&el, 299792458.0, 0.0, 0.0).unwrap().norm() - 1.0).abs() < 1e-12);
    assert!(
        array_factor(&el, 299792458.0, PI / 2.0, 0.0)
            .unwrap()
            .norm()
            < 1e-12
    );
    let vars = [AntennaOptimizationVariable {
        initial_value: 0.0,
        minimum: -5.0,
        maximum: 5.0,
        initial_step: 1.0,
    }];
    let r = optimize_antenna(
        &vars,
        |x| (x[0] - 1.234).powi(2),
        &AntennaOptimizationConfig::default(),
    )
    .unwrap();
    assert!(r.converged && (r.parameters[0] - 1.234).abs() < 1e-5);
    assert!(optimize_antenna(&vars, |_| f64::NAN, &AntennaOptimizationConfig::default()).is_err());
}

fn square(n: usize, h: f64) -> Vec<SurfaceTriangle> {
    let mut t = vec![];
    for i in 0..n {
        for j in 0..n {
            let (x, y) = (i as f64 * h, j as f64 * h);
            t.push(SurfaceTriangle {
                a: v(x, y, 0.0),
                b: v(x + h, y, 0.0),
                c: v(x + h, y + h, 0.0),
            });
            t.push(SurfaceTriangle {
                a: v(x, y, 0.0),
                b: v(x + h, y + h, 0.0),
                c: v(x, y + h, 0.0),
            });
        }
    }
    t
}

#[test]
fn physical_optics_broadside_plate_rcs() {
    let side = 0.2;
    let tri = square(4, side / 4.0);
    let f = 10e9;
    let lambda = 299792458.0 / f;
    let r = solve_physical_optics(
        &tri,
        &PhysicalOpticsConfig {
            frequency_hz: f,
            observation_distance_m: 1e3,
            ..Default::default()
        },
    )
    .unwrap();
    let expected = 4.0 * PI * (side * side).powi(2) / (lambda * lambda);
    assert!(
        (r.rcs_m2 / expected - 1.0).abs() < 1e-9,
        "{} vs {}",
        r.rcs_m2,
        expected
    );
    assert!((r.illuminated_area_m2 - side * side).abs() < 1e-15);
}

#[test]
fn dihedral_retroreflects_and_surface_mom_validates() {
    let t = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| SurfaceTriangle {
        a: v(a[0], a[1], a[2]),
        b: v(b[0], b[1], b[2]),
        c: v(c[0], c[1], c[2]),
    };
    let corner = [
        t([0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]),
        t([0.0; 3], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]),
        t([0.0; 3], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0]),
        t([0.0; 3], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]),
    ];
    let d = v(-1.0, 0.1, -1.2);
    let s = trace_shooting_bouncing_rays(
        &corner,
        Ray {
            origin: v(0.6, 0.4, 0.8),
            direction: d,
            ..Default::default()
        },
        &SbrConfig {
            frequency_hz: 1e9,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(s.escaped && s.hits.len() == 2);
    let n = d.norm();
    let e = s.escaped_ray.direction;
    assert!(
        (e.x - 1.0 / n).abs() < 1e-12
            && (e.y - 0.1 / n).abs() < 1e-12
            && (e.z - 1.2 / n).abs() < 1e-12
    );
    // Non-manifold edge and zero-contrast SIE.
    let mut bad = square(1, 0.1);
    bad.push(SurfaceTriangle {
        a: v(0.0, 0.0, 0.0),
        b: v(0.1, 0.1, 0.0),
        c: v(0.05, 0.05, 0.1),
    });
    assert!(solve_pec_surface_mom(
        &bad,
        &SurfaceMomConfig {
            frequency_hz: 1e9,
            ..Default::default()
        }
    )
    .is_err());
    let z = solve_dielectric_surface_sie(
        &square(1, 0.1),
        &DielectricSieConfig {
            frequency_hz: 1e9,
            relative_permittivity_inside: 1.0,
            ..Default::default()
        },
        v(1.0, 0.0, 0.0),
    )
    .unwrap();
    assert!(z.converged && z.normal_polarization.iter().all(|c| c.norm() == 0.0));
    let ok = solve_pec_surface_mom(
        &square(3, 0.05),
        &SurfaceMomConfig {
            frequency_hz: 1e9,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(ok.relative_residual < 1e-10);
}

#[test]
fn treecode_approaches_direct_sum() {
    let src: Vec<Vec3> = (0..200)
        .map(|i| {
            v(
                (0.37 * i as f64).sin(),
                (0.91 * i as f64).cos(),
                (0.13 * i as f64).sin() * 0.5,
            )
        })
        .collect();
    let q: Vec<C> = (0..200).map(|i| C::new(1.0, 0.01 * i as f64)).collect();
    let tgt = vec![v(20.0, 0.0, 0.0), v(0.0, -15.0, 3.0)];
    let direct = helmholtz_direct_apply(&tgt, &src, &q, 0.05).unwrap();
    let exact = helmholtz_treecode_apply(
        &tgt,
        &src,
        &q,
        0.05,
        &HelmholtzTreecodeConfig {
            opening_angle: 1e-9,
            ..Default::default()
        },
    )
    .unwrap();
    let coarse = helmholtz_treecode_apply(
        &tgt,
        &src,
        &q,
        0.05,
        &HelmholtzTreecodeConfig {
            opening_angle: 0.2,
            ..Default::default()
        },
    )
    .unwrap();
    for ((a, b), c) in direct.iter().zip(&exact).zip(&coarse) {
        assert!((a - b).norm() < 1e-12 * a.norm());
        assert!((a - c).norm() < 0.1 * a.norm());
    }
}
