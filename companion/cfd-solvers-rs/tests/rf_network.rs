//! Analytic and robustness checks for N-port algebra, Touchstone I/O, design metrics, rational
//! fitting and PEEC extraction.
use cfd_solvers::{core::Vec3, rf::*, Complex64 as C};
use std::f64::consts::PI;

fn c(re: f64, im: f64) -> C {
    C::new(re, im)
}

#[test]
fn matched_and_round_trip_conversions() {
    let refs = [50.0, 75.0];
    let mut z = ComplexMatrix::new(2);
    z[(0, 0)] = c(50.0, 0.0);
    z[(1, 1)] = c(75.0, 0.0);
    let s = nport_z_to_s(&z, &refs).unwrap();
    assert!(s.data().iter().all(|v| v.norm() < 1e-14));
    let mut z2 = ComplexMatrix::new(2);
    z2[(0, 0)] = c(30.0, 12.0);
    z2[(0, 1)] = c(5.0, -3.0);
    z2[(1, 0)] = c(5.0, -3.0);
    z2[(1, 1)] = c(90.0, -40.0);
    let back = nport_s_to_z(&nport_z_to_s(&z2, &refs).unwrap(), &refs).unwrap();
    for (a, b) in back.data().iter().zip(z2.data()) {
        assert!((a - b).norm() < 1e-10);
    }
    // Real references: Kurokawa power waves reduce to the ordinary definition.
    let pw = z_to_s_power_wave(&z2, &[c(50.0, 0.0), c(75.0, 0.0)]).unwrap();
    for (a, b) in pw
        .data()
        .iter()
        .zip(nport_z_to_s(&z2, &refs).unwrap().data())
    {
        assert!((a - b).norm() < 1e-13);
    }
    assert!(nport_z_to_s(&z2, &[50.0]).is_err());
    assert!(inverse(&ComplexMatrix::new(2)).is_err());
    assert!(z_to_s_power_wave(&z2, &[c(-1.0, 0.0), c(50.0, 0.0)]).is_err());
}

#[test]
fn touchstone_round_trip_and_errors() {
    let mut m = ComplexMatrix::new(2);
    m[(0, 0)] = c(0.25, -0.5);
    m[(1, 0)] = c(0.75, 0.125);
    m[(0, 1)] = c(-0.5, 0.0625);
    m[(1, 1)] = c(0.0, 0.375);
    let pts = vec![NPortPoint {
        frequency_hz: 1.5e9,
        s: m.clone(),
        reference_impedance: vec![50.0, 50.0],
    }];
    for v2 in [true, false] {
        let text = write_touchstone(&pts, v2).unwrap();
        let back = read_touchstone(&text, 2).unwrap();
        assert_eq!(back.len(), 1);
        assert_eq!(back[0].frequency_hz, 1.5e9);
        assert_eq!(back[0].s, m);
    }
    assert!(read_touchstone("# HZ S RI R 50\n1 0 0 0\n", 2).is_err());
    assert!(read_touchstone("# THZ S RI R 50\n1 0 0 0 0 0 0 0 0\n", 2).is_err());
    assert!(read_touchstone("# HZ Y RI R 50\n", 2).is_err());
    assert!(read_touchstone("[Number of Ports 2\n", 2).is_err());
    assert!(read_touchstone("[Reference] 50 60 70\n# HZ S RI\n1 0 0 0 0 0 0 0 0\n", 2).is_err());
    let mixed = NPortPoint {
        reference_impedance: vec![50.0, 60.0],
        ..pts[0].clone()
    };
    assert!(write_touchstone(&[mixed], false).is_err());
}

#[test]
fn design_metrics_satisfy_textbook_identities() {
    // A unilateral matched amplifier: K infinite, gain |S21|^2.
    let s = Matrix2C {
        a11: c(0.0, 0.0),
        a12: c(0.0, 0.0),
        a21: c(3.0, 0.0),
        a22: c(0.0, 0.0),
    };
    assert_eq!(stability_metrics(&s).rollett_k, f64::INFINITY);
    assert!((transducer_gain(&s, c(0.0, 0.0), c(0.0, 0.0)).unwrap() - 9.0).abs() < 1e-14);
    assert!(transducer_gain(&s, c(1.0, 0.0), c(0.0, 0.0)).is_err());
    let np = NoiseParameters {
        minimum_noise_factor: 1.5,
        optimum_source_reflection: c(0.2, 0.1),
        equivalent_noise_resistance_ohm: 10.0,
        reference_impedance_ohm: 50.0,
    };
    assert!((noise_factor(&np, np.optimum_source_reflection).unwrap() - 1.5).abs() < 1e-15);
    // Every point on the noise circle has the target noise factor.
    let circle = noise_circle(&np, 2.0).unwrap();
    for k in 0..8 {
        let g = circle.center + C::from_polar(circle.radius, f64::from(k) * PI / 4.0);
        assert!((noise_factor(&np, g).unwrap() - 2.0).abs() < 1e-12);
    }
    assert!((vswr(c(0.0, 0.0)).unwrap() - 1.0).abs() < 1e-15);
    assert_eq!(return_loss_db(c(0.0, 0.0)).unwrap(), f64::INFINITY);
    // L-match transforms the load to the source resistance.
    for (rs, rl) in [(50.0, 200.0), (300.0, 50.0)] {
        let m = synthesize_l_match(rs, rl).unwrap();
        let (x, bb) = (m.series_reactance_ohm, m.shunt_susceptance_siemens);
        let zin = if m.shunt_first {
            // Series jX next to the load, shunt jB at the source side.
            let z = c(rl, x);
            C::new(1.0, 0.0) / (C::new(1.0, 0.0) / z + c(0.0, bb))
        } else {
            c(0.0, x) + C::new(1.0, 0.0) / (c(1.0 / rl, bb))
        };
        assert!((zin.re - rs).abs() < 1e-9 * rs, "{zin} vs {rs}");
    }
}

#[test]
fn passivity_causality_and_rational_fit() {
    let mut m = ComplexMatrix::new(2);
    m[(0, 1)] = c(1.5, 0.0);
    m[(1, 0)] = c(1.5, 0.0);
    let q = network_quality(&m, 1e-10, 1e-10).unwrap();
    assert!((q.maximum_singular_value - 1.5).abs() < 1e-10 && !q.passive && q.reciprocal);
    let p = enforce_passivity(&m, 1.0).unwrap();
    assert!(
        network_quality(&p, 1e-10, 1e-10)
            .unwrap()
            .maximum_singular_value
            <= 1.0 + 1e-10
    );
    // A pure delay of a few samples is causal.
    let n = 33;
    let df = 1e8;
    let tau = 3.0 / (2.0 * (n as f64 - 1.0) * df);
    let pts: Vec<NPortPoint> = (0..n)
        .map(|k| {
            let f = df * k as f64;
            NPortPoint {
                frequency_hz: f,
                s: ComplexMatrix::from_row_major(vec![C::new(0.0, -2.0 * PI * f * tau).exp()])
                    .unwrap(),
                reference_impedance: vec![50.0],
            }
        })
        .collect();
    let causal = check_broadband_causality(&pts, 1e-6, 1e-8).unwrap();
    assert!(causal.causal, "{}", causal.negative_time_energy_ratio);
    // A model built on the same pole grid is recovered exactly.
    let truth = |f: f64, poles: &[f64]| {
        let s = C::new(0.0, 2.0 * PI * f);
        C::new(0.1, 0.0)
            + poles
                .iter()
                .enumerate()
                .map(|(k, p)| C::new(-p * 0.2 * (k as f64 + 1.0), 0.0) / (s - p))
                .sum::<C>()
    };
    let freqs: Vec<f64> = (0..15).map(|k| 1e6 * 2f64.powi(k)).collect();
    let (fmin, fmax) = (freqs[0], freqs[14]);
    let ratio = (fmax / fmin).powf(1.0 / 2.0);
    let poles: Vec<f64> = (0..3).map(|k| -2.0 * PI * fmin * ratio.powi(k)).collect();
    let pts: Vec<NPortPoint> = freqs
        .iter()
        .map(|f| NPortPoint {
            frequency_hz: *f,
            s: ComplexMatrix::from_row_major(vec![truth(*f, &poles)]).unwrap(),
            reference_impedance: vec![50.0],
        })
        .collect();
    let model = fit_rational_network(&pts, 3).unwrap();
    assert!((model.direct[(0, 0)] - c(0.1, 0.0)).norm() < 1e-8);
    for f in [3e6, 7.7e7, 2e9] {
        assert!((model.evaluate(f).unwrap()[(0, 0)] - truth(f, &poles)).norm() < 1e-8);
    }
}

#[test]
fn peec_physics_identities() {
    let seg = |y: f64| FilamentSegment {
        start: Vec3::new(0.0, y, 0.0),
        end: Vec3::new(0.1, y, 0.0),
        radius_m: 1e-3,
        conductivity_s_per_m: 5.8e7,
    };
    let mut sys = PeecFilamentSystem::default();
    sys.add_segment(seg(0.0)).unwrap();
    sys.add_segment(seg(0.01)).unwrap();
    let m = sys.extract(6).unwrap();
    let l = &m.partial_inductance_h;
    assert!(l[1] > 0.0 && l[1] == l[2] && l[1] < l[0]);
    let cap = sys.extract_capacitance(1.0, 6).unwrap();
    assert!(cap.capacitance_f[0] > 0.0 && cap.capacitance_f[1] < 0.0);
    let v = [c(1.0, 0.0), c(0.0, 0.0)];
    let i = sys.solve_currents(1e6, &v, 6).unwrap();
    let z = sys.impedance_matrix(1e6, 6).unwrap();
    for r in 0..2 {
        let lhs = z[r * 2] * i[0] + z[r * 2 + 1] * i[1];
        assert!((lhs - v[r]).norm() < 1e-12);
    }
    assert!(sys
        .add_segment(FilamentSegment {
            radius_m: 0.05,
            ..seg(0.0)
        })
        .is_err());
    assert!(sys.to_spice_subcircuit("", 6).is_err());
    // Skin effect raises AC resistance above DC.
    assert!(
        round_wire_ac_resistance(&seg(0.0), 1e8, 1.0).unwrap() > dc_resistance(&seg(0.0)).unwrap()
    );
    // A near-perfect conductor reflects TE with -1 and TM with +1; a lossless dielectric at
    // Brewster's angle has no TM reflection.
    let g = fresnel_ground_reflection(1e6, 1e9, 1.0, 0.3).unwrap();
    assert!((g.TE + c(1.0, 0.0)).norm() < 1e-3 && (g.TM - c(1.0, 0.0)).norm() < 1e-3);
    let brewster = fresnel_ground_reflection(1e6, 0.0, 4.0, 2f64.atan()).unwrap();
    assert!(brewster.TM.norm() < 1e-12 && brewster.TE.norm() > 0.1);
    // The upstream-compatible variant keeps its documented angle-independent TM.
    let u = finite_ground_reflection(1e6, 0.0, 4.0, 0.7).unwrap();
    assert!((u.TM - c(0.6, 0.0)).norm() < 1e-12 && u.TE.norm() < 1e-12);
    let s = modal_discontinuity_s(c(0.0, 0.0), c(0.0, 0.0), 50.0).unwrap();
    assert!((s.a21 - c(1.0, 0.0)).norm() < 1e-15 && s.a11.norm() < 1e-15);
}
