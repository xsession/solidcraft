//! Value-by-value comparison with the C++ upstream rf/nport, rf/peec, rf/conductor_models and
//! rf/modal_elements units. Goldens come from `parity/oracle_rf.cpp`.
use cfd_solvers::{core::Vec3, rf::*, Complex64 as C};
use serde_json::Value;
use std::f64::consts::PI;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-rf.json")).unwrap()
}

fn num(v: &Value) -> f64 {
    match v.as_str() {
        Some("inf") => f64::INFINITY,
        Some("-inf") => f64::NEG_INFINITY,
        Some(_) => f64::NAN,
        None => v.as_f64().unwrap(),
    }
}

#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(num)
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    let floor = 1e-6
        * cpp
            .iter()
            .filter(|c| c.is_finite())
            .fold(0.0_f64, |m, c| m.max(c.abs()));
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        if !c.is_finite() {
            assert!(
                r == c || (r.is_nan() && c.is_nan()),
                "{name}[{i}]: rust {r} vs C++ {c}"
            );
            continue;
        }
        assert!(
            (r - c).abs() <= rel * c.abs().max(floor) + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn text(name: &str) -> String {
    goldens()[name].as_str().unwrap().to_string()
}

fn flat(m: &ComplexMatrix) -> Vec<f64> {
    m.data().iter().flat_map(|c| [c.re, c.im]).collect()
}
fn flatc(m: &[C]) -> Vec<f64> {
    m.iter().flat_map(|c| [c.re, c.im]).collect()
}
fn cx(c: C) -> [f64; 2] {
    [c.re, c.im]
}
fn make(n: usize, seed: f64) -> ComplexMatrix {
    let mut m = ComplexMatrix::new(n);
    for r in 0..n {
        for c in 0..n {
            let (rf, cf) = (r as f64, c as f64);
            m[(r, c)] = C::new(
                (seed + 3.0 * rf + cf).sin() * if r == c { 2.0 } else { 0.3 },
                (seed * 0.7 + rf - 2.0 * cf).cos() * 0.2,
            );
        }
    }
    m
}
fn b(v: bool) -> f64 {
    f64::from(u8::from(v))
}

#[test]
fn nport_conversions_and_touchstone_match_upstream() {
    let mut z = make(3, 0.4);
    for i in 0..3 {
        z[(i, i)] += C::new(60.0 + 10.0 * i as f64, 5.0);
    }
    let z0 = [50.0, 75.0, 30.0];
    let s = nport_z_to_s(&z, &z0).unwrap();
    close("z_to_s", &flat(&s), 1e-12, 1e-15);
    close(
        "s_to_z",
        &flat(&nport_s_to_z(&s, &z0).unwrap()),
        1e-11,
        1e-12,
    );
    close(
        "renorm",
        &flat(&renormalize_s(&s, &z0, &[50.0; 3]).unwrap()),
        1e-11,
        1e-14,
    );
    let zc = [C::new(50.0, 10.0), C::new(75.0, -20.0), C::new(30.0, 0.0)];
    let sp = z_to_s_power_wave(&z, &zc).unwrap();
    close("pw_s", &flat(&sp), 1e-12, 1e-15);
    close(
        "pw_z",
        &flat(&s_to_z_power_wave(&sp, &zc).unwrap()),
        1e-11,
        1e-12,
    );
    close(
        "pw_renorm",
        &flat(&renormalize_s_power_wave(&sp, &zc, &[C::new(50.0, 0.0); 3]).unwrap()),
        1e-11,
        1e-14,
    );
    close(
        "inverse",
        &flat(&inverse(&make(4, 1.3)).unwrap()),
        1e-12,
        1e-15,
    );

    let mut tv = Vec::new();
    let t1 = "! comment\n# GHZ S MA R 50\n1.0 0.5 30 0.8 -45\n0.8 -45 0.4 10 ! row\n2.0 0.4 60 0.7 -90 0.7 -90 0.3 20\n";
    let t2 = "[Version] 2.0\n[Number of Ports] 3\n[Reference] 50 75 100\n# MHZ S DB\n[Network Data]\n100 -3 10 -20 30 -25 40 -20 30 -6 -10 -30 0 -25 40 -30 0 -9 90\n[End]\n";
    for t in [t1, t2] {
        for p in read_touchstone(t, 0).unwrap() {
            tv.push(p.frequency_hz);
            tv.extend(flat(&p.s));
            tv.extend(p.reference_impedance);
        }
    }
    close("touchstone_read", &tv, 1e-14, 1e-16);
    let pts: Vec<NPortPoint> = (0..3)
        .map(|k| NPortPoint {
            frequency_hz: 1e9 * f64::from(k + 1) + 0.123,
            s: make(2, f64::from(k) * 0.5),
            reference_impedance: vec![50.0, 50.0],
        })
        .collect();
    assert_eq!(write_touchstone(&pts, true).unwrap(), text("touchstone_v2"));
    assert_eq!(
        write_touchstone(&pts, false).unwrap(),
        text("touchstone_v1")
    );
    let gam = [C::new(0.1, 20.0), C::new(0.05, 31.0), C::new(0.0, 12.0)];
    close(
        "shift",
        &flat(&shift_reference_planes(&s, &gam, &[0.01, -0.02, 0.005]).unwrap()),
        1e-12,
        1e-15,
    );
    close(
        "rl",
        &[
            return_loss_db(C::new(0.3, 0.4)).unwrap(),
            vswr(C::new(0.3, 0.4)).unwrap(),
            vswr(C::new(0.6, 0.8)).unwrap(),
        ],
        1e-14,
        0.0,
    );
}

#[test]
fn two_port_design_metrics_match_upstream() {
    let s2 = Matrix2C {
        a11: C::new(0.6, -0.3),
        a12: C::new(0.05, 0.02),
        a21: C::new(2.5, 1.1),
        a22: C::new(0.4, -0.2),
    };
    let sm = stability_metrics(&s2);
    let (sc, lc) = (source_stability_circle(&s2), load_stability_circle(&s2));
    let v = [
        vec![sm.rollett_k, sm.mu],
        cx(sm.determinant).to_vec(),
        cx(sc.center).to_vec(),
        vec![sc.radius, b(sc.valid)],
        cx(lc.center).to_vec(),
        vec![lc.radius, b(lc.valid)],
    ]
    .concat();
    close("stability", &v, 1e-12, 0.0);
    let zero = C::new(0.0, 0.0);
    close(
        "gain",
        &[
            transducer_gain(&s2, zero, zero).unwrap(),
            transducer_gain(&s2, C::new(0.2, -0.1), C::new(0.3, 0.25)).unwrap(),
        ],
        1e-12,
        0.0,
    );
    let np = NoiseParameters {
        minimum_noise_factor: 1.4,
        optimum_source_reflection: C::new(0.3, 0.2),
        equivalent_noise_resistance_ohm: 12.0,
        reference_impedance_ohm: 50.0,
    };
    let nc = noise_circle(&np, 1.8).unwrap();
    close(
        "noise",
        &[
            vec![
                noise_factor(&np, C::new(0.1, -0.2)).unwrap(),
                nc.radius,
                nc.noise_factor,
                b(nc.valid),
            ],
            cx(nc.center).to_vec(),
        ]
        .concat(),
        1e-12,
        0.0,
    );
    let lp: Vec<f64> = sample_load_pull(&s2, C::new(0.1, 0.05), 3, 5, 0.8)
        .unwrap()
        .iter()
        .flat_map(|q| {
            [
                q.load_reflection.re,
                q.load_reflection.im,
                q.transducer_gain,
            ]
        })
        .collect();
    close("load_pull", &lp, 1e-12, 1e-15);
    let mut big = make(3, 2.0);
    for i in 0..3 {
        big[(i, 0)] *= 1.5;
    }
    let nq = network_quality(&big, 1e-10, 1e-10).unwrap();
    close(
        "quality",
        &[
            nq.maximum_singular_value,
            nq.reciprocity_error,
            b(nq.passive),
            b(nq.reciprocal),
        ],
        1e-12,
        0.0,
    );
    close(
        "passivity",
        &flat(&enforce_passivity(&big, 0.9).unwrap()),
        1e-12,
        0.0,
    );

    let sweep: Vec<NPortPoint> = (0..9)
        .map(|k| {
            let f = 1e8 * f64::from(k);
            let e = C::new(0.0, -2.0 * PI * f * 1.3e-9).exp();
            let mut m = ComplexMatrix::new(2);
            m[(0, 0)] = e * 0.2;
            m[(1, 0)] = e * 0.9;
            m[(0, 1)] = e * 0.9 * C::new(0.0, -0.2 * f64::from(k)).exp();
            m[(1, 1)] = C::new(0.1, 0.01 * f64::from(k));
            NPortPoint {
                frequency_hz: f,
                s: m,
                reference_impedance: vec![50.0, 50.0],
            }
        })
        .collect();
    let bc = check_broadband_causality(&sweep, 1e-6, 1e-8).unwrap();
    close(
        "causality",
        &[
            bc.negative_time_energy_ratio,
            b(bc.causal),
            bc.impulse_samples as f64,
        ],
        1e-10,
        0.0,
    );
    let rs: Vec<NPortPoint> = (0..12)
        .map(|k| {
            let f = 1e7 * 1.6f64.powi(k);
            let sv = C::new(0.0, 2.0 * PI * f);
            let mut m = ComplexMatrix::new(2);
            m[(0, 0)] = C::new(0.3, 0.0) + C::new(0.5 * 2e8, 0.0) / (sv + 2e8);
            m[(0, 1)] = C::new(0.1, 0.0) / (sv / 1e9 + 1.0);
            m[(1, 0)] = m[(0, 1)];
            m[(1, 1)] = C::new(0.2, 0.0) - sv * 0.1 / (sv + 5e8);
            NPortPoint {
                frequency_hz: f,
                s: m,
                reference_impedance: vec![50.0, 50.0],
            }
        })
        .collect();
    let model = fit_rational_network(&rs, 4).unwrap();
    let v = [
        model.poles_rad_per_s.clone(),
        flat(&model.direct),
        flat(&model.residues[2]),
        flat(&model.evaluate(3.3e8).unwrap()),
    ]
    .concat();
    close("rational", &v, 1e-8, 1e-12);
    let cfg = AdaptiveNetworkSweepConfig {
        initial_points: 4,
        maximum_points: 30,
        relative_tolerance: 1e-2,
        absolute_tolerance: 1e-4,
        ..Default::default()
    };
    let ad = adaptive_network_sweep(
        1e6,
        1e10,
        &[50.0],
        &|f| ComplexMatrix::from_row_major(vec![C::new(1.0, 0.0) / C::new(1.0, f / 3e8)]).unwrap(),
        &cfg,
    )
    .unwrap();
    let mut av = vec![
        ad.evaluations as f64,
        ad.refinements as f64,
        b(ad.converged),
        ad.points.len() as f64,
    ];
    for p in &ad.points {
        av.push(p.frequency_hz);
        av.extend(cx(p.s[(0, 0)]));
    }
    close("adaptive", &av, 1e-12, 1e-15);
    close(
        "mixed",
        &flat(&single_ended_to_mixed_mode(&make(4, 0.9)).unwrap()),
        1e-13,
        1e-15,
    );
    let (l1, l2) = (
        synthesize_l_match(50.0, 200.0).unwrap(),
        synthesize_l_match(300.0, 50.0).unwrap(),
    );
    close(
        "lmatch",
        &[
            l1.series_reactance_ohm,
            l1.shunt_susceptance_siemens,
            b(l1.shunt_first),
            l2.series_reactance_ohm,
            l2.shunt_susceptance_siemens,
            b(l2.shunt_first),
        ],
        1e-14,
        0.0,
    );
}

#[test]
fn peec_and_conductor_models_match_upstream() {
    let v = |x, y, z| Vec3::new(x, y, z);
    let mut sys = PeecFilamentSystem::default();
    sys.add_segment(FilamentSegment {
        start: v(0.0, 0.0, 0.0),
        end: v(0.01, 0.0, 0.0),
        radius_m: 1e-4,
        conductivity_s_per_m: 5.8e7,
    })
    .unwrap();
    sys.add_segment(FilamentSegment {
        start: v(0.01, 0.0, 0.0),
        end: v(0.01, 0.008, 0.001),
        radius_m: 1.2e-4,
        conductivity_s_per_m: 5.8e7,
    })
    .unwrap();
    sys.add_segment(FilamentSegment {
        start: v(0.0, 0.002, 0.0),
        end: v(0.012, 0.003, 0.0),
        radius_m: 0.8e-4,
        conductivity_s_per_m: 3.5e7,
    })
    .unwrap();
    let (pm, pm4) = (sys.extract(6).unwrap(), sys.extract(4).unwrap());
    close(
        "peec_rl",
        &[
            pm.resistance_ohm,
            pm.partial_inductance_h,
            pm4.partial_inductance_h,
        ]
        .concat(),
        1e-13,
        1e-25,
    );
    let pc = sys.extract_capacitance(2.2, 4).unwrap();
    close(
        "peec_c",
        &[pc.coefficient_of_potential_v_per_c, pc.capacitance_f].concat(),
        1e-12,
        0.0,
    );
    close(
        "peec_z",
        &[
            flatc(&sys.impedance_matrix(1e8, 6).unwrap()),
            flatc(&sys.impedance_matrix_skin_effect(1e8, 2).unwrap()),
        ]
        .concat(),
        1e-13,
        1e-20,
    );
    close(
        "peec_i",
        &flatc(
            &sys.solve_currents(
                5e7,
                &[C::new(1.0, 0.0), C::new(0.0, 0.5), C::new(-0.2, 0.0)],
                6,
            )
            .unwrap(),
        ),
        1e-11,
        0.0,
    );
    assert_eq!(
        sys.to_spice_subcircuit("BUS", 6).unwrap(),
        text("peec_spice")
    );
    let seg = FilamentSegment {
        start: v(0.0, 0.0, 0.0),
        end: v(0.05, 0.0, 0.0),
        radius_m: 5e-4,
        conductivity_s_per_m: 5.8e7,
    };
    close(
        "wire",
        &[
            dc_resistance(&seg).unwrap(),
            self_partial_inductance(&seg).unwrap(),
            skin_depth_m(5.8e7, 1e6, 1.0).unwrap(),
            round_wire_ac_resistance(&seg, 1e8, 1.0).unwrap(),
            round_wire_ac_resistance(&seg, 0.0, 1.0).unwrap(),
            round_wire_proximity_factor(3e-3, 5e-4, 1e8, 5.8e7, 1.0).unwrap(),
        ],
        1e-13,
        0.0,
    );
    let gr = finite_ground_reflection(1e7, 0.01, 15.0, 0.6).unwrap();
    // Re(TE) cancels to ~0; libgcc complex division differs from the textbook one at 1e-16.
    close("ground", &[cx(gr.TE), cx(gr.TM)].concat(), 1e-13, 1e-15);
    let sub: Vec<f64> = subdivide_rectangular_peec(&seg, 2e-3, 1e-3, 3, 2)
        .unwrap()
        .iter()
        .flat_map(|f| [f.start.y, f.start.z, f.end.x, f.radius_m])
        .collect();
    close("subdivide", &sub, 1e-14, 1e-18);
    let md = modal_discontinuity_s(C::new(5.0, 12.0), C::new(0.002, -0.004), 50.0).unwrap();
    close(
        "modal",
        &[cx(md.a11), cx(md.a12), cx(md.a21), cx(md.a22)].concat(),
        1e-13,
        0.0,
    );
}
