//! Value-by-value comparison with the C++ upstream rf/{thin_wire_mom, antenna, surface_mom,
//! integral_equation_extensions, asymptotic}. Goldens come from `parity/oracle_rf_mom.cpp`.
use cfd_solvers::{core::Vec3, optics::Ray, rf::*, Complex64 as C};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-rf-mom.json")).unwrap()
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
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        if c.is_infinite() {
            assert_eq!(r, c, "{name}[{i}]");
            continue;
        }
        assert!(
            (r - c).abs() <= rel * c.abs() + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn cv(z: &[C]) -> Vec<f64> {
    z.iter().flat_map(|c| [c.re, c.im]).collect()
}

fn cmat(m: &ComplexMatrix) -> Vec<f64> {
    cv(m.data())
}

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

fn wp(x: f64, y: f64, z: f64) -> WirePoint3 {
    WirePoint3::new(x, y, z)
}

const R: f64 = 1e-9;

#[test]
fn thin_wire_dipole_and_networks() {
    let r = solve_center_fed_thin_wire(&ThinWireMomConfig {
        length_m: 0.5,
        radius_m: 1e-3,
        frequency_hz: 290e6,
        segments: 21,
        ..Default::default()
    })
    .unwrap();
    close(
        "dipole",
        &[
            r.feed_impedance_ohm.re,
            r.feed_impedance_ohm.im,
            r.feed_current_a,
        ],
        R,
        0.0,
    );
    close("dipole_I", &cv(&r.current_a), R, 1e-15);
    close("dipole_z", &r.segment_center_z_m, 1e-14, 1e-16);

    let w0 = OrientedThinWire {
        start_m: wp(0.0, 0.0, 0.05),
        end_m: wp(0.0, 0.0, 0.3),
        radius_m: 1e-3,
        segments: 9,
        conductivity_s_per_m: 5.8e7,
        ..Default::default()
    };
    let w1 = OrientedThinWire {
        start_m: wp(0.0, 0.0, 0.3),
        end_m: wp(0.2, 0.05, 0.3),
        radius_m: 1.5e-3,
        segments: 7,
        dielectric_loss_ohm_per_m: 0.3,
        ..Default::default()
    };
    let w2 = OrientedThinWire {
        start_m: wp(0.1, 0.2, 0.05),
        end_m: wp(0.1, 0.2, 0.35),
        radius_m: 1e-3,
        segments: 9,
        ..Default::default()
    };
    let oc = OrientedWireMomConfig {
        frequency_hz: 400e6,
        junction_tolerance_m: 1e-6,
        ground_model: ThinWireGroundModel::PecPlane,
        wires: vec![w0, w1, w2],
        feeds: vec![ThinWireDeltaGapFeed {
            wire: 0,
            segment: Some(1),
            voltage_v: C::new(1.0, 0.0),
        }],
        loads: vec![ThinWireSeriesLoad {
            wire: 2,
            segment: 4,
            impedance_ohm: C::new(25.0, 10.0),
        }],
        transmission_line_loads: vec![ThinWireTransmissionLineLoad {
            first_wire: 1,
            first_segment: 5,
            second_wire: 2,
            second_segment: 7,
            propagation_constant_per_m: C::new(0.1, 9.0),
            length_m: 0.12,
            second_current_sign: -1,
            ..Default::default()
        }],
        ..Default::default()
    };
    let o = solve_oriented_thin_wires(&oc).unwrap();
    close(
        "oriented",
        &[
            o.accepted_power_w,
            o.conductor_loss_w,
            o.dielectric_loss_w,
            o.lumped_load_loss_w,
            o.transmission_line_loss_w,
            o.radiated_power_w,
            o.radiation_efficiency,
            o.feed_current_a[0],
            o.feed_impedance_ohm[0].re,
            o.feed_impedance_ohm[0].im,
        ],
        1e-7,
        0.0,
    );
    close("oriented_I", &cv(&o.current_a), 1e-7, 1e-16);

    let a = OrientedThinWire {
        start_m: wp(0.0, 0.0, -0.24),
        end_m: wp(0.0, 0.0, 0.24),
        radius_m: 2e-3,
        segments: 11,
        ..Default::default()
    };
    let mut b = a;
    b.start_m.x = 0.25;
    b.end_m.x = 0.25;
    let mp = OrientedWireMomConfig {
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
    close(
        "zport",
        &cmat(&oriented_wire_port_impedance_matrix(&mp).unwrap()),
        R,
        1e-9,
    );
    close(
        "sport",
        &cmat(&oriented_wire_port_s_parameters(&mp, &[50.0, 75.0]).unwrap()),
        R,
        1e-12,
    );

    let pc = ParallelWireMomConfig {
        frequency_hz: 280e6,
        wires: vec![
            ParallelThinWire {
                length_m: 0.5,
                radius_m: 1e-3,
                segments: 11,
                ..Default::default()
            },
            ParallelThinWire {
                x_m: 0.15,
                center_z_m: 0.02,
                length_m: 0.45,
                radius_m: 1e-3,
                segments: 11,
                ..Default::default()
            },
        ],
        feeds: vec![ThinWireDeltaGapFeed::default()],
        loads: vec![ThinWireSeriesLoad {
            wire: 1,
            segment: 5,
            impedance_ohm: C::new(0.0, -30.0),
        }],
        ..Default::default()
    };
    close(
        "parallel_I",
        &cv(&solve_parallel_thin_wires(&pc).unwrap().current_a),
        R,
        1e-15,
    );
}

#[test]
fn antenna_utilities() {
    let mut an = vec![];
    let hf = hertzian_dipole_far_field(0.01, 2.0, 1e9, 0.7, 10.0).unwrap();
    an.extend(cv(&[hf.e_theta, hf.h_phi]));
    an.push(shielding_effectiveness_db(3.0, 0.02).unwrap());
    let el = [
        ArrayElement {
            weight: C::new(1.0, 0.0),
            ..Default::default()
        },
        ArrayElement {
            x_m: 0.1,
            weight: C::new(0.5, 0.5),
            ..Default::default()
        },
        ArrayElement {
            x_m: 0.2,
            y_m: 0.05,
            z_m: 0.01,
            weight: C::new(0.2, -0.7),
        },
    ];
    an.extend(cv(&[array_factor(&el, 1.5e9, 0.6, 0.3).unwrap()]));
    let pm = polarization_metrics(C::new(1.0, 0.2), C::new(0.3, -0.8)).unwrap();
    an.extend([
        pm.axial_ratio,
        pm.tilt_angle_rad,
        pm.ellipticity_angle_rad,
        if pm.right_hand { 1.0 } else { 0.0 },
    ]);
    an.push(
        polarization_metrics(C::new(1.0, 0.0), C::new(0.5, 0.0))
            .unwrap()
            .axial_ratio,
    );
    let mm = antenna_match_metrics(C::new(73.0, 42.5), 50.0).unwrap();
    an.extend(cv(&[mm.reflection_coefficient]));
    an.extend([mm.return_loss_db, mm.vswr]);
    let m0 = antenna_match_metrics(C::new(50.0, 0.0), 50.0).unwrap();
    an.extend([m0.return_loss_db, m0.vswr]);
    let vars = [
        AntennaOptimizationVariable {
            initial_value: 0.0,
            minimum: -1.0,
            maximum: 2.0,
            initial_step: 0.5,
        },
        AntennaOptimizationVariable {
            initial_value: 1.0,
            minimum: 0.0,
            maximum: 3.0,
            initial_step: 0.25,
        },
    ];
    let opt = optimize_antenna(
        &vars,
        |x| (x[0] - 0.7) * (x[0] - 0.7) + 2.0 * (x[1] - 2.2) * (x[1] - 2.2) + 0.1 * x[0] * x[1],
        &AntennaOptimizationConfig::default(),
    )
    .unwrap();
    an.extend([
        opt.parameters[0],
        opt.parameters[1],
        opt.objective,
        opt.iterations as f64,
        opt.evaluations as f64,
        if opt.converged { 1.0 } else { 0.0 },
    ]);
    let sd = SinusoidalDipoleSolver::new(0.75, 300e6, 1.5).unwrap();
    an.extend([
        sd.wavelength(),
        sd.electrical_length(),
        sd.pattern_factor(0.9),
        sd.radiated_power(8192).unwrap(),
        sd.radiation_resistance(8192).unwrap(),
        sd.directivity(8192).unwrap(),
        sd.radiated_power(100).unwrap(),
    ]);
    an.extend(cv(&[sd.far_e_theta(1.1, 20.0).unwrap()]));
    for s in sd.normalized_pattern(13).unwrap() {
        an.extend([s.theta_rad, s.normalized_power]);
    }
    close("antenna", &an, 1e-12, 1e-15);
}

fn plate() -> Vec<SurfaceTriangle> {
    let h = 0.06;
    let mut t = vec![];
    for i in 0..2 {
        for j in 0..2 {
            let (fi, fj) = (i as f64, j as f64);
            let p00 = v(fi * h, fj * h, 0.0);
            let p10 = v((fi + 1.0) * h, fj * h, 0.0);
            let p01 = v(fi * h, (fj + 1.0) * h, 0.0);
            let p11 = v((fi + 1.0) * h, (fj + 1.0) * h, 0.0);
            t.push(SurfaceTriangle {
                a: p00,
                b: p10,
                c: p11,
            });
            t.push(SurfaceTriangle {
                a: p00,
                b: p11,
                c: p01,
            });
        }
    }
    t
}

#[test]
fn surface_mom_and_sie() {
    let plate = plate();
    let sc = SurfaceMomConfig {
        frequency_hz: 1.2e9,
        plane_wave_direction: v(0.2, 0.0, -1.0),
        plane_wave_electric: v(1.0, 0.3, 0.0),
        ..Default::default()
    };
    let sr = solve_pec_surface_mom(&plate, &sc).unwrap();
    let mut sm = vec![
        sr.basis.len() as f64,
        sr.relative_residual,
        sr.wave_number_per_m,
        sr.wave_impedance_ohm,
    ];
    for b in &sr.basis {
        sm.extend([
            b.edge_length_m,
            b.plus_sign as f64,
            b.minus_sign as f64,
            b.minus_triangle.map_or(-1.0, |m| m as f64),
        ]);
    }
    close("surface", &sm, 1e-10, 1e-14);
    close("surface_I", &cv(&sr.current_a), 1e-9, 1e-15);
    let ff = surface_mom_far_field(&sr, v(0.1, 0.3, 1.0), 5.0).unwrap();
    close(
        "surface_ff",
        &[
            ff.e_x.re,
            ff.e_x.im,
            ff.e_y.re,
            ff.e_y.im,
            ff.e_z.re,
            ff.e_z.im,
            ff.power_density_w_m2,
        ],
        1e-9,
        1e-15,
    );
    let dg = SurfaceMomConfig {
        frequency_hz: 0.9e9,
        excitation: SurfaceMomExcitation::DeltaGap,
        feed_edge: Some(4),
        feed_voltage_v: C::new(2.0, 0.5),
        ..Default::default()
    };
    let dr = solve_pec_surface_mom(&plate, &dg).unwrap();
    close(
        "surface_gap",
        &[
            dr.feed_current_a.re,
            dr.feed_current_a.im,
            dr.feed_impedance_ohm.re,
            dr.feed_impedance_ohm.im,
            dr.accepted_power_w,
            dr.relative_residual,
        ],
        1e-9,
        1e-14,
    );

    let (a, b, c, d) = (
        v(0.0, 0.0, 0.0),
        v(0.05, 0.0, 0.0),
        v(0.0, 0.05, 0.0),
        v(0.0, 0.0, 0.05),
    );
    let tet = [
        SurfaceTriangle { a, b: c, c: b },
        SurfaceTriangle { a, b, c: d },
        SurfaceTriangle { a, b: d, c },
        SurfaceTriangle { a: b, b: c, c: d },
    ];
    let ds = solve_dielectric_surface_sie(
        &tet,
        &DielectricSieConfig {
            frequency_hz: 2e9,
            relative_permittivity_inside: 6.0,
            ..Default::default()
        },
        v(0.3, 1.0, -0.5),
    )
    .unwrap();
    let mut dsv = cv(&ds.normal_polarization);
    dsv.extend([ds.relative_residual, if ds.converged { 1.0 } else { 0.0 }]);
    close("sie", &dsv, 1e-12, 1e-16);

    let src: Vec<Vec3> = (0..40)
        .map(|i| {
            v(
                (1.3 * i as f64).sin(),
                (0.7 * i as f64).cos() * 0.5,
                0.1 * i as f64 - 2.0,
            )
        })
        .collect();
    let q: Vec<C> = (0..40)
        .map(|i| C::new((0.3 * i as f64).cos(), 0.1 * i as f64))
        .collect();
    let mut tgt: Vec<Vec3> = (0..5)
        .map(|i| v(5.0 + i as f64, -3.0 + 0.5 * i as f64, i as f64))
        .collect();
    tgt.push(src[3]);
    close(
        "helm_direct",
        &cv(&helmholtz_direct_apply(&tgt, &src, &q, 2.5).unwrap()),
        1e-12,
        1e-15,
    );
    let tc = HelmholtzTreecodeConfig {
        leaf_size: 4,
        opening_angle: 0.4,
        ..Default::default()
    };
    close(
        "helm_tree",
        &cv(&helmholtz_treecode_apply(&tgt, &src, &q, 2.5, &tc).unwrap()),
        1e-12,
        1e-15,
    );
}

#[test]
fn physical_optics_and_sbr() {
    let plate = plate();
    let mut poc = PhysicalOpticsConfig {
        frequency_hz: 10e9,
        incident_direction: v(0.1, 0.0, -1.0),
        observation_direction: v(-0.1, 0.0, 1.0),
        observation_distance_m: 100.0,
        ..Default::default()
    };
    let por = solve_physical_optics(&plate, &poc).unwrap();
    let mut pov = vec![
        por.illuminated_area_m2,
        por.rcs_m2,
        por.wave_number_per_m,
        por.wave_impedance_ohm,
    ];
    pov.extend(cv(&por.far_field));
    for cur in &por.currents {
        pov.extend(cv(&cur.current));
    }
    poc.incident_direction = v(0.0, 0.0, 1.0);
    pov.push(solve_physical_optics(&plate, &poc).unwrap().rcs_m2);
    close("po", &pov, 1e-12, 1e-17);

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
    let ray = Ray {
        origin: v(0.6, 0.4, 0.8),
        direction: v(-1.0, 0.1, -1.2),
        ..Default::default()
    };
    let sb = trace_shooting_bouncing_rays(
        &corner,
        ray,
        &SbrConfig {
            frequency_hz: 3e9,
            ..Default::default()
        },
    )
    .unwrap();
    let e = sb.escaped_ray;
    let mut sv = vec![
        sb.hits.len() as f64,
        sb.optical_path_m,
        if sb.escaped { 1.0 } else { 0.0 },
        e.origin.x,
        e.origin.y,
        e.origin.z,
        e.direction.x,
        e.direction.y,
        e.direction.z,
        sb.phase.re,
        sb.phase.im,
    ];
    for h in &sb.hits {
        sv.extend([
            h.triangle as f64,
            h.distance_m,
            h.point.x,
            h.point.y,
            h.point.z,
            h.normal.x,
            h.normal.y,
            h.normal.z,
        ]);
    }
    close("sbr", &sv, 1e-12, 1e-14);
}
