//! Value-by-value comparison with the C++ upstream `cfd/em` units. Goldens come from
//! `parity/oracle_em.cpp`.
use cfd_solvers::{em::*, fem, rf::ComplexMatrix, Complex64 as C};
use serde_json::Value;
use std::f64::consts::PI;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-em.json")).unwrap()
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
            assert!(r == c, "{name}[{i}]: rust {r} vs C++ {c}");
            continue;
        }
        assert!(
            (r - c).abs() <= rel * c.abs().max(floor) + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn cx(c: C) -> [f64; 2] {
    [c.re, c.im]
}
fn flatc(v: &[C]) -> Vec<f64> {
    v.iter().flat_map(|c| [c.re, c.im]).collect()
}
fn b(v: bool) -> f64 {
    f64::from(u8::from(v))
}

#[test]
fn ports_cables_and_sheets_match_upstream() {
    let modes = rectangular_waveguide_modes(22.86e-3, 10.16e-3, 12e9, 1.0, 1.0, 6).unwrap();
    let mut wm = Vec::new();
    for m in &modes {
        wm.extend([
            if m.family == ModeFamily::Te { 0.0 } else { 1.0 },
            m.m as f64,
            m.n as f64,
            m.cutoff_frequency_hz,
        ]);
        wm.extend(cx(m.propagation_constant_rad_per_m));
        wm.extend(cx(m.wave_impedance_ohm));
        wm.extend([b(m.propagating), m.electric_scale]);
    }
    close("wg_modes", &wm, 1e-12, 0.0);
    let (tf, tf2) = (
        rectangular_waveguide_transverse_fields(&modes[0], 22.86e-3, 10.16e-3, 7e-3, 3e-3).unwrap(),
        rectangular_waveguide_transverse_fields(&modes[4], 22.86e-3, 10.16e-3, 7e-3, 3e-3).unwrap(),
    );
    close(
        "wg_fields",
        &[
            flatc(&tf.electric_v_per_m),
            flatc(&tf.magnetic_a_per_m),
            flatc(&tf2.electric_v_per_m),
            flatc(&tf2.magnetic_a_per_m),
        ]
        .concat(),
        1e-12,
        1e-15,
    );
    close(
        "wg_power",
        &[
            rectangular_mode_power_w(&modes[0], 22.86e-3, 10.16e-3, 16, 16).unwrap(),
            rectangular_mode_power_w(&modes[2], 22.86e-3, 10.16e-3, 8, 8).unwrap(),
        ],
        1e-12,
        1e-15,
    );

    let rl = MulticonductorRlcg {
        conductors: 2,
        resistance: vec![0.1, 0.01, 0.01, 0.12],
        inductance: vec![2.5e-7, 0.5e-7, 0.5e-7, 2.6e-7],
        conductance: vec![1e-5, 0.0, 0.0, 1e-5],
        capacitance: vec![1.0e-10, -0.2e-10, -0.2e-10, 1.1e-10],
    };
    let load = ComplexMatrix::from_row_major(vec![
        C::new(50.0, 5.0),
        C::new(2.0, 0.0),
        C::new(2.0, 0.0),
        C::new(75.0, 0.0),
    ])
    .unwrap();
    let vin = [C::new(1.0, 0.0), C::new(0.0, 0.2)];
    let cr = solve_multiconductor_cable(3.0, 5, 2e7, &rl, &vin, &load).unwrap();
    close(
        "cable",
        &[
            flatc(&cr.input_current_a),
            flatc(&cr.load_voltage_v),
            flatc(&cr.load_current_a),
            flatc(cr.transfer_matrix.data()),
        ]
        .concat(),
        1e-11,
        1e-15,
    );
    let sw = sweep_multiconductor_cable(
        1.0,
        3,
        &[1e5, 1e7],
        &|f| {
            let mut r = rl.clone();
            r.resistance[0] *= 1.0 + f * 1e-7;
            r
        },
        &vin,
        &load,
    )
    .unwrap();
    close(
        "cable_sweep",
        &[flatc(&sw[0].load_voltage_v), flatc(&sw[1].load_voltage_v)].concat(),
        1e-11,
        1e-15,
    );
    let sh = CableShieldTransferModel {
        dc_transfer_resistance_ohm_per_m: 0.01,
        transfer_inductance_h_per_m: 1e-9,
        skin_corner_hz: 2e6,
    };
    let v = [
        cx(cable_shield_transfer_impedance_ohm_per_m(5e6, &sh).unwrap()),
        cx(cable_shield_induced_voltage_v(5e6, C::new(0.3, -0.1), 2.0, &sh).unwrap()),
        cx(
            field_to_cable_open_circuit_voltage_v(C::new(10.0, 2.0), 0.05, C::new(0.8, 0.1))
                .unwrap(),
        ),
    ]
    .concat();
    close("shield", &v, 1e-14, 0.0);
    let g1 = GstcSheetConfig {
        electric_surface_admittance_s: C::new(0.002, 0.004),
        magnetic_surface_impedance_ohm: C::new(30.0, -60.0),
        ..Default::default()
    };
    let g2 = GstcSheetConfig {
        electric_surface_admittance_s: C::new(0.0, -0.01),
        wave_impedance_after_ohm: 250.0,
        ..Default::default()
    };
    let (s1, s2) = (
        solve_gstc_normal_incidence(&g1).unwrap(),
        solve_gstc_normal_incidence(&g2).unwrap(),
    );
    let v = [
        cx(s1.reflection).to_vec(),
        cx(s1.transmission).to_vec(),
        vec![s1.reflected_power, s1.transmitted_power, s1.absorbed_power],
        cx(s2.reflection).to_vec(),
        cx(s2.transmission).to_vec(),
        vec![s2.absorbed_power],
    ]
    .concat();
    close("gstc", &v, 1e-13, 1e-15);
}

#[test]
fn frequency_domain_and_system_match_upstream() {
    let fc = FrequencyDomain1DConfig {
        points: 21,
        length_m: 2.0,
        frequency_hz: 3e8,
        relative_permittivity: 2.0,
        conductivity_s_per_m: 1e-3,
        ..Default::default()
    };
    let j: Vec<C> = (0..21)
        .map(|i| {
            C::new(
                (-0.1 * f64::from((i - 7) * (i - 7))).exp(),
                0.05 * f64::from(i),
            )
        })
        .collect();
    let fd = solve_pec_driven_maxwell_1d(&fc, &j).unwrap();
    close(
        "fd1d",
        &[
            flatc(&fd.electric_v_per_m),
            flatc(&fd.magnetic_a_per_m),
            vec![fd.spacing_m],
        ]
        .concat(),
        1e-10,
        1e-14,
    );
    let eps: Vec<f64> = (0..30).map(|i| if i < 15 { 1.0 } else { 4.0 }).collect();
    let mut cav = Vec::new();
    for m in pec_cavity_eigenmodes_1d(30, 1.0, &eps, 1.0, 3).unwrap() {
        cav.extend([m.angular_frequency_rad_per_s, m.frequency_hz]);
        cav.extend(m.electric_shape);
    }
    for m in pec_cavity_eigenmodes_1d_uniform(12, 0.5, 2.5, 1.2, 2).unwrap() {
        cav.push(m.frequency_hz);
        cav.extend(m.electric_shape);
    }
    close("cavity", &cav, 1e-9, 1e-12);

    let bits = [0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0];
    let mut wave = Vec::new();
    let mut prev = 0.0;
    for (k, bit) in bits.iter().enumerate() {
        for s in 0..8 {
            let target = if *bit == 1 { 1.0 } else { 0.0 };
            prev += 0.45 * (target - prev);
            wave.push(prev + 0.02 * (0.7 * (k * 8 + s) as f64).sin());
        }
    }
    let (e, e2) = (
        analyze_nrz_eye(&wave, &bits, 8, 0).unwrap(),
        analyze_nrz_eye(&wave, &bits, 8, 6).unwrap(),
    );
    close(
        "eye",
        &[
            e.zero_mean,
            e.one_mean,
            e.zero_sigma,
            e.one_sigma,
            e.decision_threshold,
            e.eye_height,
            e.eye_width_ui,
            e.estimated_ber,
            e2.one_mean,
            e2.estimated_ber,
        ],
        1e-9,
        0.0,
    );
    let cap = |c, esr, esl, count| DecouplingCapacitor {
        capacitance_f: c,
        esr_ohm: esr,
        esl_h: esl,
        count,
    };
    let lib = [
        cap(1e-6, 0.01, 1e-9, 1),
        cap(1e-7, 0.02, 0.5e-9, 1),
        cap(1e-8, 0.03, 0.3e-9, 2),
    ];
    let zero = C::new(0.0, 0.0);
    let v = [
        cx(decoupling_impedance_ohm(&lib[0], 1e6).unwrap()),
        cx(decoupling_impedance_ohm(&lib[2], 0.0).unwrap()),
        cx(pdn_parallel_impedance_ohm(&lib, 3e7, C::new(0.001, 0.02)).unwrap()),
        cx(pdn_parallel_impedance_ohm(&[], 1e6, zero).unwrap()),
    ]
    .concat();
    close("decap", &v, 1e-13, 0.0);
    let opt = optimize_decoupling_greedy(
        &[1e5, 1e6, 1e7, 1e8],
        &[0.05, 0.02, 0.02, 0.05],
        &lib,
        6,
        zero,
    )
    .unwrap();
    let mut ov = vec![
        opt.worst_ratio_before,
        opt.worst_ratio_after,
        opt.iterations as f64,
    ];
    ov.extend(opt.selected.iter().map(|c| c.capacitance_f));
    close("decap_opt", &ov, 1e-12, 0.0);
    let pg = PdnGridConfig {
        nx: 4,
        ny: 3,
        horizontal_resistance_ohm: 2e-3,
        vertical_resistance_ohm: 3e-3,
        relative_tolerance: 1e-10,
        ..Default::default()
    };
    let mut loads = vec![0.0; 12];
    (loads[5], loads[10], loads[3]) = (1.5, 0.7, -0.2);
    let ir = solve_pdn_ir_drop_grid(
        &pg,
        &loads,
        &[
            PdnVoltageSource {
                node: 0,
                voltage_v: 1.0,
            },
            PdnVoltageSource {
                node: 11,
                voltage_v: 0.99,
            },
        ],
    )
    .unwrap();
    close(
        "pdn",
        &[
            ir.voltage_v.clone(),
            vec![
                ir.minimum_voltage_v,
                ir.maximum_drop_v,
                ir.iterations as f64,
                b(ir.converged),
            ],
        ]
        .concat(),
        1e-13,
        0.0,
    );
    let dexp = sample_double_exponential_pulse(50e-9, 2e-9, 2.0, 1e-9, 10e-9).unwrap();
    let ds = sample_damped_sine(10e-9, 0.7e-9, 1.5, 2e8, 4e-9).unwrap();
    let pr = probe_waveform(&ds, 0.7e-9).unwrap();
    close(
        "pulses",
        &[
            dexp,
            ds,
            vec![
                pr.peak,
                pr.rms,
                pr.impulse,
                pr.energy,
                normalized_double_exponential_pulse(-1.0, 1.0, 1.0, 2.0).unwrap(),
            ],
        ]
        .concat(),
        1e-13,
        1e-16,
    );
}

#[test]
fn edge_elements_match_upstream() {
    let mesh2 = fem::make_rectangle_tri_mesh(3, 2, 1.0, 0.6).unwrap();
    let c2 = EdgeMaxwell2DConfig {
        frequency_hz: 2e8,
        relative_permittivity: 2.0,
        conductivity_s_per_m: 0.01,
        ..Default::default()
    };
    let r2 = solve_driven_edge_maxwell_2d(&mesh2, &c2, &|p: fem::Node2| {
        [C::new(p.y, 0.1), C::new(-p.x, 0.0)]
    })
    .unwrap();
    let mut e2: Vec<f64> = r2
        .edges
        .iter()
        .flat_map(|e| [e[0] as f64, e[1] as f64])
        .collect();
    e2.extend(flatc(&r2.edge_voltage_v));
    for f in &r2.electric_centroid_v_per_m {
        e2.extend(flatc(f));
    }
    close("edge2d", &e2, 1e-10, 1e-14);

    let mesh3 = fem::make_box_tet_mesh(2, 2, 2, 1.0, 0.8, 0.6).unwrap();
    let mut c3 = EdgeMaxwell3DConfig {
        frequency_hz: 1.5e8,
        relative_permittivity: 2.0,
        relative_permittivity_z: 3.0,
        conductivity_s_per_m: 0.05,
        ..Default::default()
    };
    c3.linear.relative_tolerance = 1e-12;
    let ports = [EdgeCurrentPort3D {
        node_a: 13,
        node_b: 4,
        current_a: C::new(1.0, 0.0),
    }];
    let src = |p: fem::Point3| [C::new(p.z, 0.0), C::new(0.0, p.x), C::new(0.2, 0.0)];
    let r3 = solve_driven_edge_maxwell_3d(&mesh3, &c3, Some(&src), &ports).unwrap();
    let mut e3 = flatc(&r3.edge_voltage_v);
    for f in r3
        .electric_centroid_v_per_m
        .iter()
        .chain(&r3.magnetic_centroid_a_per_m)
    {
        e3.extend(flatc(f));
    }
    // GMRES to 1e-12 on both sides: agreement is limited by the iterative solve.
    close("edge3d", &e3, 1e-7, 1e-9);
    let port_edge = r3.edges.iter().position(|e| *e == [4, 13]).unwrap();
    let th = edge_port_thevenin(&r3, port_edge, C::new(1.0, 0.0)).unwrap();
    let q = resonator_quality_3d(&mesh3, &c3, &r3, 5.8e7, 1.0).unwrap();
    close(
        "edge3d_post",
        &[
            cx(th.open_circuit_voltage_v).to_vec(),
            cx(th.port_impedance_ohm).to_vec(),
            vec![
                port_edge as f64,
                q.electric_energy_j,
                q.magnetic_energy_j,
                q.dielectric_loss_w,
                q.conductor_loss_w,
                q.quality_factor,
            ],
        ]
        .concat(),
        1e-7,
        0.0,
    );
    let ind = maxwell_face_jump_indicators_3d(&mesh3, &r3, 376.730313668).unwrap();
    let mk = mark_maxwell_dorfler(&ind, 0.4).unwrap();
    let refined = fem_refine(&mesh3, &mk);
    let mut v = ind.clone();
    v.extend(mk.iter().map(|m| b(*m)));
    v.extend(refined);
    close("adapt", &v, 1e-6, 1e-9);
}

fn fem_refine(mesh: &fem::Mesh3D, marks: &[bool]) -> Vec<f64> {
    let r = refine_tet4_marked_longest_edges(mesh, marks).unwrap();
    let mut out: Vec<f64> = r.nodes.iter().flat_map(|p| [p.x, p.y, p.z]).collect();
    out.extend(r.tetrahedra.iter().flat_map(|t| t.node.map(|n| n as f64)));
    out.extend(r.boundary_node.iter().map(|x| b(*x)));
    out
}

#[test]
fn photonic_modes_match_upstream() {
    let eps: Vec<f64> = (0..16)
        .map(|i| {
            if i % 4 == 1 || i / 4 == 2 {
                9.0
            } else {
                1.0 + 0.1 * f64::from(i)
            }
        })
        .collect();
    let grid = PhotonicGrid2D {
        nx: 4,
        ny: 4,
        width_m: 1e-6,
        height_m: 1e-6,
        relative_permittivity: eps,
        dispersive_permittivity: None,
    };
    let cfg = PhotonicModeConfig {
        mode_count: 3,
        ..Default::default()
    };
    let mut pb = Vec::new();
    for m in solve_photonic_band_modes(
        &grid,
        BlochWaveVector2D {
            kx_rad_per_m: PI / 2e-6,
            ky_rad_per_m: 0.3e6,
        },
        &cfg,
    )
    .unwrap()
    {
        pb.extend([
            m.eigenvalue_per_m2,
            m.frequency_hz,
            m.relative_residual,
            m.iterations as f64,
            b(m.converged),
        ]);
        pb.extend(flatc(&m.field));
    }
    let disp = PhotonicGrid2D {
        nx: 3,
        ny: 3,
        width_m: 1e-6,
        height_m: 1e-6,
        relative_permittivity: vec![],
        dispersive_permittivity: Some(Box::new(|i, f| 2.0 + 0.1 * i as f64 + 1e-15 * f)),
    };
    let cfg2 = PhotonicModeConfig {
        mode_count: 2,
        relative_tolerance: 1e-9,
        ..Default::default()
    };
    for m in solve_photonic_band_modes(
        &disp,
        BlochWaveVector2D {
            kx_rad_per_m: 1e6,
            ky_rad_per_m: 0.0,
        },
        &cfg2,
    )
    .unwrap()
    {
        pb.extend([m.frequency_hz, m.iterations as f64, m.relative_residual]);
    }
    // Residuals are round-off level (~1e-16) and only agree in magnitude.
    close("photonic_band", &pb, 1e-8, 1e-12);
    let eps: Vec<f64> = (0..20)
        .map(|i| {
            if (1..=3).contains(&(i % 5)) && (1..=2).contains(&(i / 5)) {
                12.0
            } else {
                2.1
            }
        })
        .collect();
    let wg = PhotonicGrid2D {
        nx: 5,
        ny: 4,
        width_m: 2e-6,
        height_m: 1.5e-6,
        relative_permittivity: eps,
        dispersive_permittivity: None,
    };
    let mut pw = Vec::new();
    for m in solve_photonic_waveguide_modes(
        &wg,
        &PhotonicWaveguideConfig {
            frequency_hz: 1.9e14,
            mode_count: 3,
            ..Default::default()
        },
    )
    .unwrap()
    {
        pw.extend([
            m.transverse_eigenvalue_per_m2,
            m.propagation_constant_rad_per_m,
            m.effective_index,
            m.relative_residual,
            b(m.guided),
        ]);
        pw.extend(flatc(&m.field));
    }
    close("photonic_wg", &pw, 1e-8, 1e-12);
}
