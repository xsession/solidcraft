//! Value-by-value comparison with the C++ upstream circuit/{spice, analysis, raw, osdi_adapter,
//! external_spice_adapter, spice_compat}. Goldens come from `parity/oracle_circuit.cpp`.
#![allow(clippy::needless_range_loop, clippy::type_complexity)]
use cfd_solvers::circuit::*;
use serde_json::Value;
use std::{collections::HashMap, sync::Arc};

fn g() -> &'static Value {
    static G: std::sync::OnceLock<Value> = std::sync::OnceLock::new();
    G.get_or_init(|| serde_json::from_str(include_str!("../parity/goldens-circuit.json")).unwrap())
}

fn num(v: &Value) -> f64 {
    match v.as_str() {
        Some("inf") => f64::INFINITY,
        Some("-inf") => f64::NEG_INFINITY,
        Some(_) => f64::NAN,
        None => v.as_f64().unwrap(),
    }
}

/// Relative tolerance with an absolute floor scaled by the largest magnitude in the vector.
#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64) {
    let cpp: Vec<f64> = g()[name].as_array().unwrap().iter().map(num).collect();
    assert_eq!(rust.len(), cpp.len(), "{name} length");
    let scale = cpp
        .iter()
        .filter(|v| v.is_finite())
        .fold(0.0f64, |m, v| m.max(v.abs()));
    for (i, (a, b)) in rust.iter().zip(&cpp).enumerate() {
        let ok = a == b
            || (a.is_nan() && b.is_nan())
            || (a - b).abs() <= rel * a.abs().max(b.abs()) + rel * 1e-3 * scale;
        assert!(ok, "{name}[{i}]: rust {a} cpp {b}");
    }
}

#[track_caller]
fn s(name: &str, rust: &str) {
    assert_eq!(rust, g()[name].as_str().unwrap(), "{name}");
}

fn err<T>(r: cfd_solvers::Result<T>) -> String {
    match r {
        Ok(_) => String::new(),
        Err(e) => {
            let t = e.to_string();
            t.split_once(": ").map(|(_, m)| m.to_string()).unwrap_or(t)
        }
    }
}

fn op(g: &mut Vec<f64>, p: &OperatingPoint) {
    g.extend(&p.node_voltage);
    g.extend(&p.branch_current);
    g.push(p.iterations as f64);
    g.push(f64::from(u8::from(p.converged)));
}
fn cx(g: &mut Vec<f64>, c: Complex) {
    g.push(c.re);
    g.push(c.im);
}
fn acs(g: &mut Vec<f64>, a: &AcSolution) {
    g.push(a.frequency_hz);
    for &c in &a.node_voltage {
        cx(g, c);
    }
    for &c in &a.branch_current {
        cx(g, c);
    }
}
fn tr(g: &mut Vec<f64>, pts: &[TransientPoint]) {
    for p in pts {
        g.push(p.time);
        g.extend(&p.node_voltage);
        g.extend(&p.branch_current);
    }
}
fn cfg() -> NewtonConfig {
    NewtonConfig::default()
}
fn parse(t: &str) -> Circuit {
    Circuit::parse_spice(t).unwrap()
}
fn homotopy(c: &Circuit) -> OperatingPoint {
    c.dc_operating_point_default_homotopy(&cfg()).unwrap()
}

const NL: &str = "V1 vdd 0 5\nR1 vdd a 1k\nD1 a 0 DMOD\n.model DMOD D(IS=1e-14 N=1.5)\nVg g 0 2\nRd vdd dr 2k\nM1 dr g 0 0 NM\n.model NM NMOS(VTO=0.8 KP=2m LAMBDA=0.02)\nVb bb 0 0.7\nRb bb b 10k\nRc vdd c 1k\nQ1 c b 0 QN\n.model QN NPN(IS=1e-15 BF=150 BR=2)\nVj jg 0 -0.5\nRj vdd jd 3k\nJ1 jd jg 0 JN\n.model JN NJF(VTO=-2 BETA=1.5m LAMBDA=0.01)\nRs vdd s 100\nS1 s 0 g 0 SW1\n.model SW1 SW(RON=10 ROFF=1MEG VT=1 VH=0.2)\nVp vdd2 0 -5\nRp vdd2 pd 2k\nM2 pd pg 0 0 PM\nVpg pg 0 -2\n.model PM PMOS(VTO=-0.7 KP=1m)\n";
const LIN: &str = "* linear\n.param rv=2k gain=4\nV1 in 0 DC 10\nR1 in a {rv}\nR2 a 0 3k\nE1 b 0 a 0 {gain}\nR3 b c 1k\nVsense c 0 0\nF1 0 d Vsense 2\nR4 d 0 500\nH1 e 0 Vsense 100\nR5 e 0 1k\nG1 0 f a 0 1m\nR6 f 0 2k\nI1 0 g 1m\nR7 g 0 4.7k\n";
const RC: &str = "V1 in 0 PULSE(0 5 1u 1u 1u 5u 20u)\nR1 in out 1k\nC1 out 0 2n\nL1 out x 10m\nR2 x 0 2k\nI1 0 x SIN(0 1m 50k 2u 1e3 30)\n.ic V(out)=1\n";
const RLC: &str = "V1 in 0 0 AC 1 30\nR1 in a 50\nL1 a b 1m\nC1 b 0 1u\nL2 c 0 2m\nK1 L1 L2 0.6\nR2 c 0 75\nI2 0 b 0 AC 1m\n";

#[test]
fn numbers_and_expressions() {
    let mut out = Vec::new();
    for t in [
        "1k", "2.2MEG", "3.3m", "10u", "47n", "1p", "5f", "2T", "1.5G", "1e3", " 4.7K ", "1MIL",
        "0.5e-3k",
    ] {
        out.push(parse_spice_number(t).unwrap());
    }
    let p: HashMap<String, f64> = [("rload".to_string(), 2.0e3), ("Gain".to_string(), 3.0)]
        .into_iter()
        .collect();
    for e in [
        "1+2*3",
        "{(1+2)*3}",
        "2^3^2",
        "-2^2",
        "rload/2",
        "GAIN*1k+1MEG",
        "1e-3*2",
        "((4))",
        "2.5u*4",
    ] {
        out.push(evaluate_spice_expression(e, &p).unwrap());
    }
    close("numbers", &out, 1e-15); // serde_json may round 17-digit text by an ulp
    let e = HashMap::new();
    s(
        "expr_errors",
        &[
            err(evaluate_spice_expression("1/0", &e)),
            err(evaluate_spice_expression("foo+1", &e)),
            err(evaluate_spice_expression("(1+2", &e)),
            err(evaluate_spice_expression("1 2", &e)),
        ]
        .join("|"),
    );
}

#[test]
fn dc_linear_and_nonlinear() {
    let lin = parse(LIN);
    let mut v = Vec::new();
    op(&mut v, &lin.dc_operating_point(&cfg()).unwrap());
    op(&mut v, &homotopy(&lin));
    close("linear_dc", &v, 1e-13);
    let nl = parse(NL);
    let mut v = Vec::new();
    op(&mut v, &nl.dc_operating_point(&cfg()).unwrap());
    op(&mut v, &homotopy(&nl));
    close("nonlinear_dc", &v, 1e-11);
    let mut v = Vec::new();
    for f in [1.0, 1e3, 1e6] {
        acs(&mut v, &nl.ac(f, None, &cfg()).unwrap());
    }
    close("nonlinear_ac", &v, 1e-10);
}

#[test]
fn ac_and_transient() {
    let rlc = parse(RLC);
    let mut v = Vec::new();
    for a in rlc.ac_log_sweep(100.0, 1e5, 7, &cfg()).unwrap() {
        acs(&mut v, &a);
    }
    close("rlc_ac", &v, 1e-11);
    let rc = parse(RC);
    for (name, m) in [
        ("tr_be", TransientMethod::BackwardEuler),
        ("tr_trap", TransientMethod::Trapezoidal),
        ("tr_bdf2", TransientMethod::Bdf2),
    ] {
        let mut v = Vec::new();
        tr(&mut v, &rc.transient(0.5e-6, 40, m, &cfg()).unwrap());
        close(name, &v, 1e-11);
    }
    let xf = parse("V1 in 0 SIN(0 1 10k)\nR1 in a 10\nL1 a 0 1m\nL2 b 0 4m\nK12 L1 L2 0.9\nR2 b 0 100\nC2 b 0 10n\n");
    let mut v = Vec::new();
    tr(
        &mut v,
        &xf.transient(2e-6, 30, TransientMethod::Trapezoidal, &cfg())
            .unwrap(),
    );
    tr(
        &mut v,
        &xf.transient(2e-6, 30, TransientMethod::Bdf2, &cfg())
            .unwrap(),
    );
    close("tr_mutual", &v, 1e-10);
    let rect = parse("V1 in 0 PWL(0 0 10u 5 20u -5 30u 5)\nD1 in out DX\n.model DX D(IS=1e-12)\nR1 out 0 1k\nC1 out 0 100n\n");
    let mut v = Vec::new();
    tr(
        &mut v,
        &rect
            .transient(1e-6, 30, TransientMethod::BackwardEuler, &cfg())
            .unwrap(),
    );
    tr(
        &mut v,
        &rect
            .transient(1e-6, 30, TransientMethod::Trapezoidal, &cfg())
            .unwrap(),
    );
    close("tr_rectifier", &v, 1e-10);
    for (name, m) in [
        ("adapt_be", TransientMethod::BackwardEuler),
        ("adapt_bdf2", TransientMethod::Bdf2),
    ] {
        let ac = AdaptiveTransientConfig {
            method: m,
            initial_step: 1e-7,
            maximum_step: 2e-6,
            relative_tolerance: 1e-3,
            absolute_tolerance: 1e-5,
            ..Default::default()
        };
        let r = rc.transient_adaptive(1.5e-5, &ac, &cfg()).unwrap();
        let mut v = vec![
            r.accepted_steps as f64,
            r.rejected_steps as f64,
            r.minimum_accepted_step,
            r.maximum_accepted_step,
            r.points.len() as f64,
        ];
        tr(&mut v, &r.points);
        close(name, &v, 1e-9);
    }
    let ac = AdaptiveTransientConfig {
        method: TransientMethod::Bdf2,
        initial_step: 1e-6,
        ..Default::default()
    };
    let r = rect.transient_adaptive(3e-5, &ac, &cfg()).unwrap();
    let mut v = vec![r.accepted_steps as f64, r.rejected_steps as f64];
    tr(&mut v, &r.points);
    close("adapt_rect", &v, 1e-9);
}

#[test]
fn noise_and_hierarchy() {
    let nl = parse(NL);
    let lin = parse(LIN);
    let n1 = nl.output_noise(1e3, "a", 300.0, &cfg()).unwrap();
    let n2 = lin.output_noise(1e4, "b", 350.0, &cfg()).unwrap();
    let mut v = vec![
        n1.output_noise_v_per_sqrt_hz,
        n1.output_noise_density_v2_per_hz,
    ];
    v.extend(
        n1.contributions
            .iter()
            .map(|c| c.output_noise_density_v2_per_hz),
    );
    v.push(n2.output_noise_density_v2_per_hz);
    v.extend(
        n2.contributions
            .iter()
            .map(|c| c.output_noise_density_v2_per_hz),
    );
    close("noise", &v, 1e-9);
    let names: String = n1
        .contributions
        .iter()
        .map(|c| format!("{},", c.source))
        .collect();
    s("noise_names", &names);

    let hier = parse(".param vin=3 rbase=1k\n.func half(x) {x/2}\n.func twice(x) {2*x}\n.subckt divider top bot mid params: ra=1k rb=1k\nR1 top mid {ra}\nR2 mid bot {rb}\n.ends\n.subckt stage in out params: g=2\n.model DS D(IS=2e-14)\nX1 in 0 m divider ra={g*1k} rb=2k\nE1 out 0 m 0 {twice(g)}\nD1 out 0 DS\nR1 out 0 10k\n.ends\nV1 in 0 {vin}\nXs in out stage g=1.5\nB1 bo 0 V={half(V(in))+V(out,0)*0.1+TIME}\nRbo bo 0 1k\nGp 0 po POLY(1) in 0 1m 0.5m 0.1m\nRpo po 0 {rbase}\nEt to 0 TABLE {V(in)} = (0,0) (1,2) (4,5)\nRto to 0 1k\nBi 0 bi I={V(in)*1m}\nRbi bi 0 2k\n.nodeset V(out)=0.6\n.ic V(bo)=0.5\n");
    let mut v = Vec::new();
    op(&mut v, &hier.dc_operating_point(&cfg()).unwrap());
    op(&mut v, &homotopy(&hier));
    close("hier_dc", &v, 1e-10);
    let nodes: String = hier.node_names().iter().map(|n| format!("{n},")).collect();
    s("hier_nodes", &nodes);
    let mut v = Vec::new();
    acs(&mut v, &hier.ac(2e3, None, &cfg()).unwrap());
    close("hier_ac_tr", &v, 1e-9);
    s(
        "hier_tr_msg",
        &err(hier.transient(1e-3, 5, TransientMethod::BackwardEuler, &cfg())),
    );
    s(
        "parse_errors",
        &[
            err(Circuit::parse_spice(".include x.cir\n")),
            err(Circuit::parse_spice("X1 a b nosuch\n")),
            err(Circuit::parse_spice(".subckt a b\nR1 b 0 1\n")),
            err(Circuit::parse_spice("B1 a 0 Q=1\n")),
            err(Circuit::parse_spice(".ic V(a,b)=1\nR1 a 0 1\n")),
        ]
        .join("|"),
    );
}

#[test]
fn nports_and_compact_devices() {
    let mut c = Circuit::new();
    let (inn, out) = (c.node("in"), c.node("out"));
    c.add_voltage_source("V1", inn, 0, 0.0, Complex::new(1.0, 0.0), None)
        .unwrap();
    let p1 = c.node("p1");
    c.add_resistor("Rs", inn, p1, 50.0).unwrap();
    let freqs = [1e8, 1e9, 3e9];
    c.add_microstrip_line(
        "TL1",
        (p1, 0),
        (out, 0),
        1.9e-3,
        1e-3,
        4.4,
        0.05,
        &freqs,
        50.0,
    )
    .unwrap();
    c.add_resistor("Rl", out, 0, 75.0).unwrap();
    let mut v = Vec::new();
    for f in [1e8, 5e8, 2e9, 4e9] {
        acs(&mut v, &c.ac(f, None, &cfg()).unwrap());
    }
    let mut w = Circuit::new();
    let (a, b) = (w.node("a"), w.node("b"));
    w.add_current_source("I1", 0, a, 0.0, Complex::new(1e-3, 0.0), None)
        .unwrap();
    w.add_resistor("Ra", a, 0, 50.0).unwrap();
    w.add_resistor("Rb", b, 0, 50.0).unwrap();
    let wf = [8e9, 1e10, 1.2e10];
    w.add_rectangular_waveguide_te10_line(
        "WG",
        (a, 0),
        (b, 0),
        22.86e-3,
        10.16e-3,
        1.0,
        1.0,
        0.1,
        &wf,
        50.0,
    )
    .unwrap();
    w.add_stripline_line("SL", (a, 0), (b, 0), 1e-3, 2e-3, 2.2, 0.02, &wf, 50.0)
        .unwrap();
    w.add_coplanar_waveguide_line("CPW", (a, 0), (b, 0), 0.5e-3, 0.2e-3, 9.8, 0.01, &wf, 50.0)
        .unwrap();
    for f in [8e9, 9e9, 1.3e10] {
        acs(&mut v, &w.ac(f, None, &cfg()).unwrap());
    }
    close("nports", &v, 1e-9);

    let mut c = Circuit::new();
    let (a, b, th) = (c.node("a"), c.node("b"), c.node("th"));
    let pulse = PulseWaveform {
        initial: 0.0,
        pulsed: 2.0,
        delay: 0.0,
        rise: 1e-6,
        fall: 1e-6,
        width: 5e-6,
        period: 0.0,
    };
    c.add_voltage_source(
        "V1",
        a,
        0,
        2.0,
        Complex::new(0.1, 0.0),
        Some(pulse.waveform()),
    )
    .unwrap();
    c.add_resistor("R1", a, b, 100.0).unwrap();
    c.add_static_device(
        "NL",
        vec![b, 0],
        Arc::new(|v: &[f64]| {
            let x = v[0] - v[1];
            let i = 1e-3 * x + 2e-4 * x * x * x;
            let gd = 1e-3 + 6e-4 * x * x;
            Ok(StaticDeviceEvaluation {
                terminal_current: vec![i, -i],
                jacobian: vec![gd, -gd, -gd, gd],
            })
        }),
    )
    .unwrap();
    c.add_dynamic_device(
        "CQ",
        vec![b, 0],
        Arc::new(|v: &[f64]| {
            let x = v[0] - v[1];
            let q = 1e-9 * x + 2e-10 * x * x;
            let cq = 1e-9 + 4e-10 * x;
            Ok(DynamicDeviceEvaluation {
                terminal_current: vec![1e-6 * x, -1e-6 * x],
                current_jacobian: vec![1e-6, -1e-6, -1e-6, 1e-6],
                terminal_charge: vec![q, -q],
                charge_jacobian: vec![cq, -cq, -cq, cq],
            })
        }),
    )
    .unwrap();
    let et = ElectroThermalDeviceConfig {
        thermal_resistance_k_per_w: 50.0,
        thermal_capacitance_j_per_k: 1e-6,
        ..Default::default()
    };
    c.add_electrothermal_device(
        "ET",
        vec![b, 0],
        th,
        et,
        Arc::new(|v: &[f64], t: f64| {
            let x = v[0] - v[1];
            let g0 = 2e-3 * (1.0 - 0.004 * (t - 300.0));
            let i = g0 * x;
            let dg = -2e-3 * 0.004 * x;
            Ok(ElectroThermalDeviceEvaluation {
                electrical_current: vec![i, -i],
                electrical_jacobian: vec![g0, -g0, -g0, g0],
                electrical_temperature_derivative: vec![dg, -dg],
                dissipated_power_w: i * x,
                power_voltage_derivative: vec![2.0 * g0 * x, -2.0 * g0 * x],
                power_temperature_derivative_w_per_k: dg * x,
            })
        }),
    )
    .unwrap();
    let mut v = Vec::new();
    op(&mut v, &c.dc_operating_point(&cfg()).unwrap());
    acs(&mut v, &c.ac(1e5, None, &cfg()).unwrap());
    tr(
        &mut v,
        &c.transient(1e-6, 12, TransientMethod::Bdf2, &cfg())
            .unwrap(),
    );
    tr(
        &mut v,
        &c.transient(1e-6, 8, TransientMethod::BackwardEuler, &cfg())
            .unwrap(),
    );
    let ac = AdaptiveTransientConfig {
        method: TransientMethod::Bdf2,
        initial_step: 2e-7,
        ..Default::default()
    };
    let r = c.transient_adaptive(4e-6, &ac, &cfg()).unwrap();
    v.push(r.accepted_steps as f64);
    v.push(r.rejected_steps as f64);
    tr(&mut v, &r.points);
    close("devices", &v, 1e-9);
    s(
        "trap_dynamic_err",
        &err(c.transient(1e-6, 2, TransientMethod::Trapezoidal, &cfg())),
    );
}

#[test]
fn analysis_helpers() {
    let nl = parse(NL);
    let mut v = Vec::new();
    let mut copy = nl.clone();
    for p in copy
        .dc_sweep_voltage_source("Vg", 0.0, 3.0, 4, &cfg())
        .unwrap()
    {
        v.push(p.source_value);
        op(&mut v, &p.operating_point);
    }
    for p in nl.dc_temperature_sweep(250.0, 400.0, 3, &cfg()).unwrap() {
        v.push(p.temperature_k);
        op(&mut v, &p.operating_point);
    }
    let sens = voltage_source_dc_sensitivity(nl.clone(), "Vg", "dr", 1e-3, &cfg()).unwrap();
    v.extend([sens.baseline_output, sens.derivative]);
    let d = small_signal_distortion(nl.clone(), "Vb", "c", 1e-3, 1e-3, &cfg()).unwrap();
    v.extend([
        d.bias_output_v,
        d.first_derivative,
        d.second_derivative,
        d.third_derivative,
        d.fundamental_v,
        d.second_harmonic_v,
        d.third_harmonic_v,
        d.hd2,
        d.hd3,
    ]);
    let tol = [
        ResistorTolerance {
            element: "R1".into(),
            relative_sigma: 0.05,
        },
        ResistorTolerance {
            element: "Rd".into(),
            relative_sigma: 0.1,
        },
    ];
    let mc = monte_carlo_dc_resistors(&nl, "dr", &tol, 25, 7, &cfg()).unwrap();
    v.extend([mc.mean, mc.standard_deviation, mc.minimum, mc.maximum]);
    for x in [0.9, 0.2, -0.3] {
        v.push(limit_pn_junction_voltage(x, 0.6, 0.025852, 1e-14).unwrap());
    }
    v.push(limit_pn_junction_voltage(5.0, -1.0, 0.025852, 1e-14).unwrap());
    for p in parameter_sweep_dc(
        &nl,
        &[100.0, 200.0],
        &|c, x| c.set_resistance("R1", x),
        &cfg(),
    )
    .unwrap()
    {
        op(&mut v, &p.operating_point);
    }
    for p in temperature_sweep_dc(&nl, &[280.0, 320.0], &cfg()).unwrap() {
        op(&mut v, &p.operating_point);
    }
    // The distortion third derivative is a cancellation-dominated stencil: loose tolerance.
    close("analysis_dc", &v, 1e-6);

    let pi = parse("R1 a 0 100\nR2 a b 50\nC1 b 0 1n\nL1 b c 1u\nR3 c 0 200\n");
    let mut v = Vec::new();
    let sp = two_port_s_parameters(&pi, &AcPort::new("a"), &AcPort::new("c"), 1e7, 50.0).unwrap();
    for c in [sp.a11, sp.a12, sp.a21, sp.a22] {
        cx(&mut v, c);
    }
    let ports = [
        CircuitPort {
            positive: pi.find_node("a").unwrap(),
            negative: 0,
            reference_impedance: 50.0,
        },
        CircuitPort {
            positive: pi.find_node("c").unwrap(),
            negative: 0,
            reference_impedance: 75.0,
        },
    ];
    let e = extract_s_parameters(&pi, &ports, 1e7, &cfg()).unwrap();
    for i in 0..2 {
        for j in 0..2 {
            cx(&mut v, e[(i, j)]);
        }
    }
    let lp = parse("V1 in 0 0 AC 1\nR1 in a 1k\nC1 a 0 100n\nR2 a b 10k\nC2 b 0 1n\n");
    let dp = estimate_dominant_pole(&lp, "b", 10.0, 1e6, 41, &cfg()).unwrap();
    v.extend([
        dp.dc_gain,
        dp.pole_frequency_hz,
        f64::from(u8::from(dp.found)),
    ]);
    let pz = PoleZeroConfig {
        start_hz: 10.0,
        stop_hz: 1e7,
        samples: 40,
        numerator_order: 0,
        ..Default::default()
    };
    let r = pole_zero_analysis(&lp, "b", &pz, &cfg()).unwrap();
    for &p in &r.poles_rad_per_s {
        cx(&mut v, p);
    }
    for &z in &r.zeros_rad_per_s {
        cx(&mut v, z);
    }
    v.extend([
        r.relative_rms_fit_error,
        r.poles_rad_per_s.len() as f64,
        r.zeros_rad_per_s.len() as f64,
    ]);
    let ns = output_noise_spectrum(&nl, "a", 10.0, 1e5, 5, 300.0, &cfg()).unwrap();
    v.push(ns.integrated_rms_voltage);
    v.extend(ns.points.iter().map(|p| p.output_noise_density_v2_per_hz));
    close("analysis_ac", &v, 1e-6);

    let amp = parse("V1 in 0 SIN(0 0.5 1k)\nR1 in a 1k\nD1 a 0 DD\nD2 0 a DD\n.model DD D(IS=1e-12)\nR2 a 0 10k\nC1 a 0 10n\n");
    let t = amp
        .transient(1e-5, 200, TransientMethod::Trapezoidal, &cfg())
        .unwrap();
    let a = amp.find_node("a").unwrap();
    let mut v = Vec::new();
    for f in [0.0, 1e3, 2e3, 3e3] {
        let m = fourier_measurement(&t, a, f).unwrap();
        cx(&mut v, m.phasor);
        v.extend([m.magnitude, m.phase_rad]);
    }
    v.push(total_harmonic_distortion(&t, a, 1e3, 5).unwrap());
    let pss = PeriodicSteadyStateConfig {
        period_s: 1e-3,
        samples_per_period: 40,
        max_periods: 8,
        relative_tolerance: 1e-2,
        absolute_tolerance: 1e-4,
        ..Default::default()
    };
    let p = periodic_steady_state(&amp, &pss, &cfg()).unwrap();
    v.extend([
        p.periods as f64,
        p.normalized_residual,
        f64::from(u8::from(p.converged)),
    ]);
    tr(&mut v, &p.period);
    close("analysis_tr", &v, 1e-8);
}

#[test]
fn sparse_raw_and_plans() {
    let mut net = String::from("V1 n0 0 1\n");
    for i in 0..80 {
        net += &format!(
            "R{i} n{i} n{} {}\nRg{i} n{} 0 {}\n",
            i + 1,
            10 + i,
            i + 1,
            1000 + 7 * i
        );
    }
    let big = parse(&net);
    let mut v = Vec::new();
    op(&mut v, &big.dc_operating_point(&cfg()).unwrap());
    let dense = NewtonConfig {
        sparse_mna_threshold: 0,
        ..cfg()
    };
    op(&mut v, &big.dc_operating_point(&dense).unwrap());
    // GMRES residual stopping rules differ slightly from upstream's: compare to 1e-9.
    close("sparse", &v, 1e-9);

    let rc = parse(RC);
    let t = rc
        .transient(1e-6, 4, TransientMethod::BackwardEuler, &cfg())
        .unwrap();
    let ds = raw_from_transient(&t, rc.node_names(), "cfd_solvers transient").unwrap();
    let text = spice_raw_text(&ds).unwrap();
    s("raw_tr", &text);
    let rlc = parse(RLC);
    let a = [
        rlc.ac(1e3, None, &cfg()).unwrap(),
        rlc.ac(2e3, None, &cfg()).unwrap(),
    ];
    let da = raw_from_ac(&a, rlc.node_names(), "cfd_solvers AC analysis").unwrap();
    let text = spice_raw_text(&da).unwrap();
    // Complex division rounds differently from libstdc++ in the last bit; compare numerically.
    let cpp = g()["raw_ac"].as_str().unwrap();
    assert_eq!(text.lines().count(), cpp.lines().count());
    for (r, c) in text.lines().zip(cpp.lines()) {
        if r != c {
            let nums = |l: &str| -> Vec<f64> {
                l.trim()
                    .trim_matches(['(', ')'])
                    .split(',')
                    .map(|x| x.parse::<f64>().unwrap())
                    .collect()
            };
            for (x, y) in nums(r).iter().zip(nums(c)) {
                assert!((x - y).abs() <= 1e-12 * y.abs().max(1e-3), "{r} vs {c}");
            }
        }
    }
    let back = read_spice_raw(&text).unwrap();
    let mut v = vec![
        back.traces.len() as f64,
        f64::from(u8::from(back.complex_values)),
    ];
    for tr in &back.traces {
        for &c in &tr.values {
            cx(&mut v, c);
        }
    }
    close("raw_back", &v, 1e-12);

    let p = make_openvaf_compile_plan("bsim.va", "bsim.osdi", "openvaf").unwrap();
    let n = make_external_spice_run_plan(
        ExternalSpiceEngine::Ngspice,
        "a.cir",
        "a.raw",
        &["-x".into()],
    )
    .unwrap();
    let x = make_external_spice_run_plan(
        ExternalSpiceEngine::Xyce,
        "a.cir",
        "a.raw",
        &["-r".into(), "a.raw".into()],
    )
    .unwrap();
    let mut out = String::new();
    for a in &p.argv {
        out += &format!("{a} ");
    }
    out += "|";
    for a in &n.argv {
        out += &format!("{a} ");
    }
    out += "|";
    for a in &x.argv {
        out += &format!("{a} ");
    }
    out += "|";
    for d in spice_compatibility_profile(SpiceCompatibilityMode::Ltspice).supported_directives {
        out += &format!("{d} ");
    }
    s("plans", &out);
}
