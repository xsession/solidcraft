//! Value-by-value comparison with the C++ upstream hybrid DAE graph, power-flow, workflow
//! runner and acoustic coupling units. Goldens come from `parity/oracle_core_power.cpp`.
use cfd_solvers::{acoustics::*, core::*, fem::*, multiphysics::*, power::*, workflow::*};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-core-power.json")).unwrap()
}

#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    let floor = 1e-6 * cpp.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= rel * c.abs().max(floor) + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn flat(t: &DaeTrajectory) -> Vec<f64> {
    let mut o = vec![
        t.times.len() as f64,
        t.events.len() as f64,
        f64::from(u8::from(t.converged)),
    ];
    for (time, s) in t.times.iter().zip(&t.states) {
        o.push(*time);
        o.extend(s);
    }
    o
}

#[test]
fn hybrid_dae_graph_matches_upstream() {
    let mut g = HybridDaeGraph::new();
    let x = g.add_variable("x", 1.0, true).unwrap();
    let v = g.add_variable("v", 0.0, true).unwrap();
    let y = g.add_variable("y", 3.0, false).unwrap();
    g.add_component(
        "osc",
        &[x, v],
        Box::new(move |_, s, d, o| {
            o[0] = d[x] - s[v];
            o[1] = d[v] + 0.4 * s[v] + 3.0 * s[x];
            Ok(())
        }),
    )
    .unwrap();
    g.add_connection("mirror", y, y, x, 2.0, 1.0).unwrap();
    assert!(g.add_connection("again", y, y, x, 2.0, 1.0).is_err());
    let value = |_: f64, s: &[f64]| s[0];
    let events = [DaeEvent {
        name: "x_zero".into(),
        value: &value,
    }];
    let t = g
        .integrate(0.0, 0.05, 30, &events, MatrixFreeNewtonConfig::default())
        .unwrap();
    // Newton stops at 1e-10 absolute residual, so states agree to about 1e-9.
    close("dae_graph", &flat(&t), 1e-8, 1e-10);
    assert_eq!(t.events, vec!["x_zero".to_string()]);
}

#[test]
fn power_flow_matches_upstream() {
    let r2 = solve_two_bus_power_flow(&TwoBusPowerFlowConfig::default(), None).unwrap();
    let r3 = solve_two_bus_power_flow(
        &TwoBusPowerFlowConfig {
            line_resistance: 0.05,
            specified_active_power: -0.9,
            ..Default::default()
        },
        Some([0.95, -0.1]),
    )
    .unwrap();
    let v = [
        r2.voltage,
        r2.angle,
        r2.solve.iterations as f64,
        r2.solve.residual_rms,
        r3.voltage,
        r3.angle,
        r3.solve.iterations as f64,
    ];
    close("two_bus", &v, 1e-9, 1e-13);
    let bus = |kind, voltage, p, q| PowerBus {
        kind,
        voltage,
        angle: 0.0,
        specified_active: p,
        specified_reactive: q,
    };
    let line = |from, to, resistance, reactance| PowerLine {
        from,
        to,
        resistance,
        reactance,
    };
    let net = MultiBusPowerFlow::new(
        vec![
            bus(PowerBusType::Slack, 1.02, 0.0, 0.0),
            bus(PowerBusType::Pq, 1.0, -0.4, -0.1),
            bus(PowerBusType::Pq, 1.0, -0.3, -0.15),
            bus(PowerBusType::Pq, 1.0, 0.1, 0.02),
        ],
        vec![
            line(0, 1, 0.01, 0.08),
            line(1, 2, 0.02, 0.1),
            line(0, 2, 0.015, 0.09),
            line(2, 3, 0.01, 0.05),
        ],
    )
    .unwrap();
    let mb = net.solve(1.0, None).unwrap();
    let mut mv = mb.state.clone();
    mv.extend([
        mb.solve.iterations as f64,
        f64::from(u8::from(mb.solve.converged)),
    ]);
    close("multibus", &mv, 1e-8, 1e-10);
    let path: Vec<f64> = net
        .continuation(0.5, 2.0, 4)
        .unwrap()
        .iter()
        .flat_map(|p| std::iter::once(p.load_factor).chain(p.result.state.iter().copied()))
        .collect();
    close("continuation", &path, 1e-8, 1e-10);
    let sb = [
        SwingBusParameters {
            inertia: 2.0,
            damping: 0.3,
            mechanical_power: 0.5,
            initial_angle: 0.0,
            initial_speed: 0.0,
        },
        SwingBusParameters {
            inertia: 1.5,
            damping: 0.2,
            mechanical_power: -0.2,
            initial_angle: 0.1,
            initial_speed: 0.0,
        },
    ];
    let sl = [
        SwingLine {
            from: Some(0),
            to: Some(1),
            synchronizing_power: 1.2,
        },
        SwingLine {
            from: Some(0),
            to: None,
            synchronizing_power: 0.8,
        },
        SwingLine {
            from: Some(1),
            to: None,
            synchronizing_power: 0.5,
        },
    ];
    close(
        "swing",
        &flat(&integrate_swing_network(&sb, &sl, 0.05, 20).unwrap()),
        1e-6,
        1e-8,
    );
    let gens = [
        QuadraticDispatchGenerator {
            minimum: 0.1,
            maximum: 1.0,
            linear_cost: 2.0,
            quadratic_cost: 0.5,
        },
        QuadraticDispatchGenerator {
            minimum: 0.0,
            maximum: 0.8,
            linear_cost: 1.5,
            quadratic_cost: 1.0,
        },
        QuadraticDispatchGenerator {
            minimum: 0.2,
            maximum: 0.5,
            linear_cost: 3.0,
            quadratic_cost: 0.1,
        },
    ];
    let (d1, d2) = (
        quadratic_dispatch(&gens, 1.4).unwrap(),
        quadratic_dispatch(&gens, 5.0).unwrap(),
    );
    let mut dv = d1.generation.clone();
    dv.extend([
        d1.served,
        d1.objective,
        f64::from(u8::from(d1.feasible)),
        f64::from(u8::from(d2.feasible)),
    ]);
    close("dispatch", &dv, 1e-14, 1e-15);
}

#[test]
fn workflow_and_acoustic_coupling_match_upstream() {
    let mut w: Vec<f64> = parameter_sweep(&[0.0, 0.5, 2.0], |p| p * p - 1.0)
        .iter()
        .flat_map(|s| [s.parameter, s.objective])
        .collect();
    let inv = inverse_coordinate_search(
        |p| (p[0] - 1.3).powi(2) + 2.0 * (p[1] + 0.7).powi(2) + 0.1 * p[0] * p[1],
        vec![0.0, 0.0],
        vec![0.5, 0.5],
        200,
        1e-6,
    )
    .unwrap();
    w.extend(&inv.parameters);
    w.extend([
        inv.objective,
        inv.iterations as f64,
        f64::from(u8::from(inv.converged)),
    ]);
    close("workflow", &w, 1e-14, 1e-16);

    let mesh = make_rectangle_tri_mesh(3, 2, 1.5, 1.0).unwrap();
    let ep: Vec<f64> = (0..mesh.boundary_edges.len())
        .map(|e| 100.0 + 10.0 * e as f64)
        .collect();
    let mut loads = vec![];
    for patch in [None, Some(2)] {
        loads.extend(
            pressure_to_structural_boundary_load(&mesh, &ep, patch)
                .unwrap()
                .iter()
                .flat_map(|f| [f.x, f.y]),
        );
    }
    close("struct_load", &loads, 1e-14, 1e-12);
    let mut bar = DynamicBar1D::new(DynamicBar1DConfig {
        elements: 6,
        dt: 1e-6,
        ..Default::default()
    })
    .unwrap();
    bar.initialize(|_| 0.0, None).unwrap();
    drive_dynamic_bar_from_pressure_trace(&mut bar, &[1e5, 2e5, -5e4, 0.0], 1e-3).unwrap();
    // Each Newmark solve stops at the 1e-11 linear tolerance; the trace is compared absolutely.
    close("bar", bar.displacement(), 1e-6, 1e-18);
    let mut pen = PennesBioheat2D::new(PennesBioheat2DConfig {
        nx: 4,
        ny: 3,
        dt_s: 0.5,
        ..Default::default()
    })
    .unwrap();
    pen.initialize(310.0).unwrap();
    let inten: Vec<f64> = (0..12).map(|i| 1e4 * (1 + i % 3) as f64).collect();
    let alpha: Vec<f64> = (0..12).map(|i| 5.0 + i as f64).collect();
    apply_acoustic_heating_to_pennes(&mut pen, &inten, &alpha).unwrap();
    pen.step(3).unwrap();
    close("pennes", pen.temperature_k(), 1e-12, 0.0);
    let st: Vec<f64> = acoustic_streaming_acceleration(
        &inten,
        &alpha,
        &[1000.0; 12],
        &[1540.0; 12],
        Vec3::new(1.0, 2.0, 2.0),
    )
    .unwrap()
    .iter()
    .flat_map(|a| [a.x, a.y, a.z])
    .collect();
    close("streaming", &st, 1e-14, 0.0);
    let measured = vec![vec![0.0, 1.0, 0.5, -0.2], vec![0.3, 0.2, 0.1, 0.0]];
    let sim = |p: &[f64]| {
        vec![
            vec![0.0, p[0], 0.5 * p[0], -0.2 * p[1]],
            vec![0.3 * p[1], 0.2, 0.1 * p[0] * p[1], 0.0],
        ]
    };
    let mut fit = vec![
        sensor_trace_l2_misfit(&sim(&[0.8, 1.1]), &measured, true).unwrap(),
        sensor_trace_l2_misfit(&sim(&[0.8, 1.1]), &measured, false).unwrap(),
    ];
    let fr =
        fit_acoustic_parameters(sim, &measured, vec![0.2, 0.4], vec![0.3, 0.3], 100, 1e-7).unwrap();
    fit.extend(&fr.parameters);
    fit.extend([
        fr.misfit,
        fr.iterations as f64,
        f64::from(u8::from(fr.converged)),
    ]);
    close("fit", &fit, 1e-12, 1e-15);
}
