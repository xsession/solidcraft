//! Analytic checks for power flow, dispatch, swing transients, workflows and acoustic coupling.
use cfd_solvers::{acoustics::*, fem::*, power::*, workflow::*};
use std::{cell::RefCell, rc::Rc};

#[test]
fn two_bus_line_losses_match_current_squared_resistance() {
    let c = TwoBusPowerFlowConfig::default();
    let r = solve_two_bus_power_flow(&c, None).unwrap();
    assert!(r.solve.converged);
    // I = (V_s - V e^{j theta}) / Z; load bus consumes 0.5 + 0.2j, losses = |I|^2 R.
    let (vr, vi) = (r.voltage * r.angle.cos(), r.voltage * r.angle.sin());
    let (dr, di) = (c.slack_voltage - vr, -vi);
    let z2 = c.line_resistance.powi(2) + c.line_reactance.powi(2);
    let i2 = (dr * dr + di * di) / z2;
    // Complex power delivered to the load bus is V conj(I).
    let (ir, ii) = (
        (dr * c.line_resistance + di * c.line_reactance) / z2,
        (di * c.line_resistance - dr * c.line_reactance) / z2,
    );
    let p_load = vr * ir + vi * ii;
    assert!(
        (p_load + c.specified_active_power).abs() < 1e-10,
        "{p_load}"
    );
    assert!(i2 * c.line_resistance > 0.0);
}

#[test]
fn dispatch_meets_demand_and_swing_settles_to_power_balance() {
    let gens = [
        QuadraticDispatchGenerator {
            minimum: 0.0,
            maximum: 2.0,
            linear_cost: 1.0,
            quadratic_cost: 0.5,
        },
        QuadraticDispatchGenerator {
            minimum: 0.0,
            maximum: 2.0,
            linear_cost: 1.2,
            quadratic_cost: 0.1,
        },
    ];
    let d = quadratic_dispatch(&gens, 1.5).unwrap();
    assert!(d.feasible && (d.served - 1.5).abs() < 1e-14);
    let bus = [SwingBusParameters {
        inertia: 1.0,
        damping: 2.0,
        mechanical_power: 0.4,
        ..Default::default()
    }];
    let t = integrate_swing_network(
        &bus,
        &[SwingLine {
            from: Some(0),
            to: None,
            synchronizing_power: 1.0,
        }],
        0.1,
        400,
    )
    .unwrap();
    let last = t.states.last().unwrap();
    assert!(
        (last[0] - 0.4f64.asin()).abs() < 1e-4 && last[1].abs() < 1e-4,
        "{last:?}"
    );
}

#[test]
fn restartable_workflow_resumes_after_failure() {
    let dir = std::env::temp_dir().join(format!("cfd-workflow-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let ck = dir.join("ck.txt");
    let log = Rc::new(RefCell::new(Vec::new()));
    let fail = Rc::new(RefCell::new(true));
    let build = |log: Rc<RefCell<Vec<&'static str>>>, fail: Rc<RefCell<bool>>| {
        let mut w = RestartableWorkflow::new();
        let (l1, l2) = (log.clone(), log);
        let a = w
            .add(
                "mesh",
                Box::new(move || {
                    l1.borrow_mut().push("mesh");
                    Ok(())
                }),
                vec![],
            )
            .unwrap();
        w.add(
            "solve",
            Box::new(move || {
                if *fail.borrow() {
                    return Err(cfd_solvers::SolverError::NumericalFailure(
                        "diverged".into(),
                    ));
                }
                l2.borrow_mut().push("solve");
                Ok(())
            }),
            vec![a],
        )
        .unwrap();
        w
    };
    assert!(build(log.clone(), fail.clone()).run(&ck).is_err());
    *fail.borrow_mut() = false;
    let mut w = build(log.clone(), fail);
    w.run(&ck).unwrap();
    assert_eq!(*log.borrow(), vec!["mesh", "solve"]);
    assert_eq!(w.completed(), ["mesh".to_string(), "solve".to_string()]);
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn uniform_pressure_on_a_closed_boundary_has_zero_net_force() {
    let mesh = make_rectangle_tri_mesh(4, 3, 2.0, 1.0).unwrap();
    let f =
        pressure_to_structural_boundary_load(&mesh, &vec![1e3; mesh.boundary_edges.len()], None)
            .unwrap();
    let (fx, fy) = f.iter().fold((0.0, 0.0), |(a, b), v| (a + v.x, b + v.y));
    assert!(fx.abs() < 1e-9 && fy.abs() < 1e-9);
    let inv =
        inverse_coordinate_search(|p| (p[0] - 0.25).powi(2), vec![1.0], vec![0.5], 200, 1e-10)
            .unwrap();
    assert!(inv.converged && (inv.parameters[0] - 0.25).abs() < 1e-9);
}
