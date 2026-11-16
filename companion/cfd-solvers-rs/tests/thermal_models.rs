use cfd_solvers::{core::Vec3, fvm::*};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};
#[test]
fn thermal_and_solid_sources_satisfy_energy_balance() {
    let m = make_cartesian_hexa_mesh(5, 3, 2, 1.0, 0.6, 0.2).unwrap();
    let mut t = ThermalTransport::new(
        m.clone(),
        ThermalTransportConfig {
            dt: 1e-3,
            density: 2.0,
            default_heat_capacity: 3.0,
            default_conductivity: 0.1,
        },
    )
    .unwrap();
    t.initialize_constant(300.0).unwrap();
    t.set_volumetric_source(|_, _| 12.0);
    let initial = t.total_sensible_enthalpy().unwrap();
    t.run(10).unwrap();
    assert!(t.temperature().iter().all(|v| (v - 300.02).abs() < 1e-10));
    assert!((t.total_sensible_enthalpy().unwrap() - initial - 12.0 * 0.12 * 0.01).abs() < 1e-9);
    assert!(t.transport_stable_timestep(0.8).unwrap().unwrap() > 1e-3);
    let mut s = SolidHeatConduction::new(
        m,
        SolidHeatConfig {
            dt: 1e-3,
            density: 2.0,
            heat_capacity: 3.0,
            conductivity: 0.1,
            ..Default::default()
        },
    )
    .unwrap();
    s.initialize_constant(300.0).unwrap();
    s.set_volumetric_heat_source_constant(12.0).unwrap();
    let initial = s.total_sensible_energy();
    s.run(10).unwrap();
    assert!((s.total_sensible_energy() - initial - 12.0 * 0.12 * 0.01).abs() < 1e-8);
}
#[test]
fn temperature_dependent_material_failure_does_not_commit_state() {
    let m = make_cartesian_hexa_mesh(3, 2, 1, 1.0, 1.0, 1.0).unwrap();
    let mut t = ThermalTransport::new(m, Default::default()).unwrap();
    t.initialize_constant(300.0).unwrap();
    t.set_material(|p, _| if p.x > 0.5 { -1.0 } else { 1.0 }, |_, _| 1.0);
    assert!(t.step().is_err());
    assert_eq!(t.time(), 0.0);
    assert!(t.temperature().iter().all(|v| *v == 300.0));
}
#[test]
fn radiation_enclosure_balances_and_thermal_source_uses_actual_temperature() {
    let t = [350.0, 300.0, 280.0];
    let e = [0.7, 0.8, 1.0];
    let a = [1.0; 3];
    let f = [0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0];
    let r = solve_gray_radiosity(&t, &e, &a, &f).unwrap();
    assert!(r.net_heat[0] > 0.0 && r.net_heat[2] < 0.0);
    assert!(r.net_heat.iter().sum::<f64>().abs() < 1e-10);
    assert_eq!(
        gray_surface_exchange(300.0, 0.8, 300.0, 0.8, 1.0).unwrap(),
        0.0
    );
    let q = gray_surface_exchange(350.0, 0.7, 290.0, 0.8, 1.0).unwrap();
    assert!(q > 0.0);
    assert!(optically_thin_radiation_source(300.0, 290.0, 0.1).unwrap() < 0.0);
    let mesh = make_cartesian_hexa_mesh(2, 2, 1, 1.0, 1.0, 1.0).unwrap();
    let mut heat = ThermalTransport::new(
        mesh,
        ThermalTransportConfig {
            default_conductivity: 0.0,
            ..Default::default()
        },
    )
    .unwrap();
    heat.initialize_constant(300.0).unwrap();
    heat.set_reactive_radiative_source(
        None,
        Arc::new(OpticallyThinRadiation::new(0.1).unwrap()),
        290.0,
    )
    .unwrap();
    heat.step().unwrap();
    assert!(heat.temperature().iter().all(|v| *v < 300.0));
}
#[test]
fn closure_limits_and_turbulence_shear_are_correct() {
    assert_eq!(
        spalart_allmaras_eddy_viscosity(1e-5, 0.0, 7.1).unwrap(),
        0.0
    );
    assert!((k_epsilon_eddy_viscosity(2.0, 0.5, 0.1, 0.09).unwrap() - 0.45).abs() < 1e-14);
    assert_eq!(
        des_length_scale(0.1, 0.02, 0.65).unwrap(),
        0.013000000000000001
    );
    assert_eq!(euler_euler_drag_source(1.0, 2.0, 3.0, 4.0).unwrap(), 12.0);
    assert_eq!(euler_euler_drag_source(2.0, 1.0, 3.0, 4.0).unwrap(), -12.0);
    let p = linear_mushy_phase_change(305.0, 300.0, 310.0, 2.0, 4.0, 100.0).unwrap();
    assert_eq!(p.liquid_fraction, 0.5);
    assert_eq!(p.effective_cp, 13.0);
    let m = make_cartesian_hexa_mesh(3, 3, 2, 1.0, 1.0, 1.0).unwrap();
    let u: Vec<_> = m
        .cells()
        .iter()
        .map(|c| Vec3::new(2.0 * c.center.y, 0.0, 0.0))
        .collect();
    let s = eddy_viscosity(
        &m,
        &u,
        TurbulenceConfig {
            model: TurbulenceModel::Smagorinsky,
            ..Default::default()
        },
    )
    .unwrap();
    let expected = 0.17_f64.powi(2) * m.cells()[0].volume.cbrt().powi(2) * 2.0;
    assert!(s.iter().all(|v| (v - expected).abs() < 1e-12));
    let w = eddy_viscosity(
        &m,
        &u,
        TurbulenceConfig {
            model: TurbulenceModel::Wale,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(w.iter().all(|v| v.abs() < 1e-15));
}
#[test]
fn vortex_criteria_pressure_forces_and_probes_have_correct_units() {
    let m = make_cartesian_hexa_mesh(3, 3, 2, 1.0, 1.0, 1.0).unwrap();
    let u: Vec<_> = m
        .cells()
        .iter()
        .map(|c| Vec3::new(-2.0 * c.center.y, 2.0 * c.center.x, 0.0))
        .collect();
    let d = velocity_derived_fields(&m, &u).unwrap();
    assert!(d.vorticity.iter().all(|v| (v.z - 4.0).abs() < 1e-11));
    assert!(d.q_criterion.iter().all(|v| (v - 4.0).abs() < 1e-11));
    assert!(d.lambda2.iter().all(|v| (v + 4.0).abs() < 1e-11));
    assert!(d.enstrophy.iter().all(|v| (v - 8.0).abs() < 1e-10));
    let force = pressure_force_on_patch(
        &m,
        &vec![2.0; m.cell_count()],
        m.patch_index("left").unwrap(),
    )
    .unwrap();
    assert!((force - Vec3::new(2.0, 0.0, 0.0)).norm() < 1e-14);
    let field: Vec<_> = (0..m.cell_count()).map(|i| i as f64).collect();
    assert_eq!(
        sample_nearest_cell(&m, &field, m.cells()[5].center).unwrap(),
        5.0
    );
    assert!((wall_y_plus(0.001, 2.0, 1.0, 1e-5).unwrap() - 100.0 * 2.0_f64.sqrt()).abs() < 1e-12);
}
#[test]
fn report_callbacks_preserve_order_and_residual_history() {
    let order = RefCell::new(vec![]);
    let mut manager = FunctionObjectManager::default();
    manager.add(|_| {
        order.borrow_mut().push(1);
        Ok(())
    });
    manager.add(|_| {
        order.borrow_mut().push(2);
        Ok(())
    });
    let report = SolverReport {
        step: 3,
        time: 0.03,
        residuals: BTreeMap::from([("pressure".into(), 1e-8)]),
    };
    manager.execute(&report).unwrap();
    assert_eq!(*order.borrow(), vec![1, 2]);
    let mut history = ResidualHistory::default();
    history.observe(report).unwrap();
    assert_eq!(history.latest("pressure").unwrap(), 1e-8);
    assert!(history.latest("velocity").is_err());
}
