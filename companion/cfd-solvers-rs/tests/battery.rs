//! Analytic and invariant checks for the battery models.
use cfd_solvers::battery::*;

#[test]
fn particle_conserves_lithium_under_constant_flux() {
    let cfg = ParticleConfig::default();
    let mut p = SphericalDiffusionParticle::new(cfg).unwrap();
    let c0 = p.average_concentration();
    let (q, dt, steps) = (3e-6, 7.0, 9);
    for _ in 0..steps {
        p.step(q, dt).unwrap();
    }
    // Outward flux q over the sphere surface removes 3 q t / R per unit volume.
    let expected = c0 - 3.0 * q * dt * steps as f64 / cfg.radius_m;
    assert!((p.average_concentration() - expected).abs() < 1e-9 * c0);
    assert!(p.surface_concentration(0.0) < p.average_concentration());
}

#[test]
fn eis_limits() {
    let c = BatteryEisConfig::default();
    let hi = battery_impedance(1e7, &c).unwrap();
    assert!((hi.re - c.series_resistance_ohm).abs() < 1e-4 && hi.re.is_finite());
    let lo = battery_impedance(1e-6, &c).unwrap();
    assert!(
        lo.re
            > c.series_resistance_ohm
                + c.negative_charge_transfer_resistance_ohm
                + c.positive_charge_transfer_resistance_ohm
    );
    assert!(battery_impedance(0.0, &c).is_err());
}

#[test]
fn dfn_conserves_lithium_and_is_transactional() {
    let mut m = DoyleFullerNewmanModel::new(default_graphite_nmc_dfn_config()).unwrap();
    let li0 = m.total_lithium_mol();
    for i in 0..8 {
        m.step(if i % 2 == 0 { 4.0 } else { -2.0 }, 3.0).unwrap();
    }
    assert!((m.total_lithium_mol() - li0).abs() < 1e-9 * li0);
    let before = m.electrolyte_concentration().to_vec();
    assert!(m.step(f64::NAN, 1.0).is_err());
    assert!(m.step(1.0, 0.0).is_err());
    assert_eq!(m.electrolyte_concentration(), before.as_slice());
    let bad = DfnConfig {
        negative_nodes: 1,
        ..Default::default()
    };
    assert!(DoyleFullerNewmanModel::new(bad).is_err());
}

#[test]
fn discharge_lowers_voltage_and_rest_relaxes() {
    let mut m = SingleParticleElectrolyteModel::new(default_graphite_nmc_config()).unwrap();
    let ocv0 = m.state(0.0).unwrap().open_circuit_voltage_v;
    let load = m.step(3.0, 10.0).unwrap();
    assert!(load.terminal_voltage_v < load.open_circuit_voltage_v);
    assert!(load.open_circuit_voltage_v < ocv0);
    let rest = m.step(0.0, 10.0).unwrap();
    assert!(rest.terminal_voltage_v > load.terminal_voltage_v);
}

#[test]
fn insulated_grid_heats_uniformly_by_energy() {
    let cfg = ThermalGridConfig {
        cells: [3, 2, 2],
        ..Default::default()
    };
    let mut g = ThermalGrid3D::new(cfg, vec![ThermalVoxel::default()]).unwrap();
    let q = vec![0.5; 12];
    let total_c: f64 = g.heat_capacity_j_per_k().iter().sum();
    let r = g.step(4.0, &q).unwrap();
    let dt_expected = 6.0 * 4.0 / total_c;
    for t in g.temperature() {
        assert!((t - 298.15 - dt_expected).abs() < 1e-9);
    }
    assert!(r.energy_balance_error_j.abs() < 1e-9 && r.outward_cooling_power_w.abs() < 1e-12);
    assert!(g.step(1.0, &q[..3]).is_err());
}

#[test]
fn convective_grid_relaxes_to_ambient() {
    let mut b = [ThermalBoundary {
        convection_w_per_m2_k: 50.0,
        ambient_temperature_k: 300.0,
    }; 6];
    b[0].ambient_temperature_k = 300.0;
    let cfg = ThermalGridConfig {
        cells: [2, 2, 2],
        initial_temperature_k: 320.0,
        boundaries: b,
        ..Default::default()
    };
    let mut g = ThermalGrid3D::new(cfg, vec![ThermalVoxel::default()]).unwrap();
    for _ in 0..200 {
        g.step(500.0, &[0.0; 8]).unwrap();
    }
    assert!(g.temperature().iter().all(|t| (t - 300.0).abs() < 1e-3));
}

#[test]
fn pack_cutoff_leaves_state_unchanged() {
    let cfg = default_graphite_nmc_config();
    let pc = PackControl {
        minimum_cell_voltage_v: 3.6,
        ..Default::default()
    };
    let mut pack = SeriesBatteryPack::new(&[cfg.clone(), cfg.clone()], pc).unwrap();
    let t0 = pack.time_s();
    let v0 = pack.cells()[0].state(0.0).unwrap().terminal_voltage_v;
    assert!(pack.step(500.0, 5.0).is_err());
    assert_eq!(pack.time_s(), t0);
    assert_eq!(pack.cells()[0].state(0.0).unwrap().terminal_voltage_v, v0);
    assert!(SeriesBatteryPack::new(&[], PackControl::default()).is_err());
}

#[test]
fn electrothermal_rejects_mismatched_capacity() {
    let cfg = default_graphite_nmc_config();
    let grid = ThermalGrid3D::new(
        ThermalGridConfig {
            cells: [2, 1, 1],
            ..Default::default()
        },
        vec![ThermalVoxel::default()],
    )
    .unwrap();
    assert!(
        ElectrothermalBatteryPack::new(&[cfg], PackControl::default(), grid, vec![vec![0, 1]])
            .is_err()
    );
}
