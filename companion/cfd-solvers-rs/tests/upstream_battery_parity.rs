//! Value-by-value comparison with the C++ upstream battery/{lithium_ion, dfn, pack, thermal3d}.
//! Goldens come from `parity/oracle_battery.cpp`.
use cfd_solvers::battery::*;
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-battery.json")).unwrap()
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
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= rel * c.abs() + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn bsr(s: &BatteryStepResult) -> Vec<f64> {
    vec![
        s.time_s,
        s.current_a,
        s.terminal_voltage_v,
        s.open_circuit_voltage_v,
        s.reaction_overpotential_v,
        s.electrolyte_overpotential_v,
        s.ohmic_drop_v,
        s.negative_surface_stoichiometry,
        s.positive_surface_stoichiometry,
        s.electrolyte_min_concentration_mol_per_m3,
        s.electrolyte_max_concentration_mol_per_m3,
        s.temperature_k,
        s.heat_generation_w,
    ]
}

fn dsr(s: &DfnStepResult) -> Vec<f64> {
    vec![
        s.time_s,
        s.current_a,
        s.terminal_voltage_v,
        s.open_circuit_voltage_v,
        s.reaction_overpotential_v,
        s.ohmic_drop_v,
        s.negative_surface_stoichiometry,
        s.positive_surface_stoichiometry,
        s.electrolyte_min_concentration_mol_per_m3,
        s.electrolyte_max_concentration_mol_per_m3,
        s.temperature_k,
        s.heat_generation_w,
    ]
}

const R: f64 = 1e-10;
const A: f64 = 1e-12;

#[test]
fn particle_ocv_and_kinetics() {
    let mut p = SphericalDiffusionParticle::new(ParticleConfig::default()).unwrap();
    let mut v = vec![];
    for _ in 0..5 {
        p.step(2e-6, 10.0).unwrap();
        v.push(p.average_stoichiometry());
        v.push(p.surface_stoichiometry(2e-6));
    }
    close("particle", &v, R, A);
    close("particle_conc", p.concentrations(), R, A);
    let mut o = vec![];
    for k in 0..9 {
        let x = k as f64 * 0.125;
        o.push(graphite_ocv_v(x));
        o.push(nmc_ocv_v(x));
    }
    close("ocv", &o, R, A);
    close(
        "bv",
        &[
            symmetric_butler_volmer_overpotential_v(3.0, 2.0, 298.15).unwrap(),
            symmetric_butler_volmer_overpotential_v(-0.5, 1.5, 310.0).unwrap(),
        ],
        R,
        A,
    );
}

#[test]
fn spm_spme_drive_cycle_and_eis() {
    let cfg = default_graphite_nmc_config();
    let mut spm = SingleParticleModel::new(cfg.clone()).unwrap();
    let mut v = vec![];
    for i in 0..6 {
        v.extend(bsr(&spm.step(if i < 4 { 2.0 } else { -1.0 }, 5.0).unwrap()));
    }
    let d = spm.degradation_state();
    v.extend([
        d.sei_thickness_m,
        d.plated_lithium_mol,
        d.active_material_fraction,
        d.crack_damage,
    ]);
    close("spm", &v, R, 1e-20);
    let mut spme = SingleParticleElectrolyteModel::new(cfg.clone()).unwrap();
    let mut e = vec![];
    for i in 0..6 {
        e.extend(bsr(&spme
            .step(if i < 4 { 3.0 } else { -1.5 }, 4.0)
            .unwrap()));
    }
    close("spme", &e, R, A);
    close("spme_ce", spme.electrolyte_concentration(), R, A);
    close("spme_phie", spme.electrolyte_potential(), R, 1e-14);
    close("spme_phis", spme.solid_potential(), R, 1e-14);
    let mut spme2 = SingleParticleElectrolyteModel::new(cfg).unwrap();
    let cyc = [
        DriveCyclePoint {
            duration_s: 10.0,
            current_a: 2.0,
        },
        DriveCyclePoint {
            duration_s: 3.5,
            current_a: 0.0,
        },
        DriveCyclePoint {
            duration_s: 6.0,
            current_a: -1.0,
        },
    ];
    let dc = simulate_drive_cycle(&mut spme2, &cyc, 2.0).unwrap();
    let mut dv = vec![
        dc.samples.len() as f64,
        dc.discharged_capacity_ah,
        dc.electrical_energy_wh,
    ];
    dv.extend(bsr(dc.samples.last().unwrap()));
    close("drive", &dv, R, A);
    let z = battery_impedance_spectrum(
        &[0.01, 0.1, 1.0, 10.0, 1000.0],
        &BatteryEisConfig::default(),
    )
    .unwrap();
    close(
        "eis",
        &z.iter().flat_map(|c| [c.re, c.im]).collect::<Vec<_>>(),
        R,
        A,
    );
    let mut slab = ThermalSlab1D::new(ThermalSlabConfig::default()).unwrap();
    for _ in 0..4 {
        slab.step(0.5, &[5e5; 32]).unwrap();
    }
    close("slab", slab.temperature(), R, A);
}

#[test]
fn doyle_fuller_newman() {
    let mut dfn = DoyleFullerNewmanModel::new(default_graphite_nmc_dfn_config()).unwrap();
    let mut v = vec![dfn.total_lithium_mol()];
    for i in 0..5 {
        v.extend(dsr(&dfn.step(if i < 3 { 3.0 } else { -2.0 }, 5.0).unwrap()));
    }
    v.push(dfn.total_lithium_mol());
    close("dfn", &v, R, A);
    close("dfn_ce", dfn.electrolyte_concentration(), R, A);
    close("dfn_phie", dfn.electrolyte_potential(), R, 1e-14);
    close("dfn_phis", dfn.solid_potential(), R, 1e-14);
    close("dfn_j", dfn.local_reaction_current(), R, 1e-12);
    close("dfn_state", &dsr(&dfn.state(1.0).unwrap()), R, A);
}

#[test]
fn balanced_series_pack() {
    let cfg = default_graphite_nmc_config();
    let mut cells = vec![cfg.clone(); 3];
    cells[1].negative.particle.initial_stoichiometry = 0.72;
    cells[1].positive.particle.initial_stoichiometry = 0.45;
    let pc = PackControl {
        balancing_current_a: 0.2,
        balancing_minimum_voltage_v: 3.0,
        ..Default::default()
    };
    let mut pack = SeriesBatteryPack::new(&cells, pc).unwrap();
    let mut v = vec![];
    for _ in 0..4 {
        let r = pack.step(1.5, 5.0).unwrap();
        v.extend([
            r.time_s,
            r.terminal_voltage_v,
            r.terminal_power_w,
            r.balancing_heat_w,
            r.net_terminal_charge_ah,
            r.terminal_energy_wh,
            r.balancing_energy_wh,
        ]);
        v.extend(&r.balancing_currents_a);
        v.extend(r.cells.iter().map(|c| c.terminal_voltage_v));
    }
    close("pack", &v, R, A);
}

#[test]
fn thermal_grid_and_electrothermal_pack() {
    let mut b = [ThermalBoundary::default(); 6];
    b[0] = ThermalBoundary {
        convection_w_per_m2_k: 15.0,
        ambient_temperature_k: 290.0,
    };
    b[3] = ThermalBoundary {
        convection_w_per_m2_k: 40.0,
        ambient_temperature_k: 300.0,
    };
    b[5] = ThermalBoundary {
        convection_w_per_m2_k: 5.0,
        ambient_temperature_k: 280.0,
    };
    let tg = ThermalGridConfig {
        cells: [4, 3, 2],
        spacing_m: [0.01, 0.02, 0.015],
        boundaries: b,
        ..Default::default()
    };
    let mats = (0..24)
        .map(|i| ThermalVoxel {
            volumetric_heat_capacity_j_per_m3_k: 1.5e6 + 1e4 * i as f64,
            conductivity_w_per_m_k: [
                2.0 + 0.1 * i as f64,
                0.5,
                if i % 5 == 0 { 0.0 } else { 1.2 },
            ],
        })
        .collect();
    let mut grid = ThermalGrid3D::new(tg, mats).unwrap();
    let mut q = vec![0.0; 24];
    q[5] = 3.0;
    q[17] = 1.0;
    let mut v = vec![];
    let mut its = vec![];
    for _ in 0..3 {
        let r = grid.step(2.0, &q).unwrap();
        v.extend([
            r.supplied_power_w,
            r.outward_cooling_power_w,
            r.stored_energy_change_j,
            r.energy_balance_error_j,
            r.maximum_temperature_k,
            r.iterations as f64,
        ]);
        its.push(r.iterations);
    }
    // Energy-balance error is round-off sized; iteration counts may differ by one.
    let g = goldens();
    let cpp: Vec<f64> = g["thermal"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    for k in 0..3 {
        for j in [0, 1, 2, 4] {
            let (r, c) = (v[6 * k + j], cpp[6 * k + j]);
            assert!(
                (r - c).abs() <= 1e-9 * c.abs() + 1e-9,
                "thermal[{k},{j}] {r} vs {c}"
            );
        }
        assert!(v[6 * k + 3].abs() < 1e-6 && cpp[6 * k + 3].abs() < 1e-6);
        assert!((v[6 * k + 5] - cpp[6 * k + 5]).abs() <= 1.0);
    }
    close("thermal_T", grid.temperature(), 1e-11, 1e-9);

    let cfg = default_graphite_nmc_config();
    let mut eb = [ThermalBoundary::default(); 6];
    eb[2] = ThermalBoundary {
        convection_w_per_m2_k: 20.0,
        ambient_temperature_k: 295.0,
    };
    let eg = ThermalGridConfig {
        cells: [4, 2, 2],
        spacing_m: [0.01; 3],
        boundaries: eb,
        ..Default::default()
    };
    let cap = cfg.thermal.mass_kg * cfg.thermal.heat_capacity_j_per_kg_k;
    let vox = ThermalVoxel {
        volumetric_heat_capacity_j_per_m3_k: cap / 4.0 / 1e-6,
        conductivity_w_per_m_k: [3.0; 3],
    };
    let grid = ThermalGrid3D::new(eg, vec![vox; 16]).unwrap();
    let mut et = ElectrothermalBatteryPack::new(
        &[cfg.clone(), cfg],
        PackControl::default(),
        grid,
        vec![vec![0, 4, 8, 12], vec![1, 5, 9, 13]],
    )
    .unwrap();
    let mut ev = vec![];
    for _ in 0..3 {
        let r = et.step(4.0, 10.0).unwrap();
        ev.extend([
            r.electrical.terminal_voltage_v,
            r.thermal.supplied_power_w,
            r.thermal.maximum_temperature_k,
            r.thermal.energy_balance_error_j,
            r.electrical.cells[0].temperature_k,
            r.electrical.cells[1].temperature_k,
        ]);
    }
    let cpp: Vec<f64> = g["electrothermal"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    for (i, (r, c)) in ev.iter().zip(&cpp).enumerate() {
        if i % 6 == 3 {
            assert!(r.abs() < 1e-6 && c.abs() < 1e-6);
        } else {
            assert!(
                (r - c).abs() <= 1e-9 * c.abs() + 1e-10,
                "electrothermal[{i}] {r} vs {c}"
            );
        }
    }
    close("electrothermal_T", et.thermal().temperature(), 1e-11, 1e-9);
}
