//! End-to-end checks of the 0.7.0 domain case adapters against closed forms, conservation
//! laws and textbook values. Each test runs a (sometimes reduced) shipped example.
use cfd_solvers::case::{Case, RunManifest};
use serde_json::Value;
use std::path::PathBuf;

fn example(name: &str) -> Value {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join(format!("{name}.json"));
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn run(name: &str, edit: impl FnOnce(&mut Value)) -> (RunManifest, PathBuf) {
    let mut v = example(name);
    edit(&mut v);
    let case: Case = serde_json::from_value(v).unwrap();
    let dir = std::env::temp_dir().join(format!("cfd_domain_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let m = case.run(&dir).unwrap();
    assert!(m.execution_completed && m.capabilities.iter().any(|c| c == "domain_cases_v1"));
    for f in &m.files {
        assert!(dir.join(f).exists(), "{name}: missing {f}");
    }
    (m, dir)
}

fn metric(m: &RunManifest, k: &str) -> f64 {
    *m.metrics
        .get(k)
        .unwrap_or_else(|| panic!("missing metric {k}: {:?}", m.metrics))
}

fn csv_column(dir: &std::path::Path, file: &str, column: &str) -> Vec<f64> {
    let text = std::fs::read_to_string(dir.join(file)).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines
        .next()
        .unwrap()
        .split(',')
        .map(|h| h.trim_matches('"'))
        .collect();
    let c = header
        .iter()
        .position(|h| *h == column)
        .unwrap_or_else(|| panic!("{column} not in {header:?}"));
    lines
        .map(|l| l.split(',').nth(c).unwrap().parse().unwrap())
        .collect()
}

#[test]
fn spice_rc_filter_has_its_corner_and_settles() {
    let (_, dir) = run("spice-rc-filter", |_| {});
    // Below the 1.59 kHz corner (and below the clamp) the output follows the input.
    let f = csv_column(&dir, "ac.csv", "frequency_Hz");
    let g = csv_column(&dir, "ac.csv", "out_dB");
    let fc = 1.0 / (2.0 * std::f64::consts::PI * 1e3 * 100e-9);
    for (fi, gi) in f.iter().zip(&g) {
        let exact = -10.0 * (1.0 + (fi / fc).powi(2)).log10();
        assert!((gi - exact).abs() < 0.01, "{fi} Hz: {gi} vs {exact}");
    }
    // The transient starts at zero and the 3 V + diode (N = 1.8) clamp holds the output near 4.2 V.
    let v = csv_column(&dir, "transient.csv", "V_out");
    assert!(v[0].abs() < 1e-9 && v.iter().all(|x| *x < 4.3) && v.iter().any(|x| *x > 4.0));
}

#[test]
fn tcad_diode_is_exponential_in_forward_bias() {
    let (m, dir) = run("tcad-pn-diode", |v| v["cv"] = Value::Null);
    let vt = 0.025852;
    let i = csv_column(&dir, "iv.csv", "anode_current_A");
    let v = csv_column(&dir, "iv.csv", "anode_voltage_V");
    let n = v.len();
    let ideality = (v[n - 1] - v[n - 3]) / (vt * (i[n - 1] / i[n - 3]).ln());
    assert!(ideality > 0.95 && ideality < 1.3, "ideality {ideality}");
    // Built-in potential 2 Vt ln(N/ni) with ni = 1e16 m^-3.
    assert!((metric(&m, "built_in_potential_V") - 2.0 * vt * (1e5f64).ln()).abs() < 0.01);
}

#[test]
fn battery_discharge_stops_at_cutoff_and_pack_scales_voltage() {
    let (m, _) = run("battery-discharge", |_| {});
    assert_eq!(metric(&m, "stopped_on_cutoff"), 1.0);
    let ah = metric(&m, "net_discharge_Ah");
    assert!(ah > 0.2 && ah < 0.35, "{ah}");
    let (p, _) = run("battery-pack", |v| {
        v["protocol"] = serde_json::json!([{"current_a": 0.0, "duration_s": 10.0}])
    });
    let (c, _) = run("battery-discharge", |v| {
        v["protocol"] = serde_json::json!([{"current_a": 0.0, "duration_s": 10.0}])
    });
    assert!((metric(&p, "final_voltage_V") - 4.0 * metric(&c, "final_voltage_V")).abs() < 1e-9);
}

#[test]
fn half_wave_dipole_resonates_near_73_ohm() {
    let (m, _) = run("wire-dipole", |_| {});
    let r = metric(&m, "resonant_resistance_ohm");
    let f = metric(&m, "resonant_frequency_Hz");
    let electrical = 0.48 * f / 299_792_458.0;
    assert!(r > 62.0 && r < 80.0, "{r}");
    assert!(electrical > 0.45 && electrical < 0.49, "{electrical}");
    // A monopole over a PEC plane has half the dipole impedance.
    let (g, _) = run("wire-dipole", |v| {
        v["ground"] = "pec_plane".into();
        v["length_m"] = 0.24.into();
    });
    assert!((metric(&g, "resonant_resistance_ohm") - 0.5 * r).abs() < 1e-6 * r);
}

#[test]
fn fdtd_cavity_conserves_energy_after_the_source() {
    let (_, dir) = run("fdtd3d-cavity", |v| v["steps"] = 200.into());
    let t = csv_column(&dir, "energy.csv", "time_s");
    let e = csv_column(&dir, "energy.csv", "field_energy_J");
    let after: Vec<f64> = t
        .iter()
        .zip(&e)
        .filter(|(t, _)| **t > 1.6e-10)
        .map(|(_, e)| *e)
        .collect();
    let (lo, hi) = after
        .iter()
        .fold((f64::MAX, f64::MIN), |a, e| (a.0.min(*e), a.1.max(*e)));
    assert!(hi > 0.0 && (hi - lo) / hi < 0.02, "{lo} {hi}");
}

#[test]
fn edge_cavity_impedance_peaks_inside_the_sweep() {
    let (m, _) = run("edge-cavity", |_| {});
    assert!(m.linear_solver_converged == Some(true));
    assert!(metric(&m, "peak_impedance_ohm") > 0.0 && metric(&m, "tetrahedra") == 144.0);
    // TE101 of a 100 x 50 mm PEC box: (c/2)·sqrt(1/a² + 1/d²) = 3.352 GHz.
    let te101 = 0.5 * 299_792_458.0 * (1.0f64 / 0.01 + 1.0 / 0.0025).sqrt();
    let first = metric(&m, "first_resonance_frequency_Hz");
    assert!((first / te101 - 1.0).abs() < 0.08, "{first} vs {te101}");
}

#[test]
fn pdn_drop_is_largest_at_the_hotspot_and_sources_hold() {
    let (m, _) = run("pdn-ir-drop", |_| {});
    assert!(m.linear_solver_converged == Some(true));
    let drop = metric(&m, "maximum_ir_drop_V");
    assert!(drop > 0.0 && (metric(&m, "minimum_voltage_V") - (1.0 - drop)).abs() < 1e-12);
}

#[test]
fn galvanic_couple_balances_anodic_and_cathodic_current() {
    let (m, _) = run("galvanic-corrosion", |_| {});
    let e = metric(&m, "corrosion_potential_V");
    // Equal Tafel slopes: E = (ln(i_c A_c / i_a A_a)·b + E_a + E_c) / 2 with b = RT/(αF).
    let b = 8.314462618 * 298.15 / (0.5 * 96485.33212);
    let expected = 0.5 * (b * (1e-9 * 1e-2 / (1e-3 * 1e-3_f64)).ln() + (-0.76 + 0.4));
    assert!((e - expected).abs() < 2e-3, "{e} vs {expected}");
}

#[test]
fn cantilever_plate_and_bar_modes_match_beam_theory() {
    let (m, _) = run("elasticity-cantilever", |_| {});
    // Euler-Bernoulli tip deflection qL^4/(8EI) = 5.78e-5 m; linear triangles are stiffer.
    let u = metric(&m, "max_displacement_m");
    assert!(u > 0.6 * 5.78e-5 && u < 1.05 * 5.78e-5, "{u}");
    let (b, _) = run("bar-modes", |_| {});
    assert!((metric(&b, "mode_1_Hz") / 1261.8 - 1.0).abs() < 1e-3);
}

#[test]
fn j2_yields_at_the_yield_stress_and_hardens() {
    let (m, dir) = run("j2-cyclic", |_| {});
    assert!(metric(&m, "max_lateral_stress_residual_Pa") < 1.0);
    let s = csv_column(&dir, "stress_strain.csv", "axial_stress_Pa");
    let e = csv_column(&dir, "stress_strain.csv", "axial_strain");
    // Elastic until 250 MPa / 200 GPa = 1.25e-3.
    let k = e.iter().position(|x| *x >= 1e-3).unwrap();
    assert!((s[k] - 2e11 * e[k]).abs() < 1e-6 * s[k]);
    assert!(metric(&m, "max_axial_stress_Pa") > 2.5e8);
}

#[test]
fn phase_field_bar_softens_after_its_peak() {
    let (m, dir) = run("phase-field-bar", |_| {});
    let f = csv_column(&dir, "load_displacement.csv", "reaction_force_N");
    assert!(
        metric(&m, "peak_force_N") > *f.last().unwrap()
            && metric(&m, "final_maximum_damage") > 0.99
    );
}

#[test]
fn cosserat_rod_and_dem_stay_finite_and_settle() {
    let (c, _) = run("cosserat-cantilever", |v| v["steps"] = 1500.into());
    assert!(metric(&c, "tip_y_m") < 0.0 && metric(&c, "tip_deflection_m") < 0.2);
    let (d, _) = run("dem-settling", |v| {
        v["count"] = serde_json::json!([3, 3, 3]);
        v["steps"] = 6000.into();
        v["frames"] = 4.into();
    });
    assert_eq!(metric(&d, "escaped_particles"), 0.0);
    // Settled pile: kinetic energy small, centre of mass low.
    assert!(
        metric(&d, "final_kinetic_energy_J") < 2e-3 && metric(&d, "final_mean_height_m") < 0.02
    );
}

#[test]
fn two_stream_grows_near_the_cold_beam_rate() {
    let (m, _) = run("pic-two-stream", |_| {});
    let (g, t) = (
        metric(&m, "measured_growth_rate_per_s"),
        metric(&m, "theory_growth_rate_per_s"),
    );
    assert!((g / t - 1.0).abs() < 0.15, "{g} vs {t}");
    assert!(metric(&m, "relative_energy_change").abs() < 0.02);
}

#[test]
fn transport_deposits_most_energy_inside_the_slab() {
    let (m, _) = run("transport-depth-dose", |v| v["histories"] = 300.into());
    let deposited = metric(&m, "deposited_fraction");
    let secondaries = metric(&m, "secondary_energy_fraction");
    // The remainder is the energy left on primaries stopped at the 10 keV cut.
    assert!(
        (deposited + secondaries - 1.0).abs() < 1e-3,
        "{deposited} + {secondaries}"
    );
    assert!(metric(&m, "mean_primary_range_m") < 0.07);
}

#[test]
fn bioheat_joule_and_thmc_couplings_behave() {
    let (b, _) = run("bioheat-hotspot", |v| v["steps"] = 60.into());
    assert!(metric(&b, "max_temperature_K") > 310.15);
    let (j, _) = run("joule-busbar", |v| v["steps"] = 20.into());
    assert!(
        (metric(&j, "joule_power_W_per_m_depth") / metric(&j, "analytic_power_W_per_m_depth")
            - 1.0)
            .abs()
            < 1e-8
    );
    let (t, _) = run("thmc-column", |v| v["steps"] = 200.into());
    assert!(metric(&t, "max_temperature") > 300.0);
}

#[test]
fn microstructure_and_high_order_methods_conserve() {
    let (c, _) = run("cahn-hilliard", |v| v["steps"] = 4000.into());
    assert!(
        (metric(&c, "mean_composition") - 0.5).abs() < 1e-12
            && metric(&c, "free_energy_drop") >= 0.0
    );
    let (d, _) = run("dg-advection", |_| {});
    assert!(metric(&d, "mass_drift").abs() < 1e-12 && metric(&d, "max_error") < 5e-3);
    let unstable = run_err("cahn-hilliard", |v| v["dt"] = 1e-5.into());
    assert!(unstable.contains("unstable"), "{unstable}");
}

fn run_err(name: &str, edit: impl FnOnce(&mut Value)) -> String {
    let mut v = example(name);
    edit(&mut v);
    let case: Case = serde_json::from_value(v).unwrap();
    let dir = std::env::temp_dir().join(format!("cfd_domain_err_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    case.run(&dir).unwrap_err().to_string()
}

#[test]
fn landau_md_and_photonic_match_references() {
    let (l, _) = run("landau-damping", |_| {});
    let (g, a) = (
        metric(&l, "measured_damping_rate"),
        metric(&l, "asymptotic_damping_rate"),
    );
    assert!(g < 0.0 && (g / a - 1.0).abs() < 0.35, "{g} vs {a}");
    let (md, _) = run("lennard-jones-md", |v| v["steps"] = 3000.into());
    assert!(
        metric(&md, "nve_energy_drift_relative") < 0.01
            && (metric(&md, "mean_temperature") - 0.3).abs() < 0.1
    );
    let (p, _) = run("photonic-waveguide", |_| {});
    let n = metric(&p, "fundamental_effective_index");
    assert!(n > 2.085f64.sqrt() && n < 12.11f64.sqrt());
}

#[test]
fn hostile_domain_inputs_are_errors() {
    for (name, key, value) in [
        ("spice-rc-filter", "netlist", Value::String("R1 a\n".into())),
        ("wire-dipole", "segments", 4.into()),
        ("tcad-pn-diode", "nodes", 3.into()),
        ("pdn-ir-drop", "nx", 1.into()),
        ("dem-settling", "container_width_m", 0.01.into()),
        ("photonic-waveguide", "nx", 40.into()),
        ("lennard-jones-md", "box_length", 10.into()),
    ] {
        let e = run_err(name, |v| v[key] = value);
        assert!(!e.is_empty());
    }
}

#[test]
fn power_network_solves_traces_the_pv_curve_and_settles() {
    let (m, dir) = run("power-network", |_| {});
    assert!(m.linear_solver_converged == Some(true) && metric(&m, "base_residual_rms") < 1e-6);
    // Voltages fall monotonically with load until the nose of the PV curve.
    let v = csv_column(&dir, "pv_curve.csv", "V_elm_pu");
    assert!(v.windows(2).all(|w| w[1] < w[0]) && metric(&m, "last_converged_load_factor") > 2.0);
    assert!((metric(&m, "dispatch_served") - 1.85).abs() < 1e-12);
    // Swing settles: speeds vanish at the synchronous equilibrium.
    assert!(metric(&m, "swing_max_final_speed") < 1e-4);
}

#[test]
fn multiphase_and_spray_conserve_their_inventories() {
    let (b, _) = run("multiphase-bubble-column", |v| v["steps"] = 600.into());
    assert!(metric(&b, "volume_fraction_drift_relative").abs() < 1e-12);
    assert!(
        metric(&b, "final_centroid_y_m") > metric(&b, "initial_centroid_y_m"),
        "bubbles rise"
    );
    let (v, _) = run("multiphase-vof-stir", |v| v["steps"] = 300.into());
    assert!(
        metric(&v, "mass_drift_relative").abs() < 1e-12
            && metric(&v, "liquid_volume_drift_relative").abs() < 1e-12
    );
    let (s, _) = run("spray-evaporation", |v| v["steps"] = 600.into());
    assert!(metric(&s, "mass_balance_residual_fraction").abs() < 1e-10);
    let e = metric(&s, "evaporated_fraction");
    assert!(e > 0.0 && e < 1.0, "{e}");
}
