//! Analytic checks for the TCAD module.
use cfd_solvers::tcad::*;

fn near(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs())
}

#[test]
fn statistics_and_bernoulli() {
    let x = 0.73;
    assert!(near(
        bernoulli_function(-x).unwrap(),
        x.exp() * bernoulli_function(x).unwrap(),
        1e-12
    ));
    assert!(near(bernoulli_function(0.0).unwrap(), 1.0, 1e-15));
    assert!(near(
        fermi_dirac_half(-10.0).unwrap(),
        (-10.0f64).exp(),
        1e-12
    ));
    assert!(near(fermi_dirac_half(-7.0).unwrap(), (-7.0f64).exp(), 1e-3));
    let eta = inverse_fermi_dirac_half(fermi_dirac_half(3.2).unwrap()).unwrap();
    assert!((eta - 3.2).abs() < 1e-9);
    let rm = RecombinationModel {
        tau_n_s: 1e-6,
        tau_p_s: 2e-6,
        auger_n_m6_per_s: 2e-43,
        auger_p_m6_per_s: 1e-43,
        radiative_m3_per_s: 1e-16,
        ..Default::default()
    };
    assert!(total_recombination(1e16, 1e16, 1e16, &rm).unwrap().abs() < 1e-12);
    assert!(total_recombination(5e20, 2e20, 1e16, &rm).unwrap() > 0.0);
    assert!(bernoulli_function(f64::NAN).is_err());
}

#[test]
fn pn_built_in_voltage_and_ohmic_bar() {
    let ni = 1e16;
    let cfg = Semiconductor1DConfig {
        length_m: 2e-6,
        nodes: 101,
        max_gummel_iterations: 1000,
        relative_tolerance: 2e-7,
        under_relaxation: 0.05,
        ..Default::default()
    };
    let mut pn = SemiconductorDevice1D::new(cfg).unwrap();
    pn.set_net_doping(&pn_junction_doping(pn.x(), 1e-6, 1e21, 1e21).unwrap())
        .unwrap();
    let eq = pn.solve_equilibrium().unwrap();
    assert!(eq.converged && eq.terminal_current_a.abs() < 1e-12);
    let vt = thermal_voltage(300.0).unwrap();
    let vbi = pn.potential()[100] - pn.potential()[0];
    assert!(near(vbi, 2.0 * vt * (1e21f64 / ni).ln(), 2e-5));

    let mut rc = Semiconductor1DConfig {
        nodes: 51,
        max_gummel_iterations: 300,
        relative_tolerance: 1e-8,
        under_relaxation: 0.25,
        ..Default::default()
    };
    rc.right.voltage_v = -0.01;
    let mut bar = SemiconductorDevice1D::new(rc).unwrap();
    bar.set_net_doping_fn(|_| 1e21).unwrap();
    bar.initialize_charge_neutral().unwrap();
    let dc = bar.solve_dc().unwrap();
    assert!(dc.converged);
    let expected_j =
        ELEMENTARY_CHARGE_C * rc.material.electron_mobility_m2_per_vs * 1e21 * 0.01 / rc.length_m;
    let j = bar.total_current_density_faces();
    assert!(near(j[0], expected_j, 2e-4));
    let heat = bar.joule_heating_cells();
    let dx = rc.length_m / 50.0;
    let power: f64 = heat
        .iter()
        .enumerate()
        .map(|(i, h)| if i == 0 || i == 50 { 0.5 * h } else { *h })
        .sum::<f64>()
        * dx
        * rc.area_m2;
    assert!(near(power, (dc.terminal_current_a * 0.01).abs(), 3e-4));
    // Two insulating contacts have no potential reference.
    let mut bad = rc;
    bad.left.contact_type = ContactType::Insulating;
    bad.right.contact_type = ContactType::Insulating;
    let mut d = SemiconductorDevice1D::new(bad).unwrap();
    d.initialize_charge_neutral().unwrap();
    let before = d.clone();
    assert!(d.solve_dc().is_err());
    assert_eq!(d, before);
}

#[test]
fn schottky_boundary_and_mos_capacitor() {
    let mut sc = Semiconductor1DConfig {
        nodes: 11,
        ..Default::default()
    };
    sc.left.contact_type = ContactType::Schottky;
    sc.left.electron_barrier_ev = 0.72;
    let mut s = SemiconductorDevice1D::new(sc).unwrap();
    s.set_net_doping_fn(|_| 1e21).unwrap();
    s.initialize_charge_neutral().unwrap();
    let expected =
        sc.material.conduction_band_density_m3 * (-0.72 / thermal_voltage(300.0).unwrap()).exp();
    assert!(near(s.electron_density()[0], expected, 1e-12));

    let mc = MosCapacitor1D::new(MosCapacitorConfig::default()).unwrap();
    let cox = mc.oxide_capacitance_per_area() * 1e-10;
    let acc = mc.solve(-3.0, 1e-4).unwrap();
    let dep = mc.solve(0.3, 1e-4).unwrap();
    assert!(acc.capacitance_f > 0.9 * cox && acc.capacitance_f <= cox * 1.0001);
    assert!(dep.capacitance_f < 0.7 * cox && dep.surface_potential_v > 0.0);
}

#[test]
fn compact_transistor_identities() {
    let m = LongChannelMosfet::new(LongChannelMosfetConfig {
        channel_length_modulation_per_v: 0.0,
        ..Default::default()
    })
    .unwrap();
    let c = m.config();
    let k = c.electron_mobility_m2_per_vs * c.oxide_capacitance_f_per_m2 * c.width_m / c.length_m;
    assert!(near(
        m.drain_current(1.5, 3.0, 0.0).unwrap(),
        0.5 * k,
        1e-12
    ));
    assert_eq!(m.drain_current(0.3, 1.0, 0.0).unwrap(), 0.0);
    let ev =
        make_long_channel_mosfet_spice_evaluator(*c, 1e-6).unwrap()(&[1.0, 1.5, 0.0, 0.0]).unwrap();
    assert!(near(ev.terminal_current[0], -ev.terminal_current[2], 1e-15));

    let bjt = BipolarJunctionTransistor1D::new(BipolarJunctionConfig::default()).unwrap();
    let p = bjt.derived(0.0).unwrap();
    assert!(near(
        p.forward_alpha * p.emitter_saturation_current_a,
        p.reverse_alpha * p.collector_saturation_current_a,
        1e-9
    ));
    let op = bjt.operating_point(2.0, 0.7, 0.0, 0.0).unwrap();
    assert!(
        (op.collector_current_a + op.base_current_a + op.emitter_current_a).abs()
            < 1e-15 * op.collector_current_a.abs().max(1e-30)
    );
    assert!(near(
        op.collector_current_a / op.base_current_a,
        p.forward_beta,
        1e-6
    ));

    let igbt = IgbtPowerDevice::new(IgbtConfig::default()).unwrap();
    assert_eq!(
        igbt.operating_point(10.0, 0.5, 0.0)
            .unwrap()
            .collector_current_a,
        0.0
    );
    assert!(
        igbt.operating_point(10.0, 3.0, 0.0)
            .unwrap()
            .collector_current_a
            > 0.0
    );

    let gp = GummelPoonBjt::new(GummelPoonConfig::default())
        .unwrap()
        .operating_point(0.6, -1.0)
        .unwrap();
    assert!((gp.collector_current_a + gp.base_current_a + gp.emitter_current_a).abs() < 1e-24);
}

#[test]
fn hysteresis_loop_dissipates() {
    let mut core = HystereticCoreModel::new(HystereticCoreConfig::default()).unwrap();
    for i in 0..400 {
        core.step(
            200.0 * (2.0 * std::f64::consts::PI * i as f64 / 100.0).sin(),
            300.0,
            2e-5,
        )
        .unwrap();
    }
    let s = *core.state();
    assert!(s.loss_density_j_per_m3 > 0.0 && s.induction_t.abs() < 1.5);
    assert!(core.step(1.0, 300.0, 0.0).is_err());
}
