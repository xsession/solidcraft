//! Value-by-value comparison with the C++ upstream tcad/{semiconductor1d, mos, bipolar_power,
//! advanced_devices}. Goldens come from `parity/oracle_tcad.cpp`.
use cfd_solvers::tcad::*;
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-tcad.json")).unwrap()
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

fn b(x: bool) -> f64 {
    if x {
        1.0
    } else {
        0.0
    }
}

#[test]
fn physics_functions() {
    let rm = RecombinationModel {
        tau_n_s: 1e-6,
        tau_p_s: 2e-6,
        trap_n1_m3: 3e15,
        trap_p1_m3: 0.0,
        auger_n_m6_per_s: 2.8e-43,
        auger_p_m6_per_s: 9.9e-44,
        radiative_m3_per_s: 1.1e-20,
    };
    let f = [
        thermal_voltage(310.0).unwrap(),
        bernoulli_function(3e-6).unwrap(),
        bernoulli_function(0.7).unwrap(),
        bernoulli_function(-60.0).unwrap(),
        bernoulli_function(55.0).unwrap(),
        silicon_bandgap_ev(350.0, 1.17, 4.73e-4, 636.0).unwrap(),
        temperature_scaled_intrinsic_density(1e16, 300.0, 360.0, 1.12, 1.09).unwrap(),
        temperature_scaled_mobility(0.135, 300.0, 400.0, 2.42).unwrap(),
        srh_recombination(1e22, 1e12, 1e16, &rm).unwrap(),
        auger_recombination(1e22, 1e12, 1e16, &rm).unwrap(),
        radiative_recombination(1e22, 1e12, 1e16, &rm).unwrap(),
        total_recombination(1e18, 3e17, 1e16, &rm).unwrap(),
        caughey_thomas_mobility(3e23, 0.0068, 0.1414, 9.2e22, 0.711).unwrap(),
        high_field_mobility(0.135, 2e6, 1e5, 2.0).unwrap(),
        fermi_dirac_half(-9.0).unwrap(),
        fermi_dirac_half(0.3).unwrap(),
        fermi_dirac_half(12.0).unwrap(),
        fermi_dirac_half(45.0).unwrap(),
        inverse_fermi_dirac_half(2.5).unwrap(),
        fermi_electron_density(2.8e25, -1.0).unwrap(),
        fermi_hole_density(1.04e25, 2.0).unwrap(),
        ionized_donor_density(1e23, -3.0, 0.045, 2.0, 300.0).unwrap(),
        ionized_acceptor_density(1e23, -2.0, 0.045, 4.0, 250.0).unwrap(),
    ];
    close("funcs", &f, 1e-13, 0.0);
}

#[test]
fn drift_diffusion_device() {
    let pc = Semiconductor1DConfig {
        length_m: 2e-6,
        nodes: 101,
        max_gummel_iterations: 1000,
        relative_tolerance: 2e-7,
        under_relaxation: 0.05,
        ..Default::default()
    };
    let mut pn = SemiconductorDevice1D::new(pc).unwrap();
    pn.set_net_doping(&pn_junction_doping(pn.x(), 1e-6, 1e21, 1e21).unwrap())
        .unwrap();
    let eq = pn.solve_equilibrium().unwrap();
    assert!(eq.converged);
    close(
        "pn_eq",
        &[
            b(eq.converged),
            eq.iterations as f64,
            eq.relative_change,
            eq.terminal_current_a,
            pn.total_mobile_charge_c(),
        ],
        1e-6,
        1e-17,
    );
    close("pn_eq_phi", pn.potential(), 1e-9, 1e-12);
    close("pn_n", pn.electron_density(), 1e-8, 1e-3);
    close("pn_p", pn.hole_density(), 1e-8, 1e-3);

    let mut rc = Semiconductor1DConfig {
        nodes: 41,
        area_m2: 2e-12,
        max_gummel_iterations: 600,
        relative_tolerance: 1e-9,
        under_relaxation: 0.4,
        ..Default::default()
    };
    rc.right.voltage_v = -0.01;
    rc.recombination.radiative_m3_per_s = 1e-18;
    let mut bar = SemiconductorDevice1D::new(rc).unwrap();
    bar.set_net_doping_fn(|x| 1e21 + 0.0 * x).unwrap();
    bar.initialize_charge_neutral().unwrap();
    let dc = bar.solve_dc().unwrap();
    let mut v = vec![
        b(dc.converged),
        dc.iterations as f64,
        dc.terminal_current_a,
        bar.left_terminal_charge_c(),
    ];
    for p in bar.sweep_right_contact(-0.005, -0.02, 3).unwrap() {
        v.extend([p.voltage_v, p.current_a, b(p.converged)]);
    }
    close("bar", &v, 1e-7, 1e-24);
    close("bar_phi", bar.potential(), 1e-10, 1e-14);
    close("bar_J", &bar.total_current_density_faces(), 1e-7, 1e-3);
    close("bar_heat", &bar.joule_heating_cells(), 1e-7, 1e-1);
    close("bar_E", &bar.electric_field_faces(), 1e-7, 1e-3);
    let ss = bar.small_signal(1e6, 1e-5).unwrap();
    close(
        "bar_ac",
        &[
            ss.bias_voltage_v,
            ss.conductance_s,
            ss.capacitance_f,
            ss.admittance_s.re,
            ss.admittance_s.im,
            ss.impedance_ohm.re,
            ss.impedance_ohm.im,
            b(ss.converged),
        ],
        1e-5,
        1e-20,
    );
    let mut cv = vec![];
    for c in bar.sweep_cv(-1e-3, 1e-3, 3, 1e-5).unwrap() {
        cv.extend([
            c.voltage_v,
            c.capacitance_f,
            c.electrode_charge_c,
            b(c.converged),
        ]);
    }
    close("bar_cv", &cv, 1e-5, 1e-21);
    bar.set_contact_voltages(0.0, -0.012).unwrap();
    let mut tr = vec![];
    for _ in 0..3 {
        let t = bar.step_transient(1e-12).unwrap();
        tr.extend([
            b(t.converged),
            t.iterations as f64,
            t.time_s,
            t.conduction_current_a,
            t.displacement_current_a,
            t.terminal_current_a,
        ]);
    }
    close("bar_tr", &tr, 1e-6, 1e-16);

    let mut sc = Semiconductor1DConfig {
        nodes: 11,
        ..Default::default()
    };
    sc.left.contact_type = ContactType::Schottky;
    sc.left.electron_barrier_ev = 0.72;
    let mut sd = SemiconductorDevice1D::new(sc).unwrap();
    sd.set_net_doping_fn(|_| 1e21).unwrap();
    sd.set_temperature(340.0).unwrap();
    sd.initialize_charge_neutral().unwrap();
    let m = sd.config().material;
    close(
        "schottky",
        &[
            m.intrinsic_density_m3,
            m.electron_mobility_m2_per_vs,
            m.hole_mobility_m2_per_vs,
            m.conduction_band_density_m3,
            sd.electron_density()[0],
            sd.hole_density()[0],
        ],
        1e-12,
        0.0,
    );
    close("schottky_phi", sd.potential(), 1e-12, 1e-15);

    let mut gc = Semiconductor1DConfig {
        nodes: 31,
        length_m: 0.5e-6,
        max_gummel_iterations: 2000,
        under_relaxation: 0.1,
        relative_tolerance: 1e-7,
        ..Default::default()
    };
    gc.left = Contact1D {
        contact_type: ContactType::Gate,
        oxide_capacitance_f_per_m2: 3.45e-3,
        voltage_v: 0.2,
        ..Default::default()
    };
    let mut gd = SemiconductorDevice1D::new(gc).unwrap();
    gd.set_net_doping(&[-1e21; 31]).unwrap();
    gd.initialize_charge_neutral().unwrap();
    let gdc = gd.solve_dc().unwrap();
    close(
        "gate",
        &[
            b(gdc.converged),
            gdc.iterations as f64,
            gd.left_terminal_charge_c(),
            gd.potential()[0],
        ],
        1e-6,
        0.0,
    );
    close(
        "pin",
        &pin_diode_doping(gd.x(), 1e-7, 3e-7, 1e23, 2e23, 1e19).unwrap(),
        0.0,
        0.0,
    );
}

#[test]
fn mos_bjt_igbt_and_compact_models() {
    let mc = MosCapacitor1D::new(MosCapacitorConfig::default()).unwrap();
    let mut mv = vec![mc.oxide_capacitance_per_area()];
    for p in mc.sweep_cv(-2.0, 2.0, 7, 1e-4).unwrap() {
        mv.extend([p.surface_potential_v, p.gate_charge_c, p.capacitance_f]);
    }
    close("moscap", &mv, 1e-9, 1e-25);
    let lc = LongChannelMosfetConfig {
        body_effect_coefficient_sqrt_v: 0.4,
        ..Default::default()
    };
    let mos = LongChannelMosfet::new(lc).unwrap();
    let mut mo = vec![
        mos.threshold_voltage(-1.0),
        mos.drain_current(1.5, 0.3, 0.0).unwrap(),
        mos.drain_current(1.5, 2.0, -0.5).unwrap(),
        mos.drain_current(1.2, -0.4, 0.0).unwrap(),
        mos.drain_current(0.2, 1.0, 0.0).unwrap(),
    ];
    let op = mos.operating_point(1.4, 0.6, -0.2).unwrap();
    mo.extend([
        op.drain_current_a,
        op.transconductance_s,
        op.output_conductance_s,
        op.threshold_voltage_v,
    ]);
    mo.extend(
        mos.sweep(0.0, 2.0, 3, 0.0, 1.5, 3, 0.0)
            .unwrap()
            .iter()
            .map(|s| s.drain_current_a),
    );
    let me = make_long_channel_mosfet_spice_evaluator(lc, 1e-6).unwrap()(&[1.0, 1.6, 0.1, -0.2])
        .unwrap();
    mo.extend(&me.terminal_current);
    mo.extend(&me.jacobian);
    close("mosfet", &mo, 1e-9, 1e-15);

    let bjt = BipolarJunctionTransistor1D::new(BipolarJunctionConfig::default()).unwrap();
    let (d, d2) = (bjt.derived(0.0).unwrap(), bjt.derived(350.0).unwrap());
    let mut bv = vec![
        d.thermal_voltage_v,
        d.electron_diffusion_length_base_m,
        d.hole_diffusion_length_emitter_m,
        d.base_transport_factor,
        d.emitter_injection_efficiency,
        d.collector_injection_efficiency,
        d.forward_alpha,
        d.reverse_alpha,
        d.forward_beta,
        d.reverse_beta,
        d.emitter_saturation_current_a,
        d.collector_saturation_current_a,
        d.transport_saturation_current_a,
        d.emitter_base_built_in_v,
        d.collector_base_built_in_v,
        d2.forward_beta,
        d2.emitter_saturation_current_a,
    ];
    let bo = bjt.operating_point(3.0, 0.7, 0.0, 0.0).unwrap();
    bv.extend([
        bo.collector_current_a,
        bo.base_current_a,
        bo.emitter_current_a,
        bo.dissipated_power_w,
    ]);
    bv.extend(bjt.doping_profile(7).unwrap());
    let pnp = BipolarJunctionConfig {
        polarity: BipolarPolarity::Pnp,
        ..Default::default()
    };
    let po = BipolarJunctionTransistor1D::new(pnp)
        .unwrap()
        .operating_point(-2.0, -0.68, 0.0, 320.0)
        .unwrap();
    bv.extend([
        po.collector_current_a,
        po.base_current_a,
        po.emitter_current_a,
    ]);
    let be = make_bipolar_tcad_spice_evaluator(BipolarJunctionConfig::default(), 1e-6).unwrap()(&[
        2.0, 0.65, 0.0,
    ])
    .unwrap();
    bv.extend(&be.terminal_current);
    bv.extend(&be.jacobian);
    close("bjt", &bv, 1e-9, 1e-30);

    let mut ig = IgbtPowerDevice::new(IgbtConfig::default()).unwrap();
    let io = ig.operating_point(5.0, 3.0, 0.0).unwrap();
    let mut iv = vec![
        io.collector_current_a,
        io.channel_current_a,
        io.bipolar_current_a,
        io.effective_drift_resistance_ohm,
        io.dissipated_power_w,
    ];
    let io2 = ig.operating_point(5.0, 3.0, 360.0).unwrap();
    iv.extend([io2.collector_current_a, io2.effective_drift_resistance_ohm]);
    for i in 0..4 {
        let s = ig.step(if i < 2 { 6.0 } else { 0.0 }, 2.5, 1e-6).unwrap();
        iv.extend([
            s.collector_current_a,
            s.bipolar_current_a,
            s.junction_temperature_k,
            s.dissipated_power_w,
        ]);
    }
    let ie = make_igbt_tcad_spice_evaluator(IgbtConfig::default(), 1e-6).unwrap()(&[3.0, 5.0, 0.0])
        .unwrap();
    iv.extend(&ie.terminal_current);
    iv.extend(&ie.jacobian);
    close("igbt", &iv, 1e-9, 1e-12);

    let gp = GummelPoonBjt::new(GummelPoonConfig::default()).unwrap();
    let l23 = MosLevel23::new(MosLevel23Config::default()).unwrap();
    let vd = VdmosPowerDevice::new(VdmosConfig::default()).unwrap();
    let mut core = HystereticCoreModel::new(HystereticCoreConfig::default()).unwrap();
    let mut av = vec![];
    for (x, y) in [(0.7, -2.0), (0.85, 0.1), (-0.2, 0.6)] {
        let g = gp.operating_point(x, y).unwrap();
        av.extend([
            g.collector_current_a,
            g.base_current_a,
            g.emitter_current_a,
            g.forward_transport_current_a,
            g.reverse_transport_current_a,
        ]);
    }
    for (x, y) in [(0.3, 1.0), (1.5, 0.4), (1.5, 3.0)] {
        let m = l23.operating_point(x, y).unwrap();
        av.extend([
            m.drain_current_a,
            m.transconductance_s,
            m.output_conductance_s,
            b(m.subthreshold),
        ]);
    }
    for (x, y) in [(2.0, 1.0), (2.0, -0.7)] {
        let o = vd.operating_point(x, y).unwrap();
        av.extend([
            o.drain_current_a,
            o.channel_current_a,
            o.body_diode_current_a,
            o.drift_drop_v,
        ]);
    }
    for i in 0..12 {
        core.step(120.0 * (0.6 * i as f64).sin(), 300.0 + 5.0 * i as f64, 3e-5)
            .unwrap();
        av.extend([core.state().induction_t, core.state().loss_density_j_per_m3]);
    }
    close("advanced", &av, 1e-12, 1e-30);
}
