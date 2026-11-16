// Test oracle only. Exercises upstream tcad/{semiconductor1d, mos, bipolar_power, advanced_devices}.
#include "cfd/tcad/advanced_devices.hpp"
#include "cfd/tcad/bipolar_power.hpp"
#include "cfd/tcad/mos.hpp"
#include "cfd/tcad/semiconductor1d.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace cfd::tcad;
static bool first = true;
static void vec(const std::string& k, const std::vector<double>& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) {
        double d = v[i];
        std::cout << (i ? "," : "");
        if (std::isnan(d)) std::cout << "\"nan\""; else if (std::isinf(d)) std::cout << (d > 0 ? "\"inf\"" : "\"-inf\""); else std::cout << std::setprecision(17) << d;
    }
    std::cout << "]";
}
int main() {
    std::cout << "{\n";
    RecombinationModel rm{1e-6, 2e-6, 3e15, 0.0, 2.8e-43, 9.9e-44, 1.1e-20};
    std::vector<double> f{thermal_voltage(310.0), bernoulli_function(3e-6), bernoulli_function(0.7), bernoulli_function(-60.0), bernoulli_function(55.0),
                          silicon_bandgap_ev(350.0), temperature_scaled_intrinsic_density(1e16, 300, 360, 1.12, 1.09), temperature_scaled_mobility(0.135, 300, 400, 2.42),
                          srh_recombination(1e22, 1e12, 1e16, rm), auger_recombination(1e22, 1e12, 1e16, rm), radiative_recombination(1e22, 1e12, 1e16, rm), total_recombination(1e18, 3e17, 1e16, rm),
                          caughey_thomas_mobility(3e23, 0.0068, 0.1414, 9.2e22, 0.711), high_field_mobility(0.135, 2e6, 1e5),
                          fermi_dirac_half(-9.0), fermi_dirac_half(0.3), fermi_dirac_half(12.0), fermi_dirac_half(45.0), inverse_fermi_dirac_half(2.5),
                          fermi_electron_density(2.8e25, -1.0), fermi_hole_density(1.04e25, 2.0), ionized_donor_density(1e23, -3.0, 0.045, 2.0, 300.0), ionized_acceptor_density(1e23, -2.0, 0.045, 4.0, 250.0)};
    vec("funcs", f);
    // PN equilibrium (upstream regression configuration).
    Semiconductor1DConfig pc; pc.length_m = 2e-6; pc.nodes = 101; pc.max_gummel_iterations = 1000; pc.relative_tolerance = 2e-7; pc.under_relaxation = 0.05;
    SemiconductorDevice1D pn(pc);
    pn.set_net_doping(pn_junction_doping(pn.x(), 1e-6, 1e21, 1e21));
    auto eq = pn.solve_equilibrium();
    vec("pn_eq", {(double)eq.converged, (double)eq.iterations, eq.relative_change, eq.terminal_current_a, pn.total_mobile_charge_c()});
    vec("pn_eq_phi", pn.potential()); vec("pn_n", pn.electron_density()); vec("pn_p", pn.hole_density());
    // Doped bar with SRH: DC, IV sweep, heat, AC, C-V and transient.
    Semiconductor1DConfig rc; rc.nodes = 41; rc.area_m2 = 2e-12; rc.max_gummel_iterations = 600; rc.relative_tolerance = 1e-9; rc.under_relaxation = 0.4; rc.right.voltage_v = -0.01;
    rc.recombination.radiative_m3_per_s = 1e-18;
    SemiconductorDevice1D bar(rc);
    bar.set_net_doping([](double x) { return 1e21 + 0.0 * x; });
    bar.initialize_charge_neutral();
    auto dc = bar.solve_dc();
    std::vector<double> bv0{(double)dc.converged, (double)dc.iterations, dc.terminal_current_a, bar.left_terminal_charge_c()};
    for (auto& p : bar.sweep_right_contact(-0.005, -0.02, 3)) bv0.insert(bv0.end(), {p.voltage_v, p.current_a, (double)p.converged});
    vec("bar", bv0);
    vec("bar_phi", bar.potential()); vec("bar_J", bar.total_current_density_faces()); vec("bar_heat", bar.joule_heating_cells()); vec("bar_E", bar.electric_field_faces());
    auto ss = bar.small_signal(1e6, 1e-5);
    vec("bar_ac", {ss.bias_voltage_v, ss.conductance_s, ss.capacitance_f, ss.admittance_s.real(), ss.admittance_s.imag(), ss.impedance_ohm.real(), ss.impedance_ohm.imag(), (double)ss.converged});
    std::vector<double> cvv; for (auto& c : bar.sweep_cv(-1e-3, 1e-3, 3, 1e-5)) cvv.insert(cvv.end(), {c.voltage_v, c.capacitance_f, c.electrode_charge_c, (double)c.converged});
    vec("bar_cv", cvv);
    std::vector<double> tr;
    bar.set_contact_voltages(0.0, -0.012);
    for (int i = 0; i < 3; ++i) { auto t = bar.step_transient(1e-12); tr.insert(tr.end(), {(double)t.converged, (double)t.iterations, t.time_s, t.conduction_current_a, t.displacement_current_a, t.terminal_current_a}); }
    vec("bar_tr", tr);
    // Temperature dependence + Schottky contact.
    Semiconductor1DConfig sc; sc.nodes = 11; sc.left.type = ContactType::schottky; sc.left.electron_barrier_ev = 0.72;
    SemiconductorDevice1D sd(sc);
    sd.set_net_doping([](double) { return 1e21; });
    sd.set_temperature(340.0);
    sd.initialize_charge_neutral();
    vec("schottky", {sd.config().material.intrinsic_density_m3, sd.config().material.electron_mobility_m2_per_vs, sd.config().material.hole_mobility_m2_per_vs, sd.config().material.conduction_band_density_m3, sd.electron_density().front(), sd.hole_density().front()});
    vec("schottky_phi", sd.potential());
    // Gate contact (MOS-like) quasi-static capacitance.
    Semiconductor1DConfig gc; gc.nodes = 31; gc.length_m = 0.5e-6; gc.left.type = ContactType::gate; gc.left.oxide_capacitance_f_per_m2 = 3.45e-3; gc.left.voltage_v = 0.2; gc.max_gummel_iterations = 2000; gc.under_relaxation = 0.1; gc.relative_tolerance = 1e-7;
    SemiconductorDevice1D gd(gc);
    gd.set_net_doping(std::vector<double>(31, -1e21));
    gd.initialize_charge_neutral();
    try { auto gdc = gd.solve_dc(); vec("gate", {(double)gdc.converged, (double)gdc.iterations, gd.left_terminal_charge_c(), gd.potential().front()}); } catch (const std::exception& e) { std::cerr << "gate: " << e.what() << "\n"; }
    vec("pin", pin_diode_doping(gd.x(), 1e-7, 3e-7, 1e23, 2e23, 1e19));
    MosCapacitor1D mc;
    auto mcv = mc.sweep_cv(-2.0, 2.0, 7);
    std::vector<double> mv{mc.oxide_capacitance_per_area()}; for (auto& p : mcv) mv.insert(mv.end(), {p.surface_potential_v, p.gate_charge_c, p.capacitance_f});
    vec("moscap", mv);
    LongChannelMosfetConfig lc; lc.body_effect_coefficient_sqrt_v = 0.4;
    LongChannelMosfet mos(lc);
    std::vector<double> mo{mos.threshold_voltage(-1.0), mos.drain_current(1.5, 0.3, 0.0), mos.drain_current(1.5, 2.0, -0.5), mos.drain_current(1.2, -0.4, 0.0), mos.drain_current(0.2, 1.0, 0.0)};
    auto op = mos.operating_point(1.4, 0.6, -0.2); mo.insert(mo.end(), {op.drain_current_a, op.transconductance_s, op.output_conductance_s, op.threshold_voltage_v});
    for (auto& s : mos.sweep(0.0, 2.0, 3, 0.0, 1.5, 3)) mo.push_back(s.drain_current_a);
    auto me = make_long_channel_mosfet_spice_evaluator(lc)(std::vector<double>{1.0, 1.6, 0.1, -0.2});
    mo.insert(mo.end(), me.terminal_current.begin(), me.terminal_current.end()); mo.insert(mo.end(), me.jacobian.begin(), me.jacobian.end());
    vec("mosfet", mo);
    BipolarJunctionTransistor1D bjt;
    auto d = bjt.derived(); auto d2 = bjt.derived(350.0);
    std::vector<double> bv{d.thermal_voltage_v, d.electron_diffusion_length_base_m, d.hole_diffusion_length_emitter_m, d.base_transport_factor, d.emitter_injection_efficiency, d.collector_injection_efficiency, d.forward_alpha, d.reverse_alpha, d.forward_beta, d.reverse_beta, d.emitter_saturation_current_a, d.collector_saturation_current_a, d.transport_saturation_current_a, d.emitter_base_built_in_v, d.collector_base_built_in_v, d2.forward_beta, d2.emitter_saturation_current_a};
    auto bo = bjt.operating_point(3.0, 0.7, 0.0); bv.insert(bv.end(), {bo.collector_current_a, bo.base_current_a, bo.emitter_current_a, bo.dissipated_power_w});
    for (double x : bjt.doping_profile(7)) bv.push_back(x);
    BipolarJunctionConfig pnpc; pnpc.polarity = BipolarPolarity::pnp;
    auto po = BipolarJunctionTransistor1D(pnpc).operating_point(-2.0, -0.68, 0.0, 320.0); bv.insert(bv.end(), {po.collector_current_a, po.base_current_a, po.emitter_current_a});
    auto be = make_bipolar_tcad_spice_evaluator({})(std::vector<double>{2.0, 0.65, 0.0});
    bv.insert(bv.end(), be.terminal_current.begin(), be.terminal_current.end()); bv.insert(bv.end(), be.jacobian.begin(), be.jacobian.end());
    vec("bjt", bv);
    IgbtPowerDevice ig;
    std::vector<double> iv2;
    auto io = ig.operating_point(5.0, 3.0); iv2.insert(iv2.end(), {io.collector_current_a, io.channel_current_a, io.bipolar_current_a, io.effective_drift_resistance_ohm, io.dissipated_power_w});
    auto io2 = ig.operating_point(5.0, 3.0, 360.0); iv2.insert(iv2.end(), {io2.collector_current_a, io2.effective_drift_resistance_ohm});
    for (int i = 0; i < 4; ++i) { auto s = ig.step(i < 2 ? 6.0 : 0.0, 2.5, 1e-6); iv2.insert(iv2.end(), {s.collector_current_a, s.bipolar_current_a, s.junction_temperature_k, s.dissipated_power_w}); }
    auto ie = make_igbt_tcad_spice_evaluator({})(std::vector<double>{3.0, 5.0, 0.0});
    iv2.insert(iv2.end(), ie.terminal_current.begin(), ie.terminal_current.end()); iv2.insert(iv2.end(), ie.jacobian.begin(), ie.jacobian.end());
    vec("igbt", iv2);
    GummelPoonBjt gp; MosLevel23 l23; VdmosPowerDevice vd; HystereticCoreModel core;
    std::vector<double> av;
    for (auto [a, b] : std::vector<std::pair<double, double>>{{0.7, -2.0}, {0.85, 0.1}, {-0.2, 0.6}}) { auto g = gp.operating_point(a, b); av.insert(av.end(), {g.collector_current_a, g.base_current_a, g.emitter_current_a, g.forward_transport_current_a, g.reverse_transport_current_a}); }
    for (auto [a, b] : std::vector<std::pair<double, double>>{{0.3, 1.0}, {1.5, 0.4}, {1.5, 3.0}}) { auto m = l23.operating_point(a, b); av.insert(av.end(), {m.drain_current_a, m.transconductance_s, m.output_conductance_s, (double)m.subthreshold}); }
    for (auto [a, b] : std::vector<std::pair<double, double>>{{2.0, 1.0}, {2.0, -0.7}}) { auto v = vd.operating_point(a, b); av.insert(av.end(), {v.drain_current_a, v.channel_current_a, v.body_diode_current_a, v.drift_drop_v}); }
    for (int i = 0; i < 12; ++i) { core.step(120.0 * std::sin(0.6 * i), 300.0 + 5.0 * i, 3e-5); av.insert(av.end(), {core.state().induction_t, core.state().loss_density_j_per_m3}); }
    vec("advanced", av);
    std::cout << "\n}\n";
}
