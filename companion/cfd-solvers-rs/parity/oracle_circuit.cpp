// Test oracle only. Exercises upstream circuit/{spice, analysis, raw, osdi_adapter,
// external_spice_adapter, spice_compat}.
#include "cfd/circuit/analysis.hpp"
#include "cfd/circuit/external_spice_adapter.hpp"
#include "cfd/circuit/osdi_adapter.hpp"
#include "cfd/circuit/raw.hpp"
#include "cfd/circuit/spice.hpp"
#include "cfd/circuit/spice_compat.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace cfd::circuit;
static bool first = true;
static void key(const std::string& k) { std::cout << (first ? "" : ",\n") << "\"" << k << "\":"; first = false; }
static void vec(const std::string& k, const std::vector<double>& v) {
    key(k); std::cout << "[";
    for (std::size_t i = 0; i < v.size(); ++i) {
        std::cout << (i ? "," : "");
        if (std::isfinite(v[i])) std::cout << std::setprecision(17) << v[i];
        else std::cout << (std::isnan(v[i]) ? "\"nan\"" : (v[i] > 0 ? "\"inf\"" : "\"-inf\""));
    }
    std::cout << "]";
}
static std::string esc(const std::string& s) { std::string o; for (char c : s) { if (c == '"' || c == '\\') o += '\\'; if (c == '\n') { o += "\\n"; continue; } if (c == '\t') { o += "\\t"; continue; } o += c; } return o; }
static void str(const std::string& k, const std::string& v) { key(k); std::cout << "\"" << esc(v) << "\""; }
static void op(std::vector<double>& g, const OperatingPoint& p) { g.insert(g.end(), p.node_voltage.begin(), p.node_voltage.end()); g.insert(g.end(), p.branch_current.begin(), p.branch_current.end()); g.push_back(p.iterations); g.push_back(p.converged); }
static void cx(std::vector<double>& g, Complex c) { g.push_back(c.real()); g.push_back(c.imag()); }
static void acs(std::vector<double>& g, const AcSolution& a) { g.push_back(a.frequency_hz); for (auto c : a.node_voltage) cx(g, c); for (auto c : a.branch_current) cx(g, c); }
static void tr(std::vector<double>& g, const std::vector<TransientPoint>& pts) { for (auto& p : pts) { g.push_back(p.time); g.insert(g.end(), p.node_voltage.begin(), p.node_voltage.end()); g.insert(g.end(), p.branch_current.begin(), p.branch_current.end()); } }
static std::string err(auto f) { try { f(); return ""; } catch (const std::exception& e) { return e.what(); } }

int main() {
    std::cout << "{\n";
    NewtonConfig cfg;
    // numbers and expressions
    { std::vector<double> g; for (auto t : {"1k", "2.2MEG", "3.3m", "10u", "47n", "1p", "5f", "2T", "1.5G", "1e3", " 4.7K ", "1MIL", "0.5e-3k"}) g.push_back(parse_spice_number(t));
      std::unordered_map<std::string, double> p{{"rload", 2.0e3}, {"Gain", 3.0}};
      for (auto e : {"1+2*3", "{(1+2)*3}", "2^3^2", "-2^2", "rload/2", "GAIN*1k+1MEG", "1e-3*2", "((4))", "2.5u*4"}) g.push_back(evaluate_spice_expression(e, p));
      vec("numbers", g); }
    str("expr_errors", err([] { (void)evaluate_spice_expression("1/0"); }) + "|" + err([] { (void)evaluate_spice_expression("foo+1"); }) + "|" + err([] { (void)evaluate_spice_expression("(1+2"); }) + "|" + err([] { (void)evaluate_spice_expression("1 2"); }));
    // linear netlist with controlled sources
    auto lin = Circuit::parse_spice(
        "* linear\n.param rv=2k gain=4\nV1 in 0 DC 10\nR1 in a {rv}\nR2 a 0 3k\nE1 b 0 a 0 {gain}\nR3 b c 1k\nVsense c 0 0\nF1 0 d Vsense 2\nR4 d 0 500\nH1 e 0 Vsense 100\nR5 e 0 1k\nG1 0 f a 0 1m\nR6 f 0 2k\nI1 0 g 1m\nR7 g 0 4.7k\n");
    { std::vector<double> g; op(g, lin.dc_operating_point(cfg)); op(g, lin.dc_operating_point_homotopy(cfg)); vec("linear_dc", g); }
    // diode / MOS / BJT / JFET / switch
    auto nl = Circuit::parse_spice(
        "V1 vdd 0 5\nR1 vdd a 1k\nD1 a 0 DMOD\n.model DMOD D(IS=1e-14 N=1.5)\n"
        "Vg g 0 2\nRd vdd dr 2k\nM1 dr g 0 0 NM\n.model NM NMOS(VTO=0.8 KP=2m LAMBDA=0.02)\n"
        "Vb bb 0 0.7\nRb bb b 10k\nRc vdd c 1k\nQ1 c b 0 QN\n.model QN NPN(IS=1e-15 BF=150 BR=2)\n"
        "Vj jg 0 -0.5\nRj vdd jd 3k\nJ1 jd jg 0 JN\n.model JN NJF(VTO=-2 BETA=1.5m LAMBDA=0.01)\n"
        "Rs vdd s 100\nS1 s 0 g 0 SW1\n.model SW1 SW(RON=10 ROFF=1MEG VT=1 VH=0.2)\n"
        "Vp vdd2 0 -5\nRp vdd2 pd 2k\nM2 pd pg 0 0 PM\nVpg pg 0 -2\n.model PM PMOS(VTO=-0.7 KP=1m)\n");
    { std::vector<double> g; op(g, nl.dc_operating_point(cfg)); op(g, nl.dc_operating_point_homotopy(cfg)); vec("nonlinear_dc", g); }
    { std::vector<double> g; for (double f : {1.0, 1e3, 1e6}) acs(g, nl.ac(f)); vec("nonlinear_ac", g); }
    // AC: RLC with mutual inductance and AC sources
    auto rlc = Circuit::parse_spice("V1 in 0 0 AC 1 30\nR1 in a 50\nL1 a b 1m\nC1 b 0 1u\nL2 c 0 2m\nK1 L1 L2 0.6\nR2 c 0 75\nI2 0 b 0 AC 1m\n");
    { std::vector<double> g; for (const auto& a : rlc.ac_log_sweep(100.0, 1e5, 7)) acs(g, a); vec("rlc_ac", g); }
    // Transient: RC pulse with BE/TRAP/BDF2, RL sine with mutual
    auto rc = Circuit::parse_spice("V1 in 0 PULSE(0 5 1u 1u 1u 5u 20u)\nR1 in out 1k\nC1 out 0 2n\nL1 out x 10m\nR2 x 0 2k\nI1 0 x SIN(0 1m 50k 2u 1e3 30)\n.ic V(out)=1\n");
    for (auto [name, m] : std::vector<std::pair<std::string, TransientMethod>>{{"tr_be", TransientMethod::backward_euler}, {"tr_trap", TransientMethod::trapezoidal}, {"tr_bdf2", TransientMethod::bdf2}}) { std::vector<double> g; tr(g, rc.transient(0.5e-6, 40, m, cfg)); vec(name, g); }
    auto xf = Circuit::parse_spice("V1 in 0 SIN(0 1 10k)\nR1 in a 10\nL1 a 0 1m\nL2 b 0 4m\nK12 L1 L2 0.9\nR2 b 0 100\nC2 b 0 10n\n");
    { std::vector<double> g; tr(g, xf.transient(2e-6, 30, TransientMethod::trapezoidal, cfg)); tr(g, xf.transient(2e-6, 30, TransientMethod::bdf2, cfg)); vec("tr_mutual", g); }
    auto rect = Circuit::parse_spice("V1 in 0 PWL(0 0 10u 5 20u -5 30u 5)\nD1 in out DX\n.model DX D(IS=1e-12)\nR1 out 0 1k\nC1 out 0 100n\n");
    { std::vector<double> g; tr(g, rect.transient(1e-6, 30, TransientMethod::backward_euler, cfg)); tr(g, rect.transient(1e-6, 30, TransientMethod::trapezoidal, cfg)); vec("tr_rectifier", g); }
    // adaptive
    for (auto [name, m] : std::vector<std::pair<std::string, TransientMethod>>{{"adapt_be", TransientMethod::backward_euler}, {"adapt_bdf2", TransientMethod::bdf2}}) {
        AdaptiveTransientConfig ac; ac.method = m; ac.initial_step = 1e-7; ac.maximum_step = 2e-6; ac.relative_tolerance = 1e-3; ac.absolute_tolerance = 1e-5;
        auto r = rc.transient_adaptive(1.5e-5, ac, cfg);
        std::vector<double> g{(double)r.accepted_steps, (double)r.rejected_steps, r.minimum_accepted_step, r.maximum_accepted_step, (double)r.points.size()}; tr(g, r.points); vec(name, g);
    }
    { AdaptiveTransientConfig ac; ac.method = TransientMethod::bdf2; ac.initial_step = 1e-6; auto r = rect.transient_adaptive(3e-5, ac, cfg); std::vector<double> g{(double)r.accepted_steps, (double)r.rejected_steps}; tr(g, r.points); vec("adapt_rect", g); }
    // noise
    { auto n1 = nl.output_noise(1e3, "a", 300.0, cfg); auto n2 = lin.output_noise(1e4, "b", 350.0, cfg); std::vector<double> g{n1.output_noise_v_per_sqrt_hz, n1.output_noise_density_v2_per_hz}; for (auto& c : n1.contributions) g.push_back(c.output_noise_density_v2_per_hz); g.push_back(n2.output_noise_density_v2_per_hz); for (auto& c : n2.contributions) g.push_back(c.output_noise_density_v2_per_hz); vec("noise", g);
      std::string names; for (auto& c : n1.contributions) names += c.source + ","; str("noise_names", names); }
    // subcircuits, params, funcs, behavioral, table, poly, ic/nodeset
    auto hier = Circuit::parse_spice(
        ".param vin=3 rbase=1k\n.func half(x) {x/2}\n.func twice(x) {2*x}\n"
        ".subckt divider top bot mid params: ra=1k rb=1k\nR1 top mid {ra}\nR2 mid bot {rb}\n.ends\n"
        ".subckt stage in out params: g=2\n.model DS D(IS=2e-14)\nX1 in 0 m divider ra={g*1k} rb=2k\nE1 out 0 m 0 {twice(g)}\nD1 out 0 DS\nR1 out 0 10k\n.ends\n"
        "V1 in 0 {vin}\nXs in out stage g=1.5\nB1 bo 0 V={half(V(in))+V(out,0)*0.1+TIME}\nRbo bo 0 1k\n"
        "Gp 0 po POLY(1) in 0 1m 0.5m 0.1m\nRpo po 0 {rbase}\nEt to 0 TABLE {V(in)} = (0,0) (1,2) (4,5)\nRto to 0 1k\nBi 0 bi I={V(in)*1m}\nRbi bi 0 2k\n.nodeset V(out)=0.6\n.ic V(bo)=0.5\n");
    { std::vector<double> g; op(g, hier.dc_operating_point(cfg)); op(g, hier.dc_operating_point_homotopy(cfg)); vec("hier_dc", g); std::string n; for (auto& s : hier.node_names()) n += s + ","; str("hier_nodes", n); }
    { std::vector<double> g; acs(g, hier.ac(2e3)); vec("hier_tr_err", {0.0}); str("hier_tr_msg", err([&] { (void)hier.transient(1e-3, 5, cfg); })); vec("hier_ac_tr", g); }
    str("parse_errors", err([] { (void)Circuit::parse_spice(".include x.cir\n"); }) + "|" + err([] { (void)Circuit::parse_spice("X1 a b nosuch\n"); }) + "|" + err([] { (void)Circuit::parse_spice(".subckt a b\nR1 b 0 1\n"); }) + "|" + err([] { (void)Circuit::parse_spice("B1 a 0 Q=1\n"); }) + "|" + err([] { (void)Circuit::parse_spice(".ic V(a,b)=1\nR1 a 0 1\n"); }));
    // N-ports
    { Circuit c; auto in = c.node("in"), out = c.node("out"); c.add_voltage_source("V1", in, 0, 0.0, Complex{1.0, 0.0}); c.add_resistor("Rs", in, c.node("p1"), 50.0);
      std::vector<double> freqs{1e8, 1e9, 3e9}; c.add_microstrip_line("TL1", {c.find_node("p1"), 0}, {out, 0}, 1.9e-3, 1e-3, 4.4, 0.05, freqs);
      c.add_resistor("Rl", out, 0, 75.0); std::vector<double> g; for (double f : {1e8, 5e8, 2e9, 4e9}) acs(g, c.ac(f));
      Circuit w; auto a = w.node("a"), b = w.node("b"); w.add_current_source("I1", 0, a, 0.0, Complex{1e-3, 0.0}); w.add_resistor("Ra", a, 0, 50.0); w.add_resistor("Rb", b, 0, 50.0);
      std::vector<double> wf{8e9, 1e10, 1.2e10}; w.add_rectangular_waveguide_te10_line("WG", {a, 0}, {b, 0}, 22.86e-3, 10.16e-3, 1.0, 1.0, 0.1, wf);
      w.add_stripline_line("SL", {a, 0}, {b, 0}, 1e-3, 2e-3, 2.2, 0.02, wf); w.add_coplanar_waveguide_line("CPW", {a, 0}, {b, 0}, 0.5e-3, 0.2e-3, 9.8, 0.01, wf);
      for (double f : {8e9, 9e9, 1.3e10}) acs(g, w.ac(f)); vec("nports", g); }
    // static, dynamic and electrothermal devices
    { Circuit c; auto a = c.node("a"), b = c.node("b"), th = c.node("th");
      c.add_voltage_source("V1", a, 0, 2.0, Complex{0.1, 0.0}, PulseWaveform{0.0, 2.0, 0.0, 1e-6, 1e-6, 5e-6, 0.0});
      c.add_resistor("R1", a, b, 100.0);
      c.add_static_device("NL", {b, 0}, [](std::span<const double> v) { const double x = v[0] - v[1]; const double i = 1e-3 * x + 2e-4 * x * x * x; const double gd = 1e-3 + 6e-4 * x * x; return StaticDeviceEvaluation{{i, -i}, {gd, -gd, -gd, gd}}; });
      c.add_dynamic_device("CQ", {b, 0}, [](std::span<const double> v) { const double x = v[0] - v[1]; const double q = 1e-9 * x + 2e-10 * x * x; const double cq = 1e-9 + 4e-10 * x; return DynamicDeviceEvaluation{{1e-6 * x, -1e-6 * x}, {1e-6, -1e-6, -1e-6, 1e-6}, {q, -q}, {cq, -cq, -cq, cq}}; });
      ElectroThermalDeviceConfig et; et.thermal_resistance_k_per_w = 50.0; et.thermal_capacitance_j_per_k = 1e-6;
      c.add_electrothermal_device("ET", {b, 0}, th, et, [](std::span<const double> v, double t) { const double x = v[0] - v[1]; const double g0 = 2e-3 * (1.0 - 0.004 * (t - 300.0)); const double i = g0 * x; ElectroThermalDeviceEvaluation e; e.electrical_current = {i, -i}; e.electrical_jacobian = {g0, -g0, -g0, g0}; const double dg = -2e-3 * 0.004 * x; e.electrical_temperature_derivative = {dg, -dg}; e.dissipated_power_w = i * x; e.power_voltage_derivative = {2 * g0 * x, -2 * g0 * x}; e.power_temperature_derivative_w_per_k = dg * x; return e; });
      std::vector<double> g; op(g, c.dc_operating_point(cfg)); acs(g, c.ac(1e5)); tr(g, c.transient(1e-6, 12, TransientMethod::bdf2, cfg)); tr(g, c.transient(1e-6, 8, TransientMethod::backward_euler, cfg));
      AdaptiveTransientConfig ac; ac.method = TransientMethod::bdf2; ac.initial_step = 2e-7; auto r = c.transient_adaptive(4e-6, ac, cfg); g.push_back(r.accepted_steps); g.push_back(r.rejected_steps); tr(g, r.points);
      vec("devices", g);
      str("trap_dynamic_err", err([&] { (void)c.transient(1e-6, 2, TransientMethod::trapezoidal, cfg); })); }
    // sweeps and analysis helpers
    { std::vector<double> g; auto copy = nl; for (auto& p : copy.dc_sweep_voltage_source("Vg", 0.0, 3.0, 4, cfg)) { g.push_back(p.source_value); op(g, p.operating_point); }
      for (auto& p : nl.dc_temperature_sweep(250.0, 400.0, 3, cfg)) { g.push_back(p.temperature_k); op(g, p.operating_point); }
      auto s = voltage_source_dc_sensitivity(nl, "Vg", "dr", 1e-3, cfg); g.push_back(s.baseline_output); g.push_back(s.derivative);
      auto d = small_signal_distortion(nl, "Vb", "c", 1e-3, 1e-3, cfg); for (double v : {d.bias_output_v, d.first_derivative, d.second_derivative, d.third_derivative, d.fundamental_v, d.second_harmonic_v, d.third_harmonic_v, d.hd2, d.hd3}) g.push_back(v);
      std::vector<ResistorTolerance> tol{{"R1", 0.05}, {"Rd", 0.1}}; auto mc = monte_carlo_dc_resistors(nl, "dr", tol, 25, 7, cfg); for (double v : {mc.mean, mc.standard_deviation, mc.minimum, mc.maximum}) g.push_back(v);
      for (double v : {0.9, 0.2, -0.3}) g.push_back(limit_pn_junction_voltage(v, 0.6, 0.025852, 1e-14)); g.push_back(limit_pn_junction_voltage(5.0, -1.0, 0.025852, 1e-14));
      std::vector<double> rs{100.0, 200.0}; for (auto& p : parameter_sweep_dc(nl, rs, [](Circuit& c, double v) { c.set_resistance("R1", v); }, cfg)) op(g, p.operating_point);
      for (auto& p : temperature_sweep_dc(nl, std::vector<double>{280.0, 320.0}, cfg)) op(g, p.operating_point);
      vec("analysis_dc", g); }
    { auto pi = Circuit::parse_spice("R1 a 0 100\nR2 a b 50\nC1 b 0 1n\nL1 b c 1u\nR3 c 0 200\n");
      std::vector<double> g; auto s = two_port_s_parameters(pi, {"a"}, {"c"}, 1e7); for (auto c : {s.a11, s.a12, s.a21, s.a22}) cx(g, c);
      std::vector<CircuitPort> ports{{pi.find_node("a"), 0, 50.0}, {pi.find_node("c"), 0, 75.0}}; auto e = extract_s_parameters(pi, ports, 1e7, cfg); for (std::size_t i = 0; i < 2; ++i) for (std::size_t j = 0; j < 2; ++j) cx(g, e(i, j));
      auto lp = Circuit::parse_spice("V1 in 0 0 AC 1\nR1 in a 1k\nC1 a 0 100n\nR2 a b 10k\nC2 b 0 1n\n");
      auto dp = estimate_dominant_pole(lp, "b", 10.0, 1e6, 41, cfg); g.push_back(dp.dc_gain); g.push_back(dp.pole_frequency_hz); g.push_back(dp.found);
      PoleZeroConfig pz; pz.start_hz = 10.0; pz.stop_hz = 1e7; pz.samples = 40; pz.numerator_order = 0; auto r = pole_zero_analysis(lp, "b", pz, cfg); for (auto p : r.poles_rad_per_s) cx(g, p); for (auto z : r.zeros_rad_per_s) cx(g, z); g.push_back(r.relative_rms_fit_error); g.push_back(r.poles_rad_per_s.size()); g.push_back(r.zeros_rad_per_s.size());
      auto ns = output_noise_spectrum(nl, "a", 10.0, 1e5, 5, 300.0, cfg); g.push_back(ns.integrated_rms_voltage); for (auto& p : ns.points) g.push_back(p.output_noise_density_v2_per_hz);
      vec("analysis_ac", g); }
    { auto amp = Circuit::parse_spice("V1 in 0 SIN(0 0.5 1k)\nR1 in a 1k\nD1 a 0 DD\nD2 0 a DD\n.model DD D(IS=1e-12)\nR2 a 0 10k\nC1 a 0 10n\n");
      auto t = amp.transient(1e-5, 200, TransientMethod::trapezoidal, cfg); std::vector<double> g; auto a = amp.find_node("a");
      for (double f : {0.0, 1e3, 2e3, 3e3}) { auto m = fourier_measurement(t, a, f); cx(g, m.phasor); g.push_back(m.magnitude); g.push_back(m.phase_rad); }
      g.push_back(total_harmonic_distortion(t, a, 1e3, 5));
      PeriodicSteadyStateConfig pss; pss.period_s = 1e-3; pss.samples_per_period = 40; pss.max_periods = 8; pss.relative_tolerance = 1e-2; pss.absolute_tolerance = 1e-4;
      auto p = periodic_steady_state(amp, pss, cfg); g.push_back(p.periods); g.push_back(p.normalized_residual); g.push_back(p.converged); tr(g, p.period);
      vec("analysis_tr", g); }
    // sparse path
    { std::ostringstream net; net << "V1 n0 0 1\n"; for (int i = 0; i < 80; ++i) net << "R" << i << " n" << i << " n" << i + 1 << " " << 10 + i << "\nRg" << i << " n" << i + 1 << " 0 " << 1000 + 7 * i << "\n";
      auto big = Circuit::parse_spice(net.str()); std::vector<double> g; op(g, big.dc_operating_point(cfg)); NewtonConfig dense = cfg; dense.sparse_mna_threshold = 0; op(g, big.dc_operating_point(dense)); vec("sparse", g); }
    // raw
    { auto t = rc.transient(1e-6, 4, cfg); auto ds = raw_from_transient(t, rc.node_names()); std::ostringstream o; write_spice_raw(o, ds); str("raw_tr", o.str());
      std::vector<AcSolution> a{rlc.ac(1e3), rlc.ac(2e3)}; auto da = raw_from_ac(a, rlc.node_names()); std::ostringstream o2; write_spice_raw(o2, da); str("raw_ac", o2.str());
      std::istringstream in(o2.str()); auto back = read_spice_raw(in); std::vector<double> g{(double)back.traces.size(), (double)back.complex_values}; for (auto& tr : back.traces) for (auto c : tr.values) cx(g, c); vec("raw_back", g); }
    // plans
    { auto p = make_openvaf_compile_plan("bsim.va", "bsim.osdi"); auto n = make_external_spice_run_plan(ExternalSpiceEngine::ngspice, "a.cir", "a.raw", {"-x"}); auto x = make_external_spice_run_plan(ExternalSpiceEngine::xyce, "a.cir", "a.raw", {"-r", "a.raw"});
      std::string s; for (auto& a : p.argv) s += a + " "; s += "|"; for (auto& a : n.argv) s += a + " "; s += "|"; for (auto& a : x.argv) s += a + " "; s += "|"; for (auto d : spice_compatibility_profile(SpiceCompatibilityMode::ltspice).supported_directives) s += std::string(d) + " "; str("plans", s); }
    std::cout << "\n}\n";
}
