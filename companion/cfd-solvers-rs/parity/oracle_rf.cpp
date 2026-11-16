// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream rf/nport, rf/peec, rf/conductor_models and rf/modal_elements.
#include "cfd/rf/conductor_models.hpp"
#include "cfd/rf/modal_elements.hpp"
#include "cfd/rf/nport.hpp"
#include "cfd/rf/peec.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace cfd::rf;
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) {
        const double d = static_cast<double>(x);
        std::cout << (i++ ? "," : "");
        if (std::isfinite(d)) std::cout << std::setprecision(17) << d;
        else std::cout << (std::isnan(d) ? "\"nan\"" : d > 0 ? "\"inf\"" : "\"-inf\"");
    }
    std::cout << "]";
}
static void str(const std::string& k, const std::string& s) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":\"";
    first = false;
    for (char c : s) { if (c == '\n') std::cout << "\\n"; else std::cout << c; }
    std::cout << "\"";
}
static std::vector<double> flat(const ComplexMatrix& m) { std::vector<double> o; for (auto c : m.data()) { o.push_back(c.real()); o.push_back(c.imag()); } return o; }
static std::vector<double> flatc(const std::vector<Complex>& m) { std::vector<double> o; for (auto c : m) { o.push_back(c.real()); o.push_back(c.imag()); } return o; }
template <class T> static std::vector<double> cat(std::vector<double> a, const T& b) { a.insert(a.end(), b.begin(), b.end()); return a; }
static std::vector<double> cx(Complex c) { return {c.real(), c.imag()}; }
static ComplexMatrix make(std::size_t n, double seed) { ComplexMatrix m(n); for (std::size_t r = 0; r < n; ++r) for (std::size_t c = 0; c < n; ++c) m(r, c) = Complex{std::sin(seed + 3.0 * r + c) * (r == c ? 2.0 : 0.3), std::cos(seed * 0.7 + r - 2.0 * c) * 0.2}; return m; }

int main() {
    std::cout << "{\n";
    auto z = make(3, 0.4);
    for (std::size_t i = 0; i < 3; ++i) z(i, i) += Complex{60.0 + 10.0 * i, 5.0};
    std::vector<double> z0{50.0, 75.0, 30.0};
    auto s = z_to_s(z, z0);
    vec("z_to_s", flat(s));
    vec("s_to_z", flat(s_to_z(s, z0)));
    vec("renorm", flat(renormalize_s(s, z0, std::vector<double>{50.0, 50.0, 50.0})));
    std::vector<Complex> zc{{50.0, 10.0}, {75.0, -20.0}, {30.0, 0.0}};
    auto sp = z_to_s_power_wave(z, zc);
    vec("pw_s", flat(sp));
    vec("pw_z", flat(s_to_z_power_wave(sp, zc)));
    vec("pw_renorm", flat(renormalize_s_power_wave(sp, zc, std::vector<Complex>{{50, 0}, {50, 0}, {50, 0}})));
    vec("inverse", flat(inverse(make(4, 1.3))));

    std::istringstream t1("! comment\n# GHZ S MA R 50\n1.0 0.5 30 0.8 -45\n0.8 -45 0.4 10 ! row\n2.0 0.4 60 0.7 -90 0.7 -90 0.3 20\n");
    std::vector<double> tv;
    for (auto& p : read_touchstone(t1)) { tv.push_back(p.frequency_hz); tv = cat(tv, flat(p.s)); tv = cat(tv, p.reference_impedance); }
    std::istringstream t2("[Version] 2.0\n[Number of Ports] 3\n[Reference] 50 75 100\n# MHZ S DB\n[Network Data]\n100 -3 10 -20 30 -25 40 -20 30 -6 -10 -30 0 -25 40 -30 0 -9 90\n[End]\n");
    for (auto& p : read_touchstone(t2)) { tv.push_back(p.frequency_hz); tv = cat(tv, flat(p.s)); tv = cat(tv, p.reference_impedance); }
    vec("touchstone_read", tv);
    std::vector<NPortPoint> pts;
    for (int k = 0; k < 3; ++k) pts.push_back({1e9 * (k + 1) + 0.123, make(2, k * 0.5), {50.0, 50.0}});
    std::ostringstream w1, w2;
    write_touchstone(w2, pts, true);
    write_touchstone(w1, pts, false);
    str("touchstone_v2", w2.str());
    str("touchstone_v1", w1.str());

    std::vector<Complex> gam{{0.1, 20.0}, {0.05, 31.0}, {0.0, 12.0}};
    vec("shift", flat(shift_reference_planes(s, gam, std::vector<double>{0.01, -0.02, 0.005})));
    vec("rl", std::vector<double>{return_loss_db({0.3, 0.4}), vswr({0.3, 0.4}), vswr({0.6, 0.8})});
    Matrix2C s2{{0.6, -0.3}, {0.05, 0.02}, {2.5, 1.1}, {0.4, -0.2}};
    auto sm = stability_metrics(s2);
    auto sc = source_stability_circle(s2), lc = load_stability_circle(s2);
    vec("stability", cat(cat(std::vector<double>{sm.rollett_k, sm.mu}, cx(sm.determinant)), cat(cat(cx(sc.center), std::vector<double>{sc.radius, double(sc.valid)}), cat(cx(lc.center), std::vector<double>{lc.radius, double(lc.valid)}))));
    vec("gain", std::vector<double>{transducer_gain(s2), transducer_gain(s2, {0.2, -0.1}, {0.3, 0.25})});
    NoiseParameters np{1.4, {0.3, 0.2}, 12.0, 50.0};
    auto nc = noise_circle(np, 1.8);
    vec("noise", cat(std::vector<double>{noise_factor(np, {0.1, -0.2}), nc.radius, nc.noise_factor, double(nc.valid)}, cx(nc.center)));
    std::vector<double> lp;
    for (auto& q : sample_load_pull(s2, {0.1, 0.05}, 3, 5, 0.8)) { lp = cat(lp, cx(q.load_reflection)); lp.push_back(q.transducer_gain); }
    vec("load_pull", lp);
    auto big = make(3, 2.0); for (auto i = 0; i < 3; ++i) big(i, 0) *= 1.5;
    auto nq = network_quality(big);
    vec("quality", std::vector<double>{nq.maximum_singular_value, nq.reciprocity_error, double(nq.passive), double(nq.reciprocal)});
    vec("passivity", flat(enforce_passivity(big, 0.9)));
    std::vector<NPortPoint> sweep;
    for (int k = 0; k < 9; ++k) {
        double f = 1e8 * k;
        ComplexMatrix m(2);
        Complex e = std::exp(Complex{0.0, -2 * M_PI * f * 1.3e-9});
        m(0, 0) = 0.2 * e; m(1, 0) = 0.9 * e; m(0, 1) = 0.9 * e * std::exp(Complex{0.0, -0.2 * k}); m(1, 1) = Complex{0.1, 0.01 * k};
        sweep.push_back({f, m, {50.0, 50.0}});
    }
    auto bc = check_broadband_causality(sweep);
    vec("causality", std::vector<double>{bc.negative_time_energy_ratio, double(bc.causal), double(bc.impulse_samples)});
    std::vector<NPortPoint> rs;
    for (int k = 0; k < 12; ++k) {
        double f = 1e7 * std::pow(1.6, k);
        ComplexMatrix m(2);
        Complex sv{0.0, 2 * M_PI * f};
        m(0, 0) = 0.3 + 0.5 * 2e8 / (sv + 2e8); m(0, 1) = 0.1 / (sv / 1e9 + 1.0); m(1, 0) = m(0, 1); m(1, 1) = Complex{0.2, 0.0} - 0.1 * sv / (sv + 5e8);
        rs.push_back({f, m, {50.0, 50.0}});
    }
    auto model = fit_rational_network(rs, 4);
    vec("rational", cat(cat(model.poles_rad_per_s, flat(model.direct)), cat(flat(model.residues[2]), flat(model.evaluate(3.3e8)))));
    AdaptiveNetworkSweepConfig ac; ac.initial_points = 4; ac.maximum_points = 30; ac.relative_tolerance = 1e-2; ac.absolute_tolerance = 1e-4;
    auto ad = adaptive_network_sweep(1e6, 1e10, std::vector<double>{50.0}, [](double f) { ComplexMatrix m(1); m(0, 0) = 1.0 / (Complex{1.0, f / 3e8}); return m; }, ac);
    std::vector<double> av{double(ad.evaluations), double(ad.refinements), double(ad.converged), double(ad.points.size())};
    for (auto& p : ad.points) { av.push_back(p.frequency_hz); av = cat(av, cx(p.s(0, 0))); }
    vec("adaptive", av);
    vec("mixed", flat(single_ended_to_mixed_mode(make(4, 0.9))));
    auto l1 = synthesize_l_match(50.0, 200.0), l2 = synthesize_l_match(300.0, 50.0);
    vec("lmatch", std::vector<double>{l1.series_reactance_ohm, l1.shunt_susceptance_siemens, double(l1.shunt_first), l2.series_reactance_ohm, l2.shunt_susceptance_siemens, double(l2.shunt_first)});

    PeecFilamentSystem sys;
    sys.add_segment({{0, 0, 0}, {0.01, 0, 0}, 1e-4, 5.8e7});
    sys.add_segment({{0.01, 0, 0}, {0.01, 0.008, 0.001}, 1.2e-4, 5.8e7});
    sys.add_segment({{0.0, 0.002, 0}, {0.012, 0.003, 0}, 0.8e-4, 3.5e7});
    auto pm = sys.extract();
    auto pm4 = sys.extract(4);
    vec("peec_rl", cat(cat(pm.resistance_ohm, pm.partial_inductance_h), pm4.partial_inductance_h));
    auto pc = sys.extract_capacitance(2.2, 4);
    vec("peec_c", cat(pc.coefficient_of_potential_v_per_c, pc.capacitance_f));
    vec("peec_z", cat(flatc(sys.impedance_matrix(1e8)), flatc(sys.impedance_matrix_skin_effect(1e8, 2))));
    vec("peec_i", flatc(sys.solve_currents(5e7, std::vector<Complex>{{1, 0}, {0, 0.5}, {-0.2, 0}})));
    str("peec_spice", sys.to_spice_subcircuit("BUS"));
    FilamentSegment seg{{0, 0, 0}, {0.05, 0, 0}, 5e-4, 5.8e7};
    vec("wire", std::vector<double>{dc_resistance(seg), self_partial_inductance(seg), skin_depth_m(5.8e7, 1e6), round_wire_ac_resistance(seg, 1e8), round_wire_ac_resistance(seg, 0.0), round_wire_proximity_factor(3e-3, 5e-4, 1e8, 5.8e7)});
    auto gr = finite_ground_reflection(1e7, 0.01, 15.0, 0.6);
    vec("ground", cat(cx(gr.TE), cx(gr.TM)));
    std::vector<double> sub;
    for (auto& f : subdivide_rectangular_peec(seg, 2e-3, 1e-3, 3, 2)) { sub.push_back(f.start.y); sub.push_back(f.start.z); sub.push_back(f.end.x); sub.push_back(f.radius_m); }
    vec("subdivide", sub);
    auto md = modal_discontinuity_s({5.0, 12.0}, {0.002, -0.004}, 50.0);
    vec("modal", cat(cat(cx(md.a11), cx(md.a12)), cat(cx(md.a21), cx(md.a22))));
    std::cout << "\n}\n";
}
