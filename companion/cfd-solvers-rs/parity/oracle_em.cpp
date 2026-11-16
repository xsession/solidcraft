// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream cfd/em units.
#include "cfd/em/adaptivity3d.hpp"
#include "cfd/em/cable.hpp"
#include "cfd/em/edge_fem2d.hpp"
#include "cfd/em/edge_fem3d.hpp"
#include "cfd/em/frequency_domain.hpp"
#include "cfd/em/metasurface_gstc.hpp"
#include "cfd/em/photonic_modes.hpp"
#include "cfd/em/system.hpp"
#include "cfd/em/wave_port.hpp"
#include "cfd/fem/mesh2d.hpp"
#include "cfd/fem/mesh3d.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd::em;
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
template <class T> static std::vector<double> cat(std::vector<double> a, const T& b) { a.insert(a.end(), b.begin(), b.end()); return a; }
static std::vector<double> cx(Complex c) { return {c.real(), c.imag()}; }
template <class V> static std::vector<double> flatc(const V& v) { std::vector<double> o; for (auto c : v) { o.push_back(c.real()); o.push_back(c.imag()); } return o; }

int main() {
    std::cout << "{\n";
    auto modes = rectangular_waveguide_modes(22.86e-3, 10.16e-3, 12e9, 1.0, 1.0, 6);
    std::vector<double> wm;
    for (auto& m : modes) { wm.push_back(int(m.family)); wm.push_back(m.m); wm.push_back(m.n); wm.push_back(m.cutoff_frequency_hz); wm = cat(wm, cx(m.propagation_constant_rad_per_m)); wm = cat(wm, cx(m.wave_impedance_ohm)); wm.push_back(m.propagating); wm.push_back(m.electric_scale); }
    vec("wg_modes", wm);
    auto tf = rectangular_waveguide_transverse_fields(modes[0], 22.86e-3, 10.16e-3, 7e-3, 3e-3);
    auto tf2 = rectangular_waveguide_transverse_fields(modes[4], 22.86e-3, 10.16e-3, 7e-3, 3e-3);
    vec("wg_fields", cat(cat(flatc(tf.electric_v_per_m), flatc(tf.magnetic_a_per_m)), cat(flatc(tf2.electric_v_per_m), flatc(tf2.magnetic_a_per_m))));
    vec("wg_power", std::vector<double>{rectangular_mode_power_w(modes[0], 22.86e-3, 10.16e-3, 16, 16), rectangular_mode_power_w(modes[2], 22.86e-3, 10.16e-3, 8, 8)});

    MulticonductorRlcg rl{2, {0.1, 0.01, 0.01, 0.12}, {2.5e-7, 0.5e-7, 0.5e-7, 2.6e-7}, {1e-5, 0.0, 0.0, 1e-5}, {1.0e-10, -0.2e-10, -0.2e-10, 1.1e-10}};
    cfd::rf::ComplexMatrix load(2); load(0, 0) = Complex{50.0, 5.0}; load(1, 1) = Complex{75.0, 0.0}; load(0, 1) = load(1, 0) = Complex{2.0, 0.0};
    std::vector<Complex> vin{{1.0, 0.0}, {0.0, 0.2}};
    auto cr = solve_multiconductor_cable(3.0, 5, 2e7, rl, vin, load);
    vec("cable", cat(cat(flatc(cr.input_current_a), flatc(cr.load_voltage_v)), cat(flatc(cr.load_current_a), flatc(cr.transfer_matrix.data()))));
    std::vector<double> fs{1e5, 1e7};
    auto sw = sweep_multiconductor_cable(1.0, 3, fs, [&](double f) { auto r = rl; r.resistance[0] *= 1.0 + f * 1e-7; return r; }, vin, load);
    vec("cable_sweep", cat(flatc(sw[0].load_voltage_v), flatc(sw[1].load_voltage_v)));
    CableShieldTransferModel sh{0.01, 1e-9, 2e6};
    vec("shield", cat(cat(cx(cable_shield_transfer_impedance_ohm_per_m(5e6, sh)), cx(cable_shield_induced_voltage_v(5e6, {0.3, -0.1}, 2.0, sh))), cx(field_to_cable_open_circuit_voltage_v({10.0, 2.0}, 0.05, {0.8, 0.1}))));
    GstcSheetConfig g1; g1.electric_surface_admittance_s = {0.002, 0.004}; g1.magnetic_surface_impedance_ohm = {30.0, -60.0};
    GstcSheetConfig g2; g2.electric_surface_admittance_s = {0.0, -0.01}; g2.wave_impedance_after_ohm = 250.0;
    auto s1 = solve_gstc_normal_incidence(g1), s2 = solve_gstc_normal_incidence(g2);
    vec("gstc", cat(cat(cat(cx(s1.reflection), cx(s1.transmission)), std::vector<double>{s1.reflected_power, s1.transmitted_power, s1.absorbed_power}), cat(cat(cx(s2.reflection), cx(s2.transmission)), std::vector<double>{s2.absorbed_power})));

    FrequencyDomain1DConfig fc; fc.points = 21; fc.length_m = 2.0; fc.frequency_hz = 3e8; fc.relative_permittivity = 2.0; fc.conductivity_s_per_m = 1e-3;
    std::vector<Complex> j(21);
    for (int i = 0; i < 21; ++i) j[i] = Complex{std::exp(-0.1 * (i - 7) * (i - 7)), 0.05 * i};
    auto fd = solve_pec_driven_maxwell_1d(fc, j);
    vec("fd1d", cat(flatc(fd.electric_v_per_m), cat(flatc(fd.magnetic_a_per_m), std::vector<double>{fd.spacing_m})));
    std::vector<double> eps(30);
    for (int i = 0; i < 30; ++i) eps[i] = i < 15 ? 1.0 : 4.0;
    std::vector<double> cav;
    for (auto& m : pec_cavity_eigenmodes_1d(30, 1.0, eps, 1.0, 3)) { cav.push_back(m.angular_frequency_rad_per_s); cav.push_back(m.frequency_hz); cav = cat(cav, m.electric_shape); }
    for (auto& m : pec_cavity_eigenmodes_1d(12, 0.5, 2.5, 1.2, 2)) { cav.push_back(m.frequency_hz); cav = cat(cav, m.electric_shape); }
    vec("cavity", cav);

    std::vector<int> bits{0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0};
    std::vector<double> wave;
    double prev = 0.0;
    for (std::size_t k = 0; k < bits.size(); ++k) for (int s = 0; s < 8; ++s) { double target = bits[k] ? 1.0 : 0.0; prev += 0.45 * (target - prev); wave.push_back(prev + 0.02 * std::sin(0.7 * (k * 8 + s))); }
    auto eye = analyze_nrz_eye(wave, bits, 8, 0);
    auto eye2 = analyze_nrz_eye(wave, bits, 8, 6);
    vec("eye", std::vector<double>{eye.zero_mean, eye.one_mean, eye.zero_sigma, eye.one_sigma, eye.decision_threshold, eye.eye_height, eye.eye_width_ui, eye.estimated_ber, eye2.one_mean, eye2.estimated_ber});
    std::vector<DecouplingCapacitor> lib{{1e-6, 0.01, 1e-9, 1}, {1e-7, 0.02, 0.5e-9, 1}, {1e-8, 0.03, 0.3e-9, 2}};
    vec("decap", cat(cat(cx(decoupling_impedance_ohm(lib[0], 1e6)), cx(decoupling_impedance_ohm(lib[2], 0.0))), cat(cx(pdn_parallel_impedance_ohm(lib, 3e7, {0.001, 0.02})), cx(pdn_parallel_impedance_ohm({}, 1e6, {})))));
    std::vector<double> ff{1e5, 1e6, 1e7, 1e8}, tz{0.05, 0.02, 0.02, 0.05};
    auto opt = optimize_decoupling_greedy(ff, tz, lib, 6, {0.0, 0.0});
    std::vector<double> ov{opt.worst_ratio_before, opt.worst_ratio_after, double(opt.iterations)};
    for (auto& c : opt.selected) ov.push_back(c.capacitance_f);
    vec("decap_opt", ov);
    PdnGridConfig pg; pg.nx = 4; pg.ny = 3; pg.horizontal_resistance_ohm = 2e-3; pg.vertical_resistance_ohm = 3e-3; pg.relative_tolerance = 1e-10;
    std::vector<double> loads(12, 0.0); loads[5] = 1.5; loads[10] = 0.7; loads[3] = -0.2;
    auto ir = solve_pdn_ir_drop_grid(pg, loads, std::vector<PdnVoltageSource>{{0, 1.0}, {11, 0.99}});
    vec("pdn", cat(ir.voltage_v, std::vector<double>{ir.minimum_voltage_v, ir.maximum_drop_v, double(ir.iterations), double(ir.converged)}));
    auto dexp = sample_double_exponential_pulse(50e-9, 2e-9, 2.0, 1e-9, 10e-9);
    auto ds = sample_damped_sine(10e-9, 0.7e-9, 1.5, 2e8, 4e-9);
    auto pr = probe_waveform(ds, 0.7e-9);
    vec("pulses", cat(cat(dexp, ds), std::vector<double>{pr.peak, pr.rms, pr.impulse, pr.energy, normalized_double_exponential_pulse(-1.0, 1.0, 1.0, 2.0)}));

    auto mesh2 = cfd::fem::make_rectangle_tri_mesh(3, 2, 1.0, 0.6);
    EdgeMaxwell2DConfig c2; c2.frequency_hz = 2e8; c2.relative_permittivity = 2.0; c2.conductivity_s_per_m = 0.01;
    auto r2 = solve_driven_edge_maxwell_2d(mesh2, c2, [](cfd::fem::Node2 p) { return ComplexVec2{Complex{p.y, 0.1}, Complex{-p.x, 0.0}}; });
    std::vector<double> e2;
    for (auto& e : r2.edges) { e2.push_back(e[0]); e2.push_back(e[1]); }
    e2 = cat(e2, flatc(r2.edge_voltage_v));
    for (auto& f : r2.electric_centroid_v_per_m) e2 = cat(e2, flatc(f));
    vec("edge2d", e2);

    auto mesh3 = cfd::fem::make_box_tet_mesh(2, 2, 2, 1.0, 0.8, 0.6);
    EdgeMaxwell3DConfig c3; c3.frequency_hz = 1.5e8; c3.relative_permittivity = 2.0; c3.relative_permittivity_z = 3.0; c3.conductivity_s_per_m = 0.05; c3.linear.relative_tolerance = 1e-12;
    std::vector<EdgeCurrentPort3D> ports{{13, 4, Complex{1.0, 0.0}}};
    auto r3 = solve_driven_edge_maxwell_3d(mesh3, c3, [](cfd::fem::Point3 p) { return ComplexVec3{Complex{p.z, 0.0}, Complex{0.0, p.x}, Complex{0.2, 0.0}}; }, ports);
    std::vector<double> e3;
    e3 = cat(e3, flatc(r3.edge_voltage_v));
    for (auto& f : r3.electric_centroid_v_per_m) e3 = cat(e3, flatc(f));
    for (auto& f : r3.magnetic_centroid_a_per_m) e3 = cat(e3, flatc(f));
    vec("edge3d", e3);
    std::size_t port_edge = 0;
    for (std::size_t i = 0; i < r3.edges.size(); ++i) if (r3.edges[i][0] == 4 && r3.edges[i][1] == 13) port_edge = i;
    auto th = edge_port_thevenin(r3, port_edge, {1.0, 0.0});
    auto q = resonator_quality_3d(mesh3, c3, r3, 5.8e7, 1.0);
    vec("edge3d_post", cat(cat(cx(th.open_circuit_voltage_v), cx(th.port_impedance_ohm)), std::vector<double>{double(port_edge), q.electric_energy_j, q.magnetic_energy_j, q.dielectric_loss_w, q.conductor_loss_w, q.quality_factor}));
    auto ind = maxwell_face_jump_indicators_3d(mesh3, r3);
    auto mk = mark_maxwell_dorfler(ind, 0.4);
    auto refined = refine_tet4_marked_longest_edges(mesh3, mk);
    std::vector<double> rf;
    for (auto& p : refined.nodes) { rf.push_back(p.x); rf.push_back(p.y); rf.push_back(p.z); }
    for (auto& t : refined.tetrahedra) for (auto n : t.node) rf.push_back(n);
    for (auto b : refined.boundary_node) rf.push_back(b);
    vec("adapt", cat(cat(ind, mk), rf));

    PhotonicGrid2D pgd{4, 4, 1e-6, 1e-6, {}, {}};
    for (int i = 0; i < 16; ++i) pgd.relative_permittivity.push_back((i % 4 == 1 || i / 4 == 2) ? 9.0 : 1.0 + 0.1 * i);
    PhotonicModeConfig pc; pc.mode_count = 3;
    std::vector<double> pb;
    for (auto& m : solve_photonic_band_modes(pgd, {M_PI / 2e-6, 0.3e6}, pc)) { pb.push_back(m.eigenvalue_per_m2); pb.push_back(m.frequency_hz); pb.push_back(m.relative_residual); pb.push_back(m.iterations); pb.push_back(m.converged); pb = cat(pb, flatc(m.field)); }
    PhotonicGrid2D pdisp{3, 3, 1e-6, 1e-6, {}, [](std::size_t i, double f) { return 2.0 + 0.1 * i + 1e-15 * f; }};
    PhotonicModeConfig pc2; pc2.mode_count = 2; pc2.relative_tolerance = 1e-9;
    for (auto& m : solve_photonic_band_modes(pdisp, {1e6, 0.0}, pc2)) { pb.push_back(m.frequency_hz); pb.push_back(m.iterations); pb.push_back(m.relative_residual); }
    vec("photonic_band", pb);
    PhotonicGrid2D wg{5, 4, 2e-6, 1.5e-6, {}, {}};
    for (int i = 0; i < 20; ++i) wg.relative_permittivity.push_back((i % 5 >= 1 && i % 5 <= 3 && i / 5 >= 1 && i / 5 <= 2) ? 12.0 : 2.1);
    PhotonicWaveguideConfig wc; wc.frequency_hz = 1.9e14; wc.mode_count = 3;
    std::vector<double> pw;
    for (auto& m : solve_photonic_waveguide_modes(wg, wc)) { pw.push_back(m.transverse_eigenvalue_per_m2); pw.push_back(m.propagation_constant_rad_per_m); pw.push_back(m.effective_index); pw.push_back(m.relative_residual); pw.push_back(m.guided); pw = cat(pw, flatc(m.field)); }
    vec("photonic_wg", pw);
    std::cout << "\n}\n";
}
