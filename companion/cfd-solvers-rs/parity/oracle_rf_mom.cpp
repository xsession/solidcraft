// Test oracle only. Exercises upstream rf/{thin_wire_mom, antenna, surface_mom, integral_equation_extensions, asymptotic}.
#include "cfd/rf/antenna.hpp"
#include "cfd/rf/asymptotic.hpp"
#include "cfd/rf/integral_equation_extensions.hpp"
#include "cfd/rf/surface_mom.hpp"
#include "cfd/rf/thin_wire_mom.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace cfd::rf;
using Cx = std::complex<double>;
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
static void cv(std::vector<double>& o, Cx z) { o.push_back(z.real()); o.push_back(z.imag()); }
static std::vector<double> cvec(const std::vector<Cx>& z) { std::vector<double> o; for (auto x : z) cv(o, x); return o; }
static std::vector<double> cmat(const ComplexMatrix& m) { std::vector<double> o; for (std::size_t i = 0; i < m.size(); ++i) for (std::size_t j = 0; j < m.size(); ++j) cv(o, m(i, j)); return o; }
int main() {
    std::cout << "{\n";
    ThinWireMomConfig d; d.length_m = 0.5; d.radius_m = 1e-3; d.frequency_hz = 290e6; d.segments = 21;
    auto r = solve_center_fed_thin_wire(d);
    std::vector<double> o; cv(o, r.feed_impedance_ohm); o.push_back(r.feed_current_a);
    vec("dipole", o); vec("dipole_I", cvec(r.current_a)); vec("dipole_z", r.segment_center_z_m);
    // L-shaped junction antenna with losses and loads over PEC ground.
    OrientedWireMomConfig oc; oc.frequency_hz = 400e6; oc.junction_tolerance_m = 1e-6;
    oc.ground_model = ThinWireGroundModel::pec_plane; oc.ground_plane_z_m = 0.0;
    OrientedThinWire w0; w0.start_m = {0, 0, 0.05}; w0.end_m = {0, 0, 0.3}; w0.radius_m = 1e-3; w0.segments = 9; w0.conductivity_s_per_m = 5.8e7;
    OrientedThinWire w1; w1.start_m = {0, 0, 0.3}; w1.end_m = {0.2, 0.05, 0.3}; w1.radius_m = 1.5e-3; w1.segments = 7; w1.dielectric_loss_ohm_per_m = 0.3;
    OrientedThinWire w2; w2.start_m = {0.1, 0.2, 0.05}; w2.end_m = {0.1, 0.2, 0.35}; w2.radius_m = 1e-3; w2.segments = 9;
    oc.wires = {w0, w1, w2};
    oc.feeds = {{0, 1, Cx{1.0, 0.0}}};
    oc.loads = {{2, 4, Cx{25.0, 10.0}}};
    ThinWireTransmissionLineLoad tl; tl.first_wire = 1; tl.first_segment = 5; tl.second_wire = 2; tl.second_segment = 7;
    tl.propagation_constant_per_m = {0.1, 9.0}; tl.length_m = 0.12; tl.second_current_sign = -1;
    oc.transmission_line_loads = {tl};
    auto orr = solve_oriented_thin_wires(oc);
    std::vector<double> po{orr.accepted_power_w, orr.conductor_loss_w, orr.dielectric_loss_w, orr.lumped_load_loss_w, orr.transmission_line_loss_w, orr.radiated_power_w, orr.radiation_efficiency, orr.feed_current_a[0]};
    cv(po, orr.feed_impedance_ohm[0]);
    vec("oriented", po); vec("oriented_I", cvec(orr.current_a));
    // Two-port coupled dipoles.
    OrientedWireMomConfig mp; mp.frequency_hz = 300e6;
    OrientedThinWire a; a.start_m = {0, 0, -0.24}; a.end_m = {0, 0, 0.24}; a.radius_m = 2e-3; a.segments = 11;
    OrientedThinWire b = a; b.start_m.x = b.end_m.x = 0.25;
    mp.wires = {a, b}; mp.feeds = {{0}, {1}};
    vec("zport", cmat(oriented_wire_port_impedance_matrix(mp)));
    std::vector<double> ref{50.0, 75.0};
    vec("sport", cmat(oriented_wire_port_s_parameters(mp, ref)));
    ParallelWireMomConfig pc; pc.frequency_hz = 280e6; pc.wires = {{0, 0, 0, 0.5, 1e-3, 11}, {0.15, 0.0, 0.02, 0.45, 1e-3, 11}};
    pc.feeds = {{0}}; pc.loads = {{1, 5, Cx{0.0, -30.0}}};
    auto pr = solve_parallel_thin_wires(pc);
    vec("parallel_I", cvec(pr.current_a));
    // Antenna utilities.
    auto hf = hertzian_dipole_far_field(0.01, 2.0, 1e9, 0.7, 10.0);
    std::vector<double> an; cv(an, hf.e_theta); cv(an, hf.h_phi);
    an.push_back(shielding_effectiveness_db(3.0, 0.02));
    std::vector<ArrayElement> el{{0, 0, 0, {1, 0}}, {0.1, 0, 0, {0.5, 0.5}}, {0.2, 0.05, 0.01, {0.2, -0.7}}};
    cv(an, array_factor(el, 1.5e9, 0.6, 0.3));
    auto pm = polarization_metrics({1.0, 0.2}, {0.3, -0.8});
    an.insert(an.end(), {pm.axial_ratio, pm.tilt_angle_rad, pm.ellipticity_angle_rad, pm.right_hand ? 1.0 : 0.0});
    auto pl = polarization_metrics({1.0, 0.0}, {0.5, 0.0});
    an.push_back(pl.axial_ratio);
    auto mm = antenna_match_metrics({73.0, 42.5}, 50.0);
    cv(an, mm.reflection_coefficient); an.push_back(mm.return_loss_db); an.push_back(mm.vswr);
    auto m0 = antenna_match_metrics({50.0, 0.0}, 50.0); an.push_back(m0.return_loss_db); an.push_back(m0.vswr);
    std::vector<AntennaOptimizationVariable> vars{{0.0, -1.0, 2.0, 0.5}, {1.0, 0.0, 3.0, 0.25}};
    auto opt = optimize_antenna(vars, [](std::span<const double> x) { return (x[0] - 0.7) * (x[0] - 0.7) + 2.0 * (x[1] - 2.2) * (x[1] - 2.2) + 0.1 * x[0] * x[1]; });
    an.insert(an.end(), {opt.parameters[0], opt.parameters[1], opt.objective, (double)opt.iterations, (double)opt.evaluations, opt.converged ? 1.0 : 0.0});
    SinusoidalDipoleSolver sd(0.75, 300e6, 1.5);
    an.insert(an.end(), {sd.wavelength(), sd.electrical_length(), sd.pattern_factor(0.9), sd.radiated_power(), sd.radiation_resistance(), sd.directivity(), sd.radiated_power(100)});
    cv(an, sd.far_e_theta(1.1, 20.0));
    auto pat = sd.normalized_pattern(13); for (auto& s : pat) { an.push_back(s.theta_rad); an.push_back(s.normalized_power); }
    vec("antenna", an);
    // Surface MoM: 2x2 plate (8 triangles) with plane wave, and delta gap.
    std::vector<SurfaceTriangle> plate;
    const double h = 0.06;
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        cfd::fem::Point3 p00{i * h, j * h, 0}, p10{(i + 1) * h, j * h, 0}, p01{i * h, (j + 1) * h, 0}, p11{(i + 1) * h, (j + 1) * h, 0};
        plate.push_back({p00, p10, p11}); plate.push_back({p00, p11, p01});
    }
    SurfaceMomConfig sc; sc.frequency_hz = 1.2e9; sc.plane_wave_direction = {0.2, 0.0, -1.0}; sc.plane_wave_electric = {1.0, 0.3, 0.0};
    auto sr = solve_pec_surface_mom(plate, sc);
    std::vector<double> sm{(double)sr.basis.size(), sr.relative_residual, sr.wave_number_per_m, sr.wave_impedance_ohm};
    for (auto& bf : sr.basis) { sm.push_back(bf.edge_length_m); sm.push_back(bf.plus_sign); sm.push_back(bf.minus_sign); sm.push_back(bf.minus_triangle == std::numeric_limits<std::size_t>::max() ? -1.0 : (double)bf.minus_triangle); }
    vec("surface", sm); vec("surface_I", cvec(sr.current_a));
    auto ff = surface_mom_far_field(sr, {0.1, 0.3, 1.0}, 5.0);
    std::vector<double> fv; cv(fv, ff.e_x); cv(fv, ff.e_y); cv(fv, ff.e_z); fv.push_back(ff.power_density_w_m2);
    vec("surface_ff", fv);
    SurfaceMomConfig dg; dg.frequency_hz = 0.9e9; dg.excitation = SurfaceMomExcitation::delta_gap; dg.feed_edge = 4; dg.feed_voltage_v = {2.0, 0.5};
    auto dr = solve_pec_surface_mom(plate, dg);
    std::vector<double> dv; cv(dv, dr.feed_current_a); cv(dv, dr.feed_impedance_ohm); dv.push_back(dr.accepted_power_w); dv.push_back(dr.relative_residual);
    vec("surface_gap", dv);
    // Dielectric SIE on a tetrahedron surface.
    cfd::fem::Point3 A{0, 0, 0}, B{0.05, 0, 0}, Cc{0, 0.05, 0}, D{0, 0, 0.05};
    std::vector<SurfaceTriangle> tet{{A, Cc, B}, {A, B, D}, {A, D, Cc}, {B, Cc, D}};
    DielectricSIEConfig dc; dc.frequency_hz = 2e9; dc.relative_permittivity_inside = 6.0;
    auto ds = solve_dielectric_surface_sie(tet, dc, {0.3, 1.0, -0.5});
    auto dsv = cvec(ds.normal_polarization); dsv.push_back(ds.relative_residual); dsv.push_back(ds.converged);
    vec("sie", dsv);
    std::vector<cfd::fem::Point3> src, tgt; std::vector<Cx> q;
    for (int i = 0; i < 40; ++i) { src.push_back({std::sin(1.3 * i), std::cos(0.7 * i) * 0.5, 0.1 * i - 2.0}); q.push_back({std::cos(0.3 * i), 0.1 * i}); }
    for (int i = 0; i < 5; ++i) tgt.push_back({5.0 + i, -3.0 + 0.5 * i, 1.0 * i});
    tgt.push_back(src[3]);
    vec("helm_direct", cvec(helmholtz_direct_apply(tgt, src, q, 2.5)));
    HelmholtzTreecodeConfig tc; tc.leaf_size = 4; tc.opening_angle = 0.4;
    vec("helm_tree", cvec(helmholtz_treecode_apply(tgt, src, q, 2.5, tc)));
    // Physical optics on the plate and SBR in a dihedral corner.
    PhysicalOpticsConfig poc; poc.frequency_hz = 10e9; poc.incident_direction = {0.1, 0.0, -1.0}; poc.observation_direction = {-0.1, 0.0, 1.0}; poc.observation_distance_m = 100.0;
    auto por = solve_physical_optics(plate, poc);
    std::vector<double> pov{por.illuminated_area_m2, por.rcs_m2, por.wave_number_per_m, por.wave_impedance_ohm};
    for (auto z : por.far_field) cv(pov, z);
    for (auto& cur : por.currents) for (auto z : cur.current) cv(pov, z);
    poc.incident_direction = {0.0, 0.0, 1.0};
    pov.push_back(solve_physical_optics(plate, poc).rcs_m2);
    vec("po", pov);
    std::vector<SurfaceTriangle> corner{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}, {{0, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {{0, 0, 0}, {0, 1, 0}, {0, 1, 1}}, {{0, 0, 0}, {0, 1, 1}, {0, 0, 1}}};
    cfd::optics::Ray ray; ray.origin = {0.6, 0.4, 0.8}; ray.direction = {-1.0, 0.1, -1.2};
    SbrConfig sbc; sbc.frequency_hz = 3e9;
    auto sb = trace_shooting_bouncing_rays(corner, ray, sbc);
    std::vector<double> sv{(double)sb.hits.size(), sb.optical_path_m, sb.escaped ? 1.0 : 0.0, sb.escaped_ray.origin.x, sb.escaped_ray.origin.y, sb.escaped_ray.origin.z, sb.escaped_ray.direction.x, sb.escaped_ray.direction.y, sb.escaped_ray.direction.z};
    cv(sv, sb.phase);
    for (auto& hh : sb.hits) sv.insert(sv.end(), {(double)hh.triangle, hh.distance_m, hh.point.x, hh.point.y, hh.point.z, hh.normal.x, hh.normal.y, hh.normal.z});
    vec("sbr", sv);
    std::cout << "\n}\n";
}
