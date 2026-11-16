// Test oracle only. Exercises upstream mechanics/*.
#include "cfd/mechanics/constitutive3d.hpp"
#include "cfd/mechanics/continuum_baselines.hpp"
#include "cfd/mechanics/cosserat_rods.hpp"
#include "cfd/mechanics/finite_strain_j2.hpp"
#include "cfd/mechanics/finite_strain_j2_calibration.hpp"
#include "cfd/mechanics/finite_strain_tet4.hpp"
#include "cfd/mechanics/finite_strain_tet4_plastic.hpp"
#include "cfd/mechanics/finite_strain_tri3_plastic.hpp"
#include "cfd/mechanics/meshfree_neighbor_graph.hpp"
#include "cfd/mechanics/phase_field_fracture1d.hpp"
#include "cfd/mechanics/phase_field_fracture_operator1d.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace cfd::mechanics;
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) { double d = static_cast<double>(x); std::cout << (i++ ? "," : ""); if (std::isnan(d)) std::cout << "\"nan\""; else if (std::isinf(d)) std::cout << (d > 0 ? "\"inf\"" : "\"-inf\""); else std::cout << std::setprecision(17) << d; }
    std::cout << "]";
}
template <class A> static void app(std::vector<double>& o, const A& a) { for (auto x : a) o.push_back(static_cast<double>(x)); }
int main() {
    std::cout << "{\n";
    const Tensor3 F{1.1, 0.05, -0.02, 0.03, 0.95, 0.04, -0.01, 0.02, 1.07};
    const Tensor3 D{0.2, -0.1, 0.3, 0.05, 0.4, -0.2, 0.1, 0.0, -0.3};
    CompressibleNeoHookean3D nh({2.0, 3.0});
    auto nr = nh.evaluate(F); std::vector<double> nv{nr.determinant, nr.energy_density}; app(nv, nr.first_piola); app(nv, nh.tangent_action(F, D));
    vec("neo", nv);
    J2Plasticity3D j2({200.0, 0.28, 0.8, 15.0, 1e-7});
    std::vector<double> jv;
    J2Plasticity3DState st{};
    for (Voigt6 e : {Voigt6{1e-3, -2e-4, 0, 1e-4, 0, 0}, Voigt6{8e-3, -2e-3, -1e-3, 2e-3, -1e-3, 5e-4}, Voigt6{1.2e-2, -3e-3, -2e-3, 3e-3, -2e-3, 1e-3}}) {
        auto r = j2.update(e, st); st = r.state;
        app(jv, r.stress); app(jv, r.algorithmic_tangent); app(jv, r.state.plastic_strain); jv.insert(jv.end(), {r.state.equivalent_plastic_strain, r.trial_equivalent_stress, r.yield_function, (double)r.yielded});
    }
    vec("j2", jv);
    FiniteStrainJ2 fj({200.0, 0.3, 1.0, 5.0, 1e-7, 32});
    std::vector<double> fv;
    FiniteStrainJ2State fs{};
    const Tensor3 F2{1.03, 0.04, 0.0, 0.01, 0.99, 0.02, 0.0, -0.01, 1.0};
    for (const Tensor3& g : {F2, F}) {
        auto r = fj.evaluate(g, fs); fs = r.state;
        fv.insert(fv.end(), {r.determinant, r.energy_density, r.trial_equivalent_stress, r.yield_function, (double)r.yielded, r.state.equivalent_plastic_strain});
        app(fv, r.cauchy_stress); app(fv, r.first_piola); app(fv, r.state.plastic_deformation_gradient);
    }
    app(fv, fj.tangent_action(F, D, fs)); app(fv, fj.central_tangent_action(F, D, fs));
    vec("fsj2", fv);
    auto cal = calibrate_finite_strain_j2(fj, F, D, {});
    vec("calib", std::vector<double>{cal.determinant, cal.rotation_relative_error, cal.hydrostatic_deviatoric_ratio, cal.tangent_step_relative_change, (double)cal.positive_state, (double)cal.objective_rotation, (double)cal.hydrostatic_response, (double)cal.tangent_step_converged});
    std::array<cfd::fem::Point3, 4> ref{{{0, 0, 0}, {1, 0, 0.1}, {0.1, 1, 0}, {0, 0.2, 1}}};
    std::array<cfd::fem::Point3, 4> cur{{{0.01, 0, 0}, {1.05, 0.02, 0.1}, {0.12, 0.98, 0.01}, {0.0, 0.21, 1.04}}};
    std::array<double, 12> dir{0.1, 0, 0, -0.1, 0.2, 0, 0.05, 0, 0.3, 0, -0.2, 0.1};
    auto t4 = assemble_total_lagrangian_tet4(ref, cur, dir, nh);
    std::vector<double> tv{t4.reference_volume, t4.determinant, t4.strain_energy}; app(tv, t4.internal_force); app(tv, t4.tangent_action);
    std::array<cfd::fem::Point3, 4> cur2{{{0.0, 0, 0}, {1.08, 0.02, 0.1}, {0.12, 0.96, 0.01}, {0.0, 0.21, 1.02}}};
    auto t4p = assemble_total_lagrangian_plastic_tet4(ref, cur2, dir, fj);
    tv.insert(tv.end(), {t4p.reference_volume, t4p.determinant, t4p.strain_energy, t4p.trial_equivalent_stress, t4p.yield_function, (double)t4p.yielded}); app(tv, t4p.internal_force); app(tv, t4p.tangent_action);
    std::array<cfd::fem::Node2, 3> r3{{{0, 0}, {1, 0.1}, {0.2, 0.9}}}, c3{{{0.0, 0.0}, {1.04, 0.12}, {0.21, 0.88}}};
    auto t3 = assemble_total_lagrangian_plastic_tri3(r3, c3, {0.1, 0.0, -0.1, 0.2, 0.05, 0.3}, fj, 0.5);
    tv.insert(tv.end(), {t3.reference_area, t3.determinant, t3.strain_energy, t3.trial_equivalent_stress, t3.yield_function, (double)t3.yielded}); app(tv, t3.internal_force); app(tv, t3.tangent_action); app(tv, t3.state.plastic_deformation_gradient);
    vec("elements", tv);
    NonlinearBarConfig bc; bc.elements = 6;
    auto bar = solve_nonlinear_bar(bc, 5.0);
    std::vector<double> bv = bar.displacement; bv.insert(bv.end(), {(double)bar.solve.iterations, (double)bar.solve.converged, bar.solve.residual_rms});
    vec("bar", bv);
    PeridynamicBarConfig pc; pc.points = 21; pc.horizon = 0.16;
    PeridynamicBar1D pd(pc);
    pd.set_velocity(20, 0.8); pd.set_displacement(0, -0.01);
    std::vector<double> ext(21, 0.0); ext[10] = 0.5;
    std::vector<double> pv;
    for (int i = 0; i < 40; ++i) { auto d = pd.step(0.01, i % 2 ? std::span<const double>{} : std::span<const double>{ext}); if (i % 8 == 7) pv.insert(pv.end(), {(double)d.broken_bonds, d.kinetic_energy, d.potential_energy}); }
    app(pv, pd.displacement()); app(pv, pd.velocity());
    vec("peri", pv);
    PenaltyImpact1D im(1.0, 2.0, 0.05, 1e4, 3.0);
    im.set_velocities(1.0, -0.5);
    std::vector<double> iv;
    for (int i = 0; i < 120; ++i) { im.step(1e-3); if (i % 20 == 19) { auto s = im.state(); iv.insert(iv.end(), {s.left_position, s.right_position, s.left_velocity, s.right_velocity, s.contact_force}); } }
    vec("impact", iv);
    std::vector<cfd::fem::Point3> pts;
    for (int i = 0; i < 60; ++i) pts.push_back({std::fmod(0.37 * i, 1.0), std::fmod(0.61 * i, 0.8), std::fmod(0.13 * i, 0.3)});
    MeshfreeNeighborGraphConfig mc; mc.support_radius = 0.22;
    auto g = build_meshfree_neighbor_graph(pts, mc);
    std::vector<double> gv{(double)g.cell_counts[0], (double)g.cell_counts[1], (double)g.cell_counts[2], g.cell_size};
    app(gv, g.offsets); app(gv, g.indices);
    std::vector<int> owner(60); for (int i = 0; i < 60; ++i) owner[i] = (i * 7) % 3;
    auto h = classify_meshfree_halo(g, owner, 1); app(gv, h.owned); gv.push_back(-1.0); app(gv, h.halo);
    mc.include_self = true; auto g2 = build_meshfree_neighbor_graph(pts, mc); gv.push_back((double)g2.indices.size());
    vec("meshfree", gv);
    PhaseFieldFracture1DConfig fc; fc.cells = 24; fc.notch_cell = 12; fc.initial_notch_damage = 0.3; fc.damage_iterations = 300;
    PhaseFieldFracture1D pf(fc);
    std::vector<double> ff;
    for (double u : {0.005, 0.01, 0.02, 0.015}) { auto d = pf.step(u); ff.insert(ff.end(), {d.reaction_force, d.elastic_energy, d.fracture_energy, d.total_energy, d.maximum_damage}); }
    app(ff, pf.damage()); app(ff, pf.history()); auto dd = pf.diagnostics(); ff.push_back(dd.total_energy);
    vec("pf", ff);
    CoupledPhaseFieldFracture1DConfig cc; cc.cells = 8;
    CoupledPhaseFieldFracture1D op(cc);
    std::vector<double> state(op.size()), dirv(op.size()), res(op.size(), 0.0), out(op.size(), 0.0);
    for (std::size_t i = 0; i < op.size(); ++i) { state[i] = i < 9 ? 0.002 * i * i : 0.05 * ((i * 3) % 7); dirv[i] = std::sin(0.7 * i); }
    op.update_history(state);
    auto er = op.evaluate(state, 0.18, res);
    op.tangent_action(state, dirv, 0.18, out);
    std::vector<double> ov{er.residual_l2, er.reaction_force, er.elastic_energy, er.fracture_energy, er.total_energy, er.maximum_damage};
    app(ov, res); app(ov, out); app(ov, op.history());
    vec("pfop", ov);
    CosseratRodConfig rc; rc.dt = 2e-3;
    StableCosseratRod rod(rc);
    std::vector<cfd::multibody::Vec3> rp; for (int i = 0; i < 6; ++i) rp.push_back({0.0, 0.0, 0.1 * i});
    std::vector<cfd::multibody::Quaternion> ro(6, cfd::multibody::Quaternion{});
    ro[3] = cfd::multibody::quaternion_from_axis_angle({1, 0, 0}, 0.1);
    rod.initialize(rp, ro);
    rod.set_fixed(0);
    std::vector<double> cv;
    for (int s = 0; s < 25; ++s) { rod.add_force(5, {0.5, 0.0, 0.0}); rod.add_torque(4, {0.0, 0.0, 0.2}); rod.step(); if (s % 5 == 4) { auto d = rod.diagnostics(); cv.insert(cv.end(), {d.max_stretch_shear_error, d.max_bend_twist_error, d.kinetic_energy, d.elastic_energy, d.linear_momentum_norm, rod.checksum()}); } }
    for (auto& n : rod.nodes()) cv.insert(cv.end(), {n.position.x, n.position.y, n.position.z, n.orientation.w, n.orientation.x, n.orientation.y, n.orientation.z, n.angular_velocity.x, n.angular_velocity.y, n.angular_velocity.z});
    cv.push_back(rod.time());
    vec("cosserat", cv);
    std::cout << "\n}\n";
}
