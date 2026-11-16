// Test oracle only. Exercises upstream multibody/{math, rigid_body, constraints, contact, system,
// advanced_dynamics, cfd_dem_coupling, distributed_dem}.
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/multibody/advanced_dynamics.hpp"
#include "cfd/multibody/cfd_dem_coupling.hpp"
#include "cfd/multibody/constraints.hpp"
#include "cfd/multibody/contact.hpp"
#include "cfd/multibody/distributed_dem.hpp"
#include "cfd/multibody/system.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace cfd::multibody;
static bool first = true;
static void vec(const std::string& k, const std::vector<double>& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) { double d = v[i]; std::cout << (i ? "," : ""); if (std::isnan(d)) std::cout << "\"nan\""; else if (std::isinf(d)) std::cout << (d > 0 ? "\"inf\"" : "\"-inf\""); else std::cout << std::setprecision(17) << d; }
    std::cout << "]";
}
static void v3(std::vector<double>& o, Vec3 v) { o.insert(o.end(), {v.x, v.y, v.z}); }
static void q4(std::vector<double>& o, Quaternion q) { o.insert(o.end(), {q.w, q.x, q.y, q.z}); }
static void st(std::vector<double>& o, const RigidBody& b) { v3(o, b.state.position); q4(o, b.state.orientation); v3(o, b.state.linear_velocity); v3(o, b.state.angular_velocity); }
int main() {
    std::cout << "{\n";
    std::vector<double> m;
    auto q = quaternion_from_axis_angle({1, 2, 0.5}, 0.8);
    q4(m, q); v3(m, rotate(q, {0.3, -1, 2})); v3(m, inverse_rotate(q, {0.3, -1, 2})); v3(m, any_perpendicular({0.9, 0.1, 0.2})); v3(m, any_perpendicular({0.1, 0.9, 0.2}));
    Quaternion qq = q; integrate_orientation(qq, {0.5, -1.0, 2.0}, 0.01); q4(m, qq);
    for (auto integ : {RigidBodyIntegrator::semi_implicit_euler, RigidBodyIntegrator::velocity_verlet}) {
        RigidBody b; b.mass = 2.0; b.inertia_diagonal = {1.0, 2.0, 3.0}; b.state.orientation = q; b.state.angular_velocity = {1, 2, 0.5}; b.state.linear_velocity = {0.1, 0, 0};
        for (int i = 0; i < 5; ++i) { b.add_force_at_point({0.0, 1.0, 0.5}, b.state.position + Vec3{0.2, 0, 0}); b.add_torque({0.1, 0, 0}); integrate_rigid_body(b, 0.01, integ); }
        b.apply_impulse({0.5, 0, 0}, b.state.position + Vec3{0, 0.3, 0});
        st(m, b);
    }
    vec("rigid", m);
    // Joint rows from a 5-body chain.
    std::vector<RigidBody> bodies(5);
    for (int i = 0; i < 5; ++i) { bodies[i].mass = 1.0 + i; bodies[i].inertia_diagonal = {0.1 + 0.02 * i, 0.2, 0.15}; bodies[i].state.position = {0.5 * i, 0.1 * i, 0.0}; bodies[i].state.orientation = quaternion_from_axis_angle({0, 0, 1}, 0.1 * i); bodies[i].state.linear_velocity = {0.1 * i, -0.2, 0.05 * i}; bodies[i].state.angular_velocity = {0.0, 0.1 * i, 0.3}; }
    bodies[0].fixed = true;
    std::vector<Joint> joints{make_distance_constraint(bodies, 0, 1, {0.1, 0, 0}, {0.4, 0.1, 0}), make_spherical_joint(bodies, 1, 2, {0.75, 0.15, 0}), make_revolute_joint(bodies, 2, 3, {1.25, 0.25, 0}, {0, 0, 1}),
                              make_prismatic_joint(bodies, 3, 4, {1.75, 0.35, 0}, {1, 0.2, 0}), make_fixed_joint(bodies, 1, 3, {1.0, 0.2, 0.0}), make_gear_constraint(bodies, 2, 4, {0, 0, 1}, {0, 1, 0}, 2.0),
                              make_linear_motor(bodies, 0, 4, {1, 0, 0}, 0.3, 5.0), make_angular_motor(bodies, 1, 2, {0, 0, 1}, 1.0)};
    std::vector<double> jv;
    std::vector<ConstraintRow> rows;
    for (auto& j : joints) { auto r = build_constraint_rows(bodies, j, 0.01, 0.2); for (auto& row : r) { v3(jv, row.linear_a); v3(jv, row.angular_a); v3(jv, row.linear_b); v3(jv, row.angular_b); jv.push_back(row.bias); jv.push_back(row.lower); jv.push_back(row.upper); } rows.insert(rows.end(), r.begin(), r.end()); }
    solve_constraint_rows(bodies, rows, 12);
    for (auto& r : rows) jv.push_back(r.impulse);
    for (auto& b : bodies) st(jv, b);
    vec("joints", jv);
    // Systems.
    for (auto method : {ContactMethod::smooth_penalty, ContactMethod::nonsmooth_impulse}) {
        MultibodyConfig mc; mc.contact_method = method; mc.dt = 2e-3; mc.integrator = method == ContactMethod::smooth_penalty ? RigidBodyIntegrator::velocity_verlet : RigidBodyIntegrator::semi_implicit_euler;
        mc.impulse.restitution = 0.3;
        RigidBodySystem sys(mc);
        for (int i = 0; i < 4; ++i) { RigidBody b; b.state.position = {0.3 * i, 0.45 + 0.15 * i, 0.05 * i}; b.state.linear_velocity = {0.2 - 0.1 * i, -0.5, 0.0}; b.state.angular_velocity = {0, 0, 0.5 * i}; b.inertia_diagonal = {0.05, 0.05, 0.05}; b.mass = 1.0 + 0.5 * i; auto id = sys.add_body(b); sys.add_sphere(id, 0.2); }
        sys.add_plane({{0, 1, 0}, 0.0});
        TriangleMeshShape mesh; mesh.vertices = {{-1, 0.1, -1}, {2, 0.1, -1}, {2, 0.4, 1}, {-1, 0.4, 1}}; mesh.triangles = {{0, 1, 2}, {0, 2, 3}};
        sys.add_triangle_mesh(mesh);
        sys.add_joint(make_distance_constraint(sys.bodies(), 0, 1, sys.bodies()[0].state.position, sys.bodies()[1].state.position));
        sys.run(150);
        std::vector<double> sv{sys.time(), sys.total_kinetic_energy()}; v3(sv, sys.linear_momentum());
        for (auto& b : sys.bodies()) st(sv, b);
        vec(method == ContactMethod::smooth_penalty ? "sys_penalty" : "sys_impulse", sv);
    }
    for (bool hist : {true, false}) {
        DemConfig dc; dc.dt = 5e-5; dc.history_dependent_tangential = hist; dc.contact.cohesion_force = 0.01;
        ExplicitDemSystem dem(dc);
        for (int i = 0; i < 6; ++i) dem.add_particle({0.021 * i, 0.012 + 0.005 * (i % 2), 0.002 * (i % 3)}, 0.011, 2500.0, {0.05 * (3 - i), -0.1, 0.0});
        dem.add_plane({{0, 1, 0}, 0.0});
        BondedParticleModel bm; bm.tensile_failure_force = 0.5; bm.shear_failure_force = 0.3;
        dem.add_bond(0, 1, bm); dem.add_bond(2, 3, bm, 0.02); dem.add_bond(4, 5, {});
        std::vector<double> dv;
        for (int s = 0; s < 6; ++s) { dem.run(40); auto bs = dem.last_bond_stats(); dv.insert(dv.end(), {(double)bs.active, (double)bs.damaged, (double)bs.broken_total, (double)dem.contact_history().states().size()}); }
        for (auto& b : dem.bodies()) st(dv, b);
        for (auto& b : dem.bonds()) { dv.push_back(b.damage); dv.push_back(b.broken); v3(dv, b.tangential_displacement); }
        vec(hist ? "dem_hist" : "dem_plain", dv);
    }
    // Broad phases and GJK/EPA.
    std::vector<RigidBody> cloud(30); std::vector<SphereShape> sph;
    for (int i = 0; i < 30; ++i) { cloud[i].state.position = {std::fmod(0.37 * i, 1.0), std::fmod(0.61 * i, 0.9), std::fmod(0.13 * i, 0.4)}; sph.push_back({(std::size_t)i, 0.06 + 0.002 * (i % 5)}); }
    sph.push_back({3, 0.05});
    std::vector<double> bp;
    auto sap = broad_phase_sweep_and_prune(cloud, sph); std::sort(sap.begin(), sap.end());
    for (auto [a, b] : sap) { bp.push_back(a); bp.push_back(b); } bp.push_back(-1);
    for (auto [a, b] : broad_phase_bvh(cloud, sph)) { bp.push_back(a); bp.push_back(b); } bp.push_back(-1);
    for (auto [a, b] : broad_phase_cell_linked(cloud, sph, 0.2)) { bp.push_back(a); bp.push_back(b); } bp.push_back(-1);
    for (auto& c : sphere_contacts(cloud, sph, sap)) { v3(bp, c.point); v3(bp, c.normal); bp.push_back(c.penetration); }
    vec("broad", bp);
    std::vector<RigidBody> hulls(2);
    hulls[0].state.position = {0, 0, 0}; hulls[1].state.position = {0.8, 0.3, 0.1}; hulls[1].state.orientation = quaternion_from_axis_angle({0.3, 0.2, 1}, 0.6);
    ConvexHullShape ha{0, {{-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}, {-0.5, 0.5, -0.5}, {-0.5, -0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}}};
    ConvexHullShape hb{1, {{0, 0, 0.6}, {0.5, 0, -0.3}, {-0.4, 0.4, -0.3}, {-0.3, -0.5, -0.3}}};
    std::vector<double> gv;
    auto g = convex_gjk_epa_contact(hulls, ha, hb);
    gv.insert(gv.end(), {(double)g.intersect, (double)g.gjk_iterations, (double)g.epa_iterations, g.contact.penetration}); v3(gv, g.contact.point); v3(gv, g.contact.normal);
    hulls[1].state.position = {2.5, 0, 0};
    auto g2 = convex_gjk_epa_contact(hulls, ha, hb);
    gv.insert(gv.end(), {(double)g2.intersect, (double)g2.gjk_iterations});
    auto bx = convex_aabb(hulls[1], hb); v3(gv, bx.minimum); v3(gv, bx.maximum);
    vec("gjk", gv);
    // Implicit integrators on a Duffing oscillator pair.
    std::vector<double> iv;
    auto duffing = [](std::span<const double> x, std::span<const double> v, double t, std::span<double> f) { f[0] = -4.0 * x[0] - 0.5 * x[0] * x[0] * x[0] + 0.3 * (x[1] - x[0]) - 0.05 * v[0] + std::sin(t); f[1] = -2.0 * x[1] - 0.3 * (x[1] - x[0]) - 0.02 * v[1]; };
    for (auto cfg : {ImplicitSecondOrderConfig::newmark(0.05), ImplicitSecondOrderConfig::generalized_alpha(0.05, 0.6), ImplicitSecondOrderConfig::hht(0.05, -0.2)}) {
        ImplicitGeneralizedIntegrator in({2.0, 0.1, 0.1, 1.0}, 2, duffing, cfg);
        in.initialize({0.5, -0.2}, {0.0, 0.3});
        in.run(20);
        iv.insert(iv.end(), {in.q()[0], in.q()[1], in.v()[0], in.v()[1], in.a()[0], in.a()[1], in.time(), (double)in.last_newton_iterations()});
    }
    vec("implicit", iv);
    ArticulatedSystem arm;
    ArticulatedLink l0; l0.parent_joint_offset = {0, 1, 0}; l0.joint_to_com_local = {0.5, 0, 0}; l0.body.mass = 2.0; l0.body.inertia_diagonal = {0.05, 0.2, 0.2}; l0.q = 0.3; l0.damping = 0.05;
    ArticulatedLink l1; l1.parent = 0; l1.parent_joint_offset = {0.5, 0, 0}; l1.joint_to_com_local = {0.4, 0, 0}; l1.body.mass = 1.0; l1.body.inertia_diagonal = {0.02, 0.1, 0.1}; l1.q = -0.2; l1.qd = 0.5;
    ArticulatedLink l2; l2.parent = 1; l2.joint_type = ReducedJointType::prismatic; l2.axis_parent = {1, 0, 0}; l2.parent_joint_offset = {0.4, 0, 0}; l2.joint_to_com_local = {0.1, 0, 0}; l2.body.mass = 0.5; l2.body.inertia_diagonal = {0.01, 0.01, 0.01}; l2.q = 0.05;
    arm.add_link(l0); arm.add_link(l1); arm.add_link(l2); arm.set_effort(1, 0.4);
    std::vector<double> av; for (double x : arm.mass_matrix()) av.push_back(x); for (double x : arm.generalized_forces()) av.push_back(x);
    arm.run(50, 0.002);
    for (auto& l : arm.links()) { av.insert(av.end(), {l.q, l.qd, l.qdd}); st(av, l.body); }
    av.push_back(arm.kinetic_energy());
    vec("articulated", av);
    FlexibleBodyInterface flex({{{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}}, {FlexibleMode{{{0, 0, 0.1}, {0, 0.1, 0.2}, {0.1, 0, 0}}, 2.0, 5.0, 0.1}, FlexibleMode{{{0.05, 0, 0}, {0, 0, 0}, {0, 0.05, 0.05}}, 1.0, 3.0, 0.0}});
    flex.set_modal_state({0.3, -0.2}, {0.1, 0.05});
    RigidBody frame; frame.state.position = {1, 2, 3}; frame.state.orientation = q; frame.state.linear_velocity = {0.1, 0, 0}; frame.state.angular_velocity = {0, 0, 1};
    std::vector<double> fv; for (auto p : flex.node_positions(frame)) v3(fv, p); for (auto p : flex.node_velocities(frame)) v3(fv, p);
    std::vector<Vec3> nf{{1, 0, 0}, {0, 2, 0}, {0, 0, -1}};
    auto ld = flex.apply_nodal_forces(frame, nf); v3(fv, ld.body_force); v3(fv, ld.body_torque); for (double x : ld.modal_force) fv.push_back(x); v3(fv, frame.torque);
    vec("flex", fv);
    // CFD/DEM coupling.
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(3, 2, 2, 0.3, 0.2, 0.2);
    std::vector<double> pr(mesh.cell_count()); std::vector<SymmetricStress3> ss(mesh.cell_count()); std::vector<cfd::fvm::Vec3> uf(mesh.cell_count()), gp(mesh.cell_count()), vo(mesh.cell_count());
    for (std::size_t c = 0; c < mesh.cell_count(); ++c) { auto x = mesh.cells()[c].center; pr[c] = 100.0 + 50.0 * x.x - 20.0 * x.y; ss[c] = {0.1 * c, 0.2, -0.1, 0.05 * c, 0.01, -0.02}; uf[c] = {0.5 + x.y, 0.1 * x.z, -0.2}; gp[c] = {50.0, -20.0, 1.0 * c}; vo[c] = {0.0, 0.3 * c, 1.0}; }
    RigidBody ball; ball.state.position = {0.15, 0.1, 0.1};
    std::vector<double> cv;
    for (auto n : {std::size_t(6), std::size_t(15)}) { auto quad = sphere_surface_quadrature(ball, 0.04, n); auto load = integrate_resolved_surface_load(mesh, pr, ss, quad, ball.state.position); v3(cv, load.force); v3(cv, load.torque); cv.push_back(quad.size()); }
    std::vector<RigidBody> parts(4); std::vector<SphereShape> ps;
    for (int i = 0; i < 4; ++i) { parts[i].state.position = {0.05 + 0.07 * i, 0.04 + 0.03 * i, 0.1}; parts[i].state.linear_velocity = {0.1 * i, 0, 0.05}; ps.push_back({(std::size_t)i, 0.01 + 0.003 * i}); }
    auto ur = unresolved_cfd_dem_coupling(mesh, parts, ps, uf, gp, vo);
    for (double x : ur.solid_volume_fraction) cv.push_back(x); for (double x : ur.void_fraction) cv.push_back(x);
    for (auto f : ur.particle_force) v3(cv, f); for (auto f : ur.fluid_momentum_source) v3(cv, f);
    UnresolvedCfdDemConfig uc; uc.include_saffman_lift = false; uc.dynamic_viscosity = 1e-5;
    auto ur2 = unresolved_cfd_dem_coupling(mesh, parts, ps, uf, {}, {}, uc); for (auto f : ur2.particle_force) v3(cv, f);
    vec("coupling", cv);
    // Distributed DEM (2 emulated ranks).
    DemSlabDecomposition dd{0.0, 0.1, 2};
    std::vector<DistributedDemParticle> all;
    for (int i = 0; i < 8; ++i) { RigidBody b; b.state.position = {0.0125 * i + 0.003, 0.01 * (i % 2), 0.0}; b.state.linear_velocity = {0.02 * (4 - i), 0.0, 0.01}; b.state.angular_velocity = {0, 0, 0.3 * i}; b.mass = 0.01; b.inertia_diagonal = {1e-6, 1e-6, 1e-6}; all.push_back(make_distributed_dem_particle(100 + i, b, {0, 0.0068}, dd.owner_rank(b.state.position.x))); }
    std::vector<double> ddv{dd.slab_width(), (double)dd.owner_rank(0.07), dd.slab_min(1), dd.slab_max(0)};
    auto plan = plan_dem_slab_exchange(all, dd, 0, 0.004);
    for (auto& r : plan.migrate_indices_by_rank) { for (auto i : r) ddv.push_back(i); ddv.push_back(-1); }
    for (auto& r : plan.ghost_indices_by_rank) { for (auto i : r) ddv.push_back(i); ddv.push_back(-1); }
    std::vector<DistributedDemParticle> own0, own1, ghost0, ghost1;
    for (auto& p : all) (p.owner_rank == 0 ? own0 : own1).push_back(p);
    for (auto& p : own1) if (p.state.position.x - p.radius - 0.004 < 0.05) ghost0.push_back(p);
    for (auto& p : own0) if (p.state.position.x + p.radius + 0.004 > 0.05) ghost1.push_back(p);
    std::vector<DistributedDemContactHistory> h0, h1;
    std::vector<DistributedDemBond> b0{make_distributed_dem_bond(own0[2], own1[0], {}, 0.0), make_distributed_dem_bond(own0[0], own0[1], {}, 0.0)}, b1;
    HertzMindlinContactModel hm; hm.normal_stiffness = 2e4; hm.tangential_stiffness = 1e4;
    for (int s = 0; s < 3; ++s) {
        auto c0 = evaluate_distributed_dem_contacts(own0, ghost0, h0, 0, 1e-4, hm);
        auto c1 = evaluate_distributed_dem_contacts(own1, ghost1, h1, 1, 1e-4, hm);
        auto k0 = evaluate_distributed_dem_bonds(own0, ghost0, b0, 0, 1e-4);
        std::vector<DistributedDemForceContribution> f0 = c0.local, f1 = c1.local;
        accumulate_distributed_dem_forces(f0, k0.local);
        accumulate_distributed_dem_forces(f1, c0.remote); accumulate_distributed_dem_forces(f1, k0.remote);
        accumulate_distributed_dem_forces(f0, c1.remote);
        ddv.insert(ddv.end(), {(double)c0.contacts, (double)c1.contacts, (double)k0.active_bonds, (double)c0.remote.size(), (double)h0.size()});
        integrate_distributed_dem_particles(own0, f0, {0, -9.81, 0}, 1e-4);
        integrate_distributed_dem_particles(own1, f1, {0, -9.81, 0}, 1e-4);
        ghost0.clear(); ghost1.clear();
        for (auto& p : own1) if (p.state.position.x - p.radius - 0.004 < 0.05) ghost0.push_back(p);
        for (auto& p : own0) if (p.state.position.x + p.radius + 0.004 > 0.05) ghost1.push_back(p);
    }
    for (auto& p : own0) { v3(ddv, p.state.position); v3(ddv, p.state.linear_velocity); v3(ddv, p.state.angular_velocity); q4(ddv, p.state.orientation); }
    for (auto& p : own1) { v3(ddv, p.state.position); v3(ddv, p.state.linear_velocity); v3(ddv, p.state.angular_velocity); }
    for (auto& h : h0) { ddv.push_back(h.particle_a); ddv.push_back(h.particle_b); v3(ddv, h.tangential_displacement); ddv.push_back(h.age); }
    for (auto& b : b0) { ddv.push_back(b.rest_length); v3(ddv, b.tangential_displacement); ddv.push_back(b.damage); }
    vec("distributed", ddv);
    std::cout << "\n}\n";
}
