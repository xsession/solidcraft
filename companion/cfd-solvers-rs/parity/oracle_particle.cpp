// Test oracle only. Exercises upstream particle/{electromagnetic, cross_sections, monte_carlo,
// plasma, wakefield, transport}.
#include "cfd/particle/cross_sections.hpp"
#include "cfd/particle/electromagnetic.hpp"
#include "cfd/particle/monte_carlo.hpp"
#include "cfd/particle/plasma.hpp"
#include "cfd/particle/transport.hpp"
#include "cfd/particle/wakefield.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace cfd::particle;
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
using G = std::vector<double>;
static void v3(G& g, Vec3 v) { g.push_back(v.x); g.push_back(v.y); g.push_back(v.z); }
static void cp(G& g, const ChargedParticle& p) { v3(g, p.position_m); v3(g, p.velocity_m_per_s); g.push_back(p.weight); }
static void p2(G& g, const PicParticle2D& p) { g.push_back(p.position_m.x); g.push_back(p.position_m.y); v3(g, p.velocity_m_per_s); g.push_back(p.weight); }
static void p3(G& g, const PicParticle3D& p) { v3(g, p.position_m); v3(g, p.velocity_m_per_s); g.push_back(p.weight); }
static void app(G& g, const std::vector<double>& v) { g.insert(g.end(), v.begin(), v.end()); }
static void cs(G& g, const CollisionStatistics& s) { g.push_back(s.elastic_events); g.push_back(s.ionization_events); g.push_back(s.energy_loss_j); }
static void fld(G& g, const ElectromagneticField& f) { v3(g, f.electric_v_per_m); v3(g, f.magnetic_t); }

int main() {
    std::cout << "{\n";
    const double q = -1.602176634e-19, me = 9.1093837015e-31;
    // pushers
    {
        G g; ElectromagneticField f{{1.0e3, -2.0e3, 5.0e2}, {0.01, 0.02, -0.03}};
        ChargedParticle a{{0.1, 0.2, 0.3}, {1.0e5, -2.0e5, 3.0e5}, q, me, 2.0};
        ChargedParticle b = a, c = a; b.velocity_m_per_s = {1.0e8, 5.0e7, -2.0e7}; c.velocity_m_per_s = b.velocity_m_per_s;
        for (int i = 0; i < 7; ++i) { boris_push(a, f, 1.0e-11); relativistic_boris_push(b, f, 1.0e-11); vay_push(c, f, 1.0e-11); }
        cp(g, a); cp(g, b); cp(g, c);
        g.push_back(lorentz_gamma(a)); g.push_back(lorentz_gamma(b)); g.push_back(kinetic_energy_j(a)); g.push_back(kinetic_energy_j(c));
        ChargedParticle t{{0, 0, 0}, {1.0e4, 0, 0}, 1.0, 2.0, 1.0};
        track_particle(t, [](Vec3 p, double time) { return ElectromagneticField{{std::sin(p.x) + time, 0.5, 0.0}, {0.0, 0.0, 1.0 + p.y}}; }, 0.25, 1.0e-3, 9);
        cp(g, t);
        vec("pushers", g);
    }
    // Monte-Carlo collisions, SEY, walls
    {
        G g; std::vector<ChargedParticle> ps;
        for (int i = 0; i < 6; ++i) ps.push_back({{0.1 * i, 0.0, 0.0}, {1.0e6 * (1 + i), 2.0e5 * i, -3.0e5}, q, me, 1.0 + i});
        NeutralCollisionModel m{1.0e22, 1.0e-19, 5.0e-20, 15.0, 42ULL};
        cs(g, apply_monte_carlo_collisions(ps, 1.0e-9, m)); cs(g, apply_monte_carlo_collisions(ps, 1.0e-9, m));
        g.push_back(ps.size()); for (auto& p : ps) cp(g, p);
        SecondaryEmissionModel s;
        for (double e : {0.0, 10.0, 20.0, 50.0, 300.0, 1000.0, 5000.0}) g.push_back(secondary_electron_yield(e, s));
        AxisAlignedParticleBox box{{0, 0, 0}, {1, 2, 3}};
        for (int mode = 0; mode < 2; ++mode)
            for (Vec3 pos : {Vec3{1.2, 0.5, 0.5}, Vec3{-0.1, 2.3, 1.0}, Vec3{0.5, 1.0, -0.4}, Vec3{0.5, 0.5, 0.5}}) {
                ChargedParticle p{pos, {2.0e6, -1.0e6, 5.0e5}, q, me, 3.0};
                auto r = apply_particle_box_boundary(p, box, mode ? ParticleWallMode::specular_reflect : ParticleWallMode::absorb, &s);
                cp(g, p); g.push_back(r.alive); g.push_back(r.impacted); v3(g, r.impact_normal); g.push_back(r.secondary_macro_weight);
            }
        vec("collisions_walls", g);
    }
    // 1-D PIC
    {
        G g; std::vector<double> rho(16); for (int i = 0; i < 16; ++i) rho[i] = std::sin(2 * M_PI * i / 16.0) + 0.1 * i;
        app(g, periodic_electric_field_from_charge_density(rho, 2.0, 2.5, true));
        std::vector<double> rho0(16); for (int i = 0; i < 16; ++i) rho0[i] = (i % 2 ? 1.0 : -1.0) * (1 + i / 2);
        app(g, periodic_electric_field_from_charge_density(rho0, 2.0, 2.5, false));
        ElectrostaticPic1DConfig c; c.grid_points = 32; c.length_m = 1.0; c.dt_s = 0.01;
        ElectrostaticPic1D es(c);
        es.set_particles({{0.1, 0.5, 1.0, 1.0, 1.0}, {0.45, -0.2, -1.0, 2.0, 1.5}, {0.8, 0.0, 1.0, 1.0, 0.5}, {1.3, 0.1, 0.5, 1.0, 1.0}});
        es.step(5);
        app(g, es.charge_density()); app(g, es.electric_field()); for (auto& p : es.particles()) { g.push_back(p.position_m); g.push_back(p.velocity_m_per_s); }
        std::vector<ChargedParticle> o, n;
        for (int i = 0; i < 4; ++i) { o.push_back({{0.1 + 0.2 * i, 0, 0}, {}, (i % 2 ? -1.0 : 1.0), 1.0, 1.0 + i}); n.push_back(o.back()); n.back().position_m.x += 0.03 * (i + 1); }
        auto cc = deposit_charge_conserving_current_1d(o, n, 16, 1.0, 0.1);
        app(g, cc.charge_density_old_c_per_m3); app(g, cc.charge_density_new_c_per_m3); app(g, cc.current_x_a_per_m2); g.push_back(cc.continuity_linf_residual);
        vec("pic1d_es", g);
    }
    {
        G g; ElectromagneticPic1DConfig c; c.grid_points = 24; c.dt_s = 1.0e-10;
        ElectromagneticPic1D em(c);
        em.initialize_right_traveling_mode(100.0, 2);
        std::vector<ChargedParticle> ps; for (int i = 0; i < 5; ++i) ps.push_back({{0.13 + 0.17 * i, 0, 0}, {1.0e6 * (i - 2), 2.0e5, -1.0e5 * i}, q * 1e10, me * 1e10, 1.0});
        em.set_particles(ps);
        NeutralCollisionModel m{1.0e20, 1.0e-19, 1.0e-20, 15.0, 7ULL};
        em.step(4, &m); em.step(2);
        auto d = em.diagnostics();
        g.push_back(d.field_energy_j); g.push_back(d.particle_kinetic_energy_j); g.push_back(d.total_energy_j); g.push_back(d.charge_continuity_linf_residual); cs(g, d.collision_statistics); g.push_back(d.particle_count);
        app(g, em.charge_density()); app(g, em.current_x()); app(g, em.current_y()); app(g, em.electric_x()); app(g, em.electric_y()); app(g, em.electric_z()); app(g, em.magnetic_y()); app(g, em.magnetic_z());
        for (auto& p : em.particles()) cp(g, p); g.push_back(em.time_s()); fld(g, em.gather_field(0.37));
        vec("pic1d_em", g);
    }
    // 2-D PIC
    std::vector<PicParticle2D> ps2;
    for (int i = 0; i < 6; ++i) ps2.push_back({{0.07 + 0.15 * i, 0.9 - 0.13 * i}, {2.0e5 * (i - 3), 1.0e5 * i, 3.0e4}, (i % 2 ? q : -q) * 1e9, me * 1e9, 1.0 + 0.5 * i});
    {
        G g; app(g, deposit_cic_charge_density_2d(ps2, 8, 6, 1.0, 1.0));
        std::vector<double> rho(48); for (int i = 0; i < 48; ++i) rho[i] = std::cos(0.3 * i) + 0.01 * i;
        auto f = periodic_electric_field_from_charge_density_2d(rho, 8, 6, 1.0, 1.5, 2.0, true); app(g, f.electric_x_v_per_m); app(g, f.electric_y_v_per_m);
        auto n2 = ps2; for (auto& p : n2) { p.position_m.x += 0.02; p.position_m.y -= 0.015; }
        auto cc = deposit_charge_conserving_current_2d(ps2, n2, 8, 6, 1.0, 1.0, 1.0e-9);
        app(g, cc.charge_density_new_c_per_m3); app(g, cc.current_x_a_per_m2); app(g, cc.current_y_a_per_m2); g.push_back(cc.continuity_linf_residual);
        ElectrostaticPic2DConfig c; c.nx = 8; c.ny = 8; c.dt_s = 1.0e-9; ElectrostaticPic2D es(c); es.set_particles(ps2); es.step(3);
        app(g, es.charge_density()); app(g, es.electric_x()); for (auto& p : es.particles()) p2(g, p);
        auto e = es.gather_electric_field({0.33, 0.71}); g.push_back(e.x); g.push_back(e.y);
        ElectromagneticPic2DConfig ec; ec.nx = 10; ec.ny = 8; ec.dt_s = 1.0e-11; ElectromagneticPic2D em(ec);
        em.initialize_tm_z_mode(50.0, 1, 1); em.set_particles(ps2);
        NeutralCollisionModel m{1.0e21, 1.0e-19, 1.0e-20, 15.0, 11ULL}; em.step(3, &m); em.step(1);
        auto d = em.diagnostics(); g.push_back(d.field_energy_j); g.push_back(d.particle_kinetic_energy_j); g.push_back(d.charge_continuity_linf_residual); cs(g, d.collision_statistics); g.push_back(d.particle_count);
        app(g, em.electric_z()); app(g, em.magnetic_x()); app(g, em.magnetic_y()); app(g, em.current_x()); app(g, em.current_z());
        for (auto& p : em.particles()) p2(g, p);
        vec("pic2d", g);
    }
    {
        G g; SecondaryEmissionModel s;
        PicBoundaryBox2D box{{0, 0}, {1, 1}};
        for (int mode = 0; mode < 2; ++mode) for (Vec2 pos : {Vec2{1.1, 0.5}, Vec2{-0.2, -0.3}, Vec2{0.5, 0.5}}) {
            PicParticle2D p{pos, {1.0e6, -2.0e6, 0}, q, me, 2.0};
            auto r = apply_particle_box_boundary_2d(p, box, mode ? ParticleWallMode::specular_reflect : ParticleWallMode::absorb, &s);
            p2(g, p); g.push_back(r.alive); g.push_back(r.impacted); g.push_back(r.impact_normal.x); g.push_back(r.impact_normal.y); g.push_back(r.secondary_macro_weight);
        }
        for (auto mode : {GridBoundaryMode2D::electric_wall, GridBoundaryMode2D::absorbing_sponge}) {
            std::vector<double> ex(30, 1.0), ey(30, -2.0); GridBoundary2DConfig b; b.mode = mode; b.sponge_cells = 2; b.sponge_strength = 2.5;
            apply_electrostatic_field_boundary_2d(ex, ey, 6, 5, b); app(g, ex); app(g, ey);
        }
        for (int variant = 0; variant < 2; ++variant) {
            StaggeredElectromagneticPic2DConfig c; c.nx = 10; c.ny = 9; c.dt_s = 5.0e-11;
            if (variant) { c.field_boundary.mode = GridBoundaryMode2D::absorbing_sponge; c.periodic_particles = false; c.particle_boundary = ParticleWallMode::absorb; c.enable_secondary_yield_report = true; }
            StaggeredElectromagneticPic2D st(c); st.initialize_tm_z_mode(30.0, 1, 2);
            auto ps = ps2; for (auto& p : ps) p.velocity_m_per_s.x *= 10.0 * (variant + 1);
            ps[0].position_m.x = 0.9995; ps[0].velocity_m_per_s.x = 6.0e7; ps[1].position_m.y = 0.0004; ps[1].velocity_m_per_s.y = -5.0e7;
            st.set_particles(ps); st.step(3);
            auto d = st.diagnostics(); g.push_back(d.field_energy_j); g.push_back(d.particle_kinetic_energy_j); g.push_back(d.charge_continuity_linf_residual); g.push_back(d.particle_count); g.push_back(d.absorbed_particles); g.push_back(d.secondary_macro_weight);
            app(g, st.electric_z()); app(g, st.magnetic_x()); app(g, st.magnetic_y()); app(g, st.current_x()); app(g, st.charge_density());
            for (auto& p : st.particles()) p2(g, p); fld(g, st.gather_field({0.41, 0.27}));
        }
        vec("pic2d_staggered", g);
    }
    // 3-D PIC
    std::vector<PicParticle3D> ps3;
    for (int i = 0; i < 7; ++i) ps3.push_back({{0.05 + 0.13 * i, 0.93 - 0.11 * i, 0.2 + 0.09 * i}, {1.5e5 * (i - 3), 1.0e5 * (i % 3), -2.0e4 * i}, (i % 2 ? q : -q) * 1e9, me * 1e9, 1.0 + 0.25 * i});
    {
        G g; app(g, deposit_cic_charge_density_3d(ps3, 5, 4, 6, 1.0, 1.2, 0.8));
        std::vector<double> rho(120); for (int i = 0; i < 120; ++i) rho[i] = std::sin(0.21 * i) + 0.02 * i;
        auto f = periodic_electric_field_from_charge_density_3d(rho, 5, 4, 6, 1.0, 1.2, 0.8, 1.5, true); app(g, f.electric_x_v_per_m); app(g, f.electric_y_v_per_m); app(g, f.electric_z_v_per_m);
        auto n3 = ps3; for (auto& p : n3) { p.position_m.x += 0.03; p.position_m.y -= 0.02; p.position_m.z += 0.01; }
        for (int local = 0; local < 2; ++local) {
            auto cc = local ? deposit_charge_conserving_current_3d_local(ps3, n3, 5, 4, 6, 1.0, 1.2, 0.8, 1.0e-9) : deposit_charge_conserving_current_3d(ps3, n3, 5, 4, 6, 1.0, 1.2, 0.8, 1.0e-9);
            app(g, cc.charge_density_old_c_per_m3); app(g, cc.charge_density_new_c_per_m3); app(g, cc.current_x_a_per_m2); app(g, cc.current_y_a_per_m2); app(g, cc.current_z_a_per_m2); g.push_back(cc.continuity_linf_residual);
        }
        ElectrostaticPic3DConfig c; c.nx = 6; c.ny = 5; c.nz = 4; c.dt_s = 1.0e-9; ElectrostaticPic3D es(c); es.set_particles(ps3); es.step(3);
        app(g, es.charge_density()); app(g, es.electric_x()); app(g, es.electric_z()); for (auto& p : es.particles()) p3(g, p);
        auto e = es.gather_electric_field({0.33, 0.71, 0.52}); v3(g, e);
        vec("pic3d", g);
    }
    {
        G g; ElectromagneticPic3DConfig c; c.nx = 6; c.ny = 5; c.nz = 4; c.dt_s = 1.0e-11; c.solve_longitudinal_poisson = false;
        for (int variant = 0; variant < 2; ++variant) {
            c.solve_longitudinal_poisson = variant == 1;
            ElectromagneticPic3D em(c); em.initialize_z_polarized_mode(40.0, 1, 1, 0); em.set_particles(ps3);
            NeutralCollisionModel m{1.0e21, 1.0e-19, 1.0e-20, 15.0, 13ULL}; em.step(2, &m); em.step(1);
            auto d = em.diagnostics(); g.push_back(d.field_energy_j); g.push_back(d.particle_kinetic_energy_j); g.push_back(d.charge_continuity_linf_residual); cs(g, d.collision_statistics); g.push_back(d.particle_count);
            app(g, em.electric_x()); app(g, em.electric_z()); app(g, em.magnetic_x()); app(g, em.magnetic_y()); app(g, em.current_z());
            for (auto& p : em.particles()) p3(g, p); fld(g, em.gather_field({0.41, 0.27, 0.66}));
        }
        vec("pic3d_em", g);
    }
    {
        G g;
        for (auto mode : {GridBoundaryMode3D::electric_wall, GridBoundaryMode3D::absorbing_sponge}) {
            std::vector<double> ex(60, 1.0), ey(60, -2.0), ez(60, 3.0), bx(60, 0.5), by(60, 0.25), bz(60, -1.0);
            GridBoundary3DConfig b; b.mode = mode; b.sponge_cells = 2; b.sponge_strength = 2.0; b.sponge_polynomial_order = 1.5;
            apply_electromagnetic_field_boundary_3d(ex, ey, ez, bx, by, bz, 5, 4, 3, b); app(g, ex); app(g, ez); app(g, bx); app(g, bz);
        }
        for (int variant = 0; variant < 3; ++variant) {
            StaggeredElectromagneticPic3DConfig c; c.nx = 6; c.ny = 5; c.nz = 5; c.dt_s = 5.0e-11;
            if (variant == 1) { c.field_boundary.mode = GridBoundaryMode3D::electric_wall; c.periodic_particles = false; c.particle_boundary = ParticleWallMode::specular_reflect; c.current_deposition = CurrentDeposition3DMode::local_finite_volume; c.sort_particles_by_cell = true; c.particle_sort_interval = 2; }
            if (variant == 2) { c.field_boundary.mode = GridBoundaryMode3D::absorbing_sponge; c.periodic_particles = false; c.particle_boundary = ParticleWallMode::absorb; c.enable_secondary_yield_report = true; c.sort_particles_by_cell = true; }
            StaggeredElectromagneticPic3D st(c); st.initialize_z_polarized_mode(25.0, 1, 0, 1);
            auto ps = ps3; for (auto& p : ps) p.velocity_m_per_s.x *= 10.0 * (variant + 1);
            ps[0].position_m.x = 0.9995; ps[0].velocity_m_per_s.x = 6.0e7; ps[1].position_m.z = 0.0004; ps[1].velocity_m_per_s.z = -5.0e7;
            st.set_particles(ps); st.step(4);
            auto d = st.diagnostics(); g.push_back(d.field_energy_j); g.push_back(d.particle_kinetic_energy_j); g.push_back(d.charge_continuity_linf_residual); g.push_back(d.particle_count); g.push_back(d.absorbed_particles); g.push_back(d.secondary_macro_weight); g.push_back(d.particle_sort_passes); g.push_back(d.occupied_particle_cells);
            app(g, st.electric_x()); app(g, st.electric_z()); app(g, st.magnetic_x()); app(g, st.magnetic_y()); app(g, st.current_x()); app(g, st.charge_density());
            for (auto& p : st.particles()) p3(g, p); fld(g, st.gather_field({0.41, 0.27, 0.96}));
        }
        vec("pic3d_staggered", g);
    }
    // 3-D domain decomposition and exchange
    {
        G g; PicDomainGrid3DConfig c; c.global_nx = 6; c.global_ny = 4; c.global_nz = 4; c.domains_x = 3; c.domains_y = 2; c.domains_z = 1; c.length_x_m = 1.2; c.length_y_m = 0.8; c.length_z_m = 0.6;
        auto doms = make_pic_domain_grid_3d(c);
        for (auto& d : doms) { g.push_back(d.domain_id); g.push_back(d.first_x); g.push_back(d.cells_x); g.push_back(d.first_y); g.push_back(d.cells_y); v3(g, d.minimum_m); v3(g, d.maximum_m); }
        auto sorted = sort_particles_by_cell_3d(ps3, 3, 2, 2, 1.0, 1.0, 1.0);
        for (auto i : sorted.original_indices) g.push_back(i); for (auto i : sorted.cell_offsets) g.push_back(i); g.push_back(sorted.occupied_cells);
        ParticleGuardHalo3DConfig hc{4, 4, 4, 1.0, 1.0, 1.0, 1, false};
        auto h = classify_particle_guard_halos_3d(ps3, hc);
        for (auto* v : {&h.x_minus, &h.x_plus, &h.y_minus, &h.y_plus, &h.z_minus, &h.z_plus, &h.interior, &h.outside_domain}) { g.push_back(v->size()); for (auto i : *v) g.push_back(i); }
        std::vector<std::vector<PicParticle3D>> by(doms.size());
        for (std::size_t i = 0; i < ps3.size(); ++i) { auto p = ps3[i]; p.position_m.x *= 1.2; p.position_m.y *= 0.8; p.position_m.z *= 0.6; by[i % doms.size()].push_back(p); }
        std::vector<std::vector<double>> cells(doms.size());
        for (std::size_t d = 0; d < doms.size(); ++d) for (std::size_t k = 0; k < doms[d].cells_x * doms[d].cells_y * doms[d].cells_z; ++k) cells[d].push_back(100.0 * d + k);
        auto blocks = exchange_scalar_guard_cells_3d(cells, doms, c, 1, -1.0);
        for (auto& b : blocks) app(g, b.values);
        auto topo = make_pic_rank_topology_3d(doms, 2, 1, 1);
        for (auto r : topo.domain_to_rank) g.push_back(r);
        auto round = run_serialized_distributed_pic_exchange_round_3d(by, cells, doms, c, 1, -2.0, topo);
        for (auto& v : round.particles_by_domain) { g.push_back(v.size()); for (auto& p : v) p3(g, p); }
        for (auto& b : round.scalar_guarded_blocks) app(g, b.values);
        auto& td = round.transport_diagnostics; for (double x : {double(td.envelopes), double(td.particle_messages), double(td.scalar_guard_messages), double(td.particle_payload_count), double(td.scalar_guard_value_count), double(td.same_rank_messages), double(td.remote_rank_messages)}) g.push_back(x);
        g.push_back(round.serialized_message_count); g.push_back(round.serialized_byte_count); for (auto b : round.serialized_plan.byte_count_by_destination_rank) g.push_back(b);
        ElectromagneticFieldBlocks3D fb;
        for (std::size_t d = 0; d < doms.size(); ++d) {
            const std::size_t n = doms[d].cells_x * doms[d].cells_y * doms[d].cells_z; std::vector<double> v(n);
            for (std::size_t k = 0; k < n; ++k) v[k] = std::sin(0.3 * k + d);
            fb.electric_x_by_domain.push_back(v); fb.electric_y_by_domain.push_back(std::vector<double>(n, 10.0)); fb.electric_z_by_domain.push_back(std::vector<double>(n, -5.0));
            fb.magnetic_x_by_domain.push_back(std::vector<double>(n, 1.0e-3)); fb.magnetic_y_by_domain.push_back(std::vector<double>(n, 0.0)); fb.magnetic_z_by_domain.push_back(v);
        }
        DistributedStaggeredPicStep3DConfig sc; sc.dt_s = 1.0e-9; sc.exterior_value = 0.5;
        auto stp = run_serialized_distributed_staggered_pic_step_3d(by, fb, doms, c, sc, topo);
        for (auto& v : stp.particles_by_domain) { g.push_back(v.size()); for (auto& p : v) p3(g, p); }
        for (auto& b : stp.guarded_fields.electric_x) app(g, b.values);
        g.push_back(stp.particle_count_before); g.push_back(stp.particle_count_after_push); g.push_back(stp.particle_count_after_migration); g.push_back(stp.max_particle_displacement_m);
        g.push_back(locate_pic_domain_3d({0.9, 0.5, 0.1}, c));
        vec("pic3d_domain", g);
    }
    // nuclear / Monte Carlo slab / plasma / wakefield
    {
        G g; CrossSectionTable t({{1.0, 2.0}, {10.0, 1.0}, {100.0, 0.5}});
        for (double e : {0.5, 1.0, 3.0, 10.0, 55.0, 1000.0}) g.push_back(t.value(e));
        MultiSpeciesDepletion dep({1.0e20, 0.0, 5.0e18});
        dep.add_reaction({0, 1, ParticleKind::neutron, t, 0.9}); dep.add_reaction({1, 2, ParticleKind::photon, CrossSectionTable({{1.0, 3.0e-24}}), 1.0});
        DepletionReaction r0; r0.parent = 0; r0.daughter = 2; r0.cross_section = CrossSectionTable({{1.0, 1.0e-24}, {5.0, 2.0e-24}}); dep.add_reaction(r0);
        for (int i = 0; i < 3; ++i) { auto d = dep.advance(10.0, 1.0e18, 2.0e17, 4.0); g.push_back(d.total_inventory); g.push_back(d.absorbed_rate); g.push_back(d.photon_rate); g.push_back(d.neutron_rate); }
        for (double v : dep.inventory()) g.push_back(v);
        SlabTransport1D slab({2.0, 1.5, 0.4});
        for (std::uint64_t seed : {1ULL, 99ULL}) { auto ta = slab.run(500, seed); g.push_back(ta.absorbed); g.push_back(ta.transmitted); g.push_back(ta.reflected); g.push_back(ta.total_path); g.push_back(ta.absorption_probability()); }
        vec("nuclear", g);
    }
    {
        G g; PlasmaState st; st.species = {{"e", -1.602176634e-19, me}, {"Ar", 0.0, 6.6e-26}, {"Ar+", 1.602176634e-19, 6.6e-26}}; st.number_density_m3 = {1.0e16, 1.0e22, 1.0e16}; st.electron_temperature_ev = 3.0;
        g.push_back(plasma_charge_density_c_m3(st));
        std::vector<PlasmaReaction> rx{{{{0, 1.0}, {1, 1.0}}, {{0, 2.0}, {2, 1.0}}, 1.0e-20}, {{{0, 1.0}, {2, 1.0}}, {{1, 1.0}}, 2.0e-19}};
        auto d = advance_plasma_chemistry(st, rx, 1.0e-6, 4);
        app(g, d.reaction_rates_m3_s); g.push_back(d.charge_density_before_c_m3); g.push_back(d.charge_density_after_c_m3); g.push_back(d.total_number_density_before_m3); g.push_back(d.total_number_density_after_m3); g.push_back(d.max_relative_density_change); app(g, st.number_density_m3);
        PaschenGas gas{11.0, 180.0, 0.01};
        for (double pd : {100.0, 1000.0, 1.0e5}) g.push_back(paschen_breakdown_voltage_v(gas, pd, 0.01));
        auto b = estimate_gas_breakdown_threshold(gas, 1.0e4, 0.002, 3.0e6); g.push_back(b.breakdown_voltage_v); g.push_back(b.breakdown_field_v_m); g.push_back(b.applied_field_exceeds_threshold);
        for (std::size_t order : {1U, 3U}) { auto m = estimate_parallel_plate_multipactor(1.0e9, 0.002, order, me, -1.602176634e-19, 1.3); g.push_back(m.order); g.push_back(m.transit_time_s); g.push_back(m.resonant_field_v_m); g.push_back(m.impact_energy_ev); g.push_back(m.secondary_yield_sustains); }
        std::vector<double> w, lam;
        for (int i = 0; i < 40; ++i) { w.push_back(resonator_longitudinal_wake_v_per_c(0.001 * i, 1.0e4, 1.3e9, 50.0)); lam.push_back(std::exp(-0.5 * std::pow((i - 20) / 4.0, 2))); }
        app(g, w); app(g, bunch_wake_potential_v(w, lam, 0.001)); auto z = wake_impedance_ohm(w, 0.001, 1.0e9); g.push_back(z.real()); g.push_back(z.imag());
        vec("plasma_wake", g);
    }
    // transport
    {
        G g; TransportWorld w;
        w.particles = {{"e-", -1.6e-19, me}, {"gamma", 0.0, 1.0e-40}};
        w.materials = {{"water", 1000.0, 1.0e3}, {"lead", 11340.0, 0.0}};
        w.sensitive_detectors = {{"det", 0.0, false}, {"strict", 5.0e4, true}};
        TransportRegion r0; r0.name = "slab0"; r0.box = {{0, 0, 0}, {0.1, 1, 1}}; r0.material_index = 0; r0.sensitive_detector_index = 0;
        TransportRegion r1; r1.name = "slab1"; r1.box = {{0.1, 0, 0}, {0.25, 1, 1}}; r1.material_index = 1; r1.sensitive_detector_index = 1;
        TransportRegion r2; r2.name = "slab2"; r2.box = {{0.25, 0, 0}, {0.4, 1, 1}}; r2.material_index = 0;
        w.regions = {r0, r1, r2};
        TransportProcess cont; cont.name = "ion"; cont.stopping_power_ev_per_m = 2.0e6; cont.particle_index = 0;
        TransportProcess disc; disc.name = "brems"; disc.kind = TransportProcessKind::discrete_interaction; disc.physical_interaction_length_m = 0.02; disc.energy_loss_fraction = 0.1; disc.secondary_energy_fraction = 0.2; disc.material_index = 1;
        TransportProcess abs; abs.name = "abs"; abs.kind = TransportProcessKind::absorber; abs.absorption_threshold_ev = 2.0e4;
        TransportPhysicsList pl; pl.name = "std"; pl.processes = {cont, disc, abs}; pl.production_cut_energy_ev = 500.0;
        auto world = apply_physics_list(w, pl);
        TransportConfig cfg; cfg.max_step_m = 0.03;
        TransportTrack t; t.position_m = {0.0, 0.5, 0.5}; t.direction = {1.0, 0.2, 0.0}; t.kinetic_energy_ev = 6.0e5;
        auto run = transport_track(world, t, cfg, 200);
        auto dumpres = [&](const TransportRunResult& r) {
            g.push_back(r.steps.size()); g.push_back(r.secondaries.size()); v3(g, r.primary.position_m); g.push_back(r.primary.kinetic_energy_ev); g.push_back(r.primary.time_s); g.push_back(r.primary.track_length_m); g.push_back(int(r.primary.status));
            g.push_back(r.scoring.total_energy_deposit_ev); g.push_back(r.scoring.total_track_length_m); g.push_back(r.scoring.steps); g.push_back(r.scoring.secondaries); g.push_back(r.scoring.boundary_crossings);
            for (auto& s : r.steps) { g.push_back(s.step_length_m); g.push_back(s.energy_deposit_ev); g.push_back(int(s.limiter)); g.push_back(s.region_index); g.push_back(s.limiting_process_index == std::size_t(-1) ? -1.0 : double(s.limiting_process_index)); }
            for (auto& s : r.secondaries) { v3(g, s.track.position_m); g.push_back(s.track.kinetic_energy_ev); g.push_back(s.track.time_s); }
        };
        dumpres(run);
        TransportRandom rng; rng.state = 12345;
        auto sr = transport_track_stochastic(world, t, cfg, rng, 200); dumpres(sr); g.push_back(double(rng.state >> 11));
        TransportRandom r2n; for (int i = 0; i < 5; ++i) g.push_back(sample_exponential_interaction_length_m(0.5, r2n));
        auto hits = collect_transport_hits(world, run.steps);
        g.push_back(hits.hits.size()); g.push_back(hits.total_energy_deposit_ev); g.push_back(hits.total_weighted_energy_deposit_ev);
        for (auto& h : hits.hits) { g.push_back(h.detector_index); v3(g, h.position_m); g.push_back(h.energy_deposit_ev); }
        auto grid = make_transport_dose_grid(4, 3, 2, {{0, 0, 0}, {0.4, 1, 1}}, 1000.0);
        score_hits_to_dose_grid(hits, grid); app(g, grid.energy_deposit_ev); app(g, grid.dose_gy);
        app(g, dose_rate_to_sar_w_per_kg(grid, 2.0)); app(g, project_dose_grid_to_pennes2d_sar(grid, 2.0, 3, 2));
        g.push_back(distance_to_box_boundary(r1.box, {0.15, 0.5, 0.5}, {-1, 1, 0}));
        TransportWorld big = world; big.regions.clear();
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 3; ++j) { TransportRegion r; r.name = "c"; r.box = {{0.25 * i, 0.3 * j, 0}, {0.25 * (i + 1), 0.3 * (j + 1), 1}}; r.material_index = (i + j) % 2; big.regions.push_back(r); }
        auto bvh = build_region_bvh(big, 2);
        for (Vec3 p : {Vec3{0.1, 0.1, 0.5}, Vec3{0.6, 0.45, 0.5}, Vec3{0.99, 0.89, 0.2}, Vec3{1.5, 0.1, 0.1}}) { auto a = locate_region(big, p, 0.0), b = locate_region_bvh(big, bvh, p, 0.0); g.push_back(a == std::size_t(-1) ? -1.0 : double(a)); g.push_back(b == std::size_t(-1) ? -1.0 : double(b)); }
        g.push_back(bvh.nodes.size());
        vec("transport", g);
    }
    std::cout << "\n}\n";
}
