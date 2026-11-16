// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream FVM particle/spray/two-fluid units and prints JSON goldens.
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/solvers/fvm/dense_particle_rheology.hpp"
#include "cfd/solvers/fvm/euler_euler.hpp"
#include "cfd/solvers/fvm/particle_interactions.hpp"
#include "cfd/solvers/fvm/particles.hpp"
#include "cfd/solvers/fvm/spray.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd;
using fvm::Vec3;
static bool first = true;
static void vec(const std::string& k, const std::vector<double>& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << std::setprecision(17) << v[i];
    std::cout << "]";
}
static void push(std::vector<double>& out, Vec3 v) { out.push_back(v.x); out.push_back(v.y); out.push_back(v.z); }
static std::vector<double> flat(const std::vector<Vec3>& v) { std::vector<double> o; for (auto x : v) push(o, x); return o; }
static std::vector<double> flat(const std::vector<fvm::Particle>& v) {
    std::vector<double> o;
    for (auto& p : v) { push(o, p.position); push(o, p.velocity); o.push_back(p.diameter); o.push_back(p.density); }
    return o;
}

int main() {
    std::cout << "{\n";
    auto mesh = fvm::make_sheared_cartesian_hexa_mesh(5, 4, 1, 0.5, 0.4, 0.1, 0.0);
    std::vector<Vec3> u(mesh.cell_count());
    for (std::size_t c = 0; c < u.size(); ++c) { auto p = mesh.cells()[c].center; u[c] = {0.2 + p.y, -0.1 * p.x, 0.0}; }

    // Particle cloud, two-way coupled.
    fvm::ParticleCloud cloud;
    cloud.add({{0.05, 0.05, 0.05}, {0.0, 0.0, 0.0}, 2e-4, 1500.0});
    cloud.add({{0.31, 0.22, 0.05}, {0.5, 0.1, 0.0}, 1e-4, 1000.0});
    cloud.add({{0.44, 0.37, 0.05}, {-0.2, 0.0, 0.1}, 5e-5, 800.0});
    for (int i = 0; i < 4; ++i) cloud.advance(mesh, u, 1.8e-5, 0.01, true);
    vec("cloud", flat(cloud.particles()));
    vec("cloud_source", flat(cloud.fluid_momentum_source()));
    std::vector<double> forces;
    push(forces, fvm::saffman_lift_force({0.3, -0.1, 0.05}, {0.0, 0.0, 40.0}, 1.2, 1.5e-5, 1e-4));
    push(forces, fvm::virtual_mass_force({1.0, 2.0, 0.0}, {0.5, -1.0, 0.3}, 1000.0, 2e-3));
    vec("forces", forces);
    fvm::Particle a{{0.0, 0.0, 0.0}, {1.0, 0.2, 0.0}, 1e-3, 1000.0}, b{{8e-4, 3e-4, 0.0}, {-0.5, 0.0, 0.1}, 7e-4, 900.0};
    fvm::collide_particles_elastic(a, b, 0.8);
    auto merged = fvm::coalesce_particles(a, b);
    vec("collide_coalesce", flat(std::vector<fvm::Particle>{a, b, merged}));
    vec("breakup_evap", {fvm::breakup_by_weber(1.2, 30.0, 1e-3, 0.072) ? 1.0 : 0.0, fvm::breakup_by_weber(1.2, 5.0, 1e-3, 0.072) ? 1.0 : 0.0,
                         fvm::evaporate_d2_law(1e-4, 1e-9, 2.0), fvm::evaporate_d2_law(1e-4, 1e-8, 2.0)});
    vec("volume_fraction", fvm::particle_volume_fraction(mesh, cloud.particles()));

    // Interaction system: one coalescing pair, one bouncing/breaking pair, one isolated particle.
    std::vector<fvm::Particle> ps{
        {{0.0, 0.0, 0.0}, {0.05, 0.0, 0.0}, 1e-3, 1000.0},
        {{9e-4, 0.0, 0.0}, {-0.05, 0.0, 0.0}, 1e-3, 1000.0},
        {{0.01, 0.0, 0.0}, {15.0, 1.0, 0.0}, 1.2e-3, 1000.0},
        {{0.0109, 2e-4, 0.0}, {-15.0, 0.0, 0.0}, 1e-3, 1000.0},
        {{0.05, 0.05, 0.05}, {0.0, 0.0, 0.0}, 1e-3, 1000.0}};
    fvm::ParticleInteractionSystem sys;
    sys.resolve(ps);
    auto st = sys.last_stats();
    vec("interactions", flat(ps));
    vec("interaction_stats", {double(st.candidate_pairs), double(st.collisions), double(st.coalescences), double(st.breakups), double(st.removed_particles),
                              st.momentum_before, st.momentum_after, st.kinetic_energy_before, st.kinetic_energy_after});

    // Dense rheology.
    fvm::DenseParticleRheology dense(mesh);
    std::vector<fvm::Particle> grains;
    for (int i = 0; i < 40; ++i) grains.push_back({{0.02 + 0.011 * i, 0.03 + 0.008 * (i % 7), 0.05}, {0.1 * std::sin(i), 0.05, 0.0}, 4e-3, 2500.0});
    auto dr = dense.evaluate(grains, u);
    vec("dense_phi", dr.solid_fraction);
    vec("dense_pressure", dr.granular_pressure);
    vec("dense_viscosity", dr.effective_viscosity);
    vec("dense_particle_source", flat(dr.particle_momentum_source));

    // Spray.
    fvm::SprayConfig sc; sc.dt = 1e-3; sc.evaporation_constant = 2e-7;
    fvm::SprayInjectionEvaporation spray(mesh, sc);
    fvm::SprayInjection inj; inj.position = {0.02, 0.2, 0.05}; inj.velocity = {5.0, 0.0, 0.0}; inj.diameter = 5e-5; inj.temperature = 290.0; inj.mass_flow_rate = 1e-6; inj.duration = 4.5e-3;
    spray.set_injection(inj);
    spray.add_parcel({{0.3, 0.1, 0.05}, {0.0, 1.0, 0.0}, 2e-5, 4e-12, 310.0});
    std::vector<double> temp(mesh.cell_count());
    for (std::size_t c = 0; c < temp.size(); ++c) temp[c] = 350.0 + 100.0 * mesh.cells()[c].center.x;
    for (int i = 0; i < 8; ++i) spray.step(u, temp);
    std::vector<double> parcels;
    for (auto& p : spray.parcels()) { push(parcels, p.position); push(parcels, p.velocity); parcels.push_back(p.diameter); parcels.push_back(p.mass); parcels.push_back(p.temperature); }
    vec("spray_parcels", parcels);
    vec("spray_vapor", spray.vapor_mass_source());
    vec("spray_energy", spray.carrier_energy_source());
    vec("spray_momentum", flat(spray.carrier_momentum_source()));
    vec("spray_totals", {spray.total_liquid_mass(), spray.injected_mass(), spray.evaporated_mass(), spray.time(), double(spray.steps())});

    // Euler-Euler.
    fvm::EulerEulerConfig ec; ec.dt = 2e-3; ec.dispersed_density = 50.0; ec.drag_coefficient = 20.0; ec.quadratic_drag = true; ec.dispersed_body_acceleration = {0.0, 9.81, 0.0}; ec.primary_body_acceleration = {0.0, -0.5, 0.0};
    fvm::EulerEulerTransport ee(mesh, ec);
    ee.initialize([](Vec3 p) { return p.y < 0.2 ? 0.3 : 0.05; }, Vec3{0.1, 0.0, 0.0}, Vec3{0.0, 0.2, 0.0});
    for (int i = 0; i < 10; ++i) ee.step();
    vec("ee_alpha", ee.dispersed_fraction());
    vec("ee_primary", flat(ee.primary_velocity()));
    vec("ee_dispersed", flat(ee.dispersed_velocity()));
    std::vector<double> tot{ee.volume_fraction_integral(), ee.minimum_fraction(), ee.maximum_fraction(), ee.time()};
    push(tot, ee.total_momentum());
    push(tot, ee.integrated_drag_force());
    vec("ee_totals", tot);
    vec("ee_face_flux", ee.dispersed_face_flux());
    std::cout << "\n}\n";
}
