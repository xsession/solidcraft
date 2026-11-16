// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream FVM reacting/VOF units and prints JSON goldens.
#include "cfd/chemistry/kinetics.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/solvers/fvm/flame_models.hpp"
#include "cfd/solvers/fvm/premixed_flame1d.hpp"
#include "cfd/solvers/fvm/reacting.hpp"
#include "cfd/solvers/fvm/soot_radiation.hpp"
#include "cfd/solvers/fvm/species_transport.hpp"
#include "cfd/solvers/fvm/vof.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd;
static bool first = true;
static void vec(const std::string& k, const std::vector<double>& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << std::setprecision(17) << v[i];
    std::cout << "]";
}

int main() {
    std::cout << "{\n";
    auto mesh = fvm::make_sheared_cartesian_hexa_mesh(6, 4, 1, 0.6, 0.4, 0.1, 0.0);
    std::vector<double> flux(mesh.face_count());
    for (std::size_t f = 0; f < flux.size(); ++f) {
        const auto& face = mesh.faces()[f];
        flux[f] = face.area.x * 0.05 + face.area.y * 0.01;
    }
    // Multi-species transport.
    fvm::MultiSpeciesTransport species(mesh, {1e-3, 5e-3}, 0.01);
    species.initialize(0, [](fvm::Vec3 p) { return std::exp(-30.0 * ((p.x - 0.2) * (p.x - 0.2) + (p.y - 0.2) * (p.y - 0.2))); });
    species.initialize(1, [](fvm::Vec3 p) { return 1.0 + p.x; });
    species.set_face_flux(flux);
    species.run(5);
    vec("species0", species.values(0));
    vec("species1", species.values(1));
    vec("species_inventory", {species.inventory(0), species.inventory(1)});
    // Reacting cell sources: A + B -> C with Arrhenius kinetics.
    chemistry::ReactionNetwork net({{"A", 0.002, 0}, {"B", 0.032, 0}, {"C", 0.034, 0}});
    chemistry::ElementaryReaction r;
    r.reactants = {{0, 1.0}, {1, 1.0}};
    r.products = {{2, 1.0}};
    r.forward = {2.0e3, 0.5, 2.0e4};
    net.add_reaction(r);
    std::vector<std::vector<double>> cells{{1.0, 2.0, 0.1}, {0.5, 0.5, 0.0}};
    auto sources = fvm::reacting_cell_sources(net, cells, std::vector<double>{800.0, 1200.0});
    vec("reacting_sources", {sources[0][0], sources[0][1], sources[0][2], sources[1][0], sources[1][1], sources[1][2]});
    // Premixed flame.
    fvm::PremixedFlame1DConfig fc; fc.cells = 200; fc.length = 1.0; fc.diffusivity = 1e-3; fc.reaction_rate = 4.0; fc.dt = 0.002;
    fvm::PremixedFlame1D flame(fc);
    flame.initialize_front(0.2, 0.02);
    flame.step(500);
    vec("flame", {flame.front_location(), flame.theoretical_speed(), flame.temperature(30), flame.progress()[50], flame.progress()[120]});
    // Flame models and soot.
    fvm::PremixedFlameModel premixed({2.0, 5e5, 400.0});
    auto p1 = premixed.evaluate({0.3, 0.0, 500.0});
    auto p2 = premixed.evaluate({0.3, 0.0, 350.0});
    fvm::NonPremixedFlameModel nonpremixed({});
    auto n1 = nonpremixed.evaluate({0.2, 0.45, 600.0});
    vec("flame_models", {p1.progress_source, p1.heat_release, p2.progress_source, n1.progress_source, n1.mixture_fraction_source, n1.heat_release});
    fvm::SootRadiationModel soot({});
    auto s = soot.evaluate({0.001, 1500.0, 2.0e5}, 300.0);
    vec("soot", {s.soot_source, s.absorption_coefficient, s.radiative_source});
    // VOF.
    fvm::VofTransport vof(mesh, {0.02, 0.6});
    vof.initialize([](fvm::Vec3 p) { return p.x < 0.25 && p.y < 0.25 ? 1.0 : 0.0; });
    vof.set_face_flux(flux);
    for (int i = 0; i < 10; ++i) vof.step();
    vec("vof_alpha", vof.alpha());
    vec("vof_volume", {vof.volume()});
    std::vector<double> alpha(mesh.cell_count());
    for (std::size_t c = 0; c < alpha.size(); ++c) alpha[c] = std::tanh(10.0 * (0.3 - mesh.cells()[c].center.x));
    auto csf = fvm::csf_surface_tension(mesh, alpha, 0.07);
    std::vector<double> flat;
    for (const auto& v : csf) { flat.push_back(v.x); flat.push_back(v.y); flat.push_back(v.z); }
    vec("csf", flat);
    auto n = fvm::contact_angle_normal({0.3, 1.0, 0.0}, {0.0, 1.0, 0.0}, 1.0);
    auto n2 = fvm::contact_angle_normal({0.0, 2.0, 0.0}, {0.0, 1.0, 0.0}, 0.5);
    vec("contact_angle", {n.x, n.y, n.z, n2.x, n2.y, n2.z});
    std::cout << "\n}\n";
}
