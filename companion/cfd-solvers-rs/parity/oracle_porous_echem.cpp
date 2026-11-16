// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream porous, aqueous-equilibrium and electrochemistry units.
#include "cfd/chemistry/aqueous_equilibrium.hpp"
#include "cfd/electrochemistry/corrosion_models.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"
#include "cfd/electrochemistry/implicit_pnp.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/porous/cahn_hilliard_calphad.hpp"
#include "cfd/porous/fractured_darcy_thmc.hpp"
#include "cfd/porous/thmc1d.hpp"
#include "cfd/solvers/electrochemistry/corrosion1d.hpp"
#include "cfd/solvers/electrochemistry/mixed_potential.hpp"
#include "cfd/solvers/electrochemistry/nernst_planck1d.hpp"
#include "cfd/solvers/electrochemistry/nernst_planck_poly.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd;
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) std::cout << (i++ ? "," : "") << std::setprecision(17) << static_cast<double>(x);
    std::cout << "]";
}
template <class T> static std::vector<double> cat(std::vector<double> a, const T& b) { a.insert(a.end(), b.begin(), b.end()); return a; }

int main() {
    std::cout << "{\n";
    porous::Thmc1D th({12, 2.0, 0.02, 0.015, 0.01, 0.3, 0.05, 0.2});
    th.initialize(1.0, 300.0, 0.5);
    th.set_cell(3, 1.5, 320.0, 0.9); th.set_cell(8, 0.8, 290.0, 0.0);
    for (int i = 0; i < 30; ++i) th.step(0.05);
    auto td = th.diagnostics();
    vec("thmc", cat(cat(cat(cat(std::vector<double>{td.pressure_mean, td.temperature_mean, td.species_inventory, td.reaction_heat}, th.pressure()), th.temperature()), th.species()), th.displacement()));
    porous::AllenCahn1D ac({40, 1.0, 1.5, 0.05});
    ac.initialize_tanh(0.4, 0.2);
    double e0 = ac.free_energy();
    for (int i = 0; i < 40; ++i) ac.step(2e-4);
    vec("allen_cahn", cat(std::vector<double>{e0, ac.free_energy()}, ac.phase()));
    porous::FracturedDarcyThmcNetwork fr({5});
    fr.add_edge({0, 1, 1.0, 1.0}); fr.add_edge({1, 2, 0.5, 4.0}); fr.add_edge({2, 3, 2.0, 0.5}); fr.add_edge({3, 4, 1.0, 1.0}); fr.add_edge({4, 0, 0.3, 2.0});
    fr.initialize(1.0, 0.3, 300.0, 0.2);
    fr.set_cell(0, 3.0, 0.9, 350.0, 1.0); fr.set_cell(3, 0.2, 0.0, 280.0, 0.0);
    for (int i = 0; i < 20; ++i) fr.step(0.02);
    auto fd = fr.diagnostics();
    vec("fractured", cat(cat(cat(cat(cat(std::vector<double>{fd.pressure_mass, fd.saturation_mass, fd.temperature_mass, fd.species_mass, fd.absolute_flux, fd.reaction_consumption}, fr.pressure()), fr.saturation()), fr.temperature()), fr.species()), fr.displacement()));
    porous::BinaryCalphadTable table({{0.0, 0.0}, {0.2, -0.05}, {0.5, -0.02}, {0.8, -0.06}, {1.0, 0.0}});
    porous::CahnHilliard1D ch({16, 1.0, 0.5, 0.002, 1.2, 0.0});
    ch.initialize(0.45, 0.05); ch.seed_nucleus(4, 1.5, 0.1);
    for (int i = 0; i < 25; ++i) ch.step(1e-4);
    porous::CahnHilliard1D ct({16, 1.0, 0.5, 0.002, 1.0, 0.0}, &table);
    ct.initialize(0.5, 0.1);
    for (int i = 0; i < 25; ++i) ct.step(1e-4);
    vec("cahn_hilliard", cat(cat(cat(cat(std::vector<double>{ch.mean(), ch.free_energy(), ct.mean(), ct.free_energy(), table.free_energy(0.35), table.chemical_potential(0.9)}, ch.composition()), ch.chemical_potential()), ct.composition()), ct.chemical_potential()));

    std::vector<chemistry::AcidFamily> fams{{0.01, 0, {1.8e-5}}, {0.002, 0, {7.5e-3, 6.2e-8, 4.8e-13}}};
    auto aq = chemistry::equilibrate_acids(fams, 0.003);
    std::vector<double> aqv{aq.ph, aq.hydrogen_mol_per_litre, aq.hydroxide_mol_per_litre, aq.charge_residual_mol_per_litre};
    for (auto& f : aq.family_concentrations) aqv = cat(aqv, f);
    vec("aqueous", aqv);

    using namespace cfd::electrochemistry;
    std::vector<double> b{nernst_potential(0.77, 310.0, 1.0, 3.5), butler_volmer_current_density(2.0, 0.05, 300.0, 2.0, 0.4, 0.6), faradaic_molar_flux(10.0, 2.0), corrosion_penetration_rate(5.0, 0.0558, 7874.0, 2.0)};
    auto ip = solve_implicit_pnp_electrode_cell(0.8, 0.02, {1e-3, 1e-9, 1e-9, 1.0, 1e-6, 0.0, 1e-6, 298.15, 1.0});
    b = cat(b, std::vector<double>{ip.concentration_mol_m3, ip.potential_v, double(ip.iterations), double(ip.converged)});
    CorrosionCell1DConfig cc; cc.metal_potential = 0.12; cc.bulk_electrolyte_potential = -0.01; cc.equilibrium_potential = 0.02;
    auto cr = solve_corrosion_cell_1d(cc);
    b = cat(b, std::vector<double>{cr.current_density, cr.surface_electrolyte_potential, cr.overpotential, cr.molar_dissolution_flux, cr.penetration_rate, double(cr.iterations)});
    std::vector<ElectrodeReaction> rx{{"Fe", -0.44, 1e-3, 1.0, 2.0, 0.5, 0.5}, {"O2", 0.4, 1e-6, 2.0, 4.0, 0.5, 0.5}, {"H2", 0.0, 1e-4, 0.5, 2.0, 0.5, 0.5}};
    auto mp = solve_mixed_potential(rx);
    b = cat(cat(cat(b, std::vector<double>{mp.potential, mp.total_current, double(mp.iterations)}), mp.reaction_current_density), mp.reaction_current);
    vec("echem_basic", b);

    std::vector<double> m;
    RedoxActivityTerm prod[] = {{0.5, 2.0}}, reac[] = {{0.1, 1.0}, {2.0, 3.0}};
    m.push_back(reaction_quotient(prod, reac));
    m.push_back(concentration_dependent_nernst(0.3, 300.0, 2.0, prod, reac));
    auto orr = oxygen_reduction_cathodic(0.6, 1.229, 0.2, 1e-4, 1e-6, 5e-3);
    auto her = hydrogen_evolution_cathodic(-0.5, 1.0, 1e-3, 1e-3, 2.0);
    m = cat(m, std::vector<double>{orr.equilibrium_potential, orr.activation_current_density, orr.limited_current_density, her.equilibrium_potential, her.activation_current_density, her.limited_current_density});
    m.push_back(charge_transfer_resistance(2.0, 300.0, 2.0));
    auto z1 = randles_impedance(100.0, 1.0, 20.0, 1e-4, 5.0), z0 = randles_impedance(0.0, 1.0, 20.0, 1e-4, 5.0);
    m = cat(m, std::vector<double>{z1.real(), z1.imag(), z0.real(), z0.imag()});
    DoubleLayerState dl{0.2, 0.01}; dl.advance_current(3.0, 0.01); m.push_back(dl.potential); m.push_back(dl.charge_density());
    PassivationModel pm{-0.4, 1e-2, -0.2, 0.8, 1e-3, 0.05, 2.0};
    for (double e : {-0.45, -0.3, 0.2, 1.0}) m.push_back(passivation_current_density(e, pm));
    ProductLayerState pl{1e-6, 1e5, 0.6}; pl.advance_from_anodic_current(4.0, 100.0, 0.0558, 7874.0, 2.0); m.push_back(pl.thickness); m.push_back(pl.area_specific_resistance());
    m = cat(m, std::vector<double>{faradaic_recession_distance(-3.0, 10.0, 0.0558, 7874.0, 2.0), cathodic_protection_current_density(-0.85, -0.6, 50.0), bruggeman_effective_transport(2e-9, 0.3), solid_phase_ohmic_drop(10.0, 1e-3, 100.0),
        marcus_current_density(1.0, 0.1, 6e4), marcus_current_density(1.0, -0.2, 6e4, 320.0, 2.0), bikerman_activity_correction(1.0, 0.02, 10.0)});
    auto mf = modified_nernst_planck_flux(1.0, 1.5, 1e-9, 20.0, -1, 300.0, 1e-4, 0.01, 20.0, 1e-6);
    m.push_back(mf.molar_flux); m.push_back(mf.activity);
    auto sw = smooth_electrode_interface(0.3, 0.5); m.push_back(sw.electrolyte_fraction); m.push_back(sw.delta);
    for (auto [l, r] : std::vector<std::pair<double, double>>{{-1, -2}, {1, 2}, {-0.3, 0.7}, {0.2, -0.6}}) { auto cc2 = cut_cell_electrode_interface(l, r, 1.0); m = cat(m, std::vector<double>{cc2.electrolyte_fraction, cc2.electrode_fraction, cc2.interface_area_fraction}); }
    auto bvmt = solve_anodic_bv_mass_transfer(10.0, 1e-4, -0.44, -0.3, 1e-2);
    m = cat(m, std::vector<double>{bvmt.surface_concentration, bvmt.current_density, bvmt.molar_flux, double(bvmt.iterations)});
    m.push_back(acid_speciation_redox_potential(0.8, 298.15, 2.0, 2.0, fams, 0.003));
    m.push_back(pitting_initiates(0.5, 0.2, {0.3, 0.1})); m.push_back(pitting_initiates(0.5, 0.0, {0.3, 0.1}));
    m = cat(m, std::vector<double>{stress_assisted_exchange_current(1e-3, 2e8, 1e-6), phase_field_corrosion_step(0.3, 2.0, 0.5, 1.5, 0.01, 0.1), recessed_level_set(-0.2, 0.05),
        marcus_hush_chidsey_current_density(1.0, 0.1, 6e4), marcus_hush_chidsey_current_density(1.0, -0.05, 4e4, 310.0, 1.0, 400)});
    m = cat(m, recession_remesh_1d(std::vector<double>{0.0, 0.1, 0.3, 0.7, 1.0}, 0.2));
    auto sc = coupled_stress_corrosion_step(1e-3, 1e8, 2e-6, 300.0, 50.0, 0.0558, 7874.0, 2.0); m.push_back(sc.stress_amplified_current_density); m.push_back(sc.recession_distance);
    PhaseFieldCorrosion2D p2(6, 5, 0.1, 0.12, 1.0, 0.002);
    p2.initialize_interface(0.2, 0.25, 0.15);
    std::vector<double> drv(30); for (int i = 0; i < 30; ++i) drv[i] = 0.1 * std::sin(double(i));
    for (int i = 0; i < 4; ++i) p2.step(0.05, drv);
    PhaseFieldCorrosion1D p1(8, 0.1, 1.0, 0.003);
    p1.initialize_interface(0.35, 0.1);
    for (int i = 0; i < 4; ++i) p1.step(0.05, std::vector<double>(drv.begin(), drv.begin() + 8));
    vec("corrosion_models", cat(cat(m, p2.order_parameter()), p1.order_parameter()));

    NernstPlanck1D np({20, 1e-3, 298.15, 2e-4, 1e-4, TransportBoundary1D::no_flux});
    np.set_potential([](double x) { return 0.02 * std::sin(3000.0 * x); });
    np.add_species("Na", 1, 1.3e-9, [](double x) { return 1.0 + 0.5 * std::exp(-1e7 * (x - 3e-4) * (x - 3e-4)); });
    np.add_species("Cl", -1, 2.0e-9, 1.2);
    np.step(40);
    NernstPlanck1D npp({16, 1.0, 300.0, 1e-3, -0.05, TransportBoundary1D::periodic});
    npp.add_species("A", 2, 0.01, [](double x) { return 1.0 + 0.3 * std::cos(2 * M_PI * x); });
    npp.set_potential([](double x) { return 0.001 * x; });
    npp.step(20);
    vec("np1d", cat(cat(cat(std::vector<double>{np.total_amount(0), np.total_amount(1), np.minimum_concentration(1), np.time()}, np.species()[0].concentration), np.species()[1].concentration), npp.species()[0].concentration));

    auto mesh = fvm::make_sheared_cartesian_hexa_mesh(5, 3, 1, 1e-3, 6e-4, 1e-4, 0.1);
    std::vector<double> poly;
    for (auto scheme : {ElectromigrationFluxScheme::centered, ElectromigrationFluxScheme::scharfetterGummel}) {
        NernstPlanckPolyMesh pm2(mesh, {298.15, 2e-3, 78.5, scheme});
        auto a = pm2.add_species("Cu", 2, 7e-10, [](fvm::Vec3 p) { return 10.0 + 4000.0 * p.x; });
        auto s2 = pm2.add_species("SO4", -2, 1.0e-9, 10.5);
        pm2.set_species_butler_volmer_boundary(a, "left", 0.05, 0.0, 2.0, 2);
        pm2.set_species_boundary(a, "right", SpeciesBoundaryType::fixedConcentration, 12.0);
        pm2.set_species_boundary(s2, "top", SpeciesBoundaryType::fixedMolarFlux, -1e-6);
        pm2.set_potential_boundary("right", PotentialBoundaryType::fixedPotential, 0.0);
        pm2.set_potential([](fvm::Vec3 p) { return -5.0 * p.x; });
        std::vector<double> flux(mesh.face_count());
        for (std::size_t f = 0; f < flux.size(); ++f) flux[f] = mesh.faces()[f].boundary() ? 0.0 : 1e-12 * mesh.faces()[f].area.x;
        pm2.set_face_flux(flux);
        pm2.step(5);
        auto en = pm2.solve_electroneutral_potential();
        poly = cat(poly, std::vector<double>{double(en.iterations), pm2.total_amount(a), pm2.total_amount(s2), pm2.minimum_concentration(a), pm2.total_charge()});
        poly = cat(cat(cat(poly, pm2.species()[0].concentration), pm2.species()[1].concentration), pm2.potential());
        poly = cat(poly, pm2.current_flux_faces());
    }
    NernstPlanckPolyMesh ps(mesh, {298.15, 1e-3, 78.5, ElectromigrationFluxScheme::scharfetterGummel});
    ps.add_species("K", 1, 2e-9, [](fvm::Vec3 p) { return 1.0 + 1e-11 * std::sin(5000.0 * p.x); });
    ps.add_species("Cl", -1, 2e-9, 1.0);
    ps.set_potential_boundary("left", PotentialBoundaryType::fixedPotential, 1e-3);
    ps.set_potential_boundary("right", PotentialBoundaryType::fixedPotential, 0.0);
    auto pr = ps.solve_poisson_potential();
    ps.step_poisson_nernst_planck(2);
    poly = cat(cat(cat(poly, std::vector<double>{double(pr.iterations), ps.time()}), ps.potential()), ps.species()[0].concentration);
    vec("np_poly", poly);
    std::cout << "\n}\n";
}
