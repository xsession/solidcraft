// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream reaction_path, reactor_network, speciation_block,
// advanced_thermo and database_io chemistry units.
#include "cfd/chemistry/advanced_thermo.hpp"
#include "cfd/chemistry/database_io.hpp"
#include "cfd/chemistry/reaction_path.hpp"
#include "cfd/chemistry/reactor_network.hpp"
#include "cfd/chemistry/speciation_block.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd::chemistry;
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) std::cout << (i++ ? "," : "") << std::setprecision(17) << static_cast<double>(x);
    std::cout << "]";
}

static ReactionNetwork make_net() {
    ReactionNetwork n({{"A", 0.03, 0}, {"B", 0.06, 0}, {"C", 0.06, 0}});
    ElementaryReaction r1;
    r1.reactants = {{0, 2.0}};
    r1.products = {{1, 1.0}};
    r1.forward = {3.0, 0.0, 0.0};
    n.add_reaction(r1);
    ElementaryReaction r2;
    r2.reactants = {{1, 1.0}};
    r2.products = {{2, 1.0}};
    r2.forward = {4.0, 0.5, 1000.0};
    n.add_reaction(r2);
    return n;
}

int main() {
    std::cout << "{\n";
    auto net = make_net();
    std::vector<double> c{1.0, 0.5, 0.2};
    std::vector<double> path;
    for (auto e : reaction_path(net, c, 300.0, 0.0)) { path.push_back(e.from); path.push_back(e.to); path.push_back(e.molar_flux); }
    vec("path", path);
    std::vector<double> path_cut;
    for (auto e : reaction_path(net, c, 300.0, 3.0)) { path_cut.push_back(e.from); path_cut.push_back(e.to); path_cut.push_back(e.molar_flux); }
    vec("path_cut", path_cut);

    WellStirredReactorNetwork wsr(make_net());
    wsr.add_reactor({2.0, 300.0, {1.0, 0.0, 0.0}});
    wsr.add_reactor({1.0, 350.0, {0.0, 0.4, 0.1}});
    wsr.connect({0, 1, 0.3});
    wsr.advance(0.5, 0.05);
    std::vector<double> wc;
    for (auto& r : wsr.reactors()) for (double x : r.concentrations) wc.push_back(x);
    vec("wsr_conc", wc);
    vec("wsr_total", wsr.total_species_amounts());

    MonolithicSpeciationBlock block;
    std::vector<double> logk{std::log(10.0), std::log(0.5), 2.0};
    auto sp = block.solve(0.1, logk);
    vec("spec_c", sp.concentrations);
    vec("spec_meta", std::vector<double>{double(sp.iterations), sp.converged ? 1.0 : 0.0});
    std::vector<double> init{0.01, 0.02, 0.03, 0.04};
    auto sp2 = block.solve(2.0, std::vector<double>{-1.0, 0.0, 1.0}, init);
    vec("spec2_c", sp2.concentrations);

    Nasa7Polynomial n2{{3.53100528, -1.23660988e-4, -5.02999433e-7, 2.43530612e-9, -1.40881235e-12, -1046.97628, 2.96747038}, 200.0, 1000.0};
    std::vector<double> nasa;
    for (double t : {250.0, 298.15, 600.0, 999.0}) { nasa.push_back(n2.cp_molar(t)); nasa.push_back(n2.h_molar(t)); nasa.push_back(n2.s_molar(t)); }
    vec("nasa", nasa);
    ShomatePolynomial h2o{30.092, 6.832514, 6.793435, -2.53448, 0.082139, -250.881, 223.3967, -241.8264, 500.0, 1700.0};
    std::vector<double> sho;
    for (double t : {500.0, 1000.0, 1700.0}) { sho.push_back(h2o.cp_molar(t)); sho.push_back(h2o.h_molar_kj(t)); sho.push_back(h2o.s_molar(t)); }
    vec("shomate", sho);

    BinaryPitzer11 nacl{0.0765, 0.2664, 0.00127};
    std::vector<double> pz;
    for (double m : {0.0, 1e-10, 0.1, 1.0, 3.0}) pz.push_back(nacl.mean_activity_coefficient(m, 298.15));
    vec("pitzer11", pz);

    std::vector<int> z{1, -1, 2, -2, 0};
    std::vector<double> m{0.5, 0.6, 0.1, 0.05, 0.2};
    MulticomponentPitzer mp({{0, 1, 0.0765, 0.2664, 0.00127}, {2, 3, 0.2, 2.9, 0.01}, {0, 3, 0.02, 1.1, 0.0}});
    vec("mpitzer", mp.activity_coefficients(z, m, 298.15));
    DaviesActivity dav;
    vec("davies", dav.activity_coefficients(z, m, 298.15));
    SitActivity sit({{0, 1, 0.03}, {2, 3, -0.1}, {4, 4, 0.05}});
    vec("sit", sit.activity_coefficients(z, m, 298.15));

    std::vector<double> gibbs;
    for (auto [a, b, t] : std::vector<std::array<double, 3>>{{-1000.0, 2000.0, 300.0}, {0.0, 2.0e5, 300.0}, {3.0e5, 0.0, 400.0}, {5.0, 5.0, 1000.0}}) {
        auto g = minimize_binary_ideal_gibbs(a, b, t);
        gibbs.push_back(g.mole_fraction_a); gibbs.push_back(g.mole_fraction_b); gibbs.push_back(g.gibbs_molar);
    }
    vec("gibbs", gibbs);
    vec("sorption", std::vector<double>{langmuir_surface_coverage(2.5, 0.3), langmuir_surface_coverage(0.0, 1.0),
                                        binary_ion_exchange_fraction(3.0, 0.2, 0.1), binary_ion_exchange_fraction(3.0, 0.0, 0.1),
                                        binary_ion_exchange_fraction(3.0, 0.2, 0.0)});

    const std::string phreeqc =
        "# test db\nSOLUTION_MASTER_SPECIES\nNa  Na+  22.99  1   # sodium\nCl Cl- 35.453 -1\n  Ca   Ca+2   0.04008  2\nbad row\nSOLUTION_SPECIES\nX X 12 0\n";
    auto sp_db = load_phreeqc_master_species(phreeqc);
    std::vector<double> pdb;
    for (auto& s : sp_db) { pdb.push_back(s.molar_mass); pdb.push_back(s.charge); }
    vec("phreeqc", pdb);

    const std::string yaml =
        "units: {length: m}\nspecies:\n- name: H2\n- name: O2\n- name: H2O\n- name: OH\nreactions:\n"
        "- equation: 2 H2 + O2 => 2H2O\n  rate-constant: {A: 1.5e3, b: 0.5, Ea: 2000.0}\n"
        "- equation: H2O <=> OH + 0.5 H2   # comment\n  rate-constant:\n    A: 7.0\n    b: -0.25\n    Ea: 500\n";
    auto mech = load_cantera_yaml(yaml);
    vec("cantera_rates", mech.reaction_rates(std::vector<double>{0.4, 0.3, 0.2, 0.1}, 800.0));
    std::vector<double> cst;
    for (auto& r : mech.reactions()) {
        for (auto& t : r.reactants) { cst.push_back(t.species); cst.push_back(t.coefficient); }
        cst.push_back(-1.0);
        for (auto& t : r.products) { cst.push_back(t.species); cst.push_back(t.coefficient); }
        cst.push_back(r.forward.pre_exponential); cst.push_back(r.forward.temperature_exponent); cst.push_back(r.forward.activation_energy);
    }
    vec("cantera_terms", cst);
    std::cout << "\n}\n";
}
