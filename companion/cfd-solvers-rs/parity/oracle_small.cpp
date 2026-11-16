// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream DG, kinetic, atomistics and acoustics post-processing units.
#include "cfd/acoustics/postprocess.hpp"
#include "cfd/atomistics/md1d.hpp"
#include "cfd/atomistics/neighbor_list.hpp"
#include "cfd/dg/advection1d.hpp"
#include "cfd/dg/high_order_fem.hpp"
#include "cfd/kinetic/vlasov1d1v.hpp"
#include "cfd/kinetic/vlasov_poisson1d1v.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) std::cout << (i++ ? "," : "") << std::setprecision(17) << static_cast<double>(x);
    std::cout << "]";
}

int main() {
    std::cout << "{\n";
    for (std::size_t p : {1, 2, 4, 6}) vec("lgl" + std::to_string(p), cfd::dg::detail::gauss_lobatto_nodes(p));
    cfd::dg::DgAdvection1D dg({6, 4, 2.0, -0.7});
    dg.initialize([](double x) { return std::sin(M_PI * x) + 0.3 * (x > 1.0); });
    double m0 = dg.mass();
    for (int i = 0; i < 25; ++i) dg.step(0.01);
    vec("dg", dg.values());
    vec("dg_mass", std::vector<double>{m0, dg.mass()});
    cfd::dg::HighOrderNodalBasis1D basis(5);
    vec("cheb_nodes", basis.nodes());
    vec("cheb_d", basis.derivative_matrix());
    vec("cheb_element", basis.element_nodes(0.5, 2.0));
    std::vector<cfd::dg::OneDimensionalInterval> mesh{{0, 1, 0}, {1, 3, 1}, {3, 4, 0}};
    std::vector<std::size_t> marked{2, 0};
    std::vector<double> amr;
    for (auto c : basis.refine_marked(mesh, marked)) { amr.push_back(c.left); amr.push_back(c.right); amr.push_back(c.level); }
    vec("amr", amr);
    cfd::dg::SumFactorizedDerivative1D sf(basis);
    std::vector<double> in(6), out(6);
    for (int i = 0; i < 6; ++i) in[i] = std::exp(basis.nodes()[i]);
    sf.apply(in, out, 2.5);
    vec("sumfact", out);

    cfd::kinetic::Vlasov1D1V vl({8, 10, 2.0, 3.0});
    vl.initialize([](double x, double v) { return std::exp(-v * v) * (1.0 + 0.2 * std::cos(M_PI * x)) - 0.05; });
    for (int i = 0; i < 6; ++i) vl.step(0.02, i % 2 ? -0.4 : 0.6);
    vec("vlasov", vl.distribution());
    vec("vlasov_mass", std::vector<double>{vl.mass()});
    cfd::kinetic::VlasovPoisson1D1V vp({8, 12, 4.0 * M_PI, 4.0, 0.0});
    vp.initialize([](double x, double v) { return std::exp(-0.5 * v * v) / std::sqrt(2 * M_PI) * (1.0 + 0.1 * std::cos(0.5 * x)); });
    for (int i = 0; i < 5; ++i) vp.step(0.05);
    vec("vp_f", vp.distribution());
    vec("vp_e", vp.electric_field());
    vec("vp_rho", vp.charge_density());
    vec("vp_scalars", std::vector<double>{vp.mass(), vp.field_energy()});

    cfd::atomistics::HarmonicChain1D chain({7, 3.5, 2.0, 40.0});
    chain.set_velocity(0, 0.3); chain.set_velocity(4, -0.2);
    for (int i = 0; i < 20; ++i) chain.step(0.01);
    auto en = chain.energy();
    vec("chain_x", chain.positions());
    vec("chain_v", chain.velocities());
    vec("chain_e", std::vector<double>{en.kinetic, en.potential, en.total()});
    cfd::atomistics::PeriodicNeighborList1D nl({10.0, 1.6, 0.3});
    std::vector<double> pos{0.1, 1.2, 2.0, 3.4, 9.5, 8.6, 5.0};
    nl.rebuild(pos);
    std::vector<double> pairs; for (auto [a, b] : nl.pairs()) { pairs.push_back(a); pairs.push_back(b); }
    vec("pairs", pairs);
    cfd::atomistics::LennardJones1DForceField lj({0.8, 1.0});
    std::vector<double> f(pos.size());
    auto diag = lj.compute(pos, nl, f);
    vec("lj_f", f);
    vec("lj_diag", std::vector<double>{diag.potential_energy, diag.virial, nl.minimum_image(7.3), nl.minimum_image(-12.0)});
    std::vector<double> vel{0.5, -0.2, 0.1, 0.7};
    cfd::atomistics::VelocityRescalingThermostat1D::apply(vel, 2.0, 0.3, 0.01, 0.1);
    double L = 10.0;
    std::vector<double> bpos{1.0, 2.0, 7.0};
    double factor = cfd::atomistics::IsotropicBarostat1D::scale(bpos, L, 2.0, 0.05, 0.5);
    vel.push_back(factor); vel.push_back(L); for (double p : bpos) vel.push_back(p);
    vec("thermo_baro", vel);

    std::vector<double> p1{0, 1, 3, 2, 0, -1, 0, 0}, p2{0, 0, 1, 3, 2, 0, -1, 0}, p3{1, 2};
    std::vector<cfd::acoustics::BeamformSample> sensors{{0.0, 0.0, p1}, {0.3, 0.1, p2}, {0.1, 0.4, p3}, {1.0, 1.0, {}}};
    auto bf = cfd::acoustics::delay_and_sum_beamform(sensors, 0.2, 0.2, 1e-4, 1500.0, 6);
    bf.push_back(cfd::acoustics::plane_wave_intensity_w_m2(1e3, 1000.0, 1500.0));
    bf.push_back(cfd::acoustics::radiation_pressure_pa(3e4, 1500.0, 0.4));
    bf.push_back(cfd::acoustics::acoustic_heating_w_m3(3e4, 2.5));
    vec("acoustics", bf);
    std::cout << "\n}\n";
}
