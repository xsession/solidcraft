// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream hybrid DAE graph, power-flow, workflow runner and
// acoustic coupling units and prints JSON goldens.
#include "cfd/acoustics/coupling.hpp"
#include "cfd/core/hybrid_dae_graph.hpp"
#include "cfd/power/ac_power_flow.hpp"
#include "cfd/power/multibus_power.hpp"
#include "cfd/workflow/runner.hpp"
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
static std::vector<double> flat(const cfd::core::DaeTrajectory& t) {
    std::vector<double> o{double(t.times.size()), double(t.events.size()), t.converged ? 1.0 : 0.0};
    for (std::size_t i = 0; i < t.times.size(); ++i) { o.push_back(t.times[i]); for (double v : t.states[i]) o.push_back(v); }
    return o;
}

int main() {
    std::cout << "{\n";
    // Hybrid DAE graph: damped oscillator x'' + 0.4 x' + 3 x = 0 with an algebraic mirror y = 2x + 1.
    cfd::core::HybridDaeGraph g;
    auto x = g.add_variable("x", 1.0), v = g.add_variable("v", 0.0), y = g.add_variable("y", 3.0, false);
    std::array<std::size_t, 2> rows{x, v};
    g.add_component("osc", rows, [=](double, std::span<const double> s, std::span<const double> d, std::span<double> o) { o[0] = d[x] - s[v]; o[1] = d[v] + 0.4 * s[v] + 3.0 * s[x]; });
    g.add_connection("mirror", y, y, x, 2.0, 1.0);
    std::vector<cfd::core::DaeEvent> ev{{"x_zero", [x](double, std::span<const double> s) { return s[x]; }}};
    auto traj = g.integrate(0.0, 0.05, 30, ev);
    auto ft = flat(traj);
    vec("dae_graph", ft);

    cfd::power::TwoBusPowerFlowConfig tb;
    auto r2 = cfd::power::solve_two_bus_power_flow(tb);
    tb.line_resistance = 0.05; tb.specified_active_power = -0.9;
    auto r3 = cfd::power::solve_two_bus_power_flow(tb, std::vector<double>{0.95, -0.1});
    vec("two_bus", std::vector<double>{r2.voltage, r2.angle, double(r2.solve.iterations), r2.solve.residual_rms, r3.voltage, r3.angle, double(r3.solve.iterations)});
    using cfd::power::PowerBus; using cfd::power::PowerBusType;
    cfd::power::MultiBusPowerFlow net({PowerBus{PowerBusType::slack, 1.02, 0.0, 0, 0}, PowerBus{PowerBusType::pq, 1.0, 0.0, -0.4, -0.1}, PowerBus{PowerBusType::pq, 1.0, 0.0, -0.3, -0.15}, PowerBus{PowerBusType::pq, 1.0, 0.0, 0.1, 0.02}},
        {{0, 1, 0.01, 0.08}, {1, 2, 0.02, 0.1}, {0, 2, 0.015, 0.09}, {2, 3, 0.01, 0.05}});
    auto mb = net.solve();
    std::vector<double> mbv = mb.state; mbv.push_back(double(mb.solve.iterations)); mbv.push_back(double(mb.solve.converged));
    vec("multibus", mbv);
    std::vector<double> path;
    for (auto& p : net.continuation(0.5, 2.0, 4)) { path.push_back(p.load_factor); for (double s : p.result.state) path.push_back(s); }
    vec("continuation", path);
    std::vector<cfd::power::SwingBusParameters> sb{{2.0, 0.3, 0.5, 0.0, 0.0}, {1.5, 0.2, -0.2, 0.1, 0.0}};
    std::vector<cfd::power::SwingLine> sl{{0, 1, 1.2}, {0, cfd::power::invalid_power_index(), 0.8}, {1, cfd::power::invalid_power_index(), 0.5}};
    vec("swing", flat(cfd::power::integrate_swing_network(sb, sl, 0.05, 20)));
    std::vector<cfd::power::QuadraticDispatchGenerator> gens{{0.1, 1.0, 2.0, 0.5}, {0.0, 0.8, 1.5, 1.0}, {0.2, 0.5, 3.0, 0.1}};
    auto d1 = cfd::power::quadratic_dispatch(gens, 1.4), d2 = cfd::power::quadratic_dispatch(gens, 5.0);
    std::vector<double> dv = d1.generation; dv.push_back(d1.served); dv.push_back(d1.objective); dv.push_back(d1.feasible); dv.push_back(d2.feasible);
    vec("dispatch", dv);

    auto sweep = cfd::workflow::parameter_sweep(std::vector<double>{0.0, 0.5, 2.0}, [](double p) { return p * p - 1.0; });
    std::vector<double> sw; for (auto s : sweep) { sw.push_back(s.parameter); sw.push_back(s.objective); }
    auto inv = cfd::workflow::inverse_coordinate_search([](std::span<const double> p) { return (p[0] - 1.3) * (p[0] - 1.3) + 2.0 * (p[1] + 0.7) * (p[1] + 0.7) + 0.1 * p[0] * p[1]; }, {0.0, 0.0}, {0.5, 0.5}, 200, 1e-6);
    sw.insert(sw.end(), inv.parameters.begin(), inv.parameters.end());
    sw.push_back(inv.objective); sw.push_back(double(inv.iterations)); sw.push_back(inv.converged);
    vec("workflow", sw);

    auto mesh = cfd::fem::make_rectangle_tri_mesh(3, 2, 1.5, 1.0);
    std::vector<double> ep(mesh.boundary_edges.size());
    for (std::size_t e = 0; e < ep.size(); ++e) ep[e] = 100.0 + 10.0 * double(e);
    std::vector<double> loads;
    for (int patch : {-1, 2}) for (auto f : cfd::acoustics::pressure_to_structural_boundary_load(mesh, ep, patch)) { loads.push_back(f.x); loads.push_back(f.y); }
    vec("struct_load", loads);
    cfd::fem::DynamicBar1DConfig bc; bc.elements = 6; bc.dt = 1e-6;
    cfd::fem::DynamicBar1D bar(bc);
    bar.initialize([](double) { return 0.0; });
    cfd::acoustics::drive_dynamic_bar_from_pressure_trace(bar, std::vector<double>{1e5, 2e5, -5e4, 0.0}, 1e-3);
    vec("bar", bar.displacement());
    cfd::multiphysics::PennesBioheat2DConfig pc; pc.nx = 4; pc.ny = 3; pc.dt_s = 0.5;
    cfd::multiphysics::PennesBioheat2D pen(pc);
    pen.initialize(310.0);
    std::vector<double> inten(12), alpha(12);
    for (int i = 0; i < 12; ++i) { inten[i] = 1e4 * (1 + i % 3); alpha[i] = 5.0 + i; }
    cfd::acoustics::apply_acoustic_heating_to_pennes(pen, inten, alpha);
    pen.step(3);
    vec("pennes", pen.temperature_k());
    std::vector<double> rho(12, 1000.0), c(12, 1540.0);
    std::vector<double> st;
    for (auto a : cfd::acoustics::acoustic_streaming_acceleration(inten, alpha, rho, c, {1.0, 2.0, 2.0})) { st.push_back(a.x); st.push_back(a.y); st.push_back(a.z); }
    vec("streaming", st);
    std::vector<std::vector<double>> measured{{0.0, 1.0, 0.5, -0.2}, {0.3, 0.2, 0.1, 0.0}};
    auto sim = [](std::span<const double> p) { return std::vector<std::vector<double>>{{0.0, p[0], 0.5 * p[0], -0.2 * p[1]}, {0.3 * p[1], 0.2, 0.1 * p[0] * p[1], 0.0}}; };
    std::vector<double> fit{cfd::acoustics::sensor_trace_l2_misfit(sim(std::vector<double>{0.8, 1.1}), measured, true), cfd::acoustics::sensor_trace_l2_misfit(sim(std::vector<double>{0.8, 1.1}), measured, false)};
    auto fr = cfd::acoustics::fit_acoustic_parameters(sim, measured, {0.2, 0.4}, {0.3, 0.3}, 100, 1e-7);
    fit.insert(fit.end(), fr.parameters.begin(), fr.parameters.end());
    fit.push_back(fr.misfit); fit.push_back(double(fr.iterations)); fit.push_back(fr.converged);
    vec("fit", fit);
    std::cout << "\n}\n";
}
