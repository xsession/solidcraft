// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream FEM solver/kinematics units and prints JSON goldens.
#include "cfd/fem/mesh2d.hpp"
#include "cfd/fem/nonlinear_geometry.hpp"
#include "cfd/fem/phase_change.hpp"
#include "cfd/physics/material_bridge.hpp"
#include "cfd/solvers/fem/axisymmetric_elasticity.hpp"
#include "cfd/solvers/fem/brinkman2d.hpp"
#include "cfd/solvers/fem/darcy2d.hpp"
#include "cfd/solvers/fem/dynamic_bar1d.hpp"
#include "cfd/solvers/fem/eddy_current2d.hpp"
#include "cfd/solvers/fem/magnetostatics2d.hpp"
#include "cfd/solvers/fem/mixed_stokes.hpp"
#include "cfd/solvers/fem/nonlinear_poisson2d.hpp"
#include "cfd/solvers/fem/stokes2d.hpp"
#include "cfd/solvers/fem/thermoelastic_block.hpp"
#include "cfd/solvers/fem/transport_stabilization.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd::fem;
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << std::setprecision(17) << static_cast<double>(v[i]);
    std::cout << "]";
}
static std::vector<double> d2(const std::vector<Displacement2>& v) { std::vector<double> o; for (auto x : v) { o.push_back(x.x); o.push_back(x.y); } return o; }

int main() {
    std::cout << "{\n";
    auto mesh = make_rectangle_tri_mesh(6, 5, 1.2, 1.0);
    Brinkman2D br(mesh, {1e-3, 2e-3, 1e-4, 3000, 1e-12});
    br.solve([](Node2 p) { return Displacement2{1.0 + p.y, -0.5 * p.x}; });
    vec("brinkman", d2(br.velocity()));
    PenaltyStokes2D st(mesh, {1.0, 500.0, 5000, 1e-12});
    st.set_dirichlet([](Node2 p) { return Displacement2{p.y >= 1.0 - 1e-12 ? 1.0 : 0.0, 0.0}; });
    st.solve([](Node2 p) { return Displacement2{0.0, -p.x}; });
    vec("stokes", d2(st.velocity()));
    vec("stokes_div", std::vector<double>{st.divergence_l2()});
    auto ms = assemble_mixed_stokes_tri3(mesh, mesh.triangles[7], 0.8, 1e-3);
    std::vector<double> m(ms.matrix.begin(), ms.matrix.end()); m.push_back(ms.area);
    vec("mixed_stokes", m);
    auto tb = assemble_thermoelastic_block(0.3, 2e5, 4.0, 1e6, 0.1, 1.2e-5);
    vec("thermoelastic", tb.matrix);
    auto s1 = supg_parameters(0.01, 2.0, 1e-4, 1e-3), s2 = supg_parameters(0.02, -0.5, 0.0, 1e-2);
    vec("supg", std::vector<double>{s1.tau_s, s1.element_peclet, s2.tau_s, s2.element_peclet, dg_upwind_flux(1.0, 3.0, 2.0, 0.5), dg_upwind_flux(1.0, 3.0, -2.0, 0.5)});
    DynamicBar1DConfig bc; bc.elements = 12; bc.dt = 2e-6;
    DynamicBar1D bar(bc);
    bar.initialize([](double x) { return 1e-5 * x; }, [](double x) { return 0.1 * x * x; });
    bar.run(15, 2e3);
    std::vector<double> bv = bar.displacement(); for (double v : bar.velocity()) bv.push_back(v);
    bv.push_back(bar.time()); bv.push_back(bar.kinetic_energy()); bv.push_back(bar.strain_energy());
    vec("bar", bv);
    EddyCurrent2D ec(mesh, {5.0, 4000, 50, 1e-11});
    ec.solve([](Node2 p) { return p.x > 0.6 ? 0.05 : 1.0; }, [](Node2 p) { return p.x > 0.6 ? 1.0 : 0.0; }, [](Node2 p) { return std::complex<double>(p.x < 0.4 ? 1.0 : 0.0, p.y < 0.3 ? 0.2 : 0.0); });
    std::vector<double> ev; for (auto a : ec.vector_potential()) { ev.push_back(a.real()); ev.push_back(a.imag()); }
    vec("eddy", ev);
    Darcy2DConfig dc; dc.viscosity = 2.0; dc.density = 0.1; dc.gravity = {0.0, -9.81}; dc.linear = {4000, 1e-12};
    Darcy2D da(mesh, dc);
    da.set_boundary(0, {ScalarBoundaryType::dirichlet, [](Node2 p) { return 2.0 + 0.5 * p.y; }, {}});
    da.set_boundary(1, {ScalarBoundaryType::dirichlet, [](Node2) { return 1.0; }, {}});
    da.set_boundary(2, {ScalarBoundaryType::neumann, [](Node2) { return 0.0; }, {}});
    auto perm = [](Node2 p) { return 1.0 + p.x; };
    da.solve(perm, [](Node2 p) { return 0.3 * p.y; });
    vec("darcy_p", da.pressure());
    std::vector<double> dv; for (auto v : da.element_velocity(perm)) { dv.push_back(v.x); dv.push_back(v.y); }
    vec("darcy_v", dv);
    auto ring = make_rectangle_tri_mesh(5, 4, 1.0, 2.0);
    for (auto& n : ring.nodes) n.x += 0.5;
    AxisymmetricElasticity ax(ring, {2e5, 0.28, 5000, 1e-12});
    ax.set_dirichlet([](Node2 p) { return p.y <= 1e-12; }, [](Node2) { return Displacement2{0.0, 0.0}; }, false, true);
    ax.set_dirichlet([](Node2 p) { return p.x <= 0.5 + 1e-12; }, [](Node2) { return Displacement2{1e-3, 0.0}; }, true, false);
    ax.solve([](Node2 p) { return Displacement2{10.0 * p.x, -5.0}; });
    vec("axisym", d2(ax.displacement()));
    Magnetostatics2D mg(mesh, {4000, 1e-12});
    mg.set_boundary(0, {ScalarBoundaryType::dirichlet, [](Node2) { return 0.0; }, {}});
    mg.set_boundary(1, {ScalarBoundaryType::dirichlet, [](Node2) { return 0.0; }, {}});
    mg.set_boundary(3, {ScalarBoundaryType::dirichlet, [](Node2) { return 0.0; }, {}});
    RectangularStrandedCoil2D coil{0.3, 0.7, 0.2, 0.6, 50.0, 2.0};
    PiecewiseLinearBHCurve bh({0.0, 0.5, 1.2, 1.6}, {0.0, 100.0, 600.0, 4000.0});
    const double mu0 = 4e-7 * M_PI;
    auto nl = mg.solve_nonlinear([&](double b) { return std::min(bh.reluctivity(b), 1.0 / mu0); }, [&](Node2 p) { return coil.current_density_a_per_m2(p) * 1e3; }, 30, 1e-8, 0.6);
    vec("mag_a", mg.vector_potential());
    std::vector<double> bf; for (auto b : mg.element_flux_density()) { bf.push_back(b.x); bf.push_back(b.y); }
    vec("mag_b", bf);
    auto ft = maxwell_stress_boundary_force_2d(mesh, mg.element_flux_density(), 2, mu0, {0.6, 0.5});
    vec("mag_misc", std::vector<double>{double(nl.iterations), nl.maximum_relative_reluctivity_change, nl.converged ? 1.0 : 0.0,
        coil_flux_linkage_wb_turn(mesh, mg.vector_potential(), coil), coil_inductance_h(mesh, mg.vector_potential(), coil),
        ft.force_x_n_per_m, ft.force_y_n_per_m, ft.torque_z_n, bh.field_strength(0.8), bh.field_strength(2.0), bh.reluctivity(0.0), bh.reluctivity(1.4)});
    NonlinearPoisson2DConfig nc; nc.diffusivity = 0.7; nc.cubic_coefficient = 3.0;
    NonlinearPoisson2D np(mesh, nc);
    np.solve([](Node2 p) { return 10.0 + 5.0 * p.x * p.y; }, [](Node2 p) { return 0.2 * p.x - 0.1 * p.y; });
    std::vector<double> nps = np.solution(); nps.push_back(double(np.nonlinear_result().iterations)); nps.push_back(np.nonlinear_result().residual_rms);
    vec("nonlinear_poisson", nps);
    auto tr = green_lagrange_truss(2.0, 2.3, 1e-4, 2e11);
    auto tl = total_lagrangian_triangle(mesh, mesh.triangles[3], {Node2{0.2, 0.0}, Node2{0.45, 0.02}, Node2{0.21, 0.23}}, 1e7, 0.3);
    std::vector<double> kin{tr.stretch, tr.green_lagrange_strain, tr.second_piola_stress, tr.axial_force, tr.tangent, tl.reference_area, tl.determinant};
    for (double v : tl.deformation_gradient) kin.push_back(v);
    for (double v : tl.green_lagrange_strain) kin.push_back(v);
    for (double v : tl.internal_force) kin.push_back(v);
    vec("kinematics", kin);
    PhaseChangeMaterial pcm{273.0, 278.0, 2000.0, 4000.0, 3.3e5};
    vec("phase_change", std::vector<double>{pcm.liquid_fraction(275.5), pcm.specific_enthalpy(270.0), pcm.specific_enthalpy(276.0), pcm.specific_enthalpy(290.0),
        pcm.effective_heat_capacity(276.0), pcm.temperature_from_enthalpy(pcm.specific_enthalpy(276.0)), pcm.temperature_from_enthalpy(-50.0), porous_effective_conductivity(0.3, 2.0, 0.6)});
    std::vector<double> mat;
    cfd::physics::Material good{"steel", 7800, 500, 40, 1e6, {2, 3, 4}, {1, 1, 1}, 1.5};
    mat.push_back(cfd::physics::fem_conductivity(good));
    mat.push_back(cfd::physics::optical_refractive_index(good));
    for (auto bad : {cfd::physics::Material{"", 1, 1, 1, 1, {1, 1, 1}, {1, 1, 1}, 1}, cfd::physics::Material{"x", 1, 1, 1, 1, {1, 0, 1}, {1, 1, 1}, 1}, cfd::physics::Material{"x", 1, 1, 1, 1, {1, 1, 1}, {1, 1, 1}, 0}}) {
        try { cfd::physics::validate_material(bad); mat.push_back(1); } catch (...) { mat.push_back(0); }
    }
    vec("material", mat);
    std::cout << "\n}\n";
}
