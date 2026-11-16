// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream projection2d, incompressible2d, thin_film and
// compressible_vof units plus the matrix-free CG, and prints JSON goldens.
#include "cfd/core/conjugate_gradient.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/solvers/fvm/compressible_vof.hpp"
#include "cfd/solvers/fvm/incompressible2d.hpp"
#include "cfd/solvers/fvm/projection2d.hpp"
#include "cfd/solvers/fvm/thin_film.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd;
static bool first = true;
template <class V>
static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << std::setprecision(17) << v[i];
    std::cout << "]";
}

int main() {
    std::cout << "{\n";
    // Matrix-free CG on a small SPD tridiagonal operator.
    std::vector<double> rhs{1.0, 2.0, -1.0, 0.5, 3.0, -2.0}, x(6, 0.0);
    core::ConjugateGradientWorkspace ws;
    auto apply = [](std::span<const double> in, std::span<double> out) {
        for (std::size_t i = 0; i < in.size(); ++i) out[i] = 4.0 * in[i] - (i ? in[i - 1] : 0.0) - (i + 1 < in.size() ? in[i + 1] : 0.0);
    };
    auto cg = core::conjugate_gradient(rhs, x, apply, ws, 100, 1e-12);
    std::vector<double> cgv = x;
    cgv.push_back(double(cg.iterations)); cgv.push_back(cg.converged ? 1.0 : 0.0);
    vec("cg", cgv);

    // Projection of a divergent field and of Taylor-Green.
    fvm::Projection2DConfig pc; pc.nx = 12; pc.ny = 10; pc.lx = 1.0; pc.ly = 0.8; pc.dt = 0.01; pc.density = 2.0;
    fvm::Projection2D proj(pc);
    proj.initialize_divergent(0.05);
    double div0 = proj.divergence_l2();
    proj.project();
    vec("proj_u", proj.u_faces());
    vec("proj_v", proj.v_faces());
    vec("proj_p", proj.pressure());
    vec("proj_metrics", std::vector<double>{div0, proj.divergence_l2(), proj.kinetic_energy(), double(proj.pressure_result().iterations)});
    fvm::Projection2D tg(pc);
    tg.initialize_taylor_green(0.05);
    double ke0 = tg.kinetic_energy();
    tg.project();
    vec("proj_tg", std::vector<double>{ke0, tg.kinetic_energy(), tg.divergence_l2()});

    // Incompressible Taylor-Green.
    fvm::Incompressible2DConfig ic; ic.nx = 16; ic.ny = 16; ic.kinematic_viscosity = 0.02; ic.dt = 1e-3;
    fvm::Incompressible2D ns(ic);
    ns.initialize_taylor_green(0.1);
    ns.run(20);
    vec("ns_u", ns.u_faces());
    vec("ns_v", ns.v_faces());
    vec("ns_p", ns.pressure());
    vec("ns_metrics", std::vector<double>{ns.time(), ns.kinetic_energy(), ns.divergence_l2(), ns.taylor_green_velocity_error(0.1)});

    // Thin film: capillary levelling of a bump, closed and with boundary pressure.
    auto mesh = fvm::make_sheared_cartesian_hexa_mesh(8, 5, 1, 0.08, 0.05, 0.01, 0.0);
    fvm::ThinFilmConfig tc; tc.dt = 2e-3; tc.surface_tension = 0.07; tc.viscosity = 0.05; tc.body_acceleration = {0.0, -9.81, 0.0};
    fvm::ThinFilmTransport film(mesh, tc);
    film.initialize([](fvm::Vec3 p) { return 1e-3 + 4e-4 * std::exp(-2000.0 * ((p.x - 0.04) * (p.x - 0.04) + (p.y - 0.025) * (p.y - 0.025))); });
    for (int i = 0; i < 5; ++i) film.step();
    vec("film_h", film.thickness());
    vec("film_p", film.pressure());
    vec("film_flux", film.face_flux());
    vec("film_inventory", std::vector<double>{film.inventory()});
    std::vector<double> bp(mesh.face_count(), 0.0);
    for (std::size_t f = 0; f < bp.size(); ++f) bp[f] = 10.0 * mesh.faces()[f].center.x;
    film.set_boundary_pressure(bp);
    film.step();
    vec("film_open_h", film.thickness());

    // Compressible VOF.
    fvm::CompressibleVofConfig vc; vc.dt = 1e-3; vc.liquid_bulk_modulus = 2.0e7;
    fvm::CompressibleVofTransport cv(mesh, vc);
    cv.initialize([](fvm::Vec3 p) { return p.x < 0.03 ? 1.0 : 0.0; }, 1.2e5);
    std::vector<double> flux(mesh.face_count());
    for (std::size_t f = 0; f < flux.size(); ++f) {
        const auto& face = mesh.faces()[f];
        flux[f] = face.boundary() ? 0.0 : face.area.x * 0.4 * (1.0 + 5.0 * face.center.x) + face.area.y * 0.1;
    }
    cv.set_face_flux(flux);
    for (int i = 0; i < 6; ++i) cv.step();
    vec("cvof_alpha", cv.liquid_fraction());
    vec("cvof_density", cv.density());
    vec("cvof_pressure", cv.pressure());
    vec("cvof_phase_density", [&] { std::vector<double> o; for (std::size_t c = 0; c < mesh.cell_count(); ++c) { o.push_back(cv.liquid_density()[c]); o.push_back(cv.gas_density()[c]); } return o; }());
    vec("cvof_totals", std::vector<double>{cv.liquid_volume(), cv.total_mass(), cv.minimum_pressure(), cv.maximum_pressure(), cv.time()});
    std::cout << "\n}\n";
}
