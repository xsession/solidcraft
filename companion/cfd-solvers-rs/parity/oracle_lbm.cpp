// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream LBM CPU units and prints JSON goldens.
#include "cfd/io/mesh_io.hpp"
#include "cfd/solvers/lbm/advanced_d2q9.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"
#include "cfd/solvers/lbm/descriptors.hpp"
#include "cfd/solvers/lbm/diagnostics.hpp"
#include "cfd/solvers/lbm/free_surface.hpp"
#include "cfd/solvers/lbm/particles.hpp"
#include "cfd/solvers/lbm/passive_scalar.hpp"
#include "cfd/solvers/lbm/visualization.hpp"
#include "cfd/solvers/lbm/voxelize.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace cfd;
static bool first = true;
template <class V>
static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << std::setprecision(17) << static_cast<double>(v[i]);
    std::cout << "]";
}
static void text(const std::string& k, const std::string& path) {
    std::ifstream in(path); std::ostringstream os; os << in.rdbuf();
    std::string s = os.str(), e;
    for (unsigned char c : s) { if (c == '"' || c == '\\') { e += '\\'; e += char(c); } else if (c == '\n') e += "\\n"; else if (c < 0x20) { char b[8]; std::snprintf(b, 8, "\\u%04x", c); e += b; } else e += char(c); }
    std::cout << ",\n\"" << k << "\":\"" << e << "\"";
}
template <class D> static std::vector<double> desc() {
    std::vector<double> o{double(D::dimensions), double(D::q)};
    for (int i = 0; i < D::q; ++i) { o.push_back(D::cx(i)); o.push_back(D::cy(i)); o.push_back(D::cz(i)); o.push_back(D::weight(i)); o.push_back(D::opposite(i)); }
    return o;
}
template <class S> static std::vector<double> adv(S&& s) {
    s.initialize_taylor_green(0.04F);
    std::vector<double> o{s.mass(), s.kinetic_energy()};
    s.step(25);
    o.push_back(s.mass()); o.push_back(s.kinetic_energy());
    return o;
}

int main() {
    std::cout << "{\n";
    vec("d2q9_inplace", desc<lbm::D2Q9InPlaceDescriptor>());
    vec("d3q19", desc<lbm::D3Q19Descriptor>());
    vec("d3q27", desc<lbm::D3Q27Descriptor>());

    // Momentum exchange on a channel with a square obstacle.
    lbm::D2Q9Solver solver({24, 16, 0.8F, 1e-5F, 0.0F});
    for (std::size_t y = 6; y < 10; ++y) for (std::size_t x = 8; x < 12; ++x) solver.set_solid(x, y);
    solver.initialize_uniform(1.0F, 0.02F, 0.0F);
    solver.step(30);
    auto ft = lbm::momentum_exchange_force(solver, 10.0, 8.0);
    vec("momentum_exchange", std::vector<double>{ft.force_x, ft.force_y, ft.torque_z});
    auto dfa = lbm::darcy_forchheimer_acceleration(0.03, -0.01, 0.05, 2e-3, 0.4);
    vec("darcy", std::vector<double>{dfa[0], dfa[1]});

    // Advanced collisions.
    lbm::AdvancedD2Q9Config ac{16, 12, 0.65F};
    vec("regularized", adv(lbm::RegularizedD2Q9Solver(ac)));
    vec("cumulant", adv(lbm::CumulantD2Q9Solver(ac)));
    vec("mrt", adv(lbm::MrtD2Q9Solver(ac)));

    // Passive scalar.
    lbm::PassiveScalarD2Q5 ps({20, 10, 0.02, 0.05, -0.02});
    ps.initialize([](double x, double y) { return std::exp(-40.0 * ((x - 0.3) * (x - 0.3) + (y - 0.5) * (y - 0.5))); });
    double t0 = ps.total();
    ps.step(40);
    vec("scalar", ps.scalar());
    vec("scalar_total", std::vector<double>{t0, ps.total()});

    // Voxelize a closed tetrahedron after a rigid transform, and read it back as STL/OBJ.
    io::TriangleSurface tet;
    fem::Point3 a{0.1, 0.1, 0.1}, b{0.9, 0.15, 0.1}, c{0.2, 0.85, 0.15}, d{0.3, 0.3, 0.9};
    tet.triangles = {{a, c, b}, {a, b, d}, {b, c, d}, {c, a, d}};
    auto moved = lbm::transform_surface(tet, {0.05, -0.02, 0.0}, {0.2, 0.3, 1.0}, 0.3);
    std::vector<double> mv;
    for (auto& t : moved.triangles) for (auto p : {t.a, t.b, t.c}) { mv.push_back(p.x); mv.push_back(p.y); mv.push_back(p.z); }
    vec("transformed", mv);
    lbm::VoxelGrid3D vg{10, 9, 8, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    vec("voxels", lbm::voxelize_surface(moved, vg));
    {
        std::ofstream o("/tmp/oracle_tet.stl");
        o << "solid t\n";
        for (auto& t : tet.triangles) { o << " facet normal 0 0 0\n  outer loop\n"; for (auto p : {t.a, t.b, t.c}) o << "   vertex " << std::setprecision(17) << p.x << ' ' << p.y << ' ' << p.z << "\n"; o << "  endloop\n endfacet\n"; }
        o << "endsolid t\n";
        std::ofstream j("/tmp/oracle_tet.obj");
        j << "# tet\nv 0.1 0.1 0.1\nv 0.9 0.15 0.1\nv 0.2 0.85 0.15\nv 0.3 0.3 0.9\nf 1 3 2\nf 1/1 2/2 4/4\nf 2//2 3//3 4//4\nf 3 1 4\n";
    }
    auto stl = io::read_ascii_stl("/tmp/oracle_tet.stl");
    auto obj = io::read_obj_triangles("/tmp/oracle_tet.obj");
    vec("stl_obj_voxels", [&] { auto x = lbm::voxelize_surface(stl, vg); auto y = lbm::voxelize_surface(obj, vg); x.insert(x.end(), y.begin(), y.end()); return x; }());

    // Free surface.
    const std::size_t fx = 8, fy = 6, fz = 5;
    std::vector<double> fill(fx * fy * fz);
    std::vector<float> ux(fill.size()), uy(fill.size()), uz(fill.size());
    for (std::size_t z = 0; z < fz; ++z) for (std::size_t y = 0; y < fy; ++y) for (std::size_t x = 0; x < fx; ++x) {
        const auto n = (z * fy + y) * fx + x;
        const double r = std::sqrt((x - 3.5) * (x - 3.5) + (y - 2.5) * (y - 2.5) + (z - 2.0) * (z - 2.0));
        fill[n] = std::clamp(2.2 - r, 0.0, 1.0);
        ux[n] = 0.3F + 0.05F * float(y); uy[n] = -0.1F * std::sin(float(x)); uz[n] = 0.05F;
    }
    auto fs = lbm::make_free_surface_field(fx, fy, fz, fill, 1e-8);
    double v0 = lbm::free_surface_volume(fs);
    for (int i = 0; i < 3; ++i) lbm::advect_free_surface_vof(fs, ux, uy, uz, 0.8, 1e-8);
    vec("fs_fill", fs.fill);
    vec("fs_phase", [&] { std::vector<double> o; for (auto p : fs.phase) o.push_back(double(static_cast<int>(p))); return o; }());
    vec("fs_volume", std::vector<double>{v0, lbm::free_surface_volume(fs)});
    auto flat = [](const std::vector<lbm::ParticleVector3>& v) { std::vector<double> o; for (auto p : v) { o.push_back(p.x); o.push_back(p.y); o.push_back(p.z); } return o; };
    vec("fs_normals", flat(lbm::free_surface_normals(fs)));
    vec("fs_curvature", lbm::free_surface_curvature(fs));
    vec("fs_csf", flat(lbm::free_surface_surface_tension_acceleration(fs, 0.07, 1.3)));

    // Immersed-boundary particles.
    auto iv = lbm::interpolate_particle_velocity(fx, fy, fz, ux, uy, uz, {-0.3, 6.7, 2.25});
    std::vector<lbm::ImmersedBoundaryParticle> ib{{{1.2, 2.5, 0.4}, {0.0, 0.1, 0.0}, 2.0, 0.5}, {{7.9, 5.9, 4.9}, {0.2, 0.0, -0.1}, 1.0, 2.0}};
    lbm::ParticleCouplingResult res;
    for (int i = 0; i < 3; ++i) res = lbm::advance_immersed_boundary_particles(fx, fy, fz, ux, uy, uz, ib, 0.4, true, {0.0, 0.0, -0.01});
    std::vector<double> ibv{iv.x, iv.y, iv.z};
    for (auto& p : ib) for (double v : {p.position.x, p.position.y, p.position.z, p.velocity.x, p.velocity.y, p.velocity.z}) ibv.push_back(v);
    ibv.push_back(res.total_reaction_force.x); ibv.push_back(res.total_reaction_force.y); ibv.push_back(res.total_reaction_force.z);
    vec("ib", ibv);
    vec("ib_ax", res.acceleration_x);
    vec("ib_az", res.acceleration_z);

    // Visualization.
    const std::size_t vx = 7, vy = 5, vz = 4;
    std::vector<lbm::Vector2> u2(vx * vy);
    std::vector<double> s2(vx * vy);
    for (std::size_t y = 0; y < vy; ++y) for (std::size_t x = 0; x < vx; ++x) { u2[y * vx + x] = {std::sin(0.7 * x) * std::cos(0.4 * y) + 0.1, -0.3 * x + 0.2 * y * y}; s2[y * vx + x] = 1.0 + 0.1 * x - 1e-7 * y * y * y + (x == 3 ? 1e-30 : 0.0); }
    vec("vorticity", lbm::vorticity_2d(vx, vy, u2, 0.5, 0.25));
    std::vector<lbm::Vector3> u3(vx * vy * vz);
    for (std::size_t z = 0; z < vz; ++z) for (std::size_t y = 0; y < vy; ++y) for (std::size_t x = 0; x < vx; ++x) u3[(z * vy + y) * vx + x] = {std::sin(0.5 * y + 0.2 * z), std::cos(0.3 * x) * z, 0.1 * x * y};
    vec("q_criterion", lbm::q_criterion_3d(vx, vy, vz, u3, 1.0, 0.5, 2.0));
    vec("slice", lbm::horizontal_slice(vx, vy, s2, 3));
    std::vector<double> sl;
    for (auto p : lbm::streamline_2d(vx, vy, u2, {0.5, 0.5}, 0.4, 30)) { sl.push_back(p.x); sl.push_back(p.y); }
    vec("streamline", sl);
    std::vector<double> s3(vx * vy * vz);
    std::vector<std::uint8_t> solid(s3.size());
    for (std::size_t i = 0; i < s3.size(); ++i) { s3[i] = 0.001 * double(i) * double(i) - 3.25 + (i == 5 ? 123456789.0 : 0.0) + (i == 6 ? 1e-9 : 0.0); solid[i] = i % 3 == 0; }
    lbm::write_vtk_structured_2d("/tmp/oracle_lbm2d.vtk", vx, vy, s2, u2, "rho");
    lbm::write_vtk_structured_3d("/tmp/oracle_lbm3d.vtk", vx, vy, vz, s3, u3, solid, "field");
    lbm::write_json_state_3d("/tmp/oracle_lbm_state.json", vx, vy, vz, 42, s3, solid, "speed");
    lbm::write_pgm_scalar("/tmp/oracle_lbm.pgm", vx, vy, s2);
    lbm::write_html_scalar_viewer("/tmp/oracle_lbm.html", vx, vy, s2);
    text("vtk2d", "/tmp/oracle_lbm2d.vtk");
    text("vtk3d", "/tmp/oracle_lbm3d.vtk");
    text("json_state", "/tmp/oracle_lbm_state.json");
    text("pgm", "/tmp/oracle_lbm.pgm");
    text("html", "/tmp/oracle_lbm.html");
    std::cout << "\n}\n";
}
