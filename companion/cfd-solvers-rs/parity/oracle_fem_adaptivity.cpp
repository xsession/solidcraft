// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream fem/adaptivity unit.
#include "cfd/fem/adaptivity.hpp"
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
    std::size_t i = 0;
    for (auto x : v) std::cout << (i++ ? "," : "") << std::setprecision(17) << static_cast<double>(x);
    std::cout << "]";
}
static void dump(const std::string& p, const Mesh2D& m) {
    std::vector<double> n, t, b;
    for (auto& x : m.nodes) { n.push_back(x.x); n.push_back(x.y); }
    for (auto& x : m.triangles) for (auto i : x.node) t.push_back(i);
    for (auto& x : m.boundary_edges) { b.push_back(x.node[0]); b.push_back(x.node[1]); b.push_back(x.patch); }
    vec(p + "_nodes", n); vec(p + "_tris", t); vec(p + "_bnd", b);
}

int main() {
    std::cout << "{\n";
    auto mesh = make_rectangle_tri_mesh(4, 3, 2.0, 1.5);
    auto u_of = [](Node2 p) { return p.x * p.x + 0.3 * p.x * p.y + 0.1 * p.y * p.y * p.y; };
    std::vector<double> u;
    for (auto& p : mesh.nodes) u.push_back(u_of(p));
    auto src = [](Node2 p) { return 1.0 + p.x + 0.2 * p.y * p.y; };
    auto est = estimate_poisson_error_tri3(mesh, u, src, 1.7);
    vec("eta", est.element_indicator);
    vec("eta_global", std::vector<double>{est.global_indicator});
    auto mark = mark_dorfler(est.element_indicator, 0.5);
    vec("mark", mark);
    vec("mark_custom", mark_dorfler(std::vector<double>{0.1, 3.0, 0.5, 2.0, 0.0}, 0.9));
    vec("mark_zero", mark_dorfler(std::vector<double>{0.0, 0.0, 0.0}, 0.5));
    auto r1 = refine_tri3_longest_edges(mesh, mark);
    dump("r1", r1);
    std::vector<unsigned char> all(r1.element_count(), 0);
    for (std::size_t i = 0; i < all.size(); i += 3) all[i] = 1;
    auto r2 = refine_tri3_longest_edges(r1, all);
    dump("r2", r2);
    std::vector<double> u2;
    for (auto& p : r2.nodes) u2.push_back(u_of(p));
    auto est2 = estimate_poisson_error_tri3(r2, u2, src, 1.0);
    vec("eta2_global", std::vector<double>{est2.global_indicator});
    std::cout << "\n}\n";
}
