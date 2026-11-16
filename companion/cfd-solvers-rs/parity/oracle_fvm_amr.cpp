// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream fvm/{adaptive_mesh, combustion, mesh_motion, reordering}.
#include "cfd/fvm/adaptive_mesh.hpp"
#include "cfd/fvm/combustion.hpp"
#include "cfd/fvm/mesh_motion.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/fvm/reordering.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd::fvm;
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) std::cout << (i++ ? "," : "") << std::setprecision(17) << static_cast<double>(x);
    std::cout << "]";
}
static std::vector<double> dump(const AdaptiveHexMesh& m) {
    std::vector<double> o;
    for (auto& c : m.cells()) { o.push_back(c.lower.x); o.push_back(c.lower.y); o.push_back(c.lower.z); o.push_back(c.upper.x); o.push_back(c.upper.y); o.push_back(c.upper.z); o.push_back(c.level); o.push_back(c.root); o.push_back(double(c.lineage)); }
    return o;
}

int main() {
    std::cout << "{\n";
    std::vector<IdealGasSpeciesThermo> sp{{0.032, 918.0, 0.0}, {0.028, 1040.0, 0.0}, {0.018, 1996.0, -1.3435e7}};
    std::vector<double> y{0.2, 0.7, 0.1};
    auto st = ideal_gas_mixture_state(1200.0, 2e5, y, sp);
    vec("combustion", std::vector<double>{st.density, st.mixture_molar_mass, st.gas_constant, st.cp, st.cv, st.gamma, st.specific_enthalpy, temperature_from_mixture_enthalpy(st.specific_enthalpy, y, sp)});
    auto mesh = make_cartesian_hexa_mesh(3, 2, 2, 1.5, 1.0, 0.8);
    auto moved = rigidly_moved_mesh(mesh, {0.1, -0.2, 0.3}, {0.3, 0.4, 1.0}, 0.35);
    std::vector<double> mv;
    for (auto& c : moved.cells()) { mv.push_back(c.center.x); mv.push_back(c.center.y); mv.push_back(c.center.z); }
    for (auto& f : moved.faces()) { mv.push_back(f.area.x); mv.push_back(f.area.y); mv.push_back(f.area.z); }
    auto phi = ale_mesh_flux(mesh, moved, 0.01);
    mv.insert(mv.end(), phi.begin(), phi.end());
    std::vector<double> field; for (std::size_t i = 0; i < mesh.cell_count(); ++i) field.push_back(std::sin(1.0 + i));
    auto fine = make_cartesian_hexa_mesh(4, 3, 2, 1.5, 1.0, 0.8);
    auto rm = conservative_nearest_remap(mesh, field, fine);
    mv.insert(mv.end(), rm.begin(), rm.end());
    vec("motion", mv);
    auto sh = make_sheared_cartesian_hexa_mesh(4, 3, 2, 1.0, 1.0, 1.0, 0.3);
    auto order = morton_cell_order(sh);
    auto re = reorder_cells(sh, order);
    std::vector<double> ro(order.begin(), order.end());
    for (auto& f : re.faces()) { ro.push_back(f.owner); ro.push_back(f.boundary() ? -1.0 : double(f.neighbour)); }
    vec("reorder", ro);

    auto amr = AdaptiveHexMesh::cartesian(3, 2, 2, 1.5, 1.0, 1.0);
    std::vector<unsigned char> mk(amr.cell_count(), 0); mk[0] = 1; mk[7] = 1;
    auto r1 = MarkedHexRefinement(mk, 3).apply(amr);
    std::vector<unsigned char> mk2(r1.mesh.cell_count(), 0); mk2[1] = 1;
    auto r2 = MarkedHexRefinement(mk2, 3).apply(r1.mesh);
    auto a = dump(r2.mesh);
    a.push_back(satisfies_two_to_one(r2.mesh));
    for (auto& ws : r2.mapping.new_from_old) for (auto w : ws) { a.push_back(w.source_cell); a.push_back(w.volume); }
    vec("amr_refine", a);
    std::vector<double> f2; for (auto& c : r2.mesh.cells()) f2.push_back(std::exp(-((c.lower.x - 0.4) * (c.lower.x - 0.4) + c.lower.y * c.lower.y)) + 0.1 * c.upper.z);
    auto ind = scalar_jump_error_indicator(r2.mesh, f2);
    auto dm = mark_dorfler(ind, 0.6);
    std::vector<unsigned char> cm(r2.mesh.cell_count(), 1);
    auto co = MarkedHexCoarsening(cm, 0).apply(r2.mesh);
    auto b = ind; b.insert(b.end(), dm.begin(), dm.end());
    auto cd = dump(co.mesh); b.insert(b.end(), cd.begin(), cd.end());
    auto rem = conservative_adaptive_remap(r2.mesh, f2, co.mesh); b.insert(b.end(), rem.begin(), rem.end());
    vec("amr_coarsen", b);
    AdaptiveAmrConfig cfg; cfg.maximum_level = 2; cfg.refine_theta = 0.4; cfg.coarsen_relative_threshold = 0.3;
    auto ad = adapt_scalar_field(r2.mesh, f2, cfg);
    auto c = dump(ad.mesh); c.insert(c.end(), ad.field.begin(), ad.field.end()); c.insert(c.end(), ad.indicator.begin(), ad.indicator.end()); c.push_back(ad.refined_cell_count); c.push_back(ad.coarsened_cell_count);
    auto pm = ad.mesh.poly_mesh();
    for (auto& f : pm.faces()) { c.push_back(f.owner); c.push_back(f.area.x + 2 * f.area.y + 3 * f.area.z); c.push_back(f.center.x + f.center.y + f.center.z); }
    vec("amr_adapt", c);
    std::cout << "\n}\n";
}
