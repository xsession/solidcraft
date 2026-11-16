// Test oracle only. Exercises upstream distributed/{cartesian, halo_grid, selective_halo, device_assignment}.
#include "cfd/distributed/cartesian.hpp"
#include "cfd/distributed/device_assignment.hpp"
#include "cfd/distributed/halo_grid.hpp"
#include "cfd/distributed/selective_halo.hpp"
#include "cfd/solvers/lbm/descriptors.hpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace cfd::distributed;
using D19 = cfd::lbm::D3Q19Descriptor;
using D27 = cfd::lbm::D3Q27Descriptor;
static bool first = true;
static void vec(const std::string& k, const std::vector<double>& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << std::setprecision(17) << v[i];
    std::cout << "]";
}
int main() {
    std::cout << "{\n";
    std::vector<double> g;
    for (auto [e, r, d] : std::vector<std::tuple<Extent3, std::size_t, int>>{{{64, 64, 64}, 8, 3}, {{120, 40, 10}, 12, 3}, {{100, 50, 1}, 6, 2}, {{7, 5, 3}, 30, 3}, {{256, 32, 32}, 16, 3}}) {
        auto p = choose_process_grid(e, r, d); g.insert(g.end(), {(double)p.x, (double)p.y, (double)p.z});
    }
    ProcessGrid pg{3, 2, 2};
    for (std::size_t r = 0; r < 12; ++r) { auto c = coords_from_rank(r, pg); g.insert(g.end(), {(double)c.x, (double)c.y, (double)c.z, (double)rank_from_coords(c, pg)}); }
    for (std::size_t c = 0; c < 4; ++c) { auto [b, n] = partition_axis(10, c, 4); g.push_back(b); g.push_back(n); }
    auto brick = make_brick({10, 7, 5}, pg, 7);
    g.insert(g.end(), {(double)brick.begin.x, (double)brick.begin.y, (double)brick.begin.z, (double)brick.extent.x, (double)brick.extent.y, (double)brick.extent.z});
    for (auto o : neighbor_offsets()) { g.push_back(neighbor_tag(o)); g.push_back(neighbor_rank(brick, o)); g.push_back(neighbor_rank(brick, o, {false, true, false})); }
    vec("cartesian", g);
    HaloSoA3D<double, 2> h({3, 2, 2});
    for (std::size_t c = 0; c < 2; ++c) for (std::size_t z = 0; z < 4; ++z) for (std::size_t y = 0; y < 4; ++y) for (std::size_t x = 0; x < 5; ++x) h(c, x, y, z) = 1000.0 * c + 100.0 * z + 10.0 * y + x;
    std::vector<double> hv;
    for (auto o : neighbor_offsets()) { auto p = h.pack(o); hv.push_back(p.size()); hv.insert(hv.end(), p.begin(), p.end()); }
    HaloSoA3D<double, 2> h2({3, 2, 2});
    for (auto o : neighbor_offsets()) { auto p = h.pack(o.opposite()); for (auto& x : p) x += 0.5; h2.unpack(o, p); }
    for (std::size_t c = 0; c < 2; ++c) for (std::size_t z = 0; z < 4; ++z) for (std::size_t y = 0; y < 4; ++y) for (std::size_t x = 0; x < 5; ++x) hv.push_back(h2(c, x, y, z));
    vec("halo", hv);
    HaloSoA3D<float, 19> s({3, 2, 2});
    for (std::size_t c = 0; c < 19; ++c) for (std::size_t z = 0; z < 4; ++z) for (std::size_t y = 0; y < 4; ++y) for (std::size_t x = 0; x < 5; ++x) s(c, x, y, z) = float(c * 100 + z * 20 + y * 5 + x);
    std::vector<double> sv;
    for (auto o : neighbor_offsets()) {
        auto cp = crossing_populations<D19>(o); sv.push_back(cp.size()); for (int d : cp) sv.push_back(d);
        auto p = pack_selective<D19>(s, o); for (float x : p) sv.push_back(x);
        auto e = selective_receive_entry<D19>(s.interior_extent(), o, 1); sv.insert(sv.end(), {(double)e.direction, (double)e.x, (double)e.y, (double)e.z});
    }
    HaloSoA3D<float, 19> s2({3, 2, 2});
    for (auto o : neighbor_offsets()) unpack_selective<D19>(s2, o, pack_selective<D19>(s, o.opposite()));
    for (std::size_t c = 0; c < 19; ++c) for (std::size_t z = 0; z < 4; ++z) for (std::size_t y = 0; y < 4; ++y) for (std::size_t x = 0; x < 5; ++x) sv.push_back(s2(c, x, y, z));
    for (auto st : {selective_halo_traffic<D19>({32, 16, 8}), full_halo_traffic<D19>({32, 16, 8}), selective_halo_traffic<D27, double>({10, 10, 10})}) sv.insert(sv.end(), {(double)st.messages, (double)st.values, (double)st.bytes});
    vec("selective", sv);
    std::vector<double> dv;
    for (auto [r, n, v] : std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>{{0, 4, 2}, {3, 4, 2}, {1, 2, 8}, {2, 3, 0}}) { auto a = assign_device(r, n, v); dv.insert(dv.end(), {(double)a.device_ordinal, (double)a.oversubscribed, (double)a.visible_devices}); }
    for (auto [r, n, v, m] : std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>{{0, 3, 8, 0}, {2, 3, 8, 0}, {1, 3, 8, 2}, {3, 4, 2, 0}}) {
        auto gr = assign_device_group(r, n, v, m); dv.push_back(gr.device_ordinals.size()); for (auto o : gr.device_ordinals) dv.push_back(o); dv.push_back(gr.oversubscribed);
        for (auto sl : split_rank_work(103, gr)) dv.insert(dv.end(), {(double)sl.device_ordinal, (double)sl.begin, (double)sl.end});
    }
    vec("devices", dv);
    std::cout << "\n}\n";
}
