// Test oracle only. Exercises upstream solvers/lbm/distributed_pull (virtual exchange).
#include "cfd/solvers/lbm/distributed_pull.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace cfd::lbm;
using namespace cfd::distributed;
static bool first = true;
static void vec(const std::string& k, const std::vector<double>& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << std::setprecision(17) << v[i];
    std::cout << "]";
}
template <class D> static std::vector<double> run(Extent3 global, std::size_t ranks, bool selective, DistributedLbmConfig cfg) {
    cfg.global = global;
    auto grid = choose_process_grid(global, ranks, D::dimensions);
    std::vector<DistributedPullBlock<D>> blocks;
    for (std::size_t r = 0; r < grid.size(); ++r) { blocks.emplace_back(cfg, make_brick(global, grid, r)); blocks.back().initialize_taylor_green(0.04F); }
    std::vector<double> out{(double)grid.x, (double)grid.y, (double)grid.z};
    for (int s = 0; s < 6; ++s) { if (selective) virtual_selective_step(blocks); else virtual_step(blocks, {true, true, true}); }
    double mass = 0; for (auto& b : blocks) mass += b.local_mass(); out.push_back(mass);
    auto m = gather_virtual_macroscopic(blocks);
    for (std::size_t i = 0; i < m.rho.size(); ++i) { out.push_back(m.rho[i]); out.push_back(m.ux[i]); out.push_back(m.uy[i]); out.push_back(m.uz[i]); }
    return out;
}
int main() {
    std::cout << "{\n";
    DistributedLbmConfig c; c.tau = 0.7F; c.acceleration_x = 1e-5F; c.acceleration_y = -2e-6F;
    vec("d3q19_full", run<D3Q19Descriptor>({8, 6, 4}, 4, false, c));
    vec("d3q19_sel", run<D3Q19Descriptor>({8, 6, 4}, 4, true, c));
    vec("d3q27_sel", run<D3Q27Descriptor>({6, 5, 4}, 3, true, c));
    vec("d2q9_full", run<D2Q9InPlaceDescriptor>({12, 10, 1}, 4, false, c));
    DistributedPullBlock<D3Q19Descriptor> b(DistributedLbmConfig{{4, 3, 2}, 0.8F, 1e-4F, 0, 0}, make_brick({4, 3, 2}, {1, 1, 1}, 0));
    b.initialize_taylor_green(0.02F); b.compute_all(); b.finish_step();
    b.save_checkpoint("/tmp/oracle_lbm.chk");
    std::ifstream in("/tmp/oracle_lbm.chk", std::ios::binary); std::stringstream ss; ss << in.rdbuf(); auto bytes = ss.str();
    std::vector<double> bv; for (unsigned char ch : bytes) bv.push_back(ch);
    vec("checkpoint", bv);
    std::cout << "\n}\n";
}
