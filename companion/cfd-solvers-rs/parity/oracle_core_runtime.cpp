// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream core runtime helpers and distributed/repartition.
#include "cfd/core/autotune.hpp"
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/decomposition.hpp"
#include "cfd/core/distributed_csr.hpp"
#include "cfd/core/distributed_sparse_assembly.hpp"
#include "cfd/core/mesh_partition.hpp"
#include "cfd/core/mixed_precision.hpp"
#include "cfd/core/numa.hpp"
#include "cfd/core/reverse_ad_distributed.hpp"
#include "cfd/distributed/repartition.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using namespace cfd::core;
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) std::cout << (i++ ? "," : "") << std::setprecision(17) << static_cast<double>(x);
    std::cout << "]";
}
static void str(const std::string& k, const std::string& s) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":\"";
    first = false;
    for (char c : s) { if (c == '\n') std::cout << "\\n"; else if (c == '\t') std::cout << "\\t"; else std::cout << c; }
    std::cout << "\"";
}
static CsrMatrix lap(std::size_t n) {
    std::vector<std::size_t> off{0}, col; std::vector<double> val;
    for (std::size_t r = 0; r < n; ++r) {
        if (r > 0) { col.push_back(r - 1); val.push_back(-1.0 - 0.01 * r); }
        col.push_back(r); val.push_back(4.0 + 0.1 * r);
        if (r + 1 < n) { col.push_back(r + 1); val.push_back(-1.2); }
        if (r + 3 < n) { col.push_back(r + 3); val.push_back(0.3); }
        off.push_back(col.size());
    }
    return CsrMatrix(n, n, off, col, val);
}

int main() {
    std::cout << "{\n";
    std::vector<double> p;
    for (std::size_t r = 0; r < 4; ++r) { auto s = partition_1d(23, r, 4, 2); p.push_back(s.begin); p.push_back(s.end); p.push_back(s.halo); }
    vec("partition_1d", p);
    auto a = lap(10);
    DistributedCsrPartition part(a, 3, 7);
    std::vector<double> hc(part.halo_columns().begin(), part.halo_columns().end());
    std::vector<double> local{1, -2, 0.5, 3}, halo; for (auto c : part.halo_columns()) halo.push_back(0.1 * c + 1); std::vector<double> y(4);
    part.multiply(local, halo, y);
    vec("dcsr", [&] { auto o = hc; o.insert(o.end(), y.begin(), y.end()); return o; }());
    DistributedCsrAssembler as(6, 6, {0, 2, 4, 6});
    for (int k = 0; k < 20; ++k) as.add((k * 5) % 6, (k * 7) % 6, 0.1 * k - 0.7, k % 3);
    std::vector<double> asm_out;
    for (std::size_t r = 0; r < 3; ++r) { auto b = as.build_owned(r); asm_out.insert(asm_out.end(), b.matrix.values().begin(), b.matrix.values().end()); for (auto c : b.matrix.column_indices()) asm_out.push_back(c); }
    auto g = as.build_global_reference(); asm_out.insert(asm_out.end(), g.values().begin(), g.values().end());
    for (auto o : g.row_offsets()) asm_out.push_back(o);
    vec("assembly", asm_out);
    std::vector<double> bb(10), x(10, 0.0); for (int i = 0; i < 10; ++i) bb[i] = std::sin(i + 1.0);
    auto rr = mixed_precision_iterative_refinement(a, bb, x, 12, 30, 1e-12);
    vec("refine", [&] { auto o = x; o.push_back(rr.converged); o.push_back(rr.outer_iterations); o.push_back(rr.relative_residual); return o; }());
    std::vector<std::array<std::size_t, 3>> els{{0, 1, 4}, {1, 5, 4}, {1, 2, 5}, {2, 6, 5}, {2, 3, 6}, {3, 7, 6}, {4, 5, 8}, {5, 9, 8}};
    auto plan = partition_tri3_mesh(els, 10, 3, std::vector<double>{1, 2, 1, 1, 3, 1, 1, 2});
    std::vector<double> mp(plan.element_boundaries.begin(), plan.element_boundaries.end());
    for (auto& r : plan.partitions) { mp.push_back(-1); for (auto n : r.owned_nodes) mp.push_back(n); mp.push_back(-2); for (auto n : r.halo_nodes) mp.push_back(n); mp.push_back(-3); for (auto n : r.shared_owned_nodes) mp.push_back(n); }
    vec("mesh_partition", mp);
    auto rp = cfd::distributed::plan_contiguous_repartition(std::vector<double>{1, 5, 2, 2, 0.5, 7, 1, 1, 3, 2}, 4, 2);
    std::vector<double> rv(rp.boundaries.begin(), rp.boundaries.end()); rv.insert(rv.end(), rp.loads.begin(), rp.loads.end()); rv.push_back(rp.average_load); rv.push_back(rp.imbalance);
    auto z = cfd::distributed::plan_contiguous_repartition(std::vector<double>(7, 0.0), 3, 1);
    rv.insert(rv.end(), z.boundaries.begin(), z.boundaries.end()); rv.push_back(z.imbalance);
    for (auto s : cfd::distributed::build_migration_segments(std::vector<std::size_t>{0, 3, 7, 10}, std::vector<std::size_t>{0, 4, 5, 10})) { rv.push_back(s.source_part); rv.push_back(s.destination_part); rv.push_back(s.global_begin); rv.push_back(s.global_end); }
    rv.push_back(cfd::distributed::repartition_recommended(1.4, 1.1)); rv.push_back(cfd::distributed::repartition_recommended(1.1, 1.0)); rv.push_back(cfd::distributed::repartition_recommended(1.3, 1.28));
    vec("repartition", rv);
    DistributedReverseOperator op(5);
    op.add_partition({{0, 1, 2}, [](ReverseTape&, std::span<const ReverseScalar> v) { return v[0] * v[1] + v[2] * v[2] * 0.5; }});
    op.add_partition({{2, 3, 4, 0}, [](ReverseTape&, std::span<const ReverseScalar> v) { return (v[0] - v[1]) * (v[0] - v[1]) / (v[2] + 2.0) + v[3] * 3.0; }});
    std::vector<double> gs{0.3, -1.2, 2.0, 0.7, 1.5}, gg(5);
    auto ro = op.evaluate(gs, gg);
    gg.push_back(ro.objective); gg.push_back(ro.partitions);
    vec("reverse_ad", gg);
    vec("cpu_list", parse_cpu_list(" 0-3, 8 ,6,10-11,2"));
    KernelTuningDatabase db;
    DeviceTuningKey dev{"Acme\tGPU", "X1", "1.2"};
    std::vector<KernelTuningCandidate> cands{{64, 1, 1}, {128, 2, 1}, {256, 4, 2}};
    db.autotune(dev, "spmv", cands, [](const KernelTuningCandidate& c) { return 1e-3 / c.work_group_size + 1e-6 * c.fused_stages; });
    db.upsert({{"B", "Y", "0"}, "stencil", {32, 1, 3}, 0.00123456789});
    db.save("/tmp/oracle_tune.tsv");
    std::ifstream in("/tmp/oracle_tune.tsv"); std::stringstream ss; ss << in.rdbuf();
    str("tuning_tsv", ss.str());
    str("tuning_key", device_tuning_key_string(dev));
    std::cout << "\n}\n";
}
