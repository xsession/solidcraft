// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream io/{checkpoint, openfoam_io, xdmf_contract}.
// Run with the files of parity/openfoam-case copied to /tmp/ofcase.
#include "cfd/fem/mesh2d.hpp"
#include "cfd/io/checkpoint.hpp"
#include "cfd/io/openfoam_io.hpp"
#include "cfd/io/xdmf_contract.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace cfd::io;
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
    for (char c : s) { if (c == '\n') std::cout << "\\n"; else std::cout << c; }
    std::cout << "\"";
}
static std::string slurp(const std::string& p) { std::ifstream in(p, std::ios::binary); std::stringstream s; s << in.rdbuf(); return s.str(); }

int main() {
    std::cout << "{\n";
    std::vector<double> a{1.5, -2.25, 1e-300}, b{3.0};
    std::vector<CheckpointFieldView> views{{"U", a}, {"p", b}};
    write_checkpoint("/tmp/oracle.chk", "cavity", 0.125, 42, views);
    auto bytes = slurp("/tmp/oracle.chk");
    std::vector<double> bv; for (unsigned char c : bytes) bv.push_back(c);
    vec("checkpoint_bytes", bv);
    auto mesh = read_openfoam_polymesh("/tmp/ofcase");
    std::vector<double> mv;
    for (auto& c : mesh.cells()) { mv.push_back(c.center.x); mv.push_back(c.center.y); mv.push_back(c.center.z); mv.push_back(c.volume); }
    for (auto& f : mesh.faces()) { mv.push_back(f.owner); mv.push_back(f.boundary() ? -1.0 : double(f.neighbour)); mv.push_back(f.boundary() ? double(f.patch) : -1.0); mv.push_back(f.center.x); mv.push_back(f.center.y); mv.push_back(f.center.z); mv.push_back(f.area.x); mv.push_back(f.area.y); mv.push_back(f.area.z); }
    vec("openfoam_mesh", mv);
    std::vector<double> sv = read_openfoam_vol_scalar_field("/tmp/ofcase/T", 2);
    auto pu = read_openfoam_vol_scalar_field("/tmp/ofcase/p", 3); sv.insert(sv.end(), pu.begin(), pu.end());
    vec("openfoam_scalar", sv);
    write_openfoam_vol_scalar_field("/tmp/oracle_s", "T", {300.123456789, -1e-7, 2.5e12});
    write_openfoam_vol_vector_field("/tmp/oracle_v", "U", {{1.0, 2.0 / 3.0, -0.5}, {1e20, 0.0, 3.14159265}});
    str("openfoam_scalar_text", slurp("/tmp/oracle_s"));
    str("openfoam_vector_text", slurp("/tmp/oracle_v"));
    auto m2 = cfd::fem::make_rectangle_tri_mesh(3, 2, 1.0, 1.0);
    auto c = make_xdmf_mesh_contract(m2, std::vector<double>(m2.node_count(), 1.0), "temp.K");
    vec("xdmf", std::vector<double>{double(c.node_count), double(c.triangle_count), double(c.field_count), double(valid_xdmf_identifier("bad name"))});
    std::cout << "\n}\n";
}
