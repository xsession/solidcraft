// Test oracle only. Exercises upstream python/c_api (version string and inline XDMF writer).
#include "cfd/python/c_api.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
static std::string esc(const std::string& s) { std::string o; for (char c : s) { if (c == '"' || c == '\\') o += '\\'; if (c == '\n') { o += "\\n"; continue; } o += c; } return o; }
static std::string slurp(const char* p) { std::ifstream f(p); std::stringstream s; s << f.rdbuf(); return s.str(); }
int main() {
    const double xy[] = {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.1, 0.7, 1.0 / 3.0, -2.5e-7};
    const size_t tri[] = {0, 1, 2, 0, 2, 3, 1, 4, 0};
    const double s[] = {1.0, -0.5, 2.0 / 3.0, 1e20, 3.14159};
    char err[64] = {0};
    int a = cfd_write_xdmf_tri3("/tmp/cfd_capi_a.xmf", xy, 5, tri, 3, s, "temperature", err, sizeof err);
    int b = cfd_write_xdmf_tri3("/tmp/cfd_capi_b.xmf", xy, 5, tri, 3, nullptr, nullptr, err, sizeof err);
    int c = cfd_write_xdmf_tri3("/tmp/cfd_capi_c.xmf", xy, 0, tri, 3, nullptr, nullptr, err, sizeof err);
    std::cout << "{\"version\":\"" << cfd_solvers_version() << "\",\"codes\":[" << a << "," << b << "," << c << "],\"error\":\"" << esc(err)
              << "\",\"with_scalar\":\"" << esc(slurp("/tmp/cfd_capi_a.xmf")) << "\",\"without_scalar\":\"" << esc(slurp("/tmp/cfd_capi_b.xmf")) << "\"}\n";
}
