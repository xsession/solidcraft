// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream FDTD units (3-D, cylindrical, monitors, ports,
// post-processing, multiconductor lines, EMC sources, VTK output, material bridge).
#include "cfd/physics/material_bridge.hpp"
#include "cfd/solvers/fdtd/cylindrical_multigrid.hpp"
#include "cfd/solvers/fdtd/cylindrical_tm.hpp"
#include "cfd/solvers/fdtd/emc_sources.hpp"
#include "cfd/solvers/fdtd/hdf5_xdmf.hpp"
#include "cfd/solvers/fdtd/maxwell3d.hpp"
#include "cfd/solvers/fdtd/monitor.hpp"
#include "cfd/solvers/fdtd/multiconductor_tline.hpp"
#include "cfd/solvers/fdtd/port.hpp"
#include "cfd/solvers/fdtd/postprocess.hpp"
#include "cfd/solvers/fdtd/vtk_io.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace cfd::fdtd;
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
template <class T> static std::vector<double> cat(std::vector<double> a, const T& b) { a.insert(a.end(), b.begin(), b.end()); return a; }
static std::vector<double> cx(std::complex<double> c) { return {c.real(), c.imag()}; }

int main() {
    std::cout << "{\n";
    Maxwell3DConfig c; c.nx = 8; c.ny = 7; c.nz = 6; c.dx = 1e-3; c.dy = 1.2e-3; c.dz = 0.9e-3; c.courant = 0.9; c.epsilon_r = 2.0; c.epsilon_ry = 3.0;
    Maxwell3D m(c);
    m.initialize_gaussian_ez(1.5, 0.3);
    m.add_soft_ez_source(3, 3, 2, 0.25);
    m.step(12);
    std::vector<double> f;
    for (auto* a : {&m.ex(), &m.ey(), &m.ez(), &m.hx(), &m.hy(), &m.hz()}) f.insert(f.end(), a->begin(), a->end());
    vec("m3_fields", f);
    vec("m3_scalars", std::vector<double>{m.dt(), m.energy(), m.max_field()});
    Maxwell3DConfig cp = c; cp.boundary = Boundary3D::pmc; cp.epsilon_ry = 0.0; cp.mu_r = 1.3;
    Maxwell3D mp(cp);
    mp.initialize_gaussian_ez(1.0, 0.2);
    mp.step(9);
    vec("m3_pmc", std::vector<double>{mp.dt(), mp.energy(), mp.max_field(), mp.ez()[mp.nx() * mp.ny() * 2 + mp.nx() * 3 + 4], mp.hx()[mp.nx() * mp.ny() * 3 + mp.nx() * 2 + 1]});
    Maxwell3DConfig cv; cv.nx = 4; cv.ny = 4; cv.nz = 5; cv.dx = 2e-3; cv.dy = 1e-3; cv.dz = 1.5e-3;
    Maxwell3D mv(cv);
    mv.initialize_gaussian_ez(2.0, 0.5);
    mv.step(3);
    write_maxwell3d_vtk_ascii(mv, "/tmp/oracle_m3.vtk");
    std::ifstream in("/tmp/oracle_m3.vtk"); std::stringstream ss; ss << in.rdbuf();
    str("vtk", ss.str());
    vec("hdf5", std::vector<double>{double(hdf5_xdmf_output_available())});

    CylindricalTMConfig cc; cc.nr = 10; cc.nz = 12; cc.dr = 1e-3; cc.dz = 1.5e-3; cc.courant = 0.4; cc.epsilon_r = 2.2;
    CylindricalTM cyl(cc);
    cyl.initialize_gaussian_ez(0.0, 8e-3, 2e-3, 1.2);
    cyl.add_soft_ez_source(2, 5, 0.1);
    cyl.step(25);
    vec("cyl_ez", cyl.ez());
    vec("cyl_scalars", std::vector<double>{cyl.dt(), cyl.energy(), cyl.max_field()});
    std::vector<double> lv;
    for (auto l : make_cylindrical_multigrid_levels(40, 23, 0.02, 0.05, 6)) { lv.push_back(l.radial_cells); lv.push_back(l.axial_cells); lv.push_back(l.radial_spacing); lv.push_back(l.axial_spacing); }
    vec("levels", lv);

    DftMonitor dm(2e9);
    TimeProbe tp;
    for (int i = 0; i <= 200; ++i) { double t = i * 2e-12; double v = 0.7 * std::cos(2 * M_PI * 2e9 * t + 0.3) + 0.1 * std::sin(2 * M_PI * 5e9 * t); dm.sample(t, v); tp.sample(t, v); }
    vec("dft", cat(cx(dm.integral()), std::vector<double>{dm.amplitude(400e-12), double(dm.samples()), double(tp.time().size())}));
    WavePort1D p1(1e9, 50.0, 1), p2(1e9, 50.0, 1), p3(1e9, 75.0, -1);
    for (int i = 0; i <= 300; ++i) {
        double t = i * 5e-12, w = 2 * M_PI * 1e9;
        double a = std::sin(w * t), b = 0.3 * std::sin(w * t + 1.0), tr = 0.8 * std::sin(w * t - 0.7);
        p1.sample(t, a + b, (-a + b) / 50.0);
        p2.sample(t, tr, -tr / 50.0);
        p3.sample(t, 0.4 * a, 0.2 * b);
    }
    auto sp = s_parameters(p1, p2);
    auto inst = p3.instantaneous(1.0, 0.01);
    vec("ports", cat(cat(cat(cx(sp.s11), cx(sp.s21)), cat(cx(p3.forward_spectrum()), cx(p3.backward_spectrum()))), std::vector<double>{inst.forward, inst.backward}));

    std::vector<NearFieldSurfaceSample> surf;
    for (int i = 0; i < 6; ++i) {
        NearFieldSurfaceSample q;
        q.position = {0.01 * std::cos(i), 0.01 * std::sin(i), 0.002 * i};
        q.normal = {std::cos(i), std::sin(i), 0.1};
        q.area = 1e-5 * (1 + i);
        q.electric = {{0.1 * i, 0.2}, {0.3, -0.1 * i}, {0.05, 0.0}};
        q.magnetic = {{0.001, -0.002 * i}, {0.0005 * i, 0.001}, {-0.001, 0.0003}};
        surf.push_back(q);
    }
    auto ff = nf2ff(surf, {0.3, 0.4, 0.8}, 3e9);
    vec("nf2ff", cat(cat(cat(cx(ff.electric.x), cx(ff.electric.y)), cx(ff.electric.z)), std::vector<double>{ff.radiation_intensity}));
    SarGrid3D g; g.nx = 5; g.ny = 4; g.nz = 3; g.dx = 5e-3; g.dy = 6e-3; g.dz = 4e-3;
    for (int n = 0; n < 60; ++n) { g.ex.push_back(10.0 * std::sin(n)); g.ey.push_back(5.0 * std::cos(0.7 * n)); g.ez.push_back(1.0 + 0.1 * n); g.conductivity.push_back(0.5 + 0.01 * (n % 7)); g.density.push_back(1000.0 + 20.0 * (n % 5)); }
    vec("sar_local", local_sar(g));
    vec("sar_1g", sar_1g(g));
    vec("sar_10g", sar_10g(g));

    MulticonductorTlineConfig mc{3, {1e-7, 2e-7, 1.5e-7}, {1e-10, 2e-10, 0.5e-10}, 1e-10};
    MulticonductorTlineState st{{1.0, 0.0, -0.5}, {0.0, 0.01, 0.0}};
    for (int i = 0; i < 7; ++i) advance_multiconductor_tline(mc, st, std::vector<double>{0.2, -0.1, 0.0});
    vec("mtl", cat(st.voltage, st.current));
    vec("graded", cat(make_graded_axis(1.0, 6, 0.05, 1.3), make_graded_axis(2.0, 4, 0.1, 1.0)));
    std::vector<double> dw;
    for (auto w : deposit_thin_wire_subcell(2.3e-3, 0.7e-3, 9.9e-3, 4, 3, 5, 1e-3, 1e-3, 2e-3)) { dw.push_back(w.i); dw.push_back(w.j); dw.push_back(w.k); dw.push_back(w.weight); }
    for (auto w : deposit_thin_wire_subcell(-1.0, 1e-3, 0.0, 4, 3, 5, 1e-3, 1e-3, 2e-3)) { dw.push_back(w.i); dw.push_back(w.j); dw.push_back(w.k); dw.push_back(w.weight); }
    vec("deposit", dw);
    auto hs = huygens_equivalent_source({0.0, 2.0, 1.0}, {1.0, -2.0, 0.5}, {0.01, 0.02, -0.03});
    vec("huygens", cat(std::vector<double>(hs.electric_surface_current.begin(), hs.electric_surface_current.end()), hs.magnetic_surface_current));
    auto it = impedance_sheet_transmission({10.0, -40.0}, {377.0, 0.0}, {250.0, 5.0});
    std::vector<EmcTransferSample> ts{{1e6, {1.0, 0.0}}, {1e7, {0.5, -0.2}}, {1e8, {0.1, 0.3}}};
    vec("emc", cat(cat(cat(cx(it), cx(interpolate_emc_transfer(ts, 5e7))), cx(interpolate_emc_transfer(ts, 1e5))), cx(interpolate_emc_transfer(ts, 1e9))));
    cfd::physics::Material mat; mat.name = "x"; mat.density = 1000; mat.heat_capacity = 1000; mat.thermal_conductivity = 1; mat.relative_permittivity = {2.0, 3.0, 4.0};
    auto bc = cfd::physics::apply_to_fdtd(Maxwell3DConfig{}, mat);
    vec("bridge", std::vector<double>{bc.epsilon_rx, bc.epsilon_ry, bc.epsilon_rz, bc.epsilon_r});
    std::cout << "\n}\n";
}
