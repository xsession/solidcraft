// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream optics analysis/design units.
#include "cfd/optics/analysis.hpp"
#include "cfd/optics/birefringence.hpp"
#include "cfd/optics/diffraction.hpp"
#include "cfd/optics/freeform.hpp"
#include "cfd/optics/interchange.hpp"
#include "cfd/optics/nonlinear_polarization.hpp"
#include "cfd/optics/nonsequential.hpp"
#include "cfd/optics/optimization.hpp"
#include "cfd/optics/serialization.hpp"
#include "cfd/optics/system_models.hpp"
#include "cfd/optics/transform.hpp"
#include "cfd/optics/wavefront.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd::optics;
using cfd::optics::Vec3;
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
    for (char c : s) { if (c == '\n') std::cout << "\\n"; else if (c == '"') std::cout << "\\\""; else std::cout << c; }
    std::cout << "\"";
}
static std::vector<double> cplx(std::complex<double> c) { return {c.real(), c.imag()}; }
template <class T> static std::vector<double> cat(std::vector<double> a, const T& b) { a.insert(a.end(), b.begin(), b.end()); return a; }

static SequentialOpticalSystem lens() {
    SequentialOpticalSystem s(1.0);
    SequentialSurface a; a.type = SurfaceType::sphere; a.vertex_z = 0.0; a.radius = 50.0; a.aperture_radius = 10.0; a.refractive_index_after = 1.5;
    SequentialSurface b; b.type = SurfaceType::even_asphere; b.vertex_z = 5.0; b.radius = -50.0; b.aperture_radius = 10.0; b.refractive_index_after = 1.0; b.conic_constant = -0.5; b.even_coefficients = {1e-5, -2e-8, 0.0, 3e-12};
    SequentialSurface c; c.type = SurfaceType::plane; c.vertex_z = 60.0; c.aperture_radius = 30.0;
    s.add_surface(a); s.add_surface(b); s.add_surface(c);
    return s;
}

int main() {
    std::cout << "{\n";
    auto sys = lens();
    std::vector<Ray> rays;
    for (int i = -4; i <= 4; ++i) { Ray r; r.origin = {0.3 * i, 1.0 * i, -10.0}; r.direction = {0.0, 0.01, 1.0}; rays.push_back(r); }
    Ray bad; bad.origin = {0.0, 20.0, -10.0}; rays.push_back(bad);
    auto traced = sys.trace_many(rays);
    auto spot = spot_diagram_at_plane(traced, 70.0);
    std::vector<double> sp;
    for (auto p : spot.points) { sp.push_back(p.x); sp.push_back(p.y); }
    sp.push_back(spot.centroid_x); sp.push_back(spot.centroid_y); sp.push_back(spot.rms_radius);
    vec("spot", sp);
    std::vector<double> fan;
    for (auto p : tangential_ray_fan(traced, 70.0, 0.1)) { fan.push_back(p.pupil_coordinate); fan.push_back(p.transverse_aberration); }
    vec("fan", fan);
    vec("distortion", std::vector<double>{distortion_percent(10.3, 10.0)});

    std::vector<std::complex<double>> pupil;
    for (int i = 0; i < 16; ++i) { double x = (i - 7.5) / 8.0; pupil.push_back(std::abs(x) < 0.6 ? std::polar(1.0, 0.4 * x * x) : 0.0); }
    std::vector<double> fr;
    for (auto s : fraunhofer_1d(pupil, 1e-4, 550e-9, 0.5)) { fr.push_back(s.coordinate); fr.push_back(s.intensity); }
    vec("fraunhofer", fr);
    auto psf = normalized_psf_1d(pupil);
    vec("psf", psf);
    vec("mtf", mtf_from_psf_1d(psf));

    std::vector<double> zz;
    for (auto [n, m] : std::vector<std::pair<int, int>>{{0, 0}, {1, 1}, {1, -1}, {2, 0}, {3, -1}, {4, 2}, {5, 3}, {6, 0}, {3, 0}})
        zz.push_back(zernike(n, m, 0.73, 1.1));
    vec("zernike", zz);
    std::vector<WavefrontSample> ws;
    for (int i = 0; i < 40; ++i) {
        double rho = std::sqrt((i + 0.5) / 40.0), th = 2.399963 * i;
        ws.push_back({rho, th, 0.3 + 0.2 * rho * std::cos(th) - 0.1 * (2 * rho * rho - 1) + 0.05 * rho * rho * std::sin(2 * th) + 0.01 * std::sin(5.0 * i)});
    }
    std::vector<double> fz;
    for (auto c : fit_zernike(ws, 3)) { fz.push_back(c.n); fz.push_back(c.m); fz.push_back(c.value); }
    vec("fit_zernike", fz);
    vec("wf_rms", std::vector<double>{wavefront_rms(ws, true), wavefront_rms(ws, false), longitudinal_chromatic_shift(std::vector<double>{50.1, 49.8, 50.4})});

    double neff = uniaxial_extraordinary_index(1.658, 1.486, 0.6);
    vec("biref", std::vector<double>{neff, birefringent_retardance(1.658, neff, 2e-5, 633e-9)});
    auto j = apply_linear_retarder({{0.6, 0.1}, {0.3, -0.7}}, 1.2, 0.35);
    vec("retarder", cat(cplx(j.x), cplx(j.y)));

    NonlinearOpticalMaterial mat{2.25, {2e-22, 1e-24}, 5e12};
    auto r1 = evaluate_nonlinear_polarization({3e7, 1e7}, mat);
    auto r2 = evaluate_nonlinear_polarization({0.0, 0.0}, mat, 250.0);
    auto r3 = evaluate_nonlinear_polarization({1e8, 0.0}, NonlinearOpticalMaterial{1.5, {1e-21, 0.0}, 0.0});
    std::vector<double> nl;
    for (auto& r : {r1, r2, r3}) {
        nl.push_back(r.intensity_w_per_m2);
        nl = cat(nl, cplx(r.linear_polarization_c_per_m2)); nl = cat(nl, cplx(r.nonlinear_polarization_c_per_m2));
        nl = cat(nl, cplx(r.total_polarization_c_per_m2)); nl = cat(nl, cplx(r.effective_relative_permittivity));
    }
    vec("nonlinear", nl);

    auto tf = RigidTransform::from_euler_xyz({1.0, -2.0, 0.5}, 0.3, -0.2, 0.9);
    Ray tr; tr.origin = {0.4, 0.5, -1.0}; tr.direction = {0.1, 0.2, 0.97};
    auto ta = tf.apply(tr); auto ti = tf.inverse(ta);
    vec("transform", std::vector<double>{ta.origin.x, ta.origin.y, ta.origin.z, ta.direction.x, ta.direction.y, ta.direction.z, ti.origin.x, ti.origin.y, ti.origin.z});

    Aperture ap; std::vector<double> apv;
    for (auto k : {ApertureKind::circular, ApertureKind::rectangular, ApertureKind::elliptical}) {
        ap.kind = k; ap.x_radius = 1.0; ap.y_radius = 0.5;
        apv.push_back(ap.contains(0.8, 0.45)); apv.push_back(ap.contains(0.6, 0.3)); apv.push_back(ap.contains(0.95, 0.1));
    }
    WavelengthModel wm; wm.add({486.1, 1.0, false}); wm.add({587.6, 2.0, true}); wm.add({656.3, 1.0, false});
    apv.push_back(wm.primary_nm());
    vec("models", apv);

    PolynomialFreeformSurface ff(3.0, {{2, 0, 0.01}, {0, 2, 0.015}, {1, 2, -0.002}, {3, 1, 1e-4}}, 5.0);
    Ray fr1; fr1.origin = {0.5, -0.8, -2.0}; fr1.direction = normalized({0.05, 0.1, 1.0});
    Vec3 pt{}, nn{};
    bool hit = ff.intersect(fr1, pt, nn);
    vec("freeform", std::vector<double>{ff.sag(1.2, -0.7), double(hit), pt.x, pt.y, pt.z, nn.x, nn.y, nn.z});

    for (auto [fmt, name] : std::vector<std::pair<LensInterchangeFormat, std::string>>{{LensInterchangeFormat::zemax_zmx, "zmx"}, {LensInterchangeFormat::codev_seq, "seq"}, {LensInterchangeFormat::oslo_len, "len"}})
        str("export_" + name, export_lens_system(sys, fmt));
    auto imp = import_lens_system("! hdr\nSURF 0\nTYPE STANDARD\nCURV 0.02 # c\nDISZ 1.5\nDIAM 8\nGLAS N-BK7 1.5168\nCONI -1\nsurface 1\ntype evenasph\nRD -40\nTH 4\nAP 8\nGLASS 1\nCOEF 2 -3e-7\nASPH 9 4\nSURF 2\nTYPE STANDARD\nCURV 0\nDISZ 30\n", LensInterchangeFormat::oslo_len);
    std::vector<double> iv;
    for (auto& s : imp.surfaces()) { iv.push_back(double(int(s.type))); iv.push_back(s.vertex_z); iv.push_back(s.radius); iv.push_back(s.aperture_radius); iv.push_back(s.refractive_index_after); iv.push_back(s.conic_constant); for (double c : s.even_coefficients) iv.push_back(c); }
    vec("import", iv);
    auto js = sequential_system_to_json(sys, 1.33);
    str("json", js);
    auto back = sequential_system_from_json("junk {\"type\":\"sphere\",\"z\":1} " + js);
    std::vector<double> bv{back.object_space_index()};
    for (auto& s : back.surfaces()) { bv.push_back(double(int(s.type))); bv.push_back(s.vertex_z); bv.push_back(s.radius); bv.push_back(s.even_coefficients[3]); }
    vec("json_back", bv);

    NonSequentialScene scene;
    scene.add_plane({{0, 0, 5}, {0, 0.1, 1}, 1.0, 1.5, 0.04, 0.9, 0.01, 0.05});
    scene.add_plane({{0, 0, 10}, {0, 0, -1}, 1.5, 1.0, 0.1, 0.85, 0.05, 0.0});
    NsRay seed; seed.ray.origin = {0.1, 0.2, 0.0}; seed.ray.direction = normalized({0.02, -0.03, 1.0});
    auto out = scene.trace(seed, 4, 1e-4);
    std::vector<double> ns;
    for (auto& r : out) { ns.push_back(r.power); ns.push_back(r.bounces); ns.push_back(r.ray.origin.x); ns.push_back(r.ray.origin.y); ns.push_back(r.ray.origin.z); ns.push_back(r.ray.direction.x); ns.push_back(r.ray.direction.y); ns.push_back(r.ray.direction.z); }
    vec("ns_trace", ns);
    auto sum = summarize_stray_light(out, 1.0);
    vec("ns_summary", std::vector<double>{sum.total_power, sum.reflected_power, sum.transmitted_power, sum.scattered_power, sum.absorbed_power});
    std::vector<double> gh;
    for (auto& g : rank_ghost_paths(out, 2)) { gh.push_back(g.bounces); gh.push_back(g.power); }
    vec("ghosts", gh);

    MeritFunction rosen = [](std::span<const double> x) { return (1 - x[0]) * (1 - x[0]) + 10 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + 0.5 * x[2] * x[2]; };
    auto cd = coordinate_descent(rosen, {-1.2, 1.0, 0.7}, {0.5, 0.5, 0.5}, 300, 1e-8);
    vec("cd", cat(cd.variables, std::vector<double>{cd.merit, double(cd.iterations), double(cd.converged)}));
    std::vector<double> lo{-2, -1, -1}, hi{2, 3, 1};
    auto de = differential_evolution(rosen, lo, hi, 12, 40, 0.7, 0.85, 11);
    vec("de", cat(de.variables, std::vector<double>{de.merit, double(de.iterations)}));
    std::vector<double> p0{0.3, -0.4, 2.0};
    vec("sens", finite_difference_sensitivity(rosen, p0, 1e-6));
    std::vector<double> sig{0.01, 0.02, 0.0};
    auto mc = monte_carlo_tolerance(rosen, p0, sig, 25, 3);
    vec("mc", std::vector<double>{mc.mean, mc.standard_deviation, mc.minimum, mc.maximum, double(mc.samples)});
    vec("autodiff", std::vector<double>{autodiff_derivative([](auto x) { return square(x) * x / (x + decltype(x){2.0}) - x; }, 1.3)});
    std::cout << "\n}\n";
}
