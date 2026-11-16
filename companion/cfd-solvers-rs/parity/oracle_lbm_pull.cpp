// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream LBM pull/in-place streaming units and prints JSON goldens.
#include "cfd/solvers/lbm/compressed_pull.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"
#include "cfd/solvers/lbm/one_step_pull.hpp"
#include "cfd/solvers/lbm/precision_pull.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd::lbm;
static bool first = true;
template <class V>
static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const double x = static_cast<double>(v[i]);
        std::cout << (i ? "," : "");
        if (std::isnan(x)) std::cout << "\"nan\"";
        else if (std::isinf(x)) std::cout << (x > 0 ? "\"inf\"" : "\"-inf\"");
        else std::cout << std::setprecision(17) << x;
    }
    std::cout << "]";
}
template <class M> static std::vector<double> macro(const M& m) {
    std::vector<double> o;
    for (std::size_t i = 0; i < m.rho.size(); ++i) { o.push_back(m.rho[i]); o.push_back(m.ux[i]); o.push_back(m.uy[i]); o.push_back(m.uz[i]); }
    return o;
}

int main() {
    std::cout << "{\n";
    // D2Q9: esoteric vs two-lattice with uniform forcing.
    InPlaceLbmConfig c2{10, 7, 1, 0.7F, 2e-5F, -1e-5F, 0.0F, false};
    D2Q9InPlaceSolver e2(c2);
    e2.initialize_taylor_green(0.04F);
    e2.step(9);
    D2Q9PullSolver p2(c2);
    p2.initialize_taylor_green(0.04F);
    p2.step(9);
    vec("e2_macro", macro(e2.compute_macroscopic()));
    { std::vector<float> raw; for (std::size_t d = 0; d < 9; ++d) for (std::size_t n = 0; n < e2.cells(); ++n) raw.push_back(e2.raw_storage()(d, n)); vec("e2_raw", raw); }
    vec("e2_metrics", std::vector<double>{e2.mass(), e2.kinetic_energy(), e2.max_speed(), double(e2.time_step())});
    vec("p2_macro", macro(p2.compute_macroscopic()));
    vec("p2_mass", std::vector<double>{p2.mass()});
    D2Q9PullSolver trt(c2, 0.9F);
    trt.initialize_taylor_green(0.04F);
    trt.step(6);
    vec("trt_macro", macro(trt.compute_macroscopic()));

    // D3Q19 with solids, local acceleration and Smagorinsky LES.
    InPlaceLbmConfig c3{6, 5, 4, 0.62F, 1e-5F, 0.0F, 0.0F, true};
    D3Q19Solver e3(c3);
    std::vector<std::uint8_t> solid(e3.cells(), 0);
    for (std::size_t z = 1; z < 3; ++z) for (std::size_t y = 2; y < 3; ++y) for (std::size_t x = 2; x < 4; ++x) solid[(z * 5 + y) * 6 + x] = 1;
    std::vector<float> ax(e3.cells()), ay(e3.cells()), az(e3.cells());
    for (std::size_t n = 0; n < e3.cells(); ++n) { ax[n] = 1e-6F * float(n % 6); ay[n] = -2e-6F; az[n] = 1e-6F * float(n % 4); }
    e3.set_solid_mask(solid);
    e3.set_local_body_acceleration(ax, ay, az);
    e3.initialize_taylor_green(0.05F);
    e3.step(7);
    vec("e3_macro", macro(e3.compute_macroscopic()));
    vec("e3_logical", e3.logical_populations());
    vec("e3_metrics", std::vector<double>{e3.mass(), e3.kinetic_energy(), e3.max_speed()});
    e3.clear_solid_mask();
    e3.step(3);
    vec("e3_after_clear", std::vector<double>{e3.mass(), e3.kinetic_energy(), double(e3.time_step())});

    // D3Q27 plain and two-lattice reference.
    InPlaceLbmConfig c27{4, 4, 4, 0.8F, 0.0F, 0.0F, 3e-5F, false};
    D3Q27Solver e27(c27);
    e27.initialize_taylor_green(0.03F);
    e27.step(5);
    D3Q27PullSolver p27(c27);
    p27.initialize_taylor_green(0.03F);
    p27.step(5);
    vec("e27_macro", macro(e27.compute_macroscopic()));
    vec("p27_macro", macro(p27.compute_macroscopic()));

    // Codecs.
    std::vector<float> samples{0.0F, -0.0F, 1.0F, 0.4444444F, -0.027777F, 6.1e-5F, 5.9e-6F, 1e-9F, 65504.0F, 70000.0F, -1.5e-3F, 0.33333334F, 2.0e-8F, 0.99951172F, 0.99975586F};
    std::vector<double> half, bf, shifted, packed;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const float v = samples[i];
        const int d = int(i % 9);
        half.push_back(IeeeHalfCodec::encode(v)); half.push_back(IeeeHalfCodec::decode(IeeeHalfCodec::encode(v)));
        bf.push_back(BFloat16Codec::encode(v)); bf.push_back(BFloat16Codec::decode(BFloat16Codec::encode(v)));
        shifted.push_back(ShiftedHalfCodec<D2Q9InPlaceDescriptor>::encode(v, d)); shifted.push_back(ShiftedHalfCodec<D2Q9InPlaceDescriptor>::decode(ShiftedHalfCodec<D2Q9InPlaceDescriptor>::encode(v, d), d));
        packed.push_back(PackedDeviation16Codec<D2Q9InPlaceDescriptor>::encode(v, d)); packed.push_back(PackedDeviation16Codec<D2Q9InPlaceDescriptor>::decode(PackedDeviation16Codec<D2Q9InPlaceDescriptor>::encode(v, d), d));
    }
    vec("codec_half", half);
    vec("codec_bf16", bf);
    vec("codec_shifted", shifted);
    vec("codec_packed", packed);

    // Compressed solvers.
    CompressedPopulationPullSolver<D2Q9InPlaceDescriptor, ShiftedHalfCodec<D2Q9InPlaceDescriptor>> ch({12, 8, 1, 0.75F});
    ch.initialize_taylor_green(0.04F); ch.step(10);
    vec("comp_shifted", macro(ch.compute_macroscopic()));
    CompressedPopulationPullSolver<D2Q9InPlaceDescriptor, PackedDeviation16Codec<D2Q9InPlaceDescriptor>> cp({12, 8, 1, 0.75F});
    cp.initialize_taylor_green(0.04F); cp.step(10);
    vec("comp_packed", macro(cp.compute_macroscopic()));
    CompressedPopulationPullSolver<D3Q19Descriptor, BFloat16Codec> cb({5, 4, 3, 0.7F});
    cb.initialize_taylor_green(0.04F); cb.step(4);
    vec("comp_bf16", macro(cb.compute_macroscopic()));
    CompressedPopulationPullSolver<D3Q19Descriptor, IeeeHalfCodec> ci({5, 4, 3, 0.7F});
    ci.initialize_taylor_green(0.04F); ci.step(4);
    vec("comp_half", macro(ci.compute_macroscopic()));
    vec("comp_mass", std::vector<double>{ch.mass(), cp.mass(), cb.mass(), ci.mass(), double(ci.population_bytes())});

    // Precision-generic.
    PrecisionPullSolver<D2Q9InPlaceDescriptor, double> pd({9, 6, 1, 0.66});
    pd.initialize_taylor_green(0.05); pd.step(12);
    vec("prec_f64", macro(pd.compute_macroscopic()));
    PrecisionPullSolver<D3Q19Descriptor, float> pf({4, 5, 3, 0.66F});
    pf.initialize_taylor_green(0.05F); pf.step(5);
    vec("prec_f32", macro(pf.compute_macroscopic()));
    vec("prec_mass", std::vector<double>{double(pd.mass()), double(pf.mass())});
    std::cout << "\n}\n";
}
