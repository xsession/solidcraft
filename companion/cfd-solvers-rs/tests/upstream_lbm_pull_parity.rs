//! Value-by-value comparison with the C++ upstream in-place (Esoteric-Pull), two-lattice pull,
//! 16-bit compressed and precision-generic LBM solvers. Goldens come from
//! `parity/oracle_lbm_pull.cpp` built against the pinned upstream revision.
use cfd_solvers::lbm::*;
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-lbm-pull.json")).unwrap()
}

fn number(v: &Value) -> f64 {
    match v.as_str() {
        Some("inf") => f64::INFINITY,
        Some("-inf") => f64::NEG_INFINITY,
        Some(_) => f64::NAN,
        None => v.as_f64().unwrap(),
    }
}

/// `|rust - cpp| <= rel * max(|cpp|, 1e-6 * max|cpp|) + abs`; infinities must match exactly.
#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(number)
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    let floor = 1e-6
        * cpp
            .iter()
            .filter(|c| c.is_finite())
            .fold(0.0_f64, |m, c| m.max(c.abs()));
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        if !c.is_finite() {
            assert!(
                r == c || (r.is_nan() && c.is_nan()),
                "{name}[{i}]: rust {r} vs C++ {c}"
            );
            continue;
        }
        assert!(
            (r - c).abs() <= rel * c.abs().max(floor) + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn macro_f32(m: &MacroscopicFields) -> Vec<f64> {
    (0..m.rho.len())
        .flat_map(|i| [m.rho[i], m.ux[i], m.uy[i], m.uz[i]])
        .map(f64::from)
        .collect()
}

fn macro_typed<R: Copy + Into<f64>>(m: &TypedMacroscopicFields<R>) -> Vec<f64> {
    (0..m.rho.len())
        .flat_map(|i| [m.rho[i], m.ux[i], m.uy[i], m.uz[i]])
        .map(Into::into)
        .collect()
}

fn cfg(nx: usize, ny: usize, nz: usize, tau: f32, a: [f32; 3], les: bool) -> InPlaceLbmConfig {
    InPlaceLbmConfig {
        nx,
        ny,
        nz,
        tau,
        acceleration_x: a[0],
        acceleration_y: a[1],
        acceleration_z: a[2],
        smagorinsky_les: les,
    }
}

/// Single-precision solvers: identical operation order gives agreement to a few float ulps.
const F32: f64 = 2e-6;

#[test]
fn d2q9_in_place_two_lattice_and_trt_match_upstream() {
    let c2 = cfg(10, 7, 1, 0.7, [2e-5, -1e-5, 0.0], false);
    let mut e2 = D2Q9InPlaceSolver::new(c2).unwrap();
    e2.initialize_taylor_green(0.04);
    e2.step(9);
    let mut p2 = D2Q9PullSolver::new(c2, 0.0).unwrap();
    p2.initialize_taylor_green(0.04);
    p2.step(9);
    close("e2_macro", &macro_f32(&e2.compute_macroscopic()), F32, 0.0);
    close(
        "e2_raw",
        &e2.raw_storage()
            .iter()
            .map(|v| f64::from(*v))
            .collect::<Vec<_>>(),
        F32,
        0.0,
    );
    close(
        "e2_metrics",
        &[
            e2.mass(),
            e2.kinetic_energy(),
            f64::from(e2.max_speed()),
            e2.time_step() as f64,
        ],
        1e-6,
        0.0,
    );
    close("p2_macro", &macro_f32(&p2.compute_macroscopic()), F32, 0.0);
    close("p2_mass", &[p2.mass()], 1e-9, 0.0);
    let mut trt = D2Q9PullSolver::new(c2, 0.9).unwrap();
    trt.initialize_taylor_green(0.04);
    trt.step(6);
    close(
        "trt_macro",
        &macro_f32(&trt.compute_macroscopic()),
        F32,
        0.0,
    );
}

#[test]
fn d3q19_solids_local_forcing_and_les_match_upstream() {
    let mut e3 = D3Q19Solver::new(cfg(6, 5, 4, 0.62, [1e-5, 0.0, 0.0], true)).unwrap();
    let cells = e3.cells();
    let mut solid = vec![0u8; cells];
    for z in 1..3 {
        for x in 2..4 {
            solid[(z * 5 + 2) * 6 + x] = 1;
        }
    }
    let ax: Vec<f32> = (0..cells).map(|n| 1e-6 * (n % 6) as f32).collect();
    let ay = vec![-2e-6f32; cells];
    let az: Vec<f32> = (0..cells).map(|n| 1e-6 * (n % 4) as f32).collect();
    e3.set_solid_mask(&solid).unwrap();
    e3.set_local_body_acceleration(&ax, &ay, &az).unwrap();
    e3.initialize_taylor_green(0.05);
    e3.step(7);
    close("e3_macro", &macro_f32(&e3.compute_macroscopic()), F32, 0.0);
    close(
        "e3_logical",
        &e3.logical_populations()
            .iter()
            .map(|v| f64::from(*v))
            .collect::<Vec<_>>(),
        F32,
        0.0,
    );
    close(
        "e3_metrics",
        &[e3.mass(), e3.kinetic_energy(), f64::from(e3.max_speed())],
        1e-6,
        0.0,
    );
    e3.clear_solid_mask();
    e3.step(3);
    close(
        "e3_after_clear",
        &[e3.mass(), e3.kinetic_energy(), e3.time_step() as f64],
        1e-6,
        0.0,
    );
}

#[test]
fn d3q27_in_place_and_two_lattice_match_upstream() {
    let c = cfg(4, 4, 4, 0.8, [0.0, 0.0, 3e-5], false);
    let mut e = D3Q27Solver::new(c).unwrap();
    e.initialize_taylor_green(0.03);
    e.step(5);
    let mut p = D3Q27PullSolver::new(c, 0.0).unwrap();
    p.initialize_taylor_green(0.03);
    p.step(5);
    close("e27_macro", &macro_f32(&e.compute_macroscopic()), F32, 0.0);
    close("p27_macro", &macro_f32(&p.compute_macroscopic()), F32, 0.0);
}

#[test]
// The sample literals are the upstream driver's f32 literals verbatim.
#[allow(clippy::excessive_precision)]
fn population_codecs_match_upstream() {
    let samples = [
        0.0f32, -0.0, 1.0, 0.4444444, -0.027777, 6.1e-5, 5.9e-6, 1e-9, 65504.0, 70000.0, -1.5e-3,
        0.33333334, 2.0e-8, 0.99951172, 0.99975586,
    ];
    let (mut half, mut bf, mut shifted, mut packed) = (vec![], vec![], vec![], vec![]);
    type S = ShiftedHalfCodec<D2Q9InPlaceDescriptor>;
    type P = PackedDeviation16Codec<D2Q9InPlaceDescriptor>;
    for (i, &v) in samples.iter().enumerate() {
        let d = i % 9;
        let h = IeeeHalfCodec::encode(v, 0);
        half.extend([f64::from(h), f64::from(IeeeHalfCodec::decode(h, 0))]);
        let b = BFloat16Codec::encode(v, 0);
        bf.extend([f64::from(b), f64::from(BFloat16Codec::decode(b, 0))]);
        let s = S::encode(v, d);
        shifted.extend([f64::from(s), f64::from(S::decode(s, d))]);
        let p = P::encode(v, d);
        packed.extend([f64::from(p), f64::from(P::decode(p, d))]);
    }
    // Encoded words are integers and compare exactly; decoded floats allow one ulp of JSON
    // parsing error (serde_json is built without `float_roundtrip`).
    for (name, values) in [
        ("codec_half", &half),
        ("codec_bf16", &bf),
        ("codec_shifted", &shifted),
        ("codec_packed", &packed),
    ] {
        close(name, values, 4e-16, 0.0);
        let words: Vec<f64> = values.iter().step_by(2).copied().collect();
        let golden: Vec<f64> = goldens()[name]
            .as_array()
            .unwrap()
            .iter()
            .step_by(2)
            .map(number)
            .collect();
        assert_eq!(words, golden, "{name} encoded words");
    }
}

#[test]
fn compressed_and_precision_solvers_match_upstream() {
    let c2 = CompressedLbmConfig {
        nx: 12,
        ny: 8,
        nz: 1,
        tau: 0.75,
    };
    let mut ch = CompressedPopulationPullSolver::<
        D2Q9InPlaceDescriptor,
        ShiftedHalfCodec<D2Q9InPlaceDescriptor>,
    >::new(c2)
    .unwrap();
    ch.initialize_taylor_green(0.04);
    ch.step(10);
    close(
        "comp_shifted",
        &macro_f32(&ch.compute_macroscopic()),
        F32,
        0.0,
    );
    let mut cp = CompressedPopulationPullSolver::<
        D2Q9InPlaceDescriptor,
        PackedDeviation16Codec<D2Q9InPlaceDescriptor>,
    >::new(c2)
    .unwrap();
    cp.initialize_taylor_green(0.04);
    cp.step(10);
    close(
        "comp_packed",
        &macro_f32(&cp.compute_macroscopic()),
        F32,
        0.0,
    );
    let c3 = CompressedLbmConfig {
        nx: 5,
        ny: 4,
        nz: 3,
        tau: 0.7,
    };
    let mut cb = CompressedPopulationPullSolver::<D3Q19Descriptor, BFloat16Codec>::new(c3).unwrap();
    cb.initialize_taylor_green(0.04);
    cb.step(4);
    close("comp_bf16", &macro_f32(&cb.compute_macroscopic()), F32, 0.0);
    let mut ci = CompressedPopulationPullSolver::<D3Q19Descriptor, IeeeHalfCodec>::new(c3).unwrap();
    ci.initialize_taylor_green(0.04);
    ci.step(4);
    close("comp_half", &macro_f32(&ci.compute_macroscopic()), F32, 0.0);
    close(
        "comp_mass",
        &[
            ch.mass(),
            cp.mass(),
            cb.mass(),
            ci.mass(),
            ci.population_bytes() as f64,
        ],
        1e-12,
        0.0,
    );

    let mut pd = PrecisionPullSolver::<D2Q9InPlaceDescriptor, f64>::new(TypedLbmConfig {
        nx: 9,
        ny: 6,
        nz: 1,
        tau: 0.66,
    })
    .unwrap();
    pd.initialize_taylor_green(0.05);
    pd.step(12);
    close(
        "prec_f64",
        &macro_typed(&pd.compute_macroscopic()),
        1e-13,
        0.0,
    );
    let mut pf = PrecisionPullSolver::<D3Q19Descriptor, f32>::new(TypedLbmConfig {
        nx: 4,
        ny: 5,
        nz: 3,
        tau: 0.66,
    })
    .unwrap();
    pf.initialize_taylor_green(0.05);
    pf.step(5);
    close(
        "prec_f32",
        &macro_typed(&pf.compute_macroscopic()),
        F32,
        0.0,
    );
    // Upstream accumulates mass in long double; f64 differs only in the last bits.
    close("prec_mass", &[pd.mass(), pf.mass()], 1e-14, 0.0);
}
