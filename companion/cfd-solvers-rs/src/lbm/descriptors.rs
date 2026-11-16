//! In-place lattice descriptors (D2Q9 in-place ordering, D3Q19, D3Q27), ported from upstream
//! `cfd/solvers/lbm/descriptors` (MIT, cfd_solvers contributors). Non-rest directions come in
//! opposite pairs `(2k-1, 2k)`, which the in-place (esoteric) streaming schemes rely on.

/// A discrete velocity set with weights. `opposite` follows the paired ordering.
pub trait LatticeDescriptor {
    const DIMENSIONS: usize;
    const Q: usize;
    fn cx(i: usize) -> i32;
    fn cy(i: usize) -> i32;
    fn cz(i: usize) -> i32;
    fn weight(i: usize) -> f32;
    fn opposite(i: usize) -> usize {
        if i == 0 {
            0
        } else if i % 2 == 1 {
            i + 1
        } else {
            i - 1
        }
    }
}

macro_rules! descriptor {
    ($name:ident, $doc:literal, $d:expr, $q:expr, $cx:expr, $cy:expr, $cz:expr, $w:expr) => {
        #[doc = $doc]
        #[derive(Clone, Copy, Debug, Default)]
        pub struct $name;
        impl $name {
            pub const CX: [i32; $q] = $cx;
            pub const CY: [i32; $q] = $cy;
            pub const CZ: [i32; $q] = $cz;
            pub const WEIGHTS: [f32; $q] = $w;
        }
        impl LatticeDescriptor for $name {
            const DIMENSIONS: usize = $d;
            const Q: usize = $q;
            fn cx(i: usize) -> i32 {
                Self::CX[i]
            }
            fn cy(i: usize) -> i32 {
                Self::CY[i]
            }
            fn cz(i: usize) -> i32 {
                Self::CZ[i]
            }
            fn weight(i: usize) -> f32 {
                Self::WEIGHTS[i]
            }
        }
    };
}

const W9A: f32 = 4.0 / 9.0;
const W9B: f32 = 1.0 / 9.0;
const W9C: f32 = 1.0 / 36.0;
descriptor!(
    D2Q9InPlaceDescriptor,
    "D2Q9 with opposite-paired ordering.",
    2,
    9,
    [0, 1, -1, 0, 0, 1, -1, -1, 1],
    [0, 0, 0, 1, -1, 1, -1, 1, -1],
    [0; 9],
    [W9A, W9B, W9B, W9B, W9B, W9C, W9C, W9C, W9C]
);

const W19A: f32 = 1.0 / 3.0;
const W19B: f32 = 1.0 / 18.0;
const W19C: f32 = 1.0 / 36.0;
descriptor!(
    D3Q19Descriptor,
    "D3Q19 with opposite-paired ordering.",
    3,
    19,
    [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0],
    [0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1],
    [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1],
    [
        W19A, W19B, W19B, W19B, W19B, W19B, W19B, W19C, W19C, W19C, W19C, W19C, W19C, W19C, W19C,
        W19C, W19C, W19C, W19C
    ]
);

const W27A: f32 = 8.0 / 27.0;
const W27B: f32 = 2.0 / 27.0;
const W27C: f32 = 1.0 / 54.0;
const W27D: f32 = 1.0 / 216.0;
descriptor!(
    D3Q27Descriptor,
    "D3Q27 with opposite-paired ordering.",
    3,
    27,
    [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1],
    [0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1],
    [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1],
    [
        W27A, W27B, W27B, W27B, W27B, W27B, W27B, W27C, W27C, W27C, W27C, W27C, W27C, W27C, W27C,
        W27C, W27C, W27C, W27C, W27D, W27D, W27D, W27D, W27D, W27D, W27D, W27D
    ]
);

/// True when the rest direction is first, `Q` is odd and every pair `(2k-1, 2k)` is opposite.
pub fn descriptor_pairs_are_opposites<D: LatticeDescriptor>() -> bool {
    if D::Q < 3 || D::Q % 2 == 0 || D::cx(0) != 0 || D::cy(0) != 0 || D::cz(0) != 0 {
        return false;
    }
    (1..D::Q).step_by(2).all(|i| {
        D::cx(i) == -D::cx(i + 1) && D::cy(i) == -D::cy(i + 1) && D::cz(i) == -D::cz(i + 1)
    })
}

/// True when the weights sum to one within 1e-6.
pub fn descriptor_weights_sum_to_one<D: LatticeDescriptor>() -> bool {
    let sum: f64 = (0..D::Q).map(|i| f64::from(D::weight(i))).sum();
    sum > 0.999999 && sum < 1.000001
}
