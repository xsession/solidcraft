//! Bit-exact replicas of the libstdc++ random facilities used by upstream units, so seeded
//! stochastic algorithms reproduce the C++ sequences: `std::mt19937_64`,
//! `std::generate_canonical<double, 53>` (`uniform_real_distribution`), the Lemire
//! `uniform_int_distribution<std::size_t>` of GCC 11+, and the Marsaglia-polar
//! `normal_distribution<double>` with its cached second variate.

/// `std::mt19937_64`.
#[derive(Clone, Debug)]
pub struct Mt19937_64 {
    state: Box<[u64; 312]>,
    index: usize,
}

impl Mt19937_64 {
    pub const DEFAULT_SEED: u64 = 5489;

    pub fn new(seed: u64) -> Self {
        let mut state = Box::new([0u64; 312]);
        state[0] = seed;
        for i in 1..312 {
            let prev = state[i - 1];
            state[i] = 6364136223846793005u64
                .wrapping_mul(prev ^ (prev >> 62))
                .wrapping_add(i as u64);
        }
        Self { state, index: 312 }
    }

    fn twist(&mut self) {
        const UPPER: u64 = 0xFFFF_FFFF_8000_0000;
        const LOWER: u64 = 0x7FFF_FFFF;
        const MATRIX: u64 = 0xB502_6F5A_A966_19E9;
        for i in 0..312 {
            let y = (self.state[i] & UPPER) | (self.state[(i + 1) % 312] & LOWER);
            let mut v = self.state[(i + 156) % 312] ^ (y >> 1);
            if y & 1 != 0 {
                v ^= MATRIX;
            }
            self.state[i] = v;
        }
        self.index = 0;
    }

    pub fn next_u64(&mut self) -> u64 {
        if self.index >= 312 {
            self.twist();
        }
        let mut z = self.state[self.index];
        self.index += 1;
        z ^= (z >> 29) & 0x5555_5555_5555_5555;
        z ^= (z << 17) & 0x71D6_7FFF_EDA6_0000;
        z ^= (z << 37) & 0xFFF7_EEE0_0000_0000;
        z ^= z >> 43;
        z
    }

    /// `std::generate_canonical<double, 53>`: one 64-bit draw scaled into [0, 1).
    pub fn canonical(&mut self) -> f64 {
        let r = self.next_u64() as f64 / 18446744073709551616.0;
        if r >= 1.0 {
            1.0 - f64::EPSILON / 2.0
        } else {
            r
        }
    }

    /// `std::uniform_real_distribution<double>(a, b)`.
    pub fn uniform(&mut self, a: f64, b: f64) -> f64 {
        self.canonical() * (b - a) + a
    }

    /// `std::uniform_int_distribution<std::size_t>(a, b)` (inclusive), libstdc++ 11+ algorithm.
    pub fn uniform_index(&mut self, a: usize, b: usize) -> usize {
        let range = (b as u64).wrapping_sub(a as u64);
        if range == u64::MAX {
            return self.next_u64() as usize;
        }
        let erange = range + 1;
        let mut product = u128::from(self.next_u64()) * u128::from(erange);
        let mut low = product as u64;
        if low < erange {
            let threshold = erange.wrapping_neg() % erange;
            while low < threshold {
                product = u128::from(self.next_u64()) * u128::from(erange);
                low = product as u64;
            }
        }
        a + (product >> 64) as usize
    }
}

/// `std::normal_distribution<double>` (Marsaglia polar method, second variate cached).
#[derive(Clone, Copy, Debug)]
pub struct StdNormal {
    mean: f64,
    stddev: f64,
    saved: Option<f64>,
}

impl Default for StdNormal {
    fn default() -> Self {
        Self::new(0.0, 1.0)
    }
}

impl StdNormal {
    pub fn new(mean: f64, stddev: f64) -> Self {
        Self {
            mean,
            stddev,
            saved: None,
        }
    }
    pub fn sample(&mut self, rng: &mut Mt19937_64) -> f64 {
        let ret = if let Some(s) = self.saved.take() {
            s
        } else {
            let (mut x, mut y, mut r2);
            loop {
                x = 2.0 * rng.canonical() - 1.0;
                y = 2.0 * rng.canonical() - 1.0;
                r2 = x * x + y * y;
                if r2 <= 1.0 && r2 != 0.0 {
                    break;
                }
            }
            let mult = (-2.0 * r2.ln() / r2).sqrt();
            self.saved = Some(x * mult);
            y * mult
        };
        ret * self.stddev + self.mean
    }
}

/// Longest prefix of `s` (after leading ASCII whitespace) that reads as a decimal float, as
/// `std::stod` would consume it; returns the value and the bytes consumed including whitespace.
pub(crate) fn parse_float_prefix(s: &str) -> Option<(f64, usize)> {
    let bytes = s.as_bytes();
    let mut end = 0;
    while end < bytes.len() && bytes[end].is_ascii_whitespace() {
        end += 1;
    }
    let start = end;
    if end < bytes.len() && (bytes[end] == b'+' || bytes[end] == b'-') {
        end += 1;
    }
    let mut digits = 0;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
        digits += 1;
    }
    if end < bytes.len() && bytes[end] == b'.' {
        end += 1;
        while end < bytes.len() && bytes[end].is_ascii_digit() {
            end += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return None;
    }
    if end < bytes.len() && (bytes[end] == b'e' || bytes[end] == b'E') {
        let mut e = end + 1;
        if e < bytes.len() && (bytes[e] == b'+' || bytes[e] == b'-') {
            e += 1;
        }
        let first = e;
        while e < bytes.len() && bytes[e].is_ascii_digit() {
            e += 1;
        }
        if e > first {
            end = e;
        }
    }
    s[start..end].parse().ok().map(|v| (v, end))
}
