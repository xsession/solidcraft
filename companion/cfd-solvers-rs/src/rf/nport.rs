//! N-port network algebra: dense complex matrices, real and power-wave Z/S conversion,
//! Touchstone SnP I/O, reference-plane shifts, stability, gain, noise, load-pull, passivity,
//! broadband causality, adaptive sweeps, fixed-pole rational fitting, mixed-mode transforms and
//! L-matching, ported from upstream `cfd/rf/nport` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::Matrix2C;
use crate::{
    core::parse_float_prefix, lbm::format_general, require, Complex64 as C, Result, SolverError,
};
use std::{
    f64::consts::PI,
    ops::{Index, IndexMut},
};

const ZERO: C = C::new(0.0, 0.0);
const ONE: C = C::new(1.0, 0.0);

/// Dense row-major square complex matrix.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ComplexMatrix {
    size: usize,
    data: Vec<C>,
}

impl ComplexMatrix {
    pub fn new(size: usize) -> Self {
        Self::filled(size, ZERO)
    }
    pub fn filled(size: usize, value: C) -> Self {
        Self {
            size,
            data: vec![value; size * size],
        }
    }
    pub fn identity(size: usize) -> Self {
        let mut m = Self::new(size);
        for i in 0..size {
            m[(i, i)] = ONE;
        }
        m
    }
    /// Builds from row-major entries; the length must be a perfect square.
    pub fn from_row_major(data: Vec<C>) -> Result<Self> {
        let size = (data.len() as f64).sqrt().round() as usize;
        require(size * size == data.len(), "RF matrix data is not square")?;
        Ok(Self { size, data })
    }
    pub fn size(&self) -> usize {
        self.size
    }
    pub fn data(&self) -> &[C] {
        &self.data
    }
    /// Bounds-checked element access.
    pub fn get(&self, row: usize, column: usize) -> Result<C> {
        require(
            row < self.size && column < self.size,
            "RF matrix index out of range",
        )?;
        Ok(self.data[row * self.size + column])
    }
}

impl Index<(usize, usize)> for ComplexMatrix {
    type Output = C;
    fn index(&self, (r, c): (usize, usize)) -> &C {
        &self.data[r * self.size + c]
    }
}

impl IndexMut<(usize, usize)> for ComplexMatrix {
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut C {
        &mut self.data[r * self.size + c]
    }
}

pub fn multiply(lhs: &ComplexMatrix, rhs: &ComplexMatrix) -> Result<ComplexMatrix> {
    require(lhs.size == rhs.size, "RF matrix multiply size mismatch")?;
    let n = lhs.size;
    let mut out = ComplexMatrix::new(n);
    for row in 0..n {
        for k in 0..n {
            let a = lhs[(row, k)];
            for col in 0..n {
                out[(row, col)] += a * rhs[(k, col)];
            }
        }
    }
    Ok(out)
}

/// Gauss-Jordan inverse with partial pivoting.
pub fn inverse(matrix: &ComplexMatrix) -> Result<ComplexMatrix> {
    let n = matrix.size;
    require(n > 0, "cannot invert empty RF matrix")?;
    let stride = 2 * n;
    let mut aug = vec![ZERO; n * stride];
    for row in 0..n {
        for col in 0..n {
            aug[row * stride + col] = matrix[(row, col)];
        }
        aug[row * stride + n + row] = ONE;
    }
    for col in 0..n {
        let mut pivot = col;
        let mut best = aug[col * stride + col].norm();
        for row in col + 1..n {
            let v = aug[row * stride + col].norm();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best.is_nan() || best < 1e-30 {
            return Err(SolverError::Singular("singular RF N-port matrix".into()));
        }
        if pivot != col {
            for j in 0..stride {
                aug.swap(pivot * stride + j, col * stride + j);
            }
        }
        let d = aug[col * stride + col];
        for j in 0..stride {
            aug[col * stride + j] /= d;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = aug[row * stride + col];
            if f.norm() == 0.0 {
                continue;
            }
            for j in 0..stride {
                let v = aug[col * stride + j];
                aug[row * stride + j] -= f * v;
            }
        }
    }
    let mut out = ComplexMatrix::new(n);
    for row in 0..n {
        for col in 0..n {
            out[(row, col)] = aug[row * stride + n + col];
        }
    }
    Ok(out)
}

fn validate_reference(reference: &[f64], n: usize) -> Result<()> {
    require(
        reference.len() == n,
        "N-port reference-impedance size mismatch",
    )?;
    require(
        reference.iter().all(|z| *z > 0.0 && z.is_finite()),
        "N-port reference impedance must be finite and positive",
    )
}

fn validate_power_reference(reference: &[C], n: usize) -> Result<()> {
    require(
        reference.len() == n,
        "N-port complex reference-impedance size mismatch",
    )?;
    require(
        reference
            .iter()
            .all(|z| z.re > 0.0 && z.re.is_finite() && z.im.is_finite()),
        "power-wave reference impedance must be finite with positive real part",
    )
}

fn diagonal(values: impl Iterator<Item = C>, n: usize) -> ComplexMatrix {
    let mut m = ComplexMatrix::new(n);
    for (i, v) in values.enumerate() {
        m[(i, i)] = v;
    }
    m
}

fn diag_sqrt(reference: &[f64], reciprocal: bool) -> ComplexMatrix {
    diagonal(
        reference
            .iter()
            .map(|z| C::new(if reciprocal { 1.0 / z.sqrt() } else { z.sqrt() }, 0.0)),
        reference.len(),
    )
}

fn identity_minus(m: &ComplexMatrix) -> ComplexMatrix {
    let mut out = ComplexMatrix::identity(m.size);
    for (o, v) in out.data.iter_mut().zip(&m.data) {
        *o -= v;
    }
    out
}

/// N-port `S = (z - I)(z + I)^-1` (upstream `z_to_s` overload; the name avoids the two-port
/// [`z_to_s`]) with `z = Z0^-1/2 Z Z0^-1/2`.
pub fn nport_z_to_s(z: &ComplexMatrix, reference: &[f64]) -> Result<ComplexMatrix> {
    let n = z.size;
    validate_reference(reference, n)?;
    let inv_root = diag_sqrt(reference, true);
    let normalized = multiply(&multiply(&inv_root, z)?, &inv_root)?;
    let (mut minus, mut plus) = (normalized.clone(), normalized);
    for i in 0..n {
        minus[(i, i)] -= ONE;
        plus[(i, i)] += ONE;
    }
    multiply(&minus, &inverse(&plus)?)
}

pub fn nport_s_to_z(s: &ComplexMatrix, reference: &[f64]) -> Result<ComplexMatrix> {
    let n = s.size;
    validate_reference(reference, n)?;
    let root = diag_sqrt(reference, false);
    let mut plus = s.clone();
    for i in 0..n {
        plus[(i, i)] += ONE;
    }
    let normalized = multiply(&plus, &inverse(&identity_minus(s))?)?;
    multiply(&multiply(&root, &normalized)?, &root)
}

pub fn renormalize_s(
    s: &ComplexMatrix,
    old_reference: &[f64],
    new_reference: &[f64],
) -> Result<ComplexMatrix> {
    nport_z_to_s(&nport_s_to_z(s, old_reference)?, new_reference)
}

/// Kurokawa power waves: `S = R^-1/2 (Z - Z0*) (Z + Z0)^-1 R^1/2`.
pub fn z_to_s_power_wave(z: &ComplexMatrix, reference: &[C]) -> Result<ComplexMatrix> {
    let n = z.size;
    validate_power_reference(reference, n)?;
    let r_root = diagonal(reference.iter().map(|z| C::new(z.re.sqrt(), 0.0)), n);
    let r_inv_root = diagonal(reference.iter().map(|z| C::new(1.0 / z.re.sqrt(), 0.0)), n);
    let (mut num, mut den) = (z.clone(), z.clone());
    for i in 0..n {
        num[(i, i)] -= reference[i].conj();
        den[(i, i)] += reference[i];
    }
    multiply(
        &multiply(&multiply(&r_inv_root, &num)?, &inverse(&den)?)?,
        &r_root,
    )
}

pub fn s_to_z_power_wave(s: &ComplexMatrix, reference: &[C]) -> Result<ComplexMatrix> {
    let n = s.size;
    validate_power_reference(reference, n)?;
    let r_root = diagonal(reference.iter().map(|z| C::new(z.re.sqrt(), 0.0)), n);
    let r_inv_root = diagonal(reference.iter().map(|z| C::new(1.0 / z.re.sqrt(), 0.0)), n);
    let z0 = diagonal(reference.iter().copied(), n);
    let t = multiply(&multiply(&r_root, s)?, &r_inv_root)?;
    let mut rhs = multiply(&t, &z0)?;
    for i in 0..n {
        rhs[(i, i)] += reference[i].conj();
    }
    multiply(&inverse(&identity_minus(&t))?, &rhs)
}

pub fn renormalize_s_power_wave(
    s: &ComplexMatrix,
    old_reference: &[C],
    new_reference: &[C],
) -> Result<ComplexMatrix> {
    z_to_s_power_wave(&s_to_z_power_wave(s, old_reference)?, new_reference)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct NPortPoint {
    pub frequency_hz: f64,
    pub s: ComplexMatrix,
    pub reference_impedance: Vec<f64>,
}

fn frequency_scale(unit: &str) -> Result<f64> {
    match unit.to_ascii_uppercase().as_str() {
        "HZ" => Ok(1.0),
        "KHZ" => Ok(1e3),
        "MHZ" => Ok(1e6),
        "GHZ" => Ok(1e9),
        _ => Err(SolverError::InvalidInput(
            "unsupported Touchstone frequency unit".into(),
        )),
    }
}

fn decode_pair(first: f64, second: f64, format: &str) -> Result<C> {
    let radians = second * PI / 180.0;
    match format.to_ascii_uppercase().as_str() {
        "RI" => Ok(C::new(first, second)),
        "MA" => Ok(C::from_polar(first, radians)),
        "DB" => Ok(C::from_polar(10f64.powf(first / 20.0), radians)),
        _ => Err(SolverError::InvalidInput(
            "unsupported Touchstone data format".into(),
        )),
    }
}

/// `while (stream >> double)`: whole numeric tokens, stopping after a token whose numeric prefix
/// is followed by other characters (or that has no numeric prefix).
fn parse_numbers(line: &str) -> Vec<f64> {
    let line = line.split('!').next().unwrap_or("");
    let mut out = Vec::new();
    for token in line.split_ascii_whitespace() {
        match parse_float_prefix(token) {
            Some((v, used)) => {
                out.push(v);
                if used != token.len() {
                    break;
                }
            }
            None => break,
        }
    }
    out
}

/// Touchstone v1 option lines and the common v2 keywords (`[Number of Ports]`, `[Reference]`,
/// `[Network Data]`, `[End]`). Data are column-major (S11, S21, ..., SNN). Without a port count
/// the data are read as S2P.
pub fn read_touchstone(input: &str, port_count_hint: usize) -> Result<Vec<NPortPoint>> {
    let (mut unit, mut format) = (String::from("GHZ"), String::from("MA"));
    let mut ports = port_count_hint;
    let mut scalar_reference = 50.0;
    let mut explicit_reference: Vec<f64> = Vec::new();
    let mut pending: Vec<f64> = Vec::new();
    let mut points = Vec::new();
    for raw in input.lines() {
        let line = raw
            .split('!')
            .next()
            .unwrap_or("")
            .trim_start_matches([' ', '\t', '\r', '\n']);
        if line.trim_end_matches([' ', '\t', '\r', '\n']).is_empty() {
            continue;
        }
        if let Some(options) = line.strip_prefix('#') {
            let mut tokens = options.split_ascii_whitespace();
            if let Some(u) = tokens.next() {
                unit = u.into();
            }
            let parameter = tokens.next().unwrap_or("");
            if let Some(f) = tokens.next() {
                format = f.into();
            }
            require(
                parameter.eq_ignore_ascii_case("S"),
                "Touchstone N-port reader supports S parameters only",
            )?;
            while let Some(t) = tokens.next() {
                if t.eq_ignore_ascii_case("R") {
                    if let Some(v) = tokens
                        .next()
                        .and_then(|v| parse_float_prefix(v).map(|p| p.0))
                    {
                        scalar_reference = v;
                    }
                }
            }
            continue;
        }
        if let Some(body) = line.strip_prefix('[') {
            let close = body
                .find(']')
                .ok_or_else(|| SolverError::InvalidInput("malformed Touchstone keyword".into()))?;
            let keyword = body[..close].to_ascii_uppercase();
            let rest = &body[close + 1..];
            match keyword.as_str() {
                // A failed integer extraction stores 0, as `std::istream` does.
                "NUMBER OF PORTS" => {
                    ports = rest
                        .split_ascii_whitespace()
                        .next()
                        .and_then(|t| t.parse().ok())
                        .unwrap_or(0)
                }
                "REFERENCE" => explicit_reference = parse_numbers(rest),
                "END" => break,
                _ => {}
            }
            continue;
        }
        pending.extend(parse_numbers(line));
        if ports == 0 {
            ports = 2;
        }
        require(ports <= 1024, "Touchstone port count too large")?;
        let count = 1 + 2 * ports * ports;
        while pending.len() >= count {
            let mut point = NPortPoint {
                frequency_hz: pending[0] * frequency_scale(&unit)?,
                s: ComplexMatrix::new(ports),
                reference_impedance: vec![scalar_reference; ports],
            };
            if explicit_reference.len() == 1 {
                point.reference_impedance = vec![explicit_reference[0]; ports];
            } else if !explicit_reference.is_empty() {
                require(
                    explicit_reference.len() == ports,
                    "Touchstone [Reference] count mismatch",
                )?;
                point.reference_impedance = explicit_reference.clone();
            }
            let mut cursor = 1;
            for col in 0..ports {
                for row in 0..ports {
                    point.s[(row, col)] =
                        decode_pair(pending[cursor], pending[cursor + 1], &format)?;
                    cursor += 2;
                }
            }
            points.push(point);
            pending.drain(..count);
        }
    }
    require(pending.is_empty(), "truncated Touchstone SnP data")?;
    Ok(points)
}

/// Writes RI data in Hz with six significant digits (default `std::ostream` precision), as v2
/// (keywords) or v1 (single common reference).
pub fn write_touchstone(points: &[NPortPoint], version2: bool) -> Result<String> {
    let Some(first) = points.first() else {
        return Ok(String::new());
    };
    let ports = first.s.size();
    require(ports > 0, "cannot write empty-port Touchstone network")?;
    for p in points {
        require(
            p.s.size() == ports && p.reference_impedance.len() == ports,
            "Touchstone point topology mismatch",
        )?;
        validate_reference(&p.reference_impedance, ports)?;
    }
    let g = |v: f64| format_general(v, 6);
    let mut o = String::new();
    if version2 {
        o += &format!("[Version] 2.0\n[Number of Ports] {ports}\n[Reference]");
        for z in &first.reference_impedance {
            o += &format!(" {}", g(*z));
        }
        o += "\n# HZ S RI R 50\n[Network Data]\n";
    } else {
        let z0 = first.reference_impedance[0];
        require(
            first
                .reference_impedance
                .iter()
                .all(|v| (v - z0).abs() <= 1e-12),
            "Touchstone v1 writer requires one common reference impedance",
        )?;
        o += &format!("# HZ S RI R {}\n", g(z0));
    }
    for p in points {
        o += &g(p.frequency_hz);
        for col in 0..ports {
            for row in 0..ports {
                let v = p.s[(row, col)];
                o += &format!(" {} {}", g(v.re), g(v.im));
            }
        }
        o.push('\n');
    }
    if version2 {
        o += "[End]\n";
    }
    Ok(o)
}

/// `S'_ij = S_ij exp(gamma_i d_i) exp(gamma_j d_j)`: positive distances de-embed line length.
pub fn shift_reference_planes(
    s: &ComplexMatrix,
    gamma: &[C],
    distance: &[f64],
) -> Result<ComplexMatrix> {
    let n = s.size;
    require(
        gamma.len() == n && distance.len() == n,
        "reference-plane shift size mismatch",
    )?;
    require(
        gamma
            .iter()
            .zip(distance)
            .all(|(g, d)| d.is_finite() && g.re.is_finite() && g.im.is_finite()),
        "invalid reference-plane shift",
    )?;
    let factor: Vec<C> = gamma
        .iter()
        .zip(distance)
        .map(|(g, d)| (g * d).exp())
        .collect();
    let mut out = ComplexMatrix::new(n);
    for r in 0..n {
        for c in 0..n {
            out[(r, c)] = s[(r, c)] * factor[r] * factor[c];
        }
    }
    Ok(out)
}

pub fn return_loss_db(reflection: C) -> Result<f64> {
    let m = reflection.norm();
    require(m.is_finite(), "invalid reflection coefficient")?;
    Ok(if m == 0.0 {
        f64::INFINITY
    } else {
        -20.0 * m.log10()
    })
}

pub fn vswr(reflection: C) -> Result<f64> {
    let m = reflection.norm();
    require(m.is_finite(), "invalid reflection coefficient")?;
    Ok(if m >= 1.0 {
        f64::INFINITY
    } else {
        (1.0 + m) / (1.0 - m)
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StabilityMetrics {
    pub rollett_k: f64,
    pub mu: f64,
    pub determinant: C,
}

/// Rollett K, Edwards-Sinsky mu and `det S`.
pub fn stability_metrics(s: &Matrix2C) -> StabilityMetrics {
    let delta = s.a11 * s.a22 - s.a12 * s.a21;
    let den = 2.0 * (s.a12 * s.a21).norm();
    let k = if den > 0.0 {
        (1.0 - s.a11.norm_sqr() - s.a22.norm_sqr() + delta.norm_sqr()) / den
    } else {
        f64::INFINITY
    };
    let mu_den = (s.a22 - s.a11.conj() * delta).norm() + (s.a12 * s.a21).norm();
    let mu = if mu_den > 0.0 {
        (1.0 - s.a11.norm_sqr()) / mu_den
    } else {
        f64::INFINITY
    };
    StabilityMetrics {
        rollett_k: k,
        mu,
        determinant: delta,
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StabilityCircle {
    pub center: C,
    pub radius: f64,
    pub valid: bool,
}

fn stability_circle(sxx: C, syy: C, s: &Matrix2C) -> StabilityCircle {
    let delta = s.a11 * s.a22 - s.a12 * s.a21;
    let den = sxx.norm_sqr() - delta.norm_sqr();
    if den.abs() < 1e-30 {
        return StabilityCircle {
            center: ZERO,
            radius: f64::INFINITY,
            valid: false,
        };
    }
    let center = (sxx - delta * syy.conj()).conj() / den;
    let radius = (s.a12 * s.a21).norm() / den.abs();
    StabilityCircle {
        center,
        radius,
        valid: radius.is_finite() && center.re.is_finite() && center.im.is_finite(),
    }
}

pub fn source_stability_circle(s: &Matrix2C) -> StabilityCircle {
    stability_circle(s.a11, s.a22, s)
}

pub fn load_stability_circle(s: &Matrix2C) -> StabilityCircle {
    stability_circle(s.a22, s.a11, s)
}

/// Transducer power gain for source/load reflections inside the unit disc.
pub fn transducer_gain(s: &Matrix2C, source_reflection: C, load_reflection: C) -> Result<f64> {
    let (gs, gl) = (source_reflection, load_reflection);
    require(
        gs.norm() < 1.0 && gl.norm() < 1.0,
        "source/load reflection magnitude must be < 1",
    )?;
    let den = (ONE - s.a11 * gs) * (ONE - s.a22 * gl) - s.a12 * s.a21 * gs * gl;
    let d2 = den.norm_sqr();
    Ok(if d2 < 1e-30 {
        f64::INFINITY
    } else {
        (1.0 - gs.norm_sqr()) * s.a21.norm_sqr() * (1.0 - gl.norm_sqr()) / d2
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NoiseParameters {
    pub minimum_noise_factor: f64,
    pub optimum_source_reflection: C,
    pub equivalent_noise_resistance_ohm: f64,
    pub reference_impedance_ohm: f64,
}

impl Default for NoiseParameters {
    fn default() -> Self {
        Self {
            minimum_noise_factor: 1.0,
            optimum_source_reflection: ZERO,
            equivalent_noise_resistance_ohm: 0.0,
            reference_impedance_ohm: 50.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NoiseCircle {
    pub center: C,
    pub radius: f64,
    pub noise_factor: f64,
    pub valid: bool,
}

pub fn noise_factor(p: &NoiseParameters, source_reflection: C) -> Result<f64> {
    require(
        p.minimum_noise_factor >= 1.0
            && p.equivalent_noise_resistance_ohm >= 0.0
            && p.reference_impedance_ohm > 0.0
            && p.minimum_noise_factor.is_finite()
            && p.equivalent_noise_resistance_ohm.is_finite()
            && p.reference_impedance_ohm.is_finite()
            && p.optimum_source_reflection.norm() < 1.0
            && source_reflection.norm() < 1.0,
        "invalid RF noise parameters/reflection coefficient",
    )?;
    let den = (1.0 - source_reflection.norm_sqr()) * (ONE + p.optimum_source_reflection).norm_sqr();
    let coefficient = 4.0 * p.equivalent_noise_resistance_ohm / p.reference_impedance_ohm;
    Ok(p.minimum_noise_factor
        + coefficient * (source_reflection - p.optimum_source_reflection).norm_sqr() / den)
}

pub fn noise_circle(p: &NoiseParameters, target_noise_factor: f64) -> Result<NoiseCircle> {
    require(
        target_noise_factor >= p.minimum_noise_factor
            && target_noise_factor.is_finite()
            && p.equivalent_noise_resistance_ohm > 0.0
            && p.reference_impedance_ohm > 0.0
            && p.optimum_source_reflection.norm() < 1.0,
        "invalid RF noise-circle controls",
    )?;
    let coefficient = 4.0 * p.equivalent_noise_resistance_ohm / p.reference_impedance_ohm;
    let n = (target_noise_factor - p.minimum_noise_factor)
        * (ONE + p.optimum_source_reflection).norm_sqr()
        / coefficient;
    let den = 1.0 + n;
    let center = p.optimum_source_reflection / den;
    let radius = (n * (den - p.optimum_source_reflection.norm_sqr()) / (den * den))
        .max(0.0)
        .sqrt();
    Ok(NoiseCircle {
        center,
        radius,
        noise_factor: target_noise_factor,
        valid: radius.is_finite() && center.norm() + radius < 1.0 + 1e-12,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LoadPullSample {
    pub load_reflection: C,
    pub transducer_gain: f64,
}

/// Polar load-pull grid (origin first; upstream defaults 10, 72, 0.95).
pub fn sample_load_pull(
    s: &Matrix2C,
    source_reflection: C,
    radial_steps: usize,
    angular_steps: usize,
    maximum_radius: f64,
) -> Result<Vec<LoadPullSample>> {
    require(
        radial_steps > 0 && angular_steps > 0 && (0.0..1.0).contains(&maximum_radius),
        "invalid load-pull sampling controls",
    )?;
    let mut out = vec![LoadPullSample {
        load_reflection: ZERO,
        transducer_gain: transducer_gain(s, source_reflection, ZERO)?,
    }];
    for r in 1..=radial_steps {
        let radius = maximum_radius * r as f64 / radial_steps as f64;
        for a in 0..angular_steps {
            let gamma = C::from_polar(radius, 2.0 * PI * a as f64 / angular_steps as f64);
            out.push(LoadPullSample {
                load_reflection: gamma,
                transducer_gain: transducer_gain(s, source_reflection, gamma)?,
            });
        }
    }
    Ok(out)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkQuality {
    pub maximum_singular_value: f64,
    pub reciprocity_error: f64,
    pub passive: bool,
    pub reciprocal: bool,
}

/// Largest singular value by power iteration on `S^H S` and relative `||S - S^T||_F` (upstream
/// default tolerances 1e-10).
pub fn network_quality(
    s: &ComplexMatrix,
    passivity_tolerance: f64,
    reciprocity_tolerance: f64,
) -> Result<NetworkQuality> {
    let n = s.size;
    require(
        n > 0 && passivity_tolerance >= 0.0 && reciprocity_tolerance >= 0.0,
        "invalid network-quality controls",
    )?;
    let mut x = vec![C::new(1.0 / (n as f64).sqrt(), 0.0); n];
    let mut lambda = 0.0_f64;
    for _ in 0..80 {
        let mut y = vec![ZERO; n];
        let mut z = vec![ZERO; n];
        for r in 0..n {
            for c in 0..n {
                y[r] += s[(r, c)] * x[c];
            }
        }
        for c in 0..n {
            for r in 0..n {
                z[c] += s[(r, c)].conj() * y[r];
            }
        }
        let norm = z.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
        if norm < 1e-30 {
            lambda = 0.0;
            break;
        }
        for i in 0..n {
            x[i] = z[i] / norm;
        }
        let mut rq = ZERO;
        for i in 0..n {
            rq += x[i].conj() * z[i];
        }
        let next = rq.re.max(0.0);
        if (next - lambda).abs() < 1e-14 * next.max(1.0) {
            lambda = next;
            break;
        }
        lambda = next;
    }
    let sigma = lambda.max(0.0).sqrt();
    let (mut diff2, mut scale2) = (0.0, 0.0);
    for r in 0..n {
        for c in 0..n {
            diff2 += (s[(r, c)] - s[(c, r)]).norm_sqr();
            scale2 += s[(r, c)].norm_sqr();
        }
    }
    let reciprocity = diff2.sqrt() / scale2.sqrt().max(1e-30);
    Ok(NetworkQuality {
        maximum_singular_value: sigma,
        reciprocity_error: reciprocity,
        passive: sigma <= 1.0 + passivity_tolerance,
        reciprocal: reciprocity <= reciprocity_tolerance,
    })
}

/// Uniform scaling down to `maximum_singular_value` (upstream default 1) when exceeded.
pub fn enforce_passivity(s: &ComplexMatrix, maximum_singular_value: f64) -> Result<ComplexMatrix> {
    require(
        s.size > 0 && maximum_singular_value > 0.0 && maximum_singular_value.is_finite(),
        "invalid passivity-enforcement controls",
    )?;
    let sigma = network_quality(s, 1e-10, 1e-10)?.maximum_singular_value;
    require(sigma.is_finite(), "non-finite network singular value")?;
    if sigma <= maximum_singular_value {
        return Ok(s.clone());
    }
    let scale = maximum_singular_value / sigma;
    let mut out = s.clone();
    out.data.iter_mut().for_each(|v| *v *= scale);
    Ok(out)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BroadbandCausality {
    pub negative_time_energy_ratio: f64,
    pub causal: bool,
    pub impulse_samples: usize,
}

/// Hermitian-extended inverse DFT of a uniform DC..Nyquist sweep; reports the energy fraction in
/// wrapped negative-time samples (upstream defaults 1e-6 and 1e-8).
pub fn check_broadband_causality(
    points: &[NPortPoint],
    negative_time_energy_tolerance: f64,
    frequency_uniformity_tolerance: f64,
) -> Result<BroadbandCausality> {
    let (tol, ftol) = (
        negative_time_energy_tolerance,
        frequency_uniformity_tolerance,
    );
    require(
        points.len() >= 3 && tol >= 0.0 && tol.is_finite() && ftol >= 0.0 && ftol.is_finite(),
        "invalid broadband causality controls",
    )?;
    let ports = points[0].s.size();
    require(ports > 0, "causality check requires non-empty networks")?;
    require(
        points[0].frequency_hz.abs() <= ftol,
        "causality check requires a DC first sample",
    )?;
    let df = points[1].frequency_hz - points[0].frequency_hz;
    require(
        df > 0.0 && df.is_finite(),
        "causality check requires increasing frequencies",
    )?;
    for (k, p) in points.iter().enumerate() {
        require(p.s.size() == ports, "causality sweep port-count mismatch")?;
        let expected = df * k as f64;
        let scale = 1.0f64.max(expected.abs()).max(p.frequency_hz.abs());
        require(
            (p.frequency_hz - expected).abs() <= ftol * scale,
            "causality check requires uniformly spaced frequencies",
        )?;
    }
    let size = 2 * (points.len() - 1);
    let first_negative = points.len();
    let (mut negative, mut total) = (0.0, 0.0);
    let mut spectrum = vec![ZERO; size];
    for row in 0..ports {
        for col in 0..ports {
            spectrum.fill(ZERO);
            for (k, p) in points.iter().enumerate() {
                spectrum[k] = p.s[(row, col)];
            }
            spectrum[0] = C::new(spectrum[0].re, 0.0);
            let last = points.len() - 1;
            spectrum[last] = C::new(spectrum[last].re, 0.0);
            for k in 1..last {
                spectrum[size - k] = spectrum[k].conj();
            }
            for n in 0..size {
                let mut value = ZERO;
                for (k, sk) in spectrum.iter().enumerate() {
                    let phase = 2.0 * PI * ((k * n) as f64) / size as f64;
                    value += sk * C::new(0.0, phase).exp();
                }
                let energy = (value / size as f64).norm_sqr();
                total += energy;
                if n >= first_negative {
                    negative += energy;
                }
            }
        }
    }
    let ratio = if total > 0.0 { negative / total } else { 0.0 };
    Ok(BroadbandCausality {
        negative_time_energy_ratio: ratio,
        causal: ratio <= tol,
        impulse_samples: size,
    })
}

/// `H(s) = D + sum_k R_k (-p_k) / (s - p_k)` with stable real poles.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RationalNetworkModel {
    pub poles_rad_per_s: Vec<f64>,
    pub residues: Vec<ComplexMatrix>,
    pub direct: ComplexMatrix,
    pub reference_impedance: Vec<f64>,
}

impl RationalNetworkModel {
    pub fn evaluate(&self, frequency_hz: f64) -> Result<ComplexMatrix> {
        require(
            frequency_hz >= 0.0
                && self.residues.len() == self.poles_rad_per_s.len()
                && self.direct.size() > 0,
            "invalid rational network model/evaluation frequency",
        )?;
        require(
            self.residues.iter().all(|r| r.size() == self.direct.size()),
            "rational network residue dimension mismatch",
        )?;
        let mut value = self.direct.clone();
        let s = C::new(0.0, 2.0 * PI * frequency_hz);
        for (p, r) in self.poles_rad_per_s.iter().zip(&self.residues) {
            let scale = C::new(-p, 0.0) / (s - p);
            for (v, rv) in value.data.iter_mut().zip(&r.data) {
                *v += rv * scale;
            }
        }
        Ok(value)
    }
}

fn solve_small_complex_system(mut a: Vec<C>, mut b: Vec<C>, n: usize) -> Result<Vec<C>> {
    for k in 0..n {
        let mut pivot = k;
        let mut best = a[k * n + k].norm();
        for row in k + 1..n {
            let c = a[row * n + k].norm();
            if c > best {
                best = c;
                pivot = row;
            }
        }
        if best.is_nan() || best <= 1e-24 {
            return Err(SolverError::Singular(
                "singular rational network fit".into(),
            ));
        }
        if pivot != k {
            for col in k..n {
                a.swap(k * n + col, pivot * n + col);
            }
            b.swap(k, pivot);
        }
        for row in k + 1..n {
            let factor = a[row * n + k] / a[k * n + k];
            if factor.norm() == 0.0 {
                continue;
            }
            a[row * n + k] = ZERO;
            for col in k + 1..n {
                let v = a[k * n + col];
                a[row * n + col] -= factor * v;
            }
            let bk = b[k];
            b[row] -= factor * bk;
        }
    }
    let mut x = vec![ZERO; n];
    for i in (0..n).rev() {
        let mut v = b[i];
        for col in i + 1..n {
            v -= a[i * n + col] * x[col];
        }
        x[i] = v / a[i * n + i];
    }
    Ok(x)
}

/// Fixed log-spaced negative-real poles over the band; per-entry complex least squares (upstream
/// default 6 poles).
pub fn fit_rational_network(
    points: &[NPortPoint],
    pole_count: usize,
) -> Result<RationalNetworkModel> {
    require(
        points.len() >= 3 && pole_count > 0 && pole_count < points.len(),
        "invalid rational network fit controls",
    )?;
    let ports = points[0].s.size();
    require(
        ports > 0 && points[0].reference_impedance.len() == ports,
        "invalid rational network sweep",
    )?;
    let (mut fmin, mut fmax) = (points[0].frequency_hz, points[0].frequency_hz);
    for p in points {
        require(
            p.s.size() == ports
                && p.reference_impedance == points[0].reference_impedance
                && p.frequency_hz > 0.0,
            "inconsistent rational network sweep",
        )?;
        fmin = fmin.min(p.frequency_hz);
        fmax = fmax.max(p.frequency_hz);
    }
    require(
        fmax > fmin,
        "rational network fit requires a frequency band",
    )?;
    let ratio = (fmax / fmin).powf(1.0 / (pole_count - 1).max(1) as f64);
    let poles: Vec<f64> = (0..pole_count)
        .map(|k| -2.0 * PI * fmin * ratio.powf(k as f64))
        .collect();
    let columns = pole_count + 1;
    let design: Vec<Vec<C>> = points
        .iter()
        .map(|p| {
            let s = C::new(0.0, 2.0 * PI * p.frequency_hz);
            std::iter::once(ONE)
                .chain(poles.iter().map(|pk| C::new(-pk, 0.0) / (s - pk)))
                .collect()
        })
        .collect();
    let mut model = RationalNetworkModel {
        poles_rad_per_s: poles,
        residues: vec![ComplexMatrix::new(ports); pole_count],
        direct: ComplexMatrix::new(ports),
        reference_impedance: points[0].reference_impedance.clone(),
    };
    for row in 0..ports {
        for col in 0..ports {
            let mut normal = vec![ZERO; columns * columns];
            let mut rhs = vec![ZERO; columns];
            for (i, p) in points.iter().enumerate() {
                for a in 0..columns {
                    rhs[a] += design[i][a].conj() * p.s[(row, col)];
                    for b in 0..columns {
                        normal[a * columns + b] += design[i][a].conj() * design[i][b];
                    }
                }
            }
            let coef = solve_small_complex_system(normal, rhs, columns)?;
            model.direct[(row, col)] = coef[0];
            for k in 0..pole_count {
                model.residues[k][(row, col)] = coef[k + 1];
            }
        }
    }
    Ok(model)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveNetworkSweepConfig {
    pub initial_points: usize,
    pub maximum_points: usize,
    pub maximum_refinements: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
}

impl Default for AdaptiveNetworkSweepConfig {
    fn default() -> Self {
        Self {
            initial_points: 5,
            maximum_points: 257,
            maximum_refinements: 12,
            relative_tolerance: 1e-3,
            absolute_tolerance: 1e-6,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AdaptiveNetworkSweepResult {
    pub points: Vec<NPortPoint>,
    pub evaluations: usize,
    pub refinements: usize,
    pub converged: bool,
}

/// Log-frequency sweep refined at geometric midpoints whose linear-interpolation error exceeds
/// the mixed tolerance.
pub fn adaptive_network_sweep(
    start_hz: f64,
    stop_hz: f64,
    reference_impedance: &[f64],
    sampler: &dyn Fn(f64) -> ComplexMatrix,
    config: &AdaptiveNetworkSweepConfig,
) -> Result<AdaptiveNetworkSweepResult> {
    require(
        start_hz > 0.0
            && stop_hz > start_hz
            && stop_hz.is_finite()
            && !reference_impedance.is_empty()
            && config.initial_points >= 2
            && config.maximum_points >= config.initial_points
            && config.maximum_refinements > 0
            && config.relative_tolerance > 0.0
            && config.absolute_tolerance > 0.0,
        "invalid adaptive network sweep controls",
    )?;
    require(
        reference_impedance
            .iter()
            .all(|z| *z > 0.0 && z.is_finite()),
        "invalid network reference impedance",
    )?;
    let mut result = AdaptiveNetworkSweepResult::default();
    let sample = |f: f64, evaluations: &mut usize| -> Result<NPortPoint> {
        let s = sampler(f);
        *evaluations += 1;
        require(
            s.size() == reference_impedance.len(),
            "adaptive network sampler port-count mismatch",
        )?;
        Ok(NPortPoint {
            frequency_hz: f,
            s,
            reference_impedance: reference_impedance.to_vec(),
        })
    };
    let ratio = (stop_hz / start_hz).powf(1.0 / (config.initial_points - 1) as f64);
    for i in 0..config.initial_points {
        let f = if i + 1 == config.initial_points {
            stop_hz
        } else {
            start_hz * ratio.powf(i as f64)
        };
        let p = sample(f, &mut result.evaluations)?;
        result.points.push(p);
    }
    for _ in 0..config.maximum_refinements {
        let mut additions = Vec::new();
        for i in 0..result.points.len() - 1 {
            if result.points.len() + additions.len() >= config.maximum_points {
                break;
            }
            let (left, right) = (&result.points[i], &result.points[i + 1]);
            let mid = (left.frequency_hz * right.frequency_hz).sqrt();
            let actual = sample(mid, &mut result.evaluations)?;
            let mut error = 0.0_f64;
            for k in 0..actual.s.data.len() {
                let predicted = (left.s.data[k] + right.s.data[k]) * 0.5;
                let scale = config.absolute_tolerance
                    + config.relative_tolerance
                        * actual.s.data[k]
                            .norm()
                            .max(left.s.data[k].norm())
                            .max(right.s.data[k].norm());
                error = error.max((actual.s.data[k] - predicted).norm() / scale);
            }
            if error > 1.0 {
                additions.push(actual);
            }
        }
        if additions.is_empty() {
            result.converged = true;
            break;
        }
        result.refinements += 1;
        result.points.extend(additions);
        result
            .points
            .sort_by(|a, b| a.frequency_hz.total_cmp(&b.frequency_hz));
        if result.points.len() >= config.maximum_points {
            break;
        }
    }
    Ok(result)
}

/// Single-ended 4-port to `[d1, d2, c1, c2]` mixed-mode S for pairs (1,2) and (3,4).
pub fn single_ended_to_mixed_mode(four_port_s: &ComplexMatrix) -> Result<ComplexMatrix> {
    require(
        four_port_s.size() == 4,
        "mixed-mode conversion requires a 4-port network",
    )?;
    let h = C::new(1.0 / 2f64.sqrt(), 0.0);
    let mut c = ComplexMatrix::new(4);
    c[(0, 0)] = h;
    c[(0, 1)] = -h;
    c[(1, 2)] = h;
    c[(1, 3)] = -h;
    c[(2, 0)] = h;
    c[(2, 1)] = h;
    c[(3, 2)] = h;
    c[(3, 3)] = h;
    let mut ct = ComplexMatrix::new(4);
    for r in 0..4 {
        for col in 0..4 {
            ct[(r, col)] = c[(col, r)];
        }
    }
    multiply(&multiply(&c, four_port_s)?, &ct)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LMatch {
    pub series_reactance_ohm: f64,
    pub shunt_susceptance_siemens: f64,
    pub shunt_first: bool,
}

/// Low-pass L network between two positive resistances.
pub fn synthesize_l_match(source_resistance: f64, load_resistance: f64) -> Result<LMatch> {
    let (rs, rl) = (source_resistance, load_resistance);
    require(
        rs > 0.0 && rl > 0.0 && rs.is_finite() && rl.is_finite(),
        "L-match resistances must be finite and positive",
    )?;
    if (rs - rl).abs() < 1e-15 * rs.max(rl) {
        return Ok(LMatch::default());
    }
    if rl > rs {
        let q = (rl / rs - 1.0).sqrt();
        return Ok(LMatch {
            series_reactance_ohm: rs * q,
            shunt_susceptance_siemens: q / rl,
            shunt_first: false,
        });
    }
    let q = (rs / rl - 1.0).sqrt();
    Ok(LMatch {
        series_reactance_ohm: rl * q,
        shunt_susceptance_siemens: q / rs,
        shunt_first: true,
    })
}
