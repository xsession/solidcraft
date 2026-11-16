//! Rectangular-waveguide port modes, multiconductor cable telegrapher solves and shield transfer
//! impedance, normal-incidence GSTC metasurface sheets, 1-D driven/eigen frequency-domain
//! Maxwell, and signal/power-integrity utilities (NRZ eye, decoupling, PDN IR drop, pulses),
//! ported from upstream `cfd/em/{wave_port, cable, metasurface_gstc, frequency_domain, system}`
//! (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use crate::{
    core::erfc,
    require,
    rf::{inverse, multiply, ComplexMatrix},
    Complex64 as C, Result, SolverError,
};
use std::f64::consts::PI;

const EPS0: f64 = 8.8541878128e-12;
const MU0: f64 = 1.25663706212e-6;
const ZERO: C = C::new(0.0, 0.0);

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub enum ModeFamily {
    #[default]
    Te,
    Tm,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RectangularWaveguideMode {
    pub family: ModeFamily,
    pub m: usize,
    pub n: usize,
    pub frequency_hz: f64,
    pub cutoff_frequency_hz: f64,
    /// `beta` for `exp(-j beta z)`; `-j alpha` below cutoff.
    pub propagation_constant_rad_per_m: C,
    pub wave_impedance_ohm: C,
    pub propagating: bool,
    /// Transverse-E scale; propagating modes carry exactly 1 W forward.
    pub electric_scale: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransverseModeFields {
    pub electric_v_per_m: [C; 2],
    pub magnetic_a_per_m: [C; 2],
}

fn cos2_integral(length: f64, index: usize) -> f64 {
    if index == 0 {
        length
    } else {
        0.5 * length
    }
}
fn sin2_integral(length: f64, index: usize) -> f64 {
    if index == 0 {
        0.0
    } else {
        0.5 * length
    }
}

/// Analytic PEC rectangular-guide TE/TM modes sorted by cutoff (TE before TM, then `(m, n)`);
/// upstream defaults `eps_r = mu_r = 1`, `mode_count = 4`.
pub fn rectangular_waveguide_modes(
    width_m: f64,
    height_m: f64,
    frequency_hz: f64,
    relative_permittivity: f64,
    relative_permeability: f64,
    mode_count: usize,
) -> Result<Vec<RectangularWaveguideMode>> {
    let (a, b) = (width_m, height_m);
    require(
        a > 0.0
            && b > 0.0
            && frequency_hz > 0.0
            && relative_permittivity > 0.0
            && relative_permeability > 0.0
            && mode_count > 0
            && mode_count < 10_000,
        "invalid rectangular wave-port controls",
    )?;
    let (eps, mu) = (EPS0 * relative_permittivity, MU0 * relative_permeability);
    let omega = 2.0 * PI * frequency_hz;
    let k = omega * (mu * eps).sqrt();
    let velocity = 1.0 / (mu * eps).sqrt();
    let limit = (mode_count + 3).max(4);
    let mut modes = Vec::new();
    for m in 0..=limit {
        for n in 0..=limit {
            if m == 0 && n == 0 {
                continue;
            }
            let (kx, ky) = (PI * m as f64 / a, PI * n as f64 / b);
            let kc = kx.hypot(ky);
            let fc = velocity * kc / (2.0 * PI);
            let propagating = k > kc;
            let mut te = RectangularWaveguideMode {
                family: ModeFamily::Te,
                m,
                n,
                frequency_hz,
                cutoff_frequency_hz: fc,
                propagating,
                ..Default::default()
            };
            if propagating {
                let beta = (k * k - kc * kc).sqrt();
                te.propagation_constant_rad_per_m = C::new(beta, 0.0);
                te.wave_impedance_ohm = C::new(omega * mu / beta, 0.0);
            } else {
                let alpha = (kc * kc - k * k).sqrt();
                te.propagation_constant_rad_per_m = C::new(0.0, -alpha);
                te.wave_impedance_ohm = C::new(0.0, -omega * mu / alpha);
            }
            let integral = ky * ky * cos2_integral(a, m) * sin2_integral(b, n)
                + kx * kx * sin2_integral(a, m) * cos2_integral(b, n);
            te.electric_scale = if propagating {
                (2.0 * te.wave_impedance_ohm.re / integral).sqrt()
            } else {
                1.0 / integral.sqrt()
            };
            modes.push(te);
            if m > 0 && n > 0 {
                let mut tm = RectangularWaveguideMode {
                    family: ModeFamily::Tm,
                    ..te
                };
                if propagating {
                    let beta = (k * k - kc * kc).sqrt();
                    tm.propagation_constant_rad_per_m = C::new(beta, 0.0);
                    tm.wave_impedance_ohm = C::new(beta / (omega * eps), 0.0);
                } else {
                    let alpha = (kc * kc - k * k).sqrt();
                    tm.propagation_constant_rad_per_m = C::new(0.0, -alpha);
                    tm.wave_impedance_ohm = C::new(0.0, alpha / (omega * eps));
                }
                let ti = kx * kx * cos2_integral(a, m) * sin2_integral(b, n)
                    + ky * ky * sin2_integral(a, m) * cos2_integral(b, n);
                tm.electric_scale = if propagating {
                    (2.0 * tm.wave_impedance_ohm.re / ti).sqrt()
                } else {
                    1.0 / ti.sqrt()
                };
                modes.push(tm);
            }
        }
    }
    modes.sort_by(|l, r| {
        l.cutoff_frequency_hz
            .total_cmp(&r.cutoff_frequency_hz)
            .then(l.family.cmp(&r.family))
            .then((l.m, l.n).cmp(&(r.m, r.n)))
    });
    modes.truncate(mode_count);
    Ok(modes)
}

/// Transverse fields of a mode at `(x, y)` inside the guide.
pub fn rectangular_waveguide_transverse_fields(
    mode: &RectangularWaveguideMode,
    width_m: f64,
    height_m: f64,
    x_m: f64,
    y_m: f64,
) -> Result<TransverseModeFields> {
    let (a, b) = (width_m, height_m);
    require(
        a > 0.0 && b > 0.0 && (0.0..=a).contains(&x_m) && (0.0..=b).contains(&y_m),
        "invalid rectangular wave-port field coordinate",
    )?;
    let (kx, ky) = (PI * mode.m as f64 / a, PI * mode.n as f64 / b);
    let (ex, ey) = match mode.family {
        ModeFamily::Te => (
            ky * (kx * x_m).cos() * (ky * y_m).sin(),
            -kx * (kx * x_m).sin() * (ky * y_m).cos(),
        ),
        ModeFamily::Tm => (
            kx * (kx * x_m).cos() * (ky * y_m).sin(),
            ky * (kx * x_m).sin() * (ky * y_m).cos(),
        ),
    };
    let e = [
        C::new(mode.electric_scale * ex, 0.0),
        C::new(mode.electric_scale * ey, 0.0),
    ];
    let mut f = TransverseModeFields {
        electric_v_per_m: e,
        magnetic_a_per_m: [ZERO; 2],
    };
    if mode.wave_impedance_ohm.norm() > 0.0 {
        f.magnetic_a_per_m = [
            -e[1] / mode.wave_impedance_ohm,
            e[0] / mode.wave_impedance_ohm,
        ];
    }
    Ok(f)
}

/// Midpoint-rule time-average forward power (upstream defaults 64 x 64).
pub fn rectangular_mode_power_w(
    mode: &RectangularWaveguideMode,
    width_m: f64,
    height_m: f64,
    nx: usize,
    ny: usize,
) -> Result<f64> {
    require(
        nx > 0 && ny > 0,
        "wave-port power quadrature requires non-zero sampling",
    )?;
    let (dx, dy) = (width_m / nx as f64, height_m / ny as f64);
    let mut power = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let f = rectangular_waveguide_transverse_fields(
                mode,
                width_m,
                height_m,
                (i as f64 + 0.5) * dx,
                (j as f64 + 0.5) * dy,
            )?;
            let pz = f.electric_v_per_m[0] * f.magnetic_a_per_m[1].conj()
                - f.electric_v_per_m[1] * f.magnetic_a_per_m[0].conj();
            power += 0.5 * pz.re * dx * dy;
        }
    }
    Ok(power)
}

/// Per-unit-length row-major RLCG matrices `[ohm/m, H/m, S/m, F/m]`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MulticonductorRlcg {
    pub conductors: usize,
    pub resistance: Vec<f64>,
    pub inductance: Vec<f64>,
    pub conductance: Vec<f64>,
    pub capacitance: Vec<f64>,
}

impl MulticonductorRlcg {
    pub fn validate(&self) -> Result<()> {
        require(self.conductors > 0, "cable RLCG requires conductors")?;
        let n = self.conductors * self.conductors;
        require(
            [
                &self.resistance,
                &self.inductance,
                &self.conductance,
                &self.capacitance,
            ]
            .iter()
            .all(|v| v.len() == n),
            "cable RLCG matrix size mismatch",
        )?;
        require(
            self.resistance.iter().all(|v| v.is_finite()),
            "non-finite cable resistance",
        )?;
        require(
            self.inductance.iter().all(|v| v.is_finite()),
            "non-finite cable inductance",
        )?;
        require(
            self.conductance.iter().all(|v| v.is_finite()),
            "non-finite cable conductance",
        )?;
        require(
            self.capacitance.iter().all(|v| v.is_finite()),
            "non-finite cable capacitance",
        )
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CableFrequencyResult {
    pub frequency_hz: f64,
    pub input_voltage_v: Vec<C>,
    pub input_current_a: Vec<C>,
    pub load_voltage_v: Vec<C>,
    pub load_current_a: Vec<C>,
    pub transfer_matrix: ComplexMatrix,
}

fn block(m: &ComplexMatrix, r0: usize, c0: usize, n: usize) -> ComplexMatrix {
    let mut out = ComplexMatrix::new(n);
    for r in 0..n {
        for c in 0..n {
            out[(r, c)] = m[(r0 + r, c0 + c)];
        }
    }
    out
}

fn add_scaled(mut a: ComplexMatrix, b: &ComplexMatrix, scale: f64) -> ComplexMatrix {
    for r in 0..a.size() {
        for c in 0..a.size() {
            a[(r, c)] += b[(r, c)] * scale;
        }
    }
    a
}

fn matvec(a: &ComplexMatrix, x: &[C]) -> Vec<C> {
    let n = a.size();
    let mut y = vec![ZERO; n];
    for r in 0..n {
        for c in 0..n {
            y[r] += a[(r, c)] * x[c];
        }
    }
    y
}

/// Cascaded `[1/1]` Pade steps of the first-order RLCG system along the cable, terminated by
/// the full load matrix `V_L = Z_L I_L`.
pub fn solve_multiconductor_cable(
    length_m: f64,
    segments: usize,
    frequency_hz: f64,
    rlcg: &MulticonductorRlcg,
    imposed_input_voltage_v: &[C],
    load_impedance_ohm: &ComplexMatrix,
) -> Result<CableFrequencyResult> {
    rlcg.validate()?;
    let n = rlcg.conductors;
    require(
        length_m > 0.0
            && segments > 0
            && frequency_hz >= 0.0
            && imposed_input_voltage_v.len() == n
            && load_impedance_ohm.size() == n,
        "invalid multiconductor cable controls",
    )?;
    let omega = 2.0 * PI * frequency_hz;
    let dz = length_m / segments as f64;
    let mut system = ComplexMatrix::new(2 * n);
    for r in 0..n {
        for c in 0..n {
            let k = r * n + c;
            system[(r, n + c)] = -C::new(rlcg.resistance[k], omega * rlcg.inductance[k]);
            system[(n + r, c)] = -C::new(rlcg.conductance[k], omega * rlcg.capacitance[k]);
        }
    }
    let (mut lhs, mut rhs) = (
        ComplexMatrix::identity(2 * n),
        ComplexMatrix::identity(2 * n),
    );
    for r in 0..2 * n {
        for c in 0..2 * n {
            lhs[(r, c)] -= system[(r, c)] * (0.5 * dz);
            rhs[(r, c)] += system[(r, c)] * (0.5 * dz);
        }
    }
    let step = multiply(&inverse(&lhs)?, &rhs)?;
    let mut transfer = ComplexMatrix::identity(2 * n);
    for _ in 0..segments {
        transfer = multiply(&step, &transfer)?;
    }
    let (a, b, c, d) = (
        block(&transfer, 0, 0, n),
        block(&transfer, 0, n, n),
        block(&transfer, n, 0, n),
        block(&transfer, n, n, n),
    );
    let left = add_scaled(b, &multiply(load_impedance_ohm, &d)?, -1.0);
    let right = add_scaled(multiply(load_impedance_ohm, &c)?, &a, -1.0);
    let right_v = matvec(&right, imposed_input_voltage_v);
    let input_current = matvec(&inverse(&left)?, &right_v);
    let state: Vec<C> = imposed_input_voltage_v
        .iter()
        .chain(&input_current)
        .copied()
        .collect();
    let load_state = matvec(&transfer, &state);
    Ok(CableFrequencyResult {
        frequency_hz,
        input_voltage_v: imposed_input_voltage_v.to_vec(),
        input_current_a: input_current,
        load_voltage_v: load_state[..n].to_vec(),
        load_current_a: load_state[n..].to_vec(),
        transfer_matrix: transfer,
    })
}

/// One [`solve_multiconductor_cable`] per frequency with frequency-dependent RLCG.
pub fn sweep_multiconductor_cable(
    length_m: f64,
    segments: usize,
    frequencies_hz: &[f64],
    rlcg_sampler: &dyn Fn(f64) -> MulticonductorRlcg,
    imposed_input_voltage_v: &[C],
    load_impedance_ohm: &ComplexMatrix,
) -> Result<Vec<CableFrequencyResult>> {
    require(!frequencies_hz.is_empty(), "invalid cable frequency sweep")?;
    frequencies_hz
        .iter()
        .map(|f| {
            solve_multiconductor_cable(
                length_m,
                segments,
                *f,
                &rlcg_sampler(*f),
                imposed_input_voltage_v,
                load_impedance_ohm,
            )
        })
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CableShieldTransferModel {
    pub dc_transfer_resistance_ohm_per_m: f64,
    pub transfer_inductance_h_per_m: f64,
    pub skin_corner_hz: f64,
}

impl Default for CableShieldTransferModel {
    fn default() -> Self {
        Self {
            dc_transfer_resistance_ohm_per_m: 0.0,
            transfer_inductance_h_per_m: 0.0,
            skin_corner_hz: 1e6,
        }
    }
}

/// `Z_T = R_dc sqrt(1 + j f / f_c) + j omega L_T`.
pub fn cable_shield_transfer_impedance_ohm_per_m(
    frequency_hz: f64,
    model: &CableShieldTransferModel,
) -> Result<C> {
    require(
        frequency_hz >= 0.0
            && model.dc_transfer_resistance_ohm_per_m >= 0.0
            && model.transfer_inductance_h_per_m >= 0.0
            && model.skin_corner_hz > 0.0,
        "invalid cable shield transfer model",
    )?;
    let skin = C::new(1.0, frequency_hz / model.skin_corner_hz).sqrt();
    Ok(skin * model.dc_transfer_resistance_ohm_per_m
        + C::new(
            0.0,
            2.0 * PI * frequency_hz * model.transfer_inductance_h_per_m,
        ))
}

pub fn cable_shield_induced_voltage_v(
    frequency_hz: f64,
    shield_current_a: C,
    exposed_length_m: f64,
    model: &CableShieldTransferModel,
) -> Result<C> {
    require(exposed_length_m >= 0.0, "invalid cable exposed length")?;
    Ok(shield_current_a
        * cable_shield_transfer_impedance_ohm_per_m(frequency_hz, model)?
        * exposed_length_m)
}

/// `E_parallel h_eff k` (upstream default coupling factor 1).
pub fn field_to_cable_open_circuit_voltage_v(
    parallel_electric_field_v_per_m: C,
    effective_height_m: f64,
    coupling_factor: C,
) -> Result<C> {
    let (e, k) = (parallel_electric_field_v_per_m, coupling_factor);
    require(
        effective_height_m >= 0.0
            && e.re.is_finite()
            && e.im.is_finite()
            && k.re.is_finite()
            && k.im.is_finite(),
        "invalid field-to-cable coupling controls",
    )?;
    Ok(e * effective_height_m * k)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GstcSheetConfig {
    pub electric_surface_admittance_s: C,
    pub magnetic_surface_impedance_ohm: C,
    pub wave_impedance_before_ohm: f64,
    pub wave_impedance_after_ohm: f64,
}

impl Default for GstcSheetConfig {
    fn default() -> Self {
        Self {
            electric_surface_admittance_s: ZERO,
            magnetic_surface_impedance_ohm: ZERO,
            wave_impedance_before_ohm: 376.730313668,
            wave_impedance_after_ohm: 376.730313668,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GstcScattering {
    pub reflection: C,
    pub transmission: C,
    pub reflected_power: f64,
    pub transmitted_power: f64,
    pub absorbed_power: f64,
}

/// Normal-incidence scalar generalized sheet transition conditions with averaged sheet fields.
pub fn solve_gstc_normal_incidence(config: &GstcSheetConfig) -> Result<GstcScattering> {
    let (ye, zm) = (
        config.electric_surface_admittance_s,
        config.magnetic_surface_impedance_ohm,
    );
    require(
        ye.re.is_finite() && ye.im.is_finite(),
        "GSTC electric admittance must be finite",
    )?;
    require(
        zm.re.is_finite() && zm.im.is_finite(),
        "GSTC magnetic impedance must be finite",
    )?;
    let (eb, ea) = (
        config.wave_impedance_before_ohm,
        config.wave_impedance_after_ohm,
    );
    require(
        eb > 0.0 && ea > 0.0 && eb.is_finite() && ea.is_finite(),
        "GSTC wave impedances must be finite and positive",
    )?;
    let one = C::new(1.0, 0.0);
    let a11 = -1.0 / eb - ye / 2.0;
    let a12 = -1.0 / ea - ye / 2.0;
    let b1 = -(C::new(1.0 / eb, 0.0) - ye / 2.0);
    let a21 = -one - zm / (2.0 * eb);
    let a22 = one + zm / (2.0 * ea);
    let b2 = one - zm / (2.0 * eb);
    let det = a11 * a22 - a12 * a21;
    require(det.norm() > 1e-30, "singular GSTC sheet system")?;
    let reflection = (b1 * a22 - a12 * b2) / det;
    let transmission = (a11 * b2 - b1 * a21) / det;
    let reflected_power = reflection.norm_sqr();
    let transmitted_power = (eb / ea) * transmission.norm_sqr();
    Ok(GstcScattering {
        reflection,
        transmission,
        reflected_power,
        transmitted_power,
        absorbed_power: 1.0 - reflected_power - transmitted_power,
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrequencyDomain1DConfig {
    pub points: usize,
    pub length_m: f64,
    pub frequency_hz: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub conductivity_s_per_m: f64,
}

impl Default for FrequencyDomain1DConfig {
    fn default() -> Self {
        Self {
            points: 201,
            length_m: 1.0,
            frequency_hz: 1e8,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            conductivity_s_per_m: 0.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct FrequencyDomain1DResult {
    pub electric_v_per_m: Vec<C>,
    /// Staggered between electric samples.
    pub magnetic_a_per_m: Vec<C>,
    pub spacing_m: f64,
}

fn solve_dense_banded(
    mut m: Vec<C>,
    mut rhs: Vec<C>,
    n: usize,
    singular_floor: f64,
    what: &str,
) -> Result<Vec<C>> {
    for k in 0..n {
        let mut pivot = k;
        let mut best = m[k * n + k].norm();
        for row in k + 1..n {
            let v = m[row * n + k].norm();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best.is_nan() || best <= singular_floor {
            return Err(SolverError::Singular(what.into()));
        }
        if pivot != k {
            for col in k..n {
                m.swap(k * n + col, pivot * n + col);
            }
            rhs.swap(k, pivot);
        }
        let diag = m[k * n + k];
        for row in k + 1..n {
            let factor = m[row * n + k] / diag;
            if factor.norm() == 0.0 {
                continue;
            }
            m[row * n + k] = ZERO;
            for col in k + 1..n {
                let v = m[k * n + col];
                m[row * n + col] -= factor * v;
            }
            let r = rhs[k];
            rhs[row] -= factor * r;
        }
    }
    let mut x = vec![ZERO; n];
    for row in (0..n).rev() {
        let mut s = rhs[row];
        for col in row + 1..n {
            s -= m[row * n + col] * x[col];
        }
        x[row] = s / m[row * n + row];
    }
    Ok(x)
}

/// `E'' + k^2 E = j omega mu J` with PEC ends and `exp(+j omega t)` phasors.
pub fn solve_pec_driven_maxwell_1d(
    config: &FrequencyDomain1DConfig,
    impressed_current_density_a_per_m2: &[C],
) -> Result<FrequencyDomain1DResult> {
    let c = config;
    require(
        (3..=20_000).contains(&c.points)
            && c.length_m > 0.0
            && c.frequency_hz > 0.0
            && c.relative_permittivity > 0.0
            && c.relative_permeability > 0.0
            && c.conductivity_s_per_m >= 0.0
            && impressed_current_density_a_per_m2.len() == c.points,
        "invalid driven 1-D Maxwell configuration",
    )?;
    let interior = c.points - 2;
    let dx = c.length_m / (c.points - 1) as f64;
    let omega = 2.0 * PI * c.frequency_hz;
    let (mu, eps) = (
        MU0 * c.relative_permeability,
        EPS0 * c.relative_permittivity,
    );
    let k2 = C::new(
        omega * omega * mu * eps,
        -omega * mu * c.conductivity_s_per_m,
    );
    let inv_dx2 = 1.0 / (dx * dx);
    let mut m = vec![ZERO; interior * interior];
    let mut rhs = vec![ZERO; interior];
    for row in 0..interior {
        m[row * interior + row] = C::new(-2.0 * inv_dx2, 0.0) + k2;
        if row > 0 {
            m[row * interior + row - 1] = C::new(inv_dx2, 0.0);
        }
        if row + 1 < interior {
            m[row * interior + row + 1] = C::new(inv_dx2, 0.0);
        }
        rhs[row] = C::new(0.0, omega * mu) * impressed_current_density_a_per_m2[row + 1];
    }
    let e = solve_dense_banded(
        m,
        rhs,
        interior,
        f64::EPSILON,
        "singular frequency-domain Maxwell matrix",
    )?;
    let mut electric = vec![ZERO; c.points];
    electric[1..=interior].copy_from_slice(&e);
    let h_factor = C::new(0.0, 1.0 / (omega * mu * dx));
    let magnetic = (0..c.points - 1)
        .map(|i| h_factor * (electric[i + 1] - electric[i]))
        .collect();
    Ok(FrequencyDomain1DResult {
        electric_v_per_m: electric,
        magnetic_a_per_m: magnetic,
        spacing_m: dx,
    })
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CavityMode1D {
    pub angular_frequency_rad_per_s: f64,
    pub frequency_hz: f64,
    pub electric_shape: Vec<f64>,
}

fn jacobi_symmetric(a: &mut [f64], n: usize) -> Result<Vec<f64>> {
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _ in 0..80 * n * n {
        let (mut p, mut q, mut largest) = (0, 0, 0.0);
        for i in 0..n {
            for j in i + 1..n {
                let value = a[i * n + j].abs();
                if value > largest {
                    largest = value;
                    p = i;
                    q = j;
                }
            }
        }
        let scale = (0..n).fold(1.0_f64, |m, i| m.max(a[i * n + i].abs()));
        if largest <= 1e-13 * scale {
            return Ok(v);
        }
        let (app, aqq, apq) = (a[p * n + p], a[q * n + q], a[p * n + q]);
        let angle = 0.5 * (2.0 * apq).atan2(aqq - app);
        let (c, s) = (angle.cos(), angle.sin());
        for k in 0..n {
            if k == p || k == q {
                continue;
            }
            let (akp, akq) = (a[k * n + p], a[k * n + q]);
            let (nkp, nkq) = (c * akp - s * akq, s * akp + c * akq);
            a[k * n + p] = nkp;
            a[p * n + k] = nkp;
            a[k * n + q] = nkq;
            a[q * n + k] = nkq;
        }
        a[p * n + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q * n + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p * n + q] = 0.0;
        a[q * n + p] = 0.0;
        for k in 0..n {
            let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
            v[k * n + p] = c * vkp - s * vkq;
            v[k * n + q] = s * vkp + c * vkq;
        }
    }
    Err(SolverError::NumericalFailure(
        "cavity Jacobi eigensolver failed to converge".into(),
    ))
}

/// Lowest PEC cavity modes from the symmetrised generalized eigenproblem (Jacobi); shapes are
/// peak-normalised. `relative_permittivity` is nodal.
pub fn pec_cavity_eigenmodes_1d(
    points: usize,
    length_m: f64,
    relative_permittivity: &[f64],
    relative_permeability: f64,
    mode_count: usize,
) -> Result<Vec<CavityMode1D>> {
    require(
        (3..=2_000).contains(&points)
            && length_m > 0.0
            && relative_permittivity.len() == points
            && relative_permeability > 0.0
            && mode_count > 0,
        "invalid PEC cavity eigenmode configuration",
    )?;
    require(
        relative_permittivity
            .iter()
            .all(|v| *v > 0.0 && v.is_finite()),
        "invalid cavity permittivity",
    )?;
    let n = points - 2;
    require(mode_count <= n, "requested too many cavity modes")?;
    let dx = length_m / (points - 1) as f64;
    let inv_dx2 = 1.0 / (dx * dx);
    let mu = MU0 * relative_permeability;
    let mass: Vec<f64> = (0..n)
        .map(|i| mu * EPS0 * relative_permittivity[i + 1])
        .collect();
    let mut m = vec![0.0; n * n];
    for row in 0..n {
        m[row * n + row] = 2.0 * inv_dx2 / mass[row];
        if row + 1 < n {
            let v = -inv_dx2 / (mass[row] * mass[row + 1]).sqrt();
            m[row * n + row + 1] = v;
            m[(row + 1) * n + row] = v;
        }
    }
    let vectors = jacobi_symmetric(&mut m, n)?;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|a, b| m[a * n + a].total_cmp(&m[b * n + b]));
    let mut modes = Vec::with_capacity(mode_count);
    for &col in order.iter().take(mode_count) {
        let lambda = m[col * n + col];
        require(lambda > 0.0, "non-positive Maxwell cavity eigenvalue")?;
        let w = lambda.sqrt();
        let mut shape = vec![0.0; points];
        let mut peak = 0.0_f64;
        for i in 0..n {
            shape[i + 1] = vectors[i * n + col] / mass[i].sqrt();
            peak = peak.max(shape[i + 1].abs());
        }
        require(peak > 0.0, "empty cavity eigenvector")?;
        shape.iter_mut().for_each(|v| *v /= peak);
        modes.push(CavityMode1D {
            angular_frequency_rad_per_s: w,
            frequency_hz: w / (2.0 * PI),
            electric_shape: shape,
        });
    }
    Ok(modes)
}

/// Uniform-permittivity convenience overload.
pub fn pec_cavity_eigenmodes_1d_uniform(
    points: usize,
    length_m: f64,
    relative_permittivity: f64,
    relative_permeability: f64,
    mode_count: usize,
) -> Result<Vec<CavityMode1D>> {
    require(
        points <= 2_000,
        "invalid PEC cavity eigenmode configuration",
    )?;
    pec_cavity_eigenmodes_1d(
        points,
        length_m,
        &vec![relative_permittivity; points],
        relative_permeability,
        mode_count,
    )
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EyeDiagramMetrics {
    pub symbol_count: usize,
    pub samples_per_symbol: usize,
    pub zero_mean: f64,
    pub one_mean: f64,
    pub zero_sigma: f64,
    pub one_sigma: f64,
    pub decision_threshold: f64,
    pub eye_height: f64,
    pub eye_width_ui: f64,
    pub estimated_ber: f64,
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn sample_sigma(v: &[f64], m: f64) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() - 1) as f64).sqrt()
}

fn qfunc(x: f64) -> f64 {
    0.5 * erfc(x / 2f64.sqrt())
}

/// NRZ eye metrics against a known 0/1 bit stream; `decision_sample = 0` selects the UI centre.
pub fn analyze_nrz_eye(
    waveform_v: &[f64],
    symbols: &[i32],
    samples_per_symbol: usize,
    decision_sample: usize,
) -> Result<EyeDiagramMetrics> {
    require(
        samples_per_symbol >= 2 && !symbols.is_empty(),
        "invalid eye-diagram controls",
    )?;
    let decision = if decision_sample == 0 {
        samples_per_symbol / 2
    } else {
        decision_sample
    };
    require(decision < samples_per_symbol, "invalid eye decision sample")?;
    require(
        symbols
            .len()
            .checked_mul(samples_per_symbol)
            .is_some_and(|n| waveform_v.len() >= n),
        "waveform shorter than bit stream",
    )?;
    let (mut zeros, mut ones) = (Vec::new(), Vec::new());
    for (k, s) in symbols.iter().enumerate() {
        let v = waveform_v[k * samples_per_symbol + decision];
        require(v.is_finite(), "non-finite eye waveform")?;
        match s {
            0 => zeros.push(v),
            1 => ones.push(v),
            _ => return Err(SolverError::InvalidInput("eye symbols must be 0/1".into())),
        }
    }
    require(
        !zeros.is_empty() && !ones.is_empty(),
        "eye analysis requires both zero and one symbols",
    )?;
    let mut r = EyeDiagramMetrics {
        symbol_count: symbols.len(),
        samples_per_symbol,
        ..Default::default()
    };
    r.zero_mean = mean(&zeros);
    r.one_mean = mean(&ones);
    r.zero_sigma = sample_sigma(&zeros, r.zero_mean);
    r.one_sigma = sample_sigma(&ones, r.one_mean);
    r.decision_threshold = 0.5 * (r.zero_mean + r.one_mean);
    r.eye_height = (r.one_mean - 3.0 * r.one_sigma) - (r.zero_mean + 3.0 * r.zero_sigma);
    let (s0, s1) = (r.zero_sigma.max(1e-30), r.one_sigma.max(1e-30));
    r.estimated_ber = 0.5
        * (qfunc((r.one_mean - r.decision_threshold) / s1)
            + qfunc((r.decision_threshold - r.zero_mean) / s0));
    let mut open = 0;
    for phase in 0..samples_per_symbol {
        let (mut max_zero, mut min_one) = (f64::NEG_INFINITY, f64::INFINITY);
        for (k, s) in symbols.iter().enumerate() {
            let v = waveform_v[k * samples_per_symbol + phase];
            if *s == 0 {
                max_zero = max_zero.max(v);
            } else {
                min_one = min_one.min(v);
            }
        }
        if min_one > max_zero {
            open += 1;
        }
    }
    r.eye_width_ui = open as f64 / samples_per_symbol as f64;
    Ok(r)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecouplingCapacitor {
    pub capacitance_f: f64,
    pub esr_ohm: f64,
    pub esl_h: f64,
    pub count: usize,
}

impl Default for DecouplingCapacitor {
    fn default() -> Self {
        Self {
            capacitance_f: 0.0,
            esr_ohm: 0.0,
            esl_h: 0.0,
            count: 1,
        }
    }
}

/// `(ESR + j omega ESL + 1/(j omega C)) / count`; DC gives an infinite reactance.
pub fn decoupling_impedance_ohm(capacitor: &DecouplingCapacitor, frequency_hz: f64) -> Result<C> {
    require(
        frequency_hz >= 0.0 && frequency_hz.is_finite(),
        "invalid frequency",
    )?;
    require(
        capacitor.capacitance_f > 0.0
            && capacitor.esr_ohm >= 0.0
            && capacitor.esl_h >= 0.0
            && capacitor.count > 0,
        "invalid decoupling capacitor",
    )?;
    let omega = 2.0 * PI * frequency_hz;
    let mut z = C::new(capacitor.esr_ohm, omega * capacitor.esl_h);
    z += if omega == 0.0 {
        C::new(0.0, f64::NEG_INFINITY)
    } else {
        C::new(0.0, -1.0 / (omega * capacitor.capacitance_f))
    };
    Ok(z / capacitor.count as f64)
}

/// Parallel combination of a plane impedance (ignored when zero) and capacitor banks.
pub fn pdn_parallel_impedance_ohm(
    capacitors: &[DecouplingCapacitor],
    frequency_hz: f64,
    plane_impedance_ohm: C,
) -> Result<C> {
    require(
        frequency_hz >= 0.0 && frequency_hz.is_finite(),
        "invalid frequency",
    )?;
    let one = C::new(1.0, 0.0);
    let mut y = ZERO;
    if plane_impedance_ohm.norm() > 0.0 {
        y += one / plane_impedance_ohm;
    }
    for cap in capacitors {
        let z = decoupling_impedance_ohm(cap, frequency_hz)?;
        if z.norm() > 0.0 && z.re.is_finite() && z.im.is_finite() {
            y += one / z;
        }
    }
    Ok(if y.norm() == 0.0 {
        C::new(f64::INFINITY, 0.0)
    } else {
        one / y
    })
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DecouplingOptimizationResult {
    pub selected: Vec<DecouplingCapacitor>,
    pub worst_ratio_before: f64,
    pub worst_ratio_after: f64,
    pub iterations: usize,
}

fn worst_ratio(
    frequencies: &[f64],
    targets: &[f64],
    caps: &[DecouplingCapacitor],
    plane: C,
) -> Result<f64> {
    let mut worst = 0.0_f64;
    for (f, t) in frequencies.iter().zip(targets) {
        require(*t > 0.0 && t.is_finite(), "invalid PDN target impedance")?;
        worst = worst.max(pdn_parallel_impedance_ohm(caps, *f, plane)?.norm() / t);
    }
    Ok(worst)
}

/// Greedy part selection minimising the worst `|Z| / Z_target` over the band.
pub fn optimize_decoupling_greedy(
    frequencies_hz: &[f64],
    target_impedance_ohm: &[f64],
    library: &[DecouplingCapacitor],
    maximum_parts: usize,
    plane_impedance_ohm: C,
) -> Result<DecouplingOptimizationResult> {
    require(
        !frequencies_hz.is_empty()
            && frequencies_hz.len() == target_impedance_ohm.len()
            && !library.is_empty(),
        "invalid decoupling optimization inputs",
    )?;
    let mut r = DecouplingOptimizationResult::default();
    r.worst_ratio_before = worst_ratio(
        frequencies_hz,
        target_impedance_ohm,
        &[],
        plane_impedance_ohm,
    )?;
    r.worst_ratio_after = r.worst_ratio_before;
    for _ in 0..maximum_parts {
        let mut best = r.worst_ratio_after;
        let mut best_selection = None;
        for cap in library {
            let mut trial = r.selected.clone();
            trial.push(*cap);
            let ratio = worst_ratio(
                frequencies_hz,
                target_impedance_ohm,
                &trial,
                plane_impedance_ohm,
            )?;
            if ratio < best {
                best = ratio;
                best_selection = Some(trial);
            }
        }
        let Some(selection) = best_selection else {
            break;
        };
        r.selected = selection;
        r.worst_ratio_after = best;
        r.iterations += 1;
        if r.worst_ratio_after <= 1.0 {
            break;
        }
    }
    Ok(r)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PdnGridConfig {
    pub nx: usize,
    pub ny: usize,
    pub horizontal_resistance_ohm: f64,
    pub vertical_resistance_ohm: f64,
    pub max_iterations: usize,
    pub relative_tolerance: f64,
}

impl Default for PdnGridConfig {
    fn default() -> Self {
        Self {
            nx: 2,
            ny: 2,
            horizontal_resistance_ohm: 1e-3,
            vertical_resistance_ohm: 1e-3,
            max_iterations: 20000,
            relative_tolerance: 1e-12,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PdnVoltageSource {
    pub node: usize,
    pub voltage_v: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PdnIrDropResult {
    pub voltage_v: Vec<f64>,
    pub minimum_voltage_v: f64,
    pub maximum_drop_v: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// Gauss-Seidel resistive-grid IR drop; positive loads draw current to ground.
pub fn solve_pdn_ir_drop_grid(
    config: &PdnGridConfig,
    load_current_a: &[f64],
    voltage_sources: &[PdnVoltageSource],
) -> Result<PdnIrDropResult> {
    let c = config;
    require(
        c.nx >= 2
            && c.ny >= 2
            && c.horizontal_resistance_ohm > 0.0
            && c.vertical_resistance_ohm > 0.0
            && c.max_iterations > 0
            && c.relative_tolerance > 0.0,
        "invalid PDN grid",
    )?;
    let n =
        c.nx.checked_mul(c.ny)
            .ok_or_else(|| SolverError::InvalidInput("invalid PDN grid".into()))?;
    require(
        load_current_a.len() == n && !voltage_sources.is_empty(),
        "PDN load/source size mismatch",
    )?;
    let mut fixed: Vec<Option<f64>> = vec![None; n];
    let mut voltage = vec![0.0; n];
    for s in voltage_sources {
        require(
            s.node < n && s.voltage_v.is_finite(),
            "invalid PDN voltage source",
        )?;
        fixed[s.node] = Some(s.voltage_v);
        voltage[s.node] = s.voltage_v;
    }
    let fixed_values: Vec<f64> = fixed.iter().flatten().copied().collect();
    let source_mean = fixed_values.iter().sum::<f64>() / fixed_values.len() as f64;
    for (v, f) in voltage.iter_mut().zip(&fixed) {
        if f.is_none() {
            *v = source_mean;
        }
    }
    let (gx, gy) = (
        1.0 / c.horizontal_resistance_ohm,
        1.0 / c.vertical_resistance_ohm,
    );
    let mut residual0 = 0.0;
    let mut result = PdnIrDropResult::default();
    for iteration in 0..c.max_iterations {
        let mut residual = 0.0_f64;
        for y in 0..c.ny {
            for x in 0..c.nx {
                let node = y * c.nx + x;
                if let Some(f) = fixed[node] {
                    voltage[node] = f;
                    continue;
                }
                let (mut g, mut rhs) = (0.0, -load_current_a[node]);
                if x > 0 {
                    g += gx;
                    rhs += gx * voltage[node - 1];
                }
                if x + 1 < c.nx {
                    g += gx;
                    rhs += gx * voltage[node + 1];
                }
                if y > 0 {
                    g += gy;
                    rhs += gy * voltage[node - c.nx];
                }
                if y + 1 < c.ny {
                    g += gy;
                    rhs += gy * voltage[node + c.nx];
                }
                let next = rhs / g;
                residual = residual.max((next - voltage[node]).abs());
                voltage[node] = next;
            }
        }
        if iteration == 0 {
            residual0 = residual.max(1e-30);
        }
        result.iterations = iteration + 1;
        if residual <= c.relative_tolerance * residual0 {
            result.converged = true;
            break;
        }
    }
    result.minimum_voltage_v = voltage.iter().copied().fold(f64::INFINITY, f64::min);
    let maximum_source = fixed_values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    result.maximum_drop_v = maximum_source - result.minimum_voltage_v;
    result.voltage_v = voltage;
    Ok(result)
}

/// Double-exponential pulse normalised to `amplitude` at its peak; zero before `t = 0`.
pub fn normalized_double_exponential_pulse(
    time_s: f64,
    amplitude: f64,
    rise_tau_s: f64,
    fall_tau_s: f64,
) -> Result<f64> {
    if time_s < 0.0 {
        return Ok(0.0);
    }
    require(
        amplitude >= 0.0 && rise_tau_s > 0.0 && fall_tau_s > rise_tau_s,
        "invalid double exponential pulse",
    )?;
    let peak_time =
        rise_tau_s * fall_tau_s * (fall_tau_s / rise_tau_s).ln() / (fall_tau_s - rise_tau_s);
    let peak = (-peak_time / fall_tau_s).exp() - (-peak_time / rise_tau_s).exp();
    Ok(amplitude * ((-time_s / fall_tau_s).exp() - (-time_s / rise_tau_s).exp()) / peak)
}

fn sample_count(duration: f64, dt: f64) -> Result<usize> {
    let n = (duration / dt).floor();
    require(n.is_finite() && n < 1e8, "waveform sample count too large")?;
    Ok(n as usize + 1)
}

pub fn sample_double_exponential_pulse(
    duration_s: f64,
    dt_s: f64,
    amplitude: f64,
    rise_tau_s: f64,
    fall_tau_s: f64,
) -> Result<Vec<f64>> {
    require(duration_s > 0.0 && dt_s > 0.0, "invalid waveform sampling")?;
    (0..sample_count(duration_s, dt_s)?)
        .map(|i| {
            normalized_double_exponential_pulse(i as f64 * dt_s, amplitude, rise_tau_s, fall_tau_s)
        })
        .collect()
}

pub fn sample_damped_sine(
    duration_s: f64,
    dt_s: f64,
    amplitude: f64,
    frequency_hz: f64,
    decay_tau_s: f64,
) -> Result<Vec<f64>> {
    require(
        duration_s > 0.0 && dt_s > 0.0 && frequency_hz >= 0.0 && decay_tau_s > 0.0,
        "invalid damped sine sampling",
    )?;
    Ok((0..sample_count(duration_s, dt_s)?)
        .map(|i| {
            let t = i as f64 * dt_s;
            amplitude * (-t / decay_tau_s).exp() * (2.0 * PI * frequency_hz * t).sin()
        })
        .collect())
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ProbeMetrics {
    pub peak: f64,
    pub rms: f64,
    pub impulse: f64,
    pub energy: f64,
}

pub fn probe_waveform(samples: &[f64], dt_s: f64) -> Result<ProbeMetrics> {
    require(!samples.is_empty() && dt_s > 0.0, "invalid probe waveform")?;
    require(
        samples.iter().all(|v| v.is_finite()),
        "non-finite waveform sample",
    )?;
    let mut out = ProbeMetrics::default();
    let mut sum2 = 0.0;
    for v in samples {
        out.peak = out.peak.max(v.abs());
        sum2 += v * v;
        out.impulse += v * dt_s;
        out.energy += v * v * dt_s;
    }
    out.rms = (sum2 / samples.len() as f64).sqrt();
    Ok(out)
}
