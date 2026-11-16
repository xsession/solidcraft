//! Resonator wake functions, bunch wake potentials and wake impedance. Ported from upstream
//! `cfd/particle/wakefield` (MIT, cfd_solvers contributors).
use crate::{require, Result};
use num_complex::Complex64;
use std::f64::consts::PI;

const C0: f64 = 299_792_458.0;

/// Causal single-resonator longitudinal point-charge wake W(s) [V/C] at distance `s ≥ 0`.
pub fn resonator_longitudinal_wake_v_per_c(
    distance: f64,
    rs: f64,
    frequency: f64,
    q: f64,
) -> Result<f64> {
    require(
        distance >= 0.0 && rs >= 0.0 && frequency > 0.0 && q > 0.5,
        "invalid resonator wake controls",
    )?;
    let omega = 2.0 * PI * frequency;
    let alpha = omega / (2.0 * q);
    let omega_d = omega * (1.0 - 1.0 / (4.0 * q * q)).sqrt();
    let t = distance / C0;
    Ok((omega * rs / q)
        * (-alpha * t).exp()
        * ((omega_d * t).cos() - (alpha / omega_d) * (omega_d * t).sin()))
}

/// `V(s_i) = −Σ_lag W(lag) λ(s_i − lag) Δs` for a causal sampled wake.
pub fn bunch_wake_potential_v(wake: &[f64], line_charge: &[f64], ds: f64) -> Result<Vec<f64>> {
    require(
        !wake.is_empty() && !line_charge.is_empty() && ds > 0.0,
        "invalid bunch-wake convolution inputs",
    )?;
    Ok((0..line_charge.len())
        .map(|i| {
            let mut v = 0.0;
            for lag in 0..(i + 1).min(wake.len()) {
                v -= wake[lag] * line_charge[i - lag] * ds;
            }
            v
        })
        .collect())
}

/// Trapezoidal `Z(ω) = ∫ W(t) e^{−jωt} dt` with `t = s/c`.
pub fn wake_impedance_ohm(wake: &[f64], ds: f64, frequency: f64) -> Result<Complex64> {
    require(
        !wake.is_empty() && ds > 0.0 && frequency >= 0.0,
        "invalid wake-impedance controls",
    )?;
    let omega = 2.0 * PI * frequency;
    let dt = ds / C0;
    let mut z = Complex64::default();
    for (i, &w) in wake.iter().enumerate() {
        let t = i as f64 * dt;
        let weight = if i == 0 || i + 1 == wake.len() {
            0.5
        } else {
            1.0
        };
        z += weight * w * Complex64::new((-omega * t).cos(), (-omega * t).sin()) * dt;
    }
    Ok(z)
}
