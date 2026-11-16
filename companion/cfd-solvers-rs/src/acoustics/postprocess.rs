//! Plane-wave intensity, radiation pressure, absorption heating and delay-and-sum beamforming,
//! ported from upstream `cfd/acoustics/postprocess` (MIT, cfd_solvers contributors).
use crate::{require, Result};

/// `I = p_rms^2 / (rho c)`.
pub fn plane_wave_intensity_w_m2(
    pressure_rms_pa: f64,
    density_kg_m3: f64,
    sound_speed_m_s: f64,
) -> Result<f64> {
    let (p, rho, c) = (pressure_rms_pa, density_kg_m3, sound_speed_m_s);
    require(
        p >= 0.0 && rho > 0.0 && c > 0.0 && p.is_finite() && rho.is_finite() && c.is_finite(),
        "invalid acoustic intensity inputs",
    )?;
    Ok(p * p / (rho * c))
}

/// `(1 + R) I / c` for reflection coefficient `R` in `[0, 1]`.
pub fn radiation_pressure_pa(
    intensity_w_m2: f64,
    sound_speed_m_s: f64,
    reflection_coefficient: f64,
) -> Result<f64> {
    let (i, c, r) = (intensity_w_m2, sound_speed_m_s, reflection_coefficient);
    require(
        i >= 0.0 && c > 0.0 && (0.0..=1.0).contains(&r) && i.is_finite() && c.is_finite(),
        "invalid radiation pressure inputs",
    )?;
    Ok((1.0 + r) * i / c)
}

/// `q = 2 alpha I`.
pub fn acoustic_heating_w_m3(intensity_w_m2: f64, absorption_np_per_m: f64) -> Result<f64> {
    let (i, a) = (intensity_w_m2, absorption_np_per_m);
    require(
        i >= 0.0 && a >= 0.0 && i.is_finite() && a.is_finite(),
        "invalid acoustic heating inputs",
    )?;
    Ok(2.0 * a * i)
}

#[derive(Clone, Copy, Debug)]
pub struct BeamformSample<'a> {
    pub x_m: f64,
    pub y_m: f64,
    pub pressure: &'a [f64],
}

/// Linear-interpolated delay-and-sum toward a focus, delays relative to the nearest sensor.
pub fn delay_and_sum_beamform(
    sensors: &[BeamformSample<'_>],
    focus_x_m: f64,
    focus_y_m: f64,
    sample_dt: f64,
    sound_speed_m_s: f64,
    output_samples: usize,
) -> Result<Vec<f64>> {
    require(
        !sensors.is_empty() && sample_dt > 0.0 && sound_speed_m_s > 0.0 && output_samples > 0,
        "invalid beamforming controls",
    )?;
    let distance: Vec<f64> = sensors
        .iter()
        .map(|s| (s.x_m - focus_x_m).hypot(s.y_m - focus_y_m))
        .collect();
    let dmin = distance.iter().copied().fold(f64::INFINITY, f64::min);
    let mut out = vec![0.0; output_samples];
    for (s, d) in sensors.iter().zip(&distance) {
        if s.pressure.is_empty() {
            continue;
        }
        let delay = (d - dmin) / sound_speed_m_s / sample_dt;
        for (n, o) in out.iter_mut().enumerate() {
            let src = n as f64 + delay;
            let i0 = src.floor() as usize;
            if i0 + 1 >= s.pressure.len() {
                continue;
            }
            let frac = src - i0 as f64;
            *o += (1.0 - frac) * s.pressure[i0] + frac * s.pressure[i0 + 1];
        }
    }
    let scale = 1.0 / sensors.len() as f64;
    for v in &mut out {
        *v *= scale;
    }
    Ok(out)
}
