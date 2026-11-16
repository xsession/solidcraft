//! Analytic antennas: Hertzian and sinusoidal dipoles, array factors, polarization, matching and a
//! bounded pattern-search optimizer.
use crate::{require, Complex64 as C, Result, SolverError};
use std::f64::consts::PI;

const C0: f64 = 299792458.0;
const ETA0: f64 = 376.730313668;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HertzianFarField {
    pub e_theta: C,
    pub h_phi: C,
}

pub fn hertzian_dipole_far_field(
    length_m: f64,
    current_a: f64,
    frequency_hz: f64,
    theta_rad: f64,
    distance_m: f64,
) -> Result<HertzianFarField> {
    require(
        length_m > 0.0
            && length_m.is_finite()
            && current_a.is_finite()
            && frequency_hz > 0.0
            && frequency_hz.is_finite()
            && theta_rad.is_finite()
            && distance_m > 0.0
            && distance_m.is_finite(),
        "invalid Hertzian-dipole controls",
    )?;
    let k = 2.0 * PI * frequency_hz / C0;
    let phase = C::new(0.0, -k * distance_m).exp();
    let e = C::new(0.0, 1.0) * ETA0 * k * current_a * length_m * theta_rad.sin() * phase
        / (4.0 * PI * distance_m);
    Ok(HertzianFarField {
        e_theta: e,
        h_phi: e / ETA0,
    })
}

pub fn shielding_effectiveness_db(
    incident_magnitude: f64,
    transmitted_magnitude: f64,
) -> Result<f64> {
    require(
        incident_magnitude > 0.0
            && transmitted_magnitude > 0.0
            && incident_magnitude.is_finite()
            && transmitted_magnitude.is_finite(),
        "shielding field magnitudes must be finite and positive",
    )?;
    Ok(20.0 * (incident_magnitude / transmitted_magnitude).log10())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArrayElement {
    pub x_m: f64,
    pub y_m: f64,
    pub z_m: f64,
    pub weight: C,
}

impl Default for ArrayElement {
    fn default() -> Self {
        Self {
            x_m: 0.0,
            y_m: 0.0,
            z_m: 0.0,
            weight: C::new(1.0, 0.0),
        }
    }
}

/// Normalized by the sum of weight magnitudes.
pub fn array_factor(
    elements: &[ArrayElement],
    frequency_hz: f64,
    theta_rad: f64,
    phi_rad: f64,
) -> Result<C> {
    require(
        !elements.is_empty()
            && frequency_hz > 0.0
            && frequency_hz.is_finite()
            && theta_rad.is_finite()
            && phi_rad.is_finite(),
        "invalid phased-array controls",
    )?;
    let k = 2.0 * PI * frequency_hz / C0;
    let (sx, sy, sz) = (
        theta_rad.sin() * phi_rad.cos(),
        theta_rad.sin() * phi_rad.sin(),
        theta_rad.cos(),
    );
    let (mut result, mut norm) = (C::new(0.0, 0.0), 0.0);
    for e in elements {
        require(
            e.x_m.is_finite() && e.y_m.is_finite() && e.z_m.is_finite(),
            "non-finite array-element coordinate",
        )?;
        let phase = k * (e.x_m * sx + e.y_m * sy + e.z_m * sz);
        result += e.weight * C::new(0.0, phase).exp();
        norm += e.weight.norm();
    }
    Ok(if norm > 0.0 {
        result / norm
    } else {
        C::new(0.0, 0.0)
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PolarizationMetrics {
    pub axial_ratio: f64,
    pub tilt_angle_rad: f64,
    pub ellipticity_angle_rad: f64,
    pub right_hand: bool,
}

/// Stokes-parameter polarization ellipse; linear polarization has an infinite axial ratio.
pub fn polarization_metrics(e_theta: C, e_phi: C) -> Result<PolarizationMetrics> {
    let s0 = e_theta.norm_sqr() + e_phi.norm_sqr();
    require(
        s0 > 0.0 && s0.is_finite(),
        "polarization field must be nonzero and finite",
    )?;
    let s1 = e_theta.norm_sqr() - e_phi.norm_sqr();
    let x = e_theta * e_phi.conj();
    let (s2, s3) = (2.0 * x.re, -2.0 * x.im);
    let linear = s1.hypot(s2);
    let major2 = 0.5 * (s0 + linear);
    let minor2 = (0.5 * (s0 - linear)).max(0.0);
    Ok(PolarizationMetrics {
        axial_ratio: if minor2 > 1e-30 {
            (major2 / minor2).sqrt()
        } else {
            f64::INFINITY
        },
        tilt_angle_rad: 0.5 * s2.atan2(s1),
        ellipticity_angle_rad: 0.5 * (s3 / s0).clamp(-1.0, 1.0).asin(),
        right_hand: s3 < 0.0,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PatternSample {
    pub theta_rad: f64,
    pub normalized_power: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AntennaMatchMetrics {
    pub impedance_ohm: C,
    pub reflection_coefficient: C,
    pub return_loss_db: f64,
    pub vswr: f64,
}

/// Reference impedance defaults to 50 ohm upstream.
pub fn antenna_match_metrics(
    feed_impedance_ohm: C,
    reference_impedance_ohm: f64,
) -> Result<AntennaMatchMetrics> {
    require(
        reference_impedance_ohm > 0.0
            && reference_impedance_ohm.is_finite()
            && feed_impedance_ohm.re.is_finite()
            && feed_impedance_ohm.im.is_finite(),
        "invalid antenna match impedance",
    )?;
    let reflection = (feed_impedance_ohm - reference_impedance_ohm)
        / (feed_impedance_ohm + reference_impedance_ohm);
    let m = reflection.norm();
    Ok(AntennaMatchMetrics {
        impedance_ohm: feed_impedance_ohm,
        reflection_coefficient: reflection,
        return_loss_db: if m == 0.0 {
            f64::INFINITY
        } else {
            -20.0 * m.log10()
        },
        vswr: if m >= 1.0 || m.is_nan() {
            f64::INFINITY
        } else {
            (1.0 + m) / (1.0 - m)
        },
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AntennaOptimizationVariable {
    pub initial_value: f64,
    pub minimum: f64,
    pub maximum: f64,
    pub initial_step: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AntennaOptimizationConfig {
    pub max_iterations: usize,
    pub parameter_tolerance: f64,
    pub step_reduction: f64,
}

impl Default for AntennaOptimizationConfig {
    fn default() -> Self {
        Self {
            max_iterations: 128,
            parameter_tolerance: 1e-6,
            step_reduction: 0.5,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AntennaOptimizationResult {
    pub parameters: Vec<f64>,
    pub objective: f64,
    pub iterations: usize,
    pub evaluations: usize,
    pub converged: bool,
}

/// Bounded derivative-free coordinate/pattern search; lower objective is better.
pub fn optimize_antenna(
    variables: &[AntennaOptimizationVariable],
    mut objective: impl FnMut(&[f64]) -> f64,
    config: &AntennaOptimizationConfig,
) -> Result<AntennaOptimizationResult> {
    require(
        !variables.is_empty()
            && config.max_iterations > 0
            && config.parameter_tolerance > 0.0
            && config.parameter_tolerance.is_finite()
            && config.step_reduction > 0.0
            && config.step_reduction < 1.0,
        "invalid antenna optimization controls",
    )?;
    let mut result = AntennaOptimizationResult {
        parameters: vec![0.0; variables.len()],
        ..Default::default()
    };
    let mut step = vec![0.0; variables.len()];
    for (i, v) in variables.iter().enumerate() {
        require(
            v.initial_value.is_finite()
                && v.minimum.is_finite()
                && v.maximum.is_finite()
                && v.initial_step.is_finite()
                && v.maximum >= v.minimum
                && v.initial_step > 0.0,
            "invalid antenna optimization variable",
        )?;
        result.parameters[i] = v.initial_value.clamp(v.minimum, v.maximum);
        step[i] = v
            .initial_step
            .min((v.maximum - v.minimum).max(config.parameter_tolerance));
    }
    let mut evaluate = |values: &[f64], evaluations: &mut usize| -> Result<f64> {
        let score = objective(values);
        *evaluations += 1;
        if !score.is_finite() {
            return Err(SolverError::NumericalFailure(
                "antenna objective returned non-finite value".into(),
            ));
        }
        Ok(score)
    };
    result.objective = evaluate(&result.parameters, &mut result.evaluations)?;
    for iteration in 0..config.max_iterations {
        result.iterations = iteration + 1;
        let mut improved = false;
        for (variable, v) in variables.iter().enumerate() {
            let original = result.parameters[variable];
            let (mut best_value, mut best_score) = (original, result.objective);
            for direction in [-1.0, 1.0] {
                let candidate = (original + direction * step[variable]).clamp(v.minimum, v.maximum);
                if candidate == original {
                    continue;
                }
                result.parameters[variable] = candidate;
                let score = evaluate(&result.parameters, &mut result.evaluations)?;
                if score < best_score {
                    best_score = score;
                    best_value = candidate;
                }
            }
            result.parameters[variable] = best_value;
            if best_score < result.objective {
                result.objective = best_score;
                improved = true;
            }
        }
        if !improved {
            step.iter_mut().for_each(|s| *s *= config.step_reduction);
        }
        let largest = step.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if largest <= config.parameter_tolerance {
            result.converged = true;
            break;
        }
    }
    Ok(result)
}

/// Center-fed dipole with the ideal sinusoidal current distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SinusoidalDipoleSolver {
    length: f64,
    frequency: f64,
    current: f64,
}

impl SinusoidalDipoleSolver {
    /// Feed current defaults to 1 A upstream.
    pub fn new(length_m: f64, frequency_hz: f64, feed_current_a: f64) -> Result<Self> {
        require(
            length_m > 0.0
                && frequency_hz > 0.0
                && feed_current_a > 0.0
                && length_m.is_finite()
                && frequency_hz.is_finite()
                && feed_current_a.is_finite(),
            "invalid dipole parameters",
        )?;
        Ok(Self {
            length: length_m,
            frequency: frequency_hz,
            current: feed_current_a,
        })
    }
    pub fn wavelength(&self) -> f64 {
        C0 / self.frequency
    }
    pub fn electrical_length(&self) -> f64 {
        self.length / self.wavelength()
    }
    pub fn pattern_factor(&self, theta_rad: f64) -> f64 {
        let k = 2.0 * PI / self.wavelength();
        let s = theta_rad.sin();
        if s.abs() < 1e-10 {
            return 0.0;
        }
        ((0.5 * k * self.length * theta_rad.cos()).cos() - (0.5 * k * self.length).cos()) / s
    }
    pub fn far_e_theta(&self, theta_rad: f64, distance_m: f64) -> Result<C> {
        require(
            distance_m > 0.0 && distance_m.is_finite(),
            "far-field distance must be positive",
        )?;
        let k = 2.0 * PI / self.wavelength();
        let phase = C::new(0.0, -k * distance_m).exp();
        Ok(
            C::new(0.0, 1.0) * ETA0 * self.current / (2.0 * PI * distance_m)
                * self.pattern_factor(theta_rad)
                * phase,
        )
    }
    /// Simpson integration over theta (upstream default 8192 points).
    pub fn radiated_power(&self, integration_points: usize) -> Result<f64> {
        require(
            integration_points >= 16,
            "dipole integration requires at least 16 points",
        )?;
        let n = integration_points + integration_points % 2;
        let h = PI / n as f64;
        let mut integral = 0.0;
        for i in 0..=n {
            let theta = h * i as f64;
            let f = self.pattern_factor(theta);
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            integral += weight * f * f * theta.sin();
        }
        integral *= h / 3.0;
        Ok(ETA0 * self.current * self.current / (4.0 * PI) * integral)
    }
    pub fn radiation_resistance(&self, integration_points: usize) -> Result<f64> {
        Ok(2.0 * self.radiated_power(integration_points)? / (self.current * self.current))
    }
    pub fn directivity(&self, integration_points: usize) -> Result<f64> {
        let p = self.radiated_power(integration_points)?;
        let n = integration_points;
        let mut max_u: f64 = 0.0;
        for i in 0..=n {
            let f = self.pattern_factor(PI * i as f64 / n as f64);
            max_u = max_u.max(ETA0 * self.current * self.current / (8.0 * PI * PI) * f * f);
        }
        Ok(4.0 * PI * max_u / p)
    }
    /// Upstream default 181 samples.
    pub fn normalized_pattern(&self, samples: usize) -> Result<Vec<PatternSample>> {
        require(samples >= 2, "dipole pattern requires at least two samples")?;
        let mut out: Vec<PatternSample> = (0..samples)
            .map(|i| {
                let theta = PI * i as f64 / (samples - 1) as f64;
                let f = self.pattern_factor(theta);
                PatternSample {
                    theta_rad: theta,
                    normalized_power: f * f,
                }
            })
            .collect();
        let maximum = out.iter().map(|s| s.normalized_power).fold(0.0, f64::max);
        if maximum > 0.0 {
            out.iter_mut().for_each(|s| s.normalized_power /= maximum);
        }
        Ok(out)
    }
}
