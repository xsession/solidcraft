//! Circuit-level analyses on top of the MNA kernel: two-port and N-port S-parameter extraction,
//! DC sensitivity, quasi-static distortion, Fourier/THD measurements, resistor Monte Carlo,
//! parameter and temperature sweeps, dominant-pole and rational pole/zero fits, shooting-by-
//! settling periodic steady state, integrated noise spectra and the SPICE `pnjlim` limiter.
//! Ported from upstream `cfd/circuit/analysis` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::spice::{
    solve_dense, Circuit, Complex, NewtonConfig, NoiseResult, OperatingPoint, TransientMethod,
    TransientPoint,
};
use crate::core::{Mt19937_64, StdNormal};
use crate::rf::{self, ComplexMatrix, Matrix2C};
use crate::{require, Result, SolverError};
use std::f64::consts::PI;

fn err_runtime(msg: &str) -> SolverError {
    SolverError::NumericalFailure(msg.into())
}

/// A port given by node names (negative defaults to ground).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AcPort {
    pub positive: String,
    pub negative: String,
}

impl AcPort {
    pub fn new(positive: impl Into<String>) -> Self {
        Self {
            positive: positive.into(),
            negative: "0".into(),
        }
    }
}

/// Two-port Z by 1-A AC current injection at each port, converted to S. Independent AC sources
/// already in the circuit should be zero for a passive characterization.
pub fn two_port_s_parameters(
    circuit: &Circuit,
    first: &AcPort,
    second: &AcPort,
    frequency: f64,
    reference_impedance: f64,
) -> Result<Matrix2C> {
    require(
        frequency > 0.0 && frequency.is_finite() && reference_impedance > 0.0,
        "invalid two-port S-parameter controls",
    )?;
    let excite = |driven: &AcPort| -> Result<[Complex; 2]> {
        let mut copy = circuit.clone();
        let (dp, dn) = (
            copy.find_node(&driven.positive)?,
            copy.find_node(&driven.negative)?,
        );
        copy.add_current_source(
            "__PORT_EXCITATION__",
            dp,
            dn,
            0.0,
            Complex::new(-1.0, 0.0),
            None,
        )?;
        let ac = copy.ac(frequency, None, &NewtonConfig::default())?;
        let voltage = |port: &AcPort| -> Result<Complex> {
            let (p, n) = (
                copy.find_node(&port.positive)?,
                copy.find_node(&port.negative)?,
            );
            Ok(ac.node_voltage[p] - ac.node_voltage[n])
        };
        Ok([voltage(first)?, voltage(second)?])
    };
    let c1 = excite(first)?;
    let c2 = excite(second)?;
    rf::z_to_s(
        Matrix2C {
            a11: c1[0],
            a12: c2[0],
            a21: c1[1],
            a22: c2[1],
        },
        reference_impedance,
    )
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SensitivityResult {
    pub baseline_output: f64,
    /// Output volts per unit parameter.
    pub derivative: f64,
}

/// Central-difference `dV(output)/dV(source)` with homotopy operating points.
pub fn voltage_source_dc_sensitivity(
    mut circuit: Circuit,
    source: &str,
    output: &str,
    perturbation: f64,
    cfg: &NewtonConfig,
) -> Result<SensitivityResult> {
    require(
        perturbation > 0.0 && perturbation.is_finite(),
        "sensitivity perturbation must be positive",
    )?;
    let nominal = circuit.voltage_source_dc(source)?;
    let base = circuit.dc_operating_point_default_homotopy(cfg)?;
    if !base.converged {
        return Err(err_runtime("sensitivity baseline operating point failed"));
    }
    circuit.set_voltage_source_dc(source, nominal + perturbation)?;
    let plus = circuit.dc_operating_point_default_homotopy(cfg)?;
    circuit.set_voltage_source_dc(source, nominal - perturbation)?;
    let minus = circuit.dc_operating_point_default_homotopy(cfg)?;
    if !plus.converged || !minus.converged {
        return Err(err_runtime("sensitivity perturbed operating point failed"));
    }
    Ok(SensitivityResult {
        baseline_output: circuit.voltage(&base, output)?,
        derivative: (circuit.voltage(&plus, output)? - circuit.voltage(&minus, output)?)
            / (2.0 * perturbation),
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SmallSignalDistortion {
    pub bias_output_v: f64,
    pub first_derivative: f64,
    pub second_derivative: f64,
    pub third_derivative: f64,
    pub fundamental_v: f64,
    pub second_harmonic_v: f64,
    pub third_harmonic_v: f64,
    pub hd2: f64,
    pub hd3: f64,
}

/// Five-point DC transfer stencil → harmonics of a sinusoidal source perturbation.
pub fn small_signal_distortion(
    mut circuit: Circuit,
    source: &str,
    output: &str,
    amplitude: f64,
    h: f64,
    cfg: &NewtonConfig,
) -> Result<SmallSignalDistortion> {
    require(
        amplitude > 0.0 && h > 0.0 && amplitude.is_finite() && h.is_finite(),
        "invalid small-signal distortion controls",
    )?;
    let bias = circuit.voltage_source_dc(source)?;
    let mut evaluate = |offset: f64| -> Result<f64> {
        circuit.set_voltage_source_dc(source, bias + offset)?;
        let op = circuit.dc_operating_point_default_homotopy(cfg)?;
        if !op.converged {
            return Err(err_runtime("distortion operating point failed"));
        }
        circuit.voltage(&op, output)
    };
    let fm2 = evaluate(-2.0 * h)?;
    let fm1 = evaluate(-h)?;
    let f0 = evaluate(0.0)?;
    let fp1 = evaluate(h)?;
    let fp2 = evaluate(2.0 * h)?;
    let d1 = (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h);
    let d2 = (-fp2 + 16.0 * fp1 - 30.0 * f0 + 16.0 * fm1 - fm2) / (12.0 * h * h);
    let d3 = (fp2 - 2.0 * fp1 + 2.0 * fm1 - fm2) / (2.0 * h * h * h);
    let fundamental = (d1 * amplitude + d3 * amplitude * amplitude * amplitude / 8.0).abs();
    let second = (d2 * amplitude * amplitude / 4.0).abs();
    let third = (d3 * amplitude * amplitude * amplitude / 24.0).abs();
    let ratio = |x: f64| {
        if fundamental > 0.0 {
            x / fundamental
        } else {
            f64::INFINITY
        }
    };
    Ok(SmallSignalDistortion {
        bias_output_v: f0,
        first_derivative: d1,
        second_derivative: d2,
        third_derivative: d3,
        fundamental_v: fundamental,
        second_harmonic_v: second,
        third_harmonic_v: third,
        hd2: ratio(second),
        hd3: ratio(third),
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FourierMeasurement {
    pub frequency_hz: f64,
    /// Peak-amplitude convention for a real signal.
    pub phasor: Complex,
    pub magnitude: f64,
    pub phase_rad: f64,
}

/// Trapezoidal Fourier coefficient of one node voltage over the sampled window.
pub fn fourier_measurement(
    samples: &[TransientPoint],
    node: usize,
    frequency: f64,
) -> Result<FourierMeasurement> {
    require(
        samples.len() >= 2 && frequency >= 0.0 && frequency.is_finite(),
        "invalid Fourier measurement",
    )?;
    let duration = samples[samples.len() - 1].time - samples[0].time;
    require(duration > 0.0, "transient samples need increasing time")?;
    let mut integral = Complex::default();
    for w in samples.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let (Some(&va), Some(&vb)) = (a.node_voltage.get(node), b.node_voltage.get(node)) else {
            return Err(SolverError::InvalidInput(
                "Fourier node index out of range".into(),
            ));
        };
        let dt = b.time - a.time;
        require(dt > 0.0, "transient sample times must increase")?;
        let e0 = Complex::new(0.0, -2.0 * PI * frequency * a.time).exp();
        let e1 = Complex::new(0.0, -2.0 * PI * frequency * b.time).exp();
        integral += 0.5 * dt * (va * e0 + vb * e1);
    }
    let phasor = (if frequency == 0.0 { 1.0 } else { 2.0 }) * integral / duration;
    Ok(FourierMeasurement {
        frequency_hz: frequency,
        phasor,
        magnitude: phasor.norm(),
        phase_rad: phasor.arg(),
    })
}

pub fn total_harmonic_distortion(
    samples: &[TransientPoint],
    node: usize,
    fundamental: f64,
    harmonics: usize,
) -> Result<f64> {
    require(fundamental > 0.0 && harmonics >= 2, "invalid THD controls")?;
    let base = fourier_measurement(samples, node, fundamental)?.magnitude;
    if base.is_nan() || base <= 0.0 {
        return Ok(f64::INFINITY);
    }
    let mut sum = 0.0;
    for h in 2..=harmonics {
        let a = fourier_measurement(samples, node, fundamental * h as f64)?.magnitude;
        sum += a * a;
    }
    Ok(sum.sqrt() / base)
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResistorTolerance {
    pub element: String,
    pub relative_sigma: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MonteCarloResult {
    pub samples: usize,
    pub mean: f64,
    pub standard_deviation: f64,
    pub minimum: f64,
    pub maximum: f64,
}

/// Gaussian resistor variations drawn with `std::mt19937_64` + `std::normal_distribution`
/// (bit-identical sequence), Welford statistics of the output node voltage.
pub fn monte_carlo_dc_resistors(
    nominal: &Circuit,
    output: &str,
    tolerances: &[ResistorTolerance],
    samples: usize,
    seed: u64,
    cfg: &NewtonConfig,
) -> Result<MonteCarloResult> {
    require(samples > 0, "Monte Carlo requires at least one sample")?;
    require(
        tolerances
            .iter()
            .all(|t| t.relative_sigma >= 0.0 && t.relative_sigma.is_finite()),
        "invalid resistor tolerance",
    )?;
    let mut rng = Mt19937_64::new(seed);
    let mut normal = StdNormal::default();
    let mut result = MonteCarloResult {
        samples,
        minimum: f64::INFINITY,
        maximum: f64::NEG_INFINITY,
        ..Default::default()
    };
    let (mut mean, mut m2) = (0.0, 0.0);
    for sample in 0..samples {
        let mut circuit = nominal.clone();
        for t in tolerances {
            let r = nominal.resistance(&t.element)?;
            let mut varied = r * (1.0 + t.relative_sigma * normal.sample(&mut rng));
            if varied.is_nan() || varied <= 0.0 {
                varied = (r * 1.0e-6).max(1.0e-18);
            }
            circuit.set_resistance(&t.element, varied)?;
        }
        let op = circuit.dc_operating_point_default_homotopy(cfg)?;
        if !op.converged {
            return Err(err_runtime("Monte Carlo operating point failed"));
        }
        let value = circuit.voltage(&op, output)?;
        let delta = value - mean;
        mean += delta / (sample + 1) as f64;
        m2 += delta * (value - mean);
        result.minimum = result.minimum.min(value);
        result.maximum = result.maximum.max(value);
    }
    result.mean = mean;
    result.standard_deviation = (m2 / samples as f64).max(0.0).sqrt();
    Ok(result)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircuitPort {
    pub positive: usize,
    pub negative: usize,
    pub reference_impedance: f64,
}

impl CircuitPort {
    pub fn new(positive: usize, negative: usize) -> Self {
        Self {
            positive,
            negative,
            reference_impedance: 50.0,
        }
    }
}

/// Voltage-force each port in turn; source currents give Y, then Y → Z → S.
pub fn extract_s_parameters(
    circuit: &Circuit,
    ports: &[CircuitPort],
    frequency: f64,
    cfg: &NewtonConfig,
) -> Result<ComplexMatrix> {
    require(
        !ports.is_empty() && frequency > 0.0 && frequency.is_finite(),
        "invalid circuit S-parameter extraction controls",
    )?;
    let mut reference = Vec::with_capacity(ports.len());
    for p in ports {
        require(
            p.positive < circuit.node_count()
                && p.negative < circuit.node_count()
                && p.reference_impedance > 0.0
                && p.reference_impedance.is_finite(),
            "invalid circuit port",
        )?;
        reference.push(p.reference_impedance);
    }
    let mut y = ComplexMatrix::new(ports.len());
    for excitation in 0..ports.len() {
        let mut driven = circuit.clone();
        let names: Vec<String> = (0..ports.len())
            .map(|p| format!("__CFD_PORT_{p}"))
            .collect();
        for (p, port) in ports.iter().enumerate() {
            let ac = if p == excitation {
                Complex::new(1.0, 0.0)
            } else {
                Complex::default()
            };
            driven.add_voltage_source(
                names[p].clone(),
                port.positive,
                port.negative,
                0.0,
                ac,
                None,
            )?;
        }
        let op = driven.dc_operating_point_default_homotopy(cfg)?;
        if !op.converged {
            return Err(err_runtime("circuit S-parameter operating point failed"));
        }
        let solution = driven.ac(frequency, Some(&op), cfg)?;
        for (response, name) in names.iter().enumerate() {
            // MNA source current flows + → − through the source; current into the DUT is its negative.
            y[(response, excitation)] = -driven.branch_current(&solution, name)?;
        }
    }
    let z = rf::inverse(&y)?;
    rf::nport_z_to_s(&z, &reference)
}

/// SPICE `pnjlim`: limit a PN-junction Newton step above the critical voltage.
pub fn limit_pn_junction_voltage(
    proposed: f64,
    previous: f64,
    vt: f64,
    saturation_current: f64,
) -> Result<f64> {
    require(
        vt > 0.0 && saturation_current > 0.0 && proposed.is_finite() && previous.is_finite(),
        "invalid PN junction limiter input",
    )?;
    let vcrit = vt * (vt / (2.0f64.sqrt() * saturation_current)).ln();
    if proposed > vcrit && (proposed - previous).abs() > 2.0 * vt {
        if previous > 0.0 {
            let argument = 1.0 + (proposed - previous) / vt;
            return Ok(if argument > 0.0 {
                previous + vt * argument.ln()
            } else {
                vcrit
            });
        }
        return Ok(vt * (proposed / vt).max(1.0).ln());
    }
    Ok(proposed)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParameterSweepPoint {
    pub parameter: f64,
    pub operating_point: OperatingPoint,
}

/// Apply `setter` to a copy of `nominal` for each value and solve with homotopy.
pub fn parameter_sweep_dc(
    nominal: &Circuit,
    values: &[f64],
    setter: &dyn Fn(&mut Circuit, f64) -> Result<()>,
    cfg: &NewtonConfig,
) -> Result<Vec<ParameterSweepPoint>> {
    require(!values.is_empty(), "invalid parameter sweep controls")?;
    let mut out = Vec::with_capacity(values.len());
    for &value in values {
        require(value.is_finite(), "non-finite sweep value")?;
        let mut circuit = nominal.clone();
        setter(&mut circuit, value)?;
        let op = circuit.dc_operating_point_default_homotopy(cfg)?;
        if !op.converged {
            return Err(err_runtime("parameter sweep operating point failed"));
        }
        out.push(ParameterSweepPoint {
            parameter: value,
            operating_point: op,
        });
    }
    Ok(out)
}

pub fn temperature_sweep_dc(
    nominal: &Circuit,
    temperatures: &[f64],
    cfg: &NewtonConfig,
) -> Result<Vec<ParameterSweepPoint>> {
    parameter_sweep_dc(
        nominal,
        temperatures,
        &|c, t| c.set_device_temperature(t),
        cfg,
    )
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct NoiseSpectrum {
    pub points: Vec<NoiseResult>,
    pub integrated_rms_voltage: f64,
}

/// Log-spaced output-noise sweep; RMS from trapezoidal integration of the PSD in Hz.
pub fn output_noise_spectrum(
    circuit: &Circuit,
    output: &str,
    start: f64,
    stop: f64,
    count: usize,
    temperature: f64,
    cfg: &NewtonConfig,
) -> Result<NoiseSpectrum> {
    require(
        start > 0.0
            && stop >= start
            && count > 0
            && temperature > 0.0
            && start.is_finite()
            && stop.is_finite()
            && temperature.is_finite(),
        "invalid noise spectrum controls",
    )?;
    let mut result = NoiseSpectrum::default();
    if count == 1 {
        result
            .points
            .push(circuit.output_noise(start, output, temperature, cfg)?);
        return Ok(result);
    }
    let ratio = (stop / start).powf(1.0 / (count - 1) as f64);
    for i in 0..count {
        let f = if i + 1 == count {
            stop
        } else {
            start * ratio.powf(i as f64)
        };
        result
            .points
            .push(circuit.output_noise(f, output, temperature, cfg)?);
    }
    let mut variance = 0.0;
    for w in result.points.windows(2) {
        variance += 0.5
            * (w[0].output_noise_density_v2_per_hz + w[1].output_noise_density_v2_per_hz)
            * (w[1].frequency_hz - w[0].frequency_hz);
    }
    result.integrated_rms_voltage = variance.max(0.0).sqrt();
    Ok(result)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DominantPoleEstimate {
    pub dc_gain: f64,
    pub pole_frequency_hz: f64,
    pub found: bool,
}

/// First −3 dB crossing of the AC magnitude (log-frequency interpolation).
pub fn estimate_dominant_pole(
    circuit: &Circuit,
    output: &str,
    start: f64,
    stop: f64,
    count: usize,
    cfg: &NewtonConfig,
) -> Result<DominantPoleEstimate> {
    require(
        start > 0.0 && stop > start && count >= 3,
        "invalid pole-estimation controls",
    )?;
    let sweep = circuit.ac_log_sweep(start, stop, count, cfg)?;
    let low = circuit.ac_voltage(&sweep[0], output)?.norm();
    if !(low > 0.0 && low.is_finite()) {
        return Ok(DominantPoleEstimate {
            dc_gain: low,
            ..Default::default()
        });
    }
    let target = low / 2.0f64.sqrt();
    for i in 1..sweep.len() {
        let previous = circuit.ac_voltage(&sweep[i - 1], output)?.norm();
        let current = circuit.ac_voltage(&sweep[i], output)?.norm();
        if previous >= target && current <= target && previous != current {
            let (x0, x1) = (sweep[i - 1].frequency_hz.ln(), sweep[i].frequency_hz.ln());
            let f = ((target - previous) / (current - previous)).clamp(0.0, 1.0);
            return Ok(DominantPoleEstimate {
                dc_gain: low,
                pole_frequency_hz: (x0 + f * (x1 - x0)).exp(),
                found: true,
            });
        }
    }
    Ok(DominantPoleEstimate {
        dc_gain: low,
        ..Default::default()
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PoleZeroConfig {
    pub start_hz: f64,
    pub stop_hz: f64,
    pub samples: usize,
    pub denominator_order: usize,
    pub numerator_order: usize,
}

impl Default for PoleZeroConfig {
    fn default() -> Self {
        Self {
            start_hz: 1.0,
            stop_hz: 1.0e9,
            samples: 64,
            denominator_order: 2,
            numerator_order: 2,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PoleZeroResult {
    /// s-plane roots in rad/s, sorted by (real, imaginary).
    pub poles_rad_per_s: Vec<Complex>,
    pub zeros_rad_per_s: Vec<Complex>,
    pub relative_rms_fit_error: f64,
}

fn evaluate_polynomial(coefficients: &[Complex], x: Complex) -> Complex {
    coefficients
        .iter()
        .rev()
        .fold(Complex::default(), |v, &c| v * x + c)
}

/// Durand–Kerner roots (coefficients in ascending powers), as upstream.
fn polynomial_roots(mut c: Vec<Complex>) -> Vec<Complex> {
    let maximum = c.iter().map(|v| v.norm()).fold(0.0, f64::max);
    if maximum == 0.0 {
        return Vec::new();
    }
    while c.len() > 1 && c.last().is_some_and(|v| v.norm() < maximum * 1.0e-10) {
        c.pop();
    }
    let degree = c.len() - 1;
    if degree == 0 {
        return Vec::new();
    }
    if degree == 1 {
        return vec![-c[0] / c[1]];
    }
    let lead = c[degree];
    for v in c.iter_mut() {
        *v /= lead;
    }
    let mut radius = 1.0f64;
    for v in &c[..degree] {
        radius = radius.max(1.0 + v.norm());
    }
    let mut roots: Vec<Complex> = (0..degree)
        .map(|i| Complex::from_polar(radius, 2.0 * PI * (i as f64 + 0.37) / degree as f64))
        .collect();
    for _ in 0..300 {
        let old = roots.clone();
        let mut change = 0.0f64;
        for i in 0..degree {
            let mut denominator = Complex::new(1.0, 0.0);
            for j in 0..degree {
                if j != i {
                    denominator *= old[i] - old[j];
                }
            }
            if denominator.norm() < 1.0e-24 {
                denominator = Complex::new(1.0e-24, 0.0);
            }
            roots[i] = old[i] - evaluate_polynomial(&c, old[i]) / denominator;
            change = change.max((roots[i] - old[i]).norm());
        }
        if change < 1.0e-12 {
            break;
        }
    }
    roots
}

/// Complex rational least-squares fit (ridge-regularized normal equations) of the AC response
/// at `output`, with roots of the fitted numerator/denominator as zeros/poles.
pub fn pole_zero_analysis(
    circuit: &Circuit,
    output: &str,
    controls: &PoleZeroConfig,
    cfg: &NewtonConfig,
) -> Result<PoleZeroResult> {
    let c = controls;
    require(
        c.start_hz > 0.0
            && c.stop_hz > c.start_hz
            && c.denominator_order > 0
            && c.samples >= 3
            && c.denominator_order <= 8
            && c.numerator_order <= 8,
        "invalid pole-zero controls",
    )?;
    let unknowns = c.denominator_order + c.numerator_order + 1;
    require(
        c.samples > unknowns,
        "pole-zero fit needs more samples than coefficients",
    )?;
    let sweep = circuit.ac_log_sweep(c.start_hz, c.stop_hz, c.samples, cfg)?;
    let scale = 2.0 * PI * (c.start_hz * c.stop_hz).sqrt();
    let mut rows = Vec::with_capacity(sweep.len());
    let mut targets = Vec::with_capacity(sweep.len());
    for sample in &sweep {
        let h = circuit.ac_voltage(sample, output)?;
        let x = Complex::new(0.0, 2.0 * PI * sample.frequency_hz / scale);
        let mut row = vec![Complex::default(); unknowns];
        let mut power = x;
        for j in 0..c.denominator_order {
            row[j] = h * power;
            power *= x;
        }
        power = Complex::new(1.0, 0.0);
        for j in 0..=c.numerator_order {
            row[c.denominator_order + j] = -power;
            power *= x;
        }
        rows.push(row);
        targets.push(-h);
    }
    let mut normal = vec![Complex::default(); unknowns * unknowns];
    let mut rhs = vec![Complex::default(); unknowns];
    for (row, target) in rows.iter().zip(&targets) {
        for i in 0..unknowns {
            rhs[i] += row[i].conj() * target;
            for j in 0..unknowns {
                normal[i * unknowns + j] += row[i].conj() * row[j];
            }
        }
    }
    let mut largest = 0.0f64;
    for i in 0..unknowns {
        largest = largest.max(normal[i * unknowns + i].norm());
    }
    let ridge = largest.max(1.0) * 1.0e-14;
    for i in 0..unknowns {
        normal[i * unknowns + i] += ridge;
    }
    let coefficients = solve_dense(normal, rhs, unknowns, "singular rational-fit normal matrix")?;
    let mut denominator = vec![Complex::default(); c.denominator_order + 1];
    denominator[0] = Complex::new(1.0, 0.0);
    denominator[1..].copy_from_slice(&coefficients[..c.denominator_order]);
    let numerator = coefficients[c.denominator_order..].to_vec();
    let (mut error2, mut signal2) = (0.0, 0.0);
    for sample in &sweep {
        let x = Complex::new(0.0, 2.0 * PI * sample.frequency_hz / scale);
        let fit = evaluate_polynomial(&numerator, x) / evaluate_polynomial(&denominator, x);
        let actual = circuit.ac_voltage(sample, output)?;
        error2 += (fit - actual).norm_sqr();
        signal2 += actual.norm_sqr();
    }
    let order = |a: &Complex, b: &Complex| {
        if a.re == b.re {
            a.im.total_cmp(&b.im)
        } else {
            a.re.total_cmp(&b.re)
        }
    };
    let mut poles: Vec<Complex> = polynomial_roots(denominator)
        .into_iter()
        .map(|r| r * scale)
        .collect();
    let mut zeros: Vec<Complex> = polynomial_roots(numerator)
        .into_iter()
        .map(|r| r * scale)
        .collect();
    poles.sort_by(order);
    zeros.sort_by(order);
    Ok(PoleZeroResult {
        poles_rad_per_s: poles,
        zeros_rad_per_s: zeros,
        relative_rms_fit_error: (error2 / signal2.max(f64::MIN_POSITIVE)).sqrt(),
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeriodicSteadyStateConfig {
    pub period_s: f64,
    pub samples_per_period: usize,
    pub max_periods: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
    pub method: TransientMethod,
}

impl Default for PeriodicSteadyStateConfig {
    fn default() -> Self {
        Self {
            period_s: 0.0,
            samples_per_period: 128,
            max_periods: 64,
            relative_tolerance: 1.0e-4,
            absolute_tolerance: 1.0e-7,
            method: TransientMethod::Trapezoidal,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeriodicSteadyStateResult {
    /// Time shifted to [0, period].
    pub period: Vec<TransientPoint>,
    pub periods: usize,
    pub normalized_residual: f64,
    pub converged: bool,
}

/// Shooting-by-settling: integrate `max_periods` periods and return the first period whose
/// node voltages match the previous period at every phase sample within tolerance.
pub fn periodic_steady_state(
    circuit: &Circuit,
    controls: &PeriodicSteadyStateConfig,
    cfg: &NewtonConfig,
) -> Result<PeriodicSteadyStateResult> {
    let c = controls;
    require(
        c.period_s > 0.0
            && c.period_s.is_finite()
            && c.samples_per_period >= 4
            && c.max_periods >= 2
            && c.relative_tolerance > 0.0
            && c.absolute_tolerance > 0.0,
        "invalid periodic steady-state controls",
    )?;
    let total = c
        .max_periods
        .checked_mul(c.samples_per_period)
        .ok_or_else(|| {
            SolverError::InvalidInput("periodic steady-state step count overflow".into())
        })?;
    let dt = c.period_s / c.samples_per_period as f64;
    let transient = circuit.transient(dt, total, c.method, cfg)?;
    let spp = c.samples_per_period;
    let mut result = PeriodicSteadyStateResult {
        normalized_residual: f64::INFINITY,
        ..Default::default()
    };
    let mut selected = c.max_periods - 1;
    for period in 1..c.max_periods {
        let base = period * spp;
        let mut residual = 0.0f64;
        for phase in 0..=spp {
            let cur = &transient[base + phase];
            let prev = &transient[base - spp + phase];
            if cur.node_voltage.len() != prev.node_voltage.len() {
                return Err(err_runtime("PSS node-vector size changed"));
            }
            for node in 1..cur.node_voltage.len() {
                let (a, b) = (cur.node_voltage[node], prev.node_voltage[node]);
                let scale = c.absolute_tolerance + c.relative_tolerance * a.abs().max(b.abs());
                residual = residual.max((a - b).abs() / scale);
            }
        }
        result.normalized_residual = residual;
        if residual <= 1.0 {
            result.converged = true;
            result.periods = period + 1;
            selected = period;
            break;
        }
    }
    if !result.converged {
        result.periods = c.max_periods;
    }
    let base = selected * spp;
    let origin = transient[base].time;
    result.period = transient[base..=base + spp]
        .iter()
        .map(|p| TransientPoint {
            time: p.time - origin,
            ..p.clone()
        })
        .collect();
    Ok(result)
}
