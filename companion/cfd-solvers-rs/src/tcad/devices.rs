//! Physics-parameterized BJT, IGBT, Gummel-Poon, MOS level 2/3, VDMOS and hysteretic-core models.
use super::{
    fd_evaluation, silicon_bandgap_default, temperature_scaled_intrinsic_density,
    temperature_scaled_mobility, thermal_voltage, LongChannelMosfet, LongChannelMosfetConfig,
    SemiconductorMaterial, StaticDeviceEvaluator, ELEMENTARY_CHARGE_C,
};
use crate::{require, Result};
use std::sync::Arc;

fn safe_expm1(x: f64) -> f64 {
    x.clamp(-80.0, 80.0).exp_m1()
}

fn positive_coth(x: f64) -> Result<f64> {
    require(x > 0.0, "coth argument must be positive")?;
    Ok(if x < 1e-4 {
        1.0 / x + x / 3.0
    } else if x > 30.0 {
        1.0
    } else {
        1.0 / x.tanh()
    })
}

fn bounded_cosh_inverse(x: f64) -> f64 {
    if x > 40.0 {
        0.0
    } else {
        1.0 / x.cosh()
    }
}

fn material_at_temperature(r: &SemiconductorMaterial, t: f64) -> Result<SemiconductorMaterial> {
    require(
        t > 0.0 && t.is_finite(),
        "semiconductor temperature must be positive and finite",
    )?;
    let mut m = *r;
    if (t - r.temperature_k).abs() < 1e-12 {
        return Ok(m);
    }
    let eg_ref = silicon_bandgap_default(r.temperature_k)?;
    let eg = silicon_bandgap_default(t)?;
    m.intrinsic_density_m3 = temperature_scaled_intrinsic_density(
        r.intrinsic_density_m3,
        r.temperature_k,
        t,
        eg_ref,
        eg,
    )?;
    m.electron_mobility_m2_per_vs =
        temperature_scaled_mobility(r.electron_mobility_m2_per_vs, r.temperature_k, t, 2.42)?;
    m.hole_mobility_m2_per_vs =
        temperature_scaled_mobility(r.hole_mobility_m2_per_vs, r.temperature_k, t, 2.20)?;
    m.conduction_band_density_m3 = r.conduction_band_density_m3 * (t / r.temperature_k).powf(1.5);
    m.valence_band_density_m3 = r.valence_band_density_m3 * (t / r.temperature_k).powf(1.5);
    m.temperature_k = t;
    Ok(m)
}

fn built_in_voltage(a: f64, b: f64, ni: f64, vt: f64) -> f64 {
    vt * (a * b / (ni * ni)).max(1.0).ln()
}

fn finite_difference_step(x: f64, base: f64) -> f64 {
    base.max(x.abs() * 1e-7)
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum BipolarPolarity {
    #[default]
    Npn,
    Pnp,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BipolarJunctionConfig {
    pub polarity: BipolarPolarity,
    pub area_m2: f64,
    pub emitter_width_m: f64,
    pub base_width_m: f64,
    pub collector_width_m: f64,
    pub emitter_doping_m3: f64,
    pub base_doping_m3: f64,
    pub collector_doping_m3: f64,
    pub minority_electron_lifetime_base_s: f64,
    pub minority_hole_lifetime_emitter_s: f64,
    pub minority_hole_lifetime_collector_s: f64,
    pub material: SemiconductorMaterial,
}

impl Default for BipolarJunctionConfig {
    fn default() -> Self {
        Self {
            polarity: BipolarPolarity::Npn,
            area_m2: 1e-10,
            emitter_width_m: 0.5e-6,
            base_width_m: 0.25e-6,
            collector_width_m: 2e-6,
            emitter_doping_m3: 1e24,
            base_doping_m3: 2e22,
            collector_doping_m3: 2e21,
            minority_electron_lifetime_base_s: 2e-7,
            minority_hole_lifetime_emitter_s: 5e-8,
            minority_hole_lifetime_collector_s: 2e-7,
            material: SemiconductorMaterial::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BipolarDerivedParameters {
    pub thermal_voltage_v: f64,
    pub electron_diffusion_length_base_m: f64,
    pub hole_diffusion_length_emitter_m: f64,
    pub hole_diffusion_length_collector_m: f64,
    pub base_transport_factor: f64,
    pub emitter_injection_efficiency: f64,
    pub collector_injection_efficiency: f64,
    pub forward_alpha: f64,
    pub reverse_alpha: f64,
    pub forward_beta: f64,
    pub reverse_beta: f64,
    pub emitter_saturation_current_a: f64,
    pub collector_saturation_current_a: f64,
    pub transport_saturation_current_a: f64,
    pub emitter_base_built_in_v: f64,
    pub collector_base_built_in_v: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BipolarOperatingPoint {
    pub collector_voltage_v: f64,
    pub base_voltage_v: f64,
    pub emitter_voltage_v: f64,
    pub collector_current_a: f64,
    pub base_current_a: f64,
    pub emitter_current_a: f64,
    pub dissipated_power_w: f64,
    pub parameters: BipolarDerivedParameters,
}

/// Abrupt emitter/base/collector profile with Ebers-Moll transport derived from minority-carrier
/// diffusion. A temperature of 0 means the material temperature (upstream convention).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BipolarJunctionTransistor1D {
    config: BipolarJunctionConfig,
}

impl BipolarJunctionTransistor1D {
    pub fn new(config: BipolarJunctionConfig) -> Result<Self> {
        let c = &config;
        let m = &c.material;
        require(
            c.area_m2 > 0.0
                && c.emitter_width_m > 0.0
                && c.base_width_m > 0.0
                && c.collector_width_m > 0.0
                && c.emitter_doping_m3 > 0.0
                && c.base_doping_m3 > 0.0
                && c.collector_doping_m3 > 0.0
                && c.minority_electron_lifetime_base_s > 0.0
                && c.minority_hole_lifetime_emitter_s > 0.0
                && c.minority_hole_lifetime_collector_s > 0.0
                && m.temperature_k > 0.0
                && m.intrinsic_density_m3 > 0.0
                && m.electron_mobility_m2_per_vs > 0.0
                && m.hole_mobility_m2_per_vs > 0.0,
            "invalid bipolar-junction configuration",
        )?;
        Ok(Self { config })
    }
    pub fn config(&self) -> &BipolarJunctionConfig {
        &self.config
    }
    pub fn total_length_m(&self) -> f64 {
        self.config.emitter_width_m + self.config.base_width_m + self.config.collector_width_m
    }
    pub fn net_doping_m3(&self, x_m: f64) -> Result<f64> {
        let length = self.total_length_m();
        let tol = 16.0 * f64::EPSILON * length.max(1.0);
        require(
            x_m.is_finite() && x_m >= -tol && x_m <= length + tol,
            "BJT profile coordinate outside device",
        )?;
        let x = x_m.clamp(0.0, length);
        let c = &self.config;
        let npn = c.polarity == BipolarPolarity::Npn;
        Ok(if x < c.emitter_width_m {
            (if npn { 1.0 } else { -1.0 }) * c.emitter_doping_m3
        } else if x < c.emitter_width_m + c.base_width_m {
            (if npn { -1.0 } else { 1.0 }) * c.base_doping_m3
        } else {
            (if npn { 1.0 } else { -1.0 }) * c.collector_doping_m3
        })
    }
    pub fn doping_profile(&self, samples: usize) -> Result<Vec<f64>> {
        require(
            samples >= 3,
            "BJT doping profile requires at least three samples",
        )?;
        let l = self.total_length_m();
        (0..samples)
            .map(|i| self.net_doping_m3(l * i as f64 / (samples - 1) as f64))
            .collect()
    }
    pub fn derived(&self, temperature_k: f64) -> Result<BipolarDerivedParameters> {
        let c = &self.config;
        let t = if temperature_k == 0.0 {
            c.material.temperature_k
        } else {
            temperature_k
        };
        let m = material_at_temperature(&c.material, t)?;
        let vt = thermal_voltage(t)?;
        let dn = m.electron_mobility_m2_per_vs * vt;
        let dp_e = m.hole_mobility_m2_per_vs * vt;
        let dp_c = m.hole_mobility_m2_per_vs * vt;
        let ln = (dn * c.minority_electron_lifetime_base_s).sqrt();
        let lp_e = (dp_e * c.minority_hole_lifetime_emitter_s).sqrt();
        let lp_c = (dp_c * c.minority_hole_lifetime_collector_s).sqrt();
        let wb = c.base_width_m / ln;
        let alpha_t = bounded_cosh_inverse(wb);
        let ni2 = m.intrinsic_density_m3 * m.intrinsic_density_m3;
        let jn = ELEMENTARY_CHARGE_C * dn * ni2 / (c.base_doping_m3 * ln) * positive_coth(wb)?;
        let jp_e = ELEMENTARY_CHARGE_C * dp_e * ni2 / (c.emitter_doping_m3 * lp_e)
            * positive_coth(c.emitter_width_m / lp_e)?;
        let jp_c = ELEMENTARY_CHARGE_C * dp_c * ni2 / (c.collector_doping_m3 * lp_c)
            * positive_coth(c.collector_width_m / lp_c)?;
        let gamma_e = jn / (jn + jp_e);
        let gamma_c = jn / (jn + jp_c);
        let alpha_f = (alpha_t * gamma_e).clamp(1e-12, 1.0 - 1e-12);
        let alpha_r = (alpha_t * gamma_c).clamp(1e-12, 1.0 - 1e-12);
        Ok(BipolarDerivedParameters {
            thermal_voltage_v: vt,
            electron_diffusion_length_base_m: ln,
            hole_diffusion_length_emitter_m: lp_e,
            hole_diffusion_length_collector_m: lp_c,
            base_transport_factor: alpha_t,
            emitter_injection_efficiency: gamma_e,
            collector_injection_efficiency: gamma_c,
            forward_alpha: alpha_f,
            reverse_alpha: alpha_r,
            forward_beta: alpha_f / (1.0 - alpha_f),
            reverse_beta: alpha_r / (1.0 - alpha_r),
            emitter_saturation_current_a: c.area_m2 * (jn + jp_e),
            collector_saturation_current_a: c.area_m2 * (jn + jp_c),
            transport_saturation_current_a: c.area_m2 * alpha_t * jn,
            emitter_base_built_in_v: built_in_voltage(
                c.emitter_doping_m3,
                c.base_doping_m3,
                m.intrinsic_density_m3,
                vt,
            ),
            collector_base_built_in_v: built_in_voltage(
                c.collector_doping_m3,
                c.base_doping_m3,
                m.intrinsic_density_m3,
                vt,
            ),
        })
    }
    pub fn operating_point(
        &self,
        collector_v: f64,
        base_v: f64,
        emitter_v: f64,
        temperature_k: f64,
    ) -> Result<BipolarOperatingPoint> {
        require(
            collector_v.is_finite() && base_v.is_finite() && emitter_v.is_finite(),
            "BJT terminal voltages must be finite",
        )?;
        let p = self.derived(temperature_k)?;
        let pol = if self.config.polarity == BipolarPolarity::Npn {
            1.0
        } else {
            -1.0
        };
        let vbe = pol * (base_v - emitter_v);
        let vbc = pol * (base_v - collector_v);
        let forward = p.emitter_saturation_current_a * safe_expm1(vbe / p.thermal_voltage_v);
        let reverse = p.collector_saturation_current_a * safe_expm1(vbc / p.thermal_voltage_v);
        let ic_eff = p.forward_alpha * forward - reverse;
        let ie_eff = -forward + p.reverse_alpha * reverse;
        let ib_eff = -(ic_eff + ie_eff);
        let (ic, ib, ie) = (pol * ic_eff, pol * ib_eff, pol * ie_eff);
        Ok(BipolarOperatingPoint {
            collector_voltage_v: collector_v,
            base_voltage_v: base_v,
            emitter_voltage_v: emitter_v,
            collector_current_a: ic,
            base_current_a: ib,
            emitter_current_a: ie,
            dissipated_power_w: collector_v * ic + base_v * ib + emitter_v * ie,
            parameters: p,
        })
    }
}

/// Three-terminal `[C, B, E]` evaluator (upstream default step 1e-6 V).
pub fn make_bipolar_tcad_spice_evaluator(
    config: BipolarJunctionConfig,
    jacobian_step_v: f64,
) -> Result<StaticDeviceEvaluator> {
    require(
        jacobian_step_v > 0.0 && jacobian_step_v.is_finite(),
        "BJT Jacobian step must be positive",
    )?;
    let model = BipolarJunctionTransistor1D::new(config)?;
    let h = jacobian_step_v;
    Ok(Arc::new(move |v: &[f64]| {
        require(v.len() == 3, "TCAD BJT evaluator expects [C,B,E]")?;
        let currents = |x: &[f64]| -> Result<Vec<f64>> {
            let op = model.operating_point(x[0], x[1], x[2], 0.0)?;
            Ok(vec![
                op.collector_current_a,
                op.base_current_a,
                op.emitter_current_a,
            ])
        };
        fd_evaluation(v, 3, &currents, &|x| finite_difference_step(x, h))
    }))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IgbtConfig {
    pub channel: LongChannelMosfetConfig,
    /// Low-gain PNP path for drift-region conductivity modulation (polarity forced to PNP).
    pub bipolar: BipolarJunctionConfig,
    pub drift_resistance_ohm: f64,
    pub junction_drop_v: f64,
    pub conductivity_modulation_current_a: f64,
    pub maximum_bipolar_gain: f64,
    pub tail_lifetime_s: f64,
    pub ambient_temperature_k: f64,
    pub thermal_resistance_k_per_w: f64,
    pub thermal_capacitance_j_per_k: f64,
    pub drift_resistance_temperature_coefficient_per_k: f64,
}

impl Default for IgbtConfig {
    fn default() -> Self {
        Self {
            channel: LongChannelMosfetConfig::default(),
            bipolar: BipolarJunctionConfig::default(),
            drift_resistance_ohm: 0.08,
            junction_drop_v: 0.65,
            conductivity_modulation_current_a: 0.05,
            maximum_bipolar_gain: 12.0,
            tail_lifetime_s: 2e-6,
            ambient_temperature_k: 300.0,
            thermal_resistance_k_per_w: 0.6,
            thermal_capacitance_j_per_k: 0.02,
            drift_resistance_temperature_coefficient_per_k: 3.5e-3,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IgbtOperatingPoint {
    pub gate_emitter_voltage_v: f64,
    pub collector_emitter_voltage_v: f64,
    pub collector_current_a: f64,
    pub channel_current_a: f64,
    pub bipolar_current_a: f64,
    pub effective_drift_resistance_ohm: f64,
    pub junction_temperature_k: f64,
    pub dissipated_power_w: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IgbtTransientState {
    pub junction_temperature_k: f64,
    pub tail_current_a: f64,
    pub time_s: f64,
}

impl Default for IgbtTransientState {
    fn default() -> Self {
        Self {
            junction_temperature_k: 300.0,
            tail_current_a: 0.0,
            time_s: 0.0,
        }
    }
}

/// MOS-gated bipolar power device: long-channel MOS channel plus BJT gain, drift-region
/// conductivity modulation, turn-off tail and lumped thermal RC (not a field-resolved IGBT).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IgbtPowerDevice {
    config: IgbtConfig,
    state: IgbtTransientState,
}

impl IgbtPowerDevice {
    pub fn new(mut config: IgbtConfig) -> Result<Self> {
        let c = &config;
        require(
            c.drift_resistance_ohm > 0.0
                && c.junction_drop_v >= 0.0
                && c.conductivity_modulation_current_a > 0.0
                && c.maximum_bipolar_gain >= 0.0
                && c.tail_lifetime_s > 0.0
                && c.ambient_temperature_k > 0.0
                && c.thermal_resistance_k_per_w > 0.0
                && c.thermal_capacitance_j_per_k > 0.0
                && c.drift_resistance_temperature_coefficient_per_k >= 0.0,
            "invalid IGBT configuration",
        )?;
        config.bipolar.polarity = BipolarPolarity::Pnp;
        Ok(Self {
            config,
            state: IgbtTransientState {
                junction_temperature_k: config.ambient_temperature_k,
                ..Default::default()
            },
        })
    }
    pub fn state(&self) -> &IgbtTransientState {
        &self.state
    }
    pub fn config(&self) -> &IgbtConfig {
        &self.config
    }
    /// Returns (total, channel, bipolar, effective drift resistance).
    fn static_current(&self, vge: f64, vce: f64, t: f64) -> Result<(f64, f64, f64, f64)> {
        let c = &self.config;
        let (mut channel, mut bipolar): (f64, f64) = (0.0, 0.0);
        if vce.is_nan() || vce <= c.junction_drop_v {
            return Ok((0.0, 0.0, 0.0, c.drift_resistance_ohm));
        }
        let bjt = BipolarJunctionTransistor1D::new(c.bipolar)?.derived(t)?;
        let beta = c.maximum_bipolar_gain.min(bjt.forward_beta);
        let scale = temperature_scaled_mobility(1.0, c.ambient_temperature_k, t, 2.42)?;
        let mut ch = c.channel;
        ch.electron_mobility_m2_per_vs *= scale;
        let mos = LongChannelMosfet::new(ch)?;
        let r_temp = c.drift_resistance_ohm
            * (1.0
                + c.drift_resistance_temperature_coefficient_per_k * (t - c.ambient_temperature_k));
        let mut total: f64 = 0.0;
        let mut resistance = c.drift_resistance_ohm;
        for _ in 0..100 {
            let modulation = 1.0 + channel.abs() / c.conductivity_modulation_current_a;
            resistance = (r_temp / modulation).max(1e-6);
            let vds = (vce - c.junction_drop_v - total * resistance).max(0.0);
            let new_channel = mos.drain_current(vge, vds, 0.0)?;
            let new_bipolar = beta * new_channel;
            let next = 0.5 * total + 0.5 * (new_channel + new_bipolar);
            channel = new_channel;
            bipolar = new_bipolar;
            if (next - total).abs() <= 1e-10 * next.abs().max(1.0) {
                total = next;
                break;
            }
            total = next;
        }
        Ok((total.max(0.0), channel, bipolar, resistance))
    }
    /// Temperature 0 means ambient (upstream convention).
    pub fn operating_point(
        &self,
        vge_v: f64,
        vce_v: f64,
        temperature_k: f64,
    ) -> Result<IgbtOperatingPoint> {
        require(
            vge_v.is_finite() && vce_v.is_finite(),
            "IGBT voltages must be finite",
        )?;
        let t = if temperature_k == 0.0 {
            self.config.ambient_temperature_k
        } else {
            temperature_k
        };
        let (current, channel, bipolar, resistance) = self.static_current(vge_v, vce_v, t)?;
        Ok(IgbtOperatingPoint {
            gate_emitter_voltage_v: vge_v,
            collector_emitter_voltage_v: vce_v,
            collector_current_a: current,
            channel_current_a: channel,
            bipolar_current_a: bipolar,
            effective_drift_resistance_ohm: resistance,
            junction_temperature_k: t,
            dissipated_power_w: current * vce_v.max(0.0),
        })
    }
    pub fn step(&mut self, vge_v: f64, vce_v: f64, dt_s: f64) -> Result<IgbtOperatingPoint> {
        require(
            dt_s > 0.0 && dt_s.is_finite(),
            "IGBT timestep must be positive",
        )?;
        let mut op = self.operating_point(vge_v, vce_v, self.state.junction_temperature_k)?;
        let c = &self.config;
        let mut s = self.state;
        let target = op.bipolar_current_a;
        let decay = (-dt_s / c.tail_lifetime_s).exp();
        if op.channel_current_a > 0.0 {
            s.tail_current_a = target + (s.tail_current_a - target) * decay;
        } else {
            s.tail_current_a *= decay;
        }
        let conduction = op.channel_current_a + s.tail_current_a;
        let power = conduction * vce_v.max(0.0);
        let cooling =
            (s.junction_temperature_k - c.ambient_temperature_k) / c.thermal_resistance_k_per_w;
        s.junction_temperature_k += dt_s * (power - cooling) / c.thermal_capacitance_j_per_k;
        s.junction_temperature_k = s.junction_temperature_k.max(1.0);
        s.time_s += dt_s;
        require(
            s.junction_temperature_k.is_finite() && s.tail_current_a.is_finite(),
            "IGBT state overflow",
        )?;
        self.state = s;
        op.collector_current_a = conduction;
        op.bipolar_current_a = s.tail_current_a;
        op.junction_temperature_k = s.junction_temperature_k;
        op.dissipated_power_w = power;
        Ok(op)
    }
    /// Temperature 0 means ambient.
    pub fn reset(&mut self, temperature_k: f64) -> Result<()> {
        let t = if temperature_k == 0.0 {
            self.config.ambient_temperature_k
        } else {
            temperature_k
        };
        require(
            t > 0.0 && t.is_finite(),
            "IGBT reset temperature must be positive",
        )?;
        self.state = IgbtTransientState {
            junction_temperature_k: t,
            tail_current_a: 0.0,
            time_s: 0.0,
        };
        Ok(())
    }
}

/// Three-terminal `[C, G, E]` evaluator (upstream default step 1e-6 V).
pub fn make_igbt_tcad_spice_evaluator(
    config: IgbtConfig,
    jacobian_step_v: f64,
) -> Result<StaticDeviceEvaluator> {
    require(
        jacobian_step_v > 0.0 && jacobian_step_v.is_finite(),
        "IGBT Jacobian step must be positive",
    )?;
    let model = IgbtPowerDevice::new(config)?;
    let h = jacobian_step_v;
    Ok(Arc::new(move |v: &[f64]| {
        require(v.len() == 3, "TCAD IGBT evaluator expects [C,G,E]")?;
        let currents = |x: &[f64]| -> Result<Vec<f64>> {
            let op = model.operating_point(x[1] - x[2], x[0] - x[2], 0.0)?;
            Ok(vec![op.collector_current_a, 0.0, -op.collector_current_a])
        };
        fd_evaluation(v, 3, &currents, &|x| finite_difference_step(x, h))
    }))
}

const BOLTZMANN_OVER_Q: f64 = 8.617333262145e-5;

fn finite_positive(v: f64, msg: &str) -> Result<()> {
    require(v > 0.0 && v.is_finite(), msg)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GummelPoonConfig {
    pub saturation_current_a: f64,
    pub forward_beta: f64,
    pub reverse_beta: f64,
    pub early_voltage_forward_v: f64,
    pub early_voltage_reverse_v: f64,
    pub high_injection_current_a: f64,
    pub temperature_k: f64,
}

impl Default for GummelPoonConfig {
    fn default() -> Self {
        Self {
            saturation_current_a: 1e-15,
            forward_beta: 120.0,
            reverse_beta: 2.0,
            early_voltage_forward_v: 80.0,
            early_voltage_reverse_v: 20.0,
            high_injection_current_a: 0.02,
            temperature_k: 300.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GummelPoonOperatingPoint {
    pub collector_current_a: f64,
    pub base_current_a: f64,
    pub emitter_current_a: f64,
    pub forward_transport_current_a: f64,
    pub reverse_transport_current_a: f64,
}

/// Charge-control Gummel-Poon baseline with Early effect and smooth high-injection roll-off.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GummelPoonBjt {
    config: GummelPoonConfig,
}

impl GummelPoonBjt {
    pub fn new(config: GummelPoonConfig) -> Result<Self> {
        finite_positive(
            config.saturation_current_a,
            "invalid Gummel-Poon saturation current",
        )?;
        finite_positive(config.forward_beta, "invalid Gummel-Poon forward beta")?;
        finite_positive(config.reverse_beta, "invalid Gummel-Poon reverse beta")?;
        finite_positive(
            config.early_voltage_forward_v,
            "invalid Gummel-Poon Early voltage",
        )?;
        finite_positive(
            config.early_voltage_reverse_v,
            "invalid reverse Early voltage",
        )?;
        require(
            config.high_injection_current_a > 0.0
                && config.temperature_k > 0.0
                && config.high_injection_current_a.is_finite()
                && config.temperature_k.is_finite(),
            "invalid Gummel-Poon controls",
        )?;
        Ok(Self { config })
    }
    pub fn operating_point(&self, vbe_v: f64, vbc_v: f64) -> Result<GummelPoonOperatingPoint> {
        require(
            vbe_v.is_finite() && vbc_v.is_finite(),
            "invalid BJT terminal voltage",
        )?;
        let c = &self.config;
        let vt = BOLTZMANN_OVER_Q * c.temperature_k;
        let fr = c.saturation_current_a * (vbe_v / vt).clamp(-80.0, 80.0).exp_m1();
        let rr = c.saturation_current_a * (vbc_v / vt).clamp(-80.0, 80.0).exp_m1();
        let forward = fr / (1.0 + fr.abs() / c.high_injection_current_a);
        let reverse = rr / (1.0 + rr.abs() / c.high_injection_current_a);
        let af = c.forward_beta / (c.forward_beta + 1.0);
        let ar = c.reverse_beta / (c.reverse_beta + 1.0);
        let ef = 1.0 + (vbe_v - vbc_v) / c.early_voltage_forward_v;
        let er = 1.0 + (vbc_v - vbe_v) / c.early_voltage_reverse_v;
        let ft = af * forward * ef.max(0.05);
        let rt = ar * reverse * er.max(0.05);
        let collector = ft - reverse;
        let base = (1.0 - af) * forward + (1.0 - ar) * reverse;
        Ok(GummelPoonOperatingPoint {
            collector_current_a: collector,
            base_current_a: base,
            emitter_current_a: -collector - base,
            forward_transport_current_a: ft,
            reverse_transport_current_a: rt,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MosLevel23Config {
    pub oxide_capacitance_f_per_m2: f64,
    pub width_m: f64,
    pub length_m: f64,
    pub mobility_m2_per_vs: f64,
    pub threshold_voltage_v: f64,
    pub channel_length_modulation_per_v: f64,
    pub mobility_degradation_per_v: f64,
    pub subthreshold_factor: f64,
    pub subthreshold_current_a: f64,
    pub temperature_k: f64,
}

impl Default for MosLevel23Config {
    fn default() -> Self {
        Self {
            oxide_capacitance_f_per_m2: 3.45e-3,
            width_m: 10e-6,
            length_m: 1e-6,
            mobility_m2_per_vs: 0.05,
            threshold_voltage_v: 0.5,
            channel_length_modulation_per_v: 0.02,
            mobility_degradation_per_v: 0.08,
            subthreshold_factor: 1.5,
            subthreshold_current_a: 1e-12,
            temperature_k: 300.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MosLevel23OperatingPoint {
    pub drain_current_a: f64,
    pub transconductance_s: f64,
    pub output_conductance_s: f64,
    pub subthreshold: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MosLevel23 {
    config: MosLevel23Config,
}

impl MosLevel23 {
    pub fn new(config: MosLevel23Config) -> Result<Self> {
        finite_positive(
            config.oxide_capacitance_f_per_m2,
            "invalid MOS oxide capacitance",
        )?;
        finite_positive(config.width_m, "invalid MOS width")?;
        finite_positive(config.length_m, "invalid MOS length")?;
        finite_positive(config.mobility_m2_per_vs, "invalid MOS mobility")?;
        require(
            config.subthreshold_factor > 1.0
                && config.subthreshold_current_a > 0.0
                && config.temperature_k > 0.0
                && config.threshold_voltage_v.is_finite()
                && config.channel_length_modulation_per_v.is_finite()
                && config.mobility_degradation_per_v.is_finite(),
            "invalid MOS level 2/3 controls",
        )?;
        Ok(Self { config })
    }
    pub fn operating_point(&self, vgs_v: f64, vds_v: f64) -> Result<MosLevel23OperatingPoint> {
        require(
            vgs_v.is_finite() && vds_v.is_finite(),
            "invalid MOS terminal voltage",
        )?;
        let c = &self.config;
        let vt = BOLTZMANN_OVER_Q * c.temperature_k;
        let overdrive = vgs_v - c.threshold_voltage_v;
        if overdrive <= 0.0 {
            let e = (overdrive / (c.subthreshold_factor * vt)).clamp(-80.0, 80.0);
            let current =
                c.subthreshold_current_a * e.exp() * (1.0 - (-vds_v / vt).clamp(-80.0, 80.0).exp());
            return Ok(MosLevel23OperatingPoint {
                drain_current_a: current,
                transconductance_s: current / (c.subthreshold_factor * vt),
                output_conductance_s: c.subthreshold_current_a * e.exp() / vt,
                subthreshold: true,
            });
        }
        let mobility = c.mobility_m2_per_vs / (1.0 + c.mobility_degradation_per_v * overdrive);
        let beta = mobility * c.oxide_capacitance_f_per_m2 * c.width_m / c.length_m;
        let vds = vds_v.max(0.0);
        let sat = vds.min(overdrive);
        let lm = 1.0 + c.channel_length_modulation_per_v * vds;
        let current = beta * (overdrive * sat - 0.5 * sat * sat) * lm;
        let gm = beta * (overdrive - sat).max(0.0) + beta * sat * lm;
        let go = current * c.channel_length_modulation_per_v / lm.max(1e-30);
        Ok(MosLevel23OperatingPoint {
            drain_current_a: current,
            transconductance_s: gm,
            output_conductance_s: go,
            subthreshold: false,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VdmosConfig {
    pub channel: MosLevel23Config,
    pub drift_resistance_ohm: f64,
    pub body_diode_saturation_a: f64,
    pub body_diode_ideality: f64,
}

impl Default for VdmosConfig {
    fn default() -> Self {
        Self {
            channel: MosLevel23Config::default(),
            drift_resistance_ohm: 0.08,
            body_diode_saturation_a: 1e-12,
            body_diode_ideality: 1.4,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct VdmosOperatingPoint {
    pub drain_current_a: f64,
    pub channel_current_a: f64,
    pub body_diode_current_a: f64,
    pub drift_drop_v: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VdmosPowerDevice {
    config: VdmosConfig,
    channel: MosLevel23,
}

impl VdmosPowerDevice {
    pub fn new(config: VdmosConfig) -> Result<Self> {
        let channel = MosLevel23::new(config.channel)?;
        finite_positive(
            config.drift_resistance_ohm,
            "invalid VDMOS drift resistance",
        )?;
        finite_positive(config.body_diode_saturation_a, "invalid VDMOS body diode")?;
        finite_positive(config.body_diode_ideality, "invalid VDMOS diode ideality")?;
        Ok(Self { config, channel })
    }
    pub fn operating_point(&self, vgs_v: f64, vds_v: f64) -> Result<VdmosOperatingPoint> {
        require(
            vgs_v.is_finite() && vds_v.is_finite(),
            "invalid VDMOS terminal voltage",
        )?;
        if vds_v >= 0.0 {
            let ch = self.channel.operating_point(vgs_v, vds_v)?;
            return Ok(VdmosOperatingPoint {
                drain_current_a: ch.drain_current_a,
                channel_current_a: ch.drain_current_a,
                body_diode_current_a: 0.0,
                drift_drop_v: ch.drain_current_a * self.config.drift_resistance_ohm,
            });
        }
        let vt =
            BOLTZMANN_OVER_Q * self.config.channel.temperature_k * self.config.body_diode_ideality;
        let diode = self.config.body_diode_saturation_a * (-vds_v / vt).clamp(-80.0, 80.0).exp_m1();
        Ok(VdmosOperatingPoint {
            drain_current_a: -diode,
            channel_current_a: 0.0,
            body_diode_current_a: -diode,
            drift_drop_v: 0.0,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HystereticCoreConfig {
    pub saturation_induction_t: f64,
    pub coercive_field_a_per_m: f64,
    pub relaxation_time_s: f64,
    pub temperature_coefficient_per_k: f64,
    pub reference_temperature_k: f64,
}

impl Default for HystereticCoreConfig {
    fn default() -> Self {
        Self {
            saturation_induction_t: 1.5,
            coercive_field_a_per_m: 40.0,
            relaxation_time_s: 1e-4,
            temperature_coefficient_per_k: 0.0,
            reference_temperature_k: 300.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HystereticCoreState {
    pub induction_t: f64,
    pub loss_density_j_per_m3: f64,
    pub previous_field_a_per_m: f64,
}

/// Rate-dependent tanh hysteresis with accumulated loss density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HystereticCoreModel {
    config: HystereticCoreConfig,
    state: HystereticCoreState,
}

impl HystereticCoreModel {
    pub fn new(config: HystereticCoreConfig) -> Result<Self> {
        finite_positive(
            config.saturation_induction_t,
            "invalid core saturation induction",
        )?;
        finite_positive(config.coercive_field_a_per_m, "invalid core coercive field")?;
        finite_positive(config.relaxation_time_s, "invalid core relaxation time")?;
        finite_positive(
            config.reference_temperature_k,
            "invalid core reference temperature",
        )?;
        Ok(Self {
            config,
            state: HystereticCoreState::default(),
        })
    }
    pub fn reset(&mut self) {
        self.state = HystereticCoreState::default();
    }
    pub fn state(&self) -> &HystereticCoreState {
        &self.state
    }
    pub fn step(&mut self, field_a_per_m: f64, temperature_k: f64, dt_s: f64) -> Result<()> {
        require(
            field_a_per_m.is_finite()
                && temperature_k > 0.0
                && dt_s > 0.0
                && temperature_k.is_finite()
                && dt_s.is_finite(),
            "invalid hysteretic core state",
        )?;
        let c = &self.config;
        let tf = (1.0
            - c.temperature_coefficient_per_k * (temperature_k - c.reference_temperature_k))
            .max(0.05);
        let shift = c
            .coercive_field_a_per_m
            .copysign(field_a_per_m - self.state.previous_field_a_per_m);
        let target = c.saturation_induction_t
            * tf
            * ((field_a_per_m - shift) / c.coercive_field_a_per_m).tanh();
        let blend = 1.0 - (-dt_s / c.relaxation_time_s).exp();
        let old = self.state.induction_t;
        self.state.induction_t += blend * (target - self.state.induction_t);
        self.state.loss_density_j_per_m3 += ((self.state.induction_t - old) * field_a_per_m).abs();
        self.state.previous_field_a_per_m = field_a_per_m;
        Ok(())
    }
}
