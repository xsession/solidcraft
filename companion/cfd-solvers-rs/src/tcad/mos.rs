//! Fermi-Dirac statistics, incomplete ionization, the analytic MOS capacitor and the long-channel
//! charge-sheet MOSFET.
use super::{
    fd_evaluation, safe_exp, thermal_voltage, SemiconductorMaterial, StaticDeviceEvaluator,
    ELEMENTARY_CHARGE_C, VACUUM_PERMITTIVITY_F_PER_M,
};
use crate::{require, Result, SolverError};
use std::sync::Arc;

const PI: f64 = std::f64::consts::PI;

fn simpson_fermi_half(eta: f64) -> f64 {
    let upper = 40.0f64.max(eta + 40.0);
    const INTERVALS: usize = 4096;
    let h = upper / INTERVALS as f64;
    let f = |e: f64| {
        if e <= 0.0 {
            return 0.0;
        }
        let z = e - eta;
        let occ = if z > 50.0 {
            (-z).exp()
        } else if z < -50.0 {
            1.0
        } else {
            1.0 / (1.0 + z.exp())
        };
        e.sqrt() * occ
    };
    let mut sum = f(0.0) + f(upper);
    for i in 1..INTERVALS {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(i as f64 * h);
    }
    (2.0 / PI.sqrt()) * (h / 3.0) * sum
}

/// `F_{1/2}(eta)` normalized so that it tends to `exp(eta)` when nondegenerate.
pub fn fermi_dirac_half(eta: f64) -> Result<f64> {
    require(eta.is_finite(), "Fermi-Dirac argument must be finite")?;
    Ok(if eta < -8.0 {
        eta.exp()
    } else if eta > 40.0 {
        (4.0 / (3.0 * PI.sqrt())) * eta.powf(1.5)
    } else {
        simpson_fermi_half(eta)
    })
}

pub fn inverse_fermi_dirac_half(value: f64) -> Result<f64> {
    require(
        value > 0.0 && value.is_finite(),
        "inverse Fermi-Dirac value must be positive and finite",
    )?;
    let (mut lo, mut hi) = (-100.0, 100.0);
    for _ in 0..160 {
        let mid = 0.5 * (lo + hi);
        if fermi_dirac_half(mid)? < value {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

pub fn fermi_electron_density(
    conduction_band_density_m3: f64,
    reduced_fermi_to_conduction: f64,
) -> Result<f64> {
    require(
        conduction_band_density_m3 > 0.0,
        "conduction-band density must be positive",
    )?;
    Ok(conduction_band_density_m3 * fermi_dirac_half(reduced_fermi_to_conduction)?)
}

pub fn fermi_hole_density(
    valence_band_density_m3: f64,
    reduced_fermi_from_valence: f64,
) -> Result<f64> {
    require(
        valence_band_density_m3 > 0.0,
        "valence-band density must be positive",
    )?;
    Ok(valence_band_density_m3 * fermi_dirac_half(reduced_fermi_from_valence)?)
}

/// `reduced_fermi_to_conduction` is `(E_F - E_C)/kT`.
pub fn ionized_donor_density(
    total_donor_m3: f64,
    reduced_fermi_to_conduction: f64,
    donor_binding_ev: f64,
    degeneracy: f64,
    temperature_k: f64,
) -> Result<f64> {
    require(
        total_donor_m3 >= 0.0 && donor_binding_ev >= 0.0 && degeneracy > 0.0,
        "invalid donor-ionization parameters",
    )?;
    let exponent = reduced_fermi_to_conduction + donor_binding_ev / thermal_voltage(temperature_k)?;
    Ok(total_donor_m3 / (1.0 + degeneracy * safe_exp(exponent)))
}

/// `reduced_fermi_from_valence` is `(E_V - E_F)/kT`.
pub fn ionized_acceptor_density(
    total_acceptor_m3: f64,
    reduced_fermi_from_valence: f64,
    acceptor_binding_ev: f64,
    degeneracy: f64,
    temperature_k: f64,
) -> Result<f64> {
    require(
        total_acceptor_m3 >= 0.0 && acceptor_binding_ev >= 0.0 && degeneracy > 0.0,
        "invalid acceptor-ionization parameters",
    )?;
    let exponent =
        reduced_fermi_from_valence + acceptor_binding_ev / thermal_voltage(temperature_k)?;
    Ok(total_acceptor_m3 / (1.0 + degeneracy * safe_exp(exponent)))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MosCapacitorConfig {
    pub semiconductor_thickness_m: f64,
    pub semiconductor_nodes: usize,
    pub area_m2: f64,
    pub substrate_acceptor_m3: f64,
    pub oxide_thickness_m: f64,
    pub oxide_relative_permittivity: f64,
    pub flatband_voltage_v: f64,
    pub semiconductor: SemiconductorMaterial,
    pub max_iterations: usize,
    pub tolerance: f64,
    pub relaxation: f64,
}

impl Default for MosCapacitorConfig {
    fn default() -> Self {
        Self {
            semiconductor_thickness_m: 1e-6,
            semiconductor_nodes: 161,
            area_m2: 1e-10,
            substrate_acceptor_m3: 1e22,
            oxide_thickness_m: 10e-9,
            oxide_relative_permittivity: 3.9,
            flatband_voltage_v: 0.0,
            semiconductor: SemiconductorMaterial::default(),
            max_iterations: 1200,
            tolerance: 1e-8,
            relaxation: 0.18,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MosCapacitorPoint {
    pub gate_voltage_v: f64,
    pub surface_potential_v: f64,
    pub gate_charge_c: f64,
    pub capacitance_f: f64,
    pub converged: bool,
}

/// p-substrate MOS capacitor solved through the exact 1-D charge-surface-potential relation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MosCapacitor1D {
    config: MosCapacitorConfig,
}

impl MosCapacitor1D {
    pub fn new(config: MosCapacitorConfig) -> Result<Self> {
        let c = &config;
        require(
            c.semiconductor_thickness_m > 0.0
                && c.semiconductor_nodes >= 5
                && c.area_m2 > 0.0
                && c.substrate_acceptor_m3 > 0.0
                && c.oxide_thickness_m > 0.0
                && c.oxide_relative_permittivity > 0.0
                && c.max_iterations > 0
                && c.tolerance > 0.0
                && c.relaxation > 0.0
                && c.relaxation <= 1.0,
            "invalid MOS capacitor configuration",
        )?;
        Ok(Self { config })
    }
    pub fn oxide_capacitance_per_area(&self) -> f64 {
        VACUUM_PERMITTIVITY_F_PER_M * self.config.oxide_relative_permittivity
            / self.config.oxide_thickness_m
    }
    pub fn config(&self) -> &MosCapacitorConfig {
        &self.config
    }
    fn solve_one(&self, vg: f64) -> Result<(f64, f64)> {
        let c = &self.config;
        let ni = c.semiconductor.intrinsic_density_m3;
        let na = c.substrate_acceptor_m3;
        let vt = thermal_voltage(c.semiconductor.temperature_k)?;
        let eps = VACUUM_PERMITTIVITY_F_PER_M * c.semiconductor.relative_permittivity;
        let cox = self.oxide_capacitance_per_area();
        let root = (na * na + 4.0 * ni * ni).sqrt();
        let p0 = 0.5 * (na + root);
        let n0 = ni * ni / p0;
        let charge = |psi: f64| {
            let u = psi / vt;
            let (ep, en) = (safe_exp(-u), safe_exp(u));
            let energy = ELEMENTARY_CHARGE_C * vt * (p0 * (ep + u - 1.0) + n0 * (en - u - 1.0));
            let magnitude = (2.0 * eps * energy).max(0.0).sqrt();
            if psi.abs() < 1e-18 {
                0.0
            } else if psi > 0.0 {
                -magnitude
            } else {
                magnitude
            }
        };
        let residual = |psi: f64| c.flatband_voltage_v + psi - charge(psi) / cox - vg;
        let (mut lo, mut hi) = (-5.0, 5.0);
        if !(residual(lo) <= 0.0 && residual(hi) >= 0.0) {
            return Err(SolverError::NumericalFailure(
                "MOS capacitor surface-potential bracket failed".into(),
            ));
        }
        for _ in 0..180 {
            let mid = 0.5 * (lo + hi);
            if residual(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let psi = 0.5 * (lo + hi);
        Ok((psi, -charge(psi) * c.area_m2))
    }
    /// Upstream default perturbation 1e-4 V.
    pub fn solve(
        &self,
        gate_voltage_v: f64,
        capacitance_perturbation_v: f64,
    ) -> Result<MosCapacitorPoint> {
        let dv = capacitance_perturbation_v;
        require(
            gate_voltage_v.is_finite() && dv > 0.0 && dv.is_finite(),
            "invalid MOS capacitor bias/perturbation",
        )?;
        let (psi, q) = self.solve_one(gate_voltage_v)?;
        let (_, qp) = self.solve_one(gate_voltage_v + dv)?;
        let (_, qm) = self.solve_one(gate_voltage_v - dv)?;
        Ok(MosCapacitorPoint {
            gate_voltage_v,
            surface_potential_v: psi,
            gate_charge_c: q,
            capacitance_f: (qp - qm) / (2.0 * dv),
            converged: true,
        })
    }
    pub fn sweep_cv(
        &self,
        start_v: f64,
        stop_v: f64,
        points: usize,
        perturbation_v: f64,
    ) -> Result<Vec<MosCapacitorPoint>> {
        require(
            points >= 2 && start_v.is_finite() && stop_v.is_finite(),
            "invalid MOS capacitor C-V sweep",
        )?;
        (0..points)
            .map(|i| {
                self.solve(
                    start_v + i as f64 / (points - 1) as f64 * (stop_v - start_v),
                    perturbation_v,
                )
            })
            .collect()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongChannelMosfetConfig {
    pub width_m: f64,
    pub length_m: f64,
    pub oxide_capacitance_f_per_m2: f64,
    pub electron_mobility_m2_per_vs: f64,
    pub threshold_voltage_v: f64,
    pub channel_length_modulation_per_v: f64,
    pub body_effect_coefficient_sqrt_v: f64,
    pub surface_potential_2phi_f_v: f64,
}

impl Default for LongChannelMosfetConfig {
    fn default() -> Self {
        Self {
            width_m: 10e-6,
            length_m: 1e-6,
            oxide_capacitance_f_per_m2: 3.45e-3,
            electron_mobility_m2_per_vs: 0.05,
            threshold_voltage_v: 0.5,
            channel_length_modulation_per_v: 0.02,
            body_effect_coefficient_sqrt_v: 0.0,
            surface_potential_2phi_f_v: 0.7,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MosfetOperatingPoint {
    pub vgs_v: f64,
    pub vds_v: f64,
    pub vbs_v: f64,
    pub drain_current_a: f64,
    pub transconductance_s: f64,
    pub output_conductance_s: f64,
    pub threshold_voltage_v: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MosfetIvSample {
    pub gate_voltage_v: f64,
    pub drain_voltage_v: f64,
    pub drain_current_a: f64,
}

/// Long-channel charge-sheet nMOSFET: analytic integral of `I = W mu Qinv(V) dV/dx`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongChannelMosfet {
    config: LongChannelMosfetConfig,
}

impl LongChannelMosfet {
    pub fn new(config: LongChannelMosfetConfig) -> Result<Self> {
        let c = &config;
        require(
            c.width_m > 0.0
                && c.length_m > 0.0
                && c.oxide_capacitance_f_per_m2 > 0.0
                && c.electron_mobility_m2_per_vs > 0.0
                && c.channel_length_modulation_per_v >= 0.0
                && c.body_effect_coefficient_sqrt_v >= 0.0
                && c.surface_potential_2phi_f_v > 0.0,
            "invalid long-channel MOSFET configuration",
        )?;
        Ok(Self { config })
    }
    pub fn config(&self) -> &LongChannelMosfetConfig {
        &self.config
    }
    pub fn threshold_voltage(&self, vbs_v: f64) -> f64 {
        let c = &self.config;
        let vsb = (-vbs_v).max(0.0);
        c.threshold_voltage_v
            + c.body_effect_coefficient_sqrt_v
                * ((c.surface_potential_2phi_f_v + vsb).sqrt()
                    - c.surface_potential_2phi_f_v.sqrt())
    }
    pub fn drain_current(&self, vgs_v: f64, vds_v: f64, vbs_v: f64) -> Result<f64> {
        require(
            vgs_v.is_finite() && vds_v.is_finite() && vbs_v.is_finite(),
            "MOSFET voltages must be finite",
        )?;
        Ok(self.current(vgs_v, vds_v, vbs_v))
    }
    fn current(&self, vgs: f64, vds: f64, vbs: f64) -> f64 {
        if vds < 0.0 {
            return -self.current(vgs - vds, -vds, vbs - vds);
        }
        let c = &self.config;
        let overdrive = vgs - self.threshold_voltage(vbs);
        if overdrive.is_nan() || overdrive <= 0.0 {
            return 0.0;
        }
        let k =
            c.electron_mobility_m2_per_vs * c.oxide_capacitance_f_per_m2 * c.width_m / c.length_m;
        let veff = vds.min(overdrive);
        let mut current = k * (overdrive * veff - 0.5 * veff * veff);
        if vds > overdrive {
            current *= 1.0 + c.channel_length_modulation_per_v * (vds - overdrive);
        }
        current
    }
    pub fn operating_point(
        &self,
        vgs_v: f64,
        vds_v: f64,
        vbs_v: f64,
    ) -> Result<MosfetOperatingPoint> {
        const H: f64 = 1e-6;
        let id = self.drain_current(vgs_v, vds_v, vbs_v)?;
        let gm = (self.drain_current(vgs_v + H, vds_v, vbs_v)?
            - self.drain_current(vgs_v - H, vds_v, vbs_v)?)
            / (2.0 * H);
        let go = (self.drain_current(vgs_v, vds_v + H, vbs_v)?
            - self.drain_current(vgs_v, vds_v - H, vbs_v)?)
            / (2.0 * H);
        Ok(MosfetOperatingPoint {
            vgs_v,
            vds_v,
            vbs_v,
            drain_current_a: id,
            transconductance_s: gm,
            output_conductance_s: go,
            threshold_voltage_v: self.threshold_voltage(vbs_v),
        })
    }
    #[allow(clippy::too_many_arguments)]
    pub fn sweep(
        &self,
        gate_start_v: f64,
        gate_stop_v: f64,
        gate_points: usize,
        drain_start_v: f64,
        drain_stop_v: f64,
        drain_points: usize,
        vbs_v: f64,
    ) -> Result<Vec<MosfetIvSample>> {
        require(
            gate_points >= 2 && drain_points >= 2,
            "MOSFET sweep requires at least two points per axis",
        )?;
        let mut out = Vec::with_capacity(gate_points * drain_points);
        for ig in 0..gate_points {
            let vg =
                gate_start_v + ig as f64 / (gate_points - 1) as f64 * (gate_stop_v - gate_start_v);
            for id in 0..drain_points {
                let vd = drain_start_v
                    + id as f64 / (drain_points - 1) as f64 * (drain_stop_v - drain_start_v);
                out.push(MosfetIvSample {
                    gate_voltage_v: vg,
                    drain_voltage_v: vd,
                    drain_current_a: self.drain_current(vg, vd, vbs_v)?,
                });
            }
        }
        Ok(out)
    }
}

/// Four-terminal `[D, G, S, B]` evaluator with a centered finite-difference Jacobian (upstream
/// default step 1e-6 V).
pub fn make_long_channel_mosfet_spice_evaluator(
    config: LongChannelMosfetConfig,
    jacobian_step_v: f64,
) -> Result<StaticDeviceEvaluator> {
    require(
        jacobian_step_v > 0.0 && jacobian_step_v.is_finite(),
        "SPICE MOSFET Jacobian step must be positive",
    )?;
    let model = LongChannelMosfet::new(config)?;
    let h = jacobian_step_v;
    Ok(Arc::new(move |v: &[f64]| {
        require(
            v.len() == 4,
            "TCAD MOSFET SPICE evaluator expects [D,G,S,B]",
        )?;
        let currents = |x: &[f64]| -> Result<Vec<f64>> {
            let id = model.drain_current(x[1] - x[2], x[0] - x[2], x[3] - x[2])?;
            Ok(vec![id, 0.0, -id, 0.0])
        };
        fd_evaluation(v, 4, &currents, &|_| h)
    }))
}
