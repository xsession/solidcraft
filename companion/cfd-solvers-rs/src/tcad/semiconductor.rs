//! Scharfetter-Gummel 1-D drift-diffusion device with Gummel iteration, ohmic/Schottky/gate/
//! insulating contacts, transient, small-signal and C-V analysis.
#![allow(clippy::needless_range_loop)]
use crate::{require, Complex64 as C, Result, SolverError};

pub const ELEMENTARY_CHARGE_C: f64 = 1.602176634e-19;
pub const BOLTZMANN_J_PER_K: f64 = 1.380649e-23;
pub const VACUUM_PERMITTIVITY_F_PER_M: f64 = 8.8541878128e-12;

pub fn thermal_voltage(temperature_k: f64) -> Result<f64> {
    require(temperature_k > 0.0, "temperature must be positive")?;
    Ok(BOLTZMANN_J_PER_K * temperature_k / ELEMENTARY_CHARGE_C)
}

/// `B(x) = x / (e^x - 1)` with series and asymptotic branches.
pub fn bernoulli_function(x: f64) -> Result<f64> {
    require(x.is_finite(), "Bernoulli argument must be finite")?;
    Ok(bernoulli(x))
}

pub(crate) fn bernoulli(x: f64) -> f64 {
    if x.abs() < 1e-5 {
        let x2 = x * x;
        return 1.0 - 0.5 * x + x2 / 12.0 - x2 * x2 / 720.0;
    }
    if x > 50.0 {
        return x * (-x).exp();
    }
    if x < -50.0 {
        return -x;
    }
    x / x.exp_m1()
}

/// Varshni band gap; upstream defaults are `eg0 = 1.17`, `alpha = 4.73e-4`, `beta = 636`.
pub fn silicon_bandgap_ev(
    temperature_k: f64,
    eg0_ev: f64,
    alpha_ev_per_k: f64,
    beta_k: f64,
) -> Result<f64> {
    require(
        temperature_k > 0.0 && eg0_ev > 0.0 && alpha_ev_per_k >= 0.0 && beta_k > 0.0,
        "invalid Varshni parameters",
    )?;
    Ok(eg0_ev - alpha_ev_per_k * temperature_k * temperature_k / (temperature_k + beta_k))
}

pub(crate) fn silicon_bandgap_default(t: f64) -> Result<f64> {
    silicon_bandgap_ev(t, 1.17, 4.73e-4, 636.0)
}

pub fn temperature_scaled_intrinsic_density(
    ni_reference_m3: f64,
    reference_temperature_k: f64,
    temperature_k: f64,
    bandgap_reference_ev: f64,
    bandgap_ev: f64,
) -> Result<f64> {
    require(
        ni_reference_m3 > 0.0
            && reference_temperature_k > 0.0
            && temperature_k > 0.0
            && bandgap_reference_ev > 0.0
            && bandgap_ev > 0.0,
        "invalid intrinsic-density temperature scaling parameters",
    )?;
    let vt_ref = thermal_voltage(reference_temperature_k)?;
    let vt = thermal_voltage(temperature_k)?;
    let scale = (temperature_k / reference_temperature_k).powf(1.5);
    let exponent = bandgap_reference_ev / (2.0 * vt_ref) - bandgap_ev / (2.0 * vt);
    Ok(ni_reference_m3 * scale * exponent.clamp(-100.0, 100.0).exp())
}

pub fn temperature_scaled_mobility(
    mobility_reference_m2_per_vs: f64,
    reference_temperature_k: f64,
    temperature_k: f64,
    exponent: f64,
) -> Result<f64> {
    require(
        mobility_reference_m2_per_vs > 0.0
            && reference_temperature_k > 0.0
            && temperature_k > 0.0
            && exponent >= 0.0,
        "invalid mobility temperature scaling parameters",
    )?;
    Ok(mobility_reference_m2_per_vs * (temperature_k / reference_temperature_k).powf(-exponent))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SemiconductorMaterial {
    pub relative_permittivity: f64,
    pub intrinsic_density_m3: f64,
    pub electron_mobility_m2_per_vs: f64,
    pub hole_mobility_m2_per_vs: f64,
    pub conduction_band_density_m3: f64,
    pub valence_band_density_m3: f64,
    pub temperature_k: f64,
}

impl Default for SemiconductorMaterial {
    fn default() -> Self {
        Self {
            relative_permittivity: 11.7,
            intrinsic_density_m3: 1e16,
            electron_mobility_m2_per_vs: 0.135,
            hole_mobility_m2_per_vs: 0.048,
            conduction_band_density_m3: 2.8e25,
            valence_band_density_m3: 1.04e25,
            temperature_k: 300.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureDependence {
    pub enabled: bool,
    pub varshni_eg0_ev: f64,
    pub varshni_alpha_ev_per_k: f64,
    pub varshni_beta_k: f64,
    pub electron_mobility_exponent: f64,
    pub hole_mobility_exponent: f64,
    pub density_of_states_exponent: f64,
}

impl Default for TemperatureDependence {
    fn default() -> Self {
        Self {
            enabled: true,
            varshni_eg0_ev: 1.17,
            varshni_alpha_ev_per_k: 4.73e-4,
            varshni_beta_k: 636.0,
            electron_mobility_exponent: 2.42,
            hole_mobility_exponent: 2.20,
            density_of_states_exponent: 1.5,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RecombinationModel {
    pub tau_n_s: f64,
    pub tau_p_s: f64,
    pub trap_n1_m3: f64,
    pub trap_p1_m3: f64,
    pub auger_n_m6_per_s: f64,
    pub auger_p_m6_per_s: f64,
    pub radiative_m3_per_s: f64,
}

fn check_state(n: f64, p: f64, ni: f64, what: &str) -> Result<()> {
    require(n >= 0.0 && p >= 0.0 && ni > 0.0, what)
}

fn excess(n: f64, p: f64, ni: f64) -> f64 {
    (n - ni) * p + ni * (p - ni)
}

pub fn srh_recombination(n: f64, p: f64, ni: f64, m: &RecombinationModel) -> Result<f64> {
    check_state(n, p, ni, "invalid SRH state")?;
    if !(m.tau_n_s > 0.0 && m.tau_p_s > 0.0) {
        return Ok(0.0);
    }
    let n1 = if m.trap_n1_m3 > 0.0 { m.trap_n1_m3 } else { ni };
    let p1 = if m.trap_p1_m3 > 0.0 { m.trap_p1_m3 } else { ni };
    let d = m.tau_p_s * (n + n1) + m.tau_n_s * (p + p1);
    Ok(excess(n, p, ni) / d.max(f64::MIN_POSITIVE))
}

pub fn auger_recombination(n: f64, p: f64, ni: f64, m: &RecombinationModel) -> Result<f64> {
    check_state(n, p, ni, "invalid Auger state")?;
    Ok((m.auger_n_m6_per_s * n + m.auger_p_m6_per_s * p) * excess(n, p, ni))
}

pub fn radiative_recombination(n: f64, p: f64, ni: f64, m: &RecombinationModel) -> Result<f64> {
    check_state(n, p, ni, "invalid radiative state")?;
    Ok(m.radiative_m3_per_s * excess(n, p, ni))
}

pub fn total_recombination(n: f64, p: f64, ni: f64, m: &RecombinationModel) -> Result<f64> {
    Ok(srh_recombination(n, p, ni, m)?
        + auger_recombination(n, p, ni, m)?
        + radiative_recombination(n, p, ni, m)?)
}

pub fn caughey_thomas_mobility(
    doping_abs_m3: f64,
    mu_min: f64,
    mu_max: f64,
    reference_m3: f64,
    alpha: f64,
) -> Result<f64> {
    require(
        doping_abs_m3 >= 0.0
            && mu_min > 0.0
            && mu_max >= mu_min
            && reference_m3 > 0.0
            && alpha > 0.0,
        "invalid Caughey-Thomas mobility parameters",
    )?;
    Ok(mu_min + (mu_max - mu_min) / (1.0 + (doping_abs_m3 / reference_m3).powf(alpha)))
}

/// Upstream default `beta = 2`.
pub fn high_field_mobility(
    low_field_mobility: f64,
    electric_field_v_per_m: f64,
    saturation_velocity_m_per_s: f64,
    beta: f64,
) -> Result<f64> {
    require(
        low_field_mobility > 0.0 && saturation_velocity_m_per_s > 0.0 && beta > 0.0,
        "invalid high-field mobility parameters",
    )?;
    let ratio = low_field_mobility * electric_field_v_per_m.abs() / saturation_velocity_m_per_s;
    Ok(low_field_mobility / (1.0 + ratio.powf(beta)).powf(1.0 / beta))
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ContactType {
    #[default]
    Ohmic,
    Schottky,
    Insulating,
    Gate,
}

/// For gate contacts `oxide_capacitance_f_per_m2` drives the electrostatic Robin boundary;
/// carrier flux is zero through gate and insulating contacts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Contact1D {
    pub contact_type: ContactType,
    pub voltage_v: f64,
    pub electron_barrier_ev: f64,
    pub oxide_capacitance_f_per_m2: f64,
    pub flatband_voltage_v: f64,
}

impl Default for Contact1D {
    fn default() -> Self {
        Self {
            contact_type: ContactType::Ohmic,
            voltage_v: 0.0,
            electron_barrier_ev: 0.7,
            oxide_capacitance_f_per_m2: 0.0,
            flatband_voltage_v: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Semiconductor1DConfig {
    pub length_m: f64,
    pub nodes: usize,
    pub area_m2: f64,
    pub material: SemiconductorMaterial,
    pub temperature_dependence: TemperatureDependence,
    pub recombination: RecombinationModel,
    pub left: Contact1D,
    pub right: Contact1D,
    pub max_gummel_iterations: usize,
    pub relative_tolerance: f64,
    pub under_relaxation: f64,
    pub carrier_floor_m3: f64,
}

impl Default for Semiconductor1DConfig {
    fn default() -> Self {
        Self {
            length_m: 1e-6,
            nodes: 101,
            area_m2: 1e-12,
            material: SemiconductorMaterial::default(),
            temperature_dependence: TemperatureDependence::default(),
            recombination: RecombinationModel::default(),
            left: Contact1D::default(),
            right: Contact1D::default(),
            max_gummel_iterations: 200,
            relative_tolerance: 1e-8,
            under_relaxation: 0.35,
            carrier_floor_m3: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DcResult {
    pub converged: bool,
    pub iterations: usize,
    pub relative_change: f64,
    pub terminal_current_a: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IvPoint {
    pub voltage_v: f64,
    pub current_a: f64,
    pub converged: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransientStepResult {
    pub converged: bool,
    pub iterations: usize,
    pub relative_change: f64,
    pub time_s: f64,
    pub conduction_current_a: f64,
    pub displacement_current_a: f64,
    pub terminal_current_a: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SmallSignalResult {
    pub bias_voltage_v: f64,
    pub frequency_hz: f64,
    pub conductance_s: f64,
    pub capacitance_f: f64,
    pub admittance_s: C,
    pub impedance_ohm: C,
    pub converged: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CvPoint {
    pub voltage_v: f64,
    pub capacitance_f: f64,
    pub electrode_charge_c: f64,
    pub converged: bool,
}

fn solve_tridiagonal(
    lower: &[f64],
    mut diagonal: Vec<f64>,
    upper: &[f64],
    mut rhs: Vec<f64>,
) -> Result<Vec<f64>> {
    let n = diagonal.len();
    require(
        n > 0 && lower.len() == n && upper.len() == n && rhs.len() == n,
        "invalid tridiagonal system",
    )?;
    let singular = || SolverError::Singular("singular tridiagonal system".into());
    for i in 1..n {
        if diagonal[i - 1].is_nan() || diagonal[i - 1].abs() <= 1e-300 {
            return Err(singular());
        }
        let factor = lower[i] / diagonal[i - 1];
        diagonal[i] -= factor * upper[i - 1];
        rhs[i] -= factor * rhs[i - 1];
    }
    let mut x = vec![0.0; n];
    if diagonal[n - 1].is_nan() || diagonal[n - 1].abs() <= 1e-300 {
        return Err(singular());
    }
    x[n - 1] = rhs[n - 1] / diagonal[n - 1];
    for ii in (0..n - 1).rev() {
        if diagonal[ii].is_nan() || diagonal[ii].abs() <= 1e-300 {
            return Err(singular());
        }
        x[ii] = (rhs[ii] - upper[ii] * x[ii + 1]) / diagonal[ii];
    }
    Ok(x)
}

fn neutral_carriers(net: f64, ni: f64) -> (f64, f64) {
    let root = (net * net + 4.0 * ni * ni).sqrt();
    let n = 0.5 * (net + root);
    (n, ni * ni / n.max(f64::MIN_POSITIVE))
}

pub(crate) fn safe_exp(x: f64) -> f64 {
    x.clamp(-100.0, 100.0).exp()
}

fn validate_material(m: &SemiconductorMaterial) -> Result<()> {
    require(
        m.relative_permittivity > 0.0
            && m.intrinsic_density_m3 > 0.0
            && m.electron_mobility_m2_per_vs > 0.0
            && m.hole_mobility_m2_per_vs > 0.0
            && m.conduction_band_density_m3 > 0.0
            && m.valence_band_density_m3 > 0.0
            && m.temperature_k > 0.0,
        "invalid semiconductor material",
    )
}

fn max_relative_change(old: &[f64], new: &[f64], floor: f64) -> f64 {
    let mut result: f64 = 0.0;
    for (o, n) in old.iter().zip(new) {
        let scale = o.abs().max(n.abs()).max(floor);
        result = result.max((n - o).abs() / scale);
    }
    result
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Carrier {
    Electron,
    Hole,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SemiconductorDevice1D {
    config: Semiconductor1DConfig,
    reference_material: SemiconductorMaterial,
    x: Vec<f64>,
    potential: Vec<f64>,
    electron: Vec<f64>,
    hole: Vec<f64>,
    net_doping: Vec<f64>,
    time_s: f64,
}

impl SemiconductorDevice1D {
    pub fn new(config: Semiconductor1DConfig) -> Result<Self> {
        validate_material(&config.material)?;
        let c = &config;
        require(
            c.length_m > 0.0
                && c.nodes >= 3
                && c.nodes <= 1 << 24
                && c.area_m2 > 0.0
                && c.max_gummel_iterations > 0
                && c.relative_tolerance > 0.0
                && c.under_relaxation > 0.0
                && c.under_relaxation <= 1.0
                && c.carrier_floor_m3 > 0.0,
            "invalid semiconductor device configuration",
        )?;
        let n = c.nodes;
        let dx = c.length_m / (n - 1) as f64;
        Ok(Self {
            config,
            reference_material: config.material,
            x: (0..n).map(|i| i as f64 * dx).collect(),
            potential: vec![0.0; n],
            electron: vec![config.material.intrinsic_density_m3; n],
            hole: vec![config.material.intrinsic_density_m3; n],
            net_doping: vec![0.0; n],
            time_s: 0.0,
        })
    }
    pub fn set_net_doping(&mut self, doping: &[f64]) -> Result<()> {
        require(doping.len() == self.config.nodes, "doping size mismatch")?;
        require(
            doping.iter().all(|d| d.is_finite()),
            "doping must be finite",
        )?;
        self.net_doping = doping.to_vec();
        Ok(())
    }
    pub fn set_net_doping_fn(&mut self, doping: impl Fn(f64) -> f64) -> Result<()> {
        let values: Vec<f64> = self.x.iter().map(|x| doping(*x)).collect();
        require(
            values.iter().all(|d| d.is_finite()),
            "doping callback produced non-finite value",
        )?;
        self.net_doping = values;
        Ok(())
    }
    pub fn set_temperature(&mut self, temperature_k: f64) -> Result<()> {
        require(temperature_k > 0.0, "temperature must be positive")?;
        let td = self.config.temperature_dependence;
        let r = self.reference_material;
        let mut m = self.config.material;
        m.temperature_k = temperature_k;
        if td.enabled {
            let tref = r.temperature_k;
            let eg_ref = silicon_bandgap_ev(
                tref,
                td.varshni_eg0_ev,
                td.varshni_alpha_ev_per_k,
                td.varshni_beta_k,
            )?;
            let eg = silicon_bandgap_ev(
                temperature_k,
                td.varshni_eg0_ev,
                td.varshni_alpha_ev_per_k,
                td.varshni_beta_k,
            )?;
            m.intrinsic_density_m3 = temperature_scaled_intrinsic_density(
                r.intrinsic_density_m3,
                tref,
                temperature_k,
                eg_ref,
                eg,
            )?;
            m.electron_mobility_m2_per_vs = temperature_scaled_mobility(
                r.electron_mobility_m2_per_vs,
                tref,
                temperature_k,
                td.electron_mobility_exponent,
            )?;
            m.hole_mobility_m2_per_vs = temperature_scaled_mobility(
                r.hole_mobility_m2_per_vs,
                tref,
                temperature_k,
                td.hole_mobility_exponent,
            )?;
            let dos = (temperature_k / tref).powf(td.density_of_states_exponent);
            m.conduction_band_density_m3 = r.conduction_band_density_m3 * dos;
            m.valence_band_density_m3 = r.valence_band_density_m3 * dos;
        }
        self.config.material = m;
        Ok(())
    }
    pub fn set_contact_voltages(&mut self, left_v: f64, right_v: f64) -> Result<()> {
        require(
            left_v.is_finite() && right_v.is_finite(),
            "contact voltage must be finite",
        )?;
        self.config.left.voltage_v = left_v;
        self.config.right.voltage_v = right_v;
        Ok(())
    }
    fn vt(&self) -> f64 {
        BOLTZMANN_J_PER_K * self.config.material.temperature_k / ELEMENTARY_CHARGE_C
    }
    fn dx(&self) -> f64 {
        self.config.length_m / (self.config.nodes - 1) as f64
    }
    fn blocking(c: &Contact1D) -> bool {
        matches!(c.contact_type, ContactType::Insulating | ContactType::Gate)
    }
    fn contact_electron(&self, c: &Contact1D, doping: f64) -> Result<f64> {
        let m = &self.config.material;
        match c.contact_type {
            ContactType::Ohmic => Ok(neutral_carriers(doping, m.intrinsic_density_m3)
                .0
                .max(self.config.carrier_floor_m3)),
            ContactType::Schottky => Ok((m.conduction_band_density_m3
                * safe_exp(-c.electron_barrier_ev / self.vt()))
            .max(self.config.carrier_floor_m3)),
            _ => Err(SolverError::InvalidInput(
                "blocking contact has no prescribed electron density".into(),
            )),
        }
    }
    fn contact_hole(&self, c: &Contact1D, doping: f64) -> Result<f64> {
        let m = &self.config.material;
        match c.contact_type {
            ContactType::Ohmic => Ok(neutral_carriers(doping, m.intrinsic_density_m3)
                .1
                .max(self.config.carrier_floor_m3)),
            ContactType::Schottky => {
                let n = self.contact_electron(c, doping)?;
                Ok((m.intrinsic_density_m3 * m.intrinsic_density_m3 / n)
                    .max(self.config.carrier_floor_m3))
            }
            _ => Err(SolverError::InvalidInput(
                "blocking contact has no prescribed hole density".into(),
            )),
        }
    }
    fn contact_potential(&self, c: &Contact1D, doping: f64) -> Result<f64> {
        match c.contact_type {
            ContactType::Ohmic => Ok(self.vt()
                * (self.contact_electron(c, doping)? / self.config.material.intrinsic_density_m3)
                    .ln()
                + c.voltage_v),
            ContactType::Schottky => Ok(c.voltage_v),
            _ => Err(SolverError::InvalidInput(
                "non-Dirichlet contact has no prescribed electrostatic potential".into(),
            )),
        }
    }
    fn enforce_contact_carriers(&mut self) -> Result<()> {
        let n = self.config.nodes;
        let (l, r) = (self.config.left, self.config.right);
        if !Self::blocking(&l) {
            self.electron[0] = self.contact_electron(&l, self.net_doping[0])?;
            self.hole[0] = self.contact_hole(&l, self.net_doping[0])?;
        }
        if !Self::blocking(&r) {
            self.electron[n - 1] = self.contact_electron(&r, self.net_doping[n - 1])?;
            self.hole[n - 1] = self.contact_hole(&r, self.net_doping[n - 1])?;
        }
        Ok(())
    }
    fn apply_blocking_boundaries(&mut self) {
        let vt = self.vt();
        let floor = self.config.carrier_floor_m3;
        let n = self.config.nodes;
        if Self::blocking(&self.config.left) {
            let a = (self.potential[1] - self.potential[0]) / vt;
            self.electron[0] = (self.electron[1] * bernoulli(a) / bernoulli(-a)).max(floor);
            self.hole[0] = (self.hole[1] * bernoulli(-a) / bernoulli(a)).max(floor);
        }
        if Self::blocking(&self.config.right) {
            let a = (self.potential[n - 1] - self.potential[n - 2]) / vt;
            self.electron[n - 1] = (self.electron[n - 2] * bernoulli(-a) / bernoulli(a)).max(floor);
            self.hole[n - 1] = (self.hole[n - 2] * bernoulli(a) / bernoulli(-a)).max(floor);
        }
    }
    fn check_finite_potential(&self) -> Result<()> {
        require(
            self.potential.iter().all(|v| v.is_finite()),
            "Bernoulli argument must be finite",
        )
    }
    pub fn initialize_charge_neutral(&mut self) -> Result<()> {
        let ni = self.config.material.intrinsic_density_m3;
        let vt = self.vt();
        let floor = self.config.carrier_floor_m3;
        for i in 0..self.x.len() {
            let (n, p) = neutral_carriers(self.net_doping[i], ni);
            self.electron[i] = n.max(floor);
            self.hole[i] = p.max(floor);
            self.potential[i] = vt * (self.electron[i] / ni).ln();
        }
        let last = self.config.nodes - 1;
        let (l, r) = (self.config.left, self.config.right);
        if !Self::blocking(&l) {
            self.potential[0] = self.contact_potential(&l, self.net_doping[0])?;
        }
        if !Self::blocking(&r) {
            self.potential[last] = self.contact_potential(&r, self.net_doping[last])?;
        }
        self.enforce_contact_carriers()?;
        self.check_finite_potential()?;
        self.apply_blocking_boundaries();
        Ok(())
    }
    fn solve_poisson(&mut self) -> Result<()> {
        let n = self.config.nodes;
        let dx = self.dx();
        let eps = VACUUM_PERMITTIVITY_F_PER_M * self.config.material.relative_permittivity;
        let (l, r) = (self.config.left, self.config.right);
        require(
            !(l.contact_type == ContactType::Insulating
                && r.contact_type == ContactType::Insulating),
            "Poisson problem with two insulating boundaries has no voltage reference",
        )?;
        let (mut lo, mut diag, mut hi, mut rhs) =
            (vec![-1.0; n], vec![2.0; n], vec![-1.0; n], vec![0.0; n]);
        for i in 1..n - 1 {
            let rho = ELEMENTARY_CHARGE_C * (self.net_doping[i] + self.hole[i] - self.electron[i]);
            rhs[i] = rho * dx * dx / eps;
        }
        lo[0] = 0.0;
        match l.contact_type {
            ContactType::Ohmic | ContactType::Schottky => {
                diag[0] = 1.0;
                hi[0] = 0.0;
                rhs[0] = self.contact_potential(&l, self.net_doping[0])?;
            }
            ContactType::Insulating => {
                diag[0] = 1.0;
                hi[0] = -1.0;
                rhs[0] = 0.0;
            }
            ContactType::Gate => {
                require(
                    l.oxide_capacitance_f_per_m2 > 0.0,
                    "gate contact requires positive oxide capacitance",
                )?;
                let gamma = l.oxide_capacitance_f_per_m2 * dx / eps;
                diag[0] = 1.0 + gamma;
                hi[0] = -1.0;
                rhs[0] = gamma * (l.voltage_v - l.flatband_voltage_v);
            }
        }
        hi[n - 1] = 0.0;
        match r.contact_type {
            ContactType::Ohmic | ContactType::Schottky => {
                lo[n - 1] = 0.0;
                diag[n - 1] = 1.0;
                rhs[n - 1] = self.contact_potential(&r, self.net_doping[n - 1])?;
            }
            ContactType::Insulating => {
                lo[n - 1] = -1.0;
                diag[n - 1] = 1.0;
                rhs[n - 1] = 0.0;
            }
            ContactType::Gate => {
                require(
                    r.oxide_capacitance_f_per_m2 > 0.0,
                    "gate contact requires positive oxide capacitance",
                )?;
                let gamma = r.oxide_capacitance_f_per_m2 * dx / eps;
                lo[n - 1] = -1.0;
                diag[n - 1] = 1.0 + gamma;
                rhs[n - 1] = gamma * (r.voltage_v - r.flatband_voltage_v);
            }
        }
        self.potential = solve_tridiagonal(&lo, diag, &hi, rhs)?;
        self.check_finite_potential()
    }
    /// Scharfetter-Gummel continuity solve for one carrier; `transient = Some((old, dt))`.
    fn solve_continuity(
        &mut self,
        carrier: Carrier,
        transient: Option<(&[f64], f64)>,
    ) -> Result<()> {
        let n = self.config.nodes;
        let m = n - 2;
        let vt = self.vt();
        let dx = self.dx();
        let mobility = match carrier {
            Carrier::Electron => self.config.material.electron_mobility_m2_per_vs,
            Carrier::Hole => self.config.material.hole_mobility_m2_per_vs,
        };
        let diffusion = mobility * vt;
        let accumulation = transient.map_or(0.0, |(_, dt)| dx * dx / (diffusion * dt));
        let scale = dx * dx / diffusion;
        let s = if carrier == Carrier::Electron {
            1.0
        } else {
            -1.0
        };
        let (mut lo, mut diag, mut hi, mut rhs) =
            (vec![0.0; m], vec![0.0; m], vec![0.0; m], vec![0.0; m]);
        let ni = self.config.material.intrinsic_density_m3;
        for k in 0..m {
            let i = k + 1;
            let am = s * (self.potential[i] - self.potential[i - 1]) / vt;
            let ap = s * (self.potential[i + 1] - self.potential[i]) / vt;
            lo[k] = -bernoulli(-am);
            diag[k] = bernoulli(am) + bernoulli(-ap);
            if transient.is_some() {
                diag[k] += accumulation;
            }
            hi[k] = -bernoulli(ap);
            let recomb = total_recombination(
                self.electron[i],
                self.hole[i],
                ni,
                &self.config.recombination,
            )?;
            rhs[k] = match transient {
                Some((old, _)) => accumulation * old[i] - recomb * scale,
                None => -recomb * scale,
            };
        }
        let (l, r) = (self.config.left, self.config.right);
        if Self::blocking(&l) {
            let a = s * (self.potential[1] - self.potential[0]) / vt;
            diag[0] += lo[0] * bernoulli(a) / bernoulli(-a);
        } else {
            let edge = if carrier == Carrier::Electron {
                self.contact_electron(&l, self.net_doping[0])?
            } else {
                self.contact_hole(&l, self.net_doping[0])?
            };
            rhs[0] -= lo[0] * edge;
        }
        lo[0] = 0.0;
        if Self::blocking(&r) {
            let a = s * (self.potential[n - 1] - self.potential[n - 2]) / vt;
            diag[m - 1] += hi[m - 1] * bernoulli(-a) / bernoulli(a);
        } else {
            let edge = if carrier == Carrier::Electron {
                self.contact_electron(&r, self.net_doping[n - 1])?
            } else {
                self.contact_hole(&r, self.net_doping[n - 1])?
            };
            rhs[m - 1] -= hi[m - 1] * edge;
        }
        hi[m - 1] = 0.0;
        let interior = solve_tridiagonal(&lo, diag, &hi, rhs)?;
        let floor = self.config.carrier_floor_m3;
        let target = if carrier == Carrier::Electron {
            &mut self.electron
        } else {
            &mut self.hole
        };
        for k in 0..m {
            target[k + 1] = interior[k].max(floor);
        }
        self.enforce_contact_carriers()?;
        self.apply_blocking_boundaries();
        Ok(())
    }
    fn change(&self, op: &[f64], on: &[f64], oh: &[f64]) -> f64 {
        let f = self.config.carrier_floor_m3;
        max_relative_change(op, &self.potential, f)
            .max(max_relative_change(on, &self.electron, f))
            .max(max_relative_change(oh, &self.hole, f))
    }
    /// Gummel iteration on the equilibrium (Boltzmann) carrier relation. Errors leave the device
    /// unchanged.
    pub fn solve_equilibrium(&mut self) -> Result<DcResult> {
        let mut t = self.clone();
        let r = t.equilibrium_inner()?;
        *self = t;
        Ok(r)
    }
    fn equilibrium_inner(&mut self) -> Result<DcResult> {
        self.initialize_charge_neutral()?;
        let vt = self.vt();
        let ni = self.config.material.intrinsic_density_m3;
        let (floor, w) = (self.config.carrier_floor_m3, self.config.under_relaxation);
        let mut change = f64::INFINITY;
        for it in 0..self.config.max_gummel_iterations {
            let (op, on, oh) = (
                self.potential.clone(),
                self.electron.clone(),
                self.hole.clone(),
            );
            self.solve_poisson()?;
            for i in 0..self.x.len() {
                let nt = (ni * safe_exp(self.potential[i] / vt)).max(floor);
                let pt = (ni * safe_exp(-self.potential[i] / vt)).max(floor);
                self.electron[i] = (1.0 - w) * on[i] + w * nt;
                self.hole[i] = (1.0 - w) * oh[i] + w * pt;
            }
            self.enforce_contact_carriers()?;
            self.apply_blocking_boundaries();
            change = self.change(&op, &on, &oh);
            if change < self.config.relative_tolerance {
                return Ok(DcResult {
                    converged: true,
                    iterations: it + 1,
                    relative_change: change,
                    terminal_current_a: self.terminal_current_a(),
                });
            }
        }
        Ok(DcResult {
            converged: false,
            iterations: self.config.max_gummel_iterations,
            relative_change: change,
            terminal_current_a: self.terminal_current_a(),
        })
    }
    /// Decoupled Gummel DC solve. Errors leave the device unchanged.
    pub fn solve_dc(&mut self) -> Result<DcResult> {
        let mut t = self.clone();
        let r = t.dc_inner()?;
        *self = t;
        Ok(r)
    }
    fn dc_inner(&mut self) -> Result<DcResult> {
        self.enforce_contact_carriers()?;
        let w = self.config.under_relaxation;
        let mut change = f64::INFINITY;
        for it in 0..self.config.max_gummel_iterations {
            let (op, on, oh) = (
                self.potential.clone(),
                self.electron.clone(),
                self.hole.clone(),
            );
            self.solve_poisson()?;
            self.solve_continuity(Carrier::Electron, None)?;
            self.solve_continuity(Carrier::Hole, None)?;
            for i in 1..self.x.len() - 1 {
                self.electron[i] = (1.0 - w) * on[i] + w * self.electron[i];
                self.hole[i] = (1.0 - w) * oh[i] + w * self.hole[i];
            }
            self.enforce_contact_carriers()?;
            self.apply_blocking_boundaries();
            change = self.change(&op, &on, &oh);
            if change < self.config.relative_tolerance {
                return Ok(DcResult {
                    converged: true,
                    iterations: it + 1,
                    relative_change: change,
                    terminal_current_a: self.terminal_current_a(),
                });
            }
        }
        Ok(DcResult {
            converged: false,
            iterations: self.config.max_gummel_iterations,
            relative_change: change,
            terminal_current_a: self.terminal_current_a(),
        })
    }
    pub fn sweep_right_contact(
        &mut self,
        start_v: f64,
        stop_v: f64,
        points: usize,
    ) -> Result<Vec<IvPoint>> {
        require(
            points >= 2 && start_v.is_finite() && stop_v.is_finite(),
            "invalid IV sweep",
        )?;
        let left = self.config.left.voltage_v;
        let mut out = Vec::with_capacity(points);
        for k in 0..points {
            let t = k as f64 / (points - 1) as f64;
            self.config.right.voltage_v = start_v + t * (stop_v - start_v);
            let r = self.solve_dc()?;
            out.push(IvPoint {
                voltage_v: self.config.right.voltage_v - left,
                current_a: r.terminal_current_a,
                converged: r.converged,
            });
        }
        Ok(out)
    }
    /// Backward-Euler Gummel step; a non-converged step restores the previous state.
    pub fn step_transient(&mut self, dt_s: f64) -> Result<TransientStepResult> {
        require(
            dt_s > 0.0 && dt_s.is_finite(),
            "transient timestep must be positive and finite",
        )?;
        let mut t = self.clone();
        let r = t.transient_inner(dt_s)?;
        if r.converged {
            *self = t;
        }
        Ok(r)
    }
    fn transient_inner(&mut self, dt: f64) -> Result<TransientStepResult> {
        let (old_e, old_h) = (self.electron.clone(), self.hole.clone());
        let old_field = self.electric_field_faces();
        self.enforce_contact_carriers()?;
        let w = self.config.under_relaxation;
        let mut change = f64::INFINITY;
        let mut iterations = self.config.max_gummel_iterations;
        let mut converged = false;
        for it in 0..self.config.max_gummel_iterations {
            let (ip, ie, ih) = (
                self.potential.clone(),
                self.electron.clone(),
                self.hole.clone(),
            );
            self.solve_poisson()?;
            self.solve_continuity(Carrier::Electron, Some((&old_e, dt)))?;
            self.solve_continuity(Carrier::Hole, Some((&old_h, dt)))?;
            for i in 1..self.x.len() - 1 {
                self.electron[i] = (1.0 - w) * ie[i] + w * self.electron[i];
                self.hole[i] = (1.0 - w) * ih[i] + w * self.hole[i];
            }
            self.enforce_contact_carriers()?;
            self.apply_blocking_boundaries();
            change = self.change(&ip, &ie, &ih);
            if change < self.config.relative_tolerance {
                iterations = it + 1;
                converged = true;
                break;
            }
        }
        if !converged {
            return Ok(TransientStepResult {
                converged: false,
                iterations,
                relative_change: change,
                time_s: self.time_s,
                ..Default::default()
            });
        }
        let new_field = self.electric_field_faces();
        let conduction = self.total_current_density_faces()[0] * self.config.area_m2;
        let eps = VACUUM_PERMITTIVITY_F_PER_M * self.config.material.relative_permittivity;
        let displacement = eps * self.config.area_m2 * (new_field[0] - old_field[0]) / dt;
        self.time_s += dt;
        Ok(TransientStepResult {
            converged: true,
            iterations,
            relative_change: change,
            time_s: self.time_s,
            conduction_current_a: conduction,
            displacement_current_a: displacement,
            terminal_current_a: conduction + displacement,
        })
    }
    /// Quasi-static admittance from symmetric left-contact perturbations (default 1e-4 V).
    pub fn small_signal(
        &self,
        frequency_hz: f64,
        perturbation_v: f64,
    ) -> Result<SmallSignalResult> {
        require(
            frequency_hz >= 0.0
                && frequency_hz.is_finite()
                && perturbation_v > 0.0
                && perturbation_v.is_finite(),
            "invalid small-signal parameters",
        )?;
        let (mut plus, mut minus) = (self.clone(), self.clone());
        plus.config.left.voltage_v += perturbation_v;
        minus.config.left.voltage_v -= perturbation_v;
        let (pd, md) = (plus.solve_dc()?, minus.solve_dc()?);
        let dv = 2.0 * perturbation_v;
        let conductance = (pd.terminal_current_a - md.terminal_current_a) / dv;
        let capacitance = (plus.left_terminal_charge_c() - minus.left_terminal_charge_c()) / dv;
        let omega = 2.0 * (-1.0f64).acos() * frequency_hz;
        let admittance = C::new(conductance, omega * capacitance);
        let impedance = if admittance.norm() > 0.0 {
            1.0 / admittance
        } else {
            C::new(f64::INFINITY, 0.0)
        };
        Ok(SmallSignalResult {
            bias_voltage_v: self.config.left.voltage_v - self.config.right.voltage_v,
            frequency_hz,
            conductance_s: conductance,
            capacitance_f: capacitance,
            admittance_s: admittance,
            impedance_ohm: impedance,
            converged: pd.converged && md.converged,
        })
    }
    pub fn sweep_cv(
        &self,
        start_left_v: f64,
        stop_left_v: f64,
        points: usize,
        perturbation_v: f64,
    ) -> Result<Vec<CvPoint>> {
        require(
            points >= 2
                && start_left_v.is_finite()
                && stop_left_v.is_finite()
                && perturbation_v > 0.0
                && perturbation_v.is_finite(),
            "invalid CV sweep parameters",
        )?;
        let mut out = Vec::with_capacity(points);
        for k in 0..points {
            let t = k as f64 / (points - 1) as f64;
            let mut bias = self.clone();
            bias.config.left.voltage_v = start_left_v + t * (stop_left_v - start_left_v);
            let dc = bias.solve_dc()?;
            let ac = bias.small_signal(0.0, perturbation_v)?;
            out.push(CvPoint {
                voltage_v: bias.config.left.voltage_v - bias.config.right.voltage_v,
                capacitance_f: ac.capacitance_f,
                electrode_charge_c: bias.left_terminal_charge_c(),
                converged: dc.converged && ac.converged,
            });
        }
        Ok(out)
    }
    pub fn config(&self) -> &Semiconductor1DConfig {
        &self.config
    }
    pub fn x(&self) -> &[f64] {
        &self.x
    }
    pub fn potential(&self) -> &[f64] {
        &self.potential
    }
    pub fn electron_density(&self) -> &[f64] {
        &self.electron
    }
    pub fn hole_density(&self) -> &[f64] {
        &self.hole
    }
    pub fn net_doping(&self) -> &[f64] {
        &self.net_doping
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    fn faces(&self, carrier: Carrier) -> Vec<f64> {
        let nf = self.config.nodes - 1;
        let vt = self.vt();
        let dx = self.config.length_m / nf as f64;
        let m = &self.config.material;
        let mut out = vec![0.0; nf];
        for i in 0..nf {
            let a = (self.potential[i + 1] - self.potential[i]) / vt;
            out[i] = match carrier {
                Carrier::Electron => {
                    ELEMENTARY_CHARGE_C * m.electron_mobility_m2_per_vs * vt / dx
                        * (self.electron[i + 1] * bernoulli(a) - self.electron[i] * bernoulli(-a))
                }
                Carrier::Hole => {
                    ELEMENTARY_CHARGE_C * m.hole_mobility_m2_per_vs * vt / dx
                        * (self.hole[i] * bernoulli(a) - self.hole[i + 1] * bernoulli(-a))
                }
            };
        }
        if Self::blocking(&self.config.left) {
            out[0] = 0.0;
        }
        if Self::blocking(&self.config.right) {
            out[nf - 1] = 0.0;
        }
        out
    }
    pub fn electron_current_density_faces(&self) -> Vec<f64> {
        self.faces(Carrier::Electron)
    }
    pub fn hole_current_density_faces(&self) -> Vec<f64> {
        self.faces(Carrier::Hole)
    }
    pub fn total_current_density_faces(&self) -> Vec<f64> {
        let mut out = self.faces(Carrier::Electron);
        for (o, h) in out.iter_mut().zip(self.faces(Carrier::Hole)) {
            *o += h;
        }
        out
    }
    pub fn electric_field_faces(&self) -> Vec<f64> {
        let nf = self.config.nodes - 1;
        let dx = self.config.length_m / nf as f64;
        (0..nf)
            .map(|i| -(self.potential[i + 1] - self.potential[i]) / dx)
            .collect()
    }
    pub fn joule_heating_cells(&self) -> Vec<f64> {
        let face: Vec<f64> = self
            .total_current_density_faces()
            .iter()
            .zip(self.electric_field_faces())
            .map(|(j, e)| j * e)
            .collect();
        let n = self.config.nodes;
        let mut out = vec![0.0; n];
        out[0] = face[0];
        out[n - 1] = face[face.len() - 1];
        for i in 1..n - 1 {
            out[i] = 0.5 * (face[i - 1] + face[i]);
        }
        out
    }
    pub fn terminal_current_a(&self) -> f64 {
        let c = self.total_current_density_faces();
        0.5 * (c[0] + c[c.len() - 1]) * self.config.area_m2
    }
    pub fn left_terminal_charge_c(&self) -> f64 {
        let l = &self.config.left;
        match l.contact_type {
            ContactType::Gate => {
                l.oxide_capacitance_f_per_m2
                    * self.config.area_m2
                    * (l.voltage_v - l.flatband_voltage_v - self.potential[0])
            }
            ContactType::Insulating => 0.0,
            _ => {
                VACUUM_PERMITTIVITY_F_PER_M
                    * self.config.material.relative_permittivity
                    * self.config.area_m2
                    * self.electric_field_faces()[0]
            }
        }
    }
    pub fn total_mobile_charge_c(&self) -> f64 {
        let n = self.config.nodes;
        let mut integral = 0.0;
        for i in 0..n {
            let w = if i == 0 || i + 1 == n { 0.5 } else { 1.0 };
            integral += w * (self.hole[i] - self.electron[i]);
        }
        ELEMENTARY_CHARGE_C * self.config.area_m2 * self.dx() * integral
    }
}

pub fn pn_junction_doping(
    x: &[f64],
    junction_position_m: f64,
    acceptor_m3: f64,
    donor_m3: f64,
) -> Result<Vec<f64>> {
    require(
        !x.is_empty() && acceptor_m3 > 0.0 && donor_m3 > 0.0,
        "invalid PN doping parameters",
    )?;
    Ok(x.iter()
        .map(|xi| {
            if *xi < junction_position_m {
                -acceptor_m3
            } else {
                donor_m3
            }
        })
        .collect())
}

/// Upstream default `intrinsic_net_m3 = 0`.
pub fn pin_diode_doping(
    x: &[f64],
    p_end_m: f64,
    n_begin_m: f64,
    acceptor_m3: f64,
    donor_m3: f64,
    intrinsic_net_m3: f64,
) -> Result<Vec<f64>> {
    require(
        !x.is_empty()
            && p_end_m >= 0.0
            && n_begin_m > p_end_m
            && acceptor_m3 > 0.0
            && donor_m3 > 0.0,
        "invalid PIN doping parameters",
    )?;
    Ok(x.iter()
        .map(|xi| {
            if *xi < p_end_m {
                -acceptor_m3
            } else if *xi >= n_begin_m {
                donor_m3
            } else {
                intrinsic_net_m3
            }
        })
        .collect())
}
