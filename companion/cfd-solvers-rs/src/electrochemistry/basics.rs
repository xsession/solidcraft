use crate::{chemistry::GAS_CONSTANT, require, Result, SolverError};

/// C/mol.
pub const FARADAY_CONSTANT: f64 = 96485.33212;

/// `E = E0 - RT/(nF) ln Q`.
pub fn nernst_potential(
    standard_potential: f64,
    temperature: f64,
    electrons: f64,
    reaction_quotient: f64,
) -> Result<f64> {
    require(
        temperature > 0.0 && electrons > 0.0 && reaction_quotient > 0.0,
        "invalid Nernst parameters",
    )?;
    Ok(standard_potential
        - (GAS_CONSTANT * temperature / (electrons * FARADAY_CONSTANT)) * reaction_quotient.ln())
}

/// `i0 (exp(aa nF eta/RT) - exp(-ac nF eta/RT))`; upstream defaults n = 1, aa = ac = 0.5.
pub fn butler_volmer_current_density(
    exchange_current_density: f64,
    overpotential: f64,
    temperature: f64,
    electrons: f64,
    anodic_transfer: f64,
    cathodic_transfer: f64,
) -> Result<f64> {
    require(
        exchange_current_density >= 0.0
            && temperature > 0.0
            && electrons > 0.0
            && anodic_transfer >= 0.0
            && cathodic_transfer >= 0.0
            && anodic_transfer + cathodic_transfer > 0.0,
        "invalid Butler-Volmer parameters",
    )?;
    let scale = electrons * FARADAY_CONSTANT / (GAS_CONSTANT * temperature);
    Ok(exchange_current_density
        * ((anodic_transfer * scale * overpotential).exp()
            - (-cathodic_transfer * scale * overpotential).exp()))
}

pub(crate) fn bv(i0: f64, eta: f64, t: f64, n: f64) -> Result<f64> {
    butler_volmer_current_density(i0, eta, t, n, 0.5, 0.5)
}

pub fn faradaic_molar_flux(current_density: f64, electrons: f64) -> Result<f64> {
    require(electrons > 0.0, "electron count must be positive")?;
    Ok(current_density / (electrons * FARADAY_CONSTANT))
}

/// Metal recession speed `i M / (n F rho)` in m/s.
pub fn corrosion_penetration_rate(
    current_density: f64,
    molar_mass: f64,
    metal_density: f64,
    electrons: f64,
) -> Result<f64> {
    require(
        molar_mass > 0.0 && metal_density > 0.0,
        "invalid corrosion material parameters",
    )?;
    Ok(faradaic_molar_flux(current_density, electrons)? * molar_mass / metal_density)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitPnpElectrodeConfig {
    pub dt_s: f64,
    pub diffusivity_m2_s: f64,
    pub control_volume_m3: f64,
    pub boundary_concentration_mol_m3: f64,
    pub exchange_current_a: f64,
    pub equilibrium_potential_v: f64,
    pub double_layer_capacitance_f: f64,
    pub temperature_k: f64,
    pub electrons: f64,
}

impl Default for ImplicitPnpElectrodeConfig {
    fn default() -> Self {
        Self {
            dt_s: 1e-3,
            diffusivity_m2_s: 1e-9,
            control_volume_m3: 1e-9,
            boundary_concentration_mol_m3: 1.0,
            exchange_current_a: 1e-6,
            equilibrium_potential_v: 0.0,
            double_layer_capacitance_f: 1e-6,
            temperature_k: 298.15,
            electrons: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ImplicitPnpElectrodeState {
    pub concentration_mol_m3: f64,
    pub potential_v: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// One backward-Euler step of a lumped concentration/double-layer electrode cell solved by a
/// 2x2 Newton iteration (upstream defaults 40 iterations, tolerance 1e-10).
pub fn solve_implicit_pnp_electrode_cell(
    old_concentration_mol_m3: f64,
    old_potential_v: f64,
    config: &ImplicitPnpElectrodeConfig,
    max_iterations: usize,
    tolerance: f64,
) -> Result<ImplicitPnpElectrodeState> {
    let c0 = config;
    let (old_c, old_phi) = (old_concentration_mol_m3, old_potential_v);
    require(
        old_c >= 0.0
            && old_phi.is_finite()
            && c0.dt_s > 0.0
            && c0.diffusivity_m2_s >= 0.0
            && c0.control_volume_m3 > 0.0
            && c0.boundary_concentration_mol_m3 >= 0.0
            && c0.exchange_current_a >= 0.0
            && c0.double_layer_capacitance_f > 0.0
            && c0.temperature_k > 0.0
            && c0.electrons > 0.0
            && max_iterations > 0
            && tolerance > 0.0,
        "invalid implicit PNP/electrode controls",
    )?;
    let transport = c0.diffusivity_m2_s / c0.control_volume_m3;
    let nf = c0.electrons * FARADAY_CONSTANT;
    let (mut c, mut phi) = (old_c, old_phi);
    let mut r = ImplicitPnpElectrodeState {
        concentration_mol_m3: c,
        potential_v: phi,
        iterations: 0,
        converged: false,
    };
    for iteration in 0..max_iterations {
        let arg = nf * (phi - c0.equilibrium_potential_v) / (2.0 * GAS_CONSTANT * c0.temperature_k);
        let current = c0.exchange_current_a * (arg.exp() - (-arg).exp());
        let didphi = c0.exchange_current_a * nf / (2.0 * GAS_CONSTANT * c0.temperature_k)
            * (arg.exp() + (-arg).exp());
        let rc = c
            - old_c
            - c0.dt_s
                * (transport * (c0.boundary_concentration_mol_m3 - c)
                    - current / (nf * c0.control_volume_m3));
        let rp = c0.double_layer_capacitance_f * (phi - old_phi) + c0.dt_s * current;
        r.iterations = iteration + 1;
        if rc.abs().max(rp.abs()) <= tolerance {
            r.concentration_mol_m3 = c.max(0.0);
            r.potential_v = phi;
            r.converged = true;
            return Ok(r);
        }
        let j11 = 1.0 + c0.dt_s * transport;
        let j12 = c0.dt_s * didphi / (nf * c0.control_volume_m3);
        let j22 = c0.double_layer_capacitance_f + c0.dt_s * didphi;
        let det = j11 * j22;
        if det.abs() < 1e-30 {
            return Err(SolverError::Singular(
                "implicit PNP Jacobian is singular".into(),
            ));
        }
        let dc = (-rc * j22 + j12 * rp) / det;
        let dphi = (-j11 * rp) / det;
        c = (c + dc).max(0.0);
        phi += dphi;
    }
    r.concentration_mol_m3 = c;
    r.potential_v = phi;
    Ok(r)
}
