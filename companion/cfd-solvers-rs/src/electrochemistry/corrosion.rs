use super::{
    butler_volmer_current_density, bv, corrosion_penetration_rate, faradaic_molar_flux,
    nernst_potential, FARADAY_CONSTANT,
};
use crate::{
    chemistry::{equilibrate_acids, AcidFamily, GAS_CONSTANT},
    require, Complex64, Result, SolverError,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CorrosionCell1DConfig {
    pub electrolyte_length: f64,
    pub electrolyte_conductivity: f64,
    pub temperature: f64,
    pub metal_potential: f64,
    pub bulk_electrolyte_potential: f64,
    pub equilibrium_potential: f64,
    pub exchange_current_density: f64,
    pub electrons: f64,
    pub anodic_transfer: f64,
    pub cathodic_transfer: f64,
    pub metal_molar_mass: f64,
    pub metal_density: f64,
    pub max_iterations: usize,
    pub current_tolerance: f64,
}

impl Default for CorrosionCell1DConfig {
    fn default() -> Self {
        Self {
            electrolyte_length: 1e-3,
            electrolyte_conductivity: 5.0,
            temperature: 298.15,
            metal_potential: 0.0,
            bulk_electrolyte_potential: 0.0,
            equilibrium_potential: 0.0,
            exchange_current_density: 1.0,
            electrons: 2.0,
            anodic_transfer: 0.5,
            cathodic_transfer: 0.5,
            metal_molar_mass: 0.055845,
            metal_density: 7874.0,
            max_iterations: 200,
            current_tolerance: 1e-12,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CorrosionCell1DResult {
    pub current_density: f64,
    pub surface_electrolyte_potential: f64,
    pub overpotential: f64,
    pub molar_dissolution_flux: f64,
    pub penetration_rate: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// Butler-Volmer kinetics in series with an ohmic electrolyte layer, by bracketed bisection on
/// the current density (anodic positive).
pub fn solve_corrosion_cell_1d(c: &CorrosionCell1DConfig) -> Result<CorrosionCell1DResult> {
    require(
        c.electrolyte_length > 0.0
            && c.electrolyte_conductivity > 0.0
            && c.temperature > 0.0
            && c.exchange_current_density >= 0.0
            && c.electrons > 0.0
            && c.max_iterations > 0
            && c.current_tolerance > 0.0,
        "invalid 1-D corrosion cell configuration",
    )?;
    let surface = |i: f64| {
        c.bulk_electrolyte_potential + i * c.electrolyte_length / c.electrolyte_conductivity
    };
    let residual = |i: f64| -> Result<f64> {
        let eta = c.metal_potential - surface(i) - c.equilibrium_potential;
        Ok(i - butler_volmer_current_density(
            c.exchange_current_density,
            eta,
            c.temperature,
            c.electrons,
            c.anodic_transfer,
            c.cathodic_transfer,
        )?)
    };
    let applied = c.metal_potential - c.bulk_electrolyte_potential - c.equilibrium_potential;
    if applied.abs() <= 1e-18 || c.exchange_current_density == 0.0 {
        return Ok(CorrosionCell1DResult {
            surface_electrolyte_potential: c.bulk_electrolyte_potential,
            overpotential: applied,
            converged: true,
            ..Default::default()
        });
    }
    let scale = c
        .exchange_current_density
        .max(applied.abs() * c.electrolyte_conductivity / c.electrolyte_length)
        .max(1e-12);
    let (mut lo, mut hi) = (-scale, scale);
    let (mut flo, mut fhi) = (residual(lo)?, residual(hi)?);
    let mut expand = 0;
    while expand < 80 && !(flo <= 0.0 && fhi >= 0.0) {
        lo *= 2.0;
        hi *= 2.0;
        flo = residual(lo)?;
        fhi = residual(hi)?;
        if !flo.is_finite() || !fhi.is_finite() {
            break;
        }
        expand += 1;
    }
    if !(flo.is_finite() && fhi.is_finite() && flo <= 0.0 && fhi >= 0.0) {
        return Err(SolverError::NumericalFailure(
            "failed to bracket corrosion-cell current density".into(),
        ));
    }
    let (mut current, mut iteration, mut converged) = (0.0, 0, false);
    while iteration < c.max_iterations {
        current = 0.5 * (lo + hi);
        let f = residual(current)?;
        let tol = c.current_tolerance * current.abs().max(1.0);
        if f.abs() <= tol || 0.5 * (hi - lo).abs() <= tol {
            converged = true;
            iteration += 1;
            break;
        }
        if f > 0.0 {
            hi = current;
        } else {
            lo = current;
        }
        iteration += 1;
    }
    let s = surface(current);
    Ok(CorrosionCell1DResult {
        current_density: current,
        surface_electrolyte_potential: s,
        overpotential: c.metal_potential - s - c.equilibrium_potential,
        molar_dissolution_flux: faradaic_molar_flux(current, c.electrons)?,
        penetration_rate: corrosion_penetration_rate(
            current,
            c.metal_molar_mass,
            c.metal_density,
            c.electrons,
        )?,
        iterations: iteration,
        converged,
    })
}

#[derive(Clone, Debug, PartialEq)]
pub struct ElectrodeReaction {
    pub name: String,
    pub equilibrium_potential: f64,
    pub exchange_current_density: f64,
    pub area: f64,
    pub electrons: f64,
    pub alpha_anodic: f64,
    pub alpha_cathodic: f64,
}

impl Default for ElectrodeReaction {
    fn default() -> Self {
        Self {
            name: String::new(),
            equilibrium_potential: 0.0,
            exchange_current_density: 0.0,
            area: 1.0,
            electrons: 1.0,
            alpha_anodic: 0.5,
            alpha_cathodic: 0.5,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MixedPotentialResult {
    pub potential: f64,
    pub total_current: f64,
    pub reaction_current_density: Vec<f64>,
    pub reaction_current: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
}

/// Potential where the summed Butler-Volmer currents equal `target_current`, by bisection
/// (upstream defaults T = 298.15, target 0, 200 iterations, relative tolerance 1e-10).
pub fn solve_mixed_potential(
    reactions: &[ElectrodeReaction],
    temperature: f64,
    target_current: f64,
    max_iterations: usize,
    relative_current_tolerance: f64,
) -> Result<MixedPotentialResult> {
    require(
        !reactions.is_empty(),
        "mixed-potential solve requires reactions",
    )?;
    require(
        temperature > 0.0
            && max_iterations > 0
            && relative_current_tolerance > 0.0
            && target_current.is_finite(),
        "invalid mixed-potential controls",
    )?;
    let (mut min_eq, mut max_eq) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut scale = target_current.abs();
    for r in reactions {
        require(
            !r.name.is_empty()
                && r.exchange_current_density >= 0.0
                && r.area > 0.0
                && r.electrons > 0.0
                && r.alpha_anodic > 0.0
                && r.alpha_cathodic > 0.0
                && r.equilibrium_potential.is_finite(),
            "invalid electrode reaction",
        )?;
        min_eq = min_eq.min(r.equilibrium_potential);
        max_eq = max_eq.max(r.equilibrium_potential);
        scale += r.area * r.exchange_current_density;
    }
    let tolerance = relative_current_tolerance * scale.max(1e-30);
    let density = |r: &ElectrodeReaction, e: f64| {
        butler_volmer_current_density(
            r.exchange_current_density,
            e - r.equilibrium_potential,
            temperature,
            r.electrons,
            r.alpha_anodic,
            r.alpha_cathodic,
        )
    };
    let total = |e: f64| -> Result<f64> {
        let mut t = -target_current;
        for r in reactions {
            t += r.area * density(r, e)?;
        }
        Ok(t)
    };
    let mut margin = 0.1;
    let (mut lo, mut hi) = (min_eq - margin, max_eq + margin);
    let (mut flo, mut fhi) = (total(lo)?, total(hi)?);
    let mut expansion = 0;
    while flo * fhi > 0.0 && expansion < 12 {
        margin *= 2.0;
        lo = min_eq - margin;
        hi = max_eq + margin;
        flo = total(lo)?;
        fhi = total(hi)?;
        expansion += 1;
    }
    if !(flo <= 0.0 && fhi >= 0.0) {
        return Err(SolverError::NumericalFailure(
            "mixed-potential current balance could not be bracketed".into(),
        ));
    }
    let mut result = MixedPotentialResult::default();
    for iteration in 0..max_iterations {
        let mid = 0.5 * (lo + hi);
        let fm = total(mid)?;
        result.iterations = iteration + 1;
        result.potential = mid;
        if fm.abs() <= tolerance || (hi - lo).abs() <= 1e-14 {
            result.converged = true;
            break;
        }
        if fm > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    for r in reactions {
        let j = density(r, result.potential)?;
        result.reaction_current_density.push(j);
        result.reaction_current.push(r.area * j);
        result.total_current += r.area * j;
    }
    Ok(result)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RedoxActivityTerm {
    pub activity: f64,
    pub stoichiometric_power: f64,
}

/// `prod(a_p^nu) / prod(a_r^nu)`.
pub fn reaction_quotient(
    products: &[RedoxActivityTerm],
    reactants: &[RedoxActivityTerm],
) -> Result<f64> {
    let mut q = 1.0;
    for (terms, sign) in [(products, 1.0), (reactants, -1.0)] {
        for t in terms {
            require(
                t.activity > 0.0 && t.activity.is_finite() && t.stoichiometric_power.is_finite(),
                "invalid redox activity",
            )?;
            q *= t.activity.powf(sign * t.stoichiometric_power);
        }
    }
    require(q > 0.0 && q.is_finite(), "reaction quotient out of range")?;
    Ok(q)
}

pub fn concentration_dependent_nernst(
    standard_potential: f64,
    temperature: f64,
    electrons: f64,
    products: &[RedoxActivityTerm],
    reactants: &[RedoxActivityTerm],
) -> Result<f64> {
    nernst_potential(
        standard_potential,
        temperature,
        electrons,
        reaction_quotient(products, reactants)?,
    )
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CathodicReactionResult {
    pub equilibrium_potential: f64,
    pub activation_current_density: f64,
    pub limited_current_density: f64,
}

/// Koutecky-Levich combination of a cathodic activation current with a limiting current.
fn cathodic_limit(activation: f64, limit: f64) -> f64 {
    if !(limit > 0.0 && limit.is_finite()) {
        return activation;
    }
    let cath = (-activation).max(0.0);
    let limited = cath / (1.0 + cath / limit);
    if activation < 0.0 {
        -limited
    } else {
        activation
    }
}

/// O2 + 4H+ + 4e- -> 2H2O (upstream defaults T = 298.15, n = 4).
#[allow(clippy::too_many_arguments)]
pub fn oxygen_reduction_cathodic(
    electrode_potential: f64,
    standard_potential: f64,
    oxygen_activity: f64,
    proton_activity: f64,
    exchange_current_density: f64,
    limiting_current_density: f64,
    temperature: f64,
    electrons: f64,
) -> Result<CathodicReactionResult> {
    let react = [
        RedoxActivityTerm {
            activity: oxygen_activity,
            stoichiometric_power: 1.0,
        },
        RedoxActivityTerm {
            activity: proton_activity,
            stoichiometric_power: 4.0,
        },
    ];
    let prod = [RedoxActivityTerm {
        activity: 1.0,
        stoichiometric_power: 2.0,
    }];
    let eq =
        concentration_dependent_nernst(standard_potential, temperature, electrons, &prod, &react)?;
    let a = bv(
        exchange_current_density,
        electrode_potential - eq,
        temperature,
        electrons,
    )?;
    Ok(CathodicReactionResult {
        equilibrium_potential: eq,
        activation_current_density: a,
        limited_current_density: cathodic_limit(a, limiting_current_density),
    })
}

/// 2H+ + 2e- -> H2 against the standard hydrogen electrode (upstream defaults T = 298.15, n = 2).
pub fn hydrogen_evolution_cathodic(
    electrode_potential: f64,
    hydrogen_pressure_activity: f64,
    proton_activity: f64,
    exchange_current_density: f64,
    limiting_current_density: f64,
    temperature: f64,
    electrons: f64,
) -> Result<CathodicReactionResult> {
    let react = [RedoxActivityTerm {
        activity: proton_activity,
        stoichiometric_power: 2.0,
    }];
    let prod = [RedoxActivityTerm {
        activity: hydrogen_pressure_activity,
        stoichiometric_power: 1.0,
    }];
    let eq = concentration_dependent_nernst(0.0, temperature, electrons, &prod, &react)?;
    let a = bv(
        exchange_current_density,
        electrode_potential - eq,
        temperature,
        electrons,
    )?;
    Ok(CathodicReactionResult {
        equilibrium_potential: eq,
        activation_current_density: a,
        limited_current_density: cathodic_limit(a, limiting_current_density),
    })
}

/// Linearized `R_ct = RT / (i0 n F (aa + ac))`.
pub fn charge_transfer_resistance(
    exchange_current_density: f64,
    temperature: f64,
    electrons: f64,
    anodic_transfer: f64,
    cathodic_transfer: f64,
) -> Result<f64> {
    require(
        exchange_current_density > 0.0
            && temperature > 0.0
            && electrons > 0.0
            && anodic_transfer + cathodic_transfer > 0.0,
        "invalid charge-transfer controls",
    )?;
    Ok(GAS_CONSTANT * temperature
        / (exchange_current_density
            * electrons
            * FARADAY_CONSTANT
            * (anodic_transfer + cathodic_transfer)))
}

/// Randles cell `Rs + 1 / (1/Rct + j w Cdl + 1/Zw)` with semi-infinite Warburg `sigma/sqrt(w)(1-j)`.
pub fn randles_impedance(
    frequency_hz: f64,
    solution_resistance: f64,
    charge_transfer_resistance_ohm: f64,
    double_layer_capacitance_f: f64,
    warburg_sigma: f64,
) -> Result<Complex64> {
    require(
        frequency_hz >= 0.0
            && solution_resistance >= 0.0
            && charge_transfer_resistance_ohm > 0.0
            && double_layer_capacitance_f >= 0.0
            && warburg_sigma >= 0.0,
        "invalid EIS parameters",
    )?;
    let w = 2.0 * std::f64::consts::PI * frequency_hz;
    let mut adm = Complex64::new(
        1.0 / charge_transfer_resistance_ohm,
        w * double_layer_capacitance_f,
    );
    if warburg_sigma > 0.0 && w > 0.0 {
        let zw = warburg_sigma / w.sqrt() * Complex64::new(1.0, -1.0);
        adm += 1.0 / zw;
    }
    Ok(solution_resistance + 1.0 / adm)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DoubleLayerState {
    pub capacitance_per_area: f64,
    pub potential: f64,
}

impl DoubleLayerState {
    pub fn charge_density(&self) -> f64 {
        self.capacitance_per_area * self.potential
    }
    pub fn advance_current(&mut self, current_density: f64, dt: f64) -> Result<()> {
        require(
            self.capacitance_per_area > 0.0 && dt >= 0.0 && current_density.is_finite(),
            "invalid double-layer state",
        )?;
        self.potential += current_density * dt / self.capacitance_per_area;
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PassivationModel {
    pub active_equilibrium_potential: f64,
    pub exchange_current_density: f64,
    pub passive_onset: f64,
    pub transpassive_onset: f64,
    pub passive_current_density: f64,
    /// A/(m2 V).
    pub transpassive_slope: f64,
    pub electrons: f64,
}

impl Default for PassivationModel {
    fn default() -> Self {
        Self {
            active_equilibrium_potential: 0.0,
            exchange_current_density: 0.0,
            passive_onset: 0.0,
            transpassive_onset: 0.0,
            passive_current_density: 0.0,
            transpassive_slope: 0.0,
            electrons: 1.0,
        }
    }
}

/// Active Butler-Volmer / passive plateau / linear transpassive branches.
pub fn passivation_current_density(
    electrode_potential: f64,
    m: &PassivationModel,
    temperature: f64,
) -> Result<f64> {
    require(
        m.exchange_current_density >= 0.0
            && m.passive_current_density >= 0.0
            && m.transpassive_slope >= 0.0
            && m.electrons > 0.0
            && temperature > 0.0
            && m.transpassive_onset >= m.passive_onset,
        "invalid passivation model",
    )?;
    if electrode_potential < m.passive_onset {
        return bv(
            m.exchange_current_density,
            electrode_potential - m.active_equilibrium_potential,
            temperature,
            m.electrons,
        );
    }
    if electrode_potential <= m.transpassive_onset {
        return Ok(m.passive_current_density);
    }
    Ok(m.passive_current_density
        + m.transpassive_slope * (electrode_potential - m.transpassive_onset))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProductLayerState {
    /// m.
    pub thickness: f64,
    /// ohm m.
    pub resistivity: f64,
    pub growth_fraction: f64,
}

impl Default for ProductLayerState {
    fn default() -> Self {
        Self {
            thickness: 0.0,
            resistivity: 0.0,
            growth_fraction: 1.0,
        }
    }
}

impl ProductLayerState {
    pub fn area_specific_resistance(&self) -> f64 {
        self.thickness * self.resistivity
    }
    pub fn advance_from_anodic_current(
        &mut self,
        current_density: f64,
        dt: f64,
        molar_mass: f64,
        density: f64,
        electrons: f64,
    ) -> Result<()> {
        require(
            self.thickness >= 0.0
                && self.resistivity >= 0.0
                && (0.0..=1.0).contains(&self.growth_fraction)
                && dt >= 0.0,
            "invalid product-layer state",
        )?;
        self.thickness += self.growth_fraction
            * faradaic_recession_distance(
                current_density.max(0.0),
                dt,
                molar_mass,
                density,
                electrons,
            )?;
        Ok(())
    }
}

pub fn faradaic_recession_distance(
    current_density: f64,
    dt: f64,
    molar_mass: f64,
    metal_density: f64,
    electrons: f64,
) -> Result<f64> {
    require(dt >= 0.0 && dt.is_finite(), "invalid recession time")?;
    Ok(
        corrosion_penetration_rate(current_density.abs(), molar_mass, metal_density, electrons)?
            * dt,
    )
}

pub fn cathodic_protection_current_density(
    protected_potential: f64,
    free_corrosion_potential: f64,
    polarization_resistance: f64,
) -> Result<f64> {
    require(
        polarization_resistance > 0.0
            && protected_potential.is_finite()
            && free_corrosion_potential.is_finite(),
        "invalid cathodic-protection controls",
    )?;
    Ok((protected_potential - free_corrosion_potential) / polarization_resistance)
}

/// `bulk * porosity^exponent` (upstream default exponent 1.5).
pub fn bruggeman_effective_transport(bulk_value: f64, porosity: f64, exponent: f64) -> Result<f64> {
    require(
        bulk_value >= 0.0 && (0.0..=1.0).contains(&porosity) && exponent >= 0.0,
        "invalid Bruggeman transport controls",
    )?;
    Ok(bulk_value * porosity.powf(exponent))
}

pub fn solid_phase_ohmic_drop(current_density: f64, length: f64, conductivity: f64) -> Result<f64> {
    require(
        length >= 0.0 && conductivity > 0.0 && current_density.is_finite(),
        "invalid solid conduction controls",
    )?;
    Ok(current_density * length / conductivity)
}

/// Symmetric Marcus kinetics with reorganization energy `lambda` (J/mol).
pub fn marcus_current_density(
    exchange_current_density: f64,
    overpotential: f64,
    reorganization_energy_j_per_mol: f64,
    temperature: f64,
    electrons: f64,
) -> Result<f64> {
    let (i0, lam) = (exchange_current_density, reorganization_energy_j_per_mol);
    require(
        i0 >= 0.0 && lam > 0.0 && temperature > 0.0 && electrons > 0.0,
        "invalid Marcus kinetics controls",
    )?;
    let work = electrons * FARADAY_CONSTANT * overpotential;
    let denom = 4.0 * lam * GAS_CONSTANT * temperature;
    let forward = (-((lam - work) * (lam - work) - lam * lam) / denom).exp();
    let reverse = (-((lam + work) * (lam + work) - lam * lam) / denom).exp();
    Ok(i0 * (forward - reverse))
}

/// Marcus-Hush-Chidsey current by midpoint quadrature over the electrode Fermi distribution
/// (upstream default 1200 points, at least 100).
pub fn marcus_hush_chidsey_current_density(
    exchange_current_density: f64,
    overpotential: f64,
    reorganization_energy_j_per_mol: f64,
    temperature: f64,
    electrons: f64,
    quadrature_points: usize,
) -> Result<f64> {
    let i0 = exchange_current_density;
    require(
        i0 >= 0.0
            && reorganization_energy_j_per_mol > 0.0
            && temperature > 0.0
            && electrons > 0.0
            && quadrature_points >= 100,
        "invalid MHC kinetics controls",
    )?;
    let rt = GAS_CONSTANT * temperature;
    let lam = reorganization_energy_j_per_mol / rt;
    let v = electrons * FARADAY_CONSTANT * overpotential / rt;
    let (lo, hi) = (-20.0, 20.0);
    let h = (hi - lo) / quadrature_points as f64;
    let integral = |shift: f64| {
        let mut q = 0.0;
        for k in 0..quadrature_points {
            let x = lo + (k as f64 + 0.5) * h;
            let fermi = 1.0 / (1.0 + x.exp());
            q += (-(x - lam + shift) * (x - lam + shift) / (4.0 * lam)).exp() * fermi;
        }
        q * h
    };
    let (forward, reverse, norm) = (integral(-v), integral(v), integral(0.0));
    Ok(i0 * (forward - reverse) / norm.max(1e-300))
}

/// Finite-size (Bikerman) activity `c / (1 - v sum c)`.
pub fn bikerman_activity_correction(
    concentration: f64,
    occupied_volume_per_mole: f64,
    total_concentration: f64,
) -> Result<f64> {
    require(
        concentration >= 0.0 && occupied_volume_per_mole >= 0.0 && total_concentration >= 0.0,
        "invalid finite-size PNP state",
    )?;
    let free = 1.0 - occupied_volume_per_mole * total_concentration;
    require(free > 0.0, "finite-size PNP packing fraction >= 1")?;
    Ok(concentration / free)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModifiedPnpFlux {
    pub molar_flux: f64,
    pub activity: f64,
}

/// Face flux with steric (Bikerman) diffusion, migration and upwind-free advection.
#[allow(clippy::too_many_arguments)]
pub fn modified_nernst_planck_flux(
    concentration_left: f64,
    concentration_right: f64,
    diffusivity: f64,
    potential_gradient: f64,
    charge: i32,
    temperature: f64,
    cell_width: f64,
    occupied_volume_per_mole: f64,
    total_concentration: f64,
    advective_velocity: f64,
) -> Result<ModifiedPnpFlux> {
    let (l, r) = (concentration_left, concentration_right);
    require(
        l >= 0.0
            && r >= 0.0
            && diffusivity >= 0.0
            && temperature > 0.0
            && cell_width > 0.0
            && occupied_volume_per_mole >= 0.0
            && total_concentration >= 0.0
            && potential_gradient.is_finite()
            && advective_velocity.is_finite(),
        "invalid modified PNP flux controls",
    )?;
    let free = 1.0 - occupied_volume_per_mole * total_concentration;
    require(free > 0.0, "modified PNP packing fraction >= 1")?;
    let c = 0.5 * (l + r);
    let steric = (r - l) / cell_width / free;
    let electric = f64::from(charge) * FARADAY_CONSTANT / (GAS_CONSTANT * temperature)
        * c
        * potential_gradient;
    Ok(ModifiedPnpFlux {
        molar_flux: advective_velocity * c - diffusivity * (steric + electric),
        activity: c / free,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ImmersedInterfaceWeight {
    pub electrolyte_fraction: f64,
    pub delta: f64,
}

/// Smoothed Heaviside/delta of a signed distance (electrode negative).
pub fn smooth_electrode_interface(
    signed_distance: f64,
    thickness: f64,
) -> Result<ImmersedInterfaceWeight> {
    require(thickness > 0.0, "invalid immersed-interface thickness")?;
    let q = signed_distance / thickness;
    let c = q.cosh();
    Ok(ImmersedInterfaceWeight {
        electrolyte_fraction: 0.5 * (1.0 + q.tanh()),
        delta: 0.5 / (thickness * c * c),
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CutCellInterface {
    pub electrolyte_fraction: f64,
    pub electrode_fraction: f64,
    pub interface_area_fraction: f64,
}

/// Conservative 1-D cut cell from end-point signed distances (negative in the electrode).
pub fn cut_cell_electrode_interface(
    left_signed_distance: f64,
    right_signed_distance: f64,
    cell_width: f64,
) -> Result<CutCellInterface> {
    let (l, r) = (left_signed_distance, right_signed_distance);
    require(
        cell_width > 0.0 && l.is_finite() && r.is_finite(),
        "invalid cut-cell electrode geometry",
    )?;
    let (le, re) = (l < 0.0, r < 0.0);
    if le == re {
        return Ok(if le {
            CutCellInterface {
                electrolyte_fraction: 0.0,
                electrode_fraction: 1.0,
                interface_area_fraction: 0.0,
            }
        } else {
            CutCellInterface {
                electrolyte_fraction: 1.0,
                electrode_fraction: 0.0,
                interface_area_fraction: 0.0,
            }
        });
    }
    let fraction = (l / (l - r)).clamp(0.0, 1.0);
    let electrode = if le { fraction } else { 1.0 - fraction };
    Ok(CutCellInterface {
        electrolyte_fraction: 1.0 - electrode,
        electrode_fraction: electrode,
        interface_area_fraction: 1.0,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CoupledElectrodeResult {
    pub surface_concentration: f64,
    pub current_density: f64,
    pub molar_flux: f64,
    pub iterations: usize,
}

/// Anodic dissolution limited by film mass transfer: `k_m (c_s - c_b) = i / (nF)` with a
/// concentration-dependent Nernst equilibrium, by bisection on `c_s`.
#[allow(clippy::too_many_arguments)]
pub fn solve_anodic_bv_mass_transfer(
    bulk_concentration: f64,
    mass_transfer_coefficient: f64,
    standard_potential: f64,
    electrode_potential: f64,
    exchange_current_density: f64,
    temperature: f64,
    electrons: f64,
    reference_concentration: f64,
) -> Result<CoupledElectrodeResult> {
    let (bulk, km, cref) = (
        bulk_concentration,
        mass_transfer_coefficient,
        reference_concentration,
    );
    require(
        bulk > 0.0
            && km > 0.0
            && exchange_current_density >= 0.0
            && temperature > 0.0
            && electrons > 0.0
            && cref > 0.0,
        "invalid coupled electrode controls",
    )?;
    let nf = electrons * FARADAY_CONSTANT;
    let current = |cs: f64| -> Result<f64> {
        let eq = standard_potential + (GAS_CONSTANT * temperature / nf) * (cs / cref).ln();
        bv(
            exchange_current_density,
            electrode_potential - eq,
            temperature,
            electrons,
        )
    };
    let residual = |cs: f64| -> Result<f64> { Ok(km * (cs - bulk) - current(cs)? / nf) };
    let (mut lo, mut hi) = (bulk, (2.0 * bulk).max(bulk + 1.0));
    let (mut flo, mut fhi) = (residual(lo)?, residual(hi)?);
    let mut grow = 0;
    while grow < 80 && flo * fhi > 0.0 {
        hi *= 2.0;
        fhi = residual(hi)?;
        grow += 1;
    }
    if flo * fhi > 0.0 {
        return Err(SolverError::NumericalFailure(
            "could not bracket coupled electrode concentration".into(),
        ));
    }
    let mut it = 0;
    while it < 100 {
        let mid = 0.5 * (lo + hi);
        let fm = residual(mid)?;
        if fm.abs() < 1e-12 * (km * bulk).max(1.0) {
            lo = mid;
            hi = mid;
            break;
        }
        if flo * fm <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            flo = fm;
        }
        it += 1;
    }
    let cs = 0.5 * (lo + hi);
    let i = current(cs)?;
    Ok(CoupledElectrodeResult {
        surface_concentration: cs,
        current_density: i,
        molar_flux: i / nf,
        iterations: it + 1,
    })
}

/// Nernst potential with `Q = [H+]^-nu` from the acid speciation equilibrium.
pub fn acid_speciation_redox_potential(
    standard_potential: f64,
    temperature: f64,
    electrons: f64,
    proton_stoichiometry: f64,
    families: &[AcidFamily],
    strong_charge_mol_per_litre: f64,
) -> Result<f64> {
    require(proton_stoichiometry >= 0.0, "invalid proton stoichiometry")?;
    let eq = equilibrate_acids(families, strong_charge_mol_per_litre, 1e-14)?;
    nernst_potential(
        standard_potential,
        temperature,
        electrons,
        eq.hydrogen_mol_per_litre.powf(-proton_stoichiometry),
    )
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PittingCriterion {
    pub critical_chloride: f64,
    pub critical_potential: f64,
}

pub fn pitting_initiates(
    chloride_activity: f64,
    electrode_potential: f64,
    c: &PittingCriterion,
) -> Result<bool> {
    require(
        c.critical_chloride >= 0.0 && c.critical_potential.is_finite(),
        "invalid pitting criterion",
    )?;
    Ok(chloride_activity >= c.critical_chloride && electrode_potential >= c.critical_potential)
}

/// `i0 exp(sigma_h V / RT)`.
pub fn stress_assisted_exchange_current(
    base_exchange_current: f64,
    hydrostatic_stress: f64,
    activation_volume: f64,
    temperature: f64,
) -> Result<f64> {
    require(
        base_exchange_current >= 0.0
            && temperature > 0.0
            && hydrostatic_stress.is_finite()
            && activation_volume.is_finite(),
        "invalid stress-corrosion state",
    )?;
    Ok(base_exchange_current
        * (hydrostatic_stress * activation_volume / (GAS_CONSTANT * temperature)).exp())
}

/// One explicit Allen-Cahn corrosion update clamped to `[-1, 1]`.
pub fn phase_field_corrosion_step(
    order_parameter: f64,
    laplacian: f64,
    electrochemical_driving: f64,
    mobility: f64,
    gradient_energy: f64,
    dt: f64,
) -> Result<f64> {
    require(
        mobility >= 0.0 && gradient_energy >= 0.0 && dt >= 0.0,
        "invalid phase-field controls",
    )?;
    let p = order_parameter;
    let derivative = p * (p * p - 1.0) - gradient_energy * laplacian - electrochemical_driving;
    Ok((p - mobility * dt * derivative).clamp(-1.0, 1.0))
}

pub fn recessed_level_set(signed_distance: f64, recession_distance: f64) -> Result<f64> {
    require(recession_distance >= 0.0, "invalid recession distance")?;
    Ok(signed_distance + recession_distance)
}

/// Affinely compresses increasing coordinates toward the first one by the recession distance.
pub fn recession_remesh_1d(coordinates: &[f64], recession_distance: f64) -> Result<Vec<f64>> {
    require(
        coordinates.len() >= 2 && recession_distance >= 0.0 && recession_distance.is_finite(),
        "invalid recession remesh controls",
    )?;
    require(
        coordinates
            .windows(2)
            .all(|w| w[1] > w[0] && w[1].is_finite()),
        "remesh coordinates must be finite and increasing",
    )?;
    let length = coordinates[coordinates.len() - 1] - coordinates[0];
    require(
        recession_distance < length,
        "recession consumes the mesh interval",
    )?;
    let scale = (length - recession_distance) / length;
    Ok(coordinates
        .iter()
        .map(|x| coordinates[0] + scale * (x - coordinates[0]))
        .collect())
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StressCorrosionUpdate {
    pub stress_amplified_current_density: f64,
    pub recession_distance: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn coupled_stress_corrosion_step(
    base_exchange_current: f64,
    hydrostatic_stress: f64,
    activation_volume: f64,
    temperature: f64,
    dt: f64,
    molar_mass: f64,
    metal_density: f64,
    electrons: f64,
) -> Result<StressCorrosionUpdate> {
    require(
        base_exchange_current >= 0.0
            && hydrostatic_stress.is_finite()
            && activation_volume.is_finite()
            && temperature > 0.0
            && dt >= 0.0
            && molar_mass > 0.0
            && metal_density > 0.0
            && electrons > 0.0,
        "invalid stress-corrosion coupling controls",
    )?;
    let i = stress_assisted_exchange_current(
        base_exchange_current,
        hydrostatic_stress,
        activation_volume,
        temperature,
    )?;
    Ok(StressCorrosionUpdate {
        stress_amplified_current_density: i,
        recession_distance: faradaic_recession_distance(
            i,
            dt,
            molar_mass,
            metal_density,
            electrons,
        )?,
    })
}

/// Clamped-boundary (zero-gradient) 2-D phase-field corrosion.
#[derive(Clone, Debug)]
pub struct PhaseFieldCorrosion2D {
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    mobility: f64,
    gradient: f64,
    phi: Vec<f64>,
}

impl PhaseFieldCorrosion2D {
    pub fn new(
        nx: usize,
        ny: usize,
        dx: f64,
        dy: f64,
        mobility: f64,
        gradient_energy: f64,
    ) -> Result<Self> {
        require(
            nx >= 3 && ny >= 3 && dx > 0.0 && dy > 0.0 && mobility >= 0.0 && gradient_energy >= 0.0,
            "invalid 2-D phase-field grid",
        )?;
        Ok(Self {
            nx,
            ny,
            dx,
            dy,
            mobility,
            gradient: gradient_energy,
            phi: vec![0.0; nx * ny],
        })
    }
    /// `tanh((w - |x - x0|)/w)`: a disc of corroding phase (+1) around `(x0, y0)`.
    pub fn initialize_interface(&mut self, x0: f64, y0: f64, width: f64) -> Result<()> {
        require(
            width > 0.0 && x0.is_finite() && y0.is_finite(),
            "invalid 2-D phase-field interface",
        )?;
        for j in 0..self.ny {
            for i in 0..self.nx {
                let (x, y) = (i as f64 * self.dx, j as f64 * self.dy);
                let d = ((x - x0) * (x - x0) + (y - y0) * (y - y0)).sqrt();
                self.phi[j * self.nx + i] = ((width - d) / width).tanh();
            }
        }
        Ok(())
    }
    pub fn step(&mut self, dt: f64, driving: &[f64]) -> Result<()> {
        require(
            dt >= 0.0 && driving.len() == self.phi.len(),
            "2-D phase-field step mismatch",
        )?;
        let (nx, ny) = (self.nx, self.ny);
        let mut next = vec![0.0; self.phi.len()];
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let il = j * nx + if i > 0 { i - 1 } else { i };
                let ir = j * nx + if i + 1 < nx { i + 1 } else { i };
                let ib = if j > 0 { j - 1 } else { j } * nx + i;
                let it = if j + 1 < ny { j + 1 } else { j } * nx + i;
                let p = &self.phi;
                let lap = (p[il] - 2.0 * p[k] + p[ir]) / (self.dx * self.dx)
                    + (p[ib] - 2.0 * p[k] + p[it]) / (self.dy * self.dy);
                next[k] = phase_field_corrosion_step(
                    p[k],
                    lap,
                    driving[k],
                    self.mobility,
                    self.gradient,
                    dt,
                )?;
            }
        }
        self.phi = next;
        Ok(())
    }
    pub fn order_parameter(&self) -> &[f64] {
        &self.phi
    }
    pub fn nx(&self) -> usize {
        self.nx
    }
    pub fn ny(&self) -> usize {
        self.ny
    }
}

#[derive(Clone, Debug)]
pub struct PhaseFieldCorrosion1D {
    dx: f64,
    mobility: f64,
    gradient: f64,
    phi: Vec<f64>,
}

impl PhaseFieldCorrosion1D {
    pub fn new(cells: usize, dx: f64, mobility: f64, gradient_energy: f64) -> Result<Self> {
        require(
            cells >= 3 && dx > 0.0 && mobility >= 0.0 && gradient_energy >= 0.0,
            "invalid phase-field grid",
        )?;
        Ok(Self {
            dx,
            mobility,
            gradient: gradient_energy,
            phi: vec![0.0; cells],
        })
    }
    pub fn initialize_interface(&mut self, location: f64, width: f64) -> Result<()> {
        require(width > 0.0, "invalid phase-field width")?;
        for (i, p) in self.phi.iter_mut().enumerate() {
            *p = ((location - i as f64 * self.dx) / width).tanh();
        }
        Ok(())
    }
    pub fn step(&mut self, dt: f64, driving: &[f64]) -> Result<()> {
        let n = self.phi.len();
        require(dt >= 0.0 && driving.len() == n, "phase-field step mismatch")?;
        let mut next = vec![0.0; n];
        for i in 0..n {
            let l = if i > 0 { i - 1 } else { i };
            let r = if i + 1 < n { i + 1 } else { i };
            let lap = (self.phi[l] - 2.0 * self.phi[i] + self.phi[r]) / (self.dx * self.dx);
            next[i] = phase_field_corrosion_step(
                self.phi[i],
                lap,
                driving[i],
                self.mobility,
                self.gradient,
                dt,
            )?;
        }
        self.phi = next;
        Ok(())
    }
    pub fn order_parameter(&self) -> &[f64] {
        &self.phi
    }
}
