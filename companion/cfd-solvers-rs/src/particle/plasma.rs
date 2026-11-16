//! Plasma chemistry (mass-action explicit substeps), Paschen breakdown and parallel-plate
//! multipactor estimates. Ported from upstream `cfd/particle/plasma` (MIT, cfd_solvers
//! contributors).
use crate::{require, Result, SolverError};

const ELEMENTARY_CHARGE: f64 = 1.602176634e-19;

#[derive(Clone, Debug, PartialEq)]
pub struct PlasmaSpecies {
    pub name: String,
    pub charge_c: f64,
    pub mass_kg: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PlasmaState {
    pub species: Vec<PlasmaSpecies>,
    pub number_density_m3: Vec<f64>,
    pub electron_temperature_ev: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlasmaReactionTerm {
    pub species: usize,
    pub coefficient: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PlasmaReaction {
    pub reactants: Vec<PlasmaReactionTerm>,
    pub products: Vec<PlasmaReactionTerm>,
    pub rate_coefficient: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PlasmaChemistryDiagnostics {
    pub reaction_rates_m3_s: Vec<f64>,
    pub charge_density_before_c_m3: f64,
    pub charge_density_after_c_m3: f64,
    pub total_number_density_before_m3: f64,
    pub total_number_density_after_m3: f64,
    pub max_relative_density_change: f64,
}

fn validate_state(s: &PlasmaState) -> Result<()> {
    require(
        !s.species.is_empty() && s.species.len() == s.number_density_m3.len(),
        "invalid plasma state size",
    )?;
    for (sp, &n) in s.species.iter().zip(&s.number_density_m3) {
        require(
            sp.mass_kg > 0.0 && sp.charge_c.is_finite(),
            "invalid plasma species",
        )?;
        require(n >= 0.0 && n.is_finite(), "invalid plasma density")?;
    }
    Ok(())
}

fn reaction_rate(s: &PlasmaState, r: &PlasmaReaction) -> Result<f64> {
    require(
        r.rate_coefficient >= 0.0 && r.rate_coefficient.is_finite(),
        "invalid plasma reaction rate coefficient",
    )?;
    let mut rate = r.rate_coefficient;
    for t in &r.reactants {
        require(
            t.species < s.number_density_m3.len()
                && t.coefficient > 0.0
                && t.coefficient.is_finite(),
            "invalid plasma reactant term",
        )?;
        rate *= s.number_density_m3[t.species].powf(t.coefficient);
    }
    for t in &r.products {
        require(
            t.species < s.number_density_m3.len()
                && t.coefficient > 0.0
                && t.coefficient.is_finite(),
            "invalid plasma product term",
        )?;
    }
    Ok(rate)
}

pub fn plasma_charge_density_c_m3(s: &PlasmaState) -> Result<f64> {
    validate_state(s)?;
    Ok(s.species
        .iter()
        .zip(&s.number_density_m3)
        .fold(0.0, |acc, (sp, n)| acc + sp.charge_c * n))
}

/// Explicit mass-action substeps; a substep driving any density clearly negative is an error
/// and the state is left unchanged.
pub fn advance_plasma_chemistry(
    state: &mut PlasmaState,
    reactions: &[PlasmaReaction],
    dt: f64,
    substeps: usize,
) -> Result<PlasmaChemistryDiagnostics> {
    validate_state(state)?;
    require(
        !reactions.is_empty() && dt >= 0.0 && dt.is_finite() && substeps > 0,
        "invalid plasma chemistry controls",
    )?;
    let mut work = state.clone();
    let mut d = PlasmaChemistryDiagnostics {
        charge_density_before_c_m3: plasma_charge_density_c_m3(&work)?,
        total_number_density_before_m3: work.number_density_m3.iter().sum(),
        reaction_rates_m3_s: vec![0.0; reactions.len()],
        ..Default::default()
    };
    let initial = work.number_density_m3.clone();
    let dt_sub = dt / substeps as f64;
    for _ in 0..substeps {
        let mut delta = vec![0.0; work.number_density_m3.len()];
        for (ri, r) in reactions.iter().enumerate() {
            let rate = reaction_rate(&work, r)?;
            d.reaction_rates_m3_s[ri] = rate;
            for t in &r.reactants {
                delta[t.species] -= t.coefficient * rate * dt_sub;
            }
            for t in &r.products {
                delta[t.species] += t.coefficient * rate * dt_sub;
            }
        }
        for (n, dn) in work.number_density_m3.iter_mut().zip(&delta) {
            let updated = *n + dn;
            if updated < -1.0e-9 * n.max(1.0) {
                return Err(SolverError::NumericalFailure(
                    "plasma chemistry step produced negative density".into(),
                ));
            }
            *n = updated.max(0.0);
        }
    }
    d.charge_density_after_c_m3 = plasma_charge_density_c_m3(&work)?;
    d.total_number_density_after_m3 = work.number_density_m3.iter().sum();
    for (n, n0) in work.number_density_m3.iter().zip(&initial) {
        d.max_relative_density_change = d
            .max_relative_density_change
            .max((n - n0).abs() / n0.max(1.0));
    }
    *state = work;
    Ok(d)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PaschenGas {
    pub townsend_a_per_m_pa: f64,
    pub townsend_b_v_per_m_pa: f64,
    pub secondary_emission_yield: f64,
}

impl Default for PaschenGas {
    fn default() -> Self {
        Self {
            townsend_a_per_m_pa: 0.0,
            townsend_b_v_per_m_pa: 0.0,
            secondary_emission_yield: 0.01,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GasBreakdownEstimate {
    pub breakdown_voltage_v: f64,
    pub breakdown_field_v_m: f64,
    pub applied_field_exceeds_threshold: bool,
}

/// `V_b = B·pd / (ln(A·pd) − ln(ln(1 + 1/γ)))`.
pub fn paschen_breakdown_voltage_v(gas: &PaschenGas, pressure: f64, gap: f64) -> Result<f64> {
    require(
        gas.townsend_a_per_m_pa > 0.0
            && gas.townsend_b_v_per_m_pa > 0.0
            && gas.secondary_emission_yield > 0.0
            && pressure > 0.0
            && gap > 0.0,
        "invalid Paschen controls",
    )?;
    let pd = pressure * gap;
    let argument = gas.townsend_a_per_m_pa * pd;
    let secondary = (1.0 + 1.0 / gas.secondary_emission_yield).ln();
    require(
        argument > 1.0 && secondary > 1.0,
        "Paschen controls outside breakdown-law domain",
    )?;
    let denominator = argument.ln() - secondary.ln();
    require(denominator > 0.0, "Paschen denominator is not positive")?;
    Ok(gas.townsend_b_v_per_m_pa * pd / denominator)
}

pub fn estimate_gas_breakdown_threshold(
    gas: &PaschenGas,
    pressure: f64,
    gap: f64,
    applied_field: f64,
) -> Result<GasBreakdownEstimate> {
    require(
        applied_field >= 0.0 && applied_field.is_finite(),
        "invalid applied gas-breakdown field",
    )?;
    let v = paschen_breakdown_voltage_v(gas, pressure, gap)?;
    let e = v / gap;
    Ok(GasBreakdownEstimate {
        breakdown_voltage_v: v,
        breakdown_field_v_m: e,
        applied_field_exceeds_threshold: applied_field >= e,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MultipactorEstimate {
    pub order: usize,
    pub transit_time_s: f64,
    pub resonant_field_v_m: f64,
    pub impact_energy_ev: f64,
    pub secondary_yield_sustains: bool,
}

/// Odd-half-period resonance `t = (2N−1)/(2f)` for a parallel-plate gap.
pub fn estimate_parallel_plate_multipactor(
    frequency: f64,
    gap: f64,
    order: usize,
    mass: f64,
    charge: f64,
    secondary_yield: f64,
) -> Result<MultipactorEstimate> {
    require(
        frequency > 0.0
            && gap > 0.0
            && order > 0
            && mass > 0.0
            && charge.abs() > 0.0
            && secondary_yield >= 0.0,
        "invalid multipactor estimate controls",
    )?;
    let transit = (2 * order - 1) as f64 / (2.0 * frequency);
    let q = charge.abs();
    let impact_velocity = 2.0 * gap / transit;
    Ok(MultipactorEstimate {
        order,
        transit_time_s: transit,
        resonant_field_v_m: 2.0 * mass * gap / (q * transit * transit),
        impact_energy_ev: 0.5 * mass * impact_velocity * impact_velocity / ELEMENTARY_CHARGE,
        secondary_yield_sustains: secondary_yield > 1.0,
    })
}
