//! Ideal dilute acid/base speciation by charge-balance bisection in `ln[H+]`, ported from
//! upstream `cfd/chemistry/aqueous_equilibrium` (MIT, cfd_solvers contributors).
use crate::{require, Result, SolverError};

/// Concentrations and Ka in mol/L.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AcidFamily {
    pub total_mol_per_litre: f64,
    pub fully_protonated_charge: i32,
    pub dissociation_constants: Vec<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AqueousEquilibriumResult {
    pub ph: f64,
    pub hydrogen_mol_per_litre: f64,
    pub hydroxide_mol_per_litre: f64,
    pub charge_residual_mol_per_litre: f64,
    pub family_concentrations: Vec<Vec<f64>>,
}

fn distribution(family: &AcidFamily, log_h: f64) -> Vec<f64> {
    let mut w = vec![0.0; family.dissociation_constants.len() + 1];
    for i in 1..w.len() {
        w[i] = w[i - 1] + family.dissociation_constants[i - 1].ln() - log_h;
    }
    let maximum = w.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut sum = 0.0;
    for v in &mut w {
        *v = (*v - maximum).exp();
        sum += *v;
    }
    for v in &mut w {
        *v *= family.total_mol_per_litre / sum;
    }
    w
}

/// `strong_charge` is the signed spectator-ion charge concentration; `water_ion_product` is Kw
/// (upstream default 1e-14).
pub fn equilibrate_acids(
    families: &[AcidFamily],
    strong_charge_mol_per_litre: f64,
    water_ion_product: f64,
) -> Result<AqueousEquilibriumResult> {
    let (sc, kw) = (strong_charge_mol_per_litre, water_ion_product);
    require(
        sc.is_finite() && kw > 0.0 && kw.is_finite(),
        "invalid aqueous equilibrium controls",
    )?;
    for f in families {
        require(
            f.total_mol_per_litre >= 0.0 && f.total_mol_per_litre.is_finite(),
            "invalid analytical acid concentration",
        )?;
        require(
            f.dissociation_constants
                .iter()
                .all(|k| *k > 0.0 && k.is_finite()),
            "invalid acid dissociation constant",
        )?;
    }
    let charge = |log_h: f64| {
        let mut q = log_h.exp() - kw / log_h.exp() + sc;
        for f in families {
            for (j, c) in distribution(f, log_h).iter().enumerate() {
                q += (f64::from(f.fully_protonated_charge) - j as f64) * c;
            }
        }
        q
    };
    let (mut low, mut high) = (-690.0f64, 690.0f64);
    if !(charge(low) < 0.0 && charge(high) > 0.0) {
        return Err(SolverError::NumericalFailure(
            "aqueous charge balance cannot be bracketed".into(),
        ));
    }
    for _ in 0..200 {
        let mid = 0.5 * (low + high);
        if charge(mid) > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
        if high - low < 1e-13 {
            break;
        }
    }
    let log_h = 0.5 * (low + high);
    let h = log_h.exp();
    Ok(AqueousEquilibriumResult {
        ph: -log_h / 10f64.ln(),
        hydrogen_mol_per_litre: h,
        hydroxide_mol_per_litre: kw / h,
        charge_residual_mol_per_litre: charge(log_h),
        family_concentrations: families.iter().map(|f| distribution(f, log_h)).collect(),
    })
}
