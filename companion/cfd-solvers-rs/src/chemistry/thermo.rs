use super::GAS_CONSTANT;
use crate::{positive, require, Result, SolverError};
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ThermoSpecies {
    pub name: String,
    pub standard_chemical_potential: f64,
}
pub trait PhaseThermo {
    fn species(&self) -> &[ThermoSpecies];
    fn activities(
        &self,
        mole_fractions: &[f64],
        temperature: f64,
        pressure: f64,
    ) -> Result<Vec<f64>>;
    fn species_count(&self) -> usize {
        self.species().len()
    }
    fn chemical_potentials(&self, x: &[f64], temperature: f64, pressure: f64) -> Result<Vec<f64>> {
        let a = self.activities(x, temperature, pressure)?;
        a.iter()
            .zip(self.species())
            .map(|(&a, s)| {
                positive(a, "chemical potential activity")?;
                Ok(s.standard_chemical_potential + GAS_CONSTANT * temperature * a.ln())
            })
            .collect()
    }
}
fn normalize(x: &[f64], count: usize, temperature: f64, pressure: f64) -> Result<Vec<f64>> {
    positive(temperature, "temperature")?;
    positive(pressure, "pressure")?;
    require(
        !x.is_empty() && x.len() == count && x.iter().all(|v| v.is_finite() && *v >= 0.0),
        "invalid composition",
    )?;
    let sum: f64 = x.iter().sum();
    positive(sum, "composition sum")?;
    Ok(x.iter().map(|v| v / sum).collect())
}
#[derive(Clone, Debug)]
pub struct IdealGasPhase {
    species: Vec<ThermoSpecies>,
    reference_pressure: f64,
}
impl IdealGasPhase {
    pub fn new(species: Vec<ThermoSpecies>, reference_pressure: f64) -> Result<Self> {
        positive(reference_pressure, "reference pressure")?;
        require(
            !species.is_empty()
                && species
                    .iter()
                    .all(|s| s.standard_chemical_potential.is_finite()),
            "invalid phase species",
        )?;
        Ok(Self {
            species,
            reference_pressure,
        })
    }
}
impl PhaseThermo for IdealGasPhase {
    fn species(&self) -> &[ThermoSpecies] {
        &self.species
    }
    fn activities(&self, x: &[f64], t: f64, p: f64) -> Result<Vec<f64>> {
        let mut a = normalize(x, self.species.len(), t, p)?;
        for v in &mut a {
            *v *= p / self.reference_pressure;
        }
        Ok(a)
    }
}
#[derive(Clone, Debug)]
pub struct IdealSolutionPhase {
    species: Vec<ThermoSpecies>,
}
impl IdealSolutionPhase {
    pub fn new(species: Vec<ThermoSpecies>) -> Result<Self> {
        require(
            !species.is_empty()
                && species
                    .iter()
                    .all(|s| s.standard_chemical_potential.is_finite()),
            "invalid phase species",
        )?;
        Ok(Self { species })
    }
}
impl PhaseThermo for IdealSolutionPhase {
    fn species(&self) -> &[ThermoSpecies] {
        &self.species
    }
    fn activities(&self, x: &[f64], t: f64, p: f64) -> Result<Vec<f64>> {
        normalize(x, self.species.len(), t, p)
    }
}
pub fn parse_thermo_csv(text: &str) -> Result<Vec<ThermoSpecies>> {
    let mut out = vec![];
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (name, value) = line
            .split_once(',')
            .ok_or_else(|| SolverError::InvalidInput("invalid thermo CSV row".into()))?;
        let mu = value
            .trim()
            .parse::<f64>()
            .map_err(|_| SolverError::InvalidInput("invalid thermo CSV value".into()))?;
        require(
            !name.is_empty() && mu.is_finite(),
            "invalid thermo CSV species",
        )?;
        out.push(ThermoSpecies {
            name: name.into(),
            standard_chemical_potential: mu,
        });
    }
    require(!out.is_empty(), "empty thermo CSV")?;
    Ok(out)
}
pub fn load_thermo_csv(path: impl AsRef<std::path::Path>) -> Result<Vec<ThermoSpecies>> {
    parse_thermo_csv(&std::fs::read_to_string(path)?)
}
#[derive(Clone, Copy, Debug)]
pub struct BinaryPrecipitationResult {
    pub dissolved_a: f64,
    pub dissolved_b: f64,
    pub solid_amount: f64,
    pub reaction_extent: f64,
    pub saturation_ratio: f64,
}
pub fn equilibrate_binary_salt(
    a: f64,
    b: f64,
    solid: f64,
    ksp: f64,
) -> Result<BinaryPrecipitationResult> {
    positive(ksp, "solubility product")?;
    require(
        [a, b, solid].iter().all(|v| v.is_finite() && *v >= 0.0),
        "invalid binary-salt amounts",
    )?;
    let eq = 0.5 * (a + b - ((a - b) * (a - b) + 4.0 * ksp).sqrt());
    let extent = eq.max(-solid).min(a.min(b));
    let (af, bf, sf) = (a - extent, b - extent, solid + extent);
    Ok(BinaryPrecipitationResult {
        dissolved_a: af,
        dissolved_b: bf,
        solid_amount: sf,
        reaction_extent: extent,
        saturation_ratio: af * bf / ksp,
    })
}
