use super::GAS_CONSTANT;
use crate::{positive, require, Result};
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Species {
    pub name: String,
    pub molar_mass: f64,
    pub charge: i32,
}
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
pub struct ArrheniusRate {
    pub pre_exponential: f64,
    pub temperature_exponent: f64,
    pub activation_energy: f64,
}
impl ArrheniusRate {
    pub fn rate_constant(self, temperature: f64) -> Result<f64> {
        positive(temperature, "temperature")?;
        require(
            self.pre_exponential.is_finite()
                && self.pre_exponential >= 0.0
                && self.temperature_exponent.is_finite()
                && self.activation_energy.is_finite(),
            "invalid Arrhenius parameters",
        )?;
        let k = self.pre_exponential
            * temperature.powf(self.temperature_exponent)
            * (-self.activation_energy / (GAS_CONSTANT * temperature)).exp();
        require(k.is_finite(), "Arrhenius overflow")?;
        Ok(k)
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct StoichiometricTerm {
    pub species: usize,
    pub coefficient: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct ReactionOrder {
    pub species: usize,
    pub order: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct ThirdBodyEfficiency {
    pub species: usize,
    pub efficiency: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct EquilibriumConstant {
    pub reference_value: f64,
    pub reference_temperature: f64,
    pub reaction_enthalpy: f64,
}
impl Default for EquilibriumConstant {
    fn default() -> Self {
        Self {
            reference_value: 1.0,
            reference_temperature: 298.15,
            reaction_enthalpy: 0.0,
        }
    }
}
impl EquilibriumConstant {
    pub fn value(self, temperature: f64) -> Result<f64> {
        positive(temperature, "temperature")?;
        positive(self.reference_value, "equilibrium constant")?;
        positive(self.reference_temperature, "reference temperature")?;
        require(self.reaction_enthalpy.is_finite(), "nonfinite enthalpy")?;
        let k = self.reference_value
            * (-self.reaction_enthalpy / GAS_CONSTANT
                * (1.0 / temperature - 1.0 / self.reference_temperature))
                .exp();
        positive(k, "equilibrium constant result")?;
        Ok(k)
    }
}
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ElementaryReaction {
    pub reactants: Vec<StoichiometricTerm>,
    pub products: Vec<StoichiometricTerm>,
    pub forward: ArrheniusRate,
    pub equilibrium: Option<EquilibriumConstant>,
    pub forward_orders: Vec<ReactionOrder>,
    pub third_body: bool,
    pub third_body_efficiencies: Vec<ThirdBodyEfficiency>,
    pub low_pressure_limit: Option<ArrheniusRate>,
}
#[derive(Clone, Debug)]
pub struct ReactionNetwork {
    species: Vec<Species>,
    reactions: Vec<ElementaryReaction>,
}
impl ReactionNetwork {
    pub fn new(species: Vec<Species>) -> Self {
        Self {
            species,
            reactions: vec![],
        }
    }
    pub fn species(&self) -> &[Species] {
        &self.species
    }
    pub fn reactions(&self) -> &[ElementaryReaction] {
        &self.reactions
    }
    pub fn add_species(&mut self, species: Species) -> usize {
        self.species.push(species);
        self.species.len() - 1
    }
    pub fn add_reaction(&mut self, r: ElementaryReaction) -> Result<()> {
        require(
            !r.reactants.is_empty() || !r.products.is_empty(),
            "empty reaction",
        )?;
        for t in r.reactants.iter().chain(&r.products) {
            require(
                t.species < self.species.len(),
                "stoichiometric species index",
            )?;
            positive(t.coefficient, "stoichiometric coefficient")?;
        }
        for t in &r.forward_orders {
            require(
                t.species < self.species.len() && t.order.is_finite(),
                "invalid reaction order",
            )?;
        }
        for t in &r.third_body_efficiencies {
            require(
                t.species < self.species.len() && t.efficiency.is_finite() && t.efficiency >= 0.0,
                "invalid third-body efficiency",
            )?;
        }
        self.reactions.push(r);
        Ok(())
    }
    pub fn reaction_rates_into(
        &self,
        c: &[f64],
        temperature: f64,
        rates: &mut [f64],
    ) -> Result<()> {
        require(
            c.len() == self.species.len() && rates.len() == self.reactions.len(),
            "reaction buffer sizes",
        )?;
        positive(temperature, "temperature")?;
        require(
            c.iter().all(|v| v.is_finite() && *v >= 0.0),
            "invalid concentration",
        )?;
        for (i, r) in self.reactions.iter().enumerate() {
            let kinf = r.forward.rate_constant(temperature)?;
            let mut collider = 1.0;
            if r.third_body || r.low_pressure_limit.is_some() {
                collider = c.iter().sum();
                for e in &r.third_body_efficiencies {
                    collider += (e.efficiency - 1.0) * c[e.species];
                }
                require(
                    collider.is_finite() && collider >= 0.0,
                    "third-body concentration invalid",
                )?;
            }
            let mut k = kinf;
            if let Some(k0) = r.low_pressure_limit {
                if kinf == 0.0 {
                    k = 0.0;
                } else {
                    let pr = k0.rate_constant(temperature)? * collider / kinf;
                    k = kinf * pr / (1.0 + pr);
                }
            } else if r.third_body {
                k *= collider;
            }
            let mut rate = k;
            if !r.forward_orders.is_empty() {
                for o in &r.forward_orders {
                    rate *= c[o.species].powf(o.order);
                }
            } else {
                for t in &r.reactants {
                    rate *= c[t.species].powf(t.coefficient);
                }
            }
            if let Some(eq) = r.equilibrium {
                let mut reverse = k / eq.value(temperature)?;
                for t in &r.products {
                    reverse *= c[t.species].powf(t.coefficient);
                }
                rate -= reverse;
            }
            require(rate.is_finite(), "nonfinite reaction rate")?;
            rates[i] = rate;
        }
        Ok(())
    }
    pub fn reaction_rates(&self, c: &[f64], temperature: f64) -> Result<Vec<f64>> {
        let mut rates = vec![0.0; self.reactions.len()];
        self.reaction_rates_into(c, temperature, &mut rates)?;
        Ok(rates)
    }
    pub fn source_terms_into(
        &self,
        c: &[f64],
        temperature: f64,
        source: &mut [f64],
        rates: &mut [f64],
    ) -> Result<()> {
        require(source.len() == self.species.len(), "source buffer size")?;
        self.reaction_rates_into(c, temperature, rates)?;
        source.fill(0.0);
        for (i, r) in self.reactions.iter().enumerate() {
            for t in &r.reactants {
                source[t.species] -= t.coefficient * rates[i];
            }
            for t in &r.products {
                source[t.species] += t.coefficient * rates[i];
            }
        }
        crate::finite_slice(source, "chemical source")
    }
    pub fn source_terms(&self, c: &[f64], temperature: f64) -> Result<Vec<f64>> {
        let mut source = vec![0.0; self.species.len()];
        let mut rates = vec![0.0; self.reactions.len()];
        self.source_terms_into(c, temperature, &mut source, &mut rates)?;
        Ok(source)
    }
}
