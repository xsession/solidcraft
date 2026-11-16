//! Reaction-path analysis, well-stirred reactor networks, monolithic speciation, NASA-7/Shomate
//! thermodynamics, Pitzer/Davies/SIT activity models, binary ideal Gibbs minimization, sorption,
//! and PHREEQC/Cantera database readers, ported from upstream `cfd/chemistry/{reaction_path,
//! reactor_network, speciation_block, advanced_thermo, database_io}` (MIT, cfd_solvers
//! contributors).
#![allow(clippy::needless_range_loop)]
use super::{
    integrate_isothermal, ArrheniusRate, ElementaryReaction, ReactionNetwork, ReactorConfig,
    Species, StoichiometricTerm, GAS_CONSTANT,
};
use crate::{
    core::{
        parse_float_prefix, solve_matrix_free_newton, MatrixFreeNewtonConfig,
        MatrixFreeNonlinearSystem,
    },
    require, Result, SolverError,
};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReactionPathEdge {
    pub from: usize,
    pub to: usize,
    pub molar_flux: f64,
}

/// Element-agnostic species-to-species flux graph: every reactant-product pair of a reaction
/// carries `|rate| nu_a nu_b`, oriented by the sign of the net rate.
pub fn reaction_path(
    network: &ReactionNetwork,
    concentrations: &[f64],
    temperature: f64,
    minimum_flux: f64,
) -> Result<Vec<ReactionPathEdge>> {
    require(minimum_flux >= 0.0, "negative reaction-path cutoff")?;
    let rates = network.reaction_rates(concentrations, temperature)?;
    let mut out = Vec::new();
    for (rx, &rate) in network.reactions().iter().zip(&rates) {
        let forward = rate >= 0.0;
        let (from, to) = if forward {
            (&rx.reactants, &rx.products)
        } else {
            (&rx.products, &rx.reactants)
        };
        for a in from {
            for b in to {
                let flux = rate.abs() * a.coefficient * b.coefficient;
                if flux >= minimum_flux && flux > 0.0 {
                    out.push(ReactionPathEdge {
                        from: a.species,
                        to: b.species,
                        molar_flux: flux,
                    });
                }
            }
        }
    }
    Ok(out)
}

#[derive(Clone, Debug, PartialEq)]
pub struct WellStirredReactor {
    /// m3.
    pub volume: f64,
    pub temperature: f64,
    pub concentrations: Vec<f64>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReactorConnection {
    pub first: usize,
    pub second: usize,
    /// Symmetric mixing exchange, m3/s.
    pub exchange_flow: f64,
}

/// Operator-split network: explicit symmetric mixing, then implicit isothermal chemistry per
/// reactor, in sub-steps of at most `maximum_step`.
#[derive(Clone, Debug)]
pub struct WellStirredReactorNetwork {
    chemistry: ReactionNetwork,
    reactors: Vec<WellStirredReactor>,
    connections: Vec<ReactorConnection>,
}

impl WellStirredReactorNetwork {
    pub fn new(chemistry: ReactionNetwork) -> Self {
        Self {
            chemistry,
            reactors: Vec::new(),
            connections: Vec::new(),
        }
    }
    pub fn add_reactor(&mut self, r: WellStirredReactor) -> Result<usize> {
        require(
            r.volume > 0.0
                && r.temperature > 0.0
                && r.volume.is_finite()
                && r.temperature.is_finite()
                && r.concentrations.len() == self.chemistry.species().len(),
            "invalid well-stirred reactor",
        )?;
        require(
            r.concentrations.iter().all(|c| *c >= 0.0 && c.is_finite()),
            "invalid reactor concentration",
        )?;
        self.reactors.push(r);
        Ok(self.reactors.len() - 1)
    }
    pub fn connect(&mut self, c: ReactorConnection) -> Result<()> {
        require(
            c.first < self.reactors.len()
                && c.second < self.reactors.len()
                && c.first != c.second
                && c.exchange_flow >= 0.0
                && c.exchange_flow.is_finite(),
            "invalid reactor connection",
        )?;
        self.connections.push(c);
        Ok(())
    }
    /// Upstream default `maximum_step` 1e-3 s.
    pub fn advance(
        &mut self,
        duration: f64,
        maximum_step: f64,
        controls: ReactorConfig,
    ) -> Result<()> {
        require(
            duration >= 0.0
                && maximum_step > 0.0
                && duration.is_finite()
                && maximum_step.is_finite(),
            "invalid reactor-network step",
        )?;
        let ns = self.chemistry.species().len();
        let mut time = 0.0;
        while time < duration {
            let dt = maximum_step.min(duration - time);
            let mut delta = vec![vec![0.0; ns]; self.reactors.len()];
            for c in &self.connections {
                let (a, b) = (&self.reactors[c.first], &self.reactors[c.second]);
                for s in 0..ns {
                    let rate = c.exchange_flow * (b.concentrations[s] - a.concentrations[s]);
                    delta[c.first][s] += rate * dt / a.volume;
                    delta[c.second][s] -= rate * dt / b.volume;
                }
            }
            for (r, d) in self.reactors.iter_mut().zip(&delta) {
                for s in 0..ns {
                    let v = r.concentrations[s] + d[s];
                    require(
                        v >= -1e-14,
                        "reactor mixing produced negative concentration",
                    )?;
                    r.concentrations[s] = v.max(0.0);
                }
            }
            if ns > 0 {
                for r in &mut self.reactors {
                    integrate_isothermal(
                        &self.chemistry,
                        &mut r.concentrations,
                        r.temperature,
                        dt,
                        controls,
                    )?;
                }
            }
            time += dt;
        }
        Ok(())
    }
    pub fn reactors(&self) -> &[WellStirredReactor] {
        &self.reactors
    }
    pub fn total_species_amounts(&self) -> Vec<f64> {
        let mut amount = vec![0.0; self.chemistry.species().len()];
        for r in &self.reactors {
            for (a, c) in amount.iter_mut().zip(&r.concentrations) {
                *a += r.volume * c;
            }
        }
        amount
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeciationConfig {
    pub max_iterations: usize,
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
}

impl Default for SpeciationConfig {
    fn default() -> Self {
        Self {
            max_iterations: 40,
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-10,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpeciationResult {
    pub concentrations: Vec<f64>,
    pub iterations: usize,
    pub residual_rms: f64,
    pub converged: bool,
}

/// One conserved component over `K + 1` species in log-concentration unknowns:
/// `ln c_i - ln c_0 = ln K_i` and `sum c = total`.
#[derive(Clone, Copy, Debug)]
pub struct MonolithicSpeciationBlock {
    config: SpeciationConfig,
}

impl MonolithicSpeciationBlock {
    pub fn new(config: SpeciationConfig) -> Result<Self> {
        require(
            config.max_iterations > 0
                && config.absolute_tolerance > 0.0
                && config.relative_tolerance > 0.0,
            "invalid speciation controls",
        )?;
        Ok(Self { config })
    }
    pub fn solve(
        &self,
        total_concentration: f64,
        log_equilibrium_constants: &[f64],
        initial_concentrations: Option<&[f64]>,
    ) -> Result<SpeciationResult> {
        require(
            total_concentration > 0.0
                && total_concentration.is_finite()
                && !log_equilibrium_constants.is_empty(),
            "invalid speciation total or reaction block",
        )?;
        let n = log_equilibrium_constants.len() + 1;
        let mut x = vec![(total_concentration / n as f64).ln(); n];
        if let Some(init) = initial_concentrations {
            require(init.len() == n, "speciation initial state size mismatch")?;
            require(
                init.iter().all(|c| *c > 0.0 && c.is_finite()),
                "speciation initial concentration must be positive",
            )?;
            for (xi, c) in x.iter_mut().zip(init) {
                *xi = c.ln();
            }
        }
        let k = log_equilibrium_constants;
        let residual = |x: &[f64], out: &mut [f64]| -> Result<()> {
            let total: f64 = x.iter().map(|v| v.exp()).sum();
            for i in 1..n {
                out[i - 1] = x[i] - x[0] - k[i - 1];
            }
            out[n - 1] = total - total_concentration;
            Ok(())
        };
        let jacobian = |x: &[f64], d: &[f64], out: &mut [f64]| -> Result<()> {
            for i in 1..n {
                out[i - 1] = d[i] - d[0];
            }
            out[n - 1] = (0..n).map(|i| x[i].exp() * d[i]).sum();
            Ok(())
        };
        let system = MatrixFreeNonlinearSystem {
            size: n,
            residual: &residual,
            jacobian_action: Some(&jacobian),
            transpose_jacobian_action: None,
        };
        let c = MatrixFreeNewtonConfig {
            max_iterations: self.config.max_iterations,
            linear_max_iterations: 32,
            linear_restart: n.min(8),
            absolute_tolerance: self.config.absolute_tolerance,
            relative_tolerance: self.config.relative_tolerance,
            linear_relative_tolerance: 1e-12,
            ..Default::default()
        };
        let r = solve_matrix_free_newton(&mut x, &system, c)?;
        Ok(SpeciationResult {
            concentrations: x.iter().map(|v| v.exp()).collect(),
            iterations: r.iterations,
            residual_rms: r.residual_rms,
            converged: r.converged,
        })
    }
}

/// NASA 7-coefficient polynomial (single range) in SI molar units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Nasa7Polynomial {
    pub a: [f64; 7],
    pub minimum_temperature: f64,
    pub maximum_temperature: f64,
}

impl Default for Nasa7Polynomial {
    fn default() -> Self {
        Self {
            a: [0.0; 7],
            minimum_temperature: 200.0,
            maximum_temperature: 6000.0,
        }
    }
}

impl Nasa7Polynomial {
    fn check(&self, t: f64) -> Result<()> {
        require(
            t > 0.0
                && t.is_finite()
                && t >= self.minimum_temperature
                && t <= self.maximum_temperature,
            "temperature outside NASA polynomial range",
        )
    }
    pub fn cp_molar(&self, t: f64) -> Result<f64> {
        self.check(t)?;
        let a = self.a;
        Ok(GAS_CONSTANT
            * (a[0] + a[1] * t + a[2] * t * t + a[3] * t * t * t + a[4] * t * t * t * t))
    }
    pub fn h_molar(&self, t: f64) -> Result<f64> {
        self.check(t)?;
        let a = self.a;
        Ok(GAS_CONSTANT
            * t
            * (a[0]
                + a[1] * t / 2.0
                + a[2] * t * t / 3.0
                + a[3] * t * t * t / 4.0
                + a[4] * t * t * t * t / 5.0
                + a[5] / t))
    }
    pub fn s_molar(&self, t: f64) -> Result<f64> {
        self.check(t)?;
        let a = self.a;
        Ok(GAS_CONSTANT
            * (a[0] * t.ln()
                + a[1] * t
                + a[2] * t * t / 2.0
                + a[3] * t * t * t / 3.0
                + a[4] * t * t * t * t / 4.0
                + a[6]))
    }
}

/// NIST Shomate form with `t = T/1000`: cp, s in J/(mol K), h in kJ/mol.
#[derive(Clone, Copy, Debug, PartialEq)]
#[allow(non_snake_case)]
pub struct ShomatePolynomial {
    pub A: f64,
    pub B: f64,
    pub C: f64,
    pub D: f64,
    pub E: f64,
    pub F: f64,
    pub G: f64,
    pub H: f64,
    pub minimum_temperature: f64,
    pub maximum_temperature: f64,
}

impl Default for ShomatePolynomial {
    fn default() -> Self {
        Self {
            A: 0.0,
            B: 0.0,
            C: 0.0,
            D: 0.0,
            E: 0.0,
            F: 0.0,
            G: 0.0,
            H: 0.0,
            minimum_temperature: 200.0,
            maximum_temperature: 6000.0,
        }
    }
}

impl ShomatePolynomial {
    fn t(&self, temperature: f64) -> Result<f64> {
        require(
            temperature >= self.minimum_temperature && temperature <= self.maximum_temperature,
            "Shomate temperature out of range",
        )?;
        Ok(temperature / 1000.0)
    }
    pub fn cp_molar(&self, temperature: f64) -> Result<f64> {
        let t = self.t(temperature)?;
        Ok(self.A + self.B * t + self.C * t * t + self.D * t * t * t + self.E / (t * t))
    }
    pub fn h_molar_kj(&self, temperature: f64) -> Result<f64> {
        let t = self.t(temperature)?;
        Ok(self.A * t
            + self.B * t * t / 2.0
            + self.C * t * t * t / 3.0
            + self.D * t * t * t * t / 4.0
            - self.E / t
            + self.F
            - self.H)
    }
    pub fn s_molar(&self, temperature: f64) -> Result<f64> {
        let t = self.t(temperature)?;
        Ok(
            self.A * t.ln() + self.B * t + self.C * t * t / 2.0 + self.D * t * t * t / 3.0
                - self.E / (2.0 * t * t)
                + self.G,
        )
    }
}

/// Single 1:1 electrolyte Pitzer mean activity coefficient (A_phi = 0.392, b = 1.2, alpha = 2).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BinaryPitzer11 {
    pub beta0: f64,
    pub beta1: f64,
    pub c_phi: f64,
}

impl BinaryPitzer11 {
    pub fn mean_activity_coefficient(&self, molality: f64, temperature: f64) -> Result<f64> {
        require(molality >= 0.0 && temperature > 0.0, "invalid Pitzer state")?;
        if molality == 0.0 {
            return Ok(1.0);
        }
        let m = molality;
        let s = m.sqrt();
        let (a_phi, b) = (0.392, 1.2);
        let x = 2.0 * s;
        let g = if x < 1e-8 {
            1.0
        } else {
            2.0 * (1.0 - (1.0 + x) * (-x).exp()) / (x * x)
        };
        let f = -a_phi * (s / (1.0 + b * s) + (2.0 / b) * (1.0 + b * s).ln());
        Ok((f + m * (2.0 * self.beta0 + 2.0 * self.beta1 * g) + 1.5 * m * m * self.c_phi).exp())
    }
}

/// Activity coefficients for charged/neutral solutes given molalities.
pub trait ActivityModel {
    fn activity_coefficients(
        &self,
        charge: &[i32],
        molality: &[f64],
        temperature: f64,
    ) -> Result<Vec<f64>>;
}

fn ionic_strength(charge: &[i32], molality: &[f64]) -> Result<f64> {
    require(
        charge.len() == molality.len(),
        "activity-model size mismatch",
    )?;
    require(
        molality.iter().all(|m| *m >= 0.0 && m.is_finite()),
        "invalid molality",
    )?;
    Ok(charge
        .iter()
        .zip(molality)
        .map(|(z, m)| 0.5 * m * f64::from(*z) * f64::from(*z))
        .sum())
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PitzerInteraction {
    pub first: usize,
    pub second: usize,
    pub beta0: f64,
    pub beta1: f64,
    pub c_phi: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentPitzer {
    interactions: Vec<PitzerInteraction>,
    a_phi: f64,
    b: f64,
}

impl MulticomponentPitzer {
    /// Upstream defaults `a_phi = 0.392`, `b = 1.2`.
    pub fn new(interactions: Vec<PitzerInteraction>, a_phi: f64, b: f64) -> Self {
        Self {
            interactions,
            a_phi,
            b,
        }
    }
}

impl ActivityModel for MulticomponentPitzer {
    fn activity_coefficients(
        &self,
        charge: &[i32],
        molality: &[f64],
        temperature: f64,
    ) -> Result<Vec<f64>> {
        require(
            charge.len() == molality.len()
                && !charge.is_empty()
                && temperature > 0.0
                && self.a_phi > 0.0
                && self.b > 0.0,
            "invalid multicomponent Pitzer state",
        )?;
        require(
            molality.iter().all(|m| *m >= 0.0),
            "negative Pitzer molality",
        )?;
        let is: f64 = (0..charge.len())
            .map(|i| 0.5 * molality[i] * f64::from(charge[i] * charge[i]))
            .sum();
        let root = is.sqrt();
        let f = -self.a_phi
            * (root / (1.0 + self.b * root) + (2.0 / self.b) * (1.0 + self.b * root).ln());
        let mut out = vec![1.0; charge.len()];
        for (i, o) in out.iter_mut().enumerate() {
            let mut ln_gamma = f64::from(charge[i] * charge[i]) * f;
            for p in &self.interactions {
                require(
                    p.first < charge.len()
                        && p.second < charge.len()
                        && p.beta0 >= 0.0
                        && p.beta1 >= 0.0
                        && p.c_phi >= 0.0,
                    "invalid Pitzer interaction",
                )?;
                if i != p.first && i != p.second {
                    continue;
                }
                let partner = if i == p.first { p.second } else { p.first };
                let x = 2.0 * root;
                let g = if x < 1e-12 {
                    1.0
                } else {
                    2.0 * (1.0 - (1.0 + x) * (-x).exp()) / (x * x)
                };
                ln_gamma += molality[partner] * (2.0 * p.beta0 + 2.0 * p.beta1 * g)
                    + p.c_phi * molality[partner] * molality[partner];
            }
            *o = ln_gamma.exp();
        }
        Ok(out)
    }
}

/// Davies equation `log10 gamma = -A z^2 (sqrt I/(1 + sqrt I) - 0.3 I)` (upstream A = 0.509).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DaviesActivity {
    pub a: f64,
}

impl Default for DaviesActivity {
    fn default() -> Self {
        Self { a: 0.509 }
    }
}

impl ActivityModel for DaviesActivity {
    fn activity_coefficients(
        &self,
        charge: &[i32],
        molality: &[f64],
        temperature: f64,
    ) -> Result<Vec<f64>> {
        let i = ionic_strength(charge, molality)?;
        require(
            temperature > 0.0 && temperature.is_finite() && self.a > 0.0,
            "invalid Davies controls",
        )?;
        let root = i.sqrt();
        Ok(charge
            .iter()
            .map(|z| {
                10f64
                    .powf(-self.a * f64::from(*z) * f64::from(*z) * (root / (1.0 + root) - 0.3 * i))
            })
            .collect())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SitPair {
    pub first: usize,
    pub second: usize,
    /// kg/mol.
    pub epsilon: f64,
}

/// Specific ion interaction theory with a 1.5 sqrt(I) Debye-Huckel denominator.
#[derive(Clone, Debug, PartialEq)]
pub struct SitActivity {
    interactions: Vec<SitPair>,
    a: f64,
}

impl SitActivity {
    /// Upstream default `a = 0.509`.
    pub fn new(interactions: Vec<SitPair>, a_parameter: f64) -> Self {
        Self {
            interactions,
            a: a_parameter,
        }
    }
}

impl ActivityModel for SitActivity {
    fn activity_coefficients(
        &self,
        charge: &[i32],
        molality: &[f64],
        temperature: f64,
    ) -> Result<Vec<f64>> {
        let i = ionic_strength(charge, molality)?;
        require(
            temperature > 0.0 && temperature.is_finite() && self.a > 0.0,
            "invalid SIT controls",
        )?;
        let root = i.sqrt();
        let mut lg: Vec<f64> = charge
            .iter()
            .map(|z| -self.a * f64::from(*z) * f64::from(*z) * root / (1.0 + 1.5 * root))
            .collect();
        for p in &self.interactions {
            require(
                p.first < molality.len() && p.second < molality.len() && p.epsilon.is_finite(),
                "invalid SIT interaction",
            )?;
            lg[p.first] += p.epsilon * molality[p.second];
            if p.first != p.second {
                lg[p.second] += p.epsilon * molality[p.first];
            }
        }
        Ok(lg.iter().map(|v| 10f64.powf(*v)).collect())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BinaryIdealEquilibrium {
    pub mole_fraction_a: f64,
    pub mole_fraction_b: f64,
    pub gibbs_molar: f64,
}

/// Closed-form ideal-mixing minimum `x_b = 1/(1 + exp((mu_b - mu_a)/RT))`.
pub fn minimize_binary_ideal_gibbs(
    standard_mu_a: f64,
    standard_mu_b: f64,
    temperature: f64,
) -> Result<BinaryIdealEquilibrium> {
    let (mu_a, mu_b, t) = (standard_mu_a, standard_mu_b, temperature);
    require(
        t > 0.0 && t.is_finite() && mu_a.is_finite() && mu_b.is_finite(),
        "invalid binary Gibbs equilibrium input",
    )?;
    let delta = (mu_b - mu_a) / (GAS_CONSTANT * t);
    let xb = if delta > 50.0 {
        (-delta).exp()
    } else if delta < -50.0 {
        1.0 - delta.exp()
    } else {
        1.0 / (1.0 + delta.exp())
    }
    .clamp(1e-15, 1.0 - 1e-15);
    let xa = 1.0 - xb;
    let g = xa * mu_a + xb * mu_b + GAS_CONSTANT * t * (xa * xa.ln() + xb * xb.ln());
    Ok(BinaryIdealEquilibrium {
        mole_fraction_a: xa,
        mole_fraction_b: xb,
        gibbs_molar: g,
    })
}

pub fn langmuir_surface_coverage(equilibrium_constant: f64, activity: f64) -> Result<f64> {
    require(
        equilibrium_constant >= 0.0
            && activity >= 0.0
            && equilibrium_constant.is_finite()
            && activity.is_finite(),
        "invalid Langmuir controls",
    )?;
    let ka = equilibrium_constant * activity;
    Ok(ka / (1.0 + ka))
}

/// Exchanger fraction of B for `K = (x_B a_A)/(x_A a_B)`.
pub fn binary_ion_exchange_fraction(
    selectivity: f64,
    activity_a: f64,
    activity_b: f64,
) -> Result<f64> {
    require(
        selectivity > 0.0
            && activity_a >= 0.0
            && activity_b >= 0.0
            && selectivity.is_finite()
            && activity_a.is_finite()
            && activity_b.is_finite(),
        "invalid ion-exchange controls",
    )?;
    if activity_b == 0.0 {
        return Ok(0.0);
    }
    if activity_a == 0.0 {
        return Ok(1.0);
    }
    let ratio = selectivity * activity_b / activity_a;
    Ok(ratio / (1.0 + ratio))
}

/// Character cursor mimicking `std::istream >>` for whitespace-separated tokens and doubles.
struct Cursor<'a> {
    s: &'a str,
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(s: &'a str) -> Self {
        Self { s, pos: 0 }
    }
    fn skip_ws(&mut self) {
        while self.pos < self.s.len() && self.s.as_bytes()[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }
    fn token(&mut self) -> Option<&'a str> {
        self.skip_ws();
        let start = self.pos;
        while self.pos < self.s.len() && !self.s.as_bytes()[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
        (self.pos > start).then(|| &self.s[start..self.pos])
    }
    fn number(&mut self) -> Option<f64> {
        self.skip_ws();
        let (v, used) = parse_float_prefix(&self.s[self.pos..])?;
        self.pos += used;
        Some(v)
    }
}

fn clean_line(line: &str) -> &str {
    line.split('#')
        .next()
        .unwrap_or("")
        .trim_matches([' ', '\t', '\r'])
}

/// PHREEQC `SOLUTION_MASTER_SPECIES` rows `element name molar_mass charge`; masses above one
/// are taken as g/mol and returned in kg/mol.
pub fn load_phreeqc_master_species(database_text: &str) -> Result<Vec<Species>> {
    let mut species = Vec::new();
    let mut in_section = false;
    for line in database_text.lines() {
        let clean = clean_line(line);
        if clean.is_empty() {
            continue;
        }
        if clean.to_lowercase() == "solution_master_species" {
            in_section = true;
            continue;
        }
        if in_section && clean.as_bytes()[0].is_ascii_uppercase() && !clean.contains([' ', '\t']) {
            in_section = false;
            continue;
        }
        if !in_section {
            continue;
        }
        let mut row = Cursor::new(clean);
        let (Some(_element), Some(name)) = (row.token(), row.token()) else {
            continue;
        };
        let (Some(mut molar_mass), Some(charge)) = (row.number(), row.number()) else {
            continue;
        };
        require(
            molar_mass > 0.0 && molar_mass.is_finite() && charge.is_finite(),
            "invalid PHREEQC master species row",
        )?;
        if molar_mass > 1.0 {
            molar_mass *= 1e-3;
        }
        species.push(Species {
            name: name.into(),
            molar_mass,
            charge: charge.round() as i32,
        });
    }
    require(
        !species.is_empty(),
        "PHREEQC database has no SOLUTION_MASTER_SPECIES rows",
    )?;
    Ok(species)
}

fn numeric_after(line: &str, key: &str, fallback: f64) -> f64 {
    let lower = line.to_ascii_lowercase();
    let Some(found) = lower.find(&key.to_ascii_lowercase()) else {
        return fallback;
    };
    let begin = line[found..]
        .find(':')
        .map_or(found + key.len(), |c| found + c + 1);
    let mut c = Cursor::new(&line[begin.min(line.len())..]);
    c.number().unwrap_or(fallback)
}

fn find_species(indices: &HashMap<String, usize>, name: &str) -> Result<usize> {
    let name = name
        .trim_matches([' ', '\t', '\r'])
        .trim_start_matches(['+', '-']);
    indices.get(name).copied().ok_or_else(|| {
        SolverError::InvalidInput(format!(
            "Cantera reaction references unknown species: {name}"
        ))
    })
}

fn parse_side(side: &str, indices: &HashMap<String, usize>) -> Result<Vec<StoichiometricTerm>> {
    let mut terms = Vec::new();
    for raw in side.split('+') {
        let term = raw.trim_matches([' ', '\t', '\r']);
        if term.is_empty() {
            continue;
        }
        let mut c = Cursor::new(term);
        let Some(mut name) = c.token() else {
            return Err(SolverError::InvalidInput(
                "empty Cantera stoichiometric term".into(),
            ));
        };
        let mut coefficient = 1.0;
        if let Some((parsed, used)) = parse_float_prefix(name) {
            coefficient = parsed;
            if used == name.len() {
                name = c.token().ok_or_else(|| {
                    SolverError::InvalidInput("missing Cantera species name".into())
                })?;
            } else {
                name = &name[used..];
            }
        }
        require(
            coefficient > 0.0 && coefficient.is_finite(),
            "invalid Cantera stoichiometric coefficient",
        )?;
        terms.push(StoichiometricTerm {
            species: find_species(indices, name)?,
            coefficient,
        });
    }
    Ok(terms)
}

/// Cantera YAML subset: `species:` names and `reactions:` equations with Arrhenius `A`, `b`, `Ea`.
pub fn load_cantera_yaml(mechanism_text: &str) -> Result<ReactionNetwork> {
    let mut network = ReactionNetwork::new(Vec::new());
    let mut indices: HashMap<String, usize> = HashMap::new();
    let (mut in_species, mut in_reactions) = (false, false);
    let mut equation = String::new();
    let (mut a, mut b, mut ea) = (0.0, 0.0, 0.0);
    let flush = |network: &mut ReactionNetwork,
                 equation: &mut String,
                 a: &mut f64,
                 b: &mut f64,
                 ea: &mut f64,
                 indices: &HashMap<String, usize>|
     -> Result<()> {
        if equation.is_empty() {
            return Ok(());
        }
        let (arrow, width) = match equation.find("<=>") {
            Some(p) => (p, 3),
            None => (
                equation.find("=>").ok_or_else(|| {
                    SolverError::InvalidInput("Cantera reaction is missing an arrow".into())
                })?,
                2,
            ),
        };
        let reaction = ElementaryReaction {
            reactants: parse_side(&equation[..arrow], indices)?,
            products: parse_side(&equation[arrow + width..], indices)?,
            forward: ArrheniusRate {
                pre_exponential: *a,
                temperature_exponent: *b,
                activation_energy: *ea,
            },
            ..Default::default()
        };
        network.add_reaction(reaction)?;
        equation.clear();
        (*a, *b, *ea) = (0.0, 0.0, 0.0);
        Ok(())
    };
    for line in mechanism_text.lines() {
        let clean = clean_line(line);
        if clean.is_empty() {
            continue;
        }
        let folded = clean.to_lowercase();
        if folded == "species:" {
            (in_species, in_reactions) = (true, false);
            continue;
        }
        if folded == "reactions:" {
            (in_species, in_reactions) = (false, true);
            continue;
        }
        if in_species && folded.starts_with("- name:") {
            let name =
                clean[clean.find(':').map_or(0, |p| p + 1)..].trim_matches([' ', '\t', '\r']);
            require(!name.is_empty(), "Cantera species name is empty")?;
            require(
                !indices.contains_key(name),
                format!("duplicate Cantera species: {name}"),
            )?;
            indices.insert(name.into(), indices.len());
            network.add_species(Species {
                name: name.into(),
                molar_mass: 1.0,
                charge: 0,
            });
            continue;
        }
        if in_reactions && folded.starts_with("- equation:") {
            flush(
                &mut network,
                &mut equation,
                &mut a,
                &mut b,
                &mut ea,
                &indices,
            )?;
            equation = clean[clean.find(':').map_or(0, |p| p + 1)..]
                .trim_matches([' ', '\t', '\r'])
                .into();
            continue;
        }
        if in_reactions && !equation.is_empty() {
            a = numeric_after(clean, "a:", a);
            b = numeric_after(clean, "b:", b);
            ea = numeric_after(clean, "ea:", ea);
        }
    }
    flush(
        &mut network,
        &mut equation,
        &mut a,
        &mut b,
        &mut ea,
        &indices,
    )?;
    require(!indices.is_empty(), "Cantera YAML has no species section")?;
    Ok(network)
}

/// Optional external-equilibrium seam: an application supplies the evaluator (e.g. its own
/// Reaktoro binding); results are validated for size and non-negativity.
/// `(composition, temperature, pressure) -> equilibrium composition`.
pub type ReaktoroEvaluator = Box<dyn Fn(&[f64], f64, f64) -> Vec<f64>>;

pub struct ReaktoroDynamicAdapter {
    evaluator: ReaktoroEvaluator,
}

impl ReaktoroDynamicAdapter {
    pub fn new(evaluator: ReaktoroEvaluator) -> Self {
        Self { evaluator }
    }
    pub fn equilibrate(
        &self,
        composition: &[f64],
        temperature: f64,
        pressure: f64,
    ) -> Result<Vec<f64>> {
        require(
            temperature > 0.0 && pressure > 0.0 && temperature.is_finite() && pressure.is_finite(),
            "invalid Reaktoro adapter state",
        )?;
        let r = (self.evaluator)(composition, temperature, pressure);
        require(
            r.len() == composition.len(),
            "Reaktoro adapter size mismatch",
        )?;
        require(
            r.iter().all(|v| *v >= 0.0 && v.is_finite()),
            "Reaktoro adapter returned invalid composition",
        )?;
        Ok(r)
    }
}
