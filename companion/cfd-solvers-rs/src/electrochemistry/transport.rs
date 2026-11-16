// Stencil loops index several per-cell arrays together, mirroring upstream.
#![allow(clippy::needless_range_loop)]
use super::{butler_volmer_current_density, FARADAY_CONSTANT};
use crate::{
    chemistry::GAS_CONSTANT,
    core::{matrix_free_conjugate_gradient, ConjugateGradientResult, Vec3},
    fvm::{decompose_face_area, PolyMesh},
    require, Result, SolverError,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransportBoundary1D {
    Periodic,
    NoFlux,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NernstPlanck1DConfig {
    pub cells: usize,
    pub length: f64,
    pub temperature: f64,
    pub dt: f64,
    pub advection_velocity: f64,
    pub boundary: TransportBoundary1D,
}

impl Default for NernstPlanck1DConfig {
    fn default() -> Self {
        Self {
            cells: 128,
            length: 1.0,
            temperature: 298.15,
            dt: 1e-4,
            advection_velocity: 0.0,
            boundary: TransportBoundary1D::Periodic,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct NernstPlanckSpecies1D {
    pub name: String,
    pub charge: i32,
    pub diffusivity: f64,
    pub concentration: Vec<f64>,
}

/// Explicit conservative 1-D Nernst-Planck with upwind advection and a prescribed potential.
#[derive(Clone, Debug)]
pub struct NernstPlanck1D {
    config: NernstPlanck1DConfig,
    dx: f64,
    time: f64,
    potential: Vec<f64>,
    species: Vec<NernstPlanckSpecies1D>,
}

impl NernstPlanck1D {
    pub fn new(config: NernstPlanck1DConfig) -> Result<Self> {
        require(
            config.cells >= 3 && config.length > 0.0 && config.temperature > 0.0 && config.dt > 0.0,
            "invalid Nernst-Planck 1-D configuration",
        )?;
        Ok(Self {
            config,
            dx: config.length / config.cells as f64,
            time: 0.0,
            potential: vec![0.0; config.cells],
            species: Vec::new(),
        })
    }
    pub fn config(&self) -> &NernstPlanck1DConfig {
        &self.config
    }
    pub fn potential(&self) -> &[f64] {
        &self.potential
    }
    pub fn species(&self) -> &[NernstPlanckSpecies1D] {
        &self.species
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn dx(&self) -> f64 {
        self.dx
    }
    pub fn set_potential_fn(&mut self, f: impl Fn(f64) -> f64) {
        for (i, p) in self.potential.iter_mut().enumerate() {
            *p = f((i as f64 + 0.5) * self.dx);
        }
    }
    pub fn set_potential(&mut self, values: &[f64]) -> Result<()> {
        require(values.len() == self.config.cells, "potential size mismatch")?;
        self.potential.copy_from_slice(values);
        Ok(())
    }
    pub fn add_species_fn(
        &mut self,
        name: &str,
        charge: i32,
        diffusivity: f64,
        initial_concentration: impl Fn(f64) -> f64,
    ) -> Result<usize> {
        require(diffusivity > 0.0, "species diffusivity must be positive")?;
        let concentration: Vec<f64> = (0..self.config.cells)
            .map(|i| initial_concentration((i as f64 + 0.5) * self.dx))
            .collect();
        require(
            concentration.iter().all(|c| *c >= 0.0 && c.is_finite()),
            "invalid initial concentration",
        )?;
        self.species.push(NernstPlanckSpecies1D {
            name: name.into(),
            charge,
            diffusivity,
            concentration,
        });
        Ok(self.species.len() - 1)
    }
    pub fn add_species(
        &mut self,
        name: &str,
        charge: i32,
        diffusivity: f64,
        uniform_concentration: f64,
    ) -> Result<usize> {
        require(
            uniform_concentration >= 0.0,
            "concentration must be non-negative",
        )?;
        self.add_species_fn(name, charge, diffusivity, |_| uniform_concentration)
    }
    pub fn step(&mut self, count: usize) -> Result<()> {
        require(
            !self.species.is_empty(),
            "Nernst-Planck solver has no species",
        )?;
        for _ in 0..count {
            for s in 0..self.species.len() {
                self.step_species(s)?;
            }
            self.time += self.config.dt;
        }
        Ok(())
    }
    fn step_species(&mut self, index: usize) -> Result<()> {
        let n = self.config.cells;
        let (dx, c) = (self.dx, self.config);
        let sp = &self.species[index];
        let scale = f64::from(sp.charge) * FARADAY_CONSTANT / (GAS_CONSTANT * c.temperature);
        let face = |l: usize, r: usize| {
            let (cl, cr) = (sp.concentration[l], sp.concentration[r]);
            let adv = if c.advection_velocity >= 0.0 {
                c.advection_velocity * cl
            } else {
                c.advection_velocity * cr
            };
            adv - sp.diffusivity
                * ((cr - cl) / dx
                    + scale * 0.5 * (cl + cr) * (self.potential[r] - self.potential[l]) / dx)
        };
        let mut flux = vec![0.0; n + 1];
        for f in 1..n {
            flux[f] = face(f - 1, f);
        }
        if c.boundary == TransportBoundary1D::Periodic {
            flux[0] = face(n - 1, 0);
            flux[n] = flux[0];
        }
        let mut next = vec![0.0; n];
        for i in 0..n {
            let v = sp.concentration[i] - (c.dt / dx) * (flux[i + 1] - flux[i]);
            if !v.is_finite() || v < -1e-12 {
                return Err(SolverError::NumericalFailure(
                    "Nernst-Planck step became unstable or negative; reduce dt".into(),
                ));
            }
            next[i] = v.max(0.0);
        }
        self.species[index].concentration = next;
        Ok(())
    }
    pub fn total_amount(&self, species: usize) -> Result<f64> {
        let s = self
            .species
            .get(species)
            .ok_or_else(|| SolverError::InvalidInput("species index out of range".into()))?;
        Ok(s.concentration.iter().sum::<f64>() * self.dx)
    }
    pub fn minimum_concentration(&self, species: usize) -> Result<f64> {
        let s = self
            .species
            .get(species)
            .ok_or_else(|| SolverError::InvalidInput("species index out of range".into()))?;
        Ok(s.concentration
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum SpeciesBoundaryType {
    #[default]
    NoFlux,
    FixedConcentration,
    FixedMolarFlux,
    ButlerVolmerFlux,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum PotentialBoundaryType {
    #[default]
    Insulating,
    FixedPotential,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PotentialBoundaryCondition {
    pub kind: PotentialBoundaryType,
    pub value: f64,
}

/// `value` is a concentration (mol/m3) or outward molar flux (mol/m2/s); Butler-Volmer fields
/// describe an anodic reaction producing `stoichiometric` moles per event into the electrolyte.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeciesBoundaryCondition {
    pub kind: SpeciesBoundaryType,
    pub value: f64,
    pub electrode_potential: f64,
    pub equilibrium_potential: f64,
    pub exchange_current_density: f64,
    pub alpha_anodic: f64,
    pub alpha_cathodic: f64,
    pub electrons: i32,
    pub stoichiometric: f64,
}

impl Default for SpeciesBoundaryCondition {
    fn default() -> Self {
        Self {
            kind: SpeciesBoundaryType::NoFlux,
            value: 0.0,
            electrode_potential: 0.0,
            equilibrium_potential: 0.0,
            exchange_current_density: 0.0,
            alpha_anodic: 0.5,
            alpha_cathodic: 0.5,
            electrons: 1,
            stoichiometric: 1.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PolyTransportSpecies {
    pub name: String,
    pub charge: i32,
    pub diffusivity: f64,
    pub concentration: Vec<f64>,
    pub boundary: Vec<SpeciesBoundaryCondition>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ElectromigrationFluxScheme {
    #[default]
    Centered,
    ScharfetterGummel,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NernstPlanckPolyConfig {
    pub temperature: f64,
    pub dt: f64,
    pub relative_permittivity: f64,
    pub electromigration_scheme: ElectromigrationFluxScheme,
}

impl Default for NernstPlanckPolyConfig {
    fn default() -> Self {
        Self {
            temperature: 298.15,
            dt: 1e-5,
            relative_permittivity: 78.5,
            electromigration_scheme: ElectromigrationFluxScheme::Centered,
        }
    }
}

/// Bernoulli function `x / (exp(x) - 1)` with series and asymptotic branches.
pub fn bernoulli_exponential_fit(x: f64) -> f64 {
    if x.abs() < 1e-6 {
        1.0 - 0.5 * x + x * x / 12.0
    } else if x > 50.0 {
        x * (-x).exp()
    } else if x < -50.0 {
        -x
    } else {
        x / x.exp_m1()
    }
}

/// Explicit conservative finite-volume Nernst-Planck on a PolyMesh with prescribed or
/// electroneutral/Poisson potentials and patch boundary conditions.
#[derive(Clone, Debug)]
pub struct NernstPlanckPolyMesh {
    mesh: PolyMesh,
    config: NernstPlanckPolyConfig,
    species: Vec<PolyTransportSpecies>,
    potential: Vec<f64>,
    potential_boundary: Vec<PotentialBoundaryCondition>,
    face_flux: Vec<f64>,
    potential_result: ConjugateGradientResult,
    time: f64,
}

impl NernstPlanckPolyMesh {
    pub fn new(mesh: PolyMesh, config: NernstPlanckPolyConfig) -> Result<Self> {
        require(
            config.temperature > 0.0 && config.dt > 0.0 && config.relative_permittivity > 0.0,
            "invalid PolyMesh Nernst-Planck controls",
        )?;
        let (n, p, f) = (mesh.cell_count(), mesh.patches().len(), mesh.face_count());
        Ok(Self {
            mesh,
            config,
            species: Vec::new(),
            potential: vec![0.0; n],
            potential_boundary: vec![PotentialBoundaryCondition::default(); p],
            face_flux: vec![0.0; f],
            potential_result: ConjugateGradientResult::default(),
            time: 0.0,
        })
    }
    pub fn add_species_fn(
        &mut self,
        name: &str,
        charge: i32,
        diffusivity: f64,
        initial_concentration: impl Fn(Vec3) -> f64,
    ) -> Result<usize> {
        require(
            !name.is_empty() && diffusivity >= 0.0,
            "invalid transport species",
        )?;
        let concentration: Vec<f64> = self
            .mesh
            .cells()
            .iter()
            .map(|c| initial_concentration(c.center))
            .collect();
        require(
            concentration.iter().all(|c| *c >= 0.0 && c.is_finite()),
            "initial species concentration must be finite and non-negative",
        )?;
        self.species.push(PolyTransportSpecies {
            name: name.into(),
            charge,
            diffusivity,
            concentration,
            boundary: vec![SpeciesBoundaryCondition::default(); self.mesh.patches().len()],
        });
        Ok(self.species.len() - 1)
    }
    pub fn add_species(
        &mut self,
        name: &str,
        charge: i32,
        diffusivity: f64,
        initial_concentration: f64,
    ) -> Result<usize> {
        self.add_species_fn(name, charge, diffusivity, |_| initial_concentration)
    }
    fn species_mut(&mut self, species: usize) -> Result<&mut PolyTransportSpecies> {
        self.species
            .get_mut(species)
            .ok_or_else(|| SolverError::InvalidInput("species index out of range".into()))
    }
    pub fn set_species_boundary(
        &mut self,
        species: usize,
        patch: &str,
        kind: SpeciesBoundaryType,
        value: f64,
    ) -> Result<()> {
        require(
            !(kind == SpeciesBoundaryType::FixedConcentration
                && !(value >= 0.0 && value.is_finite())),
            "fixed concentration must be finite and non-negative",
        )?;
        require(value.is_finite(), "species boundary value must be finite")?;
        let p = self.mesh.patch_index(patch)?;
        self.species_mut(species)?.boundary[p] = SpeciesBoundaryCondition {
            kind,
            value,
            ..Default::default()
        };
        Ok(())
    }
    /// Upstream defaults: stoichiometric 1, alpha_a = alpha_c = 0.5.
    #[allow(clippy::too_many_arguments)]
    pub fn set_species_butler_volmer_boundary(
        &mut self,
        species: usize,
        patch: &str,
        electrode_potential: f64,
        equilibrium_potential: f64,
        exchange_current_density: f64,
        electrons: i32,
        stoichiometric: f64,
        alpha_anodic: f64,
        alpha_cathodic: f64,
    ) -> Result<()> {
        require(
            exchange_current_density >= 0.0
                && exchange_current_density.is_finite()
                && electrons > 0
                && electrode_potential.is_finite()
                && equilibrium_potential.is_finite()
                && stoichiometric.is_finite()
                && alpha_anodic > 0.0
                && alpha_cathodic > 0.0,
            "invalid Butler-Volmer boundary parameters",
        )?;
        let p = self.mesh.patch_index(patch)?;
        self.species_mut(species)?.boundary[p] = SpeciesBoundaryCondition {
            kind: SpeciesBoundaryType::ButlerVolmerFlux,
            value: 0.0,
            electrode_potential,
            equilibrium_potential,
            exchange_current_density,
            alpha_anodic,
            alpha_cathodic,
            electrons,
            stoichiometric,
        };
        Ok(())
    }
    pub fn set_potential(&mut self, value: f64) -> Result<()> {
        require(value.is_finite(), "potential must be finite")?;
        self.potential.fill(value);
        Ok(())
    }
    pub fn set_potential_fn(&mut self, value: impl Fn(Vec3) -> f64) -> Result<()> {
        let p: Vec<f64> = self.mesh.cells().iter().map(|c| value(c.center)).collect();
        require(p.iter().all(|v| v.is_finite()), "potential must be finite")?;
        self.potential = p;
        Ok(())
    }
    pub fn set_potential_boundary(
        &mut self,
        patch: &str,
        kind: PotentialBoundaryType,
        value: f64,
    ) -> Result<()> {
        require(value.is_finite(), "potential boundary value must be finite")?;
        let p = self.mesh.patch_index(patch)?;
        self.potential_boundary[p] = PotentialBoundaryCondition { kind, value };
        Ok(())
    }
    pub fn set_face_flux(&mut self, volumetric_flux: Vec<f64>) -> Result<()> {
        require(
            volumetric_flux.len() == self.mesh.face_count(),
            "Nernst-Planck face-flux size mismatch",
        )?;
        self.face_flux = volumetric_flux;
        Ok(())
    }
    pub fn step(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step_once()?;
        }
        Ok(())
    }
    /// Upstream defaults 1000 Poisson iterations, relative tolerance 1e-10.
    pub fn step_poisson_nernst_planck(
        &mut self,
        steps: usize,
        poisson_max_iterations: usize,
        poisson_relative_tolerance: f64,
    ) -> Result<()> {
        for _ in 0..steps {
            self.solve_poisson_potential(poisson_max_iterations, poisson_relative_tolerance)?;
            self.step_once()?;
        }
        Ok(())
    }
    fn patch_of(&self, fi: usize) -> Result<usize> {
        self.mesh.faces()[fi]
            .patch
            .ok_or_else(|| SolverError::InvalidInput("boundary face without patch".into()))
    }
    fn step_once(&mut self) -> Result<()> {
        let c = self.config;
        let scale = FARADAY_CONSTANT / (GAS_CONSTANT * c.temperature);
        let sg = c.electromigration_scheme == ElectromigrationFluxScheme::ScharfetterGummel;
        let mut updated = Vec::with_capacity(self.species.len());
        for sp in &self.species {
            let mut flux = vec![0.0; self.mesh.face_count()];
            for (fi, face) in self.mesh.faces().iter().enumerate() {
                let geometry = decompose_face_area(&self.mesh, fi)?;
                let co = sp.concentration[face.owner];
                let z = f64::from(sp.charge);
                let drift = |cb: f64, dphi: f64, adv: f64| {
                    if sg && sp.charge != 0 {
                        let psi = z * scale * dphi;
                        adv + sp.diffusivity
                            * geometry.orthogonal_metric
                            * (co * bernoulli_exponential_fit(psi)
                                - cb * bernoulli_exponential_fit(-psi))
                    } else {
                        adv - sp.diffusivity * geometry.orthogonal_metric * (cb - co)
                            - z * sp.diffusivity
                                * scale
                                * 0.5
                                * (co + cb)
                                * geometry.orthogonal_metric
                                * dphi
                    }
                };
                match face.neighbour {
                    Some(nb) => {
                        let cn = sp.concentration[nb];
                        let phi = self.face_flux[fi];
                        let cf = if phi >= 0.0 { co } else { cn };
                        flux[fi] = drift(
                            cn,
                            self.potential[nb] - self.potential[face.owner],
                            phi * cf,
                        );
                    }
                    None => {
                        let p = self.patch_of(fi)?;
                        let bc = sp.boundary[p];
                        let pbc = self.potential_boundary[p];
                        let area = face.area.norm();
                        flux[fi] = match bc.kind {
                            SpeciesBoundaryType::NoFlux => 0.0,
                            SpeciesBoundaryType::FixedMolarFlux => bc.value * area,
                            SpeciesBoundaryType::ButlerVolmerFlux => {
                                let electrolyte =
                                    if pbc.kind == PotentialBoundaryType::FixedPotential {
                                        pbc.value
                                    } else {
                                        self.potential[face.owner]
                                    };
                                let eta =
                                    bc.electrode_potential - electrolyte - bc.equilibrium_potential;
                                let current = butler_volmer_current_density(
                                    bc.exchange_current_density,
                                    eta,
                                    c.temperature,
                                    f64::from(bc.electrons),
                                    bc.alpha_anodic,
                                    bc.alpha_cathodic,
                                )?;
                                -bc.stoichiometric * current
                                    / (f64::from(bc.electrons) * FARADAY_CONSTANT)
                                    * area
                            }
                            SpeciesBoundaryType::FixedConcentration => {
                                let cb = bc.value;
                                let ff = self.face_flux[fi];
                                let adv = if ff >= 0.0 { ff * co } else { ff * cb };
                                let dphi = if pbc.kind == PotentialBoundaryType::FixedPotential {
                                    pbc.value - self.potential[face.owner]
                                } else {
                                    0.0
                                };
                                drift(cb, dphi, adv)
                            }
                        };
                    }
                }
            }
            let mut next = vec![0.0; self.mesh.cell_count()];
            for (cell, out) in next.iter_mut().enumerate() {
                let mut net = 0.0;
                for &fi in &self.mesh.cell_faces()[cell] {
                    net += if self.mesh.faces()[fi].owner == cell {
                        flux[fi]
                    } else {
                        -flux[fi]
                    };
                }
                let v = sp.concentration[cell] - c.dt * net / self.mesh.cells()[cell].volume;
                if !v.is_finite() || v < -1e-12 {
                    return Err(SolverError::NumericalFailure("Nernst-Planck explicit step produced a negative/invalid concentration; reduce dt".into()));
                }
                *out = v.max(0.0);
            }
            updated.push(next);
        }
        for (sp, next) in self.species.iter_mut().zip(updated) {
            sp.concentration = next;
        }
        self.time += c.dt;
        Ok(())
    }
    fn remove_mean(v: &mut [f64]) {
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        for x in v {
            *x -= mean;
        }
    }
    /// Current-free potential `div(F sum z D grad c + kappa grad phi) = 0` with insulating
    /// boundaries and a mean-zero gauge (upstream defaults 1000 iterations, 1e-10).
    pub fn solve_electroneutral_potential(
        &mut self,
        max_iterations: usize,
        relative_tolerance: f64,
    ) -> Result<ConjugateGradientResult> {
        require(
            !self.species.is_empty(),
            "electroneutral potential requires at least one ionic species",
        )?;
        require(
            max_iterations > 0 && relative_tolerance > 0.0,
            "invalid potential solve controls",
        )?;
        let kappa_scale =
            FARADAY_CONSTANT * FARADAY_CONSTANT / (GAS_CONSTANT * self.config.temperature);
        let nf = self.mesh.face_count();
        let (mut coefficient, mut diffusion_current) = (vec![0.0; nf], vec![0.0; nf]);
        for (fi, face) in self.mesh.faces().iter().enumerate() {
            let Some(nb) = face.neighbour else { continue };
            let g = decompose_face_area(&self.mesh, fi)?;
            let (mut k, mut d) = (0.0, 0.0);
            for sp in self
                .species
                .iter()
                .filter(|s| s.charge != 0 && s.diffusivity != 0.0)
            {
                let (co, cn) = (sp.concentration[face.owner], sp.concentration[nb]);
                let z = f64::from(sp.charge);
                k += z * z * sp.diffusivity * 0.5 * (co + cn);
                d += z * sp.diffusivity * (cn - co);
            }
            coefficient[fi] = kappa_scale * k * g.orthogonal_metric;
            diffusion_current[fi] = -FARADAY_CONSTANT * d * g.orthogonal_metric;
        }
        let mut rhs = vec![0.0; self.mesh.cell_count()];
        for (cell, r) in rhs.iter_mut().enumerate() {
            for &fi in &self.mesh.cell_faces()[cell] {
                let face = &self.mesh.faces()[fi];
                if face.neighbour.is_none() {
                    continue;
                }
                *r -= if face.owner == cell {
                    diffusion_current[fi]
                } else {
                    -diffusion_current[fi]
                };
            }
        }
        Self::remove_mean(&mut rhs);
        Self::remove_mean(&mut self.potential);
        let mesh = &self.mesh;
        let apply = |x: &[f64], out: &mut [f64]| {
            for (cell, o) in out.iter_mut().enumerate() {
                let mut sum = 0.0;
                for &fi in &mesh.cell_faces()[cell] {
                    let face = &mesh.faces()[fi];
                    let Some(nb) = face.neighbour else { continue };
                    let other = if face.owner == cell { nb } else { face.owner };
                    sum += coefficient[fi] * (x[cell] - x[other]);
                }
                *o = sum;
            }
        };
        self.potential_result = matrix_free_conjugate_gradient(
            &rhs,
            &mut self.potential,
            apply,
            max_iterations,
            relative_tolerance,
        )?;
        require(
            self.potential_result.converged,
            "electroneutral potential solve did not converge",
        )?;
        Self::remove_mean(&mut self.potential);
        Ok(self.potential_result)
    }
    /// `-div(eps grad phi) = F sum z c`, scaled by `1/eps`; fixed-potential patches are Dirichlet,
    /// a fully insulating problem needs net-zero charge and uses a mean-zero gauge.
    pub fn solve_poisson_potential(
        &mut self,
        max_iterations: usize,
        relative_tolerance: f64,
    ) -> Result<ConjugateGradientResult> {
        require(
            !self.species.is_empty(),
            "Poisson potential requires at least one charged species",
        )?;
        require(
            max_iterations > 0 && relative_tolerance > 0.0,
            "invalid Poisson solve controls",
        )?;
        let epsilon = 8.8541878128e-12 * self.config.relative_permittivity;
        let mut rhs: Vec<f64> = (0..self.mesh.cell_count())
            .map(|cell| {
                let q: f64 = self
                    .species
                    .iter()
                    .map(|s| f64::from(s.charge) * s.concentration[cell])
                    .sum();
                (FARADAY_CONSTANT / epsilon) * q * self.mesh.cells()[cell].volume
            })
            .collect();
        let has_dirichlet = self
            .potential_boundary
            .iter()
            .any(|b| b.kind == PotentialBoundaryType::FixedPotential);
        let metrics: Vec<f64> = (0..self.mesh.face_count())
            .map(|fi| decompose_face_area(&self.mesh, fi).map(|g| g.orthogonal_metric))
            .collect::<Result<_>>()?;
        if !has_dirichlet {
            let net: f64 = rhs.iter().sum();
            let scale: f64 = rhs.iter().map(|v| v.abs()).sum();
            require(
                net.abs() <= 1e-10 * scale.max(1.0),
                "pure-insulating Poisson problem requires net-zero charge",
            )?;
            let mean = net / rhs.len() as f64;
            for r in &mut rhs {
                *r -= mean;
            }
            Self::remove_mean(&mut self.potential);
        } else {
            for (fi, face) in self.mesh.faces().iter().enumerate() {
                if face.neighbour.is_some() {
                    continue;
                }
                let bc = self.potential_boundary[self.patch_of(fi)?];
                if bc.kind == PotentialBoundaryType::FixedPotential {
                    rhs[face.owner] += metrics[fi] * bc.value;
                }
            }
        }
        let fixed: Vec<bool> = self
            .mesh
            .faces()
            .iter()
            .map(|f| {
                f.neighbour.is_none()
                    && f.patch.is_some_and(|p| {
                        self.potential_boundary[p].kind == PotentialBoundaryType::FixedPotential
                    })
            })
            .collect();
        let mesh = &self.mesh;
        let apply = |x: &[f64], out: &mut [f64]| {
            for (cell, o) in out.iter_mut().enumerate() {
                let mut sum = 0.0;
                for &fi in &mesh.cell_faces()[cell] {
                    let face = &mesh.faces()[fi];
                    match face.neighbour {
                        Some(nb) => {
                            let other = if face.owner == cell { nb } else { face.owner };
                            sum += metrics[fi] * (x[cell] - x[other]);
                        }
                        None if fixed[fi] => sum += metrics[fi] * x[cell],
                        None => {}
                    }
                }
                *o = sum;
            }
        };
        self.potential_result = matrix_free_conjugate_gradient(
            &rhs,
            &mut self.potential,
            apply,
            max_iterations,
            relative_tolerance,
        )?;
        require(
            self.potential_result.converged,
            "Poisson potential solve did not converge",
        )?;
        if !has_dirichlet {
            Self::remove_mean(&mut self.potential);
        }
        Ok(self.potential_result)
    }
    /// Interior-face ionic current (owner to neighbour).
    pub fn current_flux_faces(&self) -> Result<Vec<f64>> {
        let scale = FARADAY_CONSTANT * FARADAY_CONSTANT / (GAS_CONSTANT * self.config.temperature);
        let mut current = vec![0.0; self.mesh.face_count()];
        for (fi, face) in self.mesh.faces().iter().enumerate() {
            let Some(nb) = face.neighbour else { continue };
            let g = decompose_face_area(&self.mesh, fi)?;
            let (mut d, mut k) = (0.0, 0.0);
            for sp in self
                .species
                .iter()
                .filter(|s| s.charge != 0 && s.diffusivity != 0.0)
            {
                let (co, cn) = (sp.concentration[face.owner], sp.concentration[nb]);
                let z = f64::from(sp.charge);
                d += z * sp.diffusivity * (cn - co);
                k += z * z * sp.diffusivity * 0.5 * (co + cn);
            }
            current[fi] = -FARADAY_CONSTANT * d * g.orthogonal_metric
                - scale
                    * k
                    * g.orthogonal_metric
                    * (self.potential[nb] - self.potential[face.owner]);
        }
        Ok(current)
    }
    pub fn total_charge(&self) -> f64 {
        FARADAY_CONSTANT
            * (0..self.mesh.cell_count())
                .map(|cell| {
                    self.species
                        .iter()
                        .map(|s| f64::from(s.charge) * s.concentration[cell])
                        .sum::<f64>()
                        * self.mesh.cells()[cell].volume
                })
                .sum::<f64>()
    }
    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn species(&self) -> &[PolyTransportSpecies] {
        &self.species
    }
    pub fn potential(&self) -> &[f64] {
        &self.potential
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn total_amount(&self, species: usize) -> Result<f64> {
        let s = self
            .species
            .get(species)
            .ok_or_else(|| SolverError::InvalidInput("species index out of range".into()))?;
        Ok(s.concentration
            .iter()
            .zip(self.mesh.cells())
            .map(|(c, cell)| c * cell.volume)
            .sum())
    }
    pub fn minimum_concentration(&self, species: usize) -> Result<f64> {
        let s = self
            .species
            .get(species)
            .ok_or_else(|| SolverError::InvalidInput("species index out of range".into()))?;
        Ok(s.concentration
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min))
    }
}
