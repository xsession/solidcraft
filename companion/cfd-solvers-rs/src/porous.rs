//! Periodic THMC slice, Allen-Cahn, fractured graph-FV two-phase THMC and Cahn-Hilliard with a
//! tabular CALPHAD free energy, ported from upstream `cfd/porous/{thmc1d, fractured_darcy_thmc,
//! cahn_hilliard_calphad}` (MIT, cfd_solvers contributors).
// Stencil loops index several per-cell arrays together, mirroring upstream.
#![allow(clippy::needless_range_loop)]
use crate::{require, Result};

fn left(i: usize, n: usize) -> usize {
    if i == 0 {
        n - 1
    } else {
        i - 1
    }
}
fn right(i: usize, n: usize) -> usize {
    if i + 1 == n {
        0
    } else {
        i + 1
    }
}
fn laplace(v: &[f64], i: usize, inv_h2: f64) -> f64 {
    let n = v.len();
    (v[left(i, n)] - 2.0 * v[i] + v[right(i, n)]) * inv_h2
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Thmc1DConfig {
    pub cells: usize,
    pub length: f64,
    pub pressure_diffusivity: f64,
    pub thermal_diffusivity: f64,
    pub species_diffusivity: f64,
    pub reaction_rate: f64,
    pub thermo_pressure_coupling: f64,
    pub reaction_heat: f64,
}

impl Default for Thmc1DConfig {
    fn default() -> Self {
        Self {
            cells: 64,
            length: 1.0,
            pressure_diffusivity: 0.01,
            thermal_diffusivity: 0.01,
            species_diffusivity: 0.005,
            reaction_rate: 0.1,
            thermo_pressure_coupling: 0.05,
            reaction_heat: 0.1,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThmcDiagnostics {
    pub pressure_mean: f64,
    pub temperature_mean: f64,
    pub species_inventory: f64,
    pub reaction_heat: f64,
}

/// Explicit split THMC update with thermo-pressure coupling, first-order reaction heating and a
/// displacement proxy accumulating the pressure change.
#[derive(Clone, Debug)]
pub struct Thmc1D {
    c: Thmc1DConfig,
    pressure: Vec<f64>,
    temperature: Vec<f64>,
    species: Vec<f64>,
    displacement: Vec<f64>,
}

impl Thmc1D {
    pub fn new(c: Thmc1DConfig) -> Result<Self> {
        require(
            c.cells >= 3
                && c.length > 0.0
                && c.pressure_diffusivity >= 0.0
                && c.thermal_diffusivity >= 0.0
                && c.species_diffusivity >= 0.0
                && c.reaction_rate >= 0.0
                && c.reaction_heat >= 0.0,
            "invalid THMC configuration",
        )?;
        let z = vec![0.0; c.cells];
        Ok(Self {
            c,
            pressure: z.clone(),
            temperature: z.clone(),
            species: z.clone(),
            displacement: z,
        })
    }
    pub fn initialize(&mut self, pressure: f64, temperature: f64, species: f64) -> Result<()> {
        require(
            pressure.is_finite()
                && temperature.is_finite()
                && species.is_finite()
                && species >= 0.0,
            "invalid THMC initial state",
        )?;
        self.pressure.fill(pressure);
        self.temperature.fill(temperature);
        self.species.fill(species);
        self.displacement.fill(0.0);
        Ok(())
    }
    pub fn set_cell(
        &mut self,
        cell: usize,
        pressure: f64,
        temperature: f64,
        species: f64,
    ) -> Result<()> {
        require(
            cell < self.c.cells
                && pressure.is_finite()
                && temperature.is_finite()
                && species.is_finite()
                && species >= 0.0,
            "invalid THMC cell",
        )?;
        self.pressure[cell] = pressure;
        self.temperature[cell] = temperature;
        self.species[cell] = species;
        Ok(())
    }
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(dt > 0.0 && dt.is_finite(), "invalid THMC timestep")?;
        let c = self.c;
        let h = c.length / c.cells as f64;
        let inv_h2 = 1.0 / (h * h);
        let (mut np, mut nt, mut ns) = (vec![0.0; c.cells], vec![0.0; c.cells], vec![0.0; c.cells]);
        for i in 0..c.cells {
            let (lp, lt, ls) = (
                laplace(&self.pressure, i, inv_h2),
                laplace(&self.temperature, i, inv_h2),
                laplace(&self.species, i, inv_h2),
            );
            let reaction = c.reaction_rate * self.species[i];
            np[i] = self.pressure[i]
                + dt * (c.pressure_diffusivity * lp + c.thermo_pressure_coupling * lt);
            nt[i] = self.temperature[i]
                + dt * (c.thermal_diffusivity * lt + c.reaction_heat * reaction);
            ns[i] = (self.species[i] + dt * (c.species_diffusivity * ls - reaction)).max(0.0);
        }
        for i in 0..c.cells {
            self.displacement[i] += dt * (np[i] - self.pressure[i]);
        }
        self.pressure = np;
        self.temperature = nt;
        self.species = ns;
        Ok(())
    }
    pub fn pressure(&self) -> &[f64] {
        &self.pressure
    }
    pub fn temperature(&self) -> &[f64] {
        &self.temperature
    }
    pub fn species(&self) -> &[f64] {
        &self.species
    }
    pub fn displacement(&self) -> &[f64] {
        &self.displacement
    }
    pub fn diagnostics(&self) -> ThmcDiagnostics {
        let n = self.c.cells as f64;
        let inv = 1.0 / n;
        let (ps, ts, ss) = (
            self.pressure.iter().sum::<f64>(),
            self.temperature.iter().sum::<f64>(),
            self.species.iter().sum::<f64>(),
        );
        ThmcDiagnostics {
            pressure_mean: ps * inv,
            temperature_mean: ts * inv,
            species_inventory: ss,
            reaction_heat: self.c.reaction_heat * (ss / n),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AllenCahn1DConfig {
    pub cells: usize,
    pub length: f64,
    pub mobility: f64,
    pub interface_width: f64,
}

impl Default for AllenCahn1DConfig {
    fn default() -> Self {
        Self {
            cells: 128,
            length: 1.0,
            mobility: 1.0,
            interface_width: 0.02,
        }
    }
}

/// Explicit periodic Allen-Cahn `phi_t = M (eps^2 phi_xx - phi (phi^2 - 1))`.
#[derive(Clone, Debug)]
pub struct AllenCahn1D {
    c: AllenCahn1DConfig,
    phase: Vec<f64>,
}

impl AllenCahn1D {
    pub fn new(c: AllenCahn1DConfig) -> Result<Self> {
        require(
            c.cells >= 3 && c.length > 0.0 && c.mobility > 0.0 && c.interface_width > 0.0,
            "invalid Allen-Cahn configuration",
        )?;
        Ok(Self {
            c,
            phase: vec![0.0; c.cells],
        })
    }
    pub fn initialize_tanh(&mut self, center: f64, radius: f64) -> Result<()> {
        require(
            center.is_finite() && radius > 0.0,
            "invalid phase-field interface",
        )?;
        let h = self.c.length / self.c.cells as f64;
        for (i, p) in self.phase.iter_mut().enumerate() {
            let x = h * (i as f64 + 0.5);
            *p = ((radius - (x - center).abs()) / self.c.interface_width).tanh();
        }
        Ok(())
    }
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(dt > 0.0 && dt.is_finite(), "invalid Allen-Cahn timestep")?;
        let h = self.c.length / self.c.cells as f64;
        let inv_h2 = 1.0 / (h * h);
        let e2 = self.c.interface_width * self.c.interface_width;
        let next: Vec<f64> = (0..self.phase.len())
            .map(|i| {
                let p = self.phase[i];
                p + dt
                    * self.c.mobility
                    * (e2 * laplace(&self.phase, i, inv_h2) - p * (p * p - 1.0))
            })
            .collect();
        self.phase = next;
        Ok(())
    }
    pub fn phase(&self) -> &[f64] {
        &self.phase
    }
    pub fn free_energy(&self) -> f64 {
        let n = self.phase.len();
        let h = self.c.length / n as f64;
        let e2 = self.c.interface_width * self.c.interface_width;
        (0..n)
            .map(|i| {
                let p = self.phase[i];
                let g = (self.phase[right(i, n)] - p) / h;
                (0.5 * e2 * g * g + 0.25 * (p * p - 1.0) * (p * p - 1.0)) * h
            })
            .sum()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FracturedDarcyEdge {
    pub first: usize,
    pub second: usize,
    pub transmissibility: f64,
    pub aperture: f64,
}

impl Default for FracturedDarcyEdge {
    fn default() -> Self {
        Self {
            first: 0,
            second: 0,
            transmissibility: 1.0,
            aperture: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FracturedDarcyThmcConfig {
    pub cells: usize,
    pub porosity: f64,
    pub pressure_storage: f64,
    pub thermal_capacity: f64,
    pub species_capacity: f64,
    pub thermal_conductivity: f64,
    pub species_diffusivity: f64,
    pub reaction_rate: f64,
    pub reaction_heat: f64,
    pub saturation_mobility: f64,
}

impl Default for FracturedDarcyThmcConfig {
    fn default() -> Self {
        Self {
            cells: 8,
            porosity: 0.2,
            pressure_storage: 1.0,
            thermal_capacity: 1.0,
            species_capacity: 1.0,
            thermal_conductivity: 0.01,
            species_diffusivity: 0.005,
            reaction_rate: 0.1,
            reaction_heat: 0.1,
            saturation_mobility: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FracturedDarcyThmcDiagnostics {
    pub pressure_mass: f64,
    pub saturation_mass: f64,
    pub temperature_mass: f64,
    pub species_mass: f64,
    pub absolute_flux: f64,
    pub reaction_consumption: f64,
}

/// Conservative graph finite volumes; edges with aperture above one act as fracture conduits.
#[derive(Clone, Debug)]
pub struct FracturedDarcyThmcNetwork {
    c: FracturedDarcyThmcConfig,
    edges: Vec<FracturedDarcyEdge>,
    pressure: Vec<f64>,
    saturation: Vec<f64>,
    temperature: Vec<f64>,
    species: Vec<f64>,
    displacement: Vec<f64>,
    absolute_flux: f64,
    reaction_consumption: f64,
}

impl FracturedDarcyThmcNetwork {
    pub fn new(c: FracturedDarcyThmcConfig) -> Result<Self> {
        require(
            c.cells > 0
                && c.porosity > 0.0
                && c.porosity < 1.0
                && c.pressure_storage > 0.0
                && c.thermal_capacity > 0.0
                && c.species_capacity > 0.0
                && c.thermal_conductivity >= 0.0
                && c.species_diffusivity >= 0.0
                && c.reaction_rate >= 0.0
                && c.reaction_heat >= 0.0
                && c.saturation_mobility >= 0.0,
            "invalid fractured Darcy THMC configuration",
        )?;
        let z = vec![0.0; c.cells];
        Ok(Self {
            c,
            edges: Vec::new(),
            pressure: z.clone(),
            saturation: z.clone(),
            temperature: z.clone(),
            species: z.clone(),
            displacement: z,
            absolute_flux: 0.0,
            reaction_consumption: 0.0,
        })
    }
    pub fn add_edge(&mut self, e: FracturedDarcyEdge) -> Result<()> {
        require(
            e.first < self.c.cells
                && e.second < self.c.cells
                && e.first != e.second
                && e.transmissibility > 0.0
                && e.aperture > 0.0
                && e.transmissibility.is_finite()
                && e.aperture.is_finite(),
            "invalid fractured Darcy edge",
        )?;
        self.edges.push(e);
        Ok(())
    }
    pub fn initialize(
        &mut self,
        pressure: f64,
        saturation: f64,
        temperature: f64,
        species: f64,
    ) -> Result<()> {
        require(
            pressure.is_finite()
                && (0.0..=1.0).contains(&saturation)
                && temperature.is_finite()
                && species >= 0.0
                && species.is_finite(),
            "invalid fractured Darcy initial state",
        )?;
        self.pressure.fill(pressure);
        self.saturation.fill(saturation);
        self.temperature.fill(temperature);
        self.species.fill(species);
        self.displacement.fill(0.0);
        Ok(())
    }
    pub fn set_cell(
        &mut self,
        cell: usize,
        pressure: f64,
        saturation: f64,
        temperature: f64,
        species: f64,
    ) -> Result<()> {
        require(
            cell < self.c.cells
                && pressure.is_finite()
                && (0.0..=1.0).contains(&saturation)
                && temperature.is_finite()
                && species >= 0.0
                && species.is_finite(),
            "invalid fractured Darcy cell",
        )?;
        self.pressure[cell] = pressure;
        self.saturation[cell] = saturation;
        self.temperature[cell] = temperature;
        self.species[cell] = species;
        Ok(())
    }
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(
            dt > 0.0 && dt.is_finite(),
            "invalid fractured Darcy timestep",
        )?;
        let c = self.c;
        let n = c.cells;
        let (mut dp, mut ds, mut dtemp, mut dc) =
            (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        self.absolute_flux = 0.0;
        for e in &self.edges {
            let (a, b) = (e.first, e.second);
            let fracture = e.aperture.max(1.0);
            let conductance = e.transmissibility * fracture;
            let pf = conductance * (self.pressure[a] - self.pressure[b]);
            dp[a] -= pf;
            dp[b] += pf;
            let up = if pf >= 0.0 {
                self.saturation[a]
            } else {
                self.saturation[b]
            };
            let sf = c.saturation_mobility * fracture * pf * up;
            ds[a] -= sf;
            ds[b] += sf;
            let tf =
                c.thermal_conductivity * conductance * (self.temperature[a] - self.temperature[b]);
            dtemp[a] -= tf;
            dtemp[b] += tf;
            let cf = c.species_diffusivity * conductance * (self.species[a] - self.species[b]);
            dc[a] -= cf;
            dc[b] += cf;
            self.absolute_flux += pf.abs();
        }
        self.reaction_consumption = 0.0;
        for i in 0..n {
            let reaction = c.reaction_rate * self.species[i];
            self.pressure[i] += dt * dp[i] / c.pressure_storage;
            self.saturation[i] = (self.saturation[i] + dt * ds[i] / c.porosity).clamp(0.0, 1.0);
            self.temperature[i] +=
                dt * (dtemp[i] / c.thermal_capacity + c.reaction_heat * reaction);
            self.species[i] =
                (self.species[i] + dt * (dc[i] / c.species_capacity - reaction)).max(0.0);
            self.displacement[i] += dt * dp[i];
            self.reaction_consumption += reaction;
        }
        Ok(())
    }
    pub fn pressure(&self) -> &[f64] {
        &self.pressure
    }
    pub fn saturation(&self) -> &[f64] {
        &self.saturation
    }
    pub fn temperature(&self) -> &[f64] {
        &self.temperature
    }
    pub fn species(&self) -> &[f64] {
        &self.species
    }
    pub fn displacement(&self) -> &[f64] {
        &self.displacement
    }
    pub fn diagnostics(&self) -> FracturedDarcyThmcDiagnostics {
        let c = self.c;
        let mut d = FracturedDarcyThmcDiagnostics {
            absolute_flux: self.absolute_flux,
            reaction_consumption: self.reaction_consumption,
            ..Default::default()
        };
        for i in 0..c.cells {
            d.pressure_mass += self.pressure[i] * c.pressure_storage;
            d.saturation_mass += self.saturation[i] * c.porosity;
            d.temperature_mass += self.temperature[i] * c.thermal_capacity;
            d.species_mass += self.species[i] * c.species_capacity;
        }
        d
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CalphadBinaryPoint {
    pub composition: f64,
    pub free_energy: f64,
}

/// Piecewise-linear tabular free energy on `[0, 1]`; chemical potential is the interval slope.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryCalphadTable {
    points: Vec<CalphadBinaryPoint>,
}

impl BinaryCalphadTable {
    pub fn new(points: Vec<CalphadBinaryPoint>) -> Result<Self> {
        require(points.len() >= 2, "CALPHAD table needs two points")?;
        for (i, p) in points.iter().enumerate() {
            require(
                (0.0..=1.0).contains(&p.composition)
                    && p.free_energy.is_finite()
                    && (i == 0 || p.composition > points[i - 1].composition),
                "invalid CALPHAD composition table",
            )?;
        }
        Ok(Self { points })
    }
    fn locate(&self, x: f64) -> Result<usize> {
        require(
            (0.0..=1.0).contains(&x),
            "CALPHAD composition outside [0,1]",
        )?;
        Ok((0..self.points.len() - 1)
            .find(|&i| x <= self.points[i + 1].composition)
            .unwrap_or(self.points.len() - 2))
    }
    pub fn free_energy(&self, composition: f64) -> Result<f64> {
        let i = self.locate(composition)?;
        let (l, r) = (self.points[i], self.points[i + 1]);
        let w = (composition - l.composition) / (r.composition - l.composition);
        Ok(l.free_energy + w * (r.free_energy - l.free_energy))
    }
    pub fn chemical_potential(&self, composition: f64) -> Result<f64> {
        let i = self.locate(composition)?;
        let (l, r) = (self.points[i], self.points[i + 1]);
        Ok((r.free_energy - l.free_energy) / (r.composition - l.composition))
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CahnHilliard1DConfig {
    pub cells: usize,
    pub length: f64,
    pub mobility: f64,
    pub gradient_energy: f64,
    pub anisotropy_factor: f64,
    pub nucleation_amplitude: f64,
}

impl Default for CahnHilliard1DConfig {
    fn default() -> Self {
        Self {
            cells: 64,
            length: 1.0,
            mobility: 1.0,
            gradient_energy: 0.01,
            anisotropy_factor: 1.0,
            nucleation_amplitude: 0.0,
        }
    }
}

const C_MIN: f64 = 1e-12;
const C_MAX: f64 = 1.0 - 1e-12;

/// Conservative periodic Cahn-Hilliard with a double-well or tabular CALPHAD free energy; the
/// mean is restored after clamping.
#[derive(Clone, Debug)]
pub struct CahnHilliard1D {
    c: CahnHilliard1DConfig,
    table: Option<BinaryCalphadTable>,
    composition: Vec<f64>,
    chemical_potential: Vec<f64>,
}

impl CahnHilliard1D {
    pub fn new(c: CahnHilliard1DConfig, table: Option<BinaryCalphadTable>) -> Result<Self> {
        require(
            c.cells >= 4
                && c.length > 0.0
                && c.mobility > 0.0
                && c.gradient_energy >= 0.0
                && c.anisotropy_factor > 0.0
                && c.nucleation_amplitude >= 0.0,
            "invalid Cahn-Hilliard configuration",
        )?;
        Ok(Self {
            c,
            table,
            composition: vec![0.0; c.cells],
            chemical_potential: vec![0.0; c.cells],
        })
    }
    /// `mean + perturbation cos(2 pi i / N)` clamped to the open unit interval.
    pub fn initialize(&mut self, mean: f64, perturbation: f64) -> Result<()> {
        require(
            (0.0..=1.0).contains(&mean) && perturbation >= 0.0 && perturbation.is_finite(),
            "invalid Cahn-Hilliard initial composition",
        )?;
        let n = self.composition.len() as f64;
        for (i, v) in self.composition.iter_mut().enumerate() {
            let phase = 2.0 * std::f64::consts::PI * i as f64 / n;
            *v = (mean + perturbation * phase.cos()).clamp(C_MIN, C_MAX);
        }
        Ok(())
    }
    pub fn seed_nucleus(&mut self, center: usize, radius: f64, amplitude: f64) -> Result<()> {
        require(
            center < self.composition.len()
                && radius >= 0.0
                && amplitude >= 0.0
                && radius.is_finite()
                && amplitude.is_finite(),
            "invalid Cahn-Hilliard nucleus",
        )?;
        for (i, v) in self.composition.iter_mut().enumerate() {
            if (i as f64 - center as f64).abs() <= radius {
                *v = (*v + amplitude).clamp(C_MIN, C_MAX);
            }
        }
        Ok(())
    }
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(dt > 0.0 && dt.is_finite(), "invalid Cahn-Hilliard timestep")?;
        let n = self.composition.len();
        let h = self.c.length / n as f64;
        let inv_h2 = 1.0 / (h * h);
        for i in 0..n {
            let x = self.composition[i];
            let homogeneous = match &self.table {
                None => x * (x - 0.5) * (x - 1.0),
                Some(t) => t.chemical_potential(x)?,
            };
            self.chemical_potential[i] = homogeneous
                - self.c.gradient_energy
                    * self.c.anisotropy_factor
                    * laplace(&self.composition, i, inv_h2);
        }
        let mut next: Vec<f64> = (0..n)
            .map(|i| {
                (self.composition[i]
                    + dt * self.c.mobility * laplace(&self.chemical_potential, i, inv_h2))
                .clamp(C_MIN, C_MAX)
            })
            .collect();
        let old_mean = self.mean();
        let new_mean = next.iter().sum::<f64>() / n as f64;
        let correction = old_mean - new_mean;
        for v in &mut next {
            *v = (*v + correction).clamp(C_MIN, C_MAX);
        }
        self.composition = next;
        Ok(())
    }
    pub fn composition(&self) -> &[f64] {
        &self.composition
    }
    pub fn chemical_potential(&self) -> &[f64] {
        &self.chemical_potential
    }
    pub fn mean(&self) -> f64 {
        self.composition.iter().sum::<f64>() / self.composition.len() as f64
    }
    pub fn free_energy(&self) -> Result<f64> {
        let n = self.composition.len();
        let h = self.c.length / n as f64;
        let mut total = 0.0;
        for i in 0..n {
            let x = self.composition[i];
            let homogeneous = match &self.table {
                None => 0.25 * x * x * (1.0 - x) * (1.0 - x),
                Some(t) => t.free_energy(x)?,
            };
            let g = (self.composition[right(i, n)] - x) / h;
            total +=
                (homogeneous + 0.5 * self.c.gradient_energy * self.c.anisotropy_factor * g * g) * h;
        }
        Ok(total)
    }
}
