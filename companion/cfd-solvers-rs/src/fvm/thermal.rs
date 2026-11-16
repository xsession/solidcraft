use super::{
    schemes::flux_sizes, PolyMesh, RadiationSourceModel, ScalarBoundaryType, ScalarTransport,
    ScalarTransportConfig, TemporalScheme,
};
use crate::{
    core::{IterativeSolverResult, LinearConfig, Vec3},
    finite_slice, positive, require, Result,
};
use serde::{Deserialize, Serialize};
use std::sync::Arc;
#[derive(Clone, Copy, Debug, Serialize, Deserialize, Default)]
pub struct ThermalBoundary {
    pub fixed_temperature: bool,
    pub temperature: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ThermalTransportConfig {
    pub dt: f64,
    pub density: f64,
    pub default_heat_capacity: f64,
    pub default_conductivity: f64,
}
impl Default for ThermalTransportConfig {
    fn default() -> Self {
        Self {
            dt: 1e-4,
            density: 1.0,
            default_heat_capacity: 1.0,
            default_conductivity: 1.0,
        }
    }
}
pub type HeatRelease<'a> = dyn Fn(Vec3, f64, f64) -> Result<f64> + 'a;
type MaterialProperty<'a> = dyn Fn(Vec3, f64) -> f64 + 'a;
pub struct ThermalTransport<'a> {
    mesh: PolyMesh,
    config: ThermalTransportConfig,
    temperature: Vec<f64>,
    mass_flux: Vec<f64>,
    boundary: Vec<ThermalBoundary>,
    source: Option<Box<HeatRelease<'a>>>,
    cp: Box<MaterialProperty<'a>>,
    k: Box<MaterialProperty<'a>>,
    time: f64,
}
impl<'a> ThermalTransport<'a> {
    pub fn new(mesh: PolyMesh, config: ThermalTransportConfig) -> Result<Self> {
        positive(config.dt, "thermal dt")?;
        positive(config.density, "density")?;
        positive(config.default_heat_capacity, "heat capacity")?;
        require(
            config.default_conductivity.is_finite() && config.default_conductivity >= 0.0,
            "invalid conductivity",
        )?;
        let n = mesh.cell_count();
        let nf = mesh.face_count();
        let np = mesh.patches().len();
        Ok(Self {
            mesh,
            config,
            temperature: vec![0.0; n],
            mass_flux: vec![0.0; nf],
            boundary: vec![ThermalBoundary::default(); np],
            source: None,
            cp: Box::new(move |_, _| config.default_heat_capacity),
            k: Box::new(move |_, _| config.default_conductivity),
            time: 0.0,
        })
    }
    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn temperature(&self) -> &[f64] {
        &self.temperature
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn initialize_constant(&mut self, t: f64) -> Result<()> {
        self.initialize_profile(|_| t)
    }
    pub fn initialize_profile(&mut self, f: impl Fn(Vec3) -> f64) -> Result<()> {
        let t: Vec<_> = self.mesh.cells().iter().map(|c| f(c.center)).collect();
        finite_slice(&t, "thermal initialization")?;
        self.temperature = t;
        self.time = 0.0;
        Ok(())
    }
    pub fn set_boundary(&mut self, patch: &str, b: ThermalBoundary) -> Result<()> {
        require(b.temperature.is_finite(), "nonfinite thermal boundary")?;
        let p = self.mesh.patch_index(patch)?;
        self.boundary[p] = b;
        Ok(())
    }
    pub fn set_face_mass_flux(&mut self, f: Vec<f64>) -> Result<()> {
        flux_sizes(&self.mesh, &f)?;
        self.mass_flux = f;
        Ok(())
    }
    pub fn set_volumetric_source(&mut self, f: impl Fn(Vec3, f64) -> f64 + 'a) {
        self.source = Some(Box::new(move |p, t, _| Ok(f(p, t))));
    }
    pub fn set_temperature_dependent_source(&mut self, f: impl Fn(Vec3, f64, f64) -> f64 + 'a) {
        self.source = Some(Box::new(move |p, t, u| Ok(f(p, t, u))));
    }
    pub fn clear_source(&mut self) {
        self.source = None;
    }
    pub fn set_reactive_radiative_source(
        &mut self,
        chemistry: Option<Box<HeatRelease<'a>>>,
        radiation: Arc<dyn RadiationSourceModel + Send + Sync + 'a>,
        environment: f64,
    ) -> Result<()> {
        require(
            environment.is_finite() && environment >= 0.0,
            "radiation environment temperature",
        )?;
        self.source = Some(Box::new(move |p, t, temperature| {
            let chemical = if let Some(f) = &chemistry {
                f(p, t, temperature)?
            } else {
                0.0
            };
            let radiative = radiation.volumetric_source(temperature, environment)?;
            let sum = chemical + radiative;
            require(sum.is_finite(), "reactive-radiative source overflow")?;
            Ok(sum)
        }));
        Ok(())
    }
    pub fn set_material(
        &mut self,
        cp: impl Fn(Vec3, f64) -> f64 + 'a,
        k: impl Fn(Vec3, f64) -> f64 + 'a,
    ) {
        self.cp = Box::new(cp);
        self.k = Box::new(k);
    }
    fn properties(&self) -> Result<(Vec<f64>, Vec<f64>)> {
        let mut cp = vec![];
        let mut k = vec![];
        for (c, t) in self.mesh.cells().iter().zip(&self.temperature) {
            let heat = (self.cp)(c.center, *t);
            let conductivity = (self.k)(c.center, *t);
            positive(heat, "temperature-dependent heat capacity")?;
            require(
                conductivity.is_finite() && conductivity >= 0.0,
                "invalid temperature-dependent conductivity",
            )?;
            cp.push(heat);
            k.push(conductivity);
        }
        Ok((cp, k))
    }
    /// Conduction/upwind diagonal bound only; nonlinear sources need their own time-step analysis.
    pub fn transport_stable_timestep(&self, safety: f64) -> Result<Option<f64>> {
        require(
            safety > 0.0 && safety <= 1.0,
            "thermal safety factor must be in (0,1]",
        )?;
        let (cp, k) = self.properties()?;
        let mut dt = f64::INFINITY;
        for c in 0..self.mesh.cell_count() {
            let mut diagonal = 0.0;
            for &fi in &self.mesh.cell_faces()[c] {
                let f = self.mesh.faces()[fi];
                if let Some(n) = f.neighbour {
                    let other = if f.owner == c { n } else { f.owner };
                    let kf = 2.0 * k[c] * k[other] / (k[c] + k[other]).max(1e-30);
                    let distance =
                        (self.mesh.cells()[other].center - self.mesh.cells()[c].center).norm();
                    positive(distance, "thermal neighbour distance")?;
                    diagonal += kf * f.area.norm() / distance;
                } else if self.boundary[f.patch.ok_or_else(|| {
                    crate::SolverError::InvalidInput("boundary face without patch".into())
                })?]
                .fixed_temperature
                {
                    let distance = (f.center - self.mesh.cells()[c].center).norm();
                    positive(distance, "thermal boundary distance")?;
                    diagonal += k[c] * f.area.norm() / distance;
                }
                let phi = if f.owner == c {
                    self.mass_flux[fi]
                } else {
                    -self.mass_flux[fi]
                };
                diagonal += phi.max(0.0) * cp[c];
            }
            if diagonal > 0.0 {
                dt = dt.min(
                    safety * self.config.density * cp[c] * self.mesh.cells()[c].volume / diagonal,
                );
            }
        }
        Ok(if dt.is_finite() { Some(dt) } else { None })
    }
    pub fn step(&mut self) -> Result<()> {
        let (cp, k) = self.properties()?;
        let next_time = self.time + self.config.dt;
        let mut next = self.temperature.clone();
        for c in 0..self.mesh.cell_count() {
            let cell = self.mesh.cells()[c];
            let source = if let Some(f) = &self.source {
                f(cell.center, next_time, self.temperature[c])?
            } else {
                0.0
            };
            require(source.is_finite(), "nonfinite thermal source")?;
            let mut rate = source * cell.volume;
            for &fi in &self.mesh.cell_faces()[c] {
                let face = self.mesh.faces()[fi];
                let phi = if face.owner == c {
                    self.mass_flux[fi]
                } else {
                    -self.mass_flux[fi]
                };
                if let Some(n) = face.neighbour {
                    let other = if face.owner == c { n } else { face.owner };
                    let kf = 2.0 * k[c] * k[other] / (k[c] + k[other]).max(1e-30);
                    let d = (self.mesh.cells()[other].center - cell.center).norm();
                    positive(d, "thermal centre distance")?;
                    let conduction =
                        kf * face.area.norm() * (self.temperature[other] - self.temperature[c]) / d;
                    let (up, up_cp) = if phi >= 0.0 {
                        (self.temperature[c], cp[c])
                    } else {
                        (self.temperature[other], cp[other])
                    };
                    rate += conduction - phi * up_cp * up;
                } else {
                    let b = self.boundary[face.patch.ok_or_else(|| {
                        crate::SolverError::InvalidInput("boundary face without patch".into())
                    })?];
                    if b.fixed_temperature {
                        let d = (face.center - cell.center).norm();
                        positive(d, "thermal boundary distance")?;
                        rate += k[c] * face.area.norm() * (b.temperature - self.temperature[c]) / d;
                    }
                    if phi > 0.0 {
                        rate -= phi * cp[c] * self.temperature[c];
                    } else if phi < 0.0 && b.fixed_temperature {
                        rate -= phi * cp[c] * b.temperature;
                    }
                }
            }
            next[c] = self.temperature[c]
                + self.config.dt * rate / (self.config.density * cp[c] * cell.volume);
        }
        finite_slice(&next, "thermal transport update")?;
        self.temperature = next;
        self.time = next_time;
        Ok(())
    }
    pub fn run(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }
    pub fn total_sensible_enthalpy(&self) -> Result<f64> {
        let (cp, _) = self.properties()?;
        let h = self
            .temperature
            .iter()
            .zip(cp)
            .zip(self.mesh.cells())
            .map(|((t, cp), c)| self.config.density * cp * t * c.volume)
            .sum::<f64>();
        require(h.is_finite(), "sensible heat overflow")?;
        Ok(h)
    }
}
pub fn boussinesq_acceleration(
    mesh: &PolyMesh,
    temperature: &[f64],
    reference: f64,
    beta: f64,
    gravity: Vec3,
) -> Result<Vec<Vec3>> {
    require(
        temperature.len() == mesh.cell_count()
            && reference.is_finite()
            && beta.is_finite()
            && gravity.is_finite(),
        "Boussinesq input size/parameters",
    )?;
    finite_slice(temperature, "Boussinesq temperature")?;
    let a: Vec<_> = temperature
        .iter()
        .map(|t| gravity * (-beta * (t - reference)))
        .collect();
    require(a.iter().all(|v| v.is_finite()), "Boussinesq overflow")?;
    Ok(a)
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SolidHeatConfig {
    pub dt: f64,
    pub density: f64,
    pub heat_capacity: f64,
    pub conductivity: f64,
    pub temporal_scheme: TemporalScheme,
    pub crank_nicolson_off_centering: f64,
    pub linear: LinearConfig,
}
impl Default for SolidHeatConfig {
    fn default() -> Self {
        Self {
            dt: 1e-3,
            density: 1.0,
            heat_capacity: 1.0,
            conductivity: 1.0,
            temporal_scheme: TemporalScheme::BackwardBdf2,
            crank_nicolson_off_centering: 1.0,
            linear: LinearConfig {
                max_iterations: 400,
                restart: 30,
                relative_tolerance: 1e-10,
            },
        }
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SolidThermalBoundaryType {
    Adiabatic,
    FixedTemperature,
}
#[derive(Clone, Debug)]
pub struct SolidHeatConduction {
    config: SolidHeatConfig,
    solver: ScalarTransport,
}
impl SolidHeatConduction {
    pub fn new(mesh: PolyMesh, c: SolidHeatConfig) -> Result<Self> {
        positive(c.density, "solid density")?;
        positive(c.heat_capacity, "solid heat capacity")?;
        let scalar = ScalarTransportConfig {
            dt: c.dt,
            diffusivity: c.conductivity / (c.density * c.heat_capacity),
            linear: c.linear,
            temporal_scheme: c.temporal_scheme,
            crank_nicolson_off_centering: c.crank_nicolson_off_centering,
        };
        Ok(Self {
            config: c,
            solver: ScalarTransport::new(mesh, scalar)?,
        })
    }
    pub fn initialize_constant(&mut self, t: f64) -> Result<()> {
        self.solver.initialize_constant(t)
    }
    pub fn initialize_profile(&mut self, f: impl Fn(Vec3) -> f64) -> Result<()> {
        self.solver.initialize_profile(f)
    }
    pub fn set_boundary(
        &mut self,
        patch: &str,
        kind: SolidThermalBoundaryType,
        temperature: f64,
    ) -> Result<()> {
        self.solver.set_boundary(
            patch,
            match kind {
                SolidThermalBoundaryType::Adiabatic => ScalarBoundaryType::ZeroGradient,
                SolidThermalBoundaryType::FixedTemperature => ScalarBoundaryType::FixedValue,
            },
            temperature,
        )
    }
    pub fn set_volumetric_heat_source_constant(&mut self, q: f64) -> Result<()> {
        self.set_volumetric_heat_source_profile(|_| q)
    }
    pub fn set_volumetric_heat_source_profile(&mut self, f: impl Fn(Vec3) -> f64) -> Result<()> {
        let capacity = self.config.density * self.config.heat_capacity;
        self.solver.set_source_profile(|p| f(p) / capacity)
    }
    pub fn step(&mut self) -> Result<IterativeSolverResult> {
        self.solver.step()
    }
    pub fn run(&mut self, n: usize) -> Result<Option<IterativeSolverResult>> {
        self.solver.run(n)
    }
    pub fn mesh(&self) -> &PolyMesh {
        self.solver.mesh()
    }
    pub fn temperature(&self) -> &[f64] {
        self.solver.values()
    }
    pub fn time(&self) -> f64 {
        self.solver.time()
    }
    pub fn config(&self) -> SolidHeatConfig {
        self.config
    }
    pub fn total_sensible_energy(&self) -> f64 {
        self.config.density * self.config.heat_capacity * self.solver.volume_integral()
    }
}
