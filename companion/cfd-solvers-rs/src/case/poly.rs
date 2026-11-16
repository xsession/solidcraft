//! Typed JSON access to the general FVM engine. Geometry is explicit and exportable.
use super::{write_field, RunManifest};
use crate::{
    core::*,
    fvm::*,
    io::{self, DataArray, VtuMesh},
    positive, require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, path::Path, sync::Arc};

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PolyMeshSpec {
    pub dimensions: [usize; 3],
    pub lengths: Vec3,
    pub xy_shear: f64,
}
impl Default for PolyMeshSpec {
    fn default() -> Self {
        Self {
            dimensions: [8, 6, 2],
            lengths: Vec3::new(1.0, 1.0, 0.1),
            xy_shear: 0.0,
        }
    }
}
impl PolyMeshSpec {
    pub fn build(self) -> Result<(PolyMesh, VtuMesh)> {
        let [nx, ny, nz] = self.dimensions;
        let l = self.lengths;
        let m = make_sheared_cartesian_hexa_mesh(nx, ny, nz, l.x, l.y, l.z, self.xy_shear)?;
        let export = VtuMesh::box_cells(self.dimensions, l, self.xy_shear)?;
        Ok((m, export))
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PolyOutput {
    pub steps: usize,
    pub output_interval: usize,
}
impl Default for PolyOutput {
    fn default() -> Self {
        Self {
            steps: 20,
            output_interval: 5,
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScalarPatch {
    pub patch: String,
    pub kind: ScalarBoundaryType,
    #[serde(default)]
    pub value: f64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VelocityPatch {
    pub patch: String,
    pub kind: VelocityBoundaryType,
    #[serde(default)]
    pub value: Vec3,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PressurePatch {
    pub patch: String,
    pub kind: PressureBoundaryType,
    #[serde(default)]
    pub value: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FlowMethod {
    Simple,
    Piso,
    Pimple,
}
impl Default for FlowMethod {
    fn default() -> Self {
        Self::Piso
    }
}
fn unique<'a>(patches: impl Iterator<Item = &'a str>) -> Result<()> {
    let mut seen = BTreeSet::new();
    for p in patches {
        require(seen.insert(p), format!("duplicate boundary patch {p}"))?;
    }
    Ok(())
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PolyScalarCase {
    pub mesh: PolyMeshSpec,
    pub config: ScalarTransportConfig,
    pub output: PolyOutput,
    pub initial: f64,
    pub initial_gradient: Vec3,
    pub velocity: Vec3,
    pub source: f64,
    pub boundaries: Vec<ScalarPatch>,
}
impl Default for PolyScalarCase {
    fn default() -> Self {
        Self {
            mesh: Default::default(),
            config: Default::default(),
            output: Default::default(),
            initial: 1.0,
            initial_gradient: Vec3::ZERO,
            velocity: Vec3::ZERO,
            source: 0.0,
            boundaries: vec![],
        }
    }
}
/// A grid-aligned solid, expressed as half-open integer cell-index bounds.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GridObstacle {
    pub name: String,
    pub min: [usize; 3],
    pub max: [usize; 3],
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CollocatedFlowCase {
    pub mesh: PolyMeshSpec,
    /// Solid blocks in half-open background-grid cell coordinates.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub obstacles: Vec<GridObstacle>,
    /// Solids from closed triangle surfaces in metres (`surface_bodies_v1`).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub bodies: Vec<super::bodies::SurfaceBody>,
    pub config: CollocatedIncompressibleConfig,
    pub output: PolyOutput,
    pub method: FlowMethod,
    pub initial_velocity: Vec3,
    pub initial_pressure: f64,
    pub acceleration: Vec3,
    pub velocity_boundaries: Vec<VelocityPatch>,
    pub pressure_boundaries: Vec<PressurePatch>,
    pub steady_velocity_tolerance: f64,
    pub steady_continuity_tolerance: f64,
    /// Optional external fan curve on one planar exterior inlet (`fan_curve_v1`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fan: Option<super::fan::FanInlet>,
    /// Optional isotropic Darcy-Forchheimer zones on grid cells (`porous_blocks_v1`).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub porous_blocks: Vec<super::porous::PorousBlock>,
    /// Optional steady conjugate heat transfer on the settled flow (`conjugate_heat_v1`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub thermal: Option<super::conjugate::ConjugateHeat>,
    /// Monitored engineering goals (`goals_v1`).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub goals: Vec<super::goals::Goal>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub goal_control: Option<super::goals::GoalControl>,
}
impl Default for CollocatedFlowCase {
    fn default() -> Self {
        Self {
            mesh: Default::default(),
            obstacles: vec![],
            bodies: vec![],
            config: Default::default(),
            output: Default::default(),
            method: FlowMethod::Piso,
            initial_velocity: Vec3::ZERO,
            initial_pressure: 0.0,
            acceleration: Vec3::ZERO,
            velocity_boundaries: vec![],
            pressure_boundaries: vec![],
            steady_velocity_tolerance: 1e-7,
            steady_continuity_tolerance: 1e-8,
            fan: None,
            porous_blocks: vec![],
            thermal: None,
            goals: vec![],
            goal_control: None,
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SolidHeatCase {
    pub mesh: PolyMeshSpec,
    pub config: SolidHeatConfig,
    pub output: PolyOutput,
    pub initial: f64,
    pub initial_gradient: Vec3,
    pub source: f64,
    pub boundaries: Vec<ScalarPatch>,
}
impl Default for SolidHeatCase {
    fn default() -> Self {
        Self {
            mesh: Default::default(),
            config: Default::default(),
            output: Default::default(),
            initial: 300.0,
            initial_gradient: Vec3::ZERO,
            source: 0.0,
            boundaries: vec![],
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ThermalTransportCase {
    pub mesh: PolyMeshSpec,
    pub config: ThermalTransportConfig,
    pub output: PolyOutput,
    pub initial: f64,
    pub initial_gradient: Vec3,
    pub velocity: Vec3,
    pub source: f64,
    pub absorption: f64,
    pub environment_temperature: f64,
    pub boundaries: Vec<ScalarPatch>,
}
impl Default for ThermalTransportCase {
    fn default() -> Self {
        Self {
            mesh: Default::default(),
            config: Default::default(),
            output: Default::default(),
            initial: 300.0,
            initial_gradient: Vec3::ZERO,
            velocity: Vec3::ZERO,
            source: 0.0,
            absorption: 0.0,
            environment_temperature: 300.0,
            boundaries: vec![],
        }
    }
}

// Every step gets a diagnostic record; only selected snapshots are written to VTU.
// A failed linear solve returns an error before a success manifest is written.
fn series(
    frame: impl FnMut(usize) -> Result<(f64, Vec<DataArray>, serde_json::Value, Option<bool>)>,
    out: &Path,
    m: &mut RunManifest,
    mesh: &VtuMesh,
    control: PolyOutput,
    stem: &str,
) -> Result<()> {
    series_until(
        frame,
        out,
        m,
        mesh,
        control,
        stem,
        &std::cell::Cell::new(false),
    )
}

/// Like [`series`], but a frame may request an early stop through `stop`; the stopping step
/// is then written as the final snapshot.
fn series_until(
    mut frame: impl FnMut(usize) -> Result<(f64, Vec<DataArray>, serde_json::Value, Option<bool>)>,
    out: &Path,
    m: &mut RunManifest,
    mesh: &VtuMesh,
    control: PolyOutput,
    stem: &str,
    stop: &std::cell::Cell<bool>,
) -> Result<()> {
    require(
        control.output_interval > 0,
        "output interval must be positive",
    )?;
    let mut frames = vec![];
    let mut reports = vec![];
    for i in 0..=control.steps {
        let (time, data, report, linear) = frame(i)?;
        require(time.is_finite(), "nonfinite output time")?;
        if let Some(converged) = linear {
            m.linear_solver_converged = Some(converged);
            if !converged {
                return Err(SolverError::NumericalFailure(
                    "linear solve did not converge; see case controls".into(),
                ));
            }
        }
        reports.push(serde_json::json!({"step":i,"coordinate":time,"report":report}));
        let stopping = stop.get();
        if i == 0 || i % control.output_interval == 0 || i == control.steps || stopping {
            let file = format!("{stem}_{i:06}.vtu");
            write_field(out, m, &file, mesh, vec![], data)?;
            frames.push((time, file));
        }
        if stopping {
            break;
        }
    }
    let pvd = format!("{stem}.pvd");
    io::write_pvd(out.join(&pvd), &frames)?;
    m.files.push(pvd);
    let file = format!("{stem}_diagnostics.json");
    std::fs::write(out.join(&file), serde_json::to_vec_pretty(&reports)?)?;
    m.files.push(file);
    Ok(())
}
impl PolyScalarCase {
    pub(super) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        let (mesh, export) = self.mesh.build()?;
        require(
            self.velocity.is_finite() && self.initial_gradient.is_finite(),
            "nonfinite scalar profile/velocity",
        )?;
        unique(self.boundaries.iter().map(|b| b.patch.as_str()))?;
        let mut s = ScalarTransport::new(mesh.clone(), self.config)?;
        s.initialize_profile(|p| self.initial + self.initial_gradient.dot(p))?;
        s.set_source_constant(self.source)?;
        s.set_face_flux(
            mesh.faces()
                .iter()
                .map(|f| self.velocity.dot(f.area))
                .collect(),
        )?;
        for b in &self.boundaries {
            s.set_boundary(&b.patch, b.kind, b.value)?;
        }
        series(
            |i| {
                let r = if i > 0 { Some(s.step()?) } else { None };
                Ok((
                    s.time(),
                    vec![DataArray::scalar("scalar", s.values().to_vec())],
                    serde_json::json!({"linear":r,"volume_integral":s.volume_integral(),"minimum":s.minimum(),"maximum":s.maximum()}),
                    r.map(|r| r.converged),
                ))
            },
            out,
            m,
            &export,
            self.output,
            "transport",
        )?;
        m.metrics.insert("time_s".into(), s.time());
        m.metrics
            .insert("volume_integral".into(), s.volume_integral());
        m.metrics
            .insert("operator_assemblies".into(), s.operator_assemblies() as f64);
        Ok(())
    }
}
impl CollocatedFlowCase {
    /// Build both solver and visualization topology from one retained cell set.
    /// Empty obstacles retain the original box/sheared-box behavior.
    pub fn build_mesh(&self) -> Result<(PolyMesh, VtuMesh)> {
        super::obstacles::build(self.mesh, &self.obstacles, &self.bodies)
    }

    /// Names of all solids: grid boxes first, then surface bodies.
    pub fn solid_names(&self) -> Vec<&str> {
        self.obstacles
            .iter()
            .map(|o| o.name.as_str())
            .chain(self.bodies.iter().map(|b| b.name.as_str()))
            .collect()
    }

    pub(super) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        let (mesh, export) = self.build_mesh()?;
        positive(self.steady_velocity_tolerance, "steady velocity tolerance")?;
        positive(
            self.steady_continuity_tolerance,
            "steady continuity tolerance",
        )?;
        unique(self.velocity_boundaries.iter().map(|b| b.patch.as_str()))?;
        unique(self.pressure_boundaries.iter().map(|b| b.patch.as_str()))?;
        if !self.obstacles.is_empty() || !self.bodies.is_empty() {
            for name in self
                .velocity_boundaries
                .iter()
                .map(|b| b.patch.as_str())
                .chain(self.pressure_boundaries.iter().map(|b| b.patch.as_str()))
            {
                let patch = mesh.patch_index(name)?;
                require(
                    mesh.faces().iter().any(|face| face.patch == Some(patch)),
                    format!("obstacles remove every face of configured boundary patch {name}"),
                )?;
            }
        }
        for obstacle in self.solid_names() {
            for b in self
                .velocity_boundaries
                .iter()
                .filter(|b| b.patch == obstacle)
            {
                require(
                    b.kind == VelocityBoundaryType::FixedValue && b.value == Vec3::ZERO,
                    "obstacle velocity boundaries must be stationary no-slip walls",
                )?;
            }
            for b in self
                .pressure_boundaries
                .iter()
                .filter(|b| b.patch == obstacle)
            {
                require(
                    b.kind == PressureBoundaryType::ZeroGradient && b.value == 0.0,
                    "obstacle pressure boundaries must be zero-gradient",
                )?;
            }
        }
        let mut s = CollocatedIncompressible::new(mesh.clone(), self.config)?;
        s.initialize_uniform(self.initial_velocity, self.initial_pressure)?;
        for b in &self.velocity_boundaries {
            s.set_velocity_boundary(&b.patch, b.kind, b.value)?;
        }
        for b in &self.pressure_boundaries {
            s.set_pressure_boundary(&b.patch, b.kind, b.value)?;
        }
        for obstacle in self.solid_names() {
            s.set_velocity_boundary(obstacle, VelocityBoundaryType::FixedValue, Vec3::ZERO)?;
            s.set_pressure_boundary(obstacle, PressureBoundaryType::ZeroGradient, 0.0)?;
        }
        s.set_body_acceleration(&vec![self.acceleration; mesh.cell_count()])?;
        let mut fan = match &self.fan {
            Some(f) => {
                let controller = super::fan::FanController::new(self, f, &mesh)?;
                s.set_velocity_boundary(
                    controller.patch(),
                    VelocityBoundaryType::FixedValue,
                    controller.velocity(),
                )?;
                Some(controller)
            }
            None => None,
        };
        let mut energy = match &self.thermal {
            Some(heat) => Some(super::conjugate::EnergyState::new(self, heat)?),
            None => None,
        };
        let mut monitor = super::goals::GoalMonitor::new(self, &mesh)?;
        let fixed_pressure = super::goals::fixed_pressure_lookup(self, &mesh);
        let stop = std::cell::Cell::new(false);
        let porous = if self.porous_blocks.is_empty() {
            None
        } else {
            Some(super::porous::PorousZones::new(self, &mesh)?)
        };
        let mut last = None;
        let steady = self.method == FlowMethod::Simple;
        let stem = if steady { "steady_iterations" } else { "flow" };
        series_until(
            |i| {
                let info = if i > 0 {
                    if let Some(zones) = &porous {
                        let beta = zones.drag(self, s.velocity())?;
                        s.set_implicit_drag(&beta)?;
                    }
                    let r = match self.method {
                        FlowMethod::Simple => s.iterate_simple()?,
                        FlowMethod::Piso => s.step_piso()?,
                        FlowMethod::Pimple => s.step_pimple()?,
                    };
                    let linear =
                        r.pressure.converged && s.momentum_results().iter().all(|r| r.converged);
                    last = Some((r, linear));
                    if let (Some(state), Some(heat)) = (energy.as_mut(), &self.thermal) {
                        if let Some(b) = &heat.buoyancy {
                            state.step(s.face_flux(), heat.inflow_temperature_k, self.config.dt)?;
                            let acceleration =
                                state.buoyancy(b, self.acceleration, heat.inflow_temperature_k);
                            s.set_body_acceleration(&acceleration)?;
                        }
                    }
                    if let Some(goals) = monitor.as_mut() {
                        if goals.due(i) {
                            if let (Some(state), Some(heat)) = (energy.as_mut(), &self.thermal) {
                                if heat.buoyancy.is_none() && goals.needs_thermal_refresh(i) {
                                    state.steady(s.face_flux(), heat.inflow_temperature_k)?;
                                }
                            }
                            let lookup =
                                |patch: usize| fixed_pressure.get(patch).copied().flatten();
                            if goals.record(
                                i,
                                &mesh,
                                s.face_flux(),
                                s.pressure(),
                                &lookup,
                                energy.as_ref(),
                            )? {
                                goals.stopped_at = Some(i);
                                stop.set(true);
                            }
                        }
                    }
                    if let Some(controller) = fan.as_mut() {
                        if let Some(velocity) = controller.update(i, &mesh, s.pressure())? {
                            let patch = controller.patch().to_owned();
                            s.set_velocity_boundary_value(&patch, velocity)?;
                        }
                    }
                    Some(r)
                } else {
                    None
                };
                let d = velocity_derived_fields(&mesh, s.velocity())?;
                let mut data = vec![
                    DataArray::vector("velocity_m_per_s", s.velocity().to_vec()),
                    DataArray::scalar("kinematic_pressure_m2_per_s2", s.pressure().to_vec()),
                    DataArray::vector("vorticity_per_s", d.vorticity),
                    DataArray::scalar("Q_per_s2", d.q_criterion),
                    DataArray::scalar("lambda2_per_s2", d.lambda2),
                    DataArray::scalar("enstrophy_per_s2", d.enstrophy),
                ];
                if let (Some(state), Some(heat)) = (energy.as_ref(), &self.thermal) {
                    if heat.buoyancy.is_some() {
                        data.push(DataArray::scalar(
                            "temperature_k",
                            state.fluid_temperature(heat.inflow_temperature_k),
                        ));
                    }
                }
                Ok((
                    if steady { i as f64 } else { s.time() },
                    data,
                    serde_json::json!({"physical_time_s":s.time(),"coupling":info,"momentum":s.momentum_results(),"continuity_l2":s.continuity_l2()?,
                        "fan":fan.as_ref().and_then(|f| f.history.last().copied())}),
                    last.map(|(_, linear)| linear),
                ))
            },
            out,
            m,
            &export,
            self.output,
            stem,
            &stop,
        )?;
        if let Some((r, linear)) = last {
            m.metrics
                .insert("velocity_rms_change".into(), r.velocity_rms_change);
            if steady {
                m.steady_state_converged = Some(
                    linear
                        && r.velocity_rms_change <= self.steady_velocity_tolerance
                        && r.continuity_l2 <= self.steady_continuity_tolerance,
                );
            }
        }
        m.metrics.insert("time_s".into(), s.time());
        if self.mesh.xy_shear == 0.0 {
            // The collocated PISO/PIMPLE path is reliable for small viscous diffusion numbers;
            // keep this below about 0.1 (see the README).
            let h = (0..3)
                .map(|a| self.mesh.lengths.component(a) / self.mesh.dimensions[a] as f64)
                .fold(f64::MAX, f64::min);
            m.metrics.insert(
                "max_viscous_diffusion_number".into(),
                self.config.kinematic_viscosity * self.config.dt / (h * h),
            );
        }
        m.metrics.insert("continuity_l2".into(), s.continuity_l2()?);
        m.metrics
            .insert("continuity_max".into(), s.continuity_max()?);
        m.metrics
            .insert("kinetic_energy_J".into(), s.kinetic_energy());
        if !self.obstacles.is_empty() || !self.bodies.is_empty() {
            super::obstacles::write_report(self, &mesh, s.face_flux(), out, m)?;
        }
        if let Some(zones) = &porous {
            zones.write_report(self, &mesh, s.velocity(), out, m)?;
        }
        if let Some(controller) = fan {
            super::fan::write_report(&controller, out, m)?;
        }
        if let Some(goals) = &monitor {
            goals.write(out, m)?;
        }
        if let (Some(mut state), Some(heat)) = (energy, &self.thermal) {
            if heat.buoyancy.is_none() {
                state.steady(s.face_flux(), heat.inflow_temperature_k)?;
            }
            state.finish(self, heat, s.face_flux(), out, m)?;
        }
        Ok(())
    }
}
impl SolidHeatCase {
    pub(super) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        let (mesh, export) = self.mesh.build()?;
        require(
            self.initial_gradient.is_finite(),
            "invalid temperature gradient",
        )?;
        unique(self.boundaries.iter().map(|b| b.patch.as_str()))?;
        let mut s = SolidHeatConduction::new(mesh, self.config)?;
        s.initialize_profile(|p| self.initial + self.initial_gradient.dot(p))?;
        s.set_volumetric_heat_source_constant(self.source)?;
        for b in &self.boundaries {
            s.set_boundary(
                &b.patch,
                if b.kind == ScalarBoundaryType::FixedValue {
                    SolidThermalBoundaryType::FixedTemperature
                } else {
                    SolidThermalBoundaryType::Adiabatic
                },
                b.value,
            )?;
        }
        series(
            |i| {
                let r = if i > 0 { Some(s.step()?) } else { None };
                Ok((
                    s.time(),
                    vec![DataArray::scalar("temperature_K", s.temperature().to_vec())],
                    serde_json::json!({"linear":r,"sensible_energy_J":s.total_sensible_energy()}),
                    r.map(|r| r.converged),
                ))
            },
            out,
            m,
            &export,
            self.output,
            "solid_heat",
        )?;
        m.metrics.insert("time_s".into(), s.time());
        m.metrics
            .insert("sensible_energy_J".into(), s.total_sensible_energy());
        Ok(())
    }
}
impl ThermalTransportCase {
    pub(super) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        let (mesh, export) = self.mesh.build()?;
        require(
            self.velocity.is_finite()
                && self.initial_gradient.is_finite()
                && self.source.is_finite(),
            "invalid thermal inputs",
        )?;
        unique(self.boundaries.iter().map(|b| b.patch.as_str()))?;
        let mut s = ThermalTransport::new(mesh.clone(), self.config)?;
        s.initialize_profile(|p| self.initial + self.initial_gradient.dot(p))?;
        for b in &self.boundaries {
            s.set_boundary(
                &b.patch,
                ThermalBoundary {
                    fixed_temperature: b.kind == ScalarBoundaryType::FixedValue,
                    temperature: b.value,
                },
            )?;
        }
        s.set_face_mass_flux(
            mesh.faces()
                .iter()
                .map(|f| self.config.density * self.velocity.dot(f.area))
                .collect(),
        )?;
        s.set_reactive_radiative_source(
            Some(Box::new(move |_, _, _| Ok(self.source))),
            Arc::new(OpticallyThinRadiation::new(self.absorption)?),
            self.environment_temperature,
        )?;
        if let Some(dt) = s.transport_stable_timestep(1.0)? {
            require(
                self.config.dt <= dt,
                "explicit thermal dt exceeds transport stability bound",
            )?;
            m.metrics.insert("transport_dt_bound_s".into(), dt);
        }
        series(
            |i| {
                if i > 0 {
                    s.step()?;
                }
                Ok((
                    s.time(),
                    vec![DataArray::scalar("temperature_K", s.temperature().to_vec())],
                    serde_json::json!({"sensible_energy_J":s.total_sensible_enthalpy()?}),
                    None,
                ))
            },
            out,
            m,
            &export,
            self.output,
            "thermal",
        )?;
        m.metrics.insert("time_s".into(), s.time());
        m.metrics
            .insert("sensible_energy_J".into(), s.total_sensible_enthalpy()?);
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StressPatch {
    pub patch: String,
    pub kind: ScalarBoundaryType,
    pub value: SymmetricTensor3,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "model", rename_all = "snake_case", deny_unknown_fields)]
pub enum RansModelSpec {
    SpalartAllmaras {
        #[serde(default)]
        config: SpalartAllmarasConfig,
        nu_tilde: f64,
        wall_distance: f64,
        #[serde(default)]
        boundaries: Vec<ScalarPatch>,
    },
    KEpsilon {
        #[serde(default)]
        config: KEpsilonConfig,
        k: f64,
        epsilon: f64,
        #[serde(default)]
        k_boundaries: Vec<ScalarPatch>,
        #[serde(default)]
        epsilon_boundaries: Vec<ScalarPatch>,
    },
    KOmegaSst {
        #[serde(default)]
        config: KOmegaSSTConfig,
        k: f64,
        omega: f64,
        wall_distance: f64,
        grid_scale: Option<f64>,
        #[serde(default)]
        k_boundaries: Vec<ScalarPatch>,
        #[serde(default)]
        omega_boundaries: Vec<ScalarPatch>,
    },
    ReynoldsStress {
        #[serde(default)]
        config: ReynoldsStressConfig,
        stress: SymmetricTensor3,
        epsilon: f64,
        #[serde(default)]
        stress_boundaries: Vec<StressPatch>,
        #[serde(default)]
        epsilon_boundaries: Vec<ScalarPatch>,
    },
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RansTransportCase {
    #[serde(default)]
    pub mesh: PolyMeshSpec,
    #[serde(default)]
    pub output: PolyOutput,
    pub turbulence: RansModelSpec,
    #[serde(default)]
    pub velocity: Vec3,
    #[serde(default)]
    pub velocity_gradient: VelocityGradient3,
}
enum NativeRans {
    Sa(SpalartAllmarasTransport),
    Ke(KEpsilonTransport),
    Sst(KOmegaSSTTransport),
    Rsm(ReynoldsStressTransport),
}
impl RansModelSpec {
    fn build(&self, mesh: PolyMesh) -> Result<NativeRans> {
        Ok(match self {
            Self::SpalartAllmaras {
                config,
                nu_tilde,
                wall_distance,
                boundaries,
            } => {
                let mut s = SpalartAllmarasTransport::new(mesh, *config)?;
                s.initialize(*nu_tilde)?;
                s.set_wall_distance(*wall_distance)?;
                unique(boundaries.iter().map(|b| b.patch.as_str()))?;
                for b in boundaries {
                    s.set_boundary(&b.patch, b.kind, b.value)?;
                }
                NativeRans::Sa(s)
            }
            Self::KEpsilon {
                config,
                k,
                epsilon,
                k_boundaries,
                epsilon_boundaries,
            } => {
                let mut s = KEpsilonTransport::new(mesh, *config)?;
                s.initialize(*k, *epsilon)?;
                unique(k_boundaries.iter().map(|b| b.patch.as_str()))?;
                unique(epsilon_boundaries.iter().map(|b| b.patch.as_str()))?;
                for b in k_boundaries {
                    s.set_k_boundary(&b.patch, b.kind, b.value)?;
                }
                for b in epsilon_boundaries {
                    s.set_epsilon_boundary(&b.patch, b.kind, b.value)?;
                }
                NativeRans::Ke(s)
            }
            Self::KOmegaSst {
                config,
                k,
                omega,
                wall_distance,
                grid_scale,
                k_boundaries,
                omega_boundaries,
            } => {
                let mut s = KOmegaSSTTransport::new(mesh, *config)?;
                s.initialize(*k, *omega)?;
                s.set_wall_distance(*wall_distance)?;
                if let Some(d) = grid_scale {
                    s.set_grid_scale(*d)?;
                }
                unique(k_boundaries.iter().map(|b| b.patch.as_str()))?;
                unique(omega_boundaries.iter().map(|b| b.patch.as_str()))?;
                for b in k_boundaries {
                    s.set_k_boundary(&b.patch, b.kind, b.value)?;
                }
                for b in omega_boundaries {
                    s.set_omega_boundary(&b.patch, b.kind, b.value)?;
                }
                NativeRans::Sst(s)
            }
            Self::ReynoldsStress {
                config,
                stress,
                epsilon,
                stress_boundaries,
                epsilon_boundaries,
            } => {
                let n = mesh.cell_count();
                let mut s = ReynoldsStressTransport::new(mesh, *config)?;
                s.initialize_profile(&vec![*stress; n], &vec![*epsilon; n])?;
                unique(stress_boundaries.iter().map(|b| b.patch.as_str()))?;
                unique(epsilon_boundaries.iter().map(|b| b.patch.as_str()))?;
                for b in stress_boundaries {
                    s.set_boundary(&b.patch, b.kind, b.value)?;
                }
                for b in epsilon_boundaries {
                    s.set_epsilon_boundary(&b.patch, b.kind, b.value)?;
                }
                NativeRans::Rsm(s)
            }
        })
    }
}
impl RansTransportCase {
    pub(super) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        let (mesh, export) = self.mesh.build()?;
        require(
            self.velocity.is_finite() && self.velocity_gradient.iter().all(|v| v.is_finite()),
            "nonfinite RANS velocity",
        )?;
        let velocity = |p: Vec3| {
            let g = self.velocity_gradient;
            self.velocity
                + Vec3::new(
                    g[0] * p.x + g[1] * p.y + g[2] * p.z,
                    g[3] * p.x + g[4] * p.y + g[5] * p.z,
                    g[6] * p.x + g[7] * p.y + g[8] * p.z,
                )
        };
        let u: Vec<_> = mesh.cells().iter().map(|c| velocity(c.center)).collect();
        let flux: Vec<_> = mesh
            .faces()
            .iter()
            .map(|f| velocity(f.center).dot(f.area))
            .collect();
        let mut solver = self.turbulence.build(mesh)?;
        match &mut solver {
            NativeRans::Sa(s) => {
                s.set_velocity(&u)?;
                s.set_face_flux(flux)?;
            }
            NativeRans::Ke(s) => {
                s.set_velocity(&u)?;
                s.set_face_flux(flux)?;
            }
            NativeRans::Sst(s) => {
                s.set_velocity(&u)?;
                s.set_face_flux(flux)?;
            }
            NativeRans::Rsm(s) => {
                s.set_velocity(&u)?;
                s.set_face_flux(flux)?;
            }
        }
        let mut time = 0.0;
        series(
            |i| {
                let (t, mut data, linear) = match &mut solver {
                    NativeRans::Sa(s) => {
                        let r = if i > 0 { Some(s.step()?) } else { None };
                        (
                            s.time(),
                            vec![
                                DataArray::scalar("nu_tilde_m2_per_s", s.nu_tilde().to_vec()),
                                DataArray::scalar(
                                    "eddy_viscosity_m2_per_s",
                                    s.kinematic_eddy_viscosity()?,
                                ),
                            ],
                            r.map(|r| r.converged),
                        )
                    }
                    NativeRans::Ke(s) => {
                        let r = if i > 0 { Some(s.step()?) } else { None };
                        (
                            s.time(),
                            vec![
                                DataArray::scalar("k_m2_per_s2", s.k().to_vec()),
                                DataArray::scalar("epsilon_m2_per_s3", s.epsilon().to_vec()),
                                DataArray::scalar(
                                    "eddy_viscosity_m2_per_s",
                                    s.kinematic_eddy_viscosity()?,
                                ),
                            ],
                            r.map(|r| r.first.converged && r.second.converged),
                        )
                    }
                    NativeRans::Sst(s) => {
                        let r = if i > 0 { Some(s.step()?) } else { None };
                        (
                            s.time(),
                            vec![
                                DataArray::scalar("k_m2_per_s2", s.k().to_vec()),
                                DataArray::scalar("omega_per_s", s.omega().to_vec()),
                                DataArray::scalar(
                                    "eddy_viscosity_m2_per_s",
                                    s.kinematic_eddy_viscosity()?,
                                ),
                                DataArray::scalar("SST_F1", s.blending_f1().to_vec()),
                                DataArray::scalar("SST_F2", s.blending_f2().to_vec()),
                                DataArray::scalar(
                                    "DES_dissipation_factor",
                                    s.hybrid_dissipation_factor().to_vec(),
                                ),
                            ],
                            r.map(|r| r.first.converged && r.second.converged),
                        )
                    }
                    NativeRans::Rsm(s) => {
                        let r = if i > 0 { Some(s.step()?) } else { None };
                        let stress = s.reynolds_stress();
                        let data = vec![
                            DataArray {
                                name: "Reynolds_stress_xx_yy_zz_xy_xz_yz_m2_per_s2".into(),
                                components: 6,
                                values: stress.iter().flat_map(|r| r.components()).collect(),
                            },
                            DataArray::scalar("k_m2_per_s2", s.turbulent_kinetic_energy()),
                            DataArray::scalar("epsilon_m2_per_s3", s.epsilon().to_vec()),
                            DataArray::scalar(
                                "eddy_viscosity_m2_per_s",
                                s.kinematic_eddy_viscosity()?,
                            ),
                        ];
                        (
                            s.time(),
                            data,
                            r.map(|r| r.stress.iter().all(|r| r.converged) && r.epsilon.converged),
                        )
                    }
                };
                time = t;
                data.push(DataArray::vector("prescribed_velocity_m_per_s", u.clone()));
                Ok((
                    t,
                    data,
                    serde_json::json!({"linear_converged":linear,"velocity_coupling":"prescribed"}),
                    linear,
                ))
            },
            out,
            m,
            &export,
            self.output,
            "rans",
        )?;
        m.metrics.insert("time_s".into(), time);
        Ok(())
    }
}
