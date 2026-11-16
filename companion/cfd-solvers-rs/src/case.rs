//! Typed, intentionally bounded case runner. This is a Rust API, not a parser for the upstream fixed-reference dispatcher.
use crate::{
    acoustics::*,
    chemistry::*,
    core::{LinearConfig, Vec3},
    fdtd::*,
    fem::*,
    fvm::*,
    io::{self, DataArray, VtuMesh},
    lbm::*,
    optics::*,
    require,
    rf::*,
    Result, UPSTREAM_REVISION,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, path::Path};
mod bodies;
pub mod conjugate;
pub mod domains;
mod fan;
pub mod goals;
mod obstacles;
mod poly;
mod porous;
pub mod verification;
pub use bodies::SurfaceBody;
pub use domains::{electrical::*, mechanics::*, physics::*};
pub use fan::{FanCurvePoint, FanInlet, FanPressureConvention, FanSample};
pub use poly::*;
pub use porous::PorousBlock;
fn one() -> f64 {
    1.0
}
fn half() -> f64 {
    0.5
}
fn elements() -> usize {
    32
}
fn steps() -> usize {
    100
}
fn interval() -> usize {
    10
}
fn duration() -> f64 {
    0.1
}
fn amplitude() -> f32 {
    0.03
}
fn initial_density() -> f32 {
    1.0
}
fn gaussian_center() -> f64 {
    0.3
}
fn gaussian_width() -> f64 {
    0.04
}
fn reference() -> f64 {
    50.0
}
fn temperature() -> f64 {
    300.0
}
fn sound_speed() -> f64 {
    1500.0
}
fn mass_density() -> f64 {
    1000.0
}
fn sigma() -> f64 {
    0.0003
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Mesh2DSpec {
    Rectangle {
        nx: usize,
        ny: usize,
        #[serde(default = "one")]
        width: f64,
        #[serde(default = "one")]
        height: f64,
    },
    Explicit {
        mesh: Mesh2D,
    },
}
impl Default for Mesh2DSpec {
    fn default() -> Self {
        Self::Rectangle {
            nx: 16,
            ny: 16,
            width: 1.0,
            height: 1.0,
        }
    }
}
impl Mesh2DSpec {
    pub fn build(&self) -> Result<Mesh2D> {
        match self {
            Self::Rectangle {
                nx,
                ny,
                width,
                height,
            } => make_rectangle_tri_mesh(*nx, *ny, *width, *height),
            Self::Explicit { mesh } => {
                mesh.validate()?;
                Ok(mesh.clone())
            }
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Mesh3DSpec {
    Box {
        nx: usize,
        ny: usize,
        nz: usize,
        #[serde(default = "one")]
        width: f64,
        #[serde(default = "one")]
        height: f64,
        #[serde(default = "one")]
        depth: f64,
    },
    Explicit {
        mesh: Mesh3D,
    },
}
impl Default for Mesh3DSpec {
    fn default() -> Self {
        Self::Box {
            nx: 4,
            ny: 4,
            nz: 4,
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
    }
}
impl Mesh3DSpec {
    pub fn build(&self) -> Result<Mesh3D> {
        match self {
            Self::Box {
                nx,
                ny,
                nz,
                width,
                height,
                depth,
            } => make_box_tet_mesh(*nx, *ny, *nz, *width, *height, *depth),
            Self::Explicit { mesh } => {
                mesh.validate()?;
                Ok(mesh.clone())
            }
        }
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum LbmInitial {
    Uniform {
        #[serde(default = "initial_density")]
        rho: f32,
        #[serde(default)]
        ux: f32,
        #[serde(default)]
        uy: f32,
    },
    TaylorGreen {
        #[serde(default = "amplitude")]
        amplitude: f32,
    },
}
impl Default for LbmInitial {
    fn default() -> Self {
        Self::TaylorGreen { amplitude: 0.03 }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum ConstantBoundary {
    Dirichlet { value: f64 },
    Neumann { value: f64 },
    Robin { value: f64, alpha: f64 },
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PatchBoundary {
    pub patch: i32,
    #[serde(flatten)]
    pub condition: ConstantBoundary,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "solver", rename_all = "snake_case", deny_unknown_fields)]
pub enum Case {
    PolyScalarTransport3d(PolyScalarCase),
    CollocatedFlow3d(CollocatedFlowCase),
    SolidHeat3d(SolidHeatCase),
    ThermalTransport3d(ThermalTransportCase),
    RansTransport3d(RansTransportCase),
    Poisson1d {
        #[serde(default = "elements")]
        elements: usize,
        #[serde(default = "one")]
        length: f64,
        #[serde(default = "one")]
        source: f64,
        #[serde(default)]
        left: f64,
        #[serde(default)]
        right: f64,
    },
    Poisson2d {
        #[serde(default)]
        mesh: Mesh2DSpec,
        #[serde(default = "one")]
        diffusivity: f64,
        #[serde(default = "one")]
        source: f64,
        #[serde(default)]
        boundary: f64,
    },
    Poisson3d {
        #[serde(default)]
        mesh: Mesh3DSpec,
        #[serde(default = "one")]
        diffusivity: f64,
        #[serde(default = "one")]
        source: f64,
        #[serde(default)]
        boundary: f64,
    },
    ScalarDiffusion2d {
        #[serde(default)]
        mesh: Mesh2DSpec,
        #[serde(default = "one")]
        diffusivity: f64,
        #[serde(default)]
        reaction: f64,
        #[serde(default = "one")]
        source: f64,
        boundaries: Vec<PatchBoundary>,
    },
    Heat2d {
        #[serde(default)]
        mesh: Mesh2DSpec,
        #[serde(default)]
        config: Heat2DConfig,
        #[serde(default = "steps")]
        steps: usize,
        #[serde(default = "interval")]
        output_interval: usize,
        #[serde(default = "one")]
        initial: f64,
        #[serde(default)]
        boundary: f64,
        #[serde(default)]
        source: f64,
    },
    LbmD2q9 {
        #[serde(default)]
        config: D2Q9Config,
        #[serde(default = "steps")]
        steps: usize,
        #[serde(default)]
        initial: LbmInitial,
        #[serde(default)]
        walls_y: bool,
        #[serde(default)]
        top_wall_speed: f32,
    },
    Maxwell1d {
        #[serde(default)]
        config: Maxwell1DConfig,
        #[serde(default = "steps")]
        steps: usize,
        #[serde(default = "gaussian_center")]
        center_fraction: f64,
        #[serde(default = "gaussian_width")]
        width_fraction: f64,
    },
    Euler1d {
        #[serde(default)]
        config: Compressible1DConfig,
        #[serde(default = "duration")]
        duration: f64,
        left: EulerPrimitive1D,
        right: EulerPrimitive1D,
        #[serde(default = "half")]
        interface_fraction: f64,
    },
    Diffusion2d {
        #[serde(default)]
        config: Diffusion2DConfig,
        #[serde(default = "steps")]
        iterations: usize,
        #[serde(default = "one")]
        source: f64,
        #[serde(default)]
        boundary: [f64; 4],
    },
    RfTem {
        line: TemTransmissionLine,
        frequencies_hz: Vec<f64>,
        #[serde(default = "reference")]
        reference_impedance: f64,
    },
    OpticsSequential {
        #[serde(default = "one")]
        object_space_index: f64,
        surfaces: Vec<SequentialSurface>,
        rays: Vec<Ray>,
        image_z: Option<f64>,
    },
    ChemistryBatch {
        species: Vec<Species>,
        reactions: Vec<ElementaryReaction>,
        concentrations: Vec<f64>,
        #[serde(default = "temperature")]
        temperature: f64,
        duration: f64,
        #[serde(default)]
        config: ReactorConfig,
    },
    SpiceCircuit(SpiceCircuitCase),
    TcadDiode1d(TcadDiodeCase),
    BatteryCell(BatteryCellCase),
    WireAntenna(WireAntennaCase),
    Maxwell3d(Maxwell3dCase),
    EdgeCavity3d(EdgeCavityCase),
    PdnIrDrop(PdnIrDropCase),
    GalvanicCorrosion(GalvanicCorrosionCase),
    Elasticity2d(Elasticity2dCase),
    BarModes1d(BarModes1dCase),
    J2Uniaxial(J2UniaxialCase),
    PhaseFieldFracture1d(PhaseFieldFracture1dCase),
    CosseratCantilever(CosseratCantileverCase),
    DemGranular(DemGranularCase),
    PicTwoStream1d(PicTwoStreamCase),
    ParticleTransportSlab(ParticleTransportCase),
    PennesBioheat2d(BioheatCase),
    JouleHeating2d(JouleHeatingCase),
    CahnHilliard1d(CahnHilliardCase),
    DgAdvection1d(DgAdvectionCase),
    LandauDamping1d(LandauDampingCase),
    LennardJonesMd1d(LennardJonesMdCase),
    PhotonicWaveguide(PhotonicWaveguideCase),
    ThmcColumn1d(ThmcColumnCase),
    PowerNetwork(PowerNetworkCase),
    MultiphaseBox(MultiphaseBoxCase),
    SprayInjection(SprayInjectionCase),
    AcousticKspace2d {
        #[serde(default)]
        config: KSpace2DConfig,
        #[serde(default = "sound_speed")]
        sound_speed: f64,
        #[serde(default = "mass_density")]
        density: f64,
        #[serde(default = "steps")]
        steps: usize,
        #[serde(default = "half")]
        center_x_fraction: f64,
        #[serde(default = "half")]
        center_y_fraction: f64,
        #[serde(default = "sigma")]
        sigma: f64,
        #[serde(default = "one")]
        amplitude: f64,
        #[serde(default)]
        sensor_cells: Vec<usize>,
    },
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RunManifest {
    pub schema: String,
    pub rust_port_version: String,
    pub upstream_revision: String,
    pub solver: String,
    pub backend: String,
    /// Supported protocol extensions, distinct from the upstream source baseline.
    #[serde(default)]
    pub capabilities: Vec<String>,
    pub complete_clone: bool,
    pub execution_completed: bool,
    pub linear_solver_converged: Option<bool>,
    #[serde(default)]
    pub steady_state_converged: Option<bool>,
    pub metrics: BTreeMap<String, f64>,
    pub files: Vec<String>,
    pub case: Case,
}
impl Case {
    pub fn name(&self) -> &'static str {
        match self {
            Self::PolyScalarTransport3d(..) => "poly_scalar_transport3d",
            Self::CollocatedFlow3d(..) => "collocated_flow3d",
            Self::SolidHeat3d(..) => "solid_heat3d",
            Self::ThermalTransport3d(..) => "thermal_transport3d",
            Self::RansTransport3d(..) => "rans_transport3d",
            Self::Poisson1d { .. } => "poisson1d",
            Self::Poisson2d { .. } => "poisson2d",
            Self::Poisson3d { .. } => "poisson3d",
            Self::ScalarDiffusion2d { .. } => "scalar_diffusion2d",
            Self::Heat2d { .. } => "heat2d",
            Self::LbmD2q9 { .. } => "lbm_d2q9",
            Self::Maxwell1d { .. } => "maxwell1d",
            Self::Euler1d { .. } => "euler1d",
            Self::Diffusion2d { .. } => "diffusion2d",
            Self::RfTem { .. } => "rf_tem",
            Self::OpticsSequential { .. } => "optics_sequential",
            Self::ChemistryBatch { .. } => "chemistry_batch",
            Self::AcousticKspace2d { .. } => "acoustic_kspace2d",
            Self::SpiceCircuit(..) => "spice_circuit",
            Self::TcadDiode1d(..) => "tcad_diode1d",
            Self::BatteryCell(..) => "battery_cell",
            Self::WireAntenna(..) => "wire_antenna",
            Self::Maxwell3d(..) => "maxwell3d",
            Self::EdgeCavity3d(..) => "edge_cavity3d",
            Self::PdnIrDrop(..) => "pdn_ir_drop",
            Self::GalvanicCorrosion(..) => "galvanic_corrosion",
            Self::Elasticity2d(..) => "elasticity2d",
            Self::BarModes1d(..) => "bar_modes1d",
            Self::J2Uniaxial(..) => "j2_uniaxial",
            Self::PhaseFieldFracture1d(..) => "phase_field_fracture1d",
            Self::CosseratCantilever(..) => "cosserat_cantilever",
            Self::DemGranular(..) => "dem_granular",
            Self::PicTwoStream1d(..) => "pic_two_stream1d",
            Self::ParticleTransportSlab(..) => "particle_transport_slab",
            Self::PennesBioheat2d(..) => "pennes_bioheat2d",
            Self::JouleHeating2d(..) => "joule_heating2d",
            Self::CahnHilliard1d(..) => "cahn_hilliard1d",
            Self::DgAdvection1d(..) => "dg_advection1d",
            Self::LandauDamping1d(..) => "landau_damping1d",
            Self::LennardJonesMd1d(..) => "lennard_jones_md1d",
            Self::PhotonicWaveguide(..) => "photonic_waveguide",
            Self::ThmcColumn1d(..) => "thmc_column1d",
            Self::PowerNetwork(..) => "power_network",
            Self::MultiphaseBox(..) => "multiphase_box",
            Self::SprayInjection(..) => "spray_injection",
        }
    }
    pub fn run(&self, out: &Path) -> Result<RunManifest> {
        require(
            !out.join("manifest.json").exists(),
            "output folder contains a previous run; choose a new folder",
        )?;
        std::fs::create_dir_all(out)?;
        let start = std::time::Instant::now();
        let mut m = RunManifest {
            schema: "cfd-solvers-rs/run/1".into(),
            rust_port_version: env!("CARGO_PKG_VERSION").into(),
            upstream_revision: UPSTREAM_REVISION.into(),
            solver: self.name().into(),
            backend: "native_rust_serial_cpu".into(),
            capabilities: vec![
                "grid_obstacles_v1".into(),
                "fan_curve_v1".into(),
                "porous_blocks_v1".into(),
                "conjugate_heat_v1".into(),
                "surface_bodies_v1".into(),
                "goals_v1".into(),
                "natural_convection_v1".into(),
                "domain_cases_v1".into(),
            ],
            complete_clone: false,
            execution_completed: false,
            linear_solver_converged: None,
            steady_state_converged: None,
            metrics: BTreeMap::new(),
            files: vec![],
            case: self.clone(),
        };
        let linear = LinearConfig {
            max_iterations: 4000,
            relative_tolerance: 1e-11,
            ..Default::default()
        };
        match self {
            Self::PolyScalarTransport3d(c) => c.run(out, &mut m)?,
            Self::CollocatedFlow3d(c) => c.run(out, &mut m)?,
            Self::SolidHeat3d(c) => c.run(out, &mut m)?,
            Self::ThermalTransport3d(c) => c.run(out, &mut m)?,
            Self::RansTransport3d(c) => c.run(out, &mut m)?,
            Self::SpiceCircuit(c) => c.run(out, &mut m)?,
            Self::TcadDiode1d(c) => c.run(out, &mut m)?,
            Self::BatteryCell(c) => c.run(out, &mut m)?,
            Self::WireAntenna(c) => c.run(out, &mut m)?,
            Self::Maxwell3d(c) => c.run(out, &mut m)?,
            Self::EdgeCavity3d(c) => c.run(out, &mut m)?,
            Self::PdnIrDrop(c) => c.run(out, &mut m)?,
            Self::GalvanicCorrosion(c) => c.run(out, &mut m)?,
            Self::Elasticity2d(c) => c.run(out, &mut m)?,
            Self::BarModes1d(c) => c.run(out, &mut m)?,
            Self::J2Uniaxial(c) => c.run(out, &mut m)?,
            Self::PhaseFieldFracture1d(c) => c.run(out, &mut m)?,
            Self::CosseratCantilever(c) => c.run(out, &mut m)?,
            Self::DemGranular(c) => c.run(out, &mut m)?,
            Self::PicTwoStream1d(c) => c.run(out, &mut m)?,
            Self::ParticleTransportSlab(c) => c.run(out, &mut m)?,
            Self::PennesBioheat2d(c) => c.run(out, &mut m)?,
            Self::JouleHeating2d(c) => c.run(out, &mut m)?,
            Self::CahnHilliard1d(c) => c.run(out, &mut m)?,
            Self::DgAdvection1d(c) => c.run(out, &mut m)?,
            Self::LandauDamping1d(c) => c.run(out, &mut m)?,
            Self::LennardJonesMd1d(c) => c.run(out, &mut m)?,
            Self::PhotonicWaveguide(c) => c.run(out, &mut m)?,
            Self::ThmcColumn1d(c) => c.run(out, &mut m)?,
            Self::PowerNetwork(c) => c.run(out, &mut m)?,
            Self::MultiphaseBox(c) => c.run(out, &mut m)?,
            Self::SprayInjection(c) => c.run(out, &mut m)?,
            Self::Poisson1d {
                elements,
                length,
                source,
                left,
                right,
            } => {
                require(source.is_finite(), "source must be finite")?;
                let mut s = Poisson1D::new(*elements, *length)?;
                s.solve(|_| *source, *left, *right)?;
                let mesh = VtuMesh::line(s.x.len(), length / (*elements as f64), 0.0)?;
                write_field(
                    out,
                    &mut m,
                    "solution.vtu",
                    &mesh,
                    vec![DataArray::scalar("potential", s.u.clone())],
                    vec![],
                )?;
                io::write_csv(
                    out.join("solution.csv"),
                    &[("x_m", &s.x), ("potential", &s.u)],
                )?;
                m.files.push("solution.csv".into());
            }
            Self::Poisson2d {
                mesh,
                diffusivity,
                source,
                boundary,
            } => {
                let mesh = mesh.build()?;
                let s = poisson2d(&mesh, *diffusivity, |_| *source, |_| *boundary, linear)?;
                m.linear_solver_converged = Some(s.linear.converged);
                m.metrics
                    .insert("residual_rms".into(), s.linear.residual_rms);
                write_field(
                    out,
                    &mut m,
                    "solution.vtu",
                    &VtuMesh::from_triangles(&mesh)?,
                    vec![DataArray::scalar("potential", s.values)],
                    vec![],
                )?;
            }
            Self::Poisson3d {
                mesh,
                diffusivity,
                source,
                boundary,
            } => {
                let mesh = mesh.build()?;
                let s = poisson3d(&mesh, *diffusivity, |_| *source, |_| *boundary, linear)?;
                m.linear_solver_converged = Some(s.linear.converged);
                m.metrics
                    .insert("residual_rms".into(), s.linear.residual_rms);
                write_field(
                    out,
                    &mut m,
                    "solution.vtu",
                    &VtuMesh::from_tetrahedra(&mesh)?,
                    vec![DataArray::scalar("potential", s.values)],
                    vec![],
                )?;
            }
            Self::ScalarDiffusion2d {
                mesh,
                diffusivity,
                reaction,
                source,
                boundaries,
            } => {
                let mesh = mesh.build()?;
                let mut b = BTreeMap::new();
                for bc in boundaries {
                    let condition = match bc.condition {
                        ConstantBoundary::Dirichlet { value } => {
                            ScalarBoundaryCondition2D::Dirichlet(Box::new(move |_| value))
                        }
                        ConstantBoundary::Neumann { value } => {
                            ScalarBoundaryCondition2D::Neumann(Box::new(move |_| value))
                        }
                        ConstantBoundary::Robin { value, alpha } => {
                            ScalarBoundaryCondition2D::Robin {
                                value: Box::new(move |_| value),
                                alpha: Box::new(move |_| alpha),
                            }
                        }
                    };
                    require(
                        b.insert(bc.patch, condition).is_none(),
                        "duplicate boundary patch",
                    )?;
                }
                let s = scalar_diffusion2d(
                    &mesh,
                    |_| *diffusivity,
                    |_| *reaction,
                    |_| *source,
                    &b,
                    linear,
                )?;
                m.linear_solver_converged = Some(s.linear.converged);
                m.metrics
                    .insert("residual_rms".into(), s.linear.residual_rms);
                write_field(
                    out,
                    &mut m,
                    "solution.vtu",
                    &VtuMesh::from_triangles(&mesh)?,
                    vec![DataArray::scalar("potential", s.values)],
                    vec![],
                )?;
            }
            Self::Heat2d {
                mesh,
                config,
                steps,
                output_interval,
                initial,
                boundary,
                source,
            } => {
                require(*output_interval > 0, "output interval must be positive")?;
                let mut s = Heat2D::new(mesh.build()?, *config)?;
                s.initialize(|_| *initial)?;
                let mesh = VtuMesh::from_triangles(&s.mesh)?;
                let mut frames = vec![];
                write_field(
                    out,
                    &mut m,
                    "heat_000000.vtu",
                    &mesh,
                    vec![DataArray::scalar("temperature", s.temperature.clone())],
                    vec![],
                )?;
                frames.push((0.0, "heat_000000.vtu".into()));
                for i in 1..=*steps {
                    s.step(|_, _, _| *source, |_, _| *boundary)?;
                    if i % output_interval == 0 || i == *steps {
                        let file = format!("heat_{i:06}.vtu");
                        write_field(
                            out,
                            &mut m,
                            &file,
                            &mesh,
                            vec![DataArray::scalar("temperature", s.temperature.clone())],
                            vec![],
                        )?;
                        frames.push((s.time, file));
                    }
                }
                io::write_pvd(out.join("heat.pvd"), &frames)?;
                m.files.push("heat.pvd".into());
                m.metrics.insert("time_s".into(), s.time);
                if let Some(r) = s.linear_result {
                    m.linear_solver_converged = Some(r.converged);
                    m.metrics
                        .insert("last_linear_residual_rms".into(), r.residual_rms);
                }
            }
            Self::LbmD2q9 {
                config,
                steps,
                initial,
                walls_y,
                top_wall_speed,
            } => {
                require(
                    *walls_y || *top_wall_speed == 0.0,
                    "top wall speed requires walls_y",
                )?;
                let mut s = D2Q9Solver::new(*config)?;
                if *walls_y {
                    for x in 0..config.nx {
                        s.set_solid(x, 0, true)?;
                        s.set_wall_velocity(x, config.ny - 1, *top_wall_speed, 0.0)?;
                    }
                }
                match initial {
                    LbmInitial::Uniform { rho, ux, uy } => s.initialize_uniform(*rho, *ux, *uy)?,
                    LbmInitial::TaylorGreen { amplitude } => {
                        s.initialize_taylor_green(*amplitude)?
                    }
                };
                let initial_mass = s.mass();
                s.step(*steps)?;
                if cfg!(feature = "parallel") && config.nx * config.ny >= 8192 {
                    m.backend = "native_rust_rayon_cpu".into();
                }
                m.metrics
                    .insert("initial_mass_lattice".into(), initial_mass);
                m.metrics.insert("mass_lattice".into(), s.mass());
                m.metrics
                    .insert("kinetic_energy_lattice".into(), s.kinetic_energy());
                m.metrics
                    .insert("max_speed_lattice".into(), s.max_speed() as f64);
                let velocity = (0..s.density().len())
                    .map(|i| Vec3::new(s.velocity_x()[i] as f64, s.velocity_y()[i] as f64, 0.0))
                    .collect();
                write_field(
                    out,
                    &mut m,
                    "flow.vtu",
                    &VtuMesh::rectangle_cells(config.nx, config.ny, 1.0, 1.0)?,
                    vec![],
                    vec![
                        DataArray::scalar(
                            "density_lattice",
                            s.density().iter().map(|&x| x as f64).collect(),
                        ),
                        DataArray::vector("velocity_lattice", velocity),
                        DataArray::scalar(
                            "solid",
                            s.solid_mask()
                                .iter()
                                .map(|v| if *v { 1.0 } else { 0.0 })
                                .collect(),
                        ),
                    ],
                )?;
            }
            Self::Maxwell1d {
                config,
                steps,
                center_fraction,
                width_fraction,
            } => {
                let mut s = Maxwell1D::new(*config)?;
                s.initialize_gaussian(*center_fraction, *width_fraction)?;
                let initial = s.energy();
                s.step(*steps)?;
                m.metrics
                    .insert("initial_field_energy_per_area".into(), initial);
                m.metrics.insert("field_energy_per_area".into(), s.energy());
                m.metrics.insert("time_s".into(), s.time());
                write_field(
                    out,
                    &mut m,
                    "electric.vtu",
                    &VtuMesh::line(config.cells, config.dx, 0.0)?,
                    vec![
                        DataArray::scalar("electric_V_per_m", s.electric().to_vec()),
                        DataArray::scalar("polarization_C_per_m2", s.polarization().to_vec()),
                    ],
                    vec![],
                )?;
                write_field(
                    out,
                    &mut m,
                    "magnetic.vtu",
                    &VtuMesh::line(config.cells - 1, config.dx, 0.5 * config.dx)?,
                    vec![DataArray::scalar("magnetic_A_per_m", s.magnetic().to_vec())],
                    vec![],
                )?;
                let xe: Vec<_> = (0..config.cells).map(|i| i as f64 * config.dx).collect();
                let xh: Vec<_> = (0..config.cells - 1)
                    .map(|i| (i as f64 + 0.5) * config.dx)
                    .collect();
                io::write_csv(
                    out.join("electric.csv"),
                    &[("x_m", &xe), ("E_V_per_m", s.electric())],
                )?;
                io::write_csv(
                    out.join("magnetic.csv"),
                    &[("x_m", &xh), ("H_A_per_m", s.magnetic())],
                )?;
                m.files
                    .extend(["electric.csv".into(), "magnetic.csv".into()]);
            }
            Self::Euler1d {
                config,
                duration,
                left,
                right,
                interface_fraction,
            } => {
                let mut s = CompressibleEuler1D::new(*config)?;
                s.initialize(*left, *right, *interface_fraction)?;
                s.run(*duration)?;
                let dx = config.length / config.cells as f64;
                let rho: Vec<_> = s.conserved().iter().map(|u| u.density).collect();
                let momentum: Vec<_> = s.conserved().iter().map(|u| u.momentum).collect();
                let energy: Vec<_> = s.conserved().iter().map(|u| u.energy).collect();
                let pressure = (0..config.cells)
                    .map(|i| Ok(s.primitive(i)?.pressure))
                    .collect::<Result<Vec<_>>>()?;
                let x: Vec<_> = (0..config.cells).map(|i| (i as f64 + 0.5) * dx).collect();
                write_field(
                    out,
                    &mut m,
                    "flow.vtu",
                    &VtuMesh::line(config.cells + 1, dx, 0.0)?,
                    vec![],
                    vec![
                        DataArray::scalar("density", rho.clone()),
                        DataArray::scalar("momentum", momentum.clone()),
                        DataArray::scalar("total_energy_density", energy),
                        DataArray::scalar("pressure", pressure.clone()),
                    ],
                )?;
                io::write_csv(
                    out.join("flow.csv"),
                    &[
                        ("x", &x),
                        ("density", &rho),
                        ("momentum", &momentum),
                        ("pressure", &pressure),
                    ],
                )?;
                m.files.push("flow.csv".into());
                m.metrics.insert("time".into(), s.time);
                m.metrics.insert("mass".into(), s.mass());
                m.metrics.insert("momentum".into(), s.momentum());
                m.metrics.insert("total_energy".into(), s.total_energy());
            }
            Self::Diffusion2d {
                config,
                iterations,
                source,
                boundary,
            } => {
                let mut s = Diffusion2D::new(*config)?;
                s.set_source(*source)?;
                s.set_dirichlet(boundary[0], boundary[1], boundary[2], boundary[3])?;
                let delta = s.iterate(*iterations)?;
                m.metrics.insert("last_iteration_max_delta".into(), delta);
                let mesh = VtuMesh::rectangle_cells(
                    config.nx - 1,
                    config.ny - 1,
                    config.lx / (config.nx - 1) as f64,
                    config.ly / (config.ny - 1) as f64,
                )?;
                write_field(
                    out,
                    &mut m,
                    "solution.vtu",
                    &mesh,
                    vec![DataArray::scalar("potential", s.field().to_vec())],
                    vec![],
                )?;
            }
            Self::RfTem {
                line,
                frequencies_hz,
                reference_impedance,
            } => {
                require(!frequencies_hz.is_empty(), "RF sweep needs frequencies")?;
                let points = frequencies_hz
                    .iter()
                    .map(|&f| {
                        Ok(Touchstone2PortPoint {
                            frequency_hz: f,
                            s: line.s_parameters(f, *reference_impedance)?,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                std::fs::write(
                    out.join("network.s2p"),
                    write_touchstone_s2p(&points, *reference_impedance)?,
                )?;
                m.files.push("network.s2p".into());
                let re: Vec<_> = points.iter().map(|p| p.s.a21.re).collect();
                let im: Vec<_> = points.iter().map(|p| p.s.a21.im).collect();
                io::write_csv(
                    out.join("s21.csv"),
                    &[
                        ("frequency_Hz", frequencies_hz),
                        ("S21_real", &re),
                        ("S21_imag", &im),
                    ],
                )?;
                m.files.push("s21.csv".into());
            }
            Self::OpticsSequential {
                object_space_index,
                surfaces,
                rays,
                image_z,
            } => {
                require(
                    !surfaces.is_empty() && !rays.is_empty(),
                    "optical trace needs surfaces and rays",
                )?;
                let mut s = SequentialOpticalSystem::new(*object_space_index)?;
                for surface in surfaces {
                    s.add_surface(*surface)?;
                }
                let traced = s.trace_many(rays)?;
                std::fs::write(out.join("rays.json"), serde_json::to_vec_pretty(&traced)?)?;
                m.files.push("rays.json".into());
                m.metrics.insert(
                    "valid_rays".into(),
                    traced.iter().filter(|r| r.valid).count() as f64,
                );
                if let Some(z) = image_z {
                    m.metrics.insert(
                        "rms_spot_radius_m".into(),
                        rms_spot_radius_at_plane(&traced, *z)?,
                    );
                }
            }
            Self::ChemistryBatch {
                species,
                reactions,
                concentrations,
                temperature,
                duration,
                config,
            } => {
                let mut s = ReactionNetwork::new(species.clone());
                for r in reactions {
                    s.add_reaction(r.clone())?;
                }
                let mut c = concentrations.clone();
                let r = integrate_isothermal(&s, &mut c, *temperature, *duration, *config)?;
                std::fs::write(
                    out.join("concentrations.json"),
                    serde_json::to_vec_pretty(
                        &serde_json::json!({"species":species,"initial":concentrations,"final":c,"temperature_K":temperature,"time_s":r.time_advanced}),
                    )?,
                )?;
                m.files.push("concentrations.json".into());
                m.metrics.insert("time_s".into(), r.time_advanced);
                m.metrics
                    .insert("accepted_steps".into(), r.accepted_steps as f64);
                m.metrics
                    .insert("rejected_steps".into(), r.rejected_steps as f64);
            }
            Self::AcousticKspace2d {
                config,
                sound_speed,
                density,
                steps,
                center_x_fraction,
                center_y_fraction,
                sigma,
                amplitude,
                sensor_cells,
            } => {
                let mut s = KSpaceAcoustic2D::new(*config)?;
                s.set_uniform_medium(*sound_speed, *density)?;
                s.add_sensor_array(sensor_cells)?;
                s.initialize_gaussian_pressure(
                    *center_x_fraction,
                    *center_y_fraction,
                    *sigma,
                    *amplitude,
                )?;
                let initial = s.total_energy();
                s.step(*steps)?;
                write_field(
                    out,
                    &mut m,
                    "acoustic.vtu",
                    &VtuMesh::rectangle_cells(config.nx, config.ny, config.dx, config.dy)?,
                    vec![],
                    vec![
                        DataArray::scalar("pressure_Pa", s.pressure().to_vec()),
                        DataArray::vector(
                            "velocity_m_per_s",
                            (0..s.cells())
                                .map(|q| Vec3::new(s.velocity_x()[q], s.velocity_y()[q], 0.0))
                                .collect(),
                        ),
                    ],
                )?;
                m.metrics.insert("time_s".into(), s.time());
                m.metrics.insert("initial_energy_per_depth".into(), initial);
                m.metrics
                    .insert("energy_per_depth".into(), s.total_energy());
                let time: Vec<_> = (0..=*steps).map(|i| i as f64 * s.dt()).collect();
                for (i, sensor) in s.sensors().iter().enumerate() {
                    let file = format!("sensor_{i:03}.csv");
                    io::write_csv(
                        out.join(&file),
                        &[("time_s", &time), ("pressure_Pa", &sensor.pressure)],
                    )?;
                    m.files.push(file);
                }
            }
        }
        m.metrics
            .insert("wall_time_s".into(), start.elapsed().as_secs_f64());
        require(
            m.metrics.values().all(|v| v.is_finite()),
            "nonfinite run metric",
        )?;
        m.execution_completed = true;
        m.files.push("manifest.json".into());
        std::fs::write(out.join("manifest.json"), serde_json::to_vec_pretty(&m)?)?;
        Ok(m)
    }
}
pub(crate) fn write_field(
    out: &Path,
    m: &mut RunManifest,
    file: &str,
    mesh: &VtuMesh,
    points: Vec<DataArray>,
    cells: Vec<DataArray>,
) -> Result<()> {
    io::write_vtu(out.join(file), mesh, &points, &cells)?;
    m.files.push(file.into());
    Ok(())
}
pub const SOLVERS: [&str; 45] = [
    "poly_scalar_transport3d",
    "collocated_flow3d",
    "solid_heat3d",
    "thermal_transport3d",
    "rans_transport3d",
    "poisson1d",
    "poisson2d",
    "poisson3d",
    "scalar_diffusion2d",
    "heat2d",
    "lbm_d2q9",
    "maxwell1d",
    "euler1d",
    "diffusion2d",
    "rf_tem",
    "optics_sequential",
    "chemistry_batch",
    "acoustic_kspace2d",
    "spice_circuit",
    "tcad_diode1d",
    "battery_cell",
    "wire_antenna",
    "maxwell3d",
    "edge_cavity3d",
    "pdn_ir_drop",
    "galvanic_corrosion",
    "elasticity2d",
    "bar_modes1d",
    "j2_uniaxial",
    "phase_field_fracture1d",
    "cosserat_cantilever",
    "dem_granular",
    "pic_two_stream1d",
    "particle_transport_slab",
    "pennes_bioheat2d",
    "joule_heating2d",
    "cahn_hilliard1d",
    "dg_advection1d",
    "landau_damping1d",
    "lennard_jones_md1d",
    "photonic_waveguide",
    "thmc_column1d",
    "power_network",
    "multiphase_box",
    "spray_injection",
];
