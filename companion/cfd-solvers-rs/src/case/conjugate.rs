//! Steady conjugate heat transfer on the grid-obstacle domain (`conjugate_heat_v1`).
//!
//! After the flow run, one steady energy equation is solved on every background cell:
//! fluid cells conduct and carry heat with the solver's conservative face fluxes (first-order
//! upwind); obstacle cells are solids with grid-aligned orthotropic conductivity and a
//! uniformly distributed power. Fluid-solid interfaces use the series (harmonic) face
//! conductance, so temperature and normal heat flux are continuous by construction.
//! The coupling is one-way: constant properties, no buoyancy, no radiation, and the flow
//! field is taken as already settled. The global energy balance is reported for every run.
use super::{obstacles, CollocatedFlowCase, PolyMeshSpec, RunManifest};
use crate::{
    core::{restarted_gmres, CsrBuilder, Ilu0Preconditioner, LinearConfig, Vec3},
    fvm::PolyMesh,
    io::DataArray,
    require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, path::Path};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ThermalFluid {
    pub specific_heat_j_per_kg_k: f64,
    pub conductivity_w_per_m_k: f64,
}

/// Material and power of one grid obstacle. Conductivity is per grid axis (x, y, z).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ThermalSolid {
    pub obstacle: String,
    pub conductivity_w_per_m_k: [f64; 3],
    pub power_w: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_temperature_k: Option<f64>,
    /// Needed only for physically transient (buoyancy-coupled) runs.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub density_kg_per_m3: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specific_heat_j_per_kg_k: Option<f64>,
    /// Thermal contact resistance on every face this solid shares with another solid, m^2 K/W
    /// (`thermal_devices_v1`). Faces in contact with fluid are not affected, as for a contact
    /// resistance between parts. Where both solids carry one, the larger applies.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub contact_resistance_m2_k_per_w: Option<f64>,
    /// Two-resistor component: the junction temperature reported (and checked against
    /// `max_temperature_k`) is the solid's maximum plus `power_w` times this resistance.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub junction_to_case_k_per_w: Option<f64>,
    /// Heat pipe: replaces the conductivity along one grid axis by the effective value that
    /// gives this end-to-end resistance over the solid's length.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_pipe: Option<HeatPipe>,
}

/// Effective-conductivity heat pipe along grid axis 0 (x), 1 (y) or 2 (z):
/// `k = L / (R A)` with `L` the solid's extent along the axis and `A = volume / L`.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct HeatPipe {
    pub axis: usize,
    pub effective_resistance_k_per_w: f64,
}

/// Thermoelectric cooler in fixed-pumping mode: `heat_pumped_w` is absorbed uniformly in the
/// cold solid and that heat plus `electrical_power_w` is released in the hot solid.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ThermoelectricCooler {
    pub name: String,
    pub cold_solid: String,
    pub hot_solid: String,
    pub heat_pumped_w: f64,
    pub electrical_power_w: f64,
}

/// Boussinesq natural convection: body acceleration -beta (T - T_ref) g on fluid cells.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Buoyancy {
    pub gravity_m_per_s2: [f64; 3],
    pub expansion_per_k: f64,
    pub reference_temperature_k: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ThermalCondition {
    Adiabatic,
    Temperature { temperature_k: f64 },
    HeatFlux { w_per_m2: f64 },
    Convection { h_w_per_m2_k: f64, ambient_k: f64 },
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ThermalPatch {
    pub patch: String,
    pub condition: ThermalCondition,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ConjugateHeat {
    pub fluid: ThermalFluid,
    /// Temperature carried in by every inflowing boundary face.
    pub inflow_temperature_k: f64,
    /// One entry per solid (grid obstacle or surface body); no solid material is implied.
    pub solids: Vec<ThermalSolid>,
    /// Exterior patches without an entry are adiabatic for conduction.
    #[serde(default)]
    pub boundaries: Vec<ThermalPatch>,
    /// Two-way coupling: energy advances with every flow step and drives buoyancy.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub buoyancy: Option<Buoyancy>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub thermoelectric_coolers: Vec<ThermoelectricCooler>,
    /// Heat released on parts of surface bodies (`face_heat_sources_v1`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub face_sources: Vec<FaceHeatSource>,
    /// Power, heat capacity and conductivity by what each cell holds, measured on a finer
    /// lattice (`partial_cells_v1`, `case::partial`).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub partial_cells: bool,
}

/// Power released on part of a surface body's surface: the triangles of one model face. The
/// power goes into the body's cells behind the grid faces of its staircase surface nearest those
/// triangles (exact point-triangle distance), each in proportion to its face area.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct FaceHeatSource {
    pub name: String,
    pub body: String,
    pub triangles: Vec<usize>,
    pub power_w: f64,
}

impl ConjugateHeat {
    /// Whether the case uses any `thermal_devices_v1` feature.
    pub fn uses_devices(&self) -> bool {
        !self.thermoelectric_coolers.is_empty()
            || self.solids.iter().any(|s| {
                s.contact_resistance_m2_k_per_w.is_some()
                    || s.junction_to_case_k_per_w.is_some()
                    || s.heat_pipe.is_some()
            })
    }
}

/// The solved field and its boundary heat accounting.
#[derive(Clone, Debug)]
pub struct EnergySolution {
    pub temperature_k: Vec<f64>,
    pub iterations: usize,
    pub residual_rms: f64,
    /// Conductive heat leaving through each exterior patch, W.
    pub conduction_out_w: Vec<(String, f64)>,
    /// Advected heat leaving minus entering, relative to the inflow temperature, W.
    pub advection_out_w: f64,
    /// Net volume flux out of the boundary, m^3/s (continuity error of the flow field).
    pub net_boundary_flux_m3_per_s: f64,
}

fn axis_of(area: Vec3) -> usize {
    let a = [area.x.abs(), area.y.abs(), area.z.abs()];
    if a[0] >= a[1] && a[0] >= a[2] {
        0
    } else if a[1] >= a[2] {
        1
    } else {
        2
    }
}

fn component(v: [f64; 3], axis: usize) -> f64 {
    v.get(axis).copied().unwrap_or(0.0)
}

fn conditions(mesh: &PolyMesh, patches: &[ThermalPatch]) -> Result<Vec<ThermalCondition>> {
    let mut seen = BTreeSet::new();
    let mut condition = vec![ThermalCondition::Adiabatic; mesh.patches().len()];
    for p in patches {
        require(
            seen.insert(p.patch.as_str()),
            format!("duplicate thermal patch {}", p.patch),
        )?;
        let index = mesh.patch_index(&p.patch)?;
        let ok = match p.condition {
            ThermalCondition::Adiabatic => true,
            ThermalCondition::Temperature { temperature_k } => {
                temperature_k.is_finite() && temperature_k > 0.0
            }
            ThermalCondition::HeatFlux { w_per_m2 } => w_per_m2.is_finite(),
            ThermalCondition::Convection {
                h_w_per_m2_k,
                ambient_k,
            } => {
                h_w_per_m2_k.is_finite()
                    && h_w_per_m2_k > 0.0
                    && ambient_k.is_finite()
                    && ambient_k > 0.0
            }
        };
        require(ok, format!("invalid thermal condition on {}", p.patch))?;
        if let Some(slot) = condition.get_mut(index) {
            *slot = p.condition.clone();
        }
    }
    Ok(condition)
}

fn half(mesh: &PolyMesh, cell: usize, face: Vec3, area: Vec3) -> f64 {
    (face - mesh.cells()[cell].center)
        .dot(area / area.norm())
        .abs()
}

/// Assemble and solve for the rise above the inflow temperature. With `transient`
/// (capacity / dt in W/K per cell, previous rise) one implicit Euler step is taken.
#[allow(clippy::too_many_arguments)]
fn solve_rise(
    mesh: &PolyMesh,
    conductivity: &[[f64; 3]],
    heat_capacity_flux: &[f64],
    source_w: &[f64],
    condition: &[ThermalCondition],
    inflow_temperature_k: f64,
    transient: Option<(&[f64], &[f64])>,
    contact: Option<(&[usize], &[f64])>,
    guess: &mut [f64],
) -> Result<(usize, f64)> {
    let n = mesh.cell_count();
    // The unknown is the rise above the inflow temperature. Besides scaling the system, this
    // keeps a small continuity error of the flow field from acting as a spurious heat source
    // (an absolute-temperature upwind form would add T_in times each cell's flux divergence).
    let mut a = CsrBuilder::new(n, n);
    let mut rhs: Vec<f64> = source_w.to_vec();
    if let Some((capacity, old)) = transient {
        for c in 0..n {
            a.add(c, c, capacity[c])?;
            rhs[c] += capacity[c] * old[c];
        }
    }
    for (fi, f) in mesh.faces().iter().enumerate() {
        let area = f.area.norm();
        let axis = axis_of(f.area);
        let p = f.owner;
        let kp = component(conductivity[p], axis);
        let dp = half(mesh, p, f.center, f.area);
        let m = heat_capacity_flux[fi];
        match f.neighbour {
            Some(q) => {
                let kq = component(conductivity[q], axis);
                let dq = half(mesh, q, f.center, f.area);
                // A contact resistance acts between two different solids only.
                let r = contact.map_or(0.0, |(region, resistance)| {
                    let (rp, rq) = (region[p], region[q]);
                    if rp > 0 && rq > 0 && rp != rq {
                        resistance[p].max(resistance[q])
                    } else {
                        0.0
                    }
                });
                let g = area / (dp / kp + dq / kq + r);
                a.add(p, p, g)?;
                a.add(q, q, g)?;
                a.add(p, q, -g)?;
                a.add(q, p, -g)?;
                if m > 0.0 {
                    a.add(p, p, m)?;
                    a.add(q, p, -m)?;
                } else if m < 0.0 {
                    a.add(q, q, -m)?;
                    a.add(p, q, m)?;
                }
            }
            None => {
                if m > 0.0 {
                    a.add(p, p, m)?;
                }
                // Inflow carries the inflow temperature: zero rise, so nothing on the right side.
                match f
                    .patch
                    .and_then(|i| condition.get(i))
                    .cloned()
                    .unwrap_or(ThermalCondition::Adiabatic)
                {
                    ThermalCondition::Adiabatic => {}
                    ThermalCondition::Temperature { temperature_k } => {
                        let g = kp * area / dp;
                        a.add(p, p, g)?;
                        rhs[p] += g * (temperature_k - inflow_temperature_k);
                    }
                    ThermalCondition::Convection {
                        h_w_per_m2_k,
                        ambient_k,
                    } => {
                        let g = area / (dp / kp + 1.0 / h_w_per_m2_k);
                        a.add(p, p, g)?;
                        rhs[p] += g * (ambient_k - inflow_temperature_k);
                    }
                    ThermalCondition::HeatFlux { w_per_m2 } => rhs[p] += w_per_m2 * area,
                }
            }
        }
    }
    // Row scaling by the diagonal puts residuals in kelvin, so one tolerance fits every
    // problem size and material contrast.
    let unscaled = a.build(0.0)?;
    let diagonal = unscaled.diagonal();
    require(
        diagonal.iter().all(|d| d.is_finite() && *d > 0.0),
        "energy matrix needs a positive diagonal",
    )?;
    let mut scaled = CsrBuilder::new(n, n);
    for row in 0..n {
        for k in unscaled.row_offsets()[row]..unscaled.row_offsets()[row + 1] {
            scaled.add(
                row,
                unscaled.columns()[k],
                unscaled.values()[k] / diagonal[row],
            )?;
        }
        rhs[row] /= diagonal[row];
    }
    let matrix = scaled.build(0.0)?;
    let ilu = Ilu0Preconditioner::new(&matrix, 1e-14)?;
    let result = restarted_gmres(
        &rhs,
        guess,
        |v, y| matrix.multiply(v, y),
        &ilu,
        LinearConfig {
            max_iterations: 5_000,
            restart: 60,
            relative_tolerance: 1e-10,
        },
    )?;
    if !result.converged {
        return Err(SolverError::NumericalFailure(format!(
            "conjugate energy GMRES did not converge ({} iterations, residual {:e})",
            result.iterations, result.residual_rms
        )));
    }
    Ok((result.iterations, result.residual_rms))
}

#[allow(clippy::too_many_arguments)]
fn account(
    mesh: &PolyMesh,
    conductivity: &[[f64; 3]],
    heat_capacity_flux: &[f64],
    condition: &[ThermalCondition],
    inflow_temperature_k: f64,
    rho_cp: f64,
    rise: &[f64],
    iterations: usize,
    residual_rms: f64,
) -> EnergySolution {
    let temperature_k: Vec<f64> = rise.iter().map(|r| r + inflow_temperature_k).collect();
    let mut conduction_out_w: Vec<(String, f64)> = mesh
        .patches()
        .iter()
        .map(|p| (p.name.clone(), 0.0))
        .collect();
    let mut advection_out_w = 0.0;
    let mut net_flux = 0.0;
    for (fi, f) in mesh
        .faces()
        .iter()
        .enumerate()
        .filter(|(_, f)| f.neighbour.is_none())
    {
        let p = f.owner;
        let m = heat_capacity_flux[fi];
        net_flux += m / rho_cp;
        if m > 0.0 {
            advection_out_w += m * rise[p];
        }
        let area = f.area.norm();
        let kp = component(conductivity[p], axis_of(f.area));
        let dp = half(mesh, p, f.center, f.area);
        let Some(patch) = f.patch else { continue };
        let q = match condition.get(patch) {
            Some(ThermalCondition::Temperature {
                temperature_k: wall,
            }) => kp * area / dp * (temperature_k[p] - wall),
            Some(ThermalCondition::Convection {
                h_w_per_m2_k,
                ambient_k,
            }) => area / (dp / kp + 1.0 / h_w_per_m2_k) * (temperature_k[p] - ambient_k),
            Some(ThermalCondition::HeatFlux { w_per_m2 }) => -w_per_m2 * area,
            _ => 0.0,
        };
        if let Some(slot) = conduction_out_w.get_mut(patch) {
            slot.1 += q;
        }
    }
    EnergySolution {
        temperature_k,
        iterations,
        residual_rms,
        conduction_out_w,
        advection_out_w,
        net_boundary_flux_m3_per_s: net_flux,
    }
}

/// Solve steady `div(rho cp u T) - div(K grad T) = q` on the background grid of `spec`.
/// `heat_capacity_flux` is rho*cp*volume-flux per background face (W/K, owner to neighbour),
/// `source_w` the heat generated in each cell. This is the operator used by conjugate runs
/// and by the verification suite.
pub fn solve_energy(
    spec: &PolyMeshSpec,
    conductivity: &[[f64; 3]],
    heat_capacity_flux: &[f64],
    source_w: &[f64],
    patches: &[ThermalPatch],
    inflow_temperature_k: f64,
    rho_cp: f64,
) -> Result<EnergySolution> {
    let (mesh, _) = spec.build()?;
    let n = mesh.cell_count();
    require(
        conductivity.len() == n
            && source_w.len() == n
            && heat_capacity_flux.len() == mesh.face_count(),
        "energy input sizes",
    )?;
    require(
        conductivity
            .iter()
            .flatten()
            .all(|k| k.is_finite() && *k > 0.0),
        "conductivities must be finite and positive",
    )?;
    require(
        source_w
            .iter()
            .chain(heat_capacity_flux)
            .all(|v| v.is_finite()),
        "nonfinite energy source or flux",
    )?;
    require(
        inflow_temperature_k.is_finite() && inflow_temperature_k > 0.0,
        "inflow temperature must be positive kelvin",
    )?;
    let condition = conditions(&mesh, patches)?;
    let mut rise = vec![0.0; n];
    let (iterations, residual) = solve_rise(
        &mesh,
        conductivity,
        heat_capacity_flux,
        source_w,
        &condition,
        inflow_temperature_k,
        None,
        None,
        &mut rise,
    )?;
    Ok(account(
        &mesh,
        conductivity,
        heat_capacity_flux,
        &condition,
        inflow_temperature_k,
        rho_cp,
        &rise,
        iterations,
        residual,
    ))
}

/// Adds each face heat source's power to the cells of its body behind the surface grid faces
/// nearest its triangles: faces between a cell of the body and a cell outside it, or on the
/// domain boundary.
fn face_sources(
    case: &CollocatedFlowCase,
    heat: &ConjugateHeat,
    mesh: &PolyMesh,
    cell_solid: &[Option<usize>],
    source: &mut [f64],
) -> Result<()> {
    let mut names = BTreeSet::new();
    for fs in &heat.face_sources {
        require(
            !fs.name.is_empty() && names.insert(fs.name.as_str()),
            "face heat sources need unique names",
        )?;
        require(fs.power_w.is_finite(), format!("face heat source {}: power must be finite", fs.name))?;
        let bi = case.bodies.iter().position(|b| b.name == fs.body).ok_or_else(|| {
            SolverError::InvalidInput(format!("face heat source {} names no surface body {}", fs.name, fs.body))
        })?;
        let body = &case.bodies[bi];
        require(!fs.triangles.is_empty(), format!("face heat source {} lists no triangles", fs.name))?;
        let mut wanted = vec![false; body.triangles.len()];
        for &t in &fs.triangles {
            let slot = wanted.get_mut(t).ok_or_else(|| {
                SolverError::InvalidInput(format!("face heat source {} lists triangle {t}, which {} does not have", fs.name, fs.body))
            })?;
            *slot = true;
        }
        let si = case.obstacles.len() + bi;
        let inside = |c: usize| cell_solid.get(c).copied().flatten() == Some(si);
        // Surface faces of the body: (its cell, face centre, face area).
        let surface: Vec<(usize, [f64; 3], f64)> = mesh
            .faces()
            .iter()
            .filter_map(|f| {
                let cell = match (inside(f.owner), f.neighbour.map(inside)) {
                    (true, Some(false) | None) => f.owner,
                    (false, Some(true)) => f.neighbour?,
                    _ => return None,
                };
                Some((cell, [f.center.x, f.center.y, f.center.z], f.area.norm()))
            })
            .collect();
        let centres: Vec<[f64; 3]> = surface.iter().map(|(_, c, _)| *c).collect();
        let nearest = body.nearest_triangles(&centres)?;
        let chosen: Vec<(usize, f64)> = surface
            .iter()
            .zip(&nearest)
            .filter(|(_, t)| wanted.get(**t).copied().unwrap_or(false))
            .map(|((cell, _, area), _)| (*cell, *area))
            .collect();
        let total: f64 = chosen.iter().map(|(_, a)| a).sum();
        require(
            total > 0.0,
            format!("face heat source {} covers no grid face of {}; refine the grid near it", fs.name, fs.body),
        )?;
        for (cell, area) in chosen {
            if let Some(slot) = source.get_mut(cell) {
                *slot += fs.power_w * area / total;
            }
        }
    }
    Ok(())
}

/// Per-component verdict against an optional limit. Warn when the margin is below 5% of the
/// allowed rise above the inflow temperature.
fn verdict(max_k: f64, limit: Option<f64>, inflow_k: f64) -> (&'static str, f64) {
    match limit {
        None => ("no_limit", -1.0),
        Some(limit) if max_k > limit => ("fail", 2.0),
        Some(limit) if limit - max_k < 0.05 * (limit - inflow_k).max(0.0) => ("warn", 1.0),
        Some(_) => ("pass", 0.0),
    }
}

/// Energy field of one flow run: steady at the end (one-way), or advanced with every flow
/// step and feeding buoyancy back (two-way).
pub(super) struct EnergyState {
    mesh: PolyMesh,
    region: Vec<usize>,
    names: Vec<String>,
    conductivity: Vec<[f64; 3]>,
    /// Contact resistance of each cell's solid (0 for fluid), m^2 K/W.
    contact: Vec<f64>,
    source: Vec<f64>,
    capacity_j_per_k: Vec<f64>,
    pseudo_capacity: bool,
    /// The cells each solid's temperatures are read from: its own, or for a part with none (a
    /// lumped body) the cells its volume lies in.
    readout: Vec<Vec<usize>>,
    condition: Vec<ThermalCondition>,
    origins: Vec<(usize, bool)>,
    fluid_cells: Vec<usize>,
    rho_cp: f64,
    rise: Vec<f64>,
    storage_w: f64,
    iterations: usize,
    residual: f64,
}

impl EnergyState {
    pub fn new(case: &CollocatedFlowCase, heat: &ConjugateHeat) -> Result<Self> {
        let spec = &case.mesh;
        require(
            spec.xy_shear == 0.0,
            "conjugate heat transfer requires zero xy_shear",
        )?;
        let fluid = &heat.fluid;
        require(
            fluid.specific_heat_j_per_kg_k.is_finite()
                && fluid.specific_heat_j_per_kg_k > 0.0
                && fluid.conductivity_w_per_m_k.is_finite()
                && fluid.conductivity_w_per_m_k > 0.0,
            "fluid thermal properties must be finite and positive",
        )?;
        require(
            heat.inflow_temperature_k.is_finite() && heat.inflow_temperature_k > 0.0,
            "inflow temperature must be positive kelvin",
        )?;
        let solid_names = case.solid_names();
        let mut named = BTreeSet::new();
        for s in &heat.solids {
            require(
                named.insert(s.obstacle.as_str()),
                "each solid needs exactly one material entry",
            )?;
            require(
                solid_names.contains(&s.obstacle.as_str()),
                format!("unknown solid {}", s.obstacle),
            )?;
            require(
                s.conductivity_w_per_m_k
                    .iter()
                    .all(|k| k.is_finite() && *k > 0.0)
                    && s.power_w.is_finite()
                    && s.power_w >= 0.0,
                "solid conductivity must be positive and power non-negative",
            )?;
            require(
                s.max_temperature_k.is_none_or(|t| t.is_finite() && t > 0.0),
                "temperature limits are positive kelvin",
            )?;
            require(
                s.density_kg_per_m3.is_none_or(|v| v.is_finite() && v > 0.0)
                    && s.specific_heat_j_per_kg_k
                        .is_none_or(|v| v.is_finite() && v > 0.0),
                "solid density and specific heat must be positive",
            )?;
            require(
                s.contact_resistance_m2_k_per_w
                    .is_none_or(|r| r.is_finite() && r >= 0.0)
                    && s.junction_to_case_k_per_w
                        .is_none_or(|r| r.is_finite() && r >= 0.0),
                format!("{}: contact and junction resistances must be finite and non-negative", s.obstacle),
            )?;
            if let Some(pipe) = s.heat_pipe {
                require(
                    pipe.axis < 3
                        && pipe.effective_resistance_k_per_w.is_finite()
                        && pipe.effective_resistance_k_per_w > 0.0,
                    format!("{}: a heat pipe needs axis 0, 1 or 2 and a positive resistance", s.obstacle),
                )?;
            }
        }
        let mut coolers = BTreeSet::new();
        for tec in &heat.thermoelectric_coolers {
            require(
                coolers.insert(tec.name.as_str()) && !tec.name.is_empty(),
                "thermoelectric coolers need unique names",
            )?;
            require(
                solid_names.contains(&tec.cold_solid.as_str())
                    && solid_names.contains(&tec.hot_solid.as_str())
                    && tec.cold_solid != tec.hot_solid,
                format!("{}: cold and hot sides must be two different solids", tec.name),
            )?;
            require(
                tec.heat_pumped_w.is_finite()
                    && tec.heat_pumped_w >= 0.0
                    && tec.electrical_power_w.is_finite()
                    && tec.electrical_power_w >= 0.0,
                format!("{}: pumped heat and electrical power must be finite and non-negative", tec.name),
            )?;
        }
        require(
            named.len() == solid_names.len(),
            "every solid needs a material entry",
        )?;
        for p in &heat.boundaries {
            require(
                ["left", "right", "bottom", "top", "front", "back"].contains(&p.patch.as_str()),
                "thermal boundaries apply to exterior patches",
            )?;
        }
        if let Some(b) = &heat.buoyancy {
            require(
                b.gravity_m_per_s2.iter().all(|g| g.is_finite())
                    && b.expansion_per_k.is_finite()
                    && b.expansion_per_k > 0.0
                    && b.reference_temperature_k.is_finite()
                    && b.reference_temperature_k > 0.0,
                "buoyancy needs finite gravity, positive expansion and a kelvin reference",
            )?;
        }
        let (mesh, _) = spec.build()?;
        let n = mesh.cell_count();
        let rho_cp = case.config.density * fluid.specific_heat_j_per_kg_k;
        let cell_solid = if solid_names.is_empty() {
            vec![None; n]
        } else {
            obstacles::layout(spec, &case.obstacles, &case.bodies)?.cell_solid
        };
        let mut region = vec![0usize; n];
        let mut contact = vec![0.0; n];
        let mut conductivity = vec![[fluid.conductivity_w_per_m_k; 3]; n];
        let mut source = vec![0.0; n];
        let mut capacity: Vec<f64> = mesh.cells().iter().map(|c| rho_cp * c.volume).collect();
        let mut pseudo_capacity = false;
        // Each solid's conductivity and volumetric heat capacity, for partial cells.
        let mut solid_k = Vec::with_capacity(solid_names.len());
        let mut solid_capacity = Vec::with_capacity(solid_names.len());
        for (si, name) in solid_names.iter().enumerate() {
            let solid = heat
                .solids
                .iter()
                .find(|s| s.obstacle == *name)
                .ok_or_else(|| {
                    SolverError::InvalidInput(format!("solid {name} has no material entry"))
                })?;
            let cells: Vec<usize> = (0..n).filter(|c| cell_solid[*c] == Some(si)).collect();
            let per_cell = solid.power_w / cells.len().max(1) as f64;
            let solid_rho_cp = match (solid.density_kg_per_m3, solid.specific_heat_j_per_kg_k) {
                (Some(rho), Some(cp)) => rho * cp,
                _ => {
                    pseudo_capacity = true;
                    rho_cp
                }
            };
            let mut k = solid.conductivity_w_per_m_k;
            if let Some(pipe) = solid.heat_pipe {
                // Extent along the axis from the cell centres plus one cell, area from volume.
                // Extent along the axis from the faces of the outermost cells.
                let along = |c: &usize| match pipe.axis {
                    0 => mesh.cells()[*c].center.x,
                    1 => mesh.cells()[*c].center.y,
                    _ => mesh.cells()[*c].center.z,
                };
                let nodes = spec.nodes(pipe.axis);
                let face = |x: f64, upper: bool| {
                    spec.cell_of(pipe.axis, x)
                        .and_then(|i| nodes.get(if upper { i + 1 } else { i }).copied())
                        .unwrap_or(x)
                };
                let low = face(cells.iter().map(along).fold(f64::MAX, f64::min), false);
                let high = face(cells.iter().map(along).fold(f64::MIN, f64::max), true);
                let length = high - low;
                let volume: f64 = cells.iter().map(|c| mesh.cells()[*c].volume).sum();
                require(
                    length > 0.0 && volume > 0.0,
                    format!("heat pipe {name} covers no cells"),
                )?;
                let area = volume / length;
                k[pipe.axis] = length / (pipe.effective_resistance_k_per_w * area);
            }
            solid_k.push(k);
            solid_capacity.push(solid_rho_cp);
            for c in cells {
                region[c] = si + 1;
                conductivity[c] = k;
                contact[c] = solid.contact_resistance_m2_k_per_w.unwrap_or(0.0);
                if !heat.partial_cells {
                    source[c] = per_cell;
                }
                capacity[c] = solid_rho_cp * mesh.cells()[c].volume;
            }
        }
        let mut readout: Vec<Vec<usize>> = (0..solid_names.len()).map(|si| (0..n).filter(|c| region[*c] == si + 1).collect()).collect();
        if !heat.partial_cells {
            if let Some(si) = readout.iter().position(Vec::is_empty) {
                return Err(SolverError::InvalidInput(format!(
                    "solid {} has no cell of its own; a lumped body needs thermal.partial_cells",
                    solid_names[si]
                )));
            }
        }
        if heat.partial_cells {
            let boxes = case.obstacles.len();
            let partial = super::partial::measure(spec, &cell_solid, boxes, &case.bodies, [fluid.conductivity_w_per_m_k; 3], &solid_k)?;
            // Solid volume of each part, by cell.
            let mut volume_of = vec![0.0; solid_names.len()];
            for (c, shares) in partial.shares.iter().enumerate() {
                for (s, share) in shares {
                    if let Some(v) = volume_of.get_mut(*s) {
                        *v += share * mesh.cells()[c].volume;
                    }
                }
            }
            for (c, shares) in partial.shares.iter().enumerate() {
                let volume = mesh.cells()[c].volume;
                if let Some(k) = partial.conductivity[c] {
                    conductivity[c] = k;
                }
                let solid_share: f64 = shares.iter().map(|(_, share)| share).sum();
                capacity[c] = rho_cp * volume * (1.0 - solid_share).max(0.0)
                    + shares.iter().map(|(s, share)| solid_capacity.get(*s).copied().unwrap_or(rho_cp) * share * volume).sum::<f64>();
                for (s, share) in shares {
                    let power = heat.solids.iter().find(|m| Some(m.obstacle.as_str()) == solid_names.get(*s).copied()).map_or(0.0, |m| m.power_w);
                    if volume_of[*s] > 0.0 {
                        source[c] += power * share * volume / volume_of[*s];
                    }
                }
            }
            for (si, cells) in readout.iter_mut().enumerate() {
                if !cells.is_empty() {
                    continue;
                }
                // A part with no cell of its own: the cells its volume lies in, or the cell at
                // its box's centre when it is thinner than a sub-cell (its power goes there too).
                *cells = (0..n).filter(|c| partial.shares[*c].iter().any(|(s, _)| *s == si)).collect();
                if cells.is_empty() {
                    let body = si.checked_sub(boxes).and_then(|b| case.bodies.get(b));
                    let centre = body.and_then(|b| {
                        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
                        for v in &b.vertices {
                            for a in 0..3 {
                                lo[a] = lo[a].min(v[a]);
                                hi[a] = hi[a].max(v[a]);
                            }
                        }
                        let at: Option<Vec<usize>> = (0..3).map(|a| spec.cell_of(a, 0.5 * (lo[a] + hi[a]))).collect();
                        at.map(|at| (at[2] * spec.dimensions[1] + at[1]) * spec.dimensions[0] + at[0])
                    });
                    let Some(c) = centre.filter(|c| *c < n) else {
                        return Err(SolverError::InvalidInput(format!("solid {} lies outside the grid", solid_names[si])));
                    };
                    let power = heat.solids.iter().find(|m| m.obstacle == solid_names[si]).map_or(0.0, |m| m.power_w);
                    source[c] += power;
                    *cells = vec![c];
                }
            }
        }
        for tec in &heat.thermoelectric_coolers {
            for (side, watts) in [
                (&tec.cold_solid, -tec.heat_pumped_w),
                (&tec.hot_solid, tec.heat_pumped_w + tec.electrical_power_w),
            ] {
                let index = solid_names.iter().position(|s| s == side).map(|i| i + 1);
                let cells: Vec<usize> = (0..n).filter(|c| Some(region[*c]) == index).collect();
                require(!cells.is_empty(), format!("{}: {side} covers no cells", tec.name))?;
                let per_cell = watts / cells.len() as f64;
                for c in cells {
                    source[c] += per_cell;
                }
            }
        }
        face_sources(case, heat, &mesh, &cell_solid, &mut source)?;
        let condition = conditions(&mesh, &heat.boundaries)?;
        let origins = obstacles::fluid_face_origins(spec, &case.obstacles, &case.bodies)?;
        let fluid_cells = (0..n).filter(|c| cell_solid[*c].is_none()).collect();
        Ok(Self {
            mesh,
            region,
            names: solid_names.iter().map(|s| s.to_string()).collect(),
            conductivity,
            contact,
            source,
            capacity_j_per_k: capacity,
            pseudo_capacity,
            readout,
            condition,
            origins,
            fluid_cells,
            rho_cp,
            rise: vec![0.0; n],
            storage_w: 0.0,
            iterations: 0,
            residual: 0.0,
        })
    }

    fn capacity_flux(&self, fluid_flux: &[f64]) -> Result<Vec<f64>> {
        require(
            self.origins.len() == fluid_flux.len(),
            "fluid/background face mapping mismatch",
        )?;
        let mut out = vec![0.0; self.mesh.face_count()];
        for ((face, flipped), phi) in self.origins.iter().zip(fluid_flux) {
            if let Some(slot) = out.get_mut(*face) {
                *slot = self.rho_cp * if *flipped { -phi } else { *phi };
            }
        }
        Ok(out)
    }

    pub fn steady(&mut self, fluid_flux: &[f64], inflow_k: f64) -> Result<()> {
        let flux = self.capacity_flux(fluid_flux)?;
        let mut rise = self.rise.clone();
        let (it, res) = solve_rise(
            &self.mesh,
            &self.conductivity,
            &flux,
            &self.source,
            &self.condition,
            inflow_k,
            None,
            Some((&self.region, &self.contact)),
            &mut rise,
        )?;
        self.rise = rise;
        self.storage_w = 0.0;
        self.iterations = it;
        self.residual = res;
        Ok(())
    }

    /// One implicit Euler step of length `dt` with the latest flow fluxes.
    pub fn step(&mut self, fluid_flux: &[f64], inflow_k: f64, dt: f64) -> Result<()> {
        require(
            dt.is_finite() && dt > 0.0,
            "energy time step must be positive",
        )?;
        let flux = self.capacity_flux(fluid_flux)?;
        let capacity: Vec<f64> = self.capacity_j_per_k.iter().map(|c| c / dt).collect();
        let old = self.rise.clone();
        let mut rise = old.clone();
        let (it, res) = solve_rise(
            &self.mesh,
            &self.conductivity,
            &flux,
            &self.source,
            &self.condition,
            inflow_k,
            Some((&capacity, &old)),
            Some((&self.region, &self.contact)),
            &mut rise,
        )?;
        self.storage_w = capacity
            .iter()
            .zip(rise.iter().zip(&old))
            .map(|(c, (new, old))| c * (new - old))
            .sum();
        self.rise = rise;
        self.iterations = it;
        self.residual = res;
        Ok(())
    }

    /// Temperatures of the fluid cells in flow-mesh order.
    pub fn fluid_temperature(&self, inflow_k: f64) -> Vec<f64> {
        self.fluid_cells
            .iter()
            .map(|c| self.rise.get(*c).copied().unwrap_or(0.0) + inflow_k)
            .collect()
    }

    /// Temperatures of every background-grid cell, fluid and solid, K.
    pub fn cell_temperature(&self, inflow_k: f64) -> Vec<f64> {
        self.rise.iter().map(|r| r + inflow_k).collect()
    }

    /// Body acceleration per fluid cell, `base - beta (T - T_ref) g`.
    pub fn buoyancy(&self, b: &Buoyancy, base: Vec3, inflow_k: f64) -> Vec<Vec3> {
        let g = Vec3::new(
            b.gravity_m_per_s2[0],
            b.gravity_m_per_s2[1],
            b.gravity_m_per_s2[2],
        );
        self.fluid_temperature(inflow_k)
            .iter()
            .map(|t| base - g * (b.expansion_per_k * (t - b.reference_temperature_k)))
            .collect()
    }

    pub fn finish(
        &self,
        case: &CollocatedFlowCase,
        heat: &ConjugateHeat,
        fluid_flux: &[f64],
        out: &Path,
        m: &mut RunManifest,
    ) -> Result<()> {
        let spec = &case.mesh;
        let flux = self.capacity_flux(fluid_flux)?;
        let solution = account(
            &self.mesh,
            &self.conductivity,
            &flux,
            &self.condition,
            heat.inflow_temperature_k,
            self.rho_cp,
            &self.rise,
            self.iterations,
            self.residual,
        );
        let t = &solution.temperature_k;
        let n = self.mesh.cell_count();
        let region = &self.region;
        let power: f64 = self.source.iter().sum();
        let conducted: f64 = solution.conduction_out_w.iter().map(|(_, q)| q).sum();
        let imbalance = power - solution.advection_out_w - conducted - self.storage_w;
        // Normalise by the gross heat flows, so through-flow cases (heat in one wall, out of
        // another) are judged against the heat actually transported.
        let gross = power.abs()
            + solution
                .conduction_out_w
                .iter()
                .map(|(_, q)| q.abs())
                .sum::<f64>()
            + solution.advection_out_w.abs()
            + self.storage_w.abs();
        let scale = gross;
        let relative = if scale > 0.0 {
            imbalance.abs() / scale
        } else {
            imbalance.abs()
        };
        let mut inflow_volume = 0.0;
        for (fi, f) in self.mesh.faces().iter().enumerate() {
            if f.neighbour.is_none() && flux[fi] < 0.0 {
                inflow_volume -= flux[fi] / self.rho_cp;
            }
        }
        let mut components = vec![];
        let (mut fails, mut warns) = (0.0, 0.0);
        let mut csv: [Vec<f64>; 6] = Default::default();
        for (si, name) in self.names.iter().enumerate() {
            let Some(solid) = heat.solids.iter().find(|s| s.obstacle == *name) else {
                continue;
            };
            let values: Vec<(f64, f64)> = self
                .readout
                .get(si)
                .into_iter()
                .flatten()
                .map(|c| (t[*c], self.mesh.cells()[*c].volume))
                .collect();
            let volume: f64 = values.iter().map(|(_, v)| v).sum();
            let max = values.iter().map(|(t, _)| *t).fold(f64::MIN, f64::max);
            let min = values.iter().map(|(t, _)| *t).fold(f64::MAX, f64::min);
            let mean = values.iter().map(|(t, v)| t * v).sum::<f64>() / volume;
            // A two-resistor component is judged on its junction, not on its case.
            let junction = solid.junction_to_case_k_per_w.map(|r| max + solid.power_w * r);
            let (status, code) = verdict(
                junction.unwrap_or(max),
                solid.max_temperature_k,
                heat.inflow_temperature_k,
            );
            if code == 2.0 {
                fails += 1.0;
            } else if code == 1.0 {
                warns += 1.0;
            }
            for (column, value) in
                csv.iter_mut()
                    .zip([si as f64 + 1.0, solid.power_w, min, mean, max, code])
            {
                column.push(value);
            }
            components.push(serde_json::json!({
                "name": name, "region": si + 1, "power_w": solid.power_w,
                "conductivity_w_per_m_k": solid.conductivity_w_per_m_k,
                "min_temperature_k": min, "mean_temperature_k": mean, "max_temperature_k": max,
                "max_temperature_c": max - 273.15, "limit_k": solid.max_temperature_k,
                "margin_k": solid.max_temperature_k.map(|l| l - junction.unwrap_or(max)), "status": status,
                "junction_temperature_k": junction,
                "contact_resistance_m2_k_per_w": solid.contact_resistance_m2_k_per_w,
                "heat_pipe": solid.heat_pipe.map(|p| serde_json::json!({
                    "axis": p.axis, "effective_resistance_k_per_w": p.effective_resistance_k_per_w,
                    "axial_conductivity_w_per_m_k": (0..n).find(|c| region[*c] == si + 1).map(|c| self.conductivity[c][p.axis])
                }))
            }));
        }
        let fluid_max = (0..n)
            .filter(|c| region[*c] == 0)
            .map(|c| t[c])
            .fold(f64::MIN, f64::max);
        let coupled = heat.buoyancy.is_some();
        let report = serde_json::json!({
            "schema": "cfd-solvers-rs/conjugate-heat/1",
            "model": if coupled { "transient_two_way_boussinesq_cht_upwind" } else { "steady_one_way_forced_convection_cht_upwind" },
            "assumptions": if coupled {
                serde_json::json!(["constant_properties_except_boussinesq_buoyancy", "no_radiation", "laminar", "grid_aligned_orthotropic_solids"])
            } else {
                serde_json::json!(["constant_properties", "no_buoyancy", "no_radiation", "flow_field_taken_as_settled", "grid_aligned_orthotropic_solids"])
            },
            "solid_heat_capacity": if self.pseudo_capacity && coupled { "pseudo_fluid_capacity_for_some_solids_transient_history_not_physical" } else { "as_specified" },
            "buoyancy": heat.buoyancy,
            "inflow_temperature_k": heat.inflow_temperature_k,
            "fluid": heat.fluid, "rho_cp_j_per_m3_k": self.rho_cp,
            "total_power_w": power,
            "face_sources": heat.face_sources.iter().map(|f| serde_json::json!({"name": f.name, "body": f.body, "power_w": f.power_w})).collect::<Vec<_>>(),
            "advected_heat_out_w": solution.advection_out_w,
            "conducted_heat_out_w": solution.conduction_out_w.iter().map(|(p, q)| serde_json::json!({"patch":p,"w":q})).collect::<Vec<_>>(),
            "storage_rate_w": self.storage_w,
            "energy_imbalance_w": imbalance, "relative_energy_imbalance": relative,
            "inflow_volume_m3_per_s": inflow_volume,
            "bulk_outflow_rise_k": if inflow_volume > 0.0 { solution.advection_out_w / (self.rho_cp * inflow_volume) } else { 0.0 },
            "flow_net_boundary_flux_m3_per_s": solution.net_boundary_flux_m3_per_s,
            "max_fluid_temperature_k": fluid_max,
            "components": components,
            "thermoelectric_coolers": heat.thermoelectric_coolers.iter().map(|tec| {
                let side = |name: &str| self.solid_temperature(name, heat.inflow_temperature_k).map(|(_, mean)| mean);
                let cold = side(&tec.cold_solid);
                let hot = side(&tec.hot_solid);
                serde_json::json!({
                    "name": tec.name, "cold_solid": tec.cold_solid, "hot_solid": tec.hot_solid,
                    "heat_pumped_w": tec.heat_pumped_w, "electrical_power_w": tec.electrical_power_w,
                    "heat_released_w": tec.heat_pumped_w + tec.electrical_power_w,
                    "cop": if tec.electrical_power_w > 0.0 { Some(tec.heat_pumped_w / tec.electrical_power_w) } else { None },
                    "cold_mean_temperature_k": cold, "hot_mean_temperature_k": hot,
                    "model": "fixed_heat_pumping"
                })
            }).collect::<Vec<_>>(),
            "warn_rule": "margin below 5% of the allowed rise above the inflow temperature",
            "linear": {"iterations": solution.iterations, "residual_rms": solution.residual_rms},
            "physical_validation": "not_established"
        });
        std::fs::write(
            out.join("thermal_report.json"),
            serde_json::to_vec_pretty(&report)?,
        )?;
        m.files.push("thermal_report.json".into());
        crate::io::write_csv(
            out.join("thermal_components.csv"),
            &[
                ("region", &csv[0]),
                ("power_w", &csv[1]),
                ("min_temperature_k", &csv[2]),
                ("mean_temperature_k", &csv[3]),
                ("max_temperature_k", &csv[4]),
                ("status_code", &csv[5]),
            ],
        )?;
        m.files.push("thermal_components.csv".into());
        let (_, export) = spec.build()?;
        super::write_field(
            out,
            m,
            "thermal.vtu",
            &export,
            vec![],
            vec![
                DataArray::scalar("temperature_k", t.clone()),
                DataArray::scalar("temperature_c", t.iter().map(|v| v - 273.15).collect()),
                DataArray::scalar("region_index", region.iter().map(|r| *r as f64).collect()),
            ],
        )?;
        let solid_max = (0..n)
            .filter(|c| region[*c] > 0)
            .map(|c| t[c])
            .fold(f64::MIN, f64::max);
        for (k, v) in [
            ("thermal_total_power_w", power),
            ("thermal_relative_energy_imbalance", relative),
            (
                "thermal_max_solid_temperature_k",
                if solid_max.is_finite() && solid_max > f64::MIN {
                    solid_max
                } else {
                    fluid_max
                },
            ),
            ("thermal_max_fluid_temperature_k", fluid_max),
            ("thermal_components_failing", fails),
            ("thermal_components_warning", warns),
            ("thermal_linear_iterations", solution.iterations as f64),
        ] {
            m.metrics.insert(k.into(), v);
        }
        Ok(())
    }

    /// Conductive heat leaving one exterior patch for the current field, W.
    pub fn patch_conduction(&self, patch: &str, fluid_flux: &[f64], inflow_k: f64) -> Result<f64> {
        let flux = self.capacity_flux(fluid_flux)?;
        let s = account(
            &self.mesh,
            &self.conductivity,
            &flux,
            &self.condition,
            inflow_k,
            self.rho_cp,
            &self.rise,
            0,
            0.0,
        );
        s.conduction_out_w
            .iter()
            .find(|(p, _)| p == patch)
            .map(|(_, q)| *q)
            .ok_or_else(|| SolverError::InvalidInput(format!("unknown patch {patch}")))
    }

    /// Max and volume-mean temperature of the fluid cells (`solid == false`) or of all solid
    /// cells, K; `None` when there are no such cells.
    pub fn region_temperature(&self, solid: bool, inflow_k: f64) -> Option<(f64, f64)> {
        let cells: Vec<usize> = (0..self.region.len())
            .filter(|c| (self.region[*c] > 0) == solid)
            .collect();
        let volume: f64 = cells.iter().map(|c| self.mesh.cells()[*c].volume).sum();
        if cells.is_empty() || volume <= 0.0 {
            return None;
        }
        let max = cells.iter().map(|c| self.rise[*c]).fold(f64::MIN, f64::max) + inflow_k;
        let mean = cells
            .iter()
            .map(|c| self.rise[*c] * self.mesh.cells()[*c].volume)
            .sum::<f64>()
            / volume
            + inflow_k;
        Some((max, mean))
    }

    /// Temperature of the cell whose centre is nearest to `point`, K.
    pub fn point_temperature(&self, point: [f64; 3], inflow_k: f64) -> Option<f64> {
        let at = Vec3::new(point[0], point[1], point[2]);
        let cell = (0..self.mesh.cell_count()).min_by(|a, b| {
            let da = (self.mesh.cells()[*a].center - at).norm();
            let db = (self.mesh.cells()[*b].center - at).norm();
            da.total_cmp(&db)
        })?;
        Some(self.rise.get(cell).copied()? + inflow_k)
    }

    /// Max and volume-mean temperature of a named solid, K.
    pub fn solid_temperature(&self, name: &str, inflow_k: f64) -> Option<(f64, f64)> {
        // The cells its temperature is read from: its own, or for a lumped part the cells its
        // volume lies in, as the component table reads it.
        let index = self.names.iter().position(|n| n == name)?;
        let cells = self.readout.get(index).filter(|c| !c.is_empty())?;
        let volume: f64 = cells.iter().map(|c| self.mesh.cells()[*c].volume).sum();
        let max = cells.iter().map(|c| self.rise[*c]).fold(f64::MIN, f64::max) + inflow_k;
        let mean = cells
            .iter()
            .map(|c| self.rise[*c] * self.mesh.cells()[*c].volume)
            .sum::<f64>()
            / volume
            + inflow_k;
        Some((max, mean))
    }
}
