//! Small interface laws: conjugate interface heat flux, electrochemical heat source,
//! axial deformation of a sequential optical system, and a lumped porous
//! electro-thermal-flow control volume.
use crate::{optics::SequentialOpticalSystem, require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InterfaceHeatTransfer {
    pub heat_flux: f64,
    pub interface_temperature: f64,
}

/// Series-resistance heat flux from side A to B and the interface temperature.
pub fn conjugate_interface_heat_flux(
    ta: f64,
    ka: f64,
    da: f64,
    tb: f64,
    kb: f64,
    db: f64,
) -> Result<InterfaceHeatTransfer> {
    require(
        [ka, kb, da, db].iter().all(|v| v.is_finite() && *v > 0.0),
        "invalid CHT interface",
    )?;
    let (ra, rb) = (da / ka, db / kb);
    let q = (ta - tb) / (ra + rb);
    Ok(InterfaceHeatTransfer {
        heat_flux: q,
        interface_temperature: ta - q * ra,
    })
}

/// `q = j * eta (+ ohmic)` per entry.
pub fn electrochemical_heat_source(
    current_density: &[f64],
    overpotential: &[f64],
    ohmic: &[f64],
) -> Result<Vec<f64>> {
    require(
        current_density.len() == overpotential.len()
            && (ohmic.is_empty() || ohmic.len() == current_density.len()),
        "electrochemical heat size mismatch",
    )?;
    Ok(current_density
        .iter()
        .zip(overpotential)
        .enumerate()
        .map(|(i, (j, eta))| j * eta + ohmic.get(i).copied().unwrap_or(0.0))
        .collect())
}

/// Moves each surface vertex by the axial displacement evaluated at its vertex position.
pub fn deform_optical_system_axially(
    system: &SequentialOpticalSystem,
    displacement: impl Fn(f64) -> f64,
) -> Result<SequentialOpticalSystem> {
    let mut out = SequentialOpticalSystem::new(system.object_space_index())?;
    for s in system.surfaces() {
        let mut moved = *s;
        let dz = displacement(s.vertex_z);
        require(dz.is_finite(), "non-finite optical deformation")?;
        moved.vertex_z += dz;
        out.add_surface(moved)?;
    }
    Ok(out)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousElectroThermalFlowConfig {
    pub porosity: f64,
    pub permeability_m2: f64,
    pub viscosity_pa_s: f64,
    pub reaction_rate_per_s: f64,
    pub reaction_heat_j_per_mol: f64,
    pub volumetric_heat_capacity_j_per_m3_k: f64,
}

impl Default for PorousElectroThermalFlowConfig {
    fn default() -> Self {
        Self {
            porosity: 0.4,
            permeability_m2: 1e-12,
            viscosity_pa_s: 1e-3,
            reaction_rate_per_s: 1.0,
            reaction_heat_j_per_mol: 1e4,
            volumetric_heat_capacity_j_per_m3_k: 1e6,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousElectroThermalFlowState {
    pub species_concentration_mol_m3: f64,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub darcy_velocity_m_per_s: f64,
}

impl Default for PorousElectroThermalFlowState {
    fn default() -> Self {
        Self {
            species_concentration_mol_m3: 0.0,
            temperature_k: 300.0,
            pressure_pa: 0.0,
            darcy_velocity_m_per_s: 0.0,
        }
    }
}

/// One explicit step of a single-control-volume reaction/heat/Darcy-flow baseline.
pub fn advance_porous_electrothermal_flow(
    mut state: PorousElectroThermalFlowState,
    c: &PorousElectroThermalFlowConfig,
    dt: f64,
    pressure_gradient: f64,
) -> Result<PorousElectroThermalFlowState> {
    require(
        c.porosity > 0.0
            && c.porosity <= 1.0
            && c.permeability_m2 >= 0.0
            && c.viscosity_pa_s > 0.0
            && c.reaction_rate_per_s >= 0.0
            && c.reaction_heat_j_per_mol >= 0.0
            && c.volumetric_heat_capacity_j_per_m3_k > 0.0
            && dt >= 0.0
            && state.species_concentration_mol_m3 >= 0.0
            && state.temperature_k > 0.0
            && pressure_gradient.is_finite(),
        "invalid porous electrothermal flow state",
    )?;
    let rate = c.reaction_rate_per_s * state.species_concentration_mol_m3;
    state.species_concentration_mol_m3 = (state.species_concentration_mol_m3 - dt * rate).max(0.0);
    state.temperature_k +=
        dt * rate * c.reaction_heat_j_per_mol / c.volumetric_heat_capacity_j_per_m3_k;
    state.darcy_velocity_m_per_s = -c.permeability_m2 / c.viscosity_pa_s * pressure_gradient;
    state.pressure_pa += -pressure_gradient * state.darcy_velocity_m_per_s * dt;
    Ok(state)
}
