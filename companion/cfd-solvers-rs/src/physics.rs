//! Shared material record and bridges into the FEM and optics solvers, ported from upstream
//! `cfd/physics/{material, material_bridge}` (MIT, cfd_solvers contributors).
use crate::{require, Result};

#[derive(Clone, Debug, PartialEq)]
pub struct Material {
    pub name: String,
    pub density: f64,
    pub heat_capacity: f64,
    pub thermal_conductivity: f64,
    pub electrical_conductivity: f64,
    pub relative_permittivity: [f64; 3],
    pub relative_permeability: [f64; 3],
    pub refractive_index: f64,
}

impl Default for Material {
    fn default() -> Self {
        Self {
            name: String::new(),
            density: 0.0,
            heat_capacity: 0.0,
            thermal_conductivity: 0.0,
            electrical_conductivity: 0.0,
            relative_permittivity: [1.0; 3],
            relative_permeability: [1.0; 3],
            refractive_index: 1.0,
        }
    }
}

/// Non-empty name, non-negative bulk properties, positive refractive index and positive finite
/// diagonal permittivity/permeability.
pub fn validate_material(m: &Material) -> Result<()> {
    require(
        !m.name.is_empty()
            && m.density >= 0.0
            && m.heat_capacity >= 0.0
            && m.thermal_conductivity >= 0.0
            && m.electrical_conductivity >= 0.0
            && m.refractive_index > 0.0,
        "invalid shared material",
    )?;
    require(
        m.relative_permittivity
            .iter()
            .all(|v| *v > 0.0 && v.is_finite()),
        "invalid permittivity",
    )?;
    require(
        m.relative_permeability
            .iter()
            .all(|v| *v > 0.0 && v.is_finite()),
        "invalid permeability",
    )
}

pub fn fem_conductivity(m: &Material) -> Result<f64> {
    validate_material(m)?;
    Ok(m.electrical_conductivity)
}

pub fn optical_refractive_index(m: &Material) -> Result<f64> {
    validate_material(m)?;
    Ok(m.refractive_index)
}

/// Copies the diagonal relative permittivity into the 3-D FDTD anisotropy fields.
pub fn apply_to_fdtd(
    mut config: crate::fdtd::Maxwell3DConfig,
    material: &Material,
) -> Result<crate::fdtd::Maxwell3DConfig> {
    validate_material(material)?;
    [config.epsilon_rx, config.epsilon_ry, config.epsilon_rz] = material.relative_permittivity;
    Ok(config)
}
