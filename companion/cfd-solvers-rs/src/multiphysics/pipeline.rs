//! Staged one-way multiphysics pipelines with finiteness checks between stages.
use crate::{require, Result};

pub type Field = Vec<f64>;
pub type Stage<'a> = Box<dyn Fn(&[f64]) -> Result<Field> + 'a>;

fn checked(field: Field, what: &str) -> Result<Field> {
    require(
        !field.is_empty(),
        format!("{what} stage returned an empty field"),
    )?;
    require(
        field.iter().all(|v| v.is_finite()),
        format!("{what} stage produced non-finite data"),
    )?;
    Ok(field)
}

pub struct CfdThermalStructuralOpticsPipeline<'a> {
    pub cfd_to_thermal: Stage<'a>,
    pub thermal_to_structural: Stage<'a>,
    pub structural_to_optics: Stage<'a>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PipelineResult {
    pub thermal: Field,
    pub structural: Field,
    pub optics: Field,
}

pub fn run_cfd_thermal_structural_optics(
    pipeline: &CfdThermalStructuralOpticsPipeline<'_>,
    cfd: &[f64],
) -> Result<PipelineResult> {
    require(!cfd.is_empty(), "incomplete multiphysics pipeline")?;
    let thermal = checked((pipeline.cfd_to_thermal)(cfd)?, "thermal")?;
    let structural = checked((pipeline.thermal_to_structural)(&thermal)?, "structural")?;
    let optics = checked((pipeline.structural_to_optics)(&structural)?, "optics")?;
    Ok(PipelineResult {
        thermal,
        structural,
        optics,
    })
}

pub struct CorrosionCfdStructuralPipeline<'a> {
    pub recession_to_mesh: Stage<'a>,
    pub mesh_to_cfd: Stage<'a>,
    pub cfd_to_structural: Stage<'a>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CorrosionCfdStructuralResult {
    pub mesh: Field,
    pub cfd: Field,
    pub structural: Field,
}

pub fn run_corrosion_recession_mesh_cfd_structural(
    pipeline: &CorrosionCfdStructuralPipeline<'_>,
    recession: &[f64],
) -> Result<CorrosionCfdStructuralResult> {
    require(
        !recession.is_empty(),
        "incomplete corrosion/CFD/structural pipeline",
    )?;
    let mesh = checked((pipeline.recession_to_mesh)(recession)?, "mesh")?;
    let cfd = checked((pipeline.mesh_to_cfd)(&mesh)?, "cfd")?;
    let structural = checked((pipeline.cfd_to_structural)(&cfd)?, "structural")?;
    Ok(CorrosionCfdStructuralResult {
        mesh,
        cfd,
        structural,
    })
}
