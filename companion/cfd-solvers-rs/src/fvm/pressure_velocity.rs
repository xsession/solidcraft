use super::{
    gauss_gradient_scalar,
    operators::{linear_scalar, linear_vector, scalar_sizes, vector_sizes},
    PolyMesh,
};
use crate::{core::Vec3, finite_slice, require, Result};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VelocityBoundaryType {
    FixedValue,
    ZeroGradient,
    Slip,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct VelocityBoundaryCondition {
    pub boundary_type: VelocityBoundaryType,
    pub value: Vec3,
}
impl Default for VelocityBoundaryCondition {
    fn default() -> Self {
        Self {
            boundary_type: VelocityBoundaryType::ZeroGradient,
            value: Vec3::ZERO,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PressureBoundaryType {
    FixedValue,
    ZeroGradient,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct PressureBoundaryCondition {
    pub boundary_type: PressureBoundaryType,
    pub value: f64,
}
impl Default for PressureBoundaryCondition {
    fn default() -> Self {
        Self {
            boundary_type: PressureBoundaryType::ZeroGradient,
            value: 0.0,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct FaceOrthogonalDecomposition {
    pub delta: Vec3,
    pub orthogonal_area: Vec3,
    pub nonorthogonal_area: Vec3,
    pub orthogonal_metric: f64,
}
pub fn decompose_face_area(m: &PolyMesh, fi: usize) -> Result<FaceOrthogonalDecomposition> {
    require(fi < m.face_count(), "face index out of range")?;
    let f = m.faces()[fi];
    let d = if let Some(n) = f.neighbour {
        m.cells()[n].center - m.cells()[f.owner].center
    } else {
        f.center - m.cells()[f.owner].center
    };
    let d2 = d.dot(d);
    require(
        d2.is_finite() && d2 > 0.0,
        "degenerate owner-neighbour displacement",
    )?;
    let metric = f.area.dot(d) / d2;
    require(
        metric.is_finite() && metric > 0.0,
        "face area has invalid orientation",
    )?;
    let orthogonal_area = d * metric;
    Ok(FaceOrthogonalDecomposition {
        delta: d,
        orthogonal_area,
        nonorthogonal_area: f.area - orthogonal_area,
        orthogonal_metric: metric,
    })
}
pub fn boundary_velocity_values(
    m: &PolyMesh,
    u: &[Vec3],
    bc: &[VelocityBoundaryCondition],
) -> Result<Vec<Vec3>> {
    vector_sizes(m, u, &[])?;
    require(
        bc.len() == m.patches().len() && bc.iter().all(|b| b.value.is_finite()),
        "velocity patch conditions",
    )?;
    let mut b = vec![Vec3::ZERO; m.face_count()];
    for (fi, f) in m.faces().iter().enumerate() {
        let Some(p) = f.patch else {
            continue;
        };
        let condition = bc[p];
        b[fi] = match condition.boundary_type {
            VelocityBoundaryType::FixedValue => condition.value,
            VelocityBoundaryType::ZeroGradient => u[f.owner],
            VelocityBoundaryType::Slip => {
                let n = f.area / f.area.norm();
                let relative = u[f.owner] - condition.value;
                condition.value + relative - n * relative.dot(n)
            }
        };
    }
    Ok(b)
}
pub fn boundary_pressure_values(
    m: &PolyMesh,
    p: &[f64],
    bc: &[PressureBoundaryCondition],
) -> Result<Vec<f64>> {
    scalar_sizes(m, p, &[])?;
    require(
        bc.len() == m.patches().len() && bc.iter().all(|b| b.value.is_finite()),
        "pressure patch conditions",
    )?;
    let mut b = vec![0.0; m.face_count()];
    for (fi, f) in m.faces().iter().enumerate() {
        let Some(patch) = f.patch else {
            continue;
        };
        b[fi] = if bc[patch].boundary_type == PressureBoundaryType::FixedValue {
            bc[patch].value
        } else {
            p[f.owner]
        };
    }
    Ok(b)
}
pub fn rhie_chow_face_flux(
    m: &PolyMesh,
    h: &[Vec3],
    mobility: &[f64],
    p: &[f64],
    velocity_bc: &[VelocityBoundaryCondition],
    pressure_bc: &[PressureBoundaryCondition],
    include_nonorthogonal: bool,
) -> Result<Vec<f64>> {
    vector_sizes(m, h, &[])?;
    scalar_sizes(m, mobility, &[])?;
    scalar_sizes(m, p, &[])?;
    let pb = boundary_pressure_values(m, p, pressure_bc)?;
    let hb = boundary_velocity_values(m, h, velocity_bc)?;
    let gp = gauss_gradient_scalar(m, p, &pb)?;
    let mut r = vec![0.0; m.face_count()];
    for (fi, out) in r.iter_mut().enumerate() {
        let f = m.faces()[fi];
        let geom = decompose_face_area(m, fi)?;
        if f.patch
            .is_some_and(|b| velocity_bc[b].boundary_type == VelocityBoundaryType::FixedValue)
        {
            *out = hb[fi].dot(f.area);
            continue;
        }
        let hf = if f.boundary() {
            hb[fi]
        } else {
            linear_vector(m, f, h)?
        };
        let d = linear_scalar(m, f, mobility)?;
        let mut pressure_flux = 0.0;
        if let Some(n) = f.neighbour {
            pressure_flux = d * geom.orthogonal_metric * (p[n] - p[f.owner]);
            if include_nonorthogonal {
                pressure_flux += d * linear_vector(m, f, &gp)?.dot(geom.nonorthogonal_area);
            }
        } else if pressure_bc[f.patch.ok_or_else(|| {
            crate::SolverError::InvalidInput("boundary face without patch".into())
        })?]
        .boundary_type
            == PressureBoundaryType::FixedValue
        {
            pressure_flux = d * geom.orthogonal_metric * (pb[fi] - p[f.owner]);
            if include_nonorthogonal {
                pressure_flux += d * gp[f.owner].dot(geom.nonorthogonal_area);
            }
        }
        *out = hf.dot(f.area) - pressure_flux;
    }
    finite_slice(&r, "Rhie-Chow flux")?;
    Ok(r)
}
pub fn flux_divergence(m: &PolyMesh, flux: &[f64]) -> Result<Vec<f64>> {
    require(flux.len() == m.face_count(), "face flux size")?;
    finite_slice(flux, "face flux")?;
    let mut d = vec![0.0; m.cell_count()];
    for (c, out) in d.iter_mut().enumerate() {
        *out = m.cell_faces()[c]
            .iter()
            .map(|&fi| {
                if m.faces()[fi].owner == c {
                    flux[fi]
                } else {
                    -flux[fi]
                }
            })
            .sum::<f64>()
            / m.cells()[c].volume;
    }
    finite_slice(&d, "flux divergence")?;
    Ok(d)
}
pub fn flux_divergence_l2(m: &PolyMesh, flux: &[f64]) -> Result<f64> {
    let d = flux_divergence(m, flux)?;
    Ok((d.iter().map(|v| v * v).sum::<f64>() / d.len() as f64).sqrt())
}
pub fn flux_divergence_max(m: &PolyMesh, flux: &[f64]) -> Result<f64> {
    Ok(flux_divergence(m, flux)?
        .iter()
        .map(|v| v.abs())
        .fold(0.0, f64::max))
}
