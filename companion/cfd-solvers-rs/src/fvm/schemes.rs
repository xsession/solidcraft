use super::{
    gauss_gradient_scalar_into, least_squares_gradient_scalar_into,
    operators::{linear_scalar, linear_vector, scalar_sizes, vector_sizes},
    PolyMesh,
};
use crate::{core::Vec3, finite_slice, positive, require, Result, SolverError};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FaceInterpolationScheme {
    Linear,
    Upwind,
    BoundedLinear,
    Minmod,
    VanLeer,
}
impl Default for FaceInterpolationScheme {
    fn default() -> Self {
        Self::Linear
    }
}
pub fn parse_face_interpolation_scheme(name: &str) -> Result<FaceInterpolationScheme> {
    Ok(match name {
        "linear" => FaceInterpolationScheme::Linear,
        "upwind" => FaceInterpolationScheme::Upwind,
        "boundedLinear" | "bounded_linear" => FaceInterpolationScheme::BoundedLinear,
        "minmod" => FaceInterpolationScheme::Minmod,
        "vanLeer" | "van_leer" => FaceInterpolationScheme::VanLeer,
        _ => {
            return Err(SolverError::InvalidInput(format!(
                "unknown interpolation scheme {name}"
            )))
        }
    })
}
pub fn face_interpolation_scheme_name(s: FaceInterpolationScheme) -> &'static str {
    match s {
        FaceInterpolationScheme::Linear => "linear",
        FaceInterpolationScheme::Upwind => "upwind",
        FaceInterpolationScheme::BoundedLinear => "boundedLinear",
        FaceInterpolationScheme::Minmod => "minmod",
        FaceInterpolationScheme::VanLeer => "vanLeer",
    }
}
#[derive(Clone, Debug, Default)]
pub struct ScalarSchemeWorkspace {
    pub gradient: Vec<Vec3>,
    pub minimum: Vec<f64>,
    pub maximum: Vec<f64>,
    pub limiter: Vec<f64>,
    pub face_values: Vec<f64>,
}
impl ScalarSchemeWorkspace {
    fn prepare_cells(&mut self, m: &PolyMesh) {
        self.gradient.resize(m.cell_count(), Vec3::ZERO);
        self.minimum.resize(m.cell_count(), 0.0);
        self.maximum.resize(m.cell_count(), 0.0);
        self.limiter.resize(m.cell_count(), 0.0);
    }
    pub fn resize(&mut self, m: &PolyMesh) {
        self.prepare_cells(m);
        self.face_values.resize(m.face_count(), 0.0);
    }
    pub fn allocated_bytes(&self) -> usize {
        self.gradient.capacity() * std::mem::size_of::<Vec3>()
            + (self.minimum.capacity()
                + self.maximum.capacity()
                + self.limiter.capacity()
                + self.face_values.capacity())
                * std::mem::size_of::<f64>()
    }
}
pub(crate) fn flux_sizes(m: &PolyMesh, flux: &[f64]) -> Result<()> {
    require(flux.len() == m.face_count(), "FVM face flux size")?;
    finite_slice(flux, "face flux")
}
fn bounded_linear(
    m: &PolyMesh,
    v: &[f64],
    flux: &[f64],
    boundary: &[f64],
    out: &mut [f64],
    w: &mut ScalarSchemeWorkspace,
) -> Result<()> {
    gauss_gradient_scalar_into(m, v, &mut w.gradient, boundary)?;
    for c in 0..m.cell_count() {
        let (mut lo, mut hi) = (v[c], v[c]);
        for &fi in &m.cell_faces()[c] {
            let f = m.faces()[fi];
            let other = if let Some(n) = f.neighbour {
                v[if f.owner == c { n } else { f.owner }]
            } else if !boundary.is_empty() {
                boundary[fi]
            } else {
                v[c]
            };
            lo = lo.min(other);
            hi = hi.max(other);
        }
        w.minimum[c] = lo;
        w.maximum[c] = hi;
    }
    for c in 0..m.cell_count() {
        let mut a: f64 = 1.0;
        for &fi in &m.cell_faces()[c] {
            let d = w.gradient[c].dot(m.faces()[fi].center - m.cells()[c].center);
            if d > 0.0 {
                let room = w.maximum[c] - v[c];
                a = a.min(if room <= 0.0 { 0.0 } else { room / d });
            } else if d < 0.0 {
                let room = w.minimum[c] - v[c];
                a = a.min(if room >= 0.0 { 0.0 } else { room / d });
            }
        }
        w.limiter[c] = a.clamp(0.0, 1.0);
    }
    for (fi, r) in out.iter_mut().enumerate() {
        let f = m.faces()[fi];
        if f.boundary() && flux[fi] < 0.0 && !boundary.is_empty() {
            *r = boundary[fi];
            continue;
        }
        let up = if flux[fi] >= 0.0 {
            f.owner
        } else {
            f.neighbour.unwrap_or(f.owner)
        };
        let reconstructed =
            v[up] + w.limiter[up] * w.gradient[up].dot(f.center - m.cells()[up].center);
        *r = reconstructed.clamp(w.minimum[up], w.maximum[up]);
    }
    Ok(())
}
pub fn interpolate_scalar_to_faces_into(
    m: &PolyMesh,
    v: &[f64],
    out: &mut [f64],
    w: &mut ScalarSchemeWorkspace,
    scheme: FaceInterpolationScheme,
    flux: &[f64],
    boundary: &[f64],
) -> Result<()> {
    scalar_sizes(m, v, boundary)?;
    require(
        out.len() == m.face_count(),
        "face interpolation output size",
    )?;
    if scheme != FaceInterpolationScheme::Linear {
        flux_sizes(m, flux)?;
    }
    if scheme == FaceInterpolationScheme::BoundedLinear {
        w.prepare_cells(m);
        bounded_linear(m, v, flux, boundary, out, w)?;
        return finite_slice(out, "bounded face values");
    }
    let tvd = matches!(
        scheme,
        FaceInterpolationScheme::Minmod | FaceInterpolationScheme::VanLeer
    );
    if tvd {
        w.prepare_cells(m);
        least_squares_gradient_scalar_into(m, v, &mut w.gradient, boundary)?;
    }
    for (fi, r) in out.iter_mut().enumerate() {
        let f = m.faces()[fi];
        let Some(n) = f.neighbour else {
            *r = if !boundary.is_empty()
                && (scheme == FaceInterpolationScheme::Linear || flux[fi] < 0.0)
            {
                boundary[fi]
            } else {
                v[f.owner]
            };
            continue;
        };
        if scheme == FaceInterpolationScheme::Linear {
            *r = linear_scalar(m, f, v)?;
            continue;
        }
        let (up, down) = if flux[fi] >= 0.0 {
            (f.owner, n)
        } else {
            (n, f.owner)
        };
        *r = v[up];
        let delta = v[down] - v[up];
        if !tvd || delta.abs() < 1e-30 {
            continue;
        }
        let d = m.cells()[down].center - m.cells()[up].center;
        let ratio = 2.0 * w.gradient[up].dot(d) / delta - 1.0;
        let psi = if scheme == FaceInterpolationScheme::Minmod {
            ratio.clamp(0.0, 1.0)
        } else if ratio > 0.0 {
            2.0 / (1.0 + 1.0 / ratio)
        } else {
            0.0
        };
        let d2 = d.dot(d);
        positive(d2, "TVD neighbour distance")?;
        let fraction = ((f.center - m.cells()[up].center).dot(d) / d2).clamp(0.0, 1.0);
        *r = (v[up] + fraction * psi * delta).clamp(v[up].min(v[down]), v[up].max(v[down]));
    }
    finite_slice(out, "face interpolation")
}
pub fn interpolate_scalar_to_faces(
    m: &PolyMesh,
    v: &[f64],
    scheme: FaceInterpolationScheme,
    flux: &[f64],
    boundary: &[f64],
) -> Result<Vec<f64>> {
    let mut out = vec![0.0; m.face_count()];
    interpolate_scalar_to_faces_into(
        m,
        v,
        &mut out,
        &mut ScalarSchemeWorkspace::default(),
        scheme,
        flux,
        boundary,
    )?;
    Ok(out)
}
pub fn interpolate_vector_to_faces(
    m: &PolyMesh,
    v: &[Vec3],
    scheme: FaceInterpolationScheme,
    flux: &[f64],
    boundary: &[Vec3],
) -> Result<Vec<Vec3>> {
    vector_sizes(m, v, boundary)?;
    if scheme != FaceInterpolationScheme::Linear {
        flux_sizes(m, flux)?;
    }
    let mut result = vec![Vec3::ZERO; m.face_count()];
    if matches!(
        scheme,
        FaceInterpolationScheme::Linear | FaceInterpolationScheme::Upwind
    ) {
        for (fi, out) in result.iter_mut().enumerate() {
            let f = m.faces()[fi];
            *out = if f.boundary() {
                if !boundary.is_empty()
                    && (scheme == FaceInterpolationScheme::Linear || flux[fi] < 0.0)
                {
                    boundary[fi]
                } else {
                    v[f.owner]
                }
            } else if scheme == FaceInterpolationScheme::Upwind {
                v[if flux[fi] >= 0.0 {
                    f.owner
                } else {
                    f.neighbour.unwrap_or(f.owner)
                }]
            } else {
                linear_vector(m, f, v)?
            };
        }
    } else {
        let mut components = vec![];
        for axis in 0..3 {
            let cells: Vec<_> = v.iter().map(|p| p.component(axis)).collect();
            let b: Vec<_> = boundary.iter().map(|p| p.component(axis)).collect();
            components.push(interpolate_scalar_to_faces(m, &cells, scheme, flux, &b)?);
        }
        for i in 0..result.len() {
            result[i] = Vec3::new(components[0][i], components[1][i], components[2][i]);
        }
    }
    Ok(result)
}
pub fn face_flux_from_velocity(m: &PolyMesh, u: &[Vec3], boundary: &[Vec3]) -> Result<Vec<f64>> {
    let v = interpolate_vector_to_faces(m, u, FaceInterpolationScheme::Linear, &[], boundary)?;
    Ok(v.iter()
        .zip(m.faces())
        .map(|(v, f)| v.dot(f.area))
        .collect())
}
pub fn convective_divergence_scalar_into(
    m: &PolyMesh,
    v: &[f64],
    flux: &[f64],
    out: &mut [f64],
    w: &mut ScalarSchemeWorkspace,
    scheme: FaceInterpolationScheme,
    boundary: &[f64],
) -> Result<()> {
    flux_sizes(m, flux)?;
    require(out.len() == m.cell_count(), "convection output size")?;
    let mut faces = std::mem::take(&mut w.face_values);
    faces.resize(m.face_count(), 0.0);
    let result = interpolate_scalar_to_faces_into(m, v, &mut faces, w, scheme, flux, boundary);
    if result.is_ok() {
        for (c, r) in out.iter_mut().enumerate() {
            *r = m.cell_faces()[c]
                .iter()
                .map(|&fi| {
                    (if m.faces()[fi].owner == c {
                        flux[fi]
                    } else {
                        -flux[fi]
                    }) * faces[fi]
                })
                .sum::<f64>()
                / m.cells()[c].volume;
        }
    }
    w.face_values = faces;
    result?;
    finite_slice(out, "convective divergence")
}
pub fn convective_divergence_scalar(
    m: &PolyMesh,
    v: &[f64],
    flux: &[f64],
    scheme: FaceInterpolationScheme,
    boundary: &[f64],
) -> Result<Vec<f64>> {
    let mut out = vec![0.0; m.cell_count()];
    convective_divergence_scalar_into(
        m,
        v,
        flux,
        &mut out,
        &mut ScalarSchemeWorkspace::default(),
        scheme,
        boundary,
    )?;
    Ok(out)
}
pub fn convective_divergence_vector(
    m: &PolyMesh,
    v: &[Vec3],
    flux: &[f64],
    scheme: FaceInterpolationScheme,
    boundary: &[Vec3],
) -> Result<Vec<Vec3>> {
    flux_sizes(m, flux)?;
    let faces = interpolate_vector_to_faces(m, v, scheme, flux, boundary)?;
    let mut out = vec![Vec3::ZERO; m.cell_count()];
    for (c, r) in out.iter_mut().enumerate() {
        for &fi in &m.cell_faces()[c] {
            *r += faces[fi]
                * (if m.faces()[fi].owner == c {
                    flux[fi]
                } else {
                    -flux[fi]
                });
        }
        *r = *r / m.cells()[c].volume;
    }
    Ok(out)
}
pub fn courant_limited_timestep(
    m: &PolyMesh,
    flux: &[f64],
    target: f64,
    maximum_dt: f64,
) -> Result<f64> {
    flux_sizes(m, flux)?;
    positive(target, "target Courant")?;
    positive(maximum_dt, "maximum dt")?;
    let mut dt = maximum_dt;
    for c in 0..m.cell_count() {
        let sum: f64 = m.cell_faces()[c].iter().map(|&fi| flux[fi].abs()).sum();
        if sum > 0.0 {
            dt = dt.min(target * m.cells()[c].volume / sum);
        }
    }
    positive(dt, "Courant timestep")?;
    Ok(dt)
}
