use super::{decompose_face_area, Face, PolyMesh};
use crate::{core::Vec3, finite_slice, positive, require, Result};
pub(crate) fn scalar_sizes(m: &PolyMesh, v: &[f64], b: &[f64]) -> Result<()> {
    require(
        v.len() == m.cell_count() && (b.is_empty() || b.len() == m.face_count()),
        "FVM scalar field sizes",
    )?;
    finite_slice(v, "cell scalar")?;
    finite_slice(b, "boundary scalar")
}
pub(crate) fn vector_sizes(m: &PolyMesh, v: &[Vec3], b: &[Vec3]) -> Result<()> {
    require(
        v.len() == m.cell_count()
            && (b.is_empty() || b.len() == m.face_count())
            && v.iter().chain(b).all(|v| v.is_finite()),
        "FVM vector sizes/values",
    )
}
pub(crate) fn interpolation_weight(m: &PolyMesh, f: Face) -> Result<f64> {
    let Some(n) = f.neighbour else {
        return Ok(1.0);
    };
    let dof = (f.center - m.cells()[f.owner].center).norm();
    let dnf = (m.cells()[n].center - f.center).norm();
    positive(dof + dnf, "face interpolation distance")?;
    Ok(dnf / (dof + dnf))
}
pub(crate) fn linear_scalar(m: &PolyMesh, f: Face, v: &[f64]) -> Result<f64> {
    let Some(n) = f.neighbour else {
        return Ok(v[f.owner]);
    };
    let w = interpolation_weight(m, f)?;
    Ok(w * v[f.owner] + (1.0 - w) * v[n])
}
pub(crate) fn linear_vector(m: &PolyMesh, f: Face, v: &[Vec3]) -> Result<Vec3> {
    let Some(n) = f.neighbour else {
        return Ok(v[f.owner]);
    };
    let w = interpolation_weight(m, f)?;
    Ok(v[f.owner] * w + v[n] * (1.0 - w))
}
pub fn gauss_gradient_scalar(m: &PolyMesh, v: &[f64], b: &[f64]) -> Result<Vec<Vec3>> {
    let mut r = vec![Vec3::ZERO; m.cell_count()];
    gauss_gradient_scalar_into(m, v, &mut r, b)?;
    Ok(r)
}
pub fn gauss_gradient_scalar_into(
    m: &PolyMesh,
    v: &[f64],
    r: &mut [Vec3],
    b: &[f64],
) -> Result<()> {
    scalar_sizes(m, v, b)?;
    require(r.len() == m.cell_count(), "gradient output size")?;
    for (c, out) in r.iter_mut().enumerate() {
        let mut sum = Vec3::ZERO;
        for &fi in &m.cell_faces()[c] {
            let f = m.faces()[fi];
            let value = if f.boundary() && !b.is_empty() {
                b[fi]
            } else {
                linear_scalar(m, f, v)?
            };
            sum += f.outward_area(c) * value;
        }
        *out = sum / m.cells()[c].volume;
    }
    Ok(())
}
pub fn gauss_divergence_vector(m: &PolyMesh, v: &[Vec3], b: &[Vec3]) -> Result<Vec<f64>> {
    vector_sizes(m, v, b)?;
    let mut r = vec![0.0; m.cell_count()];
    for (c, out) in r.iter_mut().enumerate() {
        let mut sum = 0.0;
        for &fi in &m.cell_faces()[c] {
            let f = m.faces()[fi];
            let value = if f.boundary() && !b.is_empty() {
                b[fi]
            } else {
                linear_vector(m, f, v)?
            };
            sum += value.dot(f.outward_area(c));
        }
        *out = sum / m.cells()[c].volume;
    }
    finite_slice(&r, "divergence")?;
    Ok(r)
}
pub fn orthogonal_laplacian_scalar(
    m: &PolyMesh,
    v: &[f64],
    diffusivity: f64,
    b: &[f64],
) -> Result<Vec<f64>> {
    scalar_sizes(m, v, b)?;
    positive(diffusivity, "diffusivity")?;
    let mut r = vec![0.0; m.cell_count()];
    for (c, out) in r.iter_mut().enumerate() {
        let mut sum = 0.0;
        for &fi in &m.cell_faces()[c] {
            let f = m.faces()[fi];
            let sf = f.outward_area(c);
            let area = sf.norm();
            let n = sf / area;
            if let Some(other) = f.neighbour.map(|n| if f.owner == c { n } else { f.owner }) {
                let d = (m.cells()[other].center - m.cells()[c].center).dot(n).abs();
                positive(d, "internal face distance")?;
                sum += diffusivity * area * (v[other] - v[c]) / d;
            } else if !b.is_empty() {
                let d = (f.center - m.cells()[c].center).dot(n).abs();
                positive(d, "boundary face distance")?;
                sum += diffusivity * area * (b[fi] - v[c]) / d;
            }
        }
        *out = sum / m.cells()[c].volume;
    }
    finite_slice(&r, "Laplacian")?;
    Ok(r)
}
fn weighted_gradient_solve(a: [[f64; 3]; 3], rhs: [f64; 3]) -> Vec3 {
    let mut aug = [[0.0; 4]; 3];
    for i in 0..3 {
        aug[i][..3].copy_from_slice(&a[i]);
        aug[i][3] = rhs[i];
    }
    let trace = a[0][0].abs() + a[1][1].abs() + a[2][2].abs();
    let reg = (trace * 1e-14).max(1e-30);
    for i in 0..3 {
        aug[i][i] += reg;
    }
    for col in 0..3 {
        let mut p = col;
        for r in col + 1..3 {
            if aug[r][col].abs() > aug[p][col].abs() {
                p = r;
            }
        }
        if aug[p][col].abs() < 1e-30 {
            continue;
        }
        aug.swap(p, col);
        let inv = 1.0 / aug[col][col];
        for c in col..4 {
            aug[col][c] *= inv;
        }
        for r in 0..3 {
            if r == col {
                continue;
            }
            let f = aug[r][col];
            for c in col..4 {
                aug[r][c] -= f * aug[col][c];
            }
        }
    }
    Vec3::new(aug[0][3], aug[1][3], aug[2][3])
}
pub fn least_squares_gradient_scalar(m: &PolyMesh, v: &[f64], b: &[f64]) -> Result<Vec<Vec3>> {
    let mut r = vec![Vec3::ZERO; m.cell_count()];
    least_squares_gradient_scalar_into(m, v, &mut r, b)?;
    Ok(r)
}
pub fn least_squares_gradient_scalar_into(
    m: &PolyMesh,
    v: &[f64],
    r: &mut [Vec3],
    b: &[f64],
) -> Result<()> {
    scalar_sizes(m, v, b)?;
    require(r.len() == m.cell_count(), "least-squares output size")?;
    for (c, out) in r.iter_mut().enumerate() {
        let mut a = [[0.0; 3]; 3];
        let mut rhs = [0.0; 3];
        let center = m.cells()[c].center;
        for &fi in &m.cell_faces()[c] {
            let f = m.faces()[fi];
            let (d, dv) = if let Some(n) = f.neighbour {
                let other = if f.owner == c { n } else { f.owner };
                (m.cells()[other].center - center, v[other] - v[c])
            } else if !b.is_empty() {
                (f.center - center, b[fi] - v[c])
            } else {
                continue;
            };
            let d2 = d.dot(d);
            if d2 <= 0.0 {
                continue;
            }
            let w = 1.0 / d2;
            let d = d.as_array();
            for i in 0..3 {
                rhs[i] += w * d[i] * dv;
                for j in 0..3 {
                    a[i][j] += w * d[i] * d[j];
                }
            }
        }
        *out = weighted_gradient_solve(a, rhs);
    }
    require(r.iter().all(|v| v.is_finite()), "gradient became nonfinite")
}
pub fn corrected_face_normal_gradient_scalar(
    m: &PolyMesh,
    v: &[f64],
    b: &[f64],
) -> Result<Vec<f64>> {
    scalar_sizes(m, v, b)?;
    let g = least_squares_gradient_scalar(m, v, b)?;
    let mut r = vec![0.0; m.face_count()];
    for (fi, out) in r.iter_mut().enumerate() {
        let f = m.faces()[fi];
        let geom = decompose_face_area(m, fi)?;
        let (delta, gf) = if let Some(n) = f.neighbour {
            (v[n] - v[f.owner], (g[f.owner] + g[n]) * 0.5)
        } else if !b.is_empty() {
            (b[fi] - v[f.owner], g[f.owner])
        } else {
            continue;
        };
        *out = (geom.orthogonal_metric * delta + gf.dot(geom.nonorthogonal_area)) / f.area.norm();
    }
    finite_slice(&r, "normal gradient")?;
    Ok(r)
}
pub fn corrected_laplacian_scalar(
    m: &PolyMesh,
    v: &[f64],
    diffusivity: f64,
    b: &[f64],
) -> Result<Vec<f64>> {
    positive(diffusivity, "diffusivity")?;
    let sn = corrected_face_normal_gradient_scalar(m, v, b)?;
    let mut r = vec![0.0; m.cell_count()];
    for (c, out) in r.iter_mut().enumerate() {
        let mut sum = 0.0;
        for &fi in &m.cell_faces()[c] {
            let f = m.faces()[fi];
            if f.boundary() && b.is_empty() {
                continue;
            }
            sum += diffusivity * (if f.owner == c { 1.0 } else { -1.0 }) * f.area.norm() * sn[fi];
        }
        *out = sum / m.cells()[c].volume;
    }
    finite_slice(&r, "corrected Laplacian")?;
    Ok(r)
}
