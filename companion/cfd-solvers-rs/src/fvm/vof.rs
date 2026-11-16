//! Volume-of-fluid transport, CSF surface tension and contact-angle normals, ported
//! from upstream `cfd/solvers/fvm/vof` (MIT, cfd_solvers contributors).
use super::{
    gauss_divergence_vector, interpolate_scalar_to_faces, least_squares_gradient_scalar,
    FaceInterpolationScheme, PolyMesh,
};
use crate::{core::Vec3, require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VofConfig {
    pub dt: f64,
    /// Blend from bounded-linear face values (0) toward bounded upwind (1).
    pub compression: f64,
}

impl Default for VofConfig {
    fn default() -> Self {
        Self {
            dt: 1e-3,
            compression: 1.0,
        }
    }
}

/// Explicit, conservative, bounded volume-fraction transport on prescribed face fluxes.
pub struct VofTransport {
    mesh: PolyMesh,
    config: VofConfig,
    alpha: Vec<f64>,
    flux: Vec<f64>,
}

impl VofTransport {
    pub fn new(mesh: PolyMesh, config: VofConfig) -> Result<Self> {
        require(
            config.dt.is_finite()
                && config.dt > 0.0
                && config.compression.is_finite()
                && config.compression >= 0.0,
            "invalid VOF controls",
        )?;
        let (n, f) = (mesh.cell_count(), mesh.face_count());
        Ok(Self {
            mesh,
            config,
            alpha: vec![0.0; n],
            flux: vec![0.0; f],
        })
    }
    pub fn initialize(&mut self, alpha: f64) -> Result<()> {
        require((0.0..=1.0).contains(&alpha), "invalid volume fraction")?;
        self.alpha.fill(alpha);
        Ok(())
    }
    pub fn initialize_profile(&mut self, f: impl Fn(Vec3) -> f64) {
        for (a, c) in self.alpha.iter_mut().zip(self.mesh.cells()) {
            *a = f(c.center).clamp(0.0, 1.0);
        }
    }
    pub fn set_face_flux(&mut self, flux: Vec<f64>) -> Result<()> {
        require(
            flux.len() == self.mesh.face_count() && flux.iter().all(|v| v.is_finite()),
            "VOF flux size",
        )?;
        self.flux = flux;
        Ok(())
    }
    pub fn step(&mut self) -> Result<()> {
        let mut face = interpolate_scalar_to_faces(
            &self.mesh,
            &self.alpha,
            FaceInterpolationScheme::BoundedLinear,
            &self.flux,
            &[],
        )?;
        if self.config.compression > 0.0 {
            let upwind = interpolate_scalar_to_faces(
                &self.mesh,
                &self.alpha,
                FaceInterpolationScheme::Upwind,
                &self.flux,
                &[],
            )?;
            let blend = self.config.compression.clamp(0.0, 1.0);
            for (f, u) in face.iter_mut().zip(upwind) {
                *f = ((1.0 - blend) * *f + blend * u).clamp(0.0, 1.0);
            }
        }
        let mut next = self.alpha.clone();
        for (c, value) in next.iter_mut().enumerate() {
            let mut sum = 0.0;
            for &fi in &self.mesh.cell_faces()[c] {
                let f = self.mesh.faces()[fi];
                let oriented = if f.owner == c {
                    self.flux[fi]
                } else {
                    -self.flux[fi]
                };
                sum += oriented * face[fi];
            }
            *value = (self.alpha[c] - self.config.dt * sum / self.mesh.cells()[c].volume)
                .clamp(0.0, 1.0);
        }
        self.alpha = next;
        Ok(())
    }
    pub fn alpha(&self) -> &[f64] {
        &self.alpha
    }
    pub fn volume(&self) -> f64 {
        self.alpha
            .iter()
            .zip(self.mesh.cells())
            .map(|(a, c)| a * c.volume)
            .sum()
    }
    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
}

/// Continuum-surface-force body force `-sigma kappa grad(alpha)`, curvature from the
/// divergence of the unit interface normal.
pub fn csf_surface_tension(mesh: &PolyMesh, alpha: &[f64], sigma: f64) -> Result<Vec<Vec3>> {
    require(
        alpha.len() == mesh.cell_count() && sigma.is_finite() && sigma >= 0.0,
        "invalid CSF inputs",
    )?;
    let grad = least_squares_gradient_scalar(mesh, alpha, &[])?;
    let normal: Vec<Vec3> = grad
        .iter()
        .map(|g| {
            if g.norm() > 1e-14 {
                *g / g.norm()
            } else {
                Vec3::ZERO
            }
        })
        .collect();
    let div = gauss_divergence_vector(mesh, &normal, &[])?;
    Ok(grad
        .iter()
        .zip(div)
        .map(|(g, d)| *g * (-sigma * d))
        .collect())
}

/// Interface normal at a wall that makes `angle` with the wall normal, keeping the
/// tangential direction of the given interface normal.
pub fn contact_angle_normal(interface: Vec3, wall: Vec3, angle: f64) -> Result<Vec3> {
    let (wn, inn) = (wall.norm(), interface.norm());
    require(
        wn > 0.0 && inn > 0.0 && (0.0..=std::f64::consts::PI).contains(&angle),
        "invalid contact angle",
    )?;
    let (nw, ni) = (wall / wn, interface / inn);
    let mut tangent = ni - nw * ni.dot(nw);
    if tangent.norm() < 1e-14 {
        let seed = if nw.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        tangent = seed - nw * seed.dot(nw);
    }
    let tangent = tangent / tangent.norm();
    Ok(nw * angle.cos() + tangent * angle.sin())
}
