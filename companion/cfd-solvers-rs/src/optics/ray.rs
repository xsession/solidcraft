use crate::{core::Vec3, positive, require, Result};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
    pub wavelength_nm: f64,
}
impl Default for Ray {
    fn default() -> Self {
        Self {
            origin: Vec3::ZERO,
            direction: Vec3::new(0.0, 0.0, 1.0),
            wavelength_nm: 550.0,
        }
    }
}
pub fn reflect(incident: Vec3, normal: Vec3) -> Result<Vec3> {
    let i = incident.normalized()?;
    let mut n = normal.normalized()?;
    if i.dot(n) > 0.0 {
        n = -n;
    }
    (i - n * (2.0 * i.dot(n))).normalized()
}
pub fn refract(incident: Vec3, normal: Vec3, n_from: f64, n_to: f64) -> Result<Option<Vec3>> {
    positive(n_from, "incident index")?;
    positive(n_to, "transmitted index")?;
    let i = incident.normalized()?;
    let mut n = normal.normalized()?;
    let mut cosine = -i.dot(n);
    if cosine < 0.0 {
        n = -n;
        cosine = -i.dot(n);
    }
    cosine = cosine.clamp(0.0, 1.0);
    let eta = n_from / n_to;
    let k = 1.0 - eta * eta * (1.0 - cosine * cosine);
    if k < 0.0 {
        return Ok(None);
    }
    Ok(Some(
        (i * eta + n * (eta * cosine - k.sqrt())).normalized()?,
    ))
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SurfaceType {
    Plane,
    Sphere,
    Conic,
    EvenAsphere,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SequentialSurface {
    pub surface_type: SurfaceType,
    pub vertex_z: f64,
    pub radius: f64,
    pub aperture_radius: f64,
    pub refractive_index_after: f64,
    pub conic_constant: f64,
    pub even_coefficients: [f64; 4],
}
impl Default for SequentialSurface {
    fn default() -> Self {
        Self {
            surface_type: SurfaceType::Plane,
            vertex_z: 0.0,
            radius: 0.0,
            aperture_radius: 1e30,
            refractive_index_after: 1.0,
            conic_constant: 0.0,
            even_coefficients: [0.0; 4],
        }
    }
}
impl SequentialSurface {
    pub fn validate(self) -> Result<()> {
        positive(self.aperture_radius, "aperture")?;
        positive(self.refractive_index_after, "surface index")?;
        require(
            self.vertex_z.is_finite()
                && self.conic_constant.is_finite()
                && self.even_coefficients.iter().all(|v| v.is_finite()),
            "nonfinite optical surface",
        )?;
        if self.surface_type != SurfaceType::Plane {
            require(
                self.radius.is_finite() && self.radius != 0.0,
                "invalid surface radius",
            )?;
        }
        Ok(())
    }
    fn sag_and_slope(self, r: f64) -> Option<(f64, f64)> {
        if self.surface_type == SurfaceType::Plane {
            return Some((0.0, 0.0));
        }
        let c = 1.0 / self.radius;
        let k = if self.surface_type == SurfaceType::Sphere {
            0.0
        } else {
            self.conic_constant
        };
        let rad = 1.0 - (1.0 + k) * c * c * r * r;
        if rad <= 0.0 {
            return None;
        }
        let root = rad.sqrt();
        let mut sag = c * r * r / (1.0 + root);
        let mut slope = c * r / root;
        if self.surface_type == SurfaceType::EvenAsphere {
            for i in 0..4 {
                let e = 4 + 2 * i as i32;
                sag += self.even_coefficients[i] * r.powi(e);
                slope += e as f64 * self.even_coefficients[i] * r.powi(e - 1);
            }
        }
        if sag.is_finite() && slope.is_finite() {
            Some((sag, slope))
        } else {
            None
        }
    }
    pub fn sag(self, r: f64) -> Result<f64> {
        self.validate()?;
        require(r.is_finite() && r >= 0.0, "invalid sag radius")?;
        self.sag_and_slope(r).map(|v| v.0).ok_or_else(|| {
            crate::SolverError::InvalidInput("surface sag outside real conic".into())
        })
    }
    fn intersect(self, ray: Ray) -> Result<Option<(Vec3, Vec3)>> {
        let eps = 1e-12;
        let (hit, normal);
        if self.surface_type == SurfaceType::Plane {
            if ray.direction.z.abs() < eps {
                return Ok(None);
            }
            let t = (self.vertex_z - ray.origin.z) / ray.direction.z;
            if t <= eps {
                return Ok(None);
            }
            hit = ray.origin + ray.direction * t;
            normal = Vec3::new(0.0, 0.0, -1.0);
        } else if matches!(
            self.surface_type,
            SurfaceType::Conic | SurfaceType::EvenAsphere
        ) {
            if ray.direction.z.abs() < eps {
                return Ok(None);
            }
            let mut t = (self.vertex_z - ray.origin.z) / ray.direction.z;
            let mut found = None;
            for _ in 0..50 {
                let h = ray.origin + ray.direction * t;
                let r = h.x.hypot(h.y);
                let Some((sag, slope)) = self.sag_and_slope(r) else {
                    return Ok(None);
                };
                let res = h.z - self.vertex_z - sag;
                let (sx, sy) = if r > 0.0 {
                    (slope * h.x / r, slope * h.y / r)
                } else {
                    (0.0, 0.0)
                };
                let n = Vec3::new(sx, sy, -1.0).normalized()?;
                if res.abs() < 1e-12 * (1.0 + sag.abs()) {
                    found = Some((h, n));
                    break;
                }
                let derivative = ray.direction.z - sx * ray.direction.x - sy * ray.direction.y;
                if derivative.abs() < eps {
                    return Ok(None);
                }
                t -= res / derivative;
                if !t.is_finite() {
                    return Ok(None);
                }
            }
            if t <= eps {
                return Ok(None);
            }
            let Some((h, n)) = found else {
                return Ok(None);
            };
            hit = h;
            normal = n;
        } else {
            if self.radius.abs() < eps {
                return Ok(None);
            }
            let center = Vec3::new(0.0, 0.0, self.vertex_z + self.radius);
            let oc = ray.origin - center;
            let b = oc.dot(ray.direction);
            let c = oc.dot(oc) - self.radius * self.radius;
            let disc = b * b - c;
            if disc < 0.0 {
                return Ok(None);
            }
            let root = disc.sqrt();
            let (t0, t1) = (-b - root, -b + root);
            let t = if t0 > eps {
                t0
            } else if t1 > eps {
                t1
            } else {
                return Ok(None);
            };
            hit = ray.origin + ray.direction * t;
            normal = (hit - center).normalized()?;
        }
        Ok(
            if hit.x * hit.x + hit.y * hit.y <= self.aperture_radius * self.aperture_radius {
                Some((hit, normal))
            } else {
                None
            },
        )
    }
}
#[derive(Clone, Copy, Debug, Serialize)]
pub struct RayTraceResult {
    pub ray: Ray,
    pub valid: bool,
    pub surfaces_traced: usize,
}
#[derive(Clone, Debug)]
pub struct SequentialOpticalSystem {
    object_space_index: f64,
    surfaces: Vec<SequentialSurface>,
}
impl SequentialOpticalSystem {
    pub fn new(index: f64) -> Result<Self> {
        positive(index, "object-space index")?;
        Ok(Self {
            object_space_index: index,
            surfaces: vec![],
        })
    }
    pub fn add_surface(&mut self, s: SequentialSurface) -> Result<()> {
        s.validate()?;
        self.surfaces.push(s);
        Ok(())
    }
    pub fn surfaces(&self) -> &[SequentialSurface] {
        &self.surfaces
    }
    pub fn object_space_index(&self) -> f64 {
        self.object_space_index
    }
    pub fn trace(&self, mut ray: Ray) -> Result<RayTraceResult> {
        require(ray.origin.is_finite(), "nonfinite ray origin")?;
        positive(ray.wavelength_nm, "wavelength")?;
        ray.direction = ray.direction.normalized()?;
        let mut n = self.object_space_index;
        let mut count = 0;
        for &s in &self.surfaces {
            let Some((hit, normal)) = s.intersect(ray)? else {
                return Ok(RayTraceResult {
                    ray,
                    valid: false,
                    surfaces_traced: count,
                });
            };
            let Some(transmitted) = refract(ray.direction, normal, n, s.refractive_index_after)?
            else {
                return Ok(RayTraceResult {
                    ray,
                    valid: false,
                    surfaces_traced: count,
                });
            };
            ray.origin = hit;
            ray.direction = transmitted;
            n = s.refractive_index_after;
            count += 1;
        }
        Ok(RayTraceResult {
            ray,
            valid: true,
            surfaces_traced: count,
        })
    }
    pub fn trace_many(&self, rays: &[Ray]) -> Result<Vec<RayTraceResult>> {
        rays.iter().map(|&r| self.trace(r)).collect()
    }
}
pub fn rms_spot_radius_at_plane(rays: &[RayTraceResult], image_z: f64) -> Result<f64> {
    require(image_z.is_finite(), "image plane must be finite")?;
    let mut xy = vec![];
    for r in rays {
        if !r.valid || r.ray.direction.z.abs() <= 1e-14 {
            continue;
        }
        let t = (image_z - r.ray.origin.z) / r.ray.direction.z;
        if t < 0.0 {
            continue;
        }
        xy.push((
            r.ray.origin.x + t * r.ray.direction.x,
            r.ray.origin.y + t * r.ray.direction.y,
        ));
    }
    require(!xy.is_empty(), "no valid rays reach image plane")?;
    let n = xy.len() as f64;
    let cx = xy.iter().map(|p| p.0).sum::<f64>() / n;
    let cy = xy.iter().map(|p| p.1).sum::<f64>() / n;
    Ok((xy
        .iter()
        .map(|p| (p.0 - cx).powi(2) + (p.1 - cy).powi(2))
        .sum::<f64>()
        / n)
        .sqrt())
}
#[derive(Clone, Copy, Debug)]
pub struct ParaxialRay {
    pub height: f64,
    pub angle: f64,
    pub z: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct ParaxialTraceResult {
    pub ray: ParaxialRay,
    pub valid: bool,
    pub surfaces_traced: usize,
}
pub fn trace_paraxial(
    system: &SequentialOpticalSystem,
    mut ray: ParaxialRay,
    mut index: f64,
) -> Result<ParaxialTraceResult> {
    positive(index, "paraxial index")?;
    require(
        [ray.height, ray.angle, ray.z].iter().all(|v| v.is_finite()),
        "nonfinite paraxial ray",
    )?;
    let mut count = 0;
    for s in system.surfaces() {
        ray.height += (s.vertex_z - ray.z) * ray.angle;
        ray.z = s.vertex_z;
        if ray.height.abs() > s.aperture_radius {
            return Ok(ParaxialTraceResult {
                ray,
                valid: false,
                surfaces_traced: count,
            });
        }
        let n2 = s.refractive_index_after;
        let curvature = if s.surface_type == SurfaceType::Sphere {
            1.0 / s.radius
        } else {
            0.0
        };
        ray.angle = index / n2 * ray.angle - (n2 - index) / n2 * curvature * ray.height;
        index = n2;
        count += 1;
    }
    Ok(ParaxialTraceResult {
        ray,
        valid: true,
        surfaces_traced: count,
    })
}
pub fn paraxial_back_focal_distance(system: &SequentialOpticalSystem, index: f64) -> Result<f64> {
    require(
        !system.surfaces.is_empty(),
        "focal distance requires surfaces",
    )?;
    let r = trace_paraxial(
        system,
        ParaxialRay {
            height: 1e-3,
            angle: 0.0,
            z: system.surfaces[0].vertex_z,
        },
        index,
    )?;
    Ok(if !r.valid || r.ray.angle.abs() < 1e-18 {
        f64::INFINITY
    } else {
        -r.ray.height / r.ray.angle
    })
}
