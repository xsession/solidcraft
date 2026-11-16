//! Asymptotic high-frequency methods: physical-optics currents and RCS, and a
//! shooting-and-bouncing-rays tracer over triangle surfaces.
use super::SurfaceTriangle;
use crate::{
    core::Vec3,
    optics::{reflect, Ray},
    require, Complex64 as C, Result, SolverError,
};
use std::f64::consts::PI;

pub type ComplexPoint3 = [C; 3];

const SPEED_OF_LIGHT: f64 = 299792458.0;
const VACUUM_IMPEDANCE: f64 = 376.730313668;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhysicalOpticsConfig {
    pub frequency_hz: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub incident_direction: Vec3,
    pub incident_electric: Vec3,
    pub observation_direction: Vec3,
    pub observation_distance_m: f64,
}

impl Default for PhysicalOpticsConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 0.0,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            incident_direction: Vec3::new(0.0, 0.0, -1.0),
            incident_electric: Vec3::new(1.0, 0.0, 0.0),
            observation_direction: Vec3::new(0.0, 0.0, 1.0),
            observation_distance_m: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PhysicalOpticsCurrent {
    pub centroid: Vec3,
    pub normal: Vec3,
    pub current: ComplexPoint3,
    pub area_m2: f64,
    pub illuminated: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PhysicalOpticsResult {
    pub currents: Vec<PhysicalOpticsCurrent>,
    pub far_field: ComplexPoint3,
    pub illuminated_area_m2: f64,
    pub rcs_m2: f64,
    pub wave_number_per_m: f64,
    pub wave_impedance_ohm: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SbrConfig {
    pub frequency_hz: f64,
    pub max_bounces: usize,
    /// Accepted for API parity; upstream does not use it.
    pub minimum_power: f64,
    pub intersection_tolerance: f64,
}

impl Default for SbrConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 0.0,
            max_bounces: 4,
            minimum_power: 1e-12,
            intersection_tolerance: 1e-10,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SbrHit {
    pub triangle: usize,
    pub distance_m: f64,
    pub point: Vec3,
    pub normal: Vec3,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SbrResult {
    pub escaped_ray: Ray,
    pub hits: Vec<SbrHit>,
    pub optical_path_m: f64,
    pub phase: C,
    pub escaped: bool,
}

fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}
fn scale(a: Vec3, s: f64) -> Vec3 {
    Vec3::new(a.x * s, a.y * s, a.z * s)
}
fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}
fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: Vec3) -> Result<Vec3> {
    let l = norm(a);
    require(l > 0.0 && l.is_finite(), "zero asymptotic direction")?;
    Ok(scale(a, 1.0 / l))
}
fn ccross_pc(a: Vec3, b: ComplexPoint3) -> ComplexPoint3 {
    [
        a.y * b[2] - a.z * b[1],
        a.z * b[0] - a.x * b[2],
        a.x * b[1] - a.y * b[0],
    ]
}
fn ccross_cp(a: ComplexPoint3, b: Vec3) -> ComplexPoint3 {
    [
        a[1] * b.z - a[2] * b.y,
        a[2] * b.x - a[0] * b.z,
        a[0] * b.y - a[1] * b.x,
    ]
}

fn triangle_normal(t: &SurfaceTriangle) -> Result<(Vec3, f64)> {
    let n = cross(sub(t.b, t.a), sub(t.c, t.a));
    let twice = norm(n);
    let area = 0.5 * twice;
    require(
        area > 0.0 && area.is_finite(),
        "degenerate asymptotic triangle",
    )?;
    Ok((scale(n, 1.0 / twice), area))
}

/// Moller-Trumbore intersection; `None` for a miss.
fn ray_triangle(ray: &Ray, t: &SurfaceTriangle, tol: f64) -> Result<Option<(f64, Vec3, Vec3)>> {
    let direction = unit(ray.direction)?;
    let (e1, e2) = (sub(t.b, t.a), sub(t.c, t.a));
    let pv = cross(direction, e2);
    let det = dot(e1, pv);
    if det.is_nan() || det.abs() <= tol {
        return Ok(None);
    }
    let inv = 1.0 / det;
    let tv = sub(ray.origin, t.a);
    let u = inv * dot(tv, pv);
    if u.is_nan() || u < -tol || u > 1.0 + tol {
        return Ok(None);
    }
    let q = cross(tv, e1);
    let v = inv * dot(direction, q);
    if v.is_nan() || v < -tol || u + v > 1.0 + tol {
        return Ok(None);
    }
    let distance = inv * dot(e2, q);
    if !(distance > tol && distance.is_finite()) {
        return Ok(None);
    }
    let point = add(ray.origin, scale(direction, distance));
    Ok(Some((distance, point, triangle_normal(t)?.0)))
}

/// Lit-face PO currents `J = 2 n x H_inc`, far field by the radiation integral and monostatic or
/// bistatic RCS.
pub fn solve_physical_optics(
    triangles: &[SurfaceTriangle],
    config: &PhysicalOpticsConfig,
) -> Result<PhysicalOpticsResult> {
    let c = config;
    require(
        c.frequency_hz > 0.0
            && c.relative_permittivity > 0.0
            && c.relative_permeability > 0.0
            && c.observation_distance_m > 0.0,
        "invalid physical-optics configuration",
    )?;
    let incident = unit(c.incident_direction)?;
    let observation = unit(c.observation_direction)?;
    let k = 2.0 * PI * c.frequency_hz * (c.relative_permittivity * c.relative_permeability).sqrt()
        / SPEED_OF_LIGHT;
    let eta = VACUUM_IMPEDANCE * (c.relative_permeability / c.relative_permittivity).sqrt();
    let transverse = sub(
        c.incident_electric,
        scale(incident, dot(c.incident_electric, incident)),
    );
    require(
        norm(transverse) > 0.0,
        "incident electric field is longitudinal",
    )?;
    let magnetic = scale(cross(incident, transverse), 1.0 / eta);
    let mut result = PhysicalOpticsResult {
        wave_number_per_m: k,
        wave_impedance_ohm: eta,
        currents: Vec::with_capacity(triangles.len()),
        ..Default::default()
    };
    let zero = C::new(0.0, 0.0);
    let mut integral = [zero; 3];
    for t in triangles {
        let (normal, area) = triangle_normal(t)?;
        let centroid = scale(add(add(t.a, t.b), t.c), 1.0 / 3.0);
        let illuminated = dot(normal, incident) < 0.0;
        let incident_phase = C::new(0.0, -k * dot(incident, centroid)).exp();
        let mut current = [zero; 3];
        if illuminated {
            let two = C::new(2.0, 0.0);
            let sc = ccross_pc(
                normal,
                [two * magnetic.x, two * magnetic.y, two * magnetic.z],
            );
            current = [
                sc[0] * incident_phase,
                sc[1] * incident_phase,
                sc[2] * incident_phase,
            ];
            result.illuminated_area_m2 += area;
            let ph = C::new(0.0, k * dot(observation, centroid)).exp();
            for a in 0..3 {
                integral[a] += current[a] * area * ph;
            }
        }
        result.currents.push(PhysicalOpticsCurrent {
            centroid,
            normal,
            current,
            area_m2: area,
            illuminated,
        });
    }
    let projected = ccross_cp(ccross_pc(observation, integral), observation);
    let prefactor = C::new(0.0, k * eta / (4.0 * PI * c.observation_distance_m));
    result.far_field = [
        prefactor * projected[0],
        prefactor * projected[1],
        prefactor * projected[2],
    ];
    let e2 = norm(transverse);
    result.rcs_m2 = 4.0
        * PI
        * c.observation_distance_m
        * c.observation_distance_m
        * result.far_field.iter().map(|z| z.norm_sqr()).sum::<f64>()
        / e2;
    Ok(result)
}

/// Specular ray tracing over triangles with accumulated optical path and phase.
pub fn trace_shooting_bouncing_rays(
    triangles: &[SurfaceTriangle],
    mut seed: Ray,
    config: &SbrConfig,
) -> Result<SbrResult> {
    require(
        config.frequency_hz > 0.0
            && config.max_bounces > 0
            && config.minimum_power >= 0.0
            && config.intersection_tolerance > 0.0,
        "invalid shooting-and-bouncing-rays configuration",
    )?;
    let l = norm(seed.direction);
    if l == 0.0 || l.is_nan() {
        return Err(SolverError::InvalidInput(
            "cannot normalize zero vector".into(),
        ));
    }
    seed.direction = Vec3::new(
        seed.direction.x / l,
        seed.direction.y / l,
        seed.direction.z / l,
    );
    let mut result = SbrResult {
        escaped_ray: seed,
        hits: vec![],
        optical_path_m: 0.0,
        phase: C::new(1.0, 0.0),
        escaped: false,
    };
    let k = 2.0 * PI * config.frequency_hz / SPEED_OF_LIGHT;
    for _ in 0..config.max_bounces {
        let mut best: Option<(usize, f64, Vec3, Vec3)> = None;
        for (i, t) in triangles.iter().enumerate() {
            if let Some((d, pnt, n)) =
                ray_triangle(&result.escaped_ray, t, config.intersection_tolerance)?
            {
                if best.is_none_or(|b| d < b.1) {
                    best = Some((i, d, pnt, n));
                }
            }
        }
        let Some((tri, d, pnt, n)) = best else {
            result.escaped = true;
            break;
        };
        result.optical_path_m += d;
        result.hits.push(SbrHit {
            triangle: tri,
            distance_m: d,
            point: pnt,
            normal: n,
        });
        let reflected = reflect(result.escaped_ray.direction, n)?;
        result.escaped_ray.origin =
            add(pnt, scale(reflected, config.intersection_tolerance * 10.0));
        result.escaped_ray.direction = reflected;
    }
    result.phase = C::new(0.0, k * result.optical_path_m).exp();
    if result.escaped_ray.wavelength_nm <= 0.0 {
        result.escaped_ray.wavelength_nm = 550.0;
    }
    Ok(result)
}
