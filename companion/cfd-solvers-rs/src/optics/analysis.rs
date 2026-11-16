//! Spot diagrams and ray fans, 1-D Fraunhofer/PSF/MTF, Zernike wavefronts, birefringence, Kerr
//! polarization, rigid transforms, field/wavelength/aperture models and polynomial freeform
//! surfaces, ported from upstream `cfd/optics/{analysis, diffraction, wavefront, birefringence,
//! nonlinear_polarization, transform, system_models, freeform}` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::{JonesVector, Ray, RayTraceResult};
use crate::{core::Vec3, require, Complex64 as C, Result, SolverError};
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpotPoint {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpotDiagram {
    pub points: Vec<SpotPoint>,
    pub centroid_x: f64,
    pub centroid_y: f64,
    pub rms_radius: f64,
}

fn plane_parameter(r: &RayTraceResult, z: f64) -> Option<f64> {
    if !r.valid || r.ray.direction.z.abs() < 1e-14 {
        return None;
    }
    let t = (z - r.ray.origin.z) / r.ray.direction.z;
    (t >= 0.0).then_some(t)
}

/// Valid rays intersected with the plane `z = image_z`, with centroid and RMS radius.
pub fn spot_diagram_at_plane(rays: &[RayTraceResult], image_z: f64) -> Result<SpotDiagram> {
    let mut d = SpotDiagram::default();
    for r in rays {
        if let Some(t) = plane_parameter(r, image_z) {
            d.points.push(SpotPoint {
                x: r.ray.origin.x + t * r.ray.direction.x,
                y: r.ray.origin.y + t * r.ray.direction.y,
            });
        }
    }
    require(!d.points.is_empty(), "no rays reach spot plane")?;
    for p in &d.points {
        d.centroid_x += p.x;
        d.centroid_y += p.y;
    }
    let n = d.points.len() as f64;
    d.centroid_x /= n;
    d.centroid_y /= n;
    for p in &d.points {
        let (x, y) = (p.x - d.centroid_x, p.y - d.centroid_y);
        d.rms_radius += x * x + y * y;
    }
    d.rms_radius = (d.rms_radius / n).sqrt();
    Ok(d)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RayFanPoint {
    pub pupil_coordinate: f64,
    pub transverse_aberration: f64,
}

/// Tangential fan: pupil coordinate `y0 / max|y0|` against image-plane `y - reference_x`.
pub fn tangential_ray_fan(
    rays: &[RayTraceResult],
    image_z: f64,
    reference_x: f64,
) -> Vec<RayFanPoint> {
    let mut maxp = rays
        .iter()
        .filter(|r| r.valid)
        .fold(0.0_f64, |m, r| m.max(r.ray.origin.y.abs()));
    if maxp == 0.0 {
        maxp = 1.0;
    }
    rays.iter()
        .filter_map(|r| {
            plane_parameter(r, image_z).map(|t| RayFanPoint {
                pupil_coordinate: r.ray.origin.y / maxp,
                transverse_aberration: r.ray.origin.y + t * r.ray.direction.y - reference_x,
            })
        })
        .collect()
}

pub fn distortion_percent(actual_image_height: f64, ideal_image_height: f64) -> Result<f64> {
    require(
        ideal_image_height != 0.0,
        "ideal image height must be nonzero",
    )?;
    Ok(100.0 * (actual_image_height - ideal_image_height) / ideal_image_height)
}

fn dft(x: &[C]) -> Vec<C> {
    let n = x.len();
    let mut out = vec![C::new(0.0, 0.0); n];
    for (k, o) in out.iter_mut().enumerate() {
        for (j, v) in x.iter().enumerate() {
            let a = -2.0 * PI * ((k * j) as f64) / n as f64;
            *o += v * C::new(a.cos(), a.sin());
        }
    }
    out
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DiffractionSample {
    pub coordinate: f64,
    pub intensity: f64,
}

/// Far-field intensity of a sampled 1-D pupil, centred and peak-normalised, at screen
/// coordinate `lambda z f`.
pub fn fraunhofer_1d(
    pupil: &[C],
    sample_pitch: f64,
    wavelength: f64,
    propagation_distance: f64,
) -> Result<Vec<DiffractionSample>> {
    require(
        pupil.len() >= 2 && sample_pitch > 0.0 && wavelength > 0.0 && propagation_distance > 0.0,
        "invalid Fraunhofer inputs",
    )?;
    let f = dft(pupil);
    let n = pupil.len();
    let max = f.iter().fold(0.0_f64, |m, q| m.max(q.norm_sqr()));
    Ok((0..n)
        .map(|k| {
            let centered = k as i64 - (n / 2) as i64;
            let src = (k + n / 2) % n;
            let freq = centered as f64 / (n as f64 * sample_pitch);
            DiffractionSample {
                coordinate: wavelength * propagation_distance * freq,
                intensity: if max > 0.0 {
                    f[src].norm_sqr() / max
                } else {
                    0.0
                },
            }
        })
        .collect())
}

/// `|DFT(pupil)|^2` normalised to unit sum (uncentred).
pub fn normalized_psf_1d(pupil: &[C]) -> Result<Vec<f64>> {
    require(!pupil.is_empty(), "empty pupil")?;
    let mut out: Vec<f64> = dft(pupil).iter().map(|v| v.norm_sqr()).collect();
    let sum: f64 = out.iter().sum();
    if sum > 0.0 {
        out.iter_mut().for_each(|v| *v /= sum);
    }
    Ok(out)
}

/// `|DFT(psf)| / |DFT(psf)[0]|`.
pub fn mtf_from_psf_1d(psf: &[f64]) -> Result<Vec<f64>> {
    require(!psf.is_empty(), "empty PSF")?;
    let x: Vec<C> = psf.iter().map(|v| C::new(*v, 0.0)).collect();
    let f = dft(&x);
    let dc = f[0].norm();
    Ok(f.iter()
        .map(|v| if dc > 0.0 { v.norm() / dc } else { 0.0 })
        .collect())
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WavefrontSample {
    pub rho: f64,
    pub theta: f64,
    pub opd: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ZernikeCoefficient {
    pub n: i32,
    pub m: i32,
    pub value: f64,
}

fn factorial(n: i32) -> f64 {
    let mut v = 1.0;
    for i in 2..=n {
        v *= f64::from(i);
    }
    v
}

fn radial(n: i32, m: i32, r: f64) -> f64 {
    let m = m.abs();
    if n < 0 || m > n || (n - m) & 1 != 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for k in 0..=(n - m) / 2 {
        let c = (-1.0f64).powf(f64::from(k)) * factorial(n - k)
            / (factorial(k) * factorial((n + m) / 2 - k) * factorial((n - m) / 2 - k));
        sum += c * r.powf(f64::from(n - 2 * k));
    }
    sum
}

/// Unnormalised Zernike polynomial: `R_n^m cos(m theta)` for `m > 0`, `R_n^|m| sin(|m| theta)`
/// for `m < 0`.
pub fn zernike(n: i32, m: i32, rho: f64, theta: f64) -> Result<f64> {
    require(
        (0.0..=1.0).contains(&rho) && theta.is_finite(),
        "invalid Zernike coordinate",
    )?;
    let r = radial(n, m, rho);
    Ok(if m > 0 {
        r * (f64::from(m) * theta).cos()
    } else if m < 0 {
        r * (f64::from(-m) * theta).sin()
    } else {
        r
    })
}

fn solve_dense(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Result<Vec<f64>> {
    for c in 0..n {
        let mut p = c;
        for r in c + 1..n {
            if a[r * n + c].abs() > a[p * n + c].abs() {
                p = r;
            }
        }
        if a[p * n + c].abs() < 1e-14 {
            return Err(SolverError::Singular("singular Zernike fit".into()));
        }
        if p != c {
            for j in c..n {
                a.swap(p * n + j, c * n + j);
            }
            b.swap(p, c);
        }
        let d = a[c * n + c];
        for j in c..n {
            a[c * n + j] /= d;
        }
        b[c] /= d;
        for r in 0..n {
            if r != c {
                let f = a[r * n + c];
                for j in c..n {
                    a[r * n + j] -= f * a[c * n + j];
                }
                b[r] -= f * b[c];
            }
        }
    }
    Ok(b)
}

/// Least-squares fit of all modes up to `max_radial_order` (normal equations, Gauss-Jordan).
pub fn fit_zernike(
    samples: &[WavefrontSample],
    max_radial_order: u32,
) -> Result<Vec<ZernikeCoefficient>> {
    require(max_radial_order <= 64, "Zernike radial order too large")?;
    let mut modes = Vec::new();
    for n in 0..=max_radial_order as i32 {
        let mut m = -n;
        while m <= n {
            modes.push((n, m));
            m += 2;
        }
    }
    require(
        samples.len() >= modes.len(),
        "insufficient wavefront samples",
    )?;
    let n = modes.len();
    let (mut a, mut b) = (vec![0.0; n * n], vec![0.0; n]);
    for q in samples {
        let v = modes
            .iter()
            .map(|(mn, mm)| zernike(*mn, *mm, q.rho, q.theta))
            .collect::<Result<Vec<_>>>()?;
        for i in 0..n {
            b[i] += v[i] * q.opd;
            for j in 0..n {
                a[i * n + j] += v[i] * v[j];
            }
        }
    }
    let x = solve_dense(a, b, n)?;
    Ok(modes
        .iter()
        .zip(x)
        .map(|((n, m), value)| ZernikeCoefficient {
            n: *n,
            m: *m,
            value,
        })
        .collect())
}

pub fn wavefront_rms(samples: &[WavefrontSample], remove_piston: bool) -> Result<f64> {
    require(!samples.is_empty(), "empty wavefront")?;
    let n = samples.len() as f64;
    let mean = if remove_piston {
        samples.iter().map(|q| q.opd).sum::<f64>() / n
    } else {
        0.0
    };
    Ok((samples
        .iter()
        .map(|q| (q.opd - mean) * (q.opd - mean))
        .sum::<f64>()
        / n)
        .sqrt())
}

pub fn longitudinal_chromatic_shift(focal_positions: &[f64]) -> Result<f64> {
    require(!focal_positions.is_empty(), "empty chromatic focal set")?;
    let lo = focal_positions
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let hi = focal_positions
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    Ok(hi - lo)
}

/// `n_e(theta) = n_o n_e / sqrt(n_e^2 cos^2 + n_o^2 sin^2)`.
pub fn uniaxial_extraordinary_index(
    ordinary_index: f64,
    extraordinary_index: f64,
    angle_to_optic_axis: f64,
) -> Result<f64> {
    let (no, ne) = (ordinary_index, extraordinary_index);
    require(
        no > 0.0 && ne > 0.0 && angle_to_optic_axis.is_finite(),
        "invalid uniaxial indices",
    )?;
    let (c, s) = (angle_to_optic_axis.cos(), angle_to_optic_axis.sin());
    Ok(no * ne / (ne * ne * c * c + no * no * s * s).sqrt())
}

pub fn birefringent_retardance(
    ordinary_index: f64,
    extraordinary_effective_index: f64,
    thickness: f64,
    vacuum_wavelength: f64,
) -> Result<f64> {
    require(
        ordinary_index > 0.0
            && extraordinary_effective_index > 0.0
            && thickness >= 0.0
            && vacuum_wavelength > 0.0,
        "invalid retardance controls",
    )?;
    Ok(2.0 * PI * (extraordinary_effective_index - ordinary_index) * thickness / vacuum_wavelength)
}

/// Linear retarder: slow-axis component delayed by `retardance` (fast axis at `fast_axis_angle`).
pub fn apply_linear_retarder(
    input: JonesVector,
    retardance: f64,
    fast_axis_angle: f64,
) -> JonesVector {
    let (c, s) = (fast_axis_angle.cos(), fast_axis_angle.sin());
    let fast = input.x * c + input.y * s;
    let slow = (input.x * -s + input.y * c) * C::new(0.0, retardance).exp();
    JonesVector {
        x: fast * c - slow * s,
        y: fast * s + slow * c,
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearOpticalMaterial {
    pub relative_permittivity: f64,
    pub chi3_m2_per_v2: C,
    pub saturation_intensity_w_per_m2: f64,
}

impl Default for NonlinearOpticalMaterial {
    fn default() -> Self {
        Self {
            relative_permittivity: 1.0,
            chi3_m2_per_v2: C::new(0.0, 0.0),
            saturation_intensity_w_per_m2: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearPolarizationResponse {
    pub intensity_w_per_m2: f64,
    pub linear_polarization_c_per_m2: C,
    pub nonlinear_polarization_c_per_m2: C,
    pub total_polarization_c_per_m2: C,
    pub effective_relative_permittivity: C,
}

/// Upstream default wave impedance, ohm.
pub const FREE_SPACE_WAVE_IMPEDANCE: f64 = 376.730313668;

/// Scalar isotropic Kerr response with optional saturation: a constitutive hook, not a
/// propagator.
pub fn evaluate_nonlinear_polarization(
    electric_field_v_per_m: C,
    material: &NonlinearOpticalMaterial,
    wave_impedance_ohm: f64,
) -> Result<NonlinearPolarizationResponse> {
    const EPS0: f64 = 8.8541878128e-12;
    let e = electric_field_v_per_m;
    require(
        e.re.is_finite() && e.im.is_finite(),
        "nonlinear optical electric field must be finite",
    )?;
    require(
        material.chi3_m2_per_v2.re.is_finite() && material.chi3_m2_per_v2.im.is_finite(),
        "nonlinear optical chi3 must be finite",
    )?;
    require(
        material.relative_permittivity > 0.0
            && material.relative_permittivity.is_finite()
            && wave_impedance_ohm > 0.0
            && wave_impedance_ohm.is_finite()
            && material.saturation_intensity_w_per_m2 >= 0.0
            && material.saturation_intensity_w_per_m2.is_finite(),
        "invalid nonlinear optical material controls",
    )?;
    let field_squared = e.norm_sqr();
    let intensity = 0.5 * field_squared / wave_impedance_ohm;
    let saturation = if material.saturation_intensity_w_per_m2 > 0.0 {
        1.0 / (1.0 + intensity / material.saturation_intensity_w_per_m2)
    } else {
        1.0
    };
    let linear = e * (EPS0 * (material.relative_permittivity - 1.0));
    let nonlinear = material.chi3_m2_per_v2 * EPS0 * field_squared * saturation * e;
    let total = linear + nonlinear;
    let effective = if e.norm() > f64::EPSILON {
        C::new(1.0, 0.0) + total / (e * EPS0)
    } else {
        C::new(material.relative_permittivity, 0.0)
    };
    Ok(NonlinearPolarizationResponse {
        intensity_w_per_m2: intensity,
        linear_polarization_c_per_m2: linear,
        nonlinear_polarization_c_per_m2: nonlinear,
        total_polarization_c_per_m2: total,
        effective_relative_permittivity: effective,
    })
}

/// Rotation `Rz Ry Rx` plus translation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidTransform {
    t: Vec3,
    r: [f64; 9],
}

impl Default for RigidTransform {
    fn default() -> Self {
        Self {
            t: Vec3::ZERO,
            r: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        }
    }
}

impl RigidTransform {
    pub fn from_euler_xyz(translation: Vec3, rx: f64, ry: f64, rz: f64) -> Self {
        let (cx, sx, cy, sy, cz, sz) = (rx.cos(), rx.sin(), ry.cos(), ry.sin(), rz.cos(), rz.sin());
        Self {
            t: translation,
            r: [
                cz * cy,
                cz * sy * sx - sz * cx,
                cz * sy * cx + sz * sx,
                sz * cy,
                sz * sy * sx + cz * cx,
                sz * sy * cx - cz * sx,
                -sy,
                cy * sx,
                cy * cx,
            ],
        }
    }
    pub fn apply_vector(&self, v: Vec3) -> Vec3 {
        let m = &self.r;
        Vec3::new(
            m[0] * v.x + m[1] * v.y + m[2] * v.z,
            m[3] * v.x + m[4] * v.y + m[5] * v.z,
            m[6] * v.x + m[7] * v.y + m[8] * v.z,
        )
    }
    pub fn apply_point(&self, p: Vec3) -> Vec3 {
        let v = self.apply_vector(p);
        Vec3::new(v.x + self.t.x, v.y + self.t.y, v.z + self.t.z)
    }
    pub fn inverse_vector(&self, v: Vec3) -> Vec3 {
        let m = &self.r;
        Vec3::new(
            m[0] * v.x + m[3] * v.y + m[6] * v.z,
            m[1] * v.x + m[4] * v.y + m[7] * v.z,
            m[2] * v.x + m[5] * v.y + m[8] * v.z,
        )
    }
    pub fn inverse_point(&self, p: Vec3) -> Vec3 {
        self.inverse_vector(Vec3::new(p.x - self.t.x, p.y - self.t.y, p.z - self.t.z))
    }
    pub fn apply(&self, mut r: Ray) -> Ray {
        r.origin = self.apply_point(r.origin);
        r.direction = self.apply_vector(r.direction);
        r
    }
    pub fn inverse(&self, mut r: Ray) -> Ray {
        r.origin = self.inverse_point(r.origin);
        r.direction = self.inverse_vector(r.direction);
        r
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FieldType {
    #[default]
    Angle,
    ObjectHeight,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldPoint {
    pub x: f64,
    pub y: f64,
    pub weight: f64,
}

impl Default for FieldPoint {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            weight: 1.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct FieldModel {
    field_type: FieldType,
    fields: Vec<FieldPoint>,
}

impl FieldModel {
    pub fn new(field_type: FieldType) -> Self {
        Self {
            field_type,
            fields: Vec::new(),
        }
    }
    pub fn add(&mut self, field: FieldPoint) -> Result<()> {
        require(field.weight >= 0.0, "invalid optical field weight")?;
        self.fields.push(field);
        Ok(())
    }
    pub fn field_type(&self) -> FieldType {
        self.field_type
    }
    pub fn fields(&self) -> &[FieldPoint] {
        &self.fields
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Wavelength {
    pub nm: f64,
    pub weight: f64,
    pub primary: bool,
}

impl Default for Wavelength {
    fn default() -> Self {
        Self {
            nm: 0.0,
            weight: 1.0,
            primary: false,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct WavelengthModel {
    wavelengths: Vec<Wavelength>,
}

impl WavelengthModel {
    /// A new primary wavelength clears the flag on earlier entries.
    pub fn add(&mut self, w: Wavelength) -> Result<()> {
        require(w.nm > 0.0 && w.weight >= 0.0, "invalid optical wavelength")?;
        if w.primary {
            self.wavelengths.iter_mut().for_each(|x| x.primary = false);
        }
        self.wavelengths.push(w);
        Ok(())
    }
    /// The primary wavelength, else the first.
    pub fn primary_nm(&self) -> Result<f64> {
        if let Some(w) = self.wavelengths.iter().find(|w| w.primary) {
            return Ok(w.nm);
        }
        self.wavelengths
            .first()
            .map(|w| w.nm)
            .ok_or_else(|| SolverError::InvalidInput("wavelength model is empty".into()))
    }
    pub fn wavelengths(&self) -> &[Wavelength] {
        &self.wavelengths
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ApertureKind {
    #[default]
    Circular,
    Rectangular,
    Elliptical,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aperture {
    pub kind: ApertureKind,
    pub x_radius: f64,
    pub y_radius: f64,
    pub stop: bool,
}

impl Default for Aperture {
    fn default() -> Self {
        Self {
            kind: ApertureKind::Circular,
            x_radius: 1.0,
            y_radius: 1.0,
            stop: false,
        }
    }
}

impl Aperture {
    pub fn contains(&self, x: f64, y: f64) -> bool {
        if self.x_radius <= 0.0
            || self.y_radius <= 0.0
            || self.x_radius.is_nan()
            || self.y_radius.is_nan()
        {
            return false;
        }
        match self.kind {
            ApertureKind::Circular => x * x + y * y <= self.x_radius * self.x_radius,
            ApertureKind::Rectangular => x.abs() <= self.x_radius && y.abs() <= self.y_radius,
            ApertureKind::Elliptical => {
                x * x / (self.x_radius * self.x_radius) + y * y / (self.y_radius * self.y_radius)
                    <= 1.0
            }
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PolynomialTerm {
    pub x_power: u32,
    pub y_power: u32,
    pub coefficient: f64,
}

fn ipow(x: f64, p: u32) -> f64 {
    let mut v = 1.0;
    for _ in 0..p {
        v *= x;
    }
    v
}

/// `z = vertex_z + sum c x^i y^j`, intersected by Newton iteration along the ray.
#[derive(Clone, Debug, PartialEq)]
pub struct PolynomialFreeformSurface {
    vertex_z: f64,
    terms: Vec<PolynomialTerm>,
    aperture_radius: f64,
}

impl PolynomialFreeformSurface {
    /// Upstream default aperture radius 1e30.
    pub fn new(vertex_z: f64, terms: Vec<PolynomialTerm>, aperture_radius: f64) -> Result<Self> {
        require(
            vertex_z.is_finite() && aperture_radius > 0.0 && aperture_radius.is_finite(),
            "invalid freeform surface",
        )?;
        require(
            terms.iter().all(|t| t.coefficient.is_finite()),
            "invalid freeform coefficient",
        )?;
        Ok(Self {
            vertex_z,
            terms,
            aperture_radius,
        })
    }
    pub fn sag(&self, x: f64, y: f64) -> f64 {
        self.terms
            .iter()
            .map(|t| t.coefficient * ipow(x, t.x_power) * ipow(y, t.y_power))
            .fold(0.0, |a, b| a + b)
    }
    /// Unit normal of `F = sag - z`, oriented towards `-z`.
    pub fn normal(&self, x: f64, y: f64) -> Vec3 {
        let (mut sx, mut sy) = (0.0, 0.0);
        for t in &self.terms {
            if t.x_power > 0 {
                sx += t.coefficient
                    * f64::from(t.x_power)
                    * ipow(x, t.x_power - 1)
                    * ipow(y, t.y_power);
            }
            if t.y_power > 0 {
                sy += t.coefficient
                    * f64::from(t.y_power)
                    * ipow(x, t.x_power)
                    * ipow(y, t.y_power - 1);
            }
        }
        let len = (sx * sx + sy * sy + 1.0).sqrt();
        Vec3::new(sx / len, sy / len, -1.0 / len)
    }
    /// `(point, normal)` of the first forward hit inside the aperture (upstream defaults 40
    /// iterations, tolerance 1e-12).
    pub fn intersect(
        &self,
        ray: &Ray,
        max_iterations: usize,
        tolerance: f64,
    ) -> Option<(Vec3, Vec3)> {
        if max_iterations == 0
            || tolerance.is_nan()
            || tolerance <= 0.0
            || ray.direction.z.abs() < 1e-15
        {
            return None;
        }
        let d = ray.direction;
        let mut t = (self.vertex_z - ray.origin.z) / d.z;
        for _ in 0..max_iterations {
            let p = Vec3::new(
                ray.origin.x + t * d.x,
                ray.origin.y + t * d.y,
                ray.origin.z + t * d.z,
            );
            if p.x * p.x + p.y * p.y > self.aperture_radius * self.aperture_radius {
                return None;
            }
            let n = self.normal(p.x, p.y);
            let f = p.z - self.vertex_z - self.sag(p.x, p.y);
            if f.abs() < tolerance {
                return (t > 0.0).then_some((p, n));
            }
            let denom = -(n.x * d.x + n.y * d.y + n.z * d.z) / (-n.z);
            if denom.abs() < 1e-15 {
                return None;
            }
            t -= f / denom;
            if !t.is_finite() {
                return None;
            }
        }
        None
    }
}
