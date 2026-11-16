//! Partial-element equivalent circuits for slender round filaments (resistance, partial and
//! mutual inductance, coefficients of potential, skin effect, SPICE export), finite-ground
//! Fresnel reflection, proximity crowding, rectangular PEEC subdivision and modal two-port
//! discontinuities, ported from upstream `cfd/rf/{peec, conductor_models, modal_elements}` (MIT,
//! cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::{abcd_to_s, cascade_abcd, Matrix2C};
use crate::{core::Vec3, lbm::format_general, require, Complex64 as C, Result, SolverError};
use std::f64::consts::PI;

const MU0: f64 = 4.0e-7 * PI;
const EPS0: f64 = 8.8541878128e-12;

/// Points are [`Vec3`]; distances use `sqrt(dot)` exactly as upstream.
pub type Point3 = Vec3;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilamentSegment {
    pub start: Vec3,
    pub end: Vec3,
    pub radius_m: f64,
    pub conductivity_s_per_m: f64,
}

impl Default for FilamentSegment {
    fn default() -> Self {
        Self {
            start: Vec3::ZERO,
            end: Vec3::ZERO,
            radius_m: 0.0,
            conductivity_s_per_m: 5.8e7,
        }
    }
}

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}
fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}
fn scale(a: Vec3, s: f64) -> Vec3 {
    Vec3::new(a.x * s, a.y * s, a.z * s)
}
fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}
fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

#[allow(clippy::excessive_precision)]
const X2: [f64; 2] = [-0.57735026918962576451, 0.57735026918962576451];
const W2: [f64; 2] = [1.0, 1.0];
#[allow(clippy::excessive_precision)]
const X4: [f64; 4] = [
    -0.86113631159405257522,
    -0.33998104358485626480,
    0.33998104358485626480,
    0.86113631159405257522,
];
#[allow(clippy::excessive_precision)]
const W4: [f64; 4] = [
    0.34785484513745385737,
    0.65214515486254614263,
    0.65214515486254614263,
    0.34785484513745385737,
];
#[allow(clippy::excessive_precision)]
const X6: [f64; 6] = [
    -0.93246951420315202781,
    -0.66120938646626451366,
    -0.23861918608319690863,
    0.23861918608319690863,
    0.66120938646626451366,
    0.93246951420315202781,
];
#[allow(clippy::excessive_precision)]
const W6: [f64; 6] = [
    0.17132449237917034504,
    0.36076157304813860757,
    0.46791393457269104739,
    0.46791393457269104739,
    0.36076157304813860757,
    0.17132449237917034504,
];

fn gauss_rule(order: usize) -> (&'static [f64], &'static [f64]) {
    if order <= 2 {
        (&X2, &W2)
    } else if order <= 4 {
        (&X4, &W4)
    } else {
        (&X6, &W6)
    }
}

pub fn filament_length(segment: &FilamentSegment) -> f64 {
    norm(sub(segment.end, segment.start))
}

fn validate(s: &FilamentSegment) -> Result<()> {
    let length = filament_length(s);
    require(
        length > 0.0
            && s.radius_m > 0.0
            && s.conductivity_s_per_m > 0.0
            && s.radius_m.is_finite()
            && s.conductivity_s_per_m.is_finite()
            && length.is_finite(),
        "invalid PEEC filament segment",
    )?;
    require(
        s.radius_m < 0.25 * length,
        "PEEC filament requires a slender conductor",
    )
}

pub fn dc_resistance(segment: &FilamentSegment) -> Result<f64> {
    validate(segment)?;
    Ok(filament_length(segment)
        / (segment.conductivity_s_per_m * PI * segment.radius_m * segment.radius_m))
}

/// `delta = sqrt(2 / (omega mu sigma))` (upstream default `relative_permeability = 1`).
pub fn skin_depth_m(
    conductivity_s_per_m: f64,
    frequency_hz: f64,
    relative_permeability: f64,
) -> Result<f64> {
    require(
        conductivity_s_per_m > 0.0
            && frequency_hz > 0.0
            && relative_permeability > 0.0
            && conductivity_s_per_m.is_finite()
            && frequency_hz.is_finite()
            && relative_permeability.is_finite(),
        "invalid skin-depth parameters",
    )?;
    let omega = 2.0 * PI * frequency_hz;
    Ok((2.0 / (omega * MU0 * relative_permeability * conductivity_s_per_m)).sqrt())
}

/// Annular skin-depth area model; DC resistance at zero frequency.
pub fn round_wire_ac_resistance(
    segment: &FilamentSegment,
    frequency_hz: f64,
    relative_permeability: f64,
) -> Result<f64> {
    validate(segment)?;
    if frequency_hz == 0.0 {
        return dc_resistance(segment);
    }
    let delta = skin_depth_m(
        segment.conductivity_s_per_m,
        frequency_hz,
        relative_permeability,
    )?;
    let inner = (segment.radius_m - delta).max(0.0);
    let area = PI * (segment.radius_m * segment.radius_m - inner * inner);
    Ok(filament_length(segment) / (segment.conductivity_s_per_m * area))
}

/// External partial self-inductance `mu0 l / (2 pi) (ln(2 l / r) - 1)`.
pub fn self_partial_inductance(segment: &FilamentSegment) -> Result<f64> {
    validate(segment)?;
    let l = filament_length(segment);
    Ok(MU0 * l / (2.0 * PI) * ((2.0 * l / segment.radius_m).ln() - 1.0))
}

/// Neumann double integral by tensor Gauss quadrature (order 2, 4 or 6; upstream default 6).
pub fn mutual_partial_inductance(
    first: &FilamentSegment,
    second: &FilamentSegment,
    quadrature_order: usize,
) -> Result<f64> {
    validate(first)?;
    validate(second)?;
    let (d1, d2) = (sub(first.end, first.start), sub(second.end, second.start));
    let tangent_dot = dot(d1, d2);
    if tangent_dot.abs() < 1e-30 {
        return Ok(0.0);
    }
    let (x, w) = gauss_rule(quadrature_order);
    let mut integral = 0.0;
    for i in 0..x.len() {
        let p1 = add(first.start, scale(d1, 0.5 * (x[i] + 1.0)));
        for j in 0..x.len() {
            let p2 = add(second.start, scale(d2, 0.5 * (x[j] + 1.0)));
            let distance = norm(sub(p1, p2));
            require(
                distance >= 0.25 * first.radius_m.min(second.radius_m),
                "PEEC mutual integral has overlapping filaments",
            )?;
            integral += 0.25 * w[i] * w[j] / distance;
        }
    }
    Ok(MU0 / (4.0 * PI) * tangent_dot * integral)
}

fn potential_coefficient(
    first: &FilamentSegment,
    second: &FilamentSegment,
    relative_permittivity: f64,
    order: usize,
    is_self: bool,
) -> Result<f64> {
    let (d1, d2) = (sub(first.end, first.start), sub(second.end, second.start));
    let (x, w) = gauss_rule(order);
    let mut kernel = 0.0;
    for i in 0..x.len() {
        let p1 = add(first.start, scale(d1, 0.5 * (x[i] + 1.0)));
        for j in 0..x.len() {
            let p2 = add(second.start, scale(d2, 0.5 * (x[j] + 1.0)));
            let delta = sub(p1, p2);
            let smoothing = if is_self {
                first.radius_m * first.radius_m
            } else {
                0.0
            };
            let distance = (dot(delta, delta) + smoothing).sqrt();
            if distance.is_nan() || distance <= 0.0 {
                return Err(SolverError::NumericalFailure(
                    "singular PEEC coefficient-of-potential integral".into(),
                ));
            }
            kernel += 0.25 * w[i] * w[j] / distance;
        }
    }
    Ok(kernel / (4.0 * PI * EPS0 * relative_permittivity))
}

fn solve_dense_complex(mut m: Vec<C>, mut rhs: Vec<C>, n: usize) -> Result<Vec<C>> {
    for col in 0..n {
        let mut pivot = col;
        let mut best = m[col * n + col].norm();
        for row in col + 1..n {
            let c = m[row * n + col].norm();
            if c > best {
                best = c;
                pivot = row;
            }
        }
        if best.is_nan() || best < 1e-24 {
            return Err(SolverError::Singular(
                "singular PEEC impedance matrix".into(),
            ));
        }
        if pivot != col {
            for j in col..n {
                m.swap(col * n + j, pivot * n + j);
            }
            rhs.swap(col, pivot);
        }
        let d = m[col * n + col];
        for j in col..n {
            m[col * n + j] /= d;
        }
        rhs[col] /= d;
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = m[row * n + col];
            for j in col..n {
                let v = m[col * n + j];
                m[row * n + j] -= f * v;
            }
            let r = rhs[col];
            rhs[row] -= f * r;
        }
    }
    Ok(rhs)
}

fn inverse_dense_real(mut m: Vec<f64>, n: usize) -> Result<Vec<f64>> {
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for col in 0..n {
        let mut pivot = col;
        let mut best = m[col * n + col].abs();
        for row in col + 1..n {
            let c = m[row * n + col].abs();
            if c > best {
                best = c;
                pivot = row;
            }
        }
        if best.is_nan() || best < 1e-30 {
            return Err(SolverError::Singular(
                "singular PEEC potential matrix".into(),
            ));
        }
        if pivot != col {
            for j in 0..n {
                m.swap(col * n + j, pivot * n + j);
                inv.swap(col * n + j, pivot * n + j);
            }
        }
        let d = m[col * n + col];
        for j in 0..n {
            m[col * n + j] /= d;
            inv[col * n + j] /= d;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = m[row * n + col];
            if f == 0.0 {
                continue;
            }
            for j in 0..n {
                m[row * n + j] -= f * m[col * n + j];
                inv[row * n + j] -= f * inv[col * n + j];
            }
        }
    }
    Ok(inv)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeecMatrices {
    pub size: usize,
    /// Diagonal branch resistances.
    pub resistance_ohm: Vec<f64>,
    /// Row-major partial-inductance matrix.
    pub partial_inductance_h: Vec<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeecCapacitanceMatrices {
    pub size: usize,
    /// Row-major `P` with `V = P q`.
    pub coefficient_of_potential_v_per_c: Vec<f64>,
    /// Row-major Maxwell capacitance `C = P^-1`.
    pub capacitance_f: Vec<f64>,
}

/// Quasi-static filament system (no radiation).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeecFilamentSystem {
    segments: Vec<FilamentSegment>,
}

impl PeecFilamentSystem {
    pub fn add_segment(&mut self, segment: FilamentSegment) -> Result<()> {
        validate(&segment)?;
        self.segments.push(segment);
        Ok(())
    }
    pub fn segments(&self) -> &[FilamentSegment] {
        &self.segments
    }
    /// Upstream default quadrature order 6.
    pub fn extract(&self, quadrature_order: usize) -> Result<PeecMatrices> {
        require(
            !self.segments.is_empty(),
            "PEEC extraction requires at least one segment",
        )?;
        let n = self.segments.len();
        let mut out = PeecMatrices {
            size: n,
            resistance_ohm: vec![0.0; n],
            partial_inductance_h: vec![0.0; n * n],
        };
        for i in 0..n {
            out.resistance_ohm[i] = dc_resistance(&self.segments[i])?;
            out.partial_inductance_h[i * n + i] = self_partial_inductance(&self.segments[i])?;
            for j in 0..i {
                let m = mutual_partial_inductance(
                    &self.segments[i],
                    &self.segments[j],
                    quadrature_order,
                )?;
                out.partial_inductance_h[i * n + j] = m;
                out.partial_inductance_h[j * n + i] = m;
            }
        }
        Ok(out)
    }
    /// Self terms use the radius as a smoothing length (upstream defaults 1, 6).
    pub fn extract_capacitance(
        &self,
        relative_permittivity: f64,
        quadrature_order: usize,
    ) -> Result<PeecCapacitanceMatrices> {
        require(
            !self.segments.is_empty(),
            "PEEC capacitance extraction requires at least one segment",
        )?;
        require(
            relative_permittivity > 0.0 && relative_permittivity.is_finite(),
            "PEEC relative permittivity must be finite and positive",
        )?;
        let n = self.segments.len();
        let mut p = vec![0.0; n * n];
        for i in 0..n {
            p[i * n + i] = potential_coefficient(
                &self.segments[i],
                &self.segments[i],
                relative_permittivity,
                quadrature_order,
                true,
            )?;
            for j in 0..i {
                let v = potential_coefficient(
                    &self.segments[i],
                    &self.segments[j],
                    relative_permittivity,
                    quadrature_order,
                    false,
                )?;
                p[i * n + j] = v;
                p[j * n + i] = v;
            }
        }
        let c = inverse_dense_real(p.clone(), n)?;
        Ok(PeecCapacitanceMatrices {
            size: n,
            coefficient_of_potential_v_per_c: p,
            capacitance_f: c,
        })
    }
    fn impedance(&self, frequency_hz: f64, quadrature_order: usize, skin: bool) -> Result<Vec<C>> {
        require(
            frequency_hz >= 0.0 && frequency_hz.is_finite(),
            "invalid PEEC frequency",
        )?;
        let d = self.extract(quadrature_order)?;
        let omega = 2.0 * PI * frequency_hz;
        let n = d.size;
        let mut z = vec![C::new(0.0, 0.0); n * n];
        for r in 0..n {
            for c in 0..n {
                let resistance = if r != c {
                    0.0
                } else if skin {
                    round_wire_ac_resistance(&self.segments[r], frequency_hz, 1.0)?
                } else {
                    d.resistance_ohm[r]
                };
                z[r * n + c] = C::new(resistance, omega * d.partial_inductance_h[r * n + c]);
            }
        }
        Ok(z)
    }
    /// Row-major `R + j omega Lp`.
    pub fn impedance_matrix(&self, frequency_hz: f64, quadrature_order: usize) -> Result<Vec<C>> {
        self.impedance(frequency_hz, quadrature_order, false)
    }
    /// As [`Self::impedance_matrix`] with skin-effect branch resistances.
    pub fn impedance_matrix_skin_effect(
        &self,
        frequency_hz: f64,
        quadrature_order: usize,
    ) -> Result<Vec<C>> {
        self.impedance(frequency_hz, quadrature_order, true)
    }
    /// Coupled-inductor SPICE subcircuit (R + L per filament, K for mutual coupling).
    pub fn to_spice_subcircuit(&self, name: &str, quadrature_order: usize) -> Result<String> {
        require(
            !self.segments.is_empty(),
            "PEEC export requires at least one segment",
        )?;
        require(!name.is_empty(), "PEEC subcircuit name must not be empty")?;
        let d = self.extract(quadrature_order)?;
        let g = |v: f64| format_general(v, 17);
        let n = d.size;
        let mut o = format!(".SUBCKT {name}");
        for i in 0..n {
            o += &format!(" P{} N{}", i + 1, i + 1);
        }
        o.push('\n');
        for i in 0..n {
            o += &format!(
                "R{} P{} X{} {}\nL{} X{} N{} {}\n",
                i + 1,
                i + 1,
                i + 1,
                g(d.resistance_ohm[i]),
                i + 1,
                i + 1,
                i + 1,
                g(d.partial_inductance_h[i * n + i])
            );
        }
        for i in 0..n {
            for j in 0..i {
                let den =
                    (d.partial_inductance_h[i * n + i] * d.partial_inductance_h[j * n + j]).sqrt();
                if den > 0.0 {
                    let k = d.partial_inductance_h[i * n + j] / den;
                    if k.abs() > 1e-15 {
                        o += &format!("K{}_{} L{} L{} {}\n", j + 1, i + 1, j + 1, i + 1, g(k));
                    }
                }
            }
        }
        o += &format!(".ENDS {name}\n");
        Ok(o)
    }
    /// Branch currents from `Z I = V`.
    pub fn solve_currents(
        &self,
        frequency_hz: f64,
        applied_branch_voltage: &[C],
        quadrature_order: usize,
    ) -> Result<Vec<C>> {
        require(
            applied_branch_voltage.len() == self.segments.len(),
            "PEEC branch-voltage size mismatch",
        )?;
        solve_dense_complex(
            self.impedance_matrix(frequency_hz, quadrature_order)?,
            applied_branch_voltage.to_vec(),
            self.segments.len(),
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
#[allow(non_snake_case)]
pub struct GroundReflection {
    pub TE: C,
    pub TM: C,
}

/// Upstream-compatible "finite ground" reflection. Upstream computes the transmitted cosine as
/// `sqrt(cos^2(theta) / eps_r)` instead of Snell's `sqrt(1 - sin^2(theta) / eps_r)`, so its TE
/// coefficient is identically ~0 and its TM coefficient is `(eps_r - 1) / (eps_r + 1)` for every
/// angle. Kept for value parity; use [`fresnel_ground_reflection`] for physical coefficients.
pub fn finite_ground_reflection(
    frequency_hz: f64,
    conductivity_s_per_m: f64,
    relative_permittivity: f64,
    incidence_angle_rad: f64,
) -> Result<GroundReflection> {
    require(
        frequency_hz > 0.0
            && conductivity_s_per_m >= 0.0
            && relative_permittivity > 0.0
            && frequency_hz.is_finite()
            && conductivity_s_per_m.is_finite()
            && relative_permittivity.is_finite()
            && incidence_angle_rad.is_finite()
            && (0.0..0.5 * PI).contains(&incidence_angle_rad),
        "invalid finite-ground controls",
    )?;
    let omega = 2.0 * PI * frequency_hz;
    let eps_r = C::new(
        relative_permittivity,
        -conductivity_s_per_m / (omega * EPS0),
    );
    let root = eps_r.sqrt();
    let (sin_i, cos_i) = (incidence_angle_rad.sin(), incidence_angle_rad.cos());
    let cos_t = (C::new(1.0 - sin_i * sin_i, 0.0) / eps_r).sqrt();
    let te = (C::new(cos_i, 0.0) - root * cos_t) / (C::new(cos_i, 0.0) + root * cos_t);
    let tm = (root * cos_i - cos_t) / (root * cos_i + cos_t);
    Ok(GroundReflection { TE: te, TM: tm })
}

/// Physical Fresnel TE/TM reflection (`E` and `H` conventions as upstream: TE -> -1 and TM -> +1
/// on a perfect conductor) from a lossy half-space with `eps_r - j sigma / (omega eps0)`.
pub fn fresnel_ground_reflection(
    frequency_hz: f64,
    conductivity_s_per_m: f64,
    relative_permittivity: f64,
    incidence_angle_rad: f64,
) -> Result<GroundReflection> {
    finite_ground_reflection(
        frequency_hz,
        conductivity_s_per_m,
        relative_permittivity,
        incidence_angle_rad,
    )?;
    let omega = 2.0 * PI * frequency_hz;
    let eps_r = C::new(
        relative_permittivity,
        -conductivity_s_per_m / (omega * EPS0),
    );
    let (sin_i, cos_i) = (incidence_angle_rad.sin(), incidence_angle_rad.cos());
    let root = (eps_r - sin_i * sin_i).sqrt();
    let te = (C::new(cos_i, 0.0) - root) / (C::new(cos_i, 0.0) + root);
    let tm = (eps_r * cos_i - root) / (eps_r * cos_i + root);
    Ok(GroundReflection { TE: te, TM: tm })
}

/// Bounded engineering crowding factor `1 + (r / gap)^2 min(10, r / delta) / 2`.
pub fn round_wire_proximity_factor(
    center_spacing_m: f64,
    radius_m: f64,
    frequency_hz: f64,
    conductivity_s_per_m: f64,
    relative_permeability: f64,
) -> Result<f64> {
    require(
        center_spacing_m > 2.0 * radius_m
            && radius_m > 0.0
            && frequency_hz > 0.0
            && conductivity_s_per_m > 0.0
            && relative_permeability > 0.0,
        "invalid round-wire proximity controls",
    )?;
    let delta = skin_depth_m(conductivity_s_per_m, frequency_hz, relative_permeability)?;
    let gap = center_spacing_m - 2.0 * radius_m;
    let crowding = (radius_m / delta).min(10.0);
    Ok(1.0 + 0.5 * (radius_m / gap).powf(2.0) * crowding)
}

/// Tiles the y-z cross-section of an x-directed segment into `ny * nz` filaments.
pub fn subdivide_rectangular_peec(
    centerline: &FilamentSegment,
    width_y_m: f64,
    height_z_m: f64,
    divisions_y: usize,
    divisions_z: usize,
) -> Result<Vec<FilamentSegment>> {
    require(
        filament_length(centerline) > 0.0
            && width_y_m > 0.0
            && height_z_m > 0.0
            && divisions_y > 0
            && divisions_z > 0
            && centerline.conductivity_s_per_m > 0.0,
        "invalid rectangular PEEC subdivision",
    )?;
    let (dy, dz) = (
        width_y_m / divisions_y as f64,
        height_z_m / divisions_z as f64,
    );
    let mut out = Vec::with_capacity(divisions_y.saturating_mul(divisions_z).min(1 << 20));
    for iy in 0..divisions_y {
        for iz in 0..divisions_z {
            let y = (iy as f64 + 0.5) * dy - 0.5 * width_y_m;
            let z = (iz as f64 + 0.5) * dz - 0.5 * height_z_m;
            out.push(FilamentSegment {
                start: Vec3::new(
                    centerline.start.x,
                    centerline.start.y + y,
                    centerline.start.z + z,
                ),
                end: Vec3::new(centerline.end.x, centerline.end.y + y, centerline.end.z + z),
                radius_m: dy.min(dz) * 0.5,
                conductivity_s_per_m: centerline.conductivity_s_per_m,
            });
        }
    }
    Ok(out)
}

fn finite(c: C) -> bool {
    c.re.is_finite() && c.im.is_finite()
}

/// ABCD of a series impedance.
pub fn series_modal_discontinuity(impedance_ohm: C) -> Result<Matrix2C> {
    require(finite(impedance_ohm), "invalid modal series impedance")?;
    Ok(Matrix2C {
        a11: C::new(1.0, 0.0),
        a12: impedance_ohm,
        a21: C::new(0.0, 0.0),
        a22: C::new(1.0, 0.0),
    })
}

/// ABCD of a shunt admittance.
pub fn shunt_modal_discontinuity(admittance_siemens: C) -> Result<Matrix2C> {
    require(finite(admittance_siemens), "invalid modal shunt admittance")?;
    Ok(Matrix2C {
        a11: C::new(1.0, 0.0),
        a12: C::new(0.0, 0.0),
        a21: admittance_siemens,
        a22: C::new(1.0, 0.0),
    })
}

/// S parameters of a series element followed by a shunt element (upstream default 50 ohm).
pub fn modal_discontinuity_s(
    series_impedance_ohm: C,
    shunt_admittance_siemens: C,
    reference_impedance: f64,
) -> Result<Matrix2C> {
    require(
        reference_impedance > 0.0 && reference_impedance.is_finite(),
        "invalid modal reference impedance",
    )?;
    abcd_to_s(
        cascade_abcd(
            series_modal_discontinuity(series_impedance_ohm)?,
            shunt_modal_discontinuity(shunt_admittance_siemens)?,
        ),
        reference_impedance,
    )
}
