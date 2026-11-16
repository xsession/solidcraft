//! Dense RWG EFIE reference for PEC triangle surfaces, a scalar dielectric SIE reference and a
//! monopole Helmholtz treecode.
#![allow(clippy::needless_range_loop)]
use super::{inverse, ComplexMatrix};
use crate::{core::Vec3, io::TriangleSurface, require, Complex64 as C, Result, SolverError};
use std::collections::BTreeMap;
use std::f64::consts::PI;

const EPS0: f64 = 8.8541878128e-12;
const MU0: f64 = 1.25663706212e-6;

type Point = Vec3;

fn p(x: f64, y: f64, z: f64) -> Point {
    Point { x, y, z }
}
fn add(a: Point, b: Point) -> Point {
    p(a.x + b.x, a.y + b.y, a.z + b.z)
}
fn sub(a: Point, b: Point) -> Point {
    p(a.x - b.x, a.y - b.y, a.z - b.z)
}
fn scale(a: Point, s: f64) -> Point {
    p(a.x * s, a.y * s, a.z * s)
}
fn dot(a: Point, b: Point) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}
fn cross(a: Point, b: Point) -> Point {
    p(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}
fn norm(a: Point) -> f64 {
    dot(a, a).sqrt()
}
fn finite(a: Point) -> bool {
    a.x.is_finite() && a.y.is_finite() && a.z.is_finite()
}
fn normalize(a: Point) -> Result<Point> {
    let l = norm(a);
    require(
        l > 0.0 && l.is_finite(),
        "surface MoM vector must be finite and nonzero",
    )?;
    Ok(scale(a, 1.0 / l))
}

/// Vertex order defines the local normal.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SurfaceTriangle {
    pub a: Point,
    pub b: Point,
    pub c: Point,
}

/// `minus_triangle: None` marks a boundary half-RWG function.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RwgBasisFunction {
    pub edge: usize,
    pub edge_start: Point,
    pub edge_end: Point,
    pub plus_free_vertex: Point,
    pub minus_free_vertex: Point,
    pub plus_triangle: usize,
    pub minus_triangle: Option<usize>,
    pub edge_length_m: f64,
    pub plus_sign: i32,
    pub minus_sign: i32,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SurfaceMomExcitation {
    #[default]
    PlaneWave,
    DeltaGap,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceMomConfig {
    pub frequency_hz: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub excitation: SurfaceMomExcitation,
    pub plane_wave_direction: Point,
    pub plane_wave_electric: Point,
    pub feed_voltage_v: C,
    pub feed_edge: Option<usize>,
    /// Three-point triangle rule with a finite self-interaction distance (not singular quadrature).
    pub singular_regularization_fraction: f64,
    pub mesh_tolerance_m: f64,
}

impl Default for SurfaceMomConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 0.0,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            excitation: SurfaceMomExcitation::PlaneWave,
            plane_wave_direction: p(0.0, 0.0, 1.0),
            plane_wave_electric: p(1.0, 0.0, 0.0),
            feed_voltage_v: C::new(1.0, 0.0),
            feed_edge: None,
            singular_regularization_fraction: 1e-3,
            mesh_tolerance_m: 1e-12,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceMomResult {
    pub triangles: Vec<SurfaceTriangle>,
    pub basis: Vec<RwgBasisFunction>,
    pub current_a: Vec<C>,
    pub system_matrix: ComplexMatrix,
    pub feed_edge: Option<usize>,
    pub feed_current_a: C,
    pub feed_impedance_ohm: C,
    pub accepted_power_w: f64,
    pub relative_residual: f64,
    pub wave_number_per_m: f64,
    pub wave_impedance_ohm: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SurfaceMomFarField {
    pub direction: Point,
    pub e_x: C,
    pub e_y: C,
    pub e_z: C,
    pub power_density_w_m2: f64,
}

#[derive(Clone, Copy)]
struct TriangleGeometry {
    a: Point,
    b: Point,
    c: Point,
    area: f64,
}

fn geometry(t: &SurfaceTriangle, what: &str) -> Result<TriangleGeometry> {
    require(
        finite(t.a) && finite(t.b) && finite(t.c),
        "surface MoM triangle contains a non-finite vertex",
    )?;
    let twice_area = norm(cross(sub(t.b, t.a), sub(t.c, t.a)));
    require(twice_area > 1e-24 && twice_area.is_finite(), what)?;
    Ok(TriangleGeometry {
        a: t.a,
        b: t.b,
        c: t.c,
        area: 0.5 * twice_area,
    })
}

fn samples(t: &TriangleGeometry) -> [(Point, f64); 3] {
    let w = t.area / 3.0;
    [
        (
            add(
                add(scale(t.a, 2.0 / 3.0), scale(t.b, 1.0 / 6.0)),
                scale(t.c, 1.0 / 6.0),
            ),
            w,
        ),
        (
            add(
                add(scale(t.a, 1.0 / 6.0), scale(t.b, 2.0 / 3.0)),
                scale(t.c, 1.0 / 6.0),
            ),
            w,
        ),
        (
            add(
                add(scale(t.a, 1.0 / 6.0), scale(t.b, 1.0 / 6.0)),
                scale(t.c, 2.0 / 3.0),
            ),
            w,
        ),
    ]
}

#[derive(Clone, Copy, Default)]
struct EdgeUse {
    triangle: usize,
    free_point: Point,
    sign: i32,
}

struct BasisBuild {
    public: RwgBasisFunction,
    plus: EdgeUse,
    minus: Option<EdgeUse>,
}

impl BasisBuild {
    fn support(&self, triangle: usize) -> Option<&EdgeUse> {
        if self.plus.triangle == triangle {
            Some(&self.plus)
        } else {
            self.minus.as_ref().filter(|m| m.triangle == triangle)
        }
    }
    fn value(&self, t: &TriangleGeometry, triangle: usize, position: Point) -> Point {
        match self.support(triangle) {
            Some(u) => scale(
                sub(position, u.free_point),
                u.sign as f64 * self.public.edge_length_m / (2.0 * t.area),
            ),
            None => Point::default(),
        }
    }
    fn divergence(&self, t: &TriangleGeometry, triangle: usize) -> f64 {
        self.support(triangle)
            .map_or(0.0, |u| u.sign as f64 * self.public.edge_length_m / t.area)
    }
    fn supports(&self) -> Vec<usize> {
        let mut s = vec![self.plus.triangle];
        if let Some(m) = self.minus {
            s.push(m.triangle);
        }
        s
    }
}

fn find_or_add_vertex(vertices: &mut Vec<Point>, point: Point, tolerance: f64) -> usize {
    let t2 = tolerance * tolerance;
    if let Some(i) = vertices
        .iter()
        .position(|v| dot(sub(*v, point), sub(*v, point)) <= t2)
    {
        return i;
    }
    vertices.push(point);
    vertices.len() - 1
}

pub fn solve_pec_surface_mom(
    input: &[SurfaceTriangle],
    config: &SurfaceMomConfig,
) -> Result<SurfaceMomResult> {
    require(
        !input.is_empty(),
        "surface MoM requires at least one triangle",
    )?;
    let c = config;
    require(
        c.frequency_hz > 0.0
            && c.frequency_hz.is_finite()
            && c.relative_permittivity > 0.0
            && c.relative_permittivity.is_finite()
            && c.relative_permeability > 0.0
            && c.relative_permeability.is_finite()
            && c.singular_regularization_fraction > 0.0
            && c.singular_regularization_fraction.is_finite()
            && c.mesh_tolerance_m > 0.0
            && c.mesh_tolerance_m.is_finite(),
        "invalid surface MoM configuration",
    )?;
    require(
        input.len() <= 1 << 14,
        "surface MoM mesh too large for the dense reference solver",
    )?;
    let direction = normalize(c.plane_wave_direction)?;
    require(
        finite(c.plane_wave_electric),
        "surface MoM plane-wave electric field is not finite",
    )?;
    let mut electric = sub(
        c.plane_wave_electric,
        scale(direction, dot(direction, c.plane_wave_electric)),
    );
    let en = norm(electric);
    require(
        en > 1e-14,
        "plane-wave electric field must have a transverse component",
    )?;
    electric = scale(electric, 1.0 / en);
    let delta_gap = c.excitation == SurfaceMomExcitation::DeltaGap;
    if delta_gap {
        require(
            c.feed_edge.is_some()
                && c.feed_voltage_v.norm() != 0.0
                && c.feed_voltage_v.re.is_finite()
                && c.feed_voltage_v.im.is_finite(),
            "delta-gap surface MoM excitation requires a finite nonzero feed",
        )?;
    }
    let mut vertices = Vec::with_capacity(input.len() * 3);
    let mut tri_vertices = Vec::with_capacity(input.len());
    let mut triangles = Vec::with_capacity(input.len());
    for t in input {
        let a = find_or_add_vertex(&mut vertices, t.a, c.mesh_tolerance_m);
        let b = find_or_add_vertex(&mut vertices, t.b, c.mesh_tolerance_m);
        let cc = find_or_add_vertex(&mut vertices, t.c, c.mesh_tolerance_m);
        require(
            a != b && b != cc && a != cc,
            "surface MoM triangle has repeated vertices",
        )?;
        tri_vertices.push([a, b, cc]);
        triangles.push(geometry(t, "surface MoM triangle is degenerate")?);
    }
    let mut edge_index: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut edges: Vec<(usize, usize, Vec<EdgeUse>)> = vec![];
    for (ti, nodes) in tri_vertices.iter().enumerate() {
        for local in 0..3 {
            let (s, e) = (nodes[local], nodes[(local + 1) % 3]);
            let key = (s.min(e), s.max(e));
            let idx = *edge_index.entry(key).or_insert_with(|| {
                edges.push((key.0, key.1, vec![]));
                edges.len() - 1
            });
            let uses = &mut edges[idx].2;
            require(uses.len() < 2, "surface MoM mesh has a non-manifold edge")?;
            let free = nodes[(local + 2) % 3];
            uses.push(EdgeUse {
                triangle: ti,
                free_point: vertices[free],
                sign: if s == key.0 { 1 } else { -1 },
            });
        }
    }
    let basis: Vec<BasisBuild> = edges
        .iter()
        .enumerate()
        .map(|(edge, (s, e, uses))| {
            let plus = uses[0];
            let minus = uses.get(1).copied();
            let (start, end) = (vertices[*s], vertices[*e]);
            BasisBuild {
                public: RwgBasisFunction {
                    edge,
                    edge_start: start,
                    edge_end: end,
                    plus_free_vertex: plus.free_point,
                    minus_free_vertex: minus.map_or(Point::default(), |m| m.free_point),
                    plus_triangle: plus.triangle,
                    minus_triangle: minus.map(|m| m.triangle),
                    edge_length_m: norm(sub(end, start)),
                    plus_sign: plus.sign,
                    minus_sign: minus.map_or(0, |m| m.sign),
                },
                plus,
                minus,
            }
        })
        .collect();
    require(
        !basis.is_empty(),
        "surface MoM mesh produced no RWG basis functions",
    )?;
    if delta_gap {
        require(
            c.feed_edge.is_some_and(|f| f < basis.len()),
            "surface MoM feed edge out of range",
        )?;
    }
    let omega = 2.0 * PI * c.frequency_hz;
    let epsilon = EPS0 * c.relative_permittivity;
    let mu = MU0 * c.relative_permeability;
    let k = omega * (mu * epsilon).sqrt();
    let eta = (mu / epsilon).sqrt();
    let j = C::new(0.0, 1.0);
    let n = basis.len();
    let mut matrix = ComplexMatrix::new(n);
    let mut rhs = vec![C::new(0.0, 0.0); n];
    for row in 0..n {
        let rb = &basis[row];
        for col in 0..n {
            let cb = &basis[col];
            let mut value = C::new(0.0, 0.0);
            for rt in rb.supports() {
                let row_tri = &triangles[rt];
                for ct in cb.supports() {
                    let col_tri = &triangles[ct];
                    let reg =
                        c.singular_regularization_fraction * row_tri.area.min(col_tri.area).sqrt();
                    for (op, ow) in samples(row_tri) {
                        let test = rb.value(row_tri, rt, op);
                        let test_div = rb.divergence(row_tri, rt);
                        for (sp, sw) in samples(col_tri) {
                            let trial = cb.value(col_tri, ct, sp);
                            let trial_div = cb.divergence(col_tri, ct);
                            let mut distance = norm(sub(op, sp));
                            if distance < reg {
                                distance = reg;
                            }
                            let green = (-j * k * distance).exp() / (4.0 * PI * distance);
                            let term = -j * k * dot(test, trial) + j / k * test_div * trial_div;
                            value += eta * green * term * ow * sw;
                        }
                    }
                }
            }
            matrix[(row, col)] = value;
        }
        if !delta_gap {
            for t in rb.supports() {
                let tri = &triangles[t];
                for (q, w) in samples(tri) {
                    let test = rb.value(tri, t, q);
                    let phase = (-j * k * dot(direction, q)).exp();
                    rhs[row] -= (test.x * (electric.x * phase)
                        + test.y * (electric.y * phase)
                        + test.z * (electric.z * phase))
                        * w;
                }
            }
        }
    }
    let feed = if delta_gap { c.feed_edge } else { None };
    if let Some(f) = feed {
        rhs[f] = c.feed_voltage_v * basis[f].public.edge_length_m;
    }
    let inv = inverse(&matrix)?;
    let mut current = vec![C::new(0.0, 0.0); n];
    for row in 0..n {
        for col in 0..n {
            current[row] += inv[(row, col)] * rhs[col];
        }
    }
    let (mut residual2, mut rhs2) = (0.0, 0.0);
    for row in 0..n {
        let mut r = C::new(0.0, 0.0);
        for col in 0..n {
            r += matrix[(row, col)] * current[col];
        }
        residual2 += (r - rhs[row]).norm_sqr();
        rhs2 += rhs[row].norm_sqr();
    }
    let mut result = SurfaceMomResult {
        triangles: input.to_vec(),
        basis: basis.iter().map(|b| b.public).collect(),
        current_a: current,
        system_matrix: matrix,
        feed_edge: feed,
        feed_current_a: C::new(0.0, 0.0),
        feed_impedance_ohm: C::new(0.0, 0.0),
        accepted_power_w: 0.0,
        relative_residual: residual2.sqrt() / rhs2.sqrt().max(1e-30),
        wave_number_per_m: k,
        wave_impedance_ohm: eta,
    };
    if let Some(f) = feed {
        result.feed_current_a = result.current_a[f];
        if result.feed_current_a.norm() > 1e-30 {
            result.feed_impedance_ohm = c.feed_voltage_v / result.feed_current_a;
            result.accepted_power_w = 0.5 * (c.feed_voltage_v * result.feed_current_a.conj()).re;
        }
    }
    Ok(result)
}

/// Bridge for the OBJ/ASCII-STL triangle readers.
pub fn solve_pec_surface_mom_from_surface(
    surface: &TriangleSurface,
    config: &SurfaceMomConfig,
) -> Result<SurfaceMomResult> {
    require(
        !surface.triangles.is_empty(),
        "surface MoM triangle surface is empty",
    )?;
    let triangles: Vec<SurfaceTriangle> = surface
        .triangles
        .iter()
        .map(|t| SurfaceTriangle {
            a: t.a,
            b: t.b,
            c: t.c,
        })
        .collect();
    solve_pec_surface_mom(&triangles, config)
}

/// Electric-current radiation integral of a solved surface current; the direction is normalized.
pub fn surface_mom_far_field(
    result: &SurfaceMomResult,
    observation_direction: Point,
    distance_m: f64,
) -> Result<SurfaceMomFarField> {
    require(
        !result.triangles.is_empty()
            && result.basis.len() == result.current_a.len()
            && result.wave_number_per_m > 0.0
            && result.wave_impedance_ohm > 0.0
            && distance_m > 0.0
            && distance_m.is_finite(),
        "invalid surface MoM far-field result or observation point",
    )?;
    let direction = normalize(observation_direction)?;
    let triangles = result
        .triangles
        .iter()
        .map(|t| geometry(t, "surface MoM triangle is degenerate"))
        .collect::<Result<Vec<_>>>()?;
    let k = result.wave_number_per_m;
    let (mut ix, mut iy, mut iz) = (C::new(0.0, 0.0), C::new(0.0, 0.0), C::new(0.0, 0.0));
    for (b, current) in result.basis.iter().zip(&result.current_a) {
        let mut integrate = |ti: usize, free: Point, sign: i32| -> Result<()> {
            let t = triangles.get(ti).ok_or_else(|| {
                SolverError::InvalidInput("RWG support triangle out of range".into())
            })?;
            for (q, w) in samples(t) {
                let v = scale(sub(q, free), sign as f64 * b.edge_length_m / (2.0 * t.area));
                let phase = (C::new(0.0, 1.0) * k * dot(direction, q)).exp();
                let wc = C::new(w, 0.0);
                ix += current * v.x * phase * wc;
                iy += current * v.y * phase * wc;
                iz += current * v.z * phase * wc;
            }
            Ok(())
        };
        integrate(b.plus_triangle, b.plus_free_vertex, b.plus_sign)?;
        if let Some(m) = b.minus_triangle {
            integrate(m, b.minus_free_vertex, b.minus_sign)?;
        }
    }
    // d x (I x d)
    let d = direction;
    let (cx, cy, cz) = (
        iy * d.z - iz * d.y,
        iz * d.x - ix * d.z,
        ix * d.y - iy * d.x,
    );
    let tx = d.y * cz - d.z * cy;
    let ty = d.z * cx - d.x * cz;
    let tz = d.x * cy - d.y * cx;
    let factor =
        C::new(0.0, 1.0) * k * result.wave_impedance_ohm * C::new(0.0, -k * distance_m).exp()
            / (4.0 * PI * distance_m);
    let (ex, ey, ez) = (tx * factor, ty * factor, tz * factor);
    Ok(SurfaceMomFarField {
        direction,
        e_x: ex,
        e_y: ey,
        e_z: ez,
        power_density_w_m2: (ex.norm_sqr() + ey.norm_sqr() + ez.norm_sqr())
            / (2.0 * result.wave_impedance_ohm),
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricSieConfig {
    pub frequency_hz: f64,
    pub relative_permittivity_inside: f64,
    pub relative_permittivity_outside: f64,
    pub singular_regularization_fraction: f64,
}

impl Default for DielectricSieConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 0.0,
            relative_permittivity_inside: 4.0,
            relative_permittivity_outside: 1.0,
            singular_regularization_fraction: 1e-3,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DielectricSieResult {
    pub normal_polarization: Vec<C>,
    pub relative_residual: f64,
    pub converged: bool,
}

fn helmholtz_kernel(target: Point, source: Point, k: f64) -> C {
    let d = norm(sub(target, source));
    if d.is_nan() || d <= 1e-15 {
        return C::new(0.0, 0.0);
    }
    C::new(0.0, -k * d).exp() / (4.0 * PI * d)
}

fn solve_gauss(mut m: Vec<C>, mut rhs: Vec<C>) -> Result<Vec<C>> {
    let n = rhs.len();
    require(m.len() == n * n, "integral-equation matrix size mismatch")?;
    for col in 0..n {
        let mut pivot = col;
        let mut mag = m[col * n + col].norm();
        for row in col + 1..n {
            let cand = m[row * n + col].norm();
            if cand > mag {
                mag = cand;
                pivot = row;
            }
        }
        if mag.is_nan() || mag <= 100.0 * f64::EPSILON {
            return Err(SolverError::Singular(
                "singular dielectric SIE matrix".into(),
            ));
        }
        if pivot != col {
            for jj in col..n {
                m.swap(col * n + jj, pivot * n + jj);
            }
            rhs.swap(col, pivot);
        }
        let diag = m[col * n + col];
        for row in col + 1..n {
            let factor = m[row * n + col] / diag;
            m[row * n + col] = C::new(0.0, 0.0);
            for jj in col + 1..n {
                let v = m[col * n + jj];
                m[row * n + jj] -= factor * v;
            }
            let v = rhs[col];
            rhs[row] -= factor * v;
        }
    }
    let mut x = vec![C::new(0.0, 0.0); n];
    for back in 0..n {
        let row = n - 1 - back;
        let mut v = rhs[row];
        for jj in row + 1..n {
            v -= m[row * n + jj] * x[jj];
        }
        x[row] = v / m[row * n + row];
    }
    Ok(x)
}

/// Scalar normal-polarization dielectric SIE with an explicit finite self regularization (not a
/// vectorial singular-quadrature SIE). Upstream default incident field is (1, 0, 0).
pub fn solve_dielectric_surface_sie(
    triangles: &[SurfaceTriangle],
    config: &DielectricSieConfig,
    incident_electric_field: Point,
) -> Result<DielectricSieResult> {
    let c = config;
    require(
        !triangles.is_empty()
            && c.frequency_hz > 0.0
            && c.relative_permittivity_inside > 0.0
            && c.relative_permittivity_outside > 0.0
            && c.singular_regularization_fraction > 0.0
            && c.frequency_hz.is_finite(),
        "invalid dielectric SIE configuration",
    )?;
    require(
        triangles.len() <= 1 << 14,
        "dielectric SIE mesh too large for the dense reference solver",
    )?;
    let contrast = (c.relative_permittivity_inside - c.relative_permittivity_outside)
        / (c.relative_permittivity_inside + 2.0 * c.relative_permittivity_outside);
    if contrast.abs() <= 1e-14 {
        return Ok(DielectricSieResult {
            normal_polarization: vec![C::new(0.0, 0.0); triangles.len()],
            relative_residual: 0.0,
            converged: true,
        });
    }
    let k = 2.0 * PI * c.frequency_hz * (MU0 * EPS0 * c.relative_permittivity_outside).sqrt();
    let mut data = Vec::with_capacity(triangles.len());
    for t in triangles {
        let ab = sub(t.b, t.a);
        let ac = sub(t.c, t.a);
        let tn = cross(ab, ac);
        let ta = norm(tn);
        require(
            ta > 1e-24 && ta.is_finite(),
            "dielectric SIE triangle is degenerate",
        )?;
        data.push((
            scale(add(add(t.a, t.b), t.c), 1.0 / 3.0),
            scale(tn, 1.0 / ta),
            0.5 * ta,
        ));
    }
    let n = data.len();
    let mut m = vec![C::new(0.0, 0.0); n * n];
    let mut rhs = vec![C::new(0.0, 0.0); n];
    let self_term =
        |row: usize| 1.0 / contrast + c.singular_regularization_fraction / data[row].2.sqrt();
    for row in 0..n {
        rhs[row] = C::new(dot(data[row].1, incident_electric_field), 0.0);
        for col in 0..n {
            m[row * n + col] = if row == col {
                C::new(self_term(row), 0.0)
            } else {
                helmholtz_kernel(data[row].0, data[col].0, k) * data[col].2
            };
        }
    }
    let solution = solve_gauss(m, rhs)?;
    let (mut res2, mut rhs2) = (0.0, 0.0);
    for row in 0..n {
        let mut v = solution[row] * self_term(row);
        for col in 0..n {
            if row != col {
                v += helmholtz_kernel(data[row].0, data[col].0, k) * data[col].2 * solution[col];
            }
        }
        let b = dot(data[row].1, incident_electric_field);
        res2 += (v - b).norm_sqr();
        rhs2 += b * b;
    }
    let residual = (res2 / rhs2.max(1e-30)).sqrt();
    Ok(DielectricSieResult {
        normal_polarization: solution,
        relative_residual: residual,
        converged: residual.is_finite() && residual < 1e-9,
    })
}

pub fn helmholtz_direct_apply(
    targets: &[Point],
    sources: &[Point],
    source_strength: &[C],
    wave_number_per_m: f64,
) -> Result<Vec<C>> {
    require(
        sources.len() == source_strength.len()
            && wave_number_per_m >= 0.0
            && wave_number_per_m.is_finite(),
        "invalid Helmholtz interaction inputs",
    )?;
    Ok(targets
        .iter()
        .map(|t| {
            let mut r = C::new(0.0, 0.0);
            for (s, q) in sources.iter().zip(source_strength) {
                r += q * helmholtz_kernel(*t, *s, wave_number_per_m);
            }
            r
        })
        .collect())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HelmholtzTreecodeConfig {
    pub leaf_size: usize,
    /// Effective depth is additionally capped at 1024 to bound recursion.
    pub maximum_depth: usize,
    pub opening_angle: f64,
}

impl Default for HelmholtzTreecodeConfig {
    fn default() -> Self {
        Self {
            leaf_size: 8,
            maximum_depth: 32,
            opening_angle: 0.5,
        }
    }
}

struct TreeNode {
    center: Point,
    radius: f64,
    indices: Vec<usize>,
    child: Option<[usize; 2]>,
    aggregate: C,
}

fn build_tree(
    tree: &mut Vec<TreeNode>,
    indices: Vec<usize>,
    sources: &[Point],
    strengths: &[C],
    leaf: usize,
    max_depth: usize,
    depth: usize,
) -> usize {
    let mut minimum = sources[indices[0]];
    let mut maximum = minimum;
    let mut aggregate = C::new(0.0, 0.0);
    for &i in &indices {
        let q = sources[i];
        minimum = p(minimum.x.min(q.x), minimum.y.min(q.y), minimum.z.min(q.z));
        maximum = p(maximum.x.max(q.x), maximum.y.max(q.y), maximum.z.max(q.z));
        aggregate += strengths[i];
    }
    let center = scale(add(minimum, maximum), 0.5);
    let radius = indices
        .iter()
        .fold(0.0, |r: f64, &i| r.max(norm(sub(sources[i], center))));
    let location = tree.len();
    let split_needed = indices.len() > leaf && depth < max_depth;
    tree.push(TreeNode {
        center,
        radius,
        indices,
        child: None,
        aggregate,
    });
    if split_needed {
        let e = sub(maximum, minimum);
        let axis = if e.x >= e.y && e.x >= e.z {
            0
        } else if e.y >= e.z {
            1
        } else {
            2
        };
        let coord = |q: Point| [q.x, q.y, q.z][axis];
        let split = coord(center);
        let (left, right): (Vec<usize>, Vec<usize>) = tree[location]
            .indices
            .iter()
            .partition(|&&i| coord(sources[i]) <= split);
        if !left.is_empty() && !right.is_empty() {
            tree[location].indices.clear();
            let a = build_tree(tree, left, sources, strengths, leaf, max_depth, depth + 1);
            let b = build_tree(tree, right, sources, strengths, leaf, max_depth, depth + 1);
            tree[location].child = Some([a, b]);
        }
    }
    location
}

fn evaluate_tree(
    tree: &[TreeNode],
    node: usize,
    target: Point,
    sources: &[Point],
    strengths: &[C],
    k: f64,
    opening: f64,
) -> C {
    let n = &tree[node];
    let distance = norm(sub(target, n.center));
    match n.child {
        Some(_) if distance > 0.0 && n.radius / distance < opening => {
            n.aggregate * helmholtz_kernel(target, n.center, k)
        }
        None => {
            let mut r = C::new(0.0, 0.0);
            for &i in &n.indices {
                r += strengths[i] * helmholtz_kernel(target, sources[i], k);
            }
            r
        }
        Some([a, b]) => {
            evaluate_tree(tree, a, target, sources, strengths, k, opening)
                + evaluate_tree(tree, b, target, sources, strengths, k, opening)
        }
    }
}

/// Deterministic monopole treecode (MLFMM-inspired seam for small studies, not production MLFMM).
pub fn helmholtz_treecode_apply(
    targets: &[Point],
    sources: &[Point],
    source_strength: &[C],
    wave_number_per_m: f64,
    config: &HelmholtzTreecodeConfig,
) -> Result<Vec<C>> {
    require(
        !sources.is_empty()
            && config.leaf_size > 0
            && config.maximum_depth > 0
            && config.opening_angle > 0.0
            && config.opening_angle < 1.0,
        "invalid Helmholtz treecode configuration",
    )?;
    require(
        sources.len() == source_strength.len(),
        "Helmholtz source size mismatch",
    )?;
    let mut tree = Vec::with_capacity(sources.len() * 2);
    build_tree(
        &mut tree,
        (0..sources.len()).collect(),
        sources,
        source_strength,
        config.leaf_size,
        config.maximum_depth.min(1024),
        0,
    );
    Ok(targets
        .iter()
        .map(|t| {
            evaluate_tree(
                &tree,
                0,
                *t,
                sources,
                source_strength,
                wave_number_per_m,
                config.opening_angle,
            )
        })
        .collect())
}
