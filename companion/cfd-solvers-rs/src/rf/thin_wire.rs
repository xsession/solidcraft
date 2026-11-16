//! Thin-wire method of moments: pulse basis, point matching and the dyadic free-space kernel,
//! with endpoint junctions, a PEC ground image, lumped and transmission-line loads, loss
//! accounting and multiport extraction.
#![allow(clippy::needless_range_loop)]
use super::{inverse, nport_z_to_s, ComplexMatrix};
use crate::{require, Complex64 as C, Result, SolverError};
use std::f64::consts::PI;

const EPS0: f64 = 8.8541878128e-12;
const MU0: f64 = 1.25663706212e-6;

fn invalid(msg: &str) -> SolverError {
    SolverError::InvalidInput(msg.into())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinWireMomConfig {
    pub length_m: f64,
    pub radius_m: f64,
    pub frequency_hz: f64,
    /// Odd: the feed is the center segment.
    pub segments: usize,
    pub feed_voltage_v: f64,
    pub quadrature_order: usize,
}

impl Default for ThinWireMomConfig {
    fn default() -> Self {
        Self {
            length_m: 0.0,
            radius_m: 0.0,
            frequency_hz: 0.0,
            segments: 31,
            feed_voltage_v: 1.0,
            quadrature_order: 8,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ThinWireMomResult {
    pub segment_center_z_m: Vec<f64>,
    pub current_a: Vec<C>,
    pub feed_impedance_ohm: C,
    pub feed_current_a: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParallelThinWire {
    pub x_m: f64,
    pub y_m: f64,
    pub center_z_m: f64,
    pub length_m: f64,
    pub radius_m: f64,
    pub segments: usize,
}

impl Default for ParallelThinWire {
    fn default() -> Self {
        Self {
            x_m: 0.0,
            y_m: 0.0,
            center_z_m: 0.0,
            length_m: 0.0,
            radius_m: 0.0,
            segments: 31,
        }
    }
}

/// `segment: None` selects the center segment (upstream `size_t::max()`), which needs odd
/// segmentation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinWireDeltaGapFeed {
    pub wire: usize,
    pub segment: Option<usize>,
    pub voltage_v: C,
}

impl Default for ThinWireDeltaGapFeed {
    fn default() -> Self {
        Self {
            wire: 0,
            segment: None,
            voltage_v: C::new(1.0, 0.0),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThinWireSeriesLoad {
    pub wire: usize,
    pub segment: usize,
    pub impedance_ohm: C,
}

/// Two-port line section between two segment-current unknowns; signs map segment orientation to
/// current entering each port and must be ±1.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinWireTransmissionLineLoad {
    pub first_wire: usize,
    pub first_segment: usize,
    pub second_wire: usize,
    pub second_segment: usize,
    pub characteristic_impedance_ohm: C,
    /// alpha + j beta.
    pub propagation_constant_per_m: C,
    pub length_m: f64,
    pub first_current_sign: i32,
    pub second_current_sign: i32,
}

impl Default for ThinWireTransmissionLineLoad {
    fn default() -> Self {
        Self {
            first_wire: 0,
            first_segment: 0,
            second_wire: 0,
            second_segment: 0,
            characteristic_impedance_ohm: C::new(50.0, 0.0),
            propagation_constant_per_m: C::new(0.0, 0.0),
            length_m: 0.0,
            first_current_sign: 1,
            second_current_sign: 1,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ParallelWireMomConfig {
    pub frequency_hz: f64,
    pub wires: Vec<ParallelThinWire>,
    pub feeds: Vec<ThinWireDeltaGapFeed>,
    pub loads: Vec<ThinWireSeriesLoad>,
    pub quadrature_order: usize,
}

impl Default for ParallelWireMomConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 0.0,
            wires: vec![],
            feeds: vec![],
            loads: vec![],
            quadrature_order: 8,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WirePoint3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl WirePoint3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

/// Zero conductivity is PEC; positive conductivity adds the round-wire surface resistance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedThinWire {
    pub start_m: WirePoint3,
    pub end_m: WirePoint3,
    pub radius_m: f64,
    pub segments: usize,
    pub conductivity_s_per_m: f64,
    pub relative_permeability: f64,
    /// Equivalent distributed dielectric loss (ohm/m).
    pub dielectric_loss_ohm_per_m: f64,
}

impl Default for OrientedThinWire {
    fn default() -> Self {
        Self {
            start_m: WirePoint3::default(),
            end_m: WirePoint3::default(),
            radius_m: 0.0,
            segments: 31,
            conductivity_s_per_m: 0.0,
            relative_permeability: 1.0,
            dielectric_loss_ohm_per_m: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct OrientedWireSegment {
    pub wire: usize,
    pub local_segment: usize,
    pub center_m: WirePoint3,
    pub tangent: WirePoint3,
    pub length_m: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ThinWireGroundModel {
    #[default]
    FreeSpace,
    PecPlane,
}

#[derive(Clone, Debug, PartialEq)]
pub struct OrientedWireMomConfig {
    pub frequency_hz: f64,
    pub wires: Vec<OrientedThinWire>,
    pub feeds: Vec<ThinWireDeltaGapFeed>,
    pub loads: Vec<ThinWireSeriesLoad>,
    pub transmission_line_loads: Vec<ThinWireTransmissionLineLoad>,
    pub quadrature_order: usize,
    /// Zero keeps wires disjoint; positive enables endpoint-only junctions.
    pub junction_tolerance_m: f64,
    pub ground_model: ThinWireGroundModel,
    pub ground_plane_z_m: f64,
}

impl Default for OrientedWireMomConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 0.0,
            wires: vec![],
            feeds: vec![],
            loads: vec![],
            transmission_line_loads: vec![],
            quadrature_order: 8,
            junction_tolerance_m: 0.0,
            ground_model: ThinWireGroundModel::FreeSpace,
            ground_plane_z_m: 0.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OrientedWireMomResult {
    pub segments: Vec<OrientedWireSegment>,
    pub current_a: Vec<C>,
    pub feed_impedance_ohm: Vec<C>,
    pub feed_current_a: Vec<f64>,
    pub accepted_power_w: f64,
    pub conductor_loss_w: f64,
    pub dielectric_loss_w: f64,
    pub lumped_load_loss_w: f64,
    pub transmission_line_loss_w: f64,
    pub radiated_power_w: f64,
    pub radiation_efficiency: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ParallelWireSegment {
    pub wire: usize,
    pub local_segment: usize,
    pub x_m: f64,
    pub y_m: f64,
    pub center_z_m: f64,
    pub length_m: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParallelWireMomResult {
    pub segments: Vec<ParallelWireSegment>,
    pub current_a: Vec<C>,
    pub feed_impedance_ohm: Vec<C>,
    pub feed_current_a: Vec<f64>,
}

pub(crate) struct Quadrature {
    pub x: Vec<f64>,
    pub w: Vec<f64>,
}

fn legendre(n: usize, z: f64) -> (f64, f64) {
    let (mut p0, mut p1) = (1.0, z);
    for k in 2..=n {
        let kd = k as f64;
        let p = ((2.0 * kd - 1.0) * z * p1 - (kd - 1.0) * p0) / kd;
        p0 = p1;
        p1 = p;
    }
    (p0, p1)
}

/// Gauss-Legendre nodes on [-1, 1] by Newton iteration (orders 2..=32).
pub(crate) fn gauss_legendre(n: usize) -> Result<Quadrature> {
    require(
        (2..=32).contains(&n),
        "thin-wire MoM quadrature order must be 2..32",
    )?;
    let mut q = Quadrature {
        x: vec![0.0; n],
        w: vec![0.0; n],
    };
    let nd = n as f64;
    for i in 0..n.div_ceil(2) {
        let mut z = (PI * (i as f64 + 0.75) / (nd + 0.5)).cos();
        for _ in 0..50 {
            let (p0, p1) = legendre(n, z);
            let derivative = nd * (z * p1 - p0) / (z * z - 1.0);
            let previous = z;
            z -= p1 / derivative;
            if (z - previous).abs() < 1e-15 {
                break;
            }
        }
        let (p0, p1) = legendre(n, z);
        let derivative = nd * (z * p1 - p0) / (z * z - 1.0);
        let weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
        q.x[i] = -z;
        q.x[n - 1 - i] = z;
        q.w[i] = weight;
        q.w[n - 1 - i] = weight;
    }
    Ok(q)
}

/// Gauss-Jordan with partial pivoting on a row-major `n x n` system (upstream `solve_dense`).
pub(crate) fn solve_dense(mut a: Vec<C>, mut b: Vec<C>, n: usize, what: &str) -> Result<Vec<C>> {
    require(
        a.len() == n * n && b.len() == n,
        "dense MoM system size mismatch",
    )?;
    for col in 0..n {
        let mut pivot = col;
        let mut best = a[col * n + col].norm();
        for row in col + 1..n {
            let v = a[row * n + col].norm();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best.is_nan() || best < 1e-24 {
            return Err(SolverError::Singular(format!("singular {what} matrix")));
        }
        if pivot != col {
            for j in col..n {
                a.swap(pivot * n + j, col * n + j);
            }
            b.swap(pivot, col);
        }
        let d = a[col * n + col];
        for j in col..n {
            a[col * n + j] /= d;
        }
        b[col] /= d;
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = a[row * n + col];
            if f.norm() == 0.0 {
                continue;
            }
            for j in col..n {
                let v = a[col * n + j];
                a[row * n + j] -= f * v;
            }
            let v = b[col];
            b[row] -= f * v;
        }
    }
    Ok(b)
}

fn add(a: WirePoint3, b: WirePoint3) -> WirePoint3 {
    WirePoint3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}
fn sub(a: WirePoint3, b: WirePoint3) -> WirePoint3 {
    WirePoint3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}
fn mul(a: WirePoint3, s: f64) -> WirePoint3 {
    WirePoint3::new(a.x * s, a.y * s, a.z * s)
}
fn dot(a: WirePoint3, b: WirePoint3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}
fn norm(a: WirePoint3) -> f64 {
    dot(a, a).sqrt()
}
fn finite(a: WirePoint3) -> bool {
    a.x.is_finite() && a.y.is_finite() && a.z.is_finite()
}

fn segment_distance(p1: WirePoint3, q1: WirePoint3, p2: WirePoint3, q2: WirePoint3) -> f64 {
    let (d1, d2, r) = (sub(q1, p1), sub(q2, p2), sub(p1, p2));
    let (a, e, f) = (dot(d1, d1), dot(d2, d2), dot(d2, r));
    let (mut s, mut t) = (0.0, 0.0);
    const TINY: f64 = 1e-30;
    if a <= TINY && e <= TINY {
        return norm(r);
    }
    if a <= TINY {
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot(d1, r);
        if e <= TINY {
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            if denom.abs() > TINY {
                s = ((b * f - c * e) / denom).clamp(0.0, 1.0);
            }
            t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
        }
    }
    norm(sub(add(p1, mul(d1, s)), add(p2, mul(d2, t))))
}

#[derive(Clone, Copy)]
struct JunctionEndpoint {
    wire: usize,
    start: bool,
}

fn endpoint_junctions(config: &OrientedWireMomConfig) -> Result<Vec<Vec<JunctionEndpoint>>> {
    let tol = config.junction_tolerance_m;
    if tol == 0.0 {
        return Ok(vec![]);
    }
    require(
        tol > 0.0 && tol.is_finite(),
        "junction tolerance must be finite and nonnegative",
    )?;
    let mut endpoints = Vec::with_capacity(2 * config.wires.len());
    let mut minimum_length = f64::INFINITY;
    for (w, wire) in config.wires.iter().enumerate() {
        endpoints.push((wire.start_m, w, true));
        endpoints.push((wire.end_m, w, false));
        minimum_length = minimum_length.min(norm(sub(wire.end_m, wire.start_m)));
    }
    require(
        tol < 0.25 * minimum_length,
        "junction tolerance is too large relative to wire length",
    )?;
    let mut parent: Vec<usize> = (0..endpoints.len()).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for a in 0..endpoints.len() {
        for b in a + 1..endpoints.len() {
            if endpoints[a].1 == endpoints[b].1 {
                continue;
            }
            if norm(sub(endpoints[a].0, endpoints[b].0)) <= tol {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    parent[rb] = ra;
                }
            }
        }
    }
    let mut groups: Vec<Vec<JunctionEndpoint>> = vec![];
    let mut roots: Vec<usize> = vec![];
    for (i, e) in endpoints.iter().enumerate() {
        let root = find(&mut parent, i);
        let g = match roots.iter().position(|r| *r == root) {
            Some(g) => g,
            None => {
                roots.push(root);
                groups.push(vec![]);
                groups.len() - 1
            }
        };
        groups[g].push(JunctionEndpoint {
            wire: e.1,
            start: e.2,
        });
    }
    groups.retain(|g| g.len() >= 2);
    Ok(groups)
}

fn endpoint_pair_is_junction(a: &OrientedThinWire, b: &OrientedThinWire, tol: f64) -> bool {
    tol > 0.0
        && (norm(sub(a.start_m, b.start_m)) <= tol
            || norm(sub(a.start_m, b.end_m)) <= tol
            || norm(sub(a.end_m, b.start_m)) <= tol
            || norm(sub(a.end_m, b.end_m)) <= tol)
}

fn endpoint_junction_overlaps(a: &OrientedThinWire, b: &OrientedThinWire, tol: f64) -> bool {
    let candidates = [
        (a.start_m, b.start_m, true, true),
        (a.start_m, b.end_m, true, false),
        (a.end_m, b.start_m, false, true),
        (a.end_m, b.end_m, false, false),
    ];
    for (ap, bp, as_, bs) in candidates {
        if norm(sub(ap, bp)) > tol {
            continue;
        }
        let ad = if as_ {
            sub(a.end_m, a.start_m)
        } else {
            sub(a.start_m, a.end_m)
        };
        let bd = if bs {
            sub(b.end_m, b.start_m)
        } else {
            sub(b.start_m, b.end_m)
        };
        if dot(ad, bd) / (norm(ad) * norm(bd)) > 1.0 - 1e-10 {
            return true;
        }
    }
    false
}

fn kernel_integral_oriented(
    observation: &OrientedWireSegment,
    source: &OrientedWireSegment,
    source_radius: f64,
    same_wire: bool,
    k: f64,
    omega: f64,
    q: &Quadrature,
) -> Result<C> {
    let j = C::new(0.0, 1.0);
    let tangent_dot = dot(observation.tangent, source.tangent);
    let mut sum = C::new(0.0, 0.0);
    for (x, w) in q.x.iter().zip(&q.w) {
        let source_point = add(
            source.center_m,
            mul(source.tangent, 0.5 * source.length_m * x),
        );
        let d = sub(observation.center_m, source_point);
        let smooth2 = if same_wire {
            source_radius * source_radius
        } else {
            0.0
        };
        let r = (dot(d, d) + smooth2).sqrt();
        if !(r > 0.0 && r.is_finite()) {
            return Err(SolverError::NumericalFailure(
                "singular oriented thin-wire interaction".into(),
            ));
        }
        let (inv_r, inv_r2) = (1.0 / r, 1.0 / (r * r));
        let inv_r3 = inv_r2 * inv_r;
        let obs_radial = dot(observation.tangent, d) * inv_r;
        let src_radial = dot(source.tangent, d) * inv_r;
        let phase = (-j * k * r).exp() / (4.0 * PI);
        let radial =
            (-k * k * inv_r + 3.0 * j * k * inv_r2 + 3.0 * inv_r3) * obs_radial * src_radial;
        let isotropic = (-j * k * inv_r2 - inv_r3) * tangent_dot;
        let vector_term = k * k * tangent_dot * inv_r;
        sum += w * phase * (vector_term + radial + isotropic);
    }
    Ok(j / (omega * EPS0) * (0.5 * source.length_m) * sum)
}

fn pec_image_segment(source: &OrientedWireSegment, plane_z: f64) -> OrientedWireSegment {
    let mut image = *source;
    image.center_m.z = 2.0 * plane_z - source.center_m.z;
    image.tangent = WirePoint3::new(-source.tangent.x, -source.tangent.y, source.tangent.z);
    image
}

fn finite_c(z: C) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn line_z_parameters(load: &ThinWireTransmissionLineLoad) -> Result<(C, C)> {
    require(
        load.length_m > 0.0
            && load.length_m.is_finite()
            && finite_c(load.characteristic_impedance_ohm)
            && load.characteristic_impedance_ohm.re > 0.0
            && finite_c(load.propagation_constant_per_m)
            && (load.first_current_sign == 1 || load.first_current_sign == -1)
            && (load.second_current_sign == 1 || load.second_current_sign == -1),
        "invalid thin-wire transmission-line load",
    )?;
    let gl = load.propagation_constant_per_m * load.length_m;
    let (sh, ch) = (gl.sinh(), gl.cosh());
    require(
        sh.norm() >= 1e-14,
        "thin-wire transmission-line section is singular at this electrical length",
    )?;
    Ok((
        load.characteristic_impedance_ohm * ch / sh,
        load.characteristic_impedance_ohm / sh,
    ))
}

fn feed_segment(segments: usize, feed: &ThinWireDeltaGapFeed) -> Result<usize> {
    match feed.segment {
        None => {
            require(
                segments % 2 == 1,
                "center delta-gap feed requires odd wire segmentation",
            )?;
            Ok(segments / 2)
        }
        Some(s) => {
            require(s < segments, "thin-wire feed segment out of range")?;
            Ok(s)
        }
    }
}

pub fn solve_oriented_thin_wires(config: &OrientedWireMomConfig) -> Result<OrientedWireMomResult> {
    require(
        config.frequency_hz > 0.0
            && config.frequency_hz.is_finite()
            && !config.wires.is_empty()
            && !config.feeds.is_empty(),
        "invalid oriented thin-wire MoM configuration",
    )?;
    let mut offsets = vec![0usize; config.wires.len() + 1];
    let mut lengths = vec![0.0; config.wires.len()];
    for (w, wire) in config.wires.iter().enumerate() {
        let length = norm(sub(wire.end_m, wire.start_m));
        lengths[w] = length;
        require(
            finite(wire.start_m)
                && finite(wire.end_m)
                && length > 0.0
                && wire.radius_m > 0.0
                && wire.radius_m.is_finite()
                && wire.radius_m < 0.1 * length
                && wire.segments >= 5
                && wire.conductivity_s_per_m >= 0.0
                && wire.conductivity_s_per_m.is_finite()
                && wire.relative_permeability > 0.0
                && wire.relative_permeability.is_finite()
                && wire.dielectric_loss_ohm_per_m >= 0.0
                && wire.dielectric_loss_ohm_per_m.is_finite(),
            "invalid oriented thin wire",
        )?;
        if config.ground_model == ThinWireGroundModel::PecPlane {
            let minimum_z = wire.start_m.z.min(wire.end_m.z);
            require(minimum_z - config.ground_plane_z_m > wire.radius_m, "PEC-ground thin wires must remain above the ground plane by more than their radius")?;
        }
        let ds = length / wire.segments as f64;
        require(
            wire.radius_m < 0.5 * ds,
            "oriented thin-wire MoM requires radius < half segment length",
        )?;
        offsets[w + 1] = offsets[w]
            .checked_add(wire.segments)
            .ok_or_else(|| invalid("thin-wire segment count overflow"))?;
    }
    for a in 0..config.wires.len() {
        for b in a + 1..config.wires.len() {
            let (wa, wb) = (&config.wires[a], &config.wires[b]);
            let separation = segment_distance(wa.start_m, wa.end_m, wb.start_m, wb.end_m);
            if separation.is_nan() || separation <= wa.radius_m + wb.radius_m {
                let junction = endpoint_pair_is_junction(wa, wb, config.junction_tolerance_m);
                if !junction || endpoint_junction_overlaps(wa, wb, config.junction_tolerance_m) {
                    return Err(invalid("oriented thin-wire axes overlap/cross outside an enabled endpoint junction"));
                }
            }
        }
    }
    let junctions = endpoint_junctions(config)?;
    let n = offsets[config.wires.len()];
    require(
        n <= 20_000,
        "thin-wire MoM system too large for the dense reference solver",
    )?;
    let mut result = OrientedWireMomResult {
        segments: Vec::with_capacity(n),
        ..Default::default()
    };
    for (w, wire) in config.wires.iter().enumerate() {
        let delta = sub(wire.end_m, wire.start_m);
        let ds = lengths[w] / wire.segments as f64;
        let tangent = mul(delta, 1.0 / lengths[w]);
        for s in 0..wire.segments {
            let fraction = (s as f64 + 0.5) / wire.segments as f64;
            result.segments.push(OrientedWireSegment {
                wire: w,
                local_segment: s,
                center_m: add(wire.start_m, mul(delta, fraction)),
                tangent,
                length_m: ds,
            });
        }
    }
    let c = 1.0 / (MU0 * EPS0).sqrt();
    let omega = 2.0 * PI * config.frequency_hz;
    let k = omega / c;
    let quadrature = gauss_legendre(config.quadrature_order)?;
    let mut matrix = vec![C::new(0.0, 0.0); n * n];
    let mut rhs = vec![C::new(0.0, 0.0); n];
    let conductor_r: Vec<f64> = config
        .wires
        .iter()
        .map(|d| {
            if d.conductivity_s_per_m > 0.0 {
                (PI * config.frequency_hz * MU0 * d.relative_permeability / d.conductivity_s_per_m)
                    .sqrt()
                    / (2.0 * PI * d.radius_m)
            } else {
                0.0
            }
        })
        .collect();
    for row in 0..n {
        let obs = result.segments[row];
        for col in 0..n {
            let src = result.segments[col];
            let radius = config.wires[src.wire].radius_m;
            let mut v = kernel_integral_oriented(
                &obs,
                &src,
                radius,
                obs.wire == src.wire,
                k,
                omega,
                &quadrature,
            )?;
            if config.ground_model == ThinWireGroundModel::PecPlane {
                v += kernel_integral_oriented(
                    &obs,
                    &pec_image_segment(&src, config.ground_plane_z_m),
                    radius,
                    false,
                    k,
                    omega,
                    &quadrature,
                )?;
            }
            matrix[row * n + col] = v;
        }
        matrix[row * n + row] += C::new(
            conductor_r[obs.wire] + config.wires[obs.wire].dielectric_loss_ohm_per_m,
            0.0,
        );
    }
    for load in &config.loads {
        require(
            load.wire < config.wires.len()
                && load.segment < config.wires[load.wire].segments
                && finite_c(load.impedance_ohm),
            "invalid thin-wire series load",
        )?;
        let i = offsets[load.wire] + load.segment;
        matrix[i * n + i] += load.impedance_ohm / result.segments[i].length_m;
    }
    let mut line_ports = Vec::with_capacity(config.transmission_line_loads.len());
    for load in &config.transmission_line_loads {
        require(
            load.first_wire < config.wires.len()
                && load.second_wire < config.wires.len()
                && load.first_segment < config.wires[load.first_wire].segments
                && load.second_segment < config.wires[load.second_wire].segments,
            "thin-wire transmission-line load segment out of range",
        )?;
        let (zd, zt) = line_z_parameters(load)?;
        let first = offsets[load.first_wire] + load.first_segment;
        let second = offsets[load.second_wire] + load.second_segment;
        require(
            first != second,
            "thin-wire transmission-line ports must use distinct segments",
        )?;
        let (s1, s2) = (
            load.first_current_sign as f64,
            load.second_current_sign as f64,
        );
        let (l1, l2) = (
            result.segments[first].length_m,
            result.segments[second].length_m,
        );
        matrix[first * n + first] += s1 * s1 * zd / l1;
        matrix[first * n + second] += s1 * s2 * zt / l1;
        matrix[second * n + first] += s2 * s1 * zt / l2;
        matrix[second * n + second] += s2 * s2 * zd / l2;
        line_ports.push((first, second, s1, s2, zd, zt));
    }
    let mut feed_indices = Vec::with_capacity(config.feeds.len());
    for feed in &config.feeds {
        require(
            feed.wire < config.wires.len()
                && feed.voltage_v.norm() != 0.0
                && finite_c(feed.voltage_v),
            "invalid thin-wire delta-gap feed",
        )?;
        let index = offsets[feed.wire] + feed_segment(config.wires[feed.wire].segments, feed)?;
        feed_indices.push(index);
        rhs[index] -= feed.voltage_v / result.segments[index].length_m;
    }
    result.current_a = if junctions.is_empty() {
        solve_dense(matrix, rhs, n, "thin-wire MoM")?
    } else {
        let dim = n + junctions.len();
        let mut a = vec![C::new(0.0, 0.0); dim * dim];
        let mut b = vec![C::new(0.0, 0.0); dim];
        let mut scale: f64 = 0.0;
        for row in 0..n {
            b[row] = rhs[row];
            scale = scale.max(matrix[row * n + row].norm());
            a[row * dim..row * dim + n].copy_from_slice(&matrix[row * n..row * n + n]);
        }
        scale = scale.max(1.0);
        for (jn, group) in junctions.iter().enumerate() {
            let row = n + jn;
            for e in group {
                let segment = if e.start {
                    offsets[e.wire]
                } else {
                    offsets[e.wire + 1] - 1
                };
                let coefficient = C::new(scale * if e.start { 1.0 } else { -1.0 }, 0.0);
                a[row * dim + segment] = coefficient;
                a[segment * dim + row] = coefficient;
            }
        }
        let mut s = solve_dense(a, b, dim, "thin-wire MoM")?;
        s.truncate(n);
        s
    };
    let first_impedance = config.feeds[0].voltage_v / result.current_a[feed_indices[0]];
    if first_impedance.re < 0.0 {
        result.current_a.iter_mut().for_each(|c| *c = -*c);
    }
    for (feed, &i) in config.feeds.iter().zip(&feed_indices) {
        let current = result.current_a[i];
        if current.norm().is_nan() || current.norm() < 1e-30 {
            return Err(SolverError::NumericalFailure(
                "thin-wire MoM zero feed current".into(),
            ));
        }
        result.feed_impedance_ohm.push(feed.voltage_v / current);
        result.feed_current_a.push(current.norm());
        result.accepted_power_w += 0.5 * (feed.voltage_v * current.conj()).re;
    }
    for (seg, current) in result.segments.iter().zip(&result.current_a) {
        let i2 = current.norm_sqr();
        result.conductor_loss_w += 0.5 * conductor_r[seg.wire] * seg.length_m * i2;
        result.dielectric_loss_w +=
            0.5 * config.wires[seg.wire].dielectric_loss_ohm_per_m * seg.length_m * i2;
    }
    for load in &config.loads {
        if load.impedance_ohm.re > 0.0 {
            result.lumped_load_loss_w += 0.5
                * load.impedance_ohm.re
                * result.current_a[offsets[load.wire] + load.segment].norm_sqr();
        }
    }
    for (first, second, s1, s2, zd, zt) in line_ports {
        let (i1, i2) = (s1 * result.current_a[first], s2 * result.current_a[second]);
        let (v1, v2) = (zd * i1 + zt * i2, zt * i1 + zd * i2);
        result.transmission_line_loss_w += (0.5 * (v1 * i1.conj() + v2 * i2.conj()).re).max(0.0);
    }
    let dissipated = result.conductor_loss_w
        + result.dielectric_loss_w
        + result.lumped_load_loss_w
        + result.transmission_line_loss_w;
    result.radiated_power_w = (result.accepted_power_w - dissipated).max(0.0);
    result.radiation_efficiency = if result.accepted_power_w > 0.0 {
        (result.radiated_power_w / result.accepted_power_w).clamp(0.0, 1.0)
    } else {
        0.0
    };
    Ok(result)
}

/// Each port is driven with 1 V alone; the port-current admittance matrix is inverted.
pub fn oriented_wire_port_impedance_matrix(
    config: &OrientedWireMomConfig,
) -> Result<ComplexMatrix> {
    require(
        !config.feeds.is_empty(),
        "wire multiport extraction requires at least one feed",
    )?;
    let ports = config.feeds.len();
    let mut offsets = vec![0usize; config.wires.len() + 1];
    for (w, wire) in config.wires.iter().enumerate() {
        offsets[w + 1] = offsets[w]
            .checked_add(wire.segments)
            .ok_or_else(|| invalid("thin-wire segment count overflow"))?;
    }
    let mut indices = Vec::with_capacity(ports);
    for feed in &config.feeds {
        require(
            feed.wire < config.wires.len(),
            "wire multiport feed wire out of range",
        )?;
        indices.push(offsets[feed.wire] + feed_segment(config.wires[feed.wire].segments, feed)?);
    }
    let mut admittance = ComplexMatrix::new(ports);
    for excited in 0..ports {
        let mut single = config.clone();
        let mut feed = config.feeds[excited];
        feed.voltage_v = C::new(1.0, 0.0);
        single.feeds = vec![feed];
        let solved = solve_oriented_thin_wires(&single)?;
        for (observed, &i) in indices.iter().enumerate() {
            admittance[(observed, excited)] = solved.current_a[i];
        }
    }
    inverse(&admittance)
}

pub fn oriented_wire_port_s_parameters(
    config: &OrientedWireMomConfig,
    reference_impedance: &[f64],
) -> Result<ComplexMatrix> {
    nport_z_to_s(
        &oriented_wire_port_impedance_matrix(config)?,
        reference_impedance,
    )
}

/// z-directed compatibility wrapper over [`solve_oriented_thin_wires`].
pub fn solve_parallel_thin_wires(config: &ParallelWireMomConfig) -> Result<ParallelWireMomResult> {
    let oriented = OrientedWireMomConfig {
        frequency_hz: config.frequency_hz,
        feeds: config.feeds.clone(),
        loads: config.loads.clone(),
        quadrature_order: config.quadrature_order,
        wires: config
            .wires
            .iter()
            .map(|w| OrientedThinWire {
                start_m: WirePoint3::new(w.x_m, w.y_m, w.center_z_m - 0.5 * w.length_m),
                end_m: WirePoint3::new(w.x_m, w.y_m, w.center_z_m + 0.5 * w.length_m),
                radius_m: w.radius_m,
                segments: w.segments,
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    };
    let solved = solve_oriented_thin_wires(&oriented)?;
    Ok(ParallelWireMomResult {
        segments: solved
            .segments
            .iter()
            .map(|s| ParallelWireSegment {
                wire: s.wire,
                local_segment: s.local_segment,
                x_m: s.center_m.x,
                y_m: s.center_m.y,
                center_z_m: s.center_m.z,
                length_m: s.length_m,
            })
            .collect(),
        current_a: solved.current_a,
        feed_impedance_ohm: solved.feed_impedance_ohm,
        feed_current_a: solved.feed_current_a,
    })
}

pub fn solve_center_fed_thin_wire(config: &ThinWireMomConfig) -> Result<ThinWireMomResult> {
    require(
        config.length_m > 0.0
            && config.radius_m > 0.0
            && config.frequency_hz > 0.0
            && config.length_m.is_finite()
            && config.radius_m.is_finite()
            && config.frequency_hz.is_finite()
            && config.radius_m < 0.1 * config.length_m
            && config.segments >= 5
            && config.segments % 2 == 1
            && config.feed_voltage_v != 0.0
            && config.feed_voltage_v.is_finite(),
        "invalid thin-wire MoM configuration",
    )?;
    let multi = ParallelWireMomConfig {
        frequency_hz: config.frequency_hz,
        quadrature_order: config.quadrature_order,
        wires: vec![ParallelThinWire {
            length_m: config.length_m,
            radius_m: config.radius_m,
            segments: config.segments,
            ..Default::default()
        }],
        feeds: vec![ThinWireDeltaGapFeed {
            wire: 0,
            segment: None,
            voltage_v: C::new(config.feed_voltage_v, 0.0),
        }],
        loads: vec![],
    };
    let solved = solve_parallel_thin_wires(&multi)?;
    Ok(ThinWireMomResult {
        segment_center_z_m: solved.segments.iter().map(|s| s.center_z_m).collect(),
        current_a: solved.current_a,
        feed_impedance_ohm: solved.feed_impedance_ohm[0],
        feed_current_a: solved.feed_current_a[0],
    })
}
