//! Quantitative post-processing of companion results: colour maps and legends, probes,
//! min/max location, line sampling, region integrals, run-to-run field comparison, goal
//! histories, component tables, metric tables across studies and unit hints.
//!
//! Every function works on already-parsed data ([`crate::vtk::Grid`], CSV or JSON text) and
//! is bounded. Interpolation is linear inside tetrahedra and triangles and uses a
//! tetrahedral split for hexahedra, wedges and pyramids; nothing is extrapolated, and
//! points outside every cell report the nearest cell explicitly.
use crate::vtk::{Association, Field, Grid};
use serde::Serialize;
use serde_json::Value;

/// Samples one line request may return.
pub const MAX_LINE_SAMPLES: usize = 10_000;
/// Goal series and component rows read from one file.
pub const MAX_TABLE_ROWS: usize = 100_000;

#[derive(Debug, thiserror::Error, PartialEq)]
#[error("results: {0}")]
pub struct Error(pub String);

pub(crate) fn bad(message: impl Into<String>) -> Error {
    Error(message.into())
}

// ---------------------------------------------------------------- colour maps

/// Perceptual and diverging maps for field plots. `Thermal` is black-red-yellow-white.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, serde::Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ColorMap {
    #[default]
    Viridis,
    Turbo,
    CoolWarm,
    Thermal,
    Grayscale,
}

impl ColorMap {
    pub const ALL: [ColorMap; 5] = [ColorMap::Viridis, ColorMap::Turbo, ColorMap::CoolWarm, ColorMap::Thermal, ColorMap::Grayscale];

    pub fn name(self) -> &'static str {
        match self {
            ColorMap::Viridis => "Viridis",
            ColorMap::Turbo => "Turbo",
            ColorMap::CoolWarm => "Cool–warm",
            ColorMap::Thermal => "Thermal",
            ColorMap::Grayscale => "Grayscale",
        }
    }

    fn stops(self) -> &'static [[u8; 3]] {
        match self {
            ColorMap::Viridis => &[[68, 1, 84], [65, 68, 135], [42, 120, 142], [34, 168, 132], [122, 209, 81], [253, 231, 37]],
            ColorMap::Turbo => &[
                [48, 18, 59],
                [70, 107, 227],
                [41, 187, 236],
                [49, 241, 152],
                [164, 252, 60],
                [237, 208, 58],
                [251, 128, 34],
                [210, 49, 5],
                [122, 4, 3],
            ],
            ColorMap::CoolWarm => &[[59, 76, 192], [124, 159, 249], [192, 212, 245], [242, 203, 183], [238, 132, 104], [180, 4, 38]],
            ColorMap::Thermal => &[[0, 0, 0], [128, 0, 0], [230, 50, 0], [255, 180, 0], [255, 255, 255]],
            ColorMap::Grayscale => &[[16, 16, 16], [240, 240, 240]],
        }
    }

    /// The colours for values below and above a marked range: ones the map itself never shows,
    /// so a cell out of range cannot be read as a value at an end of the map.
    pub fn out_of_range(self) -> ([u8; 3], [u8; 3]) {
        match self {
            ColorMap::Viridis | ColorMap::Turbo | ColorMap::CoolWarm => ([24, 24, 24], [255, 0, 255]),
            ColorMap::Thermal => ([0, 90, 255], [0, 200, 80]),
            ColorMap::Grayscale => ([0, 90, 255], [230, 40, 40]),
        }
    }

    /// Colour at `t` in `[0, 1]` (clamped; NaN maps to the lowest colour).
    pub fn color(self, t: f64) -> [u8; 3] {
        let stops = self.stops();
        let t = if t.is_nan() { 0.0 } else { t.clamp(0.0, 1.0) };
        let x = t * (stops.len() - 1) as f64;
        let i = (x.floor() as usize).min(stops.len() - 2);
        let f = x - i as f64;
        let (a, b) = (stops[i], stops[i + 1]);
        let mix = |p: u8, q: u8| (f64::from(p) + (f64::from(q) - f64::from(p)) * f).round().clamp(0.0, 255.0) as u8;
        [mix(a[0], b[0]), mix(a[1], b[1]), mix(a[2], b[2])]
    }
}

/// A colour range: automatic from data, or locked by the user so several plots or frames
/// share one scale.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, serde::Deserialize, schemars::JsonSchema)]
pub struct ColorRange {
    pub min: f64,
    pub max: f64,
    pub locked: bool,
}

impl ColorRange {
    pub fn from_values(values: impl IntoIterator<Item = f64>) -> Option<Self> {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for v in values.into_iter().filter(|v| v.is_finite()) {
            lo = lo.min(v);
            hi = hi.max(v);
        }
        (lo <= hi).then_some(Self { min: lo, max: hi, locked: false })
    }

    /// `(value - min) / (max - min)`, 0.5 for a degenerate range.
    pub fn normalize(&self, value: f64) -> f64 {
        let span = self.max - self.min;
        if span.abs() <= f64::EPSILON * self.max.abs().max(self.min.abs()).max(1.0) { 0.5 } else { (value - self.min) / span }
    }
}

/// How a plot's colour range is chosen.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, serde::Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum RangeMode {
    /// The data's lowest to highest value.
    #[default]
    Auto,
    /// The values the user typed, so several plots, frames or runs share one scale.
    Locked,
    /// The 2nd to the 98th percentile: a few extreme cells no longer wash the plot out.
    Percentile,
    /// Centred on zero, out to the largest magnitude: for signed fields on a diverging map.
    Symmetric,
}

impl RangeMode {
    pub const ALL: [RangeMode; 4] = [RangeMode::Auto, RangeMode::Locked, RangeMode::Percentile, RangeMode::Symmetric];

    pub fn key(self) -> &'static str {
        match self {
            RangeMode::Auto => "auto",
            RangeMode::Locked => "locked",
            RangeMode::Percentile => "percentile",
            RangeMode::Symmetric => "symmetric",
        }
    }
}

/// Fraction of values each percentile end leaves out.
pub const PERCENTILE_CLIP: f64 = 0.02;
/// Most values sorted for a percentile range; more are sampled evenly.
pub const MAX_PERCENTILE_VALUES: usize = 1_000_000;

/// The range `mode` gives `values`; `locked` is the user's range for [`RangeMode::Locked`].
/// `None` when there is no finite value (or the locked range is not a range).
pub fn range_for(mode: RangeMode, values: &[f64], locked: Option<(f64, f64)>) -> Option<(f64, f64)> {
    let data = ColorRange::from_values(values.iter().copied())?;
    match mode {
        RangeMode::Auto => Some((data.min, data.max)),
        RangeMode::Locked => locked.filter(|(lo, hi)| lo.is_finite() && hi.is_finite() && lo < hi).or(Some((data.min, data.max))),
        RangeMode::Symmetric => {
            let m = data.min.abs().max(data.max.abs());
            Some((-m, m))
        }
        RangeMode::Percentile => {
            let step = values.len().div_ceil(MAX_PERCENTILE_VALUES).max(1);
            let mut sorted: Vec<f64> = values.iter().copied().step_by(step).filter(|v| v.is_finite()).collect();
            sorted.sort_by(f64::total_cmp);
            let at = |q: f64| sorted.get((sorted.len().saturating_sub(1) as f64 * q).round() as usize).copied();
            Some(match (at(PERCENTILE_CLIP), at(1.0 - PERCENTILE_CLIP)) {
                (Some(lo), Some(hi)) if lo < hi => (lo, hi),
                _ => (data.min, data.max),
            })
        }
    }
}

/// Where a value falls on a scale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Shade {
    /// A colour of the map.
    In([u8; 3]),
    /// Below or above the range, with out-of-range marking on: the map's
    /// [`ColorMap::out_of_range`] colours.
    Below([u8; 3]),
    Above([u8; 3]),
}

/// A colour map as a plot uses it: reversed or not, continuous or in bands, and whether values
/// outside the range take the end colours or colours of their own.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, serde::Deserialize, schemars::JsonSchema)]
pub struct ColorScale {
    pub map: ColorMap,
    #[serde(default)]
    pub reversed: bool,
    /// 0: continuous; otherwise that many bands of one colour each (2 to [`MAX_BANDS`]).
    #[serde(default)]
    pub bands: u8,
    /// Values outside the range are shown as out of range instead of in the end colours.
    #[serde(default)]
    pub mark_out_of_range: bool,
}

/// Most colour bands a scale takes.
pub const MAX_BANDS: u8 = 32;

impl ColorScale {
    /// The colour at `t` in `[0, 1]`, reversed and banded as set.
    pub fn at(&self, t: f64) -> [u8; 3] {
        let mut t = if t.is_nan() { 0.0 } else { t.clamp(0.0, 1.0) };
        if self.bands >= 2 {
            let n = f64::from(self.bands.min(MAX_BANDS));
            // The band's middle, so n bands show n distinct colours from end to end.
            t = (((t * n).floor().min(n - 1.0)) + 0.5) / n;
        }
        self.map.color(if self.reversed { 1.0 - t } else { t })
    }

    /// The shade of `value` on `[low, high]`.
    pub fn shade(&self, value: f64, low: f64, high: f64) -> Shade {
        if self.mark_out_of_range && value.is_finite() && low < high {
            let slack = (high - low) * 1e-9;
            let (below, above) = self.map.out_of_range();
            if value < low - slack {
                return Shade::Below(below);
            }
            if value > high + slack {
                return Shade::Above(above);
            }
        }
        Shade::In(self.at(ColorRange { min: low, max: high, locked: true }.normalize(value)))
    }
}

/// "Nice" legend ticks (1, 2, 2.5, 5 × 10ⁿ steps) covering `[min, max]`, at most `target + 2`.
pub fn legend_ticks(min: f64, max: f64, target: usize) -> Vec<f64> {
    if !(min.is_finite() && max.is_finite()) || target == 0 {
        return vec![];
    }
    if max <= min {
        return vec![min];
    }
    let raw = (max - min) / target as f64;
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 2.5, 5.0, 10.0].iter().map(|m| m * mag).find(|s| *s >= raw).unwrap_or(10.0 * mag);
    let first = (min / step).ceil();
    // Round to the step's decimal places so 3 × 0.2 prints as 0.6.
    let places = (-(step.log10().floor()) + 1.0).clamp(0.0, 15.0) as i32;
    let scale = 10f64.powi(places);
    let mut ticks = Vec::new();
    let mut k = 0usize;
    while k <= target + 2 {
        let v = ((first + k as f64) * step * scale).round() / scale;
        if v > max + 1e-9 * step {
            break;
        }
        ticks.push(if v.abs() < 1e-12 * step { 0.0 } else { v });
        k += 1;
    }
    ticks
}

// ---------------------------------------------------------------- units

/// A display unit inferred from a companion field or column name, or `None` when the name
/// carries no unit. Flow `pressure` is kinematic (pressure divided by density).
pub fn unit_hint(name: &str) -> Option<&'static str> {
    let lower = name.to_ascii_lowercase();
    if lower == "pressure" || lower == "p" || lower == "kinematic_pressure" {
        return Some("m²/s² (kinematic)");
    }
    const SUFFIXES: &[(&str, &str)] = &[
        ("_w_m3", "W/m³"),
        ("_w_m2", "W/m²"),
        ("_w_kg", "W/kg"),
        ("_w_per_m_k", "W/(m·K)"),
        ("_a_m2", "A/m²"),
        ("_a_per_m2", "A/m²"),
        ("_v_per_m", "V/m"),
        ("_v_m", "V/m"),
        ("_a_per_m", "A/m"),
        ("_c_m3", "C/m³"),
        ("_c_per_m3", "C/m³"),
        ("_m_per_s", "m/s"),
        ("_m_s", "m/s"),
        ("_m3_per_s", "m³/s"),
        ("_m3", "m³"),
        ("_m2", "m²"),
        ("_kg_m3", "kg/m³"),
        ("_pa", "Pa"),
        ("_k", "K"),
        ("_c", "°C"),
        ("_w", "W"),
        ("_v", "V"),
        ("_a", "A"),
        ("_j", "J"),
        ("_n", "N"),
        ("_s", "s"),
        ("_m", "m"),
        ("_hz", "Hz"),
        ("_ohm", "Ω"),
        ("_db", "dB"),
        ("_deg", "°"),
        ("_gy", "Gy"),
        ("_ah", "Ah"),
        ("_wh", "Wh"),
        ("_f", "F"),
    ];
    if lower.ends_with("temperature") {
        return Some("K");
    }
    SUFFIXES.iter().find(|(s, _)| lower.ends_with(s)).map(|(_, u)| *u)
}

// ---------------------------------------------------------------- cell geometry

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(sub(a, b))
}

/// Simplices of a cell as point-index lists (tetrahedra for solids, triangles for faces,
/// segments for lines), using VTK node ordering.
fn simplices(kind: u8, nodes: &[usize]) -> Vec<Vec<usize>> {
    let n = |i: usize| nodes[i];
    match (kind, nodes.len()) {
        (1, 1) => vec![vec![n(0)]],
        (3, 2) => vec![vec![n(0), n(1)]],
        (5, 3) => vec![vec![n(0), n(1), n(2)]],
        (9, 4) => vec![vec![n(0), n(1), n(2)], vec![n(0), n(2), n(3)]],
        (10, 4) => vec![nodes.to_vec()],
        (12, 8) => [[0, 1, 3, 4], [1, 2, 3, 6], [1, 3, 4, 6], [1, 4, 5, 6], [3, 4, 6, 7]].iter().map(|t| t.iter().map(|&i| n(i)).collect()).collect(),
        (13, 6) => [[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]].iter().map(|t| t.iter().map(|&i| n(i)).collect()).collect(),
        (14, 5) => vec![vec![n(0), n(1), n(2), n(4)], vec![n(0), n(2), n(3), n(4)]],
        _ => vec![],
    }
}

/// Dimension (0–3) of a supported cell type, `None` otherwise.
fn cell_dimension(kind: u8) -> Option<usize> {
    match kind {
        1 => Some(0),
        3 => Some(1),
        5 | 9 => Some(2),
        10 | 12 | 13 | 14 => Some(3),
        _ => None,
    }
}

fn simplex_measure(p: &[[f64; 3]]) -> f64 {
    match p.len() {
        2 => dist(p[0], p[1]),
        3 => 0.5 * norm(cross(sub(p[1], p[0]), sub(p[2], p[0]))),
        4 => dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))).abs() / 6.0,
        _ => 0.0,
    }
}

/// Barycentric weights of `x` in a segment, triangle or tetrahedron (projected onto the
/// simplex's span) and the distance from `x` to that span.
pub(crate) fn barycentric(p: &[[f64; 3]], x: [f64; 3]) -> Option<(Vec<f64>, f64)> {
    match p.len() {
        1 => Some((vec![1.0], dist(p[0], x))),
        2 => {
            let e = sub(p[1], p[0]);
            let l2 = dot(e, e);
            if l2 <= 0.0 {
                return None;
            }
            let t = dot(sub(x, p[0]), e) / l2;
            let q = [p[0][0] + t * e[0], p[0][1] + t * e[1], p[0][2] + t * e[2]];
            Some((vec![1.0 - t, t], dist(q, x)))
        }
        3 => {
            let (e1, e2) = (sub(p[1], p[0]), sub(p[2], p[0]));
            let nrm = cross(e1, e2);
            let n2 = dot(nrm, nrm);
            if n2 <= 0.0 {
                return None;
            }
            let w = sub(x, p[0]);
            let off = dot(w, nrm) / n2.sqrt();
            let b1 = dot(cross(w, e2), nrm) / n2;
            let b2 = dot(cross(e1, w), nrm) / n2;
            Some((vec![1.0 - b1 - b2, b1, b2], off.abs()))
        }
        4 => {
            let (e1, e2, e3) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
            let det = dot(e1, cross(e2, e3));
            if det.abs() <= 0.0 {
                return None;
            }
            let w = sub(x, p[0]);
            let b1 = dot(w, cross(e2, e3)) / det;
            let b2 = dot(e1, cross(w, e3)) / det;
            let b3 = dot(e1, cross(e2, w)) / det;
            Some((vec![1.0 - b1 - b2 - b3, b1, b2, b3], 0.0))
        }
        _ => None,
    }
}

/// Precomputed cell geometry for repeated queries on one grid.
#[derive(Clone, Debug)]
pub struct CellIndex {
    pub(crate) nodes: Vec<Vec<usize>>,
    pub(crate) simplices: Vec<Vec<Vec<usize>>>,
    pub centroids: Vec<[f64; 3]>,
    /// Length, area or volume by cell dimension.
    pub measures: Vec<f64>,
    pub dimensions: Vec<Option<usize>>,
    pub(crate) bounds: Vec<[[f64; 2]; 3]>,
    /// Highest cell dimension present (integrals use only cells of this dimension).
    pub dimension: usize,
}

impl CellIndex {
    pub fn new(grid: &Grid) -> Result<Self, Error> {
        let cells = grid.cell_types.len();
        if grid.offsets.len() != cells {
            return Err(bad("offsets and cell types differ in length"));
        }
        let mut idx =
            CellIndex { nodes: vec![], simplices: vec![], centroids: vec![], measures: vec![], dimensions: vec![], bounds: vec![], dimension: 0 };
        let mut start = 0usize;
        for c in 0..cells {
            let end = grid.offsets[c];
            if end < start || end > grid.connectivity.len() {
                return Err(bad("invalid cell offsets"));
            }
            let nodes = grid.connectivity[start..end].to_vec();
            start = end;
            if nodes.iter().any(|&n| n >= grid.points.len()) || nodes.is_empty() {
                return Err(bad("cell references a missing point"));
            }
            let kind = grid.cell_types[c];
            let simp = simplices(kind, &nodes);
            let pts: Vec<[f64; 3]> = nodes.iter().map(|&n| grid.points[n]).collect();
            let k = pts.len() as f64;
            let centroid =
                [pts.iter().map(|p| p[0]).sum::<f64>() / k, pts.iter().map(|p| p[1]).sum::<f64>() / k, pts.iter().map(|p| p[2]).sum::<f64>() / k];
            let measure: f64 = simp.iter().map(|s| simplex_measure(&s.iter().map(|&n| grid.points[n]).collect::<Vec<_>>())).sum();
            let mut b = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
            for p in &pts {
                for a in 0..3 {
                    b[a][0] = b[a][0].min(p[a]);
                    b[a][1] = b[a][1].max(p[a]);
                }
            }
            let dim = cell_dimension(kind);
            idx.dimension = idx.dimension.max(dim.unwrap_or(0));
            idx.nodes.push(nodes);
            idx.simplices.push(simp);
            idx.centroids.push(centroid);
            idx.measures.push(measure);
            idx.dimensions.push(dim);
            idx.bounds.push(b);
        }
        Ok(idx)
    }

    pub fn len(&self) -> usize {
        self.centroids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.centroids.is_empty()
    }

    /// The cell containing `x` (within `tolerance`) with interpolation weights over the
    /// cell's simplex nodes, or `None`.
    fn locate(&self, grid: &Grid, x: [f64; 3], tolerance: f64) -> Option<(usize, Vec<(usize, f64)>)> {
        for c in 0..self.len() {
            let b = &self.bounds[c];
            if (0..3).any(|a| x[a] < b[a][0] - tolerance || x[a] > b[a][1] + tolerance) {
                continue;
            }
            for s in &self.simplices[c] {
                let pts: Vec<[f64; 3]> = s.iter().map(|&n| grid.points[n]).collect();
                if let Some((w, off)) = barycentric(&pts, x) {
                    let scale = tolerance / (self.measures[c].abs().powf(1.0 / self.dimensions[c].unwrap_or(1).max(1) as f64)).max(f64::MIN_POSITIVE);
                    if off <= tolerance && w.iter().all(|wi| *wi >= -1e-9 - scale) {
                        return Some((c, s.iter().copied().zip(w).collect()));
                    }
                }
            }
        }
        None
    }

    fn nearest(&self, x: [f64; 3]) -> Option<usize> {
        (0..self.len()).min_by(|&a, &b| dist(self.centroids[a], x).total_cmp(&dist(self.centroids[b], x)))
    }
}

pub(crate) fn field<'a>(grid: &'a Grid, name: &str) -> Result<&'a Field, Error> {
    grid.fields.iter().find(|f| f.name == name).ok_or_else(|| bad(format!("no field named {name:?}")))
}

/// Scalar value of a field entry: the value itself, or the magnitude of a vector.
fn scalar(f: &Field, index: usize) -> Option<f64> {
    let k = f.components.max(1);
    let slice = f.values.get(index * k..index * k + k)?;
    Some(if k == 1 { slice[0] } else { slice.iter().map(|v| v * v).sum::<f64>().sqrt() })
}

pub(crate) fn components(f: &Field, index: usize) -> Option<Vec<f64>> {
    let k = f.components.max(1);
    f.values.get(index * k..index * k + k).map(<[f64]>::to_vec)
}

pub(crate) fn check_field(grid: &Grid, f: &Field) -> Result<(), Error> {
    let k = f.components.max(1);
    let count = match f.association {
        Association::Point => grid.points.len(),
        Association::Cell => grid.cell_types.len(),
    };
    if f.values.len() != count * k {
        return Err(bad(format!("field {:?} has {} values, expected {}", f.name, f.values.len(), count * k)));
    }
    Ok(())
}

/// A cell field's mean over a cell's nodes is used for point fields where one value per
/// cell is needed.
pub(crate) fn cell_value(grid: &Grid, idx: &CellIndex, f: &Field, cell: usize) -> Option<Vec<f64>> {
    match f.association {
        Association::Cell => components(f, cell),
        Association::Point => {
            let k = f.components.max(1);
            let nodes = &idx.nodes[cell];
            let mut acc = vec![0.0; k];
            for &n in nodes {
                for (a, v) in acc.iter_mut().zip(components(f, n)?) {
                    *a += v;
                }
            }
            let _ = grid;
            Some(acc.into_iter().map(|a| a / nodes.len() as f64).collect())
        }
    }
}

pub(crate) fn magnitude(v: &[f64]) -> f64 {
    if v.len() == 1 { v[0] } else { v.iter().map(|x| x * x).sum::<f64>().sqrt() }
}

// ---------------------------------------------------------------- probe

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct Probe {
    pub point: [f64; 3],
    /// Containing cell, or the nearest cell when `inside` is false.
    pub cell: usize,
    pub inside: bool,
    /// Distance from the point to the reported cell's centroid.
    pub distance: f64,
    /// Field components at the point (one for scalars).
    pub value: Vec<f64>,
    /// Scalar value or vector magnitude.
    pub magnitude: f64,
    /// "linear" inside a cell for point data, "cell" for cell data, "nearest" outside.
    pub method: String,
}

/// Value of `field_name` at `point`. Inside the mesh, point data are interpolated linearly
/// and cell data take the containing cell's value; outside, the nearest cell is reported
/// with `inside = false` and nothing is extrapolated.
pub fn probe(grid: &Grid, index: &CellIndex, field_name: &str, point: [f64; 3]) -> Result<Probe, Error> {
    let f = field(grid, field_name)?;
    check_field(grid, f)?;
    if !point.iter().all(|v| v.is_finite()) {
        return Err(bad("probe point must be finite"));
    }
    let extent = grid.points.iter().fold(0.0f64, |m, p| m.max(norm(*p)));
    let tolerance = 1e-9 * extent.max(1e-12);
    if let Some((cell, weights)) = index.locate(grid, point, tolerance) {
        let (value, method) = match f.association {
            Association::Cell => (components(f, cell).ok_or_else(|| bad("cell value missing"))?, "cell"),
            Association::Point => {
                let k = f.components.max(1);
                let mut acc = vec![0.0; k];
                for (n, w) in &weights {
                    for (a, v) in acc.iter_mut().zip(components(f, *n).ok_or_else(|| bad("point value missing"))?) {
                        *a += w * v;
                    }
                }
                (acc, "linear")
            }
        };
        return Ok(Probe {
            point,
            cell,
            inside: true,
            distance: dist(index.centroids[cell], point),
            magnitude: magnitude(&value),
            value,
            method: method.into(),
        });
    }
    let cell = index.nearest(point).ok_or_else(|| bad("the grid has no cells"))?;
    let value = cell_value(grid, index, f, cell).ok_or_else(|| bad("value missing"))?;
    Ok(Probe {
        point,
        cell,
        inside: false,
        distance: dist(index.centroids[cell], point),
        magnitude: magnitude(&value),
        value,
        method: "nearest".into(),
    })
}

// ---------------------------------------------------------------- min / max

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct Extremum {
    pub value: f64,
    /// Point or cell index, by the field's association.
    pub index: usize,
    /// Point coordinate or cell centroid.
    pub location: [f64; 3],
}

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct MinMax {
    pub field: String,
    pub min: Extremum,
    pub max: Extremum,
}

/// Minimum and maximum (vector magnitude for vectors) with their locations.
pub fn min_max(grid: &Grid, index: &CellIndex, field_name: &str) -> Result<MinMax, Error> {
    let f = field(grid, field_name)?;
    check_field(grid, f)?;
    let count = f.values.len() / f.components.max(1);
    let location = |i: usize| match f.association {
        Association::Point => grid.points[i],
        Association::Cell => index.centroids[i],
    };
    let mut lo: Option<(usize, f64)> = None;
    let mut hi: Option<(usize, f64)> = None;
    for i in 0..count {
        let Some(v) = scalar(f, i).filter(|v| v.is_finite()) else { continue };
        if lo.is_none_or(|(_, m)| v < m) {
            lo = Some((i, v));
        }
        if hi.is_none_or(|(_, m)| v > m) {
            hi = Some((i, v));
        }
    }
    let ((li, lv), (hi_i, hv)) = lo.zip(hi).ok_or_else(|| bad("field has no finite values"))?;
    Ok(MinMax {
        field: field_name.into(),
        min: Extremum { value: lv, index: li, location: location(li) },
        max: Extremum { value: hv, index: hi_i, location: location(hi_i) },
    })
}

// ---------------------------------------------------------------- line samples

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct LineSample {
    /// Distance from the line start.
    pub s: f64,
    pub point: [f64; 3],
    /// `None` where the point lies outside the mesh.
    pub value: Option<f64>,
}

/// `samples` evenly spaced points from `a` to `b` (both included); values are scalar or
/// vector magnitude, and points outside the mesh have no value (no nearest-cell fill).
pub fn sample_line(grid: &Grid, index: &CellIndex, field_name: &str, a: [f64; 3], b: [f64; 3], samples: usize) -> Result<Vec<LineSample>, Error> {
    if !(2..=MAX_LINE_SAMPLES).contains(&samples) {
        return Err(bad(format!("use 2 to {MAX_LINE_SAMPLES} line samples")));
    }
    let length = dist(a, b);
    if !(length.is_finite() && length > 0.0) {
        return Err(bad("line end points must be distinct and finite"));
    }
    let mut out = Vec::with_capacity(samples);
    for i in 0..samples {
        let t = i as f64 / (samples - 1) as f64;
        let point = [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])];
        let p = probe(grid, index, field_name, point)?;
        out.push(LineSample { s: t * length, point, value: p.inside.then_some(p.magnitude) });
    }
    Ok(out)
}

// ---------------------------------------------------------------- integrals

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct RegionStatistics {
    pub field: String,
    /// Cells of the grid's highest dimension whose centroid lies in the region.
    pub cells: usize,
    /// Total length, area or volume of those cells.
    pub measure: f64,
    /// ∫ value dV (vector magnitude for vectors).
    pub integral: f64,
    /// Measure-weighted mean.
    pub mean: f64,
    pub min: f64,
    pub max: f64,
}

/// Measure-weighted statistics over cells of the highest dimension, optionally limited to an
/// axis-aligned box `[[xmin, xmax], [ymin, ymax], [zmin, zmax]]` by cell centroid. Point
/// data use the mean of each cell's nodes.
pub fn integrate(grid: &Grid, index: &CellIndex, field_name: &str, region: Option<[[f64; 2]; 3]>) -> Result<RegionStatistics, Error> {
    let f = field(grid, field_name)?;
    check_field(grid, f)?;
    let mut s =
        RegionStatistics { field: field_name.into(), cells: 0, measure: 0.0, integral: 0.0, mean: 0.0, min: f64::INFINITY, max: f64::NEG_INFINITY };
    for c in 0..index.len() {
        if index.dimensions[c] != Some(index.dimension) {
            continue;
        }
        let x = index.centroids[c];
        if let Some(r) = region
            && (0..3).any(|a| x[a] < r[a][0] || x[a] > r[a][1])
        {
            continue;
        }
        let v = magnitude(&cell_value(grid, index, f, c).ok_or_else(|| bad("value missing"))?);
        if !v.is_finite() {
            continue;
        }
        s.cells += 1;
        s.measure += index.measures[c];
        s.integral += v * index.measures[c];
        s.min = s.min.min(v);
        s.max = s.max.max(v);
    }
    if s.cells == 0 {
        return Err(bad("no cells lie in the region"));
    }
    s.mean = if s.measure > 0.0 { s.integral / s.measure } else { f64::NAN };
    Ok(s)
}

// ---------------------------------------------------------------- comparison

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct FieldDifference {
    pub field: String,
    pub entries: usize,
    pub max_abs: f64,
    /// Location of the largest difference.
    pub max_abs_location: [f64; 3],
    pub rms: f64,
    /// ‖b − a‖₂ / ‖a‖₂ (NaN when ‖a‖ is zero).
    pub relative_l2: f64,
    /// `b − a` per entry (scalars, or vector magnitudes), for a difference plot.
    pub difference: Vec<f64>,
}

/// Difference of one field between two results on the same mesh (`b − a`). The meshes
/// must have the same points (within 1e-9 relative) and cells; otherwise this is an error,
/// not an interpolation.
pub fn compare(a: &Grid, b: &Grid, field_name: &str) -> Result<FieldDifference, Error> {
    if a.points.len() != b.points.len() || a.cell_types != b.cell_types || a.connectivity != b.connectivity || a.offsets != b.offsets {
        return Err(bad("the two results are on different meshes; compare needs identical meshes"));
    }
    let extent = a.points.iter().fold(0.0f64, |m, p| m.max(norm(*p))).max(1e-300);
    if a.points.iter().zip(&b.points).any(|(p, q)| dist(*p, *q) > 1e-9 * extent) {
        return Err(bad("the two results have different point coordinates"));
    }
    let (fa, fb) = (field(a, field_name)?, field(b, field_name)?);
    check_field(a, fa)?;
    check_field(b, fb)?;
    if fa.association != fb.association || fa.components != fb.components {
        return Err(bad("the field differs in association or components"));
    }
    let index = CellIndex::new(a)?;
    let count = fa.values.len() / fa.components.max(1);
    let mut d = FieldDifference {
        field: field_name.into(),
        entries: count,
        max_abs: 0.0,
        max_abs_location: [0.0; 3],
        rms: 0.0,
        relative_l2: 0.0,
        difference: Vec::with_capacity(count),
    };
    let (mut sum2, mut ref2) = (0.0, 0.0);
    for i in 0..count {
        let (va, vb) = (components(fa, i).unwrap_or_default(), components(fb, i).unwrap_or_default());
        let diff: Vec<f64> = va.iter().zip(&vb).map(|(x, y)| y - x).collect();
        let m = magnitude(&diff).abs();
        let signed = if diff.len() == 1 { diff[0] } else { m };
        d.difference.push(signed);
        sum2 += diff.iter().map(|x| x * x).sum::<f64>();
        ref2 += va.iter().map(|x| x * x).sum::<f64>();
        if m > d.max_abs {
            d.max_abs = m;
            d.max_abs_location = match fa.association {
                Association::Point => a.points[i],
                Association::Cell => index.centroids[i],
            };
        }
    }
    d.rms = (sum2 / count.max(1) as f64).sqrt();
    d.relative_l2 = if ref2 > 0.0 { (sum2 / ref2).sqrt() } else { f64::NAN };
    Ok(d)
}

// ---------------------------------------------------------------- CSV tables

/// A numeric CSV table as written by the companion (`write_csv`: quoted header, numbers).
#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct Table {
    pub columns: Vec<String>,
    pub rows: Vec<Vec<f64>>,
}

impl Table {
    pub fn column(&self, name: &str) -> Option<Vec<f64>> {
        let c = self.columns.iter().position(|n| n == name)?;
        Some(self.rows.iter().map(|r| r[c]).collect())
    }
}

/// Parse a numeric CSV with one header row. Quoted header names are unquoted; a cell that
/// is not a number makes the whole table an error.
pub fn parse_csv(text: &str) -> Result<Table, Error> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = lines.next().ok_or_else(|| bad("empty CSV"))?;
    let columns: Vec<String> = header.split(',').map(|h| h.trim().trim_matches('"').replace("\"\"", "\"")).collect();
    let mut rows = Vec::new();
    for (i, line) in lines.enumerate() {
        if rows.len() >= MAX_TABLE_ROWS {
            return Err(bad("CSV exceeds the row limit"));
        }
        let row: Vec<f64> =
            line.split(',').map(|v| v.trim().parse::<f64>()).collect::<Result<_, _>>().map_err(|_| bad(format!("row {} is not numeric", i + 2)))?;
        if row.len() != columns.len() {
            return Err(bad(format!("row {} has {} values for {} columns", i + 2, row.len(), columns.len())));
        }
        rows.push(row);
    }
    Ok(Table { columns, rows })
}

// ---------------------------------------------------------------- goals

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct GoalHistory {
    pub name: String,
    pub kind: Option<String>,
    pub unit: Option<String>,
    pub steps: Vec<f64>,
    pub values: Vec<f64>,
    pub final_value: Option<f64>,
    pub last_relative_change: Option<f64>,
    pub converged: Option<bool>,
    /// Target value from the case, when one was set.
    pub target: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct GoalReport {
    pub goals: Vec<GoalHistory>,
    pub all_converged: Option<bool>,
    pub stopped_early_at_step: Option<f64>,
    pub relative_tolerance: Option<f64>,
    pub rule: Option<String>,
}

/// Combine the companion's `goals.csv` (history) and optional `goals.json` (verdicts).
pub fn goals(csv: &str, json: Option<&str>) -> Result<GoalReport, Error> {
    let table = parse_csv(csv)?;
    let steps = table.column("step").ok_or_else(|| bad("goals.csv has no step column"))?;
    let meta: Option<Value> = match json {
        Some(text) => Some(serde_json::from_str(text).map_err(|e| bad(format!("goals.json: {e}")))?),
        None => None,
    };
    let entry = |name: &str| meta.as_ref().and_then(|m| m["goals"].as_array()).and_then(|a| a.iter().find(|g| g["name"] == name)).cloned();
    let mut goals = Vec::new();
    for name in table.columns.iter().filter(|c| c.as_str() != "step") {
        let values = table.column(name).unwrap_or_default();
        let e = entry(name);
        let get = |k: &str| e.as_ref().map(|g| g[k].clone()).filter(|v| !v.is_null());
        goals.push(GoalHistory {
            name: name.clone(),
            kind: get("kind").and_then(|v| v.as_str().map(String::from)),
            unit: get("unit").and_then(|v| v.as_str().map(String::from)),
            steps: steps.clone(),
            final_value: get("final_value").and_then(|v| v.as_f64()).or_else(|| values.last().copied()),
            values,
            last_relative_change: get("last_relative_change").and_then(|v| v.as_f64()),
            converged: get("converged").and_then(|v| v.as_bool()),
            target: get("target").and_then(|v| v.as_f64()),
        });
    }
    let m = meta.as_ref();
    Ok(GoalReport {
        goals,
        all_converged: m.and_then(|m| m["all_converged"].as_bool()),
        stopped_early_at_step: m.and_then(|m| m["stopped_early_at_step"].as_f64()),
        relative_tolerance: m.and_then(|m| m["control"]["relative_tolerance"].as_f64()),
        rule: m.and_then(|m| m["rule"].as_str().map(String::from)),
    })
}

// ---------------------------------------------------------------- components

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, serde::Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ComponentStatus {
    Pass,
    Warn,
    Fail,
    NoLimit,
}

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct ComponentRow {
    pub name: String,
    pub power_w: f64,
    pub min_temperature_k: f64,
    pub mean_temperature_k: f64,
    pub max_temperature_k: f64,
    pub limit_k: Option<f64>,
    pub margin_k: Option<f64>,
    pub status: ComponentStatus,
}

/// Component temperatures and verdicts from the companion's `thermal_report.json`.
pub fn components_from_report(report: &str) -> Result<Vec<ComponentRow>, Error> {
    let v: Value = serde_json::from_str(report).map_err(|e| bad(format!("thermal_report.json: {e}")))?;
    let list = v["components"].as_array().ok_or_else(|| bad("thermal_report.json has no components"))?;
    let num = |c: &Value, k: &str| c[k].as_f64().ok_or_else(|| bad(format!("component field {k} is missing")));
    list.iter()
        .take(MAX_TABLE_ROWS)
        .map(|c| {
            let limit = c["limit_k"].as_f64();
            Ok(ComponentRow {
                name: c["name"].as_str().unwrap_or("component").into(),
                power_w: num(c, "power_w")?,
                min_temperature_k: num(c, "min_temperature_k")?,
                mean_temperature_k: num(c, "mean_temperature_k")?,
                max_temperature_k: num(c, "max_temperature_k")?,
                limit_k: limit,
                margin_k: c["margin_k"].as_f64(),
                status: match (limit, c["status"].as_str()) {
                    (None, _) => ComponentStatus::NoLimit,
                    (_, Some("fail")) => ComponentStatus::Fail,
                    (_, Some("warn")) => ComponentStatus::Warn,
                    _ => ComponentStatus::Pass,
                },
            })
        })
        .collect()
}

// ---------------------------------------------------------------- metric tables

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct MetricTable {
    /// One label per run (study name and run id).
    pub runs: Vec<String>,
    pub metrics: Vec<String>,
    /// `values[metric][run]`; `None` where a run lacks the metric.
    pub values: Vec<Vec<Option<f64>>>,
}

/// Metric table across runs from their manifests' `metrics` objects, metrics sorted by name.
pub fn metric_table(runs: &[(String, &Value)]) -> MetricTable {
    let mut names = std::collections::BTreeSet::new();
    for (_, m) in runs {
        if let Some(obj) = m["metrics"].as_object() {
            names.extend(obj.keys().cloned());
        }
    }
    let metrics: Vec<String> = names.into_iter().collect();
    let values = metrics.iter().map(|k| runs.iter().map(|(_, m)| m["metrics"][k].as_f64()).collect()).collect();
    MetricTable { runs: runs.iter().map(|(l, _)| l.clone()).collect(), metrics, values }
}

// ---------------------------------------------------------------- queries on retained runs

/// One question about a retained run. Field queries read a retained ASCII VTU artifact.
#[derive(Clone, Debug, PartialEq, serde::Deserialize, Serialize, schemars::JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Query {
    Probe {
        artifact: String,
        field: String,
        point: [f64; 3],
    },
    MinMax {
        artifact: String,
        field: String,
    },
    Line {
        artifact: String,
        field: String,
        start: [f64; 3],
        end: [f64; 3],
        samples: usize,
    },
    Integrate {
        artifact: String,
        field: String,
        #[serde(default)]
        region: Option<[[f64; 2]; 3]>,
    },
    /// `b − a` with this run as `a` and `other_run_id` (of `other_study`, default the same
    /// study) as `b`; both must hold `artifact` on the same mesh.
    Compare {
        artifact: String,
        field: String,
        other_run_id: String,
        #[serde(default)]
        other_study: Option<String>,
    },
    Goals,
    Components,
    /// Flow trajectories (or particle paths) of a vector field.
    Trajectories {
        artifact: String,
        field: String,
        seeds: crate::flowpost::Seeds,
        #[serde(default)]
        options: crate::flowpost::TrajectoryOptions,
    },
    Isosurface {
        artifact: String,
        field: String,
        level: f64,
        #[serde(default)]
        color_field: Option<String>,
    },
    /// The mesh boundary (or one side of its bounding box) with the field on it.
    Surface {
        artifact: String,
        field: String,
        #[serde(default)]
        side: Option<crate::flowpost::BoxSide>,
    },
}

#[derive(Clone, Debug, PartialEq, serde::Deserialize, Serialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct QueryRequest {
    pub id: String,
    pub run_id: String,
    pub query: Query,
}

#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum QueryResult {
    Probe(Probe),
    MinMax(MinMax),
    Line { samples: Vec<LineSample> },
    Integrate(RegionStatistics),
    Compare(FieldDifference),
    Goals(GoalReport),
    Components { rows: Vec<ComponentRow> },
    Trajectories { lines: Vec<crate::flowpost::Trajectory> },
    Isosurface(crate::flowpost::Isosurface),
    Surface(crate::flowpost::Surface),
}

/// A query answer with the provenance needed to trust it.
#[derive(Clone, Debug, PartialEq, Serialize, schemars::JsonSchema)]
pub struct QueryResponse {
    pub id: String,
    pub run_id: String,
    pub case_hash: String,
    /// The study changed after this run; the answer describes the run, not the current case.
    pub stale: bool,
    /// Unit inferred from the field name, if any.
    pub unit: Option<String>,
    pub result: QueryResult,
}

fn run_of<'a>(project: &'a crate::Project, study: &str, run: &str) -> Result<(&'a crate::Study, &'a crate::RunEvidence), Error> {
    let s = project.study(study).ok_or_else(|| bad(format!("study {study} not found")))?;
    let r = s.runs.iter().find(|r| r.id == run).ok_or_else(|| bad(format!("run {run} is not retained in study {study}")))?;
    Ok((s, r))
}

fn artifact<'a>(run: &'a crate::RunEvidence, name: &str) -> Result<&'a str, Error> {
    crate::validate_artifact_name(name).map_err(|e| bad(e.to_string()))?;
    run.artifacts
        .iter()
        .find(|a| a.name == name)
        .map(|a| a.text.as_str())
        .ok_or_else(|| bad(format!("artifact {name} is not retained in run {}", run.id)))
}

fn grid(run: &crate::RunEvidence, name: &str) -> Result<Grid, Error> {
    if !name.ends_with(".vtu") {
        return Err(bad("field queries need a VTU artifact"));
    }
    crate::vtk::parse_vtu(artifact(run, name)?).map_err(|e| bad(e.to_string()))
}

/// A zoom-in child case cut out of a retained run (see [`crate::zoom`]). The run's own input
/// case is the parent, so the child matches the results it is fed from even after later edits.
pub fn zoom_from_run(
    project: &crate::Project,
    id: &str,
    run_id: &str,
    artifact_name: &str,
    low: [f64; 3],
    high: [f64; 3],
    refine: usize,
) -> Result<crate::zoom::Zoom, Error> {
    let (_, run) = run_of(project, id, run_id)?;
    let g = grid(run, artifact_name)?;
    crate::zoom::zoom_case(&run.input_case, &g, low, high, refine).map_err(|e| bad(e.to_string()))
}

/// Answer one query about a retained run.
pub fn query(project: &crate::Project, request: &QueryRequest) -> Result<QueryResponse, Error> {
    let (study, run) = run_of(project, &request.id, &request.run_id)?;
    let hash = crate::case_hash(&study.case).map_err(|e| bad(e.to_string()))?;
    let stale = run.study_revision != study.revision || run.case_hash != hash;
    let (unit, result) = match &request.query {
        Query::Probe { artifact, field, point } => {
            let g = grid(run, artifact)?;
            let idx = CellIndex::new(&g)?;
            (unit_hint(field), QueryResult::Probe(probe(&g, &idx, field, *point)?))
        }
        Query::MinMax { artifact, field } => {
            let g = grid(run, artifact)?;
            let idx = CellIndex::new(&g)?;
            (unit_hint(field), QueryResult::MinMax(min_max(&g, &idx, field)?))
        }
        Query::Line { artifact, field, start, end, samples } => {
            let g = grid(run, artifact)?;
            let idx = CellIndex::new(&g)?;
            (unit_hint(field), QueryResult::Line { samples: sample_line(&g, &idx, field, *start, *end, *samples)? })
        }
        Query::Integrate { artifact, field, region } => {
            let g = grid(run, artifact)?;
            let idx = CellIndex::new(&g)?;
            (unit_hint(field), QueryResult::Integrate(integrate(&g, &idx, field, *region)?))
        }
        Query::Compare { artifact, field, other_run_id, other_study } => {
            let (_, other) = run_of(project, other_study.as_deref().unwrap_or(&request.id), other_run_id)?;
            let (a, b) = (grid(run, artifact)?, grid(other, artifact)?);
            (unit_hint(field), QueryResult::Compare(compare(&a, &b, field)?))
        }
        Query::Goals => {
            let csv = artifact(run, "goals.csv")?;
            let json = run.artifacts.iter().find(|a| a.name == "goals.json").map(|a| a.text.as_str());
            (None, QueryResult::Goals(goals(csv, json)?))
        }
        Query::Components => (Some("K"), QueryResult::Components { rows: components_from_report(artifact(run, "thermal_report.json")?)? }),
        Query::Trajectories { artifact, field, seeds, options } => {
            let g = grid(run, artifact)?;
            let idx = CellIndex::new(&g)?;
            let lines = crate::flowpost::trajectories(&g, &idx, field, &seeds.points(), options)?;
            (unit_hint(field), QueryResult::Trajectories { lines })
        }
        Query::Isosurface { artifact, field, level, color_field } => {
            let g = grid(run, artifact)?;
            let idx = CellIndex::new(&g)?;
            let unit = unit_hint(color_field.as_deref().unwrap_or(field));
            (unit, QueryResult::Isosurface(crate::flowpost::isosurface(&g, &idx, field, *level, color_field.as_deref())?))
        }
        Query::Surface { artifact, field, side } => {
            let g = grid(run, artifact)?;
            let idx = CellIndex::new(&g)?;
            (unit_hint(field), QueryResult::Surface(crate::flowpost::surface(&g, &idx, field, *side)?))
        }
    };
    Ok(QueryResponse {
        id: request.id.clone(),
        run_id: request.run_id.clone(),
        case_hash: run.case_hash.clone(),
        stale,
        unit: unit.map(String::from),
        result,
    })
}

/// Metric table over the most recent retained run of each listed study (all studies when
/// the list is empty), labelled `study name · run id`.
pub fn compare_studies(project: &crate::Project, ids: &[String]) -> Result<MetricTable, Error> {
    let studies: Vec<&crate::Study> = if ids.is_empty() {
        project.studies.iter().collect()
    } else {
        ids.iter().map(|id| project.study(id).ok_or_else(|| bad(format!("study {id} not found")))).collect::<Result<_, _>>()?
    };
    let runs: Vec<(String, &Value)> = studies.iter().filter_map(|s| s.runs.last().map(|r| (format!("{} · {}", s.name, r.id), &r.manifest))).collect();
    if runs.is_empty() {
        return Err(bad("none of these studies has a retained run"));
    }
    Ok(metric_table(&runs))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_cube_hex(field: Vec<f64>, association: Association) -> Grid {
        let points = vec![[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
        Grid {
            points,
            connectivity: (0..8).collect(),
            offsets: vec![8],
            cell_types: vec![12],
            fields: vec![Field { name: "T".into(), association, data_type: "Float64".into(), components: 1, values: field }],
        }
    }

    fn two_tets() -> Grid {
        // Unit cube split into a 2-cell strip along x: [0,1] and [1,2] hexes.
        let mut points = vec![];
        for z in [0.0, 1.0] {
            for y in [0.0, 1.0] {
                for x in [0.0, 1.0, 2.0] {
                    points.push([x, y, z]);
                }
            }
        }
        let id = |x: usize, y: usize, z: usize| z * 6 + y * 3 + x;
        let hex =
            |i: usize| vec![id(i, 0, 0), id(i + 1, 0, 0), id(i + 1, 1, 0), id(i, 1, 0), id(i, 0, 1), id(i + 1, 0, 1), id(i + 1, 1, 1), id(i, 1, 1)];
        let mut connectivity = hex(0);
        connectivity.extend(hex(1));
        let linear: Vec<f64> = points.iter().map(|p| 2.0 * p[0] + 3.0 * p[1] - p[2]).collect();
        Grid {
            points,
            connectivity,
            offsets: vec![8, 16],
            cell_types: vec![12, 12],
            fields: vec![
                Field { name: "phi".into(), association: Association::Point, data_type: "Float64".into(), components: 1, values: linear },
                Field { name: "q".into(), association: Association::Cell, data_type: "Float64".into(), components: 1, values: vec![1.0, 3.0] },
                Field {
                    name: "u".into(),
                    association: Association::Cell,
                    data_type: "Float64".into(),
                    components: 3,
                    values: vec![3.0, 4.0, 0.0, 0.0, 0.0, 1.0],
                },
            ],
        }
    }

    #[test]
    fn colour_maps_hit_their_end_stops_and_ranges_lock() {
        for map in ColorMap::ALL {
            assert_eq!(map.color(0.0), map.stops()[0]);
            assert_eq!(map.color(1.0), *map.stops().last().expect("test value"));
            assert_eq!(map.color(-3.0), map.color(0.0));
            assert_eq!(map.color(f64::NAN), map.color(0.0));
        }
        let r = ColorRange::from_values([3.0, f64::NAN, -1.0, 5.0]).expect("test value");
        assert_eq!((r.min, r.max), (-1.0, 5.0));
        assert!((r.normalize(2.0) - 0.5).abs() < 1e-15);
        assert_eq!(ColorRange { min: 2.0, max: 2.0, locked: true }.normalize(7.0), 0.5);
        assert!(ColorRange::from_values([f64::NAN]).is_none());
    }

    #[test]
    fn range_modes_clip_centre_and_lock() {
        let mut values: Vec<f64> = (0..=100).map(f64::from).collect();
        values.push(1e6);
        assert_eq!(range_for(RangeMode::Auto, &values, None), Some((0.0, 1e6)));
        // One extreme cell no longer sets the top of the scale.
        let (lo, hi) = range_for(RangeMode::Percentile, &values, None).unwrap();
        assert!((1.0..=3.0).contains(&lo) && (97.0..=100.0).contains(&hi), "{lo} {hi}");
        assert_eq!(range_for(RangeMode::Symmetric, &[-2.0, 5.0], None), Some((-5.0, 5.0)));
        assert_eq!(range_for(RangeMode::Locked, &values, Some((10.0, 20.0))), Some((10.0, 20.0)));
        assert_eq!(range_for(RangeMode::Locked, &values, Some((20.0, 10.0))), Some((0.0, 1e6)), "not a range: the data's");
        assert_eq!(range_for(RangeMode::Auto, &[f64::NAN], None), None);
    }

    #[test]
    fn scales_reverse_band_and_mark_values_out_of_range() {
        let plain = ColorScale::default();
        assert_eq!(plain.at(0.0), ColorMap::Viridis.color(0.0));
        let reversed = ColorScale { reversed: true, ..plain };
        assert_eq!(reversed.at(0.0), ColorMap::Viridis.color(1.0));
        // Four bands: 0.1 and 0.2 share a colour, 0.3 starts the next.
        let banded = ColorScale { bands: 4, ..plain };
        assert_eq!(banded.at(0.1), banded.at(0.2));
        assert_ne!(banded.at(0.2), banded.at(0.3));
        assert_eq!(banded.at(0.0), ColorMap::Viridis.color(0.125));
        // Clamped by default; marked when asked.
        assert_eq!(plain.shade(-5.0, 0.0, 1.0), Shade::In(ColorMap::Viridis.color(0.0)));
        let marked = ColorScale { mark_out_of_range: true, ..plain };
        assert_eq!(marked.shade(-5.0, 0.0, 1.0), Shade::Below(ColorMap::Viridis.out_of_range().0));
        assert_eq!(marked.shade(5.0, 0.0, 1.0), Shade::Above(ColorMap::Viridis.out_of_range().1));
        // No map shows its own out-of-range colours anywhere along it.
        for map in ColorMap::ALL {
            let (below, above) = map.out_of_range();
            for step in 0..=100 {
                let c = map.color(f64::from(step) / 100.0);
                let near = |o: [u8; 3]| (0..3).map(|i| (i32::from(c[i]) - i32::from(o[i])).abs()).sum::<i32>() < 60;
                assert!(!near(below) && !near(above), "{map:?} at {step}%: {c:?}");
            }
        }
        assert_eq!(marked.shade(1.0, 0.0, 1.0), Shade::In(ColorMap::Viridis.color(1.0)));
    }

    #[test]
    fn legend_ticks_are_nice_and_inside_the_range() {
        assert_eq!(legend_ticks(0.0, 1.0, 5), vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0]);
        let t = legend_ticks(293.4, 351.7, 5);
        assert_eq!(t, vec![300.0, 320.0, 340.0]);
        assert_eq!(legend_ticks(2.0, 2.0, 5), vec![2.0]);
        assert!(legend_ticks(f64::NAN, 1.0, 5).is_empty());
    }

    #[test]
    fn units_come_from_names_only() {
        assert_eq!(unit_hint("temperature_K"), Some("K"));
        assert_eq!(unit_hint("max_temperature_k"), Some("K"));
        assert_eq!(unit_hint("pressure"), Some("m²/s² (kinematic)"));
        assert_eq!(unit_hint("joule_heat_W_m3"), Some("W/m³"));
        assert_eq!(unit_hint("E_V_per_m"), Some("V/m"));
        assert_eq!(unit_hint("velocity"), None);
        assert_eq!(unit_hint("damage"), None);
    }

    #[test]
    fn probes_interpolate_linear_fields_exactly() {
        let g = two_tets();
        let idx = CellIndex::new(&g).expect("test value");
        assert_eq!(idx.dimension, 3);
        assert!((idx.measures.iter().sum::<f64>() - 2.0).abs() < 1e-12);
        for x in [[0.3, 0.7, 0.2], [1.5, 0.25, 0.9], [2.0, 1.0, 1.0], [1.0, 0.5, 0.5]] {
            let p = probe(&g, &idx, "phi", x).expect("test value");
            assert!(p.inside && p.method == "linear");
            assert!((p.magnitude - (2.0 * x[0] + 3.0 * x[1] - x[2])).abs() < 1e-12, "{x:?}");
        }
        let c = probe(&g, &idx, "q", [1.7, 0.5, 0.5]).expect("test value");
        assert_eq!((c.cell, c.value.clone(), c.method.as_str()), (1, vec![3.0], "cell"));
        let v = probe(&g, &idx, "u", [0.5, 0.5, 0.5]).expect("test value");
        assert_eq!(v.magnitude, 5.0);
        let out = probe(&g, &idx, "phi", [5.0, 0.5, 0.5]).expect("test value");
        assert!(!out.inside && out.cell == 1 && out.method == "nearest");
        assert!(probe(&g, &idx, "missing", [0.0; 3]).is_err());
        assert!(probe(&g, &idx, "phi", [f64::NAN, 0.0, 0.0]).is_err());
    }

    #[test]
    fn min_max_lines_integrals_and_comparison() {
        let g = two_tets();
        let idx = CellIndex::new(&g).expect("test value");
        let mm = min_max(&g, &idx, "phi").expect("test value");
        assert_eq!((mm.min.value, mm.max.value), (-1.0, 7.0));
        assert_eq!(mm.max.location, [2.0, 1.0, 0.0]);
        let q = min_max(&g, &idx, "q").expect("test value");
        assert_eq!(q.max.location, [1.5, 0.5, 0.5]);

        let line = sample_line(&g, &idx, "phi", [-0.5, 0.5, 0.5], [2.5, 0.5, 0.5], 7).expect("test value");
        assert_eq!(line.len(), 7);
        assert!(line[0].value.is_none() && line[6].value.is_none());
        assert!((line[3].value.expect("test value") - (2.0 + 1.5 - 0.5)).abs() < 1e-12);
        assert!(sample_line(&g, &idx, "phi", [0.0; 3], [0.0; 3], 5).is_err());
        assert!(sample_line(&g, &idx, "phi", [0.0; 3], [1.0; 3], 1).is_err());

        let all = integrate(&g, &idx, "q", None).expect("test value");
        assert_eq!(all.cells, 2);
        assert!((all.measure - 2.0).abs() < 1e-12 && (all.integral - 4.0).abs() < 1e-12 && (all.mean - 2.0).abs() < 1e-12);
        let left = integrate(&g, &idx, "q", Some([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])).expect("test value");
        assert!(left.cells == 1 && (left.integral - 1.0).abs() < 1e-12);
        // Point data: the mean of each hex's nodes of a linear field is its centroid value.
        let phi = integrate(&g, &idx, "phi", None).expect("test value");
        assert!((phi.integral - ((1.0 + 1.5 - 0.5) + (3.0 + 1.5 - 0.5))).abs() < 1e-12);
        assert!(integrate(&g, &idx, "q", Some([[5.0, 6.0], [0.0, 1.0], [0.0, 1.0]])).is_err());

        let mut h = g.clone();
        h.fields[1].values = vec![1.5, 3.0];
        let d = compare(&g, &h, "q").expect("test value");
        assert_eq!(d.difference, vec![0.5, 0.0]);
        assert_eq!(d.max_abs, 0.5);
        assert_eq!(d.max_abs_location, [0.5, 0.5, 0.5]);
        assert!((d.relative_l2 - (0.25f64 / 10.0).sqrt()).abs() < 1e-15);
        let mut moved = g.clone();
        moved.points[0][0] = -0.1;
        assert!(compare(&g, &moved, "q").is_err());
        assert!(compare(&g, &unit_cube_hex(vec![1.0], Association::Cell), "q").is_err());
    }

    #[test]
    fn triangles_and_tetrahedra_have_correct_measures() {
        let tri = Grid {
            points: vec![[0., 0., 0.], [2., 0., 0.], [0., 1., 0.]],
            connectivity: vec![0, 1, 2],
            offsets: vec![3],
            cell_types: vec![5],
            fields: vec![Field {
                name: "s".into(),
                association: Association::Point,
                data_type: "Float64".into(),
                components: 1,
                values: vec![0.0, 2.0, 1.0],
            }],
        };
        let idx = CellIndex::new(&tri).expect("test value");
        assert_eq!((idx.dimension, idx.measures[0]), (2, 1.0));
        let p = probe(&tri, &idx, "s", [0.5, 0.25, 0.0]).expect("test value");
        assert!(p.inside && (p.magnitude - 0.75).abs() < 1e-12);
        assert!(!probe(&tri, &idx, "s", [0.5, 0.25, 0.3]).expect("test value").inside, "off the plane");
        let tet = Grid {
            points: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            connectivity: vec![0, 1, 2, 3],
            offsets: vec![4],
            cell_types: vec![10],
            fields: vec![],
        };
        assert!((CellIndex::new(&tet).expect("test value").measures[0] - 1.0 / 6.0).abs() < 1e-15);
        let broken = Grid { offsets: vec![9], ..tet.clone() };
        assert!(CellIndex::new(&broken).is_err());
    }

    #[test]
    fn goal_and_component_files_parse() {
        let csv = "\"step\",\"q_out\",\"t_max\"\n1.0e0,2.0e0,3.0e2\n2.0e0,2.5e0,3.1e2\n";
        let json = r#"{"schema":"cfd-solvers-rs/goals/1","control":{"relative_tolerance":0.01},"all_converged":false,"stopped_early_at_step":null,"rule":"r","goals":[{"name":"q_out","kind":"patch_heat_out","target":null,"final_value":2.5,"last_relative_change":0.2,"converged":false,"unit":"W"}]}"#;
        let r = goals(csv, Some(json)).expect("test value");
        assert_eq!(r.goals.len(), 2);
        assert_eq!(r.goals[0].unit.as_deref(), Some("W"));
        assert_eq!(r.goals[0].converged, Some(false));
        assert_eq!(r.goals[1].values, vec![300.0, 310.0]);
        assert_eq!(r.goals[1].unit, None);
        assert_eq!(r.relative_tolerance, Some(0.01));
        assert!(goals("\"a\",\"b\"\n1,2\n", None).is_err(), "needs a step column");
        assert!(parse_csv("\"step\"\nx\n").is_err());
        assert!(parse_csv("\"step\",\"a\"\n1\n").is_err());

        let report = r#"{"components":[{"name":"cpu","power_w":5.0,"min_temperature_k":300.0,"mean_temperature_k":320.0,"max_temperature_k":340.0,"limit_k":358.15,"margin_k":18.15,"status":"pass"},{"name":"vrm","power_w":2.0,"min_temperature_k":300.0,"mean_temperature_k":350.0,"max_temperature_k":360.0,"limit_k":358.15,"margin_k":-1.85,"status":"fail"},{"name":"cap","power_w":0.1,"min_temperature_k":300.0,"mean_temperature_k":301.0,"max_temperature_k":302.0,"limit_k":null,"margin_k":null,"status":"no_limit"}]}"#;
        let rows = components_from_report(report).expect("test value");
        assert_eq!(rows.iter().map(|r| r.status).collect::<Vec<_>>(), vec![ComponentStatus::Pass, ComponentStatus::Fail, ComponentStatus::NoLimit]);
        assert!(components_from_report("{}").is_err());
    }

    #[test]
    fn metric_tables_align_runs() {
        let a = serde_json::json!({"metrics": {"t": 1.0, "p": 2.0}});
        let b = serde_json::json!({"metrics": {"t": 3.0, "q": 4.0}});
        let t = metric_table(&[("A".into(), &a), ("B".into(), &b)]);
        assert_eq!(t.metrics, vec!["p", "q", "t"]);
        assert_eq!(t.values, vec![vec![Some(2.0), None], vec![None, Some(4.0)], vec![Some(1.0), Some(3.0)]]);
    }

    #[test]
    fn queries_read_retained_runs_and_report_staleness() {
        let mut project = crate::Project::default();
        let study = project.create("flow-electronics-enclosure-cht", "Fan").expect("create");
        let goals_csv = "\"step\",\"t_max\"\n1,300\n2,301\n";
        let report = r#"{"components":[{"name":"cpu","power_w":5.0,"min_temperature_k":300.0,"mean_temperature_k":320.0,"max_temperature_k":340.0,"limit_k":358.15,"margin_k":18.15,"status":"pass"}]}"#;
        let hash = crate::case_hash(&study.case).expect("hash");
        let run = crate::RunEvidence {
            id: "run-1".into(),
            study_revision: study.revision,
            case_hash: hash,
            input_case: study.case.clone(),
            manifest: serde_json::json!({"metrics": {"t": 1.0}}),
            artifacts: vec![
                crate::Artifact { name: "goals.csv".into(), text: goals_csv.into() },
                crate::Artifact { name: "thermal_report.json".into(), text: report.into() },
            ],
        };
        project.studies[0].runs.push(run);
        let q = |query| QueryRequest { id: study.id.clone(), run_id: "run-1".into(), query };
        let g = query(&project, &q(Query::Goals)).expect("goals");
        assert!(!g.stale);
        assert!(matches!(g.result, QueryResult::Goals(ref r) if r.goals[0].values == vec![300.0, 301.0]));
        let c = query(&project, &q(Query::Components)).expect("components");
        assert!(matches!(c.result, QueryResult::Components { ref rows } if rows[0].name == "cpu"));
        assert!(query(&project, &q(Query::MinMax { artifact: "goals.csv".into(), field: "x".into() })).is_err(), "not a VTU");
        assert!(query(&project, &q(Query::MinMax { artifact: "../x.vtu".into(), field: "x".into() })).is_err());
        assert!(query(&project, &QueryRequest { id: study.id.clone(), run_id: "run-9".into(), query: Query::Goals }).is_err());
        project.studies[0].revision += 1;
        assert!(query(&project, &q(Query::Goals)).expect("goals").stale);
        let table = compare_studies(&project, &[]).expect("table");
        assert_eq!(table.runs, vec!["Fan · run-1"]);
        assert!(compare_studies(&project, &["mp-study-404".into()]).is_err());
    }

    #[test]
    fn flow_result_queries_trace_cut_and_integrate_retained_fields() {
        // One unit hexahedron: velocity (1, 0, 0) and t = x at its points.
        let vtu = r#"<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian"><UnstructuredGrid><Piece NumberOfPoints="8" NumberOfCells="1"><PointData><DataArray type="Float64" Name="velocity_m_per_s" NumberOfComponents="3" format="ascii">1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0</DataArray><DataArray type="Float64" Name="t" NumberOfComponents="1" format="ascii">0 1 1 0 0 1 1 0</DataArray></PointData><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1</DataArray></Points><Cells><DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3 4 5 6 7</DataArray><DataArray type="Int64" Name="offsets" format="ascii">8</DataArray><DataArray type="UInt8" Name="types" format="ascii">12</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>"#;
        let mut project = crate::Project::default();
        let study = project.create("flow-channel", "Flow").expect("create");
        let run = crate::RunEvidence {
            id: "run-1".into(),
            study_revision: study.revision,
            case_hash: crate::case_hash(&study.case).expect("hash"),
            input_case: study.case.clone(),
            manifest: serde_json::json!({}),
            artifacts: vec![crate::Artifact { name: "flow.vtu".into(), text: vtu.into() }],
        };
        project.studies[0].runs.push(run);
        let ask = |query| query_result(&project, &study.id, query);
        let request: Query = serde_json::from_value(serde_json::json!({
            "kind": "trajectories", "artifact": "flow.vtu", "field": "velocity_m_per_s",
            "seeds": {"kind": "line", "start": [0.1, 0.2, 0.5], "end": [0.1, 0.8, 0.5], "count": 3}
        }))
        .expect("trajectory request");
        let QueryResult::Trajectories { lines } = ask(request).expect("trajectories") else { panic!("trajectories expected") };
        assert_eq!(lines.len(), 3);
        assert!(lines.iter().all(|l| l.end == crate::flowpost::TrajectoryEnd::LeftMesh && (l.residence_time - 0.9).abs() < 1e-9));
        let iso = Query::Isosurface { artifact: "flow.vtu".into(), field: "t".into(), level: 0.5, color_field: None };
        let QueryResult::Isosurface(iso) = ask(iso).expect("isosurface") else { panic!("isosurface expected") };
        assert!((iso.area - 1.0).abs() < 1e-12);
        let outlet = Query::Surface {
            artifact: "flow.vtu".into(),
            field: "velocity_m_per_s".into(),
            side: Some(crate::flowpost::BoxSide { axis: 0, max: true }),
        };
        let QueryResult::Surface(outlet) = ask(outlet).expect("surface") else { panic!("surface expected") };
        assert_eq!(outlet.flux, Some(1.0));
        assert!(
            serde_json::from_value::<Query>(
                serde_json::json!({"kind":"surface","artifact":"flow.vtu","field":"t","side":{"axis":0,"max":true,"extra":1}})
            )
            .is_err()
        );
    }

    fn query_result(project: &crate::Project, id: &str, query: Query) -> Result<QueryResult, Error> {
        super::query(project, &QueryRequest { id: id.into(), run_id: "run-1".into(), query }).map(|r| r.result)
    }
}
