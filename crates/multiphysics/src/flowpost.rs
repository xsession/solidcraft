//! Flow-analysis result tools on a parsed results grid: flow trajectories (streamlines) and
//! particle studies, isosurfaces, and surface plots and parameters on the mesh boundary.
//!
//! Everything is bounded and works on [`crate::vtk::Grid`] with a [`CellIndex`]. Point data
//! are interpolated linearly inside cells (tetrahedral split, as for probes); cell data are
//! constant per cell. A uniform bin [`Locator`] keeps point location fast enough to trace
//! many lines.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::results::{CellIndex, Error, bad, barycentric, check_field, components, field, magnitude};
use crate::vtk::{Association, Field, Grid};

/// Seeds one trajectory request may trace.
pub const MAX_SEEDS: usize = 400;
/// Steps one trajectory may take in each direction.
pub const MAX_STEPS: usize = 20_000;
/// Points all trajectories of one request may return together.
pub const MAX_TRAJECTORY_POINTS: usize = 400_000;
/// Triangles one isosurface may return.
pub const MAX_ISO_TRIANGLES: usize = 400_000;
/// Faces one surface request may return.
pub const MAX_SURFACE_FACES: usize = 200_000;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
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

// ---------------------------------------------------------------- point location

/// Cells of the grid's highest dimension binned on a uniform grid over their bounds.
#[derive(Clone, Debug)]
pub struct Locator {
    min: [f64; 3],
    step: [f64; 3],
    bins: [usize; 3],
    cells: Vec<Vec<u32>>,
    tolerance: f64,
}

impl Locator {
    pub fn new(index: &CellIndex) -> Self {
        let solid: Vec<usize> = (0..index.len()).filter(|&c| index.dimensions[c] == Some(index.dimension)).collect();
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for &c in &solid {
            for a in 0..3 {
                min[a] = min[a].min(index.bounds[c][a][0]);
                max[a] = max[a].max(index.bounds[c][a][1]);
            }
        }
        if solid.is_empty() {
            return Self { min: [0.0; 3], step: [1.0; 3], bins: [1; 3], cells: vec![Vec::new()], tolerance: 0.0 };
        }
        let extent = [0, 1, 2].map(|a| (max[a] - min[a]).max(0.0));
        let largest = extent.iter().copied().fold(0.0, f64::max).max(f64::MIN_POSITIVE);
        let per_axis = (solid.len() as f64).cbrt().clamp(1.0, 128.0);
        let bins = extent.map(|e| ((e / largest * per_axis).round() as usize).clamp(1, 128));
        let step = [0, 1, 2].map(|a| if extent[a] > 0.0 { extent[a] / bins[a] as f64 } else { 1.0 });
        let mut cells = vec![Vec::new(); bins.iter().product()];
        let bin = |a: usize, v: f64| (((v - min[a]) / step[a]).floor().max(0.0) as usize).min(bins[a] - 1);
        for &c in &solid {
            let b = index.bounds[c];
            for i in bin(0, b[0][0])..=bin(0, b[0][1]) {
                for j in bin(1, b[1][0])..=bin(1, b[1][1]) {
                    for k in bin(2, b[2][0])..=bin(2, b[2][1]) {
                        if let Some(list) = cells.get_mut((k * bins[1] + j) * bins[0] + i) {
                            list.push(c as u32);
                        }
                    }
                }
            }
        }
        Self { min, step, bins, cells, tolerance: 1e-9 * largest }
    }

    fn candidates(&self, x: [f64; 3]) -> &[u32] {
        let mut ijk = [0usize; 3];
        for a in 0..3 {
            let t = (x[a] - self.min[a]) / self.step[a];
            if !t.is_finite() || t < -1e-9 || t > self.bins[a] as f64 + 1e-9 {
                return &[];
            }
            ijk[a] = (t.floor().max(0.0) as usize).min(self.bins[a] - 1);
        }
        self.cells.get((ijk[2] * self.bins[1] + ijk[1]) * self.bins[0] + ijk[0]).map_or(&[], Vec::as_slice)
    }

    /// The containing cell of `x` with weights over its simplex nodes.
    pub fn locate(&self, grid: &Grid, index: &CellIndex, x: [f64; 3]) -> Option<(usize, Vec<(usize, f64)>)> {
        for &c in self.candidates(x) {
            let c = c as usize;
            let b = index.bounds.get(c)?;
            let tol = self.tolerance;
            if (0..3).any(|a| x[a] < b[a][0] - tol || x[a] > b[a][1] + tol) {
                continue;
            }
            for s in index.simplices.get(c)? {
                let points: Option<Vec<[f64; 3]>> = s.iter().map(|&n| grid.points.get(n).copied()).collect();
                let points = points?;
                if points.len() != 4 {
                    continue;
                }
                if let Some((w, _)) = barycentric(&points, x)
                    && w.iter().all(|wi| *wi >= -1e-9)
                {
                    return Some((c, s.iter().copied().zip(w).collect()));
                }
            }
        }
        None
    }
}

/// Value of `f` at `x` (all components), or `None` outside the mesh.
pub(crate) fn sample(grid: &Grid, index: &CellIndex, locator: &Locator, f: &Field, x: [f64; 3]) -> Option<(usize, Vec<f64>)> {
    let (cell, weights) = locator.locate(grid, index, x)?;
    let value = match f.association {
        Association::Cell => components(f, cell)?,
        Association::Point => {
            let mut acc = vec![0.0; f.components.max(1)];
            for (n, w) in weights {
                for (a, v) in acc.iter_mut().zip(components(f, n)?) {
                    *a += w * v;
                }
            }
            acc
        }
    };
    Some((cell, value))
}

fn vector(v: &[f64]) -> Option<[f64; 3]> {
    match v {
        [x, y, z] => Some([*x, *y, *z]),
        [x, y] => Some([*x, *y, 0.0]),
        _ => None,
    }
}

/// Typical edge length of a cell.
fn cell_size(index: &CellIndex, cell: usize) -> f64 {
    let dim = index.dimensions.get(cell).copied().flatten().unwrap_or(3).max(1) as f64;
    index.measures.get(cell).map_or(0.0, |m| m.abs().powf(1.0 / dim))
}

// ---------------------------------------------------------------- trajectories

/// A particle carried by the flow with Stokes drag (and gravity), instead of a fluid tracer.
#[derive(Clone, Copy, Debug, PartialEq, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Particle {
    pub diameter_m: f64,
    pub density_kg_per_m3: f64,
    pub fluid_viscosity_pa_s: f64,
    #[serde(default)]
    pub gravity_m_per_s2: [f64; 3],
}

impl Particle {
    /// Momentum response time ρ d² / (18 μ).
    pub fn response_time(&self) -> f64 {
        self.density_kg_per_m3 * self.diameter_m * self.diameter_m / (18.0 * self.fluid_viscosity_pa_s)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TrajectoryOptions {
    /// Steps in each direction (default 2000, at most [`MAX_STEPS`]).
    #[serde(default)]
    pub max_steps: Option<usize>,
    /// Step length as a fraction of the local cell size (default 0.25).
    #[serde(default)]
    pub step_cells: Option<f64>,
    /// Also trace upstream from each seed (tracers only).
    #[serde(default)]
    pub both_directions: bool,
    #[serde(default)]
    pub particle: Option<Particle>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum TrajectoryEnd {
    /// The line reached the mesh boundary (an opening, a wall or an obstacle).
    LeftMesh,
    /// The flow speed fell to nothing.
    Stagnated,
    MaxSteps,
    /// The seed is not inside the mesh.
    SeedOutside,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Trajectory {
    pub seed: [f64; 3],
    pub points: Vec<[f64; 3]>,
    /// Time along the line (negative upstream of the seed).
    pub times: Vec<f64>,
    /// Flow speed (tracers) or particle speed at each point.
    pub speeds: Vec<f64>,
    pub end: TrajectoryEnd,
    /// Total time from the first to the last point.
    pub residence_time: f64,
}

/// Where trajectories start.
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Seeds {
    Points {
        points: Vec<[f64; 3]>,
    },
    /// `count` seeds evenly spaced from `start` to `end`.
    Line {
        start: [f64; 3],
        end: [f64; 3],
        count: usize,
    },
    /// A `nu × nv` grid on the rectangle `origin + s u + t v` (0 ≤ s, t ≤ 1).
    Rectangle {
        origin: [f64; 3],
        u: [f64; 3],
        v: [f64; 3],
        nu: usize,
        nv: usize,
    },
}

impl Seeds {
    pub fn points(&self) -> Vec<[f64; 3]> {
        match self {
            Self::Points { points } => points.clone(),
            Self::Line { start, end, count } => seeds_on_line(*start, *end, *count),
            Self::Rectangle { origin, u, v, nu, nv } => seeds_on_rectangle(*origin, *u, *v, *nu, *nv),
        }
    }
}

/// `n` seeds evenly spaced from `a` to `b` (both included; one seed is the midpoint).
pub fn seeds_on_line(a: [f64; 3], b: [f64; 3], n: usize) -> Vec<[f64; 3]> {
    let n = n.clamp(1, MAX_SEEDS);
    (0..n).map(|i| if n == 1 { 0.5 } else { i as f64 / (n - 1) as f64 }).map(|t| add(a, scale(sub(b, a), t))).collect()
}

/// Seeds at cell-centre positions of an `nu × nv` grid on the rectangle `origin + s u + t v`.
pub fn seeds_on_rectangle(origin: [f64; 3], u: [f64; 3], v: [f64; 3], nu: usize, nv: usize) -> Vec<[f64; 3]> {
    let (nu, nv) = (nu.max(1), nv.max(1));
    let mut out = Vec::new();
    for j in 0..nv {
        for i in 0..nu {
            if out.len() >= MAX_SEEDS {
                return out;
            }
            let (s, t) = ((i as f64 + 0.5) / nu as f64, (j as f64 + 0.5) / nv as f64);
            out.push(add(origin, add(scale(u, s), scale(v, t))));
        }
    }
    out
}

struct Tracer<'a> {
    grid: &'a Grid,
    index: &'a CellIndex,
    locator: &'a Locator,
    velocity: &'a Field,
    stagnant: f64,
}

impl Tracer<'_> {
    fn velocity(&self, x: [f64; 3]) -> Option<(usize, [f64; 3])> {
        let (cell, v) = sample(self.grid, self.index, self.locator, self.velocity, x)?;
        Some((cell, vector(&v)?))
    }

    /// One fluid-tracer pass from `seed` in direction `sign` (+1 downstream, −1 upstream).
    fn streamline(&self, seed: [f64; 3], sign: f64, steps: usize, step_cells: f64, budget: &mut usize) -> (Vec<([f64; 3], f64, f64)>, TrajectoryEnd) {
        let mut out = Vec::new();
        let Some((mut cell, u0)) = self.velocity(seed) else { return (out, TrajectoryEnd::SeedOutside) };
        let (mut x, mut t) = (seed, 0.0);
        out.push((x, t, norm(u0)));
        for _ in 0..steps {
            if *budget == 0 {
                return (out, TrajectoryEnd::MaxSteps);
            }
            let Some((_, k1)) = self.velocity(x) else { return (out, TrajectoryEnd::LeftMesh) };
            let speed = norm(k1);
            if speed <= self.stagnant {
                return (out, TrajectoryEnd::Stagnated);
            }
            let dt = sign * step_cells * cell_size(self.index, cell) / speed;
            // Classic Runge–Kutta; a stage outside the mesh ends the line at the boundary.
            let stage = |p: [f64; 3]| self.velocity(p).map(|(_, v)| v);
            let next = (|| {
                let k2 = stage(add(x, scale(k1, dt / 2.0)))?;
                let k3 = stage(add(x, scale(k2, dt / 2.0)))?;
                let k4 = stage(add(x, scale(k3, dt)))?;
                let slope = scale(add(add(k1, scale(k2, 2.0)), add(scale(k3, 2.0), k4)), 1.0 / 6.0);
                Some(add(x, scale(slope, dt)))
            })();
            let Some(next) = next else { return (out, TrajectoryEnd::LeftMesh) };
            let Some((next_cell, u)) = self.velocity(next) else { return (out, TrajectoryEnd::LeftMesh) };
            (x, t, cell) = (next, t + dt, next_cell);
            out.push((x, t, norm(u)));
            *budget -= 1;
        }
        (out, TrajectoryEnd::MaxSteps)
    }

    /// A particle with drag response time τ, released at the local flow velocity. Each step
    /// integrates the drag law exactly for the flow velocity at the step's start.
    fn particle(
        &self,
        seed: [f64; 3],
        particle: &Particle,
        steps: usize,
        step_cells: f64,
        budget: &mut usize,
    ) -> (Vec<([f64; 3], f64, f64)>, TrajectoryEnd) {
        let mut out = Vec::new();
        let Some((mut cell, mut v)) = self.velocity(seed) else { return (out, TrajectoryEnd::SeedOutside) };
        let tau = particle.response_time();
        let drift = scale(particle.gravity_m_per_s2, tau);
        let (mut x, mut t) = (seed, 0.0);
        out.push((x, t, norm(v)));
        for _ in 0..steps {
            if *budget == 0 {
                return (out, TrajectoryEnd::MaxSteps);
            }
            let Some((_, u)) = self.velocity(x) else { return (out, TrajectoryEnd::LeftMesh) };
            let speed = norm(v).max(norm(u));
            if speed <= self.stagnant && norm(drift) <= self.stagnant {
                return (out, TrajectoryEnd::Stagnated);
            }
            let dt = step_cells * cell_size(self.index, cell) / speed.max(norm(add(u, drift))).max(f64::MIN_POSITIVE);
            let terminal = add(u, drift);
            let decay = (-dt / tau).exp();
            let lag = sub(v, terminal);
            let moved = add(scale(terminal, dt), scale(lag, tau * (1.0 - decay)));
            let next = add(x, moved);
            v = add(terminal, scale(lag, decay));
            let Some((next_cell, _)) = self.velocity(next) else { return (out, TrajectoryEnd::LeftMesh) };
            (x, t, cell) = (next, t + dt, next_cell);
            out.push((x, t, norm(v)));
            *budget -= 1;
        }
        (out, TrajectoryEnd::MaxSteps)
    }
}

/// Trace flow trajectories of the vector field `velocity_field` from `seeds`, or particle
/// paths when `options.particle` is set.
pub fn trajectories(
    grid: &Grid,
    index: &CellIndex,
    velocity_field: &str,
    seeds: &[[f64; 3]],
    options: &TrajectoryOptions,
) -> Result<Vec<Trajectory>, Error> {
    let f = field(grid, velocity_field)?;
    check_field(grid, f)?;
    if !(2..=3).contains(&f.components) {
        return Err(bad(format!("{velocity_field} is not a vector field")));
    }
    if seeds.is_empty() || seeds.len() > MAX_SEEDS {
        return Err(bad(format!("use 1 to {MAX_SEEDS} seeds")));
    }
    if seeds.iter().flatten().any(|v| !v.is_finite()) {
        return Err(bad("seeds must be finite"));
    }
    let steps = options.max_steps.unwrap_or(2000);
    if !(1..=MAX_STEPS).contains(&steps) {
        return Err(bad(format!("use 1 to {MAX_STEPS} steps")));
    }
    let step_cells = options.step_cells.unwrap_or(0.25);
    if !(step_cells.is_finite() && step_cells > 0.0 && step_cells <= 2.0) {
        return Err(bad("step_cells must be in (0, 2]"));
    }
    if let Some(p) = &options.particle {
        let tau = p.response_time();
        if !(tau.is_finite() && tau > 0.0) || !p.gravity_m_per_s2.iter().all(|g| g.is_finite()) {
            return Err(bad("particle diameter, density and viscosity must be positive and finite"));
        }
        if options.both_directions {
            return Err(bad("particles are traced downstream only"));
        }
    }
    let largest = (0..f.values.len() / f.components).filter_map(|i| components(f, i)).map(|v| magnitude(&v)).fold(0.0, f64::max);
    let locator = Locator::new(index);
    let tracer = Tracer { grid, index, locator: &locator, velocity: f, stagnant: (largest * 1e-9).max(1e-300) };
    let mut budget = MAX_TRAJECTORY_POINTS;
    let mut out = Vec::with_capacity(seeds.len());
    for &seed in seeds {
        let (mut points, end) = match &options.particle {
            Some(p) => tracer.particle(seed, p, steps, step_cells, &mut budget),
            None => tracer.streamline(seed, 1.0, steps, step_cells, &mut budget),
        };
        if options.both_directions && options.particle.is_none() && end != TrajectoryEnd::SeedOutside {
            let (mut back, _) = tracer.streamline(seed, -1.0, steps, step_cells, &mut budget);
            back.reverse();
            back.pop();
            back.append(&mut points);
            points = back;
        }
        let residence_time = match (points.first(), points.last()) {
            (Some(a), Some(b)) => b.1 - a.1,
            _ => 0.0,
        };
        out.push(Trajectory {
            seed,
            points: points.iter().map(|p| p.0).collect(),
            times: points.iter().map(|p| p.1).collect(),
            speeds: points.iter().map(|p| p.2).collect(),
            end,
            residence_time,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------- point values

/// One scalar per grid point: point data as stored (vector magnitude), cell data averaged
/// over the cells that use the point, weighted by cell measure.
fn point_scalars(grid: &Grid, index: &CellIndex, f: &Field) -> Vec<f64> {
    match f.association {
        Association::Point => (0..grid.points.len()).map(|p| components(f, p).map_or(f64::NAN, |v| magnitude(&v))).collect(),
        Association::Cell => {
            let mut sum = vec![0.0; grid.points.len()];
            let mut weight = vec![0.0; grid.points.len()];
            for c in 0..index.len() {
                let Some(v) = components(f, c).map(|v| magnitude(&v)).filter(|v| v.is_finite()) else { continue };
                let w = index.measures.get(c).copied().unwrap_or(0.0).abs().max(f64::MIN_POSITIVE);
                for &n in index.nodes.get(c).into_iter().flatten() {
                    if let (Some(s), Some(m)) = (sum.get_mut(n), weight.get_mut(n)) {
                        *s += v * w;
                        *m += w;
                    }
                }
            }
            sum.into_iter().zip(weight).map(|(s, w)| if w > 0.0 { s / w } else { f64::NAN }).collect()
        }
    }
}

// ---------------------------------------------------------------- isosurfaces

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Isosurface {
    pub field: String,
    pub level: f64,
    pub triangles: Vec<[[f64; 3]; 3]>,
    /// `color_field` at each triangle corner (scalar or magnitude), when one was asked for.
    pub colors: Vec<[f64; 3]>,
    pub area: f64,
    /// The triangle cap was reached and the surface is incomplete.
    pub truncated: bool,
}

/// The surface where `field_name` (scalar or magnitude) equals `level`, by marching
/// tetrahedra over the 3-D cells' tetrahedral split. `color_field`, if given, is
/// interpolated onto the surface for colouring by a second quantity.
pub fn isosurface(grid: &Grid, index: &CellIndex, field_name: &str, level: f64, color_field: Option<&str>) -> Result<Isosurface, Error> {
    let f = field(grid, field_name)?;
    check_field(grid, f)?;
    if !level.is_finite() {
        return Err(bad("the iso level must be finite"));
    }
    if index.dimension != 3 {
        return Err(bad("isosurfaces need a 3-D mesh"));
    }
    let values = point_scalars(grid, index, f);
    let colors = match color_field {
        Some(name) => {
            let c = field(grid, name)?;
            check_field(grid, c)?;
            Some(point_scalars(grid, index, c))
        }
        None => None,
    };
    let mut out = Isosurface { field: field_name.into(), level, triangles: Vec::new(), colors: Vec::new(), area: 0.0, truncated: false };
    let emit = |corners: [(usize, usize); 3], out: &mut Isosurface| {
        let cut = |(a, b): (usize, usize)| -> Option<([f64; 3], f64)> {
            let (va, vb) = (*values.get(a)?, *values.get(b)?);
            let t = if (vb - va).abs() > 0.0 { ((level - va) / (vb - va)).clamp(0.0, 1.0) } else { 0.5 };
            let (pa, pb) = (*grid.points.get(a)?, *grid.points.get(b)?);
            let color = colors.as_ref().map_or(0.0, |c| {
                let (ca, cb) = (c.get(a).copied().unwrap_or(f64::NAN), c.get(b).copied().unwrap_or(f64::NAN));
                ca + t * (cb - ca)
            });
            Some((add(pa, scale(sub(pb, pa), t)), color))
        };
        let (Some(p), Some(q), Some(r)) = (cut(corners[0]), cut(corners[1]), cut(corners[2])) else { return };
        if out.triangles.len() >= MAX_ISO_TRIANGLES {
            out.truncated = true;
            return;
        }
        out.area += 0.5 * norm(cross(sub(q.0, p.0), sub(r.0, p.0)));
        out.triangles.push([p.0, q.0, r.0]);
        if colors.is_some() {
            out.colors.push([p.1, q.1, r.1]);
        }
    };
    for c in 0..index.len() {
        if index.dimensions.get(c).copied().flatten() != Some(3) {
            continue;
        }
        for tet in index.simplices.get(c).into_iter().flatten() {
            let [a, b, cc, d] = tet.as_slice() else { continue };
            let nodes = [*a, *b, *cc, *d];
            let v: Option<Vec<f64>> = nodes.iter().map(|&n| values.get(n).copied().filter(|v| v.is_finite())).collect();
            let Some(v) = v else { continue };
            let above: Vec<usize> = (0..4).filter(|&i| v[i] > level).collect();
            let below: Vec<usize> = (0..4).filter(|&i| v[i] <= level).collect();
            match (above.as_slice(), below.as_slice()) {
                ([i], [j, k, l]) | ([j, k, l], [i]) => {
                    let n = |x: usize| nodes[x];
                    emit([(n(*i), n(*j)), (n(*i), n(*k)), (n(*i), n(*l))], &mut out);
                }
                ([i, j], [k, l]) => {
                    let n = |x: usize| nodes[x];
                    // The cut is a quadrilateral i-k, i-l, j-l, j-k: two triangles.
                    emit([(n(*i), n(*k)), (n(*i), n(*l)), (n(*j), n(*l))], &mut out);
                    emit([(n(*i), n(*k)), (n(*j), n(*l)), (n(*j), n(*k))], &mut out);
                }
                _ => {}
            }
            if out.truncated {
                return Ok(out);
            }
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------- boundary surfaces

/// Faces of a 3-D cell in VTK node order (indices into the cell's node list).
fn cell_faces(kind: u8) -> &'static [&'static [usize]] {
    match kind {
        10 => &[&[0, 1, 2], &[0, 1, 3], &[1, 2, 3], &[0, 2, 3]],
        12 => &[&[0, 1, 2, 3], &[4, 5, 6, 7], &[0, 1, 5, 4], &[1, 2, 6, 5], &[2, 3, 7, 6], &[3, 0, 4, 7]],
        13 => &[&[0, 1, 2], &[3, 4, 5], &[0, 1, 4, 3], &[1, 2, 5, 4], &[2, 0, 3, 5]],
        14 => &[&[0, 1, 2, 3], &[0, 1, 4], &[1, 2, 4], &[2, 3, 4], &[3, 0, 4]],
        _ => &[],
    }
}

/// A plane of the mesh's bounding box.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoxSide {
    /// 0, 1 or 2 for x, y or z.
    pub axis: usize,
    /// The maximum side when true, else the minimum side.
    pub max: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct SurfaceFace {
    pub corners: Vec<[f64; 3]>,
    /// Outward unit normal.
    pub normal: [f64; 3],
    pub area: f64,
    /// Scalar value or vector magnitude on the face.
    pub value: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Surface {
    pub field: String,
    pub faces: Vec<SurfaceFace>,
    pub area: f64,
    /// Area-weighted mean, minimum, maximum and integral of the value.
    pub mean: f64,
    pub min: f64,
    pub max: f64,
    pub integral: f64,
    /// For vector fields, ∫ v·n dA with outward n (a volume flow for velocity).
    pub flux: Option<f64>,
    /// The face cap was reached; aggregates still cover every face.
    pub truncated: bool,
}

/// The mesh boundary (faces used by one 3-D cell), optionally only the part lying on one
/// side of the bounding box, with `field_name` on each face and area-weighted statistics.
/// Point data use the mean of a face's corners; cell data the owning cell's value.
pub fn surface(grid: &Grid, index: &CellIndex, field_name: &str, side: Option<BoxSide>) -> Result<Surface, Error> {
    let f = field(grid, field_name)?;
    check_field(grid, f)?;
    if index.dimension != 3 {
        return Err(bad("surfaces need a 3-D mesh"));
    }
    if side.is_some_and(|s| s.axis > 2) {
        return Err(bad("axis must be 0, 1 or 2"));
    }
    let mut seen: std::collections::HashMap<Vec<usize>, (usize, usize, u32)> = std::collections::HashMap::new();
    for c in 0..index.len() {
        let (Some(&kind), Some(nodes)) = (grid.cell_types.get(c), index.nodes.get(c)) else { continue };
        if index.dimensions.get(c).copied().flatten() != Some(3) {
            continue;
        }
        for (k, face) in cell_faces(kind).iter().enumerate() {
            let mut key: Vec<usize> = face.iter().filter_map(|&i| nodes.get(i).copied()).collect();
            key.sort_unstable();
            seen.entry(key).and_modify(|e| e.2 += 1).or_insert((c, k, 1));
        }
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in &grid.points {
        for a in 0..3 {
            lo[a] = lo[a].min(p[a]);
            hi[a] = hi[a].max(p[a]);
        }
    }
    let tolerance = 1e-9 * (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max).max(f64::MIN_POSITIVE);
    let mut boundary: Vec<(usize, usize)> = seen.into_values().filter(|e| e.2 == 1).map(|e| (e.0, e.1)).collect();
    boundary.sort_unstable();
    let mut out = Surface {
        field: field_name.into(),
        faces: Vec::new(),
        area: 0.0,
        mean: 0.0,
        min: f64::INFINITY,
        max: f64::NEG_INFINITY,
        integral: 0.0,
        flux: None,
        truncated: false,
    };
    let vectorial = (2..=3).contains(&f.components);
    let mut flux = 0.0;
    for (c, k) in boundary {
        let (Some(&kind), Some(nodes)) = (grid.cell_types.get(c), index.nodes.get(c)) else { continue };
        let Some(face) = cell_faces(kind).get(k) else { continue };
        let ids: Vec<usize> = face.iter().filter_map(|&i| nodes.get(i).copied()).collect();
        let corners: Option<Vec<[f64; 3]>> = ids.iter().map(|&n| grid.points.get(n).copied()).collect();
        let Some(corners) = corners else { continue };
        if let Some(s) = side {
            let plane = if s.max { hi[s.axis] } else { lo[s.axis] };
            if corners.iter().any(|p| (p[s.axis] - plane).abs() > tolerance) {
                continue;
            }
        }
        // Area vector by fan triangulation; orient it away from the cell centre.
        let mut area_vector = [0.0; 3];
        for i in 1..corners.len().saturating_sub(1) {
            area_vector = add(area_vector, scale(cross(sub(corners[i], corners[0]), sub(corners[i + 1], corners[0])), 0.5));
        }
        let area = norm(area_vector);
        if area.is_nan() || area <= 0.0 {
            continue;
        }
        let centre = scale(corners.iter().fold([0.0; 3], |a, p| add(a, *p)), 1.0 / corners.len() as f64);
        let inward = sub(index.centroids.get(c).copied().unwrap_or(centre), centre);
        let mut normal = scale(area_vector, 1.0 / area);
        if dot(normal, inward) > 0.0 {
            normal = scale(normal, -1.0);
        }
        let raw: Option<Vec<f64>> = match f.association {
            Association::Cell => components(f, c),
            Association::Point => {
                let mut acc = vec![0.0; f.components.max(1)];
                let mut ok = true;
                for &n in &ids {
                    match components(f, n) {
                        Some(v) => acc.iter_mut().zip(v).for_each(|(a, v)| *a += v / ids.len() as f64),
                        None => ok = false,
                    }
                }
                ok.then_some(acc)
            }
        };
        let Some(raw) = raw.filter(|v| v.iter().all(|x| x.is_finite())) else { continue };
        let value = magnitude(&raw);
        if vectorial && let Some(v) = vector(&raw) {
            flux += dot(v, normal) * area;
        }
        out.area += area;
        out.integral += value * area;
        out.min = out.min.min(value);
        out.max = out.max.max(value);
        if out.faces.len() < MAX_SURFACE_FACES {
            out.faces.push(SurfaceFace { corners, normal, area, value });
        } else {
            out.truncated = true;
        }
    }
    if out.area <= 0.0 {
        return Err(bad("no boundary faces lie on that side"));
    }
    out.mean = out.integral / out.area;
    out.flux = vectorial.then_some(flux);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An n × n × n hexahedral mesh of the unit cube with point fields from `point`.
    type PointField<'a> = (&'a str, usize, fn([f64; 3]) -> Vec<f64>);

    fn cube(n: usize, point: &[PointField<'_>]) -> Grid {
        let id = |i: usize, j: usize, k: usize| (k * (n + 1) + j) * (n + 1) + i;
        let h = 1.0 / n as f64;
        let mut points = Vec::new();
        for k in 0..=n {
            for j in 0..=n {
                for i in 0..=n {
                    points.push([i as f64 * h, j as f64 * h, k as f64 * h]);
                }
            }
        }
        let mut connectivity = Vec::new();
        let mut offsets = Vec::new();
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    connectivity.extend([
                        id(i, j, k),
                        id(i + 1, j, k),
                        id(i + 1, j + 1, k),
                        id(i, j + 1, k),
                        id(i, j, k + 1),
                        id(i + 1, j, k + 1),
                        id(i + 1, j + 1, k + 1),
                        id(i, j + 1, k + 1),
                    ]);
                    offsets.push(connectivity.len());
                }
            }
        }
        let cells = offsets.len();
        let fields = point
            .iter()
            .map(|(name, k, f)| Field {
                name: (*name).into(),
                association: Association::Point,
                data_type: "Float64".into(),
                components: *k,
                values: points.iter().flat_map(|p| f(*p)).collect(),
            })
            .collect();
        Grid { points, connectivity, offsets, cell_types: vec![12; cells], fields }
    }

    #[test]
    fn locator_finds_the_cell_of_any_inside_point() {
        let grid = cube(6, &[("x", 1, |p| vec![p[0]])]);
        let index = CellIndex::new(&grid).expect("index");
        let locator = Locator::new(&index);
        for x in [[0.01, 0.5, 0.99], [0.5, 0.5, 0.5], [0.999, 0.001, 0.3], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]] {
            let (cell, weights) = locator.locate(&grid, &index, x).expect("inside");
            let b = index.bounds[cell];
            assert!((0..3).all(|a| x[a] >= b[a][0] - 1e-12 && x[a] <= b[a][1] + 1e-12));
            assert!((weights.iter().map(|w| w.1).sum::<f64>() - 1.0).abs() < 1e-12);
            let (_, v) = sample(&grid, &index, &locator, &grid.fields[0], x).expect("value");
            assert!((v[0] - x[0]).abs() < 1e-12, "linear fields interpolate exactly");
        }
        assert!(locator.locate(&grid, &index, [1.1, 0.5, 0.5]).is_none());
        assert!(locator.locate(&grid, &index, [f64::NAN, 0.5, 0.5]).is_none());
    }

    #[test]
    fn streamlines_follow_uniform_and_rotating_flows() {
        // Uniform flow along x: lines are straight, take 1/u to cross, and leave the mesh.
        let grid = cube(5, &[("velocity", 3, |_| vec![2.0, 0.0, 0.0])]);
        let index = CellIndex::new(&grid).expect("index");
        let seeds = seeds_on_line([0.0, 0.2, 0.5], [0.0, 0.8, 0.5], 4);
        let lines = trajectories(
            &grid,
            &index,
            "velocity",
            &seeds,
            &TrajectoryOptions { max_steps: None, step_cells: None, both_directions: false, particle: None },
        )
        .expect("trace");
        assert_eq!(lines.len(), 4);
        for line in &lines {
            assert_eq!(line.end, TrajectoryEnd::LeftMesh);
            assert!(line.points.iter().all(|p| (p[1] - line.seed[1]).abs() < 1e-12 && (p[2] - 0.5).abs() < 1e-12));
            let last = line.points.last().copied().unwrap_or_default();
            assert!(last[0] > 0.9, "reaches the outlet side: {last:?}");
            assert!((line.residence_time - last[0] / 2.0).abs() < 1e-9);
        }
        // Solid-body rotation about the cube's z axis: RK4 keeps the radius.
        let grid = cube(8, &[("velocity", 3, |p| vec![-(p[1] - 0.5), p[0] - 0.5, 0.0])]);
        let index = CellIndex::new(&grid).expect("index");
        let options = TrajectoryOptions { max_steps: Some(400), step_cells: Some(0.5), both_directions: false, particle: None };
        let line = trajectories(&grid, &index, "velocity", &[[0.8, 0.5, 0.5]], &options).expect("trace").remove(0);
        assert_eq!(line.end, TrajectoryEnd::MaxSteps);
        let radius = |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt();
        assert!(line.points.iter().all(|p| (radius(*p) - 0.3).abs() < 1e-3), "radius drift");
        // Upstream and downstream together; the seed sits at time 0.
        let grid = cube(5, &[("velocity", 3, |_| vec![1.0, 0.0, 0.0])]);
        let index = CellIndex::new(&grid).expect("index");
        let both = TrajectoryOptions { both_directions: true, ..options.clone() };
        let line = trajectories(&grid, &index, "velocity", &[[0.5, 0.5, 0.5]], &both).expect("trace").remove(0);
        assert!(line.points.first().is_some_and(|p| p[0] < 0.1) && line.points.last().is_some_and(|p| p[0] > 0.9));
        assert!(line.times.windows(2).all(|w| w[0] < w[1]));
        assert!(line.times.contains(&0.0));
        // Seeds outside, scalar fields and empty seed lists are reported or refused.
        let outside = trajectories(&grid, &index, "velocity", &[[2.0, 0.5, 0.5]], &options).expect("trace");
        assert_eq!(outside[0].end, TrajectoryEnd::SeedOutside);
        let scalar = cube(2, &[("t", 1, |p| vec![p[0]])]);
        let scalar_index = CellIndex::new(&scalar).expect("index");
        assert!(trajectories(&scalar, &scalar_index, "t", &[[0.5; 3]], &options).is_err());
        assert!(trajectories(&grid, &index, "velocity", &[], &options).is_err());
    }

    #[test]
    fn particles_settle_at_their_stokes_terminal_velocity() {
        // Still air, gravity down y: speed approaches τ g; a 50 µm water droplet in air.
        let grid = cube(4, &[("velocity", 3, |_| vec![0.0, 0.0, 0.0])]);
        let index = CellIndex::new(&grid).expect("index");
        let particle = Particle { diameter_m: 50e-6, density_kg_per_m3: 1000.0, fluid_viscosity_pa_s: 1.8e-5, gravity_m_per_s2: [0.0, -9.81, 0.0] };
        let tau = particle.response_time();
        assert!((tau - 1000.0 * 2.5e-9 / (18.0 * 1.8e-5)).abs() < 1e-15);
        let options = TrajectoryOptions { max_steps: Some(5000), step_cells: Some(0.05), both_directions: false, particle: Some(particle) };
        let path = trajectories(&grid, &index, "velocity", &[[0.5, 0.95, 0.5]], &options).expect("trace").remove(0);
        assert_eq!(path.end, TrajectoryEnd::LeftMesh, "it falls out through the floor");
        let terminal = tau * 9.81;
        let last = path.speeds.last().copied().unwrap_or(0.0);
        assert!((last - terminal).abs() < 1e-3 * terminal, "{last} vs {terminal}");
        assert!(path.points.iter().all(|p| (p[0] - 0.5).abs() < 1e-12));
        // A heavy particle in uniform flow keeps up with the flow (released at flow speed).
        let flow = cube(4, &[("velocity", 3, |_| vec![1.0, 0.0, 0.0])]);
        let flow_index = CellIndex::new(&flow).expect("index");
        let no_gravity = TrajectoryOptions { particle: Some(Particle { gravity_m_per_s2: [0.0; 3], ..particle }), ..options.clone() };
        let path = trajectories(&flow, &flow_index, "velocity", &[[0.1, 0.5, 0.5]], &no_gravity).expect("trace").remove(0);
        assert!(path.speeds.iter().all(|s| (s - 1.0).abs() < 1e-12));
        assert!(
            trajectories(&flow, &flow_index, "velocity", &[[0.1, 0.5, 0.5]], &TrajectoryOptions { both_directions: true, ..no_gravity }).is_err()
        );
    }

    #[test]
    fn isosurfaces_of_planes_and_spheres_have_the_right_area() {
        let grid = cube(6, &[("t", 1, |p| vec![p[0] + 0.0 * p[1]]), ("c", 1, |p| vec![p[1]])]);
        let index = CellIndex::new(&grid).expect("index");
        let plane = isosurface(&grid, &index, "t", 0.42, Some("c")).expect("iso");
        assert!((plane.area - 1.0).abs() < 1e-9, "a cut through the cube: {}", plane.area);
        assert!(plane.triangles.iter().flatten().all(|p| (p[0] - 0.42).abs() < 1e-12));
        assert_eq!(plane.colors.len(), plane.triangles.len());
        for (t, c) in plane.triangles.iter().zip(&plane.colors) {
            for (p, v) in t.iter().zip(c) {
                assert!((p[1] - v).abs() < 1e-12, "colour is the second field there");
            }
        }
        // Distance from the centre: an approximate sphere of radius 0.3 (area 4πr²).
        let grid = cube(24, &[("r", 1, |p| vec![((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2) + (p[2] - 0.5).powi(2)).sqrt()])]);
        let index = CellIndex::new(&grid).expect("index");
        let sphere = isosurface(&grid, &index, "r", 0.3, None).expect("iso");
        let exact = 4.0 * std::f64::consts::PI * 0.09;
        assert!((sphere.area - exact).abs() < 0.03 * exact, "{} vs {exact}", sphere.area);
        assert!(isosurface(&grid, &index, "r", 5.0, None).expect("iso").triangles.is_empty());
        assert!(isosurface(&grid, &index, "r", f64::NAN, None).is_err());
    }

    #[test]
    fn surfaces_give_areas_means_and_outward_fluxes() {
        let grid = cube(4, &[("velocity", 3, |_| vec![3.0, 0.0, 0.0]), ("t", 1, |p| vec![10.0 * p[2]])]);
        let index = CellIndex::new(&grid).expect("index");
        let all = surface(&grid, &index, "velocity", None).expect("surface");
        assert!((all.area - 6.0).abs() < 1e-12);
        assert_eq!(all.faces.len(), 6 * 16);
        assert!(all.flux.is_some_and(|f| f.abs() < 1e-12), "what flows in flows out");
        let outlet = surface(&grid, &index, "velocity", Some(BoxSide { axis: 0, max: true })).expect("outlet");
        assert!((outlet.area - 1.0).abs() < 1e-12);
        assert!(outlet.flux.is_some_and(|f| (f - 3.0).abs() < 1e-12), "volume flow = u A");
        assert!(outlet.faces.iter().all(|f| (f.normal[0] - 1.0).abs() < 1e-12));
        let inlet = surface(&grid, &index, "velocity", Some(BoxSide { axis: 0, max: false })).expect("inlet");
        assert!(inlet.flux.is_some_and(|f| (f + 3.0).abs() < 1e-12));
        // A scalar: the top is at t = 10, the side faces average 5.
        let top = surface(&grid, &index, "t", Some(BoxSide { axis: 2, max: true })).expect("top");
        assert!((top.mean - 10.0).abs() < 1e-12 && top.flux.is_none());
        let side = surface(&grid, &index, "t", Some(BoxSide { axis: 1, max: false })).expect("side");
        assert!((side.mean - 5.0).abs() < 1e-12 && (side.min - 1.25).abs() < 1e-12 && (side.max - 8.75).abs() < 1e-12);
        assert!(surface(&grid, &index, "t", Some(BoxSide { axis: 3, max: false })).is_err());
    }
}
