//! Lens interchange (ZMX/SEQ/LEN subset), JSON serialization, a non-sequential plane scene with
//! stray-light and ghost ranking, and optical optimization (coordinate descent, seeded
//! differential evolution, sensitivities, Monte-Carlo tolerancing, forward-mode derivative),
//! ported from upstream `cfd/optics/{interchange, serialization, nonsequential, optimization}`
//! (MIT, cfd_solvers contributors).
use super::{reflect, refract, Ray, SequentialOpticalSystem, SequentialSurface, SurfaceType};
use crate::{
    core::{parse_float_prefix, Dual, Mt19937_64, StdNormal, Vec3},
    lbm::format_general,
    require, Result, SolverError,
};

fn g17(v: f64) -> String {
    format_general(v, 17)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LensInterchangeFormat {
    ZemaxZmx,
    CodevSeq,
    OsloLen,
}

fn type_word(t: SurfaceType, f: LensInterchangeFormat) -> &'static str {
    match t {
        SurfaceType::Plane if f == LensInterchangeFormat::CodevSeq => "PLANO",
        SurfaceType::Plane => "PLANE",
        SurfaceType::Sphere => "STANDARD",
        SurfaceType::Conic => "CONIC",
        SurfaceType::EvenAsphere => "EVENASPH",
    }
}

/// Text export of the common sequential subset (type, curvature/radius, spacing, aperture,
/// index, conic, four even asphere terms).
pub fn export_lens_system(
    system: &SequentialOpticalSystem,
    format: LensInterchangeFormat,
) -> String {
    use std::fmt::Write;
    let mut o = String::from("# cfd_solvers sequential lens interchange\n");
    let tag = match format {
        LensInterchangeFormat::ZemaxZmx => "ZMX",
        LensInterchangeFormat::CodevSeq => "SEQ",
        LensInterchangeFormat::OsloLen => "LEN",
    };
    let _ = writeln!(o, "# format {tag}");
    let mut previous_z = 0.0;
    for (index, s) in system.surfaces().iter().enumerate() {
        let thickness = s.vertex_z - previous_z;
        let plane = s.surface_type == SurfaceType::Plane;
        let curvature = if plane { 0.0 } else { 1.0 / s.radius };
        let radius = if plane { 0.0 } else { s.radius };
        let _ = writeln!(
            o,
            "SURF {index}\nTYPE {}",
            type_word(s.surface_type, format)
        );
        let (rk, rv, th, ap, gl, co, pa) = match format {
            LensInterchangeFormat::ZemaxZmx => {
                ("CURV", curvature, "DISZ", "DIAM", "GLAS", "CONI", "PARM")
            }
            LensInterchangeFormat::CodevSeq => (
                "RADIUS", radius, "THICK", "APERTURE", "GLASS", "CONIC", "ASPH",
            ),
            LensInterchangeFormat::OsloLen => ("RD", radius, "TH", "AP", "GLASS", "CONI", "COEF"),
        };
        let _ = write!(
            o,
            "{rk} {}\n{th} {}\n{ap} {}\n{gl} {}\n{co} {}\n",
            g17(rv),
            g17(thickness),
            g17(s.aperture_radius),
            g17(s.refractive_index_after),
            g17(s.conic_constant)
        );
        for (i, c) in s.even_coefficients.iter().enumerate() {
            let _ = writeln!(o, "{pa} {} {}", i + 1, g17(*c));
        }
        previous_z = s.vertex_z;
    }
    o
}

/// Whole-token `std::stod` (finite decimal values; out-of-range tokens are skipped as upstream
/// skips them on `std::out_of_range`).
fn whole_number(token: &str) -> Option<f64> {
    let lower = token.to_ascii_lowercase();
    let t = lower.trim_start_matches(['+', '-']);
    if t == "inf" || t == "infinity" || t == "nan" {
        return token.parse().ok();
    }
    match parse_float_prefix(token) {
        Some((v, used)) if used == token.len() && v.is_finite() => Some(v),
        _ => None,
    }
}

fn parse_last_number(tokens: &[&str]) -> Result<f64> {
    tokens
        .iter()
        .rev()
        .find_map(|t| whole_number(t))
        .ok_or_else(|| SolverError::InvalidInput("lens record is missing a numeric value".into()))
}

#[derive(Default)]
struct SurfaceRecord {
    active: bool,
    surface: SequentialSurface,
    thickness: f64,
    curvature: f64,
}

fn finish_surface(
    records: &mut Vec<SequentialSurface>,
    record: &mut SurfaceRecord,
    previous_z: &mut f64,
) {
    if !record.active {
        return;
    }
    let s = &mut record.surface;
    if s.surface_type != SurfaceType::Plane && record.curvature.abs() > 0.0 {
        s.radius = 1.0 / record.curvature;
    }
    if s.surface_type != SurfaceType::Plane && s.radius == 0.0 {
        s.surface_type = SurfaceType::Plane;
    }
    s.vertex_z = *previous_z + record.thickness;
    *previous_z = s.vertex_z;
    records.push(*s);
    *record = SurfaceRecord::default();
}

/// Reads any of the three dialects (keywords are recognised regardless of `format`); unknown
/// records are ignored and `#`/`!` start comments.
pub fn import_lens_system(
    text: &str,
    _format: LensInterchangeFormat,
) -> Result<SequentialOpticalSystem> {
    let mut records = Vec::new();
    let mut record = SurfaceRecord::default();
    let mut previous_z = 0.0;
    for raw in text.lines() {
        let line = raw.split(['#', '!']).next().unwrap_or("");
        let tokens: Vec<&str> = line.split_ascii_whitespace().collect();
        let Some(first) = tokens.first() else {
            continue;
        };
        let keyword = first.to_ascii_uppercase();
        if keyword == "SURF" || keyword == "SURFACE" {
            finish_surface(&mut records, &mut record, &mut previous_z);
            record.active = true;
            continue;
        }
        if !record.active {
            continue;
        }
        match keyword.as_str() {
            "TYPE" if tokens.len() >= 2 => {
                record.surface.surface_type = match tokens[1].to_ascii_uppercase().as_str() {
                    "PLANE" | "PLANO" => SurfaceType::Plane,
                    "CONIC" => SurfaceType::Conic,
                    "EVENASPH" | "EVEN_ASPHERE" => SurfaceType::EvenAsphere,
                    _ => SurfaceType::Sphere,
                }
            }
            "CURV" => record.curvature = parse_last_number(&tokens)?,
            "RADIUS" | "RD" => {
                let r = parse_last_number(&tokens)?;
                record.curvature = if r == 0.0 { 0.0 } else { 1.0 / r };
            }
            "DISZ" | "THICK" | "TH" => record.thickness = parse_last_number(&tokens)?,
            "DIAM" | "APERTURE" | "AP" => {
                record.surface.aperture_radius = parse_last_number(&tokens)?
            }
            "GLAS" | "GLASS" => record.surface.refractive_index_after = parse_last_number(&tokens)?,
            "CONI" => record.surface.conic_constant = parse_last_number(&tokens)?,
            "PARM" | "ASPH" | "COEF" if tokens.len() >= 3 => {
                let index = parse_last_number(&tokens[1..2])?;
                let value = parse_last_number(&tokens[2..3])?;
                if (1.0..=4.0).contains(&index) {
                    record.surface.even_coefficients[(index - 1.0) as usize] = value;
                }
            }
            _ => {}
        }
    }
    finish_surface(&mut records, &mut record, &mut previous_z);
    require(!records.is_empty(), "lens interchange contains no surfaces")?;
    let mut system = SequentialOpticalSystem::new(1.0)?;
    for s in records {
        system.add_surface(s)?;
    }
    Ok(system)
}

fn type_name(t: SurfaceType) -> &'static str {
    match t {
        SurfaceType::Plane => "plane",
        SurfaceType::Sphere => "sphere",
        SurfaceType::Conic => "conic",
        SurfaceType::EvenAsphere => "even_asphere",
    }
}

/// Compact JSON with 17 significant digits.
pub fn sequential_system_to_json(
    system: &SequentialOpticalSystem,
    object_space_index: f64,
) -> String {
    let mut o = format!(
        "{{\"object_index\":{},\"surfaces\":[",
        g17(object_space_index)
    );
    for (i, s) in system.surfaces().iter().enumerate() {
        if i > 0 {
            o.push(',');
        }
        let a: Vec<String> = s.even_coefficients.iter().map(|v| g17(*v)).collect();
        o += &format!(
            "{{\"type\":\"{}\",\"z\":{},\"radius\":{},\"aperture\":{},\"n_after\":{},\"k\":{},\"a\":[{}]}}",
            type_name(s.surface_type),
            g17(s.vertex_z),
            g17(s.radius),
            g17(s.aperture_radius),
            g17(s.refractive_index_after),
            g17(s.conic_constant),
            a.join(",")
        );
    }
    o + "]}"
}

fn is_num_char(c: u8) -> bool {
    c.is_ascii_digit() || matches!(c, b'-' | b'+' | b'.' | b'e' | b'E')
}

/// `[-+0-9.eE]+` at `pos`: returns the span end.
fn num_span(b: &[u8], pos: usize) -> Option<usize> {
    let mut e = pos;
    while e < b.len() && is_num_char(b[e]) {
        e += 1;
    }
    (e > pos).then_some(e)
}

fn stod(s: &str) -> Result<f64> {
    parse_float_prefix(s)
        .map(|(v, _)| v)
        .ok_or_else(|| SolverError::InvalidInput(format!("invalid optical JSON number: {s}")))
}

/// Matches one serialized surface record at `pos` (the exact layout written by
/// [`sequential_system_to_json`]); returns the captures and the end of the match.
fn match_surface(text: &str, pos: usize) -> Option<(Vec<&str>, usize)> {
    let b = text.as_bytes();
    let mut p = pos;
    let lit = |p: &mut usize, s: &str| -> Option<()> {
        if text.get(*p..)?.starts_with(s) {
            *p += s.len();
            Some(())
        } else {
            None
        }
    };
    let mut caps = Vec::with_capacity(11);
    lit(&mut p, "{\"type\":\"")?;
    let start = p;
    while p < b.len() && b[p] != b'"' {
        p += 1;
    }
    if p == start || p >= b.len() {
        return None;
    }
    caps.push(&text[start..p]);
    p += 1;
    for key in ["z", "radius", "aperture", "n_after", "k"] {
        lit(&mut p, &format!(",\"{key}\":"))?;
        let e = num_span(b, p)?;
        caps.push(&text[p..e]);
        p = e;
    }
    lit(&mut p, ",\"a\":[")?;
    for i in 0..4 {
        if i > 0 {
            lit(&mut p, ",")?;
        }
        let e = num_span(b, p)?;
        caps.push(&text[p..e]);
        p = e;
    }
    lit(&mut p, "]}")?;
    Some((caps, p))
}

/// Inverse of [`sequential_system_to_json`]. Like upstream's regex reader it scans for records
/// of the exact written layout and silently skips anything else.
pub fn sequential_system_from_json(json: &str) -> Result<SequentialOpticalSystem> {
    let key = "\"object_index\":";
    let mut object_index = None;
    let mut from = 0;
    while let Some(off) = json[from..].find(key) {
        let p = from + off + key.len();
        if let Some(e) = num_span(json.as_bytes(), p) {
            object_index = Some(stod(&json[p..e])?);
            break;
        }
        from = from + off + 1;
    }
    let object_index = object_index
        .ok_or_else(|| SolverError::InvalidInput("missing optical object index".into()))?;
    let mut system = SequentialOpticalSystem::new(object_index)?;
    let mut pos = 0;
    while let Some(off) = json[pos..].find("{\"type\":\"") {
        let at = pos + off;
        match match_surface(json, at) {
            Some((c, end)) => {
                let surface_type = match c[0] {
                    "plane" => SurfaceType::Plane,
                    "sphere" => SurfaceType::Sphere,
                    "conic" => SurfaceType::Conic,
                    "even_asphere" => SurfaceType::EvenAsphere,
                    _ => {
                        return Err(SolverError::InvalidInput(
                            "unknown optical surface type".into(),
                        ))
                    }
                };
                let mut s = SequentialSurface {
                    surface_type,
                    ..Default::default()
                };
                s.vertex_z = stod(c[1])?;
                s.radius = stod(c[2])?;
                s.aperture_radius = stod(c[3])?;
                s.refractive_index_after = stod(c[4])?;
                s.conic_constant = stod(c[5])?;
                for i in 0..4 {
                    s.even_coefficients[i] = stod(c[6 + i])?;
                }
                system.add_surface(s)?;
                pos = end;
            }
            None => pos = at + 1,
        }
    }
    Ok(system)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsPlane {
    pub point: Vec3,
    pub normal: Vec3,
    pub refractive_index_before: f64,
    pub refractive_index_after: f64,
    pub reflectance: f64,
    pub transmittance: f64,
    pub absorptance: f64,
    pub scatter_fraction: f64,
}

impl Default for NsPlane {
    fn default() -> Self {
        Self {
            point: Vec3::ZERO,
            normal: Vec3::new(0.0, 0.0, 1.0),
            refractive_index_before: 1.0,
            refractive_index_after: 1.5,
            reflectance: 0.0,
            transmittance: 1.0,
            absorptance: 0.0,
            scatter_fraction: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsRay {
    pub ray: Ray,
    pub power: f64,
    pub bounces: usize,
}

impl Default for NsRay {
    fn default() -> Self {
        Self {
            ray: Ray::default(),
            power: 1.0,
            bounces: 0,
        }
    }
}

/// Planar non-sequential scene: each hit splits a ray into reflected, refracted and a
/// deterministic scattered child, pruned below `min_power`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NonSequentialScene {
    planes: Vec<NsPlane>,
}

impl NonSequentialScene {
    pub fn add_plane(&mut self, mut p: NsPlane) -> Result<()> {
        p.normal = p.normal.normalized()?;
        require(
            p.reflectance >= 0.0
                && p.transmittance >= 0.0
                && p.absorptance >= 0.0
                && p.scatter_fraction >= 0.0
                && p.reflectance + p.transmittance + p.absorptance + p.scatter_fraction
                    <= 1.0000001
                && p.refractive_index_before > 0.0
                && p.refractive_index_after > 0.0
                && p.point.is_finite(),
            "invalid non-sequential plane",
        )?;
        self.planes.push(p);
        Ok(())
    }

    /// Upstream defaults `max_bounces = 4`, `min_power = 1e-6`. Returns escaped rays followed by
    /// rays still active after the last generation.
    pub fn trace(&self, seed: NsRay, max_bounces: usize, min_power: f64) -> Result<Vec<NsRay>> {
        require(
            seed.ray.origin.is_finite() && seed.ray.direction.is_finite(),
            "nonfinite non-sequential ray",
        )?;
        let mut active = vec![seed];
        let mut finished = Vec::new();
        let mut depth = 0;
        while depth < max_bounces && !active.is_empty() {
            let mut next = Vec::new();
            for r in &active {
                let mut best = f64::INFINITY;
                let mut hit = None;
                for p in &self.planes {
                    let den = r.ray.direction.dot(p.normal);
                    if den.abs() < 1e-12 {
                        continue;
                    }
                    let d = Vec3::new(
                        p.point.x - r.ray.origin.x,
                        p.point.y - r.ray.origin.y,
                        p.point.z - r.ray.origin.z,
                    );
                    let t = d.dot(p.normal) / den;
                    if t > 1e-9 && t < best {
                        best = t;
                        hit = Some(p);
                    }
                }
                let Some(h) = hit else {
                    finished.push(*r);
                    continue;
                };
                let o = r.ray.origin;
                let dir = r.ray.direction;
                let point = Vec3::new(o.x + best * dir.x, o.y + best * dir.y, o.z + best * dir.z);
                let child = |direction: Vec3, factor: f64| NsRay {
                    ray: Ray {
                        origin: point,
                        direction,
                        ..r.ray
                    },
                    power: r.power * factor,
                    bounces: r.bounces + 1,
                };
                if h.reflectance > 0.0 && r.power * h.reflectance >= min_power {
                    next.push(child(reflect(dir, h.normal)?, h.reflectance));
                }
                if h.transmittance > 0.0 && r.power * h.transmittance >= min_power {
                    if let Some(t) = refract(
                        dir,
                        h.normal,
                        h.refractive_index_before,
                        h.refractive_index_after,
                    )? {
                        next.push(child(t, h.transmittance));
                    }
                }
                if h.scatter_fraction > 0.0 && r.power * h.scatter_fraction >= min_power {
                    let side = if dir.dot(h.normal) < 0.0 { -1.0 } else { 1.0 };
                    let s = Vec3::new(
                        h.normal.x * side + 0.35,
                        h.normal.y + 0.17,
                        h.normal.z * side,
                    )
                    .normalized()?;
                    next.push(child(s, h.scatter_fraction));
                }
            }
            active = next;
            depth += 1;
        }
        finished.extend(active);
        Ok(finished)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StrayLightSummary {
    pub total_power: f64,
    pub reflected_power: f64,
    pub transmitted_power: f64,
    pub scattered_power: f64,
    pub absorbed_power: f64,
}

/// Rays with any interaction count as "reflected", unscattered ones as "transmitted" (upstream
/// classification; `scattered_power` is left zero).
pub fn summarize_stray_light(rays: &[NsRay], launched_power: f64) -> Result<StrayLightSummary> {
    require(launched_power >= 0.0, "invalid launched optical power")?;
    let mut s = StrayLightSummary::default();
    for r in rays {
        s.total_power += r.power;
        if r.bounces > 0 {
            s.reflected_power += r.power;
        } else {
            s.transmitted_power += r.power;
        }
    }
    s.absorbed_power = (launched_power - s.total_power).max(0.0);
    Ok(s)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GhostPath {
    pub bounces: usize,
    pub power: f64,
    pub ray: Ray,
}

/// Paths with at least `minimum_bounces` (upstream default 2), strongest first.
pub fn rank_ghost_paths(rays: &[NsRay], minimum_bounces: usize) -> Vec<GhostPath> {
    let mut out: Vec<GhostPath> = rays
        .iter()
        .filter(|r| r.bounces >= minimum_bounces)
        .map(|r| GhostPath {
            bounces: r.bounces,
            power: r.power,
            ray: r.ray,
        })
        .collect();
    out.sort_by(|a, b| b.power.total_cmp(&a.power));
    out
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OptimizationResult {
    pub variables: Vec<f64>,
    pub merit: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// Pattern search: `+step`, then `-step` per coordinate; halve all steps after a sweep without
/// improvement, converged once the largest step is below `tolerance` (upstream defaults 200,
/// 1e-10).
pub fn coordinate_descent(
    merit: &dyn Fn(&[f64]) -> f64,
    initial: Vec<f64>,
    mut step: Vec<f64>,
    max_iterations: usize,
    tolerance: f64,
) -> Result<OptimizationResult> {
    let mut x = initial;
    require(
        !x.is_empty() && step.len() == x.len() && max_iterations > 0 && tolerance > 0.0,
        "invalid optical optimizer controls",
    )?;
    let mut best = merit(&x);
    for it in 0..max_iterations {
        let mut improved = false;
        for j in 0..x.len() {
            let mut trial = x.clone();
            trial[j] += step[j];
            let m = merit(&trial);
            if m < best {
                x = trial;
                best = m;
                improved = true;
                continue;
            }
            let mut trial = x.clone();
            trial[j] -= step[j];
            let m = merit(&trial);
            if m < best {
                x = trial;
                best = m;
                improved = true;
            }
        }
        if !improved {
            let mut mx = 0.0_f64;
            for h in &mut step {
                *h *= 0.5;
                mx = mx.max(h.abs());
            }
            if mx < tolerance {
                return Ok(OptimizationResult {
                    variables: x,
                    merit: best,
                    iterations: it + 1,
                    converged: true,
                });
            }
        }
    }
    Ok(OptimizationResult {
        variables: x,
        merit: best,
        iterations: max_iterations,
        converged: false,
    })
}

/// DE/rand/1/bin with `std::mt19937_64(seed)` draws in upstream order, so a seed reproduces the
/// C++ population exactly (upstream defaults 32, 200, 0.8, 0.9, 1).
#[allow(clippy::too_many_arguments)]
pub fn differential_evolution(
    merit: &dyn Fn(&[f64]) -> f64,
    lower: &[f64],
    upper: &[f64],
    population: usize,
    generations: usize,
    mutation: f64,
    crossover: f64,
    seed: u64,
) -> Result<OptimizationResult> {
    require(
        !lower.is_empty()
            && lower.len() == upper.len()
            && population >= 4
            && generations > 0
            && mutation > 0.0
            && (0.0..=1.0).contains(&crossover),
        "invalid differential-evolution controls",
    )?;
    require(
        lower
            .iter()
            .zip(upper)
            .all(|(l, u)| l < u && l.is_finite() && u.is_finite()),
        "invalid differential-evolution bounds",
    )?;
    let dim = lower.len();
    let mut rng = Mt19937_64::new(seed);
    let mut pop = vec![vec![0.0; dim]; population];
    let mut score = vec![0.0; population];
    for i in 0..population {
        for j in 0..dim {
            pop[i][j] = lower[j] + rng.uniform(0.0, 1.0) * (upper[j] - lower[j]);
        }
        score[i] = merit(&pop[i]);
    }
    let mut best = 0;
    for i in 1..population {
        if score[i] < score[best] {
            best = i;
        }
    }
    let pick = |rng: &mut Mt19937_64| rng.uniform_index(0, population - 1);
    for _ in 0..generations {
        for i in 0..population {
            let a = loop {
                let v = pick(&mut rng);
                if v != i {
                    break v;
                }
            };
            let b = loop {
                let v = pick(&mut rng);
                if v != i && v != a {
                    break v;
                }
            };
            let c = loop {
                let v = pick(&mut rng);
                if v != i && v != a && v != b {
                    break v;
                }
            };
            let mut trial = pop[i].clone();
            let forced = pick(&mut rng) % dim;
            for j in 0..dim {
                if rng.uniform(0.0, 1.0) < crossover || j == forced {
                    trial[j] =
                        (pop[a][j] + mutation * (pop[b][j] - pop[c][j])).clamp(lower[j], upper[j]);
                }
            }
            let value = merit(&trial);
            if value < score[i] {
                pop[i] = trial;
                score[i] = value;
                if score[i] < score[best] {
                    best = i;
                }
            }
        }
    }
    Ok(OptimizationResult {
        variables: pop[best].clone(),
        merit: score[best],
        iterations: generations,
        converged: true,
    })
}

/// Central differences with `h = rel * max(1, |x_i|)` (upstream default 1e-6).
pub fn finite_difference_sensitivity(
    merit: &dyn Fn(&[f64]) -> f64,
    point: &[f64],
    relative_step: f64,
) -> Result<Vec<f64>> {
    require(relative_step > 0.0, "invalid sensitivity controls")?;
    let mut x = point.to_vec();
    let mut g = vec![0.0; point.len()];
    for i in 0..point.len() {
        let h = relative_step * point[i].abs().max(1.0);
        x[i] = point[i] + h;
        let a = merit(&x);
        x[i] = point[i] - h;
        let b = merit(&x);
        x[i] = point[i];
        g[i] = (a - b) / (2.0 * h);
    }
    Ok(g)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MonteCarloStats {
    pub mean: f64,
    pub standard_deviation: f64,
    pub minimum: f64,
    pub maximum: f64,
    pub samples: usize,
}

/// Gaussian tolerance sampling with libstdc++-identical draws for a given seed.
pub fn monte_carlo_tolerance(
    merit: &dyn Fn(&[f64]) -> f64,
    nominal: &[f64],
    sigma: &[f64],
    samples: usize,
    seed: u64,
) -> Result<MonteCarloStats> {
    require(
        nominal.len() == sigma.len() && samples > 0,
        "invalid tolerance controls",
    )?;
    let mut rng = Mt19937_64::new(seed);
    let mut normal = StdNormal::default();
    let mut x = nominal.to_vec();
    let (mut sum, mut sum2) = (0.0, 0.0);
    let mut s = MonteCarloStats {
        minimum: f64::INFINITY,
        maximum: f64::NEG_INFINITY,
        samples,
        ..Default::default()
    };
    for _ in 0..samples {
        for i in 0..x.len() {
            x[i] = nominal[i] + sigma[i] * normal.sample(&mut rng);
        }
        let m = merit(&x);
        sum += m;
        sum2 += m * m;
        s.minimum = s.minimum.min(m);
        s.maximum = s.maximum.max(m);
    }
    s.mean = sum / samples as f64;
    s.standard_deviation = (sum2 / samples as f64 - s.mean * s.mean).max(0.0).sqrt();
    Ok(s)
}

/// Forward-mode derivative of a scalar function written over [`Dual<1>`].
pub fn autodiff_derivative(function: impl Fn(Dual<1>) -> Dual<1>, point: f64) -> f64 {
    function(Dual {
        value: point,
        derivative: [1.0],
    })
    .derivative[0]
}
