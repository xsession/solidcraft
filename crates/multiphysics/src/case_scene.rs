//! The 3-D picture of a `collocated_flow3d` grid case, and the edits made by pointing at it:
//! drag a computational-domain face by whole cells, pick parts (obstacles, porous blocks,
//! surface bodies) to give them a solid material, pick domain faces to make them an inlet,
//! outlet, wall or symmetry plane, and close part of a face with lids.
//!
//! Everything here is pure geometry on the case JSON. Coordinates are metres in the case
//! frame (the domain spans `0..lengths`). Each edit returns a new case plus the admission
//! issues it raises; the caller decides whether to keep it, as with feature drafts.
use schemars::JsonSchema;
use serde::Serialize;
use serde_json::{Map, Value, json};

use crate::{Error, MAX_CASE_BYTES, Result, validate_case, validate_json};

/// The six exterior patches, in companion order: index / 2 is the axis, index % 2 the side
/// (0 the minimum, 1 the maximum).
pub const PATCHES: [&str; 6] = ["left", "right", "bottom", "top", "front", "back"];
/// Most triangles of surface bodies the scene carries for display and picking.
pub const MAX_SCENE_TRIANGLES: usize = 60_000;
/// Most cells along one axis that a domain drag may reach.
pub const MAX_AXIS_CELLS: usize = 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum PartKind {
    Obstacle,
    Porous,
    Body,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Part {
    pub name: String,
    pub kind: PartKind,
    /// Axis-aligned bounds in metres. A grid part is exactly this box.
    pub min: [f64; 3],
    pub max: [f64; 3],
    /// Surface-body triangles (empty for grid parts).
    pub triangles: Vec<[[f64; 3]; 3]>,
    /// The part has a solid-material (thermal solid) record.
    pub solid: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum PatchRole {
    /// Prescribed non-zero velocity.
    Inlet,
    /// Fixed (usually ambient) pressure.
    Outlet,
    /// No-slip wall.
    Wall,
    /// Slip: a symmetry plane.
    Symmetry,
    /// The external fan's patch.
    Fan,
    /// Nothing set: the solver's default.
    Default,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Scene {
    pub lengths: [f64; 3],
    pub dimensions: [usize; 3],
    pub parts: Vec<Part>,
    /// Roles of [`PATCHES`], in that order.
    pub patches: [PatchRole; 6],
    /// Surface-body triangles beyond [`MAX_SCENE_TRIANGLES`] were left out.
    pub truncated: bool,
    /// Grid node coordinates per axis (evenly spaced unless the mesh has a local refinement).
    pub nodes: [Vec<f64>; 3],
    /// Rotating regions (MRF), as boxes in metres.
    pub rotating: Vec<RotatingBox>,
    /// Internal fans, as boxes in metres with their blowing direction.
    pub fans: Vec<FanBox>,
}

/// An internal fan as drawn: its box in metres and the unit direction it blows along.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct FanBox {
    pub name: String,
    pub min: [f64; 3],
    pub max: [f64; 3],
    pub direction: [f64; 3],
}

/// Unit vector of a fan direction (`+x` … `-z`).
pub fn fan_direction(text: &str) -> Option<[f64; 3]> {
    let index = crate::airflow::FAN_DIRECTIONS.iter().position(|d| *d == text)?;
    let mut v = [0.0; 3];
    v[index / 2] = if index % 2 == 0 { 1.0 } else { -1.0 };
    Some(v)
}

/// A rotating region as drawn: its box in metres, axis direction and speed.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct RotatingBox {
    pub name: String,
    pub min: [f64; 3],
    pub max: [f64; 3],
    pub axis_point: [f64; 3],
    pub axis: [f64; 3],
    pub rpm: f64,
}

impl Scene {
    /// Mean cell size along `axis` in metres.
    pub fn cell(&self, axis: usize) -> f64 {
        match (self.lengths.get(axis), self.dimensions.get(axis)) {
            (Some(length), Some(&cells)) if cells > 0 => length / cells as f64,
            _ => 0.0,
        }
    }

    /// Width of the cell against domain face `patch` (a [`PATCHES`] index): the step a face
    /// drag moves by.
    pub fn edge_cell(&self, patch: usize) -> f64 {
        edge_width(self.nodes.get(patch / 2).map_or(&[][..], Vec::as_slice), patch % 2 == 1)
    }

    /// The mesh has explicit (locally refined) coordinates.
    pub fn graded(&self) -> bool {
        self.nodes.iter().any(|n| n.windows(3).any(|w| ((w[2] - w[1]) - (w[1] - w[0])).abs() > 1e-9 * (w[2] - w[0])))
    }

    pub fn part(&self, name: &str) -> Option<&Part> {
        self.parts.iter().find(|part| part.name == name)
    }
}

/// What a pick ray hit.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Target {
    /// A part; `face` is the box face hit (as a [`PATCHES`] index), or `None` for a surface body.
    Part { name: String, face: Option<usize> },
    /// An exterior domain face, as a [`PATCHES`] index.
    Patch { index: usize },
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Hit {
    pub target: Target,
    /// Distance along the (unit) ray.
    pub distance: f64,
}

/// A proposed case and the admission issues it raises (empty when it would be accepted).
#[derive(Clone, Debug, Serialize, serde::Deserialize, JsonSchema)]
pub struct Edit {
    pub case: Value,
    pub issues: Vec<String>,
    /// RFC 6901 pointer to the record to review next, when there is one.
    pub pointer: Option<String>,
}

fn bad(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

pub(crate) fn finish(case: Value, pointer: Option<String>) -> Result<Edit> {
    validate_json(&case, MAX_CASE_BYTES)?;
    let issues = validate_case(&case).err().map(|e| vec![e.to_string()]).unwrap_or_default();
    Ok(Edit { case, issues, pointer })
}

fn flow_case(case: &Value) -> Result<()> {
    if !case.is_object() || case["solver"].as_str() != Some("collocated_flow3d") {
        return Err(bad("the 3-D case view needs a collocated_flow3d case"));
    }
    if case.pointer("/mesh/xy_shear").is_some_and(|value| value.as_f64() != Some(0.0)) {
        return Err(bad("sheared grids are not shown in 3-D"));
    }
    Ok(())
}

pub fn grid(case: &Value) -> Result<([usize; 3], [f64; 3])> {
    flow_case(case)?;
    let dimensions: [usize; 3] = serde_json::from_value(case.pointer("/mesh/dimensions").cloned().unwrap_or(Value::Null))
        .map_err(|_| bad("the 3-D case view needs three mesh dimensions"))?;
    let mut lengths = [0.0; 3];
    for (slot, axis) in lengths.iter_mut().zip(["x", "y", "z"]) {
        *slot = case
            .pointer(&format!("/mesh/lengths/{axis}"))
            .and_then(Value::as_f64)
            .filter(|v| v.is_finite() && *v > 0.0)
            .ok_or_else(|| bad("the 3-D case view needs finite positive mesh lengths"))?;
    }
    if dimensions.iter().any(|&n| n == 0 || n > MAX_AXIS_CELLS) {
        return Err(bad(format!("mesh dimensions must be 1..{MAX_AXIS_CELLS} per axis")));
    }
    Ok((dimensions, lengths))
}

fn edge_width(nodes: &[f64], upper: bool) -> f64 {
    let n = nodes.len();
    match (upper, n) {
        (_, 0 | 1) => 0.0,
        (false, _) => nodes[1] - nodes[0],
        (true, _) => nodes[n - 1] - nodes[n - 2],
    }
}

/// Grid node coordinates per axis: the mesh's explicit coordinates, or evenly spaced nodes.
pub fn nodes(case: &Value) -> Result<[Vec<f64>; 3]> {
    let (dimensions, lengths) = grid(case)?;
    let mut out: [Vec<f64>; 3] = Default::default();
    for (axis, name) in ["x", "y", "z"].into_iter().enumerate() {
        let explicit: Option<Vec<f64>> = case.pointer(&format!("/mesh/coordinates/{name}")).and_then(|v| serde_json::from_value(v.clone()).ok());
        out[axis] = match explicit {
            Some(list) if list.len() == dimensions[axis] + 1 && list.windows(2).all(|w| w[1] > w[0]) => list,
            Some(_) => return Err(bad(format!("mesh coordinates {name} disagree with the dimensions"))),
            None => (0..=dimensions[axis]).map(|i| lengths[axis] * i as f64 / dimensions[axis] as f64).collect(),
        };
    }
    Ok(out)
}

/// Cell range `[first, last)` along an axis that covers `[low, high]` (snapped outward to
/// nodes, at least one cell, clamped to the grid).
fn snap(nodes: &[f64], low: f64, high: f64) -> (usize, usize) {
    let n = nodes.len().saturating_sub(1).max(1);
    let tolerance = 1e-6 * (nodes.last().copied().unwrap_or(1.0) / n as f64);
    let first = nodes.partition_point(|v| *v <= low + tolerance).saturating_sub(1).min(n - 1);
    let last = nodes.partition_point(|v| *v < high - tolerance).clamp(first + 1, n);
    (first, last)
}

pub(crate) fn index3(value: &Value) -> Option<[usize; 3]> {
    serde_json::from_value(value.clone()).ok()
}

fn point3(value: &Value) -> Option<[f64; 3]> {
    serde_json::from_value::<[f64; 3]>(value.clone()).ok().filter(|p| p.iter().all(|v| v.is_finite()))
}

fn solid_names(case: &Value) -> Vec<&str> {
    case.pointer("/thermal/solids").and_then(Value::as_array).into_iter().flatten().filter_map(|row| row["obstacle"].as_str()).collect()
}

/// The role each exterior patch plays, read from the boundary lists and the fan.
pub fn patch_roles(case: &Value) -> [PatchRole; 6] {
    let mut roles = [PatchRole::Default; 6];
    let rows = |key: &str| case.get(key).and_then(Value::as_array).cloned().unwrap_or_default();
    let (velocity, pressure) = (rows("velocity_boundaries"), rows("pressure_boundaries"));
    let find = |rows: &[Value], patch: &str| rows.iter().find(|row| row["patch"].as_str() == Some(patch)).cloned();
    for (role, patch) in roles.iter_mut().zip(PATCHES) {
        let fixed_pressure = find(&pressure, patch).is_some_and(|row| row["kind"].as_str() == Some("fixed_value"));
        *role = if case.pointer("/fan/patch").and_then(Value::as_str) == Some(patch) {
            PatchRole::Fan
        } else {
            match find(&velocity, patch).map(|row| (row["kind"].as_str().map(str::to_string), row["value"].clone())) {
                Some((Some(kind), value)) if kind == "fixed_value" => {
                    let moving = ["x", "y", "z"].iter().any(|axis| value.get(axis).and_then(Value::as_f64).is_some_and(|v| v != 0.0));
                    if moving { PatchRole::Inlet } else { PatchRole::Wall }
                }
                Some((Some(kind), _)) if kind == "slip" => PatchRole::Symmetry,
                _ if fixed_pressure => PatchRole::Outlet,
                _ => PatchRole::Default,
            }
        };
    }
    roles
}

/// The picture of `case`: domain, parts and patch roles.
pub fn scene(case: &Value) -> Result<Scene> {
    validate_json(case, MAX_CASE_BYTES)?;
    let (dimensions, lengths) = grid(case)?;
    let nodes = nodes(case)?;
    let solids = solid_names(case);
    let mut parts = Vec::new();
    for (key, kind) in [("obstacles", PartKind::Obstacle), ("porous_blocks", PartKind::Porous)] {
        for row in case.get(key).and_then(Value::as_array).into_iter().flatten().take(256) {
            let (Some(name), Some(min), Some(max)) = (row["name"].as_str(), index3(&row["min"]), index3(&row["max"])) else { continue };
            let at = |axis: usize, i: usize| nodes[axis].get(i).copied();
            let (Some(x0), Some(y0), Some(z0), Some(x1), Some(y1), Some(z1)) =
                (at(0, min[0]), at(1, min[1]), at(2, min[2]), at(0, max[0]), at(1, max[1]), at(2, max[2]))
            else {
                continue;
            };
            parts.push(Part { name: name.into(), kind, min: [x0, y0, z0], max: [x1, y1, z1], triangles: Vec::new(), solid: solids.contains(&name) });
        }
    }
    let mut budget = MAX_SCENE_TRIANGLES;
    let mut truncated = false;
    for row in case.get("bodies").and_then(Value::as_array).into_iter().flatten().take(256) {
        let Some(name) = row["name"].as_str() else { continue };
        let vertices: Vec<[f64; 3]> = row["vertices"].as_array().into_iter().flatten().filter_map(point3).collect();
        let mut triangles = Vec::new();
        for t in row["triangles"].as_array().into_iter().flatten() {
            if budget == 0 {
                truncated = true;
                break;
            }
            let Some([a, b, c]) = index3(t) else { continue };
            if let (Some(a), Some(b), Some(c)) = (vertices.get(a), vertices.get(b), vertices.get(c)) {
                triangles.push([*a, *b, *c]);
                budget -= 1;
            }
        }
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for v in &vertices {
            for axis in 0..3 {
                min[axis] = min[axis].min(v[axis]);
                max[axis] = max[axis].max(v[axis]);
            }
        }
        if vertices.is_empty() {
            continue;
        }
        parts.push(Part { name: name.into(), kind: PartKind::Body, min, max, triangles, solid: solids.contains(&name) });
    }
    let vector = |v: &Value| ["x", "y", "z"].map(|a| v.get(a).and_then(Value::as_f64).unwrap_or(0.0));
    let mut rotating = Vec::new();
    for row in case.get("rotating_regions").and_then(Value::as_array).into_iter().flatten().take(64) {
        let (Some(name), Some(min), Some(max)) = (row["name"].as_str(), index3(&row["min"]), index3(&row["max"])) else { continue };
        let corner = |index: [usize; 3]| -> Option<[f64; 3]> {
            Some([nodes[0].get(index[0]).copied()?, nodes[1].get(index[1]).copied()?, nodes[2].get(index[2]).copied()?])
        };
        let (Some(min), Some(max)) = (corner(min), corner(max)) else { continue };
        rotating.push(RotatingBox {
            name: name.into(),
            min,
            max,
            axis_point: vector(&row["axis_point"]),
            axis: vector(&row["axis"]),
            rpm: row["rpm"].as_f64().unwrap_or(0.0),
        });
    }
    let mut fans = Vec::new();
    for row in case.get("internal_fans").and_then(Value::as_array).into_iter().flatten().take(64) {
        let (Some(name), Some(min), Some(max)) = (row["name"].as_str(), index3(&row["min"]), index3(&row["max"])) else { continue };
        let corner = |index: [usize; 3]| -> Option<[f64; 3]> {
            Some([nodes[0].get(index[0]).copied()?, nodes[1].get(index[1]).copied()?, nodes[2].get(index[2]).copied()?])
        };
        let (Some(min), Some(max), Some(direction)) = (corner(min), corner(max), row["direction"].as_str().and_then(fan_direction)) else { continue };
        fans.push(FanBox { name: name.into(), min, max, direction });
    }
    Ok(Scene { lengths, dimensions, parts, patches: patch_roles(case), truncated, nodes, rotating, fans })
}

// ---------------------------------------------------------------- picking

/// Entry distance and entry face of a ray through a box (slab method). A ray starting inside
/// the box reports its exit instead, so a box around the eye can still be picked.
fn ray_box(origin: [f64; 3], direction: [f64; 3], min: [f64; 3], max: [f64; 3]) -> Option<(f64, usize, f64, usize)> {
    let (mut near, mut far) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut near_face, mut far_face) = (0, 0);
    for axis in 0..3 {
        let (o, d) = (origin[axis], direction[axis]);
        if d.abs() < 1e-300 {
            if o < min[axis] || o > max[axis] {
                return None;
            }
            continue;
        }
        let (t0, t1) = ((min[axis] - o) / d, (max[axis] - o) / d);
        // Entering through the minimum face when travelling up the axis.
        let (enter, exit, enter_face, exit_face) = if t0 <= t1 { (t0, t1, 2 * axis, 2 * axis + 1) } else { (t1, t0, 2 * axis + 1, 2 * axis) };
        if enter > near {
            near = enter;
            near_face = enter_face;
        }
        if exit < far {
            far = exit;
            far_face = exit_face;
        }
    }
    (near <= far && far >= 0.0 && near.is_finite() && far.is_finite()).then_some((near, near_face, far, far_face))
}

fn ray_triangle(origin: [f64; 3], direction: [f64; 3], [a, b, c]: [[f64; 3]; 3]) -> Option<f64> {
    let sub = |p: [f64; 3], q: [f64; 3]| [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
    let cross = |p: [f64; 3], q: [f64; 3]| [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
    let dot = |p: [f64; 3], q: [f64; 3]| p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
    let (e1, e2) = (sub(b, a), sub(c, a));
    let p = cross(direction, e2);
    let det = dot(e1, p);
    if det.abs() < 1e-18 {
        return None;
    }
    let s = sub(origin, a);
    let u = dot(s, p) / det;
    let q = cross(s, e1);
    let v = dot(direction, q) / det;
    let t = dot(e2, q) / det;
    (u >= 0.0 && v >= 0.0 && u + v <= 1.0 && t >= 0.0 && t.is_finite()).then_some(t)
}

/// Everything a ray hits, nearest first: parts, then the domain faces it enters and leaves
/// (so a face behind the model can still be chosen by cycling through the hits).
pub fn pick(scene: &Scene, origin: [f64; 3], direction: [f64; 3]) -> Vec<Hit> {
    let norm = direction.iter().map(|v| v * v).sum::<f64>().sqrt();
    if !(norm.is_finite() && norm > 0.0) || !origin.iter().all(|v| v.is_finite()) {
        return Vec::new();
    }
    let direction = direction.map(|v| v / norm);
    let mut hits = Vec::new();
    for part in &scene.parts {
        if part.kind == PartKind::Body {
            let nearest = part.triangles.iter().filter_map(|t| ray_triangle(origin, direction, *t)).fold(f64::INFINITY, f64::min);
            if nearest.is_finite() {
                hits.push(Hit { target: Target::Part { name: part.name.clone(), face: None }, distance: nearest });
            }
        } else if let Some((near, face, far, back)) = ray_box(origin, direction, part.min, part.max) {
            let (distance, face) = if near >= 0.0 { (near, face) } else { (far, back) };
            hits.push(Hit { target: Target::Part { name: part.name.clone(), face: Some(face) }, distance });
        }
    }
    if let Some((near, face, far, back)) = ray_box(origin, direction, [0.0; 3], scene.lengths) {
        if near >= 0.0 {
            hits.push(Hit { target: Target::Patch { index: face }, distance: near });
        }
        hits.push(Hit { target: Target::Patch { index: back }, distance: far });
    }
    hits.sort_by(|a, b| a.distance.total_cmp(&b.distance));
    hits
}

// ---------------------------------------------------------------- domain drag

/// Move exterior face `patch` (a [`PATCHES`] index) outward by `cells` whole cells (inward
/// when negative). The cell size is kept. Moving a minimum face shifts every grid part,
/// surface body and point goal so they stay where they are relative to each other.
/// Parts the moved face would cut are refused.
pub fn move_domain_face(case: &Value, patch: usize, cells: i64) -> Result<Edit> {
    validate_json(case, MAX_CASE_BYTES)?;
    let (dimensions, lengths) = grid(case)?;
    let (axis, upper) = (patch / 2, patch % 2 == 1);
    if patch >= 6 {
        return Err(bad("unknown domain face"));
    }
    if cells == 0 {
        return finish(case.clone(), None);
    }
    let old = dimensions[axis] as i64;
    let new = old + cells;
    if new < 1 || new > MAX_AXIS_CELLS as i64 {
        return Err(bad(format!("the domain must keep 1..{MAX_AXIS_CELLS} cells along this axis")));
    }
    let _ = lengths;
    let explicit = case.pointer("/mesh/coordinates").is_some_and(|c| !c.is_null());
    let old_nodes = nodes(case)?[axis].clone();
    // New cells take the width of the cell they are added beside.
    let w = edge_width(&old_nodes, upper);
    let k = cells.unsigned_abs() as usize;
    let (new_nodes, offset): (Vec<f64>, f64) = match (upper, cells > 0) {
        (true, true) => {
            let last = old_nodes.last().copied().unwrap_or(0.0);
            (old_nodes.iter().copied().chain((1..=k).map(|j| last + j as f64 * w)).collect(), 0.0)
        }
        (true, false) => (old_nodes[..old_nodes.len() - k].to_vec(), 0.0),
        (false, true) => {
            let grow = k as f64 * w;
            ((0..k).map(|j| j as f64 * w).chain(old_nodes.iter().map(|v| v + grow)).collect(), grow)
        }
        (false, false) => {
            let cut = old_nodes[k];
            (old_nodes[k..].iter().map(|v| v - cut).collect(), -cut)
        }
    };
    // Index shift of everything when the minimum face moves (outward adds cells below).
    let shift = if upper { 0 } else { cells };
    let mut proposed = case.clone();
    for key in ["obstacles", "porous_blocks"] {
        let Some(rows) = proposed.get_mut(key).and_then(Value::as_array_mut) else { continue };
        for row in rows {
            let name = row["name"].as_str().unwrap_or("a part").to_string();
            let (Some(min), Some(max)) = (index3(&row["min"]), index3(&row["max"])) else {
                return Err(bad(format!("{name} has invalid grid bounds")));
            };
            let (lo, hi) = (min[axis] as i64 + shift, max[axis] as i64 + shift);
            if lo < 0 || hi > new {
                return Err(bad(format!("moving this face would cut {name}; move the part first")));
            }
            row["min"][axis] = json!(lo);
            row["max"][axis] = json!(hi);
        }
    }
    let new_length = new_nodes.last().copied().unwrap_or(0.0);
    if let Some(rows) = proposed.get_mut("bodies").and_then(Value::as_array_mut) {
        for row in rows {
            let name = row["name"].as_str().unwrap_or("a body").to_string();
            for vertex in row.get_mut("vertices").and_then(Value::as_array_mut).into_iter().flatten() {
                let moved = vertex.get(axis).and_then(Value::as_f64).map(|v| v + offset);
                match moved {
                    Some(v) if (-1e-12..=new_length + 1e-12).contains(&v) => vertex[axis] = json!(v),
                    _ => return Err(bad(format!("moving this face would cut {name}; move the body first"))),
                }
            }
        }
    }
    if let Some(rows) = proposed.get_mut("goals").and_then(Value::as_array_mut) {
        for row in rows {
            if let Some(point) = row.get_mut("point").filter(|p| p.is_array()) {
                let moved = point.get(axis).and_then(Value::as_f64).map(|v| v + offset);
                match moved {
                    Some(v) if (0.0..=new_length).contains(&v) => point[axis] = json!(v),
                    _ => return Err(bad("moving this face would leave a point goal outside the domain")),
                }
            }
        }
    }
    proposed["mesh"]["dimensions"][axis] = json!(new);
    proposed["mesh"]["lengths"][["x", "y", "z"][axis]] = json!(new_length);
    if explicit {
        proposed["mesh"]["coordinates"][["x", "y", "z"][axis]] = json!(new_nodes);
    }
    finish(proposed, Some("/mesh".into()))
}

// ---------------------------------------------------------------- materials

/// Give every named part a solid material with isotropic `conductivity` (W/m·K). Parts that
/// already have one keep their power and other fields; only the conductivity changes. A case
/// without a thermal block gets the starter fluid, as Insert Solid Material does.
pub fn assign_solid_material(case: &Value, names: &[String], conductivity: f64) -> Result<Edit> {
    validate_json(case, MAX_CASE_BYTES)?;
    flow_case(case)?;
    if names.is_empty() || names.len() > 128 {
        return Err(bad("pick 1..128 parts in the graphics window"));
    }
    if !(conductivity.is_finite() && conductivity > 0.0 && conductivity < 1e5) {
        return Err(bad("conductivity must be positive and finite"));
    }
    let known = |name: &str| {
        ["obstacles", "bodies"]
            .iter()
            .any(|key| case.get(*key).and_then(Value::as_array).is_some_and(|rows| rows.iter().any(|row| row["name"].as_str() == Some(name))))
    };
    let mut proposed = case.clone();
    if proposed.get("thermal").is_none_or(Value::is_null) {
        proposed["thermal"] =
            json!({"fluid":{"specific_heat_j_per_kg_k":1000.0,"conductivity_w_per_m_k":0.1},"inflow_temperature_k":293.15,"solids":[]});
    }
    let solids = proposed
        .get_mut("thermal")
        .and_then(|thermal| thermal.as_object_mut())
        .ok_or_else(|| bad("thermal must be an object"))?
        .entry("solids")
        .or_insert_with(|| json!([]))
        .as_array_mut()
        .ok_or_else(|| bad("thermal.solids must be an array"))?;
    let mut first = None;
    for name in names {
        if !known(name) {
            return Err(bad(format!("{name} is not an obstacle or surface body; porous blocks cannot be solids")));
        }
        let index = match solids.iter().position(|row| row["obstacle"].as_str() == Some(name)) {
            Some(index) => index,
            None => {
                if solids.len() >= 128 {
                    return Err(bad("thermal.solids exceeds 128 records"));
                }
                solids.push(json!({"obstacle":name,"power_w":0.0,"density_kg_per_m3":1000.0,"specific_heat_j_per_kg_k":1000.0}));
                solids.len() - 1
            }
        };
        let row = solids.get_mut(index).and_then(Value::as_object_mut).ok_or_else(|| bad("thermal solids must be objects"))?;
        row.insert("conductivity_w_per_m_k".into(), json!([conductivity, conductivity, conductivity]));
        first.get_or_insert(index);
    }
    finish(proposed, first.map(|index| format!("/thermal/solids/{index}")))
}

// ---------------------------------------------------------------- boundary conditions

/// What a picked domain face should become.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RoleInput {
    /// Uniform inflow at `speed` m/s along the inward normal.
    Inlet {
        speed: f64,
    },
    /// Fixed kinematic pressure (m²/s², 0 for ambient) with zero-gradient velocity.
    Outlet {
        pressure: f64,
    },
    Wall,
    Symmetry,
}

fn replace_row(case: &mut Value, key: &str, patch: &str, row: Option<Value>) -> Result<()> {
    let object = case.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
    let rows = object.entry(key).or_insert_with(|| json!([])).as_array_mut().ok_or_else(|| bad(format!("{key} must be an array")))?;
    rows.retain(|existing| existing["patch"].as_str() != Some(patch));
    if let Some(row) = row {
        rows.push(row);
    }
    Ok(())
}

/// Make exterior face `patch` (a [`PATCHES`] index) an inlet, outlet, wall or symmetry plane.
/// Its previous velocity and pressure conditions are replaced. The fan's own patch is refused:
/// change or remove the fan instead.
pub fn set_patch_role(case: &Value, patch: usize, role: RoleInput) -> Result<Edit> {
    validate_json(case, MAX_CASE_BYTES)?;
    flow_case(case)?;
    let name = *PATCHES.get(patch).ok_or_else(|| bad("unknown domain face"))?;
    if case.pointer("/fan/patch").and_then(Value::as_str) == Some(name) {
        return Err(bad("the external fan owns this face; edit or remove the fan to change it"));
    }
    let mut proposed = case.clone();
    let (velocity, pressure) = match role {
        RoleInput::Inlet { speed } => {
            if !(speed.is_finite() && speed > 0.0 && speed < 1e3) {
                return Err(bad("inlet speed must be positive and below 1000 m/s"));
            }
            let mut value = Map::new();
            for (index, axis) in ["x", "y", "z"].into_iter().enumerate() {
                let inward = if index != patch / 2 {
                    0.0
                } else if patch.is_multiple_of(2) {
                    speed
                } else {
                    -speed
                };
                value.insert(axis.into(), json!(inward));
            }
            (json!({"patch":name,"kind":"fixed_value","value":value}), None)
        }
        RoleInput::Outlet { pressure } => {
            if !pressure.is_finite() || pressure.abs() > 1e6 {
                return Err(bad("outlet pressure must be finite"));
            }
            (json!({"patch":name,"kind":"zero_gradient"}), Some(json!({"patch":name,"kind":"fixed_value","value":pressure})))
        }
        RoleInput::Wall => (json!({"patch":name,"kind":"fixed_value"}), None),
        RoleInput::Symmetry => (json!({"patch":name,"kind":"slip"}), None),
    };
    replace_row(&mut proposed, "velocity_boundaries", name, Some(velocity))?;
    replace_row(&mut proposed, "pressure_boundaries", name, pressure)?;
    finish(proposed, Some("/velocity_boundaries".into()))
}

// ---------------------------------------------------------------- lids

/// The two in-plane axes of a face normal to `axis`, in increasing order.
pub fn face_axes(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

/// Close face `patch` with one-cell-thick lid obstacles everywhere except the rectangular
/// opening `[open_min, open_max)` (cell indices along [`face_axes`]). The opening is where
/// the face's boundary condition then acts, like a lid with an inlet face in it. Lids are
/// named `<patch>_lid_1..4`; an empty band adds no lid.
pub fn create_lids(case: &Value, patch: usize, open_min: [usize; 2], open_max: [usize; 2]) -> Result<Edit> {
    validate_json(case, MAX_CASE_BYTES)?;
    let (dimensions, _) = grid(case)?;
    let name = *PATCHES.get(patch).ok_or_else(|| bad("unknown domain face"))?;
    let axis = patch / 2;
    let [u, v] = face_axes(axis);
    if open_min[0] >= open_max[0] || open_min[1] >= open_max[1] || open_max[0] > dimensions[u] || open_max[1] > dimensions[v] {
        return Err(bad("the opening must be a non-empty rectangle inside the face"));
    }
    if dimensions[axis] < 2 {
        return Err(bad("lids need at least two cells across the domain"));
    }
    let layer = if patch.is_multiple_of(2) { (0, 1) } else { (dimensions[axis] - 1, dimensions[axis]) };
    let full = [dimensions[u], dimensions[v]];
    // Bands around the opening: below and above along u (full v), then beside it along v.
    let bands = [
        ([0, 0], [open_min[0], full[1]]),
        ([open_max[0], 0], [full[0], full[1]]),
        ([open_min[0], 0], [open_max[0], open_min[1]]),
        ([open_min[0], open_max[1]], [open_max[0], full[1]]),
    ];
    let mut proposed = case.clone();
    let existing: Vec<String> =
        case.get("obstacles").and_then(Value::as_array).into_iter().flatten().filter_map(|row| row["name"].as_str().map(str::to_string)).collect();
    let mut first = None;
    let mut added = 0;
    for (lo, hi) in bands {
        if lo[0] >= hi[0] || lo[1] >= hi[1] {
            continue;
        }
        added += 1;
        let lid = format!("{name}_lid_{added}");
        if existing.contains(&lid) {
            return Err(bad(format!("{lid} already exists; delete the old lids first")));
        }
        let (mut min, mut max) = ([0usize; 3], [0usize; 3]);
        (min[axis], max[axis]) = layer;
        (min[u], max[u], min[v], max[v]) = (lo[0], hi[0], lo[1], hi[1]);
        let object = proposed.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
        let rows = object.entry("obstacles").or_insert_with(|| json!([])).as_array_mut().ok_or_else(|| bad("obstacles must be an array"))?;
        first.get_or_insert(rows.len());
        rows.push(json!({"name":lid,"min":min,"max":max}));
    }
    if added == 0 {
        return Err(bad("the opening covers the whole face; there is nothing to close"));
    }
    finish(proposed, first.map(|index| format!("/obstacles/{index}")))
}

// ---------------------------------------------------------------- CAD components

/// One object of an OBJ mesh: its name and triangles (polygons are fanned), in file units.
#[derive(Clone, Debug, PartialEq)]
pub struct MeshObject {
    pub name: String,
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[usize; 3]>,
}

/// Most objects and triangles one import may bring in.
pub const MAX_IMPORT_OBJECTS: usize = 64;
pub const MAX_IMPORT_TRIANGLES: usize = 200_000;

/// The objects of an OBJ file (`o` groups; vertices are global as OBJ numbers them). Only the
/// geometry lines are read; normals, materials and texture indices are ignored.
pub fn parse_obj(text: &str) -> Result<Vec<MeshObject>> {
    let mut all: Vec<[f64; 3]> = Vec::new();
    let mut groups: Vec<(String, Vec<[usize; 3]>)> = Vec::new();
    let mut triangles = 0usize;
    for line in text.lines() {
        let mut words = line.split_whitespace();
        match words.next() {
            Some("v") => {
                let xyz: Vec<f64> = words.take(3).filter_map(|w| w.parse::<f64>().ok()).collect();
                let [x, y, z] = xyz.as_slice() else { return Err(bad("an OBJ vertex needs three numbers")) };
                if ![*x, *y, *z].iter().all(|v| v.is_finite() && v.abs() < 1e9) {
                    return Err(bad("OBJ vertices must be finite"));
                }
                all.push([*x, *y, *z]);
            }
            Some("o") | Some("g") => {
                let name = words.collect::<Vec<_>>().join(" ");
                if groups.len() >= MAX_IMPORT_OBJECTS {
                    return Err(bad(format!("an import takes at most {MAX_IMPORT_OBJECTS} bodies")));
                }
                groups.push((name, Vec::new()));
            }
            Some("f") => {
                let mut corners = Vec::new();
                for word in words {
                    let first = word.split('/').next().unwrap_or("");
                    let index: i64 = first.parse().map_err(|_| bad("an OBJ face index is not a number"))?;
                    let resolved = if index > 0 { index - 1 } else { all.len() as i64 + index };
                    if resolved < 0 || resolved as usize >= all.len() {
                        return Err(bad("an OBJ face refers to a missing vertex"));
                    }
                    corners.push(resolved as usize);
                }
                if corners.len() < 3 {
                    return Err(bad("an OBJ face needs three corners"));
                }
                if groups.is_empty() {
                    groups.push(("body".into(), Vec::new()));
                }
                if let Some((_, faces)) = groups.last_mut() {
                    for k in 1..corners.len() - 1 {
                        faces.push([corners[0], corners[k], corners[k + 1]]);
                    }
                    triangles += corners.len() - 2;
                }
                if triangles > MAX_IMPORT_TRIANGLES {
                    return Err(bad(format!("an import takes at most {MAX_IMPORT_TRIANGLES} triangles")));
                }
            }
            _ => {}
        }
    }
    // Re-index each object over its own vertices.
    let mut out = Vec::new();
    for (name, faces) in groups {
        if faces.is_empty() {
            continue;
        }
        let mut map = std::collections::BTreeMap::new();
        let mut vertices = Vec::new();
        let mut local = Vec::with_capacity(faces.len());
        for face in faces {
            let mut t = [0usize; 3];
            for (slot, global) in t.iter_mut().zip(face) {
                *slot = *map.entry(global).or_insert_with(|| {
                    vertices.push(all.get(global).copied().unwrap_or([0.0; 3]));
                    vertices.len() - 1
                });
            }
            local.push(t);
        }
        out.push(MeshObject { name, vertices, triangles: local });
    }
    if out.is_empty() {
        return Err(bad("the CAD export has no faces"));
    }
    Ok(out)
}

/// How a CAD body enters the case.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, serde::Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Representation {
    /// A grid block when the body fills its bounding box (within 2 %), else a surface body.
    #[default]
    Auto,
    /// Always its bounding box, snapped to whole cells.
    Obstacle,
    /// Always its closed triangle surface.
    Surface,
}

/// Choices of a component import (lengths in millimetres, the CAD unit).
#[derive(Clone, Debug, Default, serde::Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportOptions {
    #[serde(default)]
    pub representation: Representation,
    /// The CAD point at the domain origin. Required unless `fit_domain`.
    #[serde(default)]
    pub origin_mm: Option<[f64; 3]>,
    /// Size the domain around the imported bodies plus `margin_mm` on every side.
    #[serde(default)]
    pub fit_domain: bool,
    #[serde(default)]
    pub margin_mm: f64,
    /// Cell size for a fitted domain; by default the case keeps about its cell count.
    #[serde(default)]
    pub cell_mm: Option<f64>,
    /// Remove the case's other obstacles, porous blocks and surface bodies (and their solid
    /// records and goals) first.
    #[serde(default)]
    pub replace_parts: bool,
    /// With `fit_domain`, the exact domain box (CAD mm) instead of the bodies' tessellated box
    /// plus the margin. The domain then keeps that size exactly, with whole cells along each
    /// axis (their sizes may differ slightly between axes), so walls on the box stay on the
    /// domain boundary.
    #[serde(default)]
    pub domain_mm: Option<[[f64; 3]; 2]>,
}

#[derive(Clone, Debug, Serialize, JsonSchema)]
pub struct Imported {
    pub edit: Edit,
    /// The CAD point at the domain origin that was used; pass it to later imports.
    pub origin_mm: [f64; 3],
    /// Case names given to the bodies, in import order, with how each was represented.
    pub parts: Vec<(String, PartKind)>,
    /// Things to review (clipped bodies, kept parts).
    pub notes: Vec<String>,
}

/// A case identifier for a CAD body name.
pub fn part_name(cad: &str) -> String {
    let mut name: String = cad.trim().chars().map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '-' { c } else { '_' }).take(60).collect();
    if !name.as_bytes().first().is_some_and(u8::is_ascii_alphanumeric) {
        name.insert_str(0, "body_");
    }
    if PATCHES.contains(&name.as_str()) {
        name.insert_str(0, "body_");
    }
    name
}

fn bounds_of(points: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for p in points {
        for a in 0..3 {
            min[a] = min[a].min(p[a]);
            max[a] = max[a].max(p[a]);
        }
    }
    (min, max)
}

/// Signed enclosed volume of a triangle surface (positive when wound outward).
fn volume(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> f64 {
    triangles
        .iter()
        .filter_map(|t| {
            let (a, b, c) = (vertices.get(t[0])?, vertices.get(t[1])?, vertices.get(t[2])?);
            Some(a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))
        })
        .sum::<f64>()
        / 6.0
}

/// Smallest positive bounding-box extent of any body (mm), or infinity without one.
fn thinnest_extent(objects: &[MeshObject]) -> f64 {
    objects
        .iter()
        .flat_map(|object| {
            let (low, high) = bounds_of(&object.vertices);
            [0, 1, 2].map(|a| high[a] - low[a])
        })
        .filter(|extent| extent.is_finite() && *extent > 0.0)
        .fold(f64::INFINITY, f64::min)
}

/// Largest speed the case prescribes (initial field and fixed-value inlets), in m/s.
fn prescribed_speed(case: &Value) -> f64 {
    let magnitude = |v: &Value| {
        let parts = ["x", "y", "z"].map(|axis| v.get(axis).and_then(Value::as_f64).unwrap_or(0.0));
        let speed = parts.iter().map(|p| p * p).sum::<f64>().sqrt();
        if speed.is_finite() { speed } else { 0.0 }
    };
    let inlets = case.get("velocity_boundaries").and_then(Value::as_array).into_iter().flatten().map(|row| magnitude(&row["value"]));
    inlets.chain(case.get("initial_velocity").map(magnitude)).fold(0.0, f64::max)
}

/// After a refit, the case's time step was chosen for its old, usually much coarser cells; on
/// the new ones it can be far past what the PISO steps converge with (the momentum solve then
/// fails). Divide it by a whole factor until it is within the 3-D explicit limits, with a 0.8
/// safety factor: diffusion ν·dt·Σ 1/h² ≤ ½ and Courant u·dt·Σ 1/h ≤ 1 over the three axes'
/// smallest cells. Take that many more steps so the simulated time stays, as far as the step
/// limits allow.
pub(crate) fn fit_time_step(case: &mut Value, notes: &mut Vec<String>) {
    const SAFETY: f64 = 0.8;
    if case.get("solver").and_then(Value::as_str) != Some("collocated_flow3d") {
        return;
    }
    let Some(dt) = case.pointer("/config/dt").and_then(Value::as_f64).filter(|dt| dt.is_finite() && *dt > 0.0) else { return };
    let Ok(grid_nodes) = nodes(case) else { return };
    let widths = grid_nodes.clone().map(|axis| axis.windows(2).map(|w| w[1] - w[0]).fold(f64::INFINITY, f64::min));
    if widths.iter().any(|h| !(h.is_finite() && *h > 0.0)) {
        return;
    }
    let h = widths.iter().copied().fold(f64::INFINITY, f64::min);
    let nu = case.pointer("/config/kinematic_viscosity").and_then(Value::as_f64).filter(|nu| nu.is_finite() && *nu > 0.0);
    let speed = prescribed_speed(case);
    let mut limit = f64::INFINITY;
    if let Some(nu) = nu {
        limit = limit.min(SAFETY * 0.5 / (nu * widths.iter().map(|h| 1.0 / (h * h)).sum::<f64>()));
    }
    if speed > 0.0 {
        limit = limit.min(SAFETY / (speed * widths.iter().map(|h| 1.0 / h).sum::<f64>()));
    }
    if !(limit.is_finite() && dt > limit) {
        return;
    }
    // A ratio that is a whole number up to rounding needs exactly that many steps.
    let factor = (dt / limit - 1e-9).ceil().max(1.0).min(crate::MAX_STEPS as f64) as u64;
    if factor <= 1 {
        return;
    }
    let new_dt = dt / factor as f64;
    case["config"]["dt"] = json!(new_dt);
    let mut note = format!("the time step was cut from {dt:e} s to {new_dt:e} s for {:.3} mm cells", h * 1000.0);
    if let Some(steps) = case.pointer("/output/steps").and_then(Value::as_u64) {
        let cells = grid_nodes.iter().map(|axis| axis.len().saturating_sub(1) as u64).product::<u64>().max(1);
        let most = crate::MAX_STEPS.min(crate::MAX_WORK / cells).max(1);
        let new_steps = steps.saturating_mul(factor).min(most);
        case["output"]["steps"] = json!(new_steps);
        if let Some(interval) = case.pointer("/output/output_interval").and_then(Value::as_u64) {
            case["output"]["output_interval"] = json!(interval.saturating_mul(factor).clamp(1, new_steps));
        }
        if new_steps == steps.saturating_mul(factor) {
            note.push_str(&format!(", with {new_steps} steps for the same simulated time"));
        } else {
            note.push_str(&format!(
                "; the step limit allows {new_steps} steps, so it simulates {:e} s instead of {:e} s",
                new_steps as f64 * new_dt,
                steps as f64 * dt
            ));
        }
    }
    notes.push(note);
}

/// Bring CAD bodies (an OBJ export in millimetres) into the case as grid blocks or surface
/// bodies, placing them by `origin_mm` or fitting the domain around them. A body with the same
/// case name as an existing part replaces that part's geometry and keeps its material.
pub fn import_components(case: &Value, obj_text: &str, options: &ImportOptions) -> Result<Imported> {
    validate_json(case, MAX_CASE_BYTES)?;
    let (mut dimensions, mut lengths) = grid(case)?;
    let objects = parse_obj(obj_text)?;
    let all: Vec<[f64; 3]> = objects.iter().flat_map(|o| o.vertices.iter().copied()).collect();
    let (low, high) = bounds_of(&all);
    let mut proposed = case.clone();
    let mut notes = Vec::new();
    let origin = if options.fit_domain {
        if !(options.margin_mm.is_finite() && options.margin_mm >= 0.0) {
            return Err(bad("the margin must be zero or positive"));
        }
        let (origin, size) = match options.domain_mm {
            Some([min, max]) if (0..3).all(|a| min[a].is_finite() && max[a].is_finite()) => (min, [0, 1, 2].map(|a| (max[a] - min[a]) / 1000.0)),
            Some(_) => return Err(bad("the domain box must be finite")),
            None => (low.map(|v| v - options.margin_mm), [0, 1, 2].map(|a| (high[a] - low[a] + 2.0 * options.margin_mm) / 1000.0)),
        };
        if size.iter().any(|v| !(v.is_finite() && *v > 0.0)) {
            return Err(bad("the bodies and margin give the domain no size; add a margin for flat bodies"));
        }
        let budget = dimensions.iter().product::<usize>().clamp(1, crate::MAX_CELLS as usize) as f64;
        let thinnest = thinnest_extent(&objects) / 1000.0;
        let mut h = match options.cell_mm {
            Some(cell) if cell.is_finite() && cell > 0.0 => cell / 1000.0,
            Some(_) => return Err(bad("the cell size must be positive")),
            // The case's cell count alone can leave a body thinner than one cell, or the
            // margin without a fluid cell. Resolve the thinnest body with two cells and the
            // margin with one; the cell limit below coarsens again when that is too many.
            None => {
                let mut h = (size[0] * size[1] * size[2] / budget).cbrt();
                if thinnest.is_finite() && thinnest > 0.0 {
                    h = h.min(thinnest / 2.0);
                }
                if options.margin_mm > 0.0 {
                    h = h.min(options.margin_mm / 1000.0);
                }
                h
            }
        };
        let count = |h: f64| size.map(|l| ((l / h) - 1e-9).ceil().max(1.0) as usize);
        let mut fitted = count(h);
        let mut guard = 0;
        while (fitted.iter().product::<usize>() as u64 > crate::MAX_CELLS || fitted.iter().any(|&n| n > MAX_AXIS_CELLS)) && guard < 200 {
            h *= 1.05;
            fitted = count(h);
            guard += 1;
        }
        if thinnest.is_finite() && thinnest < 2.0 * h {
            notes.push(format!(
                "the thinnest body ({:.3} mm) spans fewer than two {:.3} mm cells; set a smaller cell size to resolve it",
                thinnest * 1000.0,
                h * 1000.0
            ));
        }
        // The domain is whole cells: grow it to dimensions × h, centred on the bodies, or keep
        // an exact box and divide it into whole cells.
        let (fitted_lengths, origin) = if options.domain_mm.is_some() {
            (size, origin)
        } else {
            let grown = fitted.map(|n| n as f64 * h);
            (grown, [0, 1, 2].map(|a| origin[a] - (grown[a] - size[a]) * 500.0))
        };
        let kept = ["obstacles", "porous_blocks"].iter().map(|k| case.get(*k).and_then(Value::as_array).map_or(0, Vec::len)).sum::<usize>();
        if kept > 0 && !options.replace_parts {
            notes.push(format!("{kept} existing grid parts keep their cell indices (clipped to the refitted grid); check where they are"));
        }
        dimensions = fitted;
        lengths = fitted_lengths;
        proposed["mesh"]["dimensions"] = json!(dimensions);
        proposed["mesh"]["lengths"] = json!({"x": lengths[0], "y": lengths[1], "z": lengths[2]});
        // A fitted domain is a fresh uniform grid; refine it locally afterwards.
        if let Some(mesh) = proposed.get_mut("mesh").and_then(Value::as_object_mut)
            && mesh.remove("coordinates").is_some()
        {
            notes.push("the local mesh refinement was reset by fitting the domain".into());
        }
        origin
    } else {
        options
            .origin_mm
            .filter(|o| o.iter().all(|v| v.is_finite()))
            .ok_or_else(|| bad("no CAD origin for the domain yet: fit the domain to the components first"))?
    };
    if options.fit_domain && !options.replace_parts {
        // Kept grid parts must still fit the refitted grid: clip them, drop any left empty.
        for key in ["obstacles", "porous_blocks"] {
            let Some(rows) = proposed.get_mut(key).and_then(Value::as_array_mut) else { continue };
            rows.retain_mut(|row| {
                let (Some(min), Some(max)) = (index3(&row["min"]), index3(&row["max"])) else { return false };
                let lo = [0, 1, 2].map(|a| min[a].min(dimensions[a] - 1));
                let hi = [0, 1, 2].map(|a| max[a].min(dimensions[a]).max(lo[a] + 1));
                row["min"] = json!(lo);
                row["max"] = json!(hi);
                true
            });
        }
    }
    if options.replace_parts {
        for key in ["obstacles", "bodies", "porous_blocks"] {
            if let Some(object) = proposed.as_object_mut() {
                object.remove(key);
            }
        }
    }
    let to_case = |p: [f64; 3]| [0, 1, 2].map(|a| (p[a] - origin[a]) / 1000.0);
    let grid_nodes = nodes(&proposed)?;
    let mut parts = Vec::new();
    for object in &objects {
        let name = part_name(&object.name);
        if parts.iter().any(|(n, _)| *n == name) {
            return Err(bad(format!("two bodies become {name}; rename one in the CAD model")));
        }
        let vertices: Vec<[f64; 3]> = object.vertices.iter().map(|p| to_case(*p)).collect();
        let (min, max) = bounds_of(&vertices);
        if (0..3).any(|a| max[a] <= 0.0 || min[a] >= lengths[a]) {
            return Err(bad(format!("{name} lies outside the domain; fit the domain to the components or move it")));
        }
        let box_volume: f64 = (0..3).map(|a| max[a] - min[a]).product();
        let enclosed = volume(&vertices, &object.triangles);
        let boxy = box_volume > 0.0 && (enclosed - box_volume).abs() <= 0.02 * box_volume;
        let surface = match options.representation {
            Representation::Auto => !boxy,
            Representation::Obstacle => false,
            Representation::Surface => true,
        };
        // Remove any old part of that name (obstacle or body) before adding the new one.
        for key in ["obstacles", "bodies", "porous_blocks"] {
            if let Some(rows) = proposed.get_mut(key).and_then(Value::as_array_mut) {
                rows.retain(|row| row["name"].as_str() != Some(name.as_str()));
            }
        }
        // Tessellation and exact bounds agree to rounding; a body on the domain box is inside it.
        let slack = |a: usize| 1e-9_f64.max(1e-9 * lengths[a]);
        let clipped = (0..3).any(|a| min[a] < -slack(a) || max[a] > lengths[a] + slack(a));
        if surface {
            if clipped {
                return Err(bad(format!("{name} crosses the domain boundary; a surface body must lie inside it")));
            }
            if enclosed <= 0.0 {
                return Err(bad(format!("{name} is not a closed outward surface; check the CAD body is a solid")));
            }
            // The companion refuses a body that holds no cell centre; say so before a run.
            if (0..3).any(|a| !grid_nodes[a].windows(2).any(|w| (min[a]..=max[a]).contains(&((w[0] + w[1]) / 2.0)))) {
                notes.push(format!("{name} covers no cell centre, so the solver will refuse it; set a smaller cell size"));
            }
            let row = json!({"name": name, "vertices": vertices, "triangles": object.triangles});
            let object_map = proposed.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
            object_map.entry("bodies").or_insert_with(|| json!([])).as_array_mut().ok_or_else(|| bad("bodies must be an array"))?.push(row);
            parts.push((name, PartKind::Body));
        } else {
            if clipped {
                notes.push(format!("{name} reaches past the domain and was cut at its faces"));
            }
            let mut lo = [0usize; 3];
            let mut hi = [0usize; 3];
            for a in 0..3 {
                (lo[a], hi[a]) = snap(&grid_nodes[a], min[a], max[a]);
            }
            let row = json!({"name": name, "min": lo, "max": hi});
            let object_map = proposed.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
            object_map.entry("obstacles").or_insert_with(|| json!([])).as_array_mut().ok_or_else(|| bad("obstacles must be an array"))?.push(row);
            parts.push((name, PartKind::Obstacle));
        }
    }
    // Solid records whose part is gone would make the case unreadable.
    let names: Vec<String> = ["obstacles", "bodies"]
        .iter()
        .flat_map(|k| proposed.get(*k).and_then(Value::as_array).cloned().unwrap_or_default())
        .filter_map(|row| row["name"].as_str().map(str::to_string))
        .collect();
    if let Some(solids) = proposed.pointer_mut("/thermal/solids").and_then(Value::as_array_mut) {
        let before = solids.len();
        solids.retain(|row| row["obstacle"].as_str().is_some_and(|n| names.iter().any(|m| m == n)));
        if solids.len() < before {
            notes.push(format!("{} solid materials of removed parts were dropped", before - solids.len()));
        }
        // A thermal case needs a solid for every part: new parts get the starter material.
        let mut given = 0;
        for name in &names {
            if !solids.iter().any(|row| row["obstacle"].as_str() == Some(name)) {
                solids.push(json!({"obstacle":name,"conductivity_w_per_m_k":[1.0,1.0,1.0],"power_w":0.0,"density_kg_per_m3":1000.0,"specific_heat_j_per_kg_k":1000.0}));
                given += 1;
            }
        }
        if given > 0 {
            notes.push(format!("{given} new parts got the starter solid material (1 W/m·K, no power); set their materials"));
        }
    }
    // Goals on removed parts go too, and equations that use them.
    if let Some(goals) = proposed.get_mut("goals").and_then(Value::as_array_mut) {
        let mut dropped: Vec<String> = Vec::new();
        goals.retain(|row| {
            let target = row["target"].as_str().unwrap_or("");
            let gone = !target.is_empty() && !PATCHES.contains(&target) && !names.iter().any(|n| n == target);
            if gone {
                dropped.push(row["name"].as_str().unwrap_or("").to_string());
            }
            !gone
        });
        goals.retain(|row| {
            let uses = row["expression"]
                .as_str()
                .is_some_and(|e| e.split(|c: char| !(c.is_ascii_alphanumeric() || c == '_')).any(|word| dropped.iter().any(|d| d == word)));
            if uses {
                dropped.push(row["name"].as_str().unwrap_or("").to_string());
            }
            !uses
        });
        if !dropped.is_empty() {
            notes.push(format!("goals on removed parts were dropped: {}", dropped.join(", ")));
        }
    }
    if options.fit_domain {
        fit_time_step(&mut proposed, &mut notes);
    }
    let pointer = parts.first().map(|(name, kind)| {
        let key = if *kind == PartKind::Body { "bodies" } else { "obstacles" };
        let index =
            proposed.get(key).and_then(Value::as_array).and_then(|rows| rows.iter().position(|r| r["name"].as_str() == Some(name))).unwrap_or(0);
        format!("/{key}/{index}")
    });
    Ok(Imported { edit: finish(proposed, pointer)?, origin_mm: origin, parts, notes })
}

// ---------------------------------------------------------------- local mesh

/// Largest split of one cell by a local mesh region.
pub const MAX_REFINEMENT: usize = 4;

/// Set the grid's nodes and remap every grid part's indices through `maps` (old node index to
/// new node index, per axis). Lengths are unchanged.
fn regrid(case: &Value, new_nodes: [Vec<f64>; 3], maps: [Vec<usize>; 3], uniform: bool) -> Result<Value> {
    let dimensions = [0, 1, 2].map(|a| new_nodes[a].len().saturating_sub(1));
    if dimensions.iter().any(|&n| n == 0 || n > MAX_AXIS_CELLS) {
        return Err(bad(format!("the local mesh needs 1..{MAX_AXIS_CELLS} cells per axis")));
    }
    let total = dimensions.iter().map(|&n| n as u64).product::<u64>();
    if total > crate::MAX_CELLS {
        return Err(bad(format!("the local mesh would need {total} cells; the limit is {}", crate::MAX_CELLS)));
    }
    let mut proposed = case.clone();
    for key in ["obstacles", "porous_blocks"] {
        let Some(rows) = proposed.get_mut(key).and_then(Value::as_array_mut) else { continue };
        for row in rows {
            let (Some(min), Some(max)) = (index3(&row["min"]), index3(&row["max"])) else { return Err(bad("a grid part has invalid bounds")) };
            let mut lo = [0usize; 3];
            let mut hi = [0usize; 3];
            for a in 0..3 {
                let (Some(&l), Some(&h)) = (maps[a].get(min[a]), maps[a].get(max[a])) else { return Err(bad("a grid part lies outside the grid")) };
                // A part never collapses: it keeps at least one cell.
                (lo[a], hi[a]) = (l.min(dimensions[a] - 1), h.max(l + 1).min(dimensions[a]));
            }
            row["min"] = json!(lo);
            row["max"] = json!(hi);
        }
    }
    let mesh = proposed.get_mut("mesh").and_then(Value::as_object_mut).ok_or_else(|| bad("the case has no mesh"))?;
    mesh.insert("dimensions".into(), json!(dimensions));
    if uniform {
        mesh.remove("coordinates");
    } else {
        let [x, y, z] = new_nodes;
        mesh.insert("coordinates".into(), json!({"x": x, "y": y, "z": z}));
    }
    Ok(proposed)
}

/// Local mesh: split every cell whose span along an axis overlaps `[min, max]` (metres) into
/// `factor` equal cells. The grid is a tensor product, so a refined interval runs through the
/// whole domain along the other two axes. Existing nodes stay where they are, so every grid
/// part keeps its exact position; surface bodies and point goals are in metres and unchanged.
pub fn refine_region(case: &Value, min: [f64; 3], max: [f64; 3], factor: usize) -> Result<Edit> {
    validate_json(case, MAX_CASE_BYTES)?;
    grid(case)?;
    if !(2..=MAX_REFINEMENT).contains(&factor) {
        return Err(bad(format!("split cells by 2 to {MAX_REFINEMENT}")));
    }
    if !(0..3).all(|a| min[a].is_finite() && max[a].is_finite() && min[a] < max[a]) {
        return Err(bad("the region needs a finite, non-empty box"));
    }
    let old = nodes(case)?;
    let mut new_nodes: [Vec<f64>; 3] = Default::default();
    let mut maps: [Vec<usize>; 3] = Default::default();
    let mut touched = [false; 3];
    for a in 0..3 {
        let list = &old[a];
        let mut out = vec![list[0]];
        let mut map = vec![0];
        for w in list.windows(2) {
            let (left, right) = (w[0], w[1]);
            if right > min[a] && left < max[a] {
                touched[a] = true;
                for j in 1..factor {
                    out.push(left + (right - left) * j as f64 / factor as f64);
                }
            }
            out.push(right);
            map.push(out.len() - 1);
        }
        new_nodes[a] = out;
        maps[a] = map;
    }
    if !touched.iter().all(|t| *t) {
        return Err(bad("the region lies outside the domain"));
    }
    let proposed = regrid(case, new_nodes, maps, false)?;
    finish(proposed, Some("/mesh".into()))
}

/// Refine the cells around the named parts, plus `margin_cells` cells' worth of space around
/// them (at the mean cell size), by `factor`.
pub fn refine_around(case: &Value, names: &[String], factor: usize, margin_cells: f64) -> Result<Edit> {
    let scene = scene(case)?;
    if names.is_empty() {
        return Err(bad("pick the parts to refine around"));
    }
    if !(margin_cells.is_finite() && (0.0..=64.0).contains(&margin_cells)) {
        return Err(bad("the margin must be 0 to 64 cells"));
    }
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for name in names {
        let part = scene.part(name).ok_or_else(|| bad(format!("{name} is not a part of this case")))?;
        for a in 0..3 {
            min[a] = min[a].min(part.min[a]);
            max[a] = max[a].max(part.max[a]);
        }
    }
    // Strictly inside the part's outer faces, so a margin of 0 refines only the part's cells.
    let pad = [0, 1, 2].map(|a| margin_cells * scene.cell(a));
    let inset = [0, 1, 2].map(|a| 1e-9 * scene.lengths[a]);
    let low = [0, 1, 2].map(|a| min[a] - pad[a] + inset[a]);
    let high = [0, 1, 2].map(|a| max[a] + pad[a] - inset[a]);
    refine_region(case, low, high, factor)
}

/// Remove the local mesh: back to evenly spaced cells at the grid's coarsest size along each
/// axis. Grid parts snap to the nearest new nodes (each keeps at least one cell).
pub fn uniform_mesh(case: &Value) -> Result<Edit> {
    validate_json(case, MAX_CASE_BYTES)?;
    let (_, lengths) = grid(case)?;
    if case.pointer("/mesh/coordinates").is_none_or(Value::is_null) {
        return Err(bad("the mesh has no local refinement"));
    }
    let old = nodes(case)?;
    let mut new_nodes: [Vec<f64>; 3] = Default::default();
    let mut maps: [Vec<usize>; 3] = Default::default();
    for a in 0..3 {
        let coarsest = old[a].windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max);
        let n = ((lengths[a] / coarsest).round() as usize).max(1);
        let list: Vec<f64> = (0..=n).map(|i| lengths[a] * i as f64 / n as f64).collect();
        maps[a] = old[a].iter().map(|x| ((x / lengths[a]) * n as f64).round().clamp(0.0, n as f64) as usize).collect();
        new_nodes[a] = list;
    }
    let proposed = regrid(case, new_nodes, maps, true)?;
    finish(proposed, Some("/mesh".into()))
}

/// The domain face a CAD point (millimetres) lies on, within half a cell: how a face picked on
/// the CAD model becomes an inlet, outlet or wall.
pub fn patch_at(case: &Value, origin_mm: [f64; 3], point_mm: [f64; 3]) -> Option<usize> {
    let (_, lengths) = grid(case).ok()?;
    let grid_nodes = nodes(case).ok()?;
    let p = [0, 1, 2].map(|a| (point_mm[a] - origin_mm[a]) / 1000.0);
    let mut best: Option<(f64, usize)> = None;
    for face in 0..6usize {
        let axis = face / 2;
        let plane = if face.is_multiple_of(2) { 0.0 } else { lengths[axis] };
        let tolerance = edge_width(&grid_nodes[axis], !face.is_multiple_of(2)) * 0.5;
        let distance = (p[axis] - plane).abs();
        let inside = face_axes(axis).iter().all(|&a| p[a] >= -tolerance && p[a] <= lengths[a] + tolerance);
        if distance <= tolerance && inside && best.is_none_or(|(d, _)| distance < d) {
            best = Some((distance, face));
        }
    }
    best.map(|(_, face)| face)
}

// ---------------------------------------------------------------- rotating regions

// ---------------------------------------------------------------- internal fans

/// An internal fan in place of the picked parts: the cell block covering them becomes fluid
/// that the fan drives along `direction` (`+x` … `-z`) following `curve`, a list of (volume
/// flow m³/s, pressure rise Pa) with rising flow and falling pressure. The parts are removed
/// with their solid materials and the goals on them. The case switches on the PISO
/// neighbour update and balanced body forces, which an actuator fan needs to stay accurate.
pub fn add_internal_fan(case: &Value, name: &str, parts: &[String], direction: &str, curve: &[(f64, f64)]) -> Result<Edit> {
    let scene = scene(case)?;
    let bytes = name.as_bytes();
    if bytes.is_empty() || bytes.len() > 60 || !bytes.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_') {
        return Err(bad("an internal fan needs a name of letters, digits and underscores"));
    }
    if case.get("internal_fans").and_then(Value::as_array).into_iter().flatten().any(|row| row["name"].as_str() == Some(name)) {
        return Err(bad(format!("there is already an internal fan {name}")));
    }
    if parts.is_empty() {
        return Err(bad("pick the part that stands for the fan"));
    }
    if fan_direction(direction).is_none() {
        return Err(bad("the direction must be +x, -x, +y, -y, +z or -z"));
    }
    if curve.len() < 2
        || curve.len() > 64
        || !curve.iter().all(|(q, p)| q.is_finite() && p.is_finite() && *q >= 0.0 && *p >= 0.0)
        || !curve.windows(2).all(|w| w[1].0 > w[0].0 && w[1].1 < w[0].1)
    {
        return Err(bad("a fan curve needs 2 to 64 points with rising flow and falling pressure"));
    }
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for part in parts {
        let part = scene.part(part).ok_or_else(|| bad(format!("{part} is not a part of this case")))?;
        for a in 0..3 {
            low[a] = low[a].min(part.min[a]);
            high[a] = high[a].max(part.max[a]);
        }
    }
    let mut min = [0usize; 3];
    let mut max = [0usize; 3];
    for a in 0..3 {
        (min[a], max[a]) = snap(&scene.nodes[a], low[a], high[a]);
    }
    let mut proposed = case.clone();
    let object = proposed.as_object_mut().ok_or_else(|| bad("the case must be an object"))?;
    let picked = |row: &Value| row["name"].as_str().is_some_and(|n| parts.iter().any(|p| p == n));
    for key in ["obstacles", "porous_blocks", "bodies"] {
        if let Some(rows) = object.get_mut(key).and_then(Value::as_array_mut) {
            rows.retain(|row| !picked(row));
        }
    }
    if let Some(solids) = object.get_mut("thermal").and_then(|t| t.get_mut("solids")).and_then(Value::as_array_mut) {
        solids.retain(|row| !row["obstacle"].as_str().is_some_and(|n| parts.iter().any(|p| p == n)));
    }
    if let Some(goals) = object.get_mut("goals").and_then(Value::as_array_mut) {
        goals.retain(|row| !row["target"].as_str().is_some_and(|n| parts.iter().any(|p| p == n)));
    }
    let row = json!({
        "name": name,
        "min": min,
        "max": max,
        "direction": direction,
        "curve": curve.iter().map(|(q, p)| json!({"volume_flow_m3_per_s": q, "pressure_rise_pa": p})).collect::<Vec<_>>(),
    });
    let fans = object.entry("internal_fans").or_insert_with(|| json!([]));
    let fans = fans.as_array_mut().ok_or_else(|| bad("internal_fans must be a list"))?;
    fans.push(row);
    let index = fans.len() - 1;
    let config = object.entry("config").or_insert_with(|| json!({}));
    let config = config.as_object_mut().ok_or_else(|| bad("config must be an object"))?;
    config.insert("piso_neighbour_update".into(), json!(true));
    config.insert("balanced_body_force".into(), json!(true));
    finish(proposed, Some(format!("/internal_fans/{index}")))
}

/// Remove the named internal fan; its cells stay fluid.
pub fn remove_internal_fan(case: &Value, name: &str) -> Result<Edit> {
    let mut proposed = case.clone();
    let rows = proposed.get_mut("internal_fans").and_then(Value::as_array_mut).ok_or_else(|| bad("the case has no internal fans"))?;
    let before = rows.len();
    rows.retain(|row| row["name"].as_str() != Some(name));
    if rows.len() == before {
        return Err(bad(format!("there is no internal fan {name}")));
    }
    if rows.is_empty()
        && let Some(object) = proposed.as_object_mut()
    {
        object.remove("internal_fans");
    }
    finish(proposed, None)
}

/// The most cells of margin a rotating region may take around its parts.
pub const MAX_ROTATING_MARGIN: usize = 32;
/// Rotation per time step this authoring aims for (the companion accepts up to 0.2 rad).
const TARGET_ROTATION_PER_STEP: f64 = 0.1;
/// Most time steps a lowered time step may lead to.
const MAX_ROTATING_STEPS: u64 = 1_000_000;

/// A rotating region (MRF) around the picked parts: the cell block covering them plus
/// `margin_cells`, kept at least one cell off the domain boundary so the outer walls stay
/// still. The axis runs along `axis` (0 x, 1 y, 2 z) through the block's centre; positive
/// `rpm` turns by the right-hand rule. Walls of parts inside the region turn with it. The
/// case switches on the PISO neighbour update, and a time step too long for the speed is
/// shortened to 0.1 rad per step with the step count raised to keep the simulated time.
pub fn add_rotating_region(case: &Value, name: &str, parts: &[String], axis: usize, rpm: f64, margin_cells: usize) -> Result<Edit> {
    let scene = scene(case)?;
    let bytes = name.as_bytes();
    if bytes.is_empty() || bytes.len() > 60 || !bytes.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_') {
        return Err(bad("a rotating region needs a name of letters, digits and underscores"));
    }
    let existing = case.get("rotating_regions").and_then(Value::as_array).cloned().unwrap_or_default();
    if existing.iter().any(|row| row["name"].as_str() == Some(name)) {
        return Err(bad(format!("there is already a rotating region {name}")));
    }
    if parts.is_empty() {
        return Err(bad("pick the parts that rotate"));
    }
    if axis > 2 {
        return Err(bad("the axis must be x, y or z"));
    }
    if !(rpm.is_finite() && rpm != 0.0 && rpm.abs() <= 1e6) {
        return Err(bad("the speed must be a non-zero number of revolutions per minute"));
    }
    if margin_cells > MAX_ROTATING_MARGIN {
        return Err(bad(format!("the margin must be 0 to {MAX_ROTATING_MARGIN} cells")));
    }
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for part in parts {
        let part = scene.part(part).ok_or_else(|| bad(format!("{part} is not a part of this case")))?;
        for a in 0..3 {
            low[a] = low[a].min(part.min[a]);
            high[a] = high[a].max(part.max[a]);
        }
    }
    let mut min = [0usize; 3];
    let mut max = [0usize; 3];
    for a in 0..3 {
        let n = scene.dimensions[a];
        let (first, last) = snap(&scene.nodes[a], low[a], high[a]);
        // An axis one cell thick (a 2-D case) spans the whole grid; otherwise stay off the walls.
        if n == 1 {
            (min[a], max[a]) = (0, 1);
            continue;
        }
        let first = first.saturating_sub(margin_cells).max(1);
        let last = (last + margin_cells).min(n - 1);
        if first >= last {
            return Err(bad(format!(
                "the parts touch the domain boundary along {}; a rotating region needs a cell of fluid between it and the outer walls",
                ["x", "y", "z"][a]
            )));
        }
        (min[a], max[a]) = (first, last);
    }
    for row in &existing {
        let (Some(lo), Some(hi)) = (index3(&row["min"]), index3(&row["max"])) else { continue };
        if (0..3).all(|a| min[a] < hi[a] && lo[a] < max[a]) {
            return Err(bad(format!("the region would overlap rotating region {}", row["name"].as_str().unwrap_or("?"))));
        }
    }
    let centre = [0, 1, 2].map(|a| 0.5 * (scene.nodes[a][min[a]] + scene.nodes[a][max[a]]));
    let mut direction = [0.0; 3];
    direction[axis] = 1.0;
    let mut proposed = case.clone();
    let mut rows = existing;
    rows.push(json!({
        "name": name,
        "min": min,
        "max": max,
        "axis_point": {"x": centre[0], "y": centre[1], "z": centre[2]},
        "axis": {"x": direction[0], "y": direction[1], "z": direction[2]},
        "rpm": rpm,
    }));
    let index = rows.len() - 1;
    let object = proposed.as_object_mut().ok_or_else(|| bad("the case must be an object"))?;
    object.insert("rotating_regions".into(), Value::Array(rows));
    let config = object.entry("config").or_insert_with(|| json!({}));
    let config = config.as_object_mut().ok_or_else(|| bad("config must be an object"))?;
    config.insert("piso_neighbour_update".into(), json!(true));
    let omega = rpm.abs() * std::f64::consts::TAU / 60.0;
    let dt = config.get("dt").and_then(Value::as_f64).filter(|v| v.is_finite() && *v > 0.0).unwrap_or(1e-2);
    if omega * dt > TARGET_ROTATION_PER_STEP {
        let shorter = TARGET_ROTATION_PER_STEP / omega;
        config.insert("dt".into(), json!(shorter));
        if let Some(output) = object.get_mut("output").and_then(Value::as_object_mut) {
            let scale = dt / shorter;
            for key in ["steps", "output_interval"] {
                if let Some(v) = output.get(key).and_then(Value::as_u64) {
                    output.insert(key.into(), json!(((v as f64 * scale).ceil() as u64).clamp(1, MAX_ROTATING_STEPS)));
                }
            }
        }
    }
    finish(proposed, Some(format!("/rotating_regions/{index}")))
}

/// Remove the named rotating region. The PISO neighbour update stays on.
pub fn remove_rotating_region(case: &Value, name: &str) -> Result<Edit> {
    let mut proposed = case.clone();
    let rows = proposed.get_mut("rotating_regions").and_then(Value::as_array_mut).ok_or_else(|| bad("the case has no rotating regions"))?;
    let before = rows.len();
    rows.retain(|row| row["name"].as_str() != Some(name));
    if rows.len() == before {
        return Err(bad(format!("there is no rotating region {name}")));
    }
    if rows.is_empty()
        && let Some(object) = proposed.as_object_mut()
    {
        object.remove("rotating_regions");
    }
    finish(proposed, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enclosure() -> Value {
        serde_json::from_str(include_str!("../data/templates/flow-electronics-enclosure-cht.json")).unwrap_or(Value::Null)
    }

    /// A 20 × 10 × 10 mm box at (100, 0, 0) and a right triangular prism beside it, as the OBJ
    /// export writes them (millimetres, outward, welded, 1-based indices).
    fn cad_obj() -> String {
        let mut text = String::from("# SolidCraft OBJ export (mm)\no Heat Sink\n");
        let (x0, x1, y1, z1) = (100.0, 120.0, 10.0, 10.0);
        for (x, y, z) in [(x0, 0.0, 0.0), (x1, 0.0, 0.0), (x1, y1, 0.0), (x0, y1, 0.0), (x0, 0.0, z1), (x1, 0.0, z1), (x1, y1, z1), (x0, y1, z1)] {
            text += &format!("v {x} {y} {z}\n");
        }
        for f in ["1 4 3 2", "5 6 7 8", "1 2 6 5", "2 3 7 6", "3 4 8 7", "4 1 5 8"] {
            text += &format!("f {f}\n");
        }
        // Prism: triangle (130,0)-(150,0)-(130,10) in x-y, extruded 0..10 in z. Vertices 9..14.
        text += "o wedge\n";
        for (x, y, z) in [(130.0, 0.0, 0.0), (150.0, 0.0, 0.0), (130.0, 10.0, 0.0), (130.0, 0.0, 10.0), (150.0, 0.0, 10.0), (130.0, 10.0, 10.0)] {
            text += &format!("v {x} {y} {z}\n");
        }
        for f in ["9 11 10", "12 13 14", "9 10 13 12", "10 11 14 13", "11 9 12 14"] {
            text += &format!("f {f}\n");
        }
        text
    }

    #[test]
    fn obj_objects_are_read_and_reindexed() {
        let objects = parse_obj(&cad_obj()).expect("parse");
        assert_eq!(objects.len(), 2);
        assert_eq!(objects[0].name, "Heat Sink");
        assert_eq!((objects[0].vertices.len(), objects[0].triangles.len()), (8, 12));
        assert_eq!((objects[1].vertices.len(), objects[1].triangles.len()), (6, 8));
        assert!((volume(&objects[0].vertices, &objects[0].triangles) - 2000.0).abs() < 1e-9);
        assert!((volume(&objects[1].vertices, &objects[1].triangles) - 1000.0).abs() < 1e-9);
        assert!(parse_obj("v 0 0 0\nf 1 2 3\n").is_err());
        assert!(parse_obj("v 0 0\n").is_err());
        assert!(parse_obj("# nothing\n").is_err());
        assert_eq!(part_name("Heat Sink"), "Heat_Sink");
        assert_eq!(part_name("left"), "body_left");
        assert_eq!(part_name("_x"), "body__x");
    }

    #[test]
    fn components_fit_the_domain_and_choose_blocks_or_surfaces() {
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        let options = ImportOptions { fit_domain: true, margin_mm: 10.0, cell_mm: Some(5.0), ..Default::default() };
        let imported = import_components(&channel, &cad_obj(), &options).expect("import");
        assert!(imported.edit.issues.is_empty(), "{:?}", imported.edit.issues);
        // Bodies span x 100..150, y 0..10, z 0..10 mm; plus 10 mm each side: 70 × 30 × 30 mm.
        assert_eq!(imported.edit.case["mesh"]["dimensions"], json!([14, 6, 6]));
        assert_eq!(imported.origin_mm, [90.0, -10.0, -10.0]);
        assert_eq!(imported.parts, vec![("Heat_Sink".to_string(), PartKind::Obstacle), ("wedge".to_string(), PartKind::Body)]);
        // The box covers x 10..30 mm, y and z 10..20 mm of the domain: cells 2..6, 2..4, 2..4.
        assert_eq!(imported.edit.case["obstacles"][0], json!({"name":"Heat_Sink","min":[2,2,2],"max":[6,4,4]}));
        let scene = scene(&imported.edit.case).expect("scene");
        let wedge = scene.part("wedge").expect("wedge");
        assert!((wedge.min[0] - 0.04).abs() < 1e-12 && (wedge.max[0] - 0.06).abs() < 1e-12);
        // Importing again with the origin replaces geometry and keeps materials.
        let with_material = assign_solid_material(&imported.edit.case, &["Heat_Sink".into()], 200.0).expect("material");
        let again = ImportOptions { origin_mm: Some(imported.origin_mm), representation: Representation::Surface, ..Default::default() };
        let reimported = import_components(&with_material.case, &cad_obj(), &again).expect("reimport");
        assert!(reimported.edit.issues.is_empty(), "{:?}", reimported.edit.issues);
        assert_eq!(reimported.edit.case["obstacles"].as_array().map(Vec::len), Some(0));
        assert_eq!(reimported.edit.case["bodies"].as_array().map(Vec::len), Some(2));
        assert_eq!(reimported.edit.case["thermal"]["solids"][0]["obstacle"], "Heat_Sink");
        // Without a fitted origin there is nowhere to put them; a body outside is refused.
        assert!(import_components(&channel, &cad_obj(), &ImportOptions::default()).is_err());
        let far = ImportOptions { origin_mm: Some([10_000.0, 0.0, 0.0]), ..Default::default() };
        assert!(import_components(&imported.edit.case, &cad_obj(), &far).is_err());
    }

    /// An 80 × 50 × 8 mm plate, as the OBJ export writes it.
    fn plate_obj() -> String {
        let mut text = String::from("o Plate\n");
        for (x, y, z) in [(0, 0, 0), (80, 0, 0), (80, 50, 0), (0, 50, 0), (0, 0, 8), (80, 0, 8), (80, 50, 8), (0, 50, 8)] {
            text += &format!("v {x} {y} {z}\n");
        }
        for f in ["1 4 3 2", "5 6 7 8", "1 2 6 5", "2 3 7 6", "3 4 8 7", "4 1 5 8"] {
            text += &format!("f {f}\n");
        }
        text
    }

    #[test]
    fn a_default_fit_resolves_thin_bodies_and_fits_the_time_step_to_the_cells() {
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        let options =
            ImportOptions { fit_domain: true, margin_mm: 10.0, replace_parts: true, representation: Representation::Obstacle, ..Default::default() };
        let imported = import_components(&channel, &plate_obj(), &options).expect("import");
        assert!(imported.edit.issues.is_empty(), "{:?}", imported.edit.issues);
        // The channel's 96 cells alone would give 12.7 mm cells for a 100 × 70 × 28 mm domain,
        // thicker than the plate; the 8 mm plate gets two 4 mm cells (three once snapped outward)
        // and the margin keeps fluid around it.
        let case = &imported.edit.case;
        assert_eq!(case["mesh"]["dimensions"], json!([25, 18, 7]));
        assert_eq!(case["obstacles"][0], json!({"name":"Plate","min":[2,2,2],"max":[23,16,5]}));
        // The 3-D diffusion limit with its 0.8 safety: dt ≤ 0.4 / (ν · 3/h²) = 0.4 × 0.004² /
        // (3 × 0.02) = 1.0667e-4 s, so 0.002 s is cut by ⌈18.75⌉ = 19 and the steps go 12 → 228.
        assert!((case["config"]["dt"].as_f64().unwrap_or(0.0) - 0.002 / 19.0).abs() < 1e-15);
        assert_eq!((case["output"]["steps"].as_u64(), case["output"]["output_interval"].as_u64()), (Some(228), Some(76)));
        assert!(imported.notes.iter().any(|n| n.contains("time step was cut")), "{:?}", imported.notes);
        // With no margin the block is the whole grid: the structural check says so.
        let tight = ImportOptions { margin_mm: 0.0, ..options };
        let filled = import_components(&channel, &plate_obj(), &tight).expect("import");
        assert!(filled.edit.issues.iter().any(|i| i.contains("every fluid cell")), "{:?}", filled.edit.issues);
        // A surface body thinner than its cells is reported before a run.
        let coarse = ImportOptions { cell_mm: Some(20.0), representation: Representation::Surface, ..options };
        let thin = import_components(&channel, &plate_obj(), &coarse).expect("import");
        assert!(thin.notes.iter().any(|n| n.contains("covers no cell centre")), "{:?}", thin.notes);
    }

    #[test]
    fn replacing_parts_drops_their_materials_and_cad_faces_map_to_domain_faces() {
        let options = ImportOptions {
            fit_domain: true,
            margin_mm: 0.0,
            cell_mm: Some(5.0),
            replace_parts: true,
            representation: Representation::Obstacle,
            ..Default::default()
        };
        let imported = import_components(&enclosure(), &cad_obj(), &options).expect("import");
        assert!(imported.edit.issues.is_empty(), "{:?}", imported.edit.issues);
        let names: Vec<_> = imported.edit.case["obstacles"].as_array().into_iter().flatten().filter_map(|r| r["name"].as_str()).collect();
        assert_eq!(names, ["Heat_Sink", "wedge"]);
        // The old solids went with their parts; the new parts got the starter material.
        let solids: Vec<_> =
            imported.edit.case["thermal"]["solids"].as_array().into_iter().flatten().filter_map(|r| r["obstacle"].as_str()).collect();
        assert_eq!(solids, ["Heat_Sink", "wedge"]);
        assert!(imported.notes.iter().any(|n| n.contains("dropped")));
        let goals: Vec<_> = imported.edit.case["goals"].as_array().into_iter().flatten().filter_map(|r| r["name"].as_str()).collect();
        assert_eq!(goals, ["fan_flow"], "goals on the removed processor and power blocks went");
        // The fitted domain is the bodies' box: x 100..150 mm. A face at x = 100 is "left".
        let origin = imported.origin_mm;
        let case = &imported.edit.case;
        assert_eq!(patch_at(case, origin, [100.0, 5.0, 5.0]), Some(0));
        assert_eq!(patch_at(case, origin, [150.0, 5.0, 5.0]), Some(1));
        assert_eq!(patch_at(case, origin, [125.0, 5.0, 10.0]), Some(5));
        assert_eq!(patch_at(case, origin, [125.0, 5.0, 5.0]), None, "inside, on no domain face");
    }

    #[test]
    fn local_mesh_refines_around_parts_without_moving_them() {
        let base = enclosure();
        let before = scene(&base).expect("scene");
        let refined = refine_around(&base, &["processor_block".into()], 2, 1.0).expect("refine");
        assert!(refined.issues.is_empty(), "{:?}", refined.issues);
        let after = scene(&refined.case).expect("scene");
        assert!(after.graded() && !before.graded());
        // The processor spans 3 × 3 × 3 cells; with a margin of one cell, 5 cells per axis split in two.
        assert_eq!(after.dimensions, [25, 17, 13]);
        for part in &before.parts {
            let moved = after.part(&part.name).expect("part");
            assert!((0..3).all(|a| (moved.min[a] - part.min[a]).abs() < 1e-12 && (moved.max[a] - part.max[a]).abs() < 1e-12), "{}", part.name);
        }
        let x = &after.nodes[0];
        assert!((x[4] - 0.04).abs() < 1e-12 && (x[5] - 0.045).abs() < 1e-12, "{:?}", &x[..8]);
        assert!(crate::airflow::requirements(&refined.case).iter().any(|r| r.capability == "local_mesh_v1"));
        // Moving a face of the refined grid adds cells of the edge width and keeps parts in place.
        let grown = move_domain_face(&refined.case, 0, 2).expect("grow left");
        assert!(grown.issues.is_empty(), "{:?}", grown.issues);
        let shifted = scene(&grown.case).expect("scene");
        let gap = |s: &Scene| s.part("exhaust_grille").map_or(0.0, |p| p.min[0]) - s.part("processor_block").map_or(0.0, |p| p.max[0]);
        assert!((gap(&shifted) - gap(&after)).abs() < 1e-12);
        assert!((shifted.lengths[0] - 0.22).abs() < 1e-12);
        let shrunk = move_domain_face(&grown.case, 0, -2).expect("shrink left");
        assert_eq!(shrunk.case["obstacles"], refined.case["obstacles"]);
        // Back to uniform: the original grid and indices.
        let flat = uniform_mesh(&refined.case).expect("uniform");
        assert!(flat.issues.is_empty(), "{:?}", flat.issues);
        assert_eq!(flat.case["mesh"]["dimensions"], base["mesh"]["dimensions"]);
        assert_eq!(flat.case["obstacles"], base["obstacles"]);
        assert!(flat.case.pointer("/mesh/coordinates").is_none());
        assert!(uniform_mesh(&base).is_err());
        // Refusals: factor, an empty or outside region, and the cell budget.
        assert!(refine_region(&base, [0.0; 3], [0.1; 3], 5).is_err());
        assert!(refine_region(&base, [1.0; 3], [2.0; 3], 2).is_err());
        assert!(refine_region(&base, [0.0; 3], [1.0; 3], 4).is_err(), "64 × 1920 cells is over the limit");
        // Admission catches hand-edited coordinates.
        let mut broken = refined.case.clone();
        broken["mesh"]["coordinates"]["x"][3] = json!(0.0);
        assert!(crate::validate_case(&broken).is_err());
    }

    #[test]
    fn imports_snap_to_graded_cells() {
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        let options = ImportOptions { fit_domain: true, margin_mm: 10.0, cell_mm: Some(5.0), ..Default::default() };
        let imported = import_components(&channel, &cad_obj(), &options).expect("import");
        let refined = refine_around(&imported.edit.case, &["Heat_Sink".into()], 2, 0.0).expect("refine");
        let again = ImportOptions { origin_mm: Some(imported.origin_mm), representation: Representation::Obstacle, ..Default::default() };
        let reimported = import_components(&refined.case, &cad_obj(), &again).expect("reimport");
        assert!(reimported.edit.issues.is_empty(), "{:?}", reimported.edit.issues);
        let before = scene(&imported.edit.case).expect("scene");
        let after = scene(&reimported.edit.case).expect("scene");
        let sink = |s: &Scene| s.part("Heat_Sink").map(|p| (p.min, p.max));
        assert_eq!(sink(&before).map(|b| b.0.map(|v| (v * 1e9).round())), sink(&after).map(|b| b.0.map(|v| (v * 1e9).round())));
        assert_eq!(sink(&before).map(|b| b.1.map(|v| (v * 1e9).round())), sink(&after).map(|b| b.1.map(|v| (v * 1e9).round())));
        // Fitting the domain again resets the refinement, with a note.
        let refit = import_components(&refined.case, &cad_obj(), &options).expect("refit");
        assert!(refit.edit.case.pointer("/mesh/coordinates").is_none());
        assert!(refit.notes.iter().any(|n| n.contains("refinement")));
    }

    #[test]
    fn scene_reads_parts_roles_and_cell_sizes() {
        let scene = scene(&enclosure()).expect("scene");
        assert_eq!(scene.dimensions, [20, 12, 8]);
        assert!((scene.cell(0) - 0.01).abs() < 1e-12);
        let pcb = scene.part("pcb_block").expect("pcb");
        assert_eq!(pcb.kind, PartKind::Obstacle);
        assert!(pcb.solid);
        assert!((pcb.min[0] - 0.02).abs() < 1e-12 && (pcb.max[1] - 0.11).abs() < 1e-12);
        assert_eq!(scene.part("exhaust_grille").map(|p| p.kind), Some(PartKind::Porous));
        // left: fan, right: outlet, the rest no-slip walls.
        assert_eq!(scene.patches, [PatchRole::Fan, PatchRole::Outlet, PatchRole::Wall, PatchRole::Wall, PatchRole::Wall, PatchRole::Wall]);
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        assert_eq!(patch_roles(&channel)[0], PatchRole::Inlet);
        assert_eq!(patch_roles(&channel)[4], PatchRole::Symmetry);
        let sheared: Value = serde_json::from_str(include_str!("../data/templates/flow-cavity-simple.json")).unwrap_or(Value::Null);
        assert!(super::scene(&sheared).is_err());
    }

    #[test]
    fn picking_finds_the_nearest_part_face_and_every_domain_face_behind_it() {
        let scene = scene(&enclosure()).expect("scene");
        // Straight down (-y) through the processor block at x=0.065, z=0.03.
        let hits = pick(&scene, [0.065, 1.0, 0.03], [0.0, -1.0, 0.0]);
        assert_eq!(hits.first().map(|h| &h.target), Some(&Target::Patch { index: 3 }), "enter through top");
        let parts: Vec<_> =
            hits.iter().filter_map(|h| if let Target::Part { name, face } = &h.target { Some((name.as_str(), *face)) } else { None }).collect();
        assert_eq!(parts.first(), Some(&("processor_block", Some(3))), "{hits:?}");
        assert_eq!(hits.last().map(|h| &h.target), Some(&Target::Patch { index: 2 }), "leave through bottom");
        assert!(hits.windows(2).all(|w| w[0].distance <= w[1].distance));
        // A miss hits nothing; a zero direction is refused quietly.
        assert!(pick(&scene, [5.0, 5.0, 5.0], [1.0, 0.0, 0.0]).is_empty());
        assert!(pick(&scene, [0.0; 3], [0.0; 3]).is_empty());
        // From inside a box the exit face is reported.
        let inside = pick(&scene, [0.065, 0.045, 0.03], [1.0, 0.0, 0.0]);
        assert_eq!(inside.first().map(|h| &h.target), Some(&Target::Part { name: "processor_block".into(), face: Some(1) }));
    }

    #[test]
    fn surface_bodies_are_picked_by_triangle() {
        let mut case = enclosure();
        // A small closed tetrahedron in free air near the outlet side.
        case["bodies"] = json!([{"name":"probe","vertices":[[0.15,0.05,0.03],[0.17,0.05,0.03],[0.15,0.07,0.03],[0.15,0.05,0.05]],
            "triangles":[[0,2,1],[0,1,3],[0,3,2],[1,2,3]]}]);
        let scene = scene(&case).expect("scene");
        let hits = pick(&scene, [0.155, 0.055, 1.0], [0.0, 0.0, -1.0]);
        let first_part = hits.iter().find(|h| matches!(h.target, Target::Part { .. }));
        assert_eq!(first_part.map(|h| &h.target), Some(&Target::Part { name: "probe".into(), face: None }));
        // The slanted face x + y + z = 0.25 is at z = 0.04 under the ray.
        assert!(first_part.is_some_and(|h| (h.distance - 0.96).abs() < 1e-9), "{hits:?}");
    }

    #[test]
    fn domain_faces_move_by_whole_cells_and_keep_parts_in_place() {
        let base = enclosure();
        let grown = move_domain_face(&base, 1, 5).expect("grow right");
        assert!(grown.issues.is_empty(), "{:?}", grown.issues);
        assert_eq!(grown.case["mesh"]["dimensions"], json!([25, 12, 8]));
        assert!((grown.case["mesh"]["lengths"]["x"].as_f64().unwrap_or(0.0) - 0.25).abs() < 1e-12);
        assert_eq!(grown.case["obstacles"], base["obstacles"], "a maximum face leaves indices alone");
        // Growing the minimum face shifts every part and the point goals.
        let mut with_point = base.clone();
        with_point["goals"] = json!([{"name":"p","kind":"point_temperature","point":[0.1,0.06,0.04]}]);
        let shifted = move_domain_face(&with_point, 0, 3).expect("grow left");
        assert_eq!(shifted.case["obstacles"][5]["min"], json!([8, 3, 2]));
        assert_eq!(shifted.case["porous_blocks"][0]["max"], json!([23, 12, 8]));
        assert!((shifted.case["goals"][0]["point"][0].as_f64().unwrap_or(0.0) - 0.13).abs() < 1e-12);
        // The same scene, seen from the parts: the processor stays 0.05 m from the grille.
        let before = scene(&base).expect("scene");
        let after = scene(&shifted.case).expect("scene");
        let gap = |s: &Scene| s.part("exhaust_grille").map(|p| p.min[0]).unwrap_or(0.0) - s.part("processor_block").map(|p| p.max[0]).unwrap_or(0.0);
        assert!((gap(&before) - gap(&after)).abs() < 1e-12);
        // Shrinking into a part is refused; the shroud sits at x = 0..1.
        assert!(move_domain_face(&base, 0, -1).is_err());
        assert!(move_domain_face(&base, 1, -1).is_err(), "the grille fills the last layer");
        assert!(move_domain_face(&base, 3, -100).is_err());
        assert!(move_domain_face(&base, 6, 1).is_err());
        let same = move_domain_face(&base, 2, 0).expect("no-op");
        assert_eq!(same.case, base);
    }

    #[test]
    fn materials_go_to_every_picked_part_and_keep_existing_power() {
        let base = enclosure();
        let names = vec!["processor_block".to_string(), "shroud_lower".to_string()];
        let edit = assign_solid_material(&base, &names, 200.0).expect("assign");
        assert!(edit.issues.is_empty(), "{:?}", edit.issues);
        let solids = edit.case["thermal"]["solids"].as_array().cloned().unwrap_or_default();
        assert_eq!(solids.len(), base["thermal"]["solids"].as_array().map_or(0, Vec::len));
        for row in &solids {
            if names.iter().any(|n| row["obstacle"].as_str() == Some(n)) {
                assert_eq!(row["conductivity_w_per_m_k"], json!([200.0, 200.0, 200.0]));
            }
        }
        let processor = solids.iter().find(|r| r["obstacle"] == "processor_block").cloned().unwrap_or(Value::Null);
        let original = base["thermal"]["solids"].as_array().and_then(|rows| rows.iter().find(|r| r["obstacle"] == "processor_block")).cloned();
        assert_eq!(Some(processor["power_w"].clone()), original.map(|r| r["power_w"].clone()));
        // A case without a thermal block gets one, with new records.
        let mut cold = base.clone();
        cold.as_object_mut().map(|o| o.remove("thermal"));
        cold.as_object_mut().map(|o| o.remove("goals"));
        let fresh = assign_solid_material(&cold, &names, 1.0).expect("fresh");
        assert_eq!(fresh.case["thermal"]["solids"].as_array().map(Vec::len), Some(2));
        assert_eq!(fresh.pointer.as_deref(), Some("/thermal/solids/0"));
        // Porous blocks, unknown names and bad conductivities are refused.
        assert!(assign_solid_material(&base, &["exhaust_grille".into()], 1.0).is_err());
        assert!(assign_solid_material(&base, &["nothing".into()], 1.0).is_err());
        assert!(assign_solid_material(&base, &names, 0.0).is_err());
        assert!(assign_solid_material(&base, &[], 1.0).is_err());
    }

    #[test]
    fn picked_faces_become_inlets_outlets_walls_or_symmetry() {
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        // Turn the flow around: right becomes the inlet (flowing -x), left the outlet.
        let edit = set_patch_role(&channel, 1, RoleInput::Inlet { speed: 0.2 }).expect("inlet");
        let edit = set_patch_role(&edit.case, 0, RoleInput::Outlet { pressure: 0.0 }).expect("outlet");
        assert!(edit.issues.is_empty(), "{:?}", edit.issues);
        let roles = patch_roles(&edit.case);
        assert_eq!((roles[0], roles[1]), (PatchRole::Outlet, PatchRole::Inlet));
        let inlet = edit.case["velocity_boundaries"].as_array().and_then(|r| r.iter().find(|r| r["patch"] == "right")).cloned();
        assert_eq!(inlet.map(|r| r["value"].clone()), Some(json!({"x":-0.2,"y":0.0,"z":0.0})));
        assert_eq!(edit.case["pressure_boundaries"].as_array().map(Vec::len), Some(1), "the old right pressure is gone");
        let walled = set_patch_role(&edit.case, 4, RoleInput::Wall).expect("wall");
        assert_eq!(patch_roles(&walled.case)[4], PatchRole::Wall);
        let symmetric = set_patch_role(&walled.case, 4, RoleInput::Symmetry).expect("symmetry");
        assert_eq!(patch_roles(&symmetric.case)[4], PatchRole::Symmetry);
        // The fan's patch is refused, and so are bad values.
        assert!(set_patch_role(&enclosure(), 0, RoleInput::Wall).is_err());
        assert!(set_patch_role(&channel, 0, RoleInput::Inlet { speed: -1.0 }).is_err());
        assert!(set_patch_role(&channel, 9, RoleInput::Wall).is_err());
    }

    #[test]
    fn lids_close_a_face_around_its_opening() {
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        // Left face of an 8x6x2 grid: in-plane axes y (6) and z (2). Opening y 2..4, z 0..2.
        let edit = create_lids(&channel, 0, [2, 0], [4, 2]).expect("lids");
        assert!(edit.issues.is_empty(), "{:?}", edit.issues);
        let lids = edit.case["obstacles"].as_array().cloned().unwrap_or_default();
        assert_eq!(lids.len(), 2, "z is fully open, so only the two y bands remain");
        assert_eq!(lids[0], json!({"name":"left_lid_1","min":[0,0,0],"max":[1,2,2]}));
        assert_eq!(lids[1], json!({"name":"left_lid_2","min":[0,4,0],"max":[1,6,2]}));
        // Closed cells + open cells = the whole face.
        let closed: usize = lids
            .iter()
            .map(|l| (1..3).map(|a| l["max"][a].as_u64().unwrap_or(0) - l["min"][a].as_u64().unwrap_or(0)).product::<u64>() as usize)
            .sum();
        assert_eq!(closed + 2 * 2, 6 * 2);
        let scene = scene(&edit.case).expect("scene");
        assert_eq!(scene.parts.len(), 2);
        // Top face of the enclosure, opening in the middle: four bands, on the last y layer.
        let top = create_lids(&channel, 3, [2, 0], [5, 1]).expect("top");
        let rows = top.case["obstacles"].as_array().cloned().unwrap_or_default();
        assert!(rows.iter().all(|r| r["min"][1] == 5 && r["max"][1] == 6));
        assert!(create_lids(&edit.case, 0, [2, 0], [4, 2]).is_err(), "names are not reused");
        assert!(create_lids(&channel, 0, [0, 0], [6, 2]).is_err(), "a fully open face needs no lid");
        assert!(create_lids(&channel, 0, [3, 0], [3, 2]).is_err());
        assert!(create_lids(&channel, 0, [0, 0], [7, 2]).is_err());
    }

    #[test]
    fn rotating_regions_wrap_parts_off_the_walls_and_follow_the_speed() {
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        let options = ImportOptions { fit_domain: true, margin_mm: 10.0, cell_mm: Some(5.0), ..Default::default() };
        let mut case = import_components(&channel, &cad_obj(), &options).expect("import").edit.case;
        // Heat_Sink covers cells 2..6, 2..4, 2..4 of a 14 × 6 × 6 grid.
        case["config"]["dt"] = json!(0.01);
        case["output"] = json!({"steps": 1000, "output_interval": 250});
        let edit = add_rotating_region(&case, "rotor", &["Heat_Sink".into()], 2, 600.0, 1).expect("region");
        assert!(edit.issues.is_empty(), "{:?}", edit.issues);
        assert_eq!(edit.pointer.as_deref(), Some("/rotating_regions/0"));
        let row = &edit.case["rotating_regions"][0];
        assert_eq!((row["min"].clone(), row["max"].clone()), (json!([1, 1, 1]), json!([7, 5, 5])));
        // Centre of cells 1..7, 1..5, 1..5 at 5 mm: 20, 15, 15 mm.
        assert_eq!(row["axis"], json!({"x": 0.0, "y": 0.0, "z": 1.0}));
        assert!((row["axis_point"]["x"].as_f64().unwrap_or(0.0) - 0.02).abs() < 1e-12);
        assert!((row["axis_point"]["y"].as_f64().unwrap_or(0.0) - 0.015).abs() < 1e-12);
        assert_eq!(edit.case["config"]["piso_neighbour_update"], json!(true));
        // 600 rpm = 62.8 rad/s: Δt 0.01 s would turn 0.63 rad a step. Shortened to 0.1 rad,
        // with the steps raised to keep the 10 s simulated.
        let dt = edit.case["config"]["dt"].as_f64().unwrap_or(0.0);
        assert!((dt * 600.0 * std::f64::consts::TAU / 60.0 - 0.1).abs() < 1e-12, "{dt}");
        let steps = edit.case["output"]["steps"].as_u64().unwrap_or(0);
        assert!((steps as f64 * dt - 10.0).abs() < dt, "{steps}");
        let needs: Vec<_> = crate::airflow::requirements(&edit.case).iter().map(|r| r.capability).collect();
        assert!(needs.contains(&"rotating_regions_v1") && needs.contains(&"piso_neighbour_update_v1"), "{needs:?}");
        // The scene draws it in metres.
        let scene = scene(&edit.case).expect("scene");
        assert_eq!(scene.rotating.len(), 1);
        assert!((scene.rotating[0].min[0] - 0.005).abs() < 1e-12 && (scene.rotating[0].max[0] - 0.035).abs() < 1e-12);
        assert_eq!(scene.rotating[0].rpm, 600.0);
        // A large margin stops one cell short of the walls.
        let wide = add_rotating_region(&case, "fan", &["Heat_Sink".into()], 0, -60.0, 20).expect("wide");
        assert_eq!(
            (wide.case["rotating_regions"][0]["min"].clone(), wide.case["rotating_regions"][0]["max"].clone()),
            (json!([1, 1, 1]), json!([13, 5, 5]))
        );
        assert_eq!(wide.case["config"]["dt"], json!(0.01), "slow enough: the time step stays");
        // Refusals: overlap, duplicate names, bad input, unknown parts.
        assert!(add_rotating_region(&edit.case, "other", &["Heat_Sink".into()], 2, 60.0, 0).is_err());
        assert!(add_rotating_region(&edit.case, "rotor", &["wedge".into()], 2, 60.0, 0).is_err());
        assert!(add_rotating_region(&case, "r", &[], 2, 60.0, 0).is_err());
        assert!(add_rotating_region(&case, "r", &["Heat_Sink".into()], 3, 60.0, 0).is_err());
        assert!(add_rotating_region(&case, "r", &["Heat_Sink".into()], 2, 0.0, 0).is_err());
        assert!(add_rotating_region(&case, "r s", &["Heat_Sink".into()], 2, 60.0, 0).is_err());
        assert!(add_rotating_region(&case, "r", &["nothing".into()], 2, 60.0, 0).is_err());
        // Removal.
        let removed = remove_rotating_region(&edit.case, "rotor").expect("remove");
        assert!(removed.case.get("rotating_regions").is_none());
        assert!(remove_rotating_region(&removed.case, "rotor").is_err());
    }

    #[test]
    fn admission_checks_rotating_regions() {
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        let options = ImportOptions { fit_domain: true, margin_mm: 10.0, cell_mm: Some(5.0), ..Default::default() };
        let mut case = import_components(&channel, &cad_obj(), &options).expect("import").edit.case;
        case["config"]["dt"] = json!(0.001);
        let good = json!({"name": "r", "min": [1, 1, 1], "max": [7, 5, 5], "axis_point": {"x": 0.02, "y": 0.015, "z": 0.0}, "axis": {"x": 0.0, "y": 0.0, "z": 1.0}, "rpm": 600.0});
        let check = |rows: Value| {
            let mut c = case.clone();
            c["rotating_regions"] = rows;
            validate_case(&c)
        };
        assert!(check(json!([good.clone()])).is_ok());
        let with = |key: &str, value: Value| {
            let mut row = good.clone();
            row[key] = value;
            row
        };
        assert!(check(json!([with("max", json!([15, 5, 5]))])).is_err(), "outside the grid");
        assert!(check(json!([with("max", json!([1, 5, 5]))])).is_err(), "empty");
        assert!(check(json!([with("axis", json!({"x": 0.0, "y": 0.0, "z": 0.0}))])).is_err(), "zero axis");
        assert!(check(json!([with("rpm", json!(4000.0))])).is_err(), "0.42 rad a step");
        assert!(check(json!([with("name", json!("a b"))])).is_err(), "name");
        assert!(check(json!([with("extra", json!(1))])).is_err(), "unknown field");
        let mut other = with("name", json!("s"));
        other["min"] = json!([6, 1, 1]);
        other["max"] = json!([9, 3, 3]);
        assert!(check(json!([good.clone(), other.clone()])).is_err(), "overlap");
        other["min"] = json!([7, 1, 1]);
        assert!(check(json!([good, other])).is_ok(), "touching is fine");
        let mut flag = case.clone();
        flag["config"]["piso_neighbour_update"] = json!("yes");
        assert!(validate_case(&flag).is_err());
    }

    #[test]
    fn an_internal_fan_replaces_the_picked_part_and_turns_on_balanced_forces() {
        let channel: Value = serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null);
        let options = ImportOptions { fit_domain: true, margin_mm: 10.0, cell_mm: Some(5.0), ..Default::default() };
        let imported = import_components(&channel, &cad_obj(), &options).expect("import").edit.case;
        let case = assign_solid_material(&imported, &["Heat_Sink".into(), "wedge".into()], 200.0).expect("material").case;
        assert!(validate_case(&case).is_ok());
        let curve = [(0.0, 30.0), (0.002, 12.0), (0.004, 0.0)];
        let edit = add_internal_fan(&case, "blower", &["Heat_Sink".into()], "+x", &curve).expect("fan");
        assert!(edit.issues.is_empty(), "{:?}", edit.issues);
        assert_eq!(edit.pointer.as_deref(), Some("/internal_fans/0"));
        let fan = &edit.case["internal_fans"][0];
        // Heat_Sink was cells 2..6, 2..4, 2..4: the fan takes that block.
        assert_eq!((fan["min"].clone(), fan["max"].clone()), (json!([2, 2, 2]), json!([6, 4, 4])));
        assert_eq!(fan["curve"][1], json!({"volume_flow_m3_per_s": 0.002, "pressure_rise_pa": 12.0}));
        assert!(edit.case["obstacles"].as_array().is_some_and(|rows| rows.iter().all(|r| r["name"] != "Heat_Sink")));
        assert!(edit.case.pointer("/thermal/solids").and_then(Value::as_array).is_none_or(|rows| rows.iter().all(|r| r["obstacle"] != "Heat_Sink")));
        assert_eq!(edit.case["config"]["balanced_body_force"], json!(true));
        assert_eq!(edit.case["config"]["piso_neighbour_update"], json!(true));
        let needs: Vec<_> = crate::airflow::requirements(&edit.case).iter().map(|r| r.capability).collect();
        assert!(needs.contains(&"internal_fans_v1") && needs.contains(&"balanced_body_force_v1"), "{needs:?}");
        let scene = scene(&edit.case).expect("scene");
        assert_eq!(scene.fans.len(), 1);
        assert_eq!(scene.fans[0].direction, [1.0, 0.0, 0.0]);
        assert!((scene.fans[0].min[0] - 0.01).abs() < 1e-12 && (scene.fans[0].max[0] - 0.03).abs() < 1e-12);
        assert!(scene.part("Heat_Sink").is_none());
        // Refusals.
        assert!(add_internal_fan(&edit.case, "blower", &["wedge".into()], "+x", &curve).is_err(), "duplicate name");
        assert!(add_internal_fan(&case, "f", &[], "+x", &curve).is_err());
        assert!(add_internal_fan(&case, "f", &["Heat_Sink".into()], "x", &curve).is_err());
        assert!(add_internal_fan(&case, "f", &["Heat_Sink".into()], "+x", &[(0.0, 1.0)]).is_err());
        assert!(add_internal_fan(&case, "f", &["Heat_Sink".into()], "+x", &[(0.0, 1.0), (1.0, 2.0)]).is_err(), "rising pressure");
        // Admission catches a fan that overlaps a solid.
        let mut overlapping = edit.case.clone();
        overlapping["obstacles"] = json!([{"name": "lump", "min": [3, 2, 2], "max": [4, 3, 3]}]);
        assert!(validate_case(&overlapping).is_err());
        let mut bad_direction = edit.case.clone();
        bad_direction["internal_fans"][0]["direction"] = json!("up");
        assert!(validate_case(&bad_direction).is_err());
        let removed = remove_internal_fan(&edit.case, "blower").expect("remove");
        assert!(removed.case.get("internal_fans").is_none());
        assert!(remove_internal_fan(&removed.case, "blower").is_err());
    }
}
