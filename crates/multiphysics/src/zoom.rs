//! Zooming (FloEFD's transferred boundary conditions): a child case cut out of a parent run's
//! grid, optionally refined, whose faces inside the parent take the parent's results.
//!
//! - The child box is snapped outward to the parent's grid nodes, so every child cell lies
//!   exactly inside parent cells. A refinement factor then splits each child cell.
//! - Grid parts are clipped to the box, and parts left empty are dropped. Surface bodies and
//!   point goals are shifted into the child frame; bodies that cross the box are dropped. Solid
//!   records and goals of dropped parts go too.
//! - A child face that lies on the parent's own boundary keeps the parent's conditions there.
//!   A face inside the parent gets the parent's results: if the parent flow enters the child
//!   through it on average, a fixed-value velocity with the sampled velocity profile; else a
//!   fixed pressure with the sampled pressure profile and a zero-gradient velocity.
//! - Temperatures are not transferred: the child keeps the parent's inflow temperature.
use schemars::JsonSchema;
use serde::Serialize;
use serde_json::{Value, json};

use crate::case_scene::{PATCHES, face_axes, finish, grid, index3, nodes};
use crate::flowpost::{Locator, sample};
use crate::results::{CellIndex, field};
use crate::vtk::Grid;
use crate::{Error, MAX_CASE_BYTES, Result, validate_json};

/// Velocity and pressure fields of a `collocated_flow3d` run.
pub const VELOCITY_FIELD: &str = "velocity_m_per_s";
pub const PRESSURE_FIELD: &str = "kinematic_pressure_m2_per_s2";

fn bad(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

#[derive(Clone, Debug, Serialize, JsonSchema)]
pub struct Zoom {
    pub case: Value,
    /// Where the child's origin lies in the parent frame (metres).
    pub offset_m: [f64; 3],
    /// Child faces that take the parent's results, with "velocity" or "pressure".
    pub transferred: Vec<(String, String)>,
    pub notes: Vec<String>,
    /// Admission issues of the child case (empty when it would be accepted).
    pub issues: Vec<String>,
}

/// Cut the box `[low, high]` (metres, parent frame) out of `parent` and feed its inner faces
/// from the run `grid` (the parent's VTU). `refine` (1..=4) splits each child cell.
pub fn zoom_case(parent: &Value, grid_vtu: &Grid, low: [f64; 3], high: [f64; 3], refine: usize) -> Result<Zoom> {
    validate_json(parent, MAX_CASE_BYTES)?;
    let (dimensions, _) = grid(parent)?;
    if !(1..=crate::case_scene::MAX_REFINEMENT).contains(&refine) {
        return Err(bad(format!("refine by 1 to {}", crate::case_scene::MAX_REFINEMENT)));
    }
    if !(0..3).all(|a| low[a].is_finite() && high[a].is_finite() && low[a] < high[a]) {
        return Err(bad("the zoom box needs a finite, non-empty extent"));
    }
    let parent_nodes = nodes(parent)?;
    // Snap outward to parent nodes: cells [first, last) per axis.
    let mut first = [0usize; 3];
    let mut last = [0usize; 3];
    for a in 0..3 {
        let list = &parent_nodes[a];
        let tolerance = 1e-9 * list.last().copied().unwrap_or(1.0);
        first[a] = list.partition_point(|v| *v <= low[a] + tolerance).saturating_sub(1).min(dimensions[a] - 1);
        last[a] = list.partition_point(|v| *v < high[a] - tolerance).clamp(first[a] + 1, dimensions[a]);
    }
    if (0..3).all(|a| first[a] == 0 && last[a] == dimensions[a]) {
        return Err(bad("the zoom box covers the whole parent domain"));
    }
    let offset = [0, 1, 2].map(|a| parent_nodes[a][first[a]]);
    let child_nodes: [Vec<f64>; 3] = [0, 1, 2].map(|a| {
        let mut out = vec![0.0];
        for w in parent_nodes[a][first[a]..=last[a]].windows(2) {
            for j in 1..=refine {
                out.push(w[0] + (w[1] - w[0]) * j as f64 / refine as f64 - offset[a]);
            }
        }
        out
    });
    let child_dims = [0, 1, 2].map(|a| child_nodes[a].len() - 1);
    let total = child_dims.iter().map(|&n| n as u64).product::<u64>();
    if total > crate::MAX_CELLS {
        return Err(bad(format!("the zoomed grid would need {total} cells; the limit is {}", crate::MAX_CELLS)));
    }
    let mut child = parent.clone();
    let mut notes = Vec::new();
    let mut dropped: Vec<String> = Vec::new();
    // Grid parts: clip to the box, map to child cells.
    for key in ["obstacles", "porous_blocks"] {
        let Some(rows) = child.get_mut(key).and_then(Value::as_array_mut) else { continue };
        rows.retain_mut(|row| {
            let name = row["name"].as_str().unwrap_or("").to_string();
            let (Some(min), Some(max)) = (index3(&row["min"]), index3(&row["max"])) else { return false };
            let lo = [0, 1, 2].map(|a| min[a].clamp(first[a], last[a]) - first[a]);
            let hi = [0, 1, 2].map(|a| max[a].clamp(first[a], last[a]) - first[a]);
            if (0..3).any(|a| hi[a] <= lo[a]) {
                dropped.push(name);
                return false;
            }
            row["min"] = json!(lo.map(|v| v * refine));
            row["max"] = json!(hi.map(|v| v * refine));
            true
        });
    }
    let extent = [0, 1, 2].map(|a| child_nodes[a].last().copied().unwrap_or(0.0));
    let inside = |p: [f64; 3]| (0..3).all(|a| p[a] >= -1e-12 && p[a] <= extent[a] + 1e-12);
    if let Some(rows) = child.get_mut("bodies").and_then(Value::as_array_mut) {
        rows.retain_mut(|row| {
            let name = row["name"].as_str().unwrap_or("").to_string();
            let moved: Option<Vec<[f64; 3]>> = row["vertices"].as_array().map(|list| {
                list.iter().filter_map(|v| serde_json::from_value::<[f64; 3]>(v.clone()).ok()).map(|p| [0, 1, 2].map(|a| p[a] - offset[a])).collect()
            });
            match moved {
                Some(points) if !points.is_empty() && points.iter().all(|p| inside(*p)) => {
                    row["vertices"] = json!(points);
                    true
                }
                _ => {
                    dropped.push(name);
                    false
                }
            }
        });
    }
    if !dropped.is_empty() {
        notes.push(format!("parts outside or across the zoom box were left out: {}", dropped.join(", ")));
    }
    if let Some(solids) = child.pointer_mut("/thermal/solids").and_then(Value::as_array_mut) {
        solids.retain(|row| row["obstacle"].as_str().is_none_or(|n| !dropped.iter().any(|d| d == n)));
    }
    if let Some(goals) = child.get_mut("goals").and_then(Value::as_array_mut) {
        let before = goals.len();
        goals.retain_mut(|row| {
            if row["target"].as_str().is_some_and(|t| dropped.iter().any(|d| d == t)) {
                return false;
            }
            if let Some(point) = row.get("point").and_then(|p| serde_json::from_value::<[f64; 3]>(p.clone()).ok()) {
                let moved = [0, 1, 2].map(|a| point[a] - offset[a]);
                if !inside(moved) {
                    return false;
                }
                row["point"] = json!(moved);
            }
            true
        });
        // Equation goals that lost an operand would be refused; drop them as well.
        let names: Vec<String> = goals.iter().filter_map(|g| g["name"].as_str().map(str::to_string)).collect();
        goals.retain(|row| {
            row["expression"].as_str().is_none_or(|e| {
                e.split(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
                    .filter(|w| !w.is_empty() && !w.chars().next().is_some_and(|c| c.is_ascii_digit() || c == '.'))
                    .all(|w| names.iter().any(|n| n == w))
            })
        });
        if goals.len() < before {
            notes.push(format!("{} goals outside the zoom box were left out", before - goals.len()));
        }
    }
    // The child's faces: copy the parent's conditions where they coincide, else transfer.
    let index = CellIndex::new(grid_vtu).map_err(|e| bad(e.to_string()))?;
    let locator = Locator::new(&index);
    let velocity = field(grid_vtu, VELOCITY_FIELD).map_err(|e| bad(e.to_string()))?;
    let pressure = field(grid_vtu, PRESSURE_FIELD).map_err(|e| bad(e.to_string()))?;
    let mut transferred = Vec::new();
    let mut profiles = Vec::new();
    for (patch, name) in PATCHES.iter().enumerate() {
        let axis = patch / 2;
        let upper = patch % 2 == 1;
        let on_parent_boundary = if upper { last[axis] == dimensions[axis] } else { first[axis] == 0 };
        if on_parent_boundary {
            continue;
        }
        if child.pointer("/fan/patch").and_then(Value::as_str) == Some(name) {
            return Err(bad("the fan's face would lie inside the parent; keep it on the zoom box boundary"));
        }
        let [u, v] = face_axes(axis);
        let plane = if upper { extent[axis] } else { 0.0 };
        let centre = |nodes: &[f64], i: usize| 0.5 * (nodes[i] + nodes[i + 1]);
        let inward = if upper { -1.0 } else { 1.0 };
        let mut velocities = Vec::new();
        let mut pressures = Vec::new();
        let mut inflow = 0.0;
        let mut area = 0.0;
        for jv in 0..child_dims[v] {
            for ju in 0..child_dims[u] {
                let mut point = [0.0; 3];
                point[axis] = plane + offset[axis];
                point[u] = centre(&child_nodes[u], ju) + offset[u];
                point[v] = centre(&child_nodes[v], jv) + offset[v];
                let a = (child_nodes[u][ju + 1] - child_nodes[u][ju]) * (child_nodes[v][jv + 1] - child_nodes[v][jv]);
                let vel = sample(grid_vtu, &index, &locator, velocity, point).and_then(|(_, val)| <[f64; 3]>::try_from(val.as_slice()).ok());
                let pre = sample(grid_vtu, &index, &locator, pressure, point).and_then(|(_, val)| val.first().copied());
                if let Some(vel) = vel {
                    inflow += vel[axis] * inward * a;
                    area += a;
                }
                velocities.push(vel);
                pressures.push(pre);
            }
        }
        if area <= 0.0 {
            return Err(bad(format!("the parent has no flow results on the child's {name} face (it lies in a solid)")));
        }
        let mean_pressure = {
            let valid: Vec<f64> = pressures.iter().flatten().copied().collect();
            valid.iter().sum::<f64>() / valid.len().max(1) as f64
        };
        let replace = |key: &str, row: Option<Value>, child: &mut Value| -> Result<()> {
            let object = child.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
            let rows = object.entry(key).or_insert_with(|| json!([])).as_array_mut().ok_or_else(|| bad(format!("{key} must be an array")))?;
            rows.retain(|r| r["patch"].as_str() != Some(name));
            if let Some(row) = row {
                rows.push(row);
            }
            Ok(())
        };
        if inflow > 0.0 {
            // Faces in the parent's solids are walls.
            let list: Vec<Value> =
                velocities.iter().map(|v| v.map_or(json!({"x":0.0,"y":0.0,"z":0.0}), |v| json!({"x":v[0],"y":v[1],"z":v[2]}))).collect();
            let mut mean = [0.0; 3];
            mean[axis] = inward * inflow / area;
            replace(
                "velocity_boundaries",
                Some(json!({"patch":name,"kind":"fixed_value","value":{"x":mean[0],"y":mean[1],"z":mean[2]}})),
                &mut child,
            )?;
            replace("pressure_boundaries", None, &mut child)?;
            profiles.push(json!({"patch":name,"velocity":list}));
            transferred.push(((*name).to_string(), "velocity".to_string()));
        } else {
            let list: Vec<f64> = pressures.iter().map(|p| p.unwrap_or(mean_pressure)).collect();
            replace("velocity_boundaries", Some(json!({"patch":name,"kind":"zero_gradient"})), &mut child)?;
            replace("pressure_boundaries", Some(json!({"patch":name,"kind":"fixed_value","value":mean_pressure})), &mut child)?;
            profiles.push(json!({"patch":name,"pressure":list}));
            transferred.push(((*name).to_string(), "pressure".to_string()));
        }
    }
    let object = child.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
    object.insert("transferred".into(), json!(profiles));
    let mesh = object.get_mut("mesh").and_then(Value::as_object_mut).ok_or_else(|| bad("the case has no mesh"))?;
    mesh.insert("dimensions".into(), json!(child_dims));
    mesh.insert("lengths".into(), json!({"x": extent[0], "y": extent[1], "z": extent[2]}));
    let even = child_nodes.iter().all(|list| {
        let n = list.len() - 1;
        let total = list[n];
        list.iter().enumerate().all(|(i, x)| (x - total * i as f64 / n as f64).abs() <= 1e-12 * total)
    });
    if even {
        mesh.remove("coordinates");
    } else {
        let [x, y, z] = child_nodes;
        mesh.insert("coordinates".into(), json!({"x": x, "y": y, "z": z}));
    }
    if transferred.iter().all(|(_, kind)| kind == "velocity")
        && !child.get("pressure_boundaries").and_then(Value::as_array).is_some_and(|r| r.iter().any(|row| row["kind"] == "fixed_value"))
    {
        notes.push("no face of the zoom box fixes the pressure; widen the box to include an outlet".into());
    }
    notes.push("temperatures are not transferred: the child keeps the parent's inflow temperature".into());
    let edit = finish(child, Some("/transferred".into()))?;
    Ok(Zoom { case: edit.case, offset_m: offset, transferred, notes, issues: edit.issues })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vtk::{Association, Field};

    /// The parent's results on its own grid: point data u = (a + b y, 0, 0), p = −g x.
    fn parent_results(case: &Value, a: f64, b: f64, g: f64) -> Grid {
        let [xs, ys, zs] = nodes(case).unwrap_or_default();
        let (nx, ny) = (xs.len(), ys.len());
        let mut points = Vec::new();
        for z in &zs {
            for y in &ys {
                for x in &xs {
                    points.push([*x, *y, *z]);
                }
            }
        }
        let id = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;
        let (mut connectivity, mut offsets) = (Vec::new(), Vec::new());
        for k in 0..zs.len() - 1 {
            for j in 0..ny - 1 {
                for i in 0..nx - 1 {
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
        let fields = vec![
            Field {
                name: VELOCITY_FIELD.into(),
                association: Association::Point,
                data_type: "Float64".into(),
                components: 3,
                values: points.iter().flat_map(|p| [a + b * p[1], 0.0, 0.0]).collect(),
            },
            Field {
                name: PRESSURE_FIELD.into(),
                association: Association::Point,
                data_type: "Float64".into(),
                components: 1,
                values: points.iter().map(|p| -g * p[0]).collect(),
            },
        ];
        Grid { points, connectivity, offsets, cell_types: vec![12; cells], fields }
    }

    fn channel() -> Value {
        serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null)
    }

    #[test]
    fn a_zoom_box_feeds_its_inner_faces_from_the_parent() {
        // Channel 1.0 × 0.5 × 0.1 m in 8 × 6 × 2 cells; zoom x 0.3..0.7 (snaps to 0.25..0.75).
        let parent = channel();
        let results = parent_results(&parent, 0.05, 0.2, 0.01);
        let zoom = zoom_case(&parent, &results, [0.3, 0.0, 0.0], [0.7, 0.5, 0.1], 1).expect("zoom");
        assert!(zoom.issues.is_empty(), "{:?}", zoom.issues);
        assert_eq!(zoom.offset_m, [0.25, 0.0, 0.0]);
        assert_eq!(zoom.case["mesh"]["dimensions"], json!([4, 6, 2]));
        assert!((zoom.case["mesh"]["lengths"]["x"].as_f64().unwrap_or(0.0) - 0.5).abs() < 1e-12);
        assert_eq!(zoom.transferred, vec![("left".to_string(), "velocity".to_string()), ("right".to_string(), "pressure".to_string())]);
        // Left: velocity at each (y, z) face centre is a + b y exactly (linear field).
        let rows = zoom.case["transferred"].as_array().cloned().unwrap_or_default();
        let left = rows[0]["velocity"].as_array().cloned().unwrap_or_default();
        assert_eq!(left.len(), 6 * 2);
        for (index, value) in left.iter().enumerate() {
            let y = (index % 6) as f64 * 0.5 / 6.0 + 0.5 / 12.0;
            assert!((value["x"].as_f64().unwrap_or(0.0) - (0.05 + 0.2 * y)).abs() < 1e-12, "{index}: {value}");
        }
        // Right: pressure −g x at the parent's x = 0.75.
        let right = rows[1]["pressure"].as_array().cloned().unwrap_or_default();
        assert!(right.iter().all(|p| (p.as_f64().unwrap_or(0.0) + 0.0075).abs() < 1e-12));
        // The walls and slip faces on the parent's boundary keep their conditions.
        let roles = crate::case_scene::patch_roles(&zoom.case);
        assert_eq!(roles[2], crate::case_scene::PatchRole::Wall);
        assert_eq!(roles[4], crate::case_scene::PatchRole::Symmetry);
        assert!(crate::airflow::requirements(&zoom.case).iter().any(|r| r.capability == "transferred_bc_v1"));
        // Refining doubles the cells and writes graded-free coordinates only when needed.
        let fine = zoom_case(&parent, &results, [0.3, 0.0, 0.0], [0.7, 0.5, 0.1], 2).expect("fine");
        assert!(fine.issues.is_empty(), "{:?}", fine.issues);
        assert_eq!(fine.case["mesh"]["dimensions"], json!([8, 12, 4]));
        assert!(fine.case.pointer("/mesh/coordinates").is_none());
        assert_eq!(fine.case["transferred"][0]["velocity"].as_array().map(Vec::len), Some(12 * 4));
    }

    #[test]
    fn zoom_boxes_clip_parts_and_refuse_what_cannot_work() {
        let mut parent = channel();
        parent["obstacles"] = json!([{"name":"inside","min":[3,0,0],"max":[5,2,2]},{"name":"across","min":[1,4,0],"max":[3,6,2]},{"name":"outside","min":[7,0,0],"max":[8,1,2]}]);
        let results = parent_results(&parent, 0.05, 0.0, 0.01);
        let zoom = zoom_case(&parent, &results, [0.25, 0.0, 0.0], [0.75, 0.5, 0.1], 1).expect("zoom");
        assert!(zoom.issues.is_empty(), "{:?}", zoom.issues);
        let obstacles = zoom.case["obstacles"].as_array().cloned().unwrap_or_default();
        assert_eq!(obstacles.len(), 2);
        assert_eq!(obstacles[0], json!({"name":"inside","min":[1,0,0],"max":[3,2,2]}));
        assert_eq!(obstacles[1], json!({"name":"across","min":[0,4,0],"max":[1,6,2]}));
        assert!(zoom.notes.iter().any(|n| n.contains("outside")));
        // The parent's results are a solid on a face: refused. The whole domain: refused.
        let mut walled = parent.clone();
        walled["obstacles"] = json!([{"name":"plug","min":[1,0,0],"max":[3,6,2]}]);
        let mut empty = parent_results(&walled, 0.05, 0.0, 0.01);
        empty.fields.clear();
        assert!(zoom_case(&walled, &empty, [0.25, 0.0, 0.0], [0.75, 0.5, 0.1], 1).is_err());
        assert!(zoom_case(&parent, &results, [-1.0, -1.0, -1.0], [2.0, 2.0, 2.0], 1).is_err());
        assert!(zoom_case(&parent, &results, [0.25, 0.0, 0.0], [0.75, 0.5, 0.1], 5).is_err());
    }
}
