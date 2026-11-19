//! Fluid volume check for `collocated_flow3d` cases: which grid cells the solver treats as
//! fluid and which as solid, how the fluid splits into connected regions, and which regions
//! no inlet, outlet or fan reaches (sealed pockets).
//!
//! Cells are classified as the companion meshes them. Grid obstacles take their half-open cell
//! blocks. A surface body takes the cells whose centre it encloses, with the companion's
//! winding rule (the same crossings, the same top-left edge rule), so the counts here match
//! the solver's mesh report. Porous blocks stay fluid. Regions are joined across cell faces
//! only, as the finite-volume fluxes are.
use schemars::JsonSchema;
use serde::Serialize;
use serde_json::Value;
use std::collections::BTreeMap;

use crate::case_scene::{PATCHES, PatchRole, grid, index3, nodes, patch_roles};
use crate::{Error, MAX_CASE_BYTES, Result, validate_json};

/// Most cells the check classifies (one byte and one label per cell).
pub const MAX_CHECK_CELLS: usize = 16_000_000;
/// Most blocks filling the sealed pockets may add.
pub const MAX_FILL_BLOCKS: usize = 256;
/// Most fluid regions listed one by one; the rest are summed.
pub const MAX_LISTED_REGIONS: usize = 64;

/// One connected fluid region.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct FluidRegion {
    pub cells: usize,
    pub volume_m3: f64,
    /// Bounds of its cells in metres.
    pub min: [f64; 3],
    pub max: [f64; 3],
    /// The inlet, outlet and fan patches it touches (patch names).
    pub openings: Vec<String>,
}

impl FluidRegion {
    /// No inlet, outlet or fan reaches it.
    pub fn sealed(&self) -> bool {
        self.openings.is_empty()
    }
}

/// A solid part and how many cells it takes.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct SolidPart {
    pub name: String,
    pub cells: usize,
    pub volume_m3: f64,
    /// Cells it shares with an earlier part (counted once, under the earlier part).
    pub overlapping_cells: usize,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct FluidVolume {
    pub total_cells: usize,
    pub fluid_cells: usize,
    pub solid_cells: usize,
    pub total_volume_m3: f64,
    pub fluid_volume_m3: f64,
    pub solid_volume_m3: f64,
    /// Patches that let flow in or out: inlets, outlets and the fan.
    pub openings: Vec<String>,
    /// Connected fluid regions, largest first (at most [`MAX_LISTED_REGIONS`]).
    pub regions: Vec<FluidRegion>,
    /// Regions beyond the listed ones, and their cells.
    pub unlisted_regions: usize,
    pub unlisted_cells: usize,
    /// Solid parts in case order.
    pub solids: Vec<SolidPart>,
    /// Surface bodies that enclose no cell centre (too thin for the grid).
    pub empty_bodies: Vec<String>,
}

impl FluidVolume {
    /// Sealed regions among the listed ones.
    pub fn sealed(&self) -> impl Iterator<Item = &FluidRegion> {
        self.regions.iter().filter(|r| r.sealed())
    }
}

fn bad(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

/// Cell centres along an axis.
fn centres(nodes: &[f64]) -> Vec<f64> {
    nodes.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect()
}

/// Cells whose centre a closed triangle surface encloses (non-zero winding along +x), as the
/// companion's `SurfaceBody::voxelize` decides them. Index `(k * ny + j) * nx + i`.
pub fn voxelize_body(vertices: &[[f64; 3]], triangles: &[[usize; 3]], cx: &[f64], cy: &[f64], cz: &[f64]) -> Vec<bool> {
    let (nx, ny, nz) = (cx.len(), cy.len(), cz.len());
    let span = |c: &[f64], lo: f64, hi: f64| -> Option<(usize, usize)> {
        let first = c.partition_point(|v| *v < lo);
        let end = c.partition_point(|v| *v <= hi);
        (first < end).then(|| (first, end - 1))
    };
    let mut inside = vec![false; nx * ny * nz];
    let mut columns: BTreeMap<(usize, usize), Vec<(f64, i32)>> = BTreeMap::new();
    for t in triangles {
        let (Some(a), Some(b), Some(c)) = (vertices.get(t[0]), vertices.get(t[1]), vertices.get(t[2])) else { continue };
        let (a, b, c) = (*a, *b, *c);
        let area = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]);
        if area == 0.0 {
            continue;
        }
        let (ymin, ymax) = (a[1].min(b[1]).min(c[1]), a[1].max(b[1]).max(c[1]));
        let (zmin, zmax) = (a[2].min(b[2]).min(c[2]), a[2].max(b[2]).max(c[2]));
        let (Some((j0, j1)), Some((k0, k1))) = (span(cy, ymin, ymax), span(cz, zmin, zmax)) else { continue };
        let sign = if area > 0.0 { 1.0 } else { -1.0 };
        for k in k0..=k1 {
            for j in j0..=j1 {
                let (py, pz) = (cy[j], cz[k]);
                // Evaluated with the endpoints in a fixed order, so the two triangles sharing an
                // edge get exactly opposite values and a centre on it is counted once.
                let edge = |p: [f64; 3], q: [f64; 3]| {
                    let swap = (q[1], q[2]) < (p[1], p[2]);
                    let (p, q) = if swap { (q, p) } else { (p, q) };
                    let w = (q[1] - p[1]) * (pz - p[2]) - (q[2] - p[2]) * (py - p[1]);
                    sign * if swap { -w } else { w }
                };
                let owns = |p: [f64; 3], q: [f64; 3]| {
                    let (ey, ez) = (sign * (q[1] - p[1]), sign * (q[2] - p[2]));
                    ez > 0.0 || (ez == 0.0 && ey < 0.0)
                };
                let covered = [(a, b), (b, c), (c, a)].iter().all(|(p, q)| {
                    let w = edge(*p, *q);
                    w > 0.0 || (w == 0.0 && owns(*p, *q))
                });
                if !covered {
                    continue;
                }
                let w0 = (b[1] - py) * (c[2] - pz) - (b[2] - pz) * (c[1] - py);
                let w1 = (c[1] - py) * (a[2] - pz) - (c[2] - pz) * (a[1] - py);
                let x = (w0 * a[0] + w1 * b[0] + (area - w0 - w1) * c[0]) / area;
                columns.entry((j, k)).or_default().push((x, if area > 0.0 { -1 } else { 1 }));
            }
        }
    }
    for ((j, k), mut crossings) in columns {
        crossings.sort_by(|p, q| p.0.total_cmp(&q.0));
        let mut winding = 0;
        let mut next = crossings.iter().peekable();
        for (i, &x) in cx.iter().enumerate() {
            while let Some((at, w)) = next.peek() {
                if *at <= x {
                    winding += w;
                    next.next();
                } else {
                    break;
                }
            }
            if winding != 0
                && let Some(cell) = inside.get_mut((k * ny + j) * nx + i)
            {
                *cell = true;
            }
        }
    }
    inside
}

/// Classify the cells of a `collocated_flow3d` case and split the fluid into regions.
pub fn fluid_volume(case: &Value) -> Result<FluidVolume> {
    analyse(case).map(|(report, _)| report)
}

/// The report, and per cell the position of its region in `report.regions` (`u32::MAX` for
/// solid cells and unlisted regions).
fn analyse(case: &Value) -> Result<(FluidVolume, Vec<u32>)> {
    validate_json(case, MAX_CASE_BYTES)?;
    let (dimensions, _) = grid(case)?;
    let [nx, ny, nz] = dimensions;
    let total = nx.checked_mul(ny).and_then(|v| v.checked_mul(nz)).filter(|n| *n <= MAX_CHECK_CELLS);
    let Some(total) = total else {
        return Err(Error::Limit(format!("the fluid volume check takes at most {MAX_CHECK_CELLS} cells")));
    };
    let grid_nodes = nodes(case)?;
    let width = |axis: usize, i: usize| grid_nodes[axis][i + 1] - grid_nodes[axis][i];
    let volume = |i: usize, j: usize, k: usize| width(0, i) * width(1, j) * width(2, k);
    let index = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;
    // 0 fluid, otherwise 1 + the solid part's position.
    let mut owner = vec![0u32; total];
    let mut solids = Vec::new();
    let mut empty_bodies = Vec::new();
    let claim = |owner: &mut [u32], name: &str, cells: &mut dyn Iterator<Item = (usize, usize, usize)>, solids: &mut Vec<SolidPart>| {
        let id = solids.len() as u32 + 1;
        let mut part = SolidPart { name: name.into(), cells: 0, volume_m3: 0.0, overlapping_cells: 0 };
        for (i, j, k) in cells {
            let slot = &mut owner[index(i, j, k)];
            if *slot == 0 {
                *slot = id;
                part.cells += 1;
                part.volume_m3 += volume(i, j, k);
            } else {
                part.overlapping_cells += 1;
            }
        }
        solids.push(part);
    };
    for row in case.get("obstacles").and_then(Value::as_array).into_iter().flatten() {
        let (Some(name), Some(min), Some(max)) = (row["name"].as_str(), index3(&row["min"]), index3(&row["max"])) else { continue };
        let hi = [max[0].min(nx), max[1].min(ny), max[2].min(nz)];
        let mut cells = (min[2]..hi[2]).flat_map(|k| (min[1]..hi[1]).flat_map(move |j| (min[0]..hi[0]).map(move |i| (i, j, k))));
        claim(&mut owner, name, &mut cells, &mut solids);
    }
    let (cx, cy, cz) = (centres(&grid_nodes[0]), centres(&grid_nodes[1]), centres(&grid_nodes[2]));
    for row in case.get("bodies").and_then(Value::as_array).into_iter().flatten() {
        let Some(name) = row["name"].as_str() else { continue };
        let vertices: Vec<[f64; 3]> =
            row["vertices"].as_array().into_iter().flatten().filter_map(|v| serde_json::from_value::<[f64; 3]>(v.clone()).ok()).collect();
        let triangles: Vec<[usize; 3]> = row["triangles"].as_array().into_iter().flatten().filter_map(index3).collect();
        let inside = voxelize_body(&vertices, &triangles, &cx, &cy, &cz);
        let mut cells = inside.iter().enumerate().filter(|(_, v)| **v).map(|(c, _)| (c % nx, (c / nx) % ny, c / (nx * ny)));
        claim(&mut owner, name, &mut cells, &mut solids);
        if solids.last().is_some_and(|p| p.cells + p.overlapping_cells == 0) {
            empty_bodies.push(name.to_string());
        }
    }
    // Openings: which patches let flow in or out, and the cells on each.
    let roles = patch_roles(case);
    let open: Vec<usize> = (0..6).filter(|p| matches!(roles[*p], PatchRole::Inlet | PatchRole::Outlet | PatchRole::Fan)).collect();
    let on_patch = |patch: usize, i: usize, j: usize, k: usize| {
        let at = [i, j, k];
        let axis = patch / 2;
        if patch % 2 == 0 { at[axis] == 0 } else { at[axis] + 1 == dimensions[axis] }
    };
    // Flood fill the fluid (face neighbours only).
    let mut label = vec![u32::MAX; total];
    let mut regions: Vec<FluidRegion> = Vec::new();
    let mut stack = Vec::new();
    let mut fluid_cells = 0usize;
    let mut fluid_volume_m3 = 0.0;
    for start in 0..total {
        if owner[start] != 0 || label[start] != u32::MAX {
            continue;
        }
        let id = regions.len() as u32;
        let mut region = FluidRegion { cells: 0, volume_m3: 0.0, min: [f64::INFINITY; 3], max: [f64::NEG_INFINITY; 3], openings: Vec::new() };
        let mut touched = [false; 6];
        label[start] = id;
        stack.push(start);
        while let Some(c) = stack.pop() {
            let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
            region.cells += 1;
            region.volume_m3 += volume(i, j, k);
            let at = [i, j, k];
            for a in 0..3 {
                region.min[a] = region.min[a].min(grid_nodes[a][at[a]]);
                region.max[a] = region.max[a].max(grid_nodes[a][at[a] + 1]);
            }
            for &p in &open {
                if on_patch(p, i, j, k) {
                    touched[p] = true;
                }
            }
            let mut visit = |n: usize| {
                if owner[n] == 0 && label[n] == u32::MAX {
                    label[n] = id;
                    stack.push(n);
                }
            };
            if i > 0 {
                visit(c - 1);
            }
            if i + 1 < nx {
                visit(c + 1);
            }
            if j > 0 {
                visit(c - nx);
            }
            if j + 1 < ny {
                visit(c + nx);
            }
            if k > 0 {
                visit(c - nx * ny);
            }
            if k + 1 < nz {
                visit(c + nx * ny);
            }
        }
        region.openings = (0..6).filter(|p| touched[*p]).map(|p| PATCHES[p].to_string()).collect();
        fluid_cells += region.cells;
        fluid_volume_m3 += region.volume_m3;
        regions.push(region);
    }
    let mut order: Vec<usize> = (0..regions.len()).collect();
    order.sort_by(|a, b| regions[*b].cells.cmp(&regions[*a].cells).then(regions[*a].min[0].total_cmp(&regions[*b].min[0])));
    // Old label -> position in the sorted list, if listed.
    let mut position = vec![u32::MAX; regions.len()];
    for (rank, old) in order.iter().enumerate().take(MAX_LISTED_REGIONS) {
        position[*old] = rank as u32;
    }
    for l in label.iter_mut() {
        *l = position.get(*l as usize).copied().unwrap_or(u32::MAX);
    }
    let mut sorted: Vec<FluidRegion> = order.iter().map(|i| regions[*i].clone()).collect();
    let listed = sorted.len().min(MAX_LISTED_REGIONS);
    let rest = sorted.split_off(listed);
    let regions = sorted;
    let total_volume_m3: f64 = (0..3).map(|a| grid_nodes[a].last().copied().unwrap_or(0.0) - grid_nodes[a][0]).product();
    let solid_cells = total - fluid_cells;
    let report = FluidVolume {
        total_cells: total,
        fluid_cells,
        solid_cells,
        total_volume_m3,
        fluid_volume_m3,
        solid_volume_m3: solids.iter().map(|p| p.volume_m3).sum(),
        openings: open.iter().map(|p| PATCHES[*p].to_string()).collect(),
        unlisted_regions: rest.len(),
        unlisted_cells: rest.iter().map(|r| r.cells).sum(),
        regions,
        solids,
        empty_bodies,
    };
    Ok((report, label))
}

/// Fill the sealed fluid pockets with solid: each pocket becomes grid obstacles named
/// `pocket_<n>_<m>` (its cells merged into boxes), so the solver no longer carries dead
/// fluid the inlets cannot reach. Fluid regions that reach an opening are left alone.
pub fn fill_sealed_pockets(case: &Value) -> Result<crate::case_scene::Edit> {
    let (report, label) = analyse(case)?;
    if report.openings.is_empty() {
        return Err(bad("the case has no inlet, outlet or fan, so every fluid region is sealed; set the openings first"));
    }
    if report.unlisted_regions > 0 {
        return Err(bad(format!("more than {MAX_LISTED_REGIONS} fluid regions; coarsen the parts or the grid first")));
    }
    if report.sealed().next().is_none() {
        return Err(bad("there are no sealed fluid pockets"));
    }
    let (dimensions, _) = grid(case)?;
    let [nx, ny, nz] = dimensions;
    let mut obstacles = case.get("obstacles").and_then(Value::as_array).cloned().unwrap_or_default();
    let existing: Vec<String> = obstacles.iter().filter_map(|r| r["name"].as_str().map(str::to_string)).collect();
    let mut added = 0usize;
    for (pocket, (n, _)) in report.regions.iter().enumerate().filter(|(_, r)| r.sealed()).enumerate() {
        let id = n as u32;
        // Greedy box merge: grow along x, then y, then z over the pocket's unclaimed cells.
        let mut mine: Vec<bool> = label.iter().map(|l| *l == id).collect();
        let at = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;
        let mut boxes = 0usize;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if !mine[at(i, j, k)] {
                        continue;
                    }
                    let mut i1 = i + 1;
                    while i1 < nx && mine[at(i1, j, k)] {
                        i1 += 1;
                    }
                    let mut j1 = j + 1;
                    while j1 < ny && (i..i1).all(|x| mine[at(x, j1, k)]) {
                        j1 += 1;
                    }
                    let mut k1 = k + 1;
                    while k1 < nz && (j..j1).all(|y| (i..i1).all(|x| mine[at(x, y, k1)])) {
                        k1 += 1;
                    }
                    for z in k..k1 {
                        for y in j..j1 {
                            for x in i..i1 {
                                mine[at(x, y, z)] = false;
                            }
                        }
                    }
                    boxes += 1;
                    added += 1;
                    if added > MAX_FILL_BLOCKS {
                        return Err(Error::Limit(format!("filling the pockets needs more than {MAX_FILL_BLOCKS} blocks")));
                    }
                    let mut name = format!("pocket_{}_{}", pocket + 1, boxes);
                    while existing.contains(&name) {
                        name.push('_');
                    }
                    obstacles.push(serde_json::json!({"name": name, "min": [i, j, k], "max": [i1, j1, k1]}));
                }
            }
        }
    }
    let mut proposed = case.clone();
    let object = proposed.as_object_mut().ok_or_else(|| bad("the case must be an object"))?;
    let first = obstacles.len() - added;
    object.insert("obstacles".into(), Value::Array(obstacles));
    crate::case_scene::finish(proposed, Some(format!("/obstacles/{first}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// A 10 × 6 × 1 channel, 10 mm cells: inlet on the left, outlet on the right, and a ring of
    /// four blocks around cells x 4..6, y 2..4 that seals a 2 × 2 pocket.
    fn ring_case() -> Value {
        json!({
            "solver": "collocated_flow3d",
            "mesh": {"dimensions": [10, 6, 1], "lengths": {"x": 0.1, "y": 0.06, "z": 0.01}},
            "config": {"dt": 0.001},
            "velocity_boundaries": [
                {"patch": "left", "kind": "fixed_value", "value": {"x": 0.5, "y": 0.0, "z": 0.0}},
                {"patch": "bottom", "kind": "fixed_value", "value": {"x": 0.0, "y": 0.0, "z": 0.0}},
                {"patch": "top", "kind": "fixed_value", "value": {"x": 0.0, "y": 0.0, "z": 0.0}}
            ],
            "pressure_boundaries": [{"patch": "right", "kind": "fixed_value", "value": 0.0}],
            "obstacles": [
                {"name": "ring_bottom", "min": [3, 1, 0], "max": [7, 2, 1]},
                {"name": "ring_top", "min": [3, 4, 0], "max": [7, 5, 1]},
                {"name": "ring_left", "min": [3, 2, 0], "max": [4, 4, 1]},
                {"name": "ring_right", "min": [6, 2, 0], "max": [7, 4, 1]}
            ]
        })
    }

    #[test]
    fn a_ring_of_blocks_seals_a_pocket_and_filling_it_leaves_one_region() {
        let case = ring_case();
        let report = fluid_volume(&case).expect("report");
        assert_eq!((report.total_cells, report.solid_cells, report.fluid_cells), (60, 12, 48));
        let cell = 0.01 * 0.01 * 0.01;
        assert!((report.solid_volume_m3 - 12.0 * cell).abs() < 1e-15);
        assert!((report.fluid_volume_m3 + report.solid_volume_m3 - report.total_volume_m3).abs() < 1e-15);
        assert_eq!(report.openings, vec!["left".to_string(), "right".to_string()]);
        assert_eq!(report.regions.len(), 2);
        assert_eq!((report.regions[0].cells, report.regions[0].openings.clone()), (44, vec!["left".to_string(), "right".to_string()]));
        let pocket = &report.regions[1];
        assert!(pocket.sealed());
        assert_eq!(pocket.cells, 4);
        assert!((pocket.min[0] - 0.04).abs() < 1e-12 && (pocket.max[0] - 0.06).abs() < 1e-12);
        assert!((pocket.min[1] - 0.02).abs() < 1e-12 && (pocket.max[1] - 0.04).abs() < 1e-12);
        assert_eq!(report.sealed().count(), 1);
        // Filling turns the pocket into one block, and the check then finds one open region.
        let filled = fill_sealed_pockets(&case).expect("fill");
        assert_eq!(filled.case["obstacles"][4], json!({"name": "pocket_1_1", "min": [4, 2, 0], "max": [6, 4, 1]}));
        assert_eq!(filled.pointer.as_deref(), Some("/obstacles/4"));
        let after = fluid_volume(&filled.case).expect("after");
        assert_eq!((after.regions.len(), after.solid_cells, after.sealed().count()), (1, 16, 0));
        assert!(fill_sealed_pockets(&filled.case).is_err(), "nothing left to fill");
    }

    #[test]
    fn an_l_shaped_pocket_fills_with_two_blocks() {
        let mut case = ring_case();
        // Widen the ring so the pocket is x 4..7, y 2..4 minus the corner (6, 3): an L of 5 cells.
        case["mesh"]["dimensions"] = json!([12, 6, 1]);
        case["mesh"]["lengths"]["x"] = json!(0.12);
        case["obstacles"] = json!([
            {"name": "b", "min": [3, 1, 0], "max": [8, 2, 1]},
            {"name": "t", "min": [3, 4, 0], "max": [8, 5, 1]},
            {"name": "l", "min": [3, 2, 0], "max": [4, 4, 1]},
            {"name": "r", "min": [7, 2, 0], "max": [8, 4, 1]},
            {"name": "corner", "min": [6, 3, 0], "max": [7, 4, 1]}
        ]);
        let report = fluid_volume(&case).expect("report");
        assert_eq!(report.sealed().map(|r| r.cells).collect::<Vec<_>>(), vec![5]);
        let filled = fill_sealed_pockets(&case).expect("fill");
        let added: Vec<_> = filled.case["obstacles"].as_array().map(|a| a[5..].to_vec()).unwrap_or_default();
        assert_eq!(
            added,
            vec![
                json!({"name": "pocket_1_1", "min": [4, 2, 0], "max": [7, 3, 1]}),
                json!({"name": "pocket_1_2", "min": [4, 3, 0], "max": [6, 4, 1]}),
            ]
        );
        assert_eq!(fluid_volume(&filled.case).expect("after").regions.len(), 1);
    }

    /// An axis-aligned box as an outward closed triangle surface.
    fn cube(min: [f64; 3], max: [f64; 3]) -> Value {
        let v: Vec<[f64; 3]> = (0..8).map(|i| [0, 1, 2].map(|a| if i >> a & 1 == 0 { min[a] } else { max[a] })).collect();
        let quads = [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]];
        let triangles: Vec<[usize; 3]> = quads.iter().flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]).collect();
        json!({"name": "box", "vertices": v, "triangles": triangles})
    }

    #[test]
    fn surface_bodies_take_the_cells_whose_centre_they_enclose_and_overlaps_count_once() {
        let mut case = ring_case();
        case["obstacles"] = json!([{"name": "block", "min": [1, 1, 0], "max": [3, 3, 1]}]);
        // Centres at 5, 15, 25 mm...: x 12..38 mm encloses centres 15, 25, 35 and y 12..28 mm
        // centres 15, 25: 6 cells, 4 of them already the block's (cells 1..3 × 1..3).
        let mut body = cube([0.012, 0.012, -0.001], [0.038, 0.028, 0.011]);
        body["name"] = json!("chip");
        let mut thin = cube([0.081, 0.011, 0.0], [0.084, 0.019, 0.01]);
        thin["name"] = json!("foil");
        case["bodies"] = json!([body, thin]);
        let report = fluid_volume(&case).expect("report");
        let chip = &report.solids[1];
        assert_eq!((chip.name.as_str(), chip.cells, chip.overlapping_cells), ("chip", 2, 4));
        assert_eq!(report.solid_cells, 4 + 2);
        assert_eq!(report.empty_bodies, vec!["foil".to_string()]);
        // The winding agrees with the column test for a cube covering exactly cells 5..7 in x.
        let x: Vec<f64> = (0..10).map(|i| 0.005 + 0.01 * i as f64).collect();
        let one = [0.005];
        let v: Vec<[f64; 3]> = cube([0.05, 0.0, 0.0], [0.07, 0.01, 0.01])["vertices"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|p| serde_json::from_value(p.clone()).ok())
            .collect();
        let t: Vec<[usize; 3]> = cube([0.0; 3], [1.0; 3])["triangles"].as_array().into_iter().flatten().filter_map(index3).collect();
        let inside = voxelize_body(&v, &t, &x, &one, &one);
        assert_eq!(inside.iter().enumerate().filter(|(_, b)| **b).map(|(i, _)| i).collect::<Vec<_>>(), vec![5, 6]);
    }

    #[test]
    fn closed_boxes_and_bad_cases_are_refused_or_reported() {
        let mut case = ring_case();
        case["velocity_boundaries"] = json!([]);
        case["pressure_boundaries"] = json!([]);
        let report = fluid_volume(&case).expect("report");
        assert!(report.openings.is_empty() && report.regions.iter().all(FluidRegion::sealed));
        assert!(fill_sealed_pockets(&case).is_err(), "with no openings everything is sealed");
        assert!(fluid_volume(&json!({"solver": "lbm_d2q9"})).is_err());
        let mut huge = ring_case();
        huge["mesh"]["dimensions"] = json!([1024, 1024, 1024]);
        assert!(fluid_volume(&huge).is_err());
    }

    /// Generalized winding number of a closed outward surface at `p` (1 inside, 0 outside).
    fn winding(v: &[[f64; 3]], t: &[[usize; 3]], p: [f64; 3]) -> f64 {
        let mut w = 0.0;
        for tri in t {
            let [a, b, c] = tri.map(|i| [v[i][0] - p[0], v[i][1] - p[1], v[i][2] - p[2]]);
            let n = |u: [f64; 3]| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
            let dot = |u: [f64; 3], q: [f64; 3]| u[0] * q[0] + u[1] * q[1] + u[2] * q[2];
            let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
            let den = n(a) * n(b) * n(c) + dot(a, b) * n(c) + dot(b, c) * n(a) + dot(c, a) * n(b);
            w += 2.0 * det.atan2(den);
        }
        w / (4.0 * std::f64::consts::PI)
    }

    #[test]
    fn voxels_match_the_exact_winding_number_with_centres_on_shared_edges() {
        // A prism whose side face is split along y + z = 0.01, through rows of cell centres,
        // and a tilted tetrahedron; the companion's voxelizer gives the same 540 cells.
        let mut v = vec![[0.03, 0.0, 0.0], [0.05, 0.0, 0.0], [0.03, 0.01, 0.0], [0.03, 0.0, 0.01], [0.05, 0.0, 0.01], [0.03, 0.01, 0.01]];
        let mut t = vec![[0, 2, 1], [3, 4, 5], [0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4], [2, 0, 3], [2, 3, 5]];
        let o = v.len();
        v.extend([[0.06, 0.002, 0.001], [0.078, 0.004, 0.002], [0.065, 0.018, 0.003], [0.068, 0.007, 0.0095]]);
        t.extend([[o, o + 2, o + 1], [o, o + 1, o + 3], [o + 1, o + 2, o + 3], [o, o + 3, o + 2]]);
        let c = |n: usize, l: f64| (0..n).map(|i| (i as f64 + 0.5) * l / n as f64).collect::<Vec<_>>();
        let (cx, cy, cz) = (c(40, 0.1), c(20, 0.02), c(10, 0.01));
        let inside = voxelize_body(&v, &t, &cx, &cy, &cz);
        for k in 0..10 {
            for j in 0..20 {
                for i in 0..40 {
                    let exact = winding(&v, &t, [cx[i], cy[j], cz[k]]).abs() > 0.5;
                    assert_eq!(inside[(k * 20 + j) * 40 + i], exact, "cell ({i}, {j}, {k})");
                }
            }
        }
        assert_eq!(inside.iter().filter(|c| **c).count(), 540);
    }
}
