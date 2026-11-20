//! Bounded primitive-domain meshing. No CAD or curved-surface fidelity is implied.
use super::{CollocatedFlowCase, GridObstacle, PolyMeshSpec, RunManifest};
use crate::{core::Vec3, fvm::*, io::VtuMesh, require, Result};
use std::{collections::BTreeSet, path::Path};

const EXTERIOR: [&str; 6] = ["left", "right", "bottom", "top", "front", "back"];

/// Solid cells of the background grid: for each cell, the index of the solid that owns it
/// (boxes first, then surface bodies), and the solids' names in that order.
pub(super) struct SolidLayout {
    pub cell_solid: Vec<Option<usize>>,
    pub names: Vec<String>,
}

pub(super) fn layout(
    spec: &PolyMeshSpec,
    obstacles: &[GridObstacle],
    bodies: &[super::bodies::SurfaceBody],
) -> Result<SolidLayout> {
    require(
        obstacles.len() + bodies.len() <= 128,
        "at most 128 grid obstacles and surface bodies are supported",
    )?;
    require(spec.xy_shear == 0.0, "grid obstacles require zero xy_shear")?;
    let [nx, ny, nz] = spec.dimensions;
    let count = nx.checked_mul(ny).and_then(|v| v.checked_mul(nz));
    require(
        count.is_some_and(|n| n > 0 && n <= 32_768),
        "obstacle grid must contain 1..32768 cells",
    )?;
    require(
        spec.dimensions.iter().all(|&n| n > 0),
        "positive grid dimensions required",
    )?;
    let count =
        count.ok_or_else(|| crate::SolverError::InvalidInput("grid size overflow".into()))?;
    let mut names: BTreeSet<&str> = EXTERIOR.into_iter().collect();
    let mut solid = vec![None; count];
    for (oi, obstacle) in obstacles.iter().enumerate() {
        let name = obstacle.name.as_bytes();
        require(
            !name.is_empty() && name.len() <= 80
                && name.first().is_some_and(u8::is_ascii_alphanumeric)
                && name.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_' || *c == b'-')
                && names.insert(&obstacle.name),
            "obstacle names must be unique ASCII identifiers (1..80 characters), distinct from exterior patches",
        )?;
        require(
            (0..3).all(|a| {
                obstacle.min[a] < obstacle.max[a] && obstacle.max[a] <= spec.dimensions[a]
            }),
            "obstacle bounds must be nonempty half-open grid cell indices within dimensions",
        )?;
        for k in obstacle.min[2]..obstacle.max[2] {
            for j in obstacle.min[1]..obstacle.max[1] {
                for i in obstacle.min[0]..obstacle.max[0] {
                    let slot = solid.get_mut((k * ny + j) * nx + i).ok_or_else(|| {
                        crate::SolverError::InvalidInput("grid index out of range".into())
                    })?;
                    require(slot.is_none(), "grid obstacles must not overlap")?;
                    *slot = Some(oi);
                }
            }
        }
    }
    for (bi, body) in bodies.iter().enumerate() {
        let name = body.name.as_bytes();
        require(
            !name.is_empty() && name.len() <= 80
                && name.first().is_some_and(u8::is_ascii_alphanumeric)
                && name.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_' || *c == b'-')
                && names.insert(&body.name),
            "surface body names must be unique ASCII identifiers distinct from patches and obstacles",
        )?;
        let inside = body.voxelize(spec)?;
        require(
            inside.iter().any(|v| *v),
            format!(
                "surface body {} covers no cell centre; refine the grid",
                body.name
            ),
        )?;
        for (slot, inside) in solid.iter_mut().zip(inside) {
            if inside {
                require(
                    slot.is_none(),
                    format!("surface body {} overlaps another solid", body.name),
                )?;
                *slot = Some(obstacles.len() + bi);
            }
        }
    }
    let names = obstacles
        .iter()
        .map(|o| o.name.clone())
        .chain(bodies.iter().map(|b| b.name.clone()))
        .collect();
    Ok(SolidLayout {
        cell_solid: solid,
        names,
    })
}

pub(super) fn build(
    spec: &PolyMeshSpec,
    obstacles: &[GridObstacle],
    bodies: &[super::bodies::SurfaceBody],
    conditions: &[super::bodies::BodyFaceCondition],
    cut: bool,
) -> Result<(PolyMesh, VtuMesh)> {
    if cut {
        return build_cut(spec, obstacles, bodies, conditions);
    }
    require(
        bodies.len() >= 1 || conditions.is_empty(),
        "conditions on body faces need surface bodies",
    )?;
    if obstacles.is_empty() && bodies.is_empty() {
        return spec.build();
    }
    let SolidLayout {
        cell_solid: solid,
        names,
    } = layout(spec, obstacles, bodies)?;
    let count = solid.len();
    let (base, full_export) = spec.build()?;
    let mut retained = vec![None; count];
    let mut cells = Vec::new();
    let mut export_cells = Vec::new();
    for (i, (cell, hex)) in base.cells().iter().zip(&full_export.cells).enumerate() {
        if solid.get(i) == Some(&None) {
            retained[i] = Some(cells.len());
            cells.push(*cell);
            export_cells.push(hex.clone());
        }
    }
    require(!cells.is_empty(), "obstacles remove every fluid cell")?;
    let mut patches = base.patches().to_vec();
    patches.extend(
        names
            .iter()
            .map(|name| BoundaryPatch { name: name.clone() }),
    );
    // Conditions on body faces: each one is a patch after the solids', and owns the faces of
    // its body nearest its triangles. `owner_of[c]` is the triangle owner table of body `c`.
    let solids = patches.len();
    let mut owner_of: Vec<Option<Vec<Option<usize>>>> = vec![None; bodies.len()];
    {
        let mut taken: BTreeSet<&str> = names.iter().map(String::as_str).chain(EXTERIOR).collect();
        for (ci, condition) in conditions.iter().enumerate() {
            let name = condition.name.as_bytes();
            require(
                !name.is_empty() && name.len() <= 80
                    && name.first().is_some_and(u8::is_ascii_alphanumeric)
                    && name.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_' || *c == b'-')
                    && taken.insert(&condition.name),
                "body face condition names must be unique ASCII identifiers distinct from patches and solids",
            )?;
            let bi = bodies.iter().position(|b| b.name == condition.body).ok_or_else(|| {
                crate::SolverError::InvalidInput(format!(
                    "body face condition {} names no surface body {}",
                    condition.name, condition.body
                ))
            })?;
            let count = bodies[bi].triangles.len();
            let table = owner_of[bi].get_or_insert_with(|| vec![None; count]);
            require(
                !condition.triangles.is_empty(),
                format!("body face condition {} lists no triangles", condition.name),
            )?;
            for &t in &condition.triangles {
                let slot = table.get_mut(t).ok_or_else(|| {
                    crate::SolverError::InvalidInput(format!(
                        "body face condition {} lists triangle {t}, which {} does not have",
                        condition.name, condition.body
                    ))
                })?;
                require(
                    slot.is_none(),
                    format!("triangle {t} of {} is in two body face conditions", condition.body),
                )?;
                *slot = Some(ci);
            }
            patches.push(BoundaryPatch { name: condition.name.clone() });
        }
    }
    let mut faces = Vec::new();
    for f in base.faces() {
        let owner = retained.get(f.owner).copied().flatten();
        let neighbour = f.neighbour.and_then(|n| retained.get(n).copied().flatten());
        match (owner, neighbour, f.neighbour) {
            (Some(o), Some(n), _) => faces.push(Face::internal(o, n, f.center, f.area)),
            (Some(o), None, Some(n)) => {
                let obstacle = solid.get(n).copied().flatten().ok_or_else(|| {
                    crate::SolverError::InvalidInput("missing solid neighbour".into())
                })?;
                faces.push(Face::boundary_face(
                    o,
                    f.center,
                    f.area,
                    EXTERIOR.len() + obstacle,
                ));
            }
            (None, Some(n), _) => {
                let obstacle = solid.get(f.owner).copied().flatten().ok_or_else(|| {
                    crate::SolverError::InvalidInput("missing solid owner".into())
                })?;
                faces.push(Face::boundary_face(
                    n,
                    f.center,
                    -f.area,
                    EXTERIOR.len() + obstacle,
                ));
            }
            (Some(o), None, None) => {
                let patch = f.patch.ok_or_else(|| {
                    crate::SolverError::InvalidInput("missing exterior patch".into())
                })?;
                faces.push(Face::boundary_face(o, f.center, f.area, patch));
            }
            _ => {}
        }
    }
    if !conditions.is_empty() {
        // The body's fluid-side faces nearest a condition's triangles move to its patch.
        for (bi, table) in owner_of.iter().enumerate() {
            let Some(table) = table else { continue };
            let patch = EXTERIOR.len() + obstacles.len() + bi;
            let on_body: Vec<usize> = (0..faces.len()).filter(|&f| faces[f].patch == Some(patch)).collect();
            let centres: Vec<[f64; 3]> = on_body
                .iter()
                .map(|&f| [faces[f].center.x, faces[f].center.y, faces[f].center.z])
                .collect();
            let nearest = bodies[bi].nearest_triangles(&centres)?;
            for (&f, t) in on_body.iter().zip(nearest) {
                if let Some(Some(ci)) = table.get(t) {
                    faces[f].patch = Some(solids + ci);
                }
            }
        }
    }
    let mesh = PolyMesh::new(cells, faces, patches)?;
    validate_topology(&mesh)?;
    for p in EXTERIOR.len()..mesh.patches().len() {
        let what = if p >= solids {
            format!(
                "body face condition {} covers no fluid face of the grid; refine the grid near it",
                mesh.patches()[p].name
            )
        } else {
            "an obstacle has no exposed fluid boundary".to_string()
        };
        require(mesh.faces().iter().any(|f| f.patch == Some(p)), what)?;
    }
    let export = VtuMesh {
        points: full_export.points,
        cells: export_cells,
    };
    export.validate()?;
    require(
        export.cells.len() == mesh.cell_count(),
        "export/solver cell count mismatch",
    )?;
    Ok((mesh, export))
}

/// The condition tables of `conditions`: per body, which condition (if any) owns each of its
/// triangles, after checking names and triangle lists as [`build`] does.
fn condition_tables(
    bodies: &[super::bodies::SurfaceBody],
    conditions: &[super::bodies::BodyFaceCondition],
    taken: &mut BTreeSet<String>,
) -> Result<Vec<Option<Vec<Option<usize>>>>> {
    let mut owner_of: Vec<Option<Vec<Option<usize>>>> = vec![None; bodies.len()];
    for (ci, condition) in conditions.iter().enumerate() {
        let name = condition.name.as_bytes();
        require(
            !name.is_empty() && name.len() <= 80
                && name.first().is_some_and(u8::is_ascii_alphanumeric)
                && name.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_' || *c == b'-')
                && taken.insert(condition.name.clone()),
            "body face condition names must be unique ASCII identifiers distinct from patches and solids",
        )?;
        let bi = bodies.iter().position(|b| b.name == condition.body).ok_or_else(|| {
            crate::SolverError::InvalidInput(format!("body face condition {} names no surface body {}", condition.name, condition.body))
        })?;
        let table = owner_of[bi].get_or_insert_with(|| vec![None; bodies[bi].triangles.len()]);
        require(!condition.triangles.is_empty(), format!("body face condition {} lists no triangles", condition.name))?;
        for &t in &condition.triangles {
            let slot = table.get_mut(t).ok_or_else(|| {
                crate::SolverError::InvalidInput(format!(
                    "body face condition {} lists triangle {t}, which {} does not have",
                    condition.name, condition.body
                ))
            })?;
            require(slot.is_none(), format!("triangle {t} of {} is in two body face conditions", condition.body))?;
            *slot = Some(ci);
        }
    }
    Ok(owner_of)
}

/// The mesh with cut cells (`cut_cells_v1`, see `cut_cells`): every connected piece of fluid in
/// a cell is a cell of its own, grid faces keep their open part between pieces, each piece gets a
/// wall face that closes it, and small pieces are merged.
fn build_cut(
    spec: &PolyMeshSpec,
    obstacles: &[GridObstacle],
    bodies: &[super::bodies::SurfaceBody],
    conditions: &[super::bodies::BodyFaceCondition],
) -> Result<(PolyMesh, VtuMesh)> {
    require(!bodies.is_empty(), "cut cells need surface bodies")?;
    require(spec.xy_shear == 0.0, "cut cells require zero xy_shear")?;
    for body in bodies {
        body.validate()?;
    }
    // Grid obstacles stay whole cells, exactly as without cut cells.
    let blocks = layout(spec, obstacles, &[])?;
    let blocked: Vec<bool> = blocks.cell_solid.iter().map(Option::is_some).collect();
    let mut names = blocks.names.clone();
    names.extend(bodies.iter().map(|b| b.name.clone()));
    let mut taken: BTreeSet<String> = names.iter().cloned().chain(EXTERIOR.iter().map(|s| s.to_string())).collect();
    require(taken.len() == names.len() + EXTERIOR.len(), "solid names must be unique and distinct from exterior patches")?;
    let cut = super::cut_cells::Cut::new(spec, bodies, blocked.clone())?;
    let (base, full_export) = spec.build()?;
    let pieces = cut.fragments();
    let count = pieces.len();
    require(count > 0, "obstacles remove every fluid cell")?;
    let body_patch = |b: usize| EXTERIOR.len() + obstacles.len() + b;
    let axis_of = |f: &Face| -> usize {
        let a = [f.area.x.abs(), f.area.y.abs(), f.area.z.abs()];
        if a[0] >= a[1] && a[0] >= a[2] {
            0
        } else if a[1] >= a[2] {
            1
        } else {
            2
        }
    };
    let component = |v: Vec3, a: usize| match a {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    };
    // Open grid faces between pieces, or from a piece to the exterior or an obstacle.
    let mut faces: Vec<Face> = Vec::new();
    for f in base.faces() {
        let axis = axis_of(f);
        let owner_upper = component(f.center, axis) > component(base.cells()[f.owner].center, axis);
        match f.neighbour {
            Some(n) if !blocked[f.owner] => {
                for link in cut.links(f.owner, axis, owner_upper, !blocked[n]) {
                    let area = f.area * link.fraction;
                    match link.across {
                        Some(other) => faces.push(Face::internal(link.from, other, link.centre, area)),
                        None => {
                            let obstacle = blocks.cell_solid.get(n).copied().flatten().unwrap_or(0);
                            faces.push(Face::boundary_face(link.from, link.centre, area, EXTERIOR.len() + obstacle));
                        }
                    }
                }
            }
            Some(n) if !blocked[n] => {
                let obstacle = blocks.cell_solid.get(f.owner).copied().flatten().unwrap_or(0);
                for link in cut.links(n, axis, !owner_upper, false) {
                    faces.push(Face::boundary_face(link.from, link.centre, -f.area * link.fraction, EXTERIOR.len() + obstacle));
                }
            }
            Some(_) => {}
            None if !blocked[f.owner] => {
                let patch = f.patch.ok_or_else(|| crate::SolverError::InvalidInput("missing exterior patch".into()))?;
                for link in cut.links(f.owner, axis, owner_upper, false) {
                    faces.push(Face::boundary_face(link.from, link.centre, f.area * link.fraction, patch));
                }
            }
            None => {}
        }
    }
    // Each piece's wall closes it: minus the sum of its open faces' outward area vectors.
    let mut closure = vec![Vec3::ZERO; count];
    let mut scale = vec![0.0; count];
    for f in &faces {
        closure[f.owner] += f.area;
        scale[f.owner] += f.area.norm();
        if let Some(n) = f.neighbour {
            closure[n] -= f.area;
            scale[n] += f.area.norm();
        }
    }
    let condition_of = condition_tables(bodies, conditions, &mut taken)?;
    let mut patches = base.patches().to_vec();
    patches.extend(names.iter().map(|name| BoundaryPatch { name: name.clone() }));
    let solids = patches.len();
    patches.extend(conditions.iter().map(|c| BoundaryPatch { name: c.name.clone() }));
    let measured = cut.walls();
    let volume_of = |p: usize| pieces[p].fraction * base.cells()[pieces[p].cell].volume;
    let mut walls: Vec<(usize, Vec3, Vec3, Option<usize>, usize)> = Vec::new(); // piece, centre, area, condition, body
    for p in 0..count {
        let area = -closure[p];
        if area.norm() <= 1e-9 * scale[p].max(1e-300) {
            continue;
        }
        let (at, body) = match measured[p] {
            Some(w) => w,
            None => (pieces[p].centroid + area * (0.25 / area.norm() * volume_of(p).cbrt()), 0),
        };
        walls.push((p, at, area, None, body));
    }
    if !conditions.is_empty() {
        for (bi, table) in condition_of.iter().enumerate() {
            let Some(table) = table else { continue };
            let mine: Vec<usize> = (0..walls.len()).filter(|&w| walls[w].4 == bi).collect();
            let centres: Vec<[f64; 3]> = mine.iter().map(|&w| [walls[w].1.x, walls[w].1.y, walls[w].1.z]).collect();
            let nearest = bodies[bi].nearest_triangles(&centres)?;
            for (&w, t) in mine.iter().zip(nearest) {
                if let Some(Some(ci)) = table.get(t) {
                    walls[w].3 = Some(*ci);
                }
            }
        }
    }
    for (p, at, area, condition, body) in &walls {
        let patch = condition.map_or(body_patch(*body), |ci| solids + ci);
        faces.push(Face::boundary_face(*p, *at, *area, patch));
    }
    // Small pieces join their best-connected neighbour.
    let fractions: Vec<f64> = pieces.iter().map(|p| p.fraction).collect();
    let target = super::cut_cells::merge_small(&fractions, &faces);
    let mut new_index = vec![usize::MAX; count];
    let mut kept = Vec::new();
    for p in 0..count {
        if target[p] == Some(p) {
            new_index[p] = kept.len();
            kept.push(p);
        }
    }
    require(!kept.is_empty(), "obstacles remove every fluid cell")?;
    let mut volume = vec![0.0; kept.len()];
    let mut moment = vec![Vec3::ZERO; kept.len()];
    for p in 0..count {
        let Some(root) = target[p] else { continue };
        let v = volume_of(p);
        volume[new_index[root]] += v;
        moment[new_index[root]] += pieces[p].centroid * v;
    }
    let cells: Vec<PolyCell> = (0..kept.len()).map(|i| PolyCell { center: moment[i] * (1.0 / volume[i]), volume: volume[i] }).collect();
    let to_new = |p: usize| target[p].map(|root| new_index[root]);
    let mut merged_faces = Vec::with_capacity(faces.len());
    for f in faces {
        let Some(owner) = to_new(f.owner) else { continue };
        match f.neighbour {
            Some(n) => match to_new(n) {
                Some(other) if other == owner => {}
                Some(other) => merged_faces.push(Face::internal(owner, other, f.center, f.area)),
                None => {}
            },
            None => {
                let mut face = Face::boundary_face(owner, f.center, f.area, f.patch.unwrap_or(0));
                // A wall centre measured on the lattice can fall behind a merged cell's centre:
                // keep it on the face's outward side.
                if face.area.dot(face.center - cells[owner].center) <= 0.0 {
                    face.center = cells[owner].center + face.area * (0.25 / face.area.norm() * cells[owner].volume.cbrt());
                }
                merged_faces.push(face);
            }
        }
    }
    let mesh = PolyMesh::new(cells, merged_faces, patches)?;
    validate_topology(&mesh)?;
    for p in solids..mesh.patches().len() {
        require(
            mesh.faces().iter().any(|f| f.patch == Some(p)),
            format!("body face condition {} covers no fluid face of the grid; refine the grid near it", mesh.patches()[p].name),
        )?;
    }
    // Each kept piece is shown as its cell's box (two pieces of one cell share it).
    let export = VtuMesh { points: full_export.points, cells: kept.iter().map(|&p| full_export.cells[pieces[p].cell].clone()).collect() };
    export.validate()?;
    require(export.cells.len() == mesh.cell_count(), "export/solver cell count mismatch")?;
    Ok((mesh, export))
}

fn validate_topology(mesh: &PolyMesh) -> Result<()> {
    for (ci, cell) in mesh.cells().iter().enumerate() {
        let mut closure = Vec3::ZERO;
        let mut area_sum = 0.0;
        for &fi in &mesh.cell_faces()[ci] {
            let face = mesh.faces()[fi];
            let area = face.outward_area(ci);
            closure += area;
            area_sum += area.norm();
            require(
                area.dot(face.center - cell.center) > 0.0,
                "invalid obstacle cell face orientation",
            )?;
        }
        require(
            closure.norm() <= 1e-12 * area_sum,
            "obstacle cell is not closed",
        )?;
    }
    let mut reached = vec![false; mesh.cell_count()];
    let mut pending = vec![0];
    if let Some(first) = reached.first_mut() {
        *first = true;
    }
    while let Some(c) = pending.pop() {
        for &fi in &mesh.cell_faces()[c] {
            let f = mesh.faces()[fi];
            if let Some(n) = f.neighbour {
                let other = if f.owner == c { n } else { f.owner };
                if !reached[other] {
                    reached[other] = true;
                    pending.push(other);
                }
            }
        }
    }
    require(
        reached.iter().all(|r| *r),
        "obstacles create disconnected fluid regions",
    )
}

pub(super) fn write_report(
    case: &CollocatedFlowCase,
    mesh: &PolyMesh,
    flux: &[f64],
    out: &Path,
    manifest: &mut RunManifest,
) -> Result<()> {
    require(
        flux.len() == mesh.face_count(),
        "flow report flux size mismatch",
    )?;
    let mut patch_reports = Vec::new();
    let mut inflow = 0.0;
    let mut outflow = 0.0;
    for (pi, p) in mesh.patches().iter().enumerate() {
        let mut area = 0.0;
        let mut volume_flux = 0.0;
        let mut faces = 0;
        for (fi, f) in mesh
            .faces()
            .iter()
            .enumerate()
            .filter(|(_, f)| f.patch == Some(pi))
        {
            area += f.area.norm();
            volume_flux += flux[fi];
            faces += 1;
            if flux[fi] < 0.0 {
                inflow -= flux[fi];
            } else {
                outflow += flux[fi];
            }
        }
        patch_reports.push(
            serde_json::json!({"name":p.name,"face_count":faces,"area_m2":area,
            "outward_volume_flux_m3_per_s":volume_flux,
            "outward_mass_flux_kg_per_s":volume_flux*case.config.density}),
        );
    }
    let scale = inflow.max(outflow);
    let relative = if scale > 0.0 {
        (outflow - inflow).abs() / scale
    } else {
        0.0
    };
    let background: usize = case.mesh.dimensions.iter().product();
    let mut bodies = vec![];
    // The staircase statistics describe whole-cell bodies; cut cells report their own above.
    if !case.bodies.is_empty() && !case.cut_cells {
        let layout = layout(&case.mesh, &case.obstacles, &case.bodies)?;
        let [nx, ny, _] = case.mesh.dimensions;
        let (xs, ys, zs) = (case.mesh.nodes(0), case.mesh.nodes(1), case.mesh.nodes(2));
        let width = |v: &[f64], i: usize| v.get(i + 1).zip(v.get(i)).map_or(0.0, |(b, a)| b - a);
        let cell_volume = |c: usize| width(&xs, c % nx) * width(&ys, (c / nx) % ny) * width(&zs, c / (nx * ny));
        for (bi, body) in case.bodies.iter().enumerate() {
            let index = case.obstacles.len() + bi;
            let owned: Vec<usize> = layout
                .cell_solid
                .iter()
                .enumerate()
                .filter(|(_, s)| **s == Some(index))
                .map(|(c, _)| c)
                .collect();
            let cells = owned.len();
            let surface = body.validate()?;
            let represented: f64 = owned.iter().map(|c| cell_volume(*c)).sum();
            bodies.push(serde_json::json!({"name": body.name, "triangles": body.triangles.len(), "cells": cells,
                "surface_volume_m3": surface, "represented_volume_m3": represented,
                "relative_volume_error": (represented - surface) / surface,
                "representation": "staircase_cell_centre_winding"}));
        }
    }
    // With cut cells: how much fluid the cells hold against the exact fluid volume.
    let cut = if case.cut_cells {
        let domain = case.mesh.lengths.x * case.mesh.lengths.y * case.mesh.lengths.z;
        let mut solid = 0.0;
        for body in &case.bodies {
            solid += body.validate()?;
        }
        let (xs, ys, zs) = (case.mesh.nodes(0), case.mesh.nodes(1), case.mesh.nodes(2));
        for o in &case.obstacles {
            solid += (xs[o.max[0]] - xs[o.min[0]]) * (ys[o.max[1]] - ys[o.min[1]]) * (zs[o.max[2]] - zs[o.min[2]]);
        }
        let represented: f64 = mesh.cells().iter().map(|c| c.volume).sum();
        let exact = domain - solid;
        Some(serde_json::json!({"representation":"cut_cells_sampled_sub_lattice",
            "max_sub_cells_per_axis": super::cut_cells::MAX_SUB, "merge_below_fraction": super::cut_cells::MERGE_BELOW,
            "represented_fluid_volume_m3": represented, "exact_fluid_volume_m3": exact,
            "relative_fluid_volume_error": (represented - exact) / exact}))
    } else {
        None
    };
    let report = serde_json::json!({
        "schema":"cfd-solvers-rs/grid-obstacle-mesh/1", "geometry_kind":"discretized_grid_primitives",
        "dimensions":case.mesh.dimensions,"lengths_m":case.mesh.lengths,"obstacles":case.obstacles,
        "fluid_cells":mesh.cell_count(),"excluded_solid_cells":background-mesh.cell_count(),
        "faces":mesh.face_count(),"connected_fluid_regions":1,"cell_closure_checked":true,
        "solver_and_export_cell_order":"same_retained_background_x_fastest_order",
        "patches":patch_reports,"inflow_m3_per_s":inflow,"outflow_m3_per_s":outflow,
        "relative_boundary_flux_imbalance":relative,"physics":"isothermal_incompressible_laminar",
        "physical_validation":"not_established"
    });
    let mut report = report;
    if let (Some(cut), Some(object)) = (cut, report.as_object_mut()) {
        object.insert("cut_cells".into(), cut);
    }
    if !bodies.is_empty() {
        if let Some(object) = report.as_object_mut() {
            object.insert("surface_bodies".into(), serde_json::Value::Array(bodies));
        }
    }
    std::fs::write(
        out.join("mesh_report.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    manifest.files.push("mesh_report.json".into());
    manifest
        .metrics
        .insert("fluid_cell_count".into(), mesh.cell_count() as f64);
    manifest.metrics.insert(
        "excluded_solid_cell_count".into(),
        (background - mesh.cell_count()) as f64,
    );
    manifest
        .metrics
        .insert("boundary_inflow_m3_per_s".into(), inflow);
    manifest
        .metrics
        .insert("boundary_outflow_m3_per_s".into(), outflow);
    manifest
        .metrics
        .insert("relative_boundary_flux_imbalance".into(), relative);
    Ok(())
}

/// For each face of the retained fluid mesh, the background-grid face it came from and
/// whether its orientation was flipped. Mirrors the face filter in [`build`] exactly.
pub(super) fn fluid_face_origins(
    spec: &PolyMeshSpec,
    obstacles: &[GridObstacle],
    bodies: &[super::bodies::SurfaceBody],
) -> Result<Vec<(usize, bool)>> {
    let (base, _) = spec.build()?;
    if obstacles.is_empty() && bodies.is_empty() {
        return Ok((0..base.face_count()).map(|f| (f, false)).collect());
    }
    let solid: Vec<bool> = layout(spec, obstacles, bodies)?
        .cell_solid
        .iter()
        .map(Option::is_some)
        .collect();
    let fluid = |c: usize| solid.get(c) == Some(&false);
    let mut origins = Vec::new();
    for (fi, f) in base.faces().iter().enumerate() {
        match (fluid(f.owner), f.neighbour.map(fluid)) {
            (true, _) => origins.push((fi, false)),
            (false, Some(true)) => origins.push((fi, true)),
            _ => {}
        }
    }
    Ok(origins)
}
