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
    spec: PolyMeshSpec,
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
    spec: PolyMeshSpec,
    obstacles: &[GridObstacle],
    bodies: &[super::bodies::SurfaceBody],
) -> Result<(PolyMesh, VtuMesh)> {
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
    let mesh = PolyMesh::new(cells, faces, patches)?;
    validate_topology(&mesh)?;
    for p in EXTERIOR.len()..mesh.patches().len() {
        require(
            mesh.faces().iter().any(|f| f.patch == Some(p)),
            "an obstacle has no exposed fluid boundary",
        )?;
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
    if !case.bodies.is_empty() {
        let layout = layout(case.mesh, &case.obstacles, &case.bodies)?;
        let [nx, ny, nz] = case.mesh.dimensions;
        let cell_volume =
            case.mesh.lengths.x * case.mesh.lengths.y * case.mesh.lengths.z / (nx * ny * nz) as f64;
        for (bi, body) in case.bodies.iter().enumerate() {
            let index = case.obstacles.len() + bi;
            let cells = layout
                .cell_solid
                .iter()
                .filter(|s| **s == Some(index))
                .count();
            let surface = body.validate()?;
            let represented = cells as f64 * cell_volume;
            bodies.push(serde_json::json!({"name": body.name, "triangles": body.triangles.len(), "cells": cells,
                "surface_volume_m3": surface, "represented_volume_m3": represented,
                "relative_volume_error": (represented - surface) / surface,
                "representation": "staircase_cell_centre_winding"}));
        }
    }
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
    spec: PolyMeshSpec,
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
