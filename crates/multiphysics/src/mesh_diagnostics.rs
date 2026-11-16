//! Bounded, index-based diagnostics for an explicit triangle surface.
//!
//! This is an advisory inspection of supplied mesh connectivity, not a CAD
//! geometry operation, volume mesher, mesh repairer or solver-readiness test.
//! Vertices are compared by index only. Coincident coordinates with different
//! indices, collinear triangles and self-intersections are not detected.

use std::collections::BTreeMap;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{Error, Result};

/// Limits apply before building the edge incidence table.
pub const MAX_SURFACE_VERTICES: usize = 100_000;
pub const MAX_SURFACE_TRIANGLES: usize = 100_000;

/// Explicit indexed surface in the caller's declared coordinate system and units.
/// The inspection does not infer units or merge vertices by position.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SurfaceMesh {
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[u32; 3]>,
}

/// Connectivity evidence only. A closed edge pairing cannot establish that
/// the surface bounds a valid physical domain or that it can be meshed.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct SurfaceReport {
    pub vertices: usize,
    pub triangles: usize,
    pub unused_vertices: usize,
    pub repeated_index_triangles: usize,
    pub indexed_edges: usize,
    pub boundary_edges: usize,
    pub non_manifold_edges: usize,
    pub inconsistent_winding_edges: usize,
    /// All nondegenerate indexed edges have two oppositely directed incidences,
    /// and no triangle has a repeated index. This is not watertightness proof.
    pub indexed_edge_pairing_complete: bool,
}

/// Inspect one bounded triangle surface. Invalid input never returns partial evidence.
pub fn inspect(mesh: &SurfaceMesh) -> Result<SurfaceReport> {
    if mesh.vertices.is_empty() || mesh.triangles.is_empty() {
        return Err(Error::Invalid("surface needs vertices and triangles".into()));
    }
    if mesh.vertices.len() > MAX_SURFACE_VERTICES || mesh.triangles.len() > MAX_SURFACE_TRIANGLES {
        return Err(Error::Limit("surface exceeds indexed mesh inspection limits".into()));
    }
    if mesh.vertices.iter().flatten().any(|value| !value.is_finite()) {
        return Err(Error::Invalid("surface contains a non-finite coordinate".into()));
    }
    // Each edge records incident faces and a direction sum relative to its
    // sorted index pair. Counts are bounded by MAX_SURFACE_TRIANGLES.
    let mut edges: BTreeMap<(u32, u32), (usize, i64)> = BTreeMap::new();
    let mut referenced = vec![false; mesh.vertices.len()];
    let mut repeated_index_triangles = 0;
    for triangle in &mesh.triangles {
        let [a, b, c] = *triangle;
        for vertex in [a, b, c] {
            let Some(used) = referenced.get_mut(vertex as usize) else {
                return Err(Error::Invalid("triangle references a missing vertex".into()));
            };
            *used = true;
        }
        if a == b || b == c || c == a {
            repeated_index_triangles += 1;
            continue;
        }
        for (from, to) in [(a, b), (b, c), (c, a)] {
            let (key, direction) = if from < to { ((from, to), 1) } else { ((to, from), -1) };
            let entry = edges.entry(key).or_insert((0, 0));
            entry.0 += 1;
            entry.1 += direction;
        }
    }
    let boundary_edges = edges.values().filter(|(count, _)| *count == 1).count();
    let non_manifold_edges = edges.values().filter(|(count, _)| *count > 2).count();
    let inconsistent_winding_edges = edges.values().filter(|(count, direction)| *count == 2 && *direction != 0).count();
    let indexed_edge_pairing_complete =
        !edges.is_empty() && repeated_index_triangles == 0 && boundary_edges == 0 && non_manifold_edges == 0 && inconsistent_winding_edges == 0;
    Ok(SurfaceReport {
        vertices: mesh.vertices.len(),
        triangles: mesh.triangles.len(),
        unused_vertices: referenced.iter().filter(|used| !**used).count(),
        repeated_index_triangles,
        indexed_edges: edges.len(),
        boundary_edges,
        non_manifold_edges,
        inconsistent_winding_edges,
        indexed_edge_pairing_complete,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tetrahedron() -> SurfaceMesh {
        SurfaceMesh {
            vertices: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            triangles: vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
        }
    }

    #[test]
    fn closed_indexed_tetrahedron_and_missing_face_are_distinguished() {
        let mut mesh = tetrahedron();
        let report = inspect(&mesh).unwrap();
        assert_eq!(report.indexed_edges, 6);
        assert_eq!(report.boundary_edges, 0);
        assert!(report.indexed_edge_pairing_complete);
        mesh.triangles.pop();
        let report = inspect(&mesh).unwrap();
        assert_eq!(report.boundary_edges, 3);
        assert!(!report.indexed_edge_pairing_complete);
    }

    #[test]
    fn repeated_indices_nonmanifold_and_wrong_winding_are_reported() {
        let mut mesh = tetrahedron();
        mesh.triangles[0] = [0, 1, 2];
        assert_eq!(inspect(&mesh).unwrap().inconsistent_winding_edges, 3);
        mesh.triangles.push([0, 1, 2]);
        assert_eq!(inspect(&mesh).unwrap().non_manifold_edges, 3);
        mesh.triangles.push([0, 0, 1]);
        assert_eq!(inspect(&mesh).unwrap().repeated_index_triangles, 1);
        assert!(!inspect(&mesh).unwrap().indexed_edge_pairing_complete);
    }

    #[test]
    fn invalid_and_oversized_inputs_are_rejected() {
        let mut mesh = tetrahedron();
        mesh.vertices[0][0] = f64::NAN;
        assert!(inspect(&mesh).is_err());
        mesh = tetrahedron();
        mesh.triangles[0][0] = u32::MAX;
        assert!(inspect(&mesh).is_err());
        mesh = tetrahedron();
        mesh.vertices.resize(MAX_SURFACE_VERTICES + 1, [0.0; 3]);
        assert!(matches!(inspect(&mesh), Err(Error::Limit(_))));
        assert!(inspect(&SurfaceMesh { vertices: vec![], triangles: vec![] }).is_err());
    }

    #[test]
    fn coincident_vertices_are_not_silently_welded() {
        let mut mesh = tetrahedron();
        mesh.vertices.push(mesh.vertices[0]);
        let report = inspect(&mesh).unwrap();
        assert_eq!(report.unused_vertices, 1);
        assert_eq!(report.vertices, 5);
    }
}
