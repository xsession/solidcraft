//! Residual/jump a-posteriori estimator for P1 Poisson, Dorfler marking and conforming
//! longest-edge refinement of Tri3 meshes, ported from upstream `cfd/fem/adaptivity` (MIT,
//! cfd_solvers contributors).
use super::{BoundaryEdge2, Mesh2D, Node2, Tri3};
use crate::{require, Result, SolverError};
use std::collections::{BTreeMap, BTreeSet};

type Edge = (usize, usize);

fn edge_key(a: usize, b: usize) -> Edge {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

fn signed_area(a: Node2, b: Node2, c: Node2) -> f64 {
    0.5 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x))
}

fn edge_length(a: Node2, b: Node2) -> f64 {
    (b.x - a.x).hypot(b.y - a.y)
}

fn triangle_gradient(mesh: &Mesh2D, tri: &Tri3, u: &[f64]) -> Result<[f64; 2]> {
    let (a, b, c) = (
        mesh.nodes[tri.node[0]],
        mesh.nodes[tri.node[1]],
        mesh.nodes[tri.node[2]],
    );
    let det = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    if det.abs() <= 1e-30 {
        return Err(SolverError::NumericalFailure(
            "degenerate triangle in error estimator".into(),
        ));
    }
    let inv = 1.0 / det;
    let grad = [
        [(b.y - c.y) * inv, (c.x - b.x) * inv],
        [(c.y - a.y) * inv, (a.x - c.x) * inv],
        [(a.y - b.y) * inv, (b.x - a.x) * inv],
    ];
    let mut g = [0.0; 2];
    for i in 0..3 {
        g[0] += u[tri.node[i]] * grad[i][0];
        g[1] += u[tri.node[i]] * grad[i][1];
    }
    Ok(g)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PoissonErrorEstimate2D {
    pub element_indicator: Vec<f64>,
    pub global_indicator: f64,
}

/// Estimator for `-div(k grad u) = f` on P1 triangles: `h_K^2 |K| f(centroid)^2` plus half of
/// each interior edge's `|e| [k grad u . n]^2`. An indicator, not an exact error norm.
pub fn estimate_poisson_error_tri3(
    mesh: &Mesh2D,
    nodal_solution: &[f64],
    source: &dyn Fn(Node2) -> f64,
    diffusivity: f64,
) -> Result<PoissonErrorEstimate2D> {
    mesh.validate()?;
    require(
        nodal_solution.len() == mesh.node_count() && diffusivity > 0.0 && diffusivity.is_finite(),
        "invalid Poisson estimator input",
    )?;
    let ne = mesh.triangles.len();
    let mut eta2 = vec![0.0; ne];
    let mut grad = Vec::with_capacity(ne);
    let mut adjacency: BTreeMap<Edge, Vec<usize>> = BTreeMap::new();
    for (e, tri) in mesh.triangles.iter().enumerate() {
        grad.push(triangle_gradient(mesh, tri, nodal_solution)?);
        let (a, b, c) = (
            mesh.nodes[tri.node[0]],
            mesh.nodes[tri.node[1]],
            mesh.nodes[tri.node[2]],
        );
        let area = signed_area(a, b, c).abs();
        let h = edge_length(a, b)
            .max(edge_length(b, c))
            .max(edge_length(c, a));
        let f = source(Node2 {
            x: (a.x + b.x + c.x) / 3.0,
            y: (a.y + b.y + c.y) / 3.0,
        });
        eta2[e] += h * h * area * f * f;
        for (p, q) in [(0, 1), (1, 2), (2, 0)] {
            adjacency
                .entry(edge_key(tri.node[p], tri.node[q]))
                .or_default()
                .push(e);
        }
    }
    for (&(i0, i1), adj) in &adjacency {
        if adj.len() != 2 {
            continue;
        }
        let (p0, p1) = (mesh.nodes[i0], mesh.nodes[i1]);
        let (dx, dy) = (p1.x - p0.x, p1.y - p0.y);
        let length = dx.hypot(dy);
        if length <= 0.0 || !length.is_finite() {
            continue;
        }
        let (nx, ny) = (dy / length, -dx / length);
        let (t0, t1) = (grad[adj[0]], grad[adj[1]]);
        let jump = diffusivity * ((t0[0] - t1[0]) * nx + (t0[1] - t1[1]) * ny);
        let contribution = 0.5 * length * jump * jump;
        eta2[adj[0]] += contribution;
        eta2[adj[1]] += contribution;
    }
    require(
        eta2.iter().all(|v| v.is_finite()),
        "non-finite error indicator (check the source term)",
    )?;
    let element_indicator: Vec<f64> = eta2.iter().map(|v| v.max(0.0).sqrt()).collect();
    let total: f64 = eta2.iter().sum();
    Ok(PoissonErrorEstimate2D {
        element_indicator,
        global_indicator: total.sqrt(),
    })
}

/// Dorfler bulk marking: largest indicators first until `target_fraction` of the squared total
/// is reached (`target_fraction` in (0, 1]). An all-zero indicator marks the first element.
pub fn mark_dorfler(indicator: &[f64], target_fraction: f64) -> Result<Vec<bool>> {
    require(
        !indicator.is_empty() && target_fraction > 0.0 && target_fraction <= 1.0,
        "invalid Dorfler marking input",
    )?;
    require(
        indicator.iter().all(|v| *v >= 0.0 && v.is_finite()),
        "invalid error indicator",
    )?;
    let mut order: Vec<usize> = (0..indicator.len()).collect();
    order.sort_by(|a, b| indicator[*b].total_cmp(&indicator[*a]));
    let total: f64 = indicator.iter().map(|v| v * v).sum();
    let mut marked = vec![false; indicator.len()];
    if total == 0.0 {
        marked[order[0]] = true;
        return Ok(marked);
    }
    let target = target_fraction * total;
    let mut accumulated = 0.0;
    for e in order {
        marked[e] = true;
        accumulated += indicator[e] * indicator[e];
        if accumulated >= target {
            break;
        }
    }
    Ok(marked)
}

fn add_oriented_triangle(out: &mut Mesh2D, a: usize, mut b: usize, mut c: usize) -> Result<()> {
    let area = signed_area(out.nodes[a], out.nodes[b], out.nodes[c]);
    if area.abs() <= 1e-30 {
        return Err(SolverError::NumericalFailure(
            "adaptive refinement produced degenerate triangle".into(),
        ));
    }
    if area < 0.0 {
        std::mem::swap(&mut b, &mut c);
    }
    out.triangles.push(Tri3 { node: [a, b, c] });
    Ok(())
}

/// Conforming longest-edge refinement: each marked triangle selects its longest edge (first on
/// ties), every triangle touching a selected edge is split with a 1/2/3-edge template, and
/// boundary patches are kept on split boundary edges.
pub fn refine_tri3_longest_edges(mesh: &Mesh2D, marked_element: &[bool]) -> Result<Mesh2D> {
    mesh.validate()?;
    require(
        marked_element.len() == mesh.element_count(),
        "adaptive mark size mismatch",
    )?;
    let mut split_edges = BTreeSet::new();
    for (tri, _) in mesh
        .triangles
        .iter()
        .zip(marked_element)
        .filter(|(_, m)| **m)
    {
        let edges = [
            edge_key(tri.node[0], tri.node[1]),
            edge_key(tri.node[1], tri.node[2]),
            edge_key(tri.node[2], tri.node[0]),
        ];
        let mut best = 0;
        let mut best_len = f64::NEG_INFINITY;
        for (k, e) in edges.iter().enumerate() {
            let l = edge_length(mesh.nodes[e.0], mesh.nodes[e.1]);
            // `std::max_element` keeps the first maximum.
            if l > best_len {
                best = k;
                best_len = l;
            }
        }
        split_edges.insert(edges[best]);
    }
    if split_edges.is_empty() {
        return Ok(mesh.clone());
    }
    let mut out = Mesh2D {
        nodes: mesh.nodes.clone(),
        triangles: Vec::new(),
        boundary_edges: Vec::new(),
        boundary_node: Vec::new(),
    };
    let mut midpoint = BTreeMap::new();
    for &(i, j) in &split_edges {
        let (a, b) = (mesh.nodes[i], mesh.nodes[j]);
        midpoint.insert((i, j), out.nodes.len());
        out.nodes.push(Node2 {
            x: 0.5 * (a.x + b.x),
            y: 0.5 * (a.y + b.y),
        });
    }
    for tri in &mesh.triangles {
        let [a, b, c] = tri.node;
        let m0 = midpoint.get(&edge_key(a, b)).copied();
        let m1 = midpoint.get(&edge_key(b, c)).copied();
        let m2 = midpoint.get(&edge_key(c, a)).copied();
        let o = &mut out;
        match (m0, m1, m2) {
            (None, None, None) => add_oriented_triangle(o, a, b, c)?,
            (Some(m0), None, None) => {
                add_oriented_triangle(o, a, m0, c)?;
                add_oriented_triangle(o, m0, b, c)?;
            }
            (None, Some(m1), None) => {
                add_oriented_triangle(o, b, m1, a)?;
                add_oriented_triangle(o, m1, c, a)?;
            }
            (None, None, Some(m2)) => {
                add_oriented_triangle(o, c, m2, b)?;
                add_oriented_triangle(o, m2, a, b)?;
            }
            (Some(m0), Some(m1), None) => {
                add_oriented_triangle(o, b, m1, m0)?;
                add_oriented_triangle(o, m0, m1, c)?;
                add_oriented_triangle(o, a, m0, c)?;
            }
            (None, Some(m1), Some(m2)) => {
                add_oriented_triangle(o, c, m2, m1)?;
                add_oriented_triangle(o, m1, m2, a)?;
                add_oriented_triangle(o, b, m1, a)?;
            }
            (Some(m0), None, Some(m2)) => {
                add_oriented_triangle(o, a, m0, m2)?;
                add_oriented_triangle(o, m2, m0, b)?;
                add_oriented_triangle(o, c, m2, b)?;
            }
            (Some(m0), Some(m1), Some(m2)) => {
                add_oriented_triangle(o, a, m0, m2)?;
                add_oriented_triangle(o, m0, b, m1)?;
                add_oriented_triangle(o, m2, m1, c)?;
                add_oriented_triangle(o, m0, m1, m2)?;
            }
        }
    }
    for edge in &mesh.boundary_edges {
        match midpoint.get(&edge_key(edge.node[0], edge.node[1])) {
            None => out.boundary_edges.push(*edge),
            Some(&m) => {
                out.boundary_edges.push(BoundaryEdge2 {
                    node: [edge.node[0], m],
                    patch: edge.patch,
                });
                out.boundary_edges.push(BoundaryEdge2 {
                    node: [m, edge.node[1]],
                    patch: edge.patch,
                });
            }
        }
    }
    out.boundary_node = vec![false; out.nodes.len()];
    for e in &out.boundary_edges {
        out.boundary_node[e.node[0]] = true;
        out.boundary_node[e.node[1]] = true;
    }
    out.validate()?;
    Ok(out)
}
