//! Lowest-order Nedelec edge elements: driven 2-D in-plane Maxwell on Tri3 (dense), driven 3-D
//! curl-curl Maxwell on Tet4 (sparse GMRES) with edge current ports, Thevenin extraction and
//! resonator quality, plus face-jump error indicators, Dorfler marking and conforming
//! longest-edge Tet4 refinement, ported from upstream `cfd/em/{edge_fem2d, edge_fem3d,
//! adaptivity3d}` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use crate::{
    core::{
        solve_complex_sparse, ComplexCsrBuilder, ComplexSparseSolveConfig, IterativeSolverResult,
    },
    fem::{Mesh2D, Mesh3D, Node2, Point3, Tet4, Tri3},
    require, Complex64 as C, Result, SolverError,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    f64::consts::PI,
};

const EPS0: f64 = 8.8541878128e-12;
const MU0: f64 = 1.25663706212e-6;
const ZERO: C = C::new(0.0, 0.0);

pub type ComplexVec2 = [C; 2];
pub type ComplexVec3 = [C; 3];
pub type RealVec3 = [f64; 3];

const TRI_EDGES: [[usize; 2]; 3] = [[0, 1], [1, 2], [2, 0]];
const TET_EDGES: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

fn canonical_edge(a: usize, b: usize) -> ([usize; 2], f64) {
    if a < b {
        ([a, b], 1.0)
    } else {
        ([b, a], -1.0)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NedelecTri3Geometry {
    pub area: f64,
    /// Physical gradients of the barycentric coordinates.
    pub gradient_lambda: [[f64; 2]; 3],
}

/// Local edges are oriented (0->1), (1->2), (2->0).
pub fn nedelec_tri3_geometry(mesh: &Mesh2D, triangle: &Tri3) -> Result<NedelecTri3Geometry> {
    require(
        triangle.node.iter().all(|n| *n < mesh.nodes.len()),
        "Nedelec triangle node out of range",
    )?;
    let (a, b, c) = (
        mesh.nodes[triangle.node[0]],
        mesh.nodes[triangle.node[1]],
        mesh.nodes[triangle.node[2]],
    );
    let t = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    require(t.abs() > 1e-30, "degenerate Nedelec triangle")?;
    Ok(NedelecTri3Geometry {
        area: 0.5 * t.abs(),
        gradient_lambda: [
            [(b.y - c.y) / t, (c.x - b.x) / t],
            [(c.y - a.y) / t, (a.x - c.x) / t],
            [(a.y - b.y) / t, (b.x - a.x) / t],
        ],
    })
}

pub fn nedelec_tri3_basis(
    g: &NedelecTri3Geometry,
    local_edge: usize,
    barycentric: [f64; 3],
) -> Result<[f64; 2]> {
    require(local_edge < 3, "Nedelec local edge index")?;
    let [i, j] = TRI_EDGES[local_edge];
    let l = barycentric;
    Ok([
        l[i] * g.gradient_lambda[j][0] - l[j] * g.gradient_lambda[i][0],
        l[i] * g.gradient_lambda[j][1] - l[j] * g.gradient_lambda[i][1],
    ])
}

pub fn nedelec_tri3_curl(g: &NedelecTri3Geometry, local_edge: usize) -> Result<f64> {
    require(local_edge < 3, "Nedelec local edge index")?;
    let [i, j] = TRI_EDGES[local_edge];
    let (gi, gj) = (g.gradient_lambda[i], g.gradient_lambda[j]);
    Ok(2.0 * (gi[0] * gj[1] - gi[1] * gj[0]))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeMaxwell2DConfig {
    pub frequency_hz: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub conductivity_s_per_m: f64,
}

impl Default for EdgeMaxwell2DConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 1e8,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            conductivity_s_per_m: 0.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct EdgeMaxwell2DResult {
    /// Lower node index -> higher node index.
    pub edges: Vec<[usize; 2]>,
    /// Line integral of E along each edge, V.
    pub edge_voltage_v: Vec<C>,
    /// In-plane E at triangle centroids, V/m.
    pub electric_centroid_v_per_m: Vec<ComplexVec2>,
}

fn dense_lu_solve(mut m: Vec<C>, mut rhs: Vec<C>, n: usize) -> Result<Vec<C>> {
    for k in 0..n {
        let mut pivot = k;
        let mut best = m[k * n + k].norm();
        for row in k + 1..n {
            let v = m[row * n + k].norm();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best.is_nan() || best <= 100.0 * f64::EPSILON {
            return Err(SolverError::Singular("singular edge-Maxwell system".into()));
        }
        if pivot != k {
            for col in k..n {
                m.swap(k * n + col, pivot * n + col);
            }
            rhs.swap(k, pivot);
        }
        let diag = m[k * n + k];
        for row in k + 1..n {
            let factor = m[row * n + k] / diag;
            if factor.norm() == 0.0 {
                continue;
            }
            m[row * n + k] = ZERO;
            for col in k + 1..n {
                let v = m[k * n + col];
                m[row * n + col] -= factor * v;
            }
            let r = rhs[k];
            rhs[row] -= factor * r;
        }
    }
    let mut x = vec![ZERO; n];
    for row in (0..n).rev() {
        let mut v = rhs[row];
        for col in row + 1..n {
            v -= m[row * n + col] * x[col];
        }
        x[row] = v / m[row * n + row];
    }
    Ok(x)
}

/// Driven in-plane full-wave Maxwell with PEC boundary edges (`exp(+j omega t)`), dense LU.
pub fn solve_driven_edge_maxwell_2d(
    mesh: &Mesh2D,
    config: &EdgeMaxwell2DConfig,
    impressed_current_density_a_per_m2: &dyn Fn(Node2) -> ComplexVec2,
) -> Result<EdgeMaxwell2DResult> {
    mesh.validate()?;
    require(
        config.frequency_hz > 0.0
            && config.relative_permittivity > 0.0
            && config.relative_permeability > 0.0
            && config.conductivity_s_per_m >= 0.0,
        "invalid edge-Maxwell configuration/source",
    )?;
    let mut edge_index: BTreeMap<[usize; 2], usize> = BTreeMap::new();
    let mut edges = Vec::new();
    let mut tri_edges = vec![[0usize; 3]; mesh.triangles.len()];
    let mut tri_signs = vec![[0.0f64; 3]; mesh.triangles.len()];
    for (e, tri) in mesh.triangles.iter().enumerate() {
        for local in 0..3 {
            let (edge, sign) =
                canonical_edge(tri.node[TRI_EDGES[local][0]], tri.node[TRI_EDGES[local][1]]);
            let index = *edge_index.entry(edge).or_insert_with(|| {
                edges.push(edge);
                edges.len() - 1
            });
            tri_edges[e][local] = index;
            tri_signs[e][local] = sign;
        }
    }
    let n = edges.len();
    require(
        n <= 4000,
        "2-D edge-Maxwell dense solve is limited to 4000 edges",
    )?;
    let mut m = vec![ZERO; n * n];
    let mut rhs = vec![ZERO; n];
    let omega = 2.0 * PI * config.frequency_hz;
    let (mu, eps) = (
        MU0 * config.relative_permeability,
        EPS0 * config.relative_permittivity,
    );
    let quadrature = [
        [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
    ];
    for (e, tri) in mesh.triangles.iter().enumerate() {
        let g = nedelec_tri3_geometry(mesh, tri)?;
        let mut mass = [[0.0; 3]; 3];
        let mut load = [ZERO; 3];
        for lambda in quadrature {
            let mut basis = [[0.0; 2]; 3];
            for a in 0..3 {
                basis[a] = nedelec_tri3_basis(&g, a, lambda)?;
            }
            let mut p = Node2 { x: 0.0, y: 0.0 };
            for i in 0..3 {
                p.x += lambda[i] * mesh.nodes[tri.node[i]].x;
                p.y += lambda[i] * mesh.nodes[tri.node[i]].y;
            }
            let current = impressed_current_density_a_per_m2(p);
            let w = g.area / 3.0;
            for a in 0..3 {
                load[a] += (current[0] * basis[a][0] + current[1] * basis[a][1]) * w;
                for b in 0..3 {
                    mass[a][b] += w * (basis[a][0] * basis[b][0] + basis[a][1] * basis[b][1]);
                }
            }
        }
        let mut curl = [0.0; 3];
        for a in 0..3 {
            curl[a] = nedelec_tri3_curl(&g, a)?;
        }
        for a in 0..3 {
            let (ga, sa) = (tri_edges[e][a], tri_signs[e][a]);
            rhs[ga] += C::new(0.0, -omega) * load[a] * sa;
            for b in 0..3 {
                let (gb, sb) = (tri_edges[e][b], tri_signs[e][b]);
                let local = C::new(g.area * curl[a] * curl[b] / mu, 0.0)
                    + C::new(-omega * omega * eps, omega * config.conductivity_s_per_m)
                        * mass[a][b];
                m[ga * n + gb] += local * (sa * sb);
            }
        }
    }
    let boundary: BTreeSet<[usize; 2]> = mesh
        .boundary_edges
        .iter()
        .map(|b| canonical_edge(b.node[0], b.node[1]).0)
        .collect();
    for edge in &boundary {
        let fixed = *edge_index.get(edge).ok_or_else(|| {
            SolverError::InvalidInput("boundary edge is not a triangle edge".into())
        })?;
        for col in 0..n {
            m[fixed * n + col] = ZERO;
            m[col * n + fixed] = ZERO;
        }
        m[fixed * n + fixed] = C::new(1.0, 0.0);
        rhs[fixed] = ZERO;
    }
    let voltage = dense_lu_solve(m, rhs, n)?;
    let mut centroid_e = vec![[ZERO; 2]; mesh.triangles.len()];
    for (e, tri) in mesh.triangles.iter().enumerate() {
        let g = nedelec_tri3_geometry(mesh, tri)?;
        let mut f = [ZERO; 2];
        for local in 0..3 {
            let basis = nedelec_tri3_basis(&g, local, [1.0 / 3.0; 3])?;
            let dof = voltage[tri_edges[e][local]] * tri_signs[e][local];
            f[0] += dof * basis[0];
            f[1] += dof * basis[1];
        }
        centroid_e[e] = f;
    }
    Ok(EdgeMaxwell2DResult {
        edges,
        edge_voltage_v: voltage,
        electric_centroid_v_per_m: centroid_e,
    })
}

fn sub3(a: Point3, b: Point3) -> RealVec3 {
    [a.x - b.x, a.y - b.y, a.z - b.z]
}
fn cross3(a: RealVec3, b: RealVec3) -> RealVec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot3(a: RealVec3, b: RealVec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm3(a: RealVec3) -> f64 {
    dot3(a, a).sqrt()
}
fn cnorm2(a: ComplexVec3) -> f64 {
    a[0].norm_sqr() + a[1].norm_sqr() + a[2].norm_sqr()
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NedelecTet4Geometry {
    pub volume_m3: f64,
    pub gradient_lambda: [RealVec3; 4],
}

/// Local edges (0->1), (0->2), (0->3), (1->2), (1->3), (2->3).
pub fn nedelec_tet4_geometry(mesh: &Mesh3D, tet: &Tet4) -> Result<NedelecTet4Geometry> {
    require(
        tet.node.iter().all(|n| *n < mesh.nodes.len()),
        "Nedelec tetrahedron node out of range",
    )?;
    let p = tet.node.map(|i| mesh.nodes[i]);
    let (v1, v2, v3) = (sub3(p[1], p[0]), sub3(p[2], p[0]), sub3(p[3], p[0]));
    let det = dot3(v1, cross3(v2, v3));
    require(det.abs() > 1e-30, "degenerate Nedelec tetrahedron")?;
    let (g1, g2, g3) = (cross3(v2, v3), cross3(v3, v1), cross3(v1, v2));
    let mut g = NedelecTet4Geometry {
        volume_m3: det.abs() / 6.0,
        ..Default::default()
    };
    g.gradient_lambda[1] = g1.map(|v| v / det);
    g.gradient_lambda[2] = g2.map(|v| v / det);
    g.gradient_lambda[3] = g3.map(|v| v / det);
    for d in 0..3 {
        g.gradient_lambda[0][d] =
            -(g.gradient_lambda[1][d] + g.gradient_lambda[2][d] + g.gradient_lambda[3][d]);
    }
    Ok(g)
}

pub fn nedelec_tet4_basis(
    g: &NedelecTet4Geometry,
    local_edge: usize,
    barycentric: [f64; 4],
) -> Result<RealVec3> {
    require(local_edge < 6, "Nedelec Tet4 local edge index")?;
    let [i, j] = TET_EDGES[local_edge];
    Ok(std::array::from_fn(|d| {
        barycentric[i] * g.gradient_lambda[j][d] - barycentric[j] * g.gradient_lambda[i][d]
    }))
}

pub fn nedelec_tet4_curl(g: &NedelecTet4Geometry, local_edge: usize) -> Result<RealVec3> {
    require(local_edge < 6, "Nedelec Tet4 local edge index")?;
    let [i, j] = TET_EDGES[local_edge];
    Ok(cross3(g.gradient_lambda[i], g.gradient_lambda[j]).map(|v| v * 2.0))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeMaxwell3DConfig {
    pub frequency_hz: f64,
    pub relative_permittivity: f64,
    /// Optional diagonal anisotropy; zero entries fall back to `relative_permittivity`.
    pub relative_permittivity_x: f64,
    pub relative_permittivity_y: f64,
    pub relative_permittivity_z: f64,
    pub relative_permeability: f64,
    pub conductivity_s_per_m: f64,
    pub linear: ComplexSparseSolveConfig,
    /// Upstream multiplies the conductivity term by the permittivity-weighted mass, which scales
    /// `sigma` by `eps` (~1e-11) and makes it negligible. `false` (default) keeps that behaviour
    /// for parity; `true` assembles the physical `j omega sigma` mass term.
    pub physical_conductivity_term: bool,
}

impl Default for EdgeMaxwell3DConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 1e8,
            relative_permittivity: 1.0,
            relative_permittivity_x: 0.0,
            relative_permittivity_y: 0.0,
            relative_permittivity_z: 0.0,
            relative_permeability: 1.0,
            conductivity_s_per_m: 0.0,
            linear: ComplexSparseSolveConfig::default(),
            physical_conductivity_term: false,
        }
    }
}

/// Impressed current on a directed mesh edge (positive from `node_a` to `node_b`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeCurrentPort3D {
    pub node_a: usize,
    pub node_b: usize,
    pub current_a: C,
}

#[derive(Clone, Debug, Default)]
pub struct EdgeMaxwell3DResult {
    pub edges: Vec<[usize; 2]>,
    pub boundary_edge: Vec<bool>,
    pub edge_voltage_v: Vec<C>,
    pub electric_centroid_v_per_m: Vec<ComplexVec3>,
    pub magnetic_centroid_a_per_m: Vec<ComplexVec3>,
    pub linear_result: IterativeSolverResult,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgePortThevenin {
    pub open_circuit_voltage_v: C,
    pub port_impedance_ohm: C,
}

struct Topology3D {
    edges: Vec<[usize; 2]>,
    edge_index: BTreeMap<[usize; 2], usize>,
    tet_edges: Vec<[usize; 6]>,
    tet_signs: Vec<[f64; 6]>,
    boundary_edge: Vec<bool>,
    /// Boundary face -> owning element.
    boundary_faces: BTreeMap<[usize; 3], usize>,
}

fn sorted3(mut f: [usize; 3]) -> [usize; 3] {
    f.sort_unstable();
    f
}

fn build_topology(mesh: &Mesh3D) -> Result<Topology3D> {
    let mut t = Topology3D {
        edges: Vec::new(),
        edge_index: BTreeMap::new(),
        tet_edges: vec![[0; 6]; mesh.tetrahedra.len()],
        tet_signs: vec![[0.0; 6]; mesh.tetrahedra.len()],
        boundary_edge: Vec::new(),
        boundary_faces: BTreeMap::new(),
    };
    let mut faces: BTreeMap<[usize; 3], (usize, usize)> = BTreeMap::new();
    for (e, tet) in mesh.tetrahedra.iter().enumerate() {
        for le in 0..6 {
            let (edge, sign) =
                canonical_edge(tet.node[TET_EDGES[le][0]], tet.node[TET_EDGES[le][1]]);
            let edges = &mut t.edges;
            let index = *t.edge_index.entry(edge).or_insert_with(|| {
                edges.push(edge);
                edges.len() - 1
            });
            t.tet_edges[e][le] = index;
            t.tet_signs[e][le] = sign;
        }
        for lf in TET_FACES {
            let face = sorted3([tet.node[lf[0]], tet.node[lf[1]], tet.node[lf[2]]]);
            faces.entry(face).and_modify(|v| v.1 += 1).or_insert((e, 1));
        }
    }
    t.boundary_edge = vec![false; t.edges.len()];
    for (face, (owner, count)) in faces {
        if count != 1 {
            continue;
        }
        t.boundary_faces.insert(face, owner);
        for i in 0..3 {
            for j in i + 1..3 {
                let idx = *t
                    .edge_index
                    .get(&canonical_edge(face[i], face[j]).0)
                    .ok_or_else(|| {
                        SolverError::NumericalFailure("inconsistent Tet4 topology".into())
                    })?;
                t.boundary_edge[idx] = true;
            }
        }
    }
    Ok(t)
}

fn scale_add(mut out: ComplexVec3, value: C, basis: RealVec3) -> ComplexVec3 {
    for d in 0..3 {
        out[d] += value * basis[d];
    }
    out
}

/// Driven 3-D curl-curl Maxwell with PEC exterior faces, optional volume current and edge
/// current ports; `exp(+j omega t)` phasors.
pub fn solve_driven_edge_maxwell_3d(
    mesh: &Mesh3D,
    config: &EdgeMaxwell3DConfig,
    impressed_current_density_a_per_m2: Option<&dyn Fn(Point3) -> ComplexVec3>,
    edge_ports: &[EdgeCurrentPort3D],
) -> Result<EdgeMaxwell3DResult> {
    mesh.validate()?;
    let c = config;
    require(
        c.frequency_hz > 0.0
            && c.relative_permittivity > 0.0
            && c.relative_permeability > 0.0
            && c.conductivity_s_per_m >= 0.0
            && c.relative_permittivity_x >= 0.0
            && c.relative_permittivity_y >= 0.0
            && c.relative_permittivity_z >= 0.0,
        "invalid 3-D edge-Maxwell configuration",
    )?;
    let topo = build_topology(mesh)?;
    let mut free_index = vec![usize::MAX; topo.edges.len()];
    let mut free_count = 0;
    for e in 0..topo.edges.len() {
        if !topo.boundary_edge[e] {
            free_index[e] = free_count;
            free_count += 1;
        }
    }
    require(
        free_count > 0,
        "3-D edge-Maxwell mesh has no free interior edges",
    )?;
    let mut matrix = ComplexCsrBuilder::new(free_count)?;
    let mut rhs = vec![ZERO; free_count];
    let omega = 2.0 * PI * c.frequency_hz;
    let mu = MU0 * c.relative_permeability;
    let pick = |v: f64| EPS0 * if v > 0.0 { v } else { c.relative_permittivity };
    let (ex, ey, ez) = (
        pick(c.relative_permittivity_x),
        pick(c.relative_permittivity_y),
        pick(c.relative_permittivity_z),
    );
    let (a, b) = (0.5854101966249685, 0.1381966011250105);
    let quadrature = [[a, b, b, b], [b, a, b, b], [b, b, a, b], [b, b, b, a]];
    for (element, tet) in mesh.tetrahedra.iter().enumerate() {
        let g = nedelec_tet4_geometry(mesh, tet)?;
        let mut mass = [[0.0; 6]; 6];
        let mut plain_mass = [[0.0; 6]; 6];
        let mut load = [ZERO; 6];
        for lambda in quadrature {
            let mut basis = [[0.0; 3]; 6];
            for le in 0..6 {
                basis[le] = nedelec_tet4_basis(&g, le, lambda)?;
            }
            let w = g.volume_m3 / 4.0;
            let current = match impressed_current_density_a_per_m2 {
                Some(f) => {
                    let mut p = Point3::new(0.0, 0.0, 0.0);
                    for i in 0..4 {
                        let q = mesh.nodes[tet.node[i]];
                        p.x += lambda[i] * q.x;
                        p.y += lambda[i] * q.y;
                        p.z += lambda[i] * q.z;
                    }
                    f(p)
                }
                None => [ZERO; 3],
            };
            for i in 0..6 {
                load[i] += (current[0] * basis[i][0]
                    + current[1] * basis[i][1]
                    + current[2] * basis[i][2])
                    * w;
                for j in 0..6 {
                    mass[i][j] += w
                        * (ex * basis[i][0] * basis[j][0]
                            + ey * basis[i][1] * basis[j][1]
                            + ez * basis[i][2] * basis[j][2]);
                    plain_mass[i][j] += w * dot3(basis[i], basis[j]);
                }
            }
        }
        let mut curl = [[0.0; 3]; 6];
        for le in 0..6 {
            curl[le] = nedelec_tet4_curl(&g, le)?;
        }
        for i in 0..6 {
            let gi = topo.tet_edges[element][i];
            if topo.boundary_edge[gi] {
                continue;
            }
            let ri = free_index[gi];
            let si = topo.tet_signs[element][i];
            rhs[ri] += C::new(0.0, -omega) * load[i] * si;
            for j in 0..6 {
                let gj = topo.tet_edges[element][j];
                if topo.boundary_edge[gj] {
                    continue;
                }
                let sj = topo.tet_signs[element][j];
                let local = if c.physical_conductivity_term {
                    C::new(
                        g.volume_m3 * dot3(curl[i], curl[j]) / mu - omega * omega * mass[i][j],
                        omega * c.conductivity_s_per_m * plain_mass[i][j],
                    )
                } else {
                    C::new(g.volume_m3 * dot3(curl[i], curl[j]) / mu, 0.0)
                        + C::new(-omega * omega, omega * c.conductivity_s_per_m) * mass[i][j]
                };
                matrix.add(ri, free_index[gj], local * (si * sj))?;
            }
        }
    }
    for port in edge_ports {
        require(
            port.node_a != port.node_b
                && port.node_a < mesh.node_count()
                && port.node_b < mesh.node_count(),
            "invalid 3-D edge current port",
        )?;
        let (edge, sign) = canonical_edge(port.node_a, port.node_b);
        let idx = *topo.edge_index.get(&edge).ok_or_else(|| {
            SolverError::InvalidInput("3-D edge current port is not a mesh edge".into())
        })?;
        require(
            !topo.boundary_edge[idx],
            "3-D edge current port cannot excite a PEC boundary edge",
        )?;
        rhs[free_index[idx]] += C::new(0.0, -omega) * port.current_a * sign;
    }
    let solved = solve_complex_sparse(&matrix, &rhs, c.linear)?;
    if !solved.linear_result.converged {
        return Err(SolverError::NumericalFailure(
            "3-D edge-Maxwell GMRES failed to converge".into(),
        ));
    }
    let mut voltage = vec![ZERO; topo.edges.len()];
    for e in 0..topo.edges.len() {
        if !topo.boundary_edge[e] {
            voltage[e] = solved.solution[free_index[e]];
        }
    }
    let mut electric = vec![[ZERO; 3]; mesh.tetrahedra.len()];
    let mut magnetic = vec![[ZERO; 3]; mesh.tetrahedra.len()];
    let factor = C::new(0.0, 1.0 / (omega * mu));
    for (element, tet) in mesh.tetrahedra.iter().enumerate() {
        let g = nedelec_tet4_geometry(mesh, tet)?;
        let (mut e_field, mut curl_e) = ([ZERO; 3], [ZERO; 3]);
        for le in 0..6 {
            let dof = voltage[topo.tet_edges[element][le]] * topo.tet_signs[element][le];
            e_field = scale_add(e_field, dof, nedelec_tet4_basis(&g, le, [0.25; 4])?);
            curl_e = scale_add(curl_e, dof, nedelec_tet4_curl(&g, le)?);
        }
        electric[element] = e_field;
        magnetic[element] = curl_e.map(|v| factor * v);
    }
    Ok(EdgeMaxwell3DResult {
        edges: topo.edges,
        boundary_edge: topo.boundary_edge,
        edge_voltage_v: voltage,
        electric_centroid_v_per_m: electric,
        magnetic_centroid_a_per_m: magnetic,
        linear_result: solved.linear_result,
    })
}

/// Small-signal Thevenin point `(V_edge, V_edge / I)` of a driven edge port.
pub fn edge_port_thevenin(
    field: &EdgeMaxwell3DResult,
    edge_index: usize,
    impressed_current_a: C,
) -> Result<EdgePortThevenin> {
    let i = impressed_current_a;
    require(
        edge_index < field.edge_voltage_v.len()
            && i.norm() > 1e-30
            && i.re.is_finite()
            && i.im.is_finite(),
        "invalid edge-port Thevenin request",
    )?;
    let v = field.edge_voltage_v[edge_index];
    Ok(EdgePortThevenin {
        open_circuit_voltage_v: v,
        port_impedance_ohm: v / i,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ResonatorQuality3D {
    pub electric_energy_j: f64,
    pub magnetic_energy_j: f64,
    pub dielectric_loss_w: f64,
    pub conductor_loss_w: f64,
    pub quality_factor: f64,
}

/// Stored energies, volume loss and optional wall loss `R_s = sqrt(pi f mu / sigma_wall)`
/// (upstream defaults `wall_sigma = 0`, `wall_mu_r = 1`).
pub fn resonator_quality_3d(
    mesh: &Mesh3D,
    config: &EdgeMaxwell3DConfig,
    field: &EdgeMaxwell3DResult,
    wall_conductivity_s_per_m: f64,
    wall_relative_permeability: f64,
) -> Result<ResonatorQuality3D> {
    mesh.validate()?;
    require(
        field.electric_centroid_v_per_m.len() == mesh.element_count()
            && field.magnetic_centroid_a_per_m.len() == mesh.element_count(),
        "3-D resonator field/mesh size mismatch",
    )?;
    let c = config;
    require(
        c.frequency_hz > 0.0
            && c.relative_permittivity > 0.0
            && c.relative_permeability > 0.0
            && c.conductivity_s_per_m >= 0.0
            && wall_conductivity_s_per_m >= 0.0
            && wall_relative_permeability > 0.0,
        "invalid resonator quality controls",
    )?;
    let (eps, mu, omega) = (
        EPS0 * c.relative_permittivity,
        MU0 * c.relative_permeability,
        2.0 * PI * c.frequency_hz,
    );
    let mut out = ResonatorQuality3D::default();
    for (e, tet) in mesh.tetrahedra.iter().enumerate() {
        let volume = nedelec_tet4_geometry(mesh, tet)?.volume_m3;
        let (e2, h2) = (
            cnorm2(field.electric_centroid_v_per_m[e]),
            cnorm2(field.magnetic_centroid_a_per_m[e]),
        );
        out.electric_energy_j += 0.25 * eps * e2 * volume;
        out.magnetic_energy_j += 0.25 * mu * h2 * volume;
        out.dielectric_loss_w += 0.5 * c.conductivity_s_per_m * e2 * volume;
    }
    if wall_conductivity_s_per_m > 0.0 {
        let topo = build_topology(mesh)?;
        let rs = (PI * c.frequency_hz * MU0 * wall_relative_permeability
            / wall_conductivity_s_per_m)
            .sqrt();
        for (face, owner) in &topo.boundary_faces {
            let (p0, p1, p2) = (
                mesh.nodes[face[0]],
                mesh.nodes[face[1]],
                mesh.nodes[face[2]],
            );
            let raw = cross3(sub3(p1, p0), sub3(p2, p0));
            let twice_area = norm3(raw);
            if twice_area.is_nan() || twice_area <= 0.0 {
                continue;
            }
            let n = raw.map(|v| v / twice_area);
            let h = field.magnetic_centroid_a_per_m[*owner];
            let hn = h[0] * n[0] + h[1] * n[1] + h[2] * n[2];
            let ht: ComplexVec3 = std::array::from_fn(|d| h[d] - hn * n[d]);
            out.conductor_loss_w += 0.5 * rs * cnorm2(ht) * (0.5 * twice_area);
        }
    }
    let (loss, stored) = (
        out.dielectric_loss_w + out.conductor_loss_w,
        out.electric_energy_j + out.magnetic_energy_j,
    );
    out.quality_factor = if loss > 0.0 {
        omega * stored / loss
    } else {
        f64::INFINITY
    };
    Ok(out)
}

fn face_owners(mesh: &Mesh3D) -> BTreeMap<[usize; 3], Vec<usize>> {
    let mut owners: BTreeMap<[usize; 3], Vec<usize>> = BTreeMap::new();
    for (e, t) in mesh.tetrahedra.iter().enumerate() {
        for lf in TET_FACES {
            owners
                .entry(sorted3([t.node[lf[0]], t.node[lf[1]], t.node[lf[2]]]))
                .or_default()
                .push(e);
        }
    }
    owners
}

fn tangential_jump2(a: ComplexVec3, b: ComplexVec3, n: RealVec3) -> f64 {
    let mut d: ComplexVec3 = std::array::from_fn(|k| a[k] - b[k]);
    let dn = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
    for k in 0..3 {
        d[k] -= dn * n[k];
    }
    cnorm2(d)
}

/// Per-element `sqrt(sum_faces area/2 (|[E_t]|^2 + Z^2 |[H_t]|^2) / 2)` over interior faces
/// (upstream default `Z` = 376.730313668 ohm).
pub fn maxwell_face_jump_indicators_3d(
    mesh: &Mesh3D,
    field: &EdgeMaxwell3DResult,
    impedance_scale_ohm: f64,
) -> Result<Vec<f64>> {
    mesh.validate()?;
    require(
        field.electric_centroid_v_per_m.len() == mesh.element_count()
            && field.magnetic_centroid_a_per_m.len() == mesh.element_count()
            && impedance_scale_ohm > 0.0,
        "invalid Maxwell jump-estimator inputs",
    )?;
    let mut indicator = vec![0.0; mesh.element_count()];
    for (face, owners) in face_owners(mesh) {
        if owners.len() != 2 {
            continue;
        }
        let (p0, p1, p2) = (
            mesh.nodes[face[0]],
            mesh.nodes[face[1]],
            mesh.nodes[face[2]],
        );
        let raw = cross3(sub3(p1, p0), sub3(p2, p0));
        let twice_area = norm3(raw);
        if twice_area.is_nan() || twice_area <= 0.0 {
            continue;
        }
        let n = raw.map(|v| v / twice_area);
        let e_jump = tangential_jump2(
            field.electric_centroid_v_per_m[owners[0]],
            field.electric_centroid_v_per_m[owners[1]],
            n,
        );
        let h_jump = tangential_jump2(
            field.magnetic_centroid_a_per_m[owners[0]],
            field.magnetic_centroid_a_per_m[owners[1]],
            n,
        );
        let eta = 0.5 * twice_area * (e_jump + impedance_scale_ohm * impedance_scale_ohm * h_jump);
        indicator[owners[0]] += 0.5 * eta;
        indicator[owners[1]] += 0.5 * eta;
    }
    Ok(indicator.iter().map(|v| v.max(0.0).sqrt()).collect())
}

/// Smallest descending prefix reaching `fraction` of the squared total (upstream default 0.5);
/// an all-zero indicator marks nothing.
pub fn mark_maxwell_dorfler(indicators: &[f64], fraction: f64) -> Result<Vec<bool>> {
    require(
        !indicators.is_empty() && fraction > 0.0 && fraction <= 1.0,
        "invalid Maxwell Dorfler controls",
    )?;
    require(
        indicators.iter().all(|v| *v >= 0.0 && v.is_finite()),
        "invalid Maxwell error indicator",
    )?;
    let total: f64 = indicators.iter().map(|v| v * v).sum();
    let mut marked = vec![false; indicators.len()];
    if total <= 0.0 {
        return Ok(marked);
    }
    let mut order: Vec<usize> = (0..indicators.len()).collect();
    order.sort_by(|a, b| indicators[*b].total_cmp(&indicators[*a]));
    let mut selected = 0.0;
    for i in order {
        marked[i] = true;
        selected += indicators[i] * indicators[i];
        if selected >= fraction * total {
            break;
        }
    }
    Ok(marked)
}

fn boundary_edge_set(mesh: &Mesh3D) -> BTreeSet<[usize; 2]> {
    let mut edges = BTreeSet::new();
    for (face, owners) in face_owners(mesh) {
        if owners.len() == 1 {
            for i in 0..3 {
                for j in i + 1..3 {
                    edges.insert(canonical_edge(face[i], face[j]).0);
                }
            }
        }
    }
    edges
}

/// Bisects every tetrahedron in the star of each marked element's longest edge (first longest on
/// ties), so shared faces stay conforming.
pub fn refine_tet4_marked_longest_edges(mesh: &Mesh3D, marked: &[bool]) -> Result<Mesh3D> {
    mesh.validate()?;
    require(
        marked.len() == mesh.element_count(),
        "Tet4 refinement marking size mismatch",
    )?;
    let mut selected = BTreeSet::new();
    for (tet, _) in mesh.tetrahedra.iter().zip(marked).filter(|(_, m)| **m) {
        let mut best = -1.0;
        let mut longest = [0, 0];
        for le in TET_EDGES {
            let edge = canonical_edge(tet.node[le[0]], tet.node[le[1]]).0;
            let d = sub3(mesh.nodes[edge[0]], mesh.nodes[edge[1]]);
            let l2 = dot3(d, d);
            if l2 > best {
                best = l2;
                longest = edge;
            }
        }
        selected.insert(longest);
    }
    if selected.is_empty() {
        return Ok(mesh.clone());
    }
    let mut out = mesh.clone();
    for edge in selected {
        let incident: Vec<usize> = out
            .tetrahedra
            .iter()
            .enumerate()
            .filter(|(_, t)| t.node.contains(&edge[0]) && t.node.contains(&edge[1]))
            .map(|(e, _)| e)
            .collect();
        if incident.is_empty() {
            continue;
        }
        let boundary = boundary_edge_set(&out);
        let (a, b) = (out.nodes[edge[0]], out.nodes[edge[1]]);
        let mid = out.nodes.len();
        out.nodes.push(Point3::new(
            0.5 * (a.x + b.x),
            0.5 * (a.y + b.y),
            0.5 * (a.z + b.z),
        ));
        out.boundary_node.push(boundary.contains(&edge));
        let mut rebuilt = Vec::with_capacity(out.tetrahedra.len() + incident.len());
        for (e, t) in out.tetrahedra.iter().enumerate() {
            if incident.binary_search(&e).is_err() {
                rebuilt.push(*t);
                continue;
            }
            let others: Vec<usize> = t
                .node
                .iter()
                .copied()
                .filter(|n| *n != edge[0] && *n != edge[1])
                .collect();
            require(others.len() == 2, "invalid Tet4 edge-star refinement")?;
            rebuilt.push(Tet4 {
                node: [edge[0], mid, others[0], others[1]],
            });
            rebuilt.push(Tet4 {
                node: [mid, edge[1], others[0], others[1]],
            });
        }
        out.tetrahedra = rebuilt;
    }
    out.validate()?;
    Ok(out)
}
