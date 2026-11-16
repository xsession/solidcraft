use super::*;
use crate::{
    core::{
        conjugate_gradient, solve_tridiagonal, CsrBuilder, CsrMatrix, IterativeSolverResult,
        LinearConfig,
    },
    finite_slice, positive, require, Result, SolverError,
};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct Poisson1D {
    pub x: Vec<f64>,
    pub u: Vec<f64>,
    length: f64,
}
impl Poisson1D {
    pub fn new(elements: usize, length: f64) -> Result<Self> {
        require(elements >= 2, "Poisson1D needs at least two elements")?;
        positive(length, "length")?;
        require(elements < usize::MAX, "element count overflow")?;
        Ok(Self {
            x: (0..=elements)
                .map(|i| length * i as f64 / elements as f64)
                .collect(),
            u: vec![0.0; elements + 1],
            length,
        })
    }
    pub fn solve(&mut self, source: impl Fn(f64) -> f64, left: f64, right: f64) -> Result<()> {
        let n = self.x.len() - 2;
        let h = self.length / (n + 1) as f64;
        let mut rhs = vec![0.0; n];
        for j in 0..n {
            rhs[j] = 0.5 * h * (source(self.x[j + 1] - 0.5 * h) + source(self.x[j + 1] + 0.5 * h));
        }
        rhs[0] += left / h;
        rhs[n - 1] += right / h;
        let u = solve_tridiagonal(
            &vec![-1.0 / h; n],
            &vec![2.0 / h; n],
            &vec![-1.0 / h; n],
            &rhs,
        )?;
        self.u[0] = left;
        self.u[n + 1] = right;
        self.u[1..=n].copy_from_slice(&u);
        finite_slice(&self.u, "Poisson1D result")
    }
}
#[derive(Clone, Debug)]
pub struct ScalarSolution {
    pub values: Vec<f64>,
    pub linear: IterativeSolverResult,
}
fn solve_prescribed(
    builder: &mut CsrBuilder,
    rhs: &mut [f64],
    prescribed: &[Option<f64>],
    initial: &[f64],
    config: LinearConfig,
) -> Result<ScalarSolution> {
    for (i, v) in prescribed.iter().enumerate() {
        if let Some(v) = v {
            builder.dirichlet(rhs, i, *v)?;
        }
    }
    let a = builder.build(0.0)?;
    let mut values = initial.to_vec();
    let linear = conjugate_gradient(&a, rhs, &mut values, config)?;
    if !linear.converged {
        return Err(SolverError::NumericalFailure(format!(
            "FEM PCG failed after {} iterations, RMS {}",
            linear.iterations, linear.residual_rms
        )));
    }
    finite_slice(&values, "FEM solution")?;
    Ok(ScalarSolution { values, linear })
}
pub(crate) fn tri_geometry(mesh: &Mesh2D, tri: &Tri3) -> (f64, [Point3; 3], Node2) {
    let [a, b, c] = tri.node.map(|i| mesh.nodes[i]);
    let d = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    (
        0.5 * d.abs(),
        [
            Point3::new((b.y - c.y) / d, (c.x - b.x) / d, 0.0),
            Point3::new((c.y - a.y) / d, (a.x - c.x) / d, 0.0),
            Point3::new((a.y - b.y) / d, (b.x - a.x) / d, 0.0),
        ],
        Node2 {
            x: (a.x + b.x + c.x) / 3.0,
            y: (a.y + b.y + c.y) / 3.0,
        },
    )
}
/// P1 Poisson with centroid source quadrature, as in the upstream Poisson2D.
pub fn poisson2d(
    mesh: &Mesh2D,
    diffusivity: f64,
    source: impl Fn(Node2) -> f64,
    boundary: impl Fn(Node2) -> f64,
    config: LinearConfig,
) -> Result<ScalarSolution> {
    mesh.validate()?;
    positive(diffusivity, "diffusivity")?;
    let n = mesh.node_count();
    let mut b = CsrBuilder::new(n, n);
    let mut rhs = vec![0.0; n];
    for t in &mesh.triangles {
        let (area, g, c) = tri_geometry(mesh, t);
        let f = source(c);
        require(f.is_finite(), "nonfinite source")?;
        for i in 0..3 {
            rhs[t.node[i]] += f * area / 3.0;
            for j in 0..3 {
                b.add(t.node[i], t.node[j], diffusivity * area * g[i].dot(g[j]))?;
            }
        }
    }
    let prescribed: Vec<_> = mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(i, &p)| {
            if mesh.boundary_node[i] {
                Some(boundary(p))
            } else {
                None
            }
        })
        .collect();
    solve_prescribed(&mut b, &mut rhs, &prescribed, &vec![0.0; n], config)
}
pub fn poisson3d(
    mesh: &Mesh3D,
    diffusivity: f64,
    source: impl Fn(Point3) -> f64,
    boundary: impl Fn(Point3) -> f64,
    config: LinearConfig,
) -> Result<ScalarSolution> {
    mesh.validate()?;
    positive(diffusivity, "diffusivity")?;
    let n = mesh.node_count();
    let mut b = CsrBuilder::new(n, n);
    let mut rhs = vec![0.0; n];
    for t in &mesh.tetrahedra {
        let nodes = t.node.map(|i| mesh.nodes[i]);
        for q in gaussian_quadrature(ElementType::Tet4, 2)? {
            let s = evaluate_shape(ElementType::Tet4, q.point);
            let iso = evaluate_isoparametric(ElementType::Tet4, &nodes, q.point)?;
            let f = source(iso.point);
            require(f.is_finite(), "nonfinite source")?;
            let w = q.weight * iso.measure;
            for i in 0..4 {
                rhs[t.node[i]] += f * s.value[i] * w;
                for j in 0..4 {
                    b.add(
                        t.node[i],
                        t.node[j],
                        diffusivity * iso.gradient_physical[i].dot(iso.gradient_physical[j]) * w,
                    )?;
                }
            }
        }
    }
    let prescribed: Vec<_> = mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(i, &p)| {
            if mesh.boundary_node[i] {
                Some(boundary(p))
            } else {
                None
            }
        })
        .collect();
    solve_prescribed(&mut b, &mut rhs, &prescribed, &vec![0.0; n], config)
}
pub fn assemble_tri3_laplace_matrix(
    mesh: &Mesh2D,
    diffusivity: impl Fn(Node2) -> f64,
) -> Result<CsrMatrix> {
    mesh.validate()?;
    let mut b = CsrBuilder::new(mesh.node_count(), mesh.node_count());
    for t in &mesh.triangles {
        let nodes = t.node.map(|i| mesh.nodes[i].point3());
        for q in gaussian_quadrature(ElementType::Tri3, 2)? {
            let iso = evaluate_isoparametric(ElementType::Tri3, &nodes, q.point)?;
            let k = diffusivity(Node2 {
                x: iso.point.x,
                y: iso.point.y,
            });
            positive(k, "diffusivity")?;
            for i in 0..3 {
                for j in 0..3 {
                    b.add(
                        t.node[i],
                        t.node[j],
                        k * iso.gradient_physical[i].dot(iso.gradient_physical[j])
                            * q.weight
                            * iso.measure,
                    )?;
                }
            }
        }
    }
    b.build(0.0)
}
pub fn apply_tri3_laplace_matrix_free(
    mesh: &Mesh2D,
    diffusivity: impl Fn(Node2) -> f64,
    x: &[f64],
    y: &mut [f64],
) -> Result<()> {
    mesh.validate()?;
    require(
        x.len() == mesh.node_count() && y.len() == x.len(),
        "matrix-free FEM dimensions",
    )?;
    finite_slice(x, "matrix-free FEM input")?;
    y.fill(0.0);
    for t in &mesh.triangles {
        let nodes = t.node.map(|i| mesh.nodes[i].point3());
        for q in gaussian_quadrature(ElementType::Tri3, 2)? {
            let iso = evaluate_isoparametric(ElementType::Tri3, &nodes, q.point)?;
            let k = diffusivity(Node2 {
                x: iso.point.x,
                y: iso.point.y,
            });
            positive(k, "diffusivity")?;
            for i in 0..3 {
                for j in 0..3 {
                    y[t.node[i]] += k
                        * iso.gradient_physical[i].dot(iso.gradient_physical[j])
                        * q.weight
                        * iso.measure
                        * x[t.node[j]];
                }
            }
        }
    }
    Ok(())
}
pub enum ScalarBoundaryCondition2D<'a> {
    Dirichlet(Box<dyn Fn(Node2) -> f64 + 'a>),
    Neumann(Box<dyn Fn(Node2) -> f64 + 'a>),
    Robin {
        value: Box<dyn Fn(Node2) -> f64 + 'a>,
        alpha: Box<dyn Fn(Node2) -> f64 + 'a>,
    },
}
/// General -div(k grad u) + c u = f, with per-patch conditions.
pub fn scalar_diffusion2d(
    mesh: &Mesh2D,
    diffusivity: impl Fn(Node2) -> f64,
    reaction: impl Fn(Node2) -> f64,
    source: impl Fn(Node2) -> f64,
    boundary: &BTreeMap<i32, ScalarBoundaryCondition2D<'_>>,
    config: LinearConfig,
) -> Result<ScalarSolution> {
    mesh.validate()?;
    let n = mesh.node_count();
    let mut b = CsrBuilder::new(n, n);
    let mut rhs = vec![0.0; n];
    let mut prescribed: Vec<Option<f64>> = vec![None; n];
    for e in &mesh.boundary_edges {
        if let Some(ScalarBoundaryCondition2D::Dirichlet(value)) = boundary.get(&e.patch) {
            for i in e.node {
                let v = value(mesh.nodes[i]);
                require(v.is_finite(), "nonfinite Dirichlet value")?;
                if let Some(old) = prescribed[i] {
                    require(
                        (old - v).abs() <= 1e-10 * v.abs().max(1.0),
                        "conflicting Dirichlet corner values",
                    )?;
                }
                prescribed[i] = Some(v);
            }
        }
    }
    for t in &mesh.triangles {
        let nodes = t.node.map(|i| mesh.nodes[i].point3());
        for q in gaussian_quadrature(ElementType::Tri3, 2)? {
            let iso = evaluate_isoparametric(ElementType::Tri3, &nodes, q.point)?;
            let s = evaluate_shape(ElementType::Tri3, q.point);
            let p = Node2 {
                x: iso.point.x,
                y: iso.point.y,
            };
            let (k, c, f) = (diffusivity(p), reaction(p), source(p));
            positive(k, "diffusivity")?;
            require(
                c.is_finite() && c >= 0.0 && f.is_finite(),
                "invalid reaction/source",
            )?;
            let w = q.weight * iso.measure;
            for i in 0..3 {
                rhs[t.node[i]] += f * s.value[i] * w;
                for j in 0..3 {
                    b.add(
                        t.node[i],
                        t.node[j],
                        (k * iso.gradient_physical[i].dot(iso.gradient_physical[j])
                            + c * s.value[i] * s.value[j])
                            * w,
                    )?;
                }
            }
        }
    }
    for e in &mesh.boundary_edges {
        let Some(bc) = boundary.get(&e.patch) else {
            continue;
        };
        if matches!(bc, ScalarBoundaryCondition2D::Dirichlet(_)) {
            continue;
        }
        let [a, c] = e.node.map(|i| mesh.nodes[i]);
        let length = (a.x - c.x).hypot(a.y - c.y);
        for xi in [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()] {
            let nn = [0.5 * (1.0 - xi), 0.5 * (1.0 + xi)];
            let p = Node2 {
                x: nn[0] * a.x + nn[1] * c.x,
                y: nn[0] * a.y + nn[1] * c.y,
            };
            let (beta, alpha) = match bc {
                ScalarBoundaryCondition2D::Neumann(v) => (v(p), 0.0),
                ScalarBoundaryCondition2D::Robin { value, alpha } => (value(p), alpha(p)),
                ScalarBoundaryCondition2D::Dirichlet(_) => continue,
            };
            require(
                beta.is_finite() && alpha.is_finite() && alpha >= 0.0,
                "invalid boundary coefficient",
            )?;
            for i in 0..2 {
                rhs[e.node[i]] += beta * nn[i] * 0.5 * length;
                for j in 0..2 {
                    b.add(e.node[i], e.node[j], alpha * nn[i] * nn[j] * 0.5 * length)?;
                }
            }
        }
    }
    solve_prescribed(&mut b, &mut rhs, &prescribed, &vec![0.0; n], config)
}
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Heat2DConfig {
    pub conductivity: f64,
    pub volumetric_heat_capacity: f64,
    pub dt: f64,
    pub linear: LinearConfig,
}
impl Default for Heat2DConfig {
    fn default() -> Self {
        Self {
            conductivity: 1.0,
            volumetric_heat_capacity: 1.0,
            dt: 1e-3,
            linear: LinearConfig {
                max_iterations: 2000,
                ..Default::default()
            },
        }
    }
}
#[derive(Clone, Debug)]
pub struct Heat2D {
    pub mesh: Mesh2D,
    pub temperature: Vec<f64>,
    pub time: f64,
    pub linear_result: Option<IterativeSolverResult>,
    config: Heat2DConfig,
}
impl Heat2D {
    pub fn new(mesh: Mesh2D, config: Heat2DConfig) -> Result<Self> {
        mesh.validate()?;
        positive(config.conductivity, "conductivity")?;
        positive(config.volumetric_heat_capacity, "heat capacity")?;
        positive(config.dt, "dt")?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            temperature: vec![0.0; n],
            time: 0.0,
            linear_result: None,
            config,
        })
    }
    pub fn initialize(&mut self, f: impl Fn(Node2) -> f64) -> Result<()> {
        let v: Vec<_> = self.mesh.nodes.iter().map(|&p| f(p)).collect();
        finite_slice(&v, "initial temperature")?;
        self.temperature = v;
        Ok(())
    }
    pub fn step(
        &mut self,
        source: impl Fn(usize, Node2, f64) -> f64,
        boundary: impl Fn(Node2, f64) -> f64,
    ) -> Result<()> {
        self.mesh.validate()?;
        require(
            self.temperature.len() == self.mesh.node_count(),
            "temperature size",
        )?;
        let n = self.mesh.node_count();
        let mut b = CsrBuilder::new(n, n);
        let mut rhs = vec![0.0; n];
        let t = self.time + self.config.dt;
        for (e, tri) in self.mesh.triangles.iter().enumerate() {
            let (area, g, p) = tri_geometry(&self.mesh, tri);
            let q = source(e, p, t);
            require(q.is_finite(), "nonfinite heat source")?;
            for i in 0..3 {
                rhs[tri.node[i]] += q * area / 3.0;
                for j in 0..3 {
                    let mass = self.config.volumetric_heat_capacity
                        * area
                        * (if i == j { 2.0 } else { 1.0 })
                        / 12.0;
                    let stiff = self.config.conductivity * area * g[i].dot(g[j]);
                    b.add(tri.node[i], tri.node[j], mass / self.config.dt + stiff)?;
                    rhs[tri.node[i]] += mass / self.config.dt * self.temperature[tri.node[j]];
                }
            }
        }
        let prescribed: Vec<_> = self
            .mesh
            .nodes
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                if self.mesh.boundary_node[i] {
                    Some(boundary(p, t))
                } else {
                    None
                }
            })
            .collect();
        let s = solve_prescribed(
            &mut b,
            &mut rhs,
            &prescribed,
            &self.temperature,
            self.config.linear,
        )?;
        self.temperature = s.values;
        self.linear_result = Some(s.linear);
        self.time = t;
        Ok(())
    }
    pub fn step_element_source(
        &mut self,
        values: &[f64],
        boundary: impl Fn(Node2, f64) -> f64,
    ) -> Result<()> {
        require(
            values.len() == self.mesh.element_count(),
            "element heat source size",
        )?;
        self.step(|e, _, _| values[e], boundary)
    }
}
