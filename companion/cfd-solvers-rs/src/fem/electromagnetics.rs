//! 2D electrostatics and DC conduction on Tri3 meshes, ported from upstream
//! `cfd/solvers/fem/electromagnetics2d` (MIT, cfd_solvers contributors).
//! Both solve `-div(k grad phi) = f` with the shared scalar-diffusion assembly and
//! report element-centroid gradients of the linear potential.
use super::*;
use crate::{core::LinearConfig, require, Result};
use std::collections::BTreeMap;

/// Element-wise 2D vector (field or current density).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vector2 {
    pub x: f64,
    pub y: f64,
}

fn element_gradients(mesh: &Mesh2D, potential: &[f64]) -> Result<Vec<Vector2>> {
    require(potential.len() == mesh.node_count(), "potential size")?;
    Ok(mesh
        .triangles
        .iter()
        .map(|tri| {
            let (_, g, _) = tri_geometry(mesh, tri);
            let mut out = Vector2::default();
            for (k, node) in tri.node.iter().enumerate() {
                out.x += potential[*node] * g[k].x;
                out.y += potential[*node] * g[k].y;
            }
            out
        })
        .collect())
}

fn centroid(mesh: &Mesh2D, tri: &Tri3) -> Node2 {
    let [a, b, c] = tri.node.map(|i| mesh.nodes[i]);
    Node2 {
        x: (a.x + b.x + c.x) / 3.0,
        y: (a.y + b.y + c.y) / 3.0,
    }
}

/// Shared state of the two potential solvers.
pub struct PotentialSolver2D {
    mesh: Mesh2D,
    boundary: BTreeMap<i32, ScalarBoundaryCondition2D<'static>>,
    config: LinearConfig,
    potential: Vec<f64>,
    linear: Option<crate::core::IterativeSolverResult>,
}

impl PotentialSolver2D {
    fn new(mesh: Mesh2D, config: LinearConfig) -> Result<Self> {
        mesh.validate()?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            boundary: BTreeMap::new(),
            config,
            potential: vec![0.0; n],
            linear: None,
        })
    }
    fn solve(&mut self, k: &dyn Fn(Node2) -> f64, f: &dyn Fn(Node2) -> f64) -> Result<()> {
        let s = scalar_diffusion2d(&self.mesh, k, |_| 0.0, f, &self.boundary, self.config)?;
        self.potential = s.values;
        self.linear = Some(s.linear);
        Ok(())
    }
}

/// Solves `-div(eps_r grad phi) = rho/eps0`. Supplying `rho/eps0` keeps SI systems well scaled.
pub struct Electrostatics2D(PotentialSolver2D);

impl Electrostatics2D {
    pub fn new(mesh: Mesh2D, config: LinearConfig) -> Result<Self> {
        Ok(Self(PotentialSolver2D::new(mesh, config)?))
    }
    pub fn set_boundary(&mut self, patch: i32, condition: ScalarBoundaryCondition2D<'static>) {
        self.0.boundary.insert(patch, condition);
    }
    pub fn solve(
        &mut self,
        relative_permittivity: &dyn Fn(Node2) -> f64,
        charge_over_epsilon0: &dyn Fn(Node2) -> f64,
    ) -> Result<()> {
        self.0.solve(relative_permittivity, charge_over_epsilon0)
    }
    pub fn potential(&self) -> &[f64] {
        &self.0.potential
    }
    pub fn mesh(&self) -> &Mesh2D {
        &self.0.mesh
    }
    pub fn linear_result(&self) -> Option<crate::core::IterativeSolverResult> {
        self.0.linear
    }
    /// Element electric field `E = -grad(phi)`.
    pub fn element_electric_field(&self) -> Result<Vec<Vector2>> {
        Ok(element_gradients(&self.0.mesh, &self.0.potential)?
            .into_iter()
            .map(|g| Vector2 { x: -g.x, y: -g.y })
            .collect())
    }
}

/// Solves `-div(sigma grad phi) = source` for steady current flow.
pub struct DcConduction2D(PotentialSolver2D);

impl DcConduction2D {
    pub fn new(mesh: Mesh2D, config: LinearConfig) -> Result<Self> {
        Ok(Self(PotentialSolver2D::new(mesh, config)?))
    }
    pub fn set_boundary(&mut self, patch: i32, condition: ScalarBoundaryCondition2D<'static>) {
        self.0.boundary.insert(patch, condition);
    }
    pub fn solve(
        &mut self,
        conductivity: &dyn Fn(Node2) -> f64,
        current_source: Option<&dyn Fn(Node2) -> f64>,
    ) -> Result<()> {
        let zero = |_: Node2| 0.0;
        self.0.solve(conductivity, current_source.unwrap_or(&zero))
    }
    pub fn potential(&self) -> &[f64] {
        &self.0.potential
    }
    pub fn mesh(&self) -> &Mesh2D {
        &self.0.mesh
    }
    pub fn linear_result(&self) -> Option<crate::core::IterativeSolverResult> {
        self.0.linear
    }
    /// Element current density `J = -sigma grad(phi)`, sigma sampled at the centroid.
    pub fn element_current_density(
        &self,
        conductivity: &dyn Fn(Node2) -> f64,
    ) -> Result<Vec<Vector2>> {
        let g = element_gradients(&self.0.mesh, &self.0.potential)?;
        let mut out = Vec::with_capacity(g.len());
        for (tri, g) in self.0.mesh.triangles.iter().zip(g) {
            let sigma = conductivity(centroid(&self.0.mesh, tri));
            require(sigma.is_finite() && sigma > 0.0, "invalid DC conductivity")?;
            out.push(Vector2 {
                x: -sigma * g.x,
                y: -sigma * g.y,
            });
        }
        Ok(out)
    }
}
