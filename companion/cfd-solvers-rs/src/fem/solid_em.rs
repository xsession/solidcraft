//! Axisymmetric elasticity, 2-D magnetostatics with nonlinear B-H, coils and Maxwell stress,
//! cubic nonlinear Poisson, Green-Lagrange truss/triangle kinematics and enthalpy phase change,
//! ported from upstream `cfd/solvers/fem/{axisymmetric_elasticity, magnetostatics2d,
//! nonlinear_poisson2d}` and `cfd/fem/{nonlinear_geometry, phase_change}` (MIT, cfd_solvers
//! contributors).
// Element loops index several local arrays by node, mirroring the upstream assembly.
#![allow(clippy::needless_range_loop)]
use super::{
    evaluate_isoparametric, evaluate_shape, gaussian_quadrature, scalar_diffusion2d, Displacement2,
    ElementType, Mesh2D, Node2, ReferencePoint, ScalarBoundaryCondition2D, Tri3,
};
use crate::{
    core::{
        newton_solve, preconditioned_conjugate_gradient, CsrBuilder, IterativeSolverResult,
        JacobiPreconditioner, LinearConfig, NewtonConfig, NewtonResult,
    },
    require, Result, SolverError,
};
use std::collections::BTreeMap;

const CENTROID: ReferencePoint = ReferencePoint {
    xi: 1.0 / 3.0,
    eta: 1.0 / 3.0,
    zeta: 0.0,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricElasticityConfig {
    pub youngs_modulus: f64,
    pub poisson_ratio: f64,
    pub max_iterations: usize,
    pub relative_tolerance: f64,
}

impl Default for AxisymmetricElasticityConfig {
    fn default() -> Self {
        Self {
            youngs_modulus: 1.0,
            poisson_ratio: 0.3,
            max_iterations: 5000,
            relative_tolerance: 1e-10,
        }
    }
}

/// Small-strain axisymmetric elasticity; mesh `x` is `r`, `y` is `z`, with the `2 pi r`
/// measure and hoop strain `u_r / r`.
#[derive(Clone, Debug)]
pub struct AxisymmetricElasticity {
    mesh: Mesh2D,
    config: AxisymmetricElasticityConfig,
    displacement: Vec<Displacement2>,
    prescribed_r: Vec<Option<f64>>,
    prescribed_z: Vec<Option<f64>>,
    linear: IterativeSolverResult,
}

impl AxisymmetricElasticity {
    pub fn new(mesh: Mesh2D, config: AxisymmetricElasticityConfig) -> Result<Self> {
        mesh.validate()?;
        require(
            config.youngs_modulus > 0.0
                && config.poisson_ratio > -1.0
                && config.poisson_ratio < 0.5
                && config.max_iterations > 0
                && config.relative_tolerance > 0.0,
            "invalid axisymmetric elasticity controls",
        )?;
        require(
            mesh.nodes.iter().all(|p| p.x >= -1e-14),
            "axisymmetric FEM requires non-negative radial coordinate",
        )?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            config,
            displacement: vec![Displacement2::default(); n],
            prescribed_r: vec![None; n],
            prescribed_z: vec![None; n],
            linear: IterativeSolverResult::default(),
        })
    }
    pub fn set_dirichlet(
        &mut self,
        predicate: impl Fn(Node2) -> bool,
        value: impl Fn(Node2) -> Displacement2,
        fix_radial: bool,
        fix_axial: bool,
    ) {
        for (i, p) in self.mesh.nodes.iter().enumerate() {
            if predicate(*p) {
                let v = value(*p);
                if fix_radial {
                    self.prescribed_r[i] = Some(v.x).filter(|x| x.is_finite());
                }
                if fix_axial {
                    self.prescribed_z[i] = Some(v.y).filter(|x| x.is_finite());
                }
            }
        }
    }
    pub fn solve(&mut self, body_force: Option<&dyn Fn(Node2) -> Displacement2>) -> Result<()> {
        let n = self.mesh.node_count();
        let mut map = vec![usize::MAX; 2 * n];
        let mut nf = 0;
        for i in 0..n {
            match self.prescribed_r[i] {
                Some(v) => self.displacement[i].x = v,
                None => {
                    map[2 * i] = nf;
                    nf += 1;
                }
            }
            match self.prescribed_z[i] {
                Some(v) => self.displacement[i].y = v,
                None => {
                    map[2 * i + 1] = nf;
                    nf += 1;
                }
            }
        }
        require(nf > 0, "axisymmetric elasticity system has no free DOFs")?;
        let (e, nu) = (self.config.youngs_modulus, self.config.poisson_ratio);
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let mut d = [[0.0; 4]; 4];
        for (i, row) in d.iter_mut().enumerate().take(3) {
            for v in row.iter_mut().take(3) {
                *v = lambda;
            }
            row[i] += 2.0 * mu;
        }
        d[3][3] = mu;
        let rule = gaussian_quadrature(ElementType::Tri3, 2)?;
        let mut builder = CsrBuilder::new(nf, nf);
        let mut rhs = vec![0.0; nf];
        for tri in &self.mesh.triangles {
            let nodes = tri.node.map(|i| self.mesh.nodes[i].point3());
            let mut ke = [[0.0; 6]; 6];
            let mut fe = [0.0; 6];
            for qp in &rule {
                let shape = evaluate_shape(ElementType::Tri3, qp.point);
                let iso = evaluate_isoparametric(ElementType::Tri3, &nodes, qp.point)?;
                let r = iso.point.x;
                if r.is_nan() || r <= 1e-14 {
                    continue;
                }
                let mut b = [[0.0; 6]; 4];
                for i in 0..3 {
                    let g = iso.gradient_physical[i];
                    b[0][2 * i] = g.x;
                    b[1][2 * i + 1] = g.y;
                    b[2][2 * i] = shape.value[i] / r;
                    b[3][2 * i] = g.y;
                    b[3][2 * i + 1] = g.x;
                }
                let w = 2.0 * std::f64::consts::PI * r * qp.weight * iso.measure;
                for i in 0..6 {
                    for j in 0..6 {
                        for a in 0..4 {
                            for bb in 0..4 {
                                ke[i][j] += w * b[a][i] * d[a][bb] * b[bb][j];
                            }
                        }
                    }
                }
                let bf = body_force.map_or(Displacement2::default(), |f| {
                    f(Node2 {
                        x: iso.point.x,
                        y: iso.point.y,
                    })
                });
                for i in 0..3 {
                    fe[2 * i] += w * shape.value[i] * bf.x;
                    fe[2 * i + 1] += w * shape.value[i] * bf.y;
                }
            }
            for i in 0..3 {
                for comp in 0..2 {
                    let row = map[2 * tri.node[i] + comp];
                    if row == usize::MAX {
                        continue;
                    }
                    rhs[row] += fe[2 * i + comp];
                    for j in 0..3 {
                        for cj in 0..2 {
                            let kij = ke[2 * i + comp][2 * j + cj];
                            let col = map[2 * tri.node[j] + cj];
                            if col == usize::MAX {
                                let u = self.displacement[tri.node[j]];
                                rhs[row] -= kij * if cj == 0 { u.x } else { u.y };
                            } else {
                                builder.add(row, col, kij)?;
                            }
                        }
                    }
                }
            }
        }
        let a = builder.build(0.0)?;
        let jacobi = JacobiPreconditioner::new(&a.diagonal(), 1e-300)?;
        let mut x = vec![0.0; nf];
        self.linear = preconditioned_conjugate_gradient(
            &rhs,
            &mut x,
            |v, o| a.multiply(v, o),
            &jacobi,
            LinearConfig {
                max_iterations: self.config.max_iterations,
                relative_tolerance: self.config.relative_tolerance,
                ..Default::default()
            },
        )?;
        require(
            self.linear.converged,
            "AxisymmetricElasticity PCG did not converge",
        )?;
        for i in 0..n {
            if map[2 * i] != usize::MAX {
                self.displacement[i].x = x[map[2 * i]];
            }
            if map[2 * i + 1] != usize::MAX {
                self.displacement[i].y = x[map[2 * i + 1]];
            }
        }
        Ok(())
    }
    pub fn mesh(&self) -> &Mesh2D {
        &self.mesh
    }
    pub fn displacement(&self) -> &[Displacement2] {
        &self.displacement
    }
    pub fn linear_result(&self) -> &IterativeSolverResult {
        &self.linear
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MagneticFluxDensity2 {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MagnetostaticNonlinearResult {
    pub iterations: usize,
    pub maximum_relative_reluctivity_change: f64,
    pub converged: bool,
}

/// Out-of-plane vector potential: `-div(nu grad A_z) = J_z`, `B = (dA/dy, -dA/dx)`.
pub struct Magnetostatics2D {
    mesh: Mesh2D,
    config: LinearConfig,
    boundary: BTreeMap<i32, ScalarBoundaryCondition2D<'static>>,
    potential: Vec<f64>,
    linear: IterativeSolverResult,
}

fn centroid(mesh: &Mesh2D, t: &Tri3) -> Node2 {
    let [a, b, c] = t.node.map(|i| mesh.nodes[i]);
    Node2 {
        x: (a.x + b.x + c.x) / 3.0,
        y: (a.y + b.y + c.y) / 3.0,
    }
}

impl Magnetostatics2D {
    pub fn new(mesh: Mesh2D, config: LinearConfig) -> Result<Self> {
        mesh.validate()?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            config,
            boundary: BTreeMap::new(),
            potential: vec![0.0; n],
            linear: IterativeSolverResult::default(),
        })
    }
    pub fn set_boundary(&mut self, patch: i32, condition: ScalarBoundaryCondition2D<'static>) {
        self.boundary.insert(patch, condition);
    }
    pub fn solve(
        &mut self,
        reluctivity: impl Fn(Node2) -> f64,
        current_density_z: impl Fn(Node2) -> f64,
    ) -> Result<()> {
        let s = scalar_diffusion2d(
            &self.mesh,
            reluctivity,
            |_| 0.0,
            current_density_z,
            &self.boundary,
            self.config,
        )?;
        self.potential = s.values;
        self.linear = s.linear;
        Ok(())
    }
    /// Element-wise Picard iteration on `nu(|B|)` with under-relaxation (upstream defaults 50, 1e-6, 0.7).
    pub fn solve_nonlinear(
        &mut self,
        reluctivity_law: impl Fn(f64) -> f64,
        current_density_z: impl Fn(Node2) -> f64,
        max_iterations: usize,
        relative_tolerance: f64,
        relaxation: f64,
    ) -> Result<MagnetostaticNonlinearResult> {
        require(
            max_iterations > 0 && relative_tolerance > 0.0 && relaxation > 0.0 && relaxation <= 1.0,
            "invalid nonlinear magnetostatic controls",
        )?;
        let nu0 = reluctivity_law(0.0);
        require(nu0 > 0.0 && nu0.is_finite(), "invalid B-H reluctivity law")?;
        let mut nu = vec![nu0; self.mesh.triangles.len()];
        let mesh = self.mesh.clone();
        let element_for_point = |p: Node2| {
            let mut nearest = 0;
            let mut best = f64::INFINITY;
            for (e, t) in mesh.triangles.iter().enumerate() {
                let [a, b, c] = t.node.map(|i| mesh.nodes[i]);
                let det = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
                let l0 = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / det;
                let l1 = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / det;
                let l2 = 1.0 - l0 - l1;
                if l0 >= -1e-10 && l1 >= -1e-10 && l2 >= -1e-10 {
                    return e;
                }
                let (cx, cy) = ((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0);
                let d2 = (p.x - cx) * (p.x - cx) + (p.y - cy) * (p.y - cy);
                if d2 < best {
                    best = d2;
                    nearest = e;
                }
            }
            nearest
        };
        let mut result = MagnetostaticNonlinearResult::default();
        for iteration in 0..max_iterations {
            let current = &nu;
            self.solve(|p| current[element_for_point(p)], &current_density_z)?;
            let b = self.element_flux_density()?;
            let mut max_change = 0.0f64;
            let mut updated = vec![0.0; nu.len()];
            for e in 0..nu.len() {
                let target = reluctivity_law(b[e].x.hypot(b[e].y));
                require(
                    target > 0.0 && target.is_finite(),
                    "nonlinear B-H law returned invalid reluctivity",
                )?;
                updated[e] = (1.0 - relaxation) * nu[e] + relaxation * target;
                max_change =
                    max_change.max((updated[e] - nu[e]).abs() / updated[e].abs().max(1e-30));
            }
            nu = updated;
            result = MagnetostaticNonlinearResult {
                iterations: iteration + 1,
                maximum_relative_reluctivity_change: max_change,
                converged: max_change <= relative_tolerance,
            };
            if result.converged {
                break;
            }
        }
        Ok(result)
    }
    pub fn mesh(&self) -> &Mesh2D {
        &self.mesh
    }
    pub fn vector_potential(&self) -> &[f64] {
        &self.potential
    }
    pub fn linear_result(&self) -> &IterativeSolverResult {
        &self.linear
    }
    pub fn element_flux_density(&self) -> Result<Vec<MagneticFluxDensity2>> {
        let mut out = Vec::with_capacity(self.mesh.triangles.len());
        for t in &self.mesh.triangles {
            let nodes = t.node.map(|i| self.mesh.nodes[i].point3());
            let iso = evaluate_isoparametric(ElementType::Tri3, &nodes, CENTROID)?;
            let (mut gx, mut gy) = (0.0, 0.0);
            for i in 0..3 {
                gx += self.potential[t.node[i]] * iso.gradient_physical[i].x;
                gy += self.potential[t.node[i]] * iso.gradient_physical[i].y;
            }
            out.push(MagneticFluxDensity2 { x: gy, y: -gx });
        }
        Ok(out)
    }
}

/// Piecewise-linear B-H curve from the origin with linear extrapolation beyond the last sample.
#[derive(Clone, Debug, PartialEq)]
pub struct PiecewiseLinearBHCurve {
    b: Vec<f64>,
    h: Vec<f64>,
}

impl PiecewiseLinearBHCurve {
    pub fn new(flux_density_t: Vec<f64>, field_strength_a_per_m: Vec<f64>) -> Result<Self> {
        let (b, h) = (flux_density_t, field_strength_a_per_m);
        require(
            b.len() >= 2 && b.len() == h.len() && b[0] == 0.0 && h[0] == 0.0,
            "B-H curve must start at origin and contain at least two samples",
        )?;
        require(
            b.windows(2).all(|w| w[1] > w[0]) && h.windows(2).all(|w| w[1] > w[0]),
            "B-H samples must be strictly increasing",
        )?;
        Ok(Self { b, h })
    }
    pub fn field_strength(&self, flux_density_t: f64) -> Result<f64> {
        let b = flux_density_t;
        require(b >= 0.0 && b.is_finite(), "invalid B-H query")?;
        let n = self.b.len();
        if b >= self.b[n - 1] {
            let slope = (self.h[n - 1] - self.h[n - 2]) / (self.b[n - 1] - self.b[n - 2]);
            return Ok(self.h[n - 1] + slope * (b - self.b[n - 1]));
        }
        let hi = self.b.partition_point(|v| *v <= b);
        let lo = hi - 1;
        let t = (b - self.b[lo]) / (self.b[hi] - self.b[lo]);
        Ok(self.h[lo] + t * (self.h[hi] - self.h[lo]))
    }
    /// `nu = H/B`, with the initial slope near zero field.
    pub fn reluctivity(&self, flux_density_t: f64) -> Result<f64> {
        if flux_density_t <= 1e-15 {
            return Ok(self.h[1] / self.b[1]);
        }
        Ok(self.field_strength(flux_density_t)? / flux_density_t)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularStrandedCoil2D {
    pub xmin: f64,
    pub xmax: f64,
    pub ymin: f64,
    pub ymax: f64,
    pub turns: f64,
    pub current_a: f64,
}

impl Default for RectangularStrandedCoil2D {
    fn default() -> Self {
        Self {
            xmin: 0.0,
            xmax: 0.0,
            ymin: 0.0,
            ymax: 0.0,
            turns: 1.0,
            current_a: 1.0,
        }
    }
}

impl RectangularStrandedCoil2D {
    pub fn area_m2(&self) -> Result<f64> {
        require(
            self.xmax > self.xmin
                && self.ymax > self.ymin
                && self.turns.is_finite()
                && self.current_a.is_finite(),
            "invalid rectangular stranded coil",
        )?;
        Ok((self.xmax - self.xmin) * (self.ymax - self.ymin))
    }
    pub fn contains(&self, p: Node2) -> bool {
        p.x >= self.xmin && p.x <= self.xmax && p.y >= self.ymin && p.y <= self.ymax
    }
    pub fn current_density_a_per_m2(&self, p: Node2) -> Result<f64> {
        let area = self.area_m2()?;
        Ok(if self.contains(p) {
            self.turns * self.current_a / area
        } else {
            0.0
        })
    }
}

/// `N/S * integral of A over the coil` for elements whose centroid lies in the coil.
pub fn coil_flux_linkage_wb_turn(
    mesh: &Mesh2D,
    vector_potential_a_per_m: &[f64],
    coil: &RectangularStrandedCoil2D,
) -> Result<f64> {
    mesh.validate()?;
    let a = vector_potential_a_per_m;
    require(
        a.len() == mesh.node_count(),
        "coil flux-linkage field size mismatch",
    )?;
    let area = coil.area_m2()?;
    let mut integral = 0.0;
    for t in &mesh.triangles {
        if !coil.contains(centroid(mesh, t)) {
            continue;
        }
        let [p0, p1, p2] = t.node.map(|i| mesh.nodes[i]);
        let twice_area = ((p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)).abs();
        integral += (a[t.node[0]] + a[t.node[1]] + a[t.node[2]]) / 3.0 * 0.5 * twice_area;
    }
    Ok(coil.turns * integral / area)
}

pub fn coil_inductance_h(
    mesh: &Mesh2D,
    vector_potential_a_per_m: &[f64],
    coil: &RectangularStrandedCoil2D,
) -> Result<f64> {
    require(
        coil.current_a.abs() > 1e-30,
        "coil inductance requires non-zero current",
    )?;
    Ok(coil_flux_linkage_wb_turn(mesh, vector_potential_a_per_m, coil)? / coil.current_a)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MaxwellForceTorque2D {
    pub force_x_n_per_m: f64,
    pub force_y_n_per_m: f64,
    /// Per unit out-of-plane depth.
    pub torque_z_n: f64,
}

/// Maxwell-stress traction integrated over one boundary patch with element-constant `B` and
/// uniform permeability; the normal points out of the adjacent triangle.
pub fn maxwell_stress_boundary_force_2d(
    mesh: &Mesh2D,
    element_flux_density: &[MagneticFluxDensity2],
    patch: i32,
    permeability_h_per_m: f64,
    torque_origin: Node2,
) -> Result<MaxwellForceTorque2D> {
    mesh.validate()?;
    let mu = permeability_h_per_m;
    require(
        element_flux_density.len() == mesh.element_count() && mu > 0.0,
        "invalid Maxwell-stress inputs",
    )?;
    let mut r = MaxwellForceTorque2D::default();
    for edge in mesh.boundary_edges.iter().filter(|e| e.patch == patch) {
        let owner = mesh
            .triangles
            .iter()
            .position(|t| t.node.contains(&edge.node[0]) && t.node.contains(&edge.node[1]))
            .ok_or_else(|| {
                SolverError::InvalidInput("boundary edge has no adjacent triangle".into())
            })?;
        let (p0, p1) = (mesh.nodes[edge.node[0]], mesh.nodes[edge.node[1]]);
        let (dx, dy) = (p1.x - p0.x, p1.y - p0.y);
        let length = dx.hypot(dy);
        let (mut nx, mut ny) = (dy / length, -dx / length);
        let c = centroid(mesh, &mesh.triangles[owner]);
        let mid = Node2 {
            x: 0.5 * (p0.x + p1.x),
            y: 0.5 * (p0.y + p1.y),
        };
        if nx * (c.x - mid.x) + ny * (c.y - mid.y) > 0.0 {
            nx = -nx;
            ny = -ny;
        }
        let b = element_flux_density[owner];
        let bn = b.x * nx + b.y * ny;
        let b2 = b.x * b.x + b.y * b.y;
        let fx = (b.x * bn - 0.5 * b2 * nx) / mu * length;
        let fy = (b.y * bn - 0.5 * b2 * ny) / mu * length;
        r.force_x_n_per_m += fx;
        r.force_y_n_per_m += fy;
        r.torque_z_n += (mid.x - torque_origin.x) * fy - (mid.y - torque_origin.y) * fx;
    }
    Ok(r)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearPoisson2DConfig {
    pub diffusivity: f64,
    pub cubic_coefficient: f64,
    pub newton: NewtonConfig,
}

impl Default for NonlinearPoisson2DConfig {
    fn default() -> Self {
        Self {
            diffusivity: 1.0,
            cubic_coefficient: 1.0,
            newton: NewtonConfig::default(),
        }
    }
}

/// P1 Newton solve of `-div(k grad u) + beta u^3 = f` with Dirichlet data on all boundary nodes.
#[derive(Clone, Debug)]
pub struct NonlinearPoisson2D {
    mesh: Mesh2D,
    config: NonlinearPoisson2DConfig,
    solution: Vec<f64>,
    result: Option<NewtonResult>,
}

impl NonlinearPoisson2D {
    pub fn new(mesh: Mesh2D, config: NonlinearPoisson2DConfig) -> Result<Self> {
        mesh.validate()?;
        require(
            config.diffusivity > 0.0 && config.cubic_coefficient >= 0.0,
            "invalid nonlinear Poisson coefficients",
        )?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            config,
            solution: vec![0.0; n],
            result: None,
        })
    }
    pub fn solve(
        &mut self,
        source: impl Fn(Node2) -> f64,
        dirichlet_value: impl Fn(Node2) -> f64,
    ) -> Result<()> {
        let n = self.mesh.node_count();
        let mut free_index = vec![usize::MAX; n];
        let mut free_count = 0;
        for i in 0..n {
            if self.mesh.boundary_node[i] {
                self.solution[i] = dirichlet_value(self.mesh.nodes[i]);
            } else {
                free_index[i] = free_count;
                free_count += 1;
            }
        }
        require(free_count > 0, "nonlinear Poisson mesh has no free nodes")?;
        let mut x = vec![0.0; free_count];
        for i in 0..n {
            if free_index[i] != usize::MAX {
                x[free_index[i]] = self.solution[i];
            }
        }
        let rule = gaussian_quadrature(ElementType::Tri3, 2)?;
        let (mesh, base, c) = (&self.mesh, &self.solution, self.config);
        let gather = |free: &[f64]| {
            let mut full = base.clone();
            for i in 0..n {
                if free_index[i] != usize::MAX {
                    full[i] = free[free_index[i]];
                }
            }
            full
        };
        let residual = |free: &[f64], out: &mut [f64]| -> Result<()> {
            require(out.len() == free_count, "nonlinear residual size mismatch")?;
            out.fill(0.0);
            let full = gather(free);
            for t in &mesh.triangles {
                let nodes = t.node.map(|i| mesh.nodes[i].point3());
                let mut local = [0.0; 3];
                for qp in &rule {
                    let shape = evaluate_shape(ElementType::Tri3, qp.point);
                    let iso = evaluate_isoparametric(ElementType::Tri3, &nodes, qp.point)?;
                    let (mut uq, mut gx, mut gy) = (0.0, 0.0, 0.0);
                    for j in 0..3 {
                        let u = full[t.node[j]];
                        uq += shape.value[j] * u;
                        gx += iso.gradient_physical[j].x * u;
                        gy += iso.gradient_physical[j].y * u;
                    }
                    let f = source(Node2 {
                        x: iso.point.x,
                        y: iso.point.y,
                    });
                    let w = qp.weight * iso.measure;
                    for i in 0..3 {
                        let gi = iso.gradient_physical[i];
                        local[i] += (c.diffusivity * (gi.x * gx + gi.y * gy)
                            + c.cubic_coefficient * shape.value[i] * uq * uq * uq
                            - shape.value[i] * f)
                            * w;
                    }
                }
                for i in 0..3 {
                    if free_index[t.node[i]] != usize::MAX {
                        out[free_index[t.node[i]]] += local[i];
                    }
                }
            }
            Ok(())
        };
        let jacobian = |free: &[f64]| {
            let full = gather(free);
            let mut b = CsrBuilder::new(free_count, free_count);
            for t in &mesh.triangles {
                let nodes = t.node.map(|i| mesh.nodes[i].point3());
                let mut local = [[0.0; 3]; 3];
                for qp in &rule {
                    let shape = evaluate_shape(ElementType::Tri3, qp.point);
                    let iso = evaluate_isoparametric(ElementType::Tri3, &nodes, qp.point)?;
                    let uq: f64 = (0..3).map(|j| shape.value[j] * full[t.node[j]]).sum();
                    let w = qp.weight * iso.measure;
                    for i in 0..3 {
                        for j in 0..3 {
                            let (gi, gj) = (iso.gradient_physical[i], iso.gradient_physical[j]);
                            local[i][j] += (c.diffusivity * (gi.x * gj.x + gi.y * gj.y)
                                + 3.0
                                    * c.cubic_coefficient
                                    * shape.value[i]
                                    * shape.value[j]
                                    * uq
                                    * uq)
                                * w;
                        }
                    }
                }
                for i in 0..3 {
                    let fi = free_index[t.node[i]];
                    if fi == usize::MAX {
                        continue;
                    }
                    for j in 0..3 {
                        let fj = free_index[t.node[j]];
                        if fj != usize::MAX {
                            b.add(fi, fj, local[i][j])?;
                        }
                    }
                }
            }
            b.build(0.0)
        };
        let result = newton_solve(&mut x, residual, jacobian, c.newton)?;
        require(
            result.converged,
            "NonlinearPoisson2D Newton solve did not converge",
        )?;
        self.result = Some(result);
        for i in 0..n {
            if free_index[i] != usize::MAX {
                self.solution[i] = x[free_index[i]];
            }
        }
        Ok(())
    }
    pub fn mesh(&self) -> &Mesh2D {
        &self.mesh
    }
    pub fn solution(&self) -> &[f64] {
        &self.solution
    }
    pub fn nonlinear_result(&self) -> Option<&NewtonResult> {
        self.result.as_ref()
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NonlinearTrussState {
    pub stretch: f64,
    pub green_lagrange_strain: f64,
    pub second_piola_stress: f64,
    pub axial_force: f64,
    pub tangent: f64,
}

/// St. Venant-Kirchhoff truss: `E = (lambda^2 - 1)/2`, `S = E_young E`, `N = A S lambda`.
pub fn green_lagrange_truss(
    reference_length: f64,
    current_length: f64,
    area: f64,
    youngs_modulus: f64,
) -> Result<NonlinearTrussState> {
    let (l0, l, a, e) = (reference_length, current_length, area, youngs_modulus);
    require(
        l0 > 0.0 && l > 0.0 && a > 0.0 && e > 0.0,
        "invalid nonlinear truss state",
    )?;
    let lambda = l / l0;
    let eps = 0.5 * (lambda * lambda - 1.0);
    let s = e * eps;
    Ok(NonlinearTrussState {
        stretch: lambda,
        green_lagrange_strain: eps,
        second_piola_stress: s,
        axial_force: a * s * lambda,
        tangent: a * e * (1.5 * lambda * lambda - 0.5) / l0,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NonlinearTriangleState {
    pub reference_area: f64,
    pub determinant: f64,
    /// Row-major F.
    pub deformation_gradient: [f64; 4],
    /// `[E_xx, E_yy, E_xy]`.
    pub green_lagrange_strain: [f64; 3],
    /// Nodal `x/y` internal-force residual.
    pub internal_force: [f64; 6],
}

/// Total-Lagrangian St. Venant-Kirchhoff CST: deformation gradient, Green-Lagrange strain and
/// first-Piola internal forces.
pub fn total_lagrangian_triangle(
    mesh: &Mesh2D,
    triangle: &Tri3,
    current_nodes: &[Node2; 3],
    youngs_modulus: f64,
    poisson_ratio: f64,
) -> Result<NonlinearTriangleState> {
    let (e, nu) = (youngs_modulus, poisson_ratio);
    require(
        e > 0.0 && nu > -1.0 && nu < 0.5,
        "invalid nonlinear triangle material",
    )?;
    require(
        triangle.node.iter().all(|&i| i < mesh.nodes.len()),
        "nonlinear triangle node index out of range",
    )?;
    let [a, b, c] = triangle.node.map(|i| mesh.nodes[i]);
    let (j11, j12, j21, j22) = (b.x - a.x, c.x - a.x, b.y - a.y, c.y - a.y);
    let det = j11 * j22 - j12 * j21;
    require(
        det > 0.0,
        "nonlinear triangle requires positive reference orientation",
    )?;
    let area = 0.5 * det;
    let inv = 1.0 / det;
    let grad = [
        [(b.y - c.y) * inv, (c.x - b.x) * inv],
        [(c.y - a.y) * inv, (a.x - c.x) * inv],
        [(a.y - b.y) * inv, (b.x - a.x) * inv],
    ];
    let cur = current_nodes;
    let (d11, d12, d21, d22) = (
        cur[1].x - cur[0].x,
        cur[2].x - cur[0].x,
        cur[1].y - cur[0].y,
        cur[2].y - cur[0].y,
    );
    let (i11, i12, i21, i22) = (j22 * inv, -j12 * inv, -j21 * inv, j11 * inv);
    let (f11, f12, f21, f22) = (
        d11 * i11 + d12 * i21,
        d11 * i12 + d12 * i22,
        d21 * i11 + d22 * i21,
        d21 * i12 + d22 * i22,
    );
    let exx = 0.5 * (f11 * f11 + f21 * f21 - 1.0);
    let eyy = 0.5 * (f12 * f12 + f22 * f22 - 1.0);
    let exy = 0.5 * (f11 * f12 + f21 * f22);
    let mu = e / (2.0 * (1.0 + nu));
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let sxx = 2.0 * mu * exx + lambda * (exx + eyy);
    let syy = 2.0 * mu * eyy + lambda * (exx + eyy);
    let sxy = 2.0 * mu * exy;
    let (p11, p12, p21, p22) = (
        f11 * sxx + f12 * sxy,
        f11 * sxy + f12 * syy,
        f21 * sxx + f22 * sxy,
        f21 * sxy + f22 * syy,
    );
    let mut r = NonlinearTriangleState {
        reference_area: area,
        determinant: f11 * f22 - f12 * f21,
        deformation_gradient: [f11, f12, f21, f22],
        green_lagrange_strain: [exx, eyy, exy],
        internal_force: [0.0; 6],
    };
    for i in 0..3 {
        r.internal_force[2 * i] = area * (p11 * grad[i][0] + p12 * grad[i][1]);
        r.internal_force[2 * i + 1] = area * (p21 * grad[i][0] + p22 * grad[i][1]);
    }
    Ok(r)
}

/// Linear-mushy-zone phase change with temperature-referenced enthalpy.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PhaseChangeMaterial {
    pub solidus: f64,
    pub liquidus: f64,
    pub heat_capacity_solid: f64,
    pub heat_capacity_liquid: f64,
    pub latent_heat: f64,
}

impl PhaseChangeMaterial {
    fn check(&self) -> Result<()> {
        require(
            self.liquidus > self.solidus
                && self.heat_capacity_solid > 0.0
                && self.heat_capacity_liquid > 0.0
                && self.latent_heat >= 0.0,
            "invalid phase-change material",
        )
    }
    pub fn liquid_fraction(&self, temperature: f64) -> Result<f64> {
        self.check()?;
        Ok(((temperature - self.solidus) / (self.liquidus - self.solidus)).clamp(0.0, 1.0))
    }
    pub fn specific_enthalpy(&self, temperature: f64) -> Result<f64> {
        self.check()?;
        let t = temperature;
        let (cs, cl, ts) = (
            self.heat_capacity_solid,
            self.heat_capacity_liquid,
            self.solidus,
        );
        Ok(if t <= ts {
            cs * t
        } else if t >= self.liquidus {
            cs * ts + self.latent_heat + cl * (t - ts)
        } else {
            let f = (t - ts) / (self.liquidus - ts);
            cs * ts + f * self.latent_heat + 0.5 * (cs + cl) * (t - ts)
        })
    }
    pub fn effective_heat_capacity(&self, temperature: f64) -> Result<f64> {
        self.check()?;
        Ok(if temperature <= self.solidus {
            self.heat_capacity_solid
        } else if temperature >= self.liquidus {
            self.heat_capacity_liquid
        } else {
            0.5 * (self.heat_capacity_solid + self.heat_capacity_liquid)
                + self.latent_heat / (self.liquidus - self.solidus)
        })
    }
    /// 120-step bisection inverse of [`Self::specific_enthalpy`].
    pub fn temperature_from_enthalpy(&self, enthalpy: f64) -> Result<f64> {
        self.check()?;
        let mut lo = 0.0f64.min(self.solidus - 10000.0);
        let mut hi = self.liquidus
            + 10000.0f64.max(
                enthalpy.abs() / self.heat_capacity_solid.min(self.heat_capacity_liquid) + 1.0,
            );
        for _ in 0..120 {
            let mid = 0.5 * (lo + hi);
            if self.specific_enthalpy(mid)? < enthalpy {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(0.5 * (lo + hi))
    }
}

/// Parallel (arithmetic) mixture conductivity.
pub fn porous_effective_conductivity(
    porosity: f64,
    solid_conductivity: f64,
    fluid_conductivity: f64,
) -> Result<f64> {
    require(
        (0.0..=1.0).contains(&porosity) && solid_conductivity >= 0.0 && fluid_conductivity >= 0.0,
        "invalid porous conductivity",
    )?;
    Ok((1.0 - porosity) * solid_conductivity + porosity * fluid_conductivity)
}
