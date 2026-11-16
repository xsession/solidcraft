//! P1 Brinkman, penalty Stokes, mixed Stokes local block, Newmark dynamic bar, complex
//! eddy-current, Darcy, thermoelastic block and stabilization helpers, ported from upstream
//! `cfd/solvers/fem/{brinkman2d, stokes2d, mixed_stokes, dynamic_bar1d, eddy_current2d,
//! darcy2d, thermoelastic_block, transport_stabilization}` (MIT, cfd_solvers contributors).
use super::{
    evaluate_isoparametric, scalar_diffusion2d, Displacement2, ElementType, Mesh2D, Node2,
    ScalarBoundaryCondition2D, Tri3,
};
use crate::{
    core::{
        preconditioned_conjugate_gradient, restarted_gmres, CsrBuilder, CsrMatrix,
        Ilu0Preconditioner, IterativeSolverResult, JacobiPreconditioner, LinearConfig,
    },
    require, Result, SolverError,
};
use std::collections::BTreeMap;

/// Area, constant P1 gradients (`inv = 1/det` multiplied, as upstream) and centroid.
struct Geometry {
    area: f64,
    grad: [[f64; 2]; 3],
    centroid: Node2,
}

fn geometry(mesh: &Mesh2D, t: &Tri3, what: &str) -> Result<Geometry> {
    let [a, b, c] = t.node.map(|i| mesh.nodes[i]);
    let det = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    let area = 0.5 * det.abs();
    require(area > 0.0, format!("degenerate {what} triangle"))?;
    let inv = 1.0 / det;
    Ok(Geometry {
        area,
        grad: [
            [(b.y - c.y) * inv, (c.x - b.x) * inv],
            [(c.y - a.y) * inv, (a.x - c.x) * inv],
            [(a.y - b.y) * inv, (b.x - a.x) * inv],
        ],
        centroid: Node2 {
            x: (a.x + b.x + c.x) / 3.0,
            y: (a.y + b.y + c.y) / 3.0,
        },
    })
}

fn jacobi_pcg(
    a: &CsrMatrix,
    rhs: &[f64],
    x: &mut [f64],
    max_iterations: usize,
    tolerance: f64,
) -> Result<IterativeSolverResult> {
    let jacobi = JacobiPreconditioner::new(&a.diagonal(), 1e-300)?;
    preconditioned_conjugate_gradient(
        rhs,
        x,
        |v, o| a.multiply(v, o),
        &jacobi,
        LinearConfig {
            max_iterations,
            relative_tolerance: tolerance,
            ..Default::default()
        },
    )
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Brinkman2DConfig {
    pub viscosity: f64,
    pub effective_viscosity: f64,
    pub permeability: f64,
    pub max_iterations: usize,
    pub tolerance: f64,
}

impl Default for Brinkman2DConfig {
    fn default() -> Self {
        Self {
            viscosity: 1e-3,
            effective_viscosity: 1e-3,
            permeability: 1e-8,
            max_iterations: 3000,
            tolerance: 1e-10,
        }
    }
}

/// Brinkman `-mu_e lap u + (mu/K) u = f` with component-wise P1 solves and lumped-consistent drag mass.
#[derive(Clone, Debug)]
pub struct Brinkman2D {
    mesh: Mesh2D,
    config: Brinkman2DConfig,
    velocity: Vec<Displacement2>,
    zero_boundary: bool,
}

impl Brinkman2D {
    pub fn new(mesh: Mesh2D, config: Brinkman2DConfig) -> Result<Self> {
        mesh.validate()?;
        require(
            config.viscosity > 0.0 && config.effective_viscosity > 0.0 && config.permeability > 0.0,
            "invalid Brinkman controls",
        )?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            config,
            velocity: vec![Displacement2::default(); n],
            zero_boundary: true,
        })
    }
    /// `true` (default) clamps every boundary node to zero velocity.
    pub fn set_zero_velocity_boundary(&mut self, all_boundary: bool) {
        self.zero_boundary = all_boundary;
    }
    pub fn solve(&mut self, body_force: impl Fn(Node2) -> Displacement2) -> Result<()> {
        let n = self.mesh.node_count();
        let mut map = vec![usize::MAX; n];
        let mut nf = 0;
        for (i, m) in map.iter_mut().enumerate() {
            if !(self.zero_boundary && self.mesh.boundary_node[i]) {
                *m = nf;
                nf += 1;
            }
        }
        let mut b = CsrBuilder::new(nf, nf);
        let (mut rx, mut ry) = (vec![0.0; nf], vec![0.0; nf]);
        let c = self.config;
        let drag = c.viscosity / c.permeability;
        for t in &self.mesh.triangles {
            let g = geometry(&self.mesh, t, "Brinkman")?;
            let f = body_force(g.centroid);
            for i in 0..3 {
                let row = map[t.node[i]];
                if row == usize::MAX {
                    continue;
                }
                rx[row] += f.x * g.area / 3.0;
                ry[row] += f.y * g.area / 3.0;
                for j in 0..3 {
                    let col = map[t.node[j]];
                    if col == usize::MAX {
                        continue;
                    }
                    let k = c.effective_viscosity
                        * g.area
                        * (g.grad[i][0] * g.grad[j][0] + g.grad[i][1] * g.grad[j][1])
                        + drag * g.area * if i == j { 2.0 } else { 1.0 } / 12.0;
                    b.add(row, col, k)?;
                }
            }
        }
        let a = b.build(0.0)?;
        let (mut x, mut y) = (vec![0.0; nf], vec![0.0; nf]);
        for (r, u) in [(&rx, &mut x), (&ry, &mut y)] {
            let q = jacobi_pcg(&a, r, u, c.max_iterations, c.tolerance)?;
            require(q.converged, "Brinkman solve failed")?;
        }
        for i in 0..n {
            if map[i] != usize::MAX {
                self.velocity[i] = Displacement2 {
                    x: x[map[i]],
                    y: y[map[i]],
                };
            }
        }
        Ok(())
    }
    pub fn velocity(&self) -> &[Displacement2] {
        &self.velocity
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PenaltyStokes2DConfig {
    pub viscosity: f64,
    pub divergence_penalty: f64,
    pub max_iterations: usize,
    pub relative_tolerance: f64,
}

impl Default for PenaltyStokes2DConfig {
    fn default() -> Self {
        Self {
            viscosity: 1.0,
            divergence_penalty: 1000.0,
            max_iterations: 5000,
            relative_tolerance: 1e-10,
        }
    }
}

/// Velocity-only Stokes with a divergence penalty; boundary nodes take the Dirichlet velocity.
pub struct PenaltyStokes2D {
    mesh: Mesh2D,
    config: PenaltyStokes2DConfig,
    boundary: Option<Box<dyn Fn(Node2) -> Displacement2>>,
    velocity: Vec<Displacement2>,
}

impl PenaltyStokes2D {
    pub fn new(mesh: Mesh2D, config: PenaltyStokes2DConfig) -> Result<Self> {
        mesh.validate()?;
        require(
            config.viscosity > 0.0 && config.divergence_penalty >= 0.0,
            "invalid penalty Stokes controls",
        )?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            config,
            boundary: None,
            velocity: vec![Displacement2::default(); n],
        })
    }
    pub fn set_dirichlet(&mut self, velocity: impl Fn(Node2) -> Displacement2 + 'static) {
        self.boundary = Some(Box::new(velocity));
    }
    pub fn solve(&mut self, body_force: Option<&dyn Fn(Node2) -> Displacement2>) -> Result<()> {
        let boundary = self.boundary.as_ref().ok_or_else(|| {
            SolverError::InvalidInput("PenaltyStokes2D requires boundary velocity".into())
        })?;
        let n = self.mesh.node_count();
        let mut map = vec![usize::MAX; 2 * n];
        let mut nf = 0;
        for i in 0..n {
            if self.mesh.boundary_node[i] {
                self.velocity[i] = boundary(self.mesh.nodes[i]);
            } else {
                map[2 * i] = nf;
                map[2 * i + 1] = nf + 1;
                nf += 2;
            }
        }
        let c = self.config;
        let mut b = CsrBuilder::new(nf, nf);
        let mut rhs = vec![0.0; nf];
        for t in &self.mesh.triangles {
            let g = geometry(&self.mesh, t, "Stokes")?;
            let f = body_force.map_or(Displacement2::default(), |f| f(g.centroid));
            for i in 0..3 {
                for a in 0..2 {
                    let row = map[2 * t.node[i] + a];
                    if row == usize::MAX {
                        continue;
                    }
                    rhs[row] += if a == 1 { f.y } else { f.x } * g.area / 3.0;
                    for j in 0..3 {
                        for cc in 0..2 {
                            let k = c.viscosity
                                * g.area
                                * (g.grad[i][0] * g.grad[j][0] + g.grad[i][1] * g.grad[j][1])
                                * if a == cc { 1.0 } else { 0.0 }
                                + c.divergence_penalty * g.area * g.grad[i][a] * g.grad[j][cc];
                            let col = map[2 * t.node[j] + cc];
                            if col == usize::MAX {
                                let v = self.velocity[t.node[j]];
                                rhs[row] -= k * if cc == 1 { v.y } else { v.x };
                            } else {
                                b.add(row, col, k)?;
                            }
                        }
                    }
                }
            }
        }
        let a = b.build(0.0)?;
        let mut x = vec![0.0; nf];
        let r = jacobi_pcg(&a, &rhs, &mut x, c.max_iterations, c.relative_tolerance)?;
        require(r.converged, "PenaltyStokes2D did not converge")?;
        for i in 0..n {
            if !self.mesh.boundary_node[i] {
                self.velocity[i] = Displacement2 {
                    x: x[map[2 * i]],
                    y: x[map[2 * i + 1]],
                };
            }
        }
        Ok(())
    }
    pub fn velocity(&self) -> &[Displacement2] {
        &self.velocity
    }
    /// Area-weighted RMS of the element divergence.
    pub fn divergence_l2(&self) -> Result<f64> {
        let (mut s, mut a) = (0.0, 0.0);
        for t in &self.mesh.triangles {
            let g = geometry(&self.mesh, t, "Stokes")?;
            let mut div = 0.0;
            for i in 0..3 {
                let v = self.velocity[t.node[i]];
                div += v.x * g.grad[i][0] + v.y * g.grad[i][1];
            }
            s += div * div * g.area;
            a += g.area;
        }
        Ok((s / a).sqrt())
    }
}

/// Local P1/P1 mixed Stokes block, ordered `[ux0, uy0, ux1, uy1, ux2, uy2, p0, p1, p2]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedStokesTri3LocalSystem {
    pub matrix: [f64; 81],
    pub area: f64,
}

/// Equal-order mixed Stokes element with pressure-mass stabilization (upstream default 1e-8).
pub fn assemble_mixed_stokes_tri3(
    mesh: &Mesh2D,
    triangle: &Tri3,
    viscosity: f64,
    pressure_stabilization: f64,
) -> Result<MixedStokesTri3LocalSystem> {
    require(
        viscosity > 0.0 && pressure_stabilization >= 0.0,
        "invalid mixed Stokes controls",
    )?;
    require(
        triangle.node.iter().all(|&i| i < mesh.nodes.len()),
        "mixed Stokes node index out of range",
    )?;
    let g = geometry(mesh, triangle, "mixed Stokes")?;
    let mut r = MixedStokesTri3LocalSystem {
        matrix: [0.0; 81],
        area: g.area,
    };
    let mut add = |row: usize, col: usize, v: f64| r.matrix[row * 9 + col] += v;
    for i in 0..3 {
        for j in 0..3 {
            let visc =
                viscosity * g.area * (g.grad[i][0] * g.grad[j][0] + g.grad[i][1] * g.grad[j][1]);
            add(2 * i, 2 * j, visc);
            add(2 * i + 1, 2 * j + 1, visc);
            for component in 0..2 {
                let coupling = -g.area / 3.0 * g.grad[j][component];
                add(6 + i, 2 * j + component, coupling);
                add(2 * j + component, 6 + i, coupling);
            }
            add(
                6 + i,
                6 + j,
                pressure_stabilization * g.area * if i == j { 2.0 } else { 1.0 } / 12.0,
            );
        }
    }
    Ok(r)
}

/// Local `[u_x, u_y, T]` thermoelastic block (row-major 3x3).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermoelasticBlockLocalSystem {
    pub matrix: [f64; 9],
}

pub fn assemble_thermoelastic_block(
    area: f64,
    mechanical_stiffness: f64,
    thermal_conductivity: f64,
    thermal_capacity: f64,
    time_step: f64,
    thermal_expansion: f64,
) -> Result<ThermoelasticBlockLocalSystem> {
    require(
        area > 0.0
            && mechanical_stiffness > 0.0
            && thermal_conductivity >= 0.0
            && thermal_capacity >= 0.0
            && time_step > 0.0,
        "invalid thermoelastic block controls",
    )?;
    let mut m = [0.0; 9];
    m[0] = mechanical_stiffness * area;
    m[4] = mechanical_stiffness * area;
    m[8] = (thermal_conductivity + thermal_capacity / time_step) * area;
    m[2] = -thermal_expansion * mechanical_stiffness * area;
    m[5] = m[2];
    m[6] = m[2];
    m[7] = m[2];
    Ok(ThermoelasticBlockLocalSystem { matrix: m })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SupgParameters {
    pub tau_s: f64,
    pub element_peclet: f64,
}

/// `tau = (( 2/dt)^2 + (2|u|/h)^2 + (4D/h^2)^2)^(-1/2)` and `Pe = |u| h / 2D`.
pub fn supg_parameters(
    element_length_m: f64,
    advection_speed_m_per_s: f64,
    diffusivity_m2_per_s: f64,
    timestep_s: f64,
) -> Result<SupgParameters> {
    let (h, u, d, dt) = (
        element_length_m,
        advection_speed_m_per_s,
        diffusivity_m2_per_s,
        timestep_s,
    );
    require(
        h > 0.0
            && d >= 0.0
            && dt > 0.0
            && h.is_finite()
            && u.is_finite()
            && d.is_finite()
            && dt.is_finite(),
        "invalid SUPG controls",
    )?;
    let advective = 2.0 * u.abs() / h;
    let diffusive = 4.0 * d / (h * h);
    let temporal = 2.0 / dt;
    let inverse_tau = (temporal * temporal + advective * advective + diffusive * diffusive).sqrt();
    Ok(SupgParameters {
        tau_s: 1.0 / inverse_tau,
        element_peclet: u.abs() * h / (2.0 * d).max(1e-30),
    })
}

/// Upwind plus symmetric-penalty face flux.
pub fn dg_upwind_flux(
    left_value: f64,
    right_value: f64,
    outward_speed: f64,
    penalty: f64,
) -> Result<f64> {
    require(
        left_value.is_finite()
            && right_value.is_finite()
            && outward_speed.is_finite()
            && penalty >= 0.0
            && penalty.is_finite(),
        "invalid DG flux controls",
    )?;
    let upwind = if outward_speed >= 0.0 {
        left_value
    } else {
        right_value
    };
    Ok(outward_speed * upwind + penalty * (left_value - right_value))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicBar1DConfig {
    pub elements: usize,
    pub length: f64,
    pub area: f64,
    pub young_modulus: f64,
    pub density: f64,
    pub dt: f64,
    pub newmark_beta: f64,
    pub newmark_gamma: f64,
    pub linear_iterations: usize,
    pub linear_tolerance: f64,
}

impl Default for DynamicBar1DConfig {
    fn default() -> Self {
        Self {
            elements: 80,
            length: 1.0,
            area: 1.0,
            young_modulus: 200e9,
            density: 7800.0,
            dt: 1e-5,
            newmark_beta: 0.25,
            newmark_gamma: 0.5,
            linear_iterations: 2000,
            linear_tolerance: 1e-11,
        }
    }
}

/// Fixed-free axial bar (node 0 clamped) with consistent mass and implicit Newmark integration.
#[derive(Clone, Debug)]
pub struct DynamicBar1D {
    c: DynamicBar1DConfig,
    k: CsrMatrix,
    m: CsrMatrix,
    effective: CsrMatrix,
    u: Vec<f64>,
    v: Vec<f64>,
    a: Vec<f64>,
    time: f64,
}

impl DynamicBar1D {
    pub fn new(c: DynamicBar1DConfig) -> Result<Self> {
        require(
            c.elements >= 2
                && c.length > 0.0
                && c.area > 0.0
                && c.young_modulus > 0.0
                && c.density > 0.0
                && c.dt > 0.0
                && c.newmark_beta > 0.0
                && c.newmark_gamma > 0.0,
            "invalid dynamic bar controls",
        )?;
        let n = c.elements;
        let h = c.length / n as f64;
        let ke = c.young_modulus * c.area / h;
        let me = c.density * c.area * h / 6.0;
        let (mut kb, mut mb, mut ab) = (
            CsrBuilder::new(n, n),
            CsrBuilder::new(n, n),
            CsrBuilder::new(n, n),
        );
        let kk = [[ke, -ke], [-ke, ke]];
        let mm = [[2.0 * me, me], [me, 2.0 * me]];
        for e in 0..n {
            let g = [e, e + 1];
            for i in 0..2 {
                if g[i] == 0 {
                    continue;
                }
                for j in 0..2 {
                    if g[j] == 0 {
                        continue;
                    }
                    let (r, co) = (g[i] - 1, g[j] - 1);
                    kb.add(r, co, kk[i][j])?;
                    mb.add(r, co, mm[i][j])?;
                    ab.add(r, co, mm[i][j] + c.newmark_beta * c.dt * c.dt * kk[i][j])?;
                }
            }
        }
        Ok(Self {
            c,
            k: kb.build(0.0)?,
            m: mb.build(0.0)?,
            effective: ab.build(0.0)?,
            u: vec![0.0; n],
            v: vec![0.0; n],
            a: vec![0.0; n],
            time: 0.0,
        })
    }
    /// Sets nodal displacement/velocity at `x = L (i+1)/N` and solves `M a = -K u`.
    pub fn initialize(
        &mut self,
        displacement: impl Fn(f64) -> f64,
        velocity: Option<&dyn Fn(f64) -> f64>,
    ) -> Result<()> {
        let n = self.c.elements;
        for i in 0..n {
            let x = self.c.length * (i + 1) as f64 / n as f64;
            self.u[i] = displacement(x);
            self.v[i] = velocity.map_or(0.0, |f| f(x));
        }
        let mut ku = vec![0.0; n];
        self.k.multiply(&self.u, &mut ku)?;
        let neg: Vec<f64> = ku.iter().map(|v| -v).collect();
        let mut a = vec![0.0; n];
        let r = jacobi_pcg(
            &self.m,
            &neg,
            &mut a,
            self.c.linear_iterations,
            self.c.linear_tolerance,
        )?;
        require(r.converged, "dynamic bar initial acceleration solve failed")?;
        self.a = a;
        self.time = 0.0;
        Ok(())
    }
    pub fn step(&mut self, tip_force: f64) -> Result<()> {
        let (dt, b, g) = (self.c.dt, self.c.newmark_beta, self.c.newmark_gamma);
        let n = self.u.len();
        let up: Vec<f64> = (0..n)
            .map(|i| self.u[i] + dt * self.v[i] + dt * dt * (0.5 - b) * self.a[i])
            .collect();
        let vp: Vec<f64> = (0..n)
            .map(|i| self.v[i] + dt * (1.0 - g) * self.a[i])
            .collect();
        let mut kup = vec![0.0; n];
        self.k.multiply(&up, &mut kup)?;
        let mut rhs: Vec<f64> = kup.iter().map(|v| -v).collect();
        rhs[n - 1] += tip_force;
        let mut candidate = self.a.clone();
        let r = jacobi_pcg(
            &self.effective,
            &rhs,
            &mut candidate,
            self.c.linear_iterations,
            self.c.linear_tolerance,
        )?;
        require(r.converged, "dynamic bar Newmark solve failed")?;
        self.a = candidate;
        for i in 0..n {
            self.u[i] = up[i] + b * dt * dt * self.a[i];
            self.v[i] = vp[i] + g * dt * self.a[i];
        }
        self.time += dt;
        Ok(())
    }
    pub fn run(&mut self, steps: usize, tip_force: f64) -> Result<()> {
        for _ in 0..steps {
            self.step(tip_force)?;
        }
        Ok(())
    }
    pub fn displacement(&self) -> &[f64] {
        &self.u
    }
    pub fn velocity(&self) -> &[f64] {
        &self.v
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    fn quadratic(m: &CsrMatrix, x: &[f64]) -> Result<f64> {
        let mut y = vec![0.0; x.len()];
        m.multiply(x, &mut y)?;
        Ok(x.iter().zip(&y).map(|(a, b)| 0.5 * a * b).sum())
    }
    pub fn kinetic_energy(&self) -> Result<f64> {
        Self::quadratic(&self.m, &self.v)
    }
    pub fn strain_energy(&self) -> Result<f64> {
        Self::quadratic(&self.k, &self.u)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EddyCurrent2DConfig {
    pub frequency_hz: f64,
    pub max_iterations: usize,
    pub gmres_restart: usize,
    pub relative_tolerance: f64,
}

impl Default for EddyCurrent2DConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 50.0,
            max_iterations: 4000,
            gmres_restart: 50,
            relative_tolerance: 1e-10,
        }
    }
}

/// Time-harmonic `-div(nu grad A) + j omega sigma A = J` in real block form with ILU(0) GMRES;
/// boundary nodes are `A = 0`.
#[derive(Clone, Debug)]
pub struct EddyCurrent2D {
    mesh: Mesh2D,
    config: EddyCurrent2DConfig,
    a: Vec<crate::Complex64>,
    result: IterativeSolverResult,
}

impl EddyCurrent2D {
    pub fn new(mesh: Mesh2D, config: EddyCurrent2DConfig) -> Result<Self> {
        mesh.validate()?;
        require(
            config.frequency_hz >= 0.0
                && config.max_iterations > 0
                && config.gmres_restart > 0
                && config.relative_tolerance > 0.0,
            "invalid eddy-current controls",
        )?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            config,
            a: vec![crate::Complex64::new(0.0, 0.0); n],
            result: IterativeSolverResult::default(),
        })
    }
    pub fn solve(
        &mut self,
        reluctivity: impl Fn(Node2) -> f64,
        conductivity: impl Fn(Node2) -> f64,
        source_current: impl Fn(Node2) -> crate::Complex64,
    ) -> Result<()> {
        let n = self.mesh.node_count();
        let mut map = vec![usize::MAX; n];
        let mut nf = 0;
        for (i, m) in map.iter_mut().enumerate() {
            if !self.mesh.boundary_node[i] {
                *m = nf;
                nf += 1;
            }
        }
        require(nf > 0, "no eddy-current free nodes")?;
        let mut b = CsrBuilder::new(2 * nf, 2 * nf);
        let mut rhs = vec![0.0; 2 * nf];
        let omega = 2.0 * std::f64::consts::PI * self.config.frequency_hz;
        for t in &self.mesh.triangles {
            let g = geometry(&self.mesh, t, "eddy-current")?;
            let (rel, cond, src) = (
                reluctivity(g.centroid),
                conductivity(g.centroid),
                source_current(g.centroid),
            );
            require(rel > 0.0 && cond >= 0.0, "invalid eddy-current material")?;
            for i in 0..3 {
                let ri = map[t.node[i]];
                if ri == usize::MAX {
                    continue;
                }
                rhs[ri] += src.re * g.area / 3.0;
                rhs[nf + ri] += src.im * g.area / 3.0;
                for j in 0..3 {
                    let cj = map[t.node[j]];
                    if cj == usize::MAX {
                        continue;
                    }
                    let k =
                        rel * g.area * (g.grad[i][0] * g.grad[j][0] + g.grad[i][1] * g.grad[j][1]);
                    let mass = cond * g.area * if i == j { 2.0 } else { 1.0 } / 12.0;
                    b.add(ri, cj, k)?;
                    b.add(ri, nf + cj, -omega * mass)?;
                    b.add(nf + ri, cj, omega * mass)?;
                    b.add(nf + ri, nf + cj, k)?;
                }
            }
        }
        let a = b.build(0.0)?;
        let ilu = Ilu0Preconditioner::new(&a, 1e-300)?;
        let mut x = vec![0.0; 2 * nf];
        let c = self.config;
        self.result = restarted_gmres(
            &rhs,
            &mut x,
            |v, o| a.multiply(v, o),
            &ilu,
            LinearConfig {
                max_iterations: c.max_iterations,
                relative_tolerance: c.relative_tolerance,
                restart: c.gmres_restart,
            },
        )?;
        require(self.result.converged, "eddy-current GMRES failed")?;
        for (i, m) in map.iter().enumerate() {
            self.a[i] = if *m == usize::MAX {
                crate::Complex64::new(0.0, 0.0)
            } else {
                crate::Complex64::new(x[*m], x[nf + m])
            };
        }
        Ok(())
    }
    pub fn vector_potential(&self) -> &[crate::Complex64] {
        &self.a
    }
    pub fn linear_result(&self) -> &IterativeSolverResult {
        &self.result
    }
    pub fn mesh(&self) -> &Mesh2D {
        &self.mesh
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DarcyVelocity2 {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Darcy2DConfig {
    pub viscosity: f64,
    pub density: f64,
    pub gravity: DarcyVelocity2,
    pub linear: LinearConfig,
}

impl Default for Darcy2DConfig {
    fn default() -> Self {
        Self {
            viscosity: 1e-3,
            density: 1000.0,
            gravity: DarcyVelocity2::default(),
            linear: LinearConfig::default(),
        }
    }
}

/// Saturated Darcy pressure `-div((K/mu) grad p) = q` with velocity `-(K/mu)(grad p - rho g)`.
pub struct Darcy2D {
    mesh: Mesh2D,
    config: Darcy2DConfig,
    boundary: BTreeMap<i32, ScalarBoundaryCondition2D<'static>>,
    pressure: Vec<f64>,
    linear: IterativeSolverResult,
}

impl Darcy2D {
    pub fn new(mesh: Mesh2D, config: Darcy2DConfig) -> Result<Self> {
        mesh.validate()?;
        require(
            config.viscosity > 0.0 && config.density >= 0.0,
            "invalid Darcy properties",
        )?;
        let n = mesh.node_count();
        Ok(Self {
            mesh,
            config,
            boundary: BTreeMap::new(),
            pressure: vec![0.0; n],
            linear: IterativeSolverResult::default(),
        })
    }
    pub fn set_boundary(&mut self, patch: i32, condition: ScalarBoundaryCondition2D<'static>) {
        self.boundary.insert(patch, condition);
    }
    pub fn solve(
        &mut self,
        permeability: impl Fn(Node2) -> f64,
        volumetric_source: Option<&dyn Fn(Node2) -> f64>,
    ) -> Result<()> {
        let mu = self.config.viscosity;
        let s = scalar_diffusion2d(
            &self.mesh,
            |p| permeability(p) / mu,
            |_| 0.0,
            |p| volumetric_source.map_or(0.0, |f| f(p)),
            &self.boundary,
            self.config.linear,
        )?;
        self.pressure = s.values;
        self.linear = s.linear;
        Ok(())
    }
    pub fn mesh(&self) -> &Mesh2D {
        &self.mesh
    }
    pub fn pressure(&self) -> &[f64] {
        &self.pressure
    }
    pub fn linear_result(&self) -> &IterativeSolverResult {
        &self.linear
    }
    /// Element (centroid) Darcy velocity.
    pub fn element_velocity(
        &self,
        permeability: impl Fn(Node2) -> f64,
    ) -> Result<Vec<DarcyVelocity2>> {
        let mut out = Vec::with_capacity(self.mesh.triangles.len());
        for t in &self.mesh.triangles {
            let nodes = t.node.map(|i| self.mesh.nodes[i].point3());
            let iso = evaluate_isoparametric(
                ElementType::Tri3,
                &nodes,
                super::ReferencePoint {
                    xi: 1.0 / 3.0,
                    eta: 1.0 / 3.0,
                    zeta: 0.0,
                },
            )?;
            let (mut gx, mut gy) = (0.0, 0.0);
            for i in 0..3 {
                gx += self.pressure[t.node[i]] * iso.gradient_physical[i].x;
                gy += self.pressure[t.node[i]] * iso.gradient_physical[i].y;
            }
            let [a, b, c] = t.node.map(|i| self.mesh.nodes[i]);
            let k = permeability(Node2 {
                x: (a.x + b.x + c.x) / 3.0,
                y: (a.y + b.y + c.y) / 3.0,
            });
            require(k > 0.0, "non-positive Darcy permeability")?;
            let mobility = k / self.config.viscosity;
            let (rho, g) = (self.config.density, self.config.gravity);
            out.push(DarcyVelocity2 {
                x: -mobility * (gx - rho * g.x),
                y: -mobility * (gy - rho * g.y),
            });
        }
        Ok(out)
    }
}
