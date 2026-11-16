//! Implicit Newmark/generalized-alpha/HHT integrator for generalized coordinates, a fixed-base
//! reduced-coordinate articulated tree and a floating-frame flexible-body interface.
#![allow(clippy::needless_range_loop)]
use super::{
    component_mul, cross, dot, inverse_rotate, norm, norm_squared, normalized_or_zero,
    normalized_quaternion, quaternion_from_axis_angle, rotate, Quaternion, RigidBody, Vec3,
};
use crate::{require, Result, SolverError};
use std::sync::Arc;

fn solve_dense(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Result<Vec<f64>> {
    require(
        a.len() == n * n && b.len() == n,
        "dense solve size mismatch",
    )?;
    for col in 0..n {
        let mut pivot = col;
        for row in col + 1..n {
            if a[row * n + col].abs() > a[pivot * n + col].abs() {
                pivot = row;
            }
        }
        let p = a[pivot * n + col].abs();
        if p.is_nan() || p <= 1e-14 {
            return Err(SolverError::Singular(
                "singular generalized-coordinate matrix".into(),
            ));
        }
        if pivot != col {
            for j in col..n {
                a.swap(pivot * n + j, col * n + j);
            }
            b.swap(pivot, col);
        }
        let d = a[col * n + col];
        for j in col..n {
            a[col * n + j] /= d;
        }
        b[col] /= d;
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = a[row * n + col];
            if f == 0.0 {
                continue;
            }
            for j in col..n {
                let v = a[col * n + j];
                a[row * n + j] -= f * v;
            }
            let v = b[col];
            b[row] -= f * v;
        }
    }
    Ok(b)
}

fn inf_norm(x: &[f64]) -> f64 {
    x.iter().fold(0.0, |m, v| m.max(v.abs()))
}

fn world_inertia_mul(b: &RigidBody, v: Vec3) -> Result<Vec3> {
    let l = inverse_rotate(b.state.orientation, v)?;
    rotate(b.state.orientation, component_mul(b.inertia_diagonal, l))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitSecondOrderConfig {
    pub dt: f64,
    pub beta: f64,
    pub gamma: f64,
    /// `x_{n+1-af} = (1-af) x_{n+1} + af x_n`, `a_{n+1-am} = (1-am) a_{n+1} + am a_n`.
    pub alpha_m: f64,
    pub alpha_f: f64,
    pub newton_iterations: usize,
    pub nonlinear_tolerance: f64,
    pub finite_difference_epsilon: f64,
}

impl Default for ImplicitSecondOrderConfig {
    fn default() -> Self {
        Self {
            dt: 1e-3,
            beta: 0.25,
            gamma: 0.5,
            alpha_m: 0.0,
            alpha_f: 0.0,
            newton_iterations: 16,
            nonlinear_tolerance: 1e-10,
            finite_difference_epsilon: 1e-7,
        }
    }
}

impl ImplicitSecondOrderConfig {
    /// Upstream defaults `beta = 0.25`, `gamma = 0.5`.
    pub fn newmark(dt: f64, beta: f64, gamma: f64) -> Self {
        Self {
            dt,
            beta,
            gamma,
            alpha_m: 0.0,
            alpha_f: 0.0,
            ..Default::default()
        }
    }
    pub fn generalized_alpha(dt: f64, spectral_radius_infinity: f64) -> Result<Self> {
        let rho = spectral_radius_infinity;
        require(
            (0.0..=1.0).contains(&rho),
            "generalized-alpha spectral radius must be in [0,1]",
        )?;
        let am = (2.0 * rho - 1.0) / (rho + 1.0);
        let af = rho / (rho + 1.0);
        let term = 1.0 - am + af;
        Ok(Self {
            dt,
            alpha_m: am,
            alpha_f: af,
            gamma: 0.5 - am + af,
            beta: 0.25 * term * term,
            ..Default::default()
        })
    }
    pub fn hht(dt: f64, alpha: f64) -> Result<Self> {
        require(
            (-1.0 / 3.0..=0.0).contains(&alpha),
            "HHT alpha must be in [-1/3,0]",
        )?;
        Ok(Self {
            dt,
            alpha_m: 0.0,
            alpha_f: -alpha,
            gamma: 0.5 - alpha,
            beta: 0.25 * (1.0 - alpha) * (1.0 - alpha),
            ..Default::default()
        })
    }
}

/// `force(q, v, t, out)`.
pub type GeneralizedForceFunction =
    Arc<dyn Fn(&[f64], &[f64], f64, &mut [f64]) -> Result<()> + Send + Sync>;

/// Dense implicit second-order integrator with a constant mass matrix and finite-difference
/// Newton on the accelerations.
#[derive(Clone)]
pub struct ImplicitGeneralizedIntegrator {
    mass: Vec<f64>,
    n: usize,
    force: GeneralizedForceFunction,
    config: ImplicitSecondOrderConfig,
    q: Vec<f64>,
    v: Vec<f64>,
    a: Vec<f64>,
    time: f64,
    last_newton_iterations: usize,
    last_residual_norm: f64,
}

impl ImplicitGeneralizedIntegrator {
    pub fn new(
        mass_matrix: Vec<f64>,
        dofs: usize,
        force: GeneralizedForceFunction,
        config: ImplicitSecondOrderConfig,
    ) -> Result<Self> {
        let c = &config;
        require(
            dofs > 0 && dofs <= 4096 && mass_matrix.len() == dofs * dofs,
            "implicit multibody mass matrix size mismatch",
        )?;
        require(
            c.dt > 0.0
                && c.beta > 0.0
                && c.gamma > 0.0
                && c.newton_iterations > 0
                && c.nonlinear_tolerance > 0.0
                && c.finite_difference_epsilon > 0.0,
            "invalid implicit multibody controls",
        )?;
        require(
            c.alpha_m < 1.0 && c.alpha_f < 1.0,
            "generalized-alpha interpolation must remain implicit",
        )?;
        Ok(Self {
            mass: mass_matrix,
            n: dofs,
            force,
            config,
            q: vec![0.0; dofs],
            v: vec![0.0; dofs],
            a: vec![0.0; dofs],
            time: 0.0,
            last_newton_iterations: 0,
            last_residual_norm: 0.0,
        })
    }
    /// `v` may be empty (zero velocities); the initial acceleration solves `M a = f`.
    pub fn initialize(&mut self, q: Vec<f64>, v: Vec<f64>, time: f64) -> Result<()> {
        require(q.len() == self.n, "implicit multibody q size mismatch")?;
        let v = if v.is_empty() { vec![0.0; self.n] } else { v };
        require(v.len() == self.n, "implicit multibody v size mismatch")?;
        require(time.is_finite(), "implicit multibody time must be finite")?;
        let mut f = vec![0.0; self.n];
        (self.force)(&q, &v, time, &mut f)?;
        self.a = solve_dense(self.mass.clone(), f, self.n)?;
        self.q = q;
        self.v = v;
        self.time = time;
        self.last_newton_iterations = 0;
        self.last_residual_norm = 0.0;
        Ok(())
    }
    fn predict(&self, i: usize, a1: f64, q0: &[f64], v0: &[f64], a0: &[f64]) -> (f64, f64) {
        let (dt, b, g) = (self.config.dt, self.config.beta, self.config.gamma);
        let qp = q0[i] + dt * v0[i] + dt * dt * (0.5 - b) * a0[i];
        let vp = v0[i] + dt * (1.0 - g) * a0[i];
        (qp + b * dt * dt * a1, vp + g * dt * a1)
    }
    fn residual(&self, a1: &[f64], q0: &[f64], v0: &[f64], a0: &[f64]) -> Result<Vec<f64>> {
        let n = self.n;
        let c = &self.config;
        let (mut qaf, mut vaf, mut aam) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        for i in 0..n {
            let (q1, v1) = self.predict(i, a1[i], q0, v0, a0);
            qaf[i] = (1.0 - c.alpha_f) * q1 + c.alpha_f * q0[i];
            vaf[i] = (1.0 - c.alpha_f) * v1 + c.alpha_f * v0[i];
            aam[i] = (1.0 - c.alpha_m) * a1[i] + c.alpha_m * a0[i];
        }
        let t = (1.0 - c.alpha_f) * (self.time + c.dt) + c.alpha_f * self.time;
        let mut f = vec![0.0; n];
        (self.force)(&qaf, &vaf, t, &mut f)?;
        Ok((0..n)
            .map(|i| (0..n).fold(0.0, |s, j| s + self.mass[i * n + j] * aam[j]) - f[i])
            .collect())
    }
    /// One step; on failure the state is unchanged.
    pub fn step(&mut self) -> Result<()> {
        let n = self.n;
        let c = self.config;
        let (q0, v0, a0) = (self.q.clone(), self.v.clone(), self.a.clone());
        let mut a1 = a0.clone();
        let (mut iters, mut res_norm) = (0usize, self.last_residual_norm);
        for it in 0..c.newton_iterations {
            let mut r = self.residual(&a1, &q0, &v0, &a0)?;
            res_norm = inf_norm(&r);
            if res_norm <= c.nonlinear_tolerance {
                iters = it;
                break;
            }
            let mut jac = vec![0.0; n * n];
            for col in 0..n {
                let mut p = a1.clone();
                let h = c.finite_difference_epsilon * a1[col].abs().max(1.0);
                p[col] += h;
                let rp = self.residual(&p, &q0, &v0, &a0)?;
                for row in 0..n {
                    jac[row * n + col] = (rp[row] - r[row]) / h;
                }
            }
            r.iter_mut().for_each(|v| *v = -*v);
            let delta = solve_dense(jac, r, n)?;
            for i in 0..n {
                a1[i] += delta[i];
            }
            iters = it + 1;
            if inf_norm(&delta) <= c.nonlinear_tolerance * inf_norm(&a1).max(1.0) {
                res_norm = inf_norm(&self.residual(&a1, &q0, &v0, &a0)?);
                break;
            }
        }
        self.last_newton_iterations = iters;
        self.last_residual_norm = res_norm;
        if res_norm.is_nan() || res_norm > 100.0 * c.nonlinear_tolerance {
            return Err(SolverError::NumericalFailure(
                "implicit multibody Newton solve did not converge".into(),
            ));
        }
        for i in 0..n {
            let (q1, v1) = self.predict(i, a1[i], &q0, &v0, &a0);
            self.q[i] = q1;
            self.v[i] = v1;
        }
        self.a = a1;
        self.time += c.dt;
        Ok(())
    }
    pub fn run(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }
    pub fn q(&self) -> &[f64] {
        &self.q
    }
    pub fn v(&self) -> &[f64] {
        &self.v
    }
    pub fn a(&self) -> &[f64] {
        &self.a
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn last_newton_iterations(&self) -> usize {
        self.last_newton_iterations
    }
    pub fn last_residual_norm(&self) -> f64 {
        self.last_residual_norm
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ReducedJointType {
    #[default]
    Revolute,
    Prismatic,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedLink {
    /// `None` is the fixed base.
    pub parent: Option<usize>,
    pub joint_type: ReducedJointType,
    /// Parent-local (world for a root link).
    pub parent_joint_offset: Vec3,
    /// Child-frame vector from the joint to the child COM.
    pub joint_to_com_local: Vec3,
    pub axis_parent: Vec3,
    pub reference_orientation: Quaternion,
    pub body: RigidBody,
    pub q: f64,
    pub qd: f64,
    pub qdd: f64,
    pub applied_effort: f64,
    pub damping: f64,
}

impl Default for ArticulatedLink {
    fn default() -> Self {
        Self {
            parent: None,
            joint_type: ReducedJointType::Revolute,
            parent_joint_offset: Vec3::ZERO,
            joint_to_com_local: Vec3::ZERO,
            axis_parent: Vec3::new(0.0, 0.0, 1.0),
            reference_orientation: Quaternion::default(),
            body: RigidBody::default(),
            q: 0.0,
            qd: 0.0,
            qdd: 0.0,
            applied_effort: 0.0,
            damping: 0.0,
        }
    }
}

/// Fixed-base tree: `M(q) = sum J^T M_b J` from body Jacobians, solved densely (O(n^3) baseline).
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedSystem {
    gravity: Vec3,
    links: Vec<ArticulatedLink>,
    origin: Vec<Vec3>,
    axis: Vec<Vec3>,
}

impl ArticulatedSystem {
    /// Upstream default gravity `(0, -9.81, 0)`.
    pub fn new(gravity: Vec3) -> Self {
        Self {
            gravity,
            links: vec![],
            origin: vec![],
            axis: vec![],
        }
    }
    pub fn add_link(&mut self, mut link: ArticulatedLink) -> Result<usize> {
        require(
            link.parent.is_none_or(|p| p < self.links.len()),
            "articulated parent must precede child",
        )?;
        require(
            norm(link.axis_parent) >= 1e-14,
            "articulated joint axis must be nonzero",
        )?;
        require(
            self.links.len() < 4096,
            "articulated tree too large for the dense baseline",
        )?;
        link.axis_parent = normalized_or_zero(link.axis_parent, 1e-15);
        link.body.validate()?;
        let mut t = self.clone();
        t.links.push(link);
        t.origin.push(Vec3::ZERO);
        t.axis.push(Vec3::ZERO);
        t.forward_kinematics()?;
        *self = t;
        Ok(self.links.len() - 1)
    }
    pub fn set_effort(&mut self, link: usize, effort: f64) -> Result<()> {
        self.links
            .get_mut(link)
            .ok_or_else(|| SolverError::InvalidInput("articulated link index".into()))?
            .applied_effort = effort;
        Ok(())
    }
    fn is_ancestor(&self, joint: usize, body: usize) -> bool {
        let mut current = Some(body);
        while let Some(c) = current {
            if c == joint {
                return true;
            }
            current = self.links[c].parent;
        }
        false
    }
    fn linear_jacobian(&self, joint: usize, body: usize) -> Vec3 {
        if !self.is_ancestor(joint, body) {
            return Vec3::ZERO;
        }
        if self.links[joint].joint_type == ReducedJointType::Prismatic {
            return self.axis[joint];
        }
        cross(
            self.axis[joint],
            self.links[body].body.state.position - self.origin[joint],
        )
    }
    fn angular_jacobian(&self, joint: usize, body: usize) -> Vec3 {
        if !self.is_ancestor(joint, body)
            || self.links[joint].joint_type == ReducedJointType::Prismatic
        {
            return Vec3::ZERO;
        }
        self.axis[joint]
    }
    pub fn forward_kinematics(&mut self) -> Result<()> {
        for i in 0..self.links.len() {
            let link = self.links[i];
            let (pp, po) = match link.parent {
                Some(p) => (
                    self.links[p].body.state.position,
                    self.links[p].body.state.orientation,
                ),
                None => (Vec3::ZERO, Quaternion::default()),
            };
            let axis_world = rotate(po, link.axis_parent)?;
            let mut origin = pp + rotate(po, link.parent_joint_offset)?;
            let child = if link.joint_type == ReducedJointType::Revolute {
                normalized_quaternion(
                    po * quaternion_from_axis_angle(link.axis_parent, link.q)?
                        * link.reference_orientation,
                )?
            } else {
                origin += axis_world * link.q;
                normalized_quaternion(po * link.reference_orientation)?
            };
            self.origin[i] = origin;
            self.axis[i] = normalized_or_zero(axis_world, 1e-15);
            self.links[i].body.state.orientation = child;
            self.links[i].body.state.position = origin + rotate(child, link.joint_to_com_local)?;
        }
        for b in 0..self.links.len() {
            let (mut v, mut w) = (Vec3::ZERO, Vec3::ZERO);
            for j in 0..self.links.len() {
                if !self.is_ancestor(j, b) {
                    continue;
                }
                v += self.linear_jacobian(j, b) * self.links[j].qd;
                w += self.angular_jacobian(j, b) * self.links[j].qd;
            }
            self.links[b].body.state.linear_velocity = v;
            self.links[b].body.state.angular_velocity = w;
        }
        Ok(())
    }
    pub fn mass_matrix(&self) -> Result<Vec<f64>> {
        let n = self.links.len();
        let mut m = vec![0.0; n * n];
        for b in 0..n {
            let rb = &self.links[b].body;
            for i in 0..n {
                if !self.is_ancestor(i, b) {
                    continue;
                }
                let (jvi, jwi) = (self.linear_jacobian(i, b), self.angular_jacobian(i, b));
                for j in 0..n {
                    if !self.is_ancestor(j, b) {
                        continue;
                    }
                    let (jvj, jwj) = (self.linear_jacobian(j, b), self.angular_jacobian(j, b));
                    m[i * n + j] += rb.mass * dot(jvi, jvj) + dot(jwi, world_inertia_mul(rb, jwj)?);
                }
            }
        }
        Ok(m)
    }
    pub fn generalized_forces(&self) -> Vec<f64> {
        let n = self.links.len();
        (0..n)
            .map(|i| {
                let mut f = self.links[i].applied_effort - self.links[i].damping * self.links[i].qd;
                for b in 0..n {
                    if !self.is_ancestor(i, b) {
                        continue;
                    }
                    let body = &self.links[b].body;
                    f += body.mass * dot(self.linear_jacobian(i, b), self.gravity);
                    f += dot(self.angular_jacobian(i, b), body.torque);
                    f += dot(self.linear_jacobian(i, b), body.force);
                }
                f
            })
            .collect()
    }
    /// Semi-implicit Euler on the generalized coordinates; on error the system is unchanged.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(
            dt > 0.0 && dt.is_finite(),
            "articulated dt must be finite and positive",
        )?;
        let mut t = self.clone();
        t.forward_kinematics()?;
        let qdd = solve_dense(t.mass_matrix()?, t.generalized_forces(), t.links.len())?;
        for (l, a) in t.links.iter_mut().zip(&qdd) {
            l.qdd = *a;
            l.qd += dt * a;
            l.q += dt * l.qd;
            l.body.clear_accumulators();
        }
        t.forward_kinematics()?;
        *self = t;
        Ok(())
    }
    pub fn run(&mut self, steps: usize, dt: f64) -> Result<()> {
        for _ in 0..steps {
            self.step(dt)?;
        }
        Ok(())
    }
    pub fn links(&self) -> &[ArticulatedLink] {
        &self.links
    }
    pub fn links_mut(&mut self) -> &mut [ArticulatedLink] {
        &mut self.links
    }
    pub fn kinetic_energy(&self) -> Result<f64> {
        let mut e = 0.0;
        for l in &self.links {
            let b = &l.body;
            e += 0.5 * b.mass * norm_squared(b.state.linear_velocity)
                + 0.5
                    * dot(
                        b.state.angular_velocity,
                        world_inertia_mul(b, b.state.angular_velocity)?,
                    );
        }
        Ok(e)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FlexibleInterfaceNode {
    pub reference_local: Vec3,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FlexibleMode {
    /// One displacement vector per interface node.
    pub shape_local: Vec<Vec3>,
    pub modal_mass: f64,
    pub modal_stiffness: f64,
    pub modal_damping: f64,
}

impl Default for FlexibleMode {
    fn default() -> Self {
        Self {
            shape_local: vec![],
            modal_mass: 1.0,
            modal_stiffness: 0.0,
            modal_damping: 0.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct FlexibleInterfaceLoad {
    pub body_force: Vec3,
    pub body_torque: Vec3,
    pub modal_force: Vec<f64>,
}

/// Conservative floating-frame map between FEM interface nodes and a rigid frame plus modes.
#[derive(Clone, Debug, PartialEq)]
pub struct FlexibleBodyInterface {
    nodes: Vec<FlexibleInterfaceNode>,
    modes: Vec<FlexibleMode>,
    q: Vec<f64>,
    qd: Vec<f64>,
}

impl FlexibleBodyInterface {
    pub fn new(nodes: Vec<FlexibleInterfaceNode>, modes: Vec<FlexibleMode>) -> Result<Self> {
        require(!nodes.is_empty(), "flexible interface requires nodes")?;
        for m in &modes {
            require(
                m.shape_local.len() == nodes.len(),
                "flexible mode shape size mismatch",
            )?;
            require(
                m.modal_mass > 0.0 && m.modal_mass.is_finite(),
                "flexible modal mass must be finite and positive",
            )?;
        }
        let k = modes.len();
        Ok(Self {
            nodes,
            modes,
            q: vec![0.0; k],
            qd: vec![0.0; k],
        })
    }
    /// `qd` may be empty.
    pub fn set_modal_state(&mut self, q: Vec<f64>, qd: Vec<f64>) -> Result<()> {
        require(
            q.len() == self.modes.len(),
            "flexible modal q size mismatch",
        )?;
        let qd = if qd.is_empty() {
            vec![0.0; self.modes.len()]
        } else {
            qd
        };
        require(
            qd.len() == self.modes.len(),
            "flexible modal qd size mismatch",
        )?;
        self.q = q;
        self.qd = qd;
        Ok(())
    }
    fn deformation(&self, node: usize, amp: &[f64]) -> Vec3 {
        let mut d = Vec3::ZERO;
        for (m, a) in self.modes.iter().zip(amp) {
            d += m.shape_local[node] * *a;
        }
        d
    }
    pub fn node_positions(&self, frame: &RigidBody) -> Result<Vec<Vec3>> {
        (0..self.nodes.len())
            .map(|i| {
                Ok(frame.state.position
                    + rotate(
                        frame.state.orientation,
                        self.nodes[i].reference_local + self.deformation(i, &self.q),
                    )?)
            })
            .collect()
    }
    pub fn node_velocities(&self, frame: &RigidBody) -> Result<Vec<Vec3>> {
        (0..self.nodes.len())
            .map(|i| {
                let r = rotate(
                    frame.state.orientation,
                    self.nodes[i].reference_local + self.deformation(i, &self.q),
                )?;
                Ok(frame.state.linear_velocity
                    + cross(frame.state.angular_velocity, r)
                    + rotate(frame.state.orientation, self.deformation(i, &self.qd))?)
            })
            .collect()
    }
    pub fn project_nodal_forces(
        &self,
        frame: &RigidBody,
        nodal_forces: &[Vec3],
    ) -> Result<FlexibleInterfaceLoad> {
        require(
            nodal_forces.len() == self.nodes.len(),
            "flexible nodal force size mismatch",
        )?;
        let mut load = FlexibleInterfaceLoad {
            modal_force: vec![0.0; self.modes.len()],
            ..Default::default()
        };
        let positions = self.node_positions(frame)?;
        for (i, f) in nodal_forces.iter().enumerate() {
            load.body_force += *f;
            load.body_torque += cross(positions[i] - frame.state.position, *f);
            for (m, mode) in self.modes.iter().enumerate() {
                load.modal_force[m] +=
                    dot(rotate(frame.state.orientation, mode.shape_local[i])?, *f);
            }
        }
        Ok(load)
    }
    pub fn apply_nodal_forces(
        &self,
        frame: &mut RigidBody,
        nodal_forces: &[Vec3],
    ) -> Result<FlexibleInterfaceLoad> {
        let load = self.project_nodal_forces(frame, nodal_forces)?;
        frame.add_force(load.body_force);
        frame.add_torque(load.body_torque);
        Ok(load)
    }
    pub fn modal_q(&self) -> &[f64] {
        &self.q
    }
    pub fn modal_qd(&self) -> &[f64] {
        &self.qd
    }
    pub fn modes(&self) -> &[FlexibleMode] {
        &self.modes
    }
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}
