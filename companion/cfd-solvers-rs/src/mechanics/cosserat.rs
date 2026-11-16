//! Stable Cosserat rod: alternating position (stretch/shear) and material-frame (bend/twist)
//! relaxation with XPBD-style compliance.
use crate::{
    multibody::{
        cross, div, integrate_orientation, norm, norm_squared, normalized_or_zero,
        normalized_quaternion, quaternion_from_axis_angle, rotate, Quaternion, Vec3,
    },
    require, Result, SolverError,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CosseratRodConfig {
    pub dt: f64,
    pub gravity: Vec3,
    pub linear_damping: f64,
    pub angular_damping: f64,
    pub stretch_shear_compliance: f64,
    pub bend_twist_compliance: f64,
    pub stretch_shear_stiffness: f64,
    pub bend_twist_stiffness: f64,
    pub director_alignment_strength: f64,
    pub alternating_iterations: usize,
    pub position_iterations: usize,
    pub orientation_iterations: usize,
}

impl Default for CosseratRodConfig {
    fn default() -> Self {
        Self {
            dt: 1e-3,
            gravity: Vec3::new(0.0, -9.81, 0.0),
            linear_damping: 0.02,
            angular_damping: 0.02,
            stretch_shear_compliance: 1e-8,
            bend_twist_compliance: 1e-8,
            stretch_shear_stiffness: 1e5,
            bend_twist_stiffness: 1e2,
            director_alignment_strength: 1.0,
            alternating_iterations: 4,
            position_iterations: 2,
            orientation_iterations: 2,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CosseratRodNode {
    pub position: Vec3,
    pub orientation: Quaternion,
    pub linear_velocity: Vec3,
    pub angular_velocity: Vec3,
    pub inverse_mass: f64,
    pub inverse_rotational_mass: f64,
    pub fixed: bool,
}

impl Default for CosseratRodNode {
    fn default() -> Self {
        Self {
            position: Vec3::ZERO,
            orientation: Quaternion::default(),
            linear_velocity: Vec3::ZERO,
            angular_velocity: Vec3::ZERO,
            inverse_mass: 1.0,
            inverse_rotational_mass: 1.0,
            fixed: false,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CosseratRodDiagnostics {
    pub max_stretch_shear_error: f64,
    pub max_bend_twist_error: f64,
    pub kinetic_energy: f64,
    pub elastic_energy: f64,
    pub linear_momentum_norm: f64,
    pub finite: bool,
}

fn finite_v(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

fn finite_q(q: Quaternion) -> bool {
    q.w.is_finite() && q.x.is_finite() && q.y.is_finite() && q.z.is_finite()
}

fn rotation_from_vector(r: Vec3) -> Result<Quaternion> {
    let angle = norm(r);
    if angle.is_nan() || angle <= 1e-14 {
        return normalized_quaternion(Quaternion {
            w: 1.0,
            x: 0.5 * r.x,
            y: 0.5 * r.y,
            z: 0.5 * r.z,
        });
    }
    quaternion_from_axis_angle(div(r, angle), angle)
}

fn rotation_vector(q: Quaternion) -> Result<Vec3> {
    let mut q = normalized_quaternion(q)?;
    if q.w < 0.0 {
        q = Quaternion {
            w: -q.w,
            x: -q.x,
            y: -q.y,
            z: -q.z,
        };
    }
    let v = Vec3::new(q.x, q.y, q.z);
    let s = norm(v);
    if s.is_nan() || s <= 1e-14 {
        return Ok(v * 2.0);
    }
    let angle = 2.0 * s.atan2(q.w.clamp(-1.0, 1.0));
    Ok(v * (angle / s))
}

fn safe_weight(v: f64) -> f64 {
    if v > 0.0 && v.is_finite() {
        v
    } else {
        0.0
    }
}

const EZ: Vec3 = Vec3::new(0.0, 0.0, 1.0);

#[derive(Clone, Debug, PartialEq)]
pub struct StableCosseratRod {
    config: CosseratRodConfig,
    nodes: Vec<CosseratRodNode>,
    rest_lengths: Vec<f64>,
    rest_relative: Vec<Quaternion>,
    forces: Vec<Vec3>,
    torques: Vec<Vec3>,
    diagnostics: CosseratRodDiagnostics,
    time: f64,
    steps: usize,
}

impl StableCosseratRod {
    pub fn new(config: CosseratRodConfig) -> Result<Self> {
        let s = Self {
            config,
            nodes: vec![],
            rest_lengths: vec![],
            rest_relative: vec![],
            forces: vec![],
            torques: vec![],
            diagnostics: CosseratRodDiagnostics::default(),
            time: 0.0,
            steps: 0,
        };
        s.validate_config()?;
        Ok(s)
    }
    fn validate_config(&self) -> Result<()> {
        let c = &self.config;
        require(
            c.dt > 0.0 && c.dt.is_finite(),
            "Cosserat rod dt must be positive and finite",
        )?;
        require(finite_v(c.gravity), "Cosserat rod gravity must be finite")?;
        require(
            c.linear_damping >= 0.0
                && c.angular_damping >= 0.0
                && c.linear_damping.is_finite()
                && c.angular_damping.is_finite(),
            "Cosserat rod damping must be finite and non-negative",
        )?;
        require(
            c.stretch_shear_compliance >= 0.0
                && c.bend_twist_compliance >= 0.0
                && c.stretch_shear_compliance.is_finite()
                && c.bend_twist_compliance.is_finite(),
            "Cosserat rod compliance must be finite and non-negative",
        )?;
        require(
            c.stretch_shear_stiffness > 0.0
                && c.bend_twist_stiffness > 0.0
                && c.director_alignment_strength >= 0.0
                && c.stretch_shear_stiffness.is_finite()
                && c.bend_twist_stiffness.is_finite()
                && c.director_alignment_strength.is_finite(),
            "Cosserat rod stiffness/alignment settings must be finite and valid",
        )?;
        require(
            c.alternating_iterations > 0
                && c.position_iterations > 0
                && c.orientation_iterations > 0,
            "Cosserat rod iteration counts must be non-zero",
        )
    }
    fn validate_state(&self) -> Result<()> {
        require(
            self.nodes.len() >= 2,
            "Cosserat rod needs at least two nodes",
        )?;
        require(
            self.rest_lengths.len() + 1 == self.nodes.len()
                && self.rest_relative.len() == self.rest_lengths.len(),
            "Cosserat rod rest state is inconsistent",
        )?;
        for n in &self.nodes {
            if !(finite_v(n.position)
                && finite_q(n.orientation)
                && finite_v(n.linear_velocity)
                && finite_v(n.angular_velocity)
                && n.inverse_mass >= 0.0
                && n.inverse_rotational_mass >= 0.0
                && n.inverse_mass.is_finite()
                && n.inverse_rotational_mass.is_finite())
            {
                return Err(SolverError::NumericalFailure(
                    "Cosserat rod contains a non-finite or invalid node".into(),
                ));
            }
        }
        require(
            self.rest_lengths.iter().all(|l| *l > 0.0 && l.is_finite()),
            "Cosserat rod rest length must be positive and finite",
        )
    }
    /// Orientations may be empty (identity frames); positions define rest lengths and the
    /// orientations define rest relative frames.
    pub fn initialize(&mut self, positions: &[Vec3], orientations: &[Quaternion]) -> Result<()> {
        self.validate_config()?;
        require(
            positions.len() >= 2,
            "Cosserat rod needs at least two positions",
        )?;
        require(
            orientations.is_empty() || orientations.len() == positions.len(),
            "Cosserat rod orientation count must match positions",
        )?;
        let mut nodes = vec![CosseratRodNode::default(); positions.len()];
        for (i, p) in positions.iter().enumerate() {
            require(finite_v(*p), "Cosserat rod position is non-finite")?;
            nodes[i].position = *p;
            nodes[i].orientation = if orientations.is_empty() {
                Quaternion::default()
            } else {
                normalized_quaternion(orientations[i])?
            };
        }
        let mut rl = Vec::with_capacity(nodes.len() - 1);
        let mut rr = Vec::with_capacity(nodes.len() - 1);
        for i in 0..nodes.len() - 1 {
            let l = norm(nodes[i + 1].position - nodes[i].position);
            require(
                l > 0.0 && l.is_finite(),
                "Cosserat rod positions must be distinct",
            )?;
            rl.push(l);
            rr.push(normalized_quaternion(
                nodes[i].orientation.conjugate() * nodes[i + 1].orientation,
            )?);
        }
        self.forces = vec![Vec3::ZERO; nodes.len()];
        self.torques = vec![Vec3::ZERO; nodes.len()];
        self.nodes = nodes;
        self.rest_lengths = rl;
        self.rest_relative = rr;
        self.time = 0.0;
        self.steps = 0;
        self.update_diagnostics()
    }
    pub fn set_fixed(&mut self, node: usize, fixed: bool) -> Result<()> {
        require(node < self.nodes.len(), "Cosserat rod node index")?;
        let n = &mut self.nodes[node];
        n.fixed = fixed;
        if fixed {
            n.linear_velocity = Vec3::ZERO;
            n.angular_velocity = Vec3::ZERO;
        }
        Ok(())
    }
    pub fn add_force(&mut self, node: usize, force: Vec3) -> Result<()> {
        require(node < self.nodes.len(), "Cosserat rod node index")?;
        require(finite_v(force), "Cosserat rod force is non-finite")?;
        self.forces[node] += force;
        Ok(())
    }
    pub fn add_torque(&mut self, node: usize, torque: Vec3) -> Result<()> {
        require(node < self.nodes.len(), "Cosserat rod node index")?;
        require(finite_v(torque), "Cosserat rod torque is non-finite")?;
        self.torques[node] += torque;
        Ok(())
    }
    pub fn clear_external_loads(&mut self) {
        self.forces.iter_mut().for_each(|f| *f = Vec3::ZERO);
        self.torques.iter_mut().for_each(|t| *t = Vec3::ZERO);
    }
    fn relax_positions(&mut self) -> Result<()> {
        let alpha = self.config.stretch_shear_compliance / (self.config.dt * self.config.dt);
        for _ in 0..self.config.position_iterations {
            for i in 0..self.nodes.len() - 1 {
                let j = i + 1;
                let di = rotate(self.nodes[i].orientation, EZ)?;
                let dj = rotate(self.nodes[j].orientation, EZ)?;
                let target = (di + dj) * (0.5 * self.rest_lengths[i]);
                let constraint = (self.nodes[j].position - self.nodes[i].position) - target;
                let wi = if self.nodes[i].fixed {
                    0.0
                } else {
                    safe_weight(self.nodes[i].inverse_mass)
                };
                let wj = if self.nodes[j].fixed {
                    0.0
                } else {
                    safe_weight(self.nodes[j].inverse_mass)
                };
                let den = wi + wj + alpha;
                if den.is_nan() || den <= 0.0 {
                    continue;
                }
                let corr = div(constraint, den);
                if wi > 0.0 {
                    self.nodes[i].position += corr * wi;
                }
                if wj > 0.0 {
                    self.nodes[j].position -= corr * wj;
                }
            }
        }
        Ok(())
    }
    fn relax_orientations(&mut self) -> Result<()> {
        let alpha = self.config.bend_twist_compliance / (self.config.dt * self.config.dt);
        for _ in 0..self.config.orientation_iterations {
            for i in 0..self.nodes.len() - 1 {
                let j = i + 1;
                let (oi, oj) = (self.nodes[i].orientation, self.nodes[j].orientation);
                let relative = normalized_quaternion(oi.conjugate() * oj)?;
                let error = normalized_quaternion(relative * self.rest_relative[i].conjugate())?;
                let bend_world = rotate(oi, rotation_vector(error)?)?;
                let di = rotate(oi, EZ)?;
                let dj = rotate(oj, EZ)?;
                let tangent =
                    normalized_or_zero(self.nodes[j].position - self.nodes[i].position, 1e-15);
                let avg = normalized_or_zero(di + dj, 1e-15);
                let align = cross(avg, tangent) * self.config.director_alignment_strength;
                let wi = if self.nodes[i].fixed {
                    0.0
                } else {
                    safe_weight(self.nodes[i].inverse_rotational_mass)
                };
                let wj = if self.nodes[j].fixed {
                    0.0
                } else {
                    safe_weight(self.nodes[j].inverse_rotational_mass)
                };
                let den = wi + wj + alpha;
                if den.is_nan() || den <= 0.0 {
                    continue;
                }
                let bc = div(bend_world, den);
                let ac = div(align, den);
                if wi > 0.0 {
                    self.nodes[i].orientation = normalized_quaternion(
                        rotation_from_vector((bc + ac) * wi)? * self.nodes[i].orientation,
                    )?;
                }
                if wj > 0.0 {
                    self.nodes[j].orientation = normalized_quaternion(
                        rotation_from_vector((-bc + ac) * wj)? * self.nodes[j].orientation,
                    )?;
                }
            }
        }
        Ok(())
    }
    /// One step; on error the rod is unchanged.
    pub fn step(&mut self) -> Result<()> {
        let mut t = self.clone();
        t.step_inner()?;
        *self = t;
        Ok(())
    }
    fn step_inner(&mut self) -> Result<()> {
        self.validate_config()?;
        self.validate_state()?;
        let c = self.config;
        let prev: Vec<(Vec3, Quaternion)> = self
            .nodes
            .iter()
            .map(|n| (n.position, n.orientation))
            .collect();
        for i in 0..self.nodes.len() {
            let n = &mut self.nodes[i];
            if n.fixed {
                n.linear_velocity = Vec3::ZERO;
                n.angular_velocity = Vec3::ZERO;
                continue;
            }
            n.linear_velocity += (c.gravity + self.forces[i] * n.inverse_mass) * c.dt;
            n.position += n.linear_velocity * c.dt;
            n.angular_velocity += self.torques[i] * (n.inverse_rotational_mass * c.dt);
            integrate_orientation(&mut n.orientation, n.angular_velocity, c.dt)?;
        }
        for _ in 0..c.alternating_iterations {
            self.relax_positions()?;
            self.relax_orientations()?;
        }
        let ld = (-c.linear_damping * c.dt).exp();
        let ad = (-c.angular_damping * c.dt).exp();
        for (n, (pp, po)) in self.nodes.iter_mut().zip(&prev) {
            if n.fixed {
                n.position = *pp;
                n.orientation = *po;
                n.linear_velocity = Vec3::ZERO;
                n.angular_velocity = Vec3::ZERO;
                continue;
            }
            n.linear_velocity = div(n.position - *pp, c.dt) * ld;
            let delta = normalized_quaternion(n.orientation * po.conjugate())?;
            n.angular_velocity = div(rotation_vector(delta)?, c.dt) * ad;
        }
        self.clear_external_loads();
        self.time += c.dt;
        self.steps += 1;
        self.update_diagnostics()
    }
    pub fn run(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }
    fn update_diagnostics(&mut self) -> Result<()> {
        let mut d = CosseratRodDiagnostics::default();
        let (mut mx, mut my, mut mz) = (0.0, 0.0, 0.0);
        for n in &self.nodes {
            let m = if n.inverse_mass > 0.0 {
                1.0 / n.inverse_mass
            } else {
                0.0
            };
            let r = if n.inverse_rotational_mass > 0.0 {
                1.0 / n.inverse_rotational_mass
            } else {
                0.0
            };
            d.kinetic_energy += 0.5 * m * norm_squared(n.linear_velocity);
            d.kinetic_energy += 0.5 * r * norm_squared(n.angular_velocity);
            mx += m * n.linear_velocity.x;
            my += m * n.linear_velocity.y;
            mz += m * n.linear_velocity.z;
        }
        d.linear_momentum_norm = norm(Vec3::new(mx, my, mz));
        for i in 0..self.nodes.len().saturating_sub(1) {
            let (a, b) = (&self.nodes[i], &self.nodes[i + 1]);
            let target = (rotate(a.orientation, EZ)? + rotate(b.orientation, EZ)?)
                * (0.5 * self.rest_lengths[i]);
            let se = norm((b.position - a.position) - target);
            d.max_stretch_shear_error = d.max_stretch_shear_error.max(se);
            let rel = normalized_quaternion(a.orientation.conjugate() * b.orientation)?;
            let err = normalized_quaternion(rel * self.rest_relative[i].conjugate())?;
            let be = norm(rotation_vector(err)?);
            d.max_bend_twist_error = d.max_bend_twist_error.max(be);
            d.elastic_energy += 0.5 * self.config.stretch_shear_stiffness * se * se;
            d.elastic_energy += 0.5 * self.config.bend_twist_stiffness * be * be;
        }
        d.finite = [
            d.max_stretch_shear_error,
            d.max_bend_twist_error,
            d.kinetic_energy,
            d.elastic_energy,
            d.linear_momentum_norm,
        ]
        .iter()
        .all(|v| v.is_finite());
        self.diagnostics = d;
        Ok(())
    }
    pub fn config(&self) -> &CosseratRodConfig {
        &self.config
    }
    pub fn nodes(&self) -> &[CosseratRodNode] {
        &self.nodes
    }
    pub fn nodes_mut(&mut self) -> &mut [CosseratRodNode] {
        &mut self.nodes
    }
    pub fn rest_lengths(&self) -> &[f64] {
        &self.rest_lengths
    }
    pub fn diagnostics(&self) -> &CosseratRodDiagnostics {
        &self.diagnostics
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn steps(&self) -> usize {
        self.steps
    }
    pub fn checksum(&self) -> f64 {
        let mut r = 0.0;
        for (i, n) in self.nodes.iter().enumerate() {
            let w = (i + 1) as f64;
            r += w
                * (n.position.x
                    + 3.0 * n.position.y
                    + 5.0 * n.position.z
                    + 7.0 * n.orientation.w
                    + 11.0 * n.orientation.x
                    + 13.0 * n.orientation.y
                    + 17.0 * n.orientation.z);
        }
        r
    }
}
