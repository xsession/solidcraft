//! Explicit two-fluid (Euler-Euler) transport with interphase drag, ported from upstream
//! `cfd/solvers/fvm/euler_euler` (MIT, cfd_solvers contributors).
use super::{particles::finite_vec, PolyMesh};
use crate::{core::Vec3, require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerEulerConfig {
    pub dt: f64,
    pub primary_density: f64,
    pub dispersed_density: f64,
    pub drag_coefficient: f64,
    pub quadratic_drag: bool,
    pub primary_body_acceleration: Vec3,
    pub dispersed_body_acceleration: Vec3,
    pub minimum_phase_fraction: f64,
}

impl Default for EulerEulerConfig {
    fn default() -> Self {
        Self {
            dt: 1e-3,
            primary_density: 1000.0,
            dispersed_density: 1.0,
            drag_coefficient: 1.0,
            quadratic_drag: false,
            primary_body_acceleration: Vec3::default(),
            dispersed_body_acceleration: Vec3::default(),
            minimum_phase_fraction: 1e-12,
        }
    }
}

/// Conservative upwind phase-fraction transport with phase momenta and equal-and-opposite drag.
/// Boundary faces are impermeable, as upstream.
#[derive(Clone, Debug)]
pub struct EulerEulerTransport {
    mesh: PolyMesh,
    config: EulerEulerConfig,
    alpha: Vec<f64>,
    primary_velocity: Vec<Vec3>,
    dispersed_velocity: Vec<Vec3>,
    primary_face_flux: Vec<f64>,
    dispersed_face_flux: Vec<f64>,
    drag_force: Vec<Vec3>,
    time: f64,
    steps: usize,
}

fn valid_fraction(a: f64) -> bool {
    a.is_finite() && (0.0..=1.0).contains(&a)
}

impl EulerEulerTransport {
    pub fn new(mesh: PolyMesh, config: EulerEulerConfig) -> Result<Self> {
        let c = config;
        require(
            c.dt > 0.0
                && c.dt.is_finite()
                && c.primary_density > 0.0
                && c.primary_density.is_finite()
                && c.dispersed_density > 0.0
                && c.dispersed_density.is_finite()
                && c.drag_coefficient >= 0.0
                && c.drag_coefficient.is_finite()
                && c.minimum_phase_fraction > 0.0
                && c.minimum_phase_fraction < 0.5
                && finite_vec(c.primary_body_acceleration)
                && finite_vec(c.dispersed_body_acceleration),
            "invalid Euler-Euler controls",
        )?;
        let (n, f) = (mesh.cell_count(), mesh.face_count());
        Ok(Self {
            mesh,
            config,
            alpha: vec![0.0; n],
            primary_velocity: vec![Vec3::default(); n],
            dispersed_velocity: vec![Vec3::default(); n],
            primary_face_flux: vec![0.0; f],
            dispersed_face_flux: vec![0.0; f],
            drag_force: vec![Vec3::default(); n],
            time: 0.0,
            steps: 0,
        })
    }

    pub fn initialize_uniform(
        &mut self,
        dispersed_fraction: f64,
        primary_velocity: Vec3,
        dispersed_velocity: Vec3,
    ) -> Result<()> {
        let n = self.mesh.cell_count();
        self.initialize(
            &vec![dispersed_fraction; n],
            &vec![primary_velocity; n],
            &vec![dispersed_velocity; n],
        )
    }

    pub fn initialize_profile(
        &mut self,
        fraction: impl Fn(Vec3) -> f64,
        primary_velocity: Vec3,
        dispersed_velocity: Vec3,
    ) -> Result<()> {
        let a: Vec<f64> = self
            .mesh
            .cells()
            .iter()
            .map(|c| fraction(c.center))
            .collect();
        let n = a.len();
        self.initialize(&a, &vec![primary_velocity; n], &vec![dispersed_velocity; n])
    }

    pub fn initialize(
        &mut self,
        dispersed_fraction: &[f64],
        primary_velocity: &[Vec3],
        dispersed_velocity: &[Vec3],
    ) -> Result<()> {
        let n = self.mesh.cell_count();
        require(
            dispersed_fraction.len() == n
                && primary_velocity.len() == n
                && dispersed_velocity.len() == n,
            "Euler-Euler initial field size mismatch",
        )?;
        require(
            dispersed_fraction.iter().all(|a| valid_fraction(*a)),
            "Euler-Euler volume fraction must lie in [0,1]",
        )?;
        require(
            primary_velocity
                .iter()
                .chain(dispersed_velocity)
                .all(|v| finite_vec(*v)),
            "Euler-Euler velocity must be finite",
        )?;
        self.alpha.copy_from_slice(dispersed_fraction);
        self.primary_velocity.copy_from_slice(primary_velocity);
        self.dispersed_velocity.copy_from_slice(dispersed_velocity);
        self.primary_face_flux.fill(0.0);
        self.dispersed_face_flux.fill(0.0);
        self.drag_force.fill(Vec3::default());
        self.time = 0.0;
        self.steps = 0;
        Ok(())
    }

    pub fn step(&mut self) {
        let c = self.config;
        let mesh = &self.mesh;
        let n = mesh.cell_count();
        let mut primary_momentum_flux = vec![Vec3::default(); mesh.face_count()];
        let mut dispersed_momentum_flux = vec![Vec3::default(); mesh.face_count()];
        let mut primary_momentum = vec![Vec3::default(); n];
        let mut dispersed_momentum = vec![Vec3::default(); n];
        for cell in 0..n {
            let alpha = self.alpha[cell];
            let beta = 1.0 - alpha;
            let vol = mesh.cells()[cell].volume;
            primary_momentum[cell] = self.primary_velocity[cell] * (c.primary_density * beta * vol);
            dispersed_momentum[cell] =
                self.dispersed_velocity[cell] * (c.dispersed_density * alpha * vol);
            let relative = self.dispersed_velocity[cell] - self.primary_velocity[cell];
            let mut coefficient = c.drag_coefficient * (alpha * beta);
            if c.quadratic_drag {
                coefficient *= relative.norm();
            }
            self.drag_force[cell] = relative * coefficient;
        }
        self.primary_face_flux.fill(0.0);
        self.dispersed_face_flux.fill(0.0);
        for (fi, face) in mesh.faces().iter().enumerate() {
            let Some(nb) = face.neighbour else { continue };
            let ow = face.owner;
            let pfv = (self.primary_velocity[ow] + self.primary_velocity[nb]) * 0.5;
            let dfv = (self.dispersed_velocity[ow] + self.dispersed_velocity[nb]) * 0.5;
            let pn = pfv.dot(face.area);
            let dn = dfv.dot(face.area);
            let p_owner = pn >= 0.0;
            let d_owner = dn >= 0.0;
            let alpha_face = if d_owner {
                self.alpha[ow]
            } else {
                self.alpha[nb]
            };
            let beta_face = if p_owner {
                1.0 - self.alpha[ow]
            } else {
                1.0 - self.alpha[nb]
            };
            let pf = beta_face * pn;
            let df = alpha_face * dn;
            self.primary_face_flux[fi] = pf;
            self.dispersed_face_flux[fi] = df;
            let pu = if p_owner {
                self.primary_velocity[ow]
            } else {
                self.primary_velocity[nb]
            };
            let du = if d_owner {
                self.dispersed_velocity[ow]
            } else {
                self.dispersed_velocity[nb]
            };
            primary_momentum_flux[fi] = pu * (c.primary_density * pf);
            dispersed_momentum_flux[fi] = du * (c.dispersed_density * df);
        }
        let mut next_alpha = self.alpha.clone();
        let mut next_primary = self.primary_velocity.clone();
        let mut next_dispersed = self.dispersed_velocity.clone();
        for cell in 0..n {
            let mut d_net = 0.0;
            let mut p_mom = Vec3::default();
            let mut d_mom = Vec3::default();
            for &fi in &mesh.cell_faces()[cell] {
                let o = if mesh.faces()[fi].owner == cell {
                    1.0
                } else {
                    -1.0
                };
                d_net += o * self.dispersed_face_flux[fi];
                p_mom += primary_momentum_flux[fi] * o;
                d_mom += dispersed_momentum_flux[fi] * o;
            }
            let vol = mesh.cells()[cell].volume;
            next_alpha[cell] = (self.alpha[cell] - c.dt * d_net / vol).clamp(0.0, 1.0);
            let alpha = self.alpha[cell];
            let beta = 1.0 - alpha;
            let p_mass = c.primary_density * beta * vol;
            let d_mass = c.dispersed_density * alpha * vol;
            let p_force = c.primary_body_acceleration * p_mass + self.drag_force[cell] * vol;
            let d_force = c.dispersed_body_acceleration * d_mass - self.drag_force[cell] * vol;
            primary_momentum[cell] -= p_mom * c.dt;
            dispersed_momentum[cell] -= d_mom * c.dt;
            primary_momentum[cell] += p_force * c.dt;
            dispersed_momentum[cell] += d_force * c.dt;
            let np_mass = c.primary_density * (1.0 - next_alpha[cell]) * vol;
            let nd_mass = c.dispersed_density * next_alpha[cell] * vol;
            next_primary[cell] = if np_mass > c.minimum_phase_fraction * c.primary_density * vol {
                primary_momentum[cell] / np_mass
            } else {
                Vec3::default()
            };
            next_dispersed[cell] = if nd_mass > c.minimum_phase_fraction * c.dispersed_density * vol
            {
                dispersed_momentum[cell] / nd_mass
            } else {
                Vec3::default()
            };
        }
        self.alpha = next_alpha;
        self.primary_velocity = next_primary;
        self.dispersed_velocity = next_dispersed;
        self.time += c.dt;
        self.steps += 1;
    }

    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn config(&self) -> &EulerEulerConfig {
        &self.config
    }
    pub fn dispersed_fraction(&self) -> &[f64] {
        &self.alpha
    }
    pub fn primary_velocity(&self) -> &[Vec3] {
        &self.primary_velocity
    }
    pub fn dispersed_velocity(&self) -> &[Vec3] {
        &self.dispersed_velocity
    }
    pub fn primary_face_flux(&self) -> &[f64] {
        &self.primary_face_flux
    }
    pub fn dispersed_face_flux(&self) -> &[f64] {
        &self.dispersed_face_flux
    }
    pub fn drag_force(&self) -> &[Vec3] {
        &self.drag_force
    }
    pub fn volume_fraction_integral(&self) -> f64 {
        self.alpha
            .iter()
            .zip(self.mesh.cells())
            .map(|(a, c)| a * c.volume)
            .sum()
    }
    pub fn total_momentum(&self) -> Vec3 {
        let c = self.config;
        let mut r = Vec3::default();
        for (i, cell) in self.mesh.cells().iter().enumerate() {
            r += self.primary_velocity[i]
                * (c.primary_density * (1.0 - self.alpha[i]) * cell.volume);
            r += self.dispersed_velocity[i] * (c.dispersed_density * self.alpha[i] * cell.volume);
        }
        r
    }
    pub fn integrated_drag_force(&self) -> Vec3 {
        let mut r = Vec3::default();
        for (f, cell) in self.drag_force.iter().zip(self.mesh.cells()) {
            r += *f * cell.volume;
        }
        r
    }
    pub fn minimum_fraction(&self) -> f64 {
        self.alpha.iter().copied().fold(f64::INFINITY, f64::min)
    }
    pub fn maximum_fraction(&self) -> f64 {
        self.alpha.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn steps(&self) -> usize {
        self.steps
    }
}
