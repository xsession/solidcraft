use super::minmod;
use crate::{positive, require, Result, SolverError};
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct IdealGasEquationOfState {
    pub gamma: f64,
    pub gas_constant: f64,
}
impl Default for IdealGasEquationOfState {
    fn default() -> Self {
        Self {
            gamma: 1.4,
            gas_constant: 287.05,
        }
    }
}
impl IdealGasEquationOfState {
    fn validate(self) -> Result<()> {
        require(
            self.gamma.is_finite() && self.gamma > 1.0,
            "gamma must exceed one",
        )?;
        positive(self.gas_constant, "gas constant")
    }
    pub fn pressure(self, density: f64, momentum: f64, total_energy: f64) -> Result<f64> {
        self.validate()?;
        positive(density, "density")?;
        let p = (self.gamma - 1.0) * (total_energy - 0.5 * momentum * momentum / density);
        positive(p, "pressure")?;
        Ok(p)
    }
    pub fn sound_speed(self, density: f64, pressure: f64) -> Result<f64> {
        self.validate()?;
        positive(density, "density")?;
        positive(pressure, "pressure")?;
        Ok((self.gamma * pressure / density).sqrt())
    }
    pub fn temperature(self, density: f64, pressure: f64) -> Result<f64> {
        self.validate()?;
        positive(density, "density")?;
        positive(pressure, "pressure")?;
        Ok(pressure / (density * self.gas_constant))
    }
}
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
pub struct EulerState1D {
    pub density: f64,
    pub momentum: f64,
    pub energy: f64,
}
impl EulerState1D {
    fn vector(self) -> Vector3<f64> {
        Vector3::new(self.density, self.momentum, self.energy)
    }
    fn from_vector(v: Vector3<f64>) -> Self {
        Self {
            density: v[0],
            momentum: v[1],
            energy: v[2],
        }
    }
    fn add(self, b: Self, scale: f64) -> Self {
        Self {
            density: self.density + scale * b.density,
            momentum: self.momentum + scale * b.momentum,
            energy: self.energy + scale * b.energy,
        }
    }
    fn scale(self, v: f64) -> Self {
        Self {
            density: self.density * v,
            momentum: self.momentum * v,
            energy: self.energy * v,
        }
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct EulerPrimitive1D {
    pub density: f64,
    pub velocity: f64,
    pub pressure: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CompressibleReconstruction {
    MusclMinmod,
    CharacteristicWeno5,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CompressibleTimeIntegrator {
    ForwardEuler,
    Ssprk3,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Compressible1DConfig {
    pub cells: usize,
    pub length: f64,
    pub cfl: f64,
    pub eos: IdealGasEquationOfState,
    pub periodic: bool,
    pub reconstruction: CompressibleReconstruction,
    pub time_integrator: CompressibleTimeIntegrator,
    pub weno_epsilon: f64,
    pub density_floor: f64,
    pub pressure_floor: f64,
}
impl Default for Compressible1DConfig {
    fn default() -> Self {
        Self {
            cells: 200,
            length: 1.0,
            cfl: 0.45,
            eos: Default::default(),
            periodic: false,
            reconstruction: CompressibleReconstruction::MusclMinmod,
            time_integrator: CompressibleTimeIntegrator::ForwardEuler,
            weno_epsilon: 1e-6,
            density_floor: 1e-12,
            pressure_floor: 1e-12,
        }
    }
}
pub fn weno5_left(q: [f64; 5], epsilon: f64) -> Result<f64> {
    positive(epsilon, "WENO epsilon")?;
    crate::finite_slice(&q, "WENO stencil")?;
    let p0 = (2.0 * q[0] - 7.0 * q[1] + 11.0 * q[2]) / 6.0;
    let p1 = (-q[1] + 5.0 * q[2] + 2.0 * q[3]) / 6.0;
    let p2 = (2.0 * q[2] + 5.0 * q[3] - q[4]) / 6.0;
    let b0 = (13.0 / 12.0) * (q[0] - 2.0 * q[1] + q[2]).powi(2)
        + 0.25 * (q[0] - 4.0 * q[1] + 3.0 * q[2]).powi(2);
    let b1 = (13.0 / 12.0) * (q[1] - 2.0 * q[2] + q[3]).powi(2) + 0.25 * (q[1] - q[3]).powi(2);
    let b2 = (13.0 / 12.0) * (q[2] - 2.0 * q[3] + q[4]).powi(2)
        + 0.25 * (3.0 * q[2] - 4.0 * q[3] + q[4]).powi(2);
    let a0 = 0.1 / (epsilon + b0).powi(2);
    let a1 = 0.6 / (epsilon + b1).powi(2);
    let a2 = 0.3 / (epsilon + b2).powi(2);
    Ok((a0 * p0 + a1 * p1 + a2 * p2) / (a0 + a1 + a2))
}
#[derive(Clone, Debug)]
pub struct CompressibleEuler1D {
    config: Compressible1DConfig,
    state: Vec<EulerState1D>,
    pub time: f64,
}
impl CompressibleEuler1D {
    pub fn new(c: Compressible1DConfig) -> Result<Self> {
        require(
            c.cells >= 8 && c.cells <= i64::MAX as usize && c.cfl > 0.0 && c.cfl < 1.0,
            "invalid compressible cells/CFL",
        )?;
        positive(c.length, "length")?;
        c.eos.validate()?;
        positive(c.weno_epsilon, "WENO epsilon")?;
        positive(c.density_floor, "density floor")?;
        positive(c.pressure_floor, "pressure floor")?;
        let mut s = Self {
            config: c,
            state: vec![EulerState1D::default(); c.cells],
            time: 0.0,
        };
        s.initialize_uniform(EulerPrimitive1D {
            density: 1.0,
            velocity: 0.0,
            pressure: 1.0,
        })?;
        Ok(s)
    }
    pub fn config(&self) -> Compressible1DConfig {
        self.config
    }
    pub fn conserved(&self) -> &[EulerState1D] {
        &self.state
    }
    fn conservative(&self, q: EulerPrimitive1D) -> Result<EulerState1D> {
        positive(q.density, "density")?;
        positive(q.pressure, "pressure")?;
        require(q.velocity.is_finite(), "nonfinite velocity")?;
        Ok(EulerState1D {
            density: q.density,
            momentum: q.density * q.velocity,
            energy: q.pressure / (self.config.eos.gamma - 1.0)
                + 0.5 * q.density * q.velocity * q.velocity,
        })
    }
    pub fn initialize(
        &mut self,
        left: EulerPrimitive1D,
        right: EulerPrimitive1D,
        fraction: f64,
    ) -> Result<()> {
        require(
            fraction > 0.0 && fraction < 1.0,
            "invalid Riemann interface",
        )?;
        let l = self.conservative(left)?;
        let r = self.conservative(right)?;
        let n = self.state.len();
        for i in 0..n {
            self.state[i] = if (i as f64 + 0.5) / (n as f64) < fraction {
                l
            } else {
                r
            };
        }
        self.time = 0.0;
        Ok(())
    }
    pub fn initialize_uniform(&mut self, q: EulerPrimitive1D) -> Result<()> {
        let u = self.conservative(q)?;
        self.state.fill(u);
        self.time = 0.0;
        Ok(())
    }
    pub fn initialize_profile(&mut self, q: &[EulerPrimitive1D]) -> Result<()> {
        require(q.len() == self.state.len(), "compressible profile size")?;
        let state = q
            .iter()
            .map(|&x| self.conservative(x))
            .collect::<Result<Vec<_>>>()?;
        self.state = state;
        self.time = 0.0;
        Ok(())
    }
    pub fn primitive(&self, i: usize) -> Result<EulerPrimitive1D> {
        require(i < self.state.len(), "compressible cell index")?;
        let u = self.state[i];
        Ok(EulerPrimitive1D {
            density: u.density,
            velocity: u.momentum / u.density,
            pressure: self.config.eos.pressure(u.density, u.momentum, u.energy)?,
        })
    }
    fn physical(&self, u: EulerState1D) -> bool {
        let p =
            (self.config.eos.gamma - 1.0) * (u.energy - 0.5 * u.momentum * u.momentum / u.density);
        u.density.is_finite()
            && u.momentum.is_finite()
            && u.energy.is_finite()
            && u.density >= self.config.density_floor
            && p.is_finite()
            && p >= self.config.pressure_floor
    }
    fn require_physical(&self, u: &[EulerState1D]) -> Result<()> {
        if u.iter().all(|&s| self.physical(s)) {
            Ok(())
        } else {
            Err(SolverError::NumericalFailure(
                "compressible update produced nonphysical state".into(),
            ))
        }
    }
    fn flux(&self, u: EulerState1D) -> Result<EulerState1D> {
        let p = self.config.eos.pressure(u.density, u.momentum, u.energy)?;
        let v = u.momentum / u.density;
        Ok(EulerState1D {
            density: u.momentum,
            momentum: u.momentum * v + p,
            energy: (u.energy + p) * v,
        })
    }
    fn rusanov(&self, l: EulerState1D, r: EulerState1D) -> Result<EulerState1D> {
        let e = self.config.eos;
        let pl = e.pressure(l.density, l.momentum, l.energy)?;
        let pr = e.pressure(r.density, r.momentum, r.energy)?;
        let speed = (l.momentum / l.density).abs() + e.sound_speed(l.density, pl)?;
        let speed = speed.max((r.momentum / r.density).abs() + e.sound_speed(r.density, pr)?);
        Ok(self
            .flux(l)?
            .add(self.flux(r)?, 1.0)
            .scale(0.5)
            .add(r.add(l, -1.0), -0.5 * speed))
    }
    fn roe_basis(&self, l: EulerState1D, r: EulerState1D) -> Result<(Matrix3<f64>, Matrix3<f64>)> {
        let e = self.config.eos;
        let pl = e.pressure(l.density, l.momentum, l.energy)?;
        let pr = e.pressure(r.density, r.momentum, r.energy)?;
        let (rl, rr) = (l.density.sqrt(), r.density.sqrt());
        let v = (rl * l.momentum / l.density + rr * r.momentum / r.density) / (rl + rr);
        let h = (rl * (l.energy + pl) / l.density + rr * (r.energy + pr) / r.density) / (rl + rr);
        let c2 = (e.gamma - 1.0) * (h - 0.5 * v * v);
        positive(c2, "Roe sound speed squared")?;
        let c = c2.sqrt();
        let right = Matrix3::new(
            1.0,
            1.0,
            1.0,
            v - c,
            v,
            v + c,
            h - v * c,
            0.5 * v * v,
            h + v * c,
        );
        require(
            right.determinant().abs() > 1e-14,
            "singular characteristic basis",
        )?;
        let left = right
            .try_inverse()
            .ok_or_else(|| SolverError::Singular("Euler characteristic basis".into()))?;
        Ok((right, left))
    }
    pub fn stable_timestep(&self) -> Result<f64> {
        let mut speed: f64 = 0.0;
        for i in 0..self.state.len() {
            let q = self.primitive(i)?;
            speed =
                speed.max(q.velocity.abs() + self.config.eos.sound_speed(q.density, q.pressure)?);
        }
        positive(speed, "characteristic speed")?;
        Ok(self.config.cfl * (self.config.length / self.state.len() as f64) / speed)
    }
    fn spatial_residual(&self, state: &[EulerState1D]) -> Result<Vec<EulerState1D>> {
        let n = state.len();
        let mut face = vec![EulerState1D::default(); n + 1];
        let sample = |i: i64| {
            state[if self.config.periodic {
                i.rem_euclid(n as i64) as usize
            } else {
                i.clamp(0, n as i64 - 1) as usize
            }]
        };
        if self.config.reconstruction == CompressibleReconstruction::MusclMinmod {
            let mut slope = vec![EulerState1D::default(); n];
            for i in 0..n {
                if !self.config.periodic && (i == 0 || i + 1 == n) {
                    continue;
                }
                let a = sample(i as i64 - 1).vector();
                let b = state[i].vector();
                let c = sample(i as i64 + 1).vector();
                slope[i] = EulerState1D::from_vector(Vector3::from_fn(|k, _| {
                    minmod(b[k] - a[k], c[k] - b[k])
                }));
            }
            for f in 0..=n {
                let li = if f == 0 {
                    if self.config.periodic {
                        n - 1
                    } else {
                        0
                    }
                } else {
                    f - 1
                };
                let ri = if f == n {
                    if self.config.periodic {
                        0
                    } else {
                        n - 1
                    }
                } else {
                    f
                };
                let mut l = state[li].add(slope[li], 0.5);
                let mut r = state[ri].add(slope[ri], -0.5);
                if !self.physical(l) {
                    l = state[li];
                }
                if !self.physical(r) {
                    r = state[ri];
                }
                face[f] = self.rusanov(l, r)?;
            }
        } else {
            for f in 0..=n {
                let ff = f as i64;
                let lc = sample(ff - 1);
                let rc = sample(ff);
                let (mut l, mut r) = (lc, rc);
                if let Ok((right, left)) = self.roe_basis(lc, rc) {
                    let wl: [Vector3<f64>; 5] =
                        std::array::from_fn(|s| left * sample(ff - 3 + s as i64).vector());
                    let wr: [Vector3<f64>; 5] =
                        std::array::from_fn(|s| left * sample(ff + 2 - s as i64).vector());
                    let mut a = Vector3::zeros();
                    let mut b = Vector3::zeros();
                    for k in 0..3 {
                        a[k] = weno5_left(
                            std::array::from_fn(|s| wl[s][k]),
                            self.config.weno_epsilon,
                        )?;
                        b[k] = weno5_left(
                            std::array::from_fn(|s| wr[s][k]),
                            self.config.weno_epsilon,
                        )?;
                    }
                    l = EulerState1D::from_vector(right * a);
                    r = EulerState1D::from_vector(right * b);
                }
                if !self.physical(l) {
                    l = lc;
                }
                if !self.physical(r) {
                    r = rc;
                }
                face[f] = self.rusanov(l, r)?;
            }
        }
        let dx = self.config.length / n as f64;
        Ok((0..n)
            .map(|i| face[i + 1].add(face[i], -1.0).scale(-1.0 / dx))
            .collect())
    }
    fn advance(&self, base: &[EulerState1D], dt: f64) -> Result<Vec<EulerState1D>> {
        let rhs = self.spatial_residual(base)?;
        let u: Vec<_> = base.iter().zip(rhs).map(|(&b, r)| b.add(r, dt)).collect();
        self.require_physical(&u)?;
        Ok(u)
    }
    pub fn step(&mut self, dt: f64) -> Result<f64> {
        let dt = if dt == 0.0 {
            self.stable_timestep()?
        } else {
            dt
        };
        positive(dt, "compressible dt")?;
        let next = if self.config.time_integrator == CompressibleTimeIntegrator::ForwardEuler {
            self.advance(&self.state, dt)?
        } else {
            let u0 = &self.state;
            let u1 = self.advance(u0, dt)?;
            let u1e = self.advance(&u1, dt)?;
            let u2: Vec<_> = u0
                .iter()
                .zip(u1e)
                .map(|(&a, b)| a.scale(0.75).add(b, 0.25))
                .collect();
            self.require_physical(&u2)?;
            let u2e = self.advance(&u2, dt)?;
            let next: Vec<_> = u0
                .iter()
                .zip(u2e)
                .map(|(&a, b)| a.scale(1.0 / 3.0).add(b, 2.0 / 3.0))
                .collect();
            self.require_physical(&next)?;
            next
        };
        self.state = next;
        self.time += dt;
        Ok(dt)
    }
    pub fn run(&mut self, duration: f64) -> Result<()> {
        require(duration.is_finite() && duration >= 0.0, "invalid duration")?;
        let mut t = 0.0;
        while t < duration {
            let dt = self.stable_timestep()?.min(duration - t);
            require(t + dt > t, "timestep below time resolution")?;
            t += self.step(dt)?;
        }
        Ok(())
    }
    pub fn pressure_jump_sensor(&self) -> Result<Vec<f64>> {
        let n = self.state.len();
        let p = (0..n)
            .map(|i| Ok(self.primitive(i)?.pressure))
            .collect::<Result<Vec<_>>>()?;
        let mut s = vec![0.0_f64; n];
        for i in 0..n {
            for offset in [-1, 1] {
                let j = i as i64 + offset;
                if !self.config.periodic && (j < 0 || j >= n as i64) {
                    continue;
                }
                let j = j.rem_euclid(n as i64) as usize;
                s[i] = s[i].max((p[j] - p[i]).abs() / (p[i] + p[j] + self.config.pressure_floor));
            }
        }
        Ok(s)
    }
    pub fn mass(&self) -> f64 {
        let dx = self.config.length / self.state.len() as f64;
        self.state.iter().map(|s| s.density * dx).sum()
    }
    pub fn momentum(&self) -> f64 {
        let dx = self.config.length / self.state.len() as f64;
        self.state.iter().map(|s| s.momentum * dx).sum()
    }
    pub fn total_energy(&self) -> f64 {
        let dx = self.config.length / self.state.len() as f64;
        self.state.iter().map(|s| s.energy * dx).sum()
    }
}
