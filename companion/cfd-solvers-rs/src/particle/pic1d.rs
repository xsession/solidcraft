//! 1-D periodic particle-in-cell: spectral Poisson, electrostatic PIC, charge-conserving
//! current reconstruction and the 1-D/3-velocity electromagnetic PIC with Yee-staggered
//! transverse fields and Vay pushing. Ported from upstream `cfd/particle/electromagnetic`
//! (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::pusher::*;
use crate::{require, Result};
use num_complex::Complex64 as C;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PicParticle1D {
    pub position_m: f64,
    pub velocity_m_per_s: f64,
    pub charge_c: f64,
    pub mass_kg: f64,
    pub weight: f64,
}

impl Default for PicParticle1D {
    fn default() -> Self {
        Self {
            position_m: 0.0,
            velocity_m_per_s: 0.0,
            charge_c: 1.0,
            mass_kg: 1.0,
            weight: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrostaticPic1DConfig {
    pub grid_points: usize,
    pub length_m: f64,
    pub relative_permittivity: f64,
    pub dt_s: f64,
    pub neutralize_mean_charge: bool,
}

impl Default for ElectrostaticPic1DConfig {
    fn default() -> Self {
        Self {
            grid_points: 64,
            length_m: 1.0,
            relative_permittivity: 1.0,
            dt_s: 1.0e-3,
            neutralize_mean_charge: true,
        }
    }
}

pub(crate) fn signed_mode(mode: usize, n: usize) -> f64 {
    if mode <= n / 2 {
        mode as f64
    } else {
        mode as f64 - n as f64
    }
}

pub(crate) fn phase(angle: f64) -> C {
    C::new(angle.cos(), angle.sin())
}

/// Periodic spectral Poisson field `E = −∇φ` from nodal charge density (unit cross-section).
/// With `remove_mean = false` a non-neutral density is an error.
pub fn periodic_electric_field_from_charge_density(
    rho: &[f64],
    length: f64,
    relative_permittivity: f64,
    remove_mean: bool,
) -> Result<Vec<f64>> {
    let n = rho.len();
    require(
        n >= 2 && length > 0.0 && relative_permittivity > 0.0,
        "invalid periodic Poisson grid",
    )?;
    let mut mean = 0.0;
    for &v in rho {
        require(v.is_finite(), "non-finite charge density")?;
        mean += v;
    }
    mean /= n as f64;
    if !remove_mean {
        let largest = rho
            .iter()
            .fold(0.0f64, |m, v| if m.abs() < v.abs() { *v } else { m });
        require(
            mean.abs() <= 1.0e-12 * largest.max(1.0),
            "periodic Poisson requires mean-zero charge density",
        )?;
    }
    let eps = EPSILON0 * relative_permittivity;
    let shift = if remove_mean { mean } else { 0.0 };
    let mut field_hat = vec![C::default(); n];
    for mode in 1..n {
        let k = 2.0 * PI * signed_mode(mode, n) / length;
        let mut rho_hat = C::default();
        for j in 0..n {
            let angle = -2.0 * PI * (mode * j) as f64 / n as f64;
            rho_hat += (rho[j] - shift) * phase(angle);
        }
        rho_hat /= n as f64;
        field_hat[mode] = C::new(0.0, -1.0) * rho_hat / (k * eps);
    }
    Ok((0..n)
        .map(|j| {
            let mut v = C::default();
            for mode in 1..n {
                v += field_hat[mode] * phase(2.0 * PI * (mode * j) as f64 / n as f64);
            }
            v.re
        })
        .collect())
}

fn cic(u: f64, n: usize) -> (usize, usize, f64) {
    let fl = u.floor();
    let left = (fl as usize) % n;
    (left, (left + 1) % n, u - fl)
}

/// Periodic electrostatic PIC with CIC deposition and leapfrog pushing.
#[derive(Clone, Debug)]
pub struct ElectrostaticPic1D {
    config: ElectrostaticPic1DConfig,
    particles: Vec<PicParticle1D>,
    charge_density: Vec<f64>,
    electric_field: Vec<f64>,
    time_s: f64,
}

impl ElectrostaticPic1D {
    pub fn new(config: ElectrostaticPic1DConfig) -> Result<Self> {
        require(
            config.grid_points >= 2
                && config.length_m > 0.0
                && config.relative_permittivity > 0.0
                && config.dt_s > 0.0,
            "invalid electrostatic PIC configuration",
        )?;
        Ok(Self {
            charge_density: vec![0.0; config.grid_points],
            electric_field: vec![0.0; config.grid_points],
            config,
            particles: Vec::new(),
            time_s: 0.0,
        })
    }
    fn wrap(&self, x: f64) -> f64 {
        wrap_periodic(x, self.config.length_m)
    }
    pub fn set_particles(&mut self, mut particles: Vec<PicParticle1D>) -> Result<()> {
        for p in &mut particles {
            require(
                p.mass_kg > 0.0 && p.charge_c.is_finite() && p.weight.is_finite(),
                "invalid PIC particle",
            )?;
            p.position_m = self.wrap(p.position_m);
        }
        self.particles = particles;
        Ok(())
    }
    pub fn particles(&self) -> &[PicParticle1D] {
        &self.particles
    }
    pub fn charge_density(&self) -> &[f64] {
        &self.charge_density
    }
    pub fn electric_field(&self) -> &[f64] {
        &self.electric_field
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn deposit_and_solve(&mut self) -> Result<()> {
        let n = self.config.grid_points;
        self.charge_density.iter_mut().for_each(|v| *v = 0.0);
        let dx = self.config.length_m / n as f64;
        for p in &self.particles {
            let (l, r, f) = cic(self.wrap(p.position_m) / dx, n);
            let density = p.charge_c * p.weight / dx;
            self.charge_density[l] += density * (1.0 - f);
            self.charge_density[r] += density * f;
        }
        self.electric_field = periodic_electric_field_from_charge_density(
            &self.charge_density,
            self.config.length_m,
            self.config.relative_permittivity,
            self.config.neutralize_mean_charge,
        )?;
        Ok(())
    }
    fn interpolate_field(&self, x: f64) -> f64 {
        let n = self.config.grid_points;
        let dx = self.config.length_m / n as f64;
        let (l, r, f) = cic(self.wrap(x) / dx, n);
        self.electric_field[l] * (1.0 - f) + self.electric_field[r] * f
    }
    pub fn step(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.deposit_and_solve()?;
            for i in 0..self.particles.len() {
                let p = self.particles[i];
                let v = p.velocity_m_per_s
                    + (p.charge_c / p.mass_kg)
                        * self.interpolate_field(p.position_m)
                        * self.config.dt_s;
                let x = self.wrap(p.position_m + v * self.config.dt_s);
                self.particles[i].velocity_m_per_s = v;
                self.particles[i].position_m = x;
            }
            self.time_s += self.config.dt_s;
        }
        self.deposit_and_solve()
    }
}

fn check_charged(p: &ChargedParticle, what: &str) -> Result<()> {
    require(
        p.mass_kg > 0.0 && p.charge_c.is_finite() && p.weight.is_finite(),
        what,
    )
}

pub(crate) fn deposit_cic_charge_density_1d(
    particles: &[ChargedParticle],
    n: usize,
    length: f64,
) -> Result<Vec<f64>> {
    require(n >= 2 && length > 0.0, "invalid CIC charge grid")?;
    let mut rho = vec![0.0; n];
    let cell = length / n as f64;
    for p in particles {
        check_charged(p, "invalid charged particle for deposition")?;
        let (l, r, f) = cic(wrap_periodic(p.position_m.x, length) / cell, n);
        let density = p.charge_c * p.weight / cell;
        rho[l] += density * (1.0 - f);
        rho[r] += density * f;
    }
    Ok(rho)
}

pub(crate) fn fourier_coefficient(values: &[f64], mode: usize) -> C {
    let n = values.len();
    let mut c = C::default();
    for (i, &v) in values.iter().enumerate() {
        c += v * phase(-2.0 * PI * (mode * i) as f64 / n as f64);
    }
    c / n as f64
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChargeConservingCurrent1D {
    pub charge_density_old_c_per_m3: Vec<f64>,
    pub charge_density_new_c_per_m3: Vec<f64>,
    pub current_x_a_per_m2: Vec<f64>,
    pub continuity_linf_residual: f64,
}

/// CIC charge before/after a move and the longitudinal current whose spectral divergence
/// satisfies discrete continuity.
pub fn deposit_charge_conserving_current_1d(
    old: &[ChargedParticle],
    new: &[ChargedParticle],
    n: usize,
    length: f64,
    dt: f64,
) -> Result<ChargeConservingCurrent1D> {
    require(
        old.len() == new.len(),
        "old/new particle arrays must have equal size",
    )?;
    require(
        n >= 2 && length > 0.0 && dt > 0.0 && dt.is_finite(),
        "invalid charge-conserving deposition parameters",
    )?;
    let mut r = ChargeConservingCurrent1D {
        charge_density_old_c_per_m3: deposit_cic_charge_density_1d(old, n, length)?,
        charge_density_new_c_per_m3: deposit_cic_charge_density_1d(new, n, length)?,
        current_x_a_per_m2: vec![0.0; n],
        continuity_linf_residual: 0.0,
    };
    let mut current_hat = vec![C::default(); n];
    let mut divergence_hat = vec![C::default(); n];
    for mode in 1..n {
        let k = 2.0 * PI * signed_mode(mode, n) / length;
        let delta = fourier_coefficient(&r.charge_density_new_c_per_m3, mode)
            - fourier_coefficient(&r.charge_density_old_c_per_m3, mode);
        current_hat[mode] = C::new(0.0, 1.0) * delta / (k * dt);
        divergence_hat[mode] = C::new(0.0, 1.0) * k * current_hat[mode];
    }
    let mut divergence = vec![0.0; n];
    for i in 0..n {
        let (mut cv, mut dv) = (C::default(), C::default());
        for mode in 1..n {
            let ph = phase(2.0 * PI * (mode * i) as f64 / n as f64);
            cv += current_hat[mode] * ph;
            dv += divergence_hat[mode] * ph;
        }
        r.current_x_a_per_m2[i] = cv.re;
        divergence[i] = dv.re;
    }
    for i in 0..n {
        let c = (r.charge_density_new_c_per_m3[i] - r.charge_density_old_c_per_m3[i]) / dt
            + divergence[i];
        r.continuity_linf_residual = r.continuity_linf_residual.max(c.abs());
    }
    Ok(r)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectromagneticPic1DConfig {
    pub grid_points: usize,
    pub length_m: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub dt_s: f64,
    pub solve_longitudinal_poisson: bool,
    pub neutralize_mean_charge: bool,
}

impl Default for ElectromagneticPic1DConfig {
    fn default() -> Self {
        Self {
            grid_points: 64,
            length_m: 1.0,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            dt_s: 1.0e-12,
            solve_longitudinal_poisson: true,
            neutralize_mean_charge: true,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ElectromagneticPicDiagnostics {
    pub field_energy_j: f64,
    pub particle_kinetic_energy_j: f64,
    pub total_energy_j: f64,
    pub charge_continuity_linf_residual: f64,
    pub collision_statistics: CollisionStatistics,
    pub particle_count: usize,
}

/// 1-D/3-velocity periodic electromagnetic PIC: Yee-staggered transverse fields (Ey/Ez at
/// nodes, By/Bz at half cells), spectral Ex, Vay pushing and continuity-consistent Jx.
#[derive(Clone, Debug)]
pub struct ElectromagneticPic1D {
    config: ElectromagneticPic1DConfig,
    particles: Vec<ChargedParticle>,
    charge_density: Vec<f64>,
    current_x: Vec<f64>,
    current_y: Vec<f64>,
    current_z: Vec<f64>,
    electric_x: Vec<f64>,
    electric_y: Vec<f64>,
    electric_z: Vec<f64>,
    magnetic_y: Vec<f64>,
    magnetic_z: Vec<f64>,
    time_s: f64,
    last_continuity_residual: f64,
    accumulated_collisions: CollisionStatistics,
}

impl ElectromagneticPic1D {
    pub fn new(config: ElectromagneticPic1DConfig) -> Result<Self> {
        require(
            config.grid_points >= 2
                && config.length_m > 0.0
                && config.relative_permittivity > 0.0
                && config.relative_permeability > 0.0
                && config.dt_s > 0.0,
            "invalid electromagnetic PIC configuration",
        )?;
        let n = config.grid_points;
        let s = Self {
            config,
            particles: Vec::new(),
            charge_density: vec![0.0; n],
            current_x: vec![0.0; n],
            current_y: vec![0.0; n],
            current_z: vec![0.0; n],
            electric_x: vec![0.0; n],
            electric_y: vec![0.0; n],
            electric_z: vec![0.0; n],
            magnetic_y: vec![0.0; n],
            magnetic_z: vec![0.0; n],
            time_s: 0.0,
            last_continuity_residual: 0.0,
            accumulated_collisions: CollisionStatistics::default(),
        };
        require(
            s.config.dt_s * s.medium_light_speed() / s.dx() <= 1.0,
            "1-D electromagnetic PIC CFL must be <= 1",
        )?;
        Ok(s)
    }
    fn dx(&self) -> f64 {
        self.config.length_m / self.config.grid_points as f64
    }
    fn epsilon(&self) -> f64 {
        EPSILON0 * self.config.relative_permittivity
    }
    fn mu(&self) -> f64 {
        (4.0e-7 * PI) * self.config.relative_permeability
    }
    fn medium_light_speed(&self) -> f64 {
        1.0 / (self.epsilon() * self.mu()).sqrt()
    }
    fn wrap(&self, x: f64) -> f64 {
        wrap_periodic(x, self.config.length_m)
    }
    pub fn particles(&self) -> &[ChargedParticle] {
        &self.particles
    }
    pub fn charge_density(&self) -> &[f64] {
        &self.charge_density
    }
    pub fn current_x(&self) -> &[f64] {
        &self.current_x
    }
    pub fn current_y(&self) -> &[f64] {
        &self.current_y
    }
    pub fn current_z(&self) -> &[f64] {
        &self.current_z
    }
    pub fn electric_x(&self) -> &[f64] {
        &self.electric_x
    }
    pub fn electric_y(&self) -> &[f64] {
        &self.electric_y
    }
    pub fn electric_z(&self) -> &[f64] {
        &self.electric_z
    }
    pub fn magnetic_y(&self) -> &[f64] {
        &self.magnetic_y
    }
    pub fn magnetic_z(&self) -> &[f64] {
        &self.magnetic_z
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn set_particles(&mut self, mut particles: Vec<ChargedParticle>) -> Result<()> {
        for p in &mut particles {
            check_charged(p, "invalid electromagnetic PIC particle")?;
            p.position_m.x = self.wrap(p.position_m.x);
        }
        self.particles = particles;
        self.deposit_sources()
    }
    pub fn set_transverse_fields(
        &mut self,
        ey: &[f64],
        ez: &[f64],
        by: &[f64],
        bz: &[f64],
    ) -> Result<()> {
        let n = self.config.grid_points;
        require(
            ey.len() == n && ez.len() == n && by.len() == n && bz.len() == n,
            "transverse PIC field arrays must match grid size",
        )?;
        self.electric_y = ey.to_vec();
        self.electric_z = ez.to_vec();
        self.magnetic_y = by.to_vec();
        self.magnetic_z = bz.to_vec();
        Ok(())
    }
    pub fn initialize_right_traveling_mode(
        &mut self,
        amplitude: f64,
        mode_index: usize,
    ) -> Result<()> {
        require(mode_index != 0, "mode index must be non-zero")?;
        let cell = self.dx();
        let k = 2.0 * PI * mode_index as f64 / self.config.length_m;
        let c = self.medium_light_speed();
        for i in 0..self.config.grid_points {
            let xn = i as f64 * cell;
            let xh = (i as f64 + 0.5) * cell;
            self.electric_y[i] = amplitude * (k * xn).sin();
            self.electric_z[i] = 0.0;
            self.magnetic_y[i] = 0.0;
            self.magnetic_z[i] = (amplitude / c) * (k * xh).sin();
        }
        Ok(())
    }
    fn solve_ex(&mut self) -> Result<()> {
        if self.config.solve_longitudinal_poisson {
            self.electric_x = periodic_electric_field_from_charge_density(
                &self.charge_density,
                self.config.length_m,
                self.config.relative_permittivity,
                self.config.neutralize_mean_charge,
            )?;
        } else {
            self.electric_x.iter_mut().for_each(|v| *v = 0.0);
        }
        Ok(())
    }
    pub fn deposit_sources(&mut self) -> Result<()> {
        self.charge_density = deposit_cic_charge_density_1d(
            &self.particles,
            self.config.grid_points,
            self.config.length_m,
        )?;
        for v in [
            &mut self.current_x,
            &mut self.current_y,
            &mut self.current_z,
        ] {
            v.iter_mut().for_each(|x| *x = 0.0);
        }
        self.solve_ex()
    }
    fn node_field(&self, nodal: &[f64], x: f64) -> f64 {
        let (l, r, f) = cic(self.wrap(x) / self.dx(), self.config.grid_points);
        nodal[l] * (1.0 - f) + nodal[r] * f
    }
    fn half_field(&self, half: &[f64], x: f64) -> f64 {
        let cell = self.dx();
        let (l, r, f) = cic(self.wrap(x - 0.5 * cell) / cell, self.config.grid_points);
        half[l] * (1.0 - f) + half[r] * f
    }
    pub fn gather_field(&self, x: f64) -> ElectromagneticField {
        ElectromagneticField {
            electric_v_per_m: Vec3::new(
                self.node_field(&self.electric_x, x),
                self.node_field(&self.electric_y, x),
                self.node_field(&self.electric_z, x),
            ),
            magnetic_t: Vec3::new(
                0.0,
                self.half_field(&self.magnetic_y, x),
                self.half_field(&self.magnetic_z, x),
            ),
        }
    }
    fn advance_magnetic(&mut self, dt: f64) {
        let cell = self.dx();
        let n = self.config.grid_points;
        for i in 0..n {
            let r = (i + 1) % n;
            self.magnetic_y[i] += dt * (self.electric_z[r] - self.electric_z[i]) / cell;
            self.magnetic_z[i] -= dt * (self.electric_y[r] - self.electric_y[i]) / cell;
        }
    }
    fn deposit_transverse_current(&mut self, old: &[ChargedParticle]) {
        self.current_y.iter_mut().for_each(|v| *v = 0.0);
        self.current_z.iter_mut().for_each(|v| *v = 0.0);
        let cell = self.dx();
        let len = self.config.length_m;
        for (po, pn) in old.iter().zip(&self.particles) {
            let (mut xo, mut xn) = (po.position_m.x, pn.position_m.x);
            let delta = xn - xo;
            if delta > 0.5 * len {
                xo += len;
            } else if delta < -0.5 * len {
                xn += len;
            }
            let xm = wrap_periodic(0.5 * (xo + xn), len);
            let (l, r, f) = cic(xm / cell, self.config.grid_points);
            let s = pn.charge_c * pn.weight / cell;
            let vy = 0.5 * (po.velocity_m_per_s.y + pn.velocity_m_per_s.y);
            let vz = 0.5 * (po.velocity_m_per_s.z + pn.velocity_m_per_s.z);
            self.current_y[l] += s * vy * (1.0 - f);
            self.current_y[r] += s * vy * f;
            self.current_z[l] += s * vz * (1.0 - f);
            self.current_z[r] += s * vz * f;
        }
    }
    /// Advance `steps` cycles; with a collision model, collisions follow each cycle.
    pub fn step(&mut self, steps: usize, collisions: Option<&NeutralCollisionModel>) -> Result<()> {
        for _ in 0..steps {
            let old = self.particles.clone();
            self.advance_magnetic(0.5 * self.config.dt_s);
            for i in 0..self.particles.len() {
                let f = self.gather_field(self.particles[i].position_m.x);
                vay_push(&mut self.particles[i], &f, self.config.dt_s)?;
                self.particles[i].position_m.x = self.wrap(self.particles[i].position_m.x);
            }
            let cc = deposit_charge_conserving_current_1d(
                &old,
                &self.particles,
                self.config.grid_points,
                self.config.length_m,
                self.config.dt_s,
            )?;
            self.charge_density = cc.charge_density_new_c_per_m3;
            self.current_x = cc.current_x_a_per_m2;
            self.last_continuity_residual = cc.continuity_linf_residual;
            self.deposit_transverse_current(&old);
            let (cell, eps, mu, n, dt) = (
                self.dx(),
                self.epsilon(),
                self.mu(),
                self.config.grid_points,
                self.config.dt_s,
            );
            for i in 0..n {
                let l = (i + n - 1) % n;
                let curl_y = -(self.magnetic_z[i] - self.magnetic_z[l]) / (mu * cell);
                let curl_z = (self.magnetic_y[i] - self.magnetic_y[l]) / (mu * cell);
                self.electric_y[i] += dt * (curl_y - self.current_y[i]) / eps;
                self.electric_z[i] += dt * (curl_z - self.current_z[i]) / eps;
            }
            self.advance_magnetic(0.5 * dt);
            self.solve_ex()?;
            if let Some(model) = collisions {
                let stats = collide(
                    &mut self.particles,
                    dt,
                    model,
                    "invalid collision particle mass",
                )?;
                self.accumulated_collisions.accumulate(&stats);
                let len = self.config.length_m;
                for p in &mut self.particles {
                    p.position_m.x = wrap_periodic(p.position_m.x, len);
                }
            }
            self.time_s += dt;
        }
        Ok(())
    }
    pub fn diagnostics(&self) -> Result<ElectromagneticPicDiagnostics> {
        let (cell, eps, mu) = (self.dx(), self.epsilon(), self.mu());
        let mut d = ElectromagneticPicDiagnostics::default();
        for i in 0..self.config.grid_points {
            let e2 = self.electric_x[i] * self.electric_x[i]
                + self.electric_y[i] * self.electric_y[i]
                + self.electric_z[i] * self.electric_z[i];
            let b2 =
                self.magnetic_y[i] * self.magnetic_y[i] + self.magnetic_z[i] * self.magnetic_z[i];
            d.field_energy_j += 0.5 * eps * e2 * cell + 0.5 * b2 * cell / mu;
        }
        for p in &self.particles {
            d.particle_kinetic_energy_j += p.weight * kinetic_energy_j(p)?;
        }
        d.total_energy_j = d.field_energy_j + d.particle_kinetic_energy_j;
        d.charge_continuity_linf_residual = self.last_continuity_residual;
        d.collision_statistics = self.accumulated_collisions;
        d.particle_count = self.particles.len();
        Ok(d)
    }
}
