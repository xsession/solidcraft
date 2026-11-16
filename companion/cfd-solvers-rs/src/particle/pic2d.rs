//! 2-D periodic particle-in-cell: CIC deposition, spectral Poisson, charge-conserving currents,
//! electrostatic and collocated electromagnetic PIC, wall boundaries and field boundary
//! conditions. Ported from upstream `cfd/particle/electromagnetic` (MIT, cfd_solvers
//! contributors).
#![allow(clippy::needless_range_loop)]
use super::pic1d::phase;
use super::pusher::*;
use crate::{require, Result};
use num_complex::Complex64 as C;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

impl Vec2 {
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PicParticle2D {
    pub position_m: Vec2,
    pub velocity_m_per_s: Vec3,
    pub charge_c: f64,
    pub mass_kg: f64,
    pub weight: f64,
}

impl Default for PicParticle2D {
    fn default() -> Self {
        Self {
            position_m: Vec2::default(),
            velocity_m_per_s: Vec3::default(),
            charge_c: 1.0,
            mass_kg: 1.0,
            weight: 1.0,
        }
    }
}

impl PicParticle2D {
    pub(crate) fn proxy(&self) -> ChargedParticle {
        ChargedParticle {
            position_m: Vec3::new(self.position_m.x, self.position_m.y, 0.0),
            velocity_m_per_s: self.velocity_m_per_s,
            charge_c: self.charge_c,
            mass_kg: self.mass_kg,
            weight: self.weight,
        }
    }
}

impl CollisionParticle for PicParticle2D {
    fn velocity(&self) -> Vec3 {
        self.velocity_m_per_s
    }
    fn set_velocity(&mut self, v: Vec3) {
        self.velocity_m_per_s = v;
    }
    fn mass(&self) -> f64 {
        self.mass_kg
    }
    fn weight(&self) -> f64 {
        self.weight
    }
    fn as_charged(&self) -> ChargedParticle {
        self.proxy()
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ElectrostaticField2D {
    pub electric_x_v_per_m: Vec<f64>,
    pub electric_y_v_per_m: Vec<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChargeConservingCurrent2D {
    pub charge_density_old_c_per_m3: Vec<f64>,
    pub charge_density_new_c_per_m3: Vec<f64>,
    pub current_x_a_per_m2: Vec<f64>,
    pub current_y_a_per_m2: Vec<f64>,
    pub continuity_linf_residual: f64,
}

pub(crate) fn fi2(ix: usize, iy: usize, nx: usize) -> usize {
    iy * nx + ix
}

fn fourier_2d(values: &[f64], nx: usize, ny: usize, mx: usize, my: usize) -> C {
    let mut c = C::default();
    for iy in 0..ny {
        for ix in 0..nx {
            let angle = -2.0 * PI * ((mx * ix) as f64 / nx as f64 + (my * iy) as f64 / ny as f64);
            c += values[fi2(ix, iy, nx)] * phase(angle);
        }
    }
    c / (nx * ny) as f64
}

fn signed(m: usize, n: usize) -> f64 {
    super::pic1d::signed_mode(m, n)
}

pub(crate) struct Cic2 {
    pub ix0: usize,
    pub iy0: usize,
    pub ix1: usize,
    pub iy1: usize,
    pub fx: f64,
    pub fy: f64,
}

pub(crate) fn cic2(ux: f64, uy: f64, nx: usize, ny: usize) -> Cic2 {
    let (flx, fly) = (ux.floor(), uy.floor());
    let ix0 = (flx as usize) % nx;
    let iy0 = (fly as usize) % ny;
    Cic2 {
        ix0,
        iy0,
        ix1: (ix0 + 1) % nx,
        iy1: (iy0 + 1) % ny,
        fx: ux - flx,
        fy: uy - fly,
    }
}

impl Cic2 {
    pub(crate) fn scatter(&self, grid: &mut [f64], nx: usize, value: f64) {
        grid[fi2(self.ix0, self.iy0, nx)] += value * (1.0 - self.fx) * (1.0 - self.fy);
        grid[fi2(self.ix1, self.iy0, nx)] += value * self.fx * (1.0 - self.fy);
        grid[fi2(self.ix0, self.iy1, nx)] += value * (1.0 - self.fx) * self.fy;
        grid[fi2(self.ix1, self.iy1, nx)] += value * self.fx * self.fy;
    }
    pub(crate) fn gather(&self, grid: &[f64], nx: usize) -> f64 {
        grid[fi2(self.ix0, self.iy0, nx)] * (1.0 - self.fx) * (1.0 - self.fy)
            + grid[fi2(self.ix1, self.iy0, nx)] * self.fx * (1.0 - self.fy)
            + grid[fi2(self.ix0, self.iy1, nx)] * (1.0 - self.fx) * self.fy
            + grid[fi2(self.ix1, self.iy1, nx)] * self.fx * self.fy
    }
}

pub fn deposit_cic_charge_density_2d(
    particles: &[PicParticle2D],
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
) -> Result<Vec<f64>> {
    require(
        nx >= 2 && ny >= 2 && lx > 0.0 && ly > 0.0,
        "invalid 2-D CIC grid",
    )?;
    let mut rho = vec![0.0; nx * ny];
    let (dx, dy) = (lx / nx as f64, ly / ny as f64);
    for p in particles {
        require(
            p.mass_kg > 0.0 && p.charge_c.is_finite() && p.weight.is_finite(),
            "invalid 2-D PIC particle",
        )?;
        let w = cic2(
            wrap_periodic(p.position_m.x, lx) / dx,
            wrap_periodic(p.position_m.y, ly) / dy,
            nx,
            ny,
        );
        w.scatter(&mut rho, nx, p.charge_c * p.weight / (dx * dy));
    }
    Ok(rho)
}

/// Spectral periodic Poisson in 2-D (direct DFT, O(N²)).
pub fn periodic_electric_field_from_charge_density_2d(
    rho: &[f64],
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    relative_permittivity: f64,
    remove_mean: bool,
) -> Result<ElectrostaticField2D> {
    require(
        rho.len() == nx * ny
            && nx >= 2
            && ny >= 2
            && lx > 0.0
            && ly > 0.0
            && relative_permittivity > 0.0,
        "invalid 2-D periodic Poisson grid",
    )?;
    let (mut mean, mut max_abs) = (0.0, 0.0f64);
    for &v in rho {
        require(v.is_finite(), "non-finite 2-D charge density")?;
        mean += v;
        max_abs = max_abs.max(v.abs());
    }
    mean /= rho.len() as f64;
    require(
        remove_mean || mean.abs() <= 1.0e-12 * max_abs.max(1.0),
        "2-D periodic Poisson requires mean-zero charge density",
    )?;
    let centered: Vec<f64> = rho
        .iter()
        .map(|v| if remove_mean { v - mean } else { *v })
        .collect();
    let eps = EPSILON0 * relative_permittivity;
    let mut ex_hat = vec![C::default(); nx * ny];
    let mut ey_hat = vec![C::default(); nx * ny];
    for my in 0..ny {
        for mx in 0..nx {
            if mx == 0 && my == 0 {
                continue;
            }
            let kx = 2.0 * PI * signed(mx, nx) / lx;
            let ky = 2.0 * PI * signed(my, ny) / ly;
            let k2 = kx * kx + ky * ky;
            if k2 == 0.0 {
                continue;
            }
            let factor = C::new(0.0, -1.0) * fourier_2d(&centered, nx, ny, mx, my) / (eps * k2);
            ex_hat[fi2(mx, my, nx)] = factor * kx;
            ey_hat[fi2(mx, my, nx)] = factor * ky;
        }
    }
    let mut field = ElectrostaticField2D {
        electric_x_v_per_m: vec![0.0; nx * ny],
        electric_y_v_per_m: vec![0.0; nx * ny],
    };
    for iy in 0..ny {
        for ix in 0..nx {
            let (mut ex, mut ey) = (C::default(), C::default());
            for my in 0..ny {
                for mx in 0..nx {
                    let ph = phase(
                        2.0 * PI * ((mx * ix) as f64 / nx as f64 + (my * iy) as f64 / ny as f64),
                    );
                    ex += ex_hat[fi2(mx, my, nx)] * ph;
                    ey += ey_hat[fi2(mx, my, nx)] * ph;
                }
            }
            field.electric_x_v_per_m[fi2(ix, iy, nx)] = ex.re;
            field.electric_y_v_per_m[fi2(ix, iy, nx)] = ey.re;
        }
    }
    Ok(field)
}

/// CIC charge before/after a move and the curl-free current satisfying spectral continuity.
pub fn deposit_charge_conserving_current_2d(
    old: &[PicParticle2D],
    new: &[PicParticle2D],
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    dt: f64,
) -> Result<ChargeConservingCurrent2D> {
    require(
        old.len() == new.len(),
        "old/new 2-D particle arrays must have equal size",
    )?;
    require(
        nx >= 2 && ny >= 2 && lx > 0.0 && ly > 0.0 && dt > 0.0 && dt.is_finite(),
        "invalid 2-D charge-conserving deposition parameters",
    )?;
    let mut r = ChargeConservingCurrent2D {
        charge_density_old_c_per_m3: deposit_cic_charge_density_2d(old, nx, ny, lx, ly)?,
        charge_density_new_c_per_m3: deposit_cic_charge_density_2d(new, nx, ny, lx, ly)?,
        current_x_a_per_m2: vec![0.0; nx * ny],
        current_y_a_per_m2: vec![0.0; nx * ny],
        continuity_linf_residual: 0.0,
    };
    let delta: Vec<f64> = r
        .charge_density_new_c_per_m3
        .iter()
        .zip(&r.charge_density_old_c_per_m3)
        .map(|(a, b)| a - b)
        .collect();
    let n = nx * ny;
    let (mut jx_hat, mut jy_hat, mut div_hat) = (
        vec![C::default(); n],
        vec![C::default(); n],
        vec![C::default(); n],
    );
    for my in 0..ny {
        for mx in 0..nx {
            if mx == 0 && my == 0 {
                continue;
            }
            let kx = 2.0 * PI * signed(mx, nx) / lx;
            let ky = 2.0 * PI * signed(my, ny) / ly;
            let k2 = kx * kx + ky * ky;
            if k2 == 0.0 {
                continue;
            }
            let factor = C::new(0.0, 1.0) * fourier_2d(&delta, nx, ny, mx, my) / (dt * k2);
            let k = fi2(mx, my, nx);
            jx_hat[k] = factor * kx;
            jy_hat[k] = factor * ky;
            div_hat[k] = C::new(0.0, 1.0) * (kx * jx_hat[k] + ky * jy_hat[k]);
        }
    }
    let mut divergence = vec![0.0; n];
    for iy in 0..ny {
        for ix in 0..nx {
            let (mut jx, mut jy, mut dv) = (C::default(), C::default(), C::default());
            for my in 0..ny {
                for mx in 0..nx {
                    let ph = phase(
                        2.0 * PI * ((mx * ix) as f64 / nx as f64 + (my * iy) as f64 / ny as f64),
                    );
                    let k = fi2(mx, my, nx);
                    jx += jx_hat[k] * ph;
                    jy += jy_hat[k] * ph;
                    dv += div_hat[k] * ph;
                }
            }
            let i = fi2(ix, iy, nx);
            r.current_x_a_per_m2[i] = jx.re;
            r.current_y_a_per_m2[i] = jy.re;
            divergence[i] = dv.re;
        }
    }
    for i in 0..n {
        r.continuity_linf_residual = r
            .continuity_linf_residual
            .max((delta[i] / dt + divergence[i]).abs());
    }
    Ok(r)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrostaticPic2DConfig {
    pub nx: usize,
    pub ny: usize,
    pub length_x_m: f64,
    pub length_y_m: f64,
    pub relative_permittivity: f64,
    pub dt_s: f64,
    pub neutralize_mean_charge: bool,
}

impl Default for ElectrostaticPic2DConfig {
    fn default() -> Self {
        Self {
            nx: 32,
            ny: 32,
            length_x_m: 1.0,
            length_y_m: 1.0,
            relative_permittivity: 1.0,
            dt_s: 1.0e-3,
            neutralize_mean_charge: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ElectrostaticPic2D {
    config: ElectrostaticPic2DConfig,
    particles: Vec<PicParticle2D>,
    charge_density: Vec<f64>,
    electric_x: Vec<f64>,
    electric_y: Vec<f64>,
    time_s: f64,
}

impl ElectrostaticPic2D {
    pub fn new(config: ElectrostaticPic2DConfig) -> Result<Self> {
        let c = &config;
        require(
            c.nx >= 2
                && c.ny >= 2
                && c.length_x_m > 0.0
                && c.length_y_m > 0.0
                && c.relative_permittivity > 0.0
                && c.dt_s > 0.0,
            "invalid electrostatic PIC 2-D configuration",
        )?;
        let n = c.nx * c.ny;
        Ok(Self {
            config,
            particles: Vec::new(),
            charge_density: vec![0.0; n],
            electric_x: vec![0.0; n],
            electric_y: vec![0.0; n],
            time_s: 0.0,
        })
    }
    fn wrap(&self, p: Vec2) -> Vec2 {
        Vec2::new(
            wrap_periodic(p.x, self.config.length_x_m),
            wrap_periodic(p.y, self.config.length_y_m),
        )
    }
    pub fn set_particles(&mut self, mut particles: Vec<PicParticle2D>) -> Result<()> {
        for p in &mut particles {
            require(
                p.mass_kg > 0.0 && p.charge_c.is_finite() && p.weight.is_finite(),
                "invalid electrostatic PIC 2-D particle",
            )?;
            p.position_m = self.wrap(p.position_m);
        }
        self.particles = particles;
        Ok(())
    }
    pub fn particles(&self) -> &[PicParticle2D] {
        &self.particles
    }
    pub fn charge_density(&self) -> &[f64] {
        &self.charge_density
    }
    pub fn electric_x(&self) -> &[f64] {
        &self.electric_x
    }
    pub fn electric_y(&self) -> &[f64] {
        &self.electric_y
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn deposit_and_solve(&mut self) -> Result<()> {
        let c = self.config;
        self.charge_density =
            deposit_cic_charge_density_2d(&self.particles, c.nx, c.ny, c.length_x_m, c.length_y_m)?;
        let f = periodic_electric_field_from_charge_density_2d(
            &self.charge_density,
            c.nx,
            c.ny,
            c.length_x_m,
            c.length_y_m,
            c.relative_permittivity,
            c.neutralize_mean_charge,
        )?;
        self.electric_x = f.electric_x_v_per_m;
        self.electric_y = f.electric_y_v_per_m;
        Ok(())
    }
    fn bilinear(&self, nodal: &[f64], p: Vec2) -> f64 {
        let c = &self.config;
        let w = self.wrap(p);
        cic2(
            w.x / (c.length_x_m / c.nx as f64),
            w.y / (c.length_y_m / c.ny as f64),
            c.nx,
            c.ny,
        )
        .gather(nodal, c.nx)
    }
    pub fn gather_electric_field(&self, p: Vec2) -> Vec2 {
        Vec2::new(
            self.bilinear(&self.electric_x, p),
            self.bilinear(&self.electric_y, p),
        )
    }
    pub fn step(&mut self, steps: usize) -> Result<()> {
        let dt = self.config.dt_s;
        for _ in 0..steps {
            self.deposit_and_solve()?;
            for i in 0..self.particles.len() {
                let mut p = self.particles[i];
                let f = self.gather_electric_field(p.position_m);
                p.velocity_m_per_s.x += (p.charge_c / p.mass_kg) * f.x * dt;
                p.velocity_m_per_s.y += (p.charge_c / p.mass_kg) * f.y * dt;
                p.position_m.x += p.velocity_m_per_s.x * dt;
                p.position_m.y += p.velocity_m_per_s.y * dt;
                p.position_m = self.wrap(p.position_m);
                self.particles[i] = p;
            }
            self.time_s += dt;
        }
        self.deposit_and_solve()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectromagneticPic2DConfig {
    pub nx: usize,
    pub ny: usize,
    pub length_x_m: f64,
    pub length_y_m: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub dt_s: f64,
    pub solve_longitudinal_poisson: bool,
    pub neutralize_mean_charge: bool,
}

impl Default for ElectromagneticPic2DConfig {
    fn default() -> Self {
        Self {
            nx: 32,
            ny: 32,
            length_x_m: 1.0,
            length_y_m: 1.0,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            dt_s: 1.0e-12,
            solve_longitudinal_poisson: false,
            neutralize_mean_charge: true,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ElectromagneticPic2DDiagnostics {
    pub field_energy_j: f64,
    pub particle_kinetic_energy_j: f64,
    pub total_energy_j: f64,
    pub charge_continuity_linf_residual: f64,
    pub collision_statistics: CollisionStatistics,
    pub particle_count: usize,
}

/// Six-component fields on one node grid with centred differences (collocated 2-D EM PIC).
#[derive(Clone, Debug)]
pub struct ElectromagneticPic2D {
    config: ElectromagneticPic2DConfig,
    particles: Vec<PicParticle2D>,
    charge_density: Vec<f64>,
    current_x: Vec<f64>,
    current_y: Vec<f64>,
    current_z: Vec<f64>,
    electric_x: Vec<f64>,
    electric_y: Vec<f64>,
    electric_z: Vec<f64>,
    magnetic_x: Vec<f64>,
    magnetic_y: Vec<f64>,
    magnetic_z: Vec<f64>,
    time_s: f64,
    last_continuity_residual: f64,
    accumulated_collisions: CollisionStatistics,
}

macro_rules! field_getters {
    ($($name:ident),*) => {
        $(pub fn $name(&self) -> &[f64] {
            &self.$name
        })*
    };
}
pub(crate) use field_getters;

impl ElectromagneticPic2D {
    pub fn new(config: ElectromagneticPic2DConfig) -> Result<Self> {
        let c = &config;
        require(
            c.nx >= 2
                && c.ny >= 2
                && c.length_x_m > 0.0
                && c.length_y_m > 0.0
                && c.relative_permittivity > 0.0
                && c.relative_permeability > 0.0
                && c.dt_s > 0.0
                && c.dt_s.is_finite(),
            "invalid electromagnetic PIC 2-D configuration",
        )?;
        let n = c.nx * c.ny;
        let z = vec![0.0; n];
        Ok(Self {
            config,
            particles: Vec::new(),
            charge_density: z.clone(),
            current_x: z.clone(),
            current_y: z.clone(),
            current_z: z.clone(),
            electric_x: z.clone(),
            electric_y: z.clone(),
            electric_z: z.clone(),
            magnetic_x: z.clone(),
            magnetic_y: z.clone(),
            magnetic_z: z,
            time_s: 0.0,
            last_continuity_residual: 0.0,
            accumulated_collisions: CollisionStatistics::default(),
        })
    }
    field_getters!(
        charge_density,
        current_x,
        current_y,
        current_z,
        electric_x,
        electric_y,
        electric_z,
        magnetic_x,
        magnetic_y,
        magnetic_z
    );
    pub fn particles(&self) -> &[PicParticle2D] {
        &self.particles
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    fn dx(&self) -> f64 {
        self.config.length_x_m / self.config.nx as f64
    }
    fn dy(&self) -> f64 {
        self.config.length_y_m / self.config.ny as f64
    }
    fn epsilon(&self) -> f64 {
        EPSILON0 * self.config.relative_permittivity
    }
    fn mu(&self) -> f64 {
        MU0 * self.config.relative_permeability
    }
    fn light(&self) -> f64 {
        1.0 / (self.epsilon() * self.mu()).sqrt()
    }
    fn wrap(&self, p: Vec2) -> Vec2 {
        Vec2::new(
            wrap_periodic(p.x, self.config.length_x_m),
            wrap_periodic(p.y, self.config.length_y_m),
        )
    }
    pub fn set_particles(&mut self, mut particles: Vec<PicParticle2D>) -> Result<()> {
        for p in &mut particles {
            require(
                p.mass_kg > 0.0
                    && p.mass_kg.is_finite()
                    && p.charge_c.is_finite()
                    && p.weight.is_finite(),
                "invalid electromagnetic PIC 2-D particle",
            )?;
            p.position_m = self.wrap(p.position_m);
        }
        self.particles = particles;
        Ok(())
    }
    pub fn set_fields(
        &mut self,
        ex: &[f64],
        ey: &[f64],
        ez: &[f64],
        bx: &[f64],
        by: &[f64],
        bz: &[f64],
    ) -> Result<()> {
        let n = self.config.nx * self.config.ny;
        require(
            [ex, ey, ez, bx, by, bz].iter().all(|f| f.len() == n),
            "2-D EM PIC field arrays must match grid size",
        )?;
        self.electric_x = ex.to_vec();
        self.electric_y = ey.to_vec();
        self.electric_z = ez.to_vec();
        self.magnetic_x = bx.to_vec();
        self.magnetic_y = by.to_vec();
        self.magnetic_z = bz.to_vec();
        Ok(())
    }
    /// Travelling TM_z plane wave `Ez = A sin(k·x)`, `B = k × E / ω`.
    pub fn initialize_tm_z_mode(
        &mut self,
        amplitude: f64,
        mode_x: usize,
        mode_y: usize,
    ) -> Result<()> {
        require(
            mode_x != 0 || mode_y != 0,
            "2-D EM PIC mode index must be non-zero",
        )?;
        let kx = 2.0 * PI * mode_x as f64 / self.config.length_x_m;
        let ky = 2.0 * PI * mode_y as f64 / self.config.length_y_m;
        let omega = self.light() * kx.hypot(ky);
        let (dx, dy, nx) = (self.dx(), self.dy(), self.config.nx);
        for iy in 0..self.config.ny {
            let y = iy as f64 * dy;
            for ix in 0..nx {
                let x = ix as f64 * dx;
                let i = fi2(ix, iy, nx);
                self.electric_x[i] = 0.0;
                self.electric_y[i] = 0.0;
                self.electric_z[i] = amplitude * (kx * x + ky * y).sin();
                self.magnetic_x[i] = (ky / omega) * self.electric_z[i];
                self.magnetic_y[i] = -(kx / omega) * self.electric_z[i];
                self.magnetic_z[i] = 0.0;
            }
        }
        Ok(())
    }
    fn bilinear(&self, nodal: &[f64], p: Vec2) -> f64 {
        let w = self.wrap(p);
        cic2(
            w.x / self.dx(),
            w.y / self.dy(),
            self.config.nx,
            self.config.ny,
        )
        .gather(nodal, self.config.nx)
    }
    fn ddx(&self, v: &[f64], ix: usize, iy: usize) -> f64 {
        let nx = self.config.nx;
        (v[fi2((ix + 1) % nx, iy, nx)] - v[fi2((ix + nx - 1) % nx, iy, nx)]) / (2.0 * self.dx())
    }
    fn ddy(&self, v: &[f64], ix: usize, iy: usize) -> f64 {
        let (nx, ny) = (self.config.nx, self.config.ny);
        (v[fi2(ix, (iy + 1) % ny, nx)] - v[fi2(ix, (iy + ny - 1) % ny, nx)]) / (2.0 * self.dy())
    }
    fn solve_longitudinal(&mut self) -> Result<()> {
        if self.config.solve_longitudinal_poisson {
            let c = self.config;
            let f = periodic_electric_field_from_charge_density_2d(
                &self.charge_density,
                c.nx,
                c.ny,
                c.length_x_m,
                c.length_y_m,
                c.relative_permittivity,
                c.neutralize_mean_charge,
            )?;
            self.electric_x = f.electric_x_v_per_m;
            self.electric_y = f.electric_y_v_per_m;
        }
        Ok(())
    }
    pub fn deposit_sources(&mut self) -> Result<()> {
        let c = self.config;
        self.charge_density =
            deposit_cic_charge_density_2d(&self.particles, c.nx, c.ny, c.length_x_m, c.length_y_m)?;
        for v in [
            &mut self.current_x,
            &mut self.current_y,
            &mut self.current_z,
        ] {
            v.iter_mut().for_each(|x| *x = 0.0);
        }
        self.solve_longitudinal()
    }
    pub fn gather_field(&self, p: Vec2) -> ElectromagneticField {
        ElectromagneticField {
            electric_v_per_m: Vec3::new(
                self.bilinear(&self.electric_x, p),
                self.bilinear(&self.electric_y, p),
                self.bilinear(&self.electric_z, p),
            ),
            magnetic_t: Vec3::new(
                self.bilinear(&self.magnetic_x, p),
                self.bilinear(&self.magnetic_y, p),
                self.bilinear(&self.magnetic_z, p),
            ),
        }
    }
    fn advance_magnetic(&mut self, dt: f64) {
        let (mut nbx, mut nby, mut nbz) = (
            self.magnetic_x.clone(),
            self.magnetic_y.clone(),
            self.magnetic_z.clone(),
        );
        let nx = self.config.nx;
        for iy in 0..self.config.ny {
            for ix in 0..nx {
                let i = fi2(ix, iy, nx);
                nbx[i] -= dt * self.ddy(&self.electric_z, ix, iy);
                nby[i] += dt * self.ddx(&self.electric_z, ix, iy);
                nbz[i] -=
                    dt * (self.ddx(&self.electric_y, ix, iy) - self.ddy(&self.electric_x, ix, iy));
            }
        }
        self.magnetic_x = nbx;
        self.magnetic_y = nby;
        self.magnetic_z = nbz;
    }
    fn deposit_current_z(&mut self, old: &[PicParticle2D]) {
        self.current_z.iter_mut().for_each(|v| *v = 0.0);
        let (lx, ly) = (self.config.length_x_m, self.config.length_y_m);
        let area = self.dx() * self.dy();
        for (po, pn) in old.iter().zip(&self.particles) {
            let mid = unwrap_midpoint(po.position_m, pn.position_m, lx, ly);
            let mid = self.wrap(mid);
            let w = cic2(
                mid.x / self.dx(),
                mid.y / self.dy(),
                self.config.nx,
                self.config.ny,
            );
            let vz = 0.5 * (po.velocity_m_per_s.z + pn.velocity_m_per_s.z);
            w.scatter(
                &mut self.current_z,
                self.config.nx,
                pn.charge_c * pn.weight * vz / area,
            );
        }
    }
    pub fn step(&mut self, steps: usize, collisions: Option<&NeutralCollisionModel>) -> Result<()> {
        let c = self.config;
        for _ in 0..steps {
            let old = self.particles.clone();
            self.advance_magnetic(0.5 * c.dt_s);
            for i in 0..self.particles.len() {
                let mut proxy = self.particles[i].proxy();
                let f = self.gather_field(self.particles[i].position_m);
                vay_push(&mut proxy, &f, c.dt_s)?;
                self.particles[i].velocity_m_per_s = proxy.velocity_m_per_s;
                self.particles[i].position_m =
                    self.wrap(Vec2::new(proxy.position_m.x, proxy.position_m.y));
            }
            let cc = deposit_charge_conserving_current_2d(
                &old,
                &self.particles,
                c.nx,
                c.ny,
                c.length_x_m,
                c.length_y_m,
                c.dt_s,
            )?;
            self.charge_density = cc.charge_density_new_c_per_m3;
            self.current_x = cc.current_x_a_per_m2;
            self.current_y = cc.current_y_a_per_m2;
            self.last_continuity_residual = cc.continuity_linf_residual;
            self.deposit_current_z(&old);
            let (mut nex, mut ney, mut nez) = (
                self.electric_x.clone(),
                self.electric_y.clone(),
                self.electric_z.clone(),
            );
            let (eps, mu) = (self.epsilon(), self.mu());
            for iy in 0..c.ny {
                for ix in 0..c.nx {
                    let i = fi2(ix, iy, c.nx);
                    nex[i] += c.dt_s
                        * ((self.ddy(&self.magnetic_z, ix, iy) / mu) - self.current_x[i])
                        / eps;
                    ney[i] += c.dt_s
                        * ((-self.ddx(&self.magnetic_z, ix, iy) / mu) - self.current_y[i])
                        / eps;
                    nez[i] += c.dt_s
                        * (((self.ddx(&self.magnetic_y, ix, iy)
                            - self.ddy(&self.magnetic_x, ix, iy))
                            / mu)
                            - self.current_z[i])
                        / eps;
                }
            }
            self.electric_x = nex;
            self.electric_y = ney;
            self.electric_z = nez;
            self.advance_magnetic(0.5 * c.dt_s);
            self.solve_longitudinal()?;
            if let Some(model) = collisions {
                let stats = collide(
                    &mut self.particles,
                    c.dt_s,
                    model,
                    "invalid collision particle mass",
                )?;
                self.accumulated_collisions.accumulate(&stats);
                for i in 0..self.particles.len() {
                    self.particles[i].position_m = self.wrap(self.particles[i].position_m);
                }
            }
            self.time_s += c.dt_s;
        }
        Ok(())
    }
    pub fn diagnostics(&self) -> Result<ElectromagneticPic2DDiagnostics> {
        let mut d = ElectromagneticPic2DDiagnostics::default();
        let (area, eps, mu) = (self.dx() * self.dy(), self.epsilon(), self.mu());
        for i in 0..self.electric_x.len() {
            let e2 = self.electric_x[i] * self.electric_x[i]
                + self.electric_y[i] * self.electric_y[i]
                + self.electric_z[i] * self.electric_z[i];
            let b2 = self.magnetic_x[i] * self.magnetic_x[i]
                + self.magnetic_y[i] * self.magnetic_y[i]
                + self.magnetic_z[i] * self.magnetic_z[i];
            d.field_energy_j += 0.5 * eps * e2 * area + 0.5 * b2 * area / mu;
        }
        for p in &self.particles {
            d.particle_kinetic_energy_j += p.weight * kinetic_energy_j(&p.proxy())?;
        }
        d.total_energy_j = d.field_energy_j + d.particle_kinetic_energy_j;
        d.charge_continuity_linf_residual = self.last_continuity_residual;
        d.collision_statistics = self.accumulated_collisions;
        d.particle_count = self.particles.len();
        Ok(d)
    }
}

/// Midpoint of a periodic move, unwrapping a jump larger than half the domain.
pub(crate) fn unwrap_midpoint(old: Vec2, new: Vec2, lx: f64, ly: f64) -> Vec2 {
    let (mut o, mut n) = (old, new);
    let dx = n.x - o.x;
    if dx > 0.5 * lx {
        o.x += lx;
    } else if dx < -0.5 * lx {
        n.x += lx;
    }
    let dy = n.y - o.y;
    if dy > 0.5 * ly {
        o.y += ly;
    } else if dy < -0.5 * ly {
        n.y += ly;
    }
    Vec2::new(0.5 * (o.x + n.x), 0.5 * (o.y + n.y))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PicBoundaryBox2D {
    pub minimum_m: Vec2,
    pub maximum_m: Vec2,
}

impl Default for PicBoundaryBox2D {
    fn default() -> Self {
        Self {
            minimum_m: Vec2::default(),
            maximum_m: Vec2::new(1.0, 1.0),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleBoundary2DResult {
    pub alive: bool,
    pub impacted: bool,
    pub impact_normal: Vec2,
    pub secondary_macro_weight: f64,
}

impl Default for ParticleBoundary2DResult {
    fn default() -> Self {
        Self {
            alive: true,
            impacted: false,
            impact_normal: Vec2::default(),
            secondary_macro_weight: 0.0,
        }
    }
}

pub fn apply_particle_box_boundary_2d(
    p: &mut PicParticle2D,
    b: &PicBoundaryBox2D,
    mode: ParticleWallMode,
    secondary: Option<&SecondaryEmissionModel>,
) -> Result<ParticleBoundary2DResult> {
    require(
        b.maximum_m.x > b.minimum_m.x && b.maximum_m.y > b.minimum_m.y,
        "invalid 2-D particle boundary box",
    )?;
    let mut r = ParticleBoundary2DResult::default();
    let outside = |v: f64, lo: f64, hi: f64| v < lo || v > hi;
    if !outside(p.position_m.x, b.minimum_m.x, b.maximum_m.x)
        && !outside(p.position_m.y, b.minimum_m.y, b.maximum_m.y)
    {
        return Ok(r);
    }
    r.impacted = true;
    let incident_ev = kinetic_energy_j(&p.proxy())? / ELEMENTARY_CHARGE;
    let axis = |pos: &mut f64,
                vel: &mut f64,
                lo: f64,
                hi: f64,
                low: Vec2,
                high: Vec2,
                r: &mut ParticleBoundary2DResult| {
        if *pos >= lo && *pos <= hi {
            return;
        }
        r.impact_normal = if *pos < lo { low } else { high };
        if mode == ParticleWallMode::Absorb {
            r.alive = false;
            return;
        }
        reflect_into(pos, vel, lo, hi);
    };
    axis(
        &mut p.position_m.x,
        &mut p.velocity_m_per_s.x,
        b.minimum_m.x,
        b.maximum_m.x,
        Vec2::new(-1.0, 0.0),
        Vec2::new(1.0, 0.0),
        &mut r,
    );
    if r.alive {
        axis(
            &mut p.position_m.y,
            &mut p.velocity_m_per_s.y,
            b.minimum_m.y,
            b.maximum_m.y,
            Vec2::new(0.0, -1.0),
            Vec2::new(0.0, 1.0),
            &mut r,
        );
    }
    if !r.alive {
        if let Some(s) = secondary {
            r.secondary_macro_weight = p.weight * secondary_electron_yield(incident_ev, s)?;
        }
    }
    Ok(r)
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum GridBoundaryMode2D {
    #[default]
    Periodic,
    ElectricWall,
    AbsorbingSponge,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridBoundary2DConfig {
    pub mode: GridBoundaryMode2D,
    pub sponge_cells: usize,
    pub sponge_strength: f64,
    /// Kept for API parity; upstream's sponge profile is quadratic.
    pub sponge_polynomial_order: f64,
}

impl Default for GridBoundary2DConfig {
    fn default() -> Self {
        Self {
            mode: GridBoundaryMode2D::Periodic,
            sponge_cells: 2,
            sponge_strength: 3.0,
            sponge_polynomial_order: 2.0,
        }
    }
}

/// Electric walls zero the tangential component on the boundary; the sponge damps both
/// components by `exp(−s·(1 − d/L)²)` within `L` cells of the edge.
pub fn apply_electrostatic_field_boundary_2d(
    ex: &mut [f64],
    ey: &mut [f64],
    nx: usize,
    ny: usize,
    b: &GridBoundary2DConfig,
) -> Result<()> {
    require(
        ex.len() == nx * ny && ey.len() == nx * ny && nx >= 2 && ny >= 2,
        "invalid 2-D field boundary arrays",
    )?;
    match b.mode {
        GridBoundaryMode2D::Periodic => {}
        GridBoundaryMode2D::ElectricWall => {
            for iy in 0..ny {
                ey[fi2(0, iy, nx)] = 0.0;
                ey[fi2(nx - 1, iy, nx)] = 0.0;
            }
            for ix in 0..nx {
                ex[fi2(ix, 0, nx)] = 0.0;
                ex[fi2(ix, ny - 1, nx)] = 0.0;
            }
        }
        GridBoundaryMode2D::AbsorbingSponge => {
            let layers = b.sponge_cells.min(nx.min(ny) / 2);
            if layers == 0 {
                return Ok(());
            }
            require(
                b.sponge_strength >= 0.0,
                "invalid absorbing sponge strength",
            )?;
            for iy in 0..ny {
                for ix in 0..nx {
                    let edge = ix.min(nx - 1 - ix).min(iy.min(ny - 1 - iy));
                    if edge >= layers {
                        continue;
                    }
                    let normalized = 1.0 - edge as f64 / layers as f64;
                    let damping = (-b.sponge_strength * normalized * normalized).exp();
                    let i = fi2(ix, iy, nx);
                    ex[i] *= damping;
                    ey[i] *= damping;
                }
            }
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StaggeredElectromagneticPic2DConfig {
    pub nx: usize,
    pub ny: usize,
    pub length_x_m: f64,
    pub length_y_m: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub dt_s: f64,
    pub field_boundary: GridBoundary2DConfig,
    pub particle_boundary: ParticleWallMode,
    pub periodic_particles: bool,
    pub remove_absorbed_particles: bool,
    pub enable_secondary_yield_report: bool,
    pub secondary_model: SecondaryEmissionModel,
}

impl Default for StaggeredElectromagneticPic2DConfig {
    fn default() -> Self {
        Self {
            nx: 32,
            ny: 32,
            length_x_m: 1.0,
            length_y_m: 1.0,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            dt_s: 1.0e-12,
            field_boundary: GridBoundary2DConfig::default(),
            particle_boundary: ParticleWallMode::SpecularReflect,
            periodic_particles: true,
            remove_absorbed_particles: true,
            enable_secondary_yield_report: false,
            secondary_model: SecondaryEmissionModel::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StaggeredElectromagneticPic2DDiagnostics {
    pub field_energy_j: f64,
    pub particle_kinetic_energy_j: f64,
    pub total_energy_j: f64,
    pub charge_continuity_linf_residual: f64,
    pub particle_count: usize,
    pub absorbed_particles: usize,
    pub secondary_macro_weight: f64,
    pub collision_statistics: CollisionStatistics,
}

/// Yee-staggered 2-D EM PIC (Ex at (i+½,j), Ey at (i,j+½), Ez at nodes, Bx at (i,j+½),
/// By at (i+½,j), Bz at (i+½,j+½)) with periodic, electric-wall or sponge field boundaries
/// and periodic, reflecting or absorbing particle walls.
#[derive(Clone, Debug)]
pub struct StaggeredElectromagneticPic2D {
    config: StaggeredElectromagneticPic2DConfig,
    particles: Vec<PicParticle2D>,
    charge_density: Vec<f64>,
    current_x: Vec<f64>,
    current_y: Vec<f64>,
    current_z: Vec<f64>,
    electric_x: Vec<f64>,
    electric_y: Vec<f64>,
    electric_z: Vec<f64>,
    magnetic_x: Vec<f64>,
    magnetic_y: Vec<f64>,
    magnetic_z: Vec<f64>,
    time_s: f64,
    last_continuity_residual: f64,
    absorbed_particles: usize,
    secondary_macro_weight: f64,
    accumulated_collisions: CollisionStatistics,
}

pub(crate) fn bounded_index(value: f64, size: usize) -> usize {
    if value <= 0.0 {
        return 0;
    }
    if value >= (size - 1) as f64 {
        return size - 1;
    }
    value.floor() as usize
}

impl StaggeredElectromagneticPic2D {
    pub fn new(config: StaggeredElectromagneticPic2DConfig) -> Result<Self> {
        let c = &config;
        require(
            c.nx >= 2
                && c.ny >= 2
                && c.length_x_m > 0.0
                && c.length_y_m > 0.0
                && c.relative_permittivity > 0.0
                && c.relative_permeability > 0.0
                && c.dt_s > 0.0
                && c.dt_s.is_finite(),
            "invalid staggered 2-D EM PIC configuration",
        )?;
        let n = c.nx * c.ny;
        let z = vec![0.0; n];
        let s = Self {
            config,
            particles: Vec::new(),
            charge_density: z.clone(),
            current_x: z.clone(),
            current_y: z.clone(),
            current_z: z.clone(),
            electric_x: z.clone(),
            electric_y: z.clone(),
            electric_z: z.clone(),
            magnetic_x: z.clone(),
            magnetic_y: z.clone(),
            magnetic_z: z,
            time_s: 0.0,
            last_continuity_residual: 0.0,
            absorbed_particles: 0,
            secondary_macro_weight: 0.0,
            accumulated_collisions: CollisionStatistics::default(),
        };
        let (dx, dy) = (s.dx(), s.dy());
        require(
            s.config.dt_s * s.light() * (1.0 / (dx * dx) + 1.0 / (dy * dy)).sqrt() <= 1.0,
            "staggered 2-D EM PIC CFL must be <= 1",
        )?;
        Ok(s)
    }
    field_getters!(
        charge_density,
        current_x,
        current_y,
        current_z,
        electric_x,
        electric_y,
        electric_z,
        magnetic_x,
        magnetic_y,
        magnetic_z
    );
    pub fn particles(&self) -> &[PicParticle2D] {
        &self.particles
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    fn dx(&self) -> f64 {
        self.config.length_x_m / self.config.nx as f64
    }
    fn dy(&self) -> f64 {
        self.config.length_y_m / self.config.ny as f64
    }
    fn epsilon(&self) -> f64 {
        EPSILON0 * self.config.relative_permittivity
    }
    fn mu(&self) -> f64 {
        MU0 * self.config.relative_permeability
    }
    fn light(&self) -> f64 {
        1.0 / (self.epsilon() * self.mu()).sqrt()
    }
    fn periodic_fields(&self) -> bool {
        self.config.field_boundary.mode == GridBoundaryMode2D::Periodic
    }
    fn place(&self, p: Vec2) -> Vec2 {
        if self.config.periodic_particles {
            Vec2::new(
                wrap_periodic(p.x, self.config.length_x_m),
                wrap_periodic(p.y, self.config.length_y_m),
            )
        } else {
            p
        }
    }
    pub fn set_particles(&mut self, mut particles: Vec<PicParticle2D>) -> Result<()> {
        for p in &mut particles {
            require(
                p.mass_kg > 0.0
                    && p.mass_kg.is_finite()
                    && p.charge_c.is_finite()
                    && p.weight.is_finite(),
                "invalid staggered 2-D EM PIC particle",
            )?;
            p.position_m = self.place(p.position_m);
        }
        self.particles = particles;
        Ok(())
    }
    pub fn set_fields(
        &mut self,
        ex: &[f64],
        ey: &[f64],
        ez: &[f64],
        bx: &[f64],
        by: &[f64],
        bz: &[f64],
    ) -> Result<()> {
        let n = self.config.nx * self.config.ny;
        require(
            [ex, ey, ez, bx, by, bz].iter().all(|f| f.len() == n),
            "staggered 2-D EM PIC field arrays must match grid size",
        )?;
        self.electric_x = ex.to_vec();
        self.electric_y = ey.to_vec();
        self.electric_z = ez.to_vec();
        self.magnetic_x = bx.to_vec();
        self.magnetic_y = by.to_vec();
        self.magnetic_z = bz.to_vec();
        self.apply_field_boundary()
    }
    pub fn initialize_tm_z_mode(
        &mut self,
        amplitude: f64,
        mode_x: usize,
        mode_y: usize,
    ) -> Result<()> {
        require(
            mode_x != 0 || mode_y != 0,
            "staggered 2-D EM PIC mode index must be non-zero",
        )?;
        let kx = 2.0 * PI * mode_x as f64 / self.config.length_x_m;
        let ky = 2.0 * PI * mode_y as f64 / self.config.length_y_m;
        let omega = self.light() * kx.hypot(ky);
        let (dx, dy, nx) = (self.dx(), self.dy(), self.config.nx);
        for iy in 0..self.config.ny {
            for ix in 0..nx {
                let i = fi2(ix, iy, nx);
                let (xf, yf) = (ix as f64, iy as f64);
                self.electric_x[i] = 0.0;
                self.electric_y[i] = 0.0;
                self.electric_z[i] = amplitude * (kx * (xf * dx) + ky * (yf * dy)).sin();
                let phase_bx = kx * xf * dx + ky * (yf + 0.5) * dy;
                let phase_by = kx * (xf + 0.5) * dx + ky * yf * dy;
                self.magnetic_x[i] = (ky / omega) * amplitude * phase_bx.sin();
                self.magnetic_y[i] = -(kx / omega) * amplitude * phase_by.sin();
                self.magnetic_z[i] = 0.0;
            }
        }
        self.apply_field_boundary()
    }
    fn sample(&self, values: &[f64], p: Vec2, ox: f64, oy: f64) -> f64 {
        let c = &self.config;
        let (nx, ny) = (c.nx, c.ny);
        if self.periodic_fields() {
            let px = wrap_periodic(p.x, c.length_x_m);
            let py = wrap_periodic(p.y, c.length_y_m);
            let mut ux = px / self.dx() - ox;
            let mut uy = py / self.dy() - oy;
            ux -= (ux / nx as f64).floor() * nx as f64;
            uy -= (uy / ny as f64).floor() * ny as f64;
            return cic2(ux, uy, nx, ny).gather(values, nx);
        }
        let px = p.x.max(0.0).min(c.length_x_m);
        let py = p.y.max(0.0).min(c.length_y_m);
        let ux = (px / self.dx() - ox).max(0.0).min((nx - 1) as f64);
        let uy = (py / self.dy() - oy).max(0.0).min((ny - 1) as f64);
        let ix0 = bounded_index(ux, nx);
        let iy0 = bounded_index(uy, ny);
        let w = Cic2 {
            ix0,
            iy0,
            ix1: (ix0 + 1).min(nx - 1),
            iy1: (iy0 + 1).min(ny - 1),
            fx: (ux - ix0 as f64).clamp(0.0, 1.0),
            fy: (uy - iy0 as f64).clamp(0.0, 1.0),
        };
        w.gather(values, nx)
    }
    fn dxf(&self, v: &[f64], ix: usize, iy: usize) -> f64 {
        let nx = self.config.nx;
        let r = if self.periodic_fields() {
            (ix + 1) % nx
        } else if ix + 1 < nx {
            ix + 1
        } else {
            ix
        };
        (v[fi2(r, iy, nx)] - v[fi2(ix, iy, nx)]) / self.dx()
    }
    fn dyf(&self, v: &[f64], ix: usize, iy: usize) -> f64 {
        let (nx, ny) = (self.config.nx, self.config.ny);
        let u = if self.periodic_fields() {
            (iy + 1) % ny
        } else if iy + 1 < ny {
            iy + 1
        } else {
            iy
        };
        (v[fi2(ix, u, nx)] - v[fi2(ix, iy, nx)]) / self.dy()
    }
    fn dxb(&self, v: &[f64], ix: usize, iy: usize) -> f64 {
        let nx = self.config.nx;
        let l = if self.periodic_fields() {
            (ix + nx - 1) % nx
        } else {
            ix.saturating_sub(1)
        };
        (v[fi2(ix, iy, nx)] - v[fi2(l, iy, nx)]) / self.dx()
    }
    fn dyb(&self, v: &[f64], ix: usize, iy: usize) -> f64 {
        let (nx, ny) = (self.config.nx, self.config.ny);
        let d = if self.periodic_fields() {
            (iy + ny - 1) % ny
        } else {
            iy.saturating_sub(1)
        };
        (v[fi2(ix, iy, nx)] - v[fi2(ix, d, nx)]) / self.dy()
    }
    pub fn gather_field(&self, p: Vec2) -> ElectromagneticField {
        ElectromagneticField {
            electric_v_per_m: Vec3::new(
                self.sample(&self.electric_x, p, 0.5, 0.0),
                self.sample(&self.electric_y, p, 0.0, 0.5),
                self.sample(&self.electric_z, p, 0.0, 0.0),
            ),
            magnetic_t: Vec3::new(
                self.sample(&self.magnetic_x, p, 0.0, 0.5),
                self.sample(&self.magnetic_y, p, 0.5, 0.0),
                self.sample(&self.magnetic_z, p, 0.5, 0.5),
            ),
        }
    }
    fn advance_magnetic(&mut self, dt: f64) -> Result<()> {
        let (mut nbx, mut nby, mut nbz) = (
            self.magnetic_x.clone(),
            self.magnetic_y.clone(),
            self.magnetic_z.clone(),
        );
        let nx = self.config.nx;
        for iy in 0..self.config.ny {
            for ix in 0..nx {
                let i = fi2(ix, iy, nx);
                nbx[i] -= dt * self.dyf(&self.electric_z, ix, iy);
                nby[i] += dt * self.dxf(&self.electric_z, ix, iy);
                nbz[i] -=
                    dt * (self.dxf(&self.electric_y, ix, iy) - self.dyf(&self.electric_x, ix, iy));
            }
        }
        self.magnetic_x = nbx;
        self.magnetic_y = nby;
        self.magnetic_z = nbz;
        self.apply_field_boundary()
    }
    fn deposit_current_z(&mut self, old: &[PicParticle2D]) {
        self.current_z.iter_mut().for_each(|v| *v = 0.0);
        let c = self.config;
        let area = self.dx() * self.dy();
        let periodic = self.periodic_fields();
        for (po, pn) in old.iter().zip(&self.particles) {
            let mid = if c.periodic_particles {
                unwrap_midpoint(po.position_m, pn.position_m, c.length_x_m, c.length_y_m)
            } else {
                Vec2::new(
                    0.5 * (po.position_m.x + pn.position_m.x),
                    0.5 * (po.position_m.y + pn.position_m.y),
                )
            };
            let mid = self.place(mid);
            let ux = (mid.x / self.dx()).max(0.0).min((c.nx - 1) as f64);
            let uy = (mid.y / self.dy()).max(0.0).min((c.ny - 1) as f64);
            let (flx, fly) = (ux.floor(), uy.floor());
            let ix0 = (flx as usize) % c.nx;
            let iy0 = (fly as usize) % c.ny;
            let w = Cic2 {
                ix0,
                iy0,
                ix1: if periodic {
                    (ix0 + 1) % c.nx
                } else {
                    (ix0 + 1).min(c.nx - 1)
                },
                iy1: if periodic {
                    (iy0 + 1) % c.ny
                } else {
                    (iy0 + 1).min(c.ny - 1)
                },
                fx: ux - flx,
                fy: uy - fly,
            };
            let vz = 0.5 * (po.velocity_m_per_s.z + pn.velocity_m_per_s.z);
            w.scatter(
                &mut self.current_z,
                c.nx,
                pn.charge_c * pn.weight * vz / area,
            );
        }
    }
    fn update_electric(&mut self, dt: f64) -> Result<()> {
        let (mut nex, mut ney, mut nez) = (
            self.electric_x.clone(),
            self.electric_y.clone(),
            self.electric_z.clone(),
        );
        let (eps, mu, nx) = (self.epsilon(), self.mu(), self.config.nx);
        for iy in 0..self.config.ny {
            for ix in 0..nx {
                let i = fi2(ix, iy, nx);
                nex[i] +=
                    dt * ((self.dyb(&self.magnetic_z, ix, iy) / mu) - self.current_x[i]) / eps;
                ney[i] +=
                    dt * ((-self.dxb(&self.magnetic_z, ix, iy) / mu) - self.current_y[i]) / eps;
                nez[i] += dt
                    * (((self.dxb(&self.magnetic_y, ix, iy) - self.dyb(&self.magnetic_x, ix, iy))
                        / mu)
                        - self.current_z[i])
                    / eps;
            }
        }
        self.electric_x = nex;
        self.electric_y = ney;
        self.electric_z = nez;
        self.apply_field_boundary()
    }
    fn apply_field_boundary(&mut self) -> Result<()> {
        let c = self.config;
        let (nx, ny) = (c.nx, c.ny);
        match c.field_boundary.mode {
            GridBoundaryMode2D::Periodic => {}
            GridBoundaryMode2D::ElectricWall => {
                for iy in 0..ny {
                    for ix in [0, nx - 1] {
                        let i = fi2(ix, iy, nx);
                        self.electric_y[i] = 0.0;
                        self.electric_z[i] = 0.0;
                        self.magnetic_x[i] = 0.0;
                    }
                }
                for ix in 0..nx {
                    for iy in [0, ny - 1] {
                        let i = fi2(ix, iy, nx);
                        self.electric_x[i] = 0.0;
                        self.electric_z[i] = 0.0;
                        self.magnetic_y[i] = 0.0;
                    }
                }
            }
            GridBoundaryMode2D::AbsorbingSponge => {
                let layers = c.field_boundary.sponge_cells.min(nx.min(ny) / 2);
                if layers == 0 {
                    return Ok(());
                }
                require(
                    c.field_boundary.sponge_strength >= 0.0,
                    "invalid staggered sponge strength",
                )?;
                for iy in 0..ny {
                    for ix in 0..nx {
                        let edge = ix.min(nx - 1 - ix).min(iy.min(ny - 1 - iy));
                        if edge >= layers {
                            continue;
                        }
                        let normalized = 1.0 - edge as f64 / layers as f64;
                        let damping =
                            (-c.field_boundary.sponge_strength * normalized * normalized).exp();
                        let i = fi2(ix, iy, nx);
                        for f in [
                            &mut self.electric_x,
                            &mut self.electric_y,
                            &mut self.electric_z,
                            &mut self.magnetic_x,
                            &mut self.magnetic_y,
                            &mut self.magnetic_z,
                        ] {
                            f[i] *= damping;
                        }
                    }
                }
            }
        }
        Ok(())
    }
    fn apply_particle_boundary(&mut self) -> Result<()> {
        let c = self.config;
        if c.periodic_particles {
            for i in 0..self.particles.len() {
                self.particles[i].position_m = self.place(self.particles[i].position_m);
            }
            return Ok(());
        }
        let b = PicBoundaryBox2D {
            minimum_m: Vec2::new(0.0, 0.0),
            maximum_m: Vec2::new(c.length_x_m, c.length_y_m),
        };
        let secondary = if c.enable_secondary_yield_report {
            Some(&c.secondary_model)
        } else {
            None
        };
        let mut kept = Vec::with_capacity(self.particles.len());
        for mut p in self.particles.iter().copied() {
            let r = apply_particle_box_boundary_2d(&mut p, &b, c.particle_boundary, secondary)?;
            if r.impacted {
                self.secondary_macro_weight += r.secondary_macro_weight;
            }
            if r.alive || !c.remove_absorbed_particles {
                kept.push(p);
            } else {
                self.absorbed_particles += 1;
            }
        }
        self.particles = kept;
        Ok(())
    }
    pub fn step(&mut self, steps: usize, collisions: Option<&NeutralCollisionModel>) -> Result<()> {
        let c = self.config;
        for _ in 0..steps {
            let old = self.particles.clone();
            self.advance_magnetic(0.5 * c.dt_s)?;
            for i in 0..self.particles.len() {
                let mut proxy = self.particles[i].proxy();
                let f = self.gather_field(self.particles[i].position_m);
                vay_push(&mut proxy, &f, c.dt_s)?;
                self.particles[i].velocity_m_per_s = proxy.velocity_m_per_s;
                self.particles[i].position_m = Vec2::new(proxy.position_m.x, proxy.position_m.y);
            }
            self.apply_particle_boundary()?;
            if self.particles.len() == old.len() {
                let cc = deposit_charge_conserving_current_2d(
                    &old,
                    &self.particles,
                    c.nx,
                    c.ny,
                    c.length_x_m,
                    c.length_y_m,
                    c.dt_s,
                )?;
                self.charge_density = cc.charge_density_new_c_per_m3;
                self.current_x = cc.current_x_a_per_m2;
                self.current_y = cc.current_y_a_per_m2;
                self.last_continuity_residual = cc.continuity_linf_residual;
                self.deposit_current_z(&old);
            } else {
                self.charge_density = deposit_cic_charge_density_2d(
                    &self.particles,
                    c.nx,
                    c.ny,
                    c.length_x_m,
                    c.length_y_m,
                )?;
                for v in [
                    &mut self.current_x,
                    &mut self.current_y,
                    &mut self.current_z,
                ] {
                    v.iter_mut().for_each(|x| *x = 0.0);
                }
                self.last_continuity_residual = 0.0;
            }
            self.update_electric(c.dt_s)?;
            self.advance_magnetic(0.5 * c.dt_s)?;
            if let Some(model) = collisions {
                let stats = collide(
                    &mut self.particles,
                    c.dt_s,
                    model,
                    "invalid collision particle mass",
                )?;
                self.accumulated_collisions.accumulate(&stats);
                for i in 0..self.particles.len() {
                    self.particles[i].position_m = self.place(self.particles[i].position_m);
                }
                self.apply_particle_boundary()?;
            }
            self.time_s += c.dt_s;
        }
        Ok(())
    }
    pub fn diagnostics(&self) -> Result<StaggeredElectromagneticPic2DDiagnostics> {
        let mut d = StaggeredElectromagneticPic2DDiagnostics::default();
        let (area, eps, mu) = (self.dx() * self.dy(), self.epsilon(), self.mu());
        for i in 0..self.electric_x.len() {
            let e2 = self.electric_x[i] * self.electric_x[i]
                + self.electric_y[i] * self.electric_y[i]
                + self.electric_z[i] * self.electric_z[i];
            let b2 = self.magnetic_x[i] * self.magnetic_x[i]
                + self.magnetic_y[i] * self.magnetic_y[i]
                + self.magnetic_z[i] * self.magnetic_z[i];
            d.field_energy_j += 0.5 * eps * e2 * area + 0.5 * b2 * area / mu;
        }
        for p in &self.particles {
            d.particle_kinetic_energy_j += p.weight * kinetic_energy_j(&p.proxy())?;
        }
        d.total_energy_j = d.field_energy_j + d.particle_kinetic_energy_j;
        d.charge_continuity_linf_residual = self.last_continuity_residual;
        d.particle_count = self.particles.len();
        d.absorbed_particles = self.absorbed_particles;
        d.secondary_macro_weight = self.secondary_macro_weight;
        d.collision_statistics = self.accumulated_collisions;
        Ok(d)
    }
}
