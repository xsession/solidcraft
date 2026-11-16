//! 3-D periodic particle-in-cell: CIC deposition, spectral Poisson, spectral and local
//! finite-volume charge-conserving currents, electrostatic, collocated and Yee-staggered
//! electromagnetic PIC with wall/sponge boundaries and cell sorting. Ported from upstream
//! `cfd/particle/electromagnetic` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::pic1d::{phase, signed_mode};
use super::pic3d_domain::{fi3, sort_particles_by_cell_3d, PicParticle3D};
use super::pusher::*;
use crate::{require, Result};
use num_complex::Complex64 as C;
use std::f64::consts::PI;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ElectrostaticField3D {
    pub electric_x_v_per_m: Vec<f64>,
    pub electric_y_v_per_m: Vec<f64>,
    pub electric_z_v_per_m: Vec<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChargeConservingCurrent3D {
    pub charge_density_old_c_per_m3: Vec<f64>,
    pub charge_density_new_c_per_m3: Vec<f64>,
    pub current_x_a_per_m2: Vec<f64>,
    pub current_y_a_per_m2: Vec<f64>,
    pub current_z_a_per_m2: Vec<f64>,
    pub continuity_linf_residual: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Grid3 {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub lx: f64,
    pub ly: f64,
    pub lz: f64,
}

impl Grid3 {
    pub(crate) fn n(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    fn valid(&self) -> bool {
        self.nx >= 2
            && self.ny >= 2
            && self.nz >= 2
            && self.lx > 0.0
            && self.ly > 0.0
            && self.lz > 0.0
    }
    fn angle(&self, m: (usize, usize, usize), i: (usize, usize, usize)) -> f64 {
        2.0 * PI
            * ((m.0 * i.0) as f64 / self.nx as f64
                + (m.1 * i.1) as f64 / self.ny as f64
                + (m.2 * i.2) as f64 / self.nz as f64)
    }
    fn wavevector(&self, mx: usize, my: usize, mz: usize) -> (f64, f64, f64) {
        (
            2.0 * PI * signed_mode(mx, self.nx) / self.lx,
            2.0 * PI * signed_mode(my, self.ny) / self.ly,
            2.0 * PI * signed_mode(mz, self.nz) / self.lz,
        )
    }
    fn fourier(&self, values: &[f64], m: (usize, usize, usize)) -> C {
        let mut c = C::default();
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    c += values[fi3(ix, iy, iz, self.nx, self.ny)]
                        * phase(-self.angle(m, (ix, iy, iz)));
                }
            }
        }
        c / (self.nx * self.ny * self.nz) as f64
    }
    /// Inverse DFT (real part) of several spectra at once.
    fn inverse<const K: usize>(&self, hats: [&[C]; K]) -> [Vec<f64>; K] {
        let n = self.n();
        let mut out: [Vec<f64>; K] = std::array::from_fn(|_| vec![0.0; n]);
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    let mut acc = [C::default(); K];
                    for mz in 0..self.nz {
                        for my in 0..self.ny {
                            for mx in 0..self.nx {
                                let ph = phase(self.angle((mx, my, mz), (ix, iy, iz)));
                                let k = fi3(mx, my, mz, self.nx, self.ny);
                                for c in 0..K {
                                    acc[c] += hats[c][k] * ph;
                                }
                            }
                        }
                    }
                    let i = fi3(ix, iy, iz, self.nx, self.ny);
                    for c in 0..K {
                        out[c][i] = acc[c].re;
                    }
                }
            }
        }
        out
    }
}

pub(crate) struct Cic3 {
    i0: [usize; 3],
    i1: [usize; 3],
    f: [f64; 3],
}

impl Cic3 {
    pub(crate) fn periodic(u: [f64; 3], n: [usize; 3]) -> Self {
        let mut s = Self {
            i0: [0; 3],
            i1: [0; 3],
            f: [0.0; 3],
        };
        for a in 0..3 {
            let fl = u[a].floor();
            s.i0[a] = (fl as usize) % n[a];
            s.i1[a] = (s.i0[a] + 1) % n[a];
            s.f[a] = u[a] - fl;
        }
        s
    }
    pub(crate) fn new(i0: [usize; 3], i1: [usize; 3], f: [f64; 3]) -> Self {
        Self { i0, i1, f }
    }
    fn corner(&self, k: usize) -> (usize, usize, usize, f64, f64, f64) {
        let (dx, dy, dz) = (k & 1, (k >> 1) & 1, (k >> 2) & 1);
        let pick = |a: usize, d: usize| {
            if d == 1 {
                (self.i1[a], self.f[a])
            } else {
                (self.i0[a], 1.0 - self.f[a])
            }
        };
        let ((ix, wx), (iy, wy), (iz, wz)) = (pick(0, dx), pick(1, dy), pick(2, dz));
        (ix, iy, iz, wx, wy, wz)
    }
    /// `grid[i] += value·wx·wy·wz` in upstream's (z, y, x) corner order.
    pub(crate) fn scatter(&self, grid: &mut [f64], nx: usize, ny: usize, value: f64) {
        for k in 0..8 {
            let (ix, iy, iz, wx, wy, wz) = self.corner(k);
            grid[fi3(ix, iy, iz, nx, ny)] += value * wx * wy * wz;
        }
    }
    pub(crate) fn gather(&self, grid: &[f64], nx: usize, ny: usize) -> f64 {
        let mut v = 0.0;
        for k in 0..8 {
            let (ix, iy, iz, wx, wy, wz) = self.corner(k);
            v += grid[fi3(ix, iy, iz, nx, ny)] * wx * wy * wz;
        }
        v
    }
}

pub fn deposit_cic_charge_density_3d(
    particles: &[PicParticle3D],
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<Vec<f64>> {
    let g = Grid3 {
        nx,
        ny,
        nz,
        lx,
        ly,
        lz,
    };
    require(g.valid(), "invalid 3-D CIC grid")?;
    let mut rho = vec![0.0; g.n()];
    let (dx, dy, dz) = (lx / nx as f64, ly / ny as f64, lz / nz as f64);
    for p in particles {
        require(p.valid(), "invalid 3-D PIC particle")?;
        let w = Cic3::periodic(
            [
                wrap_periodic(p.position_m.x, lx) / dx,
                wrap_periodic(p.position_m.y, ly) / dy,
                wrap_periodic(p.position_m.z, lz) / dz,
            ],
            [nx, ny, nz],
        );
        w.scatter(&mut rho, nx, ny, p.charge_c * p.weight / (dx * dy * dz));
    }
    Ok(rho)
}

/// Spectral periodic Poisson in 3-D (direct DFT; intended for small grids).
#[allow(clippy::too_many_arguments)]
pub fn periodic_electric_field_from_charge_density_3d(
    rho: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
    relative_permittivity: f64,
    remove_mean: bool,
) -> Result<ElectrostaticField3D> {
    let g = Grid3 {
        nx,
        ny,
        nz,
        lx,
        ly,
        lz,
    };
    require(
        rho.len() == g.n() && g.valid() && relative_permittivity > 0.0,
        "invalid 3-D periodic Poisson grid",
    )?;
    let (mut mean, mut max_abs) = (0.0, 0.0f64);
    for &v in rho {
        require(v.is_finite(), "non-finite 3-D charge density")?;
        mean += v;
        max_abs = max_abs.max(v.abs());
    }
    mean /= rho.len() as f64;
    require(
        remove_mean || mean.abs() <= 1.0e-12 * max_abs.max(1.0),
        "3-D periodic Poisson requires mean-zero charge density",
    )?;
    let centered: Vec<f64> = rho
        .iter()
        .map(|v| if remove_mean { v - mean } else { *v })
        .collect();
    let eps = EPSILON0 * relative_permittivity;
    let n = g.n();
    let (mut ex, mut ey, mut ez) = (
        vec![C::default(); n],
        vec![C::default(); n],
        vec![C::default(); n],
    );
    for mz in 0..nz {
        for my in 0..ny {
            for mx in 0..nx {
                if mx == 0 && my == 0 && mz == 0 {
                    continue;
                }
                let (kx, ky, kz) = g.wavevector(mx, my, mz);
                let k2 = kx * kx + ky * ky + kz * kz;
                if k2 == 0.0 {
                    continue;
                }
                let factor = C::new(0.0, -1.0) * g.fourier(&centered, (mx, my, mz)) / (eps * k2);
                let k = fi3(mx, my, mz, nx, ny);
                ex[k] = factor * kx;
                ey[k] = factor * ky;
                ez[k] = factor * kz;
            }
        }
    }
    let [x, y, z] = g.inverse([&ex, &ey, &ez]);
    Ok(ElectrostaticField3D {
        electric_x_v_per_m: x,
        electric_y_v_per_m: y,
        electric_z_v_per_m: z,
    })
}

fn cc_prelude(
    old: &[PicParticle3D],
    new: &[PicParticle3D],
    g: &Grid3,
    dt: f64,
    what: &str,
) -> Result<(ChargeConservingCurrent3D, Vec<f64>)> {
    require(
        old.len() == new.len(),
        "old/new 3-D particle arrays must have equal size",
    )?;
    require(g.valid() && dt > 0.0 && dt.is_finite(), what)?;
    let r = ChargeConservingCurrent3D {
        charge_density_old_c_per_m3: deposit_cic_charge_density_3d(
            old, g.nx, g.ny, g.nz, g.lx, g.ly, g.lz,
        )?,
        charge_density_new_c_per_m3: deposit_cic_charge_density_3d(
            new, g.nx, g.ny, g.nz, g.lx, g.ly, g.lz,
        )?,
        current_x_a_per_m2: vec![0.0; g.n()],
        current_y_a_per_m2: vec![0.0; g.n()],
        current_z_a_per_m2: vec![0.0; g.n()],
        continuity_linf_residual: 0.0,
    };
    let delta = r
        .charge_density_new_c_per_m3
        .iter()
        .zip(&r.charge_density_old_c_per_m3)
        .map(|(a, b)| a - b)
        .collect();
    Ok((r, delta))
}

/// Curl-free current whose spectral divergence matches the CIC charge change exactly.
#[allow(clippy::too_many_arguments)]
pub fn deposit_charge_conserving_current_3d(
    old: &[PicParticle3D],
    new: &[PicParticle3D],
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
    dt: f64,
) -> Result<ChargeConservingCurrent3D> {
    let g = Grid3 {
        nx,
        ny,
        nz,
        lx,
        ly,
        lz,
    };
    let (mut r, delta) = cc_prelude(
        old,
        new,
        &g,
        dt,
        "invalid 3-D charge-conserving deposition parameters",
    )?;
    let n = g.n();
    let mut hats: [Vec<C>; 4] = std::array::from_fn(|_| vec![C::default(); n]);
    for mz in 0..nz {
        for my in 0..ny {
            for mx in 0..nx {
                if mx == 0 && my == 0 && mz == 0 {
                    continue;
                }
                let (kx, ky, kz) = g.wavevector(mx, my, mz);
                let k2 = kx * kx + ky * ky + kz * kz;
                if k2 == 0.0 {
                    continue;
                }
                let factor = C::new(0.0, 1.0) * g.fourier(&delta, (mx, my, mz)) / (dt * k2);
                let k = fi3(mx, my, mz, nx, ny);
                hats[0][k] = factor * kx;
                hats[1][k] = factor * ky;
                hats[2][k] = factor * kz;
                hats[3][k] =
                    C::new(0.0, 1.0) * (kx * hats[0][k] + ky * hats[1][k] + kz * hats[2][k]);
            }
        }
    }
    let [jx, jy, jz, div] = g.inverse([&hats[0], &hats[1], &hats[2], &hats[3]]);
    r.current_x_a_per_m2 = jx;
    r.current_y_a_per_m2 = jy;
    r.current_z_a_per_m2 = jz;
    for i in 0..n {
        r.continuity_linf_residual = r
            .continuity_linf_residual
            .max((delta[i] / dt + div[i]).abs());
    }
    Ok(r)
}

/// Local finite-volume current: cumulative x sums per row, then row/plane means along y and z,
/// satisfying the discrete backward-difference continuity equation to round-off.
#[allow(clippy::too_many_arguments)]
pub fn deposit_charge_conserving_current_3d_local(
    old: &[PicParticle3D],
    new: &[PicParticle3D],
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
    dt: f64,
) -> Result<ChargeConservingCurrent3D> {
    let g = Grid3 {
        nx,
        ny,
        nz,
        lx,
        ly,
        lz,
    };
    let (mut r, delta) = cc_prelude(
        old,
        new,
        &g,
        dt,
        "invalid local 3-D current deposition parameters",
    )?;
    let (dx, dy, dz) = (lx / nx as f64, ly / ny as f64, lz / nz as f64);
    let mut mean_x = vec![0.0; ny * nz];
    let mut mean_xy = vec![0.0; nz];
    for iz in 0..nz {
        for iy in 0..ny {
            let mut sum = 0.0;
            for ix in 0..nx {
                sum += delta[fi3(ix, iy, iz, nx, ny)];
            }
            mean_x[iz * ny + iy] = sum / nx as f64;
        }
    }
    for iz in 0..nz {
        let mut sum = 0.0;
        for iy in 0..ny {
            sum += mean_x[iz * ny + iy];
        }
        mean_xy[iz] = sum / ny as f64;
    }
    let mut global_mean = 0.0;
    for v in &mean_xy {
        global_mean += v;
    }
    global_mean /= nz as f64;
    for iz in 0..nz {
        for iy in 0..ny {
            let mut cumulative = 0.0;
            let row = mean_x[iz * ny + iy];
            for ix in 0..nx {
                let i = fi3(ix, iy, iz, nx, ny);
                cumulative += -dx * (delta[i] - row) / dt;
                r.current_x_a_per_m2[i] = cumulative;
            }
        }
    }
    for iz in 0..nz {
        let mut cumulative = 0.0;
        for iy in 0..ny {
            cumulative += -dy * (mean_x[iz * ny + iy] - mean_xy[iz]) / dt;
            for ix in 0..nx {
                r.current_y_a_per_m2[fi3(ix, iy, iz, nx, ny)] = cumulative;
            }
        }
    }
    let mut cumulative = 0.0;
    for iz in 0..nz {
        cumulative += -dz * (mean_xy[iz] - global_mean) / dt;
        for iy in 0..ny {
            for ix in 0..nx {
                r.current_z_a_per_m2[fi3(ix, iy, iz, nx, ny)] = cumulative;
            }
        }
    }
    for iz in 0..nz {
        let zb = (iz + nz - 1) % nz;
        for iy in 0..ny {
            let yb = (iy + ny - 1) % ny;
            for ix in 0..nx {
                let xb = (ix + nx - 1) % nx;
                let i = fi3(ix, iy, iz, nx, ny);
                let div = (r.current_x_a_per_m2[i] - r.current_x_a_per_m2[fi3(xb, iy, iz, nx, ny)])
                    / dx
                    + (r.current_y_a_per_m2[i] - r.current_y_a_per_m2[fi3(ix, yb, iz, nx, ny)])
                        / dy
                    + (r.current_z_a_per_m2[i] - r.current_z_a_per_m2[fi3(ix, iy, zb, nx, ny)])
                        / dz;
                r.continuity_linf_residual =
                    r.continuity_linf_residual.max((delta[i] / dt + div).abs());
            }
        }
    }
    Ok(r)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrostaticPic3DConfig {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub length_x_m: f64,
    pub length_y_m: f64,
    pub length_z_m: f64,
    pub relative_permittivity: f64,
    pub dt_s: f64,
    pub neutralize_mean_charge: bool,
}

impl Default for ElectrostaticPic3DConfig {
    fn default() -> Self {
        Self {
            nx: 16,
            ny: 16,
            nz: 16,
            length_x_m: 1.0,
            length_y_m: 1.0,
            length_z_m: 1.0,
            relative_permittivity: 1.0,
            dt_s: 1.0e-3,
            neutralize_mean_charge: true,
        }
    }
}

fn wrap3(p: Vec3, g: &Grid3) -> Vec3 {
    Vec3::new(
        wrap_periodic(p.x, g.lx),
        wrap_periodic(p.y, g.ly),
        wrap_periodic(p.z, g.lz),
    )
}

fn trilinear_periodic(nodal: &[f64], p: Vec3, g: &Grid3) -> f64 {
    let w = wrap3(p, g);
    Cic3::periodic(
        [
            w.x / (g.lx / g.nx as f64),
            w.y / (g.ly / g.ny as f64),
            w.z / (g.lz / g.nz as f64),
        ],
        [g.nx, g.ny, g.nz],
    )
    .gather(nodal, g.nx, g.ny)
}

#[derive(Clone, Debug)]
pub struct ElectrostaticPic3D {
    config: ElectrostaticPic3DConfig,
    particles: Vec<PicParticle3D>,
    charge_density: Vec<f64>,
    electric_x: Vec<f64>,
    electric_y: Vec<f64>,
    electric_z: Vec<f64>,
    time_s: f64,
}

impl ElectrostaticPic3D {
    pub fn new(config: ElectrostaticPic3DConfig) -> Result<Self> {
        let s = Self {
            charge_density: Vec::new(),
            electric_x: Vec::new(),
            electric_y: Vec::new(),
            electric_z: Vec::new(),
            particles: Vec::new(),
            time_s: 0.0,
            config,
        };
        let g = s.grid();
        require(
            g.valid()
                && config.relative_permittivity > 0.0
                && config.dt_s > 0.0
                && config.dt_s.is_finite(),
            "invalid electrostatic PIC 3-D configuration",
        )?;
        let n = g.n();
        Ok(Self {
            charge_density: vec![0.0; n],
            electric_x: vec![0.0; n],
            electric_y: vec![0.0; n],
            electric_z: vec![0.0; n],
            ..s
        })
    }
    fn grid(&self) -> Grid3 {
        let c = &self.config;
        Grid3 {
            nx: c.nx,
            ny: c.ny,
            nz: c.nz,
            lx: c.length_x_m,
            ly: c.length_y_m,
            lz: c.length_z_m,
        }
    }
    super::pic2d::field_getters!(charge_density, electric_x, electric_y, electric_z);
    pub fn particles(&self) -> &[PicParticle3D] {
        &self.particles
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn set_particles(&mut self, mut particles: Vec<PicParticle3D>) -> Result<()> {
        let g = self.grid();
        for p in &mut particles {
            require(p.valid(), "invalid electrostatic PIC 3-D particle")?;
            p.position_m = wrap3(p.position_m, &g);
        }
        self.particles = particles;
        Ok(())
    }
    pub fn deposit_and_solve(&mut self) -> Result<()> {
        let (c, g) = (self.config, self.grid());
        self.charge_density =
            deposit_cic_charge_density_3d(&self.particles, g.nx, g.ny, g.nz, g.lx, g.ly, g.lz)?;
        let f = periodic_electric_field_from_charge_density_3d(
            &self.charge_density,
            g.nx,
            g.ny,
            g.nz,
            g.lx,
            g.ly,
            g.lz,
            c.relative_permittivity,
            c.neutralize_mean_charge,
        )?;
        self.electric_x = f.electric_x_v_per_m;
        self.electric_y = f.electric_y_v_per_m;
        self.electric_z = f.electric_z_v_per_m;
        Ok(())
    }
    pub fn gather_electric_field(&self, p: Vec3) -> Vec3 {
        let g = self.grid();
        Vec3::new(
            trilinear_periodic(&self.electric_x, p, &g),
            trilinear_periodic(&self.electric_y, p, &g),
            trilinear_periodic(&self.electric_z, p, &g),
        )
    }
    pub fn step(&mut self, steps: usize) -> Result<()> {
        let (dt, g) = (self.config.dt_s, self.grid());
        for _ in 0..steps {
            self.deposit_and_solve()?;
            for i in 0..self.particles.len() {
                let mut p = self.particles[i];
                let f = self.gather_electric_field(p.position_m);
                let qm = p.charge_c / p.mass_kg;
                p.velocity_m_per_s.x += qm * f.x * dt;
                p.velocity_m_per_s.y += qm * f.y * dt;
                p.velocity_m_per_s.z += qm * f.z * dt;
                p.position_m.x += p.velocity_m_per_s.x * dt;
                p.position_m.y += p.velocity_m_per_s.y * dt;
                p.position_m.z += p.velocity_m_per_s.z * dt;
                p.position_m = wrap3(p.position_m, &g);
                self.particles[i] = p;
            }
            self.time_s += dt;
        }
        self.deposit_and_solve()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectromagneticPic3DConfig {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub length_x_m: f64,
    pub length_y_m: f64,
    pub length_z_m: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub dt_s: f64,
    pub solve_longitudinal_poisson: bool,
    pub neutralize_mean_charge: bool,
}

impl Default for ElectromagneticPic3DConfig {
    fn default() -> Self {
        Self {
            nx: 12,
            ny: 12,
            nz: 12,
            length_x_m: 1.0,
            length_y_m: 1.0,
            length_z_m: 1.0,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            dt_s: 1.0e-12,
            solve_longitudinal_poisson: false,
            neutralize_mean_charge: true,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ElectromagneticPic3DDiagnostics {
    pub field_energy_j: f64,
    pub particle_kinetic_energy_j: f64,
    pub total_energy_j: f64,
    pub charge_continuity_linf_residual: f64,
    pub collision_statistics: CollisionStatistics,
    pub particle_count: usize,
}

/// Shared six-component field storage.
#[derive(Clone, Debug, Default)]
struct Fields3 {
    ex: Vec<f64>,
    ey: Vec<f64>,
    ez: Vec<f64>,
    bx: Vec<f64>,
    by: Vec<f64>,
    bz: Vec<f64>,
}

impl Fields3 {
    fn zeros(n: usize) -> Self {
        let z = vec![0.0; n];
        Self {
            ex: z.clone(),
            ey: z.clone(),
            ez: z.clone(),
            bx: z.clone(),
            by: z.clone(),
            bz: z,
        }
    }
    fn set(&mut self, f: [&[f64]; 6], n: usize, what: &str) -> Result<()> {
        require(f.iter().all(|v| v.len() == n), what)?;
        self.ex = f[0].to_vec();
        self.ey = f[1].to_vec();
        self.ez = f[2].to_vec();
        self.bx = f[3].to_vec();
        self.by = f[4].to_vec();
        self.bz = f[5].to_vec();
        Ok(())
    }
    fn energy(&self, eps: f64, mu: f64, volume: f64) -> f64 {
        let mut e = 0.0;
        for i in 0..self.ex.len() {
            let e2 = self.ex[i] * self.ex[i] + self.ey[i] * self.ey[i] + self.ez[i] * self.ez[i];
            let b2 = self.bx[i] * self.bx[i] + self.by[i] * self.by[i] + self.bz[i] * self.bz[i];
            e += 0.5 * eps * e2 * volume + 0.5 * b2 * volume / mu;
        }
        e
    }
}

/// Collocated 3-D EM PIC with centred differences (CFL guard 0.55), Vay pushing and spectral
/// charge-conserving currents.
#[derive(Clone, Debug)]
pub struct ElectromagneticPic3D {
    config: ElectromagneticPic3DConfig,
    particles: Vec<PicParticle3D>,
    charge_density: Vec<f64>,
    current_x: Vec<f64>,
    current_y: Vec<f64>,
    current_z: Vec<f64>,
    f: Fields3,
    time_s: f64,
    last_continuity_residual: f64,
    accumulated_collisions: CollisionStatistics,
}

macro_rules! fields3_getters {
    () => {
        pub fn electric_x(&self) -> &[f64] {
            &self.f.ex
        }
        pub fn electric_y(&self) -> &[f64] {
            &self.f.ey
        }
        pub fn electric_z(&self) -> &[f64] {
            &self.f.ez
        }
        pub fn magnetic_x(&self) -> &[f64] {
            &self.f.bx
        }
        pub fn magnetic_y(&self) -> &[f64] {
            &self.f.by
        }
        pub fn magnetic_z(&self) -> &[f64] {
            &self.f.bz
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
        pub fn particles(&self) -> &[PicParticle3D] {
            &self.particles
        }
        pub fn time_s(&self) -> f64 {
            self.time_s
        }
        fn grid(&self) -> Grid3 {
            let c = &self.config;
            Grid3 {
                nx: c.nx,
                ny: c.ny,
                nz: c.nz,
                lx: c.length_x_m,
                ly: c.length_y_m,
                lz: c.length_z_m,
            }
        }
        fn dx(&self) -> f64 {
            self.config.length_x_m / self.config.nx as f64
        }
        fn dy(&self) -> f64 {
            self.config.length_y_m / self.config.ny as f64
        }
        fn dz(&self) -> f64 {
            self.config.length_z_m / self.config.nz as f64
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
        fn inv_cfl(&self) -> f64 {
            let (dx, dy, dz) = (self.dx(), self.dy(), self.dz());
            (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)).sqrt()
        }
    };
}

impl ElectromagneticPic3D {
    pub fn new(config: ElectromagneticPic3DConfig) -> Result<Self> {
        let s = Self {
            config,
            particles: Vec::new(),
            charge_density: Vec::new(),
            current_x: Vec::new(),
            current_y: Vec::new(),
            current_z: Vec::new(),
            f: Fields3::default(),
            time_s: 0.0,
            last_continuity_residual: 0.0,
            accumulated_collisions: CollisionStatistics::default(),
        };
        let g = s.grid();
        require(
            g.valid()
                && config.relative_permittivity > 0.0
                && config.relative_permeability > 0.0
                && config.dt_s > 0.0
                && config.dt_s.is_finite(),
            "invalid electromagnetic PIC 3-D configuration",
        )?;
        require(
            config.dt_s * s.light() * s.inv_cfl() < 0.55,
            "3-D EM-PIC time step violates compact centered-CFL guard",
        )?;
        let n = g.n();
        Ok(Self {
            charge_density: vec![0.0; n],
            current_x: vec![0.0; n],
            current_y: vec![0.0; n],
            current_z: vec![0.0; n],
            f: Fields3::zeros(n),
            ..s
        })
    }
    fields3_getters!();
    pub fn set_particles(&mut self, mut particles: Vec<PicParticle3D>) -> Result<()> {
        let g = self.grid();
        for p in &mut particles {
            require(p.valid(), "invalid electromagnetic PIC 3-D particle")?;
            p.position_m = wrap3(p.position_m, &g);
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
        let n = self.grid().n();
        self.f.set(
            [ex, ey, ez, bx, by, bz],
            n,
            "3-D EM-PIC field arrays must match grid size",
        )
    }
    pub fn initialize_z_polarized_mode(
        &mut self,
        amplitude: f64,
        mx: usize,
        my: usize,
        mz: usize,
    ) -> Result<()> {
        require(
            mx != 0 || my != 0 || mz != 0,
            "3-D EM-PIC mode cannot be DC",
        )?;
        let c = self.config;
        self.f = Fields3::zeros(self.grid().n());
        let kx = 2.0 * PI * mx as f64 / c.length_x_m;
        let ky = 2.0 * PI * my as f64 / c.length_y_m;
        let kz = 2.0 * PI * mz as f64 / c.length_z_m;
        let k = (kx * kx + ky * ky + kz * kz).sqrt();
        let light = self.light();
        for iz in 0..c.nz {
            let z = c.length_z_m * iz as f64 / c.nz as f64;
            for iy in 0..c.ny {
                let y = c.length_y_m * iy as f64 / c.ny as f64;
                for ix in 0..c.nx {
                    let x = c.length_x_m * ix as f64 / c.nx as f64;
                    let s = (kx * x + ky * y + kz * z).sin();
                    let i = fi3(ix, iy, iz, c.nx, c.ny);
                    self.f.ez[i] = amplitude * s;
                    self.f.bx[i] = (ky / k) * (amplitude / light) * s;
                    self.f.by[i] = -(kx / k) * (amplitude / light) * s;
                }
            }
        }
        Ok(())
    }
    fn solve_longitudinal(&mut self) -> Result<()> {
        if self.config.solve_longitudinal_poisson {
            let (c, g) = (self.config, self.grid());
            let f = periodic_electric_field_from_charge_density_3d(
                &self.charge_density,
                g.nx,
                g.ny,
                g.nz,
                g.lx,
                g.ly,
                g.lz,
                c.relative_permittivity,
                c.neutralize_mean_charge,
            )?;
            self.f.ex = f.electric_x_v_per_m;
            self.f.ey = f.electric_y_v_per_m;
            self.f.ez = f.electric_z_v_per_m;
        }
        Ok(())
    }
    pub fn deposit_sources(&mut self) -> Result<()> {
        let g = self.grid();
        self.charge_density =
            deposit_cic_charge_density_3d(&self.particles, g.nx, g.ny, g.nz, g.lx, g.ly, g.lz)?;
        self.solve_longitudinal()
    }
    pub fn gather_field(&self, p: Vec3) -> ElectromagneticField {
        let g = self.grid();
        ElectromagneticField {
            electric_v_per_m: Vec3::new(
                trilinear_periodic(&self.f.ex, p, &g),
                trilinear_periodic(&self.f.ey, p, &g),
                trilinear_periodic(&self.f.ez, p, &g),
            ),
            magnetic_t: Vec3::new(
                trilinear_periodic(&self.f.bx, p, &g),
                trilinear_periodic(&self.f.by, p, &g),
                trilinear_periodic(&self.f.bz, p, &g),
            ),
        }
    }
    fn d(&self, v: &[f64], i: (usize, usize, usize), axis: usize) -> f64 {
        let c = &self.config;
        let (ix, iy, iz) = i;
        match axis {
            0 => {
                (v[fi3((ix + 1) % c.nx, iy, iz, c.nx, c.ny)]
                    - v[fi3((ix + c.nx - 1) % c.nx, iy, iz, c.nx, c.ny)])
                    / (2.0 * self.dx())
            }
            1 => {
                (v[fi3(ix, (iy + 1) % c.ny, iz, c.nx, c.ny)]
                    - v[fi3(ix, (iy + c.ny - 1) % c.ny, iz, c.nx, c.ny)])
                    / (2.0 * self.dy())
            }
            _ => {
                (v[fi3(ix, iy, (iz + 1) % c.nz, c.nx, c.ny)]
                    - v[fi3(ix, iy, (iz + c.nz - 1) % c.nz, c.nx, c.ny)])
                    / (2.0 * self.dz())
            }
        }
    }
    fn advance_magnetic(&mut self, dt: f64) {
        let c = self.config;
        let (mut bx, mut by, mut bz) = (self.f.bx.clone(), self.f.by.clone(), self.f.bz.clone());
        for iz in 0..c.nz {
            for iy in 0..c.ny {
                for ix in 0..c.nx {
                    let p = (ix, iy, iz);
                    let i = fi3(ix, iy, iz, c.nx, c.ny);
                    let cx = self.d(&self.f.ez, p, 1) - self.d(&self.f.ey, p, 2);
                    let cy = self.d(&self.f.ex, p, 2) - self.d(&self.f.ez, p, 0);
                    let cz = self.d(&self.f.ey, p, 0) - self.d(&self.f.ex, p, 1);
                    bx[i] -= dt * cx;
                    by[i] -= dt * cy;
                    bz[i] -= dt * cz;
                }
            }
        }
        self.f.bx = bx;
        self.f.by = by;
        self.f.bz = bz;
    }
    fn update_electric(&mut self, dt: f64) {
        let c = self.config;
        let (eps, mu) = (self.epsilon(), self.mu());
        let (mut ex, mut ey, mut ez) = (self.f.ex.clone(), self.f.ey.clone(), self.f.ez.clone());
        for iz in 0..c.nz {
            for iy in 0..c.ny {
                for ix in 0..c.nx {
                    let p = (ix, iy, iz);
                    let i = fi3(ix, iy, iz, c.nx, c.ny);
                    let hx = (self.d(&self.f.bz, p, 1) - self.d(&self.f.by, p, 2)) / mu;
                    let hy = (self.d(&self.f.bx, p, 2) - self.d(&self.f.bz, p, 0)) / mu;
                    let hz = (self.d(&self.f.by, p, 0) - self.d(&self.f.bx, p, 1)) / mu;
                    ex[i] += dt * (hx - self.current_x[i]) / eps;
                    ey[i] += dt * (hy - self.current_y[i]) / eps;
                    ez[i] += dt * (hz - self.current_z[i]) / eps;
                }
            }
        }
        self.f.ex = ex;
        self.f.ey = ey;
        self.f.ez = ez;
    }
    pub fn step(&mut self, steps: usize, collisions: Option<&NeutralCollisionModel>) -> Result<()> {
        let (c, g) = (self.config, self.grid());
        for _ in 0..steps {
            let old = self.particles.clone();
            self.advance_magnetic(0.5 * c.dt_s);
            for i in 0..self.particles.len() {
                let mut charged = self.particles[i].proxy();
                let f = self.gather_field(self.particles[i].position_m);
                vay_push(&mut charged, &f, c.dt_s)?;
                self.particles[i].position_m = wrap3(charged.position_m, &g);
                self.particles[i].velocity_m_per_s = charged.velocity_m_per_s;
            }
            let cur = deposit_charge_conserving_current_3d(
                &old,
                &self.particles,
                g.nx,
                g.ny,
                g.nz,
                g.lx,
                g.ly,
                g.lz,
                c.dt_s,
            )?;
            self.charge_density = cur.charge_density_new_c_per_m3;
            self.current_x = cur.current_x_a_per_m2;
            self.current_y = cur.current_y_a_per_m2;
            self.current_z = cur.current_z_a_per_m2;
            self.last_continuity_residual = cur.continuity_linf_residual;
            self.update_electric(c.dt_s);
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
                for p in &mut self.particles {
                    p.position_m = wrap3(p.position_m, &g);
                }
            }
            self.time_s += c.dt_s;
        }
        Ok(())
    }
    pub fn diagnostics(&self) -> Result<ElectromagneticPic3DDiagnostics> {
        let mut d = ElectromagneticPic3DDiagnostics {
            field_energy_j: self.f.energy(
                self.epsilon(),
                self.mu(),
                self.dx() * self.dy() * self.dz(),
            ),
            ..Default::default()
        };
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

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum GridBoundaryMode3D {
    #[default]
    Periodic,
    ElectricWall,
    AbsorbingSponge,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridBoundary3DConfig {
    pub mode: GridBoundaryMode3D,
    pub sponge_cells: usize,
    pub sponge_strength: f64,
    pub sponge_polynomial_order: f64,
}

impl Default for GridBoundary3DConfig {
    fn default() -> Self {
        Self {
            mode: GridBoundaryMode3D::Periodic,
            sponge_cells: 2,
            sponge_strength: 3.0,
            sponge_polynomial_order: 2.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum CurrentDeposition3DMode {
    #[default]
    SpectralContinuity,
    LocalFiniteVolume,
}

/// Electric walls zero all E components on boundary nodes; the sponge damps E and B by
/// `exp(−s·((L−d)/L)^p)` within `L` cells of any face.
#[allow(clippy::too_many_arguments)]
pub fn apply_electromagnetic_field_boundary_3d(
    ex: &mut [f64],
    ey: &mut [f64],
    ez: &mut [f64],
    bx: &mut [f64],
    by: &mut [f64],
    bz: &mut [f64],
    nx: usize,
    ny: usize,
    nz: usize,
    b: &GridBoundary3DConfig,
) -> Result<()> {
    let n = nx * ny * nz;
    require(
        [ex.len(), ey.len(), ez.len(), bx.len(), by.len(), bz.len()]
            .iter()
            .all(|&l| l == n),
        "invalid 3-D field boundary array size",
    )?;
    match b.mode {
        GridBoundaryMode3D::Periodic => {}
        GridBoundaryMode3D::ElectricWall => {
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        if ix == 0
                            || iy == 0
                            || iz == 0
                            || ix + 1 == nx
                            || iy + 1 == ny
                            || iz + 1 == nz
                        {
                            let i = fi3(ix, iy, iz, nx, ny);
                            ex[i] = 0.0;
                            ey[i] = 0.0;
                            ez[i] = 0.0;
                        }
                    }
                }
            }
        }
        GridBoundaryMode3D::AbsorbingSponge => {
            let cells = b.sponge_cells.max(1);
            let strength = b.sponge_strength.max(0.0);
            let order = b.sponge_polynomial_order.max(0.25);
            for iz in 0..nz {
                let ez_edge = iz.min(nz - 1 - iz);
                for iy in 0..ny {
                    let ey_edge = iy.min(ny - 1 - iy);
                    for ix in 0..nx {
                        let edge = ix.min(nx - 1 - ix).min(ey_edge.min(ez_edge));
                        if edge < cells {
                            let s = (cells - edge) as f64 / cells as f64;
                            let factor = (-strength * s.powf(order)).exp();
                            let i = fi3(ix, iy, iz, nx, ny);
                            for f in [&mut *ex, &mut *ey, &mut *ez, &mut *bx, &mut *by, &mut *bz] {
                                f[i] *= factor;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StaggeredElectromagneticPic3DConfig {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub length_x_m: f64,
    pub length_y_m: f64,
    pub length_z_m: f64,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub dt_s: f64,
    pub field_boundary: GridBoundary3DConfig,
    pub particle_boundary: ParticleWallMode,
    pub periodic_particles: bool,
    pub remove_absorbed_particles: bool,
    pub enable_secondary_yield_report: bool,
    pub secondary_model: SecondaryEmissionModel,
    pub current_deposition: CurrentDeposition3DMode,
    pub sort_particles_by_cell: bool,
    pub particle_sort_interval: usize,
}

impl Default for StaggeredElectromagneticPic3DConfig {
    fn default() -> Self {
        Self {
            nx: 10,
            ny: 10,
            nz: 10,
            length_x_m: 1.0,
            length_y_m: 1.0,
            length_z_m: 1.0,
            relative_permittivity: 1.0,
            relative_permeability: 1.0,
            dt_s: 1.0e-12,
            field_boundary: GridBoundary3DConfig::default(),
            particle_boundary: ParticleWallMode::SpecularReflect,
            periodic_particles: true,
            remove_absorbed_particles: true,
            enable_secondary_yield_report: false,
            secondary_model: SecondaryEmissionModel::default(),
            current_deposition: CurrentDeposition3DMode::SpectralContinuity,
            sort_particles_by_cell: false,
            particle_sort_interval: 1,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StaggeredElectromagneticPic3DDiagnostics {
    pub field_energy_j: f64,
    pub particle_kinetic_energy_j: f64,
    pub total_energy_j: f64,
    pub charge_continuity_linf_residual: f64,
    pub particle_count: usize,
    pub absorbed_particles: usize,
    pub secondary_macro_weight: f64,
    pub collision_statistics: CollisionStatistics,
    pub particle_sort_passes: usize,
    pub occupied_particle_cells: usize,
}

/// Yee-staggered 3-D EM PIC (CFL guard 0.98) with field and particle boundaries, spectral or
/// local finite-volume current deposition and optional periodic cell sorting.
#[derive(Clone, Debug)]
pub struct StaggeredElectromagneticPic3D {
    config: StaggeredElectromagneticPic3DConfig,
    particles: Vec<PicParticle3D>,
    charge_density: Vec<f64>,
    current_x: Vec<f64>,
    current_y: Vec<f64>,
    current_z: Vec<f64>,
    f: Fields3,
    time_s: f64,
    last_continuity_residual: f64,
    absorbed_particles: usize,
    secondary_macro_weight: f64,
    accumulated_collisions: CollisionStatistics,
    particle_sort_passes: usize,
    occupied_particle_cells: usize,
    step_counter: usize,
}

impl StaggeredElectromagneticPic3D {
    pub fn new(config: StaggeredElectromagneticPic3DConfig) -> Result<Self> {
        let s = Self {
            config,
            particles: Vec::new(),
            charge_density: Vec::new(),
            current_x: Vec::new(),
            current_y: Vec::new(),
            current_z: Vec::new(),
            f: Fields3::default(),
            time_s: 0.0,
            last_continuity_residual: 0.0,
            absorbed_particles: 0,
            secondary_macro_weight: 0.0,
            accumulated_collisions: CollisionStatistics::default(),
            particle_sort_passes: 0,
            occupied_particle_cells: 0,
            step_counter: 0,
        };
        let g = s.grid();
        require(
            g.valid()
                && config.relative_permittivity > 0.0
                && config.relative_permeability > 0.0
                && config.dt_s > 0.0
                && config.dt_s.is_finite(),
            "invalid staggered electromagnetic PIC 3-D configuration",
        )?;
        require(
            config.dt_s * s.light() * s.inv_cfl() < 0.98,
            "staggered 3-D EM-PIC time step violates Yee CFL guard",
        )?;
        let n = g.n();
        Ok(Self {
            charge_density: vec![0.0; n],
            current_x: vec![0.0; n],
            current_y: vec![0.0; n],
            current_z: vec![0.0; n],
            f: Fields3::zeros(n),
            ..s
        })
    }
    fields3_getters!();
    fn periodic_fields(&self) -> bool {
        self.config.field_boundary.mode == GridBoundaryMode3D::Periodic
    }
    fn place(&self, p: Vec3) -> Vec3 {
        let c = &self.config;
        if c.periodic_particles || self.periodic_fields() {
            wrap3(p, &self.grid())
        } else {
            Vec3::new(
                p.x.clamp(0.0, c.length_x_m),
                p.y.clamp(0.0, c.length_y_m),
                p.z.clamp(0.0, c.length_z_m),
            )
        }
    }
    pub fn set_particles(&mut self, mut particles: Vec<PicParticle3D>) -> Result<()> {
        for p in &mut particles {
            require(p.valid(), "invalid staggered 3-D PIC particle")?;
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
        let n = self.grid().n();
        self.f.set(
            [ex, ey, ez, bx, by, bz],
            n,
            "staggered 3-D field arrays must match grid size",
        )?;
        self.apply_field_boundary()
    }
    pub fn initialize_z_polarized_mode(
        &mut self,
        amplitude: f64,
        mx: usize,
        my: usize,
        mz: usize,
    ) -> Result<()> {
        self.f = Fields3::zeros(self.grid().n());
        require(
            mx != 0 || my != 0 || mz != 0,
            "staggered 3-D EM-PIC mode cannot be DC",
        )?;
        let c = self.config;
        let kx = 2.0 * PI * mx as f64 / c.length_x_m;
        let ky = 2.0 * PI * my as f64 / c.length_y_m;
        let kz = 2.0 * PI * mz as f64 / c.length_z_m;
        let k = (kx * kx + ky * ky + kz * kz).sqrt();
        let (dx, dy, dz, light) = (self.dx(), self.dy(), self.dz(), self.light());
        for iz in 0..c.nz {
            let zh = (iz as f64 + 0.5) * dz;
            for iy in 0..c.ny {
                let (y0, yh) = (iy as f64 * dy, (iy as f64 + 0.5) * dy);
                for ix in 0..c.nx {
                    let (x0, xh) = (ix as f64 * dx, (ix as f64 + 0.5) * dx);
                    let i = fi3(ix, iy, iz, c.nx, c.ny);
                    self.f.ez[i] = amplitude * (kx * x0 + ky * y0 + kz * zh).sin();
                    self.f.bx[i] =
                        (ky / k) * (amplitude / light) * (kx * x0 + ky * yh + kz * zh).sin();
                    self.f.by[i] =
                        -(kx / k) * (amplitude / light) * (kx * xh + ky * y0 + kz * zh).sin();
                }
            }
        }
        self.apply_field_boundary()
    }
    fn sample(&self, values: &[f64], p: Vec3, o: [f64; 3]) -> f64 {
        let c = &self.config;
        let periodic = self.periodic_fields();
        let wrap = |v: f64, l: f64| {
            if periodic {
                wrap_periodic(v, l)
            } else {
                v.max(0.0).min(next_down(l))
            }
        };
        let x = wrap(p.x - o[0] * self.dx(), c.length_x_m);
        let y = wrap(p.y - o[1] * self.dy(), c.length_y_m);
        let z = wrap(p.z - o[2] * self.dz(), c.length_z_m);
        let u = [x / self.dx(), y / self.dy(), z / self.dz()];
        let n = [c.nx, c.ny, c.nz];
        let mut i0 = [0; 3];
        let mut i1 = [0; 3];
        let mut f = [0.0; 3];
        for a in 0..3 {
            let fl = u[a].floor();
            i0[a] = (fl as usize) % n[a];
            i1[a] = if periodic {
                (i0[a] + 1) % n[a]
            } else {
                (i0[a] + 1).min(n[a] - 1)
            };
            f[a] = u[a] - fl;
        }
        Cic3::new(i0, i1, f).gather(values, c.nx, c.ny)
    }
    pub fn gather_field(&self, p: Vec3) -> ElectromagneticField {
        ElectromagneticField {
            electric_v_per_m: Vec3::new(
                self.sample(&self.f.ex, p, [0.5, 0.0, 0.0]),
                self.sample(&self.f.ey, p, [0.0, 0.5, 0.0]),
                self.sample(&self.f.ez, p, [0.0, 0.0, 0.5]),
            ),
            magnetic_t: Vec3::new(
                self.sample(&self.f.bx, p, [0.0, 0.5, 0.5]),
                self.sample(&self.f.by, p, [0.5, 0.0, 0.5]),
                self.sample(&self.f.bz, p, [0.5, 0.5, 0.0]),
            ),
        }
    }
    fn neighbor(&self, i: usize, n: usize, forward: bool) -> usize {
        let periodic = self.periodic_fields();
        if forward {
            if i + 1 < n {
                i + 1
            } else if periodic {
                0
            } else {
                i
            }
        } else if i > 0 {
            i - 1
        } else if periodic {
            n - 1
        } else {
            i
        }
    }
    fn diff(&self, v: &[f64], p: (usize, usize, usize), axis: usize, forward: bool) -> f64 {
        let c = &self.config;
        let (ix, iy, iz) = p;
        let here = v[fi3(ix, iy, iz, c.nx, c.ny)];
        let (other, h) = match axis {
            0 => (
                v[fi3(self.neighbor(ix, c.nx, forward), iy, iz, c.nx, c.ny)],
                self.dx(),
            ),
            1 => (
                v[fi3(ix, self.neighbor(iy, c.ny, forward), iz, c.nx, c.ny)],
                self.dy(),
            ),
            _ => (
                v[fi3(ix, iy, self.neighbor(iz, c.nz, forward), c.nx, c.ny)],
                self.dz(),
            ),
        };
        if forward {
            (other - here) / h
        } else {
            (here - other) / h
        }
    }
    fn apply_field_boundary(&mut self) -> Result<()> {
        let c = self.config;
        let f = &mut self.f;
        apply_electromagnetic_field_boundary_3d(
            &mut f.ex,
            &mut f.ey,
            &mut f.ez,
            &mut f.bx,
            &mut f.by,
            &mut f.bz,
            c.nx,
            c.ny,
            c.nz,
            &c.field_boundary,
        )
    }
    fn advance_magnetic(&mut self, dt: f64) -> Result<()> {
        let c = self.config;
        let (mut bx, mut by, mut bz) = (self.f.bx.clone(), self.f.by.clone(), self.f.bz.clone());
        for iz in 0..c.nz {
            for iy in 0..c.ny {
                for ix in 0..c.nx {
                    let p = (ix, iy, iz);
                    let i = fi3(ix, iy, iz, c.nx, c.ny);
                    bx[i] -= dt
                        * (self.diff(&self.f.ez, p, 1, true) - self.diff(&self.f.ey, p, 2, true));
                    by[i] -= dt
                        * (self.diff(&self.f.ex, p, 2, true) - self.diff(&self.f.ez, p, 0, true));
                    bz[i] -= dt
                        * (self.diff(&self.f.ey, p, 0, true) - self.diff(&self.f.ex, p, 1, true));
                }
            }
        }
        self.f.bx = bx;
        self.f.by = by;
        self.f.bz = bz;
        self.apply_field_boundary()
    }
    fn update_electric(&mut self, dt: f64) -> Result<()> {
        let c = self.config;
        let (eps, mu) = (self.epsilon(), self.mu());
        let (mut ex, mut ey, mut ez) = (self.f.ex.clone(), self.f.ey.clone(), self.f.ez.clone());
        for iz in 0..c.nz {
            for iy in 0..c.ny {
                for ix in 0..c.nx {
                    let p = (ix, iy, iz);
                    let i = fi3(ix, iy, iz, c.nx, c.ny);
                    let hx = (self.diff(&self.f.bz, p, 1, false)
                        - self.diff(&self.f.by, p, 2, false))
                        / mu;
                    let hy = (self.diff(&self.f.bx, p, 2, false)
                        - self.diff(&self.f.bz, p, 0, false))
                        / mu;
                    let hz = (self.diff(&self.f.by, p, 0, false)
                        - self.diff(&self.f.bx, p, 1, false))
                        / mu;
                    ex[i] += dt * (hx - self.current_x[i]) / eps;
                    ey[i] += dt * (hy - self.current_y[i]) / eps;
                    ez[i] += dt * (hz - self.current_z[i]) / eps;
                }
            }
        }
        self.f.ex = ex;
        self.f.ey = ey;
        self.f.ez = ez;
        self.apply_field_boundary()
    }
    fn apply_particle_boundary(&mut self) -> Result<()> {
        let c = self.config;
        if c.periodic_particles || self.periodic_fields() {
            for i in 0..self.particles.len() {
                self.particles[i].position_m = self.place(self.particles[i].position_m);
            }
            return Ok(());
        }
        let b = AxisAlignedParticleBox {
            minimum_m: Vec3::default(),
            maximum_m: Vec3::new(c.length_x_m, c.length_y_m, c.length_z_m),
        };
        let secondary = if c.enable_secondary_yield_report {
            Some(&c.secondary_model)
        } else {
            None
        };
        let mut survivors = Vec::with_capacity(self.particles.len());
        for mut p in self.particles.iter().copied() {
            let mut charged = p.proxy();
            let r = apply_particle_box_boundary(&mut charged, &b, c.particle_boundary, secondary)?;
            if r.impacted {
                self.absorbed_particles += 1;
                self.secondary_macro_weight += r.secondary_macro_weight;
            }
            p.position_m = charged.position_m;
            p.velocity_m_per_s = charged.velocity_m_per_s;
            if r.alive || !c.remove_absorbed_particles {
                survivors.push(p);
            }
        }
        self.particles = survivors;
        Ok(())
    }
    fn maybe_sort(&mut self) -> Result<()> {
        self.step_counter += 1;
        let c = self.config;
        if !c.sort_particles_by_cell
            || c.particle_sort_interval == 0
            || self.step_counter % c.particle_sort_interval != 0
        {
            return Ok(());
        }
        let sorted = sort_particles_by_cell_3d(
            &self.particles,
            c.nx,
            c.ny,
            c.nz,
            c.length_x_m,
            c.length_y_m,
            c.length_z_m,
        )?;
        self.occupied_particle_cells = sorted.occupied_cells;
        self.particles = sorted.sorted_particles;
        self.particle_sort_passes += 1;
        Ok(())
    }
    pub fn step(&mut self, steps: usize, collisions: Option<&NeutralCollisionModel>) -> Result<()> {
        let (c, g) = (self.config, self.grid());
        for _ in 0..steps {
            let old = self.particles.clone();
            self.advance_magnetic(0.5 * c.dt_s)?;
            for i in 0..self.particles.len() {
                let mut charged = self.particles[i].proxy();
                let f = self.gather_field(self.particles[i].position_m);
                vay_push(&mut charged, &f, c.dt_s)?;
                self.particles[i].position_m = charged.position_m;
                self.particles[i].velocity_m_per_s = charged.velocity_m_per_s;
            }
            self.apply_particle_boundary()?;
            if old.len() == self.particles.len() {
                let cur = match c.current_deposition {
                    CurrentDeposition3DMode::LocalFiniteVolume => {
                        deposit_charge_conserving_current_3d_local(
                            &old,
                            &self.particles,
                            g.nx,
                            g.ny,
                            g.nz,
                            g.lx,
                            g.ly,
                            g.lz,
                            c.dt_s,
                        )?
                    }
                    CurrentDeposition3DMode::SpectralContinuity => {
                        deposit_charge_conserving_current_3d(
                            &old,
                            &self.particles,
                            g.nx,
                            g.ny,
                            g.nz,
                            g.lx,
                            g.ly,
                            g.lz,
                            c.dt_s,
                        )?
                    }
                };
                self.charge_density = cur.charge_density_new_c_per_m3;
                self.current_x = cur.current_x_a_per_m2;
                self.current_y = cur.current_y_a_per_m2;
                self.current_z = cur.current_z_a_per_m2;
                self.last_continuity_residual = cur.continuity_linf_residual;
            } else {
                self.charge_density = deposit_cic_charge_density_3d(
                    &self.particles,
                    g.nx,
                    g.ny,
                    g.nz,
                    g.lx,
                    g.ly,
                    g.lz,
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
            self.maybe_sort()?;
        }
        Ok(())
    }
    pub fn diagnostics(&self) -> Result<StaggeredElectromagneticPic3DDiagnostics> {
        let mut d = StaggeredElectromagneticPic3DDiagnostics {
            field_energy_j: self.f.energy(
                self.epsilon(),
                self.mu(),
                self.dx() * self.dy() * self.dz(),
            ),
            ..Default::default()
        };
        for p in &self.particles {
            d.particle_kinetic_energy_j += p.weight * kinetic_energy_j(&p.proxy())?;
        }
        d.total_energy_j = d.field_energy_j + d.particle_kinetic_energy_j;
        d.charge_continuity_linf_residual = self.last_continuity_residual;
        d.particle_count = self.particles.len();
        d.absorbed_particles = self.absorbed_particles;
        d.secondary_macro_weight = self.secondary_macro_weight;
        d.collision_statistics = self.accumulated_collisions;
        d.particle_sort_passes = self.particle_sort_passes;
        d.occupied_particle_cells = self.occupied_particle_cells;
        Ok(d)
    }
}

/// `std::nextafter(x, 0.0)` for positive finite `x`.
fn next_down(x: f64) -> f64 {
    if x > 0.0 && x.is_finite() {
        f64::from_bits(x.to_bits() - 1)
    } else {
        x
    }
}
