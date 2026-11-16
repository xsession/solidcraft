//! Periodic 3-D VOF free-surface baseline: classification, first-order conservative advection,
//! normals, curvature and CSF acceleration, ported from upstream `cfd/solvers/lbm/free_surface`
//! (MIT, cfd_solvers contributors).
use super::ParticleVector3;
use crate::{require, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FreeSurfacePhase {
    Gas,
    Interface,
    Fluid,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceField3D {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub fill: Vec<f64>,
    pub phase: Vec<FreeSurfacePhase>,
}

impl FreeSurfaceField3D {
    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.ny + y) * self.nx + x
    }
    fn wrap(i: i64, n: usize) -> usize {
        i.rem_euclid(n as i64) as usize
    }
    /// Central-difference fill gradient at a cell (periodic).
    fn gradient(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let (xi, yi, zi) = (x as i64, y as i64, z as i64);
        let f = &self.fill;
        [
            0.5 * (f[self.idx(Self::wrap(xi + 1, self.nx), y, z)]
                - f[self.idx(Self::wrap(xi - 1, self.nx), y, z)]),
            0.5 * (f[self.idx(x, Self::wrap(yi + 1, self.ny), z)]
                - f[self.idx(x, Self::wrap(yi - 1, self.ny), z)]),
            0.5 * (f[self.idx(x, y, Self::wrap(zi + 1, self.nz))]
                - f[self.idx(x, y, Self::wrap(zi - 1, self.nz))]),
        ]
    }
    fn classify(&mut self, eps: f64) {
        for (f, p) in self.fill.iter_mut().zip(&mut self.phase) {
            *f = f.clamp(0.0, 1.0);
            *p = if *f <= eps {
                FreeSurfacePhase::Gas
            } else if *f >= 1.0 - eps {
                FreeSurfacePhase::Fluid
            } else {
                FreeSurfacePhase::Interface
            };
        }
    }
}

pub fn make_free_surface_field(
    nx: usize,
    ny: usize,
    nz: usize,
    fill: &[f64],
    interface_epsilon: f64,
) -> Result<FreeSurfaceField3D> {
    require(
        nx >= 2
            && ny >= 2
            && nz >= 2
            && Some(fill.len()) == nx.checked_mul(ny).and_then(|v| v.checked_mul(nz))
            && (0.0..0.5).contains(&interface_epsilon),
        "invalid free-surface field",
    )?;
    let mut s = FreeSurfaceField3D {
        nx,
        ny,
        nz,
        fill: fill.to_vec(),
        phase: vec![FreeSurfacePhase::Gas; fill.len()],
    };
    s.classify(interface_epsilon);
    Ok(s)
}

pub fn free_surface_volume(s: &FreeSurfaceField3D) -> f64 {
    s.fill.iter().sum()
}

/// Upwind face fluxes with face velocity `0.5 (u_a + u_b)`; requires `dt * max(|ux|+|uy|+|uz|) <= 1`.
pub fn advect_free_surface_vof(
    s: &mut FreeSurfaceField3D,
    ux: &[f32],
    uy: &[f32],
    uz: &[f32],
    dt: f64,
    interface_epsilon: f64,
) -> Result<()> {
    let n = s.fill.len();
    require(
        ux.len() == n && uy.len() == n && uz.len() == n,
        "free-surface velocity field size mismatch",
    )?;
    require(dt > 0.0, "free-surface dt must be positive")?;
    let vmax = (0..n).fold(0.0f64, |m, i| {
        m.max(f64::from(ux[i]).abs() + f64::from(uy[i]).abs() + f64::from(uz[i]).abs())
    });
    require(
        dt * vmax <= 1.0 + 1e-12,
        "free-surface CFL exceeds one cell per step",
    )?;
    let mut delta = vec![0.0; n];
    let flux = |fill: &[f64], a: usize, b: usize, vel: f64| {
        let face = 0.5 * vel;
        if face >= 0.0 {
            face * fill[a]
        } else {
            face * fill[b]
        }
    };
    for z in 0..s.nz {
        for y in 0..s.ny {
            for x in 0..s.nx {
                let a = s.idx(x, y, z);
                let xp = s.idx((x + 1) % s.nx, y, z);
                let yp = s.idx(x, (y + 1) % s.ny, z);
                let zp = s.idx(x, y, (z + 1) % s.nz);
                let fx = dt * flux(&s.fill, a, xp, f64::from(ux[a] + ux[xp]));
                let fy = dt * flux(&s.fill, a, yp, f64::from(uy[a] + uy[yp]));
                let fz = dt * flux(&s.fill, a, zp, f64::from(uz[a] + uz[zp]));
                delta[a] -= fx + fy + fz;
                delta[xp] += fx;
                delta[yp] += fy;
                delta[zp] += fz;
            }
        }
    }
    for (f, d) in s.fill.iter_mut().zip(&delta) {
        *f += d;
    }
    s.classify(interface_epsilon);
    Ok(())
}

/// Unit fill-gradient normals (zero where the gradient vanishes).
pub fn free_surface_normals(s: &FreeSurfaceField3D) -> Vec<ParticleVector3> {
    let mut out = vec![ParticleVector3::default(); s.fill.len()];
    for z in 0..s.nz {
        for y in 0..s.ny {
            for x in 0..s.nx {
                let [gx, gy, gz] = s.gradient(x, y, z);
                let m = (gx * gx + gy * gy + gz * gz).sqrt();
                if m > 1e-14 {
                    out[s.idx(x, y, z)] = ParticleVector3 {
                        x: gx / m,
                        y: gy / m,
                        z: gz / m,
                    };
                }
            }
        }
    }
    out
}

/// `kappa = -div(n)` by central differences.
pub fn free_surface_curvature(s: &FreeSurfaceField3D) -> Vec<f64> {
    let normals = free_surface_normals(s);
    let mut k = vec![0.0; s.fill.len()];
    let w = FreeSurfaceField3D::wrap;
    for z in 0..s.nz {
        for y in 0..s.ny {
            for x in 0..s.nx {
                let (xi, yi, zi) = (x as i64, y as i64, z as i64);
                let xp = s.idx(w(xi + 1, s.nx), y, z);
                let xm = s.idx(w(xi - 1, s.nx), y, z);
                let yp = s.idx(x, w(yi + 1, s.ny), z);
                let ym = s.idx(x, w(yi - 1, s.ny), z);
                let zp = s.idx(x, y, w(zi + 1, s.nz));
                let zm = s.idx(x, y, w(zi - 1, s.nz));
                k[s.idx(x, y, z)] = -0.5
                    * ((normals[xp].x - normals[xm].x)
                        + (normals[yp].y - normals[ym].y)
                        + (normals[zp].z - normals[zm].z));
            }
        }
    }
    k
}

/// Continuum-surface-force acceleration `sigma kappa grad(phi) / rho`.
pub fn free_surface_surface_tension_acceleration(
    s: &FreeSurfaceField3D,
    surface_tension: f64,
    density: f64,
) -> Result<Vec<ParticleVector3>> {
    require(
        surface_tension >= 0.0 && density > 0.0,
        "invalid free-surface material properties",
    )?;
    let k = free_surface_curvature(s);
    let mut a = vec![ParticleVector3::default(); s.fill.len()];
    for z in 0..s.nz {
        for y in 0..s.ny {
            for x in 0..s.nx {
                let n = s.idx(x, y, z);
                let [gx, gy, gz] = s.gradient(x, y, z);
                let c = surface_tension * k[n] / density;
                a[n] = ParticleVector3 {
                    x: c * gx,
                    y: c * gy,
                    z: c * gz,
                };
            }
        }
    }
    Ok(a)
}
