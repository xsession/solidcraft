//! 3-D Cartesian Yee FDTD with PEC/PMC walls and diagonal anisotropy, its legacy VTK writer and
//! the optional HDF5/XDMF seam, ported from upstream `cfd/solvers/fdtd/{maxwell3d, vtk_io,
//! hdf5_xdmf}` (MIT, cfd_solvers contributors).
use super::{EPSILON0, MU0};
use crate::{lbm::format_general, require, Result, SolverError};
use std::{fmt::Write as _, path::Path};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Boundary3D {
    #[default]
    Pec,
    Pmc,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Maxwell3DConfig {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub courant: f64,
    pub epsilon_r: f64,
    pub mu_r: f64,
    pub boundary: Boundary3D,
    /// Optional diagonal anisotropy; non-positive entries fall back to `epsilon_r`.
    pub epsilon_rx: f64,
    pub epsilon_ry: f64,
    pub epsilon_rz: f64,
}

impl Default for Maxwell3DConfig {
    fn default() -> Self {
        Self {
            nx: 32,
            ny: 32,
            nz: 32,
            dx: 1e-3,
            dy: 1e-3,
            dz: 1e-3,
            courant: 0.95,
            epsilon_r: 1.0,
            mu_r: 1.0,
            boundary: Boundary3D::Pec,
            epsilon_rx: 0.0,
            epsilon_ry: 0.0,
            epsilon_rz: 0.0,
        }
    }
}

/// Yee scheme on a common `nx * ny * nz` lattice (forward-difference H, backward-difference E).
#[derive(Clone, Debug)]
pub struct Maxwell3D {
    c: Maxwell3DConfig,
    dt: f64,
    eps: [f64; 3],
    mu: f64,
    ex: Vec<f64>,
    ey: Vec<f64>,
    ez: Vec<f64>,
    hx: Vec<f64>,
    hy: Vec<f64>,
    hz: Vec<f64>,
}

impl Maxwell3D {
    pub fn new(c: Maxwell3DConfig) -> Result<Self> {
        require(
            c.nx >= 4
                && c.ny >= 4
                && c.nz >= 4
                && c.dx > 0.0
                && c.dy > 0.0
                && c.dz > 0.0
                && c.courant > 0.0
                && c.courant < 1.0
                && c.epsilon_r > 0.0
                && c.mu_r > 0.0,
            "invalid Maxwell3D configuration",
        )?;
        require(
            [
                c.dx,
                c.dy,
                c.dz,
                c.epsilon_r,
                c.mu_r,
                c.epsilon_rx,
                c.epsilon_ry,
                c.epsilon_rz,
            ]
            .iter()
            .all(|v| v.is_finite()),
            "invalid Maxwell3D configuration",
        )?;
        let n =
            c.nx.checked_mul(c.ny)
                .and_then(|v| v.checked_mul(c.nz))
                .ok_or_else(|| SolverError::InvalidInput("Maxwell3D grid too large".into()))?;
        let pick = |v: f64| EPSILON0 * if v > 0.0 { v } else { c.epsilon_r };
        let eps = [pick(c.epsilon_rx), pick(c.epsilon_ry), pick(c.epsilon_rz)];
        let mu = MU0 * c.mu_r;
        let wave = 1.0 / (eps[0].min(eps[1]).min(eps[2]) * mu).sqrt();
        let dt = c.courant
            / (wave * (1.0 / (c.dx * c.dx) + 1.0 / (c.dy * c.dy) + 1.0 / (c.dz * c.dz)).sqrt());
        Ok(Self {
            c,
            dt,
            eps,
            mu,
            ex: vec![0.0; n],
            ey: vec![0.0; n],
            ez: vec![0.0; n],
            hx: vec![0.0; n],
            hy: vec![0.0; n],
            hz: vec![0.0; n],
        })
    }
    fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.c.ny + j) * self.c.nx + i
    }
    /// Upstream defaults `amplitude = 1`, `width_fraction = 0.12`.
    pub fn initialize_gaussian_ez(&mut self, amplitude: f64, width_fraction: f64) -> Result<()> {
        require(width_fraction > 0.0, "Gaussian width must be positive")?;
        let (nx, ny, nz) = (self.c.nx, self.c.ny, self.c.nz);
        let (cx, cy, cz) = (
            0.5 * (nx - 1) as f64,
            0.5 * (ny - 1) as f64,
            0.5 * (nz - 1) as f64,
        );
        let w = width_fraction * nx.min(ny).min(nz) as f64;
        for (n, v) in self.ez.iter_mut().enumerate() {
            let (i, j, k) = (n % nx, (n / nx) % ny, n / (nx * ny));
            let r2 = ((i as f64 - cx) * (i as f64 - cx)
                + (j as f64 - cy) * (j as f64 - cy)
                + (k as f64 - cz) * (k as f64 - cz))
                / (w * w);
            *v = amplitude * (-0.5 * r2).exp();
        }
        self.enforce_boundary();
        Ok(())
    }
    pub fn add_soft_ez_source(&mut self, i: usize, j: usize, k: usize, value: f64) -> Result<()> {
        require(
            i < self.c.nx && j < self.c.ny && k < self.c.nz && value.is_finite(),
            "invalid Maxwell3D soft source",
        )?;
        let n = self.idx(i, j, k);
        self.ez[n] += value;
        Ok(())
    }
    pub fn step(&mut self, count: usize) {
        let (nx, ny, nz) = (self.c.nx, self.c.ny, self.c.nz);
        let (dx, dy, dz) = (self.c.dx, self.c.dy, self.c.dz);
        let hdt = self.dt / self.mu;
        let (ex_dt, ey_dt, ez_dt) = (
            self.dt / self.eps[0],
            self.dt / self.eps[1],
            self.dt / self.eps[2],
        );
        let id = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;
        for _ in 0..count {
            for k in 0..nz - 1 {
                for j in 0..ny - 1 {
                    for i in 0..nx {
                        let n = id(i, j, k);
                        self.hx[n] -= hdt
                            * ((self.ez[id(i, j + 1, k)] - self.ez[n]) / dy
                                - (self.ey[id(i, j, k + 1)] - self.ey[n]) / dz);
                    }
                }
            }
            for k in 0..nz - 1 {
                for j in 0..ny {
                    for i in 0..nx - 1 {
                        let n = id(i, j, k);
                        self.hy[n] -= hdt
                            * ((self.ex[id(i, j, k + 1)] - self.ex[n]) / dz
                                - (self.ez[id(i + 1, j, k)] - self.ez[n]) / dx);
                    }
                }
            }
            for k in 0..nz {
                for j in 0..ny - 1 {
                    for i in 0..nx - 1 {
                        let n = id(i, j, k);
                        self.hz[n] -= hdt
                            * ((self.ey[id(i + 1, j, k)] - self.ey[n]) / dx
                                - (self.ex[id(i, j + 1, k)] - self.ex[n]) / dy);
                    }
                }
            }
            for k in 1..nz - 1 {
                for j in 1..ny - 1 {
                    for i in 1..nx - 1 {
                        let n = id(i, j, k);
                        self.ex[n] += ex_dt
                            * ((self.hz[n] - self.hz[id(i, j - 1, k)]) / dy
                                - (self.hy[n] - self.hy[id(i, j, k - 1)]) / dz);
                        self.ey[n] += ey_dt
                            * ((self.hx[n] - self.hx[id(i, j, k - 1)]) / dz
                                - (self.hz[n] - self.hz[id(i - 1, j, k)]) / dx);
                        self.ez[n] += ez_dt
                            * ((self.hy[n] - self.hy[id(i - 1, j, k)]) / dx
                                - (self.hx[n] - self.hx[id(i, j - 1, k)]) / dy);
                    }
                }
            }
            self.enforce_boundary();
        }
    }
    fn enforce_boundary(&mut self) {
        let (nx, ny, nz) = (self.c.nx, self.c.ny, self.c.nz);
        for n in 0..nx * ny * nz {
            let (i, j, k) = (n % nx, (n / nx) % ny, n / (nx * ny));
            if i == 0 || i + 1 == nx || j == 0 || j + 1 == ny || k == 0 || k + 1 == nz {
                if self.c.boundary == Boundary3D::Pec {
                    (self.ex[n], self.ey[n], self.ez[n]) = (0.0, 0.0, 0.0);
                } else {
                    (self.hx[n], self.hy[n], self.hz[n]) = (0.0, 0.0, 0.0);
                }
            }
        }
    }
    pub fn energy(&self) -> f64 {
        let mut e = 0.0;
        for n in 0..self.cell_count() {
            e += 0.5
                * (self.eps[0] * self.ex[n] * self.ex[n]
                    + self.eps[1] * self.ey[n] * self.ey[n]
                    + self.eps[2] * self.ez[n] * self.ez[n]
                    + self.mu
                        * (self.hx[n] * self.hx[n]
                            + self.hy[n] * self.hy[n]
                            + self.hz[n] * self.hz[n]))
                * self.c.dx
                * self.c.dy
                * self.c.dz;
        }
        e
    }
    pub fn max_field(&self) -> f64 {
        (0..self.cell_count()).fold(0.0_f64, |m, n| {
            m.max(
                (self.ex[n] * self.ex[n] + self.ey[n] * self.ey[n] + self.ez[n] * self.ez[n])
                    .sqrt(),
            )
        })
    }
    pub fn dt(&self) -> f64 {
        self.dt
    }
    pub fn config(&self) -> Maxwell3DConfig {
        self.c
    }
    pub fn cell_count(&self) -> usize {
        self.c.nx * self.c.ny * self.c.nz
    }
    pub fn nx(&self) -> usize {
        self.c.nx
    }
    pub fn ny(&self) -> usize {
        self.c.ny
    }
    pub fn nz(&self) -> usize {
        self.c.nz
    }
    pub fn dx(&self) -> f64 {
        self.c.dx
    }
    pub fn dy(&self) -> f64 {
        self.c.dy
    }
    pub fn dz(&self) -> f64 {
        self.c.dz
    }
    pub fn ex(&self) -> &[f64] {
        &self.ex
    }
    pub fn ey(&self) -> &[f64] {
        &self.ey
    }
    pub fn ez(&self) -> &[f64] {
        &self.ez
    }
    pub fn hx(&self) -> &[f64] {
        &self.hx
    }
    pub fn hy(&self) -> &[f64] {
        &self.hy
    }
    pub fn hz(&self) -> &[f64] {
        &self.hz
    }
}

/// Legacy ASCII VTK structured points with `E` and `H` vectors (17 significant digits), byte
/// for byte as upstream.
pub fn maxwell3d_vtk_ascii(s: &Maxwell3D) -> String {
    let g = |v: f64| format_general(v, 17);
    let mut o = format!(
        "# vtk DataFile Version 3.0\ncfd_solvers Maxwell3D fields\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS {} {} {}\nORIGIN 0 0 0\nSPACING {} {} {}\nPOINT_DATA {}\nVECTORS E double\n",
        s.nx(),
        s.ny(),
        s.nz(),
        g(s.dx()),
        g(s.dy()),
        g(s.dz()),
        s.cell_count()
    );
    for n in 0..s.cell_count() {
        let _ = writeln!(o, "{} {} {}", g(s.ex[n]), g(s.ey[n]), g(s.ez[n]));
    }
    o.push_str("VECTORS H double\n");
    for n in 0..s.cell_count() {
        let _ = writeln!(o, "{} {} {}", g(s.hx[n]), g(s.hy[n]), g(s.hz[n]));
    }
    o
}

pub fn write_maxwell3d_vtk_ascii(s: &Maxwell3D, path: impl AsRef<Path>) -> Result<()> {
    let path = path.as_ref();
    std::fs::write(path, maxwell3d_vtk_ascii(s)).map_err(|e| {
        SolverError::Io(format!(
            "failed to write VTK output file {}: {e}",
            path.display()
        ))
    })
}

/// HDF5 is not linked in this offline build, matching upstream builds without `CFD_ENABLE_HDF5`.
pub fn hdf5_xdmf_output_available() -> bool {
    false
}

/// Always `Unsupported` here (no HDF5 library); use [`write_maxwell3d_vtk_ascii`].
pub fn write_maxwell3d_hdf5_xdmf(
    _solver: &Maxwell3D,
    _xdmf_path: impl AsRef<Path>,
    _hdf5_path: impl AsRef<Path>,
) -> Result<()> {
    Err(SolverError::Unsupported(
        "FDTD HDF5/XDMF output requires CFD_ENABLE_HDF5 and an HDF5 library".into(),
    ))
}
