//! Single-lattice Esoteric-Pull, two-lattice pull, 16-bit compressed-population and
//! precision-generic pull solvers for D2Q9/D3Q19/D3Q27, ported from upstream
//! `cfd/solvers/lbm/{esoteric_pull, one_step_pull, compressed_pull, precision_pull}`
//! (MIT, cfd_solvers contributors). Serial CPU reference; populations are stored
//! structure-of-arrays as `f[direction * cells + cell]`.
// Index loops mirror the upstream lattice stencils direction by direction.
#![allow(clippy::needless_range_loop)]
use super::LatticeDescriptor;
use crate::{require, Result, SolverError};
use std::marker::PhantomData;

/// Largest supported velocity set (D3Q27).
const MAX_Q: usize = 27;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InPlaceLbmConfig {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub tau: f32,
    pub acceleration_x: f32,
    pub acceleration_y: f32,
    pub acceleration_z: f32,
    pub smagorinsky_les: bool,
}

impl Default for InPlaceLbmConfig {
    fn default() -> Self {
        Self {
            nx: 128,
            ny: 128,
            nz: 128,
            tau: 0.6,
            acceleration_x: 0.0,
            acceleration_y: 0.0,
            acceleration_z: 0.0,
            smagorinsky_les: false,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MacroscopicFields {
    pub rho: Vec<f32>,
    pub ux: Vec<f32>,
    pub uy: Vec<f32>,
    pub uz: Vec<f32>,
}

/// Neighbour coordinate one step in the sign of `delta` (periodic).
pub fn periodic_shift(coordinate: usize, delta: i32, extent: usize) -> usize {
    match delta.signum() {
        1 => {
            if coordinate + 1 == extent {
                0
            } else {
                coordinate + 1
            }
        }
        -1 => {
            if coordinate == 0 {
                extent - 1
            } else {
                coordinate - 1
            }
        }
        _ => coordinate,
    }
}

/// Second-order f32 equilibrium for any descriptor.
pub fn lattice_equilibrium<D: LatticeDescriptor>(
    q: usize,
    rho: f32,
    ux: f32,
    uy: f32,
    uz: f32,
) -> f32 {
    let cu = 3.0 * (D::cx(q) as f32 * ux + D::cy(q) as f32 * uy + D::cz(q) as f32 * uz);
    let uu = 1.5 * (ux * ux + uy * uy + uz * uz);
    D::weight(q) * rho * (1.0 + cu + 0.5 * cu * cu - uu)
}

/// Guo forcing term with the `(1 - omega/2)` prefactor; zero when the acceleration is zero.
#[allow(clippy::too_many_arguments)]
pub fn guo_force<D: LatticeDescriptor>(
    q: usize,
    rho: f32,
    ux: f32,
    uy: f32,
    uz: f32,
    ax: f32,
    ay: f32,
    az: f32,
    omega: f32,
) -> f32 {
    if ax == 0.0 && ay == 0.0 && az == 0.0 {
        return 0.0;
    }
    let (cx, cy, cz) = (D::cx(q) as f32, D::cy(q) as f32, D::cz(q) as f32);
    let (fx, fy, fz) = (rho * ax, rho * ay, rho * az);
    let cu = cx * ux + cy * uy + cz * uz;
    let cf = cx * fx + cy * fy + cz * fz;
    let uf = ux * fx + uy * fy + uz * fz;
    D::weight(q) * (1.0 - 0.5 * omega) * (3.0 * (cf - uf) + 9.0 * cu * cf)
}

/// Smagorinsky-Lilly effective relaxation time from the non-equilibrium momentum flux.
pub fn smagorinsky_relaxation_time<D: LatticeDescriptor>(
    populations: &[f32],
    rho: f32,
    ux: f32,
    uy: f32,
    uz: f32,
    base_tau: f32,
) -> f32 {
    if rho.is_nan() || rho <= 0.0 {
        return base_tau;
    }
    let mut pi = [[0.0f32; 3]; 3];
    for (d, f) in populations.iter().enumerate().take(D::Q) {
        let neq = f - lattice_equilibrium::<D>(d, rho, ux, uy, uz);
        let e = [D::cx(d) as f32, D::cy(d) as f32, D::cz(d) as f32];
        for a in 0..3 {
            for b in 0..3 {
                pi[a][b] += e[a] * e[b] * neq;
            }
        }
    }
    let mut q = 0.0f32;
    for row in &pi {
        for v in row {
            q += v * v;
        }
    }
    let pi_f = std::f32::consts::PI;
    let coefficient = 16.0f32 * std::f32::consts::SQRT_2 / (3.0f32 * pi_f * pi_f);
    0.5 * (base_tau + (base_tau * base_tau + coefficient * q.sqrt() / rho).sqrt())
}

#[derive(Clone, Copy, Debug)]
struct Lattice3 {
    nx: usize,
    ny: usize,
    nz: usize,
}

impl Lattice3 {
    fn normalized<D: LatticeDescriptor>(
        nx: usize,
        ny: usize,
        nz: usize,
        what: &str,
    ) -> Result<Self> {
        require(
            nx >= 2 && ny >= 2,
            format!("{what} grid requires nx, ny >= 2"),
        )?;
        let nz = if D::DIMENSIONS == 2 { 1 } else { nz };
        require(
            D::DIMENSIONS == 2 || nz >= 2,
            format!("3-D {what} grid requires nz >= 2"),
        )?;
        nx.checked_mul(ny)
            .and_then(|v| v.checked_mul(nz))
            .and_then(|v| v.checked_mul(D::Q))
            .ok_or_else(|| SolverError::InvalidInput("LBM grid too large".into()))?;
        Ok(Self { nx, ny, nz })
    }
    fn cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    fn coords(&self, n: usize) -> (usize, usize, usize) {
        let x = n % self.nx;
        let yz = n / self.nx;
        (x, yz % self.ny, yz / self.ny)
    }
    fn index(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.ny + y) * self.nx + x
    }
    /// Cell reached by moving one step along `direction` (or against it when `sign == -1`).
    fn step<D: LatticeDescriptor>(
        &self,
        x: usize,
        y: usize,
        z: usize,
        direction: usize,
        sign: i32,
    ) -> usize {
        self.index(
            periodic_shift(x, sign * D::cx(direction), self.nx),
            periodic_shift(y, sign * D::cy(direction), self.ny),
            periodic_shift(z, sign * D::cz(direction), self.nz),
        )
    }
    /// Taylor-Green state `(1, A sin cos zf, -A cos sin zf, 0)` in f32, with `zf = cos(2 pi z)` in 3-D.
    fn taylor_green<D: LatticeDescriptor>(
        &self,
        amplitude: f32,
        x: usize,
        y: usize,
        z: usize,
    ) -> [f32; 4] {
        let two_pi = 2.0 * std::f32::consts::PI;
        let xf = (x as f32 + 0.5) / self.nx as f32;
        let yf = (y as f32 + 0.5) / self.ny as f32;
        let (sx, cx, sy, cy) = (
            (two_pi * xf).sin(),
            (two_pi * xf).cos(),
            (two_pi * yf).sin(),
            (two_pi * yf).cos(),
        );
        let zf = if D::DIMENSIONS == 3 {
            (two_pi * ((z as f32 + 0.5) / self.nz as f32)).cos()
        } else {
            1.0
        };
        [
            1.0,
            amplitude * sx * cy * zf,
            -amplitude * cx * sy * zf,
            0.0,
        ]
    }
}

fn accumulate<D: LatticeDescriptor>(fin: &[f32]) -> (f32, f32, f32, f32) {
    let (mut rho, mut ux, mut uy, mut uz) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    for (d, v) in fin.iter().enumerate().take(D::Q) {
        rho += v;
        ux += v * D::cx(d) as f32;
        uy += v * D::cy(d) as f32;
        uz += v * D::cz(d) as f32;
    }
    if rho > 0.0 {
        (rho, ux / rho, uy / rho, uz / rho)
    } else {
        (rho, 0.0, 0.0, 0.0)
    }
}

/// Thread-safe single-population-grid Esoteric-Pull BGK (Lehmann, Computation 2022, 10(6), 92)
/// with Guo forcing, optional per-cell acceleration, Smagorinsky LES and mid-link bounce-back.
#[derive(Clone, Debug)]
pub struct EsotericPullSolver<D: LatticeDescriptor> {
    config: InPlaceLbmConfig,
    grid: Lattice3,
    cells: usize,
    f: Vec<f32>,
    local_a: Option<[Vec<f32>; 3]>,
    solid: Vec<u8>,
    solid_current: Vec<f32>,
    solid_next: Vec<f32>,
    time_step: u64,
    _d: PhantomData<D>,
}

pub type D2Q9InPlaceSolver = EsotericPullSolver<super::D2Q9InPlaceDescriptor>;
pub type D3Q19Solver = EsotericPullSolver<super::D3Q19Descriptor>;
pub type D3Q27Solver = EsotericPullSolver<super::D3Q27Descriptor>;

impl<D: LatticeDescriptor> EsotericPullSolver<D> {
    pub fn new(config: InPlaceLbmConfig) -> Result<Self> {
        let grid = Lattice3::normalized::<D>(config.nx, config.ny, config.nz, "LBM")?;
        require(
            config.tau > 0.5,
            "LBM tau must be > 0.5 for positive viscosity",
        )?;
        let config = InPlaceLbmConfig {
            nz: grid.nz,
            ..config
        };
        let cells = grid.cells();
        let mut s = Self {
            config,
            grid,
            cells,
            f: vec![0.0; D::Q * cells],
            local_a: None,
            solid: Vec::new(),
            solid_current: Vec::new(),
            solid_next: Vec::new(),
            time_step: 0,
            _d: PhantomData,
        };
        s.initialize_uniform(1.0, 0.0, 0.0, 0.0);
        Ok(s)
    }
    fn fi(&self, d: usize, n: usize) -> f32 {
        self.f[d * self.cells + n]
    }
    fn set(&mut self, d: usize, n: usize, v: f32) {
        self.f[d * self.cells + n] = v;
    }
    fn is_solid(&self, n: usize) -> bool {
        !self.solid.is_empty() && self.solid[n] != 0
    }
    fn initialize(&mut self, state: impl Fn(usize, usize, usize) -> [f32; 4]) {
        self.time_step = 0;
        for n in 0..self.cells {
            let (x, y, z) = self.grid.coords(n);
            let s = state(x, y, z);
            let solid = self.is_solid(n);
            let (ux, uy, uz) = if solid {
                (0.0, 0.0, 0.0)
            } else {
                (s[1], s[2], s[3])
            };
            self.set(0, n, lattice_equilibrium::<D>(0, s[0], ux, uy, uz));
            for i in (1..D::Q).step_by(2) {
                let neighbor = self.grid.step::<D>(x, y, z, i + 1, 1);
                self.set(i, neighbor, lattice_equilibrium::<D>(i, s[0], ux, uy, uz));
                self.set(i + 1, n, lattice_equilibrium::<D>(i + 1, s[0], ux, uy, uz));
            }
        }
    }
    pub fn initialize_uniform(&mut self, rho: f32, ux: f32, uy: f32, uz: f32) {
        self.initialize(|_, _, _| [rho, ux, uy, uz]);
    }
    pub fn initialize_taylor_green(&mut self, amplitude: f32) {
        let g = self.grid;
        self.initialize(|x, y, z| g.taylor_green::<D>(amplitude, x, y, z));
    }
    pub fn step(&mut self, count: usize) {
        for _ in 0..count {
            self.step_once();
        }
    }
    pub fn set_local_body_acceleration(
        &mut self,
        ax: &[f32],
        ay: &[f32],
        az: &[f32],
    ) -> Result<()> {
        require(
            ax.len() == self.cells && ay.len() == self.cells && az.len() == self.cells,
            "local LBM acceleration field size mismatch",
        )?;
        self.local_a = Some([ax.to_vec(), ay.to_vec(), az.to_vec()]);
        Ok(())
    }
    pub fn clear_local_body_acceleration(&mut self) {
        self.local_a = None;
    }
    pub fn has_local_body_acceleration(&self) -> bool {
        self.local_a.is_some()
    }
    /// Stationary mid-link bounce-back on a copied voxel mask (non-zero = solid).
    pub fn set_solid_mask(&mut self, mask: &[u8]) -> Result<()> {
        require(mask.len() == self.cells, "LBM solid mask size mismatch")?;
        self.solid = mask.to_vec();
        self.solid_current = vec![0.0; D::Q * self.cells];
        self.solid_next = vec![0.0; D::Q * self.cells];
        Ok(())
    }
    pub fn clear_solid_mask(&mut self) {
        self.solid.clear();
        self.solid_current.clear();
        self.solid_next.clear();
    }
    pub fn has_solid_mask(&self) -> bool {
        !self.solid.is_empty()
    }
    pub fn solid_mask(&self) -> &[u8] {
        &self.solid
    }
    fn acceleration(&self, n: usize) -> (f32, f32, f32) {
        let c = &self.config;
        match &self.local_a {
            Some([ax, ay, az]) => (
                c.acceleration_x + ax[n],
                c.acceleration_y + ay[n],
                c.acceleration_z + az[n],
            ),
            None => (
                c.acceleration_x + 0.0,
                c.acceleration_y + 0.0,
                c.acceleration_z + 0.0,
            ),
        }
    }
    /// Logical (post-streaming) populations of a cell and the backward neighbours of each pair.
    fn load_cell(&self, n: usize, fin: &mut [f32; MAX_Q], neighbors: &mut [usize; MAX_Q]) {
        let (x, y, z) = self.grid.coords(n);
        let odd = self.time_step & 1 == 1;
        fin[0] = self.fi(0, n);
        for (pair, i) in (1..D::Q).step_by(2).enumerate() {
            let neighbor = self.grid.step::<D>(x, y, z, i + 1, 1);
            neighbors[pair] = neighbor;
            if odd {
                fin[i] = self.fi(i + 1, neighbor);
                fin[i + 1] = self.fi(i, n);
            } else {
                fin[i] = self.fi(i, neighbor);
                fin[i + 1] = self.fi(i + 1, n);
            }
        }
    }
    fn load(&self, n: usize) -> [f32; MAX_Q] {
        let mut fin = [0.0; MAX_Q];
        let mut nb = [0; MAX_Q];
        self.load_cell(n, &mut fin, &mut nb);
        fin
    }
    fn store_cell(&mut self, n: usize, neighbors: &[usize; MAX_Q], fout: &[f32; MAX_Q]) {
        let odd = self.time_step & 1 == 1;
        self.set(0, n, fout[0]);
        for (pair, i) in (1..D::Q).step_by(2).enumerate() {
            let neighbor = neighbors[pair];
            if odd {
                self.set(i, n, fout[i]);
                self.set(i + 1, neighbor, fout[i + 1]);
            } else {
                self.set(i + 1, n, fout[i]);
                self.set(i, neighbor, fout[i + 1]);
            }
        }
    }
    fn collide_cell(&self, p: &mut [f32; MAX_Q], n: usize) {
        let (rho, mut ux, mut uy, mut uz) = accumulate::<D>(&p[..D::Q]);
        let (ax, ay, az) = self.acceleration(n);
        ux += 0.5 * ax;
        uy += 0.5 * ay;
        uz += 0.5 * az;
        let tau = if self.config.smagorinsky_les {
            smagorinsky_relaxation_time::<D>(&p[..D::Q], rho, ux, uy, uz, self.config.tau)
        } else {
            self.config.tau
        };
        let omega = 1.0 / tau;
        for (d, v) in p.iter_mut().enumerate().take(D::Q) {
            *v -= omega * (*v - lattice_equilibrium::<D>(d, rho, ux, uy, uz));
            *v += guo_force::<D>(d, rho, ux, uy, uz, ax, ay, az, omega);
        }
    }
    fn step_once(&mut self) {
        if !self.solid.is_empty() {
            self.step_once_with_solids();
            return;
        }
        // Each cell reads and writes only its own pair slots, so the in-place update is order-free.
        let mut fout = [0.0; MAX_Q];
        let mut neighbors = [0; MAX_Q];
        for n in 0..self.cells {
            self.load_cell(n, &mut fout, &mut neighbors);
            self.collide_cell(&mut fout, n);
            self.store_cell(n, &neighbors, &fout);
        }
        self.time_step += 1;
    }
    fn step_once_with_solids(&mut self) {
        let cells = self.cells;
        for n in 0..cells {
            let fin = self.load(n);
            for d in 0..D::Q {
                self.solid_current[d * cells + n] = fin[d];
            }
        }
        for n in 0..cells {
            if self.solid[n] != 0 {
                for d in 0..D::Q {
                    self.solid_next[d * cells + n] = self.solid_current[d * cells + n];
                }
                continue;
            }
            let mut fout = [0.0; MAX_Q];
            for d in 0..D::Q {
                fout[d] = self.solid_current[d * cells + n];
            }
            self.collide_cell(&mut fout, n);
            let (x, y, z) = self.grid.coords(n);
            for d in 0..D::Q {
                let target = self.grid.step::<D>(x, y, z, d, 1);
                if self.solid[target] != 0 {
                    self.solid_next[D::opposite(d) * cells + n] = fout[d];
                } else {
                    self.solid_next[d * cells + target] = fout[d];
                }
            }
        }
        let new_odd = (self.time_step + 1) & 1 == 1;
        for n in 0..cells {
            let (x, y, z) = self.grid.coords(n);
            let v0 = self.solid_next[n];
            self.set(0, n, v0);
            for i in (1..D::Q).step_by(2) {
                let backward = self.grid.step::<D>(x, y, z, i + 1, 1);
                if new_odd {
                    let forward = self.grid.step::<D>(x, y, z, i, 1);
                    let a = self.solid_next[(i + 1) * cells + forward];
                    let b = self.solid_next[i * cells + n];
                    self.set(i, n, a);
                    self.set(i + 1, backward, b);
                } else {
                    let a = self.solid_next[i * cells + n];
                    let b = self.solid_next[(i + 1) * cells + n];
                    self.set(i, backward, a);
                    self.set(i + 1, n, b);
                }
            }
        }
        self.time_step += 1;
    }
    /// Half-force-corrected velocity of a cell (solid cells report zero velocity when `zero_solid`).
    fn cell_state(&self, n: usize, zero_solid: bool) -> (f32, f32, f32, f32) {
        let fin = self.load(n);
        let (rho, mut ux, mut uy, mut uz) = accumulate::<D>(&fin[..D::Q]);
        if zero_solid && self.is_solid(n) {
            return (rho, 0.0, 0.0, 0.0);
        }
        let c = &self.config;
        let (lx, ly, lz) = match &self.local_a {
            Some([ax, ay, az]) => (ax[n], ay[n], az[n]),
            None => (0.0, 0.0, 0.0),
        };
        ux += 0.5 * (c.acceleration_x + lx);
        uy += 0.5 * (c.acceleration_y + ly);
        uz += 0.5 * (c.acceleration_z + lz);
        (rho, ux, uy, uz)
    }
    pub fn compute_macroscopic(&self) -> MacroscopicFields {
        let mut m = MacroscopicFields {
            rho: vec![0.0; self.cells],
            ux: vec![0.0; self.cells],
            uy: vec![0.0; self.cells],
            uz: vec![0.0; self.cells],
        };
        for n in 0..self.cells {
            let (rho, ux, uy, uz) = self.cell_state(n, true);
            m.rho[n] = rho;
            m.ux[n] = ux;
            m.uy[n] = uy;
            m.uz[n] = uz;
        }
        m
    }
    pub fn mass(&self) -> f64 {
        (0..self.cells)
            .map(|n| f64::from(self.load(n)[..D::Q].iter().sum::<f32>()))
            .sum()
    }
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * (0..self.cells)
            .map(|n| {
                let (rho, ux, uy, uz) = self.cell_state(n, false);
                let (ux, uy, uz) = (f64::from(ux), f64::from(uy), f64::from(uz));
                f64::from(rho) * (ux * ux + uy * uy + uz * uz)
            })
            .sum::<f64>()
    }
    pub fn max_speed(&self) -> f32 {
        (0..self.cells)
            .map(|n| {
                let (_, ux, uy, uz) = self.cell_state(n, false);
                (ux * ux + uy * uy + uz * uz).sqrt()
            })
            .fold(f32::NEG_INFINITY, f32::max)
    }
    /// Populations in logical (post-streaming) order, `[direction * cells + cell]`.
    pub fn logical_populations(&self) -> Vec<f32> {
        let mut r = vec![0.0; D::Q * self.cells];
        for n in 0..self.cells {
            let fin = self.load(n);
            for d in 0..D::Q {
                r[d * self.cells + n] = fin[d];
            }
        }
        r
    }
    pub fn config(&self) -> &InPlaceLbmConfig {
        &self.config
    }
    pub fn cells(&self) -> usize {
        self.cells
    }
    pub fn time_step(&self) -> u64 {
        self.time_step
    }
    pub fn population_bytes(&self) -> usize {
        D::Q * self.cells * std::mem::size_of::<f32>()
    }
    pub fn two_lattice_population_bytes(&self) -> usize {
        2 * self.population_bytes()
    }
    /// Physical (in-place, time-step-dependent) storage.
    pub fn raw_storage(&self) -> &[f32] {
        &self.f
    }
}

/// Conventional two-lattice pull-stream/collide with Guo forcing and an optional two-relaxation-time
/// (TRT) odd rate; the correctness oracle for the in-place solver.
#[derive(Clone, Debug)]
pub struct OneStepPullSolver<D: LatticeDescriptor> {
    config: InPlaceLbmConfig,
    odd_tau: f32,
    grid: Lattice3,
    cells: usize,
    f: Vec<f32>,
    next: Vec<f32>,
    _d: PhantomData<D>,
}

pub type D2Q9PullSolver = OneStepPullSolver<super::D2Q9InPlaceDescriptor>;
pub type D3Q19PullSolver = OneStepPullSolver<super::D3Q19Descriptor>;
pub type D3Q27PullSolver = OneStepPullSolver<super::D3Q27Descriptor>;

impl<D: LatticeDescriptor> OneStepPullSolver<D> {
    /// `odd_tau == 0` means plain BGK with `config.tau`.
    pub fn new(config: InPlaceLbmConfig, odd_tau: f32) -> Result<Self> {
        let grid = Lattice3::normalized::<D>(config.nx, config.ny, config.nz, "LBM")?;
        let config = InPlaceLbmConfig {
            nz: grid.nz,
            ..config
        };
        let odd_tau = if odd_tau == 0.0 { config.tau } else { odd_tau };
        require(
            config.tau > 0.5 && config.tau.is_finite() && odd_tau > 0.5 && odd_tau.is_finite(),
            "LBM tau must be > 0.5 for positive viscosity",
        )?;
        let cells = grid.cells();
        let mut s = Self {
            config,
            odd_tau,
            grid,
            cells,
            f: vec![0.0; D::Q * cells],
            next: vec![0.0; D::Q * cells],
            _d: PhantomData,
        };
        s.initialize_uniform(1.0, 0.0, 0.0, 0.0);
        Ok(s)
    }
    fn initialize(&mut self, state: impl Fn(usize, usize, usize) -> [f32; 4]) {
        for n in 0..self.cells {
            let (x, y, z) = self.grid.coords(n);
            let s = state(x, y, z);
            for d in 0..D::Q {
                self.f[d * self.cells + n] = lattice_equilibrium::<D>(d, s[0], s[1], s[2], s[3]);
            }
        }
    }
    pub fn initialize_uniform(&mut self, rho: f32, ux: f32, uy: f32, uz: f32) {
        self.initialize(|_, _, _| [rho, ux, uy, uz]);
    }
    pub fn initialize_taylor_green(&mut self, amplitude: f32) {
        let g = self.grid;
        self.initialize(|x, y, z| g.taylor_green::<D>(amplitude, x, y, z));
    }
    pub fn step(&mut self, count: usize) {
        for _ in 0..count {
            self.step_once();
        }
    }
    fn step_once(&mut self) {
        let c = self.config;
        let cells = self.cells;
        let omega = 1.0f32 / c.tau;
        for n in 0..cells {
            let (x, y, z) = self.grid.coords(n);
            let mut fin = [0.0f32; MAX_Q];
            for (d, v) in fin.iter_mut().enumerate().take(D::Q) {
                *v = self.f[d * cells + self.grid.step::<D>(x, y, z, d, -1)];
            }
            let (rho, mut ux, mut uy, mut uz) = accumulate::<D>(&fin[..D::Q]);
            ux += 0.5 * c.acceleration_x;
            uy += 0.5 * c.acceleration_y;
            uz += 0.5 * c.acceleration_z;
            if self.odd_tau != c.tau {
                let omega_odd = 1.0f32 / self.odd_tau;
                let mut feq = [0.0f32; MAX_Q];
                let mut force = [0.0f32; MAX_Q];
                for d in 0..D::Q {
                    feq[d] = lattice_equilibrium::<D>(d, rho, ux, uy, uz);
                    force[d] = guo_force::<D>(
                        d,
                        rho,
                        ux,
                        uy,
                        uz,
                        c.acceleration_x,
                        c.acceleration_y,
                        c.acceleration_z,
                        0.0,
                    );
                }
                for d in 0..D::Q {
                    let o = D::opposite(d);
                    let even = 0.5 * ((fin[d] - feq[d]) + (fin[o] - feq[o]));
                    let odd = 0.5 * ((fin[d] - feq[d]) - (fin[o] - feq[o]));
                    let forcing = 0.5 * (1.0 - 0.5 * omega) * (force[d] + force[o])
                        + 0.5 * (1.0 - 0.5 * omega_odd) * (force[d] - force[o]);
                    self.next[d * cells + n] = fin[d] - omega * even - omega_odd * odd + forcing;
                }
                continue;
            }
            for d in 0..D::Q {
                let feq = lattice_equilibrium::<D>(d, rho, ux, uy, uz);
                self.next[d * cells + n] = fin[d] - omega * (fin[d] - feq)
                    + guo_force::<D>(
                        d,
                        rho,
                        ux,
                        uy,
                        uz,
                        c.acceleration_x,
                        c.acceleration_y,
                        c.acceleration_z,
                        omega,
                    );
            }
        }
        std::mem::swap(&mut self.f, &mut self.next);
    }
    pub fn compute_macroscopic(&self) -> MacroscopicFields {
        let c = self.config;
        let mut m = MacroscopicFields {
            rho: vec![0.0; self.cells],
            ux: vec![0.0; self.cells],
            uy: vec![0.0; self.cells],
            uz: vec![0.0; self.cells],
        };
        for n in 0..self.cells {
            let fin: Vec<f32> = (0..D::Q).map(|d| self.f[d * self.cells + n]).collect();
            let (rho, ux, uy, uz) = accumulate::<D>(&fin);
            m.rho[n] = rho;
            m.ux[n] = ux + 0.5 * c.acceleration_x;
            m.uy[n] = uy + 0.5 * c.acceleration_y;
            m.uz[n] = uz + 0.5 * c.acceleration_z;
        }
        m
    }
    pub fn mass(&self) -> f64 {
        (0..self.cells)
            .map(|n| f64::from((0..D::Q).map(|d| self.f[d * self.cells + n]).sum::<f32>()))
            .sum()
    }
    /// Populations `[direction * cells + cell]`.
    pub fn populations(&self) -> &[f32] {
        &self.f
    }
    pub fn population_bytes(&self) -> usize {
        2 * D::Q * self.cells * std::mem::size_of::<f32>()
    }
    pub fn config(&self) -> &InPlaceLbmConfig {
        &self.config
    }
    pub fn cells(&self) -> usize {
        self.cells
    }
    pub fn odd_relaxation_time(&self) -> f32 {
        self.odd_tau
    }
}

/// A 16-bit population storage format. `direction` lets codecs store deviations from the weight.
pub trait PopulationCodec {
    fn encode(value: f32, direction: usize) -> u16;
    fn decode(bits: u16, direction: usize) -> f32;
}

/// IEEE 754 binary16 with round-half-away mantissa rounding, as upstream.
#[derive(Clone, Copy, Debug, Default)]
pub struct IeeeHalfCodec;

impl IeeeHalfCodec {
    pub fn encode_value(value: f32) -> u16 {
        if value.is_nan() {
            return 0x7e00;
        }
        let sign: u16 = if value.is_sign_negative() { 0x8000 } else { 0 };
        let a = value.abs();
        if a.is_infinite() {
            return sign | 0x7c00;
        }
        if a == 0.0 {
            return sign;
        }
        if a < 2.0f32.powi(-14) {
            let mantissa = (a * 2.0f32.powi(24)).round() as i64;
            if mantissa <= 0 {
                return sign;
            }
            return sign | mantissa.min(1023) as u16;
        }
        let mut exponent = a.log2().floor() as i32;
        if exponent > 15 {
            return sign | 0x7c00;
        }
        let normalized = a * 2.0f32.powi(-exponent);
        let mut mantissa = ((normalized - 1.0) * 1024.0).round() as i64;
        if mantissa == 1024 {
            mantissa = 0;
            exponent += 1;
            if exponent > 15 {
                return sign | 0x7c00;
            }
        }
        sign | (((exponent + 15) as u16) << 10) | mantissa as u16
    }
    pub fn decode_value(bits: u16) -> f32 {
        let sign = if bits & 0x8000 != 0 { -1.0f32 } else { 1.0 };
        let exponent = i32::from((bits >> 10) & 0x1f);
        let mantissa = bits & 0x03ff;
        if exponent == 0x1f {
            return if mantissa != 0 {
                f32::NAN
            } else {
                sign * f32::INFINITY
            };
        }
        if exponent == 0 {
            return sign * (f32::from(mantissa) * 2.0f32.powi(-24));
        }
        sign * ((1.0 + f32::from(mantissa) / 1024.0) * 2.0f32.powi(exponent - 15))
    }
}

impl PopulationCodec for IeeeHalfCodec {
    fn encode(value: f32, _: usize) -> u16 {
        Self::encode_value(value)
    }
    fn decode(bits: u16, _: usize) -> f32 {
        Self::decode_value(bits)
    }
}

/// bfloat16 with round-to-nearest-even on the dropped 16 bits.
#[derive(Clone, Copy, Debug, Default)]
pub struct BFloat16Codec;

impl PopulationCodec for BFloat16Codec {
    fn encode(value: f32, _: usize) -> u16 {
        let bits = value.to_bits();
        let lsb = (bits >> 16) & 1;
        (bits.wrapping_add(0x7fff + lsb) >> 16) as u16
    }
    fn decode(bits: u16, _: usize) -> f32 {
        f32::from_bits(u32::from(bits) << 16)
    }
}

/// binary16 of the deviation from the lattice weight.
#[derive(Clone, Copy, Debug, Default)]
pub struct ShiftedHalfCodec<D>(PhantomData<D>);

impl<D: LatticeDescriptor> PopulationCodec for ShiftedHalfCodec<D> {
    fn encode(value: f32, direction: usize) -> u16 {
        IeeeHalfCodec::encode_value(value - D::weight(direction))
    }
    fn decode(bits: u16, direction: usize) -> f32 {
        IeeeHalfCodec::decode_value(bits) + D::weight(direction)
    }
}

/// Fixed-point int16 deviation from the weight with scale 2^15.
#[derive(Clone, Copy, Debug, Default)]
pub struct PackedDeviation16Codec<D>(PhantomData<D>);

impl<D: LatticeDescriptor> PopulationCodec for PackedDeviation16Codec<D> {
    fn encode(value: f32, direction: usize) -> u16 {
        let deviation = value - D::weight(direction);
        let quantized = (deviation * 32768.0).clamp(-32768.0, 32767.0).round() as i64;
        quantized as i16 as u16
    }
    fn decode(bits: u16, direction: usize) -> f32 {
        D::weight(direction) + f32::from(bits as i16) / 32768.0
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompressedLbmConfig {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub tau: f32,
}

impl Default for CompressedLbmConfig {
    fn default() -> Self {
        Self {
            nx: 64,
            ny: 64,
            nz: 64,
            tau: 0.6,
        }
    }
}

/// Two-lattice BGK pull with 16-bit populations (half the memory traffic of f32).
#[derive(Clone, Debug)]
pub struct CompressedPopulationPullSolver<D: LatticeDescriptor, C: PopulationCodec> {
    config: CompressedLbmConfig,
    grid: Lattice3,
    cells: usize,
    f: Vec<u16>,
    next: Vec<u16>,
    _p: PhantomData<(D, C)>,
}

impl<D: LatticeDescriptor, C: PopulationCodec> CompressedPopulationPullSolver<D, C> {
    pub fn new(config: CompressedLbmConfig) -> Result<Self> {
        let grid = Lattice3::normalized::<D>(config.nx, config.ny, config.nz, "compressed LBM")?;
        require(config.tau > 0.5, "compressed LBM tau must be > 0.5")?;
        let config = CompressedLbmConfig {
            nz: grid.nz,
            ..config
        };
        let cells = grid.cells();
        let mut s = Self {
            config,
            grid,
            cells,
            f: vec![0; D::Q * cells],
            next: vec![0; D::Q * cells],
            _p: PhantomData,
        };
        s.initialize_uniform(1.0, 0.0, 0.0, 0.0);
        Ok(s)
    }
    fn load(&self, d: usize, n: usize) -> f32 {
        C::decode(self.f[d * self.cells + n], d)
    }
    fn initialize(&mut self, state: impl Fn(usize, usize, usize) -> [f32; 4]) {
        for n in 0..self.cells {
            let (x, y, z) = self.grid.coords(n);
            let s = state(x, y, z);
            for d in 0..D::Q {
                let v = C::encode(lattice_equilibrium::<D>(d, s[0], s[1], s[2], s[3]), d);
                self.f[d * self.cells + n] = v;
                self.next[d * self.cells + n] = v;
            }
        }
    }
    pub fn initialize_uniform(&mut self, rho: f32, ux: f32, uy: f32, uz: f32) {
        self.initialize(|_, _, _| [rho, ux, uy, uz]);
    }
    /// As upstream, the 3-D factor is `cos(2 pi (z + 0.5) / nz)` evaluated left to right.
    pub fn initialize_taylor_green(&mut self, amplitude: f32) {
        let g = self.grid;
        self.initialize(|x, y, z| {
            let mut s = g.taylor_green::<super::D2Q9InPlaceDescriptor>(amplitude, x, y, z);
            if D::DIMENSIONS == 3 {
                let two_pi = 2.0 * std::f32::consts::PI;
                let zf = (two_pi * (z as f32 + 0.5) / g.nz as f32).cos();
                s[1] *= zf;
                s[2] *= zf;
            }
            s
        });
    }
    pub fn step(&mut self, count: usize) {
        for _ in 0..count {
            self.step_once();
        }
    }
    fn step_once(&mut self) {
        let omega = 1.0f32 / self.config.tau;
        let cells = self.cells;
        for n in 0..cells {
            let (x, y, z) = self.grid.coords(n);
            let mut fin = [0.0f32; MAX_Q];
            for (d, v) in fin.iter_mut().enumerate().take(D::Q) {
                *v = self.load(d, self.grid.step::<D>(x, y, z, d, -1));
            }
            let (rho, ux, uy, uz) = accumulate::<D>(&fin[..D::Q]);
            for d in 0..D::Q {
                self.next[d * cells + n] = C::encode(
                    fin[d] - omega * (fin[d] - lattice_equilibrium::<D>(d, rho, ux, uy, uz)),
                    d,
                );
            }
        }
        std::mem::swap(&mut self.f, &mut self.next);
    }
    pub fn compute_macroscopic(&self) -> MacroscopicFields {
        let mut m = MacroscopicFields {
            rho: vec![0.0; self.cells],
            ux: vec![0.0; self.cells],
            uy: vec![0.0; self.cells],
            uz: vec![0.0; self.cells],
        };
        for n in 0..self.cells {
            let fin: Vec<f32> = (0..D::Q).map(|d| self.load(d, n)).collect();
            let (rho, ux, uy, uz) = accumulate::<D>(&fin);
            m.rho[n] = rho;
            m.ux[n] = ux;
            m.uy[n] = uy;
            m.uz[n] = uz;
        }
        m
    }
    pub fn mass(&self) -> f64 {
        let mut sum = 0.0f64;
        for n in 0..self.cells {
            for d in 0..D::Q {
                sum += f64::from(self.load(d, n));
            }
        }
        sum
    }
    pub fn population_bytes(&self) -> usize {
        D::Q * self.cells * std::mem::size_of::<u16>()
    }
    /// Encoded populations `[direction * cells + cell]`.
    pub fn raw_storage(&self) -> &[u16] {
        &self.f
    }
}

/// Floating-point type for [`PrecisionPullSolver`].
pub trait LbmReal:
    Copy
    + std::fmt::Debug
    + PartialOrd
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<Output = Self>
    + std::ops::Div<Output = Self>
    + std::ops::AddAssign
{
    fn from_f64(v: f64) -> Self;
    fn from_f32(v: f32) -> Self;
    fn to_f64(self) -> f64;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    const PI: Self;
}

impl LbmReal for f32 {
    fn from_f64(v: f64) -> Self {
        v as f32
    }
    fn from_f32(v: f32) -> Self {
        v
    }
    fn to_f64(self) -> f64 {
        f64::from(self)
    }
    fn sin(self) -> Self {
        f32::sin(self)
    }
    fn cos(self) -> Self {
        f32::cos(self)
    }
    const PI: Self = std::f32::consts::PI;
}

impl LbmReal for f64 {
    fn from_f64(v: f64) -> Self {
        v
    }
    fn from_f32(v: f32) -> Self {
        f64::from(v)
    }
    fn to_f64(self) -> f64 {
        self
    }
    fn sin(self) -> Self {
        f64::sin(self)
    }
    fn cos(self) -> Self {
        f64::cos(self)
    }
    const PI: Self = std::f64::consts::PI;
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TypedLbmConfig<R> {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub tau: R,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TypedMacroscopicFields<R> {
    pub rho: Vec<R>,
    pub ux: Vec<R>,
    pub uy: Vec<R>,
    pub uz: Vec<R>,
}

fn typed_equilibrium<D: LatticeDescriptor, R: LbmReal>(q: usize, rho: R, ux: R, uy: R, uz: R) -> R {
    let (cx, cy, cz) = (
        R::from_f64(f64::from(D::cx(q))),
        R::from_f64(f64::from(D::cy(q))),
        R::from_f64(f64::from(D::cz(q))),
    );
    let cu = R::from_f64(3.0) * (cx * ux + cy * uy + cz * uz);
    let uu = R::from_f64(1.5) * (ux * ux + uy * uy + uz * uz);
    R::from_f32(D::weight(q)) * rho * (R::from_f64(1.0) + cu + R::from_f64(0.5) * cu * cu - uu)
}

/// Precision-generic (f32/f64) two-lattice BGK pull; upstream's `long double` variant has no Rust
/// equivalent, so mass is accumulated in f64.
#[derive(Clone, Debug)]
pub struct PrecisionPullSolver<D: LatticeDescriptor, R: LbmReal> {
    config: TypedLbmConfig<R>,
    grid: Lattice3,
    cells: usize,
    f: Vec<R>,
    next: Vec<R>,
    _d: PhantomData<D>,
}

impl<D: LatticeDescriptor, R: LbmReal> PrecisionPullSolver<D, R> {
    pub fn new(config: TypedLbmConfig<R>) -> Result<Self> {
        let grid = Lattice3::normalized::<D>(config.nx, config.ny, config.nz, "LBM")?;
        require(
            config.tau > R::from_f64(0.5),
            "LBM tau must be > 0.5 for positive viscosity",
        )?;
        let config = TypedLbmConfig {
            nz: grid.nz,
            ..config
        };
        let cells = grid.cells();
        let zero = R::from_f64(0.0);
        let mut s = Self {
            config,
            grid,
            cells,
            f: vec![zero; D::Q * cells],
            next: vec![zero; D::Q * cells],
            _d: PhantomData,
        };
        s.initialize_uniform(R::from_f64(1.0), zero, zero, zero);
        Ok(s)
    }
    fn initialize(&mut self, state: impl Fn(usize, usize, usize) -> [R; 4]) {
        for n in 0..self.cells {
            let (x, y, z) = self.grid.coords(n);
            let s = state(x, y, z);
            for d in 0..D::Q {
                let v = typed_equilibrium::<D, R>(d, s[0], s[1], s[2], s[3]);
                self.f[d * self.cells + n] = v;
                self.next[d * self.cells + n] = v;
            }
        }
    }
    pub fn initialize_uniform(&mut self, rho: R, ux: R, uy: R, uz: R) {
        self.initialize(|_, _, _| [rho, ux, uy, uz]);
    }
    pub fn initialize_taylor_green(&mut self, amplitude: R) {
        let g = self.grid;
        let two_pi = R::from_f64(2.0) * R::PI;
        let half = R::from_f64(0.5);
        let one = R::from_f64(1.0);
        let zero = R::from_f64(0.0);
        self.initialize(|x, y, z| {
            let xf = (R::from_f64(x as f64) + half) / R::from_f64(g.nx as f64);
            let yf = (R::from_f64(y as f64) + half) / R::from_f64(g.ny as f64);
            let (sx, cx, sy, cy) = (
                (two_pi * xf).sin(),
                (two_pi * xf).cos(),
                (two_pi * yf).sin(),
                (two_pi * yf).cos(),
            );
            let zf = if D::DIMENSIONS == 3 {
                (two_pi * ((R::from_f64(z as f64) + half) / R::from_f64(g.nz as f64))).cos()
            } else {
                one
            };
            [
                one,
                amplitude * sx * cy * zf,
                zero - amplitude * cx * sy * zf,
                zero,
            ]
        });
    }
    pub fn step(&mut self, count: usize) {
        for _ in 0..count {
            self.step_once();
        }
    }
    fn moments(&self, fin: &[R]) -> (R, R, R, R) {
        let zero = R::from_f64(0.0);
        let (mut rho, mut mx, mut my, mut mz) = (zero, zero, zero, zero);
        for (d, v) in fin.iter().enumerate() {
            rho += *v;
            mx += *v * R::from_f64(f64::from(D::cx(d)));
            my += *v * R::from_f64(f64::from(D::cy(d)));
            mz += *v * R::from_f64(f64::from(D::cz(d)));
        }
        (rho, mx, my, mz)
    }
    fn step_once(&mut self) {
        let omega = R::from_f64(1.0) / self.config.tau;
        let zero = R::from_f64(0.0);
        let cells = self.cells;
        let mut fin = vec![zero; D::Q];
        for n in 0..cells {
            let (x, y, z) = self.grid.coords(n);
            for (d, v) in fin.iter_mut().enumerate() {
                *v = self.f[d * cells + self.grid.step::<D>(x, y, z, d, -1)];
            }
            let (rho, mx, my, mz) = self.moments(&fin);
            let positive = rho > zero;
            let (ux, uy, uz) = if positive {
                (mx / rho, my / rho, mz / rho)
            } else {
                (zero, zero, zero)
            };
            for d in 0..D::Q {
                let feq = typed_equilibrium::<D, R>(d, rho, ux, uy, uz);
                self.next[d * cells + n] = fin[d] - omega * (fin[d] - feq);
            }
        }
        std::mem::swap(&mut self.f, &mut self.next);
    }
    pub fn compute_macroscopic(&self) -> TypedMacroscopicFields<R> {
        let zero = R::from_f64(0.0);
        let mut m = TypedMacroscopicFields {
            rho: vec![zero; self.cells],
            ux: vec![zero; self.cells],
            uy: vec![zero; self.cells],
            uz: vec![zero; self.cells],
        };
        for n in 0..self.cells {
            let fin: Vec<R> = (0..D::Q).map(|d| self.f[d * self.cells + n]).collect();
            let (rho, mx, my, mz) = self.moments(&fin);
            m.rho[n] = rho;
            if rho > zero {
                m.ux[n] = mx / rho;
                m.uy[n] = my / rho;
                m.uz[n] = mz / rho;
            }
        }
        m
    }
    pub fn mass(&self) -> f64 {
        let mut sum = 0.0;
        for n in 0..self.cells {
            for d in 0..D::Q {
                sum += self.f[d * self.cells + n].to_f64();
            }
        }
        sum
    }
    pub fn populations(&self) -> &[R] {
        &self.f
    }
}
