//! Periodic staggered-grid pressure projection and incompressible Navier-Stokes baseline, ported
//! from upstream `cfd/solvers/fvm/{projection2d, incompressible2d}` (MIT, cfd_solvers contributors).
use crate::{
    core::{matrix_free_conjugate_gradient, ConjugateGradientResult},
    require, Result,
};
use std::f64::consts::PI;

/// Staggered periodic layout: `u` on (nx+1)*ny vertical faces, `v` on nx*(ny+1) horizontal
/// faces, `p` on nx*ny cell centres, all row-major in `j`.
#[derive(Clone, Copy, Debug)]
struct Grid {
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
}

impl Grid {
    fn u(&self, i: usize, j: usize) -> usize {
        j * (self.nx + 1) + i
    }
    fn v(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    fn p(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    fn cells(&self) -> usize {
        self.nx * self.ny
    }
    fn enforce_periodic(&self, u: &mut [f64], v: &mut [f64]) {
        for j in 0..self.ny {
            u[self.u(self.nx, j)] = u[self.u(0, j)];
        }
        for i in 0..self.nx {
            v[self.v(i, self.ny)] = v[self.v(i, 0)];
        }
    }
    fn divergence(&self, u: &[f64], v: &[f64], i: usize, j: usize) -> f64 {
        (u[self.u(i + 1, j)] - u[self.u(i, j)]) / self.dx
            + (v[self.v(i, j + 1)] - v[self.v(i, j)]) / self.dy
    }
    fn divergence_l2(&self, u: &[f64], v: &[f64]) -> f64 {
        let mut sum = 0.0;
        for n in 0..self.cells() {
            let d = self.divergence(u, v, n % self.nx, n / self.nx);
            sum += d * d;
        }
        (sum / self.cells() as f64).sqrt()
    }
    fn kinetic_energy(&self, u: &[f64], v: &[f64]) -> f64 {
        let mut sum = 0.0;
        for n in 0..self.cells() {
            let (i, j) = (n % self.nx, n / self.nx);
            let uc = 0.5 * (u[self.u(i, j)] + u[self.u(i + 1, j)]);
            let vc = 0.5 * (v[self.v(i, j)] + v[self.v(i, j + 1)]);
            sum += 0.5 * (uc * uc + vc * vc);
        }
        sum * self.dx * self.dy
    }
    fn negative_laplacian(&self, input: &[f64], out: &mut [f64]) {
        let (idx2, idy2) = (1.0 / (self.dx * self.dx), 1.0 / (self.dy * self.dy));
        for (n, o) in out.iter_mut().enumerate() {
            let (i, j) = (n % self.nx, n / self.nx);
            let im = (i + self.nx - 1) % self.nx;
            let ip = (i + 1) % self.nx;
            let jm = (j + self.ny - 1) % self.ny;
            let jp = (j + 1) % self.ny;
            *o = 2.0 * (idx2 + idy2) * input[n]
                - idx2 * (input[self.p(im, j)] + input[self.p(ip, j)])
                - idy2 * (input[self.p(i, jm)] + input[self.p(i, jp)]);
        }
    }
    /// Solves `-lap p = rhs` with zero initial guess, then removes the mean of `p`.
    fn solve_pressure(
        &self,
        rhs: &[f64],
        p: &mut [f64],
        iterations: usize,
        tolerance: f64,
    ) -> Result<ConjugateGradientResult> {
        p.fill(0.0);
        let r = matrix_free_conjugate_gradient(
            rhs,
            p,
            |a, b| self.negative_laplacian(a, b),
            iterations,
            tolerance,
        )?;
        let mean = p.iter().sum::<f64>() / p.len() as f64;
        for x in p.iter_mut() {
            *x -= mean;
        }
        Ok(r)
    }
    fn taylor_green(
        &self,
        lx: f64,
        ly: f64,
        amplitude: f64,
        u: &mut [f64],
        v: &mut [f64],
        all_u_columns: bool,
    ) {
        let (kx, ky) = (2.0 * PI / lx, 2.0 * PI / ly);
        let iu = if all_u_columns { self.nx + 1 } else { self.nx };
        for j in 0..self.ny {
            let y = (j as f64 + 0.5) * self.dy;
            for i in 0..iu {
                let x = i as f64 * self.dx;
                u[self.u(i, j)] = amplitude * (kx * x).sin() * (ky * y).cos();
            }
        }
        let jv = if all_u_columns { self.ny + 1 } else { self.ny };
        for j in 0..jv {
            let y = j as f64 * self.dy;
            for i in 0..self.nx {
                let x = (i as f64 + 0.5) * self.dx;
                v[self.v(i, j)] = -amplitude * (kx * x).cos() * (ky * y).sin();
            }
        }
        self.enforce_periodic(u, v);
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Projection2DConfig {
    pub nx: usize,
    pub ny: usize,
    pub lx: f64,
    pub ly: f64,
    pub density: f64,
    pub dt: f64,
    pub pressure_iterations: usize,
    pub pressure_tolerance: f64,
}

impl Default for Projection2DConfig {
    fn default() -> Self {
        Self {
            nx: 64,
            ny: 64,
            lx: 1.0,
            ly: 1.0,
            density: 1.0,
            dt: 1e-3,
            pressure_iterations: 600,
            pressure_tolerance: 1e-10,
        }
    }
}

/// Periodic staggered pressure projection with matrix-free CG.
#[derive(Clone, Debug)]
pub struct Projection2D {
    config: Projection2DConfig,
    grid: Grid,
    u: Vec<f64>,
    v: Vec<f64>,
    p: Vec<f64>,
    pressure_result: ConjugateGradientResult,
}

impl Projection2D {
    pub fn new(config: Projection2DConfig) -> Result<Self> {
        let c = config;
        require(
            c.nx >= 4
                && c.ny >= 4
                && c.lx > 0.0
                && c.ly > 0.0
                && c.density > 0.0
                && c.dt > 0.0
                && c.pressure_iterations > 0
                && c.pressure_tolerance > 0.0,
            "Projection2D requires nx,ny>=4 and positive physical/numerical parameters",
        )?;
        let grid = Grid {
            nx: c.nx,
            ny: c.ny,
            dx: c.lx / c.nx as f64,
            dy: c.ly / c.ny as f64,
        };
        Ok(Self {
            config,
            grid,
            u: vec![0.0; (c.nx + 1) * c.ny],
            v: vec![0.0; c.nx * (c.ny + 1)],
            p: vec![0.0; c.nx * c.ny],
            pressure_result: ConjugateGradientResult::default(),
        })
    }
    /// `u = A sin(kx x)`, `v = A sin(ky y)`: a purely divergent field.
    pub fn initialize_divergent(&mut self, amplitude: f64) {
        let g = self.grid;
        let (kx, ky) = (2.0 * PI / self.config.lx, 2.0 * PI / self.config.ly);
        for j in 0..g.ny {
            for i in 0..=g.nx {
                self.u[g.u(i, j)] = amplitude * (kx * i as f64 * g.dx).sin();
            }
        }
        for j in 0..=g.ny {
            for i in 0..g.nx {
                self.v[g.v(i, j)] = amplitude * (ky * j as f64 * g.dy).sin();
            }
        }
        self.p.fill(0.0);
        g.enforce_periodic(&mut self.u, &mut self.v);
    }
    pub fn initialize_taylor_green(&mut self, amplitude: f64) {
        let g = self.grid;
        g.taylor_green(
            self.config.lx,
            self.config.ly,
            amplitude,
            &mut self.u,
            &mut self.v,
            true,
        );
        self.p.fill(0.0);
    }
    pub fn project(&mut self) -> Result<()> {
        let g = self.grid;
        let c = self.config;
        let scale = c.density / c.dt;
        let (idx, idy) = (1.0 / g.dx, 1.0 / g.dy);
        let rhs: Vec<f64> = (0..g.cells())
            .map(|n| {
                let (i, j) = (n % g.nx, n / g.nx);
                -(scale
                    * ((self.u[g.u(i + 1, j)] - self.u[g.u(i, j)]) * idx
                        + (self.v[g.v(i, j + 1)] - self.v[g.v(i, j)]) * idy))
            })
            .collect();
        self.pressure_result = g.solve_pressure(
            &rhs,
            &mut self.p,
            c.pressure_iterations,
            c.pressure_tolerance,
        )?;
        let correction = c.dt / c.density;
        for n in 0..g.cells() {
            let (i, j) = (n % g.nx, n / g.nx);
            let im = (i + g.nx - 1) % g.nx;
            self.u[g.u(i, j)] -= correction * (self.p[g.p(i, j)] - self.p[g.p(im, j)]) * idx;
        }
        for n in 0..g.cells() {
            let (i, j) = (n % g.nx, n / g.nx);
            let jm = (j + g.ny - 1) % g.ny;
            self.v[g.v(i, j)] -= correction * (self.p[g.p(i, j)] - self.p[g.p(i, jm)]) * idy;
        }
        g.enforce_periodic(&mut self.u, &mut self.v);
        Ok(())
    }
    pub fn divergence_l2(&self) -> f64 {
        self.grid.divergence_l2(&self.u, &self.v)
    }
    pub fn kinetic_energy(&self) -> f64 {
        self.grid.kinetic_energy(&self.u, &self.v)
    }
    pub fn u_faces(&self) -> &[f64] {
        &self.u
    }
    pub fn v_faces(&self) -> &[f64] {
        &self.v
    }
    pub fn pressure(&self) -> &[f64] {
        &self.p
    }
    pub fn dx(&self) -> f64 {
        self.grid.dx
    }
    pub fn dy(&self) -> f64 {
        self.grid.dy
    }
    pub fn config(&self) -> &Projection2DConfig {
        &self.config
    }
    pub fn pressure_result(&self) -> ConjugateGradientResult {
        self.pressure_result
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Incompressible2DConfig {
    pub nx: usize,
    pub ny: usize,
    pub lx: f64,
    pub ly: f64,
    pub density: f64,
    pub kinematic_viscosity: f64,
    pub dt: f64,
    pub pressure_iterations: usize,
    pub pressure_tolerance: f64,
}

impl Default for Incompressible2DConfig {
    fn default() -> Self {
        Self {
            nx: 64,
            ny: 64,
            lx: 1.0,
            ly: 1.0,
            density: 1.0,
            kinematic_viscosity: 1e-2,
            dt: 2.5e-4,
            pressure_iterations: 800,
            pressure_tolerance: 1e-10,
        }
    }
}

/// Periodic staggered incompressible Navier-Stokes: flux-form convection, centred diffusion,
/// explicit predictor and CG pressure projection.
#[derive(Clone, Debug)]
pub struct Incompressible2D {
    config: Incompressible2DConfig,
    grid: Grid,
    time: f64,
    steps: usize,
    u: Vec<f64>,
    v: Vec<f64>,
    u_star: Vec<f64>,
    v_star: Vec<f64>,
    p: Vec<f64>,
    pressure_result: ConjugateGradientResult,
}

fn wrap(value: isize, n: usize) -> usize {
    value.rem_euclid(n as isize) as usize
}

impl Incompressible2D {
    pub fn new(config: Incompressible2DConfig) -> Result<Self> {
        let c = config;
        require(
            c.nx >= 8
                && c.ny >= 8
                && c.lx > 0.0
                && c.ly > 0.0
                && c.density > 0.0
                && c.kinematic_viscosity >= 0.0
                && c.dt > 0.0
                && c.pressure_iterations > 0
                && c.pressure_tolerance > 0.0,
            "Incompressible2D requires nx,ny>=8 and valid positive numerical parameters",
        )?;
        let grid = Grid {
            nx: c.nx,
            ny: c.ny,
            dx: c.lx / c.nx as f64,
            dy: c.ly / c.ny as f64,
        };
        let (nu, nv) = ((c.nx + 1) * c.ny, c.nx * (c.ny + 1));
        Ok(Self {
            config,
            grid,
            time: 0.0,
            steps: 0,
            u: vec![0.0; nu],
            v: vec![0.0; nv],
            u_star: vec![0.0; nu],
            v_star: vec![0.0; nv],
            p: vec![0.0; c.nx * c.ny],
            pressure_result: ConjugateGradientResult::default(),
        })
    }
    pub fn initialize_taylor_green(&mut self, amplitude: f64) {
        let g = self.grid;
        g.taylor_green(
            self.config.lx,
            self.config.ly,
            amplitude,
            &mut self.u,
            &mut self.v,
            false,
        );
        self.p.fill(0.0);
        self.time = 0.0;
        self.steps = 0;
    }
    fn up(&self, i: isize, j: isize) -> f64 {
        self.u[self.grid.u(wrap(i, self.grid.nx), wrap(j, self.grid.ny))]
    }
    fn vp(&self, i: isize, j: isize) -> f64 {
        self.v[self.grid.v(wrap(i, self.grid.nx), wrap(j, self.grid.ny))]
    }
    fn momentum_predictor(&mut self) {
        let g = self.grid;
        let (idx, idy) = (1.0 / g.dx, 1.0 / g.dy);
        let (idx2, idy2) = (idx * idx, idy * idy);
        let nu = self.config.kinematic_viscosity;
        let dt = self.config.dt;
        for n in 0..g.cells() {
            let (i, j) = ((n % g.nx) as isize, (n / g.nx) as isize);
            let uc = self.up(i, j);
            let ue = 0.5 * (uc + self.up(i + 1, j));
            let uw = 0.5 * (self.up(i - 1, j) + uc);
            let un = 0.5 * (uc + self.up(i, j + 1));
            let us = 0.5 * (self.up(i, j - 1) + uc);
            let vn = 0.5 * (self.vp(i - 1, j + 1) + self.vp(i, j + 1));
            let vs = 0.5 * (self.vp(i - 1, j) + self.vp(i, j));
            let conv = (ue * ue - uw * uw) * idx + (un * vn - us * vs) * idy;
            let lap = (self.up(i + 1, j) - 2.0 * uc + self.up(i - 1, j)) * idx2
                + (self.up(i, j + 1) - 2.0 * uc + self.up(i, j - 1)) * idy2;
            self.u_star[g.u(i as usize, j as usize)] = uc + dt * (-conv + nu * lap);
        }
        for n in 0..g.cells() {
            let (i, j) = ((n % g.nx) as isize, (n / g.nx) as isize);
            let vc = self.vp(i, j);
            let ve = 0.5 * (vc + self.vp(i + 1, j));
            let vw = 0.5 * (self.vp(i - 1, j) + vc);
            let vn = 0.5 * (vc + self.vp(i, j + 1));
            let vs = 0.5 * (self.vp(i, j - 1) + vc);
            let ue = 0.5 * (self.up(i + 1, j - 1) + self.up(i + 1, j));
            let uw = 0.5 * (self.up(i, j - 1) + self.up(i, j));
            let conv = (ue * ve - uw * vw) * idx + (vn * vn - vs * vs) * idy;
            let lap = (self.vp(i + 1, j) - 2.0 * vc + self.vp(i - 1, j)) * idx2
                + (self.vp(i, j + 1) - 2.0 * vc + self.vp(i, j - 1)) * idy2;
            self.v_star[g.v(i as usize, j as usize)] = vc + dt * (-conv + nu * lap);
        }
        g.enforce_periodic(&mut self.u_star, &mut self.v_star);
    }
    fn project_star_velocity(&mut self) -> Result<()> {
        let g = self.grid;
        let c = self.config;
        let (idx, idy) = (1.0 / g.dx, 1.0 / g.dy);
        let scale = c.density / c.dt;
        let mut rhs: Vec<f64> = (0..g.cells())
            .map(|n| {
                let (i, j) = (n % g.nx, n / g.nx);
                let div = (self.u_star[g.u(i + 1, j)] - self.u_star[g.u(i, j)]) * idx
                    + (self.v_star[g.v(i, j + 1)] - self.v_star[g.v(i, j)]) * idy;
                -scale * div
            })
            .collect();
        // Keep the right-hand side in the compatible mean-zero subspace of the periodic Poisson operator.
        let mean = rhs.iter().sum::<f64>() / rhs.len() as f64;
        for r in &mut rhs {
            *r -= mean;
        }
        self.pressure_result = g.solve_pressure(
            &rhs,
            &mut self.p,
            c.pressure_iterations,
            c.pressure_tolerance,
        )?;
        let correction = c.dt / c.density;
        for n in 0..g.cells() {
            let (i, j) = (n % g.nx, n / g.nx);
            let im = (i + g.nx - 1) % g.nx;
            let jm = (j + g.ny - 1) % g.ny;
            self.u[g.u(i, j)] = self.u_star[g.u(i, j)]
                - correction * (self.p[g.p(i, j)] - self.p[g.p(im, j)]) * idx;
            self.v[g.v(i, j)] = self.v_star[g.v(i, j)]
                - correction * (self.p[g.p(i, j)] - self.p[g.p(i, jm)]) * idy;
        }
        g.enforce_periodic(&mut self.u, &mut self.v);
        Ok(())
    }
    pub fn step(&mut self) -> Result<()> {
        self.momentum_predictor();
        self.project_star_velocity()?;
        self.time += self.config.dt;
        self.steps += 1;
        Ok(())
    }
    pub fn run(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }
    pub fn divergence_l2(&self) -> f64 {
        self.grid.divergence_l2(&self.u, &self.v)
    }
    pub fn kinetic_energy(&self) -> f64 {
        self.grid.kinetic_energy(&self.u, &self.v)
    }
    /// RMS error of the face velocities against the decaying Taylor-Green solution.
    pub fn taylor_green_velocity_error(&self, initial_amplitude: f64) -> f64 {
        let g = self.grid;
        let c = self.config;
        let (kx, ky) = (2.0 * PI / c.lx, 2.0 * PI / c.ly);
        let amplitude =
            initial_amplitude * (-c.kinematic_viscosity * (kx * kx + ky * ky) * self.time).exp();
        let mut sum = 0.0;
        for n in 0..g.cells() {
            let (i, j) = (n % g.nx, n / g.nx);
            let (xu, yu) = (i as f64 * g.dx, (j as f64 + 0.5) * g.dy);
            let (xv, yv) = ((i as f64 + 0.5) * g.dx, j as f64 * g.dy);
            let du = self.u[g.u(i, j)] - amplitude * (kx * xu).sin() * (ky * yu).cos();
            let dv = self.v[g.v(i, j)] + amplitude * (kx * xv).cos() * (ky * yv).sin();
            sum += du * du + dv * dv;
        }
        (sum / (2.0 * g.cells() as f64)).sqrt()
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn steps(&self) -> usize {
        self.steps
    }
    pub fn pressure_result(&self) -> ConjugateGradientResult {
        self.pressure_result
    }
    pub fn u_faces(&self) -> &[f64] {
        &self.u
    }
    pub fn v_faces(&self) -> &[f64] {
        &self.v
    }
    pub fn pressure(&self) -> &[f64] {
        &self.p
    }
    pub fn config(&self) -> &Incompressible2DConfig {
        &self.config
    }
}
