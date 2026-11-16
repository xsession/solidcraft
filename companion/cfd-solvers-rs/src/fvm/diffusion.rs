use crate::{positive, require, Result};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Diffusion2DConfig {
    pub nx: usize,
    pub ny: usize,
    pub lx: f64,
    pub ly: f64,
    pub diffusivity: f64,
}
impl Default for Diffusion2DConfig {
    fn default() -> Self {
        Self {
            nx: 128,
            ny: 128,
            lx: 1.0,
            ly: 1.0,
            diffusivity: 1.0,
        }
    }
}
#[derive(Clone, Debug)]
pub struct Diffusion2D {
    config: Diffusion2DConfig,
    phi: Vec<f64>,
    source: Vec<f64>,
    boundary: [f64; 4],
}
impl Diffusion2D {
    pub fn new(c: Diffusion2DConfig) -> Result<Self> {
        require(
            c.nx >= 3 && c.ny >= 3,
            "diffusion grid needs at least 3 by 3 points",
        )?;
        positive(c.lx, "lx")?;
        positive(c.ly, "ly")?;
        positive(c.diffusivity, "diffusivity")?;
        let n = c
            .nx
            .checked_mul(c.ny)
            .ok_or_else(|| crate::SolverError::InvalidInput("diffusion size overflow".into()))?;
        Ok(Self {
            config: c,
            phi: vec![0.0; n],
            source: vec![0.0; n],
            boundary: [0.0; 4],
        })
    }
    pub fn config(&self) -> Diffusion2DConfig {
        self.config
    }
    pub fn field(&self) -> &[f64] {
        &self.phi
    }
    pub fn set_source(&mut self, value: f64) -> Result<()> {
        require(value.is_finite(), "diffusion source must be finite")?;
        self.source.fill(value);
        Ok(())
    }
    pub fn set_dirichlet(&mut self, left: f64, right: f64, bottom: f64, top: f64) -> Result<()> {
        require(
            [left, right, bottom, top].iter().all(|v| v.is_finite()),
            "boundary must be finite",
        )?;
        self.boundary = [left, right, bottom, top];
        let (nx, ny) = (self.config.nx, self.config.ny);
        for y in 0..ny {
            self.phi[y * nx] = left;
            self.phi[y * nx + nx - 1] = right;
        }
        for x in 0..nx {
            self.phi[x] = bottom;
            self.phi[(ny - 1) * nx + x] = top;
        }
        Ok(())
    }
    pub fn iterate(&mut self, iterations: usize) -> Result<f64> {
        let c = self.config;
        let idx2 = ((c.nx - 1) as f64 / c.lx).powi(2);
        let idy2 = ((c.ny - 1) as f64 / c.ly).powi(2);
        let denom = 2.0 * c.diffusivity * (idx2 + idy2);
        let mut next = self.phi.clone();
        let mut delta: f64 = 0.0;
        for _ in 0..iterations {
            delta = 0.0;
            for y in 1..c.ny - 1 {
                for x in 1..c.nx - 1 {
                    let i = y * c.nx + x;
                    next[i] = (c.diffusivity
                        * (idx2 * (self.phi[i + 1] + self.phi[i - 1])
                            + idy2 * (self.phi[i + c.nx] + self.phi[i - c.nx]))
                        + self.source[i])
                        / denom;
                    require(next[i].is_finite(), "diffusion iterate became nonfinite")?;
                    delta = delta.max((next[i] - self.phi[i]).abs());
                }
            }
            std::mem::swap(&mut self.phi, &mut next);
        }
        Ok(delta)
    }
}
