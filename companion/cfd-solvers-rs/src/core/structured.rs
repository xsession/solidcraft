use crate::{finite_slice, positive, require, Result, SolverError};
#[derive(Clone, Copy, Debug)]
pub struct Grid2D {
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
}
impl Grid2D {
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Result<Self> {
        require(
            nx >= 2 && ny >= 2 && nx.checked_mul(ny).is_some(),
            "invalid 2D grid dimensions",
        )?;
        positive(dx, "dx")?;
        positive(dy, "dy")?;
        Ok(Self { nx, ny, dx, dy })
    }
    pub fn nx(self) -> usize {
        self.nx
    }
    pub fn ny(self) -> usize {
        self.ny
    }
    pub fn dx(self) -> f64 {
        self.dx
    }
    pub fn dy(self) -> f64 {
        self.dy
    }
    pub fn cells(self) -> usize {
        self.nx * self.ny
    }
    pub fn index(self, x: usize, y: usize) -> Result<usize> {
        require(x < self.nx && y < self.ny, "grid index out of range")?;
        Ok(y * self.nx + x)
    }
}
#[derive(Clone, Copy, Debug)]
pub struct StructuredSevenPointOperator {
    nx: usize,
    ny: usize,
    nz: usize,
    diagonal: f64,
    neighbor: f64,
}
impl StructuredSevenPointOperator {
    pub fn new(nx: usize, ny: usize, nz: usize, diagonal: f64, neighbor: f64) -> Result<Self> {
        require(
            nx > 0
                && ny > 0
                && nz > 0
                && nx.checked_mul(ny).and_then(|v| v.checked_mul(nz)).is_some()
                && diagonal.is_finite()
                && neighbor.is_finite(),
            "invalid stencil dimensions/coefficients",
        )?;
        Ok(Self {
            nx,
            ny,
            nz,
            diagonal,
            neighbor,
        })
    }
    pub fn nx(self) -> usize {
        self.nx
    }
    pub fn ny(self) -> usize {
        self.ny
    }
    pub fn nz(self) -> usize {
        self.nz
    }
    pub fn cells(self) -> usize {
        self.nx * self.ny * self.nz
    }
    pub fn diagonal(self) -> f64 {
        self.diagonal
    }
    pub fn neighbor(self) -> f64 {
        self.neighbor
    }
    pub fn estimated_apply_bytes(self) -> Result<usize> {
        self.cells()
            .checked_mul(64)
            .ok_or_else(|| SolverError::InvalidInput("traffic estimate overflow".into()))
    }
    pub fn estimated_fused_residual_bytes(self) -> Result<usize> {
        self.cells()
            .checked_mul(72)
            .ok_or_else(|| SolverError::InvalidInput("traffic estimate overflow".into()))
    }
    fn value(self, x: &[f64], i: usize) -> f64 {
        let plane = self.nx * self.ny;
        let cx = i % self.nx;
        let cy = (i / self.nx) % self.ny;
        let cz = i / plane;
        let mut a = self.diagonal * x[i];
        if cx > 0 {
            a += self.neighbor * x[i - 1];
        }
        if cx + 1 < self.nx {
            a += self.neighbor * x[i + 1];
        }
        if cy > 0 {
            a += self.neighbor * x[i - self.nx];
        }
        if cy + 1 < self.ny {
            a += self.neighbor * x[i + self.nx];
        }
        if cz > 0 {
            a += self.neighbor * x[i - plane];
        }
        if cz + 1 < self.nz {
            a += self.neighbor * x[i + plane];
        }
        a
    }
    pub fn apply(self, x: &[f64], y: &mut [f64]) -> Result<()> {
        require(
            x.len() == self.cells() && y.len() == x.len(),
            "stencil dimensions",
        )?;
        finite_slice(x, "stencil input")?;
        for (i, v) in y.iter_mut().enumerate() {
            *v = self.value(x, i);
        }
        finite_slice(y, "stencil output")
    }
    pub fn apply_residual_l2(self, rhs: &[f64], x: &[f64], r: &mut [f64]) -> Result<f64> {
        require(
            rhs.len() == self.cells() && x.len() == rhs.len() && r.len() == x.len(),
            "stencil residual dimensions",
        )?;
        finite_slice(rhs, "stencil RHS")?;
        finite_slice(x, "stencil vector")?;
        let mut sum = 0.0;
        for (i, v) in r.iter_mut().enumerate() {
            *v = rhs[i] - self.value(x, i);
            sum += *v * *v;
        }
        require(sum.is_finite(), "stencil residual overflow")?;
        Ok(sum.sqrt())
    }
}
