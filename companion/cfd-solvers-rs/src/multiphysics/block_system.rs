//! Block-structured coupled linear system assembled into one CSR matrix and solved by
//! ILU(0)-preconditioned restarted GMRES.
use crate::{
    core::{
        restarted_gmres, CsrBuilder, CsrMatrix, Ilu0Preconditioner, IterativeSolverResult,
        LinearConfig,
    },
    require, Result, SolverError,
};

pub struct BlockCoupledLinearSystem {
    sizes: Vec<usize>,
    offsets: Vec<usize>,
    builder: CsrBuilder,
    rhs: Vec<f64>,
}

impl BlockCoupledLinearSystem {
    pub fn new(block_sizes: Vec<usize>) -> Result<Self> {
        require(!block_sizes.is_empty(), "block system requires blocks")?;
        require(
            block_sizes.iter().all(|s| *s > 0),
            "block size must be positive",
        )?;
        let mut offsets: Vec<usize> = vec![0];
        for s in &block_sizes {
            let last = offsets.last().copied().unwrap_or(0);
            offsets.push(
                last.checked_add(*s).ok_or_else(|| {
                    SolverError::InvalidInput("block system size overflow".into())
                })?,
            );
        }
        let n = offsets.last().copied().unwrap_or(0);
        Ok(Self {
            sizes: block_sizes,
            offsets,
            builder: CsrBuilder::new(n, n),
            rhs: vec![0.0; n],
        })
    }
    fn index(&self, block: usize, local: usize) -> Result<usize> {
        require(
            block < self.sizes.len() && local < self.sizes[block],
            "block system index",
        )?;
        Ok(self.offsets[block] + local)
    }
    pub fn blocks(&self) -> usize {
        self.sizes.len()
    }
    pub fn total_size(&self) -> usize {
        self.offsets.last().copied().unwrap_or(0)
    }
    pub fn add(
        &mut self,
        row_block: usize,
        row: usize,
        col_block: usize,
        col: usize,
        value: f64,
    ) -> Result<()> {
        let (r, c) = (self.index(row_block, row)?, self.index(col_block, col)?);
        self.builder.add(r, c, value)
    }
    pub fn add_rhs(&mut self, block: usize, row: usize, value: f64) -> Result<()> {
        let i = self.index(block, row)?;
        require(value.is_finite(), "nonfinite block right-hand side")?;
        self.rhs[i] += value;
        Ok(())
    }
    pub fn rhs(&self) -> &[f64] {
        &self.rhs
    }
    pub fn matrix(&self) -> Result<CsrMatrix> {
        self.builder.build(0.0)
    }
    pub fn solve(
        &self,
        solution: &mut [f64],
        max_iterations: usize,
        restart: usize,
        relative_tolerance: f64,
    ) -> Result<IterativeSolverResult> {
        require(
            solution.len() == self.total_size(),
            "block solution size mismatch",
        )?;
        let a = self.matrix()?;
        let ilu = Ilu0Preconditioner::new(&a, 1e-12)?;
        restarted_gmres(
            &self.rhs,
            solution,
            |v, y| a.multiply(v, y),
            &ilu,
            LinearConfig {
                max_iterations,
                restart,
                relative_tolerance,
            },
        )
    }
}
