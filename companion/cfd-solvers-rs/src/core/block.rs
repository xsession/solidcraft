//! Dense blocks in CSR and block-Jacobi-preconditioned matrix-free GMRES.
use super::{
    restarted_gmres, CsrBuilder, CsrMatrix, IterativeSolverResult, LinearConfig, Preconditioner,
};
use crate::{finite_slice, positive, require, Result};
use std::collections::BTreeMap;
pub type DenseBlock<const B: usize> = [[f64; B]; B];
#[derive(Clone, Debug)]
pub struct BlockCsrMatrix<const B: usize> {
    rows: usize,
    cols: usize,
    offsets: Vec<usize>,
    columns: Vec<usize>,
    values: Vec<DenseBlock<B>>,
}
impl<const B: usize> BlockCsrMatrix<B> {
    pub fn new(
        rows: usize,
        cols: usize,
        offsets: Vec<usize>,
        columns: Vec<usize>,
        values: Vec<DenseBlock<B>>,
    ) -> Result<Self> {
        require(
            B > 0 && rows.checked_mul(B).is_some() && cols.checked_mul(B).is_some(),
            "invalid block size or dimension overflow",
        )?;
        require(
            rows.checked_add(1) == Some(offsets.len())
                && columns.len() == values.len()
                && offsets.first() == Some(&0)
                && offsets.last() == Some(&values.len()),
            "block CSR shape",
        )?;
        require(
            offsets.iter().all(|&v| v <= values.len()) && offsets.windows(2).all(|v| v[0] <= v[1]),
            "block CSR offsets",
        )?;
        for r in 0..rows {
            let c = &columns[offsets[r]..offsets[r + 1]];
            require(
                c.iter().all(|&v| v < cols) && c.windows(2).all(|v| v[0] < v[1]),
                "block CSR columns",
            )?;
        }
        require(
            values
                .iter()
                .flat_map(|b| b.iter().flatten())
                .all(|v| v.is_finite()),
            "nonfinite block CSR value",
        )?;
        Ok(Self {
            rows,
            cols,
            offsets,
            columns,
            values,
        })
    }
    pub fn block_rows(&self) -> usize {
        self.rows
    }
    pub fn block_cols(&self) -> usize {
        self.cols
    }
    pub fn block_nonzeros(&self) -> usize {
        self.values.len()
    }
    pub fn row_offsets(&self) -> &[usize] {
        &self.offsets
    }
    pub fn columns(&self) -> &[usize] {
        &self.columns
    }
    pub fn values(&self) -> &[DenseBlock<B>] {
        &self.values
    }
    pub fn multiply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
        require(
            x.len() == self.cols * B && y.len() == self.rows * B,
            "block CSR multiply size",
        )?;
        finite_slice(x, "block vector")?;
        y.fill(0.0);
        for r in 0..self.rows {
            for k in self.offsets[r]..self.offsets[r + 1] {
                for i in 0..B {
                    for j in 0..B {
                        y[r * B + i] += self.values[k][i][j] * x[self.columns[k] * B + j];
                    }
                }
            }
        }
        finite_slice(y, "block product")
    }
    pub fn apply(&self, x: &[f64]) -> Result<Vec<f64>> {
        let mut y = vec![0.0; self.rows * B];
        self.multiply(x, &mut y)?;
        Ok(y)
    }
    pub fn scalar_csr(&self) -> Result<CsrMatrix> {
        let mut b = CsrBuilder::new(self.rows * B, self.cols * B);
        for r in 0..self.rows {
            for k in self.offsets[r]..self.offsets[r + 1] {
                for i in 0..B {
                    for j in 0..B {
                        b.add(r * B + i, self.columns[k] * B + j, self.values[k][i][j])?;
                    }
                }
            }
        }
        b.build(0.0)
    }
    pub fn diagonal_blocks(&self) -> Vec<DenseBlock<B>> {
        (0..self.rows)
            .map(|r| {
                let start = self.offsets[r];
                self.columns[start..self.offsets[r + 1]]
                    .binary_search(&r)
                    .map(|i| self.values[start + i])
                    .unwrap_or([[0.0; B]; B])
            })
            .collect()
    }
}
#[derive(Clone, Debug)]
pub struct BlockCsrBuilder<const B: usize> {
    rows: usize,
    cols: usize,
    entries: Vec<BTreeMap<usize, DenseBlock<B>>>,
}
impl<const B: usize> BlockCsrBuilder<B> {
    pub fn new(rows: usize, cols: usize) -> Result<Self> {
        require(
            B > 0 && rows.checked_mul(B).is_some() && cols.checked_mul(B).is_some(),
            "invalid block dimensions",
        )?;
        Ok(Self {
            rows,
            cols,
            entries: vec![BTreeMap::new(); rows],
        })
    }
    pub fn add(&mut self, row: usize, col: usize, block: DenseBlock<B>) -> Result<()> {
        require(
            row < self.rows && col < self.cols && block.iter().flatten().all(|v| v.is_finite()),
            "invalid block add",
        )?;
        let e = self.entries[row].entry(col).or_insert([[0.0; B]; B]);
        for i in 0..B {
            for j in 0..B {
                e[i][j] += block[i][j];
            }
        }
        Ok(())
    }
    pub fn build(&self) -> Result<BlockCsrMatrix<B>> {
        let mut offsets = vec![0];
        let mut columns = vec![];
        let mut values = vec![];
        for row in &self.entries {
            for (&col, &block) in row {
                columns.push(col);
                values.push(block);
            }
            offsets.push(values.len());
        }
        BlockCsrMatrix::new(self.rows, self.cols, offsets, columns, values)
    }
}
#[derive(Clone, Debug)]
pub struct BlockJacobiPreconditioner<const B: usize> {
    inverse: Vec<DenseBlock<B>>,
}
impl<const B: usize> BlockJacobiPreconditioner<B> {
    pub fn new(diagonal: &[DenseBlock<B>], pivot_floor: f64) -> Result<Self> {
        require(B > 0 && B.checked_mul(2).is_some(), "invalid block size")?;
        positive(pivot_floor, "block pivot floor")?;
        let mut inverse = vec![];
        for block in diagonal {
            require(
                block.iter().flatten().all(|v| v.is_finite()),
                "nonfinite block diagonal",
            )?;
            let mut a = vec![vec![0.0; 2 * B]; B];
            for i in 0..B {
                for j in 0..B {
                    a[i][j] = block[i][j];
                }
                a[i][B + i] = 1.0;
            }
            for col in 0..B {
                let mut pivot = col;
                for r in col + 1..B {
                    if a[r][col].abs() > a[pivot][col].abs() {
                        pivot = r;
                    }
                }
                require(
                    a[pivot][col].abs() >= pivot_floor,
                    "singular block-Jacobi diagonal",
                )?;
                a.swap(col, pivot);
                let scale = a[col][col];
                for j in 0..2 * B {
                    a[col][j] /= scale;
                }
                for r in 0..B {
                    if r == col {
                        continue;
                    }
                    let f = a[r][col];
                    for j in 0..2 * B {
                        a[r][j] -= f * a[col][j];
                    }
                }
            }
            let mut b = [[0.0; B]; B];
            for i in 0..B {
                for j in 0..B {
                    b[i][j] = a[i][B + j];
                }
            }
            require(
                b.iter().flatten().all(|v| v.is_finite()),
                "block inverse overflow",
            )?;
            inverse.push(b);
        }
        Ok(Self { inverse })
    }
    pub fn blocks(&self) -> usize {
        self.inverse.len()
    }
}
impl<const B: usize> Preconditioner for BlockJacobiPreconditioner<B> {
    fn apply(&self, rhs: &[f64], result: &mut [f64]) -> Result<()> {
        require(
            rhs.len() == self.inverse.len() * B && result.len() == rhs.len(),
            "block-Jacobi vector dimensions",
        )?;
        for (b, a) in self.inverse.iter().enumerate() {
            for i in 0..B {
                result[b * B + i] = (0..B).map(|j| a[i][j] * rhs[b * B + j]).sum();
            }
        }
        finite_slice(result, "block preconditioned vector")
    }
}
pub fn restarted_block_gmres<const B: usize, F>(
    rhs: &[f64],
    x: &mut [f64],
    apply: F,
    diagonal: &[DenseBlock<B>],
    config: LinearConfig,
) -> Result<IterativeSolverResult>
where
    F: Fn(&[f64], &mut [f64]) -> Result<()>,
{
    require(
        !diagonal.is_empty() && rhs.len() == diagonal.len() * B && x.len() == rhs.len(),
        "block GMRES dimensions",
    )?;
    let p = BlockJacobiPreconditioner::new(diagonal, 1e-14)?;
    restarted_gmres(rhs, x, apply, &p, config)
}
