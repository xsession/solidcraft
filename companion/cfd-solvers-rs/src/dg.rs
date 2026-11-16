//! Nodal Gauss-Lobatto DG advection and high-order Chebyshev-Lobatto bases, ported from
//! upstream `cfd/dg/{advection1d, high_order_fem}` (MIT, cfd_solvers contributors).
use crate::{require, Result, SolverError};

/// `(P_n(x), P_n'(x))` by the three-term recurrence (`x = +-1` gives a non-finite derivative).
pub fn legendre_value_derivative(order: usize, x: f64) -> (f64, f64) {
    let mut previous = 1.0;
    if order == 0 {
        return (previous, 0.0);
    }
    let mut current = x;
    for degree in 2..=order {
        let d = degree as f64;
        let next = ((2.0 * d - 1.0) * x * current - (d - 1.0) * previous) / d;
        previous = current;
        current = next;
    }
    (
        current,
        order as f64 * (previous - x * current) / (1.0 - x * x),
    )
}

/// Legendre-Gauss-Lobatto nodes on `[-1, 1]` by Newton iteration from Chebyshev guesses.
pub fn gauss_lobatto_nodes(order: usize) -> Vec<f64> {
    if order == 0 {
        return vec![0.0];
    }
    let mut nodes = vec![0.0; order + 1];
    nodes[0] = -1.0;
    nodes[order] = 1.0;
    for (i, node) in nodes.iter_mut().enumerate().take(order).skip(1) {
        let mut x = -(std::f64::consts::PI * i as f64 / order as f64).cos();
        for _ in 0..30 {
            let (value, derivative) = legendre_value_derivative(order, x);
            let second =
                (2.0 * x * derivative - (order * (order + 1)) as f64 * value) / (1.0 - x * x);
            if second.abs() < 1e-14 {
                break;
            }
            let next = x - derivative / second;
            if (next - x).abs() < 1e-14 {
                x = next;
                break;
            }
            x = next.clamp(-1.0 + 1e-12, 1.0 - 1e-12);
        }
        *node = x;
    }
    nodes
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgAdvection1DConfig {
    pub elements: usize,
    pub order: usize,
    pub length: f64,
    pub speed: f64,
}

impl Default for DgAdvection1DConfig {
    fn default() -> Self {
        Self {
            elements: 16,
            order: 3,
            length: 1.0,
            speed: 1.0,
        }
    }
}

/// Periodic nodal (LGL collocated) DG for `u_t + a u_x = 0` with upwind flux and SSP-RK3.
#[derive(Clone, Debug)]
pub struct DgAdvection1D {
    config: DgAdvection1DConfig,
    nodes: Vec<f64>,
    weights: Vec<f64>,
    derivative: Vec<f64>,
    values: Vec<f64>,
}

impl DgAdvection1D {
    pub fn new(config: DgAdvection1DConfig) -> Result<Self> {
        require(
            config.elements > 0
                && config.order > 0
                && config.length > 0.0
                && config.speed.is_finite(),
            "invalid DG advection configuration",
        )?;
        let nodes = gauss_lobatto_nodes(config.order);
        let p = nodes.len();
        let order = config.order as f64;
        let mut weights = vec![0.0; p];
        let mut derivative = vec![0.0; p * p];
        for i in 0..p {
            let (vi, _) = legendre_value_derivative(config.order, nodes[i]);
            weights[i] = 2.0 / (order * (order + 1.0) * vi * vi);
            for j in 0..p {
                if i != j {
                    let (vj, _) = legendre_value_derivative(config.order, nodes[j]);
                    derivative[i * p + j] = vi / (vj * (nodes[i] - nodes[j]));
                }
            }
            let mut diagonal = 0.0;
            for column in 0..p {
                if column != i {
                    diagonal -= derivative[i * p + column];
                }
            }
            derivative[i * p + i] = diagonal;
        }
        let n = config.elements * p;
        Ok(Self {
            config,
            nodes,
            weights,
            derivative,
            values: vec![0.0; n],
        })
    }
    pub fn initialize(&mut self, f: impl Fn(f64) -> f64) {
        let h = self.config.length / self.config.elements as f64;
        let p = self.nodes.len();
        for e in 0..self.config.elements {
            for k in 0..p {
                self.values[e * p + k] = f(h * (e as f64 + 0.5 * (self.nodes[k] + 1.0)));
            }
        }
    }
    fn rhs(&self, input: &[f64], output: &mut [f64]) {
        let p = self.nodes.len();
        let (ne, a) = (self.config.elements, self.config.speed);
        let h = self.config.length / ne as f64;
        for e in 0..ne {
            let off = e * p;
            for k in 0..p {
                let mut d = 0.0;
                for c in 0..p {
                    d += self.derivative[k * p + c] * input[off + c];
                }
                output[off + k] = -a * 2.0 * d / h;
            }
            let prev = if e == 0 { ne - 1 } else { e - 1 };
            let next = if e + 1 == ne { 0 } else { e + 1 };
            let (li, ri) = (input[off], input[off + p - 1]);
            let (ln, rn) = (input[prev * p + p - 1], input[next * p]);
            let lf = if a >= 0.0 { a * ln } else { a * li };
            let rf = if a >= 0.0 { a * ri } else { a * rn };
            output[off] += 2.0 * (lf - a * li) / (h * self.weights[0]);
            output[off + p - 1] -= 2.0 * (rf - a * ri) / (h * self.weights[p - 1]);
        }
    }
    /// One SSP-RK3 step.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(dt > 0.0 && dt.is_finite(), "invalid DG timestep")?;
        let n = self.values.len();
        let mut rhs = vec![0.0; n];
        self.rhs(&self.values, &mut rhs);
        let stage: Vec<f64> = (0..n).map(|i| self.values[i] + dt * rhs[i]).collect();
        self.rhs(&stage, &mut rhs);
        let stage_two: Vec<f64> = (0..n)
            .map(|i| 0.75 * self.values[i] + 0.25 * (stage[i] + dt * rhs[i]))
            .collect();
        self.rhs(&stage_two, &mut rhs);
        for i in 0..n {
            self.values[i] = (self.values[i] + 2.0 * (stage_two[i] + dt * rhs[i])) / 3.0;
        }
        Ok(())
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    pub fn reference_nodes(&self) -> &[f64] {
        &self.nodes
    }
    /// LGL-quadrature integral of the solution.
    pub fn mass(&self) -> f64 {
        let h = self.config.length / self.config.elements as f64;
        let p = self.nodes.len();
        let mut total = 0.0;
        for e in 0..self.config.elements {
            for k in 0..p {
                total += 0.5 * h * self.weights[k] * self.values[e * p + k];
            }
        }
        total
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct OneDimensionalInterval {
    pub left: f64,
    pub right: f64,
    pub level: u32,
}

/// Chebyshev-Lobatto nodal basis with barycentric differentiation matrix and interval bisection.
#[derive(Clone, Debug)]
pub struct HighOrderNodalBasis1D {
    order: u32,
    nodes: Vec<f64>,
    derivative: Vec<f64>,
}

impl HighOrderNodalBasis1D {
    pub fn new(polynomial_order: u32) -> Result<Self> {
        require(
            polynomial_order > 0,
            "high-order basis needs positive order",
        )?;
        let p = polynomial_order as usize + 1;
        let pi = std::f64::consts::PI;
        let mut nodes = vec![0.0; p];
        let mut weights = vec![1.0; p];
        for i in 0..p {
            nodes[i] = -(pi * i as f64 / f64::from(polynomial_order)).cos();
            weights[i] = if i == 0 || i + 1 == p { 0.5 } else { 1.0 };
            if i % 2 != 0 {
                weights[i] = -weights[i];
            }
        }
        let mut derivative = vec![0.0; p * p];
        for i in 0..p {
            for j in 0..p {
                if i != j {
                    derivative[i * p + j] = weights[j] / (weights[i] * (nodes[i] - nodes[j]));
                }
            }
            let mut diagonal = 0.0;
            for j in 0..p {
                if i != j {
                    diagonal += derivative[i * p + j];
                }
            }
            derivative[i * p + i] = -diagonal;
        }
        Ok(Self {
            order: polynomial_order,
            nodes,
            derivative,
        })
    }
    pub fn order(&self) -> u32 {
        self.order
    }
    pub fn nodes(&self) -> &[f64] {
        &self.nodes
    }
    pub fn derivative_matrix(&self) -> &[f64] {
        &self.derivative
    }
    pub fn element_nodes(&self, left: f64, right: f64) -> Result<Vec<f64>> {
        require(
            right > left && left.is_finite() && right.is_finite(),
            "invalid high-order element interval",
        )?;
        Ok(self
            .nodes
            .iter()
            .map(|n| 0.5 * ((right - left) * n + right + left))
            .collect())
    }
    /// Bisects marked intervals (level + 1), keeping order.
    pub fn refine_marked(
        &self,
        mesh: &[OneDimensionalInterval],
        marked: &[usize],
    ) -> Result<Vec<OneDimensionalInterval>> {
        let mut is_marked = vec![false; mesh.len()];
        for &i in marked {
            *is_marked
                .get_mut(i)
                .ok_or_else(|| SolverError::InvalidInput("high-order AMR marker index".into()))? =
                true;
        }
        let mut refined = Vec::with_capacity(mesh.len() + marked.len());
        for (cell, m) in mesh.iter().zip(&is_marked) {
            require(
                cell.right > cell.left && cell.left.is_finite() && cell.right.is_finite(),
                "invalid AMR interval",
            )?;
            if !m {
                refined.push(*cell);
                continue;
            }
            let mid = 0.5 * (cell.left + cell.right);
            refined.push(OneDimensionalInterval {
                left: cell.left,
                right: mid,
                level: cell.level + 1,
            });
            refined.push(OneDimensionalInterval {
                left: mid,
                right: cell.right,
                level: cell.level + 1,
            });
        }
        Ok(refined)
    }
}

/// Matrix-free line derivative action for tensor-product elements.
#[derive(Clone, Debug)]
pub struct SumFactorizedDerivative1D {
    size: usize,
    derivative: Vec<f64>,
}

impl SumFactorizedDerivative1D {
    pub fn new(basis: &HighOrderNodalBasis1D) -> Self {
        Self {
            size: basis.nodes().len(),
            derivative: basis.derivative_matrix().to_vec(),
        }
    }
    pub fn apply(&self, values: &[f64], output: &mut [f64], inverse_length: f64) -> Result<()> {
        require(
            values.len() == self.size && output.len() == self.size && inverse_length.is_finite(),
            "sum-factorized derivative size mismatch",
        )?;
        for (row, o) in output.iter_mut().enumerate() {
            let mut v = 0.0;
            for (c, x) in values.iter().enumerate() {
                v += self.derivative[row * self.size + c] * x;
            }
            *o = inverse_length * v;
        }
        Ok(())
    }
}
