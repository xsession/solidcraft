use crate::{require, Result};
use std::cell::RefCell;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::rc::Rc;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dual<const N: usize> {
    pub value: f64,
    pub derivative: [f64; N],
}
impl<const N: usize> Dual<N> {
    pub fn constant(value: f64) -> Self {
        Self {
            value,
            derivative: [0.0; N],
        }
    }
    pub fn variable(index: usize, value: f64) -> Result<Self> {
        require(index < N && value.is_finite(), "dual variable index/value")?;
        let mut x = Self::constant(value);
        x.derivative[index] = 1.0;
        Ok(x)
    }
    fn unary(self, value: f64, scale: f64) -> Self {
        Self {
            value,
            derivative: self.derivative.map(|d| d * scale),
        }
    }
    pub fn sin(self) -> Self {
        self.unary(self.value.sin(), self.value.cos())
    }
    pub fn cos(self) -> Self {
        self.unary(self.value.cos(), -self.value.sin())
    }
    pub fn exp(self) -> Self {
        let v = self.value.exp();
        self.unary(v, v)
    }
    pub fn ln(self) -> Result<Self> {
        require(self.value > 0.0, "dual logarithm requires a positive value")?;
        Ok(self.unary(self.value.ln(), 1.0 / self.value))
    }
    pub fn sqrt(self) -> Result<Self> {
        require(
            self.value >= 0.0,
            "dual square root requires a nonnegative value",
        )?;
        let v = self.value.sqrt();
        Ok(self.unary(v, if v > 0.0 { 0.5 / v } else { 0.0 }))
    }
    pub fn powf(self, e: f64) -> Result<Self> {
        require(
            e.is_finite() && !(self.value == 0.0 && e < 1.0),
            "invalid dual power",
        )?;
        let v = self.value.powf(e);
        let d = if e == 0.0 {
            0.0
        } else {
            e * self.value.powf(e - 1.0)
        };
        require(v.is_finite() && d.is_finite(), "nonfinite dual power")?;
        Ok(self.unary(v, d))
    }
    pub fn checked_div(self, b: Self) -> Result<Self> {
        require(b.value != 0.0, "dual division by zero")?;
        let s = 1.0 / (b.value * b.value);
        Ok(Self {
            value: self.value / b.value,
            derivative: std::array::from_fn(|i| {
                (self.derivative[i] * b.value - self.value * b.derivative[i]) * s
            }),
        })
    }
}
impl<const N: usize> From<f64> for Dual<N> {
    fn from(v: f64) -> Self {
        Self::constant(v)
    }
}
impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            value: self.value + b.value,
            derivative: std::array::from_fn(|i| self.derivative[i] + b.derivative[i]),
        }
    }
}
impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            value: self.value - b.value,
            derivative: std::array::from_fn(|i| self.derivative[i] - b.derivative[i]),
        }
    }
}
impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            value: self.value * b.value,
            derivative: std::array::from_fn(|i| {
                self.derivative[i] * b.value + self.value * b.derivative[i]
            }),
        }
    }
}
impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        // IEEE semantics (inf/NaN) like `f64`; use `checked_div` for an error instead.
        let s = 1.0 / (b.value * b.value);
        Self {
            value: self.value / b.value,
            derivative: std::array::from_fn(|i| {
                (self.derivative[i] * b.value - self.value * b.derivative[i]) * s
            }),
        }
    }
}
impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    fn neg(self) -> Self {
        self.unary(-self.value, -1.0)
    }
}
macro_rules! scalar_op {
    ($trait:ident,$method:ident) => {
        impl<const N: usize> $trait<f64> for Dual<N> {
            type Output = Self;
            fn $method(self, b: f64) -> Self {
                self.$method(Self::constant(b))
            }
        }
        impl<const N: usize> $trait<Dual<N>> for f64 {
            type Output = Dual<N>;
            fn $method(self, b: Dual<N>) -> Dual<N> {
                Dual::constant(self).$method(b)
            }
        }
    };
}
scalar_op!(Add, add);
scalar_op!(Sub, sub);
scalar_op!(Mul, mul);
scalar_op!(Div, div);

#[derive(Clone, Debug)]
struct ReverseNode {
    value: f64,
    parents: Vec<(usize, f64)>,
}
#[derive(Clone, Debug, Default)]
pub struct ReverseTape {
    nodes: Rc<RefCell<Vec<ReverseNode>>>,
}
#[derive(Clone, Debug)]
pub struct ReverseScalar {
    index: usize,
    tape: ReverseTape,
}
impl ReverseTape {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn variable(&self, value: f64) -> ReverseScalar {
        let mut nodes = self.nodes.borrow_mut();
        let index = nodes.len();
        nodes.push(ReverseNode {
            value,
            parents: vec![],
        });
        ReverseScalar {
            index,
            tape: self.clone(),
        }
    }
    pub fn backward(&self, output: &ReverseScalar, inputs: &[ReverseScalar]) -> Result<Vec<f64>> {
        require(
            Rc::ptr_eq(&self.nodes, &output.tape.nodes)
                && inputs
                    .iter()
                    .all(|v| Rc::ptr_eq(&self.nodes, &v.tape.nodes)),
            "reverse AD tapes differ",
        )?;
        let nodes = self.nodes.borrow();
        let mut gradient = vec![0.0; nodes.len()];
        gradient[output.index] = 1.0;
        for i in (0..=output.index).rev() {
            for &(parent, partial) in &nodes[i].parents {
                gradient[parent] += gradient[i] * partial;
            }
        }
        Ok(inputs.iter().map(|v| gradient[v.index]).collect())
    }
    pub fn len(&self) -> usize {
        self.nodes.borrow().len()
    }
    pub fn is_empty(&self) -> bool {
        self.nodes.borrow().is_empty()
    }
}
impl ReverseScalar {
    pub fn value(&self) -> f64 {
        self.tape.nodes.borrow()[self.index].value
    }
    fn unary(&self, value: f64, scale: f64) -> Self {
        let mut nodes = self.tape.nodes.borrow_mut();
        let index = nodes.len();
        nodes.push(ReverseNode {
            value,
            parents: vec![(self.index, scale)],
        });
        Self {
            index,
            tape: self.tape.clone(),
        }
    }
    fn binary(&self, b: &Self, value: f64, left: f64, right: f64) -> Self {
        assert!(
            Rc::ptr_eq(&self.tape.nodes, &b.tape.nodes),
            "reverse AD tapes differ"
        );
        let mut nodes = self.tape.nodes.borrow_mut();
        let index = nodes.len();
        nodes.push(ReverseNode {
            value,
            parents: vec![(self.index, left), (b.index, right)],
        });
        Self {
            index,
            tape: self.tape.clone(),
        }
    }
    pub fn sin(&self) -> Self {
        self.unary(self.value().sin(), self.value().cos())
    }
    pub fn cos(&self) -> Self {
        self.unary(self.value().cos(), -self.value().sin())
    }
    pub fn exp(&self) -> Self {
        let v = self.value().exp();
        self.unary(v, v)
    }
    pub fn ln(&self) -> Result<Self> {
        let v = self.value();
        require(v > 0.0, "reverse AD logarithm requires a positive value")?;
        Ok(self.unary(v.ln(), 1.0 / v))
    }
    pub fn sqrt(&self) -> Result<Self> {
        let v = self.value();
        require(
            v >= 0.0,
            "reverse AD square root requires a nonnegative value",
        )?;
        let root = v.sqrt();
        Ok(self.unary(root, if root > 0.0 { 0.5 / root } else { 0.0 }))
    }
    pub fn powf(&self, e: f64) -> Result<Self> {
        let v = self.value();
        require(
            e.is_finite() && !(v == 0.0 && e < 1.0),
            "invalid reverse AD power",
        )?;
        let value = v.powf(e);
        let partial = if e == 0.0 { 0.0 } else { e * v.powf(e - 1.0) };
        require(
            value.is_finite() && partial.is_finite(),
            "nonfinite reverse AD power",
        )?;
        Ok(self.unary(value, partial))
    }
    pub fn checked_div(&self, b: &Self) -> Result<Self> {
        let value = b.value();
        require(value != 0.0, "reverse AD division by zero")?;
        Ok(self.binary(
            b,
            self.value() / value,
            1.0 / value,
            -self.value() / (value * value),
        ))
    }
}
impl Add for ReverseScalar {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        self.binary(&b, self.value() + b.value(), 1.0, 1.0)
    }
}
impl Sub for ReverseScalar {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self.binary(&b, self.value() - b.value(), 1.0, -1.0)
    }
}
impl Mul for ReverseScalar {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        self.binary(&b, self.value() * b.value(), b.value(), self.value())
    }
}
impl Div for ReverseScalar {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        // IEEE semantics (inf/NaN) like `f64`; use `checked_div` for an error instead.
        let value = b.value();
        self.binary(
            &b,
            self.value() / value,
            1.0 / value,
            -self.value() / (value * value),
        )
    }
}
impl Neg for ReverseScalar {
    type Output = Self;
    fn neg(self) -> Self {
        self.unary(-self.value(), -1.0)
    }
}
macro_rules! reverse_scalar_op {
    ($trait:ident,$method:ident) => {
        impl $trait<f64> for ReverseScalar {
            type Output = Self;
            fn $method(self, b: f64) -> Self {
                let constant = self.tape.variable(b);
                self.$method(constant)
            }
        }
        impl $trait<ReverseScalar> for f64 {
            type Output = ReverseScalar;
            fn $method(self, b: ReverseScalar) -> ReverseScalar {
                b.tape.variable(self).$method(b)
            }
        }
    };
}
reverse_scalar_op!(Add, add);
reverse_scalar_op!(Sub, sub);
reverse_scalar_op!(Mul, mul);
reverse_scalar_op!(Div, div);
