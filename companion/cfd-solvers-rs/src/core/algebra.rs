use crate::{finite_slice, positive, require, Complex64, Result, SolverError};
use serde::{Deserialize, Serialize};
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Vec3 {
    pub const ZERO: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn dot(self, b: Self) -> f64 {
        self.x * b.x + self.y * b.y + self.z * b.z
    }
    pub fn cross(self, b: Self) -> Self {
        Self::new(
            self.y * b.z - self.z * b.y,
            self.z * b.x - self.x * b.z,
            self.x * b.y - self.y * b.x,
        )
    }
    pub fn norm(self) -> f64 {
        self.x.hypot(self.y).hypot(self.z)
    }
    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }
    pub fn normalized(self) -> Result<Self> {
        let n = self.norm();
        positive(n, "vector length")?;
        Ok(self / n)
    }
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
    pub fn component(self, i: usize) -> f64 {
        match i {
            0 => self.x,
            1 => self.y,
            2 => self.z,
            // Out-of-range axes yield NaN, which every solver input check rejects.
            _ => f64::NAN,
        }
    }
    pub fn as_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}
impl From<[f64; 3]> for Vec3 {
    fn from(p: [f64; 3]) -> Self {
        Self::new(p[0], p[1], p[2])
    }
}
impl Add for Vec3 {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self::new(self.x + b.x, self.y + b.y, self.z + b.z)
    }
}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self::new(self.x - b.x, self.y - b.y, self.z - b.z)
    }
}
impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        self * (-1.0)
    }
}
impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, b: f64) -> Self {
        Self::new(self.x * b, self.y * b, self.z * b)
    }
}
impl Mul<Vec3> for f64 {
    type Output = Vec3;
    fn mul(self, b: Vec3) -> Vec3 {
        b * self
    }
}
impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, b: f64) -> Self {
        self * (1.0 / b)
    }
}
impl AddAssign for Vec3 {
    fn add_assign(&mut self, b: Self) {
        *self = *self + b;
    }
}
impl SubAssign for Vec3 {
    fn sub_assign(&mut self, b: Self) {
        *self = *self - b;
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Quaternion {
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Default for Quaternion {
    fn default() -> Self {
        Self {
            w: 1.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}
impl Quaternion {
    pub fn norm(self) -> f64 {
        self.w.hypot(self.x).hypot(self.y).hypot(self.z)
    }
    pub fn normalized(self) -> Result<Self> {
        let n = self.norm();
        positive(n, "quaternion norm")?;
        Ok(Self {
            w: self.w / n,
            x: self.x / n,
            y: self.y / n,
            z: self.z / n,
        })
    }
    pub fn conjugate(self) -> Self {
        Self {
            w: self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
    pub fn axis_angle(axis: Vec3, angle: f64) -> Result<Self> {
        require(angle.is_finite(), "angle must be finite")?;
        let a = axis.normalized()?;
        let s = (angle * 0.5).sin();
        Ok(Self {
            w: (angle * 0.5).cos(),
            x: a.x * s,
            y: a.y * s,
            z: a.z * s,
        })
    }
    pub fn rotate(self, v: Vec3) -> Vec3 {
        let u = Vec3::new(self.x, self.y, self.z);
        v * (self.w * self.w - u.dot(u)) + u * (2.0 * u.dot(v)) + u.cross(v) * (2.0 * self.w)
    }
    pub fn integrate(self, omega: Vec3, dt: f64) -> Result<Self> {
        require(
            omega.is_finite() && dt.is_finite(),
            "rotation update must be finite",
        )?;
        let q = Self {
            w: 0.0,
            x: omega.x,
            y: omega.y,
            z: omega.z,
        } * self;
        Self {
            w: self.w + dt * 0.5 * q.w,
            x: self.x + dt * 0.5 * q.x,
            y: self.y + dt * 0.5 * q.y,
            z: self.z + dt * 0.5 * q.z,
        }
        .normalized()
    }
}
impl Mul for Quaternion {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            w: self.w * b.w - self.x * b.x - self.y * b.y - self.z * b.z,
            x: self.w * b.x + self.x * b.w + self.y * b.z - self.z * b.y,
            y: self.w * b.y - self.x * b.z + self.y * b.w + self.z * b.x,
            z: self.w * b.z + self.x * b.y - self.y * b.x + self.z * b.w,
        }
    }
}

pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
pub fn rms(a: &[f64]) -> f64 {
    if a.is_empty() {
        0.0
    } else {
        // Scaling prevents overflow in x*x from turning a finite residual into
        // infinity and producing a false initial convergence decision.
        let mut scale: f64 = 0.0;
        for value in a {
            if value.is_nan() {
                return f64::NAN;
            }
            scale = scale.max(value.abs());
        }
        if scale == 0.0 {
            return 0.0;
        }
        if scale.is_infinite() {
            return f64::INFINITY;
        }
        scale * (a.iter().map(|v| (v / scale).powi(2)).sum::<f64>() / a.len() as f64).sqrt()
    }
}
pub fn l2_norm(a: &[f64]) -> f64 {
    a.iter().fold(0.0, |n, v| n.hypot(*v))
}

/// Partial-pivot dense solve. Matrix is row-major; inputs are retained on error.
pub fn dense_solve(matrix: &[f64], rhs: &[f64]) -> Result<Vec<f64>> {
    let n = rhs.len();
    require(n > 0 && matrix.len() == n * n, "dense dimensions disagree")?;
    finite_slice(matrix, "matrix")?;
    finite_slice(rhs, "rhs")?;
    let mut a = matrix.to_vec();
    let mut b = rhs.to_vec();
    let scale = a.iter().map(|x| x.abs()).fold(0.0, f64::max);
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap_or(k);
        if a[p * n + k].abs() <= f64::EPSILON * scale.max(f64::MIN_POSITIVE) * n as f64 {
            return Err(SolverError::Singular(format!("zero pivot at {k}")));
        }
        if p != k {
            for j in k..n {
                a.swap(k * n + j, p * n + j);
            }
            b.swap(k, p);
        }
        for i in k + 1..n {
            let m = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.0;
            for j in k + 1..n {
                a[i * n + j] -= m * a[k * n + j];
            }
            b[i] -= m * b[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            b[i] -= a[i * n + j] * b[j];
        }
        b[i] /= a[i * n + i];
    }
    finite_slice(&b, "dense result")?;
    Ok(b)
}

pub fn complex_dense_solve(matrix: &[Complex64], rhs: &[Complex64]) -> Result<Vec<Complex64>> {
    let n = rhs.len();
    require(
        n > 0 && matrix.len() == n * n,
        "complex dense dimensions disagree",
    )?;
    require(
        matrix
            .iter()
            .chain(rhs)
            .all(|v| v.re.is_finite() && v.im.is_finite()),
        "nonfinite complex system",
    )?;
    let mut a = matrix.to_vec();
    let mut b = rhs.to_vec();
    let scale = a.iter().map(|v| v.norm()).fold(0.0, f64::max);
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| a[i * n + k].norm_sqr().total_cmp(&a[j * n + k].norm_sqr()))
            .unwrap_or(k);
        if a[p * n + k].norm() <= f64::EPSILON * scale.max(f64::MIN_POSITIVE) * n as f64 {
            return Err(SolverError::Singular("complex pivot".into()));
        }
        if p != k {
            for j in k..n {
                a.swap(k * n + j, p * n + j);
            }
            b.swap(k, p);
        }
        for i in k + 1..n {
            let m = a[i * n + k] / a[k * n + k];
            for j in k + 1..n {
                let pivot = a[k * n + j];
                a[i * n + j] -= m * pivot;
            }
            let value = b[k];
            b[i] -= m * value;
            a[i * n + k] = Complex64::new(0.0, 0.0);
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            let v = b[j];
            b[i] -= a[i * n + j] * v;
        }
        b[i] /= a[i * n + i];
    }
    require(
        b.iter().all(|v| v.re.is_finite() && v.im.is_finite()),
        "complex solve overflow",
    )?;
    Ok(b)
}

pub fn solve_tridiagonal(
    lower: &[f64],
    diagonal: &[f64],
    upper: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>> {
    let n = diagonal.len();
    require(
        n > 0 && rhs.len() == n && lower.len() == n && upper.len() == n,
        "tridiagonal dimensions",
    )?;
    finite_slice(lower, "lower")?;
    finite_slice(diagonal, "diagonal")?;
    finite_slice(upper, "upper")?;
    finite_slice(rhs, "rhs")?;
    let mut d = diagonal.to_vec();
    let mut b = rhs.to_vec();
    for i in 1..n {
        if d[i - 1].abs() < 1e-30 {
            return Err(SolverError::Singular("tridiagonal pivot".into()));
        }
        let m = lower[i] / d[i - 1];
        d[i] -= m * upper[i - 1];
        b[i] -= m * b[i - 1];
    }
    if d[n - 1].abs() < 1e-30 {
        return Err(SolverError::Singular("tridiagonal final pivot".into()));
    }
    b[n - 1] /= d[n - 1];
    for i in (0..n - 1).rev() {
        if d[i].abs() < 1e-30 {
            return Err(SolverError::Singular("tridiagonal pivot".into()));
        }
        b[i] = (b[i] - upper[i] * b[i + 1]) / d[i];
    }
    finite_slice(&b, "tridiagonal result")?;
    Ok(b)
}

/// In-place radix-2 FFT, same forward sign and inverse normalization as upstream.
pub fn fft(values: &mut [Complex64], inverse: bool) -> Result<()> {
    let n = values.len();
    require(
        n > 0 && n.is_power_of_two(),
        "FFT length must be a positive power of two",
    )?;
    require(
        values.iter().all(|v| v.re.is_finite() && v.im.is_finite()),
        "nonfinite FFT input",
    )?;
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            values.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let angle = if inverse { 2.0 } else { -2.0 } * std::f64::consts::PI / len as f64;
        let root = Complex64::from_polar(1.0, angle);
        for start in (0..n).step_by(len) {
            let mut w = Complex64::new(1.0, 0.0);
            for k in 0..len / 2 {
                let u = values[start + k];
                let v = values[start + k + len / 2] * w;
                values[start + k] = u + v;
                values[start + k + len / 2] = u - v;
                w *= root;
            }
        }
        len = match len.checked_mul(2) {
            Some(v) => v,
            None => break,
        };
    }
    if inverse {
        for v in values.iter_mut() {
            *v /= n as f64;
        }
    }
    Ok(())
}
pub fn fft2d(values: &mut [Complex64], nx: usize, ny: usize, inverse: bool) -> Result<()> {
    require(
        nx.checked_mul(ny) == Some(values.len()) && nx.is_power_of_two() && ny.is_power_of_two(),
        "FFT2D dimensions",
    )?;
    for row in values.chunks_exact_mut(nx) {
        fft(row, inverse)?;
    }
    let mut col = vec![Complex64::default(); ny];
    for x in 0..nx {
        for y in 0..ny {
            col[y] = values[y * nx + x];
        }
        fft(&mut col, inverse)?;
        for y in 0..ny {
            values[y * nx + x] = col[y];
        }
    }
    Ok(())
}

pub fn wavenumbers(n: usize, spacing: f64) -> Result<Vec<f64>> {
    require(
        n.is_power_of_two(),
        "spectral grid must be a nonzero power of two",
    )?;
    positive(spacing, "spectral spacing")?;
    let factor = 2.0 * std::f64::consts::PI / (n as f64 * spacing);
    Ok((0..n)
        .map(|i| {
            factor
                * if i <= n / 2 {
                    i as f64
                } else {
                    i as f64 - n as f64
                }
        })
        .collect())
}
pub fn derivative_periodic_1d(field: &[f64], spacing: f64) -> Result<Vec<f64>> {
    let k = wavenumbers(field.len(), spacing)?;
    finite_slice(field, "spectral field")?;
    let mut values: Vec<_> = field.iter().map(|&v| Complex64::new(v, 0.0)).collect();
    fft(&mut values, false)?;
    let nyquist = values.len() / 2;
    for (i, v) in values.iter_mut().enumerate() {
        *v *= Complex64::new(0.0, if i == nyquist { 0.0 } else { k[i] });
    }
    fft(&mut values, true)?;
    Ok(values.iter().map(|v| v.re).collect())
}
pub fn fractional_laplacian_periodic_1d(
    field: &[f64],
    spacing: f64,
    order: f64,
) -> Result<Vec<f64>> {
    require(
        order.is_finite() && order >= 0.0,
        "invalid fractional Laplacian order",
    )?;
    let k = wavenumbers(field.len(), spacing)?;
    finite_slice(field, "spectral field")?;
    let mut values: Vec<_> = field.iter().map(|&v| Complex64::new(v, 0.0)).collect();
    fft(&mut values, false)?;
    for (i, v) in values.iter_mut().enumerate() {
        *v *= if k[i] == 0.0 {
            if order == 0.0 {
                1.0
            } else {
                0.0
            }
        } else {
            k[i].abs().powf(2.0 * order)
        };
    }
    fft(&mut values, true)?;
    Ok(values.iter().map(|v| v.re).collect())
}
