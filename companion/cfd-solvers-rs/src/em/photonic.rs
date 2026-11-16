//! Scalar photonic band (Bloch-periodic) and waveguide (Dirichlet) reference modes on small
//! rectangular grids with a dense complex Jacobi eigensolver and optional dispersive
//! permittivity, ported from upstream `cfd/em/photonic_modes` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use crate::{require, Complex64 as C, Result, SolverError};
use std::f64::consts::PI;

const SPEED_OF_LIGHT: f64 = 299792458.0;
const ZERO: C = C::new(0.0, 0.0);

/// `(cell index, trial frequency Hz) -> relative permittivity`.
pub type DispersivePermittivity = Box<dyn Fn(usize, f64) -> f64>;

/// Rectangular grid of at most 256 cells; `relative_permittivity` is ignored when a dispersive
/// callback is present.
#[derive(Default)]
pub struct PhotonicGrid2D {
    pub nx: usize,
    pub ny: usize,
    pub width_m: f64,
    pub height_m: f64,
    pub relative_permittivity: Vec<f64>,
    pub dispersive_permittivity: Option<DispersivePermittivity>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BlochWaveVector2D {
    pub kx_rad_per_m: f64,
    pub ky_rad_per_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonicModeConfig {
    pub mode_count: usize,
    pub max_jacobi_sweeps: usize,
    pub max_dispersion_iterations: usize,
    pub relative_tolerance: f64,
    pub reference_frequency_hz: f64,
}

impl Default for PhotonicModeConfig {
    fn default() -> Self {
        Self {
            mode_count: 1,
            max_jacobi_sweeps: 80,
            max_dispersion_iterations: 20,
            relative_tolerance: 1e-10,
            reference_frequency_hz: 2e14,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PhotonicBandMode {
    pub eigenvalue_per_m2: f64,
    pub frequency_hz: f64,
    pub field: Vec<C>,
    pub relative_residual: f64,
    pub iterations: usize,
    pub converged: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonicWaveguideConfig {
    pub frequency_hz: f64,
    pub mode_count: usize,
    pub max_jacobi_sweeps: usize,
    pub relative_tolerance: f64,
}

impl Default for PhotonicWaveguideConfig {
    fn default() -> Self {
        Self {
            frequency_hz: 0.0,
            mode_count: 1,
            max_jacobi_sweeps: 80,
            relative_tolerance: 1e-10,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PhotonicWaveguideMode {
    pub transverse_eigenvalue_per_m2: f64,
    pub propagation_constant_rad_per_m: f64,
    pub effective_index: f64,
    pub field: Vec<C>,
    pub relative_residual: f64,
    pub guided: bool,
    pub converged: bool,
}

struct Eigen {
    n: usize,
    matrix: Vec<C>,
    vectors: Vec<C>,
    values: Vec<f64>,
    converged: bool,
}

fn grid_size(g: &PhotonicGrid2D) -> Result<usize> {
    require(
        g.nx >= 2 && g.ny >= 2 && g.width_m > 0.0 && g.height_m > 0.0,
        "photonic grid dimensions must be positive and at least 2x2",
    )?;
    let count = g
        .nx
        .checked_mul(g.ny)
        .filter(|c| *c <= 256)
        .ok_or_else(|| SolverError::InvalidInput("photonic reference grid is too large".into()))?;
    require(
        g.relative_permittivity.len() == count || g.dispersive_permittivity.is_some(),
        "photonic permittivity profile size mismatch",
    )?;
    Ok(count)
}

fn material_profile(g: &PhotonicGrid2D, frequency_hz: f64) -> Result<Vec<f64>> {
    (0..g.nx * g.ny)
        .map(|i| {
            let v = match &g.dispersive_permittivity {
                Some(f) => f(i, frequency_hz),
                None => g.relative_permittivity[i],
            };
            require(
                v > 0.0 && v.is_finite(),
                "photonic relative permittivity must be finite and positive",
            )?;
            Ok(v)
        })
        .collect()
}

fn periodic_laplacian(g: &PhotonicGrid2D, bloch: BlochWaveVector2D) -> Vec<C> {
    let count = g.nx * g.ny;
    let (dx, dy) = (g.width_m / g.nx as f64, g.height_m / g.ny as f64);
    let (ix2, iy2) = (1.0 / (dx * dx), 1.0 / (dy * dy));
    let px = C::new(0.0, bloch.kx_rad_per_m * g.width_m).exp();
    let py = C::new(0.0, bloch.ky_rad_per_m * g.height_m).exp();
    let one = C::new(1.0, 0.0);
    let mut m = vec![ZERO; count * count];
    let idx = |i: usize, j: usize| j * g.nx + i;
    for j in 0..g.ny {
        for i in 0..g.nx {
            let row = idx(i, j);
            m[row * count + row] += 2.0 * (ix2 + iy2);
            let left = if i == 0 { g.nx - 1 } else { i - 1 };
            let right = if i + 1 == g.nx { 0 } else { i + 1 };
            let down = if j == 0 { g.ny - 1 } else { j - 1 };
            let up = if j + 1 == g.ny { 0 } else { j + 1 };
            m[row * count + idx(left, j)] += (if i == 0 { px.conj() } else { one }) * -ix2;
            m[row * count + idx(right, j)] += (if i + 1 == g.nx { px } else { one }) * -ix2;
            m[row * count + idx(i, down)] += (if j == 0 { py.conj() } else { one }) * -iy2;
            m[row * count + idx(i, up)] += (if j + 1 == g.ny { py } else { one }) * -iy2;
        }
    }
    m
}

fn dirichlet_laplacian(nx: usize, ny: usize, width_m: f64, height_m: f64) -> Vec<C> {
    let count = nx * ny;
    let (dx, dy) = (width_m / (nx + 1) as f64, height_m / (ny + 1) as f64);
    let (ix2, iy2) = (1.0 / (dx * dx), 1.0 / (dy * dy));
    let mut m = vec![ZERO; count * count];
    let idx = |i: usize, j: usize| j * nx + i;
    for j in 0..ny {
        for i in 0..nx {
            let row = idx(i, j);
            m[row * count + row] += 2.0 * (ix2 + iy2);
            if i > 0 {
                m[row * count + idx(i - 1, j)] += -ix2;
            }
            if i + 1 < nx {
                m[row * count + idx(i + 1, j)] += -ix2;
            }
            if j > 0 {
                m[row * count + idx(i, j - 1)] += -iy2;
            }
            if j + 1 < ny {
                m[row * count + idx(i, j + 1)] += -iy2;
            }
        }
    }
    m
}

/// Cyclic complex Jacobi rotations on a Hermitian matrix.
fn diagonalize(matrix: Vec<C>, n: usize, max_sweeps: usize, tolerance: f64) -> Eigen {
    let mut e = Eigen {
        n,
        matrix,
        vectors: vec![ZERO; n * n],
        values: Vec::new(),
        converged: false,
    };
    for i in 0..n {
        e.vectors[i * n + i] = C::new(1.0, 0.0);
    }
    for _ in 0..max_sweeps {
        let mut largest = 0.0_f64;
        for p in 0..n {
            for q in p + 1..n {
                let apq = e.matrix[p * n + q];
                let mag = apq.norm();
                largest = largest.max(mag);
                if mag <= tolerance {
                    continue;
                }
                let (app, aqq) = (e.matrix[p * n + p].re, e.matrix[q * n + q].re);
                let tau = (aqq - app) / (2.0 * mag);
                let sign = if tau >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (tau.abs() + (1.0 + tau * tau).sqrt());
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = apq * (c * t) / mag;
                for k in 0..n {
                    if k == p || k == q {
                        continue;
                    }
                    let (hkp, hkq) = (e.matrix[k * n + p], e.matrix[k * n + q]);
                    let nkp = hkp * c - s.conj() * hkq;
                    let nkq = s * hkp + hkq * c;
                    e.matrix[k * n + p] = nkp;
                    e.matrix[k * n + q] = nkq;
                    e.matrix[p * n + k] = nkp.conj();
                    e.matrix[q * n + k] = nkq.conj();
                }
                let ss = s.norm_sqr();
                let cross = (s.conj() * apq).re;
                e.matrix[p * n + p] = C::new(c * c * app + ss * aqq - 2.0 * c * cross, 0.0);
                e.matrix[q * n + q] = C::new(ss * app + c * c * aqq + 2.0 * c * cross, 0.0);
                e.matrix[p * n + q] = ZERO;
                e.matrix[q * n + p] = ZERO;
                for k in 0..n {
                    let (vkp, vkq) = (e.vectors[k * n + p], e.vectors[k * n + q]);
                    e.vectors[k * n + p] = vkp * c - s.conj() * vkq;
                    e.vectors[k * n + q] = s * vkp + vkq * c;
                }
            }
        }
        if largest <= tolerance {
            e.converged = true;
            break;
        }
    }
    e.values = (0..n).map(|i| e.matrix[i * n + i].re).collect();
    e
}

fn sorted_indices(values: &[f64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|a, b| values[*a].total_cmp(&values[*b]).then(a.cmp(b)));
    order
}

fn column(e: &Eigen, index: usize) -> Vec<C> {
    (0..e.n).map(|row| e.vectors[row * e.n + index]).collect()
}

fn vector_norm(v: &[C]) -> f64 {
    v.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt()
}

fn band_residual(field: &[C], laplacian: &[C], epsilon: &[f64], eigenvalue: f64) -> f64 {
    let n = field.len();
    let mut residual = vec![ZERO; n];
    let mut lhs = vec![ZERO; n];
    let mut matrix_scale = 0.0_f64;
    for row in 0..n {
        let mut value = ZERO;
        let mut row_scale = 0.0;
        for col in 0..n {
            value += laplacian[row * n + col] * field[col];
            row_scale += laplacian[row * n + col].norm();
        }
        lhs[row] = value;
        residual[row] = value - field[row] * (eigenvalue * epsilon[row]);
        matrix_scale = matrix_scale.max(row_scale);
    }
    let mass_norm2: f64 = (0..n).map(|i| (field[i] * epsilon[i]).norm_sqr()).sum();
    let scale = (vector_norm(&lhs) + eigenvalue.abs() * mass_norm2.sqrt())
        .max(matrix_scale * vector_norm(field))
        .max(f64::MIN_POSITIVE);
    vector_norm(&residual) / scale
}

fn operator_residual(field: &[C], matrix: &[C], eigenvalue: f64) -> f64 {
    let n = field.len();
    let mut residual = vec![ZERO; n];
    for row in 0..n {
        for col in 0..n {
            residual[row] += matrix[row * n + col] * field[col];
        }
        residual[row] -= field[row] * eigenvalue;
    }
    vector_norm(&residual) / (eigenvalue.abs().max(1.0) * vector_norm(field))
}

/// Bloch band modes of `-lap u = (omega/c)^2 eps u` (Hermitian FD Laplacian with Bloch phases on
/// periodic edges), with fixed-point iteration on frequency for dispersive materials.
pub fn solve_photonic_band_modes(
    grid: &PhotonicGrid2D,
    bloch: BlochWaveVector2D,
    config: &PhotonicModeConfig,
) -> Result<Vec<PhotonicBandMode>> {
    let count = grid_size(grid)?;
    require(
        config.mode_count > 0
            && config.mode_count <= count
            && config.max_jacobi_sweeps > 0
            && config.max_dispersion_iterations > 0
            && config.relative_tolerance > 0.0,
        "invalid photonic band solver controls",
    )?;
    require(
        bloch.kx_rad_per_m.is_finite() && bloch.ky_rad_per_m.is_finite(),
        "Bloch wave vector must be finite",
    )?;
    let dispersive = grid.dispersive_permittivity.is_some();
    require(
        !dispersive || config.reference_frequency_hz > 0.0,
        "dispersive photonic solve needs a positive reference frequency",
    )?;
    let laplacian = periodic_laplacian(grid, bloch);
    let mut out = Vec::with_capacity(config.mode_count);
    for mode_index in 0..config.mode_count {
        let mut frequency = if dispersive {
            config.reference_frequency_hz
        } else {
            0.0
        };
        let mut iterations = 0;
        let mut solved = None;
        for iteration in 0..config.max_dispersion_iterations {
            let epsilon = material_profile(grid, frequency)?;
            let scaled: Vec<C> = (0..count * count)
                .map(|k| laplacian[k] / (epsilon[k / count] * epsilon[k % count]).sqrt())
                .collect();
            let eigen = diagonalize(
                scaled,
                count,
                config.max_jacobi_sweeps,
                config.relative_tolerance,
            );
            let order = sorted_indices(&eigen.values);
            let value = eigen.values[order[mode_index]].max(0.0);
            let updated = SPEED_OF_LIGHT * value.sqrt() / (2.0 * PI);
            iterations = iteration + 1;
            let done = !dispersive
                || (frequency > 0.0
                    && (updated - frequency).abs() / frequency <= config.relative_tolerance);
            solved = Some((eigen, epsilon));
            if done {
                break;
            }
            frequency = if updated > 0.0 {
                updated
            } else {
                config.reference_frequency_hz
            };
        }
        let Some((eigen, epsilon)) = solved else {
            return Err(SolverError::NumericalFailure(
                "photonic band solve did not run".into(),
            ));
        };
        let order = sorted_indices(&eigen.values);
        let selected = order[mode_index];
        let eigenvalue = eigen.values[selected].max(0.0);
        let mut field = column(&eigen, selected);
        let norm2: f64 = field
            .iter()
            .zip(&epsilon)
            .map(|(f, e)| e * f.norm_sqr())
            .sum();
        require(
            norm2 > 0.0 && norm2.is_finite(),
            "invalid photonic mode norm",
        )?;
        let scale = 1.0 / norm2.sqrt();
        field.iter_mut().for_each(|v| *v *= scale);
        let relative_residual = band_residual(&field, &laplacian, &epsilon, eigenvalue);
        out.push(PhotonicBandMode {
            eigenvalue_per_m2: eigenvalue,
            frequency_hz: SPEED_OF_LIGHT * eigenvalue.sqrt() / (2.0 * PI),
            field,
            relative_residual,
            iterations,
            converged: eigen.converged,
        });
    }
    Ok(out)
}

/// Scalar guided modes from negative eigenvalues of `-lap - k0^2 eps` with Dirichlet walls.
pub fn solve_photonic_waveguide_modes(
    grid: &PhotonicGrid2D,
    config: &PhotonicWaveguideConfig,
) -> Result<Vec<PhotonicWaveguideMode>> {
    let count = grid_size(grid)?;
    require(
        config.frequency_hz > 0.0
            && config.mode_count > 0
            && config.mode_count <= count
            && config.max_jacobi_sweeps > 0
            && config.relative_tolerance > 0.0,
        "invalid photonic waveguide controls",
    )?;
    let epsilon = material_profile(grid, config.frequency_hz)?;
    let mut op = dirichlet_laplacian(grid.nx, grid.ny, grid.width_m, grid.height_m);
    let k0 = 2.0 * PI * config.frequency_hz / SPEED_OF_LIGHT;
    for i in 0..count {
        op[i * count + i] -= k0 * k0 * epsilon[i];
    }
    let original = op.clone();
    let eigen = diagonalize(
        op,
        count,
        config.max_jacobi_sweeps,
        config.relative_tolerance,
    );
    let order = sorted_indices(&eigen.values);
    let mut out = Vec::with_capacity(config.mode_count);
    for &selected in order.iter().take(config.mode_count) {
        let eigenvalue = eigen.values[selected];
        let guided = eigenvalue < 0.0;
        let beta = if guided { (-eigenvalue).sqrt() } else { 0.0 };
        let mut field = column(&eigen, selected);
        let norm = vector_norm(&field);
        require(norm > 0.0, "empty photonic waveguide eigenfield")?;
        field.iter_mut().for_each(|v| *v /= norm);
        let relative_residual = operator_residual(&field, &original, eigenvalue);
        out.push(PhotonicWaveguideMode {
            transverse_eigenvalue_per_m2: eigenvalue,
            propagation_constant_rad_per_m: beta,
            effective_index: beta / k0,
            field,
            relative_residual,
            guided,
            converged: eigen.converged,
        });
    }
    Ok(out)
}
