//! 3-D Neo-Hookean and small-strain J2 constitutive updates, objective finite-strain J2 and
//! total-Lagrangian Tet4 kernels.
#![allow(clippy::needless_range_loop)]
use crate::{core::Vec3, fem::Node2, require, Result};

/// Row-major 3x3 tensor.
pub type Tensor3 = [f64; 9];
/// Voigt order xx, yy, zz, yz, xz, xy with tensor (not engineering) shear strains.
pub type Voigt6 = [f64; 6];
pub type Matrix6 = [f64; 36];

pub const IDENTITY3: Tensor3 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

pub(crate) fn determinant(m: &Tensor3) -> f64 {
    m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
}

pub(crate) fn inverse(m: &Tensor3, det: f64) -> Tensor3 {
    let i = 1.0 / det;
    [
        (m[4] * m[8] - m[5] * m[7]) * i,
        (m[2] * m[7] - m[1] * m[8]) * i,
        (m[1] * m[5] - m[2] * m[4]) * i,
        (m[5] * m[6] - m[3] * m[8]) * i,
        (m[0] * m[8] - m[2] * m[6]) * i,
        (m[3] * m[2] - m[0] * m[5]) * i,
        (m[3] * m[7] - m[4] * m[6]) * i,
        (m[6] * m[1] - m[0] * m[7]) * i,
        (m[0] * m[4] - m[1] * m[3]) * i,
    ]
}

pub(crate) fn transpose(m: &Tensor3) -> Tensor3 {
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

pub(crate) fn multiply(a: &Tensor3, b: &Tensor3) -> Tensor3 {
    let mut r = [0.0; 9];
    for row in 0..3 {
        for col in 0..3 {
            for k in 0..3 {
                r[3 * row + col] += a[3 * row + k] * b[3 * k + col];
            }
        }
    }
    r
}

fn trace(m: &Tensor3) -> f64 {
    m[0] + m[4] + m[8]
}

fn tensor_inner(a: &Voigt6, b: &Voigt6) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5])
}

fn deviator(s: &Voigt6) -> Voigt6 {
    let mean = (s[0] + s[1] + s[2]) / 3.0;
    [s[0] - mean, s[1] - mean, s[2] - mean, s[3], s[4], s[5]]
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NeoHookean3DConfig {
    pub shear_modulus: f64,
    pub lame_lambda: f64,
}

impl Default for NeoHookean3DConfig {
    fn default() -> Self {
        Self {
            shear_modulus: 1.0,
            lame_lambda: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NeoHookean3DResult {
    pub determinant: f64,
    pub energy_density: f64,
    pub first_piola: Tensor3,
}

/// `W = mu/2 (I1 - 3 - 2 ln J) + lambda/2 (ln J)^2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompressibleNeoHookean3D {
    config: NeoHookean3DConfig,
}

impl CompressibleNeoHookean3D {
    pub fn new(config: NeoHookean3DConfig) -> Result<Self> {
        require(
            config.shear_modulus > 0.0
                && config.lame_lambda >= 0.0
                && config.shear_modulus.is_finite()
                && config.lame_lambda.is_finite(),
            "invalid 3-D Neo-Hookean material",
        )?;
        Ok(Self { config })
    }
    pub fn evaluate(&self, f: &Tensor3) -> Result<NeoHookean3DResult> {
        let det = determinant(f);
        require(
            det > 0.0 && det.is_finite(),
            "Neo-Hookean deformation gradient must have positive determinant",
        )?;
        let fit = transpose(&inverse(f, det));
        let ln = det.ln();
        let i1: f64 = f.iter().map(|v| v * v).sum();
        let (mu, lambda) = (self.config.shear_modulus, self.config.lame_lambda);
        let mut p = [0.0; 9];
        for i in 0..9 {
            p[i] = mu * (f[i] - fit[i]) + lambda * ln * fit[i];
        }
        Ok(NeoHookean3DResult {
            determinant: det,
            energy_density: 0.5 * mu * (i1 - 3.0 - 2.0 * ln) + 0.5 * lambda * ln * ln,
            first_piola: p,
        })
    }
    /// Directional derivative `dP/dF : direction`.
    pub fn tangent_action(&self, f: &Tensor3, direction: &Tensor3) -> Result<Tensor3> {
        let state = self.evaluate(f)?;
        let fit = transpose(&inverse(f, state.determinant));
        let finv = transpose(&fit);
        let a = multiply(&fit, &multiply(&transpose(direction), &fit));
        let dln = trace(&multiply(&finv, direction));
        let mut r = [0.0; 9];
        for i in 0..9 {
            r[i] = self.config.shear_modulus * (direction[i] + a[i])
                + self.config.lame_lambda * (dln * fit[i] - state.determinant.ln() * a[i]);
        }
        Ok(r)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct J2Plasticity3DConfig {
    pub youngs_modulus: f64,
    pub poisson_ratio: f64,
    pub yield_stress: f64,
    pub hardening_modulus: f64,
    pub tangent_step: f64,
}

impl Default for J2Plasticity3DConfig {
    fn default() -> Self {
        Self {
            youngs_modulus: 100.0,
            poisson_ratio: 0.3,
            yield_stress: 1.0,
            hardening_modulus: 0.0,
            tangent_step: 1e-7,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct J2Plasticity3DState {
    pub plastic_strain: Voigt6,
    pub equivalent_plastic_strain: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct J2Plasticity3DResult {
    pub stress: Voigt6,
    pub algorithmic_tangent: Matrix6,
    pub state: J2Plasticity3DState,
    pub trial_equivalent_stress: f64,
    pub yield_function: f64,
    pub yielded: bool,
}

#[derive(Clone, Copy, Default)]
struct J2Integration {
    stress: Voigt6,
    deviatoric_trial: Voigt6,
    state: J2Plasticity3DState,
    trial_equivalent_stress: f64,
    yield_function: f64,
    yield_level: f64,
    yielded: bool,
}

/// Small-strain radial-return J2 with linear isotropic hardening and the consistent tangent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct J2Plasticity3D {
    config: J2Plasticity3DConfig,
}

impl J2Plasticity3D {
    pub fn new(config: J2Plasticity3DConfig) -> Result<Self> {
        let c = &config;
        require(
            c.youngs_modulus > 0.0
                && c.poisson_ratio > -1.0
                && c.poisson_ratio < 0.5
                && c.yield_stress > 0.0
                && c.hardening_modulus >= 0.0
                && c.tangent_step > 0.0,
            "invalid 3-D J2 plasticity material",
        )?;
        Ok(Self { config })
    }
    fn moduli(&self) -> (f64, f64) {
        let (e, nu) = (self.config.youngs_modulus, self.config.poisson_ratio);
        (
            e / (2.0 * (1.0 + nu)),
            e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)),
        )
    }
    pub fn update(
        &self,
        total_strain: &Voigt6,
        previous: &J2Plasticity3DState,
    ) -> J2Plasticity3DResult {
        let base = self.integrate(total_strain, previous);
        J2Plasticity3DResult {
            stress: base.stress,
            state: base.state,
            trial_equivalent_stress: base.trial_equivalent_stress,
            yield_function: base.yield_function,
            yielded: base.yielded,
            algorithmic_tangent: self.analytic_tangent(&base),
        }
    }
    fn integrate(&self, total: &Voigt6, prev: &J2Plasticity3DState) -> J2Integration {
        let (shear, lame) = self.moduli();
        let mut e = [0.0; 6];
        for i in 0..6 {
            e[i] = total[i] - prev.plastic_strain[i];
        }
        let vol = e[0] + e[1] + e[2];
        let trial = [
            2.0 * shear * e[0] + lame * vol,
            2.0 * shear * e[1] + lame * vol,
            2.0 * shear * e[2] + lame * vol,
            2.0 * shear * e[3],
            2.0 * shear * e[4],
            2.0 * shear * e[5],
        ];
        let dev = deviator(&trial);
        let q = (1.5 * tensor_inner(&dev, &dev)).max(0.0).sqrt();
        let level = self.config.yield_stress
            + self.config.hardening_modulus * prev.equivalent_plastic_strain;
        let f = q - level;
        let mut r = J2Integration {
            stress: trial,
            deviatoric_trial: dev,
            state: *prev,
            trial_equivalent_stress: q,
            yield_function: f,
            yield_level: level,
            yielded: false,
        };
        if f <= 0.0 || q <= f64::EPSILON {
            return r;
        }
        let dg = f / (3.0 * shear + self.config.hardening_modulus);
        let corr = (1.0 - 3.0 * shear * dg / q).max(0.0);
        for i in 0..6 {
            r.stress[i] = if i < 3 {
                trial[i] - (1.0 - corr) * dev[i]
            } else {
                corr * dev[i]
            };
            r.state.plastic_strain[i] += 1.5 * dg * dev[i] / q;
        }
        r.state.equivalent_plastic_strain += dg;
        r.yielded = true;
        r
    }
    fn analytic_tangent(&self, base: &J2Integration) -> Matrix6 {
        let (shear, lame) = self.moduli();
        let bulk = lame + 2.0 * shear / 3.0;
        let mut t = [0.0; 36];
        for col in 0..6 {
            let mut dir = [0.0; 6];
            dir[col] = 1.0;
            let vol = if col < 3 { 1.0 } else { 0.0 };
            let mut ddir = [0.0; 6];
            for c in 0..3 {
                ddir[c] = (if c == col { 1.0 } else { 0.0 }) - vol / 3.0;
            }
            if col >= 3 {
                ddir[col] = 1.0;
            }
            let ea = [
                2.0 * shear * dir[0] + lame * vol,
                2.0 * shear * dir[1] + lame * vol,
                2.0 * shear * dir[2] + lame * vol,
                2.0 * shear * dir[3],
                2.0 * shear * dir[4],
                2.0 * shear * dir[5],
            ];
            if !base.yielded {
                for row in 0..6 {
                    t[row * 6 + col] = ea[row];
                }
                continue;
            }
            let q = base.trial_equivalent_stress;
            let den = 3.0 * shear + self.config.hardening_modulus;
            let alpha = (1.0 - 3.0 * shear * (q - base.yield_level) / (den * q)).max(0.0);
            let coef = 3.0 * shear * base.yield_level / (den * q * q);
            let mut da = [0.0; 6];
            for row in 0..6 {
                da[row] = 2.0 * shear * ddir[row];
            }
            let dq = 1.5 * tensor_inner(&base.deviatoric_trial, &da) / q;
            for row in 0..3 {
                t[row * 6 + col] =
                    bulk * vol + alpha * da[row] - coef * dq * base.deviatoric_trial[row];
            }
            for row in 3..6 {
                t[row * 6 + col] = alpha * da[row] - coef * dq * base.deviatoric_trial[row];
            }
        }
        t
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteStrainJ2Config {
    pub youngs_modulus: f64,
    pub poisson_ratio: f64,
    pub yield_stress: f64,
    pub hardening_modulus: f64,
    pub tangent_step: f64,
    pub spectral_iterations: u32,
}

impl Default for FiniteStrainJ2Config {
    fn default() -> Self {
        Self {
            youngs_modulus: 200.0,
            poisson_ratio: 0.3,
            yield_stress: 1.0,
            hardening_modulus: 0.0,
            tangent_step: 1e-7,
            spectral_iterations: 32,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteStrainJ2State {
    pub plastic_deformation_gradient: Tensor3,
    pub equivalent_plastic_strain: f64,
}

impl Default for FiniteStrainJ2State {
    fn default() -> Self {
        Self {
            plastic_deformation_gradient: IDENTITY3,
            equivalent_plastic_strain: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FiniteStrainJ2Result {
    pub determinant: f64,
    pub energy_density: f64,
    pub trial_equivalent_stress: f64,
    pub yield_function: f64,
    pub yielded: bool,
    pub cauchy_stress: Tensor3,
    pub first_piola: Tensor3,
    pub state: FiniteStrainJ2State,
}

/// Cyclic Jacobi diagonalization of a symmetric 3x3 matrix, eigenvalues sorted descending.
pub(crate) fn diagonalize_symmetric3(input: &Tensor3, iterations: u32) -> ([f64; 3], Tensor3) {
    let mut m = *input;
    let mut v = IDENTITY3;
    for _ in 0..iterations {
        let (mut p, mut q) = (0usize, 1usize);
        let mut largest = m[1].abs();
        if m[2].abs() > largest {
            p = 0;
            q = 2;
            largest = m[2].abs();
        }
        if m[5].abs() > largest {
            p = 1;
            q = 2;
            largest = m[5].abs();
        }
        if largest.is_nan() || largest <= 1e-14 {
            break;
        }
        let (app, aqq, apq) = (m[3 * p + p], m[3 * q + q], m[3 * p + q]);
        let tau = (aqq - app) / (2.0 * apq);
        let sign = if tau >= 0.0 { 1.0 } else { -1.0 };
        let t = sign / (tau.abs() + (1.0 + tau * tau).sqrt());
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        for k in 0..3 {
            if k == p || k == q {
                continue;
            }
            let (akp, akq) = (m[3 * k + p], m[3 * k + q]);
            m[3 * k + p] = c * akp - s * akq;
            m[3 * p + k] = m[3 * k + p];
            m[3 * k + q] = s * akp + c * akq;
            m[3 * q + k] = m[3 * k + q];
        }
        m[3 * p + p] = app - t * apq;
        m[3 * q + q] = aqq + t * apq;
        m[3 * p + q] = 0.0;
        m[3 * q + p] = 0.0;
        for k in 0..3 {
            let (vkp, vkq) = (v[3 * k + p], v[3 * k + q]);
            v[3 * k + p] = c * vkp - s * vkq;
            v[3 * k + q] = s * vkp + c * vkq;
        }
    }
    let mut e = [m[0], m[4], m[8]];
    for a in 0..3 {
        for b in a + 1..3 {
            if e[b] > e[a] {
                e.swap(a, b);
                for row in 0..3 {
                    v.swap(3 * row + a, 3 * row + b);
                }
            }
        }
    }
    (e, v)
}

pub(crate) fn reconstruct(v: &Tensor3, d: &[f64; 3]) -> Tensor3 {
    let mut r = [0.0; 9];
    for row in 0..3 {
        for col in 0..3 {
            for mode in 0..3 {
                r[3 * row + col] += v[3 * row + mode] * d[mode] * v[3 * col + mode];
            }
        }
    }
    r
}

/// Objective multiplicative (F = Fe Fp) finite-strain J2 integrated in principal logarithmic
/// elastic strain after a polar split of the trial elastic gradient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteStrainJ2 {
    config: FiniteStrainJ2Config,
}

impl FiniteStrainJ2 {
    pub fn new(config: FiniteStrainJ2Config) -> Result<Self> {
        let c = &config;
        require(
            c.youngs_modulus > 0.0
                && c.poisson_ratio > -1.0
                && c.poisson_ratio < 0.5
                && c.yield_stress > 0.0
                && c.hardening_modulus >= 0.0
                && c.tangent_step > 0.0
                && c.spectral_iterations > 0
                && c.youngs_modulus.is_finite()
                && c.poisson_ratio.is_finite()
                && c.yield_stress.is_finite()
                && c.hardening_modulus.is_finite()
                && c.tangent_step.is_finite(),
            "invalid finite-strain J2 configuration",
        )?;
        Ok(Self { config })
    }
    pub fn config(&self) -> &FiniteStrainJ2Config {
        &self.config
    }
    pub fn evaluate(
        &self,
        f: &Tensor3,
        previous: &FiniteStrainJ2State,
    ) -> Result<FiniteStrainJ2Result> {
        self.integrate(f, previous)
    }
    /// One-sided finite-difference tangent action with the configured step.
    pub fn tangent_action(
        &self,
        f: &Tensor3,
        direction: &Tensor3,
        previous: &FiniteStrainJ2State,
    ) -> Result<Tensor3> {
        self.tangent_action_with_step(f, direction, previous, self.config.tangent_step, false)
    }
    pub fn central_tangent_action(
        &self,
        f: &Tensor3,
        direction: &Tensor3,
        previous: &FiniteStrainJ2State,
    ) -> Result<Tensor3> {
        self.tangent_action_with_step(f, direction, previous, self.config.tangent_step, true)
    }
    pub fn tangent_action_with_step(
        &self,
        f: &Tensor3,
        direction: &Tensor3,
        previous: &FiniteStrainJ2State,
        step_scale: f64,
        central_difference: bool,
    ) -> Result<Tensor3> {
        require(
            direction.iter().all(|v| v.is_finite()),
            "finite-strain J2 direction must be finite",
        )?;
        require(
            step_scale > 0.0 && step_scale.is_finite(),
            "finite-strain J2 tangent step must be positive and finite",
        )?;
        let base = self.integrate(f, previous)?;
        let fro: f64 = f.iter().map(|v| v * v).sum();
        let step = step_scale * fro.sqrt().max(1.0);
        let mut plus = *f;
        for i in 0..9 {
            plus[i] += step * direction[i];
        }
        let np = self.integrate(&plus, previous)?;
        let mut r = [0.0; 9];
        if central_difference {
            let mut minus = *f;
            for i in 0..9 {
                minus[i] -= step * direction[i];
            }
            let nm = self.integrate(&minus, previous)?;
            for i in 0..9 {
                r[i] = (np.first_piola[i] - nm.first_piola[i]) / (2.0 * step);
            }
        } else {
            for i in 0..9 {
                r[i] = (np.first_piola[i] - base.first_piola[i]) / step;
            }
        }
        Ok(r)
    }
    fn integrate(&self, f: &Tensor3, prev: &FiniteStrainJ2State) -> Result<FiniteStrainJ2Result> {
        require(
            f.iter().all(|v| v.is_finite()),
            "finite-strain J2 deformation must be finite",
        )?;
        require(
            prev.plastic_deformation_gradient
                .iter()
                .all(|v| v.is_finite()),
            "finite-strain J2 plastic state must be finite",
        )?;
        let det = determinant(f);
        let pdet = determinant(&prev.plastic_deformation_gradient);
        require(
            det > 0.0
                && det.is_finite()
                && pdet > 0.0
                && pdet.is_finite()
                && prev.equivalent_plastic_strain >= 0.0
                && prev.equivalent_plastic_strain.is_finite(),
            "finite-strain J2 requires positive finite deformation states",
        )?;
        let fe = multiply(f, &inverse(&prev.plastic_deformation_gradient, pdet));
        let b = multiply(&fe, &transpose(&fe));
        let (eig, vec) = diagonalize_symmetric3(&b, self.config.spectral_iterations);
        require(
            eig.iter().all(|e| *e > f64::EPSILON && e.is_finite()),
            "finite-strain J2 spectral decomposition failed",
        )?;
        let mut ls = [0.0; 3];
        for i in 0..3 {
            ls[i] = 0.5 * eig[i].ln();
        }
        let c = &self.config;
        let shear = c.youngs_modulus / (2.0 * (1.0 + c.poisson_ratio));
        let bulk = c.youngs_modulus / (3.0 * (1.0 - 2.0 * c.poisson_ratio));
        let mean = (ls[0] + ls[1] + ls[2]) / 3.0;
        let dev = [ls[0] - mean, ls[1] - mean, ls[2] - mean];
        let dn2: f64 = dev.iter().map(|v| v * v).sum();
        let q = (6.0 * shear * shear * dn2).sqrt();
        let yf = q - (c.yield_stress + c.hardening_modulus * prev.equivalent_plastic_strain);
        let mut state = *prev;
        let mut yielded = false;
        if yf > 0.0 && q > f64::EPSILON {
            let dg = yf / (3.0 * shear + c.hardening_modulus);
            let corr = 1.5 * dg / q;
            for i in 0..3 {
                ls[i] -= corr * (2.0 * shear * dev[i]);
            }
            state.equivalent_plastic_strain += dg;
            yielded = true;
        }
        let tr = ls[0] + ls[1] + ls[2];
        let um = tr / 3.0;
        let udn2: f64 = ls.iter().map(|v| (v - um) * (v - um)).sum();
        let estretch = reconstruct(&vec, &[ls[0].exp(), ls[1].exp(), ls[2].exp()]);
        let tstretch = reconstruct(
            &vec,
            &[
                (0.5 * eig[0].ln()).exp(),
                (0.5 * eig[1].ln()).exp(),
                (0.5 * eig[2].ln()).exp(),
            ],
        );
        let rot = multiply(&fe, &inverse(&tstretch, determinant(&tstretch)));
        let eu = multiply(&rot, &estretch);
        state.plastic_deformation_gradient = multiply(&inverse(&eu, determinant(&eu)), f);
        let mut pk = [0.0; 3];
        for i in 0..3 {
            pk[i] = 2.0 * shear * (ls[i] - um) + 3.0 * bulk * um;
        }
        let tau = reconstruct(&vec, &pk);
        let fit = transpose(&inverse(f, det));
        let mut cauchy = [0.0; 9];
        for i in 0..9 {
            cauchy[i] = (1.0 / det) * tau[i];
        }
        Ok(FiniteStrainJ2Result {
            determinant: det,
            energy_density: shear * udn2
                + 0.5 * bulk * tr * tr
                + 0.5
                    * c.hardening_modulus
                    * state.equivalent_plastic_strain
                    * state.equivalent_plastic_strain,
            trial_equivalent_stress: q,
            yield_function: yf,
            yielded,
            cauchy_stress: cauchy,
            first_piola: multiply(&tau, &fit),
            state,
        })
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Tet4FiniteStrainResult {
    pub reference_volume: f64,
    pub determinant: f64,
    pub strain_energy: f64,
    pub internal_force: [f64; 12],
    pub tangent_action: [f64; 12],
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Tet4FiniteStrainPlasticResult {
    pub reference_volume: f64,
    pub determinant: f64,
    pub strain_energy: f64,
    pub trial_equivalent_stress: f64,
    pub yield_function: f64,
    pub yielded: bool,
    pub internal_force: [f64; 12],
    pub tangent_action: [f64; 12],
    pub state: FiniteStrainJ2State,
}

struct Tet4Kinematics {
    volume: f64,
    gradients: [[f64; 3]; 4],
    f: Tensor3,
    df: Tensor3,
}

fn mat_vec(m: &Tensor3, v: &[f64; 3]) -> [f64; 3] {
    [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    ]
}

fn tet4_kinematics(
    reference: &[Vec3; 4],
    current: &[Vec3; 4],
    direction: &[f64; 12],
    what: &str,
) -> Result<Tet4Kinematics> {
    for n in reference.iter().chain(current) {
        require(
            n.x.is_finite() && n.y.is_finite() && n.z.is_finite(),
            "Tet4 nodes must be finite",
        )?;
    }
    require(
        direction.iter().all(|v| v.is_finite()),
        format!("{what}Tet4 direction must be finite"),
    )?;
    let r = reference;
    let j = [
        r[1].x - r[0].x,
        r[2].x - r[0].x,
        r[3].x - r[0].x,
        r[1].y - r[0].y,
        r[2].y - r[0].y,
        r[3].y - r[0].y,
        r[1].z - r[0].z,
        r[2].z - r[0].z,
        r[3].z - r[0].z,
    ];
    let det = determinant(&j);
    require(
        det > 0.0 && det.is_finite(),
        format!("{what}Tet4 reference nodes must have positive orientation"),
    )?;
    let jit = transpose(&inverse(&j, det));
    const NAT: [[f64; 3]; 4] = [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut g = [[0.0; 3]; 4];
    for n in 0..4 {
        g[n] = mat_vec(&jit, &NAT[n]);
    }
    let (mut f, mut df) = ([0.0; 9], [0.0; 9]);
    for n in 0..4 {
        let c = [current[n].x, current[n].y, current[n].z];
        for row in 0..3 {
            for col in 0..3 {
                f[3 * row + col] += c[row] * g[n][col];
                df[3 * row + col] += direction[3 * n + row] * g[n][col];
            }
        }
    }
    Ok(Tet4Kinematics {
        volume: det / 6.0,
        gradients: g,
        f,
        df,
    })
}

fn tet4_forces(k: &Tet4Kinematics, p: &Tensor3, t: &Tensor3) -> ([f64; 12], [f64; 12]) {
    let (mut force, mut tangent) = ([0.0; 12], [0.0; 12]);
    for n in 0..4 {
        for c in 0..3 {
            let (mut a, mut b) = (0.0, 0.0);
            for d in 0..3 {
                a += p[3 * c + d] * k.gradients[n][d];
                b += t[3 * c + d] * k.gradients[n][d];
            }
            force[3 * n + c] = k.volume * a;
            tangent[3 * n + c] = k.volume * b;
        }
    }
    (force, tangent)
}

/// Single-point total-Lagrangian Tet4; `tangent_action` is `K * direction` without forming K.
pub fn assemble_total_lagrangian_tet4(
    reference_nodes: &[Vec3; 4],
    current_nodes: &[Vec3; 4],
    displacement_direction: &[f64; 12],
    material: &CompressibleNeoHookean3D,
) -> Result<Tet4FiniteStrainResult> {
    let k = tet4_kinematics(reference_nodes, current_nodes, displacement_direction, "")?;
    let s = material.evaluate(&k.f)?;
    let t = material.tangent_action(&k.f, &k.df)?;
    let (internal_force, tangent_action) = tet4_forces(&k, &s.first_piola, &t);
    Ok(Tet4FiniteStrainResult {
        reference_volume: k.volume,
        determinant: s.determinant,
        strain_energy: k.volume * s.energy_density,
        internal_force,
        tangent_action,
    })
}

/// Total-Lagrangian Tet4 bridge for [`FiniteStrainJ2`] with a matrix-free algorithmic tangent.
pub fn assemble_total_lagrangian_plastic_tet4(
    reference_nodes: &[Vec3; 4],
    current_nodes: &[Vec3; 4],
    displacement_direction: &[f64; 12],
    material: &FiniteStrainJ2,
    previous: &FiniteStrainJ2State,
) -> Result<Tet4FiniteStrainPlasticResult> {
    let k = tet4_kinematics(
        reference_nodes,
        current_nodes,
        displacement_direction,
        "plastic ",
    )?;
    let s = material.evaluate(&k.f, previous)?;
    let t = material.tangent_action(&k.f, &k.df, previous)?;
    let (internal_force, tangent_action) = tet4_forces(&k, &s.first_piola, &t);
    Ok(Tet4FiniteStrainPlasticResult {
        reference_volume: k.volume,
        determinant: s.determinant,
        strain_energy: k.volume * s.energy_density,
        trial_equivalent_stress: s.trial_equivalent_stress,
        yield_function: s.yield_function,
        yielded: s.yielded,
        internal_force,
        tangent_action,
        state: s.state,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Tri3FiniteStrainPlasticResult {
    pub reference_area: f64,
    pub determinant: f64,
    pub strain_energy: f64,
    pub trial_equivalent_stress: f64,
    pub yield_function: f64,
    pub yielded: bool,
    pub internal_force: [f64; 6],
    pub tangent_action: [f64; 6],
    pub state: FiniteStrainJ2State,
}

/// Plane-strain total-Lagrangian Tri3: the in-plane gradient is embedded with `F_zz = 1`
/// (upstream default thickness 1).
pub fn assemble_total_lagrangian_plastic_tri3(
    reference_nodes: &[Node2; 3],
    current_nodes: &[Node2; 3],
    displacement_direction: &[f64; 6],
    material: &FiniteStrainJ2,
    thickness: f64,
    previous: &FiniteStrainJ2State,
) -> Result<Tri3FiniteStrainPlasticResult> {
    require(
        thickness > 0.0 && thickness.is_finite(),
        "plastic Tri3 thickness must be positive and finite",
    )?;
    require(
        reference_nodes
            .iter()
            .all(|n| n.x.is_finite() && n.y.is_finite()),
        "plastic Tri3 reference nodes must be finite",
    )?;
    require(
        current_nodes
            .iter()
            .all(|n| n.x.is_finite() && n.y.is_finite()),
        "plastic Tri3 current nodes must be finite",
    )?;
    require(
        displacement_direction.iter().all(|v| v.is_finite()),
        "plastic Tri3 direction must be finite",
    )?;
    let r = reference_nodes;
    let (j11, j12, j21, j22) = (
        r[1].x - r[0].x,
        r[2].x - r[0].x,
        r[1].y - r[0].y,
        r[2].y - r[0].y,
    );
    let det = j11 * j22 - j12 * j21;
    require(
        det > 0.0 && det.is_finite(),
        "plastic Tri3 reference nodes must have positive orientation",
    )?;
    let id = 1.0 / det;
    let g = [
        [(r[1].y - r[2].y) * id, (r[2].x - r[1].x) * id],
        [(r[2].y - r[0].y) * id, (r[0].x - r[2].x) * id],
        [(r[0].y - r[1].y) * id, (r[1].x - r[0].x) * id],
    ];
    let c = current_nodes;
    let (c11, c12, c21, c22) = (
        c[1].x - c[0].x,
        c[2].x - c[0].x,
        c[1].y - c[0].y,
        c[2].y - c[0].y,
    );
    let d = displacement_direction;
    let (d11, d12, d21, d22) = (d[2] - d[0], d[4] - d[0], d[3] - d[1], d[5] - d[1]);
    let (r11, r12, r21, r22) = (j22 * id, -j12 * id, -j21 * id, j11 * id);
    let embed = |a11: f64, a12: f64, a21: f64, a22: f64, a33: f64| {
        [a11, a12, 0.0, a21, a22, 0.0, 0.0, 0.0, a33]
    };
    let f = embed(
        c11 * r11 + c12 * r21,
        c11 * r12 + c12 * r22,
        c21 * r11 + c22 * r21,
        c21 * r12 + c22 * r22,
        1.0,
    );
    let df = embed(
        d11 * r11 + d12 * r21,
        d11 * r12 + d12 * r22,
        d21 * r11 + d22 * r21,
        d21 * r12 + d22 * r22,
        0.0,
    );
    let s = material.evaluate(&f, previous)?;
    let t = material.tangent_action(&f, &df, previous)?;
    let area = 0.5 * det;
    let measure = area * thickness;
    let p = &s.first_piola;
    let mut out = Tri3FiniteStrainPlasticResult {
        reference_area: area,
        determinant: s.determinant,
        strain_energy: measure * s.energy_density,
        trial_equivalent_stress: s.trial_equivalent_stress,
        yield_function: s.yield_function,
        yielded: s.yielded,
        state: s.state,
        ..Default::default()
    };
    for n in 0..3 {
        out.internal_force[2 * n] = measure * (p[0] * g[n][0] + p[1] * g[n][1]);
        out.internal_force[2 * n + 1] = measure * (p[3] * g[n][0] + p[4] * g[n][1]);
        out.tangent_action[2 * n] = measure * (t[0] * g[n][0] + t[1] * g[n][1]);
        out.tangent_action[2 * n + 1] = measure * (t[3] * g[n][0] + t[4] * g[n][1]);
    }
    Ok(out)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FiniteStrainJ2CalibrationReport {
    pub determinant: f64,
    pub rotation_relative_error: f64,
    pub hydrostatic_deviatoric_ratio: f64,
    pub tangent_step_relative_change: f64,
    pub positive_state: bool,
    pub objective_rotation: bool,
    pub hydrostatic_response: bool,
    pub tangent_step_converged: bool,
}

impl FiniteStrainJ2CalibrationReport {
    pub fn passed(&self) -> bool {
        self.positive_state
            && self.objective_rotation
            && self.hydrostatic_response
            && self.tangent_step_converged
    }
}

fn frobenius(t: &Tensor3) -> f64 {
    t.iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// Material-point audit: positive state, superposed rotation, hydrostatic response and central
/// tangent convergence under step halving.
pub fn calibrate_finite_strain_j2(
    material: &FiniteStrainJ2,
    f: &Tensor3,
    direction: &Tensor3,
    previous: &FiniteStrainJ2State,
) -> Result<FiniteStrainJ2CalibrationReport> {
    let base = material.evaluate(f, previous)?;
    let base_norm = frobenius(&base.cauchy_stress).max(1.0);
    let angle: f64 = 0.37;
    let rot = [
        angle.cos(),
        -angle.sin(),
        0.0,
        angle.sin(),
        angle.cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    ];
    let rotated = material.evaluate(&multiply(&rot, f), previous)?;
    let expected = multiply(&rot, &multiply(&base.cauchy_stress, &transpose(&rot)));
    let rotation_error = rotated
        .cauchy_stress
        .iter()
        .zip(&expected)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
    let hydro = material.evaluate(&[1.05, 0.0, 0.0, 0.0, 1.05, 0.0, 0.0, 0.0, 1.05], previous)?;
    let hnorm = frobenius(&hydro.cauchy_stress).max(1.0);
    let mut dev = hydro.cauchy_stress;
    let mean = (dev[0] + dev[4] + dev[8]) / 3.0;
    dev[0] -= mean;
    dev[4] -= mean;
    dev[8] -= mean;
    let ratio = frobenius(&dev) / hnorm;
    let t = material.central_tangent_action(f, direction, previous)?;
    let tr = material.tangent_action_with_step(f, direction, previous, 0.5 * 1e-7, true)?;
    let (mut diff, mut scale): (f64, f64) = (0.0, 1.0);
    for i in 0..9 {
        diff = diff.max((t[i] - tr[i]).abs());
        scale = scale.max(tr[i].abs());
    }
    let mut r = FiniteStrainJ2CalibrationReport {
        determinant: base.determinant,
        rotation_relative_error: rotation_error / base_norm,
        hydrostatic_deviatoric_ratio: ratio,
        tangent_step_relative_change: diff / scale,
        ..Default::default()
    };
    r.positive_state = base.determinant > 0.0 && base.determinant.is_finite();
    r.objective_rotation = r.rotation_relative_error < 1e-8;
    r.hydrostatic_response = r.hydrostatic_deviatoric_ratio < 1e-10;
    r.tangent_step_converged = r.tangent_step_relative_change < 1e-3;
    Ok(r)
}
