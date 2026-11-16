//! Native segregated RANS transport from the pinned upstream implementation.
//! Coefficients are lagged per step. All equations commit together after successful solves.
use super::{models::velocity_gradients, *};
use crate::{core::*, finite_slice, positive, require, Result, SolverError};
use serde::{Deserialize, Serialize};

pub type TurbulenceScalarBoundaryType = ScalarBoundaryType;
pub type TurbulenceScalarBoundaryCondition = ScalarBoundaryCondition;
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RansTransportControls {
    pub dt: f64,
    pub linear: LinearConfig,
    pub temporal_scheme: TemporalScheme,
    pub crank_nicolson_off_centering: f64,
}
impl Default for RansTransportControls {
    fn default() -> Self {
        Self {
            dt: 1e-4,
            linear: LinearConfig {
                max_iterations: 400,
                restart: 30,
                relative_tolerance: 1e-10,
            },
            temporal_scheme: TemporalScheme::Euler,
            crank_nicolson_off_centering: 1.0,
        }
    }
}
impl RansTransportControls {
    fn validate(self) -> Result<()> {
        positive(self.dt, "RANS dt")?;
        positive(self.linear.relative_tolerance, "RANS linear tolerance")?;
        require(
            self.linear.max_iterations > 0 && self.linear.restart > 0,
            "RANS linear iterations/restart",
        )?;
        if self.temporal_scheme == TemporalScheme::CrankNicolson {
            crank_nicolson_implicit_weight(self.crank_nicolson_off_centering)?;
        }
        Ok(())
    }
}
fn checked_field(v: &[f64], n: usize, floor: f64) -> Result<Vec<f64>> {
    require(v.len() == n, "RANS field size")?;
    finite_slice(v, "RANS field")?;
    require(
        v.iter().all(|&x| x >= floor),
        "RANS field below allowed input floor",
    )?;
    Ok(v.to_vec())
}
#[derive(Clone, Debug)]
struct RansState<const N: usize> {
    mesh: PolyMesh,
    controls: RansTransportControls,
    fields: [Vec<f64>; N],
    previous: [Vec<f64>; N],
    boundary: [Vec<ScalarBoundaryCondition>; N],
    flux: Vec<f64>,
    steps: usize,
    time: f64,
    last: Option<[IterativeSolverResult; N]>,
}
impl<const N: usize> RansState<N> {
    fn new(mesh: PolyMesh, controls: RansTransportControls, initial: [f64; N]) -> Result<Self> {
        controls.validate()?;
        let fields = std::array::from_fn(|q| vec![initial[q]; mesh.cell_count()]);
        Ok(Self {
            boundary: std::array::from_fn(|_| {
                vec![ScalarBoundaryCondition::default(); mesh.patches().len()]
            }),
            flux: vec![0.0; mesh.face_count()],
            previous: fields.clone(),
            fields,
            mesh,
            controls,
            steps: 0,
            time: 0.0,
            last: None,
        })
    }
    fn initialize(&mut self, fields: [Vec<f64>; N]) {
        self.previous = fields.clone();
        self.fields = fields;
        self.steps = 0;
        self.time = 0.0;
        self.last = None;
    }
    fn set_boundary(
        &mut self,
        q: usize,
        patch: &str,
        kind: ScalarBoundaryType,
        value: f64,
        floor: f64,
    ) -> Result<()> {
        require(
            value.is_finite() && (kind != ScalarBoundaryType::FixedValue || value >= floor),
            "invalid RANS boundary value",
        )?;
        let p = self.mesh.patch_index(patch)?;
        self.boundary[q][p] = ScalarBoundaryCondition {
            boundary_type: kind,
            value,
        };
        Ok(())
    }
    fn set_flux(&mut self, flux: Vec<f64>) -> Result<()> {
        require(flux.len() == self.mesh.face_count(), "RANS face flux size")?;
        finite_slice(&flux, "RANS flux")?;
        self.flux = flux;
        Ok(())
    }
    fn advance(
        &mut self,
        gamma: [Vec<f64>; N],
        source: [Vec<f64>; N],
        sink: [Vec<f64>; N],
        floors: [f64; N],
    ) -> Result<[IterativeSolverResult; N]> {
        let mut next = self.fields.clone();
        let mut results = [IterativeSolverResult {
            iterations: 0,
            residual_rms: 0.0,
            converged: false,
        }; N];
        for q in 0..N {
            let (values, result) = advance_equation(
                &self.mesh,
                self.controls,
                &self.boundary[q],
                &self.flux,
                &gamma[q],
                &source[q],
                &sink[q],
                &self.fields[q],
                &self.previous[q],
                self.steps,
                floors[q],
            )?;
            next[q] = values;
            results[q] = result;
        }
        let time = self.time + self.controls.dt;
        require(time.is_finite(), "RANS clock overflow")?;
        self.previous = self.fields.clone();
        self.fields = next;
        self.steps += 1;
        self.time = time;
        self.last = Some(results);
        Ok(results)
    }
}
#[allow(clippy::too_many_arguments)]
fn advance_equation(
    mesh: &PolyMesh,
    c: RansTransportControls,
    bc: &[ScalarBoundaryCondition],
    flux: &[f64],
    gamma: &[f64],
    source: &[f64],
    sink: &[f64],
    old: &[f64],
    older: &[f64],
    steps: usize,
    floor: f64,
) -> Result<(Vec<f64>, IterativeSolverResult)> {
    let n = mesh.cell_count();
    let mut b = CsrBuilder::new(n, n);
    let mut boundary_rhs = vec![0.0; n];
    for i in 0..n {
        require(
            gamma[i].is_finite()
                && gamma[i] >= 0.0
                && source[i].is_finite()
                && sink[i].is_finite()
                && sink[i] >= 0.0,
            "invalid RANS equation coefficients",
        )?;
        let mut diagonal = sink[i] * mesh.cells()[i].volume;
        for &fi in &mesh.cell_faces()[i] {
            let f = mesh.faces()[fi];
            let phi = if f.owner == i { flux[fi] } else { -flux[fi] };
            if let Some(neighbour) = f.neighbour {
                let j = if f.owner == i { neighbour } else { f.owner };
                let (a, z) = (gamma[i], gamma[j]);
                let harmonic = if a > 0.0 && z > 0.0 {
                    2.0 * a * z / (a + z)
                } else {
                    (0.5 * (a + z)).max(0.0)
                };
                let d = harmonic * decompose_face_area(mesh, fi)?.orthogonal_metric;
                diagonal += d + phi.max(0.0);
                b.add(i, j, -d + phi.min(0.0))?;
            } else {
                let condition = bc[f.patch.ok_or_else(|| {
                    crate::SolverError::InvalidInput("boundary face without patch".into())
                })?];
                if condition.boundary_type == ScalarBoundaryType::FixedValue {
                    require(
                        condition.value >= floor,
                        "RANS fixed boundary below field floor",
                    )?;
                    let distance = (f.center - mesh.cells()[i].center).norm();
                    positive(distance, "RANS boundary distance")?;
                    let d = gamma[i] * f.area.norm() / distance;
                    diagonal += d + phi.max(0.0);
                    boundary_rhs[i] += d * condition.value - phi.min(0.0) * condition.value;
                } else {
                    diagonal += phi;
                }
            }
        }
        b.add(i, i, diagonal)?;
    }
    let spatial = b.build(0.0)?;
    let bdf = c.temporal_scheme == TemporalScheme::BackwardBdf2 && steps > 0;
    let alpha = if bdf { 1.5 } else { 1.0 };
    let theta = if c.temporal_scheme == TemporalScheme::CrankNicolson {
        crank_nicolson_implicit_weight(c.crank_nicolson_off_centering)?
    } else {
        1.0
    };
    let mut a = CsrBuilder::new(n, n);
    for i in 0..n {
        for p in spatial.row_offsets()[i]..spatial.row_offsets()[i + 1] {
            a.add(i, spatial.columns()[p], theta * spatial.values()[p])?;
        }
        a.add(i, i, alpha * mesh.cells()[i].volume / c.dt)?;
    }
    let a = a.build(0.0)?;
    let pre = Ilu0Preconditioner::new(&a, 1e-30)?;
    let explicit = if theta < 1.0 {
        spatial.apply(old)?
    } else {
        vec![0.0; n]
    };
    let rhs: Vec<_> = (0..n)
        .map(|i| {
            let temporal = if bdf {
                2.0 * old[i] - 0.5 * older[i]
            } else {
                old[i]
            };
            temporal * mesh.cells()[i].volume / c.dt
                + source[i] * mesh.cells()[i].volume
                + boundary_rhs[i]
                - (1.0 - theta) * explicit[i]
        })
        .collect();
    let mut next = old.to_vec();
    let r = restarted_gmres(&rhs, &mut next, |x, y| a.multiply(x, y), &pre, c.linear)?;
    if !r.converged {
        return Err(SolverError::NumericalFailure(format!(
            "RANS scalar solve did not converge: {}",
            r.residual_rms
        )));
    }
    finite_slice(&next, "RANS candidate")?;
    for x in &mut next {
        *x = x.max(floor);
    }
    Ok((next, r))
}
pub fn turbulence_strain_rate_magnitude(mesh: &PolyMesh, u: &[Vec3]) -> Result<Vec<f64>> {
    let g = velocity_gradients(mesh, u)?;
    Ok(g.iter()
        .map(|g| {
            let mut ss = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    ss += (0.5 * (g[i][j] + g[j][i])).powi(2);
                }
            }
            (2.0 * ss).sqrt()
        })
        .collect())
}
fn check_constants(positive_values: &[f64], nonnegative: &[f64]) -> Result<()> {
    require(
        positive_values.iter().all(|x| x.is_finite() && *x > 0.0)
            && nonnegative.iter().all(|x| x.is_finite() && *x >= 0.0),
        "invalid RANS model constants",
    )
}
macro_rules! state_accessors {
    () => {
        pub fn mesh(&self) -> &PolyMesh {
            &self.state.mesh
        }
        pub fn time(&self) -> f64 {
            self.state.time
        }
        pub fn set_face_flux(&mut self, flux: Vec<f64>) -> Result<()> {
            self.state.set_flux(flux)
        }
    };
}
macro_rules! strain_accessors {
    () => {
        pub fn strain_rate(&self) -> &[f64] {
            &self.strain
        }
        pub fn set_strain_rate(&mut self, value: f64) -> Result<()> {
            self.set_strain_rate_profile(&vec![value; self.state.mesh.cell_count()])
        }
        pub fn set_strain_rate_profile(&mut self, values: &[f64]) -> Result<()> {
            self.strain = checked_field(values, self.state.mesh.cell_count(), 0.0)?;
            Ok(())
        }
        pub fn set_velocity(&mut self, u: &[Vec3]) -> Result<()> {
            self.strain = turbulence_strain_rate_magnitude(&self.state.mesh, u)?;
            Ok(())
        }
    };
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SpalartAllmarasConfig {
    pub transport: RansTransportControls,
    pub molecular_viscosity: f64,
    pub cb1: f64,
    pub cb2: f64,
    pub cw2: f64,
    pub cw3: f64,
    pub cv1: f64,
    pub cs: f64,
    pub sigma_nu_tilde: f64,
    pub kappa: f64,
    pub minimum_nu_tilde: f64,
    pub minimum_wall_distance: f64,
}
impl Default for SpalartAllmarasConfig {
    fn default() -> Self {
        Self {
            transport: Default::default(),
            molecular_viscosity: 1.5e-5,
            cb1: 0.1355,
            cb2: 0.622,
            cw2: 0.3,
            cw3: 2.0,
            cv1: 7.1,
            cs: 0.3,
            sigma_nu_tilde: 2.0 / 3.0,
            kappa: 0.41,
            minimum_nu_tilde: 1e-12,
            minimum_wall_distance: 1e-9,
        }
    }
}
#[derive(Clone, Debug)]
pub struct SpalartAllmarasTransport {
    state: RansState<1>,
    config: SpalartAllmarasConfig,
    wall: Vec<f64>,
    strain: Vec<f64>,
}
impl SpalartAllmarasTransport {
    pub fn new(mesh: PolyMesh, c: SpalartAllmarasConfig) -> Result<Self> {
        check_constants(
            &[
                c.molecular_viscosity,
                c.cb1,
                c.cw3,
                c.cv1,
                c.sigma_nu_tilde,
                c.kappa,
                c.minimum_wall_distance,
            ],
            &[c.cb2, c.cw2, c.cs, c.minimum_nu_tilde],
        )?;
        let n = mesh.cell_count();
        Ok(Self {
            state: RansState::new(mesh, c.transport, [c.minimum_nu_tilde])?,
            config: c,
            wall: vec![1.0; n],
            strain: vec![0.0; n],
        })
    }
    state_accessors!();
    strain_accessors!();
    pub fn nu_tilde(&self) -> &[f64] {
        &self.state.fields[0]
    }
    pub fn wall_distance(&self) -> &[f64] {
        &self.wall
    }
    pub fn linear_result(&self) -> Option<IterativeSolverResult> {
        self.state.last.map(|r| r[0])
    }
    pub fn initialize(&mut self, value: f64) -> Result<()> {
        require(
            value.is_finite() && value >= self.config.minimum_nu_tilde,
            "SA initial value",
        )?;
        self.initialize_profile(&vec![value; self.mesh().cell_count()])
    }
    pub fn initialize_profile(&mut self, values: &[f64]) -> Result<()> {
        let mut v = checked_field(values, self.mesh().cell_count(), 0.0)?;
        for x in &mut v {
            *x = x.max(self.config.minimum_nu_tilde);
        }
        self.state.initialize([v]);
        Ok(())
    }
    pub fn set_boundary(
        &mut self,
        patch: &str,
        kind: ScalarBoundaryType,
        value: f64,
    ) -> Result<()> {
        self.state
            .set_boundary(0, patch, kind, value, self.config.minimum_nu_tilde)
    }
    pub fn set_wall_distance(&mut self, value: f64) -> Result<()> {
        require(
            value >= self.config.minimum_wall_distance,
            "SA minimum wall distance",
        )?;
        self.set_wall_distance_profile(&vec![value; self.mesh().cell_count()])
    }
    pub fn set_wall_distance_profile(&mut self, values: &[f64]) -> Result<()> {
        let mut d = checked_field(values, self.mesh().cell_count(), f64::MIN_POSITIVE)?;
        for x in &mut d {
            *x = x.max(self.config.minimum_wall_distance);
        }
        self.wall = d;
        Ok(())
    }
    pub fn kinematic_eddy_viscosity(&self) -> Result<Vec<f64>> {
        self.nu_tilde()
            .iter()
            .map(|&v| {
                spalart_allmaras_eddy_viscosity(self.config.molecular_viscosity, v, self.config.cv1)
            })
            .collect()
    }
    pub fn step(&mut self) -> Result<IterativeSolverResult> {
        let c = self.config;
        let n = self.mesh().cell_count();
        let gradient = least_squares_gradient_scalar(self.mesh(), self.nu_tilde(), &[])?;
        let (mut diffusion, mut source, mut sink) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        let cw1 = c.cb1 / c.kappa.powi(2) + (1.0 + c.cb2) / c.sigma_nu_tilde;
        let cw36 = c.cw3.powf(6.0);
        for i in 0..n {
            let nt = self.nu_tilde()[i].max(c.minimum_nu_tilde);
            let d = self.wall[i].max(c.minimum_wall_distance);
            let chi = nt / c.molecular_viscosity;
            let chi3 = chi * chi * chi;
            let fv1 = chi3 / (chi3 + c.cv1 * c.cv1 * c.cv1);
            let fv2 = 1.0 - chi / (1.0 + chi * fv1);
            let kd = c.kappa * c.kappa * d * d;
            let st = (self.strain[i] + nt * fv2 / kd).max(c.cs * self.strain[i]);
            let r = (nt / (st * kd).max(1e-30)).clamp(0.0, 10.0);
            let g = r + c.cw2 * (r.powf(6.0) - r);
            let fw = g * ((1.0 + cw36) / (g.powf(6.0) + cw36).max(1e-30)).powf(1.0 / 6.0);
            diffusion[i] = (c.molecular_viscosity + nt) / c.sigma_nu_tilde;
            source[i] = c.cb1 * st * nt + c.cb2 / c.sigma_nu_tilde * gradient[i].dot(gradient[i]);
            sink[i] = cw1 * fw * nt / (d * d);
        }
        Ok(self
            .state
            .advance([diffusion], [source], [sink], [c.minimum_nu_tilde])?[0])
    }
    pub fn run(&mut self, steps: usize) -> Result<Option<IterativeSolverResult>> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(self.linear_result())
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct KEpsilonConfig {
    pub transport: RansTransportControls,
    pub density: f64,
    pub molecular_viscosity: f64,
    pub c_mu: f64,
    pub c1: f64,
    pub c2: f64,
    pub sigma_k: f64,
    pub sigma_epsilon: f64,
    pub production_limiter: f64,
    pub minimum_k: f64,
    pub minimum_epsilon: f64,
    pub maximum_eddy_viscosity_ratio: f64,
}
impl Default for KEpsilonConfig {
    fn default() -> Self {
        Self {
            transport: Default::default(),
            density: 1.0,
            molecular_viscosity: 1.5e-5,
            c_mu: 0.09,
            c1: 1.44,
            c2: 1.92,
            sigma_k: 1.0,
            sigma_epsilon: 1.3,
            production_limiter: 20.0,
            minimum_k: 1e-12,
            minimum_epsilon: 1e-12,
            maximum_eddy_viscosity_ratio: 1e6,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct TwoEquationLinearResult {
    pub first: IterativeSolverResult,
    pub second: IterativeSolverResult,
}
impl From<[IterativeSolverResult; 2]> for TwoEquationLinearResult {
    fn from(r: [IterativeSolverResult; 2]) -> Self {
        Self {
            first: r[0],
            second: r[1],
        }
    }
}
#[derive(Clone, Debug)]
pub struct KEpsilonTransport {
    state: RansState<2>,
    config: KEpsilonConfig,
    strain: Vec<f64>,
}
impl KEpsilonTransport {
    pub fn new(mesh: PolyMesh, c: KEpsilonConfig) -> Result<Self> {
        check_constants(
            &[
                c.density,
                c.molecular_viscosity,
                c.c_mu,
                c.c1,
                c.c2,
                c.sigma_k,
                c.sigma_epsilon,
                c.production_limiter,
                c.minimum_k,
                c.minimum_epsilon,
                c.maximum_eddy_viscosity_ratio,
            ],
            &[],
        )?;
        let n = mesh.cell_count();
        Ok(Self {
            state: RansState::new(mesh, c.transport, [c.minimum_k, c.minimum_epsilon])?,
            config: c,
            strain: vec![0.0; n],
        })
    }
    state_accessors!();
    strain_accessors!();
    pub fn k(&self) -> &[f64] {
        &self.state.fields[0]
    }
    pub fn epsilon(&self) -> &[f64] {
        &self.state.fields[1]
    }
    pub fn linear_result(&self) -> Option<TwoEquationLinearResult> {
        self.state.last.map(Into::into)
    }
    pub fn initialize(&mut self, k: f64, epsilon: f64) -> Result<()> {
        require(
            k >= self.config.minimum_k && epsilon >= self.config.minimum_epsilon,
            "k-epsilon initial floors",
        )?;
        let n = self.mesh().cell_count();
        self.initialize_profile(&vec![k; n], &vec![epsilon; n])
    }
    pub fn initialize_profile(&mut self, k: &[f64], epsilon: &[f64]) -> Result<()> {
        let n = self.mesh().cell_count();
        let mut k = checked_field(k, n, f64::MIN_POSITIVE)?;
        let mut e = checked_field(epsilon, n, f64::MIN_POSITIVE)?;
        for x in &mut k {
            *x = x.max(self.config.minimum_k);
        }
        for x in &mut e {
            *x = x.max(self.config.minimum_epsilon);
        }
        self.state.initialize([k, e]);
        Ok(())
    }
    pub fn set_k_boundary(
        &mut self,
        patch: &str,
        kind: ScalarBoundaryType,
        value: f64,
    ) -> Result<()> {
        self.state
            .set_boundary(0, patch, kind, value, self.config.minimum_k)
    }
    pub fn set_epsilon_boundary(
        &mut self,
        patch: &str,
        kind: ScalarBoundaryType,
        value: f64,
    ) -> Result<()> {
        self.state
            .set_boundary(1, patch, kind, value, self.config.minimum_epsilon)
    }
    pub fn kinematic_eddy_viscosity(&self) -> Result<Vec<f64>> {
        let c = self.config;
        self.k()
            .iter()
            .zip(self.epsilon())
            .map(|(&k, &e)| {
                Ok(
                    (k_epsilon_eddy_viscosity(c.density, k, e, c.c_mu)? / c.density)
                        .min(c.maximum_eddy_viscosity_ratio * c.molecular_viscosity),
                )
            })
            .collect()
    }
    pub fn dynamic_eddy_viscosity(&self) -> Result<Vec<f64>> {
        Ok(self
            .kinematic_eddy_viscosity()?
            .iter()
            .map(|v| v * self.config.density)
            .collect())
    }
    pub fn step(&mut self) -> Result<TwoEquationLinearResult> {
        let c = self.config;
        let n = self.mesh().cell_count();
        let nut = self.kinematic_eddy_viscosity()?;
        let (mut gamma, mut source, mut sink) = (
            [vec![0.0; n], vec![0.0; n]],
            [vec![0.0; n], vec![0.0; n]],
            [vec![0.0; n], vec![0.0; n]],
        );
        for i in 0..n {
            let k = self.k()[i].max(c.minimum_k);
            let e = self.epsilon()[i].max(c.minimum_epsilon);
            let p = (nut[i] * self.strain[i] * self.strain[i]).min(c.production_limiter * e);
            gamma[0][i] = c.molecular_viscosity + nut[i] / c.sigma_k;
            gamma[1][i] = c.molecular_viscosity + nut[i] / c.sigma_epsilon;
            source[0][i] = p;
            source[1][i] = c.c1 * (e / k) * p;
            sink[0][i] = e / k;
            sink[1][i] = c.c2 * e / k;
        }
        Ok(self
            .state
            .advance(gamma, source, sink, [c.minimum_k, c.minimum_epsilon])?
            .into())
    }
    pub fn run(&mut self, steps: usize) -> Result<Option<TwoEquationLinearResult>> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(self.linear_result())
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SSTDESZonalFilter {
    None,
    F1,
    F2,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct KOmegaSSTConfig {
    pub transport: RansTransportControls,
    pub density: f64,
    pub molecular_viscosity: f64,
    pub beta_star: f64,
    pub beta1: f64,
    pub beta2: f64,
    pub gamma1: f64,
    pub gamma2: f64,
    pub alpha_k1: f64,
    pub alpha_k2: f64,
    pub alpha_omega1: f64,
    pub alpha_omega2: f64,
    pub a1: f64,
    pub production_limiter: f64,
    pub minimum_k: f64,
    pub minimum_omega: f64,
    pub minimum_wall_distance: f64,
    pub maximum_eddy_viscosity_ratio: f64,
    pub des_enabled: bool,
    pub c_des: f64,
    pub des_zonal_filter: SSTDESZonalFilter,
}
impl Default for KOmegaSSTConfig {
    fn default() -> Self {
        Self {
            transport: Default::default(),
            density: 1.0,
            molecular_viscosity: 1.5e-5,
            beta_star: 0.09,
            beta1: 0.075,
            beta2: 0.0828,
            gamma1: 5.0 / 9.0,
            gamma2: 0.44,
            alpha_k1: 0.85,
            alpha_k2: 1.0,
            alpha_omega1: 0.5,
            alpha_omega2: 0.856,
            a1: 0.31,
            production_limiter: 10.0,
            minimum_k: 1e-12,
            minimum_omega: 1e-12,
            minimum_wall_distance: 1e-9,
            maximum_eddy_viscosity_ratio: 1e6,
            des_enabled: false,
            c_des: 0.61,
            des_zonal_filter: SSTDESZonalFilter::F2,
        }
    }
}
#[derive(Clone, Debug)]
pub struct KOmegaSSTTransport {
    state: RansState<2>,
    config: KOmegaSSTConfig,
    strain: Vec<f64>,
    wall: Vec<f64>,
    grid: Vec<f64>,
    f1: Vec<f64>,
    f2: Vec<f64>,
    hybrid: Vec<f64>,
}
impl KOmegaSSTTransport {
    pub fn new(mesh: PolyMesh, c: KOmegaSSTConfig) -> Result<Self> {
        check_constants(
            &[
                c.density,
                c.molecular_viscosity,
                c.beta_star,
                c.beta1,
                c.beta2,
                c.gamma1,
                c.gamma2,
                c.alpha_k1,
                c.alpha_k2,
                c.alpha_omega1,
                c.alpha_omega2,
                c.a1,
                c.production_limiter,
                c.minimum_k,
                c.minimum_omega,
                c.minimum_wall_distance,
                c.maximum_eddy_viscosity_ratio,
                c.c_des,
            ],
            &[],
        )?;
        let n = mesh.cell_count();
        let grid = mesh.cells().iter().map(|c| c.volume.cbrt()).collect();
        let mut s = Self {
            state: RansState::new(mesh, c.transport, [c.minimum_k, c.minimum_omega])?,
            config: c,
            strain: vec![0.0; n],
            wall: vec![1.0; n],
            grid,
            f1: vec![0.0; n],
            f2: vec![0.0; n],
            hybrid: vec![1.0; n],
        };
        s.update_blending()?;
        Ok(s)
    }
    state_accessors!();
    strain_accessors!();
    pub fn k(&self) -> &[f64] {
        &self.state.fields[0]
    }
    pub fn omega(&self) -> &[f64] {
        &self.state.fields[1]
    }
    pub fn wall_distance(&self) -> &[f64] {
        &self.wall
    }
    pub fn blending_f1(&self) -> &[f64] {
        &self.f1
    }
    pub fn blending_f2(&self) -> &[f64] {
        &self.f2
    }
    pub fn hybrid_dissipation_factor(&self) -> &[f64] {
        &self.hybrid
    }
    pub fn linear_result(&self) -> Option<TwoEquationLinearResult> {
        self.state.last.map(Into::into)
    }
    pub fn initialize(&mut self, k: f64, omega: f64) -> Result<()> {
        require(
            k >= self.config.minimum_k && omega >= self.config.minimum_omega,
            "SST initial floors",
        )?;
        let n = self.mesh().cell_count();
        self.initialize_profile(&vec![k; n], &vec![omega; n])
    }
    pub fn initialize_profile(&mut self, k: &[f64], omega: &[f64]) -> Result<()> {
        let n = self.mesh().cell_count();
        let mut k = checked_field(k, n, f64::MIN_POSITIVE)?;
        let mut w = checked_field(omega, n, f64::MIN_POSITIVE)?;
        for x in &mut k {
            *x = x.max(self.config.minimum_k);
        }
        for x in &mut w {
            *x = x.max(self.config.minimum_omega);
        }
        self.state.initialize([k, w]);
        self.update_blending()
    }
    pub fn set_k_boundary(
        &mut self,
        patch: &str,
        kind: ScalarBoundaryType,
        value: f64,
    ) -> Result<()> {
        self.state
            .set_boundary(0, patch, kind, value, self.config.minimum_k)
    }
    pub fn set_omega_boundary(
        &mut self,
        patch: &str,
        kind: ScalarBoundaryType,
        value: f64,
    ) -> Result<()> {
        self.state
            .set_boundary(1, patch, kind, value, self.config.minimum_omega)
    }
    pub fn set_wall_distance(&mut self, value: f64) -> Result<()> {
        require(
            value >= self.config.minimum_wall_distance,
            "SST wall distance",
        )?;
        self.set_wall_distance_profile(&vec![value; self.mesh().cell_count()])
    }
    pub fn set_wall_distance_profile(&mut self, values: &[f64]) -> Result<()> {
        let mut d = checked_field(values, self.mesh().cell_count(), f64::MIN_POSITIVE)?;
        for x in &mut d {
            *x = x.max(self.config.minimum_wall_distance);
        }
        self.wall = d;
        self.update_blending()
    }
    pub fn set_grid_scale(&mut self, value: f64) -> Result<()> {
        self.set_grid_scale_profile(&vec![value; self.mesh().cell_count()])
    }
    pub fn set_grid_scale_profile(&mut self, values: &[f64]) -> Result<()> {
        self.grid = checked_field(values, self.mesh().cell_count(), f64::MIN_POSITIVE)?;
        self.update_blending()
    }
    fn update_blending(&mut self) -> Result<()> {
        let c = self.config;
        let gk = least_squares_gradient_scalar(self.mesh(), self.k(), &[])?;
        let gw = least_squares_gradient_scalar(self.mesh(), self.omega(), &[])?;
        for i in 0..self.mesh().cell_count() {
            let k = self.k()[i].max(c.minimum_k);
            let w = self.omega()[i].max(c.minimum_omega);
            let d = self.wall[i].max(c.minimum_wall_distance);
            let cross = gk[i].dot(gw[i]);
            let cd = (2.0 * c.alpha_omega2 * cross / w).max(1e-20);
            let a = k.sqrt() / (c.beta_star * w * d);
            let b = 500.0 * c.molecular_viscosity / (d * d * w);
            let z = 4.0 * c.alpha_omega2 * k / (cd * d * d);
            let arg = a.max(b).min(z).clamp(0.0, 100.0);
            let arg2 = arg * arg;
            self.f1[i] = (arg2 * arg2).tanh();
            let arg = (2.0 * a).max(b).clamp(0.0, 100.0);
            self.f2[i] = (arg * arg).tanh();
            self.hybrid[i] = 1.0;
            if c.des_enabled {
                let length = k.sqrt() / (c.beta_star * w);
                let les = c.c_des * self.grid[i].max(1e-30);
                let raw = (length / les).max(1.0);
                let shield = match c.des_zonal_filter {
                    SSTDESZonalFilter::None => 1.0,
                    SSTDESZonalFilter::F1 => 1.0 - self.f1[i],
                    SSTDESZonalFilter::F2 => 1.0 - self.f2[i],
                };
                self.hybrid[i] = 1.0 + (raw - 1.0) * shield.clamp(0.0, 1.0);
            }
        }
        finite_slice(&self.f1, "SST F1")?;
        finite_slice(&self.f2, "SST F2")?;
        finite_slice(&self.hybrid, "SST DES factor")
    }
    pub fn kinematic_eddy_viscosity(&self) -> Result<Vec<f64>> {
        let c = self.config;
        (0..self.mesh().cell_count())
            .map(|i| {
                Ok((k_omega_sst_eddy_viscosity(
                    c.density,
                    self.k()[i],
                    self.omega()[i],
                    self.strain[i],
                    self.f2[i],
                    c.a1,
                )? / c.density)
                    .min(c.maximum_eddy_viscosity_ratio * c.molecular_viscosity))
            })
            .collect()
    }
    pub fn dynamic_eddy_viscosity(&self) -> Result<Vec<f64>> {
        Ok(self
            .kinematic_eddy_viscosity()?
            .iter()
            .map(|v| v * self.config.density)
            .collect())
    }
    pub fn step(&mut self) -> Result<TwoEquationLinearResult> {
        self.update_blending()?;
        let c = self.config;
        let n = self.mesh().cell_count();
        let nut = self.kinematic_eddy_viscosity()?;
        let gk = least_squares_gradient_scalar(self.mesh(), self.k(), &[])?;
        let gw = least_squares_gradient_scalar(self.mesh(), self.omega(), &[])?;
        let (mut gamma, mut source, mut sink) = (
            [vec![0.0; n], vec![0.0; n]],
            [vec![0.0; n], vec![0.0; n]],
            [vec![0.0; n], vec![0.0; n]],
        );
        for i in 0..n {
            let k = self.k()[i].max(c.minimum_k);
            let w = self.omega()[i].max(c.minimum_omega);
            let f = self.f1[i];
            let ak = f * c.alpha_k1 + (1.0 - f) * c.alpha_k2;
            let aw = f * c.alpha_omega1 + (1.0 - f) * c.alpha_omega2;
            let beta = f * c.beta1 + (1.0 - f) * c.beta2;
            let gm = f * c.gamma1 + (1.0 - f) * c.gamma2;
            let p = (nut[i] * self.strain[i] * self.strain[i])
                .min(c.production_limiter * c.beta_star * k * w);
            let ratio = if nut[i] > 1e-30 { p / nut[i] } else { 0.0 };
            let cross = 2.0 * (1.0 - f) * c.alpha_omega2 * gk[i].dot(gw[i]) / w;
            gamma[0][i] = c.molecular_viscosity + ak * nut[i];
            gamma[1][i] = c.molecular_viscosity + aw * nut[i];
            source[0][i] = p;
            sink[0][i] = c.beta_star * w * self.hybrid[i];
            source[1][i] = gm * ratio + cross.max(0.0);
            sink[1][i] = beta * w + if cross < 0.0 { -cross / w } else { 0.0 };
        }
        let result = self
            .state
            .advance(gamma, source, sink, [c.minimum_k, c.minimum_omega])?;
        self.update_blending()?;
        Ok(result.into())
    }
    pub fn run(&mut self, steps: usize) -> Result<Option<TwoEquationLinearResult>> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(self.linear_result())
    }
}
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricTensor3 {
    pub xx: f64,
    pub yy: f64,
    pub zz: f64,
    pub xy: f64,
    pub xz: f64,
    pub yz: f64,
}
impl SymmetricTensor3 {
    pub fn components(self) -> [f64; 6] {
        [self.xx, self.yy, self.zz, self.xy, self.xz, self.yz]
    }
    pub fn from_components(v: [f64; 6]) -> Self {
        Self {
            xx: v[0],
            yy: v[1],
            zz: v[2],
            xy: v[3],
            xz: v[4],
            yz: v[5],
        }
    }
}
pub type VelocityGradient3 = [f64; 9];
pub fn turbulent_kinetic_energy(r: SymmetricTensor3) -> f64 {
    0.5 * (r.xx + r.yy + r.zz)
}
pub fn reynolds_stress_is_realizable(r: SymmetricTensor3, tolerance: f64) -> bool {
    if !tolerance.is_finite() || tolerance < 0.0 || !r.components().iter().all(|v| v.is_finite()) {
        return false;
    }
    if r.xx < -tolerance || r.yy < -tolerance || r.zz < -tolerance {
        return false;
    }
    if r.xx * r.yy - r.xy * r.xy < -tolerance
        || r.xx * r.zz - r.xz * r.xz < -tolerance
        || r.yy * r.zz - r.yz * r.yz < -tolerance
    {
        return false;
    }
    let d = r.xx * r.yy * r.zz + 2.0 * r.xy * r.xz * r.yz
        - r.xx * r.yz * r.yz
        - r.yy * r.xz * r.xz
        - r.zz * r.xy * r.xy;
    d >= -tolerance
}
fn project_realizable(mut r: SymmetricTensor3, floor: f64) -> SymmetricTensor3 {
    r.xx = r.xx.max(floor);
    r.yy = r.yy.max(floor);
    r.zz = r.zz.max(floor);
    let margin = 1.0 - 1e-12;
    r.xy = r.xy.clamp(
        -margin * (r.xx * r.yy).sqrt(),
        margin * (r.xx * r.yy).sqrt(),
    );
    r.xz = r.xz.clamp(
        -margin * (r.xx * r.zz).sqrt(),
        margin * (r.xx * r.zz).sqrt(),
    );
    r.yz = r.yz.clamp(
        -margin * (r.yy * r.zz).sqrt(),
        margin * (r.yy * r.zz).sqrt(),
    );
    if reynolds_stress_is_realizable(r, 1e-14) {
        return r;
    }
    let (xy, xz, yz) = (r.xy, r.xz, r.yz);
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        let mut candidate = r;
        candidate.xy = mid * xy;
        candidate.xz = mid * xz;
        candidate.yz = mid * yz;
        if reynolds_stress_is_realizable(candidate, 1e-14) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    r.xy = lo * xy;
    r.xz = lo * xz;
    r.yz = lo * yz;
    r
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ReynoldsStressConfig {
    pub transport: RansTransportControls,
    pub density: f64,
    pub molecular_viscosity: f64,
    pub c_mu: f64,
    pub pressure_strain_slow: f64,
    pub pressure_strain_rapid: f64,
    pub c_epsilon1: f64,
    pub c_epsilon2: f64,
    pub sigma_r: f64,
    pub sigma_epsilon: f64,
    pub minimum_k: f64,
    pub minimum_epsilon: f64,
    pub minimum_normal_stress: f64,
    pub maximum_eddy_viscosity_ratio: f64,
}
impl Default for ReynoldsStressConfig {
    fn default() -> Self {
        Self {
            transport: Default::default(),
            density: 1.0,
            molecular_viscosity: 1.5e-5,
            c_mu: 0.09,
            pressure_strain_slow: 1.8,
            pressure_strain_rapid: 0.6,
            c_epsilon1: 1.44,
            c_epsilon2: 1.92,
            sigma_r: 0.82,
            sigma_epsilon: 1.3,
            minimum_k: 1e-12,
            minimum_epsilon: 1e-12,
            minimum_normal_stress: 1e-14,
            maximum_eddy_viscosity_ratio: 1e6,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct ReynoldsStressLinearResult {
    pub stress: [IterativeSolverResult; 6],
    pub epsilon: IterativeSolverResult,
}
impl From<[IterativeSolverResult; 7]> for ReynoldsStressLinearResult {
    fn from(r: [IterativeSolverResult; 7]) -> Self {
        Self {
            stress: std::array::from_fn(|i| r[i]),
            epsilon: r[6],
        }
    }
}
pub type ReynoldsStressSource = dyn Fn(Vec3, f64) -> Result<SymmetricTensor3> + Send + Sync;
pub struct ReynoldsStressTransport {
    state: RansState<7>,
    config: ReynoldsStressConfig,
    gradient: Vec<VelocityGradient3>,
    additional_source: Option<Box<ReynoldsStressSource>>,
}
impl ReynoldsStressTransport {
    pub fn new(mesh: PolyMesh, c: ReynoldsStressConfig) -> Result<Self> {
        check_constants(
            &[
                c.density,
                c.molecular_viscosity,
                c.c_mu,
                c.c_epsilon1,
                c.c_epsilon2,
                c.sigma_r,
                c.sigma_epsilon,
                c.minimum_k,
                c.minimum_epsilon,
                c.minimum_normal_stress,
                c.maximum_eddy_viscosity_ratio,
            ],
            &[c.pressure_strain_slow, c.pressure_strain_rapid],
        )?;
        let n = mesh.cell_count();
        let normal = 2.0 * c.minimum_k / 3.0;
        Ok(Self {
            state: RansState::new(
                mesh,
                c.transport,
                [normal, normal, normal, 0.0, 0.0, 0.0, c.minimum_epsilon],
            )?,
            config: c,
            gradient: vec![[0.0; 9]; n],
            additional_source: None,
        })
    }
    state_accessors!();
    pub fn linear_result(&self) -> Option<ReynoldsStressLinearResult> {
        self.state.last.map(Into::into)
    }
    pub fn epsilon(&self) -> &[f64] {
        &self.state.fields[6]
    }
    pub fn reynolds_stress(&self) -> Vec<SymmetricTensor3> {
        (0..self.mesh().cell_count())
            .map(|i| {
                SymmetricTensor3::from_components(std::array::from_fn(|q| self.state.fields[q][i]))
            })
            .collect()
    }
    pub fn turbulent_kinetic_energy(&self) -> Vec<f64> {
        self.reynolds_stress()
            .into_iter()
            .map(turbulent_kinetic_energy)
            .collect()
    }
    pub fn initialize_isotropic(&mut self, k: f64, epsilon: f64) -> Result<()> {
        require(
            k.is_finite()
                && epsilon.is_finite()
                && k >= self.config.minimum_k
                && epsilon >= self.config.minimum_epsilon,
            "RSM initial values",
        )?;
        let n = self.mesh().cell_count();
        let value = 2.0 * k / 3.0;
        self.initialize_profile(
            &vec![
                SymmetricTensor3 {
                    xx: value,
                    yy: value,
                    zz: value,
                    ..Default::default()
                };
                n
            ],
            &vec![epsilon; n],
        )
    }
    pub fn initialize_profile(
        &mut self,
        stress: &[SymmetricTensor3],
        epsilon: &[f64],
    ) -> Result<()> {
        let n = self.mesh().cell_count();
        require(
            stress.len() == n
                && stress.iter().all(|r| {
                    reynolds_stress_is_realizable(*r, 1e-12)
                        && turbulent_kinetic_energy(*r) >= self.config.minimum_k
                }),
            "RSM initial stresses are not realizable",
        )?;
        let mut eps = checked_field(epsilon, n, f64::MIN_POSITIVE)?;
        for e in &mut eps {
            *e = e.max(self.config.minimum_epsilon);
        }
        let mut fields: [Vec<f64>; 7] = std::array::from_fn(|q| {
            if q < 6 {
                stress.iter().map(|r| r.components()[q]).collect()
            } else {
                vec![]
            }
        });
        fields[6] = eps;
        self.state.initialize(fields);
        Ok(())
    }
    pub fn set_boundary(
        &mut self,
        patch: &str,
        kind: ScalarBoundaryType,
        value: SymmetricTensor3,
    ) -> Result<()> {
        require(
            kind != ScalarBoundaryType::FixedValue || reynolds_stress_is_realizable(value, 1e-12),
            "RSM boundary is not realizable",
        )?;
        finite_slice(&value.components(), "RSM boundary")?;
        let p = self.mesh().patch_index(patch)?;
        for q in 0..6 {
            self.state.boundary[q][p] = ScalarBoundaryCondition {
                boundary_type: kind,
                value: value.components()[q],
            };
        }
        Ok(())
    }
    pub fn set_epsilon_boundary(
        &mut self,
        patch: &str,
        kind: ScalarBoundaryType,
        value: f64,
    ) -> Result<()> {
        self.state
            .set_boundary(6, patch, kind, value, self.config.minimum_epsilon)
    }
    pub fn set_velocity(&mut self, u: &[Vec3]) -> Result<()> {
        let g = velocity_gradients(self.mesh(), u)?;
        self.gradient = g
            .iter()
            .map(|g| std::array::from_fn(|q| g[q / 3][q % 3]))
            .collect();
        Ok(())
    }
    pub fn set_velocity_gradient(&mut self, g: &[VelocityGradient3]) -> Result<()> {
        require(
            g.len() == self.mesh().cell_count() && g.iter().flatten().all(|v| v.is_finite()),
            "RSM velocity gradient",
        )?;
        self.gradient = g.to_vec();
        Ok(())
    }
    pub fn set_additional_source(&mut self, source: Option<Box<ReynoldsStressSource>>) {
        self.additional_source = source;
    }
    pub fn kinematic_eddy_viscosity(&self) -> Result<Vec<f64>> {
        let c = self.config;
        let v: Vec<_> = self
            .turbulent_kinetic_energy()
            .iter()
            .zip(self.epsilon())
            .map(|(&k, &e)| {
                (c.c_mu * k * k / e.max(c.minimum_epsilon))
                    .min(c.maximum_eddy_viscosity_ratio * c.molecular_viscosity)
            })
            .collect();
        finite_slice(&v, "RSM eddy viscosity")?;
        Ok(v)
    }
    pub fn step(&mut self) -> Result<ReynoldsStressLinearResult> {
        let c = self.config;
        let n = self.mesh().cell_count();
        let stresses = self.reynolds_stress();
        let k = self.turbulent_kinetic_energy();
        let nut = self.kinematic_eddy_viscosity()?;
        let (mut gamma, mut source, mut sink): ([Vec<f64>; 7], [Vec<f64>; 7], [Vec<f64>; 7]) = (
            std::array::from_fn(|_| vec![0.0; n]),
            std::array::from_fn(|_| vec![0.0; n]),
            std::array::from_fn(|_| vec![0.0; n]),
        );
        for cell in 0..n {
            let r = stresses[cell];
            let rm = [[r.xx, r.xy, r.xz], [r.xy, r.yy, r.yz], [r.xz, r.yz, r.zz]];
            let g = self.gradient[cell];
            let mut p = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    for a in 0..3 {
                        p[i][j] -= rm[i][a] * g[j * 3 + a] + rm[j][a] * g[i * 3 + a];
                    }
                }
            }
            let pc = [
                p[0][0],
                p[1][1],
                p[2][2],
                0.5 * (p[0][1] + p[1][0]),
                0.5 * (p[0][2] + p[2][0]),
                0.5 * (p[1][2] + p[2][1]),
            ];
            let ki = k[cell].max(c.minimum_k);
            let eps = self.epsilon()[cell];
            let slow = c.pressure_strain_slow * eps / ki;
            let pthird = (pc[0] + pc[1] + pc[2]) / 3.0;
            let rc = r.components();
            let extra = if let Some(f) = &self.additional_source {
                f(self.mesh().cells()[cell].center, self.time())?.components()
            } else {
                [0.0; 6]
            };
            for q in 0..6 {
                source[q][cell] = pc[q]
                    - slow * (rc[q] - if q < 3 { 2.0 * ki / 3.0 } else { 0.0 })
                    - c.pressure_strain_rapid * (pc[q] - if q < 3 { pthird } else { 0.0 })
                    - if q < 3 { 2.0 * eps / 3.0 } else { 0.0 };
                source[q][cell] += extra[q];
                gamma[q][cell] = c.molecular_viscosity + nut[cell] / c.sigma_r;
            }
            gamma[6][cell] = c.molecular_viscosity + nut[cell] / c.sigma_epsilon;
            let ep = eps.max(c.minimum_epsilon);
            let pk = (0.5 * (pc[0] + pc[1] + pc[2])).max(0.0);
            source[6][cell] = c.c_epsilon1 * (ep / ki) * pk;
            sink[6][cell] = c.c_epsilon2 * ep / ki;
        }
        let result = self.state.advance(
            gamma,
            source,
            sink,
            [
                c.minimum_normal_stress,
                c.minimum_normal_stress,
                c.minimum_normal_stress,
                -f64::MAX,
                -f64::MAX,
                -f64::MAX,
                c.minimum_epsilon,
            ],
        )?;
        for i in 0..n {
            let r =
                SymmetricTensor3::from_components(std::array::from_fn(|q| self.state.fields[q][i]));
            let v = project_realizable(r, c.minimum_normal_stress).components();
            for q in 0..6 {
                self.state.fields[q][i] = v[q];
            }
        }
        Ok(result.into())
    }
    pub fn run(&mut self, steps: usize) -> Result<Option<ReynoldsStressLinearResult>> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(self.linear_result())
    }
}
