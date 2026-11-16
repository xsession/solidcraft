use super::operators::{linear_scalar, linear_vector};
use super::*;
use crate::{
    core::{
        preconditioned_conjugate_gradient, restarted_gmres, rms, CsrBuilder,
        IdentityPreconditioner, Ilu0Preconditioner, IterativeSolverResult, JacobiPreconditioner,
        LinearConfig, Preconditioner, Vec3,
    },
    finite_slice, positive, require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PressurePreconditionerKind {
    None,
    Jacobi,
    ExternalAmg,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CollocatedIncompressibleConfig {
    pub density: f64,
    pub kinematic_viscosity: f64,
    pub dt: f64,
    pub momentum_sweeps: usize,
    pub use_krylov_momentum: bool,
    pub momentum_iterations: usize,
    pub momentum_tolerance: f64,
    pub pressure_iterations: usize,
    pub pressure_tolerance: f64,
    pub nonorthogonal_correctors: usize,
    pub pressure_correctors: usize,
    pub outer_correctors: usize,
    pub velocity_relaxation: f64,
    pub pressure_relaxation: f64,
    pub convection_scheme: FaceInterpolationScheme,
    pub include_convection: bool,
    pub temporal_scheme: TemporalScheme,
    pub crank_nicolson_off_centering: f64,
    pub pressure_preconditioner: PressurePreconditionerKind,
}
impl Default for CollocatedIncompressibleConfig {
    fn default() -> Self {
        Self {
            density: 1.0,
            kinematic_viscosity: 1e-2,
            dt: 1e-2,
            momentum_sweeps: 3,
            use_krylov_momentum: true,
            momentum_iterations: 300,
            momentum_tolerance: 1e-10,
            pressure_iterations: 1000,
            pressure_tolerance: 1e-10,
            nonorthogonal_correctors: 2,
            pressure_correctors: 2,
            outer_correctors: 2,
            velocity_relaxation: 0.7,
            pressure_relaxation: 0.3,
            convection_scheme: FaceInterpolationScheme::Upwind,
            include_convection: true,
            temporal_scheme: TemporalScheme::Euler,
            crank_nicolson_off_centering: 1.0,
            pressure_preconditioner: PressurePreconditionerKind::None,
        }
    }
}
impl CollocatedIncompressibleConfig {
    fn validate(self) -> Result<()> {
        positive(self.density, "density")?;
        positive(self.dt, "dt")?;
        positive(self.momentum_tolerance, "momentum tolerance")?;
        positive(self.pressure_tolerance, "pressure tolerance")?;
        require(
            self.kinematic_viscosity.is_finite()
                && self.kinematic_viscosity >= 0.0
                && self.momentum_sweeps > 0
                && self.momentum_iterations > 0
                && self.pressure_iterations > 0
                && self.pressure_correctors > 0
                && self.outer_correctors > 0
                && self.nonorthogonal_correctors < usize::MAX
                && self.velocity_relaxation > 0.0
                && self.velocity_relaxation <= 1.0
                && self.pressure_relaxation > 0.0
                && self.pressure_relaxation <= 1.0,
            "invalid collocated controls",
        )?;
        if self.include_convection && self.convection_scheme != FaceInterpolationScheme::Upwind {
            return Err(SolverError::Unsupported("the upstream collocated momentum operator implements upwind only; higher-order face schemes are available as separate operators".into()));
        }
        if self.temporal_scheme == TemporalScheme::CrankNicolson {
            crank_nicolson_implicit_weight(self.crank_nicolson_off_centering)?;
        }
        Ok(())
    }
}
#[derive(Clone, Copy, Debug, Serialize)]
pub struct CollocatedIterationInfo {
    pub velocity_rms_change: f64,
    pub continuity_l2: f64,
    pub continuity_max: f64,
    pub pressure: IterativeSolverResult,
    pub pressure_solves: usize,
    pub steady_converged: Option<bool>,
}
type AmgCycle<'a> = dyn Fn(&[f64], &mut [f64]) -> Result<()> + 'a;
struct CyclePreconditioner<'a> {
    cycle: &'a AmgCycle<'a>,
}
impl Preconditioner for CyclePreconditioner<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<()> {
        (self.cycle)(r, z)
    }
}
pub struct CollocatedIncompressible<'a> {
    mesh: PolyMesh,
    config: CollocatedIncompressibleConfig,
    velocity_boundary: Vec<VelocityBoundaryCondition>,
    pressure_boundary: Vec<PressureBoundaryCondition>,
    velocity: Vec<Vec3>,
    old_velocity: Vec<Vec3>,
    h_by_a: Vec<Vec3>,
    pressure_mobility: Vec<f64>,
    pressure: Vec<f64>,
    face_flux: Vec<f64>,
    pressure_face_coefficient: Vec<f64>,
    body_acceleration: Vec<Vec3>,
    implicit_drag: Vec<f64>,
    amg: Option<Box<AmgCycle<'a>>>,
    momentum_results: [IterativeSolverResult; 3],
    time: f64,
    steps: usize,
}
fn empty_result() -> IterativeSolverResult {
    IterativeSolverResult {
        iterations: 0,
        residual_rms: 0.0,
        converged: false,
    }
}
impl<'a> CollocatedIncompressible<'a> {
    pub fn new(mesh: PolyMesh, config: CollocatedIncompressibleConfig) -> Result<Self> {
        config.validate()?;
        let n = mesh.cell_count();
        let nf = mesh.face_count();
        let np = mesh.patches().len();
        let mut s = Self {
            mesh,
            config,
            velocity_boundary: vec![VelocityBoundaryCondition::default(); np],
            pressure_boundary: vec![PressureBoundaryCondition::default(); np],
            velocity: vec![Vec3::ZERO; n],
            old_velocity: vec![Vec3::ZERO; n],
            h_by_a: vec![Vec3::ZERO; n],
            pressure_mobility: vec![config.dt; n],
            pressure: vec![0.0; n],
            face_flux: vec![0.0; nf],
            pressure_face_coefficient: vec![0.0; nf],
            body_acceleration: vec![Vec3::ZERO; n],
            implicit_drag: vec![0.0; n],
            amg: None,
            momentum_results: [empty_result(); 3],
            time: 0.0,
            steps: 0,
        };
        s.rebuild_flux()?;
        Ok(s)
    }
    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn config(&self) -> CollocatedIncompressibleConfig {
        self.config
    }
    pub fn velocity(&self) -> &[Vec3] {
        &self.velocity
    }
    pub fn pressure(&self) -> &[f64] {
        &self.pressure
    }
    pub fn face_flux(&self) -> &[f64] {
        &self.face_flux
    }
    pub fn h_by_a(&self) -> &[Vec3] {
        &self.h_by_a
    }
    pub fn pressure_mobility(&self) -> &[f64] {
        &self.pressure_mobility
    }
    pub fn momentum_results(&self) -> &[IterativeSolverResult; 3] {
        &self.momentum_results
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn steps(&self) -> usize {
        self.steps
    }
    pub fn set_velocity_boundary(
        &mut self,
        patch: &str,
        boundary_type: VelocityBoundaryType,
        value: Vec3,
    ) -> Result<()> {
        require(value.is_finite(), "nonfinite velocity boundary")?;
        let p = self.mesh.patch_index(patch)?;
        self.velocity_boundary[p] = VelocityBoundaryCondition {
            boundary_type,
            value,
        };
        self.rebuild_flux()
    }
    /// Change the value of an existing fixed-value velocity patch between steps.
    /// Interior Rhie-Chow fluxes are retained; only that patch's boundary fluxes change,
    /// so the next pressure correction restores continuity with the new inflow.
    pub fn set_velocity_boundary_value(&mut self, patch: &str, value: Vec3) -> Result<()> {
        require(value.is_finite(), "nonfinite velocity boundary")?;
        let p = self.mesh.patch_index(patch)?;
        let condition = self
            .velocity_boundary
            .get_mut(p)
            .ok_or_else(|| crate::SolverError::InvalidInput("velocity patch index".into()))?;
        require(
            condition.boundary_type == VelocityBoundaryType::FixedValue,
            "only fixed-value velocity patches can be updated in place",
        )?;
        condition.value = value;
        for (fi, face) in self.mesh.faces().iter().enumerate() {
            if face.patch == Some(p) {
                if let Some(flux) = self.face_flux.get_mut(fi) {
                    *flux = value.dot(face.area);
                }
            }
        }
        Ok(())
    }
    pub fn set_pressure_boundary(
        &mut self,
        patch: &str,
        boundary_type: PressureBoundaryType,
        value: f64,
    ) -> Result<()> {
        require(value.is_finite(), "nonfinite pressure boundary")?;
        let p = self.mesh.patch_index(patch)?;
        self.pressure_boundary[p] = PressureBoundaryCondition {
            boundary_type,
            value,
        };
        Ok(())
    }
    pub fn initialize_uniform(&mut self, u: Vec3, p: f64) -> Result<()> {
        self.initialize_fields(|_| u, |_| p)
    }
    pub fn initialize_fields(
        &mut self,
        u: impl Fn(Vec3) -> Vec3,
        p: impl Fn(Vec3) -> f64,
    ) -> Result<()> {
        let velocity: Vec<_> = self.mesh.cells().iter().map(|c| u(c.center)).collect();
        let pressure: Vec<_> = self.mesh.cells().iter().map(|c| p(c.center)).collect();
        require(
            velocity.iter().all(|v| v.is_finite()),
            "nonfinite initial velocity",
        )?;
        finite_slice(&pressure, "initial pressure")?;
        self.velocity = velocity;
        self.old_velocity = self.velocity.clone();
        self.h_by_a = self.velocity.clone();
        self.pressure = pressure;
        self.time = 0.0;
        self.steps = 0;
        self.rebuild_flux()
    }
    pub fn set_pressure_amg_cycle(
        &mut self,
        cycle: impl Fn(&[f64], &mut [f64]) -> Result<()> + 'a,
    ) {
        self.amg = Some(Box::new(cycle));
        self.config.pressure_preconditioner = PressurePreconditionerKind::ExternalAmg;
    }
    pub fn set_body_acceleration(&mut self, a: &[Vec3]) -> Result<()> {
        require(
            a.len() == self.mesh.cell_count() && a.iter().all(|v| v.is_finite()),
            "body acceleration size/values",
        )?;
        self.body_acceleration.copy_from_slice(a);
        Ok(())
    }
    /// Per-cell implicit linear drag coefficient beta (1/s): the momentum sink is
    /// `-beta * u` per unit mass, assembled on the momentum diagonal so it also enters
    /// the Rhie-Chow pressure mobility. Used for porous resistance zones.
    pub fn set_implicit_drag(&mut self, beta: &[f64]) -> Result<()> {
        require(beta.len() == self.mesh.cell_count(), "implicit drag size")?;
        require(
            beta.iter().all(|b| b.is_finite() && *b >= 0.0),
            "implicit drag must be finite and non-negative",
        )?;
        self.implicit_drag.clear();
        self.implicit_drag.extend_from_slice(beta);
        Ok(())
    }
    pub fn implicit_drag(&self) -> &[f64] {
        &self.implicit_drag
    }
    pub fn clear_body_acceleration(&mut self) {
        self.body_acceleration.fill(Vec3::ZERO);
    }
    fn fixed_pressure(&self) -> bool {
        self.pressure_boundary
            .iter()
            .any(|b| b.boundary_type == PressureBoundaryType::FixedValue)
    }
    fn rebuild_flux(&mut self) -> Result<()> {
        let b = boundary_velocity_values(&self.mesh, &self.velocity, &self.velocity_boundary)?;
        self.face_flux = face_flux_from_velocity(&self.mesh, &self.velocity, &b)?;
        Ok(())
    }
    fn momentum_predictor(&mut self, time_source: &[Vec3], physical: bool) -> Result<()> {
        let m = &self.mesh;
        let n = m.cell_count();
        require(time_source.len() == n, "momentum time source size")?;
        let pb = boundary_pressure_values(m, &self.pressure, &self.pressure_boundary)?;
        let gp = gauss_gradient_scalar(m, &self.pressure, &pb)?;
        let scheme = if physical {
            self.config.temporal_scheme
        } else {
            TemporalScheme::Euler
        };
        let bdf = scheme == TemporalScheme::BackwardBdf2 && self.steps > 0;
        let alpha0 = if bdf { 1.5 } else { 1.0 };
        let theta = if scheme == TemporalScheme::CrankNicolson {
            crank_nicolson_implicit_weight(self.config.crank_nicolson_off_centering)?
        } else {
            1.0
        };
        let mut spatial = CsrBuilder::new(n, n);
        let mut boundary_rhs = vec![Vec3::ZERO; n];
        for c in 0..n {
            let mut ap = 0.0;
            for &fi in &m.cell_faces()[c] {
                let f = m.faces()[fi];
                let phi = if f.owner == c {
                    self.face_flux[fi]
                } else {
                    -self.face_flux[fi]
                };
                let diff =
                    self.config.kinematic_viscosity * decompose_face_area(m, fi)?.orthogonal_metric;
                if let Some(neighbour) = f.neighbour {
                    let other = if f.owner == c { neighbour } else { f.owner };
                    let enter = if self.config.include_convection {
                        (-phi).max(0.0)
                    } else {
                        0.0
                    };
                    let leave = if self.config.include_convection {
                        phi.max(0.0)
                    } else {
                        0.0
                    };
                    ap += diff + leave;
                    let coupling = diff + enter;
                    if coupling != 0.0 {
                        spatial.add(c, other, -coupling)?;
                    }
                } else {
                    let bc = self.velocity_boundary[f.patch.ok_or_else(|| {
                        crate::SolverError::InvalidInput("boundary face without patch".into())
                    })?];
                    if bc.boundary_type == VelocityBoundaryType::FixedValue {
                        let enter = if self.config.include_convection {
                            (-phi).max(0.0)
                        } else {
                            0.0
                        };
                        let leave = if self.config.include_convection {
                            phi.max(0.0)
                        } else {
                            0.0
                        };
                        ap += diff + leave;
                        boundary_rhs[c] += bc.value * (diff + enter);
                    } else if self.config.include_convection {
                        ap += phi;
                    }
                }
            }
            if ap != 0.0 {
                spatial.add(c, c, ap)?;
            }
        }
        let spatial = spatial.build(0.0)?;
        let mut builder = CsrBuilder::new(n, n);
        for row in 0..n {
            for k in spatial.row_offsets()[row]..spatial.row_offsets()[row + 1] {
                builder.add(row, spatial.columns()[k], theta * spatial.values()[k])?;
            }
            let drag = self.implicit_drag.get(row).copied().unwrap_or(0.0);
            // Without drag keep the original expression, so existing cases stay bit-identical.
            let temporal = if drag == 0.0 {
                alpha0 * m.cells()[row].volume / self.config.dt
            } else {
                (alpha0 / self.config.dt + drag) * m.cells()[row].volume
            };
            builder.add(row, row, temporal)?;
        }
        let a = builder.build(0.0)?;
        let diagonal = a.diagonal();
        require(
            diagonal.iter().all(|d| d.is_finite() && *d > 0.0),
            "nonpositive momentum diagonal",
        )?;
        let mut components: [Vec<f64>; 3] =
            std::array::from_fn(|axis| self.velocity.iter().map(|u| u.component(axis)).collect());
        let mut rhs: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; n]);
        let mut old_action: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; n]);
        if theta < 1.0 {
            for axis in 0..3 {
                let component: Vec<_> = time_source.iter().map(|u| u.component(axis)).collect();
                old_action[axis] = spatial.apply(&component)?;
            }
        }
        for c in 0..n {
            let vol = m.cells()[c].volume;
            let temporal = if bdf {
                time_source[c] * 2.0 - self.old_velocity[c] * 0.5
            } else {
                time_source[c]
            };
            let full = temporal * (vol / self.config.dt) + boundary_rhs[c] - gp[c] * vol
                + self.body_acceleration[c] * vol;
            for axis in 0..3 {
                rhs[axis][c] = full.component(axis) - (1.0 - theta) * old_action[axis][c];
            }
        }
        if self.config.use_krylov_momentum {
            let ilu = Ilu0Preconditioner::new(&a, 1e-12)?;
            for axis in 0..3 {
                let r = restarted_gmres(
                    &rhs[axis],
                    &mut components[axis],
                    |v, y| a.multiply(v, y),
                    &ilu,
                    LinearConfig {
                        max_iterations: self.config.momentum_iterations,
                        restart: 30,
                        relative_tolerance: self.config.momentum_tolerance,
                    },
                )?;
                self.momentum_results[axis] = r;
                if !r.converged {
                    return Err(SolverError::NumericalFailure(format!(
                        "collocated momentum GMRES did not converge (axis {axis}, {} iterations, residual {:e}, rhs rms {:e}, diag min {:e} max {:e})",
                        r.iterations, r.residual_rms, rms(&rhs[axis]),
                        diagonal.iter().cloned().fold(f64::MAX, f64::min), diagonal.iter().cloned().fold(f64::MIN, f64::max)
                    )));
                }
            }
        } else {
            for axis in 0..3 {
                for _ in 0..self.config.momentum_sweeps {
                    let mut next = vec![0.0; n];
                    for row in 0..n {
                        let mut sum = rhs[axis][row];
                        for k in a.row_offsets()[row]..a.row_offsets()[row + 1] {
                            let col = a.columns()[k];
                            if col != row {
                                sum -= a.values()[k] * components[axis][col];
                            }
                        }
                        next[row] = sum / diagonal[row];
                    }
                    components[axis] = next;
                }
                let action = a.apply(&components[axis])?;
                let residual: Vec<_> = rhs[axis].iter().zip(action).map(|(b, a)| b - a).collect();
                let residual_rms = rms(&residual);
                self.momentum_results[axis] = IterativeSolverResult {
                    iterations: self.config.momentum_sweeps,
                    residual_rms,
                    converged: residual_rms
                        <= self.config.momentum_tolerance * rms(&rhs[axis]).max(1.0),
                };
            }
        }
        for c in 0..n {
            self.velocity[c] = Vec3::new(components[0][c], components[1][c], components[2][c]);
            self.pressure_mobility[c] = m.cells()[c].volume / diagonal[c];
            self.h_by_a[c] = self.velocity[c] + gp[c] * self.pressure_mobility[c];
        }
        Ok(())
    }
    fn nonorthogonal_flux(&self, p: &[f64]) -> Result<Vec<f64>> {
        let m = &self.mesh;
        let pb = boundary_pressure_values(m, p, &self.pressure_boundary)?;
        let g = gauss_gradient_scalar(m, p, &pb)?;
        let mut out = vec![0.0; m.face_count()];
        for (fi, r) in out.iter_mut().enumerate() {
            let f = m.faces()[fi];
            if f.patch.is_some_and(|p| {
                self.pressure_boundary[p].boundary_type == PressureBoundaryType::ZeroGradient
            }) {
                continue;
            }
            let geom = decompose_face_area(m, fi)?;
            let mobility = linear_scalar(m, f, &self.pressure_mobility)?;
            let gf = if f.boundary() {
                g[f.owner]
            } else {
                linear_vector(m, f, &g)?
            };
            *r = mobility * gf.dot(geom.nonorthogonal_area);
        }
        Ok(out)
    }
    fn pressure_action(&self, x: &[f64], out: &mut [f64]) -> Result<()> {
        require(
            x.len() == self.mesh.cell_count() && out.len() == x.len(),
            "pressure action size",
        )?;
        for (c, r) in out.iter_mut().enumerate() {
            let mut sum = 0.0;
            for &fi in &self.mesh.cell_faces()[c] {
                let f = self.mesh.faces()[fi];
                let coeff = self.pressure_face_coefficient[fi];
                if coeff <= 0.0 {
                    continue;
                }
                if let Some(n) = f.neighbour {
                    let other = if f.owner == c { n } else { f.owner };
                    sum += coeff * (x[c] - x[other]);
                } else {
                    sum += coeff * x[c];
                }
            }
            *r = sum;
        }
        Ok(())
    }
    fn correct_pressure(&mut self, relaxation: f64) -> Result<(IterativeSolverResult, usize)> {
        let m = &self.mesh;
        for (fi, coeff) in self.pressure_face_coefficient.iter_mut().enumerate() {
            let f = m.faces()[fi];
            *coeff = if f.patch.is_some_and(|p| {
                self.pressure_boundary[p].boundary_type == PressureBoundaryType::ZeroGradient
            }) {
                0.0
            } else {
                linear_scalar(m, f, &self.pressure_mobility)?
                    * decompose_face_area(m, fi)?.orthogonal_metric
            };
        }
        let hb = boundary_velocity_values(m, &self.h_by_a, &self.velocity_boundary)?;
        let phi_h = face_flux_from_velocity(m, &self.h_by_a, &hb)?;
        let mut p = self.pressure.clone();
        let fixed = self.fixed_pressure();
        let mut result = empty_result();
        let mut solves = 0;
        for _ in 0..self.config.nonorthogonal_correctors + 1 {
            let correction = self.nonorthogonal_flux(&p)?;
            let pb = boundary_pressure_values(m, &p, &self.pressure_boundary)?;
            let mut rhs = vec![0.0; m.cell_count()];
            for (c, r) in rhs.iter_mut().enumerate() {
                for &fi in &m.cell_faces()[c] {
                    let f = m.faces()[fi];
                    let sign = if f.owner == c { 1.0 } else { -1.0 };
                    *r += sign * (-phi_h[fi] + correction[fi]);
                    if f.patch.is_some_and(|b| {
                        self.pressure_boundary[b].boundary_type == PressureBoundaryType::FixedValue
                    }) {
                        *r += self.pressure_face_coefficient[fi] * pb[fi];
                    }
                }
            }
            if !fixed {
                let mean = rhs.iter().sum::<f64>() / rhs.len() as f64;
                let pmean = p.iter().sum::<f64>() / p.len() as f64;
                for i in 0..p.len() {
                    rhs[i] -= mean;
                    p[i] -= pmean;
                }
            }
            let controls = LinearConfig {
                max_iterations: self.config.pressure_iterations,
                restart: 30,
                relative_tolerance: self.config.pressure_tolerance,
            };
            result = match self.config.pressure_preconditioner {
                PressurePreconditionerKind::None => preconditioned_conjugate_gradient(
                    &rhs,
                    &mut p,
                    |x, out| self.pressure_action(x, out),
                    &IdentityPreconditioner,
                    controls,
                )?,
                PressurePreconditionerKind::Jacobi => {
                    let d: Vec<_> = (0..m.cell_count())
                        .map(|c| {
                            m.cell_faces()[c]
                                .iter()
                                .map(|&fi| self.pressure_face_coefficient[fi])
                                .sum()
                        })
                        .collect();
                    let jacobi = JacobiPreconditioner::new(&d, 1e-30)?;
                    preconditioned_conjugate_gradient(
                        &rhs,
                        &mut p,
                        |x, out| self.pressure_action(x, out),
                        &jacobi,
                        controls,
                    )?
                }
                PressurePreconditionerKind::ExternalAmg => {
                    let cycle = self.amg.as_ref().ok_or_else(|| {
                        SolverError::InvalidInput("external AMG cycle was not configured".into())
                    })?;
                    let pre = CyclePreconditioner {
                        cycle: cycle.as_ref(),
                    };
                    preconditioned_conjugate_gradient(
                        &rhs,
                        &mut p,
                        |x, out| self.pressure_action(x, out),
                        &pre,
                        controls,
                    )?
                }
            };
            solves += 1;
            if !result.converged {
                break;
            }
        }
        for i in 0..p.len() {
            self.pressure[i] += relaxation * (p[i] - self.pressure[i]);
        }
        if !fixed {
            let mean = self.pressure.iter().sum::<f64>() / self.pressure.len() as f64;
            for v in &mut self.pressure {
                *v -= mean;
            }
        }
        let pb = boundary_pressure_values(m, &self.pressure, &self.pressure_boundary)?;
        let gp = gauss_gradient_scalar(m, &self.pressure, &pb)?;
        for c in 0..m.cell_count() {
            let corrected = self.h_by_a[c] - gp[c] * self.pressure_mobility[c];
            let previous = self.velocity[c];
            self.velocity[c] = previous
                + (corrected - previous)
                    * (if relaxation < 1.0 {
                        self.config.velocity_relaxation
                    } else {
                        1.0
                    });
        }
        self.face_flux = rhie_chow_face_flux(
            m,
            &self.h_by_a,
            &self.pressure_mobility,
            &self.pressure,
            &self.velocity_boundary,
            &self.pressure_boundary,
            true,
        )?;
        Ok((result, solves))
    }
    fn coupled_sequence(
        &mut self,
        time_source: &[Vec3],
        correctors: usize,
        relaxation: f64,
        physical: bool,
    ) -> Result<CollocatedIterationInfo> {
        let before = self.velocity.clone();
        self.momentum_predictor(time_source, physical)?;
        let mut pressure = empty_result();
        let mut solves = 0;
        for _ in 0..correctors {
            let (r, n) = self.correct_pressure(relaxation)?;
            pressure = r;
            solves += n;
            if !pressure.converged {
                break;
            }
        }
        let velocity_rms_change = ((0..before.len())
            .map(|i| (self.velocity[i] - before[i]).norm_squared())
            .sum::<f64>()
            / before.len() as f64)
            .sqrt();
        Ok(CollocatedIterationInfo {
            velocity_rms_change,
            continuity_l2: self.continuity_l2()?,
            continuity_max: self.continuity_max()?,
            pressure,
            pressure_solves: solves,
            steady_converged: None,
        })
    }
    pub fn iterate_simple(&mut self) -> Result<CollocatedIterationInfo> {
        self.coupled_sequence(
            &self.velocity.clone(),
            1,
            self.config.pressure_relaxation,
            false,
        )
    }
    pub fn solve_simple(
        &mut self,
        max_iterations: usize,
        velocity_tolerance: f64,
        continuity_tolerance: f64,
    ) -> Result<CollocatedIterationInfo> {
        require(max_iterations > 0, "SIMPLE iterations must be positive")?;
        positive(velocity_tolerance, "velocity tolerance")?;
        positive(continuity_tolerance, "continuity tolerance")?;
        let mut last = None;
        for _ in 0..max_iterations {
            let mut info = self.iterate_simple()?;
            let converged = info.pressure.converged
                && self.momentum_results.iter().all(|r| r.converged)
                && info.velocity_rms_change <= velocity_tolerance
                && info.continuity_l2 <= continuity_tolerance;
            info.steady_converged = Some(converged);
            if converged || !info.pressure.converged {
                return Ok(info);
            }
            last = Some(info);
        }
        last.ok_or_else(|| SolverError::InvalidInput("SIMPLE iterations must be positive".into()))
    }
    pub fn step_piso(&mut self) -> Result<CollocatedIterationInfo> {
        let current = self.velocity.clone();
        let info = self.coupled_sequence(&current, self.config.pressure_correctors, 1.0, true)?;
        self.old_velocity = current;
        self.time += self.config.dt;
        self.steps += 1;
        Ok(info)
    }
    pub fn step_pimple(&mut self) -> Result<CollocatedIterationInfo> {
        let current = self.velocity.clone();
        let mut info = None;
        for _ in 0..self.config.outer_correctors {
            let r = self.coupled_sequence(&current, self.config.pressure_correctors, 1.0, true)?;
            let failed = !r.pressure.converged;
            info = Some(r);
            if failed {
                break;
            }
        }
        self.old_velocity = current;
        self.time += self.config.dt;
        self.steps += 1;
        info.ok_or_else(|| {
            SolverError::InvalidInput("PIMPLE needs at least one outer corrector".into())
        })
    }
    pub fn continuity_l2(&self) -> Result<f64> {
        flux_divergence_l2(&self.mesh, &self.face_flux)
    }
    pub fn continuity_max(&self) -> Result<f64> {
        flux_divergence_max(&self.mesh, &self.face_flux)
    }
    pub fn kinetic_energy(&self) -> f64 {
        self.config.density
            * self
                .velocity
                .iter()
                .zip(self.mesh.cells())
                .map(|(u, c)| 0.5 * u.dot(*u) * c.volume)
                .sum::<f64>()
    }
}
