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
    /// Rebuild H(u)/a_P from the corrected velocity before each further PISO corrector,
    /// and solve the pressure to a tolerance relative to the cell flux scale. Off by
    /// default for parity with the C++ cfd_solvers, whose correctors skip both, which
    /// limits the stable time step to diffusion numbers ν·Δt/h² of roughly 0.3.
    pub piso_neighbour_update: bool,
    /// Balanced-force pressure and body-force coupling: each cell's pressure-gradient-minus-
    /// body-force term is rebuilt by least squares from the face values (compact face pressure
    /// gradient minus face body force), and the face fluxes carry the face body force the same
    /// way. A force balanced by pressure then leaves no spurious velocity, also where the force
    /// jumps (a fan block, a porous edge, a buoyant plume's edge). Off by default for parity
    /// with the C++ cfd_solvers, which adds the cell body force to the cell Gauss gradient.
    #[serde(default)]
    pub balanced_body_force: bool,
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
            piso_neighbour_update: false,
            balanced_body_force: false,
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
    profiles: FaceProfiles,
    /// Subtracted from the face flux in momentum convection: the frame flux (Ω×r)·S of
    /// rotating regions (MRF), so momentum is carried by the flux relative to the frame.
    convective_offset: Vec<f64>,
    velocity: Vec<Vec3>,
    old_velocity: Vec<Vec3>,
    h_by_a: Vec<Vec3>,
    /// The last momentum operator and its source without the pressure gradient, kept
    /// so each further PISO corrector can rebuild H(u)/a_P from the corrected velocity.
    momentum_operator: Option<crate::core::CsrMatrix>,
    momentum_source: Vec<Vec3>,
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
/// Solve a 3 × 3 system by Cramer's rule (None when singular).
fn solve3(a: [[f64; 3]; 3], b: Vec3) -> Option<Vec3> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(a);
    let scale = a.iter().flatten().fold(0.0f64, |s, v| s.max(v.abs()));
    if !(d.is_finite() && d.abs() > 1e-12 * scale * scale * scale) {
        return None;
    }
    let rhs = [b.x, b.y, b.z];
    let column = |k: usize| {
        let mut m = a;
        for (row, value) in m.iter_mut().zip(rhs) {
            row[k] = value;
        }
        det(m) / d
    };
    Some(Vec3::new(column(0), column(1), column(2)))
}
/// Root-mean-square over cells of the summed face flux magnitudes; 1 without flow.
fn flux_scale(m: &PolyMesh, flux: &[f64]) -> f64 {
    let n = m.cell_count();
    if n == 0 {
        return 1.0;
    }
    let sum: f64 = (0..n)
        .map(|c| {
            let t: f64 = m.cell_faces()[c].iter().map(|&fi| flux[fi].abs()).sum();
            t * t
        })
        .sum();
    let s = (sum / n as f64).sqrt();
    if s.is_finite() && s > 0.0 { s } else { 1.0 }
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
            profiles: FaceProfiles::default(),
            convective_offset: Vec::new(),
            velocity: vec![Vec3::ZERO; n],
            old_velocity: vec![Vec3::ZERO; n],
            h_by_a: vec![Vec3::ZERO; n],
            momentum_operator: None,
            momentum_source: Vec::new(),
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
                    *flux = self.profiles.velocity(fi).unwrap_or(value).dot(face.area);
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
    /// Per-face boundary values on fixed-value patches (transferred boundary conditions).
    /// Set after the patch conditions; the face fluxes are rebuilt.
    pub fn set_face_profiles(&mut self, profiles: FaceProfiles) -> Result<()> {
        profiles.validate(&self.mesh)?;
        self.profiles = profiles;
        self.rebuild_flux()
    }
    /// Per-face flux removed from momentum convection (rotating regions); empty for none.
    pub fn set_convective_flux_offset(&mut self, offset: Vec<f64>) -> Result<()> {
        require(
            offset.is_empty() || offset.len() == self.mesh.face_count(),
            "convective flux offset size",
        )?;
        finite_slice(&offset, "convective flux offset")?;
        self.convective_offset = offset;
        Ok(())
    }
    /// The flux that convects transported quantities: the face flux, relative to the
    /// rotating frame inside rotating regions.
    pub fn convective_flux(&self) -> Vec<f64> {
        if self.convective_offset.is_empty() {
            return self.face_flux.clone();
        }
        self.face_flux
            .iter()
            .zip(&self.convective_offset)
            .map(|(phi, offset)| phi - offset)
            .collect()
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
        let b = boundary_velocity_values_with(&self.mesh, &self.velocity, &self.velocity_boundary, &self.profiles)?;
        self.face_flux = face_flux_from_velocity(&self.mesh, &self.velocity, &b)?;
        Ok(())
    }
    fn momentum_predictor(&mut self, time_source: &[Vec3], physical: bool) -> Result<()> {
        let m = &self.mesh;
        let n = m.cell_count();
        require(time_source.len() == n, "momentum time source size")?;
        let pb = boundary_pressure_values_with(m, &self.pressure, &self.pressure_boundary, &self.profiles)?;
        // With balanced body forces `gp` is the face-reconstructed (∇p − a); the body force
        // is then inside it and not added again below.
        let balanced = self.config.balanced_body_force;
        let gp = if balanced { self.balanced_residual(&self.pressure)? } else { gauss_gradient_scalar(m, &self.pressure, &pb)? };
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
                let relative = self.face_flux[fi] - self.convective_offset.get(fi).copied().unwrap_or(0.0);
                let phi = if f.owner == c { relative } else { -relative };
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
                        boundary_rhs[c] += self.profiles.velocity(fi).unwrap_or(bc.value) * (diff + enter);
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
            let body = if balanced { Vec3::ZERO } else { self.body_acceleration[c] * vol };
            let full = temporal * (vol / self.config.dt) + boundary_rhs[c] - gp[c] * vol + body;
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
        self.momentum_source = (0..n)
            .map(|c| Vec3::new(rhs[0][c], rhs[1][c], rhs[2][c]) + gp[c] * m.cells()[c].volume)
            .collect();
        for c in 0..n {
            self.velocity[c] = Vec3::new(components[0][c], components[1][c], components[2][c]);
            self.pressure_mobility[c] = m.cells()[c].volume / diagonal[c];
            self.h_by_a[c] = self.velocity[c] + gp[c] * self.pressure_mobility[c];
        }
        self.momentum_operator = Some(a);
        Ok(())
    }
    /// PISO: H(u)/a_P = (b - Σ_{N≠P} a_N u_N) / a_P with the corrected velocity, so the
    /// next pressure corrector sees the neighbour velocity corrections. Without this
    /// every corrector re-solves the same equation and a time step is one unrelaxed
    /// SIMPLE iteration, which diverges once diffusion dominates (a_P ≫ V/Δt).
    fn refresh_h_by_a(&mut self) -> Result<()> {
        let Some(a) = self.momentum_operator.as_ref() else {
            return Ok(());
        };
        let n = self.mesh.cell_count();
        require(self.momentum_source.len() == n, "momentum source size")?;
        for row in 0..n {
            let mut off = Vec3::ZERO;
            let mut diag = 0.0;
            for k in a.row_offsets()[row]..a.row_offsets()[row + 1] {
                let col = a.columns()[k];
                if col == row {
                    diag += a.values()[k];
                } else {
                    off += self.velocity[col] * a.values()[k];
                }
            }
            require(diag.is_finite() && diag > 0.0, "nonpositive momentum diagonal")?;
            self.h_by_a[row] = (self.momentum_source[row] - off) * (1.0 / diag);
        }
        Ok(())
    }
    fn nonorthogonal_flux(&self, p: &[f64]) -> Result<Vec<f64>> {
        let m = &self.mesh;
        let pb = boundary_pressure_values_with(m, p, &self.pressure_boundary, &self.profiles)?;
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
        let hb = boundary_velocity_values_with(m, &self.h_by_a, &self.velocity_boundary, &self.profiles)?;
        let mut phi_h = face_flux_from_velocity(m, &self.h_by_a, &hb)?;
        let body_flux = if self.config.balanced_body_force { Some(self.body_flux()?) } else { None };
        if let Some(extra) = &body_flux {
            for (phi, b) in phi_h.iter_mut().zip(extra) {
                *phi += b;
            }
        }
        let mut p = self.pressure.clone();
        let fixed = self.fixed_pressure();
        let mut result = empty_result();
        let mut solves = 0;
        for _ in 0..self.config.nonorthogonal_correctors + 1 {
            let correction = self.nonorthogonal_flux(&p)?;
            let pb = boundary_pressure_values_with(m, &p, &self.pressure_boundary, &self.profiles)?;
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
            // The linear solvers floor their tolerance at an absolute rms of 1, which for
            // a right-hand side in m³/s means barely solving at all. Solve for p/s with
            // the right-hand side divided by the cell flux scale s, so the tolerance is
            // relative to the flow actually passing through the cells.
            let flux_scale = if self.config.piso_neighbour_update {
                flux_scale(m, &phi_h)
            } else {
                1.0
            };
            for v in rhs.iter_mut() {
                *v /= flux_scale;
            }
            for v in p.iter_mut() {
                *v /= flux_scale;
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
            for v in p.iter_mut() {
                *v *= flux_scale;
            }
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
        let pb = boundary_pressure_values_with(m, &self.pressure, &self.pressure_boundary, &self.profiles)?;
        let gp = if body_flux.is_some() {
            self.balanced_residual(&self.pressure)?
        } else {
            gauss_gradient_scalar(m, &self.pressure, &pb)?
        };
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
        self.face_flux = rhie_chow_face_flux_with(
            m,
            &self.h_by_a,
            &self.pressure_mobility,
            &self.pressure,
            &self.velocity_boundary,
            &self.pressure_boundary,
            true,
            &self.profiles,
        )?;
        if let Some(extra) = &body_flux {
            for (phi, b) in self.face_flux.iter_mut().zip(extra) {
                *phi += b;
            }
        }
        Ok((result, solves))
    }
    /// Whether a face's Rhie-Chow flux is driven by the pressure difference across it (and so
    /// also carries the face body force): interior faces, and boundary faces with a fixed
    /// pressure whose velocity is not prescribed.
    fn pressure_driven(&self, fi: usize) -> bool {
        let f = self.mesh.faces()[fi];
        match f.patch {
            None => true,
            Some(p) => {
                self.pressure_boundary[p].boundary_type == PressureBoundaryType::FixedValue
                    && !matches!(
                        self.velocity_boundary[p].boundary_type,
                        VelocityBoundaryType::FixedValue | VelocityBoundaryType::Slip
                    )
            }
        }
    }
    /// Face body force a_f: linear between the two cells, the owner's on a boundary.
    fn face_body_force(&self, fi: usize) -> Result<Vec3> {
        let m = &self.mesh;
        let f = m.faces()[fi];
        if f.boundary() {
            Ok(self.body_acceleration[f.owner])
        } else {
            linear_vector(m, f, &self.body_acceleration)
        }
    }
    /// Balanced-force face flux of the body force, d_f (a_f · S), on pressure-driven faces.
    fn body_flux(&self) -> Result<Vec<f64>> {
        let m = &self.mesh;
        let mut out = vec![0.0; m.face_count()];
        for (fi, slot) in out.iter_mut().enumerate() {
            if !self.pressure_driven(fi) {
                continue;
            }
            let f = m.faces()[fi];
            let mobility = linear_scalar(m, f, &self.pressure_mobility)?;
            *slot = mobility * self.face_body_force(fi)?.dot(f.area);
        }
        Ok(out)
    }
    /// The cell value of (∇p − a) rebuilt by least squares from its faces' normal components,
    /// G_f |S| = (p_N − p_P) S·d/|d|² − a_f · S, so a face in balance contributes nothing.
    /// Faces with a zero-gradient pressure (walls) are taken as balanced; they still enter the
    /// normal matrix, which keeps it invertible.
    fn balanced_residual(&self, p: &[f64]) -> Result<Vec<Vec3>> {
        let m = &self.mesh;
        let n = m.cell_count();
        let pb = boundary_pressure_values_with(m, p, &self.pressure_boundary, &self.profiles)?;
        let mut matrix = vec![[[0.0f64; 3]; 3]; n];
        let mut sum = vec![Vec3::ZERO; n];
        for fi in 0..m.face_count() {
            let f = m.faces()[fi];
            let area = f.area.norm();
            if !(area > 0.0) {
                continue;
            }
            let unit = f.area * (1.0 / area);
            let u = [unit.x, unit.y, unit.z];
            let weight = {
                let mut w = [[0.0; 3]; 3];
                for (i, row) in w.iter_mut().enumerate() {
                    for (j, v) in row.iter_mut().enumerate() {
                        *v = u[i] * u[j] * area;
                    }
                }
                w
            };
            let geom = decompose_face_area(m, fi)?;
            let imbalance = match f.neighbour {
                Some(nb) => Some(geom.orthogonal_metric * (p[nb] - p[f.owner]) - self.face_body_force(fi)?.dot(f.area)),
                None if self.pressure_driven(fi)
                    || f.patch.is_some_and(|b| self.pressure_boundary[b].boundary_type == PressureBoundaryType::FixedValue) =>
                {
                    Some(geom.orthogonal_metric * (pb[fi] - p[f.owner]) - self.face_body_force(fi)?.dot(f.area))
                }
                None => None,
            };
            let mut add = |c: usize| {
                for i in 0..3 {
                    for j in 0..3 {
                        matrix[c][i][j] += weight[i][j];
                    }
                }
                if let Some(g) = imbalance {
                    sum[c] += unit * g;
                }
            };
            add(f.owner);
            if let Some(nb) = f.neighbour {
                add(nb);
            }
        }
        let mut out = Vec::with_capacity(n);
        for c in 0..n {
            let x = solve3(matrix[c], sum[c]).ok_or_else(|| {
                SolverError::NumericalFailure("singular face-normal matrix in the balanced-force reconstruction".into())
            })?;
            out.push(x);
        }
        Ok(out)
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
        for k in 0..correctors {
            if k > 0 && physical && self.config.piso_neighbour_update {
                self.refresh_h_by_a()?;
            }
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
