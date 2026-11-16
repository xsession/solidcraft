use super::{decompose_face_area, schemes::flux_sizes, PolyMesh};
use crate::{
    core::{
        restarted_gmres, CsrBuilder, CsrMatrix, Ilu0Preconditioner, IterativeSolverResult,
        LinearConfig, Vec3,
    },
    finite_slice, positive, require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TemporalScheme {
    Euler,
    BackwardBdf2,
    CrankNicolson,
}
impl Default for TemporalScheme {
    fn default() -> Self {
        Self::Euler
    }
}
pub fn bdf2_derivative(current: f64, old: f64, older: f64, dt: f64) -> Result<f64> {
    positive(dt, "BDF2 dt")?;
    let d = (1.5 * current - 2.0 * old + 0.5 * older) / dt;
    require(d.is_finite(), "nonfinite BDF2 derivative")?;
    Ok(d)
}
pub fn crank_nicolson_update(old: f64, rhs_old: f64, rhs_new: f64, dt: f64) -> Result<f64> {
    positive(dt, "CN dt")?;
    let v = old + 0.5 * dt * (rhs_old + rhs_new);
    require(v.is_finite(), "nonfinite CN update")?;
    Ok(v)
}
pub fn crank_nicolson_implicit_weight(off_centering: f64) -> Result<f64> {
    require(
        off_centering.is_finite() && (0.0..=1.0).contains(&off_centering),
        "CN off-centering must be in [0,1]",
    )?;
    Ok(1.0 / (1.0 + off_centering))
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ScalarBoundaryType {
    ZeroGradient,
    FixedValue,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct ScalarBoundaryCondition {
    pub boundary_type: ScalarBoundaryType,
    pub value: f64,
}
impl Default for ScalarBoundaryCondition {
    fn default() -> Self {
        Self {
            boundary_type: ScalarBoundaryType::ZeroGradient,
            value: 0.0,
        }
    }
}
#[derive(Clone, Debug)]
pub struct FvScalarMatrix {
    builder: CsrBuilder,
    source: Vec<f64>,
}
impl FvScalarMatrix {
    pub fn new(cells: usize) -> Result<Self> {
        require(cells > 0, "FV matrix needs cells")?;
        Ok(Self {
            builder: CsrBuilder::new(cells, cells),
            source: vec![0.0; cells],
        })
    }
    pub fn add_diagonal(&mut self, cell: usize, value: f64) -> Result<()> {
        self.builder.add(cell, cell, value)
    }
    pub fn add_coefficient(&mut self, row: usize, column: usize, value: f64) -> Result<()> {
        self.builder.add(row, column, value)
    }
    pub fn add_source(&mut self, cell: usize, value: f64) -> Result<()> {
        require(
            cell < self.source.len() && value.is_finite(),
            "FV source index/value",
        )?;
        let s = self.source[cell] + value;
        require(s.is_finite(), "FV source overflow")?;
        self.source[cell] = s;
        Ok(())
    }
    pub fn matrix(&self) -> Result<CsrMatrix> {
        self.builder.build(0.0)
    }
    pub fn source(&self) -> &[f64] {
        &self.source
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ScalarTransportConfig {
    pub dt: f64,
    pub diffusivity: f64,
    pub linear: LinearConfig,
    pub temporal_scheme: TemporalScheme,
    pub crank_nicolson_off_centering: f64,
}
impl Default for ScalarTransportConfig {
    fn default() -> Self {
        Self {
            dt: 1e-3,
            diffusivity: 1e-2,
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
#[derive(Clone, Debug)]
pub struct ScalarTransport {
    mesh: PolyMesh,
    config: ScalarTransportConfig,
    boundary: Vec<ScalarBoundaryCondition>,
    values: Vec<f64>,
    old_values: Vec<f64>,
    source: Vec<f64>,
    face_flux: Vec<f64>,
    spatial: Option<(CsrMatrix, Vec<f64>)>,
    operator_assemblies: usize,
    steps: usize,
    time: f64,
    linear_result: Option<IterativeSolverResult>,
}
impl ScalarTransport {
    pub fn new(mesh: PolyMesh, config: ScalarTransportConfig) -> Result<Self> {
        positive(config.dt, "scalar dt")?;
        require(
            config.diffusivity.is_finite()
                && config.diffusivity >= 0.0
                && config.linear.max_iterations > 0
                && config.linear.restart > 0,
            "invalid scalar controls",
        )?;
        positive(config.linear.relative_tolerance, "linear tolerance")?;
        if config.temporal_scheme == TemporalScheme::CrankNicolson {
            crank_nicolson_implicit_weight(config.crank_nicolson_off_centering)?;
        }
        let n = mesh.cell_count();
        let nf = mesh.face_count();
        let patches = mesh.patches().len();
        Ok(Self {
            mesh,
            config,
            boundary: vec![ScalarBoundaryCondition::default(); patches],
            values: vec![0.0; n],
            old_values: vec![0.0; n],
            source: vec![0.0; n],
            face_flux: vec![0.0; nf],
            spatial: None,
            operator_assemblies: 0,
            steps: 0,
            time: 0.0,
            linear_result: None,
        })
    }
    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    pub fn face_flux(&self) -> &[f64] {
        &self.face_flux
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn linear_result(&self) -> Option<IterativeSolverResult> {
        self.linear_result
    }
    pub fn operator_assemblies(&self) -> usize {
        self.operator_assemblies
    }
    pub fn set_boundary(
        &mut self,
        patch: &str,
        boundary_type: ScalarBoundaryType,
        value: f64,
    ) -> Result<()> {
        let index = self.mesh.patch_index(patch)?;
        self.set_boundary_index(index, boundary_type, value)
    }
    pub fn set_boundary_index(
        &mut self,
        index: usize,
        boundary_type: ScalarBoundaryType,
        value: f64,
    ) -> Result<()> {
        require(
            index < self.boundary.len() && value.is_finite(),
            "scalar boundary index/value",
        )?;
        self.boundary[index] = ScalarBoundaryCondition {
            boundary_type,
            value,
        };
        self.spatial = None;
        Ok(())
    }
    pub fn initialize_constant(&mut self, value: f64) -> Result<()> {
        self.initialize_profile(|_| value)
    }
    pub fn initialize_profile(&mut self, f: impl Fn(Vec3) -> f64) -> Result<()> {
        let values: Vec<_> = self.mesh.cells().iter().map(|c| f(c.center)).collect();
        finite_slice(&values, "initial scalar values")?;
        self.old_values = values.clone();
        self.values = values;
        self.time = 0.0;
        self.steps = 0;
        self.linear_result = None;
        Ok(())
    }
    pub fn set_source_constant(&mut self, value: f64) -> Result<()> {
        self.set_source_profile(|_| value)
    }
    pub fn set_source_profile(&mut self, f: impl Fn(Vec3) -> f64) -> Result<()> {
        let s: Vec<_> = self.mesh.cells().iter().map(|c| f(c.center)).collect();
        finite_slice(&s, "scalar source")?;
        self.source = s;
        Ok(())
    }
    pub fn set_face_flux(&mut self, flux: Vec<f64>) -> Result<()> {
        flux_sizes(&self.mesh, &flux)?;
        self.face_flux = flux;
        self.spatial = None;
        Ok(())
    }
    fn assemble_operator(&mut self) -> Result<()> {
        let m = &self.mesh;
        let n = m.cell_count();
        let mut b = CsrBuilder::new(n, n);
        let mut boundary_rhs = vec![0.0; n];
        for c in 0..n {
            let (mut ap, mut rhs) = (0.0, 0.0);
            for &fi in &m.cell_faces()[c] {
                let f = m.faces()[fi];
                let phi = if f.owner == c {
                    self.face_flux[fi]
                } else {
                    -self.face_flux[fi]
                };
                if let Some(neighbour) = f.neighbour {
                    let other = if f.owner == c { neighbour } else { f.owner };
                    let diff =
                        self.config.diffusivity * decompose_face_area(m, fi)?.orthogonal_metric;
                    ap += diff + phi.max(0.0);
                    let off = -diff + phi.min(0.0);
                    if off != 0.0 {
                        b.add(c, other, off)?;
                    }
                } else {
                    let bc = self.boundary[f.patch.ok_or_else(|| {
                        crate::SolverError::InvalidInput("boundary face without patch".into())
                    })?];
                    if bc.boundary_type == ScalarBoundaryType::FixedValue {
                        let distance = (f.center - m.cells()[c].center).norm();
                        positive(distance, "scalar boundary distance")?;
                        let diff = self.config.diffusivity * f.area.norm() / distance;
                        ap += diff + phi.max(0.0);
                        rhs += diff * bc.value - phi.min(0.0) * bc.value;
                    } else {
                        ap += phi;
                    }
                }
            }
            require(
                ap.is_finite() && rhs.is_finite(),
                "nonfinite transport diagonal/load",
            )?;
            if ap != 0.0 {
                b.add(c, c, ap)?;
            }
            boundary_rhs[c] = rhs;
        }
        self.spatial = Some((b.build(0.0)?, boundary_rhs));
        self.operator_assemblies += 1;
        Ok(())
    }
    pub fn step(&mut self) -> Result<IterativeSolverResult> {
        if self.spatial.is_none() {
            self.assemble_operator()?;
        }
        let (spatial, boundary_rhs) = self.spatial.as_ref().ok_or_else(|| {
            crate::SolverError::NumericalFailure("transport operator was not assembled".into())
        })?;
        let n = self.mesh.cell_count();
        let use_bdf2 =
            self.config.temporal_scheme == TemporalScheme::BackwardBdf2 && self.steps > 0;
        let alpha0 = if use_bdf2 { 1.5 } else { 1.0 };
        let theta = if self.config.temporal_scheme == TemporalScheme::CrankNicolson {
            crank_nicolson_implicit_weight(self.config.crank_nicolson_off_centering)?
        } else {
            1.0
        };
        let mut b = CsrBuilder::new(n, n);
        for row in 0..n {
            for k in spatial.row_offsets()[row]..spatial.row_offsets()[row + 1] {
                b.add(row, spatial.columns()[k], theta * spatial.values()[k])?;
            }
            b.add(
                row,
                row,
                alpha0 * self.mesh.cells()[row].volume / self.config.dt,
            )?;
        }
        let a = b.build(0.0)?;
        let ilu = Ilu0Preconditioner::new(&a, 1e-12)?;
        let action = if theta < 1.0 {
            spatial.apply(&self.values)?
        } else {
            vec![0.0; n]
        };
        let mut rhs = vec![0.0; n];
        for c in 0..n {
            let vol = self.mesh.cells()[c].volume;
            let temporal = if use_bdf2 {
                2.0 * self.values[c] - 0.5 * self.old_values[c]
            } else {
                self.values[c]
            };
            rhs[c] = temporal * vol / self.config.dt + self.source[c] * vol + boundary_rhs[c]
                - (1.0 - theta) * action[c];
        }
        let mut next = self.values.clone();
        let result = restarted_gmres(
            &rhs,
            &mut next,
            |v, y| a.multiply(v, y),
            &ilu,
            self.config.linear,
        )?;
        self.linear_result = Some(result);
        if !result.converged {
            return Err(SolverError::NumericalFailure(format!(
                "scalar transport solve failed, residual {}",
                result.residual_rms
            )));
        }
        finite_slice(&next, "scalar update")?;
        self.old_values = std::mem::replace(&mut self.values, next);
        self.steps += 1;
        self.time += self.config.dt;
        Ok(result)
    }
    pub fn run(&mut self, steps: usize) -> Result<Option<IterativeSolverResult>> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(self.linear_result)
    }
    pub fn volume_integral(&self) -> f64 {
        self.values
            .iter()
            .zip(self.mesh.cells())
            .map(|(v, c)| v * c.volume)
            .sum()
    }
    pub fn minimum(&self) -> f64 {
        self.values.iter().copied().fold(f64::INFINITY, f64::min)
    }
    pub fn maximum(&self) -> f64 {
        self.values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}
