use super::{models::velocity_gradients, operators::scalar_sizes, PolyMesh};
use crate::{core::Vec3, finite_slice, positive, require, Result, SolverError};
use std::collections::BTreeMap;
#[derive(Clone, Debug)]
pub struct DerivedFlowFields {
    pub vorticity: Vec<Vec3>,
    pub q_criterion: Vec<f64>,
    pub lambda2: Vec<f64>,
    pub enstrophy: Vec<f64>,
}
fn eigen_symmetric(mut a: [[f64; 3]; 3]) -> [f64; 3] {
    for _ in 0..20 {
        let (mut p, mut q) = (0, 1);
        let mut mx = a[0][1].abs();
        for (i, j) in [(0, 2), (1, 2)] {
            if a[i][j].abs() > mx {
                p = i;
                q = j;
                mx = a[p][q].abs();
            }
        }
        if mx < 1e-14 {
            break;
        }
        let phi = 0.5 * (2.0 * a[p][q]).atan2(a[q][q] - a[p][p]);
        let (c, s) = (phi.cos(), phi.sin());
        for k in 0..3 {
            let (ap, aq) = (a[p][k], a[q][k]);
            a[p][k] = c * ap - s * aq;
            a[q][k] = s * ap + c * aq;
        }
        for k in 0..3 {
            let (ap, aq) = (a[k][p], a[k][q]);
            a[k][p] = c * ap - s * aq;
            a[k][q] = s * ap + c * aq;
        }
    }
    let mut e = [a[0][0], a[1][1], a[2][2]];
    e.sort_by(f64::total_cmp);
    e
}
pub fn velocity_derived_fields(mesh: &PolyMesh, u: &[Vec3]) -> Result<DerivedFlowFields> {
    let gradients = velocity_gradients(mesh, u)?;
    let n = mesh.cell_count();
    let mut out = DerivedFlowFields {
        vorticity: vec![Vec3::ZERO; n],
        q_criterion: vec![0.0; n],
        lambda2: vec![0.0; n],
        enstrophy: vec![0.0; n],
    };
    for (c, g) in gradients.iter().enumerate() {
        let mut s = [[0.0; 3]; 3];
        let mut o = [[0.0; 3]; 3];
        let (mut ss, mut oo) = (0.0, 0.0);
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] = 0.5 * (g[i][j] + g[j][i]);
                o[i][j] = 0.5 * (g[i][j] - g[j][i]);
                ss += s[i][j] * s[i][j];
                oo += o[i][j] * o[i][j];
            }
        }
        out.vorticity[c] = Vec3::new(g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]);
        out.q_criterion[c] = 0.5 * (oo - ss);
        out.enstrophy[c] = 0.5 * out.vorticity[c].norm_squared();
        let mut m = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    m[i][j] += s[i][k] * s[k][j] + o[i][k] * o[k][j];
                }
            }
        }
        out.lambda2[c] = eigen_symmetric(m)[1];
    }
    finite_slice(&out.q_criterion, "Q criterion")?;
    finite_slice(&out.lambda2, "lambda2")?;
    finite_slice(&out.enstrophy, "enstrophy")?;
    Ok(out)
}
pub fn wall_y_plus(
    distance: f64,
    shear_stress: f64,
    density: f64,
    dynamic_viscosity: f64,
) -> Result<f64> {
    require(
        distance.is_finite() && distance >= 0.0 && shear_stress.is_finite() && shear_stress >= 0.0,
        "invalid wall distance/stress",
    )?;
    positive(density, "density")?;
    positive(dynamic_viscosity, "dynamic viscosity")?;
    Ok(distance * (shear_stress / density).sqrt() / (dynamic_viscosity / density))
}
pub fn pressure_force_on_patch(mesh: &PolyMesh, pressure: &[f64], patch: usize) -> Result<Vec3> {
    scalar_sizes(mesh, pressure, &[])?;
    require(patch < mesh.patches().len(), "pressure-force patch")?;
    let mut f = Vec3::ZERO;
    for face in mesh.faces() {
        if face.patch == Some(patch) {
            f -= face.area * pressure[face.owner];
        }
    }
    Ok(f)
}
pub fn sample_nearest_cell(mesh: &PolyMesh, field: &[f64], point: Vec3) -> Result<f64> {
    scalar_sizes(mesh, field, &[])?;
    require(point.is_finite(), "sample point must be finite")?;
    let mut best = f64::INFINITY;
    let mut index = 0;
    for (i, c) in mesh.cells().iter().enumerate() {
        let d = (c.center - point).norm_squared();
        if d < best {
            best = d;
            index = i;
        }
    }
    Ok(field[index])
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct SolverReport {
    pub step: usize,
    pub time: f64,
    pub residuals: BTreeMap<String, f64>,
}
pub type FunctionObject<'a> = dyn Fn(&SolverReport) -> Result<()> + 'a;
#[derive(Default)]
pub struct FunctionObjectManager<'a> {
    objects: Vec<Box<FunctionObject<'a>>>,
}
impl<'a> FunctionObjectManager<'a> {
    pub fn add(&mut self, f: impl Fn(&SolverReport) -> Result<()> + 'a) {
        self.objects.push(Box::new(f));
    }
    pub fn execute(&self, r: &SolverReport) -> Result<()> {
        require(
            r.time.is_finite() && r.residuals.values().all(|v| v.is_finite()),
            "nonfinite solver report",
        )?;
        for f in &self.objects {
            f(r)?;
        }
        Ok(())
    }
    pub fn size(&self) -> usize {
        self.objects.len()
    }
}
#[derive(Clone, Debug, Default)]
pub struct ResidualHistory {
    records: Vec<SolverReport>,
}
impl ResidualHistory {
    pub fn observe(&mut self, r: SolverReport) -> Result<()> {
        require(
            r.time.is_finite() && r.residuals.values().all(|v| v.is_finite()),
            "nonfinite residual report",
        )?;
        self.records.push(r);
        Ok(())
    }
    pub fn records(&self) -> &[SolverReport] {
        &self.records
    }
    pub fn latest(&self, field: &str) -> Result<f64> {
        self.records
            .last()
            .ok_or_else(|| SolverError::InvalidInput("empty residual history".into()))?
            .residuals
            .get(field)
            .copied()
            .ok_or_else(|| SolverError::InvalidInput(format!("missing residual {field}")))
    }
}
