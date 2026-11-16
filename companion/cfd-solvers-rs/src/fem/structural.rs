use super::{scalar::tri_geometry, Mesh2D, Node2};
use crate::{
    core::{
        conjugate_gradient, lowest_generalized_eigenpairs, CsrBuilder, GeneralizedEigenConfig,
        GeneralizedEigenpair, IterativeSolverResult, LinearConfig,
    },
    finite_slice, positive, require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ElasticityMode2D {
    PlaneStress,
    PlaneStrain,
}
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
pub struct Displacement2 {
    pub x: f64,
    pub y: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Elasticity2DConfig {
    pub youngs_modulus: f64,
    pub poisson_ratio: f64,
    pub mode: ElasticityMode2D,
    pub linear: LinearConfig,
}
impl Default for Elasticity2DConfig {
    fn default() -> Self {
        Self {
            youngs_modulus: 1.0,
            poisson_ratio: 0.3,
            mode: ElasticityMode2D::PlaneStress,
            linear: LinearConfig {
                max_iterations: 4000,
                ..Default::default()
            },
        }
    }
}
impl Elasticity2DConfig {
    fn constitutive(self) -> Result<[[f64; 3]; 3]> {
        positive(self.youngs_modulus, "Young modulus")?;
        let nu = self.poisson_ratio;
        require(
            nu.is_finite() && nu > -1.0 && nu < 0.5,
            "invalid Poisson ratio",
        )?;
        let e = self.youngs_modulus;
        Ok(if self.mode == ElasticityMode2D::PlaneStress {
            let f = e / (1.0 - nu * nu);
            [
                [f, f * nu, 0.0],
                [f * nu, f, 0.0],
                [0.0, 0.0, f * (1.0 - nu) / 2.0],
            ]
        } else {
            let f = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
            [
                [f * (1.0 - nu), f * nu, 0.0],
                [f * nu, f * (1.0 - nu), 0.0],
                [0.0, 0.0, f * (1.0 - 2.0 * nu) / 2.0],
            ]
        })
    }
}
#[derive(Clone, Debug)]
pub struct Elasticity2D {
    mesh: Mesh2D,
    config: Elasticity2DConfig,
    displacement: Vec<Displacement2>,
    prescribed: Vec<Option<f64>>,
    pub linear_result: Option<IterativeSolverResult>,
}
impl Elasticity2D {
    pub fn new(mesh: Mesh2D, config: Elasticity2DConfig) -> Result<Self> {
        mesh.validate()?;
        config.constitutive()?;
        let n = mesh.node_count();
        require(n <= usize::MAX / 2, "elasticity DOF count overflow")?;
        Ok(Self {
            mesh,
            config,
            displacement: vec![Displacement2::default(); n],
            prescribed: vec![None; 2 * n],
            linear_result: None,
        })
    }
    pub fn mesh(&self) -> &Mesh2D {
        &self.mesh
    }
    pub fn displacement(&self) -> &[Displacement2] {
        &self.displacement
    }
    pub fn set_dirichlet(
        &mut self,
        predicate: impl Fn(Node2) -> bool,
        value: impl Fn(Node2) -> Displacement2,
        fix_x: bool,
        fix_y: bool,
    ) -> Result<()> {
        for (i, &p) in self.mesh.nodes.iter().enumerate() {
            if predicate(p) {
                let v = value(p);
                require(
                    v.x.is_finite() && v.y.is_finite(),
                    "nonfinite prescribed displacement",
                )?;
                if fix_x {
                    self.prescribed[2 * i] = Some(v.x);
                }
                if fix_y {
                    self.prescribed[2 * i + 1] = Some(v.y);
                }
            }
        }
        Ok(())
    }
    pub fn solve(&mut self, body_force: impl Fn(Node2) -> Displacement2) -> Result<()> {
        self.solve_impl(body_force, None, 0.0, 0.0)
    }
    pub fn solve_thermal(
        &mut self,
        temperature: &[f64],
        expansion: f64,
        reference: f64,
        body_force: impl Fn(Node2) -> Displacement2,
    ) -> Result<()> {
        require(
            temperature.len() == self.mesh.node_count()
                && expansion.is_finite()
                && reference.is_finite(),
            "thermoelastic input",
        )?;
        finite_slice(temperature, "temperature")?;
        self.solve_impl(body_force, Some(temperature), expansion, reference)
    }
    fn solve_impl(
        &mut self,
        body_force: impl Fn(Node2) -> Displacement2,
        temperature: Option<&[f64]>,
        expansion: f64,
        reference: f64,
    ) -> Result<()> {
        let ndof = self.prescribed.len();
        let mut map = vec![None; ndof];
        let mut nf = 0;
        for i in 0..ndof {
            if self.prescribed[i].is_none() {
                map[i] = Some(nf);
                nf += 1;
            }
        }
        require(nf > 0, "elasticity system has no free DOFs")?;
        let mut builder = CsrBuilder::new(nf, nf);
        let mut rhs = vec![0.0; nf];
        let d = self.config.constitutive()?;
        for tri in &self.mesh.triangles {
            let (area, g, centroid) = tri_geometry(&self.mesh, tri);
            let mut b = [[0.0; 6]; 3];
            for i in 0..3 {
                b[0][2 * i] = g[i].x;
                b[1][2 * i + 1] = g[i].y;
                b[2][2 * i] = g[i].y;
                b[2][2 * i + 1] = g[i].x;
            }
            let force = body_force(centroid);
            require(
                force.x.is_finite() && force.y.is_finite(),
                "nonfinite body force",
            )?;
            let mut thermal = temperature.map_or(0.0, |t| {
                expansion * ((t[tri.node[0]] + t[tri.node[1]] + t[tri.node[2]]) / 3.0 - reference)
            });
            if self.config.mode == ElasticityMode2D::PlaneStrain {
                thermal *= 1.0 + self.config.poisson_ratio;
            }
            let stress: Vec<_> = (0..3).map(|a| (d[a][0] + d[a][1]) * thermal).collect();
            for i in 0..6 {
                let gi = 2 * tri.node[i / 2] + i % 2;
                let Some(row) = map[gi] else {
                    continue;
                };
                rhs[row] += if i % 2 == 0 {
                    force.x * area / 3.0
                } else {
                    force.y * area / 3.0
                };
                for a in 0..3 {
                    rhs[row] += area * b[a][i] * stress[a];
                }
                for j in 0..6 {
                    let gj = 2 * tri.node[j / 2] + j % 2;
                    let mut kij = 0.0;
                    for a in 0..3 {
                        for c in 0..3 {
                            kij += area * b[a][i] * d[a][c] * b[c][j];
                        }
                    }
                    if let Some(col) = map[gj] {
                        builder.add(row, col, kij)?;
                    } else {
                        rhs[row] -= kij * self.prescribed[gj].unwrap_or(0.0);
                    }
                }
            }
        }
        let a = builder.build(0.0)?;
        let mut x = vec![0.0; nf];
        let result = conjugate_gradient(&a, &rhs, &mut x, self.config.linear)?;
        if !result.converged {
            return Err(SolverError::NumericalFailure(
                "Elasticity2D PCG failed".into(),
            ));
        }
        for n in 0..self.mesh.node_count() {
            self.displacement[n] = Displacement2 {
                x: map[2 * n].map_or_else(|| self.prescribed[2 * n].unwrap_or(0.0), |i| x[i]),
                y: map[2 * n + 1]
                    .map_or_else(|| self.prescribed[2 * n + 1].unwrap_or(0.0), |i| x[i]),
            };
        }
        self.linear_result = Some(result);
        Ok(())
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct BarModal1DConfig {
    pub elements: usize,
    pub length: f64,
    pub area: f64,
    pub young_modulus: f64,
    pub density: f64,
    pub eigen: GeneralizedEigenConfig,
}
impl Default for BarModal1DConfig {
    fn default() -> Self {
        Self {
            elements: 80,
            length: 1.0,
            area: 1.0,
            young_modulus: 200e9,
            density: 7800.0,
            eigen: Default::default(),
        }
    }
}
#[derive(Clone, Debug, Serialize)]
pub struct BarMode1D {
    pub angular_frequency: f64,
    pub frequency_hz: f64,
    pub displacement: Vec<f64>,
    pub eigenpair: GeneralizedEigenpair,
}
pub fn bar_modal1d(config: BarModal1DConfig, mode_count: usize) -> Result<Vec<BarMode1D>> {
    require(config.elements >= 2, "modal bar needs at least 2 elements")?;
    for (v, n) in [
        (config.length, "length"),
        (config.area, "area"),
        (config.young_modulus, "Young modulus"),
        (config.density, "density"),
    ] {
        positive(v, n)?;
    }
    let n = config.elements;
    let h = config.length / n as f64;
    let ke = config.young_modulus * config.area / h;
    let me = config.density * config.area * h / 6.0;
    let mut kb = CsrBuilder::new(n, n);
    let mut mb = CsrBuilder::new(n, n);
    for e in 0..n {
        let global = [e, e + 1];
        for i in 0..2 {
            if global[i] == 0 {
                continue;
            }
            for j in 0..2 {
                if global[j] == 0 {
                    continue;
                }
                kb.add(global[i] - 1, global[j] - 1, if i == j { ke } else { -ke })?;
                mb.add(
                    global[i] - 1,
                    global[j] - 1,
                    if i == j { 2.0 * me } else { me },
                )?;
            }
        }
    }
    let pairs = lowest_generalized_eigenpairs(
        &kb.build(0.0)?,
        &mb.build(0.0)?,
        GeneralizedEigenConfig {
            mode_count,
            ..config.eigen
        },
    )?;
    require(
        pairs.len() == mode_count,
        "bar eigen solve did not converge",
    )?;
    let mut result = vec![];
    for pair in pairs {
        positive(pair.eigenvalue, "modal eigenvalue")?;
        let omega = pair.eigenvalue.sqrt();
        let mut displacement = vec![0.0; n + 1];
        displacement[1..].copy_from_slice(&pair.vector);
        result.push(BarMode1D {
            angular_frequency: omega,
            frequency_hz: omega / (2.0 * std::f64::consts::PI),
            displacement,
            eigenpair: pair,
        });
    }
    Ok(result)
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PiezoelectricBar1DConfig {
    pub elements: usize,
    pub length_m: f64,
    pub area_m2: f64,
    pub elastic_modulus_pa: f64,
    pub piezoelectric_stress_c_per_m2: f64,
    pub permittivity_f_per_m: f64,
}
impl Default for PiezoelectricBar1DConfig {
    fn default() -> Self {
        Self {
            elements: 16,
            length_m: 1e-3,
            area_m2: 1e-4,
            elastic_modulus_pa: 6e10,
            piezoelectric_stress_c_per_m2: 15.0,
            permittivity_f_per_m: 1.5e-8,
        }
    }
}
#[derive(Clone, Debug)]
pub struct PiezoelectricBar1D {
    config: PiezoelectricBar1DConfig,
    displacement: Vec<f64>,
    potential: Vec<f64>,
    electrode_charge: f64,
}
impl PiezoelectricBar1D {
    pub fn new(c: PiezoelectricBar1DConfig) -> Result<Self> {
        require(
            c.elements > 0 && c.elements < usize::MAX / 2,
            "invalid piezo element count",
        )?;
        for (v, n) in [
            (c.length_m, "length"),
            (c.area_m2, "area"),
            (c.elastic_modulus_pa, "modulus"),
            (c.permittivity_f_per_m, "permittivity"),
        ] {
            positive(v, n)?;
        }
        require(
            c.piezoelectric_stress_c_per_m2.is_finite(),
            "piezo stress coefficient",
        )?;
        Ok(Self {
            config: c,
            displacement: vec![0.0; c.elements + 1],
            potential: vec![0.0; c.elements + 1],
            electrode_charge: 0.0,
        })
    }
    pub fn config(&self) -> PiezoelectricBar1DConfig {
        self.config
    }
    pub fn displacement_m(&self) -> &[f64] {
        &self.displacement
    }
    pub fn potential_v(&self) -> &[f64] {
        &self.potential
    }
    pub fn open_circuit_voltage_v(&self) -> f64 {
        self.potential.last().copied().unwrap_or(0.0)
    }
    pub fn tip_displacement_m(&self) -> f64 {
        self.displacement.last().copied().unwrap_or(0.0)
    }
    pub fn electrode_charge_c(&self) -> f64 {
        self.electrode_charge
    }
    pub fn solve_receive(&mut self, force: f64) -> Result<()> {
        self.solve(force, None)
    }
    pub fn solve_transmit(&mut self, voltage: f64, force: f64) -> Result<()> {
        require(voltage.is_finite(), "nonfinite voltage")?;
        self.solve(force, Some(voltage))
    }
    fn solve(&mut self, force: f64, voltage: Option<f64>) -> Result<()> {
        require(force.is_finite(), "nonfinite tip force")?;
        let c = self.config;
        let nodes = c.elements + 1;
        let ndof = 2 * nodes;
        require(
            ndof.checked_mul(ndof).is_some(),
            "piezo matrix size overflow",
        )?;
        let mut a = vec![0.0; ndof * ndof];
        let mut rhs = vec![0.0; ndof];
        let h = c.length_m / c.elements as f64;
        let kuu = c.elastic_modulus_pa * c.area_m2 / h;
        let kup = c.piezoelectric_stress_c_per_m2 * c.area_m2 / h;
        let kpp = -c.permittivity_f_per_m * c.area_m2 / h;
        for e in 0..c.elements {
            for i in 0..2 {
                for j in 0..2 {
                    let sign = if i == j { 1.0 } else { -1.0 };
                    let (ui, uj, pi, pj) = (e + i, e + j, nodes + e + i, nodes + e + j);
                    a[ui * ndof + uj] += kuu * sign;
                    a[ui * ndof + pj] += kup * sign;
                    a[pi * ndof + uj] += kup * sign;
                    a[pi * ndof + pj] += kpp * sign;
                }
            }
        }
        rhs[nodes - 1] = force;
        let mut prescribed = vec![(0, 0.0), (nodes, 0.0)];
        if let Some(v) = voltage {
            prescribed.push((ndof - 1, v));
        }
        for (dof, value) in prescribed {
            for r in 0..ndof {
                if r != dof {
                    rhs[r] -= a[r * ndof + dof] * value;
                    a[r * ndof + dof] = 0.0;
                }
            }
            for j in 0..ndof {
                a[dof * ndof + j] = 0.0;
            }
            a[dof * ndof + dof] = 1.0;
            rhs[dof] = value;
        }
        // Upstream uses an absolute pivot threshold: mixed mechanical/electrical SI units span many orders of magnitude.
        for k in 0..ndof {
            let mut pivot = k;
            for r in k + 1..ndof {
                if a[r * ndof + k].abs() > a[pivot * ndof + k].abs() {
                    pivot = r;
                }
            }
            require(
                a[pivot * ndof + k].abs() > 1e-18,
                "singular piezoelectric matrix",
            )?;
            if pivot != k {
                for j in k..ndof {
                    a.swap(k * ndof + j, pivot * ndof + j);
                }
                rhs.swap(k, pivot);
            }
            for r in k + 1..ndof {
                let factor = a[r * ndof + k] / a[k * ndof + k];
                if factor == 0.0 {
                    continue;
                }
                a[r * ndof + k] = 0.0;
                for j in k + 1..ndof {
                    a[r * ndof + j] -= factor * a[k * ndof + j];
                }
                rhs[r] -= factor * rhs[k];
            }
        }
        for i in (0..ndof).rev() {
            let mut value = rhs[i];
            for j in i + 1..ndof {
                value -= a[i * ndof + j] * rhs[j];
            }
            rhs[i] = value / a[i * ndof + i];
        }
        finite_slice(&rhs, "piezo solution")?;
        self.displacement.copy_from_slice(&rhs[..nodes]);
        self.potential.copy_from_slice(&rhs[nodes..]);
        let strain = (self.displacement[nodes - 1] - self.displacement[nodes - 2]) / h;
        let efield = -(self.potential[nodes - 1] - self.potential[nodes - 2]) / h;
        self.electrode_charge = c.area_m2
            * (c.piezoelectric_stress_c_per_m2 * strain + c.permittivity_f_per_m * efield);
        Ok(())
    }
}
