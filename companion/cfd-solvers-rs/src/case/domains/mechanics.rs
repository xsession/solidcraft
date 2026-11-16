//! Solid and multibody mechanics case adapters: 2-D linear elasticity with stress recovery,
//! axial bar modes, J2 plasticity under uniaxial stress, 1-D phase-field fracture, a Cosserat
//! cantilever rod and explicit Hertz–Mindlin DEM.
use super::{csv, json, linspace, metric};
use crate::{
    case::{write_field, RunManifest},
    core::Vec3,
    fem::{
        bar_modal1d, make_rectangle_tri_mesh, BarModal1DConfig, Displacement2, Elasticity2D,
        Elasticity2DConfig, ElasticityMode2D, Node2,
    },
    io::{self, Cell, CellType, DataArray, VtuMesh},
    mechanics::{
        CosseratRodConfig, J2Plasticity3D, J2Plasticity3DConfig, J2Plasticity3DState,
        PhaseFieldFracture1D, PhaseFieldFracture1DConfig, StableCosseratRod,
    },
    multibody::{DemConfig, ExplicitDemSystem, HertzMindlinContactModel, PlaneShape},
    require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::path::Path;

// ---------------------------------------------------------------- 2-D elasticity

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Edge {
    Left,
    Right,
    Bottom,
    Top,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PlaneMode {
    #[default]
    PlaneStress,
    PlaneStrain,
}

/// Prescribed edge displacement; an omitted component stays free.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdgeDisplacement {
    pub edge: Edge,
    #[serde(default)]
    pub ux_m: Option<f64>,
    #[serde(default)]
    pub uy_m: Option<f64>,
}

fn poisson() -> f64 {
    0.3
}

/// Rectangular plate on linear triangles with clamped/prescribed edges, gravity body load
/// and an optional uniform temperature rise; reports displacement and element stresses.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Elasticity2dCase {
    pub width_m: f64,
    pub height_m: f64,
    pub nx: usize,
    pub ny: usize,
    pub youngs_modulus_pa: f64,
    #[serde(default = "poisson")]
    pub poisson_ratio: f64,
    #[serde(default)]
    pub mode: PlaneMode,
    #[serde(default)]
    pub density_kg_m3: f64,
    #[serde(default)]
    pub gravity_m_s2: [f64; 2],
    pub supports: Vec<EdgeDisplacement>,
    #[serde(default)]
    pub temperature_rise_k: f64,
    #[serde(default)]
    pub expansion_per_k: f64,
}

impl Elasticity2dCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            !self.supports.is_empty() && self.supports.len() <= 8,
            "use 1..=8 edge supports",
        )?;
        require(
            self.poisson_ratio > -1.0 && self.poisson_ratio < 0.5,
            "Poisson ratio must lie in (-1, 0.5)",
        )?;
        let mesh = make_rectangle_tri_mesh(self.nx, self.ny, self.width_m, self.height_m)?;
        let mode = match self.mode {
            PlaneMode::PlaneStress => ElasticityMode2D::PlaneStress,
            PlaneMode::PlaneStrain => ElasticityMode2D::PlaneStrain,
        };
        let mut s = Elasticity2D::new(
            mesh,
            Elasticity2DConfig {
                youngs_modulus: self.youngs_modulus_pa,
                poisson_ratio: self.poisson_ratio,
                mode,
                ..Default::default()
            },
        )?;
        let (w, h) = (self.width_m, self.height_m);
        let tol = 1e-9 * w.max(h);
        for support in &self.supports {
            let on = move |p: Node2| match support.edge {
                Edge::Left => p.x.abs() <= tol,
                Edge::Right => (p.x - w).abs() <= tol,
                Edge::Bottom => p.y.abs() <= tol,
                Edge::Top => (p.y - h).abs() <= tol,
            };
            let (ux, uy) = (support.ux_m, support.uy_m);
            require(
                ux.is_some() || uy.is_some(),
                "an edge support must fix ux_m, uy_m or both",
            )?;
            s.set_dirichlet(
                on,
                move |_| Displacement2 {
                    x: ux.unwrap_or(0.0),
                    y: uy.unwrap_or(0.0),
                },
                ux.is_some(),
                uy.is_some(),
            )?;
        }
        let rho = self.density_kg_m3;
        let g = self.gravity_m_s2;
        let body = move |_: Node2| Displacement2 {
            x: rho * g[0],
            y: rho * g[1],
        };
        if self.temperature_rise_k != 0.0 {
            let n = s.mesh().node_count();
            s.solve_thermal(
                &vec![300.0 + self.temperature_rise_k; n],
                self.expansion_per_k,
                300.0,
                body,
            )?;
        } else {
            s.solve(body)?;
        }
        // Constant-strain element stresses.
        let (e, nu) = (self.youngs_modulus_pa, self.poisson_ratio);
        let d = match self.mode {
            PlaneMode::PlaneStress => {
                let c = e / (1.0 - nu * nu);
                [
                    [c, c * nu, 0.0],
                    [c * nu, c, 0.0],
                    [0.0, 0.0, c * (1.0 - nu) / 2.0],
                ]
            }
            PlaneMode::PlaneStrain => {
                let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
                [
                    [c * (1.0 - nu), c * nu, 0.0],
                    [c * nu, c * (1.0 - nu), 0.0],
                    [0.0, 0.0, c * (1.0 - 2.0 * nu) / 2.0],
                ]
            }
        };
        let alpha_dt = self.expansion_per_k * self.temperature_rise_k;
        let thermal = match self.mode {
            PlaneMode::PlaneStress => alpha_dt,
            PlaneMode::PlaneStrain => (1.0 + nu) * alpha_dt,
        };
        let nodes = &s.mesh().nodes;
        let u = s.displacement();
        let (mut sx, mut sy, mut txy, mut vm) = (vec![], vec![], vec![], vec![]);
        for t in &s.mesh().triangles {
            let [a, b, c] = t.node;
            let (x1, y1, x2, y2, x3, y3) = (
                nodes[a].x, nodes[a].y, nodes[b].x, nodes[b].y, nodes[c].x, nodes[c].y,
            );
            let area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
            require(area2.abs() > 0.0, "degenerate triangle")?;
            let (b1, b2, b3) = (y2 - y3, y3 - y1, y1 - y2);
            let (c1, c2, c3) = (x3 - x2, x1 - x3, x2 - x1);
            let ex = (b1 * u[a].x + b2 * u[b].x + b3 * u[c].x) / area2 - thermal;
            let ey = (c1 * u[a].y + c2 * u[b].y + c3 * u[c].y) / area2 - thermal;
            let gxy =
                (c1 * u[a].x + c2 * u[b].x + c3 * u[c].x + b1 * u[a].y + b2 * u[b].y + b3 * u[c].y)
                    / area2;
            let (px, py, pxy) = (
                d[0][0] * ex + d[0][1] * ey,
                d[1][0] * ex + d[1][1] * ey,
                d[2][2] * gxy,
            );
            let pz = match self.mode {
                PlaneMode::PlaneStress => 0.0,
                PlaneMode::PlaneStrain => nu * (px + py) - e * alpha_dt,
            };
            let v = (0.5 * ((px - py).powi(2) + (py - pz).powi(2) + (pz - px).powi(2))
                + 3.0 * pxy * pxy)
                .sqrt();
            sx.push(px);
            sy.push(py);
            txy.push(pxy);
            vm.push(v);
        }
        let disp: Vec<Vec3> = u.iter().map(|d| Vec3::new(d.x, d.y, 0.0)).collect();
        let max_u = disp.iter().map(|d| d.norm()).fold(0.0, f64::max);
        let max_vm = vm.iter().copied().fold(0.0, f64::max);
        write_field(
            out,
            m,
            "elasticity.vtu",
            &VtuMesh::from_triangles(s.mesh())?,
            vec![DataArray::vector("displacement_m", disp)],
            vec![
                DataArray::scalar("stress_xx_Pa", sx),
                DataArray::scalar("stress_yy_Pa", sy),
                DataArray::scalar("stress_xy_Pa", txy),
                DataArray::scalar("von_mises_Pa", vm),
            ],
        )?;
        metric(m, "max_displacement_m", max_u);
        metric(m, "max_von_mises_Pa", max_vm);
        m.linear_solver_converged = Some(true);
        Ok(())
    }
}

// ---------------------------------------------------------------- bar modes

fn three() -> usize {
    3
}

/// Fixed–free axial bar eigenmodes by lumped/consistent FEM, compared with `(2k-1)πc/2L`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BarModes1dCase {
    pub length_m: f64,
    pub area_m2: f64,
    pub youngs_modulus_pa: f64,
    pub density_kg_m3: f64,
    pub elements: usize,
    #[serde(default = "three")]
    pub modes: usize,
}

impl BarModes1dCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.elements >= 2
                && self.elements <= 2000
                && self.modes >= 1
                && self.modes <= 10
                && self.modes < self.elements,
            "bar needs 2..=2000 elements and 1..=10 modes",
        )?;
        let modes = bar_modal1d(
            BarModal1DConfig {
                elements: self.elements,
                length: self.length_m,
                area: self.area_m2,
                young_modulus: self.youngs_modulus_pa,
                density: self.density_kg_m3,
                ..Default::default()
            },
            self.modes,
        )?;
        let c = (self.youngs_modulus_pa / self.density_kg_m3).sqrt();
        let k: Vec<f64> = (1..=modes.len()).map(|k| k as f64).collect();
        let f: Vec<f64> = modes
            .iter()
            .map(|md| md.angular_frequency / (2.0 * std::f64::consts::PI))
            .collect();
        let exact: Vec<f64> = k
            .iter()
            .map(|k| (2.0 * k - 1.0) * c / (4.0 * self.length_m))
            .collect();
        let err: Vec<f64> = f.iter().zip(&exact).map(|(a, b)| a / b - 1.0).collect();
        csv(
            out,
            m,
            "modes.csv",
            &[
                ("mode", &k),
                ("frequency_Hz", &f),
                ("analytic_Hz", &exact),
                ("relative_error", &err),
            ],
        )?;
        let x = linspace(0.0, self.length_m, self.elements + 1);
        let shapes: Vec<(String, Vec<f64>)> = modes
            .iter()
            .enumerate()
            .map(|(i, md)| (format!("mode_{}", i + 1), md.displacement.clone()))
            .collect();
        require(
            shapes.iter().all(|(_, v)| v.len() == x.len()),
            "unexpected mode-shape size",
        )?;
        let mut table: Vec<(&str, &[f64])> = vec![("x_m", &x)];
        table.extend(shapes.iter().map(|(n, v)| (n.as_str(), v.as_slice())));
        csv(out, m, "mode_shapes.csv", &table)?;
        for (i, v) in f.iter().enumerate() {
            metric(m, &format!("mode_{}_Hz", i + 1), *v);
        }
        m.linear_solver_converged = Some(modes.iter().all(|md| md.eigenpair.converged));
        Ok(())
    }
}

// ---------------------------------------------------------------- J2 plasticity

fn tangent_step() -> f64 {
    1e-7
}

/// Small-strain J2 plasticity with linear isotropic hardening driven through cyclic uniaxial
/// stress: the axial strain follows the path and the lateral strains are iterated (Newton
/// on the algorithmic tangent) until the lateral stresses vanish.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct J2UniaxialCase {
    pub youngs_modulus_pa: f64,
    #[serde(default = "poisson")]
    pub poisson_ratio: f64,
    pub yield_stress_pa: f64,
    pub hardening_modulus_pa: f64,
    /// Turning points of the axial strain, starting from zero.
    pub strain_path: Vec<f64>,
    pub increments_per_segment: usize,
}

impl J2UniaxialCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            !self.strain_path.is_empty() && self.strain_path.len() <= 64,
            "use 1..=64 strain turning points",
        )?;
        require(
            self.increments_per_segment >= 1
                && self.increments_per_segment * self.strain_path.len() <= 20_000,
            "too many strain increments",
        )?;
        let model = J2Plasticity3D::new(J2Plasticity3DConfig {
            youngs_modulus: self.youngs_modulus_pa,
            poisson_ratio: self.poisson_ratio,
            yield_stress: self.yield_stress_pa,
            hardening_modulus: self.hardening_modulus_pa,
            tangent_step: tangent_step(),
        })?;
        let mut state = J2Plasticity3DState::default();
        let mut lateral = 0.0;
        let (mut strain, mut stress, mut peeq, mut lat) =
            (vec![0.0], vec![0.0], vec![0.0], vec![0.0]);
        let mut previous = 0.0;
        let mut worst_residual: f64 = 0.0;
        for &target in &self.strain_path {
            require(
                target.is_finite() && target.abs() < 1.0,
                "strain turning points must be finite and below 1",
            )?;
            for k in 1..=self.increments_per_segment {
                let e11 =
                    previous + (target - previous) * k as f64 / self.increments_per_segment as f64;
                let mut r = model.update(&[e11, lateral, lateral, 0.0, 0.0, 0.0], &state);
                for _ in 0..30 {
                    let s22 = r.stress[1];
                    if s22.abs() <= 1e-10 * self.yield_stress_pa {
                        break;
                    }
                    let t = &r.algorithmic_tangent;
                    let slope = t[6 + 1] + t[6 + 2];
                    require(slope.abs() > 0.0, "singular lateral tangent")?;
                    lateral -= s22 / slope;
                    r = model.update(&[e11, lateral, lateral, 0.0, 0.0, 0.0], &state);
                }
                worst_residual = worst_residual.max(r.stress[1].abs());
                state = r.state;
                strain.push(e11);
                stress.push(r.stress[0]);
                peeq.push(state.equivalent_plastic_strain);
                lat.push(lateral);
            }
            previous = target;
        }
        csv(
            out,
            m,
            "stress_strain.csv",
            &[
                ("axial_strain", &strain),
                ("axial_stress_Pa", &stress),
                ("equivalent_plastic_strain", &peeq),
                ("lateral_strain", &lat),
            ],
        )?;
        metric(
            m,
            "max_axial_stress_Pa",
            stress.iter().copied().fold(f64::MIN, f64::max),
        );
        metric(
            m,
            "min_axial_stress_Pa",
            stress.iter().copied().fold(f64::MAX, f64::min),
        );
        metric(
            m,
            "final_equivalent_plastic_strain",
            state.equivalent_plastic_strain,
        );
        metric(m, "max_lateral_stress_residual_Pa", worst_residual);
        m.linear_solver_converged = Some(worst_residual <= 1e-6 * self.yield_stress_pa);
        Ok(())
    }
}

// ---------------------------------------------------------------- phase-field fracture

/// Notched 1-D bar under monotonically increasing end displacement with an AT2 phase-field
/// damage model and history irreversibility.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhaseFieldFracture1dCase {
    pub cells: usize,
    pub length_m: f64,
    #[serde(default = "unit")]
    pub area_m2: f64,
    pub youngs_modulus_pa: f64,
    pub fracture_energy_j_m2: f64,
    pub length_scale_m: f64,
    #[serde(default)]
    pub notch_cell: Option<usize>,
    #[serde(default)]
    pub notch_damage: f64,
    pub max_end_displacement_m: f64,
    pub load_steps: usize,
}

fn unit() -> f64 {
    1.0
}

impl PhaseFieldFracture1dCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.cells >= 4
                && self.cells <= 4000
                && self.load_steps >= 1
                && self.load_steps <= 5000,
            "use 4..=4000 cells and 1..=5000 load steps",
        )?;
        let mut pf = PhaseFieldFracture1D::new(PhaseFieldFracture1DConfig {
            cells: self.cells,
            length: self.length_m,
            area: self.area_m2,
            youngs_modulus: self.youngs_modulus_pa,
            fracture_energy: self.fracture_energy_j_m2,
            length_scale: self.length_scale_m,
            notch_cell: self.notch_cell,
            initial_notch_damage: self.notch_damage,
            ..Default::default()
        })?;
        let (mut u, mut f, mut d, mut el, mut fr) = (vec![], vec![], vec![], vec![], vec![]);
        for k in 1..=self.load_steps {
            let disp = self.max_end_displacement_m * k as f64 / self.load_steps as f64;
            let r = pf.step(disp)?;
            u.push(disp);
            f.push(r.reaction_force);
            d.push(r.maximum_damage);
            el.push(r.elastic_energy);
            fr.push(r.fracture_energy);
        }
        csv(
            out,
            m,
            "load_displacement.csv",
            &[
                ("end_displacement_m", &u),
                ("reaction_force_N", &f),
                ("maximum_damage", &d),
                ("elastic_energy_J", &el),
                ("fracture_energy_J", &fr),
            ],
        )?;
        let x: Vec<f64> = (0..self.cells)
            .map(|i| (i as f64 + 0.5) * self.length_m / self.cells as f64)
            .collect();
        csv(
            out,
            m,
            "damage_profile.csv",
            &[("x_m", &x), ("damage", pf.damage())],
        )?;
        if let Some((k, peak)) = f.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)) {
            metric(m, "peak_force_N", *peak);
            metric(m, "displacement_at_peak_m", u[k]);
        }
        metric(m, "final_maximum_damage", d.last().copied().unwrap_or(0.0));
        Ok(())
    }
}

// ---------------------------------------------------------------- Cosserat rod

/// Straight rod along +z (the integrator's rest director) clamped over its first segment and
/// loaded by gravity and an optional tip force, by the position-based Cosserat integrator.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CosseratCantileverCase {
    pub length_m: f64,
    pub nodes: usize,
    pub dt_s: f64,
    pub steps: usize,
    #[serde(default)]
    pub gravity_m_s2: [f64; 3],
    #[serde(default)]
    pub tip_force_n: [f64; 3],
    /// XPBD bend/twist compliance (default 1e-4 here). The integrator is position based;
    /// compliance and node masses are not calibrated to a beam's EI.
    #[serde(default)]
    pub bend_twist_compliance: Option<f64>,
    #[serde(default)]
    pub damping: Option<f64>,
    #[serde(default)]
    pub output_interval: Option<usize>,
}

impl CosseratCantileverCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.nodes >= 3 && self.nodes <= 400 && self.steps >= 1 && self.steps <= 20_000,
            "use 3..=400 nodes and 1..=20000 steps",
        )?;
        let g = self.gravity_m_s2;
        let mut cfg = CosseratRodConfig {
            dt: self.dt_s,
            gravity: Vec3::new(g[0], g[1], g[2]),
            ..Default::default()
        };
        cfg.bend_twist_compliance = self.bend_twist_compliance.unwrap_or(1e-4);
        if let Some(c) = self.damping {
            cfg.linear_damping = c;
            cfg.angular_damping = c;
        }
        let mut rod = StableCosseratRod::new(cfg)?;
        let pos: Vec<Vec3> = (0..self.nodes)
            .map(|i| Vec3::new(0.0, 0.0, self.length_m * i as f64 / (self.nodes - 1) as f64))
            .collect();
        rod.initialize(&pos, &[])?;
        // Fixing the first two nodes clamps the root segment's direction as well as its position.
        rod.set_fixed(0, true)?;
        rod.set_fixed(1, true)?;
        let tip = self.nodes - 1;
        let interval = self
            .output_interval
            .unwrap_or(self.steps.div_ceil(200).max(1))
            .max(1);
        let (mut t, mut tx, mut ty, mut tz, mut ke) = (vec![], vec![], vec![], vec![], vec![]);
        for step in 1..=self.steps {
            rod.clear_external_loads();
            if self.tip_force_n != [0.0; 3] {
                let f = self.tip_force_n;
                rod.add_force(tip, Vec3::new(f[0], f[1], f[2]))?;
            }
            rod.step()?;
            if step % interval == 0 || step == self.steps {
                let p = rod.nodes()[tip].position;
                t.push(rod.time());
                tx.push(p.x);
                ty.push(p.y);
                tz.push(p.z);
                ke.push(rod.diagnostics().kinetic_energy);
            }
        }
        require(
            rod.diagnostics().finite,
            "Cosserat rod state became non-finite",
        )?;
        csv(
            out,
            m,
            "tip_history.csv",
            &[
                ("time_s", &t),
                ("tip_x_m", &tx),
                ("tip_y_m", &ty),
                ("tip_z_m", &tz),
                ("kinetic_energy_J", &ke),
            ],
        )?;
        let (x, y, z): (Vec<f64>, Vec<f64>, Vec<f64>) = (
            rod.nodes().iter().map(|n| n.position.x).collect(),
            rod.nodes().iter().map(|n| n.position.y).collect(),
            rod.nodes().iter().map(|n| n.position.z).collect(),
        );
        csv(
            out,
            m,
            "rod_shape.csv",
            &[("x_m", &x), ("y_m", &y), ("z_m", &z)],
        )?;
        let d = *rod.diagnostics();
        let tipp = rod.nodes()[tip].position;
        metric(
            m,
            "tip_deflection_m",
            (tipp.x * tipp.x + tipp.y * tipp.y + (tipp.z - self.length_m).powi(2)).sqrt(),
        );
        metric(m, "tip_x_m", tipp.x);
        metric(m, "tip_y_m", tipp.y);
        metric(m, "tip_z_m", tipp.z);
        metric(m, "kinetic_energy_J", d.kinetic_energy);
        metric(m, "elastic_energy_J", d.elastic_energy);
        metric(m, "max_stretch_shear_error", d.max_stretch_shear_error);
        Ok(())
    }
}

// ---------------------------------------------------------------- DEM

fn friction() -> f64 {
    0.5
}

/// Spheres in a regular block above a floor inside four side walls, settling under gravity
/// with Hertz–Mindlin contact. Writes a VTU point cloud per frame and a PVD collection.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DemGranularCase {
    pub radius_m: f64,
    pub density_kg_m3: f64,
    /// Particles per axis in the initial block.
    pub count: [usize; 3],
    /// Side length of the square container footprint.
    pub container_width_m: f64,
    pub dt_s: f64,
    pub steps: usize,
    pub frames: usize,
    #[serde(default = "friction")]
    pub friction: f64,
    #[serde(default)]
    pub normal_stiffness: Option<f64>,
}

impl DemGranularCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        let n = self.count[0] * self.count[1] * self.count[2];
        require((1..=4000).contains(&n), "use 1..=4000 particles")?;
        require(
            self.steps >= 1 && self.steps <= 20_000 && self.frames >= 1 && self.frames <= 100,
            "use 1..=20000 steps and 1..=100 frames",
        )?;
        let r = self.radius_m;
        let w = self.container_width_m;
        require(
            r > 0.0 && w > 2.0 * r * self.count[0].max(self.count[2]) as f64,
            "the particle block must fit inside the container",
        )?;
        let mut contact = HertzMindlinContactModel {
            friction: self.friction,
            ..Default::default()
        };
        if let Some(k) = self.normal_stiffness {
            contact.normal_stiffness = k;
        }
        let mut dem = ExplicitDemSystem::new(DemConfig {
            dt: self.dt_s,
            contact,
            ..Default::default()
        })?;
        dem.add_plane(PlaneShape {
            normal: Vec3::new(0.0, 1.0, 0.0),
            offset: 0.0,
        });
        dem.add_plane(PlaneShape {
            normal: Vec3::new(1.0, 0.0, 0.0),
            offset: 0.0,
        });
        dem.add_plane(PlaneShape {
            normal: Vec3::new(-1.0, 0.0, 0.0),
            offset: -w,
        });
        dem.add_plane(PlaneShape {
            normal: Vec3::new(0.0, 0.0, 1.0),
            offset: 0.0,
        });
        dem.add_plane(PlaneShape {
            normal: Vec3::new(0.0, 0.0, -1.0),
            offset: -w,
        });
        let pitch = 2.2 * r;
        let x0 = 0.5 * (w - pitch * (self.count[0] - 1) as f64);
        let z0 = 0.5 * (w - pitch * (self.count[2] - 1) as f64);
        for j in 0..self.count[1] {
            for k in 0..self.count[2] {
                for i in 0..self.count[0] {
                    // A small deterministic jitter breaks the perfect stacking symmetry.
                    let jitter = 0.05 * r * (((i * 7 + j * 13 + k * 31) % 11) as f64 / 10.0 - 0.5);
                    dem.add_particle(
                        Vec3::new(
                            x0 + pitch * i as f64 + jitter,
                            1.5 * r + pitch * j as f64,
                            z0 + pitch * k as f64 - jitter,
                        ),
                        r,
                        self.density_kg_m3,
                        Vec3::default(),
                    )?;
                }
            }
        }
        let interval = (self.steps / self.frames).max(1);
        let mut frames = vec![];
        let (mut t, mut ke, mut h) = (vec![], vec![], vec![]);
        let write = |dem: &ExplicitDemSystem, file: &str, m: &mut RunManifest| -> Result<()> {
            let points: Vec<Vec3> = dem.bodies().iter().map(|b| b.state.position).collect();
            let cells = (0..points.len())
                .map(|i| Cell {
                    cell_type: CellType::Vertex,
                    nodes: vec![i],
                })
                .collect();
            let v: Vec<Vec3> = dem
                .bodies()
                .iter()
                .map(|b| b.state.linear_velocity)
                .collect();
            let radius = dem.spheres().iter().map(|s| s.radius).collect();
            write_field(
                out,
                m,
                file,
                &VtuMesh { points, cells },
                vec![
                    DataArray::vector("velocity_m_per_s", v),
                    DataArray::scalar("radius_m", radius),
                ],
                vec![],
            )
        };
        write(&dem, "particles_0000.vtu", m)?;
        frames.push((0.0, "particles_0000.vtu".to_string()));
        for step in 1..=self.steps {
            dem.step()?;
            if step % interval == 0 || step == self.steps {
                let k: f64 = dem
                    .bodies()
                    .iter()
                    .map(|b| 0.5 * b.mass * b.state.linear_velocity.norm().powi(2))
                    .sum();
                let mean_y =
                    dem.bodies().iter().map(|b| b.state.position.y).sum::<f64>() / n as f64;
                require(
                    k.is_finite() && mean_y.is_finite(),
                    "DEM state became non-finite",
                )?;
                t.push(dem.time());
                ke.push(k);
                h.push(mean_y);
                if frames.len() <= self.frames {
                    let file = format!("particles_{:04}.vtu", frames.len());
                    write(&dem, &file, m)?;
                    frames.push((dem.time(), file));
                }
            }
        }
        io::write_pvd(out.join("particles.pvd"), &frames)?;
        m.files.push("particles.pvd".into());
        csv(
            out,
            m,
            "settling.csv",
            &[
                ("time_s", &t),
                ("kinetic_energy_J", &ke),
                ("mean_height_m", &h),
            ],
        )?;
        let escaped = dem
            .bodies()
            .iter()
            .filter(|b| {
                b.state.position.y < -r || b.state.position.x < -r || b.state.position.x > w + r
            })
            .count();
        metric(m, "particles", n as f64);
        metric(
            m,
            "final_kinetic_energy_J",
            ke.last().copied().unwrap_or(0.0),
        );
        metric(m, "final_mean_height_m", h.last().copied().unwrap_or(0.0));
        metric(m, "escaped_particles", escaped as f64);
        json(
            out,
            m,
            "dem_summary.json",
            &json!({"particles": n, "frames": frames.len(), "contact": "hertz_mindlin", "friction": self.friction}),
        )?;
        if escaped > 0 {
            return Err(SolverError::NumericalFailure(format!(
                "{escaped} particles left the container; reduce dt_s"
            )));
        }
        Ok(())
    }
}
