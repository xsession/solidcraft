//! Particle, kinetic, coupled and microstructure case adapters: 1-D electrostatic PIC
//! (two-stream), condensed-history slab transport with depth dose, Pennes bioheat with an SAR
//! hot spot, 2-D Joule heating, Cahn–Hilliard spinodal decomposition, DG advection,
//! Vlasov–Poisson Landau damping, Lennard-Jones molecular dynamics, photonic waveguide
//! modes and a coupled THMC column.
use super::{csv, json, linspace, metric};
use crate::{
    atomistics::{
        velocity_rescaling_thermostat, LennardJones1DForceField, LennardJonesParameters,
        NeighborList1DConfig, PeriodicNeighborList1D,
    },
    case::{write_field, RunManifest},
    core::{LinearConfig, Vec3},
    dg::{DgAdvection1D, DgAdvection1DConfig},
    em::{solve_photonic_waveguide_modes, PhotonicGrid2D, PhotonicWaveguideConfig},
    fem::{make_rectangle_tri_mesh, Heat2DConfig, ScalarBoundaryCondition2D},
    io::{DataArray, VtuMesh},
    kinetic::{VlasovPoisson1D1V, VlasovPoisson1D1VConfig},
    multiphysics::{
        BioheatBoundary, JouleHeatingCoupler2D, PennesBioheat2D, PennesBioheat2DConfig,
    },
    particle::{
        collect_transport_hits, make_transport_dose_grid, score_hits_to_dose_grid,
        transport_track_stochastic, AxisAlignedParticleBox, ElectrostaticPic1D,
        ElectrostaticPic1DConfig, PicParticle1D, TransportConfig, TransportMaterial,
        TransportParticleDefinition, TransportProcess, TransportProcessKind, TransportRandom,
        TransportRegion, TransportSensitiveDetector, TransportTrack, TransportWorld, Vec3 as PVec3,
    },
    porous::{CahnHilliard1D, CahnHilliard1DConfig, Thmc1D, Thmc1DConfig},
    require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::f64::consts::PI;
use std::path::Path;

const EPS0: f64 = 8.8541878128e-12;

fn one() -> f64 {
    1.0
}

/// Least-squares slope of `ln(y)` against `t` over the samples where `lo <= y <= hi`.
fn log_slope(t: &[f64], y: &[f64], lo: f64, hi: f64) -> Option<f64> {
    let pts: Vec<(f64, f64)> = t
        .iter()
        .zip(y)
        .filter(|(_, v)| **v >= lo && **v <= hi && **v > 0.0)
        .map(|(a, b)| (*a, b.ln()))
        .collect();
    if pts.len() < 3 {
        return None;
    }
    let n = pts.len() as f64;
    let (sx, sy) = pts.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0, a.1 + p.1));
    let (mx, my) = (sx / n, sy / n);
    let (num, den) = pts.iter().fold((0.0, 0.0), |a, p| {
        (a.0 + (p.0 - mx) * (p.1 - my), a.1 + (p.0 - mx).powi(2))
    });
    (den > 0.0).then(|| num / den)
}

// ---------------------------------------------------------------- two-stream PIC

/// Two cold, counter-streaming electron beams in a periodic box with a neutralising
/// background. Particles have unit mass and charge −1 C; macro-weights are chosen so the
/// total plasma frequency equals `plasma_frequency_rad_s`. The seeded mode grows at the
/// cold-beam rate until the beams trap.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PicTwoStreamCase {
    pub length_m: f64,
    pub grid_points: usize,
    pub particles_per_beam: usize,
    pub drift_velocity_m_s: f64,
    #[serde(default = "one")]
    pub plasma_frequency_rad_s: f64,
    pub dt_s: f64,
    pub steps: usize,
    #[serde(default = "seed_amp")]
    pub perturbation: f64,
    #[serde(default = "mode1")]
    pub mode: usize,
}

fn seed_amp() -> f64 {
    1e-3
}
fn mode1() -> usize {
    1
}

impl PicTwoStreamCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.particles_per_beam >= 8 && self.particles_per_beam <= 50_000,
            "use 8..=50000 particles per beam",
        )?;
        require(
            self.steps >= 1
                && self.steps <= 20_000
                && self.grid_points >= 8
                && self.grid_points <= 4096,
            "use 8..=4096 grid points and 1..=20000 steps",
        )?;
        require(
            self.particles_per_beam as f64 * self.steps as f64 <= 2e8,
            "particles times steps exceeds 2e8",
        )?;
        let l = self.length_m;
        let total = 2 * self.particles_per_beam;
        let weight = EPS0 * l * self.plasma_frequency_rad_s.powi(2) / total as f64;
        let k = 2.0 * PI * self.mode as f64 / l;
        let mut particles = Vec::with_capacity(total);
        for beam in [-1.0, 1.0] {
            for i in 0..self.particles_per_beam {
                let x0 = l * (i as f64 + 0.5) / self.particles_per_beam as f64;
                let x = x0 + self.perturbation * l * (k * x0).sin();
                particles.push(PicParticle1D {
                    position_m: x,
                    velocity_m_per_s: beam * self.drift_velocity_m_s,
                    charge_c: -1.0,
                    mass_kg: 1.0,
                    weight,
                });
            }
        }
        let mut pic = ElectrostaticPic1D::new(ElectrostaticPic1DConfig {
            grid_points: self.grid_points,
            length_m: l,
            dt_s: self.dt_s,
            neutralize_mean_charge: true,
            ..Default::default()
        })?;
        pic.set_particles(particles)?;
        pic.deposit_and_solve()?;
        let dx = l / self.grid_points as f64;
        let field_energy = |p: &ElectrostaticPic1D| {
            0.5 * EPS0 * p.electric_field().iter().map(|e| e * e).sum::<f64>() * dx
        };
        let kinetic = |p: &ElectrostaticPic1D| {
            p.particles()
                .iter()
                .map(|q| 0.5 * q.mass_kg * q.weight * q.velocity_m_per_s.powi(2))
                .sum::<f64>()
        };
        let (mut t, mut we, mut wk) = (vec![0.0], vec![field_energy(&pic)], vec![kinetic(&pic)]);
        for _ in 0..self.steps {
            pic.step(1)?;
            t.push(pic.time_s());
            we.push(field_energy(&pic));
            wk.push(kinetic(&pic));
        }
        require(
            we.iter().chain(&wk).all(|v| v.is_finite()),
            "PIC energies became non-finite",
        )?;
        let total_e: Vec<f64> = we.iter().zip(&wk).map(|(a, b)| a + b).collect();
        csv(
            out,
            m,
            "energy.csv",
            &[
                ("time_s", &t),
                ("field_energy_J", &we),
                ("kinetic_energy_J", &wk),
                ("total_energy_J", &total_e),
            ],
        )?;
        let x: Vec<f64> = pic.particles().iter().map(|p| p.position_m).collect();
        let v: Vec<f64> = pic.particles().iter().map(|p| p.velocity_m_per_s).collect();
        csv(out, m, "phase_space.csv", &[("x_m", &x), ("v_m_per_s", &v)])?;
        let grid_x = linspace(0.0, l - dx, self.grid_points);
        csv(
            out,
            m,
            "fields.csv",
            &[
                ("x_m", &grid_x),
                ("charge_density_C_m3", pic.charge_density()),
                ("electric_field_V_m", pic.electric_field()),
            ],
        )?;
        let wb2 = 0.5 * self.plasma_frequency_rad_s.powi(2);
        let kv2 = (k * self.drift_velocity_m_s).powi(2);
        let unstable = kv2 + wb2 - (wb2 * (4.0 * kv2 + wb2)).sqrt();
        let theory = if unstable < 0.0 {
            (-unstable).sqrt()
        } else {
            0.0
        };
        metric(m, "theory_growth_rate_per_s", theory);
        let (ipeak, peak) =
            we.iter()
                .copied()
                .enumerate()
                .fold((0, 0.0), |a, (i, v)| if v > a.1 { (i, v) } else { a });
        // Field energy grows at twice the amplitude rate; fit the linear phase before saturation.
        if let Some(s) = log_slope(&t[..=ipeak], &we[..=ipeak], 1e-3 * peak, 0.1 * peak) {
            metric(m, "measured_growth_rate_per_s", 0.5 * s);
        }
        metric(m, "peak_field_energy_J", peak);
        metric(
            m,
            "relative_energy_change",
            (total_e[total_e.len() - 1] - total_e[0]) / total_e[0].abs().max(f64::MIN_POSITIVE),
        );
        Ok(())
    }
}

// ---------------------------------------------------------------- slab transport

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SlabLayer {
    pub name: String,
    pub thickness_m: f64,
    pub density_kg_m3: f64,
    pub stopping_power_ev_per_m: f64,
    /// Mean free path of the discrete interaction; zero disables it.
    #[serde(default)]
    pub interaction_length_m: f64,
    #[serde(default)]
    pub energy_loss_fraction: f64,
    #[serde(default)]
    pub secondary_energy_fraction: f64,
}

/// Pencil beam normally incident on stacked slab layers (each 1 m × 1 m laterally),
/// condensed-history stepping with sampled discrete interactions, depth-dose scoring.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParticleTransportCase {
    pub layers: Vec<SlabLayer>,
    pub beam_energy_ev: f64,
    pub histories: usize,
    pub depth_bins: usize,
    #[serde(default = "max_step")]
    pub max_step_m: f64,
    #[serde(default)]
    pub energy_cut_ev: f64,
    #[serde(default = "default_seed")]
    pub seed: u64,
}

fn max_step() -> f64 {
    1e-3
}
fn default_seed() -> u64 {
    12345
}

impl ParticleTransportCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            !self.layers.is_empty() && self.layers.len() <= 16,
            "use 1..=16 layers",
        )?;
        require(
            self.histories >= 1
                && self.histories <= 100_000
                && self.depth_bins >= 1
                && self.depth_bins <= 2000,
            "use 1..=100000 histories and 1..=2000 depth bins",
        )?;
        let mut world = TransportWorld {
            particles: vec![TransportParticleDefinition {
                name: "beam".into(),
                ..Default::default()
            }],
            sensitive_detectors: vec![TransportSensitiveDetector {
                name: "slab".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let mut x = 0.0;
        for (i, layer) in self.layers.iter().enumerate() {
            require(layer.thickness_m > 0.0, "layer thickness must be positive")?;
            world.materials.push(TransportMaterial {
                name: layer.name.clone(),
                density_kg_per_m3: layer.density_kg_m3,
                production_cut_energy_ev: 0.0,
            });
            world.regions.push(TransportRegion {
                name: layer.name.clone(),
                r#box: AxisAlignedParticleBox {
                    minimum_m: PVec3::new(x, 0.0, 0.0),
                    maximum_m: PVec3::new(x + layer.thickness_m, 1.0, 1.0),
                },
                material_index: i,
                sensitive_detector_index: Some(0),
            });
            world.processes.push(TransportProcess {
                name: format!("{}_dedx", layer.name),
                stopping_power_ev_per_m: layer.stopping_power_ev_per_m,
                material_index: Some(i),
                ..Default::default()
            });
            if layer.interaction_length_m > 0.0 {
                world.processes.push(TransportProcess {
                    name: format!("{}_hard", layer.name),
                    kind: TransportProcessKind::DiscreteInteraction,
                    material_index: Some(i),
                    physical_interaction_length_m: layer.interaction_length_m,
                    energy_loss_fraction: layer.energy_loss_fraction,
                    secondary_energy_fraction: layer.secondary_energy_fraction,
                    ..Default::default()
                });
            }
            x += layer.thickness_m;
        }
        let depth = x;
        let cfg = TransportConfig {
            max_step_m: self.max_step_m,
            energy_cut_ev: self.energy_cut_ev,
            ..Default::default()
        };
        let mut rng = TransportRandom { state: self.seed };
        let mut grid = make_transport_dose_grid(
            self.depth_bins,
            1,
            1,
            AxisAlignedParticleBox {
                minimum_m: PVec3::default(),
                maximum_m: PVec3::new(depth, 1.0, 1.0),
            },
            1000.0,
        )?;
        // Per-bin mass uses each bin's layer density; the grid's uniform density is replaced below.
        let (mut ranges, mut transmitted, mut secondaries, mut secondary_energy, mut deposited) =
            (vec![], 0usize, 0usize, 0.0, 0.0);
        let mut energy_bins = vec![0.0; self.depth_bins];
        for _ in 0..self.histories {
            let track = TransportTrack {
                position_m: PVec3::new(0.0, 0.5, 0.5),
                direction: PVec3::new(1.0, 0.0, 0.0),
                kinetic_energy_ev: self.beam_energy_ev,
                ..Default::default()
            };
            let r = transport_track_stochastic(&world, track, &cfg, &mut rng, 100_000)?;
            let hits = collect_transport_hits(&world, &r.steps)?;
            grid.energy_deposit_ev.iter_mut().for_each(|e| *e = 0.0);
            score_hits_to_dose_grid(&hits, &mut grid)?;
            for (a, b) in energy_bins.iter_mut().zip(&grid.energy_deposit_ev) {
                *a += b;
            }
            deposited += r.scoring.total_energy_deposit_ev;
            secondaries += r.secondaries.len();
            secondary_energy += r
                .secondaries
                .iter()
                .map(|s| s.track.kinetic_energy_ev)
                .sum::<f64>();
            if r.primary.position_m.x >= depth - 1e-9 {
                transmitted += 1;
            }
            ranges.push(r.primary.position_m.x.min(depth));
        }
        let dz = depth / self.depth_bins as f64;
        let z: Vec<f64> = (0..self.depth_bins)
            .map(|i| (i as f64 + 0.5) * dz)
            .collect();
        let density: Vec<f64> = z
            .iter()
            .map(|zc| {
                let mut edge = 0.0;
                for l in &self.layers {
                    edge += l.thickness_m;
                    if *zc < edge {
                        return l.density_kg_m3;
                    }
                }
                self.layers.last().map_or(1000.0, |l| l.density_kg_m3)
            })
            .collect();
        let per_history: Vec<f64> = energy_bins
            .iter()
            .map(|e| e / self.histories as f64)
            .collect();
        // Dose per history for a 1 m² beam footprint: E / (rho · dz · 1 m²).
        let dose: Vec<f64> = per_history
            .iter()
            .zip(&density)
            .map(|(e, rho)| e * 1.602176634e-19 / (rho.max(f64::MIN_POSITIVE) * dz))
            .collect();
        csv(
            out,
            m,
            "depth_dose.csv",
            &[
                ("depth_m", &z),
                ("energy_deposit_eV_per_history", &per_history),
                ("dose_Gy_per_history", &dose),
            ],
        )?;
        write_field(
            out,
            m,
            "depth_dose.vtu",
            &VtuMesh::box_cells([self.depth_bins, 1, 1], Vec3::new(depth, 1.0, 1.0), 0.0)?,
            vec![],
            vec![
                DataArray::scalar("dose_Gy_per_history", dose.clone()),
                DataArray::scalar("density_kg_m3", density),
            ],
        )?;
        let mean_range = ranges.iter().sum::<f64>() / ranges.len() as f64;
        metric(m, "mean_primary_range_m", mean_range);
        metric(
            m,
            "transmitted_fraction",
            transmitted as f64 / self.histories as f64,
        );
        metric(
            m,
            "deposited_fraction",
            deposited / (self.beam_energy_ev * self.histories as f64),
        );
        metric(
            m,
            "secondaries_per_history",
            secondaries as f64 / self.histories as f64,
        );
        metric(
            m,
            "secondary_energy_fraction",
            secondary_energy / (self.beam_energy_ev * self.histories as f64),
        );
        metric(
            m,
            "peak_dose_depth_m",
            dose.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map_or(0.0, |(i, _)| z[i]),
        );
        json(
            out,
            m,
            "transport_summary.json",
            &json!({"histories": self.histories, "note": "Secondaries are created and counted but not transported; their energy leaves the primary track."}),
        )?;
        Ok(())
    }
}

// ---------------------------------------------------------------- bioheat

/// Pennes bioheat in a 2-D tissue section with a Gaussian SAR hot spot and fixed or periodic
/// edges; reports temperature, the CEM43 thermal dose and the hottest point history.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BioheatCase {
    pub nx: usize,
    pub ny: usize,
    pub width_m: f64,
    pub height_m: f64,
    pub peak_sar_w_kg: f64,
    pub hotspot_center_m: [f64; 2],
    pub hotspot_radius_m: f64,
    pub dt_s: f64,
    pub steps: usize,
    #[serde(default = "perfusion")]
    pub blood_perfusion_per_s: f64,
    #[serde(default)]
    pub periodic: bool,
}

fn perfusion() -> f64 {
    0.005
}

impl BioheatCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.nx * self.ny <= 16_384 && self.steps >= 1 && self.steps <= 20_000,
            "use at most 16384 cells and 1..=20000 steps",
        )?;
        let cfg = PennesBioheat2DConfig {
            nx: self.nx,
            ny: self.ny,
            width_m: self.width_m,
            height_m: self.height_m,
            dt_s: self.dt_s,
            blood_perfusion_per_s: self.blood_perfusion_per_s,
            boundary: if self.periodic {
                BioheatBoundary::Periodic
            } else {
                BioheatBoundary::FixedTemperature
            },
            ..Default::default()
        };
        let mut p = PennesBioheat2D::new(cfg)?;
        let (dx, dy) = (
            self.width_m / self.nx as f64,
            self.height_m / self.ny as f64,
        );
        let sar: Vec<f64> = (0..self.nx * self.ny)
            .map(|i| {
                let (x, y) = (
                    ((i % self.nx) as f64 + 0.5) * dx,
                    ((i / self.nx) as f64 + 0.5) * dy,
                );
                let r2 =
                    (x - self.hotspot_center_m[0]).powi(2) + (y - self.hotspot_center_m[1]).powi(2);
                self.peak_sar_w_kg * (-r2 / (self.hotspot_radius_m * self.hotspot_radius_m)).exp()
            })
            .collect();
        p.set_sar(&sar)?;
        let mut cem43 = vec![0.0; self.nx * self.ny];
        let (mut t, mut tmax) = (
            vec![0.0],
            vec![p.temperature_k().iter().copied().fold(f64::MIN, f64::max)],
        );
        for _ in 0..self.steps {
            p.step(1)?;
            for (d, temp) in cem43.iter_mut().zip(p.temperature_k()) {
                let c = temp - 273.15;
                let r: f64 = if c >= 43.0 { 0.5 } else { 0.25 };
                *d += r.powf(43.0 - c) * self.dt_s / 60.0;
            }
            t.push(p.time_s());
            tmax.push(p.temperature_k().iter().copied().fold(f64::MIN, f64::max));
        }
        csv(
            out,
            m,
            "hotspot_history.csv",
            &[("time_s", &t), ("max_temperature_K", &tmax)],
        )?;
        let celsius: Vec<f64> = p.temperature_k().iter().map(|k| k - 273.15).collect();
        write_field(
            out,
            m,
            "bioheat.vtu",
            &VtuMesh::rectangle_cells(self.nx, self.ny, dx, dy)?,
            vec![],
            vec![
                DataArray::scalar("temperature_K", p.temperature_k().to_vec()),
                DataArray::scalar("temperature_C", celsius),
                DataArray::scalar("sar_W_kg", sar),
                DataArray::scalar("cem43_min", cem43.clone()),
            ],
        )?;
        metric(m, "max_temperature_K", tmax.last().copied().unwrap_or(0.0));
        metric(
            m,
            "max_cem43_min",
            cem43.iter().copied().fold(0.0, f64::max),
        );
        metric(m, "time_s", p.time_s());
        if let Some(r) = p.linear_result() {
            m.linear_solver_converged = Some(r.converged);
        }
        Ok(())
    }
}

// ---------------------------------------------------------------- Joule heating

/// Rectangular conductor with a voltage across its left and right edges: steady current
/// flow, Joule heating and transient FEM heat conduction with fixed-temperature edges.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JouleHeatingCase {
    pub width_m: f64,
    pub height_m: f64,
    pub nx: usize,
    pub ny: usize,
    pub voltage_v: f64,
    pub electrical_conductivity_s_m: f64,
    pub thermal_conductivity_w_mk: f64,
    pub volumetric_heat_capacity_j_m3k: f64,
    #[serde(default = "ambient")]
    pub boundary_temperature_k: f64,
    pub dt_s: f64,
    pub steps: usize,
}

fn ambient() -> f64 {
    300.0
}

impl JouleHeatingCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.steps >= 1 && self.steps <= 20_000,
            "use 1..=20000 steps",
        )?;
        require(
            self.electrical_conductivity_s_m > 0.0,
            "electrical conductivity must be positive",
        )?;
        let mesh = make_rectangle_tri_mesh(self.nx, self.ny, self.width_m, self.height_m)?;
        let mut j = JouleHeatingCoupler2D::new(
            mesh,
            Heat2DConfig {
                conductivity: self.thermal_conductivity_w_mk,
                volumetric_heat_capacity: self.volumetric_heat_capacity_j_m3k,
                dt: self.dt_s,
                ..Default::default()
            },
            LinearConfig {
                max_iterations: 4000,
                relative_tolerance: 1e-12,
                restart: 30,
            },
        )?;
        let v = self.voltage_v;
        j.set_electrical_boundary(
            0,
            ScalarBoundaryCondition2D::Dirichlet(Box::new(move |_| v)),
        );
        j.set_electrical_boundary(1, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 0.0)));
        let tb = self.boundary_temperature_k;
        j.initialize_temperature(|_| tb)?;
        j.set_thermal_dirichlet(move |_, _| tb);
        let sigma = self.electrical_conductivity_s_m;
        j.solve_electrical(&|_| sigma, None)?;
        let (mut t, mut tmax) = (vec![0.0], vec![tb]);
        for _ in 0..self.steps {
            j.thermal_step()?;
            t.push(j.thermal().time);
            tmax.push(
                j.thermal()
                    .temperature
                    .iter()
                    .copied()
                    .fold(f64::MIN, f64::max),
            );
        }
        csv(
            out,
            m,
            "temperature_history.csv",
            &[("time_s", &t), ("max_temperature_K", &tmax)],
        )?;
        let current = j.electrical().element_current_density(&|_| sigma)?;
        let jvec: Vec<Vec3> = current.iter().map(|c| Vec3::new(c.x, c.y, 0.0)).collect();
        let q = j.joule_heating_density().to_vec();
        let mesh = j.thermal().mesh.clone();
        let total_power: f64 = mesh
            .triangles
            .iter()
            .zip(&q)
            .map(|(tri, qe)| {
                let [a, b, c] = tri.node;
                let (pa, pb, pc) = (mesh.nodes[a], mesh.nodes[b], mesh.nodes[c]);
                0.5 * ((pb.x - pa.x) * (pc.y - pa.y) - (pc.x - pa.x) * (pb.y - pa.y)).abs() * qe
            })
            .sum();
        write_field(
            out,
            m,
            "joule.vtu",
            &VtuMesh::from_triangles(&mesh)?,
            vec![
                DataArray::scalar("potential_V", j.electrical().potential().to_vec()),
                DataArray::scalar("temperature_K", j.thermal().temperature.clone()),
            ],
            vec![
                DataArray::vector("current_density_A_m2", jvec),
                DataArray::scalar("joule_heat_W_m3", q),
            ],
        )?;
        metric(m, "joule_power_W_per_m_depth", total_power);
        metric(
            m,
            "analytic_power_W_per_m_depth",
            self.voltage_v.powi(2) * sigma * self.height_m / self.width_m,
        );
        metric(m, "max_temperature_K", tmax.last().copied().unwrap_or(tb));
        if let Some(r) = j.electrical().linear_result() {
            m.linear_solver_converged = Some(r.converged);
        }
        Ok(())
    }
}

// ---------------------------------------------------------------- Cahn–Hilliard

/// Periodic 1-D Cahn–Hilliard spinodal decomposition from a perturbed mean composition.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CahnHilliardCase {
    pub cells: usize,
    #[serde(default = "one")]
    pub length: f64,
    #[serde(default = "one")]
    pub mobility: f64,
    pub gradient_energy: f64,
    pub mean_composition: f64,
    pub perturbation: f64,
    pub dt: f64,
    pub steps: usize,
    pub snapshots: usize,
}

impl CahnHilliardCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.cells >= 8
                && self.cells <= 4096
                && self.steps >= 1
                && self.steps <= 20_000
                && self.snapshots >= 1
                && self.snapshots <= 20,
            "use 8..=4096 cells, 1..=20000 steps and 1..=20 snapshots",
        )?;
        let mut ch = CahnHilliard1D::new(
            CahnHilliard1DConfig {
                cells: self.cells,
                length: self.length,
                mobility: self.mobility,
                gradient_energy: self.gradient_energy,
                ..Default::default()
            },
            None,
        )?;
        ch.initialize(self.mean_composition, self.perturbation)?;
        let f0 = ch.free_energy()?;
        let interval = (self.steps / self.snapshots).max(1);
        let x: Vec<f64> = (0..self.cells)
            .map(|i| (i as f64 + 0.5) * self.length / self.cells as f64)
            .collect();
        let mut columns: Vec<(String, Vec<f64>)> = vec![("t_0".into(), ch.composition().to_vec())];
        let (mut t, mut fe) = (vec![0.0], vec![f0]);
        for s in 1..=self.steps {
            ch.step(self.dt)?;
            if s % interval == 0 || s == self.steps {
                let f = ch.free_energy()?;
                let last = fe.last().copied().unwrap_or(f0);
                // The scheme is explicit: a rising free energy means the step is unstable.
                require(
                    f <= last + 1e-9 * last.abs().max(1e-12),
                    "free energy increased; the explicit step is unstable, reduce dt",
                )?;
                t.push(s as f64 * self.dt);
                fe.push(f);
                if columns.len() <= self.snapshots {
                    columns.push((format!("t_{s}"), ch.composition().to_vec()));
                }
            }
        }
        csv(
            out,
            m,
            "free_energy.csv",
            &[("time", &t), ("free_energy", &fe)],
        )?;
        let mut table: Vec<(&str, &[f64])> = vec![("x", &x)];
        table.extend(columns.iter().map(|(n, v)| (n.as_str(), v.as_slice())));
        csv(out, m, "composition.csv", &table)?;
        let c = ch.composition();
        metric(m, "mean_composition", ch.mean());
        metric(m, "free_energy_drop", f0 - fe.last().copied().unwrap_or(f0));
        metric(
            m,
            "composition_min",
            c.iter().copied().fold(f64::MAX, f64::min),
        );
        metric(
            m,
            "composition_max",
            c.iter().copied().fold(f64::MIN, f64::max),
        );
        Ok(())
    }
}

// ---------------------------------------------------------------- DG advection

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DgProfile {
    #[default]
    Sine,
    Gaussian,
    Square,
}

/// Periodic linear advection with nodal Gauss–Lobatto DG; compares with the exact shift.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DgAdvectionCase {
    pub elements: usize,
    pub order: usize,
    #[serde(default = "one")]
    pub length: f64,
    #[serde(default = "one")]
    pub speed: f64,
    #[serde(default)]
    pub profile: DgProfile,
    pub dt: f64,
    pub steps: usize,
}

impl DgAdvectionCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.elements >= 1
                && self.elements <= 2000
                && self.order >= 1
                && self.order <= 16
                && self.steps <= 20_000,
            "use 1..=2000 elements, order 1..=16 and at most 20000 steps",
        )?;
        let l = self.length;
        let profile = self.profile;
        let f = move |x: f64| {
            let s = x.rem_euclid(l) / l;
            match profile {
                DgProfile::Sine => (2.0 * PI * s).sin(),
                DgProfile::Gaussian => (-((s - 0.5) / 0.08).powi(2)).exp(),
                DgProfile::Square => {
                    if (0.3..0.7).contains(&s) {
                        1.0
                    } else {
                        0.0
                    }
                }
            }
        };
        let mut dg = DgAdvection1D::new(DgAdvection1DConfig {
            elements: self.elements,
            order: self.order,
            length: l,
            speed: self.speed,
        })?;
        dg.initialize(f);
        let m0 = dg.mass();
        for _ in 0..self.steps {
            dg.step(self.dt)?;
        }
        let p = dg.reference_nodes().to_vec();
        let h = l / self.elements as f64;
        let shift = self.speed * self.dt * self.steps as f64;
        let (mut x, mut u, mut exact) = (vec![], vec![], vec![]);
        for e in 0..self.elements {
            for (k, xi) in p.iter().enumerate() {
                let xx = h * (e as f64 + 0.5 * (xi + 1.0));
                x.push(xx);
                u.push(dg.values()[e * p.len() + k]);
                exact.push(f(xx - shift));
            }
        }
        require(
            u.iter().all(|v| v.is_finite()),
            "DG solution became non-finite; reduce dt",
        )?;
        let err: Vec<f64> = u.iter().zip(&exact).map(|(a, b)| a - b).collect();
        csv(
            out,
            m,
            "solution.csv",
            &[("x", &x), ("dg", &u), ("exact", &exact), ("error", &err)],
        )?;
        metric(
            m,
            "max_error",
            err.iter().map(|e| e.abs()).fold(0.0, f64::max),
        );
        metric(m, "mass_drift", dg.mass() - m0);
        metric(m, "time", shift / self.speed.abs().max(f64::MIN_POSITIVE));
        Ok(())
    }
}

// ---------------------------------------------------------------- Landau damping

/// Linear Landau damping of a cosine density perturbation on a Maxwellian, in normalised
/// units (ω_p = 1, λ_D = 1) by the positivity-preserving first-order upwind Vlasov–Poisson
/// reference. Its numerical diffusion adds to the physical damping; compare with the
/// reported asymptotic rate and refine both grids to see the difference shrink.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LandauDampingCase {
    pub spatial_cells: usize,
    pub velocity_cells: usize,
    pub wavenumber: f64,
    pub amplitude: f64,
    #[serde(default = "vmax")]
    pub maximum_velocity: f64,
    pub dt: f64,
    pub steps: usize,
}

fn vmax() -> f64 {
    6.0
}

impl LandauDampingCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.spatial_cells * self.velocity_cells <= 65_536
                && self.steps >= 1
                && self.steps <= 20_000,
            "use at most 65536 phase-space cells and 1..=20000 steps",
        )?;
        require(self.wavenumber > 0.0, "wavenumber must be positive")?;
        let (k, a) = (self.wavenumber, self.amplitude);
        let mut vp = VlasovPoisson1D1V::new(VlasovPoisson1D1VConfig {
            spatial_cells: self.spatial_cells,
            velocity_cells: self.velocity_cells,
            spatial_length: 2.0 * PI / k,
            maximum_velocity: self.maximum_velocity,
            neutralizing_density: 1.0,
        })?;
        vp.initialize(move |x, v| {
            (-0.5 * v * v).exp() / (2.0 * PI).sqrt() * (1.0 + a * (k * x).cos())
        })?;
        let m0 = vp.mass();
        let (mut t, mut w) = (vec![0.0], vec![vp.field_energy()]);
        for s in 1..=self.steps {
            vp.step(self.dt)?;
            t.push(s as f64 * self.dt);
            w.push(vp.field_energy());
        }
        require(
            w.iter().all(|v| v.is_finite()),
            "field energy became non-finite",
        )?;
        let logw: Vec<f64> = w.iter().map(|v| v.max(1e-300).ln()).collect();
        csv(
            out,
            m,
            "field_energy.csv",
            &[
                ("time", &t),
                ("field_energy", &w),
                ("log_field_energy", &logw),
            ],
        )?;
        // Damping from the successive maxima of the field energy (|E|² decays at 2γ).
        let peaks: Vec<(f64, f64)> = (1..w.len() - 1)
            .filter(|&i| w[i] > w[i - 1] && w[i] >= w[i + 1])
            .map(|i| (t[i], w[i]))
            .collect();
        if peaks.len() >= 3 {
            let (pt, pw): (Vec<f64>, Vec<f64>) = peaks.iter().take(8).copied().unzip();
            if let Some(s) = log_slope(&pt, &pw, 0.0, f64::MAX) {
                metric(m, "measured_damping_rate", 0.5 * s);
            }
        }
        // Reference: γ ≈ -sqrt(π/8) k⁻³ exp(-1/(2k²) - 3/2) (weak-damping asymptote).
        metric(
            m,
            "asymptotic_damping_rate",
            -(PI / 8.0).sqrt() / k.powi(3) * (-0.5 / (k * k) - 1.5).exp(),
        );
        metric(m, "mass_drift_relative", (vp.mass() - m0) / m0);
        Ok(())
    }
}

// ---------------------------------------------------------------- Lennard-Jones MD

/// Periodic 1-D Lennard-Jones gas in reduced units (ε = σ = m = k_B = 1), velocity Verlet with
/// an optional velocity-rescaling thermostat for the first `thermostat_steps` steps.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LennardJonesMdCase {
    pub particles: usize,
    pub box_length: f64,
    pub initial_temperature: f64,
    #[serde(default)]
    pub target_temperature: Option<f64>,
    #[serde(default)]
    pub thermostat_steps: usize,
    pub dt: f64,
    pub steps: usize,
    #[serde(default = "cutoff")]
    pub cutoff: f64,
}

fn cutoff() -> f64 {
    2.5
}

impl LennardJonesMdCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.particles >= 2
                && self.particles <= 2000
                && self.steps >= 1
                && self.steps <= 20_000,
            "use 2..=2000 particles and 1..=20000 steps",
        )?;
        let n = self.particles;
        let l = self.box_length;
        require(
            l > self.particles as f64 * 0.9,
            "box too dense for a 1-D Lennard-Jones chain (spacing below 0.9 σ)",
        )?;
        let ff = LennardJones1DForceField::new(LennardJonesParameters::default())?;
        let mut nl = PeriodicNeighborList1D::new(NeighborList1DConfig {
            periodic_length: l,
            cutoff: self.cutoff,
            skin: 0.3,
        })?;
        let mut x: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) * l / n as f64).collect();
        // Deterministic velocities with zero momentum, scaled to the initial temperature.
        let mut v: Vec<f64> = (0..n)
            .map(|i| (i * 7919 % 1000) as f64 / 1000.0 - 0.5)
            .collect();
        let mean = v.iter().sum::<f64>() / n as f64;
        v.iter_mut().for_each(|u| *u -= mean);
        let t0 = v.iter().map(|u| u * u).sum::<f64>() / n as f64;
        if t0 > 0.0 {
            let s = (self.initial_temperature / t0).sqrt();
            v.iter_mut().for_each(|u| *u *= s);
        }
        let mut f = vec![0.0; n];
        nl.rebuild(&x)?;
        let mut pe = ff.compute(&x, &nl, &mut f)?.potential_energy;
        let (mut t, mut kin, mut pot, mut temp) = (vec![], vec![], vec![], vec![]);
        let record = |t: &mut Vec<f64>,
                      kin: &mut Vec<f64>,
                      pot: &mut Vec<f64>,
                      temp: &mut Vec<f64>,
                      time: f64,
                      v: &[f64],
                      pe: f64| {
            let k: f64 = v.iter().map(|u| 0.5 * u * u).sum();
            t.push(time);
            kin.push(k);
            pot.push(pe);
            temp.push(2.0 * k / n as f64);
        };
        record(&mut t, &mut kin, &mut pot, &mut temp, 0.0, &v, pe);
        let dt = self.dt;
        for s in 1..=self.steps {
            for i in 0..n {
                v[i] += 0.5 * dt * f[i];
                x[i] = (x[i] + dt * v[i]).rem_euclid(l);
            }
            nl.rebuild(&x)?;
            pe = ff.compute(&x, &nl, &mut f)?.potential_energy;
            for i in 0..n {
                v[i] += 0.5 * dt * f[i];
            }
            if let Some(target) = self.target_temperature {
                if s <= self.thermostat_steps {
                    velocity_rescaling_thermostat(&mut v, 1.0, target, dt, 10.0 * dt)?;
                }
            }
            record(&mut t, &mut kin, &mut pot, &mut temp, s as f64 * dt, &v, pe);
        }
        require(
            kin.iter().chain(&pot).all(|e| e.is_finite()),
            "MD energies became non-finite; reduce dt",
        )?;
        let total: Vec<f64> = kin.iter().zip(&pot).map(|(a, b)| a + b).collect();
        csv(
            out,
            m,
            "energies.csv",
            &[
                ("time", &t),
                ("kinetic", &kin),
                ("potential", &pot),
                ("total", &total),
                ("temperature", &temp),
            ],
        )?;
        csv(out, m, "final_state.csv", &[("x", &x), ("v", &v)])?;
        let start = self.thermostat_steps.min(self.steps);
        let nve = &total[start..];
        let e_ref = nve[0];
        metric(
            m,
            "nve_energy_drift_relative",
            nve.iter().map(|e| (e - e_ref).abs()).fold(0.0, f64::max) / e_ref.abs().max(1e-12),
        );
        metric(
            m,
            "mean_temperature",
            temp[start..].iter().sum::<f64>() / (temp.len() - start) as f64,
        );
        metric(m, "final_potential_energy_per_particle", pe / n as f64);
        Ok(())
    }
}

// ---------------------------------------------------------------- photonic waveguide

/// Scalar finite-difference modes of a rectangular-core dielectric waveguide.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhotonicWaveguideCase {
    pub nx: usize,
    pub ny: usize,
    pub width_m: f64,
    pub height_m: f64,
    /// `[x_min, x_max, y_min, y_max]` of the core.
    pub core_m: [f64; 4],
    pub core_permittivity: f64,
    pub cladding_permittivity: f64,
    pub wavelength_m: f64,
    pub modes: usize,
}

impl PhotonicWaveguideCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.nx * self.ny <= 256 && self.modes >= 1 && self.modes <= 8,
            "use at most 256 grid cells (dense eigen solve) and 1..=8 modes",
        )?;
        let (dx, dy) = (
            self.width_m / self.nx as f64,
            self.height_m / self.ny as f64,
        );
        let c = self.core_m;
        let eps: Vec<f64> = (0..self.nx * self.ny)
            .map(|i| {
                let (x, y) = (
                    ((i % self.nx) as f64 + 0.5) * dx,
                    ((i / self.nx) as f64 + 0.5) * dy,
                );
                if x >= c[0] && x <= c[1] && y >= c[2] && y <= c[3] {
                    self.core_permittivity
                } else {
                    self.cladding_permittivity
                }
            })
            .collect();
        let grid = PhotonicGrid2D {
            nx: self.nx,
            ny: self.ny,
            width_m: self.width_m,
            height_m: self.height_m,
            relative_permittivity: eps.clone(),
            dispersive_permittivity: None,
        };
        require(self.wavelength_m > 0.0, "wavelength must be positive")?;
        let modes = solve_photonic_waveguide_modes(
            &grid,
            &PhotonicWaveguideConfig {
                frequency_hz: 299_792_458.0 / self.wavelength_m,
                mode_count: self.modes,
                ..Default::default()
            },
        )?;
        let idx: Vec<f64> = (1..=modes.len()).map(|i| i as f64).collect();
        let neff: Vec<f64> = modes.iter().map(|md| md.effective_index).collect();
        let beta: Vec<f64> = modes
            .iter()
            .map(|md| md.propagation_constant_rad_per_m)
            .collect();
        let guided: Vec<f64> = modes
            .iter()
            .map(|md| if md.guided { 1.0 } else { 0.0 })
            .collect();
        let res: Vec<f64> = modes.iter().map(|md| md.relative_residual).collect();
        csv(
            out,
            m,
            "modes.csv",
            &[
                ("mode", &idx),
                ("effective_index", &neff),
                ("propagation_constant_rad_per_m", &beta),
                ("guided", &guided),
                ("relative_residual", &res),
            ],
        )?;
        let mut cells = vec![DataArray::scalar("relative_permittivity", eps)];
        for (i, md) in modes.iter().enumerate() {
            require(
                md.field.len() == self.nx * self.ny,
                "unexpected mode field size",
            )?;
            let peak = md
                .field
                .iter()
                .map(|z| z.norm())
                .fold(0.0, f64::max)
                .max(f64::MIN_POSITIVE);
            cells.push(DataArray::scalar(
                format!("mode_{}_amplitude", i + 1),
                md.field.iter().map(|z| z.norm() / peak).collect(),
            ));
        }
        write_field(
            out,
            m,
            "modes.vtu",
            &VtuMesh::rectangle_cells(self.nx, self.ny, dx, dy)?,
            vec![],
            cells,
        )?;
        metric(
            m,
            "fundamental_effective_index",
            neff.first().copied().unwrap_or(0.0),
        );
        metric(m, "guided_modes", guided.iter().sum());
        m.linear_solver_converged = Some(modes.iter().all(|md| md.converged));
        Ok(())
    }
}

// ---------------------------------------------------------------- THMC column

/// Coupled pressure, temperature, species and reaction column with one hot, pressurised
/// cell region, evolving by the periodic THMC slice reference.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThmcColumnCase {
    pub cells: usize,
    #[serde(default = "one")]
    pub length: f64,
    pub pressure_diffusivity: f64,
    pub thermal_diffusivity: f64,
    pub species_diffusivity: f64,
    #[serde(default)]
    pub reaction_rate: f64,
    #[serde(default)]
    pub thermo_pressure_coupling: f64,
    #[serde(default)]
    pub reaction_heat: f64,
    pub initial: [f64; 3],
    /// Cells `[first, last]` set to `source_state` at t = 0.
    pub source_cells: [usize; 2],
    pub source_state: [f64; 3],
    pub dt: f64,
    pub steps: usize,
}

impl ThmcColumnCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.cells >= 4 && self.cells <= 4096 && self.steps >= 1 && self.steps <= 20_000,
            "use 4..=4096 cells and 1..=20000 steps",
        )?;
        require(
            self.source_cells[0] <= self.source_cells[1] && self.source_cells[1] < self.cells,
            "source cells must lie in the column",
        )?;
        let mut th = Thmc1D::new(Thmc1DConfig {
            cells: self.cells,
            length: self.length,
            pressure_diffusivity: self.pressure_diffusivity,
            thermal_diffusivity: self.thermal_diffusivity,
            species_diffusivity: self.species_diffusivity,
            reaction_rate: self.reaction_rate,
            thermo_pressure_coupling: self.thermo_pressure_coupling,
            reaction_heat: self.reaction_heat,
        })?;
        th.initialize(self.initial[0], self.initial[1], self.initial[2])?;
        let s = self.source_state;
        for c in self.source_cells[0]..=self.source_cells[1] {
            th.set_cell(c, s[0], s[1], s[2])
                .map_err(|e| SolverError::InvalidInput(format!("source cell {c}: {e}")))?;
        }
        let d0 = th.diagnostics();
        let (mut t, mut pm, mut tm, mut inv) = (
            vec![0.0],
            vec![d0.pressure_mean],
            vec![d0.temperature_mean],
            vec![d0.species_inventory],
        );
        for k in 1..=self.steps {
            th.step(self.dt)?;
            let d = th.diagnostics();
            t.push(k as f64 * self.dt);
            pm.push(d.pressure_mean);
            tm.push(d.temperature_mean);
            inv.push(d.species_inventory);
        }
        csv(
            out,
            m,
            "history.csv",
            &[
                ("time", &t),
                ("pressure_mean", &pm),
                ("temperature_mean", &tm),
                ("species_inventory", &inv),
            ],
        )?;
        let x: Vec<f64> = (0..self.cells)
            .map(|i| (i as f64 + 0.5) * self.length / self.cells as f64)
            .collect();
        csv(
            out,
            m,
            "profiles.csv",
            &[
                ("x", &x),
                ("pressure", th.pressure()),
                ("temperature", th.temperature()),
                ("species", th.species()),
                ("displacement", th.displacement()),
            ],
        )?;
        let d = th.diagnostics();
        metric(
            m,
            "species_inventory_change",
            d.species_inventory - d0.species_inventory,
        );
        metric(m, "reaction_heat", d.reaction_heat);
        metric(
            m,
            "max_temperature",
            th.temperature().iter().copied().fold(f64::MIN, f64::max),
        );
        Ok(())
    }
}

// ---------------------------------------------------------------- multiphase

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MultiphaseModel {
    /// Two-fluid (Euler–Euler) transport with drag and per-phase body accelerations.
    #[default]
    EulerEuler,
    /// Compressible VOF slug carried by a closed recirculating face flux.
    CompressibleVof,
}

fn half_one() -> f64 {
    0.5
}

/// Closed box on hexahedral FVM cells. `euler_euler`: a dispersed phase initially below
/// `initial_height_fraction` of the box, driven by the phase body accelerations (buoyant
/// bubbles rise). `compressible_vof`: liquid initially left of `initial_height_fraction` of
/// the width, stirred by a divergence-free vortex of peak face speed `stir_speed_m_s`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MultiphaseBoxCase {
    #[serde(default)]
    pub model: MultiphaseModel,
    pub size_m: [f64; 3],
    pub cells: [usize; 3],
    pub dt_s: f64,
    pub steps: usize,
    #[serde(default = "five")]
    pub frames: usize,
    #[serde(default = "half_one")]
    pub initial_height_fraction: f64,
    #[serde(default)]
    pub dispersed_fraction: f64,
    #[serde(default)]
    pub drag_coefficient: Option<f64>,
    #[serde(default)]
    pub primary_acceleration_m_s2: [f64; 3],
    #[serde(default)]
    pub dispersed_acceleration_m_s2: [f64; 3],
    #[serde(default)]
    pub stir_speed_m_s: f64,
}

fn five() -> usize {
    5
}

impl MultiphaseBoxCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        use crate::fvm::{
            make_sheared_cartesian_hexa_mesh, CompressibleVofConfig, CompressibleVofTransport,
            EulerEulerConfig, EulerEulerTransport,
        };
        let [nx, ny, nz] = self.cells;
        require(
            nx * ny * nz >= 1
                && nx * ny * nz <= 32_768
                && self.steps >= 1
                && self.steps <= 20_000
                && (1..=50).contains(&self.frames),
            "use at most 32768 cells, 1..=20000 steps and 1..=50 frames",
        )?;
        let [lx, ly, lz] = self.size_m;
        let mesh = make_sheared_cartesian_hexa_mesh(nx, ny, nz, lx, ly, lz, 0.0)?;
        let export = VtuMesh::box_cells([nx, ny, nz], Vec3::new(lx, ly, lz), 0.0)?;
        let centres: Vec<Vec3> = mesh.cells().iter().map(|c| c.center).collect();
        let volumes: Vec<f64> = mesh.cells().iter().map(|c| c.volume).collect();
        let interval = (self.steps / self.frames).max(1);
        let h = self.initial_height_fraction;
        let mut frames = Vec::new();
        let (mut t, mut inv, mut centroid) = (vec![], vec![], vec![]);
        let centroid_of = |alpha: &[f64]| -> f64 {
            let w: f64 = alpha.iter().zip(&volumes).map(|(a, v)| a * v).sum();
            alpha
                .iter()
                .zip(&volumes)
                .zip(&centres)
                .map(|((a, v), c)| a * v * c.y)
                .sum::<f64>()
                / w.max(f64::MIN_POSITIVE)
        };
        match self.model {
            MultiphaseModel::EulerEuler => {
                require(
                    self.dispersed_fraction > 0.0 && self.dispersed_fraction < 1.0,
                    "dispersed_fraction must lie in (0, 1)",
                )?;
                let a = |v: [f64; 3]| Vec3::new(v[0], v[1], v[2]);
                let mut cfg = EulerEulerConfig {
                    dt: self.dt_s,
                    primary_body_acceleration: a(self.primary_acceleration_m_s2),
                    dispersed_body_acceleration: a(self.dispersed_acceleration_m_s2),
                    ..Default::default()
                };
                if let Some(c) = self.drag_coefficient {
                    cfg.drag_coefficient = c;
                }
                let mut ee = EulerEulerTransport::new(mesh, cfg)?;
                let alpha = self.dispersed_fraction;
                ee.initialize_profile(
                    move |p| if p.y < h * ly { alpha } else { 0.0 },
                    Vec3::default(),
                    Vec3::default(),
                )?;
                let v0 = ee.volume_fraction_integral();
                for step in 0..=self.steps {
                    if step > 0 {
                        ee.step();
                    }
                    if step % interval == 0 || step == self.steps {
                        let alpha = ee.dispersed_fraction().to_vec();
                        require(
                            alpha.iter().all(|x| x.is_finite()),
                            "phase fraction became non-finite; reduce dt_s",
                        )?;
                        t.push(ee.time());
                        inv.push(ee.volume_fraction_integral());
                        centroid.push(centroid_of(&alpha));
                        if frames.len() <= self.frames {
                            let file = format!("multiphase_{:04}.vtu", frames.len());
                            write_field(
                                out,
                                m,
                                &file,
                                &export,
                                vec![],
                                vec![
                                    DataArray::scalar("dispersed_fraction", alpha),
                                    DataArray::vector(
                                        "dispersed_velocity_m_s",
                                        ee.dispersed_velocity().to_vec(),
                                    ),
                                    DataArray::vector(
                                        "primary_velocity_m_s",
                                        ee.primary_velocity().to_vec(),
                                    ),
                                ],
                            )?;
                            frames.push((ee.time(), file));
                        }
                    }
                }
                metric(
                    m,
                    "volume_fraction_drift_relative",
                    (ee.volume_fraction_integral() - v0) / v0,
                );
                metric(m, "max_dispersed_fraction", ee.maximum_fraction());
            }
            MultiphaseModel::CompressibleVof => {
                // Discretely divergence-free vortex from psi = U·L/π · sin(πx/Lx) sin(πy/Ly).
                let s = self.stir_speed_m_s;
                let psi = |x: f64, y: f64| {
                    s * lx / std::f64::consts::PI
                        * (std::f64::consts::PI * x / lx).sin()
                        * (std::f64::consts::PI * y / ly).sin()
                };
                let (hx, hy, dz) = (0.5 * lx / nx as f64, 0.5 * ly / ny as f64, lz / nz as f64);
                let flux: Vec<f64> = mesh
                    .faces()
                    .iter()
                    .map(|f| {
                        let c = f.center;
                        if f.area.x.abs() > 1e-12 * f.area.norm() {
                            f.area.x.signum() * (psi(c.x, c.y + hy) - psi(c.x, c.y - hy)) * dz
                        } else if f.area.y.abs() > 1e-12 * f.area.norm() {
                            -f.area.y.signum() * (psi(c.x + hx, c.y) - psi(c.x - hx, c.y)) * dz
                        } else {
                            0.0
                        }
                    })
                    .collect();
                let mut vof = CompressibleVofTransport::new(
                    mesh,
                    CompressibleVofConfig {
                        dt: self.dt_s,
                        ..Default::default()
                    },
                )?;
                vof.initialize_profile(move |p| if p.x < h * lx { 1.0 } else { 0.0 }, 0.0)?;
                vof.set_face_flux(flux)?;
                let (m0, l0) = (vof.total_mass(), vof.liquid_volume());
                for step in 0..=self.steps {
                    if step > 0 {
                        vof.step()?;
                    }
                    if step % interval == 0 || step == self.steps {
                        let alpha = vof.liquid_fraction().to_vec();
                        t.push(vof.time());
                        inv.push(vof.liquid_volume());
                        centroid.push(centroid_of(&alpha));
                        if frames.len() <= self.frames {
                            let file = format!("multiphase_{:04}.vtu", frames.len());
                            write_field(
                                out,
                                m,
                                &file,
                                &export,
                                vec![],
                                vec![
                                    DataArray::scalar("liquid_fraction", alpha),
                                    DataArray::scalar("density_kg_m3", vof.density().to_vec()),
                                    DataArray::scalar("pressure_Pa", vof.pressure().to_vec()),
                                ],
                            )?;
                            frames.push((vof.time(), file));
                        }
                    }
                }
                metric(m, "mass_drift_relative", (vof.total_mass() - m0) / m0);
                metric(
                    m,
                    "liquid_volume_drift_relative",
                    (vof.liquid_volume() - l0) / l0.max(f64::MIN_POSITIVE),
                );
            }
        }
        crate::io::write_pvd(out.join("multiphase.pvd"), &frames)?;
        m.files.push("multiphase.pvd".into());
        csv(
            out,
            m,
            "phase_history.csv",
            &[
                ("time_s", &t),
                ("phase_volume_m3", &inv),
                ("phase_centroid_y_m", &centroid),
            ],
        )?;
        metric(
            m,
            "initial_centroid_y_m",
            centroid.first().copied().unwrap_or(0.0),
        );
        metric(
            m,
            "final_centroid_y_m",
            centroid.last().copied().unwrap_or(0.0),
        );
        Ok(())
    }
}

// ---------------------------------------------------------------- spray

/// Lagrangian spray injected into a uniform carrier stream in a box: drag, heating and d²-law
/// evaporation, with vapour mass sources accumulated per cell. At most the spray
/// reference's parcel limit (500 by default) is tracked.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SprayInjectionCase {
    pub size_m: [f64; 3],
    pub cells: [usize; 3],
    pub carrier_velocity_m_s: [f64; 3],
    pub carrier_temperature_k: f64,
    pub nozzle_m: [f64; 3],
    pub injection_velocity_m_s: [f64; 3],
    pub droplet_diameter_m: f64,
    pub mass_flow_rate_kg_s: f64,
    pub injection_duration_s: f64,
    #[serde(default)]
    pub evaporation_constant_m2_s: Option<f64>,
    pub dt_s: f64,
    pub steps: usize,
    #[serde(default = "five")]
    pub frames: usize,
}

impl SprayInjectionCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        use crate::fvm::{
            make_sheared_cartesian_hexa_mesh, SprayConfig, SprayInjection,
            SprayInjectionEvaporation,
        };
        let [nx, ny, nz] = self.cells;
        require(
            nx * ny * nz >= 1
                && nx * ny * nz <= 32_768
                && self.steps >= 1
                && self.steps <= 20_000
                && (1..=50).contains(&self.frames),
            "use at most 32768 cells, 1..=20000 steps and 1..=50 frames",
        )?;
        let [lx, ly, lz] = self.size_m;
        let mesh = make_sheared_cartesian_hexa_mesh(nx, ny, nz, lx, ly, lz, 0.0)?;
        let n = mesh.cell_count();
        let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
        let mut cfg = SprayConfig {
            dt: self.dt_s,
            ..Default::default()
        };
        if let Some(k) = self.evaporation_constant_m2_s {
            cfg.evaporation_constant = k;
        }
        let mut spray = SprayInjectionEvaporation::new(mesh, cfg)?;
        spray.set_injection(SprayInjection {
            position: v(self.nozzle_m),
            velocity: v(self.injection_velocity_m_s),
            diameter: self.droplet_diameter_m,
            mass_flow_rate: self.mass_flow_rate_kg_s,
            duration: self.injection_duration_s,
            ..Default::default()
        })?;
        let u = vec![v(self.carrier_velocity_m_s); n];
        let interval = (self.steps / self.frames).max(1);
        let mut frames = Vec::new();
        let (mut t, mut liquid, mut injected, mut evaporated, mut parcels) =
            (vec![], vec![], vec![], vec![], vec![]);
        let mut vapor = vec![0.0; n];
        for step in 1..=self.steps {
            spray.step_uniform(&u, self.carrier_temperature_k)?;
            for (acc, s) in vapor.iter_mut().zip(spray.vapor_mass_source()) {
                *acc += s * self.dt_s;
            }
            if step % interval == 0 || step == self.steps {
                t.push(spray.time());
                liquid.push(spray.total_liquid_mass());
                injected.push(spray.injected_mass());
                evaporated.push(spray.evaporated_mass());
                parcels.push(spray.parcels().len() as f64);
                if frames.len() < self.frames && !spray.parcels().is_empty() {
                    let p = spray.parcels();
                    let points: Vec<Vec3> = p.iter().map(|q| q.position).collect();
                    let cells = (0..points.len())
                        .map(|i| crate::io::Cell {
                            cell_type: crate::io::CellType::Vertex,
                            nodes: vec![i],
                        })
                        .collect();
                    let file = format!("parcels_{:04}.vtu", frames.len());
                    write_field(
                        out,
                        m,
                        &file,
                        &VtuMesh { points, cells },
                        vec![
                            DataArray::scalar("diameter_m", p.iter().map(|q| q.diameter).collect()),
                            DataArray::scalar(
                                "temperature_K",
                                p.iter().map(|q| q.temperature).collect(),
                            ),
                            DataArray::vector(
                                "velocity_m_s",
                                p.iter().map(|q| q.velocity).collect(),
                            ),
                        ],
                        vec![],
                    )?;
                    frames.push((spray.time(), file));
                }
            }
        }
        crate::io::write_pvd(out.join("parcels.pvd"), &frames)?;
        m.files.push("parcels.pvd".into());
        write_field(
            out,
            m,
            "vapor_source.vtu",
            &VtuMesh::box_cells([nx, ny, nz], Vec3::new(lx, ly, lz), 0.0)?,
            vec![],
            vec![DataArray::scalar("evaporated_mass_kg", vapor)],
        )?;
        csv(
            out,
            m,
            "spray_history.csv",
            &[
                ("time_s", &t),
                ("liquid_mass_kg", &liquid),
                ("injected_mass_kg", &injected),
                ("evaporated_mass_kg", &evaporated),
                ("parcels", &parcels),
            ],
        )?;
        let inj = spray.injected_mass();
        // Injected = liquid + evaporated (+ any mass removed with parcels); report the residual.
        let residual = inj - spray.total_liquid_mass() - spray.evaporated_mass();
        metric(m, "injected_mass_kg", inj);
        metric(
            m,
            "evaporated_fraction",
            spray.evaporated_mass() / inj.max(f64::MIN_POSITIVE),
        );
        metric(
            m,
            "liquid_in_box_fraction",
            spray.total_liquid_mass() / inj.max(f64::MIN_POSITIVE),
        );
        metric(
            m,
            "mass_balance_residual_fraction",
            residual / inj.max(f64::MIN_POSITIVE),
        );
        metric(m, "final_parcels", spray.parcels().len() as f64);
        Ok(())
    }
}
