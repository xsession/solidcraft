//! Dust carried by the settled flow of a `collocated_flow3d` run (`dust_deposition_v1`).
//!
//! After the flow (and heat) has been solved, particles of each size class are released and
//! followed through the final velocity field until they reach a wall, leave through an open
//! boundary or run out of time. Each particle obeys Stokes drag with the Cunningham slip
//! correction, gravity less buoyancy, thermophoresis (Talbot, Cheng, Schefer and Willis 1980)
//! and Brownian motion:
//!
//! dv/dt = (u + v_th − v)/τ + g (1 − ρ_f/ρ_p),   τ = ρ_p d² C_c / (18 μ),
//!
//! integrated exactly over each sub-step with the local air velocity held fixed. The air velocity
//! is interpolated trilinearly between cell centres, with a no-slip mirror behind walls, so it
//! falls to zero at a wall instead of carrying a particle into it at full speed.
//!
//! Particles enter with the inflow (`release: inflow`, the default): on the open boundary faces
//! where air flows in, in proportion to each face's inflow, so the dust entering is the
//! concentration times the inflowing volume. Or they start spread through the air
//! (`release: volume`), in proportion to each fluid cell's volume. A deposit counts for the
//! solid, plate or exterior wall it reaches, and carries the temperature of the surface there.
use super::{CollocatedFlowCase, RunManifest};
use crate::{
    core::{Mt19937_64, StdNormal, Vec3},
    fvm::{PolyMesh, PressureBoundaryType, VelocityBoundaryType},
    require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, f64::consts::PI, fmt::Write as _, path::Path};

/// Boltzmann's constant, J/K.
const BOLTZMANN: f64 = 1.380_649e-23;
/// Talbot's thermophoresis coefficients (thermal slip, temperature jump, momentum exchange).
const TALBOT_CS: f64 = 1.17;
const TALBOT_CT: f64 = 2.18;
const TALBOT_CM: f64 = 1.14;
/// Most particles of one run, and most sub-steps of one particle.
const MAX_PARTICLES: usize = 100_000;
const MAX_STEPS: usize = 2_000_000;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct DustClass {
    pub diameter_m: f64,
    /// Mass of this size per volume of air entering (inflow) or of the air (volume), kg/m³.
    #[serde(default)]
    pub concentration_kg_per_m3: f64,
}

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum DustRelease {
    #[default]
    Inflow,
    Volume,
    /// Where air flows into `release_box` across its sides: the dust entering an enclosure
    /// through its openings, when the box lies just inside the enclosure's walls.
    Box,
}

/// A box in the domain, m; its sides go to the nearest grid nodes.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct DustBox {
    pub min_m: [f64; 3],
    pub max_m: [f64; 3],
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Dust {
    pub classes: Vec<DustClass>,
    pub gravity_m_per_s2: [f64; 3],
    #[serde(default = "default_density")]
    pub density_kg_per_m3: f64,
    /// The particles' own conductivity, for thermophoresis, W/(m·K).
    #[serde(default = "default_conductivity")]
    pub conductivity_w_per_m_k: f64,
    #[serde(default = "default_particles")]
    pub particles_per_class: usize,
    #[serde(default)]
    pub release: DustRelease,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub release_box: Option<DustBox>,
    #[serde(default = "default_seed")]
    pub seed: u64,
    /// A particle still in the air after this long counts as airborne.
    #[serde(default = "default_max_time")]
    pub max_time_s: f64,
    #[serde(default = "yes")]
    pub brownian: bool,
    #[serde(default = "yes")]
    pub thermophoresis: bool,
    /// The air's temperature where the run has no heat transfer, K.
    #[serde(default = "default_temperature")]
    pub temperature_k: f64,
    /// The air's molecular mean free path, m.
    #[serde(default = "default_mean_free_path")]
    pub mean_free_path_m: f64,
}

fn default_density() -> f64 {
    2650.0
}
fn default_conductivity() -> f64 {
    1.0
}
fn default_particles() -> usize {
    1000
}
fn default_seed() -> u64 {
    1
}
fn default_max_time() -> f64 {
    3600.0
}
fn yes() -> bool {
    true
}
fn default_temperature() -> f64 {
    293.15
}
fn default_mean_free_path() -> f64 {
    6.8e-8
}

impl Dust {
    pub fn validate(&self) -> Result<()> {
        require(
            !self.classes.is_empty() && self.classes.len() <= 16,
            "dust needs 1..=16 size classes",
        )?;
        for c in &self.classes {
            require(
                c.diameter_m.is_finite() && (1e-8..=1e-3).contains(&c.diameter_m),
                "dust diameters must be from 10 nm to 1 mm",
            )?;
            require(
                c.concentration_kg_per_m3.is_finite() && c.concentration_kg_per_m3 >= 0.0,
                "dust concentrations must be finite and non-negative",
            )?;
        }
        require(
            self.gravity_m_per_s2
                .iter()
                .all(|g| g.is_finite() && g.abs() <= 1e3),
            "dust gravity must be finite",
        )?;
        require(
            self.density_kg_per_m3.is_finite()
                && self.density_kg_per_m3 > 0.0
                && self.density_kg_per_m3 <= 1e5,
            "dust density must be positive",
        )?;
        require(
            self.conductivity_w_per_m_k.is_finite() && self.conductivity_w_per_m_k > 0.0,
            "dust conductivity must be positive",
        )?;
        require(
            self.particles_per_class >= 1
                && self.particles_per_class.saturating_mul(self.classes.len()) <= MAX_PARTICLES,
            format!("use at least 1 particle per class and at most {MAX_PARTICLES} in all"),
        )?;
        require(
            self.max_time_s.is_finite() && self.max_time_s > 0.0 && self.max_time_s <= 1e7,
            "dust max_time_s must be positive and at most 1e7 s",
        )?;
        match (&self.release, &self.release_box) {
            (DustRelease::Box, Some(b)) => require(
                b.min_m.iter().chain(&b.max_m).all(|v| v.is_finite())
                    && b.min_m.iter().zip(&b.max_m).all(|(l, h)| l < h),
                "the dust release box needs finite corners, min below max",
            )?,
            (DustRelease::Box, None) => {
                return Err(SolverError::InvalidInput(
                    "release box needs release_box".into(),
                ))
            }
            (_, Some(_)) => {
                return Err(SolverError::InvalidInput(
                    "release_box goes with release box".into(),
                ))
            }
            _ => {}
        }
        require(
            self.temperature_k.is_finite() && self.temperature_k > 0.0,
            "dust temperature must be positive kelvin",
        )?;
        require(
            self.mean_free_path_m.is_finite()
                && self.mean_free_path_m > 0.0
                && self.mean_free_path_m < 1e-3,
            "dust mean free path must be positive",
        )?;
        Ok(())
    }
}

/// What lies beyond one side of a cell.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Side {
    Cell(usize),
    /// A wall patch, and whether the air slips along it.
    Wall(usize, bool),
    Open(usize),
}

/// The background grid, which flow cell each of its cells is, and each flow cell's sides.
struct Grid {
    nodes: [Vec<f64>; 3],
    centres: [Vec<f64>; 3],
    /// Per background cell and axis, the step a particle may take there in one go.
    reach: Vec<[f64; 3]>,
    flow_of: Vec<Option<usize>>,
    background_of: Vec<usize>,
    /// Sides per flow cell: −x, +x, −y, +y, −z, +z.
    sides: Vec<[Side; 6]>,
    targets: Vec<String>,
}

impl Grid {
    fn new(case: &CollocatedFlowCase, mesh: &PolyMesh) -> Result<Self> {
        let nodes = [case.mesh.nodes(0), case.mesh.nodes(1), case.mesh.nodes(2)];
        let n = [
            nodes[0].len().saturating_sub(1),
            nodes[1].len().saturating_sub(1),
            nodes[2].len().saturating_sub(1),
        ];
        let count = n[0] * n[1] * n[2];
        let centres = nodes.clone().map(|v| {
            v.windows(2)
                .map(|w| 0.5 * (w[0] + w[1]))
                .collect::<Vec<_>>()
        });
        let solid = if case.obstacles.is_empty() && case.bodies.is_empty() {
            vec![None; count]
        } else {
            super::obstacles::layout(&case.mesh, &case.obstacles, &case.bodies)?.cell_solid
        };
        let mut flow_of = vec![None; count];
        let mut background_of = Vec::with_capacity(mesh.cell_count());
        for (b, s) in solid.iter().enumerate() {
            if s.is_none() {
                flow_of[b] = Some(background_of.len());
                background_of.push(b);
            }
        }
        require(
            background_of.len() == mesh.cell_count(),
            "dust needs the staircase grid's cells (no cut cells)",
        )?;
        let width = |a: usize, i: usize| -> f64 {
            match (nodes[a].get(i), nodes[a].get(i + 1)) {
                (Some(l), Some(h)) => h - l,
                _ => f64::MAX,
            }
        };
        let mut reach = vec![[0.0; 3]; count];
        for (b, r) in reach.iter_mut().enumerate() {
            let ijk = unflatten(b, n);
            for a in 0..3 {
                let i = ijk[a];
                let narrow = width(a, i)
                    .min(if i > 0 { width(a, i - 1) } else { f64::MAX })
                    .min(width(a, i + 1));
                r[a] = 0.45 * narrow;
            }
        }
        // Body face conditions deposit on their body.
        let mut targets: Vec<String> = mesh.patches().iter().map(|p| p.name.clone()).collect();
        for c in &case.face_conditions {
            for t in targets.iter_mut() {
                if *t == c.name {
                    *t = c.body.clone();
                }
            }
        }
        let open: Vec<bool> = mesh
            .patches()
            .iter()
            .map(|p| {
                case.pressure_boundaries
                    .iter()
                    .any(|b| b.patch == p.name && b.kind == PressureBoundaryType::FixedValue)
                    || case.velocity_boundaries.iter().any(|b| {
                        b.patch == p.name
                            && match b.kind {
                                VelocityBoundaryType::FixedValue => b.value.norm() > 0.0,
                                VelocityBoundaryType::ZeroGradient => true,
                                VelocityBoundaryType::Slip => false,
                            }
                    })
                    || case.fan.as_ref().is_some_and(|f| f.patch == p.name)
            })
            .collect();
        let slip: Vec<bool> = mesh
            .patches()
            .iter()
            .map(|p| {
                case.velocity_boundaries
                    .iter()
                    .any(|b| b.patch == p.name && b.kind == VelocityBoundaryType::Slip)
            })
            .collect();
        let mut sides = vec![[Side::Wall(usize::MAX, false); 6]; mesh.cell_count()];
        for f in mesh.faces() {
            let axis = dominant_axis(f.area);
            let positive = f.area.component(axis) > 0.0;
            let side = 2 * axis + usize::from(positive);
            let opposite = 2 * axis + usize::from(!positive);
            match (f.neighbour, f.patch) {
                (Some(nb), _) => {
                    if let Some(s) = sides.get_mut(f.owner) {
                        s[side] = Side::Cell(nb);
                    }
                    if let Some(s) = sides.get_mut(nb) {
                        s[opposite] = Side::Cell(f.owner);
                    }
                }
                (None, Some(p)) => {
                    if let Some(s) = sides.get_mut(f.owner) {
                        s[side] = if open.get(p).copied().unwrap_or(false) {
                            Side::Open(p)
                        } else {
                            Side::Wall(p, slip.get(p).copied().unwrap_or(false))
                        };
                    }
                }
                (None, None) => {}
            }
        }
        Ok(Self {
            nodes,
            centres,
            reach,
            flow_of,
            background_of,
            sides,
            targets,
        })
    }

    fn dims(&self) -> [usize; 3] {
        [
            self.centres[0].len(),
            self.centres[1].len(),
            self.centres[2].len(),
        ]
    }

    fn width(&self, a: usize, i: usize) -> f64 {
        match (self.nodes[a].get(i), self.nodes[a].get(i + 1)) {
            (Some(l), Some(h)) => h - l,
            _ => 0.0,
        }
    }

    fn background_ijk(&self, cell: usize) -> [usize; 3] {
        unflatten(
            self.background_of.get(cell).copied().unwrap_or(0),
            self.dims(),
        )
    }

    fn side(&self, cell: usize, axis: usize, positive: bool) -> Side {
        self.sides
            .get(cell)
            .map_or(Side::Wall(usize::MAX, false), |s| {
                s[2 * axis + usize::from(positive)]
            })
    }

    fn target(&self, patch: usize) -> &str {
        self.targets.get(patch).map_or("wall", String::as_str)
    }

    /// The air velocity at `p` in flow cell `cell`: trilinear between cell centres, mirrored
    /// (no slip) behind walls and copied (zero gradient) through openings.
    fn velocity(&self, u: &[Vec3], cell: usize, p: Vec3) -> Vec3 {
        let ijk = self.background_ijk(cell);
        let mut dir = [true; 3];
        let mut weight = [0.0; 3];
        for a in 0..3 {
            let c = self.centres[a].get(ijk[a]).copied().unwrap_or(0.0);
            let x = p.component(a);
            dir[a] = x >= c;
            let distance = match self.side(cell, a, dir[a]) {
                Side::Cell(n) => {
                    let nb = self.background_ijk(n);
                    (self.centres[a].get(nb[a]).copied().unwrap_or(c) - c).abs()
                }
                _ => self.width(a, ijk[a]),
            };
            weight[a] = if distance > 0.0 {
                ((x - c).abs() / distance).clamp(0.0, 1.0)
            } else {
                0.0
            };
        }
        let own = u.get(cell).copied().unwrap_or(Vec3::ZERO);
        let mut total = Vec3::ZERO;
        for corner in 0..8usize {
            let mut value = own;
            let mut at = Some(cell);
            let mut w = 1.0;
            for a in 0..3 {
                let step = corner >> a & 1 == 1;
                w *= if step { weight[a] } else { 1.0 - weight[a] };
                if !step {
                    continue;
                }
                if let Some(c) = at {
                    match self.side(c, a, dir[a]) {
                        Side::Cell(n) => {
                            at = Some(n);
                            value = u.get(n).copied().unwrap_or(Vec3::ZERO);
                        }
                        Side::Wall(_, false) => {
                            at = None;
                            value = -value;
                        }
                        Side::Wall(_, true) => {
                            // Slip: only the normal component turns over.
                            at = None;
                            value = set(value, a, -value.component(a));
                        }
                        Side::Open(_) => at = None,
                    }
                }
            }
            if w > 0.0 {
                total = total + value * w;
            }
        }
        total
    }

    /// Temperature and its gradient in a flow cell, from the background grid's temperatures
    /// (solid cells included, so a cold wall pulls its neighbours' gradient towards it).
    fn temperature(&self, t: &[f64], cell: usize) -> (f64, Vec3) {
        let ijk = self.background_ijk(cell);
        let n = self.dims();
        let at = |q: [usize; 3]| t.get(flatten(q, n)).copied();
        let own = at(ijk).unwrap_or(0.0);
        let mut g = [0.0; 3];
        for a in 0..3 {
            let mut lo = ijk;
            let mut hi = ijk;
            if ijk[a] > 0 {
                lo[a] -= 1;
            }
            if ijk[a] + 1 < n[a] {
                hi[a] += 1;
            }
            let (Some(tl), Some(th)) = (at(lo), at(hi)) else {
                continue;
            };
            let dx = self.centres[a].get(hi[a]).copied().unwrap_or(0.0)
                - self.centres[a].get(lo[a]).copied().unwrap_or(0.0);
            if dx > 0.0 {
                g[a] = (th - tl) / dx;
            }
        }
        (own, Vec3::from(g))
    }

    /// The flow cell holding `p`, if `p` is in a fluid cell of the grid.
    fn locate(&self, p: Vec3) -> Option<usize> {
        let mut ijk = [0; 3];
        for (a, q) in ijk.iter_mut().enumerate() {
            let nodes = &self.nodes[a];
            let x = p.component(a);
            let last = *nodes.last()?;
            if !(x >= *nodes.first()? && x <= last) {
                return None;
            }
            *q = nodes
                .partition_point(|v| *v <= x)
                .saturating_sub(1)
                .min(nodes.len().saturating_sub(2));
        }
        self.flow_of
            .get(flatten(ijk, self.dims()))
            .copied()
            .flatten()
    }
}

fn dominant_axis(v: Vec3) -> usize {
    let a = [v.x.abs(), v.y.abs(), v.z.abs()];
    if a[0] >= a[1] && a[0] >= a[2] {
        0
    } else if a[1] >= a[2] {
        1
    } else {
        2
    }
}

fn flatten(q: [usize; 3], n: [usize; 3]) -> usize {
    q[0] + n[0] * (q[1] + n[1] * q[2])
}

fn unflatten(b: usize, n: [usize; 3]) -> [usize; 3] {
    let nx = n[0].max(1);
    let ny = n[1].max(1);
    [b % nx, (b / nx) % ny, b / (nx * ny)]
}

/// One size class's particle physics.
struct Physics {
    tau: f64,
    gravity: Vec3,
    thermophoretic: f64,
    diffusivity_per_k: f64,
}

impl Physics {
    fn new(dust: &Dust, diameter: f64, rho_f: f64, nu: f64, k_gas: f64) -> Self {
        let mu = rho_f * nu;
        let kn = 2.0 * dust.mean_free_path_m / diameter;
        // Cunningham's slip correction (Allen and Raabe's coefficients).
        let slip = 1.0 + kn * (1.257 + 0.4 * (-1.1 / kn).exp());
        let tau = dust.density_kg_per_m3 * diameter * diameter * slip / (18.0 * mu);
        let ratio = k_gas / dust.conductivity_w_per_m_k;
        let thermophoretic = 2.0 * TALBOT_CS * slip * (ratio + TALBOT_CT * kn)
            / ((1.0 + 3.0 * TALBOT_CM * kn) * (1.0 + 2.0 * ratio + 2.0 * TALBOT_CT * kn));
        let g = Vec3::from(dust.gravity_m_per_s2);
        Self {
            tau,
            gravity: g * (1.0 - rho_f / dust.density_kg_per_m3),
            thermophoretic: thermophoretic * nu,
            diffusivity_per_k: if dust.brownian {
                BOLTZMANN * slip / (3.0 * PI * mu * diameter)
            } else {
                0.0
            },
        }
    }
}

enum Fate {
    Deposited {
        target: String,
        at: Vec3,
        temperature_k: f64,
    },
    Escaped,
    Airborne,
}

struct Tracked {
    fate: Fate,
    time_s: f64,
}

/// What the background grid's temperatures are (none without heat transfer).
pub(super) struct DustFields<'a> {
    pub velocity: &'a [Vec3],
    pub face_flux: &'a [f64],
    pub background_temperature_k: Option<Vec<f64>>,
    pub gas_conductivity_w_per_m_k: f64,
}

#[allow(clippy::too_many_arguments)]
fn track(
    grid: &Grid,
    dust: &Dust,
    physics: &Physics,
    fields: &DustFields<'_>,
    rng: &mut Mt19937_64,
    normal: &mut StdNormal,
    start: Vec3,
    cell: usize,
) -> Tracked {
    let mut p = start;
    let mut cell = cell;
    let mut t = 0.0;
    let temperature = |c: usize| match &fields.background_temperature_k {
        Some(bt) => grid.temperature(bt, c),
        None => (dust.temperature_k, Vec3::ZERO),
    };
    let air = grid.velocity(fields.velocity, cell, p);
    let mut v = air + physics.gravity * physics.tau;
    for _ in 0..MAX_STEPS {
        if t >= dust.max_time_s {
            return Tracked {
                fate: Fate::Airborne,
                time_s: t,
            };
        }
        let air = grid.velocity(fields.velocity, cell, p);
        let (temp, grad) = temperature(cell);
        let drift =
            if dust.thermophoresis && fields.background_temperature_k.is_some() && temp > 0.0 {
                grad * (-physics.thermophoretic / temp)
            } else {
                Vec3::ZERO
            };
        let terminal = air + drift + physics.gravity * physics.tau;
        let diffusivity = physics.diffusivity_per_k * temp.max(1.0);
        let reach = grid
            .reach
            .get(grid.background_of.get(cell).copied().unwrap_or(0))
            .copied()
            .unwrap_or([0.0; 3]);
        let mut dt = dust.max_time_s - t;
        for a in 0..3 {
            let speed = v.component(a).abs().max(terminal.component(a).abs());
            if speed > 0.0 {
                dt = dt.min(reach[a] / speed);
            }
            if diffusivity > 0.0 {
                // Keep a Brownian step (three standard deviations) within the reach too.
                dt = dt.min((reach[a] / 3.0).powi(2) / (2.0 * diffusivity));
            }
        }
        if !(dt > 0.0) || !dt.is_finite() {
            return Tracked {
                fate: Fate::Airborne,
                time_s: t,
            };
        }
        let decay = if physics.tau > 0.0 {
            (-dt / physics.tau).exp()
        } else {
            0.0
        };
        let mut step = terminal * dt + (v - terminal) * (physics.tau * (1.0 - decay));
        v = terminal + (v - terminal) * decay;
        if diffusivity > 0.0 {
            let s = (2.0 * diffusivity * dt).sqrt();
            let jitter = Vec3::new(normal.sample(rng), normal.sample(rng), normal.sample(rng));
            // Clip to the reach so the step stays within the neighbouring cells.
            let mut j = jitter * s;
            for a in 0..3 {
                let limit = reach[a];
                let c = j.component(a).clamp(-limit, limit);
                j = set(j, a, c);
            }
            step = step + j;
        }
        for a in 0..3 {
            let c = step.component(a).clamp(-reach[a], reach[a]);
            step = set(step, a, c);
        }
        t += dt;
        // Cross the faces the step passes, in the order it reaches them.
        let ijk = grid.background_ijk(cell);
        let mut crossings: Vec<(f64, usize, bool)> = Vec::with_capacity(3);
        for a in 0..3 {
            let d = step.component(a);
            if d == 0.0 {
                continue;
            }
            let face = if d > 0.0 {
                grid.nodes[a].get(ijk[a] + 1).copied().unwrap_or(f64::MAX)
            } else {
                grid.nodes[a].get(ijk[a]).copied().unwrap_or(f64::MIN)
            };
            let f = (face - p.component(a)) / d;
            if (0.0..=1.0).contains(&f) {
                crossings.push((f, a, d > 0.0));
            }
        }
        crossings.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut moved_to = cell;
        for (f, a, positive) in crossings {
            match grid.side(moved_to, a, positive) {
                Side::Cell(n) => moved_to = n,
                Side::Wall(patch, _) => {
                    let at = p + step * f;
                    let temperature_k = wall_temperature(grid, fields, dust, moved_to, a, positive);
                    // The time it reached the wall: the share of the step taken so far.
                    return Tracked {
                        fate: Fate::Deposited {
                            target: grid.target(patch).to_owned(),
                            at,
                            temperature_k,
                        },
                        time_s: t - dt * (1.0 - f),
                    };
                }
                Side::Open(_) => {
                    return Tracked {
                        fate: Fate::Escaped,
                        time_s: t - dt * (1.0 - f),
                    }
                }
            }
        }
        p = p + step;
        cell = moved_to;
        // Keep the position inside the cell it is counted in (rounding at a face).
        let q = grid.background_ijk(cell);
        for a in 0..3 {
            let lo = grid.nodes[a].get(q[a]).copied().unwrap_or(0.0);
            let hi = grid.nodes[a].get(q[a] + 1).copied().unwrap_or(0.0);
            let w = hi - lo;
            let c = p.component(a).clamp(lo + 1e-9 * w, hi - 1e-9 * w);
            p = set(p, a, c);
        }
    }
    Tracked {
        fate: Fate::Airborne,
        time_s: t,
    }
}

fn set(v: Vec3, axis: usize, value: f64) -> Vec3 {
    match axis {
        0 => Vec3::new(value, v.y, v.z),
        1 => Vec3::new(v.x, value, v.z),
        _ => Vec3::new(v.x, v.y, value),
    }
}

/// The temperature of the surface a particle reached: the solid cell beyond the face, or the
/// fluid cell's own where the wall is a plate or the domain's side.
fn wall_temperature(
    grid: &Grid,
    fields: &DustFields<'_>,
    dust: &Dust,
    cell: usize,
    axis: usize,
    positive: bool,
) -> f64 {
    let Some(bt) = &fields.background_temperature_k else {
        return dust.temperature_k;
    };
    let n = grid.dims();
    let ijk = grid.background_ijk(cell);
    let own = bt
        .get(flatten(ijk, n))
        .copied()
        .unwrap_or(dust.temperature_k);
    let mut q = ijk;
    if positive {
        if q[axis] + 1 >= n[axis] {
            return own;
        }
        q[axis] += 1;
    } else {
        if q[axis] == 0 {
            return own;
        }
        q[axis] -= 1;
    }
    let b = flatten(q, n);
    match grid.flow_of.get(b) {
        Some(None) => bt.get(b).copied().unwrap_or(own),
        _ => own,
    }
}

/// Where particles start: points on inflow faces, or in fluid cells, with each start's weight.
fn starts(
    grid: &Grid,
    mesh: &PolyMesh,
    dust: &Dust,
    fields: &DustFields<'_>,
) -> Result<(Vec<(usize, f64)>, f64)> {
    match dust.release {
        DustRelease::Inflow => {
            let mut faces = Vec::new();
            let mut total = 0.0;
            for (fi, f) in mesh.faces().iter().enumerate() {
                if f.neighbour.is_some() {
                    continue;
                }
                let Some(flux) = fields.face_flux.get(fi) else {
                    continue;
                };
                if *flux < 0.0 {
                    let axis = dominant_axis(f.area);
                    let positive = f.area.component(axis) > 0.0;
                    if matches!(grid.side(f.owner, axis, positive), Side::Open(_)) {
                        faces.push((fi, -flux));
                        total += -flux;
                    }
                }
            }
            Ok((faces, total))
        }
        DustRelease::Box => {
            let Some(b) = &dust.release_box else {
                return Ok((Vec::new(), 0.0));
            };
            let snap = |a: usize, x: f64| -> f64 {
                let nodes = &grid.nodes[a];
                nodes
                    .iter()
                    .copied()
                    .min_by(|p, q| (p - x).abs().total_cmp(&(q - x).abs()))
                    .unwrap_or(x)
            };
            let lo = [
                snap(0, b.min_m[0]),
                snap(1, b.min_m[1]),
                snap(2, b.min_m[2]),
            ];
            let hi = [
                snap(0, b.max_m[0]),
                snap(1, b.max_m[1]),
                snap(2, b.max_m[2]),
            ];
            let inside = |p: Vec3| (0..3).all(|a| p.component(a) > lo[a] && p.component(a) < hi[a]);
            let mut faces = Vec::new();
            let mut total = 0.0;
            for (fi, f) in mesh.faces().iter().enumerate() {
                let Some(nb) = f.neighbour else { continue };
                let (Some(o), Some(n)) = (mesh.cells().get(f.owner), mesh.cells().get(nb)) else {
                    continue;
                };
                let (oi, ni) = (inside(o.center), inside(n.center));
                if oi == ni {
                    continue;
                }
                let Some(flux) = fields.face_flux.get(fi).copied() else {
                    continue;
                };
                // The flux runs from owner to neighbour; into the box is towards the inside cell.
                let inward = if ni { flux } else { -flux };
                if inward > 0.0 {
                    faces.push((fi, inward));
                    total += inward;
                }
            }
            Ok((faces, total))
        }
        DustRelease::Volume => {
            let cells: Vec<(usize, f64)> = mesh
                .cells()
                .iter()
                .enumerate()
                .map(|(c, cell)| (c, cell.volume))
                .collect();
            let total = cells.iter().map(|c| c.1).sum();
            Ok((cells, total))
        }
    }
}

fn pick(cumulative: &[f64], x: f64) -> usize {
    cumulative
        .partition_point(|c| *c <= x)
        .min(cumulative.len().saturating_sub(1))
}

/// Follow the dust through the settled flow and write `dust_deposition.json` and
/// `dust_deposits.csv`.
pub(super) fn run(
    case: &CollocatedFlowCase,
    dust: &Dust,
    mesh: &PolyMesh,
    fields: &DustFields<'_>,
    out: &Path,
    m: &mut RunManifest,
) -> Result<()> {
    dust.validate()?;
    require(!case.cut_cells, "dust is not supported with cut cells yet")?;
    require(case.mesh.xy_shear == 0.0, "dust needs an unsheared grid")?;
    require(
        case.rotating_regions.is_empty(),
        "dust is not supported with rotating regions yet",
    )?;
    let grid = Grid::new(case, mesh)?;
    require(
        fields.velocity.len() == mesh.cell_count(),
        "dust: velocity does not match the mesh",
    )?;
    let (sources, total) = starts(&grid, mesh, dust, fields)?;
    let unit = match dust.release {
        DustRelease::Inflow | DustRelease::Box => "kg/s",
        DustRelease::Volume => "kg",
    };
    let mut cumulative = Vec::with_capacity(sources.len());
    let mut sum = 0.0;
    for (_, w) in &sources {
        sum += w;
        cumulative.push(sum);
    }
    let rho_f = case.config.density;
    let nu = case.config.kinematic_viscosity;
    require(
        rho_f > 0.0 && nu > 0.0,
        "dust needs a positive air density and viscosity",
    )?;
    let mut rng = Mt19937_64::new(dust.seed);
    let mut normal = StdNormal::new(0.0, 1.0);
    let mut csv =
        String::from("class,diameter_m,x_m,y_m,z_m,target,amount,surface_temperature_k,time_s\n");
    let mut classes = Vec::new();
    let mut totals: BTreeMap<String, f64> = BTreeMap::new();
    let (mut entering, mut deposited_amount) = (0.0, 0.0);
    let (mut deposited_n, mut escaped_n, mut airborne_n, mut all_n) =
        (0usize, 0usize, 0usize, 0usize);
    for (ci, class) in dust.classes.iter().enumerate() {
        let physics = Physics::new(
            dust,
            class.diameter_m,
            rho_f,
            nu,
            fields.gas_conductivity_w_per_m_k,
        );
        let amount = class.concentration_kg_per_m3 * total;
        let each = amount / dust.particles_per_class as f64;
        entering += amount;
        let mut per_target: BTreeMap<String, (usize, f64)> = BTreeMap::new();
        let (mut escaped, mut airborne, mut times) = (0usize, 0usize, Vec::new());
        if sources.is_empty() || !(sum > 0.0) {
            classes.push(serde_json::json!({
                "diameter_m": class.diameter_m, "released": 0, "deposited": 0, "escaped": 0, "airborne": 0,
                "amount": 0.0, "relaxation_time_s": physics.tau,
                "settling_speed_m_per_s": (physics.gravity * physics.tau).norm(), "targets": {}
            }));
            continue;
        }
        for _ in 0..dust.particles_per_class {
            let r = rng.uniform(0.0, sum);
            let (index, _) = sources
                .get(pick(&cumulative, r))
                .copied()
                .unwrap_or((0, 0.0));
            let (start, cell) = match dust.release {
                DustRelease::Inflow | DustRelease::Box => {
                    let Some(face) = mesh.faces().get(index) else {
                        continue;
                    };
                    let axis = dominant_axis(face.area);
                    // The cell the dust enters: the owner of an inflow boundary face, or the
                    // box's inside cell.
                    let into = match (dust.release, face.neighbour) {
                        (DustRelease::Box, Some(nb)) => {
                            let flux = fields.face_flux.get(index).copied().unwrap_or(0.0);
                            if flux > 0.0 {
                                nb
                            } else {
                                face.owner
                            }
                        }
                        _ => face.owner,
                    };
                    let ijk = grid.background_ijk(into);
                    let mut s = face.center;
                    for a in 0..3 {
                        let lo = grid.nodes[a].get(ijk[a]).copied().unwrap_or(0.0);
                        let hi = grid.nodes[a].get(ijk[a] + 1).copied().unwrap_or(0.0);
                        let w = hi - lo;
                        let c = if a == axis {
                            // Just inside the cell.
                            let centre = grid.centres[a].get(ijk[a]).copied().unwrap_or(0.0);
                            let inward = if centre < face.center.component(a) {
                                -1.0
                            } else {
                                1.0
                            };
                            face.center.component(a) + inward * 1e-6 * w
                        } else {
                            lo + w * rng.uniform(1e-6, 1.0 - 1e-6)
                        };
                        s = set(s, a, c);
                    }
                    (s, into)
                }
                DustRelease::Volume => {
                    let ijk = grid.background_ijk(index);
                    let mut s = Vec3::ZERO;
                    for a in 0..3 {
                        let lo = grid.nodes[a].get(ijk[a]).copied().unwrap_or(0.0);
                        let hi = grid.nodes[a].get(ijk[a] + 1).copied().unwrap_or(0.0);
                        s = set(s, a, lo + (hi - lo) * rng.uniform(1e-6, 1.0 - 1e-6));
                    }
                    (s, index)
                }
            };
            let cell = grid.locate(start).unwrap_or(cell);
            let tracked = track(
                &grid,
                dust,
                &physics,
                fields,
                &mut rng,
                &mut normal,
                start,
                cell,
            );
            all_n += 1;
            match tracked.fate {
                Fate::Deposited {
                    target,
                    at,
                    temperature_k,
                } => {
                    let _ = writeln!(
                        csv,
                        "{ci},{:e},{:e},{:e},{:e},{target},{:e},{:.4},{:.4}",
                        class.diameter_m, at.x, at.y, at.z, each, temperature_k, tracked.time_s
                    );
                    let e = per_target.entry(target.clone()).or_insert((0, 0.0));
                    e.0 += 1;
                    e.1 += each;
                    *totals.entry(target).or_insert(0.0) += each;
                    deposited_amount += each;
                    deposited_n += 1;
                    times.push(tracked.time_s);
                }
                Fate::Escaped => escaped += 1,
                Fate::Airborne => airborne += 1,
            }
        }
        escaped_n += escaped;
        airborne_n += airborne;
        let deposited: usize = per_target.values().map(|v| v.0).sum();
        times.sort_by(f64::total_cmp);
        let median = times.get(times.len() / 2).copied();
        classes.push(serde_json::json!({
            "diameter_m": class.diameter_m,
            "released": dust.particles_per_class,
            "deposited": deposited,
            "escaped": escaped,
            "airborne": airborne,
            "amount": amount,
            "relaxation_time_s": physics.tau,
            "settling_speed_m_per_s": (physics.gravity * physics.tau).norm(),
            "median_deposit_time_s": median,
            "longest_deposit_time_s": times.last().copied(),
            "targets": per_target.iter().map(|(k, v)| (k.clone(), serde_json::json!({"count": v.0, "amount": v.1}))).collect::<BTreeMap<_, _>>(),
        }));
    }
    let report = serde_json::json!({
        "schema": "cfd-solvers-rs/dust/1",
        "release": dust.release,
        "amount_unit": unit,
        "source_m3_per_s_or_m3": total,
        "entering": entering,
        "deposited": deposited_amount,
        "targets": totals,
        "classes": classes,
        "model": "Stokes drag with Cunningham slip, gravity less buoyancy, Talbot thermophoresis, Brownian motion; air velocity trilinear between cell centres, no-slip mirror behind walls",
    });
    if csv.len() > 64 << 20 {
        return Err(SolverError::InvalidInput(
            "dust deposits exceed 64 MiB".into(),
        ));
    }
    std::fs::write(
        out.join("dust_deposition.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    std::fs::write(out.join("dust_deposits.csv"), csv)?;
    m.files.push("dust_deposition.json".into());
    m.files.push("dust_deposits.csv".into());
    let n = all_n.max(1) as f64;
    match dust.release {
        DustRelease::Inflow | DustRelease::Box => {
            m.metrics.insert("dust_inflow_m3_per_s".into(), total);
            m.metrics.insert("dust_entering_kg_per_s".into(), entering);
            m.metrics
                .insert("dust_deposited_kg_per_s".into(), deposited_amount);
        }
        DustRelease::Volume => {
            m.metrics.insert("dust_air_volume_m3".into(), total);
            m.metrics.insert("dust_suspended_kg".into(), entering);
            m.metrics
                .insert("dust_deposited_kg".into(), deposited_amount);
        }
    }
    m.metrics
        .insert("dust_deposited_fraction".into(), deposited_n as f64 / n);
    m.metrics
        .insert("dust_escaped_fraction".into(), escaped_n as f64 / n);
    m.metrics
        .insert("dust_airborne_fraction".into(), airborne_n as f64 / n);
    Ok(())
}
