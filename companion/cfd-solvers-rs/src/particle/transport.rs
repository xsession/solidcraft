//! Geant4-style condensed-history particle transport: worlds of box regions with materials,
//! continuous / discrete / absorber processes, physics lists, a region BVH, sensitive
//! detectors, dose grids and SAR projection. Ported from upstream `cfd/particle/transport`
//! (MIT, cfd_solvers contributors). Index fields that are `npos` upstream are `Option<usize>`.
#![allow(clippy::needless_range_loop)]
use super::pusher::{AxisAlignedParticleBox, Vec3};
use crate::{require, Result};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum TrackStatus {
    #[default]
    Alive,
    Stopped,
    Escaped,
    Absorbed,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum StepLimiter {
    #[default]
    User,
    GeometryBoundary,
    ContinuousProcess,
    DiscreteProcess,
    EnergyCut,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum TransportProcessKind {
    #[default]
    ContinuousEnergyLoss,
    DiscreteInteraction,
    Absorber,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TransportParticleDefinition {
    pub name: String,
    pub charge_c: f64,
    pub mass_kg: f64,
}
impl Default for TransportParticleDefinition {
    fn default() -> Self {
        Self {
            name: "particle".into(),
            charge_c: 0.0,
            mass_kg: 1.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TransportMaterial {
    pub name: String,
    pub density_kg_per_m3: f64,
    pub production_cut_energy_ev: f64,
}
impl Default for TransportMaterial {
    fn default() -> Self {
        Self {
            name: "material".into(),
            density_kg_per_m3: 1.0,
            production_cut_energy_ev: 0.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TransportRegion {
    pub name: String,
    pub r#box: AxisAlignedParticleBox,
    pub material_index: usize,
    pub sensitive_detector_index: Option<usize>,
}
impl Default for TransportRegion {
    fn default() -> Self {
        Self {
            name: "region".into(),
            r#box: AxisAlignedParticleBox::default(),
            material_index: 0,
            sensitive_detector_index: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TransportProcess {
    pub name: String,
    pub kind: TransportProcessKind,
    /// `None` applies to all particles.
    pub particle_index: Option<usize>,
    /// `None` applies to all materials.
    pub material_index: Option<usize>,
    pub stopping_power_ev_per_m: f64,
    pub physical_interaction_length_m: f64,
    pub energy_loss_fraction: f64,
    pub secondary_energy_fraction: f64,
    pub absorption_threshold_ev: f64,
}
impl Default for TransportProcess {
    fn default() -> Self {
        Self {
            name: "process".into(),
            kind: TransportProcessKind::ContinuousEnergyLoss,
            particle_index: None,
            material_index: None,
            stopping_power_ev_per_m: 0.0,
            physical_interaction_length_m: 0.0,
            energy_loss_fraction: 0.0,
            secondary_energy_fraction: 0.0,
            absorption_threshold_ev: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransportTrack {
    pub position_m: Vec3,
    pub direction: Vec3,
    pub kinetic_energy_ev: f64,
    pub time_s: f64,
    pub track_length_m: f64,
    pub weight: f64,
    pub particle_index: usize,
    pub status: TrackStatus,
}
impl Default for TransportTrack {
    fn default() -> Self {
        Self {
            position_m: Vec3::default(),
            direction: Vec3::default(),
            kinetic_energy_ev: 0.0,
            time_s: 0.0,
            track_length_m: 0.0,
            weight: 1.0,
            particle_index: 0,
            status: TrackStatus::Alive,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransportSecondary {
    pub track: TransportTrack,
    pub creator_process: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransportStep {
    pub pre_position_m: Vec3,
    pub post_position_m: Vec3,
    pub step_length_m: f64,
    pub pre_energy_ev: f64,
    pub post_energy_ev: f64,
    pub energy_deposit_ev: f64,
    pub global_time_s: f64,
    pub material_index: Option<usize>,
    pub region_index: Option<usize>,
    pub limiter: StepLimiter,
    pub limiting_process_index: Option<usize>,
    pub secondaries: Vec<TransportSecondary>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransportScoring {
    pub total_energy_deposit_ev: f64,
    pub total_track_length_m: f64,
    pub steps: usize,
    pub secondaries: usize,
    pub boundary_crossings: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TransportSensitiveDetector {
    pub name: String,
    pub minimum_energy_deposit_ev: f64,
    pub record_zero_deposit_steps: bool,
}
impl Default for TransportSensitiveDetector {
    fn default() -> Self {
        Self {
            name: "detector".into(),
            minimum_energy_deposit_ev: 0.0,
            record_zero_deposit_steps: false,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransportWorld {
    pub particles: Vec<TransportParticleDefinition>,
    pub materials: Vec<TransportMaterial>,
    pub regions: Vec<TransportRegion>,
    pub processes: Vec<TransportProcess>,
    pub sensitive_detectors: Vec<TransportSensitiveDetector>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TransportPhysicsList {
    pub name: String,
    pub processes: Vec<TransportProcess>,
    /// Negative leaves material cuts unchanged.
    pub production_cut_energy_ev: f64,
}
impl Default for TransportPhysicsList {
    fn default() -> Self {
        Self {
            name: "physics-list".into(),
            processes: Vec::new(),
            production_cut_energy_ev: -1.0,
        }
    }
}

/// 64-bit LCG with open-interval uniform output (bit-identical to upstream).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TransportRandom {
    pub state: u64,
}
impl Default for TransportRandom {
    fn default() -> Self {
        Self {
            state: 0x9e37_79b9_7f4a_7c15,
        }
    }
}
impl TransportRandom {
    pub fn uniform_open01(&mut self) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let mantissa = (self.state >> 11) | 1;
        let u = mantissa as f64 * (1.0 / 9_007_199_254_740_992.0);
        u.max(f64::MIN_POSITIVE)
            .min(f64::from_bits(1.0f64.to_bits() - 1))
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransportBvhNode {
    pub r#box: AxisAlignedParticleBox,
    pub left: Option<usize>,
    pub right: Option<usize>,
    pub first: usize,
    pub count: usize,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransportRegionBvh {
    pub nodes: Vec<TransportBvhNode>,
    pub region_indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TransportHit {
    pub detector_index: usize,
    pub region_index: usize,
    pub material_index: usize,
    pub particle_index: usize,
    pub position_m: Vec3,
    pub energy_deposit_ev: f64,
    pub step_length_m: f64,
    pub global_time_s: f64,
    pub weight: f64,
    pub limiter: StepLimiter,
    pub process_name: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransportHitCollection {
    pub hits: Vec<TransportHit>,
    pub total_energy_deposit_ev: f64,
    pub total_weighted_energy_deposit_ev: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TransportDoseGrid3D {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub r#box: AxisAlignedParticleBox,
    pub density_kg_per_m3: f64,
    pub energy_deposit_ev: Vec<f64>,
    pub dose_gy: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransportConfig {
    pub max_step_m: f64,
    pub energy_cut_ev: f64,
    pub boundary_epsilon_m: f64,
    pub light_speed_m_per_s: f64,
}
impl Default for TransportConfig {
    fn default() -> Self {
        Self {
            max_step_m: 1.0,
            energy_cut_ev: 0.0,
            boundary_epsilon_m: 1.0e-12,
            light_speed_m_per_s: 299_792_458.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransportRunResult {
    pub primary: TransportTrack,
    pub steps: Vec<TransportStep>,
    pub secondaries: Vec<TransportSecondary>,
    pub scoring: TransportScoring,
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}
fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}
fn scale(a: Vec3, s: f64) -> Vec3 {
    Vec3::new(a.x * s, a.y * s, a.z * s)
}
fn finite(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

fn validate_world(w: &TransportWorld) -> Result<()> {
    require(
        !w.particles.is_empty(),
        "transport world requires at least one particle definition",
    )?;
    require(
        !w.materials.is_empty(),
        "transport world requires at least one material",
    )?;
    require(
        !w.regions.is_empty(),
        "transport world requires at least one region",
    )?;
    for p in &w.particles {
        require(
            !p.name.is_empty() && p.mass_kg > 0.0,
            "invalid transport particle definition",
        )?;
    }
    for m in &w.materials {
        require(
            !(m.name.is_empty() || m.density_kg_per_m3 < 0.0 || m.production_cut_energy_ev < 0.0),
            "invalid transport material",
        )?;
    }
    for r in &w.regions {
        require(
            r.material_index < w.materials.len(),
            "transport region references invalid material",
        )?;
        if let Some(d) = r.sensitive_detector_index {
            require(
                d < w.sensitive_detectors.len(),
                "transport region references invalid sensitive detector",
            )?;
        }
        let b = &r.r#box;
        require(
            b.maximum_m.x > b.minimum_m.x
                && b.maximum_m.y > b.minimum_m.y
                && b.maximum_m.z > b.minimum_m.z,
            "invalid transport region box",
        )?;
    }
    for d in &w.sensitive_detectors {
        require(
            !d.name.is_empty()
                && (d.minimum_energy_deposit_ev >= 0.0 || d.minimum_energy_deposit_ev.is_nan()),
            "invalid transport sensitive detector",
        )?;
    }
    for p in &w.processes {
        require(!p.name.is_empty(), "transport process requires a name")?;
        if let Some(i) = p.particle_index {
            require(
                i < w.particles.len(),
                "transport process references invalid particle",
            )?;
        }
        if let Some(i) = p.material_index {
            require(
                i < w.materials.len(),
                "transport process references invalid material",
            )?;
        }
        require(
            !(p.stopping_power_ev_per_m < 0.0
                || p.physical_interaction_length_m < 0.0
                || p.energy_loss_fraction < 0.0
                || p.secondary_energy_fraction < 0.0
                || p.absorption_threshold_ev < 0.0),
            "invalid transport process coefficient",
        )?;
    }
    Ok(())
}

fn validate_config(c: &TransportConfig) -> Result<()> {
    require(
        c.max_step_m > 0.0
            && c.boundary_epsilon_m >= 0.0
            && c.light_speed_m_per_s > 0.0
            && (c.energy_cut_ev >= 0.0 || c.energy_cut_ev.is_nan()),
        "invalid transport configuration",
    )
}

pub fn normalize_direction(d: Vec3) -> Result<Vec3> {
    require(finite(d), "invalid direction vector")?;
    let n = dot(d, d).sqrt();
    require(n > 0.0, "zero transport direction")?;
    Ok(Vec3::new(d.x / n, d.y / n, d.z / n))
}

pub fn contains(b: &AxisAlignedParticleBox, p: Vec3, tolerance_m: f64) -> bool {
    p.x >= b.minimum_m.x - tolerance_m
        && p.x <= b.maximum_m.x + tolerance_m
        && p.y >= b.minimum_m.y - tolerance_m
        && p.y <= b.maximum_m.y + tolerance_m
        && p.z >= b.minimum_m.z - tolerance_m
        && p.z <= b.maximum_m.z + tolerance_m
}

/// First region (declaration order) containing the point.
pub fn locate_region(world: &TransportWorld, position: Vec3, tolerance_m: f64) -> Option<usize> {
    world
        .regions
        .iter()
        .position(|r| contains(&r.r#box, position, tolerance_m))
}

pub fn distance_to_box_boundary(b: &AxisAlignedParticleBox, p: Vec3, d: Vec3) -> Result<f64> {
    let d = normalize_direction(d)?;
    let mut distance = f64::INFINITY;
    let mut candidate = |dir: f64, pos: f64, plane: f64| {
        if dir != 0.0 {
            let t = (plane - pos) / dir;
            if t >= 0.0 {
                distance = distance.min(t);
            }
        }
    };
    candidate(
        d.x,
        p.x,
        if d.x > 0.0 {
            b.maximum_m.x
        } else {
            b.minimum_m.x
        },
    );
    candidate(
        d.y,
        p.y,
        if d.y > 0.0 {
            b.maximum_m.y
        } else {
            b.minimum_m.y
        },
    );
    candidate(
        d.z,
        p.z,
        if d.z > 0.0 {
            b.maximum_m.z
        } else {
            b.minimum_m.z
        },
    );
    Ok(distance)
}

pub fn process_applies(
    p: &TransportProcess,
    track: &TransportTrack,
    material_index: usize,
) -> bool {
    p.particle_index.is_none_or(|i| i == track.particle_index)
        && p.material_index.is_none_or(|i| i == material_index)
}

pub fn apply_physics_list(
    mut world: TransportWorld,
    list: &TransportPhysicsList,
) -> Result<TransportWorld> {
    require(
        !list.name.is_empty(),
        "transport physics list requires a name",
    )?;
    if list.production_cut_energy_ev >= 0.0 {
        for m in &mut world.materials {
            m.production_cut_energy_ev = list.production_cut_energy_ev;
        }
    }
    world.processes.extend(list.processes.iter().cloned());
    validate_world(&world)?;
    Ok(world)
}

pub fn sample_exponential_interaction_length_m(
    mean_free_path_m: f64,
    rng: &mut TransportRandom,
) -> Result<f64> {
    require(
        mean_free_path_m > 0.0 && mean_free_path_m.is_finite(),
        "invalid mean free path",
    )?;
    Ok(-mean_free_path_m * rng.uniform_open01().ln())
}

fn union_box(a: &AxisAlignedParticleBox, b: &AxisAlignedParticleBox) -> AxisAlignedParticleBox {
    AxisAlignedParticleBox {
        minimum_m: Vec3::new(
            a.minimum_m.x.min(b.minimum_m.x),
            a.minimum_m.y.min(b.minimum_m.y),
            a.minimum_m.z.min(b.minimum_m.z),
        ),
        maximum_m: Vec3::new(
            a.maximum_m.x.max(b.maximum_m.x),
            a.maximum_m.y.max(b.maximum_m.y),
            a.maximum_m.z.max(b.maximum_m.z),
        ),
    }
}

fn center(b: &AxisAlignedParticleBox) -> Vec3 {
    scale(add(b.minimum_m, b.maximum_m), 0.5)
}

/// Median-split BVH over region boxes. The split uses `select_nth_unstable_by`, so the
/// leaf ordering may differ from upstream's `std::nth_element`; queries agree for
/// non-overlapping regions.
pub fn build_region_bvh(world: &TransportWorld, leaf_size: usize) -> Result<TransportRegionBvh> {
    validate_world(world)?;
    require(leaf_size != 0, "BVH leaf size must be positive")?;
    let mut bvh = TransportRegionBvh {
        nodes: Vec::new(),
        region_indices: (0..world.regions.len()).collect(),
    };
    fn build(
        world: &TransportWorld,
        bvh: &mut TransportRegionBvh,
        first: usize,
        count: usize,
        leaf: usize,
    ) -> usize {
        let mut bounds = world.regions[bvh.region_indices[first]].r#box;
        for i in 1..count {
            bounds = union_box(&bounds, &world.regions[bvh.region_indices[first + i]].r#box);
        }
        let node_index = bvh.nodes.len();
        bvh.nodes.push(TransportBvhNode {
            r#box: bounds,
            left: None,
            right: None,
            first,
            count,
        });
        if count > leaf {
            let sx = bounds.maximum_m.x - bounds.minimum_m.x;
            let sy = bounds.maximum_m.y - bounds.minimum_m.y;
            let sz = bounds.maximum_m.z - bounds.minimum_m.z;
            let axis = if sx >= sy && sx >= sz {
                0
            } else if sy >= sz {
                1
            } else {
                2
            };
            let coord = |r: usize| {
                let c = center(&world.regions[r].r#box);
                match axis {
                    0 => c.x,
                    1 => c.y,
                    _ => c.z,
                }
            };
            let mid = first + count / 2;
            bvh.region_indices[first..first + count]
                .select_nth_unstable_by(mid - first, |&a, &b| coord(a).total_cmp(&coord(b)));
            let l = build(world, bvh, first, mid - first, leaf);
            let r = build(world, bvh, mid, first + count - mid, leaf);
            let n = &mut bvh.nodes[node_index];
            n.left = Some(l);
            n.right = Some(r);
            n.count = 0;
        }
        node_index
    }
    if !bvh.region_indices.is_empty() {
        let n = bvh.region_indices.len();
        build(world, &mut bvh, 0, n, leaf_size);
    }
    Ok(bvh)
}

pub fn locate_region_bvh(
    world: &TransportWorld,
    bvh: &TransportRegionBvh,
    position: Vec3,
    tolerance_m: f64,
) -> Option<usize> {
    fn query(
        world: &TransportWorld,
        bvh: &TransportRegionBvh,
        node: usize,
        p: Vec3,
        tol: f64,
    ) -> Option<usize> {
        let n = bvh.nodes.get(node)?;
        if !contains(&n.r#box, p, tol) {
            return None;
        }
        if n.count > 0 {
            return (0..n.count)
                .filter_map(|i| bvh.region_indices.get(n.first + i).copied())
                .find(|&r| {
                    world
                        .regions
                        .get(r)
                        .is_some_and(|reg| contains(&reg.r#box, p, tol))
                });
        }
        n.left
            .and_then(|l| query(world, bvh, l, p, tol))
            .or_else(|| n.right.and_then(|r| query(world, bvh, r, p, tol)))
    }
    if bvh.nodes.is_empty() {
        return None;
    }
    query(world, bvh, 0, position, tolerance_m)
}

pub fn collect_transport_hits(
    world: &TransportWorld,
    steps: &[TransportStep],
) -> Result<TransportHitCollection> {
    validate_world(world)?;
    let mut c = TransportHitCollection::default();
    for s in steps {
        let (Some(region), Some(material)) = (s.region_index, s.material_index) else {
            continue;
        };
        if region >= world.regions.len() || material >= world.materials.len() {
            continue;
        }
        let Some(detector) = world.regions[region].sensitive_detector_index else {
            continue;
        };
        let Some(def) = world.sensitive_detectors.get(detector) else {
            continue;
        };
        if s.energy_deposit_ev < def.minimum_energy_deposit_ev
            && !(def.record_zero_deposit_steps && s.step_length_m > 0.0)
        {
            continue;
        }
        let hit = TransportHit {
            detector_index: detector,
            region_index: region,
            material_index: material,
            particle_index: 0,
            position_m: scale(add(s.pre_position_m, s.post_position_m), 0.5),
            energy_deposit_ev: s.energy_deposit_ev,
            step_length_m: s.step_length_m,
            global_time_s: s.global_time_s,
            weight: 1.0,
            limiter: s.limiter,
            process_name: s
                .limiting_process_index
                .and_then(|i| world.processes.get(i))
                .map(|p| p.name.clone())
                .unwrap_or_default(),
        };
        c.total_energy_deposit_ev += hit.energy_deposit_ev;
        c.total_weighted_energy_deposit_ev += hit.energy_deposit_ev * hit.weight;
        c.hits.push(hit);
    }
    Ok(c)
}

pub fn make_transport_dose_grid(
    nx: usize,
    ny: usize,
    nz: usize,
    b: AxisAlignedParticleBox,
    density_kg_per_m3: f64,
) -> Result<TransportDoseGrid3D> {
    require(
        nx != 0
            && ny != 0
            && nz != 0
            && density_kg_per_m3 > 0.0
            && contains(&b, b.minimum_m, 0.0)
            && b.maximum_m.x > b.minimum_m.x
            && b.maximum_m.y > b.minimum_m.y
            && b.maximum_m.z > b.minimum_m.z,
        "invalid transport dose grid",
    )?;
    let n = nx * ny * nz;
    Ok(TransportDoseGrid3D {
        nx,
        ny,
        nz,
        r#box: b,
        density_kg_per_m3,
        energy_deposit_ev: vec![0.0; n],
        dose_gy: vec![0.0; n],
    })
}

pub fn score_hits_to_dose_grid(
    hits: &TransportHitCollection,
    g: &mut TransportDoseGrid3D,
) -> Result<()> {
    require(
        g.nx != 0
            && g.ny != 0
            && g.nz != 0
            && g.energy_deposit_ev.len() == g.nx * g.ny * g.nz
            && g.dose_gy.len() == g.energy_deposit_ev.len()
            && g.density_kg_per_m3 > 0.0,
        "invalid transport dose grid storage",
    )?;
    let lx = g.r#box.maximum_m.x - g.r#box.minimum_m.x;
    let ly = g.r#box.maximum_m.y - g.r#box.minimum_m.y;
    let lz = g.r#box.maximum_m.z - g.r#box.minimum_m.z;
    require(
        lx > 0.0 && ly > 0.0 && lz > 0.0,
        "invalid transport dose grid box",
    )?;
    for h in &hits.hits {
        if !contains(&g.r#box, h.position_m, 0.0) {
            continue;
        }
        let fx = (h.position_m.x - g.r#box.minimum_m.x) / lx;
        let fy = (h.position_m.y - g.r#box.minimum_m.y) / ly;
        let fz = (h.position_m.z - g.r#box.minimum_m.z) / lz;
        let ix = (g.nx - 1).min((fx * g.nx as f64).floor() as usize);
        let iy = (g.ny - 1).min((fy * g.ny as f64).floor() as usize);
        let iz = (g.nz - 1).min((fz * g.nz as f64).floor() as usize);
        g.energy_deposit_ev[(iz * g.ny + iy) * g.nx + ix] += h.energy_deposit_ev * h.weight;
    }
    const ELECTRON_VOLT_J: f64 = 1.602_176_634e-19;
    let voxel_mass =
        g.density_kg_per_m3 * (lx / g.nx as f64) * (ly / g.ny as f64) * (lz / g.nz as f64);
    for i in 0..g.energy_deposit_ev.len() {
        g.dose_gy[i] = g.energy_deposit_ev[i] * ELECTRON_VOLT_J / voxel_mass;
    }
    Ok(())
}

pub fn dose_rate_to_sar_w_per_kg(
    g: &TransportDoseGrid3D,
    exposure_time_s: f64,
) -> Result<Vec<f64>> {
    require(
        exposure_time_s > 0.0 && !g.dose_gy.is_empty(),
        "invalid dose-rate conversion",
    )?;
    g.dose_gy
        .iter()
        .map(|&d| {
            require(d >= 0.0 && d.is_finite(), "invalid dose grid value")?;
            Ok(d / exposure_time_s)
        })
        .collect()
}

/// Averages voxel SAR columns onto an `nx × ny` Pennes-2D grid.
pub fn project_dose_grid_to_pennes2d_sar(
    g: &TransportDoseGrid3D,
    exposure_time_s: f64,
    nx: usize,
    ny: usize,
) -> Result<Vec<f64>> {
    require(nx != 0 && ny != 0, "invalid Pennes projection size")?;
    let sar3 = dose_rate_to_sar_w_per_kg(g, exposure_time_s)?;
    require(
        sar3.len() == g.nx * g.ny * g.nz,
        "invalid transport dose grid storage",
    )?;
    let mut out = vec![0.0; nx * ny];
    let mut count = vec![0.0; nx * ny];
    for iz in 0..g.nz {
        for iy in 0..g.ny {
            for ix in 0..g.nx {
                let src = (iz * g.ny + iy) * g.nx + ix;
                let tx = (nx - 1).min(ix * nx / g.nx);
                let ty = (ny - 1).min(iy * ny / g.ny);
                out[ty * nx + tx] += sar3[src];
                count[ty * nx + tx] += 1.0;
            }
        }
    }
    for i in 0..out.len() {
        if count[i] > 0.0 {
            out[i] /= count[i];
        }
    }
    Ok(out)
}

/// One condensed-history step. Returns a default step when the track is not alive or has
/// left every region (the latter marks it escaped).
pub fn transport_one_step(
    world: &TransportWorld,
    track: &mut TransportTrack,
    config: &TransportConfig,
) -> Result<TransportStep> {
    validate_world(world)?;
    validate_config(config)?;
    require(
        track.particle_index < world.particles.len(),
        "track references invalid particle",
    )?;
    if track.status != TrackStatus::Alive {
        return Ok(TransportStep::default());
    }
    track.direction = normalize_direction(track.direction)?;
    require(
        track.kinetic_energy_ev >= 0.0 && finite(track.position_m) && track.weight.is_finite(),
        "invalid transport track",
    )?;
    let Some(region_index) = locate_region(world, track.position_m, config.boundary_epsilon_m)
    else {
        track.status = TrackStatus::Escaped;
        return Ok(TransportStep::default());
    };
    let region = &world.regions[region_index];
    let material_index = region.material_index;
    let material = &world.materials[material_index];

    let mut step = TransportStep {
        pre_position_m: track.position_m,
        pre_energy_ev: track.kinetic_energy_ev,
        material_index: Some(material_index),
        region_index: Some(region_index),
        global_time_s: track.time_s,
        limiter: StepLimiter::User,
        ..Default::default()
    };
    let mut selected = config.max_step_m;
    let boundary = distance_to_box_boundary(&region.r#box, track.position_m, track.direction)?;
    if boundary < selected {
        selected = boundary;
        step.limiter = StepLimiter::GeometryBoundary;
    }
    for (i, p) in world.processes.iter().enumerate() {
        if !process_applies(p, track, material_index) {
            continue;
        }
        if p.kind == TransportProcessKind::DiscreteInteraction
            && p.physical_interaction_length_m > 0.0
            && p.physical_interaction_length_m < selected
        {
            selected = p.physical_interaction_length_m;
            step.limiter = StepLimiter::DiscreteProcess;
            step.limiting_process_index = Some(i);
        }
        if p.kind == TransportProcessKind::Absorber
            && track.kinetic_energy_ev <= p.absorption_threshold_ev
        {
            selected = 0.0;
            step.limiter = StepLimiter::EnergyCut;
            step.limiting_process_index = Some(i);
        }
    }
    selected = selected.max(0.0);
    step.step_length_m = selected;
    step.post_position_m = add(track.position_m, scale(track.direction, selected));

    let mut continuous = 0.0;
    for p in &world.processes {
        if p.kind == TransportProcessKind::ContinuousEnergyLoss
            && process_applies(p, track, material_index)
        {
            continuous += p.stopping_power_ev_per_m * selected;
        }
    }
    continuous = track.kinetic_energy_ev.min(continuous);
    track.kinetic_energy_ev -= continuous;
    step.energy_deposit_ev += continuous;

    if step.limiter == StepLimiter::DiscreteProcess {
        if let Some(p) = step
            .limiting_process_index
            .and_then(|i| world.processes.get(i))
        {
            let loss = track
                .kinetic_energy_ev
                .min(track.kinetic_energy_ev * p.energy_loss_fraction.min(1.0));
            track.kinetic_energy_ev -= loss;
            step.energy_deposit_ev += loss;
            let secondary_energy = track.kinetic_energy_ev * p.secondary_energy_fraction;
            if secondary_energy > material.production_cut_energy_ev && secondary_energy > 0.0 {
                let mut st = *track;
                st.position_m = step.post_position_m;
                st.kinetic_energy_ev = secondary_energy;
                st.track_length_m = 0.0;
                st.time_s = track.time_s + selected / config.light_speed_m_per_s;
                st.weight = track.weight;
                st.status = TrackStatus::Alive;
                step.secondaries.push(TransportSecondary {
                    track: st,
                    creator_process: p.name.clone(),
                });
                track.kinetic_energy_ev = (track.kinetic_energy_ev - secondary_energy).max(0.0);
            }
        }
    }

    track.position_m = step.post_position_m;
    track.track_length_m += selected;
    track.time_s += selected / config.light_speed_m_per_s;
    step.global_time_s = track.time_s;
    step.post_energy_ev = track.kinetic_energy_ev;
    if step.limiter == StepLimiter::EnergyCut
        || track.kinetic_energy_ev <= config.energy_cut_ev.max(material.production_cut_energy_ev)
    {
        track.status = TrackStatus::Stopped;
    }
    if step.limiter == StepLimiter::GeometryBoundary {
        let nudged = add(
            track.position_m,
            scale(track.direction, config.boundary_epsilon_m * 4.0),
        );
        if locate_region(world, nudged, config.boundary_epsilon_m).is_none() {
            track.status = TrackStatus::Escaped;
        } else {
            track.position_m = nudged;
        }
    }
    Ok(step)
}

fn run(
    world: &TransportWorld,
    mut primary: TransportTrack,
    config: &TransportConfig,
    max_steps: usize,
    mut one: impl FnMut(&mut TransportTrack) -> Result<TransportStep>,
) -> Result<TransportRunResult> {
    validate_world(world)?;
    validate_config(config)?;
    require(max_steps != 0, "transport requires positive max steps")?;
    primary.direction = normalize_direction(primary.direction)?;
    let mut r = TransportRunResult {
        primary,
        ..Default::default()
    };
    for _ in 0..max_steps {
        if r.primary.status != TrackStatus::Alive {
            break;
        }
        let step = one(&mut r.primary)?;
        if step.step_length_m == 0.0
            && step.secondaries.is_empty()
            && r.primary.status == TrackStatus::Alive
        {
            break;
        }
        r.scoring.total_energy_deposit_ev += step.energy_deposit_ev;
        r.scoring.total_track_length_m += step.step_length_m;
        r.scoring.secondaries += step.secondaries.len();
        if step.limiter == StepLimiter::GeometryBoundary {
            r.scoring.boundary_crossings += 1;
        }
        r.scoring.steps += 1;
        r.secondaries.extend(step.secondaries.iter().cloned());
        r.steps.push(step);
    }
    Ok(r)
}

pub fn transport_track(
    world: &TransportWorld,
    primary: TransportTrack,
    config: &TransportConfig,
    max_steps: usize,
) -> Result<TransportRunResult> {
    run(world, primary, config, max_steps, |t| {
        transport_one_step(world, t, config)
    })
}

/// Resamples every discrete interaction length from an exponential distribution, then steps.
pub fn transport_one_step_stochastic(
    world: &TransportWorld,
    track: &mut TransportTrack,
    config: &TransportConfig,
    rng: &mut TransportRandom,
) -> Result<TransportStep> {
    let mut sampled = world.clone();
    for p in &mut sampled.processes {
        if p.kind == TransportProcessKind::DiscreteInteraction
            && p.physical_interaction_length_m > 0.0
        {
            p.physical_interaction_length_m =
                sample_exponential_interaction_length_m(p.physical_interaction_length_m, rng)?;
        }
    }
    transport_one_step(&sampled, track, config)
}

pub fn transport_track_stochastic(
    world: &TransportWorld,
    primary: TransportTrack,
    config: &TransportConfig,
    rng: &mut TransportRandom,
    max_steps: usize,
) -> Result<TransportRunResult> {
    run(world, primary, config, max_steps, |t| {
        transport_one_step_stochastic(world, t, config, rng)
    })
}
