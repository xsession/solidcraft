//! Geometry-based simulation setup records: persistent named selections of CAD faces or
//! bodies, typed boundary assignments, thermal contacts and oriented conduction, and the
//! coverage audit that decides whether a setup may be applied or run.
//!
//! This module is portable. It never evaluates geometry. The document/kernel layer supplies
//! a [`TopologySnapshot`] of the exact current faces and bodies; references are re-resolved
//! against it by stable identity only. A missing or ambiguously split reference is reported
//! as unresolved and blocks Apply; it is never reattached to a nearby face. Display
//! triangles, screen positions and face ordinals are not identities.
//!
//! These records are not yet consumed by a numerical adapter. A complete, conflict-free
//! audit is preparation evidence, not mesh validity, solver admission or physical validation.
use std::collections::{BTreeMap, BTreeSet};

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{Error, Result};

pub const MAX_SELECTIONS: usize = 256;
pub const MAX_MEMBERS: usize = 4096;
pub const MAX_ASSIGNMENTS: usize = 256;
pub const MAX_SNAPSHOT_ENTITIES: usize = 200_000;
/// The only coordinate frame defined so far. Named frames remain planned.
pub const GLOBAL_FRAME: &str = "global";

// ---------------------------------------------------------------- references and selections

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum EntityKind {
    Face,
    Body,
}

/// A stable CAD identity: occurrence path, body identity and, for faces, the persistent
/// face name (a lineage name, or an exact split-piece name chosen during review).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EntityRef {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub occurrence: Vec<String>,
    pub body: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub face: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NamedSelection {
    pub id: String,
    pub name: String,
    pub kind: EntityKind,
    pub members: Vec<EntityRef>,
}

// ---------------------------------------------------------------- boundary conditions

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "spec", rename_all = "snake_case", deny_unknown_fields)]
pub enum InletSpec {
    /// Positive speed along each face's inward fluid normal: u(x) = -U n(x).
    NormalSpeed { inward_m_per_s: f64 },
    /// A fixed vector in a named frame. It does not rotate with a curved face's normal.
    Velocity { vector_m_per_s: [f64; 3], frame: String },
    /// Volume flow integrated over the whole selected area.
    VolumeFlow { m3_per_s: f64 },
    /// Mass flow needs an explicit density to become a volume flow.
    MassFlow { kg_per_s: f64, density_kg_per_m3: f64 },
}

impl InletSpec {
    /// Implied volume flow over a resolved area; `None` when it depends on local normals.
    pub fn volume_flow(&self, area_m2: f64) -> Option<f64> {
        match self {
            Self::NormalSpeed { inward_m_per_s } => Some(inward_m_per_s * area_m2),
            Self::VolumeFlow { m3_per_s } => Some(*m3_per_s),
            Self::MassFlow { kg_per_s, density_kg_per_m3 } => Some(kg_per_s / density_kg_per_m3),
            Self::Velocity { .. } => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum PressureKind {
    Static,
    Total,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Backflow {
    Prevent,
    Allow,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum WallSlip {
    NoSlip,
    Slip,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "motion", rename_all = "snake_case", deny_unknown_fields)]
pub enum WallMotion {
    Stationary,
    Translating {
        velocity_m_per_s: [f64; 3],
        frame: String,
    },
    /// Surface rotation about an axis. This is not a rotating fluid reference frame.
    Rotating {
        origin_m: [f64; 3],
        axis: [f64; 3],
        rad_per_s: f64,
        frame: String,
    },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "thermal", rename_all = "snake_case", deny_unknown_fields)]
pub enum WallThermal {
    Adiabatic,
    Temperature { kelvin: f64 },
    HeatFlux { w_per_m2: f64 },
    Convection { w_per_m2_k: f64, ambient_kelvin: f64 },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WallCondition {
    pub slip: WallSlip,
    pub motion: WallMotion,
    pub thermal: WallThermal,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum BoundaryCondition {
    Inlet {
        inflow: InletSpec,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        temperature_kelvin: Option<f64>,
    },
    /// Predominantly outflowing pressure boundary.
    Outlet {
        gauge_pressure_pa: f64,
        pressure: PressureKind,
        backflow: Backflow,
    },
    /// Explicitly two-way pressure boundary; backflow carries its own temperature.
    Opening {
        gauge_pressure_pa: f64,
        pressure: PressureKind,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        backflow_temperature_kelvin: Option<f64>,
    },
    Wall(WallCondition),
    Symmetry,
}

impl BoundaryCondition {
    fn label(&self) -> &'static str {
        match self {
            Self::Inlet { .. } => "inlet",
            Self::Outlet { .. } => "outlet",
            Self::Opening { .. } => "opening",
            Self::Wall(_) => "wall",
            Self::Symmetry => "symmetry",
        }
    }
    fn sets_pressure(&self) -> bool {
        matches!(self, Self::Outlet { .. } | Self::Opening { .. })
    }
    fn admits_inflow(&self) -> bool {
        matches!(self, Self::Inlet { .. } | Self::Opening { .. })
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoundaryAssignment {
    pub id: String,
    pub selection: String,
    pub condition: BoundaryCondition,
}

// ---------------------------------------------------------------- contact and conduction

/// One interface resistance law. Exactly one is stored per physical contact.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "law", rename_all = "snake_case", deny_unknown_fields)]
pub enum ContactResistance {
    /// Resistance of the whole contact, K/W.
    Total { kelvin_per_watt: f64 },
    /// Area-normalized resistance, m²·K/W.
    AreaSpecific { m2_kelvin_per_watt: f64 },
    /// Area conductance, W/(m²·K): the reciprocal of area-specific resistance.
    Conductance { w_per_m2_kelvin: f64 },
    /// Homogeneous one-dimensional layer: thickness / conductivity.
    ThinLayer { thickness_m: f64, conductivity_w_per_m_kelvin: f64 },
}

impl ContactResistance {
    fn check(&self) -> Result<()> {
        let ok = match self {
            Self::Total { kelvin_per_watt: v } | Self::AreaSpecific { m2_kelvin_per_watt: v } | Self::Conductance { w_per_m2_kelvin: v } => {
                positive(*v)
            }
            Self::ThinLayer { thickness_m, conductivity_w_per_m_kelvin } => positive(*thickness_m) && positive(*conductivity_w_per_m_kelvin),
        };
        if ok { Ok(()) } else { Err(Error::Invalid("contact resistance values must be finite and positive".into())) }
    }

    /// Area-specific resistance (m²·K/W) for a contact of the given effective area.
    pub fn area_specific(&self, area_m2: f64) -> Result<f64> {
        self.check()?;
        if !positive(area_m2) {
            return Err(Error::Invalid("contact area must be finite and positive".into()));
        }
        Ok(match self {
            Self::Total { kelvin_per_watt } => kelvin_per_watt * area_m2,
            Self::AreaSpecific { m2_kelvin_per_watt } => *m2_kelvin_per_watt,
            Self::Conductance { w_per_m2_kelvin } => 1.0 / w_per_m2_kelvin,
            Self::ThinLayer { thickness_m, conductivity_w_per_m_kelvin } => thickness_m / conductivity_w_per_m_kelvin,
        })
    }

    /// Total resistance (K/W) for a contact of the given effective area.
    pub fn total(&self, area_m2: f64) -> Result<f64> {
        Ok(self.area_specific(area_m2)? / area_m2)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ContactInterface {
    pub id: String,
    pub side_a: String,
    pub side_b: String,
    pub resistance: ContactResistance,
}

/// Conductivity of a solid region. Orthotropic axes are the material's local principal
/// axes expressed in the global frame (one row per axis), so several occurrences can share
/// a material while each assignment keeps its own orientation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "model", rename_all = "snake_case", deny_unknown_fields)]
pub enum Conductivity {
    Isotropic { w_per_m_kelvin: f64 },
    Orthotropic { principal_w_per_m_kelvin: [f64; 3], axes: [[f64; 3]; 3] },
}

pub type Tensor = [[f64; 3]; 3];

impl Conductivity {
    /// World-frame tensor K_world = R diag(k) Rᵀ, where R's columns are the principal axes.
    pub fn world_tensor(&self) -> Result<Tensor> {
        match self {
            Self::Isotropic { w_per_m_kelvin: k } => {
                if !positive(*k) {
                    return Err(Error::Invalid("conductivity must be finite and positive".into()));
                }
                Ok([[*k, 0.0, 0.0], [0.0, *k, 0.0], [0.0, 0.0, *k]])
            }
            Self::Orthotropic { principal_w_per_m_kelvin, axes } => {
                if !principal_w_per_m_kelvin.iter().all(|k| positive(*k)) {
                    return Err(Error::Invalid("principal conductivities must be finite and positive".into()));
                }
                check_axes(axes)?;
                let mut tensor = [[0.0; 3]; 3];
                for (k, axis) in principal_w_per_m_kelvin.iter().zip(axes) {
                    for (row, a) in tensor.iter_mut().zip(axis) {
                        for (cell, b) in row.iter_mut().zip(axis) {
                            *cell += k * a * b;
                        }
                    }
                }
                Ok(tensor)
            }
        }
    }
}

/// Whether a tensor is symmetric positive definite (Sylvester's criterion).
pub fn is_symmetric_positive_definite(t: &Tensor) -> bool {
    let [[a, b, c], [d, e, f], [g, h, i]] = *t;
    let finite = t.iter().flatten().all(|v| v.is_finite());
    let scale = t.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs())).max(f64::MIN_POSITIVE);
    let symmetric = (b - d).abs() <= 1e-12 * scale && (c - g).abs() <= 1e-12 * scale && (f - h).abs() <= 1e-12 * scale;
    let minor2 = a * e - b * d;
    let det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    finite && symmetric && a > 0.0 && minor2 > 0.0 && det > 0.0
}

fn check_axes(axes: &[[f64; 3]; 3]) -> Result<()> {
    let dot = |u: &[f64; 3], v: &[f64; 3]| u.iter().zip(v).map(|(a, b)| a * b).sum::<f64>();
    let [x, y, z] = axes;
    let cross = [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]];
    let unit = axes.iter().all(|a| a.iter().all(|v| v.is_finite()) && (dot(a, a) - 1.0).abs() <= 1e-9);
    let orthogonal = dot(x, y).abs() <= 1e-9 && dot(y, z).abs() <= 1e-9 && dot(x, z).abs() <= 1e-9;
    if !(unit && orthogonal && dot(&cross, z) > 0.0) {
        return Err(Error::Invalid("material axes must be a right-handed orthonormal set".into()));
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ConductionAssignment {
    pub id: String,
    /// A body selection.
    pub selection: String,
    /// Optional design material name the values were taken from; not a live binding.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub material: Option<String>,
    pub conductivity: Conductivity,
}

// ---------------------------------------------------------------- the saved setup

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Setup {
    /// Geometry revision the setup was last applied against.
    pub geometry_revision: u64,
    pub selections: Vec<NamedSelection>,
    pub boundaries: Vec<BoundaryAssignment>,
    pub contacts: Vec<ContactInterface>,
    pub conduction: Vec<ConductionAssignment>,
    /// Explicitly chosen wall for otherwise unassigned fluid boundary faces. Faces using it
    /// are reported as inherited, never as reviewed assignments.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub default_wall: Option<WallCondition>,
}

impl Setup {
    pub fn is_empty(&self) -> bool {
        *self == Self::default()
    }

    /// Structural limits and internal references only; geometry is checked by [`audit`].
    pub fn validate(&self) -> Result<()> {
        if self.selections.len() > MAX_SELECTIONS
            || [self.boundaries.len(), self.contacts.len(), self.conduction.len()].iter().any(|n| *n > MAX_ASSIGNMENTS)
        {
            return Err(Error::Limit("too many setup records".into()));
        }
        let mut ids = BTreeSet::new();
        let mut kinds = BTreeMap::new();
        for selection in &self.selections {
            text(&selection.id, "selection identifier", 80)?;
            text(&selection.name, "selection name", 160)?;
            if !ids.insert(selection.id.as_str()) {
                return Err(Error::Invalid("duplicate setup identifier".into()));
            }
            if selection.members.is_empty() || selection.members.len() > MAX_MEMBERS {
                return Err(Error::Invalid(format!("a selection needs 1..{MAX_MEMBERS} members")));
            }
            let mut members = BTreeSet::new();
            for member in &selection.members {
                validate_ref(member, selection.kind)?;
                if !members.insert(member) {
                    return Err(Error::Invalid(format!("selection {} lists a member twice", selection.id)));
                }
            }
            kinds.insert(selection.id.as_str(), selection.kind);
        }
        let expect = |id: &str, kind: EntityKind, what: &str| match kinds.get(id) {
            Some(found) if *found == kind => Ok(()),
            Some(_) => Err(Error::Invalid(format!("{what} {id} needs a {} selection", if kind == EntityKind::Face { "face" } else { "body" }))),
            None => Err(Error::Invalid(format!("{what} refers to unknown selection {id}"))),
        };
        for boundary in &self.boundaries {
            text(&boundary.id, "assignment identifier", 80)?;
            if !ids.insert(boundary.id.as_str()) {
                return Err(Error::Invalid("duplicate setup identifier".into()));
            }
            expect(&boundary.selection, EntityKind::Face, "boundary")?;
            validate_condition(&boundary.condition)?;
        }
        if let Some(wall) = &self.default_wall {
            validate_wall(wall)?;
        }
        let mut pairs = BTreeSet::new();
        for contact in &self.contacts {
            text(&contact.id, "contact identifier", 80)?;
            if !ids.insert(contact.id.as_str()) {
                return Err(Error::Invalid("duplicate setup identifier".into()));
            }
            expect(&contact.side_a, EntityKind::Face, "contact")?;
            expect(&contact.side_b, EntityKind::Face, "contact")?;
            if contact.side_a == contact.side_b {
                return Err(Error::Invalid("a contact needs two different sides".into()));
            }
            // One physical interface is counted once, whichever side was picked first.
            let pair = if contact.side_a < contact.side_b { (&contact.side_a, &contact.side_b) } else { (&contact.side_b, &contact.side_a) };
            if !pairs.insert(pair) {
                return Err(Error::Invalid("the same contact interface is defined twice; count each contact once".into()));
            }
            contact.resistance.check()?;
        }
        for conduction in &self.conduction {
            text(&conduction.id, "conduction identifier", 80)?;
            if !ids.insert(conduction.id.as_str()) {
                return Err(Error::Invalid("duplicate setup identifier".into()));
            }
            expect(&conduction.selection, EntityKind::Body, "conduction")?;
            if let Some(material) = &conduction.material {
                text(material, "material name", 160)?;
            }
            let tensor = conduction.conductivity.world_tensor()?;
            if !is_symmetric_positive_definite(&tensor) {
                return Err(Error::Invalid("conductivity must be symmetric positive definite".into()));
            }
        }
        Ok(())
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn finite3(v: &[f64; 3]) -> bool {
    v.iter().all(|x| x.is_finite())
}

fn text(value: &str, label: &str, limit: usize) -> Result<()> {
    if value.trim().is_empty() || value.len() > limit || value.chars().any(char::is_control) {
        return Err(Error::Invalid(format!("{label} must be 1..{limit} bytes without control characters")));
    }
    Ok(())
}

fn validate_ref(member: &EntityRef, kind: EntityKind) -> Result<()> {
    if member.occurrence.len() > 32 {
        return Err(Error::Limit("occurrence paths exceed 32 levels".into()));
    }
    for part in &member.occurrence {
        text(part, "occurrence identifier", 128)?;
    }
    text(&member.body, "body identifier", 128)?;
    match (kind, &member.face) {
        (EntityKind::Face, Some(face)) => text(face, "face name", 512),
        (EntityKind::Face, None) => Err(Error::Invalid("a face selection member needs a face name".into())),
        (EntityKind::Body, Some(_)) => Err(Error::Invalid("a body selection member cannot name a face".into())),
        (EntityKind::Body, None) => Ok(()),
    }
}

fn frame(name: &str) -> Result<()> {
    if name != GLOBAL_FRAME {
        return Err(Error::Invalid(format!("unknown coordinate frame {name}; only \"{GLOBAL_FRAME}\" is defined")));
    }
    Ok(())
}

fn temperature(value: Option<f64>) -> Result<()> {
    if value.is_some_and(|k| !positive(k)) {
        return Err(Error::Invalid("temperatures are absolute, in kelvin, and must be positive".into()));
    }
    Ok(())
}

fn validate_wall(wall: &WallCondition) -> Result<()> {
    match &wall.motion {
        WallMotion::Stationary => {}
        WallMotion::Translating { velocity_m_per_s, frame: name } => {
            frame(name)?;
            if !finite3(velocity_m_per_s) {
                return Err(Error::Invalid("wall velocity must be finite".into()));
            }
        }
        WallMotion::Rotating { origin_m, axis, rad_per_s, frame: name } => {
            frame(name)?;
            let length = axis.iter().map(|v| v * v).sum::<f64>().sqrt();
            if !finite3(origin_m) || !positive(length) || !rad_per_s.is_finite() {
                return Err(Error::Invalid("a rotating wall needs a finite origin, a nonzero axis and a finite speed".into()));
            }
        }
    }
    match wall.thermal {
        WallThermal::Adiabatic => Ok(()),
        WallThermal::Temperature { kelvin } => temperature(Some(kelvin)),
        WallThermal::HeatFlux { w_per_m2 } if w_per_m2.is_finite() => Ok(()),
        WallThermal::HeatFlux { .. } => Err(Error::Invalid("wall heat flux must be finite".into())),
        WallThermal::Convection { w_per_m2_k, ambient_kelvin } if positive(w_per_m2_k) => temperature(Some(ambient_kelvin)),
        WallThermal::Convection { .. } => Err(Error::Invalid("a convection coefficient must be finite and positive".into())),
    }
}

fn validate_condition(condition: &BoundaryCondition) -> Result<()> {
    match condition {
        BoundaryCondition::Inlet { inflow, temperature_kelvin } => {
            temperature(*temperature_kelvin)?;
            let ok = match inflow {
                InletSpec::NormalSpeed { inward_m_per_s } => positive(*inward_m_per_s),
                InletSpec::VolumeFlow { m3_per_s } => positive(*m3_per_s),
                InletSpec::MassFlow { kg_per_s, density_kg_per_m3 } => positive(*kg_per_s) && positive(*density_kg_per_m3),
                InletSpec::Velocity { vector_m_per_s, frame: name } => {
                    frame(name)?;
                    finite3(vector_m_per_s) && vector_m_per_s.iter().any(|v| *v != 0.0)
                }
            };
            if ok { Ok(()) } else { Err(Error::Invalid("an inlet needs a finite positive inflow".into())) }
        }
        BoundaryCondition::Outlet { gauge_pressure_pa, .. } if gauge_pressure_pa.is_finite() => Ok(()),
        BoundaryCondition::Opening { gauge_pressure_pa, backflow_temperature_kelvin, .. } if gauge_pressure_pa.is_finite() => {
            temperature(*backflow_temperature_kelvin)
        }
        BoundaryCondition::Outlet { .. } | BoundaryCondition::Opening { .. } => Err(Error::Invalid("boundary pressure must be finite".into())),
        BoundaryCondition::Wall(wall) => validate_wall(wall),
        BoundaryCondition::Symmetry => Ok(()),
    }
}

// ---------------------------------------------------------------- geometry snapshot

/// One exact face of the current evaluated geometry, provided by the document/kernel layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SnapshotFace {
    #[serde(default)]
    pub occurrence: Vec<String>,
    pub body: String,
    /// Exact persistent name, including any split-piece suffix.
    pub name: String,
    /// The lineage name the piece descends from (equal to `name` when unsplit).
    pub lineage: String,
    pub area_m2: f64,
    /// Whether this face bounds the prepared fluid region (and on which side it is wetted).
    #[serde(default)]
    pub fluid_boundary: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SnapshotBody {
    #[serde(default)]
    pub occurrence: Vec<String>,
    pub body: String,
    pub volume_m3: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct TopologySnapshot {
    pub geometry_revision: u64,
    pub faces: Vec<SnapshotFace>,
    pub bodies: Vec<SnapshotBody>,
}

/// Stable key of one resolved snapshot face.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
pub struct FaceKey {
    pub occurrence: Vec<String>,
    pub body: String,
    pub name: String,
}

impl SnapshotFace {
    fn key(&self) -> FaceKey {
        FaceKey { occurrence: self.occurrence.clone(), body: self.body.clone(), name: self.name.clone() }
    }
}

impl TopologySnapshot {
    pub fn validate(&self) -> Result<()> {
        if self.faces.len() + self.bodies.len() > MAX_SNAPSHOT_ENTITIES {
            return Err(Error::Limit("geometry snapshot is too large".into()));
        }
        let mut keys = BTreeSet::new();
        for face in &self.faces {
            if !(face.area_m2.is_finite() && face.area_m2 >= 0.0) || face.name.is_empty() || face.lineage.is_empty() {
                return Err(Error::Invalid("snapshot faces need names and finite non-negative areas".into()));
            }
            if !keys.insert(face.key()) {
                return Err(Error::Invalid("snapshot face identities must be unique".into()));
            }
        }
        let mut bodies = BTreeSet::new();
        for body in &self.bodies {
            if !(body.volume_m3.is_finite() && body.volume_m3 >= 0.0) || !bodies.insert((&body.occurrence, &body.body)) {
                return Err(Error::Invalid("snapshot bodies need unique identities and finite volumes".into()));
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------- resolution and audit

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "state", rename_all = "snake_case")]
pub enum MemberState {
    Resolved {
        faces: Vec<FaceKey>,
    },
    ResolvedBody,
    Missing,
    /// The lineage now has several pieces. Review must name the intended pieces explicitly.
    Ambiguous {
        candidates: Vec<String>,
    },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct MemberReport {
    pub reference: EntityRef,
    pub resolution: MemberState,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SelectionReport {
    pub id: String,
    pub kind: EntityKind,
    pub resolved: bool,
    pub face_count: usize,
    pub body_count: usize,
    pub area_m2: f64,
    pub members: Vec<MemberReport>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    Error,
    Warning,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Issue {
    pub severity: Severity,
    pub code: String,
    pub message: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct BoundarySummary {
    pub id: String,
    pub kind: String,
    pub face_count: usize,
    pub area_m2: f64,
    /// Volume flow implied by the inlet data and resolved area, when it does not depend on normals.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub implied_volume_flow_m3_per_s: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct ContactSummary {
    pub id: String,
    pub area_a_m2: f64,
    pub area_b_m2: f64,
    /// Area-specific and total resistance using side A's area as the effective contact area.
    pub area_specific_m2_kelvin_per_watt: Option<f64>,
    pub total_kelvin_per_watt: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct FaceConflict {
    pub face: FaceKey,
    pub assignments: Vec<String>,
}

/// Coverage and readiness of a setup against the current exact geometry. `ready` means only
/// that preparation is complete and unambiguous; meshing and solver admission are separate.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SetupReport {
    pub geometry_revision: u64,
    pub geometry_changed: bool,
    pub selections: Vec<SelectionReport>,
    pub fluid_boundary_faces: usize,
    pub explicit_faces: usize,
    pub inherited_default_faces: usize,
    pub uncovered: Vec<FaceKey>,
    pub conflicts: Vec<FaceConflict>,
    pub boundaries: Vec<BoundarySummary>,
    pub contacts: Vec<ContactSummary>,
    pub issues: Vec<Issue>,
    /// No unresolved reference and no structural error: the setup may be applied.
    pub applicable: bool,
    /// Applicable, every fluid boundary face covered once, and a pressure reference present.
    pub ready: bool,
}

fn resolve(selection: &NamedSelection, snapshot: &TopologySnapshot) -> SelectionReport {
    let mut members = vec![];
    let mut faces = BTreeSet::new();
    let mut bodies = 0;
    let mut resolved = true;
    for reference in &selection.members {
        let owner = |occurrence: &Vec<String>, body: &String| *occurrence == reference.occurrence && *body == reference.body;
        let state = match &reference.face {
            None => {
                if snapshot.bodies.iter().any(|b| owner(&b.occurrence, &b.body)) {
                    bodies += 1;
                    MemberState::ResolvedBody
                } else {
                    MemberState::Missing
                }
            }
            Some(name) => {
                let exact: Vec<&SnapshotFace> = snapshot.faces.iter().filter(|f| owner(&f.occurrence, &f.body) && f.name == *name).collect();
                let pieces: Vec<&SnapshotFace> = snapshot.faces.iter().filter(|f| owner(&f.occurrence, &f.body) && f.lineage == *name).collect();
                match (exact.as_slice(), pieces.as_slice()) {
                    ([face], _) => MemberState::Resolved { faces: vec![face.key()] },
                    ([], [face]) => MemberState::Resolved { faces: vec![face.key()] },
                    ([], []) => MemberState::Missing,
                    (_, candidates) => MemberState::Ambiguous { candidates: candidates.iter().map(|f| f.name.clone()).collect() },
                }
            }
        };
        match &state {
            MemberState::Resolved { faces: found } => faces.extend(found.iter().cloned()),
            MemberState::ResolvedBody => {}
            MemberState::Missing | MemberState::Ambiguous { .. } => resolved = false,
        }
        members.push(MemberReport { reference: reference.clone(), resolution: state });
    }
    let area = snapshot.faces.iter().filter(|f| faces.contains(&f.key())).map(|f| f.area_m2).sum();
    SelectionReport { id: selection.id.clone(), kind: selection.kind, resolved, face_count: faces.len(), body_count: bodies, area_m2: area, members }
}

fn issue(issues: &mut Vec<Issue>, severity: Severity, code: &str, message: String) {
    issues.push(Issue { severity, code: code.into(), message });
}

/// Re-resolve every selection against the current geometry and audit boundary coverage,
/// conflicts, contacts and conduction. Structural errors are reported, not returned.
pub fn audit(setup: &Setup, snapshot: &TopologySnapshot) -> SetupReport {
    let mut issues = vec![];
    if let Err(error) = setup.validate().and_then(|()| snapshot.validate()) {
        issue(&mut issues, Severity::Error, "invalid", error.to_string());
    }
    let selections: Vec<SelectionReport> = setup.selections.iter().map(|s| resolve(s, snapshot)).collect();
    let by_id: BTreeMap<&str, &SelectionReport> = selections.iter().map(|r| (r.id.as_str(), r)).collect();
    for report in selections.iter().filter(|r| !r.resolved) {
        issue(
            &mut issues,
            Severity::Error,
            "unresolved_reference",
            format!("selection {} has missing or ambiguous members; review it before applying", report.id),
        );
    }
    let faces_of = |id: &str| -> Vec<FaceKey> {
        by_id
            .get(id)
            .map(|r| {
                r.members
                    .iter()
                    .filter_map(|m| match &m.resolution {
                        MemberState::Resolved { faces } => Some(faces.clone()),
                        _ => None,
                    })
                    .flatten()
                    .collect()
            })
            .unwrap_or_default()
    };
    let wetted: BTreeMap<FaceKey, f64> = snapshot.faces.iter().filter(|f| f.fluid_boundary).map(|f| (f.key(), f.area_m2)).collect();
    let mut owners: BTreeMap<FaceKey, Vec<String>> = BTreeMap::new();
    let mut boundaries = vec![];
    for assignment in &setup.boundaries {
        let faces: BTreeSet<FaceKey> = faces_of(&assignment.selection).into_iter().collect();
        let mut area = 0.0;
        for face in &faces {
            match wetted.get(face) {
                Some(a) => area += a,
                None => issue(
                    &mut issues,
                    Severity::Error,
                    "not_fluid_boundary",
                    format!("{} assigns face {} that does not bound the fluid region", assignment.id, face.name),
                ),
            }
            owners.entry(face.clone()).or_default().push(assignment.id.clone());
        }
        let implied = match &assignment.condition {
            BoundaryCondition::Inlet { inflow, .. } => inflow.volume_flow(area),
            _ => None,
        };
        if let BoundaryCondition::Inlet { inflow: InletSpec::Velocity { .. }, .. } = &assignment.condition {
            issue(
                &mut issues,
                Severity::Warning,
                "fixed_vector_inlet",
                format!("{} uses a fixed velocity vector; it does not follow a curved face's normal", assignment.id),
            );
        }
        boundaries.push(BoundarySummary {
            id: assignment.id.clone(),
            kind: assignment.condition.label().into(),
            face_count: faces.len(),
            area_m2: area,
            implied_volume_flow_m3_per_s: implied,
        });
    }
    let conflicts: Vec<FaceConflict> =
        owners.iter().filter(|(_, ids)| ids.len() > 1).map(|(face, ids)| FaceConflict { face: face.clone(), assignments: ids.clone() }).collect();
    for conflict in &conflicts {
        issue(
            &mut issues,
            Severity::Error,
            "conflicting_assignments",
            format!("face {} is assigned by {}; one assignment must own it", conflict.face.name, conflict.assignments.join(", ")),
        );
    }
    let mut uncovered = vec![];
    let mut inherited = 0;
    for face in wetted.keys().filter(|face| !owners.contains_key(*face)) {
        if setup.default_wall.is_some() {
            inherited += 1;
        } else {
            uncovered.push(face.clone());
        }
    }
    if inherited > 0 {
        issue(
            &mut issues,
            Severity::Warning,
            "inherited_default_wall",
            format!("{inherited} fluid boundary faces use the default wall; that is not an engineer-reviewed assignment"),
        );
    }
    if !uncovered.is_empty() {
        issue(&mut issues, Severity::Error, "uncovered_faces", format!("{} fluid boundary faces have no boundary condition", uncovered.len()));
    }
    if !wetted.is_empty() && !setup.boundaries.iter().any(|b| b.condition.sets_pressure()) {
        issue(&mut issues, Severity::Error, "no_pressure_reference", "incompressible flow needs an outlet or opening as pressure reference".into());
    }
    if !wetted.is_empty() && !setup.boundaries.iter().any(|b| b.condition.admits_inflow()) {
        issue(&mut issues, Severity::Warning, "no_inflow", "no inlet or opening admits flow into the region".into());
    }
    let area_of = |id: &str| by_id.get(id).map_or(0.0, |r| r.area_m2);
    let mut contacts = vec![];
    for contact in &setup.contacts {
        let (a, b) = (area_of(&contact.side_a), area_of(&contact.side_b));
        let side_a: BTreeSet<FaceKey> = faces_of(&contact.side_a).into_iter().collect();
        if faces_of(&contact.side_b).iter().any(|face| side_a.contains(face)) {
            issue(&mut issues, Severity::Error, "contact_sides_overlap", format!("{} uses the same face on both sides", contact.id));
        }
        if a > 0.0 && b > 0.0 && (a - b).abs() > 0.05 * a.max(b) {
            issue(
                &mut issues,
                Severity::Warning,
                "contact_area_mismatch",
                format!("{} sides differ in area ({a:.3e} vs {b:.3e} m²); total resistance uses side A", contact.id),
            );
        }
        contacts.push(ContactSummary {
            id: contact.id.clone(),
            area_a_m2: a,
            area_b_m2: b,
            area_specific_m2_kelvin_per_watt: contact.resistance.area_specific(a).ok(),
            total_kelvin_per_watt: contact.resistance.total(a).ok(),
        });
    }
    let mut conducted: BTreeMap<String, String> = BTreeMap::new();
    for conduction in &setup.conduction {
        for member in by_id.get(conduction.selection.as_str()).map(|r| r.members.as_slice()).unwrap_or_default() {
            let key = format!("{}/{}", member.reference.occurrence.join("/"), member.reference.body);
            if let Some(previous) = conducted.insert(key.clone(), conduction.id.clone()) {
                issue(
                    &mut issues,
                    Severity::Error,
                    "conflicting_conduction",
                    format!("body {key} has conductivity from both {previous} and {}", conduction.id),
                );
            }
        }
    }
    let geometry_changed = setup.geometry_revision != 0 && setup.geometry_revision != snapshot.geometry_revision;
    let applicable = !issues.iter().any(|i| i.severity == Severity::Error && matches!(i.code.as_str(), "invalid" | "unresolved_reference"));
    let ready = !issues.iter().any(|i| i.severity == Severity::Error);
    SetupReport {
        geometry_revision: snapshot.geometry_revision,
        geometry_changed,
        selections,
        fluid_boundary_faces: wetted.len(),
        explicit_faces: owners.keys().filter(|face| wetted.contains_key(*face)).count(),
        inherited_default_faces: inherited,
        uncovered,
        conflicts,
        boundaries,
        contacts,
        issues,
        applicable,
        ready,
    }
}
