//! Studies of the loaded model: a study whose geometry is the design itself, not a case copy.
//!
//! A [`ModelLink`] names the design's bodies by name and the faces that carry boundary
//! conditions by their persistent face names. The case is built from the design as it is when
//! the user builds or runs: the engine tessellates the bodies (the solver's input only) and
//! resolves every named face with exact geometry into a [`ResolvedFace`]; [`build`] then makes
//! the case. Nothing here evaluates geometry, so it stays portable.
//!
//! Phase 1 uses only what the companion solver takes today: bodies become voxelised surface
//! bodies or grid blocks, and a condition sits on one of the six faces of the domain. A model
//! face that lies on the domain boundary gives its role to that domain face; any other face is
//! refused with the reason, never moved to the nearest domain face.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::case_scene::{self, ImportOptions, PATCHES, Representation, RoleInput};
use crate::{Error, Result};

/// The finest mesh level, and the one that gives the automatic cell size.
pub const MAX_MESH_LEVEL: u8 = 7;
pub const AUTOMATIC_MESH_LEVEL: u8 = 3;

/// How much a mesh level scales the automatic cell size.
pub fn cell_scale(level: u8) -> f64 {
    2f64.powf((f64::from(AUTOMATIC_MESH_LEVEL) - f64::from(level)) / 2.0)
}

/// Most bodies a study may name.
pub const MAX_BODIES: usize = 256;
/// Most faces with a condition.
pub const MAX_OPENINGS: usize = 64;

/// Where the fluid is.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Analysis {
    /// Inside the model (an enclosure or duct): the domain is the bodies' box, and the bodies'
    /// walls bound the fluid. Domain faces without a condition are walls.
    #[default]
    Internal,
    /// Around the model: the domain is the bodies' box plus `margin_mm` on every side, and
    /// the domain faces keep the roles the case gives them (a wind-tunnel inlet and outlet).
    External,
}

/// A face of the design by its persistent name, so it follows the face through model edits.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FaceRef {
    pub body: String,
    pub face: String,
}

/// What a face with a condition does.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum OpeningRole {
    /// Uniform inflow at `speed_m_per_s` along the inward normal.
    Inlet {
        speed_m_per_s: f64,
    },
    /// Fixed kinematic pressure (m²/s², 0 for ambient).
    Outlet {
        pressure_m2_per_s2: f64,
    },
    Wall,
    Symmetry,
}

impl OpeningRole {
    fn role(self) -> RoleInput {
        match self {
            Self::Inlet { speed_m_per_s } => RoleInput::Inlet { speed: speed_m_per_s },
            Self::Outlet { pressure_m2_per_s2 } => RoleInput::Outlet { pressure: pressure_m2_per_s2 },
            Self::Wall => RoleInput::Wall,
            Self::Symmetry => RoleInput::Symmetry,
        }
    }
}

/// A named condition on a face of the model.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Opening {
    pub name: String,
    pub face: FaceRef,
    pub role: OpeningRole,
}

/// What a local refinement splits the cells around.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum RefineTarget {
    Body { body: String },
    Face { face: FaceRef },
    /// A box of the design (mm), such as one solution-adaptive refinement gives.
    Region { low_mm: [f64; 3], high_mm: [f64; 3] },
}

/// Smaller cells around a body or a face of the model: every cell whose span overlaps the
/// target's box, plus `margin_cells` cells, is split into `factor` along each axis.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Refinement {
    pub name: String,
    pub target: RefineTarget,
    pub factor: usize,
    #[serde(default)]
    pub margin_cells: f64,
}

/// What a goal of a model study monitors.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ModelGoalKind {
    /// Volume flow out through a face of the model, m³/s.
    FaceVolumeFlow,
    /// Area-weighted static pressure on a face of the model, Pa.
    FaceMeanPressure,
    /// Static pressure at a point of the fluid, Pa.
    PointPressure,
    /// Arithmetic of other goals (`+ - * /`, parentheses, numbers and goal names).
    Equation,
    /// Highest temperature of a body, K (heat transfer).
    BodyMaxTemperature,
    /// Volume-mean temperature of a body, K (heat transfer).
    BodyMeanTemperature,
    /// Temperature at a point, K (heat transfer).
    PointTemperature,
    /// Over the whole fluid: largest and mean speed, m/s, and largest, smallest and mean
    /// static pressure, Pa (companion 0.16.0, `goals_v3`).
    GlobalMaxSpeed,
    GlobalMeanSpeed,
    GlobalMaxPressure,
    GlobalMinPressure,
    GlobalMeanPressure,
    /// Over the whole fluid or every solid: largest and mean fluid temperature and largest
    /// solid temperature, K (heat transfer).
    GlobalMaxFluidTemperature,
    GlobalMeanFluidTemperature,
    GlobalMaxSolidTemperature,
}

impl ModelGoalKind {
    /// Whether the goal reads temperatures.
    pub fn thermal(self) -> bool {
        matches!(
            self,
            Self::BodyMaxTemperature
                | Self::BodyMeanTemperature
                | Self::PointTemperature
                | Self::GlobalMaxFluidTemperature
                | Self::GlobalMeanFluidTemperature
                | Self::GlobalMaxSolidTemperature
        )
    }
}

/// A monitored quantity of a model study.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelGoal {
    /// Letters, digits and underscores, starting with a letter (equations name goals).
    pub name: String,
    pub kind: ModelGoalKind,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub face: Option<FaceRef>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub point_mm: Option<[f64; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expression: Option<String>,
    /// The body of a body goal.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub body: Option<String>,
}

/// A fan blowing into the domain through a face of the model on the domain boundary, following
/// its pressure-flow curve against an outlet (the companion's external fan, `fan_curve_v1`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelFan {
    pub name: String,
    pub face: FaceRef,
    /// `[volume flow m³/s, pressure rise Pa]`, flow strictly increasing, pressure strictly
    /// decreasing, 2 to 64 points.
    pub curve: Vec<[f64; 2]>,
    /// The curve's pressure is total (static plus dynamic at the fan face), not static.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub total_pressure: bool,
}

impl ModelFan {
    fn validate(&self) -> Result<()> {
        text(&self.name, "a fan name")?;
        text(&self.face.body, "a body name")?;
        text(&self.face.face, "a face name")?;
        let ok = (2..=64).contains(&self.curve.len())
            && self.curve.iter().all(|[q, p]| q.is_finite() && p.is_finite() && *q >= 0.0 && *p >= 0.0)
            && self.curve.windows(2).all(|w| w[1][0] > w[0][0] && w[1][1] < w[0][1]);
        if !ok {
            return Err(bad(format!(
                "{}: a fan curve has 2 to 64 points of flow (m³/s) and pressure rise (Pa), at least 0, with flow rising and pressure falling",
                self.name
            )));
        }
        Ok(())
    }
}

/// What a thermal condition on a face on the domain boundary does.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ThermalFace {
    Adiabatic,
    Temperature { temperature_k: f64 },
    HeatFlux { w_per_m2: f64 },
    Convection { h_w_per_m2_k: f64, ambient_k: f64 },
}

/// The material and the power of a body of the model.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BodyHeat {
    pub body: String,
    /// Isotropic thermal conductivity, W/(m·K).
    pub conductivity_w_per_m_k: f64,
    /// Heat released evenly in the body, W.
    #[serde(default)]
    pub power_w: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub density_kg_per_m3: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specific_heat_j_per_kg_k: Option<f64>,
    /// The temperature the body must stay below, K (reported as a margin).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_temperature_k: Option<f64>,
}

/// Heat released on a face of the model, W (companion 0.15.0, `face_heat_sources_v1`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FaceHeat {
    pub name: String,
    pub face: FaceRef,
    pub power_w: f64,
}

/// A thermal condition on a face of the model that lies on the domain boundary.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FaceThermal {
    pub name: String,
    pub face: FaceRef,
    pub condition: ThermalFace,
}

/// The fluid's thermal properties; air at room temperature by default.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FluidHeat {
    pub specific_heat_j_per_kg_k: f64,
    pub conductivity_w_per_m_k: f64,
}

impl Default for FluidHeat {
    fn default() -> Self {
        Self { specific_heat_j_per_kg_k: 1006.0, conductivity_w_per_m_k: 0.0257 }
    }
}

fn room_temperature() -> f64 {
    293.15
}

/// Heat transfer in a study of the model: the bodies conduct heat (conjugate heat transfer),
/// with their materials and power, heat sources on faces, and thermal conditions on the
/// faces on the domain boundary. Bodies without a material get the starter one (1 W/(m·K)).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelHeat {
    #[serde(default)]
    pub fluid: FluidHeat,
    /// The temperature of the fluid flowing in, K.
    #[serde(default = "room_temperature")]
    pub inflow_temperature_k: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub materials: Vec<BodyHeat>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub face_sources: Vec<FaceHeat>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub face_conditions: Vec<FaceThermal>,
}

impl Default for ModelHeat {
    fn default() -> Self {
        Self {
            fluid: FluidHeat::default(),
            inflow_temperature_k: room_temperature(),
            materials: Vec::new(),
            face_sources: Vec::new(),
            face_conditions: Vec::new(),
        }
    }
}

/// The starter solid material: what a body without a material conducts, W/(m·K).
pub const STARTER_CONDUCTIVITY: f64 = 1.0;

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl ModelHeat {
    fn validate(&self) -> Result<()> {
        if !(positive(self.fluid.specific_heat_j_per_kg_k) && positive(self.fluid.conductivity_w_per_m_k)) {
            return Err(bad("the fluid's specific heat and conductivity must be positive"));
        }
        if !(positive(self.inflow_temperature_k) && self.inflow_temperature_k < 1e4) {
            return Err(bad("the inflow temperature is in kelvin, above 0 and below 10 000"));
        }
        if self.materials.len() > MAX_BODIES || self.face_sources.len() > MAX_OPENINGS || self.face_conditions.len() > MAX_OPENINGS {
            return Err(Error::Limit(format!("a model study has at most {MAX_BODIES} materials and {MAX_OPENINGS} heat sources and thermal conditions")));
        }
        let mut bodies = std::collections::BTreeSet::new();
        for m in &self.materials {
            text(&m.body, "a body name")?;
            if !bodies.insert(m.body.as_str()) {
                return Err(bad(format!("{} has two materials", m.body)));
            }
            let optional = |v: Option<f64>| v.is_none_or(positive);
            if !(positive(m.conductivity_w_per_m_k) && m.conductivity_w_per_m_k < 1e5)
                || !(m.power_w.is_finite() && m.power_w >= 0.0)
                || !optional(m.density_kg_per_m3)
                || !optional(m.specific_heat_j_per_kg_k)
                || !optional(m.max_temperature_k)
            {
                return Err(bad(format!(
                    "{}: conductivity must be positive (below 1e5), power at least 0, density, specific heat and the limit positive",
                    m.body
                )));
            }
        }
        let mut names = std::collections::BTreeSet::new();
        for source in &self.face_sources {
            text(&source.name, "a heat source name")?;
            text(&source.face.body, "a body name")?;
            text(&source.face.face, "a face name")?;
            if !names.insert(source.name.as_str()) {
                return Err(bad(format!("two heat sources or thermal conditions are named {}", source.name)));
            }
            if !(source.power_w.is_finite() && source.power_w.abs() <= 1e9) {
                return Err(bad(format!("{}: the power must be finite", source.name)));
            }
        }
        for condition in &self.face_conditions {
            text(&condition.name, "a thermal condition name")?;
            text(&condition.face.body, "a body name")?;
            text(&condition.face.face, "a face name")?;
            if !names.insert(condition.name.as_str()) {
                return Err(bad(format!("two heat sources or thermal conditions are named {}", condition.name)));
            }
            let ok = match condition.condition {
                ThermalFace::Adiabatic => true,
                ThermalFace::Temperature { temperature_k } => positive(temperature_k),
                ThermalFace::HeatFlux { w_per_m2 } => w_per_m2.is_finite(),
                ThermalFace::Convection { h_w_per_m2_k, ambient_k } => positive(h_w_per_m2_k) && positive(ambient_k),
            };
            if !ok {
                return Err(bad(format!("{}: temperatures are positive kelvin, the flux finite and the coefficient positive", condition.name)));
            }
        }
        Ok(())
    }
}

/// Most refinements and goals a model study holds.
pub const MAX_REFINEMENTS: usize = 16;
pub const MAX_GOALS: usize = 32;

/// What a model study takes from the design.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelLink {
    #[serde(default)]
    pub analysis: Analysis,
    /// Bodies of the design in the study; empty means every body.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bodies: Vec<String>,
    /// A part selection of the project (`parts`) whose bodies the study uses, worked out again
    /// at every build; with it, `bodies` is empty.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub selection: Option<String>,
    /// Space around the bodies on every side, for an external analysis.
    #[serde(default)]
    pub margin_mm: f64,
    /// Cell size; by default chosen from the bodies (see `case_scene::import_components`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_mm: Option<f64>,
    /// Without `cell_mm`: how fine the automatic grid is, 1 (coarsest) to 7; 3 is the automatic
    /// size, and each level makes the cells √2 times smaller (`cell_scale`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mesh_level: Option<u8>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub openings: Vec<Opening>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub refinements: Vec<Refinement>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub goals: Vec<ModelGoal>,
    /// Keep the fluid part of the cells the bodies cut, instead of a staircase of whole cells
    /// (companion 0.14.0, `cut_cells_v1`).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub cut_cells: bool,
    /// Heat transfer; none for a flow-only study.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat: Option<ModelHeat>,
    /// A fan on a face of the model on the domain boundary.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fan: Option<ModelFan>,
    /// The digest of the model input the case was last built from. Another digest now means
    /// the design changed since, and the case must be built again before a run.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub built_from: Option<String>,
}

fn text(value: &str, what: &str) -> Result<()> {
    if value.trim().is_empty() || value.len() > 200 || value.chars().any(char::is_control) {
        return Err(Error::Invalid(format!("{what} must be 1 to 200 characters without control characters")));
    }
    Ok(())
}

impl ModelLink {
    /// Every face the link names, in one order: the conditions', the refinements', the goals',
    /// the heat sources', the thermal conditions' and the fan's. [`ModelInput::faces`] resolves them in this order.
    pub fn face_refs(&self) -> Vec<&FaceRef> {
        let openings = self.openings.iter().map(|o| &o.face);
        let refined = self.refinements.iter().filter_map(|r| match &r.target {
            RefineTarget::Face { face } => Some(face),
            RefineTarget::Body { .. } | RefineTarget::Region { .. } => None,
        });
        let heat = self.heat.iter().flat_map(|h| h.face_sources.iter().map(|f| &f.face).chain(h.face_conditions.iter().map(|f| &f.face)));
        let fan = self.fan.iter().map(|f| &f.face);
        openings.chain(refined).chain(self.goals.iter().filter_map(|g| g.face.as_ref())).chain(heat).chain(fan).collect()
    }

    /// Bounded structure only; whether the bodies and faces exist is decided when building.
    pub fn validate(&self) -> Result<()> {
        if self.bodies.len() > MAX_BODIES {
            return Err(Error::Limit(format!("a model study names at most {MAX_BODIES} bodies")));
        }
        for body in &self.bodies {
            text(body, "a body name")?;
        }
        if let Some(selection) = &self.selection {
            text(selection, "a selection name")?;
            if !self.bodies.is_empty() {
                return Err(Error::Invalid("a model study takes its bodies from a list or from a part selection, not both".into()));
            }
        }
        if !(self.margin_mm.is_finite() && (0.0..=1e6).contains(&self.margin_mm)) {
            return Err(Error::Invalid("the margin must be 0 to 1e6 mm".into()));
        }
        if self.mesh_level.is_some_and(|level| !(1..=MAX_MESH_LEVEL).contains(&level)) {
            return Err(Error::Invalid(format!("the mesh level is 1 to {MAX_MESH_LEVEL}")));
        }
        if self.cell_mm.is_some_and(|cell| !(cell.is_finite() && cell > 0.0 && cell <= 1e6)) {
            return Err(Error::Invalid("the cell size must be positive and at most 1e6 mm".into()));
        }
        if self.openings.len() > MAX_OPENINGS {
            return Err(Error::Limit(format!("a model study has at most {MAX_OPENINGS} face conditions")));
        }
        let mut names = std::collections::BTreeSet::new();
        for opening in &self.openings {
            text(&opening.name, "a condition name")?;
            text(&opening.face.body, "a body name")?;
            text(&opening.face.face, "a face name")?;
            if !names.insert(opening.name.as_str()) {
                return Err(Error::Invalid(format!("two face conditions are named {}", opening.name)));
            }
            match opening.role {
                OpeningRole::Inlet { speed_m_per_s } if !(speed_m_per_s.is_finite() && speed_m_per_s > 0.0 && speed_m_per_s < 1e3) => {
                    return Err(Error::Invalid("inlet speed must be positive and below 1000 m/s".into()));
                }
                OpeningRole::Outlet { pressure_m2_per_s2 } if !(pressure_m2_per_s2.is_finite() && pressure_m2_per_s2.abs() <= 1e6) => {
                    return Err(Error::Invalid("outlet pressure must be finite and at most 1e6".into()));
                }
                _ => {}
            }
        }
        if self.refinements.len() > MAX_REFINEMENTS {
            return Err(Error::Limit(format!("a model study has at most {MAX_REFINEMENTS} local refinements")));
        }
        let mut refined = std::collections::BTreeSet::new();
        for r in &self.refinements {
            text(&r.name, "a refinement name")?;
            if !refined.insert(r.name.as_str()) {
                return Err(Error::Invalid(format!("two refinements are named {}", r.name)));
            }
            if !(2..=case_scene::MAX_REFINEMENT).contains(&r.factor) {
                return Err(Error::Invalid(format!("{}: split cells by 2 to {}", r.name, case_scene::MAX_REFINEMENT)));
            }
            if !(r.margin_cells.is_finite() && (0.0..=64.0).contains(&r.margin_cells)) {
                return Err(Error::Invalid(format!("{}: the margin must be 0 to 64 cells", r.name)));
            }
            match &r.target {
                RefineTarget::Body { body } => text(body, "a body name")?,
                RefineTarget::Face { face } => {
                    text(&face.body, "a body name")?;
                    text(&face.face, "a face name")?;
                }
                RefineTarget::Region { low_mm, high_mm } => {
                    if !(0..3).all(|a| low_mm[a].is_finite() && high_mm[a].is_finite() && low_mm[a] < high_mm[a]) {
                        return Err(bad(format!("{}: the region needs a finite box with low below high", r.name)));
                    }
                }
            }
        }
        if self.goals.len() > MAX_GOALS {
            return Err(Error::Limit(format!("a model study has at most {MAX_GOALS} goals")));
        }
        let mut goals = std::collections::BTreeSet::new();
        for g in &self.goals {
            let identifier = g.name.len() <= 60
                && g.name.chars().next().is_some_and(|c| c.is_ascii_alphabetic())
                && g.name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_');
            if !identifier {
                return Err(Error::Invalid(format!("goal name {:?}: use letters, digits and underscores, starting with a letter", g.name)));
            }
            if !goals.insert(g.name.as_str()) {
                return Err(Error::Invalid(format!("two goals are named {}", g.name)));
            }
            let ok = match g.kind {
                ModelGoalKind::FaceVolumeFlow | ModelGoalKind::FaceMeanPressure => g.face.is_some(),
                ModelGoalKind::PointPressure => g.point_mm.is_some_and(|p| p.iter().all(|v| v.is_finite())),
                ModelGoalKind::Equation => g.expression.as_ref().is_some_and(|e| !e.trim().is_empty() && e.len() <= 200),
                ModelGoalKind::BodyMaxTemperature | ModelGoalKind::BodyMeanTemperature => {
                    g.body.as_ref().is_some_and(|b| text(b, "a body name").is_ok())
                }
                ModelGoalKind::PointTemperature => g.point_mm.is_some_and(|p| p.iter().all(|v| v.is_finite())),
                ModelGoalKind::GlobalMaxSpeed
                | ModelGoalKind::GlobalMeanSpeed
                | ModelGoalKind::GlobalMaxPressure
                | ModelGoalKind::GlobalMinPressure
                | ModelGoalKind::GlobalMeanPressure
                | ModelGoalKind::GlobalMaxFluidTemperature
                | ModelGoalKind::GlobalMeanFluidTemperature
                | ModelGoalKind::GlobalMaxSolidTemperature => true,
            };
            if g.kind.thermal() && self.heat.is_none() {
                return Err(Error::Invalid(format!("goal {}: temperature goals need heat transfer on", g.name)));
            }
            if !ok {
                return Err(Error::Invalid(format!(
                    "goal {}: a face goal needs a face, a body goal a body, a point goal a finite point, an equation its expression",
                    g.name
                )));
            }
        }
        if let Some(heat) = &self.heat {
            heat.validate()?;
        }
        if let Some(fan) = &self.fan {
            fan.validate()?;
            if self.cut_cells {
                return Err(bad("cut cells do not work with a fan yet; turn cut cells off"));
            }
        }
        if self.built_from.as_ref().is_some_and(|hash| hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit())) {
            return Err(Error::Invalid("the model digest must be 64 hexadecimal digits".into()));
        }
        Ok(())
    }
}

/// A named face of the model as exact geometry resolved it now.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ResolvedFace {
    /// The condition's name.
    pub name: String,
    /// The exact plane of a planar face (outward unit normal, a point on it, mm); `None` for a
    /// curved face.
    pub plane: Option<([f64; 3], [f64; 3])>,
    /// The design body the face is on.
    pub body: String,
    /// Which of that body's triangles in [`ModelInput::obj`] lie on the face (by exact nearest
    /// face of each triangle's centroid), in the body's own triangle order.
    pub triangles: Vec<usize>,
}

/// The centroid of each triangle of body `body` in the solver's tessellation `obj`, in its
/// triangle order: the points whose exact nearest model face says which face each triangle
/// belongs to (the engine answers that with the B-rep, never with this mesh).
pub fn triangle_centroids(obj: &str, body: &str) -> Result<Vec<[f64; 3]>> {
    let objects = case_scene::parse_obj(obj)?;
    let object = objects.iter().find(|o| o.name == body).ok_or_else(|| bad(format!("{body} was not tessellated")))?;
    Ok(object
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| object.vertices[i]);
            [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0)
        })
        .collect())
}

/// What the engine reads from the design for one build.
#[derive(Clone, Debug)]
pub struct ModelInput {
    /// The bodies, as the OBJ export writes them (mm).
    pub obj: String,
    /// The bodies' exact bounding box (mm).
    pub bounds_mm: [[f64; 3]; 2],
    /// Their exact total volume (mm³), to tell how much of it the grid holds.
    pub volume_mm3: f64,
    /// Every face the link names, in the order of [`ModelLink::face_refs`].
    pub faces: Vec<ResolvedFace>,
}

impl ModelInput {
    /// The digest of everything the case is built from: the link's settings, the bodies'
    /// triangles and the faces' exact planes. The stored `built_from` is left out.
    pub fn digest(&self, link: &ModelLink) -> Result<String> {
        let mut settings = link.clone();
        settings.built_from = None;
        let mut hasher = Sha256::new();
        hasher.update(serde_json::to_vec(&settings).map_err(|e| Error::Json(e.to_string()))?);
        hasher.update(self.obj.as_bytes());
        for face in &self.faces {
            hasher.update(face.name.as_bytes());
            if let Some((normal, point)) = face.plane {
                for v in normal.iter().chain(&point) {
                    hasher.update(v.to_bits().to_le_bytes());
                }
            }
        }
        Ok(format!("{:x}", hasher.finalize()))
    }
}

/// A case built from the model, and what the user should know about it.
#[derive(Clone, Debug, Serialize, JsonSchema)]
pub struct Built {
    pub case: Value,
    /// The input digest to store in [`ModelLink::built_from`].
    pub digest: String,
    /// The CAD point at the domain origin (mm).
    pub origin_mm: [f64; 3],
    /// Which domain face each condition landed on, in order.
    pub placed: Vec<(String, String)>,
    /// Fraction of the bodies' exact volume that the grid's solid cells hold.
    pub solid_fraction: f64,
    /// Separate fluid regions on the grid.
    pub fluid_regions: usize,
    pub notes: Vec<String>,
    /// Admission problems of the built case; empty when it would be accepted.
    pub issues: Vec<String>,
}

fn bad(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

/// The domain face a planar model face lies on: its normal along an axis, pointing out of the
/// domain, and its plane within half a cell of that domain face.
fn patch_of(plane: ([f64; 3], [f64; 3]), domain: [[f64; 3]; 2], cell_mm: [f64; 3]) -> Option<usize> {
    let (normal, point) = plane;
    let axis = (0..3).find(|&a| normal[a].abs() > 1.0 - 1e-9)?;
    let side = usize::from(normal[axis] > 0.0);
    let boundary = domain[side][axis];
    ((point[axis] - boundary).abs() <= 0.5 * cell_mm[axis]).then_some(2 * axis + side)
}

/// The sub-cells per cell the companion's cut cells measure on (`cut_cells::MAX_SUB` there).
const CUT_SUB: usize = 8;

/// The bodies' volume as the companion's cut cells hold it: sub-cell centres inside a body on the
/// lattice `CUT_SUB` times finer than the grid, times their volume (m³).
fn cut_solid_volume(case: &Value) -> Result<f64> {
    let nodes = case_scene::nodes(case)?;
    let fine: [Vec<f64>; 3] = std::array::from_fn(|axis| {
        nodes[axis].windows(2).flat_map(|w| (0..CUT_SUB).map(move |s| w[0] + (w[1] - w[0]) * (s as f64 + 0.5) / CUT_SUB as f64)).collect()
    });
    let width: [Vec<f64>; 3] =
        std::array::from_fn(|axis| nodes[axis].windows(2).flat_map(|w| std::iter::repeat_n((w[1] - w[0]) / CUT_SUB as f64, CUT_SUB)).collect());
    let mut solid = vec![false; fine[0].len() * fine[1].len() * fine[2].len()];
    for b in case.get("bodies").and_then(Value::as_array).into_iter().flatten() {
        let vertices: Vec<[f64; 3]> = serde_json::from_value(b["vertices"].clone()).unwrap_or_default();
        let triangles: Vec<[usize; 3]> = serde_json::from_value(b["triangles"].clone()).unwrap_or_default();
        for (s, i) in solid.iter_mut().zip(crate::fluid_volume::voxelize_body(&vertices, &triangles, &fine[0], &fine[1], &fine[2])) {
            *s |= i;
        }
    }
    let (nx, ny) = (fine[0].len(), fine[1].len());
    Ok(solid.iter().enumerate().filter(|(_, s)| **s).map(|(c, _)| width[0][c % nx] * width[1][(c / nx) % ny] * width[2][c / (nx * ny)]).sum())
}

/// For each domain face, whether fluid touches it in the thin layer the companion's cut cells
/// measure next to it (an eighth of a cell): sample points there, against every surface body.
fn wetted_layers(case: &Value) -> Result<[bool; 6]> {
    let nodes = case_scene::nodes(case)?;
    let fine = |axis: usize| -> Vec<f64> {
        nodes[axis].windows(2).flat_map(|w| (0..CUT_SUB).map(move |s| w[0] + (w[1] - w[0]) * (s as f64 + 0.5) / CUT_SUB as f64)).collect()
    };
    let bodies: Vec<(Vec<[f64; 3]>, Vec<[usize; 3]>)> = case
        .get("bodies")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .map(|b| {
            (serde_json::from_value(b["vertices"].clone()).unwrap_or_default(), serde_json::from_value(b["triangles"].clone()).unwrap_or_default())
        })
        .collect();
    let mut wetted = [false; 6];
    for (patch, wet) in wetted.iter_mut().enumerate() {
        let axis = patch / 2;
        let list = &nodes[axis];
        let width = if patch % 2 == 0 { list[1] - list[0] } else { list[list.len() - 1] - list[list.len() - 2] };
        let layer = if patch % 2 == 0 { list[0] + 0.5 * width / CUT_SUB as f64 } else { list[list.len() - 1] - 0.5 * width / CUT_SUB as f64 };
        let mut axes: [Vec<f64>; 3] = [fine(0), fine(1), fine(2)];
        axes[axis] = vec![layer];
        let mut solid = vec![false; axes[0].len() * axes[1].len() * axes[2].len()];
        for (vertices, triangles) in &bodies {
            let inside = crate::fluid_volume::voxelize_body(vertices, triangles, &axes[0], &axes[1], &axes[2]);
            for (s, i) in solid.iter_mut().zip(inside) {
                *s |= i;
            }
        }
        *wet = solid.iter().any(|s| !s);
    }
    Ok(wetted)
}

/// The box of a face's triangles in case coordinates (metres), from the body as the case holds it.
fn face_box(case: &Value, face: &ResolvedFace) -> Option<([f64; 3], [f64; 3])> {
    let name = case_scene::part_name(&face.body);
    let body = case.get("bodies")?.as_array()?.iter().find(|b| b["name"].as_str() == Some(name.as_str()))?;
    let vertices: Vec<[f64; 3]> = serde_json::from_value(body["vertices"].clone()).ok()?;
    let triangles: Vec<[usize; 3]> = serde_json::from_value(body["triangles"].clone()).ok()?;
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for &t in &face.triangles {
        for &v in triangles.get(t)? {
            let p = vertices.get(v)?;
            for a in 0..3 {
                low[a] = low[a].min(p[a]);
                high[a] = high[a].max(p[a]);
            }
        }
    }
    low[0].is_finite().then_some((low, high))
}

/// A case identifier for a condition's patch on a body: its name made ASCII and unique.
fn patch_name(name: &str, taken: &mut std::collections::BTreeSet<String>) -> String {
    let mut base: String = name.trim().chars().map(|c| if c.is_ascii_alphanumeric() || c == '-' { c } else { '_' }).take(70).collect();
    if !base.as_bytes().first().is_some_and(u8::is_ascii_alphanumeric) {
        base.insert_str(0, "face_");
    }
    let mut candidate = base.clone();
    let mut n = 2;
    while taken.contains(&candidate) {
        candidate = format!("{base}_{n}");
        n += 1;
    }
    taken.insert(candidate.clone());
    candidate
}

/// Put `opening` on its own body's surface: a solver patch over the body's triangles on the
/// face (`body_face_conditions_v1`), with the condition's velocity and pressure. Returns the patch.
fn body_face_condition(
    case: &mut Value,
    opening: &Opening,
    face: &ResolvedFace,
    objects: &[case_scene::MeshObject],
    taken: &mut std::collections::BTreeSet<String>,
) -> Result<String> {
    if case.get("thermal").is_some_and(|t| !t.is_null()) && matches!(opening.role, OpeningRole::Inlet { .. } | OpeningRole::Outlet { .. }) {
        return Err(bad(format!(
            "{}: inlets and outlets on faces inside the domain do not work with heat transfer yet; put it on a face on the domain boundary",
            opening.name
        )));
    }
    let object = objects
        .iter()
        .find(|o| o.name == face.body)
        .ok_or_else(|| bad(format!("{}: the body {} is not in the model's tessellation", opening.name, face.body)))?;
    if face.triangles.is_empty() || face.triangles.iter().any(|t| *t >= object.triangles.len()) {
        return Err(bad(format!("{}: the face has no triangles in the model's tessellation", opening.name)));
    }
    // The face's outward normal, area-weighted over its triangles (exact for a flat face).
    let mut normal = [0.0; 3];
    for &t in &face.triangles {
        let [a, b, c] = object.triangles[t].map(|i| object.vertices[i]);
        let (u, v) = ([0, 1, 2].map(|k| b[k] - a[k]), [0, 1, 2].map(|k| c[k] - a[k]));
        let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        for k in 0..3 {
            normal[k] += n[k];
        }
    }
    let length = normal.iter().map(|v| v * v).sum::<f64>().sqrt();
    if !(length.is_finite() && length > 0.0) {
        return Err(bad(format!("{}: the face's triangles enclose no area", opening.name)));
    }
    let normal = normal.map(|v| v / length);
    let patch = patch_name(&opening.name, taken);
    let push = |case: &mut Value, key: &str, row: Value| -> Result<()> {
        let object = case.as_object_mut().ok_or_else(|| bad("case must be an object"))?;
        object.entry(key).or_insert_with(|| serde_json::json!([])).as_array_mut().ok_or_else(|| bad(format!("{key} must be an array")))?.push(row);
        Ok(())
    };
    push(case, "face_conditions", serde_json::json!({"name": patch, "body": case_scene::part_name(&face.body), "triangles": face.triangles}))?;
    // The fluid is on the outward side of the body's face, so inflow runs along the normal.
    match opening.role {
        OpeningRole::Inlet { speed_m_per_s } => {
            let value = serde_json::json!({"x": normal[0] * speed_m_per_s, "y": normal[1] * speed_m_per_s, "z": normal[2] * speed_m_per_s});
            push(case, "velocity_boundaries", serde_json::json!({"patch": patch, "kind": "fixed_value", "value": value}))?;
            push(case, "pressure_boundaries", serde_json::json!({"patch": patch, "kind": "zero_gradient", "value": 0.0}))?;
        }
        OpeningRole::Outlet { pressure_m2_per_s2 } => {
            push(case, "velocity_boundaries", serde_json::json!({"patch": patch, "kind": "zero_gradient", "value": {"x": 0.0, "y": 0.0, "z": 0.0}}))?;
            push(case, "pressure_boundaries", serde_json::json!({"patch": patch, "kind": "fixed_value", "value": pressure_m2_per_s2}))?;
        }
        OpeningRole::Wall => {
            push(case, "velocity_boundaries", serde_json::json!({"patch": patch, "kind": "fixed_value", "value": {"x": 0.0, "y": 0.0, "z": 0.0}}))?;
        }
        OpeningRole::Symmetry => {
            push(case, "velocity_boundaries", serde_json::json!({"patch": patch, "kind": "slip", "value": {"x": 0.0, "y": 0.0, "z": 0.0}}))?;
        }
    }
    Ok(patch)
}

/// Put `fan` on the domain face its model face lies on, blowing against a domain face with a
/// fixed pressure (an outlet).
fn put_fan(
    case: &mut Value,
    fan: &ModelFan,
    face: &ResolvedFace,
    domain: [[f64; 3]; 2],
    cell_mm: [f64; 3],
    wetted: [bool; 6],
    given: &[Option<(String, OpeningRole)>; 6],
) -> Result<()> {
    let Some(side) = face.plane.and_then(|plane| patch_of(plane, domain, cell_mm)).filter(|p| wetted[*p]) else {
        return Err(bad(format!("{}: a fan goes on a flat face of the model on the domain boundary, such as a duct's open end", fan.name)));
    };
    if let Some((other, _)) = &given[side] {
        return Err(bad(format!("{} and {other} are on the same domain face ({}); a fan takes the whole face", fan.name, PATCHES[side])));
    }
    let patch = PATCHES[side];
    let reference = case
        .get("pressure_boundaries")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter(|row| row["kind"] == "fixed_value")
        .filter_map(|row| row["patch"].as_str())
        .find(|p| PATCHES.contains(p) && *p != patch)
        .map(str::to_string)
        .ok_or_else(|| bad(format!("{}: the fan needs an outlet on the domain boundary to blow against", fan.name)))?;
    if let Some(object) = case.as_object_mut() {
        // The fan owns its face's velocity and pressure.
        for key in ["velocity_boundaries", "pressure_boundaries"] {
            if let Some(rows) = object.get_mut(key).and_then(Value::as_array_mut) {
                rows.retain(|row| row["patch"].as_str() != Some(patch));
            }
        }
        let curve: Vec<Value> = fan.curve.iter().map(|[q, p]| serde_json::json!({"volume_flow_m3_per_s": q, "pressure_rise_pa": p})).collect();
        let middle = fan.curve.first().zip(fan.curve.last()).map_or(0.0, |(a, b)| 0.5 * (a[0] + b[0]));
        object.insert(
            "fan".into(),
            serde_json::json!({
                "patch": patch,
                "reference_patch": reference,
                "pressure_convention": if fan.total_pressure { "total" } else { "static" },
                "curve": curve,
                "initial_volume_flow_m3_per_s": middle,
            }),
        );
    }
    Ok(())
}

/// The case's thermal block for `heat`: the fluid, the inflow temperature and a solid for
/// every body of the case (its material, or the starter one with a note).
fn thermal_block(case: &Value, heat: &ModelHeat, notes: &mut Vec<String>) -> Result<Value> {
    let parts: Vec<String> = case
        .get("bodies")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .chain(case.get("obstacles").and_then(Value::as_array).into_iter().flatten())
        .filter_map(|b| b["name"].as_str().map(str::to_string))
        .collect();
    // A part the study leaves out keeps its material for when it is back: a live part selection
    // or a size limit takes parts in and out of a study as the design changes.
    let left_out: Vec<&str> = heat.materials.iter().filter(|m| !parts.contains(&case_scene::part_name(&m.body))).map(|m| m.body.as_str()).collect();
    if !left_out.is_empty() {
        notes.push(format!("{} have a material but the study leaves them out; it is kept for when they are back", left_out.join(", ")));
    }
    let mut starter = Vec::new();
    let solids: Vec<Value> = parts
        .iter()
        .map(|part| match heat.materials.iter().find(|m| case_scene::part_name(&m.body) == *part) {
            Some(m) => {
                let mut row = serde_json::json!({
                    "obstacle": part,
                    "conductivity_w_per_m_k": vec![m.conductivity_w_per_m_k; 3],
                    "power_w": m.power_w,
                });
                for (key, value) in [
                    ("density_kg_per_m3", m.density_kg_per_m3),
                    ("specific_heat_j_per_kg_k", m.specific_heat_j_per_kg_k),
                    ("max_temperature_k", m.max_temperature_k),
                ] {
                    if let Some(value) = value {
                        row[key] = serde_json::json!(value);
                    }
                }
                row
            }
            None => {
                starter.push(part.clone());
                serde_json::json!({"obstacle": part, "conductivity_w_per_m_k": vec![STARTER_CONDUCTIVITY; 3], "power_w": 0.0})
            }
        })
        .collect();
    if !starter.is_empty() {
        notes.push(format!("{} have no material and conduct like the starter material ({STARTER_CONDUCTIVITY} W/(m·K), no power); set their materials", starter.join(", ")));
    }
    Ok(serde_json::json!({
        "fluid": {"specific_heat_j_per_kg_k": heat.fluid.specific_heat_j_per_kg_k, "conductivity_w_per_m_k": heat.fluid.conductivity_w_per_m_k},
        "inflow_temperature_k": heat.inflow_temperature_k,
        "solids": solids,
    }))
}

/// Heat sources on body faces (`face_heat_sources_v1`) and thermal conditions on the faces on
/// the domain boundary, into the case's thermal block.
fn heat_on_faces(
    case: &mut Value,
    heat: &ModelHeat,
    source_faces: &[ResolvedFace],
    thermal_faces: &[ResolvedFace],
    domain: [[f64; 3]; 2],
    cell_mm: [f64; 3],
    wetted: [bool; 6],
) -> Result<()> {
    let mut sources = Vec::new();
    for (source, face) in heat.face_sources.iter().zip(source_faces) {
        if face.triangles.is_empty() {
            return Err(bad(format!("{}: the face has no triangles in the model's tessellation", source.name)));
        }
        sources.push(serde_json::json!({"name": source.name, "body": case_scene::part_name(&face.body), "triangles": face.triangles, "power_w": source.power_w}));
    }
    let mut boundaries: Vec<Value> = Vec::new();
    let mut given: [Option<&str>; 6] = Default::default();
    for (condition, face) in heat.face_conditions.iter().zip(thermal_faces) {
        let Some(side) = face.plane.and_then(|plane| patch_of(plane, domain, cell_mm)).filter(|p| wetted[*p]) else {
            return Err(bad(format!(
                "{}: a thermal condition goes on a flat face on the domain boundary; on a face inside the domain, give a heat source",
                condition.name
            )));
        };
        if let Some(other) = given[side] {
            return Err(bad(format!("{} and {other} are on the same domain face ({}); keep one", condition.name, PATCHES[side])));
        }
        given[side] = Some(&condition.name);
        let row = match condition.condition {
            ThermalFace::Adiabatic => serde_json::json!({"kind": "adiabatic"}),
            ThermalFace::Temperature { temperature_k } => serde_json::json!({"kind": "temperature", "temperature_k": temperature_k}),
            ThermalFace::HeatFlux { w_per_m2 } => serde_json::json!({"kind": "heat_flux", "w_per_m2": w_per_m2}),
            ThermalFace::Convection { h_w_per_m2_k, ambient_k } => {
                serde_json::json!({"kind": "convection", "h_w_per_m2_k": h_w_per_m2_k, "ambient_k": ambient_k})
            }
        };
        boundaries.push(serde_json::json!({"patch": PATCHES[side], "condition": row}));
    }
    if let Some(thermal) = case.get_mut("thermal").and_then(Value::as_object_mut) {
        if !sources.is_empty() {
            thermal.insert("face_sources".into(), Value::Array(sources));
        }
        if !boundaries.is_empty() {
            thermal.insert("boundaries".into(), Value::Array(boundaries));
        }
    }
    Ok(())
}

/// Build the solver case for `link` from `base` (the study's settings: fluid, time step,
/// thermal records) and the model as `input` describes it now.
pub fn build(base: &Value, link: &ModelLink, input: &ModelInput) -> Result<Built> {
    link.validate()?;
    if base.get("solver").and_then(Value::as_str) != Some("collocated_flow3d") {
        return Err(bad("a model study needs a 3-D flow case (collocated_flow3d)"));
    }
    if input.faces.len() != link.face_refs().len() {
        return Err(bad("every face the study names needs its resolved face"));
    }
    let (opening_faces, rest) = input.faces.split_at(link.openings.len());
    let refine_count = link.refinements.iter().filter(|r| matches!(r.target, RefineTarget::Face { .. })).count();
    let (refine_faces, rest) = rest.split_at(refine_count);
    let goal_count = link.goals.iter().filter(|g| g.face.is_some()).count();
    let (goal_faces, rest) = rest.split_at(goal_count.min(rest.len()));
    let source_count = link.heat.as_ref().map_or(0, |h| h.face_sources.len());
    let (source_faces, rest) = rest.split_at(source_count.min(rest.len()));
    let thermal_count = link.heat.as_ref().map_or(0, |h| h.face_conditions.len());
    let (thermal_faces, fan_faces) = rest.split_at(thermal_count.min(rest.len()));
    let margin = match link.analysis {
        Analysis::Internal => 0.0,
        Analysis::External => link.margin_mm,
    };
    let [low, high] = input.bounds_mm;
    if (0..3).any(|a| !(low[a].is_finite() && high[a].is_finite() && high[a] > low[a])) {
        return Err(bad("the bodies have no extent along every axis; a flat body cannot bound a fluid domain"));
    }
    let domain = [low.map(|v| v - margin), high.map(|v| v + margin)];
    let options = ImportOptions {
        // Every body is its exact triangle surface, so a condition can name any of its faces.
        representation: Representation::Surface,
        fit_domain: true,
        margin_mm: margin,
        cell_mm: link.cell_mm,
        replace_parts: true,
        domain_mm: Some(domain),
        cell_scale: link.mesh_level.filter(|_| link.cell_mm.is_none()).map(cell_scale),
        ..Default::default()
    };
    // Conditions on body faces belong to the last build: drop them and their patches' rows,
    // so a removed condition leaves nothing behind. The domain faces' rows stay.
    let mut base = base.clone();
    if let Some(object) = base.as_object_mut() {
        object.remove("face_conditions");
        // Goals, heat transfer and the fan of a study of the model are the link's.
        object.remove("goals");
        object.remove("thermal");
        object.remove("fan");
        for key in ["velocity_boundaries", "pressure_boundaries"] {
            if let Some(rows) = object.get_mut(key).and_then(Value::as_array_mut) {
                rows.retain(|row| row["patch"].as_str().is_some_and(|patch| PATCHES.contains(&patch)));
            }
        }
    }
    // Local refinements on the imported grid.
    let refine = |imported: &case_scene::Imported| -> Result<Value> {
        let mut case = imported.edit.case.clone();
        let mut refine_face = refine_faces.iter();
        for r in &link.refinements {
            let edit = match &r.target {
                RefineTarget::Body { body } => case_scene::refine_around(&case, &[case_scene::part_name(body)], r.factor, r.margin_cells),
                RefineTarget::Region { low_mm, high_mm } => {
                    let (dimensions, lengths) = case_scene::grid(&case)?;
                    let pad = [0, 1, 2].map(|a| r.margin_cells * lengths[a] / dimensions[a] as f64);
                    let origin = imported.origin_mm;
                    // A box edge on a grid node (an adaptive box is whole cells) must not reach
                    // into the next cell through rounding: a nanometre inside.
                    const INSIDE_M: f64 = 1e-9;
                    case_scene::refine_region(
                        &case,
                        [0, 1, 2].map(|a| (low_mm[a] - origin[a]) / 1000.0 - pad[a] + INSIDE_M),
                        [0, 1, 2].map(|a| (high_mm[a] - origin[a]) / 1000.0 + pad[a] - INSIDE_M),
                        r.factor,
                    )
                }
                RefineTarget::Face { .. } => {
                    let face = refine_face.next().ok_or_else(|| bad("a refinement's face is missing"))?;
                    let (low, high) =
                        face_box(&case, face).ok_or_else(|| bad(format!("{}: the face has no triangles in the model's tessellation", r.name)))?;
                    let (dimensions, lengths) = case_scene::grid(&case)?;
                    let pad = [0, 1, 2].map(|a| r.margin_cells * lengths[a] / dimensions[a] as f64);
                    // A flat face has no thickness along its normal: give the box a sliver there
                    // so the cells on both sides of it count, and only there.
                    let sliver = [0, 1, 2].map(|a| if high[a] - low[a] <= 1e-9 * lengths[a] { 1e-9 * lengths[a] } else { 0.0 });
                    case_scene::refine_region(
                        &case,
                        [0, 1, 2].map(|a| low[a] - pad[a] - sliver[a]),
                        [0, 1, 2].map(|a| high[a] + pad[a] + sliver[a]),
                        r.factor,
                    )
                }
            }
            .map_err(|e| match e {
                // Kept apart: an automatic grid makes room for it below.
                Error::Limit(m) => Error::Limit(format!("{}: {m}", r.name)),
                e => bad(format!("{}: {e}", r.name)),
            })?;
            case = edit.case;
        }
        Ok(case)
    };
    // An automatic grid is the finest that fits: when the refinements would take it past the
    // cell limit, its base cells grow (15 % a time, to at most 8 times the size) until they fit.
    // A typed cell size is the user's: the refinement is refused instead.
    let automatic = link.cell_mm.is_none();
    let mut grow = 1.0f64;
    let (imported, mut case) = loop {
        let mut options = options.clone();
        if automatic {
            options.cell_scale = Some(options.cell_scale.unwrap_or(1.0) * grow);
        }
        let imported = case_scene::import_components(&base, &input.obj, &options)?;
        match refine(&imported) {
            Ok(case) => break (imported, case),
            Err(Error::Limit(_)) if automatic && grow * 1.15 <= 8.0 => grow *= 1.15,
            Err(e) => return Err(e),
        }
    };
    let mut notes = imported.notes.clone();
    if grow > 1.0 {
        notes.push(format!(
            "the local refinements needed room under the cell limit: the automatic cells are {grow:.2} times the size the bodies ask for"
        ));
    }
    if !link.refinements.is_empty() {
        case_scene::fit_time_step(&mut case, &mut notes);
    }
    if let Some(heat) = &link.heat {
        let thermal = thermal_block(&case, heat, &mut notes)?;
        if let Some(object) = case.as_object_mut() {
            object.insert("thermal".into(), thermal);
        }
    }
    let grid_nodes = case_scene::nodes(&case)?;
    let cell_mm = [0, 1, 2].map(|a| grid_nodes[a].windows(2).map(|w| w[1] - w[0]).fold(f64::INFINITY, f64::min) * 1000.0);
    if link.cut_cells {
        for (key, what) in [
            ("thermal", "heat transfer"),
            ("transferred", "transferred boundary conditions"),
            ("rotating_regions", "rotating regions"),
            ("internal_fans", "internal fans"),
            ("porous_blocks", "porous blocks"),
            ("fan", "an external fan"),
        ] {
            if case.get(key).is_some_and(|v| !v.is_null() && v.as_array().is_none_or(|rows| !rows.is_empty())) {
                return Err(bad(format!("cut cells do not work with {what} yet; turn cut cells off")));
            }
        }
        if let Some(object) = case.as_object_mut() {
            object.insert("cut_cells".into(), Value::Bool(true));
        }
    } else if let Some(object) = case.as_object_mut() {
        object.remove("cut_cells");
    }
    // A domain face that the bodies cover entirely has no fluid next to it: it takes no
    // condition at all (the solver refuses one there). With cut cells the solver looks at the
    // thin layer next to the face, so this does too.
    let wetted = if link.cut_cells { wetted_layers(&case)? } else { crate::fluid_volume::fluid_volume(&case)?.wetted };
    for (patch, name) in PATCHES.iter().enumerate() {
        if wetted[patch] {
            continue;
        }
        if let Some(object) = case.as_object_mut() {
            for key in ["velocity_boundaries", "pressure_boundaries"] {
                if let Some(rows) = object.get_mut(key).and_then(Value::as_array_mut) {
                    rows.retain(|row| row["patch"].as_str() != Some(name));
                }
            }
        }
    }
    if link.analysis == Analysis::Internal {
        // Inside an enclosure every domain face with fluid at it is a wall unless a condition
        // opens it.
        for (patch, name) in PATCHES.iter().enumerate() {
            if wetted[patch] && case.pointer("/fan/patch").and_then(Value::as_str) != Some(*name) {
                case = case_scene::set_patch_role(&case, patch, RoleInput::Wall)?.case;
            }
        }
    }
    let mut placed: Vec<(String, String)> = Vec::new();
    let mut given: [Option<(String, OpeningRole)>; 6] = Default::default();
    let objects = case_scene::parse_obj(&input.obj)?;
    let mut taken: std::collections::BTreeSet<String> = PATCHES.iter().map(|p| (*p).to_string()).collect();
    taken.extend(objects.iter().map(|o| case_scene::part_name(&o.name)));
    for (opening, face) in link.openings.iter().zip(opening_faces) {
        // A flat face on the domain boundary opens that domain face; any other face becomes a
        // patch on its own body.
        let Some(patch) = face.plane.and_then(|plane| patch_of(plane, domain, cell_mm)).filter(|p| wetted[*p]) else {
            let patch = body_face_condition(&mut case, opening, face, &objects, &mut taken)?;
            placed.push((opening.name.clone(), patch));
            continue;
        };
        if let Some((other, role)) = &given[patch]
            && *role != opening.role
        {
            return Err(bad(format!(
                "{} and {other} lie on the same domain face ({}) with different conditions; give them the same condition",
                opening.name, PATCHES[patch]
            )));
        }
        given[patch] = Some((opening.name.clone(), opening.role));
        case = case_scene::set_patch_role(&case, patch, opening.role.role())?.case;
        placed.push((opening.name.clone(), PATCHES[patch].to_string()));
    }
    if let Some(heat) = &link.heat {
        heat_on_faces(&mut case, heat, source_faces, thermal_faces, domain, cell_mm, wetted)?;
    }
    if let (Some(fan), Some(face)) = (&link.fan, fan_faces.first()) {
        put_fan(&mut case, fan, face, domain, cell_mm, wetted, &given)?;
    }
    // Goals: a face goal reads the patch the face is on, or a wall patch made for it.
    let mut goal_face = goal_faces.iter();
    let mut goals = Vec::new();
    for goal in &link.goals {
        let row = match goal.kind {
            ModelGoalKind::FaceVolumeFlow | ModelGoalKind::FaceMeanPressure => {
                let face = goal_face.next().ok_or_else(|| bad("a goal's face is missing"))?;
                let wanted = goal.face.as_ref();
                let patch = match link.openings.iter().position(|o| Some(&o.face) == wanted) {
                    Some(i) => placed[i].1.clone(),
                    None => match face.plane.and_then(|plane| patch_of(plane, domain, cell_mm)).filter(|p| wetted[*p]) {
                        Some(side) => PATCHES[side].to_string(),
                        None => {
                            let wall = Opening {
                                name: format!("goal_{}", goal.name),
                                face: wanted.cloned().ok_or_else(|| bad("a face goal needs its face"))?,
                                role: OpeningRole::Wall,
                            };
                            body_face_condition(&mut case, &wall, face, &objects, &mut taken)?
                        }
                    },
                };
                let kind = if goal.kind == ModelGoalKind::FaceVolumeFlow { "patch_volume_flow" } else { "patch_mean_pressure" };
                serde_json::json!({"name": goal.name, "kind": kind, "target": patch})
            }
            ModelGoalKind::PointPressure => {
                let p = goal.point_mm.unwrap_or_default();
                let at = [0, 1, 2].map(|a| (p[a] - imported.origin_mm[a]) / 1000.0);
                serde_json::json!({"name": goal.name, "kind": "point_pressure", "point": at})
            }
            ModelGoalKind::Equation => serde_json::json!({"name": goal.name, "kind": "equation", "expression": goal.expression}),
            ModelGoalKind::BodyMaxTemperature | ModelGoalKind::BodyMeanTemperature => {
                let body = goal.body.as_deref().unwrap_or_default();
                let part = case_scene::part_name(body);
                if !objects.iter().any(|o| o.name == body) {
                    return Err(bad(format!("goal {}: {body} is not a body of the study", goal.name)));
                }
                let kind = if goal.kind == ModelGoalKind::BodyMaxTemperature { "solid_max_temperature" } else { "solid_mean_temperature" };
                serde_json::json!({"name": goal.name, "kind": kind, "target": part})
            }
            ModelGoalKind::PointTemperature => {
                let p = goal.point_mm.unwrap_or_default();
                let at = [0, 1, 2].map(|a| (p[a] - imported.origin_mm[a]) / 1000.0);
                serde_json::json!({"name": goal.name, "kind": "point_temperature", "point": at})
            }
            // The global kinds have the companion's names.
            global => serde_json::json!({"name": goal.name, "kind": global}),
        };
        goals.push(row);
    }
    if !goals.is_empty()
        && let Some(object) = case.as_object_mut()
    {
        object.insert("goals".into(), Value::Array(goals));
    }
    // How much of the model the grid keeps, and where the fluid is.
    let volume = crate::fluid_volume::fluid_volume(&case)?;
    // What the solver will hold: whole cells for the staircase, the cut-cell sub-lattice with cut cells.
    let solid_m3 = if link.cut_cells { cut_solid_volume(&case)? } else { volume.solid_volume_m3 };
    let solid_fraction = if input.volume_mm3 > 0.0 { solid_m3 * 1e9 / input.volume_mm3 } else { 0.0 };
    if !link.cut_cells && (solid_fraction < 0.8 || solid_fraction > 1.25) {
        notes.push(format!(
            "the grid's solid cells hold {:.0} % of the bodies' volume, so walls thinner than a cell are lost or thickened; set a smaller cell size",
            solid_fraction * 100.0
        ));
    }
    if volume.fluid_cells == 0 {
        notes.push("no fluid cell is left: the bodies fill the domain".into());
    }
    if link.analysis == Analysis::Internal && placed.is_empty() {
        notes.push("no face has a condition, so the fluid is sealed in; give an inlet and an outlet".into());
    }
    let issues = crate::validate_case(&case).err().map(|e| vec![e.to_string()]).unwrap_or_default();
    Ok(Built {
        digest: input.digest(link)?,
        case,
        origin_mm: imported.origin_mm,
        placed,
        solid_fraction,
        fluid_regions: volume.regions.len() + volume.unlisted_regions,
        notes,
        issues,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Welded OBJ vertices and faces.
    #[derive(Default)]
    struct Obj {
        points: Vec<[f64; 3]>,
        faces: String,
    }

    impl Obj {
        fn index(&mut self, p: [f64; 3]) -> usize {
            match self.points.iter().position(|q| *q == p) {
                Some(i) => i + 1,
                None => {
                    self.points.push(p);
                    self.points.len()
                }
            }
        }
        /// A planar convex quad wound so its normal points along `outward`.
        fn quad(&mut self, mut corners: [[f64; 3]; 4], outward: [f64; 3]) {
            let [a, b, c, _] = corners;
            let (u, v) = ([0, 1, 2].map(|i| b[i] - a[i]), [0, 1, 2].map(|i| c[i] - a[i]));
            let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            if n[0] * outward[0] + n[1] * outward[1] + n[2] * outward[2] < 0.0 {
                corners.reverse();
            }
            let ids = corners.map(|p| self.index(p));
            self.faces.push_str(&format!("f {} {} {} {}\n", ids[0], ids[1], ids[2], ids[3]));
        }
        fn text(&self, name: &str) -> String {
            let mut text = format!("o {name}\n");
            for p in &self.points {
                text.push_str(&format!("v {} {} {}\n", p[0], p[1], p[2]));
            }
            text + &self.faces
        }
    }

    /// A square duct along x, 100 mm long, 40 × 40 mm outside with 4 mm walls, open at both
    /// ends: volume 100 × (40² − 32²) = 57 600 mm³.
    fn duct() -> ModelInput {
        let mut obj = Obj::default();
        let (l, o, i0, i1) = (100.0, 40.0, 4.0, 36.0);
        // Outer walls, then inner walls facing the channel.
        for (lo, hi, sign) in [(0.0, o, 1.0), (i0, i1, -1.0)] {
            obj.quad([[0.0, lo, lo], [l, lo, lo], [l, hi, lo], [0.0, hi, lo]], [0.0, 0.0, -sign]);
            obj.quad([[0.0, lo, hi], [l, lo, hi], [l, hi, hi], [0.0, hi, hi]], [0.0, 0.0, sign]);
            obj.quad([[0.0, lo, lo], [l, lo, lo], [l, lo, hi], [0.0, lo, hi]], [0.0, -sign, 0.0]);
            obj.quad([[0.0, hi, lo], [l, hi, lo], [l, hi, hi], [0.0, hi, hi]], [0.0, sign, 0.0]);
        }
        // The two end frames, mitred like a picture frame so every edge is shared by two faces.
        for (x, nx) in [(0.0, -1.0), (l, 1.0)] {
            let p = |y: f64, z: f64| [x, y, z];
            obj.quad([p(0.0, 0.0), p(o, 0.0), p(i1, i0), p(i0, i0)], [nx, 0.0, 0.0]);
            obj.quad([p(0.0, o), p(i0, i1), p(i1, i1), p(o, o)], [nx, 0.0, 0.0]);
            obj.quad([p(0.0, 0.0), p(i0, i0), p(i0, i1), p(0.0, o)], [nx, 0.0, 0.0]);
            obj.quad([p(o, 0.0), p(o, o), p(i1, i1), p(i1, i0)], [nx, 0.0, 0.0]);
        }
        ModelInput { obj: obj.text("Duct"), bounds_mm: [[0.0; 3], [l, o, o]], volume_mm3: l * (o * o - 32.0 * 32.0), faces: vec![] }
    }

    fn channel() -> Value {
        serde_json::from_str(include_str!("../data/templates/flow-channel.json")).unwrap_or(Value::Null)
    }

    fn end(x: f64, nx: f64) -> ResolvedFace {
        ResolvedFace { plane: Some(([nx, 0.0, 0.0], [x, 20.0, 2.0])), body: "Duct".into(), ..Default::default() }
    }

    fn link(openings: &[(&str, OpeningRole)]) -> ModelLink {
        ModelLink {
            cell_mm: Some(4.0),
            openings: openings
                .iter()
                .map(|(name, role)| Opening {
                    name: (*name).into(),
                    face: FaceRef { body: "Duct".into(), face: format!("{name} face") },
                    role: *role,
                })
                .collect(),
            ..Default::default()
        }
    }

    /// Heat transfer in a study of the model: the duct conducts with its material and power,
    /// one inner wall releases heat, the inlet end holds a temperature, a wall goal on an inner
    /// face still works with heat on, and body temperature goals read the duct.
    #[test]
    fn heat_goes_on_the_bodies_their_faces_and_the_domain_boundary() {
        let inlet = OpeningRole::Inlet { speed_m_per_s: 0.5 };
        let outlet = OpeningRole::Outlet { pressure_m2_per_s2: 0.0 };
        let mut heated = link(&[("in", inlet), ("out", outlet)]);
        let face = |name: &str| FaceRef { body: "Duct".into(), face: name.into() };
        let goal = |name: &str, kind, face: Option<FaceRef>, body: Option<&str>| ModelGoal {
            name: name.into(),
            kind,
            face,
            point_mm: (kind == ModelGoalKind::PointTemperature).then_some([50.0, 20.0, 20.0]),
            expression: None,
            body: body.map(str::to_string),
        };
        heated.goals = vec![
            goal("wall_p", ModelGoalKind::FaceMeanPressure, Some(face("inner y")), None),
            goal("duct_max", ModelGoalKind::BodyMaxTemperature, None, Some("Duct")),
            goal("duct_mean", ModelGoalKind::BodyMeanTemperature, None, Some("Duct")),
            goal("mid_t", ModelGoalKind::PointTemperature, None, None),
        ];
        heated.heat = Some(ModelHeat {
            inflow_temperature_k: 300.0,
            materials: vec![BodyHeat {
                body: "Duct".into(),
                conductivity_w_per_m_k: 200.0,
                power_w: 1.0,
                density_kg_per_m3: Some(2700.0),
                specific_heat_j_per_kg_k: Some(900.0),
                max_temperature_k: Some(360.0),
            }],
            face_sources: vec![FaceHeat { name: "heater".into(), face: face("inner y"), power_w: 5.0 }],
            face_conditions: vec![FaceThermal { name: "hot inlet".into(), face: face("in face"), condition: ThermalFace::Temperature { temperature_k: 310.0 } }],
            ..Default::default()
        });
        let mut input = duct();
        let inner = ResolvedFace { name: "inner y".into(), body: "Duct".into(), triangles: vec![12, 13], plane: Some(([0.0, 1.0, 0.0], [50.0, 4.0, 20.0])) };
        input.faces = vec![
            ResolvedFace { name: "in".into(), ..end(0.0, -1.0) },
            ResolvedFace { name: "out".into(), ..end(100.0, 1.0) },
            inner.clone(),
            ResolvedFace { name: "heater".into(), ..inner.clone() },
            ResolvedFace { name: "hot inlet".into(), ..end(0.0, -1.0) },
        ];
        let built = build(&channel(), &heated, &input).expect("build");
        assert!(built.issues.is_empty(), "{:?}", built.issues);
        if let Ok(path) = std::env::var("SC_DUMP_HEATED") { let _ = std::fs::write(path, built.case.to_string()); }
        let thermal = &built.case["thermal"];
        assert_eq!(thermal["inflow_temperature_k"], serde_json::json!(300.0));
        assert_eq!(thermal["solids"], serde_json::json!([{
            "obstacle": "Duct", "conductivity_w_per_m_k": [200.0, 200.0, 200.0], "power_w": 1.0,
            "density_kg_per_m3": 2700.0, "specific_heat_j_per_kg_k": 900.0, "max_temperature_k": 360.0
        }]));
        assert_eq!(thermal["face_sources"], serde_json::json!([{"name": "heater", "body": "Duct", "triangles": [12, 13], "power_w": 5.0}]));
        assert_eq!(thermal["boundaries"], serde_json::json!([{"patch": "left", "condition": {"kind": "temperature", "temperature_k": 310.0}}]));
        let goals = built.case["goals"].as_array().cloned().unwrap_or_default();
        let kinds: Vec<&str> = goals.iter().filter_map(|g| g["kind"].as_str()).collect();
        assert_eq!(kinds, ["patch_mean_pressure", "solid_max_temperature", "solid_mean_temperature", "point_temperature"]);
        assert_eq!(goals[1]["target"], "Duct");
        let needs: Vec<&str> = crate::airflow::requirements(&built.case).iter().map(|r| r.capability).collect();
        assert!(needs.contains(&"face_heat_sources_v1") && needs.contains(&"conjugate_heat_v1"), "{needs:?}");
        // Heat settings are part of the digest.
        let mut hotter = heated.clone();
        if let Some(heat) = hotter.heat.as_mut() {
            heat.inflow_temperature_k = 320.0;
        }
        assert_ne!(built.digest, input.digest(&hotter).expect("digest"));
        // Without a material the duct conducts like the starter material, and says so.
        let mut bare = heated.clone();
        if let Some(heat) = bare.heat.as_mut() {
            heat.materials.clear();
        }
        let starter = build(&channel(), &bare, &input).expect("starter");
        assert_eq!(starter.case["thermal"]["solids"][0]["conductivity_w_per_m_k"], serde_json::json!([1.0, 1.0, 1.0]));
        assert!(starter.notes.iter().any(|n| n.contains("starter material")), "{:?}", starter.notes);
        // A thermal condition on a face inside the domain is refused with the reason; a material
        // on a body the study does not hold (a part a selection or a size limit left out) is kept
        // and noted.
        let mut inside = input.clone();
        inside.faces[4] = ResolvedFace { name: "hot inlet".into(), ..inner.clone() };
        let why = build(&channel(), &heated, &inside).expect_err("inside").to_string();
        assert!(why.contains("domain boundary"), "{why}");
        let mut stranger = heated.clone();
        if let Some(heat) = stranger.heat.as_mut() {
            heat.materials[0].body = "Lid".into();
        }
        let left = build(&channel(), &stranger, &input).expect("a left-out part's material is kept");
        assert!(left.notes.iter().any(|n| n.starts_with("Lid have a material but the study leaves them out")), "{:?}", left.notes);
        assert_eq!(left.case["thermal"]["solids"][0]["conductivity_w_per_m_k"], serde_json::json!([1.0, 1.0, 1.0]), "the duct itself has none now");
        // Temperature goals need heat on.
        let mut cold = heated.clone();
        cold.heat = None;
        assert!(cold.validate().expect_err("cold").to_string().contains("heat transfer on"));
    }

    /// The mesh level scales the automatic cells (each level by √2, 3 is automatic), a typed
    /// cell size wins over it, and a region of the design refines the cells it overlaps.
    #[test]
    fn mesh_levels_scale_the_automatic_cells_and_regions_refine_a_box() {
        let mut auto = link(&[]);
        auto.cell_mm = None;
        let cells = |link: &ModelLink| -> usize {
            let built = build(&channel(), link, &duct()).expect("build");
            built.case["mesh"]["dimensions"].as_array().into_iter().flatten().filter_map(Value::as_u64).product::<u64>() as usize
        };
        let automatic = cells(&auto);
        let mut level = auto.clone();
        level.mesh_level = Some(AUTOMATIC_MESH_LEVEL);
        assert_eq!(cells(&level), automatic);
        level.mesh_level = Some(5);
        let finer = cells(&level);
        level.mesh_level = Some(1);
        let coarser = cells(&level);
        // The automatic cells make the duct's 4 mm walls two cells thick: 2 mm, 50 × 20 × 20.
        assert_eq!(automatic, 50 * 20 * 20);
        // Two levels double the cell size: about an eighth of the cells (whole cells per axis).
        // Two levels finer would be eight times as many, past the cell limit, which caps them.
        assert!(coarser * 5 < automatic, "{coarser} {automatic}");
        assert!(finer > automatic && finer as u64 <= crate::MAX_CELLS, "{automatic} {finer}");
        level.cell_mm = Some(4.0);
        assert_eq!(cells(&level), 25 * 10 * 10, "a typed cell size wins");
        level.mesh_level = Some(9);
        assert!(level.validate().is_err());
        assert!((cell_scale(5) - 0.5).abs() < 1e-12 && (cell_scale(1) - 2.0).abs() < 1e-12);
        // A region from x = 40 to 60 mm splits the five 4 mm cells it overlaps along x, and
        // every cell along y and z (the box spans the duct's section).
        let mut region = link(&[]);
        region.refinements = vec![Refinement {
            name: "adapt_1".into(),
            target: RefineTarget::Region { low_mm: [41.0, 1.0, 1.0], high_mm: [59.0, 39.0, 39.0] },
            factor: 2,
            margin_cells: 0.0,
        }];
        let built = build(&channel(), &region, &duct()).expect("region");
        assert_eq!(built.case["mesh"]["dimensions"], serde_json::json!([30, 20, 20]));
        region.refinements[0].target = RefineTarget::Region { low_mm: [60.0, 1.0, 1.0], high_mm: [40.0, 39.0, 39.0] };
        assert!(region.validate().is_err());
    }

    /// A fan on the duct's inlet end blows against the outlet, and global goals read the whole
    /// fluid.
    #[test]
    fn a_fan_on_an_end_face_and_global_goals() {
        let outlet = OpeningRole::Outlet { pressure_m2_per_s2: 0.0 };
        let mut fanned = link(&[("out", outlet)]);
        let curve = vec![[0.0, 4.0], [4e-4, 3.0], [1e-3, 0.0]];
        fanned.fan = Some(ModelFan { name: "blower".into(), face: FaceRef { body: "Duct".into(), face: "left end".into() }, curve, total_pressure: false });
        let global = |name: &str, kind| ModelGoal { name: name.into(), kind, face: None, point_mm: None, expression: None, body: None };
        fanned.goals = vec![global("v_max", ModelGoalKind::GlobalMaxSpeed), global("p_mean", ModelGoalKind::GlobalMeanPressure)];
        let mut input = duct();
        input.faces = vec![ResolvedFace { name: "out".into(), ..end(100.0, 1.0) }, ResolvedFace { name: "left end".into(), ..end(0.0, -1.0) }];
        let built = build(&channel(), &fanned, &input).expect("fan");
        assert!(built.issues.is_empty(), "{:?}", built.issues);
        let fan = &built.case["fan"];
        assert_eq!((fan["patch"].clone(), fan["reference_patch"].clone()), (serde_json::json!("left"), serde_json::json!("right")));
        assert_eq!(fan["curve"][1], serde_json::json!({"volume_flow_m3_per_s": 4e-4, "pressure_rise_pa": 3.0}));
        let rows = |key: &str| built.case[key].as_array().into_iter().flatten().filter(|r| r["patch"] == "left").count();
        assert_eq!(rows("velocity_boundaries") + rows("pressure_boundaries"), 0, "the fan owns its face");
        assert_eq!(built.case["goals"], serde_json::json!([{"name": "v_max", "kind": "global_max_speed"}, {"name": "p_mean", "kind": "global_mean_pressure"}]));
        let needs: Vec<&str> = crate::airflow::requirements(&built.case).iter().map(|r| r.capability).collect();
        assert!(needs.contains(&"fan_curve_v1") && needs.contains(&"goals_v3"), "{needs:?}");
        // Without an outlet the fan has nothing to blow against; a rising curve is refused.
        let mut alone = fanned.clone();
        alone.openings.clear();
        let mut lone = input.clone();
        lone.faces.remove(0);
        assert!(build(&channel(), &alone, &lone).expect_err("alone").to_string().contains("outlet"));
        if let Some(fan) = fanned.fan.as_mut() {
            fan.curve = vec![[0.0, 1.0], [1e-3, 2.0]];
        }
        assert!(fanned.validate().is_err());
        // Global temperatures need heat on.
        let mut hot = link(&[]);
        hot.goals = vec![global("t_max", ModelGoalKind::GlobalMaxFluidTemperature)];
        assert!(hot.validate().is_err());
    }

    #[test]
    fn an_internal_duct_takes_its_box_and_opens_its_end_faces() {
        let inlet = OpeningRole::Inlet { speed_m_per_s: 0.5 };
        let outlet = OpeningRole::Outlet { pressure_m2_per_s2: 0.0 };
        let link = link(&[("in", inlet), ("out", outlet)]);
        let mut input = duct();
        input.faces = vec![ResolvedFace { name: "in".into(), ..end(0.0, -1.0) }, ResolvedFace { name: "out".into(), ..end(100.0, 1.0) }];
        let built = build(&channel(), &link, &input).expect("build");
        assert!(built.issues.is_empty(), "{:?}", built.issues);
        // The domain is the duct's box exactly, in 4 mm cells: 25 × 10 × 10.
        assert_eq!(built.case["mesh"]["dimensions"], serde_json::json!([25, 10, 10]));
        assert_eq!(built.origin_mm, [0.0; 3]);
        assert_eq!(built.placed, vec![("in".to_string(), "left".to_string()), ("out".to_string(), "right".to_string())]);
        // The 4 mm walls are one cell each, so the grid holds the duct's volume exactly.
        assert!((built.solid_fraction - 1.0).abs() < 1e-9, "{}", built.solid_fraction);
        assert_eq!(built.fluid_regions, 1);
        // Inflow along +x on the left, ambient pressure on the right. The duct's walls cover
        // the other four domain faces entirely, so no fluid touches them and they take no
        // condition (the solver refuses one on a face without fluid).
        let velocity = |patch: &str| built.case["velocity_boundaries"].as_array().into_iter().flatten().find(|r| r["patch"] == patch).cloned();
        assert_eq!(velocity("left").map(|r| r["value"]["x"].clone()), Some(serde_json::json!(0.5)));
        assert_eq!(velocity("right").map(|r| r["kind"].clone()), Some(serde_json::json!("zero_gradient")));
        for patch in ["bottom", "top", "front", "back"] {
            assert_eq!(velocity(patch), None, "{patch}");
            let pressure = built.case["pressure_boundaries"].as_array().into_iter().flatten().find(|r| r["patch"] == patch).cloned();
            assert_eq!(pressure, None, "{patch}");
        }
        // The same input gives the same digest; another setting gives another.
        assert_eq!(built.digest, input.digest(&link).expect("digest"));
        let mut finer = link.clone();
        finer.cell_mm = Some(2.0);
        assert_ne!(built.digest, input.digest(&finer).expect("digest"));
    }

    #[test]
    fn a_face_inside_the_domain_becomes_a_patch_on_its_body() {
        let outlet = OpeningRole::Outlet { pressure_m2_per_s2: 0.0 };
        let mut input = duct();
        // The inner channel wall at y = 4 mm faces into the fluid, inside the domain: quad 6 of
        // the duct, its triangles 12 and 13, outward along +y.
        input.faces = vec![ResolvedFace {
            name: "Side vent".into(),
            plane: Some(([0.0, 1.0, 0.0], [50.0, 4.0, 20.0])),
            body: "Duct".into(),
            triangles: vec![12, 13],
        }];
        let built = build(&channel(), &link(&[("Side vent", outlet)]), &input).expect("inside");
        assert!(built.issues.is_empty(), "{:?}", built.issues);
        assert_eq!(built.placed, vec![("Side vent".to_string(), "Side_vent".to_string())]);
        assert_eq!(built.case["face_conditions"], serde_json::json!([{"name":"Side_vent","body":"Duct","triangles":[12,13]}]));
        let pressure = built.case["pressure_boundaries"].as_array().into_iter().flatten().find(|r| r["patch"] == "Side_vent").cloned();
        assert_eq!(pressure, Some(serde_json::json!({"patch":"Side_vent","kind":"fixed_value","value":0.0})));
        // An inlet there blows along the face's outward normal, into the channel.
        let inlet = OpeningRole::Inlet { speed_m_per_s: 0.5 };
        let built = build(&channel(), &link(&[("Side vent", inlet)]), &input).expect("inlet");
        let velocity = built.case["velocity_boundaries"].as_array().into_iter().flatten().find(|r| r["patch"] == "Side_vent").cloned();
        assert_eq!(velocity.map(|r| r["value"].clone()), Some(serde_json::json!({"x":0.0,"y":0.5,"z":0.0})));
        // Building again without the condition leaves none of it behind.
        let mut without = duct();
        without.faces = vec![];
        let rebuilt = build(&built.case, &link(&[]), &without).expect("rebuild");
        assert_eq!(rebuilt.case.get("face_conditions"), None);
        let rows = |key: &str| rebuilt.case[key].as_array().into_iter().flatten().filter(|r| r["patch"] == "Side_vent").count();
        assert_eq!(rows("velocity_boundaries") + rows("pressure_boundaries"), 0);
        // The same outlet with heat transfer is refused, not left uncoupled.
        let mut heated = link(&[("Side vent", outlet)]);
        heated.heat = Some(ModelHeat::default());
        let why = build(&channel(), &heated, &input).expect_err("thermal").to_string();
        assert!(why.contains("heat transfer"), "{why}");
    }

    #[test]
    fn faces_without_triangles_and_clashes_are_refused() {
        let inlet = OpeningRole::Inlet { speed_m_per_s: 0.5 };
        let mut input = duct();
        // A face the tessellation does not hold (a curved face too small to keep, say).
        input.faces = vec![ResolvedFace { name: "in".into(), body: "Duct".into(), ..Default::default() }];
        let why = build(&channel(), &link(&[("in", inlet)]), &input).expect_err("no triangles").to_string();
        assert!(why.contains("no triangles"), "{why}");
        // Two faces on the left domain face with different conditions.
        let wall = OpeningRole::Wall;
        input.faces = vec![ResolvedFace { name: "a".into(), ..end(0.0, -1.0) }, ResolvedFace { name: "b".into(), ..end(0.0, -1.0) }];
        let why = build(&channel(), &link(&[("a", inlet), ("b", wall)]), &input).expect_err("clash").to_string();
        assert!(why.contains("same domain face"), "{why}");
        // A study with no conditions builds, sealed, and says so.
        input.faces = vec![];
        let sealed = build(&channel(), &link(&[]), &input).expect("sealed");
        assert!(sealed.notes.iter().any(|n| n.contains("sealed")), "{:?}", sealed.notes);
    }

    #[test]
    fn refinements_split_the_cells_around_a_body_or_a_face_and_refit_the_time_step() {
        let mut around_body = link(&[]);
        around_body.refinements =
            vec![Refinement { name: "all".into(), target: RefineTarget::Body { body: "Duct".into() }, factor: 2, margin_cells: 0.0 }];
        let built = build(&channel(), &around_body, &duct()).expect("body");
        // The duct fills the domain's box, so every 4 mm cell is split in two along each axis.
        assert_eq!(built.case["mesh"]["dimensions"], serde_json::json!([50, 20, 20]));
        assert!((built.solid_fraction - 1.0).abs() < 1e-9, "{}", built.solid_fraction);
        // The step is fitted to the 4 mm cells (0.002 s ÷ 19, see case_scene) and then to the
        // 2 mm ones: their limit 0.4 × 0.002² / (3 × 0.02) = 2.667e-5 s divides it by
        // ⌈1.0526e-4 / 2.667e-5⌉ = 4.
        assert!((built.case["config"]["dt"].as_f64().unwrap_or(0.0) - 0.002 / 76.0).abs() < 1e-15, "{}", built.case["config"]["dt"]);
        // Around the inner wall at y = 4 mm only (margin 0): the two y cells meeting it (0..4 and
        // 4..8 mm) split, 10 → 12; along x the wall runs the whole length (25 → 50), along z it
        // covers the 8 cells from 4 to 36 mm (10 → 18).
        let mut around_face = link(&[]);
        let wall = FaceRef { body: "Duct".into(), face: "inner y".into() };
        around_face.refinements = vec![Refinement { name: "wall".into(), target: RefineTarget::Face { face: wall }, factor: 2, margin_cells: 0.0 }];
        let mut input = duct();
        input.faces = vec![ResolvedFace { name: "inner y".into(), body: "Duct".into(), triangles: vec![12, 13], ..Default::default() }];
        let built = build(&channel(), &around_face, &input).expect("face");
        let dims = built.case["mesh"]["dimensions"].clone();
        assert_eq!(dims[1], serde_json::json!(12), "{dims}");
        assert_eq!((dims[0].clone(), dims[2].clone()), (serde_json::json!(50), serde_json::json!(18)), "{dims}");
    }

    /// The automatic grid (2 mm, the duct's walls two cells thick: 20 000 cells) split in two
    /// everywhere would be 160 000 cells: its base cells grow until the refinement fits. A
    /// typed cell size is refused instead.
    #[test]
    fn an_automatic_grid_makes_room_for_its_refinements() {
        let mut auto = link(&[]);
        auto.cell_mm = None;
        auto.refinements = vec![Refinement { name: "all".into(), target: RefineTarget::Body { body: "Duct".into() }, factor: 2, margin_cells: 0.0 }];
        let built = build(&channel(), &auto, &duct()).expect("room made");
        let cells: u64 = built.case["mesh"]["dimensions"].as_array().into_iter().flatten().filter_map(Value::as_u64).product();
        assert!(cells <= crate::MAX_CELLS && cells > crate::MAX_CELLS / 4, "{cells}");
        assert!(built.notes.iter().any(|n| n.contains("needed room under the cell limit")), "{:?}", built.notes);
        let mut typed = auto.clone();
        typed.cell_mm = Some(2.0);
        let refused = build(&channel(), &typed, &duct()).expect_err("over the limit");
        assert!(refused.to_string().contains("the local mesh would need"), "{refused}");
    }

    #[test]
    fn goals_read_condition_patches_domain_sides_and_wall_patches_made_for_them() {
        let inlet = OpeningRole::Inlet { speed_m_per_s: 0.5 };
        let mut goals = link(&[("in", inlet)]);
        let face = |name: &str| FaceRef { body: "Duct".into(), face: name.into() };
        goals.openings[0].face = face("left end");
        goals.goals = vec![
            ModelGoal { name: "inflow".into(), kind: ModelGoalKind::FaceVolumeFlow, face: Some(face("left end")), point_mm: None, expression: None, body: None },
            ModelGoal {
                name: "outflow".into(),
                kind: ModelGoalKind::FaceVolumeFlow,
                face: Some(face("right end")),
                point_mm: None,
                expression: None,
                body: None,
            },
            ModelGoal { name: "wall_p".into(), kind: ModelGoalKind::FaceMeanPressure, face: Some(face("inner y")), point_mm: None, expression: None, body: None },
            ModelGoal { name: "mid_p".into(), kind: ModelGoalKind::PointPressure, face: None, point_mm: Some([50.0, 20.0, 20.0]), expression: None, body: None },
            ModelGoal {
                name: "balance".into(),
                kind: ModelGoalKind::Equation,
                face: None,
                point_mm: None,
                expression: Some("inflow + outflow".into()),
                body: None,
            },
        ];
        let mut input = duct();
        let left = ResolvedFace { name: "left end".into(), ..end(0.0, -1.0) };
        let right = ResolvedFace { name: "right end".into(), ..end(100.0, 1.0) };
        let inner =
            ResolvedFace { name: "inner y".into(), plane: Some(([0.0, 1.0, 0.0], [50.0, 4.0, 20.0])), body: "Duct".into(), triangles: vec![12, 13] };
        // In face_refs order: the condition's face, then the goals' faces.
        input.faces = vec![left.clone(), left, right, inner];
        let built = build(&channel(), &goals, &input).expect("goals");
        assert!(built.issues.is_empty(), "{:?}", built.issues);
        let rows = built.case["goals"].as_array().cloned().unwrap_or_default();
        let target = |name: &str| rows.iter().find(|g| g["name"] == name).map(|g| g["target"].clone());
        assert_eq!(target("inflow"), Some(serde_json::json!("left")));
        assert_eq!(target("outflow"), Some(serde_json::json!("right")));
        // A face without a condition gets a wall patch of its own to be measured on.
        assert_eq!(target("wall_p"), Some(serde_json::json!("goal_wall_p")));
        assert_eq!(built.case["face_conditions"][0]["name"], "goal_wall_p");
        let mid = rows.iter().find(|g| g["name"] == "mid_p").cloned().unwrap_or_default();
        assert_eq!(mid["point"], serde_json::json!([0.05, 0.02, 0.02]));
        assert_eq!(rows.iter().find(|g| g["name"] == "balance").map(|g| g["expression"].clone()), Some(serde_json::json!("inflow + outflow")));
        // Goal names must be identifiers an equation can use.
        let mut bad_name = goals.clone();
        bad_name.goals[0].name = "in flow".into();
        assert!(bad_name.validate().is_err());
    }

    #[test]
    fn cut_cells_are_asked_for_and_keep_thin_walls_without_the_staircase_note() {
        // 10 mm cells and 4 mm walls: the staircase loses them, cut cells keep them.
        let mut coarse = link(&[]);
        coarse.cell_mm = Some(10.0);
        coarse.cut_cells = true;
        let built = build(&channel(), &coarse, &duct()).expect("cut");
        assert_eq!(built.case["cut_cells"], serde_json::json!(true));
        // The 4 mm walls read as 3.75 mm on the 1.25 mm sub-lattice: 100 × (40² − 32.5²) of the
        // duct's 100 × (40² − 32²) mm³.
        let expected = (40.0f64 * 40.0 - 32.5 * 32.5) / (40.0 * 40.0 - 32.0 * 32.0);
        assert!((built.solid_fraction - expected).abs() < 1e-9, "{} against {expected}", built.solid_fraction);
        assert!(!built.notes.iter().any(|n| n.contains("walls thinner than a cell")), "{:?}", built.notes);
        // The duct's walls cover the boundary layer of the four side faces of the domain.
        let velocity = |case: &Value, patch: &str| case["velocity_boundaries"].as_array().into_iter().flatten().any(|r| r["patch"] == patch);
        for patch in ["bottom", "top", "front", "back"] {
            assert!(!velocity(&built.case, patch), "{patch}");
        }
        assert!(velocity(&built.case, "left"), "the open ends have fluid at them");
        // Off again: the flag goes.
        coarse.cut_cells = false;
        assert_eq!(build(&built.case, &coarse, &duct()).expect("staircase").case.get("cut_cells"), None);
        // Not with heat transfer yet.
        coarse.cut_cells = true;
        coarse.heat = Some(ModelHeat::default());
        let why = build(&channel(), &coarse, &duct()).expect_err("thermal").to_string();
        assert!(why.contains("cut cells do not work with heat transfer"), "{why}");
    }

    #[test]
    fn coarse_cells_lose_thin_walls_and_say_so() {
        let mut coarse = link(&[]);
        coarse.cell_mm = Some(10.0);
        let built = build(&channel(), &coarse, &duct()).expect("build");
        assert!(built.notes.iter().any(|n| n.contains("walls thinner than a cell")), "{:?}", built.notes);
    }
}
