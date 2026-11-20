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
    /// Space around the bodies on every side, for an external analysis.
    #[serde(default)]
    pub margin_mm: f64,
    /// Cell size; by default chosen from the bodies (see `case_scene::import_components`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell_mm: Option<f64>,
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
    /// Every face the link names, in one order: the conditions', the refinements' and the
    /// goals'. [`ModelInput::faces`] resolves them in this order.
    pub fn face_refs(&self) -> Vec<&FaceRef> {
        let openings = self.openings.iter().map(|o| &o.face);
        let refined = self.refinements.iter().filter_map(|r| match &r.target {
            RefineTarget::Face { face } => Some(face),
            RefineTarget::Body { .. } => None,
        });
        openings.chain(refined).chain(self.goals.iter().filter_map(|g| g.face.as_ref())).collect()
    }

    /// Bounded structure only; whether the bodies and faces exist is decided when building.
    pub fn validate(&self) -> Result<()> {
        if self.bodies.len() > MAX_BODIES {
            return Err(Error::Limit(format!("a model study names at most {MAX_BODIES} bodies")));
        }
        for body in &self.bodies {
            text(body, "a body name")?;
        }
        if !(self.margin_mm.is_finite() && (0.0..=1e6).contains(&self.margin_mm)) {
            return Err(Error::Invalid("the margin must be 0 to 1e6 mm".into()));
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
            };
            if !ok {
                return Err(Error::Invalid(format!(
                    "goal {}: a face goal needs a face, a point goal a finite point, an equation its expression",
                    g.name
                )));
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
    if case.get("thermal").is_some_and(|t| !t.is_null()) {
        return Err(bad(format!(
            "{}: conditions on faces inside the domain do not work with heat transfer yet; put it on a face on the domain boundary",
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
    let (refine_faces, goal_faces) = rest.split_at(refine_count);
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
        ..Default::default()
    };
    // Conditions on body faces belong to the last build: drop them and their patches' rows,
    // so a removed condition leaves nothing behind. The domain faces' rows stay.
    let mut base = base.clone();
    if let Some(object) = base.as_object_mut() {
        object.remove("face_conditions");
        // Goals of a study of the model are the link's.
        object.remove("goals");
        for key in ["velocity_boundaries", "pressure_boundaries"] {
            if let Some(rows) = object.get_mut(key).and_then(Value::as_array_mut) {
                rows.retain(|row| row["patch"].as_str().is_some_and(|patch| PATCHES.contains(&patch)));
            }
        }
    }
    let imported = case_scene::import_components(&base, &input.obj, &options)?;
    let mut case = imported.edit.case;
    let mut notes = imported.notes;
    // Local refinements, then a time step that fits the smallest cells.
    let mut refine_face = refine_faces.iter();
    for r in &link.refinements {
        let edit = match &r.target {
            RefineTarget::Body { body } => case_scene::refine_around(&case, &[case_scene::part_name(body)], r.factor, r.margin_cells),
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
        .map_err(|e| bad(format!("{}: {e}", r.name)))?;
        case = edit.case;
    }
    if !link.refinements.is_empty() {
        case_scene::fit_time_step(&mut case, &mut notes);
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
        // The same face in a case with heat transfer is refused, not left uncoupled.
        let mut thermal = channel();
        thermal["thermal"] = serde_json::json!({"solids": []});
        let why = build(&thermal, &link(&[("Side vent", outlet)]), &input).expect_err("thermal").to_string();
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

    #[test]
    fn goals_read_condition_patches_domain_sides_and_wall_patches_made_for_them() {
        let inlet = OpeningRole::Inlet { speed_m_per_s: 0.5 };
        let mut goals = link(&[("in", inlet)]);
        let face = |name: &str| FaceRef { body: "Duct".into(), face: name.into() };
        goals.openings[0].face = face("left end");
        goals.goals = vec![
            ModelGoal { name: "inflow".into(), kind: ModelGoalKind::FaceVolumeFlow, face: Some(face("left end")), point_mm: None, expression: None },
            ModelGoal {
                name: "outflow".into(),
                kind: ModelGoalKind::FaceVolumeFlow,
                face: Some(face("right end")),
                point_mm: None,
                expression: None,
            },
            ModelGoal { name: "wall_p".into(), kind: ModelGoalKind::FaceMeanPressure, face: Some(face("inner y")), point_mm: None, expression: None },
            ModelGoal { name: "mid_p".into(), kind: ModelGoalKind::PointPressure, face: None, point_mm: Some([50.0, 20.0, 20.0]), expression: None },
            ModelGoal {
                name: "balance".into(),
                kind: ModelGoalKind::Equation,
                face: None,
                point_mm: None,
                expression: Some("inflow + outflow".into()),
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
        let mut thermal = channel();
        thermal["thermal"] = serde_json::json!({"solids": []});
        coarse.cut_cells = true;
        let why = build(&thermal, &coarse, &duct()).expect_err("thermal").to_string();
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
