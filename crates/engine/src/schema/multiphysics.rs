//! Typed contracts for the independent Multiphysics commands.
use crate::multiphysics_jobs::JobSnapshot;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use solidcraft_multiphysics::{Project, Study, Template, Validation, case_hash};

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SurfaceInput {
    pub mesh: solidcraft_multiphysics::mesh_diagnostics::SurfaceMesh,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Empty {}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Create {
    pub template: String,
    pub name: String,
    #[serde(default)]
    pub case: Option<Value>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct StudyId {
    pub id: String,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Update {
    pub id: String,
    pub expected_revision: u64,
    pub name: String,
    pub question: String,
    pub case: Value,
}
/// One study per value, replacing the existing scalar at an RFC 6901 pointer in a copy of
/// the source study's case.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Sweep {
    pub id: String,
    pub pointer: String,
    pub values: Vec<Value>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CompareStudies {
    #[serde(default)]
    pub ids: Vec<String>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Wizard {
    pub name: String,
    #[serde(default)]
    pub question: String,
    pub spec: solidcraft_multiphysics::wizard::EnclosureSpec,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct JobId {
    pub job_id: String,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExportCase {
    pub id: String,
    pub path: String,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExportResult {
    pub id: String,
    pub run_id: String,
    pub directory: String,
}
/// One retained original VTU or PVD artifact, addressed by saved study and run identity.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InspectResult {
    pub id: String,
    pub run_id: String,
    pub artifact: String,
}

#[derive(Serialize, JsonSchema)]
pub struct ResultInspection {
    pub id: String,
    pub run_id: String,
    pub artifact: String,
    pub captured_study_revision: u64,
    pub current_study_revision: u64,
    pub case_hash: String,
    pub stale: bool,
    pub content: ResultInspectionContent,
}

/// Metadata only: no field arrays, guessed units, or implied physical validation.
#[derive(Serialize, JsonSchema)]
#[serde(tag = "format", rename_all = "lowercase")]
pub enum ResultInspectionContent {
    Vtu { points: usize, cells: usize, connectivity_indices: usize, cell_types: Vec<u8>, fields: Vec<ResultField> },
    Pvd { frames: Vec<ResultFrame> },
}

#[derive(Serialize, JsonSchema)]
pub struct ResultField {
    pub name: String,
    pub association: String,
    pub data_type: String,
    pub components: usize,
    pub tuples: usize,
}

#[derive(Serialize, JsonSchema)]
pub struct ResultFrame {
    /// Producer coordinate, which may be physical time or iteration; no unit is inferred.
    pub coordinate: f64,
    pub file: String,
    pub part: usize,
    /// True when the referenced file is retained with this same run.
    pub retained: bool,
}
#[derive(Serialize, JsonSchema)]
pub struct StudyReply {
    pub study: StudyView,
}
/// A new study with a copy of another study's committed case.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CloneStudy {
    pub id: String,
    #[serde(default)]
    pub name: String,
}
#[derive(Serialize, JsonSchema)]
pub struct CloneReply {
    pub study: StudyView,
}
#[derive(Serialize, JsonSchema)]
pub struct SweepReply {
    pub studies: Vec<StudyView>,
}
#[derive(Serialize, JsonSchema)]
pub struct WizardReply {
    pub study: StudyView,
    pub report: solidcraft_multiphysics::wizard::WizardReport,
}
#[derive(Serialize, JsonSchema)]
pub struct ProjectReply {
    pub project: ProjectView,
}

/// One selected study's authoring state. Retained dataset bytes never enter command replies.
#[derive(Serialize, JsonSchema)]
pub struct StudyView {
    pub id: String,
    pub name: String,
    pub question: String,
    pub profile: String,
    pub revision: u64,
    pub case: Value,
    pub runs: Vec<RunSummary>,
    /// For a study of the loaded model: what it takes from the design.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub model: Option<solidcraft_multiphysics::model::ModelLink>,
}

/// A project listing contains identity and run metadata, without numerical case payloads.
#[derive(Serialize, JsonSchema)]
pub struct ProjectView {
    pub next_id: u64,
    pub studies: Vec<StudySummary>,
}

#[derive(Serialize, JsonSchema)]
pub struct StudySummary {
    pub id: String,
    pub name: String,
    pub question: String,
    pub profile: String,
    pub revision: u64,
    pub case_hash: String,
    /// The authored solver id when it is text of at most 80 bytes; otherwise absent.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub solver: Option<String>,
    pub runs: Vec<RunSummary>,
}

#[derive(Serialize, JsonSchema)]
pub struct RunSummary {
    pub id: String,
    pub study_revision: u64,
    pub case_hash: String,
    pub stale: bool,
    pub artifact_count: usize,
    pub execution_completed: bool,
    pub linear_solver_converged: Option<bool>,
    pub steady_state_converged: Option<bool>,
}

fn run_summaries(study: &Study, hash: &str) -> Vec<RunSummary> {
    study
        .runs
        .iter()
        .map(|run| RunSummary {
            id: run.id.clone(),
            study_revision: run.study_revision,
            case_hash: run.case_hash.clone(),
            stale: run.study_revision != study.revision || run.case_hash != hash,
            artifact_count: run.artifacts.len(),
            execution_completed: run.manifest.get("execution_completed").and_then(Value::as_bool).unwrap_or(false),
            linear_solver_converged: run.manifest.get("linear_solver_converged").and_then(Value::as_bool),
            steady_state_converged: run.manifest.get("steady_state_converged").and_then(Value::as_bool),
        })
        .collect()
}

impl TryFrom<&Study> for StudyView {
    type Error = solidcraft_multiphysics::Error;
    fn try_from(study: &Study) -> Result<Self, Self::Error> {
        let hash = case_hash(&study.case)?;
        Ok(Self {
            id: study.id.clone(),
            name: study.name.clone(),
            question: study.question.clone(),
            profile: study.profile.clone(),
            revision: study.revision,
            case: study.case.clone(),
            runs: run_summaries(study, &hash),
            model: study.model.clone(),
        })
    }
}

impl TryFrom<&Study> for StudySummary {
    type Error = solidcraft_multiphysics::Error;
    fn try_from(study: &Study) -> Result<Self, Self::Error> {
        let hash = case_hash(&study.case)?;
        let solver = study.case.get("solver").and_then(Value::as_str).filter(|solver| solver.len() <= 80).map(str::to_owned);
        Ok(Self {
            id: study.id.clone(),
            name: study.name.clone(),
            question: study.question.clone(),
            profile: study.profile.clone(),
            revision: study.revision,
            solver,
            runs: run_summaries(study, &hash),
            case_hash: hash,
        })
    }
}

impl TryFrom<&Project> for ProjectView {
    type Error = solidcraft_multiphysics::Error;
    fn try_from(project: &Project) -> Result<Self, Self::Error> {
        Ok(Self { next_id: project.next_id, studies: project.studies.iter().map(StudySummary::try_from).collect::<Result<_, _>>()? })
    }
}
#[derive(Serialize, JsonSchema)]
pub struct CatalogReply {
    pub catalog: Value,
    pub templates: Vec<Template>,
    pub host_available: bool,
}
#[derive(Serialize, JsonSchema)]
pub struct PreflightReply {
    pub validation: Validation,
    pub host_available: bool,
    pub execution_scope: String,
}
#[derive(Serialize, JsonSchema)]
pub struct JobsReply {
    pub jobs: Vec<JobSnapshot>,
}
#[derive(Serialize, JsonSchema)]
pub struct Removed {
    pub removed: String,
}
#[derive(Serialize, JsonSchema)]
pub struct Exported {
    pub files: Vec<String>,
}

/// CAD bodies to bring into a 3-D airflow case (see `case_scene::import_components`).
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ComponentImport {
    pub case: Value,
    /// Body names of the current model; empty means every body.
    #[serde(default)]
    pub bodies: Vec<String>,
    /// Tessellation quality: standard (default), fine or finest.
    #[serde(default)]
    pub quality: Option<String>,
    #[serde(default)]
    pub options: solidcraft_multiphysics::case_scene::ImportOptions,
}

/// The bodies of the current model, by name.
#[derive(Serialize, JsonSchema)]
pub struct CadBodies {
    pub bodies: Vec<String>,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FeaturePrepare {
    pub case: Value,
    pub kind: solidcraft_multiphysics::authoring::FeatureKind,
    pub name: String,
    #[serde(default)]
    pub replace_fan_patch_conditions: bool,
}

// ---------------------------------------------------------------- studies of the loaded model

/// `multiphysics.model.create`: a study of the loaded model.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelCreate {
    pub name: String,
    /// The 3-D flow template whose fluid and solver settings the study starts from
    /// (default `flow-channel`).
    #[serde(default)]
    pub template: Option<String>,
    #[serde(default)]
    pub analysis: solidcraft_multiphysics::model::Analysis,
    /// Bodies of the design; empty or absent means every body.
    #[serde(default)]
    pub bodies: Vec<String>,
    #[serde(default)]
    pub margin_mm: f64,
    #[serde(default)]
    pub cell_mm: Option<f64>,
}

/// `multiphysics.model.set`: replace a model study's settings and build its case again.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelSet {
    pub id: String,
    pub expected_revision: u64,
    #[serde(default)]
    pub analysis: Option<solidcraft_multiphysics::model::Analysis>,
    #[serde(default)]
    pub bodies: Option<Vec<String>>,
    #[serde(default)]
    pub margin_mm: Option<f64>,
    /// A cell size in mm; `0` goes back to the automatic size.
    #[serde(default)]
    pub cell_mm: Option<f64>,
    /// Keep the fluid part of the cells the bodies cut (companion 0.14.0).
    #[serde(default)]
    pub cut_cells: Option<bool>,
}

/// A face of the model, by its name or by a point on it (mm).
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelFacePick {
    pub body: String,
    #[serde(default)]
    pub face: Option<String>,
    #[serde(default)]
    pub point: Option<[f64; 3]>,
}

/// `multiphysics.model.condition`: put a condition on a face of the model, or remove one.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelCondition {
    pub id: String,
    pub expected_revision: u64,
    pub name: String,
    #[serde(default)]
    pub face: Option<ModelFacePick>,
    #[serde(default)]
    pub role: Option<solidcraft_multiphysics::model::OpeningRole>,
    /// Remove the condition called `name` instead.
    #[serde(default)]
    pub remove: bool,
}

/// What building a model study's case gave.
#[derive(Serialize, JsonSchema)]
pub struct ModelReply {
    pub study: StudyView,
    pub link: solidcraft_multiphysics::model::ModelLink,
    /// The CAD point at the domain origin (mm).
    pub origin_mm: [f64; 3],
    /// Which domain face each condition lies on: `[condition, face]`.
    pub placed: Vec<(String, String)>,
    pub solid_fraction: f64,
    pub fluid_regions: usize,
    pub notes: Vec<String>,
    pub issues: Vec<String>,
}

/// `multiphysics.model.status`: whether the case still matches the model.
#[derive(Serialize, JsonSchema)]
pub struct ModelStatus {
    pub id: String,
    /// The case was built from the model as it is now.
    pub current: bool,
    /// Why the case cannot be built now (a body or face is gone), if so.
    pub problem: Option<String>,
    pub link: solidcraft_multiphysics::model::ModelLink,
    /// The computational domain the case is built in, in CAD coordinates (mm): the bodies'
    /// exact box, plus the margin for an external study.
    pub domain_mm: Option<[[f64; 3]; 2]>,
    /// Where each condition's face is, for drawing it: in the order of `link.openings`.
    pub faces: Vec<ModelFaceMark>,
}

/// A point on a condition's face and the face's outward normal there (exact geometry).
#[derive(Serialize, JsonSchema)]
pub struct ModelFaceMark {
    pub name: String,
    pub point_mm: Option<[f64; 3]>,
    pub normal: Option<[f64; 3]>,
}

/// `multiphysics.model.leak`: a path through the fluid between two points of the model (mm).
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelLeak {
    pub id: String,
    /// A face, by a point on it (mm), on the side the path starts.
    pub from: ModelFacePick,
    pub to: ModelFacePick,
}

/// The fluid path between the two points, if fluid connects them.
#[derive(Serialize, JsonSchema)]
pub struct ModelLeakReply {
    pub connected: bool,
    /// Fluid cell centres along the shortest path, in CAD coordinates (mm).
    pub path_mm: Vec<[f64; 3]>,
}

/// `multiphysics.model.refine`: smaller cells around a body or a face of the model, or remove.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelRefine {
    pub id: String,
    pub expected_revision: u64,
    pub name: String,
    #[serde(default)]
    pub body: Option<String>,
    #[serde(default)]
    pub face: Option<ModelFacePick>,
    /// Split cells by 2 to 4 (default 2).
    #[serde(default)]
    pub factor: Option<usize>,
    /// Cells of space around the target that are split too (default 1).
    #[serde(default)]
    pub margin_cells: Option<f64>,
    #[serde(default)]
    pub remove: bool,
}

/// `multiphysics.model.goal`: a monitored quantity on the model, or remove one.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelGoalInput {
    pub id: String,
    pub expected_revision: u64,
    pub name: String,
    #[serde(default)]
    pub kind: Option<solidcraft_multiphysics::model::ModelGoalKind>,
    #[serde(default)]
    pub face: Option<ModelFacePick>,
    #[serde(default)]
    pub point_mm: Option<[f64; 3]>,
    #[serde(default)]
    pub expression: Option<String>,
    #[serde(default)]
    pub remove: bool,
}
