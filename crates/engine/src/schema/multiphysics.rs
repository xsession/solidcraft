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
/// `multiphysics.result.report`: one study's component table over several runs.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReportRuns {
    pub id: String,
    /// The runs, in this order (default: every run with a thermal report, oldest first).
    #[serde(default)]
    pub run_ids: Vec<String>,
    /// A new directory to write report.csv and report.html into.
    #[serde(default)]
    pub directory: Option<String>,
    /// Also return the CSV and HTML texts in the reply.
    #[serde(default)]
    pub inline: bool,
}

#[derive(Serialize, JsonSchema)]
pub struct ReportReply {
    pub runs: Vec<solidcraft_multiphysics::report::RunSummary>,
    /// The files written into `directory`.
    pub files: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub csv: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub html: Option<String>,
}

/// One run to compare: a study and one of its runs (default: its latest run with a thermal
/// report).
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CompareSide {
    pub id: String,
    #[serde(default)]
    pub run_id: Option<String>,
}

/// `multiphysics.result.compare`: two runs' component tables part by part.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CompareRuns {
    pub first: CompareSide,
    pub second: CompareSide,
    /// Also return the comparison as CSV.
    #[serde(default)]
    pub csv: bool,
}

#[derive(Serialize, JsonSchema)]
pub struct CompareReply {
    pub comparison: solidcraft_multiphysics::variants::Comparison,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub csv: Option<String>,
}

/// `multiphysics.model.mesh_study` makes the two coarser grids of a mesh study.
#[derive(Serialize, JsonSchema)]
pub struct MeshStudyReply {
    /// The fine study and its two clones, coarser by √2 and by 2: fine, medium, coarse.
    pub studies: [String; 3],
    /// Their base cell sizes (mm) and cell counts.
    pub cell_mm: [f64; 3],
    pub cells: [f64; 3],
}

/// `multiphysics.result.convergence`: the same study on three grids.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Convergence {
    pub fine: CompareSide,
    pub medium: CompareSide,
    pub coarse: CompareSide,
}

#[derive(Serialize, JsonSchema)]
pub struct ConvergenceReply {
    /// Cells of the three runs' grids, fine to coarse.
    pub cells: [f64; 3],
    /// Ratios of the grids' representative sizes, fine to medium and medium to coarse.
    pub ratios: [f64; 2],
    /// The inflow the temperature rises are measured from, °C.
    pub inflow_c: f64,
    pub parts: Vec<solidcraft_multiphysics::convergence::PartEstimate>,
    /// Notes on the grids (a ratio under 1.3, a run of an earlier revision).
    pub notes: Vec<String>,
}

/// `multiphysics.view.save`: keep a result view in a study, by name.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ViewSave {
    pub id: String,
    /// The view: its name, and what it shows. Without `ghosted`, the parts ghosted now.
    pub view: solidcraft_multiphysics::result_views::ResultView,
    /// Take the parts ghosted now in place of the view's `ghosted`.
    #[serde(default)]
    pub ghosted_now: bool,
}

/// `multiphysics.view.apply`, `multiphysics.view.remove`: a study's view by name.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ViewName {
    pub id: String,
    pub name: String,
}

#[derive(Serialize, JsonSchema)]
pub struct ViewsReply {
    /// The study's views, in the order they were first saved.
    pub views: Vec<solidcraft_multiphysics::result_views::ResultView>,
    /// For a save, whether it replaced a view of the same name.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub replaced: Option<bool>,
}

#[derive(Serialize, JsonSchema)]
pub struct ViewApplied {
    /// The view, for the host to set its cut plot, colour scale and range from.
    pub view: solidcraft_multiphysics::result_views::ResultView,
    /// The parts drawn see-through now.
    pub ghosted: Vec<String>,
    /// Parts the view names that the design no longer has (skipped).
    pub missing: Vec<String>,
    /// What colour-by gave, when the view colours the parts and the study has what it needs;
    /// otherwise why not, in `notes`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub colouring: Option<solidcraft_multiphysics::setup_check::Colouring>,
    pub notes: Vec<String>,
}

/// `multiphysics.module.capture`: keep a study's setup of one placement of a component.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModuleCapture {
    /// The study of the model whose setup is kept.
    pub id: String,
    /// The component definition, by name.
    pub component: String,
    /// The module's name (default "<component> setup"); one of the same name is replaced.
    #[serde(default)]
    pub name: Option<String>,
    /// Which placement to read, by its occurrences' names joined with `/` (default the first).
    #[serde(default)]
    pub placement: Option<String>,
}

/// `multiphysics.module.apply`: give a module to every placement of its component.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModuleApply {
    pub id: String,
    pub expected_revision: u64,
    pub module: String,
}

/// `multiphysics.module.remove`: a module by name.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModuleName {
    pub module: String,
}

/// `multiphysics.module.export` (the module to write, a new file) and `.import` (a file; module
/// renames what it holds).
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModuleFile {
    #[serde(default)]
    pub module: Option<String>,
    pub path: String,
}

#[derive(Serialize, JsonSchema)]
pub struct ModuleReply {
    pub module: solidcraft_multiphysics::modules::Module,
    /// Whether a module of the same name was replaced.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub replaced: Option<bool>,
    /// The component's placements, by label (for a capture).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub placements: Vec<String>,
}

#[derive(Serialize, JsonSchema)]
pub struct ModulesReply {
    pub modules: Vec<solidcraft_multiphysics::modules::Module>,
}

#[derive(Serialize, JsonSchema)]
pub struct ModuleApplied {
    /// Placements set up, and bodies set in them.
    pub placements: usize,
    pub bodies: usize,
    /// `placement/body` of the module's bodies a placement does not have.
    pub missing: Vec<String>,
    /// The study as built again (as `multiphysics.model.set` replies).
    pub build: Value,
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
    /// A part selection whose bodies the study uses, worked out again at every build
    /// (`multiphysics.parts.selection`); not with `bodies`.
    #[serde(default)]
    pub selection: Option<String>,
    #[serde(default)]
    pub margin_mm: f64,
    #[serde(default)]
    pub cell_mm: Option<f64>,
    /// The fluid and run length (default air for two seconds, stopping once the goals settle).
    #[serde(default)]
    pub flow: Option<solidcraft_multiphysics::model::ModelFlow>,
}

/// `multiphysics.model.leave_out_small`: leave the parts smaller than a size out of a study.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelLeaveOutSmall {
    pub id: String,
    pub expected_revision: u64,
    /// Parts whose largest bounding-box side is under this many millimetres are left out.
    pub below_mm: f64,
}

/// `multiphysics.model.set`: replace a model study's settings and build its case again.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelSet {
    pub id: String,
    pub expected_revision: u64,
    #[serde(default)]
    pub analysis: Option<solidcraft_multiphysics::model::Analysis>,
    /// The study's bodies by name; replaces a part selection.
    #[serde(default)]
    pub bodies: Option<Vec<String>>,
    /// A part selection whose bodies the study uses from now on; `""` goes back to the list.
    #[serde(default)]
    pub selection: Option<String>,
    #[serde(default)]
    pub margin_mm: Option<f64>,
    /// A cell size in mm; `0` goes back to the automatic size.
    #[serde(default)]
    pub cell_mm: Option<f64>,
    /// Keep the fluid part of the cells the bodies cut (companion 0.14.0).
    #[serde(default)]
    pub cut_cells: Option<bool>,
    /// How fine the automatic grid is, 1 to 7 (3 is the automatic size); 0 goes back to it.
    /// Setting a level clears a typed cell size.
    #[serde(default)]
    pub mesh_level: Option<u8>,
    /// The fluid and how long the run simulates.
    #[serde(default)]
    pub flow: Option<solidcraft_multiphysics::model::ModelFlow>,
    /// The bodies that are porous zones, with their resistances (replaces the list; [] for none).
    #[serde(default)]
    pub porous: Option<Vec<solidcraft_multiphysics::model::ModelPorous>>,
    /// Bodies that stand for fans inside the domain (replaces the list; [] for none).
    #[serde(default)]
    pub internal_fans: Option<Vec<solidcraft_multiphysics::model::ModelInternalFan>>,
    /// Regions turning about an axis around bodies (replaces the list; [] for none).
    #[serde(default)]
    pub rotating: Option<Vec<solidcraft_multiphysics::model::ModelRotating>>,
}

/// `multiphysics.model.fan`: a fan on a face of the model on the domain boundary, or remove it.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelFanInput {
    pub id: String,
    pub expected_revision: u64,
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub face: Option<ModelFacePick>,
    /// `[volume flow m³/s, pressure rise Pa]` points, flow rising and pressure falling.
    #[serde(default)]
    pub curve: Option<Vec<[f64; 2]>>,
    /// The curve's pressure is total, not static.
    #[serde(default)]
    pub total_pressure: bool,
    #[serde(default)]
    pub remove: bool,
}

/// `multiphysics.model.adapt`: refine where a kept run's field changes most from cell to cell.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelAdapt {
    pub id: String,
    pub expected_revision: u64,
    pub run_id: String,
    /// The run's field file (default its last flow_*.vtu).
    #[serde(default)]
    pub artifact: Option<String>,
    /// The field (default the velocity).
    #[serde(default)]
    pub field: Option<String>,
    /// The share of cells with the largest jumps to refine, above 0 and at most 0.5 (default 0.1).
    #[serde(default)]
    pub fraction: Option<f64>,
    /// Most regions, 1 to 16 (default 4).
    #[serde(default)]
    pub regions: Option<usize>,
    /// Split cells by 2 to 4 (default 2).
    #[serde(default)]
    pub factor: Option<usize>,
    /// Remove the adaptive refinements instead.
    #[serde(default)]
    pub remove: bool,
}

/// What a solution-adaptive refinement found.
#[derive(Serialize, JsonSchema)]
pub struct ModelAdapted {
    pub run_id: String,
    pub artifact: String,
    pub field: String,
    pub cells: usize,
    pub flagged: usize,
    /// The refined boxes in CAD coordinates (mm), as `[low, high]`.
    pub regions_mm: Vec<[[f64; 3]; 2]>,
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
    /// What `multiphysics.model.adapt` found, for its reply.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub adapted: Option<ModelAdapted>,
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
    /// The body of a body temperature goal.
    #[serde(default)]
    pub body: Option<String>,
    #[serde(default)]
    pub remove: bool,
}

/// `multiphysics.model.heat`: heat transfer in a study of the model, on or off, and the fluid's
/// thermal settings.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelHeatSet {
    pub id: String,
    pub expected_revision: u64,
    #[serde(default)]
    pub enabled: Option<bool>,
    /// The temperature of the fluid flowing in, K.
    #[serde(default)]
    pub inflow_temperature_k: Option<f64>,
    #[serde(default)]
    pub fluid: Option<solidcraft_multiphysics::model::FluidHeat>,
    /// Natural convection: true (air under the design's gravity, Z up), settings of its own,
    /// or false to solve heat on the settled flow again.
    #[serde(default)]
    pub buoyancy: Option<BuoyancyInput>,
}

/// Natural convection on (air's defaults), off, or given.
#[derive(Deserialize, JsonSchema)]
#[serde(untagged)]
pub enum BuoyancyInput {
    On(bool),
    Given(solidcraft_multiphysics::model::ModelBuoyancy),
}

/// `multiphysics.model.material`: the material and power of bodies of the model, or remove.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelMaterial {
    pub id: String,
    pub expected_revision: u64,
    /// The bodies, by name (1 to 256).
    #[serde(default)]
    pub bodies: Vec<String>,
    /// Or a part selection: the material goes on the bodies it holds now.
    #[serde(default)]
    pub selection: Option<String>,
    /// W/(m·K); kept when a body already has a material and this is absent.
    #[serde(default)]
    pub conductivity_w_per_m_k: Option<f64>,
    /// Heat released evenly in each body, W.
    #[serde(default)]
    pub power_w: Option<f64>,
    #[serde(default)]
    pub density_kg_per_m3: Option<f64>,
    #[serde(default)]
    pub specific_heat_j_per_kg_k: Option<f64>,
    #[serde(default)]
    pub max_temperature_k: Option<f64>,
    /// Remove the bodies' materials instead (they conduct like the starter material).
    #[serde(default)]
    pub remove: bool,
}

/// `multiphysics.model.heat_source`: heat released on a face of the model, or remove one.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelHeatSource {
    pub id: String,
    pub expected_revision: u64,
    pub name: String,
    #[serde(default)]
    pub face: Option<ModelFacePick>,
    #[serde(default)]
    pub power_w: Option<f64>,
    #[serde(default)]
    pub remove: bool,
}

/// `multiphysics.model.thermal_face`: a thermal condition on a face of the model on the domain
/// boundary, or remove one.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelThermalFace {
    pub id: String,
    pub expected_revision: u64,
    pub name: String,
    #[serde(default)]
    pub face: Option<ModelFacePick>,
    #[serde(default)]
    pub condition: Option<solidcraft_multiphysics::model::ThermalFace>,
    #[serde(default)]
    pub remove: bool,
}

/// `multiphysics.parts.tag`: tags and attributes on some of the design's parts.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsTag {
    /// The bodies, by name; or `selection`.
    #[serde(default)]
    pub bodies: Vec<String>,
    /// A part selection whose bodies, as the design is now, are changed.
    #[serde(default)]
    pub selection: Option<String>,
    #[serde(default)]
    pub add_tags: Vec<String>,
    #[serde(default)]
    pub remove_tags: Vec<String>,
    /// Attributes to set (a number or a text), or `null` to remove one.
    #[serde(default)]
    pub attributes: std::collections::BTreeMap<String, Option<solidcraft_multiphysics::parts::AttrValue>>,
}

/// `multiphysics.parts.selection`: add, replace or remove a named part selection.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsSelection {
    pub name: String,
    /// What it holds: saved bodies or a live query. Absent with `remove`.
    #[serde(default)]
    pub members: Option<solidcraft_multiphysics::parts::Members>,
    #[serde(default)]
    pub remove: bool,
}

/// One body of the design with what the registry records about it.
#[derive(Serialize, JsonSchema)]
pub struct PartView {
    pub body: String,
    pub tags: Vec<String>,
    pub attributes: std::collections::BTreeMap<String, solidcraft_multiphysics::parts::AttrValue>,
}

/// A selection and the bodies it holds now.
#[derive(Serialize, JsonSchema)]
pub struct PartSelectionView {
    pub selection: solidcraft_multiphysics::parts::PartSelection,
    pub bodies: Vec<String>,
    /// Bodies a saved selection names that the design no longer has.
    pub missing: Vec<String>,
    /// Why the selection cannot be worked out now, if so.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub problem: Option<String>,
}

/// `multiphysics.parts.list` and the reply of the commands that change parts.
#[derive(Serialize, JsonSchema)]
pub struct PartsReply {
    /// Every body of the design, in the design's order.
    pub parts: Vec<PartView>,
    pub selections: Vec<PartSelectionView>,
    /// Parts the registry records that the design no longer has.
    pub missing: Vec<String>,
}

/// `multiphysics.parts.table`: the component table of a study.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsTable {
    /// The study.
    pub id: String,
    /// The run the results come from; by default the latest run with a thermal report.
    #[serde(default)]
    pub run_id: Option<String>,
    /// Only the bodies of this part selection.
    #[serde(default)]
    pub selection: Option<String>,
    /// Show the design's bodies the study leaves out too.
    #[serde(default)]
    pub all_bodies: bool,
    /// Add the table as CSV to the reply.
    #[serde(default)]
    pub csv: bool,
}

#[derive(Serialize, JsonSchema)]
pub struct PartsTableReply {
    pub table: solidcraft_multiphysics::component_table::ComponentTable,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub csv: Option<String>,
}

/// `multiphysics.parts.select`: put some of the design's parts in its selection, the one the
/// tree, the component table and the 3-D view share.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsSelect {
    /// The bodies, by name; or `selection`. None (and no selection): the selection is cleared.
    #[serde(default)]
    pub bodies: Vec<String>,
    /// A part selection whose bodies, as the design is now, are selected.
    #[serde(default)]
    pub selection: Option<String>,
    /// Add to what is selected instead of replacing it.
    #[serde(default)]
    pub add: bool,
}

#[derive(Serialize, JsonSchema)]
pub struct PartsSelectReply {
    /// The bodies selected now, in the design's order.
    pub selected: Vec<String>,
    /// Bodies a saved selection names that the design no longer has.
    pub missing: Vec<String>,
}

/// `multiphysics.parts.ghost`: draw some parts see-through, letting clicks reach what is behind.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsGhost {
    /// The bodies, by name; or `selection`.
    #[serde(default)]
    pub bodies: Vec<String>,
    /// A part selection whose bodies, as the design is now, are changed.
    #[serde(default)]
    pub selection: Option<String>,
    /// False: draw them solid again.
    #[serde(default = "yes")]
    pub ghost: bool,
    /// Ghost only these, drawing every other part solid (with `ghost`); or, with no bodies or
    /// selection, draw every part solid.
    #[serde(default)]
    pub only: bool,
}

fn yes() -> bool {
    true
}

/// `multiphysics.parts.find`: the study's parts that meet every rule, on setup and results.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsFind {
    pub id: String,
    pub rules: Vec<solidcraft_multiphysics::find::FindRule>,
    /// The run results come from (default: the latest with a thermal report).
    #[serde(default)]
    pub run_id: Option<String>,
    /// Also look at the design's bodies the study leaves out.
    #[serde(default)]
    pub all_bodies: bool,
    /// Put what is found in the design's selection.
    #[serde(default)]
    pub select: bool,
}

#[derive(Serialize, JsonSchema)]
pub struct PartsFindReply {
    pub found: Vec<String>,
    /// The run the result rules read; none before a run with a thermal report.
    pub run: Option<String>,
}

/// `multiphysics.parts.group`: tag parts by the kind of part their reference designators name.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsGroup {
    /// The bodies, by name; or `selection`; neither: every body of the design.
    #[serde(default)]
    pub bodies: Vec<String>,
    /// A part selection whose bodies, as the design is now, are grouped.
    #[serde(default)]
    pub selection: Option<String>,
}

#[derive(Serialize, JsonSchema)]
pub struct PartsGroupReply {
    /// The parts tagged, each with its designator and the kind of part (its new tag).
    pub grouped: Vec<solidcraft_multiphysics::grouping::Grouped>,
    /// The parts whose names carry no designator of a known kind, left as they were.
    pub other: Vec<String>,
}

/// `multiphysics.parts.similar`: the parts shaped like one part.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsSimilar {
    /// The part the others are compared with.
    pub body: String,
    /// How far, in per cent of the larger, a volume or a side may differ (default 1, at most 50).
    #[serde(default)]
    pub tolerance_pct: Option<f64>,
    /// Put the similar parts in the design's selection.
    #[serde(default)]
    pub select: bool,
}

#[derive(Serialize, JsonSchema)]
pub struct PartsSimilarReply {
    /// The parts alike in volume and bounding box, `body` among them, in the design's order.
    pub similar: Vec<String>,
}

/// `multiphysics.parts.check`: a study of the model checked before it runs.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsCheck {
    pub id: String,
}

/// `multiphysics.parts.colour`: colour the design's parts by their setup in a study.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartsColour {
    /// The study; not needed with `by: none`.
    #[serde(default)]
    pub id: Option<String>,
    pub by: solidcraft_multiphysics::setup_check::ColourBy,
    /// For verdicts: the run (default the latest with a thermal report).
    #[serde(default)]
    pub run_id: Option<String>,
}

#[derive(Serialize, JsonSchema)]
pub struct PartsGhostReply {
    /// The parts drawn see-through now, in the design's order.
    pub ghosted: Vec<String>,
}
