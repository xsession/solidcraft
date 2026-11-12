//! Typed contracts for the independent Multiphysics commands.
use crate::multiphysics_jobs::JobSnapshot;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use solidcraft_multiphysics::{Project, Study, Template, Validation, case_hash};

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Empty {}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Create {
    pub template: String,
    pub name: String,
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
#[derive(Serialize, JsonSchema)]
pub struct StudyReply {
    pub study: StudyView,
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
