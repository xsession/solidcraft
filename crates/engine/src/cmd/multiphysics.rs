//! Independent study authoring and asynchronous numerical work. No exact CAD geometry changes.
use super::CommandSpec;
use crate::schema::multiphysics::*;
use crate::schema::{TypedContract, output_schema_for};
use crate::{EngineError, Result, Session, multiphysics_jobs::JobSnapshot, params::bad};
use serde::{Serialize, de::DeserializeOwned};
use serde_json::Value;
use solidcraft_multiphysics::{Study, catalog as solver_catalog, preflight, templates};
use std::sync::Arc;

pub static COMMANDS: &[CommandSpec] = &[
    CommandSpec::new("multiphysics.catalog", "Multiphysics Catalogue", crate::checked!("multiphysics.catalog", catalog))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<CatalogReply>)
        .params("No parameters. Original 34-profile inventory, 23 bounded Rust case templates and separate host status. Inventory is not execution qualification."),
    CommandSpec::new("multiphysics.study.list", "List Multiphysics Studies", crate::checked!("multiphysics.study.list", list))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<ProjectReply>)
        .params("No parameters. Document-owned study and run summaries with case hashes; case payloads, manifests and artifact bytes are omitted."),
    CommandSpec::new("multiphysics.study.create", "Create Multiphysics Study", crate::checked!("multiphysics.study.create", create))
        .doc_only().icon("new").contract(&TypedContract::<Create>::new()).output(output_schema_for::<StudyReply>)
        .params("template: known case id; name: nonblank name. Creates one saved study and undo step without CAD-to-mesh conversion or execution."),
    CommandSpec::new("multiphysics.study.update", "Apply Multiphysics Study", crate::checked!("multiphysics.study.update", update))
        .doc_only().icon("settings").contract(&TypedContract::<Update>::new()).output(output_schema_for::<StudyReply>)
        .params("id, expected_revision, name, question, case: bounded JSON object. Atomic Apply; detailed solver schema and physical prerequisites are checked at execution."),
    CommandSpec::new("multiphysics.study.remove", "Remove Multiphysics Study", crate::checked!("multiphysics.study.remove", remove))
        .doc_only().icon("delete").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<Removed>)
        .params("id. Removes authored study in one undo step. External jobs are retained separately."),
    CommandSpec::new("multiphysics.preflight", "Check Multiphysics Study", crate::checked!("multiphysics.preflight", check))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<PreflightReply>)
        .params("id. Bounded structural check and separate worker status; not numerical or physical validation."),
    CommandSpec::new("multiphysics.run.start", "Run Multiphysics Case", crate::checked!("multiphysics.run.start", start))
        .doc_only().noundo().icon("play").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<JobSnapshot>)
        .params("id: committed study. Freezes revision and original input hash; starts injected asynchronous native companion. Refuses missing host and scratch sessions."),
    CommandSpec::new("multiphysics.run.status", "Multiphysics Run Status", crate::checked!("multiphysics.run.status", status))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<JobsReply>)
        .params("No parameters. Polls up to eight recent jobs; dataset bytes omitted to keep status replies bounded. Native UI reads shared admitted evidence."),
    CommandSpec::new("multiphysics.run.cancel", "Cancel Multiphysics Run", crate::checked!("multiphysics.run.cancel", cancel))
        .doc_only().noundo().icon("stop").contract(&TypedContract::<JobId>::new()).output(output_schema_for::<JobSnapshot>)
        .params("job_id. Cancelling persists until worker acknowledgement. Late completion after cancellation is discarded."),
    CommandSpec::new("multiphysics.run.keep", "Keep Multiphysics Result", crate::checked!("multiphysics.run.keep", keep))
        .doc_only().icon("save").contract(&TypedContract::<JobId>::new()).output(output_schema_for::<StudyReply>)
        .params("job_id: completed current job. Rechecks document, revision, original case and artifacts; retains portable evidence in one undo step."),
    CommandSpec::new("multiphysics.case.export", "Export Multiphysics Case", crate::checked!("multiphysics.case.export", export_case))
        .doc_only().noundo().icon("export").contract(&TypedContract::<ExportCase>::new()).output(output_schema_for::<Exported>)
        .params("id, path: new filename. Writes committed numerical case through atomic no-replace writer."),
    CommandSpec::new("multiphysics.result.export", "Export Multiphysics Results", crate::checked!("multiphysics.result.export", export_result))
        .doc_only().noundo().icon("export").contract(&TypedContract::<ExportResult>::new()).output(output_schema_for::<Exported>)
        .params("id, run_id, directory: new directory. Original VTU/PVD/CSV/JSON/S2P, case and manifest for ParaView; associations unchanged, no invented .pvsm or units."),
];

fn parse<T: DeserializeOwned>(params: &Value) -> Result<T> {
    serde_json::from_value(params.clone()).map_err(|e| bad("", e.to_string()))
}
fn reply(value: impl Serialize) -> Result<Value> {
    serde_json::to_value(value).map_err(|e| EngineError::Other(e.to_string()))
}
fn study_reply(study: &Study) -> Result<Value> {
    reply(StudyReply { study: StudyView::try_from(study).map_err(|e| bad("", e.to_string()))? })
}
fn selected(session: &Session, id: &str) -> Result<Study> {
    session.doc.multiphysics.study(id).cloned().ok_or_else(|| bad("", "Multiphysics study not found in this document"))
}
fn catalog(session: &mut Session, params: &Value) -> Result<Value> {
    let _: Empty = parse(params)?;
    reply(CatalogReply {
        catalog: solver_catalog().map_err(|e| bad("", e.to_string()))?,
        templates: templates().map_err(|e| bad("", e.to_string()))?,
        host_available: session.multiphysics_host_available(),
    })
}
fn list(session: &mut Session, params: &Value) -> Result<Value> {
    let _: Empty = parse(params)?;
    reply(ProjectReply { project: ProjectView::try_from(session.doc.multiphysics.as_ref()).map_err(|e| bad("", e.to_string()))? })
}
fn create(session: &mut Session, params: &Value) -> Result<Value> {
    let input: Create = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    // Undo must never allocate an old live identity to a different study.
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let study = project.create(&input.template, &input.name).map_err(|e| bad("", e.to_string()))?;
    session.multiphysics_jobs.next_study = project.next_id;
    session.doc_mut().multiphysics = Arc::new(project);
    study_reply(&study)
}
fn update(session: &mut Session, params: &Value) -> Result<Value> {
    let input: Update = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    let study = project.update(&input.id, input.expected_revision, &input.name, &input.question, input.case).map_err(|e| bad("", e.to_string()))?;
    if project != *session.doc.multiphysics {
        session.doc_mut().multiphysics = Arc::new(project);
    }
    study_reply(&study)
}
fn remove(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.remove(&input.id).map_err(|e| bad("", e.to_string()))?;
    session.doc_mut().multiphysics = Arc::new(project);
    reply(Removed { removed: input.id })
}
fn check(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    let study = selected(session, &input.id)?;
    reply(PreflightReply {
        validation: preflight(&study),
        host_available: session.multiphysics_host_available(),
        execution_scope: "Bounded cfd-solvers-rs 0.2.0 case, independent of CAD; full case schema checked at execution".into(),
    })
}
fn start(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    reply(session.start_multiphysics_run(&input.id)?)
}
fn status(session: &mut Session, params: &Value) -> Result<Value> {
    let _: Empty = parse(params)?;
    session.poll_multiphysics_jobs();
    let mut jobs = session.multiphysics_jobs();
    for job in &mut jobs {
        job.evidence = None;
    }
    reply(JobsReply { jobs })
}
fn cancel(session: &mut Session, params: &Value) -> Result<Value> {
    let input: JobId = parse(params)?;
    let mut job = session.cancel_multiphysics_run(&input.job_id)?;
    job.evidence = None;
    reply(job)
}
fn keep(session: &mut Session, params: &Value) -> Result<Value> {
    let input: JobId = parse(params)?;
    session.poll_multiphysics_jobs();
    let job =
        session.multiphysics_jobs().into_iter().find(|j| j.id == input.job_id).ok_or_else(|| bad("", "Numerical job not found in this document"))?;
    if job.stale {
        return Err(bad("", "The study changed; this historical result cannot become current evidence"));
    }
    let evidence = job.evidence.ok_or_else(|| bad("", "Wait for an admitted completed result"))?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    let study = project.attach_evidence(&job.study_id, evidence.as_ref().clone()).map_err(|e| bad("", e.to_string()))?;
    if project != *session.doc.multiphysics {
        session.doc_mut().multiphysics = Arc::new(project);
    }
    study_reply(&study)
}
fn export_case(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ExportCase = parse(params)?;
    let study = selected(session, &input.id)?;
    solidcraft_io::multiphysics::write_case(std::path::Path::new(&input.path), &study).map_err(|e| bad("", e.to_string()))?;
    reply(Exported { files: vec![input.path] })
}
fn export_result(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ExportResult = parse(params)?;
    let study = selected(session, &input.id)?;
    let run = study.runs.iter().find(|r| r.id == input.run_id).ok_or_else(|| bad("", "Select a retained result"))?;
    let files =
        solidcraft_io::multiphysics::write_results(std::path::Path::new(&input.directory), &study, run).map_err(|e| bad("", e.to_string()))?;
    reply(Exported { files })
}
