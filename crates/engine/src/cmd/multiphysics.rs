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
    CommandSpec::new("multiphysics.surface.inspect", "Inspect Multiphysics Surface", crate::checked!("multiphysics.surface.inspect", inspect_surface))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<SurfaceInput>::new()).output(output_schema_for::<solidcraft_multiphysics::mesh_diagnostics::SurfaceReport>)
        .params("mesh: bounded indexed triangle surface. Reports index connectivity only; does not establish watertight geometry or solver readiness."),
    CommandSpec::new("multiphysics.catalog", "Multiphysics Catalogue", crate::checked!("multiphysics.catalog", catalog))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<CatalogReply>)
        .params("No parameters. Original 34-profile inventory, 56 bounded Rust case templates and separate host status. Inventory is not execution qualification."),
    CommandSpec::new("multiphysics.study.list", "List Multiphysics Studies", crate::checked!("multiphysics.study.list", list))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<ProjectReply>)
        .params("No parameters. Document-owned study and run summaries with case hashes; case payloads, manifests and artifact bytes are omitted."),
    CommandSpec::new("multiphysics.study.create", "Create Multiphysics Study", crate::checked!("multiphysics.study.create", create))
        .doc_only().icon("new").contract(&TypedContract::<Create>::new()).output(output_schema_for::<StudyReply>)
        .params("template: known case id; name: nonblank name. Creates one saved study and undo step without CAD-to-mesh conversion or execution."),
    CommandSpec::new("multiphysics.study.update", "Apply Multiphysics Study", crate::checked!("multiphysics.study.update", update))
        .doc_only().icon("settings").contract(&TypedContract::<Update>::new()).output(output_schema_for::<StudyReply>)
        .params("id, expected_revision, name, question, case: bounded JSON object. Atomic Apply; detailed solver schema and physical prerequisites are checked at execution."),
    CommandSpec::new("multiphysics.study.sweep", "Sweep Multiphysics Study", crate::checked!("multiphysics.study.sweep", sweep))
        .doc_only().icon("new").contract(&TypedContract::<Sweep>::new()).output(output_schema_for::<SweepReply>)
        .params("id: source study; pointer: JSON pointer to an existing scalar in its case (not /solver); values: 1..16 values of the same type. Creates one checked study per value in one undo step; nothing runs."),
    CommandSpec::new("multiphysics.study.clone", "Clone Multiphysics Study", crate::checked!("multiphysics.study.clone", clone_study))
        .doc_only().icon("copy").contract(&TypedContract::<CloneStudy>::new()).output(output_schema_for::<CloneReply>)
        .params("id: source study; name: optional (default: \"<name> (copy)\"). Copies the committed case, question and setup into a new study in one undo step; retained runs stay with the source."),
    CommandSpec::new("multiphysics.wizard.enclosure", "Multiphysics Enclosure Wizard", crate::checked!("multiphysics.wizard.enclosure", wizard))
        .doc_only().icon("new").contract(&TypedContract::<Wizard>::new()).output(output_schema_for::<WizardReply>)
        .params("name, question, spec: enclosure size (mm), cooling, ambient (°C), wall convection, components (mm, W) and mesh level 1..5. Creates one study and reports every automatic choice; nothing runs."),
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
    CommandSpec::new("multiphysics.result.inspect", "Inspect Multiphysics Result", crate::checked!("multiphysics.result.inspect", inspect_result))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<InspectResult>::new()).output(output_schema_for::<ResultInspection>)
        .params("id, run_id, artifact: one retained ASCII VTU/PVD file. Returns bounded geometry, field or collection metadata, without dataset bytes or inferred units."),
    CommandSpec::new("multiphysics.result.query", "Query Multiphysics Result", crate::checked!("multiphysics.result.query", query_result))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<solidcraft_multiphysics::results::QueryRequest>::new()).output(output_schema_for::<solidcraft_multiphysics::results::QueryResponse>)
        .params("id, run_id, query: {kind: probe|min_max|line|integrate|compare|goals|components, ...}. Probes interpolate inside cells and report the nearest cell outside; line samples outside the mesh have no value; integrals are measure weighted; compare needs identical meshes. Every answer carries the run's case hash and a stale flag."),
    CommandSpec::new("multiphysics.study.compare", "Compare Multiphysics Studies", crate::checked!("multiphysics.study.compare", compare_studies))
        .doc_only().noundo().icon("list").contract(&TypedContract::<CompareStudies>::new()).output(output_schema_for::<solidcraft_multiphysics::results::MetricTable>)
        .params("ids: studies to tabulate (empty for all). Returns every manifest metric of each study's latest retained run, aligned by name, for sweep and what-if comparison."),
];

fn parse<T: DeserializeOwned>(params: &Value) -> Result<T> {
    serde_json::from_value(params.clone()).map_err(|e| bad("", e.to_string()))
}
fn reply(value: impl Serialize) -> Result<Value> {
    serde_json::to_value(value).map_err(|e| EngineError::Other(e.to_string()))
}
fn inspect_surface(_session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::mesh_diagnostics::{self, MAX_SURFACE_TRIANGLES, MAX_SURFACE_VERTICES};
    let mesh = params.get("mesh").ok_or_else(|| bad("mesh", "surface mesh is required"))?;
    for (field, limit) in [("vertices", MAX_SURFACE_VERTICES), ("triangles", MAX_SURFACE_TRIANGLES)] {
        let values = mesh.get(field).and_then(Value::as_array).ok_or_else(|| bad(field, "expected an array"))?;
        if values.len() > limit {
            return Err(bad(field, "surface exceeds inspection limits"));
        }
    }
    let input: SurfaceInput = parse(params)?;
    reply(mesh_diagnostics::inspect(&input.mesh).map_err(|error| bad("mesh", error.to_string()))?)
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
fn sweep(session: &mut Session, params: &Value) -> Result<Value> {
    let input: Sweep = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let studies = project.sweep(&input.id, &input.pointer, &input.values).map_err(|e| bad("", e.to_string()))?;
    let views = studies.iter().map(StudyView::try_from).collect::<std::result::Result<Vec<_>, _>>().map_err(|e| bad("", e.to_string()))?;
    session.multiphysics_jobs.next_study = project.next_id;
    session.doc_mut().multiphysics = Arc::new(project);
    reply(SweepReply { studies: views })
}
fn clone_study(session: &mut Session, params: &Value) -> Result<Value> {
    let input: CloneStudy = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let study = project.clone_study(&input.id, &input.name).map_err(|e| bad("id", e.to_string()))?;
    let view = StudyView::try_from(&study).map_err(|e| bad("", e.to_string()))?;
    session.multiphysics_jobs.next_study = project.next_id;
    session.doc_mut().multiphysics = Arc::new(project);
    reply(CloneReply { study: view })
}
fn wizard(session: &mut Session, params: &Value) -> Result<Value> {
    let input: Wizard = parse(params)?;
    let (case, report) = solidcraft_multiphysics::wizard::enclosure_case(&input.spec).map_err(|e| bad("spec", e.to_string()))?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let study = project.create_from_case(&input.name, &input.question, case).map_err(|e| bad("", e.to_string()))?;
    let view = StudyView::try_from(&study).map_err(|e| bad("", e.to_string()))?;
    session.multiphysics_jobs.next_study = project.next_id;
    session.doc_mut().multiphysics = Arc::new(project);
    reply(WizardReply { study: view, report })
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
fn query_result(session: &mut Session, params: &Value) -> Result<Value> {
    let input: solidcraft_multiphysics::results::QueryRequest = parse(params)?;
    let answer = solidcraft_multiphysics::results::query(&session.doc.multiphysics, &input).map_err(|error| bad("query", error.to_string()))?;
    reply(answer)
}
fn compare_studies(session: &mut Session, params: &Value) -> Result<Value> {
    let input: CompareStudies = parse(params)?;
    let table =
        solidcraft_multiphysics::results::compare_studies(&session.doc.multiphysics, &input.ids).map_err(|error| bad("ids", error.to_string()))?;
    reply(table)
}
fn inspect_result(session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::vtk::{self, Association};
    let input: InspectResult = parse(params)?;
    solidcraft_multiphysics::validate_artifact_name(&input.artifact).map_err(|error| bad("artifact", error.to_string()))?;
    let study = session.doc.multiphysics.study(&input.id).ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    let run = study.runs.iter().find(|run| run.id == input.run_id).ok_or_else(|| bad("run_id", "Select a retained result"))?;
    let artifact = run
        .artifacts
        .iter()
        .find(|artifact| artifact.name == input.artifact)
        .ok_or_else(|| bad("artifact", "Artifact is not retained in this run"))?;
    let content = if input.artifact.ends_with(".vtu") {
        let grid = vtk::parse_vtu(&artifact.text).map_err(|error| bad("artifact", error.to_string()))?;
        let cell_types = grid.cell_types.iter().copied().collect::<std::collections::BTreeSet<_>>().into_iter().collect();
        ResultInspectionContent::Vtu {
            points: grid.points.len(),
            cells: grid.offsets.len(),
            connectivity_indices: grid.connectivity.len(),
            cell_types,
            fields: grid
                .fields
                .into_iter()
                .map(|field| ResultField {
                    name: field.name,
                    association: match field.association {
                        Association::Point => "point",
                        Association::Cell => "cell",
                    }
                    .into(),
                    data_type: field.data_type,
                    components: field.components,
                    tuples: field.values.len() / field.components,
                })
                .collect(),
        }
    } else if input.artifact.ends_with(".pvd") {
        let frames = vtk::parse_pvd(&artifact.text)
            .map_err(|error| bad("artifact", error.to_string()))?
            .into_iter()
            .map(|frame| {
                let retained = run.artifacts.iter().any(|item| item.name == frame.file);
                ResultFrame { coordinate: frame.coordinate, file: frame.file, part: frame.part, retained }
            })
            .collect();
        ResultInspectionContent::Pvd { frames }
    } else {
        return Err(bad("artifact", "Choose a retained ASCII VTU or PVD artifact"));
    };
    let hash = solidcraft_multiphysics::case_hash(&study.case).map_err(|error| bad("id", error.to_string()))?;
    reply(ResultInspection {
        id: input.id,
        run_id: input.run_id,
        artifact: input.artifact,
        captured_study_revision: run.study_revision,
        current_study_revision: study.revision,
        case_hash: run.case_hash.clone(),
        stale: run.study_revision != study.revision || run.case_hash != hash,
        content,
    })
}
