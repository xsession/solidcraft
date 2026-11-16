use super::*;
use serde_json::{Value, json};
use solidcraft_multiphysics::{Artifact, RunEvidence, SOLVER_REVISION, SOLVER_VERSION};
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicUsize, Ordering},
};

#[test]
fn surface_inspection_command_reports_connectivity_and_rejects_bad_indices() {
    let mut session = Session::default();
    let mesh = json!({"vertices":[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],"triangles":[[0,1,2]]});
    let report = session.execute("multiphysics.surface.inspect", &json!({"mesh":mesh})).unwrap();
    assert_eq!(report["boundary_edges"], 3);
    assert_eq!(report["indexed_edge_pairing_complete"], false);
    let bad = json!({"vertices":[[0.,0.,0.]],"triangles":[[0,1,2]]});
    assert!(session.execute("multiphysics.surface.inspect", &json!({"mesh":bad})).is_err());
}

#[derive(Default)]
struct FakeHost {
    requests: Arc<Mutex<Vec<FrozenRun>>>,
    outcome: Arc<Mutex<Option<Outcome>>>,
    cancellations: Arc<AtomicUsize>,
}
struct FakeJob {
    outcome: Arc<Mutex<Option<Outcome>>>,
    cancellations: Arc<AtomicUsize>,
}
impl Host for FakeHost {
    fn start(&self, request: FrozenRun) -> std::result::Result<Box<dyn Job>, String> {
        self.requests.lock().unwrap().push(request);
        Ok(Box::new(FakeJob { outcome: self.outcome.clone(), cancellations: self.cancellations.clone() }))
    }
}
impl Job for FakeJob {
    fn poll(&mut self) -> Option<Outcome> {
        self.outcome.lock().unwrap().take()
    }
    fn cancel(&mut self) {
        self.cancellations.fetch_add(1, Ordering::SeqCst);
    }
}
fn evidence(request: &FrozenRun) -> RunEvidence {
    let manifest = json!({"schema":"cfd-solvers-rs/run/1","rust_port_version":SOLVER_VERSION,"upstream_revision":SOLVER_REVISION,
        "backend":"native_rust_serial_cpu","complete_clone":false,"execution_completed":true,"solver":request.study.case["solver"],
        "linear_solver_converged":true,"steady_state_converged":null,"metrics":{"residual_rms":0.0},"files":["solution.csv","manifest.json"],"case":request.study.case});
    RunEvidence {
        id: request.id.clone(),
        study_revision: request.study.revision,
        case_hash: request.case_hash.clone(),
        input_case: request.study.case.clone(),
        artifacts: vec![
            Artifact { name: "solution.csv".into(), text: "x_m,potential\n0,0\n0.5,0.25\n1,0\n".into() },
            Artifact { name: "manifest.json".into(), text: serde_json::to_string_pretty(&manifest).unwrap() },
        ],
        manifest,
    }
}
fn create(s: &mut Session) -> String {
    s.execute("multiphysics.study.create", &json!({"template":"poisson1d","name":"Potential study"})).unwrap()["study"]["id"].as_str().unwrap().into()
}
fn update(s: &mut Session, id: &str) {
    let study = s.doc.multiphysics.study(id).unwrap().clone();
    let mut case = study.case;
    case["source"] = json!(3.0);
    s.execute(
        "multiphysics.study.update",
        &json!({"id":id,"expected_revision":study.revision,"name":study.name,"question":"What is the maximum potential?","case":case}),
    )
    .unwrap();
}
fn setup() -> (Session, Arc<FakeHost>, String) {
    let mut s = Session::default();
    let host = Arc::new(FakeHost::default());
    s.set_multiphysics_host(Some(host.clone()));
    let id = create(&mut s);
    (s, host, id)
}
fn complete(s: &mut Session, host: &FakeHost) {
    let request = host.requests.lock().unwrap().last().unwrap().clone();
    *host.outcome.lock().unwrap() = Some(Outcome::Completed(evidence(&request)));
    s.poll_multiphysics_jobs();
}

#[test]
fn authoring_is_atomic_undoable_and_portable_without_geometry_evaluation() {
    let mut s = Session::default();
    let id = create(&mut s);
    let initial = s.doc.clone();
    assert_eq!(s.undo.len(), 1);
    let evaluations = s.evaluations;
    update(&mut s, &id);
    assert_eq!(s.doc.multiphysics.study(&id).unwrap().revision, 2);
    assert_eq!(s.undo.len(), 2);
    assert_eq!(s.evaluations, evaluations);
    let after = s.doc.clone();
    assert!(s.execute("multiphysics.study.update", &json!({"id":id,"expected_revision":1,"name":"stale","question":"","case":{}})).is_err());
    assert_eq!(s.doc, after);
    s.execute("edit.undo", &json!({})).unwrap();
    assert_eq!(s.doc.multiphysics, initial.multiphysics);
    s.execute("edit.redo", &json!({})).unwrap();
    assert_eq!(s.doc.multiphysics, after.multiphysics);
    let packed = solidcraft_io::write_design_container(&s.doc).unwrap();
    assert_eq!(solidcraft_io::read_design(&packed).unwrap().multiphysics, s.doc.multiphysics);
    assert!(solidcraft_doc::Document::from_json(&solidcraft_doc::Document::new("Old").to_json()).unwrap().multiphysics.is_empty());
}
#[test]
fn undo_does_not_reuse_study_identity_and_failed_apply_does_not_create_history() {
    let mut s = Session::default();
    let first = create(&mut s);
    s.execute("edit.undo", &json!({})).unwrap();
    let second = create(&mut s);
    assert_ne!(first, second);
    let before = s.doc.clone();
    let history = s.undo.len();
    assert!(
        s.execute(
            "multiphysics.study.update",
            &json!({"id":second,"expected_revision":1,"name":"valid","question":"","case":{"solver":"poisson1d","elements":10000000}})
        )
        .is_err()
    );
    assert_eq!(s.doc, before);
    assert_eq!(s.undo.len(), history);
}
#[test]
fn results_keep_only_current_evidence_and_artifacts_survive_reopen() {
    let (mut s, host, id) = setup();
    let job = s.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    assert_eq!(job["state"], "running");
    complete(&mut s, &host);
    assert_eq!(s.multiphysics_jobs()[0].state, "succeeded");
    s.execute("multiphysics.run.keep", &json!({"job_id":job["id"]})).unwrap();
    assert_eq!(s.doc.multiphysics.study(&id).unwrap().runs.len(), 1);
    assert_eq!(solidcraft_io::read_design(&solidcraft_io::write_design_container(&s.doc).unwrap()).unwrap().multiphysics, s.doc.multiphysics);
    update(&mut s, &id);
    assert!(s.multiphysics_jobs()[0].stale);
    assert!(s.execute("multiphysics.run.keep", &json!({"job_id":job["id"]})).is_err());
    assert_eq!(s.doc.multiphysics.study(&id).unwrap().runs.len(), 1);
}
#[test]
fn cancellation_rejects_late_completion_and_document_replacement_retires_workers() {
    let (mut s, host, id) = setup();
    let job = s.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    s.execute("multiphysics.run.cancel", &json!({"job_id":job["id"]})).unwrap();
    assert_eq!(s.multiphysics_jobs()[0].state, "cancelling");
    complete(&mut s, &host);
    assert_eq!(s.multiphysics_jobs()[0].state, "cancelled");
    assert!(s.multiphysics_jobs()[0].evidence.is_none());
    s.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    let cancelled = host.cancellations.load(Ordering::SeqCst);
    s.execute("file.new", &json!({"name":"Another design"})).unwrap();
    assert!(s.multiphysics_jobs().is_empty());
    assert!(s.multiphysics_host_available());
    assert!(host.cancellations.load(Ordering::SeqCst) > cancelled);
}
#[test]
fn forged_worker_identity_and_changed_study_never_attach() {
    let (mut s, host, id) = setup();
    s.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    let request = host.requests.lock().unwrap().last().unwrap().clone();
    let mut forged = evidence(&request);
    forged.id = "other-worker".into();
    *host.outcome.lock().unwrap() = Some(Outcome::Completed(forged));
    s.poll_multiphysics_jobs();
    assert_eq!(s.multiphysics_jobs()[0].state, "failed");
    assert!(s.multiphysics_jobs()[0].evidence.is_none());
    s.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    update(&mut s, &id);
    complete(&mut s, &host);
    let last = s.multiphysics_jobs().pop().unwrap();
    assert!(last.stale);
    assert!(s.execute("multiphysics.run.keep", &json!({"job_id":last.id})).is_err());
}
#[test]
fn exports_respect_agent_scope_and_preserve_historical_inputs() {
    let (mut s, host, id) = setup();
    let job = s.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    let original = host.requests.lock().unwrap().last().unwrap().study.case.clone();
    complete(&mut s, &host);
    s.execute("multiphysics.run.keep", &json!({"job_id":job["id"]})).unwrap();
    update(&mut s, &id);
    let base = std::env::temp_dir().join(format!("solidcraft-mp-export-{}-{}", std::process::id(), s.uid));
    std::fs::create_dir_all(base.join("allowed")).unwrap();
    solidcraft_io::vfs::scoped(Some(solidcraft_io::vfs::PathScope::new([base.join("allowed")])), || {
        assert!(s.execute("multiphysics.case.export", &json!({"id":id,"path":base.join("outside.json")})).is_err());
        let output = base.join("allowed/result");
        s.execute("multiphysics.result.export", &json!({"id":id,"run_id":job["id"],"directory":output})).unwrap();
        let back: Value = serde_json::from_slice(&std::fs::read(output.join("case.json")).unwrap()).unwrap();
        assert_eq!(back, original);
        assert!(output.join("manifest.json").is_file());
        assert!(output.join("solution.csv").is_file());
        assert!(!output.join("state.pvsm").exists());
        assert!(s.execute("multiphysics.result.export", &json!({"id":id,"run_id":job["id"],"directory":output})).is_err());
    });
    std::fs::remove_dir_all(base).unwrap();
}
#[test]
fn missing_companion_and_scratch_sessions_cannot_launch_work() {
    let mut s = Session::default();
    let id = create(&mut s);
    assert!(s.execute("multiphysics.run.start", &json!({"id":id})).is_err());
    let host = Arc::new(FakeHost::default());
    s.set_multiphysics_host(Some(host.clone()));
    assert!(s.scratch().execute("multiphysics.run.start", &json!({"id":id})).is_err());
    assert!(host.requests.lock().unwrap().is_empty());
}

#[test]
fn command_replies_keep_authoring_and_run_metadata_without_dataset_bytes() {
    let (mut session, host, id) = setup();
    let job = session.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    complete(&mut session, &host);
    let kept = session.execute("multiphysics.run.keep", &json!({"job_id":job["id"]})).unwrap();
    assert_eq!(kept["study"]["id"], id);
    assert_eq!(kept["study"]["case"]["solver"], "poisson1d");
    let retained = &kept["study"]["runs"][0];
    assert_eq!(retained["id"], job["id"]);
    assert_eq!(retained["artifact_count"], 2);
    assert_eq!(retained["execution_completed"], true);
    assert_eq!(retained["linear_solver_converged"], true);
    assert!(retained["steady_state_converged"].is_null());
    assert_eq!(retained["stale"], false);
    assert!(retained.get("manifest").is_none());
    assert!(retained.get("artifacts").is_none());
    assert!(!serde_json::to_string(&kept).unwrap().contains("0.5,0.25"));
    let listing = session.execute("multiphysics.study.list", &json!({})).unwrap();
    let listed = &listing["project"]["studies"][0];
    assert_eq!(listed["name"], "Potential study");
    assert_eq!(listed["solver"], "poisson1d");
    assert_eq!(listed["case_hash"].as_str().unwrap().len(), 64);
    assert!(listed.get("case").is_none());
    assert_eq!(listed["runs"][0]["artifact_count"], 2);
    assert!(!serde_json::to_string(&listing).unwrap().contains("0.5,0.25"));
    update(&mut session, &id);
    let listing = session.execute("multiphysics.study.list", &json!({})).unwrap();
    assert_eq!(listing["project"]["studies"][0]["runs"][0]["stale"], true);
}

#[test]
fn retained_result_inspection_reports_original_field_associations_without_dataset_bytes() {
    let (mut session, host, id) = setup();
    let job = session.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    let request = host.requests.lock().unwrap().last().unwrap().clone();
    let mut run = evidence(&request);
    let vtu = r#"<VTKFile type="UnstructuredGrid"><UnstructuredGrid><Piece NumberOfPoints="2" NumberOfCells="1"><PointData><DataArray type="Float64" Name="potential" format="ascii">0 0.25</DataArray></PointData><CellData><DataArray type="Float64" Name="flux" NumberOfComponents="2" format="ascii">2 3</DataArray></CellData><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0</DataArray></Points><Cells><DataArray type="Int64" Name="connectivity" format="ascii">0 1</DataArray><DataArray type="Int64" Name="offsets" format="ascii">2</DataArray><DataArray type="UInt8" Name="types" format="ascii">3</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>"#;
    let pvd = r#"<VTKFile type="Collection"><Collection><DataSet timestep="2.5" part="0" file="solution.vtu"/><DataSet timestep="3.5" part="1" file="missing.vtu"/></Collection></VTKFile>"#;
    run.artifacts.push(Artifact { name: "solution.vtu".into(), text: vtu.into() });
    run.artifacts.push(Artifact { name: "series.pvd".into(), text: pvd.into() });
    run.artifacts.push(Artifact { name: "invalid.vtu".into(), text: "<invalid/>".into() });
    run.manifest["files"] = json!(["solution.csv", "manifest.json", "solution.vtu", "series.pvd", "invalid.vtu"]);
    run.artifacts.iter_mut().find(|artifact| artifact.name == "manifest.json").unwrap().text = serde_json::to_string_pretty(&run.manifest).unwrap();
    *host.outcome.lock().unwrap() = Some(Outcome::Completed(run));
    session.poll_multiphysics_jobs();
    session.execute("multiphysics.run.keep", &json!({"job_id":job["id"]})).unwrap();
    let history = session.undo.len();
    let inspect = |artifact: &str| json!({"id":id,"run_id":job["id"],"artifact":artifact});
    let grid = session.execute("multiphysics.result.inspect", &inspect("solution.vtu")).unwrap();
    assert_eq!(grid["content"]["format"], "vtu");
    assert_eq!(grid["content"]["points"], 2);
    assert_eq!(grid["content"]["cells"], 1);
    assert_eq!(grid["content"]["cell_types"], json!([3]));
    assert_eq!(grid["content"]["fields"][0]["association"], "point");
    assert_eq!(grid["content"]["fields"][1]["association"], "cell");
    assert_eq!(grid["content"]["fields"][1]["components"], 2);
    assert_eq!(grid["content"]["fields"][1]["tuples"], 1);
    assert!(grid["content"]["fields"][0].get("values").is_none());
    assert!(grid["content"]["fields"][0].get("units").is_none());
    let collection = session.execute("multiphysics.result.inspect", &inspect("series.pvd")).unwrap();
    assert_eq!(collection["content"]["frames"][0]["coordinate"], 2.5);
    assert_eq!(collection["content"]["frames"][0]["retained"], true);
    assert_eq!(collection["content"]["frames"][1]["retained"], false);
    assert!(session.execute("multiphysics.result.inspect", &inspect("invalid.vtu")).is_err());
    assert!(session.execute("multiphysics.result.inspect", &inspect("solution.csv")).is_err());
    assert!(session.execute("multiphysics.result.inspect", &inspect("missing.vtu")).is_err());
    assert_eq!(session.undo.len(), history, "inspection is read-only");
    update(&mut session, &id);
    let historical = session.execute("multiphysics.result.inspect", &inspect("solution.vtu")).unwrap();
    assert_eq!(historical["stale"], true);
    assert_eq!(historical["captured_study_revision"], request.study.revision);
}

#[test]
fn newer_design_keeps_unknown_numerical_evidence_while_current_format_rejects_it() {
    let (mut session, host, id) = setup();
    session.execute("multiphysics.run.start", &json!({"id":id})).unwrap();
    let request = host.requests.lock().unwrap().last().unwrap().clone();
    let mut future = evidence(&request);
    future.manifest["rust_port_version"] = json!("99.0.0");
    let manifest = serde_json::to_string_pretty(&future.manifest).unwrap();
    future.artifacts.iter_mut().find(|artifact| artifact.name == "manifest.json").unwrap().text = manifest;
    let mut raw = serde_json::to_value(session.doc.as_ref()).unwrap();
    raw["name"] = json!("CAD with future numerical evidence");
    raw["params"] = json!([{"name":"part_length","expr":"10","unit":"mm","model":false}]);
    raw["multiphysics"]["studies"][0]["runs"] = json!([future]);
    let original_section = raw["multiphysics"].clone();
    assert!(solidcraft_doc::Document::from_json(&serde_json::to_string(&raw).unwrap()).is_err());
    raw["format"] = json!("solidcraft/9999");
    let opened = solidcraft_doc::Document::from_json(&serde_json::to_string(&raw).unwrap()).unwrap();
    assert!(opened.multiphysics.is_empty());
    assert_eq!(opened.name, "CAD with future numerical evidence");
    assert_eq!(opened.params[0].name, "part_length");
    assert_eq!(opened.params[0].expr, "10");
    assert!(opened.newer_file.as_ref().unwrap().warnings().iter().any(|warning| warning.contains("multiphysics")));
    assert_eq!(serde_json::from_str::<Value>(&opened.to_json()).unwrap()["multiphysics"], original_section);
    let container = solidcraft_io::write_design_container(&opened).unwrap();
    let reopened = solidcraft_io::read_design(&container).unwrap();
    assert!(reopened.multiphysics.is_empty());
    assert_eq!(reopened.name, opened.name);
    assert_eq!(reopened.params, opened.params);
    assert_eq!(serde_json::from_str::<Value>(&reopened.to_json()).unwrap()["multiphysics"], original_section);
}
