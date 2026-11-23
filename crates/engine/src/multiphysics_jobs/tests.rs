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

#[test]
fn reference_feature_prepare_and_wizard_creation_preserve_atomic_history() {
    let mut s = Session::default();
    let case = solidcraft_multiphysics::Project::default().create("flow-channel", "Base").unwrap().case;
    let document = s.doc.clone();
    let undo = s.undo.len();
    let prepared = s.execute("multiphysics.case.feature.prepare", &json!({"case":case,"kind":"porous","name":"filter"})).unwrap();
    assert!(Arc::ptr_eq(&document, &s.doc));
    assert_eq!(s.undo.len(), undo);
    assert!(s.multiphysics_jobs().is_empty());
    let mut invalid = case.clone();
    invalid["mesh"]["dimensions"] = json!([u64::MAX, u64::MAX, 2]);
    for rejected in [invalid, json!({"solver":"poisson1d","cells":8})] {
        assert!(s.execute("multiphysics.study.create", &json!({"template":"flow-channel","name":"Rejected","case":rejected})).is_err());
        assert_eq!(*s.doc, *document);
        assert_eq!(s.undo.len(), undo);
    }
    let reply = s.execute("multiphysics.study.create", &json!({"template":"flow-channel","name":"Reviewed","case":prepared["case"]})).unwrap();
    assert_eq!(reply["study"]["id"], "mp-study-1", "failed creation must not allocate identities");
    assert_eq!(s.doc.multiphysics.studies[0].case, prepared["case"]);
    assert_eq!(s.undo.len(), undo + 1);
    s.execute("edit.undo", &json!({})).unwrap();
    assert!(s.doc.multiphysics.studies.is_empty());
    s.execute("edit.redo", &json!({})).unwrap();
    assert_eq!(s.doc.multiphysics.studies[0].case, prepared["case"]);
    let saved = solidcraft_io::write_design_container(&s.doc).unwrap();
    let reopened = solidcraft_io::read_design(&saved).unwrap();
    assert_eq!(reopened.multiphysics, s.doc.multiphysics);
}

#[test]
fn cad_bodies_become_case_parts_without_touching_the_document() {
    let mut s = Session::default();
    s.execute("solid.box", &json!({"length":20,"width":10,"height":10})).unwrap();
    let names = s.execute("multiphysics.cad.bodies", &json!({})).unwrap()["bodies"].clone();
    let body = names[0].as_str().unwrap().to_string();
    let case = solidcraft_multiphysics::Project::default().create("flow-channel", "Base").unwrap().case;
    let document = s.doc.clone();
    let undo = s.undo.len();
    let reply = s
        .execute(
            "multiphysics.case.component.import",
            &json!({"case":case,"bodies":[body],"options":{"fit_domain":true,"margin_mm":10.0,"cell_mm":5.0}}),
        )
        .unwrap();
    assert!(Arc::ptr_eq(&document, &s.doc));
    assert_eq!(s.undo.len(), undo);
    assert_eq!(reply["edit"]["issues"], json!([]), "{reply}");
    // A box is a grid block; 20 × 10 × 10 mm plus 10 mm each side in 5 mm cells.
    assert_eq!(reply["parts"][0][1], "obstacle");
    let dimensions: Vec<u64> = serde_json::from_value(reply["edit"]["case"]["mesh"]["dimensions"].clone()).unwrap();
    let mut sorted = dimensions.clone();
    sorted.sort();
    assert_eq!(sorted, [6, 6, 8]);
    assert_eq!(reply["edit"]["case"]["obstacles"][0]["min"], json!([2, 2, 2]));
    assert_eq!(reply["edit"]["case"]["obstacles"][0]["max"], json!(dimensions.iter().map(|n| n - 2).collect::<Vec<_>>()));
    assert!(s.execute("multiphysics.case.component.import", &json!({"case":case,"bodies":["missing"]})).is_err());
}

/// A square duct along x: 100 mm long, 40 × 40 mm outside, 4 mm walls, open at both ends.
fn duct_session() -> Session {
    let mut s = Session::default();
    s.execute("solid.box", &json!({"length":100,"width":40,"height":40,"center":[50,20,20],"body_name":"Duct"})).unwrap();
    s.execute("solid.box", &json!({"length":120,"width":32,"height":32,"center":[50,20,20],"operation":"cut"})).unwrap();
    s
}

#[test]
fn a_study_of_the_model_builds_from_the_design_and_follows_its_edits() {
    let mut s = duct_session();
    let reply = s.execute("multiphysics.model.create", &json!({"name":"Duct flow","cell_mm":4.0})).unwrap();
    let id = reply["study"]["id"].as_str().unwrap().to_string();
    // The domain is the duct's exact box in 4 mm cells, and its 4 mm walls are one cell each,
    // so the grid holds the duct's 100 × (40² − 32²) = 57 600 mm³ exactly.
    assert_eq!(reply["study"]["case"]["mesh"]["dimensions"], json!([25, 10, 10]), "{reply}");
    assert!((reply["solid_fraction"].as_f64().unwrap() - 1.0).abs() < 1e-6, "{reply}");
    assert_eq!(reply["study"]["model"]["analysis"], "internal");
    let revision = |s: &Session| s.doc.multiphysics.study(&id).unwrap().revision;
    // Conditions go on the duct's end faces, picked by a point on each.
    let inlet = json!({"id":id,"expected_revision":revision(&s),"name":"Inlet","face":{"body":"Duct","point":[0.0,20.0,2.0]},"role":{"kind":"inlet","speed_m_per_s":0.5}});
    s.execute("multiphysics.model.condition", &inlet).unwrap();
    let outlet = json!({"id":id,"expected_revision":revision(&s),"name":"Outlet","face":{"body":"Duct","point":[100.0,20.0,2.0]},"role":{"kind":"outlet","pressure_m2_per_s2":0.0}});
    let reply = s.execute("multiphysics.model.condition", &outlet).unwrap();
    assert_eq!(reply["placed"], json!([["Inlet", "left"], ["Outlet", "right"]]), "{reply}");
    assert_eq!(reply["notes"], json!([]), "{reply}");
    // A face inside the domain (the channel's wall at y = 4 mm) becomes a patch on the duct
    // itself: exactly its triangles, every one on that plane.
    let inside = json!({"id":id,"expected_revision":revision(&s),"name":"Side","face":{"body":"Duct","point":[50.0,4.0,20.0]},"role":{"kind":"outlet","pressure_m2_per_s2":0.0}});
    let side = s.execute("multiphysics.model.condition", &inside).unwrap();
    assert_eq!(side["placed"][2], json!(["Side", "Side"]), "{side}");
    let condition = &side["study"]["case"]["face_conditions"][0];
    assert_eq!(condition["body"], "Duct");
    let body = side["study"]["case"]["bodies"].as_array().unwrap().iter().find(|b| b["name"] == "Duct").unwrap().clone();
    let triangles = condition["triangles"].as_array().unwrap();
    assert!(!triangles.is_empty());
    for t in triangles {
        for v in body["triangles"][t.as_u64().unwrap() as usize].as_array().unwrap() {
            // Case coordinates are metres from the domain origin at the duct's corner.
            let y = body["vertices"][v.as_u64().unwrap() as usize][1].as_f64().unwrap();
            assert!((y - 0.004).abs() < 1e-9, "triangle {t} leaves the face: y = {y}");
        }
    }
    // Remove it again: the two end conditions stay.
    let remove = json!({"id":id,"expected_revision":revision(&s),"name":"Side","remove":true});
    assert_eq!(s.execute("multiphysics.model.condition", &remove).unwrap()["placed"], reply["placed"]);
    // The case matches the model; editing the design makes it out of date and Run refuses.
    assert_eq!(s.execute("multiphysics.model.status", &json!({"id":id})).unwrap()["current"], true);
    s.execute("solid.box", &json!({"length":10,"width":10,"height":10,"center":[50,60,20],"operation":"new","body_name":"Block"})).unwrap();
    let status = s.execute("multiphysics.model.status", &json!({"id":id})).unwrap();
    assert_eq!(status["current"], false, "{status}");
    let why = s.execute("multiphysics.run.start", &json!({"id":id})).unwrap_err().to_string();
    assert!(why.contains("model changed"), "{why}");
    // Building again takes the new body in (the box now reaches y = 65 mm) and is undoable.
    let rebuilt = s.execute("multiphysics.model.build", &json!({"id":id})).unwrap();
    assert_eq!(s.execute("multiphysics.model.status", &json!({"id":id})).unwrap()["current"], true);
    assert_ne!(rebuilt["study"]["case"]["mesh"], reply["study"]["case"]["mesh"]);
    s.execute("edit.undo", &json!({})).unwrap();
    assert_eq!(s.doc.multiphysics.study(&id).unwrap().case, reply["study"]["case"]);
    // The link and its conditions survive saving and reopening.
    let reopened = solidcraft_doc::Document::from_json(&s.doc.to_json()).unwrap();
    assert_eq!(reopened.multiphysics, s.doc.multiphysics);
}

#[test]
fn a_study_of_the_model_needs_solid_bodies_and_named_faces() {
    let mut s = Session::default();
    let why = s.execute("multiphysics.model.create", &json!({"name":"Empty"})).unwrap_err().to_string();
    assert!(why.contains("no bodies"), "{why}");
    assert!(s.doc.multiphysics.studies.is_empty());
    let mut s = duct_session();
    let why = s.execute("multiphysics.model.create", &json!({"name":"Missing","bodies":["Nope"]})).unwrap_err().to_string();
    assert!(why.contains("no body called Nope"), "{why}");
    assert!(s.doc.multiphysics.studies.is_empty());
    let id = s.execute("multiphysics.model.create", &json!({"name":"Duct","cell_mm":4.0})).unwrap()["study"]["id"].as_str().unwrap().to_string();
    let revision = s.doc.multiphysics.study(&id).unwrap().revision;
    let far = json!({"id":id,"expected_revision":revision,"name":"In","face":{"body":"Duct","point":[0.0,20.0,500.0]},"role":{"kind":"wall"}});
    let why = s.execute("multiphysics.model.condition", &far).unwrap_err().to_string();
    assert!(why.contains("pick a point on the face"), "{why}");
    let plain = s.execute("multiphysics.study.create", &json!({"template":"flow-channel","name":"Plain"})).unwrap()["study"]["id"]
        .as_str()
        .unwrap()
        .to_string();
    let why = s.execute("multiphysics.model.build", &json!({"id":plain})).unwrap_err().to_string();
    assert!(why.contains("not a study of the loaded model"), "{why}");
}

#[test]
fn a_leak_path_finds_the_open_ends_of_a_duct_and_none_once_they_are_capped() {
    let mut s = duct_session();
    let reply = s.execute("multiphysics.model.create", &json!({"name":"Leak","analysis":"external","margin_mm":12.0,"cell_mm":4.0})).unwrap();
    let id = reply["study"]["id"].as_str().unwrap().to_string();
    // From the channel's wall (y = 4 mm, inside) to the duct's outer wall (y = 0, outside):
    // the fluid joins them round the open ends, so the path leaves through x < 0 or x > 100.
    let leak = json!({"id":id,"from":{"body":"Duct","point":[50.0,4.0,20.0]},"to":{"body":"Duct","point":[50.0,0.0,20.0]}});
    let found = s.execute("multiphysics.model.leak", &leak).unwrap();
    assert_eq!(found["connected"], true, "{found}");
    let path = found["path_mm"].as_array().unwrap();
    assert!(
        path.iter().any(|p| {
            let x = p[0].as_f64().unwrap();
            !(0.0..=100.0).contains(&x)
        }),
        "{found}"
    );
    // The start is half a cell off the inner wall, on the channel side, not through the wall.
    assert!(path[0][1].as_f64().unwrap() > 4.0, "{found}");
    // Cap both ends with 4 mm plates: the inside is sealed.
    for x in [-2.0, 102.0] {
        s.execute("solid.box", &json!({"length":4,"width":40,"height":40,"center":[x,20,20],"operation":"new"})).unwrap();
    }
    let reply = s.execute("multiphysics.model.build", &json!({"id":id})).unwrap();
    assert_eq!(reply["issues"], json!([]), "{reply}");
    let sealed = s.execute("multiphysics.model.leak", &leak).unwrap();
    assert_eq!(sealed["connected"], false, "{sealed}");
}

#[test]
fn a_study_of_the_model_refines_around_faces_and_monitors_goals_on_them() {
    let mut s = duct_session();
    let id = s.execute("multiphysics.model.create", &json!({"name":"Goals","cell_mm":4.0})).unwrap()["study"]["id"].as_str().unwrap().to_string();
    let revision = |s: &Session| s.doc.multiphysics.study(&id).unwrap().revision;
    let face = |point: [f64; 3]| json!({"body":"Duct","point":point});
    for (name, point, role) in [
        ("Inlet", [0.0, 20.0, 2.0], json!({"kind":"inlet","speed_m_per_s":0.05})),
        ("Outlet", [100.0, 20.0, 2.0], json!({"kind":"outlet","pressure_m2_per_s2":0.0})),
    ] {
        s.execute("multiphysics.model.condition", &json!({"id":id,"expected_revision":revision(&s),"name":name,"face":face(point),"role":role}))
            .unwrap();
    }
    // Smaller cells along the channel's wall at y = 4 mm: its two y cells split (10 → 12), the
    // wall's whole length (25 → 50) and its 8 z cells from 4 to 36 mm (10 → 18).
    let refined = s
        .execute(
            "multiphysics.model.refine",
            &json!({"id":id,"expected_revision":revision(&s),"name":"Wall","face":face([50.0, 4.0, 20.0]),"factor":2,"margin_cells":0.0}),
        )
        .unwrap();
    assert_eq!(refined["study"]["case"]["mesh"]["dimensions"], json!([50, 12, 18]), "{refined}");
    // Goals on the inlet's face, on the wall (no condition there) and at a point, and their sum.
    for goal in [
        json!({"name":"inflow","kind":"face_volume_flow","face":face([0.0, 20.0, 2.0])}),
        json!({"name":"wall_p","kind":"face_mean_pressure","face":face([50.0, 4.0, 20.0])}),
        json!({"name":"mid_p","kind":"point_pressure","point_mm":[50.0, 20.0, 20.0]}),
        json!({"name":"both","kind":"equation","expression":"wall_p + mid_p"}),
    ] {
        let mut params = goal;
        params["id"] = json!(id);
        params["expected_revision"] = json!(revision(&s));
        s.execute("multiphysics.model.goal", &params).unwrap();
    }
    let case = s.doc.multiphysics.study(&id).unwrap().case.clone();
    let target = |name: &str| case["goals"].as_array().unwrap().iter().find(|g| g["name"] == name).map(|g| g["target"].clone());
    assert_eq!(target("inflow"), Some(json!("left")));
    assert_eq!(target("wall_p"), Some(json!("goal_wall_p")));
    assert_eq!(case["face_conditions"][0]["name"], "goal_wall_p");
    // The point is in case metres from the domain origin at the duct's corner.
    let mid = case["goals"].as_array().unwrap().iter().find(|g| g["name"] == "mid_p").unwrap()["point"].clone();
    assert_eq!(mid, json!([0.05, 0.02, 0.02]));
    // Removing the refinement goes back to the uniform 4 mm grid.
    let back = s.execute("multiphysics.model.refine", &json!({"id":id,"expected_revision":revision(&s),"name":"Wall","remove":true})).unwrap();
    assert_eq!(back["study"]["case"]["mesh"]["dimensions"], json!([25, 10, 10]));
    assert_eq!(back["study"]["case"]["mesh"].get("coordinates"), None);
}

#[test]
fn a_study_of_the_model_carries_heat_on_its_bodies_and_faces() {
    let mut s = duct_session();
    let id = s.execute("multiphysics.model.create", &json!({"name":"Heat","cell_mm":4.0})).unwrap()["study"]["id"].as_str().unwrap().to_string();
    let revision = |s: &Session| s.doc.multiphysics.study(&id).unwrap().revision;
    let face = |point: [f64; 3]| json!({"body":"Duct","point":point});
    for (name, point, role) in [
        ("Inlet", [0.0, 20.0, 2.0], json!({"kind":"inlet","speed_m_per_s":0.05})),
        ("Outlet", [100.0, 20.0, 2.0], json!({"kind":"outlet","pressure_m2_per_s2":0.0})),
    ] {
        s.execute("multiphysics.model.condition", &json!({"id":id,"expected_revision":revision(&s),"name":name,"face":face(point),"role":role}))
            .unwrap();
    }
    // Materials need heat on.
    let early = json!({"id":id,"expected_revision":revision(&s),"bodies":["Duct"],"conductivity_w_per_m_k":200.0});
    assert!(s.execute("multiphysics.model.material", &early).unwrap_err().to_string().contains("heat transfer on"));
    let on =
        s.execute("multiphysics.model.heat", &json!({"id":id,"expected_revision":revision(&s),"enabled":true,"inflow_temperature_k":300.0})).unwrap();
    // The duct has no material yet: the starter one, with a note.
    assert_eq!(on["study"]["case"]["thermal"]["solids"][0]["conductivity_w_per_m_k"], json!([1.0, 1.0, 1.0]), "{on}");
    assert!(on["notes"].as_array().unwrap().iter().any(|n| n.as_str().unwrap().contains("starter material")), "{on}");
    let material = json!({"id":id,"expected_revision":revision(&s),"bodies":["Duct"],"conductivity_w_per_m_k":200.0,"power_w":1.0});
    s.execute("multiphysics.model.material", &material).unwrap();
    // 5 W on the channel's wall at y = 4 mm, and the inlet's end face held at 310 K.
    let source = json!({"id":id,"expected_revision":revision(&s),"name":"Heater","face":face([50.0, 4.0, 20.0]),"power_w":5.0});
    s.execute("multiphysics.model.heat_source", &source).unwrap();
    let hot = json!({"id":id,"expected_revision":revision(&s),"name":"Hot inlet","face":face([0.0, 20.0, 2.0]),"condition":{"kind":"temperature","temperature_k":310.0}});
    s.execute("multiphysics.model.thermal_face", &hot).unwrap();
    let goal = json!({"id":id,"expected_revision":revision(&s),"name":"duct_max","kind":"body_max_temperature","body":"Duct"});
    let reply = s.execute("multiphysics.model.goal", &goal).unwrap();
    let case = &reply["study"]["case"];
    assert_eq!(case["thermal"]["inflow_temperature_k"], json!(300.0));
    assert_eq!(case["thermal"]["solids"], json!([{"obstacle":"Duct","conductivity_w_per_m_k":[200.0,200.0,200.0],"power_w":1.0}]));
    assert_eq!(case["thermal"]["boundaries"], json!([{"patch":"left","condition":{"kind":"temperature","temperature_k":310.0}}]));
    let heater = &case["thermal"]["face_sources"][0];
    assert_eq!((heater["name"].clone(), heater["body"].clone(), heater["power_w"].clone()), (json!("Heater"), json!("Duct"), json!(5.0)));
    // The source covers exactly the wall's triangles: every one on the plane y = 4 mm.
    let body = case["bodies"].as_array().unwrap().iter().find(|b| b["name"] == "Duct").unwrap().clone();
    for t in heater["triangles"].as_array().unwrap() {
        for v in body["triangles"][t.as_u64().unwrap() as usize].as_array().unwrap() {
            let y = body["vertices"][v.as_u64().unwrap() as usize][1].as_f64().unwrap();
            assert!((y - 0.004).abs() < 1e-9, "triangle {t} leaves the face: y = {y}");
        }
    }
    assert_eq!(case["goals"][0], json!({"name":"duct_max","kind":"solid_max_temperature","target":"Duct"}));
    if let Ok(path) = std::env::var("SC_DUMP_HEATED_ENGINE") {
        let _ = std::fs::write(path, case.to_string());
    }
    // A thermal condition on a face inside the domain is refused with the reason.
    let inside = json!({"id":id,"expected_revision":revision(&s),"name":"Wall","face":face([50.0, 4.0, 20.0]),"condition":{"kind":"adiabatic"}});
    assert!(s.execute("multiphysics.model.thermal_face", &inside).unwrap_err().to_string().contains("domain boundary"));
    // Heat off: the thermal block and the temperature goal go.
    let off = s.execute("multiphysics.model.heat", &json!({"id":id,"expected_revision":revision(&s),"enabled":false})).unwrap();
    assert_eq!(off["study"]["case"].get("thermal"), None);
    assert_eq!(off["study"]["case"].get("goals"), None);
    assert_eq!(off["link"].get("heat"), None);
}

/// A grid of `n` cells of `cell` metres with a velocity along x of 0 below cell `step` and 1
/// m/s from it on, as the companion writes a flow field.
fn step_vtu(n: [usize; 3], cell: f64, step: usize) -> String {
    let (mut points, mut connectivity, mut offsets, mut velocity) = (String::new(), String::new(), String::new(), String::new());
    let mut count = 0;
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                for (dx, dy, dz) in [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)] {
                    points.push_str(&format!("{} {} {} ", (i + dx) as f64 * cell, (j + dy) as f64 * cell, (k + dz) as f64 * cell));
                }
                for v in 0..8 {
                    connectivity.push_str(&format!("{} ", count * 8 + v));
                }
                count += 1;
                offsets.push_str(&format!("{} ", count * 8));
                velocity.push_str(if i < step { "0 0 0 " } else { "1 0 0 " });
            }
        }
    }
    format!(
        r#"<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian"><UnstructuredGrid><Piece NumberOfPoints="{}" NumberOfCells="{count}"><CellData><DataArray type="Float64" Name="velocity_m_per_s" NumberOfComponents="3" format="ascii">{velocity}</DataArray></CellData><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{points}</DataArray></Points><Cells><DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray><DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray><DataArray type="UInt8" Name="types" format="ascii">{}</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>"#,
        count * 8,
        "12 ".repeat(count)
    )
}

#[test]
fn a_study_of_the_model_takes_a_mesh_level_and_adapts_to_a_run() {
    let mut s = duct_session();
    let id = s.execute("multiphysics.model.create", &json!({"name":"Adapt","cell_mm":4.0})).unwrap()["study"]["id"].as_str().unwrap().to_string();
    let revision = |s: &Session| s.doc.multiphysics.study(&id).unwrap().revision;
    let dimensions = |reply: &Value| reply["study"]["case"]["mesh"]["dimensions"].clone();
    // A level clears the typed 4 mm cells: level 5 is finer than the automatic level 3.
    let level3 = s.execute("multiphysics.model.set", &json!({"id":id,"expected_revision":revision(&s),"mesh_level":3})).unwrap();
    assert_eq!(level3["link"].get("cell_mm"), None);
    let level5 = s.execute("multiphysics.model.set", &json!({"id":id,"expected_revision":revision(&s),"mesh_level":5})).unwrap();
    let count = |d: Value| d.as_array().unwrap().iter().map(|n| n.as_u64().unwrap()).product::<u64>();
    assert!(count(dimensions(&level5)) > 5 * count(dimensions(&level3)), "{} {}", dimensions(&level3), dimensions(&level5));
    // Back to 4 mm cells, then keep a run whose velocity steps up at x = 48 mm.
    let reply = s.execute("multiphysics.model.set", &json!({"id":id,"expected_revision":revision(&s),"cell_mm":4.0})).unwrap();
    assert_eq!(dimensions(&reply), json!([25, 10, 10]));
    let study = s.doc.multiphysics.study(&id).unwrap().clone();
    let capabilities: Vec<&str> = solidcraft_multiphysics::airflow::requirements(&study.case).iter().map(|r| r.capability).collect();
    let manifest = json!({"schema":"cfd-solvers-rs/run/1","rust_port_version":SOLVER_VERSION,"upstream_revision":SOLVER_REVISION,
        "backend":"native_rust_serial_cpu","complete_clone":false,"execution_completed":true,"solver":"collocated_flow3d",
        "linear_solver_converged":null,"steady_state_converged":null,"capabilities":capabilities,"metrics":{"wall_time_s":0.01},
        "files":["flow_000010.vtu","manifest.json"],"case":study.case});
    let run = RunEvidence {
        id: "first".into(),
        study_revision: study.revision,
        case_hash: solidcraft_multiphysics::case_hash(&study.case).unwrap(),
        input_case: study.case.clone(),
        artifacts: vec![
            Artifact { name: "flow_000010.vtu".into(), text: step_vtu([25, 10, 10], 0.004, 12) },
            Artifact { name: "manifest.json".into(), text: manifest.to_string() },
        ],
        manifest,
    };
    let mut project = s.doc.multiphysics.as_ref().clone();
    project.attach_evidence(&id, run).unwrap();
    s.doc_mut().multiphysics = Arc::new(project);
    let adapted = s.execute("multiphysics.model.adapt", &json!({"id":id,"expected_revision":revision(&s),"run_id":"first"})).unwrap();
    // The cells on both sides of the step, x from 44 to 52 mm, through the whole section.
    assert_eq!(adapted["adapted"]["flagged"], json!(200), "{}", adapted["adapted"]);
    let regions: Vec<[[f64; 3]; 2]> = serde_json::from_value(adapted["adapted"]["regions_mm"].clone()).unwrap();
    assert_eq!(regions.len(), 1);
    for (got, want) in regions[0].iter().flatten().zip([44.0, 0.0, 0.0, 52.0, 40.0, 40.0]) {
        assert!((got - want).abs() < 1e-9, "{regions:?}");
    }
    assert_eq!(adapted["link"]["refinements"][0]["name"], "adapt_1");
    // Those two x cells split in two; y and z split everywhere the box reaches (all of them).
    assert_eq!(dimensions(&adapted), json!([27, 20, 20]));
    // Adapting again replaces the earlier boxes; remove takes them away.
    let again = s.execute("multiphysics.model.adapt", &json!({"id":id,"expected_revision":revision(&s),"run_id":"first","factor":3})).unwrap();
    assert_eq!(again["link"]["refinements"].as_array().unwrap().len(), 1);
    let removed = s.execute("multiphysics.model.adapt", &json!({"id":id,"expected_revision":revision(&s),"run_id":"first","remove":true})).unwrap();
    assert_eq!(dimensions(&removed), json!([25, 10, 10]));
    assert!(s.execute("multiphysics.model.adapt", &json!({"id":id,"expected_revision":revision(&s),"run_id":"nope"})).is_err());
}

#[test]
fn a_study_of_the_model_takes_a_fan_and_global_goals() {
    let mut s = duct_session();
    let id = s.execute("multiphysics.model.create", &json!({"name":"Fan","cell_mm":4.0})).unwrap()["study"]["id"].as_str().unwrap().to_string();
    let revision = |s: &Session| s.doc.multiphysics.study(&id).unwrap().revision;
    let face = |point: [f64; 3]| json!({"body":"Duct","point":point});
    let outlet = json!({"id":id,"expected_revision":revision(&s),"name":"Outlet","face":face([100.0, 20.0, 2.0]),"role":{"kind":"outlet","pressure_m2_per_s2":0.0}});
    s.execute("multiphysics.model.condition", &outlet).unwrap();
    // The fan on the duct's other end blows against the outlet.
    let fan = json!({"id":id,"expected_revision":revision(&s),"face":face([0.0, 20.0, 2.0]),"curve":[[0.0, 4.0], [4e-4, 3.0], [1e-3, 0.0]]});
    let reply = s.execute("multiphysics.model.fan", &fan).unwrap();
    let case = &reply["study"]["case"];
    assert_eq!((case["fan"]["patch"].clone(), case["fan"]["reference_patch"].clone()), (json!("left"), json!("right")), "{case}");
    for (name, kind) in [("v_max", "global_max_speed"), ("p_mean", "global_mean_pressure")] {
        s.execute("multiphysics.model.goal", &json!({"id":id,"expected_revision":revision(&s),"name":name,"kind":kind})).unwrap();
    }
    let case = s.doc.multiphysics.study(&id).unwrap().case.clone();
    assert_eq!(case["goals"], json!([{"name":"v_max","kind":"global_max_speed"},{"name":"p_mean","kind":"global_mean_pressure"}]));
    if let Ok(path) = std::env::var("SC_DUMP_FAN") {
        let _ = std::fs::write(path, case.to_string());
    }
    // A fan on a face inside the domain is refused; remove takes it away.
    let inside = json!({"id":id,"expected_revision":revision(&s),"face":face([50.0, 4.0, 20.0]),"curve":[[0.0, 4.0], [1e-3, 0.0]]});
    assert!(s.execute("multiphysics.model.fan", &inside).unwrap_err().to_string().contains("domain boundary"));
    let removed = s.execute("multiphysics.model.fan", &json!({"id":id,"expected_revision":revision(&s),"remove":true})).unwrap();
    assert_eq!(removed["study"]["case"].get("fan"), None);
}

#[test]
fn part_selections_follow_tags_and_a_study_of_one_follows_the_design() {
    let mut s = Session::default();
    // Two 10 mm cubes and a 4 mm one, apart.
    for (name, x, side) in [("U1", 0.0, 10.0), ("U2", 30.0, 10.0), ("C1", 60.0, 4.0)] {
        s.execute("solid.box", &json!({"length":side,"width":side,"height":side,"center":[x,0,0],"operation":"new","body_name":name})).unwrap();
    }
    let undo = s.undo.len();
    let tagged = s.execute("multiphysics.parts.tag", &json!({"bodies":["U1","U2"],"add_tags":["power"],"attributes":{"power_w":2.5}})).unwrap();
    assert_eq!(s.undo.len(), undo + 1);
    let u1 = tagged["parts"].as_array().unwrap().iter().find(|p| p["body"] == "U1").unwrap().clone();
    assert_eq!((u1["tags"].clone(), u1["attributes"].clone()), (json!(["power"]), json!({"power_w":2.5})), "{tagged}");
    // A live query on the tag, and one on exact size: the 4 mm cube is under 5 mm.
    let query = |rules: Value| json!({"name":"q","members":{"kind":"query","rules":rules}});
    let mut power = query(json!([{"rule":"tag","tag":"power"}]));
    power["name"] = json!("Power parts");
    s.execute("multiphysics.parts.selection", &power).unwrap();
    let mut small = query(json!([{"rule":"attribute","key":"size_mm","compare":"lt","value":5.0}]));
    small["name"] = json!("Small");
    let listed = s.execute("multiphysics.parts.selection", &small).unwrap();
    assert_eq!(listed["selections"][0]["bodies"], json!(["U1", "U2"]), "{listed}");
    assert_eq!(listed["selections"][1]["bodies"], json!(["C1"]), "{listed}");
    // A study of the selection builds from its two bodies, and follows a retag on rebuild.
    let reply = s
        .execute("multiphysics.model.create", &json!({"name":"Power","analysis":"external","margin_mm":10.0,"cell_mm":5.0,"selection":"Power parts"}))
        .unwrap();
    let id = reply["study"]["id"].as_str().unwrap().to_string();
    s.execute("multiphysics.parts.tag", &json!({"bodies":["C1"],"add_tags":["power"]})).unwrap();
    assert_eq!(s.execute("multiphysics.model.status", &json!({"id":id})).unwrap()["current"], false);
    let rebuilt = s.execute("multiphysics.model.build", &json!({"id":id})).unwrap();
    assert_eq!(rebuilt["study"]["model"]["selection"], "Power parts");
    // C1 at x = 60 mm now widens the domain.
    assert_ne!(rebuilt["study"]["case"]["mesh"], reply["study"]["case"]["mesh"], "{rebuilt}");
    // A selection a study uses cannot be removed; size_mm cannot be set; unknown bodies are refused.
    assert!(s.execute("multiphysics.parts.selection", &json!({"name":"Power parts","remove":true})).is_err());
    assert!(s.execute("multiphysics.parts.tag", &json!({"bodies":["U1"],"attributes":{"size_mm":1.0}})).is_err());
    assert!(s.execute("multiphysics.parts.tag", &json!({"bodies":["nope"],"add_tags":["x"]})).is_err());
    // The component table: 2 W in each of the three power parts, before any run.
    let revision = s.doc.multiphysics.study(&id).unwrap().revision;
    s.execute("multiphysics.model.heat", &json!({"id":id,"expected_revision":revision,"enabled":true})).unwrap();
    let revision = s.doc.multiphysics.study(&id).unwrap().revision;
    let material = json!({"id":id,"expected_revision":revision,"selection":"Power parts","conductivity_w_per_m_k":150.0,"power_w":2.0,"max_temperature_k":358.0});
    s.execute("multiphysics.model.material", &material).unwrap();
    let table = s.execute("multiphysics.parts.table", &json!({"id":id,"csv":true})).unwrap();
    assert_eq!(table["table"]["total_power_w"], json!(6.0), "{table}");
    assert_eq!(table["table"]["run"], Value::Null);
    let rows = table["table"]["rows"].as_array().unwrap();
    assert_eq!(rows.iter().map(|r| r["body"].as_str().unwrap()).collect::<Vec<_>>(), ["U1", "U2", "C1"]);
    assert_eq!((rows[0]["attributes"]["power_w"].clone(), rows[0]["limit_k"].clone()), (json!(2.5), json!(358.0)));
    assert!(table["csv"].as_str().unwrap().starts_with("body,in_study,tags,"), "{table}");
    let small = s.execute("multiphysics.parts.table", &json!({"id":id,"selection":"Small"})).unwrap();
    assert_eq!(small["table"]["rows"].as_array().unwrap().len(), 1);
    assert!(s.execute("multiphysics.parts.table", &json!({"id":id,"run_id":"missing"})).is_err());
    // Linked selection: a part selection becomes the design's selection, the one the tree, the
    // table and the 3-D view show; selecting is not an undo step.
    let undo = s.undo.len();
    let picked = s.execute("multiphysics.parts.select", &json!({"selection":"Small"})).unwrap();
    assert_eq!(picked["selected"], json!(["C1"]), "{picked}");
    let picked = s.execute("multiphysics.parts.select", &json!({"bodies":["U2"],"add":true})).unwrap();
    assert_eq!(picked["selected"], json!(["U2", "C1"]), "in the design's order: {picked}");
    assert_eq!(s.selection, vec![crate::Sel::Body { name: "C1".into() }, crate::Sel::Body { name: "U2".into() }]);
    assert!(s.execute("multiphysics.parts.select", &json!({"bodies":["nope"]})).is_err());
    assert!(s.execute("multiphysics.parts.select", &json!({"bodies":["U1"],"selection":"Small"})).is_err());
    assert_eq!(s.execute("multiphysics.parts.select", &json!({})).unwrap()["selected"], json!([]));
    // Ghosting is view state: see-through parts, not saved and not undone.
    let ghosted = s.execute("multiphysics.parts.ghost", &json!({"selection":"Power parts"})).unwrap();
    assert_eq!(ghosted["ghosted"], json!(["U1", "U2", "C1"]), "{ghosted}");
    let ghosted = s.execute("multiphysics.parts.ghost", &json!({"bodies":["U2"],"ghost":false})).unwrap();
    assert_eq!(ghosted["ghosted"], json!(["U1", "C1"]));
    let ghosted = s.execute("multiphysics.parts.ghost", &json!({"bodies":["U2"],"only":true})).unwrap();
    assert_eq!(ghosted["ghosted"], json!(["U2"]), "only: every other part solid");
    assert!(s.execute("multiphysics.parts.ghost", &json!({})).is_err(), "nothing named and not only");
    assert_eq!(s.undo.len(), undo);
    assert!(!s.doc.to_json().contains("ghost"), "not saved");
    assert_eq!(s.execute("multiphysics.parts.ghost", &json!({"only":true})).unwrap()["ghosted"], json!([]));
    // Everything survives saving and reopening.
    let reopened = solidcraft_doc::Document::from_json(&s.doc.to_json()).unwrap();
    assert_eq!(reopened.multiphysics, s.doc.multiphysics);
}
