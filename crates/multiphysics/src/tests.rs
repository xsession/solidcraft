use super::*;
use serde_json::json;

fn study() -> Study {
    Project::default().create("poisson1d", "Potential").expect("template")
}

fn evidence(study: &Study) -> RunEvidence {
    let manifest = json!({
        "schema":"cfd-solvers-rs/run/1", "rust_port_version":SOLVER_VERSION,
        "upstream_revision":SOLVER_REVISION, "backend":"native_rust_serial_cpu",
        "solver":"poisson1d", "complete_clone":false, "execution_completed":true,
        "linear_solver_converged":null, "steady_state_converged":null,
        "case":study.case, "metrics":{"wall_time_s":0.01}, "files":["solution.csv","manifest.json"]
    });
    RunEvidence {
        id: "run-1".into(),
        study_revision: study.revision,
        case_hash: case_hash(&study.case).expect("hash"),
        input_case: study.case.clone(),
        artifacts: vec![
            Artifact { name: "solution.csv".into(), text: "x_m,potential\n0,0\n1,0\n".into() },
            Artifact { name: "manifest.json".into(), text: serde_json::to_string(&manifest).expect("manifest") },
        ],
        manifest,
    }
}

fn refresh_manifest(run: &mut RunEvidence) {
    run.artifacts.iter_mut().find(|artifact| artifact.name == "manifest.json").expect("manifest artifact").text =
        serde_json::to_string(&run.manifest).expect("serialize");
}

#[test]
fn all_templates_have_a_catalogue_profile_and_bounded_preflight() {
    let templates = templates().expect("templates");
    assert_eq!(templates.len(), 23);
    let profiles = catalog().expect("catalogue")["profiles"].as_array().expect("profiles").clone();
    assert_eq!(profiles.len(), 34);
    let mut solvers = BTreeSet::new();
    let mut ids = BTreeSet::new();
    for template in templates {
        assert!(ids.insert(template.id.clone()));
        assert!(profiles.iter().any(|profile| profile["id"] == template.profile));
        solvers.insert(template.case["solver"].as_str().expect("solver").to_owned());
        let mut project = Project::default();
        let study = project.create(&template.id, &template.name).expect("create");
        assert!(preflight(&study).errors.is_empty());
        assert!(!preflight(&study).warnings.is_empty());
    }
    assert_eq!(solvers.len(), 18);
}

#[test]
fn create_remove_never_reuses_a_durable_id() {
    let mut project = Project::default();
    assert!(project.is_empty());
    let a = project.create("poisson1d", "  A  ").expect("create");
    assert_eq!(a.id, "mp-study-1");
    assert_eq!(a.name, "A");
    project.remove(&a.id).expect("remove");
    let b = project.create("poisson1d", "B").expect("create second");
    assert_eq!(b.id, "mp-study-2");
    assert_ne!(a.id, b.id);
    project.validate().expect("project");
}

#[test]
fn stale_apply_and_invalid_apply_leave_the_project_unchanged() {
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    let before = project.clone();
    assert_eq!(project.update(&a.id, a.revision + 1, "B", "", a.case.clone()), Err(Error::Stale));
    assert!(project.update(&a.id, a.revision, "", "", a.case.clone()).is_err());
    assert!(project.update(&a.id, a.revision, "B", "", json!([])).is_err());
    assert_eq!(project, before);
}

#[test]
fn incomplete_cases_remain_savable_and_noop_apply_preserves_revision() {
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    let same = project.update(&a.id, a.revision, &a.name, &a.question, a.case.clone()).expect("noop");
    assert_eq!(same.revision, a.revision);
    let draft = project.update(&a.id, a.revision, "Draft", "Find potential", json!({})).expect("save draft");
    assert_eq!(draft.revision, a.revision + 1);
    assert!(!preflight(&draft).errors.is_empty());
    let bytes = serde_json::to_vec(&project).expect("save");
    let reopened: Project = serde_json::from_slice(&bytes).expect("reopen");
    reopened.validate().expect("valid draft");
    assert_eq!(project, reopened);
}

#[test]
fn missing_or_unconnected_solver_does_not_report_executable_readiness() {
    let mut a = study();
    a.case = json!({"solver":"spice_circuit"});
    assert!(preflight(&a).errors.iter().any(|error| error.contains("no native case adapter")));
    a.case = json!({"solver":"maxwell1d"});
    assert!(preflight(&a).errors.iter().any(|error| error.contains("selected profile")));
}

#[test]
fn hostile_mesh_steps_nesting_and_text_are_bounded() {
    assert!(validate_case(&json!({"mesh":{"dimensions":[u64::MAX,2,2]}})).is_err());
    assert!(validate_case(&json!({"mesh":{"nx":512,"ny":512}})).is_err());
    assert!(validate_case(&json!({"steps":MAX_STEPS + 1})).is_err());
    assert!(validate_case(&json!({"mesh":{"dimensions":[128,128,2]},"steps":2000})).is_err());
    assert!(validate_case(&json!({"note":"x".repeat(MAX_CASE_BYTES + 1)})).is_err());
    let mut deep = Value::Null;
    for _ in 0..40 {
        deep = json!([deep]);
    }
    assert!(validate_case(&json!({"nested":deep})).is_err());
}

#[test]
fn invalid_id_sequences_duplicates_and_capacity_are_rejected() {
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    project.next_id = 1;
    assert!(project.validate().is_err());
    project.next_id = 2;
    project.studies.push(a);
    assert!(project.validate().is_err());
    let mut full = Project::default();
    for _ in 0..MAX_STUDIES {
        full.create("poisson1d", "A").expect("create bounded study");
    }
    let before = full.clone();
    assert!(full.create("poisson1d", "Overflow").is_err());
    assert_eq!(full, before);
}

#[test]
fn raw_case_hash_is_deterministic_and_input_sensitive() {
    let a = json!({"a":1,"b":2});
    let b: Value = serde_json::from_str("{\"b\":2,\"a\":1}").expect("case");
    assert_eq!(case_hash(&a).expect("hash"), case_hash(&b).expect("hash"));
    assert_ne!(case_hash(&a).expect("hash"), case_hash(&json!({"a":2,"b":2})).expect("hash"));
    assert_eq!(case_hash(&json!({})).expect("hash"), "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a");
}

#[test]
fn completion_retains_convergence_outcome_without_claiming_physical_validation() {
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    let mut run = evidence(&a);
    run.manifest["linear_solver_converged"] = json!(false);
    refresh_manifest(&mut run);
    let admitted = project.attach_evidence(&a.id, run).expect("completed evidence");
    assert_eq!(admitted.runs[0].manifest["linear_solver_converged"], false);
    assert_eq!(admitted.revision, a.revision);
    assert!(admitted.runs[0].manifest.get("physical_validation").is_none());
}

#[test]
fn normalized_solver_defaults_are_allowed_but_authored_changes_are_rejected() {
    let a = study();
    let mut run = evidence(&a);
    run.manifest["case"]["unused_default"] = json!(0);
    refresh_manifest(&mut run);
    validate_evidence(&a, &run).expect("extra defaults");
    run.manifest["case"]["source"] = json!(42.0);
    refresh_manifest(&mut run);
    assert!(validate_evidence(&a, &run).is_err());
}

#[test]
fn normalized_solver_resource_defaults_are_not_new_authored_requests() {
    let study = Project::default().create("chemistry-batch", "Batch reactor").expect("template");
    assert!(study.case.pointer("/config/max_steps").is_none());
    let mut run = evidence(&study);
    run.manifest["solver"] = json!("chemistry_batch");
    run.manifest["case"]["config"]["max_steps"] = json!(100_000);
    refresh_manifest(&mut run);
    validate_evidence(&study, &run).expect("completed record retains omitted producer default");

    let mut explicit = study.clone();
    explicit.case["config"]["max_steps"] = json!(100_000);
    assert!(validate_case(&explicit.case).is_err());
    run.input_case = explicit.case.clone();
    run.case_hash = case_hash(&explicit.case).expect("digest still encodes raw data");
    assert!(validate_evidence(&explicit, &run).is_err());
}

#[test]
fn old_runs_survive_edits_and_reopen_as_historical_evidence() {
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    let run = evidence(&a);
    project.attach_evidence(&a.id, run.clone()).expect("attach");
    let b = project.update(&a.id, a.revision, "B", "Changed question", json!({"solver":"poisson1d","source":2.0})).expect("edit");
    assert_eq!(b.runs.len(), 1);
    validate_record(&b.runs[0]).expect("historical export remains valid");
    assert!(validate_evidence(&b, &run).is_err());
    let before = project.clone();
    assert!(project.attach_evidence(&a.id, run).is_err());
    assert_eq!(project, before);
    let saved = serde_json::to_string(&project).expect("save");
    let reopened: Project = serde_json::from_str(&saved).expect("reopen");
    reopened.validate().expect("historical run");
    assert_eq!(reopened, project);
}

#[test]
fn retained_historical_input_digest_and_manifest_are_checked_after_edits() {
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    project.attach_evidence(&a.id, evidence(&a)).expect("attach");
    project.update(&a.id, a.revision, "B", "", json!({"solver":"poisson1d","source":4})).expect("edit");
    let mut tampered = project.clone();
    tampered.studies[0].runs[0].input_case["source"] = json!(99.0);
    assert!(tampered.validate().is_err());
    tampered.studies[0].runs[0].case_hash = case_hash(&tampered.studies[0].runs[0].input_case).expect("rehash");
    assert!(tampered.validate().is_err());
}

#[test]
fn normalized_numeric_encodings_match_without_rounding_large_exact_integers() {
    let mut a = study();
    a.case["source"] = json!(2);
    let mut run = evidence(&a);
    run.manifest["case"]["source"] = json!(2.0);
    refresh_manifest(&mut run);
    validate_evidence(&a, &run).expect("normalized float");
    assert!(!authored_subset(&json!(9_007_199_254_740_993u64), &json!(9_007_199_254_740_992f64)));
}

#[test]
fn only_declared_d2q9_single_precision_scalars_are_rounded_for_matching() {
    let authored = json!({"solver":"lbm_d2q9","config":{"nx":32,"ny":32,"tau":0.7},"initial":{"type":"taylor_green","amplitude":0.03},"steps":100});
    let produced = json!({"solver":"lbm_d2q9","config":{"nx":32,"ny":32,"tau":0.7f32 as f64},"initial":{"type":"taylor_green","amplitude":0.03f32 as f64},"steps":100,"walls_y":false});
    assert!(normalized_authored_matches(&authored, &produced));
    assert!(normalized_authored_matches(&authored, &authored));
    let mut changed = produced.clone();
    changed["steps"] = json!(101);
    assert!(!normalized_authored_matches(&authored, &changed));
    assert!(!normalized_authored_matches(&json!({"solver":"poisson1d","source":0.03}), &json!({"solver":"poisson1d","source":0.03f32 as f64})));
}

#[test]
fn admission_rejects_bad_origin_paths_files_and_metric_structure() {
    let a = study();
    for (key, value) in [
        ("schema", json!("other")),
        ("rust_port_version", json!("9.0")),
        ("backend", json!("sycl")),
        ("execution_completed", json!(false)),
        ("complete_clone", json!(true)),
        ("solver", json!("maxwell1d")),
        ("metrics", json!({"bad":"0"})),
    ] {
        let mut run = evidence(&a);
        run.manifest[key] = value;
        refresh_manifest(&mut run);
        assert!(validate_evidence(&a, &run).is_err(), "{key}");
    }
    for name in ["../escape.csv", "/tmp/result", "x\\y", "C:drive", ".hidden", "", "CON.csv", "NUL", "Com1.vtu", "LPT9.txt", "alias."] {
        assert!(validate_artifact_name(name).is_err(), "{name}");
    }
    let mut run = evidence(&a);
    run.artifacts.pop();
    assert!(validate_evidence(&a, &run).is_err());
    let mut run = evidence(&a);
    run.artifacts[0].name = "manifest.json".into();
    assert!(validate_evidence(&a, &run).is_err());
    let mut run = evidence(&a);
    run.manifest["metrics"]["different"] = json!(1);
    assert!(validate_evidence(&a, &run).is_err());
}

#[test]
fn oversized_evidence_and_future_revisions_are_rejected_atomically() {
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    let before = project.clone();
    let mut run = evidence(&a);
    run.artifacts[0].text = "x".repeat(MAX_ARTIFACT_BYTES + 1);
    assert!(project.attach_evidence(&a.id, run).is_err());
    let mut run = evidence(&a);
    run.study_revision += 1;
    assert!(project.attach_evidence(&a.id, run).is_err());
    assert_eq!(project, before);
}

#[test]
fn final_project_budget_rejects_create_and_update_without_mutating_state() {
    let mut empty = Project::default();
    let before = empty.clone();
    let empty_bytes = serde_json::to_vec(&empty).expect("empty size").len();
    assert!(empty.create_with_budget("poisson1d", "A", empty_bytes + 1).is_err());
    assert_eq!(empty, before);
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    let before = project.clone();
    let current_bytes = serde_json::to_vec(&project).expect("size").len();
    let unchanged = project.update_with_budget(&a.id, a.revision, &a.name, &a.question, a.case.clone(), current_bytes).expect("exact-budget noop");
    assert_eq!(unchanged.revision, a.revision);
    assert_eq!(project, before);
    assert!(project.update_with_budget(&a.id, a.revision, "A", "", json!({"note":"x".repeat(2048)}), current_bytes + 100).is_err());
    assert_eq!(project, before);
}

#[test]
fn aggregate_budget_counts_frozen_cases_manifests_and_artifact_encoding() {
    let mut project = Project::default();
    let a = project.create("poisson1d", "A").expect("create");
    let run = evidence(&a);
    let mut admitted = project.clone();
    admitted.attach_evidence(&a.id, run.clone()).expect("attach");
    let serialized_bytes = serde_json::to_vec(&admitted).expect("serialized size").len();
    let artifact_bytes = run.artifacts.iter().map(|artifact| artifact.text.len()).sum::<usize>();
    let base_bytes = serde_json::to_vec(&project).expect("base size").len();
    assert!(serialized_bytes > base_bytes + artifact_bytes);
    let before = project.clone();
    assert!(project.attach_with_budget(&a.id, run, base_bytes + artifact_bytes).is_err());
    assert_eq!(project, before);
    admitted.validate_with_budget(serialized_bytes).expect("exact payload fits");
    assert!(admitted.validate_with_budget(serialized_bytes - 1).is_err());
}
