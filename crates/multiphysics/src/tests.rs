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
    assert_eq!(templates.len(), 56);
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
    assert_eq!(solvers.len(), 45);
    for solver in crate::airflow::DOMAIN_SOLVERS {
        assert!(solvers.contains(*solver), "{solver} has no template");
    }
}

#[test]
fn create_remove_never_reuses_a_durable_id() {
    #[derive(serde::Serialize, serde::Deserialize)]
    struct SavedProject {
        #[serde(default, skip_serializing_if = "Project::is_default")]
        project: Project,
    }
    let mut project = Project::default();
    assert!(project.is_empty());
    assert!(project.is_default());
    assert_eq!(serde_json::to_value(SavedProject { project: project.clone() }).unwrap(), json!({}));
    let a = project.create("poisson1d", "  A  ").expect("create");
    assert_eq!(a.id, "mp-study-1");
    assert_eq!(a.name, "A");
    project.remove(&a.id).expect("remove");
    assert!(project.is_empty());
    assert!(!project.is_default(), "an empty project must retain its identity history");
    let saved = serde_json::to_value(SavedProject { project }).expect("save empty project");
    assert_eq!(saved["project"]["next_id"], 2);
    let mut project = serde_json::from_value::<SavedProject>(saved).expect("reopen empty project").project;
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
    a.case = json!({"solver":"spray_fvm"});
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
fn electronics_obstacles_admit_touching_blocks_but_reject_ambiguous_geometry() {
    let enclosure = templates().expect("templates").into_iter().find(|template| template.id == "flow-electronics-enclosure").expect("enclosure");
    validate_case(&enclosure.case).expect("touching PCB and components are valid");
    let mut overlapping = enclosure.case.clone();
    overlapping["obstacles"][1]["min"][2] = json!(1);
    assert!(validate_case(&overlapping).expect_err("overlap").to_string().contains("overlap"));
    let mut renamed = enclosure.case.clone();
    renamed["obstacles"][1]["name"] = json!("left");
    assert!(validate_case(&renamed).is_err());
    renamed["obstacles"][1]["name"] = json!("pcb_block");
    assert!(validate_case(&renamed).is_err());
    let mut shear = enclosure.case.clone();
    shear["mesh"]["xy_shear"] = json!(0.1);
    assert!(validate_case(&shear).is_err());
    let mut outside = enclosure.case.clone();
    outside["obstacles"][0]["max"][0] = json!(21);
    assert!(validate_case(&outside).is_err());
    let mut malformed = enclosure.case.clone();
    malformed["obstacles"][0]["min"] = json!([1, 2]);
    assert!(validate_case(&malformed).is_err());
    malformed = enclosure.case.clone();
    malformed["obstacles"][0]["min"][0] = json!(-1);
    assert!(validate_case(&malformed).is_err());
    malformed = enclosure.case.clone();
    malformed["mesh"]["lengths"]["x"] = json!(0);
    assert!(validate_case(&malformed).is_err());
    malformed = enclosure.case.clone();
    malformed["obstacles"] = json!(vec![enclosure.case["obstacles"][0].clone(); 129]);
    assert!(validate_case(&malformed).is_err());
    malformed = enclosure.case;
    malformed["mesh"]["dimensions"] = json!([0, 12, 8]);
    assert!(validate_case(&malformed).is_err());
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
fn producer_admission_preserves_legacy_runs_and_requires_obstacle_capability() {
    let a = study();
    let mut run = evidence(&a);
    run.manifest["rust_port_version"] = json!("0.2.0");
    refresh_manifest(&mut run);
    validate_evidence(&a, &run).expect("legacy scalar evidence remains valid");
    let enclosure = templates().expect("templates").into_iter().find(|template| template.id == "flow-electronics-enclosure").expect("enclosure");
    assert!(validate_producer(&run.manifest, &enclosure.case).is_err());
    run.manifest["rust_port_version"] = json!(SOLVER_VERSION);
    assert!(validate_producer(&run.manifest, &enclosure.case).is_err());
    run.manifest["capabilities"] = json!(["grid_obstacles_v1"]);
    validate_producer(&run.manifest, &enclosure.case).expect("new obstacle capability");
    run.manifest["rust_port_version"] = json!("0.3.0");
    validate_producer(&run.manifest, &enclosure.case).expect("0.3.0 obstacle runs remain admissible");
    run.manifest["rust_port_version"] = json!("0.2.0");
    assert!(validate_producer(&run.manifest, &enclosure.case).is_err());
    run.manifest["rust_port_version"] = json!(SOLVER_VERSION);
    run.manifest["capabilities"] = json!(["grid_obstacles_v1", "grid_obstacles_v1"]);
    assert!(validate_producer(&run.manifest, &enclosure.case).is_err());
    run.manifest["capabilities"] = json!([7]);
    assert!(validate_producer(&run.manifest, &a.case).is_err());
    run.manifest["capabilities"] = json!(["grid_obstacles_v1"]);
    run.manifest["rust_port_version"] = json!("999.0.0");
    assert!(validate_producer(&run.manifest, &a.case).is_err());
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

fn fan_enclosure() -> Value {
    templates().expect("templates").into_iter().find(|template| template.id == "flow-electronics-enclosure-fan").expect("fan template").case
}

#[test]
fn fan_and_porous_cases_require_their_own_producer_capabilities() {
    let case = fan_enclosure();
    let names: Vec<_> = airflow::requirements(&case).iter().map(|r| r.capability).collect();
    assert_eq!(names, ["grid_obstacles_v1", "fan_curve_v1", "porous_blocks_v1"]);
    let mut manifest = evidence(&study()).manifest;
    manifest["rust_port_version"] = json!("0.3.0");
    manifest["capabilities"] = json!(["grid_obstacles_v1", "fan_curve_v1", "porous_blocks_v1"]);
    assert!(validate_producer(&manifest, &case).is_err(), "fan and porous zones are 0.4.0 extensions");
    manifest["rust_port_version"] = json!(SOLVER_VERSION);
    validate_producer(&manifest, &case).expect("0.4.0 with all three capabilities");
    for missing in ["grid_obstacles_v1", "fan_curve_v1", "porous_blocks_v1"] {
        let names: Vec<_> = ["grid_obstacles_v1", "fan_curve_v1", "porous_blocks_v1"].into_iter().filter(|name| *name != missing).collect();
        manifest["capabilities"] = json!(names);
        let error = validate_producer(&manifest, &case).expect_err("missing capability");
        assert!(error.to_string().contains(missing), "{error}");
    }
    let mut fan_only = case.clone();
    fan_only["obstacles"] = json!([]);
    fan_only["porous_blocks"] = json!([]);
    manifest["capabilities"] = json!(["fan_curve_v1"]);
    validate_producer(&manifest, &fan_only).expect("a fan alone needs only its own capability");
}

#[test]
fn fan_curves_and_porous_blocks_are_bounded_and_unambiguous() {
    let base = fan_enclosure();
    validate_case(&base).expect("template admits");
    let rejected = |edit: &dyn Fn(&mut Value)| {
        let mut case = base.clone();
        edit(&mut case);
        validate_case(&case).is_err()
    };
    assert!(rejected(&|c| c["solver"] = json!("solid_heat3d")));
    assert!(rejected(&|c| c["fan"]["curve"] = json!([{"volume_flow_m3_per_s":0.0,"pressure_rise_pa":1.0}])));
    assert!(rejected(&|c| c["fan"]["curve"][1]["pressure_rise_pa"] = json!(5.0)), "rising pressure is not invertible");
    assert!(rejected(&|c| c["fan"]["curve"][1]["volume_flow_m3_per_s"] = json!(0.0)));
    assert!(rejected(&|c| c["fan"]["curve"][0]["pressure_rise_pa"] = json!(-1.0)));
    assert!(rejected(
        &|c| c["fan"]["curve"] = json!((0..65).map(|i| json!({"volume_flow_m3_per_s":i,"pressure_rise_pa":100 - i})).collect::<Vec<_>>())
    ));
    assert!(
        rejected(&|c| {
            c["fan"].as_object_mut().expect("fan").remove("pressure_convention");
        }),
        "convention must be explicit"
    );
    assert!(rejected(&|c| c["fan"]["pressure_convention"] = json!("dynamic")));
    assert!(rejected(&|c| c["fan"]["typo"] = json!(1)));
    assert!(rejected(&|c| c["fan"]["patch"] = json!("pcb_block")));
    assert!(rejected(&|c| c["fan"]["reference_patch"] = json!("left")));
    assert!(rejected(&|c| c["fan"]["initial_volume_flow_m3_per_s"] = json!(1.0)));
    assert!(rejected(&|c| c["fan"]["relaxation"] = json!(0.0)));
    assert!(rejected(&|c| c["fan"]["update_interval"] = json!(0)));
    assert!(rejected(&|c| c["fan"]["settle_tolerance"] = json!(2.0)));
    assert!(rejected(&|c| c["velocity_boundaries"].as_array_mut().expect("list").push(json!({"patch":"left","kind":"fixed_value"}))));
    assert!(rejected(&|c| c["pressure_boundaries"] = json!([{"patch":"right","kind":"zero_gradient"}])));
    assert!(rejected(&|c| c["pressure_boundaries"].as_array_mut().expect("list").push(json!({"patch":"left","kind":"fixed_value","value":0.0}))));
    assert!(rejected(&|c| c["porous_blocks"][0]["viscous_resistance_per_m2"] = json!(-1.0)));
    assert!(rejected(&|c| {
        c["porous_blocks"][0]["viscous_resistance_per_m2"] = json!(0.0);
        c["porous_blocks"][0]["inertial_resistance_per_m"] = json!(0.0);
    }));
    assert!(rejected(&|c| c["porous_blocks"][0]["max"] = json!([21, 12, 8])));
    assert!(rejected(&|c| c["porous_blocks"][0]["min"] = json!([17, 0, 0])), "overlaps the PCB block");
    assert!(rejected(&|c| c["porous_blocks"][0]["name"] = json!("pcb_block")));
    assert!(rejected(&|c| c["porous_blocks"][0]["name"] = json!("right")));
    assert!(rejected(&|c| {
        let first = c["porous_blocks"][0].clone();
        let mut second = first.clone();
        second["name"] = json!("second_grille");
        c["porous_blocks"] = json!([first, second]);
    }));
    assert!(rejected(&|c| c["mesh"]["xy_shear"] = json!(0.1)));
    assert!(rejected(&|c| c["porous_blocks"] = json!((0..33).map(|_| json!({})).collect::<Vec<_>>())));
    // Removing the optional features leaves an ordinary admitted case.
    let mut plain = base.clone();
    plain.as_object_mut().expect("case").remove("fan");
    plain["porous_blocks"] = json!([]);
    validate_case(&plain).expect("features are optional");
}

#[test]
fn conjugate_heat_cases_need_0_5_0_and_describe_every_solid() {
    let case = templates().expect("templates").into_iter().find(|t| t.id == "flow-electronics-enclosure-cht").expect("cht").case;
    let names: Vec<_> = airflow::requirements(&case).iter().map(|r| r.capability).collect();
    assert_eq!(names, ["grid_obstacles_v1", "fan_curve_v1", "porous_blocks_v1", "conjugate_heat_v1", "goals_v1"]);
    let mut manifest = evidence(&study()).manifest;
    manifest["capabilities"] = json!(["grid_obstacles_v1", "fan_curve_v1", "porous_blocks_v1", "conjugate_heat_v1", "goals_v1"]);
    manifest["rust_port_version"] = json!("0.4.0");
    assert!(validate_producer(&manifest, &case).is_err(), "conjugate heat is a 0.5.0 extension");
    manifest["rust_port_version"] = json!("0.5.0");
    assert!(validate_producer(&manifest, &case).is_err(), "goals are a 0.6.0 extension");
    let mut no_goals = case.clone();
    if let Some(object) = no_goals.as_object_mut() {
        object.remove("goals");
        object.remove("goal_control");
    }
    validate_producer(&manifest, &no_goals).expect("0.5.0 admits conjugate heat without goals");
    manifest["rust_port_version"] = json!(SOLVER_VERSION);
    validate_producer(&manifest, &case).expect("the current producer admits it");
    let rejected = |edit: &dyn Fn(&mut Value)| {
        let mut changed = case.clone();
        edit(&mut changed);
        validate_case(&changed).is_err()
    };
    assert!(
        rejected(&|c| {
            c["thermal"]["solids"].as_array_mut().expect("solids").pop();
        }),
        "an obstacle without a material"
    );
    assert!(rejected(&|c| c["thermal"]["solids"][0]["obstacle"] = json!("ghost")));
    assert!(rejected(&|c| c["thermal"]["solids"][0]["conductivity_w_per_m_k"] = json!([1.0, 0.0, 1.0])));
    assert!(rejected(&|c| c["thermal"]["solids"][0]["conductivity_w_per_m_k"] = json!(1.0)), "per-axis values are explicit");
    assert!(rejected(&|c| c["thermal"]["solids"][0]["power_w"] = json!(-1.0)));
    assert!(rejected(&|c| c["thermal"]["solids"][0]["max_temperature_k"] = json!(-40.0)), "limits are kelvin");
    assert!(rejected(&|c| c["thermal"]["inflow_temperature_k"] = json!(0.0)));
    assert!(rejected(&|c| c["thermal"]["fluid"]["conductivity_w_per_m_k"] = json!(0.0)));
    assert!(rejected(&|c| c["thermal"]["boundaries"][0]["patch"] = json!("top")), "duplicate patch");
    assert!(rejected(&|c| c["thermal"]["boundaries"][0]["patch"] = json!("pcb_block")));
    assert!(rejected(&|c| c["thermal"]["boundaries"][0]["condition"] = json!({"kind":"convection","h_w_per_m2_k":0.0,"ambient_k":298.15})));
    assert!(rejected(&|c| c["thermal"]["typo"] = json!(1)));
    assert!(rejected(&|c| c["solver"] = json!("solid_heat3d")));
}

// ---------------------------------------------------------------- geometry setup records

mod geometry_setup {
    use super::*;
    use crate::setup::*;

    fn face(body: &str, name: &str, lineage: &str, area: f64, wetted: bool) -> SnapshotFace {
        SnapshotFace {
            occurrence: vec!["enclosure-1".into()],
            body: body.into(),
            name: name.into(),
            lineage: lineage.into(),
            area_m2: area,
            fluid_boundary: wetted,
        }
    }

    /// An enclosure with a curved inlet fillet, a hidden internal rib and an exterior face
    /// that does not touch the fluid. Names follow persistent lineage, not face ordinals.
    fn snapshot(revision: u64, split_inlet: bool) -> TopologySnapshot {
        let mut faces = vec![
            face("case", "F1:-x", "F1:-x", 0.004, true),
            face("case", "F3:fillet", "F3:fillet", 0.0015, true),
            face("case", "F1:+x", "F1:+x", 0.006, true),
            face("case", "F1:walls", "F1:walls", 0.05, true),
            face("case", "F5:rib", "F5:rib", 0.002, true),
            face("case", "F9:outside", "F9:outside", 0.08, false),
            face("chip", "F2:top", "F2:top", 0.0004, true),
            face("heatsink", "F7:base", "F7:base", 0.0004, false),
        ];
        if split_inlet {
            faces.retain(|f| f.name != "F1:-x");
            faces.push(face("case", "F1:-x#[F8:slot]", "F1:-x", 0.0025, true));
            faces.push(face("case", "F1:-x#[F8:slot,1]", "F1:-x", 0.0015, true));
        }
        let bodies = ["case", "chip", "heatsink"]
            .into_iter()
            .map(|body| SnapshotBody { occurrence: vec!["enclosure-1".into()], body: body.into(), volume_m3: 1e-5 })
            .collect();
        TopologySnapshot { geometry_revision: revision, faces, bodies }
    }

    fn member(body: &str, face: Option<&str>) -> EntityRef {
        EntityRef { occurrence: vec!["enclosure-1".into()], body: body.into(), face: face.map(String::from) }
    }

    fn selection(id: &str, kind: EntityKind, members: Vec<EntityRef>) -> NamedSelection {
        NamedSelection { id: id.into(), name: id.replace('-', " "), kind, members }
    }

    fn wall() -> WallCondition {
        WallCondition { slip: WallSlip::NoSlip, motion: WallMotion::Stationary, thermal: WallThermal::Adiabatic }
    }

    fn complete_setup() -> Setup {
        Setup {
            selections: vec![
                selection("inlet-faces", EntityKind::Face, vec![member("case", Some("F1:-x")), member("case", Some("F3:fillet"))]),
                selection("outlet", EntityKind::Face, vec![member("case", Some("F1:+x"))]),
                selection(
                    "walls",
                    EntityKind::Face,
                    vec![member("case", Some("F1:walls")), member("case", Some("F5:rib")), member("chip", Some("F2:top"))],
                ),
                selection("chip-top", EntityKind::Face, vec![member("chip", Some("F2:top"))]),
                selection("heatsink-base", EntityKind::Face, vec![member("heatsink", Some("F7:base"))]),
                selection("board", EntityKind::Body, vec![member("chip", None)]),
            ],
            boundaries: vec![
                BoundaryAssignment {
                    id: "inlet".into(),
                    selection: "inlet-faces".into(),
                    condition: BoundaryCondition::Inlet { inflow: InletSpec::NormalSpeed { inward_m_per_s: 0.5 }, temperature_kelvin: Some(298.15) },
                },
                BoundaryAssignment {
                    id: "exhaust".into(),
                    selection: "outlet".into(),
                    condition: BoundaryCondition::Outlet { gauge_pressure_pa: 0.0, pressure: PressureKind::Static, backflow: Backflow::Prevent },
                },
                BoundaryAssignment { id: "no-slip".into(), selection: "walls".into(), condition: BoundaryCondition::Wall(wall()) },
            ],
            contacts: vec![ContactInterface {
                id: "tim".into(),
                side_a: "chip-top".into(),
                side_b: "heatsink-base".into(),
                resistance: ContactResistance::ThinLayer { thickness_m: 1e-4, conductivity_w_per_m_kelvin: 2.0 },
            }],
            conduction: vec![ConductionAssignment {
                id: "pcb-k".into(),
                selection: "board".into(),
                material: Some("FR-4 laminate".into()),
                conductivity: Conductivity::Orthotropic {
                    principal_w_per_m_kelvin: [30.0, 30.0, 0.3],
                    axes: [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]],
                },
            }],
            ..Default::default()
        }
    }

    fn codes(report: &SetupReport) -> Vec<&str> {
        report.issues.iter().map(|issue| issue.code.as_str()).collect()
    }

    #[test]
    fn complete_setup_covers_curved_and_hidden_faces_once() {
        let report = audit(&complete_setup(), &snapshot(4, false));
        assert!(report.applicable && report.ready, "{:?}", report.issues);
        assert_eq!(report.fluid_boundary_faces, 6);
        assert_eq!(report.explicit_faces, 6);
        assert_eq!(report.inherited_default_faces, 0);
        let inlet = report.boundaries.iter().find(|b| b.id == "inlet").expect("inlet");
        // Several faces, one of them curved: the area and implied flow integrate all of them.
        assert_eq!(inlet.face_count, 2);
        assert!((inlet.area_m2 - 0.0055).abs() < 1e-15);
        assert!((inlet.implied_volume_flow_m3_per_s.expect("flow") - 0.00275).abs() < 1e-15);
        let contact = &report.contacts[0];
        assert!((contact.area_specific_m2_kelvin_per_watt.expect("r") - 5e-5).abs() < 1e-18);
        assert!((contact.total_kelvin_per_watt.expect("R") - 0.125).abs() < 1e-12);
    }

    #[test]
    fn a_split_or_deleted_face_blocks_apply_until_reviewed() {
        let mut project = Project::default();
        let study = project.create("flow-electronics-enclosure", "Setup").expect("study");
        let (applied, _) = project.apply_setup(&study.id, study.revision, complete_setup(), &snapshot(4, false)).expect("apply");
        assert_eq!(applied.revision, study.revision + 1);
        assert_eq!(applied.setup.geometry_revision, 4);
        // A later feature splits the inlet face. Its lineage now has two pieces.
        let report = audit(&applied.setup, &snapshot(5, true));
        assert!(report.geometry_changed && !report.applicable && !report.ready);
        let inlet = report.selections.iter().find(|s| s.id == "inlet-faces").expect("selection");
        assert!(matches!(&inlet.members[0].resolution, MemberState::Ambiguous { candidates } if candidates.len() == 2));
        let before = project.clone();
        let error = project.apply_setup(&study.id, applied.revision, applied.setup.clone(), &snapshot(5, true)).expect_err("blocked");
        assert!(error.to_string().contains("inlet-faces"), "{error}");
        assert_eq!(project, before, "a refused Apply changes nothing");
        // Review: name both pieces explicitly. Nothing is reattached by proximity.
        let mut reviewed = applied.setup.clone();
        reviewed.selections[0].members =
            vec![member("case", Some("F1:-x#[F8:slot]")), member("case", Some("F1:-x#[F8:slot,1]")), member("case", Some("F3:fillet"))];
        let (study, report) = project.apply_setup(&study.id, applied.revision, reviewed, &snapshot(5, true)).expect("reviewed");
        assert!(report.ready, "{:?}", report.issues);
        assert_eq!(study.setup.geometry_revision, 5);
        // Deleting the fillet leaves a missing member, which blocks Apply as well.
        let mut gone = snapshot(6, true);
        gone.faces.retain(|f| f.name != "F3:fillet");
        let report = audit(&study.setup, &gone);
        assert!(!report.applicable && codes(&report).contains(&"unresolved_reference"));
        assert!(project.apply_setup(&study.id, study.revision, study.setup.clone(), &gone).is_err());
        assert!(matches!(project.apply_setup(&study.id, study.revision - 1, study.setup.clone(), &snapshot(5, true)), Err(Error::Stale)));
    }

    #[test]
    fn coverage_reports_gaps_defaults_conflicts_and_non_fluid_faces() {
        let geometry = snapshot(1, false);
        let mut setup = complete_setup();
        setup.boundaries.retain(|b| b.id != "no-slip");
        let report = audit(&setup, &geometry);
        assert!(report.applicable && !report.ready);
        assert_eq!(report.uncovered.len(), 3);
        setup.default_wall = Some(wall());
        let report = audit(&setup, &geometry);
        assert!(report.ready, "an explicit default is allowed");
        assert_eq!(report.inherited_default_faces, 3);
        assert!(codes(&report).contains(&"inherited_default_wall"));
        // Overlap never silently wins: the chip face assigned twice is a conflict.
        let mut setup = complete_setup();
        setup.boundaries.push(BoundaryAssignment {
            id: "hot-chip".into(),
            selection: "chip-top".into(),
            condition: BoundaryCondition::Wall(WallCondition { thermal: WallThermal::HeatFlux { w_per_m2: 2000.0 }, ..wall() }),
        });
        let report = audit(&setup, &geometry);
        assert!(report.applicable && !report.ready);
        assert_eq!(report.conflicts.len(), 1);
        assert_eq!(report.conflicts[0].assignments, ["no-slip", "hot-chip"]);
        // An exterior face is not a fluid boundary and cannot become an inlet.
        let mut setup = complete_setup();
        setup.selections[1].members.push(member("case", Some("F9:outside")));
        assert!(codes(&audit(&setup, &geometry)).contains(&"not_fluid_boundary"));
        // Without an outlet or opening there is no pressure reference.
        let mut setup = complete_setup();
        setup.boundaries.retain(|b| b.id != "exhaust");
        setup.default_wall = Some(wall());
        assert!(codes(&audit(&setup, &geometry)).contains(&"no_pressure_reference"));
    }

    #[test]
    fn inlet_meanings_stay_distinct() {
        assert_eq!(InletSpec::NormalSpeed { inward_m_per_s: 2.0 }.volume_flow(0.01), Some(0.02));
        assert_eq!(InletSpec::VolumeFlow { m3_per_s: 0.003 }.volume_flow(99.0), Some(0.003));
        assert_eq!(InletSpec::MassFlow { kg_per_s: 0.012, density_kg_per_m3: 1.2 }.volume_flow(1.0), Some(0.01));
        let vector = InletSpec::Velocity { vector_m_per_s: [1.0, 0.0, 0.0], frame: GLOBAL_FRAME.into() };
        assert_eq!(vector.volume_flow(1.0), None, "a fixed vector's flux depends on each face normal");
        let mut setup = complete_setup();
        setup.boundaries[0].condition = BoundaryCondition::Inlet { inflow: vector, temperature_kelvin: None };
        assert!(codes(&audit(&setup, &snapshot(1, false))).contains(&"fixed_vector_inlet"));
        for invalid in [
            InletSpec::NormalSpeed { inward_m_per_s: -1.0 },
            InletSpec::MassFlow { kg_per_s: 1.0, density_kg_per_m3: 0.0 },
            InletSpec::Velocity { vector_m_per_s: [1.0, 0.0, 0.0], frame: "rotor".into() },
            InletSpec::Velocity { vector_m_per_s: [0.0; 3], frame: GLOBAL_FRAME.into() },
        ] {
            setup.boundaries[0].condition = BoundaryCondition::Inlet { inflow: invalid, temperature_kelvin: None };
            assert!(setup.validate().is_err());
        }
        setup.boundaries[0].condition =
            BoundaryCondition::Inlet { inflow: InletSpec::NormalSpeed { inward_m_per_s: 1.0 }, temperature_kelvin: Some(-5.0) };
        assert!(setup.validate().is_err(), "temperatures are absolute");
    }

    #[test]
    fn contact_laws_convert_once_and_are_counted_once() {
        let area = 1e-4;
        let laws = [
            ContactResistance::Total { kelvin_per_watt: 0.5 },
            ContactResistance::AreaSpecific { m2_kelvin_per_watt: 5e-5 },
            ContactResistance::Conductance { w_per_m2_kelvin: 2e4 },
            ContactResistance::ThinLayer { thickness_m: 1e-4, conductivity_w_per_m_kelvin: 2.0 },
        ];
        for law in &laws {
            assert!((law.area_specific(area).expect("r") - 5e-5).abs() < 1e-18, "{law:?}");
            assert!((law.total(area).expect("R") - 0.5).abs() < 1e-12, "{law:?}");
        }
        assert!(laws[0].total(0.0).is_err());
        assert!(ContactResistance::Conductance { w_per_m2_kelvin: 0.0 }.area_specific(area).is_err());
        let mut setup = complete_setup();
        let mut reversed = setup.contacts[0].clone();
        reversed.id = "tim-again".into();
        std::mem::swap(&mut reversed.side_a, &mut reversed.side_b);
        setup.contacts.push(reversed);
        assert!(setup.validate().expect_err("double counted").to_string().contains("once"));
        let mut setup = complete_setup();
        setup.contacts[0].side_b = "chip-top".into();
        assert!(setup.validate().is_err());
        let mut geometry = snapshot(1, false);
        if let Some(base) = geometry.faces.iter_mut().find(|f| f.name == "F7:base") {
            base.area_m2 = 0.001;
        }
        assert!(codes(&audit(&complete_setup(), &geometry)).contains(&"contact_area_mismatch"));
    }

    #[test]
    fn orthotropic_conductivity_rotates_with_its_assignment() {
        let pcb = complete_setup().conduction[0].conductivity.clone();
        let tensor = pcb.world_tensor().expect("tensor");
        // Through-thickness (local z) lies along global y for this placement.
        let expected = [[30.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 30.0]];
        for (row, want) in tensor.iter().zip(expected) {
            for (value, w) in row.iter().zip(want) {
                assert!((value - w).abs() < 1e-12, "{tensor:?}");
            }
        }
        let (c, s) = (std::f64::consts::FRAC_PI_4.cos(), std::f64::consts::FRAC_PI_4.sin());
        let rotated = Conductivity::Orthotropic { principal_w_per_m_kelvin: [10.0, 2.0, 1.0], axes: [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]] };
        let t = rotated.world_tensor().expect("tensor");
        assert!((t[0][1] - 4.0).abs() < 1e-12 && (t[1][0] - t[0][1]).abs() < 1e-15, "K_xy = (k1 - k2) c s");
        assert!((t[0][0] + t[1][1] - 12.0).abs() < 1e-12, "the trace is invariant");
        assert!(is_symmetric_positive_definite(&t));
        for axes in [
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
            [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        ] {
            assert!(Conductivity::Orthotropic { principal_w_per_m_kelvin: [1.0; 3], axes }.world_tensor().is_err(), "{axes:?}");
        }
        assert!(Conductivity::Isotropic { w_per_m_kelvin: 0.0 }.world_tensor().is_err());
        assert!(!is_symmetric_positive_definite(&[[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]]));
        let mut setup = complete_setup();
        let mut second = setup.conduction[0].clone();
        second.id = "pcb-k-2".into();
        setup.conduction.push(second);
        assert!(codes(&audit(&setup, &snapshot(1, false))).contains(&"conflicting_conduction"));
    }

    #[test]
    fn setup_records_are_structurally_bounded_and_optional_in_saved_studies() {
        let study = Project::default().create("flow-electronics-enclosure", "Format").expect("study");
        let text = serde_json::to_string(&study).expect("serialize");
        assert!(!text.contains("\"setup\""), "an empty setup leaves the saved study unchanged");
        let mut with_setup = study.clone();
        with_setup.setup = complete_setup();
        let round: Study = serde_json::from_str(&serde_json::to_string(&with_setup).expect("serialize")).expect("parse");
        assert_eq!(round, with_setup);
        let invalid = [
            |s: &mut Setup| s.boundaries[0].selection = "nowhere".into(),
            |s: &mut Setup| s.boundaries[0].selection = "board".into(),
            |s: &mut Setup| s.conduction[0].selection = "outlet".into(),
            |s: &mut Setup| {
                let first = s.selections[0].members[0].clone();
                s.selections[0].members.push(first);
            },
            |s: &mut Setup| s.selections[0].members[0].face = None,
            |s: &mut Setup| s.selections[5].members[0].face = Some("F2:top".into()),
            |s: &mut Setup| s.selections[0].members.clear(),
            |s: &mut Setup| s.boundaries[1].id = "inlet".into(),
            |s: &mut Setup| s.selections[0].id = "bad\nid".into(),
            |s: &mut Setup| {
                s.default_wall = Some(WallCondition { thermal: WallThermal::Convection { w_per_m2_k: 0.0, ambient_kelvin: 300.0 }, ..wall() })
            },
            |s: &mut Setup| {
                s.default_wall = Some(WallCondition {
                    motion: WallMotion::Rotating { origin_m: [0.0; 3], axis: [0.0; 3], rad_per_s: 1.0, frame: GLOBAL_FRAME.into() },
                    ..wall()
                })
            },
        ];
        for (i, edit) in invalid.iter().enumerate() {
            let mut setup = complete_setup();
            edit(&mut setup);
            assert!(setup.validate().is_err(), "edit {i}");
            let mut study = study.clone();
            study.setup = setup;
            assert!(validate_study(&study).is_err(), "edit {i}");
        }
        let json = serde_json::json!({"kind":"inlet","inflow":{"spec":"normal_speed","inward_m_per_s":1.0},"speed":2});
        assert!(serde_json::from_value::<BoundaryCondition>(json).is_err(), "unknown members are refused");
    }
}

// ---------------------------------------------------------------- enclosure wizard

mod enclosure_wizard {
    use super::*;
    use crate::wizard::*;

    pub fn component(name: &str, min: [f64; 3], max: [f64; 3], power: f64, limit: Option<f64>) -> Component {
        Component { name: name.into(), min_mm: min, max_mm: max, power_w: power, conductivity_w_per_m_k: [15.0; 3], max_temperature_c: limit }
    }

    pub fn fan_spec() -> EnclosureSpec {
        EnclosureSpec {
            size_mm: [200.0, 120.0, 80.0],
            cooling: Cooling::Fan {
                inlet: Face::MinusX,
                opening_center_mm: [60.0, 40.0],
                opening_size_mm: [40.0, 40.0],
                curve: vec![
                    CurvePoint { litres_per_s: 0.0, pressure_pa: 4.0 },
                    CurvePoint { litres_per_s: 0.4, pressure_pa: 3.0 },
                    CurvePoint { litres_per_s: 0.8, pressure_pa: 1.5 },
                    CurvePoint { litres_per_s: 1.0, pressure_pa: 0.0 },
                ],
                static_pressure_curve: true,
                exhaust: Face::PlusX,
                exhaust_filter: Some(Filter { viscous_resistance_per_m2: 1e8, inertial_resistance_per_m: 2000.0 }),
            },
            ambient_c: 25.0,
            wall_h_w_per_m2_k: 5.0,
            components: vec![
                Component {
                    conductivity_w_per_m_k: [20.0, 20.0, 0.4],
                    ..component("pcb", [20.0, 10.0, 10.0], [180.0, 110.0, 20.0], 0.5, Some(105.0))
                },
                component("cpu", [50.0, 30.0, 20.0], [80.0, 60.0, 50.0], 3.0, Some(85.0)),
                component("regulator", [110.0, 70.0, 20.0], [150.0, 100.0, 40.0], 2.0, Some(100.0)),
            ],
            mesh_level: 2,
        }
    }

    pub fn natural_spec() -> EnclosureSpec {
        EnclosureSpec { cooling: Cooling::Natural { gravity_toward: Face::MinusZ }, mesh_level: 1, ..fan_spec() }
    }

    #[test]
    fn fan_wizard_builds_an_admissible_case_and_reports_every_choice() {
        let (case, report) = enclosure_case(&fan_spec()).expect("wizard");
        validate_case(&case).expect("admissible");
        let names: Vec<_> = airflow::requirements(&case).iter().map(|r| r.capability).collect();
        assert_eq!(names, ["grid_obstacles_v1", "fan_curve_v1", "porous_blocks_v1", "conjugate_heat_v1", "goals_v1"]);
        let cells: u64 = report.dimensions.iter().product();
        assert!((2_000..=8_000).contains(&cells), "{report:?}");
        // Time step respects both Courant 0.5 at the fan's free-delivery speed and diffusion 0.05.
        let h = report.cell_size_mm.iter().cloned().fold(f64::MAX, f64::min) / 1000.0;
        assert!(report.time_step_s <= 0.05 * h * h / 1.5e-5 + 1e-15);
        assert_eq!(report.components.len(), 3);
        for snapped in &report.components {
            assert!(snapped.max_shift_mm <= report.cell_size_mm.iter().cloned().fold(0.0, f64::max) / 2.0 + 1e-9, "{snapped:?}");
        }
        // Units: litres per second become m^3/s, Celsius limits become kelvin.
        assert_eq!(case["fan"]["curve"][3]["volume_flow_m3_per_s"], json!(0.001));
        let cpu = case["thermal"]["solids"].as_array().expect("solids").iter().find(|s| s["obstacle"] == "cpu").expect("cpu");
        assert!((cpu["max_temperature_k"].as_f64().expect("limit") - 358.15).abs() < 1e-9);
        // The shroud leaves exactly the opening exposed on the inlet face.
        let shrouds =
            case["obstacles"].as_array().expect("obstacles").iter().filter(|o| o["name"].as_str().is_some_and(|n| n.starts_with("shroud_"))).count();
        assert_eq!(shrouds, 4);
        assert_eq!(case["goal_control"]["stop_when_converged"], true);
        assert!(case["goals"].as_array().expect("goals").iter().any(|g| g["name"] == "cpu_max"));
    }

    #[test]
    fn natural_wizard_and_refusals() {
        let (case, report) = enclosure_case(&natural_spec()).expect("wizard");
        validate_case(&case).expect("admissible");
        assert_eq!(case["thermal"]["buoyancy"]["gravity_m_per_s2"], json!([0.0, 0.0, -9.81]));
        assert!(case.get("fan").is_none() && case["pressure_boundaries"] == json!([]));
        assert!(report.warnings.iter().any(|w| w.contains("pcb")), "a thin board is flagged: {:?}", report.warnings);
        let refused = |edit: &dyn Fn(&mut EnclosureSpec)| {
            let mut spec = fan_spec();
            edit(&mut spec);
            enclosure_case(&spec).is_err()
        };
        assert!(refused(&|s| s.mesh_level = 0));
        assert!(refused(&|s| s.mesh_level = 6));
        assert!(refused(&|s| s.size_mm[0] = -1.0));
        assert!(refused(&|s| s.components.clear()));
        assert!(refused(&|s| s.components[1].max_mm[0] = 500.0), "outside the box");
        assert!(refused(&|s| s.components[1].name = "pcb".into()), "duplicate name");
        assert!(refused(&|s| s.components[1].name = "top".into()), "patch name");
        assert!(refused(&|s| s.components[1].power_w = -1.0));
        assert!(refused(&|s| s.components[1].max_temperature_c = Some(20.0)), "limit below ambient");
        assert!(refused(&|s| s.components[1].min_mm[0] = 0.0), "touches the fan face");
        assert!(refused(&|s| s.wall_h_w_per_m2_k = 0.0));
        assert!(refused(&|s| if let Cooling::Fan { exhaust, .. } = &mut s.cooling {
            *exhaust = Face::MinusX
        }));
        assert!(
            refused(&|s| if let Cooling::Fan { curve, .. } = &mut s.cooling {
                curve[1].pressure_pa = 9.0
            }),
            "non-monotone curve"
        );
        assert!(refused(&|s| if let Cooling::Fan { opening_size_mm, .. } = &mut s.cooling {
            opening_size_mm[0] = 500.0
        }));
        assert!(refused(&|s| s.components[2].max_mm[0] = 200.0), "touches the filtered exhaust");
    }
}

#[test]
fn domain_cases_need_the_0_7_0_capability() {
    for template in templates().expect("templates") {
        let solver = template.case["solver"].as_str().expect("solver").to_owned();
        let names: Vec<_> = airflow::requirements(&template.case).iter().map(|r| r.capability).collect();
        let domain = airflow::DOMAIN_SOLVERS.contains(&solver.as_str());
        assert_eq!(names.contains(&"domain_cases_v1"), domain, "{}", template.id);
        if !domain {
            continue;
        }
        let mut manifest = evidence(&study()).manifest;
        manifest["solver"] = json!(solver);
        manifest["case"] = template.case.clone();
        manifest["capabilities"] = json!(["domain_cases_v1"]);
        manifest["rust_port_version"] = json!("0.6.0");
        assert!(validate_producer(&manifest, &template.case).is_err(), "{} needs 0.7.0", template.id);
        manifest["rust_port_version"] = json!(SOLVER_VERSION);
        manifest["capabilities"] = json!([]);
        assert!(validate_producer(&manifest, &template.case).is_err(), "{} needs the capability", template.id);
        manifest["capabilities"] = json!(["domain_cases_v1"]);
        validate_producer(&manifest, &template.case).unwrap_or_else(|e| panic!("{}: {e}", template.id));
    }
}

#[test]
fn sweeps_create_checked_variants_atomically() {
    let mut project = Project::default();
    let base = project.create("battery-discharge", "Cell").expect("create");
    let values = [json!(0.1), json!(0.25), json!(0.5)];
    let made = project.sweep(&base.id, "/protocol/0/current_a", &values).expect("sweep");
    assert_eq!(made.len(), 3);
    assert_eq!(project.studies.len(), 4);
    for (study, value) in made.iter().zip(&values) {
        assert_eq!(study.case["protocol"][0]["current_a"], *value);
        assert_eq!(study.profile, base.profile);
        assert!(study.name.contains("/protocol/0/current_a"));
        assert!(preflight(study).errors.is_empty());
    }
    assert_eq!(project.study(&base.id).expect("source").case, base.case, "source unchanged");
    let before = project.clone();
    for (pointer, values) in [
        ("/protocol/0/current_a", vec![json!("fast")]),
        ("/solver", vec![json!("spm")]),
        ("/no_such_key", vec![json!(1)]),
        ("protocol", vec![json!(1)]),
        ("/protocol/0/current_a", vec![]),
        ("/protocol/0/current_a", (0..17).map(|i| json!(i)).collect()),
        // The first value is valid; the second is not, so neither is committed.
        ("/dt_s", vec![json!(1.0), json!("x")]),
    ] {
        assert!(project.sweep(&base.id, pointer, &values).is_err(), "{pointer} {values:?}");
        assert_eq!(project, before, "a rejected sweep leaves the project unchanged");
    }
    assert!(project.sweep("mp-study-99", "/dt_s", &[json!(1.0)]).is_err());
}

#[test]
fn wizard_cases_become_studies_with_their_solver_profile() {
    let (case, _) = crate::wizard::enclosure_case(&enclosure_wizard::fan_spec()).expect("wizard");
    let mut project = Project::default();
    let study = project.create_from_case("Wizard enclosure", "Hot spot below 70 C?", case.clone()).expect("create");
    assert_eq!(study.case, case);
    assert_eq!(Some(study.profile.as_str()), profile_for_solver(case["solver"].as_str().expect("solver")));
    assert!(project.create_from_case("Bad", "", json!({"solver":"spray_fvm"})).is_err());
    assert_eq!(project.studies.len(), 1);
}

#[test]
fn cloned_projects_copy_the_case_but_not_its_evidence() {
    let mut project = Project::default();
    let base = project.create("battery-discharge", "Cell").expect("create");
    let copy = project.clone_study(&base.id, "").expect("clone");
    assert_ne!(copy.id, base.id);
    assert_eq!(copy.name, "Cell (copy)");
    assert_eq!(copy.case, base.case);
    assert_eq!(copy.profile, base.profile);
    assert!(copy.runs.is_empty());
    assert_eq!(project.studies.len(), 2);
    let named = project.clone_study(&base.id, "  Cell B ").expect("named clone");
    assert_eq!(named.name, "Cell B");
    let before = project.clone();
    assert!(project.clone_study("mp-study-missing", "x").is_err());
    assert_eq!(project, before, "a failed clone changes nothing");
}
