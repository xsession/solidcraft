use super::*;
use serde_json::json;
use solidcraft_engine::multiphysics::{Project, case_hash, templates};

fn request() -> FrozenRun {
    let study = Project::default().create("poisson1d", "Host fixture").unwrap();
    FrozenRun { id: "run-host-fixture".into(), case_hash: case_hash(&study.case).unwrap(), study }
}

fn manifest(request: &FrozenRun) -> Value {
    json!({
        "schema":"cfd-solvers-rs/run/1",
        "rust_port_version":"0.2.0",
        "upstream_revision":solidcraft_engine::multiphysics::SOLVER_REVISION,
        "solver":request.study.case.get("solver").and_then(Value::as_str).unwrap(),
        "backend":"native_rust_serial_cpu",
        "complete_clone":false,
        "execution_completed":true,
        "linear_solver_converged":true,
        "steady_state_converged":null,
        "metrics":{"residual":1e-12},
        "files":["field.csv","manifest.json"],
        "case":request.study.case,
    })
}

fn await_outcome(job: &mut dyn Job) -> Outcome {
    await_outcome_for(job, Duration::from_secs(20))
}

fn await_outcome_for(job: &mut dyn Job, timeout: Duration) -> Outcome {
    // The deadline bounds a broken test; no assertion depends on how long a run takes.
    let deadline = Instant::now() + timeout;
    loop {
        if let Some(outcome) = job.poll() {
            return outcome;
        }
        assert!(Instant::now() < deadline, "solver job never reached a terminal state");
        std::thread::sleep(Duration::from_millis(5));
    }
}

fn preserve_evidence(name: &str, evidence: &RunEvidence) {
    let Some(root) = std::env::var_os("SOLIDCRAFT_TEST_OUTPUT_DIR") else { return };
    let root = PathBuf::from(root);
    fs::create_dir_all(&root).unwrap();
    let directory = root.join(name);
    fs::create_dir(&directory).expect("choose a fresh SOLIDCRAFT_TEST_OUTPUT_DIR for each verification run");
    for artifact in &evidence.artifacts {
        fs::write(directory.join(&artifact.name), &artifact.text).unwrap();
    }
    fs::write(directory.join("case.json"), serde_json::to_vec_pretty(&evidence.input_case).unwrap()).unwrap();
}

#[test]
fn artifact_paths_reject_cross_platform_escapes() {
    for path in ["", "../outside", "a/../../outside", "/absolute", "C:\\absolute", "\\\\server\\file", "a\\b", "file\0x", ".", "a/./b"] {
        assert!(!ordinary_relative(path), "{path:?}");
    }
    assert!(ordinary_relative("field.vtu"));
    assert!(ordinary_relative("frames/field-000.vtu"));
}

#[test]
fn bounded_reads_reject_excessive_artifacts() {
    let directory = WorkDirectory::create().unwrap();
    fs::write(directory.0.join("field.csv"), b"123456789").unwrap();
    let root = directory.0.canonicalize().unwrap();
    assert!(read_file(&root, "field.csv", 8).unwrap_err().contains("bounded regular file"));
    assert_eq!(read_file(&root, "field.csv", 9).unwrap(), "123456789");
}

#[test]
fn output_growth_is_refused_before_admission() {
    let directory = WorkDirectory::create().unwrap();
    let output = directory.0.join("results");
    fs::create_dir(&output).unwrap();
    File::create(output.join("large.vtu")).unwrap().set_len(MAX_ARTIFACT_BYTES as u64 + 1).unwrap();
    assert!(output_budget(&output).unwrap_err().contains("16 MiB"));
}

#[cfg(unix)]
fn quote(path: &Path) -> String {
    format!("'{}'", path.to_string_lossy().replace('\'', "'\\''"))
}

#[cfg(unix)]
struct FakeCompanion {
    directory: WorkDirectory,
    executable: PathBuf,
}

#[cfg(unix)]
impl FakeCompanion {
    fn new(request: &FrozenRun, manifest: Value, body: &str) -> Self {
        use std::os::unix::fs::PermissionsExt as _;
        let directory = WorkDirectory::create().unwrap();
        let source = directory.0.join("manifest-source.json");
        fs::write(&source, serde_json::to_vec(&manifest).unwrap()).unwrap();
        let executable = directory.0.join("companion");
        let script = format!(
            "#!/bin/sh\nset -eu\n[ \"$1\" = run ]\n[ \"$3\" = --out ]\nmkdir \"$4\"\ncp {} \"$4/manifest.json\"\nprintf 'x,value\\n0,1\\n1,2\\n' > \"$4/field.csv\"\n{}\n",
            quote(&source),
            body,
        );
        fs::write(&executable, script).unwrap();
        fs::set_permissions(&executable, fs::Permissions::from_mode(0o700)).unwrap();
        // Keep the request borrow explicit: a fixture always agrees with its frozen case.
        assert_eq!(request.case_hash, case_hash(&request.study.case).unwrap());
        Self { directory, executable }
    }

    fn start(&self, request: FrozenRun) -> Box<dyn Job> {
        ProcessHost::new(self.executable.clone()).start(request).unwrap()
    }
}

#[test]
#[cfg(unix)]
fn process_result_admission_preserves_frozen_case_and_artifacts() {
    let request = request();
    let companion = FakeCompanion::new(&request, manifest(&request), "");
    let mut job = companion.start(request.clone());
    match await_outcome(job.as_mut()) {
        Outcome::Completed(evidence) => {
            assert_eq!(evidence.id, request.id);
            assert_eq!(evidence.case_hash, request.case_hash);
            assert_eq!(evidence.study_revision, request.study.revision);
            assert_eq!(evidence.artifacts.len(), 2);
            assert_eq!(evidence.artifacts.first().unwrap().text, "x,value\n0,1\n1,2\n");
        }
        _ => panic!("a conforming completed process was refused"),
    }
    assert!(job.poll().is_none(), "terminal results are delivered once");
}

#[test]
#[cfg(unix)]
fn authored_case_larger_than_one_mib_reaches_worker_and_retains_manifest() {
    let mut request = request();
    request.study.case["extension"] = json!("x".repeat((1 << 20) + 1024));
    request.case_hash = case_hash(&request.study.case).unwrap();
    let companion = FakeCompanion::new(&request, manifest(&request), "");
    let mut job = companion.start(request.clone());
    match await_outcome(job.as_mut()) {
        Outcome::Completed(evidence) => assert_eq!(evidence.input_case, request.study.case),
        Outcome::Failed(error) => panic!("bounded case was refused: {error}"),
        Outcome::Cancelled => panic!("bounded case was cancelled"),
    }
}

#[test]
#[cfg(unix)]
fn incompatible_companion_version_never_becomes_a_result() {
    let request = request();
    let mut value = manifest(&request);
    value["rust_port_version"] = json!("0.1.0");
    let companion = FakeCompanion::new(&request, value, "");
    let mut job = companion.start(request);
    assert!(matches!(await_outcome(job.as_mut()), Outcome::Failed(message) if message.contains("rust_port_version")));
}

#[test]
#[cfg(unix)]
fn current_companion_accepts_legacy_case() {
    let request = request();
    let mut value = manifest(&request);
    value["rust_port_version"] = json!("0.3.0");
    value["capabilities"] = json!(["grid_obstacles_v1"]);
    let companion = FakeCompanion::new(&request, value, "");
    let mut job = companion.start(request);
    assert!(matches!(await_outcome(job.as_mut()), Outcome::Completed(_)));
}

#[test]
#[cfg(unix)]
fn obstacle_case_requires_current_companion_and_capability() {
    let study = Project::default().create("flow-electronics-enclosure", "Enclosure admission").unwrap();
    let request = FrozenRun { id: "run-enclosure-admission".into(), case_hash: case_hash(&study.case).unwrap(), study };
    for (version, capabilities) in [("0.2.0", json!(["grid_obstacles_v1"])), ("0.3.0", json!([]))] {
        let mut value = manifest(&request);
        value["rust_port_version"] = json!(version);
        value["capabilities"] = capabilities;
        let companion = FakeCompanion::new(&request, value, "");
        let mut job = companion.start(request.clone());
        assert!(matches!(await_outcome(job.as_mut()), Outcome::Failed(_)), "version {version}");
    }
}

#[test]
#[cfg(unix)]
fn fan_and_porous_case_requires_0_4_0_and_each_capability() {
    let study = Project::default().create("flow-electronics-enclosure-fan", "Fan admission").unwrap();
    let request = FrozenRun { id: "run-fan-admission".into(), case_hash: case_hash(&study.case).unwrap(), study };
    let all = json!(["grid_obstacles_v1", "fan_curve_v1", "porous_blocks_v1"]);
    for (version, capabilities) in
        [("0.3.0", all.clone()), ("0.4.0", json!(["grid_obstacles_v1", "porous_blocks_v1"])), ("0.4.0", json!(["grid_obstacles_v1", "fan_curve_v1"]))]
    {
        let mut value = manifest(&request);
        value["rust_port_version"] = json!(version);
        value["capabilities"] = capabilities.clone();
        let companion = FakeCompanion::new(&request, value, "");
        let mut job = companion.start(request.clone());
        assert!(matches!(await_outcome(job.as_mut()), Outcome::Failed(_)), "version {version} {capabilities}");
    }
    let mut value = manifest(&request);
    value["rust_port_version"] = json!("0.4.0");
    value["capabilities"] = all;
    let companion = FakeCompanion::new(&request, value, "");
    let mut job = companion.start(request);
    assert!(matches!(await_outcome(job.as_mut()), Outcome::Completed(_)));
}

#[test]
#[cfg(unix)]
fn manifest_case_substitution_never_becomes_a_result() {
    let request = request();
    let mut value = manifest(&request);
    value["case"]["solver"] = json!("poisson3d");
    let companion = FakeCompanion::new(&request, value, "");
    let mut job = companion.start(request);
    assert!(matches!(await_outcome(job.as_mut()), Outcome::Failed(_)));
}

#[test]
#[cfg(unix)]
fn nonzero_exit_never_admits_even_a_complete_manifest() {
    let request = request();
    let companion = FakeCompanion::new(&request, manifest(&request), "exit 7");
    let mut job = companion.start(request);
    assert!(matches!(await_outcome(job.as_mut()), Outcome::Failed(message) if message.contains("failed with")));
}

#[test]
#[cfg(unix)]
fn artifact_parent_escape_never_reads_outside_the_run() {
    let request = request();
    let mut value = manifest(&request);
    value["files"] = json!(["../case.json"]);
    let companion = FakeCompanion::new(&request, value, "");
    let mut job = companion.start(request);
    assert!(matches!(await_outcome(job.as_mut()), Outcome::Failed(message) if message.contains("relative path")));
}

#[test]
#[cfg(unix)]
fn symlink_artifact_escape_never_reads_outside_the_run() {
    let request = request();
    let companion = FakeCompanion::new(&request, manifest(&request), "rm \"$4/field.csv\"\nln -s \"$2\" \"$4/field.csv\"");
    let mut job = companion.start(request);
    assert!(matches!(await_outcome(job.as_mut()), Outcome::Failed(message) if message.contains("leaves its task directory")));
}

#[test]
#[cfg(unix)]
fn cancellation_waits_for_a_worker_acknowledgement() {
    let request = request();
    let companion = FakeCompanion::new(&request, manifest(&request), "");
    let marker = companion.directory.0.join("started");
    use std::os::unix::fs::PermissionsExt as _;
    let script = format!("#!/bin/sh\nprintf started > {}\nwhile :; do :; done\n", quote(&marker));
    fs::write(&companion.executable, script).unwrap();
    fs::set_permissions(&companion.executable, fs::Permissions::from_mode(0o700)).unwrap();
    let mut job = companion.start(request);
    let deadline = Instant::now() + Duration::from_secs(20);
    while !marker.exists() {
        assert!(Instant::now() < deadline, "fake companion never started");
        std::thread::sleep(Duration::from_millis(5));
    }
    job.cancel();
    assert!(matches!(await_outcome(job.as_mut()), Outcome::Cancelled));
}

#[test]
#[ignore = "requires SOLIDCRAFT_TEST_CFD_SOLVE pointing to a compatible native Rust companion"]
fn real_reference_runner_round_trip() {
    let executable = std::env::var_os("SOLIDCRAFT_TEST_CFD_SOLVE").expect("set SOLIDCRAFT_TEST_CFD_SOLVE");
    let request = request();
    let mut job = ProcessHost::new(PathBuf::from(executable)).start(request).unwrap();
    match await_outcome(job.as_mut()) {
        Outcome::Completed(evidence) => {
            assert!(!evidence.artifacts.is_empty());
            assert_eq!(evidence.manifest["execution_completed"], true);
        }
        Outcome::Failed(error) => panic!("real companion was refused: {error}"),
        Outcome::Cancelled => panic!("real companion was cancelled"),
    }
}

#[test]
#[ignore = "requires SOLIDCRAFT_TEST_CFD_SOLVE pointing to a 0.4.0 native companion"]
fn real_electronics_enclosure_retains_fluid_only_fields() {
    use solidcraft_engine::multiphysics::vtk::{Association, parse_pvd, parse_vtu};

    let executable = std::env::var_os("SOLIDCRAFT_TEST_CFD_SOLVE").expect("set SOLIDCRAFT_TEST_CFD_SOLVE");
    let study = Project::default().create("flow-electronics-enclosure", "Prescribed-flow isothermal enclosure").unwrap();
    let input = study.case.clone();
    let request = FrozenRun { id: "run-enclosure-reference".into(), case_hash: case_hash(&input).unwrap(), study };
    let mut job = ProcessHost::new(PathBuf::from(executable)).start(request).unwrap();
    // The deadline detects a broken test/companion; acceptance never depends on speed.
    let deadline = Instant::now() + Duration::from_secs(300);
    let outcome = loop {
        if let Some(outcome) = job.poll() {
            break outcome;
        }
        assert!(Instant::now() < deadline, "enclosure companion never reached a terminal state");
        std::thread::sleep(Duration::from_millis(20));
    };
    let evidence = match outcome {
        Outcome::Completed(evidence) => evidence,
        Outcome::Failed(error) => panic!("enclosure companion refused: {error}"),
        Outcome::Cancelled => panic!("enclosure companion was cancelled"),
    };
    assert_eq!(evidence.input_case, input);
    assert_eq!(evidence.manifest["rust_port_version"], solidcraft_engine::multiphysics::SOLVER_VERSION);
    assert!(evidence.manifest["capabilities"].as_array().unwrap().contains(&json!("grid_obstacles_v1")));
    assert_eq!(evidence.manifest["linear_solver_converged"], true);
    let pvd = evidence.artifacts.iter().find(|artifact| artifact.name.ends_with(".pvd")).unwrap();
    let frames = parse_pvd(&pvd.text).unwrap();
    let frame = frames.last().unwrap();
    assert!(frame.coordinate > 0.0);
    let vtu = evidence.artifacts.iter().find(|artifact| artifact.name == frame.file).unwrap();
    let grid = parse_vtu(&vtu.text).unwrap();
    let dimensions: Vec<usize> = input["mesh"]["dimensions"].as_array().unwrap().iter().map(|value| value.as_u64().unwrap() as usize).collect();
    let total_cells: usize = dimensions.iter().product();
    let obstacles = input["obstacles"].as_array().unwrap();
    let removed: usize = obstacles
        .iter()
        .map(|obstacle| {
            (0..3).map(|axis| (obstacle["max"][axis].as_u64().unwrap() - obstacle["min"][axis].as_u64().unwrap()) as usize).product::<usize>()
        })
        .sum();
    assert_eq!(grid.cell_types.len(), total_cells - removed);
    assert!(grid.cell_types.iter().all(|kind| *kind == 12));
    let velocity = grid.fields.iter().find(|field| field.name == "velocity_m_per_s").unwrap();
    assert_eq!(velocity.association, Association::Cell);
    assert_eq!(velocity.components, 3);
    assert_eq!(velocity.values.len(), grid.cell_types.len() * 3);
    assert!(velocity.values.iter().all(|value| value.is_finite()));
    assert!(velocity.values.iter().any(|value| value.abs() > 1e-8));
    let lengths: Vec<f64> = ["x", "y", "z"].iter().map(|axis| input["mesh"]["lengths"][axis].as_f64().unwrap()).collect();
    let mut start = 0;
    for end in &grid.offsets {
        let cell_points: Vec<_> = grid.connectivity[start..*end].iter().map(|index| grid.points[*index]).collect();
        let center: Vec<f64> = (0..3).map(|axis| cell_points.iter().map(|point| point[axis]).sum::<f64>() / cell_points.len() as f64).collect();
        for obstacle in obstacles {
            let inside = (0..3).all(|axis| {
                let lo = obstacle["min"][axis].as_u64().unwrap() as f64 * lengths[axis] / dimensions[axis] as f64;
                let hi = obstacle["max"][axis].as_u64().unwrap() as f64 * lengths[axis] / dimensions[axis] as f64;
                center[axis] > lo && center[axis] < hi
            });
            assert!(!inside, "fluid export cell occupies a solid obstacle");
        }
        start = *end;
    }
    for name in ["continuity_l2", "continuity_max"] {
        let value = evidence.manifest["metrics"][name].as_f64().unwrap();
        assert!(value.is_finite() && value >= 0.0);
    }
    assert_eq!(evidence.manifest["metrics"]["fluid_cell_count"].as_f64(), Some(grid.cell_types.len() as f64));
    assert_eq!(evidence.manifest["metrics"]["excluded_solid_cell_count"].as_f64(), Some(removed as f64));
    for name in ["boundary_inflow_m3_per_s", "boundary_outflow_m3_per_s"] {
        assert!(evidence.manifest["metrics"][name].as_f64().unwrap() > 0.0);
    }
    let imbalance = evidence.manifest["metrics"]["relative_boundary_flux_imbalance"].as_f64().unwrap();
    assert!(imbalance.is_finite() && imbalance < 1e-4, "relative flux imbalance {imbalance}");
    let report = evidence.artifacts.iter().find(|artifact| artifact.name == "mesh_report.json").unwrap();
    let report: Value = serde_json::from_str(&report.text).unwrap();
    assert_eq!(report["schema"], "cfd-solvers-rs/grid-obstacle-mesh/1");
    for obstacle in obstacles {
        assert!(report["patches"].as_array().unwrap().iter().any(|patch| patch["name"] == obstacle["name"]));
    }
    preserve_evidence("flow-electronics-enclosure", &evidence);
    let temporary = WorkDirectory::create().unwrap();
    let root = if let Some(output) = std::env::var_os("SOLIDCRAFT_TEST_OUTPUT_DIR") {
        let root = PathBuf::from(output).join("enclosure-document");
        fs::create_dir(&root).unwrap();
        root
    } else {
        temporary.0.clone()
    };
    let mut session = solidcraft_engine::Session::default();
    let created = session
        .execute("multiphysics.study.create", &json!({"template":"flow-electronics-enclosure","name":"Electronics enclosure airflow"}))
        .unwrap();
    let id = created["study"]["id"].as_str().unwrap().to_string();
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.attach_evidence(&id, evidence.clone()).unwrap();
    session.doc_mut().multiphysics = Arc::new(project);
    let design = root.join("electronics-enclosure-airflow.sc");
    session.execute("file.save_as", &json!({"path":design})).unwrap();
    let mut reopened = solidcraft_engine::Session::default();
    reopened.execute("doc.open", &json!({"path":design})).unwrap();
    let saved = reopened.doc.multiphysics.study(&id).unwrap();
    assert_eq!(saved.case, input);
    assert_eq!(saved.runs.first().unwrap(), &evidence);
    let inspected = reopened.execute("multiphysics.result.inspect", &json!({"id":id,"run_id":evidence.id,"artifact":frame.file})).unwrap();
    assert_eq!(inspected["content"]["cells"].as_u64(), Some(grid.cell_types.len() as u64));
    let exported = root.join("exported-results");
    reopened.execute("multiphysics.result.export", &json!({"id":id,"run_id":evidence.id,"directory":exported})).unwrap();
    assert_eq!(fs::read_to_string(exported.join(&frame.file)).unwrap(), vtu.text);
}

#[test]
#[ignore = "requires SOLIDCRAFT_TEST_CFD_SOLVE pointing to a native Rust companion with all template capabilities"]
fn real_reference_runner_all_templates() {
    let executable = std::env::var_os("SOLIDCRAFT_TEST_CFD_SOLVE").expect("set SOLIDCRAFT_TEST_CFD_SOLVE");
    let host = ProcessHost::new(PathBuf::from(executable));
    let templates = templates().unwrap();
    assert!(templates.iter().any(|template| template.id == "flow-electronics-enclosure"));
    for template in templates {
        let study = Project::default().create(&template.id, "Reference round trip").unwrap();
        let request = FrozenRun { id: format!("run-{}", template.id), case_hash: case_hash(&study.case).unwrap(), study };
        let mut job = host.start(request).unwrap();
        match await_outcome_for(job.as_mut(), Duration::from_secs(300)) {
            Outcome::Completed(evidence) => {
                assert!(!evidence.artifacts.is_empty(), "{}", template.id);
                assert_eq!(evidence.manifest["execution_completed"], true, "{}", template.id);
                preserve_evidence(&template.id, &evidence);
            }
            Outcome::Failed(error) => panic!("{} companion result refused: {error}", template.id),
            Outcome::Cancelled => panic!("{} companion was cancelled", template.id),
        }
    }
}

#[test]
#[ignore = "requires SOLIDCRAFT_TEST_CFD_SOLVE pointing to a compatible native Rust companion"]
fn real_engine_run_keep_reopen_and_export() {
    use solidcraft_engine::Session;

    let executable = std::env::var_os("SOLIDCRAFT_TEST_CFD_SOLVE").expect("set SOLIDCRAFT_TEST_CFD_SOLVE");
    let temporary = WorkDirectory::create().unwrap();
    let root = if let Some(output) = std::env::var_os("SOLIDCRAFT_TEST_OUTPUT_DIR") {
        let root = PathBuf::from(output);
        fs::create_dir_all(&root).unwrap();
        let root = root.join("engine-workflow");
        fs::create_dir(&root).expect("choose a fresh output directory for each workflow verification");
        root
    } else {
        temporary.0.clone()
    };
    let mut session = Session::default();
    session.set_multiphysics_host(Some(Arc::new(ProcessHost::new(PathBuf::from(executable)))));
    let created = session.execute("multiphysics.study.create", &json!({"template":"poisson1d","name":"Analytical potential"})).unwrap();
    let study_id = created["study"]["id"].as_str().unwrap().to_string();
    let original_case = created["study"]["case"].clone();
    assert_eq!(original_case["source"].as_f64(), Some(2.0));
    let started = session.execute("multiphysics.run.start", &json!({"id":study_id})).unwrap();
    let job_id = started["id"].as_str().unwrap().to_string();
    let deadline = Instant::now() + Duration::from_secs(20);
    loop {
        let status = session.execute("multiphysics.run.status", &json!({})).unwrap();
        let job = status["jobs"].as_array().unwrap().iter().find(|job| job["id"] == job_id).unwrap();
        match job["state"].as_str().unwrap() {
            "running" | "cancelling" => {}
            "succeeded" => break,
            other => panic!("real engine job reached {other}: {}", job["message"]),
        }
        assert!(Instant::now() < deadline, "engine command workflow never completed");
        std::thread::sleep(Duration::from_millis(5));
    }
    session.execute("multiphysics.run.keep", &json!({"job_id":job_id})).unwrap();
    let design = root.join("analytical-potential.sc");
    session.execute("file.save_as", &json!({"path":design.to_string_lossy()})).unwrap();
    let mut reopened = Session::default();
    reopened.execute("doc.open", &json!({"path":design.to_string_lossy()})).unwrap();
    let saved = reopened.doc.multiphysics.study(&study_id).unwrap();
    let evidence = saved.runs.first().unwrap();
    assert_eq!(evidence.input_case, original_case);
    let csv = evidence.artifacts.iter().find(|artifact| artifact.name == "solution.csv").unwrap();
    let samples: Vec<(f64, f64)> = csv
        .text
        .lines()
        .skip(1)
        .map(|line| {
            let (x, value) = line.split_once(',').unwrap();
            (x.parse().unwrap(), value.parse().unwrap())
        })
        .collect();
    assert_eq!(samples.len(), 65);
    for (x, value) in &samples {
        assert!((value - x * (1.0 - x)).abs() < 1e-10, "solution violates -u''=2 at x={x}");
    }
    let maximum = samples.iter().map(|(_, value)| *value).fold(f64::NEG_INFINITY, f64::max);
    assert!((maximum - 0.25).abs() < 1e-10);
    let retained_manifest = evidence.manifest.clone();
    let export = root.join("paraview");
    reopened.execute("multiphysics.result.export", &json!({"id":study_id,"run_id":job_id,"directory":export.to_string_lossy()})).unwrap();
    let exported_case: Value = serde_json::from_slice(&fs::read(export.join("case.json")).unwrap()).unwrap();
    let exported_manifest: Value = serde_json::from_slice(&fs::read(export.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(exported_case, original_case);
    assert_eq!(exported_manifest, retained_manifest);
    assert!(export.join("solution.vtu").is_file());
    assert!(export.join("solidcraft-study.json").is_file());
}
