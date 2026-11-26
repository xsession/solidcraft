//! The example cases in `example_cases/model`: each a list of SolidCraft commands that builds a
//! study of the model, the solver capabilities it shows, and the checks its run must pass.
//!
//! Every example is built in a headless session and must ask the solver for the capabilities it
//! says it shows (always run). With `SOLIDCRAFT_TEST_CFD_SOLVE` set, the ignored test runs each
//! through that companion and judges its checks; `SOLIDCRAFT_TEST_OUTPUT_DIR` keeps what each
//! measured, as `model-examples.json`.
use super::*;
use serde_json::{Value, json};
use solidcraft_engine::Session;

fn folder() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../example_cases/model")
}

/// Every example, by file name, in order.
fn examples() -> Vec<(String, Value)> {
    let mut files: Vec<PathBuf> = fs::read_dir(folder())
        .expect("example_cases/model")
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|path| path.extension().is_some_and(|e| e == "json"))
        .collect();
    files.sort();
    files
        .into_iter()
        .map(|path| {
            let text = fs::read_to_string(&path).unwrap();
            let name = path.file_stem().unwrap().to_string_lossy().into_owned();
            (name, serde_json::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display())))
        })
        .collect()
}

/// `params` with `"$study"` as the study's id and `"$revision"` as its revision now.
fn fill(params: &Value, study: &str, revision: u64) -> Value {
    match params {
        Value::String(s) if s == "$study" => json!(study),
        Value::String(s) if s == "$revision" => json!(revision),
        Value::Array(items) => Value::Array(items.iter().map(|v| fill(v, study, revision)).collect()),
        Value::Object(map) => Value::Object(map.iter().map(|(k, v)| (k.clone(), fill(v, study, revision))).collect()),
        other => other.clone(),
    }
}

/// Run an example's steps in `session`; the study the steps made.
fn build(session: &mut Session, name: &str, example: &Value) -> String {
    let mut study = String::new();
    for (index, step) in example["steps"].as_array().expect("steps").iter().enumerate() {
        let command = step["command"].as_str().expect("command");
        let revision = session.doc.multiphysics.study(&study).map_or(0, |s| s.revision);
        let params = fill(&step["params"], &study, revision);
        let reply = session.execute(command, &params).unwrap_or_else(|e| panic!("{name} step {} ({command}): {e}", index + 1));
        if let Some(id) = reply["study"]["id"].as_str() {
            study = id.to_string();
        }
    }
    assert!(!study.is_empty(), "{name}: no step made a study");
    study
}

/// A check's value from a run's metrics: `metric`, less `minus`, over `divide`, its size with `abs`.
fn value_of(check: &Value, metrics: &Value) -> Option<f64> {
    let spec = &check["value"];
    let get = |key: &str| spec[key].as_str().map(|m| metrics[m].as_f64());
    let mut v = get("metric")??;
    if let Some(minus) = get("minus") {
        v -= minus?;
    }
    if let Some(divide) = get("divide") {
        v /= divide?;
    }
    if spec["abs"].as_bool() == Some(true) {
        v = v.abs();
    }
    Some(v)
}

/// Whether `value` passes `check`: within `tolerance` (a share) of `expect`, and between
/// `at_least` and `at_most`.
pub(super) fn judge(check: &Value, value: f64) -> Result<(), String> {
    if let Some(expect) = check["expect"].as_f64() {
        let tolerance = check["tolerance"].as_f64().unwrap_or(0.05);
        if !((value - expect).abs() <= tolerance * expect.abs().max(f64::MIN_POSITIVE)) {
            return Err(format!("{value:e} is not within {:.1} % of {expect:e}", tolerance * 100.0));
        }
    }
    if let Some(least) = check["at_least"].as_f64()
        && !(value >= least)
    {
        return Err(format!("{value:e} is below {least:e}"));
    }
    if let Some(most) = check["at_most"].as_f64()
        && !(value <= most)
    {
        return Err(format!("{value:e} is above {most:e}"));
    }
    Ok(())
}

/// Each example builds from its commands alone, and the case it builds asks the solver for
/// every capability it says it shows: the examples stay true as SolidCraft changes.
#[test]
fn every_model_example_builds_and_asks_for_what_it_shows() {
    let all = examples();
    assert!(all.len() >= 11, "{} examples", all.len());
    for (name, example) in &all {
        assert!(example["title"].as_str().is_some_and(|t| !t.is_empty()) && example["shows"].as_str().is_some(), "{name}: title and shows");
        assert!(example["checks"].as_array().is_some_and(|c| !c.is_empty()), "{name}: no checks");
        let mut session = Session::default();
        let study = build(&mut session, name, example);
        let case = session.doc.multiphysics.study(&study).map(|s| s.case.clone()).unwrap();
        let needs: Vec<&str> = solidcraft_engine::multiphysics::airflow::requirements(&case).iter().map(|r| r.capability).collect();
        for capability in example["capabilities"].as_array().expect("capabilities") {
            let capability = capability.as_str().unwrap();
            assert!(needs.contains(&capability), "{name} says it shows {capability}, but its case asks for {needs:?}");
        }
        assert!(solidcraft_engine::multiphysics::validate_case(&case).is_ok(), "{name}: the case is refused");
    }
}

/// Every example run by a real companion, judged by its checks; all failures are reported
/// together, with what each example measured.
#[test]
#[ignore = "requires SOLIDCRAFT_TEST_CFD_SOLVE pointing to a native Rust companion (0.18.0 or later)"]
fn real_model_examples_meet_their_checks() {
    let executable = std::env::var_os("SOLIDCRAFT_TEST_CFD_SOLVE").expect("set SOLIDCRAFT_TEST_CFD_SOLVE");
    let only = std::env::var("SOLIDCRAFT_EXAMPLES").ok();
    let mut failures = Vec::new();
    let mut measured = serde_json::Map::new();
    for (name, example) in examples() {
        if only.as_deref().is_some_and(|o| !o.split(',').any(|n| name.starts_with(n))) {
            continue;
        }
        let mut session = Session::default();
        session.set_multiphysics_host(Some(Arc::new(ProcessHost::new(PathBuf::from(&executable)))));
        let study = build(&mut session, &name, &example);
        let started = Instant::now();
        let job = session.execute("multiphysics.run.start", &json!({"id": study})).unwrap();
        let job_id = job["id"].as_str().unwrap().to_string();
        let deadline = Instant::now() + Duration::from_secs(2 * 60 * 60);
        let state = loop {
            let status = session.execute("multiphysics.run.status", &json!({})).unwrap();
            let entry = status["jobs"].as_array().unwrap().iter().find(|j| j["id"] == job_id).cloned().unwrap();
            match entry["state"].as_str().unwrap_or_default() {
                "running" | "cancelling" | "queued" => {}
                done => break (done.to_string(), entry["message"].clone()),
            }
            assert!(Instant::now() < deadline, "{name}: the run never ended");
            std::thread::sleep(Duration::from_millis(200));
        };
        if state.0 != "succeeded" {
            failures.push(format!("{name}: the run {}: {}", state.0, state.1));
            continue;
        }
        session.execute("multiphysics.run.keep", &json!({"job_id": job_id})).unwrap();
        let run = session.doc.multiphysics.study(&study).and_then(|s| s.runs.last().cloned()).unwrap();
        let metrics = run.manifest["metrics"].clone();
        let reported: Vec<&str> = run.manifest["capabilities"].as_array().into_iter().flatten().filter_map(Value::as_str).collect();
        let mut results = Vec::new();
        for capability in example["capabilities"].as_array().unwrap().iter().filter_map(Value::as_str) {
            if !reported.contains(&capability) {
                failures.push(format!("{name}: the companion did not report {capability}"));
            }
        }
        for check in example["checks"].as_array().unwrap() {
            let what = check["what"].as_str().unwrap_or_default();
            match value_of(check, &metrics) {
                Some(value) => {
                    let verdict = judge(check, value);
                    if let Err(why) = &verdict {
                        failures.push(format!("{name}: {what}: {why}"));
                    }
                    results.push(json!({"what": what, "value": value, "expect": check["expect"], "pass": verdict.is_ok()}));
                }
                None => failures.push(format!("{name}: {what}: the run reported no such metric ({})", check["value"])),
            }
        }
        measured.insert(name.clone(), json!({"title": example["title"], "seconds": started.elapsed().as_secs_f64(), "checks": results}));
    }
    if let Some(root) = std::env::var_os("SOLIDCRAFT_TEST_OUTPUT_DIR") {
        fs::create_dir_all(&root).unwrap();
        fs::write(PathBuf::from(root).join("model-examples.json"), serde_json::to_vec_pretty(&Value::Object(measured)).unwrap()).unwrap();
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn checks_take_differences_ratios_and_bounds() {
    let metrics = json!({"a": 5.0, "b": 2.0, "q": -0.5});
    let check = json!({"value": {"metric": "a", "minus": "b", "divide": "q", "abs": true}, "expect": 6.0, "tolerance": 0.01});
    assert_eq!(value_of(&check, &metrics), Some(6.0));
    assert!(judge(&check, 6.0).is_ok() && judge(&check, 6.1).is_err());
    assert!(judge(&json!({"at_least": 1.0, "at_most": 2.0}), 1.5).is_ok());
    assert!(judge(&json!({"at_least": 1.0}), 0.5).is_err());
    assert!(judge(&json!({"at_most": 1.0}), f64::NAN).is_err());
    assert_eq!(value_of(&json!({"value": {"metric": "missing"}}), &metrics), None);
    assert_eq!(fill(&json!({"id": "$study", "n": ["$revision"]}), "s1", 3), json!({"id": "s1", "n": [3]}));
}
