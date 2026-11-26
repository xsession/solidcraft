//! The benchmarks in `example_cases/benchmarks`: published problems with known answers, one or
//! more for each of the companion's domains, each a complete case and the checks its run must
//! pass. `make_benchmarks.py` there writes them and works out every expected value.
//!
//! Every benchmark must be a case SolidCraft accepts, with its source named (always run). With
//! `SOLIDCRAFT_TEST_CFD_SOLVE` set, the ignored test runs each through that companion and
//! judges its checks; `SOLIDCRAFT_TEST_OUTPUT_DIR` keeps what each measured, as
//! `benchmarks.json`.
use super::*;
use serde_json::{Value, json};
use solidcraft_engine::multiphysics::{Project, case_hash};

fn folder() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../example_cases/benchmarks")
}

/// Every benchmark, by file name, in order.
fn benchmarks() -> Vec<(String, Value)> {
    let mut files: Vec<PathBuf> = fs::read_dir(folder())
        .expect("example_cases/benchmarks")
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

/// A column of a CSV artifact, by its header.
fn column(text: &str, name: &str) -> Option<Vec<f64>> {
    let mut lines = text.lines();
    let index = lines.next()?.split(',').position(|h| h.trim() == name)?;
    lines.filter(|l| !l.trim().is_empty()).map(|l| l.split(',').nth(index).and_then(|v| v.trim().parse().ok())).collect()
}

/// A check's value from a run: a metric, a JSON artifact's value by path, or a CSV artifact's
/// column in a given row or interpolated at a value of another column; then over another such
/// value (`over`), then `add`, then `divide_by`.
fn value_of(spec: &Value, metrics: &Value, artifacts: &[Artifact]) -> Option<f64> {
    let artifact = |name: &str| artifacts.iter().find(|a| a.name == name).map(|a| a.text.as_str());
    let mut v = if let Some(metric) = spec["metric"].as_str() {
        metrics[metric].as_f64()?
    } else if let Some(file) = spec["json"].as_str() {
        let mut node: Value = serde_json::from_str(artifact(file)?).ok()?;
        for key in spec["path"].as_array()? {
            node = match key {
                Value::Number(i) => node.get(usize::try_from(i.as_u64()?).ok()?)?.clone(),
                Value::String(k) => node.get(k)?.clone(),
                _ => return None,
            };
        }
        node.as_f64()?
    } else {
        let text = artifact(spec["file"].as_str()?)?;
        let ys = column(text, spec["column"].as_str()?)?;
        if let Some(which) = spec["row"].as_str() {
            match which {
                "first" => *ys.first()?,
                "last" => *ys.last()?,
                "max" => ys.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                "min" => ys.iter().copied().fold(f64::INFINITY, f64::min),
                _ => return None,
            }
        } else {
            let xs = column(text, spec["at"]["column"].as_str()?)?;
            let x = spec["at"]["value"].as_f64()?;
            let mut pairs: Vec<(f64, f64)> = xs.into_iter().zip(ys).collect();
            pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
            pairs.windows(2).find(|w| w[0].0 <= x && x <= w[1].0).map(|w| {
                let ((x0, y0), (x1, y1)) = (w[0], w[1]);
                if x1 == x0 { y0 } else { y0 + (y1 - y0) * (x - x0) / (x1 - x0) }
            })?
        }
    };
    if spec.get("over").is_some() {
        v /= value_of(&spec["over"], metrics, artifacts)?;
    }
    if let Some(add) = spec["add"].as_f64() {
        v += add;
    }
    if let Some(divide) = spec["divide_by"].as_f64() {
        v /= divide;
    }
    Some(v)
}

/// Each benchmark names its source and what it shows, has checks with an expected value or a
/// bound, and is a case SolidCraft accepts as it stands.
#[test]
fn every_benchmark_is_an_accepted_case_with_its_source() {
    let all = benchmarks();
    assert!(all.len() >= 20, "{} benchmarks", all.len());
    for (name, benchmark) in &all {
        for key in ["title", "domain", "shows", "source"] {
            assert!(benchmark[key].as_str().is_some_and(|t| !t.is_empty()), "{name}: {key}");
        }
        let checks = benchmark["checks"].as_array().unwrap_or_else(|| panic!("{name}: checks"));
        assert!(!checks.is_empty(), "{name}: no checks");
        for check in checks {
            assert!(check["expect"].is_number() || check["at_least"].is_number() || check["at_most"].is_number(), "{name}: {check}");
        }
        let solver = benchmark["case"]["solver"].as_str().unwrap_or_else(|| panic!("{name}: the case names no solver"));
        assert!(solidcraft_engine::multiphysics::profile_for_solver(solver).is_some(), "{name}: no profile runs {solver}");
        if let Err(error) = solidcraft_engine::multiphysics::validate_case(&benchmark["case"]) {
            panic!("{name}: the case is refused: {error}");
        }
    }
}

/// Every benchmark run by a real companion, judged by its checks; all failures are reported
/// together, with what each benchmark measured.
#[test]
#[ignore = "requires SOLIDCRAFT_TEST_CFD_SOLVE pointing to a native Rust companion (0.18.0 or later)"]
fn real_benchmarks_meet_their_references() {
    let executable = std::env::var_os("SOLIDCRAFT_TEST_CFD_SOLVE").expect("set SOLIDCRAFT_TEST_CFD_SOLVE");
    let only = std::env::var("SOLIDCRAFT_BENCHMARKS").ok();
    let host = ProcessHost::new(PathBuf::from(executable));
    let mut failures = Vec::new();
    let mut measured = serde_json::Map::new();
    for (name, benchmark) in benchmarks() {
        if only.as_deref().is_some_and(|o| !o.split(',').any(|n| name.starts_with(n))) {
            continue;
        }
        // A run must belong to its study's profile: the profile of the benchmark's solver.
        let mut study = Project::default().create("poisson1d", &name).unwrap();
        let solver = benchmark["case"]["solver"].as_str().unwrap_or_default();
        study.profile =
            solidcraft_engine::multiphysics::profile_for_solver(solver).unwrap_or_else(|| panic!("{name}: no profile for {solver}")).into();
        study.case = benchmark["case"].clone();
        let request = FrozenRun { id: format!("run-{name}"), case_hash: case_hash(&study.case).unwrap(), study };
        let started = Instant::now();
        let mut job = host.start(request).unwrap();
        let deadline = Instant::now() + Duration::from_secs(30 * 60);
        let outcome = loop {
            if let Some(outcome) = job.poll() {
                break Some(outcome);
            }
            if Instant::now() > deadline {
                job.cancel();
                break None;
            }
            std::thread::sleep(Duration::from_millis(20));
        };
        let evidence = match outcome {
            Some(Outcome::Completed(evidence)) => evidence,
            Some(Outcome::Failed(error)) => {
                failures.push(format!("{name}: the run failed: {error}"));
                continue;
            }
            other => {
                failures.push(format!("{name}: the run did not finish ({})", if other.is_none() { "over 30 minutes" } else { "cancelled" }));
                continue;
            }
        };
        let metrics = evidence.manifest["metrics"].clone();
        let mut results = Vec::new();
        for check in benchmark["checks"].as_array().unwrap() {
            let what = check["what"].as_str().unwrap_or_default();
            match value_of(&check["value"], &metrics, &evidence.artifacts) {
                Some(value) => {
                    let verdict = super::example_cases::judge(check, value);
                    if let Err(why) = &verdict {
                        failures.push(format!("{name}: {what}: {why}"));
                    }
                    results.push(json!({"what": what, "value": value, "expect": check["expect"], "pass": verdict.is_ok()}));
                }
                None => failures.push(format!("{name}: {what}: the run gave no such value ({})", check["value"])),
            }
        }
        measured.insert(
            name.clone(),
            json!({"title": benchmark["title"], "domain": benchmark["domain"], "seconds": started.elapsed().as_secs_f64(), "checks": results}),
        );
    }
    if let Some(root) = std::env::var_os("SOLIDCRAFT_TEST_OUTPUT_DIR") {
        fs::create_dir_all(&root).unwrap();
        fs::write(PathBuf::from(root).join("benchmarks.json"), serde_json::to_vec_pretty(&Value::Object(measured)).unwrap()).unwrap();
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn values_come_from_metrics_json_and_csv_artifacts() {
    let artifacts = vec![
        Artifact { name: "a.csv".into(), text: "x,y\n0,1\n2,5\n".into() },
        Artifact { name: "c.json".into(), text: r#"{"final": [0.25, 0.75]}"#.into() },
    ];
    let metrics = json!({"m": 10.0});
    let ask = |spec: Value| value_of(&spec, &metrics, &artifacts);
    assert_eq!(ask(json!({"metric": "m", "add": -4.0, "divide_by": 2.0})), Some(3.0));
    assert_eq!(ask(json!({"json": "c.json", "path": ["final", 1]})), Some(0.75));
    assert_eq!(ask(json!({"file": "a.csv", "column": "y", "at": {"column": "x", "value": 1.0}})), Some(3.0));
    assert_eq!(ask(json!({"file": "a.csv", "column": "y", "row": "max"})), Some(5.0));
    assert_eq!(ask(json!({"file": "a.csv", "column": "y", "at": {"column": "x", "value": 3.0}})), None);
    assert_eq!(ask(json!({"metric": "missing"})), None);
    let ratio = json!({"file": "a.csv", "column": "y", "row": "max", "over": {"metric": "m"}});
    assert_eq!(ask(ratio), Some(0.5));
}
