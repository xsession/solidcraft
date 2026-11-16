//! Analytic and behavioural checks of the workflow port: Latin-hypercube stratification,
//! timestamps, template round trips, a real local campaign run through a shell script, control
//! directive round trips and deploy-file permissions.
use cfd_solvers::workflow::*;
use std::path::PathBuf;

fn tmp(name: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("cfdwf_analytic_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

#[test]
fn latin_hypercube_hits_every_stratum_once() {
    let params = vec![
        CampaignParameter::continuous("a", 0.0, 1.0),
        CampaignParameter::continuous("b", -10.0, 10.0),
    ];
    for n in [1usize, 4, 17, 64] {
        let cases = generate_latin_hypercube_campaign(
            &params,
            n,
            CampaignRandom {
                state: 99 + n as u64,
            },
        )
        .unwrap();
        for (name, lo, hi) in [("a", 0.0, 1.0), ("b", -10.0, 10.0)] {
            let mut seen = vec![false; n];
            for c in &cases {
                let v: f64 = c.values[name].parse().unwrap();
                let t = (v - lo) / (hi - lo);
                assert!((0.0..1.0).contains(&t));
                let k = ((t * n as f64) as usize).min(n - 1);
                assert!(!seen[k], "stratum {k} of {name} hit twice (n={n})");
                seen[k] = true;
            }
        }
    }
}

#[test]
fn utc_timestamps_and_provenance_file() {
    assert_eq!(format_utc_timestamp(0), "1970-01-01T00:00:00Z");
    assert_eq!(format_utc_timestamp(951_782_400), "2000-02-29T00:00:00Z");
    assert_eq!(format_utc_timestamp(1_794_659_696), "2026-11-14T12:34:56Z");
    assert_eq!(format_utc_timestamp(-1), "1969-12-31T23:59:59Z");
    let p = make_provenance("fvm", "1.0", "unknown");
    assert_eq!(p.timestamp_utc.len(), 20);
    assert!(p.timestamp_utc.starts_with("20") && p.timestamp_utc.ends_with('Z'));
    let dir = tmp("prov");
    write_provenance(&dir.join("p.json"), &p).unwrap();
    let text = std::fs::read_to_string(dir.join("p.json")).unwrap();
    let v: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(v["solver"], "fvm");
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn case_config_json_round_trip() {
    let mut c = CaseConfig::default();
    c.set("u", UnitValue::new(0.1 + 0.2, "m/s"));
    c.set("nu", UnitValue::new(1.5e-5, "m^2/s"));
    let back = case_config_from_json(&c.to_json()).unwrap();
    assert_eq!(back, c);
    assert!(case_config_from_json("{\"x\":{\"value\":--,\"unit\":\"m\"}}").is_err());
}

#[test]
fn gradient_descent_reaches_quadratic_minimum() {
    let q = |x: &[f64]| (x[0] - 1.5).powi(2) + 3.0 * (x[1] + 0.25).powi(2);
    let r = run_gradient_descent_campaign(
        &q,
        vec![4.0, 2.0],
        &CampaignOptimizationConfig {
            max_iterations: 500,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        (r.parameters[0] - 1.5).abs() < 1e-4 && (r.parameters[1] + 0.25).abs() < 1e-4,
        "{r:?}"
    );
    let gr = finite_difference_gradient(&q, &[0.0, 0.0], &Default::default()).unwrap();
    assert!((gr[0] + 3.0).abs() < 1e-8 && (gr[1] - 1.5).abs() < 1e-8);
}

#[cfg(unix)]
#[test]
fn local_campaign_runs_a_real_solver_script() {
    use std::os::unix::fs::PermissionsExt;
    let root = tmp("local");
    let script = root.join("fake-solver.sh");
    std::fs::write(
        &script,
        "#!/bin/sh\nfor i in 1 2 3; do echo \"iter $i res 0.$i\"; done\necho \"perf.wall_time 1.5 s\"\ncase \"$PWD\" in *case0002) echo FATAL >&2; exit 3;; esac\necho DONE\n",
    )
    .unwrap();
    std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
    let camp = root.join("camp");
    let cases =
        generate_factorial_campaign(&[CampaignParameter::discrete("mesh", ["a", "b"])]).unwrap();
    write_campaign_case_folders(
        &camp,
        &cases,
        &[CampaignTemplateFile {
            relative_path: "case.txt".into(),
            contents: "mesh {mesh}\n".into(),
        }],
        &Default::default(),
    )
    .unwrap();
    let adapter = SolverAdapterDescriptor {
        name: "fake".into(),
        native_binary: script.to_string_lossy().into_owned(),
        container_binary: "fake".into(),
        capabilities: vec![SolverCapability::Residuals, SolverCapability::Control],
        control_actions: vec!["stop".into()],
    };
    let out = run_local_campaign(
        &camp,
        &LocalCampaignRunOptions {
            adapter: adapter.clone(),
            ..Default::default()
        },
        "registry.tsv",
    )
    .unwrap();
    assert_eq!(out.launched, 2);
    assert_eq!((out.summary.done, out.summary.failed), (1, 1));
    assert_eq!(out.summary.best_case_id, "case0001");
    assert!((out.summary.best_objective - 0.3).abs() < 1e-15);
    let r0 = &out.case_results[0];
    assert_eq!(r0.residuals.len(), 3);
    assert_eq!(r0.performance[0].key, "wall_time");
    assert_eq!(out.case_results[1].exit_code, 3);
    let reg = load_campaign_registry(&camp.join("registry.tsv")).unwrap();
    assert_eq!(reg[1].status, CampaignCaseStatus::Failed);
    assert_eq!(reg[0].message, "exit=0 residuals=3 perf=1");
    // A second pass has nothing pending.
    let again = run_local_campaign(
        &camp,
        &LocalCampaignRunOptions {
            adapter,
            ..Default::default()
        },
        "registry.tsv",
    )
    .unwrap();
    assert_eq!(again.launched, 0);
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn missing_launcher_marks_case_failed_and_registry_reloads() {
    let root = tmp("missing");
    let camp = root.join("camp");
    let cases = generate_factorial_campaign(&[CampaignParameter::discrete("m", ["x"])]).unwrap();
    write_campaign_case_folders(
        &camp,
        &cases,
        &[CampaignTemplateFile {
            relative_path: "f".into(),
            contents: "{m}\n".into(),
        }],
        &Default::default(),
    )
    .unwrap();
    let adapter = SolverAdapterDescriptor {
        name: "n".into(),
        native_binary: "/nonexistent/solver-binary".into(),
        container_binary: "n".into(),
        ..Default::default()
    };
    let out = run_local_campaign(
        &camp,
        &LocalCampaignRunOptions {
            adapter,
            ..Default::default()
        },
        "registry.tsv",
    )
    .unwrap();
    assert_eq!((out.launched, out.summary.failed), (0, 1));
    let reg = load_campaign_registry(&camp.join("registry.tsv")).unwrap();
    assert!(reg[0].objective.is_infinite());
    assert!(
        reg[0].message.starts_with("runtime doctor failed"),
        "{}",
        reg[0].message
    );
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn control_directive_round_trip_and_capability_gate() {
    let root = tmp("ctl");
    std::fs::create_dir_all(root.join("case0001")).unwrap();
    let mut adapter = SolverAdapterDescriptor {
        name: "s".into(),
        native_binary: "s".into(),
        container_binary: "s".into(),
        capabilities: vec![],
        control_actions: vec!["extend".into()],
    };
    let req = CampaignControlRequest {
        case_id: "case0001".into(),
        action: CampaignControlAction::Extend,
        value: "a=b\\c\nd".into(),
        reason: "r\r".into(),
    };
    assert!(write_campaign_control_directive(&root, &adapter, &req, &Default::default()).is_err());
    adapter.capabilities.push(SolverCapability::Control);
    write_campaign_control_directive(&root, &adapter, &req, &Default::default()).unwrap();
    let found =
        discover_campaign_control_directives(&root.join("case0001"), &Default::default()).unwrap();
    assert_eq!(found.len(), 2); // extend.directive + latest.directive
    for d in &found {
        assert_eq!((d.value.as_str(), d.reason.as_str()), ("a=b\\c\nd", "r\r"));
    }
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn deploy_scripts_are_executable_and_paths_stay_inside_root() {
    let root = tmp("docker");
    let r = write_docker_deploy_system(&root, &Default::default(), true).unwrap();
    assert_eq!(r.files_written, 9);
    assert!(write_docker_deploy_system(&root, &Default::default(), false).is_err());
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mode = std::fs::metadata(root.join("deploy/deploy.sh"))
            .unwrap()
            .permissions()
            .mode();
        assert!(mode & 0o100 != 0);
    }
    assert_eq!(shell_quote("a b'c"), "'a b'\\''c'");
    assert_eq!(shell_quote(""), "''");
    assert_eq!(shell_quote("/x/y:z=1,2"), "/x/y:z=1,2");
    let _ = std::fs::remove_dir_all(root);
}
