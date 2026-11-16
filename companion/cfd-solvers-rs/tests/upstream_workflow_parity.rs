//! Text-exact comparison with the C++ upstream workflow/{schema, provenance, validation_catalog,
//! campaign, campaign_control, deploy}. Goldens come from `parity/oracle_workflow.cpp`, which was
//! run with output root `/tmp/cfdwf_oracle`.
#![allow(clippy::needless_range_loop, clippy::type_complexity)]
use cfd_solvers::workflow::*;
use serde_json::Value;
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

const ORACLE_ROOT: &str = "/tmp/cfdwf_oracle";

fn g() -> &'static Value {
    static G: std::sync::OnceLock<Value> = std::sync::OnceLock::new();
    G.get_or_init(|| serde_json::from_str(include_str!("../parity/goldens-workflow.json")).unwrap())
}

#[track_caller]
fn s(name: &str, rust: &str) {
    assert_eq!(rust, g()[name].as_str().unwrap(), "{name}");
}

#[track_caller]
fn ss(name: &str, rust: &[String]) {
    let cpp: Vec<String> = g()[name]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect();
    assert_eq!(rust, cpp.as_slice(), "{name}");
}

fn num(v: &Value) -> f64 {
    match v.as_str() {
        Some("inf") => f64::INFINITY,
        Some("-inf") => f64::NEG_INFINITY,
        Some(_) => f64::NAN,
        None => v.as_f64().unwrap(),
    }
}

#[track_caller]
fn v(name: &str, rust: &[f64], rel: f64) {
    let cpp: Vec<f64> = g()[name].as_array().unwrap().iter().map(num).collect();
    assert_eq!(rust.len(), cpp.len(), "{name} length");
    for (i, (a, b)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            a == b || (a - b).abs() <= rel * b.abs().max(1e-300),
            "{name}[{i}]: rust {a} cpp {b}"
        );
    }
}

fn err<T>(r: cfd_solvers::Result<T>) -> String {
    match r {
        Ok(_) => String::new(),
        Err(e) => {
            let t = e.to_string();
            t.split_once(": ").map(|(_, m)| m.to_string()).unwrap_or(t)
        }
    }
}

fn g17(x: f64) -> String {
    cfd_solvers::lbm::format_general(x, 17)
}

fn tmp(name: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("cfdwf_rust_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&p);
    p
}

fn slurp(p: &Path) -> String {
    std::fs::read_to_string(p).unwrap()
}

fn cases_text(cs: &[CampaignCase]) -> Vec<String> {
    cs.iter()
        .map(|c| {
            let mut s = c.case_id.clone();
            for (k, v) in &c.values {
                s += &format!("|{k}={v}");
            }
            s
        })
        .collect()
}

fn registry_text(r: &[CampaignRegistryEntry]) -> Vec<String> {
    r.iter()
        .map(|e| {
            format!(
                "{}|{}|{}|{}|{}",
                e.case_id,
                campaign_status_name(e.status),
                g17(e.objective),
                e.iterations,
                e.message
            )
        })
        .collect()
}

fn summary_text(s: &CampaignSummary) -> String {
    format!(
        "{} {} {} {} {} {} {}",
        s.pending,
        s.running,
        s.done,
        s.failed,
        s.stopped,
        g17(s.best_objective),
        s.best_case_id
    )
}

fn cmds_text(cs: &[&RemoteCommandPlan]) -> Vec<String> {
    cs.iter()
        .map(|c| {
            let mut s = format!(
                "{}|{}|{}|{}",
                c.server_name,
                c.description,
                if c.uses_ssh { "ssh" } else { "local" },
                c.display_command
            );
            for a in &c.argv {
                s += &format!("|{a}");
            }
            s
        })
        .collect()
}

fn rule(
    key: &str,
    unit: &str,
    min: Option<f64>,
    max: Option<f64>,
    required: bool,
) -> ParameterRule {
    ParameterRule {
        key: key.into(),
        unit: unit.into(),
        minimum: min,
        maximum: max,
        required,
    }
}

#[test]
fn schema_provenance_catalog() {
    let mut cfg = CaseConfig::default();
    cfg.set("velocity", UnitValue::new(1.25, "m/s"));
    cfg.set("density", UnitValue::new(998.2, "kg/m^3"));
    cfg.set("a\"b", UnitValue::new(1e-7, "x\\y"));
    cfg.set("temperature", UnitValue::new(300.0 / 7.0, "K"));
    s("case_json", &cfg.to_json());
    s(
        "validate_ok",
        &err(cfg.validate(&[
            rule("velocity", "m/s", Some(0.0), Some(10.0), true),
            rule("missing", "", None, None, false),
        ])),
    );
    s(
        "validate_missing",
        &err(cfg.validate(&[rule("missing", "", None, None, true)])),
    );
    s(
        "validate_unit",
        &err(cfg.validate(&[rule("velocity", "km/h", None, None, true)])),
    );
    s(
        "validate_min",
        &err(cfg.validate(&[rule("velocity", "", Some(2.0), None, true)])),
    );
    s(
        "validate_max",
        &err(cfg.validate(&[rule("velocity", "", None, Some(1.0), true)])),
    );
    let parsed = case_config_from_json("{\"parameters\":{\"p\":{\"value\":101325,\"unit\":\"Pa\"}, \"q\" : { \"value\" : -1.5e-3 , \"unit\" : \"\" }, \"bad\":{\"value\":\"x\",\"unit\":\"m\"}, \"r\":{\"value\":2.5e,\"unit\":\"s\"}}}").unwrap();
    s("from_json", &parsed.to_json());
    let yaml = case_config_from_yaml("# comment\nvelocity: 2.5 m/s\n  pressure : -1e5 Pa  \nbad: 3\nunit_hash: 4 m#x\nT: 300 K\r\nlast:7 s").unwrap();
    s("from_yaml", &yaml.to_json());
    let mut db = MaterialDatabase::default();
    let mat = |name: &str, props: &[(&str, f64, &str)]| Material {
        name: name.into(),
        property: props
            .iter()
            .map(|(k, v, u)| (k.to_string(), UnitValue::new(*v, *u)))
            .collect::<BTreeMap<_, _>>(),
    };
    db.add(mat(
        "water",
        &[("rho", 998.2, "kg/m^3"), ("k", 0.6, "W/m/K")],
    ))
    .unwrap();
    db.add(mat("air", &[("rho", 1.225, "kg/m^3")])).unwrap();
    db.add(mat("empty", &[])).unwrap();
    s("materials", &db.to_json());
    s("material_empty_name", &err(db.add(mat("", &[]))));
    let prov = ResultProvenance {
        solver: "fvm".into(),
        version: "0.4.0".into(),
        git_commit: "abc123".into(),
        timestamp_utc: "2026-11-14T12:34:56Z".into(),
        metadata: [("zeta", "z"), ("alpha", "q\"uote")]
            .iter()
            .map(|(a, b)| (a.to_string(), b.to_string()))
            .collect(),
    };
    s("provenance", &provenance_json(&prov));
    let cat: Vec<String> = validation_catalog()
        .iter()
        .map(|c| format!("{}|{}|{}|{}", c.name, c.family, c.command, c.criterion))
        .collect();
    ss("catalog", &cat);
}

fn mixed() -> Vec<CampaignParameter> {
    vec![
        CampaignParameter::continuous("u", 0.5, 2.5),
        CampaignParameter::continuous("T", 280.0, 320.0),
        CampaignParameter::discrete("mesh", ["a", "b", "c"]),
    ]
}

fn adapter() -> SolverAdapterDescriptor {
    SolverAdapterDescriptor {
        name: "cfd-solve".into(),
        native_binary: "cfd-solve".into(),
        container_binary: "/opt/cfd/bin/cfd-solve".into(),
        capabilities: vec![SolverCapability::Residuals, SolverCapability::Control],
        control_actions: vec!["stop".into(), "checkpoint".into()],
    }
}

fn base_registry() -> Vec<CampaignRegistryEntry> {
    let mut reg = Vec::new();
    for e in [
        CampaignRegistryEntry::new("case0001", CampaignCaseStatus::Done, 0.25, 12, "ok"),
        CampaignRegistryEntry::new(
            "case0002",
            CampaignCaseStatus::Done,
            0.125,
            7,
            "tab\there\nnew\\line",
        ),
        CampaignRegistryEntry::new("case0003", CampaignCaseStatus::Failed, 1.0 / 3.0, 0, ""),
        CampaignRegistryEntry::new("case0001", CampaignCaseStatus::Done, 0.5, 13, "again"),
        CampaignRegistryEntry::new("case0004", CampaignCaseStatus::Running, 0.0, 0, "r"),
    ] {
        update_case_status(&mut reg, e).unwrap();
    }
    reg
}

#[test]
fn campaign_generation_rendering_registry() {
    let mut rng = CampaignRandom::default();
    let u: Vec<f64> = (0..8).map(|_| rng.uniform_open01()).collect();
    v("uniform", &u, 1e-15); // serde_json may round the 17-digit text by one ulp
    let disc = vec![
        CampaignParameter::discrete("mesh", ["coarse", "fine"]),
        CampaignParameter::discrete("model", ["kEps", "kOmega", "SA"]),
    ];
    ss(
        "factorial",
        &cases_text(&generate_factorial_campaign(&disc).unwrap()),
    );
    ss(
        "lhs",
        &cases_text(
            &generate_latin_hypercube_campaign(&mixed(), 7, CampaignRandom::default()).unwrap(),
        ),
    );
    ss(
        "lhs_seed",
        &cases_text(
            &generate_latin_hypercube_campaign(&mixed(), 5, CampaignRandom { state: 12345 })
                .unwrap(),
        ),
    );
    s("factorial_err", &err(generate_factorial_campaign(&mixed())));
    s(
        "reserved_err",
        &err(generate_factorial_campaign(&[CampaignParameter::discrete(
            "case_id",
            ["x"],
        )])),
    );
    let row = CampaignCase {
        case_id: "case0001".into(),
        values: [
            ("mesh", "fine"),
            ("u", "1.5"),
            ("turb", "1"),
            ("off", "false"),
        ]
        .iter()
        .map(|(a, b)| (a.to_string(), b.to_string()))
        .collect(),
    };
    let templ = "u={u} mesh={mesh} {unknown_brace\n<!-- IF turb -->\nturbulence on\n<!-- IF mesh=fine -->\nfine {mesh}{u}\n<!-- ENDIF -->\n<!-- ENDIF -->\n<!-- IF off -->\nhidden\n<!-- ENDIF -->\n<!-- IF mesh != coarse -->\nnot coarse\n<!-- ENDIF -->\n{{u}}\nend";
    s(
        "render",
        &render_campaign_template(templ, &row, false).unwrap(),
    );
    s(
        "render_strict_err",
        &err(render_campaign_template("x {nope}", &row, true)),
    );
    s(
        "render_unclosed",
        &err(render_campaign_template("<!-- IF u -->\nx", &row, true)),
    );
    s(
        "render_endif",
        &err(render_campaign_template("<!-- ENDIF -->", &row, true)),
    );

    let reg = base_registry();
    s("summary", &summary_text(&summarize_campaign_registry(&reg)));
    let root = tmp("reg");
    let path = root.join("reg").join("registry.tsv");
    save_campaign_registry(&path, &reg).unwrap();
    s("registry_tsv", &slurp(&path));
    ss(
        "registry_loaded",
        &registry_text(&load_campaign_registry(&path).unwrap()),
    );
    let names: Vec<String> = [
        "pending",
        "RUNNING",
        "Success",
        "completed",
        "error",
        "cancelled",
        "STOPPED",
    ]
    .iter()
    .map(|n| campaign_status_name(campaign_status_from_name(n).unwrap()).to_string())
    .collect();
    ss("status_names", &names);

    let mut plans = Vec::new();
    for rt in [
        SolverRuntime::Native,
        SolverRuntime::Docker,
        SolverRuntime::Singularity,
        SolverRuntime::Slurm,
    ] {
        let req = SolverRunRequest {
            adapter: adapter(),
            case_directory: "/data/my case".into(),
            runtime: rt,
            mpi_ranks: 4,
            threads: 2,
            extra_args: vec!["--flag".into(), "it's".into()],
            docker_image: "img:1".into(),
            singularity_image: "img.sif".into(),
            slurm_partition: "gpu".into(),
        };
        let p = build_solver_command_plan(&req).unwrap();
        let mut t = format!(
            "{}|{}|{}",
            p.runtime.ordinal(),
            u8::from(p.requires_scheduler),
            p.display_command
        );
        for (k, val) in &p.environment {
            t += &format!("|{k}={val}");
        }
        for a in &p.argv {
            t += &format!("|{a}");
        }
        plans.push(t);
    }
    ss("plans", &plans);
    s(
        "plan_docker_err",
        &err(build_solver_command_plan(&SolverRunRequest {
            adapter: adapter(),
            case_directory: "/x".into(),
            runtime: SolverRuntime::Docker,
            ..Default::default()
        })),
    );

    let log = "iter 1 res 1e-1\nIteration=2, residual=5e-2 field=U\n3 0.025 extra\nEQUATION p ITER 4 RES 1.25e-2\nnoise 5\n7 inf\nperf.wall_time 12.5 s\nWALL_TIME: 3\nCELLS_PER_S=1e6 cells/s\nperf.mem 1.5\nother 1\n";
    let res = parse_residual_history(log);
    let flat: Vec<f64> = res
        .iter()
        .flat_map(|r| [r.iteration as f64, r.residual])
        .collect();
    v("residuals", &flat, 0.0);
    ss(
        "residual_eq",
        &res.iter().map(|r| r.equation.clone()).collect::<Vec<_>>(),
    );
    let perf: Vec<String> = parse_performance_metrics(log)
        .iter()
        .map(|m| format!("{}|{}|{}", m.key, g17(m.value), m.units))
        .collect();
    ss("perf", &perf);
    let outcomes: Vec<String> = [
        "all done",
        "Fatal: x",
        "job cancelled",
        "iter 3",
        "SUCCESS but error",
    ]
    .iter()
    .map(|l| campaign_status_name(detect_solver_outcome_from_log(l)).to_string())
    .collect();
    ss("outcomes", &outcomes);
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn case_folders_control_and_refresh() {
    let root = tmp("camp");
    let camp = root.join("camp");
    let tf = vec![
        CampaignTemplateFile {
            relative_path: "system/controlDict".into(),
            contents: "u {u}\nmesh {mesh}\n".into(),
        },
        CampaignTemplateFile {
            relative_path: "run.sh".into(),
            contents: "#!/bin/sh\necho {mesh}\n<!-- IF mesh=b -->\necho b\n<!-- ENDIF -->\n".into(),
        },
    ];
    let lcases = generate_latin_hypercube_campaign(&mixed(), 3, CampaignRandom::default()).unwrap();
    let pr =
        write_campaign_case_folders(&camp, &lcases, &tf, &CampaignWriteOptions::default()).unwrap();
    v(
        "persist",
        &[
            pr.cases_written as f64,
            pr.files_written as f64,
            pr.registry_entries as f64,
        ],
        0.0,
    );
    s("persist_registry", &slurp(&camp.join("registry.tsv")));
    s(
        "persist_manifest",
        &slurp(&camp.join("template_manifest.tsv")),
    );
    s("persist_doe", &slurp(&camp.join("case0002/doe_row.csv")));
    s(
        "persist_ctrl",
        &slurp(&camp.join("case0002/system/controlDict")),
    );
    s("persist_run", &slurp(&camp.join("case0003/run.sh")));
    let e = err(write_campaign_case_folders(
        &camp,
        &lcases,
        &tf,
        &CampaignWriteOptions::default(),
    ));
    assert!(
        e.starts_with(g()["persist_exists_err"].as_str().unwrap()),
        "{e}"
    );
    s(
        "persist_escape_err",
        &err(write_campaign_case_folders(
            &root.join("camp2"),
            &lcases,
            &[CampaignTemplateFile {
                relative_path: "../x".into(),
                contents: String::new(),
            }],
            &CampaignWriteOptions::default(),
        )),
    );

    let c1 = camp.join("case0001");
    std::fs::create_dir_all(c1.join("logs")).unwrap();
    std::fs::write(
        c1.join("logs/solver.log"),
        "iter 1 res 0.5\niter 2 res 0.25\n",
    )
    .unwrap();
    std::fs::write(c1.join("residual_history.dat"), "3 0.125\n").unwrap();
    std::fs::write(c1.join("timing.txt"), "perf.wall 2 s\n").unwrap();
    std::fs::write(c1.join("notes.md"), "DONE\n").unwrap();
    std::fs::write(camp.join("case0002/run.out"), "iter 9 res 1e-3\nFAILED\n").unwrap();
    let d = discover_campaign_outputs(&c1);
    let rel = |p: &PathBuf| {
        p.strip_prefix(&root)
            .unwrap()
            .to_string_lossy()
            .into_owned()
    };
    let mut found = Vec::new();
    found.extend(d.log_files.iter().map(|p| format!("L{}", rel(p))));
    found.extend(d.residual_files.iter().map(|p| format!("R{}", rel(p))));
    found.extend(d.performance_files.iter().map(|p| format!("P{}", rel(p))));
    ss("discover", &found);
    s(
        "refresh_summary",
        &summary_text(&refresh_campaign_registry_from_outputs(&camp, "registry.tsv").unwrap()),
    );
    ss(
        "refresh_registry",
        &registry_text(&load_campaign_registry(&camp.join("registry.tsv")).unwrap()),
    );

    let creq = CampaignControlRequest {
        case_id: "case0002".into(),
        action: CampaignControlAction::Checkpoint,
        value: "100\nsteps".into(),
        reason: "user\\req".into(),
    };
    let dir = write_campaign_control_directive(
        &camp,
        &adapter(),
        &creq,
        &CampaignControlOptions::default(),
    )
    .unwrap();
    s("directive", &slurp(&dir.path));
    s("directive_name", &rel(&dir.path));
    s(
        "control_err",
        &err(write_campaign_control_directive(
            &camp,
            &adapter(),
            &CampaignControlRequest {
                case_id: "case0002".into(),
                action: CampaignControlAction::Extend,
                ..Default::default()
            },
            &CampaignControlOptions::default(),
        )),
    );
    let ds: Vec<String> = discover_campaign_control_directives(
        &camp.join("case0002"),
        &CampaignControlOptions::default(),
    )
    .unwrap()
    .iter()
    .map(|d| {
        format!(
            "{}|{}|{}|{}|{}",
            d.case_id,
            campaign_control_action_name(d.action),
            d.value,
            d.reason,
            d.path.file_name().unwrap().to_string_lossy()
        )
    })
    .collect();
    ss("directives", &ds);
    let jobs = parse_slurm_queue_table(
        "JOBID NAME ST TIME\n# c\n101,case0001,R,12.5\n102 case0005 PD\n103\tcase0004\tCD\t7\n104 case0003 WEIRD x\nshort row\n",
    );
    let jt: Vec<String> = jobs
        .iter()
        .map(|j| {
            format!(
                "{}|{}|{}|{}|{}",
                j.job_id,
                j.case_id,
                scheduler_job_state_name(j.state),
                g17(j.elapsed_seconds),
                j.raw_state
            )
        })
        .collect();
    ss("slurm", &jt);
    let mut r2 = base_registry();
    s(
        "sched_summary",
        &summary_text(&apply_scheduler_records_to_registry(&mut r2, &jobs)),
    );
    ss("sched_registry", &registry_text(&r2));
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn finite_difference_and_descent() {
    let rosen = |x: &[f64]| {
        (1.0 - x[0]) * (1.0 - x[0]) + 10.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0])
    };
    v(
        "fd_central",
        &finite_difference_gradient(&rosen, &[0.3, -0.2], &Default::default()).unwrap(),
        1e-12,
    );
    v(
        "fd_forward",
        &finite_difference_gradient(
            &rosen,
            &[30.0, -0.2],
            &FiniteDifferenceGradientConfig {
                absolute_step: 1e-6,
                relative_step: 1e-4,
                central: false,
            },
        )
        .unwrap(),
        1e-12,
    );
    let r = run_gradient_descent_campaign(
        &rosen,
        vec![-0.5, 0.5],
        &CampaignOptimizationConfig {
            max_iterations: 200,
            step_size: 0.05,
            tolerance: 1e-8,
            finite_difference: Default::default(),
        },
    )
    .unwrap();
    let mut out = r.parameters.clone();
    out.extend([
        r.objective,
        r.iterations as f64,
        f64::from(u8::from(r.converged)),
    ]);
    out.extend(&r.last_gradient);
    v("gd", &out, 1e-9);
    let q = |x: &[f64]| (x[0] - 3.0) * (x[0] - 3.0) + 2.0 * (x[1] + 1.0) * (x[1] + 1.0);
    let r = run_gradient_descent_campaign(&q, vec![0.0, 0.0], &Default::default()).unwrap();
    let mut out = r.parameters.clone();
    out.extend([
        r.objective,
        r.iterations as f64,
        f64::from(u8::from(r.converged)),
    ]);
    v("gd_quad", &out, 1e-9);
}

fn servers() -> Vec<CampaignServerDescriptor> {
    let sv = |name: &str,
              host: &str,
              user: &str,
              root: &str,
              rt,
              slots,
              tags: &[&str],
              img: &str,
              online| CampaignServerDescriptor {
        name: name.into(),
        host: host.into(),
        user: user.into(),
        remote_root: root.into(),
        runtime: rt,
        slots,
        tags: tags.iter().map(|t| t.to_string()).collect(),
        docker_image: img.into(),
        online,
    };
    vec![
        sv(
            "gpu-a",
            "10.0.0.11",
            "cfd",
            "/srv/cfd/campaigns",
            CampaignServerRuntime::Docker,
            2,
            &["gpu", "linux"],
            "cfd-solvers:latest",
            true,
        ),
        sv(
            "cpu-b",
            "localhost",
            "",
            "/scratch/my runs",
            CampaignServerRuntime::Native,
            1,
            &["linux"],
            "",
            true,
        ),
        sv(
            "cpu-a",
            "cluster.example",
            "",
            "/srv/cfd",
            CampaignServerRuntime::Native,
            3,
            &["linux"],
            "",
            true,
        ),
        sv(
            "off",
            "10.0.0.99",
            "",
            "/srv",
            CampaignServerRuntime::Native,
            4,
            &["linux"],
            "",
            false,
        ),
    ]
}

#[test]
fn multiserver_deploy_controller_docker() {
    let servers = servers();
    let ad = adapter();
    let preg: Vec<CampaignRegistryEntry> = (1..=8)
        .map(|i| {
            CampaignRegistryEntry::new(
                format!("case{i:04}"),
                if i == 3 {
                    CampaignCaseStatus::Done
                } else {
                    CampaignCaseStatus::Pending
                },
                0.0,
                0,
                "",
            )
        })
        .collect();
    let so = MultiServerSchedulingOptions {
        max_cases_per_server: 2,
        ..Default::default()
    };
    let mplan =
        plan_multi_server_campaign(Path::new("/campaigns/local"), &preg, &ad, &servers, &so)
            .unwrap();
    let at: Vec<String> = mplan
        .assignments
        .iter()
        .map(|a| {
            format!(
                "{}|{}|{}|{}|{}",
                a.case_id,
                a.server_name,
                a.local_case_dir.display(),
                a.remote_case_dir.display(),
                a.solver_plan.display_command
            )
        })
        .collect();
    ss("ms_assign", &at);
    ss(
        "ms_cmds",
        &cmds_text(&mplan.commands.iter().collect::<Vec<_>>()),
    );
    v(
        "ms_counts",
        &[mplan.unassigned_cases as f64, mplan.active_servers as f64],
        0.0,
    );
    let t = MultiServerSchedulingOptions {
        required_tag: "gpu".into(),
        ..Default::default()
    };
    let p = plan_multi_server_campaign(Path::new("/c"), &preg, &ad, &servers, &t).unwrap();
    v(
        "ms_tag",
        &[
            p.assignments.len() as f64,
            p.unassigned_cases as f64,
            p.active_servers as f64,
        ],
        0.0,
    );
    let t = MultiServerSchedulingOptions {
        required_tag: "nope".into(),
        ..Default::default()
    };
    s(
        "ms_err",
        &err(plan_multi_server_campaign(
            Path::new("/c"),
            &preg,
            &ad,
            &servers,
            &t,
        )),
    );

    let root = tmp("ms");
    let ms = root.join("ms");
    let pw = write_multiserver_plan_files(&ms, &mplan, &Default::default()).unwrap();
    s("ms_assign_tsv", &slurp(&pw.assignments_path));
    s("ms_cmd_sh", &slurp(&pw.commands_path));
    let ex = plan_multiserver_execution(&mplan, &servers, ".cfd_run").unwrap();
    let ht: Vec<String> = ex
        .health_checks
        .iter()
        .map(|h| {
            format!(
                "{}|{}|{}|{}",
                h.server_name,
                h.check_name,
                u8::from(h.required),
                h.command.display_command
            )
        })
        .collect();
    ss("ex_health", &ht);
    let jc: Vec<&RemoteCommandPlan> = ex
        .jobs
        .iter()
        .flat_map(|j| {
            [
                &j.launch_command,
                &j.status_command,
                &j.cancel_command,
                &j.fetch_logs_command,
            ]
        })
        .collect();
    ss("ex_jobs", &cmds_text(&jc));
    v(
        "ex_counts",
        &[
            ex.native_jobs as f64,
            ex.docker_jobs as f64,
            ex.remote_jobs as f64,
        ],
        0.0,
    );
    let ew = write_multiserver_execution_files(&ms, &ex, &Default::default()).unwrap();
    s("ex_health_sh", &slurp(&ew.health_script_path));
    s("ex_launch_sh", &slurp(&ew.launch_script_path));
    s("ex_status_sh", &slurp(&ew.status_script_path));
    s("ex_fetch_sh", &slurp(&ew.fetch_logs_script_path));
    s("ex_jobs_tsv", &slurp(&ew.jobs_path));
    let rs = parse_remote_job_status_table(
        "case_id\tserver\tstate\texit\ncase0001\tgpu-a\tdone\t0\ncase0002\tgpu-a\tfailed\t2x\ncase0004\tcpu-b\tmissing\ncase0009\tcpu-a\tqueued\tzz\n",
    )
    .unwrap();
    let rt: Vec<String> = rs
        .iter()
        .map(|r| {
            format!(
                "{}|{}|{}|{}|{}",
                r.case_id,
                r.server_name,
                remote_job_state_name(r.state),
                r.exit_code,
                r.raw_state
            )
        })
        .collect();
    ss("remote_status", &rt);
    let mut r2 = preg.clone();
    s(
        "remote_summary",
        &summary_text(&apply_remote_job_status_to_registry(&mut r2, &rs).unwrap()),
    );
    ss("remote_registry", &registry_text(&r2));
    s(
        "remote_err",
        &err(parse_remote_job_status_table("a\tb\tweird\n")),
    );
    let mk: Vec<String> = [
        "TOKEN=",
        "PASSWORD=",
        "SECRET=",
        "KEY=",
        "AWS_",
        "GITHUB_",
        "",
    ]
    .map(String::from)
    .to_vec();
    ss(
        "redact",
        &[
            redact_sensitive_command_display(
                "env TOKEN=abc123 PASSWORD='p w' AWS_SECRET=xyz;run GITHUB_TOKEN=t KEY=",
                &mk,
            ),
            redact_sensitive_command_display("plain command", &mk),
        ],
    );

    let mut ex2 = ex.clone();
    ex2.jobs[0].launch_command.display_command += " TOKEN=supersecret";
    let sp = plan_multiserver_supervision(&ex2, &servers, &Default::default()).unwrap();
    let pt: Vec<String> = sp
        .access_probes
        .iter()
        .map(|p| {
            format!(
                "{}|{}|{}",
                p.server_name, p.probe_name, p.command.display_command
            )
        })
        .collect();
    ss("sup_probes", &pt);
    let mut sc: Vec<&RemoteCommandPlan> = sp.retry_launches.iter().map(|r| &r.command).collect();
    sc.extend(
        sp.log_tails
            .iter()
            .flat_map(|t| [&t.stdout_tail_command, &t.stderr_tail_command]),
    );
    ss("sup_cmds", &cmds_text(&sc));
    v(
        "sup_counts",
        &[
            sp.sensitive_commands_redacted as f64,
            sp.dashboard_cases.len() as f64,
        ],
        0.0,
    );
    s(
        "dashboard_json",
        &build_multiserver_dashboard_json(&sp.dashboard_cases),
    );
    let sw = write_multiserver_supervision_files(&ms, &sp, &Default::default()).unwrap();
    s("sup_access_sh", &slurp(&sw.access_script_path));
    s("sup_retry_sh", &slurp(&sw.retry_launch_script_path));
    s("sup_tail_sh", &slurp(&sw.tail_logs_script_path));
    s("sup_tsv", &slurp(&sw.summary_path));
    let mut cases = sp.dashboard_cases.clone();
    cases[1].state = "Running".into();
    cases[2].state = "stopped".into();
    cases[3].state = "weird\u{1}".into();
    cases[0].exit_code = -3;
    s(
        "ctl_status",
        &build_multiserver_controller_status_json(&cases),
    );
    s(
        "ctl_events",
        &build_multiserver_controller_events_ndjson(&cases),
    );
    let cc = RemoteControllerConfig {
        dashboard_title: "My \"Campaign\"".into(),
        control_actions: vec!["stop".into(), "flush".into()],
        port: 9443,
        ..Default::default()
    };
    let cp = plan_multiserver_controller(&sp, &cc).unwrap();
    s("ctl_openapi", &cp.openapi_json);
    s("ctl_script", &cp.controller_script);
    s("ctl_html", &cp.index_html);
    s("ctl_js", &cp.dashboard_js);
    s("ctl_css", &cp.dashboard_css);
    let c2 = RemoteControllerConfig {
        enable_event_history: false,
        allow_mutating_actions: false,
        require_authentication: false,
        enable_tls: false,
        ..Default::default()
    };
    let p2 = plan_multiserver_controller(&sp, &c2).unwrap();
    s("ctl_openapi2", &p2.openapi_json);
    s("ctl_script2", &p2.controller_script);
    let cw = write_multiserver_controller_files(&ms, &cp, &Default::default()).unwrap();
    for (k, p) in [
        ("routes", &cw.routes_path),
        ("env", &cw.env_path),
        ("readme", &cw.readme_path),
        ("tokens", &cw.tokens_example_path),
        ("systemd", &cw.systemd_path),
        ("caddy", &cw.reverse_proxy_path),
        ("start", &cw.start_script_path),
    ] {
        let text = slurp(p).replace(
            &ms.to_string_lossy().into_owned(),
            &format!("{ORACLE_ROOT}/ms"),
        );
        s(&format!("ctlw_{k}"), &text);
    }
    s(
        "ctl_err",
        &err(plan_multiserver_controller(
            &sp,
            &RemoteControllerConfig {
                port: 0,
                ..Default::default()
            },
        )),
    );

    let dc = DockerDeployConfig {
        worker_replicas: 3,
        worker_cpus: 4,
        worker_memory: "8g".into(),
        enable_gpu_runtime: true,
        environment: vec![("OMP_NUM_THREADS".into(), "4".into())],
        ..Default::default()
    };
    let flat = |files: Vec<DockerDeployFile>| -> Vec<String> {
        files
            .iter()
            .map(|f| format!("{}\u{1f}{}", f.relative_path.display(), f.contents))
            .collect()
    };
    ss(
        "docker_full",
        &flat(generate_docker_deploy_files(&dc).unwrap()),
    );
    ss(
        "docker_min",
        &flat(
            generate_docker_deploy_files(&DockerDeployConfig {
                include_manager_service: false,
                ..Default::default()
            })
            .unwrap(),
        ),
    );
    let r = write_docker_deploy_system(&root.join("docker"), &dc, true).unwrap();
    v("docker_written", &[r.files_written as f64], 0.0);
    let _ = std::fs::remove_dir_all(&root);
}
