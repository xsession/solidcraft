//! Multi-server campaign deployment planning: server descriptors, case placement, SSH/rsync
//! command plans, supervised remote launch/status/cancel scripts, retry and log-tail plans,
//! dashboard JSON, a generated dependency-free controller (Python HTTP server, dashboard
//! HTML/JS/CSS, OpenAPI, systemd and Caddy examples) and a Docker/Compose scaffold.
//!
//! Everything here *generates text and files*; nothing is executed or contacted over the
//! network. Ported from upstream `cfd/workflow/deploy` (MIT, cfd_solvers contributors).
use super::campaign::{
    build_solver_command_plan, summarize_campaign_registry, update_case_status,
    validate_solver_adapter_descriptor, CampaignCaseStatus, CampaignRegistryEntry, CampaignSummary,
    SolverAdapterDescriptor, SolverCommandPlan, SolverRunRequest, SolverRuntime,
};
use super::text::{
    add_owner_exec, c_isspace, create_dirs, getline_split, lower, path_text, stoi_prefix,
    write_file,
};
use crate::{require, Result, SolverError};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum CampaignServerRuntime {
    #[default]
    Native,
    Docker,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CampaignServerDescriptor {
    pub name: String,
    pub host: String,
    pub user: String,
    pub remote_root: PathBuf,
    pub runtime: CampaignServerRuntime,
    pub slots: usize,
    pub tags: Vec<String>,
    pub docker_image: String,
    pub online: bool,
}

impl Default for CampaignServerDescriptor {
    fn default() -> Self {
        Self {
            name: String::new(),
            host: "localhost".into(),
            user: String::new(),
            remote_root: PathBuf::new(),
            runtime: CampaignServerRuntime::Native,
            slots: 1,
            tags: Vec::new(),
            docker_image: String::new(),
            online: true,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiServerSchedulingOptions {
    pub pending_only: bool,
    pub require_online: bool,
    pub required_tag: String,
    /// 0 = no extra cap beyond the server slots.
    pub max_cases_per_server: usize,
}

impl Default for MultiServerSchedulingOptions {
    fn default() -> Self {
        Self {
            pending_only: true,
            require_online: true,
            required_tag: String::new(),
            max_cases_per_server: 0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteCommandPlan {
    pub server_name: String,
    pub description: String,
    pub argv: Vec<String>,
    pub display_command: String,
    pub uses_ssh: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ServerCaseAssignment {
    pub case_id: String,
    pub server_name: String,
    pub local_case_dir: PathBuf,
    pub remote_case_dir: PathBuf,
    pub solver_plan: SolverCommandPlan,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MultiServerCampaignPlan {
    pub assignments: Vec<ServerCaseAssignment>,
    pub commands: Vec<RemoteCommandPlan>,
    pub unassigned_cases: usize,
    pub active_servers: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiServerPlanWriteOptions {
    pub overwrite: bool,
    pub assignments_filename: String,
    pub commands_filename: String,
}

impl Default for MultiServerPlanWriteOptions {
    fn default() -> Self {
        Self {
            overwrite: true,
            assignments_filename: "multiserver_assignments.tsv".into(),
            commands_filename: "multiserver_commands.sh".into(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MultiServerPlanWriteResult {
    pub assignments_path: PathBuf,
    pub commands_path: PathBuf,
    pub assignments_written: usize,
    pub commands_written: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum RemoteJobState {
    #[default]
    Planned,
    Submitted,
    Running,
    Done,
    Failed,
    Cancelled,
    Lost,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteHostHealthCheck {
    pub server_name: String,
    pub check_name: String,
    pub command: RemoteCommandPlan,
    pub required: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteSupervisedJobPlan {
    pub case_id: String,
    pub server_name: String,
    pub remote_case_dir: PathBuf,
    pub run_dir: PathBuf,
    pub stdout_path: PathBuf,
    pub stderr_path: PathBuf,
    pub pid_path: PathBuf,
    pub exit_code_path: PathBuf,
    pub launch_command: RemoteCommandPlan,
    pub status_command: RemoteCommandPlan,
    pub cancel_command: RemoteCommandPlan,
    pub fetch_logs_command: RemoteCommandPlan,
    pub solver_plan: SolverCommandPlan,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MultiServerExecutionPlan {
    pub health_checks: Vec<RemoteHostHealthCheck>,
    pub jobs: Vec<RemoteSupervisedJobPlan>,
    pub native_jobs: usize,
    pub docker_jobs: usize,
    pub remote_jobs: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiServerExecutionWriteOptions {
    pub overwrite: bool,
    pub health_filename: String,
    pub launch_filename: String,
    pub status_filename: String,
    pub fetch_logs_filename: String,
    pub jobs_filename: String,
}

impl Default for MultiServerExecutionWriteOptions {
    fn default() -> Self {
        Self {
            overwrite: true,
            health_filename: "multiserver_health.sh".into(),
            launch_filename: "multiserver_launch.sh".into(),
            status_filename: "multiserver_status.sh".into(),
            fetch_logs_filename: "multiserver_fetch_logs.sh".into(),
            jobs_filename: "multiserver_jobs.tsv".into(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MultiServerExecutionWriteResult {
    pub health_script_path: PathBuf,
    pub launch_script_path: PathBuf,
    pub status_script_path: PathBuf,
    pub fetch_logs_script_path: PathBuf,
    pub jobs_path: PathBuf,
    pub health_checks_written: usize,
    pub jobs_written: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteJobStatusRecord {
    pub case_id: String,
    pub server_name: String,
    pub state: RemoteJobState,
    pub exit_code: i32,
    pub raw_state: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemoteSupervisionOptions {
    pub max_attempts: usize,
    pub retry_delay_seconds: usize,
    pub log_tail_lines: usize,
    pub include_ssh_probe: bool,
    pub include_docker_probe: bool,
    pub redact_sensitive_values: bool,
    pub sensitive_markers: Vec<String>,
}

impl Default for RemoteSupervisionOptions {
    fn default() -> Self {
        Self {
            max_attempts: 3,
            retry_delay_seconds: 5,
            log_tail_lines: 120,
            include_ssh_probe: true,
            include_docker_probe: true,
            redact_sensitive_values: true,
            sensitive_markers: ["TOKEN=", "PASSWORD=", "SECRET=", "KEY=", "AWS_", "GITHUB_"]
                .map(String::from)
                .to_vec(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteAccessProbe {
    pub server_name: String,
    pub probe_name: String,
    pub command: RemoteCommandPlan,
    pub required: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteRetryCommandPlan {
    pub case_id: String,
    pub server_name: String,
    pub max_attempts: usize,
    pub retry_delay_seconds: usize,
    pub command: RemoteCommandPlan,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteLogTailPlan {
    pub case_id: String,
    pub server_name: String,
    pub stdout_tail_command: RemoteCommandPlan,
    pub stderr_tail_command: RemoteCommandPlan,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteDashboardCase {
    pub case_id: String,
    pub server_name: String,
    pub state: String,
    pub exit_code: i32,
    pub stdout_path: PathBuf,
    pub stderr_path: PathBuf,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MultiServerSupervisionPlan {
    pub access_probes: Vec<RemoteAccessProbe>,
    pub retry_launches: Vec<RemoteRetryCommandPlan>,
    pub log_tails: Vec<RemoteLogTailPlan>,
    pub dashboard_cases: Vec<RemoteDashboardCase>,
    pub sensitive_commands_redacted: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiServerSupervisionWriteOptions {
    pub overwrite: bool,
    pub access_filename: String,
    pub retry_launch_filename: String,
    pub tail_logs_filename: String,
    pub dashboard_filename: String,
    pub summary_filename: String,
}

impl Default for MultiServerSupervisionWriteOptions {
    fn default() -> Self {
        Self {
            overwrite: true,
            access_filename: "multiserver_access_checks.sh".into(),
            retry_launch_filename: "multiserver_retry_launch.sh".into(),
            tail_logs_filename: "multiserver_tail_logs.sh".into(),
            dashboard_filename: "multiserver_dashboard.json".into(),
            summary_filename: "multiserver_supervision.tsv".into(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MultiServerSupervisionWriteResult {
    pub access_script_path: PathBuf,
    pub retry_launch_script_path: PathBuf,
    pub tail_logs_script_path: PathBuf,
    pub dashboard_json_path: PathBuf,
    pub summary_path: PathBuf,
    pub access_checks_written: usize,
    pub retry_launches_written: usize,
    pub log_tail_commands_written: usize,
    pub dashboard_cases_written: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemoteControllerConfig {
    pub bind_host: String,
    pub port: usize,
    pub token_env_var: String,
    pub campaign_root_env_var: String,
    pub log_tail_lines: usize,
    pub dashboard_refresh_seconds: usize,
    pub dashboard_title: String,
    pub expose_event_snapshot: bool,
    pub expose_event_stream: bool,
    pub event_stream_heartbeat_seconds: usize,
    pub require_authentication: bool,
    pub token_file_env_var: String,
    pub enable_tls: bool,
    pub tls_cert_env_var: String,
    pub tls_key_env_var: String,
    pub enable_event_history: bool,
    pub max_event_history_bytes: usize,
    pub enable_live_probes: bool,
    pub probe_timeout_seconds: usize,
    pub allow_mutating_actions: bool,
    pub require_token_for_mutation: bool,
    pub control_actions: Vec<String>,
}

impl Default for RemoteControllerConfig {
    fn default() -> Self {
        Self {
            bind_host: "127.0.0.1".into(),
            port: 8080,
            token_env_var: "CFD_CONTROLLER_TOKEN".into(),
            campaign_root_env_var: "CFD_CAMPAIGN_ROOT".into(),
            log_tail_lines: 200,
            dashboard_refresh_seconds: 5,
            dashboard_title: "cfd_solvers campaign dashboard".into(),
            expose_event_snapshot: true,
            expose_event_stream: true,
            event_stream_heartbeat_seconds: 15,
            require_authentication: true,
            token_file_env_var: "CFD_CONTROLLER_TOKEN_FILE".into(),
            enable_tls: true,
            tls_cert_env_var: "CFD_CONTROLLER_TLS_CERT".into(),
            tls_key_env_var: "CFD_CONTROLLER_TLS_KEY".into(),
            enable_event_history: true,
            max_event_history_bytes: 10_485_760,
            enable_live_probes: true,
            probe_timeout_seconds: 15,
            allow_mutating_actions: true,
            require_token_for_mutation: true,
            control_actions: ["stop", "extend", "checkpoint", "flush"]
                .map(String::from)
                .to_vec(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoteControllerRoute {
    pub method: String,
    pub path: String,
    pub purpose: String,
    pub requires_token: bool,
    pub mutates: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MultiServerControllerPlan {
    pub config: RemoteControllerConfig,
    pub routes: Vec<RemoteControllerRoute>,
    pub cases: Vec<RemoteDashboardCase>,
    pub openapi_json: String,
    pub status_json: String,
    pub controller_script: String,
    pub start_script: String,
    pub index_html: String,
    pub dashboard_js: String,
    pub dashboard_css: String,
    pub events_ndjson: String,
    pub tokens_example_json: String,
    pub systemd_unit: String,
    pub reverse_proxy_config: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiServerControllerWriteOptions {
    pub overwrite: bool,
    pub controller_dirname: String,
    pub controller_script_filename: String,
    pub start_script_filename: String,
    pub openapi_filename: String,
    pub status_filename: String,
    pub routes_filename: String,
    pub env_filename: String,
    pub readme_filename: String,
    pub index_filename: String,
    pub dashboard_js_filename: String,
    pub dashboard_css_filename: String,
    pub events_filename: String,
    pub tokens_example_filename: String,
    pub systemd_filename: String,
    pub reverse_proxy_filename: String,
}

impl Default for MultiServerControllerWriteOptions {
    fn default() -> Self {
        Self {
            overwrite: true,
            controller_dirname: "controller".into(),
            controller_script_filename: "cfd_controller.py".into(),
            start_script_filename: "multiserver_controller.sh".into(),
            openapi_filename: "openapi.json".into(),
            status_filename: "status.json".into(),
            routes_filename: "routes.tsv".into(),
            env_filename: "env.example".into(),
            readme_filename: "README.md".into(),
            index_filename: "index.html".into(),
            dashboard_js_filename: "dashboard.js".into(),
            dashboard_css_filename: "dashboard.css".into(),
            events_filename: "events.ndjson".into(),
            tokens_example_filename: "controller.tokens.example.json".into(),
            systemd_filename: "cfd-controller.service".into(),
            reverse_proxy_filename: "Caddyfile.example".into(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MultiServerControllerWriteResult {
    pub controller_script_path: PathBuf,
    pub start_script_path: PathBuf,
    pub openapi_path: PathBuf,
    pub status_path: PathBuf,
    pub routes_path: PathBuf,
    pub env_path: PathBuf,
    pub readme_path: PathBuf,
    pub index_path: PathBuf,
    pub dashboard_js_path: PathBuf,
    pub dashboard_css_path: PathBuf,
    pub events_path: PathBuf,
    pub tokens_example_path: PathBuf,
    pub systemd_path: PathBuf,
    pub reverse_proxy_path: PathBuf,
    pub routes_written: usize,
    pub cases_written: usize,
    pub static_assets_written: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DockerDeployConfig {
    pub image_name: String,
    pub service_name: String,
    pub compose_project_name: String,
    pub container_campaign_dir: String,
    pub host_campaign_dir: String,
    pub worker_command: String,
    pub worker_replicas: usize,
    pub worker_cpus: usize,
    pub worker_memory: String,
    pub enable_gpu_runtime: bool,
    pub include_manager_service: bool,
    pub manager_command: String,
    pub dashboard_port: usize,
    pub environment: Vec<(String, String)>,
}

impl Default for DockerDeployConfig {
    fn default() -> Self {
        Self {
            image_name: "cfd-solvers:latest".into(),
            service_name: "cfd-campaign-worker".into(),
            compose_project_name: "cfd_solvers".into(),
            container_campaign_dir: "/work/campaigns".into(),
            host_campaign_dir: "./campaigns".into(),
            worker_command: "cfd-solve particle-campaign-local-runner".into(),
            worker_replicas: 1,
            worker_cpus: 0,
            worker_memory: String::new(),
            enable_gpu_runtime: false,
            include_manager_service: true,
            manager_command: "cfd-solve particle-campaign-control".into(),
            dashboard_port: 8080,
            environment: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DockerDeployFile {
    pub relative_path: PathBuf,
    pub contents: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DockerDeployResult {
    pub root: PathBuf,
    pub files: Vec<DockerDeployFile>,
    pub files_written: usize,
}

/// POSIX shell quoting: arguments made only of `[A-Za-z0-9_\-./:=,]` stay bare.
pub fn shell_quote(text: &str) -> String {
    if text.is_empty() {
        return "''".into();
    }
    let safe = text.bytes().all(|c| {
        c.is_ascii_alphanumeric() || matches!(c, b'_' | b'-' | b'.' | b'/' | b':' | b'=' | b',')
    });
    if safe {
        return text.to_string();
    }
    let mut out = String::from("'");
    for c in text.chars() {
        if c == '\'' {
            out.push_str("'\\''");
        } else {
            out.push(c);
        }
    }
    out.push('\'');
    out
}

fn join_display(argv: &[String]) -> String {
    argv.iter()
        .map(|a| shell_quote(a))
        .collect::<Vec<_>>()
        .join(" ")
}

fn has_tag(server: &CampaignServerDescriptor, tag: &str) -> bool {
    tag.is_empty() || server.tags.iter().any(|t| t == tag)
}

fn ssh_target(server: &CampaignServerDescriptor) -> String {
    if server.user.is_empty() {
        server.host.clone()
    } else {
        format!("{}@{}", server.user, server.host)
    }
}

fn is_localhost(host: &str) -> bool {
    matches!(host, "" | "localhost" | "127.0.0.1" | "::1")
}

fn make_command(
    server: &CampaignServerDescriptor,
    description: String,
    remote_argv: Vec<String>,
) -> RemoteCommandPlan {
    let (argv, uses_ssh) = if is_localhost(&server.host) {
        (remote_argv, false)
    } else {
        (
            vec!["ssh".into(), ssh_target(server), join_display(&remote_argv)],
            true,
        )
    };
    RemoteCommandPlan {
        server_name: server.name.clone(),
        description,
        display_command: join_display(&argv),
        argv,
        uses_ssh,
    }
}

fn os(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn make_rsync_command(
    server: &CampaignServerDescriptor,
    local_case: &Path,
    remote_case: &Path,
) -> RemoteCommandPlan {
    let name = local_case
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let source = format!("{}/", os(local_case));
    let local = is_localhost(&server.host);
    let target = if local {
        format!("{}/", os(remote_case))
    } else {
        format!("{}:{}/", ssh_target(server), os(remote_case))
    };
    let argv = vec!["rsync".into(), "-az".into(), source, target];
    RemoteCommandPlan {
        server_name: server.name.clone(),
        description: format!("sync {name} to {}", server.name),
        display_command: join_display(&argv),
        argv,
        uses_ssh: !local,
    }
}

fn ensure_relative(path: &Path) -> Result<()> {
    require(
        !path.as_os_str().is_empty() && !path.is_absolute() && !path.has_root(),
        "deploy file path must be relative",
    )?;
    require(
        !path
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir)),
        "deploy file path must not escape root",
    )
}

fn docker_env_block(config: &DockerDeployConfig, indent: &str) -> String {
    let mut out = format!(
        "{indent}CFD_CAMPAIGN_ROOT: {}\n{indent}CFD_DEPLOY_PROJECT: {}\n",
        config.container_campaign_dir, config.compose_project_name
    );
    for (k, v) in &config.environment {
        out.push_str(&format!("{indent}{k}: {v}\n"));
    }
    out
}

pub fn campaign_server_runtime_name(runtime: CampaignServerRuntime) -> &'static str {
    match runtime {
        CampaignServerRuntime::Native => "native",
        CampaignServerRuntime::Docker => "docker",
    }
}

pub fn campaign_server_runtime_from_name(name: &str) -> Result<CampaignServerRuntime> {
    match lower(name).as_str() {
        "native" => Ok(CampaignServerRuntime::Native),
        "docker" => Ok(CampaignServerRuntime::Docker),
        _ => Err(SolverError::InvalidInput(
            "unknown campaign server runtime".into(),
        )),
    }
}

pub fn validate_campaign_server_descriptor(server: &CampaignServerDescriptor) -> Result<()> {
    require(!server.name.is_empty(), "campaign server requires a name")?;
    require(!server.host.is_empty(), "campaign server requires a host")?;
    require(
        !server.remote_root.as_os_str().is_empty(),
        "campaign server requires remote root",
    )?;
    require(
        server.slots > 0,
        "campaign server requires at least one slot",
    )?;
    require(
        server.runtime != CampaignServerRuntime::Docker || !server.docker_image.is_empty(),
        "docker campaign server requires image",
    )
}

/// Greedy least-loaded placement (ties broken by server name) of registry cases onto active
/// servers within their slot limits; three commands per case: mkdir, rsync, run.
pub fn plan_multi_server_campaign(
    local_root: &Path,
    registry: &[CampaignRegistryEntry],
    adapter: &SolverAdapterDescriptor,
    servers: &[CampaignServerDescriptor],
    options: &MultiServerSchedulingOptions,
) -> Result<MultiServerCampaignPlan> {
    require(
        !local_root.as_os_str().is_empty(),
        "multi-server campaign planning requires local root",
    )?;
    validate_solver_adapter_descriptor(adapter)?;
    require(
        !servers.is_empty(),
        "multi-server campaign planning requires servers",
    )?;
    let mut active: Vec<&CampaignServerDescriptor> = Vec::new();
    for server in servers {
        validate_campaign_server_descriptor(server)?;
        if options.require_online && !server.online {
            continue;
        }
        if !has_tag(server, &options.required_tag) {
            continue;
        }
        active.push(server);
    }
    require(
        !active.is_empty(),
        "no active campaign servers match scheduling options",
    )?;
    let limits: Vec<usize> = active
        .iter()
        .map(|s| {
            if options.max_cases_per_server > 0 {
                s.slots.min(options.max_cases_per_server)
            } else {
                s.slots
            }
        })
        .collect();
    let mut assigned = vec![0usize; active.len()];
    let mut plan = MultiServerCampaignPlan {
        active_servers: active.len(),
        ..Default::default()
    };
    for entry in registry {
        if options.pending_only && entry.status != CampaignCaseStatus::Pending {
            continue;
        }
        let mut best: Option<usize> = None;
        for i in 0..active.len() {
            if assigned[i] >= limits[i] {
                continue;
            }
            let better = match best {
                None => true,
                Some(b) => (assigned[i], &active[i].name) < (assigned[b], &active[b].name),
            };
            if better {
                best = Some(i);
            }
        }
        let Some(b) = best else {
            plan.unassigned_cases += 1;
            continue;
        };
        let server = active[b];
        assigned[b] += 1;
        let local_case = local_root.join(&entry.case_id);
        let remote_case = server.remote_root.join(&entry.case_id);
        let mut request = SolverRunRequest {
            adapter: adapter.clone(),
            case_directory: remote_case.clone(),
            mpi_ranks: 1,
            threads: 1,
            ..Default::default()
        };
        if server.runtime == CampaignServerRuntime::Docker {
            request.runtime = SolverRuntime::Docker;
            request.docker_image = server.docker_image.clone();
        }
        let solver_plan = build_solver_command_plan(&request)?;
        plan.commands.push(make_command(
            server,
            "create remote case directory".into(),
            vec!["mkdir".into(), "-p".into(), os(&remote_case)],
        ));
        plan.commands
            .push(make_rsync_command(server, &local_case, &remote_case));
        plan.commands.push(make_command(
            server,
            format!("run {} on {}", entry.case_id, server.name),
            solver_plan.argv.clone(),
        ));
        plan.assignments.push(ServerCaseAssignment {
            case_id: entry.case_id.clone(),
            server_name: server.name.clone(),
            local_case_dir: local_case,
            remote_case_dir: remote_case,
            solver_plan,
        });
    }
    Ok(plan)
}

pub fn multiserver_assignments_tsv(plan: &MultiServerCampaignPlan) -> String {
    let mut out =
        String::from("case_id\tserver\tlocal_case_dir\tremote_case_dir\truntime\tcommand\n");
    for a in &plan.assignments {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            a.case_id,
            a.server_name,
            path_text(&a.local_case_dir),
            path_text(&a.remote_case_dir),
            a.solver_plan.runtime.ordinal(),
            a.solver_plan.display_command
        ));
    }
    out
}

pub fn multiserver_commands_script(plan: &MultiServerCampaignPlan) -> String {
    let mut out = String::from("#!/usr/bin/env sh\nset -eu\n");
    for c in &plan.commands {
        out.push_str(&format!("\n# {}\n{}\n", c.description, c.display_command));
    }
    out
}

fn ensure_absent(paths: &[&Path], what: &str) -> Result<()> {
    for p in paths {
        if p.exists() {
            return Err(SolverError::Io(format!(
                "{what} already exists: {}",
                p.display()
            )));
        }
    }
    Ok(())
}

pub fn write_multiserver_plan_files(
    root: &Path,
    plan: &MultiServerCampaignPlan,
    options: &MultiServerPlanWriteOptions,
) -> Result<MultiServerPlanWriteResult> {
    require(
        !root.as_os_str().is_empty()
            && !options.assignments_filename.is_empty()
            && !options.commands_filename.is_empty(),
        "multi-server plan writer requires root and filenames",
    )?;
    create_dirs(root)?;
    let assignments_path = root.join(&options.assignments_filename);
    let commands_path = root.join(&options.commands_filename);
    if !options.overwrite && (assignments_path.exists() || commands_path.exists()) {
        return Err(SolverError::Io(
            "multi-server plan files already exist".into(),
        ));
    }
    write_file(
        &assignments_path,
        &multiserver_assignments_tsv(plan),
        "failed to write multi-server assignments",
    )?;
    write_file(
        &commands_path,
        &multiserver_commands_script(plan),
        "failed to write multi-server commands",
    )?;
    add_owner_exec(&commands_path)?;
    Ok(MultiServerPlanWriteResult {
        assignments_path,
        commands_path,
        assignments_written: plan.assignments.len(),
        commands_written: plan.commands.len(),
    })
}

pub fn remote_job_state_name(state: RemoteJobState) -> &'static str {
    match state {
        RemoteJobState::Planned => "planned",
        RemoteJobState::Submitted => "submitted",
        RemoteJobState::Running => "running",
        RemoteJobState::Done => "done",
        RemoteJobState::Failed => "failed",
        RemoteJobState::Cancelled => "cancelled",
        RemoteJobState::Lost => "lost",
    }
}

pub fn remote_job_state_from_name(name: &str) -> Result<RemoteJobState> {
    Ok(match lower(name).as_str() {
        "planned" => RemoteJobState::Planned,
        "submitted" | "queued" => RemoteJobState::Submitted,
        "running" => RemoteJobState::Running,
        "done" | "complete" | "completed" | "success" => RemoteJobState::Done,
        "failed" | "error" => RemoteJobState::Failed,
        "cancelled" | "canceled" | "stopped" => RemoteJobState::Cancelled,
        "lost" | "missing" => RemoteJobState::Lost,
        _ => return Err(SolverError::InvalidInput("unknown remote job state".into())),
    })
}

fn find_server_by_name<'a>(
    servers: &'a [CampaignServerDescriptor],
    name: &str,
) -> Result<&'a CampaignServerDescriptor> {
    servers.iter().find(|s| s.name == name).ok_or_else(|| {
        SolverError::InvalidInput(format!(
            "assignment references unknown campaign server: {name}"
        ))
    })
}

fn make_fetch_logs_command(
    server: &CampaignServerDescriptor,
    remote_run: &Path,
    case_id: &str,
) -> RemoteCommandPlan {
    let destination = format!(
        "{}/",
        path_text(&Path::new("remote_logs").join(&server.name).join(case_id))
    );
    let source = format!("{}/", path_text(remote_run));
    let local = is_localhost(&server.host);
    let argv = if local {
        vec![
            "sh".into(),
            "-lc".into(),
            format!(
                "mkdir -p {} && cp -a {}. {}",
                shell_quote(&destination),
                shell_quote(&source),
                shell_quote(&destination)
            ),
        ]
    } else {
        vec![
            "rsync".into(),
            "-az".into(),
            format!("{}:{source}", ssh_target(server)),
            destination,
        ]
    };
    RemoteCommandPlan {
        server_name: server.name.clone(),
        description: format!("fetch logs for {case_id} from {}", server.name),
        display_command: join_display(&argv),
        argv,
        uses_ssh: !local,
    }
}

fn sh_lc(script: String) -> Vec<String> {
    vec!["sh".into(), "-lc".into(), script]
}

/// Health checks per server (root, runtime, disk) and supervised launch/status/cancel/fetch
/// commands per assignment, with logs, pid and exit code under `<case>/<run_dirname>`.
pub fn plan_multiserver_execution(
    campaign_plan: &MultiServerCampaignPlan,
    servers: &[CampaignServerDescriptor],
    run_dirname: &str,
) -> Result<MultiServerExecutionPlan> {
    require(
        !run_dirname.is_empty(),
        "multi-server execution requires a run directory name",
    )?;
    require(
        !campaign_plan.assignments.is_empty(),
        "multi-server execution requires assignments",
    )?;
    require(
        !servers.is_empty(),
        "multi-server execution requires server descriptors",
    )?;
    let mut plan = MultiServerExecutionPlan::default();
    let mut health_added: Vec<String> = Vec::new();
    for a in &campaign_plan.assignments {
        let server = find_server_by_name(servers, &a.server_name)?;
        validate_campaign_server_descriptor(server)?;
        let root_q = shell_quote(&path_text(&server.remote_root));
        if !health_added.contains(&server.name) {
            health_added.push(server.name.clone());
            plan.health_checks.push(RemoteHostHealthCheck {
                server_name: server.name.clone(),
                check_name: "remote-root".into(),
                command: make_command(
                    server,
                    "ensure remote campaign root".into(),
                    sh_lc(format!("mkdir -p {root_q} && test -d {root_q}")),
                ),
                required: true,
            });
            let docker = server.runtime == CampaignServerRuntime::Docker;
            let command = if docker {
                make_command(
                    server,
                    "verify docker runtime".into(),
                    sh_lc("command -v docker >/dev/null && docker --version >/dev/null".into()),
                )
            } else {
                let bin = a.solver_plan.argv.first().cloned().unwrap_or_default();
                make_command(
                    server,
                    "verify native solver binary".into(),
                    sh_lc(format!("command -v {} >/dev/null", shell_quote(&bin))),
                )
            };
            plan.health_checks.push(RemoteHostHealthCheck {
                server_name: server.name.clone(),
                check_name: if docker {
                    "docker-runtime"
                } else {
                    "native-runtime"
                }
                .into(),
                command,
                required: true,
            });
            plan.health_checks.push(RemoteHostHealthCheck {
                server_name: server.name.clone(),
                check_name: "disk-space".into(),
                command: make_command(
                    server,
                    "show remote disk space".into(),
                    sh_lc(format!("df -P {root_q} | tail -n 1")),
                ),
                required: false,
            });
        }
        let run_dir = a.remote_case_dir.join(run_dirname);
        let stdout_path = run_dir.join("solver.stdout.log");
        let stderr_path = run_dir.join("solver.stderr.log");
        let pid_path = run_dir.join("solver.pid");
        let exit_code_path = run_dir.join("exit_code");
        let q = |p: &Path| shell_quote(&path_text(p));
        let solver_cmd = &a.solver_plan.display_command;
        let launch = format!(
            "mkdir -p {} && rm -f {} && cd {} && nohup sh -lc {} > {} 2> {} & echo $! > {}",
            q(&run_dir),
            q(&exit_code_path),
            q(&a.remote_case_dir),
            shell_quote(&format!("{solver_cmd}; echo $? > {}", q(&exit_code_path))),
            q(&stdout_path),
            q(&stderr_path),
            q(&pid_path)
        );
        let status = format!(
            "state=lost; ec=0; if [ -f {ec} ]; then ec=$(cat {ec}); if [ \"$ec\" = 0 ]; then state=done; else state=failed; fi; elif [ -f {pid} ] && kill -0 $(cat {pid}) 2>/dev/null; then state=running; else state=lost; fi; printf '%s\\t%s\\t%s\\t%s\\n' {case} {srv} \"$state\" \"$ec\"",
            ec = q(&exit_code_path),
            pid = q(&pid_path),
            case = shell_quote(&a.case_id),
            srv = shell_quote(&server.name)
        );
        let cancel = format!(
            "if [ -f {pid} ]; then kill $(cat {pid}) 2>/dev/null || true; fi; echo cancelled > {} ; true",
            shell_quote(&format!("{}/cancelled", path_text(&run_dir))),
            pid = q(&pid_path)
        );
        let job = RemoteSupervisedJobPlan {
            case_id: a.case_id.clone(),
            server_name: server.name.clone(),
            remote_case_dir: a.remote_case_dir.clone(),
            launch_command: make_command(
                server,
                format!("launch supervised {}", a.case_id),
                sh_lc(launch),
            ),
            status_command: make_command(server, format!("status {}", a.case_id), sh_lc(status)),
            cancel_command: make_command(server, format!("cancel {}", a.case_id), sh_lc(cancel)),
            fetch_logs_command: make_fetch_logs_command(server, &run_dir, &a.case_id),
            solver_plan: a.solver_plan.clone(),
            run_dir,
            stdout_path,
            stderr_path,
            pid_path,
            exit_code_path,
        };
        if server.runtime == CampaignServerRuntime::Docker {
            plan.docker_jobs += 1;
        } else {
            plan.native_jobs += 1;
        }
        if !is_localhost(&server.host) {
            plan.remote_jobs += 1;
        }
        plan.jobs.push(job);
    }
    Ok(plan)
}

/// `#!/usr/bin/env sh` script with one `# [server] description` comment per command.
pub fn command_script(title: &str, commands: &[&RemoteCommandPlan]) -> String {
    let mut out = format!("#!/usr/bin/env sh\nset -eu\n\n# {title}\n");
    for c in commands {
        out.push_str(&format!(
            "\n# [{}] {}\n{}\n",
            c.server_name, c.description, c.display_command
        ));
    }
    out
}

fn write_command_script(
    path: &Path,
    title: &str,
    commands: &[&RemoteCommandPlan],
    what: &str,
) -> Result<()> {
    write_file(path, &command_script(title, commands), what)?;
    add_owner_exec(path)
}

pub fn multiserver_jobs_tsv(plan: &MultiServerExecutionPlan) -> String {
    let mut out = String::from(
        "case_id\tserver\tremote_case_dir\trun_dir\tstdout\tstderr\tpid\texit_code\truntime\tlaunch\tstatus\tcancel\tfetch_logs\n",
    );
    for j in &plan.jobs {
        let fields = [
            j.case_id.clone(),
            j.server_name.clone(),
            path_text(&j.remote_case_dir),
            path_text(&j.run_dir),
            path_text(&j.stdout_path),
            path_text(&j.stderr_path),
            path_text(&j.pid_path),
            path_text(&j.exit_code_path),
            j.solver_plan.runtime.ordinal().to_string(),
            j.launch_command.display_command.clone(),
            j.status_command.display_command.clone(),
            j.cancel_command.display_command.clone(),
            j.fetch_logs_command.display_command.clone(),
        ];
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    out
}

pub fn write_multiserver_execution_files(
    root: &Path,
    plan: &MultiServerExecutionPlan,
    options: &MultiServerExecutionWriteOptions,
) -> Result<MultiServerExecutionWriteResult> {
    require(
        !root.as_os_str().is_empty(),
        "execution writer requires a root path",
    )?;
    require(
        !options.health_filename.is_empty()
            && !options.launch_filename.is_empty()
            && !options.status_filename.is_empty()
            && !options.fetch_logs_filename.is_empty()
            && !options.jobs_filename.is_empty(),
        "execution writer requires filenames",
    )?;
    create_dirs(root)?;
    let result = MultiServerExecutionWriteResult {
        health_script_path: root.join(&options.health_filename),
        launch_script_path: root.join(&options.launch_filename),
        status_script_path: root.join(&options.status_filename),
        fetch_logs_script_path: root.join(&options.fetch_logs_filename),
        jobs_path: root.join(&options.jobs_filename),
        health_checks_written: plan.health_checks.len(),
        jobs_written: plan.jobs.len(),
    };
    if !options.overwrite {
        ensure_absent(
            &[
                &result.health_script_path,
                &result.launch_script_path,
                &result.status_script_path,
                &result.fetch_logs_script_path,
                &result.jobs_path,
            ],
            "execution file",
        )?;
    }
    let health: Vec<&RemoteCommandPlan> = plan.health_checks.iter().map(|c| &c.command).collect();
    let launch: Vec<&RemoteCommandPlan> = plan.jobs.iter().map(|j| &j.launch_command).collect();
    let status: Vec<&RemoteCommandPlan> = plan.jobs.iter().map(|j| &j.status_command).collect();
    let fetch: Vec<&RemoteCommandPlan> = plan.jobs.iter().map(|j| &j.fetch_logs_command).collect();
    let what = "failed to write execution script";
    write_command_script(
        &result.health_script_path,
        "multi-server health checks",
        &health,
        what,
    )?;
    write_command_script(
        &result.launch_script_path,
        "multi-server supervised launch",
        &launch,
        what,
    )?;
    write_command_script(
        &result.status_script_path,
        "multi-server status probes",
        &status,
        what,
    )?;
    write_command_script(
        &result.fetch_logs_script_path,
        "multi-server log fetch",
        &fetch,
        what,
    )?;
    write_file(
        &result.jobs_path,
        &multiserver_jobs_tsv(plan),
        "failed to write multi-server jobs table",
    )?;
    Ok(result)
}

/// Tab-separated `case_id server state [exit_code]` rows (as printed by the status script);
/// an optional `case_id…` header row is skipped.
pub fn parse_remote_job_status_table(text: &str) -> Result<Vec<RemoteJobStatusRecord>> {
    let mut out = Vec::new();
    let mut first = true;
    for line in getline_split(text, b'\n') {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if first && fields.len() >= 3 && fields[0] == "case_id" {
            first = false;
            continue;
        }
        first = false;
        require(fields.len() >= 3, "malformed remote job status row")?;
        out.push(RemoteJobStatusRecord {
            case_id: fields[0].to_string(),
            server_name: fields[1].to_string(),
            state: remote_job_state_from_name(fields[2])?,
            exit_code: fields.get(3).and_then(|f| stoi_prefix(f)).unwrap_or(0),
            raw_state: fields[2].to_string(),
        });
    }
    Ok(out)
}

pub fn apply_remote_job_status_to_registry(
    registry: &mut Vec<CampaignRegistryEntry>,
    records: &[RemoteJobStatusRecord],
) -> Result<CampaignSummary> {
    for r in records {
        let status = match r.state {
            RemoteJobState::Planned | RemoteJobState::Submitted | RemoteJobState::Running => {
                CampaignCaseStatus::Running
            }
            RemoteJobState::Done => CampaignCaseStatus::Done,
            RemoteJobState::Failed | RemoteJobState::Lost => CampaignCaseStatus::Failed,
            RemoteJobState::Cancelled => CampaignCaseStatus::Stopped,
        };
        update_case_status(
            registry,
            CampaignRegistryEntry::new(
                r.case_id.clone(),
                status,
                0.0,
                0,
                format!(
                    "remote {} on {}",
                    remote_job_state_name(r.state),
                    r.server_name
                ),
            ),
        )?;
    }
    Ok(summarize_campaign_registry(registry))
}

/// Replace the value following each marker (up to whitespace, a quote or `;`) by `<redacted>`.
pub fn redact_sensitive_command_display(command: &str, sensitive_markers: &[String]) -> String {
    let mut out: Vec<u8> = command.as_bytes().to_vec();
    const REDACTED: &[u8] = b"<redacted>";
    for marker in sensitive_markers {
        let m = marker.as_bytes();
        if m.is_empty() {
            continue;
        }
        let mut pos = 0;
        while let Some(found) = find_bytes(&out, m, pos) {
            let start = found + m.len();
            let mut end = start;
            while end < out.len()
                && !c_isspace(out[end])
                && !matches!(out[end], b'\'' | b'"' | b';')
            {
                end += 1;
            }
            out.splice(start..end, REDACTED.iter().copied());
            pos = start + REDACTED.len();
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn find_bytes(hay: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    if from > hay.len() || needle.len() > hay.len() - from {
        return None;
    }
    (from..=hay.len() - needle.len()).find(|&i| &hay[i..i + needle.len()] == needle)
}

/// JSON string escaping as upstream: `\\ \" \n \r \t`, other control characters become `?`.
pub fn json_escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 8);
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push('?'),
            c => out.push(c),
        }
    }
    out
}

fn make_retry_command(
    server: &CampaignServerDescriptor,
    job: &RemoteSupervisedJobPlan,
    options: &RemoteSupervisionOptions,
    launch_display: &str,
) -> RemoteCommandPlan {
    let attempts = options.max_attempts.max(1);
    let script = format!(
        "attempt=1; ec=0; while :; do {launch_display} && exit 0; ec=$?; if [ $attempt -ge {attempts} ]; then exit $ec; fi; sleep {}; attempt=$((attempt+1)); done",
        options.retry_delay_seconds
    );
    make_command(
        server,
        format!("retry launch {}", job.case_id),
        sh_lc(script),
    )
}

fn make_tail_command(
    server: &CampaignServerDescriptor,
    description: String,
    log_path: &Path,
    lines: usize,
) -> RemoteCommandPlan {
    make_command(
        server,
        description,
        sh_lc(format!(
            "tail -n {lines} {} 2>/dev/null || true",
            shell_quote(&path_text(log_path))
        )),
    )
}

/// Access probes per server, retry-wrapped (optionally redacted) launches, log tails and
/// dashboard rows (state `planned`) per job.
pub fn plan_multiserver_supervision(
    execution_plan: &MultiServerExecutionPlan,
    servers: &[CampaignServerDescriptor],
    options: &RemoteSupervisionOptions,
) -> Result<MultiServerSupervisionPlan> {
    require(
        !execution_plan.jobs.is_empty(),
        "multi-server supervision requires jobs",
    )?;
    require(
        !servers.is_empty(),
        "multi-server supervision requires server descriptors",
    )?;
    require(options.max_attempts > 0, "retry attempts must be non-zero")?;
    require(
        options.log_tail_lines > 0,
        "log tail line count must be non-zero",
    )?;
    let mut plan = MultiServerSupervisionPlan::default();
    let mut probed: Vec<String> = Vec::new();
    for job in &execution_plan.jobs {
        let server = find_server_by_name(servers, &job.server_name)?;
        validate_campaign_server_descriptor(server)?;
        if !probed.contains(&server.name) {
            probed.push(server.name.clone());
            let root_q = shell_quote(&path_text(&server.remote_root));
            let local = is_localhost(&server.host);
            if options.include_ssh_probe {
                let command = if local {
                    make_command(
                        server,
                        "verify local execution shell".into(),
                        sh_lc("test -d / && echo local-ok".into()),
                    )
                } else {
                    make_command(
                        server,
                        "verify SSH batch access".into(),
                        sh_lc(format!("echo ssh-ok && test -d {root_q}")),
                    )
                };
                plan.access_probes.push(RemoteAccessProbe {
                    server_name: server.name.clone(),
                    probe_name: if local {
                        "local-access"
                    } else {
                        "ssh-batch-access"
                    }
                    .into(),
                    command,
                    required: true,
                });
            }
            if options.include_docker_probe && server.runtime == CampaignServerRuntime::Docker {
                plan.access_probes.push(RemoteAccessProbe {
                    server_name: server.name.clone(),
                    probe_name: "docker-daemon".into(),
                    command: make_command(
                        server,
                        "probe Docker daemon".into(),
                        sh_lc("docker info >/dev/null && docker ps --format '{{.ID}}' | head -n 1 >/dev/null".into()),
                    ),
                    required: true,
                });
            }
            plan.access_probes.push(RemoteAccessProbe {
                server_name: server.name.clone(),
                probe_name: "campaign-root-writable".into(),
                command: make_command(
                    server,
                    "verify campaign root is writable".into(),
                    sh_lc(format!("mkdir -p {root_q} && test -w {root_q}")),
                ),
                required: true,
            });
        }
        let mut launch_display = job.launch_command.display_command.clone();
        if options.redact_sensitive_values {
            let redacted =
                redact_sensitive_command_display(&launch_display, &options.sensitive_markers);
            if redacted != launch_display {
                plan.sensitive_commands_redacted += 1;
            }
            launch_display = redacted;
        }
        plan.retry_launches.push(RemoteRetryCommandPlan {
            case_id: job.case_id.clone(),
            server_name: job.server_name.clone(),
            max_attempts: options.max_attempts,
            retry_delay_seconds: options.retry_delay_seconds,
            command: make_retry_command(server, job, options, &launch_display),
        });
        plan.log_tails.push(RemoteLogTailPlan {
            case_id: job.case_id.clone(),
            server_name: job.server_name.clone(),
            stdout_tail_command: make_tail_command(
                server,
                format!("tail stdout for {}", job.case_id),
                &job.stdout_path,
                options.log_tail_lines,
            ),
            stderr_tail_command: make_tail_command(
                server,
                format!("tail stderr for {}", job.case_id),
                &job.stderr_path,
                options.log_tail_lines,
            ),
        });
        plan.dashboard_cases.push(RemoteDashboardCase {
            case_id: job.case_id.clone(),
            server_name: job.server_name.clone(),
            state: "planned".into(),
            exit_code: 0,
            stdout_path: job.stdout_path.clone(),
            stderr_path: job.stderr_path.clone(),
        });
    }
    Ok(plan)
}

fn dashboard_case_rows(cases: &[RemoteDashboardCase]) -> String {
    let mut out = String::new();
    for (i, c) in cases.iter().enumerate() {
        out.push_str(&format!(
            "    {{\"case_id\": \"{}\", \"server\": \"{}\", \"state\": \"{}\", \"exit_code\": {}, \"stdout\": \"{}\", \"stderr\": \"{}\"}}",
            json_escape(&c.case_id),
            json_escape(&c.server_name),
            json_escape(&c.state),
            c.exit_code,
            json_escape(&path_text(&c.stdout_path)),
            json_escape(&path_text(&c.stderr_path))
        ));
        if i + 1 < cases.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out
}

pub fn build_multiserver_dashboard_json(cases: &[RemoteDashboardCase]) -> String {
    format!(
        "{{\n  \"schema\": \"cfd.multiserver.dashboard.v1\",\n  \"cases\": [\n{}  ]\n}}\n",
        dashboard_case_rows(cases)
    )
}

pub fn multiserver_supervision_tsv(plan: &MultiServerSupervisionPlan) -> String {
    let mut out = String::from("case_id\tserver\tstate\tstdout\tstderr\n");
    for c in &plan.dashboard_cases {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            c.case_id,
            c.server_name,
            c.state,
            path_text(&c.stdout_path),
            path_text(&c.stderr_path)
        ));
    }
    out
}

pub fn write_multiserver_supervision_files(
    root: &Path,
    plan: &MultiServerSupervisionPlan,
    options: &MultiServerSupervisionWriteOptions,
) -> Result<MultiServerSupervisionWriteResult> {
    require(
        !root.as_os_str().is_empty(),
        "supervision writer requires a root path",
    )?;
    require(
        !options.access_filename.is_empty()
            && !options.retry_launch_filename.is_empty()
            && !options.tail_logs_filename.is_empty()
            && !options.dashboard_filename.is_empty()
            && !options.summary_filename.is_empty(),
        "supervision writer requires filenames",
    )?;
    create_dirs(root)?;
    let result = MultiServerSupervisionWriteResult {
        access_script_path: root.join(&options.access_filename),
        retry_launch_script_path: root.join(&options.retry_launch_filename),
        tail_logs_script_path: root.join(&options.tail_logs_filename),
        dashboard_json_path: root.join(&options.dashboard_filename),
        summary_path: root.join(&options.summary_filename),
        access_checks_written: plan.access_probes.len(),
        retry_launches_written: plan.retry_launches.len(),
        log_tail_commands_written: plan.log_tails.len() * 2,
        dashboard_cases_written: plan.dashboard_cases.len(),
    };
    if !options.overwrite {
        ensure_absent(
            &[
                &result.access_script_path,
                &result.retry_launch_script_path,
                &result.tail_logs_script_path,
                &result.dashboard_json_path,
                &result.summary_path,
            ],
            "supervision file",
        )?;
    }
    let access: Vec<&RemoteCommandPlan> = plan.access_probes.iter().map(|p| &p.command).collect();
    let retry: Vec<&RemoteCommandPlan> = plan.retry_launches.iter().map(|r| &r.command).collect();
    let tail: Vec<&RemoteCommandPlan> = plan
        .log_tails
        .iter()
        .flat_map(|t| [&t.stdout_tail_command, &t.stderr_tail_command])
        .collect();
    let what = "failed to write supervision script";
    write_command_script(
        &result.access_script_path,
        "multi-server access and runtime probes",
        &access,
        what,
    )?;
    write_command_script(
        &result.retry_launch_script_path,
        "multi-server retry launch commands",
        &retry,
        what,
    )?;
    write_command_script(
        &result.tail_logs_script_path,
        "multi-server log tail commands",
        &tail,
        what,
    )?;
    write_file(
        &result.dashboard_json_path,
        &build_multiserver_dashboard_json(&plan.dashboard_cases),
        "failed to write dashboard json",
    )?;
    write_file(
        &result.summary_path,
        &multiserver_supervision_tsv(plan),
        "failed to write supervision summary",
    )?;
    Ok(result)
}

pub fn build_multiserver_controller_status_json(cases: &[RemoteDashboardCase]) -> String {
    let (mut planned, mut running, mut done, mut failed, mut cancelled, mut other) =
        (0, 0, 0, 0, 0, 0);
    for c in cases {
        match lower(&c.state).as_str() {
            "planned" => planned += 1,
            "running" => running += 1,
            "done" => done += 1,
            "failed" => failed += 1,
            "cancelled" | "stopped" => cancelled += 1,
            _ => other += 1,
        }
    }
    format!(
        "{{\n  \"schema\": \"cfd.multiserver.controller.status.v1\",\n  \"summary\": {{\"cases\": {}, \"planned\": {planned}, \"running\": {running}, \"done\": {done}, \"failed\": {failed}, \"cancelled\": {cancelled}, \"other\": {other}}},\n  \"cases\": [\n{}  ]\n}}\n",
        cases.len(),
        dashboard_case_rows(cases)
    )
}

pub fn build_multiserver_controller_events_ndjson(cases: &[RemoteDashboardCase]) -> String {
    let mut out = format!(
        "{{\"event\":\"summary\",\"schema\":\"cfd.multiserver.controller.events.v1\",\"cases\":{}}}\n",
        cases.len()
    );
    for c in cases {
        out.push_str(&format!(
            "{{\"event\":\"case\",\"case_id\":\"{}\",\"server\":\"{}\",\"state\":\"{}\",\"exit_code\":{}}}\n",
            json_escape(&c.case_id),
            json_escape(&c.server_name),
            json_escape(&c.state),
            c.exit_code
        ));
    }
    out
}

pub fn build_multiserver_controller_dashboard_html(config: &RemoteControllerConfig) -> String {
    let mut html = String::new();
    html.push_str("<!doctype html>\n<html lang=\"en\">\n<head>\n");
    html.push_str("  <meta charset=\"utf-8\">\n  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
    html.push_str("  <title>");
    html.push_str(&json_escape(&config.dashboard_title));
    html.push_str("</title>\n");
    html.push_str("  <link rel=\"stylesheet\" href=\"/dashboard.css\">\n</head>\n<body>\n");
    html.push_str("  <main data-refresh-seconds=\"");
    html.push_str(&config.dashboard_refresh_seconds.to_string());
    html.push_str("\">\n");
    html.push_str("    <header><div><p class=\"eyebrow\">cfd_solvers controller</p><h1>");
    html.push_str(&json_escape(&config.dashboard_title));
    html.push_str("</h1></div><span id=\"health\">offline</span></header>\n");
    html.push_str("    <section id=\"summary\" class=\"cards\"></section>\n");
    html.push_str("    <section><h2>Cases</h2><table><thead><tr><th>Case</th><th>Server</th><th>State</th><th>Exit</th><th>Logs</th><th>Control</th></tr></thead><tbody id=\"cases\"></tbody></table></section>\n");
    html.push_str("    <section><h2>Log preview</h2><pre id=\"log\">Select stdout or stderr.</pre></section>\n");
    html.push_str("  </main>\n  <script src=\"/dashboard.js\"></script>\n</body>\n</html>\n");
    html
}

const DASHBOARD_JS_BODY: &str = r##"
async function getJson(path) {
  const r = await fetch(path, {cache: 'no-store', headers: authHeaders()});
  if (!r.ok) throw new Error(path + ' -> ' + r.status);
  return await r.json();
}
let controllerToken = sessionStorage.getItem('cfd-controller-token') || '';
function authHeaders() { return controllerToken ? {'Authorization': 'Bearer ' + controllerToken} : {}; }
function ensureToken() {
  if (!controllerToken) controllerToken = window.prompt('Controller access token:', '') || '';
  if (controllerToken) sessionStorage.setItem('cfd-controller-token', controllerToken);
  return controllerToken;
}
function setHealth(text, ok) {
  const el = document.getElementById('health');
  el.textContent = text;
  el.className = ok ? 'ok' : 'bad';
}
let eventSource;
function connectEvents() {
  if (!window.EventSource) return;
  const token = ensureToken();
  eventSource = new EventSource('/api/events/stream' + (token ? '?token=' + encodeURIComponent(token) : ''));
  eventSource.addEventListener('status', event => {
    const status = JSON.parse(event.data);
    setHealth('live', true);
    renderSummary(status.summary || {});
    renderCases(status.cases || []);
  });
  eventSource.onerror = () => setHealth('reconnecting', false);
}
function renderSummary(summary) {
  const root = document.getElementById('summary');
  const keys = ['cases','planned','running','done','failed','cancelled','other'];
  root.innerHTML = keys.map(k => `<article><span>${k}</span><strong>${summary[k] ?? 0}</strong></article>`).join('');
}
async function showLog(caseId, stream) {
  ensureToken();
  const r = await fetch(`/api/cases/${encodeURIComponent(caseId)}/logs/${stream}`, {cache: 'no-store', headers: authHeaders()});
  document.getElementById('log').textContent = await r.text();
}
async function sendControl(caseId, action) {
  ensureToken();
  const r = await fetch(`/api/cases/${encodeURIComponent(caseId)}/control/${action}`, {method: 'POST', headers: authHeaders()});
  const text = await r.text();
  if (!r.ok) window.alert(text);
  await refresh();
}
function renderCases(cases) {
  const body = document.getElementById('cases');
  body.innerHTML = cases.map(c => {
    const logs = `<button onclick="showLog('${c.case_id}','stdout')">stdout</button> <button onclick="showLog('${c.case_id}','stderr')">stderr</button>`;
    const control = controlActions.map(a => `<button onclick="sendControl('${c.case_id}','${a}')">${a}</button>`).join(' ');
    return `<tr><td>${c.case_id}</td><td>${c.server}</td><td><span class="state ${c.state}">${c.state}</span></td><td>${c.exit_code}</td><td>${logs}</td><td>${control}</td></tr>`;
  }).join('');
}
async function refresh() {
  try {
    const health = await getJson('/api/health');
    const status = await getJson('/api/status');
    setHealth(health.ok ? 'online' : 'degraded', !!health.ok);
    renderSummary(status.summary || {});
    renderCases(status.cases || []);
  } catch (err) {
    setHealth('offline: ' + err.message, false);
  }
}
refresh();
setInterval(refresh, refreshSeconds * 1000);
connectEvents();
"##;

pub fn build_multiserver_controller_dashboard_js(config: &RemoteControllerConfig) -> String {
    let actions: Vec<String> = config
        .control_actions
        .iter()
        .map(|a| format!("'{}'", json_escape(a)))
        .collect();
    format!(
        "const refreshSeconds = {};\nconst controlActions = [{}];\n{DASHBOARD_JS_BODY}",
        config.dashboard_refresh_seconds,
        actions.join(",")
    )
}

const DASHBOARD_CSS: &str = r##":root { color-scheme: light dark; font-family: Inter, ui-sans-serif, system-ui, sans-serif; }
body { margin: 0; background: #0f172a; color: #e2e8f0; }
main { max-width: 1180px; margin: 0 auto; padding: 28px; }
header { display: flex; justify-content: space-between; align-items: center; gap: 16px; margin-bottom: 24px; }
h1 { margin: 0; font-size: clamp(1.8rem, 4vw, 3rem); }
.eyebrow { margin: 0 0 4px; color: #93c5fd; text-transform: uppercase; letter-spacing: .14em; font-size: .78rem; }
#health { padding: 8px 12px; border-radius: 999px; background: #334155; }
#health.ok { background: #14532d; }
#health.bad { background: #7f1d1d; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(130px, 1fr)); gap: 12px; margin-bottom: 24px; }
article { background: #1e293b; border: 1px solid #334155; border-radius: 14px; padding: 14px; }
article span { color: #94a3b8; display: block; }
article strong { display: block; font-size: 2rem; margin-top: 6px; }
section { background: #111827; border: 1px solid #334155; border-radius: 16px; padding: 18px; margin-bottom: 20px; overflow-x: auto; }
table { width: 100%; border-collapse: collapse; }
th, td { padding: 10px; border-bottom: 1px solid #334155; text-align: left; }
button { border: 1px solid #64748b; background: #1e293b; color: #e2e8f0; border-radius: 8px; padding: 6px 10px; margin: 2px; cursor: pointer; }
button:hover { background: #334155; }
.state { padding: 3px 8px; border-radius: 999px; background: #334155; }
.state.running { background: #1d4ed8; }
.state.done { background: #166534; }
.state.failed { background: #991b1b; }
pre { min-height: 160px; white-space: pre-wrap; background: #020617; padding: 14px; border-radius: 12px; }
"##;

pub fn build_multiserver_controller_dashboard_css() -> String {
    DASHBOARD_CSS.to_string()
}

pub fn build_multiserver_controller_openapi_json(
    routes: &[RemoteControllerRoute],
    config: &RemoteControllerConfig,
) -> String {
    let mut out = format!(
        "{{\n  \"openapi\": \"3.0.3\",\n  \"info\": {{\"title\": \"cfd_solvers multi-server controller\", \"version\": \"0.12.0\"}},\n  \"servers\": [{{\"url\": \"http://{}:{}\"}}],\n  \"paths\": {{\n",
        json_escape(&config.bind_host),
        config.port
    );
    for (i, r) in routes.iter().enumerate() {
        out.push_str(&format!(
            "    \"{}\": {{\"{}\": {{\"summary\": \"{}\", \"x-mutates\": {}, \"x-requires-token\": {}}}}}",
            json_escape(&r.path),
            lower(&r.method),
            json_escape(&r.purpose),
            r.mutates,
            r.requires_token
        ));
        if i + 1 < routes.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  }\n}\n");
    out
}

/// The dependency-free Python controller (`http.server`, token roles, TLS, SSE, probes).
pub fn build_multiserver_controller_script(config: &RemoteControllerConfig) -> String {
    let mut py = String::new();
    py.push_str("#!/usr/bin/env python3\n");
    py.push_str(
        "import hashlib, hmac, json, os, pathlib, ssl, subprocess, sys, time, urllib.parse\n",
    );
    py.push_str("from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer\n\n");
    py.push_str("HERE = pathlib.Path(__file__).resolve().parent\n");
    py.push_str("CAMPAIGN_ROOT = pathlib.Path(os.environ.get('");
    py.push_str(&config.campaign_root_env_var);
    py.push_str("', HERE.parent)).resolve()\n");
    py.push_str("TOKEN_ENV = '");
    py.push_str(&config.token_env_var);
    py.push_str("'\n");
    py.push_str("TAIL_LINES = ");
    py.push_str(&config.log_tail_lines.to_string());
    py.push('\n');
    py.push_str("EVENT_HEARTBEAT_SECONDS = ");
    py.push_str(&config.event_stream_heartbeat_seconds.to_string());
    py.push('\n');
    py.push_str("TOKEN_FILE_ENV = '");
    py.push_str(&config.token_file_env_var);
    py.push_str("'\n");
    py.push_str("REQUIRE_AUTH = ");
    py.push_str(if config.require_authentication {
        "True"
    } else {
        "False"
    });
    py.push('\n');
    py.push_str("TLS_ENABLED = ");
    py.push_str(if config.enable_tls { "True" } else { "False" });
    py.push('\n');
    py.push_str("TLS_CERT_ENV = '");
    py.push_str(&config.tls_cert_env_var);
    py.push_str("'\n");
    py.push_str("TLS_KEY_ENV = '");
    py.push_str(&config.tls_key_env_var);
    py.push_str("'\n");
    py.push_str("HISTORY_ENABLED = ");
    py.push_str(if config.enable_event_history {
        "True"
    } else {
        "False"
    });
    py.push('\n');
    py.push_str("MAX_HISTORY_BYTES = ");
    py.push_str(&config.max_event_history_bytes.to_string());
    py.push('\n');
    py.push_str("LIVE_PROBES_ENABLED = ");
    py.push_str(if config.enable_live_probes {
        "True"
    } else {
        "False"
    });
    py.push('\n');
    py.push_str("PROBE_TIMEOUT_SECONDS = ");
    py.push_str(&config.probe_timeout_seconds.to_string());
    py.push('\n');
    py.push_str("ALLOW_MUTATION = ");
    py.push_str(if config.allow_mutating_actions {
        "True"
    } else {
        "False"
    });
    py.push('\n');
    py.push_str("REQUIRE_TOKEN = ");
    py.push_str(if config.require_token_for_mutation {
        "True"
    } else {
        "False"
    });
    py.push('\n');
    py.push_str("ACTIONS = set(");
    py.push('[');
    py.push_str(
        &config
            .control_actions
            .iter()
            .map(|a| format!("'{a}'"))
            .collect::<Vec<_>>()
            .join(", "),
    );
    py.push(']');
    py.push_str(")\n\n");
    py.push_str("def read_status():\n    path = HERE / 'status.json'\n    if not path.exists():\n        return {'schema': 'cfd.multiserver.controller.status.v1', 'summary': {'cases': 0}, 'cases': []}\n    return json.loads(path.read_text(encoding='utf-8'))\n\ndef append_history(event, payload):\n    if not HISTORY_ENABLED:\n        return\n    path = HERE / 'events-history.ndjson'\n    record = json.dumps({'time_ns': time.time_ns(), 'event': event, 'payload': payload}, separators=(',', ':')) + '\\n'\n    if path.exists() and path.stat().st_size + len(record.encode('utf-8')) > MAX_HISTORY_BYTES:\n        rotated = HERE / 'events-history.ndjson.1'\n        path.replace(rotated)\n    with path.open('a', encoding='utf-8') as stream:\n        stream.write(record)\n\n");
    py.push_str("def reply(handler, code, obj, content_type='application/json'):\n    data = json.dumps(obj, indent=2).encode('utf-8') if not isinstance(obj, (bytes, bytearray)) else obj\n    handler.send_response(code)\n    handler.send_header('Content-Type', content_type)\n    handler.send_header('Content-Length', str(len(data)))\n    handler.end_headers()\n    handler.wfile.write(data)\n\n");
    py.push_str("def authorized(handler, required_role='viewer'):\n    if not REQUIRE_AUTH and not REQUIRE_TOKEN:\n        return True\n    header = handler.headers.get('Authorization', '')\n    query = urllib.parse.parse_qs(urllib.parse.urlparse(handler.path).query)\n    raw = header[7:] if header.startswith('Bearer ') else query.get('token', [''])[0]\n    legacy = os.environ.get(TOKEN_ENV, '')\n    if legacy and hmac.compare_digest(raw, legacy):\n        return True\n    token_path = pathlib.Path(os.environ.get(TOKEN_FILE_ENV, HERE / 'controller.tokens.json'))\n    try:\n        entries = json.loads(token_path.read_text(encoding='utf-8')).get('tokens', [])\n    except (OSError, ValueError):\n        return False\n    digest = hashlib.sha256(raw.encode('utf-8')).hexdigest()\n    levels = {'viewer': 0, 'operator': 1, 'admin': 2}\n    for entry in entries:\n        if hmac.compare_digest(digest, str(entry.get('sha256', ''))) and levels.get(entry.get('role'), -1) >= levels[required_role]:\n            return True\n    return False\n\n");
    py.push_str("def tail(path_text, lines):\n    path = pathlib.Path(path_text)\n    if not path.is_absolute():\n        path = CAMPAIGN_ROOT / path\n    try:\n        buf = path.read_text(errors='replace').splitlines()\n        return '\\n'.join(buf[-lines:]) + ('\\n' if buf else '')\n    except OSError as exc:\n        return 'log unavailable: %s\\n' % exc\n\n");
    py.push_str("class Handler(BaseHTTPRequestHandler):\n    def do_GET(self):\n        parsed = urllib.parse.urlparse(self.path)\n        parts = [p for p in parsed.path.split('/') if p]\n        status = read_status()\n        if parsed.path in ('/', '/index.html'):\n            data = (HERE / 'index.html').read_bytes()\n            reply(self, 200, data, 'text/html; charset=utf-8')\n            return\n        if parsed.path == '/dashboard.js':\n            data = (HERE / 'dashboard.js').read_bytes()\n            reply(self, 200, data, 'application/javascript; charset=utf-8')\n            return\n        if parsed.path == '/dashboard.css':\n            data = (HERE / 'dashboard.css').read_bytes()\n            reply(self, 200, data, 'text/css; charset=utf-8')\n            return\n        if parsed.path == '/api/events':\n            path = HERE / 'events.ndjson'\n            data = path.read_bytes() if path.exists() else b''\n            reply(self, 200, data, 'application/x-ndjson; charset=utf-8')\n            return\n        if parsed.path == '/api/health':\n            reply(self, 200, {'ok': True, 'schema': 'cfd.multiserver.controller.health.v1'})\n            return\n        if parsed.path == '/api/status':\n            reply(self, 200, status)\n            return\n        if parsed.path == '/api/cases':\n            reply(self, 200, {'cases': status.get('cases', [])})\n            return\n        if len(parts) == 5 and parts[:2] == ['api', 'cases'] and parts[3] == 'logs':\n            case_id, stream = parts[2], parts[4]\n            for case in status.get('cases', []):\n                if case.get('case_id') == case_id and stream in ('stdout', 'stderr'):\n                    data = tail(case.get(stream, ''), TAIL_LINES).encode('utf-8')\n                    reply(self, 200, data, 'text/plain; charset=utf-8')\n                    return\n            reply(self, 404, {'error': 'case or log stream not found'})\n            return\n        reply(self, 404, {'error': 'not found'})\n\n    def do_POST(self):\n        parsed = urllib.parse.urlparse(self.path)\n        parts = [p for p in parsed.path.split('/') if p]\n        if len(parts) == 5 and parts[:2] == ['api', 'cases'] and parts[3] == 'control':\n            if not ALLOW_MUTATION:\n                reply(self, 403, {'error': 'mutations disabled'})\n                return\n            if not authorized(self):\n                reply(self, 401, {'error': 'missing or invalid controller token'})\n                return\n            case_id, action = parts[2], parts[4]\n            if action not in ACTIONS:\n                reply(self, 400, {'error': 'unsupported action'})\n                return\n            control_dir = CAMPAIGN_ROOT / case_id / '.cfd_control'\n            control_dir.mkdir(parents=True, exist_ok=True)\n            body = 'case_id=%s\\naction=%s\\nreason=controller-api\\n' % (case_id, action)\n            (control_dir / (action + '.directive')).write_text(body, encoding='utf-8')\n            (control_dir / 'latest.directive').write_text(body, encoding='utf-8')\n            reply(self, 202, {'ok': True, 'case_id': case_id, 'action': action})\n            return\n        reply(self, 404, {'error': 'not found'})\n\n    def log_message(self, fmt, *args):\n        sys.stderr.write('controller: ' + fmt % args + '\\n')\n\n");
    py.push_str("def sse_frame(event, payload, event_id):\n    data = json.dumps(payload, separators=(',', ':'))\n    return ('id: %s\\nevent: %s\\ndata: %s\\n\\n' % (event_id, event, data)).encode('utf-8')\n\n");
    py.push_str("class StreamingHandler(Handler):\n    protocol_version = 'HTTP/1.1'\n\n    def do_GET(self):\n        if urllib.parse.urlparse(self.path).path != '/api/events/stream':\n            return super().do_GET()\n        self.send_response(200)\n        self.send_header('Content-Type', 'text/event-stream; charset=utf-8')\n        self.send_header('Cache-Control', 'no-cache')\n        self.send_header('Connection', 'keep-alive')\n        self.end_headers()\n        event_id = 0\n        last_payload = None\n        try:\n            while True:\n                payload = read_status()\n                encoded = json.dumps(payload, sort_keys=True)\n                if encoded != last_payload:\n                    event_id += 1\n                    self.wfile.write(sse_frame('status', payload, event_id))\n                    last_payload = encoded\n                else:\n                    self.wfile.write(b': heartbeat\\n\\n')\n                self.wfile.flush()\n                time.sleep(EVENT_HEARTBEAT_SECONDS)\n        except (BrokenPipeError, ConnectionResetError):\n            return\n\n");
    py.push_str("class ProductionHandler(StreamingHandler):\n    def reject_auth(self):\n        reply(self, 401, {'error': 'authentication required'})\n\n    def do_GET(self):\n        path = urllib.parse.urlparse(self.path).path\n        if path.startswith('/api/') and path != '/api/health' and not authorized(self, 'viewer'):\n            return self.reject_auth()\n        if path == '/api/events/history':\n            history = HERE / 'events-history.ndjson'\n            return reply(self, 200, history.read_bytes() if history.exists() else b'', 'application/x-ndjson; charset=utf-8')\n        if path == '/api/probes':\n            if not LIVE_PROBES_ENABLED:\n                return reply(self, 403, {'error': 'live probes disabled'})\n            script = HERE.parent / 'multiserver_access_checks.sh'\n            try:\n                run = subprocess.run([str(script)], cwd=HERE.parent, capture_output=True, text=True, timeout=PROBE_TIMEOUT_SECONDS, check=False)\n                payload = {'ok': run.returncode == 0, 'exit_code': run.returncode, 'stdout': run.stdout[-8192:], 'stderr': run.stderr[-8192:]}\n            except (OSError, subprocess.TimeoutExpired) as exc:\n                payload = {'ok': False, 'error': str(exc)}\n            append_history('probe', payload)\n            return reply(self, 200 if payload.get('ok') else 503, payload)\n        return super().do_GET()\n\n    def do_POST(self):\n        if not authorized(self, 'operator'):\n            return self.reject_auth()\n        append_history('control-request', {'path': urllib.parse.urlparse(self.path).path})\n        return super().do_POST()\n\n");
    py.push_str("if __name__ == '__main__':\n    host = os.environ.get('CFD_CONTROLLER_HOST', '");
    py.push_str(&config.bind_host);
    py.push_str("')\n    port = int(os.environ.get('CFD_CONTROLLER_PORT', '");
    py.push_str(&config.port.to_string());
    py.push_str("'))\n    server = ThreadingHTTPServer((host, port), ProductionHandler)\n    if TLS_ENABLED:\n        cert = os.environ.get(TLS_CERT_ENV, '')\n        key = os.environ.get(TLS_KEY_ENV, '')\n        if not cert or not key:\n            raise SystemExit('TLS enabled but certificate/key environment variables are missing')\n        context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)\n        context.minimum_version = ssl.TLSVersion.TLSv1_2\n        context.load_cert_chain(cert, key)\n        server.socket = context.wrap_socket(server.socket, server_side=True)\n    print('serving cfd controller on %s:%d campaign_root=%s tls=%s' % (host, port, CAMPAIGN_ROOT, TLS_ENABLED))\n    server.serve_forever()\n");
    py
}

fn route(
    method: &str,
    path: &str,
    purpose: &str,
    requires_token: bool,
    mutates: bool,
) -> RemoteControllerRoute {
    RemoteControllerRoute {
        method: method.into(),
        path: path.into(),
        purpose: purpose.into(),
        requires_token,
        mutates,
    }
}

const TOKENS_EXAMPLE_JSON: &str = "{\n  \"tokens\": [\n    {\"name\": \"dashboard-viewer\", \"role\": \"viewer\", \"sha256\": \"replace-with-sha256-of-random-token\"},\n    {\"name\": \"campaign-operator\", \"role\": \"operator\", \"sha256\": \"replace-with-sha256-of-another-random-token\"}\n  ]\n}\n";
const SYSTEMD_UNIT: &str = "[Unit]\nDescription=cfd_solvers campaign controller\nAfter=network-online.target\nWants=network-online.target\n\n[Service]\nType=simple\nUser=cfd-controller\nGroup=cfd-controller\nWorkingDirectory=/var/lib/cfd-controller\nEnvironmentFile=/etc/cfd-controller/controller.env\nExecStart=/usr/bin/python3 /var/lib/cfd-controller/controller/cfd_controller.py\nRestart=on-failure\nRestartSec=3\nNoNewPrivileges=true\nPrivateTmp=true\nProtectSystem=strict\nProtectHome=true\nReadWritePaths=/var/lib/cfd-controller\nCapabilityBoundingSet=\nLockPersonality=true\nMemoryDenyWriteExecute=true\n\n[Install]\nWantedBy=multi-user.target\n";
const REVERSE_PROXY_CONFIG: &str = ":443 {\n  tls /etc/cfd-controller/tls.crt /etc/cfd-controller/tls.key\n  encode zstd gzip\n  reverse_proxy 127.0.0.1:8080 {\n    flush_interval -1\n  }\n  header {\n    Strict-Transport-Security \"max-age=31536000; includeSubDomains\"\n    X-Content-Type-Options nosniff\n    X-Frame-Options DENY\n    Referrer-Policy no-referrer\n  }\n}\n";
const START_SCRIPT: &str = "#!/usr/bin/env sh\nset -eu\nSCRIPT_DIR=$(CDPATH= cd -- \"$(dirname -- \"$0\")\" && pwd)\ncd \"$SCRIPT_DIR/controller\"\nexec python3 cfd_controller.py\n";

pub fn plan_multiserver_controller(
    supervision_plan: &MultiServerSupervisionPlan,
    config: &RemoteControllerConfig,
) -> Result<MultiServerControllerPlan> {
    require(config.port != 0, "controller port must be non-zero")?;
    require(
        !config.token_env_var.is_empty(),
        "controller token env var must not be empty",
    )?;
    require(
        config.log_tail_lines != 0,
        "controller log tail line count must be non-zero",
    )?;
    require(
        config.dashboard_refresh_seconds != 0,
        "controller dashboard refresh interval must be non-zero",
    )?;
    require(
        config.event_stream_heartbeat_seconds != 0,
        "controller event heartbeat interval must be non-zero",
    )?;
    require(
        !config.require_authentication || !config.token_file_env_var.is_empty(),
        "controller token file env var must not be empty",
    )?;
    require(
        !config.enable_tls
            || (!config.tls_cert_env_var.is_empty() && !config.tls_key_env_var.is_empty()),
        "controller TLS env vars must not be empty",
    )?;
    require(
        config.max_event_history_bytes != 0,
        "controller event history limit must be non-zero",
    )?;
    require(
        config.probe_timeout_seconds != 0,
        "controller probe timeout must be non-zero",
    )?;
    let mut routes = vec![
        route("GET", "/", "browser dashboard", false, false),
        route(
            "GET",
            "/dashboard.js",
            "dashboard JavaScript asset",
            false,
            false,
        ),
        route(
            "GET",
            "/dashboard.css",
            "dashboard stylesheet asset",
            false,
            false,
        ),
        route(
            "GET",
            "/api/health",
            "controller health check",
            false,
            false,
        ),
        route(
            "GET",
            "/api/status",
            "summary plus case status",
            false,
            false,
        ),
        route("GET", "/api/cases", "case table", false, false),
        route(
            "GET",
            "/api/cases/{case_id}/logs/{stream}",
            "tail stdout or stderr metadata path",
            false,
            false,
        ),
    ];
    if config.expose_event_snapshot {
        routes.push(route(
            "GET",
            "/api/events",
            "newline-delimited dashboard event snapshot",
            false,
            false,
        ));
    }
    if config.expose_event_stream {
        routes.push(route(
            "GET",
            "/api/events/stream",
            "server-sent status event stream",
            false,
            false,
        ));
    }
    if config.enable_event_history {
        routes.push(route(
            "GET",
            "/api/events/history",
            "durable bounded event history",
            config.require_authentication,
            false,
        ));
    }
    if config.enable_live_probes {
        routes.push(route(
            "GET",
            "/api/probes",
            "run supervised SSH Docker and filesystem probes",
            config.require_authentication,
            false,
        ));
    }
    if config.allow_mutating_actions {
        routes.push(route(
            "POST",
            "/api/cases/{case_id}/control/{action}",
            "write a campaign control directive",
            config.require_token_for_mutation,
            true,
        ));
    }
    let cases = supervision_plan.dashboard_cases.clone();
    Ok(MultiServerControllerPlan {
        openapi_json: build_multiserver_controller_openapi_json(&routes, config),
        status_json: build_multiserver_controller_status_json(&cases),
        controller_script: build_multiserver_controller_script(config),
        index_html: build_multiserver_controller_dashboard_html(config),
        dashboard_js: build_multiserver_controller_dashboard_js(config),
        dashboard_css: build_multiserver_controller_dashboard_css(),
        events_ndjson: build_multiserver_controller_events_ndjson(&cases),
        tokens_example_json: TOKENS_EXAMPLE_JSON.into(),
        systemd_unit: SYSTEMD_UNIT.into(),
        reverse_proxy_config: REVERSE_PROXY_CONFIG.into(),
        start_script: START_SCRIPT.into(),
        config: config.clone(),
        routes,
        cases,
    })
}

pub fn controller_routes_tsv(plan: &MultiServerControllerPlan) -> String {
    let mut out = String::from("method\tpath\trequires_token\tmutates\tpurpose\n");
    for r in &plan.routes {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            r.method,
            r.path,
            u8::from(r.requires_token),
            u8::from(r.mutates),
            r.purpose
        ));
    }
    out
}

pub fn controller_env_example(root: &Path, plan: &MultiServerControllerPlan) -> String {
    let c = &plan.config;
    format!(
        "# Store real credentials outside source control.\n{}=/etc/cfd-controller/controller.tokens.json\nCFD_CONTROLLER_HOST={}\nCFD_CONTROLLER_PORT={}\n{}={}\n{}=/etc/cfd-controller/tls.crt\n{}=/etc/cfd-controller/tls.key\n",
        c.token_file_env_var,
        c.bind_host,
        c.port,
        c.campaign_root_env_var,
        path_text(root),
        c.tls_cert_env_var,
        c.tls_key_env_var
    )
}

pub fn controller_readme(
    plan: &MultiServerControllerPlan,
    options: &MultiServerControllerWriteOptions,
) -> String {
    format!(
        "# Multi-server controller\n\nThis directory contains a dependency-free Python controller for generated multi-server campaign supervision artifacts.\n\n\
Start it from the campaign root with `./{start}`. Copy `{tokens}` outside the source tree, replace the placeholder hashes with SHA-256 digests of random tokens, and set `{tokenvar}`. Configure the TLS certificate/key variables before startup.\n\n\
Routes are documented in `{routes}` and `{openapi}`. Status is served from `{status}`; remote workers still use the generated supervision scripts.\n\
Open `{index}` or `http://127.0.0.1:{port}/` for the generated live dashboard. It prefers the SSE status stream and retains polling as a reconnect fallback. `{events}` stores a newline-delimited event snapshot for downstream dashboards.\n",
        start = options.start_script_filename,
        tokens = options.tokens_example_filename,
        tokenvar = plan.config.token_file_env_var,
        routes = options.routes_filename,
        openapi = options.openapi_filename,
        status = options.status_filename,
        index = options.index_filename,
        port = plan.config.port,
        events = options.events_filename
    )
}

pub fn write_multiserver_controller_files(
    root: &Path,
    plan: &MultiServerControllerPlan,
    options: &MultiServerControllerWriteOptions,
) -> Result<MultiServerControllerWriteResult> {
    require(
        !root.as_os_str().is_empty(),
        "controller writer requires root",
    )?;
    let names = [
        &options.controller_dirname,
        &options.controller_script_filename,
        &options.start_script_filename,
        &options.openapi_filename,
        &options.status_filename,
        &options.routes_filename,
        &options.index_filename,
        &options.dashboard_js_filename,
        &options.dashboard_css_filename,
        &options.events_filename,
        &options.tokens_example_filename,
        &options.systemd_filename,
        &options.reverse_proxy_filename,
    ];
    require(
        names.iter().all(|n| !n.is_empty()),
        "controller writer requires filenames",
    )?;
    let dir = root.join(&options.controller_dirname);
    create_dirs(&dir)?;
    let r = MultiServerControllerWriteResult {
        controller_script_path: dir.join(&options.controller_script_filename),
        start_script_path: root.join(&options.start_script_filename),
        openapi_path: dir.join(&options.openapi_filename),
        status_path: dir.join(&options.status_filename),
        routes_path: dir.join(&options.routes_filename),
        env_path: dir.join(&options.env_filename),
        readme_path: dir.join(&options.readme_filename),
        index_path: dir.join(&options.index_filename),
        dashboard_js_path: dir.join(&options.dashboard_js_filename),
        dashboard_css_path: dir.join(&options.dashboard_css_filename),
        events_path: dir.join(&options.events_filename),
        tokens_example_path: dir.join(&options.tokens_example_filename),
        systemd_path: dir.join(&options.systemd_filename),
        reverse_proxy_path: dir.join(&options.reverse_proxy_filename),
        routes_written: plan.routes.len(),
        cases_written: plan.cases.len(),
        static_assets_written: 7,
    };
    if !options.overwrite {
        ensure_absent(
            &[
                &r.controller_script_path,
                &r.start_script_path,
                &r.openapi_path,
                &r.status_path,
                &r.routes_path,
                &r.env_path,
                &r.readme_path,
                &r.index_path,
                &r.dashboard_js_path,
                &r.dashboard_css_path,
                &r.events_path,
            ],
            "controller file",
        )?;
    }
    let what = "failed to write controller file";
    write_file(&r.controller_script_path, &plan.controller_script, what)?;
    write_file(&r.start_script_path, &plan.start_script, what)?;
    write_file(&r.openapi_path, &plan.openapi_json, what)?;
    write_file(&r.status_path, &plan.status_json, what)?;
    write_file(&r.index_path, &plan.index_html, what)?;
    write_file(&r.dashboard_js_path, &plan.dashboard_js, what)?;
    write_file(&r.dashboard_css_path, &plan.dashboard_css, what)?;
    write_file(&r.events_path, &plan.events_ndjson, what)?;
    write_file(&r.tokens_example_path, &plan.tokens_example_json, what)?;
    write_file(&r.systemd_path, &plan.systemd_unit, what)?;
    write_file(&r.reverse_proxy_path, &plan.reverse_proxy_config, what)?;
    write_file(&r.routes_path, &controller_routes_tsv(plan), what)?;
    write_file(&r.env_path, &controller_env_example(root, plan), what)?;
    write_file(&r.readme_path, &controller_readme(plan, options), what)?;
    add_owner_exec(&r.controller_script_path)?;
    add_owner_exec(&r.start_script_path)?;
    Ok(r)
}

fn deploy_file(path: &str, contents: String) -> DockerDeployFile {
    DockerDeployFile {
        relative_path: PathBuf::from(path),
        contents,
    }
}

/// Dockerfile, compose.yaml (optional manager service, replicas, CPU/memory limits, GPU
/// reservation), helper scripts, env and server examples, and a README.
pub fn generate_docker_deploy_files(config: &DockerDeployConfig) -> Result<Vec<DockerDeployFile>> {
    require(
        !config.image_name.is_empty()
            && !config.service_name.is_empty()
            && !config.container_campaign_dir.is_empty()
            && !config.worker_command.is_empty(),
        "docker deploy config requires image, service, campaign dir and command",
    )?;
    let c = config;
    let mut files = Vec::new();
    files.push(deploy_file(
        "Dockerfile",
        format!(
            "FROM debian:bookworm-slim\nRUN apt-get update && apt-get install -y --no-install-recommends ca-certificates libgomp1 && rm -rf /var/lib/apt/lists/*\nWORKDIR /opt/cfd_solvers\nCOPY . /opt/cfd_solvers\nENV PATH=/opt/cfd_solvers/bin:$PATH\nRUN mkdir -p {}\nCMD [\"/bin/sh\",\"-lc\",\"{}\"]\n",
            c.container_campaign_dir, c.worker_command
        ),
    ));
    let env_block = docker_env_block(c, "      ");
    let healthcheck = format!(
        "    healthcheck:\n      test: [\"CMD-SHELL\", \"test -d {}\"]\n      interval: 30s\n      timeout: 5s\n      retries: 3\n",
        c.container_campaign_dir
    );
    let mut compose = String::from("services:\n");
    if c.include_manager_service {
        compose.push_str(&format!(
            "  campaign-manager:\n    image: {}\n    build: .\n    command: [\"/bin/sh\", \"-lc\", \"{}\"]\n    environment:\n{env_block}    volumes:\n      - {}:{}\n    ports:\n      - \"{port}:{port}\"\n{healthcheck}",
            c.image_name,
            c.manager_command,
            c.host_campaign_dir,
            c.container_campaign_dir,
            port = c.dashboard_port
        ));
    }
    compose.push_str(&format!(
        "  {}:\n    image: {}\n    build: .\n    command: [\"/bin/sh\", \"-lc\", \"{}\"]\n    environment:\n{env_block}    volumes:\n      - {}:{}\n{healthcheck}",
        c.service_name, c.image_name, c.worker_command, c.host_campaign_dir, c.container_campaign_dir
    ));
    let limits = c.worker_cpus > 0 || !c.worker_memory.is_empty();
    if c.worker_replicas > 1 || limits || c.enable_gpu_runtime {
        compose.push_str("    deploy:\n");
        if c.worker_replicas > 1 {
            compose.push_str(&format!("      replicas: {}\n", c.worker_replicas));
        }
        if limits || c.enable_gpu_runtime {
            compose.push_str("      resources:\n");
            if limits {
                compose.push_str("        limits:\n");
                if c.worker_cpus > 0 {
                    compose.push_str(&format!("          cpus: '{}'\n", c.worker_cpus));
                }
                if !c.worker_memory.is_empty() {
                    compose.push_str(&format!("          memory: {}\n", c.worker_memory));
                }
            }
            if c.enable_gpu_runtime {
                compose.push_str("        reservations:\n          devices:\n            - capabilities: [gpu]\n");
            }
        }
    }
    files.push(deploy_file("compose.yaml", compose));
    files.push(deploy_file(
        ".dockerignore",
        "build*/\n.git/\n*.zip\n*.patch\n*.log\n/mnt/\n".into(),
    ));
    files.push(deploy_file(
        "deploy/healthcheck.sh",
        "#!/usr/bin/env sh\nset -eu\ntest -d \"${CFD_CAMPAIGN_ROOT:-/work/campaigns}\"\n".into(),
    ));
    files.push(deploy_file(
        "deploy/worker-entrypoint.sh",
        "#!/usr/bin/env sh\nset -eu\n: \"${CFD_CAMPAIGN_ROOT:=/work/campaigns}\"\nmkdir -p \"$CFD_CAMPAIGN_ROOT\"\nexec /bin/sh -lc \"${CFD_WORKER_COMMAND:-cfd-solve particle-campaign-local-runner}\"\n".into(),
    ));
    files.push(deploy_file(
        "deploy/env.example",
        format!(
            "CFD_IMAGE={}\nCFD_CAMPAIGN_DIR={}\nCFD_WORKER_REPLICAS={}\n",
            c.image_name, c.host_campaign_dir, c.worker_replicas
        ),
    ));
    files.push(deploy_file(
        "deploy/servers.example.tsv",
        format!(
            "name\thost\tuser\tslots\truntime\tremote_root\tdocker_image\ttags\ngpu-a\t10.0.0.11\tcfd\t2\tdocker\t/srv/cfd/campaigns\t{}\tgpu,linux\ncpu-a\t10.0.0.12\tcfd\t4\tnative\t/srv/cfd/campaigns\t\tcpu,linux\n",
            c.image_name
        ),
    ));
    files.push(deploy_file(
        "deploy/deploy.sh",
        format!(
            "#!/usr/bin/env sh\nset -eu\ndocker compose -p {p} build\ndocker compose -p {p} up -d --scale {}={}\n",
            c.service_name,
            c.worker_replicas,
            p = c.compose_project_name
        ),
    ));
    files.push(deploy_file(
        "deploy/README.md",
        format!(
            "# Docker deployment\n\nThis directory contains a deterministic Docker/Compose deployment scaffold for cfd_solvers campaigns.\n\n1. Copy or generate campaign folders under `{}`.\n2. Build the image with `docker compose build`.\n3. Start workers with `docker compose up -d --scale {}={}`.\n4. Use `deploy/servers.example.tsv` as the shape for multi-server placement.\n\nThe generated files are an execution scaffold, not a security boundary. Review images, mounts and network exposure before production use.\n",
            c.host_campaign_dir, c.service_name, c.worker_replicas
        ),
    ));
    Ok(files)
}

pub fn write_docker_deploy_system(
    root: &Path,
    config: &DockerDeployConfig,
    overwrite: bool,
) -> Result<DockerDeployResult> {
    require(
        !root.as_os_str().is_empty(),
        "docker deploy writer requires root",
    )?;
    let files = generate_docker_deploy_files(config)?;
    for f in &files {
        ensure_relative(&f.relative_path)?;
    }
    create_dirs(root)?;
    let mut result = DockerDeployResult {
        root: root.to_path_buf(),
        files: files.clone(),
        files_written: 0,
    };
    for f in &files {
        let target = root.join(&f.relative_path);
        if target.exists() && !overwrite {
            return Err(SolverError::Io(format!(
                "docker deploy file already exists: {}",
                target.display()
            )));
        }
        if let Some(parent) = target.parent() {
            create_dirs(parent)?;
        }
        write_file(&target, &f.contents, "failed to write docker deploy file")?;
        result.files_written += 1;
        if target.extension().is_some_and(|e| e == "sh") {
            add_owner_exec(&target)?;
        }
    }
    Ok(result)
}
