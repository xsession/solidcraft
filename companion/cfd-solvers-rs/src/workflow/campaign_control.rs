//! File-based campaign control directives (stop / extend / checkpoint / flush), Slurm queue
//! parsing and registry refresh from solver outputs. Ported from upstream
//! `cfd/workflow/campaign_control` (MIT, cfd_solvers contributors).
use super::campaign::{
    discover_campaign_outputs, load_campaign_registry, max_iteration, parse_performance_metrics,
    parse_residual_history, read_text_full, save_campaign_registry, solver_has_capability,
    sort_residuals, summarize_campaign_registry, validate_solver_adapter_descriptor,
    CampaignCaseStatus, CampaignRegistryEntry, CampaignSummary, SolverAdapterDescriptor,
    SolverCapability,
};
use super::text::{
    create_dirs, getline_split, lower, stod_full, trim_ascii, upper, write_file, ws_tokens,
};
use crate::{require, Result, SolverError};
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum CampaignControlAction {
    #[default]
    Stop,
    Extend,
    Checkpoint,
    Flush,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CampaignControlRequest {
    pub case_id: String,
    pub action: CampaignControlAction,
    pub value: String,
    pub reason: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CampaignControlOptions {
    pub control_dirname: String,
    pub overwrite_latest: bool,
}

impl Default for CampaignControlOptions {
    fn default() -> Self {
        Self {
            control_dirname: ".cfd_control".into(),
            overwrite_latest: true,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CampaignControlDirective {
    pub case_id: String,
    pub action: CampaignControlAction,
    pub value: String,
    pub reason: String,
    pub path: PathBuf,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SchedulerJobState {
    #[default]
    Unknown,
    Queued,
    Running,
    Done,
    Failed,
    Cancelled,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SchedulerJobRecord {
    pub job_id: String,
    pub case_id: String,
    pub state: SchedulerJobState,
    pub elapsed_seconds: f64,
    pub raw_state: String,
}

fn escape_value(v: &str) -> String {
    let mut o = String::with_capacity(v.len());
    for c in v.chars() {
        match c {
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            _ => o.push(c),
        }
    }
    o
}

fn unescape_value(v: &str) -> String {
    let mut o = String::with_capacity(v.len());
    let mut it = v.chars();
    while let Some(c) = it.next() {
        if c == '\\' {
            if let Some(n) = it.next() {
                o.push(match n {
                    'n' => '\n',
                    'r' => '\r',
                    other => other,
                });
                continue;
            }
        }
        o.push(c);
    }
    o
}

fn parse_key_values(text: &str) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for line in getline_split(text, b'\n') {
        if let Some(pos) = line.find('=') {
            out.insert(
                trim_ascii(&line[..pos]).to_string(),
                unescape_value(trim_ascii(&line[pos + 1..])),
            );
        }
    }
    out
}

/// Directive file body: `case_id=`, `action=`, `value=`, `reason=` lines (backslash-escaped).
pub fn campaign_control_payload(request: &CampaignControlRequest) -> String {
    format!(
        "case_id={}\naction={}\nvalue={}\nreason={}\n",
        escape_value(&request.case_id),
        campaign_control_action_name(request.action),
        escape_value(&request.value),
        escape_value(&request.reason)
    )
}

pub fn campaign_control_action_name(action: CampaignControlAction) -> &'static str {
    match action {
        CampaignControlAction::Stop => "stop",
        CampaignControlAction::Extend => "extend",
        CampaignControlAction::Checkpoint => "checkpoint",
        CampaignControlAction::Flush => "flush",
    }
}

pub fn campaign_control_action_from_name(name: &str) -> Result<CampaignControlAction> {
    Ok(match lower(name).as_str() {
        "stop" => CampaignControlAction::Stop,
        "extend" => CampaignControlAction::Extend,
        "checkpoint" => CampaignControlAction::Checkpoint,
        "flush" => CampaignControlAction::Flush,
        _ => {
            return Err(SolverError::InvalidInput(
                "unknown campaign control action".into(),
            ))
        }
    })
}

/// The adapter must advertise the `control` capability and list the action by name.
pub fn validate_campaign_control_request(
    adapter: &SolverAdapterDescriptor,
    request: &CampaignControlRequest,
) -> Result<()> {
    validate_solver_adapter_descriptor(adapter)?;
    require(
        !request.case_id.is_empty(),
        "campaign control requires case_id",
    )?;
    let action = campaign_control_action_name(request.action);
    require(
        solver_has_capability(adapter, SolverCapability::Control)
            && adapter.control_actions.iter().any(|a| a == action),
        format!("solver adapter does not support control action: {action}"),
    )
}

/// Write `<root>/<case>/<control dir>/<action>.directive` (and `latest.directive`).
pub fn write_campaign_control_directive(
    root: &Path,
    adapter: &SolverAdapterDescriptor,
    request: &CampaignControlRequest,
    options: &CampaignControlOptions,
) -> Result<CampaignControlDirective> {
    require(
        !root.as_os_str().is_empty() && !options.control_dirname.is_empty(),
        "control directive requires root and control directory",
    )?;
    validate_campaign_control_request(adapter, request)?;
    let case_dir = root.join(&request.case_id);
    if !case_dir.is_dir() {
        return Err(SolverError::Io(format!(
            "campaign control case directory is missing: {}",
            request.case_id
        )));
    }
    let control_dir = case_dir.join(&options.control_dirname);
    create_dirs(&control_dir)?;
    let action = campaign_control_action_name(request.action);
    let payload = campaign_control_payload(request);
    let target = control_dir.join(format!("{action}.directive"));
    if target.exists() && !options.overwrite_latest {
        return Err(SolverError::Io(format!(
            "control directive already exists: {}",
            target.display()
        )));
    }
    write_file(&target, &payload, "failed to write control directive")?;
    if options.overwrite_latest {
        // Upstream ignores a failure to update the convenience copy.
        let _ = std::fs::write(control_dir.join("latest.directive"), &payload);
    }
    Ok(CampaignControlDirective {
        case_id: request.case_id.clone(),
        action: request.action,
        value: request.value.clone(),
        reason: request.reason.clone(),
        path: target,
    })
}

/// Every `*.directive` file with `case_id` and `action` keys, sorted by (case id, path).
pub fn discover_campaign_control_directives(
    case_directory: &Path,
    options: &CampaignControlOptions,
) -> Result<Vec<CampaignControlDirective>> {
    let mut out = Vec::new();
    if case_directory.as_os_str().is_empty() || options.control_dirname.is_empty() {
        return Ok(out);
    }
    let control_dir = case_directory.join(&options.control_dirname);
    let Ok(entries) = std::fs::read_dir(&control_dir) else {
        return Ok(out);
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_file() || path.extension().is_none_or(|e| e != "directive") {
            continue;
        }
        let kv = parse_key_values(&read_text_full(&path));
        let (Some(case_id), Some(action)) = (kv.get("case_id"), kv.get("action")) else {
            continue;
        };
        out.push(CampaignControlDirective {
            case_id: case_id.clone(),
            action: campaign_control_action_from_name(action)?,
            value: kv.get("value").cloned().unwrap_or_default(),
            reason: kv.get("reason").cloned().unwrap_or_default(),
            path,
        });
    }
    out.sort_by(|a, b| a.case_id.cmp(&b.case_id).then_with(|| a.path.cmp(&b.path)));
    Ok(out)
}

pub fn scheduler_job_state_name(state: SchedulerJobState) -> &'static str {
    match state {
        SchedulerJobState::Unknown => "unknown",
        SchedulerJobState::Queued => "queued",
        SchedulerJobState::Running => "running",
        SchedulerJobState::Done => "done",
        SchedulerJobState::Failed => "failed",
        SchedulerJobState::Cancelled => "cancelled",
    }
}

/// Slurm short/long state codes and generic names (case-insensitive); anything else is unknown.
pub fn scheduler_job_state_from_name(name: &str) -> SchedulerJobState {
    match upper(name).as_str() {
        "PD" | "PENDING" | "QUEUED" | "CONFIGURING" => SchedulerJobState::Queued,
        "R" | "RUNNING" | "COMPLETING" => SchedulerJobState::Running,
        "CD" | "COMPLETED" | "DONE" | "SUCCESS" => SchedulerJobState::Done,
        "F" | "FAILED" | "TIMEOUT" | "NODE_FAIL" | "OUT_OF_MEMORY" | "ERROR" => {
            SchedulerJobState::Failed
        }
        "CA" | "CANCELLED" | "CANCELED" | "STOPPED" => SchedulerJobState::Cancelled,
        _ => SchedulerJobState::Unknown,
    }
}

/// Rows `job_id case_id state [elapsed_seconds]` separated by whitespace, `,`, `;` or tabs;
/// `JOBID` headers and `#` comments are skipped.
pub fn parse_slurm_queue_table(text: &str) -> Vec<SchedulerJobRecord> {
    let mut out = Vec::new();
    for line in getline_split(text, b'\n') {
        let mapped: String = line
            .chars()
            .map(|c| {
                if matches!(c, ',' | ';' | '\t') {
                    ' '
                } else {
                    c
                }
            })
            .collect();
        let row = ws_tokens(&mapped);
        let Some(first) = row.first() else {
            continue;
        };
        if upper(first) == "JOBID" || first.starts_with('#') || row.len() < 3 {
            continue;
        }
        out.push(SchedulerJobRecord {
            job_id: row[0].clone(),
            case_id: row[1].clone(),
            state: scheduler_job_state_from_name(&row[2]),
            elapsed_seconds: row.get(3).and_then(|t| stod_full(t)).unwrap_or(0.0),
            raw_state: row[2].clone(),
        });
    }
    out
}

pub fn apply_scheduler_records_to_registry(
    registry: &mut Vec<CampaignRegistryEntry>,
    records: &[SchedulerJobRecord],
) -> CampaignSummary {
    for record in records {
        if record.case_id.is_empty() {
            continue;
        }
        let idx = match registry.iter().position(|e| e.case_id == record.case_id) {
            Some(i) => i,
            None => {
                registry.push(CampaignRegistryEntry::new(
                    record.case_id.clone(),
                    CampaignCaseStatus::Pending,
                    0.0,
                    0,
                    format!("scheduler discovered job={}", record.job_id),
                ));
                registry.len() - 1
            }
        };
        let e = &mut registry[idx];
        match record.state {
            SchedulerJobState::Queued | SchedulerJobState::Running => {
                e.status = CampaignCaseStatus::Running
            }
            SchedulerJobState::Done => e.status = CampaignCaseStatus::Done,
            SchedulerJobState::Failed => e.status = CampaignCaseStatus::Failed,
            SchedulerJobState::Cancelled => e.status = CampaignCaseStatus::Stopped,
            SchedulerJobState::Unknown => {}
        }
        e.message = format!("scheduler job={} state={}", record.job_id, record.raw_state);
    }
    summarize_campaign_registry(registry)
}

/// Re-derive status, iterations and objective of every case that has log files.
pub fn refresh_campaign_status_from_outputs(
    root: &Path,
    registry: &mut [CampaignRegistryEntry],
) -> Result<CampaignSummary> {
    require(
        !root.as_os_str().is_empty(),
        "campaign output refresh requires root",
    )?;
    for entry in registry.iter_mut() {
        let case_dir = root.join(&entry.case_id);
        if !case_dir.exists() {
            continue;
        }
        let discovery = discover_campaign_outputs(&case_dir);
        let mut combined = String::new();
        for p in &discovery.log_files {
            combined.push_str(&read_text_full(p));
            combined.push('\n');
        }
        if combined.is_empty() {
            continue;
        }
        let mut residuals = parse_residual_history(&combined);
        let mut perf = parse_performance_metrics(&combined);
        for p in &discovery.residual_files {
            residuals.extend(parse_residual_history(&read_text_full(p)));
        }
        for p in &discovery.performance_files {
            perf.extend(parse_performance_metrics(&read_text_full(p)));
        }
        sort_residuals(&mut residuals);
        let outcome = super::campaign::detect_solver_outcome_from_log(&combined);
        if outcome != CampaignCaseStatus::Running {
            entry.status = outcome;
        } else if entry.status == CampaignCaseStatus::Pending {
            entry.status = CampaignCaseStatus::Running;
        }
        entry.iterations = max_iteration(&residuals);
        if let Some(last) = residuals.last() {
            entry.objective = last.residual;
        }
        entry.message = format!(
            "refreshed residuals={} perf={}",
            residuals.len(),
            perf.len()
        );
    }
    Ok(summarize_campaign_registry(registry))
}

pub fn refresh_campaign_registry_from_outputs(
    root: &Path,
    registry_filename: &str,
) -> Result<CampaignSummary> {
    require(
        !root.as_os_str().is_empty() && !registry_filename.is_empty(),
        "campaign registry refresh requires root and registry filename",
    )?;
    let path = root.join(registry_filename);
    let mut registry = load_campaign_registry(&path)?;
    let summary = refresh_campaign_status_from_outputs(root, &mut registry)?;
    save_campaign_registry(&path, &registry)?;
    Ok(summary)
}
