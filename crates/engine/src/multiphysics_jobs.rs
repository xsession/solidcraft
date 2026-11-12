//! Session-owned numerical work; saved authoring and immutable retained evidence are separate.
use crate::{EngineError, Result, Session};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use solidcraft_multiphysics::{RunEvidence, Study, case_hash, preflight, validate_evidence};
use std::sync::Arc;

const MAX_JOBS: usize = 8;
const MAX_LIVE: usize = 2;

#[derive(Clone)]
pub struct FrozenRun {
    pub id: String,
    pub study: Study,
    pub case_hash: String,
}
pub enum Outcome {
    Completed(RunEvidence),
    Failed(String),
    Cancelled,
}
pub trait Host: Send + Sync {
    fn start(&self, request: FrozenRun) -> std::result::Result<Box<dyn Job>, String>;
}
pub trait Job: Send + Sync {
    fn poll(&mut self) -> Option<Outcome>;
    fn cancel(&mut self);
}

#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct JobSnapshot {
    pub id: String,
    pub study_id: String,
    pub study_revision: u64,
    pub case_hash: String,
    pub state: String,
    pub message: String,
    pub stale: bool,
    pub evidence: Option<Arc<RunEvidence>>,
}
impl JobSnapshot {
    pub fn active(&self) -> bool {
        matches!(self.state.as_str(), "running" | "cancelling")
    }
}
struct Entry {
    snapshot: JobSnapshot,
    study: Study,
    session_uid: u64,
    document_id: String,
    handle: Option<Box<dyn Job>>,
}
impl Drop for Entry {
    fn drop(&mut self) {
        if let Some(handle) = self.handle.as_mut() {
            handle.cancel();
        }
    }
}
#[derive(Default)]
pub(crate) struct Jobs {
    host: Option<Arc<dyn Host>>,
    entries: Vec<Entry>,
    next_run: u64,
    pub(crate) next_study: u64,
}
impl Jobs {
    pub(crate) fn scratch(&self) -> Self {
        Self { host: self.host.clone(), ..Self::default() }
    }
}

impl Session {
    pub fn set_multiphysics_host(&mut self, host: Option<Arc<dyn Host>>) {
        self.cancel_multiphysics_jobs();
        self.multiphysics_jobs.host = host;
    }
    pub fn multiphysics_host_available(&self) -> bool {
        self.multiphysics_jobs.host.is_some()
    }
    pub fn cancel_multiphysics_jobs(&mut self) {
        for entry in &mut self.multiphysics_jobs.entries {
            if let Some(job) = entry.handle.as_mut() {
                job.cancel();
                entry.snapshot.state = "cancelling".into();
                entry.snapshot.message = "Waiting for the numerical worker to stop".into();
            }
        }
    }
    pub fn poll_multiphysics_jobs(&mut self) {
        for entry in &mut self.multiphysics_jobs.entries {
            let Some(outcome) = entry.handle.as_mut().and_then(|job| job.poll()) else { continue };
            let cancelling = entry.snapshot.state == "cancelling";
            entry.handle = None;
            match outcome {
                Outcome::Cancelled => {
                    entry.snapshot.state = "cancelled".into();
                    entry.snapshot.message = "The numerical worker stopped".into();
                }
                Outcome::Failed(message) => {
                    entry.snapshot.state = if cancelling { "cancelled" } else { "failed" }.into();
                    entry.snapshot.message = message.chars().take(4096).collect();
                }
                Outcome::Completed(_) if cancelling => {
                    entry.snapshot.state = "cancelled".into();
                    entry.snapshot.message = "Completion arrived after cancellation; result discarded".into();
                }
                Outcome::Completed(evidence) => {
                    let checked = if evidence.id != entry.snapshot.id {
                        Err("Worker returned another job identity".into())
                    } else {
                        validate_evidence(&entry.study, &evidence).map_err(|e| e.to_string())
                    };
                    match checked {
                        Err(message) => {
                            entry.snapshot.state = "failed".into();
                            entry.snapshot.message = message;
                        }
                        Ok(()) => {
                            entry.snapshot.state = "succeeded".into();
                            entry.snapshot.message = "Execution completed; inspect convergence and model assumptions".into();
                            entry.snapshot.evidence = Some(Arc::new(evidence));
                        }
                    }
                }
            }
        }
    }
    pub fn multiphysics_jobs(&self) -> Vec<JobSnapshot> {
        self.multiphysics_jobs
            .entries
            .iter()
            .map(|entry| {
                let mut snapshot = entry.snapshot.clone();
                snapshot.stale = entry.session_uid != self.uid
                    || entry.document_id != self.document_id
                    || self.doc.multiphysics.study(&entry.study.id).is_none_or(|study| {
                        study.revision != entry.study.revision || case_hash(&study.case).ok().as_ref() != Some(&snapshot.case_hash)
                    });
                snapshot
            })
            .collect()
    }
    pub(crate) fn start_multiphysics_run(&mut self, id: &str) -> Result<JobSnapshot> {
        if self.is_copy {
            return Err(EngineError::Other("Numerical work starts only in a live document session".into()));
        }
        self.poll_multiphysics_jobs();
        let mut study = self.doc.multiphysics.study(id).cloned().ok_or_else(|| EngineError::Other("Select an existing Multiphysics study".into()))?;
        let report = preflight(&study);
        if !report.errors.is_empty() {
            return Err(EngineError::Other(report.errors.join("; ")));
        }
        let host = self.multiphysics_jobs.host.clone().ok_or_else(|| {
            EngineError::Other(
                "Configure SOLIDCRAFT_CFD_SOLVE with the cfd-solvers-rs 0.2.0 companion; this host currently supports preparation only".into(),
            )
        })?;
        if self.multiphysics_jobs.entries.iter().filter(|entry| entry.snapshot.active()).count() >= MAX_LIVE
            || self.multiphysics_jobs.entries.iter().any(|entry| entry.snapshot.study_id == id && entry.snapshot.active())
        {
            return Err(EngineError::Other("Wait for the running study or cancel it; at most two numerical workers run at once".into()));
        }
        if self.multiphysics_jobs.entries.len() >= MAX_JOBS {
            if let Some(index) = self.multiphysics_jobs.entries.iter().position(|entry| !entry.snapshot.active()) {
                self.multiphysics_jobs.entries.remove(index);
            } else {
                return Err(EngineError::Other("Numerical job history is full".into()));
            }
        }
        let next = self.multiphysics_jobs.next_run.checked_add(1).ok_or_else(|| EngineError::Other("Numerical run identities exhausted".into()))?;
        let job_id = format!("mp-{}-{next}", self.uid);
        let hash = case_hash(&study.case).map_err(|e| EngineError::Other(e.to_string()))?;
        study.runs.clear();
        let handle = host.start(FrozenRun { id: job_id.clone(), study: study.clone(), case_hash: hash.clone() }).map_err(EngineError::Other)?;
        self.multiphysics_jobs.next_run = next;
        let snapshot = JobSnapshot {
            id: job_id,
            study_id: id.into(),
            study_revision: study.revision,
            case_hash: hash,
            state: "running".into(),
            message: "Numerical worker started".into(),
            stale: false,
            evidence: None,
        };
        self.multiphysics_jobs.entries.push(Entry {
            snapshot: snapshot.clone(),
            study,
            session_uid: self.uid,
            document_id: self.document_id.clone(),
            handle: Some(handle),
        });
        Ok(snapshot)
    }
    pub(crate) fn cancel_multiphysics_run(&mut self, id: &str) -> Result<JobSnapshot> {
        self.poll_multiphysics_jobs();
        let entry = self
            .multiphysics_jobs
            .entries
            .iter_mut()
            .find(|entry| entry.snapshot.id == id)
            .ok_or_else(|| EngineError::Other("Numerical job not found in this document".into()))?;
        if let Some(job) = entry.handle.as_mut() {
            job.cancel();
            entry.snapshot.state = "cancelling".into();
            entry.snapshot.message = "Waiting for the numerical worker to stop".into();
        }
        Ok(entry.snapshot.clone())
    }
}
#[cfg(test)]
mod tests;
