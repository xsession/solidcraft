use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::Read;
use std::path::{Component, Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, mpsc};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use serde_json::Value;
use solidcraft_engine::multiphysics::{Artifact, RunEvidence, validate_artifact_name, validate_evidence};
use solidcraft_engine::multiphysics_jobs::{FrozenRun, Host, Job, Outcome};

const MAX_MANIFEST_BYTES: usize = 1 << 20;
const MAX_ARTIFACT_BYTES: usize = 16 << 20;
const MAX_ARTIFACTS: usize = 128;
const MAX_STDERR_BYTES: usize = 16 << 10;
const PORT_VERSION: &str = "0.2.0";
const UPSTREAM_REVISION: &str = "7cceaaaddafaf2515e01a86406e13d6037fdbf9a";
static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(1);

/// Host configuration is application state, never a field in an imported study.
#[derive(Clone)]
pub struct ProcessHost {
    executable: PathBuf,
    timeout: Duration,
}

impl ProcessHost {
    /// Use an explicitly chosen companion. No subprocess is started by this constructor.
    pub fn new(executable: PathBuf) -> Self {
        Self { executable, timeout: Duration::from_secs(300) }
    }
}

/// Discover an explicit environment configuration or an installed `cfd-solve` on PATH.
/// Discovery checks files only; a slow or broken executable cannot block application startup.
pub fn configured_host() -> Option<Arc<dyn Host>> {
    if let Some(value) = std::env::var_os("SOLIDCRAFT_CFD_SOLVE") {
        if value.is_empty() {
            return None;
        }
        // Keep an explicit invalid path configured: starting a run reports the actual OS error.
        return Some(Arc::new(ProcessHost::new(PathBuf::from(value))));
    }
    let path = std::env::var_os("PATH")?;
    let name = if cfg!(windows) { "cfd-solve.exe" } else { "cfd-solve" };
    std::env::split_paths(&path).map(|folder| folder.join(name)).find(|candidate| candidate.is_file()).map(|candidate| {
        let host: Arc<dyn Host> = Arc::new(ProcessHost::new(candidate));
        host
    })
}

struct ProcessJob {
    receiver: Mutex<mpsc::Receiver<Outcome>>,
    cancelled: Arc<AtomicBool>,
    finished: bool,
}

impl Host for ProcessHost {
    fn start(&self, request: FrozenRun) -> Result<Box<dyn Job>, String> {
        let host = self.clone();
        let cancelled = Arc::new(AtomicBool::new(false));
        let worker_cancelled = cancelled.clone();
        let (sender, receiver) = mpsc::channel();
        std::thread::Builder::new()
            .name("solidcraft-multiphysics".into())
            .spawn(move || {
                let outcome = host.run(request, &worker_cancelled);
                let _ = sender.send(outcome);
            })
            .map_err(|error| format!("cannot start solver worker: {error}"))?;
        Ok(Box::new(ProcessJob { receiver: Mutex::new(receiver), cancelled, finished: false }))
    }
}

impl Job for ProcessJob {
    fn poll(&mut self) -> Option<Outcome> {
        if self.finished {
            return None;
        }
        let outcome = match self.receiver.lock() {
            Ok(receiver) => match receiver.try_recv() {
                Ok(outcome) => Some(outcome),
                Err(mpsc::TryRecvError::Empty) => None,
                Err(mpsc::TryRecvError::Disconnected) => Some(Outcome::Failed("solver worker ended without a terminal result".into())),
            },
            Err(_) => Some(Outcome::Failed("solver result channel is unavailable".into())),
        };
        self.finished = outcome.is_some();
        outcome
    }

    fn cancel(&mut self) {
        self.cancelled.store(true, Ordering::Release);
    }
}

impl Drop for ProcessJob {
    fn drop(&mut self) {
        self.cancelled.store(true, Ordering::Release);
    }
}

impl ProcessHost {
    fn run(&self, request: FrozenRun, cancelled: &AtomicBool) -> Outcome {
        match self.run_inner(request, cancelled) {
            Ok(outcome) => outcome,
            Err(error) => Outcome::Failed(error),
        }
    }

    fn run_inner(&self, request: FrozenRun, cancelled: &AtomicBool) -> Result<Outcome, String> {
        if cancelled.load(Ordering::Acquire) {
            return Ok(Outcome::Cancelled);
        }
        let directory = WorkDirectory::create()?;
        let case_path = directory.0.join("case.json");
        let output = directory.0.join("results");
        let case_bytes = serde_json::to_vec(&request.study.case).map_err(|error| format!("cannot encode frozen case: {error}"))?;
        if case_bytes.len() > MAX_MANIFEST_BYTES {
            return Err("frozen case exceeds the 1 MiB companion input limit".into());
        }
        solidcraft_engine::io::safe_write::write_atomic_no_follow(&case_path, &case_bytes)
            .map_err(|error| format!("cannot write frozen case: {error}"))?;
        let mut child = Command::new(&self.executable)
            .arg("run")
            .arg(&case_path)
            .arg("--out")
            .arg(&output)
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|error| format!("cannot start cfd-solve ({}): {error}", self.executable.display()))?;
        // Drain stderr while the solver runs, retaining only a bounded prefix. The pipe cannot
        // fill and deadlock a numerical run, and a verbose failure cannot allocate unbounded RAM.
        let stderr = child.stderr.take();
        let diagnostic = Arc::new(Mutex::new(Vec::new()));
        if let Some(mut stderr) = stderr {
            let diagnostic = diagnostic.clone();
            if let Err(error) = std::thread::Builder::new().name("solidcraft-solver-log".into()).spawn(move || {
                let mut buffer = [0_u8; 4096];
                loop {
                    let count = match stderr.read(&mut buffer) {
                        Ok(0) | Err(_) => break,
                        Ok(count) => count,
                    };
                    if let Ok(mut retained) = diagnostic.lock() {
                        let keep = count.min(MAX_STDERR_BYTES.saturating_sub(retained.len()));
                        if let Some(bytes) = buffer.get(..keep) {
                            retained.extend_from_slice(bytes);
                        }
                    }
                }
            }) {
                terminate(&mut child);
                return Err(format!("cannot drain companion diagnostics: {error}"));
            }
        }
        let started = Instant::now();
        let status = loop {
            if cancelled.load(Ordering::Acquire) {
                terminate(&mut child);
                return Ok(Outcome::Cancelled);
            }
            if started.elapsed() > self.timeout {
                terminate(&mut child);
                return Err("cfd-solve exceeded the 300-second runtime limit".into());
            }
            if let Err(error) = output_budget(&output) {
                terminate(&mut child);
                return Err(error);
            }
            match child.try_wait() {
                Ok(Some(status)) => break status,
                Ok(None) => std::thread::sleep(Duration::from_millis(20)),
                Err(error) => {
                    terminate(&mut child);
                    return Err(format!("cannot read companion status: {error}"));
                }
            }
        };
        if cancelled.load(Ordering::Acquire) {
            return Ok(Outcome::Cancelled);
        }
        if !status.success() {
            let log = diagnostic.lock().ok().map(|bytes| String::from_utf8_lossy(&bytes).trim().to_string()).unwrap_or_default();
            return Err(if log.is_empty() { format!("cfd-solve failed with {status}") } else { format!("cfd-solve failed with {status}: {log}") });
        }
        let evidence = read_evidence(&output, &request)?;
        if cancelled.load(Ordering::Acquire) {
            return Ok(Outcome::Cancelled);
        }
        Ok(Outcome::Completed(evidence))
    }
}

fn terminate(child: &mut Child) {
    let _ = child.kill();
    let _ = child.wait();
}

/// The admitted format is flat. Check growth during execution as well as after exit, so a
/// faulty companion cannot fill the disk for the whole runtime before artifact admission.
fn output_budget(output: &Path) -> Result<(), String> {
    let metadata = match fs::symlink_metadata(output) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(format!("cannot inspect solver output: {error}")),
    };
    if !metadata.is_dir() || metadata.file_type().is_symlink() {
        return Err("companion output must be its own ordinary directory".into());
    }
    let entries = fs::read_dir(output).map_err(|error| format!("cannot inspect solver output: {error}"))?;
    let mut total = 0_u64;
    for (index, entry) in entries.enumerate() {
        if index >= MAX_ARTIFACTS {
            return Err("companion produced more than 128 files".into());
        }
        let entry = entry.map_err(|error| format!("cannot inspect companion artifact: {error}"))?;
        let metadata = match entry.metadata() {
            Ok(metadata) => metadata,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => continue,
            Err(error) => return Err(format!("cannot inspect companion artifact: {error}")),
        };
        // read_evidence rejects links which leave the result directory. Never walk a subtree.
        if metadata.is_dir() {
            return Err("companion output contains an unsupported artifact directory".into());
        }
        total = total.saturating_add(metadata.len());
        if total > MAX_ARTIFACT_BYTES as u64 {
            return Err("companion result exceeds the 16 MiB output limit".into());
        }
    }
    Ok(())
}

struct WorkDirectory(PathBuf);

impl WorkDirectory {
    fn create() -> Result<Self, String> {
        let stamp = SystemTime::now().duration_since(UNIX_EPOCH).map_err(|error| format!("cannot name solver task: {error}"))?.as_nanos();
        for _ in 0..16 {
            let serial = NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!("solidcraft-multiphysics-{}-{stamp}-{serial}", std::process::id()));
            let mut builder = fs::DirBuilder::new();
            builder.recursive(false);
            #[cfg(unix)]
            {
                use std::os::unix::fs::DirBuilderExt as _;
                builder.mode(0o700);
            }
            match builder.create(&path) {
                Ok(()) => return Ok(Self(path)),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(format!("cannot create private solver task: {error}")),
            }
        }
        Err("cannot allocate a unique solver task directory".into())
    }
}

impl Drop for WorkDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn ordinary_relative(name: &str) -> bool {
    !name.is_empty()
        && name.len() <= 240
        && !name.contains(['\\', ':', '\0'])
        && name.split('/').all(|part| !part.is_empty() && part != "." && part != "..")
        && Path::new(name).components().all(|component| matches!(component, Component::Normal(_)))
}

fn read_file(root: &Path, name: &str, limit: usize) -> Result<String, String> {
    if !ordinary_relative(name) {
        return Err(format!("companion artifact path is not a portable relative path: {name}"));
    }
    let path = root.join(name).canonicalize().map_err(|error| format!("cannot resolve artifact {name}: {error}"))?;
    if !path.starts_with(root) {
        return Err(format!("companion artifact leaves its task directory: {name}"));
    }
    let metadata = fs::metadata(&path).map_err(|error| format!("cannot inspect artifact {name}: {error}"))?;
    if !metadata.is_file() || metadata.len() > limit as u64 {
        return Err(format!("companion artifact is not a bounded regular file: {name}"));
    }
    let mut bytes = Vec::new();
    File::open(&path)
        .and_then(|file| file.take(limit as u64 + 1).read_to_end(&mut bytes))
        .map_err(|error| format!("cannot read artifact {name}: {error}"))?;
    if bytes.len() > limit {
        return Err(format!("companion artifact grew beyond its byte limit: {name}"));
    }
    String::from_utf8(bytes).map_err(|error| format!("companion artifact is not UTF-8 text ({name}): {error}"))
}

fn read_evidence(output: &Path, request: &FrozenRun) -> Result<RunEvidence, String> {
    let root = output.canonicalize().map_err(|error| format!("companion produced no results directory: {error}"))?;
    // The directory is freshly allocated by this host; a solver may not replace it with a link.
    if fs::symlink_metadata(output).map_err(|error| error.to_string())?.file_type().is_symlink() {
        return Err("companion results directory is a symbolic link".into());
    }
    let text = read_file(&root, "manifest.json", MAX_MANIFEST_BYTES)?;
    let manifest: Value = serde_json::from_str(&text).map_err(|error| format!("invalid companion manifest: {error}"))?;
    for (key, expected) in [
        ("schema", "cfd-solvers-rs/run/1"),
        ("rust_port_version", PORT_VERSION),
        ("upstream_revision", UPSTREAM_REVISION),
        ("backend", "native_rust_serial_cpu"),
    ] {
        if manifest.get(key).and_then(Value::as_str) != Some(expected) {
            return Err(format!("unsupported companion {key}: expected {expected}"));
        }
    }
    let names = manifest.get("files").and_then(Value::as_array).ok_or("companion manifest has no artifact list")?;
    if names.is_empty() || names.len() > MAX_ARTIFACTS {
        return Err("companion artifact count is outside 1..128".into());
    }
    let mut seen = BTreeSet::new();
    let mut artifacts = Vec::new();
    let mut remaining = MAX_ARTIFACT_BYTES;
    for value in names {
        let name = value.as_str().ok_or("companion artifact name is not a string")?;
        validate_artifact_name(name).map_err(|error| format!("invalid relative path for companion artifact: {error}"))?;
        if !seen.insert(name) {
            return Err("companion artifact names repeat".into());
        }
        let artifact_text = if name == "manifest.json" { text.clone() } else { read_file(&root, name, remaining.min(8 << 20))? };
        if artifact_text.len() > remaining {
            return Err("companion result exceeds 16 MiB".into());
        }
        remaining -= artifact_text.len();
        artifacts.push(Artifact { name: name.to_string(), text: artifact_text });
    }
    let evidence = RunEvidence {
        id: request.id.clone(),
        study_revision: request.study.revision,
        case_hash: request.case_hash.clone(),
        input_case: request.study.case.clone(),
        manifest,
        artifacts,
    };
    validate_evidence(&request.study, &evidence).map_err(|error| format!("companion evidence was refused: {error}"))?;
    Ok(evidence)
}

#[cfg(test)]
mod tests;
