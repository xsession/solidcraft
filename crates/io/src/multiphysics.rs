//! Portable numerical cases and original ParaView-compatible result handoffs.
use crate::{IoError, Result, safe_write};
use serde_json::json;
use solidcraft_multiphysics::{RunEvidence, Study, validate_record};
use std::path::Path;

fn invalid(error: impl ToString) -> IoError {
    IoError::Invalid(error.to_string())
}
pub fn write_case(path: &Path, study: &Study) -> Result<()> {
    crate::vfs::check(&path.to_string_lossy()).map_err(invalid)?;
    solidcraft_multiphysics::case_hash(&study.case).map_err(invalid)?;
    if path.as_os_str().is_empty() {
        return Err(invalid("Choose a new case filename"));
    }
    safe_write::create_atomic(path, &serde_json::to_vec_pretty(&study.case).map_err(invalid)?).map_err(invalid)
}
/// A module of a component's setup as a new file of its own (`solidcraft-multiphysics-module/1`).
pub fn write_module(path: &Path, module: &solidcraft_multiphysics::modules::Module) -> Result<()> {
    crate::vfs::check(&path.to_string_lossy()).map_err(invalid)?;
    if path.as_os_str().is_empty() {
        return Err(invalid("Choose a new module filename"));
    }
    let text = solidcraft_multiphysics::modules::to_file(module).map_err(invalid)?;
    safe_write::create_atomic(path, text.as_bytes()).map_err(invalid)
}
/// A module file [`write_module`] wrote, read and checked.
pub fn read_module(path: &Path) -> Result<solidcraft_multiphysics::modules::Module> {
    crate::vfs::check(&path.to_string_lossy()).map_err(invalid)?;
    let size = std::fs::metadata(path).map_err(invalid)?.len();
    if size > solidcraft_multiphysics::modules::MAX_FILE_BYTES as u64 {
        return Err(invalid("A module file is at most 4 MiB"));
    }
    let text = std::fs::read_to_string(path).map_err(invalid)?;
    solidcraft_multiphysics::modules::from_file(&text).map_err(invalid)
}
/// A batch report into a new directory: `report.csv`, `report.html` and its pictures (plain
/// file names), written atomically; the directory is removed again if any file fails.
pub fn write_report(directory: &Path, csv: &str, html: &str, pictures: &[(String, Vec<u8>)]) -> Result<Vec<String>> {
    crate::vfs::check(&directory.to_string_lossy()).map_err(invalid)?;
    if cfg!(target_arch = "wasm32") {
        return Err(invalid("Report directories require a native host"));
    }
    if directory.as_os_str().is_empty() {
        return Err(invalid("Choose a new output directory"));
    }
    if pictures.iter().any(|(name, _)| !name.ends_with(".png") || name.contains(['/', '\\']) || name.starts_with('.') || name.len() > 64) {
        return Err(invalid("A report picture needs a plain .png file name"));
    }
    std::fs::create_dir(directory).map_err(invalid)?;
    let write = || -> Result<Vec<String>> {
        let mut files = Vec::new();
        for (name, bytes) in [("report.csv", csv.as_bytes()), ("report.html", html.as_bytes())] {
            safe_write::create_atomic(&directory.join(name), bytes).map_err(invalid)?;
            files.push(name.to_string());
        }
        for (name, bytes) in pictures {
            safe_write::create_atomic(&directory.join(name), bytes).map_err(invalid)?;
            files.push(name.clone());
        }
        Ok(files)
    };
    let result = write();
    if result.is_err() {
        let _ = std::fs::remove_dir_all(directory);
    }
    result
}

/// The named text files of a run folder that exist there (each at most 64 MiB), for reading
/// results the design did not keep.
pub fn read_run_texts(directory: &Path, names: &[&str]) -> Result<Vec<(String, String)>> {
    crate::vfs::check(&directory.to_string_lossy()).map_err(invalid)?;
    if !directory.join("manifest.json").is_file() {
        return Err(invalid(format!("{} holds no run (no manifest.json)", directory.display())));
    }
    let mut out = Vec::new();
    for name in names {
        let path = directory.join(name);
        let Ok(meta) = std::fs::metadata(&path) else { continue };
        if !meta.is_file() || meta.len() > 64 << 20 {
            continue;
        }
        out.push(((*name).to_string(), std::fs::read_to_string(&path).map_err(invalid)?));
    }
    Ok(out)
}

/// Text files (plain names) into a new directory, written atomically; the directory is removed
/// again if any file fails.
pub fn write_texts(directory: &Path, files: &[(&str, &str)]) -> Result<Vec<String>> {
    crate::vfs::check(&directory.to_string_lossy()).map_err(invalid)?;
    if cfg!(target_arch = "wasm32") {
        return Err(invalid("Output directories require a native host"));
    }
    if directory.as_os_str().is_empty() {
        return Err(invalid("Choose a new output directory"));
    }
    if files.iter().any(|(name, _)| name.contains(['/', '\\']) || name.starts_with('.') || name.is_empty() || name.len() > 64) {
        return Err(invalid("An output file needs a plain file name"));
    }
    std::fs::create_dir(directory).map_err(invalid)?;
    let write = || -> Result<Vec<String>> {
        let mut written = Vec::new();
        for (name, text) in files {
            safe_write::create_atomic(&directory.join(name), text.as_bytes()).map_err(invalid)?;
            written.push((*name).to_string());
        }
        Ok(written)
    };
    let result = write();
    if result.is_err() {
        let _ = std::fs::remove_dir_all(directory);
    }
    result
}

pub fn write_results(directory: &Path, study: &Study, run: &RunEvidence) -> Result<Vec<String>> {
    validate_record(run).map_err(invalid)?;
    crate::vfs::check(&directory.to_string_lossy()).map_err(invalid)?;
    if cfg!(target_arch = "wasm32") {
        return Err(invalid("Portable result directories require a native host"));
    }
    if directory.as_os_str().is_empty() {
        return Err(invalid("Choose a new output directory"));
    }
    for reserved in ["case.json", "solidcraft-study.json", "README.txt"] {
        if run.artifacts.iter().any(|a| a.name == reserved) {
            return Err(invalid("Artifact collides with a reserved handoff filename"));
        }
    }
    std::fs::create_dir(directory).map_err(invalid)?;
    let write = || -> Result<Vec<String>> {
        let mut files = Vec::new();
        for artifact in &run.artifacts {
            safe_write::create_atomic(&directory.join(&artifact.name), artifact.text.as_bytes()).map_err(invalid)?;
            files.push(artifact.name.clone());
        }
        for (name, value) in [
            ("case.json", run.input_case.clone()),
            (
                "solidcraft-study.json",
                json!({"schema":"solidcraft-multiphysics-handoff/1","study_id":study.id,
                "current_study_metadata":{"name":study.name,"profile":study.profile,"question":study.question,"revision":study.revision},
                "run_id":run.id,"captured_study_revision":run.study_revision,"case_sha256":run.case_hash,
                "representation":"explicit companion case; independent of CAD geometry","physical_validation":"not_assessed",
                "field_semantics":"producer names; no inferred SI units","time_axis":"steady SIMPLE series use iteration, not seconds"}),
            ),
        ] {
            safe_write::create_atomic(&directory.join(name), &serde_json::to_vec_pretty(&value).map_err(invalid)?).map_err(invalid)?;
            files.push(name.into());
        }
        let guide = "SolidCraft Multiphysics handoff\n\nOpen the .pvd collection or .vtu files in ParaView. CSV, JSON and Touchstone\n.s2p artifacts remain in their original producer format. PointData/CellData\nassociations and values are unchanged. Save a ParaView state after opening;\nno .pvsm was fabricated.\n\nRead manifest.json and the original case.json before comparisons. Flow pressure\nis kinematic pressure (m^2/s^2), not Pa. SIMPLE collection axes are iterations.\nMaxwell E/H may use separate staggered grids. Generic potential has no inferred\nphysical unit. Completion, linear convergence, steady convergence and physical\nvalidation are distinct. Current-study metadata may differ from the captured\nrun's original case; captured_study_revision identifies that historical input.\n";
        safe_write::create_atomic(&directory.join("README.txt"), guide.as_bytes()).map_err(invalid)?;
        files.push("README.txt".into());
        Ok(files)
    };
    let result = write();
    if result.is_err() {
        let _ = std::fs::remove_dir_all(directory);
    }
    result
}
