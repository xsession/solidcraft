use cfd_solvers::{
    case::{Case, RunManifest, SOLVERS},
    core::Vec3,
    io::*,
};
use std::{
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "cfd-rust-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
#[test]
fn every_embedded_example_executes_and_exports() {
    let folder = Scratch::new();
    for name in [
        "poisson1d",
        "poisson2d",
        "poisson3d",
        "scalar-mixed",
        "heat2d",
        "lbm-taylor-green",
        "maxwell1d",
        "sod-weno",
        "diffusion2d",
        "rf-tem",
        "optics-lens",
        "chemistry-batch",
        "acoustic-kspace2d",
        "poly-scalar",
        "flow-channel",
        "flow-cavity-simple",
        "flow-cavity-pimple",
        "solid-heat",
        "thermal-radiation",
        "rans-sa",
        "rans-kepsilon",
        "rans-sst-des",
        "rans-reynolds-stress",
    ] {
        let out = folder.0.join(name);
        let result = Command::new(env!("CARGO_BIN_EXE_cfd-solve"))
            .args(["example", name, "--out"])
            .arg(&out)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let printed: RunManifest = serde_json::from_slice(&result.stdout).unwrap();
        let manifest: RunManifest =
            serde_json::from_slice(&std::fs::read(out.join("manifest.json")).unwrap()).unwrap();
        assert!(manifest.execution_completed);
        assert!(!manifest.complete_clone);
        assert_eq!(printed.solver, manifest.solver);
        assert!(SOLVERS.contains(&manifest.solver.as_str()));
        for f in manifest.files {
            assert!(out.join(f).is_file());
        }
    }
}
#[test]
fn invalid_cases_fail_and_previous_runs_are_retained() {
    let folder = Scratch::new();
    let invalid = folder.0.join("invalid.json");
    std::fs::write(&invalid, r#"{"solver":"not_a_solver","dummy":true}"#).unwrap();
    let out = folder.0.join("out");
    let failed = Command::new(env!("CARGO_BIN_EXE_cfd-solve"))
        .arg("run")
        .arg(&invalid)
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    assert_eq!(failed.status.code(), Some(2));
    assert!(!out.join("manifest.json").exists());
    assert!(serde_json::from_str::<Case>(r#"{"solver":"poisson1d","diffusvity":1.0}"#).is_err());
    let c: Case = serde_json::from_str(r#"{"solver":"poisson1d"}"#).unwrap();
    c.run(&out).unwrap();
    let old = std::fs::read(out.join("manifest.json")).unwrap();
    assert!(c.run(&out).is_err());
    assert_eq!(old, std::fs::read(out.join("manifest.json")).unwrap());
}
#[test]
fn vtu_writer_validates_associations_and_escapes_names() {
    let folder = Scratch::new();
    let m = VtuMesh {
        points: vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)],
        cells: vec![Cell {
            cell_type: CellType::Line,
            nodes: vec![0, 1],
        }],
    };
    let path = folder.0.join("mesh.vtu");
    assert!(write_vtu(&path, &m, &[DataArray::scalar("bad", vec![1.0])], &[]).is_err());
    assert!(!path.exists());
    write_vtu(
        &path,
        &m,
        &[DataArray::scalar("T & \"wall\"", vec![2.0, 3.0])],
        &[DataArray::scalar("material", vec![1.0])],
    )
    .unwrap();
    let text = std::fs::read_to_string(path).unwrap();
    assert!(text.contains("T &amp; &quot;wall&quot;"));
    assert!(text.contains("NumberOfPoints=\"2\" NumberOfCells=\"1\""));
    assert!(write_csv(folder.0.join("bad.csv"), &[("x", &[0.0]), ("y", &[])]).is_err());
}
#[cfg(feature = "parallel")]
#[test]
fn rayon_lbm_matches_serial_bit_for_bit() {
    use cfd_solvers::lbm::*;
    let mut serial = D2Q9Solver::new(D2Q9Config {
        nx: 96,
        ny: 96,
        tau: 0.8,
        ..Default::default()
    })
    .unwrap();
    serial.initialize_taylor_green(0.02).unwrap();
    let mut parallel = serial.clone();
    serial.step_serial(4).unwrap();
    parallel.step(4).unwrap();
    assert_eq!(serial.populations(), parallel.populations());
    assert_eq!(serial.density(), parallel.density());
    assert_eq!(serial.velocity_x(), parallel.velocity_x());
}
