//! Byte comparison with the C++ upstream python/c_api (`cfd_solvers_version`,
//! `cfd_write_xdmf_tri3`). Goldens come from `parity/oracle_python.cpp`.
#![allow(clippy::approx_constant)]
use cfd_solvers::python::*;
use serde_json::Value;

fn g() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-python.json")).unwrap()
}

const XY: [f64; 10] = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.1, 0.7, 1.0 / 3.0, -2.5e-7];
const TRI: [usize; 9] = [0, 1, 2, 0, 2, 3, 1, 4, 0];

#[test]
fn xdmf_writer_matches_upstream_bytes() {
    let g = g();
    let dir = std::env::temp_dir().join(format!("cfd_capi_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let s = [1.0, -0.5, 2.0 / 3.0, 1e20, 3.14159];
    let a = dir.join("a.xmf");
    write_xdmf_tri3(&a, &XY, 5, &TRI, 3, Some(&s), Some("temperature")).unwrap();
    assert_eq!(
        std::fs::read_to_string(&a).unwrap(),
        g["with_scalar"].as_str().unwrap()
    );
    let b = dir.join("b.xmf");
    write_xdmf_tri3(&b, &XY, 5, &TRI, 3, None, None).unwrap();
    assert_eq!(
        std::fs::read_to_string(&b).unwrap(),
        g["without_scalar"].as_str().unwrap()
    );
    let e = write_xdmf_tri3(&dir.join("c.xmf"), &XY, 0, &TRI, 3, None, None).unwrap_err();
    assert!(e.to_string().contains(g["error"].as_str().unwrap()), "{e}");
    assert_eq!(UPSTREAM_VERSION, g["version"].as_str().unwrap());
    std::fs::remove_dir_all(&dir).ok();
}
