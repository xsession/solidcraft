//! Checks for the safe counterpart of upstream `python/c_api.h`.
use cfd_solvers::python::*;

#[test]
fn short_buffers_and_bad_indices_are_errors() {
    let dir = std::env::temp_dir().join(format!("cfd_capi_api_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("x.xmf");
    let xy = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
    assert!(write_xdmf_tri3(&p, &xy[..4], 3, &[0, 1, 2], 1, None, None).is_err());
    assert!(write_xdmf_tri3(&p, &xy, 3, &[0, 1], 1, None, None).is_err());
    assert!(write_xdmf_tri3(&p, &xy, 3, &[0, 1, 2], 1, Some(&[1.0]), None).is_err());
    assert!(write_xdmf_tri3(&p, &xy, 3, &[0, 1, 7], 1, None, None).is_err());
    assert!(write_xdmf_tri3(
        &dir.join("missing").join("x.xmf"),
        &xy,
        3,
        &[0, 1, 2],
        1,
        None,
        None
    )
    .is_err());
    write_xdmf_tri3(&p, &xy, 3, &[0, 1, 2], 1, Some(&[1.0, 2.0, 3.0]), None).unwrap();
    let text = std::fs::read_to_string(&p).unwrap();
    assert!(text.contains("Attribute Name=\"field\"") && text.contains("1 2 3 \n"));
    assert_eq!(cfd_solvers_version(), env!("CARGO_PKG_VERSION"));
    std::fs::remove_dir_all(&dir).ok();
}
