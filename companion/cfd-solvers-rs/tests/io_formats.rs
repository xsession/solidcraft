//! Round-trip and robustness checks for checkpoints, OpenFOAM I/O and the XDMF contract.
use cfd_solvers::{core::Vec3, fem, fvm, io::*, SolverError};

#[test]
fn checkpoint_round_trip_and_corruption() {
    let dir = std::env::temp_dir().join(format!("cfd-io-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("state.chk");
    write_checkpoint(
        &path,
        "run",
        2.5,
        7,
        &[("rho", &[1.0, 2.0]), ("T", &[300.0])],
    )
    .unwrap();
    let c = read_checkpoint(&path).unwrap();
    assert_eq!(
        (c.case_name.as_str(), c.time, c.step, c.fields.len()),
        ("run", 2.5, 7, 2)
    );
    assert_eq!(c.fields[0].values, vec![1.0, 2.0]);
    let bytes = std::fs::read(&path).unwrap();
    for cut in [0, 7, 20, bytes.len() - 1] {
        assert!(decode_checkpoint(&bytes[..cut]).is_err());
    }
    let mut bad = bytes.clone();
    bad[0] = b'X';
    assert!(decode_checkpoint(&bad).is_err());
    // A huge declared field length is rejected before allocating.
    let mut huge = encode_checkpoint("a", 0.0, 0, &[("f", &[])]).unwrap();
    let n = huge.len();
    huge[n - 8..].copy_from_slice(&(u64::MAX / 16).to_le_bytes());
    assert!(decode_checkpoint(&huge).is_err());
    assert!(encode_checkpoint("a/b", 0.0, 0, &[]).is_err());
    assert!(encode_checkpoint("a", -1.0, 0, &[]).is_err());
    assert!(encode_checkpoint("a", 0.0, 0, &[("x", &[f64::NAN])]).is_err());
    assert!(encode_checkpoint("a", 0.0, 0, &[("x", &[]), ("x", &[])]).is_err());
    assert!(read_checkpoint(dir.join("missing")).is_err());
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn openfoam_fields_round_trip_and_mesh_is_closed() {
    let values = [1.25, -3.5, 1e-3];
    let text =
        openfoam_vol_scalar_field_text("T", &values, "[0 0 0 1 0 0 0]", "boundaryField {}\n");
    assert_eq!(parse_openfoam_vol_scalar_field(&text, 3).unwrap(), values);
    assert!(parse_openfoam_vol_scalar_field(&text, 4).is_err());
    assert!(parse_openfoam_vol_scalar_field(
        "internalField nonuniform List<vector> 1 ((0 0 0));",
        1
    )
    .is_err());
    let v =
        openfoam_vol_vector_field_text("U", &[Vec3::new(1.0, 0.0, 0.0)], "[0 1 -1 0 0 0 0]", "");
    assert!(v.contains("(1 0 0)\n"));
    let mesh = read_openfoam_polymesh("parity/openfoam-case").unwrap();
    assert_eq!(mesh.cells().len(), 2);
    // Gauss: the outward face areas of every cell sum to zero.
    for cell in 0..2 {
        let mut s = Vec3::new(0.0, 0.0, 0.0);
        for f in mesh.faces() {
            if f.owner == cell || f.neighbour == Some(cell) {
                let a = f.outward_area(cell);
                s = Vec3::new(s.x + a.x, s.y + a.y, s.z + a.z);
            }
        }
        assert!(s.x.abs() + s.y.abs() + s.z.abs() < 1e-12);
    }
    assert!((mesh.cells()[0].volume - 1.0).abs() < 1e-14);
    assert!(parse_openfoam_polymesh("1((0 0 0))", "1(3(0 0 0))", "1(0)", "0()", "0()").is_err());
    assert!(parse_openfoam_polymesh("no list", "", "", "", "").is_err());
    assert!(matches!(
        read_openfoam_polymesh("/nonexistent"),
        Err(SolverError::Io(_))
    ));
}

#[test]
fn hdf5_seams_and_xdmf_contract() {
    assert!(!fvm_hdf5_checkpoint_available());
    let mesh = fvm::make_cartesian_hexa_mesh(2, 1, 1, 1.0, 1.0, 1.0).unwrap();
    assert!(matches!(
        write_fvm_hdf5_checkpoint("/tmp/x.h5", &mesh, &[("T", &[1.0, 2.0])], 0.0, 0),
        Err(SolverError::Unsupported(_))
    ));
    assert!(matches!(
        write_fvm_hdf5_checkpoint("/tmp/x.h5", &mesh, &[("bad name", &[1.0, 2.0])], 0.0, 0),
        Err(SolverError::InvalidInput(_))
    ));
    assert!(read_fvm_hdf5_checkpoint("/tmp/x.h5").is_err());
    let m = fem::make_rectangle_tri_mesh(2, 2, 1.0, 1.0).unwrap();
    assert!(make_xdmf_mesh_contract(&m, &[], "anything goes").is_ok());
    assert!(make_xdmf_mesh_contract(&m, &[1.0], "T").is_err());
    assert!(make_xdmf_mesh_contract(&m, &vec![1.0; m.node_count()], "bad/name").is_err());
    assert!(!valid_xdmf_identifier(&"x".repeat(129)));
}
