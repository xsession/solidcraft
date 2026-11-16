//! Value-by-value comparison with the C++ upstream io/{checkpoint, openfoam_io, xdmf_contract}.
//! Goldens come from `parity/oracle_io.cpp`; the OpenFOAM case is `parity/openfoam-case`.
#![allow(clippy::approx_constant)]
use cfd_solvers::{core::Vec3, fem, io::*};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-io.json")).unwrap()
}

#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= rel * c.abs() + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn text(name: &str) -> String {
    goldens()[name].as_str().unwrap().to_string()
}

#[test]
fn checkpoint_openfoam_and_xdmf_match_upstream() {
    let bytes = encode_checkpoint(
        "cavity",
        0.125,
        42,
        &[("U", &[1.5, -2.25, 1e-300]), ("p", &[3.0])],
    )
    .unwrap();
    close(
        "checkpoint_bytes",
        &bytes.iter().map(|b| f64::from(*b)).collect::<Vec<_>>(),
        0.0,
        0.0,
    );
    let mesh = parse_openfoam_polymesh(
        include_str!("../parity/openfoam-case/points"),
        include_str!("../parity/openfoam-case/faces"),
        include_str!("../parity/openfoam-case/owner"),
        include_str!("../parity/openfoam-case/neighbour"),
        include_str!("../parity/openfoam-case/boundary"),
    )
    .unwrap();
    let mut mv: Vec<f64> = mesh
        .cells()
        .iter()
        .flat_map(|c| [c.center.x, c.center.y, c.center.z, c.volume])
        .collect();
    for f in mesh.faces() {
        mv.extend([
            f.owner as f64,
            f.neighbour.map_or(-1.0, |n| n as f64),
            f.patch.map_or(-1.0, |p| p as f64),
            f.center.x,
            f.center.y,
            f.center.z,
            f.area.x,
            f.area.y,
            f.area.z,
        ]);
    }
    close("openfoam_mesh", &mv, 1e-15, 1e-16);
    assert_eq!(
        mesh.patches()
            .iter()
            .map(|p| p.name.as_str())
            .collect::<Vec<_>>(),
        ["left", "right", "walls", "frontBack"]
    );
    let mut sv =
        parse_openfoam_vol_scalar_field(include_str!("../parity/openfoam-case/T"), 2).unwrap();
    sv.extend(
        parse_openfoam_vol_scalar_field(include_str!("../parity/openfoam-case/p"), 3).unwrap(),
    );
    close("openfoam_scalar", &sv, 0.0, 0.0);
    assert_eq!(
        openfoam_vol_scalar_field_text(
            "T",
            &[300.123456789, -1e-7, 2.5e12],
            "[0 0 0 0 0 0 0]",
            "boundaryField {}\n"
        ),
        text("openfoam_scalar_text")
    );
    assert_eq!(
        openfoam_vol_vector_field_text(
            "U",
            &[
                Vec3::new(1.0, 2.0 / 3.0, -0.5),
                Vec3::new(1e20, 0.0, 3.14159265)
            ],
            "[0 1 -1 0 0 0 0]",
            "boundaryField {}\n"
        ),
        text("openfoam_vector_text")
    );
    let m2 = fem::make_rectangle_tri_mesh(3, 2, 1.0, 1.0).unwrap();
    let c = make_xdmf_mesh_contract(&m2, &vec![1.0; m2.node_count()], "temp.K").unwrap();
    close(
        "xdmf",
        &[
            c.node_count as f64,
            c.triangle_count as f64,
            c.field_count as f64,
            f64::from(u8::from(valid_xdmf_identifier("bad name"))),
        ],
        0.0,
        0.0,
    );
}
