//! Value-by-value comparison with the C++ upstream `fem/adaptivity` unit. Goldens come from
//! `parity/oracle_fem_adaptivity.cpp`.
use cfd_solvers::fem::*;
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-fem-adaptivity.json")).unwrap()
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
    let floor = 1e-6 * cpp.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= rel * c.abs().max(floor) + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn b(v: &[bool]) -> Vec<f64> {
    v.iter().map(|x| f64::from(u8::from(*x))).collect()
}

fn dump(p: &str, m: &Mesh2D) {
    let n: Vec<f64> = m.nodes.iter().flat_map(|x| [x.x, x.y]).collect();
    let t: Vec<f64> = m
        .triangles
        .iter()
        .flat_map(|x| x.node.map(|i| i as f64))
        .collect();
    let e: Vec<f64> = m
        .boundary_edges
        .iter()
        .flat_map(|x| [x.node[0] as f64, x.node[1] as f64, f64::from(x.patch)])
        .collect();
    close(&format!("{p}_nodes"), &n, 1e-15, 0.0);
    close(&format!("{p}_tris"), &t, 0.0, 0.0);
    close(&format!("{p}_bnd"), &e, 0.0, 0.0);
}

#[test]
fn adaptivity_matches_upstream() {
    let mesh = make_rectangle_tri_mesh(4, 3, 2.0, 1.5).unwrap();
    let u_of = |p: Node2| p.x * p.x + 0.3 * p.x * p.y + 0.1 * p.y * p.y * p.y;
    let u: Vec<f64> = mesh.nodes.iter().map(|p| u_of(*p)).collect();
    let src = |p: Node2| 1.0 + p.x + 0.2 * p.y * p.y;
    let est = estimate_poisson_error_tri3(&mesh, &u, &src, 1.7).unwrap();
    close("eta", &est.element_indicator, 1e-12, 0.0);
    close("eta_global", &[est.global_indicator], 1e-12, 0.0);
    let mark = mark_dorfler(&est.element_indicator, 0.5).unwrap();
    close("mark", &b(&mark), 0.0, 0.0);
    close(
        "mark_custom",
        &b(&mark_dorfler(&[0.1, 3.0, 0.5, 2.0, 0.0], 0.9).unwrap()),
        0.0,
        0.0,
    );
    close(
        "mark_zero",
        &b(&mark_dorfler(&[0.0, 0.0, 0.0], 0.5).unwrap()),
        0.0,
        0.0,
    );
    let r1 = refine_tri3_longest_edges(&mesh, &mark).unwrap();
    dump("r1", &r1);
    let all: Vec<bool> = (0..r1.element_count()).map(|i| i % 3 == 0).collect();
    let r2 = refine_tri3_longest_edges(&r1, &all).unwrap();
    dump("r2", &r2);
    let u2: Vec<f64> = r2.nodes.iter().map(|p| u_of(*p)).collect();
    close(
        "eta2_global",
        &[estimate_poisson_error_tri3(&r2, &u2, &src, 1.0)
            .unwrap()
            .global_indicator],
        1e-12,
        0.0,
    );
}
