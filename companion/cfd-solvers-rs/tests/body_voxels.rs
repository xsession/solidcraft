//! Surface-body voxelization against the exact generalized winding number (solid angles).
//! Cell centres lying exactly on an edge shared by two triangles must count once: this prism's
//! side face is split along the diagonal y + z = 0.01, which passes through 50 cell centres.

use cfd_solvers::{case::*, core::Vec3};

fn prism_and_tetrahedron() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut v = vec![[0.03, 0.0, 0.0], [0.05, 0.0, 0.0], [0.03, 0.01, 0.0], [0.03, 0.0, 0.01], [0.05, 0.0, 0.01], [0.03, 0.01, 0.01]];
    let mut t = vec![[0, 2, 1], [3, 4, 5], [0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4], [2, 0, 3], [2, 3, 5]];
    let o = v.len();
    v.extend([[0.06, 0.002, 0.001], [0.078, 0.004, 0.002], [0.065, 0.018, 0.003], [0.068, 0.007, 0.0095]]);
    t.extend([[o, o + 2, o + 1], [o, o + 1, o + 3], [o + 1, o + 2, o + 3], [o, o + 3, o + 2]]);
    (v, t)
}

/// Generalized winding number of a closed outward surface at `p` (1 inside, 0 outside).
fn winding(v: &[[f64; 3]], t: &[[usize; 3]], p: [f64; 3]) -> f64 {
    let mut w = 0.0;
    for tri in t {
        let [a, b, c] = tri.map(|i| [v[i][0] - p[0], v[i][1] - p[1], v[i][2] - p[2]]);
        let n = |u: [f64; 3]| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let dot = |u: [f64; 3], q: [f64; 3]| u[0] * q[0] + u[1] * q[1] + u[2] * q[2];
        let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
        let den = n(a) * n(b) * n(c) + dot(a, b) * n(c) + dot(b, c) * n(a) + dot(c, a) * n(b);
        w += 2.0 * det.atan2(den);
    }
    w / (4.0 * std::f64::consts::PI)
}

#[test]
fn voxels_match_the_exact_winding_number_with_centres_on_shared_edges() {
    let (v, t) = prism_and_tetrahedron();
    let body: SurfaceBody = serde_json::from_value(serde_json::json!({"name": "w", "vertices": v, "triangles": t})).unwrap();
    let [nx, ny, nz] = [40, 20, 10];
    let spec = PolyMeshSpec { dimensions: [nx, ny, nz], lengths: Vec3::new(0.1, 0.02, 0.01), xy_shear: 0.0, coordinates: None };
    let inside = body.voxelize(&spec).unwrap();
    let mut on_diagonal = 0;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let p = [(i as f64 + 0.5) * 0.1 / nx as f64, (j as f64 + 0.5) * 0.02 / ny as f64, (k as f64 + 0.5) * 0.01 / nz as f64];
                let exact = winding(&v, &t, p).abs() > 0.5;
                assert_eq!(inside[(k * ny + j) * nx + i], exact, "cell ({i}, {j}, {k})");
                if (p[1] + p[2] - 0.01).abs() < 1e-15 && (0.0..=0.01).contains(&p[1]) {
                    on_diagonal += 1;
                }
            }
        }
    }
    assert!(on_diagonal >= 50, "{on_diagonal}");
    assert_eq!(inside.iter().filter(|c| **c).count(), 540);
}
