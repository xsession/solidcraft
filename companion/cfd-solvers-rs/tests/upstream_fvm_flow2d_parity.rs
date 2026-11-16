//! Value-by-value comparison with the C++ upstream matrix-free CG, `projection2d`,
//! `incompressible2d`, `thin_film` and `compressible_vof` units.
//! Goldens come from `parity/oracle_fvm_flow2d.cpp`.
use cfd_solvers::{core::*, fvm::*};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-fvm-flow2d.json")).unwrap()
}

/// `|rust - cpp| <= rel * max(|cpp|, 1e-6 * max|cpp|) + abs` entry by entry.
#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let n = goldens()[name].as_array().map_or(0, |a| a.len());
    assert_eq!(rust.len(), n, "{name}: length");
    close_prefix(name, rust, rel, abs);
}

/// Like [`close`], but compares only the first `rust.len()` golden entries.
#[track_caller]
fn close_prefix(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let mut cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert!(rust.len() <= cpp.len(), "{name}: length");
    cpp.truncate(rust.len());
    let floor = 1e-6 * cpp.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= rel * c.abs().max(floor) + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

#[test]
fn matrix_free_cg_matches_upstream() {
    let rhs = [1.0, 2.0, -1.0, 0.5, 3.0, -2.0];
    let mut x = [0.0; 6];
    let r = matrix_free_conjugate_gradient(
        &rhs,
        &mut x,
        |a, b| {
            for i in 0..a.len() {
                b[i] = 4.0 * a[i]
                    - if i > 0 { a[i - 1] } else { 0.0 }
                    - if i + 1 < a.len() { a[i + 1] } else { 0.0 };
            }
        },
        100,
        1e-12,
    )
    .unwrap();
    let mut v = x.to_vec();
    v.extend([r.iterations as f64, if r.converged { 1.0 } else { 0.0 }]);
    close("cg", &v, 1e-12, 0.0);
}

#[test]
fn projection_matches_upstream() {
    let pc = Projection2DConfig {
        nx: 12,
        ny: 10,
        lx: 1.0,
        ly: 0.8,
        dt: 0.01,
        density: 2.0,
        ..Default::default()
    };
    let mut proj = Projection2D::new(pc).unwrap();
    proj.initialize_divergent(0.05);
    let div0 = proj.divergence_l2();
    proj.project().unwrap();
    // The projected divergent field is zero to round-off: compare absolutely.
    close("proj_u", proj.u_faces(), 0.0, 1e-13);
    close("proj_v", proj.v_faces(), 0.0, 1e-13);
    close("proj_p", proj.pressure(), 1e-10, 1e-12);
    let m = [
        div0,
        proj.divergence_l2(),
        proj.kinetic_energy(),
        proj.pressure_result().iterations as f64,
    ];
    close_prefix("proj_metrics", &m[..1], 1e-13, 0.0);
    assert!(m[1] < 1e-12 && m[2] < 1e-28);
    let g = goldens();
    assert_eq!(m[3], g["proj_metrics"][3].as_f64().unwrap());
    let mut tg = Projection2D::new(pc).unwrap();
    tg.initialize_taylor_green(0.05);
    let ke0 = tg.kinetic_energy();
    tg.project().unwrap();
    close_prefix("proj_tg", &[ke0, tg.kinetic_energy()], 1e-11, 0.0);
    assert!(tg.divergence_l2() < 1e-12);
}

#[test]
fn incompressible_taylor_green_matches_upstream() {
    let mut ns = Incompressible2D::new(Incompressible2DConfig {
        nx: 16,
        ny: 16,
        kinematic_viscosity: 0.02,
        dt: 1e-3,
        ..Default::default()
    })
    .unwrap();
    ns.initialize_taylor_green(0.1);
    ns.run(20).unwrap();
    close("ns_u", ns.u_faces(), 1e-9, 1e-13);
    close("ns_v", ns.v_faces(), 1e-9, 1e-13);
    close("ns_p", ns.pressure(), 1e-7, 1e-11);
    close_prefix("ns_metrics", &[ns.time(), ns.kinetic_energy()], 1e-10, 0.0);
    assert!(ns.divergence_l2() < 1e-12);
    let g = goldens();
    let err = g["ns_metrics"][3].as_f64().unwrap();
    assert!((ns.taylor_green_velocity_error(0.1) - err).abs() < 1e-9 * err.max(1e-12) + 1e-14);
}

#[test]
fn thin_film_matches_upstream() {
    let mesh = make_sheared_cartesian_hexa_mesh(8, 5, 1, 0.08, 0.05, 0.01, 0.0).unwrap();
    let cfg = ThinFilmConfig {
        dt: 2e-3,
        surface_tension: 0.07,
        viscosity: 0.05,
        body_acceleration: Vec3::new(0.0, -9.81, 0.0),
        ..Default::default()
    };
    let mut film = ThinFilmTransport::new(mesh.clone(), cfg).unwrap();
    film.initialize_profile(|p| {
        1e-3 + 4e-4 * (-2000.0 * ((p.x - 0.04).powi(2) + (p.y - 0.025).powi(2))).exp()
    });
    for _ in 0..5 {
        film.step().unwrap();
    }
    close("film_h", film.thickness(), 1e-11, 0.0);
    close("film_p", film.pressure(), 1e-8, 0.0);
    close("film_flux", film.face_flux(), 1e-7, 0.0);
    close("film_inventory", &[film.inventory()], 1e-12, 0.0);
    let bp: Vec<f64> = mesh.faces().iter().map(|f| 10.0 * f.center.x).collect();
    film.set_boundary_pressure(bp).unwrap();
    film.step().unwrap();
    close("film_open_h", film.thickness(), 1e-11, 0.0);
}

#[test]
fn compressible_vof_matches_upstream() {
    let mesh = make_sheared_cartesian_hexa_mesh(8, 5, 1, 0.08, 0.05, 0.01, 0.0).unwrap();
    let mut cv = CompressibleVofTransport::new(
        mesh.clone(),
        CompressibleVofConfig {
            dt: 1e-3,
            liquid_bulk_modulus: 2.0e7,
            ..Default::default()
        },
    )
    .unwrap();
    cv.initialize_profile(|p| if p.x < 0.03 { 1.0 } else { 0.0 }, 1.2e5)
        .unwrap();
    let flux = mesh
        .faces()
        .iter()
        .map(|f| {
            if f.neighbour.is_none() {
                0.0
            } else {
                f.area.x * 0.4 * (1.0 + 5.0 * f.center.x) + f.area.y * 0.1
            }
        })
        .collect();
    cv.set_face_flux(flux).unwrap();
    for _ in 0..6 {
        cv.step().unwrap();
    }
    close("cvof_alpha", cv.liquid_fraction(), 1e-10, 0.0);
    close("cvof_density", cv.density(), 1e-10, 0.0);
    close("cvof_pressure", cv.pressure(), 1e-9, 0.0);
    let phase: Vec<f64> = cv
        .liquid_density()
        .iter()
        .zip(cv.gas_density())
        .flat_map(|(l, g)| [*l, *g])
        .collect();
    close("cvof_phase_density", &phase, 1e-10, 0.0);
    close(
        "cvof_totals",
        &[
            cv.liquid_volume(),
            cv.total_mass(),
            cv.minimum_pressure(),
            cv.maximum_pressure(),
            cv.time(),
        ],
        1e-10,
        0.0,
    );
}
