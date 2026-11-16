//! Value-by-value (and, for the writers, byte-by-byte) comparison with the C++ upstream LBM CPU
//! units. Goldens come from `parity/oracle_lbm.cpp` built against the pinned upstream revision.
use cfd_solvers::{core::Vec3, io::*, lbm::*};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-lbm.json")).unwrap()
}

/// `|rust - cpp| <= rel * max(|cpp|, 1e-6 * max|cpp|) + abs` entry by entry.
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

#[track_caller]
fn same_text(name: &str, path: &std::path::Path) {
    let g = goldens();
    let expected = g[name].as_str().unwrap();
    let actual = std::fs::read_to_string(path).unwrap();
    if actual != expected {
        let at = actual
            .chars()
            .zip(expected.chars())
            .position(|(a, b)| a != b)
            .unwrap_or(actual.len().min(expected.len()));
        let lo = at.saturating_sub(40);
        panic!(
            "{name} differs at char {at}:\nrust: {:?}\nC++:  {:?}",
            &actual[lo..(at + 40).min(actual.len())],
            &expected[lo..(at + 40).min(expected.len())]
        );
    }
}

fn desc<D: LatticeDescriptor>() -> Vec<f64> {
    let mut o = vec![D::DIMENSIONS as f64, D::Q as f64];
    for i in 0..D::Q {
        o.extend([
            D::cx(i) as f64,
            D::cy(i) as f64,
            D::cz(i) as f64,
            f64::from(D::weight(i)),
            D::opposite(i) as f64,
        ]);
    }
    o
}

#[test]
fn descriptors_match_upstream() {
    close("d2q9_inplace", &desc::<D2Q9InPlaceDescriptor>(), 0.0, 0.0);
    close("d3q19", &desc::<D3Q19Descriptor>(), 0.0, 0.0);
    close("d3q27", &desc::<D3Q27Descriptor>(), 0.0, 0.0);
    assert!(
        descriptor_pairs_are_opposites::<D2Q9InPlaceDescriptor>()
            && descriptor_pairs_are_opposites::<D3Q19Descriptor>()
            && descriptor_pairs_are_opposites::<D3Q27Descriptor>()
    );
    assert!(
        descriptor_weights_sum_to_one::<D2Q9InPlaceDescriptor>()
            && descriptor_weights_sum_to_one::<D3Q19Descriptor>()
            && descriptor_weights_sum_to_one::<D3Q27Descriptor>()
    );
}

#[test]
fn diagnostics_match_upstream() {
    let mut solver = D2Q9Solver::new(D2Q9Config {
        nx: 24,
        ny: 16,
        tau: 0.8,
        acceleration_x: 1e-5,
        acceleration_y: 0.0,
    })
    .unwrap();
    for y in 6..10 {
        for x in 8..12 {
            solver.set_solid(x, y, true).unwrap();
        }
    }
    solver.initialize_uniform(1.0, 0.02, 0.0).unwrap();
    solver.step(30).unwrap();
    let ft = momentum_exchange_force(&solver, 10.0, 8.0);
    // Drag is compared relatively; lift and torque are float round-off around zero.
    let drag = goldens()["momentum_exchange"][0].as_f64().unwrap();
    assert!(
        (ft.force_x - drag).abs() <= 1e-5 * drag.abs(),
        "drag {} vs C++ {drag}",
        ft.force_x
    );
    assert!(ft.force_y.abs() < 1e-6 && ft.torque_z.abs() < 1e-6);
    let a = darcy_forchheimer_acceleration(0.03, -0.01, 0.05, 2e-3, 0.4).unwrap();
    close("darcy", &a, 1e-15, 0.0);
}

fn advanced<S>(
    mut s: S,
    tg: impl Fn(&mut S),
    step: impl Fn(&mut S),
    mass: impl Fn(&S) -> f64,
    ke: impl Fn(&S) -> f64,
) -> Vec<f64> {
    tg(&mut s);
    let mut o = vec![mass(&s), ke(&s)];
    step(&mut s);
    o.extend([mass(&s), ke(&s)]);
    o
}

#[test]
fn advanced_collisions_match_upstream() {
    let ac = AdvancedD2Q9Config {
        nx: 16,
        ny: 12,
        tau: 0.65,
    };
    let r = advanced(
        RegularizedD2Q9Solver::new(ac).unwrap(),
        |s| s.initialize_taylor_green(0.04),
        |s| s.step(25),
        |s| s.mass(),
        |s| s.kinetic_energy(),
    );
    close("regularized", &r, 1e-6, 0.0);
    let c = advanced(
        CumulantD2Q9Solver::new(ac).unwrap(),
        |s| s.initialize_taylor_green(0.04),
        |s| s.step(25),
        |s| s.mass(),
        |s| s.kinetic_energy(),
    );
    close("cumulant", &c, 1e-6, 0.0);
    let m = advanced(
        MrtD2Q9Solver::new(ac).unwrap(),
        |s| s.initialize_taylor_green(0.04),
        |s| s.step(25),
        |s| s.mass(),
        |s| s.kinetic_energy(),
    );
    close("mrt", &m, 1e-6, 0.0);
}

#[test]
fn passive_scalar_matches_upstream() {
    let mut ps = PassiveScalarD2Q5::new(PassiveScalarConfig {
        nx: 20,
        ny: 10,
        diffusivity: 0.02,
        ux: 0.05,
        uy: -0.02,
    })
    .unwrap();
    ps.initialize(|x, y| (-40.0 * ((x - 0.3).powi(2) + (y - 0.5).powi(2))).exp());
    let t0 = ps.total();
    ps.step(40);
    close("scalar", &ps.scalar(), 1e-12, 0.0);
    close("scalar_total", &[t0, ps.total()], 1e-13, 0.0);
}

fn tet() -> TriangleSurface {
    let (a, b, c, d) = (
        Vec3::new(0.1, 0.1, 0.1),
        Vec3::new(0.9, 0.15, 0.1),
        Vec3::new(0.2, 0.85, 0.15),
        Vec3::new(0.3, 0.3, 0.9),
    );
    let t = |a, b, c| Triangle3 { a, b, c };
    TriangleSurface {
        triangles: vec![t(a, c, b), t(a, b, d), t(b, c, d), t(c, a, d)],
    }
}

#[test]
fn voxelization_and_surface_readers_match_upstream() {
    let moved = transform_surface(
        &tet(),
        Vec3::new(0.05, -0.02, 0.0),
        Vec3::new(0.2, 0.3, 1.0),
        0.3,
    );
    let mv: Vec<f64> = moved
        .triangles
        .iter()
        .flat_map(|t| [t.a, t.b, t.c])
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();
    close("transformed", &mv, 1e-14, 1e-16);
    let vg = VoxelGrid3D {
        nx: 10,
        ny: 9,
        nz: 8,
        minimum: Vec3::new(0.0, 0.0, 0.0),
        maximum: Vec3::new(1.0, 1.0, 1.0),
    };
    let as_f = |m: Vec<u8>| m.into_iter().map(f64::from).collect::<Vec<_>>();
    close(
        "voxels",
        &as_f(voxelize_surface(&moved, &vg).unwrap()),
        0.0,
        0.0,
    );

    let dir = std::env::temp_dir().join(format!("cfd-lbm-parity-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut stl = String::from("solid t\n");
    for t in &tet().triangles {
        stl.push_str(" facet normal 0 0 0\n  outer loop\n");
        for p in [t.a, t.b, t.c] {
            stl.push_str(&format!("   vertex {} {} {}\n", p.x, p.y, p.z));
        }
        stl.push_str("  endloop\n endfacet\n");
    }
    stl.push_str("endsolid t\n");
    std::fs::write(dir.join("t.stl"), stl).unwrap();
    std::fs::write(dir.join("t.obj"), "# tet\nv 0.1 0.1 0.1\nv 0.9 0.15 0.1\nv 0.2 0.85 0.15\nv 0.3 0.3 0.9\nf 1 3 2\nf 1/1 2/2 4/4\nf 2//2 3//3 4//4\nf 3 1 4\n").unwrap();
    let mut both =
        as_f(voxelize_surface(&read_ascii_stl(dir.join("t.stl")).unwrap(), &vg).unwrap());
    both.extend(as_f(
        voxelize_surface(&read_obj_triangles(dir.join("t.obj")).unwrap(), &vg).unwrap(),
    ));
    close("stl_obj_voxels", &both, 0.0, 0.0);
    std::fs::write(dir.join("bad.obj"), "v 0 0 0\nf 1 2 3\n").unwrap();
    assert!(read_obj_triangles(dir.join("bad.obj")).is_err());
    std::fs::write(
        dir.join("bad.stl"),
        "solid x\nvertex 0 0 0\nvertex 1 0 0\nendsolid\n",
    )
    .unwrap();
    assert!(read_ascii_stl(dir.join("bad.stl")).is_err());
    let _ = std::fs::remove_dir_all(dir);
}

type Fields = (Vec<f64>, Vec<f32>, Vec<f32>, Vec<f32>);

fn free_surface_fields() -> Fields {
    let (fx, fy, fz) = (8usize, 6usize, 5usize);
    let n = fx * fy * fz;
    let (mut fill, mut ux, mut uy, mut uz) = (
        vec![0.0; n],
        vec![0.0f32; n],
        vec![0.0f32; n],
        vec![0.0f32; n],
    );
    for z in 0..fz {
        for y in 0..fy {
            for x in 0..fx {
                let i = (z * fy + y) * fx + x;
                let r = ((x as f64 - 3.5).powi(2)
                    + (y as f64 - 2.5).powi(2)
                    + (z as f64 - 2.0).powi(2))
                .sqrt();
                fill[i] = (2.2 - r).clamp(0.0, 1.0);
                ux[i] = 0.3 + 0.05 * y as f32;
                uy[i] = -0.1 * (x as f32).sin();
                uz[i] = 0.05;
            }
        }
    }
    (fill, ux, uy, uz)
}

fn flat(v: &[ParticleVector3]) -> Vec<f64> {
    v.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

#[test]
fn free_surface_matches_upstream() {
    let (fill, ux, uy, uz) = free_surface_fields();
    let mut fs = make_free_surface_field(8, 6, 5, &fill, 1e-8).unwrap();
    let v0 = free_surface_volume(&fs);
    for _ in 0..3 {
        advect_free_surface_vof(&mut fs, &ux, &uy, &uz, 0.8, 1e-8).unwrap();
    }
    close("fs_fill", &fs.fill, 1e-12, 1e-15);
    let phase: Vec<f64> = fs
        .phase
        .iter()
        .map(|p| match p {
            FreeSurfacePhase::Gas => 0.0,
            FreeSurfacePhase::Interface => 1.0,
            FreeSurfacePhase::Fluid => 2.0,
        })
        .collect();
    close("fs_phase", &phase, 0.0, 0.0);
    close("fs_volume", &[v0, free_surface_volume(&fs)], 1e-14, 0.0);
    close(
        "fs_normals",
        &flat(&free_surface_normals(&fs)),
        1e-10,
        1e-13,
    );
    close("fs_curvature", &free_surface_curvature(&fs), 1e-10, 1e-13);
    close(
        "fs_csf",
        &flat(&free_surface_surface_tension_acceleration(&fs, 0.07, 1.3).unwrap()),
        1e-10,
        1e-15,
    );
}

#[test]
fn immersed_boundary_particles_match_upstream() {
    let (_, ux, uy, uz) = free_surface_fields();
    let iv = interpolate_particle_velocity(
        8,
        6,
        5,
        &ux,
        &uy,
        &uz,
        ParticleVector3 {
            x: -0.3,
            y: 6.7,
            z: 2.25,
        },
    )
    .unwrap();
    let p = |pos: [f64; 3], vel: [f64; 3], mass, response_time| ImmersedBoundaryParticle {
        position: ParticleVector3 {
            x: pos[0],
            y: pos[1],
            z: pos[2],
        },
        velocity: ParticleVector3 {
            x: vel[0],
            y: vel[1],
            z: vel[2],
        },
        mass,
        response_time,
    };
    let mut ib = vec![
        p([1.2, 2.5, 0.4], [0.0, 0.1, 0.0], 2.0, 0.5),
        p([7.9, 5.9, 4.9], [0.2, 0.0, -0.1], 1.0, 2.0),
    ];
    let mut res = ParticleCouplingResult::default();
    for _ in 0..3 {
        res = advance_immersed_boundary_particles(
            8,
            6,
            5,
            &ux,
            &uy,
            &uz,
            &mut ib,
            0.4,
            true,
            ParticleVector3 {
                x: 0.0,
                y: 0.0,
                z: -0.01,
            },
        )
        .unwrap();
    }
    let mut v = vec![iv.x, iv.y, iv.z];
    for q in &ib {
        v.extend([
            q.position.x,
            q.position.y,
            q.position.z,
            q.velocity.x,
            q.velocity.y,
            q.velocity.z,
        ]);
    }
    v.extend([
        res.total_reaction_force.x,
        res.total_reaction_force.y,
        res.total_reaction_force.z,
    ]);
    close("ib", &v, 1e-13, 1e-16);
    let f32s = |a: &[f32]| a.iter().map(|x| f64::from(*x)).collect::<Vec<_>>();
    close("ib_ax", &f32s(&res.acceleration_x), 1e-7, 0.0);
    close("ib_az", &f32s(&res.acceleration_z), 1e-7, 0.0);
}

#[test]
fn visualization_matches_upstream_byte_for_byte() {
    let (vx, vy, vz) = (7usize, 5usize, 4usize);
    let mut u2 = vec![Vector2::default(); vx * vy];
    let mut s2 = vec![0.0; vx * vy];
    for y in 0..vy {
        for x in 0..vx {
            let (xf, yf) = (x as f64, y as f64);
            u2[y * vx + x] = Vector2 {
                x: (0.7 * xf).sin() * (0.4 * yf).cos() + 0.1,
                y: -0.3 * xf + 0.2 * yf * yf,
            };
            s2[y * vx + x] =
                1.0 + 0.1 * xf - 1e-7 * yf * yf * yf + if x == 3 { 1e-30 } else { 0.0 };
        }
    }
    close(
        "vorticity",
        &vorticity_2d(vx, vy, &u2, 0.5, 0.25).unwrap(),
        1e-14,
        0.0,
    );
    let mut u3 = vec![Vector3::default(); vx * vy * vz];
    for z in 0..vz {
        for y in 0..vy {
            for x in 0..vx {
                let (xf, yf, zf) = (x as f64, y as f64, z as f64);
                u3[(z * vy + y) * vx + x] = Vector3 {
                    x: (0.5 * yf + 0.2 * zf).sin(),
                    y: (0.3 * xf).cos() * zf,
                    z: 0.1 * xf * yf,
                };
            }
        }
    }
    close(
        "q_criterion",
        &q_criterion_3d(vx, vy, vz, &u3, 1.0, 0.5, 2.0).unwrap(),
        1e-13,
        1e-15,
    );
    close(
        "slice",
        &horizontal_slice(vx, vy, &s2, 3).unwrap(),
        0.0,
        0.0,
    );
    let sl: Vec<f64> = streamline_2d(vx, vy, &u2, Vector2 { x: 0.5, y: 0.5 }, 0.4, 30)
        .unwrap()
        .iter()
        .flat_map(|p| [p.x, p.y])
        .collect();
    close("streamline", &sl, 1e-14, 0.0);

    let n3 = vx * vy * vz;
    let s3: Vec<f64> = (0..n3)
        .map(|i| {
            0.001 * (i as f64) * (i as f64) - 3.25
                + if i == 5 { 123456789.0 } else { 0.0 }
                + if i == 6 { 1e-9 } else { 0.0 }
        })
        .collect();
    let solid: Vec<u8> = (0..n3).map(|i| u8::from(i % 3 == 0)).collect();
    let dir = std::env::temp_dir().join(format!("cfd-lbm-viz-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    write_vtk_structured_2d(dir.join("a.vtk"), vx, vy, &s2, &u2, "rho").unwrap();
    write_vtk_structured_3d(dir.join("b.vtk"), vx, vy, vz, &s3, &u3, &solid, "field").unwrap();
    write_json_state_3d(
        dir.join("nested/state.json"),
        vx,
        vy,
        vz,
        42,
        &s3,
        &solid,
        "speed",
    )
    .unwrap();
    write_pgm_scalar(dir.join("c.pgm"), vx, vy, &s2).unwrap();
    write_html_scalar_viewer(dir.join("d.html"), vx, vy, &s2).unwrap();
    same_text("vtk2d", &dir.join("a.vtk"));
    same_text("vtk3d", &dir.join("b.vtk"));
    same_text("json_state", &dir.join("nested/state.json"));
    same_text("pgm", &dir.join("c.pgm"));
    same_text("html", &dir.join("d.html"));
    assert!(!dir.join("nested/state.json.tmp").exists());
    let _ = std::fs::remove_dir_all(dir);
}

/// Expected strings are C `printf("%.*g")` output (glibc), which iostreams reproduce.
#[test]
fn general_format_matches_printf() {
    let cases = [
        (0.1, 17, "0.10000000000000001"),
        (1e-5, 9, "1e-05"),
        (123456789.0, 9, "123456789"),
        (1234567890.0, 9, "1.23456789e+09"),
        (-0.0, 6, "-0"),
        (2.5e-30, 17, "2.4999999999999999e-30"),
        (100.0, 3, "100"),
        (1000.0, 3, "1e+03"),
        (0.0001, 6, "0.0001"),
        (9.9999999e-5, 3, "0.0001"),
        (5e-324, 17, "4.9406564584124654e-324"),
        (1.7976931348623157e308, 17, "1.7976931348623157e+308"),
    ];
    for (v, p, expected) in cases {
        assert_eq!(format_general(v, p), expected, "{v:e} at {p}");
    }
}
