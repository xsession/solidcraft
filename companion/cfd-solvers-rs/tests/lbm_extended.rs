//! Analytic and conservation checks for the extended LBM units.
use cfd_solvers::{core::Vec3, io::*, lbm::*};
use std::f64::consts::PI;

#[test]
fn advanced_collisions_conserve_mass_and_decay_taylor_green_viscously() {
    let ac = AdvancedD2Q9Config {
        nx: 32,
        ny: 32,
        tau: 0.8,
    };
    let nu = (f64::from(ac.tau) - 0.5) / 3.0;
    let k2 = 2.0 * (2.0 * PI / 32.0).powi(2);
    let steps = 200;
    let expected = (-2.0 * nu * k2 * steps as f64).exp();
    let check = |name: &str, m0: f64, e0: f64, m1: f64, e1: f64| {
        // Populations are f32, so mass is conserved to single-precision round-off.
        assert!((m1 - m0).abs() < 1e-5 * m0, "{name} mass drift {}", m1 - m0);
        let ratio = e1 / e0;
        assert!(
            (ratio - expected).abs() < 0.03 * expected,
            "{name} energy ratio {ratio} vs {expected}"
        );
    };
    let mut r = RegularizedD2Q9Solver::new(ac).unwrap();
    r.initialize_taylor_green(0.02);
    let (m0, e0) = (r.mass(), r.kinetic_energy());
    r.step(steps);
    check("regularized", m0, e0, r.mass(), r.kinetic_energy());
    let mut c = CumulantD2Q9Solver::new(ac).unwrap();
    c.initialize_taylor_green(0.02);
    c.step(steps);
    check("cumulant", m0, e0, c.mass(), c.kinetic_energy());
    let mut m = MrtD2Q9Solver::new(ac).unwrap();
    m.initialize_taylor_green(0.02);
    m.step(steps);
    check("mrt", m0, e0, m.mass(), m.kinetic_energy());
    assert!(MrtD2Q9Solver::new(AdvancedD2Q9Config { tau: 0.5, ..ac }).is_err());
}

#[test]
fn passive_scalar_sine_mode_decays_at_diffusive_rate() {
    let (nx, d) = (64usize, 0.05);
    let mut s = PassiveScalarD2Q5::new(PassiveScalarConfig {
        nx,
        ny: 4,
        diffusivity: d,
        ux: 0.0,
        uy: 0.0,
    })
    .unwrap();
    s.initialize(|x, _| 1.0 + 0.1 * (2.0 * PI * x).sin());
    let total = s.total();
    let steps = 400;
    s.step(steps);
    let k = 2.0 * PI / nx as f64;
    let amplitude = s.scalar()[..nx]
        .iter()
        .enumerate()
        .map(|(i, v)| (v - 1.0) * (2.0 * PI * (i as f64 + 0.5) / nx as f64).sin())
        .sum::<f64>()
        * 2.0
        / nx as f64;
    let expected = 0.1 * (-d * k * k * steps as f64).exp();
    assert!(
        (amplitude - expected).abs() < 0.01 * expected,
        "{amplitude} vs {expected}"
    );
    assert!((s.total() - total).abs() < 1e-10 * total);
}

fn unit_cube(lo: f64, hi: f64) -> TriangleSurface {
    let p = |x: f64, y: f64, z: f64| {
        Vec3::new(
            if x > 0.5 { hi } else { lo },
            if y > 0.5 { hi } else { lo },
            if z > 0.5 { hi } else { lo },
        )
    };
    let quad =
        |a: Vec3, b: Vec3, c: Vec3, d: Vec3| [Triangle3 { a, b, c }, Triangle3 { a, b: c, c: d }];
    let mut t = Vec::new();
    t.extend(quad(
        p(0., 0., 0.),
        p(0., 1., 0.),
        p(1., 1., 0.),
        p(1., 0., 0.),
    ));
    t.extend(quad(
        p(0., 0., 1.),
        p(1., 0., 1.),
        p(1., 1., 1.),
        p(0., 1., 1.),
    ));
    t.extend(quad(
        p(0., 0., 0.),
        p(1., 0., 0.),
        p(1., 0., 1.),
        p(0., 0., 1.),
    ));
    t.extend(quad(
        p(0., 1., 0.),
        p(0., 1., 1.),
        p(1., 1., 1.),
        p(1., 1., 0.),
    ));
    t.extend(quad(
        p(0., 0., 0.),
        p(0., 0., 1.),
        p(0., 1., 1.),
        p(0., 1., 0.),
    ));
    t.extend(quad(
        p(1., 0., 0.),
        p(1., 1., 0.),
        p(1., 1., 1.),
        p(1., 0., 1.),
    ));
    TriangleSurface { triangles: t }
}

#[test]
fn voxelized_cube_has_exact_cell_count_and_rotation_by_full_turn_is_identity() {
    let grid = VoxelGrid3D {
        nx: 10,
        ny: 10,
        nz: 10,
        minimum: Vec3::new(0.0, 0.0, 0.0),
        maximum: Vec3::new(1.0, 1.0, 1.0),
    };
    let mask = voxelize_surface(&unit_cube(0.2, 0.7), &grid).unwrap();
    assert_eq!(mask.iter().filter(|v| **v == 1).count(), 125);
    let turned = transform_surface(
        &unit_cube(0.2, 0.7),
        Vec3::default(),
        Vec3::new(1.0, 2.0, 3.0),
        2.0 * PI,
    );
    assert_eq!(voxelize_surface(&turned, &grid).unwrap(), mask);
    assert!(voxelize_surface(&TriangleSurface::default(), &grid).is_err());
}

#[test]
fn free_surface_volume_is_conserved_and_planar_normals_point_into_liquid() {
    let (nx, ny, nz) = (10, 4, 4);
    let n = nx * ny * nz;
    let fill: Vec<f64> = (0..n)
        .map(|i| {
            if i % nx < 4 {
                1.0
            } else if i % nx == 4 {
                0.5
            } else {
                0.0
            }
        })
        .collect();
    let mut fs = make_free_surface_field(nx, ny, nz, &fill, 1e-8).unwrap();
    let normals = free_surface_normals(&fs);
    assert!((normals[4].x + 1.0).abs() < 1e-14 && normals[4].y == 0.0);
    assert_eq!(fs.phase[4], FreeSurfacePhase::Interface);
    let v0 = free_surface_volume(&fs);
    let (ux, uy, uz) = (vec![0.4f32; n], vec![0.1f32; n], vec![-0.2f32; n]);
    for _ in 0..20 {
        advect_free_surface_vof(&mut fs, &ux, &uy, &uz, 1.0, 1e-8).unwrap();
    }
    assert!((free_surface_volume(&fs) - v0).abs() < 1e-12);
    assert!(fs.fill.iter().all(|f| (0.0..=1.0).contains(f)));
    assert!(advect_free_surface_vof(&mut fs, &ux, &uy, &uz, 2.0, 1e-8).is_err());
}

#[test]
fn immersed_boundary_reaction_spreads_exactly_and_particle_relaxes() {
    let (nx, ny, nz) = (4, 4, 4);
    let n = nx * ny * nz;
    let (ux, uy, uz) = (vec![0.1f32; n], vec![0.0f32; n], vec![0.0f32; n]);
    let mut ps = vec![ImmersedBoundaryParticle {
        position: ParticleVector3 {
            x: 1.3,
            y: 2.7,
            z: 0.2,
        },
        ..Default::default()
    }];
    let r = advance_immersed_boundary_particles(
        nx,
        ny,
        nz,
        &ux,
        &uy,
        &uz,
        &mut ps,
        0.5,
        true,
        ParticleVector3::default(),
    )
    .unwrap();
    let spread: f64 = r.acceleration_x.iter().map(|v| f64::from(*v)).sum();
    assert!((spread - r.total_reaction_force.x).abs() < 1e-7);
    assert!((ps[0].velocity.x - 0.05).abs() < 1e-8);
    let one_way = advance_immersed_boundary_particles(
        nx,
        ny,
        nz,
        &ux,
        &uy,
        &uz,
        &mut ps,
        0.5,
        false,
        ParticleVector3::default(),
    )
    .unwrap();
    assert!(one_way.acceleration_x.is_empty());
    let mut bad = vec![ImmersedBoundaryParticle {
        mass: 0.0,
        ..Default::default()
    }];
    assert!(advance_immersed_boundary_particles(
        nx,
        ny,
        nz,
        &ux,
        &uy,
        &uz,
        &mut bad,
        0.5,
        true,
        ParticleVector3::default()
    )
    .is_err());
}

#[test]
fn darcy_limit_and_rigid_rotation_vorticity() {
    let a = darcy_forchheimer_acceleration(0.01, 0.0, 0.1, 0.5, 0.0).unwrap();
    assert!((a[0] + 0.002).abs() < 1e-15 && a[1] == 0.0);
    assert!(darcy_forchheimer_acceleration(0.0, 0.0, 0.1, 0.0, 0.0).is_err());
    let (nx, ny) = (6, 6);
    let u: Vec<Vector2> = (0..nx * ny)
        .map(|i| Vector2 {
            x: -((i / nx) as f64),
            y: (i % nx) as f64,
        })
        .collect();
    let w = vorticity_2d(nx, ny, &u, 1.0, 1.0).unwrap();
    assert!((w[2 * nx + 2] - 2.0).abs() < 1e-14);
}

fn in_place(nx: usize, ny: usize, nz: usize, tau: f32, ax: f32) -> InPlaceLbmConfig {
    InPlaceLbmConfig {
        nx,
        ny,
        nz,
        tau,
        acceleration_x: ax,
        ..Default::default()
    }
}

/// Starting from the same equilibrium, the in-place solver streams after colliding while the
/// two-lattice solver streams before, so the two agree on viscous decay rather than cell by cell.
#[test]
fn in_place_and_two_lattice_solvers_both_decay_taylor_green_viscously() {
    let c = in_place(32, 32, 1, 0.8, 0.0);
    let nu = (f64::from(c.tau) - 0.5) / 3.0;
    let k2 = 2.0 * (2.0 * PI / 32.0).powi(2);
    let steps = 200;
    let expected = (-2.0 * nu * k2 * steps as f64).exp();
    let mut e = D2Q9InPlaceSolver::new(c).unwrap();
    e.initialize_taylor_green(0.02);
    let (m0, e0) = (e.mass(), e.kinetic_energy());
    e.step(steps);
    assert!((e.mass() - m0).abs() < 1e-5 * m0);
    let ratio = e.kinetic_energy() / e0;
    assert!(
        (ratio - expected).abs() < 0.03 * expected,
        "in-place ratio {ratio} vs {expected}"
    );
    let mut p = D2Q9PullSolver::new(c, 0.0).unwrap();
    p.initialize_taylor_green(0.02);
    p.step(steps);
    let (me, mp) = (e.compute_macroscopic(), p.compute_macroscopic());
    let worst = me
        .ux
        .iter()
        .zip(&mp.ux)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    assert!(
        worst < 2e-4,
        "in-place vs two-lattice velocity difference {worst}"
    );
    assert_eq!(e.population_bytes() * 2, p.population_bytes());
}

#[test]
fn bounce_back_conserves_mass_and_channel_flow_develops_downstream() {
    let (nx, ny, nz) = (6, 9, 3);
    let mut s = D3Q19Solver::new(in_place(nx, ny, nz, 0.8, 1e-5)).unwrap();
    let mut solid = vec![0u8; s.cells()];
    for z in 0..nz {
        for x in 0..nx {
            solid[(z * ny) * nx + x] = 1;
            solid[(z * ny + ny - 1) * nx + x] = 1;
        }
    }
    s.set_solid_mask(&solid).unwrap();
    s.initialize_uniform(1.0, 0.0, 0.0, 0.0);
    let m0 = s.mass();
    s.step(300);
    assert!(
        (s.mass() - m0).abs() < 1e-4 * m0,
        "mass drift {}",
        s.mass() - m0
    );
    let m = s.compute_macroscopic();
    let centre = m.ux[(ny / 2) * nx + 2];
    let near_wall = m.ux[nx + 2];
    assert!(
        centre > 0.0 && centre > 2.0 * near_wall,
        "centre {centre}, near wall {near_wall}"
    );
    assert!(m
        .ux
        .iter()
        .zip(&solid)
        .filter(|(_, s)| **s != 0)
        .all(|(u, _)| *u == 0.0));
}

#[test]
fn compressed_codecs_round_trip_lattice_weights_and_track_the_f32_solver() {
    type S = ShiftedHalfCodec<D2Q9InPlaceDescriptor>;
    for d in 0..9 {
        let w = D2Q9InPlaceDescriptor::WEIGHTS[d];
        assert_eq!(S::decode(S::encode(w, d), d), w);
        assert_eq!(
            PackedDeviation16Codec::<D2Q9InPlaceDescriptor>::decode(
                PackedDeviation16Codec::<D2Q9InPlaceDescriptor>::encode(w, d),
                d
            ),
            w
        );
    }
    assert_eq!(
        IeeeHalfCodec::decode_value(IeeeHalfCodec::encode_value(1.0e6)),
        f32::INFINITY
    );
    let mut c =
        CompressedPopulationPullSolver::<D2Q9InPlaceDescriptor, S>::new(CompressedLbmConfig {
            nx: 16,
            ny: 16,
            nz: 1,
            tau: 0.8,
        })
        .unwrap();
    let mut r = PrecisionPullSolver::<D2Q9InPlaceDescriptor, f64>::new(TypedLbmConfig {
        nx: 16,
        ny: 16,
        nz: 1,
        tau: 0.8,
    })
    .unwrap();
    c.initialize_taylor_green(0.03);
    r.initialize_taylor_green(0.03);
    c.step(50);
    r.step(50);
    let (a, b) = (c.compute_macroscopic(), r.compute_macroscopic());
    let worst =
        a.ux.iter()
            .zip(&b.ux)
            .map(|(x, y)| (f64::from(*x) - y).abs())
            .fold(0.0, f64::max);
    assert!(worst < 2e-4, "shifted-half drift {worst}");
    assert_eq!(c.population_bytes(), 9 * 256 * 2);
}

#[test]
fn smagorinsky_raises_the_relaxation_time_only_under_shear() {
    let mut f = [0.0f32; 19];
    for (d, v) in f.iter_mut().enumerate() {
        *v = lattice_equilibrium::<D3Q19Descriptor>(d, 1.0, 0.01, 0.0, 0.0);
    }
    assert_eq!(
        smagorinsky_relaxation_time::<D3Q19Descriptor>(&f, 1.0, 0.01, 0.0, 0.0, 0.6),
        0.6
    );
    f[7] += 1e-3;
    f[8] += 1e-3;
    assert!(smagorinsky_relaxation_time::<D3Q19Descriptor>(&f, 1.0, 0.01, 0.0, 0.0, 0.6) > 0.6);
    assert!(D3Q19Solver::new(in_place(4, 4, 1, 0.6, 0.0)).is_err());
}
