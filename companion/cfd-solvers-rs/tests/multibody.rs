//! Invariant and analytic checks for the multibody module.
#![allow(clippy::needless_range_loop, clippy::type_complexity)]
use cfd_solvers::{fvm::make_cartesian_hexa_mesh, multibody::*};
use std::sync::Arc;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

#[test]
fn torque_free_body_conserves_energy_and_momentum() {
    let mut b = RigidBody {
        mass: 2.0,
        inertia_diagonal: v(1.0, 2.0, 3.0),
        ..Default::default()
    };
    b.state.angular_velocity = v(0.3, 0.01, 0.2);
    b.state.linear_velocity = v(1.0, 0.0, 0.0);
    let energy = |b: &RigidBody| {
        let w = inverse_rotate(b.state.orientation, b.state.angular_velocity).unwrap();
        0.5 * (w.x * w.x + 2.0 * w.y * w.y + 3.0 * w.z * w.z)
    };
    let e0 = energy(&b);
    for _ in 0..2000 {
        integrate_rigid_body(&mut b, 1e-3, RigidBodyIntegrator::SemiImplicitEuler).unwrap();
    }
    assert!((energy(&b) - e0).abs() < 2e-3 * e0);
    assert!(
        (b.state.position.x - 2.0).abs() < 1e-12 && b.state.linear_velocity == v(1.0, 0.0, 0.0)
    );
}

#[test]
fn sphere_settles_on_plane_and_pendulum_keeps_length() {
    for method in [
        ContactMethod::SmoothPenalty,
        ContactMethod::NonsmoothImpulse,
    ] {
        let mut sys = RigidBodySystem::new(MultibodyConfig {
            contact_method: method,
            ..Default::default()
        })
        .unwrap();
        let mut b = RigidBody::default();
        b.state.position = v(0.0, 0.3, 0.0);
        let id = sys.add_body(b).unwrap();
        sys.add_sphere(id, 0.1).unwrap();
        sys.add_plane(PlaneShape::default());
        sys.run(3000).unwrap();
        let y = sys.bodies()[0].state.position.y;
        assert!(y > 0.08 && y < 0.105, "{method:?} {y}");
        assert!(sys.bodies()[0].state.linear_velocity.norm() < 0.05);
    }
    let mut sys = RigidBodySystem::new(MultibodyConfig {
        constraint_iterations: 30,
        ..Default::default()
    })
    .unwrap();
    let anchor = sys
        .add_body(RigidBody {
            fixed: true,
            ..Default::default()
        })
        .unwrap();
    let mut bob = RigidBody::default();
    bob.state.position = v(1.0, 0.0, 0.0);
    let id = sys.add_body(bob).unwrap();
    let j = make_spherical_joint(sys.bodies(), anchor, id, Vec3::ZERO).unwrap();
    sys.add_joint(Joint::Spherical(j));
    sys.run(500).unwrap();
    let r = sys.bodies()[1].state.position.norm();
    assert!((r - 1.0).abs() < 0.02, "{r}");
    assert!(sys.bodies()[1].state.position.y < -0.3);
}

#[test]
fn broad_phases_agree_with_brute_force() {
    let mut bodies = vec![RigidBody::default(); 40];
    let mut spheres = vec![];
    for i in 0..40 {
        let f = i as f64;
        bodies[i].state.position = v((0.731 * f).sin(), (0.413 * f).cos(), 0.3 * (0.17 * f).sin());
        spheres.push(SphereShape {
            body: i,
            radius: 0.12,
        });
    }
    let mut brute = vec![];
    for i in 0..40 {
        for j in i + 1..40 {
            let (a, b) = (
                sphere_aabb(&bodies[i], &spheres[i]).unwrap(),
                sphere_aabb(&bodies[j], &spheres[j]).unwrap(),
            );
            if a.minimum.x <= b.maximum.x
                && b.minimum.x <= a.maximum.x
                && a.minimum.y <= b.maximum.y
                && b.minimum.y <= a.maximum.y
                && a.minimum.z <= b.maximum.z
                && b.minimum.z <= a.maximum.z
            {
                brute.push((i, j));
            }
        }
    }
    let mut sap: Vec<(usize, usize)> = broad_phase_sweep_and_prune(&bodies, &spheres)
        .unwrap()
        .into_iter()
        .map(|(a, b)| (a.min(b), a.max(b)))
        .collect();
    sap.sort_unstable();
    assert_eq!(sap, brute);
    assert_eq!(broad_phase_bvh(&bodies, &spheres).unwrap(), brute);
    assert_eq!(
        broad_phase_cell_linked(&bodies, &spheres, 0.25).unwrap(),
        brute
    );
}

#[test]
fn gjk_epa_cube_overlap_depth() {
    let mut bodies = vec![RigidBody::default(); 2];
    bodies[1].state.position = v(0.9, 0.0, 0.0);
    let cube = |b| ConvexHullShape {
        body: b,
        vertices: vec![
            v(-0.5, -0.5, -0.5),
            v(0.5, -0.5, -0.5),
            v(0.5, 0.5, -0.5),
            v(-0.5, 0.5, -0.5),
            v(-0.5, -0.5, 0.5),
            v(0.5, -0.5, 0.5),
            v(0.5, 0.5, 0.5),
            v(-0.5, 0.5, 0.5),
        ],
    };
    let r = convex_gjk_epa_contact(&bodies, &cube(0), &cube(1), 32, 64, 1e-9).unwrap();
    assert!(r.intersect && (r.contact.penetration - 0.1).abs() < 1e-9);
    assert!((r.contact.normal.x + 1.0).abs() < 1e-9);
    bodies[1].state.position = v(1.2, 0.0, 0.0);
    assert!(
        !convex_gjk_epa_contact(&bodies, &cube(0), &cube(1), 32, 64, 1e-9)
            .unwrap()
            .intersect
    );
}

#[test]
fn newmark_conserves_linear_oscillator_energy() {
    let force: GeneralizedForceFunction =
        Arc::new(|q: &[f64], _: &[f64], _: f64, f: &mut [f64]| {
            f[0] = -9.0 * q[0];
            Ok(())
        });
    let mut it = ImplicitGeneralizedIntegrator::new(
        vec![1.0],
        1,
        force,
        ImplicitSecondOrderConfig::newmark(0.05, 0.25, 0.5),
    )
    .unwrap();
    it.initialize(vec![1.0], vec![], 0.0).unwrap();
    it.run(400).unwrap();
    let e = 0.5 * it.v()[0] * it.v()[0] + 4.5 * it.q()[0] * it.q()[0];
    assert!((e - 4.5).abs() < 1e-9);
    assert!(ImplicitSecondOrderConfig::generalized_alpha(0.1, 1.5).is_err());
}

#[test]
fn articulated_pendulum_small_angle_period() {
    let mut arm = ArticulatedSystem::new(v(0.0, -9.81, 0.0));
    let rot = quaternion_from_axis_angle(v(0.0, 0.0, 1.0), -std::f64::consts::FRAC_PI_2).unwrap();
    arm.add_link(ArticulatedLink {
        joint_to_com_local: v(1.0, 0.0, 0.0),
        reference_orientation: rot,
        body: RigidBody {
            mass: 1.0,
            inertia_diagonal: v(1e-9, 1e-9, 1e-9),
            ..Default::default()
        },
        q: 0.05,
        ..Default::default()
    })
    .unwrap();
    // Point mass on a unit rod: hanging straight down at q = 0.
    let p = arm.links()[0].body.state.position;
    assert!((p.y + (0.05f64).cos()).abs() < 1e-9);
    let dt = 1e-4;
    let mut t = 0.0;
    let mut prev = arm.links()[0].q;
    let mut crossings = vec![];
    while crossings.len() < 3 && t < 10.0 {
        arm.step(dt).unwrap();
        t += dt;
        let q = arm.links()[0].q;
        if prev > 0.0 && q <= 0.0 {
            crossings.push(t);
        }
        prev = q;
    }
    let period = crossings[2] - crossings[1];
    assert!(
        (period - 2.0 * std::f64::consts::PI / 9.81f64.sqrt()).abs() < 5e-3,
        "{period}"
    );
}

#[test]
fn coupling_and_distributed_conservation() {
    let mesh = make_cartesian_hexa_mesh(2, 2, 2, 0.2, 0.2, 0.2).unwrap();
    let n = mesh.cell_count();
    let uf = vec![v(0.3, 0.0, 0.1); n];
    let mut parts = vec![RigidBody::default(); 3];
    let mut ps = vec![];
    for i in 0..3 {
        parts[i].state.position = v(0.03 + 0.06 * i as f64, 0.05, 0.12);
        ps.push(SphereShape {
            body: i,
            radius: 0.01,
        });
    }
    let r = unresolved_cfd_dem_coupling(
        &mesh,
        &parts,
        &ps,
        &uf,
        &vec![v(1.0, 2.0, 3.0); n],
        &vec![v(0.0, 0.0, 2.0); n],
        &UnresolvedCfdDemConfig::default(),
    )
    .unwrap();
    let mut total = Vec3::ZERO;
    for f in &r.particle_force {
        total += *f;
    }
    for (c, s) in r.fluid_momentum_source.iter().enumerate() {
        total += *s * mesh.cells()[c].volume;
    }
    assert!(total.norm() < 1e-15);
    // Uniform pressure gives zero net resolved force with the octahedral rule.
    let load = integrate_resolved_surface_load(
        &mesh,
        &vec![5.0; n],
        &[],
        &sphere_surface_quadrature(&parts[0], 0.02, 6).unwrap(),
        parts[0].state.position,
    )
    .unwrap();
    assert!(load.force.norm() < 1e-15);

    // Two emulated ranks reproduce the single-rank result.
    let make = |owner: &dyn Fn(f64) -> i32| -> Vec<DistributedDemParticle> {
        (0..6)
            .map(|i| {
                let mut b = RigidBody {
                    mass: 0.01,
                    inertia_diagonal: v(1e-6, 1e-6, 1e-6),
                    ..Default::default()
                };
                b.state.position = v(0.0125 * i as f64 + 0.01, 0.0, 0.0);
                b.state.linear_velocity = v(0.02 * (3.0 - i as f64), 0.0, 0.0);
                make_distributed_dem_particle(
                    i as u64,
                    &b,
                    &SphereShape {
                        body: 0,
                        radius: 0.0068,
                    },
                    owner(b.state.position.x),
                )
                .unwrap()
            })
            .collect()
    };
    let model = HertzMindlinContactModel {
        normal_stiffness: 2e4,
        ..Default::default()
    };
    let mut single = make(&|_| 0);
    let mut hs = vec![];
    for _ in 0..5 {
        let c = evaluate_distributed_dem_contacts(&single, &[], &mut hs, 0, 1e-4, &model).unwrap();
        integrate_distributed_dem_particles(&mut single, &c.local, Vec3::ZERO, 1e-4).unwrap();
    }
    let all = make(&|x| if x < 0.05 { 0 } else { 1 });
    let (mut o0, mut o1): (Vec<_>, Vec<_>) = all.into_iter().partition(|p| p.owner_rank == 0);
    let (mut h0, mut h1) = (vec![], vec![]);
    for _ in 0..5 {
        let (g0, g1) = (o1.clone(), o0.clone());
        let c0 = evaluate_distributed_dem_contacts(&o0, &g0, &mut h0, 0, 1e-4, &model).unwrap();
        let c1 = evaluate_distributed_dem_contacts(&o1, &g1, &mut h1, 1, 1e-4, &model).unwrap();
        let (mut f0, mut f1) = (c0.local.clone(), c1.local.clone());
        accumulate_distributed_dem_forces(&mut f1, &c0.remote).unwrap();
        accumulate_distributed_dem_forces(&mut f0, &c1.remote).unwrap();
        integrate_distributed_dem_particles(&mut o0, &f0, Vec3::ZERO, 1e-4).unwrap();
        integrate_distributed_dem_particles(&mut o1, &f1, Vec3::ZERO, 1e-4).unwrap();
    }
    for p in o0.iter().chain(&o1) {
        let s = single.iter().find(|q| q.global_id == p.global_id).unwrap();
        assert!(
            (p.state.position - s.state.position).norm() < 1e-15
                && (p.state.linear_velocity - s.state.linear_velocity).norm() < 1e-12
        );
    }
}
