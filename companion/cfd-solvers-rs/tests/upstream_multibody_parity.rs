//! Value-by-value comparison with the C++ upstream multibody/*. Goldens come from
//! `parity/oracle_multibody.cpp`.
#![allow(clippy::needless_range_loop, clippy::type_complexity)]
use cfd_solvers::{fvm::make_cartesian_hexa_mesh, multibody::*};
use serde_json::Value;
use std::sync::Arc;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-multibody.json")).unwrap()
}

fn num(v: &Value) -> f64 {
    match v.as_str() {
        Some("inf") => f64::INFINITY,
        Some("-inf") => f64::NEG_INFINITY,
        Some(_) => f64::NAN,
        None => v.as_f64().unwrap(),
    }
}

#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(num)
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        if c.is_infinite() {
            assert_eq!(r, c, "{name}[{i}]");
            continue;
        }
        assert!(
            (r - c).abs() <= rel * c.abs() + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn v3(o: &mut Vec<f64>, a: Vec3) {
    o.extend([a.x, a.y, a.z]);
}
fn q4(o: &mut Vec<f64>, q: Quaternion) {
    o.extend([q.w, q.x, q.y, q.z]);
}
fn st(o: &mut Vec<f64>, b: &RigidBody) {
    v3(o, b.state.position);
    q4(o, b.state.orientation);
    v3(o, b.state.linear_velocity);
    v3(o, b.state.angular_velocity);
}
fn b(x: bool) -> f64 {
    if x {
        1.0
    } else {
        0.0
    }
}

#[test]
fn rigid_bodies_and_joints() {
    let mut m = vec![];
    let q = quaternion_from_axis_angle(v(1.0, 2.0, 0.5), 0.8).unwrap();
    q4(&mut m, q);
    v3(&mut m, rotate(q, v(0.3, -1.0, 2.0)).unwrap());
    v3(&mut m, inverse_rotate(q, v(0.3, -1.0, 2.0)).unwrap());
    v3(&mut m, any_perpendicular(v(0.9, 0.1, 0.2)));
    v3(&mut m, any_perpendicular(v(0.1, 0.9, 0.2)));
    let mut qq = q;
    integrate_orientation(&mut qq, v(0.5, -1.0, 2.0), 0.01).unwrap();
    q4(&mut m, qq);
    for integ in [
        RigidBodyIntegrator::SemiImplicitEuler,
        RigidBodyIntegrator::VelocityVerlet,
    ] {
        let mut bd = RigidBody {
            mass: 2.0,
            inertia_diagonal: v(1.0, 2.0, 3.0),
            ..Default::default()
        };
        bd.state.orientation = q;
        bd.state.angular_velocity = v(1.0, 2.0, 0.5);
        bd.state.linear_velocity = v(0.1, 0.0, 0.0);
        for _ in 0..5 {
            let p = bd.state.position + v(0.2, 0.0, 0.0);
            bd.add_force_at_point(v(0.0, 1.0, 0.5), p);
            bd.add_torque(v(0.1, 0.0, 0.0));
            integrate_rigid_body(&mut bd, 0.01, integ).unwrap();
        }
        let p = bd.state.position + v(0.0, 0.3, 0.0);
        bd.apply_impulse(v(0.5, 0.0, 0.0), p).unwrap();
        st(&mut m, &bd);
    }
    close("rigid", &m, 1e-13, 1e-15);

    let mut bodies: Vec<RigidBody> = (0..5)
        .map(|i| {
            let f = i as f64;
            let mut bd = RigidBody {
                mass: 1.0 + f,
                inertia_diagonal: v(0.1 + 0.02 * f, 0.2, 0.15),
                ..Default::default()
            };
            bd.state.position = v(0.5 * f, 0.1 * f, 0.0);
            bd.state.orientation = quaternion_from_axis_angle(v(0.0, 0.0, 1.0), 0.1 * f).unwrap();
            bd.state.linear_velocity = v(0.1 * f, -0.2, 0.05 * f);
            bd.state.angular_velocity = v(0.0, 0.1 * f, 0.3);
            bd
        })
        .collect();
    bodies[0].fixed = true;
    let joints = vec![
        Joint::Distance(
            make_distance_constraint(&bodies, 0, 1, v(0.1, 0.0, 0.0), v(0.4, 0.1, 0.0), None)
                .unwrap(),
        ),
        Joint::Spherical(make_spherical_joint(&bodies, 1, 2, v(0.75, 0.15, 0.0)).unwrap()),
        Joint::Revolute(
            make_revolute_joint(&bodies, 2, 3, v(1.25, 0.25, 0.0), v(0.0, 0.0, 1.0)).unwrap(),
        ),
        Joint::Prismatic(
            make_prismatic_joint(&bodies, 3, 4, v(1.75, 0.35, 0.0), v(1.0, 0.2, 0.0)).unwrap(),
        ),
        Joint::Fixed(make_fixed_joint(&bodies, 1, 3, v(1.0, 0.2, 0.0)).unwrap()),
        Joint::Gear(
            make_gear_constraint(&bodies, 2, 4, v(0.0, 0.0, 1.0), v(0.0, 1.0, 0.0), 2.0).unwrap(),
        ),
        Joint::LinearMotor(make_linear_motor(&bodies, 0, 4, v(1.0, 0.0, 0.0), 0.3, 5.0).unwrap()),
        Joint::AngularMotor(
            make_angular_motor(&bodies, 1, 2, v(0.0, 0.0, 1.0), 1.0, f64::INFINITY).unwrap(),
        ),
    ];
    let mut jv = vec![];
    let mut rows = vec![];
    for j in &joints {
        for r in build_constraint_rows(&bodies, j, 0.01, 0.2).unwrap() {
            v3(&mut jv, r.linear_a);
            v3(&mut jv, r.angular_a);
            v3(&mut jv, r.linear_b);
            v3(&mut jv, r.angular_b);
            jv.extend([r.bias, r.lower, r.upper]);
            rows.push(r);
        }
    }
    solve_constraint_rows(&mut bodies, &mut rows, 12).unwrap();
    jv.extend(rows.iter().map(|r| r.impulse));
    for bd in &bodies {
        st(&mut jv, bd);
    }
    close("joints", &jv, 1e-11, 1e-13);
}

#[test]
fn systems_and_dem() {
    for (method, name) in [
        (ContactMethod::SmoothPenalty, "sys_penalty"),
        (ContactMethod::NonsmoothImpulse, "sys_impulse"),
    ] {
        let mut mc = MultibodyConfig {
            contact_method: method,
            dt: 2e-3,
            ..Default::default()
        };
        mc.integrator = if method == ContactMethod::SmoothPenalty {
            RigidBodyIntegrator::VelocityVerlet
        } else {
            RigidBodyIntegrator::SemiImplicitEuler
        };
        mc.impulse.restitution = 0.3;
        let mut sys = RigidBodySystem::new(mc).unwrap();
        for i in 0..4 {
            let f = i as f64;
            let mut bd = RigidBody {
                inertia_diagonal: v(0.05, 0.05, 0.05),
                mass: 1.0 + 0.5 * f,
                ..Default::default()
            };
            bd.state.position = v(0.3 * f, 0.45 + 0.15 * f, 0.05 * f);
            bd.state.linear_velocity = v(0.2 - 0.1 * f, -0.5, 0.0);
            bd.state.angular_velocity = v(0.0, 0.0, 0.5 * f);
            let id = sys.add_body(bd).unwrap();
            sys.add_sphere(id, 0.2).unwrap();
        }
        sys.add_plane(PlaneShape {
            normal: v(0.0, 1.0, 0.0),
            offset: 0.0,
        });
        sys.add_triangle_mesh(TriangleMeshShape {
            body: None,
            vertices: vec![
                v(-1.0, 0.1, -1.0),
                v(2.0, 0.1, -1.0),
                v(2.0, 0.4, 1.0),
                v(-1.0, 0.4, 1.0),
            ],
            triangles: vec![[0, 1, 2], [0, 2, 3]],
            feature_namespace: 0,
        })
        .unwrap();
        let (p0, p1) = (
            sys.bodies()[0].state.position,
            sys.bodies()[1].state.position,
        );
        let dcon = make_distance_constraint(sys.bodies(), 0, 1, p0, p1, None).unwrap();
        sys.add_joint(Joint::Distance(dcon));
        sys.run(150).unwrap();
        let mut sv = vec![sys.time(), sys.total_kinetic_energy().unwrap()];
        v3(&mut sv, sys.linear_momentum());
        for bd in sys.bodies() {
            st(&mut sv, bd);
        }
        close(name, &sv, 1e-8, 1e-11);
    }
    for (hist, name) in [(true, "dem_hist"), (false, "dem_plain")] {
        let mut dc = DemConfig {
            dt: 5e-5,
            history_dependent_tangential: hist,
            ..Default::default()
        };
        dc.contact.cohesion_force = 0.01;
        let mut dem = ExplicitDemSystem::new(dc).unwrap();
        for i in 0..6 {
            let f = i as f64;
            dem.add_particle(
                v(
                    0.021 * f,
                    0.012 + 0.005 * (i % 2) as f64,
                    0.002 * (i % 3) as f64,
                ),
                0.011,
                2500.0,
                v(0.05 * (3.0 - f), -0.1, 0.0),
            )
            .unwrap();
        }
        dem.add_plane(PlaneShape {
            normal: v(0.0, 1.0, 0.0),
            offset: 0.0,
        });
        let bm = BondedParticleModel {
            tensile_failure_force: 0.5,
            shear_failure_force: 0.3,
            ..Default::default()
        };
        dem.add_bond(0, 1, bm, 0.0).unwrap();
        dem.add_bond(2, 3, bm, 0.02).unwrap();
        dem.add_bond(4, 5, BondedParticleModel::default(), 0.0)
            .unwrap();
        let mut dv = vec![];
        for _ in 0..6 {
            dem.run(40).unwrap();
            let bs = dem.last_bond_stats();
            dv.extend([
                bs.active as f64,
                bs.damaged as f64,
                bs.broken_total as f64,
                dem.contact_history().states().len() as f64,
            ]);
        }
        for bd in dem.bodies() {
            st(&mut dv, bd);
        }
        for bo in dem.bonds() {
            dv.extend([bo.damage, b(bo.broken)]);
            v3(&mut dv, bo.tangential_displacement);
        }
        close(name, &dv, 1e-8, 1e-12);
    }
}

#[test]
fn broad_phase_and_gjk() {
    let mut cloud = vec![RigidBody::default(); 30];
    let mut sph = vec![];
    for i in 0..30 {
        let f = i as f64;
        cloud[i].state.position = v((0.37 * f) % 1.0, (0.61 * f) % 0.9, (0.13 * f) % 0.4);
        sph.push(SphereShape {
            body: i,
            radius: 0.06 + 0.002 * (i % 5) as f64,
        });
    }
    sph.push(SphereShape {
        body: 3,
        radius: 0.05,
    });
    // Upstream's BVH traversal drops pairs whose lower index lands in the right subtree (and its
    // pair set depends on libstdc++ nth_element); the port emits every overlapping pair, so the
    // upstream BVH list must be a subset of ours, and ours must equal sweep-and-prune.
    let g = goldens();
    let cpp: Vec<f64> = g["broad"].as_array().unwrap().iter().map(num).collect();
    let mut segments = cpp.split(|x| *x == -1.0);
    let to_pairs = |s: &[f64]| -> Vec<(usize, usize)> {
        s.chunks(2)
            .map(|c| (c[0] as usize, c[1] as usize))
            .collect()
    };
    let (seg_sap, seg_bvh, seg_cell, seg_contacts) = (
        segments.next().unwrap(),
        segments.next().unwrap(),
        segments.next().unwrap(),
        segments.next().unwrap(),
    );
    let mut sap = broad_phase_sweep_and_prune(&cloud, &sph).unwrap();
    sap.sort_unstable();
    assert_eq!(to_pairs(seg_sap), sap);
    let bvh = broad_phase_bvh(&cloud, &sph).unwrap();
    let mut canonical: Vec<(usize, usize)> =
        sap.iter().map(|&(a, c)| (a.min(c), a.max(c))).collect();
    canonical.sort_unstable();
    assert_eq!(bvh, canonical);
    assert!(to_pairs(seg_bvh).iter().all(|p| bvh.contains(p)));
    assert_eq!(
        to_pairs(seg_cell),
        broad_phase_cell_linked(&cloud, &sph, 0.2).unwrap()
    );
    let mut cv = vec![];
    for c in sphere_contacts(&cloud, &sph, &sap).unwrap() {
        v3(&mut cv, c.point);
        v3(&mut cv, c.normal);
        cv.push(c.penetration);
    }
    assert_eq!(cv.len(), seg_contacts.len());
    for (r, c) in cv.iter().zip(seg_contacts) {
        assert!((r - c).abs() <= 1e-13 * c.abs() + 1e-15);
    }

    let mut hulls = vec![RigidBody::default(); 2];
    hulls[1].state.position = v(0.8, 0.3, 0.1);
    hulls[1].state.orientation = quaternion_from_axis_angle(v(0.3, 0.2, 1.0), 0.6).unwrap();
    let ha = ConvexHullShape {
        body: 0,
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
    let hb = ConvexHullShape {
        body: 1,
        vertices: vec![
            v(0.0, 0.0, 0.6),
            v(0.5, 0.0, -0.3),
            v(-0.4, 0.4, -0.3),
            v(-0.3, -0.5, -0.3),
        ],
    };
    let g = convex_gjk_epa_contact(&hulls, &ha, &hb, 32, 64, 1e-9).unwrap();
    let mut gv = vec![
        b(g.intersect),
        g.gjk_iterations as f64,
        g.epa_iterations as f64,
        g.contact.penetration,
    ];
    v3(&mut gv, g.contact.point);
    v3(&mut gv, g.contact.normal);
    hulls[1].state.position = v(2.5, 0.0, 0.0);
    let g2 = convex_gjk_epa_contact(&hulls, &ha, &hb, 32, 64, 1e-9).unwrap();
    gv.extend([b(g2.intersect), g2.gjk_iterations as f64]);
    let bx = convex_aabb(&hulls[1], &hb).unwrap();
    v3(&mut gv, bx.minimum);
    v3(&mut gv, bx.maximum);
    close("gjk", &gv, 1e-12, 1e-14);
}

#[test]
fn implicit_articulated_and_flexible() {
    let duffing: GeneralizedForceFunction =
        Arc::new(|x: &[f64], vv: &[f64], t: f64, f: &mut [f64]| {
            f[0] = -4.0 * x[0] - 0.5 * x[0] * x[0] * x[0] + 0.3 * (x[1] - x[0]) - 0.05 * vv[0]
                + t.sin();
            f[1] = -2.0 * x[1] - 0.3 * (x[1] - x[0]) - 0.02 * vv[1];
            Ok(())
        });
    let mut iv = vec![];
    for cfg in [
        ImplicitSecondOrderConfig::newmark(0.05, 0.25, 0.5),
        ImplicitSecondOrderConfig::generalized_alpha(0.05, 0.6).unwrap(),
        ImplicitSecondOrderConfig::hht(0.05, -0.2).unwrap(),
    ] {
        let mut it =
            ImplicitGeneralizedIntegrator::new(vec![2.0, 0.1, 0.1, 1.0], 2, duffing.clone(), cfg)
                .unwrap();
        it.initialize(vec![0.5, -0.2], vec![0.0, 0.3], 0.0).unwrap();
        it.run(20).unwrap();
        iv.extend([
            it.q()[0],
            it.q()[1],
            it.v()[0],
            it.v()[1],
            it.a()[0],
            it.a()[1],
            it.time(),
            it.last_newton_iterations() as f64,
        ]);
    }
    close("implicit", &iv, 1e-8, 1e-10);

    let mut arm = ArticulatedSystem::new(v(0.0, -9.81, 0.0));
    let l0 = ArticulatedLink {
        parent_joint_offset: v(0.0, 1.0, 0.0),
        joint_to_com_local: v(0.5, 0.0, 0.0),
        body: RigidBody {
            mass: 2.0,
            inertia_diagonal: v(0.05, 0.2, 0.2),
            ..Default::default()
        },
        q: 0.3,
        damping: 0.05,
        ..Default::default()
    };
    let l1 = ArticulatedLink {
        parent: Some(0),
        parent_joint_offset: v(0.5, 0.0, 0.0),
        joint_to_com_local: v(0.4, 0.0, 0.0),
        body: RigidBody {
            mass: 1.0,
            inertia_diagonal: v(0.02, 0.1, 0.1),
            ..Default::default()
        },
        q: -0.2,
        qd: 0.5,
        ..Default::default()
    };
    let l2 = ArticulatedLink {
        parent: Some(1),
        joint_type: ReducedJointType::Prismatic,
        axis_parent: v(1.0, 0.0, 0.0),
        parent_joint_offset: v(0.4, 0.0, 0.0),
        joint_to_com_local: v(0.1, 0.0, 0.0),
        body: RigidBody {
            mass: 0.5,
            inertia_diagonal: v(0.01, 0.01, 0.01),
            ..Default::default()
        },
        q: 0.05,
        ..Default::default()
    };
    arm.add_link(l0).unwrap();
    arm.add_link(l1).unwrap();
    arm.add_link(l2).unwrap();
    arm.set_effort(1, 0.4).unwrap();
    let mut av = arm.mass_matrix().unwrap();
    av.extend(arm.generalized_forces());
    arm.run(50, 0.002).unwrap();
    for l in arm.links() {
        av.extend([l.q, l.qd, l.qdd]);
        st(&mut av, &l.body);
    }
    av.push(arm.kinetic_energy().unwrap());
    close("articulated", &av, 1e-10, 1e-13);

    let node = |x, y, z| FlexibleInterfaceNode {
        reference_local: v(x, y, z),
    };
    let mut flex = FlexibleBodyInterface::new(
        vec![
            node(0.0, 0.0, 0.0),
            node(1.0, 0.0, 0.0),
            node(1.0, 1.0, 0.0),
        ],
        vec![
            FlexibleMode {
                shape_local: vec![v(0.0, 0.0, 0.1), v(0.0, 0.1, 0.2), v(0.1, 0.0, 0.0)],
                modal_mass: 2.0,
                modal_stiffness: 5.0,
                modal_damping: 0.1,
            },
            FlexibleMode {
                shape_local: vec![v(0.05, 0.0, 0.0), v(0.0, 0.0, 0.0), v(0.0, 0.05, 0.05)],
                modal_mass: 1.0,
                modal_stiffness: 3.0,
                modal_damping: 0.0,
            },
        ],
    )
    .unwrap();
    flex.set_modal_state(vec![0.3, -0.2], vec![0.1, 0.05])
        .unwrap();
    let mut frame = RigidBody::default();
    frame.state.position = v(1.0, 2.0, 3.0);
    frame.state.orientation = quaternion_from_axis_angle(v(1.0, 2.0, 0.5), 0.8).unwrap();
    frame.state.linear_velocity = v(0.1, 0.0, 0.0);
    frame.state.angular_velocity = v(0.0, 0.0, 1.0);
    let mut fv = vec![];
    for p in flex.node_positions(&frame).unwrap() {
        v3(&mut fv, p);
    }
    for p in flex.node_velocities(&frame).unwrap() {
        v3(&mut fv, p);
    }
    let ld = flex
        .apply_nodal_forces(
            &mut frame,
            &[v(1.0, 0.0, 0.0), v(0.0, 2.0, 0.0), v(0.0, 0.0, -1.0)],
        )
        .unwrap();
    v3(&mut fv, ld.body_force);
    v3(&mut fv, ld.body_torque);
    fv.extend(&ld.modal_force);
    v3(&mut fv, frame.torque);
    close("flex", &fv, 1e-13, 1e-15);
}

#[test]
fn cfd_dem_coupling_and_distributed_dem() {
    let mesh = make_cartesian_hexa_mesh(3, 2, 2, 0.3, 0.2, 0.2).unwrap();
    let n = mesh.cell_count();
    let (mut pr, mut ss, mut uf, mut gp, mut vo) = (
        vec![0.0; n],
        vec![SymmetricStress3::default(); n],
        vec![Vec3::ZERO; n],
        vec![Vec3::ZERO; n],
        vec![Vec3::ZERO; n],
    );
    for c in 0..n {
        let x = mesh.cells()[c].center;
        let cf = c as f64;
        pr[c] = 100.0 + 50.0 * x.x - 20.0 * x.y;
        ss[c] = SymmetricStress3 {
            xx: 0.1 * cf,
            yy: 0.2,
            zz: -0.1,
            xy: 0.05 * cf,
            xz: 0.01,
            yz: -0.02,
        };
        uf[c] = v(0.5 + x.y, 0.1 * x.z, -0.2);
        gp[c] = v(50.0, -20.0, cf);
        vo[c] = v(0.0, 0.3 * cf, 1.0);
    }
    let mut ball = RigidBody::default();
    ball.state.position = v(0.15, 0.1, 0.1);
    let mut cv = vec![];
    for cnt in [6, 15] {
        let quad = sphere_surface_quadrature(&ball, 0.04, cnt).unwrap();
        let load =
            integrate_resolved_surface_load(&mesh, &pr, &ss, &quad, ball.state.position).unwrap();
        v3(&mut cv, load.force);
        v3(&mut cv, load.torque);
        cv.push(quad.len() as f64);
    }
    let mut parts = vec![RigidBody::default(); 4];
    let mut ps = vec![];
    for i in 0..4 {
        let f = i as f64;
        parts[i].state.position = v(0.05 + 0.07 * f, 0.04 + 0.03 * f, 0.1);
        parts[i].state.linear_velocity = v(0.1 * f, 0.0, 0.05);
        ps.push(SphereShape {
            body: i,
            radius: 0.01 + 0.003 * f,
        });
    }
    let ur = unresolved_cfd_dem_coupling(
        &mesh,
        &parts,
        &ps,
        &uf,
        &gp,
        &vo,
        &UnresolvedCfdDemConfig::default(),
    )
    .unwrap();
    cv.extend(&ur.solid_volume_fraction);
    cv.extend(&ur.void_fraction);
    for f in ur.particle_force.iter().chain(&ur.fluid_momentum_source) {
        v3(&mut cv, *f);
    }
    let uc = UnresolvedCfdDemConfig {
        include_saffman_lift: false,
        dynamic_viscosity: 1e-5,
        ..Default::default()
    };
    for f in unresolved_cfd_dem_coupling(&mesh, &parts, &ps, &uf, &[], &[], &uc)
        .unwrap()
        .particle_force
    {
        v3(&mut cv, f);
    }
    close("coupling", &cv, 1e-12, 1e-14);

    let dd = DemSlabDecomposition {
        x_min: 0.0,
        x_max: 0.1,
        ranks: 2,
    };
    let mut all = vec![];
    for i in 0..8 {
        let f = i as f64;
        let mut bd = RigidBody {
            mass: 0.01,
            inertia_diagonal: v(1e-6, 1e-6, 1e-6),
            ..Default::default()
        };
        bd.state.position = v(0.0125 * f + 0.003, 0.01 * (i % 2) as f64, 0.0);
        bd.state.linear_velocity = v(0.02 * (4.0 - f), 0.0, 0.01);
        bd.state.angular_velocity = v(0.0, 0.0, 0.3 * f);
        let owner = dd.owner_rank(bd.state.position.x).unwrap();
        all.push(
            make_distributed_dem_particle(
                100 + i as u64,
                &bd,
                &SphereShape {
                    body: 0,
                    radius: 0.0068,
                },
                owner,
            )
            .unwrap(),
        );
    }
    let mut ddv = vec![
        dd.slab_width().unwrap(),
        dd.owner_rank(0.07).unwrap() as f64,
        dd.slab_min(1).unwrap(),
        dd.slab_max(0).unwrap(),
    ];
    let plan = plan_dem_slab_exchange(&all, &dd, 0, 0.004).unwrap();
    for r in plan
        .migrate_indices_by_rank
        .iter()
        .chain(&plan.ghost_indices_by_rank)
    {
        ddv.extend(r.iter().map(|i| *i as f64));
        ddv.push(-1.0);
    }
    let (mut own0, mut own1): (Vec<_>, Vec<_>) = all.iter().partition(|p| p.owner_rank == 0);
    let ghosts = |own: &[DistributedDemParticle], low: bool| -> Vec<DistributedDemParticle> {
        own.iter()
            .filter(|p| {
                if low {
                    p.state.position.x - p.radius - 0.004 < 0.05
                } else {
                    p.state.position.x + p.radius + 0.004 > 0.05
                }
            })
            .copied()
            .collect()
    };
    let mut ghost0 = ghosts(&own1, true);
    let mut ghost1 = ghosts(&own0, false);
    let (mut h0, mut h1) = (vec![], vec![]);
    let mut b0 = vec![
        make_distributed_dem_bond(&own0[2], &own1[0], BondedParticleModel::default(), 0.0).unwrap(),
        make_distributed_dem_bond(&own0[0], &own0[1], BondedParticleModel::default(), 0.0).unwrap(),
    ];
    let hm = HertzMindlinContactModel {
        normal_stiffness: 2e4,
        tangential_stiffness: 1e4,
        ..Default::default()
    };
    for _ in 0..3 {
        let c0 = evaluate_distributed_dem_contacts(&own0, &ghost0, &mut h0, 0, 1e-4, &hm).unwrap();
        let c1 = evaluate_distributed_dem_contacts(&own1, &ghost1, &mut h1, 1, 1e-4, &hm).unwrap();
        let k0 = evaluate_distributed_dem_bonds(&own0, &ghost0, &mut b0, 0, 1e-4).unwrap();
        let (mut f0, mut f1) = (c0.local.clone(), c1.local.clone());
        accumulate_distributed_dem_forces(&mut f0, &k0.local).unwrap();
        accumulate_distributed_dem_forces(&mut f1, &c0.remote).unwrap();
        accumulate_distributed_dem_forces(&mut f1, &k0.remote).unwrap();
        accumulate_distributed_dem_forces(&mut f0, &c1.remote).unwrap();
        ddv.extend([
            c0.contacts as f64,
            c1.contacts as f64,
            k0.active_bonds as f64,
            c0.remote.len() as f64,
            h0.len() as f64,
        ]);
        integrate_distributed_dem_particles(&mut own0, &f0, v(0.0, -9.81, 0.0), 1e-4).unwrap();
        integrate_distributed_dem_particles(&mut own1, &f1, v(0.0, -9.81, 0.0), 1e-4).unwrap();
        ghost0 = ghosts(&own1, true);
        ghost1 = ghosts(&own0, false);
    }
    for p in &own0 {
        v3(&mut ddv, p.state.position);
        v3(&mut ddv, p.state.linear_velocity);
        v3(&mut ddv, p.state.angular_velocity);
        q4(&mut ddv, p.state.orientation);
    }
    for p in &own1 {
        v3(&mut ddv, p.state.position);
        v3(&mut ddv, p.state.linear_velocity);
        v3(&mut ddv, p.state.angular_velocity);
    }
    for h in &h0 {
        ddv.extend([h.particle_a as f64, h.particle_b as f64]);
        v3(&mut ddv, h.tangential_displacement);
        ddv.push(h.age as f64);
    }
    for bo in &b0 {
        ddv.push(bo.rest_length);
        v3(&mut ddv, bo.tangential_displacement);
        ddv.push(bo.damage);
    }
    close("distributed", &ddv, 1e-10, 1e-13);
}
