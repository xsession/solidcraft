//! Value-by-value comparison with the C++ upstream mechanics/*. Goldens come from
//! `parity/oracle_mechanics.cpp`.
use cfd_solvers::{core::Vec3, fem::Node2, mechanics::*, multibody};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-mechanics.json")).unwrap()
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

fn b(x: bool) -> f64 {
    if x {
        1.0
    } else {
        0.0
    }
}

const F: Tensor3 = [1.1, 0.05, -0.02, 0.03, 0.95, 0.04, -0.01, 0.02, 1.07];
const D: Tensor3 = [0.2, -0.1, 0.3, 0.05, 0.4, -0.2, 0.1, 0.0, -0.3];

fn fsj2() -> FiniteStrainJ2 {
    FiniteStrainJ2::new(FiniteStrainJ2Config {
        youngs_modulus: 200.0,
        poisson_ratio: 0.3,
        yield_stress: 1.0,
        hardening_modulus: 5.0,
        tangent_step: 1e-7,
        spectral_iterations: 32,
    })
    .unwrap()
}

#[test]
fn constitutive_updates() {
    let nh = CompressibleNeoHookean3D::new(NeoHookean3DConfig {
        shear_modulus: 2.0,
        lame_lambda: 3.0,
    })
    .unwrap();
    let r = nh.evaluate(&F).unwrap();
    let mut v = vec![r.determinant, r.energy_density];
    v.extend(r.first_piola);
    v.extend(nh.tangent_action(&F, &D).unwrap());
    close("neo", &v, 1e-13, 1e-15);

    let j2 = J2Plasticity3D::new(J2Plasticity3DConfig {
        youngs_modulus: 200.0,
        poisson_ratio: 0.28,
        yield_stress: 0.8,
        hardening_modulus: 15.0,
        tangent_step: 1e-7,
    })
    .unwrap();
    let mut st = J2Plasticity3DState::default();
    let mut jv = vec![];
    for e in [
        [1e-3, -2e-4, 0.0, 1e-4, 0.0, 0.0],
        [8e-3, -2e-3, -1e-3, 2e-3, -1e-3, 5e-4],
        [1.2e-2, -3e-3, -2e-3, 3e-3, -2e-3, 1e-3],
    ] {
        let r = j2.update(&e, &st);
        st = r.state;
        jv.extend(r.stress);
        jv.extend(r.algorithmic_tangent);
        jv.extend(r.state.plastic_strain);
        jv.extend([
            r.state.equivalent_plastic_strain,
            r.trial_equivalent_stress,
            r.yield_function,
            b(r.yielded),
        ]);
    }
    close("j2", &jv, 1e-12, 1e-15);

    let fj = fsj2();
    let mut fs = FiniteStrainJ2State::default();
    let f2: Tensor3 = [1.03, 0.04, 0.0, 0.01, 0.99, 0.02, 0.0, -0.01, 1.0];
    let mut fv = vec![];
    for g in [f2, F] {
        let r = fj.evaluate(&g, &fs).unwrap();
        fs = r.state;
        fv.extend([
            r.determinant,
            r.energy_density,
            r.trial_equivalent_stress,
            r.yield_function,
            b(r.yielded),
            r.state.equivalent_plastic_strain,
        ]);
        fv.extend(r.cauchy_stress);
        fv.extend(r.first_piola);
        fv.extend(r.state.plastic_deformation_gradient);
    }
    fv.extend(fj.tangent_action(&F, &D, &fs).unwrap());
    fv.extend(fj.central_tangent_action(&F, &D, &fs).unwrap());
    close("fsj2", &fv, 1e-6, 1e-9);
    let cal = calibrate_finite_strain_j2(&fj, &F, &D, &FiniteStrainJ2State::default()).unwrap();
    assert!(cal.passed());
    close(
        "calib",
        &[
            cal.determinant,
            cal.rotation_relative_error,
            cal.hydrostatic_deviatoric_ratio,
            cal.tangent_step_relative_change,
            b(cal.positive_state),
            b(cal.objective_rotation),
            b(cal.hydrostatic_response),
            b(cal.tangent_step_converged),
        ],
        1e-3,
        1e-12,
    );
}

#[test]
fn finite_strain_elements() {
    let nh = CompressibleNeoHookean3D::new(NeoHookean3DConfig {
        shear_modulus: 2.0,
        lame_lambda: 3.0,
    })
    .unwrap();
    let v = Vec3::new;
    let reference = [
        v(0.0, 0.0, 0.0),
        v(1.0, 0.0, 0.1),
        v(0.1, 1.0, 0.0),
        v(0.0, 0.2, 1.0),
    ];
    let current = [
        v(0.01, 0.0, 0.0),
        v(1.05, 0.02, 0.1),
        v(0.12, 0.98, 0.01),
        v(0.0, 0.21, 1.04),
    ];
    let dir = [
        0.1, 0.0, 0.0, -0.1, 0.2, 0.0, 0.05, 0.0, 0.3, 0.0, -0.2, 0.1,
    ];
    let t4 = assemble_total_lagrangian_tet4(&reference, &current, &dir, &nh).unwrap();
    let mut out = vec![t4.reference_volume, t4.determinant, t4.strain_energy];
    out.extend(t4.internal_force);
    out.extend(t4.tangent_action);
    let current2 = [
        v(0.0, 0.0, 0.0),
        v(1.08, 0.02, 0.1),
        v(0.12, 0.96, 0.01),
        v(0.0, 0.21, 1.02),
    ];
    let fj = fsj2();
    let p = assemble_total_lagrangian_plastic_tet4(
        &reference,
        &current2,
        &dir,
        &fj,
        &FiniteStrainJ2State::default(),
    )
    .unwrap();
    out.extend([
        p.reference_volume,
        p.determinant,
        p.strain_energy,
        p.trial_equivalent_stress,
        p.yield_function,
        b(p.yielded),
    ]);
    out.extend(p.internal_force);
    out.extend(p.tangent_action);
    let n = |x, y| Node2 { x, y };
    let t3 = assemble_total_lagrangian_plastic_tri3(
        &[n(0.0, 0.0), n(1.0, 0.1), n(0.2, 0.9)],
        &[n(0.0, 0.0), n(1.04, 0.12), n(0.21, 0.88)],
        &[0.1, 0.0, -0.1, 0.2, 0.05, 0.3],
        &fj,
        0.5,
        &FiniteStrainJ2State::default(),
    )
    .unwrap();
    out.extend([
        t3.reference_area,
        t3.determinant,
        t3.strain_energy,
        t3.trial_equivalent_stress,
        t3.yield_function,
        b(t3.yielded),
    ]);
    out.extend(t3.internal_force);
    out.extend(t3.tangent_action);
    out.extend(t3.state.plastic_deformation_gradient);
    close("elements", &out, 1e-6, 1e-9);
}

#[test]
fn continuum_meshfree_and_fracture() {
    let bar = solve_nonlinear_bar(
        &NonlinearBarConfig {
            elements: 6,
            ..Default::default()
        },
        5.0,
        &[],
    )
    .unwrap();
    let mut bv = bar.displacement.clone();
    bv.extend([
        bar.solve.iterations as f64,
        b(bar.solve.converged),
        bar.solve.residual_rms,
    ]);
    close("bar", &bv, 1e-10, 1e-12);

    let mut pd = PeridynamicBar1D::new(PeridynamicBarConfig {
        points: 21,
        horizon: 0.16,
        ..Default::default()
    })
    .unwrap();
    pd.set_velocity(20, 0.8).unwrap();
    pd.set_displacement(0, -0.01).unwrap();
    let mut ext = vec![0.0; 21];
    ext[10] = 0.5;
    let mut pv = vec![];
    for i in 0..40 {
        let d = pd.step(0.01, if i % 2 == 1 { &[] } else { &ext }).unwrap();
        if i % 8 == 7 {
            pv.extend([d.broken_bonds as f64, d.kinetic_energy, d.potential_energy]);
        }
    }
    pv.extend(pd.displacement());
    pv.extend(pd.velocity());
    close("peri", &pv, 1e-12, 1e-15);

    let mut im = PenaltyImpact1D::new(1.0, 2.0, 0.05, 1e4, 3.0).unwrap();
    im.set_velocities(1.0, -0.5).unwrap();
    let mut iv = vec![];
    for i in 0..120 {
        im.step(1e-3).unwrap();
        if i % 20 == 19 {
            let s = im.state();
            iv.extend([
                s.left_position,
                s.right_position,
                s.left_velocity,
                s.right_velocity,
                s.contact_force,
            ]);
        }
    }
    close("impact", &iv, 1e-12, 1e-14);

    let pts: Vec<Vec3> = (0..60)
        .map(|i| {
            Vec3::new(
                (0.37 * i as f64) % 1.0,
                (0.61 * i as f64) % 0.8,
                (0.13 * i as f64) % 0.3,
            )
        })
        .collect();
    let mut mc = MeshfreeNeighborGraphConfig {
        support_radius: 0.22,
        ..Default::default()
    };
    let g = build_meshfree_neighbor_graph(&pts, &mc).unwrap();
    let mut gv = vec![
        g.cell_counts[0] as f64,
        g.cell_counts[1] as f64,
        g.cell_counts[2] as f64,
        g.cell_size,
    ];
    gv.extend(g.offsets.iter().map(|v| *v as f64));
    gv.extend(g.indices.iter().map(|v| *v as f64));
    let owner: Vec<i32> = (0..60).map(|i| (i * 7) % 3).collect();
    let h = classify_meshfree_halo(&g, &owner, 1).unwrap();
    gv.extend(h.owned.iter().map(|v| *v as f64));
    gv.push(-1.0);
    gv.extend(h.halo.iter().map(|v| *v as f64));
    mc.include_self = true;
    gv.push(
        build_meshfree_neighbor_graph(&pts, &mc)
            .unwrap()
            .indices
            .len() as f64,
    );
    close("meshfree", &gv, 0.0, 0.0);

    let mut pf = PhaseFieldFracture1D::new(PhaseFieldFracture1DConfig {
        cells: 24,
        notch_cell: Some(12),
        initial_notch_damage: 0.3,
        damage_iterations: 300,
        ..Default::default()
    })
    .unwrap();
    let mut ff = vec![];
    for u in [0.005, 0.01, 0.02, 0.015] {
        let d = pf.step(u).unwrap();
        ff.extend([
            d.reaction_force,
            d.elastic_energy,
            d.fracture_energy,
            d.total_energy,
            d.maximum_damage,
        ]);
    }
    ff.extend(pf.damage());
    ff.extend(pf.history());
    ff.push(pf.diagnostics().total_energy);
    close("pf", &ff, 1e-12, 1e-16);

    let mut op = CoupledPhaseFieldFracture1D::new(CoupledPhaseFieldFracture1DConfig {
        cells: 8,
        ..Default::default()
    })
    .unwrap();
    let n = op.size();
    let state: Vec<f64> = (0..n)
        .map(|i| {
            if i < 9 {
                0.002 * (i * i) as f64
            } else {
                0.05 * ((i * 3) % 7) as f64
            }
        })
        .collect();
    let dir: Vec<f64> = (0..n).map(|i| (0.7 * i as f64).sin()).collect();
    let (mut res, mut out) = (vec![0.0; n], vec![0.0; n]);
    op.update_history(&state).unwrap();
    let er = op.evaluate(&state, 0.18, &mut res).unwrap();
    op.tangent_action(&state, &dir, 0.18, &mut out).unwrap();
    let mut ov = vec![
        er.residual_l2,
        er.reaction_force,
        er.elastic_energy,
        er.fracture_energy,
        er.total_energy,
        er.maximum_damage,
    ];
    ov.extend(&res);
    ov.extend(&out);
    ov.extend(op.history());
    close("pfop", &ov, 1e-12, 1e-14);
}

#[test]
fn cosserat_rod() {
    let mut rod = StableCosseratRod::new(CosseratRodConfig {
        dt: 2e-3,
        ..Default::default()
    })
    .unwrap();
    let rp: Vec<Vec3> = (0..6)
        .map(|i| Vec3::new(0.0, 0.0, 0.1 * i as f64))
        .collect();
    let mut ro = vec![multibody::Quaternion::default(); 6];
    ro[3] = multibody::quaternion_from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.1).unwrap();
    rod.initialize(&rp, &ro).unwrap();
    rod.set_fixed(0, true).unwrap();
    let mut cv = vec![];
    for s in 0..25 {
        rod.add_force(5, Vec3::new(0.5, 0.0, 0.0)).unwrap();
        rod.add_torque(4, Vec3::new(0.0, 0.0, 0.2)).unwrap();
        rod.step().unwrap();
        if s % 5 == 4 {
            let d = rod.diagnostics();
            cv.extend([
                d.max_stretch_shear_error,
                d.max_bend_twist_error,
                d.kinetic_energy,
                d.elastic_energy,
                d.linear_momentum_norm,
                rod.checksum(),
            ]);
        }
    }
    for n in rod.nodes() {
        cv.extend([
            n.position.x,
            n.position.y,
            n.position.z,
            n.orientation.w,
            n.orientation.x,
            n.orientation.y,
            n.orientation.z,
            n.angular_velocity.x,
            n.angular_velocity.y,
            n.angular_velocity.z,
        ]);
    }
    cv.push(rod.time());
    close("cosserat", &cv, 1e-9, 1e-12);
}
