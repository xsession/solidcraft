//! Analytic and invariant checks for the mechanics module.
use cfd_solvers::{core::Vec3, mechanics::*, multibody::Quaternion};

const F: Tensor3 = [1.1, 0.05, -0.02, 0.03, 0.95, 0.04, -0.01, 0.02, 1.07];
const D: Tensor3 = [0.2, -0.1, 0.3, 0.05, 0.4, -0.2, 0.1, 0.0, -0.3];

#[test]
fn neo_hookean_identity_and_tangent() {
    let m = CompressibleNeoHookean3D::new(NeoHookean3DConfig {
        shear_modulus: 2.0,
        lame_lambda: 3.0,
    })
    .unwrap();
    let r = m.evaluate(&IDENTITY3).unwrap();
    assert!(r.first_piola.iter().all(|v| v.abs() < 1e-15) && r.energy_density.abs() < 1e-15);
    let h = 1e-6;
    let (mut fp, mut fm) = (F, F);
    for i in 0..9 {
        fp[i] += h * D[i];
        fm[i] -= h * D[i];
    }
    let (pp, pm) = (
        m.evaluate(&fp).unwrap().first_piola,
        m.evaluate(&fm).unwrap().first_piola,
    );
    let t = m.tangent_action(&F, &D).unwrap();
    for i in 0..9 {
        assert!(((pp[i] - pm[i]) / (2.0 * h) - t[i]).abs() < 1e-7);
    }
    assert!(m
        .evaluate(&[-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
        .is_err());
}

#[test]
fn j2_consistent_tangent_matches_finite_difference() {
    let m = J2Plasticity3D::new(J2Plasticity3DConfig {
        youngs_modulus: 200.0,
        poisson_ratio: 0.28,
        yield_stress: 0.8,
        hardening_modulus: 15.0,
        tangent_step: 1e-7,
    })
    .unwrap();
    let e = [8e-3, -2e-3, -1e-3, 2e-3, -1e-3, 5e-4];
    let s0 = J2Plasticity3DState::default();
    let r = m.update(&e, &s0);
    assert!(r.yielded);
    let h = 1e-8;
    for col in 0..6 {
        let (mut ep, mut em) = (e, e);
        ep[col] += h;
        em[col] -= h;
        let (sp, sm) = (m.update(&ep, &s0).stress, m.update(&em, &s0).stress);
        for row in 0..6 {
            let fd = (sp[row] - sm[row]) / (2.0 * h);
            assert!(
                (fd - r.algorithmic_tangent[row * 6 + col]).abs() < 1e-4 * 200.0,
                "({row},{col})"
            );
        }
    }
    // Consistency: the returned stress lies on the hardened yield surface.
    let mean = (r.stress[0] + r.stress[1] + r.stress[2]) / 3.0;
    let s = [
        r.stress[0] - mean,
        r.stress[1] - mean,
        r.stress[2] - mean,
        r.stress[3],
        r.stress[4],
        r.stress[5],
    ];
    let q = (1.5
        * (s[0] * s[0]
            + s[1] * s[1]
            + s[2] * s[2]
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])))
        .sqrt();
    assert!((q - (0.8 + 15.0 * r.state.equivalent_plastic_strain)).abs() < 1e-10);
}

#[test]
fn finite_strain_j2_objective_and_hydrostatic() {
    let m = FiniteStrainJ2::new(FiniteStrainJ2Config {
        hardening_modulus: 5.0,
        ..Default::default()
    })
    .unwrap();
    let report = calibrate_finite_strain_j2(&m, &F, &D, &FiniteStrainJ2State::default()).unwrap();
    assert!(report.passed(), "{report:?}");
    let r = m
        .evaluate(&IDENTITY3, &FiniteStrainJ2State::default())
        .unwrap();
    assert!(r.cauchy_stress.iter().all(|v| v.abs() < 1e-14) && !r.yielded);
}

#[test]
fn tet4_rigid_translation_and_equilibrium() {
    let m = CompressibleNeoHookean3D::new(NeoHookean3DConfig::default()).unwrap();
    let v = Vec3::new;
    let reference = [
        v(0.0, 0.0, 0.0),
        v(1.0, 0.0, 0.0),
        v(0.0, 1.0, 0.0),
        v(0.0, 0.0, 1.0),
    ];
    let shifted = reference.map(|p| p + v(0.3, -0.2, 0.7));
    let r = assemble_total_lagrangian_tet4(&reference, &shifted, &[0.0; 12], &m).unwrap();
    assert!(
        r.internal_force.iter().all(|f| f.abs() < 1e-14)
            && (r.reference_volume - 1.0 / 6.0).abs() < 1e-15
    );
    let deformed = [
        v(0.0, 0.0, 0.0),
        v(1.2, 0.1, 0.0),
        v(0.0, 0.9, 0.05),
        v(0.02, 0.0, 1.1),
    ];
    let r = assemble_total_lagrangian_tet4(&reference, &deformed, &[0.0; 12], &m).unwrap();
    for c in 0..3 {
        let sum: f64 = (0..4).map(|n| r.internal_force[3 * n + c]).sum();
        assert!(sum.abs() < 1e-13);
    }
    let flipped = [reference[0], reference[2], reference[1], reference[3]];
    assert!(assemble_total_lagrangian_tet4(&flipped, &deformed, &[0.0; 12], &m).is_err());
}

#[test]
fn linear_bar_and_conservative_particles() {
    let cfg = NonlinearBarConfig {
        elements: 5,
        cubic_modulus: 0.0,
        ..Default::default()
    };
    let r = solve_nonlinear_bar(&cfg, 3.0, &[]).unwrap();
    assert!(r.solve.converged);
    for (i, u) in r.displacement.iter().enumerate() {
        assert!((u - 3.0 / 100.0 * (i + 1) as f64 / 5.0).abs() < 1e-12);
    }
    let mut pd = PeridynamicBar1D::new(PeridynamicBarConfig {
        points: 15,
        horizon: 0.2,
        ..Default::default()
    })
    .unwrap();
    pd.set_velocity(3, 0.4).unwrap();
    pd.set_velocity(9, -0.1).unwrap();
    for _ in 0..50 {
        pd.step(0.01, &[]).unwrap();
    }
    let momentum: f64 = pd.velocity().iter().sum();
    assert!((momentum - 0.3).abs() < 1e-12);
    let mut im = PenaltyImpact1D::new(1.0, 2.0, 0.01, 1e4, 0.0).unwrap();
    im.set_velocities(1.0, -0.5).unwrap();
    for _ in 0..500 {
        im.step(1e-4).unwrap();
        let s = im.state();
        assert!((s.left_velocity + 2.0 * s.right_velocity).abs() < 1e-12);
    }
}

#[test]
fn meshfree_graph_matches_brute_force() {
    let pts: Vec<Vec3> = (0..80)
        .map(|i| {
            Vec3::new(
                (0.731 * i as f64).sin(),
                (0.413 * i as f64).cos(),
                0.01 * i as f64,
            )
        })
        .collect();
    let cfg = MeshfreeNeighborGraphConfig {
        support_radius: 0.3,
        ..Default::default()
    };
    let g = build_meshfree_neighbor_graph(&pts, &cfg).unwrap();
    for i in 0..pts.len() {
        let brute: Vec<usize> = (0..pts.len())
            .filter(|&j| j != i && (pts[i] - pts[j]).norm_squared() <= 0.09)
            .collect();
        assert_eq!(g.neighbors(i).unwrap(), brute.as_slice());
    }
    assert!(g.neighbors(80).is_err());
}

#[test]
fn phase_field_irreversibility_and_operator_tangent() {
    let mut pf = PhaseFieldFracture1D::new(PhaseFieldFracture1DConfig {
        cells: 16,
        ..Default::default()
    })
    .unwrap();
    pf.step(0.03).unwrap();
    let loaded = pf.damage().to_vec();
    pf.step(0.0).unwrap();
    assert!(pf.damage().iter().zip(&loaded).all(|(a, b)| a >= b));
    let mut op = CoupledPhaseFieldFracture1D::new(CoupledPhaseFieldFracture1DConfig {
        cells: 6,
        ..Default::default()
    })
    .unwrap();
    let n = op.size();
    let state: Vec<f64> = (0..n)
        .map(|i| {
            if i < 7 {
                0.01 * i as f64
            } else {
                0.1 + 0.05 * (i - 7) as f64
            }
        })
        .collect();
    op.update_history(&state).unwrap();
    let dir: Vec<f64> = (0..n).map(|i| (1.3 * i as f64).cos()).collect();
    let mut t = vec![0.0; n];
    op.tangent_action(&state, &dir, 0.06, &mut t).unwrap();
    let h = 1e-7;
    let (sp, sm): (Vec<f64>, Vec<f64>) = (
        state.iter().zip(&dir).map(|(s, d)| s + h * d).collect(),
        state.iter().zip(&dir).map(|(s, d)| s - h * d).collect(),
    );
    let (mut rp, mut rm) = (vec![0.0; n], vec![0.0; n]);
    op.evaluate(&sp, 0.06, &mut rp).unwrap();
    op.evaluate(&sm, 0.06, &mut rm).unwrap();
    for i in 0..n {
        assert!(
            ((rp[i] - rm[i]) / (2.0 * h) - t[i]).abs() < 1e-5 * (1.0 + t[i].abs()),
            "{i}"
        );
    }
}

#[test]
fn cosserat_rest_state_is_stationary() {
    let mut rod = StableCosseratRod::new(CosseratRodConfig {
        gravity: Vec3::ZERO,
        ..Default::default()
    })
    .unwrap();
    let pos: Vec<Vec3> = (0..5)
        .map(|i| Vec3::new(0.0, 0.0, 0.2 * i as f64))
        .collect();
    rod.initialize(&pos, &[]).unwrap();
    rod.run(20).unwrap();
    let d = *rod.diagnostics();
    assert!(d.finite && d.kinetic_energy < 1e-20 && d.max_stretch_shear_error < 1e-12);
    for (n, p) in rod.nodes().iter().zip(&pos) {
        assert!((n.position - *p).norm() < 1e-12 && n.orientation == Quaternion::default());
    }
    assert!(rod.set_fixed(9, true).is_err());
}
