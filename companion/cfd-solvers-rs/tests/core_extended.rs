use cfd_solvers::{core::*, Complex64 as C};

#[test]
fn finite_large_residuals_do_not_overflow_to_false_convergence() {
    assert!(rms(&[f64::NAN, 0.0]).is_nan());
    assert!((rms(&[1e300, 1e300]) / 1e300 - 1.0).abs() < 1e-14);
    let mut x = [0.0, 0.0];
    let r = preconditioned_conjugate_gradient(
        &[1e300, 1e300],
        &mut x,
        |v, y| {
            y.copy_from_slice(v);
            Ok(())
        },
        &IdentityPreconditioner,
        LinearConfig::default(),
    )
    .unwrap();
    assert!(!r.converged);
    assert!(r.residual_rms.is_finite());
}

#[test]
fn geometric_multigrid_solves_a_known_discrete_mode() {
    let n = 15;
    let h = 1.0 / (n + 1) as f64;
    let exact: Vec<_> = (0..n)
        .flat_map(|j| {
            (0..n).map(move |i| {
                (std::f64::consts::PI * (i + 1) as f64 * h).sin()
                    * (std::f64::consts::PI * (j + 1) as f64 * h).sin()
            })
        })
        .collect();
    let lambda = 8.0 * (0.5 * std::f64::consts::PI * h).sin().powi(2) / (h * h);
    let rhs: Vec<_> = exact.iter().map(|v| v * lambda).collect();
    let mut x = vec![0.0; n * n];
    let r = PoissonMultigrid2D::new(n)
        .unwrap()
        .solve(&rhs, &mut x, 200, 1e-10)
        .unwrap();
    assert!(
        r.converged,
        "{} cycles, residual {}",
        r.cycles, r.relative_residual
    );
    assert!(x.iter().zip(exact).all(|(a, b)| (a - b).abs() < 1e-8));
    assert!(PoissonMultigrid2D::new(14).is_err());
}
#[test]
fn block_csr_and_block_preconditioner_solve_coupled_system() {
    let mut b = BlockCsrBuilder::<2>::new(3, 3).unwrap();
    for i in 0..3 {
        b.add(i, i, [[4.0, 0.7], [-0.2, 3.0]]).unwrap();
        if i + 1 < 3 {
            b.add(i, i + 1, [[-1.0, 0.1], [0.2, -0.5]]).unwrap();
            b.add(i + 1, i, [[-0.8, 0.0], [0.1, -0.6]]).unwrap();
        }
    }
    b.add(0, 0, [[1.0, 0.0], [0.0, 1.0]]).unwrap();
    let a = b.build().unwrap();
    let exact = [1.0, -2.0, 0.3, 0.7, -1.0, 2.0];
    let rhs = a.apply(&exact).unwrap();
    let scalar = a.scalar_csr().unwrap();
    assert_eq!(rhs, scalar.apply(&exact).unwrap());
    let mut x = [0.0; 6];
    let r = restarted_block_gmres(
        &rhs,
        &mut x,
        |v, out| a.multiply(v, out),
        &a.diagonal_blocks(),
        LinearConfig {
            relative_tolerance: 1e-12,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(r.converged);
    assert!(x.iter().zip(exact).all(|(a, b)| (a - b).abs() < 1e-10));
}
#[test]
fn invalid_block_shapes_and_singular_diagonals_are_rejected() {
    assert!(BlockCsrBuilder::<0>::new(1, 1).is_err());
    assert!(
        BlockCsrMatrix::<2>::new(2, 2, vec![0, 2, 1], vec![0], vec![[[1.0, 0.0], [0.0, 1.0]]])
            .is_err()
    );
    assert!(BlockJacobiPreconditioner::<2>::new(&[[[1.0, 2.0], [2.0, 4.0]]], 1e-14).is_err());
    let p = BlockJacobiPreconditioner::<2>::new(&[[[0.0, 2.0], [1.0, 3.0]]], 1e-14).unwrap();
    let mut out = [0.0; 2];
    p.apply(&[4.0, 7.0], &mut out).unwrap();
    assert!((out[0] - 1.0).abs() < 1e-14 && out[1] == 2.0);
}
#[test]
fn complex_sparse_solves_real_and_imaginary_coupling() {
    let mut a = ComplexCsrBuilder::new(3).unwrap();
    a.add(0, 0, C::new(3.0, 1.0)).unwrap();
    a.add(0, 1, C::new(0.2, -0.3)).unwrap();
    a.add(1, 0, C::new(-0.4, 0.1)).unwrap();
    a.add(1, 1, C::new(2.0, 0.4)).unwrap();
    a.add(1, 2, C::new(0.0, 0.2)).unwrap();
    a.add(2, 1, C::new(0.1, 0.0)).unwrap();
    a.add(2, 2, C::new(4.0, -0.5)).unwrap();
    let exact = [C::new(1.0, 0.2), C::new(-0.3, 0.8), C::new(0.7, -0.1)];
    let rhs = a.apply(&exact).unwrap();
    let r = solve_complex_sparse(
        &a,
        &rhs,
        ComplexSparseSolveConfig {
            relative_tolerance: 1e-12,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(r.linear_result.converged);
    for (x, e) in r.solution.iter().zip(exact) {
        assert!((*x - e).norm() < 1e-10);
    }
    let mut a = ComplexCsrBuilder::new(1).unwrap();
    a.add(0, 0, C::new(0.0, 2.0)).unwrap();
    let r = solve_complex_sparse(&a, &[C::new(2.0, 4.0)], Default::default()).unwrap();
    assert!(r.linear_result.converged);
    assert!((r.solution[0] - C::new(2.0, -1.0)).norm() < 1e-8);
}
#[test]
fn stencil_and_fused_residual_agree_with_assembled_operator() {
    let s = StructuredSevenPointOperator::new(4, 3, 2, 7.0, -0.8).unwrap();
    let mut b = CsrBuilder::new(24, 24);
    for z in 0..2 {
        for y in 0..3 {
            for x in 0..4 {
                let i = (z * 3 + y) * 4 + x;
                b.add(i, i, 7.0).unwrap();
                for (dx, dy, dz) in [
                    (-1, 0, 0),
                    (1, 0, 0),
                    (0, -1, 0),
                    (0, 1, 0),
                    (0, 0, -1),
                    (0, 0, 1),
                ] {
                    let (xx, yy, zz) = (x as i32 + dx, y as i32 + dy, z as i32 + dz);
                    if (0..4).contains(&xx) && (0..3).contains(&yy) && (0..2).contains(&zz) {
                        b.add(i, ((zz * 3 + yy) * 4 + xx) as usize, -0.8).unwrap();
                    }
                }
            }
        }
    }
    let a = b.build(0.0).unwrap();
    let x: Vec<_> = (0..24).map(|i| (i as f64).sin()).collect();
    let mut y = vec![0.0; 24];
    s.apply(&x, &mut y).unwrap();
    let exact = a.apply(&x).unwrap();
    for (x, y) in y.iter().zip(exact) {
        assert!((x - y).abs() < 1e-14);
    }
    let rhs: Vec<_> = y.iter().map(|v| v + 0.2).collect();
    let norm = s.apply_residual_l2(&rhs, &x, &mut y).unwrap();
    assert!((norm - (24.0_f64).sqrt() * 0.2).abs() < 1e-14);
    assert_eq!(s.estimated_apply_bytes().unwrap(), 24 * 64);
    assert!(Grid2D::new(1, 2, 1.0, 1.0).is_err());
    let g = Grid2D::new(4, 3, 0.1, 0.2).unwrap();
    assert_eq!(g.index(3, 2).unwrap(), 11);
    assert!(g.index(4, 2).is_err());
}
