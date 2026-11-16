use cfd_solvers::{core::*, Complex64, Result};
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * (1.0 + b.abs()), "{a} != {b}");
}
fn matrix() -> CsrMatrix {
    CsrMatrix::from_dense(&[4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0], 3, 3, 0.0).unwrap()
}
#[test]
fn sparse_duplicate_assembly_and_transpose() {
    let mut b = CsrBuilder::new(2, 3);
    b.add(0, 2, 3.0).unwrap();
    b.add(0, 2, -1.0).unwrap();
    b.add(1, 1, 4.0).unwrap();
    let a = b.build(0.0).unwrap();
    assert_eq!(a.nonzeros(), 2);
    assert_eq!(a.apply(&[1.0, 2.0, 3.0]).unwrap(), vec![6.0, 8.0]);
    assert_eq!(
        a.transpose().unwrap().apply(&[3.0, 4.0]).unwrap(),
        vec![0.0, 16.0, 6.0]
    );
}
#[test]
fn malformed_sparse_inputs_are_rejected_without_panics() {
    assert!(CsrMatrix::new(2, 2, vec![0, 10, 1], vec![0], vec![1.0]).is_err());
    assert!(CsrMatrix::new(1, 1, vec![0, 2], vec![0, 0], vec![1.0, 2.0]).is_err());
    let j = r#"{"rows":2,"cols":2,"row_offsets":[0,10,1],"column_indices":[0],"values":[1]}"#;
    assert!(serde_json::from_str::<CsrMatrix>(j).is_err());
}
#[test]
fn pcg_solves_spd_system_with_actual_residual() {
    let a = matrix();
    let b = vec![6.0, 10.0, 8.0];
    let mut x = vec![0.0; 3];
    let r = conjugate_gradient(&a, &b, &mut x, LinearConfig::default()).unwrap();
    assert!(r.converged);
    for (i, &v) in x.iter().enumerate() {
        near(v, (i + 1) as f64, 1e-10);
    }
    let residual: Vec<_> = a
        .apply(&x)
        .unwrap()
        .iter()
        .zip(&b)
        .map(|(v, b)| v - b)
        .collect();
    assert!(rms(&residual) < 1e-9);
}
#[test]
fn gmres_and_bicgstab_solve_nonsymmetric_system() {
    let a =
        CsrMatrix::from_dense(&[4.0, 1.0, 0.0, -2.0, 5.0, 1.0, 1.0, 0.0, 3.0], 3, 3, 0.0).unwrap();
    let b = a.apply(&[1.0, -2.0, 3.0]).unwrap();
    for method in 0..2 {
        let mut x = vec![0.0; 3];
        let r = if method == 0 {
            restarted_gmres(
                &b,
                &mut x,
                |v, y| a.multiply(v, y),
                &IdentityPreconditioner,
                LinearConfig::default(),
            )
        } else {
            bicgstab(
                &b,
                &mut x,
                |v, y| a.multiply(v, y),
                &IdentityPreconditioner,
                LinearConfig::default(),
            )
        }
        .unwrap();
        assert!(r.converged);
        for (a, b) in x.iter().zip([1.0, -2.0, 3.0]) {
            near(*a, b, 1e-9);
        }
    }
}
#[test]
fn ilu0_preserves_tridiagonal_factor_solve() {
    let a = matrix();
    let ilu = Ilu0Preconditioner::new(&a, 1e-12).unwrap();
    let mut x = vec![0.0; 3];
    ilu.apply(&[6.0, 10.0, 8.0], &mut x).unwrap();
    for (i, v) in x.iter().enumerate() {
        near(*v, (i + 1) as f64, 1e-12);
    }
}
#[test]
fn failure_does_not_report_false_convergence() {
    let mut x = vec![0.0, 0.0];
    let a = CsrMatrix::from_dense(&[-1.0, 0.0, 0.0, -1.0], 2, 2, 0.0).unwrap();
    let r = preconditioned_conjugate_gradient(
        &[1.0, 1.0],
        &mut x,
        |v, y| a.multiply(v, y),
        &IdentityPreconditioner,
        LinearConfig::default(),
    )
    .unwrap();
    assert!(!r.converged);
}
#[test]
fn dense_real_and_complex_solutions_match_known_vectors() {
    let a = matrix();
    let x = dense_solve(&a.to_dense(), &[6.0, 10.0, 8.0]).unwrap();
    for (i, v) in x.iter().enumerate() {
        near(*v, (i + 1) as f64, 1e-12);
    }
    let z = Complex64::new;
    let a = [z(2.0, 1.0), z(1.0, 0.0), z(0.0, -1.0), z(3.0, 0.5)];
    let expected = [z(1.0, -2.0), z(0.5, 1.0)];
    let b = [
        a[0] * expected[0] + a[1] * expected[1],
        a[2] * expected[0] + a[3] * expected[1],
    ];
    let x = complex_dense_solve(&a, &b).unwrap();
    for i in 0..2 {
        near((x[i] - expected[i]).norm(), 0.0, 1e-12);
    }
    assert!(dense_solve(&[1.0, 2.0, 2.0, 4.0], &[1.0, 2.0]).is_err());
}
#[test]
fn newton_solves_nonlinear_system() {
    let mut x = [1.0, 1.0];
    let residual = |v: &[f64], r: &mut [f64]| -> Result<()> {
        r[0] = v[0] * v[0] - 2.0;
        r[1] = v[1] * v[1] - 3.0;
        Ok(())
    };
    let result = finite_difference_newton(&mut x, residual, NewtonConfig::default()).unwrap();
    assert!(result.converged);
    near(x[0], 2.0_f64.sqrt(), 1e-9);
    near(x[1], 3.0_f64.sqrt(), 1e-9);
}
#[test]
fn fft_inverse_and_spectral_derivative_match_analytical_modes() {
    let n = 64;
    let dx = 2.0 * std::f64::consts::PI / n as f64;
    let f: Vec<_> = (0..n).map(|i| (3.0 * i as f64 * dx).sin()).collect();
    let d = derivative_periodic_1d(&f, dx).unwrap();
    let lap = fractional_laplacian_periodic_1d(&f, dx, 1.0).unwrap();
    for i in 0..n {
        near(d[i], 3.0 * (3.0 * i as f64 * dx).cos(), 1e-10);
        near(lap[i], 9.0 * f[i], 1e-10);
    }
    let mut z: Vec<_> = f.iter().map(|&x| Complex64::new(x, 0.0)).collect();
    fft(&mut z, false).unwrap();
    fft(&mut z, true).unwrap();
    for i in 0..n {
        near(z[i].re, f[i], 1e-12);
    }
    assert!(fft(&mut [z[0]; 3], false).is_err());
}
#[test]
fn quaternion_rotation_and_normalization() {
    let q = Quaternion::axis_angle(Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2).unwrap();
    let p = q.rotate(Vec3::new(1.0, 0.0, 0.0));
    near(p.x, 0.0, 1e-12);
    near(p.y, 1.0, 1e-12);
    near(
        q.integrate(Vec3::new(0.0, 0.0, 2.0), 0.01).unwrap().norm(),
        1.0,
        1e-12,
    );
}
