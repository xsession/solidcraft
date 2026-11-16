use cfd_solvers::{core::*, Result};
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() < t * (1.0 + b.abs()), "{a} != {b}");
}
#[test]
fn forward_and_reverse_ad_match_analytical_derivatives() {
    let x = Dual::<2>::variable(0, 0.7).unwrap();
    let y = Dual::<2>::variable(1, 1.2).unwrap();
    let f = (x * y).sin() + x * x * y.exp();
    let dx = 1.2 * (0.84_f64).cos() + 1.4 * 1.2_f64.exp();
    let dy = 0.7 * 0.84_f64.cos() + 0.49 * 1.2_f64.exp();
    near(f.derivative[0], dx, 1e-12);
    near(f.derivative[1], dy, 1e-12);
    let tape = ReverseTape::new();
    let x = tape.variable(0.7);
    let y = tape.variable(1.2);
    let r = (x.clone() * y.clone()).sin() + x.clone() * x.clone() * y.exp();
    let g = tape.backward(&r, &[x, y]).unwrap();
    near(g[0], dx, 1e-12);
    near(g[1], dy, 1e-12);
    assert!(Dual::<1>::constant(-1.0).ln().is_err());
}
#[test]
fn matrix_free_jacobian_and_adjoint_are_consistent() {
    let residual = |x: &[f64], out: &mut [f64]| -> Result<()> {
        out[0] = 2.0 * x[0] + 3.0 * x[1];
        out[1] = -x[0] + 4.0 * x[1];
        Ok(())
    };
    let x = [0.3, 0.7];
    let mut base = [0.0; 2];
    residual(&x, &mut base).unwrap();
    let u = [0.2, -0.8];
    let v = [0.9, 0.4];
    let mut jv = [0.0; 2];
    let mut jtu = [0.0; 2];
    finite_difference_jacobian_action(&residual, &x, &base, &v, &mut jv, 1e-7).unwrap();
    finite_difference_transpose_jacobian_action(&residual, &x, &base, &u, &mut jtu, 1e-7).unwrap();
    near(dot(&u, &jv), dot(&jtu, &v), 1e-7);
    let system = MatrixFreeNonlinearSystem {
        size: 2,
        residual: &residual,
        jacobian_action: None,
        transpose_jacobian_action: None,
    };
    let mut adjoint = [0.0; 2];
    assert!(
        solve_matrix_free_adjoint(
            &x,
            &system,
            &[2.0, 3.0],
            &mut adjoint,
            MatrixFreeNewtonConfig::default()
        )
        .unwrap()
        .converged
    );
    near(adjoint[0], 1.0, 1e-6);
    near(adjoint[1], 0.0, 1e-6);
}
#[test]
fn jfnk_converges_without_forming_a_matrix() {
    let residual = |x: &[f64], r: &mut [f64]| -> Result<()> {
        r[0] = x[0] * x[0] - 2.0;
        r[1] = x[1] * x[1] - 3.0;
        Ok(())
    };
    let system = MatrixFreeNonlinearSystem {
        size: 2,
        residual: &residual,
        jacobian_action: None,
        transpose_jacobian_action: None,
    };
    let mut x = [1.0, 1.0];
    let result =
        solve_matrix_free_newton(&mut x, &system, MatrixFreeNewtonConfig::default()).unwrap();
    assert!(result.converged);
    near(x[0], 2.0_f64.sqrt(), 1e-7);
    near(x[1], 3.0_f64.sqrt(), 1e-7);
}
#[test]
fn dae_implicit_euler_retains_events_and_known_decay() {
    let residual = |_: f64, x: &[f64], dx: &[f64], r: &mut [f64]| -> Result<()> {
        r[0] = dx[0] + x[0];
        Ok(())
    };
    let event = |_: f64, x: &[f64]| x[0] - 0.5;
    let system = DaeSystem {
        size: 1,
        residual: &residual,
    };
    let result = integrate_implicit_euler(
        &system,
        0.0,
        0.1,
        &[1.0],
        10,
        &[DaeEvent {
            name: "half".into(),
            value: &event,
        }],
        MatrixFreeNewtonConfig::default(),
    )
    .unwrap();
    assert!(result.converged);
    assert_eq!(result.times.len(), 11);
    near(result.states[10][0], 1.1_f64.powi(-10), 1e-7);
    assert_eq!(result.events, vec!["half"]);
}
#[test]
fn generalized_eigenpairs_have_correct_values_and_mass_orthogonality() {
    let k = CsrMatrix::from_dense(&[2.0, 0.0, 0.0, 8.0], 2, 2, 0.0).unwrap();
    let m = CsrMatrix::from_dense(&[2.0, 0.0, 0.0, 2.0], 2, 2, 0.0).unwrap();
    let modes = lowest_generalized_eigenpairs(
        &k,
        &m,
        GeneralizedEigenConfig {
            mode_count: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(modes.len(), 2);
    near(modes[0].eigenvalue, 1.0, 1e-8);
    near(modes[1].eigenvalue, 4.0, 1e-8);
    near(
        csr_bilinear(&m, &modes[0].vector, &modes[1].vector).unwrap(),
        0.0,
        1e-9,
    );
}
