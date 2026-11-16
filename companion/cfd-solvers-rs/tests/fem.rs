use cfd_solvers::{core::*, fem::*};
fn config() -> LinearConfig {
    LinearConfig {
        max_iterations: 4000,
        relative_tolerance: 1e-12,
        ..Default::default()
    }
}
#[test]
fn shapes_form_a_partition_of_unity() {
    for t in [
        ElementType::Line2,
        ElementType::Tri3,
        ElementType::Quad4,
        ElementType::Tet4,
        ElementType::Hex8,
        ElementType::Prism6,
        ElementType::Pyramid5,
    ] {
        let s = evaluate_shape(t, ReferencePoint::new(0.17, 0.21, 0.29));
        assert!((s.value.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        let sum = s
            .gradient_reference
            .iter()
            .copied()
            .fold(Vec3::ZERO, |a, b| a + b);
        assert!(sum.norm() < 1e-14);
    }
}
#[test]
fn quadrature_integrates_reference_volumes() {
    for (t, volume) in [
        (ElementType::Line2, 2.0),
        (ElementType::Tri3, 0.5),
        (ElementType::Quad4, 4.0),
        (ElementType::Tet4, 1.0 / 6.0),
        (ElementType::Hex8, 8.0),
        (ElementType::Prism6, 1.0),
        (ElementType::Pyramid5, 4.0 / 3.0),
    ] {
        for order in 1..=2 {
            let sum: f64 = gaussian_quadrature(t, order)
                .unwrap()
                .iter()
                .map(|q| q.weight)
                .sum();
            assert!((sum - volume).abs() < 1e-14);
        }
    }
}
#[test]
fn poisson1d_matches_exact_constant_load() {
    let mut p = Poisson1D::new(32, 2.0).unwrap();
    p.solve(|_| 3.0, 1.0, 2.0).unwrap();
    for (x, u) in p.x.iter().zip(&p.u) {
        assert!((u - (1.0 + 0.5 * x + 1.5 * x * (2.0 - x))).abs() < 1e-12);
    }
}
#[test]
fn poisson_tri_and_tet_pass_affine_patch_tests() {
    let mesh = make_rectangle_tri_mesh(7, 5, 2.0, 1.0).unwrap();
    let exact = |p: Node2| 2.0 + 3.0 * p.x - 0.7 * p.y;
    let u = poisson2d(&mesh, 1.4, |_| 0.0, exact, config()).unwrap();
    for (p, v) in mesh.nodes.iter().zip(&u.values) {
        assert!((v - exact(*p)).abs() < 1e-10);
    }
    let mesh = make_box_tet_mesh(3, 3, 3, 1.0, 1.0, 1.0).unwrap();
    let exact = |p: Point3| 1.0 + 2.0 * p.x - 0.3 * p.y + 0.4 * p.z;
    let u = poisson3d(&mesh, 1.2, |_| 0.0, exact, config()).unwrap();
    for (p, v) in mesh.nodes.iter().zip(&u.values) {
        assert!((v - exact(*p)).abs() < 1e-10);
    }
}
#[test]
fn matrix_free_matches_csr_and_is_symmetric() {
    let mesh = make_rectangle_tri_mesh(4, 3, 2.0, 1.0).unwrap();
    let k = |p: Node2| 1.0 + p.x * p.x;
    let a = assemble_tri3_laplace_matrix(&mesh, k).unwrap();
    let x: Vec<_> = (0..mesh.node_count()).map(|i| (i as f64).sin()).collect();
    let expected = a.apply(&x).unwrap();
    let mut actual = vec![0.0; x.len()];
    apply_tri3_laplace_matrix_free(&mesh, k, &x, &mut actual).unwrap();
    for (a, b) in actual.iter().zip(expected) {
        assert!((a - b).abs() < 1e-12);
    }
}
#[test]
fn mixed_robin_and_neumann_reproduce_linear_field() {
    use std::collections::BTreeMap;
    let mesh = make_rectangle_tri_mesh(7, 5, 1.0, 1.0).unwrap();
    let mut bc = BTreeMap::new();
    bc.insert(0, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 1.0)));
    bc.insert(
        1,
        ScalarBoundaryCondition2D::Robin {
            alpha: Box::new(|_| 2.0),
            value: Box::new(|_| 5.0),
        },
    );
    bc.insert(2, ScalarBoundaryCondition2D::Neumann(Box::new(|_| 0.0)));
    bc.insert(3, ScalarBoundaryCondition2D::Neumann(Box::new(|_| 0.0)));
    let s = scalar_diffusion2d(&mesh, |_| 1.0, |_| 0.0, |_| 0.0, &bc, config()).unwrap();
    for (p, u) in mesh.nodes.iter().zip(s.values) {
        assert!((u - (1.0 + p.x)).abs() < 1e-10);
    }
}
#[test]
fn heat_preserves_uniform_field_and_element_sources() {
    let mesh = make_rectangle_tri_mesh(5, 5, 1.0, 1.0).unwrap();
    let n = mesh.element_count();
    let mut heat = Heat2D::new(mesh, Heat2DConfig::default()).unwrap();
    heat.initialize(|_| 7.0).unwrap();
    heat.step_element_source(&vec![0.0; n], |_, _| 7.0).unwrap();
    assert!(heat.temperature.iter().all(|t| (t - 7.0).abs() < 1e-12));
    assert_eq!(heat.time, 0.001);
}
