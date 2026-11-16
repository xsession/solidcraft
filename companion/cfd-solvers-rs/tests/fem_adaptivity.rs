//! Analytic checks for the P1 Poisson estimator, Dorfler marking and conforming refinement.
use cfd_solvers::fem::*;

fn area(m: &Mesh2D) -> f64 {
    m.triangles
        .iter()
        .map(|t| {
            let (a, b, c) = (m.nodes[t.node[0]], m.nodes[t.node[1]], m.nodes[t.node[2]]);
            0.5 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x))
        })
        .sum()
}

#[test]
fn linear_solution_with_zero_source_has_zero_indicator() {
    let m = make_rectangle_tri_mesh(5, 4, 1.0, 1.0).unwrap();
    let u: Vec<f64> = m
        .nodes
        .iter()
        .map(|p| 2.0 * p.x - 3.0 * p.y + 1.0)
        .collect();
    let e = estimate_poisson_error_tri3(&m, &u, &|_| 0.0, 1.0).unwrap();
    assert!(e.global_indicator < 1e-12);
    assert!(estimate_poisson_error_tri3(&m, &u[1..], &|_| 0.0, 1.0).is_err());
    assert!(estimate_poisson_error_tri3(&m, &u, &|_| 0.0, 0.0).is_err());
    assert!(estimate_poisson_error_tri3(&m, &u, &|_| f64::NAN, 1.0).is_err());
}

#[test]
fn refinement_is_conforming_counter_clockwise_and_area_preserving() {
    let mut m = make_rectangle_tri_mesh(3, 3, 2.0, 1.0).unwrap();
    for round in 0..4 {
        let marks: Vec<bool> = (0..m.element_count())
            .map(|i| (i + round) % 2 == 0)
            .collect();
        let r = refine_tri3_longest_edges(&m, &marks).unwrap();
        assert!(r.element_count() > m.element_count());
        assert!((area(&r) - 2.0).abs() < 1e-12);
        // Every triangle is CCW and every interior edge is shared by exactly two triangles.
        let mut count = std::collections::BTreeMap::new();
        for t in &r.triangles {
            let (a, b, c) = (r.nodes[t.node[0]], r.nodes[t.node[1]], r.nodes[t.node[2]]);
            assert!((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 0.0);
            for (p, q) in [(0, 1), (1, 2), (2, 0)] {
                let k = (t.node[p].min(t.node[q]), t.node[p].max(t.node[q]));
                *count.entry(k).or_insert(0) += 1;
            }
        }
        let boundary = count.values().filter(|c| **c == 1).count();
        assert!(count.values().all(|c| *c <= 2));
        assert_eq!(boundary, r.boundary_edges.len(), "no hanging nodes");
        let perimeter: f64 = r
            .boundary_edges
            .iter()
            .map(|e| {
                (r.nodes[e.node[1]].x - r.nodes[e.node[0]].x)
                    .hypot(r.nodes[e.node[1]].y - r.nodes[e.node[0]].y)
            })
            .sum();
        assert!((perimeter - 6.0).abs() < 1e-12);
        m = r;
    }
    let none = vec![false; m.element_count()];
    assert_eq!(
        refine_tri3_longest_edges(&m, &none)
            .unwrap()
            .element_count(),
        m.element_count()
    );
    assert!(refine_tri3_longest_edges(&m, &none[1..]).is_err());
}

#[test]
fn dorfler_marks_minimal_bulk() {
    let marks = mark_dorfler(&[1.0, 2.0, 3.0, 4.0], 0.6).unwrap();
    // 16 / 30 < 0.6, 25 / 30 >= 0.6.
    assert_eq!(marks, vec![false, false, true, true]);
    assert_eq!(mark_dorfler(&[1.0, 2.0], 1.0).unwrap(), vec![true, true]);
    assert!(mark_dorfler(&[], 0.5).is_err());
    assert!(mark_dorfler(&[1.0, f64::NAN], 0.5).is_err());
    assert!(mark_dorfler(&[1.0], 1.5).is_err());
}

#[test]
fn adaptive_loop_reduces_estimator_on_a_corner_source() {
    // Solve -lap u = f with a localized source, refine adaptively, and check the estimator drops.
    let mut m = make_rectangle_tri_mesh(4, 4, 1.0, 1.0).unwrap();
    let f = |p: Node2| if p.x < 0.25 && p.y < 0.25 { 50.0 } else { 0.0 };
    let mut last = f64::INFINITY;
    for _ in 0..4 {
        let bc: std::collections::BTreeMap<i32, ScalarBoundaryCondition2D> = (0..4)
            .map(|p| (p, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 0.0))))
            .collect();
        let sol = scalar_diffusion2d(&m, |_| 1.0, |_| 0.0, f, &bc, Default::default()).unwrap();
        let est = estimate_poisson_error_tri3(&m, &sol.values, &f, 1.0).unwrap();
        assert!(est.global_indicator < last * 1.05);
        last = est.global_indicator;
        let marks = mark_dorfler(&est.element_indicator, 0.5).unwrap();
        m = refine_tri3_longest_edges(&m, &marks).unwrap();
    }
}
