use cfd_solvers::{core::LinearConfig, fem::*};
#[test]
fn elasticity_recovers_affine_displacement_in_both_modes() {
    for mode in [ElasticityMode2D::PlaneStress, ElasticityMode2D::PlaneStrain] {
        let mesh = make_rectangle_tri_mesh(6, 5, 1.0, 1.0).unwrap();
        let mut s = Elasticity2D::new(
            mesh,
            Elasticity2DConfig {
                youngs_modulus: 3.0,
                mode,
                linear: LinearConfig {
                    relative_tolerance: 1e-12,
                    ..Default::default()
                },
                ..Default::default()
            },
        )
        .unwrap();
        let exact = |p: Node2| Displacement2 {
            x: 0.01 * p.x + 0.003 * p.y,
            y: -0.004 * p.x + 0.002 * p.y,
        };
        s.set_dirichlet(
            |p| p.x == 0.0 || p.x == 1.0 || p.y == 0.0 || p.y == 1.0,
            exact,
            true,
            true,
        )
        .unwrap();
        s.solve(|_| Displacement2::default()).unwrap();
        for (p, u) in s.mesh().nodes.iter().zip(s.displacement()) {
            let e = exact(*p);
            assert!((u.x - e.x).abs() < 1e-10);
            assert!((u.y - e.y).abs() < 1e-10);
        }
    }
}
#[test]
fn free_thermal_expansion_and_fixed_free_bar_modes() {
    let mesh = make_rectangle_tri_mesh(5, 4, 1.0, 1.0).unwrap();
    let n = mesh.node_count();
    let mut s = Elasticity2D::new(
        mesh,
        Elasticity2DConfig {
            youngs_modulus: 10.0,
            ..Default::default()
        },
    )
    .unwrap();
    s.set_dirichlet(
        |p| p.x == 0.0 && p.y == 0.0,
        |_| Displacement2::default(),
        true,
        true,
    )
    .unwrap();
    s.set_dirichlet(
        |p| p.x == 1.0 && p.y == 0.0,
        |_| Displacement2::default(),
        false,
        true,
    )
    .unwrap();
    s.solve_thermal(&vec![310.0; n], 1e-5, 300.0, |_| Displacement2::default())
        .unwrap();
    for (p, u) in s.mesh().nodes.iter().zip(s.displacement()) {
        assert!((u.x - 1e-4 * p.x).abs() < 1e-8);
        assert!((u.y - 1e-4 * p.y).abs() < 1e-8);
    }
    let config = BarModal1DConfig {
        elements: 40,
        young_modulus: 100.0,
        density: 4.0,
        ..Default::default()
    };
    let modes = bar_modal1d(config, 2).unwrap();
    for (i, m) in modes.iter().enumerate() {
        let exact = (2 * i + 1) as f64 * std::f64::consts::PI / 2.0 * 5.0;
        assert!((m.angular_frequency / exact - 1.0).abs() < 0.001);
        assert!(m.eigenpair.converged);
        assert_eq!(m.displacement[0], 0.0);
    }
}
#[test]
fn piezo_bar_recovers_uniform_coupled_solution() {
    let c = PiezoelectricBar1DConfig::default();
    let mut s = PiezoelectricBar1D::new(c).unwrap();
    s.solve_receive(1.0).unwrap();
    let effective =
        c.elastic_modulus_pa + c.piezoelectric_stress_c_per_m2.powi(2) / c.permittivity_f_per_m;
    let strain = 1.0 / (c.area_m2 * effective);
    let expected = strain * c.length_m;
    assert!((s.tip_displacement_m() / expected - 1.0).abs() < 1e-10);
    assert!(
        (s.open_circuit_voltage_v()
            / (c.piezoelectric_stress_c_per_m2 / c.permittivity_f_per_m * expected)
            - 1.0)
            .abs()
            < 1e-10
    );
    assert!(s.electrode_charge_c().abs() < 1e-18);
    s.solve_transmit(10.0, 0.0).unwrap();
    let expected = -c.piezoelectric_stress_c_per_m2 / c.elastic_modulus_pa * 10.0;
    assert!((s.tip_displacement_m() / expected - 1.0).abs() < 1e-9);
}
#[test]
fn plasticity_has_elastic_unloading_and_neo_hookean_gradient() {
    let mut m = BilinearPlasticity1D::new(200.0, 2.0, 10.0).unwrap();
    let elastic = m.update(0.005).unwrap();
    assert_eq!(elastic.stress, 1.0);
    let plastic = m.update(0.03).unwrap();
    assert!(plastic.plastic_strain > 0.0);
    assert!((plastic.tangent - 200.0 * 10.0 / 210.0).abs() < 1e-12);
    let unload = m.update(0.025).unwrap();
    assert_eq!(unload.tangent, 200.0);
    assert_eq!(
        unload.equivalent_plastic_strain,
        plastic.equivalent_plastic_strain
    );
    let stretch = [1.1, 0.95, 1.03];
    let r = compressible_neo_hookean(stretch, 7.0, 11.0).unwrap();
    for i in 0..3 {
        let mut a = stretch;
        let mut b = stretch;
        a[i] += 1e-6;
        b[i] -= 1e-6;
        let derivative = (compressible_neo_hookean(a, 7.0, 11.0)
            .unwrap()
            .energy_density
            - compressible_neo_hookean(b, 7.0, 11.0)
                .unwrap()
                .energy_density)
            / 2e-6;
        assert!((derivative - r.nominal_stress[i]).abs() < 1e-8);
    }
}
