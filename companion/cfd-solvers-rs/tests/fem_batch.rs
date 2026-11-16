//! Analytic checks for the FEM flow, dynamics, magnetics and kinematics units.
use cfd_solvers::{core::LinearConfig, fem::*};

fn lin() -> LinearConfig {
    LinearConfig {
        max_iterations: 4000,
        relative_tolerance: 1e-12,
        ..Default::default()
    }
}

#[test]
fn newmark_average_acceleration_conserves_energy_without_load() {
    let mut bar = DynamicBar1D::new(DynamicBar1DConfig {
        elements: 20,
        dt: 1e-6,
        ..Default::default()
    })
    .unwrap();
    bar.initialize(|x| 1e-4 * (std::f64::consts::FRAC_PI_2 * x).sin(), None)
        .unwrap();
    let e0 = bar.kinetic_energy().unwrap() + bar.strain_energy().unwrap();
    bar.run(400, 0.0).unwrap();
    let e1 = bar.kinetic_energy().unwrap() + bar.strain_energy().unwrap();
    assert!(
        (e1 - e0).abs() < 1e-8 * e0,
        "energy drift {}",
        (e1 - e0) / e0
    );
}

#[test]
fn darcy_recovers_one_dimensional_linear_pressure_and_velocity() {
    let mesh = make_rectangle_tri_mesh(8, 3, 2.0, 1.0).unwrap();
    let mut d = Darcy2D::new(
        mesh,
        Darcy2DConfig {
            viscosity: 0.5,
            density: 0.0,
            linear: lin(),
            ..Default::default()
        },
    )
    .unwrap();
    d.set_boundary(0, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 3.0)));
    d.set_boundary(1, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 1.0)));
    d.solve(|_| 2.0, None).unwrap();
    for (p, node) in d.pressure().iter().zip(&d.mesh().nodes) {
        assert!((p - (3.0 - node.x)).abs() < 1e-9);
    }
    for v in d.element_velocity(|_| 2.0).unwrap() {
        assert!((v.x - 4.0).abs() < 1e-8 && v.y.abs() < 1e-8);
    }
}

#[test]
fn axisymmetric_rigid_axial_shift_is_stress_free() {
    let mut mesh = make_rectangle_tri_mesh(4, 4, 1.0, 1.0).unwrap();
    for p in &mut mesh.nodes {
        p.x += 1.0;
    }
    let mut ax = AxisymmetricElasticity::new(
        mesh,
        AxisymmetricElasticityConfig {
            relative_tolerance: 1e-12,
            ..Default::default()
        },
    )
    .unwrap();
    ax.set_dirichlet(
        |p| p.y <= 1e-12,
        |_| Displacement2 { x: 0.0, y: 0.25 },
        true,
        true,
    );
    ax.solve(None).unwrap();
    assert!(ax
        .displacement()
        .iter()
        .all(|d| d.x.abs() < 1e-9 && (d.y - 0.25).abs() < 1e-9));
    assert!(AxisymmetricElasticity::new(
        make_rectangle_tri_mesh(2, 2, 1.0, 1.0)
            .map(|mut m| {
                m.nodes[0].x = -1.0;
                m
            })
            .unwrap(),
        AxisymmetricElasticityConfig::default()
    )
    .is_err());
}

#[test]
fn penalty_stokes_drives_divergence_down_with_penalty() {
    let mesh = make_rectangle_tri_mesh(8, 8, 1.0, 1.0).unwrap();
    let run = |penalty: f64| {
        let mut st = PenaltyStokes2D::new(
            mesh.clone(),
            PenaltyStokes2DConfig {
                divergence_penalty: penalty,
                relative_tolerance: 1e-12,
                max_iterations: 20000,
                ..Default::default()
            },
        )
        .unwrap();
        st.set_dirichlet(|p| Displacement2 {
            x: if p.y > 1.0 - 1e-12 { 1.0 } else { 0.0 },
            y: 0.0,
        });
        st.solve(None).unwrap();
        st.divergence_l2().unwrap()
    };
    assert!(run(1e3) < 0.2 * run(10.0));
}

#[test]
fn magnetostatics_linear_material_converges_in_one_picard_step_and_coil_is_consistent() {
    let mesh = make_rectangle_tri_mesh(6, 6, 1.0, 1.0).unwrap();
    let mut m = Magnetostatics2D::new(mesh.clone(), lin()).unwrap();
    for patch in 0..4 {
        m.set_boundary(
            patch,
            ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 0.0)),
        );
    }
    let coil = RectangularStrandedCoil2D {
        xmin: 0.3,
        xmax: 0.7,
        ymin: 0.3,
        ymax: 0.7,
        turns: 10.0,
        current_a: 1.0,
    };
    let r = m
        .solve_nonlinear(
            |_| 2.0,
            |p| coil.current_density_a_per_m2(p).unwrap(),
            10,
            1e-12,
            1.0,
        )
        .unwrap();
    assert!(r.converged && r.iterations == 1);
    let l = coil_inductance_h(&mesh, m.vector_potential(), &coil).unwrap();
    // Energy identity: L I^2 = integral A J, positive for a positive coil current.
    assert!(l > 0.0);
    let curve = PiecewiseLinearBHCurve::new(vec![0.0, 1.0], vec![0.0, 2.0]).unwrap();
    assert_eq!(curve.reluctivity(0.5).unwrap(), 2.0);
    assert!(PiecewiseLinearBHCurve::new(vec![0.0, 1.0, 0.5], vec![0.0, 1.0, 2.0]).is_err());
}

#[test]
fn nonlinear_poisson_reproduces_a_manufactured_linear_solution() {
    let mesh = make_rectangle_tri_mesh(6, 6, 1.0, 1.0).unwrap();
    let mut np = NonlinearPoisson2D::new(
        mesh,
        NonlinearPoisson2DConfig {
            cubic_coefficient: 0.0,
            ..Default::default()
        },
    )
    .unwrap();
    np.solve(|_| 0.0, |p| 1.0 + 2.0 * p.x - p.y).unwrap();
    for (u, p) in np.solution().iter().zip(&np.mesh().nodes) {
        assert!((u - (1.0 + 2.0 * p.x - p.y)).abs() < 1e-8);
    }
}

#[test]
fn kinematics_reduce_to_identity_for_rigid_motion() {
    let mesh = make_rectangle_tri_mesh(1, 1, 1.0, 1.0).unwrap();
    let t = mesh.triangles[0];
    let (c, s) = (0.6f64.cos(), 0.6f64.sin());
    let cur = t.node.map(|i| {
        let p = mesh.nodes[i];
        Node2 {
            x: c * p.x - s * p.y + 3.0,
            y: s * p.x + c * p.y - 1.0,
        }
    });
    let st = total_lagrangian_triangle(&mesh, &t, &cur, 1e6, 0.3).unwrap();
    assert!(st
        .green_lagrange_strain
        .iter()
        .chain(&st.internal_force)
        .all(|v| v.abs() < 1e-9));
    assert!((st.determinant - 1.0).abs() < 1e-12);
    let truss = green_lagrange_truss(1.0, 1.0, 1.0, 1.0).unwrap();
    assert_eq!(truss.axial_force, 0.0);
    let pcm = PhaseChangeMaterial {
        solidus: 0.0,
        liquidus: 1.0,
        heat_capacity_solid: 1.0,
        heat_capacity_liquid: 1.0,
        latent_heat: 10.0,
    };
    for t in [-3.0, 0.4, 7.0] {
        let h = pcm.specific_enthalpy(t).unwrap();
        assert!((pcm.temperature_from_enthalpy(h).unwrap() - t).abs() < 1e-9);
    }
}
