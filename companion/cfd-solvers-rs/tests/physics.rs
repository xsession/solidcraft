use cfd_solvers::{fdtd::*, fvm::*, lbm::*};
#[test]
fn lbm_preserves_uniform_flow_and_mass() {
    let mut s = D2Q9Solver::new(D2Q9Config {
        nx: 12,
        ny: 10,
        ..Default::default()
    })
    .unwrap();
    s.initialize_uniform(1.2, 0.02, -0.01).unwrap();
    let mass = s.mass();
    s.step(30).unwrap();
    assert!((s.mass() - mass).abs() < 1e-3);
    assert!(s.velocity_x().iter().all(|x| (*x - 0.02).abs() < 1e-6));
    assert!(s.velocity_y().iter().all(|x| (*x + 0.01).abs() < 1e-6));
}
#[test]
fn taylor_green_decays_and_walls_drive_flow() {
    let mut s = D2Q9Solver::new(D2Q9Config {
        nx: 20,
        ny: 20,
        tau: 0.8,
        ..Default::default()
    })
    .unwrap();
    s.initialize_taylor_green(0.03).unwrap();
    let e = s.kinetic_energy();
    let m = s.mass();
    s.step(50).unwrap();
    assert!(s.kinetic_energy() < 0.5 * e);
    assert!((s.mass() - m).abs() < 0.01);
    let mut s = D2Q9Solver::new(D2Q9Config {
        nx: 12,
        ny: 12,
        ..Default::default()
    })
    .unwrap();
    for x in 0..12 {
        s.set_solid(x, 0, true).unwrap();
        s.set_wall_velocity(x, 11, 0.03, 0.0).unwrap();
    }
    s.initialize_uniform(1.0, 0.0, 0.0).unwrap();
    s.step(100).unwrap();
    assert!(s.average_ux(0, 12, 1, 11).unwrap() > 0.0);
}
#[test]
fn bouzidi_branches_match_formula() {
    assert_eq!(
        bouzidi_interpolated_bounce_back(0.25, 2.0, 4.0, 7.0).unwrap(),
        3.0
    );
    assert!(
        (bouzidi_interpolated_bounce_back(0.75, 2.0, 4.0, 7.0).unwrap()
            - (2.0 * 2.0 / 3.0 + 7.0 / 3.0))
            .abs()
            < 1e-6
    );
    assert!(bouzidi_interpolated_bounce_back(0.0, 1.0, 1.0, 1.0).is_err());
}
#[test]
fn absorbing_boundaries_remove_gaussian_energy() {
    for b in [Boundary1D::Mur1, Boundary1D::Pml, Boundary1D::Cpml] {
        let mut s = Maxwell1D::new(Maxwell1DConfig {
            cells: 128,
            boundary: b,
            pml_cells: 12,
            ..Default::default()
        })
        .unwrap();
        s.initialize_gaussian(0.4, 0.04).unwrap();
        let e = s.energy();
        s.step(400).unwrap();
        assert!(s.energy() < e * 0.01, "{b:?}: {}", s.energy() / e);
    }
}
#[test]
fn zero_dispersion_limits_agree_with_dielectric() {
    for model in 0..3 {
        let c = Maxwell1DConfig {
            cells: 80,
            courant: 0.5,
            ..Default::default()
        };
        let mut a = Maxwell1D::new(c).unwrap();
        let mut b = Maxwell1D::new(c).unwrap();
        a.set_material(25, 55, 2.0, 0.0).unwrap();
        match model {
            0 => b.set_debye_material(25..55, 2.0, 0.0, 1e-11, 0.0).unwrap(),
            1 => b.set_drude_material(25..55, 2.0, 0.0, 1e9, 0.0).unwrap(),
            _ => b
                .set_lorentz_material(25..55, 2.0, 0.0, 1e10, 1e9, 0.0)
                .unwrap(),
        };
        a.initialize_gaussian(0.2, 0.04).unwrap();
        b.initialize_gaussian(0.2, 0.04).unwrap();
        a.step(60).unwrap();
        b.step(60).unwrap();
        for (x, y) in a.electric().iter().zip(b.electric()) {
            assert!((x - y).abs() < 1e-12);
        }
    }
}
#[test]
fn lumped_and_dispersive_models_cannot_overlap() {
    let mut s = Maxwell1D::new(Maxwell1DConfig::default()).unwrap();
    s.set_debye_material(20..30, 2.0, 0.2, 1e-11, 0.0).unwrap();
    assert!(s
        .set_parallel_lumped_rlc(22, 1e-3, 1e-6, 50.0, f64::INFINITY, 0.0)
        .is_err());
    assert!(s.set_tfsf_source(22, |_| 1.0, 0.0).is_err());
}
#[test]
fn euler_uniform_state_and_periodic_conservation() {
    for reconstruction in [
        CompressibleReconstruction::MusclMinmod,
        CompressibleReconstruction::CharacteristicWeno5,
    ] {
        let mut s = CompressibleEuler1D::new(Compressible1DConfig {
            cells: 32,
            periodic: true,
            reconstruction,
            time_integrator: CompressibleTimeIntegrator::Ssprk3,
            ..Default::default()
        })
        .unwrap();
        s.initialize_uniform(EulerPrimitive1D {
            density: 1.2,
            velocity: 0.3,
            pressure: 2.0,
        })
        .unwrap();
        s.run(0.02).unwrap();
        assert!(s
            .conserved()
            .iter()
            .all(|u| (u.density - 1.2).abs() < 1e-13));
        let q: Vec<_> = (0..32)
            .map(|i| EulerPrimitive1D {
                density: 1.0 + 0.1 * (2.0 * std::f64::consts::PI * (i as f64 + 0.5) / 32.0).sin(),
                velocity: 0.5,
                pressure: 1.0,
            })
            .collect();
        s.initialize_profile(&q).unwrap();
        let (m, p, e) = (s.mass(), s.momentum(), s.total_energy());
        s.run(0.03).unwrap();
        assert!((s.mass() - m).abs() < 1e-12);
        assert!((s.momentum() - p).abs() < 1e-12);
        assert!((s.total_energy() - e).abs() < 1e-12);
    }
}
#[test]
fn shock_tube_is_physical_and_diffusion_converges() {
    let mut s = CompressibleEuler1D::new(Compressible1DConfig {
        cells: 100,
        ..Default::default()
    })
    .unwrap();
    s.initialize(
        EulerPrimitive1D {
            density: 1.0,
            velocity: 0.0,
            pressure: 1.0,
        },
        EulerPrimitive1D {
            density: 0.125,
            velocity: 0.0,
            pressure: 0.1,
        },
        0.5,
    )
    .unwrap();
    s.run(0.1).unwrap();
    assert!(s
        .pressure_jump_sensor()
        .unwrap()
        .iter()
        .all(|v| *v >= 0.0 && *v < 1.0));
    assert!((s.mass() - 0.5625).abs() < 1e-8);
    let mut d = Diffusion2D::new(Diffusion2DConfig {
        nx: 12,
        ny: 9,
        ..Default::default()
    })
    .unwrap();
    d.set_dirichlet(3.0, 3.0, 3.0, 3.0).unwrap();
    assert!(d.iterate(1500).unwrap() < 1e-10);
    assert!(d.field().iter().all(|u| (*u - 3.0).abs() < 1e-9));
}
