//! Analytic and conservation checks for the particle, spray, two-fluid, film, compressible VOF
//! and periodic 2-D flow units.
use cfd_solvers::{core::Vec3, fvm::*};

fn box_mesh() -> PolyMesh {
    make_sheared_cartesian_hexa_mesh(6, 5, 1, 0.6, 0.5, 0.1, 0.0).unwrap()
}

/// Discretely divergence-free recirculating face flux from the streamfunction
/// `psi = sin(pi x / 0.6) sin(pi y / 0.5)`, which vanishes on the box walls (0.1 m cells).
fn closed_flux(mesh: &PolyMesh) -> Vec<f64> {
    let psi = |x: f64, y: f64| {
        (std::f64::consts::PI * x / 0.6).sin() * (std::f64::consts::PI * y / 0.5).sin()
    };
    let (h, dz) = (0.05, 0.1);
    mesh.faces()
        .iter()
        .map(|f| {
            let c = f.center;
            if f.area.x.abs() > 1e-12 {
                f.area.x.signum() * (psi(c.x, c.y + h) - psi(c.x, c.y - h)) * dz
            } else if f.area.y.abs() > 1e-12 {
                -f.area.y.signum() * (psi(c.x + h, c.y) - psi(c.x - h, c.y)) * dz
            } else {
                0.0
            }
        })
        .collect()
}

#[test]
fn stokes_particle_relaxes_exponentially() {
    let mesh = box_mesh();
    let u = vec![Vec3::new(1.0, 0.0, 0.0); mesh.cell_count()];
    let mut cloud = ParticleCloud::new();
    let p = Particle {
        position: Vec3::new(0.3, 0.25, 0.05),
        ..Default::default()
    };
    cloud.add(p).unwrap();
    let mu = 1.8e-5;
    let tau = p.density * p.diameter * p.diameter / (18.0 * mu);
    cloud.advance(&mesh, &u, mu, 0.01, false).unwrap();
    let exact = 1.0 - (-0.01 / tau).exp();
    assert!((cloud.particles()[0].velocity.x - exact).abs() < 1e-14);
}

#[test]
fn interactions_conserve_momentum() {
    let mut ps = vec![
        Particle {
            position: Vec3::new(0.0, 0.0, 0.0),
            velocity: Vec3::new(3.0, 0.0, 0.0),
            diameter: 1e-3,
            density: 1000.0,
        },
        Particle {
            position: Vec3::new(8.5e-4, 1e-4, 0.0),
            velocity: Vec3::new(-2.0, 0.5, 0.0),
            diameter: 8e-4,
            density: 900.0,
        },
        Particle {
            position: Vec3::new(0.0, 0.01, 0.0),
            velocity: Vec3::new(20.0, 0.0, 0.0),
            diameter: 1e-3,
            density: 1000.0,
        },
        Particle {
            position: Vec3::new(9e-4, 0.01, 0.0),
            velocity: Vec3::new(-20.0, 0.0, 0.0),
            diameter: 1e-3,
            density: 1000.0,
        },
    ];
    let mut sys = ParticleInteractionSystem::new(ParticleInteractionConfig::default()).unwrap();
    sys.resolve(&mut ps).unwrap();
    let s = sys.last_stats();
    assert_eq!(s.collisions, 2);
    assert_eq!(s.breakups, 1);
    assert!((s.momentum_after - s.momentum_before).abs() < 1e-12 * s.momentum_before);
    assert!(s.kinetic_energy_after <= s.kinetic_energy_before);
}

#[test]
fn spray_mass_balance_closes() {
    let mesh = box_mesh();
    let u = vec![Vec3::new(0.5, 0.0, 0.0); mesh.cell_count()];
    let mut spray = SprayInjectionEvaporation::new(
        mesh,
        SprayConfig {
            dt: 1e-3,
            evaporation_constant: 5e-7,
            ..Default::default()
        },
    )
    .unwrap();
    spray
        .set_injection(SprayInjection {
            position: Vec3::new(0.05, 0.25, 0.05),
            diameter: 4e-5,
            duration: 5e-3,
            ..Default::default()
        })
        .unwrap();
    for _ in 0..20 {
        spray.step_uniform(&u, 400.0).unwrap();
    }
    let balance = spray.injected_mass() - spray.total_liquid_mass() - spray.evaporated_mass();
    assert!(balance.abs() < 1e-12 * spray.injected_mass());
    assert!(spray.evaporated_mass() > 0.0);
}

#[test]
fn euler_euler_conserves_phase_volume_and_momentum_without_body_forces() {
    let mut ee = EulerEulerTransport::new(
        box_mesh(),
        EulerEulerConfig {
            dt: 1e-3,
            drag_coefficient: 50.0,
            ..Default::default()
        },
    )
    .unwrap();
    ee.initialize_profile(
        |p| if p.x < 0.3 { 0.4 } else { 0.1 },
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.3, 0.1, 0.0),
    )
    .unwrap();
    let (v0, m0) = (ee.volume_fraction_integral(), ee.total_momentum());
    for _ in 0..20 {
        ee.step();
    }
    assert!((ee.volume_fraction_integral() - v0).abs() < 1e-13);
    let dm = ee.total_momentum() - m0;
    // Boundary faces are impermeable but carry no pressure, so momentum is conserved only
    // up to the momentum flux that would cross them; with zero boundary flux it is exact.
    assert!(dm.norm() < 1e-9, "momentum drift {dm:?}");
    assert!(ee.minimum_fraction() >= 0.0 && ee.maximum_fraction() <= 1.0);
}

#[test]
fn closed_thin_film_and_compressible_vof_conserve_inventory() {
    let mesh = box_mesh();
    let mut film = ThinFilmTransport::new(
        mesh.clone(),
        ThinFilmConfig {
            dt: 1e-3,
            surface_tension: 0.07,
            viscosity: 0.1,
            ..Default::default()
        },
    )
    .unwrap();
    film.initialize_profile(|p| 1e-3 * (1.0 + 0.2 * (10.0 * p.x).sin()));
    let h0 = film.inventory();
    for _ in 0..10 {
        film.step().unwrap();
    }
    assert!((film.inventory() - h0).abs() < 1e-12 * h0);

    let mut cv = CompressibleVofTransport::new(
        mesh.clone(),
        CompressibleVofConfig {
            dt: 1e-3,
            ..Default::default()
        },
    )
    .unwrap();
    cv.initialize_profile(|p| if p.x < 0.2 { 1.0 } else { 0.0 }, 0.0)
        .unwrap();
    cv.set_face_flux(closed_flux(&mesh)).unwrap();
    let (m0, v0) = (cv.total_mass(), cv.liquid_volume());
    for _ in 0..10 {
        cv.step().unwrap();
    }
    assert!((cv.total_mass() - m0).abs() < 1e-12 * m0);
    assert!((cv.liquid_volume() - v0).abs() < 1e-12 * v0.max(1.0));
}

#[test]
fn projection_removes_divergence_and_taylor_green_decays_at_viscous_rate() {
    let mut proj = Projection2D::new(Projection2DConfig {
        nx: 16,
        ny: 16,
        ..Default::default()
    })
    .unwrap();
    proj.initialize_divergent(0.05);
    assert!(proj.divergence_l2() > 0.1);
    proj.project().unwrap();
    assert!(proj.divergence_l2() < 1e-10);
    assert!(proj.pressure_result().converged);

    let cfg = Incompressible2DConfig {
        nx: 32,
        ny: 32,
        kinematic_viscosity: 0.01,
        dt: 1e-3,
        ..Default::default()
    };
    let mut ns = Incompressible2D::new(cfg).unwrap();
    ns.initialize_taylor_green(0.1);
    let e0 = ns.kinetic_energy();
    ns.run(100).unwrap();
    let k2 = 2.0 * (2.0 * std::f64::consts::PI).powi(2);
    let exact = e0 * (-2.0 * cfg.kinematic_viscosity * k2 * ns.time()).exp();
    assert!(
        (ns.kinetic_energy() - exact).abs() < 0.02 * exact,
        "{} vs {exact}",
        ns.kinetic_energy()
    );
    assert!(ns.taylor_green_velocity_error(0.1) < 2e-3);
    assert!(ns.divergence_l2() < 1e-10);
}

#[test]
fn invalid_controls_are_errors() {
    assert!(Projection2D::new(Projection2DConfig {
        nx: 2,
        ..Default::default()
    })
    .is_err());
    assert!(Incompressible2D::new(Incompressible2DConfig {
        dt: 0.0,
        ..Default::default()
    })
    .is_err());
    assert!(CompressibleVofTransport::new(
        box_mesh(),
        CompressibleVofConfig {
            gas_gamma: 1.0,
            ..Default::default()
        }
    )
    .is_err());
    assert!(EulerEulerTransport::new(
        box_mesh(),
        EulerEulerConfig {
            minimum_phase_fraction: 0.6,
            ..Default::default()
        }
    )
    .is_err());
    assert!(ParticleInteractionSystem::new(ParticleInteractionConfig {
        breakup_weber: 1.0,
        ..Default::default()
    })
    .is_err());
    assert!(DenseParticleRheology::new(
        box_mesh(),
        DenseRheologyConfig {
            maximum_solid_fraction: 1.0,
            ..Default::default()
        }
    )
    .is_err());
    let mut cloud = ParticleCloud::new();
    assert!(cloud
        .add(Particle {
            diameter: 0.0,
            ..Default::default()
        })
        .is_err());
}
