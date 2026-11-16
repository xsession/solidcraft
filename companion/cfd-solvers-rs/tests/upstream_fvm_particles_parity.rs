//! Value-by-value comparison with the C++ upstream FVM particle, interaction, dense-rheology,
//! spray and Euler-Euler units. Goldens come from `parity/oracle_fvm_particles.cpp`.
use cfd_solvers::{core::Vec3, fvm::*};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-fvm-particles.json")).unwrap()
}

#[track_caller]
fn same(name: &str, rust: &[f64], tol: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    // Relative per entry, with a floor at a millionth of the array's largest magnitude so that
    // round-off around exact zeros does not count as a mismatch.
    let floor = 1e-6 * cpp.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= tol * c.abs().max(floor),
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn flat3(v: &[Vec3]) -> Vec<f64> {
    v.iter().flat_map(|a| [a.x, a.y, a.z]).collect()
}
fn flatp(v: &[Particle]) -> Vec<f64> {
    v.iter()
        .flat_map(|p| {
            [
                p.position.x,
                p.position.y,
                p.position.z,
                p.velocity.x,
                p.velocity.y,
                p.velocity.z,
                p.diameter,
                p.density,
            ]
        })
        .collect()
}
fn part(position: Vec3, velocity: Vec3, diameter: f64, density: f64) -> Particle {
    Particle {
        position,
        velocity,
        diameter,
        density,
    }
}

fn setup() -> (PolyMesh, Vec<Vec3>) {
    let mesh = make_sheared_cartesian_hexa_mesh(5, 4, 1, 0.5, 0.4, 0.1, 0.0).unwrap();
    let u = mesh
        .cells()
        .iter()
        .map(|c| v(0.2 + c.center.y, -0.1 * c.center.x, 0.0))
        .collect();
    (mesh, u)
}

#[test]
fn lagrangian_particles_match_upstream() {
    let (mesh, u) = setup();
    let mut cloud = ParticleCloud::new();
    cloud
        .add(part(v(0.05, 0.05, 0.05), v(0.0, 0.0, 0.0), 2e-4, 1500.0))
        .unwrap();
    cloud
        .add(part(v(0.31, 0.22, 0.05), v(0.5, 0.1, 0.0), 1e-4, 1000.0))
        .unwrap();
    cloud
        .add(part(v(0.44, 0.37, 0.05), v(-0.2, 0.0, 0.1), 5e-5, 800.0))
        .unwrap();
    for _ in 0..4 {
        cloud.advance(&mesh, &u, 1.8e-5, 0.01, true).unwrap();
    }
    same("cloud", &flatp(cloud.particles()), 1e-12);
    same("cloud_source", &flat3(cloud.fluid_momentum_source()), 1e-10);
    let mut forces = flat3(&[saffman_lift_force(
        v(0.3, -0.1, 0.05),
        v(0.0, 0.0, 40.0),
        1.2,
        1.5e-5,
        1e-4,
        1.615,
    )
    .unwrap()]);
    forces.extend(flat3(&[virtual_mass_force(
        v(1.0, 2.0, 0.0),
        v(0.5, -1.0, 0.3),
        1000.0,
        2e-3,
        0.5,
    )
    .unwrap()]));
    same("forces", &forces, 1e-12);
    let mut a = part(v(0.0, 0.0, 0.0), v(1.0, 0.2, 0.0), 1e-3, 1000.0);
    let mut b = part(v(8e-4, 3e-4, 0.0), v(-0.5, 0.0, 0.1), 7e-4, 900.0);
    collide_particles_elastic(&mut a, &mut b, 0.8).unwrap();
    let merged = coalesce_particles(&a, &b);
    same("collide_coalesce", &flatp(&[a, b, merged]), 1e-12);
    let flag = |x: bool| if x { 1.0 } else { 0.0 };
    same(
        "breakup_evap",
        &[
            flag(breakup_by_weber(1.2, 30.0, 1e-3, 0.072, 12.0).unwrap()),
            flag(breakup_by_weber(1.2, 5.0, 1e-3, 0.072, 12.0).unwrap()),
            evaporate_d2_law(1e-4, 1e-9, 2.0).unwrap(),
            evaporate_d2_law(1e-4, 1e-8, 2.0).unwrap(),
        ],
        1e-14,
    );
    same(
        "volume_fraction",
        &particle_volume_fraction(&mesh, cloud.particles()),
        1e-12,
    );
}

#[test]
fn interactions_match_upstream() {
    let mut ps = vec![
        part(v(0.0, 0.0, 0.0), v(0.05, 0.0, 0.0), 1e-3, 1000.0),
        part(v(9e-4, 0.0, 0.0), v(-0.05, 0.0, 0.0), 1e-3, 1000.0),
        part(v(0.01, 0.0, 0.0), v(15.0, 1.0, 0.0), 1.2e-3, 1000.0),
        part(v(0.0109, 2e-4, 0.0), v(-15.0, 0.0, 0.0), 1e-3, 1000.0),
        part(v(0.05, 0.05, 0.05), v(0.0, 0.0, 0.0), 1e-3, 1000.0),
    ];
    let mut sys = ParticleInteractionSystem::new(ParticleInteractionConfig::default()).unwrap();
    sys.resolve(&mut ps).unwrap();
    let s = *sys.last_stats();
    same("interactions", &flatp(&ps), 1e-12);
    same(
        "interaction_stats",
        &[
            s.candidate_pairs as f64,
            s.collisions as f64,
            s.coalescences as f64,
            s.breakups as f64,
            s.removed_particles as f64,
            s.momentum_before,
            s.momentum_after,
            s.kinetic_energy_before,
            s.kinetic_energy_after,
        ],
        1e-12,
    );
}

#[test]
fn dense_rheology_matches_upstream() {
    let (mesh, u) = setup();
    let dense = DenseParticleRheology::new(mesh, DenseRheologyConfig::default()).unwrap();
    let grains: Vec<Particle> = (0..40)
        .map(|i| {
            part(
                v(0.02 + 0.011 * i as f64, 0.03 + 0.008 * (i % 7) as f64, 0.05),
                v(0.1 * (i as f64).sin(), 0.05, 0.0),
                4e-3,
                2500.0,
            )
        })
        .collect();
    let r = dense.evaluate(&grains, &u).unwrap();
    same("dense_phi", &r.solid_fraction, 1e-12);
    same("dense_pressure", &r.granular_pressure, 1e-10);
    same("dense_viscosity", &r.effective_viscosity, 1e-10);
    same(
        "dense_particle_source",
        &flat3(&r.particle_momentum_source),
        1e-10,
    );
    assert_eq!(
        flat3(&r.carrier_momentum_source),
        flat3(&r.particle_momentum_source)
            .iter()
            .map(|x| -x)
            .collect::<Vec<_>>()
    );
}

#[test]
fn spray_matches_upstream() {
    let (mesh, u) = setup();
    let temp: Vec<f64> = mesh
        .cells()
        .iter()
        .map(|c| 350.0 + 100.0 * c.center.x)
        .collect();
    let mut spray = SprayInjectionEvaporation::new(
        mesh,
        SprayConfig {
            dt: 1e-3,
            evaporation_constant: 2e-7,
            ..Default::default()
        },
    )
    .unwrap();
    spray
        .set_injection(SprayInjection {
            position: v(0.02, 0.2, 0.05),
            velocity: v(5.0, 0.0, 0.0),
            diameter: 5e-5,
            temperature: 290.0,
            mass_flow_rate: 1e-6,
            duration: 4.5e-3,
            ..Default::default()
        })
        .unwrap();
    spray
        .add_parcel(SprayParcel {
            position: v(0.3, 0.1, 0.05),
            velocity: v(0.0, 1.0, 0.0),
            diameter: 2e-5,
            mass: 4e-12,
            temperature: 310.0,
        })
        .unwrap();
    for _ in 0..8 {
        spray.step(&u, &temp).unwrap();
    }
    let parcels: Vec<f64> = spray
        .parcels()
        .iter()
        .flat_map(|p| {
            [
                p.position.x,
                p.position.y,
                p.position.z,
                p.velocity.x,
                p.velocity.y,
                p.velocity.z,
                p.diameter,
                p.mass,
                p.temperature,
            ]
        })
        .collect();
    same("spray_parcels", &parcels, 1e-11);
    same("spray_vapor", spray.vapor_mass_source(), 1e-10);
    same("spray_energy", spray.carrier_energy_source(), 1e-10);
    same(
        "spray_momentum",
        &flat3(spray.carrier_momentum_source()),
        1e-10,
    );
    same(
        "spray_totals",
        &[
            spray.total_liquid_mass(),
            spray.injected_mass(),
            spray.evaporated_mass(),
            spray.time(),
            spray.steps() as f64,
        ],
        1e-12,
    );
}

#[test]
fn euler_euler_matches_upstream() {
    let (mesh, _) = setup();
    let cfg = EulerEulerConfig {
        dt: 2e-3,
        dispersed_density: 50.0,
        drag_coefficient: 20.0,
        quadratic_drag: true,
        dispersed_body_acceleration: v(0.0, 9.81, 0.0),
        primary_body_acceleration: v(0.0, -0.5, 0.0),
        ..Default::default()
    };
    let mut ee = EulerEulerTransport::new(mesh, cfg).unwrap();
    ee.initialize_profile(
        |p| if p.y < 0.2 { 0.3 } else { 0.05 },
        v(0.1, 0.0, 0.0),
        v(0.0, 0.2, 0.0),
    )
    .unwrap();
    for _ in 0..10 {
        ee.step();
    }
    same("ee_alpha", ee.dispersed_fraction(), 1e-11);
    same("ee_primary", &flat3(ee.primary_velocity()), 1e-10);
    same("ee_dispersed", &flat3(ee.dispersed_velocity()), 1e-10);
    let mut tot = vec![
        ee.volume_fraction_integral(),
        ee.minimum_fraction(),
        ee.maximum_fraction(),
        ee.time(),
    ];
    tot.extend(flat3(&[ee.total_momentum(), ee.integrated_drag_force()]));
    same("ee_totals", &tot, 1e-10);
    same("ee_face_flux", ee.dispersed_face_flux(), 1e-10);
}
