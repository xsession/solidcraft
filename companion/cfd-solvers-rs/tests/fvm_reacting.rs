//! Analytic checks for the species, flame and VOF units.
use cfd_solvers::fvm::*;

#[test]
fn fisher_kpp_front_moves_near_theoretical_speed() {
    let cfg = PremixedFlame1DConfig {
        cells: 400,
        length: 1.0,
        diffusivity: 1e-3,
        reaction_rate: 4.0,
        dt: 0.001,
        ..Default::default()
    };
    let mut flame = PremixedFlame1D::new(cfg).unwrap();
    flame.initialize_front(0.2, 0.02).unwrap();
    flame.step(1000);
    let x1 = flame.front_location();
    flame.step(1000);
    let speed = (flame.front_location() - x1) / 1.0;
    let theory = flame.theoretical_speed();
    assert!(
        (speed - theory).abs() < 0.1 * theory,
        "speed {speed} vs {theory}"
    );
}

#[test]
fn vof_stays_bounded_and_closed_box_conserves_volume() {
    let mesh = make_sheared_cartesian_hexa_mesh(8, 8, 1, 1.0, 1.0, 0.1, 0.0).unwrap();
    let flux: Vec<f64> = mesh
        .faces()
        .iter()
        .map(|f| {
            if f.neighbour.is_some() {
                f.area.x * 0.2
            } else {
                0.0
            }
        })
        .collect();
    let mut vof = VofTransport::new(
        mesh,
        VofConfig {
            dt: 0.05,
            compression: 0.5,
        },
    )
    .unwrap();
    vof.initialize_profile(|p| if p.x < 0.4 { 1.0 } else { 0.0 });
    vof.set_face_flux(flux).unwrap();
    let v0 = vof.volume();
    for _ in 0..20 {
        vof.step().unwrap();
    }
    assert!(vof
        .alpha()
        .iter()
        .all(|a| (-1e-12..=1.0 + 1e-12).contains(a)));
    assert!((vof.volume() - v0).abs() < 1e-12 * v0.max(1.0));
}

#[test]
fn invalid_inputs_are_errors_not_panics() {
    assert!(PremixedFlame1D::new(PremixedFlame1DConfig {
        dt: 10.0,
        ..Default::default()
    })
    .is_err());
    assert!(contact_angle_normal(Default::default(), Default::default(), 1.0).is_err());
    assert!(SootRadiationModel::new(SootRadiationConfig {
        soot_yield: -1.0,
        ..Default::default()
    })
    .is_err());
}
