use cfd_solvers::{acoustics::*, chemistry::*, core::Vec3, optics::*, rf::*, Complex64 as C};
#[test]
fn matched_quarter_wave_line_and_network_conversions() {
    let line = TransmissionLine {
        characteristic_impedance: 50.0,
        propagation_constant: C::new(0.0, std::f64::consts::FRAC_PI_2),
        length: 1.0,
    };
    let s = line.s_parameters(50.0).unwrap();
    assert!(s.a11.norm() < 1e-14 && s.a22.norm() < 1e-14);
    assert!((s.a21 - C::new(0.0, -1.0)).norm() < 1e-14);
    let a = s_to_abcd(s, 50.0).unwrap();
    let round = abcd_to_s(a, 50.0).unwrap();
    for (x, y) in s.values().iter().zip(round.values()) {
        assert!((*x - y).norm() < 1e-14);
    }
    let z = Matrix2C {
        a11: C::new(100.0, 5.0),
        a12: C::new(7.0, 1.0),
        a21: C::new(7.0, 1.0),
        a22: C::new(75.0, -3.0),
    };
    let round = s_to_z(z_to_s(z, 50.0).unwrap(), 50.0).unwrap();
    for (x, y) in z.values().iter().zip(round.values()) {
        assert!((*x - y).norm() < 1e-12);
    }
}
#[test]
fn touchstone_round_trip_preserves_order_and_phase() {
    let input = "! wrapped record\n# MHz S MA R 75\n 100 0.1 180 0.8 -90\n 0.7 45 0.2 0 ! end\n";
    let (p, z) = read_touchstone_s2p(input).unwrap();
    assert_eq!(z, 75.0);
    assert_eq!(p[0].frequency_hz, 1e8);
    assert!((p[0].s.a21 - C::new(0.0, -0.8)).norm() < 1e-14);
    let text = write_touchstone_s2p(&p, z).unwrap();
    let (round, _) = read_touchstone_s2p(&text).unwrap();
    for (x, y) in p[0].s.values().iter().zip(round[0].s.values()) {
        assert!((*x - y).norm() < 1e-15);
    }
    assert!(read_touchstone_s2p("# HZ S RI R 50\n1 0 0\n").is_err());
    assert!(read_touchstone_s2p("# HZ Y RI R 50\n").is_err());
}
#[test]
fn snell_total_internal_reflection_and_coating() {
    let transmitted = refract(
        Vec3::new(0.5, 0.0, 0.75_f64.sqrt()),
        Vec3::new(0.0, 0.0, -1.0),
        1.0,
        1.5,
    )
    .unwrap()
    .unwrap();
    assert!((transmitted.x - 1.0 / 3.0).abs() < 1e-14);
    assert!(refract(
        Vec3::new(0.9, 0.0, 0.19_f64.sqrt()),
        Vec3::new(0.0, 0.0, -1.0),
        1.5,
        1.0
    )
    .unwrap()
    .is_none());
    let f = fresnel_dielectric(1.0, 1.5, 1.0).unwrap();
    assert!((f.reflectance_s - 0.04).abs() < 1e-14);
    let n = 1.5_f64.sqrt();
    let r = thin_film_reflectance_normal(
        1.0,
        1.5,
        550.0,
        &[ThinFilmLayer {
            refractive_index: n,
            thickness_nm: 550.0 / (4.0 * n),
        }],
    )
    .unwrap();
    assert!(r < 1e-28);
}
#[test]
fn beam_and_paraxial_focus() {
    let mut b = GaussianBeam::new(500e-9, 1e-3, 0.0, 1.0).unwrap();
    let z = b.q().im;
    b.propagate(z).unwrap();
    assert!((b.radius() - 2.0_f64.sqrt() * 1e-3).abs() < 1e-15);
    assert!(b.transform(1.0, 0.0, 0.0, 2.0, 1.0).is_err());
    let mut sys = SequentialOpticalSystem::new(1.0).unwrap();
    sys.add_surface(SequentialSurface {
        surface_type: SurfaceType::Sphere,
        vertex_z: 0.01,
        radius: 0.05,
        aperture_radius: 0.005,
        refractive_index_after: 1.5,
        ..Default::default()
    })
    .unwrap();
    let focal = paraxial_back_focal_distance(&sys, 1.0).unwrap();
    assert!((focal - 0.15).abs() < 1e-14);
    let result = sys
        .trace(Ray {
            origin: Vec3::new(1e-5, 0.0, 0.0),
            ..Default::default()
        })
        .unwrap();
    assert!(result.valid);
    let t = -result.ray.origin.x / result.ray.direction.x;
    let focus = result.ray.origin.z + t * result.ray.direction.z;
    assert!((focus - 0.16).abs() < 1e-7);
}
fn network() -> ReactionNetwork {
    let mut n = ReactionNetwork::new(vec![
        Species {
            name: "A".into(),
            molar_mass: 0.01,
            charge: 0,
        },
        Species {
            name: "B".into(),
            molar_mass: 0.01,
            charge: 0,
        },
    ]);
    n.add_reaction(ElementaryReaction {
        reactants: vec![StoichiometricTerm {
            species: 0,
            coefficient: 1.0,
        }],
        products: vec![StoichiometricTerm {
            species: 1,
            coefficient: 1.0,
        }],
        forward: ArrheniusRate {
            pre_exponential: 12.0,
            ..Default::default()
        },
        ..Default::default()
    })
    .unwrap();
    n
}
#[test]
fn reaction_rates_third_body_and_falloff() {
    let mut n = network();
    assert_eq!(
        n.source_terms(&[2.0, 0.0], 300.0).unwrap(),
        vec![-24.0, 24.0]
    );
    let mut r = n.reactions()[0].clone();
    r.third_body = true;
    r.third_body_efficiencies = vec![ThirdBodyEfficiency {
        species: 1,
        efficiency: 2.0,
    }];
    r.low_pressure_limit = Some(ArrheniusRate {
        pre_exponential: 3.0,
        ..Default::default()
    });
    n.add_reaction(r).unwrap();
    let rates = n.reaction_rates(&[2.0, 1.0], 300.0).unwrap();
    assert!((rates[1] - 12.0).abs() < 1e-14);
}
#[test]
fn implicit_reactor_matches_decay_and_rolls_back_failure() {
    let n = network();
    let mut c = [1.0, 0.0];
    let result = integrate_isothermal(
        &n,
        &mut c,
        300.0,
        0.2,
        ReactorConfig {
            relative_tolerance: 1e-5,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(result.accepted_steps > 0);
    assert!((c[0] - (-2.4_f64).exp()).abs() < 4e-4);
    assert!((c.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    let old = c;
    assert!(integrate_isothermal(
        &n,
        &mut c,
        300.0,
        1.0,
        ReactorConfig {
            max_steps: 1,
            ..Default::default()
        }
    )
    .is_err());
    assert_eq!(c, old);
}
#[test]
fn phase_potentials_and_precipitation() {
    let species = vec![
        ThermoSpecies {
            name: "A".into(),
            standard_chemical_potential: 10.0,
        },
        ThermoSpecies {
            name: "B".into(),
            standard_chemical_potential: 20.0,
        },
    ];
    let phase = IdealGasPhase::new(species, 1e5).unwrap();
    let a = phase.activities(&[2.0, 2.0], 300.0, 2e5).unwrap();
    assert_eq!(a, vec![1.0, 1.0]);
    assert_eq!(
        phase.chemical_potentials(&[2.0, 2.0], 300.0, 2e5).unwrap(),
        vec![10.0, 20.0]
    );
    assert!(phase.chemical_potentials(&[1.0, 0.0], 300.0, 1e5).is_err());
    let eq = equilibrate_binary_salt(3.0, 2.0, 0.0, 0.25).unwrap();
    assert!((eq.saturation_ratio - 1.0).abs() < 1e-13);
    assert!((eq.dissolved_a + eq.solid_amount - 3.0).abs() < 1e-14);
}
#[test]
fn acoustic_fractional_derivative_and_absorption_match_mode() {
    let (nx, ny) = (16, 16);
    let h = 2.0 * std::f64::consts::PI / 16.0;
    let p: Vec<_> = (0..ny)
        .flat_map(|_| (0..nx).map(|i| (2.0 * i as f64 * h).sin()))
        .collect();
    let f = fractional_laplacian_periodic_2d(&p, nx, ny, h, h, 0.5).unwrap();
    for (x, y) in f.iter().zip(&p) {
        assert!((x - 2.0 * y).abs() < 1e-13);
    }
    let a = PowerLawAbsorption2D {
        alpha0_np_per_m: 0.1,
        exponent: 1.5,
        reference_frequency_hz: 1.0,
    };
    let c = 2.0;
    let dt = 0.3;
    let atten = power_law_attenuate_periodic_2d(&p, nx, ny, h, h, c, dt, a).unwrap();
    let factor = (-0.1 * (c * 2.0 / (2.0 * std::f64::consts::PI)).powf(1.5) * c * dt).exp();
    for (x, y) in atten.iter().zip(&p) {
        assert!((x - factor * y).abs() < 1e-13);
    }
}
#[test]
fn acoustic_uniform_field_and_sensor_timing() {
    let mut s = KSpaceAcoustic2D::new(KSpace2DConfig {
        nx: 16,
        ny: 16,
        pml_cells: 0,
        ..Default::default()
    })
    .unwrap();
    s.add_sensor(77).unwrap();
    s.initialize_pressure(&vec![3.0; 256]).unwrap();
    s.step(8).unwrap();
    assert!(s.pressure().iter().all(|p| (*p - 3.0).abs() < 1e-13));
    assert_eq!(s.sensors()[0].pressure.len(), 9);
    assert!((s.time() - 8.0 * s.dt()).abs() < 1e-20);
    let reconstructed = s
        .time_reversal_reconstruct(&[77], &[vec![1.0, 2.0, 3.0]], true)
        .unwrap();
    assert_eq!(reconstructed[77], 1.0);
    assert!(reconstructed.iter().all(|p| *p >= 0.0));
}
