//! Comparisons with separately compiled, source-pinned C++ numerical output.
use cfd_solvers::{core::LinearConfig, fdtd::*, fem::*, fvm::*, lbm::*};
use serde_json::Value;
fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens.json")).unwrap()
}

#[test]
fn structural_elastic_modal_piezo_and_material_results_match_cpp() {
    let g = goldens();
    for mode in 0..2 {
        let mesh = make_rectangle_tri_mesh(5, 4, 1.0, 1.0).unwrap();
        let temperature: Vec<_> = mesh.nodes.iter().map(|p| 310.0 + 5.0 * p.x).collect();
        let mut s = Elasticity2D::new(
            mesh,
            Elasticity2DConfig {
                youngs_modulus: 3.0,
                mode: if mode == 0 {
                    ElasticityMode2D::PlaneStress
                } else {
                    ElasticityMode2D::PlaneStrain
                },
                ..Default::default()
            },
        )
        .unwrap();
        s.set_dirichlet(|p| p.x == 0.0, |_| Displacement2::default(), true, true)
            .unwrap();
        s.solve_thermal(&temperature, 0.002, 300.0, |_| Displacement2 {
            x: 0.1,
            y: -0.2,
        })
        .unwrap();
        let out: Vec<_> = s.displacement().iter().flat_map(|p| [p.x, p.y]).collect();
        compare(&g, &format!("elasticity_{mode}"), &out, 1e-9, 1e-9);
    }
    let modes = bar_modal1d(
        BarModal1DConfig {
            elements: 16,
            young_modulus: 100.0,
            density: 4.0,
            ..Default::default()
        },
        2,
    )
    .unwrap();
    let mut out = vec![];
    for mode in modes {
        out.push(mode.angular_frequency);
        out.extend(mode.displacement);
    }
    compare(&g, "modal_bar", &out, 1e-8, 1e-8);
    let mut p = PiezoelectricBar1D::new(PiezoelectricBar1DConfig::default()).unwrap();
    p.solve_receive(1.2).unwrap();
    let mut out = p.displacement_m().to_vec();
    out.extend_from_slice(p.potential_v());
    out.push(p.electrode_charge_c());
    p.solve_transmit(10.0, 0.0).unwrap();
    out.extend_from_slice(p.displacement_m());
    out.extend_from_slice(p.potential_v());
    out.push(p.electrode_charge_c());
    compare(&g, "piezo_bar", &out, 1e-15, 1e-11);
    let mut law = BilinearPlasticity1D::new(200.0, 2.0, 10.0).unwrap();
    let mut out = vec![];
    for strain in [0.005, 0.03, 0.025, -0.03, 0.0] {
        let r = law.update(strain).unwrap();
        out.extend([
            r.stress,
            r.tangent,
            r.plastic_strain,
            r.equivalent_plastic_strain,
        ]);
    }
    let r = compressible_neo_hookean([1.1, 0.95, 1.03], 7.0, 11.0).unwrap();
    out.push(r.energy_density);
    out.extend(r.nominal_stress);
    compare(&g, "material_laws", &out, 1e-13, 1e-12);
}
fn compare(g: &Value, key: &str, actual: &[f64], absolute: f64, relative: f64) {
    let expected = g[key].as_array().unwrap();
    assert_eq!(actual.len(), expected.len(), "{key} output size");
    for (i, (a, b)) in actual.iter().zip(expected).enumerate() {
        let b = b.as_f64().unwrap();
        let bound = absolute + relative * b.abs();
        assert!(
            a.is_finite() && (a - b).abs() <= bound,
            "{key}[{i}]: Rust={a:.17e}, C++={b:.17e}, error={:.3e}, bound={bound:.3e}",
            (a - b).abs()
        );
    }
}
#[test]
fn poisson_and_diffusion_match_cpp() {
    let g = goldens();
    let mut p = Poisson1D::new(19, 1.7).unwrap();
    p.solve(|x| 2.0 + x * x, 0.2, 0.7).unwrap();
    compare(&g, "poisson1d", &p.u, 1e-13, 1e-12);
    let mut s = Diffusion2D::new(Diffusion2DConfig {
        nx: 9,
        ny: 8,
        lx: 1.5,
        ly: 0.7,
        diffusivity: 0.3,
    })
    .unwrap();
    s.set_source(1.3).unwrap();
    s.set_dirichlet(1.0, 2.0, 0.0, 0.5).unwrap();
    s.iterate(45).unwrap();
    compare(&g, "diffusion2d", s.field(), 1e-13, 1e-12);
}
#[test]
fn all_lbm_boundary_scenarios_match_cpp() {
    let g = goldens();
    for scenario in 0..4 {
        let mut c = D2Q9Config {
            nx: 12,
            ny: 10,
            tau: 0.73,
            ..Default::default()
        };
        if scenario == 1 {
            c.acceleration_x = 1e-5;
        }
        let mut s = D2Q9Solver::new(c).unwrap();
        if scenario == 0 {
            s.initialize_taylor_green(0.023).unwrap();
        } else {
            for x in 0..12 {
                s.set_solid(x, 0, true).unwrap();
                s.set_wall_velocity(x, 9, 0.017, 0.0).unwrap();
            }
            s.initialize_uniform(1.0, 0.0, 0.0).unwrap();
        }
        if scenario == 2 {
            s.set_velocity_inlet_left(0.01, 0.0).unwrap();
            s.set_pressure_outlet_right(1.0, 0.0).unwrap();
        }
        if scenario == 3 {
            for x in 0..12 {
                s.set_interpolated_wall_link(x, 1, 2, 0.3).unwrap();
                s.set_interpolated_wall_link(x, 8, 4, 0.7).unwrap();
            }
        }
        s.step(17).unwrap();
        let out: Vec<_> = (0..120)
            .flat_map(|i| {
                [
                    s.density()[i] as f64,
                    s.velocity_x()[i] as f64,
                    s.velocity_y()[i] as f64,
                ]
            })
            .collect();
        compare(&g, &format!("lbm_{scenario}"), &out, 8e-7, 8e-7);
    }
}
#[test]
fn yee_absorption_ade_tfsf_and_rlc_match_cpp() {
    let g = goldens();
    for scenario in 0..9 {
        let boundary = match scenario {
            1 => Boundary1D::Mur1,
            2 => Boundary1D::Pml,
            3 | 7 => Boundary1D::Cpml,
            _ => Boundary1D::Pec,
        };
        let mut s = Maxwell1D::new(Maxwell1DConfig {
            cells: 96,
            courant: 0.7,
            pml_cells: 10,
            boundary,
            ..Default::default()
        })
        .unwrap();
        match scenario {
            4 => s
                .set_debye_material(38..62, 2.1, 0.4, 2e-11, 0.001)
                .unwrap(),
            5 => s
                .set_drude_material(38..62, 2.1, 1.4e10, 4e8, 0.001)
                .unwrap(),
            6 => s
                .set_lorentz_material(38..62, 2.1, 0.4, 1.4e10, 4e8, 0.001)
                .unwrap(),
            _ => (),
        };
        if scenario == 7 {
            s.set_tfsf_source(
                22,
                |t| {
                    let x = (t - 8e-11) / 2e-11;
                    (-0.5 * x * x).exp()
                },
                0.0,
            )
            .unwrap();
        } else {
            s.initialize_gaussian(0.3, 0.045).unwrap();
        }
        if scenario == 8 {
            s.set_parallel_lumped_rlc(45, 1e-3, 1e-6, 100.0, 1e-9, 1e-13)
                .unwrap();
        }
        s.step(105).unwrap();
        let mut out = s.electric().to_vec();
        out.extend_from_slice(s.magnetic());
        out.extend_from_slice(s.polarization());
        out.push(s.energy());
        compare(&g, &format!("maxwell_{scenario}"), &out, 2e-13, 2e-11);
    }
}
#[test]
fn muscl_and_characteristic_weno_with_time_integrators_match_cpp() {
    let g = goldens();
    for scenario in 0..4 {
        let mut s = CompressibleEuler1D::new(Compressible1DConfig {
            cells: 40,
            periodic: scenario >= 2,
            reconstruction: if scenario % 2 == 1 {
                CompressibleReconstruction::CharacteristicWeno5
            } else {
                CompressibleReconstruction::MusclMinmod
            },
            time_integrator: if scenario % 2 == 1 {
                CompressibleTimeIntegrator::Ssprk3
            } else {
                CompressibleTimeIntegrator::ForwardEuler
            },
            ..Default::default()
        })
        .unwrap();
        if scenario < 2 {
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
        } else {
            let q: Vec<_> = (0..40)
                .map(|i| EulerPrimitive1D {
                    density: 1.0
                        + 0.1 * (2.0 * std::f64::consts::PI * (i as f64 + 0.5) / 40.0).sin(),
                    velocity: 0.5,
                    pressure: 1.0,
                })
                .collect();
            s.initialize_profile(&q).unwrap();
        }
        for _ in 0..9 {
            s.step(0.001).unwrap();
        }
        let out: Vec<_> = s
            .conserved()
            .iter()
            .flat_map(|u| [u.density, u.momentum, u.energy])
            .collect();
        compare(&g, &format!("euler_{scenario}"), &out, 3e-13, 3e-12);
    }
}
#[test]
fn triangle_tetrahedron_and_heat_solutions_match_cpp() {
    let g = goldens();
    let mesh = make_rectangle_tri_mesh(5, 4, 1.3, 0.8).unwrap();
    let s = poisson2d(
        &mesh,
        1.7,
        |x| 2.0 + x.x * x.y,
        |x| 0.4 + x.x - 0.2 * x.y,
        LinearConfig::default(),
    )
    .unwrap();
    compare(&g, "poisson2d", &s.values, 1e-9, 1e-9);
    let mesh = make_box_tet_mesh(3, 3, 3, 1.0, 1.0, 1.0).unwrap();
    let s = poisson3d(
        &mesh,
        1.0,
        |x| 1.0 + x.x,
        |x| 0.2 + x.x - x.y,
        LinearConfig::default(),
    )
    .unwrap();
    compare(&g, "poisson3d", &s.values, 1e-9, 1e-9);
    let mesh = make_rectangle_tri_mesh(5, 5, 1.0, 1.0).unwrap();
    let mut s = Heat2D::new(
        mesh,
        Heat2DConfig {
            conductivity: 0.8,
            volumetric_heat_capacity: 4.0,
            dt: 0.002,
            ..Default::default()
        },
    )
    .unwrap();
    s.initialize(|p| (std::f64::consts::PI * p.x).sin() * (std::f64::consts::PI * p.y).sin())
        .unwrap();
    for _ in 0..4 {
        s.step(|_, _, _| 0.3, |_, _| 0.0).unwrap();
    }
    compare(&g, "heat2d", &s.temperature, 1e-9, 1e-9);
}
#[test]
fn rf_network_and_guided_wave_models_match_cpp() {
    use cfd_solvers::{rf::*, Complex64 as C};
    let g = goldens();
    let line = TransmissionLine {
        characteristic_impedance: 63.0,
        propagation_constant: C::new(0.07, 15.2),
        length: 0.08,
    };
    let s = line.s_parameters(50.0).unwrap();
    let mut out = vec![];
    for m in [
        s,
        s_to_z(s, 50.0).unwrap(),
        s_to_abcd(s, 50.0).unwrap(),
        cascade_s(s, s, 50.0).unwrap(),
    ] {
        for z in m.values() {
            out.extend([z.re, z.im]);
        }
    }
    let m = microstrip_quasi_static(0.001, 0.0008, 4.2, 1.7e9).unwrap();
    out.extend([
        m.characteristic_impedance,
        m.effective_permittivity,
        m.phase_velocity,
        m.guided_wavelength,
    ]);
    let c = coax_quasi_static(0.0005, 0.0017, 2.1, 1.7e9).unwrap();
    out.extend([
        c.characteristic_impedance,
        c.phase_velocity,
        c.propagation_constant_rad_per_m,
        c.guided_wavelength,
    ]);
    let w = rectangular_waveguide_te10(0.02286, 0.01016, 1.0, 1.0, 10e9).unwrap();
    out.extend([
        w.cutoff_frequency_hz,
        w.propagation_constant_rad_per_m,
        w.guide_wavelength,
        w.wave_impedance,
    ]);
    let cp = coplanar_waveguide_quasi_static(0.001, 0.0002, 4.2, 1.7e9).unwrap();
    out.extend([
        cp.characteristic_impedance,
        cp.effective_permittivity,
        cp.phase_velocity,
        cp.guided_wavelength,
    ]);
    compare(&g, "rf_network", &out, 1e-12, 1e-12);
}
#[test]
fn geometric_polarization_material_and_gaussian_optics_match_cpp() {
    use cfd_solvers::{core::Vec3, optics::*};
    let g = goldens();
    let mut b = GaussianBeam::new(532e-9, 0.0007, 0.012, 1.0).unwrap();
    b.propagate(0.03).unwrap();
    b.thin_lens(0.08).unwrap();
    let mut out = vec![b.q().re, b.q().im, b.radius(), b.curvature_radius()];
    let materials = MaterialCatalog::default();
    out.push(
        materials
            .at("N-BK7")
            .unwrap()
            .refractive_index_nm(532.0)
            .unwrap(),
    );
    let f = fresnel_dielectric(1.0, 1.52, 0.7).unwrap();
    out.extend([
        f.cos_transmitted,
        f.rs,
        f.rp,
        f.ts,
        f.tp,
        f.reflectance_s,
        f.reflectance_p,
    ]);
    out.push(
        thin_film_reflectance_normal(
            1.0,
            1.52,
            532.0,
            &[
                ThinFilmLayer {
                    refractive_index: 1.38,
                    thickness_nm: 93.0,
                },
                ThinFilmLayer {
                    refractive_index: 2.1,
                    thickness_nm: 63.0,
                },
            ],
        )
        .unwrap(),
    );
    let mut s = SequentialOpticalSystem::new(1.0).unwrap();
    s.add_surface(SequentialSurface {
        surface_type: SurfaceType::Sphere,
        vertex_z: 0.01,
        radius: 0.05,
        aperture_radius: 0.005,
        refractive_index_after: 1.5,
        ..Default::default()
    })
    .unwrap();
    for x in [-0.002, 0.0, 0.002] {
        let r = s
            .trace(Ray {
                origin: Vec3::new(x, 0.0, 0.0),
                wavelength_nm: 532.0,
                ..Default::default()
            })
            .unwrap();
        out.extend([
            r.ray.origin.x,
            r.ray.origin.z,
            r.ray.direction.x,
            r.ray.direction.z,
            if r.valid { 1.0 } else { 0.0 },
        ]);
    }
    compare(&g, "optics", &out, 2e-13, 2e-12);
}
#[test]
fn kinetic_reactor_and_thermodynamic_equilibrium_match_cpp() {
    use cfd_solvers::chemistry::*;
    let g = goldens();
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
    let mut out = n.source_terms(&[2.0, 0.0], 300.0).unwrap();
    let mut c = [1.0, 0.0];
    let r = integrate_isothermal(&n, &mut c, 300.0, 0.2, ReactorConfig::default()).unwrap();
    out.extend(c);
    out.push(r.time_advanced);
    let phase = IdealGasPhase::new(
        vec![
            ThermoSpecies {
                name: "A".into(),
                standard_chemical_potential: 10.0,
            },
            ThermoSpecies {
                name: "B".into(),
                standard_chemical_potential: 20.0,
            },
        ],
        1e5,
    )
    .unwrap();
    out.extend(phase.chemical_potentials(&[0.3, 0.7], 310.0, 2e5).unwrap());
    let p = equilibrate_binary_salt(3.0, 2.0, 0.0, 0.25).unwrap();
    out.extend([
        p.dissolved_a,
        p.dissolved_b,
        p.solid_amount,
        p.reaction_extent,
        p.saturation_ratio,
    ]);
    compare(&g, "chemistry", &out, 2e-9, 2e-11);
}
#[test]
fn nonlinear_lossy_acoustics_sensors_and_time_reversal_match_cpp() {
    use cfd_solvers::acoustics::*;
    let g = goldens();
    let mut s = KSpaceAcoustic2D::new(KSpace2DConfig {
        nx: 16,
        ny: 16,
        pml_cells: 2,
        ..Default::default()
    })
    .unwrap();
    s.set_uniform_medium(1200.0, 998.0).unwrap();
    s.set_uniform_nonlinearity(3.0).unwrap();
    s.set_power_law_absorption(PowerLawAbsorption2D {
        alpha0_np_per_m: 0.3,
        exponent: 1.5,
        reference_frequency_hz: 1e6,
    })
    .unwrap();
    s.add_sensor(136).unwrap();
    s.add_sensor(140).unwrap();
    s.initialize_gaussian_pressure(0.5, 0.5, 0.0002, 700.0)
        .unwrap();
    s.add_pressure_source(100, |t| 2.0 * (1e6 * t).sin())
        .unwrap();
    s.step(4).unwrap();
    let mut out = s.pressure().to_vec();
    out.extend_from_slice(s.velocity_x());
    out.extend_from_slice(s.velocity_y());
    for sensor in s.sensors() {
        out.extend_from_slice(&sensor.pressure);
    }
    compare(&g, "acoustics", &out, 2e-8, 2e-10);
    let traces = vec![
        s.sensors()[0].pressure.clone(),
        s.sensors()[1].pressure.clone(),
    ];
    let reverse = s
        .time_reversal_reconstruct(&[136, 140], &traces, true)
        .unwrap();
    compare(&g, "acoustic_time_reversal", &reverse, 2e-8, 2e-10);
}
