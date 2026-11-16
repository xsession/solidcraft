//! Analytic and robustness checks for the `em` module.
use cfd_solvers::{em::*, fem, rf::ComplexMatrix, Complex64 as C};

const C0: f64 = 299_792_458.0;

#[test]
fn waveguide_modes_have_textbook_cutoffs_and_unit_power() {
    let (a, b) = (22.86e-3, 10.16e-3);
    let modes = rectangular_waveguide_modes(a, b, 10e9, 1.0, 1.0, 3).unwrap();
    assert_eq!(
        (modes[0].family, modes[0].m, modes[0].n),
        (ModeFamily::Te, 1, 0)
    );
    // c / 2a with the eps0/mu0 constants used upstream.
    assert!((modes[0].cutoff_frequency_hz / (C0 / (2.0 * a)) - 1.0).abs() < 1e-6);
    assert!(modes[0].propagating && !modes[2].propagating);
    let p = rectangular_mode_power_w(&modes[0], a, b, 200, 50).unwrap();
    assert!((p - 1.0).abs() < 1e-4, "{p}");
    assert!(rectangular_waveguide_transverse_fields(&modes[0], a, b, -1e-3, 0.0).is_err());
}

#[test]
fn gstc_limits_and_power_balance() {
    let t = solve_gstc_normal_incidence(&GstcSheetConfig::default()).unwrap();
    assert!(t.reflection.norm() < 1e-15 && (t.transmission - C::new(1.0, 0.0)).norm() < 1e-15);
    // A purely reactive sheet between equal media is lossless.
    let r = solve_gstc_normal_incidence(&GstcSheetConfig {
        electric_surface_admittance_s: C::new(0.0, 0.02),
        magnetic_surface_impedance_ohm: C::new(0.0, 40.0),
        ..Default::default()
    })
    .unwrap();
    assert!(r.absorbed_power.abs() < 1e-12, "{}", r.absorbed_power);
    assert!(solve_gstc_normal_incidence(&GstcSheetConfig {
        wave_impedance_after_ohm: 0.0,
        ..Default::default()
    })
    .is_err());
}

#[test]
fn one_dimensional_cavity_and_driven_solutions() {
    let modes = pec_cavity_eigenmodes_1d_uniform(61, 1.0, 4.0, 1.0, 2).unwrap();
    let f1 = C0 / (2.0 * 1.0 * 2.0);
    assert!((modes[0].frequency_hz / f1 - 1.0).abs() < 1e-3);
    assert!((modes[1].frequency_hz / (2.0 * f1) - 1.0).abs() < 1e-3);
    assert!(pec_cavity_eigenmodes_1d_uniform(5, 1.0, 1.0, 1.0, 4).is_err());
    let cfg = FrequencyDomain1DConfig {
        points: 51,
        ..Default::default()
    };
    let j = vec![C::new(1.0, 0.0); 51];
    let r = solve_pec_driven_maxwell_1d(&cfg, &j).unwrap();
    assert_eq!(r.electric_v_per_m[0], C::new(0.0, 0.0));
    assert_eq!(r.electric_v_per_m[50], C::new(0.0, 0.0));
    assert!(r.electric_v_per_m[25].norm() > 0.0);
    assert!(solve_pec_driven_maxwell_1d(&cfg, &j[1..]).is_err());
}

#[test]
fn matched_lossless_line_preserves_amplitude() {
    let (l, c): (f64, f64) = (2.5e-7, 1e-10);
    let z0 = (l / c).sqrt();
    let rlcg = MulticonductorRlcg {
        conductors: 1,
        resistance: vec![0.0],
        inductance: vec![l],
        conductance: vec![0.0],
        capacitance: vec![c],
    };
    let load = ComplexMatrix::from_row_major(vec![C::new(z0, 0.0)]).unwrap();
    let r = solve_multiconductor_cable(2.0, 400, 5e7, &rlcg, &[C::new(1.0, 0.0)], &load).unwrap();
    assert!((r.load_voltage_v[0].norm() - 1.0).abs() < 1e-9);
    assert!((r.input_current_a[0] - C::new(1.0 / z0, 0.0)).norm() < 1e-9);
    assert!(solve_multiconductor_cable(-1.0, 4, 1e6, &rlcg, &[C::new(1.0, 0.0)], &load).is_err());
}

#[test]
fn signal_and_power_integrity_identities() {
    let bits = [0, 1, 0, 1, 1, 0];
    let wave: Vec<f64> = bits.iter().flat_map(|b| [f64::from(*b); 4]).collect();
    let e = analyze_nrz_eye(&wave, &bits, 4, 0).unwrap();
    assert_eq!(e.eye_width_ui, 1.0);
    assert!(e.estimated_ber < 1e-300 && (e.eye_height - 1.0).abs() < 1e-15);
    assert!(analyze_nrz_eye(&wave, &[0, 2, 1, 0, 1, 0], 4, 0).is_err());
    let cfg = PdnGridConfig {
        nx: 3,
        ny: 3,
        ..Default::default()
    };
    let ir = solve_pdn_ir_drop_grid(
        &cfg,
        &[0.0; 9],
        &[PdnVoltageSource {
            node: 4,
            voltage_v: 1.2,
        }],
    )
    .unwrap();
    assert!(ir.converged && ir.voltage_v.iter().all(|v| (v - 1.2).abs() < 1e-12));
    // Series resonance of a decoupling capacitor: |Z| = ESR.
    let cap = DecouplingCapacitor {
        capacitance_f: 1e-6,
        esr_ohm: 0.01,
        esl_h: 1e-9,
        count: 1,
    };
    let f0 = 1.0 / (2.0 * std::f64::consts::PI * (1e-6f64 * 1e-9).sqrt());
    assert!((decoupling_impedance_ohm(&cap, f0).unwrap().norm() - 0.01).abs() < 1e-9);
    let p = probe_waveform(&[1.0, -1.0, 1.0, -1.0], 0.5).unwrap();
    assert_eq!((p.peak, p.rms, p.impulse, p.energy), (1.0, 1.0, 0.0, 2.0));
    let peak = sample_double_exponential_pulse(100e-9, 0.01e-9, 3.0, 1e-9, 20e-9)
        .unwrap()
        .into_iter()
        .fold(0.0, f64::max);
    assert!((peak - 3.0).abs() < 1e-4);
}

#[test]
fn edge_fem_post_processing_and_refinement() {
    let mesh = fem::make_box_tet_mesh(2, 2, 2, 1.0, 1.0, 1.0).unwrap();
    let cfg = EdgeMaxwell3DConfig {
        frequency_hz: 1e8,
        ..Default::default()
    };
    let r = solve_driven_edge_maxwell_3d(
        &mesh,
        &cfg,
        None,
        &[EdgeCurrentPort3D {
            node_a: 13,
            node_b: 4,
            current_a: C::new(1.0, 0.0),
        }],
    )
    .unwrap();
    let q = resonator_quality_3d(&mesh, &cfg, &r, 0.0, 1.0).unwrap();
    assert_eq!(q.quality_factor, f64::INFINITY);
    assert!(q.electric_energy_j > 0.0 && q.magnetic_energy_j > 0.0);
    // Boundary (PEC) edges carry no voltage.
    for (v, bnd) in r.edge_voltage_v.iter().zip(&r.boundary_edge) {
        if *bnd {
            assert_eq!(*v, C::new(0.0, 0.0));
        }
    }
    assert!(edge_port_thevenin(&r, r.edges.len(), C::new(1.0, 0.0)).is_err());
    assert!(solve_driven_edge_maxwell_3d(
        &mesh,
        &cfg,
        None,
        &[EdgeCurrentPort3D {
            node_a: 0,
            node_b: 1,
            current_a: C::new(1.0, 0.0)
        }]
    )
    .is_err());
    // Physical conductivity damps the response much more than the upstream-compatible term.
    let lossy = EdgeMaxwell3DConfig {
        conductivity_s_per_m: 1.0,
        ..cfg
    };
    let physical = EdgeMaxwell3DConfig {
        physical_conductivity_term: true,
        ..lossy
    };
    let port = [EdgeCurrentPort3D {
        node_a: 13,
        node_b: 4,
        current_a: C::new(1.0, 0.0),
    }];
    let (a, b) = (
        solve_driven_edge_maxwell_3d(&mesh, &lossy, None, &port).unwrap(),
        solve_driven_edge_maxwell_3d(&mesh, &physical, None, &port).unwrap(),
    );
    let peak = |r: &EdgeMaxwell3DResult| {
        r.edge_voltage_v
            .iter()
            .fold(0.0_f64, |m, v| m.max(v.norm()))
    };
    assert!(peak(&b) < peak(&a));

    let marks: Vec<bool> = (0..mesh.element_count()).map(|i| i % 5 == 0).collect();
    let refined = refine_tet4_marked_longest_edges(&mesh, &marks).unwrap();
    let volume = |m: &fem::Mesh3D| {
        m.tetrahedra
            .iter()
            .map(|t| nedelec_tet4_geometry(m, t).unwrap().volume_m3)
            .sum::<f64>()
    };
    assert!(refined.element_count() > mesh.element_count());
    assert!((volume(&refined) - 1.0).abs() < 1e-12);
    assert!(mark_maxwell_dorfler(&[0.0, 0.0], 0.5)
        .unwrap()
        .iter()
        .all(|m| !m));
    assert!(mark_maxwell_dorfler(&[1.0, f64::NAN], 0.5).is_err());
}

#[test]
fn photonic_modes_reduce_to_free_space_and_guide_light() {
    let free = PhotonicGrid2D {
        nx: 6,
        ny: 6,
        width_m: 1e-6,
        height_m: 1e-6,
        relative_permittivity: vec![1.0; 36],
        dispersive_permittivity: None,
    };
    let m = solve_photonic_band_modes(
        &free,
        BlochWaveVector2D::default(),
        &PhotonicModeConfig::default(),
    )
    .unwrap();
    assert!(m[0].eigenvalue_per_m2.abs() < 1e-3 * 1e12 && m[0].converged);
    let eps: Vec<f64> = (0..64)
        .map(|i| {
            if (2..=5).contains(&(i % 8)) && (2..=5).contains(&(i / 8)) {
                12.0
            } else {
                2.0
            }
        })
        .collect();
    let wg = PhotonicGrid2D {
        nx: 8,
        ny: 8,
        width_m: 3e-6,
        height_m: 3e-6,
        relative_permittivity: eps,
        dispersive_permittivity: None,
    };
    let modes = solve_photonic_waveguide_modes(
        &wg,
        &PhotonicWaveguideConfig {
            frequency_hz: 2e14,
            mode_count: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(modes[0].guided);
    assert!(modes[0].effective_index > 2f64.sqrt() && modes[0].effective_index < 12f64.sqrt());
    assert!(modes[0].relative_residual < 1e-8);
    let big = PhotonicGrid2D {
        nx: 20,
        ny: 20,
        width_m: 1.0,
        height_m: 1.0,
        relative_permittivity: vec![1.0; 400],
        dispersive_permittivity: None,
    };
    assert!(solve_photonic_waveguide_modes(
        &big,
        &PhotonicWaveguideConfig {
            frequency_hz: 1e9,
            ..Default::default()
        }
    )
    .is_err());
}
