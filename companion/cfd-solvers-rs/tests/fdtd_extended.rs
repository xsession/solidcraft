//! Analytic and robustness checks for 3-D and axisymmetric FDTD, monitors, ports, NF2FF, SAR and
//! EMC helpers.
use cfd_solvers::{core::Vec3, fdtd::*, Complex64 as C, SolverError};
use std::f64::consts::PI;

#[test]
fn pec_cavity_energy_stays_bounded_and_walls_stay_conducting() {
    let mut m = Maxwell3D::new(Maxwell3DConfig {
        nx: 11,
        ny: 11,
        nz: 11,
        ..Default::default()
    })
    .unwrap();
    m.initialize_gaussian_ez(1.0, 0.15).unwrap();
    let e0 = m.energy();
    let mut lo = f64::INFINITY;
    let mut hi = 0.0_f64;
    for _ in 0..60 {
        m.step(5);
        lo = lo.min(m.energy());
        hi = hi.max(m.energy());
    }
    // Leapfrog energy oscillates (E and H are staggered in time) but has no secular growth.
    assert!(hi < 1.6 * e0 && lo > 0.3 * e0, "{lo} {hi} {e0}");
    // PEC walls: every E component vanishes on the boundary layer.
    let idx = |i: usize, j: usize, k: usize| (k * 11 + j) * 11 + i;
    for a in 0..11 {
        for b in 0..11 {
            for n in [
                idx(0, a, b),
                idx(10, a, b),
                idx(a, 0, b),
                idx(a, 10, b),
                idx(a, b, 0),
                idx(a, b, 10),
            ] {
                assert_eq!((m.ex()[n], m.ey()[n], m.ez()[n]), (0.0, 0.0, 0.0));
            }
        }
    }
    assert!(m.max_field() > 0.0 && m.max_field().is_finite());
    assert!(Maxwell3D::new(Maxwell3DConfig {
        courant: 1.0,
        ..Default::default()
    })
    .is_err());
    assert!(Maxwell3D::new(Maxwell3DConfig {
        nx: 3,
        ..Default::default()
    })
    .is_err());
    assert!(m.add_soft_ez_source(11, 0, 0, 1.0).is_err());
    assert!(matches!(
        write_maxwell3d_hdf5_xdmf(&m, "/tmp/a.xdmf", "/tmp/a.h5"),
        Err(SolverError::Unsupported(_))
    ));
    assert!(write_maxwell3d_vtk_ascii(&m, "/nonexistent-dir/x.vtk").is_err());
}

#[test]
fn cylindrical_cavity_is_stable_and_validates_indices() {
    let mut c = CylindricalTM::new(CylindricalTMConfig {
        nr: 20,
        nz: 24,
        ..Default::default()
    })
    .unwrap();
    c.initialize_gaussian_ez(0.0, 12e-3, 3e-3, 1.0).unwrap();
    let e0 = c.energy();
    c.step(400);
    let e1 = c.energy();
    assert!(
        e1.is_finite() && e1 < 1.6 * e0 && e1 > 0.3 * e0,
        "{e0} {e1}"
    );
    assert!(c.add_soft_ez_source(0, 23, 1.0).is_err());
    assert!(CylindricalTM::new(CylindricalTMConfig {
        courant: 0.6,
        ..Default::default()
    })
    .is_err());
    let lv = make_cylindrical_multigrid_levels(3, 64, 1.0, 1.0, 10).unwrap();
    assert_eq!(lv.len(), 2);
    assert_eq!(lv[1].radial_cells, 2);
}

#[test]
fn dft_monitor_recovers_sinusoid_amplitude_and_matched_port_has_no_reflection() {
    let f = 1e9;
    let mut d = DftMonitor::new(f).unwrap();
    let n = 4000;
    let duration = 10.0 / f;
    for i in 0..=n {
        let t = duration * f64::from(i) / f64::from(n);
        d.sample(t, 0.8 * (2.0 * PI * f * t + 0.4).cos()).unwrap();
    }
    assert!((d.amplitude(duration).unwrap() - 0.8).abs() < 1e-4);
    assert!(d.sample(0.0, 1.0).is_err());

    let mut a = WavePort1D::new(f, 50.0, 1).unwrap();
    let mut b = WavePort1D::new(f, 50.0, 1).unwrap();
    for i in 0..=n {
        let t = duration * f64::from(i) / f64::from(n);
        let e = (2.0 * PI * f * t).sin();
        a.sample(t, e, -e / 50.0).unwrap();
        b.sample(t, e, -e / 50.0).unwrap();
    }
    let s = s_parameters(&a, &b).unwrap();
    assert!(s.s11.norm() < 1e-12 && (s.s21 - C::new(1.0, 0.0)).norm() < 1e-12);
    assert!(WavePort1D::new(f, 50.0, 0).is_err());
    let empty = WavePort1D::new(f, 50.0, 1).unwrap();
    assert!(s_parameters(&empty, &b).is_err());
}

#[test]
fn sar_of_uniform_tissue_is_uniform_and_mass_averaging_preserves_it() {
    let n = 4 * 4 * 4;
    let g = SarGrid3D {
        nx: 4,
        ny: 4,
        nz: 4,
        dx: 4e-3,
        dy: 4e-3,
        dz: 4e-3,
        ex: vec![10.0; n],
        ey: vec![0.0; n],
        ez: vec![0.0; n],
        conductivity: vec![0.8; n],
        density: vec![1000.0; n],
    };
    let expected = 0.8 * 100.0 / 1000.0;
    for v in local_sar(&g)
        .unwrap()
        .iter()
        .chain(&sar_1g(&g).unwrap())
        .chain(&sar_10g(&g).unwrap())
    {
        assert!((v - expected).abs() < 1e-14);
    }
    let bad = SarGrid3D {
        density: vec![0.0; n],
        ..g.clone()
    };
    assert!(local_sar(&bad).is_err());
    assert!(mass_averaged_sar(&g, 0.0).is_err());
}

#[test]
fn nf2ff_of_a_single_patch_is_finite_and_scales_with_area() {
    let s = NearFieldSurfaceSample {
        normal: Vec3::new(0.0, 0.0, 1.0),
        area: 1e-4,
        magnetic: ComplexVec3 {
            x: C::new(1.0, 0.0),
            ..Default::default()
        },
        ..Default::default()
    };
    let a = nf2ff(&[s], Vec3::new(0.0, 0.0, 1.0), 1e9, NF2FF_WAVE_IMPEDANCE).unwrap();
    let b = nf2ff(
        &[NearFieldSurfaceSample { area: 2e-4, ..s }],
        Vec3::new(0.0, 0.0, 1.0),
        1e9,
        NF2FF_WAVE_IMPEDANCE,
    )
    .unwrap();
    assert!((b.radiation_intensity / a.radiation_intensity - 4.0).abs() < 1e-12);
    assert!(nf2ff(&[s], Vec3::ZERO, 1e9, NF2FF_WAVE_IMPEDANCE).is_err());
    assert!(nf2ff(&[], Vec3::new(0.0, 0.0, 1.0), 1e9, NF2FF_WAVE_IMPEDANCE).is_err());
}

#[test]
fn emc_helpers_satisfy_identities() {
    let x = make_graded_axis(1.0, 8, 0.05, 1.2).unwrap();
    assert_eq!((x[0], x[8]), (0.0, 1.0));
    assert!(x.windows(3).all(|w| w[2] - w[1] > w[1] - w[0]));
    let w = deposit_thin_wire_subcell(1.3e-3, 2.2e-3, 0.4e-3, 5, 5, 5, 1e-3, 1e-3, 1e-3).unwrap();
    assert!((w.iter().map(|w| w.weight).sum::<f64>() - 1.0).abs() < 1e-12);
    let h = huygens_equivalent_source([0.0, 0.0, 3.0], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]).unwrap();
    assert_eq!(h.electric_surface_current[2], 0.0);
    assert!(huygens_equivalent_source([0.0; 3], [0.0; 3], [0.0; 3]).is_err());
    // A vanishing sheet between equal media transmits fully.
    let t = impedance_sheet_transmission(C::new(0.0, 0.0), C::new(377.0, 0.0), C::new(377.0, 0.0))
        .unwrap();
    assert!((t - C::new(1.0, 0.0)).norm() < 1e-15);
    let bad = [
        EmcTransferSample {
            frequency_hz: 2.0,
            gain: C::new(1.0, 0.0),
        },
        EmcTransferSample {
            frequency_hz: 1.0,
            gain: C::new(0.0, 0.0),
        },
        EmcTransferSample {
            frequency_hz: 3.0,
            gain: C::new(0.0, 0.0),
        },
    ];
    assert!(interpolate_emc_transfer(&bad, 2.5).is_err());
    assert!(interpolate_emc_transfer(&[], 1.0).is_err());
    // An LC line with zero source conserves 0.5 C V^2 + 0.5 L I^2 to first order in dt.
    let cfg = MulticonductorTlineConfig {
        conductors: 1,
        inductance_h: vec![1e-6],
        capacitance_f: vec![1e-9],
        dt_s: 1e-10,
    };
    let mut st = MulticonductorTlineState {
        voltage: vec![1.0],
        current: vec![0.0],
    };
    let energy = |s: &MulticonductorTlineState| {
        0.5 * 1e-9 * s.voltage[0] * s.voltage[0] + 0.5 * 1e-6 * s.current[0] * s.current[0]
    };
    let e0 = energy(&st);
    for _ in 0..1000 {
        advance_multiconductor_tline(&cfg, &mut st, &[0.0]).unwrap();
    }
    assert!((energy(&st) / e0 - 1.0).abs() < 0.05);
    assert!(advance_multiconductor_tline(&cfg, &mut st, &[0.0, 1.0]).is_err());
}
