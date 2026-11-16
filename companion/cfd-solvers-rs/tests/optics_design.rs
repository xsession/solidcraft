//! Analytic and robustness checks for optics analysis, interchange, non-sequential tracing and
//! optimization.
use cfd_solvers::{core::Vec3, optics::*, Complex64 as C};
use std::f64::consts::PI;

#[test]
fn zernike_fit_recovers_exact_modes_and_retarder_acts_as_quarter_wave() {
    let ws: Vec<WavefrontSample> = (0..60)
        .map(|i| {
            let rho = ((f64::from(i) + 0.5) / 60.0).sqrt();
            let th = 2.4 * f64::from(i);
            WavefrontSample {
                rho,
                theta: th,
                opd: 0.5 - 0.2 * rho * th.sin() + 0.3 * (2.0 * rho * rho - 1.0),
            }
        })
        .collect();
    let fit = fit_zernike(&ws, 2).unwrap();
    let get = |n, m| fit.iter().find(|c| c.n == n && c.m == m).unwrap().value;
    assert!(
        (get(0, 0) - 0.5).abs() < 1e-12
            && (get(1, -1) + 0.2).abs() < 1e-12
            && (get(2, 0) - 0.3).abs() < 1e-12
    );
    assert!(get(2, 2).abs() < 1e-12);
    assert!(fit_zernike(&ws[..3], 2).is_err());
    assert!(zernike(2, 0, 1.5, 0.0).is_err());

    // 45-degree linear light through a quarter-wave plate becomes circular (S3 = +-1).
    let j = apply_linear_retarder(
        JonesVector {
            x: C::new(1.0 / 2f64.sqrt(), 0.0),
            y: C::new(1.0 / 2f64.sqrt(), 0.0),
        },
        PI / 2.0,
        0.0,
    );
    let s = stokes_from_jones(j);
    assert!((s.v.abs() - s.i).abs() < 1e-12);
    // Along the optic axis the extraordinary wave sees n_o; perpendicular to it, n_e.
    assert!((uniaxial_extraordinary_index(1.6, 1.5, 0.0).unwrap() - 1.6).abs() < 1e-15);
    assert!((uniaxial_extraordinary_index(1.6, 1.5, PI / 2.0).unwrap() - 1.5).abs() < 1e-12);
}

#[test]
fn diffraction_psf_sums_to_one_and_mtf_starts_at_one() {
    let pupil = vec![C::new(1.0, 0.0); 8];
    let psf = normalized_psf_1d(&pupil).unwrap();
    assert!((psf.iter().sum::<f64>() - 1.0).abs() < 1e-14);
    // A uniform pupil puts all energy in the DC bin.
    assert!((psf[0] - 1.0).abs() < 1e-12);
    let mtf = mtf_from_psf_1d(&psf).unwrap();
    assert!((mtf[0] - 1.0).abs() < 1e-15);
    let ff = fraunhofer_1d(&pupil, 1e-3, 1e-6, 1.0).unwrap();
    assert_eq!(ff.iter().map(|s| s.intensity).fold(0.0, f64::max), 1.0);
    assert!(fraunhofer_1d(&pupil[..1], 1e-3, 1e-6, 1.0).is_err());
}

#[test]
fn kerr_response_reduces_to_linear_permittivity_without_chi3() {
    let m = NonlinearOpticalMaterial {
        relative_permittivity: 4.0,
        ..Default::default()
    };
    let r =
        evaluate_nonlinear_polarization(C::new(1e5, -2e5), &m, FREE_SPACE_WAVE_IMPEDANCE).unwrap();
    assert!((r.effective_relative_permittivity - C::new(4.0, 0.0)).norm() < 1e-12);
    assert!(evaluate_nonlinear_polarization(C::new(f64::NAN, 0.0), &m, 377.0).is_err());
    assert!(evaluate_nonlinear_polarization(C::new(1.0, 0.0), &m, 0.0).is_err());
}

#[test]
fn rigid_transform_round_trips_and_preserves_length() {
    let t = RigidTransform::from_euler_xyz(Vec3::new(3.0, -1.0, 2.0), 0.7, -1.1, 2.3);
    let p = Vec3::new(0.3, -0.4, 1.7);
    let q = t.inverse_point(t.apply_point(p));
    assert!((q - p).norm() < 1e-14);
    assert!((t.apply_vector(p).norm() - p.norm()).abs() < 1e-14);
}

#[test]
fn freeform_plane_intersection_and_misses() {
    let flat = PolynomialFreeformSurface::new(2.0, vec![], 1.0).unwrap();
    let r = Ray {
        origin: Vec3::new(0.1, 0.2, 0.0),
        ..Default::default()
    };
    let (p, n) = flat.intersect(&r, 40, 1e-12).unwrap();
    assert!((p.z - 2.0).abs() < 1e-15 && (n.z + 1.0).abs() < 1e-15);
    let outside = Ray {
        origin: Vec3::new(5.0, 0.0, 0.0),
        ..Default::default()
    };
    assert!(flat.intersect(&outside, 40, 1e-12).is_none());
    let behind = Ray {
        origin: Vec3::new(0.0, 0.0, 5.0),
        ..Default::default()
    };
    assert!(flat.intersect(&behind, 40, 1e-12).is_none());
    assert!(PolynomialFreeformSurface::new(
        0.0,
        vec![PolynomialTerm {
            x_power: 1,
            y_power: 0,
            coefficient: f64::NAN
        }],
        1.0
    )
    .is_err());
}

#[test]
fn json_and_interchange_round_trip_and_reject_garbage() {
    let mut s = SequentialOpticalSystem::new(1.2).unwrap();
    s.add_surface(SequentialSurface {
        surface_type: SurfaceType::Conic,
        vertex_z: 1.0,
        radius: 20.0,
        conic_constant: -1.0,
        aperture_radius: 4.0,
        refractive_index_after: 1.6,
        even_coefficients: [1e-3, 0.0, 0.0, 0.0],
    })
    .unwrap();
    let back = sequential_system_from_json(&sequential_system_to_json(&s, 1.2)).unwrap();
    assert_eq!(back.surfaces()[0].radius, 20.0);
    assert_eq!(back.object_space_index(), 1.2);
    assert!(sequential_system_from_json("{}").is_err());
    assert!(sequential_system_from_json("{\"object_index\":1,\"surfaces\":[{\"type\":\"torus\",\"z\":0,\"radius\":0,\"aperture\":1,\"n_after\":1,\"k\":0,\"a\":[0,0,0,0]}]}").is_err());
    assert!(import_lens_system("nothing here", LensInterchangeFormat::ZemaxZmx).is_err());
    assert!(import_lens_system("SURF\nCURV abc\n", LensInterchangeFormat::ZemaxZmx).is_err());
    assert!(import_lens_system("SURF\nPARM x 1\n", LensInterchangeFormat::ZemaxZmx).is_err());
    // Overflowing tokens are skipped, as std::stod's out_of_range is caught upstream.
    let s = import_lens_system("SURF\nDISZ 1e999 2.5\n", LensInterchangeFormat::ZemaxZmx).unwrap();
    assert_eq!(s.surfaces()[0].vertex_z, 2.5);
    for fmt in [
        LensInterchangeFormat::ZemaxZmx,
        LensInterchangeFormat::CodevSeq,
        LensInterchangeFormat::OsloLen,
    ] {
        let t = import_lens_system(&export_lens_system(&back, fmt), fmt).unwrap();
        let (a, b) = (t.surfaces()[0], back.surfaces()[0]);
        assert_eq!(
            (a.surface_type, a.vertex_z, a.radius, a.even_coefficients),
            (b.surface_type, b.vertex_z, b.radius, b.even_coefficients)
        );
    }
}

#[test]
fn non_sequential_power_is_bounded_by_launch() {
    let mut scene = NonSequentialScene::default();
    scene
        .add_plane(NsPlane {
            point: Vec3::new(0.0, 0.0, 1.0),
            reflectance: 0.3,
            transmittance: 0.6,
            ..Default::default()
        })
        .unwrap();
    let out = scene.trace(NsRay::default(), 4, 1e-9).unwrap();
    let s = summarize_stray_light(&out, 1.0).unwrap();
    assert!((s.total_power - 0.9).abs() < 1e-12 && (s.absorbed_power - 0.1).abs() < 1e-12);
    assert!(scene
        .add_plane(NsPlane {
            reflectance: 0.7,
            transmittance: 0.6,
            ..Default::default()
        })
        .is_err());
    assert!(scene
        .add_plane(NsPlane {
            normal: Vec3::ZERO,
            ..Default::default()
        })
        .is_err());
    let ghosts = rank_ghost_paths(&out, 1);
    assert!(ghosts.windows(2).all(|w| w[0].power >= w[1].power));
}

#[test]
fn optimizers_find_quadratic_minimum() {
    let f = |x: &[f64]| (x[0] - 1.5).powi(2) + 2.0 * (x[1] + 0.5).powi(2);
    let cd = coordinate_descent(&f, vec![0.0, 0.0], vec![1.0, 1.0], 500, 1e-10).unwrap();
    assert!(
        cd.converged
            && (cd.variables[0] - 1.5).abs() < 1e-8
            && (cd.variables[1] + 0.5).abs() < 1e-8
    );
    let de = differential_evolution(&f, &[-3.0, -3.0], &[3.0, 3.0], 20, 150, 0.8, 0.9, 5).unwrap();
    assert!(de.merit < 1e-6);
    let g = finite_difference_sensitivity(&f, &[0.0, 0.0], 1e-6).unwrap();
    assert!((g[0] + 3.0).abs() < 1e-6 && (g[1] - 2.0).abs() < 1e-6);
    let mc = monte_carlo_tolerance(&f, &[1.5, -0.5], &[0.0, 0.0], 10, 1).unwrap();
    assert_eq!(mc.standard_deviation, 0.0);
    assert!(differential_evolution(&f, &[1.0], &[0.0], 20, 10, 0.8, 0.9, 1).is_err());
    assert!(coordinate_descent(&f, vec![], vec![], 10, 1e-6).is_err());
}
