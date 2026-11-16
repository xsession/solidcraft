//! Value-by-value comparison with the C++ upstream optics analysis/design units. Goldens come
//! from `parity/oracle_optics.cpp`.
use cfd_solvers::{core::Dual, core::Vec3, optics::*, Complex64 as C};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-optics.json")).unwrap()
}

#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    let floor = 1e-6 * cpp.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= rel * c.abs().max(floor) + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn text(name: &str) -> String {
    goldens()[name].as_str().unwrap().to_string()
}

fn cx(c: C) -> [f64; 2] {
    [c.re, c.im]
}

fn lens() -> SequentialOpticalSystem {
    let mut s = SequentialOpticalSystem::new(1.0).unwrap();
    s.add_surface(SequentialSurface {
        surface_type: SurfaceType::Sphere,
        vertex_z: 0.0,
        radius: 50.0,
        aperture_radius: 10.0,
        refractive_index_after: 1.5,
        ..Default::default()
    })
    .unwrap();
    s.add_surface(SequentialSurface {
        surface_type: SurfaceType::EvenAsphere,
        vertex_z: 5.0,
        radius: -50.0,
        aperture_radius: 10.0,
        refractive_index_after: 1.0,
        conic_constant: -0.5,
        even_coefficients: [1e-5, -2e-8, 0.0, 3e-12],
    })
    .unwrap();
    s.add_surface(SequentialSurface {
        surface_type: SurfaceType::Plane,
        vertex_z: 60.0,
        aperture_radius: 30.0,
        ..Default::default()
    })
    .unwrap();
    s
}

fn surf_type(t: SurfaceType) -> f64 {
    match t {
        SurfaceType::Plane => 0.0,
        SurfaceType::Sphere => 1.0,
        SurfaceType::Conic => 2.0,
        SurfaceType::EvenAsphere => 3.0,
    }
}

#[test]
fn analysis_matches_upstream() {
    let sys = lens();
    let mut rays: Vec<Ray> = (-4..=4)
        .map(|i| Ray {
            origin: Vec3::new(0.3 * f64::from(i), f64::from(i), -10.0),
            direction: Vec3::new(0.0, 0.01, 1.0),
            ..Default::default()
        })
        .collect();
    rays.push(Ray {
        origin: Vec3::new(0.0, 20.0, -10.0),
        ..Default::default()
    });
    let traced = sys.trace_many(&rays).unwrap();
    let spot = spot_diagram_at_plane(&traced, 70.0).unwrap();
    let mut sp: Vec<f64> = spot.points.iter().flat_map(|p| [p.x, p.y]).collect();
    sp.extend([spot.centroid_x, spot.centroid_y, spot.rms_radius]);
    close("spot", &sp, 1e-10, 1e-12);
    let fan: Vec<f64> = tangential_ray_fan(&traced, 70.0, 0.1)
        .iter()
        .flat_map(|p| [p.pupil_coordinate, p.transverse_aberration])
        .collect();
    close("fan", &fan, 1e-10, 1e-12);
    close(
        "distortion",
        &[distortion_percent(10.3, 10.0).unwrap()],
        1e-14,
        0.0,
    );

    let pupil: Vec<C> = (0..16)
        .map(|i| {
            let x = (f64::from(i) - 7.5) / 8.0;
            if x.abs() < 0.6 {
                C::from_polar(1.0, 0.4 * x * x)
            } else {
                C::new(0.0, 0.0)
            }
        })
        .collect();
    let fr: Vec<f64> = fraunhofer_1d(&pupil, 1e-4, 550e-9, 0.5)
        .unwrap()
        .iter()
        .flat_map(|s| [s.coordinate, s.intensity])
        .collect();
    close("fraunhofer", &fr, 1e-12, 1e-14);
    let psf = normalized_psf_1d(&pupil).unwrap();
    close("psf", &psf, 1e-12, 1e-15);
    close("mtf", &mtf_from_psf_1d(&psf).unwrap(), 1e-11, 1e-14);

    let zz: Vec<f64> = [
        (0, 0),
        (1, 1),
        (1, -1),
        (2, 0),
        (3, -1),
        (4, 2),
        (5, 3),
        (6, 0),
        (3, 0),
    ]
    .iter()
    .map(|(n, m)| zernike(*n, *m, 0.73, 1.1).unwrap())
    .collect();
    close("zernike", &zz, 1e-13, 0.0);
    let ws: Vec<WavefrontSample> = (0..40)
        .map(|i| {
            let fi = f64::from(i);
            let (rho, th) = (((fi + 0.5) / 40.0).sqrt(), 2.399963 * fi);
            WavefrontSample {
                rho,
                theta: th,
                opd: 0.3 + 0.2 * rho * th.cos() - 0.1 * (2.0 * rho * rho - 1.0)
                    + 0.05 * rho * rho * (2.0 * th).sin()
                    + 0.01 * (5.0 * fi).sin(),
            }
        })
        .collect();
    let fz: Vec<f64> = fit_zernike(&ws, 3)
        .unwrap()
        .iter()
        .flat_map(|c| [f64::from(c.n), f64::from(c.m), c.value])
        .collect();
    close("fit_zernike", &fz, 1e-10, 1e-13);
    close(
        "wf_rms",
        &[
            wavefront_rms(&ws, true).unwrap(),
            wavefront_rms(&ws, false).unwrap(),
            longitudinal_chromatic_shift(&[50.1, 49.8, 50.4]).unwrap(),
        ],
        1e-13,
        0.0,
    );

    let neff = uniaxial_extraordinary_index(1.658, 1.486, 0.6).unwrap();
    close(
        "biref",
        &[
            neff,
            birefringent_retardance(1.658, neff, 2e-5, 633e-9).unwrap(),
        ],
        1e-13,
        0.0,
    );
    let j = apply_linear_retarder(
        JonesVector {
            x: C::new(0.6, 0.1),
            y: C::new(0.3, -0.7),
        },
        1.2,
        0.35,
    );
    close("retarder", &[cx(j.x), cx(j.y)].concat(), 1e-13, 1e-15);

    let mat = NonlinearOpticalMaterial {
        relative_permittivity: 2.25,
        chi3_m2_per_v2: C::new(2e-22, 1e-24),
        saturation_intensity_w_per_m2: 5e12,
    };
    let r1 =
        evaluate_nonlinear_polarization(C::new(3e7, 1e7), &mat, FREE_SPACE_WAVE_IMPEDANCE).unwrap();
    let r2 = evaluate_nonlinear_polarization(C::new(0.0, 0.0), &mat, 250.0).unwrap();
    let m3 = NonlinearOpticalMaterial {
        relative_permittivity: 1.5,
        chi3_m2_per_v2: C::new(1e-21, 0.0),
        saturation_intensity_w_per_m2: 0.0,
    };
    let r3 =
        evaluate_nonlinear_polarization(C::new(1e8, 0.0), &m3, FREE_SPACE_WAVE_IMPEDANCE).unwrap();
    let mut nl = Vec::new();
    for r in [r1, r2, r3] {
        nl.push(r.intensity_w_per_m2);
        for c in [
            r.linear_polarization_c_per_m2,
            r.nonlinear_polarization_c_per_m2,
            r.total_polarization_c_per_m2,
            r.effective_relative_permittivity,
        ] {
            nl.extend(cx(c));
        }
    }
    // Complex division: libgcc `__divdc3` scales differently from the textbook formula.
    close("nonlinear", &nl, 1e-12, 0.0);

    let tf = RigidTransform::from_euler_xyz(Vec3::new(1.0, -2.0, 0.5), 0.3, -0.2, 0.9);
    let tr = Ray {
        origin: Vec3::new(0.4, 0.5, -1.0),
        direction: Vec3::new(0.1, 0.2, 0.97),
        ..Default::default()
    };
    let (ta, ti) = (tf.apply(tr), tf.inverse(tf.apply(tr)));
    close(
        "transform",
        &[
            ta.origin.x,
            ta.origin.y,
            ta.origin.z,
            ta.direction.x,
            ta.direction.y,
            ta.direction.z,
            ti.origin.x,
            ti.origin.y,
            ti.origin.z,
        ],
        1e-14,
        1e-15,
    );

    let mut apv = Vec::new();
    for kind in [
        ApertureKind::Circular,
        ApertureKind::Rectangular,
        ApertureKind::Elliptical,
    ] {
        let ap = Aperture {
            kind,
            x_radius: 1.0,
            y_radius: 0.5,
            stop: false,
        };
        apv.extend(
            [
                ap.contains(0.8, 0.45),
                ap.contains(0.6, 0.3),
                ap.contains(0.95, 0.1),
            ]
            .map(|b| f64::from(u8::from(b))),
        );
    }
    let mut wm = WavelengthModel::default();
    for (nm, weight, primary) in [(486.1, 1.0, false), (587.6, 2.0, true), (656.3, 1.0, false)] {
        wm.add(Wavelength {
            nm,
            weight,
            primary,
        })
        .unwrap();
    }
    apv.push(wm.primary_nm().unwrap());
    close("models", &apv, 0.0, 0.0);

    let pt = |x, y, c| PolynomialTerm {
        x_power: x,
        y_power: y,
        coefficient: c,
    };
    let ff = PolynomialFreeformSurface::new(
        3.0,
        vec![
            pt(2, 0, 0.01),
            pt(0, 2, 0.015),
            pt(1, 2, -0.002),
            pt(3, 1, 1e-4),
        ],
        5.0,
    )
    .unwrap();
    let fr1 = Ray {
        origin: Vec3::new(0.5, -0.8, -2.0),
        direction: Vec3::new(0.05, 0.1, 1.0).normalized().unwrap(),
        ..Default::default()
    };
    let (p, n) = ff.intersect(&fr1, 40, 1e-12).unwrap();
    close(
        "freeform",
        &[ff.sag(1.2, -0.7), 1.0, p.x, p.y, p.z, n.x, n.y, n.z],
        1e-12,
        1e-15,
    );
}

#[test]
fn design_matches_upstream() {
    let sys = lens();
    for (fmt, name) in [
        (LensInterchangeFormat::ZemaxZmx, "zmx"),
        (LensInterchangeFormat::CodevSeq, "seq"),
        (LensInterchangeFormat::OsloLen, "len"),
    ] {
        assert_eq!(
            export_lens_system(&sys, fmt),
            text(&format!("export_{name}")),
            "{name}"
        );
        // Round trip through each dialect.
        let back = import_lens_system(&export_lens_system(&sys, fmt), fmt).unwrap();
        assert_eq!(back.surfaces().len(), 3);
    }
    let imp = import_lens_system(
        "! hdr\nSURF 0\nTYPE STANDARD\nCURV 0.02 # c\nDISZ 1.5\nDIAM 8\nGLAS N-BK7 1.5168\nCONI -1\nsurface 1\ntype evenasph\nRD -40\nTH 4\nAP 8\nGLASS 1\nCOEF 2 -3e-7\nASPH 9 4\nSURF 2\nTYPE STANDARD\nCURV 0\nDISZ 30\n",
        LensInterchangeFormat::OsloLen,
    )
    .unwrap();
    let iv: Vec<f64> = imp
        .surfaces()
        .iter()
        .flat_map(|s| {
            let mut v = vec![
                surf_type(s.surface_type),
                s.vertex_z,
                s.radius,
                s.aperture_radius,
                s.refractive_index_after,
                s.conic_constant,
            ];
            v.extend(s.even_coefficients);
            v
        })
        .collect();
    close("import", &iv, 1e-15, 0.0);
    let js = sequential_system_to_json(&sys, 1.33);
    assert_eq!(js, text("json"));
    let back =
        sequential_system_from_json(&format!("junk {{\"type\":\"sphere\",\"z\":1}} {js}")).unwrap();
    let mut bv = vec![back.object_space_index()];
    for s in back.surfaces() {
        bv.extend([
            surf_type(s.surface_type),
            s.vertex_z,
            s.radius,
            s.even_coefficients[3],
        ]);
    }
    close("json_back", &bv, 1e-15, 0.0);

    let mut scene = NonSequentialScene::default();
    scene
        .add_plane(NsPlane {
            point: Vec3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.1, 1.0),
            refractive_index_before: 1.0,
            refractive_index_after: 1.5,
            reflectance: 0.04,
            transmittance: 0.9,
            absorptance: 0.01,
            scatter_fraction: 0.05,
        })
        .unwrap();
    scene
        .add_plane(NsPlane {
            point: Vec3::new(0.0, 0.0, 10.0),
            normal: Vec3::new(0.0, 0.0, -1.0),
            refractive_index_before: 1.5,
            refractive_index_after: 1.0,
            reflectance: 0.1,
            transmittance: 0.85,
            absorptance: 0.05,
            scatter_fraction: 0.0,
        })
        .unwrap();
    let seed = NsRay {
        ray: Ray {
            origin: Vec3::new(0.1, 0.2, 0.0),
            direction: Vec3::new(0.02, -0.03, 1.0).normalized().unwrap(),
            ..Default::default()
        },
        ..Default::default()
    };
    let out = scene.trace(seed, 4, 1e-4).unwrap();
    let ns: Vec<f64> = out
        .iter()
        .flat_map(|r| {
            [
                r.power,
                r.bounces as f64,
                r.ray.origin.x,
                r.ray.origin.y,
                r.ray.origin.z,
                r.ray.direction.x,
                r.ray.direction.y,
                r.ray.direction.z,
            ]
        })
        .collect();
    close("ns_trace", &ns, 1e-12, 1e-14);
    let s = summarize_stray_light(&out, 1.0).unwrap();
    close(
        "ns_summary",
        &[
            s.total_power,
            s.reflected_power,
            s.transmitted_power,
            s.scattered_power,
            s.absorbed_power,
        ],
        1e-13,
        1e-16,
    );
    let gh: Vec<f64> = rank_ghost_paths(&out, 2)
        .iter()
        .flat_map(|g| [g.bounces as f64, g.power])
        .collect();
    close("ghosts", &gh, 1e-13, 0.0);

    let rosen = |x: &[f64]| {
        (1.0 - x[0]) * (1.0 - x[0])
            + 10.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0])
            + 0.5 * x[2] * x[2]
    };
    let cd =
        coordinate_descent(&rosen, vec![-1.2, 1.0, 0.7], vec![0.5, 0.5, 0.5], 300, 1e-8).unwrap();
    let mut v = cd.variables.clone();
    v.extend([
        cd.merit,
        cd.iterations as f64,
        f64::from(u8::from(cd.converged)),
    ]);
    close("cd", &v, 1e-12, 1e-15);
    let de = differential_evolution(
        &rosen,
        &[-2.0, -1.0, -1.0],
        &[2.0, 3.0, 1.0],
        12,
        40,
        0.7,
        0.85,
        11,
    )
    .unwrap();
    let mut v = de.variables.clone();
    v.extend([de.merit, de.iterations as f64]);
    close("de", &v, 1e-12, 1e-15);
    let p0 = [0.3, -0.4, 2.0];
    close(
        "sens",
        &finite_difference_sensitivity(&rosen, &p0, 1e-6).unwrap(),
        1e-9,
        0.0,
    );
    let mc = monte_carlo_tolerance(&rosen, &p0, &[0.01, 0.02, 0.0], 25, 3).unwrap();
    close(
        "mc",
        &[
            mc.mean,
            mc.standard_deviation,
            mc.minimum,
            mc.maximum,
            mc.samples as f64,
        ],
        1e-12,
        0.0,
    );
    let f = |x: Dual<1>| x * x * x / (x + Dual::constant(2.0)) - x;
    close("autodiff", &[autodiff_derivative(f, 1.3)], 1e-14, 0.0);
}
