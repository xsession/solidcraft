//! Value-by-value comparison with the C++ upstream FDTD units. Goldens come from
//! `parity/oracle_fdtd.cpp`.
use cfd_solvers::{core::Vec3, fdtd::*, physics::*, Complex64 as C};
use serde_json::Value;
use std::f64::consts::PI;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-fdtd.json")).unwrap()
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

fn cx(c: C) -> [f64; 2] {
    [c.re, c.im]
}

#[test]
fn maxwell3d_matches_upstream() {
    let c = Maxwell3DConfig {
        nx: 8,
        ny: 7,
        nz: 6,
        dx: 1e-3,
        dy: 1.2e-3,
        dz: 0.9e-3,
        courant: 0.9,
        epsilon_r: 2.0,
        epsilon_ry: 3.0,
        ..Default::default()
    };
    let mut m = Maxwell3D::new(c).unwrap();
    m.initialize_gaussian_ez(1.5, 0.3).unwrap();
    m.add_soft_ez_source(3, 3, 2, 0.25).unwrap();
    m.step(12);
    let f: Vec<f64> = [m.ex(), m.ey(), m.ez(), m.hx(), m.hy(), m.hz()].concat();
    close("m3_fields", &f, 1e-12, 0.0);
    close(
        "m3_scalars",
        &[m.dt(), m.energy(), m.max_field()],
        1e-12,
        0.0,
    );
    let cp = Maxwell3DConfig {
        boundary: Boundary3D::Pmc,
        epsilon_ry: 0.0,
        mu_r: 1.3,
        ..c
    };
    let mut mp = Maxwell3D::new(cp).unwrap();
    mp.initialize_gaussian_ez(1.0, 0.2).unwrap();
    mp.step(9);
    let (nx, ny) = (mp.nx(), mp.ny());
    close(
        "m3_pmc",
        &[
            mp.dt(),
            mp.energy(),
            mp.max_field(),
            mp.ez()[nx * ny * 2 + nx * 3 + 4],
            mp.hx()[nx * ny * 3 + nx * 2 + 1],
        ],
        1e-12,
        0.0,
    );
    let cv = Maxwell3DConfig {
        nx: 4,
        ny: 4,
        nz: 5,
        dx: 2e-3,
        dy: 1e-3,
        dz: 1.5e-3,
        ..Default::default()
    };
    let mut mv = Maxwell3D::new(cv).unwrap();
    mv.initialize_gaussian_ez(2.0, 0.5).unwrap();
    mv.step(3);
    assert_eq!(maxwell3d_vtk_ascii(&mv), goldens()["vtk"].as_str().unwrap());
    close(
        "hdf5",
        &[f64::from(u8::from(hdf5_xdmf_output_available()))],
        0.0,
        0.0,
    );
    let mat = Material {
        name: "x".into(),
        density: 1000.0,
        heat_capacity: 1000.0,
        thermal_conductivity: 1.0,
        relative_permittivity: [2.0, 3.0, 4.0],
        ..Default::default()
    };
    let b = apply_to_fdtd(Maxwell3DConfig::default(), &mat).unwrap();
    close(
        "bridge",
        &[b.epsilon_rx, b.epsilon_ry, b.epsilon_rz, b.epsilon_r],
        0.0,
        0.0,
    );
}

#[test]
fn cylindrical_matches_upstream() {
    let cc = CylindricalTMConfig {
        nr: 10,
        nz: 12,
        dr: 1e-3,
        dz: 1.5e-3,
        courant: 0.4,
        epsilon_r: 2.2,
        mu_r: 1.0,
    };
    let mut cyl = CylindricalTM::new(cc).unwrap();
    cyl.initialize_gaussian_ez(0.0, 8e-3, 2e-3, 1.2).unwrap();
    cyl.add_soft_ez_source(2, 5, 0.1).unwrap();
    cyl.step(25);
    close("cyl_ez", cyl.ez(), 1e-11, 0.0);
    close(
        "cyl_scalars",
        &[cyl.dt(), cyl.energy(), cyl.max_field()],
        1e-11,
        0.0,
    );
    let lv: Vec<f64> = make_cylindrical_multigrid_levels(40, 23, 0.02, 0.05, 6)
        .unwrap()
        .iter()
        .flat_map(|l| {
            [
                l.radial_cells as f64,
                l.axial_cells as f64,
                l.radial_spacing,
                l.axial_spacing,
            ]
        })
        .collect();
    close("levels", &lv, 1e-15, 0.0);
}

#[test]
fn monitors_and_ports_match_upstream() {
    let mut dm = DftMonitor::new(2e9).unwrap();
    let mut tp = TimeProbe::default();
    for i in 0..=200 {
        let t = f64::from(i) * 2e-12;
        let v = 0.7 * (2.0 * PI * 2e9 * t + 0.3).cos() + 0.1 * (2.0 * PI * 5e9 * t).sin();
        dm.sample(t, v).unwrap();
        tp.sample(t, v).unwrap();
    }
    let mut d = cx(dm.integral()).to_vec();
    d.extend([
        dm.amplitude(400e-12).unwrap(),
        dm.samples() as f64,
        tp.time().len() as f64,
    ]);
    close("dft", &d, 1e-11, 0.0);
    let (mut p1, mut p2, mut p3) = (
        WavePort1D::new(1e9, 50.0, 1).unwrap(),
        WavePort1D::new(1e9, 50.0, 1).unwrap(),
        WavePort1D::new(1e9, 75.0, -1).unwrap(),
    );
    for i in 0..=300 {
        let (t, w) = (f64::from(i) * 5e-12, 2.0 * PI * 1e9);
        let (a, b, tr) = (
            (w * t).sin(),
            0.3 * (w * t + 1.0).sin(),
            0.8 * (w * t - 0.7).sin(),
        );
        p1.sample(t, a + b, (-a + b) / 50.0).unwrap();
        p2.sample(t, tr, -tr / 50.0).unwrap();
        p3.sample(t, 0.4 * a, 0.2 * b).unwrap();
    }
    let sp = s_parameters(&p1, &p2).unwrap();
    let inst = p3.instantaneous(1.0, 0.01);
    let v = [
        cx(sp.s11),
        cx(sp.s21),
        cx(p3.forward_spectrum()),
        cx(p3.backward_spectrum()),
        [inst.forward, inst.backward],
    ]
    .concat();
    close("ports", &v, 1e-10, 0.0);
}

#[test]
fn postprocessing_and_emc_match_upstream() {
    let surf: Vec<NearFieldSurfaceSample> = (0..6)
        .map(|i| {
            let fi = f64::from(i);
            NearFieldSurfaceSample {
                position: Vec3::new(0.01 * fi.cos(), 0.01 * fi.sin(), 0.002 * fi),
                normal: Vec3::new(fi.cos(), fi.sin(), 0.1),
                area: 1e-5 * (1.0 + fi),
                electric: ComplexVec3 {
                    x: C::new(0.1 * fi, 0.2),
                    y: C::new(0.3, -0.1 * fi),
                    z: C::new(0.05, 0.0),
                },
                magnetic: ComplexVec3 {
                    x: C::new(0.001, -0.002 * fi),
                    y: C::new(0.0005 * fi, 0.001),
                    z: C::new(-0.001, 0.0003),
                },
            }
        })
        .collect();
    let ff = nf2ff(&surf, Vec3::new(0.3, 0.4, 0.8), 3e9, NF2FF_WAVE_IMPEDANCE).unwrap();
    close(
        "nf2ff",
        &[
            cx(ff.electric.x),
            cx(ff.electric.y),
            cx(ff.electric.z),
            [ff.radiation_intensity, 0.0],
        ]
        .concat()[..7],
        1e-12,
        0.0,
    );
    let mut g = SarGrid3D {
        nx: 5,
        ny: 4,
        nz: 3,
        dx: 5e-3,
        dy: 6e-3,
        dz: 4e-3,
        ..Default::default()
    };
    for n in 0..60 {
        let f = f64::from(n);
        g.ex.push(10.0 * f.sin());
        g.ey.push(5.0 * (0.7 * f).cos());
        g.ez.push(1.0 + 0.1 * f);
        g.conductivity.push(0.5 + 0.01 * f64::from(n % 7));
        g.density.push(1000.0 + 20.0 * f64::from(n % 5));
    }
    close("sar_local", &local_sar(&g).unwrap(), 1e-13, 0.0);
    close("sar_1g", &sar_1g(&g).unwrap(), 1e-12, 0.0);
    close("sar_10g", &sar_10g(&g).unwrap(), 1e-12, 0.0);

    let mc = MulticonductorTlineConfig {
        conductors: 3,
        inductance_h: vec![1e-7, 2e-7, 1.5e-7],
        capacitance_f: vec![1e-10, 2e-10, 0.5e-10],
        dt_s: 1e-10,
    };
    let mut st = MulticonductorTlineState {
        voltage: vec![1.0, 0.0, -0.5],
        current: vec![0.0, 0.01, 0.0],
    };
    for _ in 0..7 {
        advance_multiconductor_tline(&mc, &mut st, &[0.2, -0.1, 0.0]).unwrap();
    }
    close(
        "mtl",
        &[st.voltage.clone(), st.current.clone()].concat(),
        1e-13,
        0.0,
    );
    close(
        "graded",
        &[
            make_graded_axis(1.0, 6, 0.05, 1.3).unwrap(),
            make_graded_axis(2.0, 4, 0.1, 1.0).unwrap(),
        ]
        .concat(),
        1e-14,
        0.0,
    );
    let mut dw = Vec::new();
    for w in deposit_thin_wire_subcell(2.3e-3, 0.7e-3, 9.9e-3, 4, 3, 5, 1e-3, 1e-3, 2e-3)
        .unwrap()
        .iter()
        .chain(&deposit_thin_wire_subcell(-1.0, 1e-3, 0.0, 4, 3, 5, 1e-3, 1e-3, 2e-3).unwrap())
    {
        dw.extend([w.i as f64, w.j as f64, w.k as f64, w.weight]);
    }
    close("deposit", &dw, 1e-12, 0.0);
    let hs =
        huygens_equivalent_source([0.0, 2.0, 1.0], [1.0, -2.0, 0.5], [0.01, 0.02, -0.03]).unwrap();
    close(
        "huygens",
        &[hs.electric_surface_current, hs.magnetic_surface_current].concat(),
        1e-14,
        0.0,
    );
    let it =
        impedance_sheet_transmission(C::new(10.0, -40.0), C::new(377.0, 0.0), C::new(250.0, 5.0))
            .unwrap();
    let ts = [
        EmcTransferSample {
            frequency_hz: 1e6,
            gain: C::new(1.0, 0.0),
        },
        EmcTransferSample {
            frequency_hz: 1e7,
            gain: C::new(0.5, -0.2),
        },
        EmcTransferSample {
            frequency_hz: 1e8,
            gain: C::new(0.1, 0.3),
        },
    ];
    let e = [
        cx(it),
        cx(interpolate_emc_transfer(&ts, 5e7).unwrap()),
        cx(interpolate_emc_transfer(&ts, 1e5).unwrap()),
        cx(interpolate_emc_transfer(&ts, 1e9).unwrap()),
    ]
    .concat();
    close("emc", &e, 1e-14, 0.0);
}
