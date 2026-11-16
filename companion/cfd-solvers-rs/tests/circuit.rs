//! Analytic checks of the circuit port: closed-form DC/AC/transient responses, KCL at a
//! nonlinear operating point, transformer and S-parameter identities, pole/zero recovery,
//! netlist files with `.INCLUDE` / `.LIB`, RAW round trips and OSDI symbol discovery.
use cfd_solvers::circuit::*;
use std::{f64::consts::PI, sync::Arc};

fn cfg() -> NewtonConfig {
    NewtonConfig::default()
}

#[test]
fn divider_and_controlled_sources_close_form() {
    let c = Circuit::parse_spice(
        "V1 in 0 12\nR1 in mid 3k\nR2 mid 0 1k\nE1 out 0 mid 0 -2\nRl out 0 1k\n",
    )
    .unwrap();
    let op = c.dc_operating_point(&cfg()).unwrap();
    assert!(op.converged);
    let mid = c.voltage(&op, "mid").unwrap();
    assert!((mid - 3.0).abs() < 1e-6); // gmin = 1e-12 S on every node
    assert!((c.voltage(&op, "out").unwrap() + 6.0).abs() < 1e-6);
    // Source current flows + → − through the source: −12 V / 4 kΩ.
    assert!((op.branch_current[0] + 3e-3).abs() < 1e-9);
}

#[test]
fn diode_operating_point_satisfies_kcl() {
    let c = Circuit::parse_spice("V1 in 0 5\nR1 in a 1k\nD1 a 0 DX\n.model DX D(IS=1e-14 N=1)\n")
        .unwrap();
    let op = c.dc_operating_point_default_homotopy(&cfg()).unwrap();
    let va = c.voltage(&op, "a").unwrap();
    let vt = 1.380649e-23 * 300.0 / 1.602176634e-19;
    let id = 1e-14 * ((va / vt).exp() - 1.0);
    let ir = (5.0 - va) / 1e3;
    assert!((id - ir).abs() < 1e-9 * ir, "diode {id} resistor {ir}");
    assert!(va > 0.6 && va < 0.75);
}

#[test]
fn rc_step_matches_exponential() {
    let c = Circuit::parse_spice("V1 in 0 PWL(0 1 1 1)\nR1 in out 1k\nC1 out 0 1u\n").unwrap();
    let pts = c
        .transient(1e-5, 300, TransientMethod::Bdf2, &cfg())
        .unwrap();
    let out = c.find_node("out").unwrap();
    for p in pts.iter().skip(20).step_by(40) {
        let exact = 1.0 - (-p.time / 1e-3).exp();
        assert!(
            (p.node_voltage[out] - exact).abs() < 5e-4,
            "t={} {} vs {exact}",
            p.time,
            p.node_voltage[out]
        );
    }
    let ad = c
        .transient_adaptive(
            3e-3,
            &AdaptiveTransientConfig {
                method: TransientMethod::Bdf2,
                initial_step: 1e-6,
                maximum_step: 1e-4,
                relative_tolerance: 1e-5,
                absolute_tolerance: 1e-8,
                ..Default::default()
            },
            &cfg(),
        )
        .unwrap();
    let last = ad.points.last().unwrap();
    assert!((last.time - 3e-3).abs() < 1e-15);
    assert!((last.node_voltage[out] - (1.0 - (-3.0f64).exp())).abs() < 1e-3);
    assert!(ad.accepted_steps > 10);
}

#[test]
fn series_rlc_resonance_and_pole_zero() {
    let c =
        Circuit::parse_spice("V1 in 0 0 AC 1\nR1 in a 10\nL1 a b 1m\nC1 b out 1u\nRo out 0 10\n")
            .unwrap();
    let f0 = 1.0 / (2.0 * PI * (1e-3f64 * 1e-6).sqrt());
    let at = |f: f64| {
        c.ac_voltage(&c.ac(f, None, &cfg()).unwrap(), "out")
            .unwrap()
            .norm()
    };
    assert!((at(f0) - 0.5).abs() < 1e-9, "resonance gain {}", at(f0));
    assert!(at(f0 * 0.5) < at(f0) && at(f0 * 2.0) < at(f0));

    let lp = Circuit::parse_spice("V1 in 0 0 AC 1\nR1 in out 1k\nC1 out 0 1u\n").unwrap();
    let r = pole_zero_analysis(
        &lp,
        "out",
        &PoleZeroConfig {
            start_hz: 1.0,
            stop_hz: 1e5,
            samples: 30,
            denominator_order: 1,
            numerator_order: 0,
        },
        &cfg(),
    )
    .unwrap();
    assert_eq!(r.poles_rad_per_s.len(), 1);
    assert!(
        (r.poles_rad_per_s[0].re + 1000.0).abs() < 1e-6 && r.poles_rad_per_s[0].im.abs() < 1e-6
    );
    let dp = estimate_dominant_pole(&lp, "out", 1.0, 1e5, 200, &cfg()).unwrap();
    assert!(dp.found && (dp.pole_frequency_hz / (1000.0 / (2.0 * PI)) - 1.0).abs() < 1e-2);
}

#[test]
fn transformer_and_matched_s_parameters() {
    let mut c = Circuit::new();
    let (p, s) = (c.node("p"), c.node("s"));
    c.add_voltage_source("V1", p, 0, 10.0, Complex::default(), None)
        .unwrap();
    c.add_ideal_transformer("T1", p, 0, s, 0, 5.0).unwrap();
    c.add_resistor("RL", s, 0, 4.0).unwrap();
    let op = c.dc_operating_point(&cfg()).unwrap();
    assert!((c.voltage(&op, "s").unwrap() - 2.0).abs() < 1e-9);
    // Power conservation: primary delivers V·I = secondary V²/R.
    assert!((10.0 * -op.branch_current[0] - 1.0).abs() < 1e-9);

    // 6.02 dB T attenuator matched to 50 Ω: S11 = 0, S21 = 1/2.
    let pad = Circuit::parse_spice(
        "R1 a m 16.666666666666667\nR2 m b 16.666666666666667\nR3 m 0 66.666666666666667\n",
    )
    .unwrap();
    let sp = two_port_s_parameters(&pad, &AcPort::new("a"), &AcPort::new("b"), 1e6, 50.0).unwrap();
    assert!(
        sp.a11.norm() < 1e-9 && (sp.a21.norm() - 0.5).abs() < 1e-9,
        "{sp:?}"
    );
    let ports = [CircuitPort::new(1, 0), CircuitPort::new(3, 0)];
    let s = extract_s_parameters(&pad, &ports, 1e6, &cfg()).unwrap();
    assert!(s[(0, 0)].norm() < 1e-9 && (s[(1, 0)].norm() - 0.5).abs() < 1e-9);
}

#[test]
fn netlist_files_with_include_and_library_sections() {
    let dir = std::env::temp_dir().join(format!("cfd_circuit_files_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("models.lib"),
        ".lib fast\n.model DF D(IS=1e-12)\n.endl\n.lib slow\n.model DF D(IS=1e-16)\n.endl\n",
    )
    .unwrap();
    std::fs::write(dir.join("load.inc"), "R2 out 0 2k\n").unwrap();
    std::fs::write(
        dir.join("top.cir"),
        "* top\nV1 in 0 5\nR1 in out 1k\n.include \"load.inc\"\nD1 out 0 DF\n.lib models.lib slow\n",
    )
    .unwrap();
    let c = Circuit::parse_spice_file(dir.join("top.cir"), "").unwrap();
    let op = c.dc_operating_point_default_homotopy(&cfg()).unwrap();
    let v = c.voltage(&op, "out").unwrap();
    assert!(v > 0.6 && v < 1.0, "{v}");
    std::fs::write(dir.join("loop.cir"), ".include loop.cir\n").unwrap();
    assert!(Circuit::parse_spice_file(dir.join("loop.cir"), "").is_err());
    assert!(Circuit::parse_spice_file(dir.join("models.lib"), "missing").is_err());
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn raw_round_trip_preserves_values() {
    let c = Circuit::parse_spice("V1 in 0 SIN(0 1 1k)\nR1 in out 1k\nC1 out 0 1u\n").unwrap();
    let t = c
        .transient(1e-4, 10, TransientMethod::BackwardEuler, &cfg())
        .unwrap();
    let ds = raw_from_transient(&t, c.node_names(), "rt").unwrap();
    let back = read_spice_raw(&spice_raw_text(&ds).unwrap()).unwrap();
    assert_eq!(back, ds);
    let mut buf = Vec::new();
    write_spice_raw(&mut buf, &ds).unwrap();
    assert!(String::from_utf8(buf)
        .unwrap()
        .starts_with("Title: rt\nPlotname: Transient Analysis\n"));
}

#[test]
fn osdi_adapter_drives_a_static_device() {
    let adapter = OsdiDescriptorAdapter::new(
        2,
        Arc::new(|v: &[f64], i: &mut [f64], j: &mut [f64]| {
            let g = 2e-3;
            i[0] = g * (v[0] - v[1]);
            i[1] = -i[0];
            j.copy_from_slice(&[g, -g, -g, g]);
            Ok(())
        }),
    )
    .unwrap();
    let mut c = Circuit::new();
    let a = c.node("a");
    c.add_voltage_source("V1", a, 0, 1.0, Complex::default(), None)
        .unwrap();
    c.add_static_device("X", vec![a, 0], adapter.static_evaluator())
        .unwrap();
    let op = c.dc_operating_point(&cfg()).unwrap();
    assert!((op.branch_current[0] + 2e-3).abs() < 1e-12);
    assert!(adapter.stamp(&[0.0], &mut [0.0], &mut [0.0]).is_err());
}

/// Minimal ELF64 LE image with .dynsym/.dynstr/.data sections exporting the OSDI symbols.
fn synthetic_osdi(major: u32, minor: u32, count: u32, with_descriptors: bool) -> Vec<u8> {
    let mut names = vec![
        "OSDI_VERSION_MAJOR",
        "OSDI_VERSION_MINOR",
        "OSDI_NUM_DESCRIPTORS",
    ];
    if with_descriptors {
        names.push("OSDI_DESCRIPTORS");
    }
    let mut dynstr = vec![0u8];
    let mut name_offsets = Vec::new();
    for n in &names {
        name_offsets.push(dynstr.len() as u32);
        dynstr.extend_from_slice(n.as_bytes());
        dynstr.push(0);
    }
    let data: Vec<u8> = [major, minor, count, 0]
        .iter()
        .flat_map(|v| v.to_le_bytes())
        .collect();
    let data_addr = 0x4000u64;
    let mut dynsym = vec![0u8; 24];
    for (k, off) in name_offsets.iter().enumerate() {
        let mut e = vec![0u8; 24];
        e[0..4].copy_from_slice(&off.to_le_bytes());
        e[4] = 0x11;
        e[6..8].copy_from_slice(&3u16.to_le_bytes());
        e[8..16].copy_from_slice(&(data_addr + 4 * k as u64).to_le_bytes());
        e[16..24].copy_from_slice(&4u64.to_le_bytes());
        dynsym.extend(e);
    }
    let mut file = vec![0u8; 64];
    let dynsym_off = file.len() as u64;
    file.extend(&dynsym);
    let dynstr_off = file.len() as u64;
    file.extend(&dynstr);
    let data_off = file.len() as u64;
    file.extend(&data);
    let shoff = file.len() as u64;
    let section = |kind: u32, addr: u64, off: u64, size: u64, link: u32, entsize: u64| {
        let mut h = vec![0u8; 64];
        h[4..8].copy_from_slice(&kind.to_le_bytes());
        h[16..24].copy_from_slice(&addr.to_le_bytes());
        h[24..32].copy_from_slice(&off.to_le_bytes());
        h[32..40].copy_from_slice(&size.to_le_bytes());
        h[40..44].copy_from_slice(&link.to_le_bytes());
        h[56..64].copy_from_slice(&entsize.to_le_bytes());
        h
    };
    file.extend(section(0, 0, 0, 0, 0, 0));
    file.extend(section(11, 0, dynsym_off, dynsym.len() as u64, 2, 24));
    file.extend(section(3, 0, dynstr_off, dynstr.len() as u64, 0, 0));
    file.extend(section(1, data_addr, data_off, data.len() as u64, 0, 0));
    file[0..4].copy_from_slice(b"\x7fELF");
    file[4] = 2;
    file[5] = 1;
    file[0x28..0x30].copy_from_slice(&shoff.to_le_bytes());
    file[0x3a..0x3c].copy_from_slice(&64u16.to_le_bytes());
    file[0x3c..0x3e].copy_from_slice(&4u16.to_le_bytes());
    file
}

#[test]
fn osdi_library_reads_exported_abi_symbols() {
    let lib = OsdiLibrary::from_bytes("bsim.osdi", &synthetic_osdi(0, 4, 2, true)).unwrap();
    assert_eq!(
        (
            lib.version_major(),
            lib.version_minor(),
            lib.descriptor_count()
        ),
        (0, 4, 2)
    );
    assert!(lib.descriptors_file_offset().is_some());
    assert!(OsdiLibrary::from_bytes("x", &synthetic_osdi(0, 4, 2, false)).is_err());
    assert!(OsdiLibrary::from_bytes("x", &synthetic_osdi(0, 3, 0, false)).is_ok());
    assert!(OsdiLibrary::from_bytes("x", b"not an elf").is_err());
    assert!(OsdiLibrary::open("/nonexistent/module.osdi").is_err());
}
