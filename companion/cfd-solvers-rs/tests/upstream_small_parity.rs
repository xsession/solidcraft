//! Value-by-value comparison with the C++ upstream DG, kinetic, atomistics and acoustics
//! post-processing units. Goldens come from `parity/oracle_small.cpp`.
use cfd_solvers::{acoustics::*, atomistics::*, dg::*, kinetic::*};
use serde_json::Value;
use std::f64::consts::PI;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-small.json")).unwrap()
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

const TIGHT: f64 = 1e-12;

#[test]
fn dg_matches_upstream() {
    for p in [1, 2, 4, 6] {
        close(&format!("lgl{p}"), &gauss_lobatto_nodes(p), TIGHT, 1e-15);
    }
    let mut dg = DgAdvection1D::new(DgAdvection1DConfig {
        elements: 6,
        order: 4,
        length: 2.0,
        speed: -0.7,
    })
    .unwrap();
    dg.initialize(|x| (PI * x).sin() + if x > 1.0 { 0.3 } else { 0.0 });
    let m0 = dg.mass();
    for _ in 0..25 {
        dg.step(0.01).unwrap();
    }
    close("dg", dg.values(), 1e-11, 1e-14);
    close("dg_mass", &[m0, dg.mass()], 1e-12, 1e-15);
    let basis = HighOrderNodalBasis1D::new(5).unwrap();
    close("cheb_nodes", basis.nodes(), TIGHT, 1e-15);
    close("cheb_d", basis.derivative_matrix(), TIGHT, 1e-13);
    close(
        "cheb_element",
        &basis.element_nodes(0.5, 2.0).unwrap(),
        TIGHT,
        0.0,
    );
    let mesh = [
        OneDimensionalInterval {
            left: 0.0,
            right: 1.0,
            level: 0,
        },
        OneDimensionalInterval {
            left: 1.0,
            right: 3.0,
            level: 1,
        },
        OneDimensionalInterval {
            left: 3.0,
            right: 4.0,
            level: 0,
        },
    ];
    let amr: Vec<f64> = basis
        .refine_marked(&mesh, &[2, 0])
        .unwrap()
        .iter()
        .flat_map(|c| [c.left, c.right, f64::from(c.level)])
        .collect();
    close("amr", &amr, 0.0, 0.0);
    let sf = SumFactorizedDerivative1D::new(&basis);
    let input: Vec<f64> = basis.nodes().iter().map(|x| x.exp()).collect();
    let mut out = vec![0.0; 6];
    sf.apply(&input, &mut out, 2.5).unwrap();
    close("sumfact", &out, 1e-12, 0.0);
    assert!(basis.refine_marked(&mesh, &[3]).is_err());
}

#[test]
fn kinetic_matches_upstream() {
    let mut vl = Vlasov1D1V::new(Vlasov1D1VConfig {
        spatial_cells: 8,
        velocity_cells: 10,
        spatial_length: 2.0,
        maximum_velocity: 3.0,
    })
    .unwrap();
    vl.initialize(|x, v| (-v * v).exp() * (1.0 + 0.2 * (PI * x).cos()) - 0.05)
        .unwrap();
    for i in 0..6 {
        vl.step(0.02, if i % 2 == 1 { -0.4 } else { 0.6 }).unwrap();
    }
    close("vlasov", vl.distribution(), TIGHT, 1e-16);
    close("vlasov_mass", &[vl.mass()], TIGHT, 0.0);
    let mut vp = VlasovPoisson1D1V::new(VlasovPoisson1D1VConfig {
        spatial_cells: 8,
        velocity_cells: 12,
        spatial_length: 4.0 * PI,
        maximum_velocity: 4.0,
        neutralizing_density: 0.0,
    })
    .unwrap();
    vp.initialize(|x, v| (-0.5 * v * v).exp() / (2.0 * PI).sqrt() * (1.0 + 0.1 * (0.5 * x).cos()))
        .unwrap();
    for _ in 0..5 {
        vp.step(0.05).unwrap();
    }
    close("vp_f", vp.distribution(), 1e-11, 1e-16);
    close("vp_e", vp.electric_field(), 1e-10, 1e-14);
    close("vp_rho", vp.charge_density(), TIGHT, 0.0);
    close("vp_scalars", &[vp.mass(), vp.field_energy()], 1e-10, 0.0);
}

#[test]
fn atomistics_matches_upstream() {
    let mut chain = HarmonicChain1D::new(HarmonicChainConfig {
        particles: 7,
        length: 3.5,
        mass: 2.0,
        spring_constant: 40.0,
    })
    .unwrap();
    chain.set_velocity(0, 0.3).unwrap();
    chain.set_velocity(4, -0.2).unwrap();
    for _ in 0..20 {
        chain.step(0.01).unwrap();
    }
    let e = chain.energy();
    close("chain_x", chain.positions(), TIGHT, 1e-15);
    close("chain_v", chain.velocities(), 1e-11, 1e-15);
    close("chain_e", &[e.kinetic, e.potential, e.total()], 1e-11, 0.0);
    let mut nl = PeriodicNeighborList1D::new(NeighborList1DConfig {
        periodic_length: 10.0,
        cutoff: 1.6,
        skin: 0.3,
    })
    .unwrap();
    let pos = [0.1, 1.2, 2.0, 3.4, 9.5, 8.6, 5.0];
    nl.rebuild(&pos).unwrap();
    let pairs: Vec<f64> = nl
        .pairs()
        .iter()
        .flat_map(|&(a, b)| [a as f64, b as f64])
        .collect();
    close("pairs", &pairs, 0.0, 0.0);
    let lj = LennardJones1DForceField::new(LennardJonesParameters {
        epsilon: 0.8,
        sigma: 1.0,
    })
    .unwrap();
    let mut f = vec![0.0; pos.len()];
    let d = lj.compute(&pos, &nl, &mut f).unwrap();
    close("lj_f", &f, TIGHT, 0.0);
    close(
        "lj_diag",
        &[
            d.potential_energy,
            d.virial,
            nl.minimum_image(7.3),
            nl.minimum_image(-12.0),
        ],
        TIGHT,
        1e-15,
    );
    let mut vel = vec![0.5, -0.2, 0.1, 0.7];
    velocity_rescaling_thermostat(&mut vel, 2.0, 0.3, 0.01, 0.1).unwrap();
    let mut length = 10.0;
    let mut bpos = [1.0, 2.0, 7.0];
    let factor = isotropic_barostat_scale(&mut bpos, &mut length, 2.0, 0.05, 0.5).unwrap();
    vel.extend([factor, length]);
    vel.extend(bpos);
    close("thermo_baro", &vel, TIGHT, 0.0);
}

#[test]
fn acoustics_postprocess_matches_upstream() {
    let (p1, p2, p3) = (
        [0.0, 1.0, 3.0, 2.0, 0.0, -1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 3.0, 2.0, 0.0, -1.0, 0.0],
        [1.0, 2.0],
    );
    let sensors = [
        BeamformSample {
            x_m: 0.0,
            y_m: 0.0,
            pressure: &p1,
        },
        BeamformSample {
            x_m: 0.3,
            y_m: 0.1,
            pressure: &p2,
        },
        BeamformSample {
            x_m: 0.1,
            y_m: 0.4,
            pressure: &p3,
        },
        BeamformSample {
            x_m: 1.0,
            y_m: 1.0,
            pressure: &[],
        },
    ];
    let mut v = delay_and_sum_beamform(&sensors, 0.2, 0.2, 1e-4, 1500.0, 6).unwrap();
    v.extend([
        plane_wave_intensity_w_m2(1e3, 1000.0, 1500.0).unwrap(),
        radiation_pressure_pa(3e4, 1500.0, 0.4).unwrap(),
        acoustic_heating_w_m3(3e4, 2.5).unwrap(),
    ]);
    close("acoustics", &v, TIGHT, 1e-15);
}
