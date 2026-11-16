//! Analytic checks for DG advection, Vlasov solvers, molecular dynamics and beamforming.
use cfd_solvers::{acoustics::*, atomistics::*, dg::*, kinetic::*};
use std::f64::consts::PI;

#[test]
fn dg_advects_a_sine_one_period_with_spectral_accuracy_and_conserves_mass() {
    let mut dg = DgAdvection1D::new(DgAdvection1DConfig {
        elements: 8,
        order: 5,
        length: 1.0,
        speed: 1.0,
    })
    .unwrap();
    dg.initialize(|x| (2.0 * PI * x).sin());
    let m0 = dg.mass();
    for _ in 0..1000 {
        dg.step(1e-3).unwrap();
    }
    let p = dg.reference_nodes().to_vec();
    let h = 1.0 / 8.0;
    let mut worst = 0.0f64;
    for e in 0..8 {
        for (k, xi) in p.iter().enumerate() {
            let x = h * (e as f64 + 0.5 * (xi + 1.0));
            worst = worst.max((dg.values()[e * p.len() + k] - (2.0 * PI * x).sin()).abs());
        }
    }
    assert!(worst < 1e-5, "error {worst}");
    assert!((dg.mass() - m0).abs() < 1e-13);
}

#[test]
fn lobatto_weights_integrate_polynomials_and_chebyshev_derivative_is_exact() {
    let basis = HighOrderNodalBasis1D::new(6).unwrap();
    let d = SumFactorizedDerivative1D::new(&basis);
    let u: Vec<f64> = basis.nodes().iter().map(|x| x.powi(5) - 2.0 * x).collect();
    let mut du = vec![0.0; 7];
    d.apply(&u, &mut du, 1.0).unwrap();
    for (x, g) in basis.nodes().iter().zip(&du) {
        assert!((g - (5.0 * x.powi(4) - 2.0)).abs() < 1e-11);
    }
    let nodes = gauss_lobatto_nodes(4);
    assert!(
        nodes.windows(2).all(|w| w[1] > w[0]) && (nodes[1] + (3.0f64 / 7.0).sqrt()).abs() < 1e-14
    );
}

#[test]
fn vlasov_conserves_mass_and_poisson_field_matches_cosine_perturbation() {
    let mut vl = Vlasov1D1V::new(Vlasov1D1VConfig::default()).unwrap();
    vl.initialize(|x, v| (-v * v).exp() * (1.0 + 0.1 * (2.0 * PI * x).cos()))
        .unwrap();
    let m0 = vl.mass();
    for _ in 0..50 {
        vl.step(1e-3, 0.0).unwrap();
    }
    assert!((vl.mass() - m0).abs() < 1e-10 * m0);
    let n = 64;
    let mut vp = VlasovPoisson1D1V::new(VlasovPoisson1D1VConfig {
        spatial_cells: n,
        velocity_cells: 64,
        spatial_length: 2.0 * PI,
        maximum_velocity: 6.0,
        neutralizing_density: 1.0,
    })
    .unwrap();
    vp.initialize(|x, v| (-0.5 * v * v).exp() / (2.0 * PI).sqrt() * (1.0 + 0.05 * x.cos()))
        .unwrap();
    // dE/dx = rho with rho = 0.05 cos x gives E = 0.05 sin x (density discretization error ~1e-9).
    for (i, e) in vp.electric_field().iter().enumerate() {
        let x = (i as f64 + 0.5) * 2.0 * PI / n as f64;
        assert!((e - 0.05 * x.sin()).abs() < 1e-6, "E[{i}] = {e}");
    }
}

#[test]
fn verlet_chain_conserves_energy_and_lennard_jones_minimum_is_force_free() {
    let mut chain = HarmonicChain1D::new(HarmonicChainConfig::default()).unwrap();
    chain.set_velocity(3, 0.4).unwrap();
    let e0 = chain.energy().total();
    for _ in 0..2000 {
        chain.step(1e-3).unwrap();
    }
    assert!((chain.energy().total() - e0).abs() < 1e-4 * e0);
    let r = 2f64.powf(1.0 / 6.0);
    let pos = [1.0, 1.0 + r];
    let mut nl = PeriodicNeighborList1D::new(NeighborList1DConfig {
        periodic_length: 10.0,
        cutoff: 2.5,
        skin: 0.2,
    })
    .unwrap();
    nl.rebuild(&pos).unwrap();
    let mut f = [0.0; 2];
    let d = LennardJones1DForceField::new(LennardJonesParameters::default())
        .unwrap()
        .compute(&pos, &nl, &mut f)
        .unwrap();
    assert!(f.iter().all(|v| v.abs() < 1e-12) && (d.potential_energy + 1.0).abs() < 1e-12);
}

#[test]
fn beamforming_aligns_delayed_pulses() {
    let pulse = |shift: usize| {
        (0..32)
            .map(|i| if i == 5 + shift { 1.0 } else { 0.0 })
            .collect::<Vec<f64>>()
    };
    let (a, b) = (pulse(0), pulse(2));
    // Sensor b is 2 samples (0.3 m at c = 1500, dt = 1e-4) farther from the focus.
    let s = [
        BeamformSample {
            x_m: 0.0,
            y_m: 0.0,
            pressure: &a,
        },
        BeamformSample {
            x_m: 0.3,
            y_m: 0.0,
            pressure: &b,
        },
    ];
    let out = delay_and_sum_beamform(&s, -1e-9, 0.0, 1e-4, 1500.0, 16).unwrap();
    assert!((out[5] - 1.0).abs() < 1e-6);
}
