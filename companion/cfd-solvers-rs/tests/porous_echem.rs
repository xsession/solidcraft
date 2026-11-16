//! Analytic and conservation checks for porous media and electrochemistry.
use cfd_solvers::{
    chemistry::*, electrochemistry::*, fvm::make_sheared_cartesian_hexa_mesh, porous::*,
};

#[test]
fn fractured_network_and_cahn_hilliard_conserve_their_inventories() {
    let mut fr = FracturedDarcyThmcNetwork::new(FracturedDarcyThmcConfig {
        cells: 4,
        reaction_rate: 0.0,
        ..Default::default()
    })
    .unwrap();
    for (a, b) in [(0, 1), (1, 2), (2, 3)] {
        fr.add_edge(FracturedDarcyEdge {
            first: a,
            second: b,
            transmissibility: 1.0,
            aperture: 3.0,
        })
        .unwrap();
    }
    fr.initialize(1.0, 0.5, 300.0, 0.2).unwrap();
    fr.set_cell(0, 2.0, 0.5, 320.0, 0.6).unwrap();
    let d0 = fr.diagnostics();
    for _ in 0..50 {
        fr.step(0.01).unwrap();
    }
    let d1 = fr.diagnostics();
    for (a, b) in [
        (d0.pressure_mass, d1.pressure_mass),
        (d0.temperature_mass, d1.temperature_mass),
        (d0.species_mass, d1.species_mass),
    ] {
        assert!((a - b).abs() < 1e-12 * a.abs().max(1.0));
    }
    let mut ch = CahnHilliard1D::new(CahnHilliard1DConfig::default(), None).unwrap();
    ch.initialize(0.4, 0.05).unwrap();
    let m0 = ch.mean();
    let f0 = ch.free_energy().unwrap();
    for _ in 0..200 {
        ch.step(2e-7).unwrap();
    }
    assert!((ch.mean() - m0).abs() < 1e-14 && ch.free_energy().unwrap() <= f0);
}

#[test]
fn acetic_acid_ph_matches_the_textbook_value() {
    // 0.1 M acetic acid, Ka = 1.8e-5: pH ~ 2.88.
    let r = equilibrate_acids(
        &[AcidFamily {
            total_mol_per_litre: 0.1,
            fully_protonated_charge: 0,
            dissociation_constants: vec![1.8e-5],
        }],
        0.0,
        1e-14,
    )
    .unwrap();
    assert!((r.ph - 2.88).abs() < 0.01, "{}", r.ph);
    assert!(r.charge_residual_mol_per_litre.abs() < 1e-14);
}

#[test]
fn butler_volmer_is_odd_symmetric_and_mixed_potential_balances_currents() {
    let up = butler_volmer_current_density(1.0, 0.03, 298.15, 1.0, 0.5, 0.5).unwrap();
    let down = butler_volmer_current_density(1.0, -0.03, 298.15, 1.0, 0.5, 0.5).unwrap();
    assert!((up + down).abs() < 1e-14);
    let r = |e: f64, i0: f64| ElectrodeReaction {
        name: "r".into(),
        equilibrium_potential: e,
        exchange_current_density: i0,
        ..Default::default()
    };
    let mp =
        solve_mixed_potential(&[r(-0.4, 1e-3), r(0.3, 1e-3)], 298.15, 0.0, 200, 1e-12).unwrap();
    // Equal exchange currents and transfer coefficients put the mixed potential at the midpoint.
    assert!((mp.potential + 0.05).abs() < 1e-9 && mp.total_current.abs() < 1e-12);
    let cell = solve_corrosion_cell_1d(&CorrosionCell1DConfig {
        metal_potential: 0.1,
        ..Default::default()
    })
    .unwrap();
    assert!(
        cell.converged
            && cell.current_density > 0.0
            && (cell.overpotential + cell.current_density * 1e-3 / 5.0 - 0.1).abs() < 1e-12
    );
}

#[test]
fn randles_low_frequency_limit_and_cut_cells_partition_unity() {
    let z = randles_impedance(0.0, 2.0, 30.0, 1e-5, 0.0).unwrap();
    assert!((z.re - 32.0).abs() < 1e-12 && z.im.abs() < 1e-12);
    for (l, r) in [(-0.2, 0.8), (0.5, -0.5), (1.0, 1.0)] {
        let c = cut_cell_electrode_interface(l, r, 1.0).unwrap();
        assert!((c.electrolyte_fraction + c.electrode_fraction - 1.0).abs() < 1e-15);
    }
}

#[test]
fn nernst_planck_conserves_mass_with_closed_boundaries() {
    let mut np = NernstPlanck1D::new(NernstPlanck1DConfig {
        cells: 30,
        length: 1.0,
        dt: 1e-3,
        boundary: TransportBoundary1D::NoFlux,
        ..Default::default()
    })
    .unwrap();
    np.set_potential_fn(|x| 0.01 * x);
    np.add_species_fn("A", 1, 0.01, |x| 1.0 + x).unwrap();
    let m0 = np.total_amount(0).unwrap();
    np.step(100).unwrap();
    assert!((np.total_amount(0).unwrap() - m0).abs() < 1e-12);
    let mesh = make_sheared_cartesian_hexa_mesh(4, 3, 1, 1.0, 1.0, 0.1, 0.0).unwrap();
    let mut pm = NernstPlanckPolyMesh::new(
        mesh,
        NernstPlanckPolyConfig {
            dt: 1e-3,
            ..Default::default()
        },
    )
    .unwrap();
    pm.add_species_fn("K", 1, 0.01, |p| 1.0 + 0.5 * p.x)
        .unwrap();
    pm.add_species_fn("Cl", -1, 0.02, |p| 1.0 + 0.5 * p.x)
        .unwrap();
    let (a0, b0) = (pm.total_amount(0).unwrap(), pm.total_amount(1).unwrap());
    for _ in 0..5 {
        pm.solve_electroneutral_potential(1000, 1e-12).unwrap();
        pm.step(10).unwrap();
    }
    assert!(
        (pm.total_amount(0).unwrap() - a0).abs() < 1e-12
            && (pm.total_amount(1).unwrap() - b0).abs() < 1e-12
    );
    assert!(pm.total_charge().abs() < 1e-6);
}
