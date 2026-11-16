//! Value-by-value comparison with the C++ upstream porous, aqueous-equilibrium and
//! electrochemistry units. Goldens come from `parity/oracle_porous_echem.cpp`.
use cfd_solvers::{
    chemistry::*, core::Vec3, electrochemistry::*, fvm::make_sheared_cartesian_hexa_mesh, porous::*,
};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-porous-echem.json")).unwrap()
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

fn cat(mut a: Vec<f64>, b: &[f64]) -> Vec<f64> {
    a.extend_from_slice(b);
    a
}
fn f(b: bool) -> f64 {
    f64::from(u8::from(b))
}

#[test]
fn porous_matches_upstream() {
    let mut th = Thmc1D::new(Thmc1DConfig {
        cells: 12,
        length: 2.0,
        pressure_diffusivity: 0.02,
        thermal_diffusivity: 0.015,
        species_diffusivity: 0.01,
        reaction_rate: 0.3,
        thermo_pressure_coupling: 0.05,
        reaction_heat: 0.2,
    })
    .unwrap();
    th.initialize(1.0, 300.0, 0.5).unwrap();
    th.set_cell(3, 1.5, 320.0, 0.9).unwrap();
    th.set_cell(8, 0.8, 290.0, 0.0).unwrap();
    for _ in 0..30 {
        th.step(0.05).unwrap();
    }
    let d = th.diagnostics();
    let v = cat(
        cat(
            cat(
                cat(
                    vec![
                        d.pressure_mean,
                        d.temperature_mean,
                        d.species_inventory,
                        d.reaction_heat,
                    ],
                    th.pressure(),
                ),
                th.temperature(),
            ),
            th.species(),
        ),
        th.displacement(),
    );
    close("thmc", &v, 1e-12, 1e-15);
    let mut ac = AllenCahn1D::new(AllenCahn1DConfig {
        cells: 40,
        length: 1.0,
        mobility: 1.5,
        interface_width: 0.05,
    })
    .unwrap();
    ac.initialize_tanh(0.4, 0.2).unwrap();
    let e0 = ac.free_energy();
    for _ in 0..40 {
        ac.step(2e-4).unwrap();
    }
    close(
        "allen_cahn",
        &cat(vec![e0, ac.free_energy()], ac.phase()),
        1e-12,
        1e-15,
    );
    let mut fr = FracturedDarcyThmcNetwork::new(FracturedDarcyThmcConfig {
        cells: 5,
        ..Default::default()
    })
    .unwrap();
    for (a, b, t, ap) in [
        (0, 1, 1.0, 1.0),
        (1, 2, 0.5, 4.0),
        (2, 3, 2.0, 0.5),
        (3, 4, 1.0, 1.0),
        (4, 0, 0.3, 2.0),
    ] {
        fr.add_edge(FracturedDarcyEdge {
            first: a,
            second: b,
            transmissibility: t,
            aperture: ap,
        })
        .unwrap();
    }
    fr.initialize(1.0, 0.3, 300.0, 0.2).unwrap();
    fr.set_cell(0, 3.0, 0.9, 350.0, 1.0).unwrap();
    fr.set_cell(3, 0.2, 0.0, 280.0, 0.0).unwrap();
    for _ in 0..20 {
        fr.step(0.02).unwrap();
    }
    let fd = fr.diagnostics();
    let v = cat(
        cat(
            cat(
                cat(
                    cat(
                        vec![
                            fd.pressure_mass,
                            fd.saturation_mass,
                            fd.temperature_mass,
                            fd.species_mass,
                            fd.absolute_flux,
                            fd.reaction_consumption,
                        ],
                        fr.pressure(),
                    ),
                    fr.saturation(),
                ),
                fr.temperature(),
            ),
            fr.species(),
        ),
        fr.displacement(),
    );
    close("fractured", &v, 1e-12, 1e-15);
    let pts = [
        (0.0, 0.0),
        (0.2, -0.05),
        (0.5, -0.02),
        (0.8, -0.06),
        (1.0, 0.0),
    ]
    .map(|(composition, free_energy)| CalphadBinaryPoint {
        composition,
        free_energy,
    });
    let table = BinaryCalphadTable::new(pts.to_vec()).unwrap();
    let mut ch = CahnHilliard1D::new(
        CahnHilliard1DConfig {
            cells: 16,
            length: 1.0,
            mobility: 0.5,
            gradient_energy: 0.002,
            anisotropy_factor: 1.2,
            nucleation_amplitude: 0.0,
        },
        None,
    )
    .unwrap();
    ch.initialize(0.45, 0.05).unwrap();
    ch.seed_nucleus(4, 1.5, 0.1).unwrap();
    for _ in 0..25 {
        ch.step(1e-4).unwrap();
    }
    let mut ct = CahnHilliard1D::new(
        CahnHilliard1DConfig {
            cells: 16,
            length: 1.0,
            mobility: 0.5,
            gradient_energy: 0.002,
            anisotropy_factor: 1.0,
            nucleation_amplitude: 0.0,
        },
        Some(table.clone()),
    )
    .unwrap();
    ct.initialize(0.5, 0.1).unwrap();
    for _ in 0..25 {
        ct.step(1e-4).unwrap();
    }
    let head = vec![
        ch.mean(),
        ch.free_energy().unwrap(),
        ct.mean(),
        ct.free_energy().unwrap(),
        table.free_energy(0.35).unwrap(),
        table.chemical_potential(0.9).unwrap(),
    ];
    close(
        "cahn_hilliard",
        &cat(
            cat(
                cat(cat(head, ch.composition()), ch.chemical_potential()),
                ct.composition(),
            ),
            ct.chemical_potential(),
        ),
        1e-11,
        1e-14,
    );
}

fn families() -> Vec<AcidFamily> {
    vec![
        AcidFamily {
            total_mol_per_litre: 0.01,
            fully_protonated_charge: 0,
            dissociation_constants: vec![1.8e-5],
        },
        AcidFamily {
            total_mol_per_litre: 0.002,
            fully_protonated_charge: 0,
            dissociation_constants: vec![7.5e-3, 6.2e-8, 4.8e-13],
        },
    ]
}

#[test]
fn aqueous_and_electrochemistry_basics_match_upstream() {
    let aq = equilibrate_acids(&families(), 0.003, 1e-14).unwrap();
    let mut v = vec![aq.ph, aq.hydrogen_mol_per_litre, aq.hydroxide_mol_per_litre];
    v.push(aq.charge_residual_mol_per_litre);
    for fam in &aq.family_concentrations {
        v.extend(fam);
    }
    close("aqueous", &v, 1e-11, 1e-16);

    let mut b = vec![
        nernst_potential(0.77, 310.0, 1.0, 3.5).unwrap(),
        butler_volmer_current_density(2.0, 0.05, 300.0, 2.0, 0.4, 0.6).unwrap(),
        faradaic_molar_flux(10.0, 2.0).unwrap(),
        corrosion_penetration_rate(5.0, 0.0558, 7874.0, 2.0).unwrap(),
    ];
    let ip = solve_implicit_pnp_electrode_cell(
        0.8,
        0.02,
        &ImplicitPnpElectrodeConfig::default(),
        40,
        1e-10,
    )
    .unwrap();
    b.extend([
        ip.concentration_mol_m3,
        ip.potential_v,
        ip.iterations as f64,
        f(ip.converged),
    ]);
    let cr = solve_corrosion_cell_1d(&CorrosionCell1DConfig {
        metal_potential: 0.12,
        bulk_electrolyte_potential: -0.01,
        equilibrium_potential: 0.02,
        ..Default::default()
    })
    .unwrap();
    b.extend([
        cr.current_density,
        cr.surface_electrolyte_potential,
        cr.overpotential,
        cr.molar_dissolution_flux,
        cr.penetration_rate,
        cr.iterations as f64,
    ]);
    let r = |name: &str, e: f64, i0: f64, area: f64, n: f64| ElectrodeReaction {
        name: name.into(),
        equilibrium_potential: e,
        exchange_current_density: i0,
        area,
        electrons: n,
        ..Default::default()
    };
    let mp = solve_mixed_potential(
        &[
            r("Fe", -0.44, 1e-3, 1.0, 2.0),
            r("O2", 0.4, 1e-6, 2.0, 4.0),
            r("H2", 0.0, 1e-4, 0.5, 2.0),
        ],
        298.15,
        0.0,
        200,
        1e-10,
    )
    .unwrap();
    b.extend([mp.potential, mp.total_current, mp.iterations as f64]);
    b.extend(&mp.reaction_current_density);
    b.extend(&mp.reaction_current);
    close("echem_basic", &b, 1e-12, 1e-18);
}

#[test]
fn corrosion_models_match_upstream() {
    let t = |a: f64, p: f64| RedoxActivityTerm {
        activity: a,
        stoichiometric_power: p,
    };
    let (prod, reac) = ([t(0.5, 2.0)], [t(0.1, 1.0), t(2.0, 3.0)]);
    let mut m = vec![
        reaction_quotient(&prod, &reac).unwrap(),
        concentration_dependent_nernst(0.3, 300.0, 2.0, &prod, &reac).unwrap(),
    ];
    let orr = oxygen_reduction_cathodic(0.6, 1.229, 0.2, 1e-4, 1e-6, 5e-3, 298.15, 4.0).unwrap();
    let her = hydrogen_evolution_cathodic(-0.5, 1.0, 1e-3, 1e-3, 2.0, 298.15, 2.0).unwrap();
    m.extend([
        orr.equilibrium_potential,
        orr.activation_current_density,
        orr.limited_current_density,
        her.equilibrium_potential,
        her.activation_current_density,
        her.limited_current_density,
    ]);
    m.push(charge_transfer_resistance(2.0, 300.0, 2.0, 0.5, 0.5).unwrap());
    let (z1, z0) = (
        randles_impedance(100.0, 1.0, 20.0, 1e-4, 5.0).unwrap(),
        randles_impedance(0.0, 1.0, 20.0, 1e-4, 5.0).unwrap(),
    );
    m.extend([z1.re, z1.im, z0.re, z0.im]);
    let mut dl = DoubleLayerState {
        capacitance_per_area: 0.2,
        potential: 0.01,
    };
    dl.advance_current(3.0, 0.01).unwrap();
    m.extend([dl.potential, dl.charge_density()]);
    let pm = PassivationModel {
        active_equilibrium_potential: -0.4,
        exchange_current_density: 1e-2,
        passive_onset: -0.2,
        transpassive_onset: 0.8,
        passive_current_density: 1e-3,
        transpassive_slope: 0.05,
        electrons: 2.0,
    };
    for e in [-0.45, -0.3, 0.2, 1.0] {
        m.push(passivation_current_density(e, &pm, 298.15).unwrap());
    }
    let mut pl = ProductLayerState {
        thickness: 1e-6,
        resistivity: 1e5,
        growth_fraction: 0.6,
    };
    pl.advance_from_anodic_current(4.0, 100.0, 0.0558, 7874.0, 2.0)
        .unwrap();
    m.extend([pl.thickness, pl.area_specific_resistance()]);
    m.extend([
        faradaic_recession_distance(-3.0, 10.0, 0.0558, 7874.0, 2.0).unwrap(),
        cathodic_protection_current_density(-0.85, -0.6, 50.0).unwrap(),
        bruggeman_effective_transport(2e-9, 0.3, 1.5).unwrap(),
        solid_phase_ohmic_drop(10.0, 1e-3, 100.0).unwrap(),
        marcus_current_density(1.0, 0.1, 6e4, 298.15, 1.0).unwrap(),
        marcus_current_density(1.0, -0.2, 6e4, 320.0, 2.0).unwrap(),
        bikerman_activity_correction(1.0, 0.02, 10.0).unwrap(),
    ]);
    let mf = modified_nernst_planck_flux(1.0, 1.5, 1e-9, 20.0, -1, 300.0, 1e-4, 0.01, 20.0, 1e-6)
        .unwrap();
    m.extend([mf.molar_flux, mf.activity]);
    let sw = smooth_electrode_interface(0.3, 0.5).unwrap();
    m.extend([sw.electrolyte_fraction, sw.delta]);
    for (l, r) in [(-1.0, -2.0), (1.0, 2.0), (-0.3, 0.7), (0.2, -0.6)] {
        let c = cut_cell_electrode_interface(l, r, 1.0).unwrap();
        m.extend([
            c.electrolyte_fraction,
            c.electrode_fraction,
            c.interface_area_fraction,
        ]);
    }
    let bv =
        solve_anodic_bv_mass_transfer(10.0, 1e-4, -0.44, -0.3, 1e-2, 298.15, 2.0, 1.0).unwrap();
    m.extend([
        bv.surface_concentration,
        bv.current_density,
        bv.molar_flux,
        bv.iterations as f64,
    ]);
    m.push(acid_speciation_redox_potential(0.8, 298.15, 2.0, 2.0, &families(), 0.003).unwrap());
    let crit = PittingCriterion {
        critical_chloride: 0.3,
        critical_potential: 0.1,
    };
    m.extend([
        f(pitting_initiates(0.5, 0.2, &crit).unwrap()),
        f(pitting_initiates(0.5, 0.0, &crit).unwrap()),
    ]);
    m.extend([
        stress_assisted_exchange_current(1e-3, 2e8, 1e-6, 298.15).unwrap(),
        phase_field_corrosion_step(0.3, 2.0, 0.5, 1.5, 0.01, 0.1).unwrap(),
        recessed_level_set(-0.2, 0.05).unwrap(),
        marcus_hush_chidsey_current_density(1.0, 0.1, 6e4, 298.15, 1.0, 1200).unwrap(),
        marcus_hush_chidsey_current_density(1.0, -0.05, 4e4, 310.0, 1.0, 400).unwrap(),
    ]);
    m.extend(recession_remesh_1d(&[0.0, 0.1, 0.3, 0.7, 1.0], 0.2).unwrap());
    let sc =
        coupled_stress_corrosion_step(1e-3, 1e8, 2e-6, 300.0, 50.0, 0.0558, 7874.0, 2.0).unwrap();
    m.extend([sc.stress_amplified_current_density, sc.recession_distance]);
    let mut p2 = PhaseFieldCorrosion2D::new(6, 5, 0.1, 0.12, 1.0, 0.002).unwrap();
    p2.initialize_interface(0.2, 0.25, 0.15).unwrap();
    let drv: Vec<f64> = (0..30).map(|i| 0.1 * (i as f64).sin()).collect();
    for _ in 0..4 {
        p2.step(0.05, &drv).unwrap();
    }
    let mut p1 = PhaseFieldCorrosion1D::new(8, 0.1, 1.0, 0.003).unwrap();
    p1.initialize_interface(0.35, 0.1).unwrap();
    for _ in 0..4 {
        p1.step(0.05, &drv[..8]).unwrap();
    }
    close(
        "corrosion_models",
        &cat(cat(m, p2.order_parameter()), p1.order_parameter()),
        1e-11,
        1e-16,
    );
}

#[test]
fn nernst_planck_matches_upstream() {
    let mut np = NernstPlanck1D::new(NernstPlanck1DConfig {
        cells: 20,
        length: 1e-3,
        temperature: 298.15,
        dt: 2e-4,
        advection_velocity: 1e-4,
        boundary: TransportBoundary1D::NoFlux,
    })
    .unwrap();
    np.set_potential_fn(|x| 0.02 * (3000.0 * x).sin());
    np.add_species_fn("Na", 1, 1.3e-9, |x| {
        1.0 + 0.5 * (-1e7 * (x - 3e-4) * (x - 3e-4)).exp()
    })
    .unwrap();
    np.add_species("Cl", -1, 2.0e-9, 1.2).unwrap();
    np.step(40).unwrap();
    let mut npp = NernstPlanck1D::new(NernstPlanck1DConfig {
        cells: 16,
        length: 1.0,
        temperature: 300.0,
        dt: 1e-3,
        advection_velocity: -0.05,
        boundary: TransportBoundary1D::Periodic,
    })
    .unwrap();
    npp.add_species_fn("A", 2, 0.01, |x| {
        1.0 + 0.3 * (2.0 * std::f64::consts::PI * x).cos()
    })
    .unwrap();
    npp.set_potential_fn(|x| 0.001 * x);
    npp.step(20).unwrap();
    let head = vec![
        np.total_amount(0).unwrap(),
        np.total_amount(1).unwrap(),
        np.minimum_concentration(1).unwrap(),
        np.time(),
    ];
    close(
        "np1d",
        &cat(
            cat(
                cat(head, &np.species()[0].concentration),
                &np.species()[1].concentration,
            ),
            &npp.species()[0].concentration,
        ),
        1e-12,
        0.0,
    );

    let mesh = make_sheared_cartesian_hexa_mesh(5, 3, 1, 1e-3, 6e-4, 1e-4, 0.1).unwrap();
    let mut poly = vec![];
    for scheme in [
        ElectromigrationFluxScheme::Centered,
        ElectromigrationFluxScheme::ScharfetterGummel,
    ] {
        let mut pm = NernstPlanckPolyMesh::new(
            mesh.clone(),
            NernstPlanckPolyConfig {
                temperature: 298.15,
                dt: 2e-3,
                relative_permittivity: 78.5,
                electromigration_scheme: scheme,
            },
        )
        .unwrap();
        let a = pm
            .add_species_fn("Cu", 2, 7e-10, |p: Vec3| 10.0 + 4000.0 * p.x)
            .unwrap();
        let s2 = pm.add_species("SO4", -2, 1.0e-9, 10.5).unwrap();
        pm.set_species_butler_volmer_boundary(a, "left", 0.05, 0.0, 2.0, 2, 1.0, 0.5, 0.5)
            .unwrap();
        pm.set_species_boundary(a, "right", SpeciesBoundaryType::FixedConcentration, 12.0)
            .unwrap();
        pm.set_species_boundary(s2, "top", SpeciesBoundaryType::FixedMolarFlux, -1e-6)
            .unwrap();
        pm.set_potential_boundary("right", PotentialBoundaryType::FixedPotential, 0.0)
            .unwrap();
        pm.set_potential_fn(|p| -5.0 * p.x).unwrap();
        let flux = mesh
            .faces()
            .iter()
            .map(|f| {
                if f.neighbour.is_none() {
                    0.0
                } else {
                    1e-12 * f.area.x
                }
            })
            .collect();
        pm.set_face_flux(flux).unwrap();
        pm.step(5).unwrap();
        let en = pm.solve_electroneutral_potential(1000, 1e-10).unwrap();
        poly.extend([
            en.iterations as f64,
            pm.total_amount(a).unwrap(),
            pm.total_amount(s2).unwrap(),
            pm.minimum_concentration(a).unwrap(),
            pm.total_charge(),
        ]);
        poly.extend(&pm.species()[0].concentration);
        poly.extend(&pm.species()[1].concentration);
        poly.extend(pm.potential());
        poly.extend(pm.current_flux_faces().unwrap());
    }
    let mut ps = NernstPlanckPolyMesh::new(
        mesh,
        NernstPlanckPolyConfig {
            temperature: 298.15,
            dt: 1e-3,
            relative_permittivity: 78.5,
            electromigration_scheme: ElectromigrationFluxScheme::ScharfetterGummel,
        },
    )
    .unwrap();
    ps.add_species_fn("K", 1, 2e-9, |p| 1.0 + 1e-11 * (5000.0 * p.x).sin())
        .unwrap();
    ps.add_species("Cl", -1, 2e-9, 1.0).unwrap();
    ps.set_potential_boundary("left", PotentialBoundaryType::FixedPotential, 1e-3)
        .unwrap();
    ps.set_potential_boundary("right", PotentialBoundaryType::FixedPotential, 0.0)
        .unwrap();
    let pr = ps.solve_poisson_potential(1000, 1e-10).unwrap();
    ps.step_poisson_nernst_planck(2, 1000, 1e-10).unwrap();
    poly.extend([pr.iterations as f64, ps.time()]);
    poly.extend(ps.potential());
    poly.extend(&ps.species()[0].concentration);
    // CG potentials converge to 1e-10 relative; concentrations agree to the same order.
    close("np_poly", &poly, 1e-7, 1e-14);
}
