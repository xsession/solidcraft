//! Analytic and robustness checks for reaction paths, reactor networks, speciation, advanced
//! thermodynamics and database readers.
use cfd_solvers::chemistry::*;

fn sp(name: &str) -> Species {
    Species {
        name: name.into(),
        molar_mass: 0.03,
        charge: 0,
    }
}

#[test]
fn pure_mixing_conserves_amounts_and_relaxes_to_uniform() {
    let net = ReactionNetwork::new(vec![sp("A"), sp("B")]);
    let mut w = WellStirredReactorNetwork::new(net);
    w.add_reactor(WellStirredReactor {
        volume: 1.0,
        temperature: 300.0,
        concentrations: vec![2.0, 0.0],
    })
    .unwrap();
    w.add_reactor(WellStirredReactor {
        volume: 3.0,
        temperature: 300.0,
        concentrations: vec![0.0, 1.0],
    })
    .unwrap();
    w.connect(ReactorConnection {
        first: 0,
        second: 1,
        exchange_flow: 0.5,
    })
    .unwrap();
    let before = w.total_species_amounts();
    w.advance(40.0, 0.1, ReactorConfig::default()).unwrap();
    let after = w.total_species_amounts();
    for (a, b) in before.iter().zip(&after) {
        assert!((a - b).abs() < 1e-12);
    }
    // Uniform concentration = total amount / total volume.
    for r in w.reactors() {
        assert!(
            (r.concentrations[0] - 0.5).abs() < 1e-6 && (r.concentrations[1] - 0.75).abs() < 1e-6
        );
    }
    assert!(w
        .connect(ReactorConnection {
            first: 0,
            second: 0,
            exchange_flow: 1.0
        })
        .is_err());
    assert!(w
        .add_reactor(WellStirredReactor {
            volume: 1.0,
            temperature: 300.0,
            concentrations: vec![1.0]
        })
        .is_err());
}

#[test]
fn speciation_satisfies_mass_action_and_balance() {
    let b = MonolithicSpeciationBlock::new(SpeciationConfig::default()).unwrap();
    let logk = [1.5, -0.7];
    let r = b.solve(3.0, &logk, None).unwrap();
    assert!(r.converged);
    let c = &r.concentrations;
    assert!((c.iter().sum::<f64>() - 3.0).abs() < 1e-10);
    for (i, k) in logk.iter().enumerate() {
        assert!(((c[i + 1] / c[0]).ln() - k).abs() < 1e-9);
    }
    assert!(b.solve(-1.0, &logk, None).is_err());
    assert!(b.solve(1.0, &[], None).is_err());
    assert!(b.solve(1.0, &logk, Some(&[1.0, 0.0, 1.0])).is_err());
}

#[test]
fn reaction_path_flux_is_rate_times_stoichiometry() {
    let mut n = ReactionNetwork::new(vec![sp("A"), sp("B")]);
    n.add_reaction(ElementaryReaction {
        reactants: vec![StoichiometricTerm {
            species: 0,
            coefficient: 2.0,
        }],
        products: vec![StoichiometricTerm {
            species: 1,
            coefficient: 3.0,
        }],
        forward: ArrheniusRate {
            pre_exponential: 1.0,
            temperature_exponent: 0.0,
            activation_energy: 0.0,
        },
        ..Default::default()
    })
    .unwrap();
    let c = [0.5, 0.1];
    let rate = n.reaction_rates(&c, 300.0).unwrap()[0];
    let e = reaction_path(&n, &c, 300.0, 0.0).unwrap();
    assert_eq!(e.len(), 1);
    assert!((e[0].molar_flux - 6.0 * rate).abs() < 1e-15);
    assert!(reaction_path(&n, &c, 300.0, 10.0).unwrap().is_empty());
    assert!(reaction_path(&n, &c, 300.0, -1.0).is_err());
}

#[test]
fn thermo_polynomials_reduce_to_constant_cp() {
    let mut a = [0.0; 7];
    a[0] = 3.5;
    a[5] = -1000.0;
    a[6] = 4.0;
    let p = Nasa7Polynomial {
        a,
        ..Default::default()
    };
    let t = 500.0;
    assert!((p.cp_molar(t).unwrap() - 3.5 * GAS_CONSTANT).abs() < 1e-12);
    assert!((p.h_molar(t).unwrap() - GAS_CONSTANT * (3.5 * t - 1000.0)).abs() < 1e-9);
    assert!((p.s_molar(t).unwrap() - GAS_CONSTANT * (3.5 * t.ln() + 4.0)).abs() < 1e-12);
    assert!(p.cp_molar(100.0).is_err());
    let s = ShomatePolynomial {
        A: 29.0,
        ..Default::default()
    };
    assert!((s.cp_molar(800.0).unwrap() - 29.0).abs() < 1e-14);
    assert!(s.cp_molar(10.0).is_err());
}

#[test]
fn activity_and_equilibrium_limits() {
    assert_eq!(
        BinaryPitzer11 {
            beta0: 0.1,
            beta1: 0.2,
            c_phi: 0.0
        }
        .mean_activity_coefficient(0.0, 298.0)
        .unwrap(),
        1.0
    );
    let d = DaviesActivity::default()
        .activity_coefficients(&[1, -1, 0], &[0.0, 0.0, 1.0], 298.0)
        .unwrap();
    assert!(d.iter().all(|g| (g - 1.0).abs() < 1e-15));
    assert!(DaviesActivity::default()
        .activity_coefficients(&[1], &[0.1, 0.2], 298.0)
        .is_err());
    assert!(MulticomponentPitzer::new(
        vec![PitzerInteraction {
            first: 5,
            second: 0,
            ..Default::default()
        }],
        0.392,
        1.2
    )
    .activity_coefficients(&[1, -1], &[0.1, 0.1], 298.0)
    .is_err());
    let g = minimize_binary_ideal_gibbs(100.0, 100.0, 300.0).unwrap();
    assert!((g.mole_fraction_a - 0.5).abs() < 1e-15);
    assert!((g.gibbs_molar - (100.0 - GAS_CONSTANT * 300.0 * 2f64.ln())).abs() < 1e-9);
    assert!(minimize_binary_ideal_gibbs(0.0, 0.0, -1.0).is_err());
    assert!((langmuir_surface_coverage(1.0, 1.0).unwrap() - 0.5).abs() < 1e-15);
    assert!(binary_ion_exchange_fraction(0.0, 1.0, 1.0).is_err());
}

#[test]
fn database_readers_reject_malformed_input_without_panicking() {
    assert!(load_phreeqc_master_species("").is_err());
    assert!(load_phreeqc_master_species("SOLUTION_MASTER_SPECIES\nNa Na+ -1 1\n").is_err());
    assert!(load_cantera_yaml("reactions:\n- equation: A => B\n").is_err());
    assert!(load_cantera_yaml("species:\n- name: A\n- name: A\n").is_err());
    assert!(load_cantera_yaml("species:\n- name: A\nreactions:\n- equation: A = A\n").is_err());
    assert!(load_cantera_yaml("species:\n- name: A\nreactions:\n- equation: -2 A => A\n").is_err());
    let n = load_cantera_yaml("species:\n- name: A\n- name: B\nreactions:\n- equation: A <=> 2B\n")
        .unwrap();
    assert_eq!(n.reactions()[0].products[0].coefficient, 2.0);
    // Non-ASCII content is handled without byte-slicing panics.
    assert!(
        load_cantera_yaml("species:\n- name: Ä\nreactions:\n- equation: Ä => Ä\n  A: é5\n").is_ok()
    );

    let adapter = ReaktoroDynamicAdapter::new(Box::new(|c: &[f64], _, _| {
        c.iter().map(|v| v * 0.5).collect()
    }));
    assert_eq!(
        adapter.equilibrate(&[2.0, 4.0], 298.0, 1e5).unwrap(),
        vec![1.0, 2.0]
    );
    let bad = ReaktoroDynamicAdapter::new(Box::new(|_: &[f64], _, _| vec![-1.0]));
    assert!(bad.equilibrate(&[1.0], 298.0, 1e5).is_err());
    assert!(adapter.equilibrate(&[1.0], 0.0, 1e5).is_err());
}

#[test]
fn case_folding_with_length_changing_characters_does_not_panic() {
    // 'İ' lowercases to two chars under Unicode rules; slicing must stay on char boundaries.
    assert!(load_cantera_yaml(
        "species:\n- name: A\nreactions:\n- equation: A => A\n  İİİ b: 2 A: 3\n"
    )
    .is_ok());
}
