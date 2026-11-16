//! Value-by-value comparison with the C++ upstream reaction-path, reactor-network, speciation,
//! advanced-thermo and database-IO chemistry units. Goldens come from `parity/oracle_chemistry.cpp`.
use cfd_solvers::chemistry::*;
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-chemistry.json")).unwrap()
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

fn sp(name: &str, m: f64) -> Species {
    Species {
        name: name.into(),
        molar_mass: m,
        charge: 0,
    }
}
fn term(s: usize, c: f64) -> StoichiometricTerm {
    StoichiometricTerm {
        species: s,
        coefficient: c,
    }
}

fn make_net() -> ReactionNetwork {
    let mut n = ReactionNetwork::new(vec![sp("A", 0.03), sp("B", 0.06), sp("C", 0.06)]);
    n.add_reaction(ElementaryReaction {
        reactants: vec![term(0, 2.0)],
        products: vec![term(1, 1.0)],
        forward: ArrheniusRate {
            pre_exponential: 3.0,
            temperature_exponent: 0.0,
            activation_energy: 0.0,
        },
        ..Default::default()
    })
    .unwrap();
    n.add_reaction(ElementaryReaction {
        reactants: vec![term(1, 1.0)],
        products: vec![term(2, 1.0)],
        forward: ArrheniusRate {
            pre_exponential: 4.0,
            temperature_exponent: 0.5,
            activation_energy: 1000.0,
        },
        ..Default::default()
    })
    .unwrap();
    n
}

fn flat(e: &[ReactionPathEdge]) -> Vec<f64> {
    e.iter()
        .flat_map(|e| [e.from as f64, e.to as f64, e.molar_flux])
        .collect()
}

#[test]
fn reaction_path_and_reactor_network_match_upstream() {
    let net = make_net();
    let c = [1.0, 0.5, 0.2];
    close(
        "path",
        &flat(&reaction_path(&net, &c, 300.0, 0.0).unwrap()),
        1e-12,
        0.0,
    );
    close(
        "path_cut",
        &flat(&reaction_path(&net, &c, 300.0, 3.0).unwrap()),
        1e-12,
        0.0,
    );

    let mut w = WellStirredReactorNetwork::new(make_net());
    w.add_reactor(WellStirredReactor {
        volume: 2.0,
        temperature: 300.0,
        concentrations: vec![1.0, 0.0, 0.0],
    })
    .unwrap();
    w.add_reactor(WellStirredReactor {
        volume: 1.0,
        temperature: 350.0,
        concentrations: vec![0.0, 0.4, 0.1],
    })
    .unwrap();
    w.connect(ReactorConnection {
        first: 0,
        second: 1,
        exchange_flow: 0.3,
    })
    .unwrap();
    w.advance(0.5, 0.05, ReactorConfig::default()).unwrap();
    let wc: Vec<f64> = w
        .reactors()
        .iter()
        .flat_map(|r| r.concentrations.clone())
        .collect();
    close("wsr_conc", &wc, 1e-9, 1e-14);
    close("wsr_total", &w.total_species_amounts(), 1e-9, 1e-14);
}

#[test]
fn speciation_block_matches_upstream() {
    let b = MonolithicSpeciationBlock::new(SpeciationConfig::default()).unwrap();
    let r = b.solve(0.1, &[10f64.ln(), 0.5f64.ln(), 2.0], None).unwrap();
    close("spec_c", &r.concentrations, 1e-10, 0.0);
    close(
        "spec_meta",
        &[r.iterations as f64, f64::from(u8::from(r.converged))],
        0.0,
        0.0,
    );
    let r2 = b
        .solve(2.0, &[-1.0, 0.0, 1.0], Some(&[0.01, 0.02, 0.03, 0.04]))
        .unwrap();
    close("spec2_c", &r2.concentrations, 1e-10, 0.0);
}

#[test]
fn advanced_thermo_matches_upstream() {
    let n2 = Nasa7Polynomial {
        a: [
            3.53100528,
            -1.23660988e-4,
            -5.02999433e-7,
            2.43530612e-9,
            -1.40881235e-12,
            -1046.97628,
            2.96747038,
        ],
        minimum_temperature: 200.0,
        maximum_temperature: 1000.0,
    };
    let mut nasa = Vec::new();
    for t in [250.0, 298.15, 600.0, 999.0] {
        nasa.extend([
            n2.cp_molar(t).unwrap(),
            n2.h_molar(t).unwrap(),
            n2.s_molar(t).unwrap(),
        ]);
    }
    // h at 298.15 K is a near-cancellation (~4e-6 J/mol): absolute tolerance only.
    close("nasa", &nasa, 1e-12, 1e-8);
    let h2o = ShomatePolynomial {
        A: 30.092,
        B: 6.832514,
        C: 6.793435,
        D: -2.53448,
        E: 0.082139,
        F: -250.881,
        G: 223.3967,
        H: -241.8264,
        minimum_temperature: 500.0,
        maximum_temperature: 1700.0,
    };
    let mut s = Vec::new();
    for t in [500.0, 1000.0, 1700.0] {
        s.extend([
            h2o.cp_molar(t).unwrap(),
            h2o.h_molar_kj(t).unwrap(),
            h2o.s_molar(t).unwrap(),
        ]);
    }
    close("shomate", &s, 1e-12, 0.0);
    let nacl = BinaryPitzer11 {
        beta0: 0.0765,
        beta1: 0.2664,
        c_phi: 0.00127,
    };
    let pz: Vec<f64> = [0.0, 1e-10, 0.1, 1.0, 3.0]
        .iter()
        .map(|m| nacl.mean_activity_coefficient(*m, 298.15).unwrap())
        .collect();
    close("pitzer11", &pz, 1e-12, 0.0);

    let z = [1, -1, 2, -2, 0];
    let m = [0.5, 0.6, 0.1, 0.05, 0.2];
    let pi = |first, second, beta0, beta1, c_phi| PitzerInteraction {
        first,
        second,
        beta0,
        beta1,
        c_phi,
    };
    let mp = MulticomponentPitzer::new(
        vec![
            pi(0, 1, 0.0765, 0.2664, 0.00127),
            pi(2, 3, 0.2, 2.9, 0.01),
            pi(0, 3, 0.02, 1.1, 0.0),
        ],
        0.392,
        1.2,
    );
    close(
        "mpitzer",
        &mp.activity_coefficients(&z, &m, 298.15).unwrap(),
        1e-12,
        0.0,
    );
    close(
        "davies",
        &DaviesActivity::default()
            .activity_coefficients(&z, &m, 298.15)
            .unwrap(),
        1e-12,
        0.0,
    );
    let sit = SitActivity::new(
        vec![
            SitPair {
                first: 0,
                second: 1,
                epsilon: 0.03,
            },
            SitPair {
                first: 2,
                second: 3,
                epsilon: -0.1,
            },
            SitPair {
                first: 4,
                second: 4,
                epsilon: 0.05,
            },
        ],
        0.509,
    );
    close(
        "sit",
        &sit.activity_coefficients(&z, &m, 298.15).unwrap(),
        1e-12,
        0.0,
    );

    let mut g = Vec::new();
    for (a, b, t) in [
        (-1000.0, 2000.0, 300.0),
        (0.0, 2.0e5, 300.0),
        (3.0e5, 0.0, 400.0),
        (5.0, 5.0, 1000.0),
    ] {
        let r = minimize_binary_ideal_gibbs(a, b, t).unwrap();
        g.extend([r.mole_fraction_a, r.mole_fraction_b, r.gibbs_molar]);
    }
    close("gibbs", &g, 1e-12, 1e-14);
    close(
        "sorption",
        &[
            langmuir_surface_coverage(2.5, 0.3).unwrap(),
            langmuir_surface_coverage(0.0, 1.0).unwrap(),
            binary_ion_exchange_fraction(3.0, 0.2, 0.1).unwrap(),
            binary_ion_exchange_fraction(3.0, 0.0, 0.1).unwrap(),
            binary_ion_exchange_fraction(3.0, 0.2, 0.0).unwrap(),
        ],
        1e-14,
        0.0,
    );
}

#[test]
fn database_readers_match_upstream() {
    let phreeqc = "# test db\nSOLUTION_MASTER_SPECIES\nNa  Na+  22.99  1   # sodium\nCl Cl- 35.453 -1\n  Ca   Ca+2   0.04008  2\nbad row\nSOLUTION_SPECIES\nX X 12 0\n";
    let s = load_phreeqc_master_species(phreeqc).unwrap();
    assert_eq!(
        s.iter().map(|s| s.name.as_str()).collect::<Vec<_>>(),
        ["Na+", "Cl-", "Ca+2"]
    );
    let pdb: Vec<f64> = s
        .iter()
        .flat_map(|s| [s.molar_mass, f64::from(s.charge)])
        .collect();
    close("phreeqc", &pdb, 1e-14, 0.0);

    let yaml = "units: {length: m}\nspecies:\n- name: H2\n- name: O2\n- name: H2O\n- name: OH\nreactions:\n- equation: 2 H2 + O2 => 2H2O\n  rate-constant: {A: 1.5e3, b: 0.5, Ea: 2000.0}\n- equation: H2O <=> OH + 0.5 H2   # comment\n  rate-constant:\n    A: 7.0\n    b: -0.25\n    Ea: 500\n";
    let mech = load_cantera_yaml(yaml).unwrap();
    close(
        "cantera_rates",
        &mech.reaction_rates(&[0.4, 0.3, 0.2, 0.1], 800.0).unwrap(),
        1e-12,
        0.0,
    );
    let mut cst = Vec::new();
    for r in mech.reactions() {
        for t in &r.reactants {
            cst.extend([t.species as f64, t.coefficient]);
        }
        cst.push(-1.0);
        for t in &r.products {
            cst.extend([t.species as f64, t.coefficient]);
        }
        cst.extend([
            r.forward.pre_exponential,
            r.forward.temperature_exponent,
            r.forward.activation_energy,
        ]);
    }
    close("cantera_terms", &cst, 1e-14, 0.0);
}
