//! Value-by-value comparison with the C++ upstream FVM species/reacting/flame/soot/VOF units.
//! Goldens come from `parity/oracle_fvm_reacting.cpp` built against the pinned upstream revision.
use cfd_solvers::{
    chemistry::{ArrheniusRate, ElementaryReaction, ReactionNetwork, Species, StoichiometricTerm},
    core::Vec3,
    fvm::*,
};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-fvm-reacting.json")).unwrap()
}

#[track_caller]
fn same(name: &str, rust: &[f64], tol: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= tol * c.abs().max(1.0),
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn mesh_and_flux() -> (PolyMesh, Vec<f64>) {
    let mesh = make_sheared_cartesian_hexa_mesh(6, 4, 1, 0.6, 0.4, 0.1, 0.0).unwrap();
    let flux = mesh
        .faces()
        .iter()
        .map(|f| f.area.x * 0.05 + f.area.y * 0.01)
        .collect();
    (mesh, flux)
}

#[test]
fn species_transport_matches_upstream() {
    let (mesh, flux) = mesh_and_flux();
    let mut s = MultiSpeciesTransport::new(mesh, &[1e-3, 5e-3], 0.01).unwrap();
    s.initialize(0, |p| {
        (-30.0 * ((p.x - 0.2).powi(2) + (p.y - 0.2).powi(2))).exp()
    })
    .unwrap();
    s.initialize(1, |p| 1.0 + p.x).unwrap();
    s.set_face_flux(&flux).unwrap();
    s.run(5).unwrap();
    same("species0", s.values(0).unwrap(), 1e-9);
    same("species1", s.values(1).unwrap(), 1e-9);
    same(
        "species_inventory",
        &[s.inventory(0).unwrap(), s.inventory(1).unwrap()],
        1e-9,
    );
}

#[test]
fn reacting_sources_match_upstream() {
    let sp = |n: &str, m: f64| Species {
        name: n.into(),
        molar_mass: m,
        charge: 0,
    };
    let mut net = ReactionNetwork::new(vec![sp("A", 0.002), sp("B", 0.032), sp("C", 0.034)]);
    net.add_reaction(ElementaryReaction {
        reactants: vec![
            StoichiometricTerm {
                species: 0,
                coefficient: 1.0,
            },
            StoichiometricTerm {
                species: 1,
                coefficient: 1.0,
            },
        ],
        products: vec![StoichiometricTerm {
            species: 2,
            coefficient: 1.0,
        }],
        forward: ArrheniusRate {
            pre_exponential: 2.0e3,
            temperature_exponent: 0.5,
            activation_energy: 2.0e4,
        },
        ..Default::default()
    })
    .unwrap();
    let src = reacting_cell_sources(
        &net,
        &[vec![1.0, 2.0, 0.1], vec![0.5, 0.5, 0.0]],
        &[800.0, 1200.0],
    )
    .unwrap();
    let flat: Vec<f64> = src.iter().flatten().copied().collect();
    same("reacting_sources", &flat, 1e-12);
}

#[test]
fn premixed_flame_and_models_match_upstream() {
    let mut flame = PremixedFlame1D::new(PremixedFlame1DConfig {
        cells: 200,
        length: 1.0,
        diffusivity: 1e-3,
        reaction_rate: 4.0,
        dt: 0.002,
        ..Default::default()
    })
    .unwrap();
    flame.initialize_front(0.2, 0.02).unwrap();
    flame.step(500);
    same(
        "flame",
        &[
            flame.front_location(),
            flame.theoretical_speed(),
            flame.temperature(30).unwrap(),
            flame.progress()[50],
            flame.progress()[120],
        ],
        1e-9,
    );

    let premixed = PremixedFlameModel::new(PremixedFlameModelConfig {
        reaction_rate: 2.0,
        heat_release: 5e5,
        ignition_temperature: 400.0,
    })
    .unwrap();
    let p1 = premixed
        .evaluate(FlameModelState {
            progress: 0.3,
            mixture_fraction: 0.0,
            temperature: 500.0,
        })
        .unwrap();
    let p2 = premixed
        .evaluate(FlameModelState {
            progress: 0.3,
            mixture_fraction: 0.0,
            temperature: 350.0,
        })
        .unwrap();
    let nonpremixed = NonPremixedFlameModel::new(NonPremixedFlameModelConfig::default()).unwrap();
    let n1 = nonpremixed
        .evaluate(FlameModelState {
            progress: 0.2,
            mixture_fraction: 0.45,
            temperature: 600.0,
        })
        .unwrap();
    same(
        "flame_models",
        &[
            p1.progress_source,
            p1.heat_release,
            p2.progress_source,
            n1.progress_source,
            n1.mixture_fraction_source,
            n1.heat_release,
        ],
        1e-12,
    );

    let soot = SootRadiationModel::new(SootRadiationConfig::default()).unwrap();
    let s = soot
        .evaluate(
            SootRadiationState {
                soot_mass_fraction: 0.001,
                temperature: 1500.0,
                heat_release: 2.0e5,
            },
            300.0,
        )
        .unwrap();
    same(
        "soot",
        &[s.soot_source, s.absorption_coefficient, s.radiative_source],
        1e-12,
    );
}

#[test]
fn vof_and_surface_tension_match_upstream() {
    let (mesh, flux) = mesh_and_flux();
    let mut vof = VofTransport::new(
        mesh.clone(),
        VofConfig {
            dt: 0.02,
            compression: 0.6,
        },
    )
    .unwrap();
    vof.initialize_profile(|p| if p.x < 0.25 && p.y < 0.25 { 1.0 } else { 0.0 });
    vof.set_face_flux(flux).unwrap();
    for _ in 0..10 {
        vof.step().unwrap();
    }
    same("vof_alpha", vof.alpha(), 1e-10);
    same("vof_volume", &[vof.volume()], 1e-10);
    let alpha: Vec<f64> = mesh
        .cells()
        .iter()
        .map(|c| (10.0 * (0.3 - c.center.x)).tanh())
        .collect();
    let csf: Vec<f64> = csf_surface_tension(&mesh, &alpha, 0.07)
        .unwrap()
        .iter()
        .flat_map(|v| [v.x, v.y, v.z])
        .collect();
    same("csf", &csf, 1e-10);
    let n = contact_angle_normal(Vec3::new(0.3, 1.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 1.0).unwrap();
    let n2 = contact_angle_normal(Vec3::new(0.0, 2.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 0.5).unwrap();
    same("contact_angle", &[n.x, n.y, n.z, n2.x, n2.y, n2.z], 1e-14);
}
