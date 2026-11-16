//! Value-by-value comparison with the C++ upstream fvm adaptive mesh, combustion thermo, mesh
//! motion and reordering units. Goldens come from `parity/oracle_fvm_amr.cpp`.
use cfd_solvers::{core::Vec3, fvm::*};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-fvm-amr.json")).unwrap()
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
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= rel * c.abs() + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn dump(m: &AdaptiveHexMesh) -> Vec<f64> {
    m.cells()
        .iter()
        .flat_map(|c| {
            [
                c.lower.x,
                c.lower.y,
                c.lower.z,
                c.upper.x,
                c.upper.y,
                c.upper.z,
                c.level as f64,
                c.root as f64,
                c.lineage as f64,
            ]
        })
        .collect()
}
fn b(v: bool) -> f64 {
    f64::from(u8::from(v))
}

#[test]
fn thermo_motion_and_reordering_match_upstream() {
    let sp = [
        IdealGasSpeciesThermo {
            molar_mass: 0.032,
            cp_mass: 918.0,
            formation_enthalpy_mass: 0.0,
        },
        IdealGasSpeciesThermo {
            molar_mass: 0.028,
            cp_mass: 1040.0,
            formation_enthalpy_mass: 0.0,
        },
        IdealGasSpeciesThermo {
            molar_mass: 0.018,
            cp_mass: 1996.0,
            formation_enthalpy_mass: -1.3435e7,
        },
    ];
    let y = [0.2, 0.7, 0.1];
    let st = ideal_gas_mixture_state(1200.0, 2e5, &y, &sp).unwrap();
    close(
        "combustion",
        &[
            st.density,
            st.mixture_molar_mass,
            st.gas_constant,
            st.cp,
            st.cv,
            st.gamma,
            st.specific_enthalpy,
            temperature_from_mixture_enthalpy(st.specific_enthalpy, &y, &sp).unwrap(),
        ],
        1e-14,
        0.0,
    );
    let mesh = make_cartesian_hexa_mesh(3, 2, 2, 1.5, 1.0, 0.8).unwrap();
    let moved = rigidly_moved_mesh(
        &mesh,
        Vec3::new(0.1, -0.2, 0.3),
        Vec3::new(0.3, 0.4, 1.0),
        0.35,
    )
    .unwrap();
    let mut mv: Vec<f64> = moved
        .cells()
        .iter()
        .flat_map(|c| [c.center.x, c.center.y, c.center.z])
        .collect();
    mv.extend(
        moved
            .faces()
            .iter()
            .flat_map(|f| [f.area.x, f.area.y, f.area.z]),
    );
    mv.extend(ale_mesh_flux(&mesh, &moved, 0.01).unwrap());
    let field: Vec<f64> = (0..mesh.cells().len())
        .map(|i| (1.0 + i as f64).sin())
        .collect();
    let fine = make_cartesian_hexa_mesh(4, 3, 2, 1.5, 1.0, 0.8).unwrap();
    mv.extend(conservative_nearest_remap(&mesh, &field, &fine).unwrap());
    close("motion", &mv, 1e-13, 1e-15);
    let sh = make_sheared_cartesian_hexa_mesh(4, 3, 2, 1.0, 1.0, 1.0, 0.3).unwrap();
    let order = morton_cell_order(&sh);
    let re = reorder_cells(&sh, &order).unwrap();
    let mut ro: Vec<f64> = order.iter().map(|v| *v as f64).collect();
    for f in re.faces() {
        ro.extend([f.owner as f64, f.neighbour.map_or(-1.0, |n| n as f64)]);
    }
    close("reorder", &ro, 0.0, 0.0);
}

#[test]
fn adaptive_mesh_matches_upstream() {
    let amr = AdaptiveHexMesh::cartesian(3, 2, 2, 1.5, 1.0, 1.0).unwrap();
    let mut mk = vec![false; amr.cell_count()];
    (mk[0], mk[7]) = (true, true);
    let r1 = MarkedHexRefinement::new(mk, 3, true)
        .unwrap()
        .apply(&amr)
        .unwrap();
    let mut mk2 = vec![false; r1.mesh.cell_count()];
    mk2[1] = true;
    let r2 = MarkedHexRefinement::new(mk2, 3, true)
        .unwrap()
        .apply(&r1.mesh)
        .unwrap();
    let mut a = dump(&r2.mesh);
    a.push(b(satisfies_two_to_one(&r2.mesh)));
    for ws in &r2.mapping.new_from_old {
        for w in ws {
            a.extend([w.source_cell as f64, w.volume]);
        }
    }
    close("amr_refine", &a, 1e-15, 1e-17);
    let f2: Vec<f64> = r2
        .mesh
        .cells()
        .iter()
        .map(|c| {
            (-((c.lower.x - 0.4) * (c.lower.x - 0.4) + c.lower.y * c.lower.y)).exp()
                + 0.1 * c.upper.z
        })
        .collect();
    let ind = scalar_jump_error_indicator(&r2.mesh, &f2).unwrap();
    let dm = mark_dorfler(&ind, 0.6).unwrap();
    let co = MarkedHexCoarsening::new(vec![true; r2.mesh.cell_count()], 0, true)
        .apply(&r2.mesh)
        .unwrap();
    let mut v = ind.clone();
    v.extend(dm.iter().map(|m| b(*m)));
    v.extend(dump(&co.mesh));
    v.extend(conservative_adaptive_remap(&r2.mesh, &f2, &co.mesh).unwrap());
    close("amr_coarsen", &v, 1e-13, 1e-16);
    let cfg = AdaptiveAmrConfig {
        maximum_level: 2,
        refine_theta: 0.4,
        coarsen_relative_threshold: 0.3,
        ..Default::default()
    };
    let ad = adapt_scalar_field(&r2.mesh, &f2, &cfg).unwrap();
    let mut c = dump(&ad.mesh);
    c.extend(&ad.field);
    c.extend(&ad.indicator);
    c.extend([ad.refined_cell_count as f64, ad.coarsened_cell_count as f64]);
    for f in ad.mesh.poly_mesh().unwrap().faces() {
        c.extend([
            f.owner as f64,
            f.area.x + 2.0 * f.area.y + 3.0 * f.area.z,
            f.center.x + f.center.y + f.center.z,
        ]);
    }
    close("amr_adapt", &c, 1e-13, 1e-16);
}
