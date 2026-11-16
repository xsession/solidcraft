//! Analytic checks for hexahedral AMR, mixture thermo, mesh motion, reordering and fields.
use cfd_solvers::{core::Vec3, fvm::*};

#[test]
fn refinement_conserves_volume_balance_and_integrals() {
    let mesh = AdaptiveHexMesh::cartesian(2, 2, 2, 1.0, 1.0, 1.0).unwrap();
    let mut marks = vec![false; 8];
    marks[0] = true;
    let r1 = MarkedHexRefinement::new(marks, 4, true)
        .unwrap()
        .apply(&mesh)
        .unwrap();
    let mut marks = vec![false; r1.mesh.cell_count()];
    marks[0] = true;
    let r2 = MarkedHexRefinement::new(marks, 4, true)
        .unwrap()
        .apply(&r1.mesh)
        .unwrap();
    assert!(satisfies_two_to_one(&r2.mesh));
    assert!((r2.mesh.total_volume() - 1.0).abs() < 1e-14);
    let field: Vec<f64> = r1
        .mesh
        .cells()
        .iter()
        .enumerate()
        .map(|(i, _)| i as f64 * 0.3 - 1.0)
        .collect();
    let integral = |m: &AdaptiveHexMesh, f: &[f64]| {
        m.cells()
            .iter()
            .zip(f)
            .map(|(c, v)| {
                v * (c.upper.x - c.lower.x) * (c.upper.y - c.lower.y) * (c.upper.z - c.lower.z)
            })
            .sum::<f64>()
    };
    let fine = conservative_adaptive_remap(&r1.mesh, &field, &r2.mesh).unwrap();
    assert!((integral(&r1.mesh, &field) - integral(&r2.mesh, &fine)).abs() < 1e-14);
    // Coarsening everything returns to the root level with the same integral.
    let back = MarkedHexCoarsening::new(vec![true; r2.mesh.cell_count()], 0, true)
        .apply(&r2.mesh)
        .unwrap();
    let back2 = MarkedHexCoarsening::new(vec![true; back.mesh.cell_count()], 0, true)
        .apply(&back.mesh)
        .unwrap();
    assert_eq!(back2.mesh.cell_count(), 8);
    // The poly mesh closes: outward areas of every cell sum to zero.
    let poly = r2.mesh.poly_mesh().unwrap();
    let mut sum = vec![Vec3::new(0.0, 0.0, 0.0); poly.cells().len()];
    for f in poly.faces() {
        let s = &mut sum[f.owner];
        *s = Vec3::new(s.x + f.area.x, s.y + f.area.y, s.z + f.area.z);
        if let Some(n) = f.neighbour {
            let s = &mut sum[n];
            *s = Vec3::new(s.x - f.area.x, s.y - f.area.y, s.z - f.area.z);
        }
    }
    assert!(sum
        .iter()
        .all(|s| s.x.abs() + s.y.abs() + s.z.abs() < 1e-14));
    assert!(MarkedHexRefinement::new(vec![], 21, true).is_err());
    assert!(MarkedHexRefinement::new(vec![true], 4, true)
        .unwrap()
        .apply(&mesh)
        .is_err());
    let bad = AdaptiveHexCell {
        lower: Vec3::new(0.0, 0.0, 0.0),
        upper: Vec3::new(2.0, 1.0, 1.0),
        ..Default::default()
    };
    assert!(AdaptiveHexMesh::new(
        vec![bad],
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 1.0, 1.0)
    )
    .is_err());
}

#[test]
fn adapt_step_refines_a_front_and_stays_balanced() {
    let mesh = AdaptiveHexMesh::cartesian(4, 4, 1, 1.0, 1.0, 0.25).unwrap();
    let field: Vec<f64> = mesh
        .cells()
        .iter()
        .map(|c| if c.lower.x < 0.5 { 1.0 } else { 0.0 })
        .collect();
    let out = adapt_scalar_field(&mesh, &field, &AdaptiveAmrConfig::default()).unwrap();
    assert!(out.refined_cell_count > 0 && satisfies_two_to_one(&out.mesh));
    assert!(mark_dorfler(&[0.0, 0.0], 0.5).unwrap().iter().all(|m| !m));
    assert!(mark_dorfler(&[1.0], 0.0).is_err());
    assert!(scalar_jump_error_indicator(&mesh, &field)
        .unwrap()
        .iter()
        .all(|v| v.is_finite()));
}

#[test]
fn mixture_thermo_motion_reordering_and_fields() {
    let sp = [IdealGasSpeciesThermo {
        molar_mass: 0.029,
        cp_mass: 1005.0,
        formation_enthalpy_mass: 0.0,
    }];
    let s = ideal_gas_mixture_state(300.0, 101325.0, &[1.0], &sp).unwrap();
    assert!((s.density - 101325.0 * 0.029 / (8.31446261815324 * 300.0)).abs() < 1e-12);
    assert!(
        (temperature_from_mixture_enthalpy(s.specific_enthalpy, &[1.0], &sp).unwrap() - 300.0)
            .abs()
            < 1e-12
    );
    assert!(ideal_gas_mixture_state(300.0, 1e5, &[0.5], &sp).is_err());
    let mesh = make_cartesian_hexa_mesh(3, 3, 1, 1.0, 1.0, 1.0).unwrap();
    // A pure translation sweeps flux U.S through each face.
    let moved = rigidly_moved_mesh(
        &mesh,
        Vec3::new(0.2, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        0.0,
    )
    .unwrap();
    let phi = ale_mesh_flux(&mesh, &moved, 0.1).unwrap();
    for (p, f) in phi.iter().zip(mesh.faces()) {
        assert!((p - 2.0 * f.area.x).abs() < 1e-12);
    }
    // A rotation preserves face-area magnitudes.
    let rot = rigidly_moved_mesh(
        &mesh,
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 1.0, 0.0),
        1.0,
    )
    .unwrap();
    for (a, b) in mesh.faces().iter().zip(rot.faces()) {
        assert!((a.area.norm() - b.area.norm()).abs() < 1e-14);
    }
    let order = morton_cell_order(&mesh);
    let mut sorted = order.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, (0..9).collect::<Vec<_>>());
    let re = reorder_cells(&mesh, &order).unwrap();
    assert_eq!(re.cells().len(), 9);
    assert!(reorder_cells(&mesh, &[0; 9]).is_err());
    assert_eq!(
        reorder_cell_field(&[10, 20, 30], &[2, 0, 1]).unwrap(),
        vec![30, 10, 20]
    );
    let f = vec![1.0; 9];
    let remapped = conservative_nearest_remap(
        &mesh,
        &f,
        &make_cartesian_hexa_mesh(5, 5, 1, 1.0, 1.0, 1.0).unwrap(),
    )
    .unwrap();
    assert!(remapped.iter().all(|v| (v - 1.0).abs() < 1e-14));
    let mut h = FieldHistory::new(2).unwrap();
    h.push(vec![1.0]).unwrap();
    h.push(vec![2.0]).unwrap();
    h.push(vec![3.0]).unwrap();
    assert_eq!(
        (
            h.current().unwrap()[0],
            h.old(1).unwrap()[0],
            h.stored_levels()
        ),
        (3.0, 2.0, 2)
    );
    assert!(h.old(2).is_err() && h.push(vec![1.0, 2.0]).is_err());
    let velocity = Dimensions::LENGTH - Dimensions::TIME;
    assert_eq!(velocity.exponent[..3], [1, 0, -1]);
    let field = DimensionedField::new("U", velocity, FieldLocation::Cell, 4, 0.0);
    assert_eq!((field.size(), field.name()), (4, "U"));
    let bdf = Bdf2Coefficients::default();
    assert_eq!(bdf.current + bdf.old + bdf.older, 0.0);
}
