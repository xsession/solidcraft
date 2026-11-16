//! Value-by-value comparison with the C++ upstream FEM solver/kinematics units and the shared
//! material record. Goldens come from `parity/oracle_fem_batch.cpp`.
use cfd_solvers::{core::*, fem::*, physics::*, Complex64};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-fem-batch.json")).unwrap()
}

/// `|rust - cpp| <= rel * max(|cpp|, 1e-6 * max|cpp|) + abs` entry by entry.
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

fn d2(v: &[Displacement2]) -> Vec<f64> {
    v.iter().flat_map(|d| [d.x, d.y]).collect()
}
fn n(x: f64, y: f64) -> Node2 {
    Node2 { x, y }
}
fn lin(max_iterations: usize, relative_tolerance: f64) -> LinearConfig {
    LinearConfig {
        max_iterations,
        relative_tolerance,
        ..Default::default()
    }
}

/// Iterative solves stop at 1e-12 relative residual, so solutions agree to about 1e-9.
const SOLVE: f64 = 1e-8;

#[test]
fn flow_solvers_match_upstream() {
    let mesh = make_rectangle_tri_mesh(6, 5, 1.2, 1.0).unwrap();
    let mut br = Brinkman2D::new(
        mesh.clone(),
        Brinkman2DConfig {
            viscosity: 1e-3,
            effective_viscosity: 2e-3,
            permeability: 1e-4,
            max_iterations: 3000,
            tolerance: 1e-12,
        },
    )
    .unwrap();
    br.solve(|p| Displacement2 {
        x: 1.0 + p.y,
        y: -0.5 * p.x,
    })
    .unwrap();
    close("brinkman", &d2(br.velocity()), SOLVE, 0.0);
    let mut st = PenaltyStokes2D::new(
        mesh.clone(),
        PenaltyStokes2DConfig {
            viscosity: 1.0,
            divergence_penalty: 500.0,
            max_iterations: 5000,
            relative_tolerance: 1e-12,
        },
    )
    .unwrap();
    st.set_dirichlet(|p| Displacement2 {
        x: if p.y >= 1.0 - 1e-12 { 1.0 } else { 0.0 },
        y: 0.0,
    });
    st.solve(Some(&|p: Node2| Displacement2 { x: 0.0, y: -p.x }))
        .unwrap();
    close("stokes", &d2(st.velocity()), SOLVE, 0.0);
    close("stokes_div", &[st.divergence_l2().unwrap()], SOLVE, 0.0);
    let ms = assemble_mixed_stokes_tri3(&mesh, &mesh.triangles[7], 0.8, 1e-3).unwrap();
    let mut m = ms.matrix.to_vec();
    m.push(ms.area);
    close("mixed_stokes", &m, 1e-14, 1e-16);
    close(
        "thermoelastic",
        &assemble_thermoelastic_block(0.3, 2e5, 4.0, 1e6, 0.1, 1.2e-5)
            .unwrap()
            .matrix,
        1e-15,
        0.0,
    );
    let (s1, s2) = (
        supg_parameters(0.01, 2.0, 1e-4, 1e-3).unwrap(),
        supg_parameters(0.02, -0.5, 0.0, 1e-2).unwrap(),
    );
    close(
        "supg",
        &[
            s1.tau_s,
            s1.element_peclet,
            s2.tau_s,
            s2.element_peclet,
            dg_upwind_flux(1.0, 3.0, 2.0, 0.5).unwrap(),
            dg_upwind_flux(1.0, 3.0, -2.0, 0.5).unwrap(),
        ],
        1e-15,
        0.0,
    );

    let mut bar = DynamicBar1D::new(DynamicBar1DConfig {
        elements: 12,
        dt: 2e-6,
        ..Default::default()
    })
    .unwrap();
    bar.initialize(|x| 1e-5 * x, Some(&|x: f64| 0.1 * x * x))
        .unwrap();
    bar.run(15, 2e3).unwrap();
    let mut bv = bar.displacement().to_vec();
    bv.extend_from_slice(bar.velocity());
    bv.extend([
        bar.time(),
        bar.kinetic_energy().unwrap(),
        bar.strain_energy().unwrap(),
    ]);
    close("bar", &bv, 1e-8, 0.0);

    let mut ec = EddyCurrent2D::new(
        mesh.clone(),
        EddyCurrent2DConfig {
            frequency_hz: 5.0,
            max_iterations: 4000,
            gmres_restart: 50,
            relative_tolerance: 1e-11,
        },
    )
    .unwrap();
    ec.solve(
        |p| if p.x > 0.6 { 0.05 } else { 1.0 },
        |p| if p.x > 0.6 { 1.0 } else { 0.0 },
        |p| {
            Complex64::new(
                if p.x < 0.4 { 1.0 } else { 0.0 },
                if p.y < 0.3 { 0.2 } else { 0.0 },
            )
        },
    )
    .unwrap();
    let ev: Vec<f64> = ec
        .vector_potential()
        .iter()
        .flat_map(|a| [a.re, a.im])
        .collect();
    close("eddy", &ev, 1e-7, 0.0);

    let mut da = Darcy2D::new(
        mesh,
        Darcy2DConfig {
            viscosity: 2.0,
            density: 0.1,
            gravity: DarcyVelocity2 { x: 0.0, y: -9.81 },
            linear: lin(4000, 1e-12),
        },
    )
    .unwrap();
    da.set_boundary(
        0,
        ScalarBoundaryCondition2D::Dirichlet(Box::new(|p| 2.0 + 0.5 * p.y)),
    );
    da.set_boundary(1, ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 1.0)));
    da.set_boundary(2, ScalarBoundaryCondition2D::Neumann(Box::new(|_| 0.0)));
    let perm = |p: Node2| 1.0 + p.x;
    da.solve(perm, Some(&|p: Node2| 0.3 * p.y)).unwrap();
    close("darcy_p", da.pressure(), SOLVE, 0.0);
    let dv: Vec<f64> = da
        .element_velocity(perm)
        .unwrap()
        .iter()
        .flat_map(|v| [v.x, v.y])
        .collect();
    close("darcy_v", &dv, 1e-6, 0.0);
}

#[test]
fn axisymmetric_magnetostatic_and_nonlinear_poisson_match_upstream() {
    let mut ring = make_rectangle_tri_mesh(5, 4, 1.0, 2.0).unwrap();
    for p in &mut ring.nodes {
        p.x += 0.5;
    }
    let mut ax = AxisymmetricElasticity::new(
        ring,
        AxisymmetricElasticityConfig {
            youngs_modulus: 2e5,
            poisson_ratio: 0.28,
            max_iterations: 5000,
            relative_tolerance: 1e-12,
        },
    )
    .unwrap();
    ax.set_dirichlet(|p| p.y <= 1e-12, |_| Displacement2::default(), false, true);
    ax.set_dirichlet(
        |p| p.x <= 0.5 + 1e-12,
        |_| Displacement2 { x: 1e-3, y: 0.0 },
        true,
        false,
    );
    ax.solve(Some(&|p: Node2| Displacement2 {
        x: 10.0 * p.x,
        y: -5.0,
    }))
    .unwrap();
    close("axisym", &d2(ax.displacement()), SOLVE, 0.0);

    let mesh = make_rectangle_tri_mesh(6, 5, 1.2, 1.0).unwrap();
    let mut mg = Magnetostatics2D::new(mesh.clone(), lin(4000, 1e-12)).unwrap();
    for patch in [0, 1, 3] {
        mg.set_boundary(
            patch,
            ScalarBoundaryCondition2D::Dirichlet(Box::new(|_| 0.0)),
        );
    }
    let coil = RectangularStrandedCoil2D {
        xmin: 0.3,
        xmax: 0.7,
        ymin: 0.2,
        ymax: 0.6,
        turns: 50.0,
        current_a: 2.0,
    };
    let bh = PiecewiseLinearBHCurve::new(vec![0.0, 0.5, 1.2, 1.6], vec![0.0, 100.0, 600.0, 4000.0])
        .unwrap();
    let mu0 = 4e-7 * std::f64::consts::PI;
    let nl = mg
        .solve_nonlinear(
            |b| bh.reluctivity(b).unwrap().min(1.0 / mu0),
            |p| coil.current_density_a_per_m2(p).unwrap() * 1e3,
            30,
            1e-8,
            0.6,
        )
        .unwrap();
    close("mag_a", mg.vector_potential(), 1e-7, 0.0);
    let b = mg.element_flux_density().unwrap();
    close(
        "mag_b",
        &b.iter().flat_map(|v| [v.x, v.y]).collect::<Vec<_>>(),
        1e-7,
        0.0,
    );
    let ft = maxwell_stress_boundary_force_2d(&mesh, &b, 2, mu0, n(0.6, 0.5)).unwrap();
    let misc = [
        nl.iterations as f64,
        nl.maximum_relative_reluctivity_change,
        f64::from(u8::from(nl.converged)),
        coil_flux_linkage_wb_turn(&mesh, mg.vector_potential(), &coil).unwrap(),
        coil_inductance_h(&mesh, mg.vector_potential(), &coil).unwrap(),
        ft.force_x_n_per_m,
        ft.force_y_n_per_m,
        ft.torque_z_n,
        bh.field_strength(0.8).unwrap(),
        bh.field_strength(2.0).unwrap(),
        bh.reluctivity(0.0).unwrap(),
        bh.reluctivity(1.4).unwrap(),
    ];
    // The final Picard change is at the 1e-9 level, so it is compared only loosely.
    close("mag_misc", &misc, 1e-6, 1e-8);

    let mut np = NonlinearPoisson2D::new(
        mesh,
        NonlinearPoisson2DConfig {
            diffusivity: 0.7,
            cubic_coefficient: 3.0,
            ..Default::default()
        },
    )
    .unwrap();
    np.solve(|p| 10.0 + 5.0 * p.x * p.y, |p| 0.2 * p.x - 0.1 * p.y)
        .unwrap();
    let mut s = np.solution().to_vec();
    let r = np.nonlinear_result().unwrap();
    s.extend([r.iterations as f64]);
    close(
        "nonlinear_poisson",
        &[s.as_slice(), &[r.residual_rms]].concat(),
        1e-9,
        1e-14,
    );
}

#[test]
fn kinematics_phase_change_and_material_match_upstream() {
    let mesh = make_rectangle_tri_mesh(6, 5, 1.2, 1.0).unwrap();
    let tr = green_lagrange_truss(2.0, 2.3, 1e-4, 2e11).unwrap();
    let tl = total_lagrangian_triangle(
        &mesh,
        &mesh.triangles[3],
        &[n(0.2, 0.0), n(0.45, 0.02), n(0.21, 0.23)],
        1e7,
        0.3,
    )
    .unwrap();
    let mut k = vec![
        tr.stretch,
        tr.green_lagrange_strain,
        tr.second_piola_stress,
        tr.axial_force,
        tr.tangent,
        tl.reference_area,
        tl.determinant,
    ];
    k.extend(tl.deformation_gradient);
    k.extend(tl.green_lagrange_strain);
    k.extend(tl.internal_force);
    close("kinematics", &k, 1e-13, 1e-16);
    let pcm = PhaseChangeMaterial {
        solidus: 273.0,
        liquidus: 278.0,
        heat_capacity_solid: 2000.0,
        heat_capacity_liquid: 4000.0,
        latent_heat: 3.3e5,
    };
    close(
        "phase_change",
        &[
            pcm.liquid_fraction(275.5).unwrap(),
            pcm.specific_enthalpy(270.0).unwrap(),
            pcm.specific_enthalpy(276.0).unwrap(),
            pcm.specific_enthalpy(290.0).unwrap(),
            pcm.effective_heat_capacity(276.0).unwrap(),
            pcm.temperature_from_enthalpy(pcm.specific_enthalpy(276.0).unwrap())
                .unwrap(),
            pcm.temperature_from_enthalpy(-50.0).unwrap(),
            porous_effective_conductivity(0.3, 2.0, 0.6).unwrap(),
        ],
        1e-14,
        0.0,
    );
    let good = Material {
        name: "steel".into(),
        density: 7800.0,
        heat_capacity: 500.0,
        thermal_conductivity: 40.0,
        electrical_conductivity: 1e6,
        relative_permittivity: [2.0, 3.0, 4.0],
        relative_permeability: [1.0; 3],
        refractive_index: 1.5,
    };
    let base = Material {
        name: "x".into(),
        density: 1.0,
        heat_capacity: 1.0,
        thermal_conductivity: 1.0,
        electrical_conductivity: 1.0,
        ..Default::default()
    };
    let bad = [
        Material {
            name: String::new(),
            ..base.clone()
        },
        Material {
            relative_permittivity: [1.0, 0.0, 1.0],
            ..base.clone()
        },
        Material {
            refractive_index: 0.0,
            ..base
        },
    ];
    let mut v = vec![
        fem_conductivity(&good).unwrap(),
        optical_refractive_index(&good).unwrap(),
    ];
    v.extend(
        bad.iter()
            .map(|m| f64::from(u8::from(validate_material(m).is_ok()))),
    );
    close("material", &v, 0.0, 0.0);
}
