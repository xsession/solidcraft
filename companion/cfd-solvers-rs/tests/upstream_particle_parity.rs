//! Value-by-value comparison with the C++ upstream particle/{electromagnetic, cross_sections,
//! monte_carlo, plasma, wakefield, transport}. Goldens come from `parity/oracle_particle.cpp`.
#![allow(clippy::needless_range_loop)]
use cfd_solvers::particle::*;
use serde_json::Value;

fn g() -> &'static Value {
    static G: std::sync::OnceLock<Value> = std::sync::OnceLock::new();
    G.get_or_init(|| serde_json::from_str(include_str!("../parity/goldens-particle.json")).unwrap())
}

fn num(v: &Value) -> f64 {
    match v.as_str() {
        Some("inf") => f64::INFINITY,
        Some("-inf") => f64::NEG_INFINITY,
        Some(_) => f64::NAN,
        None => v.as_f64().unwrap(),
    }
}

/// Relative tolerance with an absolute floor scaled by the largest magnitude in the vector.
#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64) {
    let cpp: Vec<f64> = g()[name].as_array().unwrap().iter().map(num).collect();
    assert_eq!(rust.len(), cpp.len(), "{name} length");
    let scale = cpp
        .iter()
        .filter(|v| v.is_finite())
        .fold(0.0f64, |m, v| m.max(v.abs()));
    for (i, (a, b)) in rust.iter().zip(&cpp).enumerate() {
        let ok = a == b
            || (a.is_nan() && b.is_nan())
            || (a - b).abs() <= rel * a.abs().max(b.abs()) + rel * 1e-3 * scale;
        assert!(ok, "{name}[{i}]: rust {a} cpp {b}");
    }
}

type G = Vec<f64>;
fn v3(g: &mut G, v: Vec3) {
    g.extend([v.x, v.y, v.z]);
}
fn cp(g: &mut G, p: &ChargedParticle) {
    v3(g, p.position_m);
    v3(g, p.velocity_m_per_s);
    g.push(p.weight);
}
fn p2(g: &mut G, p: &PicParticle2D) {
    g.extend([p.position_m.x, p.position_m.y]);
    v3(g, p.velocity_m_per_s);
    g.push(p.weight);
}
fn p3(g: &mut G, p: &PicParticle3D) {
    v3(g, p.position_m);
    v3(g, p.velocity_m_per_s);
    g.push(p.weight);
}
fn cs(g: &mut G, s: &CollisionStatistics) {
    g.extend([
        s.elastic_events as f64,
        s.ionization_events as f64,
        s.energy_loss_j,
    ]);
}
fn fld(g: &mut G, f: &ElectromagneticField) {
    v3(g, f.electric_v_per_m);
    v3(g, f.magnetic_t);
}
fn b(x: bool) -> f64 {
    if x {
        1.0
    } else {
        0.0
    }
}
fn v2(x: f64, y: f64) -> Vec2 {
    Vec2 { x, y }
}
const Q: f64 = -1.602176634e-19;
const ME: f64 = 9.1093837015e-31;
fn cpart(pos: Vec3, vel: Vec3, charge: f64, mass: f64, weight: f64) -> ChargedParticle {
    ChargedParticle {
        position_m: pos,
        velocity_m_per_s: vel,
        charge_c: charge,
        mass_kg: mass,
        weight,
    }
}

#[test]
fn pushers_collisions_walls() {
    let mut g = G::new();
    let f = ElectromagneticField {
        electric_v_per_m: Vec3::new(1.0e3, -2.0e3, 5.0e2),
        magnetic_t: Vec3::new(0.01, 0.02, -0.03),
    };
    let mut a = cpart(
        Vec3::new(0.1, 0.2, 0.3),
        Vec3::new(1.0e5, -2.0e5, 3.0e5),
        Q,
        ME,
        2.0,
    );
    let mut bb = a;
    bb.velocity_m_per_s = Vec3::new(1.0e8, 5.0e7, -2.0e7);
    let mut c = bb;
    for _ in 0..7 {
        boris_push(&mut a, &f, 1.0e-11).unwrap();
        relativistic_boris_push(&mut bb, &f, 1.0e-11).unwrap();
        vay_push(&mut c, &f, 1.0e-11).unwrap();
    }
    cp(&mut g, &a);
    cp(&mut g, &bb);
    cp(&mut g, &c);
    g.extend([
        lorentz_gamma(&a).unwrap(),
        lorentz_gamma(&bb).unwrap(),
        kinetic_energy_j(&a).unwrap(),
        kinetic_energy_j(&c).unwrap(),
    ]);
    let mut t = cpart(Vec3::default(), Vec3::new(1.0e4, 0.0, 0.0), 1.0, 2.0, 1.0);
    let sampler = |p: Vec3, time: f64| {
        Ok(ElectromagneticField {
            electric_v_per_m: Vec3::new(p.x.sin() + time, 0.5, 0.0),
            magnetic_t: Vec3::new(0.0, 0.0, 1.0 + p.y),
        })
    };
    track_particle(&mut t, &sampler, 0.25, 1.0e-3, 9).unwrap();
    cp(&mut g, &t);
    close("pushers", &g, 1e-12);

    let mut g = G::new();
    let mut ps: Vec<ChargedParticle> = (0..6)
        .map(|i| {
            cpart(
                Vec3::new(0.1 * i as f64, 0.0, 0.0),
                Vec3::new(1.0e6 * (1 + i) as f64, 2.0e5 * i as f64, -3.0e5),
                Q,
                ME,
                1.0 + i as f64,
            )
        })
        .collect();
    let m = NeutralCollisionModel {
        neutral_density_per_m3: 1.0e22,
        elastic_cross_section_m2: 1.0e-19,
        ionization_cross_section_m2: 5.0e-20,
        ionization_energy_ev: 15.0,
        random_seed: 42,
    };
    cs(
        &mut g,
        &apply_monte_carlo_collisions(&mut ps, 1.0e-9, &m).unwrap(),
    );
    cs(
        &mut g,
        &apply_monte_carlo_collisions(&mut ps, 1.0e-9, &m).unwrap(),
    );
    g.push(ps.len() as f64);
    for p in &ps {
        cp(&mut g, p);
    }
    let s = SecondaryEmissionModel::default();
    for e in [0.0, 10.0, 20.0, 50.0, 300.0, 1000.0, 5000.0] {
        g.push(secondary_electron_yield(e, &s).unwrap());
    }
    let bx = AxisAlignedParticleBox {
        minimum_m: Vec3::default(),
        maximum_m: Vec3::new(1.0, 2.0, 3.0),
    };
    for mode in [ParticleWallMode::Absorb, ParticleWallMode::SpecularReflect] {
        for pos in [
            Vec3::new(1.2, 0.5, 0.5),
            Vec3::new(-0.1, 2.3, 1.0),
            Vec3::new(0.5, 1.0, -0.4),
            Vec3::new(0.5, 0.5, 0.5),
        ] {
            let mut p = cpart(pos, Vec3::new(2.0e6, -1.0e6, 5.0e5), Q, ME, 3.0);
            let r = apply_particle_box_boundary(&mut p, &bx, mode, Some(&s)).unwrap();
            cp(&mut g, &p);
            g.extend([b(r.alive), b(r.impacted)]);
            v3(&mut g, r.impact_normal);
            g.push(r.secondary_macro_weight);
        }
    }
    close("collisions_walls", &g, 1e-12);
}

#[test]
fn pic_1d() {
    let mut g = G::new();
    let rho: Vec<f64> = (0..16)
        .map(|i| (2.0 * std::f64::consts::PI * i as f64 / 16.0).sin() + 0.1 * i as f64)
        .collect();
    g.extend(periodic_electric_field_from_charge_density(&rho, 2.0, 2.5, true).unwrap());
    let rho0: Vec<f64> = (0..16)
        .map(|i| if i % 2 == 1 { 1.0 } else { -1.0 } * (1 + i / 2) as f64)
        .collect();
    g.extend(periodic_electric_field_from_charge_density(&rho0, 2.0, 2.5, false).unwrap());
    let c = ElectrostaticPic1DConfig {
        grid_points: 32,
        length_m: 1.0,
        dt_s: 0.01,
        ..Default::default()
    };
    let mut es = ElectrostaticPic1D::new(c).unwrap();
    let pp = |x, v, q, m, w| PicParticle1D {
        position_m: x,
        velocity_m_per_s: v,
        charge_c: q,
        mass_kg: m,
        weight: w,
    };
    es.set_particles(vec![
        pp(0.1, 0.5, 1.0, 1.0, 1.0),
        pp(0.45, -0.2, -1.0, 2.0, 1.5),
        pp(0.8, 0.0, 1.0, 1.0, 0.5),
        pp(1.3, 0.1, 0.5, 1.0, 1.0),
    ])
    .unwrap();
    es.step(5).unwrap();
    g.extend_from_slice(es.charge_density());
    g.extend_from_slice(es.electric_field());
    for p in es.particles() {
        g.extend([p.position_m, p.velocity_m_per_s]);
    }
    let mut o = Vec::new();
    let mut n = Vec::new();
    for i in 0..4 {
        o.push(cpart(
            Vec3::new(0.1 + 0.2 * i as f64, 0.0, 0.0),
            Vec3::default(),
            if i % 2 == 1 { -1.0 } else { 1.0 },
            1.0,
            1.0 + i as f64,
        ));
        let mut q = o[i];
        q.position_m.x += 0.03 * (i + 1) as f64;
        n.push(q);
    }
    let cc = deposit_charge_conserving_current_1d(&o, &n, 16, 1.0, 0.1).unwrap();
    g.extend(&cc.charge_density_old_c_per_m3);
    g.extend(&cc.charge_density_new_c_per_m3);
    g.extend(&cc.current_x_a_per_m2);
    g.push(cc.continuity_linf_residual);
    close("pic1d_es", &g, 1e-11);

    let mut g = G::new();
    let c = ElectromagneticPic1DConfig {
        grid_points: 24,
        dt_s: 1.0e-10,
        ..Default::default()
    };
    let mut em = ElectromagneticPic1D::new(c).unwrap();
    em.initialize_right_traveling_mode(100.0, 2).unwrap();
    let ps = (0..5)
        .map(|i| {
            cpart(
                Vec3::new(0.13 + 0.17 * i as f64, 0.0, 0.0),
                Vec3::new(1.0e6 * (i - 2) as f64, 2.0e5, -1.0e5 * i as f64),
                Q * 1e10,
                ME * 1e10,
                1.0,
            )
        })
        .collect();
    em.set_particles(ps).unwrap();
    let m = NeutralCollisionModel {
        neutral_density_per_m3: 1.0e20,
        elastic_cross_section_m2: 1.0e-19,
        ionization_cross_section_m2: 1.0e-20,
        ionization_energy_ev: 15.0,
        random_seed: 7,
    };
    em.step(4, Some(&m)).unwrap();
    em.step(2, None).unwrap();
    let d = em.diagnostics().unwrap();
    g.extend([
        d.field_energy_j,
        d.particle_kinetic_energy_j,
        d.total_energy_j,
        d.charge_continuity_linf_residual,
    ]);
    cs(&mut g, &d.collision_statistics);
    g.push(d.particle_count as f64);
    for v in [
        em.charge_density(),
        em.current_x(),
        em.current_y(),
        em.electric_x(),
        em.electric_y(),
        em.electric_z(),
        em.magnetic_y(),
        em.magnetic_z(),
    ] {
        g.extend_from_slice(v);
    }
    for p in em.particles() {
        cp(&mut g, p);
    }
    g.push(em.time_s());
    fld(&mut g, &em.gather_field(0.37));
    close("pic1d_em", &g, 1e-10);
}

fn ps2() -> Vec<PicParticle2D> {
    (0..6)
        .map(|i| PicParticle2D {
            position_m: v2(0.07 + 0.15 * i as f64, 0.9 - 0.13 * i as f64),
            velocity_m_per_s: Vec3::new(2.0e5 * (i - 3) as f64, 1.0e5 * i as f64, 3.0e4),
            charge_c: if i % 2 == 1 { Q } else { -Q } * 1e9,
            mass_kg: ME * 1e9,
            weight: 1.0 + 0.5 * i as f64,
        })
        .collect()
}

#[test]
fn pic_2d() {
    let ps2 = ps2();
    let mut g = G::new();
    g.extend(deposit_cic_charge_density_2d(&ps2, 8, 6, 1.0, 1.0).unwrap());
    let rho: Vec<f64> = (0..48)
        .map(|i| (0.3 * i as f64).cos() + 0.01 * i as f64)
        .collect();
    let f =
        periodic_electric_field_from_charge_density_2d(&rho, 8, 6, 1.0, 1.5, 2.0, true).unwrap();
    g.extend(&f.electric_x_v_per_m);
    g.extend(&f.electric_y_v_per_m);
    let mut n2 = ps2.clone();
    for p in &mut n2 {
        p.position_m.x += 0.02;
        p.position_m.y -= 0.015;
    }
    let cc = deposit_charge_conserving_current_2d(&ps2, &n2, 8, 6, 1.0, 1.0, 1.0e-9).unwrap();
    g.extend(&cc.charge_density_new_c_per_m3);
    g.extend(&cc.current_x_a_per_m2);
    g.extend(&cc.current_y_a_per_m2);
    g.push(cc.continuity_linf_residual);
    let mut es = ElectrostaticPic2D::new(ElectrostaticPic2DConfig {
        nx: 8,
        ny: 8,
        dt_s: 1.0e-9,
        ..Default::default()
    })
    .unwrap();
    es.set_particles(ps2.clone()).unwrap();
    es.step(3).unwrap();
    g.extend_from_slice(es.charge_density());
    g.extend_from_slice(es.electric_x());
    for p in es.particles() {
        p2(&mut g, p);
    }
    let e = es.gather_electric_field(v2(0.33, 0.71));
    g.extend([e.x, e.y]);
    let mut em = ElectromagneticPic2D::new(ElectromagneticPic2DConfig {
        nx: 10,
        ny: 8,
        dt_s: 1.0e-11,
        ..Default::default()
    })
    .unwrap();
    em.initialize_tm_z_mode(50.0, 1, 1).unwrap();
    em.set_particles(ps2.clone()).unwrap();
    let m = NeutralCollisionModel {
        neutral_density_per_m3: 1.0e21,
        elastic_cross_section_m2: 1.0e-19,
        ionization_cross_section_m2: 1.0e-20,
        ionization_energy_ev: 15.0,
        random_seed: 11,
    };
    em.step(3, Some(&m)).unwrap();
    em.step(1, None).unwrap();
    let d = em.diagnostics().unwrap();
    g.extend([
        d.field_energy_j,
        d.particle_kinetic_energy_j,
        d.charge_continuity_linf_residual,
    ]);
    cs(&mut g, &d.collision_statistics);
    g.push(d.particle_count as f64);
    for v in [
        em.electric_z(),
        em.magnetic_x(),
        em.magnetic_y(),
        em.current_x(),
        em.current_z(),
    ] {
        g.extend_from_slice(v);
    }
    for p in em.particles() {
        p2(&mut g, p);
    }
    close("pic2d", &g, 1e-10);

    let mut g = G::new();
    let s = SecondaryEmissionModel::default();
    let bx = PicBoundaryBox2D {
        minimum_m: v2(0.0, 0.0),
        maximum_m: v2(1.0, 1.0),
    };
    for mode in [ParticleWallMode::Absorb, ParticleWallMode::SpecularReflect] {
        for pos in [v2(1.1, 0.5), v2(-0.2, -0.3), v2(0.5, 0.5)] {
            let mut p = PicParticle2D {
                position_m: pos,
                velocity_m_per_s: Vec3::new(1.0e6, -2.0e6, 0.0),
                charge_c: Q,
                mass_kg: ME,
                weight: 2.0,
            };
            let r = apply_particle_box_boundary_2d(&mut p, &bx, mode, Some(&s)).unwrap();
            p2(&mut g, &p);
            g.extend([
                b(r.alive),
                b(r.impacted),
                r.impact_normal.x,
                r.impact_normal.y,
                r.secondary_macro_weight,
            ]);
        }
    }
    for mode in [
        GridBoundaryMode2D::ElectricWall,
        GridBoundaryMode2D::AbsorbingSponge,
    ] {
        let (mut ex, mut ey) = (vec![1.0; 30], vec![-2.0; 30]);
        let bc = GridBoundary2DConfig {
            mode,
            sponge_cells: 2,
            sponge_strength: 2.5,
            ..Default::default()
        };
        apply_electrostatic_field_boundary_2d(&mut ex, &mut ey, 6, 5, &bc).unwrap();
        g.extend(ex);
        g.extend(ey);
    }
    for variant in 0..2 {
        let mut c = StaggeredElectromagneticPic2DConfig {
            nx: 10,
            ny: 9,
            dt_s: 5.0e-11,
            ..Default::default()
        };
        if variant == 1 {
            c.field_boundary.mode = GridBoundaryMode2D::AbsorbingSponge;
            c.periodic_particles = false;
            c.particle_boundary = ParticleWallMode::Absorb;
            c.enable_secondary_yield_report = true;
        }
        let mut st = StaggeredElectromagneticPic2D::new(c).unwrap();
        st.initialize_tm_z_mode(30.0, 1, 2).unwrap();
        let mut ps = ps2.clone();
        for p in &mut ps {
            p.velocity_m_per_s.x *= 10.0 * (variant + 1) as f64;
        }
        ps[0].position_m.x = 0.9995;
        ps[0].velocity_m_per_s.x = 6.0e7;
        ps[1].position_m.y = 0.0004;
        ps[1].velocity_m_per_s.y = -5.0e7;
        st.set_particles(ps).unwrap();
        st.step(3, None).unwrap();
        let d = st.diagnostics().unwrap();
        g.extend([
            d.field_energy_j,
            d.particle_kinetic_energy_j,
            d.charge_continuity_linf_residual,
            d.particle_count as f64,
            d.absorbed_particles as f64,
            d.secondary_macro_weight,
        ]);
        for v in [
            st.electric_z(),
            st.magnetic_x(),
            st.magnetic_y(),
            st.current_x(),
            st.charge_density(),
        ] {
            g.extend_from_slice(v);
        }
        for p in st.particles() {
            p2(&mut g, p);
        }
        fld(&mut g, &st.gather_field(v2(0.41, 0.27)));
    }
    close("pic2d_staggered", &g, 1e-10);
}

fn ps3() -> Vec<PicParticle3D> {
    (0..7)
        .map(|i| PicParticle3D {
            position_m: Vec3::new(
                0.05 + 0.13 * i as f64,
                0.93 - 0.11 * i as f64,
                0.2 + 0.09 * i as f64,
            ),
            velocity_m_per_s: Vec3::new(
                1.5e5 * (i - 3) as f64,
                1.0e5 * (i % 3) as f64,
                -2.0e4 * i as f64,
            ),
            charge_c: if i % 2 == 1 { Q } else { -Q } * 1e9,
            mass_kg: ME * 1e9,
            weight: 1.0 + 0.25 * i as f64,
        })
        .collect()
}

#[test]
fn pic_3d() {
    let ps3 = ps3();
    let mut g = G::new();
    g.extend(deposit_cic_charge_density_3d(&ps3, 5, 4, 6, 1.0, 1.2, 0.8).unwrap());
    let rho: Vec<f64> = (0..120)
        .map(|i| (0.21 * i as f64).sin() + 0.02 * i as f64)
        .collect();
    let f = periodic_electric_field_from_charge_density_3d(&rho, 5, 4, 6, 1.0, 1.2, 0.8, 1.5, true)
        .unwrap();
    g.extend(&f.electric_x_v_per_m);
    g.extend(&f.electric_y_v_per_m);
    g.extend(&f.electric_z_v_per_m);
    let mut n3 = ps3.clone();
    for p in &mut n3 {
        p.position_m.x += 0.03;
        p.position_m.y -= 0.02;
        p.position_m.z += 0.01;
    }
    for local in [false, true] {
        let cc = if local {
            deposit_charge_conserving_current_3d_local(&ps3, &n3, 5, 4, 6, 1.0, 1.2, 0.8, 1.0e-9)
        } else {
            deposit_charge_conserving_current_3d(&ps3, &n3, 5, 4, 6, 1.0, 1.2, 0.8, 1.0e-9)
        }
        .unwrap();
        for v in [
            &cc.charge_density_old_c_per_m3,
            &cc.charge_density_new_c_per_m3,
            &cc.current_x_a_per_m2,
            &cc.current_y_a_per_m2,
            &cc.current_z_a_per_m2,
        ] {
            g.extend(v);
        }
        g.push(cc.continuity_linf_residual);
    }
    let mut es = ElectrostaticPic3D::new(ElectrostaticPic3DConfig {
        nx: 6,
        ny: 5,
        nz: 4,
        dt_s: 1.0e-9,
        ..Default::default()
    })
    .unwrap();
    es.set_particles(ps3.clone()).unwrap();
    es.step(3).unwrap();
    g.extend_from_slice(es.charge_density());
    g.extend_from_slice(es.electric_x());
    g.extend_from_slice(es.electric_z());
    for p in es.particles() {
        p3(&mut g, p);
    }
    v3(
        &mut g,
        es.gather_electric_field(Vec3::new(0.33, 0.71, 0.52)),
    );
    close("pic3d", &g, 1e-10);

    let mut g = G::new();
    for variant in 0..2 {
        let c = ElectromagneticPic3DConfig {
            nx: 6,
            ny: 5,
            nz: 4,
            dt_s: 1.0e-11,
            solve_longitudinal_poisson: variant == 1,
            ..Default::default()
        };
        let mut em = ElectromagneticPic3D::new(c).unwrap();
        em.initialize_z_polarized_mode(40.0, 1, 1, 0).unwrap();
        em.set_particles(ps3.clone()).unwrap();
        let m = NeutralCollisionModel {
            neutral_density_per_m3: 1.0e21,
            elastic_cross_section_m2: 1.0e-19,
            ionization_cross_section_m2: 1.0e-20,
            ionization_energy_ev: 15.0,
            random_seed: 13,
        };
        em.step(2, Some(&m)).unwrap();
        em.step(1, None).unwrap();
        let d = em.diagnostics().unwrap();
        g.extend([
            d.field_energy_j,
            d.particle_kinetic_energy_j,
            d.charge_continuity_linf_residual,
        ]);
        cs(&mut g, &d.collision_statistics);
        g.push(d.particle_count as f64);
        for v in [
            em.electric_x(),
            em.electric_z(),
            em.magnetic_x(),
            em.magnetic_y(),
            em.current_z(),
        ] {
            g.extend_from_slice(v);
        }
        for p in em.particles() {
            p3(&mut g, p);
        }
        fld(&mut g, &em.gather_field(Vec3::new(0.41, 0.27, 0.66)));
    }
    close("pic3d_em", &g, 1e-10);

    let mut g = G::new();
    for mode in [
        GridBoundaryMode3D::ElectricWall,
        GridBoundaryMode3D::AbsorbingSponge,
    ] {
        let (mut ex, mut ey, mut ez, mut bx, mut by, mut bz) = (
            vec![1.0; 60],
            vec![-2.0; 60],
            vec![3.0; 60],
            vec![0.5; 60],
            vec![0.25; 60],
            vec![-1.0; 60],
        );
        let bc = GridBoundary3DConfig {
            mode,
            sponge_cells: 2,
            sponge_strength: 2.0,
            sponge_polynomial_order: 1.5,
        };
        apply_electromagnetic_field_boundary_3d(
            &mut ex, &mut ey, &mut ez, &mut bx, &mut by, &mut bz, 5, 4, 3, &bc,
        )
        .unwrap();
        for v in [ex, ez, bx, bz] {
            g.extend(v);
        }
    }
    for variant in 0..3 {
        let mut c = StaggeredElectromagneticPic3DConfig {
            nx: 6,
            ny: 5,
            nz: 5,
            dt_s: 5.0e-11,
            ..Default::default()
        };
        if variant == 1 {
            c.field_boundary.mode = GridBoundaryMode3D::ElectricWall;
            c.periodic_particles = false;
            c.particle_boundary = ParticleWallMode::SpecularReflect;
            c.current_deposition = CurrentDeposition3DMode::LocalFiniteVolume;
            c.sort_particles_by_cell = true;
            c.particle_sort_interval = 2;
        }
        if variant == 2 {
            c.field_boundary.mode = GridBoundaryMode3D::AbsorbingSponge;
            c.periodic_particles = false;
            c.particle_boundary = ParticleWallMode::Absorb;
            c.enable_secondary_yield_report = true;
            c.sort_particles_by_cell = true;
        }
        let mut st = StaggeredElectromagneticPic3D::new(c).unwrap();
        st.initialize_z_polarized_mode(25.0, 1, 0, 1).unwrap();
        let mut ps = ps3.clone();
        for p in &mut ps {
            p.velocity_m_per_s.x *= 10.0 * (variant + 1) as f64;
        }
        ps[0].position_m.x = 0.9995;
        ps[0].velocity_m_per_s.x = 6.0e7;
        ps[1].position_m.z = 0.0004;
        ps[1].velocity_m_per_s.z = -5.0e7;
        st.set_particles(ps).unwrap();
        st.step(4, None).unwrap();
        let d = st.diagnostics().unwrap();
        g.extend([
            d.field_energy_j,
            d.particle_kinetic_energy_j,
            d.charge_continuity_linf_residual,
            d.particle_count as f64,
            d.absorbed_particles as f64,
            d.secondary_macro_weight,
            d.particle_sort_passes as f64,
            d.occupied_particle_cells as f64,
        ]);
        for v in [
            st.electric_x(),
            st.electric_z(),
            st.magnetic_x(),
            st.magnetic_y(),
            st.current_x(),
            st.charge_density(),
        ] {
            g.extend_from_slice(v);
        }
        for p in st.particles() {
            p3(&mut g, p);
        }
        fld(&mut g, &st.gather_field(Vec3::new(0.41, 0.27, 0.96)));
    }
    close("pic3d_staggered", &g, 1e-10);
}

#[test]
fn pic_3d_domains() {
    let ps3 = ps3();
    let mut g = G::new();
    let c = PicDomainGrid3DConfig {
        global_nx: 6,
        global_ny: 4,
        global_nz: 4,
        domains_x: 3,
        domains_y: 2,
        domains_z: 1,
        length_x_m: 1.2,
        length_y_m: 0.8,
        length_z_m: 0.6,
        periodic: true,
    };
    let doms = make_pic_domain_grid_3d(&c).unwrap();
    for d in &doms {
        g.extend([
            d.domain_id as f64,
            d.first_x as f64,
            d.cells_x as f64,
            d.first_y as f64,
            d.cells_y as f64,
        ]);
        v3(&mut g, d.minimum_m);
        v3(&mut g, d.maximum_m);
    }
    let sorted = sort_particles_by_cell_3d(&ps3, 3, 2, 2, 1.0, 1.0, 1.0).unwrap();
    g.extend(sorted.original_indices.iter().map(|&i| i as f64));
    g.extend(sorted.cell_offsets.iter().map(|&i| i as f64));
    g.push(sorted.occupied_cells as f64);
    let hc = ParticleGuardHalo3DConfig {
        nx: 4,
        ny: 4,
        nz: 4,
        length_x_m: 1.0,
        length_y_m: 1.0,
        length_z_m: 1.0,
        guard_cells: 1,
        periodic: false,
    };
    let h = classify_particle_guard_halos_3d(&ps3, &hc).unwrap();
    for v in [
        &h.x_minus,
        &h.x_plus,
        &h.y_minus,
        &h.y_plus,
        &h.z_minus,
        &h.z_plus,
        &h.interior,
        &h.outside_domain,
    ] {
        g.push(v.len() as f64);
        g.extend(v.iter().map(|&i| i as f64));
    }
    let mut by = vec![Vec::new(); doms.len()];
    for (i, p) in ps3.iter().enumerate() {
        let mut p = *p;
        p.position_m.x *= 1.2;
        p.position_m.y *= 0.8;
        p.position_m.z *= 0.6;
        by[i % doms.len()].push(p);
    }
    let cells: Vec<Vec<f64>> = doms
        .iter()
        .enumerate()
        .map(|(d, dm)| {
            (0..dm.cells_x * dm.cells_y * dm.cells_z)
                .map(|k| 100.0 * d as f64 + k as f64)
                .collect()
        })
        .collect();
    for bl in exchange_scalar_guard_cells_3d(&cells, &doms, &c, 1, -1.0).unwrap() {
        g.extend(bl.values);
    }
    let topo = make_pic_rank_topology_3d(&doms, 2, 1, 1).unwrap();
    g.extend(topo.domain_to_rank.iter().map(|&r| r as f64));
    let round =
        run_serialized_distributed_pic_exchange_round_3d(&by, &cells, &doms, &c, 1, -2.0, &topo)
            .unwrap();
    for v in &round.particles_by_domain {
        g.push(v.len() as f64);
        for p in v {
            p3(&mut g, p);
        }
    }
    for bl in &round.scalar_guarded_blocks {
        g.extend(&bl.values);
    }
    let td = &round.transport_diagnostics;
    g.extend(
        [
            td.envelopes,
            td.particle_messages,
            td.scalar_guard_messages,
            td.particle_payload_count,
            td.scalar_guard_value_count,
            td.same_rank_messages,
            td.remote_rank_messages,
        ]
        .map(|x| x as f64),
    );
    g.extend([
        round.serialized_message_count as f64,
        round.serialized_byte_count as f64,
    ]);
    g.extend(
        round
            .serialized_plan
            .byte_count_by_destination_rank
            .iter()
            .map(|&x| x as f64),
    );
    let mut fb = ElectromagneticFieldBlocks3D::default();
    for (d, dm) in doms.iter().enumerate() {
        let n = dm.cells_x * dm.cells_y * dm.cells_z;
        let v: Vec<f64> = (0..n).map(|k| (0.3 * k as f64 + d as f64).sin()).collect();
        fb.electric_x_by_domain.push(v.clone());
        fb.electric_y_by_domain.push(vec![10.0; n]);
        fb.electric_z_by_domain.push(vec![-5.0; n]);
        fb.magnetic_x_by_domain.push(vec![1.0e-3; n]);
        fb.magnetic_y_by_domain.push(vec![0.0; n]);
        fb.magnetic_z_by_domain.push(v);
    }
    let sc = DistributedStaggeredPicStep3DConfig {
        dt_s: 1.0e-9,
        exterior_value: 0.5,
        ..Default::default()
    };
    let stp =
        run_serialized_distributed_staggered_pic_step_3d(&by, &fb, &doms, &c, &sc, &topo).unwrap();
    for v in &stp.particles_by_domain {
        g.push(v.len() as f64);
        for p in v {
            p3(&mut g, p);
        }
    }
    for bl in &stp.guarded_fields.electric_x {
        g.extend(&bl.values);
    }
    g.extend([
        stp.particle_count_before as f64,
        stp.particle_count_after_push as f64,
        stp.particle_count_after_migration as f64,
        stp.max_particle_displacement_m,
    ]);
    g.push(locate_pic_domain_3d(Vec3::new(0.9, 0.5, 0.1), &c).unwrap() as f64);
    close("pic3d_domain", &g, 1e-12);
}

#[test]
fn nuclear_plasma_wakefield() {
    let mut g = G::new();
    let pt = |energy, area| CrossSectionPoint { energy, area };
    let t = CrossSectionTable::new(vec![pt(1.0, 2.0), pt(10.0, 1.0), pt(100.0, 0.5)]).unwrap();
    for e in [0.5, 1.0, 3.0, 10.0, 55.0, 1000.0] {
        g.push(t.value(e).unwrap());
    }
    let mut dep = MultiSpeciesDepletion::new(vec![1.0e20, 0.0, 5.0e18]).unwrap();
    dep.add_reaction(DepletionReaction {
        parent: 0,
        daughter: 1,
        kind: ParticleKind::Neutron,
        cross_section: t.clone(),
        branching_fraction: 0.9,
    })
    .unwrap();
    dep.add_reaction(DepletionReaction {
        parent: 1,
        daughter: 2,
        kind: ParticleKind::Photon,
        cross_section: CrossSectionTable::new(vec![pt(1.0, 3.0e-24)]).unwrap(),
        branching_fraction: 1.0,
    })
    .unwrap();
    dep.add_reaction(DepletionReaction {
        parent: 0,
        daughter: 2,
        kind: ParticleKind::Neutron,
        cross_section: CrossSectionTable::new(vec![pt(1.0, 1.0e-24), pt(5.0, 2.0e-24)]).unwrap(),
        branching_fraction: 1.0,
    })
    .unwrap();
    for _ in 0..3 {
        let d = dep.advance(10.0, 1.0e18, 2.0e17, 4.0).unwrap();
        g.extend([
            d.total_inventory,
            d.absorbed_rate,
            d.photon_rate,
            d.neutron_rate,
        ]);
    }
    g.extend_from_slice(dep.inventory());
    let slab = SlabTransport1D::new(SlabTransportConfig {
        width: 2.0,
        total_cross_section: 1.5,
        absorption_cross_section: 0.4,
    })
    .unwrap();
    for seed in [1u64, 99] {
        let ta = slab.run(500, seed).unwrap();
        g.extend([
            ta.absorbed as f64,
            ta.transmitted as f64,
            ta.reflected as f64,
            ta.total_path,
            ta.absorption_probability(),
        ]);
    }
    close("nuclear", &g, 1e-12);

    let mut g = G::new();
    let sp = |n: &str, q, m| PlasmaSpecies {
        name: n.into(),
        charge_c: q,
        mass_kg: m,
    };
    let mut st = PlasmaState {
        species: vec![
            sp("e", -1.602176634e-19, ME),
            sp("Ar", 0.0, 6.6e-26),
            sp("Ar+", 1.602176634e-19, 6.6e-26),
        ],
        number_density_m3: vec![1.0e16, 1.0e22, 1.0e16],
        electron_temperature_ev: 3.0,
    };
    g.push(plasma_charge_density_c_m3(&st).unwrap());
    let term = |species, coefficient| PlasmaReactionTerm {
        species,
        coefficient,
    };
    let rx = vec![
        PlasmaReaction {
            reactants: vec![term(0, 1.0), term(1, 1.0)],
            products: vec![term(0, 2.0), term(2, 1.0)],
            rate_coefficient: 1.0e-20,
        },
        PlasmaReaction {
            reactants: vec![term(0, 1.0), term(2, 1.0)],
            products: vec![term(1, 1.0)],
            rate_coefficient: 2.0e-19,
        },
    ];
    let d = advance_plasma_chemistry(&mut st, &rx, 1.0e-6, 4).unwrap();
    g.extend(&d.reaction_rates_m3_s);
    g.extend([
        d.charge_density_before_c_m3,
        d.charge_density_after_c_m3,
        d.total_number_density_before_m3,
        d.total_number_density_after_m3,
        d.max_relative_density_change,
    ]);
    g.extend(&st.number_density_m3);
    let gas = PaschenGas {
        townsend_a_per_m_pa: 11.0,
        townsend_b_v_per_m_pa: 180.0,
        secondary_emission_yield: 0.01,
    };
    for pd in [100.0, 1000.0, 1.0e5] {
        g.push(paschen_breakdown_voltage_v(&gas, pd, 0.01).unwrap());
    }
    let bd = estimate_gas_breakdown_threshold(&gas, 1.0e4, 0.002, 3.0e6).unwrap();
    g.extend([
        bd.breakdown_voltage_v,
        bd.breakdown_field_v_m,
        b(bd.applied_field_exceeds_threshold),
    ]);
    for order in [1usize, 3] {
        let m = estimate_parallel_plate_multipactor(1.0e9, 0.002, order, ME, -1.602176634e-19, 1.3)
            .unwrap();
        g.extend([
            m.order as f64,
            m.transit_time_s,
            m.resonant_field_v_m,
            m.impact_energy_ev,
            b(m.secondary_yield_sustains),
        ]);
    }
    let w: Vec<f64> = (0..40)
        .map(|i| resonator_longitudinal_wake_v_per_c(0.001 * i as f64, 1.0e4, 1.3e9, 50.0).unwrap())
        .collect();
    let lam: Vec<f64> = (0..40)
        .map(|i| (-0.5 * ((i as f64 - 20.0) / 4.0).powi(2)).exp())
        .collect();
    g.extend(&w);
    g.extend(bunch_wake_potential_v(&w, &lam, 0.001).unwrap());
    let z = wake_impedance_ohm(&w, 0.001, 1.0e9).unwrap();
    g.extend([z.re, z.im]);
    close("plasma_wake", &g, 1e-11);
}

fn bx(a: [f64; 3], c: [f64; 3]) -> AxisAlignedParticleBox {
    AxisAlignedParticleBox {
        minimum_m: Vec3::new(a[0], a[1], a[2]),
        maximum_m: Vec3::new(c[0], c[1], c[2]),
    }
}

#[test]
fn transport() {
    let mut g = G::new();
    let mut w = TransportWorld {
        particles: vec![
            TransportParticleDefinition {
                name: "e-".into(),
                charge_c: -1.6e-19,
                mass_kg: ME,
            },
            TransportParticleDefinition {
                name: "gamma".into(),
                charge_c: 0.0,
                mass_kg: 1.0e-40,
            },
        ],
        materials: vec![
            TransportMaterial {
                name: "water".into(),
                density_kg_per_m3: 1000.0,
                production_cut_energy_ev: 1.0e3,
            },
            TransportMaterial {
                name: "lead".into(),
                density_kg_per_m3: 11340.0,
                production_cut_energy_ev: 0.0,
            },
        ],
        sensitive_detectors: vec![
            TransportSensitiveDetector {
                name: "det".into(),
                minimum_energy_deposit_ev: 0.0,
                record_zero_deposit_steps: false,
            },
            TransportSensitiveDetector {
                name: "strict".into(),
                minimum_energy_deposit_ev: 5.0e4,
                record_zero_deposit_steps: true,
            },
        ],
        ..Default::default()
    };
    let r1box = bx([0.1, 0.0, 0.0], [0.25, 1.0, 1.0]);
    w.regions = vec![
        TransportRegion {
            name: "slab0".into(),
            r#box: bx([0.0; 3], [0.1, 1.0, 1.0]),
            material_index: 0,
            sensitive_detector_index: Some(0),
        },
        TransportRegion {
            name: "slab1".into(),
            r#box: r1box,
            material_index: 1,
            sensitive_detector_index: Some(1),
        },
        TransportRegion {
            name: "slab2".into(),
            r#box: bx([0.25, 0.0, 0.0], [0.4, 1.0, 1.0]),
            material_index: 0,
            sensitive_detector_index: None,
        },
    ];
    let cont = TransportProcess {
        name: "ion".into(),
        stopping_power_ev_per_m: 2.0e6,
        particle_index: Some(0),
        ..Default::default()
    };
    let disc = TransportProcess {
        name: "brems".into(),
        kind: TransportProcessKind::DiscreteInteraction,
        physical_interaction_length_m: 0.02,
        energy_loss_fraction: 0.1,
        secondary_energy_fraction: 0.2,
        material_index: Some(1),
        ..Default::default()
    };
    let abs = TransportProcess {
        name: "abs".into(),
        kind: TransportProcessKind::Absorber,
        absorption_threshold_ev: 2.0e4,
        ..Default::default()
    };
    let pl = TransportPhysicsList {
        name: "std".into(),
        processes: vec![cont, disc, abs],
        production_cut_energy_ev: 500.0,
    };
    let world = apply_physics_list(w, &pl).unwrap();
    let cfg = TransportConfig {
        max_step_m: 0.03,
        ..Default::default()
    };
    let t = TransportTrack {
        position_m: Vec3::new(0.0, 0.5, 0.5),
        direction: Vec3::new(1.0, 0.2, 0.0),
        kinetic_energy_ev: 6.0e5,
        ..Default::default()
    };
    let run = transport_track(&world, t, &cfg, 200).unwrap();
    let opt = |o: Option<usize>| o.map_or(-1.0, |i| i as f64);
    let dump = |g: &mut G, r: &TransportRunResult| {
        g.extend([r.steps.len() as f64, r.secondaries.len() as f64]);
        v3(g, r.primary.position_m);
        g.extend([
            r.primary.kinetic_energy_ev,
            r.primary.time_s,
            r.primary.track_length_m,
            r.primary.status as usize as f64,
        ]);
        g.extend([
            r.scoring.total_energy_deposit_ev,
            r.scoring.total_track_length_m,
            r.scoring.steps as f64,
            r.scoring.secondaries as f64,
            r.scoring.boundary_crossings as f64,
        ]);
        for s in &r.steps {
            g.extend([
                s.step_length_m,
                s.energy_deposit_ev,
                s.limiter as usize as f64,
                s.region_index.map_or(u64::MAX as f64, |i| i as f64),
                opt(s.limiting_process_index),
            ]);
        }
        for s in &r.secondaries {
            v3(g, s.track.position_m);
            g.extend([s.track.kinetic_energy_ev, s.track.time_s]);
        }
    };
    dump(&mut g, &run);
    let mut rng = TransportRandom { state: 12345 };
    let sr = transport_track_stochastic(&world, t, &cfg, &mut rng, 200).unwrap();
    dump(&mut g, &sr);
    g.push((rng.state >> 11) as f64);
    let mut r2n = TransportRandom::default();
    for _ in 0..5 {
        g.push(sample_exponential_interaction_length_m(0.5, &mut r2n).unwrap());
    }
    let hits = collect_transport_hits(&world, &run.steps).unwrap();
    g.extend([
        hits.hits.len() as f64,
        hits.total_energy_deposit_ev,
        hits.total_weighted_energy_deposit_ev,
    ]);
    for h in &hits.hits {
        g.push(h.detector_index as f64);
        v3(&mut g, h.position_m);
        g.push(h.energy_deposit_ev);
    }
    let mut grid =
        make_transport_dose_grid(4, 3, 2, bx([0.0; 3], [0.4, 1.0, 1.0]), 1000.0).unwrap();
    score_hits_to_dose_grid(&hits, &mut grid).unwrap();
    g.extend(&grid.energy_deposit_ev);
    g.extend(&grid.dose_gy);
    g.extend(dose_rate_to_sar_w_per_kg(&grid, 2.0).unwrap());
    g.extend(project_dose_grid_to_pennes2d_sar(&grid, 2.0, 3, 2).unwrap());
    g.push(
        distance_to_box_boundary(&r1box, Vec3::new(0.15, 0.5, 0.5), Vec3::new(-1.0, 1.0, 0.0))
            .unwrap(),
    );
    let mut big = world.clone();
    big.regions.clear();
    for i in 0..4 {
        for j in 0..3 {
            let (fi, fj) = (i as f64, j as f64);
            big.regions.push(TransportRegion {
                name: "c".into(),
                r#box: bx(
                    [0.25 * fi, 0.3 * fj, 0.0],
                    [0.25 * (fi + 1.0), 0.3 * (fj + 1.0), 1.0],
                ),
                material_index: (i + j) % 2,
                sensitive_detector_index: None,
            });
        }
    }
    let bvh = build_region_bvh(&big, 2).unwrap();
    for p in [
        Vec3::new(0.1, 0.1, 0.5),
        Vec3::new(0.6, 0.45, 0.5),
        Vec3::new(0.99, 0.89, 0.2),
        Vec3::new(1.5, 0.1, 0.1),
    ] {
        g.push(opt(locate_region(&big, p, 0.0)));
        g.push(opt(locate_region_bvh(&big, &bvh, p, 0.0)));
    }
    g.push(bvh.nodes.len() as f64);
    close("transport", &g, 1e-12);
}
