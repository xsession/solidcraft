//! cargo run --release --example thermoelastic_plate -- results/plate
//! One-way thermal-to-structural transfer on a shared unstructured mesh.
use cfd_solvers::{
    core::Vec3,
    fem::*,
    io::{self, DataArray, VtuMesh},
    Result,
};
fn main() -> Result<()> {
    let out = std::path::PathBuf::from(
        std::env::args()
            .nth(1)
            .unwrap_or_else(|| "results/plate".into()),
    );
    std::fs::create_dir_all(&out)?;
    let mesh = make_rectangle_tri_mesh(16, 12, 0.1, 0.06)?;
    let mut heat = Heat2D::new(
        mesh.clone(),
        Heat2DConfig {
            conductivity: 205.0,
            volumetric_heat_capacity: 2.43e6,
            dt: 0.5,
            ..Default::default()
        },
    )?;
    heat.initialize(|_| 293.15)?;
    for _ in 0..5 {
        heat.step(|_, _, _| 1e5, |_, _| 343.15)?;
    }
    let mut structure = Elasticity2D::new(
        mesh,
        Elasticity2DConfig {
            youngs_modulus: 70e9,
            poisson_ratio: 0.33,
            ..Default::default()
        },
    )?;
    structure.set_dirichlet(
        |p| p.x == 0.0 && p.y == 0.0,
        |_| Displacement2::default(),
        true,
        true,
    )?;
    structure.set_dirichlet(
        |p| p.x == 0.1 && p.y == 0.0,
        |_| Displacement2::default(),
        false,
        true,
    )?;
    structure.solve_thermal(&heat.temperature, 23e-6, 293.15, |_| {
        Displacement2::default()
    })?;
    let displacement: Vec<_> = structure
        .displacement()
        .iter()
        .map(|u| Vec3::new(u.x, u.y, 0.0))
        .collect();
    io::write_vtu(
        out.join("thermoelastic.vtu"),
        &VtuMesh::from_triangles(structure.mesh())?,
        &[
            DataArray::scalar("temperature_K", heat.temperature),
            DataArray::vector("displacement_m", displacement),
        ],
        &[],
    )?;
    println!("{}", out.join("thermoelastic.vtu").display());
    Ok(())
}
