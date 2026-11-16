//! Acoustic-to-structure loads, acoustic heating of Pennes tissue, streaming body force, sensor
//! misfit and coordinate-search parameter fitting, ported from upstream
//! `cfd/acoustics/coupling` (MIT, cfd_solvers contributors).
use crate::{
    core::Vec3,
    fem::{DynamicBar1D, Mesh2D},
    multiphysics::PennesBioheat2D,
    require,
    workflow::inverse_coordinate_search,
    Result,
};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StructuralNodalForce2D {
    pub x: f64,
    pub y: f64,
}

/// Traction `-p n` over boundary edges (all patches when `patch` is `None`), split equally to
/// the two edge nodes; `n` is the right-hand normal of the stored edge.
pub fn pressure_to_structural_boundary_load(
    mesh: &Mesh2D,
    edge_pressure_pa: &[f64],
    patch: Option<i32>,
) -> Result<Vec<StructuralNodalForce2D>> {
    mesh.validate()?;
    require(
        edge_pressure_pa.len() == mesh.boundary_edges.len(),
        "boundary pressure size mismatch",
    )?;
    let mut out = vec![StructuralNodalForce2D::default(); mesh.node_count()];
    for (edge, &p) in mesh.boundary_edges.iter().zip(edge_pressure_pa) {
        if patch.is_some_and(|q| edge.patch != q) {
            continue;
        }
        require(p.is_finite(), "non-finite acoustic boundary pressure")?;
        let (a, b) = (mesh.nodes[edge.node[0]], mesh.nodes[edge.node[1]]);
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let length = dx.hypot(dy);
        require(length > 0.0, "degenerate FEM boundary edge")?;
        let (nx, ny) = (dy / length, -dx / length);
        let (fx, fy) = (-p * nx * length * 0.5, -p * ny * length * 0.5);
        for n in edge.node {
            out[n].x += fx;
            out[n].y += fy;
        }
    }
    Ok(out)
}

pub fn acoustic_pressure_to_bar_tip_force(pressure_pa: f64, loaded_area_m2: f64) -> Result<f64> {
    require(
        pressure_pa.is_finite() && loaded_area_m2 >= 0.0 && loaded_area_m2.is_finite(),
        "invalid acoustic bar loading",
    )?;
    Ok(pressure_pa * loaded_area_m2)
}

/// One Newmark step per pressure sample.
pub fn drive_dynamic_bar_from_pressure_trace(
    bar: &mut DynamicBar1D,
    pressure_pa: &[f64],
    loaded_area_m2: f64,
) -> Result<()> {
    require(
        loaded_area_m2 >= 0.0 && loaded_area_m2.is_finite(),
        "invalid bar loading area",
    )?;
    for &p in pressure_pa {
        bar.step(acoustic_pressure_to_bar_tip_force(p, loaded_area_m2)?)?;
    }
    Ok(())
}

/// Sets `q = 2 alpha I` as the Pennes volumetric heating.
pub fn apply_acoustic_heating_to_pennes(
    solver: &mut PennesBioheat2D,
    intensity_w_m2: &[f64],
    absorption_np_per_m: &[f64],
) -> Result<()> {
    require(
        intensity_w_m2.len() == absorption_np_per_m.len()
            && intensity_w_m2.len() == solver.temperature_k().len(),
        "acoustic/Pennes grid size mismatch",
    )?;
    let mut heating = Vec::with_capacity(intensity_w_m2.len());
    for (&i, &a) in intensity_w_m2.iter().zip(absorption_np_per_m) {
        require(
            i >= 0.0 && a >= 0.0 && i.is_finite() && a.is_finite(),
            "invalid acoustic heating input",
        )?;
        heating.push(2.0 * a * i);
    }
    solver.set_volumetric_heating(&heating)
}

/// Eckart streaming acceleration `2 alpha I / (rho c)` along the propagation direction.
pub fn acoustic_streaming_acceleration(
    intensity_w_m2: &[f64],
    absorption_np_per_m: &[f64],
    density_kg_m3: &[f64],
    sound_speed_m_s: &[f64],
    propagation_direction: Vec3,
) -> Result<Vec<Vec3>> {
    let n = intensity_w_m2.len();
    require(
        absorption_np_per_m.len() == n && density_kg_m3.len() == n && sound_speed_m_s.len() == n,
        "acoustic streaming field size mismatch",
    )?;
    let mag = propagation_direction.norm();
    require(
        mag > 0.0 && mag.is_finite(),
        "invalid acoustic propagation direction",
    )?;
    let dir = propagation_direction / mag;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let (it, a, rho, c) = (
            intensity_w_m2[i],
            absorption_np_per_m[i],
            density_kg_m3[i],
            sound_speed_m_s[i],
        );
        require(
            it >= 0.0
                && a >= 0.0
                && rho > 0.0
                && c > 0.0
                && it.is_finite()
                && a.is_finite()
                && rho.is_finite()
                && c.is_finite(),
            "invalid acoustic streaming material/field value",
        )?;
        out.push(dir * (2.0 * a * it / (rho * c)));
    }
    Ok(out)
}

/// `sqrt(sum (p - m)^2)`, normalized by `sum m^2` when requested and non-zero.
pub fn sensor_trace_l2_misfit(
    predicted: &[Vec<f64>],
    measured: &[Vec<f64>],
    normalize: bool,
) -> Result<f64> {
    require(
        predicted.len() == measured.len() && !predicted.is_empty(),
        "sensor trace set mismatch",
    )?;
    let (mut error, mut reference) = (0.0, 0.0);
    for (p, m) in predicted.iter().zip(measured) {
        require(
            p.len() == m.len() && !p.is_empty(),
            "sensor trace length mismatch",
        )?;
        for (a, b) in p.iter().zip(m) {
            require(a.is_finite() && b.is_finite(), "non-finite sensor trace")?;
            error += (a - b) * (a - b);
            reference += b * b;
        }
    }
    if normalize && reference > 0.0 {
        error /= reference;
    }
    Ok(error.sqrt())
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AcousticInverseResult {
    pub parameters: Vec<f64>,
    pub misfit: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// Fits simulator parameters to measured traces by compass search on the normalized misfit
/// (upstream defaults 100 iterations, tolerance 1e-8). Simulator or misfit errors count as an
/// infinite objective.
pub fn fit_acoustic_parameters(
    simulator: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    measured: &[Vec<f64>],
    initial: Vec<f64>,
    step: Vec<f64>,
    max_iterations: usize,
    tolerance: f64,
) -> Result<AcousticInverseResult> {
    let r = inverse_coordinate_search(
        |x| sensor_trace_l2_misfit(&simulator(x), measured, true).unwrap_or(f64::INFINITY),
        initial,
        step,
        max_iterations,
        tolerance,
    )?;
    Ok(AcousticInverseResult {
        parameters: r.parameters,
        misfit: r.objective,
        iterations: r.iterations,
        converged: r.converged,
    })
}
