//! Local model closures. A closure alone is not a coupled RANS/LES/multiphase solver.
use super::{least_squares_gradient_scalar, operators::vector_sizes, PolyMesh};
use crate::{
    core::{dense_solve, Vec3},
    finite_slice, positive, require, Result, SolverError,
};
use serde::{Deserialize, Serialize};
fn nonnegative(v: f64, name: &str) -> Result<()> {
    require(
        v.is_finite() && v >= 0.0,
        format!("{name} must be nonnegative and finite"),
    )
}
pub fn spalart_allmaras_eddy_viscosity(nu: f64, nu_tilde: f64, cv1: f64) -> Result<f64> {
    positive(nu, "molecular viscosity")?;
    nonnegative(nu_tilde, "nu tilde")?;
    positive(cv1, "cv1")?;
    let chi = nu_tilde / nu;
    let c3 = chi.powi(3);
    let v = nu_tilde * c3 / (c3 + cv1.powi(3));
    require(v.is_finite(), "SA closure overflow")?;
    Ok(v)
}
pub fn k_epsilon_eddy_viscosity(density: f64, k: f64, epsilon: f64, c_mu: f64) -> Result<f64> {
    positive(density, "density")?;
    nonnegative(k, "k")?;
    positive(epsilon, "epsilon")?;
    positive(c_mu, "c_mu")?;
    let v = density * c_mu * k * k / epsilon;
    require(v.is_finite(), "k-epsilon closure overflow")?;
    Ok(v)
}
pub fn k_omega_sst_eddy_viscosity(
    density: f64,
    k: f64,
    omega: f64,
    strain: f64,
    f2: f64,
    a1: f64,
) -> Result<f64> {
    positive(density, "density")?;
    nonnegative(k, "k")?;
    positive(omega, "omega")?;
    nonnegative(strain, "strain")?;
    nonnegative(f2, "f2")?;
    positive(a1, "a1")?;
    let v = density * a1 * k / (a1 * omega).max(strain * f2);
    require(v.is_finite(), "SST closure overflow")?;
    Ok(v)
}
pub fn des_length_scale(wall: f64, grid: f64, c_des: f64) -> Result<f64> {
    nonnegative(wall, "wall distance")?;
    positive(grid, "grid length")?;
    positive(c_des, "DES coefficient")?;
    Ok(wall.min(c_des * grid))
}
#[derive(Clone, Copy, Debug)]
pub struct ReynoldsStress2D {
    pub xx: f64,
    pub yy: f64,
    pub xy: f64,
}
pub fn boussinesq_reynolds_stress(
    density: f64,
    k: f64,
    eddy_viscosity: f64,
    du_dx: f64,
    dv_dy: f64,
    shear: f64,
) -> Result<ReynoldsStress2D> {
    positive(density, "density")?;
    nonnegative(k, "k")?;
    nonnegative(eddy_viscosity, "eddy viscosity")?;
    finite_slice(&[du_dx, dv_dy, shear], "velocity gradient")?;
    let iso = 2.0 * density * k / 3.0;
    Ok(ReynoldsStress2D {
        xx: iso - 2.0 * eddy_viscosity * du_dx,
        yy: iso - 2.0 * eddy_viscosity * dv_dy,
        xy: -eddy_viscosity * shear,
    })
}
#[derive(Clone, Copy, Debug)]
pub struct PhaseChangeState {
    pub liquid_fraction: f64,
    pub effective_cp: f64,
    pub enthalpy: f64,
}
pub fn linear_mushy_phase_change(
    t: f64,
    solidus: f64,
    liquidus: f64,
    cp_solid: f64,
    cp_liquid: f64,
    latent: f64,
) -> Result<PhaseChangeState> {
    require(
        t.is_finite() && solidus.is_finite() && liquidus.is_finite() && liquidus > solidus,
        "invalid phase interval/temperature",
    )?;
    positive(cp_solid, "solid cp")?;
    positive(cp_liquid, "liquid cp")?;
    nonnegative(latent, "latent heat")?;
    let f = ((t - solidus) / (liquidus - solidus)).clamp(0.0, 1.0);
    let cp = (1.0 - f) * cp_solid
        + f * cp_liquid
        + if t > solidus && t < liquidus {
            latent / (liquidus - solidus)
        } else {
            0.0
        };
    let mut h = cp_solid * t.min(solidus);
    if t > solidus {
        let tm = t.min(liquidus) - solidus;
        h += 0.5 * (cp_solid + cp_liquid) * tm
            + latent * (tm / (liquidus - solidus)).clamp(0.0, 1.0);
    }
    if t > liquidus {
        h += cp_liquid * (t - liquidus);
    }
    Ok(PhaseChangeState {
        liquid_fraction: f,
        effective_cp: cp,
        enthalpy: h,
    })
}
pub fn drift_flux_slip_velocity(
    rho_c: f64,
    rho_d: f64,
    diameter: f64,
    viscosity: f64,
    gravity: f64,
) -> Result<f64> {
    for (v, n) in [
        (rho_c, "continuous density"),
        (rho_d, "dispersed density"),
        (diameter, "diameter"),
        (viscosity, "viscosity"),
    ] {
        positive(v, n)?;
    }
    require(gravity.is_finite(), "nonfinite gravity")?;
    Ok((rho_c - rho_d) * gravity * diameter * diameter / (18.0 * viscosity))
}
pub fn euler_euler_drag_source(
    continuous: f64,
    dispersed: f64,
    coefficient: f64,
    area: f64,
) -> Result<f64> {
    nonnegative(coefficient, "drag coefficient")?;
    nonnegative(area, "interface area")?;
    finite_slice(&[continuous, dispersed], "velocities")?;
    let d = dispersed - continuous;
    Ok(coefficient * area * d * d.abs())
}
pub fn lubrication_thin_film_flux(
    thickness: f64,
    pressure_gradient: f64,
    viscosity: f64,
    surface_velocity: f64,
) -> Result<f64> {
    nonnegative(thickness, "film thickness")?;
    positive(viscosity, "viscosity")?;
    finite_slice(&[pressure_gradient, surface_velocity], "thin-film inputs")?;
    Ok(-thickness.powi(3) * pressure_gradient / (12.0 * viscosity)
        + 0.5 * thickness * surface_velocity)
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TurbulenceModel {
    Laminar,
    MixingLength,
    Smagorinsky,
    Wale,
}
pub fn parse_turbulence_model(name: &str) -> Result<TurbulenceModel> {
    Ok(match name {
        "laminar" => TurbulenceModel::Laminar,
        "mixingLength" | "mixing_length" => TurbulenceModel::MixingLength,
        "Smagorinsky" | "smagorinsky" => TurbulenceModel::Smagorinsky,
        "WALE" | "wale" => TurbulenceModel::Wale,
        _ => {
            return Err(SolverError::InvalidInput(format!(
                "unknown turbulence model {name}"
            )))
        }
    })
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct TurbulenceConfig {
    pub model: TurbulenceModel,
    pub mixing_length: f64,
    pub smagorinsky_constant: f64,
    pub wale_constant: f64,
}
impl Default for TurbulenceConfig {
    fn default() -> Self {
        Self {
            model: TurbulenceModel::Laminar,
            mixing_length: 0.01,
            smagorinsky_constant: 0.17,
            wale_constant: 0.325,
        }
    }
}
pub(crate) fn velocity_gradients(mesh: &PolyMesh, u: &[Vec3]) -> Result<Vec<[[f64; 3]; 3]>> {
    vector_sizes(mesh, u, &[])?;
    let mut g: Vec<Vec<Vec3>> = vec![];
    for a in 0..3 {
        let v: Vec<_> = u.iter().map(|v| v.component(a)).collect();
        g.push(least_squares_gradient_scalar(mesh, &v, &[])?);
    }
    Ok((0..mesh.cell_count())
        .map(|i| [g[0][i].as_array(), g[1][i].as_array(), g[2][i].as_array()])
        .collect())
}
pub fn eddy_viscosity(mesh: &PolyMesh, u: &[Vec3], c: TurbulenceConfig) -> Result<Vec<f64>> {
    nonnegative(c.mixing_length, "mixing length")?;
    nonnegative(c.smagorinsky_constant, "Smagorinsky coefficient")?;
    nonnegative(c.wale_constant, "WALE coefficient")?;
    vector_sizes(mesh, u, &[])?;
    let mut out = vec![0.0; mesh.cell_count()];
    if c.model == TurbulenceModel::Laminar {
        return Ok(out);
    }
    let gradients = velocity_gradients(mesh, u)?;
    for (n, g) in gradients.iter().enumerate() {
        let mut ss = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                ss += (0.5 * (g[i][j] + g[j][i])).powi(2);
            }
        }
        let mag_s = (2.0 * ss).sqrt();
        let delta = mesh.cells()[n].volume.cbrt();
        out[n] = match c.model {
            TurbulenceModel::MixingLength => c.mixing_length.powi(2) * mag_s,
            TurbulenceModel::Smagorinsky => (c.smagorinsky_constant * delta).powi(2) * mag_s,
            TurbulenceModel::Wale => {
                let mut g2 = [[0.0; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        for k in 0..3 {
                            g2[i][j] += g[i][k] * g[k][j];
                        }
                    }
                }
                let tr = (0..3).map(|i| g2[i][i]).sum::<f64>();
                let mut sd2 = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        let sd = 0.5 * (g2[i][j] + g2[j][i]) - if i == j { tr / 3.0 } else { 0.0 };
                        sd2 += sd * sd;
                    }
                }
                (c.wale_constant * delta).powi(2) * sd2.powf(1.5)
                    / (ss.powf(2.5) + sd2.powf(1.25) + 1e-30)
            }
            TurbulenceModel::Laminar => 0.0,
        };
    }
    finite_slice(&out, "eddy viscosity")?;
    Ok(out)
}
pub const STEFAN_BOLTZMANN: f64 = 5.670374419e-8;
pub fn gray_surface_exchange(ta: f64, ea: f64, tb: f64, eb: f64, view_factor: f64) -> Result<f64> {
    nonnegative(ta, "temperature A")?;
    nonnegative(tb, "temperature B")?;
    require(
        ea > 0.0 && ea <= 1.0 && eb > 0.0 && eb <= 1.0 && (0.0..=1.0).contains(&view_factor),
        "invalid emissivity/view factor",
    )?;
    let q =
        view_factor / (1.0 / ea + 1.0 / eb - 1.0) * STEFAN_BOLTZMANN * (ta.powi(4) - tb.powi(4));
    require(q.is_finite(), "radiation overflow")?;
    Ok(q)
}
pub fn optically_thin_radiation_source(t: f64, environment: f64, absorption: f64) -> Result<f64> {
    nonnegative(t, "temperature")?;
    nonnegative(environment, "environment temperature")?;
    nonnegative(absorption, "absorption")?;
    let q = -4.0 * absorption * STEFAN_BOLTZMANN * (t.powi(4) - environment.powi(4));
    require(q.is_finite(), "radiation source overflow")?;
    Ok(q)
}
pub trait RadiationSourceModel {
    fn volumetric_source(&self, temperature: f64, environment_temperature: f64) -> Result<f64>;
}
#[derive(Clone, Copy, Debug)]
pub struct OpticallyThinRadiation {
    absorption: f64,
}
impl OpticallyThinRadiation {
    pub fn new(absorption: f64) -> Result<Self> {
        nonnegative(absorption, "absorption")?;
        Ok(Self { absorption })
    }
}
impl RadiationSourceModel for OpticallyThinRadiation {
    fn volumetric_source(&self, t: f64, environment: f64) -> Result<f64> {
        optically_thin_radiation_source(t, environment, self.absorption)
    }
}
#[derive(Clone, Debug)]
pub struct RadiosityResult {
    pub radiosity: Vec<f64>,
    pub net_heat: Vec<f64>,
}
pub fn solve_gray_radiosity(
    t: &[f64],
    emissivity: &[f64],
    area: &[f64],
    view_factors: &[f64],
) -> Result<RadiosityResult> {
    let n = t.len();
    require(
        n > 0
            && emissivity.len() == n
            && area.len() == n
            && n.checked_mul(n) == Some(view_factors.len()),
        "radiosity dimensions",
    )?;
    finite_slice(view_factors, "view factors")?;
    let mut matrix = vec![0.0; n * n];
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        nonnegative(t[i], "radiation temperature")?;
        positive(area[i], "radiation area")?;
        require(
            emissivity[i] > 0.0 && emissivity[i] <= 1.0,
            "surface emissivity",
        )?;
        let mut sum = 0.0;
        for j in 0..n {
            let f = view_factors[i * n + j];
            nonnegative(f, "view factor")?;
            sum += f;
            matrix[i * n + j] = (if i == j { 1.0 } else { 0.0 }) - (1.0 - emissivity[i]) * f;
        }
        require((sum - 1.0).abs() <= 1e-8, "view-factor row must sum to one")?;
        rhs[i] = emissivity[i] * STEFAN_BOLTZMANN * t[i].powi(4);
    }
    let radiosity = dense_solve(&matrix, &rhs)?;
    let net_heat = (0..n)
        .map(|i| {
            area[i]
                * (radiosity[i]
                    - (0..n)
                        .map(|j| view_factors[i * n + j] * radiosity[j])
                        .sum::<f64>())
        })
        .collect();
    Ok(RadiosityResult {
        radiosity,
        net_heat,
    })
}
