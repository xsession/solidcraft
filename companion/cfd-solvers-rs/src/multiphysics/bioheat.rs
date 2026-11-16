//! Bioelectromagnetic heating: SAR from RMS fields, mass-averaged peak SAR on voxel
//! tissue grids, temperature feedback, voxel-to-2D projection and an implicit Pennes
//! bioheat solver (periodic or fixed-temperature boundary).
use crate::{
    core::{
        preconditioned_conjugate_gradient, IdentityPreconditioner, IterativeSolverResult,
        LinearConfig,
    },
    require, Result, SolverError,
};

fn finite_nonneg(v: f64) -> bool {
    v.is_finite() && v >= 0.0
}

/// `SAR = sigma |E_rms|^2 / rho` (W/kg).
pub fn sar_from_rms_electric_field(
    conductivity: f64,
    density: f64,
    electric_rms: f64,
) -> Result<f64> {
    require(
        finite_nonneg(conductivity)
            && density.is_finite()
            && density > 0.0
            && finite_nonneg(electric_rms),
        "invalid SAR material/field data",
    )?;
    Ok(conductivity * electric_rms * electric_rms / density)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VoxelTissueMaterial {
    pub density_kg_per_m3: f64,
    pub specific_heat_j_per_kg_k: f64,
    pub thermal_conductivity_w_per_m_k: f64,
    pub electrical_conductivity_s_per_m: f64,
    pub blood_perfusion_per_s: f64,
    pub metabolic_heat_w_per_m3: f64,
}

impl Default for VoxelTissueMaterial {
    fn default() -> Self {
        Self {
            density_kg_per_m3: 1000.0,
            specific_heat_j_per_kg_k: 3600.0,
            thermal_conductivity_w_per_m_k: 0.5,
            electrical_conductivity_s_per_m: 1.0,
            blood_perfusion_per_s: 0.0,
            metabolic_heat_w_per_m3: 0.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct VoxelTissueGrid3D {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx_m: f64,
    pub dy_m: f64,
    pub dz_m: f64,
    pub materials: Vec<VoxelTissueMaterial>,
    pub material_index: Vec<usize>,
    pub electric_rms_v_per_m: Vec<f64>,
}

impl VoxelTissueGrid3D {
    fn count(&self) -> Option<usize> {
        self.nx.checked_mul(self.ny)?.checked_mul(self.nz)
    }
    pub fn validate(&self) -> Result<()> {
        require(
            self.nx > 0
                && self.ny > 0
                && self.nz > 0
                && [self.dx_m, self.dy_m, self.dz_m]
                    .iter()
                    .all(|d| d.is_finite() && *d > 0.0)
                && !self.materials.is_empty(),
            "invalid voxel tissue grid geometry",
        )?;
        let n = self
            .count()
            .ok_or_else(|| SolverError::InvalidInput("voxel grid overflow".into()))?;
        require(
            self.material_index.len() == n && self.electric_rms_v_per_m.len() == n,
            "voxel tissue field size mismatch",
        )?;
        require(
            self.materials.iter().all(|m| {
                m.density_kg_per_m3 > 0.0
                    && m.specific_heat_j_per_kg_k > 0.0
                    && m.thermal_conductivity_w_per_m_k >= 0.0
                    && m.electrical_conductivity_s_per_m >= 0.0
                    && m.blood_perfusion_per_s >= 0.0
            }),
            "invalid voxel tissue material",
        )?;
        require(
            self.material_index
                .iter()
                .all(|i| *i < self.materials.len())
                && self.electric_rms_v_per_m.iter().all(|e| finite_nonneg(*e)),
            "invalid voxel tissue index/field",
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureFeedbackConfig {
    pub reference_temperature_k: f64,
    pub conductivity_reference_s_per_m: f64,
    pub conductivity_coefficient_per_k: f64,
    pub perfusion_reference_per_s: f64,
    pub perfusion_coefficient_per_k: f64,
}

impl Default for TemperatureFeedbackConfig {
    fn default() -> Self {
        Self {
            reference_temperature_k: 310.15,
            conductivity_reference_s_per_m: 1.0,
            conductivity_coefficient_per_k: 0.01,
            perfusion_reference_per_s: 0.01,
            perfusion_coefficient_per_k: 0.02,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureFeedbackValues {
    pub conductivity_s_per_m: f64,
    pub perfusion_per_s: f64,
}

/// Linear temperature feedback on conductivity and perfusion, clamped at zero.
pub fn temperature_feedback(
    c: &TemperatureFeedbackConfig,
    temperature: f64,
) -> Result<TemperatureFeedbackValues> {
    require(
        c.reference_temperature_k.is_finite()
            && c.reference_temperature_k > 0.0
            && finite_nonneg(c.conductivity_reference_s_per_m)
            && finite_nonneg(c.perfusion_reference_per_s)
            && c.conductivity_coefficient_per_k.is_finite()
            && c.perfusion_coefficient_per_k.is_finite()
            && finite_nonneg(temperature),
        "invalid bioelectromagnetic temperature feedback",
    )?;
    let delta = temperature - c.reference_temperature_k;
    Ok(TemperatureFeedbackValues {
        conductivity_s_per_m: (c.conductivity_reference_s_per_m
            * (1.0 + c.conductivity_coefficient_per_k * delta))
            .max(0.0),
        perfusion_per_s: (c.perfusion_reference_per_s
            * (1.0 + c.perfusion_coefficient_per_k * delta))
            .max(0.0),
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExposureScenario {
    Implant,
    Wearable,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BioelectromagneticExposureCase {
    pub name: String,
    pub scenario: ExposureScenario,
    pub grid: VoxelTissueGrid3D,
    pub averaging_mass_kg: f64,
    pub limit_w_per_kg: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureValidationResult {
    pub maximum_mass_averaged_sar_w_per_kg: f64,
    pub within_limit: bool,
}

pub fn validate_exposure_case(
    exposure: &BioelectromagneticExposureCase,
) -> Result<ExposureValidationResult> {
    require(
        !exposure.name.is_empty()
            && exposure.averaging_mass_kg.is_finite()
            && exposure.averaging_mass_kg > 0.0
            && exposure.limit_w_per_kg.is_finite()
            && exposure.limit_w_per_kg > 0.0,
        "invalid bioelectromagnetic exposure case",
    )?;
    let sar = voxel_sar_w_per_kg(&exposure.grid)?;
    let maximum = max_mass_averaged_sar_w_per_kg(&exposure.grid, &sar, exposure.averaging_mass_kg)?;
    Ok(ExposureValidationResult {
        maximum_mass_averaged_sar_w_per_kg: maximum,
        within_limit: maximum <= exposure.limit_w_per_kg,
    })
}

pub fn voxel_sar_w_per_kg(grid: &VoxelTissueGrid3D) -> Result<Vec<f64>> {
    grid.validate()?;
    grid.material_index
        .iter()
        .zip(&grid.electric_rms_v_per_m)
        .map(|(m, e)| {
            let m = grid.materials[*m];
            sar_from_rms_electric_field(m.electrical_conductivity_s_per_m, m.density_kg_per_m3, *e)
        })
        .collect()
}

/// Peak SAR averaged over the nearest voxels (by centre distance) that make up the
/// averaging mass around every voxel, taking a fraction of the last voxel's mass.
pub fn max_mass_averaged_sar_w_per_kg(
    grid: &VoxelTissueGrid3D,
    sar: &[f64],
    averaging_mass_kg: f64,
) -> Result<f64> {
    grid.validate()?;
    let n = grid.material_index.len();
    require(
        sar.len() == n && averaging_mass_kg.is_finite() && averaging_mass_kg > 0.0,
        "invalid SAR averaging request",
    )?;
    require(sar.iter().all(|s| finite_nonneg(*s)), "invalid SAR value")?;
    let volume = grid.dx_m * grid.dy_m * grid.dz_m;
    let plane = grid.nx * grid.ny;
    let coords = |i: usize| (i % grid.nx, (i % plane) / grid.nx, i / plane);
    let mut best = 0.0f64;
    let mut order: Vec<(f64, usize)> = Vec::with_capacity(n);
    for centre in 0..n {
        let (cx, cy, cz) = coords(centre);
        order.clear();
        for idx in 0..n {
            let (x, y, z) = coords(idx);
            let dx = (x as f64 - cx as f64) * grid.dx_m;
            let dy = (y as f64 - cy as f64) * grid.dy_m;
            let dz = (z as f64 - cz as f64) * grid.dz_m;
            order.push((dx * dx + dy * dy + dz * dz, idx));
        }
        order.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (mut mass, mut weighted) = (0.0, 0.0);
        for (_, idx) in &order {
            let voxel_mass = grid.materials[grid.material_index[*idx]].density_kg_per_m3 * volume;
            let take = voxel_mass.min(averaging_mass_kg - mass);
            weighted += take * sar[*idx];
            mass += take;
            if mass >= averaging_mass_kg - 1e-18 {
                break;
            }
        }
        if mass > 0.0 {
            best = best.max(weighted / mass);
        }
    }
    Ok(best)
}

/// Mass-weighted projection of voxel SAR onto an `nx` x `ny` grid (z collapsed).
pub fn project_voxel_sar_to_pennes2d(
    grid: &VoxelTissueGrid3D,
    sar: &[f64],
    nx: usize,
    ny: usize,
) -> Result<Vec<f64>> {
    grid.validate()?;
    require(
        sar.len() == grid.material_index.len() && nx > 0 && ny > 0,
        "invalid voxel-to-Pennes projection",
    )?;
    let mut projected = vec![0.0; nx * ny];
    let mut mass = vec![0.0; nx * ny];
    let volume = grid.dx_m * grid.dy_m * grid.dz_m;
    for z in 0..grid.nz {
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                let src = (z * grid.ny + y) * grid.nx + x;
                let (tx, ty) = (
                    (x * nx / grid.nx).min(nx - 1),
                    (y * ny / grid.ny).min(ny - 1),
                );
                let m = grid.materials[grid.material_index[src]].density_kg_per_m3 * volume;
                projected[ty * nx + tx] += m * sar[src];
                mass[ty * nx + tx] += m;
            }
        }
    }
    for (p, m) in projected.iter_mut().zip(&mass) {
        if *m > 0.0 {
            *p /= m;
        }
    }
    Ok(projected)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BioheatBoundary {
    Periodic,
    FixedTemperature,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PennesBioheat2DConfig {
    pub nx: usize,
    pub ny: usize,
    pub width_m: f64,
    pub height_m: f64,
    pub tissue_density_kg_per_m3: f64,
    pub tissue_specific_heat_j_per_kg_k: f64,
    pub thermal_conductivity_w_per_m_k: f64,
    pub blood_density_kg_per_m3: f64,
    pub blood_specific_heat_j_per_kg_k: f64,
    pub blood_perfusion_per_s: f64,
    pub blood_temperature_k: f64,
    pub metabolic_heat_w_per_m3: f64,
    pub dt_s: f64,
    pub boundary: BioheatBoundary,
    pub fixed_boundary_temperature_k: f64,
    pub max_iterations: usize,
    pub relative_tolerance: f64,
}

impl Default for PennesBioheat2DConfig {
    fn default() -> Self {
        Self {
            nx: 32,
            ny: 32,
            width_m: 0.1,
            height_m: 0.1,
            tissue_density_kg_per_m3: 1000.0,
            tissue_specific_heat_j_per_kg_k: 3600.0,
            thermal_conductivity_w_per_m_k: 0.5,
            blood_density_kg_per_m3: 1060.0,
            blood_specific_heat_j_per_kg_k: 3770.0,
            blood_perfusion_per_s: 0.0,
            blood_temperature_k: 310.15,
            metabolic_heat_w_per_m3: 0.0,
            dt_s: 0.1,
            boundary: BioheatBoundary::Periodic,
            fixed_boundary_temperature_k: 310.15,
            max_iterations: 4000,
            relative_tolerance: 1e-11,
        }
    }
}

/// Implicit-Euler Pennes bioheat equation on a uniform node grid with a 5-point stencil.
pub struct PennesBioheat2D {
    config: PennesBioheat2DConfig,
    temperature: Vec<f64>,
    sar: Vec<f64>,
    heating: Vec<f64>,
    linear: Option<IterativeSolverResult>,
    time_s: f64,
}

impl PennesBioheat2D {
    pub fn new(config: PennesBioheat2DConfig) -> Result<Self> {
        let c = &config;
        require(
            c.nx >= 3
                && c.ny >= 3
                && c.width_m > 0.0
                && c.height_m > 0.0
                && c.tissue_density_kg_per_m3 > 0.0
                && c.tissue_specific_heat_j_per_kg_k > 0.0
                && c.thermal_conductivity_w_per_m_k >= 0.0
                && c.blood_density_kg_per_m3 > 0.0
                && c.blood_specific_heat_j_per_kg_k > 0.0
                && c.blood_perfusion_per_s >= 0.0
                && c.dt_s > 0.0
                && c.max_iterations > 0
                && c.relative_tolerance > 0.0,
            "invalid Pennes bioheat configuration",
        )?;
        let n =
            c.nx.checked_mul(c.ny)
                .ok_or_else(|| SolverError::InvalidInput("bioheat grid overflow".into()))?;
        Ok(Self {
            config,
            temperature: vec![config.blood_temperature_k; n],
            sar: vec![0.0; n],
            heating: vec![0.0; n],
            linear: None,
            time_s: 0.0,
        })
    }
    pub fn initialize(&mut self, temperature: f64) -> Result<()> {
        require(
            temperature.is_finite(),
            "invalid initial tissue temperature",
        )?;
        self.temperature.fill(temperature);
        self.time_s = 0.0;
        Ok(())
    }
    pub fn set_sar(&mut self, values: &[f64]) -> Result<()> {
        require(
            values.len() == self.sar.len() && values.iter().all(|v| finite_nonneg(*v)),
            "invalid SAR map",
        )?;
        self.sar.copy_from_slice(values);
        Ok(())
    }
    pub fn set_uniform_sar(&mut self, value: f64) -> Result<()> {
        require(finite_nonneg(value), "invalid SAR")?;
        self.sar.fill(value);
        Ok(())
    }
    pub fn set_volumetric_heating(&mut self, values: &[f64]) -> Result<()> {
        require(
            values.len() == self.heating.len() && values.iter().all(|v| finite_nonneg(*v)),
            "invalid volumetric heating map",
        )?;
        self.heating.copy_from_slice(values);
        Ok(())
    }
    pub fn set_uniform_volumetric_heating(&mut self, value: f64) -> Result<()> {
        require(finite_nonneg(value), "invalid volumetric heating")?;
        self.heating.fill(value);
        Ok(())
    }
    pub fn config(&self) -> &PennesBioheat2DConfig {
        &self.config
    }
    pub fn temperature_k(&self) -> &[f64] {
        &self.temperature
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn linear_result(&self) -> Option<IterativeSolverResult> {
        self.linear
    }
    pub fn step(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.step_once()?;
        }
        Ok(())
    }
    fn step_once(&mut self) -> Result<()> {
        let c = self.config;
        let (nx, ny) = (c.nx, c.ny);
        let periodic = c.boundary == BioheatBoundary::Periodic;
        let dx = c.width_m / (if periodic { nx } else { nx - 1 }) as f64;
        let dy = c.height_m / (if periodic { ny } else { ny - 1 }) as f64;
        let kx = c.thermal_conductivity_w_per_m_k / (dx * dx);
        let ky = c.thermal_conductivity_w_per_m_k / (dy * dy);
        let transient = c.tissue_density_kg_per_m3 * c.tissue_specific_heat_j_per_kg_k / c.dt_s;
        let perfusion =
            c.blood_density_kg_per_m3 * c.blood_specific_heat_j_per_kg_k * c.blood_perfusion_per_s;
        let diagonal = transient + perfusion + 2.0 * (kx + ky);
        let index = |x: usize, y: usize| y * nx + x;
        let boundary = |x: usize, y: usize| x == 0 || y == 0 || x + 1 == nx || y + 1 == ny;
        let tb = c.fixed_boundary_temperature_k;
        let mut rhs = vec![0.0; nx * ny];
        let mut next = self.temperature.clone();
        for y in 0..ny {
            for x in 0..nx {
                let i = index(x, y);
                if !periodic && boundary(x, y) {
                    rhs[i] = tb;
                    next[i] = tb;
                    continue;
                }
                rhs[i] = transient * self.temperature[i]
                    + perfusion * c.blood_temperature_k
                    + c.metabolic_heat_w_per_m3
                    + c.tissue_density_kg_per_m3 * self.sar[i]
                    + self.heating[i];
                if !periodic {
                    if x == 1 {
                        rhs[i] += kx * tb;
                    }
                    if x + 2 == nx {
                        rhs[i] += kx * tb;
                    }
                    if y == 1 {
                        rhs[i] += ky * tb;
                    }
                    if y + 2 == ny {
                        rhs[i] += ky * tb;
                    }
                }
            }
        }
        let apply = |input: &[f64], out: &mut [f64]| -> Result<()> {
            for y in 0..ny {
                for x in 0..nx {
                    let i = index(x, y);
                    if !periodic && boundary(x, y) {
                        out[i] = input[i];
                        continue;
                    }
                    let mut v = diagonal * input[i];
                    if periodic {
                        let (xl, xr) = ((x + nx - 1) % nx, (x + 1) % nx);
                        let (yd, yu) = ((y + ny - 1) % ny, (y + 1) % ny);
                        v -= kx * (input[index(xl, y)] + input[index(xr, y)])
                            + ky * (input[index(x, yd)] + input[index(x, yu)]);
                    } else {
                        if x > 1 {
                            v -= kx * input[index(x - 1, y)];
                        }
                        if x + 2 < nx {
                            v -= kx * input[index(x + 1, y)];
                        }
                        if y > 1 {
                            v -= ky * input[index(x, y - 1)];
                        }
                        if y + 2 < ny {
                            v -= ky * input[index(x, y + 1)];
                        }
                    }
                    out[i] = v;
                }
            }
            Ok(())
        };
        let result = preconditioned_conjugate_gradient(
            &rhs,
            &mut next,
            apply,
            &IdentityPreconditioner,
            LinearConfig {
                max_iterations: c.max_iterations,
                relative_tolerance: c.relative_tolerance,
                restart: 30,
            },
        )?;
        if !result.converged {
            return Err(SolverError::NumericalFailure(
                "Pennes bioheat linear solve failed".into(),
            ));
        }
        self.linear = Some(result);
        self.temperature = next;
        self.time_s += c.dt_s;
        Ok(())
    }
}
