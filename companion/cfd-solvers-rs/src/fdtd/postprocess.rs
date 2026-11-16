//! Near-to-far-field transform, local and mass-averaged SAR, multiconductor line stepping and
//! EMC helpers (graded axes, thin-wire deposition, Huygens sources, impedance sheets, transfer
//! interpolation), ported from upstream `cfd/solvers/fdtd/{postprocess, multiconductor_tline,
//! emc_sources}` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use crate::{core::Vec3, require, Complex64 as C, Result};
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexVec3 {
    pub x: C,
    pub y: C,
    pub z: C,
}

impl Default for ComplexVec3 {
    fn default() -> Self {
        let z = C::new(0.0, 0.0);
        Self { x: z, y: z, z }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NearFieldSurfaceSample {
    pub position: Vec3,
    pub normal: Vec3,
    pub area: f64,
    pub electric: ComplexVec3,
    pub magnetic: ComplexVec3,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FarFieldSample {
    pub electric: ComplexVec3,
    pub radiation_intensity: f64,
}

fn normed(a: Vec3) -> Result<Vec3> {
    let n = (a.x * a.x + a.y * a.y + a.z * a.z).sqrt();
    require(n > 0.0 && n.is_finite(), "zero far-field direction")?;
    Ok(Vec3::new(a.x / n, a.y / n, a.z / n))
}
fn cross_pc(a: Vec3, b: ComplexVec3) -> ComplexVec3 {
    ComplexVec3 {
        x: b.z * a.y - b.y * a.z,
        y: b.x * a.z - b.z * a.x,
        z: b.y * a.x - b.x * a.y,
    }
}
fn cross_cp(a: ComplexVec3, b: Vec3) -> ComplexVec3 {
    ComplexVec3 {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x,
    }
}
fn add(a: ComplexVec3, b: ComplexVec3) -> ComplexVec3 {
    ComplexVec3 {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
}
fn scale(a: ComplexVec3, s: C) -> ComplexVec3 {
    ComplexVec3 {
        x: a.x * s,
        y: a.y * s,
        z: a.z * s,
    }
}
fn abs2(a: ComplexVec3) -> f64 {
    a.x.norm_sqr() + a.y.norm_sqr() + a.z.norm_sqr()
}

/// Upstream default wave impedance, ohm.
pub const NF2FF_WAVE_IMPEDANCE: f64 = 376.730313668;

/// Far-zone `E` (without the `exp(-ikr)/r` factor) from equivalent currents `J = n x H`,
/// `M = -n x E` on a closed sampled surface, and radiation intensity `|E|^2 / (2 eta)`.
pub fn nf2ff(
    surface: &[NearFieldSurfaceSample],
    direction: Vec3,
    frequency_hz: f64,
    wave_impedance: f64,
) -> Result<FarFieldSample> {
    require(
        !surface.is_empty() && frequency_hz > 0.0 && wave_impedance > 0.0,
        "invalid NF2FF inputs",
    )?;
    let r = normed(direction)?;
    let k = 2.0 * PI * frequency_hz / 299792458.0;
    let eta = C::new(wave_impedance, 0.0);
    let mut sum = ComplexVec3::default();
    for q in surface {
        require(q.area > 0.0, "invalid NF2FF surface area")?;
        let n = normed(q.normal)?;
        let j = cross_pc(n, q.magnetic);
        let m = scale(cross_pc(n, q.electric), C::new(-1.0, 0.0));
        let transverse_j = cross_pc(r, cross_cp(j, r));
        let term = add(scale(transverse_j, eta), cross_cp(m, r));
        let phase = k * (r.x * q.position.x + r.y * q.position.y + r.z * q.position.z);
        sum = add(sum, scale(term, C::new(0.0, phase).exp() * q.area));
    }
    sum = scale(sum, C::new(0.0, k / (4.0 * PI)));
    Ok(FarFieldSample {
        electric: sum,
        radiation_intensity: abs2(sum) / (2.0 * wave_impedance),
    })
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct SarGrid3D {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub ex: Vec<f64>,
    pub ey: Vec<f64>,
    pub ez: Vec<f64>,
    pub conductivity: Vec<f64>,
    pub density: Vec<f64>,
}

impl SarGrid3D {
    pub fn validate(&self) -> Result<()> {
        let n = self
            .nx
            .checked_mul(self.ny)
            .and_then(|v| v.checked_mul(self.nz))
            .unwrap_or(0);
        require(
            n > 0
                && self.dx > 0.0
                && self.dy > 0.0
                && self.dz > 0.0
                && [
                    &self.ex,
                    &self.ey,
                    &self.ez,
                    &self.conductivity,
                    &self.density,
                ]
                .iter()
                .all(|v| v.len() == n),
            "invalid SAR grid",
        )?;
        require(
            self.conductivity
                .iter()
                .zip(&self.density)
                .all(|(s, d)| *s >= 0.0 && *d > 0.0),
            "invalid SAR material",
        )
    }
    fn id(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.ny + j) * self.nx + i
    }
}

/// `sigma |E|^2 / rho` per cell (peak-field convention as upstream: no 1/2 factor).
pub fn local_sar(g: &SarGrid3D) -> Result<Vec<f64>> {
    g.validate()?;
    Ok((0..g.ex.len())
        .map(|i| {
            g.conductivity[i] * (g.ex[i] * g.ex[i] + g.ey[i] * g.ey[i] + g.ez[i] * g.ez[i])
                / g.density[i]
        })
        .collect())
}

/// Mass-weighted SAR over the smallest cubic neighbourhood (growing by one cell per side, clipped
/// to the grid) whose mass reaches `target_mass_kg`.
pub fn mass_averaged_sar(g: &SarGrid3D, target_mass_kg: f64) -> Result<Vec<f64>> {
    g.validate()?;
    require(target_mass_kg > 0.0, "SAR target mass must be positive")?;
    let sar = local_sar(g)?;
    let mut out = vec![0.0; sar.len()];
    let vol = g.dx * g.dy * g.dz;
    let maxr = g.nx.max(g.ny).max(g.nz);
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let (mut mass, mut weighted) = (0.0, 0.0);
                let mut r = 0;
                while r < maxr && mass < target_mass_kg {
                    (mass, weighted) = (0.0, 0.0);
                    for kk in k.saturating_sub(r)..=(k + r).min(g.nz - 1) {
                        for jj in j.saturating_sub(r)..=(j + r).min(g.ny - 1) {
                            for ii in i.saturating_sub(r)..=(i + r).min(g.nx - 1) {
                                let n = g.id(ii, jj, kk);
                                let m = g.density[n] * vol;
                                mass += m;
                                weighted += m * sar[n];
                            }
                        }
                    }
                    r += 1;
                }
                out[g.id(i, j, k)] = weighted / mass.max(1e-30);
            }
        }
    }
    Ok(out)
}

pub fn sar_1g(g: &SarGrid3D) -> Result<Vec<f64>> {
    mass_averaged_sar(g, 0.001)
}

pub fn sar_10g(g: &SarGrid3D) -> Result<Vec<f64>> {
    mass_averaged_sar(g, 0.010)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MulticonductorTlineConfig {
    pub conductors: usize,
    /// Diagonal per-unit-length baselines.
    pub inductance_h: Vec<f64>,
    pub capacitance_f: Vec<f64>,
    pub dt_s: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MulticonductorTlineState {
    pub voltage: Vec<f64>,
    pub current: Vec<f64>,
}

/// Symplectic-Euler step of decoupled LC conductors driven by `source_voltage`. Inputs are fully
/// validated before the state is touched.
pub fn advance_multiconductor_tline(
    config: &MulticonductorTlineConfig,
    state: &mut MulticonductorTlineState,
    source_voltage: &[f64],
) -> Result<()> {
    let n = config.conductors;
    require(
        n > 0
            && config.inductance_h.len() == n
            && config.capacitance_f.len() == n
            && state.voltage.len() == n
            && state.current.len() == n
            && source_voltage.len() == n
            && config.dt_s > 0.0,
        "invalid multiconductor line controls",
    )?;
    require(
        (0..n).all(|i| {
            config.inductance_h[i] > 0.0
                && config.capacitance_f[i] > 0.0
                && source_voltage[i].is_finite()
        }),
        "invalid multiconductor line state",
    )?;
    for i in 0..n {
        let current_old = state.current[i];
        state.current[i] +=
            config.dt_s * (state.voltage[i] - source_voltage[i]) / config.inductance_h[i];
        state.voltage[i] -= config.dt_s * current_old / config.capacitance_f[i];
    }
    Ok(())
}

/// Geometric node coordinates `0..=length` whose widths grow by `growth_ratio`.
pub fn make_graded_axis(
    length_m: f64,
    cells: usize,
    first_cell_width_m: f64,
    growth_ratio: f64,
) -> Result<Vec<f64>> {
    require(
        length_m > 0.0
            && cells > 0
            && first_cell_width_m > 0.0
            && growth_ratio > 0.0
            && length_m.is_finite()
            && first_cell_width_m.is_finite()
            && growth_ratio.is_finite(),
        "invalid graded-axis controls",
    )?;
    let sum = if (growth_ratio - 1.0).abs() < 1e-12 {
        cells as f64
    } else {
        (growth_ratio.powf(cells as f64) - 1.0) / (growth_ratio - 1.0)
    };
    let scale = length_m / (first_cell_width_m * sum);
    let mut x = vec![0.0; cells + 1];
    for i in 0..cells {
        let width = first_cell_width_m * scale * growth_ratio.powf(i as f64);
        x[i + 1] = x[i] + width;
    }
    x[cells] = length_m;
    require(x.iter().all(|v| v.is_finite()), "graded axis overflowed")?;
    Ok(x)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SubcellWeight {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub weight: f64,
}

/// Trilinear deposition weights of a point onto the eight surrounding nodes (zero weights
/// dropped; the point is clamped into the grid).
#[allow(clippy::too_many_arguments)]
pub fn deposit_thin_wire_subcell(
    x_m: f64,
    y_m: f64,
    z_m: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    dx_m: f64,
    dy_m: f64,
    dz_m: f64,
) -> Result<Vec<SubcellWeight>> {
    require(
        nx > 0
            && ny > 0
            && nz > 0
            && dx_m > 0.0
            && dy_m > 0.0
            && dz_m > 0.0
            && x_m.is_finite()
            && y_m.is_finite()
            && z_m.is_finite(),
        "invalid thin-wire subcell controls",
    )?;
    let axis = |c: f64, h: f64, n: usize| {
        let q = (c / h).clamp(0.0, n as f64 - 1.0);
        let lo = q.floor() as usize;
        [lo, (lo + 1).min(n - 1)]
    };
    let (ix, iy, iz) = (
        axis(x_m, dx_m, nx),
        axis(y_m, dy_m, ny),
        axis(z_m, dz_m, nz),
    );
    let fx = (x_m / dx_m - ix[0] as f64).clamp(0.0, 1.0);
    let fy = (y_m / dy_m - iy[0] as f64).clamp(0.0, 1.0);
    let fz = (z_m / dz_m - iz[0] as f64).clamp(0.0, 1.0);
    let mut out = Vec::new();
    for a in 0..2 {
        for b in 0..2 {
            for c in 0..2 {
                let w = (if a == 1 { fx } else { 1.0 - fx })
                    * (if b == 1 { fy } else { 1.0 - fy })
                    * (if c == 1 { fz } else { 1.0 - fz });
                if w > 0.0 {
                    out.push(SubcellWeight {
                        i: ix[a],
                        j: iy[b],
                        k: iz[c],
                        weight: w,
                    });
                }
            }
        }
    }
    Ok(out)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HuygensEquivalentSource {
    /// `n x H`.
    pub electric_surface_current: [f64; 3],
    /// `-n x E`.
    pub magnetic_surface_current: [f64; 3],
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn huygens_equivalent_source(
    normal: [f64; 3],
    electric_field: [f64; 3],
    magnetic_field: [f64; 3],
) -> Result<HuygensEquivalentSource> {
    let norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    require(norm > 0.0 && norm.is_finite(), "invalid Huygens normal")?;
    require(
        electric_field.iter().all(|v| v.is_finite()),
        "invalid Huygens electric field",
    )?;
    require(
        magnetic_field.iter().all(|v| v.is_finite()),
        "invalid Huygens magnetic field",
    )?;
    let n = normal.map(|v| v / norm);
    Ok(HuygensEquivalentSource {
        electric_surface_current: cross3(n, magnetic_field),
        magnetic_surface_current: cross3(n, electric_field).map(|v| -v),
    })
}

/// `2 Z2 / (Z1 + Z2 + Zs)` transmission through a series impedance sheet.
pub fn impedance_sheet_transmission(
    sheet_impedance_ohm: C,
    medium1_impedance_ohm: C,
    medium2_impedance_ohm: C,
) -> Result<C> {
    let d = medium1_impedance_ohm + medium2_impedance_ohm + sheet_impedance_ohm;
    require(d.norm() >= 1e-30, "singular impedance-sheet transmission")?;
    Ok(medium2_impedance_ohm * 2.0 / d)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EmcTransferSample {
    pub frequency_hz: f64,
    pub gain: C,
}

/// Linear interpolation in frequency, clamped at both ends.
pub fn interpolate_emc_transfer(samples: &[EmcTransferSample], frequency_hz: f64) -> Result<C> {
    let (Some(first), Some(last)) = (samples.first(), samples.last()) else {
        return Err(crate::SolverError::InvalidInput(
            "invalid EMC transfer request".into(),
        ));
    };
    require(
        frequency_hz >= 0.0 && frequency_hz.is_finite(),
        "invalid EMC transfer request",
    )?;
    if frequency_hz <= first.frequency_hz {
        return Ok(first.gain);
    }
    if frequency_hz >= last.frequency_hz {
        return Ok(last.gain);
    }
    for w in samples.windows(2) {
        require(
            w[1].frequency_hz > w[0].frequency_hz,
            "EMC transfer frequencies must be strictly increasing",
        )?;
        if frequency_hz <= w[1].frequency_hz {
            let t = (frequency_hz - w[0].frequency_hz) / (w[1].frequency_hz - w[0].frequency_hz);
            return Ok(w[0].gain + (w[1].gain - w[0].gain) * t);
        }
    }
    Ok(last.gain)
}
