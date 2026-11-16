//! Ideal-gas mixture thermodynamics, dimensioned fields with time-level history, rigid mesh
//! motion with ALE face fluxes and conservative nearest remapping, and Morton cell reordering,
//! ported from upstream `cfd/fvm/{combustion, field, mesh_motion, reordering}` (MIT, cfd_solvers
//! contributors).
use super::{PolyCell, PolyMesh};
use crate::{chemistry::GAS_CONSTANT, core::Vec3, require, Result, SolverError};
use std::ops::{Add, Sub};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IdealGasSpeciesThermo {
    pub molar_mass: f64,
    pub cp_mass: f64,
    pub formation_enthalpy_mass: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CombustionThermoState {
    pub temperature: f64,
    pub pressure: f64,
    pub density: f64,
    pub mixture_molar_mass: f64,
    pub gas_constant: f64,
    pub cp: f64,
    pub cv: f64,
    pub gamma: f64,
    pub specific_enthalpy: f64,
}

fn check_mixture(y: &[f64], s: &[IdealGasSpeciesThermo]) -> Result<()> {
    require(!y.is_empty() && y.len() == s.len(), "mixture size mismatch")?;
    let mut sum = 0.0;
    for (yi, si) in y.iter().zip(s) {
        require(
            *yi >= 0.0 && yi.is_finite() && si.molar_mass > 0.0 && si.cp_mass > 0.0,
            "invalid mixture thermo",
        )?;
        sum += yi;
    }
    require((sum - 1.0).abs() <= 1e-10, "mass fractions must sum to one")
}

/// Constant-cp ideal-gas mixture: `h = sum Y_i hf_i + cp T`.
pub fn ideal_gas_mixture_state(
    temperature: f64,
    pressure: f64,
    mass_fractions: &[f64],
    species: &[IdealGasSpeciesThermo],
) -> Result<CombustionThermoState> {
    check_mixture(mass_fractions, species)?;
    require(
        temperature > 0.0 && pressure > 0.0,
        "invalid thermodynamic state",
    )?;
    let (mut inv_m, mut cp, mut hf) = (0.0, 0.0, 0.0);
    for (y, s) in mass_fractions.iter().zip(species) {
        inv_m += y / s.molar_mass;
        cp += y * s.cp_mass;
        hf += y * s.formation_enthalpy_mass;
    }
    let m = 1.0 / inv_m;
    let r = GAS_CONSTANT / m;
    let cv = cp - r;
    require(cv > 0.0, "mixture cp must exceed gas constant")?;
    Ok(CombustionThermoState {
        temperature,
        pressure,
        density: pressure / (r * temperature),
        mixture_molar_mass: m,
        gas_constant: r,
        cp,
        cv,
        gamma: cp / cv,
        specific_enthalpy: hf + cp * temperature,
    })
}

pub fn temperature_from_mixture_enthalpy(
    specific_enthalpy: f64,
    mass_fractions: &[f64],
    species: &[IdealGasSpeciesThermo],
) -> Result<f64> {
    check_mixture(mass_fractions, species)?;
    let (mut cp, mut hf) = (0.0, 0.0);
    for (y, s) in mass_fractions.iter().zip(species) {
        cp += y * s.cp_mass;
        hf += y * s.formation_enthalpy_mass;
    }
    let t = (specific_enthalpy - hf) / cp;
    require(
        t > 0.0 && t.is_finite(),
        "mixture enthalpy implies invalid temperature",
    )?;
    Ok(t)
}

/// SI base-unit exponents (length, mass, time, temperature, current, amount, luminosity).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Dimensions {
    pub exponent: [i32; 7],
}

impl Dimensions {
    pub const DIMLESS: Self = Self { exponent: [0; 7] };
    pub const LENGTH: Self = Self {
        exponent: [1, 0, 0, 0, 0, 0, 0],
    };
    pub const MASS: Self = Self {
        exponent: [0, 1, 0, 0, 0, 0, 0],
    };
    pub const TIME: Self = Self {
        exponent: [0, 0, 1, 0, 0, 0, 0],
    };
    pub const TEMPERATURE: Self = Self {
        exponent: [0, 0, 0, 1, 0, 0, 0],
    };
}

/// Dimensions of a product.
impl Add for Dimensions {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            exponent: std::array::from_fn(|i| self.exponent[i] + b.exponent[i]),
        }
    }
}

/// Dimensions of a quotient.
impl Sub for Dimensions {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            exponent: std::array::from_fn(|i| self.exponent[i] - b.exponent[i]),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FieldLocation {
    #[default]
    Cell,
    Face,
    Point,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DimensionedField<T: Clone> {
    name: String,
    dimensions: Dimensions,
    location: FieldLocation,
    values: Vec<T>,
}

impl<T: Clone> DimensionedField<T> {
    pub fn new(
        name: &str,
        dimensions: Dimensions,
        location: FieldLocation,
        size: usize,
        initial: T,
    ) -> Self {
        Self {
            name: name.into(),
            dimensions,
            location,
            values: vec![initial; size],
        }
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn dimensions(&self) -> Dimensions {
        self.dimensions
    }
    pub fn location(&self) -> FieldLocation {
        self.location
    }
    pub fn values(&self) -> &[T] {
        &self.values
    }
    pub fn values_mut(&mut self) -> &mut Vec<T> {
        &mut self.values
    }
    pub fn size(&self) -> usize {
        self.values.len()
    }
}

/// Newest-first time levels, capped at `levels`.
#[derive(Clone, Debug, PartialEq)]
pub struct FieldHistory<T: Clone> {
    levels: usize,
    states: Vec<Vec<T>>,
}

impl<T: Clone> FieldHistory<T> {
    /// Upstream default 2 levels.
    pub fn new(levels: usize) -> Result<Self> {
        require(levels > 0, "field history requires at least one level")?;
        Ok(Self {
            levels,
            states: Vec::new(),
        })
    }
    pub fn push(&mut self, state: Vec<T>) -> Result<()> {
        require(
            self.states.first().is_none_or(|s| s.len() == state.len()),
            "field history size mismatch",
        )?;
        self.states.insert(0, state);
        self.states.truncate(self.levels);
        Ok(())
    }
    pub fn current(&self) -> Result<&[T]> {
        self.states
            .first()
            .map(Vec::as_slice)
            .ok_or_else(|| SolverError::InvalidInput("empty field history".into()))
    }
    /// Upstream default level 1.
    pub fn old(&self, level: usize) -> Result<&[T]> {
        self.states
            .get(level)
            .map(Vec::as_slice)
            .ok_or_else(|| SolverError::InvalidInput("field history level unavailable".into()))
    }
    pub fn stored_levels(&self) -> usize {
        self.states.len()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bdf2Coefficients {
    pub current: f64,
    pub old: f64,
    pub older: f64,
}

impl Default for Bdf2Coefficients {
    fn default() -> Self {
        Self {
            current: 1.5,
            old: -2.0,
            older: 0.5,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrankNicolsonWeights {
    pub implicit: f64,
    pub explicit_old: f64,
}

impl Default for CrankNicolsonWeights {
    fn default() -> Self {
        Self {
            implicit: 0.5,
            explicit_old: 0.5,
        }
    }
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

fn rotate(p: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let m = dot(axis, axis).sqrt();
    if m < 1e-30 {
        return p;
    }
    let a = Vec3::new(axis.x / m, axis.y / m, axis.z / m);
    let (c, s) = (angle.cos(), angle.sin());
    let cross = Vec3::new(
        a.y * p.z - a.z * p.y,
        a.z * p.x - a.x * p.z,
        a.x * p.y - a.y * p.x,
    );
    let k = dot(a, p) * (1.0 - c);
    Vec3::new(
        p.x * c + cross.x * s + a.x * k,
        p.y * c + cross.y * s + a.y * k,
        p.z * c + cross.z * s + a.z * k,
    )
}

/// Rodrigues rotation about the origin followed by translation of cell/face centres and areas.
pub fn rigidly_moved_mesh(
    mesh: &PolyMesh,
    translation: Vec3,
    rotation_axis: Vec3,
    angle_radians: f64,
) -> Result<PolyMesh> {
    let shift = |p: Vec3| {
        let r = rotate(p, rotation_axis, angle_radians);
        Vec3::new(
            r.x + translation.x,
            r.y + translation.y,
            r.z + translation.z,
        )
    };
    let cells = mesh
        .cells()
        .iter()
        .map(|c| PolyCell {
            center: shift(c.center),
            volume: c.volume,
        })
        .collect();
    let faces = mesh
        .faces()
        .iter()
        .map(|f| {
            let mut f = *f;
            f.center = shift(f.center);
            f.area = rotate(f.area, rotation_axis, angle_radians);
            f
        })
        .collect();
    PolyMesh::new(cells, faces, mesh.patches().to_vec())
}

/// Face-centre swept flux `((x_new - x_old) / dt) . S_old`.
pub fn ale_mesh_flux(old_mesh: &PolyMesh, new_mesh: &PolyMesh, dt: f64) -> Result<Vec<f64>> {
    require(
        old_mesh.faces().len() == new_mesh.faces().len() && dt > 0.0,
        "ALE mesh mismatch",
    )?;
    Ok(old_mesh
        .faces()
        .iter()
        .zip(new_mesh.faces())
        .map(|(a, b)| {
            let v = Vec3::new(
                (b.center.x - a.center.x) / dt,
                (b.center.y - a.center.y) / dt,
                (b.center.z - a.center.z) / dt,
            );
            dot(v, a.area)
        })
        .collect())
}

/// Nearest-centre injection rescaled so the volume integral is conserved.
pub fn conservative_nearest_remap(
    old_mesh: &PolyMesh,
    old_field: &[f64],
    new_mesh: &PolyMesh,
) -> Result<Vec<f64>> {
    require(
        old_field.len() == old_mesh.cells().len(),
        "remap field size",
    )?;
    let mut out: Vec<f64> = new_mesh
        .cells()
        .iter()
        .map(|nc| {
            let mut best = 0;
            let mut d2 = 1e300;
            for (i, oc) in old_mesh.cells().iter().enumerate() {
                let d = Vec3::new(
                    oc.center.x - nc.center.x,
                    oc.center.y - nc.center.y,
                    oc.center.z - nc.center.z,
                );
                let q = dot(d, d);
                if q < d2 {
                    d2 = q;
                    best = i;
                }
            }
            old_field[best]
        })
        .collect();
    let old_integral: f64 = old_field
        .iter()
        .zip(old_mesh.cells())
        .map(|(v, c)| v * c.volume)
        .fold(0.0, |a, b| a + b);
    let new_integral: f64 = out
        .iter()
        .zip(new_mesh.cells())
        .map(|(v, c)| v * c.volume)
        .fold(0.0, |a, b| a + b);
    if new_integral.abs() > 1e-30 {
        let s = old_integral / new_integral;
        out.iter_mut().for_each(|v| *v *= s);
    }
    Ok(out)
}

fn part1(x: u32) -> u64 {
    let mut v = u64::from(x & 0x1f_ffff);
    v = (v | (v << 32)) & 0x1f_0000_0000_ffff;
    v = (v | (v << 16)) & 0x1f_0000_ff00_00ff;
    v = (v | (v << 8)) & 0x100f_00f0_0f00_f00f;
    v = (v | (v << 4)) & 0x10c3_0c30_c30c_30c3;
    v = (v | (v << 2)) & 0x1249_2492_4924_9249;
    v
}

fn morton(x: u32, y: u32, z: u32) -> u64 {
    part1(x) | (part1(y) << 1) | (part1(z) << 2)
}

/// Cells sorted by the 63-bit Morton code of their 21-bit-quantised centres.
pub fn morton_cell_order(mesh: &PolyMesh) -> Vec<usize> {
    let cells = mesh.cells();
    let Some(first) = cells.first() else {
        return Vec::new();
    };
    let (mut lo, mut hi) = (first.center, first.center);
    for c in cells {
        lo = Vec3::new(
            lo.x.min(c.center.x),
            lo.y.min(c.center.y),
            lo.z.min(c.center.z),
        );
        hi = Vec3::new(
            hi.x.max(c.center.x),
            hi.y.max(c.center.y),
            hi.z.max(c.center.z),
        );
    }
    let quant = |v: f64, a: f64, b: f64| {
        if b > a {
            (((v - a) / (b - a)).clamp(0.0, 1.0) * 2_097_151.0) as u32
        } else {
            0
        }
    };
    let key = |c: &PolyCell| {
        morton(
            quant(c.center.x, lo.x, hi.x),
            quant(c.center.y, lo.y, hi.y),
            quant(c.center.z, lo.z, hi.z),
        )
    };
    let mut order: Vec<usize> = (0..cells.len()).collect();
    order.sort_by_key(|i| key(&cells[*i]));
    order
}

/// New mesh whose cell `k` is old cell `order[k]`; faces are renumbered, not reordered.
pub fn reorder_cells(mesh: &PolyMesh, order: &[usize]) -> Result<PolyMesh> {
    let n = mesh.cells().len();
    require(order.len() == n, "cell reorder permutation size mismatch")?;
    let mut inverse: Vec<Option<usize>> = vec![None; n];
    let mut cells = Vec::with_capacity(n);
    for (next, &old) in order.iter().enumerate() {
        require(
            old < n && inverse[old].is_none(),
            "cell reorder must be a permutation",
        )?;
        inverse[old] = Some(next);
        cells.push(mesh.cells()[old]);
    }
    let mut faces = mesh.faces().to_vec();
    for f in &mut faces {
        f.owner = inverse[f.owner]
            .ok_or_else(|| SolverError::InvalidInput("face owner outside cell reorder".into()))?;
        if let Some(nb) = f.neighbour {
            f.neighbour = Some(inverse[nb].ok_or_else(|| {
                SolverError::InvalidInput("face neighbour outside cell reorder".into())
            })?);
        }
    }
    PolyMesh::new(cells, faces, mesh.patches().to_vec())
}

pub fn reorder_cell_field<T: Clone>(field: &[T], order: &[usize]) -> Result<Vec<T>> {
    require(field.len() == order.len(), "cell reorder size mismatch")?;
    order
        .iter()
        .map(|i| {
            field
                .get(*i)
                .cloned()
                .ok_or_else(|| SolverError::InvalidInput("cell reorder index out of range".into()))
        })
        .collect()
}
