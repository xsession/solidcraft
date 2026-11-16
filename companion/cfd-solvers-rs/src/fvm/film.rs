//! Lubrication thin-film transport and barotropic compressible VOF, ported from upstream
//! `cfd/solvers/fvm/{thin_film, compressible_vof}` (MIT, cfd_solvers contributors).
use super::{
    gauss_divergence_vector, interpolate_scalar_to_faces, least_squares_gradient_scalar,
    lubrication_thin_film_flux, particles::finite_vec, Face, FaceInterpolationScheme, PolyMesh,
};
use crate::{core::Vec3, require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinFilmConfig {
    pub dt: f64,
    pub density: f64,
    pub viscosity: f64,
    pub surface_tension: f64,
    pub body_acceleration: Vec3,
    pub minimum_thickness: f64,
}

impl Default for ThinFilmConfig {
    fn default() -> Self {
        Self {
            dt: 1e-4,
            density: 1000.0,
            viscosity: 1e-3,
            surface_tension: 0.0,
            body_acceleration: Vec3::default(),
            minimum_thickness: 1e-9,
        }
    }
}

/// Explicit lubrication film on a unit-depth mesh: capillary pressure from interface curvature
/// plus a body-force potential drive the `h^3` flux. Boundaries are impermeable unless a
/// boundary pressure is set.
#[derive(Clone, Debug)]
pub struct ThinFilmTransport {
    mesh: PolyMesh,
    config: ThinFilmConfig,
    thickness: Vec<f64>,
    pressure: Vec<f64>,
    face_flux: Vec<f64>,
    boundary_pressure: Vec<f64>,
}

fn face_distance(mesh: &PolyMesh, face: &Face) -> Result<f64> {
    let delta = match face.neighbour {
        None => face.center - mesh.cells()[face.owner].center,
        Some(nb) => mesh.cells()[nb].center - mesh.cells()[face.owner].center,
    };
    let area = face.area.norm();
    require(area > 0.0, "thin-film face has zero area")?;
    let projected = delta.dot(face.area / area).abs();
    let distance = if projected > 1e-14 {
        projected
    } else {
        delta.norm()
    };
    require(distance > 0.0, "thin-film face has zero center distance")?;
    Ok(distance)
}

impl ThinFilmTransport {
    pub fn new(mesh: PolyMesh, config: ThinFilmConfig) -> Result<Self> {
        let c = config;
        require(
            c.dt > 0.0
                && c.density > 0.0
                && c.viscosity > 0.0
                && c.surface_tension >= 0.0
                && c.minimum_thickness >= 0.0
                && [
                    c.dt,
                    c.density,
                    c.viscosity,
                    c.surface_tension,
                    c.minimum_thickness,
                ]
                .iter()
                .all(|v| v.is_finite()),
            "invalid thin-film controls",
        )?;
        require(
            finite_vec(c.body_acceleration),
            "invalid thin-film body acceleration",
        )?;
        let (n, f) = (mesh.cell_count(), mesh.face_count());
        Ok(Self {
            mesh,
            config,
            thickness: vec![c.minimum_thickness; n],
            pressure: vec![0.0; n],
            face_flux: vec![0.0; f],
            boundary_pressure: Vec::new(),
        })
    }
    pub fn initialize(&mut self, thickness: f64) -> Result<()> {
        require(
            thickness.is_finite() && thickness >= 0.0,
            "invalid thin-film thickness",
        )?;
        self.thickness
            .fill(self.config.minimum_thickness.max(thickness));
        self.pressure.fill(0.0);
        self.face_flux.fill(0.0);
        Ok(())
    }
    pub fn initialize_profile(&mut self, f: impl Fn(Vec3) -> f64) {
        let min = self.config.minimum_thickness;
        for (h, c) in self.thickness.iter_mut().zip(self.mesh.cells()) {
            *h = min.max(f(c.center));
        }
        self.pressure.fill(0.0);
        self.face_flux.fill(0.0);
    }
    /// Fixed pressure per face (only boundary entries are read).
    pub fn set_boundary_pressure(&mut self, pressure: Vec<f64>) -> Result<()> {
        require(
            pressure.len() == self.mesh.face_count(),
            "thin-film boundary pressure size mismatch",
        )?;
        require(
            pressure.iter().all(|p| p.is_finite()),
            "invalid thin-film boundary pressure",
        )?;
        self.boundary_pressure = pressure;
        Ok(())
    }
    pub fn clear_boundary_pressure(&mut self) {
        self.boundary_pressure.clear();
    }
    fn update_pressure(&mut self) -> Result<()> {
        let gradient = least_squares_gradient_scalar(&self.mesh, &self.thickness, &[])?;
        let normal: Vec<Vec3> = gradient
            .iter()
            .map(|g| {
                let m = g.norm();
                if m > 1e-14 {
                    *g / m
                } else {
                    Vec3::default()
                }
            })
            .collect();
        let curvature = gauss_divergence_vector(&self.mesh, &normal, &[])?;
        for (p, k) in self.pressure.iter_mut().zip(&curvature) {
            *p = -self.config.surface_tension * k;
        }
        Ok(())
    }
    fn update_face_flux(&mut self) -> Result<()> {
        self.face_flux.fill(0.0);
        let c = self.config;
        for (fi, face) in self.mesh.faces().iter().enumerate() {
            if face.neighbour.is_none() && self.boundary_pressure.is_empty() {
                continue;
            }
            let area = face.area.norm();
            let distance = face_distance(&self.mesh, face)?;
            let po = self.pressure[face.owner];
            let (pn, h) = match face.neighbour {
                None => (self.boundary_pressure[fi], self.thickness[face.owner]),
                Some(nb) => (
                    self.pressure[nb],
                    0.5 * (self.thickness[face.owner] + self.thickness[nb]),
                ),
            };
            let normal = face.area / area;
            let gradient = (pn - po) / distance - c.density * c.body_acceleration.dot(normal);
            self.face_flux[fi] = lubrication_thin_film_flux(h, gradient, c.viscosity, 0.0)? * area;
        }
        Ok(())
    }
    pub fn step(&mut self) -> Result<()> {
        self.update_pressure()?;
        self.update_face_flux()?;
        let mut next = self.thickness.clone();
        for (cell, out) in next.iter_mut().enumerate() {
            let mut net = 0.0;
            for &fi in &self.mesh.cell_faces()[cell] {
                net += if self.mesh.faces()[fi].owner == cell {
                    self.face_flux[fi]
                } else {
                    -self.face_flux[fi]
                };
            }
            *out = self
                .config
                .minimum_thickness
                .max(self.thickness[cell] - self.config.dt * net / self.mesh.cells()[cell].volume);
        }
        self.thickness = next;
        Ok(())
    }
    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn thickness(&self) -> &[f64] {
        &self.thickness
    }
    pub fn pressure(&self) -> &[f64] {
        &self.pressure
    }
    pub fn face_flux(&self) -> &[f64] {
        &self.face_flux
    }
    pub fn inventory(&self) -> f64 {
        self.thickness
            .iter()
            .zip(self.mesh.cells())
            .map(|(h, c)| h * c.volume)
            .sum()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompressibleVofConfig {
    pub dt: f64,
    pub liquid_density: f64,
    pub gas_density: f64,
    pub liquid_bulk_modulus: f64,
    pub gas_gamma: f64,
    pub reference_pressure: f64,
    pub minimum_pressure: f64,
    pub minimum_density: f64,
}

impl Default for CompressibleVofConfig {
    fn default() -> Self {
        Self {
            dt: 1e-4,
            liquid_density: 1000.0,
            gas_density: 1.2,
            liquid_bulk_modulus: 2.2e9,
            gas_gamma: 1.4,
            reference_pressure: 1e5,
            minimum_pressure: 1.0,
            minimum_density: 1e-12,
        }
    }
}

/// Isothermal/barotropic two-phase VOF: liquid fraction and mixture mass are transported
/// conservatively with upwind face values; pressure is recovered by bisection on the mixture EOS
/// (linear liquid compressibility, isentropic gas).
#[derive(Clone, Debug)]
pub struct CompressibleVofTransport {
    mesh: PolyMesh,
    config: CompressibleVofConfig,
    liquid_fraction: Vec<f64>,
    mass: Vec<f64>,
    density: Vec<f64>,
    pressure: Vec<f64>,
    liquid_density: Vec<f64>,
    gas_density: Vec<f64>,
    face_flux: Vec<f64>,
    time: f64,
    steps: usize,
}

impl CompressibleVofTransport {
    pub fn new(mesh: PolyMesh, config: CompressibleVofConfig) -> Result<Self> {
        let c = config;
        require(
            [
                c.dt,
                c.liquid_density,
                c.gas_density,
                c.liquid_bulk_modulus,
                c.gas_gamma,
                c.reference_pressure,
                c.minimum_pressure,
                c.minimum_density,
            ]
            .iter()
            .all(|v| v.is_finite())
                && c.dt > 0.0
                && c.liquid_density > 0.0
                && c.gas_density > 0.0
                && c.liquid_bulk_modulus > 0.0
                && c.gas_gamma > 1.0
                && c.reference_pressure > 0.0
                && c.minimum_pressure > 0.0
                && c.minimum_density > 0.0
                && c.minimum_pressure <= c.reference_pressure,
            "invalid compressible VOF controls",
        )?;
        let (n, f) = (mesh.cell_count(), mesh.face_count());
        let mut s = Self {
            mesh,
            config,
            liquid_fraction: vec![0.0; n],
            mass: vec![0.0; n],
            density: vec![c.gas_density; n],
            pressure: vec![c.reference_pressure; n],
            liquid_density: vec![c.liquid_density; n],
            gas_density: vec![c.gas_density; n],
            face_flux: vec![0.0; f],
            time: 0.0,
            steps: 0,
        };
        s.initialize_uniform(0.0, c.reference_pressure)?;
        Ok(s)
    }
    fn phase_densities(&self, pressure: f64) -> (f64, f64) {
        let c = self.config;
        (
            c.liquid_density * (1.0 + (pressure - c.reference_pressure) / c.liquid_bulk_modulus),
            c.gas_density * (pressure / c.reference_pressure).powf(1.0 / c.gas_gamma),
        )
    }
    fn mixture_density(&self, alpha: f64, pressure: f64) -> f64 {
        let (l, g) = self.phase_densities(pressure);
        alpha * l + (1.0 - alpha) * g
    }
    fn pressure_for_density(&self, alpha: f64, density: f64) -> Result<f64> {
        let c = self.config;
        require(
            (0.0..=1.0).contains(&alpha) && density >= c.minimum_density && density.is_finite(),
            "invalid compressible VOF mixture state",
        )?;
        let mut lower = c.minimum_pressure;
        if density <= self.mixture_density(alpha, lower) {
            return Ok(lower);
        }
        let mut upper = (c.reference_pressure * 2.0).max(lower * 2.0);
        for _ in 0..80 {
            if self.mixture_density(alpha, upper) >= density {
                break;
            }
            upper *= 2.0;
        }
        require(
            self.mixture_density(alpha, upper) >= density && upper.is_finite(),
            "compressible VOF pressure bracket failed",
        )?;
        for _ in 0..80 {
            let middle = 0.5 * (lower + upper);
            if self.mixture_density(alpha, middle) < density {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        Ok(0.5 * (lower + upper))
    }
    pub fn initialize_uniform(&mut self, liquid_fraction: f64, pressure: f64) -> Result<()> {
        let n = self.mesh.cell_count();
        self.initialize(&vec![liquid_fraction; n], &vec![pressure; n])
    }
    /// `pressure == 0` selects the reference pressure, as upstream.
    pub fn initialize_profile(
        &mut self,
        fraction: impl Fn(Vec3) -> f64,
        pressure: f64,
    ) -> Result<()> {
        let p = if pressure == 0.0 {
            self.config.reference_pressure
        } else {
            pressure
        };
        let a: Vec<f64> = self
            .mesh
            .cells()
            .iter()
            .map(|c| fraction(c.center))
            .collect();
        let n = a.len();
        self.initialize(&a, &vec![p; n])
    }
    pub fn initialize(&mut self, liquid_fraction: &[f64], pressure: &[f64]) -> Result<()> {
        let n = self.mesh.cell_count();
        require(
            liquid_fraction.len() == n && pressure.len() == n,
            "compressible VOF initial field size mismatch",
        )?;
        require(
            liquid_fraction
                .iter()
                .all(|a| a.is_finite() && (0.0..=1.0).contains(a)),
            "compressible VOF fraction must lie in [0,1]",
        )?;
        require(
            pressure
                .iter()
                .all(|p| p.is_finite() && *p >= self.config.minimum_pressure),
            "compressible VOF pressure must be positive",
        )?;
        for c in 0..n {
            self.mass[c] =
                self.mixture_density(liquid_fraction[c], pressure[c]) * self.mesh.cells()[c].volume;
        }
        self.liquid_fraction.copy_from_slice(liquid_fraction);
        self.pressure.copy_from_slice(pressure);
        self.face_flux.fill(0.0);
        self.update_properties()?;
        self.time = 0.0;
        self.steps = 0;
        Ok(())
    }
    pub fn set_face_flux(&mut self, volumetric_flux: Vec<f64>) -> Result<()> {
        require(
            volumetric_flux.len() == self.mesh.face_count(),
            "compressible VOF face-flux size mismatch",
        )?;
        require(
            volumetric_flux.iter().all(|f| f.is_finite()),
            "compressible VOF face flux must be finite",
        )?;
        self.face_flux = volumetric_flux;
        Ok(())
    }
    fn update_properties(&mut self) -> Result<()> {
        for c in 0..self.mesh.cell_count() {
            let rho = self.mass[c] / self.mesh.cells()[c].volume;
            require(
                rho >= self.config.minimum_density && rho.is_finite(),
                "compressible VOF update produced non-positive density",
            )?;
            self.density[c] = rho;
            let p = self.pressure_for_density(self.liquid_fraction[c], rho)?;
            self.pressure[c] = p;
            (self.liquid_density[c], self.gas_density[c]) = self.phase_densities(p);
        }
        Ok(())
    }
    pub fn step(&mut self) -> Result<()> {
        let face_fraction = interpolate_scalar_to_faces(
            &self.mesh,
            &self.liquid_fraction,
            FaceInterpolationScheme::Upwind,
            &self.face_flux,
            &[],
        )?;
        let face_density = interpolate_scalar_to_faces(
            &self.mesh,
            &self.density,
            FaceInterpolationScheme::Upwind,
            &self.face_flux,
            &[],
        )?;
        let mut next_fraction = self.liquid_fraction.clone();
        let mut next_mass = self.mass.clone();
        let dt = self.config.dt;
        for cell in 0..self.mesh.cell_count() {
            let (mut vflux, mut mflux) = (0.0, 0.0);
            for &fi in &self.mesh.cell_faces()[cell] {
                let o = if self.mesh.faces()[fi].owner == cell {
                    1.0
                } else {
                    -1.0
                };
                vflux += o * self.face_flux[fi] * face_fraction[fi];
                mflux += o * self.face_flux[fi] * face_density[fi];
            }
            let volume = self.mesh.cells()[cell].volume;
            next_fraction[cell] =
                (self.liquid_fraction[cell] - dt * vflux / volume).clamp(0.0, 1.0);
            next_mass[cell] = self.mass[cell] - dt * mflux;
            require(
                next_mass[cell] >= self.config.minimum_density * volume
                    && next_mass[cell].is_finite(),
                "compressible VOF update produced non-positive mass",
            )?;
        }
        self.liquid_fraction = next_fraction;
        self.mass = next_mass;
        self.update_properties()?;
        self.time += dt;
        self.steps += 1;
        Ok(())
    }
    pub fn mesh(&self) -> &PolyMesh {
        &self.mesh
    }
    pub fn config(&self) -> &CompressibleVofConfig {
        &self.config
    }
    pub fn liquid_fraction(&self) -> &[f64] {
        &self.liquid_fraction
    }
    pub fn density(&self) -> &[f64] {
        &self.density
    }
    pub fn pressure(&self) -> &[f64] {
        &self.pressure
    }
    pub fn liquid_density(&self) -> &[f64] {
        &self.liquid_density
    }
    pub fn gas_density(&self) -> &[f64] {
        &self.gas_density
    }
    pub fn face_flux(&self) -> &[f64] {
        &self.face_flux
    }
    pub fn liquid_volume(&self) -> f64 {
        self.liquid_fraction
            .iter()
            .zip(self.mesh.cells())
            .map(|(a, c)| a * c.volume)
            .sum()
    }
    pub fn total_mass(&self) -> f64 {
        self.mass.iter().sum()
    }
    pub fn minimum_pressure(&self) -> f64 {
        self.pressure.iter().copied().fold(f64::INFINITY, f64::min)
    }
    pub fn maximum_pressure(&self) -> f64 {
        self.pressure
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn steps(&self) -> usize {
        self.steps
    }
}
