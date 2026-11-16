//! 3-D particle bookkeeping for domain-decomposed PIC: cell sorting, guard-halo classification,
//! domain grids, particle migration plans/messages, scalar and field guard-cell exchange, rank
//! topologies, in-memory and serialized transports. Ported from upstream
//! `cfd/particle/electromagnetic` (MIT, cfd_solvers contributors).
#![allow(clippy::needless_range_loop)]
use super::pusher::*;
use crate::{require, Result, SolverError};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PicParticle3D {
    pub position_m: Vec3,
    pub velocity_m_per_s: Vec3,
    pub charge_c: f64,
    pub mass_kg: f64,
    pub weight: f64,
}

impl Default for PicParticle3D {
    fn default() -> Self {
        Self {
            position_m: Vec3::default(),
            velocity_m_per_s: Vec3::default(),
            charge_c: 1.0,
            mass_kg: 1.0,
            weight: 1.0,
        }
    }
}

impl PicParticle3D {
    pub(crate) fn proxy(&self) -> ChargedParticle {
        ChargedParticle {
            position_m: self.position_m,
            velocity_m_per_s: self.velocity_m_per_s,
            charge_c: self.charge_c,
            mass_kg: self.mass_kg,
            weight: self.weight,
        }
    }
    pub(crate) fn valid(&self) -> bool {
        self.mass_kg > 0.0 && self.charge_c.is_finite() && self.weight.is_finite()
    }
}

impl CollisionParticle for PicParticle3D {
    fn velocity(&self) -> Vec3 {
        self.velocity_m_per_s
    }
    fn set_velocity(&mut self, v: Vec3) {
        self.velocity_m_per_s = v;
    }
    fn mass(&self) -> f64 {
        self.mass_kg
    }
    fn weight(&self) -> f64 {
        self.weight
    }
    fn as_charged(&self) -> ChargedParticle {
        self.proxy()
    }
}

pub(crate) fn fi3(ix: usize, iy: usize, iz: usize, nx: usize, ny: usize) -> usize {
    (iz * ny + iy) * nx + ix
}

pub fn particle_cell_index_3d(
    p: Vec3,
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<usize> {
    require(
        nx >= 1 && ny >= 1 && nz >= 1 && lx > 0.0 && ly > 0.0 && lz > 0.0,
        "invalid 3-D particle cell grid",
    )?;
    let idx = |v: f64, l: f64, n: usize| {
        ((wrap_periodic(v, l) / (l / n as f64)).floor() as usize).min(n - 1)
    };
    Ok(fi3(
        idx(p.x, lx, nx),
        idx(p.y, ly, ny),
        idx(p.z, lz, nz),
        nx,
        ny,
    ))
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParticleCellSort3DResult {
    pub sorted_particles: Vec<PicParticle3D>,
    pub original_indices: Vec<usize>,
    pub cell_offsets: Vec<usize>,
    pub cell_counts: Vec<usize>,
    pub occupied_cells: usize,
}

/// Stable counting sort of particles by periodic cell index.
pub fn sort_particles_by_cell_3d(
    particles: &[PicParticle3D],
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<ParticleCellSort3DResult> {
    require(
        nx >= 1 && ny >= 1 && nz >= 1 && lx > 0.0 && ly > 0.0 && lz > 0.0,
        "invalid 3-D particle sorting grid",
    )?;
    let cells = nx * ny * nz;
    let mut r = ParticleCellSort3DResult {
        cell_counts: vec![0; cells],
        cell_offsets: vec![0; cells + 1],
        ..Default::default()
    };
    let mut particle_cells = Vec::with_capacity(particles.len());
    for p in particles {
        require(p.valid(), "invalid particle for 3-D sorting")?;
        let c = particle_cell_index_3d(p.position_m, nx, ny, nz, lx, ly, lz)?;
        particle_cells.push(c);
        r.cell_counts[c] += 1;
    }
    for c in 0..cells {
        r.cell_offsets[c + 1] = r.cell_offsets[c] + r.cell_counts[c];
        if r.cell_counts[c] > 0 {
            r.occupied_cells += 1;
        }
    }
    r.sorted_particles = vec![PicParticle3D::default(); particles.len()];
    r.original_indices = vec![0; particles.len()];
    let mut write = r.cell_offsets.clone();
    for (i, p) in particles.iter().enumerate() {
        let d = write[particle_cells[i]];
        write[particle_cells[i]] += 1;
        r.sorted_particles[d] = *p;
        r.original_indices[d] = i;
    }
    Ok(r)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleGuardHalo3DConfig {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub length_x_m: f64,
    pub length_y_m: f64,
    pub length_z_m: f64,
    pub guard_cells: usize,
    pub periodic: bool,
}

impl Default for ParticleGuardHalo3DConfig {
    fn default() -> Self {
        Self {
            nx: 1,
            ny: 1,
            nz: 1,
            length_x_m: 1.0,
            length_y_m: 1.0,
            length_z_m: 1.0,
            guard_cells: 1,
            periodic: true,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ParticleGuardHalo3DResult {
    pub x_minus: Vec<usize>,
    pub x_plus: Vec<usize>,
    pub y_minus: Vec<usize>,
    pub y_plus: Vec<usize>,
    pub z_minus: Vec<usize>,
    pub z_plus: Vec<usize>,
    pub interior: Vec<usize>,
    pub outside_domain: Vec<usize>,
}

/// Particles within `guard_cells` cells of each face (a particle may be listed on several faces).
pub fn classify_particle_guard_halos_3d(
    particles: &[PicParticle3D],
    c: &ParticleGuardHalo3DConfig,
) -> Result<ParticleGuardHalo3DResult> {
    require(
        c.nx >= 1
            && c.ny >= 1
            && c.nz >= 1
            && c.guard_cells >= 1
            && c.length_x_m > 0.0
            && c.length_y_m > 0.0
            && c.length_z_m > 0.0,
        "invalid 3-D particle guard halo config",
    )?;
    let mut r = ParticleGuardHalo3DResult::default();
    let g = c.guard_cells as f64;
    let tx = c.length_x_m.min(g * c.length_x_m / c.nx as f64);
    let ty = c.length_y_m.min(g * c.length_y_m / c.ny as f64);
    let tz = c.length_z_m.min(g * c.length_z_m / c.nz as f64);
    for (i, particle) in particles.iter().enumerate() {
        let mut p = particle.position_m;
        let outside = p.x < 0.0
            || p.x >= c.length_x_m
            || p.y < 0.0
            || p.y >= c.length_y_m
            || p.z < 0.0
            || p.z >= c.length_z_m;
        if outside && !c.periodic {
            r.outside_domain.push(i);
            continue;
        }
        if c.periodic {
            p = Vec3::new(
                wrap_periodic(p.x, c.length_x_m),
                wrap_periodic(p.y, c.length_y_m),
                wrap_periodic(p.z, c.length_z_m),
            );
        }
        let mut halo = false;
        for (cond, list) in [
            (p.x < tx, &mut r.x_minus),
            (p.x >= c.length_x_m - tx, &mut r.x_plus),
            (p.y < ty, &mut r.y_minus),
            (p.y >= c.length_y_m - ty, &mut r.y_plus),
            (p.z < tz, &mut r.z_minus),
            (p.z >= c.length_z_m - tz, &mut r.z_plus),
        ] {
            if cond {
                list.push(i);
                halo = true;
            }
        }
        if !halo {
            r.interior.push(i);
        }
    }
    Ok(r)
}

/// Marker for "no domain".
pub const INVALID_PIC_DOMAIN_ID_3D: usize = usize::MAX;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PicDomainGrid3DConfig {
    pub global_nx: usize,
    pub global_ny: usize,
    pub global_nz: usize,
    pub domains_x: usize,
    pub domains_y: usize,
    pub domains_z: usize,
    pub length_x_m: f64,
    pub length_y_m: f64,
    pub length_z_m: f64,
    pub periodic: bool,
}

impl Default for PicDomainGrid3DConfig {
    fn default() -> Self {
        Self {
            global_nx: 1,
            global_ny: 1,
            global_nz: 1,
            domains_x: 1,
            domains_y: 1,
            domains_z: 1,
            length_x_m: 1.0,
            length_y_m: 1.0,
            length_z_m: 1.0,
            periodic: true,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PicDomain3D {
    pub domain_id: usize,
    pub ix: usize,
    pub iy: usize,
    pub iz: usize,
    pub first_x: usize,
    pub first_y: usize,
    pub first_z: usize,
    pub cells_x: usize,
    pub cells_y: usize,
    pub cells_z: usize,
    pub minimum_m: Vec3,
    pub maximum_m: Vec3,
}

pub(crate) fn validate_domain_grid(c: &PicDomainGrid3DConfig) -> Result<()> {
    require(
        c.global_nx >= 1
            && c.global_ny >= 1
            && c.global_nz >= 1
            && c.domains_x >= 1
            && c.domains_y >= 1
            && c.domains_z >= 1
            && c.domains_x <= c.global_nx
            && c.domains_y <= c.global_ny
            && c.domains_z <= c.global_nz
            && c.length_x_m > 0.0
            && c.length_y_m > 0.0
            && c.length_z_m > 0.0,
        "invalid 3-D PIC domain grid",
    )
}

fn block_start(b: usize, blocks: usize, cells: usize) -> usize {
    (b * cells) / blocks
}

pub(crate) fn wrap_cell(cell: i64, cells: usize) -> usize {
    cell.rem_euclid(cells as i64) as usize
}

pub(crate) fn locate_domain_from_cell(
    gx: usize,
    gy: usize,
    gz: usize,
    domains: &[PicDomain3D],
) -> usize {
    domains
        .iter()
        .find(|d| {
            gx >= d.first_x
                && gx < d.first_x + d.cells_x
                && gy >= d.first_y
                && gy < d.first_y + d.cells_y
                && gz >= d.first_z
                && gz < d.first_z + d.cells_z
        })
        .map_or(INVALID_PIC_DOMAIN_ID_3D, |d| d.domain_id)
}

/// Balanced block decomposition, domain id `(iz·Dy + iy)·Dx + ix`.
pub fn make_pic_domain_grid_3d(c: &PicDomainGrid3DConfig) -> Result<Vec<PicDomain3D>> {
    validate_domain_grid(c)?;
    let (dx, dy, dz) = (
        c.length_x_m / c.global_nx as f64,
        c.length_y_m / c.global_ny as f64,
        c.length_z_m / c.global_nz as f64,
    );
    let mut out = Vec::with_capacity(c.domains_x * c.domains_y * c.domains_z);
    for iz in 0..c.domains_z {
        for iy in 0..c.domains_y {
            for ix in 0..c.domains_x {
                let fx = block_start(ix, c.domains_x, c.global_nx);
                let fy = block_start(iy, c.domains_y, c.global_ny);
                let fz = block_start(iz, c.domains_z, c.global_nz);
                let cx = block_start(ix + 1, c.domains_x, c.global_nx) - fx;
                let cy = block_start(iy + 1, c.domains_y, c.global_ny) - fy;
                let cz = block_start(iz + 1, c.domains_z, c.global_nz) - fz;
                out.push(PicDomain3D {
                    domain_id: fi3(ix, iy, iz, c.domains_x, c.domains_y),
                    ix,
                    iy,
                    iz,
                    first_x: fx,
                    first_y: fy,
                    first_z: fz,
                    cells_x: cx,
                    cells_y: cy,
                    cells_z: cz,
                    minimum_m: Vec3::new(fx as f64 * dx, fy as f64 * dy, fz as f64 * dz),
                    maximum_m: Vec3::new(
                        (fx + cx) as f64 * dx,
                        (fy + cy) as f64 * dy,
                        (fz + cz) as f64 * dz,
                    ),
                });
            }
        }
    }
    Ok(out)
}

fn wrap3(p: Vec3, c: &PicDomainGrid3DConfig) -> Vec3 {
    Vec3::new(
        wrap_periodic(p.x, c.length_x_m),
        wrap_periodic(p.y, c.length_y_m),
        wrap_periodic(p.z, c.length_z_m),
    )
}

/// Owning domain of a position ([`INVALID_PIC_DOMAIN_ID_3D`] outside a non-periodic domain).
pub fn locate_pic_domain_3d(p: Vec3, c: &PicDomainGrid3DConfig) -> Result<usize> {
    validate_domain_grid(c)?;
    let p = if c.periodic {
        wrap3(p, c)
    } else {
        if p.x < 0.0
            || p.x >= c.length_x_m
            || p.y < 0.0
            || p.y >= c.length_y_m
            || p.z < 0.0
            || p.z >= c.length_z_m
        {
            return Ok(INVALID_PIC_DOMAIN_ID_3D);
        }
        p
    };
    let cell = |v: f64, l: f64, n: usize| ((v / (l / n as f64)).floor() as usize).min(n - 1);
    let domains = make_pic_domain_grid_3d(c)?;
    Ok(locate_domain_from_cell(
        cell(p.x, c.length_x_m, c.global_nx),
        cell(p.y, c.length_y_m, c.global_ny),
        cell(p.z, c.length_z_m, c.global_nz),
        &domains,
    ))
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParticleDomainMigration3D {
    pub particles_by_domain: Vec<Vec<PicParticle3D>>,
    pub source_indices_by_domain: Vec<Vec<usize>>,
    pub destination_domain_by_particle: Vec<usize>,
    pub outside_indices: Vec<usize>,
}

pub fn plan_particle_domain_migration_3d(
    particles: &[PicParticle3D],
    c: &PicDomainGrid3DConfig,
) -> Result<ParticleDomainMigration3D> {
    let domains = make_pic_domain_grid_3d(c)?;
    let mut r = ParticleDomainMigration3D {
        particles_by_domain: vec![Vec::new(); domains.len()],
        source_indices_by_domain: vec![Vec::new(); domains.len()],
        destination_domain_by_particle: vec![INVALID_PIC_DOMAIN_ID_3D; particles.len()],
        outside_indices: Vec::new(),
    };
    for (i, particle) in particles.iter().enumerate() {
        let d = locate_pic_domain_3d(particle.position_m, c)?;
        r.destination_domain_by_particle[i] = d;
        if d == INVALID_PIC_DOMAIN_ID_3D {
            r.outside_indices.push(i);
            continue;
        }
        let mut p = *particle;
        if c.periodic {
            p.position_m = wrap3(p.position_m, c);
        }
        r.particles_by_domain[d].push(p);
        r.source_indices_by_domain[d].push(i);
    }
    Ok(r)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScalarGuardedBlock3D {
    pub domain_id: usize,
    pub interior_nx: usize,
    pub interior_ny: usize,
    pub interior_nz: usize,
    pub guard_cells: usize,
    pub values: Vec<f64>,
}

impl ScalarGuardedBlock3D {
    pub fn padded_nx(&self) -> usize {
        self.interior_nx + 2 * self.guard_cells
    }
    pub fn padded_ny(&self) -> usize {
        self.interior_ny + 2 * self.guard_cells
    }
    pub fn padded_nz(&self) -> usize {
        self.interior_nz + 2 * self.guard_cells
    }
}

pub(crate) fn validate_domain_blocks(values: &[Vec<f64>], domains: &[PicDomain3D]) -> Result<()> {
    require(
        values.len() == domains.len(),
        "domain value count must match domains",
    )?;
    for d in domains {
        require(d.domain_id < domains.len(), "domain id out of range")?;
        require(
            values[d.domain_id].len() == d.cells_x * d.cells_y * d.cells_z,
            "scalar domain block size mismatch",
        )?;
    }
    Ok(())
}

/// Global cell owning padded index `(ix,iy,iz)` of `domain`, after periodic wrapping.
pub(crate) fn padded_to_global(
    d: &PicDomain3D,
    ix: usize,
    iy: usize,
    iz: usize,
    g: usize,
    c: &PicDomainGrid3DConfig,
) -> Option<(usize, usize, usize)> {
    let mut gx = d.first_x as i64 + ix as i64 - g as i64;
    let mut gy = d.first_y as i64 + iy as i64 - g as i64;
    let mut gz = d.first_z as i64 + iz as i64 - g as i64;
    if c.periodic {
        gx = wrap_cell(gx, c.global_nx) as i64;
        gy = wrap_cell(gy, c.global_ny) as i64;
        gz = wrap_cell(gz, c.global_nz) as i64;
    }
    if gx < 0
        || gx >= c.global_nx as i64
        || gy < 0
        || gy >= c.global_ny as i64
        || gz < 0
        || gz >= c.global_nz as i64
    {
        return None;
    }
    Some((gx as usize, gy as usize, gz as usize))
}

/// Padded blocks (guard width `g`) filled from owning domains; cells outside a non-periodic
/// domain keep `exterior_value`.
pub fn exchange_scalar_guard_cells_3d(
    values: &[Vec<f64>],
    domains: &[PicDomain3D],
    c: &PicDomainGrid3DConfig,
    g: usize,
    exterior_value: f64,
) -> Result<Vec<ScalarGuardedBlock3D>> {
    validate_domain_grid(c)?;
    require(
        g >= 1 && g <= c.global_nx && g <= c.global_ny && g <= c.global_nz,
        "invalid scalar guard-cell width",
    )?;
    validate_domain_blocks(values, domains)?;
    let mut out = vec![ScalarGuardedBlock3D::default(); domains.len()];
    for d in domains {
        let mut b = ScalarGuardedBlock3D {
            domain_id: d.domain_id,
            interior_nx: d.cells_x,
            interior_ny: d.cells_y,
            interior_nz: d.cells_z,
            guard_cells: g,
            values: Vec::new(),
        };
        let (px, py, pz) = (b.padded_nx(), b.padded_ny(), b.padded_nz());
        b.values = vec![exterior_value; px * py * pz];
        for iz in 0..pz {
            for iy in 0..py {
                for ix in 0..px {
                    let Some((gx, gy, gz)) = padded_to_global(d, ix, iy, iz, g, c) else {
                        continue;
                    };
                    let owner = locate_domain_from_cell(gx, gy, gz, domains);
                    if owner == INVALID_PIC_DOMAIN_ID_3D {
                        continue;
                    }
                    let od = &domains[owner];
                    b.values[fi3(ix, iy, iz, px, py)] = values[owner][fi3(
                        gx - od.first_x,
                        gy - od.first_y,
                        gz - od.first_z,
                        od.cells_x,
                        od.cells_y,
                    )];
                }
            }
        }
        out[d.domain_id] = b;
    }
    Ok(out)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ElectromagneticFieldBlocks3D {
    pub electric_x_by_domain: Vec<Vec<f64>>,
    pub electric_y_by_domain: Vec<Vec<f64>>,
    pub electric_z_by_domain: Vec<Vec<f64>>,
    pub magnetic_x_by_domain: Vec<Vec<f64>>,
    pub magnetic_y_by_domain: Vec<Vec<f64>>,
    pub magnetic_z_by_domain: Vec<Vec<f64>>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ElectromagneticGuardedFieldBlocks3D {
    pub electric_x: Vec<ScalarGuardedBlock3D>,
    pub electric_y: Vec<ScalarGuardedBlock3D>,
    pub electric_z: Vec<ScalarGuardedBlock3D>,
    pub magnetic_x: Vec<ScalarGuardedBlock3D>,
    pub magnetic_y: Vec<ScalarGuardedBlock3D>,
    pub magnetic_z: Vec<ScalarGuardedBlock3D>,
}

/// Guard-cell exchange of all six field components.
pub fn exchange_electromagnetic_field_guard_cells_3d(
    f: &ElectromagneticFieldBlocks3D,
    domains: &[PicDomain3D],
    c: &PicDomainGrid3DConfig,
    g: usize,
    exterior_value: f64,
) -> Result<ElectromagneticGuardedFieldBlocks3D> {
    let comps = [
        (&f.electric_x_by_domain, "Ex"),
        (&f.electric_y_by_domain, "Ey"),
        (&f.electric_z_by_domain, "Ez"),
        (&f.magnetic_x_by_domain, "Bx"),
        (&f.magnetic_y_by_domain, "By"),
        (&f.magnetic_z_by_domain, "Bz"),
    ];
    for (comp, name) in comps {
        require(
            comp.len() == domains.len(),
            format!("{name} field domain count mismatch"),
        )?;
        for d in domains {
            require(
                comp.get(d.domain_id)
                    .is_some_and(|v| v.len() == d.cells_x * d.cells_y * d.cells_z),
                format!("{name} field block size mismatch"),
            )?;
        }
    }
    Ok(ElectromagneticGuardedFieldBlocks3D {
        electric_x: exchange_scalar_guard_cells_3d(
            &f.electric_x_by_domain,
            domains,
            c,
            g,
            exterior_value,
        )?,
        electric_y: exchange_scalar_guard_cells_3d(
            &f.electric_y_by_domain,
            domains,
            c,
            g,
            exterior_value,
        )?,
        electric_z: exchange_scalar_guard_cells_3d(
            &f.electric_z_by_domain,
            domains,
            c,
            g,
            exterior_value,
        )?,
        magnetic_x: exchange_scalar_guard_cells_3d(
            &f.magnetic_x_by_domain,
            domains,
            c,
            g,
            exterior_value,
        )?,
        magnetic_y: exchange_scalar_guard_cells_3d(
            &f.magnetic_y_by_domain,
            domains,
            c,
            g,
            exterior_value,
        )?,
        magnetic_z: exchange_scalar_guard_cells_3d(
            &f.magnetic_z_by_domain,
            domains,
            c,
            g,
            exterior_value,
        )?,
    })
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ParticleMigrationRecord3D {
    pub source_domain_id: usize,
    pub source_index: usize,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParticleMigrationMessage3D {
    pub source_domain_id: usize,
    pub destination_domain_id: usize,
    pub particles: Vec<PicParticle3D>,
    pub source_indices: Vec<usize>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PackedParticleMigration3D {
    pub retained_by_domain: Vec<Vec<PicParticle3D>>,
    pub retained_source_indices_by_domain: Vec<Vec<usize>>,
    pub messages: Vec<ParticleMigrationMessage3D>,
    pub outside_particles: Vec<ParticleMigrationRecord3D>,
}

fn validate_domain_span(domains: &[PicDomain3D], c: &PicDomainGrid3DConfig) -> Result<()> {
    validate_domain_grid(c)?;
    require(
        domains.len() == c.domains_x * c.domains_y * c.domains_z,
        "PIC domain span size does not match grid",
    )?;
    require(
        domains.iter().all(|d| d.domain_id < domains.len()),
        "PIC domain id out of range",
    )
}

/// Split each source domain's particles into retained ones and per-(source, destination)
/// migration messages (first-seen order), recording particles that left the domain.
pub fn pack_particle_migration_messages_3d(
    by_source: &[Vec<PicParticle3D>],
    domains: &[PicDomain3D],
    c: &PicDomainGrid3DConfig,
) -> Result<PackedParticleMigration3D> {
    validate_domain_span(domains, c)?;
    require(
        by_source.len() == domains.len(),
        "particle source bucket count must match PIC domains",
    )?;
    let mut r = PackedParticleMigration3D {
        retained_by_domain: vec![Vec::new(); domains.len()],
        retained_source_indices_by_domain: vec![Vec::new(); domains.len()],
        ..Default::default()
    };
    for d in domains {
        let src = d.domain_id;
        for (si, particle) in by_source[src].iter().enumerate() {
            let dst = locate_pic_domain_3d(particle.position_m, c)?;
            if dst == INVALID_PIC_DOMAIN_ID_3D {
                r.outside_particles.push(ParticleMigrationRecord3D {
                    source_domain_id: src,
                    source_index: si,
                });
                continue;
            }
            let mut p = *particle;
            if c.periodic {
                p.position_m = wrap3(p.position_m, c);
            }
            if dst == src {
                r.retained_by_domain[src].push(p);
                r.retained_source_indices_by_domain[src].push(si);
            } else {
                let pos = match r
                    .messages
                    .iter()
                    .position(|m| m.source_domain_id == src && m.destination_domain_id == dst)
                {
                    Some(i) => i,
                    None => {
                        r.messages.push(ParticleMigrationMessage3D {
                            source_domain_id: src,
                            destination_domain_id: dst,
                            ..Default::default()
                        });
                        r.messages.len() - 1
                    }
                };
                r.messages[pos].particles.push(p);
                r.messages[pos].source_indices.push(si);
            }
        }
    }
    Ok(r)
}

pub fn apply_particle_migration_messages_3d(
    mut retained: Vec<Vec<PicParticle3D>>,
    messages: &[ParticleMigrationMessage3D],
) -> Result<Vec<Vec<PicParticle3D>>> {
    for m in messages {
        require(
            m.destination_domain_id < retained.len(),
            "particle migration message destination out of range",
        )?;
        require(
            m.particles.len() == m.source_indices.len(),
            "particle migration message source-index size mismatch",
        )?;
        retained[m.destination_domain_id].extend_from_slice(&m.particles);
    }
    Ok(retained)
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ScalarGuardCellValue3D {
    pub padded_x: usize,
    pub padded_y: usize,
    pub padded_z: usize,
    pub value: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScalarGuardCellMessage3D {
    pub source_domain_id: usize,
    pub destination_domain_id: usize,
    pub offset_x: i32,
    pub offset_y: i32,
    pub offset_z: i32,
    pub values: Vec<ScalarGuardCellValue3D>,
}

fn offset_sign(padded: usize, g: usize, interior: usize) -> i32 {
    if padded < g {
        -1
    } else if padded >= g + interior {
        1
    } else {
        0
    }
}

/// One message per (source, destination, face/edge/corner offset) carrying guard values.
pub fn pack_scalar_guard_cell_messages_3d(
    values: &[Vec<f64>],
    domains: &[PicDomain3D],
    c: &PicDomainGrid3DConfig,
    g: usize,
) -> Result<Vec<ScalarGuardCellMessage3D>> {
    validate_domain_span(domains, c)?;
    require(
        g >= 1 && g <= c.global_nx && g <= c.global_ny && g <= c.global_nz,
        "invalid scalar guard-cell width",
    )?;
    validate_domain_blocks(values, domains)?;
    let mut messages: Vec<ScalarGuardCellMessage3D> = Vec::new();
    for d in domains {
        let (px, py, pz) = (d.cells_x + 2 * g, d.cells_y + 2 * g, d.cells_z + 2 * g);
        for iz in 0..pz {
            for iy in 0..py {
                for ix in 0..px {
                    let interior = ix >= g
                        && ix < g + d.cells_x
                        && iy >= g
                        && iy < g + d.cells_y
                        && iz >= g
                        && iz < g + d.cells_z;
                    if interior {
                        continue;
                    }
                    let Some((gx, gy, gz)) = padded_to_global(d, ix, iy, iz, g, c) else {
                        continue;
                    };
                    let src = locate_domain_from_cell(gx, gy, gz, domains);
                    if src == INVALID_PIC_DOMAIN_ID_3D {
                        continue;
                    }
                    let sd = &domains[src];
                    let value = values[src][fi3(
                        gx - sd.first_x,
                        gy - sd.first_y,
                        gz - sd.first_z,
                        sd.cells_x,
                        sd.cells_y,
                    )];
                    let (ox, oy, oz) = (
                        offset_sign(ix, g, d.cells_x),
                        offset_sign(iy, g, d.cells_y),
                        offset_sign(iz, g, d.cells_z),
                    );
                    let pos = match messages.iter().position(|m| {
                        m.source_domain_id == src
                            && m.destination_domain_id == d.domain_id
                            && m.offset_x == ox
                            && m.offset_y == oy
                            && m.offset_z == oz
                    }) {
                        Some(i) => i,
                        None => {
                            messages.push(ScalarGuardCellMessage3D {
                                source_domain_id: src,
                                destination_domain_id: d.domain_id,
                                offset_x: ox,
                                offset_y: oy,
                                offset_z: oz,
                                values: Vec::new(),
                            });
                            messages.len() - 1
                        }
                    };
                    messages[pos].values.push(ScalarGuardCellValue3D {
                        padded_x: ix,
                        padded_y: iy,
                        padded_z: iz,
                        value,
                    });
                }
            }
        }
    }
    Ok(messages)
}

/// Rebuild padded blocks from interiors plus received guard messages.
pub fn apply_scalar_guard_cell_messages_3d(
    values: &[Vec<f64>],
    domains: &[PicDomain3D],
    g: usize,
    exterior_value: f64,
    messages: &[ScalarGuardCellMessage3D],
) -> Result<Vec<ScalarGuardedBlock3D>> {
    require(g >= 1, "invalid scalar guard-cell width")?;
    validate_domain_blocks(values, domains)?;
    let mut blocks = vec![ScalarGuardedBlock3D::default(); domains.len()];
    for d in domains {
        let mut b = ScalarGuardedBlock3D {
            domain_id: d.domain_id,
            interior_nx: d.cells_x,
            interior_ny: d.cells_y,
            interior_nz: d.cells_z,
            guard_cells: g,
            values: Vec::new(),
        };
        let (px, py, pz) = (b.padded_nx(), b.padded_ny(), b.padded_nz());
        b.values = vec![exterior_value; px * py * pz];
        for iz in 0..d.cells_z {
            for iy in 0..d.cells_y {
                for ix in 0..d.cells_x {
                    b.values[fi3(ix + g, iy + g, iz + g, px, py)] =
                        values[d.domain_id][fi3(ix, iy, iz, d.cells_x, d.cells_y)];
                }
            }
        }
        blocks[d.domain_id] = b;
    }
    for m in messages {
        let b = blocks.get_mut(m.destination_domain_id).ok_or_else(|| {
            SolverError::InvalidInput("scalar guard message destination out of range".into())
        })?;
        let (px, py, pz) = (b.padded_nx(), b.padded_ny(), b.padded_nz());
        for e in &m.values {
            require(
                e.padded_x < px && e.padded_y < py && e.padded_z < pz,
                "scalar guard message entry out of range",
            )?;
            b.values[fi3(e.padded_x, e.padded_y, e.padded_z, px, py)] = e.value;
        }
    }
    Ok(blocks)
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum PicTransportPayloadKind3D {
    #[default]
    ParticleMigration,
    ScalarGuardCells,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PicRankTopology3D {
    pub rank_count: usize,
    pub domain_to_rank: Vec<usize>,
}

/// Block-map domains onto a `ranks_x × ranks_y × ranks_z` rank grid.
pub fn make_pic_rank_topology_3d(
    domains: &[PicDomain3D],
    rx: usize,
    ry: usize,
    rz: usize,
) -> Result<PicRankTopology3D> {
    require(
        !domains.is_empty() && rx > 0 && ry > 0 && rz > 0,
        "invalid PIC rank topology dimensions",
    )?;
    let (mut cx, mut cy, mut cz) = (0usize, 0usize, 0usize);
    for d in domains {
        require(
            d.domain_id < domains.len(),
            "PIC domain id out of range for rank topology",
        )?;
        cx = cx.max(d.ix + 1);
        cy = cy.max(d.iy + 1);
        cz = cz.max(d.iz + 1);
    }
    require(
        rx <= cx && ry <= cy && rz <= cz,
        "PIC rank topology cannot exceed domain grid dimensions",
    )?;
    let mut t = PicRankTopology3D {
        rank_count: rx * ry * rz,
        domain_to_rank: vec![0; domains.len()],
    };
    for d in domains {
        let a = (rx - 1).min(d.ix * rx / cx);
        let b = (ry - 1).min(d.iy * ry / cy);
        let e = (rz - 1).min(d.iz * rz / cz);
        t.domain_to_rank[d.domain_id] = (e * ry + b) * rx + a;
    }
    Ok(t)
}

pub fn rank_for_pic_domain_3d(t: &PicRankTopology3D, domain_id: usize) -> Result<usize> {
    let rank = *t.domain_to_rank.get(domain_id).ok_or_else(|| {
        SolverError::InvalidInput("PIC domain id out of range for rank lookup".into())
    })?;
    require(rank < t.rank_count, "PIC topology contains invalid rank id")?;
    Ok(rank)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PicTransportEnvelope3D {
    pub kind: PicTransportPayloadKind3D,
    pub source_domain_id: usize,
    pub destination_domain_id: usize,
    pub source_rank: usize,
    pub destination_rank: usize,
    pub particle_migration: ParticleMigrationMessage3D,
    pub scalar_guard: ScalarGuardCellMessage3D,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PicTransportDiagnostics3D {
    pub envelopes: usize,
    pub particle_messages: usize,
    pub scalar_guard_messages: usize,
    pub particle_payload_count: usize,
    pub scalar_guard_value_count: usize,
    pub same_rank_messages: usize,
    pub remote_rank_messages: usize,
}

/// Per-rank inboxes standing in for a message-passing transport.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InMemoryPicTransport3D {
    rank_count: usize,
    inbox_by_rank: Vec<Vec<PicTransportEnvelope3D>>,
}

impl InMemoryPicTransport3D {
    pub fn new(rank_count: usize) -> Result<Self> {
        require(
            rank_count > 0,
            "PIC in-memory transport requires at least one rank",
        )?;
        Ok(Self {
            rank_count,
            inbox_by_rank: vec![Vec::new(); rank_count],
        })
    }
    /// Posts all envelopes, or none when any rank is out of range.
    pub fn post(&mut self, envelopes: &[PicTransportEnvelope3D]) -> Result<()> {
        for e in envelopes {
            require(
                e.destination_rank < self.rank_count,
                "PIC transport destination rank out of range",
            )?;
            require(
                e.source_rank < self.rank_count,
                "PIC transport source rank out of range",
            )?;
        }
        for e in envelopes {
            self.inbox_by_rank[e.destination_rank].push(e.clone());
        }
        Ok(())
    }
    pub fn receive_rank(&mut self, rank: usize) -> Result<Vec<PicTransportEnvelope3D>> {
        require(
            rank < self.rank_count,
            "PIC transport receive rank out of range",
        )?;
        Ok(std::mem::take(&mut self.inbox_by_rank[rank]))
    }
    pub fn pending_count(&self) -> usize {
        self.inbox_by_rank.iter().map(Vec::len).sum()
    }
    pub fn clear(&mut self) {
        self.inbox_by_rank.iter_mut().for_each(Vec::clear);
    }
}

pub fn make_pic_transport_envelopes_3d(
    particle_messages: &[ParticleMigrationMessage3D],
    scalar_messages: &[ScalarGuardCellMessage3D],
    t: &PicRankTopology3D,
) -> Result<Vec<PicTransportEnvelope3D>> {
    require(t.rank_count > 0, "PIC transport topology has no ranks")?;
    let mut out = Vec::with_capacity(particle_messages.len() + scalar_messages.len());
    for m in particle_messages {
        out.push(PicTransportEnvelope3D {
            kind: PicTransportPayloadKind3D::ParticleMigration,
            source_domain_id: m.source_domain_id,
            destination_domain_id: m.destination_domain_id,
            source_rank: rank_for_pic_domain_3d(t, m.source_domain_id)?,
            destination_rank: rank_for_pic_domain_3d(t, m.destination_domain_id)?,
            particle_migration: m.clone(),
            scalar_guard: ScalarGuardCellMessage3D::default(),
        });
    }
    for m in scalar_messages {
        out.push(PicTransportEnvelope3D {
            kind: PicTransportPayloadKind3D::ScalarGuardCells,
            source_domain_id: m.source_domain_id,
            destination_domain_id: m.destination_domain_id,
            source_rank: rank_for_pic_domain_3d(t, m.source_domain_id)?,
            destination_rank: rank_for_pic_domain_3d(t, m.destination_domain_id)?,
            particle_migration: ParticleMigrationMessage3D::default(),
            scalar_guard: m.clone(),
        });
    }
    Ok(out)
}

pub fn summarize_pic_transport_3d(
    envelopes: &[PicTransportEnvelope3D],
) -> PicTransportDiagnostics3D {
    let mut d = PicTransportDiagnostics3D {
        envelopes: envelopes.len(),
        ..Default::default()
    };
    for e in envelopes {
        if e.source_rank == e.destination_rank {
            d.same_rank_messages += 1;
        } else {
            d.remote_rank_messages += 1;
        }
        if e.kind == PicTransportPayloadKind3D::ParticleMigration {
            d.particle_messages += 1;
            d.particle_payload_count += e.particle_migration.particles.len();
        } else {
            d.scalar_guard_messages += 1;
            d.scalar_guard_value_count += e.scalar_guard.values.len();
        }
    }
    d
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PicTransportedExchange3D {
    pub particles_by_domain: Vec<Vec<PicParticle3D>>,
    pub guarded_blocks: Vec<ScalarGuardedBlock3D>,
    pub diagnostics: PicTransportDiagnostics3D,
}

fn split_delivered(
    envelopes: Vec<PicTransportEnvelope3D>,
) -> (
    Vec<ParticleMigrationMessage3D>,
    Vec<ScalarGuardCellMessage3D>,
) {
    let (mut p, mut s) = (Vec::new(), Vec::new());
    for e in envelopes {
        if e.kind == PicTransportPayloadKind3D::ParticleMigration {
            p.push(e.particle_migration);
        } else {
            s.push(e.scalar_guard);
        }
    }
    (p, s)
}

/// Route messages through an [`InMemoryPicTransport3D`] and apply them.
#[allow(clippy::too_many_arguments)]
pub fn exchange_pic_messages_in_memory_3d(
    retained: Vec<Vec<PicParticle3D>>,
    particle_messages: &[ParticleMigrationMessage3D],
    values: &[Vec<f64>],
    scalar_messages: &[ScalarGuardCellMessage3D],
    domains: &[PicDomain3D],
    g: usize,
    exterior_value: f64,
    t: &PicRankTopology3D,
) -> Result<PicTransportedExchange3D> {
    let envelopes = make_pic_transport_envelopes_3d(particle_messages, scalar_messages, t)?;
    let mut transport = InMemoryPicTransport3D::new(t.rank_count)?;
    transport.post(&envelopes)?;
    let mut delivered = Vec::new();
    for rank in 0..t.rank_count {
        delivered.extend(transport.receive_rank(rank)?);
    }
    let (pm, sm) = split_delivered(delivered);
    Ok(PicTransportedExchange3D {
        particles_by_domain: apply_particle_migration_messages_3d(retained, &pm)?,
        guarded_blocks: apply_scalar_guard_cell_messages_3d(
            values,
            domains,
            g,
            exterior_value,
            &sm,
        )?,
        diagnostics: summarize_pic_transport_3d(&envelopes),
    })
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PicSerializedEnvelope3D {
    pub bytes: Vec<u8>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PicSerializedExchangePlan3D {
    pub rank_count: usize,
    pub envelopes_by_destination_rank: Vec<Vec<PicSerializedEnvelope3D>>,
    pub message_count_by_destination_rank: Vec<usize>,
    pub byte_count_by_destination_rank: Vec<usize>,
}

const PIC_TRANSPORT_MAGIC: u32 = 0x5049_4333;
const PIC_TRANSPORT_VERSION: u32 = 1;

#[derive(Default)]
struct Writer(Vec<u8>);
impl Writer {
    fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn i32(&mut self, v: i32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn size(&mut self, v: usize) {
        self.0.extend_from_slice(&(v as u64).to_le_bytes());
    }
    fn f64(&mut self, v: f64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn vec3(&mut self, v: Vec3) {
        self.f64(v.x);
        self.f64(v.y);
        self.f64(v.z);
    }
    fn particle(&mut self, p: &PicParticle3D) {
        self.vec3(p.position_m);
        self.vec3(p.velocity_m_per_s);
        self.f64(p.charge_c);
        self.f64(p.mass_kg);
        self.f64(p.weight);
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    offset: usize,
}
impl Reader<'_> {
    fn take<const N: usize>(&mut self) -> Result<[u8; N]> {
        let end = self
            .offset
            .checked_add(N)
            .filter(|&e| e <= self.bytes.len());
        let Some(end) = end else {
            return Err(SolverError::InvalidInput(
                "truncated PIC transport envelope".into(),
            ));
        };
        let mut out = [0u8; N];
        out.copy_from_slice(&self.bytes[self.offset..end]);
        self.offset = end;
        Ok(out)
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take()?))
    }
    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take()?))
    }
    fn size(&mut self) -> Result<usize> {
        usize::try_from(u64::from_le_bytes(self.take()?))
            .map_err(|_| SolverError::InvalidInput("PIC transport size overflow".into()))
    }
    fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.take()?))
    }
    fn vec3(&mut self) -> Result<Vec3> {
        Ok(Vec3::new(self.f64()?, self.f64()?, self.f64()?))
    }
    fn particle(&mut self) -> Result<PicParticle3D> {
        Ok(PicParticle3D {
            position_m: self.vec3()?,
            velocity_m_per_s: self.vec3()?,
            charge_c: self.f64()?,
            mass_kg: self.f64()?,
            weight: self.f64()?,
        })
    }
}

/// Little-endian wire format identical to upstream on x86-64 (`PIC3` magic, version 1,
/// 64-bit sizes).
pub fn serialize_pic_transport_envelope_3d(
    e: &PicTransportEnvelope3D,
) -> Result<PicSerializedEnvelope3D> {
    let mut w = Writer::default();
    w.u32(PIC_TRANSPORT_MAGIC);
    w.u32(PIC_TRANSPORT_VERSION);
    w.u32(match e.kind {
        PicTransportPayloadKind3D::ParticleMigration => 0,
        PicTransportPayloadKind3D::ScalarGuardCells => 1,
    });
    w.size(e.source_domain_id);
    w.size(e.destination_domain_id);
    w.size(e.source_rank);
    w.size(e.destination_rank);
    if e.kind == PicTransportPayloadKind3D::ParticleMigration {
        let m = &e.particle_migration;
        require(
            m.particles.len() == m.source_indices.len(),
            "PIC particle message size mismatch before serialization",
        )?;
        require(
            m.source_domain_id == e.source_domain_id
                && m.destination_domain_id == e.destination_domain_id,
            "PIC particle envelope domain mismatch",
        )?;
        w.size(m.source_domain_id);
        w.size(m.destination_domain_id);
        w.size(m.particles.len());
        for (p, &si) in m.particles.iter().zip(&m.source_indices) {
            w.particle(p);
            w.size(si);
        }
    } else {
        let m = &e.scalar_guard;
        require(
            m.source_domain_id == e.source_domain_id
                && m.destination_domain_id == e.destination_domain_id,
            "PIC scalar guard envelope domain mismatch",
        )?;
        w.size(m.source_domain_id);
        w.size(m.destination_domain_id);
        w.i32(m.offset_x);
        w.i32(m.offset_y);
        w.i32(m.offset_z);
        w.size(m.values.len());
        for v in &m.values {
            w.size(v.padded_x);
            w.size(v.padded_y);
            w.size(v.padded_z);
            w.f64(v.value);
        }
    }
    Ok(PicSerializedEnvelope3D { bytes: w.0 })
}

pub fn deserialize_pic_transport_envelope_3d(
    s: &PicSerializedEnvelope3D,
) -> Result<PicTransportEnvelope3D> {
    require(!s.bytes.is_empty(), "empty PIC transport envelope")?;
    let mut r = Reader {
        bytes: &s.bytes,
        offset: 0,
    };
    let magic = r.u32()?;
    let version = r.u32()?;
    require(
        magic == PIC_TRANSPORT_MAGIC && version == PIC_TRANSPORT_VERSION,
        "unsupported PIC transport envelope version",
    )?;
    let kind = match r.u32()? {
        0 => PicTransportPayloadKind3D::ParticleMigration,
        1 => PicTransportPayloadKind3D::ScalarGuardCells,
        _ => {
            return Err(SolverError::InvalidInput(
                "invalid PIC transport payload kind".into(),
            ))
        }
    };
    let mut e = PicTransportEnvelope3D {
        kind,
        source_domain_id: r.size()?,
        destination_domain_id: r.size()?,
        source_rank: r.size()?,
        destination_rank: r.size()?,
        ..Default::default()
    };
    if kind == PicTransportPayloadKind3D::ParticleMigration {
        let mut m = ParticleMigrationMessage3D {
            source_domain_id: r.size()?,
            destination_domain_id: r.size()?,
            ..Default::default()
        };
        let count = r.size()?;
        require(
            count <= 1 << 28,
            "unreasonable PIC particle message payload",
        )?;
        for _ in 0..count {
            m.particles.push(r.particle()?);
            m.source_indices.push(r.size()?);
        }
        require(
            m.source_domain_id == e.source_domain_id
                && m.destination_domain_id == e.destination_domain_id,
            "PIC particle message/envelope domain mismatch after deserialize",
        )?;
        e.particle_migration = m;
    } else {
        let mut m = ScalarGuardCellMessage3D {
            source_domain_id: r.size()?,
            destination_domain_id: r.size()?,
            offset_x: r.i32()?,
            offset_y: r.i32()?,
            offset_z: r.i32()?,
            values: Vec::new(),
        };
        let count = r.size()?;
        require(
            count <= 1 << 28,
            "unreasonable PIC scalar guard message payload",
        )?;
        for _ in 0..count {
            m.values.push(ScalarGuardCellValue3D {
                padded_x: r.size()?,
                padded_y: r.size()?,
                padded_z: r.size()?,
                value: r.f64()?,
            });
        }
        require(
            m.source_domain_id == e.source_domain_id
                && m.destination_domain_id == e.destination_domain_id,
            "PIC scalar guard message/envelope domain mismatch after deserialize",
        )?;
        e.scalar_guard = m;
    }
    require(
        r.offset == s.bytes.len(),
        "PIC transport envelope has trailing bytes",
    )?;
    Ok(e)
}

pub fn serialize_pic_transport_envelopes_3d(
    envelopes: &[PicTransportEnvelope3D],
) -> Result<Vec<PicSerializedEnvelope3D>> {
    envelopes
        .iter()
        .map(serialize_pic_transport_envelope_3d)
        .collect()
}

pub fn deserialize_pic_transport_envelopes_3d(
    serialized: &[PicSerializedEnvelope3D],
) -> Result<Vec<PicTransportEnvelope3D>> {
    serialized
        .iter()
        .map(deserialize_pic_transport_envelope_3d)
        .collect()
}

/// Serialize envelopes into per-destination-rank buckets with message and byte counts.
pub fn plan_serialized_pic_exchange_3d(
    envelopes: &[PicTransportEnvelope3D],
    rank_count: usize,
) -> Result<PicSerializedExchangePlan3D> {
    require(
        rank_count > 0,
        "PIC serialized exchange plan requires ranks",
    )?;
    let mut plan = PicSerializedExchangePlan3D {
        rank_count,
        envelopes_by_destination_rank: vec![Vec::new(); rank_count],
        message_count_by_destination_rank: vec![0; rank_count],
        byte_count_by_destination_rank: vec![0; rank_count],
    };
    for e in envelopes {
        require(
            e.source_rank < rank_count,
            "PIC serialized exchange source rank out of range",
        )?;
        require(
            e.destination_rank < rank_count,
            "PIC serialized exchange destination rank out of range",
        )?;
        let s = serialize_pic_transport_envelope_3d(e)?;
        plan.byte_count_by_destination_rank[e.destination_rank] += s.bytes.len();
        plan.message_count_by_destination_rank[e.destination_rank] += 1;
        plan.envelopes_by_destination_rank[e.destination_rank].push(s);
    }
    Ok(plan)
}

pub fn summarize_pic_serialized_transport_3d(
    serialized: &[PicSerializedEnvelope3D],
) -> Result<PicTransportDiagnostics3D> {
    Ok(summarize_pic_transport_3d(
        &deserialize_pic_transport_envelopes_3d(serialized)?,
    ))
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DistributedPicExchangeRound3D {
    pub particles_by_domain: Vec<Vec<PicParticle3D>>,
    pub scalar_guarded_blocks: Vec<ScalarGuardedBlock3D>,
    pub transport_diagnostics: PicTransportDiagnostics3D,
    pub serialized_plan: PicSerializedExchangePlan3D,
    pub serialized_message_count: usize,
    pub serialized_byte_count: usize,
    pub outside_particles: Vec<ParticleMigrationRecord3D>,
}

/// Pack particle migration and scalar guard messages, route them through the serialized
/// wire format by destination rank, decode and apply them.
pub fn run_serialized_distributed_pic_exchange_round_3d(
    particles_by_domain: &[Vec<PicParticle3D>],
    values: &[Vec<f64>],
    domains: &[PicDomain3D],
    c: &PicDomainGrid3DConfig,
    g: usize,
    exterior_value: f64,
    t: &PicRankTopology3D,
) -> Result<DistributedPicExchangeRound3D> {
    validate_domain_span(domains, c)?;
    require(
        particles_by_domain.len() == domains.len(),
        "distributed PIC exchange particle domain count mismatch",
    )?;
    require(
        values.len() == domains.len(),
        "distributed PIC exchange scalar domain count mismatch",
    )?;
    require(
        t.domain_to_rank.len() == domains.len(),
        "distributed PIC exchange topology/domain mismatch",
    )?;
    let packed = pack_particle_migration_messages_3d(particles_by_domain, domains, c)?;
    let guards = pack_scalar_guard_cell_messages_3d(values, domains, c, g)?;
    let envelopes = make_pic_transport_envelopes_3d(&packed.messages, &guards, t)?;
    let plan = plan_serialized_pic_exchange_3d(&envelopes, t.rank_count)?;
    let mut delivered = Vec::new();
    for bucket in &plan.envelopes_by_destination_rank {
        delivered.extend(deserialize_pic_transport_envelopes_3d(bucket)?);
    }
    let (pm, sm) = split_delivered(delivered);
    let mut r = DistributedPicExchangeRound3D {
        particles_by_domain: apply_particle_migration_messages_3d(packed.retained_by_domain, &pm)?,
        scalar_guarded_blocks: apply_scalar_guard_cell_messages_3d(
            values,
            domains,
            g,
            exterior_value,
            &sm,
        )?,
        transport_diagnostics: summarize_pic_transport_3d(&envelopes),
        serialized_message_count: plan.message_count_by_destination_rank.iter().sum(),
        serialized_byte_count: plan.byte_count_by_destination_rank.iter().sum(),
        serialized_plan: plan,
        outside_particles: packed.outside_particles,
    };
    r.serialized_plan.rank_count = t.rank_count;
    Ok(r)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistributedStaggeredPicStep3DConfig {
    pub guard_cells: usize,
    pub exterior_value: f64,
    pub dt_s: f64,
    pub exchange_em_field_guards: bool,
}

impl Default for DistributedStaggeredPicStep3DConfig {
    fn default() -> Self {
        Self {
            guard_cells: 1,
            exterior_value: 0.0,
            dt_s: 1.0e-12,
            exchange_em_field_guards: true,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DistributedStaggeredPicStep3D {
    pub particles_by_domain: Vec<Vec<PicParticle3D>>,
    pub guarded_fields: ElectromagneticGuardedFieldBlocks3D,
    pub exchange_round: DistributedPicExchangeRound3D,
    pub particle_count_before: usize,
    pub particle_count_after_push: usize,
    pub particle_count_after_migration: usize,
    pub max_particle_displacement_m: f64,
}

fn local_sample(values: &[f64], d: &PicDomain3D, p: Vec3) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let dx = (d.maximum_m.x - d.minimum_m.x) / d.cells_x as f64;
    let dy = (d.maximum_m.y - d.minimum_m.y) / d.cells_y as f64;
    let dz = (d.maximum_m.z - d.minimum_m.z) / d.cells_z as f64;
    let cell = |u: f64, n: usize| -> (usize, f64) {
        if u <= 0.0 {
            return (0, 0.0);
        }
        if u >= n as f64 - 1.0e-12 {
            return (n - 1, 0.0);
        }
        let f = u.floor();
        (f as usize, u - f)
    };
    let (ix0, fx) = cell((p.x - d.minimum_m.x) / dx, d.cells_x);
    let (iy0, fy) = cell((p.y - d.minimum_m.y) / dy, d.cells_y);
    let (iz0, fz) = cell((p.z - d.minimum_m.z) / dz, d.cells_z);
    let (ix1, iy1, iz1) = (
        (ix0 + 1).min(d.cells_x - 1),
        (iy0 + 1).min(d.cells_y - 1),
        (iz0 + 1).min(d.cells_z - 1),
    );
    let mut v = 0.0;
    for (iz, wz) in [(iz0, 1.0 - fz), (iz1, fz)] {
        for (iy, wy) in [(iy0, 1.0 - fy), (iy1, fy)] {
            for (ix, wx) in [(ix0, 1.0 - fx), (ix1, fx)] {
                v += values[fi3(ix, iy, iz, d.cells_x, d.cells_y)] * wx * wy * wz;
            }
        }
    }
    v
}

/// One distributed step: push every domain's particles in its local (cell-centred, clamped)
/// field with Vay, optionally exchange field guards, then run a serialized exchange round
/// (the Ex blocks ride along as the scalar guard payload, as upstream).
pub fn run_serialized_distributed_staggered_pic_step_3d(
    particles_by_domain: &[Vec<PicParticle3D>],
    f: &ElectromagneticFieldBlocks3D,
    domains: &[PicDomain3D],
    c: &PicDomainGrid3DConfig,
    s: &DistributedStaggeredPicStep3DConfig,
    t: &PicRankTopology3D,
) -> Result<DistributedStaggeredPicStep3D> {
    validate_domain_span(domains, c)?;
    require(
        s.dt_s > 0.0 && s.dt_s.is_finite(),
        "distributed PIC step requires positive dt",
    )?;
    require(
        particles_by_domain.len() == domains.len(),
        "distributed PIC step particle domain count mismatch",
    )?;
    require(
        t.domain_to_rank.len() == domains.len(),
        "distributed PIC step topology/domain mismatch",
    )?;
    let n = domains.len();
    let comps = [
        &f.electric_x_by_domain,
        &f.electric_y_by_domain,
        &f.electric_z_by_domain,
        &f.magnetic_x_by_domain,
        &f.magnetic_y_by_domain,
        &f.magnetic_z_by_domain,
    ];
    require(
        comps.iter().all(|c| c.len() == n),
        "distributed PIC step field domain count mismatch",
    )?;
    for d in domains {
        let cells = d.cells_x * d.cells_y * d.cells_z;
        require(
            comps.iter().all(|c| c[d.domain_id].len() == cells),
            "distributed PIC step field array size mismatch",
        )?;
    }
    let mut advanced = particles_by_domain.to_vec();
    let mut r = DistributedStaggeredPicStep3D {
        particle_count_before: particles_by_domain.iter().map(Vec::len).sum(),
        ..Default::default()
    };
    for d in domains {
        let id = d.domain_id;
        for particle in advanced[id].iter_mut() {
            let old = particle.position_m;
            let field = ElectromagneticField {
                electric_v_per_m: Vec3::new(
                    local_sample(&f.electric_x_by_domain[id], d, old),
                    local_sample(&f.electric_y_by_domain[id], d, old),
                    local_sample(&f.electric_z_by_domain[id], d, old),
                ),
                magnetic_t: Vec3::new(
                    local_sample(&f.magnetic_x_by_domain[id], d, old),
                    local_sample(&f.magnetic_y_by_domain[id], d, old),
                    local_sample(&f.magnetic_z_by_domain[id], d, old),
                ),
            };
            let mut charged = particle.proxy();
            vay_push(&mut charged, &field, s.dt_s)?;
            particle.position_m = charged.position_m;
            particle.velocity_m_per_s = charged.velocity_m_per_s;
            let delta = Vec3::new(
                particle.position_m.x - old.x,
                particle.position_m.y - old.y,
                particle.position_m.z - old.z,
            );
            r.max_particle_displacement_m = r.max_particle_displacement_m.max(norm2(delta).sqrt());
        }
    }
    r.particle_count_after_push = advanced.iter().map(Vec::len).sum();
    if s.exchange_em_field_guards {
        r.guarded_fields = exchange_electromagnetic_field_guard_cells_3d(
            f,
            domains,
            c,
            s.guard_cells,
            s.exterior_value,
        )?;
    }
    r.exchange_round = run_serialized_distributed_pic_exchange_round_3d(
        &advanced,
        &f.electric_x_by_domain,
        domains,
        c,
        s.guard_cells,
        s.exterior_value,
        t,
    )?;
    r.particles_by_domain = r.exchange_round.particles_by_domain.clone();
    r.particle_count_after_migration = r.particles_by_domain.iter().map(Vec::len).sum();
    Ok(r)
}
