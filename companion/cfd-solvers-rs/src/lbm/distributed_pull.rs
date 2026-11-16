//! Brick-decomposed BGK/Guo pull LBM on halo grids with full or lattice-selective in-process
//! ("virtual") exchange, macroscopic gather and binary per-brick checkpoints.
#![allow(clippy::needless_range_loop)]
use super::{guo_force, lattice_equilibrium, LatticeDescriptor, MacroscopicFields};
use crate::{
    distributed::{
        neighbor_offsets, neighbor_rank, pack_selective, selective_halo_value_count,
        unpack_selective, Brick, Extent3, HaloSoA3D, NeighborOffset,
    },
    require, Result, SolverError,
};
use std::marker::PhantomData;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistributedLbmConfig {
    pub global: Extent3,
    pub tau: f32,
    pub acceleration_x: f32,
    pub acceleration_y: f32,
    pub acceleration_z: f32,
}

impl Default for DistributedLbmConfig {
    fn default() -> Self {
        Self {
            global: Extent3::new(128, 128, 128),
            tau: 0.6,
            acceleration_x: 0.0,
            acceleration_y: 0.0,
            acceleration_z: 0.0,
        }
    }
}

/// One brick of a decomposed lattice; `Q` must equal `D::Q`.
pub struct DistributedPullBlock<D: LatticeDescriptor, const Q: usize> {
    config: DistributedLbmConfig,
    brick: Brick,
    current: HaloSoA3D<f32, Q>,
    next: HaloSoA3D<f32, Q>,
    time_step: u64,
    _d: PhantomData<fn() -> D>,
}

impl<D: LatticeDescriptor, const Q: usize> Clone for DistributedPullBlock<D, Q> {
    fn clone(&self) -> Self {
        Self {
            config: self.config,
            brick: self.brick,
            current: self.current.clone(),
            next: self.next.clone(),
            time_step: self.time_step,
            _d: PhantomData,
        }
    }
}

impl<D: LatticeDescriptor, const Q: usize> std::fmt::Debug for DistributedPullBlock<D, Q> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DistributedPullBlock")
            .field("brick", &self.brick)
            .field("time_step", &self.time_step)
            .finish()
    }
}

const MAGIC: [u8; 8] = *b"CFDCHK3\0";

impl<D: LatticeDescriptor, const Q: usize> DistributedPullBlock<D, Q> {
    pub fn new(config: DistributedLbmConfig, brick: Brick) -> Result<Self> {
        require(
            Q == D::Q,
            "descriptor population count must match halo-grid components",
        )?;
        require(config.tau > 0.5, "LBM tau must be > 0.5")?;
        require(
            brick.global == config.global,
            "brick/global LBM extent mismatch",
        )?;
        if D::DIMENSIONS == 2 {
            require(
                config.global.z == 1 && brick.extent.z == 1,
                "2-D distributed LBM requires z=1",
            )?;
        }
        let mut b = Self {
            config,
            brick,
            current: HaloSoA3D::new(brick.extent, 1)?,
            next: HaloSoA3D::new(brick.extent, 1)?,
            time_step: 0,
            _d: PhantomData,
        };
        b.initialize_uniform(1.0, 0.0, 0.0, 0.0)?;
        Ok(b)
    }
    pub fn initialize(&mut self, init: impl Fn(usize, usize, usize) -> [f32; 4]) -> Result<()> {
        self.time_step = 0;
        let e = self.brick.extent;
        for z in 1..=e.z {
            for y in 1..=e.y {
                for x in 1..=e.x {
                    let s = init(
                        self.brick.begin.x + x - 1,
                        self.brick.begin.y + y - 1,
                        self.brick.begin.z + z - 1,
                    );
                    for d in 0..Q {
                        let v = lattice_equilibrium::<D>(d, s[0], s[1], s[2], s[3]);
                        self.current.set(d, x, y, z, v)?;
                        self.next.set(d, x, y, z, v)?;
                    }
                }
            }
        }
        Ok(())
    }
    /// Upstream defaults `(1, 0, 0, 0)`.
    pub fn initialize_uniform(&mut self, rho: f32, ux: f32, uy: f32, uz: f32) -> Result<()> {
        self.initialize(move |_, _, _| [rho, ux, uy, uz])
    }
    /// Upstream default amplitude 0.03.
    pub fn initialize_taylor_green(&mut self, amplitude: f32) -> Result<()> {
        let g = self.config.global;
        let two_pi = 2.0 * std::f32::consts::PI;
        self.initialize(move |gx, gy, gz| {
            let xf = (gx as f32 + 0.5) / g.x as f32;
            let yf = (gy as f32 + 0.5) / g.y as f32;
            let (sx, cx, sy, cy) = (
                (two_pi * xf).sin(),
                (two_pi * xf).cos(),
                (two_pi * yf).sin(),
                (two_pi * yf).cos(),
            );
            let zf = if D::DIMENSIONS == 3 {
                (two_pi * ((gz as f32 + 0.5) / g.z as f32)).cos()
            } else {
                1.0
            };
            [
                1.0,
                amplitude * sx * cy * zf,
                -amplitude * cx * sy * zf,
                0.0,
            ]
        })
    }
    fn local(&self, n: usize) -> (usize, usize, usize) {
        let e = self.brick.extent;
        let (x, yz) = (n % e.x, n / e.x);
        (x + 1, yz % e.y + 1, yz / e.y + 1)
    }
    fn is_strict_interior(&self, x: usize, y: usize, z: usize) -> bool {
        let e = self.brick.extent;
        let xy = x > 1 && x < e.x && y > 1 && y < e.y;
        xy && (D::DIMENSIONS == 2 || (z > 1 && z < e.z))
    }
    fn collide_pull(&mut self, x: usize, y: usize, z: usize) -> Result<()> {
        let mut fin = [0.0f32; Q];
        let (mut rho, mut mx, mut my, mut mz) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
        for d in 0..Q {
            let xs = (x as i64 - D::cx(d) as i64) as usize;
            let ys = (y as i64 - D::cy(d) as i64) as usize;
            let zs = (z as i64 - D::cz(d) as i64) as usize;
            let v = self.current.get(d, xs, ys, zs)?;
            fin[d] = v;
            rho += v;
            mx += v * D::cx(d) as f32;
            my += v * D::cy(d) as f32;
            mz += v * D::cz(d) as f32;
        }
        let c = &self.config;
        let ux = if rho > 0.0 { mx / rho } else { 0.0 } + 0.5 * c.acceleration_x;
        let uy = if rho > 0.0 { my / rho } else { 0.0 } + 0.5 * c.acceleration_y;
        let uz = if rho > 0.0 { mz / rho } else { 0.0 } + 0.5 * c.acceleration_z;
        let omega = 1.0 / c.tau;
        for d in 0..Q {
            let feq = lattice_equilibrium::<D>(d, rho, ux, uy, uz);
            let v = fin[d] - omega * (fin[d] - feq)
                + guo_force::<D>(
                    d,
                    rho,
                    ux,
                    uy,
                    uz,
                    c.acceleration_x,
                    c.acceleration_y,
                    c.acceleration_z,
                    omega,
                );
            self.next.set(d, x, y, z, v)?;
        }
        Ok(())
    }
    fn compute_where(&mut self, keep: impl Fn(&Self, usize, usize, usize) -> bool) -> Result<()> {
        for n in 0..self.brick.extent.cells() {
            let (x, y, z) = self.local(n);
            if keep(self, x, y, z) {
                self.collide_pull(x, y, z)?;
            }
        }
        Ok(())
    }
    /// Cells whose pull stencil does not touch the halo (skipped for thin bricks).
    pub fn compute_strict_interior(&mut self) -> Result<()> {
        let e = self.brick.extent;
        if e.x <= 2 || e.y <= 2 || (D::DIMENSIONS == 3 && e.z <= 2) {
            return Ok(());
        }
        self.compute_where(|s, x, y, z| s.is_strict_interior(x, y, z))
    }
    pub fn compute_boundary(&mut self) -> Result<()> {
        self.compute_where(|s, x, y, z| !s.is_strict_interior(x, y, z))
    }
    pub fn compute_all(&mut self) -> Result<()> {
        self.compute_where(|_, _, _, _| true)
    }
    pub fn finish_step(&mut self) {
        std::mem::swap(&mut self.current, &mut self.next);
        self.time_step += 1;
    }
    pub fn compute_macroscopic(&self) -> Result<MacroscopicFields> {
        let n = self.brick.extent.cells();
        let mut f = MacroscopicFields {
            rho: vec![0.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            uz: vec![0.0; n],
        };
        let c = &self.config;
        for i in 0..n {
            let (x, y, z) = self.local(i);
            let (mut rho, mut mx, mut my, mut mz) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
            for d in 0..Q {
                let v = self.current.get(d, x, y, z)?;
                rho += v;
                mx += v * D::cx(d) as f32;
                my += v * D::cy(d) as f32;
                mz += v * D::cz(d) as f32;
            }
            f.rho[i] = rho;
            if rho > 0.0 {
                f.ux[i] = mx / rho + 0.5 * c.acceleration_x;
                f.uy[i] = my / rho + 0.5 * c.acceleration_y;
                f.uz[i] = mz / rho + 0.5 * c.acceleration_z;
            }
        }
        Ok(f)
    }
    pub fn local_mass(&self) -> Result<f64> {
        let mut total = 0.0;
        for i in 0..self.brick.extent.cells() {
            let (x, y, z) = self.local(i);
            let mut rho = 0.0f64;
            for d in 0..Q {
                rho += self.current.get(d, x, y, z)? as f64;
            }
            total += rho;
        }
        Ok(total)
    }
    pub fn config(&self) -> &DistributedLbmConfig {
        &self.config
    }
    pub fn brick(&self) -> &Brick {
        &self.brick
    }
    pub fn current_grid(&self) -> &HaloSoA3D<f32, Q> {
        &self.current
    }
    pub fn current_grid_mut(&mut self) -> &mut HaloSoA3D<f32, Q> {
        &mut self.current
    }
    pub fn time_step(&self) -> u64 {
        self.time_step
    }
    /// Upstream `CFDCHK3` layout in little-endian byte order.
    pub fn save_checkpoint(&self, path: impl AsRef<Path>) -> Result<()> {
        let mut out = Vec::new();
        out.extend_from_slice(&MAGIC);
        for v in [1u32, 0x0102_0304, Q as u32] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        for e in [self.config.global, self.brick.begin, self.brick.extent] {
            for v in [e.x, e.y, e.z] {
                out.extend_from_slice(&(v as u64).to_le_bytes());
            }
        }
        out.extend_from_slice(&self.time_step.to_le_bytes());
        for v in [
            self.config.tau,
            self.config.acceleration_x,
            self.config.acceleration_y,
            self.config.acceleration_z,
        ] {
            out.extend_from_slice(&v.to_le_bytes());
        }
        let e = self.brick.extent;
        for d in 0..Q {
            for z in 1..=e.z {
                for y in 1..=e.y {
                    for x in 1..=e.x {
                        out.extend_from_slice(&self.current.get(d, x, y, z)?.to_le_bytes());
                    }
                }
            }
        }
        std::fs::write(path.as_ref(), out).map_err(|e| {
            SolverError::Io(format!(
                "cannot write checkpoint {}: {e}",
                path.as_ref().display()
            ))
        })
    }
    /// Restores a checkpoint written for the same brick and configuration; the block is unchanged
    /// on error.
    pub fn load_checkpoint(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let bytes = std::fs::read(path.as_ref()).map_err(|e| {
            SolverError::Io(format!(
                "cannot open checkpoint {}: {e}",
                path.as_ref().display()
            ))
        })?;
        let mut pos = 0usize;
        let mut take = |n: usize| -> Result<&[u8]> {
            let s = bytes
                .get(pos..pos + n)
                .ok_or_else(|| SolverError::Io("checkpoint read failed".into()))?;
            pos += n;
            Ok(s)
        };
        let fail = |m: &str| SolverError::Io(m.into());
        if take(8)? != MAGIC {
            return Err(fail("invalid checkpoint magic"));
        }
        let u32v = |s: &[u8]| u32::from_le_bytes([s[0], s[1], s[2], s[3]]);
        let (version, endian, q) = (u32v(take(4)?), u32v(take(4)?), u32v(take(4)?));
        if version != 1 || endian != 0x0102_0304 || q != Q as u32 {
            return Err(fail("unsupported checkpoint format/descriptor"));
        }
        let mut ext = [Extent3::default(); 3];
        for e in ext.iter_mut() {
            let mut v = [0usize; 3];
            for c in v.iter_mut() {
                let s = take(8)?;
                *c = u64::from_le_bytes([s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]]) as usize;
            }
            *e = Extent3::new(v[0], v[1], v[2]);
        }
        if ext[0] != self.config.global || ext[1] != self.brick.begin || ext[2] != self.brick.extent
        {
            return Err(fail("checkpoint decomposition does not match block"));
        }
        let s = take(8)?;
        let step = u64::from_le_bytes([s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]]);
        let mut f = [0f32; 4];
        for v in f.iter_mut() {
            *v = f32::from_le_bytes(
                take(4)?
                    .try_into()
                    .map_err(|_| fail("checkpoint read failed"))?,
            );
        }
        let c = &self.config;
        if f != [c.tau, c.acceleration_x, c.acceleration_y, c.acceleration_z] {
            return Err(fail("checkpoint LBM configuration mismatch"));
        }
        let e = self.brick.extent;
        let (mut cur, mut nxt) = (self.current.clone(), self.next.clone());
        for d in 0..Q {
            for z in 1..=e.z {
                for y in 1..=e.y {
                    for x in 1..=e.x {
                        let v = f32::from_le_bytes(
                            take(4)?
                                .try_into()
                                .map_err(|_| fail("checkpoint read failed"))?,
                        );
                        cur.set(d, x, y, z, v)?;
                        nxt.set(d, x, y, z, v)?;
                    }
                }
            }
        }
        self.current = cur;
        self.next = nxt;
        self.time_step = step;
        Ok(())
    }
}

/// Full 26-neighbor halo exchange between in-process blocks ordered by rank.
pub fn virtual_exchange<D: LatticeDescriptor, const Q: usize>(
    blocks: &mut [DistributedPullBlock<D, Q>],
    periodic: [bool; 3],
) -> Result<()> {
    let mut msgs: Vec<(usize, NeighborOffset, Vec<f32>)> = vec![];
    for (rank, b) in blocks.iter().enumerate() {
        require(
            b.brick.rank == rank,
            "virtual blocks must be ordered by rank",
        )?;
        for o in neighbor_offsets() {
            if let Some(nb) = neighbor_rank(&b.brick, o, periodic)? {
                msgs.push((nb, o.opposite(), b.current.pack(o)));
            }
        }
    }
    for (dest, from, payload) in msgs {
        blocks
            .get_mut(dest)
            .ok_or_else(|| SolverError::InvalidInput("neighbor rank outside the block set".into()))?
            .current
            .unpack(from, &payload)?;
    }
    Ok(())
}

/// Exchange of only the populations that cross each face/edge/corner.
pub fn virtual_selective_exchange<D: LatticeDescriptor, const Q: usize>(
    blocks: &mut [DistributedPullBlock<D, Q>],
    periodic: [bool; 3],
) -> Result<()> {
    let mut msgs: Vec<(usize, NeighborOffset, Vec<f32>)> = vec![];
    for (rank, b) in blocks.iter().enumerate() {
        require(
            b.brick.rank == rank,
            "virtual blocks must be ordered by rank",
        )?;
        for o in neighbor_offsets() {
            if selective_halo_value_count::<D>(b.current.interior_extent(), o) == 0 {
                continue;
            }
            if let Some(nb) = neighbor_rank(&b.brick, o, periodic)? {
                msgs.push((
                    nb,
                    o.opposite(),
                    pack_selective::<D, f32, Q>(&b.current, o)?,
                ));
            }
        }
    }
    for (dest, from, payload) in msgs {
        let b = blocks.get_mut(dest).ok_or_else(|| {
            SolverError::InvalidInput("neighbor rank outside the block set".into())
        })?;
        unpack_selective::<D, f32, Q>(&mut b.current, from, &payload)?;
    }
    Ok(())
}

fn virtual_step_with<D: LatticeDescriptor, const Q: usize>(
    blocks: &mut [DistributedPullBlock<D, Q>],
    periodic: [bool; 3],
    selective: bool,
) -> Result<()> {
    let mut trial = blocks.to_vec();
    for b in trial.iter_mut() {
        b.compute_strict_interior()?;
    }
    if selective {
        virtual_selective_exchange(&mut trial, periodic)?;
    } else {
        virtual_exchange(&mut trial, periodic)?;
    }
    for b in trial.iter_mut() {
        b.compute_boundary()?;
        b.finish_step();
    }
    blocks.clone_from_slice(&trial);
    Ok(())
}

/// Interior collide-pull, full exchange, boundary collide-pull and swap (transactional).
pub fn virtual_step<D: LatticeDescriptor, const Q: usize>(
    blocks: &mut [DistributedPullBlock<D, Q>],
    periodic: [bool; 3],
) -> Result<()> {
    virtual_step_with(blocks, periodic, false)
}

pub fn virtual_selective_step<D: LatticeDescriptor, const Q: usize>(
    blocks: &mut [DistributedPullBlock<D, Q>],
    periodic: [bool; 3],
) -> Result<()> {
    virtual_step_with(blocks, periodic, true)
}

pub fn gather_virtual_macroscopic<D: LatticeDescriptor, const Q: usize>(
    blocks: &[DistributedPullBlock<D, Q>],
) -> Result<MacroscopicFields> {
    let Some(first) = blocks.first() else {
        return Ok(MacroscopicFields::default());
    };
    let g = first.config.global;
    let n = g.cells();
    let mut r = MacroscopicFields {
        rho: vec![0.0; n],
        ux: vec![0.0; n],
        uy: vec![0.0; n],
        uz: vec![0.0; n],
    };
    for b in blocks {
        let m = b.compute_macroscopic()?;
        let br = &b.brick;
        let mut i = 0;
        for z in 0..br.extent.z {
            for y in 0..br.extent.y {
                for x in 0..br.extent.x {
                    let gi = ((br.begin.z + z) * g.y + br.begin.y + y) * g.x + br.begin.x + x;
                    require(gi < n, "brick outside the global lattice")?;
                    r.rho[gi] = m.rho[i];
                    r.ux[gi] = m.ux[i];
                    r.uy[gi] = m.uy[i];
                    r.uz[gi] = m.uz[i];
                    i += 1;
                }
            }
        }
    }
    Ok(r)
}
