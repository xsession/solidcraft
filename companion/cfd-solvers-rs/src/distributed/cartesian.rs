//! Cartesian bricks, neighbor topology, width-1 SoA halo grids with full and lattice-selective
//! packing, accelerator assignment arithmetic and an in-process halo exchange that emulates the
//! MPI neighbor exchange.
use crate::{lbm::LatticeDescriptor, require, Result, SolverError};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Extent3 {
    pub x: usize,
    pub y: usize,
    pub z: usize,
}

impl Default for Extent3 {
    fn default() -> Self {
        Self { x: 1, y: 1, z: 1 }
    }
}

impl Extent3 {
    pub const fn new(x: usize, y: usize, z: usize) -> Self {
        Self { x, y, z }
    }
    pub fn cells(&self) -> usize {
        self.x * self.y * self.z
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Coord3 {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

pub type ProcessGrid = Extent3;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Brick {
    pub global: Extent3,
    pub begin: Extent3,
    pub extent: Extent3,
    pub processes: ProcessGrid,
    pub coordinate: Coord3,
    pub rank: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NeighborOffset {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

impl NeighborOffset {
    pub const fn new(x: i32, y: i32, z: i32) -> Self {
        Self { x, y, z }
    }
    pub fn opposite(&self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
    pub fn is_center(&self) -> bool {
        self.x == 0 && self.y == 0 && self.z == 0
    }
}

/// The 26 non-center offsets, x fastest.
pub fn neighbor_offsets() -> [NeighborOffset; 26] {
    let mut out = [NeighborOffset::default(); 26];
    let mut n = 0;
    for z in -1..=1 {
        for y in -1..=1 {
            for x in -1..=1 {
                if x == 0 && y == 0 && z == 0 {
                    continue;
                }
                out[n] = NeighborOffset::new(x, y, z);
                n += 1;
            }
        }
    }
    out
}

pub fn neighbor_tag(o: NeighborOffset) -> i32 {
    (o.z + 1) * 9 + (o.y + 1) * 3 + (o.x + 1)
}

/// Minimizes halo surface with a small aspect-ratio penalty (upstream default 3 dimensions).
pub fn choose_process_grid(global: Extent3, ranks: usize, dimensions: u32) -> Result<ProcessGrid> {
    require(ranks > 0, "process count must be positive")?;
    require(
        global.x > 0 && global.y > 0 && global.z > 0,
        "global extents must be positive",
    )?;
    require(
        dimensions == 2 || dimensions == 3,
        "dimensions must be 2 or 3",
    )?;
    let (mut best, mut best_score, mut found) = (Extent3::default(), f64::INFINITY, false);
    for px in 1..=ranks.min(global.x) {
        if ranks % px != 0 {
            continue;
        }
        let rem = ranks / px;
        for py in 1..=rem.min(global.y) {
            if rem % py != 0 {
                continue;
            }
            let pz = rem / py;
            if (dimensions == 2 && pz != 1) || pz > global.z {
                continue;
            }
            let (lx, ly, lz) = (
                global.x as f64 / px as f64,
                global.y as f64 / py as f64,
                global.z as f64 / pz as f64,
            );
            let surface = if dimensions == 2 {
                2.0 * (lx + ly)
            } else {
                2.0 * (lx * ly + lx * lz + ly * lz)
            };
            let max_l = lx.max(ly).max(lz);
            let min_l = lx.min(ly).min(lz);
            let score = surface * (1.0 + 0.02 * (max_l / min_l.max(1e-12)));
            if score < best_score {
                best_score = score;
                best = Extent3::new(px, py, pz);
                found = true;
            }
        }
    }
    if !found {
        return Err(SolverError::InvalidInput(
            "cannot map ranks onto non-empty Cartesian bricks".into(),
        ));
    }
    Ok(best)
}

pub fn coords_from_rank(rank: usize, grid: ProcessGrid) -> Result<Coord3> {
    require(rank < grid.cells(), "rank outside process grid")?;
    Ok(Coord3 {
        x: (rank % grid.x) as i32,
        y: ((rank / grid.x) % grid.y) as i32,
        z: (rank / (grid.x * grid.y)) as i32,
    })
}

fn in_grid(c: Coord3, g: ProcessGrid) -> bool {
    c.x >= 0
        && c.y >= 0
        && c.z >= 0
        && (c.x as usize) < g.x
        && (c.y as usize) < g.y
        && (c.z as usize) < g.z
}

pub fn rank_from_coords(c: Coord3, grid: ProcessGrid) -> Result<usize> {
    require(in_grid(c, grid), "coordinates outside process grid")?;
    Ok((c.z as usize * grid.y + c.y as usize) * grid.x + c.x as usize)
}

/// `(begin, extent)` of a balanced 1-D split.
pub fn partition_axis(
    global: usize,
    coordinate: usize,
    partitions: usize,
) -> Result<(usize, usize)> {
    require(
        partitions > 0 && coordinate < partitions,
        "invalid axis partition",
    )?;
    let (base, extra) = (global / partitions, global % partitions);
    Ok((
        coordinate * base + coordinate.min(extra),
        base + usize::from(coordinate < extra),
    ))
}

pub fn make_brick_at(
    global: Extent3,
    grid: ProcessGrid,
    coord: Coord3,
    rank: usize,
) -> Result<Brick> {
    require(
        grid.cells() > 0 && rank < grid.cells(),
        "invalid process grid/rank",
    )?;
    require(in_grid(coord, grid), "invalid process-grid coordinates")?;
    let (bx, nx) = partition_axis(global.x, coord.x as usize, grid.x)?;
    let (by, ny) = partition_axis(global.y, coord.y as usize, grid.y)?;
    let (bz, nz) = partition_axis(global.z, coord.z as usize, grid.z)?;
    Ok(Brick {
        global,
        begin: Extent3::new(bx, by, bz),
        extent: Extent3::new(nx, ny, nz),
        processes: grid,
        coordinate: coord,
        rank,
    })
}

pub fn make_brick(global: Extent3, grid: ProcessGrid, rank: usize) -> Result<Brick> {
    make_brick_at(global, grid, coords_from_rank(rank, grid)?, rank)
}

/// Neighbor rank, or `None` across a non-periodic boundary (upstream `-1`).
pub fn neighbor_rank(
    brick: &Brick,
    offset: NeighborOffset,
    periodic: [bool; 3],
) -> Result<Option<usize>> {
    if offset.is_center() {
        return Ok(Some(brick.rank));
    }
    let dims = [
        brick.processes.x as i32,
        brick.processes.y as i32,
        brick.processes.z as i32,
    ];
    let mut c = [
        brick.coordinate.x + offset.x,
        brick.coordinate.y + offset.y,
        brick.coordinate.z + offset.z,
    ];
    for a in 0..3 {
        if c[a] < 0 || c[a] >= dims[a] {
            if !periodic[a] {
                return Ok(None);
            }
            c[a] = c[a].rem_euclid(dims[a]);
        }
    }
    rank_from_coords(
        Coord3 {
            x: c[0],
            y: c[1],
            z: c[2],
        },
        brick.processes,
    )
    .map(Some)
}

/// Component-major SoA field with a one-cell halo on every side.
#[derive(Clone, Debug, PartialEq)]
pub struct HaloSoA3D<T: Copy + Default, const C: usize> {
    interior: Extent3,
    total: Extent3,
    data: Vec<T>,
}

impl<T: Copy + Default, const C: usize> HaloSoA3D<T, C> {
    /// Only halo width 1 is supported (as upstream).
    pub fn new(interior: Extent3, halo: usize) -> Result<Self> {
        require(
            halo == 1,
            "Phase-3 halo grid currently supports halo width 1",
        )?;
        require(
            interior.x > 0 && interior.y > 0 && interior.z > 0,
            "halo grid interior must be non-empty",
        )?;
        let total = Extent3::new(interior.x + 2, interior.y + 2, interior.z + 2);
        let n = total
            .cells()
            .checked_mul(C)
            .filter(|n| *n <= 1 << 30)
            .ok_or_else(|| SolverError::InvalidInput("halo grid too large".into()))?;
        Ok(Self {
            interior,
            total,
            data: vec![T::default(); n],
        })
    }
    pub fn interior_extent(&self) -> Extent3 {
        self.interior
    }
    pub fn total_extent(&self) -> Extent3 {
        self.total
    }
    pub fn halo(&self) -> usize {
        1
    }
    pub fn interior_cells(&self) -> usize {
        self.interior.cells()
    }
    pub fn index(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.total.y + y) * self.total.x + x
    }
    fn slot(&self, c: usize, x: usize, y: usize, z: usize) -> Result<usize> {
        require(
            c < C && x < self.total.x && y < self.total.y && z < self.total.z,
            "halo grid index out of range",
        )?;
        Ok(c * self.total.cells() + self.index(x, y, z))
    }
    pub fn get(&self, c: usize, x: usize, y: usize, z: usize) -> Result<T> {
        Ok(self.data[self.slot(c, x, y, z)?])
    }
    pub fn set(&mut self, c: usize, x: usize, y: usize, z: usize, value: T) -> Result<()> {
        let i = self.slot(c, x, y, z)?;
        self.data[i] = value;
        Ok(())
    }
    pub fn component(&self, c: usize) -> Result<&[T]> {
        require(c < C, "StaticSoA component")?;
        let n = self.total.cells();
        Ok(&self.data[c * n..(c + 1) * n])
    }
    pub fn halo_value_count(&self, o: NeighborOffset) -> usize {
        C * boundary_cells(self.interior, o)
    }
    fn send_range(n: usize, o: i32) -> (usize, usize) {
        match o.signum() {
            -1 => (1, 2),
            1 => (n, n + 1),
            _ => (1, n + 1),
        }
    }
    fn receive_range(n: usize, o: i32) -> (usize, usize) {
        match o.signum() {
            -1 => (0, 1),
            1 => (n + 1, n + 2),
            _ => (1, n + 1),
        }
    }
    pub fn pack(&self, offset: NeighborOffset) -> Vec<T> {
        let mut out = Vec::with_capacity(self.halo_value_count(offset));
        let (xb, xe) = Self::send_range(self.interior.x, offset.x);
        let (yb, ye) = Self::send_range(self.interior.y, offset.y);
        let (zb, ze) = Self::send_range(self.interior.z, offset.z);
        let n = self.total.cells();
        for c in 0..C {
            for z in zb..ze {
                for y in yb..ye {
                    for x in xb..xe {
                        out.push(self.data[c * n + self.index(x, y, z)]);
                    }
                }
            }
        }
        out
    }
    pub fn unpack(&mut self, incoming_from: NeighborOffset, values: &[T]) -> Result<()> {
        require(
            values.len() == self.halo_value_count(incoming_from),
            "halo payload size mismatch",
        )?;
        let (xb, xe) = Self::receive_range(self.interior.x, incoming_from.x);
        let (yb, ye) = Self::receive_range(self.interior.y, incoming_from.y);
        let (zb, ze) = Self::receive_range(self.interior.z, incoming_from.z);
        let n = self.total.cells();
        let mut p = 0;
        for c in 0..C {
            for z in zb..ze {
                for y in yb..ye {
                    for x in xb..xe {
                        let i = c * n + self.index(x, y, z);
                        self.data[i] = values[p];
                        p += 1;
                    }
                }
            }
        }
        Ok(())
    }
}

pub fn boundary_cells(interior: Extent3, o: NeighborOffset) -> usize {
    let nx = if o.x == 0 { interior.x } else { 1 };
    let ny = if o.y == 0 { interior.y } else { 1 };
    let nz = if o.z == 0 { interior.z } else { 1 };
    nx * ny * nz
}

pub fn population_crosses<D: LatticeDescriptor>(direction: usize, o: NeighborOffset) -> bool {
    !o.is_center()
        && (o.x == 0 || D::cx(direction) == o.x)
        && (o.y == 0 || D::cy(direction) == o.y)
        && (o.z == 0 || D::cz(direction) == o.z)
}

pub fn crossing_populations<D: LatticeDescriptor>(o: NeighborOffset) -> Vec<usize> {
    (0..D::Q)
        .filter(|d| population_crosses::<D>(*d, o))
        .collect()
}

pub fn selective_halo_value_count<D: LatticeDescriptor>(
    interior: Extent3,
    o: NeighborOffset,
) -> usize {
    crossing_populations::<D>(o).len() * boundary_cells(interior, o)
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SelectiveHaloEntry {
    pub direction: usize,
    pub x: usize,
    pub y: usize,
    pub z: usize,
}

fn selective_entry(
    directions: &[usize],
    interior: Extent3,
    o: NeighborOffset,
    linear: usize,
    receive: bool,
) -> SelectiveHaloEntry {
    let cells = boundary_cells(interior, o);
    if cells == 0 || linear >= directions.len() * cells {
        return SelectiveHaloEntry::default();
    }
    let di = linear / cells;
    let local = linear - di * cells;
    let rx = if o.x == 0 { interior.x } else { 1 };
    let ry = if o.y == 0 { interior.y } else { 1 };
    let (lx, yz) = (local % rx, local / rx);
    let (ly, lz) = (yz % ry, yz / ry);
    let coord = |off: i32, n: usize, l: usize| match (off.signum(), receive) {
        (-1, false) => 1,
        (1, false) => n,
        (-1, true) => 0,
        (1, true) => n + 1,
        _ => l + 1,
    };
    SelectiveHaloEntry {
        direction: directions[di],
        x: coord(o.x, interior.x, lx),
        y: coord(o.y, interior.y, ly),
        z: coord(o.z, interior.z, lz),
    }
}

pub fn selective_send_entry<D: LatticeDescriptor>(
    interior: Extent3,
    offset: NeighborOffset,
    linear: usize,
) -> SelectiveHaloEntry {
    selective_entry(
        &crossing_populations::<D>(offset),
        interior,
        offset,
        linear,
        false,
    )
}

pub fn selective_receive_entry<D: LatticeDescriptor>(
    interior: Extent3,
    incoming_from: NeighborOffset,
    linear: usize,
) -> SelectiveHaloEntry {
    selective_entry(
        &crossing_populations::<D>(incoming_from.opposite()),
        interior,
        incoming_from,
        linear,
        true,
    )
}

/// Packs only the populations that stream across the face/edge/corner toward `offset`.
pub fn pack_selective<D: LatticeDescriptor, T: Copy + Default, const C: usize>(
    grid: &HaloSoA3D<T, C>,
    offset: NeighborOffset,
) -> Result<Vec<T>> {
    require(
        C == D::Q,
        "descriptor population count must match halo-grid components",
    )?;
    let e = grid.interior_extent();
    let dirs = crossing_populations::<D>(offset);
    let n = dirs.len() * boundary_cells(e, offset);
    (0..n)
        .map(|p| {
            let en = selective_entry(&dirs, e, offset, p, false);
            grid.get(en.direction, en.x, en.y, en.z)
        })
        .collect()
}

pub fn unpack_selective<D: LatticeDescriptor, T: Copy + Default, const C: usize>(
    grid: &mut HaloSoA3D<T, C>,
    incoming_from: NeighborOffset,
    values: &[T],
) -> Result<()> {
    require(
        C == D::Q,
        "descriptor population count must match halo-grid components",
    )?;
    let e = grid.interior_extent();
    require(
        values.len() == selective_halo_value_count::<D>(e, incoming_from),
        "selective halo payload size mismatch",
    )?;
    let dirs = crossing_populations::<D>(incoming_from.opposite());
    for (p, v) in values.iter().enumerate() {
        let en = selective_entry(&dirs, e, incoming_from, p, true);
        grid.set(en.direction, en.x, en.y, en.z, *v)?;
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HaloTrafficStats {
    pub messages: usize,
    pub values: usize,
    pub bytes: usize,
}

/// Messages/values/bytes per step for selective exchange (`value_bytes` = upstream sizeof(T)).
pub fn selective_halo_traffic<D: LatticeDescriptor>(
    interior: Extent3,
    value_bytes: usize,
) -> HaloTrafficStats {
    let mut s = HaloTrafficStats::default();
    for o in neighbor_offsets() {
        let v = selective_halo_value_count::<D>(interior, o);
        if v == 0 {
            continue;
        }
        s.messages += 1;
        s.values += v;
    }
    s.bytes = s.values * value_bytes;
    s
}

pub fn full_halo_traffic<D: LatticeDescriptor>(
    interior: Extent3,
    value_bytes: usize,
) -> HaloTrafficStats {
    let mut s = HaloTrafficStats::default();
    for o in neighbor_offsets() {
        s.messages += 1;
        s.values += D::Q * boundary_cells(interior, o);
    }
    s.bytes = s.values * value_bytes;
    s
}

/// In-process emulation of the MPI 26-neighbor halo exchange: each rank's interior boundary layer
/// toward `offset` is unpacked into the neighbor's halo (`grids[rank]` must match `bricks[rank]`).
pub fn exchange_halos_in_process<T: Copy + Default, const C: usize>(
    grids: &mut [HaloSoA3D<T, C>],
    bricks: &[Brick],
    periodic: [bool; 3],
) -> Result<()> {
    require(
        grids.len() == bricks.len(),
        "one halo grid per brick required",
    )?;
    for (g, b) in grids.iter().zip(bricks) {
        require(
            g.interior_extent() == b.extent,
            "halo grid extent must match its brick",
        )?;
    }
    let mut messages = vec![];
    for (rank, b) in bricks.iter().enumerate() {
        require(b.rank == rank, "bricks must be ordered by rank")?;
        for o in neighbor_offsets() {
            if let Some(dest) = neighbor_rank(b, o, periodic)? {
                messages.push((dest, o.opposite(), grids[rank].pack(o)));
            }
        }
    }
    for (dest, from, payload) in messages {
        let g = grids.get_mut(dest).ok_or_else(|| {
            SolverError::InvalidInput("neighbor rank outside the brick set".into())
        })?;
        g.unpack(from, &payload)?;
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DeviceAssignment {
    pub local_rank: usize,
    pub local_size: usize,
    pub visible_devices: usize,
    pub device_ordinal: usize,
    pub oversubscribed: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RankDeviceGroup {
    pub local_rank: usize,
    pub local_size: usize,
    pub visible_devices: usize,
    pub device_ordinals: Vec<usize>,
    pub oversubscribed: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DeviceWorkSlice {
    pub device_ordinal: usize,
    pub begin: usize,
    pub end: usize,
}

pub fn assign_device(
    local_rank: usize,
    local_size: usize,
    visible_devices: usize,
) -> Result<DeviceAssignment> {
    require(
        local_size > 0 && local_rank < local_size,
        "invalid local rank topology",
    )?;
    if visible_devices == 0 {
        return Ok(DeviceAssignment {
            local_rank,
            local_size,
            ..Default::default()
        });
    }
    Ok(DeviceAssignment {
        local_rank,
        local_size,
        visible_devices,
        device_ordinal: local_rank % visible_devices,
        oversubscribed: local_size > visible_devices,
    })
}

/// `max_devices_per_rank = 0` means no cap.
pub fn assign_device_group(
    local_rank: usize,
    local_size: usize,
    visible_devices: usize,
    max_devices_per_rank: usize,
) -> Result<RankDeviceGroup> {
    require(
        local_size > 0 && local_rank < local_size,
        "invalid local rank topology",
    )?;
    let mut g = RankDeviceGroup {
        local_rank,
        local_size,
        visible_devices,
        ..Default::default()
    };
    if visible_devices == 0 {
        return Ok(g);
    }
    if visible_devices < local_size {
        g.device_ordinals.push(local_rank % visible_devices);
        g.oversubscribed = true;
        return Ok(g);
    }
    let (base, rem) = (visible_devices / local_size, visible_devices % local_size);
    let natural = base + usize::from(local_rank < rem);
    let count = if max_devices_per_rank == 0 {
        natural
    } else {
        natural.min(max_devices_per_rank)
    };
    let first = local_rank * base + local_rank.min(rem);
    g.device_ordinals = (first..first + count).collect();
    Ok(g)
}

pub fn split_rank_work(work_items: usize, group: &RankDeviceGroup) -> Vec<DeviceWorkSlice> {
    let n = group.device_ordinals.len();
    if n == 0 {
        return vec![];
    }
    let (base, rem) = (work_items / n, work_items % n);
    let mut begin = 0;
    group
        .device_ordinals
        .iter()
        .enumerate()
        .map(|(i, d)| {
            let count = base + usize::from(i < rem);
            let s = DeviceWorkSlice {
                device_ordinal: *d,
                begin,
                end: begin + count,
            };
            begin += count;
            s
        })
        .collect()
}
