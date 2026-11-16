use crate::{core::Vec3, positive, require, Result, SolverError};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct PolyCell {
    pub center: Vec3,
    pub volume: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Face {
    pub owner: usize,
    pub neighbour: Option<usize>,
    pub center: Vec3,
    pub area: Vec3,
    pub patch: Option<usize>,
}
impl Face {
    pub fn boundary(self) -> bool {
        self.neighbour.is_none()
    }
    pub fn internal(owner: usize, neighbour: usize, center: Vec3, area: Vec3) -> Self {
        Self {
            owner,
            neighbour: Some(neighbour),
            center,
            area,
            patch: None,
        }
    }
    pub fn boundary_face(owner: usize, center: Vec3, area: Vec3, patch: usize) -> Self {
        Self {
            owner,
            neighbour: None,
            center,
            area,
            patch: Some(patch),
        }
    }
    pub fn outward_area(self, cell: usize) -> Vec3 {
        if self.owner == cell {
            self.area
        } else {
            -self.area
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BoundaryPatch {
    pub name: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(try_from = "RawPolyMesh")]
pub struct PolyMesh {
    cells: Vec<PolyCell>,
    faces: Vec<Face>,
    patches: Vec<BoundaryPatch>,
    #[serde(skip)]
    cell_faces: Vec<Vec<usize>>,
}
#[derive(Deserialize)]
struct RawPolyMesh {
    cells: Vec<PolyCell>,
    faces: Vec<Face>,
    patches: Vec<BoundaryPatch>,
}
impl TryFrom<RawPolyMesh> for PolyMesh {
    type Error = SolverError;
    fn try_from(r: RawPolyMesh) -> Result<Self> {
        Self::new(r.cells, r.faces, r.patches)
    }
}
impl PolyMesh {
    pub fn new(
        cells: Vec<PolyCell>,
        faces: Vec<Face>,
        patches: Vec<BoundaryPatch>,
    ) -> Result<Self> {
        require(!cells.is_empty(), "poly mesh needs cells")?;
        for c in &cells {
            positive(c.volume, "cell volume")?;
            require(c.center.is_finite(), "nonfinite cell centre")?;
        }
        let mut cell_faces = vec![vec![]; cells.len()];
        for (i, f) in faces.iter().enumerate() {
            require(
                f.owner < cells.len()
                    && f.center.is_finite()
                    && f.area.is_finite()
                    && f.area.norm() > 0.0,
                "invalid poly face geometry/owner",
            )?;
            cell_faces[f.owner].push(i);
            if let Some(n) = f.neighbour {
                require(
                    n < cells.len() && n != f.owner && f.patch.is_none(),
                    "invalid internal neighbour/patch",
                )?;
                cell_faces[n].push(i);
            } else {
                require(
                    f.patch.is_some_and(|p| p < patches.len()),
                    "boundary face patch missing/out of range",
                )?;
            }
        }
        Ok(Self {
            cells,
            faces,
            patches,
            cell_faces,
        })
    }
    pub fn cells(&self) -> &[PolyCell] {
        &self.cells
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn patches(&self) -> &[BoundaryPatch] {
        &self.patches
    }
    pub fn cell_faces(&self) -> &[Vec<usize>] {
        &self.cell_faces
    }
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }
    pub fn boundary_face_count(&self) -> usize {
        self.faces.iter().filter(|f| f.boundary()).count()
    }
    pub fn patch_index(&self, name: &str) -> Result<usize> {
        self.patches
            .iter()
            .position(|p| p.name == name)
            .ok_or_else(|| SolverError::InvalidInput(format!("unknown boundary patch {name}")))
    }
}
pub fn make_cartesian_hexa_mesh(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
) -> Result<PolyMesh> {
    require(
        nx > 0 && ny > 0 && nz > 0,
        "positive poly mesh subdivisions required",
    )?;
    positive(lx, "lx")?;
    positive(ly, "ly")?;
    positive(lz, "lz")?;
    let count = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(nz))
        .ok_or_else(|| SolverError::InvalidInput("poly mesh size overflow".into()))?;
    let face_count = || {
        let a = nx.checked_add(1)?.checked_mul(ny)?.checked_mul(nz)?;
        let b = ny.checked_add(1)?.checked_mul(nx)?.checked_mul(nz)?;
        let c = nz.checked_add(1)?.checked_mul(nx)?.checked_mul(ny)?;
        a.checked_add(b)?.checked_add(c)
    };
    let face_count =
        face_count().ok_or_else(|| SolverError::InvalidInput("face count overflow".into()))?;
    let (dx, dy, dz) = (lx / nx as f64, ly / ny as f64, lz / nz as f64);
    let volume = dx * dy * dz;
    positive(volume, "poly cell volume")?;
    let id = |i, j, k| (k * ny + j) * nx + i;
    let mut cells = Vec::with_capacity(count);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                cells.push(PolyCell {
                    center: Vec3::new(
                        (i as f64 + 0.5) * dx,
                        (j as f64 + 0.5) * dy,
                        (k as f64 + 0.5) * dz,
                    ),
                    volume,
                });
            }
        }
    }
    let patches = ["left", "right", "bottom", "top", "front", "back"]
        .into_iter()
        .map(|name| BoundaryPatch { name: name.into() })
        .collect();
    let mut faces = Vec::with_capacity(face_count);
    for k in 0..nz {
        for j in 0..ny {
            let (y, z) = ((j as f64 + 0.5) * dy, (k as f64 + 0.5) * dz);
            faces.push(Face::boundary_face(
                id(0, j, k),
                Vec3::new(0.0, y, z),
                Vec3::new(-dy * dz, 0.0, 0.0),
                0,
            ));
            for i in 1..nx {
                faces.push(Face::internal(
                    id(i - 1, j, k),
                    id(i, j, k),
                    Vec3::new(i as f64 * dx, y, z),
                    Vec3::new(dy * dz, 0.0, 0.0),
                ));
            }
            faces.push(Face::boundary_face(
                id(nx - 1, j, k),
                Vec3::new(lx, y, z),
                Vec3::new(dy * dz, 0.0, 0.0),
                1,
            ));
        }
    }
    for k in 0..nz {
        for i in 0..nx {
            let (x, z) = ((i as f64 + 0.5) * dx, (k as f64 + 0.5) * dz);
            faces.push(Face::boundary_face(
                id(i, 0, k),
                Vec3::new(x, 0.0, z),
                Vec3::new(0.0, -dx * dz, 0.0),
                2,
            ));
            for j in 1..ny {
                faces.push(Face::internal(
                    id(i, j - 1, k),
                    id(i, j, k),
                    Vec3::new(x, j as f64 * dy, z),
                    Vec3::new(0.0, dx * dz, 0.0),
                ));
            }
            faces.push(Face::boundary_face(
                id(i, ny - 1, k),
                Vec3::new(x, ly, z),
                Vec3::new(0.0, dx * dz, 0.0),
                3,
            ));
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            let (x, y) = ((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy);
            faces.push(Face::boundary_face(
                id(i, j, 0),
                Vec3::new(x, y, 0.0),
                Vec3::new(0.0, 0.0, -dx * dy),
                4,
            ));
            for k in 1..nz {
                faces.push(Face::internal(
                    id(i, j, k - 1),
                    id(i, j, k),
                    Vec3::new(x, y, k as f64 * dz),
                    Vec3::new(0.0, 0.0, dx * dy),
                ));
            }
            faces.push(Face::boundary_face(
                id(i, j, nz - 1),
                Vec3::new(x, y, lz),
                Vec3::new(0.0, 0.0, dx * dy),
                5,
            ));
        }
    }
    PolyMesh::new(cells, faces, patches)
}
pub fn make_sheared_cartesian_hexa_mesh(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
    shear_xy: f64,
) -> Result<PolyMesh> {
    require(shear_xy.is_finite(), "nonfinite mesh shear")?;
    let b = make_cartesian_hexa_mesh(nx, ny, nz, lx, ly, lz)?;
    let p = |v: Vec3| Vec3::new(v.x + shear_xy * v.y, v.y, v.z);
    let cells = b
        .cells
        .into_iter()
        .map(|mut c| {
            c.center = p(c.center);
            c
        })
        .collect();
    let faces = b
        .faces
        .into_iter()
        .map(|mut f| {
            f.center = p(f.center);
            f.area = Vec3::new(f.area.x, f.area.y - shear_xy * f.area.x, f.area.z);
            f
        })
        .collect();
    PolyMesh::new(cells, faces, b.patches)
}
