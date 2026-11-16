//! Validated ASCII VTU/PVD and CSV, OBJ and ASCII STL triangle-surface readers, binary
//! checkpoints, OpenFOAM polyMesh/field I/O and the XDMF mesh contract. HDF5 is not linked.
mod formats;
use crate::{
    core::Vec3,
    fem::{Mesh2D, Mesh3D},
    finite_slice, positive, require, Result, SolverError,
};
pub use formats::*;
use std::{
    io::{BufWriter, Write},
    path::Path,
};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CellType {
    Vertex,
    Line,
    Triangle,
    Quad,
    Tetrahedron,
    Hexahedron,
    Wedge,
    Pyramid,
}
impl CellType {
    pub fn vtk_id(self) -> u8 {
        match self {
            Self::Vertex => 1,
            Self::Line => 3,
            Self::Triangle => 5,
            Self::Quad => 9,
            Self::Tetrahedron => 10,
            Self::Hexahedron => 12,
            Self::Wedge => 13,
            Self::Pyramid => 14,
        }
    }
    pub fn nodes(self) -> usize {
        match self {
            Self::Vertex => 1,
            Self::Line => 2,
            Self::Triangle => 3,
            Self::Quad | Self::Tetrahedron => 4,
            Self::Hexahedron => 8,
            Self::Wedge => 6,
            Self::Pyramid => 5,
        }
    }
}
#[derive(Clone, Debug)]
pub struct Cell {
    pub cell_type: CellType,
    pub nodes: Vec<usize>,
}
#[derive(Clone, Debug)]
pub struct VtuMesh {
    pub points: Vec<Vec3>,
    pub cells: Vec<Cell>,
}
impl VtuMesh {
    pub fn from_triangles(m: &Mesh2D) -> Result<Self> {
        m.validate()?;
        Ok(Self {
            points: m.nodes.iter().map(|p| p.point3()).collect(),
            cells: m
                .triangles
                .iter()
                .map(|t| Cell {
                    cell_type: CellType::Triangle,
                    nodes: t.node.to_vec(),
                })
                .collect(),
        })
    }
    pub fn from_tetrahedra(m: &Mesh3D) -> Result<Self> {
        m.validate()?;
        Ok(Self {
            points: m.nodes.clone(),
            cells: m
                .tetrahedra
                .iter()
                .map(|t| Cell {
                    cell_type: CellType::Tetrahedron,
                    nodes: t.node.to_vec(),
                })
                .collect(),
        })
    }
    pub fn line(nodes: usize, spacing: f64, offset: f64) -> Result<Self> {
        require(nodes >= 2 && offset.is_finite(), "line mesh size/offset")?;
        positive(spacing, "line spacing")?;
        Ok(Self {
            points: (0..nodes)
                .map(|i| Vec3::new(offset + i as f64 * spacing, 0.0, 0.0))
                .collect(),
            cells: (0..nodes - 1)
                .map(|i| Cell {
                    cell_type: CellType::Line,
                    nodes: vec![i, i + 1],
                })
                .collect(),
        })
    }
    /// Corner points and quad cells; suitable for cell-centred lattice fields.
    pub fn rectangle_cells(nx: usize, ny: usize, dx: f64, dy: f64) -> Result<Self> {
        require(nx > 0 && ny > 0, "rectangle mesh dimensions")?;
        positive(dx, "dx")?;
        positive(dy, "dy")?;
        let count = nx
            .checked_add(1)
            .and_then(|x| ny.checked_add(1).and_then(|y| x.checked_mul(y)))
            .ok_or_else(|| crate::SolverError::InvalidInput("mesh size overflow".into()))?;
        let mut points = Vec::with_capacity(count);
        for j in 0..=ny {
            for i in 0..=nx {
                points.push(Vec3::new(i as f64 * dx, j as f64 * dy, 0.0));
            }
        }
        let mut cells = vec![];
        for j in 0..ny {
            for i in 0..nx {
                let a = j * (nx + 1) + i;
                cells.push(Cell {
                    cell_type: CellType::Quad,
                    nodes: vec![a, a + 1, a + nx + 2, a + nx + 1],
                });
            }
        }
        Ok(Self { points, cells })
    }
    /// Real corner connectivity in the same x-fastest cell order as PolyMesh.
    /// The affine transform is x' = x + shear*y and has unit determinant.
    pub fn box_cells(dimensions: [usize; 3], lengths: Vec3, shear: f64) -> Result<Self> {
        let [nx, ny, nz] = dimensions;
        require(
            nx > 0 && ny > 0 && nz > 0 && shear.is_finite(),
            "box dimensions/shear",
        )?;
        positive(lengths.x, "box width")?;
        positive(lengths.y, "box height")?;
        positive(lengths.z, "box depth")?;
        let count = nx
            .checked_add(1)
            .and_then(|x| ny.checked_add(1).and_then(|y| x.checked_mul(y)))
            .and_then(|xy| nz.checked_add(1).and_then(|z| xy.checked_mul(z)))
            .ok_or_else(|| crate::SolverError::InvalidInput("box size overflow".into()))?;
        let mut points = Vec::with_capacity(count);
        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    let y = j as f64 * lengths.y / ny as f64;
                    points.push(Vec3::new(
                        i as f64 * lengths.x / nx as f64 + shear * y,
                        y,
                        k as f64 * lengths.z / nz as f64,
                    ));
                }
            }
        }
        let id = |i, j, k| (k * (ny + 1) + j) * (nx + 1) + i;
        let mut cells = vec![];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    cells.push(Cell {
                        cell_type: CellType::Hexahedron,
                        nodes: vec![
                            id(i, j, k),
                            id(i + 1, j, k),
                            id(i + 1, j + 1, k),
                            id(i, j + 1, k),
                            id(i, j, k + 1),
                            id(i + 1, j, k + 1),
                            id(i + 1, j + 1, k + 1),
                            id(i, j + 1, k + 1),
                        ],
                    });
                }
            }
        }
        let m = Self { points, cells };
        m.validate()?;
        Ok(m)
    }
    pub fn validate(&self) -> Result<()> {
        require(
            !self.points.is_empty() && self.points.iter().all(|p| p.is_finite()),
            "invalid VTU points",
        )?;
        for c in &self.cells {
            require(
                c.nodes.len() == c.cell_type.nodes()
                    && c.nodes.iter().all(|&n| n < self.points.len()),
                "invalid VTU cell connectivity",
            )?;
        }
        Ok(())
    }
}
#[derive(Clone, Debug)]
pub struct DataArray {
    pub name: String,
    pub components: usize,
    pub values: Vec<f64>,
}
impl DataArray {
    pub fn scalar(name: impl Into<String>, values: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            components: 1,
            values,
        }
    }
    pub fn vector(name: impl Into<String>, values: Vec<Vec3>) -> Self {
        Self {
            name: name.into(),
            components: 3,
            values: values.into_iter().flat_map(|p| p.as_array()).collect(),
        }
    }
    fn validate(&self, count: usize) -> Result<()> {
        require(
            !self.name.is_empty()
                && !self
                    .name
                    .chars()
                    .any(|c| c < ' ' && c != '\t' && c != '\n' && c != '\r')
                && self.components > 0
                && count.checked_mul(self.components) == Some(self.values.len()),
            "VTU field name/components/length",
        )?;
        finite_slice(&self.values, &self.name)
    }
}
fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}
fn arrays(w: &mut impl Write, tag: &str, data: &[DataArray]) -> Result<()> {
    writeln!(w, "<{tag}>")?;
    for a in data {
        writeln!(
            w,
            "<DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"{}\" format=\"ascii\">",
            xml(&a.name),
            a.components
        )?;
        for v in &a.values {
            write!(w, "{v:.17e} ")?;
        }
        writeln!(w, "\n</DataArray>")?;
    }
    writeln!(w, "</{tag}>")?;
    Ok(())
}
pub fn write_vtu(
    path: impl AsRef<Path>,
    mesh: &VtuMesh,
    point_data: &[DataArray],
    cell_data: &[DataArray],
) -> Result<()> {
    mesh.validate()?;
    for a in point_data {
        a.validate(mesh.points.len())?;
    }
    for a in cell_data {
        a.validate(mesh.cells.len())?;
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(w,"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n<UnstructuredGrid><Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",mesh.points.len(),mesh.cells.len())?;
    arrays(&mut w, "PointData", point_data)?;
    arrays(&mut w, "CellData", cell_data)?;
    writeln!(
        w,
        "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
    )?;
    for p in &mesh.points {
        writeln!(w, "{:.17e} {:.17e} {:.17e}", p.x, p.y, p.z)?;
    }
    writeln!(w,"</DataArray></Points>\n<Cells><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">")?;
    for c in &mesh.cells {
        for n in &c.nodes {
            write!(w, "{n} ")?;
        }
        writeln!(w)?;
    }
    writeln!(
        w,
        "</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
    )?;
    let mut offset = 0;
    for c in &mesh.cells {
        offset += c.nodes.len();
        write!(w, "{offset} ")?;
    }
    writeln!(
        w,
        "\n</DataArray><DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    )?;
    for c in &mesh.cells {
        write!(w, "{} ", c.cell_type.vtk_id())?;
    }
    writeln!(
        w,
        "\n</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>"
    )?;
    w.flush()?;
    Ok(())
}
pub fn write_pvd(path: impl AsRef<Path>, frames: &[(f64, String)]) -> Result<()> {
    for (time, file) in frames {
        require(
            time.is_finite() && !file.is_empty() && !file.chars().any(|c| c < ' '),
            "invalid PVD frame",
        )?;
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(w,"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\"><Collection>")?;
    for (time, file) in frames {
        writeln!(
            w,
            "<DataSet timestep=\"{time:.17e}\" group=\"\" part=\"0\" file=\"{}\"/>",
            xml(file)
        )?;
    }
    writeln!(w, "</Collection></VTKFile>")?;
    w.flush()?;
    Ok(())
}
fn csv_quote(s: &str) -> String {
    if s.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}
pub fn write_csv(path: impl AsRef<Path>, columns: &[(&str, &[f64])]) -> Result<()> {
    require(!columns.is_empty(), "CSV needs columns")?;
    let rows = columns[0].1.len();
    for (_, v) in columns {
        require(v.len() == rows, "CSV column length mismatch")?;
        finite_slice(v, "CSV values")?;
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "{}",
        columns
            .iter()
            .map(|c| csv_quote(c.0))
            .collect::<Vec<_>>()
            .join(",")
    )?;
    for i in 0..rows {
        for (j, (_, v)) in columns.iter().enumerate() {
            if j > 0 {
                write!(w, ",")?;
            }
            write!(w, "{:.17e}", v[i])?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}

/// A triangle in 3-D, as upstream `cfd::io::Triangle3`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Triangle3 {
    pub a: Vec3,
    pub b: Vec3,
    pub c: Vec3,
}

/// Triangle soup, as upstream `cfd::io::TriangleSurface`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TriangleSurface {
    pub triangles: Vec<Triangle3>,
}

fn parse_f64(token: Option<&str>, what: &str) -> Result<f64> {
    token
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| SolverError::InvalidInput(format!("malformed {what} coordinate")))
}

/// Wavefront OBJ: `v x y z` vertices and the first three indices of each `f` record (1-based,
/// `i/j/k` forms accepted), ported from upstream `read_obj_triangles` (MIT, cfd_solvers contributors).
pub fn read_obj_triangles(path: impl AsRef<Path>) -> Result<TriangleSurface> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| SolverError::Io(format!("cannot open OBJ: {e}")))?;
    let mut vertices = Vec::new();
    let mut surface = TriangleSurface::default();
    for line in text.lines() {
        let mut tokens = line.split_whitespace();
        match tokens.next() {
            Some("v") => {
                let x = parse_f64(tokens.next(), "OBJ")?;
                let y = parse_f64(tokens.next(), "OBJ")?;
                let z = parse_f64(tokens.next(), "OBJ")?;
                vertices.push(Vec3::new(x, y, z));
            }
            Some("f") => {
                let mut corner = [Vec3::default(); 3];
                for c in &mut corner {
                    let token = tokens.next().unwrap_or("");
                    let id: i64 = token
                        .split('/')
                        .next()
                        .unwrap_or("")
                        .parse()
                        .map_err(|_| SolverError::InvalidInput("malformed OBJ face".into()))?;
                    if id <= 0 || id as usize > vertices.len() {
                        return Err(SolverError::InvalidInput("OBJ index out of range".into()));
                    }
                    *c = vertices[id as usize - 1];
                }
                surface.triangles.push(Triangle3 {
                    a: corner[0],
                    b: corner[1],
                    c: corner[2],
                });
            }
            _ => {}
        }
    }
    Ok(surface)
}

/// ASCII STL: every three `vertex` records form a triangle, ported from upstream
/// `read_ascii_stl` (MIT, cfd_solvers contributors).
pub fn read_ascii_stl(path: impl AsRef<Path>) -> Result<TriangleSurface> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| SolverError::Io(format!("cannot open STL: {e}")))?;
    let mut surface = TriangleSurface::default();
    let mut pending = Vec::with_capacity(3);
    let mut tokens = text.split_whitespace();
    while let Some(word) = tokens.next() {
        if word == "vertex" {
            let x = parse_f64(tokens.next(), "STL")?;
            let y = parse_f64(tokens.next(), "STL")?;
            let z = parse_f64(tokens.next(), "STL")?;
            pending.push(Vec3::new(x, y, z));
            if pending.len() == 3 {
                surface.triangles.push(Triangle3 {
                    a: pending[0],
                    b: pending[1],
                    c: pending[2],
                });
                pending.clear();
            }
        }
    }
    require(pending.is_empty(), "malformed ASCII STL")?;
    Ok(surface)
}

/// Inline-XML XDMF 3 writer for a triangle mesh with an optional nodal scalar; byte-identical
/// to upstream `cfd::io::write_xdmf_inline` (`%.17g` numbers).
pub fn write_xdmf_inline(
    path: impl AsRef<Path>,
    mesh: &Mesh2D,
    nodal_scalar: &[f64],
    scalar_name: &str,
) -> Result<()> {
    use crate::workflow::text::format_g17 as g;
    use std::fmt::Write as _;
    mesh.validate()?;
    require(
        nodal_scalar.is_empty() || nodal_scalar.len() == mesh.node_count(),
        "XDMF scalar size mismatch",
    )?;
    let mut out = String::from("<?xml version=\"1.0\" ?>\n<Xdmf Version=\"3.0\">\n <Domain>\n  <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    let ne = mesh.element_count();
    let _ = write!(
        out,
        "   <Topology TopologyType=\"Triangle\" NumberOfElements=\"{ne}\">\n    <DataItem Dimensions=\"{ne} 3\" NumberType=\"UInt\" Precision=\"8\" Format=\"XML\">\n"
    );
    for t in &mesh.triangles {
        let _ = writeln!(out, "{} {} {}", t.node[0], t.node[1], t.node[2]);
    }
    let _ = write!(
        out,
        "    </DataItem>\n   </Topology>\n   <Geometry GeometryType=\"XY\">\n    <DataItem Dimensions=\"{} 2\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n",
        mesh.node_count()
    );
    for n in &mesh.nodes {
        let _ = writeln!(out, "{} {}", g(n.x), g(n.y));
    }
    out.push_str("    </DataItem>\n   </Geometry>\n");
    if !nodal_scalar.is_empty() {
        let _ = write!(
            out,
            "   <Attribute Name=\"{scalar_name}\" AttributeType=\"Scalar\" Center=\"Node\">\n    <DataItem Dimensions=\"{}\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n",
            nodal_scalar.len()
        );
        for v in nodal_scalar {
            out.push_str(&g(*v));
            out.push(' ');
        }
        out.push_str("\n    </DataItem>\n   </Attribute>\n");
    }
    out.push_str("  </Grid>\n </Domain>\n</Xdmf>\n");
    std::fs::write(path, out).map_err(|e| SolverError::Io(format!("cannot create XDMF: {e}")))
}
