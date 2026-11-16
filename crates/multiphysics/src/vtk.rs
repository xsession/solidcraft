//! Bounded inspection of ASCII VTK XML companion output.
//! No field units or physical meaning are inferred. This is not a general VTK reader.
use quick_xml::{
    Reader,
    events::{BytesStart, Event},
};
use std::collections::BTreeSet;

const MAX_XML: usize = 8 * 1024 * 1024;
const MAX_POINTS: usize = 1_000_000;
const MAX_CELLS: usize = 1_000_000;
const MAX_VALUES: usize = 2_000_000;
const MAX_ARRAYS: usize = 128;
const MAX_FRAMES: usize = 4096;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Association {
    Point,
    Cell,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    pub name: String,
    pub association: Association,
    pub data_type: String,
    pub components: usize,
    /// Producer values in their original component order; units are unspecified.
    pub values: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Grid {
    pub points: Vec<[f64; 3]>,
    pub connectivity: Vec<usize>,
    pub offsets: Vec<usize>,
    pub cell_types: Vec<u8>,
    pub fields: Vec<Field>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Frame {
    pub coordinate: f64,
    pub file: String,
    pub part: usize,
}

/// A bounded preview representation, not a general hexahedral slicing algorithm.
/// Only canonical, axis-aligned VTK HEX8 cells are admitted. Original cell values
/// and ordering are retained; no point interpolation or physical units are inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct BoxGrid {
    pub bounds: [[f64; 2]; 3],
    pub cells: Vec<[[f64; 2]; 3]>,
    pub planes: [Vec<f64>; 3],
    pub fields: Vec<Field>,
}

/// One original cell's section polygon, in unchanged producer coordinates.
/// Corners wind toward the positive slice normal; no values are interpolated.
#[derive(Clone, Debug, PartialEq)]
pub struct SliceCell3d {
    pub cell: usize,
    pub corners: [[f64; 3]; 4],
    pub center: [f64; 3],
    pub value: f64,
    /// Original vector tuple, when the selected field has three components.
    pub vector: Option<[f64; 3]>,
}

impl BoxGrid {
    pub fn from_grid(grid: Grid) -> Result<Self, Error> {
        if grid.offsets.is_empty() || grid.offsets.len() > 32_768 || grid.cell_types.len() != grid.offsets.len() {
            return Err(bad("slice preview requires 1–32768 axis-aligned HEX8 cells"));
        }
        let mut cells = Vec::with_capacity(grid.offsets.len());
        let mut bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
        let mut planes: [Vec<f64>; 3] = Default::default();
        let mut start = 0;
        for (index, &end) in grid.offsets.iter().enumerate() {
            if grid.cell_types.get(index) != Some(&12) || end.checked_sub(start) != Some(8) {
                return Err(bad("slice preview supports canonical axis-aligned HEX8 cells only"));
            }
            let indices = grid.connectivity.get(start..end).ok_or_else(|| bad("invalid slice cell connectivity"))?;
            let mut cell = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
            for &vertex in indices {
                let point = grid.points.get(vertex).ok_or_else(|| bad("invalid slice point index"))?;
                for (axis, &value) in point.iter().enumerate() {
                    if !value.is_finite() {
                        return Err(bad("non-finite slice geometry"));
                    }
                    if let Some(range) = cell.get_mut(axis) {
                        range[0] = range[0].min(value);
                        range[1] = range[1].max(value);
                    }
                }
            }
            for axis in 0..3 {
                let Some(&[low, high]) = cell.get(axis) else { return Err(bad("invalid slice axis")) };
                if high <= low || !(high - low).is_finite() {
                    return Err(bad("degenerate or oversized slice cell"));
                }
                if let Some(range) = bounds.get_mut(axis) {
                    range[0] = range[0].min(low);
                    range[1] = range[1].max(high);
                }
                if let Some(levels) = planes.get_mut(axis) {
                    levels.push(low + (high - low) * 0.5);
                }
            }
            // Requiring the VTK corner order also rejects twisted connectivity.
            for (&vertex, mask) in indices.iter().zip([0_u8, 1, 3, 2, 4, 5, 7, 6]) {
                let point = grid.points.get(vertex).ok_or_else(|| bad("invalid slice point index"))?;
                for (axis, value) in point.iter().enumerate() {
                    let range = cell.get(axis).ok_or_else(|| bad("invalid slice axis"))?;
                    let expected = range[usize::from(mask & (1 << axis) != 0)];
                    if *value != expected {
                        return Err(bad("slice preview does not support sheared or noncanonical HEX8 cells"));
                    }
                }
            }
            cells.push(cell);
            start = end;
        }
        if start != grid.connectivity.len() {
            return Err(bad("unused slice connectivity"));
        }
        for axis in 0..3 {
            let range = bounds.get(axis).ok_or_else(|| bad("invalid slice axis"))?;
            if !(range[1] - range[0]).is_finite() {
                return Err(bad("slice extent exceeds numeric range"));
            }
        }
        for levels in &mut planes {
            levels.sort_by(f64::total_cmp);
            levels.dedup();
        }
        let fields = grid.fields.into_iter().filter(|field| field.association == Association::Cell).collect::<Vec<_>>();
        for field in &fields {
            if !matches!(field.components, 1 | 3)
                || field.components.checked_mul(cells.len()) != Some(field.values.len())
                || field.values.iter().any(|value| !value.is_finite())
            {
                return Err(bad("slice fields require finite scalar or three-component cell tuples"));
            }
        }
        if fields.is_empty() {
            return Err(bad("slice preview requires cell-associated scalar or vector fields"));
        }
        Ok(Self { bounds, cells, planes, fields })
    }

    /// Scalar tuple or Euclidean vector magnitude; original arrays are unchanged.
    pub fn value(&self, field: usize, cell: usize) -> Option<f64> {
        let field = self.fields.get(field)?;
        let start = cell.checked_mul(field.components)?;
        if field.components == 1 {
            field.values.get(start).copied()
        } else {
            let tuple = field.values.get(start..start.checked_add(3)?)?;
            let value = tuple.iter().fold(0.0_f64, |length, component| length.hypot(*component));
            value.is_finite().then_some(value)
        }
    }

    /// Half-open intervals avoid drawing two layers when the plane lies on a face.
    pub fn intersects(&self, cell: usize, normal: usize, coordinate: f64) -> bool {
        let Some(range) = self.cells.get(cell).and_then(|cell| cell.get(normal)) else { return false };
        coordinate.is_finite() && coordinate >= range[0] && coordinate < range[1]
    }

    /// Bounded, axis-aligned cell sections. A plane on a cell's upper face excludes
    /// that cell, matching [`Self::intersects`]. Empty sections remain empty.
    /// Tangential axes are (Y,Z), (Z,X), (X,Y) for X, Y, Z normals respectively.
    pub fn slice3d(&self, field: usize, normal: usize, coordinate: f64) -> Result<Vec<SliceCell3d>, Error> {
        let (u, v) = match normal {
            0 => (1, 2),
            1 => (2, 0),
            2 => (0, 1),
            _ => return Err(bad("3D slice normal must be X, Y or Z")),
        };
        if !coordinate.is_finite() {
            return Err(bad("3D slice coordinate must be finite"));
        }
        if self.cells.is_empty() || self.cells.len() > 32_768 {
            return Err(bad("3D slice requires 1–32768 cells"));
        }
        let selected = self.fields.get(field).ok_or_else(|| bad("3D slice field is unavailable"))?;
        if selected.association != Association::Cell
            || !matches!(selected.components, 1 | 3)
            || self.cells.len().checked_mul(selected.components) != Some(selected.values.len())
            || selected.values.iter().any(|value| !value.is_finite())
        {
            return Err(bad("3D slice requires finite scalar or vector cell tuples"));
        }
        let mut section = Vec::new();
        for (index, cell) in self.cells.iter().enumerate() {
            if cell.iter().any(|[low, high]| !low.is_finite() || !high.is_finite() || high <= low || !(high - low).is_finite()) {
                return Err(bad("3D slice has invalid cell bounds"));
            }
            if !self.intersects(index, normal, coordinate) {
                continue;
            }
            let mut center = cell.map(|[low, high]| low + (high - low) * 0.5);
            *center.get_mut(normal).ok_or_else(|| bad("invalid slice normal"))? = coordinate;
            let mut corners = [center; 4];
            for (corner, (upper_u, upper_v)) in corners.iter_mut().zip([(false, false), (true, false), (true, true), (false, true)]) {
                *corner.get_mut(u).ok_or_else(|| bad("invalid tangent axis"))? =
                    *cell.get(u).and_then(|range| range.get(usize::from(upper_u))).ok_or_else(|| bad("invalid tangent bounds"))?;
                *corner.get_mut(v).ok_or_else(|| bad("invalid tangent axis"))? =
                    *cell.get(v).and_then(|range| range.get(usize::from(upper_v))).ok_or_else(|| bad("invalid tangent bounds"))?;
            }
            let value = self.value(field, index).filter(|value| value.is_finite()).ok_or_else(|| bad("3D slice magnitude exceeds numeric range"))?;
            let vector = if selected.components == 3 {
                let start = index.checked_mul(3).ok_or_else(|| bad("3D slice tuple offset overflow"))?;
                let tuple = selected
                    .values
                    .get(start..start.checked_add(3).ok_or_else(|| bad("3D slice tuple offset overflow"))?)
                    .ok_or_else(|| bad("3D slice vector tuple is unavailable"))?;
                Some(<[f64; 3]>::try_from(tuple).map_err(|_| bad("invalid 3D slice vector tuple"))?)
            } else {
                None
            };
            section.push(SliceCell3d { cell: index, corners, center, value, vector });
        }
        Ok(section)
    }
}

#[derive(Debug, thiserror::Error)]
#[error("VTK XML: {0}")]
pub struct Error(pub String);

fn bad(message: impl Into<String>) -> Error {
    Error(message.into())
}

fn attr(element: &BytesStart<'_>, name: &[u8]) -> Result<Option<String>, Error> {
    let mut found = None;
    for attribute in element.attributes() {
        let attribute = attribute.map_err(|e| bad(e.to_string()))?;
        if attribute.key.as_ref() == name {
            if found.is_some() {
                return Err(bad("duplicate attribute"));
            }
            let value = std::str::from_utf8(attribute.value.as_ref()).map_err(|e| bad(e.to_string()))?;
            if value.contains('&') || value.len() > 512 {
                return Err(bad("escaped or oversized attributes are unsupported"));
            }
            found = Some(value.to_string());
        }
    }
    Ok(found)
}

fn required(element: &BytesStart<'_>, name: &[u8]) -> Result<String, Error> {
    attr(element, name)?.ok_or_else(|| bad(format!("missing {}", String::from_utf8_lossy(name))))
}
fn natural(value: &str, label: &str, max: usize) -> Result<usize, Error> {
    let n: usize = value.parse().map_err(|_| bad(format!("invalid {label}")))?;
    if n > max {
        return Err(bad(format!("{label} exceeds limit")));
    }
    Ok(n)
}
fn checked_input(xml: &str) -> Result<(), Error> {
    if xml.len() > MAX_XML {
        return Err(bad("document exceeds 8 MiB"));
    }
    Ok(())
}

struct Array {
    parent: String,
    name: String,
    data_type: String,
    components: usize,
    text: String,
}
fn numeric(text: &str, data_type: &str) -> Result<Vec<f64>, Error> {
    let mut values = Vec::new();
    for word in text.split_whitespace() {
        if values.len() >= MAX_VALUES {
            return Err(bad("too many array values"));
        }
        let value: f64 = if data_type.starts_with("Int") || data_type.starts_with("UInt") {
            let n: i128 = word.parse().map_err(|_| bad("invalid integer array"))?;
            if n.unsigned_abs() > (1_u128 << 53) || !integer_fits(n, data_type) {
                return Err(bad("integer field cannot be represented exactly"));
            }
            n as f64
        } else if data_type == "Float32" {
            f64::from(word.parse::<f32>().map_err(|_| bad("invalid Float32 array"))?)
        } else {
            word.parse().map_err(|_| bad("invalid numeric array"))?
        };
        if !value.is_finite() {
            return Err(bad("non-finite array value"));
        }
        values.push(value);
    }
    Ok(values)
}
fn integer_fits(n: i128, kind: &str) -> bool {
    match kind {
        "Int8" => i8::try_from(n).is_ok(),
        "UInt8" => u8::try_from(n).is_ok(),
        "Int16" => i16::try_from(n).is_ok(),
        "UInt16" => u16::try_from(n).is_ok(),
        "Int32" => i32::try_from(n).is_ok(),
        "UInt32" => u32::try_from(n).is_ok(),
        "Int64" => i64::try_from(n).is_ok(),
        "UInt64" => u64::try_from(n).is_ok(),
        _ => false,
    }
}
fn integers(text: &str, kind: &str) -> Result<Vec<usize>, Error> {
    let mut values = Vec::new();
    for word in text.split_whitespace() {
        if values.len() >= MAX_VALUES {
            return Err(bad("too many integer values"));
        }
        let n: i128 = word.parse().map_err(|_| bad("invalid cell integer"))?;
        if !integer_fits(n, kind) {
            return Err(bad("cell integer outside its declared type"));
        }
        values.push(usize::try_from(n).map_err(|_| bad("cell integer outside platform range"))?);
    }
    Ok(values)
}

/// Decode one ASCII unstructured-grid piece. Unsupported encodings fail explicitly.
pub fn parse_vtu(xml: &str) -> Result<Grid, Error> {
    checked_input(xml)?;
    let mut reader = Reader::from_str(xml);
    let mut stack: Vec<String> = Vec::new();
    let mut kind = false;
    let mut piece = None;
    let mut current: Option<Array> = None;
    let mut points: Option<Vec<[f64; 3]>> = None;
    let mut connectivity = None;
    let mut offsets = None;
    let mut types = None;
    let mut fields = Vec::new();
    let mut names = BTreeSet::new();
    let mut total_values = 0_usize;
    loop {
        match reader.read_event().map_err(|e| bad(e.to_string()))? {
            Event::Start(e) => {
                let tag = std::str::from_utf8(e.name().as_ref()).map_err(|_| bad("invalid tag"))?.to_string();
                if stack.len() >= 32 {
                    return Err(bad("XML nesting exceeds limit"));
                }
                match tag.as_str() {
                    "VTKFile" => {
                        if !stack.is_empty() || kind || required(&e, b"type")? != "UnstructuredGrid" {
                            return Err(bad("expected one unstructured grid"));
                        }
                        kind = true;
                    }
                    "Piece" => {
                        if piece.is_some() || stack.last().map(String::as_str) != Some("UnstructuredGrid") {
                            return Err(bad("expected one grid piece"));
                        }
                        piece = Some((
                            natural(&required(&e, b"NumberOfPoints")?, "point count", MAX_POINTS)?,
                            natural(&required(&e, b"NumberOfCells")?, "cell count", MAX_CELLS)?,
                        ));
                    }
                    "DataArray" => {
                        if current.is_some() || fields.len() >= MAX_ARRAYS {
                            return Err(bad("nested or excessive arrays"));
                        }
                        if attr(&e, b"format")?.as_deref().unwrap_or("ascii") != "ascii" {
                            return Err(bad("only ASCII DataArray is supported"));
                        }
                        let components = natural(attr(&e, b"NumberOfComponents")?.as_deref().unwrap_or("1"), "components", 16)?;
                        if components == 0 {
                            return Err(bad("zero components"));
                        }
                        let data_type = attr(&e, b"type")?.unwrap_or_else(|| "Float64".into());
                        if !matches!(
                            data_type.as_str(),
                            "Float32" | "Float64" | "Int8" | "UInt8" | "Int16" | "UInt16" | "Int32" | "UInt32" | "Int64" | "UInt64"
                        ) {
                            return Err(bad("unsupported array type"));
                        }
                        current = Some(Array {
                            parent: stack.last().cloned().unwrap_or_default(),
                            name: attr(&e, b"Name")?.unwrap_or_default(),
                            data_type,
                            components,
                            text: String::new(),
                        });
                    }
                    "UnstructuredGrid" if stack.last().map(String::as_str) == Some("VTKFile") => {}
                    "PointData" | "CellData" | "Points" | "Cells" if stack.last().map(String::as_str) == Some("Piece") => {}
                    "AppendedData" | "FieldData" => return Err(bad("appended and FieldData arrays are unsupported")),
                    _ => return Err(bad(format!("unsupported or misplaced element {tag}"))),
                }
                stack.push(tag);
            }
            Event::Text(e) => {
                if let Some(array) = current.as_mut() {
                    let bytes = e.as_ref();
                    if bytes.contains(&b'&') {
                        return Err(bad("entity references in arrays are unsupported"));
                    }
                    array.text.push_str(std::str::from_utf8(bytes).map_err(|_| bad("array is not UTF-8"))?);
                }
            }
            Event::End(e) => {
                let tag = std::str::from_utf8(e.name().as_ref()).map_err(|_| bad("invalid tag"))?.to_string();
                if stack.pop().as_deref() != Some(tag.as_str()) {
                    return Err(bad("mismatched tags"));
                }
                if tag == "DataArray" {
                    let array = current.take().ok_or_else(|| bad("unexpected DataArray end"))?;
                    let (np, nc) = piece.ok_or_else(|| bad("array outside piece"))?;
                    match array.parent.as_str() {
                        "Points" => {
                            if points.is_some() || array.components != 3 {
                                return Err(bad("expected one three-component points array"));
                            }
                            let values = numeric(&array.text, &array.data_type)?;
                            total_values = total_values.saturating_add(values.len());
                            if values.len() != np.saturating_mul(3) {
                                return Err(bad("point count mismatch"));
                            }
                            points = Some(values.chunks_exact(3).map(|v| [v[0], v[1], v[2]]).collect());
                        }
                        "Cells" => {
                            if array.components != 1 {
                                return Err(bad("cell arrays must be scalar"));
                            }
                            let values = integers(&array.text, &array.data_type)?;
                            total_values = total_values.saturating_add(values.len());
                            let target = match array.name.as_str() {
                                "connectivity" => &mut connectivity,
                                "offsets" => &mut offsets,
                                "types" => &mut types,
                                _ => return Err(bad("unsupported cell array")),
                            };
                            if target.replace(values).is_some() {
                                return Err(bad("duplicate cell array"));
                            }
                        }
                        "PointData" | "CellData" => {
                            let association = if array.parent == "PointData" { Association::Point } else { Association::Cell };
                            if array.name.is_empty() || !names.insert((array.parent.clone(), array.name.clone())) {
                                return Err(bad("missing or duplicate field name"));
                            }
                            let values = numeric(&array.text, &array.data_type)?;
                            total_values = total_values.saturating_add(values.len());
                            let expected = if association == Association::Point { np } else { nc };
                            if values.len() != expected.saturating_mul(array.components) {
                                return Err(bad("field tuple count mismatch"));
                            }
                            fields.push(Field { name: array.name, association, data_type: array.data_type, components: array.components, values });
                        }
                        _ => return Err(bad("DataArray in unsupported location")),
                    }
                    if total_values > MAX_VALUES {
                        return Err(bad("total array values exceed limit"));
                    }
                }
            }
            Event::Empty(e) => {
                let tag = std::str::from_utf8(e.name().as_ref()).map_err(|_| bad("invalid tag"))?.to_string();
                if !matches!(tag.as_str(), "PointData" | "CellData") {
                    return Err(bad(format!("unsupported empty element {tag}")));
                }
            }
            Event::DocType(_) | Event::GeneralRef(_) | Event::CData(_) => return Err(bad("DTD, entities and CDATA are unsupported")),
            Event::Eof => break,
            Event::Decl(_) | Event::Comment(_) | Event::PI(_) => {}
        }
    }
    if !kind || !stack.is_empty() || current.is_some() {
        return Err(bad("incomplete XML"));
    }
    let (np, nc) = piece.ok_or_else(|| bad("missing piece"))?;
    let points = points.ok_or_else(|| bad("missing points"))?;
    let connectivity: Vec<usize> = connectivity.ok_or_else(|| bad("missing connectivity"))?;
    let offsets: Vec<usize> = offsets.ok_or_else(|| bad("missing offsets"))?;
    let raw_types: Vec<usize> = types.ok_or_else(|| bad("missing types"))?;
    if offsets.len() != nc
        || raw_types.len() != nc
        || offsets.last().copied().unwrap_or(0) != connectivity.len()
        || offsets.windows(2).any(|w| w[0] >= w[1])
        || offsets.first().is_some_and(|n| *n == 0)
        || connectivity.iter().any(|i| *i >= np)
    {
        return Err(bad("invalid cell connectivity or offsets"));
    }
    let cell_types = raw_types.into_iter().map(|t| u8::try_from(t).map_err(|_| bad("cell type out of range"))).collect::<Result<Vec<_>, _>>()?;
    Ok(Grid { points, connectivity, offsets, cell_types, fields })
}

/// Parse a PVD collection without opening paths. Coordinate meaning follows the producer.
pub fn parse_pvd(xml: &str) -> Result<Vec<Frame>, Error> {
    checked_input(xml)?;
    let mut reader = Reader::from_str(xml);
    let mut stack: Vec<String> = Vec::new();
    let mut collection = false;
    let mut collection_body = false;
    let mut frames = Vec::new();
    loop {
        match reader.read_event().map_err(|e| bad(e.to_string()))? {
            Event::Start(e) => {
                let tag = std::str::from_utf8(e.name().as_ref()).map_err(|_| bad("invalid tag"))?.to_string();
                if stack.len() >= 32 {
                    return Err(bad("XML nesting exceeds limit"));
                }
                match tag.as_str() {
                    "VTKFile" if stack.is_empty() && !collection && required(&e, b"type")? == "Collection" => collection = true,
                    "Collection" if stack.last().map(String::as_str) == Some("VTKFile") && !collection_body => collection_body = true,
                    _ => return Err(bad("unsupported or misplaced collection element")),
                }
                stack.push(tag);
            }
            Event::Empty(e) if e.name().as_ref() == b"DataSet" => {
                if stack.last().map(String::as_str) != Some("Collection") || frames.len() >= MAX_FRAMES {
                    return Err(bad("invalid or excessive datasets"));
                }
                let coordinate: f64 = required(&e, b"timestep")?.parse().map_err(|_| bad("invalid collection coordinate"))?;
                if !coordinate.is_finite() {
                    return Err(bad("non-finite collection coordinate"));
                }
                let file = required(&e, b"file")?;
                crate::validate_artifact_name(&file).map_err(|e| bad(e.to_string()))?;
                let part = natural(attr(&e, b"part")?.as_deref().unwrap_or("0"), "part", usize::MAX)?;
                frames.push(Frame { coordinate, file, part });
            }
            Event::End(e) => {
                let tag = std::str::from_utf8(e.name().as_ref()).map_err(|_| bad("invalid tag"))?.to_string();
                if stack.pop().as_deref() != Some(tag.as_str()) {
                    return Err(bad("mismatched tags"));
                }
            }
            Event::DocType(_) | Event::GeneralRef(_) | Event::CData(_) => return Err(bad("DTD, entities and CDATA are unsupported")),
            Event::Eof => break,
            Event::Empty(_) => return Err(bad("unsupported empty collection element")),
            Event::Text(_) | Event::Decl(_) | Event::Comment(_) | Event::PI(_) => {}
        }
    }
    if !collection || !collection_body || !stack.is_empty() || frames.is_empty() {
        return Err(bad("empty or incomplete collection"));
    }
    Ok(frames)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn boxes() -> Grid {
        let mut points = Vec::new();
        for x in [0.0, 2.0] {
            points.extend([
                [x, 0., 0.],
                [x + 1., 0., 0.],
                [x + 1., 1., 0.],
                [x, 1., 0.],
                [x, 0., 1.],
                [x + 1., 0., 1.],
                [x + 1., 1., 1.],
                [x, 1., 1.],
            ]);
        }
        Grid {
            points,
            connectivity: (0..16).collect(),
            offsets: vec![8, 16],
            cell_types: vec![12, 12],
            fields: vec![Field {
                name: "velocity".into(),
                association: Association::Cell,
                data_type: "Float64".into(),
                components: 3,
                values: vec![3., 4., 0., 0., 0., 2.],
            }],
        }
    }
    #[test]
    fn box_slice_preserves_obstacle_gap_and_cell_values() {
        let preview = BoxGrid::from_grid(boxes()).unwrap();
        assert_eq!(preview.bounds, [[0., 3.], [0., 1.], [0., 1.]]);
        assert_eq!(preview.planes[0], [0.5, 2.5]);
        assert_eq!(preview.cells.len(), 2, "the gap is not filled by display topology");
        assert_eq!(preview.value(0, 0), Some(5.));
        assert_eq!(preview.value(0, 1), Some(2.));
        assert!(preview.intersects(0, 2, 0.5));
        assert!(!preview.intersects(0, 0, 1.));
        assert!(!preview.intersects(0, 3, 0.5));
        assert_eq!(preview.value(0, 2), None);
    }
    #[test]
    fn box_slice_refuses_sheared_twisted_and_point_only_fields() {
        let mut grid = boxes();
        grid.points[2][0] += 0.1;
        assert!(BoxGrid::from_grid(grid).is_err());
        let mut grid = boxes();
        grid.connectivity.swap(0, 1);
        assert!(BoxGrid::from_grid(grid).is_err());
        let mut grid = boxes();
        grid.fields[0].association = Association::Point;
        assert!(BoxGrid::from_grid(grid).is_err());
        let mut grid = boxes();
        grid.offsets[0] = usize::MAX;
        assert!(BoxGrid::from_grid(grid).is_err());
    }
    #[test]
    fn slice3d_preserves_planes_winding_holes_and_original_tuples() {
        let preview = BoxGrid::from_grid(boxes()).unwrap();
        let xy = preview.slice3d(0, 2, 0.25).unwrap();
        assert_eq!(xy.len(), 2);
        assert_eq!(
            xy[0],
            SliceCell3d {
                cell: 0,
                corners: [[0., 0., 0.25], [1., 0., 0.25], [1., 1., 0.25], [0., 1., 0.25]],
                center: [0.5, 0.5, 0.25],
                value: 5.,
                vector: Some([3., 4., 0.]),
            }
        );
        assert_eq!(xy[1].cell, 1);
        assert_eq!(xy[1].center, [2.5, 0.5, 0.25]);
        assert_eq!(xy[1].value, 2.);
        assert!(xy.iter().all(|cell| cell.corners.iter().all(|point| point[0] <= 1. || point[0] >= 2.)));
        let xz = preview.slice3d(0, 1, 0.75).unwrap();
        assert_eq!(xz[0].corners, [[0., 0.75, 0.], [0., 0.75, 1.], [1., 0.75, 1.], [1., 0.75, 0.]]);
        let yz = preview.slice3d(0, 0, 2.25).unwrap();
        assert_eq!(yz.len(), 1);
        assert_eq!(yz[0].cell, 1);
        assert_eq!(yz[0].corners, [[2.25, 0., 0.], [2.25, 1., 0.], [2.25, 1., 1.], [2.25, 0., 1.]]);
        assert_eq!(yz[0].vector, Some([0., 0., 2.]));
        assert!(preview.slice3d(0, 0, 1.).unwrap().is_empty(), "upper faces and the hole remain excluded");
        assert_eq!(preview.slice3d(0, 0, 0.).unwrap().len(), 1, "lower face is included");
        assert!(preview.slice3d(0, 2, 1.).unwrap().is_empty());
        let mut scalar = preview;
        scalar.fields =
            vec![Field { name: "p".into(), association: Association::Cell, data_type: "Float64".into(), components: 1, values: vec![-2., 7.] }];
        let section = scalar.slice3d(0, 2, 0.5).unwrap();
        assert_eq!(section[0].value, -2.);
        assert_eq!(section[0].vector, None);
        assert_eq!(section[1].value, 7.);
    }
    #[test]
    fn slice3d_rejects_invalid_arguments_and_bounded_hostile_grids() {
        let mut preview = BoxGrid::from_grid(boxes()).unwrap();
        assert!(preview.slice3d(0, 3, 0.5).is_err());
        assert!(preview.slice3d(1, 2, 0.5).is_err());
        for coordinate in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(preview.slice3d(0, 2, coordinate).is_err());
        }
        preview.fields[0].values.pop();
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        preview = BoxGrid::from_grid(boxes()).unwrap();
        preview.fields[0].values[0] = f64::NAN;
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        preview.fields[0].values = vec![f64::MAX; 6];
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        preview = BoxGrid::from_grid(boxes()).unwrap();
        preview.fields[0].association = Association::Point;
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        preview.fields[0].association = Association::Cell;
        preview.fields[0].components = 2;
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        preview = BoxGrid::from_grid(boxes()).unwrap();
        preview.cells[0][0] = [1., 1.];
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        preview.cells[0][0] = [-f64::MAX, f64::MAX];
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        preview.cells.clear();
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        preview.cells.resize(32_769, [[0., 1.]; 3]);
        assert!(preview.slice3d(0, 2, 0.5).is_err());
        let maximum = BoxGrid {
            bounds: [[0., 32_768.], [0., 1.], [0., 1.]],
            cells: (0..32_768).map(|i| [[i as f64, i as f64 + 1.], [0., 1.], [0., 1.]]).collect(),
            planes: Default::default(),
            fields: vec![Field {
                name: "scalar".into(),
                association: Association::Cell,
                data_type: "Float64".into(),
                components: 1,
                values: vec![7.; 32_768],
            }],
        };
        let section = maximum.slice3d(0, 2, 0.5).unwrap();
        assert_eq!(section.len(), 32_768);
        assert_eq!(section.last().unwrap().cell, 32_767);
    }
    const LINE: &str = r#"<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian"><UnstructuredGrid><Piece NumberOfPoints="2" NumberOfCells="1"><PointData><DataArray type="Float64" Name="potential" NumberOfComponents="1" format="ascii">0 1</DataArray></PointData><CellData><DataArray type="Float64" Name="flux" NumberOfComponents="2" format="ascii">2 3</DataArray></CellData><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0</DataArray></Points><Cells><DataArray type="Int64" Name="connectivity" format="ascii">0 1</DataArray><DataArray type="Int64" Name="offsets" format="ascii">2</DataArray><DataArray type="UInt8" Name="types" format="ascii">3</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>"#;
    #[test]
    fn reads_topology_and_keeps_point_and_cell_fields_distinct() {
        let grid = parse_vtu(LINE).unwrap();
        assert_eq!(grid.points, [[0., 0., 0.], [1., 0., 0.]]);
        assert_eq!(grid.connectivity, [0, 1]);
        assert_eq!(grid.offsets, [2]);
        assert_eq!(grid.cell_types, [3]);
        assert_eq!(grid.fields[0].association, Association::Point);
        assert_eq!(grid.fields[1].association, Association::Cell);
        assert_eq!(grid.fields[1].values, [2., 3.]);
    }
    #[test]
    fn rejects_unbounded_untrusted_and_malformed_data() {
        for xml in [
            LINE.replace("0 1</DataArray>", "NaN 1</DataArray>"),
            LINE.replace("format=\"ascii\"", "format=\"binary\""),
            format!("<!DOCTYPE VTKFile [<!ENTITY x SYSTEM 'file:///etc/passwd'>]>{LINE}"),
        ] {
            assert!(parse_vtu(&xml).is_err());
        }
        assert!(parse_vtu(&" ".repeat(MAX_XML + 1)).is_err());
        assert!(
            parse_vtu(
                &LINE.replace(
                    "0 1</DataArray><DataArray type=\"Int64\" Name=\"offsets\"",
                    "0 9</DataArray><DataArray type=\"Int64\" Name=\"offsets\""
                )
            )
            .is_err()
        );
        assert!(
            parse_vtu(
                &LINE
                    .replace("type=\"Float64\" Name=\"potential\"", "type=\"UInt64\" Name=\"potential\"")
                    .replace(">0 1</DataArray></PointData>", ">9007199254740993 1</DataArray></PointData>")
            )
            .is_err()
        );
    }
    #[test]
    fn collection_preserves_coordinate_and_rejects_unsafe_paths() {
        let xml = r#"<VTKFile type="Collection"><Collection><DataSet timestep="2.5" part="1" file="step.vtu"/></Collection></VTKFile>"#;
        assert_eq!(parse_pvd(xml).unwrap(), vec![Frame { coordinate: 2.5, file: "step.vtu".into(), part: 1 }]);
        assert!(parse_pvd(&xml.replace("step.vtu", "../step.vtu")).is_err());
        assert!(parse_pvd(&xml.replace("2.5", "NaN")).is_err());
    }
    #[test]
    #[ignore = "requires SOLIDCRAFT_TEST_VTU pointing to a real companion output"]
    fn real_companion_poisson_result_preserves_topology_and_solution() {
        let path = std::env::var_os("SOLIDCRAFT_TEST_VTU").expect("set SOLIDCRAFT_TEST_VTU");
        let grid = parse_vtu(&std::fs::read_to_string(path).unwrap()).unwrap();
        assert_eq!(grid.points.len(), 65);
        assert_eq!(grid.offsets.len(), 64);
        let field = grid.fields.iter().find(|f| f.name == "potential").unwrap();
        assert_eq!(field.association, Association::Point);
        for (point, value) in grid.points.iter().zip(&field.values) {
            assert!((value - point[0] * (1.0 - point[0])).abs() < 1e-10);
        }
    }
}
