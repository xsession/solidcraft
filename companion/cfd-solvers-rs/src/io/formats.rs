//! Binary field checkpoints, the FVM HDF5 checkpoint seam, OpenFOAM polyMesh and
//! volScalar/volVector field I/O, and the XDMF mesh contract, ported from upstream
//! `cfd/io/{checkpoint, fvm_checkpoint, openfoam_io, xdmf_contract}` (MIT, cfd_solvers
//! contributors).
#![allow(clippy::needless_range_loop)]
use crate::{
    core::{parse_float_prefix, Vec3},
    fem::Mesh2D,
    fvm::{BoundaryPatch, Face, PolyCell, PolyMesh},
    lbm::format_general,
    require, Result, SolverError,
};
use std::path::Path;

const MAGIC: [u8; 8] = *b"CFDCHK1\0";
const MAX_ITEMS: u64 = 100_000_000;
const MAX_STRING: u64 = 1 << 20;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CheckpointField {
    pub name: String,
    pub values: Vec<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Checkpoint {
    pub case_name: String,
    pub time: f64,
    pub step: usize,
    pub fields: Vec<CheckpointField>,
}

fn validate_name(name: &str) -> Result<()> {
    require(
        !name.is_empty() && name.len() as u64 <= MAX_STRING,
        "invalid checkpoint name",
    )?;
    require(
        !name.contains(['\0', '/', '\\']),
        "checkpoint names may not contain path separators",
    )
}

fn put_u64(o: &mut Vec<u8>, v: u64) {
    o.extend_from_slice(&v.to_le_bytes());
}

fn put_string(o: &mut Vec<u8>, s: &str) {
    put_u64(o, s.len() as u64);
    o.extend_from_slice(s.as_bytes());
}

/// Little-endian `CFDCHK1` byte image (version 1).
pub fn encode_checkpoint(
    case_name: &str,
    time: f64,
    step: usize,
    fields: &[(&str, &[f64])],
) -> Result<Vec<u8>> {
    require(
        !case_name.is_empty()
            && time.is_finite()
            && time >= 0.0
            && fields.len() as u64 <= MAX_ITEMS,
        "invalid checkpoint metadata",
    )?;
    validate_name(case_name)?;
    let mut names: Vec<&str> = Vec::with_capacity(fields.len());
    for (name, values) in fields {
        validate_name(name)?;
        require(
            values.len() as u64 <= MAX_ITEMS && !names.contains(name),
            "invalid or duplicate checkpoint field",
        )?;
        require(
            values.iter().all(|v| v.is_finite()),
            "checkpoint field contains a non-finite value",
        )?;
        names.push(name);
    }
    let mut o = MAGIC.to_vec();
    put_u64(&mut o, 1);
    put_string(&mut o, case_name);
    put_u64(&mut o, time.to_bits());
    put_u64(&mut o, step as u64);
    put_u64(&mut o, fields.len() as u64);
    for (name, values) in fields {
        put_string(&mut o, name);
        put_u64(&mut o, values.len() as u64);
        for v in *values {
            put_u64(&mut o, v.to_bits());
        }
    }
    Ok(o)
}

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8]> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|e| *e <= self.bytes.len())
            .ok_or_else(|| SolverError::Io("cannot read checkpoint".into()))?;
        let s = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(s)
    }
    fn u64(&mut self) -> Result<u64> {
        let mut b = [0u8; 8];
        b.copy_from_slice(self.take(8)?);
        Ok(u64::from_le_bytes(b))
    }
    fn string(&mut self) -> Result<String> {
        let n = self.u64()?;
        require(n <= MAX_STRING, "checkpoint string is too long")?;
        String::from_utf8(self.take(n as usize)?.to_vec())
            .map_err(|_| SolverError::Io("checkpoint string is not UTF-8".into()))
    }
}

pub fn decode_checkpoint(bytes: &[u8]) -> Result<Checkpoint> {
    let mut r = Reader { bytes, pos: 0 };
    require(
        r.take(8)? == MAGIC && r.u64()? == 1,
        "unsupported checkpoint format",
    )?;
    let case_name = r.string()?;
    validate_name(&case_name)?;
    let time = f64::from_bits(r.u64()?);
    let step = usize::try_from(r.u64()?)
        .map_err(|_| SolverError::Io("checkpoint step overflows".into()))?;
    require(time.is_finite() && time >= 0.0, "invalid checkpoint time")?;
    let count = r.u64()?;
    require(count <= MAX_ITEMS, "checkpoint has too many fields")?;
    let mut fields: Vec<CheckpointField> = Vec::new();
    for _ in 0..count {
        let name = r.string()?;
        validate_name(&name)?;
        let n = r.u64()?;
        require(
            n <= MAX_ITEMS && n.saturating_mul(8) <= (bytes.len() - r.pos) as u64,
            "checkpoint field is too large",
        )?;
        let mut values = Vec::with_capacity(n as usize);
        for _ in 0..n {
            let v = f64::from_bits(r.u64()?);
            require(
                v.is_finite(),
                "checkpoint field contains a non-finite value",
            )?;
            values.push(v);
        }
        require(
            !fields.iter().any(|f| f.name == name),
            "checkpoint contains duplicate fields",
        )?;
        fields.push(CheckpointField { name, values });
    }
    Ok(Checkpoint {
        case_name,
        time,
        step,
        fields,
    })
}

pub fn write_checkpoint(
    path: impl AsRef<Path>,
    case_name: &str,
    time: f64,
    step: usize,
    fields: &[(&str, &[f64])],
) -> Result<()> {
    let bytes = encode_checkpoint(case_name, time, step, fields)?;
    std::fs::write(path, bytes)
        .map_err(|e| SolverError::Io(format!("cannot open checkpoint for writing: {e}")))
}

pub fn read_checkpoint(path: impl AsRef<Path>) -> Result<Checkpoint> {
    let bytes = std::fs::read(path)
        .map_err(|e| SolverError::Io(format!("cannot open checkpoint for reading: {e}")))?;
    decode_checkpoint(&bytes)
}

#[derive(Clone, Debug)]
pub struct FvmCheckpoint {
    pub mesh: PolyMesh,
    pub time: f64,
    pub step: usize,
    pub scalar_fields: Vec<CheckpointField>,
}

/// No HDF5 library is linked, matching upstream builds without `CFD_ENABLE_HDF5`.
pub fn fvm_hdf5_checkpoint_available() -> bool {
    false
}

fn valid_identifier(name: &str, max: usize) -> bool {
    !name.is_empty()
        && name.len() <= max
        && name
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-'))
}

/// Validates like upstream, then reports `Unsupported` (no HDF5 here); use [`write_checkpoint`]
/// with a serialized mesh instead.
pub fn write_fvm_hdf5_checkpoint(
    _path: impl AsRef<Path>,
    mesh: &PolyMesh,
    scalar_fields: &[(&str, &[f64])],
    time: f64,
    _step: usize,
) -> Result<()> {
    require(
        time.is_finite() && time >= 0.0,
        "invalid FVM checkpoint time",
    )?;
    let mut names: Vec<&str> = Vec::new();
    for (name, values) in scalar_fields {
        require(
            valid_identifier(name, usize::MAX),
            "FVM HDF5 field names are restricted to [A-Za-z0-9_.-]+",
        )?;
        require(
            values.len() == mesh.cells().len(),
            "FVM checkpoint scalar field size mismatch",
        )?;
        require(
            values.iter().all(|v| v.is_finite()),
            "FVM checkpoint scalar field contains a non-finite value",
        )?;
        require(!names.contains(name), "duplicate FVM checkpoint field name")?;
        names.push(name);
    }
    Err(SolverError::Unsupported(
        "FVM HDF5 checkpoint I/O requires CFD_ENABLE_HDF5 and an HDF5 library".into(),
    ))
}

pub fn read_fvm_hdf5_checkpoint(_path: impl AsRef<Path>) -> Result<FvmCheckpoint> {
    Err(SolverError::Unsupported(
        "FVM HDF5 checkpoint I/O requires CFD_ENABLE_HDF5 and an HDF5 library".into(),
    ))
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct XdmfMeshContract {
    pub node_count: usize,
    pub triangle_count: usize,
    pub field_count: usize,
    pub field_name: String,
}

pub fn valid_xdmf_identifier(name: &str) -> bool {
    valid_identifier(name, 128)
}

pub fn make_xdmf_mesh_contract(
    mesh: &Mesh2D,
    scalar: &[f64],
    field_name: &str,
) -> Result<XdmfMeshContract> {
    mesh.validate()?;
    require(
        scalar.is_empty() || scalar.len() == mesh.node_count(),
        "XDMF scalar size mismatch",
    )?;
    require(
        scalar.is_empty() || valid_xdmf_identifier(field_name),
        "XDMF field name contains unsupported characters",
    )?;
    require(
        scalar.iter().all(|v| v.is_finite()),
        "XDMF scalar contains a non-finite value",
    )?;
    Ok(XdmfMeshContract {
        node_count: mesh.node_count(),
        triangle_count: mesh.element_count(),
        field_count: scalar.len(),
        field_name: field_name.into(),
    })
}

/// Removes `//` line comments, then `/* */` block comments.
fn strip_comments(s: &str) -> String {
    let mut no_line = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(p) = rest.find("//") {
        no_line.push_str(&rest[..p]);
        rest = &rest[p..];
        let end = rest.find(['\n', '\r']).unwrap_or(rest.len());
        rest = &rest[end..];
    }
    no_line.push_str(rest);
    let mut out = String::with_capacity(no_line.len());
    let mut rest = no_line.as_str();
    while let Some(p) = rest.find("/*") {
        match rest[p + 2..].find("*/") {
            Some(e) => {
                out.push_str(&rest[..p]);
                rest = &rest[p + 2 + e + 2..];
            }
            None => break,
        }
    }
    out.push_str(rest);
    out
}

fn skip_ws(b: &[u8], mut p: usize) -> usize {
    while p < b.len() && b[p].is_ascii_whitespace() {
        p += 1;
    }
    p
}

/// Text inside the first `N (` ... `)` list, matching nested parentheses.
fn list_body(text: &str) -> Result<&str> {
    let b = text.as_bytes();
    let mut open = None;
    let mut p = 0;
    while p < b.len() {
        if b[p].is_ascii_digit() {
            let mut q = p;
            while q < b.len() && b[q].is_ascii_digit() {
                q += 1;
            }
            let r = skip_ws(b, q);
            if r < b.len() && b[r] == b'(' {
                open = Some(r);
                break;
            }
            p = q;
            continue;
        }
        p += 1;
    }
    let open = open.ok_or_else(|| SolverError::Io("OpenFOAM list header not found".into()))?;
    let mut depth = 0;
    for (i, c) in b.iter().enumerate().skip(open) {
        if *c == b'(' {
            depth += 1;
        } else if *c == b')' {
            depth -= 1;
            if depth == 0 {
                return Ok(&text[open + 1..i]);
            }
        }
    }
    Err(SolverError::Io("unterminated OpenFOAM list".into()))
}

fn stod(s: &str) -> Result<f64> {
    parse_float_prefix(s)
        .map(|v| v.0)
        .ok_or_else(|| SolverError::Io(format!("invalid OpenFOAM number: {s}")))
}

fn num_run(b: &[u8], p: usize) -> usize {
    let mut q = p;
    while q < b.len() && (b[q].is_ascii_digit() || matches!(b[q], b'-' | b'+' | b'.' | b'e' | b'E'))
    {
        q += 1;
    }
    q
}

/// `\(\s*N\s+N\s+N\s*\)` tuples.
fn parse_points(body: &str) -> Result<Vec<Vec3>> {
    let b = body.as_bytes();
    let mut out = Vec::new();
    let mut p = 0;
    'scan: while p < b.len() {
        if b[p] != b'(' {
            p += 1;
            continue;
        }
        let mut q = skip_ws(b, p + 1);
        let mut v = [0.0; 3];
        for (k, slot) in v.iter_mut().enumerate() {
            let e = num_run(b, q);
            if e == q || (k < 2 && (e >= b.len() || !b[e].is_ascii_whitespace())) {
                p += 1;
                continue 'scan;
            }
            *slot = stod(&body[q..e])?;
            q = if k < 2 { skip_ws(b, e) } else { e };
        }
        q = skip_ws(b, q);
        if q < b.len() && b[q] == b')' {
            out.push(Vec3::new(v[0], v[1], v[2]));
            p = q + 1;
        } else {
            p += 1;
        }
    }
    require(!out.is_empty(), "OpenFOAM points list empty")?;
    Ok(out)
}

fn leading_unsigned(text: &str) -> Vec<usize> {
    let mut out = Vec::new();
    for t in text.split_ascii_whitespace() {
        match t.parse::<usize>() {
            Ok(v) => out.push(v),
            Err(_) => {
                let digits: String = t.chars().take_while(char::is_ascii_digit).collect();
                if let Ok(v) = digits.parse() {
                    out.push(v);
                }
                break;
            }
        }
    }
    out
}

/// `N(i j k ...)` faces.
fn parse_faces(body: &str) -> Result<Vec<Vec<usize>>> {
    let b = body.as_bytes();
    let mut out = Vec::new();
    let mut p = 0;
    while p < b.len() {
        if !b[p].is_ascii_digit() {
            p += 1;
            continue;
        }
        let mut q = p;
        while q < b.len() && b[q].is_ascii_digit() {
            q += 1;
        }
        let r = skip_ws(b, q);
        if r >= b.len() || b[r] != b'(' {
            p = q;
            continue;
        }
        let Some(close) = body[r + 1..].find(')') else {
            p = q;
            continue;
        };
        let n: usize = body[p..q]
            .parse()
            .map_err(|_| SolverError::Io("OpenFOAM face size overflows".into()))?;
        let f = leading_unsigned(&body[r + 1..r + 1 + close]);
        require(f.len() == n && n >= 3, "malformed OpenFOAM face")?;
        out.push(f);
        p = r + 1 + close + 1;
    }
    require(!out.is_empty(), "OpenFOAM faces list empty")?;
    Ok(out)
}

fn keyword_number(block: &str, key: &str) -> Option<usize> {
    let b = block.as_bytes();
    let mut from = 0;
    while let Some(off) = block[from..].find(key) {
        let start = from + off + key.len();
        let q = skip_ws(b, start);
        if q > start {
            let mut e = q;
            while e < b.len() && b[e].is_ascii_digit() {
                e += 1;
            }
            let s = skip_ws(b, e);
            if e > q && s < b.len() && b[s] == b';' {
                if let Ok(v) = block[q..e].parse() {
                    return Some(v);
                }
            }
        }
        from = from + off + 1;
    }
    None
}

fn parse_boundary(body: &str) -> Result<Vec<(String, usize, usize)>> {
    let b = body.as_bytes();
    let is_name = |c: u8| c.is_ascii_alphanumeric() || matches!(c, b'_' | b'.' | b'-');
    let mut out = Vec::new();
    let mut p = 0;
    while p < b.len() {
        if !is_name(b[p]) {
            p += 1;
            continue;
        }
        let mut q = p;
        while q < b.len() && is_name(b[q]) {
            q += 1;
        }
        let r = skip_ws(b, q);
        let close = if r < b.len() && b[r] == b'{' {
            body[r + 1..].find('}')
        } else {
            None
        };
        let Some(close) = close else {
            p += 1;
            continue;
        };
        let block = &body[r + 1..r + 1 + close];
        let (Some(n), Some(start)) = (
            keyword_number(block, "nFaces"),
            keyword_number(block, "startFace"),
        ) else {
            return Err(SolverError::Io(
                "OpenFOAM boundary patch missing nFaces/startFace".into(),
            ));
        };
        out.push((body[p..q].to_string(), n, start));
        p = r + 1 + close + 1;
    }
    Ok(out)
}

fn read_text(path: &Path) -> Result<String> {
    std::fs::read_to_string(path).map_err(|e| {
        SolverError::Io(format!(
            "cannot open OpenFOAM file: {} ({e})",
            path.display()
        ))
    })
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// Builds a [`PolyMesh`] from the five polyMesh text files (`points`, `faces`, `owner`,
/// `neighbour`, `boundary`), already loaded.
pub fn parse_openfoam_polymesh(
    points: &str,
    faces: &str,
    owner: &str,
    neighbour: &str,
    boundary: &str,
) -> Result<PolyMesh> {
    let points = parse_points(list_body(&strip_comments(points))?)?;
    let raw_faces = parse_faces(list_body(&strip_comments(faces))?)?;
    let owner = leading_unsigned(list_body(&strip_comments(owner))?);
    let neighbour = leading_unsigned(list_body(&strip_comments(neighbour))?);
    require(
        owner.len() == raw_faces.len(),
        "OpenFOAM owner count mismatch",
    )?;
    let cells = owner
        .iter()
        .chain(&neighbour)
        .map(|c| c + 1)
        .max()
        .unwrap_or(0);
    require(cells <= 1 << 26, "OpenFOAM mesh too large")?;
    let boundary_text = strip_comments(boundary);
    let raw_patches = parse_boundary(list_body(&boundary_text)?)?;
    let mut face_patch = vec![None; raw_faces.len()];
    for (p, (_, n, start)) in raw_patches.iter().enumerate() {
        for f in *start..start.saturating_add(*n) {
            require(f < raw_faces.len(), "OpenFOAM boundary face range invalid")?;
            face_patch[f] = Some(p);
        }
    }
    let mut out_faces = Vec::with_capacity(raw_faces.len());
    let mut center_sum = vec![Vec3::new(0.0, 0.0, 0.0); cells];
    let mut center_weight = vec![0.0; cells];
    for (f, rf) in raw_faces.iter().enumerate() {
        require(
            rf.iter().all(|i| *i < points.len()),
            "OpenFOAM face point index out of range",
        )?;
        let mut c = Vec3::new(0.0, 0.0, 0.0);
        for i in rf {
            c = Vec3::new(c.x + points[*i].x, c.y + points[*i].y, c.z + points[*i].z);
        }
        let n = rf.len() as f64;
        let fc = Vec3::new(c.x / n, c.y / n, c.z / n);
        let o = points[rf[0]];
        let mut sf = Vec3::new(0.0, 0.0, 0.0);
        for i in 1..rf.len() - 1 {
            let (a, b) = (points[rf[i]], points[rf[i + 1]]);
            let (a, b) = (
                Vec3::new(a.x - o.x, a.y - o.y, a.z - o.z),
                Vec3::new(b.x - o.x, b.y - o.y, b.z - o.z),
            );
            let cr = Vec3::new(
                a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x,
            );
            sf = Vec3::new(sf.x + cr.x * 0.5, sf.y + cr.y * 0.5, sf.z + cr.z * 0.5);
        }
        let nb = neighbour.get(f).copied();
        let w = dot(sf, sf).sqrt().max(1e-30);
        for cell in std::iter::once(owner[f]).chain(nb) {
            let s = &mut center_sum[cell];
            *s = Vec3::new(s.x + fc.x * w, s.y + fc.y * w, s.z + fc.z * w);
            center_weight[cell] += w;
        }
        out_faces.push(Face {
            owner: owner[f],
            neighbour: nb,
            center: fc,
            area: sf,
            patch: if nb.is_some() { None } else { face_patch[f] },
        });
    }
    let mut vol = vec![0.0; cells];
    for f in &out_faces {
        vol[f.owner] += dot(f.area, f.center) / 3.0;
        if let Some(n) = f.neighbour {
            vol[n] -= dot(f.area, f.center) / 3.0;
        }
    }
    let mut out_cells = Vec::with_capacity(cells);
    for c in 0..cells {
        require(center_weight[c] > 0.0, "OpenFOAM cell has no faces")?;
        let s = center_sum[c];
        let w = center_weight[c];
        let volume = vol[c].abs();
        require(volume > 0.0, "OpenFOAM cell volume non-positive")?;
        out_cells.push(PolyCell {
            center: Vec3::new(s.x / w, s.y / w, s.z / w),
            volume,
        });
    }
    let patches = raw_patches
        .into_iter()
        .map(|(name, _, _)| BoundaryPatch { name })
        .collect();
    PolyMesh::new(out_cells, out_faces, patches)
}

pub fn read_openfoam_polymesh(poly_mesh_directory: impl AsRef<Path>) -> Result<PolyMesh> {
    let d = poly_mesh_directory.as_ref();
    parse_openfoam_polymesh(
        &read_text(&d.join("points"))?,
        &read_text(&d.join("faces"))?,
        &read_text(&d.join("owner"))?,
        &read_text(&d.join("neighbour"))?,
        &read_text(&d.join("boundary"))?,
    )
}

fn leading_doubles(text: &str) -> Vec<f64> {
    let mut out = Vec::new();
    for t in text.split_ascii_whitespace() {
        match parse_float_prefix(t) {
            Some((v, used)) => {
                out.push(v);
                if used != t.len() {
                    break;
                }
            }
            None => break,
        }
    }
    out
}

/// `internalField uniform v;` or `internalField nonuniform List<scalar> N ( ... );`.
pub fn parse_openfoam_vol_scalar_field(text: &str, expected_cells: usize) -> Result<Vec<f64>> {
    let text = strip_comments(text);
    let b = text.as_bytes();
    let mut from = 0;
    while let Some(off) = text[from..].find("internalField") {
        let start = from + off + "internalField".len();
        let q = skip_ws(b, start);
        if q > start && text[q..].starts_with("uniform") {
            let r = skip_ws(b, q + "uniform".len());
            let e = num_run(b, r);
            if r > q + "uniform".len()
                && e > r
                && skip_ws(b, e) < b.len()
                && b[skip_ws(b, e)] == b';'
            {
                require(expected_cells <= 1 << 28, "OpenFOAM field too large")?;
                return Ok(vec![stod(&text[r..e])?; expected_cells]);
            }
        }
        from = from + off + 1;
    }
    let key = "nonuniform";
    let mut from = 0;
    while let Some(off) = text[from..].find("internalField") {
        let start = from + off + "internalField".len();
        from = from + off + 1;
        let q = skip_ws(b, start);
        if q == start || !text[q..].starts_with(key) {
            continue;
        }
        let r = skip_ws(b, q + key.len());
        if r == q + key.len() || !text[r..].starts_with("List<scalar>") {
            continue;
        }
        let s = skip_ws(b, r + "List<scalar>".len());
        let mut e = s;
        while e < b.len() && b[e].is_ascii_digit() {
            e += 1;
        }
        let o = skip_ws(b, e);
        if s == r + "List<scalar>".len() || e == s || o >= b.len() || b[o] != b'(' {
            continue;
        }
        // Lazy body: the first ')' followed by optional whitespace and ';'.
        let mut k = o + 1;
        let close = loop {
            match text[k..].find(')') {
                Some(c) => {
                    let after = skip_ws(b, k + c + 1);
                    if after < b.len() && b[after] == b';' {
                        break Some(k + c);
                    }
                    k += c + 1;
                }
                None => break None,
            }
        };
        let Some(close) = close else { continue };
        let n: usize = text[s..e]
            .parse()
            .map_err(|_| SolverError::Io("OpenFOAM scalar count overflows".into()))?;
        require(n == expected_cells, "OpenFOAM scalar field size mismatch")?;
        let v = leading_doubles(&text[o + 1..close]);
        require(v.len() == n, "malformed OpenFOAM scalar field")?;
        return Ok(v);
    }
    Err(SolverError::Io(
        "unsupported OpenFOAM scalar internalField".into(),
    ))
}

pub fn read_openfoam_vol_scalar_field(
    field_file: impl AsRef<Path>,
    expected_cells: usize,
) -> Result<Vec<f64>> {
    parse_openfoam_vol_scalar_field(&read_text(field_file.as_ref())?, expected_cells)
}

/// ASCII volScalarField text with six significant digits (upstream defaults: dimensionless,
/// `boundaryField {}`).
pub fn openfoam_vol_scalar_field_text(
    object_name: &str,
    values: &[f64],
    dimensions: &str,
    boundary_field: &str,
) -> String {
    let mut o = format!("FoamFile\n{{ version 2.0; format ascii; class volScalarField; object {object_name}; }}\ndimensions {dimensions};\ninternalField nonuniform List<scalar>\n{}\n(\n", values.len());
    for v in values {
        o += &format_general(*v, 6);
        o.push('\n');
    }
    o + ")\n;\n" + boundary_field
}

/// ASCII volVectorField text (upstream defaults: `[0 1 -1 0 0 0 0]`, `boundaryField {}`).
pub fn openfoam_vol_vector_field_text(
    object_name: &str,
    values: &[Vec3],
    dimensions: &str,
    boundary_field: &str,
) -> String {
    let mut o = format!("FoamFile\n{{ version 2.0; format ascii; class volVectorField; object {object_name}; }}\ndimensions {dimensions};\ninternalField nonuniform List<vector>\n{}\n(\n", values.len());
    for v in values {
        o += &format!(
            "({} {} {})\n",
            format_general(v.x, 6),
            format_general(v.y, 6),
            format_general(v.z, 6)
        );
    }
    o + ")\n;\n" + boundary_field
}

pub fn write_openfoam_vol_scalar_field(
    field_file: impl AsRef<Path>,
    object_name: &str,
    values: &[f64],
    dimensions: &str,
    boundary_field: &str,
) -> Result<()> {
    std::fs::write(
        field_file,
        openfoam_vol_scalar_field_text(object_name, values, dimensions, boundary_field),
    )
    .map_err(|e| SolverError::Io(format!("cannot create OpenFOAM scalar field: {e}")))
}

pub fn write_openfoam_vol_vector_field(
    field_file: impl AsRef<Path>,
    object_name: &str,
    values: &[Vec3],
    dimensions: &str,
    boundary_field: &str,
) -> Result<()> {
    std::fs::write(
        field_file,
        openfoam_vol_vector_field_text(object_name, values, dimensions, boundary_field),
    )
    .map_err(|e| SolverError::Io(format!("cannot create OpenFOAM vector field: {e}")))
}
