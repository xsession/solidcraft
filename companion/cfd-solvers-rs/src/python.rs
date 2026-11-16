//! Safe Rust counterpart of upstream `cfd/python/c_api.h` (MIT, cfd_solvers contributors).
//!
//! Upstream exports two `extern "C"` symbols for the Python bindings: a version string and an
//! inline-XDMF writer for triangle meshes taking raw pointers, element counts and a caller
//! error buffer. Exporting a C ABI needs `unsafe` pointer handling and `#[no_mangle]`, which
//! this crate does not use, so the functions here take slices and return `Result`; a thin FFI
//! shim can wrap them where a C ABI is required.
use crate::{
    fem::{Mesh2D, Node2, Tri3},
    io::write_xdmf_inline,
    require, Result,
};
use std::path::Path;

/// Upstream library version reported by `cfd_solvers_version()` at the pinned revision.
pub const UPSTREAM_VERSION: &str = "0.19.2";

/// Version of this Rust port (the crate version).
pub fn cfd_solvers_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Mirrors `cfd_write_xdmf_tri3`: `xy` holds `node_count` interleaved coordinates, `tri`
/// holds `triangle_count` node triples, `scalar` (when given) holds `node_count` values and
/// `scalar_name` defaults to `"field"`. Slices shorter than the counts are an error rather
/// than an out-of-bounds read.
pub fn write_xdmf_tri3(
    path: &Path,
    xy: &[f64],
    node_count: usize,
    tri: &[usize],
    triangle_count: usize,
    scalar: Option<&[f64]>,
    scalar_name: Option<&str>,
) -> Result<()> {
    require(
        !path.as_os_str().is_empty() && node_count != 0 && triangle_count != 0,
        "path, coordinates and triangles are required",
    )?;
    let coordinates = xy.get(..2 * node_count).ok_or_else(|| {
        crate::SolverError::InvalidInput("coordinate buffer shorter than node count".into())
    })?;
    let connectivity = tri.get(..3 * triangle_count).ok_or_else(|| {
        crate::SolverError::InvalidInput("triangle buffer shorter than triangle count".into())
    })?;
    let values = match scalar {
        Some(s) => s
            .get(..node_count)
            .ok_or_else(|| {
                crate::SolverError::InvalidInput("scalar buffer shorter than node count".into())
            })?
            .to_vec(),
        None => Vec::new(),
    };
    let mesh = Mesh2D {
        nodes: coordinates
            .chunks_exact(2)
            .map(|c| Node2 { x: c[0], y: c[1] })
            .collect(),
        triangles: connectivity
            .chunks_exact(3)
            .map(|c| Tri3 {
                node: [c[0], c[1], c[2]],
            })
            .collect(),
        boundary_edges: Vec::new(),
        boundary_node: vec![false; node_count],
    };
    write_xdmf_inline(path, &mesh, &values, scalar_name.unwrap_or("field"))
}
