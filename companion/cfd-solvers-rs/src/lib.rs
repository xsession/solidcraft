//! Native Rust solver implementations. See parity/status.json for exact coverage.
#![forbid(unsafe_code)]

pub mod acoustics;
pub mod atomistics;
pub mod battery;
pub mod case;
pub mod chemistry;
pub mod circuit;
pub mod core;
pub mod dg;
pub mod distributed;
pub mod electrochemistry;
pub mod em;
pub mod fdtd;
pub mod fem;
pub mod fvm;
pub mod io;
pub mod kinetic;
pub mod lbm;
pub mod mechanics;
pub mod multibody;
pub mod multiphysics;
pub mod optics;
pub mod particle;
pub mod physics;
pub mod porous;
pub mod power;
pub mod python;
pub mod rf;
pub mod tcad;
pub mod workflow;
pub use num_complex::Complex64;

#[derive(Debug, Clone, PartialEq)]
pub enum SolverError {
    InvalidInput(String),
    Singular(String),
    NumericalFailure(String),
    Unsupported(String),
    Io(String),
}
impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (kind, text) = match self {
            Self::InvalidInput(s) => ("invalid input", s),
            Self::Singular(s) => ("singular system", s),
            Self::NumericalFailure(s) => ("numerical failure", s),
            Self::Unsupported(s) => ("unsupported", s),
            Self::Io(s) => ("I/O", s),
        };
        write!(f, "{kind}: {text}")
    }
}
impl std::error::Error for SolverError {}
impl From<std::io::Error> for SolverError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e.to_string())
    }
}
impl From<serde_json::Error> for SolverError {
    fn from(e: serde_json::Error) -> Self {
        Self::InvalidInput(e.to_string())
    }
}
pub type Result<T> = std::result::Result<T, SolverError>;
pub(crate) fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(SolverError::InvalidInput(message.into()))
    }
}
pub(crate) fn positive(x: f64, name: &str) -> Result<()> {
    require(
        x.is_finite() && x > 0.0,
        format!("{name} must be positive and finite"),
    )
}
pub(crate) fn finite_slice(x: &[f64], name: &str) -> Result<()> {
    require(
        x.iter().all(|v| v.is_finite()),
        format!("{name} contains nonfinite values"),
    )
}
pub const UPSTREAM_REVISION: &str = "7cceaaaddafaf2515e01a86406e13d6037fdbf9a";
