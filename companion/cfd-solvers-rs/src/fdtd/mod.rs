//! Yee time-domain electromagnetics: 1-D with ADE dispersion and absorbing layers, 3-D Cartesian,
//! axisymmetric TM, monitors and ports, near-to-far-field, SAR and EMC helpers.
mod cylindrical;
mod maxwell1d;
mod maxwell3d;
mod monitor;
mod postprocess;
pub use cylindrical::*;
pub use maxwell1d::*;
pub use maxwell3d::*;
pub use monitor::*;
pub use postprocess::*;
pub const EPSILON0: f64 = 8.8541878128e-12;
pub const MU0: f64 = 1.25663706212e-6;
