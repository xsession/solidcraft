//! Electrochemical kinetics, corrosion, mixed potentials, implicit electrode cells and
//! Nernst-Planck transport, ported from upstream `cfd/electrochemistry` and
//! `cfd/solvers/electrochemistry` (MIT, cfd_solvers contributors).
mod basics;
mod corrosion;
mod transport;
pub use basics::*;
pub use corrosion::*;
pub use transport::*;
