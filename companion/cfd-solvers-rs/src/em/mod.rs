//! Frequency-domain electromagnetics: waveguide ports, cables, metasurfaces, 1-D driven and
//! eigen solves, signal/power integrity, Nedelec edge elements in 2-D and 3-D with resonator
//! quality, a-posteriori estimation and conforming refinement, and photonic modes.
mod circuits;
mod edge;
mod photonic;
pub use circuits::*;
pub use edge::*;
pub use photonic::*;
