//! Mass-action kinetics, ideal phases, binary precipitation and adaptive implicit batch reactors.
mod advanced;
mod aqueous;
mod kinetics;
mod reactor;
mod thermo;
pub use advanced::*;
pub use aqueous::*;
pub use kinetics::*;
pub use reactor::*;
pub use thermo::*;
pub const GAS_CONSTANT: f64 = 8.31446261815324;
