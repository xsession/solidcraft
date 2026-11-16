//! Solid mechanics references: constitutive updates, finite-strain elements, continuum
//! baselines, meshfree graphs, phase-field fracture and Cosserat rods.
mod constitutive;
mod continuum;
mod cosserat;
mod meshfree;
mod phase_field;
pub use constitutive::*;
pub use continuum::*;
pub use cosserat::*;
pub use meshfree::*;
pub use phase_field::*;
