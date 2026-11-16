//! Multiphysics coupling framework, ported from upstream `cfd/multiphysics`
//! (MIT, cfd_solvers contributors): field registry with SI unit signatures,
//! conservative transfers, mesh projections, partitioned fixed-point coupling with
//! Aitken relaxation, block-coupled linear systems, staged pipelines, interface laws,
//! Joule heating, Pennes bioheat/SAR, porous electrothermal and EM/circuit port coupling.
mod bioheat;
mod block_system;
mod coupling;
mod electro_thermal;
mod fields;
mod interface;
mod pipeline;
mod projection;
mod transfer;
pub use bioheat::*;
pub use block_system::*;
pub use coupling::*;
pub use electro_thermal::*;
pub use fields::*;
pub use interface::*;
pub use pipeline::*;
pub use projection::*;
pub use transfer::*;
