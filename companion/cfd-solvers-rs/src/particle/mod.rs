//! Particle physics: charged-particle pushers and Monte-Carlo collisions, 1-D/2-D/3-D
//! particle-in-cell (electrostatic, collocated and Yee-staggered electromagnetic), 3-D PIC
//! domain decomposition with in-memory/serialized message exchange, nuclear cross sections,
//! depletion and slab transport, plasma chemistry and breakdown estimates, wakefields and
//! Geant4-style condensed-history transport with dose scoring. Ported from upstream
//! `cfd/particle` (MIT, cfd_solvers contributors).
//!
//! The upstream MPI exchange entry point is not ported (no MPI offline); the serialized
//! exchange path it wraps is available as `run_serialized_distributed_pic_exchange_round_3d`.
pub mod nuclear;
pub mod pic1d;
pub mod pic2d;
pub mod pic3d;
pub mod pic3d_domain;
pub mod plasma;
pub mod pusher;
pub mod transport;
pub mod wakefield;

pub use nuclear::*;
pub use pic1d::*;
pub use pic2d::*;
pub use pic3d::*;
pub use pic3d_domain::*;
pub use plasma::*;
pub use pusher::*;
pub use transport::*;
pub use wakefield::*;
