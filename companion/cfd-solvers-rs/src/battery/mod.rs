//! Lithium-ion cell models: spherical-particle diffusion, SPM and SPMe with degradation and
//! lumped thermal, the pseudo-2D Doyle-Fuller-Newman model, EIS, drive cycles, balanced series
//! packs and a 3-D voxel electrothermal pack, ported from upstream `cfd/battery/{lithium_ion,
//! dfn, pack, thermal3d}` (MIT, cfd_solvers contributors).
mod cell;
mod dfn;
mod pack;
pub use cell::*;
pub use dfn::*;
pub use pack::*;
