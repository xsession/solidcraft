//! First-order finite elements, matching the upstream reference-element conventions.
mod adaptivity;
mod electromagnetics;
mod flow;
mod materials;
mod mesh;
mod reference;
mod scalar;
mod solid_em;
mod structural;
pub use adaptivity::*;
pub use electromagnetics::*;
pub use flow::*;
pub use materials::*;
pub use mesh::*;
pub use reference::*;
pub use scalar::*;
pub use solid_em::*;
pub use structural::*;
