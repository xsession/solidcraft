//! Finite-volume reference solvers. General polyhedral and multiphase solvers are tracked separately.
mod adaptive_mesh;
mod collocated;
mod combustion;
mod compressible;
mod diffusion;
mod euler_euler;
mod film;
mod models;
mod operators;
mod particles;
mod periodic2d;
mod poly_mesh;
mod postprocess;
mod pressure_velocity;
mod rans;
mod schemes;
mod spray;
mod thermal;
mod transport;
mod utilities;
mod vof;
pub use adaptive_mesh::*;
pub use collocated::*;
pub use combustion::*;
pub use compressible::*;
pub use diffusion::*;
pub use euler_euler::*;
pub use film::*;
pub use models::*;
pub use operators::*;
pub use particles::*;
pub use periodic2d::*;
pub use poly_mesh::*;
pub use postprocess::*;
pub use pressure_velocity::*;
pub use rans::*;
pub use schemes::*;
pub use spray::*;
pub use thermal::*;
pub use transport::*;
pub use utilities::*;
pub use vof::*;
pub fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}
