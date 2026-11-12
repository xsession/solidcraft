//! Optional native companion host. A solver executable never runs on the window thread.
//!
//! This adapter accepts only the bounded native Rust reference runner's versioned output.
//! A successful process exit is followed by independent evidence admission in the domain crate.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

#[cfg(not(target_arch = "wasm32"))]
mod process;
#[cfg(not(target_arch = "wasm32"))]
pub use process::{ProcessHost, configured_host};

/// Browser builds can prepare cases, but cannot start local executable companions.
#[cfg(target_arch = "wasm32")]
pub fn configured_host() -> Option<std::sync::Arc<dyn solidcraft_engine::multiphysics_jobs::Host>> {
    None
}
