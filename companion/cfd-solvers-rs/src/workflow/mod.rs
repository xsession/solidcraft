//! Workflow tooling ported from upstream `cfd/workflow` (MIT, cfd_solvers contributors):
//! parameter sweeps, coordinate-search inversion and restartable checkpointed workflows
//! (`runner`); unit-tagged case configs, materials, provenance and the validation catalogue
//! (`schema`); design-of-experiments campaigns with registries, command plans, local launching
//! and log parsing (`campaign`); file-based control directives and scheduler parsing
//! (`campaign_control`); and multi-server / Docker deployment script generation (`deploy`).
mod campaign;
mod campaign_control;
mod deploy;
mod runner;
mod schema;
pub(crate) mod text;
pub use campaign::*;
pub use campaign_control::*;
pub use deploy::*;
pub use runner::*;
pub use schema::*;
