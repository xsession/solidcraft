//! SolidCraft's Multiphysics workspace: studies and their case, the model studies, runs, the
//! solver monitor and results. A workspace of the app in `solidcraft-ui-egui`, kept in its own
//! crate so editing it rebuilds only this crate. The app calls it through the hooks [`register`]
//! installs.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

mod multiphysics_case_view;
mod multiphysics_field_view;
mod multiphysics_material_library;
mod multiphysics_property_grid;
pub mod multiphysics_workspace;

// The app's modules the workspace draws with, under the names it has always used.
pub(crate) use solidcraft_ui_egui::{
    SolveApp, context_menu, icons, material_editor, navigation, prefs, render_workspace, ribbon, scenario, ui_kit, viewport, viewport_layout,
    workspace,
};
#[cfg(test)]
pub(crate) use solidcraft_ui_egui::{documents, fonts};
#[cfg(not(target_arch = "wasm32"))]
pub(crate) use solidcraft_ui_egui::open_meshing;

/// Installs the workspace in the app. Call once at start-up, before the first frame.
pub fn register() {
    solidcraft_ui_egui::multiphysics_workspace::register(multiphysics_workspace::hooks());
}
