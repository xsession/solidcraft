//! The Multiphysics workspace's place in the app: its tab, when it may be entered, and the hooks
//! its crate (`solidcraft-multiphysics-ui`) registers at start-up. The workspace itself lives in
//! that crate, so editing it rebuilds only that crate and not the whole interface. Without the
//! hooks (a host that does not register them) the workspace cannot be entered.

use std::any::Any;
use std::sync::OnceLock;

use crate::SolveApp;

/// The tab the workspace occupies in `app.ui.tab`.
pub const TAB: &str = "MULTIPHYSICS";

/// The workspace's own state, owned by its crate. The app only keeps it, moves it between
/// documents and hands it back.
#[derive(Default)]
pub struct Slot(pub Option<Box<dyn Any>>);

impl std::fmt::Debug for Slot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(if self.0.is_some() { "Slot(..)" } else { "Slot(empty)" })
    }
}

/// What the workspace's crate provides.
pub struct Hooks {
    /// Draws the workspace in the window's central area.
    pub show: fn(&mut SolveApp, &mut egui::Ui),
    /// Polls its jobs and workers once a frame, whether or not it is shown.
    pub tick: fn(&mut SolveApp, &egui::Context),
    /// Forgets the document's workspace state on New and Open.
    pub reset_document: fn(&mut SolveApp),
    /// Finishes its work before the host tears down its resources.
    pub close: fn(&mut SolveApp),
    /// Takes the shown document's view state (to keep it while another document is shown).
    pub take_view: fn(&mut SolveApp) -> Slot,
    /// Puts back a document's view state taken by `take_view`.
    pub put_view: fn(&mut SolveApp, Slot),
}

static HOOKS: OnceLock<Hooks> = OnceLock::new();

/// Installs the workspace. The first registration wins; later ones are ignored.
pub fn register(hooks: Hooks) {
    let _ = HOOKS.set(hooks);
}

fn hooks() -> Option<&'static Hooks> {
    HOOKS.get()
}

pub fn active(app: &SolveApp) -> bool {
    app.ui.tab == TAB
}

pub fn workspace_name(app: &SolveApp) -> &'static str {
    if active(app) {
        "multiphysics"
    } else if crate::render_workspace::active(app) {
        "render"
    } else {
        "design"
    }
}

pub fn can_enter(app: &SolveApp) -> bool {
    hooks().is_some()
        && app.session.active_sketch.is_none()
        && app.tool.is_none()
        && app.dialog.is_none()
        && !app.file_load.busy()
        && !app.file_picker_busy()
}

/// UI workspace routing is separate from numerical document commands.
pub fn switch(app: &mut SolveApp, workspace: &str) -> Result<(), String> {
    match workspace {
        "multiphysics" if hooks().is_some() => {
            if !can_enter(app) {
                return Err(solidcraft_ui_text::tr("multiphysics.finish_operation").into());
            }
            if crate::render_workspace::active(app) {
                crate::render_workspace::enter(app, false)?;
            }
            crate::timeline::stop_playback(app);
            app.ui.tab = TAB.into();
            app.home.open = false;
            app.viewport.hover = None;
            app.cam_anim = None;
            Ok(())
        }
        "design" | "render" => crate::render_workspace::enter(app, workspace == "render"),
        _ => Err(solidcraft_ui_text::tr("multiphysics.unknown_workspace").into()),
    }
}

pub fn show(app: &mut SolveApp, ui: &mut egui::Ui) {
    if let Some(hooks) = hooks() {
        (hooks.show)(app, ui);
    }
}

pub fn tick(app: &mut SolveApp, ctx: &egui::Context) {
    match hooks() {
        Some(hooks) => (hooks.tick)(app, ctx),
        None => app.session.poll_multiphysics_jobs(),
    }
}

pub fn reset_document(app: &mut SolveApp) {
    match hooks() {
        Some(hooks) => (hooks.reset_document)(app),
        None => app.multiphysics_workspace = Slot::default(),
    }
}

/// Finish the workspace's work before the native host tears down its resources.
pub fn close(app: &mut SolveApp) {
    if let Some(hooks) = hooks() {
        (hooks.close)(app);
    }
}

pub(crate) fn take_view(app: &mut SolveApp) -> Slot {
    hooks().map_or_else(Slot::default, |hooks| (hooks.take_view)(app))
}

pub(crate) fn put_view(app: &mut SolveApp, view: Slot) {
    if let Some(hooks) = hooks() {
        (hooks.put_view)(app, view);
    }
}
