//! Numerical studies with independent authoring, jobs and post-processing.
//! The ordinary Simulation placeholder and Design/Render workspace state are not reused.

use std::collections::BTreeMap;

use egui::{RichText, Sense};
use serde_json::{Value, json};
use solidcraft_engine::multiphysics::results as post;
use solidcraft_engine::multiphysics::{Study, Template};
use solidcraft_theme::Tokens;

use crate::SolveApp;
use crate::ui_kit::{BackingTooltipExt as _, Container, Edge, LayoutMode, Surface};

#[path = "multiphysics_shell.rs"]
mod shell;

const TAB: &str = "MULTIPHYSICS";
const MAX_TABLE_ROWS: usize = 2048;
const MAX_TABLE_COLUMNS: usize = 24;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Stage {
    #[default]
    Prepare,
    Physics,
    Mesh,
    Couple,
    Run,
    Results,
    Verify,
}

#[derive(Clone, Debug)]
pub struct Draft {
    pub name: String,
    pub question: String,
    pub case_text: String,
    pub expected_revision: u64,
    base_name: String,
    base_question: String,
    base_case: Value,
    // Edits can bypass widgets, so the complete text is the cache key.
    case_dirty_cache: std::cell::RefCell<Option<(String, bool)>>,
}
impl Draft {
    fn from_study(study: &Study) -> Self {
        let case_text = pretty(&study.case);
        Self {
            name: study.name.clone(),
            question: study.question.clone(),
            case_dirty_cache: std::cell::RefCell::new(Some((case_text.clone(), false))),
            case_text,
            expected_revision: study.revision,
            base_name: study.name.clone(),
            base_question: study.question.clone(),
            base_case: study.case.clone(),
        }
    }
    pub fn dirty(&self) -> bool {
        if self.name != self.base_name || self.question != self.base_question {
            return true;
        }
        if let Ok(cache) = self.case_dirty_cache.try_borrow()
            && let Some((text, dirty)) = cache.as_ref()
            && text == &self.case_text
        {
            return *dirty;
        }
        let dirty = serde_json::from_str::<Value>(&self.case_text).ok().as_ref() != Some(&self.base_case);
        if let Ok(mut cache) = self.case_dirty_cache.try_borrow_mut() {
            *cache = Some((self.case_text.clone(), dirty));
        }
        dirty
    }
}

/// Per-document presentation and unapplied edits; the document-tab owner retains it.
#[derive(Clone, Debug)]
pub struct ViewState {
    pub selected: Option<String>,
    pub stage: Stage,
    pub drafts: BTreeMap<String, Draft>,
    json_editor: bool,
    case_path: String,
    export_directory: String,
    selected_run: Option<String>,
    shown_run: Option<String>,
    selected_artifact: usize,
    x_column: usize,
    y_column: usize,
    slice_axis: usize,
    slice_layer: usize,
    slice_field: usize,
    color_map: post::ColorMap,
    /// A user-locked colour range shared by 2D and 3D field views and by PVD frames.
    locked_range: Option<(f64, f64)>,
    probe_field: String,
    probe_point: [f64; 3],
    line_end: [f64; 3],
    line_samples: usize,
    compare_run: Option<String>,
    goal_index: Option<usize>,
    #[cfg(not(target_arch = "wasm32"))]
    show_3d: bool,
    #[cfg(not(target_arch = "wasm32"))]
    field_3d: crate::multiphysics_field_view::State3d,
    error: Option<String>,
    notice: Option<String>,
}
impl Default for ViewState {
    fn default() -> Self {
        Self {
            selected: None,
            stage: Stage::default(),
            drafts: BTreeMap::new(),
            json_editor: false,
            case_path: String::new(),
            export_directory: String::new(),
            selected_run: None,
            shown_run: None,
            selected_artifact: 0,
            x_column: 0,
            y_column: 1,
            slice_axis: 2,
            slice_layer: 0,
            slice_field: 0,
            color_map: post::ColorMap::default(),
            locked_range: None,
            probe_field: String::new(),
            probe_point: [0.0; 3],
            line_end: [1.0, 0.0, 0.0],
            line_samples: 100,
            compare_run: None,
            goal_index: None,
            #[cfg(not(target_arch = "wasm32"))]
            show_3d: false,
            #[cfg(not(target_arch = "wasm32"))]
            field_3d: Default::default(),
            error: None,
            notice: None,
        }
    }
}
/// Window-wide catalogue caches and a single export slot stay outside document-tab state.
#[derive(Default)]
pub struct State {
    pub view: ViewState,
    preflight: Option<PreflightCache>,
    templates: Option<Result<Vec<Template>, String>>,
    catalog: Option<Result<Value, String>>,
    template_index: usize,
    search: String,
    show_catalog: bool,
    show_materials: bool,
    materials: crate::multiphysics_material_library::State,
    show_browser: bool,
    #[cfg(not(target_arch = "wasm32"))]
    export: Option<ExportJob>,
    #[cfg(not(target_arch = "wasm32"))]
    inspection: Option<InspectionJob>,
    #[cfg(not(target_arch = "wasm32"))]
    inspected: Option<Inspected>,
    /// Post-processing answers for the shown run; dropped when the project or run changes.
    post: Option<PostCache>,
    /// The Analysis tree node picked for `stage` (see [`shell::current`]).
    node: Option<(Stage, shell::Node)>,
    /// The picked node's feature dialog or result window is open over the graphics window.
    dialog_open: bool,
    show_monitor: bool,
    wizard: Option<shell::Wizard>,
    ribbon_tab: usize,
    /// Input Data and Results folders folded in the Analysis tree.
    folders_closed: [bool; 2],
    /// Other tree folders folded, by node key.
    tree_folded: std::collections::BTreeSet<&'static str>,
    /// Customize Tree's choices, by node key; absent means the folder's default.
    tree_custom: BTreeMap<&'static str, bool>,
    show_customize: bool,
    /// The case entry (JSON pointer) picked in the tree, edited alone in its feature dialog.
    entry: Option<String>,
}

/// Results derived from one retained run. The run is immutable evidence, so answers stay
/// valid while the project `Arc` and the run id are unchanged.
struct PostCache {
    project: std::sync::Arc<solidcraft_engine::multiphysics::Project>,
    study: String,
    run: String,
    goals: Option<Result<post::GoalReport, String>>,
    components: Option<Result<Vec<post::ComponentRow>, String>>,
    comparison: Option<Result<post::MetricTable, String>>,
    field: Option<Result<post::QueryResponse, String>>,
}
struct PreflightCache {
    project: std::sync::Weak<solidcraft_engine::multiphysics::Project>,
    study_id: String,
    validation: std::sync::Arc<solidcraft_engine::multiphysics::Validation>,
}
impl State {
    fn study_preflight(
        &mut self,
        project: &std::sync::Arc<solidcraft_engine::multiphysics::Project>,
        study: &Study,
    ) -> std::sync::Arc<solidcraft_engine::multiphysics::Validation> {
        if let Some(cached) = &self.preflight
            && cached.study_id == study.id
            && let Some(previous) = cached.project.upgrade()
            && std::sync::Arc::ptr_eq(&previous, project)
        {
            return cached.validation.clone();
        }
        let validation = std::sync::Arc::new(solidcraft_engine::multiphysics::preflight(study));
        self.preflight =
            Some(PreflightCache { project: std::sync::Arc::downgrade(project), study_id: study.id.clone(), validation: validation.clone() });
        validation
    }
}
impl std::ops::Deref for State {
    type Target = ViewState;
    fn deref(&self) -> &Self::Target {
        &self.view
    }
}
impl std::ops::DerefMut for State {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.view
    }
}
#[cfg(not(target_arch = "wasm32"))]
struct ExportJob {
    uid: u64,
    command: &'static str,
    done: std::sync::mpsc::Receiver<Result<Value, String>>,
}

/// The project pointer is an immutable snapshot. A changed project retires its inspection.
#[cfg(not(target_arch = "wasm32"))]
struct InspectionKey {
    uid: u64,
    project: std::sync::Arc<solidcraft_engine::multiphysics::Project>,
    study: String,
    run: String,
    artifact: String,
}
#[cfg(not(target_arch = "wasm32"))]
impl InspectionKey {
    fn matches(&self, app: &SolveApp, study: &str, run: &str, artifact: &str) -> bool {
        self.uid == app.session.uid
            && std::sync::Arc::ptr_eq(&self.project, &app.session.doc.multiphysics)
            && self.study == study
            && self.run == run
            && self.artifact == artifact
    }
}
#[cfg(not(target_arch = "wasm32"))]
struct InspectionJob {
    key: InspectionKey,
    cancelled: std::sync::Arc<std::sync::atomic::AtomicBool>,
    done: std::sync::mpsc::Receiver<Result<InspectionOutput, String>>,
}
#[cfg(not(target_arch = "wasm32"))]
struct InspectionOutput {
    metadata: Value,
    spatial: Option<Result<solidcraft_engine::multiphysics::vtk::BoxGrid, String>>,
    frame_coordinate: Option<f64>,
}
#[cfg(not(target_arch = "wasm32"))]
impl From<Value> for InspectionOutput {
    fn from(metadata: Value) -> Self {
        Self { metadata, spatial: None, frame_coordinate: None }
    }
}
#[cfg(not(target_arch = "wasm32"))]
struct Inspected {
    key: InspectionKey,
    outcome: Result<InspectionOutput, String>,
}
impl State {
    fn export_busy(&self) -> bool {
        #[cfg(not(target_arch = "wasm32"))]
        {
            self.export.is_some()
        }
        #[cfg(target_arch = "wasm32")]
        {
            false
        }
    }
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
    app.session.active_sketch.is_none() && app.tool.is_none() && app.dialog.is_none() && !app.file_load.busy() && !app.file_picker_busy()
}

/// UI workspace routing is separate from numerical document commands.
pub fn switch(app: &mut SolveApp, workspace: &str) -> Result<(), String> {
    match workspace {
        "multiphysics" => {
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

pub fn reset_document(app: &mut SolveApp) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        // Keep the one export slot occupied until retirement, even across New/Open.
        // Its captured UID prevents a completion from becoming the new document's notice.
        let export = app.multiphysics_workspace.export.take();
        let inspection = app.multiphysics_workspace.inspection.take();
        if let Some(job) = &inspection {
            job.cancelled.store(true, std::sync::atomic::Ordering::Release);
        }
        app.multiphysics_workspace = State { export, inspection, ..Default::default() };
    }
    #[cfg(target_arch = "wasm32")]
    {
        app.multiphysics_workspace = Default::default();
    }
}

/// Finish the frozen export before the native host tears down its resources.
pub fn close(app: &mut SolveApp) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let mut state = std::mem::take(&mut app.multiphysics_workspace);
        if let Some(export) = state.export.take() {
            let outcome = export.done.recv().unwrap_or_else(|error| Err(error.to_string()));
            finish_export(app, &mut state, export.uid, export.command, outcome);
        }
        if let Some(job) = state.inspection.take() {
            // Parsing is bounded, but cannot be interrupted inside the parser yet. Wait until
            // its scratch session has retired before the host tears down its resources.
            job.cancelled.store(true, std::sync::atomic::Ordering::Release);
            let _ = job.done.recv();
        }
        app.multiphysics_workspace = state;
    }
    #[cfg(target_arch = "wasm32")]
    {
        let _ = app;
    }
}

pub fn tick(app: &mut SolveApp, ctx: &egui::Context) {
    app.session.poll_multiphysics_jobs();
    #[cfg(not(target_arch = "wasm32"))]
    {
        let mut state = std::mem::take(&mut app.multiphysics_workspace);
        poll_export(app, &mut state);
        poll_inspection(app, &mut state);
        app.multiphysics_workspace = state;
    }
    if app.multiphysics_workspace.export_busy()
        || inspection_busy(&app.multiphysics_workspace)
        || app.session.multiphysics_jobs().iter().any(|job| matches!(job.state.as_str(), "queued" | "running" | "cancelling"))
    {
        ctx.request_repaint_after(std::time::Duration::from_millis(100));
    }
}

fn tr(key: &str) -> &str {
    solidcraft_ui_text::tr(key)
}
fn pretty(value: &Value) -> String {
    serde_json::to_string_pretty(value).unwrap_or_default()
}
fn button(ui: &mut egui::Ui, key: &str, handle: &str, enabled: bool) -> egui::Response {
    let response = crate::ui_kit::button(ui, tr(key), crate::ui_kit::Button::Enabled(enabled));
    crate::scenario::publish_target(handle, response.rect);
    response
}
fn command(app: &mut SolveApp, state: &mut State, id: &str, params: Value) -> Option<Value> {
    #[cfg(not(target_arch = "wasm32"))]
    if let Some(id) = match id {
        "multiphysics.case.export" => Some("multiphysics.case.export"),
        "multiphysics.result.export" => Some("multiphysics.result.export"),
        _ => None,
    } {
        return start_export(app, state, id, params);
    }
    match app.run(id, params) {
        Ok(reply) => {
            state.error = None;
            state.notice = if matches!(id, "multiphysics.case.export" | "multiphysics.result.export") {
                let key = format!("cmd.{id}");
                Some(solidcraft_ui_text::format("ui.command_done", &[("command", solidcraft_ui_text::tr(&key))]))
            } else {
                None
            };
            Some(reply)
        }
        Err(error) => {
            state.error = Some(error);
            state.notice = None;
            None
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn start_export(app: &SolveApp, state: &mut State, id: &'static str, params: Value) -> Option<Value> {
    if app.file_load.busy() || app.file_picker_busy() {
        state.error = Some(tr("a file operation is in progress").into());
        return None;
    }
    if state.export_busy() {
        state.error = Some(tr("multiphysics.export_busy").into());
        return None;
    }
    let uid = app.session.uid;
    let mut copy = app.session.scratch();
    let scope = solidcraft_engine::io::vfs::current_scope();
    let (send, done) = std::sync::mpsc::channel();
    let worker = std::thread::Builder::new().name("multiphysics-export".into()).spawn(move || {
        let outcome = solidcraft_engine::io::vfs::scoped(scope, || copy.execute(id, &params).map_err(|error| error.to_string()));
        // Retire engine/device references before acknowledging completion or host shutdown.
        drop(copy);
        let _ = send.send(outcome);
    });
    match worker {
        Ok(_) => {
            state.export = Some(ExportJob { uid, command: id, done });
            state.error = None;
            state.notice = None;
            Some(json!({"queued":true}))
        }
        Err(error) => {
            state.error = Some(error.to_string());
            state.notice = None;
            None
        }
    }
}
#[cfg(not(target_arch = "wasm32"))]
fn poll_export(app: &mut SolveApp, state: &mut State) {
    let Some(export) = &state.export else { return };
    let outcome = match export.done.try_recv() {
        Ok(outcome) => outcome,
        Err(std::sync::mpsc::TryRecvError::Empty) => return,
        Err(error) => Err(error.to_string()),
    };
    let Some(export) = state.export.take() else { return };
    finish_export(app, state, export.uid, export.command, outcome);
}
#[cfg(not(target_arch = "wasm32"))]
fn finish_export(app: &mut SolveApp, state: &mut State, uid: u64, id: &'static str, outcome: Result<Value, String>) {
    if app.session.uid != uid {
        return;
    }
    app.finish_command(id, &outcome);
    match outcome {
        Ok(_) => {
            let key = format!("cmd.{id}");
            state.notice = Some(solidcraft_ui_text::format("ui.command_done", &[("command", tr(&key))]));
            state.error = None;
        }
        Err(error) => {
            state.error = Some(error);
            state.notice = None;
        }
    }
}

fn inspection_busy(state: &State) -> bool {
    #[cfg(not(target_arch = "wasm32"))]
    {
        state.inspection.is_some()
    }
    #[cfg(target_arch = "wasm32")]
    {
        let _ = state;
        false
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn start_inspection(app: &SolveApp, state: &mut State, study: &str, run: &str, artifact: &str) {
    use std::sync::atomic::{AtomicBool, Ordering};
    if state.inspection.is_some() {
        return;
    }
    let key = InspectionKey {
        uid: app.session.uid,
        project: app.session.doc.multiphysics.clone(),
        study: study.into(),
        run: run.into(),
        artifact: artifact.into(),
    };
    let params = json!({"id":study,"run_id":run,"artifact":artifact});
    let source = key.project.clone();
    let source_study = study.to_string();
    let source_run = run.to_string();
    let source_artifact = artifact.to_string();
    let mut scratch = app.session.scratch();
    let cancelled = std::sync::Arc::new(AtomicBool::new(false));
    let worker_cancelled = cancelled.clone();
    let (send, done) = std::sync::mpsc::channel();
    let worker = std::thread::Builder::new().name("multiphysics-inspection".into()).spawn(move || {
        let result = if worker_cancelled.load(Ordering::Acquire) {
            Err(String::new())
        } else {
            scratch.execute("multiphysics.result.inspect", &params).map_err(|error| error.to_string()).map(|metadata| {
                let spatial = if source_artifact.ends_with(".vtu") && !worker_cancelled.load(Ordering::Acquire) {
                    Some((|| {
                        let text = source
                            .study(&source_study)
                            .and_then(|study| study.runs.iter().find(|run| run.id == source_run))
                            .and_then(|run| run.artifacts.iter().find(|artifact| artifact.name == source_artifact))
                            .map(|artifact| artifact.text.as_str())
                            .ok_or_else(|| "Retained field artifact is unavailable".to_string())?;
                        let grid = solidcraft_engine::multiphysics::vtk::parse_vtu(text).map_err(|error| error.to_string())?;
                        solidcraft_engine::multiphysics::vtk::BoxGrid::from_grid(grid).map_err(|error| error.to_string())
                    })())
                } else {
                    None
                };
                let frame_coordinate = source
                    .study(&source_study)
                    .and_then(|study| study.runs.iter().find(|run| run.id == source_run))
                    .and_then(|run| retained_frame_coordinate(run, &source_artifact));
                InspectionOutput { metadata, spatial, frame_coordinate }
            })
        };
        drop(scratch);
        if !worker_cancelled.load(Ordering::Acquire) {
            let _ = send.send(result);
        }
    });
    match worker {
        Ok(_) => {
            // The previous answer stays until this one lands (it is shown only for its own
            // artifact), so frame stepping can keep the chosen field and camera.
            state.inspection = Some(InspectionJob { key, cancelled, done });
        }
        Err(error) => state.error = Some(error.to_string()),
    }
}

/// A PVD coordinate may mean physical time or a steady iteration, depending on the producer.
#[cfg(not(target_arch = "wasm32"))]
fn retained_frame_coordinate(run: &solidcraft_engine::multiphysics::RunEvidence, file: &str) -> Option<f64> {
    let mut coordinate = None;
    for artifact in run.artifacts.iter().filter(|artifact| artifact.name.ends_with(".pvd")) {
        for frame in solidcraft_engine::multiphysics::vtk::parse_pvd(&artifact.text).ok()? {
            if frame.file == file {
                if coordinate.is_some_and(|old| old != frame.coordinate) {
                    return None;
                }
                coordinate = Some(frame.coordinate);
            }
        }
    }
    coordinate
}

#[cfg(not(target_arch = "wasm32"))]
fn inspection_is_selected(app: &SolveApp, state: &State, key: &InspectionKey) -> bool {
    if state.stage != Stage::Results || state.show_catalog {
        return false;
    }
    let Some(study) = state.selected.as_ref().and_then(|id| app.session.doc.multiphysics.study(id)) else { return false };
    let run = study.runs.iter().find(|run| Some(run.id.as_str()) == state.selected_run.as_deref()).or_else(|| study.runs.last());
    let Some(run) = run else { return false };
    let Some(artifact) = run.artifacts.get(state.selected_artifact) else { return false };
    key.matches(app, &study.id, &run.id, &artifact.name)
}

#[cfg(not(target_arch = "wasm32"))]
fn poll_inspection(app: &SolveApp, state: &mut State) {
    if state
        .inspected
        .as_ref()
        .is_some_and(|cached| cached.key.uid != app.session.uid || !std::sync::Arc::ptr_eq(&cached.key.project, &app.session.doc.multiphysics))
    {
        state.inspected = None;
    }
    let Some(job) = &state.inspection else { return };
    let outcome = match job.done.try_recv() {
        Ok(outcome) => outcome,
        Err(std::sync::mpsc::TryRecvError::Empty) => return,
        Err(error) => Err(error.to_string()),
    };
    let Some(job) = state.inspection.take() else { return };
    if !job.cancelled.load(std::sync::atomic::Ordering::Acquire) && inspection_is_selected(app, state, &job.key) {
        // Stepping through frames of one run keeps the chosen field, plane and camera.
        let same_run = state.inspected.as_ref().is_some_and(|old| old.key.study == job.key.study && old.key.run == job.key.run);
        if !same_run {
            state.slice_field = 0;
            state.slice_layer = 0;
            state.field_3d = Default::default();
        }
        state.inspected = Some(Inspected { key: job.key, outcome });
    }
}

/// Own docked composition. No Design Browser, History, Render inspector or camera takeover.
pub fn show(app: &mut SolveApp, ui: &mut egui::Ui) {
    let mut state = std::mem::take(&mut app.multiphysics_workspace);
    if state.templates.is_none() {
        state.templates = Some(solidcraft_engine::multiphysics::templates().map_err(|error| error.to_string()));
        state.show_browser = true;
        state.y_column = 1;
    }
    if state.catalog.is_none() {
        state.catalog = Some(solidcraft_engine::multiphysics::catalog().map_err(|error| error.to_string()));
    }
    let document = app.session.doc.clone();
    let studies = &document.multiphysics.studies;
    if state.selected.as_ref().is_none_or(|id| !studies.iter().any(|study| &study.id == id)) {
        state.selected = studies.first().map(|study| study.id.clone());
    }
    for study in studies {
        let draft = state.drafts.entry(study.id.clone()).or_insert_with(|| Draft::from_study(study));
        if draft.expected_revision != study.revision && !draft.dirty() {
            *draft = Draft::from_study(study);
        }
    }
    state.drafts.retain(|id, _| studies.iter().any(|study| &study.id == id));
    let t = Tokens::get();
    Container::for_layout(LayoutMode::Docked, Surface::Chrome)
        .edge("multiphysics_ribbon", Edge::Top)
        .exact_size(shell::RIBBON)
        .resizable(false)
        .show(ui, |ui| shell::ribbon_bar(app, &mut state, studies, ui));
    let selected = state.selected.as_ref().and_then(|id| studies.iter().find(|study| &study.id == id));
    Container::for_layout(LayoutMode::Docked, Surface::Panel)
        .edge("multiphysics_status", Edge::Bottom)
        .exact_size(t.control_height + t.space_4)
        .resizable(false)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                let viewing_retained = state.stage == Stage::Results && selected.is_some_and(|study| !study.runs.is_empty());
                ui.label(tr(if state.show_materials {
                    "multiphysics.materials_status"
                } else if app.session.multiphysics_host_available() {
                    "multiphysics.worker_ready"
                } else if viewing_retained {
                    "multiphysics.retained_view"
                } else {
                    "multiphysics.preparation_only"
                }));
                if state.export_busy() {
                    ui.spinner();
                    ui.add(egui::Label::new(tr("multiphysics.export_busy")).truncate()).backing_on_hover_text(tr("multiphysics.export_busy"));
                }
                if let Some(error) = &state.error {
                    ui.add(egui::Label::new(RichText::new(error).color(t.form_error)).truncate()).backing_on_hover_text(error);
                } else if let Some(notice) = &state.notice {
                    ui.add(egui::Label::new(RichText::new(notice).color(t.form_success)).truncate());
                }
            });
        });
    Container::for_layout(LayoutMode::Docked, Surface::Panel)
        .edge("multiphysics_resource_bar", Edge::Left)
        .exact_size(shell::RESOURCE_BAR)
        .resizable(false)
        .show(ui, |ui| shell::resource_bar(&mut state, ui));
    if state.show_browser && !state.show_materials {
        Container::for_layout(LayoutMode::Docked, Surface::Panel)
            .edge("multiphysics_browser", Edge::Left)
            .default_size(264.0)
            .size_range(210.0..=420.0)
            .resizable(true)
            .show(ui, |ui| shell::tree(app, &mut state, studies, ui));
    }
    // The ribbon and the tree may have changed the active project this frame.
    let selected = state.selected.as_ref().and_then(|id| studies.iter().find(|study| &study.id == id));
    let report = selected.map(|study| state.study_preflight(&document.multiphysics, study));
    let mut graphics = egui::Rect::NOTHING;
    Container::for_layout(LayoutMode::Docked, Surface::Viewport).central().show(ui, |ui| {
        graphics = ui.max_rect();
        shell::graphics(app, &mut state, selected, ui);
    });
    let ctx = ui.ctx().clone();
    shell::windows(app, &mut state, selected, report.as_deref(), graphics, &ctx);
    app.multiphysics_workspace = state;
}

#[derive(Debug, PartialEq, Eq)]
struct Readiness {
    case_valid: bool,
    draft_dirty: bool,
    draft_conflict: bool,
    worker_available: bool,
    current_result: bool,
    historical_result: bool,
}

impl Readiness {
    fn for_study(study: &Study, draft: Option<&Draft>, worker_available: bool, case_valid: bool) -> Self {
        Self {
            case_valid,
            draft_dirty: draft.is_some_and(Draft::dirty),
            draft_conflict: draft.is_some_and(|draft| draft.expected_revision != study.revision),
            worker_available,
            current_result: study.runs.iter().any(|run| run.study_revision == study.revision),
            historical_result: study.runs.iter().any(|run| run.study_revision != study.revision),
        }
    }

    fn can_run(&self) -> bool {
        self.case_valid && !self.draft_dirty && !self.draft_conflict && self.worker_available
    }
}

fn prepare(app: &mut SolveApp, state: &mut State, study: &Study, report: &solidcraft_engine::multiphysics::Validation, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical().id_salt("multiphysics_prepare").show(ui, |ui| {
        let readiness = Readiness::for_study(study, state.drafts.get(&study.id), app.session.multiphysics_host_available(), report.errors.is_empty());
        ui.heading(tr("multiphysics.committed_preflight"));
        if readiness.case_valid {
            ui.colored_label(Tokens::get().form_success, tr("multiphysics.case_ready"));
        } else {
            for error in report.errors.iter().take(12) {
                ui.colored_label(Tokens::get().form_error, error);
            }
        }
        if readiness.draft_conflict {
            ui.colored_label(Tokens::get().warning, tr("multiphysics.conflict"));
        } else if readiness.draft_dirty {
            ui.colored_label(Tokens::get().warning, tr("multiphysics.apply_before_run"));
        }
        ui.label(tr(if readiness.worker_available { "multiphysics.worker_ready" } else { "multiphysics.preparation_only" }));
        if let Some(run) = study.runs.iter().rev().find(|run| run.study_revision == study.revision) {
            ui.horizontal(|ui| {
                ui.label(tr("multiphysics.stage_results"));
                ui.monospace(&run.id);
            });
        }
        if readiness.historical_result {
            ui.colored_label(Tokens::get().warning, tr("multiphysics.stale_run"));
        }
        if !readiness.current_result && !readiness.historical_result {
            ui.label(tr("multiphysics.no_results"));
        }
        ui.horizontal_wrapped(|ui| {
            for (key, target) in [
                ("multiphysics.stage_physics", Stage::Physics),
                ("multiphysics.stage_run", Stage::Run),
                ("multiphysics.stage_results", Stage::Results),
                ("multiphysics.stage_verify", Stage::Verify),
            ] {
                let response = crate::ui_kit::button(ui, tr(key), crate::ui_kit::Button::Normal);
                crate::scenario::publish_target(&format!("multiphysics:readiness:{key}"), response.rect);
                if response.clicked() {
                    state.stage = target;
                }
            }
        });
        ui.separator();
        ui.heading(tr("multiphysics.intent"));
        if let Some(draft) = state.drafts.get_mut(&study.id) {
            ui.label(tr("multiphysics.question"));
            let response = crate::ui_kit::text_area(
                ui,
                tr("multiphysics.question"),
                egui::Id::new(("multiphysics_question", &study.id)),
                &mut draft.question,
                3,
                4096,
            );
            crate::scenario::publish_target("multiphysics:question", response.rect);
        }
        ui.add_space(Tokens::get().space_12);
        ui.strong(tr("multiphysics.representation"));
        ui.label(tr("multiphysics.case_scope"));
        ui.label(tr("multiphysics.cad_limit"));
        ui.separator();
        ui.label(tr("multiphysics.units"));
        ui.label(tr("multiphysics.reproducible"));
        if button(ui, "multiphysics.edit_settings", "multiphysics:edit_settings", true).clicked() {
            state.stage = Stage::Physics;
        }
        ui.separator();
        ui.strong(tr("multiphysics.case_export"));
        ui.label(tr("multiphysics.case_path"));
        crate::ui_kit::text_row(ui, tr("multiphysics.path"), &mut state.case_path, "");
        if button(ui, "multiphysics.export_case", "multiphysics:export_case", !state.case_path.trim().is_empty() && !state.export_busy())
            .backing_on_hover_text(tr(if state.export_busy() { "multiphysics.export_busy" } else { "multiphysics.case_path" }))
            .clicked()
        {
            let path = state.case_path.clone();
            command(app, state, "multiphysics.case.export", json!({"id":study.id,"path":path}));
        }
    });
}

fn inspector(
    app: &mut SolveApp,
    state: &mut State,
    study: Option<&Study>,
    report: Option<&solidcraft_engine::multiphysics::Validation>,
    ui: &mut egui::Ui,
) {
    let (Some(study), Some(report)) = (study, report) else { return };
    if state.stage == Stage::Results && !study.runs.is_empty() {
        result_summary(state, study, ui);
        return;
    }
    ui.heading(tr("multiphysics.inspector"));
    if let Some(draft) = state.drafts.get_mut(&study.id) {
        crate::ui_kit::text_row(ui, tr("multiphysics.name"), &mut draft.name, "");
        ui.label(RichText::new(&study.id).small().color(Tokens::get().text_dim));
        ui.label(solidcraft_ui_text::format("multiphysics.revision", &[("revision", &study.revision.to_string())]));
        if draft.expected_revision != study.revision {
            ui.colored_label(Tokens::get().warning, tr("multiphysics.conflict"));
        }
    }
    let draft_dirty = state.drafts.get(&study.id).is_some_and(Draft::dirty);
    ui.horizontal(|ui| {
        if button(ui, "multiphysics.apply", "multiphysics:apply", draft_dirty).clicked() {
            apply(app, state, study);
        }
        if button(ui, "multiphysics.cancel_draft", "multiphysics:cancel_draft", draft_dirty).clicked() {
            state.drafts.insert(study.id.clone(), Draft::from_study(study));
            state.error = None;
        }
    });
    ui.separator();
    ui.strong(tr("multiphysics.committed_preflight"));
    if report.errors.is_empty() {
        ui.colored_label(Tokens::get().form_success, tr("multiphysics.case_ready"));
    }
    for message in report.errors.iter().take(12) {
        ui.colored_label(Tokens::get().form_error, message);
    }
    for message in report.warnings.iter().take(8) {
        ui.colored_label(Tokens::get().warning, message);
    }
    if draft_dirty {
        ui.label(tr("multiphysics.apply_before_run"));
        ui.label(RichText::new(tr("multiphysics.draft_session_only")).small().color(Tokens::get().text_dim));
    }
    let available = app.session.multiphysics_host_available();
    let ready = Readiness::for_study(study, state.drafts.get(&study.id), available, report.errors.is_empty()).can_run();
    if button(ui, "multiphysics.run", "multiphysics:run", ready)
        .backing_on_hover_text(tr(if available { "multiphysics.run_hint" } else { "multiphysics.preparation_only" }))
        .clicked()
    {
        if command(app, state, "multiphysics.run.start", json!({"id":study.id})).is_some() {
            state.stage = Stage::Run;
        }
    }
    ui.separator();
    ui.label(RichText::new(tr("multiphysics.artifact_hint")).small().color(Tokens::get().text_dim));
}

fn result_summary(state: &State, study: &Study, ui: &mut egui::Ui) {
    let run = study.runs.iter().find(|run| Some(run.id.as_str()) == state.selected_run.as_deref()).or_else(|| study.runs.last());
    let Some(run) = run else { return };
    ui.heading(tr("multiphysics.result_summary"));
    ui.monospace(&run.id);
    ui.label(solidcraft_ui_text::format("multiphysics.revision", &[("revision", &run.study_revision.to_string())]));
    if run.study_revision != study.revision {
        ui.colored_label(Tokens::get().warning, tr("multiphysics.stale_run"));
    }
    ui.separator();
    // Every numeric metric the solver reported, not a fixed subset: domain cases each
    // report their own quantities (forces, eigenvalues, component temperatures, …).
    let known = [
        ("time_s", "multiphysics.result_time"),
        ("fluid_cell_count", "multiphysics.result_fluid_cells"),
        ("excluded_solid_cell_count", "multiphysics.result_excluded_cells"),
        ("relative_boundary_flux_imbalance", "multiphysics.result_flux_imbalance"),
        ("continuity_l2", "multiphysics.result_continuity"),
    ];
    if let Some(metrics) = run.manifest.get("metrics").and_then(Value::as_object) {
        egui::Grid::new("multiphysics_summary_metrics").striped(true).show(ui, |ui| {
            for (metric, value) in metrics.iter().take(MAX_TABLE_ROWS) {
                let Some(value) = value.as_f64() else { continue };
                match known.iter().find(|(name, _)| name == metric) {
                    Some((_, label)) => ui.label(tr(label)),
                    None => ui.monospace(metric),
                };
                let text = if metric.ends_with("_count") && value.fract() == 0.0 { format!("{value:.0}") } else { format!("{value:.6e}") };
                ui.monospace(text);
                ui.label(post::unit_hint(metric).unwrap_or(""));
                ui.end_row();
            }
        });
    }
    ui.separator();
    for (key, label) in [("linear_solver_converged", "multiphysics.result_linear"), ("steady_state_converged", "multiphysics.result_steady")] {
        ui.label(tr(label));
        ui.label(tr(match run.manifest.get(key).and_then(Value::as_bool) {
            Some(true) => "multiphysics.result_converged",
            Some(false) => "multiphysics.result_not_converged",
            None => "multiphysics.result_unreported",
        }));
    }
    ui.separator();
    ui.label(tr("multiphysics.result_scope"));
}

fn apply(app: &mut SolveApp, state: &mut State, study: &Study) {
    let Some(draft) = state.drafts.get(&study.id).cloned() else { return };
    let case: Value = match serde_json::from_str(&draft.case_text) {
        Ok(case) => case,
        Err(error) => {
            state.error = Some(error.to_string());
            return;
        }
    };
    if command(
        app,
        state,
        "multiphysics.study.update",
        json!({"id":study.id,"expected_revision":draft.expected_revision,"name":draft.name,"question":draft.question,"case":case}),
    )
    .is_some()
        && let Some(saved) = app.session.doc.multiphysics.study(&study.id)
    {
        // Command replies contain bounded run summaries, not retained dataset bytes.
        // The committed document is the authority for rebuilding the applied draft.
        state.drafts.insert(saved.id.clone(), Draft::from_study(saved));
    }
}

/// The case editor for one Analysis tree node: its top-level case entries, or every entry for
/// Input Data and in the JSON editor (which always edits the whole case).
fn settings(state: &mut State, study: &Study, node: shell::Node, entry: Option<&str>, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        crate::ui_kit::checkbox(ui, tr("multiphysics.json_editor"), &mut state.json_editor);
    });
    ui.label(
        RichText::new(tr(if node == shell::Node::Mesh { "multiphysics.mesh_limit" } else { "multiphysics.settings_hint" }))
            .small()
            .color(Tokens::get().text_dim),
    );
    let json_editor = state.json_editor;
    let Some(draft) = state.drafts.get_mut(&study.id) else { return };
    let parsed = serde_json::from_str::<Value>(&draft.case_text);
    let raw = json_editor || parsed.is_err();
    egui::ScrollArea::both().id_salt(("multiphysics_settings", &study.id)).max_height(520.0).show(ui, |ui| {
        if raw {
            let response = crate::ui_kit::text_area(
                ui,
                tr("multiphysics.case_settings"),
                egui::Id::new(("multiphysics_case_json", &study.id)),
                &mut draft.case_text,
                16,
                1024 * 1024,
            );
            crate::scenario::publish_target("multiphysics:case_json", response.rect);
            if let Err(error) = &parsed {
                ui.colored_label(Tokens::get().form_error, error.to_string());
            }
        } else if let Ok(mut value) = parsed {
            let changed = match entry {
                Some(pointer) => match value.pointer_mut(pointer) {
                    Some(slot) => field_entry(ui, &shell::entry_label(pointer), slot, 0),
                    None => {
                        ui.label(RichText::new(tr("multiphysics.tree_revert_missing")).small().color(Tokens::get().text_dim));
                        false
                    }
                },
                None if node == shell::Node::Input => fields(ui, &mut value, 0),
                None => fields_group(ui, &mut value, node),
            };
            if changed {
                draft.case_text = pretty(&value);
            }
        }
    });
}

/// The case entries of one tree folder (top-level entries and the children of mixed blocks
/// such as `config` filed under it).
fn fields_group(ui: &mut egui::Ui, value: &mut Value, node: shell::Node) -> bool {
    let pointers = shell::entries(value, node);
    if pointers.is_empty() {
        ui.label(RichText::new(tr("multiphysics.tree_folder_empty")).small().color(Tokens::get().text_dim));
        return false;
    }
    let mut changed = false;
    for pointer in pointers.iter().take(80) {
        if let Some(slot) = value.pointer_mut(pointer) {
            changed |= field_entry(ui, &shell::entry_label(pointer), slot, 0);
        }
    }
    changed
}

/// One named case entry: a folder for objects and arrays, a kit control for a value.
fn field_entry(ui: &mut egui::Ui, name: &str, value: &mut Value, depth: usize) -> bool {
    let mut changed = false;
    ui.push_id(name, |ui| {
        if value.is_object() || value.is_array() {
            egui::CollapsingHeader::new(name).default_open(depth == 0).show(ui, |ui| changed |= fields(ui, value, depth + 1));
        } else {
            crate::ui_kit::row(ui, name, |ui| match value {
                Value::Bool(boolean) => changed |= crate::ui_kit::checkbox(ui, "", boolean).changed(),
                Value::Number(number) => {
                    if let Some(mut numeric) = number.as_u64() {
                        if crate::ui_kit::number_field(ui, &mut numeric, 0..=u64::MAX, 1.0, "").changed() {
                            *number = numeric.into();
                            changed = true;
                        }
                    } else if let Some(mut numeric) = number.as_i64() {
                        if crate::ui_kit::number_field(ui, &mut numeric, i64::MIN..=i64::MAX, 1.0, "").changed() {
                            *number = numeric.into();
                            changed = true;
                        }
                    } else if let Some(mut numeric) = number.as_f64() {
                        let response = crate::ui_kit::number_field(ui, &mut numeric, f64::MIN..=f64::MAX, 0.001, "");
                        if response.changed()
                            && numeric.is_finite()
                            && let Some(updated) = serde_json::Number::from_f64(numeric)
                        {
                            *number = updated;
                            changed = true;
                        }
                    }
                }
                Value::String(string) if name == "solver" || name == "schema" => {
                    ui.monospace(string.as_str());
                }
                Value::String(string) => changed |= crate::ui_kit::text_field(ui, string, "").1,
                _ => {
                    ui.label(tr("multiphysics.not_specified"));
                }
            });
        }
    });
    changed
}

/// The case is the numerical worker's explicit schema. Strings/enum discriminators stay
/// visible, and arbitrary nested arrays use the JSON editor rather than implicit defaults.
fn fields(ui: &mut egui::Ui, value: &mut Value, depth: usize) -> bool {
    if depth > 5 {
        ui.monospace(pretty(value));
        return false;
    }
    let mut changed = false;
    match value {
        Value::Object(object) => {
            for (name, value) in object.iter_mut().take(80) {
                changed |= field_entry(ui, name, value, depth);
            }
        }
        Value::Array(values) => {
            if values.len() <= 16 {
                for (index, value) in values.iter_mut().enumerate() {
                    ui.push_id(index, |ui| changed |= fields(ui, value, depth + 1));
                }
            } else {
                ui.label(solidcraft_ui_text::format("multiphysics.array_entries", &[("count", &values.len().to_string())]));
            }
        }
        _ => {
            ui.monospace(pretty(value));
        }
    }
    changed
}

fn jobs(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    ui.heading(tr("multiphysics.jobs"));
    ui.label(tr("multiphysics.run_identity"));
    let jobs = app.session.multiphysics_jobs();
    egui::ScrollArea::vertical().id_salt("multiphysics_jobs").show(ui, |ui| {
        for job in jobs.iter().filter(|job| job.study_id == study.id) {
            ui.separator();
            ui.horizontal(|ui| {
                ui.strong(&job.state);
                ui.monospace(&job.id);
            });
            ui.label(&job.message);
            ui.label(solidcraft_ui_text::format("multiphysics.revision", &[("revision", &job.study_revision.to_string())]));
            if job.study_revision != study.revision || job.stale {
                ui.colored_label(Tokens::get().warning, tr("multiphysics.stale_run"));
            }
            if matches!(job.state.as_str(), "queued" | "running" | "cancelling") {
                ui.spinner();
                if button(ui, "multiphysics.cancel_job", &format!("multiphysics:cancel_job:{}", job.id), job.state != "cancelling").clicked() {
                    command(app, state, "multiphysics.run.cancel", json!({"job_id":job.id}));
                }
            }
            let retained = study.runs.iter().any(|run| run.id == job.id);
            if job.evidence.is_some()
                && button(ui, "multiphysics.retain_result", &format!("multiphysics:keep:{}", job.id), !retained && !job.stale).clicked()
                && command(app, state, "multiphysics.run.keep", json!({"job_id":job.id})).is_some()
            {
                state.selected_run = Some(job.id.clone());
                state.selected_artifact = 0;
                state.stage = Stage::Results;
            }
        }
        if !jobs.iter().any(|job| job.study_id == study.id) {
            ui.label(tr("multiphysics.no_jobs"));
        }
    });
}

/// The category of a template name such as `"Battery · four-cell series pack cycle"`.
fn template_category(name: &str) -> &str {
    name.split_once(" · ").map_or(name, |(category, _)| category)
}

fn catalog(state: &mut State, ui: &mut egui::Ui) {
    ui.heading(tr("multiphysics.catalog"));
    ui.label(tr("multiphysics.catalog_hint"));
    crate::ui_kit::text_row(ui, tr("multiphysics.search"), &mut state.search, "");
    let search = state.search.to_lowercase();
    let templates = state.templates.as_ref().and_then(|result| result.as_ref().ok()).cloned().unwrap_or_default();
    let mut selected_template = None;
    egui::ScrollArea::vertical().id_salt("multiphysics_catalog").show(ui, |ui| match state.catalog.as_ref() {
        Some(Ok(catalog)) => {
            if let Some(profiles) = catalog.get("profiles").and_then(Value::as_array) {
                for profile in profiles.iter().filter(|profile| profile.to_string().to_lowercase().contains(&search)) {
                    let name = profile.get("name").and_then(Value::as_str).unwrap_or_default();
                    egui::CollapsingHeader::new(name).id_salt(profile.get("id").and_then(Value::as_str).unwrap_or(name)).show(ui, |ui| {
                        for key in ["lens", "source_scope", "representation", "limitation"] {
                            if let Some(text) = profile.get(key).and_then(Value::as_str) {
                                ui.label(text);
                            }
                        }
                        ui.colored_label(Tokens::get().warning, tr("multiphysics.catalog_execution"));
                        let id = profile.get("id").and_then(Value::as_str).unwrap_or_default();
                        let runnable: Vec<(usize, String)> = templates
                            .iter()
                            .enumerate()
                            .filter(|(_, template)| template.profile == id)
                            .map(|(index, template)| (index, template.name.clone()))
                            .collect();
                        if !runnable.is_empty() {
                            ui.strong(tr("multiphysics.catalog_templates"));
                            for (index, name) in runnable {
                                ui.horizontal(|ui| {
                                    ui.label(&name);
                                    let response = crate::ui_kit::button(ui, tr("multiphysics.catalog_use_template"), crate::ui_kit::Button::Small);
                                    crate::scenario::publish_target(&format!("multiphysics:catalog:template:{index}"), response.rect);
                                    if response.clicked() {
                                        selected_template = Some((index, template_category(&name).to_owned()));
                                    }
                                });
                            }
                        }
                        if let Some(url) = profile.get("source_url").and_then(Value::as_str) {
                            ui.hyperlink_to(tr("multiphysics.source"), url);
                        }
                    });
                }
            }
        }
        Some(Err(error)) => {
            ui.colored_label(Tokens::get().form_error, error);
        }
        None => {
            ui.spinner();
        }
    });
    // Choosing a runnable template opens the wizard on its summary, ready to finish.
    if let Some((index, category)) = selected_template {
        state.template_index = index;
        state.wizard = Some(shell::Wizard::from_template(index, category));
        state.show_catalog = false;
    }
}

fn verify(study: &Study, ui: &mut egui::Ui) {
    ui.heading(tr("multiphysics.verification"));
    ui.label(tr("multiphysics.verification_hint"));
    ui.label(tr("multiphysics.validation_boundary"));
    egui::ScrollArea::vertical().id_salt("multiphysics_verification").show(ui, |ui| {
        for run in &study.runs {
            ui.separator();
            ui.monospace(&run.id);
            ui.label(&run.case_hash);
            if run.study_revision != study.revision {
                ui.colored_label(Tokens::get().warning, tr("multiphysics.stale_run"));
            }
            for key in [
                "solver",
                "rust_port_version",
                "upstream_revision",
                "backend",
                "execution_completed",
                "linear_solver_converged",
                "steady_state_converged",
            ] {
                if let Some(value) = run.manifest.get(key) {
                    ui.horizontal(|ui| {
                        ui.label(key);
                        ui.monospace(value.to_string());
                    });
                }
            }
        }
    });
}

fn results(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    if study.runs.is_empty() {
        ui.label(tr("multiphysics.no_results"));
        return;
    }
    let selected = state.selected_run.as_deref();
    let run = study.runs.iter().find(|run| Some(run.id.as_str()) == selected).or_else(|| study.runs.last());
    let Some(run) = run else { return };
    if state.shown_run.as_deref() != Some(&run.id) {
        state.shown_run = Some(run.id.clone());
        state.selected_artifact = run.artifacts.iter().position(|artifact| artifact.name.ends_with(".csv")).unwrap_or(0);
    }
    let current_id = run.id.clone();
    crate::ui_kit::dropdown_menu(ui, "multiphysics_run", &current_id, |ui| {
        for run in &study.runs {
            if crate::ui_kit::list_picker(
                ui,
                &format!("multiphysics_run_{}", run.id),
                &run.id,
                state.selected_run.as_deref() == Some(run.id.as_str()),
            )
            .clicked()
            {
                state.selected_run = Some(run.id.clone());
            }
        }
    });
    if run.study_revision != study.revision {
        ui.colored_label(Tokens::get().warning, tr("multiphysics.stale_run"));
    }
    // The Analysis tree node decides which of the run's tools show: the Results folder shows
    // them all, a feature (cut plot, goal plot, report, …) only its own.
    let focus = shell::current(state);
    let all = focus == shell::Node::Results;
    if focus == shell::Node::Report {
        result_summary(state, study, ui);
        ui.separator();
    }
    post_processing(app, state, study, run, focus, ui);
    if !(all || focus == shell::Node::Report) {
        return results_artifacts(app, state, study, run, focus, ui);
    }
    egui::CollapsingHeader::new(tr("multiphysics.manifest")).show(ui, |ui| {
        ui.monospace(pretty(&run.manifest).chars().take(16384).collect::<String>());
    });
    ui.horizontal(|ui| {
        ui.label(tr("multiphysics.export_directory"));
        crate::ui_kit::text_field(ui, &mut state.export_directory, "");
        if button(
            ui,
            "multiphysics.export_result",
            "multiphysics:export_result",
            !state.export_directory.trim().is_empty() && !state.export_busy() && !cfg!(target_arch = "wasm32"),
        )
        .backing_on_hover_text(tr(if cfg!(target_arch = "wasm32") {
            "prefs.appearance.desktop_only"
        } else if state.export_busy() {
            "multiphysics.export_busy"
        } else {
            "multiphysics.artifact_hint"
        }))
        .clicked()
        {
            let directory = state.export_directory.clone();
            command(app, state, "multiphysics.result.export", json!({"id":study.id,"run_id":run.id,"directory":directory}));
        }
    });
    if all {
        results_artifacts(app, state, study, run, focus, ui);
    }
}

/// The run's artifacts as the focused node needs them: every artifact for the Results folder,
/// field files for cut plots, parameters and animations, tables for XY plots.
fn results_artifacts(
    app: &mut SolveApp,
    state: &mut State,
    study: &Study,
    run: &solidcraft_engine::multiphysics::RunEvidence,
    focus: shell::Node,
    ui: &mut egui::Ui,
) {
    use shell::Node;
    let all = focus == Node::Results;
    let wanted = |name: &str| match focus {
        Node::CutPlots | Node::Animation | Node::PointParams | Node::VolumeParams | Node::Compare => name.ends_with(".vtu"),
        Node::XyPlots => name.ends_with(".csv"),
        Node::Results => true,
        _ => false,
    };
    if !run.artifacts.iter().any(|artifact| wanted(&artifact.name)) {
        if !matches!(focus, Node::GoalPlots | Node::Report | Node::Compare) {
            ui.label(tr("multiphysics.no_artifacts"));
        }
        return;
    }
    if !run.artifacts.get(state.selected_artifact).is_some_and(|artifact| wanted(&artifact.name))
        && let Some(index) = run.artifacts.iter().position(|artifact| wanted(&artifact.name))
    {
        state.selected_artifact = index;
    }
    if run.artifacts.is_empty() {
        ui.label(tr("multiphysics.no_artifacts"));
        return;
    }
    state.selected_artifact = state.selected_artifact.min(run.artifacts.len().saturating_sub(1));
    let Some(artifact) = run.artifacts.get(state.selected_artifact) else { return };
    crate::ui_kit::dropdown_menu(ui, "multiphysics_artifact", &artifact.name, |ui| {
        for (index, artifact) in run.artifacts.iter().enumerate().filter(|(_, artifact)| wanted(&artifact.name)) {
            if crate::ui_kit::list_picker(ui, &format!("multiphysics_artifact_{index}"), &artifact.name, state.selected_artifact == index).clicked() {
                state.selected_artifact = index;
            }
        }
    });
    let frames = if all || matches!(focus, Node::CutPlots | Node::Animation) { pvd_frames(run) } else { Vec::new() };
    if frame_stepper(ui, &mut state.view, &frames) {
        ui.ctx().request_repaint();
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(next) = run.artifacts.get(state.selected_artifact)
            && next.name.ends_with(".vtu")
        {
            // Keep the field view on screen while stepping through frames.
            start_inspection(app, state, &study.id, &run.id, &next.name);
        }
    }
    if artifact.name.ends_with(".vtu") || artifact.name.ends_with(".pvd") {
        #[cfg(not(target_arch = "wasm32"))]
        {
            if button(ui, "cmd.multiphysics.result.inspect", "multiphysics:inspect_result", !inspection_busy(state)).clicked() {
                start_inspection(app, state, &study.id, &run.id, &artifact.name);
            }
            if inspection_busy(state) {
                ui.spinner();
            }
            if let Some(cached) = &state.inspected
                && cached.key.matches(app, &study.id, &run.id, &artifact.name)
            {
                match &cached.outcome {
                    Ok(value) => {
                        if let Some(coordinate) = value.frame_coordinate {
                            ui.monospace(format!("{}: {coordinate:.6e}", tr("multiphysics.slice_frame_coordinate")));
                        }
                        if let Some(spatial) = &value.spatial {
                            match spatial {
                                Ok(grid) => slice_view(app, ui, &mut state.view, grid),
                                Err(error) => {
                                    ui.label(tr("multiphysics.slice_unsupported"));
                                    ui.colored_label(Tokens::get().warning, error);
                                }
                            }
                        }
                        egui::CollapsingHeader::new(tr("multiphysics.slice_metadata")).show(ui, |ui| {
                            ui.monospace(bounded_metadata(&value.metadata, 16384));
                        });
                        ui.label(tr("multiphysics.preview_bounded"));
                    }
                    Err(error) => {
                        ui.colored_label(Tokens::get().form_error, error.chars().take(4096).collect::<String>());
                    }
                }
            }
        }
        #[cfg(target_arch = "wasm32")]
        ui.label(tr("prefs.appearance.desktop_only"));
    }
    if artifact.name.ends_with(".vtu") && (all || matches!(focus, Node::PointParams | Node::VolumeParams | Node::Compare)) {
        let fields = inspected_fields(app, state, &study.id, &run.id, &artifact.name);
        field_tools(app, state, study, run, &artifact.name, &fields, ui);
    }
    if artifact.name.ends_with(".csv") {
        let table = NumericTable::parse(&artifact.text);
        table_view(state, &table, ui);
    } else if all {
        egui::CollapsingHeader::new(tr("multiphysics.source_artifact")).id_salt("multiphysics_artifact_source").default_open(false).show(ui, |ui| {
            egui::ScrollArea::both().id_salt("multiphysics_artifact_text").max_height(180.0).show(ui, |ui| {
                ui.monospace(artifact.text.chars().take(32768).collect::<String>());
                ui.label(tr("multiphysics.preview_bounded"));
            });
        });
    }
}

/// A fixed orthographic data chart; it does not install another geometry camera.
#[cfg(not(target_arch = "wasm32"))]
fn slice_view(app: &mut SolveApp, ui: &mut egui::Ui, view: &mut ViewState, grid: &solidcraft_engine::multiphysics::vtk::BoxGrid) {
    ui.heading(tr("multiphysics.slice_title"));
    let previous_axis = view.slice_axis;
    ui.horizontal_wrapped(|ui| {
        for (enabled, label) in [(false, "multiphysics.field_2d"), (true, "multiphysics.field_3d")] {
            let response = crate::ui_kit::choice_chip(ui, label, tr(label), view.show_3d == enabled);
            if response.clicked() {
                view.show_3d = enabled;
            }
            crate::scenario::publish_target(if enabled { "multiphysics:view_3d" } else { "multiphysics:view_2d" }, response.rect);
        }
        ui.label(tr("multiphysics.slice_plane"));
        for (axis, label) in [(2, "XY"), (1, "XZ"), (0, "YZ")] {
            let response = crate::ui_kit::choice_chip(ui, &format!("multiphysics_plane_{axis}"), label, view.slice_axis == axis);
            if response.clicked() {
                view.slice_axis = axis;
            }
            crate::scenario::publish_target(&format!("multiphysics:slice_plane:{axis}"), response.rect);
        }
        if previous_axis != view.slice_axis {
            view.slice_layer = 0;
        }
        let layer_count = grid.planes.get(view.slice_axis).map_or(0, Vec::len);
        view.slice_layer = view.slice_layer.min(layer_count.saturating_sub(1));
        ui.label(tr("multiphysics.slice_layer"));
        let response = crate::ui_kit::number_field(ui, &mut view.slice_layer, 0..=layer_count.saturating_sub(1), 1.0, "");
        crate::scenario::publish_target("multiphysics:slice_layer", response.rect);
    });
    let Some(planes) = grid.planes.get(view.slice_axis) else { return };
    view.slice_field = view.slice_field.min(grid.fields.len().saturating_sub(1));
    let Some(field) = grid.fields.get(view.slice_field) else { return };
    let field_menu = ui
        .horizontal(|ui| {
            ui.label(tr(if field.components == 3 { "multiphysics.slice_magnitude" } else { "multiphysics.slice_scalar" }));
            crate::ui_kit::dropdown_menu(ui, "multiphysics_slice_field", &field.name, |ui| {
                for (index, field) in grid.fields.iter().enumerate() {
                    let response =
                        crate::ui_kit::list_picker(ui, &format!("multiphysics_slice_field_{index}"), &field.name, view.slice_field == index);
                    if response.clicked() {
                        view.slice_field = index;
                    }
                    crate::scenario::publish_target(&format!("multiphysics:slice_field:{index}"), response.rect);
                }
            })
        })
        .inner;
    crate::scenario::publish_target("multiphysics:slice_field_menu", field_menu.rect);
    field_menu.backing_on_hover_text(tr("multiphysics.slice_scope"));
    let Some(&coordinate) = planes.get(view.slice_layer) else { return };
    let (horizontal, vertical) = match view.slice_axis {
        0 => (1, 2),
        1 => (0, 2),
        2 => (0, 1),
        _ => return,
    };
    let Some(&[xmin, xmax]) = grid.bounds.get(horizontal) else { return };
    let Some(&[ymin, ymax]) = grid.bounds.get(vertical) else { return };
    ui.monospace(format!(
        "{}: {xmin:.3e} … {xmax:.3e}   {}: {ymin:.3e} … {ymax:.3e}   {}={coordinate:.3e}",
        ["X", "Y", "Z"].get(horizontal).copied().unwrap_or("?"),
        ["X", "Y", "Z"].get(vertical).copied().unwrap_or("?"),
        ["X", "Y", "Z"].get(view.slice_axis).copied().unwrap_or("?")
    ));
    let dx = xmax - xmin;
    let dy = ymax - ymin;
    if dx <= 0.0 || dy <= 0.0 {
        return;
    }
    let extrema = (0..grid.cells.len()).filter_map(|cell| grid.value(view.slice_field, cell)).fold(None, |range, value| {
        Some(match range {
            None => (value, value),
            Some((low, high)) => (value.min(low), value.max(high)),
        })
    });
    let Some(data_range) = extrema else { return };
    let (low, high) = color_controls(ui, view, data_range);
    let map = view.color_map;
    let unit = post::unit_hint(&field.name);
    let t = Tokens::get();
    let visible_width = (ui.clip_rect().right() - ui.next_widget_position().x).max(1.0);
    let plot_height = (ui.clip_rect().bottom() - ui.next_widget_position().y - 92.0).clamp(100.0, 600.0);
    if view.show_3d {
        view.field_3d.color_map = map;
        view.field_3d.color_range = Some((low, high));
        match crate::multiphysics_field_view::show(app, ui, &mut view.field_3d, grid, view.slice_field, view.slice_axis, coordinate, plot_height) {
            Ok(visible) => field_legend(ui, map, low, high, unit, visible),
            Err(error) => {
                ui.colored_label(t.form_error, error);
            }
        }
        return;
    }
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width().min(visible_width).max(1.0), plot_height), Sense::hover());
    crate::scenario::publish_rect("multiphysics:field_slice", rect);
    ui.painter().rect_filled(rect, t.control_radius, t.control_fill);
    let available = rect.shrink(t.space_12);
    let aspect = dx / dy;
    let width = available.width().min(available.height() * aspect as f32);
    let height = available.height().min(available.width() / aspect as f32);
    if !width.is_finite() || !height.is_finite() || width <= 0.0 || height <= 0.0 {
        return;
    }
    let area = egui::Rect::from_center_size(available.center(), egui::vec2(width, height));
    crate::scenario::publish_rect("multiphysics:field_data", area);
    let position = |x: f64, y: f64| egui::pos2(area.left() + ((x - xmin) / dx) as f32 * width, area.bottom() - ((y - ymin) / dy) as f32 * height);
    let mut visible = 0;
    for (cell_index, cell) in grid.cells.iter().enumerate() {
        if !grid.intersects(cell_index, view.slice_axis, coordinate) {
            continue;
        }
        let Some(value) = grid.value(view.slice_field, cell_index) else { continue };
        let Some(x) = cell.get(horizontal) else { continue };
        let Some(y) = cell.get(vertical) else { continue };
        let color = map_color(map, post::ColorRange { min: low, max: high, locked: true }.normalize(value));
        ui.painter().rect_filled(egui::Rect::from_two_pos(position(x[0], y[0]), position(x[1], y[1])), 0.0, color);
        visible += 1;
    }
    field_legend(ui, map, low, high, unit, visible);
}

#[cfg(not(target_arch = "wasm32"))]
fn tick_label(value: f64) -> String {
    if value != 0.0 && !(1e-3..1e4).contains(&value.abs()) { format!("{value:.2e}") } else { format!("{value}") }
}

/// A continuous gradient bar with nice tick values, the range ends and the unit.
#[cfg(not(target_arch = "wasm32"))]
fn field_legend(ui: &mut egui::Ui, map: post::ColorMap, low: f64, high: f64, unit: Option<&str>, visible: usize) {
    let t = Tokens::get();
    let legend = ui.vertical(|ui| {
        ui.horizontal_wrapped(|ui| {
            ui.label(tr("multiphysics.slice_value_range")).backing_on_hover_text(tr("multiphysics.slice_units"));
            ui.monospace(format!("{low:.4e} … {high:.4e}"));
            if let Some(unit) = unit {
                ui.monospace(format!("[{unit}]"));
            }
            ui.monospace(format!("{}: {visible}", tr("multiphysics.slice_cells")));
        });
        let width = ui.available_width().clamp(120.0, 480.0);
        let (rect, _) = ui.allocate_exact_size(egui::vec2(width, 14.0 + t.font_small + 4.0), Sense::hover());
        let bar = egui::Rect::from_min_size(rect.min, egui::vec2(width, 12.0));
        let steps = 64;
        for step in 0..steps {
            let left = bar.left() + bar.width() * step as f32 / steps as f32;
            let right = bar.left() + bar.width() * (step + 1) as f32 / steps as f32;
            let color = map_color(map, (f64::from(step) + 0.5) / f64::from(steps));
            ui.painter().rect_filled(egui::Rect::from_min_max(egui::pos2(left, bar.top()), egui::pos2(right, bar.bottom())), 0.0, color);
        }
        let span = high - low;
        if span > 0.0 && span.is_finite() {
            for tick in post::legend_ticks(low, high, 5) {
                let fraction = ((tick - low) / span) as f32;
                if !(0.0..=1.0).contains(&fraction) {
                    continue;
                }
                let x = bar.left() + bar.width() * fraction;
                ui.painter().line_segment([egui::pos2(x, bar.bottom()), egui::pos2(x, bar.bottom() + 3.0)], egui::Stroke::new(1.0, t.text_dim));
                ui.painter().text(
                    egui::pos2(x, bar.bottom() + 3.0),
                    egui::Align2::CENTER_TOP,
                    tick_label(tick),
                    egui::FontId::proportional(t.font_small),
                    t.text_dim,
                );
            }
        }
    });
    crate::scenario::publish_rect("multiphysics:field_legend", legend.response.rect);
}

/// Limit serialization work as well as the number of displayed bytes.
#[cfg(not(target_arch = "wasm32"))]
fn bounded_metadata(value: &Value, limit: usize) -> String {
    struct Bounded {
        bytes: Vec<u8>,
        limit: usize,
    }
    impl std::io::Write for Bounded {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            let remaining = self.limit.saturating_sub(self.bytes.len());
            if remaining == 0 {
                return Err(std::io::Error::other("metadata preview limit"));
            }
            let count = remaining.min(bytes.len());
            self.bytes.extend_from_slice(&bytes[..count]);
            Ok(count)
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut output = Bounded { bytes: Vec::new(), limit };
    let _ = serde_json::to_writer_pretty(&mut output, value);
    String::from_utf8_lossy(&output.bytes).into_owned()
}

struct NumericTable {
    headers: Vec<String>,
    rows: Vec<Vec<Option<f64>>>,
}
impl NumericTable {
    fn parse(text: &str) -> Self {
        let mut lines = text.lines().filter(|line| !line.trim().is_empty() && !line.starts_with('#'));
        let headers = lines.next().unwrap_or_default().split(',').take(MAX_TABLE_COLUMNS).map(str::to_string).collect::<Vec<_>>();
        let rows = lines
            .take(MAX_TABLE_ROWS)
            .map(|line| {
                line.split(',').take(headers.len()).map(|value| value.trim().parse::<f64>().ok().filter(|number| number.is_finite())).collect()
            })
            .collect();
        Self { headers, rows }
    }
}

fn table_view(state: &mut State, table: &NumericTable, ui: &mut egui::Ui) {
    if table.headers.is_empty() {
        ui.label(tr("multiphysics.no_numeric_data"));
        return;
    }
    state.x_column = state.x_column.min(table.headers.len().saturating_sub(1));
    state.y_column = state.y_column.min(table.headers.len().saturating_sub(1));
    ui.columns(2, |columns| {
        let ui = &mut columns[0];
        ui.label(tr("multiphysics.x_axis"));
        let x = table.headers.get(state.x_column).cloned().unwrap_or_default();
        let width = ui.available_width();
        crate::ui_kit::bounded_dropdown_menu(ui, "multiphysics_x", &x, width, |ui| {
            for (index, name) in table.headers.iter().enumerate() {
                if crate::ui_kit::list_picker(ui, &format!("multiphysics_x_{index}"), name, state.x_column == index).clicked() {
                    state.x_column = index;
                }
            }
        });
        crate::scenario::publish_rect("multiphysics:x_axis", ui.min_rect());
        let ui = &mut columns[1];
        ui.label(tr("multiphysics.y_axis"));
        let y = table.headers.get(state.y_column).cloned().unwrap_or_default();
        let width = ui.available_width();
        crate::ui_kit::bounded_dropdown_menu(ui, "multiphysics_y", &y, width, |ui| {
            for (index, name) in table.headers.iter().enumerate() {
                if crate::ui_kit::list_picker(ui, &format!("multiphysics_y_{index}"), name, state.y_column == index).clicked() {
                    state.y_column = index;
                }
            }
        });
        crate::scenario::publish_rect("multiphysics:y_axis", ui.min_rect());
    });
    let points: Vec<_> =
        table.rows.iter().map(|row| Some((row.get(state.x_column).copied().flatten()?, row.get(state.y_column).copied().flatten()?))).collect();
    plot(ui, &points);
    ui.label(RichText::new(tr("multiphysics.numeric_preview")).small().color(Tokens::get().text_dim));
    let row_height = ui.text_style_height(&egui::TextStyle::Monospace) + ui.spacing().item_spacing.y;
    egui::ScrollArea::both().id_salt("multiphysics_table").max_height(220.0).show_rows(ui, row_height, table.rows.len() + 1, |ui, range| {
        egui::Grid::new("multiphysics_data").striped(true).show(ui, |ui| {
            for index in range {
                if index == 0 {
                    for header in &table.headers {
                        ui.strong(header);
                    }
                } else if let Some(row) = table.rows.get(index - 1) {
                    for column in 0..table.headers.len() {
                        if let Some(number) = row.get(column).copied().flatten() {
                            ui.monospace(format!("{number:.6e}"));
                        } else {
                            ui.label(tr("multiphysics.missing_sample"));
                        }
                    }
                }
                ui.end_row();
            }
        });
    });
}

/// A bounded data chart, not a second geometry camera. Gaps remain gaps.
fn plot(ui: &mut egui::Ui, points: &[Option<(f64, f64)>]) {
    plot_series(ui, &[(Tokens::get().accent, points.to_vec())]);
}

/// Several series on one shared x/y scale (callers normalise series of different units).
fn plot_series(ui: &mut egui::Ui, series: &[(egui::Color32, Vec<Option<(f64, f64)>>)]) {
    let t = Tokens::get();
    let visible_width = (ui.clip_rect().right() - ui.next_widget_position().x).max(1.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width().min(visible_width).max(1.0), 220.0), Sense::hover());
    crate::scenario::publish_rect("multiphysics:plot", rect);
    let bounds = series.iter().flat_map(|(_, points)| points.iter().flatten()).fold(None, |bounds, &(x, y)| {
        Some(match bounds {
            None => (x, x, y, y),
            Some((xmin, xmax, ymin, ymax)) => (x.min(xmin), x.max(xmax), y.min(ymin), y.max(ymax)),
        })
    });
    ui.painter().rect_filled(rect, t.control_radius, t.control_fill);
    let Some((xmin, xmax, ymin, ymax)) = bounds else { return };
    let area = rect.shrink2(egui::vec2(t.space_12, t.space_12));
    let dx = if xmax == xmin { 1.0 } else { xmax - xmin };
    let dy = if ymax == ymin { 1.0 } else { ymax - ymin };
    if !dx.is_finite() || !dy.is_finite() {
        return;
    }
    let position = |(x, y): (f64, f64)| {
        egui::pos2(area.left() + ((x - xmin) / dx) as f32 * area.width(), area.bottom() - ((y - ymin) / dy) as f32 * area.height())
    };
    // Recessive horizontal grid at nice values.
    for tick in post::legend_ticks(ymin, ymax, 4) {
        let y = position((xmin, tick)).y;
        if (area.top()..=area.bottom()).contains(&y) {
            ui.painter()
                .line_segment([egui::pos2(area.left(), y), egui::pos2(area.right(), y)], egui::Stroke::new(0.5, t.text_dim.gamma_multiply(0.35)));
        }
    }
    if let Some(point) = series.first().and_then(|(_, points)| points.last().copied().flatten()) {
        crate::scenario::publish_rect("multiphysics:plot_last_sample", egui::Rect::from_center_size(position(point), egui::Vec2::ZERO));
    }
    for (color, points) in series {
        let mut previous = None;
        for point in points {
            let current = point.map(position);
            if let (Some(a), Some(b)) = (previous, current) {
                ui.painter().line_segment([a, b], egui::Stroke::new(1.5, *color));
            }
            previous = current;
        }
    }
    ui.painter().text(area.left_top(), egui::Align2::LEFT_TOP, format!("{ymax:.4e}"), egui::FontId::proportional(t.font_small), t.text_dim);
    ui.painter().text(area.left_bottom(), egui::Align2::LEFT_BOTTOM, format!("{xmin:.4e}"), egui::FontId::proportional(t.font_small), t.text_dim);
    ui.painter().text(area.right_bottom(), egui::Align2::RIGHT_BOTTOM, format!("{xmax:.4e}"), egui::FontId::proportional(t.font_small), t.text_dim);
}

/// Field names of the inspected VTU, when the shown artifact has been inspected.
fn inspected_fields(app: &SolveApp, state: &State, study: &str, run: &str, artifact: &str) -> Vec<String> {
    #[cfg(not(target_arch = "wasm32"))]
    if let Some(cached) = &state.inspected
        && cached.key.matches(app, study, run, artifact)
        && let Ok(output) = &cached.outcome
        && let Some(Ok(grid)) = &output.spatial
    {
        return grid.fields.iter().map(|field| field.name.clone()).collect();
    }
    let _ = (app, state, study, run, artifact);
    Vec::new()
}

/// Categorical series colours (Okabe–Ito order: distinguishable with common colour-vision
/// deficiencies). Identity is also given by the goal table beside each swatch.
const SERIES_COLORS: [egui::Color32; 8] = [
    egui::Color32::from_rgb(0, 114, 178),
    egui::Color32::from_rgb(230, 159, 0),
    egui::Color32::from_rgb(0, 158, 115),
    egui::Color32::from_rgb(213, 94, 0),
    egui::Color32::from_rgb(86, 180, 233),
    egui::Color32::from_rgb(204, 121, 167),
    egui::Color32::from_rgb(240, 228, 66),
    egui::Color32::from_rgb(153, 153, 153),
];

#[cfg(not(target_arch = "wasm32"))]
fn map_color(map: post::ColorMap, t: f64) -> egui::Color32 {
    let [r, g, b] = map.color(t);
    egui::Color32::from_rgb(r, g, b)
}

fn artifact_text<'a>(run: &'a solidcraft_engine::multiphysics::RunEvidence, name: &str) -> Option<&'a str> {
    run.artifacts.iter().find(|artifact| artifact.name == name).map(|artifact| artifact.text.as_str())
}

fn post_cache<'a>(
    slot: &'a mut Option<PostCache>,
    project: &std::sync::Arc<solidcraft_engine::multiphysics::Project>,
    study: &str,
    run: &str,
) -> &'a mut PostCache {
    if !slot.as_ref().is_some_and(|cache| std::sync::Arc::ptr_eq(&cache.project, project) && cache.study == study && cache.run == run) {
        *slot = None;
    }
    slot.get_or_insert_with(|| PostCache {
        project: project.clone(),
        study: study.into(),
        run: run.into(),
        goals: None,
        components: None,
        comparison: None,
        field: None,
    })
}

/// Goals, component compliance and study comparison for the shown run (FloEFD-style goal
/// plots and results summary). Everything is read from retained evidence; nothing re-solves.
fn post_processing(
    app: &SolveApp,
    state: &mut State,
    study: &Study,
    run: &solidcraft_engine::multiphysics::RunEvidence,
    focus: shell::Node,
    ui: &mut egui::Ui,
) {
    let all = focus == shell::Node::Results;
    let project = app.session.doc.multiphysics.clone();
    let cache = post_cache(&mut state.post, &project, &study.id, &run.id);
    if (all || focus == shell::Node::GoalPlots)
        && let Some(csv) = artifact_text(run, "goals.csv")
    {
        let report = cache.goals.get_or_insert_with(|| post::goals(csv, artifact_text(run, "goals.json")).map_err(|error| error.to_string()));
        egui::CollapsingHeader::new(tr("multiphysics.post_goals")).id_salt("multiphysics_post_goals").default_open(true).show(ui, |ui| {
            goals_view(ui, report, &mut state.view.goal_index);
        });
    }
    if (all || focus == shell::Node::Report)
        && let Some(text) = artifact_text(run, "thermal_report.json")
    {
        let rows = cache.components.get_or_insert_with(|| post::components_from_report(text).map_err(|error| error.to_string()));
        egui::CollapsingHeader::new(tr("multiphysics.post_components")).id_salt("multiphysics_post_components").default_open(true).show(ui, |ui| {
            components_view(ui, rows);
        });
    }
    if !(all || focus == shell::Node::Compare) {
        return;
    }
    egui::CollapsingHeader::new(tr("multiphysics.post_compare_studies"))
        .id_salt("multiphysics_post_compare")
        .default_open(focus == shell::Node::Compare)
        .show(ui, |ui| {
            let response = crate::ui_kit::button(ui, tr("cmd.multiphysics.study.compare"), crate::ui_kit::Button::Normal);
            crate::scenario::publish_target("multiphysics:compare_studies", response.rect);
            if response.clicked() {
                cache.comparison = Some(post::compare_studies(&project, &[]).map_err(|error| error.to_string()));
            }
            match &cache.comparison {
                Some(Ok(table)) => metric_table_view(ui, table),
                Some(Err(error)) => {
                    ui.colored_label(Tokens::get().form_error, error);
                }
                None => {
                    ui.label(RichText::new(tr("multiphysics.post_compare_hint")).small().color(Tokens::get().text_dim));
                }
            }
        });
}

fn goals_view(ui: &mut egui::Ui, report: &Result<post::GoalReport, String>, selected: &mut Option<usize>) {
    let t = Tokens::get();
    let report = match report {
        Ok(report) => report,
        Err(error) => {
            ui.colored_label(t.form_error, error);
            return;
        }
    };
    match report.all_converged {
        Some(true) => ui.label(format!("✓ {}", tr("multiphysics.post_goals_converged"))),
        Some(false) => ui.colored_label(t.warning, format!("⚠ {}", tr("multiphysics.post_goals_not_converged"))),
        None => ui.label(tr("multiphysics.result_unreported")),
    };
    if let Some(step) = report.stopped_early_at_step {
        ui.label(solidcraft_ui_text::format("multiphysics.post_goals_stopped", &[("step", &format!("{step:.0}"))]));
    }
    if let Some(tolerance) = report.relative_tolerance {
        ui.monospace(format!("{}: {tolerance:.1e}", tr("multiphysics.post_goals_tolerance")));
    }
    egui::Grid::new("multiphysics_goals").striped(true).show(ui, |ui| {
        ui.label("");
        for header in ["multiphysics.post_goal", "multiphysics.post_goal_final", "multiphysics.post_goal_change", "multiphysics.post_goal_status"] {
            ui.strong(tr(header));
        }
        ui.end_row();
        for (index, goal) in report.goals.iter().enumerate() {
            let color = SERIES_COLORS[index % SERIES_COLORS.len()];
            let shown = selected.is_none_or(|selected| selected == index);
            // The swatch ties the row to its line in the chart; the name selects it.
            let (swatch, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), Sense::hover());
            ui.painter().rect_filled(swatch, 2.0, color);
            let unit = goal.unit.as_deref().or_else(|| post::unit_hint(&goal.name)).unwrap_or("");
            let name = if unit.is_empty() { goal.name.clone() } else { format!("{} [{unit}]", goal.name) };
            if crate::ui_kit::choice_chip(ui, &format!("multiphysics_goal_{index}"), &name, shown && selected.is_some()).clicked() {
                *selected = if *selected == Some(index) { None } else { Some(index) };
            }
            ui.monospace(goal.final_value.map_or_else(|| tr("multiphysics.missing_sample").to_string(), |value| format!("{value:.6e}")));
            ui.monospace(goal.last_relative_change.map_or_else(String::new, |value| format!("{value:.2e}")));
            match goal.converged {
                Some(true) => ui.label("✓"),
                Some(false) => ui.colored_label(t.warning, "⚠"),
                None => ui.label("–"),
            };
            ui.end_row();
        }
    });
    // Each goal is scaled to its own range so goals of different units share one chart
    // without a second axis; the table above carries the values.
    let series: Vec<(egui::Color32, Vec<Option<(f64, f64)>>)> = report
        .goals
        .iter()
        .enumerate()
        .filter(|(index, _)| selected.is_none_or(|selected| selected == *index))
        .map(|(index, goal)| {
            let range = post::ColorRange::from_values(goal.values.iter().copied());
            let points = goal
                .steps
                .iter()
                .zip(&goal.values)
                .map(|(&step, &value)| {
                    let y = if selected.is_some() { value } else { range.map_or(0.5, |range| range.normalize(value)) };
                    (step.is_finite() && y.is_finite()).then_some((step, y))
                })
                .collect();
            (SERIES_COLORS[index % SERIES_COLORS.len()], points)
        })
        .collect();
    plot_series(ui, &series);
    if selected.is_none() && report.goals.len() > 1 {
        ui.label(RichText::new(tr("multiphysics.post_goals_normalized")).small().color(t.text_dim));
    }
}

fn components_view(ui: &mut egui::Ui, rows: &Result<Vec<post::ComponentRow>, String>) {
    let t = Tokens::get();
    let rows = match rows {
        Ok(rows) => rows,
        Err(error) => {
            ui.colored_label(t.form_error, error);
            return;
        }
    };
    egui::ScrollArea::horizontal().id_salt("multiphysics_components_scroll").show(ui, |ui| {
        egui::Grid::new("multiphysics_components").striped(true).show(ui, |ui| {
            for header in [
                "multiphysics.post_component",
                "multiphysics.post_component_power",
                "multiphysics.post_component_mean",
                "multiphysics.post_component_max",
                "multiphysics.post_component_limit",
                "multiphysics.post_component_margin",
                "multiphysics.post_goal_status",
            ] {
                ui.strong(tr(header));
            }
            ui.end_row();
            for row in rows {
                ui.monospace(&row.name);
                ui.monospace(format!("{:.3}", row.power_w));
                ui.monospace(format!("{:.2}", row.mean_temperature_k));
                ui.monospace(format!("{:.2}", row.max_temperature_k));
                ui.monospace(row.limit_k.map_or_else(String::new, |value| format!("{value:.2}")));
                ui.monospace(row.margin_k.map_or_else(String::new, |value| format!("{value:+.2}")));
                match row.status {
                    post::ComponentStatus::Pass => ui.label(format!("✓ {}", tr("multiphysics.post_status_pass"))),
                    post::ComponentStatus::Warn => ui.colored_label(t.warning, format!("⚠ {}", tr("multiphysics.post_status_warn"))),
                    post::ComponentStatus::Fail => ui.colored_label(t.form_error, format!("✗ {}", tr("multiphysics.post_status_fail"))),
                    post::ComponentStatus::NoLimit => ui.colored_label(t.text_dim, format!("– {}", tr("multiphysics.post_status_no_limit"))),
                };
                ui.end_row();
            }
        });
    });
    ui.label(RichText::new(tr("multiphysics.post_component_units")).small().color(t.text_dim));
}

fn metric_table_view(ui: &mut egui::Ui, table: &post::MetricTable) {
    egui::ScrollArea::both().id_salt("multiphysics_compare_scroll").max_height(260.0).show(ui, |ui| {
        egui::Grid::new("multiphysics_compare").striped(true).show(ui, |ui| {
            ui.strong(tr("multiphysics.post_metric"));
            for run in &table.runs {
                ui.strong(run);
            }
            ui.end_row();
            for (metric, values) in table.metrics.iter().zip(&table.values).take(MAX_TABLE_ROWS) {
                ui.monospace(metric);
                for value in values {
                    ui.monospace(value.map_or_else(|| tr("multiphysics.missing_sample").to_string(), |value| format!("{value:.6e}")));
                }
                ui.end_row();
            }
        });
    });
}

/// Frames of the run's PVD collection, as `(artifact index, coordinate)` in file order.
fn pvd_frames(run: &solidcraft_engine::multiphysics::RunEvidence) -> Vec<(usize, f64)> {
    let mut frames = Vec::new();
    for artifact in run.artifacts.iter().filter(|artifact| artifact.name.ends_with(".pvd")) {
        let Ok(list) = solidcraft_engine::multiphysics::vtk::parse_pvd(&artifact.text) else { continue };
        for frame in list {
            if let Some(index) = run.artifacts.iter().position(|artifact| artifact.name == frame.file)
                && !frames.iter().any(|&(known, _)| known == index)
            {
                frames.push((index, frame.coordinate));
            }
        }
    }
    frames
}

/// Step through a transient or iteration series without opening the artifact list.
fn frame_stepper(ui: &mut egui::Ui, view: &mut ViewState, frames: &[(usize, f64)]) -> bool {
    if frames.len() < 2 {
        return false;
    }
    let mut position = frames.iter().position(|&(index, _)| index == view.selected_artifact).unwrap_or(0);
    let before = position;
    ui.horizontal(|ui| {
        ui.label(tr("multiphysics.post_frame"));
        if crate::ui_kit::button(ui, "◀", crate::ui_kit::Button::Small).clicked() {
            position = position.saturating_sub(1);
        }
        let response = crate::ui_kit::slider(
            ui,
            &mut position,
            0..=frames.len() - 1,
            crate::ui_kit::Slider { integer: true, show_value: false, ..Default::default() },
        );
        crate::scenario::publish_target("multiphysics:frame_slider", response.rect);
        if crate::ui_kit::button(ui, "▶", crate::ui_kit::Button::Small).clicked() {
            position = (position + 1).min(frames.len() - 1);
        }
        if let Some(&(_, coordinate)) = frames.get(position) {
            ui.monospace(format!("{}/{}  {coordinate:.6e}", position + 1, frames.len()));
        }
    });
    if position != before
        && let Some(&(index, _)) = frames.get(position)
    {
        view.selected_artifact = index;
        return true;
    }
    false
}

/// Probe, extrema, line sample, volume integral and run-to-run difference on a retained VTU.
fn field_tools(
    app: &SolveApp,
    state: &mut State,
    study: &Study,
    run: &solidcraft_engine::multiphysics::RunEvidence,
    artifact: &str,
    fields: &[String],
    ui: &mut egui::Ui,
) {
    let project = app.session.doc.multiphysics.clone();
    let view = &mut state.view;
    let cache = post_cache(&mut state.post, &project, &study.id, &run.id);
    egui::CollapsingHeader::new(tr("multiphysics.post_field_tools")).id_salt("multiphysics_field_tools").default_open(true).show(ui, |ui| {
        if view.probe_field.is_empty()
            && let Some(first) = fields.first()
        {
            view.probe_field = first.clone();
        }
        ui.horizontal_wrapped(|ui| {
            ui.label(tr("multiphysics.post_field"));
            crate::ui_kit::text_field(ui, &mut view.probe_field, "");
            if !fields.is_empty() {
                crate::ui_kit::dropdown_menu(ui, "multiphysics_probe_field", "▾", |ui| {
                    for name in fields {
                        if crate::ui_kit::list_picker(ui, &format!("multiphysics_probe_field_{name}"), name, &view.probe_field == name).clicked() {
                            view.probe_field = name.clone();
                        }
                    }
                });
            }
        });
        let coordinates = |ui: &mut egui::Ui, label: &str, point: &mut [f64; 3]| {
            ui.horizontal(|ui| {
                ui.label(tr(label));
                for (axis, value) in ["X", "Y", "Z"].into_iter().zip(point.iter_mut()) {
                    ui.label(axis);
                    crate::ui_kit::number_field(ui, value, f64::MIN..=f64::MAX, 1e-3, "");
                }
            });
        };
        coordinates(ui, "multiphysics.post_point", &mut view.probe_point);
        coordinates(ui, "multiphysics.post_line_end", &mut view.line_end);
        let field = view.probe_field.trim().to_string();
        let mut asked = None;
        ui.horizontal_wrapped(|ui| {
            let ready = !field.is_empty();
            for (key, target, query) in [
                (
                    "multiphysics.post_probe",
                    "multiphysics:probe",
                    post::Query::Probe { artifact: artifact.into(), field: field.clone(), point: view.probe_point },
                ),
                ("multiphysics.post_min_max", "multiphysics:min_max", post::Query::MinMax { artifact: artifact.into(), field: field.clone() }),
                (
                    "multiphysics.post_line",
                    "multiphysics:line",
                    post::Query::Line {
                        artifact: artifact.into(),
                        field: field.clone(),
                        start: view.probe_point,
                        end: view.line_end,
                        samples: view.line_samples,
                    },
                ),
                (
                    "multiphysics.post_integrate",
                    "multiphysics:integrate",
                    post::Query::Integrate { artifact: artifact.into(), field: field.clone(), region: None },
                ),
            ] {
                let response = crate::ui_kit::button(ui, tr(key), crate::ui_kit::Button::Enabled(ready));
                crate::scenario::publish_target(target, response.rect);
                if response.clicked() {
                    asked = Some(query);
                }
            }
            ui.label("n");
            crate::ui_kit::number_field(ui, &mut view.line_samples, 2..=10_000, 1.0, "");
        });
        let others: Vec<&str> = study
            .runs
            .iter()
            .filter(|other| other.id != run.id && other.artifacts.iter().any(|a| a.name == artifact))
            .map(|other| other.id.as_str())
            .collect();
        if !others.is_empty() {
            ui.horizontal_wrapped(|ui| {
                ui.label(tr("multiphysics.post_compare_run"));
                let first = others.first().map_or_else(String::new, |id| (*id).to_string());
                let shown = view.compare_run.clone().filter(|id| others.contains(&id.as_str())).unwrap_or(first);
                crate::ui_kit::dropdown_menu(ui, "multiphysics_compare_run", &shown, |ui| {
                    for id in &others {
                        if crate::ui_kit::list_picker(ui, &format!("multiphysics_compare_run_{id}"), id, view.compare_run.as_deref() == Some(*id))
                            .clicked()
                        {
                            view.compare_run = Some((*id).to_string());
                        }
                    }
                });
                let response = crate::ui_kit::button(ui, tr("multiphysics.post_difference"), crate::ui_kit::Button::Enabled(!field.is_empty()));
                crate::scenario::publish_target("multiphysics:difference", response.rect);
                if response.clicked() {
                    asked = Some(post::Query::Compare { artifact: artifact.into(), field: field.clone(), other_run_id: shown, other_study: None });
                }
            });
        }
        if let Some(query) = asked {
            let request = post::QueryRequest { id: study.id.clone(), run_id: run.id.clone(), query };
            cache.field = Some(post::query(&project, &request).map_err(|error| error.to_string()));
        }
        match &cache.field {
            Some(Ok(answer)) => field_answer(ui, answer),
            Some(Err(error)) => {
                ui.colored_label(Tokens::get().form_error, error);
            }
            None => {}
        }
    });
}

fn field_answer(ui: &mut egui::Ui, answer: &post::QueryResponse) {
    let t = Tokens::get();
    if answer.stale {
        ui.colored_label(t.warning, tr("multiphysics.stale_run"));
    }
    let unit = answer.unit.as_deref().map(|unit| format!(" {unit}")).unwrap_or_default();
    let point = |p: &[f64; 3]| format!("({:.4e}, {:.4e}, {:.4e})", p[0], p[1], p[2]);
    match &answer.result {
        post::QueryResult::Probe(probe) => {
            ui.monospace(format!("{} = {:.6e}{unit}  @ {}", tr("multiphysics.post_value"), probe.magnitude, point(&probe.point)));
            if probe.value.len() > 1 {
                ui.monospace(format!("[{}]", probe.value.iter().map(|v| format!("{v:.4e}")).collect::<Vec<_>>().join(", ")));
            }
            if !probe.inside {
                ui.colored_label(
                    t.warning,
                    solidcraft_ui_text::format("multiphysics.post_probe_outside", &[("distance", &format!("{:.3e}", probe.distance))]),
                );
            }
            ui.label(RichText::new(format!("{} {} · {}", tr("multiphysics.slice_cells"), probe.cell, probe.method)).small().color(t.text_dim));
        }
        post::QueryResult::MinMax(extrema) => {
            for (label, extremum) in [("multiphysics.post_minimum", &extrema.min), ("multiphysics.post_maximum", &extrema.max)] {
                ui.monospace(format!("{}: {:.6e}{unit}  @ {}", tr(label), extremum.value, point(&extremum.location)));
            }
        }
        post::QueryResult::Line { samples } => {
            let points: Vec<_> = samples.iter().map(|sample| sample.value.map(|value| (sample.s, value))).collect();
            plot(ui, &points);
            let missing = samples.iter().filter(|sample| sample.value.is_none()).count();
            if missing > 0 {
                ui.label(solidcraft_ui_text::format("multiphysics.post_line_gaps", &[("count", &missing.to_string())]));
            }
        }
        post::QueryResult::Integrate(stats) => {
            egui::Grid::new("multiphysics_integral").show(ui, |ui| {
                for (label, value) in [
                    ("multiphysics.post_integral", stats.integral),
                    ("multiphysics.post_mean", stats.mean),
                    ("multiphysics.post_minimum", stats.min),
                    ("multiphysics.post_maximum", stats.max),
                    ("multiphysics.post_measure", stats.measure),
                ] {
                    ui.label(tr(label));
                    ui.monospace(format!("{value:.6e}"));
                    ui.end_row();
                }
            });
            ui.label(RichText::new(format!("{} {}", stats.cells, tr("multiphysics.slice_cells"))).small().color(t.text_dim));
        }
        post::QueryResult::Compare(difference) => {
            egui::Grid::new("multiphysics_difference").show(ui, |ui| {
                ui.label(tr("multiphysics.post_max_difference"));
                ui.monospace(format!("{:.6e}{unit}  @ {}", difference.max_abs, point(&difference.max_abs_location)));
                ui.end_row();
                ui.label(tr("multiphysics.post_rms_difference"));
                ui.monospace(format!("{:.6e}{unit}", difference.rms));
                ui.end_row();
                ui.label(tr("multiphysics.post_relative_l2"));
                ui.monospace(if difference.relative_l2.is_finite() { format!("{:.3e}", difference.relative_l2) } else { "—".into() });
                ui.end_row();
            });
        }
        post::QueryResult::Goals(_) | post::QueryResult::Components { .. } => {}
    }
}

/// Colour map and an optional locked range shared by frames and views.
#[cfg(not(target_arch = "wasm32"))]
fn color_controls(ui: &mut egui::Ui, view: &mut ViewState, data: (f64, f64)) -> (f64, f64) {
    ui.horizontal_wrapped(|ui| {
        ui.label(tr("multiphysics.post_color_map"));
        crate::ui_kit::dropdown_menu(ui, "multiphysics_color_map", view.color_map.name(), |ui| {
            for map in post::ColorMap::ALL {
                if crate::ui_kit::list_picker(ui, &format!("multiphysics_color_map_{}", map.name()), map.name(), view.color_map == map).clicked() {
                    view.color_map = map;
                }
            }
        });
        let mut locked = view.locked_range.is_some();
        let response = crate::ui_kit::checkbox(ui, tr("multiphysics.post_lock_range"), &mut locked)
            .backing_on_hover_text(tr("multiphysics.post_lock_range_hint"));
        crate::scenario::publish_target("multiphysics:lock_range", response.rect);
        if locked != view.locked_range.is_some() {
            view.locked_range = locked.then_some(data);
        }
        if let Some((low, high)) = &mut view.locked_range {
            crate::ui_kit::number_field(ui, low, f64::MIN..=f64::MAX, (data.1 - data.0).abs().max(1e-12) * 1e-3, "");
            ui.label("…");
            crate::ui_kit::number_field(ui, high, f64::MIN..=f64::MAX, (data.1 - data.0).abs().max(1e-12) * 1e-3, "");
        }
    });
    view.locked_range.filter(|(low, high)| low.is_finite() && high.is_finite() && low < high).unwrap_or(data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scenario::Harness;

    #[test]
    fn dirty_cache_preserves_semantic_edits_and_invalid_drafts() {
        let study = solidcraft_engine::multiphysics::Project::default().create("poisson1d", "Cache").unwrap();
        let mut draft = Draft::from_study(&study);
        assert!(!draft.dirty());
        draft.case_text = serde_json::to_string(&study.case).unwrap();
        assert!(!draft.dirty(), "formatting alone is not an authored change");
        let cached = draft.case_dirty_cache.borrow().clone();
        assert!(!draft.dirty());
        assert_eq!(*draft.case_dirty_cache.borrow(), cached);
        draft.case_text = "{invalid".into();
        assert!(draft.dirty());
        assert!(draft.dirty());
        draft.case_text = "{}".into();
        assert!(draft.dirty());
        draft.case_text = pretty(&study.case);
        assert!(!draft.dirty());
        draft.question.push_str("new question");
        assert!(draft.dirty());
    }

    #[test]
    fn preflight_cache_follows_project_identity_not_revision_numbers() {
        use std::sync::Arc;
        let mut project = Arc::new(solidcraft_engine::multiphysics::Project::default());
        let study = Arc::make_mut(&mut project).create("poisson1d", "Cache").unwrap();
        let mut state = State::default();
        let original = state.study_preflight(&project, project.study(&study.id).unwrap());
        let repeat = state.study_preflight(&project, project.study(&study.id).unwrap());
        assert!(Arc::ptr_eq(&original, &repeat));
        assert!(original.errors.is_empty());
        // No revision increment: exercise replacement/undo identity rather than a revision key.
        let mut replacement = Arc::new(project.as_ref().clone());
        Arc::make_mut(&mut replacement).studies[0].case["solver"] = json!("unsupported");
        let changed = state.study_preflight(&replacement, replacement.study(&study.id).unwrap());
        assert!(!changed.errors.is_empty());
        assert!(!Arc::ptr_eq(&original, &changed));
        let restored = state.study_preflight(&project, project.study(&study.id).unwrap());
        assert!(restored.errors.is_empty());
        // A weak cache prevents in-place Arc mutation from retaining its old identity.
        Arc::make_mut(&mut project).studies[0].case["solver"] = json!("unsupported");
        let mutated = state.study_preflight(&project, project.study(&study.id).unwrap());
        assert!(!mutated.errors.is_empty());
        assert!(!Arc::ptr_eq(&restored, &mutated));
    }

    #[test]
    fn readiness_distinguishes_drafts_workers_and_historical_results() {
        fn readiness(study: &Study, draft: Option<&Draft>, worker: bool) -> Readiness {
            let valid = solidcraft_engine::multiphysics::preflight(study).errors.is_empty();
            Readiness::for_study(study, draft, worker, valid)
        }
        let mut study = solidcraft_engine::multiphysics::Project::default().create("poisson1d", "Readiness").unwrap();
        let mut draft = Draft::from_study(&study);
        assert!(!readiness(&study, Some(&draft), false).can_run());
        assert!(readiness(&study, Some(&draft), true).can_run());
        draft.question = "Unapplied question".into();
        assert!(!readiness(&study, Some(&draft), true).can_run());
        draft = Draft::from_study(&study);
        draft.expected_revision += 1;
        assert!(!readiness(&study, Some(&draft), true).can_run());
        let old_revision = study.revision;
        study.runs.push(solidcraft_engine::multiphysics::RunEvidence {
            id: "run-fixture".into(),
            study_revision: old_revision,
            case_hash: "fixture".into(),
            input_case: study.case.clone(),
            manifest: json!({}),
            artifacts: Vec::new(),
        });
        let current = readiness(&study, None, true);
        assert!(current.current_result && !current.historical_result);
        study.revision += 1;
        let historical = readiness(&study, None, true);
        assert!(!historical.current_result && historical.historical_result);
        study.case["solver"] = json!("unsupported-solver");
        assert!(!readiness(&study, None, true).can_run());
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn inspection_discards_changed_document_and_selection() {
        use std::sync::{Arc, atomic::AtomicBool, mpsc};
        let mut app = SolveApp::new(Default::default(), Default::default());
        let reply = app.run("multiphysics.study.create", json!({"template":"poisson1d","name":"Inspection"})).unwrap();
        let id = reply["study"]["id"].as_str().unwrap().to_string();
        let project = &mut std::sync::Arc::make_mut(&mut app.session.doc).multiphysics;
        let study = Arc::make_mut(project).studies.iter_mut().find(|study| study.id == id).unwrap();
        study.runs.push(solidcraft_engine::multiphysics::RunEvidence {
            id: "run-fixture".into(),
            study_revision: study.revision,
            case_hash: "fixture".into(),
            input_case: study.case.clone(),
            manifest: json!({}),
            artifacts: vec![solidcraft_engine::multiphysics::Artifact { name: "field.vtu".into(), text: "fixture".into() }],
        });
        let key = || InspectionKey {
            uid: app.session.uid,
            project: app.session.doc.multiphysics.clone(),
            study: id.clone(),
            run: "run-fixture".into(),
            artifact: "field.vtu".into(),
        };
        let mut state = State::default();
        state.selected = Some(id.clone());
        state.stage = Stage::Results;
        let (send, done) = mpsc::channel();
        state.inspection = Some(InspectionJob { key: key(), cancelled: Arc::new(AtomicBool::new(false)), done });
        send.send(Ok(json!({"content":{"format":"vtu"}}).into())).unwrap();
        poll_inspection(&app, &mut state);
        assert!(state.inspected.is_some());

        let (send, done) = mpsc::channel();
        state.inspection = Some(InspectionJob { key: key(), cancelled: Arc::new(AtomicBool::new(false)), done });
        state.selected_artifact = 1;
        send.send(Ok(json!({"content":{"format":"vtu"}}).into())).unwrap();
        state.inspected = None;
        poll_inspection(&app, &mut state);
        assert!(state.inspected.is_none(), "a completion for another artifact must be discarded");

        state.selected_artifact = 0;
        let (send, done) = mpsc::channel();
        state.inspection = Some(InspectionJob { key: key(), cancelled: Arc::new(AtomicBool::new(false)), done });
        let project = &mut std::sync::Arc::make_mut(&mut app.session.doc).multiphysics;
        Arc::make_mut(project).studies[0].question = "Changed".into();
        send.send(Ok(json!({"content":{"format":"vtu"}}).into())).unwrap();
        poll_inspection(&app, &mut state);
        assert!(state.inspected.is_none(), "a completion for a changed project must be discarded");

        let (send, done) = mpsc::channel();
        let wrong_document = InspectionKey {
            uid: app.session.uid.wrapping_add(1),
            project: app.session.doc.multiphysics.clone(),
            study: id.clone(),
            run: "run-fixture".into(),
            artifact: "field.vtu".into(),
        };
        state.inspection = Some(InspectionJob { key: wrong_document, cancelled: Arc::new(AtomicBool::new(false)), done });
        send.send(Ok(json!({"content":{"format":"vtu"}}).into())).unwrap();
        poll_inspection(&app, &mut state);
        assert!(state.inspected.is_none(), "a completion from another document must be discarded");

        let (send, done) = mpsc::channel();
        let cancelled = Arc::new(AtomicBool::new(false));
        let key = InspectionKey {
            uid: app.session.uid,
            project: app.session.doc.multiphysics.clone(),
            study: id,
            run: "run-fixture".into(),
            artifact: "field.vtu".into(),
        };
        app.multiphysics_workspace.inspection = Some(InspectionJob { key, cancelled: cancelled.clone(), done });
        send.send(Ok(json!({}).into())).unwrap();
        close(&mut app);
        assert!(cancelled.load(std::sync::atomic::Ordering::Acquire), "closing requests inspection cancellation");
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn inspection_preview_limits_serialization_work() {
        let value = json!({"fields": ["x".repeat(100_000)]});
        assert!(bounded_metadata(&value, 512).len() <= 512);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn scientific_navigation_preserves_design_pending_fit() {
        use solidcraft_engine::multiphysics::vtk::{Association, BoxGrid, Field};
        let grid = BoxGrid {
            bounds: [[0.0, 1.0]; 3],
            cells: vec![[[0.0, 1.0]; 3]],
            planes: [vec![0.5], vec![0.5], vec![0.5]],
            fields: vec![Field {
                name: "test_vector".into(),
                association: Association::Cell,
                data_type: "Float64".into(),
                components: 3,
                values: vec![1.0, 0.0, 0.0],
            }],
        };
        let mut app = SolveApp::new(Default::default(), Default::default());
        let camera = app.cam;
        app.pending_fit = Some((
            camera,
            solidcraft_engine::geom::Aabb3 { min: solidcore_record::Vec3::ZERO, max: solidcore_record::Vec3::new(1.0, 1.0, 1.0) },
            false,
            false,
        ));
        let pending = app.pending_fit;
        let mut state = crate::multiphysics_field_view::State3d::default();
        let ctx = crate::fonts::test_context();
        let bounds = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(800.0, 600.0));
        let mut before_zoom = 0.0;
        for frame in 0..3 {
            let mut events = vec![egui::Event::PointerMoved(egui::pos2(400.0, 200.0))];
            if frame == 1 {
                before_zoom = state.camera.distance;
                events.push(egui::Event::MouseWheel {
                    source: egui::MouseWheelSource::Unknown,
                    unit: egui::MouseWheelUnit::Point,
                    delta: egui::vec2(0.0, 120.0),
                    modifiers: egui::Modifiers::NONE,
                    phase: egui::TouchPhase::Move,
                });
            }
            let _ = ctx.run_ui(egui::RawInput { screen_rect: Some(bounds), time: Some(frame as f64 / 60.0), events, ..Default::default() }, |ui| {
                Container::for_layout(LayoutMode::Docked, Surface::Viewport).central().show(ui, |ui| {
                    crate::multiphysics_field_view::show(&mut app, ui, &mut state, &grid, 0, 2, 0.5, 400.0).unwrap();
                });
            });
            assert_eq!(app.cam, camera);
            assert_eq!(app.pending_fit, pending);
        }
        assert!(state.camera.distance < before_zoom, "the preserved pending fit must coexist with actual field zoom");
    }

    #[test]
    fn material_library_mouse_import_edit_undo_and_reopen_preserve_study() {
        let mut h = Harness::new();
        h.app.run("multiphysics.study.create", json!({"template":"poisson1d","name":"Unchanged numerical study"})).unwrap();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        let project = h.app.session.doc.multiphysics.clone();
        let camera = h.app.cam;
        let assignments = h.app.session.doc.materials.clone();
        let selection = h.app.session.selection.clone();
        let before = h.app.session.undo.len();
        h.step(&json!({"click":{"handle":"multiphysics:materials"}})).unwrap();
        material_library_text(&mut h, "multiphysics:materials_search", "C11000");
        h.step(&json!({"click":{"handle":"multiphysics:materials_profile:c11000"}})).unwrap();
        material_library_text(&mut h, "multiphysics:materials_name", "Enclosure copper reference");
        let plan = h
            .app
            .run("render.library.physical.prepare", json!({"id":"c11000","name":"Enclosure copper reference","acknowledge_reference":true}))
            .unwrap();
        assert_eq!(h.app.session.undo.len(), before, "browsing and preparation must be read-only");
        h.step(&json!({"click":{"handle":"multiphysics:materials_create"}})).unwrap();
        assert!(h.app.session.doc.user_materials.is_empty(), "an unchecked review must not import");
        h.step(&json!({"click":{"handle":"multiphysics:materials_acknowledge"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:materials_create"}})).unwrap();
        let material = h.app.session.doc.user_material("Enclosure copper reference").unwrap();
        assert_eq!(serde_json::to_value(material.physical_card.as_ref().unwrap()).unwrap(), plan["command"]["params"]["physical_card"]);
        assert_eq!(material.density, plan["command"]["params"]["density_g_cm3"].as_f64().unwrap());
        assert!(material.physical_card.as_ref().unwrap().card().properties.iter().any(|p| p.data.kind() == "thermal_conductivity"));
        assert_eq!(h.app.session.undo.len(), before + 1);
        assert_eq!(h.app.multiphysics_workspace.materials.selected_design.as_deref(), Some("Enclosure copper reference"));
        h.step(&json!({"click":{"handle":"multiphysics:undo"}})).unwrap();
        assert!(h.app.session.doc.user_materials.is_empty());
        assert!(h.app.multiphysics_workspace.materials.selected_design.is_none());
        h.step(&json!({"click":{"handle":"multiphysics:redo"}})).unwrap();
        assert!(h.app.session.doc.user_material("Enclosure copper reference").is_some());
        h.step(&json!({"click":{"handle":"material:edit"}})).unwrap();
        material_library_text(&mut h, "material:field:name", "Reviewed copper reference");
        h.step(&json!({"click":{"handle":"material:save"}})).unwrap();
        assert!(h.app.session.doc.user_material("Enclosure copper reference").is_none());
        assert!(h.app.session.doc.user_material("Reviewed copper reference").is_some());
        let saved = h.app.session.doc.to_json();
        let reopened = solidcraft_engine::doc::Document::from_json(&saved).unwrap();
        assert_eq!(reopened.user_materials, h.app.session.doc.user_materials);
        assert_eq!(reopened.multiphysics, h.app.session.doc.multiphysics);
        assert_eq!(h.app.cam, camera);
        assert_eq!(h.app.session.selection, selection);
        assert_eq!(h.app.session.doc.materials, assignments);
        assert!(std::sync::Arc::ptr_eq(&project, &h.app.session.doc.multiphysics));
        let evidence = std::env::var("SOLIDCRAFT_MULTIPHYSICS_EVIDENCE").ok();
        for (width, height) in [(1440., 900.), (1280., 800.)] {
            render_frame(&mut h.app, width, height, evidence.as_deref(), "materials-design");
        }
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:materials_references"}})).unwrap();
        for (width, height) in [(1440., 900.), (1280., 800.)] {
            render_frame(&mut h.app, width, height, evidence.as_deref(), "materials-reference");
            let content = crate::scenario::rect_of("multiphysics:content").unwrap();
            for target in [
                "multiphysics:materials_search",
                "multiphysics:materials_name",
                "multiphysics:materials_acknowledge",
                "multiphysics:materials_create",
            ] {
                let rect = crate::scenario::rect_of(target).unwrap();
                assert!(content.expand(1.).contains_rect(rect), "{target}: {rect:?} outside {content:?}");
            }
        }
        h.frames(3);
        material_library_text(&mut h, "multiphysics:materials_name", "Reviewed copper reference");
        h.step(&json!({"click":{"handle":"multiphysics:materials_acknowledge"}})).unwrap();
        let document = h.app.session.doc.clone();
        let history = h.app.session.undo.len();
        h.step(&json!({"click":{"handle":"multiphysics:materials_create"}})).unwrap();
        assert!(std::sync::Arc::ptr_eq(&document, &h.app.session.doc), "duplicate refusal must be atomic");
        assert_eq!(history, h.app.session.undo.len());
    }

    #[test]
    fn material_library_replacing_document_clears_pending_editor() {
        let mut h = Harness::new();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:materials"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:materials_design"}})).unwrap();
        h.step(&json!({"click":{"handle":"material:new"}})).unwrap();
        assert!(crate::scenario::handle_at("material:save").is_some());
        material_library_text(&mut h, "material:field:name", "Unsaved old-document material");
        h.app.session = solidcraft_engine::Session::default();
        h.frames(3);
        assert!(h.app.multiphysics_workspace.materials.selected_design.is_none());
        assert!(h.app.session.doc.user_materials.is_empty());
        assert!(crate::scenario::handle_at("material:save").is_none());
        assert!(crate::scenario::rect_of("multiphysics:materials_create").is_none());
        assert!(h.app.session.doc.user_materials.is_empty());
    }

    fn material_library_text(h: &mut Harness, handle: &str, value: &str) {
        h.step(&json!({"click":{"handle":handle},"double":true})).unwrap();
        h.step(&json!({"key":"A","ctrl":true})).unwrap();
        h.step(&json!({"text":value})).unwrap();
        h.step(&json!({"key":"Enter"})).unwrap();
    }

    #[test]
    fn numeric_preview_preserves_missing_samples_and_bounds_work() {
        let table = NumericTable::parse("x,T_K\n0,300\n1,NaN\n2,\n3,305\n");
        assert_eq!(table.headers, ["x", "T_K"]);
        assert_eq!(table.rows.len(), 4);
        assert_eq!(table.rows[1][1], None);
        assert_eq!(table.rows[2][1], None);
        assert_eq!(table.rows[3][1], Some(305.0));
        let huge = format!("x,y\n{}", "0,1\n".repeat(MAX_TABLE_ROWS + 1));
        assert_eq!(NumericTable::parse(&huge).rows.len(), MAX_TABLE_ROWS);
    }

    #[test]
    fn workspace_renders_bounded_native_frames_at_desktop_sizes() {
        let mut app = SolveApp::new(Default::default(), Default::default());
        let reply = app.run("multiphysics.study.create", json!({"template":"poisson1d","name":"Poisson 1D · source = 2"})).unwrap();
        let study = app.session.doc.multiphysics.study(reply["study"]["id"].as_str().unwrap()).unwrap().clone();
        switch(&mut app, "multiphysics").unwrap();
        app.multiphysics_workspace.selected = Some(study.id.clone());
        app.multiphysics_workspace.stage = Stage::Physics;
        app.multiphysics_workspace.drafts.insert(study.id.clone(), Draft::from_study(&study));
        app.multiphysics_workspace.drafts.get_mut(&study.id).unwrap().question = "Check the scalar field under a uniform source".into();
        let evidence_dir = std::env::var("SOLIDCRAFT_MULTIPHYSICS_EVIDENCE").ok();
        for (width, height) in [(1440., 900.), (1280., 800.)] {
            render_frame(&mut app, width, height, evidence_dir.as_deref(), "authoring");
        }
        if let Ok(folder) = std::env::var("SOLIDCRAFT_MULTIPHYSICS_RESULT") {
            let folder = std::path::Path::new(&folder);
            let manifest: Value = serde_json::from_str(&std::fs::read_to_string(folder.join("manifest.json")).unwrap()).unwrap();
            let artifacts = manifest["files"]
                .as_array()
                .unwrap()
                .iter()
                .map(|file| {
                    let name = file.as_str().unwrap().to_string();
                    solidcraft_engine::multiphysics::Artifact { text: std::fs::read_to_string(folder.join(&name)).unwrap(), name }
                })
                .collect();
            let evidence = solidcraft_engine::multiphysics::RunEvidence {
                id: "rendered-poisson1d-run".into(),
                study_revision: study.revision,
                case_hash: solidcraft_engine::multiphysics::case_hash(&study.case).unwrap(),
                input_case: study.case.clone(),
                manifest,
                artifacts,
            };
            let project = &mut std::sync::Arc::make_mut(&mut app.session.doc).multiphysics;
            std::sync::Arc::make_mut(project).attach_evidence(&study.id, evidence).unwrap();
            app.multiphysics_workspace.stage = Stage::Results;
            app.multiphysics_workspace.drafts.insert(study.id.clone(), Draft::from_study(&study));
            app.multiphysics_workspace.selected_artifact =
                app.session.doc.multiphysics.studies[0].runs[0].artifacts.iter().position(|artifact| artifact.name.ends_with(".csv")).unwrap_or(0);
            for (width, height) in [(1440., 900.), (1280., 800.)] {
                render_frame(&mut app, width, height, evidence_dir.as_deref(), "results");
            }
            let committed = app.session.doc.multiphysics.study(&study.id).unwrap().clone();
            let mut state = std::mem::take(&mut app.multiphysics_workspace);
            #[cfg(not(target_arch = "wasm32"))]
            {
                state.selected_artifact = committed.runs[0].artifacts.iter().position(|artifact| artifact.name.ends_with(".vtu")).unwrap();
                let artifact = committed.runs[0].artifacts[state.selected_artifact].name.clone();
                start_inspection(&app, &mut state, &committed.id, &committed.runs[0].id, &artifact);
                assert!(state.inspection.is_some(), "inspection should start off the window thread");
                let deadline = std::time::Instant::now() + std::time::Duration::from_secs(20);
                while state.inspection.is_some() {
                    poll_inspection(&app, &mut state);
                    assert!(std::time::Instant::now() < deadline, "inspection worker did not finish");
                    std::thread::sleep(std::time::Duration::from_millis(5));
                }
                let inspected = state.inspected.as_ref().expect("current inspection must be delivered");
                let result = inspected.outcome.as_ref().unwrap();
                assert_eq!(result.metadata["content"]["format"], "vtu");
                assert_eq!(result.metadata["content"]["points"], 65);
            }
            state.drafts.get_mut(&study.id).unwrap().question = "Assess the retained solution".into();
            apply(&mut app, &mut state, &committed);
            assert!(state.error.is_none(), "{:?}", state.error);
            let saved = app.session.doc.multiphysics.study(&study.id).unwrap();
            let draft = state.drafts.get(&study.id).unwrap();
            assert!(saved.revision > committed.revision);
            assert_eq!(draft.expected_revision, saved.revision);
            assert!(!draft.dirty(), "successful Apply after a retained run must clear its draft");
            assert_eq!(saved.runs.len(), 1, "Apply must retain existing evidence");
            app.multiphysics_workspace = state;
        }
    }

    /// Render genuine retained solver arrays through the actual workspace widgets.
    /// This is raster-harness evidence, not a mouse-driven native application capture.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    #[ignore = "requires SOLIDCRAFT_ENCLOSURE_RESULT containing a real solver run"]
    fn real_enclosure_result_renders_spatial_slices() {
        let folder = std::env::var("SOLIDCRAFT_ENCLOSURE_RESULT").expect("set SOLIDCRAFT_ENCLOSURE_RESULT");
        let folder = std::path::Path::new(&folder);
        let manifest: Value = serde_json::from_str(&std::fs::read_to_string(folder.join("manifest.json")).unwrap()).unwrap();
        let artifacts = manifest["files"]
            .as_array()
            .unwrap()
            .iter()
            .map(|file| {
                let name = file.as_str().unwrap().to_string();
                solidcraft_engine::multiphysics::Artifact { text: std::fs::read_to_string(folder.join(&name)).unwrap(), name }
            })
            .collect();
        let mut app = SolveApp::new(Default::default(), Default::default());
        let reply = app
            .run(
                "multiphysics.study.create",
                json!({"template":"flow-electronics-enclosure","name":"Electronics enclosure · prescribed-flow demonstration"}),
            )
            .unwrap();
        let id = reply["study"]["id"].as_str().unwrap().to_string();
        let revision = app.session.doc.multiphysics.study(&id).unwrap().revision;
        app.run(
            "multiphysics.study.update",
            json!({"id":id,"expected_revision":revision,
            "name":"Electronics enclosure · prescribed-flow demonstration",
            "question":"Isothermal primitive blocks; inspect conservation before engineering use.",
            "case":manifest["case"]}),
        )
        .unwrap();
        let study = app.session.doc.multiphysics.study(&id).unwrap().clone();
        let evidence = solidcraft_engine::multiphysics::RunEvidence {
            id: "real-enclosure-run".into(),
            study_revision: study.revision,
            case_hash: solidcraft_engine::multiphysics::case_hash(&study.case).unwrap(),
            input_case: study.case.clone(),
            manifest,
            artifacts,
        };
        solidcraft_engine::multiphysics::validate_evidence(&study, &evidence).unwrap();
        let project = &mut std::sync::Arc::make_mut(&mut app.session.doc).multiphysics;
        std::sync::Arc::make_mut(project).attach_evidence(&id, evidence).unwrap();
        switch(&mut app, "multiphysics").unwrap();
        let retained = app.session.doc.multiphysics.study(&id).unwrap().runs[0].clone();
        let selected_artifact = retained.artifacts.iter().rposition(|artifact| artifact.name.ends_with(".vtu")).unwrap();
        let artifact = &retained.artifacts[selected_artifact].name;
        let mut state = std::mem::take(&mut app.multiphysics_workspace);
        state.selected = Some(id.clone());
        state.stage = Stage::Results;
        state.selected_run = Some(retained.id.clone());
        state.shown_run = Some(retained.id.clone());
        state.selected_artifact = selected_artifact;
        start_inspection(&app, &mut state, &id, &retained.id, artifact);
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(20);
        while state.inspection.is_some() {
            poll_inspection(&app, &mut state);
            assert!(std::time::Instant::now() < deadline, "inspection worker did not finish");
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        let data = state.inspected.as_ref().unwrap().outcome.as_ref().unwrap();
        let grid = data.spatial.as_ref().unwrap().as_ref().unwrap();
        assert!(grid.cells.len() > 1);
        let velocity_field = grid.fields.iter().position(|field| field.name == "velocity_m_per_s").unwrap();
        let pressure_field = grid.fields.iter().position(|field| field.name == "kinematic_pressure_m2_per_s2").unwrap();
        let pvd = retained.artifacts.iter().find(|artifact| artifact.name.ends_with(".pvd")).unwrap();
        let expected_coordinate =
            solidcraft_engine::multiphysics::vtk::parse_pvd(&pvd.text).unwrap().into_iter().find(|frame| &frame.file == artifact).unwrap().coordinate;
        assert_eq!(data.frame_coordinate, Some(expected_coordinate));
        state.view.slice_field = velocity_field;
        state.view.slice_axis = 2;
        state.view.slice_layer = 2.min(grid.planes[2].len().saturating_sub(1));
        app.multiphysics_workspace = state;
        let before_case = app.session.doc.multiphysics.study(&id).unwrap().case.clone();
        let mut harness = Harness::new();
        harness.app = app;
        harness.frames(3);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_plane:1"}})).unwrap();
        assert_eq!(harness.app.multiphysics_workspace.view.slice_axis, 1);
        assert_eq!(harness.app.multiphysics_workspace.view.slice_layer, 0);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_layer"},"double":true})).unwrap();
        harness.step(&json!({"key":"A","ctrl":true})).unwrap();
        harness.step(&json!({"text":"2"})).unwrap();
        harness.step(&json!({"key":"Enter"})).unwrap();
        assert_eq!(harness.app.multiphysics_workspace.view.slice_layer, 2);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_field_menu"}})).unwrap();
        harness.step(&json!({"click":{"handle":format!("multiphysics:slice_field:{pressure_field}")}})).unwrap();
        assert_eq!(harness.app.multiphysics_workspace.view.slice_field, pressure_field);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_field_menu"}})).unwrap();
        harness.step(&json!({"click":{"handle":format!("multiphysics:slice_field:{velocity_field}")}})).unwrap();
        assert_eq!(harness.app.multiphysics_workspace.view.slice_field, velocity_field);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_plane:0"}})).unwrap();
        assert_eq!(harness.app.multiphysics_workspace.view.slice_axis, 0);
        assert_eq!(harness.app.multiphysics_workspace.view.slice_layer, 0);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_plane:2"}})).unwrap();
        assert_eq!(harness.app.multiphysics_workspace.view.slice_axis, 2);
        assert_eq!(harness.app.session.doc.multiphysics.study(&id).unwrap().case, before_case, "view controls must not edit solver inputs");
        harness.app.multiphysics_workspace.view.slice_layer = 2;
        let evidence_dir = std::env::var("SOLIDCRAFT_MULTIPHYSICS_EVIDENCE").ok();
        for (width, height) in [(1440., 900.), (1280., 800.)] {
            render_frame(&mut harness.app, width, height, evidence_dir.as_deref(), "enclosure-slice");
            let slice = crate::scenario::rect_of("multiphysics:field_slice").unwrap();
            let legend = crate::scenario::rect_of("multiphysics:field_legend").unwrap();
            let data = crate::scenario::rect_of("multiphysics:field_data").unwrap();
            let content = crate::scenario::rect_of("multiphysics:content").unwrap();
            assert!(slice.is_finite() && slice.width() > 400. && slice.height() >= 220.);
            assert!(data.width() >= 350. && data.height() >= 200., "field view is too small: {data:?}");
            for (name, rect) in [("slice", slice), ("legend", legend), ("data", data)] {
                assert!(content.expand(1.0).contains_rect(rect), "{name}: {rect:?} outside {content:?}");
            }
        }
        harness.frames(3);
        harness.app.preferences.nav = "fusion".into();
        harness.app.ui.view_angles_lock_yaw = true;
        harness.app.ui.view_angles_lock_pitch = true;
        harness.app.ui.view_angles_lock_roll = true;
        let design_camera = harness.app.cam;
        let document = harness.app.session.doc.clone();
        let original_pivot = harness.app.preferences.orbit_around;
        let original_middle = harness.app.preferences.middle_double_click.clone();
        let original_zoom_cursor = harness.app.preferences.zoom_to_cursor;
        let original_zoom_surface = harness.app.preferences.zoom_stop_at_model;
        harness.step(&json!({"click":{"handle":"multiphysics:view_3d"}})).unwrap();
        assert!(harness.app.multiphysics_workspace.view.show_3d);
        let before_orbit = harness.app.multiphysics_workspace.view.field_3d.camera;
        let at = crate::scenario::rect_of("multiphysics:field_slice").unwrap().center();
        let orbit = harness.call("ui.drag", json!({"x0":at.x,"y0":at.y,"x1":at.x+50.0,"y1":at.y+24.0,"button":"middle","shift":true,"steps":8}));
        assert_eq!(orbit["ok"], true, "{orbit}");
        let after_orbit = harness.app.multiphysics_workspace.view.field_3d.camera;
        assert!(before_orbit.yaw != after_orbit.yaw || before_orbit.pitch != after_orbit.pitch, "3D orbit must move its own camera");
        let pan = harness.call("ui.drag", json!({"x0":at.x,"y0":at.y,"x1":at.x+30.0,"y1":at.y-14.0,"button":"middle","steps":8}));
        assert_eq!(pan["ok"], true, "{pan}");
        let after_pan = harness.app.multiphysics_workspace.view.field_3d.camera;
        assert_ne!(after_pan.target, after_orbit.target, "3D pan must move the field view");
        let zoom = harness.call("ui.scroll", json!({"x":at.x,"y":at.y,"delta":120.0}));
        assert_eq!(zoom["ok"], true, "{zoom}");
        assert!(harness.app.multiphysics_workspace.view.field_3d.camera.distance < after_pan.distance, "3D wheel must zoom its field camera");
        let before_glyphs = harness.app.multiphysics_workspace.view.field_3d.glyphs;
        harness.step(&json!({"click":{"handle":"multiphysics:3d_glyphs"}})).unwrap();
        assert_ne!(harness.app.multiphysics_workspace.view.field_3d.glyphs, before_glyphs);
        let glyph_count = crate::scenario::count_of("multiphysics:3d_glyph_count").unwrap();
        assert!(glyph_count > 0.0 && glyph_count <= 128.0, "vector sampling must be visible and bounded: {glyph_count}");
        let before_fit = harness.app.multiphysics_workspace.view.field_3d.camera;
        harness.step(&json!({"click":{"handle":"multiphysics:3d_fit"}})).unwrap();
        let fitted = harness.app.multiphysics_workspace.view.field_3d.camera;
        assert_eq!(fitted.target, before_orbit.target, "Fit must restore the domain centre");
        assert_ne!(fitted.target, after_pan.target, "Fit must undo the displayed pan");
        assert_eq!((fitted.yaw, fitted.pitch), (before_fit.yaw, before_fit.pitch), "Fit must preserve the chosen orientation");
        assert_eq!(harness.app.cam, design_camera, "scientific navigation must preserve the Design camera");
        assert!(std::sync::Arc::ptr_eq(&harness.app.session.doc, &document), "scientific navigation must not modify the document");
        assert_eq!(harness.app.preferences.orbit_around, original_pivot);
        assert_eq!(harness.app.preferences.middle_double_click, original_middle);
        assert_eq!(harness.app.preferences.zoom_to_cursor, original_zoom_cursor);
        assert_eq!(harness.app.preferences.zoom_stop_at_model, original_zoom_surface);
        assert!(harness.app.ui.view_angles_lock_yaw && harness.app.ui.view_angles_lock_pitch && harness.app.ui.view_angles_lock_roll);
        assert_eq!(harness.app.session.doc.multiphysics.study(&id).unwrap().case, before_case);
        for (width, height) in [(1440., 900.), (1280., 800.)] {
            render_frame(&mut harness.app, width, height, evidence_dir.as_deref(), "enclosure-3d");
            let canvas = crate::scenario::rect_of("multiphysics:field_slice").unwrap();
            let legend = crate::scenario::rect_of("multiphysics:field_legend").unwrap();
            let content = crate::scenario::rect_of("multiphysics:content").unwrap();
            assert!(canvas.width() > 400. && canvas.height() >= 180.);
            assert!(content.expand(1.0).contains_rect(canvas));
            assert!(content.expand(1.0).contains_rect(legend));
            assert_eq!(harness.app.cam, design_camera);
        }
    }

    fn render_frame(app: &mut SolveApp, width: f32, height: f32, evidence: Option<&str>, stage: &str) {
        let ctx = crate::fonts::test_context();
        let mut textures = std::collections::HashMap::new();
        let bounds = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(width, height));
        for frame in 0..5 {
            let mut output =
                ctx.run_ui(egui::RawInput { screen_rect: Some(bounds), time: Some(frame as f64 / 60.), ..Default::default() }, |ui| app.ui(ui));
            crate::ui_kit::test_textures_update(&mut textures, &output.textures_delta);
            output.textures_delta.clear();
            let content = crate::scenario::rect_of("multiphysics:content").unwrap();
            assert!(content.is_finite() && content.width() > 300. && content.height() > 400., "{content:?}");
            if stage == "results" {
                for target in ["multiphysics:x_axis", "multiphysics:y_axis", "multiphysics:plot", "multiphysics:plot_last_sample"] {
                    let rect = crate::scenario::rect_of(target).unwrap();
                    assert!(
                        rect.is_finite() && rect.left() >= content.left() - 1. && rect.right() <= content.right() + 1.,
                        "{target}: {rect:?} outside {content:?}"
                    );
                }
            }
            if frame == 4
                && let Some(folder) = evidence
            {
                std::fs::create_dir_all(folder).unwrap();
                let image = crate::ui_kit::test_raster(&ctx, output.shapes, output.pixels_per_point, bounds, &textures);
                image.save(std::path::Path::new(folder).join(format!("multiphysics-{stage}-{}x{}.png", width as u32, height as u32))).unwrap();
            }
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn export_uses_frozen_case_without_blocking_edits_and_preserves_scope() {
        let mut h = Harness::new();
        let reply = h.app.run("multiphysics.study.create", json!({"template":"poisson1d","name":"Original case"})).unwrap();
        let study = h.app.session.doc.multiphysics.study(reply["study"]["id"].as_str().unwrap()).unwrap().clone();
        let path = h.dir.join("frozen-case.json");
        let duplicate = h.dir.join("duplicate-case.json");
        let mut state = State::default();
        let started = command(&mut h.app, &mut state, "multiphysics.case.export", json!({"id":study.id,"path":path})).unwrap();
        assert_eq!(started["queued"], true);
        assert!(state.export_busy());
        assert!(command(&mut h.app, &mut state, "multiphysics.case.export", json!({"id":study.id,"path":duplicate})).is_none());
        let mut changed_case = study.case.clone();
        changed_case["source"] = json!(9.0);
        h.app
            .run(
                "multiphysics.study.update",
                json!({"id":study.id,"expected_revision":study.revision,"name":"Edited while exporting","question":"","case":changed_case}),
            )
            .unwrap();
        let history = h.app.session.undo.len();
        h.app.multiphysics_workspace = state;
        close(&mut h.app);
        let mut state = std::mem::take(&mut h.app.multiphysics_workspace);
        assert!(!state.export_busy());
        assert!(state.error.is_none(), "{:?}", state.error);
        let written: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        assert_eq!(written, study.case, "export freezes the original case");
        assert_eq!(h.app.session.doc.multiphysics.study(&study.id).unwrap().case["source"], 9.0);
        assert_eq!(h.app.session.undo.len(), history, "export completion never adopts its scratch session");
        assert!(!duplicate.exists());
        assert!(state.notice.is_some());
        let outside = h.dir.with_extension("outside.json");
        let scope = solidcraft_engine::io::vfs::PathScope::new([h.dir.clone()]);
        solidcraft_engine::io::vfs::scoped(Some(scope), || {
            command(&mut h.app, &mut state, "multiphysics.case.export", json!({"id":study.id,"path":outside})).unwrap();
        });
        let export = state.export.take().unwrap();
        let outcome = export.done.recv().unwrap();
        assert!(outcome.is_err(), "worker must retain the caller's VFS scope");
        finish_export(&mut h.app, &mut state, export.uid, export.command, outcome);
        assert!(state.error.is_some());
        assert!(!outside.exists());
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn replacing_document_keeps_export_slot_until_retirement_and_discards_old_notice() {
        let mut app = SolveApp::new(Default::default(), Default::default());
        let old_uid = app.session.uid;
        let (send, done) = std::sync::mpsc::channel();
        app.multiphysics_workspace.export = Some(ExportJob { uid: old_uid, command: "multiphysics.case.export", done });
        app.run("file.new", json!({"name":"New document"})).unwrap();
        assert_ne!(app.session.uid, old_uid);
        assert!(app.multiphysics_workspace.export_busy());
        let log = app.session.log.clone();
        send.send(Ok(json!({"files":["old-case.json"]}))).unwrap();
        tick(&mut app, &crate::fonts::test_context());
        assert!(!app.multiphysics_workspace.export_busy());
        assert!(app.multiphysics_workspace.notice.is_none());
        assert!(app.multiphysics_workspace.error.is_none());
        assert_eq!(app.session.log, log, "an old completion must not echo into another document");
        assert!(app.session.doc.multiphysics.is_empty());
    }

    #[test]
    fn document_tabs_keep_same_id_study_drafts_with_their_document() {
        let mut h = Harness::new();
        let first = h.app.run("multiphysics.study.create", json!({"template":"poisson1d","name":"First document study"})).unwrap();
        let first_id = first["study"]["id"].as_str().unwrap().to_string();
        let first_uid = h.app.session.uid;
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(1);
        h.app.multiphysics_workspace.drafts.get_mut(&first_id).unwrap().case_text = "{unfinished first document".into();
        let mut second = solidcraft_engine::Session::default();
        let added = second.execute("multiphysics.study.create", &json!({"template":"poisson2d","name":"Second document study"})).unwrap();
        assert_eq!(added["study"]["id"], first_id, "the identity collision exercises document ownership");
        crate::documents::open_tab(&mut h.app, second).unwrap();
        let second_uid = h.app.session.uid;
        assert_ne!(first_uid, second_uid);
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(1);
        let draft = h.app.multiphysics_workspace.drafts.get_mut(&first_id).unwrap();
        assert!(serde_json::from_str::<Value>(&draft.case_text).is_ok());
        draft.question = "Second document draft".into();
        assert_eq!(h.app.multiphysics_workspace.y_column, 1);
        crate::documents::switch(&mut h.app, 0);
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(1);
        assert_eq!(h.app.session.uid, first_uid);
        assert_eq!(h.app.multiphysics_workspace.drafts.get(&first_id).unwrap().case_text, "{unfinished first document");
        crate::documents::switch(&mut h.app, 1);
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(1);
        assert_eq!(h.app.session.uid, second_uid);
        assert_eq!(h.app.multiphysics_workspace.drafts.get(&first_id).unwrap().question, "Second document draft");
        let clean = solidcraft_engine::Session::new(h.app.session.doc.as_ref().clone());
        let app = SolveApp::new(clean, Default::default());
        assert!(!crate::documents::active_is_blank(&app), "a study-only document contains authored work even without CAD features");
    }

    #[test]
    fn workspace_pointer_creates_study_and_retains_invalid_draft() {
        let mut h = Harness::new();
        h.call("ui.resize", json!({"width":1440,"height":900}));
        h.step(&json!({"click":{"handle":"ribbon:workspace"}})).unwrap();
        h.step(&json!({"click":{"handle":"workspace:menu:multiphysics"}})).unwrap();
        assert!(active(&h.app));
        assert!(crate::workspace::visible_tabs(&h.app).is_empty());
        h.step(&json!({"click":{"handle":"multiphysics:create"}})).unwrap();
        h.frames(2);
        assert_eq!(h.app.session.doc.multiphysics.studies.len(), 1);
        let id = h.app.session.doc.multiphysics.studies[0].id.clone();
        h.step(&json!({"click":{"handle":"multiphysics.stage_physics"}})).unwrap();
        h.step(&json!({"widget":"JSON editor"})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:case_json"}})).unwrap();
        h.step(&json!({"key":"A","cmd":true})).unwrap();
        h.step(&json!({"text":"{invalid json"})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:apply"}})).unwrap();
        assert_eq!(h.app.multiphysics_workspace.drafts.get(&id).unwrap().case_text, "{invalid json");
        assert!(h.app.multiphysics_workspace.error.is_some());
        h.step(&json!({"click":{"handle":"multiphysics:cancel_draft"}})).unwrap();
        assert!(!h.app.multiphysics_workspace.drafts.get(&id).unwrap().dirty());
        let camera = h.app.cam;
        switch(&mut h.app, "design").unwrap();
        assert_eq!(h.app.cam, camera);
        switch(&mut h.app, "multiphysics").unwrap();
        assert_eq!(h.app.session.doc.multiphysics.studies.len(), 1);
    }
}
