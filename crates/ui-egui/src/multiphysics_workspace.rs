//! Numerical studies with independent authoring, jobs and post-processing.
//! The ordinary Simulation placeholder and Design/Render workspace state are not reused.

use std::collections::BTreeMap;

use egui::{RichText, Sense};
use serde_json::{Value, json};
use solidcraft_engine::multiphysics::{Study, Template};
use solidcraft_theme::Tokens;

use crate::SolveApp;
use crate::ui_kit::{BackingTooltipExt as _, Container, Edge, LayoutMode, Surface};

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
impl Stage {
    const ALL: [Self; 7] = [Self::Prepare, Self::Physics, Self::Mesh, Self::Couple, Self::Run, Self::Results, Self::Verify];
    fn key(self) -> &'static str {
        match self {
            Self::Prepare => "multiphysics.stage_prepare",
            Self::Physics => "multiphysics.stage_physics",
            Self::Mesh => "multiphysics.stage_mesh",
            Self::Couple => "multiphysics.stage_couple",
            Self::Run => "multiphysics.stage_run",
            Self::Results => "multiphysics.stage_results",
            Self::Verify => "multiphysics.stage_verify",
        }
    }
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
}
impl Draft {
    fn from_study(study: &Study) -> Self {
        Self {
            name: study.name.clone(),
            question: study.question.clone(),
            case_text: pretty(&study.case),
            expected_revision: study.revision,
            base_name: study.name.clone(),
            base_question: study.question.clone(),
            base_case: study.case.clone(),
        }
    }
    pub fn dirty(&self) -> bool {
        self.name != self.base_name
            || self.question != self.base_question
            || serde_json::from_str::<Value>(&self.case_text).ok().as_ref() != Some(&self.base_case)
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
            error: None,
            notice: None,
        }
    }
}
/// Window-wide catalogue caches and a single export slot stay outside document-tab state.
#[derive(Default)]
pub struct State {
    pub view: ViewState,
    templates: Option<Result<Vec<Template>, String>>,
    catalog: Option<Result<Value, String>>,
    template_index: usize,
    search: String,
    show_catalog: bool,
    show_browser: bool,
    #[cfg(not(target_arch = "wasm32"))]
    export: Option<ExportJob>,
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
        app.multiphysics_workspace = State { export, ..Default::default() };
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
        app.multiphysics_workspace = state;
    }
    if app.multiphysics_workspace.export_busy()
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
    let response = ui.add_enabled(enabled, egui::Button::new(tr(key)));
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
    let selected = state.selected.as_ref().and_then(|id| studies.iter().find(|study| &study.id == id));
    let t = Tokens::get();
    let host = Container::for_layout(LayoutMode::Docked, Surface::Chrome);
    host.edge("multiphysics_header", Edge::Top).exact_size(76.0).resizable(false).show(ui, |ui| header(app, &mut state, ui));
    Container::for_layout(LayoutMode::Docked, Surface::Panel)
        .edge("multiphysics_status", Edge::Bottom)
        .exact_size(t.control_height + t.space_4)
        .resizable(false)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(tr(if app.session.multiphysics_host_available() { "multiphysics.worker_ready" } else { "multiphysics.preparation_only" }));
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
    if state.show_browser {
        Container::for_layout(LayoutMode::Docked, Surface::Panel)
            .edge("multiphysics_browser", Edge::Left)
            .default_size(248.0)
            .size_range(210.0..=420.0)
            .resizable(true)
            .show(ui, |ui| browser(app, &mut state, studies, ui));
    }
    let wide = ui.available_width() >= 720.0;
    if wide {
        Container::for_layout(LayoutMode::Docked, Surface::Panel)
            .edge("multiphysics_inspector", Edge::Right)
            .default_size(290.0)
            .size_range(240.0..=440.0)
            .resizable(true)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().id_salt("multiphysics_inspector_scroll").show(ui, |ui| inspector(app, &mut state, selected, ui));
            });
    }
    Container::for_layout(LayoutMode::Docked, Surface::Viewport).central().show(ui, |ui| {
        crate::scenario::publish_rect("multiphysics:content", ui.max_rect());
        if !wide {
            egui::CollapsingHeader::new(tr("multiphysics.inspector")).default_open(false).show(ui, |ui| inspector(app, &mut state, selected, ui));
        }
        if state.show_catalog {
            catalog(&mut state, ui);
        } else if let Some(study) = selected {
            ui.horizontal(|ui| {
                ui.heading(&study.name);
                if state.drafts.get(&study.id).is_some_and(Draft::dirty) {
                    ui.colored_label(t.warning, tr("multiphysics.unapplied"));
                }
            });
            ui.label(RichText::new(&study.profile).color(t.text_dim));
            ui.separator();
            match state.stage {
                Stage::Prepare => prepare(app, &mut state, study, ui),
                Stage::Physics | Stage::Mesh => settings(&mut state, study, ui),
                Stage::Couple => {
                    ui.heading(tr("multiphysics.coupling"));
                    ui.label(tr("multiphysics.coupling_unavailable"));
                    ui.label(tr("multiphysics.coupling_preserved"));
                }
                Stage::Run => jobs(app, &mut state, study, ui),
                Stage::Results => results(app, &mut state, study, ui),
                Stage::Verify => verify(study, ui),
            }
        } else {
            ui.heading(tr("multiphysics.welcome"));
            ui.label(tr("multiphysics.create_hint"));
            ui.label(tr("multiphysics.cad_limit"));
        }
    });
    app.multiphysics_workspace = state;
}

fn header(app: &mut SolveApp, state: &mut State, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        let response = ui.button(tr("multiphysics.workspace"));
        crate::scenario::publish_target("ribbon:workspace", response.rect);
        if response.clicked() {
            crate::context_menu::toggle_for(app, response.rect, response.rect.left_bottom(), crate::context_menu::Target::Workspace);
        }
        ui.separator();
        if ui.selectable_label(state.show_browser, tr("multiphysics.studies")).clicked() {
            state.show_browser = !state.show_browser;
        }
        if ui.selectable_label(state.show_catalog, tr("multiphysics.catalog")).clicked() {
            state.show_catalog = !state.show_catalog;
        }
        if button(ui, "multiphysics.undo", "multiphysics:undo", true).clicked() {
            command(app, state, "edit.undo", json!({}));
        }
        if button(ui, "multiphysics.redo", "multiphysics:redo", true).clicked() {
            command(app, state, "edit.redo", json!({}));
        }
        if button(ui, "multiphysics.save", "multiphysics:save", true).clicked() {
            app.start("file.save");
        }
    });
    egui::ScrollArea::horizontal().id_salt("multiphysics_stages").show(ui, |ui| {
        ui.horizontal(|ui| {
            for stage in Stage::ALL {
                let response = ui.selectable_label(state.stage == stage && !state.show_catalog, tr(stage.key()));
                crate::scenario::publish_target(stage.key(), response.rect);
                if response.clicked() {
                    state.stage = stage;
                    state.show_catalog = false;
                }
            }
        });
    });
}

fn browser(app: &mut SolveApp, state: &mut State, studies: &[Study], ui: &mut egui::Ui) {
    ui.heading(tr("multiphysics.studies"));
    egui::ScrollArea::vertical().id_salt("multiphysics_study_browser").show(ui, |ui| {
        for study in studies {
            let dirty = state.drafts.get(&study.id).is_some_and(Draft::dirty);
            let title = if dirty { format!("{} •", study.name) } else { study.name.clone() };
            let response = ui.selectable_label(state.selected.as_ref() == Some(&study.id), title);
            crate::scenario::publish_target(&format!("multiphysics:study:{}", study.id), response.rect);
            if response.clicked() {
                state.selected = Some(study.id.clone());
                state.selected_run = None;
                state.selected_artifact = 0;
            }
        }
        ui.separator();
        ui.strong(tr("multiphysics.new_study"));
        let templates = state.templates.as_ref().and_then(|result| result.as_ref().ok()).cloned().unwrap_or_default();
        let names: Vec<_> = templates.iter().map(|template| template.name.clone()).collect();
        if !names.is_empty() {
            state.template_index = state.template_index.min(names.len().saturating_sub(1));
            let selected = names.get(state.template_index).map(String::as_str).unwrap_or_default();
            crate::ui_kit::dropdown_menu(ui, "multiphysics_template", selected, |ui| {
                for (index, name) in names.iter().enumerate() {
                    ui.selectable_value(&mut state.template_index, index, name);
                }
            });
            if button(ui, "multiphysics.create", "multiphysics:create", true).clicked()
                && let Some(template) = templates.get(state.template_index)
                && let Some(reply) = command(app, state, "multiphysics.study.create", json!({"template":template.id,"name":template.name}))
            {
                state.selected = reply.get("study").and_then(|study| study.get("id")).and_then(Value::as_str).map(str::to_string);
                state.stage = Stage::Prepare;
                state.show_catalog = false;
            }
        }
        if let Some(Err(error)) = &state.templates {
            ui.colored_label(Tokens::get().form_error, error);
        }
        ui.add_space(Tokens::get().space_8);
        ui.label(RichText::new(tr("multiphysics.case_scope")).small().color(Tokens::get().text_dim));
        ui.separator();
        if button(ui, "multiphysics.remove", "multiphysics:remove", state.selected.is_some()).clicked()
            && let Some(id) = state.selected.clone()
        {
            command(app, state, "multiphysics.study.remove", json!({"id":id}));
        }
    });
}

fn prepare(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical().id_salt("multiphysics_prepare").show(ui, |ui| {
        ui.heading(tr("multiphysics.intent"));
        if let Some(draft) = state.drafts.get_mut(&study.id) {
            ui.label(tr("multiphysics.question"));
            let response = ui.add(egui::TextEdit::multiline(&mut draft.question).desired_width(f32::INFINITY).desired_rows(3));
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

fn inspector(app: &mut SolveApp, state: &mut State, study: Option<&Study>, ui: &mut egui::Ui) {
    let Some(study) = study else { return };
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
    let report = solidcraft_engine::multiphysics::preflight(study);
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
    let ready = !draft_dirty && report.errors.is_empty() && available;
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

fn settings(state: &mut State, study: &Study, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.heading(tr(if state.stage == Stage::Mesh { "multiphysics.discretisation" } else { "multiphysics.case_settings" }));
        ui.checkbox(&mut state.json_editor, tr("multiphysics.json_editor"));
    });
    ui.label(tr(if state.stage == Stage::Mesh { "multiphysics.mesh_limit" } else { "multiphysics.settings_hint" }));
    let json_editor = state.json_editor;
    let Some(draft) = state.drafts.get_mut(&study.id) else { return };
    let parsed = serde_json::from_str::<Value>(&draft.case_text);
    let raw = json_editor || parsed.is_err();
    let editor_width = ui.available_width().max(240.0);
    egui::ScrollArea::both().id_salt(("multiphysics_settings", &study.id)).show(ui, |ui| {
        if raw {
            let response = ui.add(
                egui::TextEdit::multiline(&mut draft.case_text).code_editor().desired_width(editor_width).desired_rows(26).char_limit(1024 * 1024),
            );
            crate::scenario::publish_target("multiphysics:case_json", response.rect);
            if let Err(error) = &parsed {
                ui.colored_label(Tokens::get().form_error, error.to_string());
            }
        } else if let Ok(mut value) = parsed
            && fields(ui, &mut value, 0)
        {
            draft.case_text = pretty(&value);
        }
    });
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
                ui.push_id(name, |ui| {
                    if value.is_object() || value.is_array() {
                        egui::CollapsingHeader::new(name).default_open(depth == 0).show(ui, |ui| changed |= fields(ui, value, depth + 1));
                    } else {
                        crate::ui_kit::row(ui, name, |ui| match value {
                            Value::Bool(boolean) => changed |= ui.checkbox(boolean, "").changed(),
                            Value::Number(number) => {
                                if let Some(mut numeric) = number.as_u64() {
                                    if ui.add(egui::DragValue::new(&mut numeric)).changed() {
                                        *number = numeric.into();
                                        changed = true;
                                    }
                                } else if let Some(mut numeric) = number.as_i64() {
                                    if ui.add(egui::DragValue::new(&mut numeric)).changed() {
                                        *number = numeric.into();
                                        changed = true;
                                    }
                                } else if let Some(mut numeric) = number.as_f64() {
                                    let response = ui.add(egui::DragValue::new(&mut numeric).speed(0.001));
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
                            Value::String(string) => changed |= ui.text_edit_singleline(string).changed(),
                            _ => {
                                ui.label(tr("multiphysics.not_specified"));
                            }
                        });
                    }
                });
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

fn catalog(state: &mut State, ui: &mut egui::Ui) {
    ui.heading(tr("multiphysics.catalog"));
    ui.label(tr("multiphysics.catalog_hint"));
    crate::ui_kit::text_row(ui, tr("multiphysics.search"), &mut state.search, "");
    let search = state.search.to_lowercase();
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
            ui.selectable_value(&mut state.selected_run, Some(run.id.clone()), &run.id);
        }
    });
    if run.study_revision != study.revision {
        ui.colored_label(Tokens::get().warning, tr("multiphysics.stale_run"));
    }
    egui::CollapsingHeader::new(tr("multiphysics.manifest")).show(ui, |ui| {
        ui.monospace(pretty(&run.manifest).chars().take(16384).collect::<String>());
    });
    ui.horizontal(|ui| {
        ui.label(tr("multiphysics.export_directory"));
        ui.text_edit_singleline(&mut state.export_directory);
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
    if run.artifacts.is_empty() {
        ui.label(tr("multiphysics.no_artifacts"));
        return;
    }
    state.selected_artifact = state.selected_artifact.min(run.artifacts.len().saturating_sub(1));
    let Some(artifact) = run.artifacts.get(state.selected_artifact) else { return };
    crate::ui_kit::dropdown_menu(ui, "multiphysics_artifact", &artifact.name, |ui| {
        for (index, artifact) in run.artifacts.iter().enumerate() {
            ui.selectable_value(&mut state.selected_artifact, index, &artifact.name);
        }
    });
    if artifact.name.ends_with(".csv") {
        let table = NumericTable::parse(&artifact.text);
        table_view(state, &table, ui);
    } else {
        egui::ScrollArea::both().id_salt("multiphysics_artifact_text").show(ui, |ui| {
            ui.monospace(artifact.text.chars().take(32768).collect::<String>());
            ui.label(tr("multiphysics.preview_bounded"));
        });
    }
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
                ui.selectable_value(&mut state.x_column, index, name);
            }
        });
        crate::scenario::publish_rect("multiphysics:x_axis", ui.min_rect());
        let ui = &mut columns[1];
        ui.label(tr("multiphysics.y_axis"));
        let y = table.headers.get(state.y_column).cloned().unwrap_or_default();
        let width = ui.available_width();
        crate::ui_kit::bounded_dropdown_menu(ui, "multiphysics_y", &y, width, |ui| {
            for (index, name) in table.headers.iter().enumerate() {
                ui.selectable_value(&mut state.y_column, index, name);
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
    let t = Tokens::get();
    let visible_width = (ui.clip_rect().right() - ui.next_widget_position().x).max(1.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width().min(visible_width).max(1.0), 220.0), Sense::hover());
    crate::scenario::publish_rect("multiphysics:plot", rect);
    let bounds = points.iter().flatten().fold(None, |bounds, &(x, y)| {
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
    if let Some(point) = points.last().copied().flatten() {
        crate::scenario::publish_rect("multiphysics:plot_last_sample", egui::Rect::from_center_size(position(point), egui::Vec2::ZERO));
    }
    let mut previous = None;
    for point in points {
        let current = point.map(position);
        if let (Some(a), Some(b)) = (previous, current) {
            ui.painter().line_segment([a, b], egui::Stroke::new(1.5, t.accent));
        }
        previous = current;
    }
    ui.painter().text(area.left_top(), egui::Align2::LEFT_TOP, format!("{ymax:.4e}"), egui::FontId::proportional(t.font_small), t.text_dim);
    ui.painter().text(area.left_bottom(), egui::Align2::LEFT_BOTTOM, format!("{xmin:.4e}"), egui::FontId::proportional(t.font_small), t.text_dim);
    ui.painter().text(area.right_bottom(), egui::Align2::RIGHT_BOTTOM, format!("{xmax:.4e}"), egui::FontId::proportional(t.font_small), t.text_dim);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scenario::Harness;

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
