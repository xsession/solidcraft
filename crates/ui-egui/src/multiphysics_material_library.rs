//! Multiphysics access to the existing physical reference and design material stores.
//! Only presentation state lives here; definitions, validation and Undo use material commands.
use std::sync::Arc;

use serde_json::{Value, json};
use solidcraft_theme::Tokens;
use solidcraft_ui_text::tr;

use crate::{SolveApp, material_editor, ui_kit};

#[derive(Clone, Debug)]
pub(crate) struct State {
    uid: Option<u64>,
    design: bool,
    query: String,
    loaded_query: Option<String>,
    entries: Vec<Value>,
    total: usize,
    selected_profile: Option<String>,
    detail: Option<Arc<Value>>,
    name: String,
    acknowledged: bool,
    error: Option<String>,
    note: Option<String>,
    pub(crate) selected_design: Option<String>,
    design_index: usize,
    editor: material_editor::Editor,
}

impl Default for State {
    fn default() -> Self {
        Self {
            uid: None,
            design: false,
            query: String::new(),
            loaded_query: None,
            entries: Vec::new(),
            total: 0,
            selected_profile: None,
            detail: None,
            name: String::new(),
            acknowledged: false,
            error: None,
            note: None,
            selected_design: None,
            design_index: usize::MAX,
            editor: Default::default(),
        }
    }
}

fn list(app: &mut SolveApp, state: &mut State) {
    state.loaded_query = Some(state.query.clone());
    match app.run("render.library.list", json!({"kind":"physical_profiles","query":state.query,"offset":0,"limit":128,"availability":"all"})) {
        Ok(reply) => {
            state.entries = reply["items"].as_array().map(|items| items.iter().take(128).cloned().collect()).unwrap_or_default();
            state.total = reply["total"].as_u64().and_then(|n| usize::try_from(n).ok()).unwrap_or(state.entries.len());
            if state.selected_profile.as_ref().is_some_and(|id| !state.entries.iter().any(|item| item["id"].as_str() == Some(id.as_str()))) {
                state.selected_profile = None;
                state.detail = None;
                state.name.clear();
                state.acknowledged = false;
            }
            state.error = None;
        }
        Err(error) => {
            state.entries.clear();
            state.total = 0;
            state.selected_profile = None;
            state.detail = None;
            state.name.clear();
            state.acknowledged = false;
            state.error = Some(error);
        }
    }
}

fn choose(app: &mut SolveApp, state: &mut State, id: &str) {
    state.selected_profile = Some(id.to_string());
    state.detail = None;
    state.name.clear();
    state.acknowledged = false;
    state.error = None;
    state.note = None;
    match app.run("render.library.describe", json!({"kind":"physical_profiles","id":id})) {
        Ok(detail) => {
            state.name = detail["resource"]["designation"].as_str().filter(|name| name.len() <= 128).unwrap_or(id).to_string();
            state.detail = Some(Arc::new(detail));
        }
        Err(error) => state.error = Some(error),
    }
}

fn create(app: &mut SolveApp, state: &mut State) -> Result<(), String> {
    let id = state.selected_profile.as_deref().ok_or_else(|| tr("multiphysics.materials_choose").to_string())?;
    let name = state.name.trim();
    if !state.acknowledged || name.is_empty() || name.len() > 128 || name.chars().any(char::is_control) || state.editor.editing() {
        return Err(tr("multiphysics.materials_create_requirements").into());
    }
    let plan = app.run("render.library.physical.prepare", json!({"id":id,"name":name,"acknowledge_reference":true}))?;
    // The prepared record is a read-only proposal. Never execute an arbitrary
    // command supplied by a catalog or infer an assignment from a reference.
    if plan["command"]["command"].as_str() != Some("material.create") || !plan["command"]["params"].is_object() {
        return Err(tr("multiphysics.materials_invalid_plan").into());
    }
    let reply = app.run("material.create", plan["command"]["params"].clone())?;
    let name = reply["name"].as_str().unwrap_or(name).to_string();
    state.selected_design = Some(name.clone());
    state.editor = Default::default();
    state.design_index = usize::MAX;
    // Initialise the existing editor's observation, then make an explicit name
    // selection. Its own stable-name tracking handles later insertions and Undo.
    state.editor.sync_selection(&app.session.doc, &[], &mut state.design_index);
    state.design_index = app.session.doc.material_names().iter().position(|candidate| candidate == &name).unwrap_or(usize::MAX);
    state.editor.sync_selection(&app.session.doc, &[], &mut state.design_index);
    state.note = Some(solidcraft_ui_text::format("multiphysics.materials_created", &[("name", &name)]));
    state.error = None;
    state.acknowledged = false;
    state.design = true;
    Ok(())
}

/// The workspace supplies the outer shared container. No camera or body selection is changed.
pub(crate) fn show(app: &mut SolveApp, state: &mut State, ui: &mut egui::Ui) {
    if state.uid != Some(app.session.uid) {
        *state = State { uid: Some(app.session.uid), ..Default::default() };
    }
    ui.heading(tr("multiphysics.materials_title"));
    ui.small(tr("multiphysics.materials_scope"));
    ui.horizontal_wrapped(|ui| {
        let response = crate::ui_kit::choice_chip(ui, "multiphysics_materials_references", tr("multiphysics.materials_references"), !state.design);
        crate::scenario::publish_target("multiphysics:materials_references", response.rect);
        if response.clicked() {
            state.design = false;
            state.error = None;
        }
        let response = crate::ui_kit::choice_chip(ui, "multiphysics_materials_design", tr("multiphysics.materials_design"), state.design);
        crate::scenario::publish_target("multiphysics:materials_design", response.rect);
        if response.clicked() {
            state.design = true;
            state.error = None;
        }
    });
    if let Some(error) = &state.error {
        ui.colored_label(Tokens::get().form_error, error);
    }
    if let Some(note) = &state.note {
        ui.colored_label(Tokens::get().form_success, note);
    }
    if state.design {
        ui.small(tr("multiphysics.materials_design_hint"));
        egui::ScrollArea::vertical().id_salt("multiphysics_material_design_scroll").show(ui, |ui| {
            state.editor.sync_selection(&app.session.doc, &[], &mut state.design_index);
            material_editor::rows(app, ui, &mut state.design_index, &mut state.editor);
            if !state.editor.editing() {
                let names = app.session.doc.material_names();
                state.selected_design = names.get(state.design_index).cloned();
            }
        });
    } else {
        references(app, state, ui);
    }
}

fn references(app: &mut SolveApp, state: &mut State, ui: &mut egui::Ui) {
    let response = ui_kit::row(ui, tr("multiphysics.materials_search"), |ui| {
        ui_kit::bounded_cell_editor(ui, ui.scope_id().with("multiphysics_material_query"), &mut state.query, "", false, 256)
    });
    crate::scenario::publish_target("multiphysics:materials_search", response.rect);
    if state.loaded_query.as_deref() != Some(state.query.as_str()) {
        list(app, state);
    }
    ui.small(solidcraft_ui_text::format("multiphysics.materials_count", &[("count", &state.total.to_string())]));
    let mut selected = None;
    egui::ScrollArea::vertical().id_salt("multiphysics_material_reference_list").max_height(Tokens::get().control_height * 4.0).show(ui, |ui| {
        for item in &state.entries {
            let Some(id) = item["id"].as_str() else { continue };
            let response = crate::ui_kit::list_picker(
                ui,
                &format!("multiphysics_materials_profile_{id}"),
                item["name"].as_str().unwrap_or(id),
                state.selected_profile.as_deref() == Some(id),
            );
            crate::scenario::publish_target(&format!("multiphysics:materials_profile:{id}"), response.rect);
            if response.clicked() {
                selected = Some(id.to_string());
            }
        }
    });
    if let Some(id) = selected {
        choose(app, state, &id);
    }
    if state.entries.is_empty() {
        ui.label(tr("multiphysics.materials_no_matches"));
    }
    let Some(detail) = state.detail.clone() else {
        ui.label(tr("multiphysics.materials_choose"));
        return;
    };
    ui.separator();
    ui.strong(detail["resource"]["designation"].as_str().unwrap_or(""));
    let response = ui_kit::row(ui, tr("multiphysics.materials_name"), |ui| {
        ui_kit::bounded_cell_editor(ui, ui.scope_id().with("multiphysics_material_name"), &mut state.name, "", false, 128)
    });
    crate::scenario::publish_target("multiphysics:materials_name", response.rect);
    if response.changed() {
        state.error = None;
        state.note = None;
    }
    let response = ui_kit::checkbox(ui, tr("multiphysics.materials_acknowledge"), &mut state.acknowledged);
    crate::scenario::publish_target("multiphysics:materials_acknowledge", response.rect);
    if state.editor.editing() {
        ui.small(tr("multiphysics.materials_finish_edit"));
    }
    let valid_name = !state.name.trim().is_empty() && state.name.trim().len() <= 128 && !state.name.chars().any(char::is_control);
    let response =
        ui_kit::button(ui, tr("multiphysics.materials_create"), ui_kit::Button::Enabled(state.acknowledged && valid_name && !state.editor.editing()));
    crate::scenario::publish_target("multiphysics:materials_create", response.rect);
    if response.clicked()
        && let Err(error) = create(app, state)
    {
        state.error = Some(error);
    }
    let height = ui.available_height().max(Tokens::get().control_height * 3.0);
    egui::ScrollArea::vertical().id_salt("multiphysics_material_detail_scroll").max_height(height).show(ui, |ui| details(ui, &detail));
}

fn json_line(ui: &mut egui::Ui, label: &str, value: &Value) {
    ui.add(egui::Label::new(egui::RichText::new(format!("{label}: {value}")).small().monospace()).wrap());
}

fn details(ui: &mut egui::Ui, detail: &Value) {
    let resource = &detail["resource"];
    ui.small(resource["reference_state"].as_str().unwrap_or(""));
    ui_kit::section(ui, "multiphysics_material_properties", tr("multiphysics.materials_properties"), true, |ui| {
        if let Some(properties) = resource["properties"].as_array() {
            for (index, property) in properties.iter().take(128).enumerate() {
                // Quantity identifiers, SI units and source metadata are original
                // engineering data, not translated or inferred constitutive laws.
                let title =
                    format!("{}: {} {}", property["quantity"].as_str().unwrap_or(""), property["value"], property["unit"].as_str().unwrap_or(""));
                ui_kit::section(ui, ("multiphysics_material_property", index), &title, false, |ui| {
                    json_line(ui, tr("multiphysics.materials_conditions"), &property["conditions"]);
                    json_line(ui, tr("multiphysics.materials_qualifier"), &property["qualifier"]);
                    json_line(ui, tr("multiphysics.materials_uncertainty"), &property["uncertainty"]);
                    json_line(ui, tr("multiphysics.materials_original"), &property["original"]);
                    ui.small(format!(
                        "{}: {} · {}",
                        tr("multiphysics.materials_source"),
                        property["source_id"].as_str().unwrap_or(""),
                        property["source_locator"].as_str().unwrap_or("")
                    ));
                    ui.small(property["note"].as_str().unwrap_or(""));
                });
            }
        }
    });
    ui_kit::section(ui, "multiphysics_material_sources", tr("multiphysics.materials_sources"), true, |ui| {
        if let Some(sources) = detail["sources"].as_array() {
            for source in sources.iter().take(128) {
                let title = source["title"].as_str().unwrap_or("");
                if let Some(url) = source["url"].as_str() {
                    ui.hyperlink_to(title, url);
                } else {
                    ui.label(title);
                }
                json_line(ui, tr("multiphysics.materials_source"), source);
            }
        }
    });
    ui_kit::section(ui, "multiphysics_material_omissions", tr("multiphysics.materials_omissions"), true, |ui| {
        if let Some(items) = detail["native_projection"]["not_projected"].as_array() {
            for item in items.iter().take(128) {
                ui.small(format!("{}: {}", item["quantity"].as_str().unwrap_or(""), item["reason"].as_str().unwrap_or("")));
            }
        }
        if let Some(items) = resource["limitations"].as_array() {
            for item in items.iter().take(128) {
                ui.small(item.as_str().unwrap_or(""));
            }
        }
    });
    ui_kit::section(ui, "multiphysics_material_native_card", tr("multiphysics.materials_native_card"), false, |ui| {
        json_line(ui, tr("multiphysics.materials_native_card"), &detail["native_projection"]["card"]);
        json_line(ui, tr("multiphysics.materials_receipts"), &detail["receipts"]);
    });
}
