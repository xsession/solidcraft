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
    /// The reference family chosen in the database tree; `None` lists every family.
    family: Option<String>,
    /// The right pane's page: the item list, or the chosen item's properties.
    properties: bool,
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
            family: None,
            properties: false,
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

/// Reference families in the order the database tree lists them; others follow by name.
const FAMILIES: [&str; 5] = ["metal", "polymer", "glass", "wood", "other"];

fn family_label(family: &str) -> String {
    match family {
        "metal" => tr("multiphysics.materials_family_metal").to_string(),
        "polymer" => tr("multiphysics.materials_family_polymer").to_string(),
        "glass" => tr("multiphysics.materials_family_glass").to_string(),
        "wood" => tr("multiphysics.materials_family_wood").to_string(),
        "other" => tr("multiphysics.materials_family_other").to_string(),
        _ => family.replace('_', " "),
    }
}

/// The families present in the listed entries with their counts, known families first.
fn families(entries: &[Value]) -> Vec<(String, usize)> {
    let mut found: Vec<(String, usize)> = Vec::new();
    for item in entries {
        let family = item["family"].as_str().unwrap_or("other");
        match found.iter_mut().find(|(name, _)| name == family) {
            Some((_, count)) => *count += 1,
            None => found.push((family.to_string(), 1)),
        }
    }
    found.sort_by_key(|(name, _)| (FAMILIES.iter().position(|known| known == name).unwrap_or(FAMILIES.len()), name.clone()));
    found
}

/// The workspace supplies the outer shared container. No camera or body selection is changed.
/// Laid out as an engineering database: the database tree on the left, and on the right the
/// items of the chosen folder or the properties of the chosen item.
pub(crate) fn show(app: &mut SolveApp, state: &mut State, ui: &mut egui::Ui) {
    if state.uid != Some(app.session.uid) {
        *state = State { uid: Some(app.session.uid), ..Default::default() };
    }
    if !state.design && state.loaded_query.as_deref() != Some(state.query.as_str()) {
        list(app, state);
    }
    let t = Tokens::get();
    let outer = ui.max_rect().shrink(t.space_8);
    ui.scope_builder(egui::UiBuilder::new().max_rect(outer).layout(egui::Layout::top_down(egui::Align::Min)), |ui| {
        ui.horizontal(|ui| {
            ui.heading(tr("multiphysics.materials_title"));
            ui.add_space(t.space_12);
            ui.small(tr("multiphysics.materials_scope"));
        });
        if let Some(error) = &state.error {
            ui.colored_label(t.form_error, error);
        }
        if let Some(note) = &state.note {
            ui.colored_label(t.form_success, note);
        }
        ui.add_space(t.space_4);
        // The footer that copies a chosen reference into the design stays under both pages.
        let footer = if !state.design && state.detail.is_some() { t.control_height * 2.0 + t.space_12 } else { 0.0 };
        let body = ui.available_rect_before_wrap();
        let body = egui::Rect::from_min_max(body.min, egui::pos2(body.right(), (body.bottom() - footer).max(body.top())));
        let tree_width = (body.width() * 0.26).clamp(200.0, 300.0);
        let tree = egui::Rect::from_min_max(body.min, egui::pos2(body.left() + tree_width, body.bottom()));
        let pane = egui::Rect::from_min_max(egui::pos2(tree.right() + t.space_8, body.top()), body.max);
        let response = ui.interact(body, ui.scope_id().with("multiphysics_materials_body"), egui::Sense::hover());
        crate::scenario::publish_rect("multiphysics:materials_body", response.rect);
        database_tree(app, state, ui, tree);
        ui.scope_builder(egui::UiBuilder::new().max_rect(pane).layout(egui::Layout::top_down(egui::Align::Min)), |ui| {
            ui.set_clip_rect(pane.intersect(ui.clip_rect()));
            if state.design {
                crate::multiphysics_property_grid::band(ui, tr("multiphysics.materials_design"));
                ui.small(tr("multiphysics.materials_design_hint"));
                // While a material is being edited the editor gets the whole pane.
                if !state.editor.editing() {
                    design_items(app, state, ui);
                    ui.add_space(t.space_8);
                }
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
        });
        if footer > 0.0 {
            let rect = egui::Rect::from_min_max(egui::pos2(outer.left(), body.bottom() + t.space_8), outer.max);
            ui.scope_builder(egui::UiBuilder::new().max_rect(rect).layout(egui::Layout::top_down(egui::Align::Min)), |ui| copy_bar(app, state, ui));
        }
    });
}

/// The database tree: the reference library by family, and the design's own materials.
fn database_tree(app: &SolveApp, state: &mut State, ui: &mut egui::Ui, rect: egui::Rect) {
    let grid = crate::multiphysics_property_grid::tree_row;
    crate::multiphysics_property_grid::pane(ui, rect, tr("multiphysics.materials_tree"), |ui| {
        egui::ScrollArea::vertical().id_salt("multiphysics_material_tree").show(ui, |ui| {
            let all = solidcraft_ui_text::format(
                "multiphysics.materials_tree_count",
                &[("name", tr("multiphysics.materials_references")), ("count", &state.total.to_string())],
            );
            let response = grid(ui, 0, Some(true), "material", &all, !state.design && state.family.is_none(), true, false);
            crate::scenario::publish_target("multiphysics:materials_references", response.rect);
            if response.clicked() {
                state.design = false;
                state.family = None;
                state.properties = false;
                state.error = None;
            }
            for (family, count) in families(&state.entries) {
                let text = solidcraft_ui_text::format(
                    "multiphysics.materials_tree_count",
                    &[("name", &family_label(&family)), ("count", &count.to_string())],
                );
                let selected = !state.design && state.family.as_deref() == Some(family.as_str());
                let response = grid(ui, 1, None, "material", &text, selected, true, false);
                crate::scenario::publish_target(&format!("multiphysics:materials_family:{family}"), response.rect);
                if response.clicked() {
                    state.design = false;
                    state.family = Some(family);
                    state.properties = false;
                    state.error = None;
                }
            }
            let names = app.session.doc.material_names().len();
            let text = solidcraft_ui_text::format(
                "multiphysics.materials_tree_count",
                &[("name", tr("multiphysics.materials_design")), ("count", &names.to_string())],
            );
            let response = grid(ui, 0, None, "material", &text, state.design, true, false);
            crate::scenario::publish_target("multiphysics:materials_design", response.rect);
            if response.clicked() {
                state.design = true;
                state.error = None;
            }
        });
    });
}

/// The right pane for references: a search row, the Items | Item properties pages.
fn references(app: &mut SolveApp, state: &mut State, ui: &mut egui::Ui) {
    let t = Tokens::get();
    ui.horizontal(|ui| {
        let response = crate::ui_kit::choice_chip(ui, "multiphysics_materials_items", tr("multiphysics.materials_items"), !state.properties);
        crate::scenario::publish_target("multiphysics:materials_items", response.rect);
        if response.clicked() {
            state.properties = false;
        }
        let response =
            crate::ui_kit::choice_chip(ui, "multiphysics_materials_properties", tr("multiphysics.materials_item_properties"), state.properties);
        crate::scenario::publish_target("multiphysics:materials_properties", response.rect);
        if response.clicked() && state.detail.is_some() {
            state.properties = true;
        }
    });
    ui.add_space(t.space_4);
    if state.properties
        && let Some(detail) = state.detail.clone()
    {
        egui::ScrollArea::vertical().id_salt("multiphysics_material_detail_scroll").show(ui, |ui| details(ui, &detail));
        return;
    }
    state.properties = false;
    let response = ui_kit::row(ui, tr("multiphysics.materials_search"), |ui| {
        ui_kit::bounded_cell_editor(ui, ui.scope_id().with("multiphysics_material_query"), &mut state.query, "", false, 256)
    });
    crate::scenario::publish_target("multiphysics:materials_search", response.rect);
    if state.loaded_query.as_deref() != Some(state.query.as_str()) {
        list(app, state);
    }
    let shares = [0.62, 0.2, 0.18];
    crate::multiphysics_property_grid::header(
        ui,
        &[tr("multiphysics.materials_col_name"), tr("multiphysics.materials_col_family"), tr("multiphysics.materials_col_properties")],
        &shares,
    );
    let mut selected = None;
    let mut open = false;
    egui::ScrollArea::vertical().id_salt("multiphysics_material_reference_list").auto_shrink([false, true]).show(ui, |ui| {
        let mut shown = 0usize;
        for item in &state.entries {
            let Some(id) = item["id"].as_str() else { continue };
            let family = item["family"].as_str().unwrap_or("other");
            if state.family.as_deref().is_some_and(|chosen| chosen != family) {
                continue;
            }
            let count = item["property_count"].as_u64().map(|n| n.to_string()).unwrap_or_default();
            let response = crate::multiphysics_property_grid::item(
                ui,
                &[item["name"].as_str().unwrap_or(id), &family_label(family), &count],
                &shares,
                state.selected_profile.as_deref() == Some(id),
                !shown.is_multiple_of(2),
            );
            shown += 1;
            if ui.clip_rect().contains_rect(response.rect) {
                crate::scenario::publish_target(&format!("multiphysics:materials_profile:{id}"), response.rect);
            }
            if response.clicked() || response.double_clicked() {
                selected = Some(id.to_string());
                open = response.double_clicked();
            }
        }
        if shown == 0 {
            ui.add_space(t.space_8);
            ui.label(tr("multiphysics.materials_no_matches"));
        }
    });
    if let Some(id) = selected {
        if state.selected_profile.as_deref() != Some(id.as_str()) {
            choose(app, state, &id);
        }
        state.properties = open && state.detail.is_some();
    }
    if state.detail.is_none() {
        ui.small(tr("multiphysics.materials_choose"));
    }
}

/// The design's materials as a table, built-in first: choosing a row opens it in the editor below.
fn design_items(app: &SolveApp, state: &mut State, ui: &mut egui::Ui) {
    let doc = &app.session.doc;
    let names = doc.material_names();
    let shares = [0.5, 0.25, 0.25];
    crate::multiphysics_property_grid::header(
        ui,
        &[tr("multiphysics.materials_col_name"), tr("multiphysics.materials_col_origin"), tr("multiphysics.materials_col_density")],
        &shares,
    );
    let height = (ui.available_height() * 0.35).clamp(crate::ui_kit::row_height() * 3.0, crate::ui_kit::row_height() * 8.0);
    egui::ScrollArea::vertical().id_salt("multiphysics_material_design_items").max_height(height).auto_shrink([false, true]).show(ui, |ui| {
        for (index, name) in names.iter().enumerate().take(512) {
            let origin = if doc.user_material(name).is_some() {
                tr("multiphysics.materials_origin_design")
            } else {
                tr("multiphysics.materials_origin_builtin")
            };
            let density = doc.material_density(name).map(|d| format!("{d:.3}")).unwrap_or_default();
            let response = crate::multiphysics_property_grid::item(
                ui,
                &[name, origin, &density],
                &shares,
                state.design_index == index,
                !index.is_multiple_of(2),
            );
            if ui.clip_rect().contains_rect(response.rect) {
                crate::scenario::publish_target(&format!("multiphysics:materials_design_item:{name}"), response.rect);
            }
            if response.clicked() {
                state.design_index = index;
            }
        }
    });
}

/// Copy the chosen reference into the design: its name there, the review, the action.
fn copy_bar(app: &mut SolveApp, state: &mut State, ui: &mut egui::Ui) {
    let Some(detail) = state.detail.clone() else { return };
    let t = Tokens::get();
    ui.painter().hline(ui.max_rect().x_range(), ui.max_rect().top(), egui::Stroke::new(1.0, t.dialog_rule));
    ui.add_space(t.space_4);
    ui.horizontal(|ui| {
        ui.strong(detail["resource"]["designation"].as_str().unwrap_or(""));
        if state.editor.editing() {
            ui.small(tr("multiphysics.materials_finish_edit"));
        }
    });
    ui.horizontal(|ui| {
        ui.label(tr("multiphysics.materials_name"));
        let response = ui
            .scope(|ui| {
                ui.set_width((ui.available_width() * 0.4).clamp(160.0, 360.0));
                ui_kit::bounded_cell_editor(ui, ui.scope_id().with("multiphysics_material_name"), &mut state.name, "", false, 128)
            })
            .inner;
        crate::scenario::publish_target("multiphysics:materials_name", response.rect);
        if response.changed() {
            state.error = None;
            state.note = None;
        }
        let response = ui_kit::checkbox(ui, tr("multiphysics.materials_acknowledge"), &mut state.acknowledged);
        crate::scenario::publish_target("multiphysics:materials_acknowledge", response.rect);
        let valid_name = !state.name.trim().is_empty() && state.name.trim().len() <= 128 && !state.name.chars().any(char::is_control);
        let response = ui_kit::button(
            ui,
            tr("multiphysics.materials_create"),
            ui_kit::Button::Enabled(state.acknowledged && valid_name && !state.editor.editing()),
        );
        crate::scenario::publish_target("multiphysics:materials_create", response.rect);
        if response.clicked()
            && let Err(error) = create(app, state)
        {
            state.error = Some(error);
        }
    });
}

fn json_line(ui: &mut egui::Ui, label: &str, value: &Value) {
    ui.add(egui::Label::new(egui::RichText::new(format!("{label}: {value}")).small().monospace()).wrap());
}

fn details(ui: &mut egui::Ui, detail: &Value) {
    let resource = &detail["resource"];
    ui.strong(resource["designation"].as_str().unwrap_or(""));
    ui.small(resource["reference_state"].as_str().unwrap_or(""));
    let grid = crate::multiphysics_property_grid::item;
    let shares = [0.34, 0.24, 0.16, 0.26];
    crate::multiphysics_property_grid::header(
        ui,
        &[
            tr("multiphysics.materials_col_property"),
            tr("multiphysics.materials_col_value"),
            tr("multiphysics.materials_col_unit"),
            tr("multiphysics.materials_qualifier"),
        ],
        &shares,
    );
    if let Some(properties) = resource["properties"].as_array() {
        for (index, property) in properties.iter().take(128).enumerate() {
            // Quantity identifiers, SI units and source metadata are original
            // engineering data, not translated or inferred constitutive laws.
            let quantity = property["quantity"].as_str().unwrap_or("").replace('_', " ");
            let value = match &property["value"] {
                Value::Number(n) => n.as_f64().map(|v| format!("{v}")).unwrap_or_default(),
                Value::Null => tr("multiphysics.not_specified").to_string(),
                other => other.to_string(),
            };
            let row = ui.push_id(("multiphysics_material_property", index), |ui| {
                grid(
                    ui,
                    &[&quantity, &value, property["unit"].as_str().unwrap_or(""), property["qualifier"].as_str().unwrap_or("")],
                    &shares,
                    false,
                    !index.is_multiple_of(2),
                )
            });
            row.inner.on_hover_ui(|ui| {
                json_line(ui, tr("multiphysics.materials_conditions"), &property["conditions"]);
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
    ui.add_space(Tokens::get().space_8);
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
