//! The Parts folder of a study of the model: its parts in the Analysis tree and in the component
//! table, both showing the design's own selection, so a part picked in the tree, the table or the
//! 3-D view is picked in all three. Hide, Isolate and Show All are the design's
//! `browser.visibility` (saved, undone); Ghost draws parts see-through (`multiphysics.parts.ghost`,
//! view state) so the parts inside an assembly show and can be picked.
use egui::RichText;
use serde_json::{Value, json};
use solidcraft_theme::Tokens;

use super::{State, command, tr};
use crate::SolveApp;
use solidcraft_engine::multiphysics::Study;

/// What the panel remembers between frames.
#[derive(Default)]
pub(super) struct PartsUi {
    /// The last `multiphysics.parts.table` reply: for the document (its address) and study.
    table: Option<(usize, String, Result<Value, String>)>,
    /// The part selection Select picks from the list.
    named: String,
}

/// Something done to parts from the tree, the table or a menu.
#[derive(Clone, Debug, PartialEq)]
pub(super) enum PartAction {
    /// Pick a part; `add`: toggle it in the selection (Ctrl or Cmd held).
    Pick {
        body: String,
        add: bool,
    },
    /// The bodies of a named part selection.
    PickNamed(String),
    Hide(Vec<String>),
    Isolate(Vec<String>),
    /// See-through (`true`) or solid again.
    Ghost(Vec<String>, bool),
    ShowAll,
    SolidAll,
}

/// The study's component table, worked out again when the design or the study changes.
fn table(app: &mut SolveApp, state: &mut State, study: &Study) -> Result<Value, String> {
    let key = std::sync::Arc::as_ptr(&app.session.doc) as usize;
    if let Some((k, id, reply)) = &state.parts_ui.table
        && *k == key
        && *id == study.id
    {
        return reply.clone();
    }
    let reply = app.run("multiphysics.parts.table", json!({"id": study.id}));
    state.parts_ui.table = Some((key, study.id.clone(), reply.clone()));
    reply
}

/// The study's parts, in the design's order.
pub(super) fn bodies(app: &mut SolveApp, state: &mut State, study: &Study) -> Vec<String> {
    match table(app, state, study) {
        Ok(reply) => reply["table"]["rows"].as_array().into_iter().flatten().filter_map(|r| r["body"].as_str().map(str::to_string)).collect(),
        Err(_) => Vec::new(),
    }
}

/// The bodies in the design's selection.
pub(super) fn selected(app: &SolveApp) -> Vec<String> {
    app.session
        .selection
        .iter()
        .filter_map(|sel| match sel {
            solidcraft_engine::Sel::Body { name } => Some(name.clone()),
            _ => None,
        })
        .collect()
}

/// Is `body` drawn see-through?
pub(super) fn ghosted(app: &SolveApp, body: &str) -> bool {
    app.session.ghosted_bodies.iter().any(|g| g == body)
}

/// A part's name as its rows show it: with "hidden" or "ghosted" after it when it is.
pub(super) fn label(app: &SolveApp, body: &str) -> String {
    let key = match (app.body_visible(body), ghosted(app, body)) {
        (false, _) => "multiphysics.parts_hidden",
        (true, true) => "multiphysics.parts_ghosted",
        (true, false) => return body.to_string(),
    };
    solidcraft_ui_text::format(key, &[("part", body)])
}

/// What a click on a part's row does: with Ctrl or Cmd it toggles the part in the selection.
pub(super) fn pick(ui: &egui::Ui, body: &str) -> PartAction {
    PartAction::Pick { body: body.to_string(), add: ui.input(|i| i.modifiers.command) }
}

/// The part's context menu: on the selection when the part is in it, otherwise on the part.
pub(super) fn menu(app: &SolveApp, ui: &mut egui::Ui, body: &str, actions: &mut Vec<PartAction>) {
    let chosen = selected(app);
    let bodies = if chosen.iter().any(|b| b == body) { chosen } else { vec![body.to_string()] };
    let item = |ui: &mut egui::Ui, key: &str| {
        let picked = crate::ui_kit::button(ui, tr(key), crate::ui_kit::Button::Ghost).clicked();
        if picked {
            ui.close();
        }
        picked
    };
    if item(ui, "multiphysics.parts_select") {
        actions.push(PartAction::Pick { body: body.to_string(), add: false });
    }
    if item(ui, "multiphysics.parts_hide") {
        actions.push(PartAction::Hide(bodies.clone()));
    }
    if item(ui, "multiphysics.parts_isolate") {
        actions.push(PartAction::Isolate(bodies.clone()));
    }
    let all_ghosted = bodies.iter().all(|b| ghosted(app, b));
    if item(ui, if all_ghosted { "multiphysics.parts_solid" } else { "multiphysics.parts_ghost" }) {
        actions.push(PartAction::Ghost(bodies, !all_ghosted));
    }
    ui.separator();
    if item(ui, "multiphysics.parts_show_all") {
        actions.push(PartAction::ShowAll);
    }
    if item(ui, "multiphysics.parts_solid_all") {
        actions.push(PartAction::SolidAll);
    }
}

/// Do `action` through the commands, so it is the same from the tree, the table and a script.
pub(super) fn apply(app: &mut SolveApp, state: &mut State, action: PartAction) {
    let (id, params) = match action {
        PartAction::Pick { body, add: false } => ("multiphysics.parts.select", json!({"bodies": [body]})),
        PartAction::Pick { body, add: true } => {
            let mut bodies = selected(app);
            match bodies.iter().position(|b| *b == body) {
                Some(at) => {
                    bodies.remove(at);
                }
                None => bodies.push(body),
            }
            ("multiphysics.parts.select", json!({"bodies": bodies}))
        }
        PartAction::PickNamed(name) => ("multiphysics.parts.select", json!({"selection": name})),
        PartAction::Hide(bodies) => ("browser.visibility", json!({"items": bodies, "visible": false})),
        PartAction::Isolate(bodies) => ("browser.visibility", json!({"items": bodies, "visible": true, "isolate": true})),
        PartAction::Ghost(bodies, ghost) => ("multiphysics.parts.ghost", json!({"bodies": bodies, "ghost": ghost})),
        PartAction::ShowAll => ("browser.visibility", json!({"folder": "bodies", "visible": true})),
        PartAction::SolidAll => ("multiphysics.parts.ghost", json!({"only": true})),
    };
    let _ = command(app, state, id, params);
}

/// Temperature in °C, from kelvin.
fn celsius(k: Option<f64>) -> String {
    k.map_or_else(|| "–".into(), |k| format!("{:.1}", k - 273.15))
}

/// The Parts node's panel: the component table, its rows the design's selection, and what the
/// selection can be shown as.
pub(super) fn panel(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let t = Tokens::get();
    ui.label(RichText::new(tr("multiphysics.parts_hint")).small().color(t.text_dim));
    let reply = match table(app, state, study) {
        Ok(reply) => reply,
        Err(error) => {
            ui.colored_label(t.form_error, error);
            return;
        }
    };
    let chosen = selected(app);
    let mut actions = Vec::new();
    // Select a named part selection.
    let names: Vec<String> = app.session.doc.multiphysics.parts.selections.iter().map(|s| s.name.clone()).collect();
    if !names.is_empty() {
        if !names.contains(&state.parts_ui.named) {
            state.parts_ui.named = names[0].clone();
        }
        crate::ui_kit::row(ui, tr("multiphysics.parts_named"), |ui| {
            let current = state.parts_ui.named.clone();
            crate::ui_kit::dropdown_menu(ui, "multiphysics_parts_named", &current, |ui| {
                for name in &names {
                    let response = crate::ui_kit::list_picker(ui, &format!("multiphysics_parts_named_{name}"), name, *name == current);
                    crate::scenario::publish_target(&format!("multiphysics:parts_named:{name}"), response.rect);
                    if response.clicked() {
                        state.parts_ui.named = name.clone();
                        ui.close();
                    }
                }
            });
        });
        // Its own line: beside the list, which takes the row's width, it would push the panel wide.
        let response = crate::ui_kit::button(ui, tr("multiphysics.parts_select"), crate::ui_kit::Button::Small);
        crate::scenario::publish_target("multiphysics:parts_select_named", response.rect);
        if response.clicked() {
            actions.push(PartAction::PickNamed(state.parts_ui.named.clone()));
        }
    }
    // What the selection is shown as.
    ui.horizontal_wrapped(|ui| {
        let some = !chosen.is_empty();
        let all_ghosted = some && chosen.iter().all(|b| ghosted(app, b));
        for (key, handle, enabled, action) in [
            ("multiphysics.parts_hide", "multiphysics:parts_hide", some, PartAction::Hide(chosen.clone())),
            ("multiphysics.parts_isolate", "multiphysics:parts_isolate", some, PartAction::Isolate(chosen.clone())),
            (
                if all_ghosted { "multiphysics.parts_solid" } else { "multiphysics.parts_ghost" },
                "multiphysics:parts_ghost",
                some,
                PartAction::Ghost(chosen.clone(), !all_ghosted),
            ),
            ("multiphysics.parts_show_all", "multiphysics:parts_show_all", true, PartAction::ShowAll),
            ("multiphysics.parts_solid_all", "multiphysics:parts_solid_all", !app.session.ghosted_bodies.is_empty(), PartAction::SolidAll),
        ] {
            let response = crate::ui_kit::button(ui, tr(key), crate::ui_kit::Button::Enabled(enabled));
            crate::scenario::publish_target(handle, response.rect);
            if response.clicked() {
                actions.push(action);
            }
        }
    });
    ui.separator();
    let table = &reply["table"];
    if table["run"].is_null() {
        ui.label(RichText::new(tr("multiphysics.parts_no_run")).small().color(t.text_dim));
    } else if table["stale"] == true {
        ui.colored_label(t.warning, tr("multiphysics.parts_stale"));
    }
    let rows = table["rows"].as_array().cloned().unwrap_or_default();
    egui::ScrollArea::vertical().id_salt("multiphysics_parts_table").max_height(360.0).show(ui, |ui| {
        egui::Grid::new("multiphysics_parts_grid").striped(true).show(ui, |ui| {
            for key in [
                "multiphysics.parts_column_part",
                "multiphysics.parts_column_power",
                "multiphysics.parts_column_limit",
                "multiphysics.parts_column_max",
                "multiphysics.parts_column_margin",
                "multiphysics.parts_column_status",
            ] {
                ui.strong(tr(key));
            }
            ui.end_row();
            for row in rows.iter().take(super::MAX_TABLE_ROWS) {
                let Some(body) = row["body"].as_str() else { continue };
                let on = chosen.iter().any(|b| b == body);
                let label = label(app, body);
                let response = crate::ui_kit::choice_chip(ui, &format!("multiphysics_part_{body}"), &label, on);
                crate::scenario::publish_target(&format!("multiphysics:parts_row:{body}"), response.rect);
                if response.clicked() {
                    actions.push(pick(ui, body));
                }
                crate::ui_kit::Menu::context(&response, |ui| menu(app, ui, body, &mut actions));
                ui.monospace(format!("{:.2}", row["power_w"].as_f64().unwrap_or(0.0)));
                ui.monospace(celsius(row["limit_k"].as_f64()));
                ui.monospace(celsius(row["result"]["max_temperature_k"].as_f64()));
                ui.monospace(row["margin_k"].as_f64().map_or_else(|| "–".into(), |m| format!("{m:.1}")));
                // No verdict until a run reports the part.
                let (key, color) = match (row["result"].is_null(), row["status"].as_str()) {
                    (true, _) => ("multiphysics.parts_no_result", t.text_dim),
                    (false, Some("pass")) => ("multiphysics.post_status_pass", t.form_success),
                    (false, Some("warn")) => ("multiphysics.post_status_warn", t.warning),
                    (false, Some("fail")) => ("multiphysics.post_status_fail", t.form_error),
                    _ => ("multiphysics.post_status_no_limit", t.text_dim),
                };
                ui.colored_label(color, tr(key));
                ui.end_row();
            }
        });
    });
    let total = table["total_power_w"].as_f64().unwrap_or(0.0);
    ui.label(RichText::new(solidcraft_ui_text::format("multiphysics.parts_total_power", &[("watts", &format!("{total:.2}"))])).small());
    for action in actions {
        apply(app, state, action);
    }
}
