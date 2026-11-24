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
use solidcraft_engine::multiphysics::setup_check::{ColourBy, Swatch};

/// What the panel remembers between frames.
#[derive(Default)]
pub(super) struct PartsUi {
    /// The last `multiphysics.parts.table` reply: for the document (its address) and study.
    table: Option<(usize, String, Result<Value, String>)>,
    /// The part selection Select picks from the list.
    named: String,
    /// What the parts are coloured by, and the last colouring: for the document, study and mode.
    colour_by: ColourBy,
    colouring: Option<(usize, String, ColourBy, Result<Value, String>)>,
    /// The last checklist: for the document and study.
    checks: Option<(usize, String, Result<Value, String>)>,
    /// The Find section's stacked rules, and its last answer.
    find: Vec<FindRow>,
    found: Option<Result<Value, String>>,
    /// The example rule was put in once; removing every rule leaves none.
    find_seeded: bool,
    /// The Group section: its last answer, and the size under which parts are left out, as typed.
    grouped: Option<Result<String, String>>,
    small_mm: String,
    /// The batch report's new directory, as typed, and its last answer (native only).
    #[cfg(not(target_arch = "wasm32"))]
    report_dir: String,
    #[cfg(not(target_arch = "wasm32"))]
    reported: Option<Result<String, String>>,
    /// Compare runs: the two runs picked (study id, run id), and the last comparison of them.
    pub(super) compare: [Option<(String, String)>; 2],
    pub(super) compared: Option<Result<Value, String>>,
    /// The mesh study's last answer: making its grids, or their convergence.
    pub(super) mesh_study: Option<Result<Value, String>>,
}

/// One rule of the Find section: a field of `FIND_FIELDS`, a comparison of `COMPARES` for
/// numbers, and the value as typed (for a verdict, its index in `VERDICTS`).
#[derive(Clone, Debug, PartialEq)]
pub(super) struct FindRow {
    field: usize,
    compare: usize,
    value: String,
}

impl Default for FindRow {
    /// The Find panel's own example: hotter than 75 °C.
    fn default() -> Self {
        Self { field: 5, compare: 2, value: "75".into() }
    }
}

/// The fields a Find rule reads: `(wire name, label key, takes a number)`.
const FIND_FIELDS: &[(&str, &str, bool)] = &[
    ("tag", "multiphysics.find_tag", false),
    ("not_tag", "multiphysics.find_not_tag", false),
    ("name", "multiphysics.find_name", false),
    ("power_w", "multiphysics.find_power", true),
    ("limit_c", "multiphysics.find_limit", true),
    ("max_temperature_c", "multiphysics.find_max", true),
    ("mean_temperature_c", "multiphysics.find_mean", true),
    ("margin_k", "multiphysics.find_margin", true),
    ("verdict", "multiphysics.find_verdict", false),
];
/// Comparisons for numbers: `(wire name, shown as)`.
const COMPARES: &[(&str, &str)] = &[("lt", "<"), ("le", "≤"), ("gt", ">"), ("ge", "≥"), ("eq", "=")];
/// Verdicts a rule can ask for: `(wire name, label key)`.
const VERDICTS: &[(&str, &str)] = &[("fail", "multiphysics.colour_fail"), ("warn", "multiphysics.colour_warn"), ("pass", "multiphysics.colour_pass")];

impl FindRow {
    /// The rule as `multiphysics.parts.find` takes it; `None` while it is incomplete.
    fn rule(&self) -> Option<Value> {
        let (field, _, numeric) = FIND_FIELDS.get(self.field)?;
        let text = self.value.trim();
        Some(match *field {
            "tag" | "not_tag" if !text.is_empty() => json!({"rule": field, "tag": text}),
            "name" if !text.is_empty() => json!({"rule": "name", "pattern": text}),
            "verdict" => json!({"rule": "verdict", "status": VERDICTS.get(text.parse::<usize>().unwrap_or(0))?.0}),
            quantity if *numeric => {
                let value = text.replace(',', ".").parse::<f64>().ok().filter(|v| v.is_finite())?;
                json!({"rule": "result", "quantity": quantity, "compare": COMPARES.get(self.compare)?.0, "value": value})
            }
            _ => return None,
        })
    }
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
    /// These bodies (a checklist line's).
    PickMany(Vec<String>),
    /// Colour the parts by their setup.
    ColourBy(ColourBy),
    /// A legend entry's parts on their own: `keep` selected and solid, `others` ghosted; or,
    /// when exactly `others` are ghosted already, every part solid again.
    Focus {
        keep: Vec<String>,
        others: Vec<String>,
    },
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
    app.session.ghosted_bodies.contains(body)
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
        PartAction::PickMany(bodies) => ("multiphysics.parts.select", json!({"bodies": bodies})),
        PartAction::Focus { keep, others } => {
            let mut ghosted: Vec<&String> = app.session.ghosted_bodies.iter().collect();
            let mut wanted: Vec<&String> = others.iter().collect();
            ghosted.sort();
            wanted.sort();
            if !others.is_empty() && ghosted == wanted {
                let _ = command(app, state, "multiphysics.parts.ghost", json!({"only": true}));
                return;
            }
            let _ = command(app, state, "multiphysics.parts.select", json!({"bodies": keep}));
            ("multiphysics.parts.ghost", json!({"bodies": others, "only": true}))
        }
        PartAction::ColourBy(by) => {
            state.parts_ui.colour_by = by;
            state.parts_ui.colouring = None;
            if by == ColourBy::None {
                let _ = command(app, state, "multiphysics.parts.colour", json!({"by": "none"}));
            }
            return;
        }
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
    group_section(app, state, study, ui);
    ui.separator();
    find_section(app, state, study, ui, &mut actions);
    ui.separator();
    colour_section(app, state, study, ui, &mut actions);
    ui.separator();
    checklist_section(app, state, study, ui, &mut actions);
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
    #[cfg(not(target_arch = "wasm32"))]
    report_section(app, state, study, ui);
    compare_section(app, state, study, ui);
    mesh_study_section(app, state, study, ui);
    for action in actions {
        apply(app, state, action);
    }
}

/// The colour-by modes in the order the panel offers them, with their wire names.
const COLOUR_BY: &[(ColourBy, &str)] =
    &[(ColourBy::None, "none"), (ColourBy::Material, "material"), (ColourBy::Power, "power"), (ColourBy::Tag, "tag"), (ColourBy::Verdict, "verdict")];

/// Keep the parts coloured as picked: worked out again when the design or the study changes,
/// or when leaving the workspace dropped the colours. Called each frame the Parts folder shows.
pub(super) fn keep_colour(app: &mut SolveApp, state: &mut State, study: &Study) -> Option<Value> {
    let by = state.parts_ui.colour_by;
    if by == ColourBy::None {
        return None;
    }
    let key = std::sync::Arc::as_ptr(&app.session.doc) as usize;
    let dropped = app.session.part_swatches.is_empty();
    if let Some((k, id, was, reply)) = &state.parts_ui.colouring
        && *k == key
        && *id == study.id
        && *was == by
        && !(dropped && reply.as_ref().is_ok_and(|r| r["parts"].as_array().is_some_and(|p| !p.is_empty())))
    {
        return reply.as_ref().ok().cloned();
    }
    let name = COLOUR_BY.iter().find(|(b, _)| *b == by).map_or("none", |(_, n)| *n);
    let reply = app.run("multiphysics.parts.colour", json!({"id": study.id, "by": name}));
    state.parts_ui.colouring = Some((key, study.id.clone(), by, reply.clone()));
    reply.ok()
}

/// The study's checklist, worked out again when the design or the study changes.
fn checks(app: &mut SolveApp, state: &mut State, study: &Study) -> Result<Value, String> {
    let key = std::sync::Arc::as_ptr(&app.session.doc) as usize;
    if let Some((k, id, reply)) = &state.parts_ui.checks
        && *k == key
        && *id == study.id
    {
        return reply.clone();
    }
    let reply = app.run("multiphysics.parts.check", json!({"id": study.id}));
    state.parts_ui.checks = Some((key, study.id.clone(), reply.clone()));
    reply
}

/// A small square in a swatch's colour, as the 3-D view draws it.
fn swatch_square(ui: &mut egui::Ui, swatch: Swatch) {
    let t = Tokens::get();
    let (rect, _) = ui.allocate_exact_size(egui::vec2(t.space_14, t.space_14), egui::Sense::hover());
    ui.painter().rect_filled(rect, t.control_radius, crate::viewport::swatch_colour(swatch));
}

/// What a legend entry's label says: a known id worded, otherwise the label itself (a tag).
fn legend_label(label: &str) -> String {
    let key = format!("multiphysics.colour_{label}");
    let worded = tr(&key);
    if worded == key { label.to_string() } else { worded.to_string() }
}

/// Colour-by and its legend.
pub(super) fn colour_section(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui, actions: &mut Vec<PartAction>) {
    let t = Tokens::get();
    ui.strong(tr("multiphysics.parts_colour_by"));
    ui.horizontal_wrapped(|ui| {
        for (by, name) in COLOUR_BY {
            let on = state.parts_ui.colour_by == *by;
            let response =
                crate::ui_kit::choice_chip(ui, &format!("multiphysics_parts_colour_{name}"), tr(&format!("multiphysics.colour_by_{name}")), on);
            crate::scenario::publish_target(&format!("multiphysics:parts_colour:{name}"), response.rect);
            if response.clicked() && !on {
                actions.push(PartAction::ColourBy(*by));
            }
        }
    });
    let Some(colouring) = keep_colour(app, state, study) else {
        if let Some((_, _, _, Err(error))) = &state.parts_ui.colouring {
            ui.colored_label(t.form_error, error);
        }
        return;
    };
    if let Some(max) = colouring["max"].as_f64() {
        // The power map, low to high.
        ui.horizontal(|ui| {
            ui.label(RichText::new("0 W").small());
            for step in 0..=8 {
                swatch_square(ui, Swatch::Scale { t: f64::from(step) / 8.0 });
            }
            ui.label(RichText::new(format!("{max:.2} W")).small());
        });
    }
    let parts: Vec<(String, String)> = colouring["parts"]
        .as_array()
        .into_iter()
        .flatten()
        .filter_map(|p| Some((p["body"].as_str()?.to_string(), p["label"].as_str()?.to_string())))
        .collect();
    let legend = colouring["legend"].as_array().cloned().unwrap_or_default();
    if legend.len() > 1 {
        ui.label(RichText::new(tr("multiphysics.colour_focus_hint")).small().color(t.text_dim));
    }
    for entry in &legend {
        let Ok(swatch) = serde_json::from_value::<Swatch>(entry["swatch"].clone()) else { continue };
        let id = entry["label"].as_str().unwrap_or("");
        let label = legend_label(id);
        ui.horizontal(|ui| {
            swatch_square(ui, swatch);
            // A click shows only this entry's parts: the grade plot's way to see what fails.
            let response = crate::ui_kit::button(ui, format!("{label} ({})", entry["count"].as_u64().unwrap_or(0)), crate::ui_kit::Button::Ghost);
            crate::scenario::publish_target(&format!("multiphysics:parts_legend:{id}"), response.rect);
            if response.clicked() {
                let (keep, others): (Vec<_>, Vec<_>) = parts.iter().partition(|(_, l)| l == id);
                actions.push(PartAction::Focus {
                    keep: keep.into_iter().map(|(b, _)| b.clone()).collect(),
                    others: others.into_iter().map(|(b, _)| b.clone()).collect(),
                });
            }
        });
    }
}

/// The checklist: one line per check, problems first worded, with the parts they are about.
pub(super) fn checklist_section(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui, actions: &mut Vec<PartAction>) {
    let t = Tokens::get();
    let reply = match checks(app, state, study) {
        Ok(reply) => reply,
        Err(error) => {
            ui.colored_label(t.form_error, error);
            return;
        }
    };
    let (errors, warnings) = (reply["errors"].as_u64().unwrap_or(0), reply["warnings"].as_u64().unwrap_or(0));
    ui.horizontal(|ui| {
        ui.strong(tr("multiphysics.parts_checklist"));
        let summary = solidcraft_ui_text::format(
            "multiphysics.parts_checklist_summary",
            &[("errors", &errors.to_string()), ("warnings", &warnings.to_string())],
        );
        let colour = if errors > 0 {
            t.form_error
        } else if warnings > 0 {
            t.warning
        } else {
            t.form_success
        };
        ui.colored_label(colour, summary);
    });
    for check in reply["checks"].as_array().into_iter().flatten() {
        let (id, status) = (check["check"].as_str().unwrap_or(""), check["status"].as_str().unwrap_or("ok"));
        let (mark, colour) = match status {
            "error" => ("✖", t.form_error),
            "warn" => ("⚠", t.warning),
            _ => ("✓", t.form_success),
        };
        let key = format!("multiphysics.check_{id}_{status}");
        let count = check["count"].as_u64().unwrap_or(0).to_string();
        let worded = solidcraft_ui_text::format(&key, &[("count", &count)]);
        let text = if worded == key { check["message"].as_str().unwrap_or(id).to_string() } else { worded };
        let bodies: Vec<String> = check["bodies"].as_array().into_iter().flatten().filter_map(|b| b.as_str().map(str::to_string)).collect();
        ui.horizontal_wrapped(|ui| {
            ui.colored_label(colour, mark);
            let response = ui.label(text);
            crate::scenario::publish_target(&format!("multiphysics:parts_check:{id}"), response.rect);
            if !bodies.is_empty() {
                let response = crate::ui_kit::button(ui, tr("multiphysics.parts_select"), crate::ui_kit::Button::Small);
                crate::scenario::publish_target(&format!("multiphysics:parts_check_select:{id}"), response.rect);
                if response.clicked() {
                    actions.push(PartAction::PickMany(bodies));
                }
            }
        });
    }
    capture_table(ui, &reply["capture"], actions);
}

/// How the grid holds each part: its cells, how many thick it is, and its cells' volume as a
/// share of its exact volume; the worst held first.
fn capture_table(ui: &mut egui::Ui, capture: &Value, actions: &mut Vec<PartAction>) {
    let t = Tokens::get();
    let mut rows: Vec<&Value> = capture.as_array().into_iter().flatten().collect();
    if rows.is_empty() {
        return;
    }
    // Lost first, then the thinnest, then the farthest from their own volume.
    let off = |r: &Value| r["captured"].as_f64().map_or(0.0, |c| (c - 1.0).abs());
    rows.sort_by(|a, b| {
        let key = |r: &Value| (r["cells"].as_u64() != Some(0), r["cells_across"].as_u64().unwrap_or(0));
        key(a).cmp(&key(b)).then(off(b).total_cmp(&off(a)))
    });
    let min_across = solidcraft_engine::multiphysics::setup_check::MIN_CELLS_ACROSS as u64;
    let tolerance = solidcraft_engine::multiphysics::setup_check::CAPTURE_TOLERANCE;
    egui::CollapsingHeader::new(tr("multiphysics.capture_title")).id_salt("multiphysics_capture").show(ui, |ui| {
        ui.label(RichText::new(tr("multiphysics.capture_hint")).small().color(t.text_dim));
        ui.scope(|ui| {
            ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
            egui::Grid::new("multiphysics_capture_grid").striped(true).show(ui, |ui| {
                for key in
                    ["multiphysics.parts_column_part", "multiphysics.capture_cells", "multiphysics.capture_across", "multiphysics.capture_held"]
                {
                    ui.strong(tr(key));
                }
                ui.end_row();
                for row in rows.iter().take(super::MAX_TABLE_ROWS) {
                    let body = row["body"].as_str().unwrap_or_default();
                    let response = ui.add(egui::Label::new(body).sense(egui::Sense::click()));
                    crate::scenario::publish_target(&format!("multiphysics:capture_row:{body}"), response.rect);
                    if response.clicked() {
                        actions.push(PartAction::PickMany(vec![body.to_string()]));
                    }
                    let cells = row["cells"].as_u64().unwrap_or(0);
                    let across = row["cells_across"].as_u64().unwrap_or(0);
                    ui.monospace(cells.to_string());
                    let colour = if cells == 0 {
                        t.form_error
                    } else if across < min_across {
                        t.warning
                    } else {
                        t.text
                    };
                    ui.colored_label(colour, RichText::new(across.to_string()).monospace());
                    match row["captured"].as_f64() {
                        Some(held) => {
                            let colour = if (held - 1.0).abs() > tolerance { t.warning } else { t.text };
                            ui.colored_label(colour, RichText::new(format!("{:.0} %", held * 100.0)).monospace());
                        }
                        None => {
                            ui.monospace("–");
                        }
                    }
                    ui.end_row();
                }
            });
        });
    });
}

/// Find: stacked rules on the parts' setup and the run's results; Find selects what meets them.
pub(super) fn find_section(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui, actions: &mut Vec<PartAction>) {
    let t = Tokens::get();
    ui.strong(tr("multiphysics.find_title"));
    if !state.parts_ui.find_seeded {
        state.parts_ui.find_seeded = true;
        state.parts_ui.find.push(FindRow::default());
    }
    let mut remove = None;
    for (index, row) in state.parts_ui.find.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            let (_, label, numeric) = FIND_FIELDS.get(row.field).copied().unwrap_or(FIND_FIELDS[0]);
            let field_menu = crate::ui_kit::bounded_dropdown_menu(ui, format!("multiphysics_find_field_{index}"), tr(label), 170.0, |ui| {
                for (i, (name, label, _)) in FIND_FIELDS.iter().enumerate() {
                    let response = crate::ui_kit::list_picker(ui, &format!("multiphysics_find_field_{index}_{name}"), tr(label), i == row.field);
                    crate::scenario::publish_target(&format!("multiphysics:find_field:{index}:{name}"), response.rect);
                    if response.clicked() {
                        if row.field != i {
                            row.value = if *name == "verdict" { "0".into() } else { String::new() };
                        }
                        row.field = i;
                        ui.close();
                    }
                }
            });
            crate::scenario::publish_target(&format!("multiphysics:find_field:{index}"), field_menu.rect);
            if numeric {
                let shown = COMPARES.get(row.compare).map_or("<", |c| c.1);
                crate::ui_kit::bounded_dropdown_menu(ui, format!("multiphysics_find_compare_{index}"), shown, 44.0, |ui| {
                    for (i, (name, shown)) in COMPARES.iter().enumerate() {
                        let response = crate::ui_kit::list_picker(ui, &format!("multiphysics_find_compare_{index}_{name}"), shown, i == row.compare);
                        if response.clicked() {
                            row.compare = i;
                            ui.close();
                        }
                    }
                });
            }
            if FIND_FIELDS.get(row.field).is_some_and(|f| f.0 == "verdict") {
                let current = row.value.parse::<usize>().unwrap_or(0);
                let shown = VERDICTS.get(current).map_or("", |v| tr(v.1));
                crate::ui_kit::bounded_dropdown_menu(ui, format!("multiphysics_find_verdict_{index}"), shown, 120.0, |ui| {
                    for (i, (_, label)) in VERDICTS.iter().enumerate() {
                        if crate::ui_kit::list_picker(ui, &format!("multiphysics_find_verdict_{index}_{i}"), tr(label), i == current).clicked() {
                            row.value = i.to_string();
                            ui.close();
                        }
                    }
                });
            } else {
                let width = 110.0_f32.min(ui.available_width());
                let response = ui
                    .allocate_ui(egui::vec2(width, crate::ui_kit::row_height()), |ui| {
                        crate::ui_kit::text_field_reserving(ui, &mut row.value, "", 0.0).1
                    })
                    .inner;
                crate::scenario::publish_target(&format!("multiphysics:find_value:{index}"), response.rect);
            }
            let response = crate::ui_kit::button(ui, "×", crate::ui_kit::Button::Small);
            crate::scenario::publish_target(&format!("multiphysics:find_remove:{index}"), response.rect);
            if response.clicked() {
                remove = Some(index);
            }
        });
    }
    if let Some(index) = remove {
        state.parts_ui.find.remove(index);
    }
    let rules: Option<Vec<Value>> = state.parts_ui.find.iter().map(FindRow::rule).collect();
    ui.horizontal_wrapped(|ui| {
        let more = state.parts_ui.find.len() < solidcraft_engine::multiphysics::find::MAX_FIND_RULES;
        let response = crate::ui_kit::button(ui, tr("multiphysics.find_add"), crate::ui_kit::Button::Enabled(more));
        crate::scenario::publish_target("multiphysics:find_add", response.rect);
        if response.clicked() {
            state.parts_ui.find.push(FindRow { field: 0, compare: 2, value: String::new() });
        }
        let ready = rules.as_ref().is_some_and(|r| !r.is_empty());
        let response = crate::ui_kit::button(ui, tr("multiphysics.find_run"), crate::ui_kit::Button::Enabled(ready));
        crate::scenario::publish_target("multiphysics:find_run", response.rect);
        if response.clicked()
            && let Some(rules) = &rules
        {
            state.parts_ui.found = Some(app.run("multiphysics.parts.find", json!({"id": study.id, "rules": rules, "select": true})));
        }
        // What the last Find found, beside it, and a way to see only those parts.
        if let Some(Ok(reply)) = &state.parts_ui.found {
            let found: Vec<String> = reply["found"].as_array().into_iter().flatten().filter_map(|b| b.as_str().map(str::to_string)).collect();
            let count = found.len().to_string();
            ui.label(RichText::new(solidcraft_ui_text::format("multiphysics.find_found", &[("count", &count)])).small());
            if !found.is_empty() {
                let response = crate::ui_kit::button(ui, tr("multiphysics.find_show_only"), crate::ui_kit::Button::Small);
                crate::scenario::publish_target("multiphysics:find_show_only", response.rect);
                if response.clicked() {
                    let others = bodies(app, state, study).into_iter().filter(|b| !found.contains(b)).collect();
                    actions.push(PartAction::Focus { keep: found, others });
                }
            }
        }
    });
    if let Some(Err(error)) = &state.parts_ui.found {
        ui.colored_label(t.form_error, error);
    }
}

/// Group: making sense of an imported assembly. Tag its parts by their reference designators,
/// select every part shaped like the one picked, and leave the small parts out of the study.
fn group_section(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let t = Tokens::get();
    ui.strong(tr("multiphysics.group_title"));
    let chosen = selected(app);
    ui.horizontal_wrapped(|ui| {
        let response = crate::ui_kit::button(ui, tr("multiphysics.group_designators"), crate::ui_kit::Button::Small);
        crate::scenario::publish_target("multiphysics:group_designators", response.rect);
        if response.clicked() {
            state.parts_ui.grouped = Some(app.run("multiphysics.parts.group", json!({})).map(|reply| {
                let count = reply["grouped"].as_array().map_or(0, Vec::len).to_string();
                solidcraft_ui_text::format("multiphysics.group_tagged", &[("count", &count)])
            }));
        }
        let one = chosen.len() == 1;
        let response = crate::ui_kit::button(ui, tr("multiphysics.group_similar"), crate::ui_kit::Button::Enabled(one))
            .on_hover_text(tr("multiphysics.group_similar_hint"));
        crate::scenario::publish_target("multiphysics:group_similar", response.rect);
        if response.clicked()
            && let Some(body) = chosen.first()
        {
            state.parts_ui.grouped = Some(app.run("multiphysics.parts.similar", json!({"body": body, "select": true})).map(|reply| {
                let count = reply["similar"].as_array().map_or(0, Vec::len).to_string();
                solidcraft_ui_text::format("multiphysics.group_similar_found", &[("count", &count)])
            }));
        }
    });
    ui.horizontal(|ui| {
        ui.label(tr("multiphysics.group_small"));
        let width = 80.0_f32.min(ui.available_width());
        let response = ui
            .allocate_ui(egui::vec2(width, crate::ui_kit::row_height()), |ui| {
                crate::ui_kit::text_field_reserving(ui, &mut state.parts_ui.small_mm, "", 0.0).1
            })
            .inner;
        crate::scenario::publish_target("multiphysics:group_small_mm", response.rect);
        let size = state.parts_ui.small_mm.trim().parse::<f64>().ok().filter(|v| v.is_finite() && *v > 0.0);
        let response = crate::ui_kit::button(ui, tr("multiphysics.group_small_apply"), crate::ui_kit::Button::Enabled(size.is_some()));
        crate::scenario::publish_target("multiphysics:group_small_apply", response.rect);
        if response.clicked()
            && let Some(size) = size
        {
            let params = json!({"id": study.id, "expected_revision": study.revision, "below_mm": size});
            state.parts_ui.grouped = Some(
                app.run("multiphysics.model.leave_out_small", params)
                    .map(|_| solidcraft_ui_text::format("multiphysics.group_small_done", &[("size", &format!("{size}"))])),
            );
        }
    });
    match &state.parts_ui.grouped {
        Some(Ok(done)) => {
            ui.label(RichText::new(done).small());
        }
        Some(Err(error)) => {
            ui.colored_label(t.form_error, error);
        }
        None => {}
    }
}

/// Compare runs: two runs that kept a thermal report, of this study or of any study of the
/// model (a clone, a sweep's point), part by part, the largest rise first
/// (`multiphysics.result.compare`).
fn compare_section(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let t = Tokens::get();
    // Every run with a thermal report of a study of the model: (study id, run id, its label).
    let runs: Vec<(String, String, String)> = app
        .session
        .doc
        .multiphysics
        .studies
        .iter()
        .filter(|s| s.model.is_some())
        .flat_map(|s| {
            s.runs
                .iter()
                .filter(|r| r.artifacts.iter().any(|a| a.name == "thermal_report.json"))
                .map(move |r| (s.id.clone(), r.id.clone(), format!("{} · {}", s.name, r.id)))
        })
        .collect();
    ui.separator();
    ui.strong(tr("multiphysics.compare_title"));
    if runs.len() < 2 {
        ui.label(RichText::new(tr("multiphysics.compare_needs_two")).small().color(t.text_dim));
        return;
    }
    // By default this study's last run against the run before it.
    let mine: Vec<usize> = runs.iter().enumerate().filter(|(_, r)| r.0 == study.id).map(|(i, _)| i).collect();
    let latest = mine.last().copied().unwrap_or(runs.len() - 1);
    let defaults = [mine.iter().rev().nth(1).copied().unwrap_or(if latest == 0 { 1 } else { latest - 1 }), latest];
    let ui_state = &mut state.parts_ui;
    for (slot, default) in ui_state.compare.iter_mut().zip(defaults) {
        if slot.as_ref().is_none_or(|(s, r)| !runs.iter().any(|x| &x.0 == s && &x.1 == r)) {
            *slot = runs.get(default).map(|r| (r.0.clone(), r.1.clone()));
        }
    }
    for (index, key) in ["multiphysics.compare_first", "multiphysics.compare_second"].into_iter().enumerate() {
        let current = ui_state.compare[index].clone();
        let shown = runs.iter().find(|r| current.as_ref().is_some_and(|(s, id)| &r.0 == s && &r.1 == id)).map_or("", |r| r.2.as_str());
        crate::ui_kit::row(ui, tr(key), |ui| {
            crate::ui_kit::dropdown_menu(ui, format!("multiphysics_compare_{index}"), shown, |ui| {
                for (study_id, run, label) in &runs {
                    let on = current.as_ref().is_some_and(|(s, id)| s == study_id && id == run);
                    let response = crate::ui_kit::list_picker(ui, &format!("multiphysics_compare_{index}_{study_id}_{run}"), label, on);
                    crate::scenario::publish_target(&format!("multiphysics:compare_{index}:{label}"), response.rect);
                    if response.clicked() {
                        ui_state.compare[index] = Some((study_id.clone(), run.clone()));
                        ui_state.compared = None;
                        ui.close();
                    }
                }
            });
        });
    }
    let response =
        crate::ui_kit::button(ui, tr("multiphysics.compare_run"), crate::ui_kit::Button::Enabled(ui_state.compare[0] != ui_state.compare[1]));
    crate::scenario::publish_target("multiphysics:compare_run", response.rect);
    if response.clicked()
        && let [Some((first_id, first_run)), Some((second_id, second_run))] = &ui_state.compare
    {
        let params = json!({"first": {"id": first_id, "run_id": first_run}, "second": {"id": second_id, "run_id": second_run}});
        ui_state.compared = Some(app.run("multiphysics.result.compare", params).map_err(|e| e.to_string()));
    }
    let reply = match &state.parts_ui.compared {
        Some(Ok(reply)) => reply,
        Some(Err(error)) => {
            ui.colored_label(t.form_error, error);
            return;
        }
        None => return,
    };
    let comparison = &reply["comparison"];
    let count = |key: &str| comparison[key].as_u64().unwrap_or(0).to_string();
    let power = format!("{:+.2}", comparison["power_delta_w"].as_f64().unwrap_or(0.0));
    let summary = solidcraft_ui_text::format(
        "multiphysics.compare_summary",
        &[("warmer", &count("warmer")), ("cooler", &count("cooler")), ("worse", &count("worse")), ("better", &count("better")), ("watts", &power)],
    );
    ui.label(RichText::new(summary).small());
    ui.label(RichText::new(tr("multiphysics.compare_units")).small().color(t.text_dim));
    let rows = comparison["rows"].as_array().cloned().unwrap_or_default();
    let celsius = |v: &Value| v.as_f64().map_or_else(|| "–".to_string(), |c| format!("{c:.1}"));
    let verdict = |v: &Value| match v.as_str() {
        Some("pass") => tr("multiphysics.post_status_pass"),
        Some("warn") => tr("multiphysics.post_status_warn"),
        Some("fail") => tr("multiphysics.post_status_fail"),
        Some(_) => tr("multiphysics.post_status_no_limit"),
        None => tr("multiphysics.parts_no_result"),
    };
    // In the panel's own scroll (a table scrolling inside it shows a row or two), unwrapped:
    // the headings stay whole however narrow the panel.
    ui.scope(|ui| {
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
        egui::Grid::new("multiphysics_compare_grid").striped(true).show(ui, |ui| {
            for key in [
                "multiphysics.parts_column_part",
                "multiphysics.compare_column_first",
                "multiphysics.compare_column_second",
                "multiphysics.compare_column_delta",
                "multiphysics.parts_column_status",
            ] {
                ui.strong(tr(key));
            }
            ui.end_row();
            for row in rows.iter().take(super::MAX_TABLE_ROWS) {
                ui.label(row["body"].as_str().unwrap_or_default());
                ui.monospace(celsius(&row["first_c"]));
                ui.monospace(celsius(&row["second_c"]));
                let delta = row["delta_k"].as_f64();
                let colour = match delta {
                    Some(d) if d >= solidcraft_engine::multiphysics::variants::SAME_K => t.warning,
                    Some(d) if d <= -solidcraft_engine::multiphysics::variants::SAME_K => t.form_success,
                    _ => t.text_dim,
                };
                ui.colored_label(colour, RichText::new(delta.map_or_else(|| "–".to_string(), |d| format!("{d:+.1}"))).monospace());
                let colour = match row["change"].as_str() {
                    Some("worse") => t.form_error,
                    Some("better") => t.form_success,
                    _ => t.text,
                };
                // One verdict where it stayed, both where it moved.
                let text = if row["first_status"] == row["second_status"] {
                    verdict(&row["second_status"]).to_string()
                } else {
                    format!("{} → {}", verdict(&row["first_status"]), verdict(&row["second_status"]))
                };
                ui.colored_label(colour, text);
                ui.end_row();
            }
        });
    });
}

/// Mesh study: the study on two coarser grids, and once all three have run, how far each part's
/// temperature is from settling (`multiphysics.model.mesh_study`, `multiphysics.result.convergence`).
fn mesh_study_section(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let t = Tokens::get();
    if study.model.is_none() {
        return;
    }
    ui.separator();
    ui.strong(tr("multiphysics.mesh_study_title"));
    let studies = &app.session.doc.multiphysics.studies;
    let clone = |label: &str| studies.iter().find(|s| s.name == format!("{} · mesh ×{label}", study.name)).map(|s| s.id.clone());
    let thermal = |id: &str| {
        studies.iter().find(|s| s.id == id).is_some_and(|s| s.runs.iter().any(|r| r.artifacts.iter().any(|a| a.name == "thermal_report.json")))
    };
    match (clone("1.41"), clone("2")) {
        (Some(medium), Some(coarse)) => {
            let ready = [study.id.as_str(), &medium, &coarse].iter().all(|id| thermal(id));
            let key = if ready { "multiphysics.mesh_study_ready" } else { "multiphysics.mesh_study_run_all" };
            ui.label(RichText::new(tr(key)).small().color(t.text_dim));
            let response = crate::ui_kit::button(ui, tr("multiphysics.mesh_study_check"), crate::ui_kit::Button::Enabled(ready));
            crate::scenario::publish_target("multiphysics:mesh_study_check", response.rect);
            if response.clicked() {
                let params = json!({"fine": {"id": study.id}, "medium": {"id": medium}, "coarse": {"id": coarse}});
                state.parts_ui.mesh_study = Some(app.run("multiphysics.result.convergence", params).map_err(|e| e.to_string()));
            }
        }
        _ => {
            ui.label(RichText::new(tr("multiphysics.mesh_study_hint")).small().color(t.text_dim));
            let response = crate::ui_kit::button(ui, tr("multiphysics.mesh_study_make"), crate::ui_kit::Button::Small);
            crate::scenario::publish_target("multiphysics:mesh_study_make", response.rect);
            if response.clicked() {
                state.parts_ui.mesh_study = Some(app.run("multiphysics.model.mesh_study", json!({"id": study.id})).map_err(|e| e.to_string()));
            }
        }
    }
    let reply = match &state.parts_ui.mesh_study {
        Some(Ok(reply)) => reply,
        Some(Err(error)) => {
            ui.colored_label(t.form_error, error);
            return;
        }
        None => return,
    };
    // Grids made: say so; convergence: the table.
    let Some(parts) = reply["parts"].as_array() else {
        if let Some(cells) = reply["cells"].as_array() {
            let cells: Vec<String> = cells.iter().map(|c| format!("{:.0}", c.as_f64().unwrap_or(0.0))).collect();
            ui.label(RichText::new(solidcraft_ui_text::format("multiphysics.mesh_study_made", &[("cells", &cells.join(" / "))])).small());
        }
        return;
    };
    for note in reply["notes"].as_array().into_iter().flatten().filter_map(Value::as_str) {
        ui.colored_label(t.warning, RichText::new(note).small());
    }
    ui.label(RichText::new(tr("multiphysics.mesh_study_units")).small().color(t.text_dim));
    let number = |v: &Value, digits: usize| v.as_f64().map_or_else(|| "–".to_string(), |x| format!("{x:.digits$}"));
    let inflow = reply["inflow_c"].as_f64().unwrap_or(0.0);
    ui.scope(|ui| {
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
        egui::Grid::new("multiphysics_mesh_study_grid").striped(true).show(ui, |ui| {
            for key in [
                "multiphysics.parts_column_part",
                "multiphysics.mesh_study_column_fine",
                "multiphysics.mesh_study_column_settled",
                "multiphysics.mesh_study_column_gci",
                "multiphysics.mesh_study_column_order",
            ] {
                ui.strong(tr(key));
            }
            ui.end_row();
            for part in parts.iter().take(super::MAX_TABLE_ROWS) {
                let estimate = &part["estimate"];
                ui.label(part["body"].as_str().unwrap_or_default());
                ui.monospace(number(&part["fine_c"], 1));
                ui.monospace(estimate["extrapolated"].as_f64().map_or_else(|| "–".to_string(), |rise| format!("{:.1}", rise + inflow)));
                let gci = estimate["gci_percent"].as_f64();
                let colour = match (estimate["behaviour"].as_str(), gci) {
                    (Some("diverging"), _) => t.form_error,
                    (Some("oscillatory"), _) => t.warning,
                    (_, Some(g)) if g > 5.0 => t.warning,
                    _ => t.text,
                };
                let text = match estimate["behaviour"].as_str() {
                    Some("diverging") => tr("multiphysics.mesh_study_diverging").to_string(),
                    _ => gci.map_or_else(|| "–".to_string(), |g| format!("{g:.1} %")),
                };
                ui.colored_label(colour, RichText::new(text).monospace());
                ui.monospace(number(&estimate["order"], 2));
                ui.end_row();
            }
        });
    });
}

/// Batch report: the component table of every run that kept a thermal report, side by side, as
/// CSV and a self-contained HTML page in a new directory.
#[cfg(not(target_arch = "wasm32"))]
fn report_section(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let t = Tokens::get();
    let runs = study.runs.iter().filter(|r| r.artifacts.iter().any(|a| a.name == "thermal_report.json")).count();
    ui.separator();
    ui.strong(tr("multiphysics.report_title"));
    ui.label(RichText::new(solidcraft_ui_text::format("multiphysics.report_runs", &[("count", &runs.to_string())])).small().color(t.text_dim));
    ui.horizontal(|ui| {
        ui.label(tr("multiphysics.report_directory"));
        let width = 220.0_f32.min(ui.available_width());
        let response = ui
            .allocate_ui(egui::vec2(width, crate::ui_kit::row_height()), |ui| {
                crate::ui_kit::text_field_reserving(ui, &mut state.parts_ui.report_dir, "", 0.0).1
            })
            .inner;
        crate::scenario::publish_target("multiphysics:report_dir", response.rect);
        let ready = runs > 0 && !state.parts_ui.report_dir.trim().is_empty();
        let response = crate::ui_kit::button(ui, tr("multiphysics.report_write"), crate::ui_kit::Button::Enabled(ready));
        crate::scenario::publish_target("multiphysics:report_write", response.rect);
        if response.clicked() {
            let directory = state.parts_ui.report_dir.trim().to_string();
            let params = json!({"id": study.id, "directory": directory});
            state.parts_ui.reported = Some(app.run("multiphysics.result.report", params).map(|reply| {
                let count = reply["runs"].as_array().map_or(0, Vec::len).to_string();
                solidcraft_ui_text::format("multiphysics.report_written", &[("count", &count), ("directory", &directory)])
            }));
        }
    });
    match &state.parts_ui.reported {
        Some(Ok(done)) => {
            ui.label(RichText::new(done).small());
        }
        Some(Err(error)) => {
            ui.colored_label(t.form_error, error);
        }
        None => {}
    }
}
