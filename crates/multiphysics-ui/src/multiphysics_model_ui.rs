//! The Multiphysics workspace for a study of the loaded model.
//!
//! The graphics window stays the design's own 3-D view (one camera, the normal selection), with
//! the computational domain and the conditions drawn over the model. The panel sets whether the
//! fluid is inside or around the bodies, the cell size, and conditions on faces picked in the
//! view. Every change is a `multiphysics.model.*` command that builds the case from the design
//! again, so the case never drifts from the model.
use egui::{Align2, Color32, FontId, RichText, Stroke};
use serde_json::{Value, json};
use solidcraft_theme::Tokens;

use super::{State, button, command, tr};
use crate::SolveApp;
use crate::ui_kit::BackingTooltipExt as _;
use solidcraft_engine::multiphysics::Study;

/// The condition kinds the panel offers, in order, as `multiphysics.model.condition` names them.
const ROLES: [&str; 4] = ["inlet", "outlet", "wall", "symmetry"];

/// What the panel remembers between frames.
#[derive(Default)]
pub(super) struct ModelUi {
    /// The last `multiphysics.model.status` reply, for the document it was read from.
    status: Option<(usize, String, Result<Value, String>)>,
    /// The last build's reply (placed faces, solid fraction, notes).
    pub(super) last_build: Option<(String, Value)>,
    /// The condition being added.
    name: String,
    role: usize,
    value: String,
    cell: String,
    margin: String,
    /// The end time as typed, for the study it was typed for.
    end_time: Option<(String, String)>,
    /// A porous body's Darcy (1/m²) and Forchheimer (1/m) coefficients, as typed.
    porous_d: String,
    porous_c: String,
    /// The last leak search: the study, whether fluid connects the faces, and the path (mm).
    leak: Option<(String, bool, Vec<[f64; 3]>)>,
    /// Local mesh: the split factor offered, 2 to 4 (0 means 2).
    refine_factor: usize,
    /// The goal being added: its kind (index in [`GOAL_KINDS`]), name and equation.
    goal_kind: usize,
    goal_name: String,
    goal_expression: String,
    /// Heat: the inflow temperature (°C), a material's conductivity and power, a heat source's
    /// name and power, and a thermal condition's kind (index in [`THERMAL_KINDS`]) and values.
    inflow_c: String,
    conductivity: String,
    power: String,
    source_name: String,
    source_power: String,
    thermal_kind: usize,
    thermal_value: String,
    thermal_ambient: String,
    /// The mesh level being dragged (applied when the drag ends), and the grid drawn on the model.
    level: Option<u8>,
    pub(super) show_grid: bool,
    /// What the last adaptation found, for the study.
    adapted: Option<(String, Value)>,
    /// The fan curve being typed: `flow pressure; flow pressure; …` (m³/s, Pa).
    fan_curve: String,
    /// Environment: the breeze (m/s), the dust (mg/m³), the wall inset (mm) of a dust box made
    /// from a body, the seasons being typed (low, high, dew point °C, dust mg/m³ each), the
    /// folder a year summary goes to, and the last summary, for the study.
    breeze: String,
    dust_mg: String,
    dust_inset: String,
    seasons: Vec<[String; 4]>,
    climate_folder: String,
    climate: Option<(String, Value)>,
}

/// A temperate year to start from: winter, spring, summer, autumn (low, high, dew point °C,
/// dust mg/m³ by day).
const SEASONS: [(&str, [f64; 4]); 4] =
    [("winter", [-4.0, 3.0, -5.0, 0.1]), ("spring", [4.0, 16.0, 3.0, 0.4]), ("summer", [14.0, 30.0, 13.0, 0.8]), ("autumn", [6.0, 16.0, 7.0, 0.5])];

/// The goal kinds the panel offers, as `multiphysics.model.goal` names them; the temperature
/// goals only with heat transfer on.
const GOAL_KINDS: &[&str] = &[
    "face_volume_flow",
    "face_mean_pressure",
    "point_pressure",
    "equation",
    "global_max_speed",
    "global_mean_speed",
    "global_max_pressure",
    "global_min_pressure",
    "global_mean_pressure",
    "body_max_temperature",
    "body_mean_temperature",
    "point_temperature",
    "global_max_fluid_temperature",
    "global_mean_fluid_temperature",
    "global_max_solid_temperature",
];
/// The goal kinds without heat transfer: the first nine of [`GOAL_KINDS`].
const FLOW_GOAL_KINDS: usize = 9;

/// The thermal conditions on faces on the domain boundary, as `multiphysics.model.thermal_face` names them.
const THERMAL_KINDS: [&str; 4] = ["temperature", "heat_flux", "convection", "adiabatic"];

/// A row's remove button: the delete icon, named for screen readers and its tooltip.
fn remove_button(ui: &mut egui::Ui) -> egui::Response {
    let size = Tokens::get().control_height;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(size, size), egui::Sense::hover());
    crate::ui_kit::quiet_icon_button_at(ui, rect, "delete", tr("multiphysics.model_remove")).backing_on_hover_text(tr("multiphysics.model_remove"))
}

/// A text field published under `handle`, so scripts and agents can type into it.
fn field(ui: &mut egui::Ui, value: &mut String, hint: &str, handle: &str) {
    let (_, response) = crate::ui_kit::text_field_reserving(ui, value, hint, 0.0);
    crate::scenario::publish_target(handle, response.rect);
}

/// Kelvin from a typed Celsius value.
fn kelvin(text: &str) -> Option<f64> {
    text.trim().parse::<f64>().ok().filter(|c| c.is_finite()).map(|c| c + 273.15)
}

/// Heat transfer in a study of the model (the Solid Materials node): on or off, the inflow
/// temperature, the bodies' materials and power, heat sources on faces and thermal conditions
/// on the faces on the domain boundary.
pub(super) fn heat_panel(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let Some(link) = study.model.clone() else { return };
    let t = Tokens::get();
    let base = || json!({"id":study.id,"expected_revision":study.revision});
    let mut on = link.heat.is_some();
    let response = crate::ui_kit::checkbox(ui, tr("multiphysics.model_heat"), &mut on);
    crate::scenario::publish_target("multiphysics:model_heat", response.rect);
    ui.label(RichText::new(tr("multiphysics.model_heat_hint")).small().color(t.text_dim));
    if on != link.heat.is_some() {
        let mut params = base();
        params["enabled"] = json!(on);
        model_command(app, state, study, "multiphysics.model.heat", params);
        return;
    }
    let Some(heat) = link.heat.clone() else { return };
    if let Some((id, reply)) = &state.model_ui.last_build
        && *id == study.id
    {
        for note in reply["notes"].as_array().into_iter().flatten().filter_map(Value::as_str) {
            ui.colored_label(t.warning, note);
        }
    }
    // The fluid flowing in.
    if state.model_ui.inflow_c.is_empty() {
        state.model_ui.inflow_c = format!("{}", ((heat.inflow_temperature_k - 273.15) * 100.0).round() / 100.0);
    }
    let mut apply_inflow = false;
    crate::ui_kit::row(ui, tr("multiphysics.model_inflow_temperature"), |ui| {
        field(ui, &mut state.model_ui.inflow_c, "20", "multiphysics:model_inflow");
        let response = crate::ui_kit::button(ui, tr("multiphysics.model_apply"), crate::ui_kit::Button::Small);
        crate::scenario::publish_target("multiphysics:model_apply_inflow", response.rect);
        apply_inflow = response.clicked();
    });
    if apply_inflow && let Some(k) = kelvin(&state.model_ui.inflow_c) {
        let mut params = base();
        params["inflow_temperature_k"] = json!(k);
        model_command(app, state, study, "multiphysics.model.heat", params);
        return;
    }
    // Natural convection: the temperature drives the flow.
    let mut buoyant = heat.buoyancy.is_some();
    let response = crate::ui_kit::checkbox(ui, tr("multiphysics.model_buoyancy"), &mut buoyant);
    crate::scenario::publish_target("multiphysics:model_buoyancy", response.rect);
    ui.label(RichText::new(tr("multiphysics.model_buoyancy_hint")).small().color(t.text_dim));
    if buoyant != heat.buoyancy.is_some() {
        let mut params = base();
        params["buoyancy"] = json!(buoyant);
        model_command(app, state, study, "multiphysics.model.heat", params);
        return;
    }
    ui.separator();
    // Materials and power of the bodies.
    ui.strong(tr("multiphysics.model_materials"));
    let mut remove: Option<(&str, Value)> = None;
    for m in &heat.materials {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {} W/(m·K) · {} W", m.body, m.conductivity_w_per_m_k, m.power_w));
            let response = remove_button(ui);
            crate::scenario::publish_target(&format!("multiphysics:model_remove_material:{}", m.body), response.rect);
            if response.clicked() {
                remove = Some(("multiphysics.model.material", json!({"bodies":[m.body],"remove":true})));
            }
        });
    }
    if heat.materials.is_empty() {
        ui.label(RichText::new(tr("multiphysics.model_material_starter")).small().color(t.text_dim));
    }
    crate::ui_kit::row(ui, tr("multiphysics.model_conductivity"), |ui| {
        field(ui, &mut state.model_ui.conductivity, "200", "multiphysics:model_conductivity")
    });
    crate::ui_kit::row(ui, tr("multiphysics.model_power"), |ui| field(ui, &mut state.model_ui.power, "0", "multiphysics:model_power"));
    let bodies = selected_bodies(app);
    let conductivity = state.model_ui.conductivity.trim().parse::<f64>().ok();
    if button(ui, "multiphysics.model_apply_material", "multiphysics:model_apply_material", !bodies.is_empty() && conductivity.is_some()).clicked() {
        let mut params = base();
        params["bodies"] = json!(bodies);
        params["conductivity_w_per_m_k"] = json!(conductivity);
        params["power_w"] = json!(state.model_ui.power.trim().parse::<f64>().unwrap_or(0.0));
        model_command(app, state, study, "multiphysics.model.material", params);
        return;
    }
    if bodies.is_empty() {
        ui.label(RichText::new(tr("multiphysics.model_pick_bodies")).small().color(t.text_dim));
    }
    ui.separator();
    // Heat sources on faces.
    ui.strong(tr("multiphysics.model_heat_sources"));
    for source in &heat.face_sources {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {} W · {}", source.name, source.power_w, source.face.body));
            let response = remove_button(ui);
            crate::scenario::publish_target(&format!("multiphysics:model_remove_heat_source:{}", source.name), response.rect);
            if response.clicked() {
                remove = Some(("multiphysics.model.heat_source", json!({"name":source.name,"remove":true})));
            }
        });
    }
    crate::ui_kit::row(ui, tr("multiphysics.name"), |ui| {
        crate::ui_kit::text_field(ui, &mut state.model_ui.source_name, &format!("source {}", heat.face_sources.len() + 1))
    });
    crate::ui_kit::row(ui, tr("multiphysics.model_power"), |ui| field(ui, &mut state.model_ui.source_power, "1", "multiphysics:model_source_power"));
    let faces = selected_faces(app);
    let power = state.model_ui.source_power.trim().parse::<f64>().ok();
    if button(ui, "multiphysics.model_add_heat_source", "multiphysics:model_add_heat_source", faces.len() == 1 && power.is_some()).clicked()
        && let Some(face) = faces.first()
    {
        let name = if state.model_ui.source_name.trim().is_empty() {
            format!("source {}", heat.face_sources.len() + 1)
        } else {
            state.model_ui.source_name.trim().to_string()
        };
        let mut params = base();
        params["name"] = json!(name);
        params["face"] = face.clone();
        params["power_w"] = json!(power);
        if model_command(app, state, study, "multiphysics.model.heat_source", params) {
            state.model_ui.source_name.clear();
        }
        return;
    }
    ui.separator();
    // Thermal conditions on the domain boundary.
    ui.strong(tr("multiphysics.model_thermal_faces"));
    for condition in &heat.face_conditions {
        ui.horizontal(|ui| {
            let kind = serde_json::to_value(condition.condition).ok().and_then(|v| v["kind"].as_str().map(str::to_string)).unwrap_or_default();
            ui.label(format!("{} · {}", condition.name, tr(&format!("multiphysics.model_thermal_{kind}"))));
            let response = remove_button(ui);
            crate::scenario::publish_target(&format!("multiphysics:model_remove_thermal_face:{}", condition.name), response.rect);
            if response.clicked() {
                remove = Some(("multiphysics.model.thermal_face", json!({"name":condition.name,"remove":true})));
            }
        });
    }
    if let Some((command_id, extra)) = remove {
        let mut params = base();
        if let (Some(target), Some(extra)) = (params.as_object_mut(), extra.as_object()) {
            target.extend(extra.clone());
        }
        model_command(app, state, study, command_id, params);
        return;
    }
    let kind = THERMAL_KINDS.get(state.model_ui.thermal_kind).copied().unwrap_or(THERMAL_KINDS[0]);
    crate::ui_kit::row(ui, tr("multiphysics.model_thermal_kind"), |ui| {
        crate::ui_kit::dropdown_menu(ui, "multiphysics_model_thermal_kind", tr(&format!("multiphysics.model_thermal_{kind}")), |ui| {
            for (index, k) in THERMAL_KINDS.iter().enumerate() {
                let response = crate::ui_kit::list_picker(
                    ui,
                    &format!("multiphysics_model_thermal_kind_{k}"),
                    tr(&format!("multiphysics.model_thermal_{k}")),
                    index == state.model_ui.thermal_kind,
                );
                crate::scenario::publish_target(&format!("multiphysics:model_thermal_kind:{k}"), response.rect);
                if response.clicked() {
                    state.model_ui.thermal_kind = index;
                    ui.close();
                }
            }
        });
    });
    let value_label = match kind {
        "temperature" => Some("multiphysics.model_temperature_c"),
        "heat_flux" => Some("multiphysics.model_heat_flux"),
        "convection" => Some("multiphysics.model_htc"),
        _ => None,
    };
    if let Some(label) = value_label {
        crate::ui_kit::row(ui, tr(label), |ui| field(ui, &mut state.model_ui.thermal_value, "", "multiphysics:model_thermal_value"));
    }
    if kind == "convection" {
        crate::ui_kit::row(ui, tr("multiphysics.model_ambient_c"), |ui| {
            field(ui, &mut state.model_ui.thermal_ambient, "20", "multiphysics:model_thermal_ambient")
        });
    }
    let value = state.model_ui.thermal_value.trim().parse::<f64>().ok().filter(|v| v.is_finite());
    let condition = match kind {
        "temperature" => kelvin(&state.model_ui.thermal_value).map(|k| json!({"kind":"temperature","temperature_k":k})),
        "heat_flux" => value.map(|w| json!({"kind":"heat_flux","w_per_m2":w})),
        "convection" => value.zip(kelvin(&state.model_ui.thermal_ambient)).map(|(h, k)| json!({"kind":"convection","h_w_per_m2_k":h,"ambient_k":k})),
        _ => Some(json!({"kind":"adiabatic"})),
    };
    if button(ui, "multiphysics.model_add_thermal_face", "multiphysics:model_add_thermal_face", faces.len() == 1 && condition.is_some()).clicked()
        && let (Some(face), Some(condition)) = (faces.first(), condition)
    {
        let mut params = base();
        params["name"] = json!(format!("{} {}", tr(&format!("multiphysics.model_thermal_{kind}")), heat.face_conditions.len() + 1));
        params["face"] = face.clone();
        params["condition"] = condition;
        model_command(app, state, study, "multiphysics.model.thermal_face", params);
        return;
    }
    ui.label(RichText::new(tr("multiphysics.model_thermal_face_hint")).small().color(t.text_dim));
}

/// The bodies in the design's selection.
fn selected_bodies(app: &SolveApp) -> Vec<String> {
    app.session
        .selection
        .iter()
        .filter_map(|sel| match sel {
            solidcraft_engine::Sel::Body { name } => Some(name.clone()),
            _ => None,
        })
        .collect()
}

/// The Mesh node of a study of the model: local refinements around its bodies and faces.
pub(super) fn mesh_panel(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let Some(link) = study.model.clone() else { return };
    let t = Tokens::get();
    let faces = selected_faces(app);
    let bodies = selected_bodies(app);
    // The level of the initial mesh: the automatic cells scaled by √2 per level.
    ui.strong(tr("multiphysics.model_mesh_level"));
    let current = link.mesh_level.unwrap_or(solidcraft_engine::multiphysics::model::AUTOMATIC_MESH_LEVEL);
    let mut level = state.model_ui.level.unwrap_or(current);
    let options = crate::ui_kit::Slider { integer: true, format: |v| format!("{v:.0}"), ..Default::default() };
    let response = crate::ui_kit::slider(ui, &mut level, 1..=solidcraft_engine::multiphysics::model::MAX_MESH_LEVEL, options);
    crate::scenario::publish_target("multiphysics:model_mesh_level", response.rect);
    state.model_ui.level = Some(level);
    let settled = response.drag_stopped() || (response.changed() && !response.dragged());
    if settled && (level != current || link.cell_mm.is_some()) {
        state.model_ui.level = None;
        model_command(app, state, study, "multiphysics.model.set", json!({"id":study.id,"expected_revision":study.revision,"mesh_level":level}));
        return;
    }
    if !response.dragged() && !response.has_focus() {
        state.model_ui.level = None;
    }
    let hint = if link.cell_mm.is_some() { "multiphysics.model_mesh_level_typed" } else { "multiphysics.model_mesh_level_hint" };
    ui.label(RichText::new(tr(hint)).small().color(t.text_dim));
    if let Ok((dimensions, _)) = solidcraft_engine::multiphysics::case_scene::grid(&study.case) {
        ui.monospace(format!("{} × {} × {} = {}", dimensions[0], dimensions[1], dimensions[2], dimensions.iter().product::<usize>()));
    }
    let response = crate::ui_kit::checkbox(ui, tr("multiphysics.model_show_grid"), &mut state.model_ui.show_grid);
    crate::scenario::publish_target("multiphysics:model_show_grid", response.rect);
    ui.separator();
    // Solution-adaptive refinement from a kept run.
    ui.strong(tr("multiphysics.model_adaptive"));
    let run = state.view.shown_run.clone().filter(|id| study.runs.iter().any(|r| &r.id == id)).or_else(|| study.runs.last().map(|r| r.id.clone()));
    if button(ui, "multiphysics.model_adapt", "multiphysics:model_adapt", run.is_some()).clicked()
        && let Some(run) = run
    {
        if let Some(reply) = command(app, state, "multiphysics.model.adapt", json!({"id":study.id,"expected_revision":study.revision,"run_id":run})) {
            state.model_ui.adapted = Some((study.id.clone(), reply["adapted"].clone()));
            state.model_ui.last_build = Some((study.id.clone(), reply));
        }
        return;
    }
    let adaptive = link.refinements.iter().any(|r| r.name.starts_with("adapt_"));
    if adaptive && button(ui, "multiphysics.model_adapt_remove", "multiphysics:model_adapt_remove", true).clicked() {
        model_command(
            app,
            state,
            study,
            "multiphysics.model.adapt",
            json!({"id":study.id,"expected_revision":study.revision,"run_id":"","remove":true}),
        );
        state.model_ui.adapted = None;
        return;
    }
    match &state.model_ui.adapted {
        Some((id, found)) if *id == study.id => {
            let regions = found["regions_mm"].as_array().map_or(0, Vec::len).to_string();
            ui.label(solidcraft_ui_text::format(
                "multiphysics.model_adapt_result",
                &[("flagged", &found["flagged"].to_string()), ("cells", &found["cells"].to_string()), ("regions", &regions)],
            ));
        }
        _ => {
            ui.label(
                RichText::new(tr(if run.is_some() { "multiphysics.model_adapt_hint" } else { "multiphysics.model_adapt_no_run" }))
                    .small()
                    .color(t.text_dim),
            );
        }
    }
    ui.separator();
    ui.strong(tr("multiphysics.model_local_mesh"));
    let mut remove: Option<(&str, String)> = None;
    for r in &link.refinements {
        ui.horizontal(|ui| {
            let what = match &r.target {
                solidcraft_engine::multiphysics::model::RefineTarget::Body { body } => body.clone(),
                solidcraft_engine::multiphysics::model::RefineTarget::Face { face } => format!("{} · {}", face.body, face.face),
                solidcraft_engine::multiphysics::model::RefineTarget::Region { .. } => tr("multiphysics.model_adaptive_box").to_string(),
            };
            ui.label(format!("{} · ×{} · {what}", r.name, r.factor));
            let response = remove_button(ui);
            crate::scenario::publish_target(&format!("multiphysics:model_remove_refinement:{}", r.name), response.rect);
            if response.clicked() {
                remove = Some(("multiphysics.model.refine", r.name.clone()));
            }
        });
    }
    let factor = state.model_ui.refine_factor.clamp(2, 4);
    crate::ui_kit::row(ui, tr("multiphysics.model_refine_factor"), |ui| {
        crate::ui_kit::dropdown_menu(ui, "multiphysics_model_refine_factor", &format!("×{factor}"), |ui| {
            for f in 2..=4 {
                if crate::ui_kit::list_picker(ui, &format!("multiphysics_model_refine_factor_{f}"), &format!("×{f}"), f == factor).clicked() {
                    state.model_ui.refine_factor = f;
                    ui.close();
                }
            }
        });
    });
    let target = match (faces.as_slice(), bodies.as_slice()) {
        ([face], []) => Some(json!({"face":face})),
        ([], [body]) => Some(json!({"body":body})),
        _ => None,
    };
    if button(ui, "multiphysics.model_refine", "multiphysics:model_refine", target.is_some()).clicked()
        && let Some(target) = target
    {
        let mut params = json!({"id":study.id,"expected_revision":study.revision,"name":format!("refine {}", link.refinements.len() + 1),"factor":factor,"margin_cells":1.0});
        if let (Some(object), Some(extra)) = (params.as_object_mut(), target.as_object()) {
            object.extend(extra.clone());
        }
        model_command(app, state, study, "multiphysics.model.refine", params);
        return;
    }
    ui.label(RichText::new(tr("multiphysics.model_refine_hint")).small().color(t.text_dim));
    if let Some((command_id, name)) = remove {
        model_command(app, state, study, command_id, json!({"id":study.id,"expected_revision":study.revision,"name":name,"remove":true}));
        return;
    }
    if porous_section(app, state, study, &link, ui) {
        return;
    }
    fluid_section(app, state, study, &link, ui);
}

/// The Goals node of a study of the model: quantities monitored on its faces and in its fluid.
pub(super) fn goals_panel(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let Some(link) = study.model.clone() else { return };
    let t = Tokens::get();
    let faces = selected_faces(app);
    let mut remove: Option<(&str, String)> = None;
    ui.strong(tr("multiphysics.model_goals"));
    for g in &link.goals {
        ui.horizontal(|ui| {
            let kind = serde_json::to_value(g.kind).ok().and_then(|v| v.as_str().map(str::to_string)).unwrap_or_default();
            ui.label(format!("{} · {}", g.name, tr(&format!("multiphysics.model_goal_{kind}"))));
            let response = remove_button(ui);
            crate::scenario::publish_target(&format!("multiphysics:model_remove_goal:{}", g.name), response.rect);
            if response.clicked() {
                remove = Some(("multiphysics.model.goal", g.name.clone()));
            }
        });
    }
    if let Some((command_id, name)) = remove {
        model_command(app, state, study, command_id, json!({"id":study.id,"expected_revision":study.revision,"name":name,"remove":true}));
        return;
    }
    let offered = if link.heat.is_some() { GOAL_KINDS.len() } else { FLOW_GOAL_KINDS };
    if state.model_ui.goal_kind >= offered {
        state.model_ui.goal_kind = 0;
    }
    let kind = GOAL_KINDS.get(state.model_ui.goal_kind).copied().unwrap_or(GOAL_KINDS[0]);
    let bodies = selected_bodies(app);
    crate::ui_kit::row(ui, tr("multiphysics.name"), |ui| {
        crate::ui_kit::text_field(ui, &mut state.model_ui.goal_name, &format!("goal_{}", link.goals.len() + 1))
    });
    crate::ui_kit::row(ui, tr("multiphysics.model_goal_kind"), |ui| {
        let menu = crate::ui_kit::dropdown_menu(ui, "multiphysics_model_goal_kind", tr(&format!("multiphysics.model_goal_{kind}")), |ui| {
            for (index, k) in GOAL_KINDS.iter().enumerate().take(offered) {
                let response = crate::ui_kit::list_picker(
                    ui,
                    &format!("multiphysics_model_goal_kind_{k}"),
                    tr(&format!("multiphysics.model_goal_{k}")),
                    index == state.model_ui.goal_kind,
                );
                crate::scenario::publish_target(&format!("multiphysics:model_goal_kind:{k}"), response.rect);
                if response.clicked() {
                    state.model_ui.goal_kind = index;
                    ui.close();
                }
            }
        });
        crate::scenario::publish_target("multiphysics:model_goal_kind_menu", menu.rect);
    });
    if kind == "equation" {
        crate::ui_kit::row(ui, tr("multiphysics.model_goal_expression"), |ui| {
            crate::ui_kit::text_field(ui, &mut state.model_ui.goal_expression, "a - b")
        });
    }
    let ready = match kind {
        "equation" => !state.model_ui.goal_expression.trim().is_empty(),
        "body_max_temperature" | "body_mean_temperature" => bodies.len() == 1,
        global if global.starts_with("global_") => true,
        _ => faces.len() == 1,
    };
    if button(ui, "multiphysics.model_add_goal", "multiphysics:model_add_goal", ready).clicked() {
        let name = if state.model_ui.goal_name.trim().is_empty() {
            format!("goal_{}", link.goals.len() + 1)
        } else {
            state.model_ui.goal_name.trim().to_string()
        };
        let mut params = json!({"id":study.id,"expected_revision":study.revision,"name":name,"kind":kind});
        match kind {
            "equation" => params["expression"] = json!(state.model_ui.goal_expression.trim()),
            "point_pressure" | "point_temperature" => params["point_mm"] = faces[0]["point"].clone(),
            "body_max_temperature" | "body_mean_temperature" => params["body"] = json!(bodies[0]),
            global if global.starts_with("global_") => {}
            _ => params["face"] = faces[0].clone(),
        }
        if model_command(app, state, study, "multiphysics.model.goal", params) {
            state.model_ui.goal_name.clear();
        }
    }
    ui.label(RichText::new(tr("multiphysics.model_goal_hint")).small().color(t.text_dim));
}

/// The study's status, read again only when the document changes.
fn status(app: &mut SolveApp, state: &mut State, study: &Study) -> Result<Value, String> {
    let key = std::sync::Arc::as_ptr(&app.session.doc) as usize;
    if let Some((k, id, reply)) = &state.model_ui.status
        && *k == key
        && *id == study.id
    {
        return reply.clone();
    }
    let reply = app.run("multiphysics.model.status", json!({"id":study.id}));
    state.model_ui.status = Some((key, study.id.clone(), reply.clone()));
    reply
}

/// Run a model command; keep its reply as the last build.
fn model_command(app: &mut SolveApp, state: &mut State, study: &Study, id: &str, params: Value) -> bool {
    match command(app, state, id, params) {
        Some(reply) => {
            state.model_ui.last_build = Some((study.id.clone(), reply));
            true
        }
        None => false,
    }
}

/// Build the case again when the model changed; `false` when it could not be built.
pub(super) fn ensure_current(app: &mut SolveApp, state: &mut State, study: &Study) -> bool {
    if study.model.is_none() {
        return true;
    }
    match status(app, state, study) {
        Ok(reply) if reply["current"] == true => true,
        _ => model_command(app, state, study, "multiphysics.model.build", json!({"id":study.id})),
    }
}

fn role_label(role: &Value) -> String {
    match role["kind"].as_str() {
        Some("inlet") => solidcraft_ui_text::format("multiphysics.model_inlet_value", &[("speed", &role["speed_m_per_s"].to_string())]),
        Some("outlet") => solidcraft_ui_text::format("multiphysics.model_outlet_value", &[("pressure", &role["pressure_m2_per_s2"].to_string())]),
        Some("wall") => tr("multiphysics.model_role_wall").to_string(),
        Some("symmetry") => tr("multiphysics.model_role_symmetry").to_string(),
        _ => String::new(),
    }
}

/// The faces in the design's selection, as `{body, point}` picks.
fn selected_faces(app: &SolveApp) -> Vec<Value> {
    app.session
        .selection
        .iter()
        .filter_map(|sel| match sel {
            solidcraft_engine::Sel::Face { body, point, .. } => Some(json!({"body":body,"point":[point.x, point.y, point.z]})),
            _ => None,
        })
        .collect()
}

/// The panel of a study of the loaded model.
pub(super) fn panel(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    let t = Tokens::get();
    let Some(link) = study.model.clone() else { return };
    ui.label(RichText::new(tr("multiphysics.model_hint")).small().color(t.text_dim));
    // Is the case the model as it is now?
    match status(app, state, study) {
        Ok(reply) if reply["current"] == true => {
            ui.colored_label(t.form_success, tr("multiphysics.model_current"));
        }
        Ok(reply) => {
            match reply["problem"].as_str() {
                Some(problem) => ui.colored_label(t.form_error, problem),
                None => ui.colored_label(t.warning, tr("multiphysics.model_changed")),
            };
        }
        Err(error) => {
            ui.colored_label(t.form_error, error);
        }
    }
    if button(ui, "multiphysics.model_rebuild", "multiphysics:model_rebuild", true).clicked() {
        model_command(app, state, study, "multiphysics.model.build", json!({"id":study.id}));
    }
    ui.separator();
    // Where the fluid is.
    ui.strong(tr("multiphysics.model_analysis"));
    let mut analysis = None;
    ui.horizontal(|ui| {
        for (kind, key) in [("internal", "multiphysics.model_internal"), ("external", "multiphysics.model_external")] {
            let on = serde_json::to_value(link.analysis).ok().as_ref().and_then(Value::as_str) == Some(kind);
            let response = crate::ui_kit::choice_chip(ui, &format!("multiphysics_model_{kind}"), tr(key), on);
            crate::scenario::publish_target(&format!("multiphysics:model_{kind}"), response.rect);
            if response.clicked() && !on {
                analysis = Some(kind);
            }
        }
    });
    ui.label(
        RichText::new(tr(if analysis.unwrap_or("") == "external" || link.analysis == solidcraft_engine::multiphysics::model::Analysis::External {
            "multiphysics.model_external_hint"
        } else {
            "multiphysics.model_internal_hint"
        }))
        .small()
        .color(t.text_dim),
    );
    if let Some(kind) = analysis {
        model_command(app, state, study, "multiphysics.model.set", json!({"id":study.id,"expected_revision":study.revision,"analysis":kind}));
        return;
    }
    // Margin and cell size.
    if state.model_ui.cell.is_empty() {
        state.model_ui.cell = link.cell_mm.map_or_else(String::new, |c| c.to_string());
    }
    if state.model_ui.margin.is_empty() {
        state.model_ui.margin = link.margin_mm.to_string();
    }
    crate::ui_kit::row(ui, tr("multiphysics.model_cell"), |ui| {
        crate::ui_kit::text_field(ui, &mut state.model_ui.cell, tr("multiphysics.model_cell_auto"))
    });
    if link.analysis == solidcraft_engine::multiphysics::model::Analysis::External {
        crate::ui_kit::row(ui, tr("multiphysics.model_margin"), |ui| crate::ui_kit::text_field(ui, &mut state.model_ui.margin, ""));
    }
    let mut cut = link.cut_cells;
    let response = crate::ui_kit::checkbox(ui, tr("multiphysics.model_cut_cells"), &mut cut);
    crate::scenario::publish_target("multiphysics:model_cut_cells", response.rect);
    if cut != link.cut_cells {
        model_command(app, state, study, "multiphysics.model.set", json!({"id":study.id,"expected_revision":study.revision,"cut_cells":cut}));
        return;
    }
    ui.label(RichText::new(tr("multiphysics.model_cut_cells_hint")).small().color(t.text_dim));
    let mut close = link.close_thin_walls;
    let response = crate::ui_kit::checkbox(ui, tr("multiphysics.model_thin_walls"), &mut close);
    crate::scenario::publish_target("multiphysics:model_thin_walls", response.rect);
    if close != link.close_thin_walls {
        model_command(
            app,
            state,
            study,
            "multiphysics.model.set",
            json!({"id":study.id,"expected_revision":study.revision,"close_thin_walls":close}),
        );
        return;
    }
    ui.label(RichText::new(tr("multiphysics.model_thin_walls_hint")).small().color(t.text_dim));
    if button(ui, "multiphysics.model_apply_grid", "multiphysics:model_apply_grid", true).clicked() {
        let cell = state.model_ui.cell.trim().parse::<f64>().unwrap_or(0.0);
        let mut params = json!({"id":study.id,"expected_revision":study.revision,"cell_mm":cell});
        if let Ok(margin) = state.model_ui.margin.trim().parse::<f64>() {
            params["margin_mm"] = json!(margin);
        }
        model_command(app, state, study, "multiphysics.model.set", params);
    }
    // What the last build gave.
    if let Some((id, reply)) = &state.model_ui.last_build
        && *id == study.id
    {
        if let Some(fraction) = reply["solid_fraction"].as_f64() {
            let regions = reply["fluid_regions"].as_u64().unwrap_or(0).to_string();
            ui.label(solidcraft_ui_text::format(
                "multiphysics.model_built",
                &[("percent", &format!("{:.0}", fraction * 100.0)), ("regions", &regions)],
            ));
        }
        for note in reply["notes"].as_array().into_iter().flatten().filter_map(Value::as_str) {
            ui.colored_label(t.warning, note);
        }
    }
    ui.separator();
    // Conditions on model faces.
    ui.strong(tr("multiphysics.model_conditions"));
    let mut remove = None;
    for opening in &link.openings {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {} · {}", opening.name, role_label(&serde_json::to_value(opening.role).unwrap_or_default()), opening.face.body));
            let response = remove_button(ui);
            crate::scenario::publish_target(&format!("multiphysics:model_remove:{}", opening.name), response.rect);
            if response.clicked() {
                remove = Some(opening.name.clone());
            }
        });
    }
    if link.openings.is_empty() {
        ui.label(RichText::new(tr("multiphysics.model_no_conditions")).small().color(t.text_dim));
    }
    if let Some(name) = remove {
        model_command(
            app,
            state,
            study,
            "multiphysics.model.condition",
            json!({"id":study.id,"expected_revision":study.revision,"name":name,"remove":true}),
        );
        return;
    }
    // Add one on the face picked in the view.
    ui.label(RichText::new(tr("multiphysics.model_pick_hint")).small().color(t.text_dim));
    let role = ROLES.get(state.model_ui.role).copied().unwrap_or("inlet");
    crate::ui_kit::row(ui, tr("multiphysics.name"), |ui| crate::ui_kit::text_field(ui, &mut state.model_ui.name, &format!("{role} 1")));
    crate::ui_kit::row(ui, tr("multiphysics.model_role"), |ui| {
        crate::ui_kit::dropdown_menu(ui, "multiphysics_model_role", tr(&format!("multiphysics.model_role_{role}")), |ui| {
            for (index, kind) in ROLES.iter().enumerate() {
                let response = crate::ui_kit::list_picker(
                    ui,
                    &format!("multiphysics_model_role_{kind}"),
                    tr(&format!("multiphysics.model_role_{kind}")),
                    index == state.model_ui.role,
                );
                crate::scenario::publish_target(&format!("multiphysics:model_role:{kind}"), response.rect);
                if response.clicked() {
                    state.model_ui.role = index;
                    ui.close();
                }
            }
        });
    });
    let value_key = match role {
        "inlet" => Some("multiphysics.model_speed"),
        "outlet" => Some("multiphysics.model_pressure"),
        _ => None,
    };
    if let Some(key) = value_key {
        crate::ui_kit::row(ui, tr(key), |ui| crate::ui_kit::text_field(ui, &mut state.model_ui.value, if role == "inlet" { "1" } else { "0" }));
    }
    let faces = selected_faces(app);
    let add = button(ui, "multiphysics.model_add_condition", "multiphysics:model_add_condition", faces.len() == 1);
    if faces.len() != 1 {
        ui.label(RichText::new(tr("multiphysics.model_pick_one")).small().color(t.text_dim));
    }
    if add.clicked()
        && let Some(face) = faces.clone().into_iter().next()
    {
        let typed = state.model_ui.value.trim().parse::<f64>().ok();
        let role_value = match role {
            "inlet" => json!({"kind":"inlet","speed_m_per_s":typed.unwrap_or(1.0)}),
            "outlet" => json!({"kind":"outlet","pressure_m2_per_s2":typed.unwrap_or(0.0)}),
            kind => json!({"kind":kind}),
        };
        let name = if state.model_ui.name.trim().is_empty() {
            let used = link.openings.iter().filter(|o| matches!(serde_json::to_value(o.role).ok(), Some(v) if v["kind"] == role)).count();
            format!("{role} {}", used + 1)
        } else {
            state.model_ui.name.trim().to_string()
        };
        let params = json!({"id":study.id,"expected_revision":study.revision,"name":name,"face":face,"role":role_value});
        if model_command(app, state, study, "multiphysics.model.condition", params) {
            state.model_ui.name.clear();
        }
    }
    ui.label(RichText::new(tr("multiphysics.model_inside_face_hint")).small().color(t.text_dim));
    ui.separator();
    // A fan on a face on the domain boundary, following its curve.
    ui.strong(tr("multiphysics.model_fan"));
    if let Some(fan) = &link.fan {
        let mut removed = false;
        ui.horizontal(|ui| {
            ui.label(format!("{} · {} · {} {}", fan.name, fan.face.body, fan.curve.len(), tr("multiphysics.model_fan_points")));
            let response = remove_button(ui);
            crate::scenario::publish_target("multiphysics:model_remove_fan", response.rect);
            removed = response.clicked();
        });
        if removed {
            model_command(app, state, study, "multiphysics.model.fan", json!({"id":study.id,"expected_revision":study.revision,"remove":true}));
            return;
        }
    }
    crate::ui_kit::row(ui, tr("multiphysics.model_fan_curve"), |ui| {
        field(ui, &mut state.model_ui.fan_curve, "0 4; 0.0004 3; 0.001 0", "multiphysics:model_fan_curve")
    });
    let curve: Option<Vec<[f64; 2]>> = state
        .model_ui
        .fan_curve
        .split(';')
        .filter(|point| !point.trim().is_empty())
        .map(|point| {
            let values: Vec<f64> = point.split_whitespace().filter_map(|v| v.parse().ok()).collect();
            match values.as_slice() {
                [q, p] => Some([*q, *p]),
                _ => None,
            }
        })
        .collect::<Option<Vec<_>>>()
        .filter(|points| points.len() >= 2);
    if button(ui, "multiphysics.model_add_fan", "multiphysics:model_add_fan", faces.len() == 1 && curve.is_some()).clicked()
        && let (Some(face), Some(curve)) = (faces.first(), curve)
    {
        model_command(
            app,
            state,
            study,
            "multiphysics.model.fan",
            json!({"id":study.id,"expected_revision":study.revision,"face":face,"curve":curve}),
        );
        return;
    }
    ui.label(RichText::new(tr("multiphysics.model_fan_hint")).small().color(t.text_dim));
    ui.separator();
    // Leak tracking: the fluid path between two picked faces.
    ui.label(RichText::new(tr("multiphysics.model_leak_hint")).small().color(t.text_dim));
    if button(ui, "multiphysics.model_leak", "multiphysics:model_leak", faces.len() == 2).clicked() {
        let params = json!({"id":study.id,"from":faces[0],"to":faces[1]});
        if let Some(reply) = command(app, state, "multiphysics.model.leak", params) {
            let path: Vec<[f64; 3]> = serde_json::from_value(reply["path_mm"].clone()).unwrap_or_default();
            state.model_ui.leak = Some((study.id.clone(), reply["connected"] == true, path));
        }
    }
    if let Some((id, connected, _)) = &state.model_ui.leak
        && *id == study.id
    {
        if *connected {
            ui.colored_label(t.warning, tr("multiphysics.model_leak_found"));
        } else {
            ui.colored_label(t.form_success, tr("multiphysics.model_leak_sealed"));
        }
    }
    ui.separator();
    environment_section(app, state, study, &link, ui);
}

/// Where the model stands and what it breathes: the air around it (external studies), the dust
/// in that air, and the year of conditions it works through, with the year's summary.
fn environment_section(
    app: &mut SolveApp,
    state: &mut State,
    study: &Study,
    link: &solidcraft_engine::multiphysics::model::ModelLink,
    ui: &mut egui::Ui,
) {
    let t = Tokens::get();
    let base = || json!({"id":study.id,"expected_revision":study.revision});
    ui.strong(tr("multiphysics.model_environment"));
    ui.label(RichText::new(tr("multiphysics.model_environment_hint")).small().color(t.text_dim));
    // The air around the model.
    if link.analysis == solidcraft_engine::multiphysics::model::Analysis::External {
        if state.model_ui.breeze.is_empty() {
            state.model_ui.breeze = link.surroundings.map_or_else(|| "0.1".into(), |s| s.breeze_m_per_s.to_string());
        }
        let mut apply = false;
        crate::ui_kit::row(ui, tr("multiphysics.model_breeze"), |ui| {
            field(ui, &mut state.model_ui.breeze, "0.1", "multiphysics:model_breeze");
            let response = crate::ui_kit::button(ui, tr("multiphysics.model_apply"), crate::ui_kit::Button::Small);
            crate::scenario::publish_target("multiphysics:model_apply_breeze", response.rect);
            apply = response.clicked();
        });
        let mut floor = link.surroundings.is_none_or(|s| s.floor);
        let floor_changed = crate::ui_kit::checkbox(ui, tr("multiphysics.model_floor"), &mut floor).changed();
        ui.label(RichText::new(tr("multiphysics.model_breeze_hint")).small().color(t.text_dim));
        if (apply || floor_changed)
            && let Ok(breeze) = state.model_ui.breeze.trim().parse::<f64>()
        {
            let mut params = base();
            params["breeze_m_per_s"] = json!(breeze);
            params["floor"] = json!(floor);
            model_command(app, state, study, "multiphysics.model.surroundings", params);
            return;
        }
    }
    // Dust.
    if state.model_ui.dust_mg.is_empty() {
        state.model_ui.dust_mg = link.dust.as_ref().map_or_else(|| "0.5".into(), |d| d.concentration_mg_per_m3.to_string());
    }
    if state.model_ui.dust_inset.is_empty() {
        state.model_ui.dust_inset = "2".into();
    }
    let mut apply_dust = false;
    crate::ui_kit::row(ui, tr("multiphysics.model_dust"), |ui| {
        field(ui, &mut state.model_ui.dust_mg, "0.5", "multiphysics:model_dust");
        let response = crate::ui_kit::button(ui, tr("multiphysics.model_apply"), crate::ui_kit::Button::Small);
        crate::scenario::publish_target("multiphysics:model_apply_dust", response.rect);
        apply_dust = response.clicked();
    });
    if apply_dust && let Ok(mg) = state.model_ui.dust_mg.trim().parse::<f64>() {
        let mut params = base();
        params["concentration_mg_per_m3"] = json!(mg);
        model_command(app, state, study, "multiphysics.model.dust", params);
        return;
    }
    if let Some(dust) = &link.dust {
        let entry = match &dust.entry {
            solidcraft_engine::multiphysics::environment::DustEntry::Air => tr("multiphysics.model_dust_from_air").to_string(),
            solidcraft_engine::multiphysics::environment::DustEntry::Box { low_mm, high_mm } => solidcraft_ui_text::format(
                "multiphysics.model_dust_from_box",
                &[
                    ("low", &format!("{:.1}, {:.1}, {:.1}", low_mm[0], low_mm[1], low_mm[2])),
                    ("high", &format!("{:.1}, {:.1}, {:.1}", high_mm[0], high_mm[1], high_mm[2])),
                ],
            ),
        };
        ui.label(RichText::new(entry).small());
        // Count the dust that comes into a picked body (an enclosure) through its openings.
        let picked = selected_bodies(app);
        crate::ui_kit::row(ui, tr("multiphysics.model_dust_inset"), |ui| {
            field(ui, &mut state.model_ui.dust_inset, "2", "multiphysics:model_dust_inset");
        });
        if button(ui, "multiphysics.model_dust_inside", "multiphysics:model_dust_inside", picked.len() == 1).clicked()
            && let Some(body) = picked.first()
        {
            let mut params = base();
            params["inside_of"] = json!(body);
            params["inset_mm"] = json!(state.model_ui.dust_inset.trim().parse::<f64>().unwrap_or(0.0));
            model_command(app, state, study, "multiphysics.model.dust", params);
            return;
        }
        let mut remove = false;
        ui.horizontal(|ui| {
            ui.label(RichText::new(tr("multiphysics.model_dust_hint")).small().color(t.text_dim));
            remove = remove_button(ui).clicked();
        });
        if remove {
            let mut params = base();
            params["remove"] = json!(true);
            model_command(app, state, study, "multiphysics.model.dust", params);
            return;
        }
    }
    ui.separator();
    // The year: four seasons, each a night, morning, afternoon and evening.
    ui.strong(tr("multiphysics.model_climate"));
    ui.label(RichText::new(tr("multiphysics.model_climate_hint")).small().color(t.text_dim));
    if state.model_ui.seasons.len() != SEASONS.len() {
        state.model_ui.seasons = SEASONS.iter().map(|(_, v)| v.map(|x| x.to_string())).collect();
    }
    egui::Grid::new(("multiphysics_seasons", &study.id)).num_columns(5).show(ui, |ui| {
        for key in [
            "multiphysics.model_season",
            "multiphysics.model_season_low",
            "multiphysics.model_season_high",
            "multiphysics.model_season_dew",
            "multiphysics.model_season_dust",
        ] {
            ui.label(RichText::new(tr(key)).small());
        }
        ui.end_row();
        for (row, (name, _)) in state.model_ui.seasons.iter_mut().zip(SEASONS) {
            ui.label(name);
            for (k, value) in row.iter_mut().enumerate() {
                field(ui, value, "", &format!("multiphysics:model_season:{name}:{k}"));
            }
            ui.end_row();
        }
    });
    if button(ui, "multiphysics.model_climate_set", "multiphysics:model_climate_set", true).clicked() {
        let seasons: Option<Vec<Value>> = state
            .model_ui
            .seasons
            .iter()
            .zip(SEASONS)
            .map(|(row, (name, _))| {
                let v: Option<Vec<f64>> = row.iter().map(|x| x.trim().parse::<f64>().ok()).collect();
                v.map(|v| json!({"name": name, "low_c": v[0], "high_c": v[1], "dew_point_c": v[2], "dust_mg_per_m3": v[3]}))
            })
            .collect();
        if let Some(seasons) = seasons {
            let mut params = base();
            params["seasons"] = json!(seasons);
            model_command(app, state, study, "multiphysics.model.climate", params);
            return;
        }
    }
    let Some(climate) = &link.climate else { return };
    egui::Grid::new(("multiphysics_conditions", &study.id)).num_columns(4).striped(true).show(ui, |ui| {
        for c in &climate.conditions {
            ui.label(&c.name);
            ui.label(format!("{:.0} h", c.hours_per_year));
            ui.label(format!("{:.1} °C", c.ambient_c));
            ui.label(format!("{:.0} % RH", c.relative_humidity_pct()));
            ui.end_row();
        }
    });
    if button(ui, "multiphysics.model_climate_studies", "multiphysics:model_climate_studies", true).clicked() {
        command(app, state, "multiphysics.model.climate_studies", json!({"id":study.id}));
        return;
    }
    ui.label(RichText::new(tr("multiphysics.model_climate_studies_hint")).small().color(t.text_dim));
    crate::ui_kit::row(ui, tr("multiphysics.model_climate_folder"), |ui| {
        field(ui, &mut state.model_ui.climate_folder, "", "multiphysics:model_climate_folder");
    });
    if button(ui, "multiphysics.model_climate_summary", "multiphysics:model_climate_summary", true).clicked() {
        let mut params = json!({"id":study.id});
        let folder = state.model_ui.climate_folder.trim();
        if !folder.is_empty() {
            params["directory"] = json!(folder);
        }
        if let Some(reply) = command(app, state, "multiphysics.result.climate", params) {
            state.model_ui.climate = Some((study.id.clone(), reply));
        }
    }
    if let Some((id, summary)) = &state.model_ui.climate
        && *id == study.id
    {
        ui.label(solidcraft_ui_text::format(
            "multiphysics.model_climate_totals",
            &[
                ("in", &format!("{:.3}", summary["dust_in_g_per_year"].as_f64().unwrap_or(0.0))),
                ("deposited", &format!("{:.3}", summary["deposited_g_per_year"].as_f64().unwrap_or(0.0))),
                ("hours", &format!("{:.0}", summary["hours_covered"].as_f64().unwrap_or(0.0))),
            ],
        ));
        for missing in summary["missing"].as_array().into_iter().flatten().filter_map(Value::as_str) {
            ui.colored_label(t.warning, solidcraft_ui_text::format("multiphysics.model_climate_missing", &[("condition", missing)]));
        }
        egui::Grid::new(("multiphysics_climate_parts", &study.id)).num_columns(4).striped(true).show(ui, |ui| {
            for key in [
                "multiphysics.model_climate_part",
                "multiphysics.model_climate_dust",
                "multiphysics.model_climate_humid",
                "multiphysics.model_climate_wet",
            ] {
                ui.label(RichText::new(tr(key)).small());
            }
            ui.end_row();
            for part in summary["parts"].as_array().into_iter().flatten().take(24) {
                ui.label(part["body"].as_str().unwrap_or_default());
                ui.label(format!("{:.4} g", part["dust_g_per_year"].as_f64().unwrap_or(0.0)));
                ui.label(format!("{:.0} h", part["humid_dust_hours"].as_f64().unwrap_or(0.0)));
                ui.label(format!("{:.0} h", part["condensing_hours"].as_f64().unwrap_or(0.0)));
                ui.end_row();
            }
        });
    }
}

/// Draw the computational domain and the conditions over the model in the 3-D view.
pub(super) fn overlay(app: &mut SolveApp, state: &mut State, study: &Study, ui: &egui::Ui, rect: egui::Rect) {
    let Ok(reply) = status(app, state, study) else { return };
    let t = Tokens::get();
    let proj = crate::viewport::projection(app, rect);
    let painter = ui.painter_at(rect);
    let to_screen = |p: [f64; 3]| proj.to_screen(solidcore_record::Vec3::new(p[0], p[1], p[2]));
    let domain: Option<[[f64; 3]; 2]> = serde_json::from_value(reply["domain_mm"].clone()).ok();
    if let Some([lo, hi]) = domain {
        let corner =
            |i: usize| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }];
        let stroke = Stroke::new(1.5, t.accent);
        for a in 0..8usize {
            for bit in [1usize, 2, 4] {
                if a & bit == 0
                    && let (Some(p), Some(q)) = (to_screen(corner(a)), to_screen(corner(a | bit)))
                {
                    painter.line_segment([p, q], stroke);
                }
            }
        }
        if let Some(p) = to_screen(hi) {
            painter.text(p, Align2::LEFT_BOTTOM, tr("multiphysics.model_domain"), FontId::proportional(t.font_small), t.accent);
        }
    }
    // The grid's planes on the domain's three low faces.
    if state.model_ui.show_grid
        && let Some([lo, hi]) = domain
        && let Ok(nodes) = solidcraft_engine::multiphysics::case_scene::nodes(&study.case)
    {
        let stroke = Stroke::new(1.0, t.text_dim.gamma_multiply(0.6));
        let mut lines = 0usize;
        for face in 0..3 {
            for along in (0..3).filter(|a| *a != face) {
                let across = 3 - face - along;
                for x in &nodes[along] {
                    let mut p = lo;
                    p[along] = lo[along] + x * 1000.0;
                    let mut q = p;
                    q[across] = hi[across];
                    if let (Some(a), Some(b)) = (to_screen(p), to_screen(q)) {
                        painter.line_segment([a, b], stroke);
                        lines += 1;
                    }
                }
            }
        }
        crate::scenario::publish_count("multiphysics:model_grid_lines", lines as f64);
    }
    // The last leak path, through the fluid.
    if let Some((id, true, path)) = &state.model_ui.leak
        && *id == study.id
    {
        let points: Vec<egui::Pos2> = path.iter().filter_map(|p| to_screen(*p)).collect();
        if points.len() >= 2 {
            painter.add(egui::Shape::line(points, Stroke::new(2.5, t.error)));
        }
    }
    let size = domain.map_or(10.0, |[lo, hi]| (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max) * 0.12);
    let openings = study.model.as_ref().map(|link| link.openings.clone()).unwrap_or_default();
    for (mark, opening) in reply["faces"].as_array().into_iter().flatten().zip(&openings) {
        let (Some(point), normal) =
            (serde_json::from_value::<[f64; 3]>(mark["point_mm"].clone()).ok(), serde_json::from_value::<[f64; 3]>(mark["normal"].clone()).ok())
        else {
            continue;
        };
        let role = serde_json::to_value(opening.role).unwrap_or_default();
        let color: Color32 = match role["kind"].as_str() {
            Some("inlet") => t.accent,
            Some("outlet") => t.form_success,
            _ => t.text_dim,
        };
        let Some(at) = to_screen(point) else { continue };
        painter.circle_filled(at, 4.0, color);
        if let Some(n) = normal {
            // Inflow points into the fluid (against the face's outward normal), outflow out of it.
            let sign = if role["kind"] == "inlet" { -1.0 } else { 1.0 };
            let (tail, head) =
                if sign < 0.0 { ([0, 1, 2].map(|a| point[a] + n[a] * size), point) } else { (point, [0, 1, 2].map(|a| point[a] + n[a] * size)) };
            if let (Some(p), Some(q)) = (to_screen(tail), to_screen(head)) {
                painter.arrow(p, q - p, Stroke::new(2.0, color));
            }
        }
        painter.text(
            at + egui::vec2(6.0, -6.0),
            Align2::LEFT_BOTTOM,
            format!("{} · {}", opening.name, role_label(&role)),
            FontId::proportional(t.font_small),
            color,
        );
    }
}

/// The fluids offered: their name key, kinematic viscosity (m²/s) and density (kg/m³), at 20 °C.
const FLUIDS: [(&str, f64, f64); 2] = [("multiphysics.model_fluid_air", 1.5e-5, 1.2), ("multiphysics.model_fluid_water", 1.0e-6, 998.0)];

/// The fluid and how long the run simulates. A study made before these existed runs with its
/// template's until one is picked. True when a change was sent (the study is built again).
fn fluid_section(
    app: &mut SolveApp,
    state: &mut State,
    study: &solidcraft_engine::multiphysics::Study,
    link: &solidcraft_engine::multiphysics::model::ModelLink,
    ui: &mut egui::Ui,
) -> bool {
    use solidcraft_engine::multiphysics::model::ModelFlow;
    let t = Tokens::get();
    ui.separator();
    ui.strong(tr("multiphysics.model_fluid"));
    let current = link.flow;
    let mut wanted = current;
    ui.horizontal_wrapped(|ui| {
        for (index, (key, nu, rho)) in FLUIDS.iter().enumerate() {
            let on = current.is_some_and(|f| f.kinematic_viscosity_m2_per_s == *nu && f.density_kg_per_m3 == *rho);
            let response = crate::ui_kit::choice_chip(ui, key, tr(key), on);
            crate::scenario::publish_target(&format!("multiphysics:model_fluid:{index}"), response.rect);
            if response.clicked() && !on {
                let base = current.unwrap_or(ModelFlow::AIR);
                wanted = Some(ModelFlow { kinematic_viscosity_m2_per_s: *nu, density_kg_per_m3: *rho, ..base });
            }
        }
    });
    match current {
        None => {
            ui.label(RichText::new(tr("multiphysics.model_fluid_template")).small().color(t.warning));
        }
        Some(flow) => {
            ui.label(
                RichText::new(solidcraft_ui_text::format(
                    "multiphysics.model_fluid_values",
                    &[("nu", &format!("{:.3e}", flow.kinematic_viscosity_m2_per_s)), ("rho", &format!("{}", flow.density_kg_per_m3))],
                ))
                .small()
                .color(t.text_dim),
            );
            // The end time as typed: applied when the field is left with a new value.
            if state.model_ui.end_time.as_ref().is_none_or(|(id, _)| *id != study.id) {
                state.model_ui.end_time = Some((study.id.clone(), flow.end_time_s.to_string()));
            }
            let mut text = state.model_ui.end_time.as_ref().map(|(_, t)| t.clone()).unwrap_or_default();
            let (enter, response) =
                crate::ui_kit::row(ui, tr("multiphysics.model_end_time"), |ui| crate::ui_kit::text_field_reserving(ui, &mut text, "s", 0.0));
            crate::scenario::publish_target("multiphysics:model_end_time", response.rect);
            let typed = text.trim().parse::<f64>().ok().filter(|v| v.is_finite() && *v > 0.0);
            let left = enter || response.lost_focus();
            state.model_ui.end_time = Some((study.id.clone(), text));
            let mut settle = flow.converge_tolerance.is_some();
            let toggled = crate::ui_kit::checkbox(ui, tr("multiphysics.model_settle"), &mut settle);
            crate::scenario::publish_target("multiphysics:model_settle", toggled.rect);
            let end = match typed {
                Some(v) if left => v,
                _ => flow.end_time_s,
            };
            if toggled.changed() || (end - flow.end_time_s).abs() > 1e-12 {
                wanted = Some(ModelFlow { end_time_s: end, converge_tolerance: settle.then_some(flow.converge_tolerance.unwrap_or(1e-3)), ..flow });
                state.model_ui.end_time = None;
            }
            ui.label(RichText::new(tr("multiphysics.model_end_time_hint")).small().color(t.text_dim));
        }
    }
    if wanted != current
        && let Some(flow) = wanted
    {
        model_command(app, state, study, "multiphysics.model.set", json!({"id":study.id,"expected_revision":study.revision,"flow":flow}));
        return true;
    }
    false
}

/// Bodies taken as porous zones (a catalyst brick, a filter): the picked body becomes one with
/// the typed resistances. True when a change was sent (the study is built again).
fn porous_section(
    app: &mut SolveApp,
    state: &mut State,
    study: &solidcraft_engine::multiphysics::Study,
    link: &solidcraft_engine::multiphysics::model::ModelLink,
    ui: &mut egui::Ui,
) -> bool {
    let t = Tokens::get();
    ui.separator();
    ui.strong(tr("multiphysics.model_porous"));
    let mut list = link.porous.clone();
    let mut changed = false;
    for (index, p) in link.porous.iter().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("{} · D {:e} 1/m² · C {} 1/m", p.body, p.viscous_resistance_per_m2, p.inertial_resistance_per_m));
            let response = remove_button(ui);
            crate::scenario::publish_target(&format!("multiphysics:model_porous_remove:{index}"), response.rect);
            if response.clicked() {
                list.retain(|q| q.body != p.body);
                changed = true;
            }
        });
    }
    let bodies = selected_bodies(app);
    crate::ui_kit::row(ui, tr("multiphysics.model_porous_d"), |ui| crate::ui_kit::text_field(ui, &mut state.model_ui.porous_d, "1e6"));
    crate::ui_kit::row(ui, tr("multiphysics.model_porous_c"), |ui| crate::ui_kit::text_field(ui, &mut state.model_ui.porous_c, "0"));
    let d = state.model_ui.porous_d.trim().parse::<f64>().ok().filter(|v| v.is_finite() && *v >= 0.0);
    let c = if state.model_ui.porous_c.trim().is_empty() {
        Some(0.0)
    } else {
        state.model_ui.porous_c.trim().parse::<f64>().ok().filter(|v| v.is_finite() && *v >= 0.0)
    };
    let ready = bodies.len() == 1 && d.is_some() && c.is_some();
    if button(ui, "multiphysics.model_porous_make", "multiphysics:model_porous_make", ready).clicked()
        && let (Some(body), Some(d), Some(c)) = (bodies.first(), d, c)
    {
        list.retain(|q| &q.body != body);
        list.push(solidcraft_engine::multiphysics::model::ModelPorous {
            body: body.clone(),
            viscous_resistance_per_m2: d,
            inertial_resistance_per_m: c,
        });
        changed = true;
    }
    ui.label(RichText::new(tr("multiphysics.model_porous_hint")).small().color(t.text_dim));
    if changed {
        model_command(app, state, study, "multiphysics.model.set", json!({"id":study.id,"expected_revision":study.revision,"porous":list}));
    }
    changed
}
