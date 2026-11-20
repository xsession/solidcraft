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
    /// The last leak search: the study, whether fluid connects the faces, and the path (mm).
    leak: Option<(String, bool, Vec<[f64; 3]>)>,
    /// Local mesh: the split factor offered, 2 to 4 (0 means 2).
    refine_factor: usize,
    /// The goal being added: its kind (index in [`GOAL_KINDS`]), name and equation.
    goal_kind: usize,
    goal_name: String,
    goal_expression: String,
}

/// The goal kinds the panel offers, as `multiphysics.model.goal` names them.
const GOAL_KINDS: [&str; 4] = ["face_volume_flow", "face_mean_pressure", "point_pressure", "equation"];

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
    ui.strong(tr("multiphysics.model_local_mesh"));
    let mut remove: Option<(&str, String)> = None;
    for r in &link.refinements {
        ui.horizontal(|ui| {
            let what = match &r.target {
                solidcraft_engine::multiphysics::model::RefineTarget::Body { body } => body.clone(),
                solidcraft_engine::multiphysics::model::RefineTarget::Face { face } => format!("{} · {}", face.body, face.face),
            };
            ui.label(format!("{} · ×{} · {what}", r.name, r.factor));
            let response = crate::ui_kit::button(ui, "✕", crate::ui_kit::Button::Small).backing_on_hover_text(tr("multiphysics.model_remove"));
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
    }
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
            let response = crate::ui_kit::button(ui, "✕", crate::ui_kit::Button::Small).backing_on_hover_text(tr("multiphysics.model_remove"));
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
    let kind = GOAL_KINDS.get(state.model_ui.goal_kind).copied().unwrap_or(GOAL_KINDS[0]);
    crate::ui_kit::row(ui, tr("multiphysics.name"), |ui| {
        crate::ui_kit::text_field(ui, &mut state.model_ui.goal_name, &format!("goal_{}", link.goals.len() + 1))
    });
    crate::ui_kit::row(ui, tr("multiphysics.model_goal_kind"), |ui| {
        crate::ui_kit::dropdown_menu(ui, "multiphysics_model_goal_kind", tr(&format!("multiphysics.model_goal_{kind}")), |ui| {
            for (index, k) in GOAL_KINDS.iter().enumerate() {
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
    });
    if kind == "equation" {
        crate::ui_kit::row(ui, tr("multiphysics.model_goal_expression"), |ui| {
            crate::ui_kit::text_field(ui, &mut state.model_ui.goal_expression, "a - b")
        });
    }
    let ready = match kind {
        "equation" => !state.model_ui.goal_expression.trim().is_empty(),
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
            "point_pressure" => params["point_mm"] = faces[0]["point"].clone(),
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
            let response = crate::ui_kit::list_picker(ui, &format!("multiphysics_model_{kind}"), tr(key), on);
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
            let response = crate::ui_kit::button(ui, "✕", crate::ui_kit::Button::Small).backing_on_hover_text(tr("multiphysics.model_remove"));
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
