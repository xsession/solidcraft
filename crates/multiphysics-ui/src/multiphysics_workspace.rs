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

#[path = "multiphysics_model_ui.rs"]
mod model_ui;
#[cfg(not(target_arch = "wasm32"))]
#[path = "multiphysics_on_model.rs"]
mod on_model;
#[path = "multiphysics_parts_ui.rs"]
mod parts_ui;
#[path = "multiphysics_shell.rs"]
mod shell;

pub use solidcraft_ui_egui::multiphysics_workspace::{active, can_enter, switch, workspace_name};
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
// The slice and colour fields are read by the result views, which are native only.
#[cfg_attr(target_arch = "wasm32", allow(dead_code))]
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
    /// What the cut plot's plane goes through: the typed layer, the picked part or the field's
    /// maximum, worked out again for every run, field and pick.
    #[cfg(not(target_arch = "wasm32"))]
    slice_anchor: SliceAnchor,
    /// The picked part's centre in the design (mm), by document and body, for a plot through it.
    #[cfg(not(target_arch = "wasm32"))]
    anchor_part: Option<(usize, String, Option<[f64; 3]>)>,
    /// The colour map as plots use it: reversed, banded, out-of-range marking.
    color_scale: post::ColorScale,
    /// How the plots' colour range is chosen; `locked_range` holds the Locked values.
    range_mode: post::RangeMode,
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
    /// For a study of the loaded model: the cut plot drawn over the model in the 3-D view.
    #[cfg(not(target_arch = "wasm32"))]
    on_model: bool,
    /// Where the study's domain origin is in the design (mm), read once per document.
    #[cfg(not(target_arch = "wasm32"))]
    model_origin: Option<(usize, String, Option<[f64; 3]>)>,
    /// The Surface Plots node is shown: a study of the model colours its wetted faces.
    #[cfg(not(target_arch = "wasm32"))]
    model_surfaces: bool,
    /// A flow tool's node is shown (its heading key): a study of the model draws the tool's
    /// trajectories or isosurface on the model instead of a cut plot.
    #[cfg(not(target_arch = "wasm32"))]
    model_flow: Option<&'static str>,
    /// Velocity arrows on the cut plot on the model.
    #[cfg(not(target_arch = "wasm32"))]
    model_vectors: bool,
    /// The study whose model the view last framed, and the flow lines it last scrolled to.
    #[cfg(not(target_arch = "wasm32"))]
    model_framed: Option<String>,
    #[cfg(not(target_arch = "wasm32"))]
    model_flow_seen: usize,
    /// The flow tool is in the docked property panel: its lines go on the graphics window's
    /// model (one 3-D model view a frame), on this colour range.
    #[cfg(not(target_arch = "wasm32"))]
    model_flow_docked: bool,
    #[cfg(not(target_arch = "wasm32"))]
    model_flow_range: (f64, f64),
    /// The model's faces sorted into wetted and dry for the shown result.
    #[cfg(not(target_arch = "wasm32"))]
    model_surface: Option<std::sync::Arc<on_model::Surface>>,
    /// The run whose default field file was chosen when it was first shown.
    defaulted_run: Option<String>,
    /// Where the user last probed the result on the model.
    #[cfg(not(target_arch = "wasm32"))]
    model_probe: Option<on_model::Probe>,
    #[cfg(not(target_arch = "wasm32"))]
    field_3d: crate::multiphysics_field_view::State3d,
    error: Option<String>,
    notice: Option<String>,
}
/// What a cut plot's plane goes through.
#[cfg(not(target_arch = "wasm32"))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum SliceAnchor {
    /// The layer typed in.
    #[default]
    Layer,
    /// The centre of the part picked in the design (a study of the model).
    Part,
    /// The cell where the plotted field is largest (the hottest point, the fastest flow).
    Maximum,
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
            #[cfg(not(target_arch = "wasm32"))]
            slice_anchor: SliceAnchor::Layer,
            #[cfg(not(target_arch = "wasm32"))]
            anchor_part: None,
            color_scale: post::ColorScale::default(),
            range_mode: post::RangeMode::default(),
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
            on_model: false,
            #[cfg(not(target_arch = "wasm32"))]
            model_origin: None,
            #[cfg(not(target_arch = "wasm32"))]
            model_surfaces: false,
            #[cfg(not(target_arch = "wasm32"))]
            model_flow: None,
            #[cfg(not(target_arch = "wasm32"))]
            model_vectors: false,
            #[cfg(not(target_arch = "wasm32"))]
            model_framed: None,
            #[cfg(not(target_arch = "wasm32"))]
            model_flow_seen: 0,
            #[cfg(not(target_arch = "wasm32"))]
            model_flow_docked: false,
            #[cfg(not(target_arch = "wasm32"))]
            model_flow_range: (0.0, 1.0),
            #[cfg(not(target_arch = "wasm32"))]
            model_surface: None,
            defaulted_run: None,
            #[cfg(not(target_arch = "wasm32"))]
            model_probe: None,
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
    /// Trajectory, particle, isosurface and surface tools.
    flow: FlowTool,
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
    feature_name: String,
    replace_fan_patch_conditions: bool,
    feature_owner: Option<(u64, String, &'static str)>,
    /// The goal kind Add creates in the Goals panel (`authoring::FeatureKind`, snake case).
    goal_kind: &'static str,
    /// The 3-D case view in the graphics window and what is picked in it.
    case_view: crate::multiphysics_case_view::CaseView,
    /// Values typed in the panel for the picked parts and faces.
    case_tool: CaseTool,
    /// Add from Components: the model's bodies and the import choices.
    components: ComponentTool,
    /// Per project, the CAD point (mm) at the domain origin, from the last fitted import.
    /// Kept for the session only: the solver case has no place for it.
    cad_origin: BTreeMap<String, [f64; 3]>,
    /// A study of the loaded model: its status, last build and the condition being added.
    model_ui: model_ui::ModelUi,
    /// A study of the model's Parts folder: its component table and the picked part selection.
    parts_ui: parts_ui::PartsUi,
}

#[derive(Default)]
struct ComponentTool {
    /// The model's body names, read when the panel opens or on Refresh.
    bodies: Option<Result<Vec<String>, String>>,
    chosen: std::collections::BTreeSet<String>,
    representation: solidcraft_engine::multiphysics::case_scene::Representation,
    keep_domain: bool,
    margin: String,
    cell: String,
    replace: bool,
    notes: Vec<String>,
}

const REPRESENTATIONS: [(solidcraft_engine::multiphysics::case_scene::Representation, &str); 3] = [
    (solidcraft_engine::multiphysics::case_scene::Representation::Auto, "multiphysics.components_auto"),
    (solidcraft_engine::multiphysics::case_scene::Representation::Obstacle, "multiphysics.components_obstacle"),
    (solidcraft_engine::multiphysics::case_scene::Representation::Surface, "multiphysics.components_surface"),
];

/// The bodies picked in the CAD model.
fn cad_bodies(app: &SolveApp) -> Vec<String> {
    app.session
        .selection
        .iter()
        .filter_map(|sel| match sel {
            solidcraft_engine::Sel::Body { name } => Some(name.clone()),
            solidcraft_engine::Sel::Face { body, .. } => Some(body.clone()),
            _ => None,
        })
        .collect()
}

/// Points (mm) of the faces picked in the CAD model.
fn cad_face_points(app: &SolveApp) -> Vec<[f64; 3]> {
    app.session
        .selection
        .iter()
        .filter_map(|sel| match sel {
            solidcraft_engine::Sel::Face { point, .. } => Some([point.x, point.y, point.z]),
            _ => None,
        })
        .collect()
}

/// Add from Components: turn model bodies into case parts (grid blocks or surface bodies).
fn components_panel(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    use solidcraft_engine::multiphysics::case_scene::Representation;
    let t = Tokens::get();
    ui.small(tr("multiphysics.components_hint"));
    if state.components.bodies.is_none() || button(ui, "multiphysics.components_refresh", "multiphysics:components_refresh", true).clicked() {
        let listed = app
            .run("multiphysics.cad.bodies", json!({}))
            .map(|reply| reply["bodies"].as_array().into_iter().flatten().filter_map(|b| b.as_str().map(str::to_string)).collect::<Vec<_>>());
        if let Ok(bodies) = &listed {
            state.components.chosen.retain(|name| bodies.contains(name));
        }
        state.components.bodies = Some(listed);
    }
    let picked = cad_bodies(app);
    if !picked.is_empty() && button(ui, "multiphysics.components_use_selection", "multiphysics:components_use_selection", true).clicked() {
        state.components.chosen = picked.into_iter().collect();
    }
    match state.components.bodies.clone() {
        Some(Ok(bodies)) if bodies.is_empty() => {
            ui.label(RichText::new(tr("multiphysics.components_none")).color(t.text_dim));
        }
        Some(Ok(bodies)) => {
            egui::ScrollArea::vertical().id_salt("multiphysics_components").max_height(180.0).show(ui, |ui| {
                for name in bodies {
                    let mut on = state.components.chosen.contains(&name);
                    let response = crate::ui_kit::checkbox(ui, &name, &mut on);
                    crate::scenario::publish_target(&format!("multiphysics:component:{name}"), response.rect);
                    if response.changed() {
                        if on {
                            state.components.chosen.insert(name);
                        } else {
                            state.components.chosen.remove(&name);
                        }
                    }
                }
            });
        }
        Some(Err(error)) => {
            ui.label(RichText::new(error).color(t.form_error));
        }
        None => {}
    }
    let current = REPRESENTATIONS.iter().find(|(r, _)| *r == state.components.representation).map_or("", |(_, key)| key);
    let r = crate::ui_kit::row(ui, tr("multiphysics.components_as"), |ui| {
        crate::ui_kit::dropdown_menu(ui, "multiphysics_components_as", tr(current), |ui| {
            for (representation, key) in REPRESENTATIONS {
                let response = crate::ui_kit::list_picker(
                    ui,
                    &format!("multiphysics_components_as_{key}"),
                    tr(key),
                    representation == state.components.representation,
                );
                if response.clicked() {
                    state.components.representation = representation;
                    ui.close();
                }
            }
        })
    });
    crate::scenario::publish_target("multiphysics:components_as", r.rect);
    let origin = state.cad_origin.get(&study.id).copied();
    // Without a CAD origin the first import has to size the domain around the bodies.
    let fit = !state.components.keep_domain || origin.is_none();
    let r = ui
        .add_enabled_ui(origin.is_some(), |ui| {
            crate::ui_kit::checkbox(ui, tr("multiphysics.components_keep_domain"), &mut state.components.keep_domain)
        })
        .inner;
    crate::scenario::publish_target("multiphysics:components_keep_domain", r.rect);
    match origin {
        Some(o) => {
            ui.small(solidcraft_ui_text::format("multiphysics.components_origin", &[("origin", &format!("{:.2}, {:.2}, {:.2}", o[0], o[1], o[2]))]))
        }
        None => ui.small(tr("multiphysics.components_no_origin")),
    };
    if fit {
        if state.components.margin.is_empty() {
            state.components.margin = "10".into();
        }
        let r = crate::ui_kit::row(ui, tr("multiphysics.components_margin"), |ui| {
            crate::ui_kit::bounded_cell_editor(ui, ui.scope_id().with("components_margin"), &mut state.components.margin, "", false, 16)
        });
        crate::scenario::publish_target("multiphysics:components_margin", r.rect);
        let r = crate::ui_kit::row(ui, tr("multiphysics.components_cell"), |ui| {
            crate::ui_kit::bounded_cell_editor(
                ui,
                ui.scope_id().with("components_cell"),
                &mut state.components.cell,
                tr("multiphysics.components_cell_auto"),
                false,
                16,
            )
        });
        crate::scenario::publish_target("multiphysics:components_cell", r.rect);
    }
    let r = crate::ui_kit::checkbox(ui, tr("multiphysics.components_replace"), &mut state.components.replace);
    crate::scenario::publish_target("multiphysics:components_replace", r.rect);
    let margin = if fit { number(&state.components.margin) } else { Some(0.0) };
    let cell = if state.components.cell.trim().is_empty() { Ok(None) } else { number(&state.components.cell).map(Some).ok_or(()) };
    let ready = !state.components.chosen.is_empty() && margin.is_some() && cell.is_ok();
    if button(ui, "multiphysics.components_add", "multiphysics:components_add", ready).clicked()
        && let (Some(margin), Ok(cell)) = (margin, cell)
        && let Some(case) = state.drafts.get(&study.id).and_then(|draft| serde_json::from_str::<Value>(&draft.case_text).ok())
    {
        let options = json!({
            "representation": match state.components.representation {
                Representation::Auto => "auto",
                Representation::Obstacle => "obstacle",
                Representation::Surface => "surface",
            },
            "origin_mm": if fit { None } else { origin },
            "fit_domain": fit,
            "margin_mm": margin,
            "cell_mm": if fit { cell } else { None },
            "replace_parts": state.components.replace,
        });
        let bodies: Vec<String> = state.components.chosen.iter().cloned().collect();
        match app.run("multiphysics.case.component.import", json!({"case": case, "bodies": bodies, "options": options})) {
            Ok(reply) => match serde_json::from_value::<solidcraft_engine::multiphysics::case_scene::Edit>(reply["edit"].clone()) {
                Ok(edit) => {
                    if let Ok(origin) = serde_json::from_value::<[f64; 3]>(reply["origin_mm"].clone()) {
                        state.cad_origin.insert(study.id.clone(), origin);
                    }
                    state.components.notes = reply["notes"].as_array().into_iter().flatten().filter_map(|n| n.as_str().map(str::to_string)).collect();
                    state.case_view.reset();
                    apply_case_edit(state, study, Ok(edit));
                }
                Err(error) => state.error = Some(error.to_string()),
            },
            Err(error) => state.error = Some(error),
        }
    }
    for note in &state.components.notes {
        ui.small(note.as_str());
    }
}

/// Text fields of the panel's selection tools (kept as typed until used).
#[derive(Default)]
struct CaseTool {
    conductivity: String,
    speed: String,
    pressure: String,
    /// Lid opening along the face's two in-plane axes: [u from, u to, v from, v to] in cells.
    opening: [String; 4],
    /// The face the opening fields were filled for.
    opening_face: Option<usize>,
    /// Device inputs, as typed (filled with typical values on first use).
    devices_ready: bool,
    pcb_thickness: String,
    /// `thickness_um:coverage` per copper layer, comma separated.
    pcb_layers: String,
    led_power: String,
    led_efficiency: String,
    led_rjc: String,
    battery_current: String,
    battery_resistance: String,
    battery_entropic: String,
    battery_cells: String,
    plate_name: String,
    plate_open: String,
    /// Local mesh: cell split factor (2..4) and margin around the picked parts, in cells.
    refine_factor: usize,
    refine_margin: String,
    /// Zoom-in: refinement (1..4), margin in parent cells and the child project's name.
    zoom_refine: usize,
    zoom_margin: String,
    zoom_name: String,
    /// Rotating region: name, speed in rpm, axis (0 x, 1 y, 2 z) and margin in cells.
    rotating_ready: bool,
    rotating_name: String,
    rotating_rpm: String,
    rotating_axis: usize,
    rotating_margin: String,
    /// Fluid volume check: the report for the case text it was made from (by hash), and a
    /// request to run it (Tools → Fluid Volume).
    fluid: Option<(u64, Result<solidcraft_engine::multiphysics::fluid_volume::FluidVolume, String>)>,
    fluid_requested: bool,
    /// Fans: an internal fan (picked part) instead of an external one (domain face).
    fan_internal: bool,
    /// Internal fan inputs: name, direction (index into `FAN_DIRECTIONS`) and the curve as
    /// `flow L/s : pressure Pa` pairs, comma separated.
    fan_ready: bool,
    fan_name: String,
    fan_direction: usize,
    fan_curve: String,
}

/// Results Mesh: the computational grid of a retained field file, drawn as cell lines on the
/// faces of its bounding box in the 3-D field view, with its size.
fn mesh_view(
    app: &mut SolveApp,
    state: &mut State,
    study: &Study,
    run: &solidcraft_engine::multiphysics::RunEvidence,
    artifact: &str,
    ui: &mut egui::Ui,
) {
    ui.label(RichText::new(tr("multiphysics.tree_results_mesh")).strong());
    #[cfg(not(target_arch = "wasm32"))]
    {
        let grid = state.inspected.as_ref().and_then(|cached| {
            if !cached.key.matches(app, &study.id, &run.id, artifact) {
                return None;
            }
            match &cached.outcome {
                Ok(output) => match &output.spatial {
                    Some(Ok(grid)) => {
                        // Node positions per axis: every cell face level.
                        let nodes: [Vec<f64>; 3] = [0, 1, 2].map(|a| {
                            let mut list: Vec<f64> = grid.cells.iter().flat_map(|c| [c[a][0], c[a][1]]).collect();
                            list.sort_by(f64::total_cmp);
                            let tolerance = 1e-12 * (grid.bounds[a][1] - grid.bounds[a][0]).abs().max(1e-300);
                            list.dedup_by(|x, y| (*x - *y).abs() <= tolerance);
                            list
                        });
                        Some((grid.cells.len(), grid.bounds, nodes))
                    }
                    _ => None,
                },
                Err(_) => None,
            }
        });
        let Some((cell_count, bounds, nodes)) = grid else {
            if !inspection_busy(state) {
                start_inspection(app, state, &study.id, &run.id, artifact);
            }
            ui.small(tr("multiphysics.results_mesh_loading"));
            return;
        };
        let widths = |a: usize| nodes[a].windows(2).map(|w| w[1] - w[0]).fold((f64::INFINITY, 0.0_f64), |(lo, hi), w| (lo.min(w), hi.max(w)));
        egui::Grid::new("multiphysics_results_mesh").show(ui, |ui| {
            ui.label(tr("multiphysics.results_mesh_cells"));
            ui.monospace(cell_count.to_string());
            ui.end_row();
            for (a, name) in ["x", "y", "z"].into_iter().enumerate() {
                let (lo, hi) = widths(a);
                ui.label(name);
                ui.monospace(format!("{} · {:.4e} – {:.4e} m", nodes[a].len().saturating_sub(1), lo, hi));
                ui.end_row();
            }
        });
        let mut lines = Vec::new();
        for face in 0..6usize {
            let axis = face / 2;
            let fixed = bounds[axis][face % 2];
            let [u, v] = solidcraft_engine::multiphysics::case_scene::face_axes(axis);
            for (along, across) in [(u, v), (v, u)] {
                let stride = nodes[along].len().div_ceil(64).max(1);
                for &at in nodes[along].iter().step_by(stride) {
                    let mut a = [0.0; 3];
                    a[axis] = fixed;
                    a[along] = at;
                    a[across] = bounds[across][0];
                    let mut b = a;
                    b[across] = bounds[across][1];
                    lines.push(vec![(a, 0.5), (b, 0.5)]);
                }
            }
        }
        if state.view.field_3d.overlay.as_ref().is_none_or(|o| o.lines.len() != lines.len()) {
            state.view.field_3d.overlay =
                Some(std::sync::Arc::new(crate::multiphysics_field_view::Overlay { lines, range: (0.0, 1.0), ..Default::default() }));
            state.view.show_3d = true;
        }
        ui.small(tr("multiphysics.flow_in_view"));
    }
    #[cfg(target_arch = "wasm32")]
    {
        let _ = (app, state, study, run, artifact);
    }
}

/// The newest flow field of a run: its last VTU frame that is not the thermal export.
fn flow_artifact(run: &solidcraft_engine::multiphysics::RunEvidence) -> Option<String> {
    run.artifacts
        .iter()
        .map(|artifact| artifact.name.as_str())
        .filter(|name| name.ends_with(".vtu") && !name.starts_with("thermal"))
        .max()
        .map(str::to_string)
}

/// Transferred Boundary Conditions: a zoom-in child project around the picked parts, fed by
/// this project's latest run.
fn zoom_tools(app: &mut SolveApp, state: &mut State, study: &Study, parts: &[String], ui: &mut egui::Ui) {
    use solidcraft_engine::multiphysics::case_scene as cs;
    let t = Tokens::get();
    ui.separator();
    ui.label(RichText::new(tr("multiphysics.zoom_title")).strong());
    ui.small(tr("multiphysics.zoom_hint"));
    let Some((run, artifact)) = study.runs.iter().rev().find_map(|run| flow_artifact(run).map(|name| (run, name))) else {
        ui.colored_label(t.text_dim, tr("multiphysics.zoom_needs_run"));
        return;
    };
    ui.small(solidcraft_ui_text::format("multiphysics.zoom_source", &[("run", &run.id), ("file", &artifact)]));
    let Ok(parent) = cs::scene(&run.input_case) else {
        ui.colored_label(t.text_dim, tr("multiphysics.zoom_needs_run"));
        return;
    };
    let tool = &mut state.case_tool;
    if tool.zoom_refine == 0 {
        tool.zoom_refine = 2;
        tool.zoom_margin = "2".into();
        tool.zoom_name = format!("{} zoom", study.name);
    }
    ui.horizontal_wrapped(|ui| {
        ui.label(tr("multiphysics.zoom_refine"));
        for factor in 1..=cs::MAX_REFINEMENT {
            let r = crate::ui_kit::choice_chip(ui, &format!("multiphysics_zoom_{factor}"), &format!("× {factor}"), tool.zoom_refine == factor);
            crate::scenario::publish_target(&format!("multiphysics:zoom_refine:{factor}"), r.rect);
            if r.clicked() {
                tool.zoom_refine = factor;
            }
        }
    });
    tool_field(ui, "multiphysics.zoom_margin", "multiphysics:zoom_margin", &mut tool.zoom_margin);
    tool_field(ui, "multiphysics.zoom_name", "multiphysics:zoom_name", &mut tool.zoom_name);
    // The box: the picked parts (as the run's grid places them) plus the margin.
    let margin = number(&tool.zoom_margin).filter(|m| (0.0..=64.0).contains(m));
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for name in parts {
        if let Some(part) = parent.part(name) {
            for a in 0..3 {
                low[a] = low[a].min(part.min[a]);
                high[a] = high[a].max(part.max[a]);
            }
        }
    }
    let boxed = low.iter().all(|v| v.is_finite());
    if let (true, Some(m)) = (boxed, margin) {
        for a in 0..3 {
            low[a] = (low[a] - m * parent.cell(a)).max(0.0);
            high[a] = (high[a] + m * parent.cell(a)).min(parent.lengths[a]);
        }
        ui.small(format!("x {:.4}–{:.4} · y {:.4}–{:.4} · z {:.4}–{:.4} m", low[0], high[0], low[1], high[1], low[2], high[2]));
    } else {
        ui.colored_label(t.text_dim, tr("multiphysics.zoom_pick"));
    }
    let (refine, name) = (tool.zoom_refine, tool.zoom_name.trim().to_string());
    let ready = boxed && margin.is_some() && !name.is_empty();
    if !button(ui, "multiphysics.zoom_create", "multiphysics:zoom_create", ready).clicked() {
        return;
    }
    let project = app.session.doc.multiphysics.clone();
    let zoom = match post::zoom_from_run(&project, &study.id, &run.id, &artifact, low, high, refine) {
        Ok(zoom) => zoom,
        Err(error) => {
            state.error = Some(error.to_string());
            return;
        }
    };
    let Some(reply) = command(app, state, "multiphysics.study.clone", json!({"id": study.id, "name": name})) else { return };
    let Some(child_id) = reply["study"]["id"].as_str().map(str::to_string) else { return };
    let Some(child) = app.session.doc.multiphysics.studies.iter().find(|s| s.id == child_id).cloned() else { return };
    let mut draft = Draft::from_study(&child);
    draft.case_text = pretty(&zoom.case);
    state.drafts.insert(child_id.clone(), draft);
    state.selected = Some(child_id);
    state.case_view.reset();
    state.error = zoom.issues.first().cloned();
    let faces: Vec<String> = zoom.transferred.iter().map(|(face, kind)| format!("{face} ({kind})")).collect();
    let mut notice = solidcraft_ui_text::format("multiphysics.zoom_created", &[("faces", &faces.join(", "))]);
    for note in &zoom.notes {
        notice.push('\n');
        notice.push_str(note);
    }
    state.notice = Some(notice);
}

/// Global Mesh panel: the grid size, and local refinement around the parts picked in the
/// graphics window.
fn local_mesh_tools(state: &mut State, study: &Study, case: &Value, parts: &[String], ui: &mut egui::Ui) {
    use solidcraft_engine::multiphysics::case_scene as cs;
    let Ok(scene) = cs::scene(case) else { return };
    let t = Tokens::get();
    ui.separator();
    ui.label(RichText::new(tr("multiphysics.local_mesh_title")).strong());
    let [nx, ny, nz] = scene.dimensions;
    let finest = (0..3).flat_map(|a| scene.nodes[a].windows(2).map(|w| w[1] - w[0]).collect::<Vec<_>>()).fold(f64::INFINITY, f64::min);
    ui.label(solidcraft_ui_text::format(
        if scene.graded() { "multiphysics.local_mesh_graded" } else { "multiphysics.local_mesh_uniform" },
        &[("cells", &format!("{nx} × {ny} × {nz} = {}", nx * ny * nz)), ("finest", &format!("{:.4}", finest * 1000.0))],
    ));
    let tool = &mut state.case_tool;
    if tool.refine_factor == 0 {
        tool.refine_factor = 2;
    }
    if tool.refine_margin.is_empty() {
        tool.refine_margin = "1".into();
    }
    ui.horizontal_wrapped(|ui| {
        ui.label(tr("multiphysics.local_mesh_factor"));
        for factor in 2..=solidcraft_engine::multiphysics::case_scene::MAX_REFINEMENT {
            let r = crate::ui_kit::choice_chip(ui, &format!("multiphysics_refine_{factor}"), &format!("× {factor}"), tool.refine_factor == factor);
            crate::scenario::publish_target(&format!("multiphysics:refine_factor:{factor}"), r.rect);
            if r.clicked() {
                tool.refine_factor = factor;
            }
        }
    });
    tool_field(ui, "multiphysics.local_mesh_margin", "multiphysics:refine_margin", &mut tool.refine_margin);
    let margin = number(&tool.refine_margin);
    let factor = tool.refine_factor;
    if parts.is_empty() {
        ui.colored_label(t.text_dim, tr("multiphysics.local_mesh_pick"));
    }
    let mut result = None;
    ui.horizontal_wrapped(|ui| {
        if button(ui, "multiphysics.local_mesh_refine", "multiphysics:local_mesh_refine", !parts.is_empty() && margin.is_some()).clicked()
            && let Some(margin) = margin
        {
            result = Some(cs::refine_around(case, parts, factor, margin));
        }
        let graded = case.pointer("/mesh/coordinates").is_some_and(|c| !c.is_null());
        if button(ui, "multiphysics.local_mesh_reset", "multiphysics:local_mesh_reset", graded).clicked() {
            result = Some(cs::uniform_mesh(case));
        }
    });
    if let Some(result) = result {
        apply_case_edit(state, study, result);
    }
}

/// Rotating regions (MRF): a region around the picked parts, turning about x, y or z through
/// its centre; the case's regions are listed with their speed and a remove button.
fn rotating_tools(state: &mut State, study: &Study, case: &Value, parts: &[String], ui: &mut egui::Ui) {
    use solidcraft_engine::multiphysics::case_scene as cs;
    let Ok(scene) = cs::scene(case) else { return };
    let t = Tokens::get();
    ui.separator();
    ui.label(RichText::new(tr("multiphysics.rotating_title")).strong());
    let tool = &mut state.case_tool;
    if !tool.rotating_ready {
        tool.rotating_ready = true;
        tool.rotating_name = format!("rotor_{}", scene.rotating.len() + 1);
        tool.rotating_rpm = "1500".into();
        tool.rotating_axis = 2;
        tool.rotating_margin = "1".into();
    }
    tool_field(ui, "multiphysics.rotating_name", "multiphysics:rotating_name", &mut tool.rotating_name);
    tool_field(ui, "multiphysics.rotating_rpm", "multiphysics:rotating_rpm", &mut tool.rotating_rpm);
    ui.horizontal_wrapped(|ui| {
        ui.label(tr("multiphysics.rotating_axis"));
        for (axis, label) in ["X", "Y", "Z"].into_iter().enumerate() {
            let r = crate::ui_kit::choice_chip(ui, &format!("multiphysics_rotating_axis_{axis}"), label, tool.rotating_axis == axis);
            crate::scenario::publish_target(&format!("multiphysics:rotating_axis:{axis}"), r.rect);
            if r.clicked() {
                tool.rotating_axis = axis;
            }
        }
    });
    tool_field(ui, "multiphysics.rotating_margin", "multiphysics:rotating_margin", &mut tool.rotating_margin);
    let rpm = number(&tool.rotating_rpm);
    let margin = number(&tool.rotating_margin).filter(|m| m.fract() == 0.0 && *m >= 0.0).map(|m| m as usize);
    if parts.is_empty() {
        ui.colored_label(t.text_dim, tr("multiphysics.rotating_pick"));
    }
    let mut result = None;
    let ready = !parts.is_empty() && rpm.is_some() && margin.is_some() && !tool.rotating_name.trim().is_empty();
    if button(ui, "multiphysics.rotating_add", "multiphysics:rotating_add", ready).clicked()
        && let (Some(rpm), Some(margin)) = (rpm, margin)
    {
        result = Some(cs::add_rotating_region(case, tool.rotating_name.trim(), parts, tool.rotating_axis, rpm, margin));
        tool.rotating_ready = false;
    }
    for region in &scene.rotating {
        ui.horizontal_wrapped(|ui| {
            let axis = ["X", "Y", "Z"].into_iter().zip(region.axis).find(|(_, c)| c.abs() > 0.5).map_or("?", |(name, _)| name);
            ui.label(solidcraft_ui_text::format(
                "multiphysics.rotating_row",
                &[("name", &region.name), ("rpm", &format!("{}", region.rpm)), ("axis", axis)],
            ));
            let handle = format!("multiphysics:rotating_remove:{}", region.name);
            if button(ui, "multiphysics.rotating_remove", &handle, true).clicked() {
                result = Some(cs::remove_rotating_region(case, &region.name));
            }
        });
    }
    if let Some(result) = result {
        apply_case_edit(state, study, result);
    }
}

/// Fluid volume: fluid and solid cells and volumes, the connected fluid regions and the
/// openings each reaches. Sealed pockets are drawn in the graphics window and can be filled
/// with solid blocks.
fn fluid_volume_tools(state: &mut State, study: &Study, case: &Value, ui: &mut egui::Ui) {
    use solidcraft_engine::multiphysics::fluid_volume as fv;
    use std::hash::{Hash, Hasher};
    let t = Tokens::get();
    ui.separator();
    ui.label(RichText::new(tr("multiphysics.fluid_title")).strong());
    let key = {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        case.to_string().hash(&mut hasher);
        hasher.finish()
    };
    let tool = &mut state.case_tool;
    // A report for an older case is redone, so the figures always match the draft.
    let stale = tool.fluid.as_ref().is_some_and(|(k, _)| *k != key);
    let pressed = button(ui, "multiphysics.fluid_check", "multiphysics:fluid_check", true).clicked();
    if pressed || stale || std::mem::take(&mut tool.fluid_requested) {
        tool.fluid = Some((key, fv::fluid_volume(case).map_err(|e| e.to_string())));
    }
    let litres = |m3: f64| format!("{:.4}", m3 * 1000.0);
    let mut fill = false;
    state.case_view.pockets.clear();
    match tool.fluid.as_ref().map(|(_, report)| report) {
        None => {
            ui.colored_label(t.text_dim, tr("multiphysics.fluid_hint"));
        }
        Some(Err(error)) => {
            ui.colored_label(t.error, error);
        }
        Some(Ok(report)) => {
            let percent = |part: f64| if report.total_volume_m3 > 0.0 { 100.0 * part / report.total_volume_m3 } else { 0.0 };
            ui.label(solidcraft_ui_text::format(
                "multiphysics.fluid_fluid",
                &[
                    ("cells", &report.fluid_cells.to_string()),
                    ("volume", &litres(report.fluid_volume_m3)),
                    ("percent", &format!("{:.1}", percent(report.fluid_volume_m3))),
                ],
            ));
            ui.label(solidcraft_ui_text::format(
                "multiphysics.fluid_solid",
                &[
                    ("cells", &report.solid_cells.to_string()),
                    ("volume", &litres(report.solid_volume_m3)),
                    ("percent", &format!("{:.1}", percent(report.solid_volume_m3))),
                ],
            ));
            if report.openings.is_empty() {
                ui.colored_label(t.warning, tr("multiphysics.fluid_no_openings"));
            } else {
                ui.label(solidcraft_ui_text::format("multiphysics.fluid_openings", &[("openings", &report.openings.join(", "))]));
            }
            for (n, region) in report.regions.iter().enumerate() {
                let number = (n + 1).to_string();
                let cells = region.cells.to_string();
                let volume = litres(region.volume_m3);
                if region.sealed() {
                    state.case_view.pockets.push((region.min, region.max));
                    let text =
                        solidcraft_ui_text::format("multiphysics.fluid_region_sealed", &[("n", &number), ("cells", &cells), ("volume", &volume)]);
                    let r = ui.colored_label(t.warning, text);
                    crate::scenario::publish_target(&format!("multiphysics:fluid_region:{}", n + 1), r.rect);
                } else {
                    let text = solidcraft_ui_text::format(
                        "multiphysics.fluid_region",
                        &[("n", &number), ("cells", &cells), ("volume", &volume), ("openings", &region.openings.join(", "))],
                    );
                    let r = ui.label(text);
                    crate::scenario::publish_target(&format!("multiphysics:fluid_region:{}", n + 1), r.rect);
                }
            }
            if report.unlisted_regions > 0 {
                ui.label(solidcraft_ui_text::format(
                    "multiphysics.fluid_more",
                    &[("regions", &report.unlisted_regions.to_string()), ("cells", &report.unlisted_cells.to_string())],
                ));
            }
            for part in report.solids.iter().filter(|p| p.overlapping_cells > 0) {
                ui.colored_label(
                    t.warning,
                    solidcraft_ui_text::format("multiphysics.fluid_overlap", &[("name", &part.name), ("cells", &part.overlapping_cells.to_string())]),
                );
            }
            for name in &report.empty_bodies {
                ui.colored_label(t.warning, solidcraft_ui_text::format("multiphysics.fluid_empty_body", &[("name", name)]));
            }
            let sealed = report.sealed().count();
            if sealed == 0 && !report.openings.is_empty() {
                ui.colored_label(t.text_dim, tr("multiphysics.fluid_no_pockets"));
            }
            fill = button(ui, "multiphysics.fluid_fill", "multiphysics:fluid_fill", sealed > 0 && !report.openings.is_empty()).clicked();
            crate::scenario::publish_count("multiphysics:fluid_sealed", sealed as f64);
        }
    }
    if fill {
        apply_case_edit(state, study, fv::fill_sealed_pockets(case));
    }
}

/// `flow L/s : pressure Pa` pairs, comma separated, as (m³/s, Pa) curve points.
fn parse_fan_curve(text: &str) -> Option<Vec<(f64, f64)>> {
    text.split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(|item| {
            let (flow, pressure) = item.split_once(':')?;
            Some((number(flow)? / 1000.0, number(pressure)?))
        })
        .collect()
}

/// Internal fans: the picked part becomes a fan block blowing along a grid direction with a
/// pressure-flow curve; the case's internal fans are listed with a remove button.
fn internal_fan_tools(state: &mut State, study: &Study, case: &Value, parts: &[String], ui: &mut egui::Ui) {
    use solidcraft_engine::multiphysics::{airflow::FAN_DIRECTIONS, case_scene as cs};
    let Ok(scene) = cs::scene(case) else { return };
    let t = Tokens::get();
    ui.separator();
    ui.label(RichText::new(tr("multiphysics.internal_fan_title")).strong());
    let tool = &mut state.case_tool;
    if !tool.fan_ready {
        tool.fan_ready = true;
        tool.fan_name = format!("fan_{}", scene.fans.len() + 1);
        tool.fan_direction = 0;
        // A typical 60 mm axial fan: 0 L/s at 40 Pa, 10 L/s at 20 Pa, 16 L/s free delivery.
        tool.fan_curve = "0:40, 10:20, 16:0".into();
    }
    tool_field(ui, "multiphysics.internal_fan_name", "multiphysics:internal_fan_name", &mut tool.fan_name);
    ui.horizontal_wrapped(|ui| {
        ui.label(tr("multiphysics.internal_fan_direction"));
        for (index, label) in FAN_DIRECTIONS.iter().enumerate() {
            let r = crate::ui_kit::choice_chip(ui, &format!("multiphysics_internal_fan_direction_{index}"), label, tool.fan_direction == index);
            crate::scenario::publish_target(&format!("multiphysics:internal_fan_direction:{label}"), r.rect);
            if r.clicked() {
                tool.fan_direction = index;
            }
        }
    });
    tool_field(ui, "multiphysics.internal_fan_curve", "multiphysics:internal_fan_curve", &mut tool.fan_curve);
    let curve = parse_fan_curve(&tool.fan_curve);
    if curve.is_none() {
        ui.colored_label(t.error, tr("multiphysics.internal_fan_curve_bad"));
    }
    if parts.is_empty() {
        ui.colored_label(t.text_dim, tr("multiphysics.internal_fan_pick"));
    }
    let mut result = None;
    let ready = !parts.is_empty() && curve.is_some() && !tool.fan_name.trim().is_empty();
    if button(ui, "multiphysics.internal_fan_add", "multiphysics:internal_fan_add", ready).clicked()
        && let Some(curve) = &curve
    {
        let direction = FAN_DIRECTIONS.get(tool.fan_direction).copied().unwrap_or("+x");
        result = Some(cs::add_internal_fan(case, tool.fan_name.trim(), parts, direction, curve));
        tool.fan_ready = false;
    }
    for fan in &scene.fans {
        ui.horizontal_wrapped(|ui| {
            let index = fan.direction.iter().position(|c| c.abs() > 0.5).unwrap_or(0);
            let direction = FAN_DIRECTIONS.get(2 * index + usize::from(fan.direction[index] < 0.0)).copied().unwrap_or("?");
            ui.label(solidcraft_ui_text::format("multiphysics.internal_fan_row", &[("name", &fan.name), ("direction", direction)]));
            let handle = format!("multiphysics:internal_fan_remove:{}", fan.name);
            if button(ui, "multiphysics.internal_fan_remove", &handle, true).clicked() {
                result = Some(cs::remove_internal_fan(case, &fan.name));
            }
        });
    }
    if let Some(result) = result {
        if result.is_ok() {
            state.case_view.selection.clear();
        }
        apply_case_edit(state, study, result);
    }
}

/// One labelled text field of the selection tools.
fn tool_field(ui: &mut egui::Ui, label: &str, handle: &str, value: &mut String) {
    let r = crate::ui_kit::row(ui, tr(label), |ui| crate::ui_kit::bounded_cell_editor(ui, ui.scope_id().with(handle), value, "", false, 128));
    crate::scenario::publish_target(handle, r.rect);
}

/// PCB, LED and battery inputs for the one picked part, compiled into its solid record.
fn device_tools(state: &mut State, study: &Study, case: &Value, node: shell::Node, parts: &[String], ui: &mut egui::Ui) {
    use solidcraft_engine::multiphysics::devices;
    let tool = &mut state.case_tool;
    if !tool.devices_ready {
        tool.devices_ready = true;
        tool.pcb_thickness = "1.6".into();
        tool.pcb_layers = "35:0.9, 35:0.3, 35:0.3, 35:0.9".into();
        tool.led_power = "1".into();
        tool.led_efficiency = "0.35".into();
        tool.battery_current = "10".into();
        tool.battery_resistance = "0.02".into();
        tool.battery_entropic = "0".into();
        tool.battery_cells = "1".into();
        tool.plate_name = "perforated_plate".into();
        tool.plate_open = "0.5".into();
    }
    ui.separator();
    let [part] = parts else {
        ui.small(tr("multiphysics.device_pick_one"));
        return;
    };
    let result = match node {
        shell::Node::Pcb => {
            tool_field(ui, "multiphysics.device_pcb_thickness", "multiphysics:device_pcb_thickness", &mut tool.pcb_thickness);
            tool_field(ui, "multiphysics.device_pcb_layers", "multiphysics:device_pcb_layers", &mut tool.pcb_layers);
            let layers: Option<Vec<devices::CopperLayer>> = tool
                .pcb_layers
                .split(',')
                .map(str::trim)
                .filter(|item| !item.is_empty())
                .map(|item| {
                    let (thickness, coverage) = item.split_once(':')?;
                    Some(devices::CopperLayer { thickness_um: number(thickness)?, coverage: number(coverage)? })
                })
                .collect();
            let stack = layers.zip(number(&tool.pcb_thickness)).map(|(layers, thickness)| devices::PcbStack {
                board_thickness_mm: thickness,
                layers,
                dielectric_w_per_m_k: devices::FR4_W_PER_M_K,
                normal_axis: None,
            });
            if !button(ui, "multiphysics.device_apply", "multiphysics:device_apply", stack.is_some()).clicked() {
                return;
            }
            let Some(stack) = stack else { return };
            devices::apply_pcb(case, part, &stack).map(|(edit, k)| {
                let note = solidcraft_ui_text::format(
                    "multiphysics.device_pcb_result",
                    &[("in_plane", &format!("{:.2}", k.in_plane_w_per_m_k)), ("through", &format!("{:.3}", k.through_plane_w_per_m_k))],
                );
                (edit, note)
            })
        }
        shell::Node::Leds => {
            tool_field(ui, "multiphysics.device_led_power", "multiphysics:device_led_power", &mut tool.led_power);
            tool_field(ui, "multiphysics.device_led_efficiency", "multiphysics:device_led_efficiency", &mut tool.led_efficiency);
            tool_field(ui, "multiphysics.device_led_rjc", "multiphysics:device_led_rjc", &mut tool.led_rjc);
            let rjc = if tool.led_rjc.trim().is_empty() { Some(None) } else { number(&tool.led_rjc).map(Some) };
            let led = number(&tool.led_power).zip(number(&tool.led_efficiency)).zip(rjc).map(|((power, efficiency), rjc)| devices::Led {
                electrical_power_w: power,
                optical_efficiency: efficiency,
                junction_to_case_k_per_w: rjc,
                max_junction_temperature_k: None,
            });
            if !button(ui, "multiphysics.device_apply", "multiphysics:device_apply", led.is_some()).clicked() {
                return;
            }
            let Some(led) = led else { return };
            devices::apply_led(case, part, &led)
                .map(|(edit, heat)| (edit, solidcraft_ui_text::format("multiphysics.device_heat_result", &[("watts", &format!("{heat:.3}"))])))
        }
        _ => {
            tool_field(ui, "multiphysics.device_battery_current", "multiphysics:device_battery_current", &mut tool.battery_current);
            tool_field(ui, "multiphysics.device_battery_resistance", "multiphysics:device_battery_resistance", &mut tool.battery_resistance);
            tool_field(ui, "multiphysics.device_battery_entropic", "multiphysics:device_battery_entropic", &mut tool.battery_entropic);
            tool_field(ui, "multiphysics.device_battery_cells", "multiphysics:device_battery_cells", &mut tool.battery_cells);
            let cells = tool.battery_cells.trim().parse::<u32>().ok();
            let battery = number(&tool.battery_current).zip(number(&tool.battery_resistance)).zip(number(&tool.battery_entropic)).zip(cells).map(
                |(((current, resistance), entropic), cells)| devices::Battery {
                    current_a: current,
                    internal_resistance_ohm: resistance,
                    entropic_coefficient_v_per_k: entropic,
                    temperature_k: None,
                    cells,
                },
            );
            if !button(ui, "multiphysics.device_apply", "multiphysics:device_apply", battery.is_some()).clicked() {
                return;
            }
            let Some(battery) = battery else { return };
            devices::apply_battery(case, part, &battery)
                .map(|(edit, heat)| (edit, solidcraft_ui_text::format("multiphysics.device_heat_result", &[("watts", &format!("{heat:.3}"))])))
        }
    };
    match result {
        Ok((edit, note)) => {
            if apply_case_edit(state, study, Ok(edit)) {
                state.notice = Some(note);
            }
        }
        Err(error) => state.error = Some(error.to_string()),
    }
}

/// A perforated plate on each picked domain face.
fn plate_tools(state: &mut State, study: &Study, case: &Value, faces: &[usize], ui: &mut egui::Ui) {
    use solidcraft_engine::multiphysics::devices;
    let tool = &mut state.case_tool;
    if tool.plate_name.is_empty() {
        tool.plate_name = "perforated_plate".into();
    }
    if tool.plate_open.is_empty() {
        tool.plate_open = "0.5".into();
    }
    tool_field(ui, "multiphysics.device_plate_name", "multiphysics:device_plate_name", &mut tool.plate_name);
    tool_field(ui, "multiphysics.device_plate_open", "multiphysics:device_plate_open", &mut tool.plate_open);
    let open = number(&tool.plate_open);
    if let Some(open) = open
        && let Ok(k) = devices::perforated_loss(open)
    {
        ui.small(solidcraft_ui_text::format("multiphysics.device_plate_loss", &[("k", &format!("{k:.2}"))]));
    }
    let name = tool.plate_name.trim().to_string();
    if !button(ui, "multiphysics.device_apply", "multiphysics:device_apply", !faces.is_empty() && open.is_some() && !name.is_empty()).clicked() {
        return;
    }
    let Some(open) = open else { return };
    let mut current = case.clone();
    let mut last = None;
    for (index, face) in faces.iter().enumerate() {
        let label = if faces.len() == 1 { name.clone() } else { format!("{name}_{}", index + 1) };
        match devices::apply_perforated_plate(&current, &label, *face, open, None) {
            Ok((edit, _)) => {
                current = edit.case.clone();
                last = Some(edit);
            }
            Err(error) => {
                state.error = Some(error.to_string());
                return;
            }
        }
    }
    if let Some(edit) = last {
        apply_case_edit(state, study, Ok(edit));
    }
}

/// What clicking in the 3-D case view picks for a tree node. Fans pick a domain face for an
/// external fan, or a part for an internal one.
fn case_mode(node: shell::Node, internal_fan: bool) -> crate::multiphysics_case_view::Mode {
    use crate::multiphysics_case_view::Mode;
    match node {
        shell::Node::Domain => Mode::Domain,
        shell::Node::Fans if internal_fan => Mode::Parts,
        shell::Node::Boundary | shell::Node::Fans | shell::Node::Perforated => Mode::Faces,
        _ => Mode::Parts,
    }
}

/// Put a case edit made by pointing into the draft (Apply still commits it).
fn apply_case_edit(
    state: &mut State,
    study: &Study,
    result: solidcraft_engine::multiphysics::Result<solidcraft_engine::multiphysics::case_scene::Edit>,
) -> bool {
    match result {
        Ok(edit) => {
            if let Some(draft) = state.drafts.get_mut(&study.id) {
                draft.case_text = pretty(&edit.case);
            }
            state.entry = edit.pointer;
            state.error = edit.issues.first().cloned();
            state.notice = Some(tr("multiphysics.case_edit_drafted").to_string());
            true
        }
        Err(error) => {
            state.error = Some(error.to_string());
            false
        }
    }
}

fn number(text: &str) -> Option<f64> {
    text.trim().parse::<f64>().ok().filter(|v| v.is_finite())
}

/// The panel's half of picking in the graphics window: what is picked, and what to do with it.
fn case_tools(app: &mut SolveApp, state: &mut State, study: &Study, node: shell::Node, ui: &mut egui::Ui) {
    use crate::multiphysics_case_view::Mode;
    use solidcraft_engine::multiphysics::case_scene as cs;
    let Some(case) = state.drafts.get(&study.id).and_then(|draft| serde_json::from_str::<Value>(&draft.case_text).ok()) else { return };
    let Ok(scene) = cs::scene(&case) else { return };
    state.case_view.own(&study.id);
    if node == shell::Node::FromComponents {
        components_panel(app, state, study, ui);
        ui.separator();
        return;
    }
    let origin = state.cad_origin.get(&study.id).copied();
    if node == shell::Node::Fans {
        ui.horizontal_wrapped(|ui| {
            for (internal, key) in [(false, "multiphysics.fan_kind_external"), (true, "multiphysics.fan_kind_internal")] {
                let r =
                    crate::ui_kit::choice_chip(ui, &format!("multiphysics_fan_kind_{internal}"), tr(key), state.case_tool.fan_internal == internal);
                crate::scenario::publish_target(&format!("multiphysics:fan_kind:{}", if internal { "internal" } else { "external" }), r.rect);
                if r.clicked() && state.case_tool.fan_internal != internal {
                    state.case_tool.fan_internal = internal;
                    state.case_view.selection.clear();
                }
            }
        });
    }
    let mode = case_mode(node, state.case_tool.fan_internal);
    let t = Tokens::get();
    ui.label(RichText::new(tr("multiphysics.case_from_graphics")).strong());
    match mode {
        Mode::Domain => {
            ui.small(tr("multiphysics.case_hint_domain"));
            let size = format!(
                "{} × {} × {} · {:.4} × {:.4} × {:.4} m",
                scene.dimensions[0], scene.dimensions[1], scene.dimensions[2], scene.lengths[0], scene.lengths[1], scene.lengths[2]
            );
            ui.label(size);
            for (face, name) in cs::PATCHES.iter().enumerate() {
                let mut step = 0;
                ui.horizontal(|ui| {
                    ui.label(*name);
                    if button(ui, "multiphysics.case_face_in", &format!("multiphysics:domain_in:{name}"), true).clicked() {
                        step = -1;
                    }
                    if button(ui, "multiphysics.case_face_out", &format!("multiphysics:domain_out:{name}"), true).clicked() {
                        step = 1;
                    }
                });
                if step != 0 && apply_case_edit(state, study, cs::move_domain_face(&case, face, step)) {
                    state.case_view.moved(face, step, scene.edge_cell(face));
                }
            }
            fluid_volume_tools(state, study, &case, ui);
        }
        Mode::Parts => {
            // Bodies picked on the CAD model stand for the case parts brought in from them.
            let from_cad: Vec<String> = cad_bodies(app).iter().map(|name| cs::part_name(name)).filter(|name| scene.part(name).is_some()).collect();
            if !from_cad.is_empty() && button(ui, "multiphysics.case_use_cad_bodies", "multiphysics:case_use_cad_bodies", true).clicked() {
                state.case_view.selection = from_cad.into_iter().map(crate::multiphysics_case_view::Pick::Part).collect();
            }
            let parts = state.case_view.parts();
            if parts.is_empty() {
                ui.small(tr("multiphysics.case_hint_parts"));
            } else {
                ui.label(RichText::new(parts.join(", ")).color(t.text));
            }
            // A device goes on one picked solid: the picked part fills the name field.
            if matches!(node, shell::Node::Contact | shell::Node::ThermalJoints | shell::Node::TwoResistor | shell::Node::HeatPipes)
                && let [only] = parts.as_slice()
                && state.feature_name != *only
            {
                state.feature_name = only.clone();
            }
            if node == shell::Node::SolidMaterials {
                if state.case_tool.conductivity.is_empty() {
                    state.case_tool.conductivity = "1".into();
                }
                let r = crate::ui_kit::row(ui, tr("multiphysics.case_conductivity"), |ui| {
                    crate::ui_kit::bounded_cell_editor(ui, ui.scope_id().with("case_conductivity"), &mut state.case_tool.conductivity, "", false, 24)
                });
                crate::scenario::publish_target("multiphysics:case_conductivity", r.rect);
                let k = number(&state.case_tool.conductivity);
                if button(ui, "multiphysics.case_apply_material", "multiphysics:case_apply_material", !parts.is_empty() && k.is_some()).clicked()
                    && let Some(k) = k
                {
                    apply_case_edit(state, study, cs::assign_solid_material(&case, &parts, k));
                }
            }
            if matches!(node, shell::Node::Pcb | shell::Node::Leds | shell::Node::Batteries) {
                device_tools(state, study, &case, node, &parts, ui);
            }
            if node == shell::Node::Mesh {
                local_mesh_tools(state, study, &case, &parts, ui);
            }
            if node == shell::Node::Transferred {
                zoom_tools(app, state, study, &parts, ui);
            }
            if node == shell::Node::Rotating {
                rotating_tools(state, study, &case, &parts, ui);
            }
            if node == shell::Node::Fans {
                internal_fan_tools(state, study, &case, &parts, ui);
            }
        }
        Mode::Faces => {
            // Faces picked on the CAD model that lie on a domain face pick that domain face.
            if let Some(origin) = origin {
                let mut from_cad: Vec<usize> = cad_face_points(app).into_iter().filter_map(|p| cs::patch_at(&case, origin, p)).collect();
                from_cad.sort_unstable();
                from_cad.dedup();
                if !from_cad.is_empty() && button(ui, "multiphysics.case_use_cad_faces", "multiphysics:case_use_cad_faces", true).clicked() {
                    state.case_view.selection = from_cad.into_iter().map(crate::multiphysics_case_view::Pick::Patch).collect();
                }
            }
            let faces = state.case_view.patches();
            if faces.is_empty() {
                ui.small(tr("multiphysics.case_hint_faces"));
            } else {
                let names: Vec<String> = faces
                    .iter()
                    .map(|face| format!("{} ({})", cs::PATCHES.get(*face).copied().unwrap_or("?"), tr(role_key(scene.patches.get(*face).copied()))))
                    .collect();
                ui.label(RichText::new(names.join(", ")).color(t.text));
            }
            if node == shell::Node::Perforated {
                ui.separator();
                plate_tools(state, study, &case, &faces, ui);
                ui.separator();
                return;
            }
            if node != shell::Node::Boundary {
                return;
            }
            if state.case_tool.speed.is_empty() {
                state.case_tool.speed = "0.1".into();
            }
            if state.case_tool.pressure.is_empty() {
                state.case_tool.pressure = "0".into();
            }
            let r = crate::ui_kit::row(ui, tr("multiphysics.case_inlet_speed"), |ui| {
                crate::ui_kit::bounded_cell_editor(ui, ui.scope_id().with("case_speed"), &mut state.case_tool.speed, "", false, 24)
            });
            crate::scenario::publish_target("multiphysics:case_speed", r.rect);
            let r = crate::ui_kit::row(ui, tr("multiphysics.case_outlet_pressure"), |ui| {
                crate::ui_kit::bounded_cell_editor(ui, ui.scope_id().with("case_pressure"), &mut state.case_tool.pressure, "", false, 24)
            });
            crate::scenario::publish_target("multiphysics:case_pressure", r.rect);
            let mut role = None;
            ui.horizontal_wrapped(|ui| {
                let any = !faces.is_empty();
                let speed = number(&state.case_tool.speed);
                let pressure = number(&state.case_tool.pressure);
                if button(ui, "multiphysics.case_role_inlet", "multiphysics:case_set_inlet", any && speed.is_some()).clicked()
                    && let Some(speed) = speed
                {
                    role = Some(cs::RoleInput::Inlet { speed });
                }
                if button(ui, "multiphysics.case_role_outlet", "multiphysics:case_set_outlet", any && pressure.is_some()).clicked()
                    && let Some(pressure) = pressure
                {
                    role = Some(cs::RoleInput::Outlet { pressure });
                }
                if button(ui, "multiphysics.case_role_wall", "multiphysics:case_set_wall", any).clicked() {
                    role = Some(cs::RoleInput::Wall);
                }
                if button(ui, "multiphysics.case_role_symmetry", "multiphysics:case_set_symmetry", any).clicked() {
                    role = Some(cs::RoleInput::Symmetry);
                }
            });
            if let Some(role) = role {
                let mut result = Ok(cs::Edit { case: case.clone(), issues: Vec::new(), pointer: None });
                for face in &faces {
                    result = result.and_then(|edit| cs::set_patch_role(&edit.case, *face, role));
                }
                apply_case_edit(state, study, result);
            }
            // Lids: one picked face, an opening in cells along its two in-plane axes.
            ui.separator();
            ui.label(RichText::new(tr("multiphysics.case_lids")).strong());
            let [face] = faces.as_slice() else {
                ui.small(tr("multiphysics.case_lids_hint"));
                return;
            };
            let [u, v] = cs::face_axes(face / 2);
            let (nu, nv) = (scene.dimensions[u], scene.dimensions[v]);
            if state.case_tool.opening_face != Some(*face) {
                state.case_tool.opening_face = Some(*face);
                state.case_tool.opening = [nu / 4, nu - nu / 4, nv / 4, nv - nv / 4].map(|n| n.to_string());
            }
            let axes = ["x", "y", "z"];
            for (slot, (axis, label)) in [(u, "from"), (u, "to"), (v, "from"), (v, "to")].into_iter().enumerate() {
                let text = solidcraft_ui_text::format(
                    if label == "from" { "multiphysics.case_opening_from" } else { "multiphysics.case_opening_to" },
                    &[("axis", axes.get(axis).copied().unwrap_or("?"))],
                );
                let r = crate::ui_kit::row(ui, &text, |ui| {
                    crate::ui_kit::bounded_cell_editor(
                        ui,
                        ui.scope_id().with(("case_opening", slot)),
                        &mut state.case_tool.opening[slot],
                        "",
                        false,
                        8,
                    )
                });
                crate::scenario::publish_target(&format!("multiphysics:case_opening:{slot}"), r.rect);
            }
            let cells: Option<Vec<usize>> = state.case_tool.opening.iter().map(|text| text.trim().parse::<usize>().ok()).collect();
            let opening = cells.and_then(|c| match c.as_slice() {
                [a, b, c, d] => Some(([*a, *c], [*b, *d])),
                _ => None,
            });
            if button(ui, "multiphysics.case_create_lids", "multiphysics:case_create_lids", opening.is_some()).clicked()
                && let Some((low, high)) = opening
            {
                apply_case_edit(state, study, cs::create_lids(&case, *face, low, high));
            }
        }
    }
    ui.separator();
}

fn role_key(role: Option<solidcraft_engine::multiphysics::case_scene::PatchRole>) -> &'static str {
    use solidcraft_engine::multiphysics::case_scene::PatchRole;
    match role {
        Some(PatchRole::Inlet) => "multiphysics.case_role_inlet",
        Some(PatchRole::Outlet) => "multiphysics.case_role_outlet",
        Some(PatchRole::Wall) => "multiphysics.case_role_wall",
        Some(PatchRole::Symmetry) => "multiphysics.case_role_symmetry",
        Some(PatchRole::Fan) => "multiphysics.case_role_fan",
        _ => "multiphysics.case_role_default",
    }
}

/// Goal kinds offered in the Goals panel, with their labels.
const GOAL_KINDS: [(&str, &str); 5] = [
    ("goal", "multiphysics.goal_kind_surface"),
    ("global_goal", "multiphysics.goal_kind_global"),
    ("point_goal", "multiphysics.goal_kind_point"),
    ("volume_goal", "multiphysics.goal_kind_volume"),
    ("equation_goal", "multiphysics.goal_kind_equation"),
];

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

/// The workspace's state, out of the app's slot (a fresh one when there is none yet).
pub fn take_state(app: &mut SolveApp) -> State {
    app.multiphysics_workspace.0.take().and_then(|state| state.downcast::<State>().ok()).map(|state| *state).unwrap_or_default()
}

/// Puts the workspace's state back into the app's slot.
pub fn put_state(app: &mut SolveApp, state: State) {
    app.multiphysics_workspace.0 = Some(Box::new(state));
}

fn peek(app: &SolveApp) -> Option<&State> {
    app.multiphysics_workspace.0.as_ref().and_then(|state| state.downcast_ref::<State>())
}

/// The hooks the app calls the workspace through.
pub fn hooks() -> solidcraft_ui_egui::multiphysics_workspace::Hooks {
    solidcraft_ui_egui::multiphysics_workspace::Hooks { show, tick, reset_document, close, take_view, put_view }
}

fn take_view(app: &mut SolveApp) -> solidcraft_ui_egui::multiphysics_workspace::Slot {
    let mut state = take_state(app);
    let view = std::mem::take(&mut state.view);
    put_state(app, state);
    solidcraft_ui_egui::multiphysics_workspace::Slot(Some(Box::new(view)))
}

fn put_view(app: &mut SolveApp, view: solidcraft_ui_egui::multiphysics_workspace::Slot) {
    let mut state = take_state(app);
    state.view = view.0.and_then(|view| view.downcast::<ViewState>().ok()).map(|view| *view).unwrap_or_default();
    put_state(app, state);
}

pub fn reset_document(app: &mut SolveApp) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        // Keep the one export slot occupied until retirement, even across New/Open.
        // Its captured UID prevents a completion from becoming the new document's notice.
        let mut old = take_state(app);
        let (export, inspection) = (old.export.take(), old.inspection.take());
        if let Some(job) = &inspection {
            job.cancelled.store(true, std::sync::atomic::Ordering::Release);
        }
        put_state(app, State { export, inspection, ..Default::default() });
    }
    #[cfg(target_arch = "wasm32")]
    {
        put_state(app, State::default());
    }
}

/// Finish the frozen export before the native host tears down its resources.
pub fn close(app: &mut SolveApp) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let mut state = take_state(app);
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
        put_state(app, state);
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
        let mut state = take_state(app);
        poll_export(app, &mut state);
        poll_inspection(app, &mut state);
        put_state(app, state);
    }
    if peek(app).is_some_and(|state| state.export_busy() || inspection_busy(state))
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
    let mut state = take_state(app);
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
    // Feature dialogs dock at the right of the graphics window, as a flow-analysis property
    // panel does, rather than floating over the ribbon or the model.
    if let Some(study) = selected
        && shell::feature_docked(&state)
    {
        Container::for_layout(LayoutMode::Docked, Surface::Panel)
            .edge("multiphysics_feature_panel", Edge::Right)
            .default_size(380.0)
            .size_range(300.0..=560.0)
            .resizable(true)
            .show(ui, |ui| shell::feature_panel(app, &mut state, study, report.as_deref(), ui));
    }
    let mut graphics = egui::Rect::NOTHING;
    Container::for_layout(LayoutMode::Docked, Surface::Viewport).central().show(ui, |ui| {
        graphics = ui.max_rect();
        shell::graphics(app, &mut state, selected, ui);
    });
    let ctx = ui.ctx().clone();
    shell::windows(app, &mut state, selected, report.as_deref(), graphics, &ctx);
    put_state(app, state);
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
                egui::Id::unique(("multiphysics_question", &study.id)),
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
        && command(app, state, "multiphysics.run.start", json!({"id":study.id})).is_some()
    {
        state.stage = Stage::Run;
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
fn feature_creation(app: &mut SolveApp, state: &mut State, study: &Study, node: shell::Node, ui: &mut egui::Ui) {
    let kind = match node {
        shell::Node::Fans => "fan",
        shell::Node::Porous => "porous",
        shell::Node::SolidMaterials => "solid_material",
        shell::Node::Goals => GOAL_KINDS.iter().find(|(kind, _)| *kind == state.goal_kind).map_or("goal", |(kind, _)| *kind),
        shell::Node::Contact | shell::Node::ThermalJoints => "contact_resistance",
        shell::Node::TwoResistor => "two_resistor",
        shell::Node::HeatPipes => "heat_pipe",
        shell::Node::Tec => "thermoelectric_cooler",
        _ => return,
    };
    // Devices placed on an existing thermal solid take that solid's name.
    let on_solid = matches!(kind, "solid_material" | "contact_resistance" | "two_resistor" | "heat_pipe");
    let owner = (app.session.uid, study.id.clone(), kind);
    if state.feature_owner.as_ref() != Some(&owner) {
        state.feature_owner = Some(owner);
        state.feature_name = if kind == "fan" { "external_fan".into() } else { String::new() };
        state.replace_fan_patch_conditions = false;
    }
    let compatible = state
        .drafts
        .get(&study.id)
        .and_then(|d| serde_json::from_str::<Value>(&d.case_text).ok())
        .is_some_and(|c| c["solver"].as_str() == Some("collocated_flow3d"));
    ui.small(tr("multiphysics.feature_draft_scope"));
    if node == shell::Node::Goals {
        let current = GOAL_KINDS.iter().find(|(id, _)| *id == kind).map_or("", |(_, label)| label);
        let r = crate::ui_kit::row(ui, tr("multiphysics.goal_kind"), |ui| {
            crate::ui_kit::dropdown_menu(ui, "multiphysics_goal_kind", tr(current), |ui| {
                for (id, label) in GOAL_KINDS {
                    let response = crate::ui_kit::list_picker(ui, &format!("multiphysics_goal_kind_{id}"), tr(label), id == kind);
                    crate::scenario::publish_target(&format!("multiphysics:goal_kind:{id}"), response.rect);
                    if response.clicked() {
                        state.goal_kind = id;
                        ui.close();
                    }
                }
            })
        });
        crate::scenario::publish_target("multiphysics:goal_kind", r.rect);
    }
    if kind != "fan" {
        let r = crate::ui_kit::row(ui, tr(if on_solid { "multiphysics.feature_body" } else { "multiphysics.feature_name" }), |ui| {
            crate::ui_kit::bounded_cell_editor(ui, ui.scope_id().with("new_feature_name"), &mut state.feature_name, "", false, 60)
        });
        crate::scenario::publish_target("multiphysics:feature_name", r.rect);
    } else {
        ui.small(tr("multiphysics.feature_fan_singleton"));
    }
    if kind == "fan" {
        let r = crate::ui_kit::checkbox(ui, tr("multiphysics.feature_fan_replace"), &mut state.replace_fan_patch_conditions);
        crate::scenario::publish_target("multiphysics:feature_fan_replace", r.rect);
    }
    let r = button(ui, "multiphysics.feature_add", "multiphysics:feature_add", compatible && !state.feature_name.trim().is_empty());
    if !compatible {
        r.clone().backing_on_hover_text(tr("multiphysics.feature_solver_required"));
    }
    if r.clicked() {
        let result = state
            .drafts
            .get(&study.id)
            .and_then(|draft| serde_json::from_str::<Value>(&draft.case_text).ok())
            .ok_or_else(|| tr("multiphysics.feature_solver_required").to_string())
            .and_then(|case| {
                app.run(
                    "multiphysics.case.feature.prepare",
                    json!({
                        "case": case,
                        "kind": kind,
                        "name": state.feature_name.trim(),
                        "replace_fan_patch_conditions": kind == "fan" && state.replace_fan_patch_conditions,
                    }),
                )
            });
        match result {
            Ok(prepared) => {
                if let Some(case) = prepared.get("case")
                    && let Some(draft) = state.drafts.get_mut(&study.id)
                {
                    draft.case_text = pretty(case);
                    state.entry = prepared["pointer"].as_str().map(str::to_string);
                    state.notice = Some(tr("multiphysics.feature_draft_scope").to_string());
                    state.error = prepared["issues"].as_array().and_then(|v| v.first()).and_then(Value::as_str).map(str::to_string);
                }
            }
            Err(e) => state.error = Some(e),
        }
    }
    ui.separator();
}
fn settings(app: &mut SolveApp, state: &mut State, study: &Study, node: shell::Node, entry: Option<&str>, ui: &mut egui::Ui) {
    if shell::Node::INPUT.contains(&node) || node == shell::Node::Mesh {
        case_tools(app, state, study, node, ui);
    }
    feature_creation(app, state, study, node, ui);
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
    egui::ScrollArea::vertical().id_salt(("multiphysics_settings", &study.id)).max_height(520.0).show(ui, |ui| {
        if raw {
            let response = crate::ui_kit::text_area(
                ui,
                tr("multiphysics.case_settings"),
                egui::Id::unique(("multiphysics_case_json", &study.id)),
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
                    Some(slot) => {
                        crate::multiphysics_property_grid::header(ui, &[tr("multiphysics.wizard_col_parameter"), tr("multiphysics.wizard_col_value")], &[0.46, 0.54]);
                        case_row(ui, &entry_title(pointer, slot), slot, 0, pointer, &mut 0)
                    }
                    None => {
                        ui.label(RichText::new(tr("multiphysics.tree_revert_missing")).small().color(Tokens::get().text_dim));
                        false
                    }
                },
                None => settings_table(ui, &mut value, node),
            };
            if changed {
                draft.case_text = pretty(&value);
            }
        }
    });
}

/// The case entries of one tree folder (top-level entries and the children of mixed blocks
/// such as `config` filed under it), or the whole case for Input data, as a Parameter | Value
/// table.
fn settings_table(ui: &mut egui::Ui, value: &mut Value, node: shell::Node) -> bool {
    let pointers: Vec<String> = if node == shell::Node::Input {
        value.as_object().map(|object| object.keys().map(|key| field_pointer("", key)).collect()).unwrap_or_default()
    } else {
        shell::entries(value, node)
    };
    if pointers.is_empty() {
        ui.label(RichText::new(tr("multiphysics.tree_folder_empty")).small().color(Tokens::get().text_dim));
        return false;
    }
    crate::multiphysics_property_grid::header(ui, &[tr("multiphysics.wizard_col_parameter"), tr("multiphysics.wizard_col_value")], &[0.46, 0.54]);
    let mut rows = 0;
    let mut changed = false;
    for pointer in pointers.iter().take(80) {
        if let Some(slot) = value.pointer_mut(pointer) {
            changed |= case_row(ui, &entry_title(pointer, slot), slot, 0, pointer, &mut rows);
        }
    }
    changed
}

/// A case entry's title: a list item's own `name`, `id`, `patch` or `type`, as the tree names it,
/// otherwise the last segment of its pointer.
fn entry_title(pointer: &str, value: &Value) -> String {
    let own = ["name", "id", "patch", "type"].iter().find_map(|field| value.get(*field).and_then(Value::as_str));
    match own {
        Some(own) if pointer.rsplit('/').next().is_some_and(|last| last.parse::<usize>().is_ok()) => own.to_string(),
        _ => shell::entry_label(pointer),
    }
}

/// The JSON pointer of `key` under `parent`, escaped as RFC 6901 asks.
fn field_pointer(parent: &str, key: &str) -> String {
    format!("{parent}/{}", key.replace('~', "~0").replace('/', "~1"))
}
/// The control for one case value, published as `multiphysics:case_field:<pointer>`.
fn field_editor(ui: &mut egui::Ui, name: &str, value: &mut Value, pointer: &str) -> bool {
    let mut changed = false;
    match value {
        Value::Bool(v) => {
            let r = crate::ui_kit::checkbox(ui, "", v);
            crate::scenario::publish_target(&format!("multiphysics:case_field:{pointer}"), r.rect);
            changed |= r.changed();
        }
        Value::Number(n) => {
            if let Some(mut v) = n.as_u64() {
                let r = crate::ui_kit::number_field(ui, &mut v, 0..=u64::MAX, 1.0, "");
                crate::scenario::publish_target(&format!("multiphysics:case_field:{pointer}"), r.rect);
                if r.changed() {
                    *n = v.into();
                    changed = true;
                }
            } else if let Some(mut v) = n.as_i64() {
                let r = crate::ui_kit::number_field(ui, &mut v, i64::MIN..=i64::MAX, 1.0, "");
                crate::scenario::publish_target(&format!("multiphysics:case_field:{pointer}"), r.rect);
                if r.changed() {
                    *n = v.into();
                    changed = true;
                }
            } else if let Some(mut v) = n.as_f64() {
                let r = crate::ui_kit::number_field(ui, &mut v, f64::MIN..=f64::MAX, 0.001, "");
                crate::scenario::publish_target(&format!("multiphysics:case_field:{pointer}"), r.rect);
                if r.changed()
                    && v.is_finite()
                    && let Some(updated) = serde_json::Number::from_f64(v)
                {
                    *n = updated;
                    changed = true;
                }
            }
        }
        Value::String(v) if name == "solver" || name == "schema" => {
            ui.monospace(v.as_str());
        }
        Value::String(v) => {
            let r = crate::ui_kit::text_field_reserving(ui, v, "", 0.0).1;
            crate::scenario::publish_target(&format!("multiphysics:case_field:{pointer}"), r.rect);
            changed |= r.changed();
        }
        _ => {
            ui.label(tr("multiphysics.not_specified"));
        }
    }
    changed
}

/// Most rows one case table draws, so a hostile case cannot stall the frame.
const CASE_TABLE_ROWS: usize = 400;

/// The entries of `node` in `case` as Parameter | Value rows under a band for the group, as
/// the wizard pages and the summary show them. Returns whether a value changed.
fn case_table(ui: &mut egui::Ui, case: &mut Value, node: shell::Node, rows: &mut usize) -> bool {
    let pointers = shell::entries(case, node);
    if pointers.is_empty() {
        return false;
    }
    crate::multiphysics_property_grid::band(ui, tr(node.key()));
    let mut changed = false;
    for pointer in pointers.iter().take(80) {
        if let Some(slot) = case.pointer_mut(pointer) {
            changed |= case_row(ui, &entry_title(pointer, slot), slot, 0, pointer, rows);
        }
    }
    changed
}

fn case_row(ui: &mut egui::Ui, name: &str, value: &mut Value, depth: usize, pointer: &str, rows: &mut usize) -> bool {
    if depth > 5 || *rows >= CASE_TABLE_ROWS {
        return false;
    }
    *rows += 1;
    let label = name.replace('_', " ");
    let mut changed = false;
    ui.push_id(pointer, |ui| match value {
        Value::Object(object) => {
            crate::multiphysics_property_grid::heading(ui, depth, &label);
            for (key, child) in object.iter_mut().take(80) {
                changed |= case_row(ui, key, child, depth + 1, &field_pointer(pointer, key), rows);
            }
        }
        Value::Array(values) if values.len() <= 4 && values.iter().all(|v| !v.is_object() && !v.is_array()) => {
            // A short tuple (a point, a direction, a range) is one row of fields.
            crate::multiphysics_property_grid::row(ui, depth, &label, rows.is_multiple_of(2), |ui| {
                let width = (ui.available_width() / values.len().max(1) as f32 - Tokens::get().space_4).max(32.0);
                for (index, item) in values.iter_mut().enumerate() {
                    ui.scope(|ui| {
                        ui.set_max_width(width);
                        changed |= field_editor(ui, &index.to_string(), item, &field_pointer(pointer, &index.to_string()));
                    });
                }
            });
        }
        Value::Array(values) if values.len() <= 16 => {
            crate::multiphysics_property_grid::heading(ui, depth, &label);
            for (index, item) in values.iter_mut().enumerate() {
                let own = ["name", "id", "patch", "type"].iter().find_map(|field| item.get(*field).and_then(Value::as_str).map(str::to_string));
                let title = own.unwrap_or_else(|| format!("{label} {}", index + 1));
                changed |= case_row(ui, &title, item, depth + 1, &field_pointer(pointer, &index.to_string()), rows);
            }
        }
        Value::Array(values) => {
            let count = values.len().to_string();
            crate::multiphysics_property_grid::row(ui, depth, &label, rows.is_multiple_of(2), |ui| {
                ui.label(solidcraft_ui_text::format("multiphysics.array_entries", &[("count", &count)]));
            });
        }
        _ => {
            crate::multiphysics_property_grid::row(ui, depth, &label, rows.is_multiple_of(2), |ui| changed |= field_editor(ui, name, value, pointer));
        }
    });
    changed
}
fn jobs(app: &mut SolveApp, state: &mut State, study: &Study, ui: &mut egui::Ui) {
    ui.heading(tr("multiphysics.jobs"));
    ui.add(egui::Label::new(tr("multiphysics.run_identity")).wrap());
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
                state.view.shown_run = Some(job.id.clone());
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
        Node::CutPlots
        | Node::Animation
        | Node::PointParams
        | Node::VolumeParams
        | Node::Compare
        | Node::Trajectories
        | Node::Particles
        | Node::Isosurfaces
        | Node::SurfacePlots
        | Node::SurfaceParams
        | Node::ResultsMesh => name.ends_with(".vtu"),
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
    // A run shown for the first time opens on its final flow field, not on the field at time 0.
    if state.view.defaulted_run.as_deref() != Some(run.id.as_str()) {
        state.view.defaulted_run = Some(run.id.clone());
        if let Some(index) = run
            .artifacts
            .iter()
            .enumerate()
            .filter(|(_, artifact)| artifact.name.starts_with("flow_") && artifact.name.ends_with(".vtu") && wanted(&artifact.name))
            .max_by(|a, b| a.1.name.cmp(&b.1.name))
            .map(|(index, _)| index)
        {
            state.selected_artifact = index;
        }
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
    let tools_first =
        study.model.is_some() && artifact.name.ends_with(".vtu") && matches!(focus, Node::Trajectories | Node::Particles | Node::Isosurfaces);
    if artifact.name.ends_with(".vtu") || artifact.name.ends_with(".pvd") {
        // A study of the model draws a flow tool's lines in a tall view of the model, so the
        // tool's seeds and its Trace button come first, where they stay on screen.
        if tools_first {
            let fields = inspected_fields(app, state, &study.id, &run.id, &artifact.name);
            let name = artifact.name.clone();
            flow_tools(app, state, study, run, &name, &fields, focus, ui);
            ui.separator();
        }
        #[cfg(not(target_arch = "wasm32"))]
        {
            if button(ui, "cmd.multiphysics.result.inspect", "multiphysics:inspect_result", !inspection_busy(state)).clicked() {
                start_inspection(app, state, &study.id, &run.id, &artifact.name);
            }
            // A field file is read as soon as it is shown (on a worker thread), so its plots
            // appear without a separate Inspect step.
            let read = state.inspected.as_ref().is_some_and(|cached| cached.key.matches(app, &study.id, &run.id, &artifact.name));
            if !read && !inspection_busy(state) && artifact.name.ends_with(".vtu") {
                start_inspection(app, state, &study.id, &run.id, &artifact.name);
                ui.ctx().request_repaint();
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
                                Ok(grid) => {
                                    state.view.model_surfaces = focus == Node::SurfacePlots;
                                    state.view.model_flow =
                                        matches!(focus, Node::Trajectories | Node::Particles | Node::Isosurfaces).then(|| focus.key());
                                    let docked = shell::feature_docked(state);
                                    state.view.model_flow_docked = docked;
                                    slice_view(app, ui, &mut state.view, grid)
                                }
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
    let flow = matches!(focus, Node::Trajectories | Node::Particles | Node::Isosurfaces | Node::SurfacePlots | Node::SurfaceParams);
    if artifact.name.ends_with(".vtu") && flow && !tools_first {
        let fields = inspected_fields(app, state, &study.id, &run.id, &artifact.name);
        let name = artifact.name.clone();
        flow_tools(app, state, study, run, &name, &fields, focus, ui);
    }
    if artifact.name.ends_with(".vtu") && focus == Node::ResultsMesh {
        mesh_view(app, state, study, run, &artifact.name, ui);
    }
    // The overlay belongs to the flow tools and the mesh view; other nodes show the plain field.
    #[cfg(not(target_arch = "wasm32"))]
    if !flow && focus != Node::ResultsMesh {
        state.view.field_3d.overlay = None;
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

/// The centre of the part picked in the design, in the shown study's result coordinates (m from
/// its domain's corner), from the part's exact bounding box; none without one picked part.
#[cfg(not(target_arch = "wasm32"))]
fn part_centre(app: &mut SolveApp, view: &mut ViewState) -> Option<[f64; 3]> {
    let id = view.selected.clone()?;
    let body = app.session.selection.iter().find_map(|sel| match sel {
        solidcraft_engine::Sel::Body { name } => Some(name.clone()),
        _ => None,
    })?;
    let key = std::sync::Arc::as_ptr(&app.session.doc) as usize;
    if view.anchor_part.as_ref().is_none_or(|(k, b, _)| *k != key || *b != body) {
        let centre = app.session.world_state().ok().and_then(|state| {
            let part = state.bodies.iter().find(|b| b.name == body)?;
            let b = part.body.exact_bounds().ok()?;
            Some([(b.min.x + b.max.x) * 0.5, (b.min.y + b.max.y) * 0.5, (b.min.z + b.max.z) * 0.5])
        });
        view.anchor_part = Some((key, body, centre));
    }
    let centre = view.anchor_part.as_ref().and_then(|(_, _, c)| *c)?;
    let origin = on_model::origin(app, view, &id)?;
    Some([0, 1, 2].map(|a| (centre[a] - origin[a]) / 1000.0))
}

/// A fixed orthographic data chart; it does not install another geometry camera.
#[cfg(not(target_arch = "wasm32"))]
fn slice_view(app: &mut SolveApp, ui: &mut egui::Ui, view: &mut ViewState, grid: &solidcraft_engine::multiphysics::vtk::BoxGrid) {
    let previous_axis = view.slice_axis;
    let model_study = view.selected.as_ref().and_then(|id| app.session.doc.multiphysics.study(id)).is_some_and(|s| s.model.is_some());
    if !model_study {
        view.on_model = false;
    }
    // A surface plot of a study of the model is always on the model.
    let surfaces = model_study && view.model_surfaces;
    if surfaces {
        view.on_model = true;
    }
    // A flow tool of a study of the model draws its trajectories or isosurface on the model.
    if let Some(heading) = view.model_flow.filter(|_| model_study) {
        view.on_model = true;
        ui.heading(tr(heading));
        let overlay = view.field_3d.overlay.clone();
        if overlay.is_none() {
            ui.label(tr("multiphysics.model_flow_hint"));
        }
        let (low, high) = color_controls(ui, view, &overlay.as_ref().map_or([0.0, 1.0], |o| [o.range.0, o.range.1]));
        view.model_flow_range = (low, high);
        if view.model_flow_docked {
            // The graphics window beside this panel shows the model, and the lines go on it.
            if let Some(overlay) = overlay {
                ui.label(tr("multiphysics.model_flow_in_graphics"));
                let segments = overlay.lines.iter().map(|line| line.len().saturating_sub(1)).sum::<usize>() + overlay.triangles.len();
                let unit = (!overlay.lines.is_empty()).then_some("m/s");
                field_legend(ui, view, low, high, unit, segments);
            }
            return;
        }
        let plot_height = (ui.clip_rect().bottom() - ui.next_widget_position().y - 92.0).clamp(320.0, 600.0);
        let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width().max(1.0), plot_height), Sense::hover());
        crate::scenario::publish_rect("multiphysics:field_on_model", rect);
        // Fresh lines are scrolled into sight (the seeds and Trace sit above the view).
        let fresh = overlay.as_ref().map_or(0, |o| std::sync::Arc::as_ptr(o) as usize);
        if fresh != 0 && fresh != view.model_flow_seen {
            view.model_flow_seen = fresh;
            ui.scroll_to_rect(rect, Some(egui::Align::Center));
        }
        let visible = on_model::show(app, ui, view, grid, rect, on_model::Plot::Flow, (low, high));
        if let Some(overlay) = overlay {
            let unit = (!overlay.lines.is_empty()).then_some("m/s");
            field_legend(ui, view, low, high, unit, visible);
        }
        return;
    }
    ui.heading(tr(if surfaces { "multiphysics.tree_surface_plots" } else { "multiphysics.slice_title" }));
    ui.horizontal_wrapped(|ui| {
        if surfaces {
            return;
        }
        for (enabled, label) in [(false, "multiphysics.field_2d"), (true, "multiphysics.field_3d")] {
            let response = crate::ui_kit::choice_chip(ui, label, tr(label), view.show_3d == enabled && !view.on_model);
            if response.clicked() {
                view.show_3d = enabled;
                view.on_model = false;
            }
            crate::scenario::publish_target(if enabled { "multiphysics:view_3d" } else { "multiphysics:view_2d" }, response.rect);
        }
        if model_study {
            let response = crate::ui_kit::choice_chip(ui, "multiphysics.field_on_model", tr("multiphysics.field_on_model"), view.on_model);
            if response.clicked() {
                view.on_model = true;
            }
            crate::scenario::publish_target("multiphysics:view_on_model", response.rect);
            if view.on_model && on_model::vector_field(grid).is_some() {
                let response = crate::ui_kit::checkbox(ui, tr("multiphysics.model_vectors"), &mut view.model_vectors);
                crate::scenario::publish_target("multiphysics:model_vectors", response.rect);
            }
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
        ui.label(tr("multiphysics.slice_through"));
        for (anchor, key, handle, enabled) in [
            (SliceAnchor::Layer, "multiphysics.slice_through_layer", "multiphysics:slice_through:layer", true),
            (SliceAnchor::Part, "multiphysics.slice_through_part", "multiphysics:slice_through:part", model_study),
            (SliceAnchor::Maximum, "multiphysics.slice_through_maximum", "multiphysics:slice_through:maximum", true),
        ] {
            let response = ui.add_enabled_ui(enabled, |ui| crate::ui_kit::choice_chip(ui, key, tr(key), view.slice_anchor == anchor)).inner;
            crate::scenario::publish_target(handle, response.rect);
            if response.clicked() {
                view.slice_anchor = anchor;
            }
        }
        if !model_study && view.slice_anchor == SliceAnchor::Part {
            view.slice_anchor = SliceAnchor::Layer;
        }
        ui.label(tr("multiphysics.slice_layer"));
        let typed = view.slice_anchor == SliceAnchor::Layer;
        let response =
            ui.add_enabled_ui(typed, |ui| crate::ui_kit::number_field(ui, &mut view.slice_layer, 0..=layer_count.saturating_sub(1), 1.0, "")).inner;
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
    // A plot tied to a point goes through the layer nearest it, worked out again each frame so
    // it follows a new run, field or pick.
    let through = match view.slice_anchor {
        SliceAnchor::Layer => None,
        SliceAnchor::Maximum => grid.maximum_at(view.slice_field),
        SliceAnchor::Part => part_centre(app, view),
    };
    match through.and_then(|p| grid.nearest_layer(view.slice_axis, p[view.slice_axis])) {
        Some(layer) => view.slice_layer = layer,
        None if view.slice_anchor == SliceAnchor::Part => {
            ui.label(RichText::new(tr("multiphysics.slice_through_part_hint")).small().color(Tokens::get().text_dim));
        }
        None => {}
    }
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
    let values: Vec<f64> = (0..grid.cells.len()).filter_map(|cell| grid.value(view.slice_field, cell)).collect();
    if !values.iter().any(|v| v.is_finite()) {
        return;
    }
    let (low, high) = color_controls(ui, view, &values);
    let scale = view.color_scale;
    let unit = post::unit_hint(&field.name);
    let t = Tokens::get();
    let visible_width = (ui.clip_rect().right() - ui.next_widget_position().x).max(1.0);
    let plot_height = (ui.clip_rect().bottom() - ui.next_widget_position().y - 92.0).clamp(100.0, 600.0);
    if view.on_model {
        let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width().min(visible_width).max(1.0), plot_height), Sense::hover());
        crate::scenario::publish_rect("multiphysics:field_on_model", rect);
        let plot = if surfaces {
            on_model::Plot::Surfaces
        } else {
            on_model::Plot::Slice { axis: view.slice_axis, along: (horizontal, vertical), coordinate }
        };
        let visible = on_model::show(app, ui, view, grid, rect, plot, (low, high));
        field_legend(ui, view, low, high, unit, visible);
        return;
    }
    if view.show_3d {
        view.field_3d.color_scale = scale;
        view.field_3d.color_range = Some((low, high));
        match crate::multiphysics_field_view::show(app, ui, &mut view.field_3d, grid, view.slice_field, view.slice_axis, coordinate, plot_height) {
            Ok(visible) => field_legend(ui, view, low, high, unit, visible),
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
        let color = scale_color(scale, value, low, high);
        ui.painter().rect_filled(egui::Rect::from_two_pos(position(x[0], y[0]), position(x[1], y[1])), 0.0, color);
        visible += 1;
    }
    field_legend(ui, view, low, high, unit, visible);
}

#[cfg(not(target_arch = "wasm32"))]
fn tick_label(value: f64) -> String {
    if value != 0.0 && !(1e-3..1e4).contains(&value.abs()) { format!("{value:.2e}") } else { format!("{value}") }
}

/// A continuous gradient bar with nice tick values, the range ends and the unit.
#[cfg(not(target_arch = "wasm32"))]
fn field_legend(ui: &mut egui::Ui, view: &ViewState, low: f64, high: f64, unit: Option<&str>, visible: usize) {
    let t = Tokens::get();
    let scale = view.color_scale;
    let legend = ui.vertical(|ui| {
        ui.horizontal_wrapped(|ui| {
            ui.label(tr("multiphysics.slice_value_range")).backing_on_hover_text(tr("multiphysics.slice_units"));
            ui.monospace(format!("{low:.4e} … {high:.4e}"));
            if let Some(unit) = unit {
                ui.monospace(format!("[{unit}]"));
            }
            ui.monospace(format!("{}: {visible}", tr("multiphysics.slice_cells")));
            ui.label(RichText::new(tr(&format!("multiphysics.range_{}", view.range_mode.key()))).small().color(t.text_dim));
        });
        let width = ui.available_width().clamp(120.0, 480.0);
        let (rect, _) = ui.allocate_exact_size(egui::vec2(width, 14.0 + t.font_small + 4.0), Sense::hover());
        let bar = egui::Rect::from_min_size(rect.min, egui::vec2(width, 12.0));
        let steps = 64;
        for step in 0..steps {
            let left = bar.left() + bar.width() * step as f32 / steps as f32;
            let right = bar.left() + bar.width() * (step + 1) as f32 / steps as f32;
            let color = gradient_color(scale, (f64::from(step) + 0.5) / f64::from(steps));
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
    if scale.mark_out_of_range {
        // What the plot draws for values outside the range.
        ui.horizontal(|ui| {
            let (below, above) = scale.map.out_of_range();
            for (shade, key) in [(post::Shade::Below(below), "multiphysics.range_below"), (post::Shade::Above(above), "multiphysics.range_above")] {
                let (rect, _) = ui.allocate_exact_size(egui::vec2(t.space_14, t.space_14), Sense::hover());
                ui.painter().rect_filled(rect, t.control_radius, shade_color(shade));
                ui.label(RichText::new(tr(key)).small());
            }
        });
    }
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
fn gradient_color(scale: post::ColorScale, t: f64) -> egui::Color32 {
    let [r, g, b] = scale.at(t);
    egui::Color32::from_rgb(r, g, b)
}

/// A value's colour on the plots' scale; values out of range, when marked, in the map's own
/// out-of-range colours, which it never shows for a value.
#[cfg(not(target_arch = "wasm32"))]
fn scale_color(scale: post::ColorScale, value: f64, low: f64, high: f64) -> egui::Color32 {
    shade_color(scale.shade(value, low, high))
}

/// Also used by the 3-D field view, which the web build has too.
pub(crate) fn shade_color(shade: post::Shade) -> egui::Color32 {
    let (post::Shade::In([r, g, b]) | post::Shade::Below([r, g, b]) | post::Shade::Above([r, g, b])) = shade;
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

/// Inputs and the last answer of the flow result tools (Flow Trajectories, Particle Studies,
/// Isosurfaces, Surface Plots and Surface Parameters).
#[derive(Default)]
struct FlowTool {
    /// Inputs were given their defaults (from the inspected grid when there is one).
    ready: bool,
    field: String,
    color_field: String,
    start: [f64; 3],
    end: [f64; 3],
    count: usize,
    both: bool,
    max_steps: usize,
    diameter_um: f64,
    density: f64,
    viscosity: f64,
    gravity: bool,
    level: f64,
    /// `None` for the whole boundary, else (axis, maximum side).
    side: Option<(usize, bool)>,
    /// (study, run, artifact, node) the answer belongs to, and the answer.
    answer: Option<((String, String, String, shell::Node), Result<post::QueryResponse, String>)>,
}

const SIDES: [(Option<(usize, bool)>, &str); 7] = [
    (None, "multiphysics.flow_side_all"),
    (Some((0, false)), "X−"),
    (Some((0, true)), "X+"),
    (Some((1, false)), "Y−"),
    (Some((1, true)), "Y+"),
    (Some((2, false)), "Z−"),
    (Some((2, true)), "Z+"),
];

/// The inspected grid's bounds and vector fields, when this artifact has been inspected.
fn inspected_bounds(app: &SolveApp, state: &State, study: &str, run: &str, artifact: &str) -> Option<([[f64; 2]; 3], Vec<String>)> {
    #[cfg(not(target_arch = "wasm32"))]
    if let Some(cached) = &state.inspected
        && cached.key.matches(app, study, run, artifact)
        && let Ok(output) = &cached.outcome
        && let Some(Ok(grid)) = &output.spatial
    {
        let vectors = grid.fields.iter().filter(|field| field.components == 3).map(|field| field.name.clone()).collect();
        return Some((grid.bounds, vectors));
    }
    let _ = (app, state, study, run, artifact);
    None
}

/// The answer as geometry for the 3-D field view.
#[cfg(not(target_arch = "wasm32"))]
fn flow_overlay(result: &post::QueryResult) -> Option<crate::multiphysics_field_view::Overlay> {
    use crate::multiphysics_field_view::Overlay;
    let range = |values: &mut dyn Iterator<Item = f64>| {
        values
            .filter(|v| v.is_finite())
            .fold(None, |r: Option<(f64, f64)>, v| Some(r.map_or((v, v), |(a, b)| (a.min(v), b.max(v)))))
            .unwrap_or((0.0, 1.0))
    };
    match result {
        post::QueryResult::Trajectories { lines } => Some(Overlay {
            lines: lines.iter().map(|line| line.points.iter().copied().zip(line.speeds.iter().copied()).collect()).collect(),
            range: range(&mut lines.iter().flat_map(|line| line.speeds.iter().copied())),
            ..Default::default()
        }),
        post::QueryResult::Isosurface(iso) => {
            let colors = |i: usize| iso.colors.get(i).copied().unwrap_or([iso.level; 3]);
            Some(Overlay {
                triangles: iso.triangles.iter().enumerate().map(|(i, t)| (*t, colors(i))).collect(),
                range: if iso.colors.is_empty() { (iso.level - 1.0, iso.level + 1.0) } else { range(&mut iso.colors.iter().flatten().copied()) },
                ..Default::default()
            })
        }
        post::QueryResult::Surface(surface) => {
            let mut triangles = Vec::new();
            for face in &surface.faces {
                for i in 1..face.corners.len().saturating_sub(1) {
                    if let (Some(a), Some(b), Some(c)) = (face.corners.first(), face.corners.get(i), face.corners.get(i + 1)) {
                        triangles.push(([*a, *b, *c], [face.value; 3]));
                    }
                }
            }
            Some(Overlay { triangles, range: (surface.min, surface.max), hide_slice: true, ..Default::default() })
        }
        _ => None,
    }
}

/// Flow trajectories, particle studies, isosurfaces, surface plots and surface parameters on a
/// retained VTU. The geometry is drawn in the 3-D field view; numbers are listed here.
fn flow_tools(
    app: &mut SolveApp,
    state: &mut State,
    study: &Study,
    run: &solidcraft_engine::multiphysics::RunEvidence,
    artifact: &str,
    fields: &[String],
    focus: shell::Node,
    ui: &mut egui::Ui,
) {
    use shell::Node;
    let inspected = inspected_bounds(app, state, &study.id, &run.id, artifact);
    #[cfg(not(target_arch = "wasm32"))]
    if inspected.is_none() && !inspection_busy(state) {
        start_inspection(app, state, &study.id, &run.id, artifact);
    }
    let tool = &mut state.flow;
    let tracing = matches!(focus, Node::Trajectories | Node::Particles);
    if !tool.ready {
        tool.ready = true;
        tool.count = 10;
        tool.max_steps = 2000;
        tool.diameter_um = 10.0;
        tool.density = 1000.0;
        tool.viscosity = 1.8e-5;
        tool.gravity = true;
        // Seeds across the minimum-x face, where an inlet usually is.
        if let Some((b, _)) = &inspected {
            let inset = |axis: usize, f: f64| b[axis][0] + (b[axis][1] - b[axis][0]) * f;
            tool.start = [inset(0, 0.02), inset(1, 0.1), inset(2, 0.5)];
            tool.end = [inset(0, 0.02), inset(1, 0.9), inset(2, 0.5)];
        }
    }
    let vectors = inspected.as_ref().map(|(_, v)| v.clone()).unwrap_or_default();
    let choices: Vec<String> = if tracing && !vectors.is_empty() { vectors } else { fields.to_vec() };
    if !choices.contains(&tool.field)
        && let Some(first) = choices.first()
    {
        tool.field = first.clone();
    }
    let t = Tokens::get();
    ui.label(RichText::new(tr(focus.key())).strong());
    let current = tool.field.clone();
    let r = crate::ui_kit::row(ui, tr("multiphysics.post_field"), |ui| {
        crate::ui_kit::dropdown_menu(ui, "multiphysics_flow_field", &current, |ui| {
            for name in &choices {
                if crate::ui_kit::list_picker(ui, &format!("multiphysics_flow_field_{name}"), name, &tool.field == name).clicked() {
                    tool.field = name.clone();
                }
            }
        })
    });
    crate::scenario::publish_target("multiphysics:flow_field", r.rect);
    let point = |ui: &mut egui::Ui, label: &str, handle: &str, p: &mut [f64; 3]| {
        let r = crate::ui_kit::row(ui, tr(label), |ui| {
            ui.horizontal(|ui| {
                for value in p.iter_mut() {
                    crate::ui_kit::number_field(ui, value, f64::MIN..=f64::MAX, 1e-3, "");
                }
            })
            .response
        });
        crate::scenario::publish_target(handle, r.rect);
    };
    let query = match focus {
        Node::Trajectories | Node::Particles => {
            point(ui, "multiphysics.flow_seed_start", "multiphysics:flow_seed_start", &mut tool.start);
            point(ui, "multiphysics.flow_seed_end", "multiphysics:flow_seed_end", &mut tool.end);
            ui.horizontal_wrapped(|ui| {
                ui.label(tr("multiphysics.flow_seed_count"));
                crate::ui_kit::number_field(ui, &mut tool.count, 1..=solidcraft_engine::multiphysics::flowpost::MAX_SEEDS, 1.0, "");
                ui.label(tr("multiphysics.flow_max_steps"));
                crate::ui_kit::number_field(ui, &mut tool.max_steps, 1..=solidcraft_engine::multiphysics::flowpost::MAX_STEPS, 10.0, "");
            });
            let particle = if focus == Node::Particles {
                ui.horizontal_wrapped(|ui| {
                    ui.label(tr("multiphysics.flow_diameter"));
                    crate::ui_kit::number_field(ui, &mut tool.diameter_um, 0.001..=1e6, 0.1, " µm");
                    ui.label(tr("multiphysics.flow_density"));
                    crate::ui_kit::number_field(ui, &mut tool.density, 1e-3..=1e5, 1.0, " kg/m³");
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label(tr("multiphysics.flow_viscosity"));
                    crate::ui_kit::number_field(ui, &mut tool.viscosity, 1e-7..=10.0, 1e-6, " Pa·s");
                    let r = crate::ui_kit::checkbox(ui, tr("multiphysics.flow_gravity"), &mut tool.gravity);
                    crate::scenario::publish_target("multiphysics:flow_gravity", r.rect);
                });
                Some(json!({
                    "diameter_m": tool.diameter_um * 1e-6,
                    "density_kg_per_m3": tool.density,
                    "fluid_viscosity_pa_s": tool.viscosity,
                    "gravity_m_per_s2": if tool.gravity { [0.0, -9.81, 0.0] } else { [0.0; 3] },
                }))
            } else {
                let r = crate::ui_kit::checkbox(ui, tr("multiphysics.flow_both"), &mut tool.both);
                crate::scenario::publish_target("multiphysics:flow_both", r.rect);
                None
            };
            json!({
                "kind": "trajectories", "artifact": artifact, "field": tool.field,
                "seeds": {"kind": "line", "start": tool.start, "end": tool.end, "count": tool.count},
                "options": {"max_steps": tool.max_steps, "both_directions": tool.both && particle.is_none(), "particle": particle},
            })
        }
        Node::Isosurfaces => {
            ui.horizontal_wrapped(|ui| {
                ui.label(tr("multiphysics.flow_level"));
                let r = crate::ui_kit::number_field(ui, &mut tool.level, f64::MIN..=f64::MAX, 1e-2, "");
                crate::scenario::publish_target("multiphysics:flow_level", r.rect);
                ui.label(tr("multiphysics.flow_color_by"));
                let shown = if tool.color_field.is_empty() { tr("multiphysics.flow_color_none").to_string() } else { tool.color_field.clone() };
                crate::ui_kit::dropdown_menu(ui, "multiphysics_flow_color", &shown, |ui| {
                    if crate::ui_kit::list_picker(ui, "multiphysics_flow_color_none", tr("multiphysics.flow_color_none"), tool.color_field.is_empty())
                        .clicked()
                    {
                        tool.color_field.clear();
                    }
                    for name in fields {
                        if crate::ui_kit::list_picker(ui, &format!("multiphysics_flow_color_{name}"), name, &tool.color_field == name).clicked() {
                            tool.color_field = name.clone();
                        }
                    }
                });
            });
            let color = (!tool.color_field.is_empty()).then(|| tool.color_field.clone());
            json!({"kind": "isosurface", "artifact": artifact, "field": tool.field, "level": tool.level, "color_field": color})
        }
        _ => {
            ui.horizontal_wrapped(|ui| {
                ui.label(tr("multiphysics.flow_side"));
                for (index, (side, label)) in SIDES.iter().enumerate() {
                    let text = if label.starts_with("multiphysics.") { tr(label) } else { *label };
                    let r = crate::ui_kit::choice_chip(ui, &format!("multiphysics_flow_side_{index}"), text, tool.side == *side);
                    crate::scenario::publish_target(&format!("multiphysics:flow_side:{index}"), r.rect);
                    if r.clicked() {
                        tool.side = *side;
                    }
                }
            });
            let side = tool.side.map(|(axis, max)| json!({"axis": axis, "max": max}));
            json!({"kind": "surface", "artifact": artifact, "field": tool.field, "side": side})
        }
    };
    let key = (study.id.clone(), run.id.clone(), artifact.to_string(), focus);
    let response = crate::ui_kit::button(
        ui,
        tr(if tracing { "multiphysics.flow_trace" } else { "multiphysics.flow_compute" }),
        crate::ui_kit::Button::Enabled(!tool.field.is_empty()),
    );
    crate::scenario::publish_target("multiphysics:flow_run", response.rect);
    if response.clicked() {
        let project = app.session.doc.multiphysics.clone();
        let answer = serde_json::from_value::<post::Query>(query).map_err(|error| error.to_string()).and_then(|query| {
            post::query(&project, &post::QueryRequest { id: study.id.clone(), run_id: run.id.clone(), query }).map_err(|e| e.to_string())
        });
        #[cfg(not(target_arch = "wasm32"))]
        {
            state.view.field_3d.overlay = answer.as_ref().ok().and_then(|a| flow_overlay(&a.result)).map(std::sync::Arc::new);
            if state.view.field_3d.overlay.is_some() {
                state.view.show_3d = true;
                // A study of the model shows the lines on the model itself.
                state.view.on_model = study.model.is_some();
            }
        }
        state.flow.answer = Some((key.clone(), answer));
    }
    let tool = &state.flow;
    match &tool.answer {
        Some((owner, Ok(answer))) if *owner == key => {
            field_answer(ui, answer);
            // A study of the model draws the answer on the model, not in the 3-D field view.
            if study.model.is_none() {
                ui.small(tr("multiphysics.flow_in_view"));
            }
        }
        Some((owner, Err(error))) if *owner == key => {
            ui.colored_label(t.form_error, error);
        }
        _ => {}
    }
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
                    // A narrow property panel must not break the names mid-word.
                    ui.add(egui::Label::new(tr(label)).wrap_mode(egui::TextWrapMode::Extend));
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
        post::QueryResult::Trajectories { lines } => {
            use solidcraft_engine::multiphysics::flowpost::TrajectoryEnd;
            let count = |end: TrajectoryEnd| lines.iter().filter(|line| line.end == end).count();
            let times: Vec<f64> = lines.iter().filter(|line| line.end == TrajectoryEnd::LeftMesh).map(|line| line.residence_time).collect();
            egui::Grid::new("multiphysics_trajectories").show(ui, |ui| {
                for (label, value) in [
                    ("multiphysics.flow_lines", lines.len().to_string()),
                    ("multiphysics.flow_left_mesh", count(TrajectoryEnd::LeftMesh).to_string()),
                    ("multiphysics.flow_stagnated", count(TrajectoryEnd::Stagnated).to_string()),
                    ("multiphysics.flow_max_steps_reached", count(TrajectoryEnd::MaxSteps).to_string()),
                    ("multiphysics.flow_seed_outside", count(TrajectoryEnd::SeedOutside).to_string()),
                ] {
                    // A narrow property panel must not break the names mid-word.
                    ui.add(egui::Label::new(tr(label)).wrap_mode(egui::TextWrapMode::Extend));
                    ui.monospace(value);
                    ui.end_row();
                }
                if !times.is_empty() {
                    let mean = times.iter().sum::<f64>() / times.len() as f64;
                    let max = times.iter().copied().fold(0.0, f64::max);
                    ui.label(tr("multiphysics.flow_residence"));
                    ui.add(
                        egui::Label::new(RichText::new(format!("{mean:.4e} s (max {max:.4e} s)")).monospace()).wrap_mode(egui::TextWrapMode::Extend),
                    );
                    ui.end_row();
                }
            });
        }
        post::QueryResult::Isosurface(iso) => {
            ui.monospace(format!("{} = {:.6e}{unit}", iso.field, iso.level));
            ui.monospace(format!(
                "{}: {}   {}: {:.6e}",
                tr("multiphysics.flow_triangles"),
                iso.triangles.len(),
                tr("multiphysics.flow_area"),
                iso.area
            ));
            if iso.truncated {
                ui.colored_label(t.warning, tr("multiphysics.flow_truncated"));
            }
        }
        post::QueryResult::Surface(surface) => {
            egui::Grid::new("multiphysics_surface").show(ui, |ui| {
                let mut rows = vec![
                    ("multiphysics.flow_area", surface.area),
                    ("multiphysics.post_mean", surface.mean),
                    ("multiphysics.post_minimum", surface.min),
                    ("multiphysics.post_maximum", surface.max),
                    ("multiphysics.post_integral", surface.integral),
                ];
                if let Some(flux) = surface.flux {
                    rows.push(("multiphysics.flow_flux", flux));
                }
                for (label, value) in rows {
                    // A narrow property panel must not break the names mid-word.
                    ui.add(egui::Label::new(tr(label)).wrap_mode(egui::TextWrapMode::Extend));
                    ui.monospace(format!("{value:.6e}"));
                    ui.end_row();
                }
            });
            ui.label(RichText::new(format!("{} {}", surface.faces.len(), tr("multiphysics.flow_faces"))).small().color(t.text_dim));
        }
        post::QueryResult::Goals(_) | post::QueryResult::Components { .. } => {}
    }
}

/// The colour map (reversed, in bands, out-of-range marking) and how its range is chosen:
/// from the data, locked so frames, views and runs share one scale, clipped to percentiles, or
/// centred on zero. `values` are the plotted field's.
#[cfg(not(target_arch = "wasm32"))]
fn color_controls(ui: &mut egui::Ui, view: &mut ViewState, values: &[f64]) -> (f64, f64) {
    let data = post::range_for(post::RangeMode::Auto, values, None).unwrap_or((0.0, 1.0));
    ui.horizontal_wrapped(|ui| {
        ui.label(tr("multiphysics.post_color_map"));
        let current = view.color_scale.map;
        crate::ui_kit::dropdown_menu(ui, "multiphysics_color_map", current.name(), |ui| {
            for map in post::ColorMap::ALL {
                if crate::ui_kit::list_picker(ui, &format!("multiphysics_color_map_{}", map.name()), map.name(), current == map).clicked() {
                    view.color_scale.map = map;
                }
            }
        });
        let response = crate::ui_kit::checkbox(ui, tr("multiphysics.range_reversed"), &mut view.color_scale.reversed);
        crate::scenario::publish_target("multiphysics:color_reversed", response.rect);
        ui.label(tr("multiphysics.range_bands"));
        let response = crate::ui_kit::number_field(ui, &mut view.color_scale.bands, 0..=post::MAX_BANDS, 1.0, "")
            .backing_on_hover_text(tr("multiphysics.range_bands_hint"));
        crate::scenario::publish_target("multiphysics:color_bands", response.rect);
        if view.color_scale.bands == 1 {
            view.color_scale.bands = 2;
        }
        let response = crate::ui_kit::checkbox(ui, tr("multiphysics.range_mark"), &mut view.color_scale.mark_out_of_range)
            .backing_on_hover_text(tr("multiphysics.range_mark_hint"));
        crate::scenario::publish_target("multiphysics:color_mark_out_of_range", response.rect);
    });
    ui.horizontal_wrapped(|ui| {
        ui.label(tr("multiphysics.range_mode"));
        for mode in post::RangeMode::ALL {
            let key = mode.key();
            let response = crate::ui_kit::choice_chip(
                ui,
                &format!("multiphysics_range_{key}"),
                tr(&format!("multiphysics.range_{key}")),
                view.range_mode == mode,
            );
            crate::scenario::publish_target(&format!("multiphysics:range_mode:{key}"), response.rect);
            if response.clicked() && view.range_mode != mode {
                // Locking starts from the range shown now.
                if mode == post::RangeMode::Locked && view.locked_range.is_none() {
                    view.locked_range = post::range_for(view.range_mode, values, None);
                }
                view.range_mode = mode;
            }
        }
        if view.range_mode == post::RangeMode::Locked {
            let (low, high) = view.locked_range.get_or_insert(data);
            let step = (data.1 - data.0).abs().max(1e-12) * 1e-3;
            let response = crate::ui_kit::number_field(ui, low, f64::MIN..=f64::MAX, step, "");
            crate::scenario::publish_target("multiphysics:range_low", response.rect);
            ui.label("…");
            let response = crate::ui_kit::number_field(ui, high, f64::MIN..=f64::MAX, step, "");
            crate::scenario::publish_target("multiphysics:range_high", response.rect);
        }
    });
    post::range_for(view.range_mode, values, view.locked_range).unwrap_or(data)
}

#[cfg(test)]
mod tests {
    /// A harness with the workspace registered, as the app registers it at start-up.
    fn harness() -> Harness {
        crate::register();
        Harness::new()
    }

    /// The workspace's state in the app (a fresh one when it has none yet).
    fn ws(app: &mut SolveApp) -> &mut State {
        if peek(app).is_none() {
            put_state(app, State::default());
        }
        app.multiphysics_workspace.0.as_mut().and_then(|state| state.downcast_mut::<State>()).expect("the workspace's state")
    }

    fn wsr(app: &SolveApp) -> &State {
        peek(app).expect("the workspace's state")
    }

    use super::*;
    use crate::scenario::Harness;
    use std::sync::Arc;

    #[test]
    fn case_view_drags_domain_faces_and_sets_picked_faces_and_parts() {
        use solidcraft_engine::multiphysics::case_scene::{PatchRole, patch_roles};
        let mut h = harness();
        let reply = h.app.run("multiphysics.study.create", json!({"template":"flow-electronics-enclosure-cht","name":"Picked"})).unwrap();
        let id = reply["study"]["id"].as_str().unwrap().to_string();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        let draft = |h: &Harness| serde_json::from_str::<Value>(&wsr(&h.app).drafts[&id].case_text).unwrap();
        // Computational domain: the panel's + button and a drag on the arrow both move the
        // right face outward by whole cells, and nothing reaches the document until Apply.
        let document = h.app.session.doc.clone();
        h.step(&json!({"click":{"handle":"multiphysics:menu:project"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:menu:project:domain"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:domain_out:right"}})).unwrap();
        assert_eq!(draft(&h)["mesh"]["dimensions"][0], 21);
        h.frames(2);
        let base = h.at(&json!({"handle":"multiphysics:domain_arrow_base:right"})).unwrap();
        let tip = h.at(&json!({"handle":"multiphysics:domain_arrow_tip:right"})).unwrap();
        let (dx, dy) = (tip[0] - base[0], tip[1] - base[1]);
        let length = (dx * dx + dy * dy).sqrt().max(1.0);
        let end = [tip[0] + dx / length * 150.0, tip[1] + dy / length * 150.0];
        h.step(&json!({"drag":[tip, end], "steps": 12})).unwrap();
        let grown = draft(&h)["mesh"]["dimensions"][0].as_u64().unwrap();
        assert!(grown > 21, "dragging the arrow outward adds cells: {grown}");
        assert!(Arc::ptr_eq(&document, &h.app.session.doc));
        // Boundary conditions: pick a domain face in the graphics window and make it a wall.
        h.step(&json!({"click":{"handle":"multiphysics:menu:insert"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:menu:insert:boundary"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:case_face:top"}})).unwrap();
        let picked = ws(&mut h.app).case_view.patches();
        assert_eq!(picked.len(), 1, "one face picked");
        h.step(&json!({"click":{"handle":"multiphysics:case_set_symmetry"}})).unwrap();
        assert_eq!(patch_roles(&draft(&h))[picked[0]], PatchRole::Symmetry);
        // Solid materials: pick a part and give it a conductivity.
        h.step(&json!({"click":{"handle":"multiphysics:menu:insert"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:menu:insert:solid_materials"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:case_part:processor_block"}})).unwrap();
        let parts = ws(&mut h.app).case_view.parts();
        assert_eq!(parts.len(), 1, "one part picked");
        material_library_text(&mut h, "multiphysics:case_conductivity", "150");
        h.step(&json!({"click":{"handle":"multiphysics:case_apply_material"}})).unwrap();
        let case = draft(&h);
        let row = case["thermal"]["solids"].as_array().unwrap().iter().find(|row| row["obstacle"] == parts[0].as_str()).cloned().unwrap();
        assert_eq!(row["conductivity_w_per_m_k"], json!([150.0, 150.0, 150.0]));
        assert!(Arc::ptr_eq(&document, &h.app.session.doc));
    }

    #[test]
    fn reference_menus_insert_apply_and_clone_without_running() {
        let mut h = harness();
        let reply = h.app.run("multiphysics.study.create", json!({"template":"flow-channel","name":"Reviewed flow"})).unwrap();
        let id = reply["study"]["id"].as_str().unwrap().to_string();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        let document = h.app.session.doc.clone();
        let camera = h.app.cam;
        let history = h.app.session.undo.len();
        h.step(&json!({"click":{"handle":"multiphysics:menu:insert"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:menu:insert:leds"}})).unwrap();
        assert!(Arc::ptr_eq(&document, &h.app.session.doc), "unavailable model must be inert");
        h.step(&json!({"key":"Escape"})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:menu:insert"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:menu:insert:porous"}})).unwrap();
        material_library_text(&mut h, "multiphysics:feature_name", "reviewed_filter");
        h.step(&json!({"click":{"handle":"multiphysics:feature_add"}})).unwrap();
        assert!(Arc::ptr_eq(&document, &h.app.session.doc));
        assert_eq!(history, h.app.session.undo.len());
        material_library_text(&mut h, "multiphysics:case_field:/porous_blocks/0/max/0", "2");
        let case: Value = serde_json::from_str(&wsr(&h.app).drafts[&id].case_text).unwrap();
        assert_eq!(case["porous_blocks"][0]["max"][0], 2, "scalar array fields must be editable");
        assert_eq!(case["porous_blocks"][0]["name"], "reviewed_filter");
        h.step(&json!({"click":{"handle":"multiphysics:apply"}})).unwrap();
        assert_eq!(h.app.session.doc.multiphysics.study(&id).unwrap().case, case);
        assert_eq!(h.app.session.undo.len(), history + 1);
        assert!(h.app.session.multiphysics_jobs().is_empty());
        assert_eq!(h.app.cam, camera);
        let evidence = std::env::var("SOLIDCRAFT_MULTIPHYSICS_EVIDENCE").ok();
        for (w, height) in [(1440., 900.), (1280., 800.)] {
            render_frame(&mut h.app, w, height, evidence.as_deref(), "reference-porous");
        }
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:menu:project"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:menu:project:clone"}})).unwrap();
        assert_eq!(h.app.session.doc.multiphysics.studies.len(), 2);
        assert_eq!(h.app.session.doc.multiphysics.studies[1].case, case);
        assert!(h.app.session.doc.multiphysics.studies[1].runs.is_empty());
        assert_eq!(h.app.session.undo.len(), history + 2);
        h.app.run("edit.undo", json!({})).unwrap();
        assert_eq!(h.app.session.doc.multiphysics.studies.len(), 1);
    }

    #[test]
    fn reference_seven_step_wizard_preserves_edits_and_creates_atomically() {
        let mut h = harness();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        let document = h.app.session.doc.clone();
        let history = h.app.session.undo.len();
        h.step(&json!({"click":{"handle":"multiphysics:menu:project"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:menu:project:wizard"}})).unwrap();
        material_library_text(&mut h, "multiphysics:wizard_name", "Reviewed seven step flow");
        h.step(&json!({"click":{"handle":"multiphysics:wizard_next"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:wizard_category:Flow"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:wizard_template:flow-channel"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:wizard_next"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:wizard_next"}})).unwrap();
        material_library_text(&mut h, "multiphysics:case_field:/method", "simple");
        h.step(&json!({"click":{"handle":"multiphysics:wizard_back"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:wizard_next"}})).unwrap();
        for _ in 0..3 {
            h.step(&json!({"click":{"handle":"multiphysics:wizard_next"}})).unwrap();
        }
        assert_eq!(crate::scenario::count_of("multiphysics:wizard_step"), Some(6.0));
        assert!(Arc::ptr_eq(&document, &h.app.session.doc));
        assert_eq!(history, h.app.session.undo.len());
        let evidence = std::env::var("SOLIDCRAFT_MULTIPHYSICS_EVIDENCE").ok();
        for (w, height) in [(1440., 900.), (1280., 800.)] {
            render_frame(&mut h.app, w, height, evidence.as_deref(), "reference-wizard");
        }
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:wizard_finish"}})).unwrap();
        assert!(ws(&mut h.app).wizard.is_none());
        assert_eq!(h.app.session.doc.multiphysics.studies.len(), 1);
        assert_eq!(h.app.session.doc.multiphysics.studies[0].case["method"], "simple");
        assert_eq!(h.app.session.doc.multiphysics.studies[0].name, "Reviewed seven step flow");
        assert_eq!(h.app.session.undo.len(), history + 1);
        h.app.run("edit.undo", json!({})).unwrap();
        assert!(h.app.session.doc.multiphysics.studies.is_empty());
    }

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
        crate::register();
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
        state.view.selected = Some(id.clone());
        state.view.stage = Stage::Results;
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
        ws(&mut app).inspection = Some(InspectionJob { key, cancelled: cancelled.clone(), done });
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
        crate::register();
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
            let mut output =
                ctx.run_ui(egui::RawInput { screen_rect: Some(bounds), time: Some(frame as f64 / 60.0), events, ..Default::default() }, |ui| {
                    Container::for_layout(LayoutMode::Docked, Surface::Viewport).central().show(ui, |ui| {
                        crate::multiphysics_field_view::show(&mut app, ui, &mut state, &grid, 0, 2, 0.5, 400.0).unwrap();
                    });
                });
            // The frame's font textures are not uploaded anywhere here; discard them explicitly.
            output.textures_delta.clear();
            assert_eq!(app.cam, camera);
            assert_eq!(app.pending_fit, pending);
        }
        assert!(state.camera.distance < before_zoom, "the preserved pending fit must coexist with actual field zoom");
    }

    #[test]
    fn material_library_mouse_import_edit_undo_and_reopen_preserve_study() {
        let mut h = harness();
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
        assert_eq!(ws(&mut h.app).materials.selected_design.as_deref(), Some("Enclosure copper reference"));
        h.step(&json!({"click":{"handle":"multiphysics:undo"}})).unwrap();
        assert!(h.app.session.doc.user_materials.is_empty());
        assert!(ws(&mut h.app).materials.selected_design.is_none());
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
        let mut h = harness();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:materials"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:materials_design"}})).unwrap();
        h.step(&json!({"click":{"handle":"material:new"}})).unwrap();
        assert!(crate::scenario::handle_at("material:save").is_some());
        material_library_text(&mut h, "material:field:name", "Unsaved old-document material");
        h.app.session = solidcraft_engine::Session::default();
        h.frames(3);
        assert!(ws(&mut h.app).materials.selected_design.is_none());
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
    fn a_study_of_the_model_is_set_up_on_the_model_in_the_3d_view() {
        let mut h = harness();
        // A square duct along x: 100 mm long, 40 × 40 mm outside, 4 mm walls, open at both ends.
        h.app.run("solid.box", json!({"length":100,"width":40,"height":40,"center":[50,20,20],"body_name":"Duct"})).unwrap();
        h.app.run("solid.box", json!({"length":120,"width":32,"height":32,"center":[50,20,20],"operation":"cut"})).unwrap();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:model_study"}})).unwrap();
        h.frames(3);
        let study = h.app.session.doc.multiphysics.studies.last().cloned().expect("a study of the model");
        assert!(study.model.is_some());
        // The graphics window is the design's 3-D view, not the case drawing.
        assert!(crate::scenario::rect_of("multiphysics:case_show_cad").is_none());
        assert!(crate::scenario::rect_of("multiphysics:model_rebuild").is_some());
        // Select the duct's end face at x = 0 in the view and make it an inlet.
        h.app.session.selection =
            vec![solidcraft_engine::Sel::Face { body: "Duct".into(), index: 0, point: solidcore_record::Vec3::new(0.0, 20.0, 2.0) }];
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:model_add_condition"}})).unwrap();
        h.frames(3);
        let study = h.app.session.doc.multiphysics.study(&study.id).cloned().unwrap();
        let link = study.model.clone().unwrap();
        assert_eq!(link.openings.len(), 1, "{:?}", ws(&mut h.app).error);
        assert_eq!(link.openings[0].name, "inlet 1");
        let left = study.case["velocity_boundaries"].as_array().into_iter().flatten().find(|r| r["patch"] == "left").cloned().unwrap();
        assert_eq!(left["value"]["x"], json!(1.0));
        // With the same face selected: smaller cells around it (the Mesh node), and a goal on its
        // flow (the Goals node).
        h.step(&json!({"click":{"handle":"multiphysics:mesh"}})).unwrap();
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:model_refine"}})).unwrap();
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:goals"}})).unwrap();
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:model_add_goal"}})).unwrap();
        h.frames(3);
        let link = h.app.session.doc.multiphysics.study(&study.id).unwrap().model.clone().unwrap();
        assert_eq!(link.refinements.len(), 1, "{:?}", ws(&mut h.app).error);
        assert_eq!(link.goals.len(), 1, "{:?}", ws(&mut h.app).error);
        assert_eq!(link.goals[0].name, "goal_1");
        let study = h.app.session.doc.multiphysics.study(&study.id).cloned().unwrap();
        assert_eq!(study.case["goals"][0]["target"], json!("left"));
        // Editing the model makes the case out of date, and Rebuild brings it back.
        h.app.run("solid.box", json!({"length":10,"width":10,"height":10,"center":[50,60,20],"operation":"new"})).unwrap();
        h.frames(2);
        assert_eq!(h.app.run("multiphysics.model.status", json!({"id":study.id})).unwrap()["current"], false);
        h.step(&json!({"click":{"handle":"multiphysics:general"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:model_rebuild"}})).unwrap();
        h.frames(2);
        assert_eq!(h.app.run("multiphysics.model.status", json!({"id":study.id})).unwrap()["current"], true);
    }

    /// Heat transfer in a study of the model, set up by clicks: heat on in Solid Materials, the
    /// duct's material and power on the selected body, a heat source on the selected face, a
    /// temperature on the inlet's end face, and the duct's highest temperature as a goal.
    #[test]
    fn heat_in_a_study_of_the_model_is_set_on_its_bodies_and_faces() {
        let mut h = harness();
        h.app.run("solid.box", json!({"length":100,"width":40,"height":40,"center":[50,20,20],"body_name":"Duct"})).unwrap();
        h.app.run("solid.box", json!({"length":120,"width":32,"height":32,"center":[50,20,20],"operation":"cut"})).unwrap();
        let reply = h.app.run("multiphysics.model.create", json!({"name":"Duct heat","cell_mm":4.0})).unwrap();
        let id = reply["study"]["id"].as_str().unwrap().to_string();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        let link = |h: &Harness| h.app.session.doc.multiphysics.study(&id).and_then(|s| s.model.clone()).unwrap();
        let face = |point: [f64; 3]| {
            vec![solidcraft_engine::Sel::Face { body: "Duct".into(), index: 0, point: solidcore_record::Vec3::new(point[0], point[1], point[2]) }]
        };
        // Solid Materials is listed for a study of the model, and turns heat on.
        h.step(&json!({"click":{"handle":"multiphysics:tree:tree_solid_materials"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:model_heat"}})).unwrap();
        h.frames(2);
        assert!(link(&h).heat.is_some(), "{:?}", ws(&mut h.app).error);
        // The duct's material and power, on the selected body.
        h.app.session.selection = vec![solidcraft_engine::Sel::Body { name: "Duct".into() }];
        h.frames(2);
        material_library_text(&mut h, "multiphysics:model_conductivity", "200");
        material_library_text(&mut h, "multiphysics:model_power", "1");
        h.step(&json!({"click":{"handle":"multiphysics:model_apply_material"}})).unwrap();
        h.frames(2);
        let materials = link(&h).heat.unwrap().materials;
        assert_eq!((materials.len(), materials[0].conductivity_w_per_m_k, materials[0].power_w), (1, 200.0, 1.0), "{:?}", ws(&mut h.app).error);
        // 5 W on the channel's wall at y = 4 mm.
        h.app.session.selection = face([50.0, 4.0, 20.0]);
        h.frames(2);
        material_library_text(&mut h, "multiphysics:model_source_power", "5");
        h.step(&json!({"click":{"handle":"multiphysics:model_add_heat_source"}})).unwrap();
        h.frames(2);
        let sources = link(&h).heat.unwrap().face_sources;
        assert_eq!((sources.len(), sources[0].power_w), (1, 5.0), "{:?}", ws(&mut h.app).error);
        // The inlet's end face held at 37 °C.
        h.app.session.selection = face([0.0, 20.0, 2.0]);
        h.frames(2);
        material_library_text(&mut h, "multiphysics:model_thermal_value", "37");
        h.step(&json!({"click":{"handle":"multiphysics:model_add_thermal_face"}})).unwrap();
        h.frames(2);
        let conditions = link(&h).heat.unwrap().face_conditions;
        assert_eq!(conditions.len(), 1, "{:?}", ws(&mut h.app).error);
        let case = h.app.session.doc.multiphysics.study(&id).unwrap().case.clone();
        let held = case["thermal"]["boundaries"][0]["condition"]["temperature_k"].as_f64().unwrap();
        assert!((held - 310.15).abs() < 1e-9, "{case}");
        assert_eq!(case["thermal"]["boundaries"][0]["patch"], "left");
        // The duct's highest temperature, on the selected body, in the Goals node.
        h.app.session.selection = vec![solidcraft_engine::Sel::Body { name: "Duct".into() }];
        h.step(&json!({"click":{"handle":"multiphysics:goals"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:model_goal_kind_menu"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:model_goal_kind:body_max_temperature"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:model_add_goal"}})).unwrap();
        h.frames(2);
        let goals = link(&h).goals;
        assert_eq!(goals.len(), 1, "{:?}", ws(&mut h.app).error);
        assert_eq!(goals[0].body.as_deref(), Some("Duct"));
        let case = h.app.session.doc.multiphysics.study(&id).unwrap().case.clone();
        assert_eq!(case["goals"][0]["kind"], "solid_max_temperature");
    }

    /// The plots' range follows its mode: the data's, 2–98 %, around zero, or locked (starting
    /// from the range shown when it was locked).
    #[test]
    fn color_controls_follow_the_range_mode() {
        let ctx = crate::fonts::test_context();
        let mut view = ViewState::default();
        let values: Vec<f64> = (0..=100).map(f64::from).chain([1e6]).collect();
        let range = |view: &mut ViewState| {
            let mut out = (0.0, 0.0);
            let mut output = ctx.run_ui(egui::RawInput::default(), |ui| out = color_controls(ui, view, &values));
            output.textures_delta.clear();
            out
        };
        assert_eq!(range(&mut view), (0.0, 1e6));
        view.range_mode = post::RangeMode::Percentile;
        let (low, high) = range(&mut view);
        assert!(low >= 1.0 && high <= 100.0, "{low} {high}");
        view.range_mode = post::RangeMode::Symmetric;
        assert_eq!(range(&mut view), (-1e6, 1e6));
        view.range_mode = post::RangeMode::Locked;
        view.locked_range = Some((10.0, 20.0));
        assert_eq!(range(&mut view), (10.0, 20.0));
        // Marked out-of-range values take the map's out-of-range colours, not its end colours.
        view.color_scale.mark_out_of_range = true;
        let (below, above) = view.color_scale.map.out_of_range();
        let rgb = |[r, g, b]: [u8; 3]| egui::Color32::from_rgb(r, g, b);
        assert_eq!(scale_color(view.color_scale, 5.0, 10.0, 20.0), rgb(below));
        assert_eq!(scale_color(view.color_scale, 25.0, 10.0, 20.0), rgb(above));
        assert_ne!(scale_color(view.color_scale, 15.0, 10.0, 20.0), rgb(above));
    }

    /// Evidence: the colour controls, a field drawn in 8 bands on a locked range with values
    /// outside it marked, and its legend, rendered offscreen. Written only with
    /// SOLIDCRAFT_MULTIPHYSICS_EVIDENCE set.
    #[test]
    fn color_scale_controls_and_legend_render() {
        let ctx = crate::fonts::test_context();
        let mut view = ViewState {
            color_scale: post::ColorScale { map: post::ColorMap::CoolWarm, reversed: false, bands: 8, mark_out_of_range: true },
            range_mode: post::RangeMode::Locked,
            locked_range: Some((-0.6, 0.6)),
            ..Default::default()
        };
        let values: Vec<f64> = (0..400).map(|i| ((i % 20) as f64 / 19.0 * 2.0 - 1.0) * (1.0 - (i / 20) as f64 / 40.0)).collect();
        let bounds = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(720.0, 420.0));
        let mut textures = std::collections::HashMap::new();
        let mut last = None;
        for frame in 0..3 {
            let mut output =
                ctx.run_ui(egui::RawInput { screen_rect: Some(bounds), time: Some(f64::from(frame) / 60.0), ..Default::default() }, |ui| {
                    let (low, high) = color_controls(ui, &mut view, &values);
                    let (rect, _) = ui.allocate_exact_size(egui::vec2(400.0, 200.0), egui::Sense::hover());
                    for (i, value) in values.iter().enumerate() {
                        let (x, y) = ((i % 20) as f32, (i / 20) as f32);
                        let cell = egui::Rect::from_min_size(rect.min + egui::vec2(x * 20.0, y * 10.0), egui::vec2(20.0, 10.0));
                        ui.painter().rect_filled(cell, 0.0, scale_color(view.color_scale, *value, low, high));
                    }
                    field_legend(ui, &view, low, high, Some("K"), values.len());
                });
            crate::ui_kit::test_textures_update(&mut textures, &output.textures_delta);
            output.textures_delta.clear();
            last = Some(output);
        }
        let output = last.unwrap();
        // Out-of-range cells are drawn, in the neutral colours, and the legend says what they are.
        assert!(crate::scenario::rect_of("multiphysics:field_legend").is_some());
        if let Ok(folder) = std::env::var("SOLIDCRAFT_MULTIPHYSICS_EVIDENCE") {
            std::fs::create_dir_all(&folder).unwrap();
            let image = crate::ui_kit::test_raster(&ctx, output.shapes, output.pixels_per_point, bounds, &textures);
            image.save(std::path::Path::new(&folder).join("multiphysics-color-scale.png")).unwrap();
        }
    }

    /// The Parts folder: the tree's part rows, the component table and the 3-D view share the
    /// design's selection, and the selection hides, isolates and ghosts through commands.
    #[test]
    fn parts_in_the_tree_and_table_share_the_design_selection() {
        let mut h = harness();
        for (name, x) in [("U1", 0.0), ("U2", 30.0)] {
            h.app.run("solid.box", json!({"length":10,"width":10,"height":10,"center":[x,0,0],"operation":"new","body_name":name})).unwrap();
        }
        h.app.run("multiphysics.model.create", json!({"name":"Board","analysis":"external","margin_mm":10.0,"cell_mm":5.0})).unwrap();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        let selected = |h: &Harness| super::parts_ui::selected(&h.app);
        // A part row in the tree picks the body in the design, as the 3-D view would.
        h.step(&json!({"click":{"handle":"multiphysics:tree_part:U2"}})).unwrap();
        h.frames(2);
        assert_eq!(selected(&h), ["U2"], "{:?}", ws(&mut h.app).error);
        // The Parts panel's table rows show and change the same selection.
        h.step(&json!({"click":{"handle":"multiphysics:tree:tree_parts"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:parts_row:U1"}})).unwrap();
        h.frames(2);
        assert_eq!(selected(&h), ["U1"]);
        // A pick in the 3-D view is what the tree and the table show.
        h.app.session.selection = vec![solidcraft_engine::Sel::Body { name: "U2".into() }];
        h.frames(2);
        // Ghost draws U2 see-through, All Solid brings it back; neither is an undo step.
        let undo = h.app.session.undo.len();
        h.step(&json!({"click":{"handle":"multiphysics:parts_ghost"}})).unwrap();
        h.frames(2);
        let ghosted: Vec<String> = h.app.session.ghosted_bodies.iter().cloned().collect();
        assert_eq!(ghosted, ["U2"], "{:?}", ws(&mut h.app).error);
        h.step(&json!({"click":{"handle":"multiphysics:parts_solid_all"}})).unwrap();
        h.frames(2);
        assert!(h.app.session.ghosted_bodies.is_empty());
        assert_eq!(h.app.session.undo.len(), undo);
        // Isolate hides the other part in the design; Show All shows it again.
        h.step(&json!({"click":{"handle":"multiphysics:parts_isolate"}})).unwrap();
        h.frames(2);
        assert!(!h.app.body_visible("U1") && h.app.body_visible("U2"), "{:?}", ws(&mut h.app).error);
        h.step(&json!({"click":{"handle":"multiphysics:parts_show_all"}})).unwrap();
        h.frames(2);
        assert!(h.app.body_visible("U1"));
        // Colour by material: neither part has one yet, so both stand out.
        h.step(&json!({"click":{"handle":"multiphysics:parts_colour:material"}})).unwrap();
        h.frames(2);
        let swatches = h.app.session.part_swatches.clone();
        assert_eq!(swatches.len(), 2, "{:?}", ws(&mut h.app).error);
        assert!(swatches.iter().all(|(_, s)| *s == solidcraft_engine::multiphysics::setup_check::Swatch::Muted), "heat is off: {swatches:?}");
        // Grade plot: a legend entry shows only its parts. By tag: U1 is an "ic", U2 has none.
        h.app.run("multiphysics.parts.tag", json!({"bodies":["U1"],"add_tags":["ic"]})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:parts_colour:tag"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:parts_legend:ic"}})).unwrap();
        h.frames(2);
        let ghosted: Vec<String> = h.app.session.ghosted_bodies.iter().cloned().collect();
        assert_eq!((ghosted, selected(&h)), (vec!["U2".to_string()], vec!["U1".to_string()]), "{:?}", ws(&mut h.app).error);
        // Again: every part solid.
        h.step(&json!({"click":{"handle":"multiphysics:parts_legend:ic"}})).unwrap();
        h.frames(2);
        assert!(h.app.session.ghosted_bodies.is_empty());
        // Find: the example rule reads results, which this study has none of; a tag rule finds U1,
        // and Show only these ghosts the rest.
        h.step(&json!({"click":{"handle":"multiphysics:find_run"}})).unwrap();
        h.frames(2);
        assert!(crate::scenario::handles("multiphysics:find_show_only").is_empty(), "no run: nothing found");
        h.step(&json!({"click":{"handle":"multiphysics:find_remove:0"}})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:find_add"}})).unwrap();
        h.frames(2);
        material_library_text(&mut h, "multiphysics:find_value:0", "ic");
        assert!(crate::scenario::handles("multiphysics:find_value:1").is_empty(), "the example rule does not come back");
        h.app.run("multiphysics.parts.select", json!({})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:find_run"}})).unwrap();
        h.frames(2);
        assert_eq!(selected(&h), ["U1"], "{:?}", ws(&mut h.app).error);
        h.step(&json!({"click":{"handle":"multiphysics:find_show_only"}})).unwrap();
        h.frames(2);
        assert_eq!(h.app.session.ghosted_bodies.iter().collect::<Vec<_>>(), ["U2"]);
        h.app.run("multiphysics.parts.ghost", json!({"only":true})).unwrap();
        // Group: U1 and U2 are tagged ic by their designators; Select similar picks both cubes
        // from one; a size limit leaves out what is smaller (nothing here) by the study's own query.
        h.step(&json!({"click":{"handle":"multiphysics:group_designators"}})).unwrap();
        h.frames(2);
        let tags = |h: &Harness, body: &str| h.app.session.doc.multiphysics.parts.parts.get(body).map(|p| p.tags.iter().cloned().collect::<Vec<_>>());
        assert_eq!((tags(&h, "U1"), tags(&h, "U2")), (Some(vec!["ic".to_string()]), Some(vec!["ic".to_string()])), "{:?}", ws(&mut h.app).error);
        h.app.run("multiphysics.parts.select", json!({"bodies":["U2"]})).unwrap();
        h.frames(2);
        h.step(&json!({"click":{"handle":"multiphysics:group_similar"}})).unwrap();
        h.frames(2);
        assert_eq!(selected(&h), ["U1", "U2"]);
        material_library_text(&mut h, "multiphysics:group_small_mm", "5");
        h.step(&json!({"click":{"handle":"multiphysics:group_small_apply"}})).unwrap();
        h.frames(2);
        let link = h.app.session.doc.multiphysics.studies.iter().find_map(|s| s.model.clone()).unwrap();
        assert_eq!(link.selection.as_deref(), Some("Board without small parts"), "{:?}", ws(&mut h.app).error);
        h.step(&json!({"click":{"handle":"multiphysics:parts_colour:material"}})).unwrap();
        h.frames(2);
        // The batch report waits for a run with a thermal report.
        assert!(crate::scenario::rect_of("multiphysics:report_write").is_some());
        // The checklist is listed, and its goals line has no parts to select.
        assert!(crate::scenario::handles("multiphysics:parts_check:goals").len() == 1);
        // Close shuts the Parts panel; the colours stay.
        h.step(&json!({"click":{"handle":"multiphysics:panel_close"}})).unwrap();
        h.frames(2);
        assert!(!ws(&mut h.app).dialog_open);
        assert_eq!(h.app.session.part_swatches.len(), 2);
        // Leaving the workspace draws ghosted and coloured parts as before.
        h.app.run("multiphysics.parts.ghost", json!({"bodies":["U1"]})).unwrap();
        switch(&mut h.app, "design").unwrap();
        assert!(h.app.session.ghosted_bodies.is_empty() && h.app.session.part_swatches.is_empty());
    }

    /// The Mesh node of a study of the model: the grid drawn on the model, the level of the
    /// initial mesh, and refinement where a kept run's velocity steps.
    #[test]
    fn the_mesh_of_a_study_of_the_model_takes_a_level_and_adapts_to_a_run() {
        use solidcraft_engine::multiphysics as mp;
        let mut h = harness();
        h.app.run("solid.box", json!({"length":100,"width":40,"height":40,"center":[50,20,20],"body_name":"Duct"})).unwrap();
        h.app.run("solid.box", json!({"length":120,"width":32,"height":32,"center":[50,20,20],"operation":"cut"})).unwrap();
        let reply = h.app.run("multiphysics.model.create", json!({"name":"Duct mesh","cell_mm":4.0})).unwrap();
        let id = reply["study"]["id"].as_str().unwrap().to_string();
        // A kept run whose velocity steps up from x cell 12 on (4 mm cells, 25 × 10 × 10).
        let study = h.app.session.doc.multiphysics.study(&id).unwrap().clone();
        let (mut points, mut connectivity, mut offsets, mut velocity) = (String::new(), String::new(), String::new(), String::new());
        let mut count = 0;
        for k in 0..10 {
            for j in 0..10 {
                for i in 0..25 {
                    for (dx, dy, dz) in [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)] {
                        points.push_str(&format!("{} {} {} ", (i + dx) as f64 * 0.004, (j + dy) as f64 * 0.004, (k + dz) as f64 * 0.004));
                    }
                    for v in 0..8 {
                        connectivity.push_str(&format!("{} ", count * 8 + v));
                    }
                    count += 1;
                    offsets.push_str(&format!("{} ", count * 8));
                    velocity.push_str(if i < 12 { "0 0 0 " } else { "1 0 0 " });
                }
            }
        }
        let vtu = format!(
            r#"<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian"><UnstructuredGrid><Piece NumberOfPoints="{}" NumberOfCells="{count}"><CellData><DataArray type="Float64" Name="velocity_m_per_s" NumberOfComponents="3" format="ascii">{velocity}</DataArray></CellData><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{points}</DataArray></Points><Cells><DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray><DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray><DataArray type="UInt8" Name="types" format="ascii">{}</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>"#,
            count * 8,
            "12 ".repeat(count)
        );
        let capabilities: Vec<&str> = mp::airflow::requirements(&study.case).iter().map(|r| r.capability).collect();
        let manifest = json!({
            "schema":"cfd-solvers-rs/run/1", "rust_port_version":mp::SOLVER_VERSION, "upstream_revision":mp::SOLVER_REVISION,
            "backend":"native_rust_serial_cpu", "solver":"collocated_flow3d", "complete_clone":false, "execution_completed":true,
            "linear_solver_converged":null, "steady_state_converged":null, "capabilities":capabilities, "case":study.case,
            "metrics":{"wall_time_s":0.01}, "files":["flow_000010.vtu","manifest.json"]
        });
        let evidence = mp::RunEvidence {
            id: "first".into(),
            study_revision: study.revision,
            case_hash: mp::case_hash(&study.case).unwrap(),
            input_case: study.case.clone(),
            artifacts: vec![
                mp::Artifact { name: "flow_000010.vtu".into(), text: vtu },
                mp::Artifact { name: "manifest.json".into(), text: manifest.to_string() },
            ],
            manifest,
        };
        let project = &mut Arc::make_mut(&mut h.app.session.doc).multiphysics;
        Arc::make_mut(project).attach_evidence(&id, evidence).unwrap();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        ws(&mut h.app).selected = Some(id.clone());
        h.step(&json!({"click":{"handle":"multiphysics:mesh"}})).unwrap();
        h.frames(2);
        let link = |h: &Harness| h.app.session.doc.multiphysics.study(&id).and_then(|s| s.model.clone()).unwrap();
        // The grid on the domain's three low faces: 26 + 11 + 11 node planes, each on two faces.
        h.step(&json!({"click":{"handle":"multiphysics:model_show_grid"}})).unwrap();
        h.frames(2);
        assert_eq!(crate::scenario::count_of("multiphysics:model_grid_lines"), Some(2.0 * (26.0 + 11.0 + 11.0)));
        // Refinement where the run's velocity steps: one box around x cells 11 and 12.
        h.step(&json!({"click":{"handle":"multiphysics:model_adapt"}})).unwrap();
        h.frames(2);
        let refinements = link(&h).refinements;
        assert_eq!(refinements.len(), 1, "{:?}", ws(&mut h.app).error);
        assert_eq!(refinements[0].name, "adapt_1");
        h.step(&json!({"click":{"handle":"multiphysics:model_adapt_remove"}})).unwrap();
        h.frames(2);
        assert!(link(&h).refinements.is_empty(), "{:?}", ws(&mut h.app).error);
        // A click on the level slider sets a level and drops the typed 4 mm cells.
        h.step(&json!({"click":{"handle":"multiphysics:model_mesh_level"}})).unwrap();
        h.frames(2);
        let link = link(&h);
        assert!(link.mesh_level.is_some() && link.cell_mm.is_none(), "{:?} {:?}", link, ws(&mut h.app).error);
    }

    /// A cut plot and a surface plot of a study of the model drawn on the model, and a click probing it.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn results_of_a_study_of_the_model_are_drawn_and_probed_on_the_model() {
        use solidcraft_engine::multiphysics as mp;
        let mut h = harness();
        h.app.run("solid.box", json!({"length":100,"width":40,"height":40,"center":[50,20,20],"body_name":"Duct"})).unwrap();
        h.app.run("solid.box", json!({"length":120,"width":32,"height":32,"center":[50,20,20],"operation":"cut"})).unwrap();
        let reply = h.app.run("multiphysics.model.create", json!({"name":"Duct","cell_mm":4.0})).unwrap();
        let study = h.app.session.doc.multiphysics.study(reply["study"]["id"].as_str().unwrap()).unwrap().clone();
        // Two cells across the duct's box (case metres), each with one speed.
        let vtu = r#"<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian"><UnstructuredGrid><Piece NumberOfPoints="12" NumberOfCells="2"><CellData><DataArray type="Float64" Name="speed_m_per_s" NumberOfComponents="1" format="ascii">0.1 0.2</DataArray><DataArray type="Float64" Name="velocity_m_per_s" NumberOfComponents="3" format="ascii">0.1 0 0 0.2 0.02 0</DataArray></CellData><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 0.05 0 0 0.1 0 0 0 0.04 0 0.05 0.04 0 0.1 0.04 0 0 0 0.04 0.05 0 0.04 0.1 0 0.04 0 0.04 0.04 0.05 0.04 0.04 0.1 0.04 0.04</DataArray></Points><Cells><DataArray type="Int64" Name="connectivity" format="ascii">0 1 4 3 6 7 10 9 1 2 5 4 7 8 11 10</DataArray><DataArray type="Int64" Name="offsets" format="ascii">8 16</DataArray><DataArray type="UInt8" Name="types" format="ascii">12 12</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>"#;
        let capabilities: Vec<&str> = mp::airflow::requirements(&study.case).iter().map(|r| r.capability).collect();
        let manifest = json!({
            "schema":"cfd-solvers-rs/run/1", "rust_port_version":mp::SOLVER_VERSION, "upstream_revision":mp::SOLVER_REVISION,
            "backend":"native_rust_serial_cpu", "solver":"collocated_flow3d", "complete_clone":false, "execution_completed":true,
            "linear_solver_converged":null, "steady_state_converged":null, "capabilities":capabilities, "case":study.case,
            "metrics":{"wall_time_s":0.01}, "files":["flow.vtu","manifest.json"]
        });
        let evidence = mp::RunEvidence {
            id: "duct-run".into(),
            study_revision: study.revision,
            case_hash: mp::case_hash(&study.case).unwrap(),
            input_case: study.case.clone(),
            artifacts: vec![
                mp::Artifact { name: "flow.vtu".into(), text: vtu.into() },
                mp::Artifact { name: "manifest.json".into(), text: manifest.to_string() },
            ],
            manifest,
        };
        let project = &mut Arc::make_mut(&mut h.app.session.doc).multiphysics;
        Arc::make_mut(project).attach_evidence(&study.id, evidence).unwrap();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(2);
        let mut state = take_state(&mut h.app);
        state.selected = Some(study.id.clone());
        state.selected_run = Some("duct-run".into());
        state.shown_run = Some("duct-run".into());
        state.selected_artifact = 0;
        shell::select(&mut state, shell::Node::CutPlots);
        start_inspection(&h.app, &mut state, &study.id, "duct-run", "flow.vtu");
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(20);
        while state.inspection.is_some() {
            poll_inspection(&h.app, &mut state);
            assert!(std::time::Instant::now() < deadline, "inspection worker did not finish");
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        state.view.on_model = true;
        put_state(&mut h.app, state);
        h.frames(3);
        // The plot area is the 3-D view, and both cells of the slice are drawn on the model.
        assert!(crate::scenario::rect_of("multiphysics:field_on_model").is_some(), "{:?}", ws(&mut h.app).error);
        assert_eq!(crate::scenario::count_of("multiphysics:slice_on_model_cells"), Some(2.0));
        assert!(crate::scenario::rect_of("multiphysics:view_on_model").is_some());

        // Velocity arrows go on the cut plot, one per cell of the slice.
        h.step(&json!({"click":{"handle":"multiphysics:model_vectors"}})).unwrap();
        h.frames(2);
        assert!(ws(&mut h.app).view.model_vectors);
        assert_eq!(crate::scenario::count_of("multiphysics:model_vectors_drawn"), Some(2.0));

        // Through the picked part: the plane goes through the duct's centre, taken from its
        // exact box and placed in the case's metres; with nothing picked it says how.
        h.app.run("multiphysics.parts.select", json!({})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:slice_through:part"}})).unwrap();
        h.frames(2);
        assert_eq!(ws(&mut h.app).view.slice_anchor, SliceAnchor::Part);
        let mut state = take_state(&mut h.app);
        assert_eq!(part_centre(&mut h.app, &mut state.view), None, "nothing picked");
        put_state(&mut h.app, state);
        h.app.run("multiphysics.parts.select", json!({"bodies":["Duct"]})).unwrap();
        h.frames(2);
        let mut state = take_state(&mut h.app);
        let centre = part_centre(&mut h.app, &mut state.view).expect("the duct's centre");
        put_state(&mut h.app, state);
        for (got, want) in centre.iter().zip([0.05, 0.02, 0.02]) {
            assert!((got - want).abs() < 1e-9, "{centre:?}");
        }
        // Through the field's maximum: the faster of the two cells along X is the second, and the
        // plot moves there; it stays tied to it as the plane changes back.
        h.step(&json!({"click":{"handle":"multiphysics:slice_plane:0"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:slice_through:maximum"}})).unwrap();
        h.frames(2);
        assert_eq!((ws(&mut h.app).view.slice_anchor, ws(&mut h.app).view.slice_layer), (SliceAnchor::Maximum, 1));
        h.step(&json!({"click":{"handle":"multiphysics:slice_through:layer"}})).unwrap();
        h.step(&json!({"click":{"handle":"multiphysics:slice_plane:2"}})).unwrap();
        h.frames(2);
        assert_eq!(ws(&mut h.app).view.slice_anchor, SliceAnchor::Layer);

        // A click on the slice probes the cell under it.
        h.step(&json!({"click":{"handle":"multiphysics:slice_probe_target"}})).unwrap();
        h.frames(2);
        let probe = ws(&mut h.app).view.model_probe.clone().expect("a probe from the click on the slice");
        assert!(probe.cell.is_some(), "{probe:?}");
        assert!(crate::scenario::count_of("multiphysics:model_probe_value").is_some());
        ws(&mut h.app).view.model_probe = None;

        // The surface plot colours the duct's wetted inner walls, and a click on one probes it.
        shell::select(ws(&mut h.app), shell::Node::SurfacePlots);
        h.frames(3);
        let faces = crate::scenario::count_of("multiphysics:surface_on_model_faces").unwrap_or(0.0);
        assert!(faces > 0.0, "{:?}", ws(&mut h.app).error);
        let surface = ws(&mut h.app).view.model_surface.clone().unwrap();
        // The outer faces and the open ends are dry; only the four inner walls are wetted.
        assert!(surface.wet_count() > 0 && surface.dry_count() > 0, "{surface:?}");
        h.step(&json!({"click":{"handle":"multiphysics:model_probe_target"}})).unwrap();
        h.frames(2);
        let probe = ws(&mut h.app).view.model_probe.clone().expect("a probe from the click");
        assert!(probe.cell.is_some(), "{probe:?}");
        let value = crate::scenario::count_of("multiphysics:model_probe_value").unwrap();
        assert!(value == 0.1 || value == 0.2, "{value}");
        assert_eq!(ws(&mut h.app).view.probe_point, probe.case_m);

        // Flow trajectories traced from the default seeds are drawn on the model (in the graphics
        // window when the tool is in the property panel beside it), not the cut plot.
        shell::select(ws(&mut h.app), shell::Node::Trajectories);
        h.frames(3);
        h.step(&json!({"click":{"handle":"multiphysics:flow_run"}})).unwrap();
        h.frames(3);
        let view = &wsr(&h.app).view;
        assert!(view.on_model && view.model_flow.is_some());
        let overlay = view.field_3d.overlay.clone().expect("trajectories from the trace");
        assert!(!overlay.lines.is_empty());
        let drawn = crate::scenario::count_of("multiphysics:model_flow_drawn").unwrap_or(0.0);
        assert!(drawn > 0.0, "{drawn}");
    }

    #[test]
    fn a_result_row_opens_its_own_file_when_the_run_is_first_shown() {
        use solidcraft_engine::multiphysics as mp;
        let mut h = harness();
        let reply = h.app.run("multiphysics.study.create", json!({"template":"poisson1d","name":"Rows"})).unwrap();
        let study = h.app.session.doc.multiphysics.study(reply["study"]["id"].as_str().unwrap()).unwrap().clone();
        // A PVD series, as flow runs write, puts the Frame stepper on the cut plot too.
        let files = ["series.csv", "field_000000.vtu", "field_000010.vtu", "field.pvd", "manifest.json"];
        let pvd = r#"<?xml version="1.0"?><VTKFile type="Collection" version="0.1"><Collection><DataSet timestep="0" part="0" file="field_000000.vtu"/><DataSet timestep="1" part="0" file="field_000010.vtu"/></Collection></VTKFile>"#;
        let manifest = json!({
            "schema":"cfd-solvers-rs/run/1", "rust_port_version":mp::SOLVER_VERSION, "upstream_revision":mp::SOLVER_REVISION,
            "backend":"native_rust_serial_cpu", "solver":"poisson1d", "complete_clone":false, "execution_completed":true,
            "linear_solver_converged":null, "steady_state_converged":null, "case":study.case, "metrics":{"wall_time_s":0.01}, "files":files
        });
        let artifacts = files
            .iter()
            .map(|name| mp::Artifact {
                name: (*name).into(),
                text: match *name {
                    "manifest.json" => manifest.to_string(),
                    "field.pvd" => pvd.into(),
                    _ => "x,phi\n0,0\n1,0\n".into(),
                },
            })
            .collect();
        let evidence = mp::RunEvidence {
            id: "rows-run".into(),
            study_revision: study.revision,
            case_hash: mp::case_hash(&study.case).unwrap(),
            input_case: study.case.clone(),
            manifest,
            artifacts,
        };
        let project = &mut Arc::make_mut(&mut h.app.session.doc).multiphysics;
        Arc::make_mut(project).attach_evidence(&study.id, evidence).unwrap();
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(3);
        // Clicking the second VTU row of a run the view has not drawn yet must show that file,
        // not the run's first CSV or VTU (the view resets its file when it first shows a run).
        let row = "multiphysics:tree_result:field_000010.vtu";
        if crate::scenario::rect_of(row).is_none() {
            h.step(&json!({"click":{"handle":format!("multiphysics:study:{}", study.id)}})).unwrap();
            h.frames(3);
        }
        h.step(&json!({"click":{"handle":row}})).unwrap();
        h.frames(3);
        assert_eq!(ws(&mut h.app).selected_artifact, 2);
    }

    #[test]
    fn workspace_renders_bounded_native_frames_at_desktop_sizes() {
        crate::register();
        let mut app = SolveApp::new(Default::default(), Default::default());
        let reply = app.run("multiphysics.study.create", json!({"template":"poisson1d","name":"Poisson 1D · source = 2"})).unwrap();
        let study = app.session.doc.multiphysics.study(reply["study"]["id"].as_str().unwrap()).unwrap().clone();
        switch(&mut app, "multiphysics").unwrap();
        ws(&mut app).selected = Some(study.id.clone());
        ws(&mut app).stage = Stage::Physics;
        ws(&mut app).drafts.insert(study.id.clone(), Draft::from_study(&study));
        ws(&mut app).drafts.get_mut(&study.id).unwrap().question = "Check the scalar field under a uniform source".into();
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
            ws(&mut app).stage = Stage::Results;
            ws(&mut app).drafts.insert(study.id.clone(), Draft::from_study(&study));
            ws(&mut app).selected_artifact =
                app.session.doc.multiphysics.studies[0].runs[0].artifacts.iter().position(|artifact| artifact.name.ends_with(".csv")).unwrap_or(0);
            for (width, height) in [(1440., 900.), (1280., 800.)] {
                render_frame(&mut app, width, height, evidence_dir.as_deref(), "results");
            }
            let committed = app.session.doc.multiphysics.study(&study.id).unwrap().clone();
            let mut state = take_state(&mut app);
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
            put_state(&mut app, state);
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
        crate::register();
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
        let mut state = take_state(&mut app);
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
        put_state(&mut app, state);
        let before_case = app.session.doc.multiphysics.study(&id).unwrap().case.clone();
        let mut harness = harness();
        harness.app = app;
        harness.frames(3);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_plane:1"}})).unwrap();
        assert_eq!(ws(&mut harness.app).view.slice_axis, 1);
        assert_eq!(ws(&mut harness.app).view.slice_layer, 0);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_layer"},"double":true})).unwrap();
        harness.step(&json!({"key":"A","ctrl":true})).unwrap();
        harness.step(&json!({"text":"2"})).unwrap();
        harness.step(&json!({"key":"Enter"})).unwrap();
        assert_eq!(ws(&mut harness.app).view.slice_layer, 2);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_field_menu"}})).unwrap();
        harness.step(&json!({"click":{"handle":format!("multiphysics:slice_field:{pressure_field}")}})).unwrap();
        assert_eq!(ws(&mut harness.app).view.slice_field, pressure_field);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_field_menu"}})).unwrap();
        harness.step(&json!({"click":{"handle":format!("multiphysics:slice_field:{velocity_field}")}})).unwrap();
        assert_eq!(ws(&mut harness.app).view.slice_field, velocity_field);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_plane:0"}})).unwrap();
        assert_eq!(ws(&mut harness.app).view.slice_axis, 0);
        assert_eq!(ws(&mut harness.app).view.slice_layer, 0);
        harness.step(&json!({"click":{"handle":"multiphysics:slice_plane:2"}})).unwrap();
        assert_eq!(ws(&mut harness.app).view.slice_axis, 2);
        assert_eq!(harness.app.session.doc.multiphysics.study(&id).unwrap().case, before_case, "view controls must not edit solver inputs");
        ws(&mut harness.app).view.slice_layer = 2;
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
        assert!(ws(&mut harness.app).view.show_3d);
        let before_orbit = ws(&mut harness.app).view.field_3d.camera;
        let at = crate::scenario::rect_of("multiphysics:field_slice").unwrap().center();
        let orbit = harness.call("ui.drag", json!({"x0":at.x,"y0":at.y,"x1":at.x+50.0,"y1":at.y+24.0,"button":"middle","shift":true,"steps":8}));
        assert_eq!(orbit["ok"], true, "{orbit}");
        let after_orbit = ws(&mut harness.app).view.field_3d.camera;
        assert!(before_orbit.yaw != after_orbit.yaw || before_orbit.pitch != after_orbit.pitch, "3D orbit must move its own camera");
        let pan = harness.call("ui.drag", json!({"x0":at.x,"y0":at.y,"x1":at.x+30.0,"y1":at.y-14.0,"button":"middle","steps":8}));
        assert_eq!(pan["ok"], true, "{pan}");
        let after_pan = ws(&mut harness.app).view.field_3d.camera;
        assert_ne!(after_pan.target, after_orbit.target, "3D pan must move the field view");
        let zoom = harness.call("ui.scroll", json!({"x":at.x,"y":at.y,"delta":120.0}));
        assert_eq!(zoom["ok"], true, "{zoom}");
        assert!(ws(&mut harness.app).view.field_3d.camera.distance < after_pan.distance, "3D wheel must zoom its field camera");
        let before_glyphs = ws(&mut harness.app).view.field_3d.glyphs;
        harness.step(&json!({"click":{"handle":"multiphysics:3d_glyphs"}})).unwrap();
        assert_ne!(ws(&mut harness.app).view.field_3d.glyphs, before_glyphs);
        let glyph_count = crate::scenario::count_of("multiphysics:3d_glyph_count").unwrap();
        assert!(glyph_count > 0.0 && glyph_count <= 128.0, "vector sampling must be visible and bounded: {glyph_count}");
        let before_fit = ws(&mut harness.app).view.field_3d.camera;
        harness.step(&json!({"click":{"handle":"multiphysics:3d_fit"}})).unwrap();
        let fitted = ws(&mut harness.app).view.field_3d.camera;
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
        let mut h = harness();
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
        put_state(&mut h.app, state);
        close(&mut h.app);
        let mut state = take_state(&mut h.app);
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
        crate::register();
        let mut app = SolveApp::new(Default::default(), Default::default());
        let old_uid = app.session.uid;
        let (send, done) = std::sync::mpsc::channel();
        ws(&mut app).export = Some(ExportJob { uid: old_uid, command: "multiphysics.case.export", done });
        app.run("file.new", json!({"name":"New document"})).unwrap();
        assert_ne!(app.session.uid, old_uid);
        assert!(ws(&mut app).export_busy());
        let log = app.session.log.clone();
        send.send(Ok(json!({"files":["old-case.json"]}))).unwrap();
        tick(&mut app, &crate::fonts::test_context());
        assert!(!ws(&mut app).export_busy());
        assert!(ws(&mut app).notice.is_none());
        assert!(ws(&mut app).error.is_none());
        assert_eq!(app.session.log, log, "an old completion must not echo into another document");
        assert!(app.session.doc.multiphysics.is_empty());
    }

    #[test]
    fn document_tabs_keep_same_id_study_drafts_with_their_document() {
        let mut h = harness();
        let first = h.app.run("multiphysics.study.create", json!({"template":"poisson1d","name":"First document study"})).unwrap();
        let first_id = first["study"]["id"].as_str().unwrap().to_string();
        let first_uid = h.app.session.uid;
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(1);
        ws(&mut h.app).drafts.get_mut(&first_id).unwrap().case_text = "{unfinished first document".into();
        let mut second = solidcraft_engine::Session::default();
        let added = second.execute("multiphysics.study.create", &json!({"template":"poisson2d","name":"Second document study"})).unwrap();
        assert_eq!(added["study"]["id"], first_id, "the identity collision exercises document ownership");
        crate::documents::open_tab(&mut h.app, second).unwrap();
        let second_uid = h.app.session.uid;
        assert_ne!(first_uid, second_uid);
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(1);
        let draft = ws(&mut h.app).drafts.get_mut(&first_id).unwrap();
        assert!(serde_json::from_str::<Value>(&draft.case_text).is_ok());
        draft.question = "Second document draft".into();
        assert_eq!(ws(&mut h.app).y_column, 1);
        crate::documents::switch(&mut h.app, 0);
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(1);
        assert_eq!(h.app.session.uid, first_uid);
        assert_eq!(ws(&mut h.app).drafts.get(&first_id).unwrap().case_text, "{unfinished first document");
        crate::documents::switch(&mut h.app, 1);
        switch(&mut h.app, "multiphysics").unwrap();
        h.frames(1);
        assert_eq!(h.app.session.uid, second_uid);
        assert_eq!(ws(&mut h.app).drafts.get(&first_id).unwrap().question, "Second document draft");
        let clean = solidcraft_engine::Session::new(h.app.session.doc.as_ref().clone());
        crate::register();
        let app = SolveApp::new(clean, Default::default());
        assert!(!crate::documents::active_is_blank(&app), "a study-only document contains authored work even without CAD features");
    }

    #[test]
    fn workspace_pointer_creates_study_and_retains_invalid_draft() {
        let mut h = harness();
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
        assert_eq!(ws(&mut h.app).drafts.get(&id).unwrap().case_text, "{invalid json");
        assert!(ws(&mut h.app).error.is_some());
        h.step(&json!({"click":{"handle":"multiphysics:cancel_draft"}})).unwrap();
        assert!(!ws(&mut h.app).drafts.get(&id).unwrap().dirty());
        let camera = h.app.cam;
        switch(&mut h.app, "design").unwrap();
        assert_eq!(h.app.cam, camera);
        switch(&mut h.app, "multiphysics").unwrap();
        assert_eq!(h.app.session.doc.multiphysics.studies.len(), 1);
    }
}
