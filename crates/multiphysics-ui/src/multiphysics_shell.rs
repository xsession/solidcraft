//! The Multiphysics workspace laid out as a CAD-embedded flow-analysis environment: a ribbon of
//! grouped commands across the top, a resource bar with the Analysis tree (Input Data and
//! Results folders of the active project) on the left, the CAD graphics window in the centre,
//! and feature dialogs, result windows, the solver monitor and the project wizard floating over
//! the graphics window. Everything here routes to the workspace's existing editors and commands;
//! it adds no numerical behaviour of its own.

use super::*;
use crate::ribbon;
use egui::{Rect, pos2, vec2};

#[path = "multiphysics_menus.rs"]
mod menus;

/// Width of the resource bar's icon strip.
pub(super) const RESOURCE_BAR: f32 = 34.0;
/// Width kept for the workspace switcher at the left of the ribbon.
const SWITCHER: f32 = 112.0;
/// The ribbon: its tab row and one row of labelled groups.
pub(super) const RIBBON: f32 = menus::HEIGHT + ribbon::TAB_ROW_HEIGHT + ribbon::GROUP_HEIGHT + 1.0;

/// One entry of the Analysis tree. Folders and features map onto the workspace stages, so
/// scripts that select a stage and people who click the tree reach the same editors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Node {
    Project,
    Input,
    General,
    Control,
    Domain,
    ComponentControl,
    FluidSubdomains,
    SolidMaterials,
    Boundary,
    Fans,
    Sources,
    Radiation,
    RadiationSources,
    HeatSinks,
    ThermalJoints,
    FromComponents,
    Porous,
    Perforated,
    Contact,
    Initial,
    Rotating,
    TwoResistor,
    HeatPipes,
    Pcb,
    Tec,
    Leds,
    Electrical,
    Batteries,
    Transferred,
    Goals,
    Mesh,
    /// A study of the model's parts, linked to the design's selection.
    Parts,
    Coupling,
    Results,
    ResultsMesh,
    CutPlots,
    SurfacePlots,
    Isosurfaces,
    Trajectories,
    Particles,
    PointParams,
    SurfaceParams,
    VolumeParams,
    XyPlots,
    GoalPlots,
    Report,
    ExportResults,
    Animation,
    Compare,
    Verify,
}

impl Node {
    /// The Input data folder's features in the reference flow-analysis tree's order
    /// (Computational Domain … Initial Conditions, Goals, Mesh, From Components), then the
    /// SolidCraft-only folders. Computational domain, component control, fluid subdomains,
    /// boundary conditions, goals and mesh show by default; the others once the case has
    /// entries for them or Customize Tree turns them on.
    pub(super) const INPUT: &[Self] = &[
        Self::Domain,
        Self::ComponentControl,
        Self::FluidSubdomains,
        Self::SolidMaterials,
        Self::Boundary,
        Self::Transferred,
        Self::Fans,
        Self::Sources,
        Self::Radiation,
        Self::RadiationSources,
        Self::Contact,
        Self::Tec,
        Self::HeatSinks,
        Self::Porous,
        Self::Perforated,
        Self::TwoResistor,
        Self::HeatPipes,
        Self::ThermalJoints,
        Self::Pcb,
        Self::Electrical,
        Self::Initial,
        Self::Rotating,
        Self::Leds,
        Self::Batteries,
        Self::Goals,
        Self::Mesh,
        Self::Parts,
        Self::FromComponents,
        Self::Coupling,
    ];
    /// The Results folder's features in the reference tree's order (Mesh … Report, Export
    /// Results), then SolidCraft's Animation and Compare. Unbuilt tools are listed greyed when
    /// Customize Tree shows them, so the tree says what is missing.
    const RESULTS: &[Self] = &[
        Self::ResultsMesh,
        Self::CutPlots,
        Self::SurfacePlots,
        Self::Isosurfaces,
        Self::Trajectories,
        Self::Particles,
        Self::PointParams,
        Self::SurfaceParams,
        Self::VolumeParams,
        Self::XyPlots,
        Self::GoalPlots,
        Self::Report,
        Self::ExportResults,
        Self::Animation,
        Self::Compare,
    ];

    pub(super) fn key(self) -> &'static str {
        match self {
            Self::Project => "multiphysics.tree_project",
            Self::Input => "multiphysics.tree_input",
            Self::General => "multiphysics.tree_general",
            Self::Control => "multiphysics.tree_control",
            Self::Domain => "multiphysics.tree_domain",
            Self::ComponentControl => "multiphysics.tree_component_control",
            Self::RadiationSources => "multiphysics.tree_radiation_sources",
            Self::HeatSinks => "multiphysics.tree_heat_sinks",
            Self::ThermalJoints => "multiphysics.tree_thermal_joints",
            Self::FromComponents => "multiphysics.tree_from_components",
            Self::ResultsMesh => "multiphysics.tree_results_mesh",
            Self::ExportResults => "multiphysics.tree_export_results",
            Self::FluidSubdomains => "multiphysics.tree_fluid_subdomains",
            Self::SolidMaterials => "multiphysics.tree_solid_materials",
            Self::Boundary => "multiphysics.tree_boundary",
            Self::Fans => "multiphysics.tree_fans",
            Self::Sources => "multiphysics.tree_sources",
            Self::Radiation => "multiphysics.tree_radiation",
            Self::Porous => "multiphysics.tree_porous",
            Self::Perforated => "multiphysics.tree_perforated",
            Self::Contact => "multiphysics.tree_contact",
            Self::Initial => "multiphysics.tree_initial",
            Self::Rotating => "multiphysics.tree_rotating",
            Self::TwoResistor => "multiphysics.tree_two_resistor",
            Self::HeatPipes => "multiphysics.tree_heat_pipes",
            Self::Pcb => "multiphysics.tree_pcb",
            Self::Tec => "multiphysics.tree_tec",
            Self::Leds => "multiphysics.tree_leds",
            Self::Electrical => "multiphysics.tree_electrical",
            Self::Batteries => "multiphysics.tree_batteries",
            Self::Transferred => "multiphysics.tree_transferred",
            Self::Goals => "multiphysics.tree_goals",
            Self::Mesh => "multiphysics.tree_mesh",
            Self::Parts => "multiphysics.tree_parts",
            Self::Coupling => "multiphysics.tree_coupling",
            Self::Results => "multiphysics.tree_results",
            Self::CutPlots => "multiphysics.tree_cut_plots",
            Self::SurfacePlots => "multiphysics.tree_surface_plots",
            Self::Isosurfaces => "multiphysics.tree_isosurfaces",
            Self::Trajectories => "multiphysics.tree_trajectories",
            Self::Particles => "multiphysics.tree_particles",
            Self::PointParams => "multiphysics.tree_point_params",
            Self::SurfaceParams => "multiphysics.tree_surface_params",
            Self::VolumeParams => "multiphysics.tree_volume_params",
            Self::XyPlots => "multiphysics.tree_xy_plots",
            Self::GoalPlots => "multiphysics.tree_goal_plots",
            Self::Report => "multiphysics.tree_report",
            Self::Animation => "multiphysics.tree_animation",
            Self::Compare => "multiphysics.tree_compare",
            Self::Verify => "multiphysics.tree_verify",
        }
    }

    fn icon(self) -> &'static str {
        match self {
            Self::Project => "configuration_table",
            Self::Input | Self::Results => "folder",
            Self::General => "settings",
            Self::Control => "timeline_play",
            Self::Domain => "box",
            Self::ComponentControl => "component",
            Self::RadiationSources => "render_light_point",
            Self::HeatSinks => "pattern",
            Self::ThermalJoints => "contact_sets",
            Self::FromComponents => "import",
            Self::ResultsMesh => "mesh_body",
            Self::ExportResults => "export",
            Self::FluidSubdomains => "render_environment",
            Self::SolidMaterials => "material",
            Self::Boundary => "surface",
            Self::Fans => "spun",
            Self::Sources => "render_light_point",
            Self::Radiation => "render_sun",
            Self::Porous => "perforate",
            Self::Perforated => "dots",
            Self::Contact => "contact_sets",
            Self::Initial => "restore",
            Self::Rotating => "revolve",
            Self::TwoResistor => "component",
            Self::HeatPipes => "pipe",
            Self::Pcb => "sheet",
            Self::Tec => "render_light_area",
            Self::Leds => "render_light_spot",
            Self::Electrical => "drive",
            Self::Batteries => "rect",
            Self::Transferred => "import",
            Self::Goals => "checkpoint",
            Self::Mesh => "mesh_body",
            Self::Parts => "body",
            Self::Coupling => "joint",
            Self::CutPlots => "section",
            Self::SurfacePlots => "surface",
            Self::Isosurfaces => "iso_analysis",
            Self::Trajectories => "spline",
            Self::Particles => "dots",
            Self::PointParams => "point",
            Self::SurfaceParams => "plane",
            Self::VolumeParams => "body",
            Self::XyPlots => "line",
            Self::GoalPlots => "graph",
            Self::Report => "list",
            Self::Animation => "ws_animate",
            Self::Compare => "compare",
            Self::Verify => "inspect",
        }
    }

    pub(super) fn stage(self) -> Stage {
        match self {
            Self::Project => Stage::Prepare,
            Self::Mesh => Stage::Mesh,
            Self::Parts => Stage::Physics,
            Self::Coupling => Stage::Couple,
            Self::Verify => Stage::Verify,
            node if node.edits_case() => Stage::Physics,
            _ => Stage::Results,
        }
    }

    fn from_stage(stage: Stage) -> Self {
        match stage {
            Stage::Prepare | Stage::Run => Self::Project,
            Stage::Physics => Self::Input,
            Stage::Mesh => Self::Mesh,
            Stage::Couple => Self::Coupling,
            Stage::Results => Self::Results,
            Stage::Verify => Self::Verify,
        }
    }

    /// Shown in a new project's tree before Customize Tree changes anything (the flow-analysis
    /// defaults); other folders appear once the case has entries for them.
    fn default_visible(self) -> bool {
        matches!(self, Self::Domain | Self::ComponentControl | Self::FluidSubdomains | Self::Boundary | Self::Goals | Self::Mesh)
            || (Self::RESULTS.contains(&self) && self.available())
    }

    /// Whether the tree lists this folder now: Customize Tree's choice, or the default, and
    /// always when the case has entries for it.
    /// A study of the loaded model (`model`) turns heat transfer on in Solid Materials, so that
    /// folder is always listed for it.
    fn visible(self, state: &State, case: &Value, model: bool) -> bool {
        // The Parts folder lists the design's bodies a study of the model holds; a case has none.
        if self == Self::Parts {
            return model && state.tree_custom.get(self.key()).copied().unwrap_or(true);
        }
        state.tree_custom.get(self.key()).copied().unwrap_or_else(|| self.default_visible())
            || case_has(case, self)
            || (model && self == Self::SolidMaterials)
    }

    /// Whether the node has an implementation; the others are listed for orientation only.
    /// Every tree node has a working panel now; kept so menus can still grey out items.
    pub(super) fn available(self) -> bool {
        let _ = self;
        true
    }

    /// Shown in the graphics window rather than in a floating window over it.
    fn in_graphics(self) -> bool {
        matches!(self, Self::Results | Self::CutPlots | Self::Animation | Self::Verify)
    }

    /// The case editor's view of this node: which top-level case entries it lists.
    pub(super) fn edits_case(self) -> bool {
        matches!(self, Self::Input | Self::General | Self::Control | Self::Mesh)
            || (Self::INPUT.contains(&self) && !matches!(self, Self::Coupling | Self::Parts))
    }

    /// The stage handle this node is also published under (scripts that drove the old stage
    /// bar keep working).
    fn stage_handle(self) -> Option<&'static str> {
        match self {
            Self::Project => Some("multiphysics.stage_prepare"),
            Self::Input => Some("multiphysics.stage_physics"),
            Self::Mesh => Some("multiphysics.stage_mesh"),
            Self::Coupling => Some("multiphysics.stage_couple"),
            Self::Results => Some("multiphysics.stage_results"),
            Self::Verify => Some("multiphysics.stage_verify"),
            _ => None,
        }
    }
}

/// The node shown now: the one picked in the tree while its stage is still current, otherwise
/// the stage's own node (code that switches stages directly lands on the matching node).
pub(super) fn current(state: &State) -> Node {
    match state.node {
        Some((stage, node)) if stage == state.stage => node,
        _ => Node::from_stage(state.stage),
    }
}

pub(super) fn select(state: &mut State, node: Node) {
    state.stage = node.stage();
    state.node = Some((state.stage, node));
    state.show_catalog = false;
    state.show_materials = false;
    state.dialog_open = true;
    state.entry = None;
}

/// The tree folder a top-level case entry belongs to. Case schemas differ per solver, so
/// entries are sorted by what their names say; anything unrecognised is a general setting.
/// Short tokens match whole `_`-separated words so `width` is not read as `dt`.
pub(super) fn group_of(key: &str) -> Node {
    let key = key.to_ascii_lowercase();
    let words: Vec<&str> = key.split('_').collect();
    let has = |parts: &[&str]| parts.iter().any(|part| key.contains(part));
    let word = |tokens: &[&str]| tokens.iter().any(|token| words.contains(token));
    if key == "solver" || key == "schema" {
        Node::General
    } else if key == "mesh"
        || key.ends_with("cells")
        || word(&["nx", "ny", "nz", "dx", "dy", "dz", "elements", "nodes", "segments", "bins"])
        || key == "grid_points"
        || has(&["mesh", "resolution", "refine"])
    {
        Node::Mesh
    } else if has(&["goal", "probe", "sensor", "monitor"]) {
        Node::Goals
    } else if has(&["injection"]) {
        Node::Sources
    } else if word(&[
        "step",
        "steps",
        "dt",
        "iteration",
        "iterations",
        "duration",
        "output",
        "snapshots",
        "frames",
        "interval",
        "tolerance",
        "continuation",
        "sweep",
        "increments",
        "points",
    ]) {
        Node::Control
    } else if word(&["battery", "soc", "cutoff", "protocol"]) {
        Node::Batteries
    } else if has(&["heat_pipe"]) {
        Node::HeatPipes
    } else if has(&["heat_sink"]) {
        Node::HeatSinks
    } else if has(&["thermal_joint"]) {
        Node::ThermalJoints
    } else if has(&["radiation_source"]) {
        Node::RadiationSources
    } else if word(&["tec", "thermoelectric", "peltier"]) {
        Node::Tec
    } else if word(&["led", "leds", "phosphor"]) {
        Node::Leds
    } else if word(&["pcb", "board", "via", "vias"]) {
        Node::Pcb
    } else if has(&["two_resistor", "junction"]) {
        Node::TwoResistor
    } else if has(&["transferred", "zoom"]) {
        Node::Transferred
    } else if has(&["electr", "voltage", "netlist", "impedance", "potential"])
        || word(&["current", "bus", "buses", "line", "lines", "resistance", "ohm", "charge", "donor", "acceptor", "carrier", "iv", "cv"])
    {
        Node::Electrical
    } else if has(&["initial"]) {
        Node::Initial
    } else if has(&["rotat", "swirl", "angular", "rpm", "stir"]) {
        Node::Rotating
    } else if has(&["fan"]) {
        Node::Fans
    } else if has(&["porous"]) {
        Node::Porous
    } else if has(&["perforat"]) {
        Node::Perforated
    } else if has(&["contact", "interface", "joint"]) {
        Node::Contact
    } else if has(&["radiat", "emissiv", "solar", "absorption", "scatter"]) {
        Node::Radiation
    } else if has(&["bound", "inlet", "outlet", "wall", "ambient", "environment", "support", "flow_rate"])
        || word(&["left", "right", "top", "bottom", "inflow", "outflow"])
    {
        Node::Boundary
    } else if has(&[
        "material",
        "solid",
        "conductivity",
        "density",
        "modulus",
        "poisson",
        "permittivity",
        "yield",
        "hardening",
        "capacity",
        "diffusivity",
        "sound_speed",
    ]) {
        Node::SolidMaterials
    } else if has(&["source", "heat", "power", "load", "excitation", "reaction", "beam", "flux", "wavelength", "frequenc"]) {
        Node::Sources
    } else if has(&["fluid", "subdomain", "phase", "species", "concentration", "dispersed", "gas", "liquid", "viscosity"]) {
        Node::FluidSubdomains
    } else if has(&[
        "length",
        "width",
        "height",
        "size",
        "radius",
        "depth",
        "box",
        "domain",
        "obstacle",
        "bounds",
        "geometry",
        "layer",
        "core",
        "cladding",
        "container",
        "nozzle",
        "pitch",
    ]) {
        Node::Domain
    } else {
        Node::General
    }
}

/// One segment of an RFC 6901 JSON pointer.
fn escape(segment: &str) -> String {
    segment.replace('~', "~0").replace('/', "~1")
}

/// Every editable case entry as `(JSON pointer, folder)`: the top-level entries, and the
/// children of general-purpose blocks such as `config` or `thermal`, which mix settings of
/// several folders (fluid, solids, walls, time step) and so are filed child by child.
pub(super) fn case_entries(case: &Value) -> Vec<(String, Node)> {
    let mut out = Vec::new();
    let Some(object) = case.as_object() else { return out };
    for (key, value) in object.iter().take(256) {
        let group = group_of(key);
        match value.as_object() {
            Some(children) if group == Node::General && !children.is_empty() => {
                for child in children.keys().take(256) {
                    out.push((format!("/{}/{}", escape(key), escape(child)), group_of(child)));
                }
            }
            _ => out.push((format!("/{}", escape(key)), group)),
        }
    }
    out
}

/// The entries filed under `node`, in case order.
pub(super) fn entries(case: &Value, node: Node) -> Vec<String> {
    case_entries(case).into_iter().filter(|(_, group)| *group == node).map(|(pointer, _)| pointer).collect()
}

/// Whether `case` has an entry for `node`'s folder (empty folders stay out of the tree).
fn case_has(case: &Value, node: Node) -> bool {
    case_entries(case).iter().any(|(_, group)| *group == node)
}

/// A readable name for the entry at `pointer`: its last segment with spaces for underscores.
pub(super) fn entry_label(pointer: &str) -> String {
    pointer.rsplit('/').next().unwrap_or(pointer).replace("~1", "/").replace("~0", "~").replace('_', " ")
}

/// The named items of a list entry (one row per fan, solid, patch or goal), as `(pointer,
/// label)`: an element's own `name`, `id`, `patch` or `type`, otherwise the list's name and
/// a number, as flow-analysis trees number their features.
fn items(case: &Value, pointer: &str) -> Vec<(String, String)> {
    let Some(list) = case.pointer(pointer).and_then(Value::as_array) else { return Vec::new() };
    if list.len() > 64 || !list.iter().all(Value::is_object) {
        return Vec::new();
    }
    let base = entry_label(pointer);
    list.iter()
        .enumerate()
        .map(|(index, item)| {
            let own = ["name", "id", "patch", "type"].iter().find_map(|field| item.get(*field).and_then(Value::as_str));
            (format!("{pointer}/{index}"), own.map_or_else(|| format!("{base} {}", index + 1), str::to_string))
        })
        .collect()
}

/// Seven-page study preparation wizard; edits remain transient until Finish.
#[derive(Clone, Debug, Default)]
pub(super) struct Wizard {
    step: usize,
    name: String,
    category: Option<String>,
    template: Option<usize>,
    candidate: Option<Value>,
    case_template: Option<String>,
    error: Option<String>,
}

impl Wizard {
    /// The wizard on its summary page for a template picked in the template list.
    pub(super) fn from_template(index: usize, category: String) -> Self {
        Self { step: WIZARD_STEPS.len() - 1, category: Some(category), template: Some(index), ..Default::default() }
    }
}

static WIZARD_STEPS: [&str; 7] = [
    "multiphysics.wizard_step_name",
    "multiphysics.wizard_setup_analysis",
    "multiphysics.wizard_setup_units",
    "multiphysics.wizard_setup_general",
    "multiphysics.wizard_setup_materials",
    "multiphysics.wizard_setup_initial",
    "multiphysics.wizard_step_summary",
];

/// The unit system the wizard states: every case is solved in SI. Symbols are not translated.
static SI_UNITS: [(&str, &str); 8] = [
    ("multiphysics.wizard_unit_length", "m"),
    ("multiphysics.wizard_unit_time", "s"),
    ("multiphysics.wizard_unit_mass", "kg"),
    ("multiphysics.wizard_unit_temperature", "K"),
    ("multiphysics.wizard_unit_pressure", "Pa"),
    ("multiphysics.wizard_unit_velocity", "m/s"),
    ("multiphysics.wizard_unit_heat", "W"),
    ("multiphysics.wizard_unit_conductivity", "W/(m·K)"),
];

/// Every group the summary page lists, in the pages' order.
const SUMMARY_GROUPS: &[Node] = &[
    Node::General,
    Node::Control,
    Node::Domain,
    Node::Mesh,
    Node::FluidSubdomains,
    Node::SolidMaterials,
    Node::Boundary,
    Node::Radiation,
    Node::Sources,
    Node::Fans,
    Node::Porous,
    Node::Initial,
];

/// A list row's click target, published only while the row is in view: a row scrolled out of
/// its list has no place on screen to click.
fn visible_target(ui: &egui::Ui, handle: &str, response: &egui::Response) {
    if ui.clip_rect().contains_rect(response.rect) {
        crate::scenario::publish_target(handle, response.rect);
    }
}

/// One step in the wizard's step list: a numbered disc (a tick once passed) and its name.
fn wizard_step_row(ui: &mut egui::Ui, index: usize, text: &str, current: usize, enabled: bool) -> egui::Response {
    let t = Tokens::get();
    let height = crate::ui_kit::row_height() + t.space_4;
    let (rect, response) = ui.allocate_exact_size(vec2(ui.available_width(), height), if enabled { Sense::click() } else { Sense::hover() });
    let painter = ui.painter();
    let here = index == current;
    if here {
        painter.rect_filled(rect, t.control_radius, t.accent_soft);
    } else if response.hovered() && enabled {
        painter.rect_filled(rect, t.control_radius, t.control_fill);
    }
    let radius = (height * 0.32).min(t.space_10);
    let center = pos2(rect.left() + t.space_8 + radius, rect.center().y);
    if index < current {
        painter.circle_filled(center, radius, t.accent);
        let r = radius * 0.5;
        painter
            .line(vec![center + vec2(-r, 0.0), center + vec2(-r * 0.25, r * 0.7), center + vec2(r, -r * 0.6)], egui::Stroke::new(1.8, t.on_accent));
    } else {
        if here {
            painter.circle_filled(center, radius, t.accent);
        } else {
            painter.circle_stroke(center, radius, egui::Stroke::new(1.0, t.border));
        }
        let ink = if here { t.on_accent } else { t.text_dim };
        painter.text(center, egui::Align2::CENTER_CENTER, (index + 1).to_string(), egui::FontId::proportional(t.font_small), ink);
    }
    let font = if here { crate::ui_kit::type_face::semibold(t.font_body) } else { egui::FontId::proportional(t.font_body) };
    let ink = if here {
        t.selected_ink
    } else if index < current {
        t.text
    } else {
        t.text_dim
    };
    let x = center.x + radius + t.space_8;
    let galley = painter.layout(text.to_owned(), font, ink, (rect.right() - x - t.space_4).max(0.0));
    let y = rect.center().y - galley.size().y / 2.0;
    painter.with_clip_rect(rect).galley(pos2(x, y), galley, ink);
    response.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, enabled, text));
    response
}

// ── Ribbon ───────────────────────────────────────────────────────────────────────────────────

struct Command {
    id: &'static str,
    key: &'static str,
    icon: &'static str,
    enabled: bool,
    pressed: bool,
}

fn cmd(id: &'static str, key: &'static str, icon: &'static str, enabled: bool) -> Command {
    Command { id, key, icon, enabled, pressed: false }
}

/// The ribbon's groups for tab `tab` (0: analysis, 1: results).
fn groups(app: &SolveApp, state: &State, study: Option<&Study>, tab: usize) -> Vec<(&'static str, &'static str, Vec<Command>)> {
    let has_study = study.is_some();
    let has_runs = study.is_some_and(|study| !study.runs.is_empty());
    // A finished run waiting to be kept: Load results keeps and opens it.
    let waiting = study.is_some_and(|study| keepable_job(app, study).is_some());
    let node = current(state);
    let running = study.is_some_and(|study| {
        app.session.multiphysics_jobs().iter().any(|job| job.study_id == study.id && matches!(job.state.as_str(), "queued" | "running"))
    });
    // One ribbon tab laid out like the reference flow-analysis tab: Project | Insert | Mesh |
    // Solve | Results | Tools | Help. `tab` is kept for callers; there is a single tab.
    let _ = tab;
    let mut out = vec![
        (
            "project",
            "multiphysics.group_project",
            vec![
                cmd("wizard", "multiphysics.cmd_wizard", "new", true),
                cmd("create", "multiphysics.cmd_new_project", "plus", true),
                cmd("model_study", "multiphysics.cmd_model_study", "box", true),
                cmd("clone", "multiphysics.menu_clone", "copy", has_study),
                cmd("general", "multiphysics.cmd_general_settings", "settings", has_study),
                cmd("component_control", "multiphysics.tree_component_control", "component", has_study),
                cmd("domain", "multiphysics.cmd_domain", "box", has_study),
            ],
        ),
        (
            "input",
            "multiphysics.menu_insert",
            vec![
                cmd("add_from_components", "multiphysics.cmd_add_components", "import", has_study),
                cmd("boundary", "multiphysics.cmd_boundary", "surface", has_study),
                cmd("goals", "multiphysics.cmd_goals", "checkpoint", has_study),
            ],
        ),
        (
            "mesh",
            "multiphysics.group_mesh",
            vec![
                cmd("mesh", "multiphysics.tree_global_mesh", "mesh_body", has_study),
                cmd("local_mesh", "multiphysics.menu_local_mesh", "mesh_body", has_study),
            ],
        ),
        (
            "solve",
            "multiphysics.group_solve",
            vec![
                cmd("run_solver", "multiphysics.cmd_run", "play", has_study && !running),
                cmd("batch_run", "multiphysics.cmd_batch_run", "timeline_play", has_study),
                cmd("control", "multiphysics.cmd_calculation_control", "settings_monitor", has_study),
                cmd("monitor", "multiphysics.cmd_monitor", "settings_monitor", has_study),
                cmd("stop_solver", "multiphysics.cmd_stop", "stop", running),
            ],
        ),
        (
            "results",
            "multiphysics.group_results",
            vec![
                cmd("load_results", "multiphysics.cmd_load_results", "open", has_runs || waiting),
                cmd("cut_plot", "multiphysics.tree_cut_plots", "section", has_runs),
                cmd("surface_plot", "multiphysics.tree_surface_plots", "surface", has_runs),
                cmd("isosurface", "multiphysics.tree_isosurfaces", "iso_analysis", has_runs),
                cmd("trajectories", "multiphysics.tree_trajectories", "spline", has_runs),
                cmd("particles", "multiphysics.tree_particles", "dots", has_runs),
                cmd("point_params", "multiphysics.tree_point_params", "point", has_runs),
                cmd("surface_params", "multiphysics.tree_surface_params", "plane", has_runs),
                cmd("xy_plot", "multiphysics.tree_xy_plots", "line", has_runs),
                cmd("goal_plot", "multiphysics.tree_goal_plots", "graph", has_runs),
                cmd("report", "multiphysics.tree_report", "list", has_runs),
                cmd("result_export", "multiphysics.cmd_export_results", "export", has_runs),
            ],
        ),
        (
            "tools",
            "multiphysics.group_tools",
            vec![
                cmd("check", "multiphysics.cmd_check", "render_check", has_study),
                cmd("materials", "multiphysics.cmd_engineering_database", "material", true),
                cmd("create_lids", "multiphysics.menu_create_lids", "surface", has_study),
                cmd("templates", "multiphysics.cmd_templates", "library", true),
            ],
        ),
        ("help", "multiphysics.menu_help", vec![cmd("verify", "multiphysics.menu_validation", "inspect", has_study)]),
    ];
    for (_, _, commands) in &mut out {
        for command in commands.iter_mut() {
            command.pressed = match command.id {
                "materials" => state.show_materials,
                "templates" => state.show_catalog,
                "monitor" => state.show_monitor || state.stage == Stage::Run,
                id => node_for(id).is_some_and(|target| target == node),
            };
        }
    }
    out
}

/// The latest finished run of `study` that has its result and is not kept yet, for the study's
/// current revision.
pub(super) fn keepable_job(app: &SolveApp, study: &Study) -> Option<String> {
    app.session
        .multiphysics_jobs()
        .iter()
        .rev()
        .find(|job| {
            job.study_id == study.id
                && job.evidence.is_some()
                && !job.stale
                && job.study_revision == study.revision
                && !study.runs.iter().any(|run| run.id == job.id)
        })
        .map(|job| job.id.clone())
}

/// The tree node a ribbon command opens.
fn node_for(id: &str) -> Option<Node> {
    Some(match id {
        "general" => Node::General,
        "control" => Node::Control,
        "domain" => Node::Domain,
        "boundary" => Node::Boundary,
        "fluid_subdomains" => Node::FluidSubdomains,
        "solid_materials" => Node::SolidMaterials,
        "fans" => Node::Fans,
        "radiation" => Node::Radiation,
        "radiation_sources" => Node::RadiationSources,
        "heat_sinks" => Node::HeatSinks,
        "thermal_joint" => Node::ThermalJoints,
        "transferred" => Node::Transferred,
        "rotating" => Node::Rotating,
        "fluid_volume" => Node::Domain,
        "contact" => Node::Contact,
        "create_lids" => Node::Boundary,
        "local_mesh" => Node::Mesh,
        "pcb" => Node::Pcb,
        "leds" => Node::Leds,
        "batteries" => Node::Batteries,
        "perforated" => Node::Perforated,
        "tec" => Node::Tec,
        "two_resistor" => Node::TwoResistor,
        "heat_pipes" => Node::HeatPipes,
        "component_control" => Node::ComponentControl,
        "add_from_components" => Node::FromComponents,
        "porous" => Node::Porous,
        "initial" => Node::Initial,
        "electrical" => Node::Electrical,
        "sources" => Node::Sources,
        "goals" => Node::Goals,
        "coupling" => Node::Coupling,
        "mesh" => Node::Mesh,
        "check" | "case_export" => Node::Project,
        "verify" => Node::Verify,
        "load_results" => Node::Results,
        "cut_plot" => Node::CutPlots,
        "point_params" => Node::PointParams,
        "surface_plot" => Node::SurfacePlots,
        "isosurface" => Node::Isosurfaces,
        "trajectories" => Node::Trajectories,
        "particles" => Node::Particles,
        "surface_params" => Node::SurfaceParams,
        "volume_params" => Node::VolumeParams,
        "xy_plot" => Node::XyPlots,
        "goal_plot" => Node::GoalPlots,
        "report" | "result_export" => Node::Report,
        "animation" => Node::Animation,
        "compare" => Node::Compare,
        _ => return None,
    })
}

/// Commands that stay in their group's menu (and the menu bar) instead of the ribbon row: one
/// row of every command squeezed each button to its narrowest and cut most labels short.
const MENU_ONLY: &[&str] = &[
    "clone",
    "component_control",
    "add_from_components",
    "batch_run",
    "control",
    "isosurface",
    "particles",
    "point_params",
    "surface_params",
    "create_lids",
    "templates",
];

/// A button wide enough for its label's longest word, which cannot wrap: the kit's width caps a
/// word such as "Boundary" and cuts it short.
fn label_room(painter: &egui::Painter, label: &str) -> f32 {
    let word = label.split_whitespace().map(|w| ribbon::label_width(painter, w)).fold(0.0_f32, f32::max);
    ribbon::button_width(painter, label).max(word + ribbon::BUTTON_LABEL_PADDING)
}

/// A ribbon label short enough for the button (the full name stays its accessible name and the
/// group menu's row).
fn ribbon_label(id: &str) -> Option<&'static str> {
    Some(match id {
        "model_study" => "multiphysics.ribbon_model_study",
        "general" => "multiphysics.ribbon_general",
        "domain" => "multiphysics.ribbon_domain",
        "boundary" => "multiphysics.ribbon_boundary",
        "mesh" => "multiphysics.ribbon_global_mesh",
        "local_mesh" => "multiphysics.ribbon_local_mesh",
        "monitor" => "multiphysics.ribbon_monitor",
        "load_results" => "multiphysics.ribbon_load_results",
        "trajectories" => "multiphysics.ribbon_trajectories",
        "xy_plot" => "multiphysics.ribbon_xy_plot",
        "goal_plot" => "multiphysics.ribbon_goal_plot",
        "result_export" => "multiphysics.ribbon_export",
        "check" => "multiphysics.ribbon_check",
        "materials" => "multiphysics.ribbon_materials",
        "verify" => "multiphysics.ribbon_verify",
        _ => return None,
    })
}

const RIBBON_HANDLES: ribbon::Handles<'static> = ribbon::Handles {
    button: &["multiphysics:"],
    button_rect: "multiphysics_button:",
    split: "multiphysics_split:",
    cell: &["multiphysics:"],
    collapsed: "multiphysics_collapsed:",
    section: "multiphysics_group:",
    section_target: "multiphysics_group_target:",
    group: "multiphysics_group:",
    running: "multiphysics_running:",
    reveal: "multiphysics_reveal",
    more: "multiphysics:more",
    more_row: "multiphysics:",
    collapsed_count: "multiphysics_collapsed_count",
    more_count: "multiphysics_more_count",
    trimmed_count: "multiphysics_trimmed_count",
    drag: "multiphysics_drag",
    drag_remove: "multiphysics_drag:remove",
};

const TAB_HANDLES: ribbon::TabHandles<'static> =
    ribbon::TabHandles { tab: "multiphysics_tab:", tab_rect: "multiphysics_tab_rect:", add: "multiphysics_tab:add", action: "multiphysics:" };

static TABS: [(&str, &str); 1] = [("FLOW", "multiphysics.ribbon_tab_flow")];

/// The workspace switcher, the tab row (with Save, Undo and Redo at its right, as a quick
/// access toolbar) and the active tab's groups.
pub(super) fn ribbon_bar(app: &mut SolveApp, state: &mut State, studies: &[Study], ui: &mut egui::Ui) {
    let t = Tokens::get();
    let study = state.selected.as_ref().and_then(|id| studies.iter().find(|study| &study.id == id));
    let whole = ui.max_rect();
    let switcher = Rect::from_min_size(whole.min + vec2(t.space_4, t.space_4), vec2(SWITCHER - 2.0 * t.space_4, whole.height() - 2.0 * t.space_4));
    let response =
        crate::ui_kit::button_at(ui, "multiphysics_workspace_switch", switcher, tr("multiphysics.workspace"), crate::ui_kit::Button::Normal);
    crate::scenario::publish_target("ribbon:workspace", response.rect);
    if response.clicked() {
        crate::context_menu::toggle_for(app, response.rect, response.rect.left_bottom(), crate::context_menu::Target::Workspace);
    }
    let left = whole.left() + SWITCHER;
    let menu_row = Rect::from_min_max(pos2(left, whole.top()), pos2(whole.right(), whole.top() + menus::HEIGHT));
    ui.scope_builder(egui::UiBuilder::new().max_rect(menu_row), |ui| menus::show(app, state, studies, ui));
    let tab_row = Rect::from_min_max(pos2(left, menu_row.bottom()), pos2(whole.right(), menu_row.bottom() + ribbon::TAB_ROW_HEIGHT));
    let labels: Vec<&str> = TABS.iter().map(|(_, key)| tr(key)).collect();
    let quick = [("save", tr("multiphysics.save")), ("undo", tr("multiphysics.undo")), ("redo", tr("multiphysics.redo"))];
    let accent = crate::workspace::tab_accent(&t, "MODEL");
    let active = TABS.get(state.ribbon_tab).map_or("FLOW", |&(id, _)| id);
    let out = {
        let row = ribbon::Tabs {
            tabs: TABS.iter().zip(&labels).map(|(&(id, _), &label)| ribbon::Tab { id, label, accent, contextual: false, aliases: &[] }).collect(),
            active,
            add: None,
            menu_open_from: &|_| false,
            search_width: 0.0,
            buttons: quick.iter().map(|&(id, label)| ribbon::TabButton { id, label, ready: Ok(()) }).collect(),
            floating: false,
            underline: None,
            handles: TAB_HANDLES,
        };
        ribbon::tabs(ui, &row, tab_row)
    };
    for action in out.actions {
        match action {
            ribbon::TabAction::Select { index, .. } => state.ribbon_tab = index.min(TABS.len() - 1),
            ribbon::TabAction::Run { id } => match id.as_str() {
                "save" => app.start("file.save"),
                "undo" => {
                    command(app, state, "edit.undo", json!({}));
                }
                "redo" => {
                    command(app, state, "edit.redo", json!({}));
                }
                _ => {}
            },
            ribbon::TabAction::Menu { .. } | ribbon::TabAction::Add { .. } => {}
        }
    }
    let groups = groups(app, state, study, state.ribbon_tab);
    let painter = ui.painter().clone();
    let texts: Vec<Vec<&str>> = groups.iter().map(|(_, _, commands)| commands.iter().map(|command| tr(command.key)).collect()).collect();
    // The ribbon row's buttons: the group's commands without the menu-only ones.
    let shown: Vec<Vec<usize>> =
        groups.iter().map(|(_, _, commands)| (0..commands.len()).filter(|&i| !MENU_ONLY.contains(&commands[i].id)).collect()).collect();
    let group_labels: Vec<&str> = groups.iter().map(|(_, key, _)| tr(key)).collect();
    let ids: Vec<Vec<&str>> = groups.iter().map(|(_, _, commands)| commands.iter().map(|command| command.id).collect()).collect();
    let more_tip = tr("multiphysics.ribbon_more");
    let row = ribbon::Row {
        groups: groups
            .iter()
            .enumerate()
            .map(|(gi, (key, _, commands))| ribbon::Group {
                key,
                name: group_labels.get(gi).copied().unwrap_or_default(),
                label: group_labels.get(gi).copied().unwrap_or_default(),
                tip: group_labels.get(gi).copied().unwrap_or_default(),
                buttons: shown
                    .get(gi)
                    .into_iter()
                    .flatten()
                    .filter_map(|&i| Some((commands.get(i)?, *texts.get(gi)?.get(i)?)))
                    .map(|(command, name)| (command, ribbon_label(command.id).map_or(name, tr), name))
                    .map(|(command, label, name)| ribbon::Button {
                        id: command.id,
                        run: command.id,
                        label,
                        name,
                        icon: command.icon,
                        size: ribbon::ButtonSize { width: label_room(&painter, label), min_width: ribbon::BUTTON_MIN_WIDTH },
                        enabled: command.enabled,
                        planned: false,
                        pressed: command.pressed,
                        exit: false,
                        split: None,
                        app_paints: false,
                        covers: &[],
                    })
                    .collect(),
                cells: Vec::new(),
                commands: ids.get(gi).cloned().unwrap_or_default(),
                icon: commands.first().map(|command| command.icon),
                enabled: false,
                open: false,
                running: None,
                more: commands
                    .iter()
                    .zip(texts.get(gi).into_iter().flatten())
                    .map(|(command, &label)| ribbon::MoreRow {
                        id: command.id,
                        label,
                        icon: command.icon,
                        shortcut: "",
                        enabled: command.enabled,
                        checked: None,
                        hint: None,
                    })
                    .collect(),
            })
            .collect(),
        labels: true,
        accent,
        exit: false,
        reveal: None,
        more_tip,
        handles: RIBBON_HANDLES,
    };
    let panel_rect = Rect::from_min_max(pos2(left, tab_row.bottom()), whole.max);
    let drawn = ribbon::show(ui, &row, panel_rect, &mut |_, _| false);
    let mut chosen = Vec::new();
    for button in &drawn.buttons {
        if let Some(command) = shown
            .get(button.group)
            .and_then(|row| row.get(button.index))
            .and_then(|&i| groups.get(button.group).and_then(|(_, _, commands)| commands.get(i)))
            && command.id == "monitor"
        {
            crate::scenario::publish_target("multiphysics.stage_run", button.whole);
        }
    }
    for action in drawn.actions {
        if let ribbon::Action::Run { id, .. } = action {
            chosen.push(id);
        }
    }
    for id in chosen {
        run(app, state, studies, &id);
    }
}

/// What a ribbon command does.
fn run(app: &mut SolveApp, state: &mut State, studies: &[Study], id: &str) {
    if id == "fluid_volume" {
        state.case_tool.fluid_requested = true;
    }
    let study = state.selected.as_ref().and_then(|selected| studies.iter().find(|study| &study.id == selected)).cloned();
    // Load results keeps a finished run first, so its result can be shown at once.
    if id == "load_results"
        && let Some(study) = &study
        && let Some(job) = keepable_job(app, study)
        && command(app, state, "multiphysics.run.keep", json!({"job_id":job})).is_some()
    {
        state.selected_run = Some(job.clone());
        state.shown_run = Some(job);
        state.selected_artifact = 0;
    }
    if let Some(node) = node_for(id) {
        if node.stage() == Stage::Results && id != "load_results" {
            state.ribbon_tab = 1;
        }
        select(state, node);
        return;
    }
    match id {
        "wizard" => {
            state.wizard = Some(Wizard::default());
        }
        "create" => create_from_template(app, state, None),
        "model_study" => {
            let name = solidcraft_ui_text::format("multiphysics.model_study_name", &[("design", &app.session.doc.name)]);
            if let Some(reply) = command(app, state, "multiphysics.model.create", json!({"name":name})) {
                state.selected = reply["study"]["id"].as_str().map(str::to_string);
                state.selected_run = None;
                state.model_ui.last_build = state.selected.clone().map(|id| (id, reply));
                select(state, Node::Project);
            }
        }
        "clone" => {
            if let Some(study) = study
                && let Some(reply) = command(app, state, "multiphysics.study.clone", json!({"id":study.id}))
            {
                state.selected = reply["study"]["id"].as_str().map(str::to_string);
                state.selected_run = None;
                select(state, Node::Project);
            }
        }
        "unload_results" => {
            select(state, Node::Project);
            state.dialog_open = false;
        }
        "customize_tree" => state.show_customize = true,
        "remove" => {
            if let Some(study) = study {
                command(app, state, "multiphysics.study.remove", json!({"id":study.id}));
            }
        }
        "materials" => {
            state.show_materials = !state.show_materials;
            state.show_catalog = false;
        }
        "templates" => state.show_catalog = !state.show_catalog,
        "monitor" => state.show_monitor = !(state.show_monitor || state.stage == Stage::Run),
        "run_solver" => {
            if let Some(study) = study {
                start_run(app, state, &study);
            }
        }
        "batch_run" => batch_run(app, state, studies),
        "stop_solver" => {
            if let Some(study) = study {
                let jobs: Vec<String> = app
                    .session
                    .multiphysics_jobs()
                    .iter()
                    .filter(|job| job.study_id == study.id && matches!(job.state.as_str(), "queued" | "running"))
                    .map(|job| job.id.clone())
                    .collect();
                for job in jobs {
                    command(app, state, "multiphysics.run.cancel", json!({"job_id":job}));
                }
            }
        }
        _ => {}
    }
}

/// Run the committed case when it can run; the solver monitor opens either way, so the reason it
/// did not start is in front of the person.
pub(super) fn start_run(app: &mut SolveApp, state: &mut State, study: &Study) {
    // A study of the loaded model runs on the model as it is now.
    if !super::model_ui::ensure_current(app, state, study) {
        state.show_monitor = true;
        return;
    }
    let study = &app.session.doc.multiphysics.study(&study.id).cloned().unwrap_or_else(|| study.clone());
    let report = state.study_preflight(&app.session.doc.multiphysics.clone(), study);
    let ready =
        Readiness::for_study(study, state.drafts.get(&study.id), app.session.multiphysics_host_available(), report.errors.is_empty()).can_run();
    state.show_monitor = true;
    if !ready {
        state.error = Some(
            tr(if app.session.multiphysics_host_available() { "multiphysics.apply_before_run" } else { "multiphysics.preparation_only" }).to_string(),
        );
        return;
    }
    command(app, state, "multiphysics.run.start", json!({"id":study.id}));
}

/// Batch Run: start every project that is ready and not already queued or running. Projects
/// that are not ready (unapplied edits, failing checks, no solver host) are skipped and counted.
fn batch_run(app: &mut SolveApp, state: &mut State, studies: &[Study]) {
    let host = app.session.multiphysics_host_available();
    let mut started = 0usize;
    for study in studies {
        let busy = app.session.multiphysics_jobs().iter().any(|job| job.study_id == study.id && matches!(job.state.as_str(), "queued" | "running"));
        let report = state.study_preflight(&app.session.doc.multiphysics.clone(), study);
        let ready = Readiness::for_study(study, state.drafts.get(&study.id), host, report.errors.is_empty()).can_run();
        if !busy && ready && command(app, state, "multiphysics.run.start", json!({"id":study.id})).is_some() {
            started += 1;
        }
    }
    state.show_monitor = started > 0;
    state.notice = Some(format!("{} {started}/{}", tr("multiphysics.batch_run_started"), studies.len()));
}

/// A new project from template `index` (the wizard's choice), or from the template last used.
fn create_from_template(app: &mut SolveApp, state: &mut State, chosen: Option<(usize, String)>) {
    let templates = state.templates.as_ref().and_then(|result| result.as_ref().ok()).cloned().unwrap_or_default();
    let index = chosen.as_ref().map_or(state.template_index, |(index, _)| *index);
    let Some(template) = templates.get(index).or_else(|| templates.first()) else { return };
    let name = chosen.map(|(_, name)| name).filter(|name| !name.trim().is_empty()).unwrap_or_else(|| template.name.clone());
    if let Some(reply) = command(app, state, "multiphysics.study.create", json!({"template":template.id,"name":name})) {
        state.template_index = index;
        state.selected = reply.get("study").and_then(|study| study.get("id")).and_then(Value::as_str).map(str::to_string);
        state.selected_run = None;
        state.selected_artifact = 0;
        select(state, Node::Project);
    }
}

// ── Resource bar and Analysis tree ───────────────────────────────────────────────────────────

/// The resource bar's icon strip: Analysis tree, engineering database, templates.
pub(super) fn resource_bar(state: &mut State, ui: &mut egui::Ui) {
    let t = Tokens::get();
    let mut top = ui.max_rect().top() + t.space_4;
    let side = RESOURCE_BAR - 2.0 * t.space_4;
    for (icon, key, target) in [
        ("list", "multiphysics.resource_tree", "multiphysics:resource_tree"),
        ("material", "multiphysics.cmd_engineering_database", "multiphysics:resource_materials"),
        ("library", "multiphysics.cmd_templates", "multiphysics:resource_templates"),
    ] {
        let rect = Rect::from_min_size(pos2(ui.max_rect().left() + t.space_4, top), vec2(side, side));
        let on = match target {
            "multiphysics:resource_tree" => state.show_browser && !state.show_materials,
            "multiphysics:resource_materials" => state.show_materials,
            _ => state.show_catalog,
        };
        if on {
            ui.painter().rect_filled(rect.expand(t.space_2), t.control_radius, t.accent_soft);
        }
        let response = crate::ui_kit::icon_button_at(ui, rect, icon, tr(key));
        crate::scenario::publish_target(target, response.rect);
        if response.clicked() {
            match target {
                "multiphysics:resource_tree" => {
                    state.show_browser = !state.show_browser || state.show_materials;
                    state.show_materials = false;
                }
                "multiphysics:resource_materials" => {
                    state.show_materials = !state.show_materials;
                    state.show_catalog = false;
                }
                _ => state.show_catalog = !state.show_catalog,
            }
        }
        top += side + t.space_6;
    }
}

/// What a tree row asks for.
pub(super) enum TreeAction {
    Select(String, Node),
    /// Open the Goals panel with the given goal kind picked for Add.
    InsertGoal(String, &'static str),
    /// Toggle the basic mesh (grid lines) in the graphics window.
    ShowBasicMesh(String),
    Entry(String, Node, String),
    /// Study, folder, run and the artifact's index in that run.
    Artifact(String, Node, String, usize),
    Fold(&'static str),
    Run(String),
    Clone(String),
    Remove(String),
    LoadRun(String, String),
    Unload,
    Hide(Node),
    Customize,
    Revert(String, String),
    Copy(String),
    /// Duplicate a list item in the draft, right after itself.
    CloneItem(String, String),
    /// Remove a list item from the draft.
    DeleteItem(String, String),
    /// Copy an entry's draft value into another study's draft at the same pointer.
    CopyTo(String, String, Value),
    /// Pick, hide, isolate or ghost parts (the Parts folder's rows).
    Part(super::parts_ui::PartAction),
}

/// A context-menu entry (kit ghost button) that closes the menu when chosen.
fn menu_item(ui: &mut egui::Ui, key: &str) -> bool {
    let chosen = crate::ui_kit::button(ui, tr(key), crate::ui_kit::Button::Ghost).clicked();
    if chosen {
        ui.close();
    }
    chosen
}

/// A context-menu entry listed for orientation but not available, with the reason on hover.
fn menu_item_off(ui: &mut egui::Ui, key: &str, why: &str) {
    let response = ui.add_enabled_ui(false, |ui| crate::ui_kit::button(ui, tr(key), crate::ui_kit::Button::Ghost)).inner;
    let _ = response.backing_on_hover_text(tr(why));
}

/// The right-click menu of an Input data folder. Goals and Mesh carry the reference tree's own
/// commands; other folders insert a feature of their type through the same dialog as the
/// Insert menu.
fn folder_menu(ui: &mut egui::Ui, study: &Study, folder: Node, actions: &mut Vec<TreeAction>) {
    match folder {
        Node::Goals => {
            for (key, kind) in [
                ("multiphysics.tree_insert_global_goals", "global_goal"),
                ("multiphysics.tree_insert_point_goals", "point_goal"),
                ("multiphysics.tree_insert_surface_goals", "goal"),
                ("multiphysics.tree_insert_volume_goals", "volume_goal"),
                ("multiphysics.tree_insert_equation_goal", "equation_goal"),
            ] {
                if menu_item(ui, key) {
                    actions.push(TreeAction::InsertGoal(study.id.clone(), kind));
                }
            }
        }
        Node::Mesh => {
            if menu_item(ui, "multiphysics.tree_global_mesh") {
                actions.push(TreeAction::Select(study.id.clone(), Node::Mesh));
            }
            if menu_item(ui, "multiphysics.tree_insert_local_mesh") {
                actions.push(TreeAction::Select(study.id.clone(), Node::Mesh));
            }
            if menu_item(ui, "multiphysics.tree_show_basic_mesh") {
                actions.push(TreeAction::ShowBasicMesh(study.id.clone()));
            }
            ui.separator();
            menu_item_off(ui, "multiphysics.tree_create_mesh", "multiphysics.tree_mesh_at_run");
            ui.separator();
            if menu_item(ui, "multiphysics.tree_delete_local_meshes") {
                actions.push(TreeAction::Select(study.id.clone(), Node::Mesh));
            }
            return;
        }
        _ => {
            if menu_item(ui, "multiphysics.tree_insert") {
                actions.push(TreeAction::Select(study.id.clone(), folder));
            }
        }
    }
    ui.separator();
    if menu_item(ui, "multiphysics.tree_hide_folder") {
        actions.push(TreeAction::Hide(folder));
    }
    if menu_item(ui, "multiphysics.tree_customize") {
        actions.push(TreeAction::Customize);
    }
}

/// The list and index of a list item's pointer (`/list/3` → `/list`, 3).
fn list_item(pointer: &str) -> Option<(&str, usize)> {
    let (parent, index) = pointer.rsplit_once('/')?;
    Some((parent, index.parse().ok()?))
}

/// Parse a study's draft case, change it and write it back; a draft that does not parse is left
/// alone with a message (fix the JSON first).
fn edit_draft(state: &mut State, id: &str, change: impl FnOnce(&mut Value) -> bool) {
    let Some(draft) = state.drafts.get_mut(id) else { return };
    let Ok(mut case) = serde_json::from_str::<Value>(&draft.case_text) else {
        state.error = Some(tr("multiphysics.tree_revert_invalid").to_string());
        return;
    };
    if change(&mut case) {
        draft.case_text = pretty(&case);
    } else {
        state.error = Some(tr("multiphysics.tree_revert_missing").to_string());
    }
}

/// The Analysis tree: every project of the design, the active one expanded into General
/// settings, Calculation control, the Input data folders with one row per case entry (and per
/// named item of a list), the Results folders with one row per result file, and Verification.
/// Right-click any row for its commands.
pub(super) fn tree(app: &mut SolveApp, state: &mut State, studies: &[Study], ui: &mut egui::Ui) {
    let t = Tokens::get();
    ui.add_space(t.space_4);
    ui.horizontal(|ui| {
        ui.add_space(t.space_6);
        ui.strong(tr("multiphysics.analysis_tree"));
    });
    ui.add_space(t.space_4);
    // The active project's structural check, for its error state in the tree.
    let errors = state
        .selected
        .as_ref()
        .and_then(|id| studies.iter().find(|study| &study.id == id))
        .map_or(0, |study| state.study_preflight(&app.session.doc.multiphysics.clone(), study).errors.len());
    let mut actions = Vec::new();
    egui::ScrollArea::vertical().id_salt("multiphysics_analysis_tree").show(ui, |ui| {
        if studies.is_empty() {
            ui.add_space(t.space_8);
            ui.label(RichText::new(tr("multiphysics.tree_empty")).small().color(t.text_dim));
        }
        let node = current(state);
        for study in studies {
            let active = state.selected.as_ref() == Some(&study.id);
            let dirty = state.drafts.get(&study.id).is_some_and(Draft::dirty);
            // The draft's case when it parses, so new entries show before Apply.
            let case = state
                .drafts
                .get(&study.id)
                .and_then(|draft| serde_json::from_str::<Value>(&draft.case_text).ok())
                .unwrap_or_else(|| study.case.clone());
            let title = if dirty { format!("{} •", study.name) } else { study.name.clone() };
            let icon = if active && errors > 0 { "status_error" } else { Node::Project.icon() };
            let response = row(ui, 0, Some(active), icon, &title, active && node == Node::Project, true, active);
            crate::scenario::publish_target(&format!("multiphysics:study:{}", study.id), response.rect);
            if active && let Some(handle) = Node::Project.stage_handle() {
                crate::scenario::publish_target(handle, response.rect);
            }
            let response = if active && errors > 0 { response.backing_on_hover_text(tr("multiphysics.tree_project_errors")) } else { response };
            if response.clicked() {
                actions.push(TreeAction::Select(study.id.clone(), Node::Project));
            }
            crate::ui_kit::Menu::context(&response, |ui| {
                if menu_item(ui, "multiphysics.tree_activate") {
                    actions.push(TreeAction::Select(study.id.clone(), Node::Project));
                }
                if menu_item(ui, "multiphysics.cmd_run") {
                    actions.push(TreeAction::Run(study.id.clone()));
                }
                if menu_item(ui, "multiphysics.tree_clone_project") {
                    actions.push(TreeAction::Clone(study.id.clone()));
                }
                ui.separator();
                if menu_item(ui, "multiphysics.cmd_general_settings") {
                    actions.push(TreeAction::Select(study.id.clone(), Node::General));
                }
                if menu_item(ui, "multiphysics.cmd_calculation_control") {
                    actions.push(TreeAction::Select(study.id.clone(), Node::Control));
                }
                if !study.runs.is_empty() && menu_item(ui, "multiphysics.cmd_load_results") {
                    actions.push(TreeAction::Select(study.id.clone(), Node::Results));
                }
                ui.separator();
                if menu_item(ui, "multiphysics.tree_customize") {
                    actions.push(TreeAction::Customize);
                }
                if menu_item(ui, "multiphysics.cmd_delete_project") {
                    actions.push(TreeAction::Remove(study.id.clone()));
                }
            });
            if !active {
                continue;
            }
            // Input data.
            let open = !state.folders_closed[0];
            let label = if dirty { format!("{} •", tr(Node::Input.key())) } else { tr(Node::Input.key()).to_string() };
            let response = row(ui, 1, Some(open), Node::Input.icon(), &label, node == Node::Input, true, false);
            publish_node(Node::Input, response.rect);
            if disclosure_clicked(ui, &response, 1) || response.double_clicked() {
                actions.push(TreeAction::Fold("input"));
            } else if response.clicked() {
                actions.push(TreeAction::Select(study.id.clone(), Node::Input));
            }
            crate::ui_kit::Menu::context(&response, |ui| {
                if menu_item(ui, "multiphysics.tree_edit_definition") {
                    actions.push(TreeAction::Select(study.id.clone(), Node::Input));
                }
                if menu_item(ui, "multiphysics.tree_customize") {
                    actions.push(TreeAction::Customize);
                }
            });
            if open {
                for &folder in Node::INPUT {
                    if !folder.visible(state, &case, study.model.is_some()) {
                        continue;
                    }
                    // A study of the loaded model rebuilds its case from the model, so the
                    // case's own entries are not listed for editing.
                    let pointers = if study.model.is_some() { Vec::new() } else { entries(&case, folder) };
                    let folded = state.tree_folded.contains(folder.key());
                    let parts = if folder == Node::Parts {
                        // Parts coloured by their setup stay so while the design changes.
                        super::parts_ui::keep_colour(app, state, study);
                        super::parts_ui::bodies(app, state, study)
                    } else {
                        Vec::new()
                    };
                    let has_rows = !pointers.is_empty() || folder == Node::Mesh || !parts.is_empty();
                    let response = row(
                        ui,
                        2,
                        has_rows.then_some(!folded),
                        folder.icon(),
                        tr(folder.key()),
                        node == folder && state.entry.is_none(),
                        true,
                        false,
                    );
                    publish_node(folder, response.rect);
                    if has_rows && disclosure_clicked(ui, &response, 2) {
                        actions.push(TreeAction::Fold(folder.key()));
                    } else if response.clicked() {
                        actions.push(TreeAction::Select(study.id.clone(), folder));
                    }
                    crate::ui_kit::Menu::context(&response, |ui| folder_menu(ui, study, folder, &mut actions));
                    if folded || folder == Node::Coupling {
                        continue;
                    }
                    if folder == Node::Parts {
                        // One row per part, picked when the design's selection holds it.
                        let chosen = super::parts_ui::selected(app);
                        for body in &parts {
                            let on = chosen.iter().any(|b| b == body);
                            let label = super::parts_ui::label(app, body);
                            let response = row(ui, 3, None, "body", &label, on, true, false);
                            crate::scenario::publish_target(&format!("multiphysics:tree_part:{body}"), response.rect);
                            if response.clicked() {
                                actions.push(TreeAction::Part(super::parts_ui::pick(ui, body)));
                            }
                            let mut picked = Vec::new();
                            crate::ui_kit::Menu::context(&response, |ui| super::parts_ui::menu(app, ui, body, &mut picked));
                            actions.extend(picked.into_iter().map(TreeAction::Part));
                        }
                        continue;
                    }
                    if folder == Node::Mesh {
                        // The project-wide mesh settings, as the reference tree's Global Mesh.
                        let picked = node == Node::Mesh && state.entry.is_none();
                        let response = row(ui, 3, None, Node::Mesh.icon(), tr("multiphysics.tree_global_mesh"), picked, true, false);
                        crate::scenario::publish_target("multiphysics:tree:global_mesh", response.rect);
                        if response.clicked() {
                            actions.push(TreeAction::Select(study.id.clone(), Node::Mesh));
                        }
                        crate::ui_kit::Menu::context(&response, |ui| {
                            if menu_item(ui, "multiphysics.tree_edit_definition") {
                                actions.push(TreeAction::Select(study.id.clone(), Node::Mesh));
                            }
                            if menu_item(ui, "multiphysics.tree_show_basic_mesh") {
                                actions.push(TreeAction::ShowBasicMesh(study.id.clone()));
                            }
                        });
                    }
                    for pointer in pointers {
                        let picked = node == folder && state.entry.as_deref() == Some(pointer.as_str());
                        let response = row(ui, 3, None, "params", &entry_label(&pointer), picked, true, false);
                        crate::scenario::publish_target(&format!("multiphysics:tree_entry:{pointer}"), response.rect);
                        if response.clicked() {
                            actions.push(TreeAction::Entry(study.id.clone(), folder, pointer.clone()));
                        }
                        entry_menu(&response, study, studies, folder, &pointer, &case, &mut actions);
                        for (item, label) in items(&case, &pointer) {
                            let picked = node == folder && state.entry.as_deref() == Some(item.as_str());
                            let response = row(ui, 4, None, "component", &label, picked, true, false);
                            crate::scenario::publish_target(&format!("multiphysics:tree_entry:{item}"), response.rect);
                            if response.clicked() {
                                actions.push(TreeAction::Entry(study.id.clone(), folder, item.clone()));
                            }
                            entry_menu(&response, study, studies, folder, &item, &case, &mut actions);
                        }
                    }
                }
            }
            // Results.
            let open = !state.folders_closed[1];
            let loaded = state.selected_run.as_deref().and_then(|id| study.runs.iter().find(|run| run.id == id)).or_else(|| study.runs.last());
            let stale = loaded.is_some_and(|run| run.study_revision != study.revision);
            let label = match loaded {
                Some(run) => format!("{} ({})", tr(Node::Results.key()), run.id),
                None => tr("multiphysics.tree_results_not_loaded").to_string(),
            };
            let icon = if stale { "status_warning" } else { Node::Results.icon() };
            let response = row(ui, 1, Some(open), icon, &label, node == Node::Results, loaded.is_some(), false);
            publish_node(Node::Results, response.rect);
            let response = if stale { response.backing_on_hover_text(tr("multiphysics.stale_run")) } else { response };
            if disclosure_clicked(ui, &response, 1) || response.double_clicked() {
                actions.push(TreeAction::Fold("results"));
            } else if response.clicked() && loaded.is_some() {
                actions.push(TreeAction::Select(study.id.clone(), Node::Results));
            }
            crate::ui_kit::Menu::context(&response, |ui| {
                for run in study.runs.iter().rev().take(16) {
                    let text = solidcraft_ui_text::format("multiphysics.tree_load_run", &[("run", &run.id)]);
                    if crate::ui_kit::button(ui, text, crate::ui_kit::Button::Ghost).clicked() {
                        actions.push(TreeAction::LoadRun(study.id.clone(), run.id.clone()));
                        ui.close();
                    }
                }
                if loaded.is_some() {
                    ui.separator();
                    if menu_item(ui, "multiphysics.tree_unload_results") {
                        actions.push(TreeAction::Unload);
                    }
                    if menu_item(ui, "multiphysics.tree_summary") {
                        actions.push(TreeAction::Select(study.id.clone(), Node::Report));
                    }
                    if menu_item(ui, "multiphysics.tree_compare") {
                        actions.push(TreeAction::Select(study.id.clone(), Node::Compare));
                    }
                }
                if menu_item(ui, "multiphysics.tree_customize") {
                    actions.push(TreeAction::Customize);
                }
            });
            if open && let Some(run) = loaded {
                for &folder in Node::RESULTS {
                    if !folder.visible(state, &Value::Null, false) {
                        continue;
                    }
                    let files: Vec<(usize, &str)> = run
                        .artifacts
                        .iter()
                        .enumerate()
                        .filter(|(_, artifact)| match folder {
                            Node::CutPlots => artifact.name.ends_with(".vtu"),
                            Node::XyPlots => artifact.name.ends_with(".csv"),
                            Node::Animation => artifact.name.ends_with(".pvd"),
                            _ => false,
                        })
                        .map(|(index, artifact)| (index, artifact.name.as_str()))
                        .take(64)
                        .collect();
                    let folded = state.tree_folded.contains(folder.key());
                    let response = row(
                        ui,
                        2,
                        (!files.is_empty()).then_some(!folded),
                        folder.icon(),
                        tr(folder.key()),
                        node == folder,
                        folder.available(),
                        false,
                    );
                    publish_node(folder, response.rect);
                    if !folder.available() {
                        let _ = response.clone().backing_on_hover_text(tr("multiphysics.tree_unavailable"));
                    } else if !files.is_empty() && disclosure_clicked(ui, &response, 2) {
                        actions.push(TreeAction::Fold(folder.key()));
                    } else if response.clicked() {
                        actions.push(TreeAction::Select(study.id.clone(), opens(folder)));
                    }
                    crate::ui_kit::Menu::context(&response, |ui| {
                        if folder.available() && menu_item(ui, "multiphysics.tree_insert") {
                            actions.push(TreeAction::Select(study.id.clone(), opens(folder)));
                        }
                        if menu_item(ui, "multiphysics.tree_hide_folder") {
                            actions.push(TreeAction::Hide(folder));
                        }
                    });
                    if folded {
                        continue;
                    }
                    for (index, name) in files {
                        let picked = node == folder && state.selected_artifact == index;
                        let response = row(ui, 3, None, folder.icon(), name, picked, true, false);
                        crate::scenario::publish_target(&format!("multiphysics:tree_result:{name}"), response.rect);
                        if response.clicked() {
                            actions.push(TreeAction::Artifact(study.id.clone(), folder, run.id.clone(), index));
                        }
                        crate::ui_kit::Menu::context(&response, |ui| {
                            if menu_item(ui, "multiphysics.tree_show") {
                                actions.push(TreeAction::Artifact(study.id.clone(), folder, run.id.clone(), index));
                            }
                            if menu_item(ui, "multiphysics.cmd_export_results") {
                                actions.push(TreeAction::Select(study.id.clone(), Node::Report));
                            }
                        });
                    }
                }
            }
            let response = row(ui, 1, None, Node::Verify.icon(), tr(Node::Verify.key()), node == Node::Verify, true, false);
            publish_node(Node::Verify, response.rect);
            if response.clicked() {
                actions.push(TreeAction::Select(study.id.clone(), Node::Verify));
            }
        }
    });
    for action in actions {
        apply_tree_action(app, state, studies, action, ui.ctx());
    }
}

/// The node a Results row opens: Export Results is the report's export section.
fn opens(folder: Node) -> Node {
    if folder == Node::ExportResults { Node::Report } else { folder }
}

/// The right-click menu of a case entry or list item, in the reference tree's order: Edit
/// Definition, Clone, Delete, Copy to Project, then SolidCraft's Revert and Copy as JSON. Clone
/// and Delete act on list items only; removing a whole case entry could leave the case without
/// a setting its solver requires.
fn entry_menu(response: &egui::Response, study: &Study, studies: &[Study], folder: Node, pointer: &str, case: &Value, actions: &mut Vec<TreeAction>) {
    crate::ui_kit::Menu::context(response, |ui| {
        if menu_item(ui, "multiphysics.tree_edit_definition") {
            actions.push(TreeAction::Entry(study.id.clone(), folder, pointer.to_string()));
        }
        ui.separator();
        if list_item(pointer).is_some() {
            if menu_item(ui, "multiphysics.tree_clone") {
                actions.push(TreeAction::CloneItem(study.id.clone(), pointer.to_string()));
            }
            if menu_item(ui, "multiphysics.tree_delete") {
                actions.push(TreeAction::DeleteItem(study.id.clone(), pointer.to_string()));
            }
        } else {
            menu_item_off(ui, "multiphysics.tree_clone", "multiphysics.tree_items_only");
            menu_item_off(ui, "multiphysics.tree_delete", "multiphysics.tree_items_only");
        }
        ui.separator();
        let others: Vec<&Study> = studies.iter().filter(|other| other.id != study.id).take(16).collect();
        if others.is_empty() {
            menu_item_off(ui, "multiphysics.tree_copy_to_project", "multiphysics.tree_no_other_project");
        } else if let Some(value) = case.pointer(pointer) {
            let _ = crate::ui_kit::Menu::button(ui, tr("multiphysics.tree_copy_to_project"), |ui| {
                for other in others {
                    if crate::ui_kit::button(ui, &other.name, crate::ui_kit::Button::Ghost).clicked() {
                        actions.push(TreeAction::CopyTo(other.id.clone(), pointer.to_string(), value.clone()));
                        ui.close();
                    }
                }
            });
        }
        ui.separator();
        if menu_item(ui, "multiphysics.tree_revert") {
            actions.push(TreeAction::Revert(study.id.clone(), pointer.to_string()));
        }
        if menu_item(ui, "multiphysics.tree_copy_json") {
            actions.push(TreeAction::Copy(case.pointer(pointer).map(pretty).unwrap_or_default()));
        }
    });
}

pub(super) fn apply_tree_action(app: &mut SolveApp, state: &mut State, studies: &[Study], action: TreeAction, ctx: &egui::Context) {
    let activate = |state: &mut State, id: String| {
        if state.selected.as_ref() != Some(&id) {
            state.selected_run = None;
            state.selected_artifact = 0;
        }
        state.selected = Some(id);
    };
    match action {
        TreeAction::ShowBasicMesh(id) => {
            activate(state, id);
            select(state, Node::Mesh);
            state.case_view.mesh = !state.case_view.mesh;
            state.case_view.cad = false;
        }
        TreeAction::InsertGoal(id, kind) => {
            activate(state, id);
            select(state, Node::Goals);
            state.goal_kind = kind;
        }
        TreeAction::Select(id, node) => {
            activate(state, id);
            if node.stage() == Stage::Results && node != Node::Results {
                state.ribbon_tab = 1;
            }
            select(state, node);
        }
        TreeAction::Entry(id, node, pointer) => {
            activate(state, id);
            select(state, node);
            state.entry = Some(pointer);
        }
        TreeAction::Artifact(id, node, run, index) => {
            activate(state, id);
            state.ribbon_tab = 1;
            select(state, node);
            // The row named this run's file: show that run as already shown, or the Results
            // view would reset the file to its first when it first draws the run.
            state.selected_run = Some(run.clone());
            state.shown_run = Some(run);
            state.selected_artifact = index;
        }
        TreeAction::Fold(key) => match key {
            "input" => state.folders_closed[0] = !state.folders_closed[0],
            "results" => state.folders_closed[1] = !state.folders_closed[1],
            key => {
                if !state.tree_folded.remove(key) {
                    state.tree_folded.insert(key);
                }
            }
        },
        TreeAction::Run(id) => {
            if let Some(study) = studies.iter().find(|study| study.id == id) {
                state.selected = Some(id);
                start_run(app, state, study);
            }
        }
        TreeAction::Clone(id) => {
            if let Some(reply) = command(app, state, "multiphysics.study.clone", json!({"id":id})) {
                state.selected = reply.get("study").and_then(|study| study.get("id")).and_then(Value::as_str).map(str::to_string);
                state.selected_run = None;
                select(state, Node::Project);
            }
        }
        TreeAction::Remove(id) => {
            command(app, state, "multiphysics.study.remove", json!({"id":id}));
        }
        TreeAction::LoadRun(id, run) => {
            activate(state, id);
            state.selected_run = Some(run);
            state.selected_artifact = 0;
            select(state, Node::Results);
        }
        TreeAction::Unload => {
            select(state, Node::Project);
            state.dialog_open = false;
        }
        TreeAction::Hide(node) => {
            state.tree_custom.insert(node.key(), false);
        }
        TreeAction::Customize => state.show_customize = true,
        TreeAction::Part(action) => super::parts_ui::apply(app, state, action),
        TreeAction::Revert(id, pointer) => {
            let Some(draft) = state.drafts.get_mut(&id) else { return };
            let Ok(mut case) = serde_json::from_str::<Value>(&draft.case_text) else {
                state.error = Some(tr("multiphysics.tree_revert_invalid").to_string());
                return;
            };
            match (case.pointer_mut(&pointer), draft.base_case.pointer(&pointer)) {
                (Some(slot), Some(applied)) => {
                    *slot = applied.clone();
                    draft.case_text = pretty(&case);
                }
                _ => state.error = Some(tr("multiphysics.tree_revert_missing").to_string()),
            }
        }
        TreeAction::Copy(text) => ctx.copy_text(text),
        TreeAction::CloneItem(id, pointer) => edit_draft(state, &id, |case| {
            let Some((parent, index)) = list_item(&pointer) else { return false };
            let Some(list) = case.pointer_mut(parent).and_then(Value::as_array_mut) else { return false };
            let Some(item) = list.get(index).cloned() else { return false };
            if list.len() >= 64 {
                return false;
            }
            list.insert(index + 1, item);
            true
        }),
        TreeAction::DeleteItem(id, pointer) => {
            edit_draft(state, &id, |case| {
                let Some((parent, index)) = list_item(&pointer) else { return false };
                let Some(list) = case.pointer_mut(parent).and_then(Value::as_array_mut) else { return false };
                if index >= list.len() {
                    return false;
                }
                list.remove(index);
                true
            });
            state.entry = None;
        }
        TreeAction::CopyTo(id, pointer, value) => edit_draft(state, &id, |case| match case.pointer_mut(&pointer) {
            Some(slot) => {
                *slot = value;
                true
            }
            None => false,
        }),
    }
}

/// The Customize Tree window: which Input data and Results folders the tree lists. A folder
/// the case has entries for is always listed, whatever is chosen here.
pub(super) fn customize(state: &mut State, ctx: &egui::Context) {
    if !state.show_customize {
        return;
    }
    let closed = crate::ui_kit::Window::new(egui::Id::unique("multiphysics_customize_tree"), "list", tr("multiphysics.tree_customize"))
        .cards(1.6)
        .page_rows(22.0)
        .show_body(ctx, |ui| {
            crate::scenario::publish_rect("multiphysics:customize_tree", ui.max_rect());
            ui.label(RichText::new(tr("multiphysics.tree_customize_hint")).small().color(Tokens::get().text_dim));
            for (heading, folders) in [("multiphysics.tree_input", Node::INPUT), ("multiphysics.tree_results", Node::RESULTS)] {
                ui.add_space(Tokens::get().space_6);
                ui.strong(tr(heading));
                for &folder in folders {
                    let mut shown = state.tree_custom.get(folder.key()).copied().unwrap_or_else(|| folder.default_visible());
                    let label = if folder.available() {
                        tr(folder.key()).to_string()
                    } else {
                        format!("{} — {}", tr(folder.key()), tr("multiphysics.tree_not_built"))
                    };
                    if crate::ui_kit::checkbox(ui, &label, &mut shown).changed() {
                        state.tree_custom.insert(folder.key(), shown);
                    }
                }
            }
            ui.add_space(Tokens::get().space_6);
            if button(ui, "multiphysics.tree_customize_defaults", "multiphysics:tree_customize_defaults", !state.tree_custom.is_empty()).clicked() {
                state.tree_custom.clear();
            }
        });
    if closed {
        state.show_customize = false;
    }
}

fn publish_node(node: Node, rect: Rect) {
    crate::scenario::publish_target(&format!("multiphysics:tree:{}", node.key().rsplit('.').next().unwrap_or_default()), rect);
    if let Some(handle) = node.stage_handle() {
        crate::scenario::publish_target(handle, rect);
    }
}

/// Whether the click landed on the disclosure triangle of a row at `depth`.
fn disclosure_clicked(ui: &egui::Ui, response: &egui::Response, depth: usize) -> bool {
    response.clicked()
        && ui.ctx().pointer_interact_pos().is_some_and(|at| at.x < response.rect.left() + disclosure_left(depth) + Tokens::get().space_16)
}

use crate::multiphysics_property_grid::{disclosure_left, tree_row as row};

// ── Floating windows ─────────────────────────────────────────────────────────────────────────

/// Feature dialogs and result windows over the graphics window, the solver monitor, the
/// template list and the wizard.
pub(super) fn windows(
    app: &mut SolveApp,
    state: &mut State,
    study: Option<&Study>,
    _report: Option<&solidcraft_engine::multiphysics::Validation>,
    graphics: Rect,
    ctx: &egui::Context,
) {
    let t = Tokens::get();
    let at = graphics.left_top() + vec2(t.space_12, t.space_12);
    if let Some(study) = study
        && (state.show_monitor || state.stage == Stage::Run)
    {
        let title = format!("{} — {}", tr("multiphysics.cmd_monitor"), study.name);
        let closed = crate::ui_kit::Window::new(egui::Id::unique("multiphysics_solver_monitor"), "settings_monitor", &title)
            .placed(egui::Align2::RIGHT_TOP, graphics.right_top() + vec2(-t.space_12, t.space_12))
            .cards(1.6)
            .show_body(ctx, |ui| {
                crate::scenario::publish_rect("multiphysics:solver_monitor", ui.max_rect());
                jobs(app, state, study, ui);
            });
        if closed {
            state.show_monitor = false;
            if state.stage == Stage::Run {
                state.stage = Stage::Prepare;
            }
        }
    }
    if state.show_catalog {
        let closed = crate::ui_kit::Window::new(egui::Id::unique("multiphysics_templates"), "library", tr("multiphysics.cmd_templates"))
            .placed(egui::Align2::LEFT_TOP, at)
            .cards(2.2)
            .page_rows(24.0)
            .show_body(ctx, |ui| catalog(state, ui));
        if closed {
            state.show_catalog = false;
        }
    }
    if state.wizard.is_some() {
        wizard(app, state, ctx);
    }
    customize(state, ctx);
}

/// Whether the selected node's feature dialog is shown in the right-hand property panel.
pub(super) fn feature_docked(state: &State) -> bool {
    state.dialog_open && !current(state).in_graphics() && !state.show_materials
}

/// The right-hand property panel: the feature's title, its dialog, and the ✓ / ✕ footer of a
/// flow-analysis property panel. ✓ applies pending case edits (when this node edits the case)
/// and closes; ✕ closes and keeps the draft, so nothing typed is lost.
pub(super) fn feature_panel(
    app: &mut SolveApp,
    state: &mut State,
    study: &Study,
    report: Option<&solidcraft_engine::multiphysics::Validation>,
    ui: &mut egui::Ui,
) {
    let t = Tokens::get();
    let node = current(state);
    crate::scenario::publish_rect("multiphysics:feature_window", ui.max_rect());
    let title = match &state.entry {
        Some(pointer) if node.edits_case() => format!("{} — {}", tr(node.key()), entry_label(pointer)),
        _ => tr(node.key()).to_string(),
    };
    ui.add_space(t.space_4);
    ui.horizontal(|ui| {
        let (rect, _) = ui.allocate_exact_size(vec2(t.space_16, t.space_16), Sense::hover());
        crate::icons::paint_icon(ui.painter(), rect, node.icon(), 1.0);
        ui.strong(&title);
    });
    ui.label(RichText::new(&study.name).small().color(t.text_dim));
    ui.separator();
    let footer = t.control_height + t.space_12;
    egui::ScrollArea::vertical().id_salt("multiphysics_feature_panel").max_height((ui.available_height() - footer).max(80.0)).show(ui, |ui| {
        // Hints and notes wrap at the panel's width instead of widening it to their length.
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Wrap);
        ui.set_max_width(ui.available_width());
        feature(app, state, study, report, node, ui);
    });
    ui.separator();
    ui.horizontal(|ui| {
        let dirty = state.drafts.get(&study.id).is_some_and(Draft::dirty);
        let ok = crate::ui_kit::button(ui, tr("multiphysics.panel_ok_label"), crate::ui_kit::Button::Primary)
            .backing_on_hover_text(tr("multiphysics.panel_ok"));
        crate::scenario::publish_target("multiphysics:panel_ok", ok.rect);
        let close = crate::ui_kit::button(ui, tr("multiphysics.panel_close_label"), crate::ui_kit::Button::Normal)
            .backing_on_hover_text(tr("multiphysics.panel_close"));
        crate::scenario::publish_target("multiphysics:panel_close", close.rect);
        if ok.clicked() {
            if node.edits_case() && dirty {
                apply(app, state, study);
            }
            if state.error.is_none() {
                state.dialog_open = false;
            }
        }
        if close.clicked() {
            state.dialog_open = false;
        }
    });
}

/// The body of the property panel for `node`.
fn feature(
    app: &mut SolveApp,
    state: &mut State,
    study: &Study,
    report: Option<&solidcraft_engine::multiphysics::Validation>,
    node: Node,
    ui: &mut egui::Ui,
) {
    // A study of the loaded model is set up on the model itself: where the fluid is, the
    // grid and the conditions on its faces.
    // Its mesh and goals have their own nodes, as in the analysis tree of a flow tool.
    if study.model.is_some() && node == Node::Mesh {
        super::model_ui::mesh_panel(app, state, study, ui);
        return;
    }
    if study.model.is_some() && node == Node::SolidMaterials {
        super::model_ui::heat_panel(app, state, study, ui);
        return;
    }
    if study.model.is_some() && node == Node::Parts {
        super::parts_ui::panel(app, state, study, ui);
        return;
    }
    if study.model.is_some() && node == Node::Goals {
        super::model_ui::goals_panel(app, state, study, ui);
        return;
    }
    if study.model.is_some() && matches!(node, Node::Project | Node::General | Node::Domain | Node::Boundary) {
        super::model_ui::panel(app, state, study, ui);
        if node != Node::Project {
            return;
        }
        ui.separator();
    }
    match node {
        Node::Project => {
            inspector(app, state, Some(study), report, ui);
            if let Some(report) = report {
                ui.separator();
                prepare(app, state, study, report, ui);
            }
        }
        Node::Coupling => {
            ui.label(tr("multiphysics.coupling_unavailable"));
            ui.label(tr("multiphysics.coupling_preserved"));
        }
        node if node.edits_case() => {
            let entry = state.entry.clone();
            if let Some(pointer) = &entry {
                ui.strong(entry_label(pointer));
                ui.label(RichText::new(pointer).small().monospace().color(Tokens::get().text_dim));
            }
            settings(app, state, study, node, entry.as_deref(), ui);
            ui.separator();
            let dirty = state.drafts.get(&study.id).is_some_and(Draft::dirty);
            ui.horizontal(|ui| {
                if button(ui, "multiphysics.apply", "multiphysics:apply", dirty).clicked() {
                    apply(app, state, study);
                }
                if button(ui, "multiphysics.cancel_draft", "multiphysics:cancel_draft", dirty).clicked() {
                    state.drafts.insert(study.id.clone(), Draft::from_study(study));
                    state.error = None;
                }
            });
            if let Some(draft) = state.drafts.get(&study.id)
                && draft.expected_revision != study.revision
            {
                ui.colored_label(Tokens::get().warning, tr("multiphysics.conflict"));
            }
        }
        _ => {
            egui::ScrollArea::vertical().id_salt("multiphysics_result_window").show(ui, |ui| results(app, state, study, ui));
        }
    }
}

/// The graphics window: the design itself, or the plot that the selected result node draws.
pub(super) fn graphics(app: &mut SolveApp, state: &mut State, study: Option<&Study>, ui: &mut egui::Ui) {
    let t = Tokens::get();
    let rect = ui.max_rect();
    crate::scenario::publish_rect("multiphysics:content", rect);
    let node = current(state);
    if state.show_materials {
        crate::multiphysics_material_library::show(app, &mut state.materials, ui);
        return;
    }
    match study {
        Some(study) if node.in_graphics() => {
            ui.add_space(t.space_8);
            ui.horizontal(|ui| {
                ui.add_space(t.space_8);
                ui.heading(format!("{} — {}", tr(node.key()), study.name));
            });
            match node {
                Node::Verify => verify(study, ui),
                _ => {
                    egui::ScrollArea::vertical().id_salt("multiphysics_results_scroll").show(ui, |ui| results(app, state, study, ui));
                }
            }
        }
        _ => {
            let mut offer_case = false;
            // A 3-D airflow case is drawn as itself, with its domain, parts and face roles,
            // and picked in place; anything else (and "Show CAD model") is the CAD window.
            if let Some(study) = study
                && study.model.is_none()
                && node.stage() != Stage::Results
                && let Some(text) = state.drafts.get(&study.id).map(|draft| draft.case_text.clone())
                && crate::multiphysics_case_view::drawable(&mut state.case_view, &text)
            {
                state.case_view.own(&study.id);
                if !state.case_view.cad {
                    let mode = super::case_mode(node, state.case_tool.fan_internal);
                    state.case_view.force_mesh = node == Node::Mesh;
                    state.case_view.show_pockets = node == Node::Domain;
                    if let Some(outcome) = crate::multiphysics_case_view::show(app, ui, &mut state.case_view, &text, mode) {
                        if let Some((face, cells)) = outcome.move_face
                            && let Ok(case) = serde_json::from_str::<Value>(&text)
                        {
                            let cell = solidcraft_engine::multiphysics::case_scene::scene(&case).map(|scene| scene.edge_cell(face)).unwrap_or(0.0);
                            let result = solidcraft_engine::multiphysics::case_scene::move_domain_face(&case, face, cells);
                            if super::apply_case_edit(state, study, result) {
                                state.case_view.moved(face, cells, cell);
                            }
                        }
                        // Picking in the graphics window opens the matching panel.
                        if outcome.picked && !state.case_view.selection.is_empty() && node.edits_case() && !feature_docked(state) {
                            select(state, node);
                        }
                    }
                    return;
                }
                offer_case = true;
            }
            crate::viewport::show(app, ui);
            // A study of the loaded model draws its domain and conditions over the model.
            if let Some(study) = study
                && study.model.is_some()
            {
                super::model_ui::overlay(app, state, study, ui, rect);
                // A flow tool open in the property panel draws its lines on this model.
                #[cfg(not(target_arch = "wasm32"))]
                if state.view.model_flow_docked && matches!(current(state), Node::Trajectories | Node::Particles | Node::Isosurfaces) {
                    super::on_model::flow_on_graphics(app, ui, &mut state.view, rect);
                }
            }
            // With parts ghosted, the whole model in the corner, over what is drawn on it.
            crate::overview_inset::paint(app, ui, rect);
            // Drawn after the CAD window so it stays on top of it.
            if offer_case {
                let corner = Rect::from_min_size(rect.min + vec2(t.space_8, t.space_8), vec2(260.0, t.control_height + t.space_4));
                ui.scope_builder(egui::UiBuilder::new().max_rect(corner), |ui| {
                    let response = crate::ui_kit::checkbox(ui, tr("multiphysics.case_show_cad"), &mut state.case_view.cad);
                    crate::scenario::publish_target("multiphysics:case_show_cad", response.rect);
                });
            }
            // The cue line: what the graphics window is showing and what to do next.
            let cue = match study {
                None => tr("multiphysics.cue_no_project"),
                Some(study) if study.runs.is_empty() => tr("multiphysics.cue_setup"),
                Some(_) => tr("multiphysics.cue_results"),
            };
            let galley = ui.painter().layout(cue.to_owned(), egui::FontId::proportional(t.font_small), t.text_dim, rect.width() - 2.0 * t.space_12);
            let at = pos2(rect.left() + t.space_12, rect.bottom() - t.space_8 - galley.size().y);
            ui.painter().galley(at, galley, t.text_dim);
        }
    }
}

/// The New Project wizard, one page per step.
fn wizard(app: &mut SolveApp, state: &mut State, ctx: &egui::Context) {
    let Some(mut wizard) = state.wizard.take() else { return };
    let templates = state.templates.as_ref().and_then(|r| r.as_ref().ok()).cloned().unwrap_or_default();
    let mut categories: Vec<String> = templates.iter().map(|t| template_category(&t.name).to_string()).collect();
    categories.sort();
    categories.dedup();
    let synchronize = |wizard: &mut Wizard| {
        let chosen = wizard.template.and_then(|i| templates.get(i));
        let id = chosen.map(|t| t.id.as_str());
        if wizard.case_template.as_deref() != id {
            wizard.case_template = id.map(str::to_string);
            wizard.candidate = chosen.map(|t| t.case.clone());
            wizard.error = None;
        }
    };
    synchronize(&mut wizard);
    wizard.step = wizard.step.min(6);
    let title = solidcraft_ui_text::format(
        "multiphysics.wizard_title",
        &[("step", &(wizard.step + 1).to_string()), ("count", "7"), ("name", tr(WIZARD_STEPS[wizard.step]))],
    );
    let mut finish = false;
    let mut cancel = false;
    let closed = crate::ui_kit::Window::new(egui::Id::unique("multiphysics_wizard"), "new", &title).cards(2.8).page_rows(24.0).show_body(ctx, |ui| {
        crate::scenario::publish_rect("multiphysics:wizard", ui.max_rect());
        crate::scenario::publish_count("multiphysics:wizard_step", wizard.step as f64);
        let t = Tokens::get();
        // Fill the window's body down to the footer, whatever its rows make it.
        let footer = t.control_height + t.space_16 + t.space_8;
        let room = ui.available_height() - footer;
        let height = if room.is_finite() && room >= 360.0 { room } else { (t.control_height * 16.5).max(360.0) };
        let top = ui.available_rect_before_wrap().min;
        let width = ui.available_width();
        let rail_width = (width * 0.27).clamp(170.0, 230.0);
        let rail = Rect::from_min_size(top, vec2(rail_width, height));
        let page = Rect::from_min_max(pos2(rail.right() + t.space_12, top.y), pos2(top.x + width, top.y + height));
        ui.allocate_rect(Rect::from_min_size(top, vec2(width, height)), Sense::hover());
        let reachable = wizard.candidate.is_some();
        crate::multiphysics_property_grid::pane(ui, rail, tr("multiphysics.wizard_steps"), |ui| {
            for (index, key) in WIZARD_STEPS.iter().enumerate() {
                let enabled = index <= wizard.step || reachable;
                let response = wizard_step_row(ui, index, tr(key), wizard.step, enabled);
                crate::scenario::publish_target(&format!("multiphysics:wizard_rail:{index}"), response.rect);
                if response.clicked() && enabled {
                    wizard.step = index;
                }
            }
        });
        ui.scope_builder(egui::UiBuilder::new().max_rect(page).layout(egui::Layout::top_down(egui::Align::Min)), |ui| {
            ui.set_clip_rect(page.intersect(ui.clip_rect()));
            ui.label(egui::RichText::new(tr(WIZARD_STEPS[wizard.step])).font(crate::ui_kit::type_face::semibold(t.font_subtitle)));
            let hint = match wizard.step {
                0 => "multiphysics.wizard_name_hint",
                1 => "multiphysics.wizard_template_hint",
                2 => "multiphysics.wizard_si_scope",
                3..=5 => "multiphysics.wizard_case_scope",
                _ => "multiphysics.wizard_summary_hint",
            };
            ui.add(egui::Label::new(egui::RichText::new(tr(hint)).small().color(t.text_dim)).wrap());
            ui.add_space(t.space_8);
            let columns = [tr("multiphysics.wizard_col_parameter"), tr("multiphysics.wizard_col_value")];
            match wizard.step {
                0 => {
                    crate::multiphysics_property_grid::header(ui, &columns, &[0.46, 0.54]);
                    let (_, r) = crate::multiphysics_property_grid::row(ui, 0, tr("multiphysics.name"), false, |ui: &mut egui::Ui| {
                        crate::ui_kit::bounded_cell_editor(ui, ui.scope_id().with("wizard_name"), &mut wizard.name, "", false, 160)
                    });
                    crate::scenario::publish_target("multiphysics:wizard_name", r.rect);
                }
                1 => {
                    // Every category as a filter chip, then the chosen category's templates.
                    ui.horizontal_wrapped(|ui| {
                        let response = crate::ui_kit::choice_chip(
                            ui,
                            "multiphysics_wizard_category_all",
                            tr("multiphysics.template_category_all"),
                            wizard.category.is_none(),
                        );
                        crate::scenario::publish_target("multiphysics:wizard_category:all", response.rect);
                        if response.clicked() {
                            wizard.category = None;
                        }
                        for category in &categories {
                            let selected = wizard.category.as_ref() == Some(category);
                            let response = crate::ui_kit::choice_chip(ui, &format!("multiphysics_wizard_category_{category}"), category, selected);
                            crate::scenario::publish_target(&format!("multiphysics:wizard_category:{category}"), response.rect);
                            if response.clicked() {
                                wizard.category = Some(category.clone());
                            }
                        }
                    });
                    ui.add_space(t.space_4);
                    let right = ui.available_rect_before_wrap();
                    ui.allocate_rect(right, Sense::hover());
                    ui.scope_builder(egui::UiBuilder::new().max_rect(right).layout(egui::Layout::top_down(egui::Align::Min)), |ui| {
                        let shares = [0.56, 0.44];
                        crate::multiphysics_property_grid::header(
                            ui,
                            &[tr("multiphysics.wizard_step_template"), tr("multiphysics.wizard_profile")],
                            &shares,
                        );
                        egui::ScrollArea::vertical().id_salt("multiphysics_wizard_templates").auto_shrink([false, false]).show(ui, |ui| {
                            let mut shown = 0usize;
                            for (i, template) in templates.iter().enumerate() {
                                if wizard.category.as_deref().is_some_and(|c| template_category(&template.name) != c) {
                                    continue;
                                }
                                let response = crate::multiphysics_property_grid::item(
                                    ui,
                                    &[&template.name, &template.profile],
                                    &shares,
                                    wizard.template == Some(i),
                                    !shown.is_multiple_of(2),
                                );
                                shown += 1;
                                visible_target(ui, &format!("multiphysics:wizard_template:{}", template.id), &response);
                                if response.clicked() {
                                    wizard.template = Some(i);
                                }
                                if response.double_clicked() {
                                    wizard.template = Some(i);
                                    wizard.step += 1;
                                }
                            }
                        });
                    });
                    synchronize(&mut wizard);
                }
                2 => {
                    let shares = [0.46, 0.54];
                    crate::multiphysics_property_grid::header(
                        ui,
                        &[tr("multiphysics.wizard_col_parameter"), tr("multiphysics.wizard_col_unit")],
                        &shares,
                    );
                    for (index, (key, unit)) in SI_UNITS.iter().enumerate() {
                        crate::multiphysics_property_grid::item(ui, &[tr(key), unit], &shares, false, !index.is_multiple_of(2));
                    }
                }
                3..=5 => {
                    let groups: &[Node] = match wizard.step {
                        3 => &[Node::General, Node::Control, Node::Domain, Node::Mesh],
                        4 => &[Node::FluidSubdomains, Node::SolidMaterials, Node::Boundary, Node::Radiation, Node::Sources, Node::Fans, Node::Porous],
                        _ => &[Node::Initial],
                    };
                    crate::multiphysics_property_grid::header(ui, &columns, &[0.46, 0.54]);
                    let r = ui
                        .scope(|ui| {
                            egui::ScrollArea::vertical().id_salt(("wizard_case", wizard.step)).auto_shrink([false, true]).show(ui, |ui| {
                                let mut rows = 0;
                                if let Some(case) = &mut wizard.candidate {
                                    for &group in groups {
                                        super::case_table(ui, case, group, &mut rows);
                                    }
                                }
                                if rows == 0 {
                                    ui.add_space(t.space_8);
                                    ui.label(tr("multiphysics.wizard_no_case_fields"));
                                }
                            });
                        })
                        .response;
                    crate::scenario::publish_rect(
                        match wizard.step {
                            3 => "multiphysics:wizard_case_general",
                            4 => "multiphysics:wizard_case_materials",
                            _ => "multiphysics:wizard_case_initial",
                        },
                        r.rect,
                    );
                }
                _ => {
                    crate::multiphysics_property_grid::header(ui, &columns, &[0.46, 0.54]);
                    egui::ScrollArea::vertical().id_salt("multiphysics_wizard_summary").auto_shrink([false, true]).show(ui, |ui| {
                        if let Some(t) = wizard.template.and_then(|i| templates.get(i)) {
                            let name = if wizard.name.trim().is_empty() { t.name.clone() } else { wizard.name.clone() };
                            crate::multiphysics_property_grid::row(ui, 0, tr("multiphysics.name"), false, |ui: &mut egui::Ui| ui.label(name));
                            crate::multiphysics_property_grid::row(ui, 0, tr("multiphysics.wizard_step_template"), true, |ui: &mut egui::Ui| {
                                ui.label(&t.name)
                            });
                            crate::multiphysics_property_grid::row(ui, 0, tr("multiphysics.wizard_profile"), false, |ui: &mut egui::Ui| {
                                ui.label(&t.profile)
                            });
                        }
                        // The whole case, read only: edits belong to their own pages.
                        if let Some(case) = &wizard.candidate {
                            let mut copy = case.clone();
                            let mut rows = 0;
                            ui.add_enabled_ui(false, |ui| {
                                for &group in SUMMARY_GROUPS {
                                    super::case_table(ui, &mut copy, group, &mut rows);
                                }
                            });
                        }
                    });
                }
            }
        });
        let valid = wizard.candidate.as_ref().map(solidcraft_engine::multiphysics::validate_case);
        if let Some(Err(e)) = &valid {
            ui.colored_label(Tokens::get().form_error, e.to_string());
        }
        if let Some(e) = &wizard.error {
            ui.colored_label(Tokens::get().form_error, e);
        }
        ui.separator();
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if button(ui, "multiphysics.cancel_draft", "multiphysics:wizard_cancel", true).clicked() {
                cancel = true;
            }
            if wizard.step >= 6 && button(ui, "multiphysics.wizard_finish", "multiphysics:wizard_finish", matches!(valid, Some(Ok(())))).clicked() {
                finish = true;
            }
            if wizard.step < 6
                && button(ui, "multiphysics.wizard_next", "multiphysics:wizard_next", wizard.step == 0 || wizard.candidate.is_some()).clicked()
            {
                wizard.step += 1;
            }
            if button(ui, "multiphysics.wizard_back", "multiphysics:wizard_back", wizard.step > 0).clicked() {
                wizard.step = wizard.step.saturating_sub(1);
            }
        });
    });
    if closed || cancel {
        return;
    }
    if finish
        && let Some(t) = wizard.template.and_then(|i| templates.get(i))
        && let Some(case) = wizard.candidate.as_ref()
    {
        let name = if wizard.name.trim().is_empty() { &t.name } else { &wizard.name };
        if let Some(reply) = command(app, state, "multiphysics.study.create", json!({"template":t.id,"name":name,"case":case})) {
            state.selected = reply["study"]["id"].as_str().map(str::to_string);
            state.selected_run = None;
            state.selected_artifact = 0;
            select(state, Node::Project);
            return;
        }
        wizard.error = state.error.clone();
    }
    state.wizard = Some(wizard);
}
