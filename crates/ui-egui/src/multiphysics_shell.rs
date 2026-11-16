//! The Multiphysics workspace laid out as a CAD-embedded flow-analysis environment: a ribbon of
//! grouped commands across the top, a resource bar with the Analysis tree (Input Data and
//! Results folders of the active project) on the left, the CAD graphics window in the centre,
//! and feature dialogs, result windows, the solver monitor and the project wizard floating over
//! the graphics window. Everything here routes to the workspace's existing editors and commands;
//! it adds no numerical behaviour of its own.

use super::*;
use crate::ribbon;
use egui::{Rect, pos2, vec2};

/// Width of the resource bar's icon strip.
pub(super) const RESOURCE_BAR: f32 = 34.0;
/// Width kept for the workspace switcher at the left of the ribbon.
const SWITCHER: f32 = 88.0;
/// The ribbon: its tab row and one row of labelled groups.
pub(super) const RIBBON: f32 = ribbon::TAB_ROW_HEIGHT + ribbon::GROUP_HEIGHT + 1.0;

/// One entry of the Analysis tree. Folders and features map onto the workspace stages, so
/// scripts that select a stage and people who click the tree reach the same editors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Node {
    Project,
    Input,
    General,
    Control,
    Domain,
    FluidSubdomains,
    SolidMaterials,
    Boundary,
    Fans,
    Sources,
    Radiation,
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
    Coupling,
    Results,
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
    Animation,
    Compare,
    Verify,
}

impl Node {
    /// The Input Data folder's features, in tree order.
    /// The Input data folder's features, in tree order: the flow-analysis defaults
    /// (computational domain, fluid subdomains, boundary conditions, goals, mesh) and the
    /// folders that show once the case has entries for them or Customize Tree turns them on.
    const INPUT: [Self; 23] = [
        Self::Domain,
        Self::FluidSubdomains,
        Self::SolidMaterials,
        Self::Boundary,
        Self::Fans,
        Self::Sources,
        Self::Radiation,
        Self::Porous,
        Self::Perforated,
        Self::Contact,
        Self::Initial,
        Self::Rotating,
        Self::TwoResistor,
        Self::HeatPipes,
        Self::Pcb,
        Self::Tec,
        Self::Leds,
        Self::Electrical,
        Self::Batteries,
        Self::Transferred,
        Self::Goals,
        Self::Mesh,
        Self::Coupling,
    ];
    /// The Results folder's features, in tree order. The unavailable ones stay listed (greyed)
    /// so the tree reads the same as a full flow-analysis tool and says what is missing.
    const RESULTS: [Self; 13] = [
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
        matches!(self, Self::Domain | Self::FluidSubdomains | Self::Boundary | Self::Goals | Self::Mesh)
            || (Self::RESULTS.contains(&self) && self.available())
    }

    /// Whether the tree lists this folder now: Customize Tree's choice, or the default, and
    /// always when the case has entries for it.
    fn visible(self, state: &State, case: &Value) -> bool {
        state.tree_custom.get(self.key()).copied().unwrap_or_else(|| self.default_visible()) || case_has(case, self)
    }

    /// Whether the node has an implementation; the others are listed for orientation only.
    pub(super) fn available(self) -> bool {
        !matches!(self, Self::SurfacePlots | Self::Isosurfaces | Self::Trajectories | Self::Particles | Self::SurfaceParams)
    }

    /// Shown in the graphics window rather than in a floating window over it.
    fn in_graphics(self) -> bool {
        matches!(self, Self::Results | Self::CutPlots | Self::Animation | Self::Verify)
    }

    /// The case editor's view of this node: which top-level case entries it lists.
    pub(super) fn edits_case(self) -> bool {
        matches!(self, Self::Input | Self::General | Self::Control | Self::Mesh) || (Self::INPUT.contains(&self) && self != Self::Coupling)
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

/// The New Project wizard: name, analysis type, template, summary.
#[derive(Clone, Debug, Default)]
pub(super) struct Wizard {
    step: usize,
    name: String,
    category: Option<String>,
    template: Option<usize>,
}

impl Wizard {
    /// The wizard on its summary page for a template picked in the template list.
    pub(super) fn from_template(index: usize, category: String) -> Self {
        Self { step: WIZARD_STEPS.len() - 1, name: String::new(), category: Some(category), template: Some(index) }
    }
}

static WIZARD_STEPS: [&str; 4] =
    ["multiphysics.wizard_step_name", "multiphysics.wizard_step_type", "multiphysics.wizard_step_template", "multiphysics.wizard_step_summary"];

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
    let node = current(state);
    let running = study.is_some_and(|study| {
        app.session.multiphysics_jobs().iter().any(|job| job.study_id == study.id && matches!(job.state.as_str(), "queued" | "running"))
    });
    let mut out = if tab == 0 {
        vec![
            (
                "project",
                "multiphysics.group_project",
                vec![
                    cmd("wizard", "multiphysics.cmd_wizard", "new", true),
                    cmd("create", "multiphysics.cmd_new_project", "plus", true),
                    cmd("general", "multiphysics.cmd_general_settings", "settings", has_study),
                    cmd("control", "multiphysics.cmd_calculation_control", "timeline_play", has_study),
                    cmd("remove", "multiphysics.cmd_delete_project", "delete", has_study),
                ],
            ),
            (
                "input",
                "multiphysics.group_input",
                vec![
                    cmd("domain", "multiphysics.cmd_domain", "box", has_study),
                    cmd("boundary", "multiphysics.cmd_boundary", "surface", has_study),
                    cmd("sources", "multiphysics.cmd_sources", "render_light_point", has_study),
                    cmd("materials", "multiphysics.cmd_engineering_database", "material", true),
                    cmd("goals", "multiphysics.cmd_goals", "checkpoint", has_study),
                    cmd("coupling", "multiphysics.cmd_coupling", "joint", has_study),
                ],
            ),
            (
                "mesh",
                "multiphysics.group_mesh",
                vec![cmd("mesh", "multiphysics.cmd_mesh", "mesh_body", has_study), cmd("check", "multiphysics.cmd_check", "render_check", has_study)],
            ),
            (
                "solve",
                "multiphysics.group_solve",
                vec![
                    cmd("run_solver", "multiphysics.cmd_run", "play", has_study && !running),
                    cmd("monitor", "multiphysics.cmd_monitor", "settings_monitor", has_study),
                    cmd("stop_solver", "multiphysics.cmd_stop", "stop", running),
                ],
            ),
            (
                "tools",
                "multiphysics.group_tools",
                vec![
                    cmd("templates", "multiphysics.cmd_templates", "library", true),
                    cmd("verify", "multiphysics.cmd_verify", "inspect", has_study),
                    cmd("case_export", "multiphysics.cmd_export_case", "export_tray", has_study),
                ],
            ),
        ]
    } else {
        vec![
            ("load", "multiphysics.group_results", vec![cmd("load_results", "multiphysics.cmd_load_results", "open", has_runs)]),
            (
                "plots",
                "multiphysics.group_plots",
                vec![
                    cmd("cut_plot", "multiphysics.tree_cut_plots", "section", has_runs),
                    cmd("surface_plot", "multiphysics.tree_surface_plots", "surface", false),
                    cmd("isosurface", "multiphysics.tree_isosurfaces", "iso_analysis", false),
                    cmd("trajectories", "multiphysics.tree_trajectories", "spline", false),
                    cmd("particles", "multiphysics.tree_particles", "dots", false),
                ],
            ),
            (
                "parameters",
                "multiphysics.group_parameters",
                vec![
                    cmd("point_params", "multiphysics.tree_point_params", "point", has_runs),
                    cmd("surface_params", "multiphysics.tree_surface_params", "plane", false),
                    cmd("volume_params", "multiphysics.tree_volume_params", "body", has_runs),
                ],
            ),
            (
                "charts",
                "multiphysics.group_charts",
                vec![
                    cmd("xy_plot", "multiphysics.tree_xy_plots", "line", has_runs),
                    cmd("goal_plot", "multiphysics.tree_goal_plots", "graph", has_runs),
                ],
            ),
            (
                "output",
                "multiphysics.group_output",
                vec![
                    cmd("report", "multiphysics.tree_report", "list", has_runs),
                    cmd("animation", "multiphysics.tree_animation", "ws_animate", has_runs),
                    cmd("compare", "multiphysics.tree_compare", "compare", has_runs),
                    cmd("result_export", "multiphysics.cmd_export_results", "export", has_runs),
                ],
            ),
        ]
    };
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

/// The tree node a ribbon command opens.
fn node_for(id: &str) -> Option<Node> {
    Some(match id {
        "general" => Node::General,
        "control" => Node::Control,
        "domain" => Node::Domain,
        "boundary" => Node::Boundary,
        "sources" => Node::Sources,
        "goals" => Node::Goals,
        "coupling" => Node::Coupling,
        "mesh" => Node::Mesh,
        "check" | "case_export" => Node::Project,
        "verify" => Node::Verify,
        "load_results" => Node::Results,
        "cut_plot" => Node::CutPlots,
        "point_params" => Node::PointParams,
        "volume_params" => Node::VolumeParams,
        "xy_plot" => Node::XyPlots,
        "goal_plot" => Node::GoalPlots,
        "report" | "result_export" => Node::Report,
        "animation" => Node::Animation,
        "compare" => Node::Compare,
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

static TABS: [(&str, &str); 2] = [("ANALYSIS", "multiphysics.ribbon_tab_analysis"), ("RESULTS", "multiphysics.ribbon_tab_results")];

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
    let tab_row = Rect::from_min_max(pos2(left, whole.top()), pos2(whole.right(), whole.top() + ribbon::TAB_ROW_HEIGHT));
    let labels: Vec<&str> = TABS.iter().map(|(_, key)| tr(key)).collect();
    let quick = [("save", tr("multiphysics.save")), ("undo", tr("multiphysics.undo")), ("redo", tr("multiphysics.redo"))];
    let accent = crate::workspace::tab_accent(&t, "MODEL");
    let active = TABS.get(state.ribbon_tab).map_or("ANALYSIS", |&(id, _)| id);
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
    let group_labels: Vec<&str> = groups.iter().map(|(_, key, _)| tr(key)).collect();
    let ids: Vec<Vec<&str>> = groups.iter().map(|(_, _, commands)| commands.iter().map(|command| command.id).collect()).collect();
    let more_tip = tr("multiphysics.ribbon_more");
    let row = ribbon::Row {
        groups: groups
            .iter()
            .enumerate()
            .map(|(gi, (key, _, commands))| ribbon::Group {
                key: *key,
                name: group_labels.get(gi).copied().unwrap_or_default(),
                label: group_labels.get(gi).copied().unwrap_or_default(),
                tip: group_labels.get(gi).copied().unwrap_or_default(),
                buttons: commands
                    .iter()
                    .zip(texts.get(gi).into_iter().flatten())
                    .map(|(command, &label)| ribbon::Button {
                        id: command.id,
                        run: command.id,
                        label,
                        name: label,
                        icon: command.icon,
                        size: ribbon::ButtonSize { width: ribbon::button_width(&painter, label), min_width: ribbon::BUTTON_MIN_WIDTH },
                        enabled: command.enabled,
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
        if let Some(command) = groups.get(button.group).and_then(|(_, _, commands)| commands.get(button.index)) {
            if command.id == "monitor" {
                crate::scenario::publish_target("multiphysics.stage_run", button.whole);
            }
            if !command.enabled && unavailable(command.id) {
                let _ = button.response.clone().backing_on_hover_text(tr("multiphysics.tree_unavailable"));
            }
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

/// The result tools that are listed but not built yet.
fn unavailable(id: &str) -> bool {
    matches!(id, "surface_plot" | "isosurface" | "trajectories" | "particles" | "surface_params")
}

/// What a ribbon command does.
fn run(app: &mut SolveApp, state: &mut State, studies: &[Study], id: &str) {
    let study = state.selected.as_ref().and_then(|selected| studies.iter().find(|study| &study.id == selected)).cloned();
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
enum TreeAction {
    Select(String, Node),
    Entry(String, Node, String),
    Artifact(String, Node, usize),
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
}

/// A context-menu entry (kit ghost button) that closes the menu when chosen.
fn menu_item(ui: &mut egui::Ui, key: &str) -> bool {
    let chosen = crate::ui_kit::button(ui, tr(key), crate::ui_kit::Button::Ghost).clicked();
    if chosen {
        ui.close();
    }
    chosen
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
            response.context_menu(|ui| {
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
            for setting in [Node::General, Node::Control] {
                let response = row(ui, 1, None, setting.icon(), tr(setting.key()), node == setting, true, false);
                publish_node(setting, response.rect);
                if response.clicked() {
                    actions.push(TreeAction::Select(study.id.clone(), setting));
                }
                response.context_menu(|ui| {
                    if menu_item(ui, "multiphysics.tree_edit_definition") {
                        actions.push(TreeAction::Select(study.id.clone(), setting));
                    }
                });
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
            response.context_menu(|ui| {
                if menu_item(ui, "multiphysics.tree_edit_definition") {
                    actions.push(TreeAction::Select(study.id.clone(), Node::Input));
                }
                if menu_item(ui, "multiphysics.tree_customize") {
                    actions.push(TreeAction::Customize);
                }
            });
            if open {
                for folder in Node::INPUT {
                    if !folder.visible(state, &case) {
                        continue;
                    }
                    let pointers = entries(&case, folder);
                    let folded = state.tree_folded.contains(folder.key());
                    let has_rows = !pointers.is_empty();
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
                    response.context_menu(|ui| {
                        if menu_item(ui, "multiphysics.tree_edit_definition") {
                            actions.push(TreeAction::Select(study.id.clone(), folder));
                        }
                        if menu_item(ui, "multiphysics.tree_hide_folder") {
                            actions.push(TreeAction::Hide(folder));
                        }
                        if menu_item(ui, "multiphysics.tree_customize") {
                            actions.push(TreeAction::Customize);
                        }
                    });
                    if folded || folder == Node::Coupling {
                        continue;
                    }
                    for pointer in pointers {
                        let picked = node == folder && state.entry.as_deref() == Some(pointer.as_str());
                        let response = row(ui, 3, None, "params", &entry_label(&pointer), picked, true, false);
                        crate::scenario::publish_target(&format!("multiphysics:tree_entry:{pointer}"), response.rect);
                        if response.clicked() {
                            actions.push(TreeAction::Entry(study.id.clone(), folder, pointer.clone()));
                        }
                        entry_menu(&response, study, folder, &pointer, &case, &mut actions);
                        for (item, label) in items(&case, &pointer) {
                            let picked = node == folder && state.entry.as_deref() == Some(item.as_str());
                            let response = row(ui, 4, None, "component", &label, picked, true, false);
                            crate::scenario::publish_target(&format!("multiphysics:tree_entry:{item}"), response.rect);
                            if response.clicked() {
                                actions.push(TreeAction::Entry(study.id.clone(), folder, item.clone()));
                            }
                            entry_menu(&response, study, folder, &item, &case, &mut actions);
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
                None => tr(Node::Results.key()).to_string(),
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
            response.context_menu(|ui| {
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
                for folder in Node::RESULTS {
                    if !folder.visible(state, &Value::Null) {
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
                        actions.push(TreeAction::Select(study.id.clone(), folder));
                    }
                    response.context_menu(|ui| {
                        if folder.available() && menu_item(ui, "multiphysics.tree_insert") {
                            actions.push(TreeAction::Select(study.id.clone(), folder));
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
                            actions.push(TreeAction::Artifact(study.id.clone(), folder, index));
                        }
                        response.context_menu(|ui| {
                            if menu_item(ui, "multiphysics.tree_show") {
                                actions.push(TreeAction::Artifact(study.id.clone(), folder, index));
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

/// The right-click menu of a case entry or list item.
fn entry_menu(response: &egui::Response, study: &Study, folder: Node, pointer: &str, case: &Value, actions: &mut Vec<TreeAction>) {
    response.context_menu(|ui| {
        if menu_item(ui, "multiphysics.tree_edit_definition") {
            actions.push(TreeAction::Entry(study.id.clone(), folder, pointer.to_string()));
        }
        if menu_item(ui, "multiphysics.tree_revert") {
            actions.push(TreeAction::Revert(study.id.clone(), pointer.to_string()));
        }
        if menu_item(ui, "multiphysics.tree_copy_json") {
            actions.push(TreeAction::Copy(case.pointer(pointer).map(pretty).unwrap_or_default()));
        }
    });
}

fn apply_tree_action(app: &mut SolveApp, state: &mut State, studies: &[Study], action: TreeAction, ctx: &egui::Context) {
    let activate = |state: &mut State, id: String| {
        if state.selected.as_ref() != Some(&id) {
            state.selected_run = None;
            state.selected_artifact = 0;
        }
        state.selected = Some(id);
    };
    match action {
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
        TreeAction::Artifact(id, node, index) => {
            activate(state, id);
            state.ribbon_tab = 1;
            select(state, node);
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
            for (heading, folders) in [("multiphysics.tree_input", &Node::INPUT[..]), ("multiphysics.tree_results", &Node::RESULTS[..])] {
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

fn disclosure_left(depth: usize) -> f32 {
    Tokens::get().space_4 + depth as f32 * Tokens::get().space_16
}

/// One tree row: indent, an optional disclosure triangle, the icon and the label. `bold` marks
/// the active project, as flow-analysis trees do.
#[allow(clippy::too_many_arguments)]
fn row(ui: &mut egui::Ui, depth: usize, open: Option<bool>, icon: &str, text: &str, selected: bool, enabled: bool, bold: bool) -> egui::Response {
    let t = Tokens::get();
    let height = crate::ui_kit::row_height();
    let (rect, response) = ui.allocate_exact_size(vec2(ui.available_width(), height), if enabled { Sense::click() } else { Sense::hover() });
    let painter = ui.painter();
    if selected {
        painter.rect_filled(rect, t.control_radius, t.accent_soft);
    } else if response.hovered() && enabled {
        painter.rect_filled(rect, t.control_radius, t.control_fill);
    }
    let mut x = rect.left() + disclosure_left(depth);
    let glyph = t.space_12;
    if let Some(open) = open {
        let disclosure = Rect::from_center_size(pos2(x + glyph / 2.0, rect.center().y), vec2(glyph, glyph));
        crate::icons::paint_ink(painter, disclosure, if open { "disclosure_open" } else { "disclosure_closed" }, t.text_dim);
    }
    x += glyph + t.space_4;
    let icon_size = t.space_16;
    let icon_rect = Rect::from_center_size(pos2(x + icon_size / 2.0, rect.center().y), vec2(icon_size, icon_size));
    crate::icons::paint_icon(painter, icon_rect, icon, if enabled { 1.0 } else { 0.4 });
    x += icon_size + t.space_6;
    let font = if bold { crate::ui_kit::type_face::semibold(t.font_body) } else { egui::FontId::proportional(t.font_body) };
    let color = if !enabled {
        t.text_dim
    } else if selected {
        t.selected_ink
    } else {
        t.text
    };
    let room = (rect.right() - x - t.space_4).max(0.0);
    let galley = painter.layout(text.to_owned(), font, color, room);
    let galley_height = galley.size().y;
    painter.with_clip_rect(rect).galley(pos2(x, rect.center().y - galley_height / 2.0), galley, color);
    response.widget_info(|| egui::WidgetInfo::labeled(egui::WidgetType::Button, enabled, text));
    response
}

// ── Floating windows ─────────────────────────────────────────────────────────────────────────

/// Feature dialogs and result windows over the graphics window, the solver monitor, the
/// template list and the wizard.
pub(super) fn windows(
    app: &mut SolveApp,
    state: &mut State,
    study: Option<&Study>,
    report: Option<&solidcraft_engine::multiphysics::Validation>,
    graphics: Rect,
    ctx: &egui::Context,
) {
    let t = Tokens::get();
    let node = current(state);
    let at = graphics.left_top() + vec2(t.space_12, t.space_12);
    if let Some(study) = study
        && state.dialog_open
        && !node.in_graphics()
        && !state.show_materials
    {
        let title = match &state.entry {
            Some(pointer) if node.edits_case() => format!("{} — {} — {}", tr(node.key()), entry_label(pointer), study.name),
            _ => format!("{} — {}", tr(node.key()), study.name),
        };
        let closed = crate::ui_kit::Window::new(egui::Id::unique("multiphysics_feature_window"), node.icon(), &title)
            .placed(egui::Align2::LEFT_TOP, at)
            .cards(if node.stage() == Stage::Results { 2.4 } else { 1.8 })
            .page_rows(24.0)
            .keep_open_on_escape()
            .show_body(ctx, |ui| {
                crate::scenario::publish_rect("multiphysics:feature_window", ui.max_rect());
                feature(app, state, study, report, node, ui);
            });
        if closed {
            state.dialog_open = false;
        }
    }
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

/// The body of the floating window for `node`.
fn feature(
    app: &mut SolveApp,
    state: &mut State,
    study: &Study,
    report: Option<&solidcraft_engine::multiphysics::Validation>,
    node: Node,
    ui: &mut egui::Ui,
) {
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
            settings(state, study, node, entry.as_deref(), ui);
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
            crate::viewport::show(app, ui);
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
    let templates = state.templates.as_ref().and_then(|result| result.as_ref().ok()).cloned().unwrap_or_default();
    let mut categories: Vec<String> = templates.iter().map(|template| template_category(&template.name).to_owned()).collect();
    categories.sort();
    categories.dedup();
    let in_category: Vec<usize> = (0..templates.len())
        .filter(|&index| {
            wizard
                .category
                .as_deref()
                .is_none_or(|category| templates.get(index).is_some_and(|template| template_category(&template.name) == category))
        })
        .collect();
    let title = solidcraft_ui_text::format(
        "multiphysics.wizard_title",
        &[("step", &(wizard.step + 1).to_string()), ("count", &WIZARD_STEPS.len().to_string()), ("name", tr(WIZARD_STEPS[wizard.step.min(3)]))],
    );
    let mut finish = false;
    let mut cancel = false;
    let closed = crate::ui_kit::Window::new(egui::Id::unique("multiphysics_wizard"), "new", &title).cards(2.0).page_rows(18.0).show_body(ctx, |ui| {
        crate::scenario::publish_rect("multiphysics:wizard", ui.max_rect());
        match wizard.step {
            0 => {
                ui.label(tr("multiphysics.wizard_name_hint"));
                crate::ui_kit::text_row(ui, tr("multiphysics.name"), &mut wizard.name, tr("multiphysics.wizard_name_placeholder"));
            }
            1 => {
                ui.label(tr("multiphysics.wizard_type_hint"));
                egui::ScrollArea::vertical().id_salt("multiphysics_wizard_types").max_height(320.0).show(ui, |ui| {
                    if crate::ui_kit::list_picker(
                        ui,
                        "multiphysics_wizard_type_all",
                        tr("multiphysics.template_category_all"),
                        wizard.category.is_none(),
                    )
                    .clicked()
                    {
                        wizard.category = None;
                        wizard.template = None;
                    }
                    for category in &categories {
                        let selected = wizard.category.as_deref() == Some(category.as_str());
                        if crate::ui_kit::list_picker(ui, &format!("multiphysics_wizard_type_{category}"), category, selected).clicked() {
                            wizard.category = Some(category.clone());
                            wizard.template = None;
                        }
                    }
                });
            }
            2 => {
                ui.label(tr("multiphysics.wizard_template_hint"));
                egui::ScrollArea::vertical().id_salt("multiphysics_wizard_templates").max_height(320.0).show(ui, |ui| {
                    for &index in &in_category {
                        let Some(template) = templates.get(index) else { continue };
                        let selected = wizard.template == Some(index);
                        if crate::ui_kit::list_picker(ui, &format!("multiphysics_wizard_template_{}", template.id), &template.name, selected)
                            .clicked()
                        {
                            wizard.template = Some(index);
                        }
                    }
                });
                if let Some(template) = wizard.template.and_then(|index| templates.get(index)) {
                    ui.label(RichText::new(&template.profile).small().color(Tokens::get().text_dim));
                }
            }
            _ => {
                let template = wizard.template.and_then(|index| templates.get(index));
                egui::Grid::new("multiphysics_wizard_summary").striped(true).show(ui, |ui| {
                    ui.label(tr("multiphysics.name"));
                    ui.monospace(if wizard.name.trim().is_empty() {
                        template.map_or("", |template| template.name.as_str())
                    } else {
                        wizard.name.as_str()
                    });
                    ui.end_row();
                    ui.label(tr(WIZARD_STEPS[1]));
                    ui.monospace(wizard.category.as_deref().unwrap_or(tr("multiphysics.template_category_all")));
                    ui.end_row();
                    ui.label(tr(WIZARD_STEPS[2]));
                    ui.monospace(template.map_or("", |template| template.name.as_str()));
                    ui.end_row();
                    ui.label(tr("multiphysics.wizard_profile"));
                    ui.monospace(template.map_or("", |template| template.profile.as_str()));
                    ui.end_row();
                    let entries = template.and_then(|template| template.case.as_object()).map_or(0, serde_json::Map::len);
                    ui.label(tr("multiphysics.wizard_settings"));
                    ui.monospace(entries.to_string());
                    ui.end_row();
                });
                ui.label(RichText::new(tr("multiphysics.wizard_summary_hint")).small().color(Tokens::get().text_dim));
            }
        }
        ui.separator();
        ui.horizontal(|ui| {
            let ready = match wizard.step {
                2 => wizard.template.is_some(),
                3 => wizard.template.is_some(),
                _ => true,
            };
            if button(ui, "multiphysics.wizard_back", "multiphysics:wizard_back", wizard.step > 0).clicked() {
                wizard.step = wizard.step.saturating_sub(1);
            }
            if wizard.step < WIZARD_STEPS.len() - 1 {
                if button(ui, "multiphysics.wizard_next", "multiphysics:wizard_next", ready).clicked() {
                    wizard.step += 1;
                }
            } else if button(ui, "multiphysics.wizard_finish", "multiphysics:wizard_finish", ready).clicked() {
                finish = true;
            }
            if button(ui, "multiphysics.cancel_draft", "multiphysics:wizard_cancel", true).clicked() {
                cancel = true;
            }
        });
    });
    if finish && let Some(index) = wizard.template {
        create_from_template(app, state, Some((index, wizard.name.clone())));
        return;
    }
    if !closed && !cancel {
        state.wizard = Some(wizard);
    }
}
