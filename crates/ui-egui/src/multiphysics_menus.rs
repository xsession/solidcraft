//! Original functional menu structure. Active actions share the ribbon dispatcher.
use super::*;
pub(super) const HEIGHT: f32 = ribbon::TAB_ROW_HEIGHT;
const GROUPS: &[(&str, &str)] = &[
    ("project", "multiphysics.menu_project"),
    ("insert", "multiphysics.menu_insert"),
    ("mesh", "multiphysics.group_mesh"),
    ("solve", "multiphysics.group_solve"),
    ("results", "multiphysics.group_results"),
    ("tools", "multiphysics.group_tools"),
    ("help", "multiphysics.menu_help"),
];
fn items(group: &str) -> &'static [(&'static str, &'static str)] {
    match group {
        "project" => &[
            ("wizard", "multiphysics.cmd_wizard"),
            ("create", "multiphysics.cmd_new_project"),
            ("clone", "multiphysics.menu_clone"),
            ("general", "multiphysics.cmd_general_settings"),
            ("control", "multiphysics.cmd_calculation_control"),
            ("domain", "multiphysics.cmd_domain"),
            ("case_export", "multiphysics.cmd_export_case"),
            ("remove", "multiphysics.cmd_delete_project"),
        ],
        "insert" => &[
            ("add_from_components", "multiphysics.cmd_add_components"),
            ("initial", "multiphysics.tree_initial"),
            ("porous", "multiphysics.tree_porous"),
            ("solid_materials", "multiphysics.tree_solid_materials"),
            ("boundary", "multiphysics.cmd_boundary"),
            ("fans", "multiphysics.tree_fans"),
            ("sources", "multiphysics.cmd_sources"),
            ("radiation", "multiphysics.tree_radiation"),
            ("fluid_subdomains", "multiphysics.tree_fluid_subdomains"),
            ("contact", "multiphysics.tree_contact"),
            ("tec", "multiphysics.tree_tec"),
            ("leds", "multiphysics.tree_leds"),
            ("batteries", "multiphysics.tree_batteries"),
            ("perforated", "multiphysics.tree_perforated"),
            ("two_resistor", "multiphysics.tree_two_resistor"),
            ("electrical", "multiphysics.tree_electrical"),
            ("heat_pipes", "multiphysics.tree_heat_pipes"),
            ("thermal_joint", "multiphysics.menu_thermal_joint"),
            ("pcb", "multiphysics.tree_pcb"),
            ("rotating", "multiphysics.tree_rotating"),
            ("transferred", "multiphysics.tree_transferred"),
            ("goals", "multiphysics.cmd_goals"),
            ("local_mesh", "multiphysics.menu_local_mesh"),
        ],
        "mesh" => &[("mesh", "multiphysics.cmd_mesh"), ("local_mesh", "multiphysics.menu_local_mesh"), ("check", "multiphysics.cmd_check")],
        "solve" => &[
            ("run_solver", "multiphysics.cmd_run"),
            ("batch_run", "multiphysics.cmd_batch_run"),
            ("monitor", "multiphysics.cmd_monitor"),
            ("stop_solver", "multiphysics.cmd_stop"),
        ],
        "results" => &[
            ("load_results", "multiphysics.cmd_load_results"),
            ("unload_results", "multiphysics.menu_unload"),
            ("cut_plot", "multiphysics.tree_cut_plots"),
            ("surface_plot", "multiphysics.tree_surface_plots"),
            ("isosurface", "multiphysics.tree_isosurfaces"),
            ("trajectories", "multiphysics.tree_trajectories"),
            ("particles", "multiphysics.tree_particles"),
            ("point_params", "multiphysics.tree_point_params"),
            ("surface_params", "multiphysics.tree_surface_params"),
            ("volume_params", "multiphysics.tree_volume_params"),
            ("xy_plot", "multiphysics.tree_xy_plots"),
            ("goal_plot", "multiphysics.tree_goal_plots"),
            ("animation", "multiphysics.tree_animation"),
            ("compare", "multiphysics.tree_compare"),
            ("report", "multiphysics.tree_report"),
            ("result_export", "multiphysics.cmd_export_results"),
        ],
        "tools" => &[
            ("materials", "multiphysics.cmd_engineering_database"),
            ("templates", "multiphysics.cmd_templates"),
            ("customize_tree", "multiphysics.menu_customize"),
            ("fluid_volume", "multiphysics.menu_fluid_volume"),
            ("create_lids", "multiphysics.menu_create_lids"),
            ("verify", "multiphysics.cmd_verify"),
        ],
        "help" => &[("templates", "multiphysics.menu_capabilities"), ("verify", "multiphysics.menu_validation")],
        _ => &[],
    }
}
fn reason(app: &SolveApp, state: &State, study: Option<&Study>, id: &str) -> Option<&'static str> {
    if matches!(id, "wizard" | "create" | "materials" | "templates" | "customize_tree" | "batch_run") {
        return None;
    }
    let Some(study) = study else { return Some("multiphysics.menu_needs_project") };
    if matches!(
        id,
        "load_results"
            | "unload_results"
            | "cut_plot"
            | "point_params"
            | "volume_params"
            | "xy_plot"
            | "goal_plot"
            | "animation"
            | "compare"
            | "report"
            | "result_export"
            | "transferred"
    ) && study.runs.is_empty()
    {
        return Some("multiphysics.menu_needs_results");
    }
    let running = app.session.multiphysics_jobs().iter().any(|job| job.study_id == study.id && matches!(job.state.as_str(), "queued" | "running"));
    if (id == "run_solver" && running) || (id == "stop_solver" && !running) {
        return Some("multiphysics.menu_job_state");
    }
    if matches!(
        id,
        "fans"
            | "porous"
            | "create_lids"
            | "add_from_components"
            | "local_mesh"
            | "transferred"
            | "perforated"
            | "pcb"
            | "leds"
            | "batteries"
            | "contact"
            | "tec"
            | "two_resistor"
            | "heat_pipes"
            | "thermal_joint"
            | "rotating"
            | "fluid_volume"
    ) {
        let case = state.drafts.get(&study.id).and_then(|draft| serde_json::from_str::<Value>(&draft.case_text).ok());
        let case = case.as_ref().unwrap_or(&study.case);
        if case["solver"].as_str() != Some("collocated_flow3d") {
            return Some("multiphysics.feature_solver_required");
        }
        // Thermal devices act on the conjugate-heat block's solids.
        if !matches!(
            id,
            "fans" | "porous" | "create_lids" | "add_from_components" | "local_mesh" | "transferred" | "perforated" | "rotating" | "fluid_volume"
        ) && case.get("thermal").is_none_or(Value::is_null)
        {
            return Some("multiphysics.feature_thermal_required");
        }
    }
    None
}
pub(super) fn show(app: &mut SolveApp, state: &mut State, studies: &[Study], ui: &mut egui::Ui) {
    let study = state.selected.as_ref().and_then(|id| studies.iter().find(|study| &study.id == id));
    let mut chosen = None;
    egui::ScrollArea::horizontal().id_salt("multiphysics_menus").show(ui, |ui| {
        ui.horizontal(|ui| {
            for &(group, key) in GROUPS {
                let response = crate::ui_kit::Menu::button(ui, tr(key), |ui| {
                    egui::ScrollArea::vertical().max_height(ui.ctx().content_rect().height() * 0.75).show(ui, |ui| {
                        for &(id, label) in items(group) {
                            let why = reason(app, state, study, id);
                            let identity = format!("multiphysics_menu_{group}_{id}");
                            let mut item = crate::ui_kit::MenuItem::new(&identity, tr(label), node_for(id).map_or("settings", Node::icon));
                            item.enabled = why.is_none();
                            let response = item.show(ui);
                            crate::scenario::publish_target(&format!("multiphysics:menu:{group}:{id}"), response.rect);
                            if let Some(why) = why {
                                response.clone().backing_on_hover_text(tr(why));
                            }
                            if response.clicked() {
                                chosen = Some(id);
                                ui.close();
                            }
                        }
                    });
                });
                crate::scenario::publish_target(&format!("multiphysics:menu:{group}"), response.response.rect);
            }
        });
    });
    if let Some(id) = chosen {
        run(app, state, studies, id);
    }
}
