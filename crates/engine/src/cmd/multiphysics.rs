//! Independent study authoring and asynchronous numerical work. No exact CAD geometry changes.
use super::CommandSpec;
use crate::schema::multiphysics::*;
use crate::schema::{TypedContract, output_schema_for};
use crate::{EngineError, Result, Session, multiphysics_jobs::JobSnapshot, params::bad};
use serde::{Serialize, de::DeserializeOwned};
use serde_json::Value;
use solidcraft_multiphysics::{Study, catalog as solver_catalog, preflight, templates};
use std::sync::Arc;

#[path = "multiphysics_model.rs"]
mod model;
#[path = "multiphysics_parts.rs"]
mod parts;

pub static COMMANDS: &[CommandSpec] = &[
 CommandSpec::new("multiphysics.case.feature.prepare", "Prepare Multiphysics Feature", crate::checked!("multiphysics.case.feature.prepare", prepare_feature))
 .doc_only().noundo().icon("settings").contract(&TypedContract::<FeaturePrepare>::new()).output(output_schema_for::<solidcraft_multiphysics::authoring::PreparedFeature>)
 .params("case, kind: fan|porous|solid_material|goal|global_goal|point_goal|volume_goal|equation_goal|contact_resistance|two_resistor|heat_pipe|thermoelectric_cooler, name (new porous/goal/cooler identifier or existing thermal solid; reserved for singleton fan), replace_fan_patch_conditions?: false. Returns an unapplied reviewable candidate; no document mutation or solver run."),
    CommandSpec::new("multiphysics.case.component.import", "Add Components to Multiphysics Case", crate::checked!("multiphysics.case.component.import", import_components))
        .doc_only().noundo().icon("import").contract(&TypedContract::<ComponentImport>::new()).output(output_schema_for::<solidcraft_multiphysics::case_scene::Imported>)
        .params("case: collocated_flow3d case; bodies: model body names (empty: all); quality?: standard|fine|finest; options: {representation: auto|obstacle|surface, origin_mm?, fit_domain, margin_mm, cell_mm?, replace_parts}. Tessellates the bodies (mm) into grid blocks or surface bodies. Returns an unapplied candidate and the CAD origin used; no document mutation."),
    CommandSpec::new("multiphysics.cad.bodies", "List Bodies for Multiphysics", crate::checked!("multiphysics.cad.bodies", cad_bodies))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<CadBodies>)
        .params("No parameters. Names of the current model's bodies, as multiphysics.case.component.import takes them."),
    CommandSpec::new("multiphysics.surface.inspect", "Inspect Multiphysics Surface", crate::checked!("multiphysics.surface.inspect", inspect_surface))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<SurfaceInput>::new()).output(output_schema_for::<solidcraft_multiphysics::mesh_diagnostics::SurfaceReport>)
        .params("mesh: bounded indexed triangle surface. Reports index connectivity only; does not establish watertight geometry or solver readiness."),
    CommandSpec::new("multiphysics.catalog", "Multiphysics Catalogue", crate::checked!("multiphysics.catalog", catalog))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<CatalogReply>)
        .params("No parameters. Original 34-profile inventory, 56 bounded Rust case templates and separate host status. Inventory is not execution qualification."),
    CommandSpec::new("multiphysics.study.list", "List Multiphysics Studies", crate::checked!("multiphysics.study.list", list))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<ProjectReply>)
        .params("No parameters. Document-owned study and run summaries with case hashes; case payloads, manifests and artifact bytes are omitted."),
    CommandSpec::new("multiphysics.study.create", "Create Multiphysics Study", crate::checked!("multiphysics.study.create", create))
        .doc_only().icon("new").contract(&TypedContract::<Create>::new()).output(output_schema_for::<StudyReply>)
        .params("template: known case id; name: nonblank name; case?: reviewed template-compatible candidate. Creates one saved study atomically in one undo step without CAD-to-mesh conversion or execution."),
    CommandSpec::new("multiphysics.study.update", "Apply Multiphysics Study", crate::checked!("multiphysics.study.update", update))
        .doc_only().icon("settings").contract(&TypedContract::<Update>::new()).output(output_schema_for::<StudyReply>)
        .params("id, expected_revision, name, question, case: bounded JSON object. Atomic Apply; detailed solver schema and physical prerequisites are checked at execution."),
    CommandSpec::new("multiphysics.study.sweep", "Sweep Multiphysics Study", crate::checked!("multiphysics.study.sweep", sweep))
        .doc_only().icon("new").contract(&TypedContract::<Sweep>::new()).output(output_schema_for::<SweepReply>)
        .params("id: source study; pointer: JSON pointer to an existing scalar in its case (not /solver); values: 1..16 values of the same type. Creates one checked study per value in one undo step; nothing runs."),
    CommandSpec::new("multiphysics.study.clone", "Clone Multiphysics Study", crate::checked!("multiphysics.study.clone", clone_study))
        .doc_only().icon("copy").contract(&TypedContract::<CloneStudy>::new()).output(output_schema_for::<CloneReply>)
        .params("id: source study; name: optional (default: \"<name> (copy)\"). Copies the committed case, question and setup into a new study in one undo step; retained runs stay with the source."),
    CommandSpec::new("multiphysics.wizard.enclosure", "Multiphysics Enclosure Wizard", crate::checked!("multiphysics.wizard.enclosure", wizard))
        .doc_only().icon("new").contract(&TypedContract::<Wizard>::new()).output(output_schema_for::<WizardReply>)
        .params("name, question, spec: enclosure size (mm), cooling, ambient (°C), wall convection, components (mm, W) and mesh level 1..5. Creates one study and reports every automatic choice; nothing runs."),
    CommandSpec::new("multiphysics.study.remove", "Remove Multiphysics Study", crate::checked!("multiphysics.study.remove", remove))
        .doc_only().icon("delete").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<Removed>)
        .params("id. Removes authored study in one undo step. External jobs are retained separately."),
    CommandSpec::new("multiphysics.preflight", "Check Multiphysics Study", crate::checked!("multiphysics.preflight", check))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<PreflightReply>)
        .params("id. Bounded structural check and separate worker status; not numerical or physical validation."),
    CommandSpec::new("multiphysics.run.start", "Run Multiphysics Case", crate::checked!("multiphysics.run.start", start))
        .doc_only().noundo().icon("play").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<JobSnapshot>)
        .params("id: committed study. Freezes revision and original input hash; starts injected asynchronous native companion. Refuses missing host and scratch sessions."),
    CommandSpec::new("multiphysics.run.status", "Multiphysics Run Status", crate::checked!("multiphysics.run.status", status))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<JobsReply>)
        .params("No parameters. Polls up to eight recent jobs; dataset bytes omitted to keep status replies bounded. Native UI reads shared admitted evidence."),
    CommandSpec::new("multiphysics.run.cancel", "Cancel Multiphysics Run", crate::checked!("multiphysics.run.cancel", cancel))
        .doc_only().noundo().icon("stop").contract(&TypedContract::<JobId>::new()).output(output_schema_for::<JobSnapshot>)
        .params("job_id. Cancelling persists until worker acknowledgement. Late completion after cancellation is discarded."),
    CommandSpec::new("multiphysics.run.keep", "Keep Multiphysics Result", crate::checked!("multiphysics.run.keep", keep))
        .doc_only().icon("save").contract(&TypedContract::<JobId>::new()).output(output_schema_for::<StudyReply>)
        .params("job_id: completed current job. Rechecks document, revision, original case and artifacts; retains portable evidence in one undo step."),
    CommandSpec::new("multiphysics.case.export", "Export Multiphysics Case", crate::checked!("multiphysics.case.export", export_case))
        .doc_only().noundo().icon("export").contract(&TypedContract::<ExportCase>::new()).output(output_schema_for::<Exported>)
        .params("id, path: new filename. Writes committed numerical case through atomic no-replace writer."),
    CommandSpec::new("multiphysics.result.export", "Export Multiphysics Results", crate::checked!("multiphysics.result.export", export_result))
        .doc_only().noundo().icon("export").contract(&TypedContract::<ExportResult>::new()).output(output_schema_for::<Exported>)
        .params("id, run_id, directory: new directory. Original VTU/PVD/CSV/JSON/S2P, case and manifest for ParaView; associations unchanged, no invented .pvsm or units."),
    CommandSpec::new("multiphysics.result.compare", "Compare Two Runs", crate::checked!("multiphysics.result.compare", compare_runs))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<CompareRuns>::new()).output(output_schema_for::<CompareReply>)
        .params("first, second: {id (a study of the model), run_id? (default its latest run with a thermal report)}: two runs of one study (operating points) or of two studies (design variants, such as a clone or a sweep's points); csv? (add the comparison as CSV). Each part of either study, matched by name: its highest temperature in both runs (°C), the change (K), both margins and verdicts, and whether the verdict got worse or better; the largest rise first. Counts of parts warmer, cooler, worse and better (changes under 0.05 K count as none), and the change in total power. Changes nothing."),
    CommandSpec::new("multiphysics.result.convergence", "Mesh Study Convergence", crate::checked!("multiphysics.result.convergence", convergence))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<Convergence>::new()).output(output_schema_for::<ConvergenceReply>)
        .params("fine, medium, coarse: {id (a study of the model), run_id? (default its latest run with a thermal report)}: the same study on three grids, finest first (multiphysics.model.mesh_study makes the two coarser ones). Per part all three report: its highest temperature's rise above the inflow on each grid, the observed order of accuracy, the extrapolated grid-independent value and the fine grid's convergence index (GCI, % of the fine value; Celik et al. 2008), and whether the grids converge monotonically, oscillate, have converged or are not settling yet. Notes a size ratio under 1.3 between neighbouring grids. Changes nothing."),
    CommandSpec::new("multiphysics.module.capture", "Keep Component Setup as Module", crate::checked!("multiphysics.module.capture", model::module_capture))
        .doc_only().icon("save").contract(&TypedContract::<ModuleCapture>::new()).output(output_schema_for::<ModuleReply>)
        .params("id (a study of the model), component (a component definition's name), name? (default \"<component> setup\"; replaces one of the same name), placement? (its occurrences' names joined by /; default the first). Keeps what the study and the part registry give that placement's bodies (material, power, limit, tags, attributes, a highest-temperature goal) and the heat sources on its faces, by the component's own body names, as a module of the design. One undo step."),
    CommandSpec::new("multiphysics.module.apply", "Apply Module to Every Copy", crate::checked!("multiphysics.module.apply", model::module_apply))
        .doc_only().icon("copy").contract(&TypedContract::<ModuleApply>::new()).output(output_schema_for::<ModuleApplied>)
        .params("id, expected_revision, module (name). Gives the module's setup to every placement of its component in the study: materials, tags, attributes, goals and face heat sources, replacing those of the same body or name, so applying twice changes nothing more. Lists bodies a placement does not have. Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.module.list", "List Modules", crate::checked!("multiphysics.module.list", model::module_list))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<ModulesReply>)
        .params("No parameters. The design's modules. Changes nothing."),
    CommandSpec::new("multiphysics.module.remove", "Remove Module", crate::checked!("multiphysics.module.remove", model::module_remove))
        .doc_only().icon("delete").contract(&TypedContract::<ModuleName>::new()).output(output_schema_for::<ModulesReply>)
        .params("module (name). Drops it from the design in one undo step; studies it was applied to keep their setup."),
    CommandSpec::new("multiphysics.module.export", "Export Module", crate::checked!("multiphysics.module.export", model::module_export))
        .doc_only().noundo().icon("export").contract(&TypedContract::<ModuleFile>::new()).output(output_schema_for::<Exported>)
        .params("module (name), path (a new .json file): writes the module (solidcraft-multiphysics-module/1) for another design to import. Never overwrites."),
    CommandSpec::new("multiphysics.module.import", "Import Module", crate::checked!("multiphysics.module.import", model::module_import))
        .doc_only().icon("import").contract(&TypedContract::<ModuleFile>::new()).output(output_schema_for::<ModuleReply>)
        .params("path (a module file), module? (a name to keep it under). Checks and keeps it in the design, replacing one of the same name, in one undo step; apply it to a study with multiphysics.module.apply."),
    CommandSpec::new("multiphysics.view.save", "Save Multiphysics Result View", crate::checked!("multiphysics.view.save", view_save))
        .doc_only().icon("camera").contract(&TypedContract::<ViewSave>::new()).output(output_schema_for::<ViewsReply>)
        .params("id (a study); view: {name, camera? {target [mm], yaw, pitch, distance, fov?, roll?} (none: the camera stays where it is when applied), colour_by?: none|material|power|tag|verdict, ghosted?: [body names], plot?: {field (name in the run's field file), axis: 0|1|2, anchor?: layer|part|maximum, layer?, part?, vectors?}, surfaces?, scale?: {map, reversed?, bands?, mark_out_of_range?}, range?: auto|locked|percentile|symmetric, locked_range?: [low, high]}; ghosted_now? (take the parts ghosted now). Keeps the view in the study by name, replacing one of the same name; at most 64. One undo step; runs do not go stale."),
    CommandSpec::new("multiphysics.view.apply", "Apply Multiphysics Result View", crate::checked!("multiphysics.view.apply", view_apply))
        .doc_only().noundo().icon("camera").contract(&TypedContract::<ViewName>::new()).output(output_schema_for::<ViewApplied>)
        .params("id, name. Moves the camera to the view's, ghosts exactly its parts (those the design still has; the rest are listed as missing) and colours the parts by its colour_by from the study's latest run; returns the view so the host sets its cut plot, colour scale and range. View state: not saved, not an undo step."),
    CommandSpec::new("multiphysics.view.list", "List Multiphysics Result Views", crate::checked!("multiphysics.view.list", view_list))
        .doc_only().noundo().icon("list").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<ViewsReply>)
        .params("id. The study's saved result views. Changes nothing."),
    CommandSpec::new("multiphysics.view.remove", "Remove Multiphysics Result View", crate::checked!("multiphysics.view.remove", view_remove))
        .doc_only().icon("delete").contract(&TypedContract::<ViewName>::new()).output(output_schema_for::<ViewsReply>)
        .params("id, name. Drops the study's view of that name in one undo step."),
    CommandSpec::new("multiphysics.result.report", "Batch Report over Runs", crate::checked!("multiphysics.result.report", report_runs))
        .doc_only().noundo().icon("export").contract(&TypedContract::<ReportRuns>::new()).output(output_schema_for::<ReportReply>)
        .params("id (a study of the model); run_ids? (default every run with a thermal report, oldest first; at most 64); directory? (a new directory: writes report.csv, every part's numbers per run, and report.html, a self-contained page of the runs' summary and each part's highest temperature per run, and run-N.png, each run's parts by verdict, when a picture renderer is installed); inline? (also return the texts). Give directory or inline. A summary per run: total power, hottest part, least margin, pass/warn/fail counts."),
    CommandSpec::new("multiphysics.result.inspect", "Inspect Multiphysics Result", crate::checked!("multiphysics.result.inspect", inspect_result))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<InspectResult>::new()).output(output_schema_for::<ResultInspection>)
        .params("id, run_id, artifact: one retained ASCII VTU/PVD file. Returns bounded geometry, field or collection metadata, without dataset bytes or inferred units."),
    CommandSpec::new("multiphysics.result.query", "Query Multiphysics Result", crate::checked!("multiphysics.result.query", query_result))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<solidcraft_multiphysics::results::QueryRequest>::new()).output(output_schema_for::<solidcraft_multiphysics::results::QueryResponse>)
        .params("id, run_id, query: {kind: probe|min_max|line|integrate|compare|goals|components, ...}. Probes interpolate inside cells and report the nearest cell outside; line samples outside the mesh have no value; integrals are measure weighted; compare needs identical meshes. Every answer carries the run's case hash and a stale flag."),
    CommandSpec::new("multiphysics.parts.group", "Group Parts by Designator", crate::checked!("multiphysics.parts.group", parts::group))
        .doc_only().icon("settings").contract(&TypedContract::<PartsGroup>::new()).output(output_schema_for::<PartsGroupReply>)
        .params("bodies (names) or selection (a part selection's bodies as they are now); neither: every body. A part named by a reference designator (U12, C3, R10_1: one to three letters, one to five digits) gets the kind its letters name as a tag (U, IC: ic; C: capacitor; R: resistor; L: inductor; D: diode; Q: transistor; J, P: connector; Y: crystal; F: fuse; T: transformer; K: relay; SW: switch; LED: led; TP: test_point; FB: ferrite_bead; BT: battery) and its designator as the attribute designator. Other parts are left as they are. One undo step."),
    CommandSpec::new("multiphysics.parts.list", "List Parts and Selections", crate::checked!("multiphysics.parts.list", parts::list))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<PartsReply>)
        .params("No parameters. Every body of the design with its tags and attributes, every part selection with the bodies it holds now, and recorded parts the design no longer has."),
    CommandSpec::new("multiphysics.parts.similar", "Find Similar Parts", crate::checked!("multiphysics.parts.similar", parts::similar))
        .noundo().icon("inspect").contract(&TypedContract::<PartsSimilar>::new()).output(output_schema_for::<PartsSimilarReply>)
        .params("body (a part); tolerance_pct? (default 1, at most 50); select? (put them in the design's selection). The parts whose exact volume and bounding-box sides (in any order, so a part turned a quarter turn about an axis counts) each differ from the body's by at most the tolerance, the body among them. Mesh bodies have no exact geometry and are never similar. Not an undo step."),
    CommandSpec::new("multiphysics.parts.table", "Component Table", crate::checked!("multiphysics.parts.table", parts::table))
        .doc_only().noundo().icon("list").contract(&TypedContract::<PartsTable>::new()).output(output_schema_for::<PartsTableReply>)
        .params("id (a study); run_id? (default: the latest run with a thermal report); selection? (only its bodies); all_bodies? (also the design's bodies the study leaves out); csv? (add the table as CSV). One row per body: tags, attributes, conductivity, power (in the body and on its faces), limit, the run's temperatures, margin and pass/warn/fail, with the study's total power and the part with the least margin."),
    CommandSpec::new("multiphysics.parts.tag", "Tag Parts", crate::checked!("multiphysics.parts.tag", parts::tag))
        .doc_only().icon("settings").contract(&TypedContract::<PartsTag>::new()).output(output_schema_for::<PartsReply>)
        .params("bodies (names) or selection (a part selection's bodies as they are now); add_tags?, remove_tags? (single words); attributes? {name: number | text | null (remove)}. volume_mm3 and size_mm are measured and cannot be set. One undo step."),
    CommandSpec::new("multiphysics.parts.check", "Check Study Setup", crate::checked!("multiphysics.parts.check", parts::check))
        .noundo().icon("inspect").contract(&TypedContract::<PartsCheck>::new()).output(output_schema_for::<solidcraft_multiphysics::setup_check::Checklist>)
        .params("id (a study of the model). One line per check (built, fluid, thin_parts, overlap, openings, materials, limits, power, goals) with ok|warn|error, a count, a message and the parts it is about; the fluid lines are the fluid volume check of the built case. Changes nothing."),
    CommandSpec::new("multiphysics.parts.colour", "Colour Parts by Setup", crate::checked!("multiphysics.parts.colour", parts::colour))
        .noundo().icon("appearance").contract(&TypedContract::<PartsColour>::new()).output(output_schema_for::<solidcraft_multiphysics::setup_check::Colouring>)
        .params("by: material (unassigned parts stand out) | power (on a colour map) | tag (first tag) | verdict (the run's pass/warn/fail) | none (their own looks again); id (the study; not with none); run_id? (for verdict). Draws the design's parts in swatches the app colours from its theme, with a legend. View state: not saved, not an undo step."),
    CommandSpec::new("multiphysics.parts.find", "Find Parts", crate::checked!("multiphysics.parts.find", parts::find))
        .noundo().icon("inspect").contract(&TypedContract::<PartsFind>::new()).output(output_schema_for::<PartsFindReply>)
        .params("id (a study); rules: 1 to 16, each a part rule ({rule: tag|not_tag, tag} | {rule: name, pattern} | {rule: attribute, key, compare, value?} | {rule: in_selection, name}) or a result rule ({rule: result, quantity: power_w|limit_c|max_temperature_c|mean_temperature_c|margin_k, compare: eq|ne|lt|le|gt|ge, value} | {rule: verdict, status: pass|warn|fail|no_limit}); a part is found when it meets every rule. run_id? (default the latest with a thermal report); all_bodies?; select? (put them in the design's selection). Not an undo step."),
    CommandSpec::new("multiphysics.parts.ghost", "Ghost Parts", crate::checked!("multiphysics.parts.ghost", parts::ghost))
        .noundo().icon("eye").contract(&TypedContract::<PartsGhost>::new()).output(output_schema_for::<PartsGhostReply>)
        .params("bodies (names) or selection (a part selection's bodies as they are now); ghost?: true (false: solid again); only?: false (true: every other part solid; with no bodies or selection, every part solid). Ghosted parts are drawn see-through and clicks reach what is behind them. View state: not saved, not an undo step."),
    CommandSpec::new("multiphysics.parts.select", "Select Parts", crate::checked!("multiphysics.parts.select", parts::select))
        .noundo().icon("select").contract(&TypedContract::<PartsSelect>::new()).output(output_schema_for::<PartsSelectReply>)
        .params("bodies (names) or selection (a part selection's bodies as they are now; none: clear); add?: false (add to what is selected). Sets the design's selection, the one the Multiphysics tree, the component table and the 3-D view share; hide or isolate it with browser.visibility. Not an undo step."),
    CommandSpec::new("multiphysics.parts.selection", "Part Selection", crate::checked!("multiphysics.parts.selection", parts::selection))
        .doc_only().icon("settings").contract(&TypedContract::<PartsSelection>::new()).output(output_schema_for::<PartsReply>)
        .params("name; members: {kind: saved, bodies: [names]} | {kind: query, rules: [{rule: tag|not_tag, tag} | {rule: name, pattern (* and ?)} | {rule: attribute, key (own or volume_mm3, size_mm), compare: exists|eq|ne|lt|le|gt|ge, value?} | {rule: in_selection, name (an earlier selection)}]} (a part meets every rule); or remove: true. A query is worked out again each time it is used; model studies use one with multiphysics.model.set selection. One undo step."),
    CommandSpec::new("multiphysics.model.create", "New Study of the Model", crate::checked!("multiphysics.model.create", model::create))
        .doc_only().icon("new").contract(&TypedContract::<ModelCreate>::new()).output(output_schema_for::<ModelReply>)
        .params("name; template? (a 3-D flow template whose fluid and solver settings it starts from, default flow-channel); analysis?: internal (fluid inside the bodies; the domain is their box and its faces are walls) | external (around them, with margin_mm); bodies? (names; default every body) or selection? (a part selection, worked out again at every build); margin_mm?; cell_mm? (default from the bodies); flow?: {kinematic_viscosity_m2_per_s, density_kg_per_m3, end_time_s, converge_tolerance?} (default air, 1.5e-5 m²/s and 1.2 kg/m³, for 2 s, stopping once the goals settle to 0.1 %; the time step follows from the smallest cell, the inflow speed and the viscosity). Builds the case from the design's bodies in one undo step; nothing runs."),
    CommandSpec::new("multiphysics.model.set", "Change Study of the Model", crate::checked!("multiphysics.model.set", model::set))
        .doc_only().icon("settings").contract(&TypedContract::<ModelSet>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision; any of analysis, bodies, selection (a part selection; \"\" back to bodies), margin_mm, cell_mm (0: automatic), cut_cells (keep the fluid part of the cells the bodies cut; companion 0.14.0), mesh_level (1 coarsest to 7 finest, 3 the automatic size; 0 back to automatic; clears a typed cell_mm), flow ({kinematic_viscosity_m2_per_s, density_kg_per_m3, end_time_s, converge_tolerance?}: the fluid and how long the run simulates), porous ([{body, viscous_resistance_per_m2, inertial_resistance_per_m}]: bodies taken as porous zones over their exact boxes instead of walls, Darcy and Forchheimer; companion 0.4.0; replaces the list), internal_fans ([{name, body, direction: +x|-x|+y|-y|+z|-z, curve: [[m3_per_s, pa], ...]}]: bodies whose cells become a fan blowing along the direction by its curve; companion 0.12.0; replaces the list), rotating ([{name, bodies, axis: 0|1|2, rpm, margin_cells?}]: a rotating frame around bodies whose walls turn with it, about the axis through the middle of their box; companion 0.11.0; replaces the list). Builds the case again from the design in one undo step."),
    CommandSpec::new("multiphysics.model.leave_out_small", "Leave Out Small Parts", crate::checked!("multiphysics.model.leave_out_small", model::leave_out_small))
        .doc_only().icon("filter").contract(&TypedContract::<ModelLeaveOutSmall>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, below_mm (> 0). The study keeps its parts whose largest bounding-box side, measured from exact geometry, is at least below_mm: a live part selection \"<study> without small parts\" (its earlier parts and size_mm >= below_mm) that the study then uses; a list of bodies becomes the saved selection \"<study> parts\" first. Calling again with another size replaces the size. Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.condition", "Condition on a Model Face", crate::checked!("multiphysics.model.condition", model::condition))
        .doc_only().icon("settings").contract(&TypedContract::<ModelCondition>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, name; face: {body, face (its name) | point [x,y,z] mm on it}; role: {kind: inlet, speed_m_per_s} | {kind: outlet, pressure_m2_per_s2} | {kind: wall} | {kind: symmetry}; or remove: true. The face must be planar and lie on the computational domain's boundary (a vent's outer face in an internal study). Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.mesh_study", "Make a Mesh Study", crate::checked!("multiphysics.model.mesh_study", model::mesh_study))
        .doc_only().icon("copy").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<MeshStudyReply>)
        .params("id (a built study of the model): makes the two coarser grids of a mesh study, clones of it named \"<name> · mesh ×1.41\" and \"· mesh ×2\" with its base cells √2 and 2 times as large (local refinements kept), and builds them, in one undo step; nothing runs. Run all three, then multiphysics.result.convergence."),
    CommandSpec::new("multiphysics.model.build", "Rebuild Study from the Model", crate::checked!("multiphysics.model.build", model::build))
        .doc_only().icon("refresh").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<ModelReply>)
        .params("id. Builds the case again from the design as it is now (after a model edit), in one undo step."),
    CommandSpec::new("multiphysics.model.status", "Study of the Model Status", crate::checked!("multiphysics.model.status", model::status))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<ModelStatus>)
        .params("id. Whether the study's case was built from the design as it is now, or why it cannot be built."),
    CommandSpec::new("multiphysics.model.refine", "Refine Mesh Around the Model", crate::checked!("multiphysics.model.refine", model::refine))
        .doc_only().icon("mesh_body").contract(&TypedContract::<ModelRefine>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, name; body (a body's name) or face: {body, face | point}; factor? 2..4 (default 2); margin_cells? 0..64 (default 1); or remove: true. Splits the cells around the target (local mesh, companion 0.9.0 or later), then fits the time step to the smallest cells. Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.goal", "Goal on the Model", crate::checked!("multiphysics.model.goal", model::goal))
        .doc_only().icon("checkpoint").contract(&TypedContract::<ModelGoalInput>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, name (letters, digits, underscores); kind: face_volume_flow | face_mean_pressure (with face: {body, face | point}) | point_pressure (with point_mm) | equation (with expression naming other goals) | body_max_temperature | body_mean_temperature (with body; heat on) | point_temperature (with point_mm; heat on) | global_max_speed | global_mean_speed | global_max_pressure | global_min_pressure | global_mean_pressure (whole fluid; companion 0.16.0) | global_max_fluid_temperature | global_mean_fluid_temperature | global_max_solid_temperature (heat on); or remove: true. A face goal reads the face's condition patch, the domain side it lies on, or a wall patch made on its body. Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.heat", "Heat Transfer in the Model Study", crate::checked!("multiphysics.model.heat", model::heat))
        .doc_only().icon("settings").contract(&TypedContract::<ModelHeatSet>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision; any of enabled (heat transfer on or off), inflow_temperature_k, fluid: {specific_heat_j_per_kg_k, conductivity_w_per_m_k} (air by default). With heat on, the bodies conduct heat (conjugate heat transfer) and temperature is solved after the flow; cut cells and inlets or outlets on faces inside the domain are not available with it yet. buoyancy? (natural convection: true for air under the design's gravity, Z up; {gravity_m_per_s2, expansion_per_k, reference_temperature_k?}; false to solve heat on the settled flow; companion 0.6.0). Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.material", "Material of Model Bodies", crate::checked!("multiphysics.model.material", model::material))
        .doc_only().icon("material").contract(&TypedContract::<ModelMaterial>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, bodies (names) or selection (its bodies as they are now); conductivity_w_per_m_k (W/(m·K)), power_w? (heat released evenly in each body), density_kg_per_m3?, specific_heat_j_per_kg_k?, max_temperature_k? (a limit reported as a margin); or remove: true. Needs heat on. Bodies without a material conduct like the starter material (1 W/(m·K)). Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.heat_source", "Heat Source on a Model Face", crate::checked!("multiphysics.model.heat_source", model::heat_source))
        .doc_only().icon("settings").contract(&TypedContract::<ModelHeatSource>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, name; face: {body, face | point}; power_w; or remove: true. Needs heat on. The power goes into the body behind the face (companion 0.15.0). Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.thermal_face", "Thermal Condition on a Model Face", crate::checked!("multiphysics.model.thermal_face", model::thermal_face))
        .doc_only().icon("settings").contract(&TypedContract::<ModelThermalFace>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, name; face: {body, face | point} (a flat face on the domain boundary); condition: {kind: adiabatic} | {kind: temperature, temperature_k} | {kind: heat_flux, w_per_m2} | {kind: convection, h_w_per_m2_k, ambient_k}; or remove: true. Needs heat on. Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.fan", "Fan on a Model Face", crate::checked!("multiphysics.model.fan", model::fan))
        .doc_only().icon("settings").contract(&TypedContract::<ModelFanInput>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision; face: {body, face | point} (a flat face of the model on the domain boundary, such as a duct's open end); curve: [[volume flow m³/s, pressure rise Pa], …] (2 to 64 points, flow rising, pressure falling); name?; total_pressure? (the curve is total pressure); or remove: true. The fan blows into the domain through that face against an outlet on the domain boundary (fan curve, companion 0.4.0). One fan per study; not with cut cells. Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.surroundings", "Surroundings of the Model Study", crate::checked!("multiphysics.model.surroundings", model::surroundings))
        .doc_only().icon("settings").contract(&TypedContract::<ModelSurroundingsInput>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision; any of breeze_m_per_s (0 for still air), breeze_from: -x|+x|-y|+y (default -x), floor (a wall under the model, default true); or remove: true. For a study of the air around the model (external), Z up: the breeze comes in at its side, leaves at the opposite one, slips along the other two; the top is open, and so is the bottom without a floor. Conditions on faces of the model still win. Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.dust", "Dust in the Model Study", crate::checked!("multiphysics.model.dust", model::dust))
        .doc_only().icon("settings").contract(&TypedContract::<ModelDustInput>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision; any of concentration_mg_per_m3 (all sizes), sizes: [{diameter_um, share}] (default coarse outdoor dust: 1, 2.5, 5, 10, 20 µm), density_kg_per_m3 (default 2000), particles_per_size (default 1000), max_time_s (default 3600), entry: {kind: air} (the air coming in at the domain's sides) | {kind: box, low_mm, high_mm} (air crossing into a box: put its sides inside an enclosure's walls to count what comes in through the openings), or inside_of: a body and inset_mm (that body's box shrunk by the wall, as a box entry); or remove: true. The companion follows the dust through the settled flow to where it lands (companion 0.19.0). Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.climate", "Climate of the Model Study", crate::checked!("multiphysics.model.climate", model::climate))
        .doc_only().icon("settings").contract(&TypedContract::<ModelClimateInput>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision; conditions: [{name, hours_per_year, ambient_c, dew_point_c, breeze_m_per_s?, dust_mg_per_m3?}] or seasons: [{name, low_c, high_c, dew_point_c, dust_mg_per_m3}] (each a night, morning, afternoon and evening, the warmest at 15:00, a quarter of the season's hours each, breeze and dust following the day); humid_limit_pct? (default 60: the surface humidity at which dust counts as humid); or remove: true. At most 16 conditions, at most 8760 hours. Builds the case again in one undo step."),
    CommandSpec::new("multiphysics.model.climate_studies", "Make the Climate's Studies", crate::checked!("multiphysics.model.climate_studies", model::climate_studies))
        .doc_only().icon("copy").contract(&TypedContract::<StudyId>::new()).output(output_schema_for::<ClimateStudiesReply>)
        .params("id (a study of the model with a climate): makes one study per condition, named \"<study> · <condition>\", each with the air at that condition's temperature (inflow, and natural convection's reference), its breeze and its dust, and builds them, in one undo step; nothing runs. Studies made before for this study are replaced. Run them, then multiphysics.result.climate."),
    CommandSpec::new("multiphysics.result.climate", "Climate Year Summary", crate::checked!("multiphysics.result.climate", model::climate_result))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<ClimateResultInput>::new()).output(output_schema_for::<ClimateResultReply>)
        .params("id (the climate's study); runs?: [{study (id or name), directory (a cfd-solve output folder)}] for conditions not run in this session; directory?: a new folder for climate-parts.csv, climate-deposits.csv (every deposit: place in mm, part, size, g a year, humid hours) and climate.json. Adds the condition studies' runs (each one's kept run, else its latest finished run in this session, else its folder) over the year by their hours: dust carried in and deposited, each part's dust (g a year), the hours its coldest point is humid or below the dew point, and the hours its dust is humid. Changes nothing."),
    CommandSpec::new("multiphysics.model.adapt", "Adapt Mesh to a Model Result", crate::checked!("multiphysics.model.adapt", model::adapt))
        .doc_only().icon("mesh_body").contract(&TypedContract::<ModelAdapt>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, run_id (a kept run of the study); artifact? (default its last flow_*.vtu); field? (default the velocity); fraction? (share of cells with the largest jumps to neighbours, above 0 to 0.5, default 0.1); regions? 1..16 (default 4); factor? 2..4 (default 2); or remove: true. Replaces the study's earlier adaptive refinements (adapt_1, adapt_2, …) with boxes around the clusters of flagged cells, builds the case again in one undo step, and replies with what it found (adapted). Run again to use the finer grid."),
    CommandSpec::new("multiphysics.model.leak", "Find a Leak in the Model", crate::checked!("multiphysics.model.leak", model::leak))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<ModelLeak>::new()).output(output_schema_for::<ModelLeakReply>)
        .params("id, from, to: two faces of the model, each {body, point: [x,y,z] mm on it}, such as a face inside an enclosure and one outside it. Each end steps half a cell off its face to the fluid side. Returns the shortest path through the fluid cells of the study's built case, or connected: false when no fluid joins them. Use External with a margin so the outside has fluid."),
    CommandSpec::new("multiphysics.study.compare", "Compare Multiphysics Studies", crate::checked!("multiphysics.study.compare", compare_studies))
        .doc_only().noundo().icon("list").contract(&TypedContract::<CompareStudies>::new()).output(output_schema_for::<solidcraft_multiphysics::results::MetricTable>)
        .params("ids: studies to tabulate (empty for all). Returns every manifest metric of each study's latest retained run, aligned by name, for sweep and what-if comparison."),
];

fn parse<T: DeserializeOwned>(params: &Value) -> Result<T> {
    serde_json::from_value(params.clone()).map_err(|e| bad("", e.to_string()))
}
fn reply(value: impl Serialize) -> Result<Value> {
    serde_json::to_value(value).map_err(|e| EngineError::Other(e.to_string()))
}
fn inspect_surface(_session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::mesh_diagnostics::{self, MAX_SURFACE_TRIANGLES, MAX_SURFACE_VERTICES};
    let mesh = params.get("mesh").ok_or_else(|| bad("mesh", "surface mesh is required"))?;
    for (field, limit) in [("vertices", MAX_SURFACE_VERTICES), ("triangles", MAX_SURFACE_TRIANGLES)] {
        let values = mesh.get(field).and_then(Value::as_array).ok_or_else(|| bad(field, "expected an array"))?;
        if values.len() > limit {
            return Err(bad(field, "surface exceeds inspection limits"));
        }
    }
    let input: SurfaceInput = parse(params)?;
    reply(mesh_diagnostics::inspect(&input.mesh).map_err(|error| bad("mesh", error.to_string()))?)
}
fn study_reply(study: &Study) -> Result<Value> {
    reply(StudyReply { study: StudyView::try_from(study).map_err(|e| bad("", e.to_string()))? })
}
fn selected(session: &Session, id: &str) -> Result<Study> {
    session.doc.multiphysics.study(id).cloned().ok_or_else(|| bad("", "Multiphysics study not found in this document"))
}
fn catalog(session: &mut Session, params: &Value) -> Result<Value> {
    let _: Empty = parse(params)?;
    reply(CatalogReply {
        catalog: solver_catalog().map_err(|e| bad("", e.to_string()))?,
        templates: templates().map_err(|e| bad("", e.to_string()))?,
        host_available: session.multiphysics_host_available(),
    })
}
fn list(session: &mut Session, params: &Value) -> Result<Value> {
    let _: Empty = parse(params)?;
    reply(ProjectReply { project: ProjectView::try_from(session.doc.multiphysics.as_ref()).map_err(|e| bad("", e.to_string()))? })
}
fn create(session: &mut Session, params: &Value) -> Result<Value> {
    let input: Create = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    // Undo must never allocate an old live identity to a different study.
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let mut study = project.create(&input.template, &input.name).map_err(|e| bad("", e.to_string()))?;
    if let Some(case) = input.case {
        if case["solver"] != study.case["solver"] {
            return Err(bad("", "wizard case must use the selected template solver"));
        }
        solidcraft_multiphysics::validate_case(&case).map_err(|e| bad("", e.to_string()))?;
        study = project.update(&study.id, study.revision, &study.name, &study.question, case).map_err(|e| bad("", e.to_string()))?;
    }
    let reply = study_reply(&study)?;
    session.multiphysics_jobs.next_study = project.next_id;
    session.doc_mut().multiphysics = Arc::new(project);
    Ok(reply)
}
fn update(session: &mut Session, params: &Value) -> Result<Value> {
    let input: Update = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    let study = project.update(&input.id, input.expected_revision, &input.name, &input.question, input.case).map_err(|e| bad("", e.to_string()))?;
    if project != *session.doc.multiphysics {
        session.doc_mut().multiphysics = Arc::new(project);
    }
    study_reply(&study)
}
fn sweep(session: &mut Session, params: &Value) -> Result<Value> {
    let input: Sweep = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let studies = project.sweep(&input.id, &input.pointer, &input.values).map_err(|e| bad("", e.to_string()))?;
    let views = studies.iter().map(StudyView::try_from).collect::<std::result::Result<Vec<_>, _>>().map_err(|e| bad("", e.to_string()))?;
    session.multiphysics_jobs.next_study = project.next_id;
    session.doc_mut().multiphysics = Arc::new(project);
    reply(SweepReply { studies: views })
}
fn clone_study(session: &mut Session, params: &Value) -> Result<Value> {
    let input: CloneStudy = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let study = project.clone_study(&input.id, &input.name).map_err(|e| bad("id", e.to_string()))?;
    let view = StudyView::try_from(&study).map_err(|e| bad("", e.to_string()))?;
    session.multiphysics_jobs.next_study = project.next_id;
    session.doc_mut().multiphysics = Arc::new(project);
    reply(CloneReply { study: view })
}
fn wizard(session: &mut Session, params: &Value) -> Result<Value> {
    let input: Wizard = parse(params)?;
    let (case, report) = solidcraft_multiphysics::wizard::enclosure_case(&input.spec).map_err(|e| bad("spec", e.to_string()))?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let study = project.create_from_case(&input.name, &input.question, case).map_err(|e| bad("", e.to_string()))?;
    let view = StudyView::try_from(&study).map_err(|e| bad("", e.to_string()))?;
    session.multiphysics_jobs.next_study = project.next_id;
    session.doc_mut().multiphysics = Arc::new(project);
    reply(WizardReply { study: view, report })
}
fn remove(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.remove(&input.id).map_err(|e| bad("", e.to_string()))?;
    session.doc_mut().multiphysics = Arc::new(project);
    reply(Removed { removed: input.id })
}
fn check(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    let study = selected(session, &input.id)?;
    reply(PreflightReply {
        validation: preflight(&study),
        host_available: session.multiphysics_host_available(),
        execution_scope: "Bounded cfd-solvers-rs 0.2.0 case, independent of CAD; full case schema checked at execution".into(),
    })
}
fn start(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    if selected(session, &input.id)?.model.is_some() {
        match model::current(session, &input.id)?.1 {
            Ok(true) => {}
            Ok(false) => {
                return Err(bad("id", "the model changed since this study's case was built; rebuild it (multiphysics.model.build) before running"));
            }
            Err(problem) => return Err(bad("id", format!("the study cannot be built from the model: {problem}"))),
        }
    }
    reply(session.start_multiphysics_run(&input.id)?)
}
fn status(session: &mut Session, params: &Value) -> Result<Value> {
    let _: Empty = parse(params)?;
    session.poll_multiphysics_jobs();
    let mut jobs = session.multiphysics_jobs();
    for job in &mut jobs {
        job.evidence = None;
    }
    reply(JobsReply { jobs })
}
fn cancel(session: &mut Session, params: &Value) -> Result<Value> {
    let input: JobId = parse(params)?;
    let mut job = session.cancel_multiphysics_run(&input.job_id)?;
    job.evidence = None;
    reply(job)
}
fn keep(session: &mut Session, params: &Value) -> Result<Value> {
    let input: JobId = parse(params)?;
    session.poll_multiphysics_jobs();
    let job =
        session.multiphysics_jobs().into_iter().find(|j| j.id == input.job_id).ok_or_else(|| bad("", "Numerical job not found in this document"))?;
    if job.stale {
        return Err(bad("", "The study changed; this historical result cannot become current evidence"));
    }
    let evidence = job.evidence.ok_or_else(|| bad("", "Wait for an admitted completed result"))?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    let study = project.attach_evidence(&job.study_id, evidence.as_ref().clone()).map_err(|e| bad("", e.to_string()))?;
    if project != *session.doc.multiphysics {
        session.doc_mut().multiphysics = Arc::new(project);
    }
    study_reply(&study)
}
fn export_case(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ExportCase = parse(params)?;
    let study = selected(session, &input.id)?;
    solidcraft_io::multiphysics::write_case(std::path::Path::new(&input.path), &study).map_err(|e| bad("", e.to_string()))?;
    reply(Exported { files: vec![input.path] })
}
fn export_result(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ExportResult = parse(params)?;
    let study = selected(session, &input.id)?;
    let run = study.runs.iter().find(|r| r.id == input.run_id).ok_or_else(|| bad("", "Select a retained result"))?;
    let files =
        solidcraft_io::multiphysics::write_results(std::path::Path::new(&input.directory), &study, run).map_err(|e| bad("", e.to_string()))?;
    reply(Exported { files })
}
fn report_runs(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ReportRuns = parse(params)?;
    if input.directory.is_none() && !input.inline {
        return Err(bad("directory", "give a new directory to write the report into, or inline: true"));
    }
    let study = selected(session, &input.id)?;
    let state = session.world_state()?;
    let parts = &session.doc.multiphysics.parts;
    let in_study: Vec<String> = match &study.model {
        Some(link) => model::bodies(&state, parts, link)?.iter().map(|b| b.name.clone()).collect(),
        None => vec![],
    };
    let report = solidcraft_multiphysics::report::build(&study, parts, &in_study, &input.run_ids).map_err(|e| bad("run_ids", e.to_string()))?;
    // A picture of each run's parts by verdict, beside a written report, where a renderer is installed.
    let pictures: Vec<Option<(String, Vec<u8>)>> = match (&input.directory, crate::capture::installed()) {
        (Some(_), Some(renderer)) => report
            .tables
            .iter()
            .enumerate()
            .map(|(i, table)| verdict_picture(session, table, &*renderer).map(|png| (format!("run-{}.png", i + 1), png)))
            .collect(),
        _ => vec![None; report.tables.len()],
    };
    let names: Vec<Option<String>> = pictures.iter().map(|p| p.as_ref().map(|(name, _)| name.clone())).collect();
    let (csv, html) = (solidcraft_multiphysics::report::to_csv(&report), solidcraft_multiphysics::report::to_html(&report, &names));
    let files = match &input.directory {
        Some(directory) => {
            let pictures: Vec<(String, Vec<u8>)> = pictures.into_iter().flatten().collect();
            solidcraft_io::multiphysics::write_report(std::path::Path::new(directory), &csv, &html, &pictures)
                .map_err(|e| bad("directory", e.to_string()))?
        }
        None => Vec::new(),
    };
    let (csv, html) = if input.inline { (Some(csv), Some(html)) } else { (None, None) };
    reply(ReportReply { runs: report.runs, files, csv, html })
}
/// Keeps a result view in a study by name (one undo step; the study's revision is unchanged, so
/// its runs stay current).
fn view_save(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ViewSave = parse(params)?;
    let mut view = input.view;
    if input.ghosted_now {
        view.ghosted = session.ghosted_bodies.iter().cloned().collect();
    }
    let mut project = session.doc.multiphysics.as_ref().clone();
    let study = project.studies.iter_mut().find(|s| s.id == input.id).ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    let replaced = solidcraft_multiphysics::result_views::save(study, view).map_err(|e| bad("view", e.to_string()))?;
    let views = study.views.clone();
    session.doc_mut().multiphysics = Arc::new(project);
    reply(ViewsReply { views, replaced: Some(replaced) })
}
/// Puts the session where a saved view was: camera, ghosted parts, colour-by.
fn view_apply(session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::setup_check::ColourBy;
    let input: ViewName = parse(params)?;
    let study = session.doc.multiphysics.study(&input.id).ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    let view = solidcraft_multiphysics::result_views::find(study, &input.name).map_err(|e| bad("name", e.to_string()))?.clone();
    let of_model = study.model.is_some();
    let state = session.world_state()?;
    let names: std::collections::HashSet<&str> = state.bodies.iter().map(|b| b.name.as_str()).collect();
    let (kept, missing): (Vec<String>, Vec<String>) = view.ghosted.iter().cloned().partition(|g| names.contains(g.as_str()));
    let mut notes = Vec::new();
    if let Some(c) = &view.camera {
        let up = crate::view::home_camera(session).map(|h| h.up).unwrap_or_default();
        let camera = solidcraft_viewport::Camera {
            target: solidcore_record::Vec3::new(c.target[0], c.target[1], c.target[2]),
            yaw: c.yaw,
            pitch: c.pitch,
            distance: c.distance,
            fov: c.fov,
            roll: c.roll,
            up,
        };
        session.camera = Some(camera);
    }
    session.ghosted_bodies = kept.into_iter().collect();
    let ghosted = state.bodies.iter().filter(|b| session.ghosted_bodies.contains(&b.name)).map(|b| b.name.clone()).collect();
    let colouring = if view.colour_by == ColourBy::None {
        session.part_swatches.clear();
        None
    } else if !of_model {
        session.part_swatches.clear();
        notes.push("colour-by needs a study of the model; the parts keep their own looks".into());
        None
    } else {
        match parts::colour_parts(session, &input.id, view.colour_by, None) {
            Ok(colouring) => Some(colouring),
            Err(e) => {
                session.part_swatches.clear();
                notes.push(format!("the parts keep their own looks: {e}"));
                None
            }
        }
    };
    if let Some(plot) = view.plot.as_ref().filter(|p| p.part.as_ref().is_some_and(|part| !names.contains(part.as_str()))) {
        notes.push(format!("the cut plot's part {} is no longer in the design", plot.part.as_deref().unwrap_or_default()));
    }
    session.refresh_display();
    reply(ViewApplied { view, ghosted, missing, colouring, notes })
}
fn view_list(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    let study = session.doc.multiphysics.study(&input.id).ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    reply(ViewsReply { views: study.views.clone(), replaced: None })
}
fn view_remove(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ViewName = parse(params)?;
    let mut project = session.doc.multiphysics.as_ref().clone();
    let study = project.studies.iter_mut().find(|s| s.id == input.id).ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    solidcraft_multiphysics::result_views::remove(study, &input.name).map_err(|e| bad("name", e.to_string()))?;
    let views = study.views.clone();
    session.doc_mut().multiphysics = Arc::new(project);
    reply(ViewsReply { views, replaced: None })
}
/// Two runs, of one study or of two studies of the model, part by part.
/// One run with a thermal report of a study of the model (`side.run_id`, or the study's latest
/// such run), its study and its component table; errors name `field`.
fn run_table(
    session: &Session,
    side: &CompareSide,
    field: &str,
) -> Result<(solidcraft_multiphysics::Study, solidcraft_multiphysics::RunEvidence, solidcraft_multiphysics::component_table::ComponentTable)> {
    let state = session.world_state()?;
    let parts = &session.doc.multiphysics.parts;
    let study = session.doc.multiphysics.study(&side.id).ok_or_else(|| bad(field, "Multiphysics study not found in this document"))?;
    let Some(link) = &study.model else { return Err(bad(field, "a study of the model has a component table; this study is not one")) };
    let run = match &side.run_id {
        Some(id) => study.runs.iter().find(|r| &r.id == id),
        None => study.runs.iter().rev().find(|r| r.artifacts.iter().any(|a| a.name == "thermal_report.json")),
    };
    let Some(run) = run.filter(|r| r.artifacts.iter().any(|a| a.name == "thermal_report.json")) else {
        return Err(bad(field, format!("study {} has no such run with a thermal report", study.name)));
    };
    let in_study: Vec<String> = model::bodies(&state, parts, link)?.iter().map(|b| b.name.clone()).collect();
    let table =
        solidcraft_multiphysics::component_table::build(study, parts, &in_study, &in_study, Some(&run.id)).map_err(|e| bad(field, e.to_string()))?;
    Ok((study.clone(), run.clone(), table))
}

/// Two runs, of one study or of two studies of the model, part by part.
fn compare_runs(session: &mut Session, params: &Value) -> Result<Value> {
    let input: CompareRuns = parse(params)?;
    let (_, _, first) = run_table(session, &input.first, "first")?;
    let (_, _, second) = run_table(session, &input.second, "second")?;
    let comparison = solidcraft_multiphysics::variants::compare(&first, &second);
    let csv = input.csv.then(|| solidcraft_multiphysics::variants::to_csv(&comparison));
    reply(CompareReply { comparison, csv })
}

/// The same study on three grids: per part, the grid convergence index of its temperature rise.
fn convergence(session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::convergence;
    let input: Convergence = parse(params)?;
    let fine = run_table(session, &input.fine, "fine")?;
    let medium = run_table(session, &input.medium, "medium")?;
    let coarse = run_table(session, &input.coarse, "coarse")?;
    let cells_of = |run: &solidcraft_multiphysics::RunEvidence, field: &str| {
        convergence::grid_of(&run.input_case).map(|(_, cells)| cells).ok_or_else(|| bad(field, "the run's case has no grid to count"))
    };
    let cells = [cells_of(&fine.1, "fine")?, cells_of(&medium.1, "medium")?, cells_of(&coarse.1, "coarse")?];
    if !(cells[0] > cells[1] && cells[1] > cells[2]) {
        return Err(bad("fine", "give the runs finest grid first: each needs more cells than the next"));
    }
    let ratios = convergence::ratios(cells).ok_or_else(|| bad("fine", "the grids' cell counts are not usable"))?;
    let inflow_k = fine
        .1
        .input_case
        .pointer("/thermal/inflow_temperature_k")
        .and_then(Value::as_f64)
        .ok_or_else(|| bad("fine", "the fine run has no heat transfer: the mesh study compares the parts' temperatures"))?;
    let mut notes = Vec::new();
    for (pair, ratio) in ["fine to medium", "medium to coarse"].iter().zip(ratios) {
        if ratio < convergence::MIN_RATIO {
            notes.push(format!("the {pair} grids differ by {ratio:.2} in cell size, under the {} the method asks for", convergence::MIN_RATIO));
        }
    }
    for (name, (study, run, _)) in [("fine", &fine), ("medium", &medium), ("coarse", &coarse)] {
        if run.study_revision != study.revision {
            notes.push(format!("the {name} run was made from an earlier setup of {}", study.name));
        }
    }
    let parts = convergence::parts([&fine.2, &medium.2, &coarse.2], cells, inflow_k);
    reply(ConvergenceReply { cells, ratios, inflow_c: inflow_k - 273.15, parts, notes })
}
/// One run's parts in their verdict colours (the study's parts; the others see-through), iso and
/// fitted, as a PNG; none for a design too large to picture or a failed drawing.
fn verdict_picture(
    session: &Session,
    table: &solidcraft_multiphysics::component_table::ComponentTable,
    renderer: &crate::capture::Renderer,
) -> Option<Vec<u8>> {
    use solidcraft_geom::{Aabb3, MeshGeom};
    use solidcraft_viewport::{Camera, Rgb, Scene, SceneMesh, StandardView};
    let world = session.world_state().ok()?;
    let verdicts: std::collections::HashMap<String, Option<solidcraft_multiphysics::results::ComponentStatus>> =
        solidcraft_multiphysics::report::verdicts(table).into_iter().collect();
    let mut scene = Scene::default();
    let mut bounds = Aabb3::EMPTY;
    let mut triangles = 0usize;
    for body in world.bodies.iter().filter(|b| session.body_visible(&b.name)) {
        let mesh = body.mesh();
        triangles = triangles.saturating_add(mesh.triangles.len());
        if triangles > 400_000 {
            return None;
        }
        bounds = bounds.union(&mesh.bounds());
        // A part with a verdict is solid in its colour; one without (no limit, no result, not
        // in the study) is see-through, so an enclosure does not hide the parts inside it.
        let (color, opacity) = match verdicts.get(&body.name) {
            Some(status) => {
                let [r, g, b] = solidcraft_multiphysics::report::verdict_rgb(*status);
                let judged = matches!(status, Some(s) if *s != solidcraft_multiphysics::results::ComponentStatus::NoLimit);
                (Rgb(r, g, b), if judged { 1.0 } else { 0.25 })
            }
            None => (crate::view::colors::BODY, 0.25),
        };
        scene.meshes.push(SceneMesh {
            mesh,
            color,
            face_colors: Vec::new(),
            opacity,
            closed: body.is_closed(),
            hatch: Default::default(),
            values: None,
        });
    }
    if bounds.is_empty() {
        return None;
    }
    let mut camera = Camera::default();
    camera.set_view(StandardView::Iso);
    camera.fit(&solidcraft_render::view_bounds(&bounds));
    scene.radius = (bounds.diagonal() * 0.5).max(camera.half_height()).max(10.);
    scene.background = Some((crate::view::colors::BG_TOP, crate::view::colors::BG_BOTTOM));
    crate::capture::png_with(Some(renderer), &scene, &camera, 960, 640).ok().map(|png| png.bytes)
}
fn query_result(session: &mut Session, params: &Value) -> Result<Value> {
    let input: solidcraft_multiphysics::results::QueryRequest = parse(params)?;
    let answer = solidcraft_multiphysics::results::query(&session.doc.multiphysics, &input).map_err(|error| bad("query", error.to_string()))?;
    reply(answer)
}
fn compare_studies(session: &mut Session, params: &Value) -> Result<Value> {
    let input: CompareStudies = parse(params)?;
    let table =
        solidcraft_multiphysics::results::compare_studies(&session.doc.multiphysics, &input.ids).map_err(|error| bad("ids", error.to_string()))?;
    reply(table)
}
fn inspect_result(session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::vtk::{self, Association};
    let input: InspectResult = parse(params)?;
    solidcraft_multiphysics::validate_artifact_name(&input.artifact).map_err(|error| bad("artifact", error.to_string()))?;
    let study = session.doc.multiphysics.study(&input.id).ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    let run = study.runs.iter().find(|run| run.id == input.run_id).ok_or_else(|| bad("run_id", "Select a retained result"))?;
    let artifact = run
        .artifacts
        .iter()
        .find(|artifact| artifact.name == input.artifact)
        .ok_or_else(|| bad("artifact", "Artifact is not retained in this run"))?;
    let content = if input.artifact.ends_with(".vtu") {
        let grid = vtk::parse_vtu(&artifact.text).map_err(|error| bad("artifact", error.to_string()))?;
        let cell_types = grid.cell_types.iter().copied().collect::<std::collections::BTreeSet<_>>().into_iter().collect();
        ResultInspectionContent::Vtu {
            points: grid.points.len(),
            cells: grid.offsets.len(),
            connectivity_indices: grid.connectivity.len(),
            cell_types,
            fields: grid
                .fields
                .into_iter()
                .map(|field| ResultField {
                    name: field.name,
                    association: match field.association {
                        Association::Point => "point",
                        Association::Cell => "cell",
                    }
                    .into(),
                    data_type: field.data_type,
                    components: field.components,
                    tuples: field.values.len() / field.components,
                })
                .collect(),
        }
    } else if input.artifact.ends_with(".pvd") {
        let frames = vtk::parse_pvd(&artifact.text)
            .map_err(|error| bad("artifact", error.to_string()))?
            .into_iter()
            .map(|frame| {
                let retained = run.artifacts.iter().any(|item| item.name == frame.file);
                ResultFrame { coordinate: frame.coordinate, file: frame.file, part: frame.part, retained }
            })
            .collect();
        ResultInspectionContent::Pvd { frames }
    } else {
        return Err(bad("artifact", "Choose a retained ASCII VTU or PVD artifact"));
    };
    let hash = solidcraft_multiphysics::case_hash(&study.case).map_err(|error| bad("id", error.to_string()))?;
    reply(ResultInspection {
        id: input.id,
        run_id: input.run_id,
        artifact: input.artifact,
        captured_study_revision: run.study_revision,
        current_study_revision: study.revision,
        case_hash: run.case_hash.clone(),
        stale: run.study_revision != study.revision || run.case_hash != hash,
        content,
    })
}

fn import_components(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ComponentImport = parse(params)?;
    let quality = match input.quality.as_deref() {
        None => solidcraft_io::MeshQuality::default(),
        Some(name) => solidcraft_io::MeshQuality::from_name(name).ok_or_else(|| bad("quality", "quality must be standard, fine or finest"))?,
    };
    let state = session.world_state()?;
    let painted = session.doc.painted(&state);
    let files = solidcraft_io::export_obj(&painted, &input.bodies, "components.mtl", quality).map_err(|e| bad("bodies", e.to_string()))?;
    reply(solidcraft_multiphysics::case_scene::import_components(&input.case, &files.obj, &input.options).map_err(|e| bad("", e.to_string()))?)
}

fn cad_bodies(session: &mut Session, _: &Value) -> Result<Value> {
    let state = session.world_state()?;
    reply(CadBodies { bodies: state.bodies.iter().map(|body| body.name.clone()).collect() })
}

fn prepare_feature(_: &mut Session, params: &Value) -> Result<Value> {
    let input: FeaturePrepare = parse(params)?;
    reply(
        solidcraft_multiphysics::authoring::prepare_feature_with_options(&input.case, input.kind, &input.name, input.replace_fan_patch_conditions)
            .map_err(|e| bad("", e.to_string()))?,
    )
}
