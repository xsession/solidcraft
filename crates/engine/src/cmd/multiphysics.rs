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
    CommandSpec::new("multiphysics.result.inspect", "Inspect Multiphysics Result", crate::checked!("multiphysics.result.inspect", inspect_result))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<InspectResult>::new()).output(output_schema_for::<ResultInspection>)
        .params("id, run_id, artifact: one retained ASCII VTU/PVD file. Returns bounded geometry, field or collection metadata, without dataset bytes or inferred units."),
    CommandSpec::new("multiphysics.result.query", "Query Multiphysics Result", crate::checked!("multiphysics.result.query", query_result))
        .doc_only().noundo().icon("inspect").contract(&TypedContract::<solidcraft_multiphysics::results::QueryRequest>::new()).output(output_schema_for::<solidcraft_multiphysics::results::QueryResponse>)
        .params("id, run_id, query: {kind: probe|min_max|line|integrate|compare|goals|components, ...}. Probes interpolate inside cells and report the nearest cell outside; line samples outside the mesh have no value; integrals are measure weighted; compare needs identical meshes. Every answer carries the run's case hash and a stale flag."),
    CommandSpec::new("multiphysics.parts.list", "List Parts and Selections", crate::checked!("multiphysics.parts.list", parts::list))
        .doc_only().noundo().icon("list").contract(&TypedContract::<Empty>::new()).output(output_schema_for::<PartsReply>)
        .params("No parameters. Every body of the design with its tags and attributes, every part selection with the bodies it holds now, and recorded parts the design no longer has."),
    CommandSpec::new("multiphysics.parts.table", "Component Table", crate::checked!("multiphysics.parts.table", parts::table))
        .doc_only().noundo().icon("list").contract(&TypedContract::<PartsTable>::new()).output(output_schema_for::<PartsTableReply>)
        .params("id (a study); run_id? (default: the latest run with a thermal report); selection? (only its bodies); all_bodies? (also the design's bodies the study leaves out); csv? (add the table as CSV). One row per body: tags, attributes, conductivity, power (in the body and on its faces), limit, the run's temperatures, margin and pass/warn/fail, with the study's total power and the part with the least margin."),
    CommandSpec::new("multiphysics.parts.tag", "Tag Parts", crate::checked!("multiphysics.parts.tag", parts::tag))
        .doc_only().icon("settings").contract(&TypedContract::<PartsTag>::new()).output(output_schema_for::<PartsReply>)
        .params("bodies (names) or selection (a part selection's bodies as they are now); add_tags?, remove_tags? (single words); attributes? {name: number | text | null (remove)}. volume_mm3 and size_mm are measured and cannot be set. One undo step."),
    CommandSpec::new("multiphysics.parts.selection", "Part Selection", crate::checked!("multiphysics.parts.selection", parts::selection))
        .doc_only().icon("settings").contract(&TypedContract::<PartsSelection>::new()).output(output_schema_for::<PartsReply>)
        .params("name; members: {kind: saved, bodies: [names]} | {kind: query, rules: [{rule: tag|not_tag, tag} | {rule: name, pattern (* and ?)} | {rule: attribute, key (own or volume_mm3, size_mm), compare: exists|eq|ne|lt|le|gt|ge, value?} | {rule: in_selection, name (an earlier selection)}]} (a part meets every rule); or remove: true. A query is worked out again each time it is used; model studies use one with multiphysics.model.set selection. One undo step."),
    CommandSpec::new("multiphysics.model.create", "New Study of the Model", crate::checked!("multiphysics.model.create", model::create))
        .doc_only().icon("new").contract(&TypedContract::<ModelCreate>::new()).output(output_schema_for::<ModelReply>)
        .params("name; template? (a 3-D flow template whose fluid and solver settings it starts from, default flow-channel); analysis?: internal (fluid inside the bodies; the domain is their box and its faces are walls) | external (around them, with margin_mm); bodies? (names; default every body) or selection? (a part selection, worked out again at every build); margin_mm?; cell_mm? (default from the bodies). Builds the case from the design's bodies in one undo step; nothing runs."),
    CommandSpec::new("multiphysics.model.set", "Change Study of the Model", crate::checked!("multiphysics.model.set", model::set))
        .doc_only().icon("settings").contract(&TypedContract::<ModelSet>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision; any of analysis, bodies, selection (a part selection; \"\" back to bodies), margin_mm, cell_mm (0: automatic), cut_cells (keep the fluid part of the cells the bodies cut; companion 0.14.0), mesh_level (1 coarsest to 7 finest, 3 the automatic size; 0 back to automatic; clears a typed cell_mm). Builds the case again from the design in one undo step."),
    CommandSpec::new("multiphysics.model.condition", "Condition on a Model Face", crate::checked!("multiphysics.model.condition", model::condition))
        .doc_only().icon("settings").contract(&TypedContract::<ModelCondition>::new()).output(output_schema_for::<ModelReply>)
        .params("id, expected_revision, name; face: {body, face (its name) | point [x,y,z] mm on it}; role: {kind: inlet, speed_m_per_s} | {kind: outlet, pressure_m2_per_s2} | {kind: wall} | {kind: symmetry}; or remove: true. The face must be planar and lie on the computational domain's boundary (a vent's outer face in an internal study). Builds the case again in one undo step."),
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
        .params("id, expected_revision; any of enabled (heat transfer on or off), inflow_temperature_k, fluid: {specific_heat_j_per_kg_k, conductivity_w_per_m_k} (air by default). With heat on, the bodies conduct heat (conjugate heat transfer) and temperature is solved after the flow; cut cells and inlets or outlets on faces inside the domain are not available with it yet. Builds the case again in one undo step."),
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
