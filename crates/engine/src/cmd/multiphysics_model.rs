//! Studies of the loaded model: the case is built from the design's bodies and faces.
//!
//! Exact geometry decides everything about the model here: the bodies' bounding box and volume,
//! which face a picked point means and each face's plane. The tessellation that goes into the
//! case is the solver's input only (`solidcraft_multiphysics::model`).
use super::{parse, reply, selected};
use crate::schema::multiphysics::*;
use crate::{Result, Session, params::bad};
use serde_json::Value;
use solidcraft_multiphysics::model::{self, FaceRef, ModelInput, ModelLink, Opening, ResolvedFace};
use std::sync::Arc;

/// The template a model study starts from when none is named.
const DEFAULT_TEMPLATE: &str = "flow-channel";

/// The bodies of the design a link uses, by name, each one existing and solid.
fn bodies<'a>(state: &'a solidcraft_doc::ModelState, link: &ModelLink) -> Result<Vec<&'a solidcraft_doc::ModelBody>> {
    let chosen: Vec<&solidcraft_doc::ModelBody> = if link.bodies.is_empty() {
        state.bodies.iter().collect()
    } else {
        link.bodies
            .iter()
            .map(|name| state.bodies.iter().find(|b| &b.name == name).ok_or_else(|| bad("bodies", format!("the design has no body called {name}"))))
            .collect::<Result<_>>()?
    };
    if chosen.is_empty() {
        return Err(bad("bodies", "the design has no bodies to study"));
    }
    if let Some(mesh) = chosen.iter().find(|b| b.body.is_mesh()) {
        return Err(bad("bodies", format!("{} is a mesh body; a model study needs solid bodies", mesh.name)));
    }
    Ok(chosen)
}

/// A face of a body by its persistent name.
fn face_index(body: &solidcraft_doc::ModelBody, name: &str) -> Option<usize> {
    solidcraft_doc::naming::face_names(body).iter().position(|n| n == name)
}

/// Read what a build needs from the design as it is now.
fn input(session: &Session, link: &ModelLink) -> Result<ModelInput> {
    let state = session.world_state()?;
    let chosen = bodies(&state, link)?;
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    let mut volume = 0.0;
    for body in &chosen {
        let bounds = body.body.exact_bounds()?;
        for (a, (lo, hi)) in [(bounds.min.x, bounds.max.x), (bounds.min.y, bounds.max.y), (bounds.min.z, bounds.max.z)].into_iter().enumerate() {
            low[a] = low[a].min(lo);
            high[a] = high[a].max(hi);
        }
        volume += body.measure()?.volume;
    }
    let names: Vec<String> = chosen.iter().map(|b| b.name.clone()).collect();
    let painted = session.doc.painted(&state);
    let files = solidcraft_io::export_obj(&painted, &names, "components.mtl", solidcraft_io::MeshQuality::default())
        .map_err(|e| bad("bodies", e.to_string()))?;
    // Which exact face each triangle of a body lies on, read once per body that needs it.
    let mut labels: std::collections::BTreeMap<String, Vec<Option<usize>>> = std::collections::BTreeMap::new();
    let refs = link.face_refs();
    let mut faces = Vec::with_capacity(refs.len());
    for face_ref in refs {
        let body =
            chosen.iter().find(|b| b.name == face_ref.body).ok_or_else(|| bad("face", format!("the body {} is not in the study", face_ref.body)))?;
        let index = face_index(body, &face_ref.face)
            .ok_or_else(|| bad("face", format!("{} has no face called {} any more; pick the face again", body.name, face_ref.face)))?;
        let plane = match solidcraft_kernel::face_surface(&body.body, index) {
            Some(solidcraft_kernel::FaceSurface::Plane { normal, point }) => Some(([normal.x, normal.y, normal.z], [point.x, point.y, point.z])),
            _ => None,
        };
        if !labels.contains_key(&body.name) {
            // Which face each of the solver's triangles is on: the exact face nearest its centroid.
            let points: Vec<solidcore_record::Vec3> = model::triangle_centroids(&files.obj, &body.name)
                .map_err(|e| bad("bodies", e.to_string()))?
                .into_iter()
                .map(|p| solidcore_record::Vec3::new(p[0], p[1], p[2]))
                .collect();
            let nearest = body.body.nearest_faces(&points).ok_or_else(|| bad("bodies", format!("{} has no exact faces", body.name)))?;
            labels.insert(body.name.clone(), nearest);
        }
        let triangles =
            labels.get(&body.name).map(|l| l.iter().enumerate().filter(|(_, f)| **f == Some(index)).map(|(t, _)| t).collect()).unwrap_or_default();
        faces.push(ResolvedFace { name: face_ref.face.clone(), plane, body: body.name.clone(), triangles });
    }
    Ok(ModelInput { obj: files.obj, bounds_mm: [low, high], volume_mm3: volume, faces })
}

/// Build `link`'s case for the study and store both in one revision; reply with what it gave.
fn commit(session: &mut Session, id: &str, expected_revision: u64, mut link: ModelLink) -> Result<Value> {
    let study = selected(session, id)?;
    let input = input(session, &link)?;
    let built = model::build(&study.case, &link, &input).map_err(|e| bad("", e.to_string()))?;
    if let Some(issue) = built.issues.first() {
        return Err(bad("", format!("the case built from the model would not be accepted: {issue}")));
    }
    link.built_from = Some(built.digest.clone());
    let mut project = session.doc.multiphysics.as_ref().clone();
    let study = project.set_model(id, expected_revision, Some(link.clone()), built.case).map_err(|e| bad("", e.to_string()))?;
    if project != *session.doc.multiphysics {
        session.doc_mut().multiphysics = Arc::new(project);
    }
    reply(ModelReply {
        study: StudyView::try_from(&study).map_err(|e| bad("", e.to_string()))?,
        link,
        origin_mm: built.origin_mm,
        placed: built.placed,
        solid_fraction: built.solid_fraction,
        fluid_regions: built.fluid_regions,
        notes: built.notes,
        issues: built.issues,
        adapted: None,
    })
}

fn linked(session: &Session, id: &str) -> Result<(solidcraft_multiphysics::Study, ModelLink)> {
    let study = selected(session, id)?;
    let link = study.model.clone().ok_or_else(|| bad("id", "this study is not a study of the loaded model"))?;
    Ok((study, link))
}

pub(super) fn create(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelCreate = parse(params)?;
    let link = ModelLink { analysis: input.analysis, bodies: input.bodies, margin_mm: input.margin_mm, cell_mm: input.cell_mm, ..Default::default() };
    link.validate().map_err(|e| bad("", e.to_string()))?;
    let template = input.template.as_deref().unwrap_or(DEFAULT_TEMPLATE);
    let mut project = session.doc.multiphysics.as_ref().clone();
    project.next_id = project.next_id.max(session.multiphysics_jobs.next_study);
    let study = project.create(template, &input.name).map_err(|e| bad("template", e.to_string()))?;
    // Build against the new study before anything reaches the document, so a model that
    // cannot be studied leaves no half-made study behind.
    let previous = session.doc.multiphysics.clone();
    session.doc_mut().multiphysics = Arc::new(project.clone());
    let built = commit(session, &study.id, study.revision, link);
    match built {
        Ok(reply) => {
            session.multiphysics_jobs.next_study = project.next_id;
            Ok(reply)
        }
        Err(error) => {
            session.doc_mut().multiphysics = previous;
            Err(error)
        }
    }
}

pub(super) fn set(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelSet = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    if let Some(analysis) = input.analysis {
        link.analysis = analysis;
    }
    if let Some(bodies) = input.bodies {
        link.bodies = bodies;
    }
    if let Some(margin) = input.margin_mm {
        link.margin_mm = margin;
    }
    if let Some(cell) = input.cell_mm {
        link.cell_mm = (cell != 0.0).then_some(cell);
    }
    if let Some(cut) = input.cut_cells {
        link.cut_cells = cut;
    }
    if let Some(level) = input.mesh_level {
        link.mesh_level = (level != 0).then_some(level);
        link.cell_mm = None;
    }
    commit(session, &input.id, input.expected_revision, link)
}

pub(super) fn fan(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelFanInput = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    if input.remove {
        link.fan = None;
        return commit(session, &input.id, input.expected_revision, link);
    }
    let pick = input.face.ok_or_else(|| bad("face", "give the face: {body, face} or {body, point}"))?;
    let face = pick_face(session, &link, pick)?;
    let curve = input.curve.or_else(|| link.fan.as_ref().map(|f| f.curve.clone())).ok_or_else(|| bad("curve", "give the fan curve"))?;
    let name = input.name.or_else(|| link.fan.as_ref().map(|f| f.name.clone())).unwrap_or_else(|| "Fan".into());
    link.fan = Some(model::ModelFan { name, face, curve, total_pressure: input.total_pressure });
    commit(session, &input.id, input.expected_revision, link)
}

/// The name every adaptive refinement starts with.
const ADAPT: &str = "adapt_";

pub(super) fn adapt(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelAdapt = parse(params)?;
    let (study, mut link) = linked(session, &input.id)?;
    link.refinements.retain(|r| !r.name.starts_with(ADAPT));
    if input.remove {
        return commit(session, &input.id, input.expected_revision, link);
    }
    let run = study.runs.iter().find(|r| r.id == input.run_id).ok_or_else(|| bad("run_id", format!("the study has no kept run {}", input.run_id)))?;
    let artifact = match &input.artifact {
        Some(name) => run.artifacts.iter().find(|a| &a.name == name),
        None => run.artifacts.iter().filter(|a| a.name.starts_with("flow_") && a.name.ends_with(".vtu")).max_by(|a, b| a.name.cmp(&b.name)),
    }
    .ok_or_else(|| bad("artifact", "the run kept no such field file (flow_*.vtu)"))?;
    let found = solidcraft_multiphysics::adapt::regions(
        &artifact.text,
        input.field.as_deref(),
        input.fraction.unwrap_or(0.1),
        input.regions.unwrap_or(4),
    )
    .map_err(|e| bad("", e.to_string()))?;
    // The run's coordinates are case metres from the domain's low corner.
    let (domain, _) = marks(session, &link)?;
    let origin = domain.map(|[low, _]| low).ok_or_else(|| bad("id", "the study's bodies have no extent"))?;
    let regions_mm: Vec<[[f64; 3]; 2]> =
        found.regions_m.iter().map(|(low, high)| [[0, 1, 2].map(|a| low[a] * 1000.0 + origin[a]), [0, 1, 2].map(|a| high[a] * 1000.0 + origin[a])]).collect();
    for (index, [low, high]) in regions_mm.iter().enumerate() {
        link.refinements.push(model::Refinement {
            name: format!("{ADAPT}{}", index + 1),
            target: model::RefineTarget::Region { low_mm: *low, high_mm: *high },
            factor: input.factor.unwrap_or(2),
            margin_cells: 0.0,
        });
    }
    let mut reply = commit(session, &input.id, input.expected_revision, link)?;
    let adapted = ModelAdapted {
        run_id: run.id.clone(),
        artifact: artifact.name.clone(),
        field: found.field,
        cells: found.cells,
        flagged: found.flagged,
        regions_mm,
    };
    reply["adapted"] = serde_json::to_value(adapted).map_err(|e| bad("", e.to_string()))?;
    Ok(reply)
}

/// The persistent name of a picked face (by its name or by a point on it, exact nearest face),
/// on one of the study's bodies.
fn pick_face(session: &Session, link: &ModelLink, pick: ModelFacePick) -> Result<FaceRef> {
    let state = session.world_state()?;
    let body =
        state.bodies.iter().find(|b| b.name == pick.body).ok_or_else(|| bad("face", format!("the design has no body called {}", pick.body)))?;
    if !link.bodies.is_empty() && !link.bodies.contains(&body.name) {
        return Err(bad("face", format!("{} is not one of the study's bodies", body.name)));
    }
    let face = match (pick.face, pick.point) {
        (Some(name), None) => face_index(body, &name).map(|_| name).ok_or_else(|| bad("face", format!("{} has no face called that", body.name)))?,
        (None, Some(p)) => {
            if p.iter().any(|v| !v.is_finite()) {
                return Err(bad("face", "the point must be finite"));
            }
            let (distance, index) = body
                .body
                .nearest_face(solidcore_record::Vec3::new(p[0], p[1], p[2]))
                .ok_or_else(|| bad("face", format!("{} has no faces to pick", body.name)))?;
            let reach = (body.body.size() * 0.01).max(0.01);
            if distance > reach {
                return Err(bad("face", format!("the point is {distance:.3} mm from {}'s nearest face; pick a point on the face", body.name)));
            }
            solidcraft_doc::naming::face_names(body).get(index).cloned().ok_or_else(|| bad("face", "the picked face has no name"))?
        }
        _ => return Err(bad("face", "give the face by name or by a point, not both")),
    };
    Ok(FaceRef { body: body.name.clone(), face })
}

pub(super) fn condition(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelCondition = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    if input.remove {
        let before = link.openings.len();
        link.openings.retain(|opening| opening.name != input.name);
        if link.openings.len() == before {
            return Err(bad("name", format!("there is no condition called {}", input.name)));
        }
        return commit(session, &input.id, input.expected_revision, link);
    }
    let role = input.role.ok_or_else(|| bad("role", "give the condition's role, or remove: true"))?;
    let pick = input.face.ok_or_else(|| bad("face", "give the face: {body, face} or {body, point}"))?;
    let face = pick_face(session, &link, pick)?;
    let opening = Opening { name: input.name.clone(), face, role };
    match link.openings.iter_mut().find(|o| o.name == input.name) {
        Some(existing) => *existing = opening,
        None => link.openings.push(opening),
    }
    commit(session, &input.id, input.expected_revision, link)
}

pub(super) fn refine(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelRefine = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    if input.remove {
        let before = link.refinements.len();
        link.refinements.retain(|r| r.name != input.name);
        if link.refinements.len() == before {
            return Err(bad("name", format!("there is no refinement called {}", input.name)));
        }
        return commit(session, &input.id, input.expected_revision, link);
    }
    let target = match (input.body, input.face) {
        (Some(body), None) => {
            let state = session.world_state()?;
            if !state.bodies.iter().any(|b| b.name == body) || (!link.bodies.is_empty() && !link.bodies.contains(&body)) {
                return Err(bad("body", format!("{body} is not one of the study's bodies")));
            }
            model::RefineTarget::Body { body }
        }
        (None, Some(pick)) => model::RefineTarget::Face { face: pick_face(session, &link, pick)? },
        _ => return Err(bad("body", "refine around a body or a face, one of them")),
    };
    let refinement =
        model::Refinement { name: input.name.clone(), target, factor: input.factor.unwrap_or(2), margin_cells: input.margin_cells.unwrap_or(1.0) };
    match link.refinements.iter_mut().find(|r| r.name == input.name) {
        Some(existing) => *existing = refinement,
        None => link.refinements.push(refinement),
    }
    commit(session, &input.id, input.expected_revision, link)
}

pub(super) fn goal(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelGoalInput = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    if input.remove {
        let before = link.goals.len();
        link.goals.retain(|g| g.name != input.name);
        if link.goals.len() == before {
            return Err(bad("name", format!("there is no goal called {}", input.name)));
        }
        return commit(session, &input.id, input.expected_revision, link);
    }
    let kind = input.kind.ok_or_else(|| bad("kind", "give the goal's kind, or remove: true"))?;
    let face = input.face.map(|pick| pick_face(session, &link, pick)).transpose()?;
    if let Some(body) = &input.body {
        known_body(session, &link, body)?;
    }
    let goal = model::ModelGoal { name: input.name.clone(), kind, face, point_mm: input.point_mm, expression: input.expression, body: input.body };
    match link.goals.iter_mut().find(|g| g.name == input.name) {
        Some(existing) => *existing = goal,
        None => link.goals.push(goal),
    }
    commit(session, &input.id, input.expected_revision, link)
}

/// That `body` is one of the study's bodies in the design now.
fn known_body(session: &Session, link: &ModelLink, body: &str) -> Result<()> {
    let state = session.world_state()?;
    if !state.bodies.iter().any(|b| b.name == body) || (!link.bodies.is_empty() && !link.bodies.iter().any(|b| b == body)) {
        return Err(bad("body", format!("{body} is not one of the study's bodies")));
    }
    Ok(())
}

/// The study's heat settings, which must be on.
fn heat_of(link: &mut ModelLink) -> Result<&mut model::ModelHeat> {
    link.heat.as_mut().ok_or_else(|| bad("id", "turn heat transfer on first (multiphysics.model.heat with enabled: true)"))
}

pub(super) fn heat(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelHeatSet = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    match input.enabled {
        Some(false) => {
            link.heat = None;
            // Temperature goals go with it.
            link.goals.retain(|g| !g.kind.thermal());
        }
        Some(true) if link.heat.is_none() => link.heat = Some(model::ModelHeat::default()),
        _ => {}
    }
    if input.inflow_temperature_k.is_some() || input.fluid.is_some() {
        let heat = heat_of(&mut link)?;
        if let Some(t) = input.inflow_temperature_k {
            heat.inflow_temperature_k = t;
        }
        if let Some(fluid) = input.fluid {
            heat.fluid = fluid;
        }
    }
    commit(session, &input.id, input.expected_revision, link)
}

pub(super) fn material(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelMaterial = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    if input.bodies.is_empty() || input.bodies.len() > model::MAX_BODIES {
        return Err(bad("bodies", format!("name 1 to {} bodies", model::MAX_BODIES)));
    }
    for body in &input.bodies {
        known_body(session, &link, body)?;
    }
    let heat = heat_of(&mut link)?;
    if input.remove {
        heat.materials.retain(|m| !input.bodies.contains(&m.body));
    } else {
        for body in &input.bodies {
            let index = match heat.materials.iter().position(|m| &m.body == body) {
                Some(index) => index,
                None => {
                    let conductivity = input.conductivity_w_per_m_k.ok_or_else(|| bad("conductivity_w_per_m_k", format!("give {body}'s conductivity")))?;
                    heat.materials.push(model::BodyHeat {
                        body: body.clone(),
                        conductivity_w_per_m_k: conductivity,
                        power_w: 0.0,
                        density_kg_per_m3: None,
                        specific_heat_j_per_kg_k: None,
                        max_temperature_k: None,
                    });
                    heat.materials.len() - 1
                }
            };
            let Some(m) = heat.materials.get_mut(index) else { continue };
            if let Some(k) = input.conductivity_w_per_m_k {
                m.conductivity_w_per_m_k = k;
            }
            if let Some(p) = input.power_w {
                m.power_w = p;
            }
            if input.density_kg_per_m3.is_some() {
                m.density_kg_per_m3 = input.density_kg_per_m3;
            }
            if input.specific_heat_j_per_kg_k.is_some() {
                m.specific_heat_j_per_kg_k = input.specific_heat_j_per_kg_k;
            }
            if input.max_temperature_k.is_some() {
                m.max_temperature_k = input.max_temperature_k;
            }
        }
    }
    commit(session, &input.id, input.expected_revision, link)
}

pub(super) fn heat_source(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelHeatSource = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    if input.remove {
        let heat = heat_of(&mut link)?;
        let before = heat.face_sources.len();
        heat.face_sources.retain(|s| s.name != input.name);
        if heat.face_sources.len() == before {
            return Err(bad("name", format!("there is no heat source called {}", input.name)));
        }
        return commit(session, &input.id, input.expected_revision, link);
    }
    let pick = input.face.ok_or_else(|| bad("face", "give the face: {body, face} or {body, point}"))?;
    let face = pick_face(session, &link, pick)?;
    let heat = heat_of(&mut link)?;
    let power_w = match (input.power_w, heat.face_sources.iter().find(|s| s.name == input.name)) {
        (Some(p), _) => p,
        (None, Some(existing)) => existing.power_w,
        (None, None) => return Err(bad("power_w", "give the heat source's power")),
    };
    let source = model::FaceHeat { name: input.name.clone(), face, power_w };
    match heat.face_sources.iter_mut().find(|s| s.name == input.name) {
        Some(existing) => *existing = source,
        None => heat.face_sources.push(source),
    }
    commit(session, &input.id, input.expected_revision, link)
}

pub(super) fn thermal_face(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelThermalFace = parse(params)?;
    let (_, mut link) = linked(session, &input.id)?;
    if input.remove {
        let heat = heat_of(&mut link)?;
        let before = heat.face_conditions.len();
        heat.face_conditions.retain(|c| c.name != input.name);
        if heat.face_conditions.len() == before {
            return Err(bad("name", format!("there is no thermal condition called {}", input.name)));
        }
        return commit(session, &input.id, input.expected_revision, link);
    }
    let condition = input.condition.ok_or_else(|| bad("condition", "give the thermal condition, or remove: true"))?;
    let pick = input.face.ok_or_else(|| bad("face", "give the face: {body, face} or {body, point}"))?;
    let face = pick_face(session, &link, pick)?;
    let heat = heat_of(&mut link)?;
    let thermal = model::FaceThermal { name: input.name.clone(), face, condition };
    match heat.face_conditions.iter_mut().find(|c| c.name == input.name) {
        Some(existing) => *existing = thermal,
        None => heat.face_conditions.push(thermal),
    }
    commit(session, &input.id, input.expected_revision, link)
}

pub(super) fn build(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    let (study, link) = linked(session, &input.id)?;
    commit(session, &input.id, study.revision, link)
}

/// Whether the study's case was built from the model as it is now; `Err` text when it cannot
/// even be read (a body or face is gone).
pub(super) fn current(session: &Session, id: &str) -> Result<(ModelLink, std::result::Result<bool, String>)> {
    let (_, link) = linked(session, id)?;
    let now = input(session, &link).and_then(|input| input.digest(&link).map_err(|e| bad("", e.to_string())));
    Ok(match now {
        Ok(digest) => (link.clone(), Ok(link.built_from.as_deref() == Some(digest.as_str()))),
        Err(error) => (link, Err(error.to_string())),
    })
}

/// The domain box in CAD mm and a mark on each condition's face, as the design is now.
fn marks(session: &Session, link: &ModelLink) -> Result<(Option<[[f64; 3]; 2]>, Vec<ModelFaceMark>)> {
    let state = session.world_state()?;
    let chosen = bodies(&state, link)?;
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for body in &chosen {
        let b = body.body.exact_bounds()?;
        for (a, (lo, hi)) in [(b.min.x, b.max.x), (b.min.y, b.max.y), (b.min.z, b.max.z)].into_iter().enumerate() {
            low[a] = low[a].min(lo);
            high[a] = high[a].max(hi);
        }
    }
    let margin = if link.analysis == model::Analysis::External { link.margin_mm } else { 0.0 };
    let domain = (0..3).all(|a| low[a].is_finite() && high[a].is_finite()).then(|| [low.map(|v| v - margin), high.map(|v| v + margin)]);
    let faces = link
        .openings
        .iter()
        .map(|opening| {
            let found = chosen.iter().find(|b| b.name == opening.face.body).and_then(|b| face_index(b, &opening.face.face).map(|i| (b, i)));
            let point = found.and_then(|(b, i)| solidcraft_doc::naming::face_spot(b, i)).map(|spot| [spot.point.x, spot.point.y, spot.point.z]);
            let normal = found.and_then(|(b, i)| match solidcraft_kernel::face_surface(&b.body, i) {
                Some(solidcraft_kernel::FaceSurface::Plane { normal, .. }) => Some([normal.x, normal.y, normal.z]),
                _ => None,
            });
            ModelFaceMark { name: opening.name.clone(), point_mm: point, normal }
        })
        .collect();
    Ok((domain, faces))
}

pub(super) fn status(session: &mut Session, params: &Value) -> Result<Value> {
    let input: StudyId = parse(params)?;
    let (link, now) = current(session, &input.id)?;
    let (domain_mm, faces) = marks(session, &link).unwrap_or((None, Vec::new()));
    reply(ModelStatus { id: input.id, current: now.as_ref().is_ok_and(|current| *current), problem: now.err(), link, domain_mm, faces })
}

/// A point half a cell off a picked face, on its fluid side: the face's exact outward side,
/// found by classifying points just off it along the axes (no mesh decides it).
fn off_face(session: &Session, pick: &ModelFacePick, step_mm: f64) -> Result<[f64; 3]> {
    let state = session.world_state()?;
    let body =
        state.bodies.iter().find(|b| b.name == pick.body).ok_or_else(|| bad("face", format!("the design has no body called {}", pick.body)))?;
    let p = pick.point.ok_or_else(|| bad("face", "give a point on the face"))?;
    if p.iter().any(|v| !v.is_finite()) {
        return Err(bad("face", "the point must be finite"));
    }
    let probe = (body.body.size() * 1e-4).max(1e-4);
    let outside =
        |q: [f64; 3]| matches!(body.body.classify_point(solidcore_record::Vec3::new(q[0], q[1], q[2])), Ok(solidcraft_kernel::PointClass::Outside));
    let mut direction = None;
    for (axis, sign) in [(0, 1.0), (0, -1.0), (1, 1.0), (1, -1.0), (2, 1.0), (2, -1.0)] {
        let mut ahead = p;
        let mut behind = p;
        ahead[axis] += sign * probe;
        behind[axis] -= sign * probe;
        if outside(ahead) && !outside(behind) {
            direction = Some((axis, sign));
            break;
        }
    }
    let (axis, sign) = direction.ok_or_else(|| bad("face", format!("the point is not on a face of {}", body.name)))?;
    let mut q = p;
    q[axis] += sign * step_mm;
    Ok(q)
}

pub(super) fn leak(session: &mut Session, params: &Value) -> Result<Value> {
    let input: ModelLeak = parse(params)?;
    let (study, link) = linked(session, &input.id)?;
    let (domain, _) = marks(session, &link)?;
    let [origin, _] = domain.ok_or_else(|| bad("id", "the study's bodies have no extent"))?;
    let (dimensions, lengths) = solidcraft_multiphysics::case_scene::grid(&study.case).map_err(|e| bad("", e.to_string()))?;
    let half_cell_mm = (0..3).map(|a| lengths[a] * 1000.0 / dimensions[a] as f64).fold(f64::INFINITY, f64::min) * 0.5;
    let (from, to) = (off_face(session, &input.from, half_cell_mm)?, off_face(session, &input.to, half_cell_mm)?);
    let to_case = |p: [f64; 3]| [0, 1, 2].map(|a| (p[a] - origin[a]) / 1000.0);
    let path = solidcraft_multiphysics::fluid_volume::fluid_path(&study.case, to_case(from), to_case(to)).map_err(|e| bad("", e.to_string()))?;
    let path_mm: Vec<[f64; 3]> = path.iter().flatten().map(|p| [0, 1, 2].map(|a| p[a] * 1000.0 + origin[a])).collect();
    reply(ModelLeakReply { connected: path.is_some(), path_mm })
}
