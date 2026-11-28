//! The part registry: tags and attributes on the design's bodies, and named selections of them
//! that model studies point at (`solidcraft_multiphysics::parts`).
//!
//! A query on size or volume reads them from each body's exact geometry.
use super::{parse, reply};
use crate::schema::multiphysics::*;
use crate::{Result, Session, params::bad};
use serde_json::Value;
use solidcraft_multiphysics::parts::{PartFacts, Registry, Retag};
use std::sync::Arc;

/// Every body of the design as a query sees it, measured when `measure` says a query needs it.
fn facts(state: &solidcraft_doc::ModelState, measure: bool) -> Result<Vec<PartFacts>> {
    state
        .bodies
        .iter()
        .map(|body| {
            let mut fact = PartFacts { name: body.name.clone(), ..PartFacts::default() };
            if measure && !body.body.is_mesh() {
                let b = body.body.exact_bounds()?;
                fact.size_mm = Some((b.max.x - b.min.x).max(b.max.y - b.min.y).max(b.max.z - b.min.z));
                fact.volume_mm3 = Some(body.measure()?.volume);
            }
            Ok(fact)
        })
        .collect()
}

/// The bodies of the design the selection called `name` holds now, in the design's order, and
/// the bodies a saved selection names that the design no longer has.
pub(super) fn resolve(state: &solidcraft_doc::ModelState, parts: &Registry, name: &str) -> Result<(Vec<String>, Vec<String>)> {
    let facts = facts(state, parts.needs_measure(name))?;
    parts.resolve(name, &facts).map_err(|e| bad("selection", e.to_string()))
}

/// The bodies a command names, by name or by a part selection (not both), and the bodies a
/// saved selection names that the design no longer has.
fn named(state: &solidcraft_doc::ModelState, parts: &Registry, bodies: Vec<String>, selection: Option<&str>) -> Result<(Vec<String>, Vec<String>)> {
    match (selection, bodies.is_empty()) {
        (Some(_), false) => Err(bad("selection", "name bodies or a part selection, not both")),
        (Some(selection), true) => resolve(state, parts, selection),
        (None, _) => {
            if let Some(unknown) = bodies.iter().find(|n| !state.bodies.iter().any(|b| &&b.name == n)) {
                return Err(bad("bodies", format!("the design has no body called {unknown}")));
            }
            Ok((bodies, Vec::new()))
        }
    }
}

fn view(session: &Session) -> Result<PartsReply> {
    let state = session.world_state()?;
    let parts = &session.doc.multiphysics.parts;
    let measure = parts.selections.iter().any(|s| parts.needs_measure(&s.name));
    let facts = facts(&state, measure)?;
    let selections = parts
        .selections
        .iter()
        .map(|s| match parts.resolve(&s.name, &facts) {
            Ok((bodies, missing)) => PartSelectionView { selection: s.clone(), bodies, missing, problem: None },
            Err(e) => PartSelectionView { selection: s.clone(), bodies: vec![], missing: vec![], problem: Some(e.to_string()) },
        })
        .collect();
    let views = state
        .bodies
        .iter()
        .map(|b| {
            let info = parts.parts.get(&b.name);
            PartView {
                body: b.name.clone(),
                tags: info.map(|i| i.tags.iter().cloned().collect()).unwrap_or_default(),
                attributes: info.map(|i| i.attributes.clone()).unwrap_or_default(),
            }
        })
        .collect();
    let missing = parts.parts.keys().filter(|k| !state.bodies.iter().any(|b| &&b.name == k)).cloned().collect();
    Ok(PartsReply { parts: views, selections, missing })
}

/// Put `parts` in the document when it changed, and reply with the registry as it is now.
fn store(session: &mut Session, parts: Registry) -> Result<Value> {
    if parts != session.doc.multiphysics.parts {
        let mut project = session.doc.multiphysics.as_ref().clone();
        project.parts = parts;
        project.validate().map_err(|e| bad("", e.to_string()))?;
        session.doc_mut().multiphysics = Arc::new(project);
    }
    reply(view(session)?)
}

pub(super) fn list(session: &mut Session, _: &Value) -> Result<Value> {
    reply(view(session)?)
}

pub(super) fn tag(session: &mut Session, params: &Value) -> Result<Value> {
    let input: PartsTag = parse(params)?;
    let mut parts = session.doc.multiphysics.parts.clone();
    let state = session.world_state()?;
    let bodies = named(&state, &parts, input.bodies, input.selection.as_deref())?.0;
    if bodies.is_empty() {
        return Err(bad("bodies", "no bodies to tag"));
    }
    if input.add_tags.is_empty() && input.remove_tags.is_empty() && input.attributes.is_empty() {
        return Err(bad("add_tags", "give tags to add or remove, or attributes to set"));
    }
    let change = Retag { add_tags: input.add_tags, remove_tags: input.remove_tags, attributes: input.attributes };
    parts.retag(&bodies, &change).map_err(|e| bad("", e.to_string()))?;
    store(session, parts)
}

/// Tag parts named by reference designators with the kind of part their letters name.
pub(super) fn group(session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::parts::AttrValue;
    let input: PartsGroup = parse(params)?;
    let mut parts = session.doc.multiphysics.parts.clone();
    let state = session.world_state()?;
    let bodies = if input.bodies.is_empty() && input.selection.is_none() {
        state.bodies.iter().map(|b| b.name.clone()).collect()
    } else {
        named(&state, &parts, input.bodies, input.selection.as_deref())?.0
    };
    let (grouped, other) = solidcraft_multiphysics::grouping::group(&bodies);
    for part in &grouped {
        let designator = Some(AttrValue::Text(part.designator.clone()));
        let change = Retag {
            add_tags: vec![part.class.clone()],
            attributes: [(solidcraft_multiphysics::grouping::DESIGNATOR.to_string(), designator)].into_iter().collect(),
            ..Retag::default()
        };
        parts.retag(std::slice::from_ref(&part.body), &change).map_err(|e| bad("bodies", e.to_string()))?;
    }
    store(session, parts)?;
    reply(PartsGroupReply { grouped, other })
}

/// A body's shape from its exact geometry, or none for a mesh body (it has none).
fn shape(body: &solidcraft_doc::ModelBody) -> Result<Option<solidcraft_multiphysics::grouping::Shape>> {
    if body.body.is_mesh() {
        return Ok(None);
    }
    let b = body.body.exact_bounds()?;
    let sides = [b.max.x - b.min.x, b.max.y - b.min.y, b.max.z - b.min.z];
    let shape = solidcraft_multiphysics::grouping::Shape::new(body.measure()?.volume, sides).map_err(|e| bad("body", e.to_string()))?;
    Ok(Some(shape))
}

/// The parts shaped like one part, from their exact volumes and bounding boxes.
pub(super) fn similar(session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::grouping::{DEFAULT_TOLERANCE, MAX_TOLERANCE, Shape};
    let input: PartsSimilar = parse(params)?;
    let tolerance = match input.tolerance_pct {
        Some(pct) if pct.is_finite() && (0.0..=MAX_TOLERANCE * 100.0).contains(&pct) => pct / 100.0,
        Some(_) => return Err(bad("tolerance_pct", format!("give 0 to {} per cent", MAX_TOLERANCE * 100.0))),
        None => DEFAULT_TOLERANCE,
    };
    let state = session.world_state()?;
    let target =
        state.bodies.iter().find(|b| b.name == input.body).ok_or_else(|| bad("body", format!("the design has no body called {}", input.body)))?;
    let like = shape(target)?.ok_or_else(|| bad("body", "a mesh body has no exact geometry to compare"))?;
    // The bounding box first: only parts whose sides already match are measured for volume.
    let near = |a: f64, b: f64| (a - b).abs() <= tolerance * a.abs().max(b.abs()) + 1e-9;
    let mut shapes: Vec<(String, Shape)> = Vec::new();
    for body in state.bodies.iter().filter(|b| !b.body.is_mesh()) {
        let b = body.body.exact_bounds()?;
        let mut sides = [b.max.x - b.min.x, b.max.y - b.min.y, b.max.z - b.min.z];
        sides.sort_by(f64::total_cmp);
        if sides.iter().zip(like.sides_mm).all(|(a, b)| near(*a, b))
            && let Some(shape) = shape(body)?
        {
            shapes.push((body.name.clone(), shape));
        }
    }
    let similar = solidcraft_multiphysics::grouping::similar(&like, &shapes, tolerance).map_err(|e| bad("tolerance_pct", e.to_string()))?;
    if input.select {
        session.selection = similar.iter().map(|name| crate::Sel::Body { name: name.clone() }).collect();
        session.refresh_display();
    }
    reply(PartsSimilarReply { similar })
}

pub(super) fn selection(session: &mut Session, params: &Value) -> Result<Value> {
    let input: PartsSelection = parse(params)?;
    let members = match (input.remove, input.members) {
        (true, Some(_)) => return Err(bad("members", "give members or remove, not both")),
        (true, None) => {
            if let Some(study) =
                session.doc.multiphysics.studies.iter().find(|s| s.model.as_ref().and_then(|m| m.selection.as_ref()) == Some(&input.name))
            {
                return Err(bad("name", format!("{} uses this selection; give that study its own bodies first", study.name)));
            }
            None
        }
        (false, None) => return Err(bad("members", "give the selection's members")),
        (false, Some(members)) => Some(members),
    };
    let mut parts = session.doc.multiphysics.parts.clone();
    parts.set_selection(&input.name, members).map_err(|e| bad("", e.to_string()))?;
    store(session, parts)
}

pub(super) fn table(session: &mut Session, params: &Value) -> Result<Value> {
    let input: PartsTable = parse(params)?;
    let study = session.doc.multiphysics.study(&input.id).cloned().ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    let state = session.world_state()?;
    let parts = &session.doc.multiphysics.parts;
    let in_study: Vec<String> = match &study.model {
        Some(link) => super::model::bodies(&state, parts, link)?.iter().map(|b| b.name.clone()).collect(),
        None => vec![],
    };
    let only = match &input.selection {
        Some(selection) => Some(resolve(&state, parts, selection)?.0),
        None => None,
    };
    let shown: Vec<String> = state
        .bodies
        .iter()
        .map(|b| b.name.clone())
        .filter(|n| input.all_bodies || in_study.contains(n))
        .filter(|n| only.as_ref().is_none_or(|o| o.contains(n)))
        .collect();
    let table = solidcraft_multiphysics::component_table::build(&study, parts, &shown, &in_study, input.run_id.as_deref())
        .map_err(|e| bad("run_id", e.to_string()))?;
    let csv = input.csv.then(|| solidcraft_multiphysics::component_table::to_csv(&table));
    reply(PartsTableReply { table, csv })
}

/// The design's selection: the one the Multiphysics tree, the component table and the 3-D view
/// all show, so picking a part in any of them picks it in the others.
pub(super) fn select(session: &mut Session, params: &Value) -> Result<Value> {
    let input: PartsSelect = parse(params)?;
    let state = session.world_state()?;
    let (bodies, missing) = named(&state, &session.doc.multiphysics.parts, input.bodies, input.selection.as_deref())?;
    if !input.add {
        session.selection.clear();
    }
    for name in bodies {
        let pick = crate::Sel::Body { name };
        if !session.selection.contains(&pick) {
            session.selection.push(pick);
        }
    }
    session.refresh_display();
    let selected = state
        .bodies
        .iter()
        .filter(|b| session.selection.iter().any(|s| matches!(s, crate::Sel::Body { name } if name == &b.name)))
        .map(|b| b.name.clone())
        .collect();
    reply(PartsSelectReply { selected, missing })
}

/// Parts drawn see-through, so the parts inside an assembly show and can be picked.
pub(super) fn ghost(session: &mut Session, params: &Value) -> Result<Value> {
    let input: PartsGhost = parse(params)?;
    let state = session.world_state()?;
    let named_any = input.selection.is_some() || !input.bodies.is_empty();
    let bodies = if named_any { named(&state, &session.doc.multiphysics.parts, input.bodies, input.selection.as_deref())?.0 } else { Vec::new() };
    if !named_any && !input.only {
        return Err(bad("bodies", "name bodies or a part selection, or give only: true to draw every part solid"));
    }
    if input.only {
        session.ghosted_bodies.clear();
    }
    for name in bodies {
        if input.ghost {
            session.ghosted_bodies.insert(name);
        } else {
            session.ghosted_bodies.remove(&name);
        }
    }
    // Names the design no longer has are dropped, so the set stays as large as the design.
    let names: std::collections::HashSet<&str> = state.bodies.iter().map(|b| b.name.as_str()).collect();
    session.ghosted_bodies.retain(|g| names.contains(g.as_str()));
    session.refresh_display();
    let ghosted = state.bodies.iter().filter(|b| session.ghosted_bodies.contains(&b.name)).map(|b| b.name.clone()).collect();
    reply(PartsGhostReply { ghosted })
}

/// The study's bodies (in the design's order) and the design's.
fn study_bodies(session: &Session, study: &solidcraft_multiphysics::Study) -> Result<(Vec<String>, Vec<String>)> {
    let state = session.world_state()?;
    let parts = &session.doc.multiphysics.parts;
    let in_study: Vec<String> = match &study.model {
        Some(link) => super::model::bodies(&state, parts, link)?.iter().map(|b| b.name.clone()).collect(),
        None => vec![],
    };
    Ok((state.bodies.iter().map(|b| b.name.clone()).collect(), in_study))
}

/// The checklist of a study of the model: is it built from the design as it is now, where the
/// fluid is, openings, materials, limits, power and goals.
pub(super) fn check(session: &mut Session, params: &Value) -> Result<Value> {
    let input: PartsCheck = parse(params)?;
    let study = session.doc.multiphysics.study(&input.id).cloned().ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    if study.model.is_none() {
        return Err(bad("id", "the study is not a study of the model"));
    }
    let (_, in_study) = study_bodies(session, &study)?;
    let current = super::model::current(session, &study.id)?.1.unwrap_or(false);
    // The study's parts' exact volumes, for how much of each the grid holds. A mesh body has
    // no exact volume; a part that fails to measure is left without one.
    let state = session.world_state()?;
    let exact: std::collections::BTreeMap<String, f64> = state
        .bodies
        .iter()
        .filter(|b| in_study.contains(&b.name) && !b.body.is_mesh())
        .filter_map(|b| b.measure().ok().map(|m| (b.name.clone(), m.volume)))
        .collect();
    reply(solidcraft_multiphysics::setup_check::checklist_measured(&study, &in_study, current, &exact))
}

/// Colour the design's parts by their setup in a study, or (`by: none`) by their own looks.
pub(super) fn colour(session: &mut Session, params: &Value) -> Result<Value> {
    use solidcraft_multiphysics::setup_check::{ColourBy, Colouring};
    let input: PartsColour = parse(params)?;
    if input.by == ColourBy::None {
        session.part_swatches.clear();
        session.refresh_display();
        return reply(Colouring { by: ColourBy::None, parts: vec![], legend: vec![], max: None });
    }
    let id = input.id.ok_or_else(|| bad("id", "give the study whose setup colours the parts"))?;
    reply(colour_parts(session, &id, input.by, input.run_id.as_deref())?)
}

/// Colours the parts by `by` from study `id` (and `run_id`, default its latest thermal run).
pub(super) fn colour_parts(
    session: &mut Session,
    id: &str,
    by: solidcraft_multiphysics::setup_check::ColourBy,
    run_id: Option<&str>,
) -> Result<solidcraft_multiphysics::setup_check::Colouring> {
    let study = session.doc.multiphysics.study(id).cloned().ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    let (bodies, in_study) = study_bodies(session, &study)?;
    let parts = &session.doc.multiphysics.parts;
    let table =
        solidcraft_multiphysics::component_table::build(&study, parts, &bodies, &in_study, run_id).map_err(|e| bad("run_id", e.to_string()))?;
    let colouring = solidcraft_multiphysics::setup_check::colour(by, &study, parts, &bodies, &in_study, &table);
    session.part_swatches = colouring.parts.iter().map(|p| (p.body.clone(), p.swatch)).collect();
    session.refresh_display();
    Ok(colouring)
}

/// The study's parts that meet every rule of a find, on their setup and on a run's results.
pub(super) fn find(session: &mut Session, params: &Value) -> Result<Value> {
    let input: PartsFind = parse(params)?;
    let study = session.doc.multiphysics.study(&input.id).cloned().ok_or_else(|| bad("id", "Multiphysics study not found in this document"))?;
    let (bodies, in_study) = study_bodies(session, &study)?;
    let shown: Vec<String> = bodies.iter().filter(|b| input.all_bodies || in_study.contains(b)).cloned().collect();
    let parts = &session.doc.multiphysics.parts;
    let table = solidcraft_multiphysics::component_table::build(&study, parts, &shown, &in_study, input.run_id.as_deref())
        .map_err(|e| bad("run_id", e.to_string()))?;
    if table.run.is_none() && solidcraft_multiphysics::find::needs_run(&input.rules) {
        return Err(bad("rules", "a rule reads results, but the study has no run with a thermal report"));
    }
    let state = session.world_state()?;
    let measured = facts(&state, solidcraft_multiphysics::find::needs_measure(&input.rules, parts))?;
    let found = solidcraft_multiphysics::find::find(&input.rules, parts, &table, &measured).map_err(|e| bad("rules", e.to_string()))?;
    if input.select {
        session.selection = found.iter().map(|name| crate::Sel::Body { name: name.clone() }).collect();
        session.refresh_display();
    }
    reply(PartsFindReply { found, run: table.run })
}
