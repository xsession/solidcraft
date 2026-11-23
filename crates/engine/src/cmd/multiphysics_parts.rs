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
        let there = session.ghosted_bodies.iter().position(|g| g == &name);
        match (input.ghost, there) {
            (true, None) => session.ghosted_bodies.push(name),
            (false, Some(at)) => {
                session.ghosted_bodies.remove(at);
            }
            _ => {}
        }
    }
    // Names the design no longer has are dropped, so the list stays as long as the design.
    session.ghosted_bodies.retain(|g| state.bodies.iter().any(|b| &b.name == g));
    session.refresh_display();
    let ghosted = state.bodies.iter().filter(|b| session.ghosted_bodies.contains(&b.name)).map(|b| b.name.clone()).collect();
    reply(PartsGhostReply { ghosted })
}
