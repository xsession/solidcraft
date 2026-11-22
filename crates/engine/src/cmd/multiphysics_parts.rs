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
    let bodies = match (&input.selection, input.bodies.is_empty()) {
        (Some(_), false) => return Err(bad("selection", "name bodies or a part selection, not both")),
        (Some(selection), true) => resolve(&state, &parts, selection)?.0,
        (None, _) => {
            if let Some(unknown) = input.bodies.iter().find(|n| !state.bodies.iter().any(|b| &&b.name == n)) {
                return Err(bad("bodies", format!("the design has no body called {unknown}")));
            }
            input.bodies
        }
    };
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
