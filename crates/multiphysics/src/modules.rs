//! Modules: the Multiphysics setup of one component, kept to be given again to every copy of it
//! (a board placed four times in a rack) and, through a file, to other designs.
//!
//! A module names its bodies as the component knows them (its own body names, and for a body of
//! a component inside it, that component's name first: `Sub/U1`), never as one placement names
//! them, so it fits every placement. It carries, per body, the material, power and temperature
//! limit; the tags and attributes; and whether a goal watches its highest temperature; and heat
//! sources on the component's faces, by face name.
//!
//! Applying a module to a study writes those into the study and the part registry for each
//! placement: a body's material and power, a heat source of the same name, a goal on the same
//! body are replaced, so applying twice changes nothing more. Bodies a placement does not have
//! are skipped and listed.
//!
//! Portable: the engine says which placed body is which body of the component.
use std::collections::{BTreeMap, BTreeSet};

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::model::{BodyHeat, FaceHeat, FaceRef, ModelGoal, ModelGoalKind, ModelHeat, ModelLink};
use crate::parts::{AttrValue, PartInfo, Registry};
use crate::{Error, Result};

/// What a module file says it is.
pub const FILE_FORMAT: &str = "solidcraft-multiphysics-module/1";
/// Most modules a design keeps.
pub const MAX_MODULES: usize = 64;
/// Most bodies, and most face heat sources, one module carries.
pub const MAX_BODIES: usize = 2048;
/// Largest module file read.
pub const MAX_FILE_BYTES: usize = 4 * 1024 * 1024;

/// A body's solid material, as a study gives it (see [`BodyHeat`]).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Material {
    pub conductivity_w_per_m_k: f64,
    #[serde(default)]
    pub power_w: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub density_kg_per_m3: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specific_heat_j_per_kg_k: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_temperature_k: Option<f64>,
}

/// One body of the component and its setup.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModuleBody {
    /// The body as the component names it (`Sub/U1` for a body of a component inside it).
    pub body: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub material: Option<Material>,
    #[serde(default, skip_serializing_if = "BTreeSet::is_empty")]
    pub tags: BTreeSet<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub attributes: BTreeMap<String, AttrValue>,
    /// A goal watches the body's highest temperature.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub max_temperature_goal: bool,
}

/// Heat released on a face of one of the component's bodies.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModuleFaceHeat {
    pub name: String,
    pub body: String,
    pub face: String,
    pub power_w: f64,
}

/// The setup of a component, by its definition's name.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Module {
    pub name: String,
    /// The component definition it fits, by name.
    pub component: String,
    #[serde(default)]
    pub bodies: Vec<ModuleBody>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub face_sources: Vec<ModuleFaceHeat>,
}

/// One placement of the component: a label for it (its occurrences' names) and its bodies,
/// each as the design names it and as the component does.
#[derive(Clone, Debug, PartialEq)]
pub struct Placement {
    pub label: String,
    /// (placed body name, the body as the component names it).
    pub bodies: Vec<(String, String)>,
}

fn invalid(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

fn name(value: &str, what: &str) -> Result<()> {
    if value.trim().is_empty() || value.chars().count() > 200 || value.chars().any(char::is_control) {
        return Err(invalid(format!("{what} is 1 to 200 characters without control characters")));
    }
    Ok(())
}

/// Refuses a module that could not be applied as written.
pub fn check(module: &Module) -> Result<()> {
    name(&module.name, "a module's name")?;
    name(&module.component, "a module's component")?;
    if module.bodies.len() > MAX_BODIES || module.face_sources.len() > MAX_BODIES {
        return Err(Error::Limit(format!("a module carries at most {MAX_BODIES} bodies and {MAX_BODIES} heat sources")));
    }
    let mut seen = BTreeSet::new();
    let mut registry = Registry::default();
    for body in &module.bodies {
        name(&body.body, "a module body")?;
        if !seen.insert(body.body.as_str()) {
            return Err(invalid(format!("{} is in the module twice", body.body)));
        }
        if let Some(m) = &body.material {
            let positive = |v: f64| v.is_finite() && v > 0.0;
            if !positive(m.conductivity_w_per_m_k) || !m.power_w.is_finite() || m.power_w < 0.0 {
                return Err(invalid(format!("{}: a positive conductivity and a power of 0 W or more", body.body)));
            }
            if [m.density_kg_per_m3, m.specific_heat_j_per_kg_k, m.max_temperature_k].into_iter().flatten().any(|v| !positive(v)) {
                return Err(invalid(format!("{}: density, specific heat and limit are positive", body.body)));
            }
        }
        // Tags and attributes follow the registry's own rules.
        let info = PartInfo { tags: body.tags.clone(), attributes: body.attributes.clone() };
        if !info.is_empty() {
            registry.parts.insert(body.body.clone(), info);
        }
    }
    registry.validate()?;
    let mut sources = BTreeSet::new();
    for source in &module.face_sources {
        name(&source.name, "a heat source's name")?;
        name(&source.body, "a heat source's body")?;
        name(&source.face, "a heat source's face")?;
        if !(source.power_w.is_finite() && source.power_w >= 0.0) {
            return Err(invalid(format!("{}: a power of 0 W or more", source.name)));
        }
        if !sources.insert(source.name.as_str()) {
            return Err(invalid(format!("two heat sources are named {}", source.name)));
        }
    }
    Ok(())
}

/// The module of `component` as `placement` of it is set up in a study's `link` and the part
/// `registry`: every body with something set, and the heat sources on its faces.
pub fn capture(module_name: &str, component: &str, placement: &Placement, link: &ModelLink, registry: &Registry) -> Result<Module> {
    let relative: BTreeMap<&str, &str> = placement.bodies.iter().map(|(world, own)| (world.as_str(), own.as_str())).collect();
    let mut bodies = Vec::new();
    for (world, own) in &placement.bodies {
        let material = link.heat.as_ref().and_then(|h| h.materials.iter().find(|m| &m.body == world)).map(|m| Material {
            conductivity_w_per_m_k: m.conductivity_w_per_m_k,
            power_w: m.power_w,
            density_kg_per_m3: m.density_kg_per_m3,
            specific_heat_j_per_kg_k: m.specific_heat_j_per_kg_k,
            max_temperature_k: m.max_temperature_k,
        });
        let info = registry.parts.get(world).cloned().unwrap_or_default();
        let goal = link.goals.iter().any(|g| g.kind == ModelGoalKind::BodyMaxTemperature && g.body.as_deref() == Some(world.as_str()));
        if material.is_some() || !info.is_empty() || goal {
            bodies.push(ModuleBody { body: own.clone(), material, tags: info.tags, attributes: info.attributes, max_temperature_goal: goal });
        }
    }
    let face_sources = link
        .heat
        .iter()
        .flat_map(|h| &h.face_sources)
        .filter_map(|s| {
            let own = relative.get(s.face.body.as_str())?;
            Some(ModuleFaceHeat { name: s.name.clone(), body: (*own).to_string(), face: s.face.face.clone(), power_w: s.power_w })
        })
        .collect();
    let module = Module { name: module_name.trim().to_string(), component: component.to_string(), bodies, face_sources };
    check(&module)?;
    Ok(module)
}

/// What applying a module gave.
#[derive(Clone, Debug, PartialEq)]
pub struct Applied {
    pub link: ModelLink,
    pub registry: Registry,
    /// Placements it was applied to, and bodies set in them.
    pub placements: usize,
    pub bodies: usize,
    /// `placement/body` of the module's bodies a placement does not have.
    pub missing: Vec<String>,
}

/// A goal name (letters, digits, underscores, a letter first) from free text.
fn goal_name(parts: &[&str]) -> String {
    let joined: String = parts.join("_").chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect();
    let squeezed = joined.split('_').filter(|s| !s.is_empty()).collect::<Vec<_>>().join("_");
    let start = if squeezed.starts_with(|c: char| c.is_ascii_alphabetic()) { "" } else { "m_" };
    format!("{start}{squeezed}").chars().take(64).collect()
}

/// `module` written into `link` and `registry` for each of `placements`.
pub fn apply(module: &Module, placements: &[Placement], link: &ModelLink, registry: &Registry) -> Result<Applied> {
    check(module)?;
    if placements.is_empty() {
        return Err(invalid(format!("the design has no placement of {}", module.component)));
    }
    let mut link = link.clone();
    let mut registry = registry.clone();
    let mut missing = Vec::new();
    let mut set = 0usize;
    let needs_heat = module.bodies.iter().any(|b| b.material.is_some() || b.max_temperature_goal) || !module.face_sources.is_empty();
    if needs_heat && link.heat.is_none() {
        link.heat = Some(ModelHeat::default());
    }
    for placement in placements {
        let world: BTreeMap<&str, &str> = placement.bodies.iter().map(|(world, own)| (own.as_str(), world.as_str())).collect();
        for body in &module.bodies {
            let Some(target) = world.get(body.body.as_str()).map(|w| (*w).to_string()) else {
                missing.push(format!("{}/{}", placement.label, body.body));
                continue;
            };
            set += 1;
            if let (Some(m), Some(heat)) = (&body.material, link.heat.as_mut()) {
                let row = BodyHeat {
                    body: target.clone(),
                    conductivity_w_per_m_k: m.conductivity_w_per_m_k,
                    power_w: m.power_w,
                    density_kg_per_m3: m.density_kg_per_m3,
                    specific_heat_j_per_kg_k: m.specific_heat_j_per_kg_k,
                    max_temperature_k: m.max_temperature_k,
                };
                match heat.materials.iter_mut().find(|r| r.body == target) {
                    Some(slot) => *slot = row,
                    None => heat.materials.push(row),
                }
            }
            if !body.tags.is_empty() || !body.attributes.is_empty() {
                let info = registry.parts.entry(target.clone()).or_default();
                info.tags.extend(body.tags.iter().cloned());
                info.attributes.extend(body.attributes.iter().map(|(k, v)| (k.clone(), v.clone())));
            }
            if body.max_temperature_goal {
                let goal = ModelGoal {
                    name: goal_name(&[&module.name, &placement.label, &body.body, "max"]),
                    kind: ModelGoalKind::BodyMaxTemperature,
                    face: None,
                    point_mm: None,
                    expression: None,
                    body: Some(target.clone()),
                };
                match link.goals.iter_mut().find(|g| g.name == goal.name) {
                    Some(slot) => *slot = goal,
                    None => link.goals.push(goal),
                }
            }
        }
        for source in &module.face_sources {
            let Some(target) = world.get(source.body.as_str()) else {
                missing.push(format!("{}/{}", placement.label, source.body));
                continue;
            };
            let row = FaceHeat {
                name: format!("{} · {}", source.name, placement.label),
                face: FaceRef { body: (*target).to_string(), face: source.face.clone() },
                power_w: source.power_w,
            };
            if let Some(heat) = link.heat.as_mut() {
                match heat.face_sources.iter_mut().find(|s| s.name == row.name) {
                    Some(slot) => *slot = row,
                    None => heat.face_sources.push(row),
                }
            }
        }
    }
    registry.validate()?;
    link.validate()?;
    missing.sort();
    missing.dedup();
    Ok(Applied { link, registry, placements: placements.len(), bodies: set, missing })
}

/// A module as a file of its own: JSON that says what it is.
pub fn to_file(module: &Module) -> Result<String> {
    check(module)?;
    serde_json::to_string_pretty(&serde_json::json!({"format": FILE_FORMAT, "module": module})).map_err(|e| Error::Json(e.to_string()))
}

/// A module read from a file [`to_file`] wrote; refused if it is not one or does not check.
pub fn from_file(text: &str) -> Result<Module> {
    if text.len() > MAX_FILE_BYTES {
        return Err(Error::Limit(format!("a module file is at most {MAX_FILE_BYTES} bytes")));
    }
    let value: serde_json::Value = serde_json::from_str(text).map_err(|e| Error::Json(e.to_string()))?;
    if value.get("format").and_then(serde_json::Value::as_str) != Some(FILE_FORMAT) {
        return Err(invalid(format!("not a module file (format {FILE_FORMAT})")));
    }
    let module: Module = serde_json::from_value(value.get("module").cloned().unwrap_or_default()).map_err(|e| Error::Json(e.to_string()))?;
    check(&module)?;
    Ok(module)
}

/// Keeps `module` in `modules`, in place of one of the same name; refused past [`MAX_MODULES`].
pub fn save(modules: &mut Vec<Module>, module: Module) -> Result<bool> {
    check(&module)?;
    if let Some(slot) = modules.iter_mut().find(|m| m.name == module.name) {
        *slot = module;
        return Ok(true);
    }
    if modules.len() >= MAX_MODULES {
        return Err(Error::Limit(format!("a design keeps at most {MAX_MODULES} modules")));
    }
    modules.push(module);
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn placement(label: &str) -> Placement {
        Placement { label: label.into(), bodies: ["PCB", "U1", "Sub/R1"].iter().map(|own| (format!("{label}/{own}"), (*own).to_string())).collect() }
    }

    fn study_of_first() -> (ModelLink, Registry) {
        let link = ModelLink {
            heat: Some(ModelHeat {
                materials: vec![BodyHeat {
                    body: "Board 1/U1".into(),
                    conductivity_w_per_m_k: 150.0,
                    power_w: 0.8,
                    density_kg_per_m3: None,
                    specific_heat_j_per_kg_k: None,
                    max_temperature_k: Some(358.15),
                }],
                face_sources: vec![FaceHeat { name: "pad".into(), face: FaceRef { body: "Board 1/PCB".into(), face: "F3:+z".into() }, power_w: 0.2 }],
                ..Default::default()
            }),
            goals: vec![ModelGoal {
                name: "u1_hot".into(),
                kind: ModelGoalKind::BodyMaxTemperature,
                face: None,
                point_mm: None,
                expression: None,
                body: Some("Board 1/U1".into()),
            }],
            ..Default::default()
        };
        let mut registry = Registry::default();
        registry.parts.insert("Board 1/Sub/R1".into(), PartInfo { tags: ["resistor".to_string()].into(), attributes: BTreeMap::new() });
        (link, registry)
    }

    /// Set up on one board, captured as the board's module, given to the other three: each
    /// gets the material, the tag, the goal and the face heat, on its own bodies.
    #[test]
    fn a_module_captured_from_one_placement_sets_up_every_other() {
        let (link, registry) = study_of_first();
        let module = capture("Board setup", "Board", &placement("Board 1"), &link, &registry).unwrap();
        assert_eq!(module.bodies.len(), 2, "{module:?}");
        let u1 = module.bodies.iter().find(|b| b.body == "U1").unwrap();
        assert_eq!(u1.material.map(|m| m.power_w), Some(0.8));
        assert!(u1.max_temperature_goal);
        assert_eq!(module.bodies.iter().find(|b| b.body == "Sub/R1").unwrap().tags, ["resistor".to_string()].into());
        assert_eq!(module.face_sources, vec![ModuleFaceHeat { name: "pad".into(), body: "PCB".into(), face: "F3:+z".into(), power_w: 0.2 }]);
        // Through a file, as another design would read it.
        let module = from_file(&to_file(&module).unwrap()).unwrap();
        let placements: Vec<Placement> = (1..=4).map(|i| placement(&format!("Board {i}"))).collect();
        let applied = apply(&module, &placements, &link, &registry).unwrap();
        assert_eq!((applied.placements, applied.bodies, applied.missing.len()), (4, 8, 0));
        let heat = applied.link.heat.as_ref().unwrap();
        for i in 1..=4 {
            let u1 = heat.materials.iter().find(|m| m.body == format!("Board {i}/U1")).unwrap();
            assert_eq!((u1.power_w, u1.max_temperature_k), (0.8, Some(358.15)));
            assert!(applied.registry.parts[&format!("Board {i}/Sub/R1")].tags.contains("resistor"));
            assert!(heat.face_sources.iter().any(|s| s.face.body == format!("Board {i}/PCB") && s.power_w == 0.2));
        }
        assert_eq!(heat.materials.len(), 4, "the first board's own row is replaced, not doubled");
        assert_eq!(
            applied.link.goals.iter().filter(|g| g.kind == ModelGoalKind::BodyMaxTemperature).count(),
            5,
            "the study's own goal and one per board"
        );
        // Applying again changes nothing.
        let again = apply(&module, &placements, &applied.link, &applied.registry).unwrap();
        assert_eq!((again.link, again.registry), (applied.link, applied.registry));
    }

    #[test]
    fn missing_bodies_are_listed_and_bad_modules_refused() {
        let (link, registry) = study_of_first();
        let module = capture("Board setup", "Board", &placement("Board 1"), &link, &registry).unwrap();
        let short = Placement { label: "Board 9".into(), bodies: vec![("Board 9/PCB".into(), "PCB".into())] };
        let applied = apply(&module, &[short], &link, &registry).unwrap();
        assert_eq!(applied.missing, ["Board 9/Sub/R1", "Board 9/U1"]);
        assert!(apply(&module, &[], &link, &registry).is_err(), "no placement");
        let mut bad = module.clone();
        bad.bodies[0].material = Some(Material {
            conductivity_w_per_m_k: -1.0,
            power_w: 0.0,
            density_kg_per_m3: None,
            specific_heat_j_per_kg_k: None,
            max_temperature_k: None,
        });
        assert!(check(&bad).is_err());
        let mut twice = module.clone();
        twice.bodies.push(twice.bodies[0].clone());
        assert!(check(&twice).is_err());
        assert!(from_file("{\"format\":\"something else\",\"module\":{}}").is_err());
        assert!(from_file("not json").is_err());
        let mut modules = Vec::new();
        assert!(!save(&mut modules, module.clone()).unwrap());
        assert!(save(&mut modules, module).unwrap(), "same name replaces");
        assert_eq!(goal_name(&["Board setup", "Rack/Board 1", "U1", "max"]), "Board_setup_Rack_Board_1_U1_max");
        assert_eq!(goal_name(&["1st"]), "m_1st");
    }
}
