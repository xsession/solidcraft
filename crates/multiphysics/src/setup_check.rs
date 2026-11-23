//! A study of the model checked before it runs, and its parts coloured by their setup.
//!
//! The checklist reads the study's built case and its link to the model: whether the case is
//! the design as it is now, where the fluid is (the fluid volume check: no fluid left, sealed
//! pockets, parts too thin for the grid, parts sharing cells), the openings an internal flow
//! needs, and with heat on the parts without a material, power without a limit and goals.
//!
//! Colour-by sorts the parts into swatches (assigned or not, power, tag, the run's verdict);
//! the host draws each swatch in its own theme's colours, numeric ones from a colour map.
//!
//! Portable: the engine says which bodies the study holds and whether its case is current.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::Study;
use crate::component_table::ComponentTable;
use crate::model::{Analysis, OpeningRole};
use crate::parts::Registry;
use crate::results::ComponentStatus;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum CheckStatus {
    Ok,
    Warn,
    Error,
}

/// One line of the checklist.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Check {
    /// What is checked: built, fluid, openings, thin_parts, overlap, materials, limits, power, goals.
    pub check: String,
    pub status: CheckStatus,
    /// How many things the line is about (pockets, parts, watts rounded), for the host's wording.
    pub count: usize,
    /// Plain English for scripts and agents; hosts word it from `check`, `status` and `count`.
    pub message: String,
    /// The parts it is about, to select or colour.
    pub bodies: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Checklist {
    pub study: String,
    pub checks: Vec<Check>,
    pub errors: usize,
    pub warnings: usize,
}

fn line(check: &str, status: CheckStatus, count: usize, message: impl Into<String>, bodies: Vec<String>) -> Check {
    Check { check: check.into(), status, count, message: message.into(), bodies }
}

/// The bodies among `in_study` whose case part is one of `parts` (case names).
fn bodies_of(in_study: &[String], parts: &[String]) -> Vec<String> {
    in_study.iter().filter(|b| parts.contains(&crate::case_scene::part_name(b))).cloned().collect()
}

/// The checklist of `study`, a study of the model holding `in_study`; `current`: its case is
/// the design as it is now.
pub fn checklist(study: &Study, in_study: &[String], current: bool) -> Checklist {
    use CheckStatus::*;
    let mut checks = Vec::new();
    checks.push(if current {
        line("built", Ok, 0, "the case is built from the design as it is now", vec![])
    } else {
        line("built", Warn, 0, "the design changed since the case was built; build it again", vec![])
    });
    // Where the fluid is, as the solver's mesh will see it.
    match crate::fluid_volume::fluid_volume(&study.case) {
        Err(e) => checks.push(line("fluid", Warn, 0, format!("the fluid volume could not be checked: {e}"), vec![])),
        Result::Ok(volume) => {
            let sealed = volume.sealed().count();
            checks.push(if volume.fluid_cells == 0 {
                line("fluid", Error, 0, "no fluid cell is left: the parts fill the domain", vec![])
            } else if sealed > 0 {
                line("fluid", Warn, sealed, format!("{sealed} sealed pocket(s) of fluid that no opening reaches"), vec![])
            } else {
                let regions = volume.regions.len() + volume.unlisted_regions;
                line("fluid", Ok, regions, format!("{regions} fluid region(s), all reached by an opening"), vec![])
            });
            let mut thin: Vec<String> = volume.empty_bodies.clone();
            thin.extend(volume.solids.iter().filter(|s| s.cells == 0).map(|s| s.name.clone()));
            let thin = bodies_of(in_study, &thin);
            checks.push(if thin.is_empty() {
                line("thin_parts", Ok, 0, "every part takes grid cells", vec![])
            } else {
                line("thin_parts", Warn, thin.len(), "parts thinner than a cell are lost from the grid; refine the mesh around them", thin)
            });
            let shared: Vec<String> = volume.solids.iter().filter(|s| s.overlapping_cells > 0).map(|s| s.name.clone()).collect();
            let shared = bodies_of(in_study, &shared);
            checks.push(if shared.is_empty() {
                line("overlap", Ok, 0, "no two parts share grid cells", vec![])
            } else {
                line("overlap", Warn, shared.len(), "parts share grid cells with an earlier part; their heat goes to the earlier one", shared)
            });
        }
    }
    let Some(link) = &study.model else {
        return finish(study, checks);
    };
    if link.analysis == Analysis::Internal {
        let inlets = link.openings.iter().filter(|o| matches!(o.role, OpeningRole::Inlet { .. })).count();
        let outlets = link.openings.iter().filter(|o| matches!(o.role, OpeningRole::Outlet { .. })).count();
        let fan = usize::from(link.fan.is_some());
        checks.push(if inlets + fan == 0 || outlets == 0 {
            line("openings", Error, inlets + outlets, "an internal flow needs an inlet (or a fan) and an outlet", vec![])
        } else {
            line("openings", Ok, inlets + outlets, "the flow has an inlet and an outlet", vec![])
        });
    }
    if let Some(heat) = &link.heat {
        let unassigned: Vec<String> = in_study.iter().filter(|b| !heat.materials.iter().any(|m| &m.body == *b)).cloned().collect();
        checks.push(if unassigned.is_empty() {
            line("materials", Ok, 0, "every part has a material", vec![])
        } else {
            line("materials", Warn, unassigned.len(), "parts without a material take the default solid", unassigned)
        });
        let powered = |body: &str| {
            heat.materials.iter().any(|m| m.body == body && m.power_w > 0.0)
                || heat.face_sources.iter().any(|s| s.face.body == body && s.power_w > 0.0)
        };
        let unlimited: Vec<String> = in_study
            .iter()
            .filter(|b| powered(b) && heat.materials.iter().find(|m| &m.body == *b).is_none_or(|m| m.max_temperature_k.is_none()))
            .cloned()
            .collect();
        checks.push(if unlimited.is_empty() {
            line("limits", Ok, 0, "every part with power has a temperature limit", vec![])
        } else {
            line("limits", Warn, unlimited.len(), "parts with power but no temperature limit get no verdict", unlimited)
        });
        let total: f64 = heat.materials.iter().map(|m| m.power_w).sum::<f64>() + heat.face_sources.iter().map(|s| s.power_w).sum::<f64>();
        checks.push(if total > 0.0 {
            line("power", Ok, total.round() as usize, format!("{total} W in the study"), vec![])
        } else {
            line("power", Warn, 0, "heat transfer is on but no part gives off heat", vec![])
        });
    }
    checks.push(if link.goals.is_empty() {
        line("goals", Warn, 0, "no goals: the run reports no values to watch or compare", vec![])
    } else {
        line("goals", Ok, link.goals.len(), format!("{} goal(s)", link.goals.len()), vec![])
    });
    finish(study, checks)
}

fn finish(study: &Study, checks: Vec<Check>) -> Checklist {
    let errors = checks.iter().filter(|c| c.status == CheckStatus::Error).count();
    let warnings = checks.iter().filter(|c| c.status == CheckStatus::Warn).count();
    Checklist { study: study.id.clone(), checks, errors, warnings }
}

/// What the parts are coloured by.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ColourBy {
    /// Their own appearance.
    #[default]
    None,
    /// With a material or without (unassigned parts stand out).
    Material,
    /// Their power, on a colour map.
    Power,
    /// Their first tag.
    Tag,
    /// The run's verdict against their limit.
    Verdict,
}

/// How a host draws a part: its theme's status colours, a muted body colour, a colour map
/// position, or one of its categorical colours.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Swatch {
    Good,
    Warn,
    Bad,
    Muted,
    /// Position on the colour map, 0 to 1.
    Scale {
        t: f64,
    },
    Category {
        index: u8,
    },
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct PartSwatch {
    pub body: String,
    pub swatch: Swatch,
    /// The legend entry it belongs to.
    pub label: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct LegendEntry {
    pub swatch: Swatch,
    /// An id the host words (assigned, unassigned, outside, no_heat, pass, warn, fail,
    /// no_limit, no_result, untagged), a tag, or for power the value in W.
    pub label: String,
    pub count: usize,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Colouring {
    pub by: ColourBy,
    pub parts: Vec<PartSwatch>,
    pub legend: Vec<LegendEntry>,
    /// For power: the highest value, W, at the top of the map.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max: Option<f64>,
}

/// Most tags given a colour of their own; later ones share the last.
pub const MAX_TAG_COLOURS: usize = 6;

/// The design's `bodies` coloured `by` for `study`, which holds `in_study`; `table` is its
/// component table (for power and verdicts).
pub fn colour(by: ColourBy, study: &Study, registry: &Registry, bodies: &[String], in_study: &[String], table: &ComponentTable) -> Colouring {
    let heat = study.model.as_ref().and_then(|l| l.heat.as_ref());
    let mut max = None;
    let mut tags: Vec<String> = Vec::new();
    let parts: Vec<PartSwatch> = bodies
        .iter()
        .filter(|_| by != ColourBy::None)
        .map(|body| {
            let (swatch, label) = if !in_study.contains(body) {
                (Swatch::Muted, "outside".to_string())
            } else {
                match by {
                    ColourBy::None | ColourBy::Material => match heat {
                        None => (Swatch::Muted, "no_heat".into()),
                        Some(h) if h.materials.iter().any(|m| &m.body == body) => (Swatch::Good, "assigned".into()),
                        Some(_) => (Swatch::Bad, "unassigned".into()),
                    },
                    ColourBy::Power => {
                        let top = table.rows.iter().map(|r| r.power_w).fold(0.0_f64, f64::max);
                        max = Some(top);
                        let power = table.rows.iter().find(|r| &r.body == body).map_or(0.0, |r| r.power_w);
                        let t = if top > 0.0 { power / top } else { 0.0 };
                        (Swatch::Scale { t }, format!("{power}"))
                    }
                    ColourBy::Tag => match registry.parts.get(body).and_then(|p| p.tags.iter().next()) {
                        None => (Swatch::Muted, "untagged".into()),
                        Some(tag) => {
                            let index = tags.iter().position(|t| t == tag).unwrap_or_else(|| {
                                tags.push(tag.clone());
                                tags.len() - 1
                            });
                            (Swatch::Category { index: index.min(MAX_TAG_COLOURS - 1) as u8 }, tag.clone())
                        }
                    },
                    ColourBy::Verdict => match table.rows.iter().find(|r| &r.body == body) {
                        Some(row) if row.result.is_some() => match row.status {
                            ComponentStatus::Pass => (Swatch::Good, "pass".into()),
                            ComponentStatus::Warn => (Swatch::Warn, "warn".into()),
                            ComponentStatus::Fail => (Swatch::Bad, "fail".into()),
                            ComponentStatus::NoLimit => (Swatch::Muted, "no_limit".into()),
                        },
                        _ => (Swatch::Muted, "no_result".into()),
                    },
                }
            };
            PartSwatch { body: body.clone(), swatch, label }
        })
        .collect();
    // One legend entry per label, in the order parts first show it; power is a range instead.
    let mut legend: Vec<LegendEntry> = Vec::new();
    for part in &parts {
        if matches!(part.swatch, Swatch::Scale { .. }) {
            continue;
        }
        match legend.iter_mut().find(|e| e.label == part.label) {
            Some(entry) => entry.count += 1,
            None => legend.push(LegendEntry { swatch: part.swatch, label: part.label.clone(), count: 1 }),
        }
    }
    Colouring { by, parts, legend, max }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{BodyHeat, ModelHeat, ModelLink};

    fn study(heat: Option<ModelHeat>) -> Study {
        let mut project = crate::Project::default();
        let mut study = project.create("flow-channel", "S").expect("study");
        study.model = Some(ModelLink { analysis: Analysis::External, heat, ..ModelLink::default() });
        study
    }

    fn heat(materials: &[(&str, f64, Option<f64>)]) -> ModelHeat {
        let materials = materials
            .iter()
            .map(|(body, power, limit)| BodyHeat {
                body: (*body).into(),
                conductivity_w_per_m_k: 1.0,
                power_w: *power,
                density_kg_per_m3: None,
                specific_heat_j_per_kg_k: None,
                max_temperature_k: *limit,
            })
            .collect();
        ModelHeat { materials, ..ModelHeat::default() }
    }

    fn table(study: &Study, registry: &Registry, bodies: &[String]) -> ComponentTable {
        crate::component_table::build(study, registry, bodies, bodies, None).expect("table")
    }

    #[test]
    fn the_checklist_names_unassigned_parts_and_power_without_a_limit() {
        let bodies: Vec<String> = ["U1", "U2", "Board"].map(String::from).to_vec();
        let s = study(Some(heat(&[("U1", 2.0, Some(358.0)), ("U2", 1.0, None)])));
        let list = checklist(&s, &bodies, false);
        let get = |id: &str| list.checks.iter().find(|c| c.check == id).cloned().expect(id);
        assert_eq!(get("built").status, CheckStatus::Warn);
        assert_eq!(get("materials").bodies, ["Board"]);
        assert_eq!(get("limits").bodies, ["U2"]);
        assert_eq!((get("power").status, get("power").count), (CheckStatus::Ok, 3));
        assert_eq!(get("goals").status, CheckStatus::Warn);
        // The template's channel has fluid that its openings reach.
        assert_eq!(get("fluid").status, CheckStatus::Ok, "{:?}", get("fluid"));
        assert_eq!(list.errors, 0);
        assert!(list.warnings >= 4, "{list:?}");
    }

    #[test]
    fn colour_by_material_makes_unassigned_parts_stand_out() {
        let bodies: Vec<String> = ["U1", "U2", "Board", "Lid"].map(String::from).to_vec();
        let in_study = &bodies[..3];
        let s = study(Some(heat(&[("U1", 2.0, None), ("U2", 1.0, None)])));
        let registry = Registry::default();
        let t = table(&s, &registry, &bodies);
        let c = colour(ColourBy::Material, &s, &registry, &bodies, in_study, &t);
        let swatches: Vec<Swatch> = c.parts.iter().map(|p| p.swatch).collect();
        assert_eq!(swatches, [Swatch::Good, Swatch::Good, Swatch::Bad, Swatch::Muted]);
        let legend: Vec<(&str, usize)> = c.legend.iter().map(|e| (e.label.as_str(), e.count)).collect();
        assert_eq!(legend, [("assigned", 2), ("unassigned", 1), ("outside", 1)]);
        // Power on the map: the hottest part at the top, an unpowered one at the bottom.
        let p = colour(ColourBy::Power, &s, &registry, &bodies, in_study, &t);
        assert_eq!(p.max, Some(2.0));
        assert_eq!(p.parts[1].swatch, Swatch::Scale { t: 0.5 });
        assert_eq!(p.parts[2].swatch, Swatch::Scale { t: 0.0 });
        // Nothing before a run: every verdict waits for one.
        let v = colour(ColourBy::Verdict, &s, &registry, &bodies, in_study, &t);
        assert!(v.parts[..3].iter().all(|p| p.label == "no_result"));
        assert!(colour(ColourBy::None, &s, &registry, &bodies, in_study, &t).parts.is_empty());
    }

    #[test]
    fn colour_by_tag_gives_each_tag_its_colour() {
        let bodies: Vec<String> = ["U1", "U2", "C1"].map(String::from).to_vec();
        let s = study(None);
        let mut registry = Registry::default();
        let retag = |tag: &str| crate::parts::Retag { add_tags: vec![tag.into()], remove_tags: vec![], attributes: Default::default() };
        registry.retag(&bodies[..2], &retag("ic")).expect("tag");
        let t = table(&s, &registry, &bodies);
        let c = colour(ColourBy::Tag, &s, &registry, &bodies, &bodies, &t);
        let swatches: Vec<Swatch> = c.parts.iter().map(|p| p.swatch).collect();
        assert_eq!(swatches, [Swatch::Category { index: 0 }, Swatch::Category { index: 0 }, Swatch::Muted]);
        assert_eq!(c.legend.iter().map(|e| e.label.as_str()).collect::<Vec<_>>(), ["ic", "untagged"]);
    }
}
