//! The component table: one row per part of a study, its setup (tags, attributes, material,
//! power, temperature limit) beside its result from a run (temperatures, margin, verdict), with
//! the study's total power, and as CSV.
//!
//! Portable: the engine says which bodies the study holds; results come from the run's
//! `thermal_report.json`, matched to bodies by their case part name (`case_scene::part_name`).
use std::collections::BTreeMap;

use schemars::JsonSchema;
use serde::Serialize;

use crate::model::ModelLink;
use crate::parts::{AttrValue, Registry};
use crate::results::{ComponentRow, ComponentStatus, components_from_report};
use crate::{Error, Result, RunEvidence, Study};

/// One part's row.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct TableRow {
    pub body: String,
    pub tags: Vec<String>,
    pub attributes: BTreeMap<String, AttrValue>,
    /// The body is one of the study's bodies (a row may show a design body the study leaves out).
    pub in_study: bool,
    /// W/(m·K), where the study gives the body a material.
    pub conductivity_w_per_m_k: Option<f64>,
    /// Heat released in the body and on its faces, W.
    pub power_w: f64,
    /// The temperature the body must stay below, K.
    pub limit_k: Option<f64>,
    /// The run's temperatures for the body, where the run reports it.
    pub result: Option<ComponentRow>,
    /// `limit_k` less the run's highest temperature, K.
    pub margin_k: Option<f64>,
    pub status: ComponentStatus,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct ComponentTable {
    pub study: String,
    /// The run the results come from; none before a run with a thermal report.
    pub run: Option<String>,
    /// The run was made from an earlier revision of the study, so its results may not match the setup.
    pub stale: bool,
    pub rows: Vec<TableRow>,
    /// Power of the rows in the study, W.
    pub total_power_w: f64,
    /// The row in the study with the least margin, if any has one.
    pub worst: Option<String>,
}

/// The run results come from: `run_id`, or the latest run that kept a thermal report.
fn pick_run<'a>(study: &'a Study, run_id: Option<&str>) -> Result<Option<&'a RunEvidence>> {
    match run_id {
        Some(id) => {
            study.runs.iter().find(|r| r.id == id).map(Some).ok_or_else(|| Error::NotFound(format!("run {id} is not retained in this study")))
        }
        None => Ok(study.runs.iter().rev().find(|r| r.artifacts.iter().any(|a| a.name == "thermal_report.json"))),
    }
}

/// The table for `study` over `bodies` (the design's bodies, in its order) of which the study
/// holds `in_study`.
pub fn build(study: &Study, registry: &Registry, bodies: &[String], in_study: &[String], run_id: Option<&str>) -> Result<ComponentTable> {
    let run = pick_run(study, run_id)?;
    let results: Vec<ComponentRow> = match run {
        Some(run) => match run.artifacts.iter().find(|a| a.name == "thermal_report.json") {
            Some(report) => components_from_report(&report.text).map_err(|e| Error::Invalid(e.to_string()))?,
            None => vec![],
        },
        None => vec![],
    };
    let link: Option<&ModelLink> = study.model.as_ref();
    let heat = link.and_then(|l| l.heat.as_ref());
    let mut rows = Vec::with_capacity(bodies.len());
    for body in bodies {
        let info = registry.parts.get(body);
        let material = heat.and_then(|h| h.materials.iter().find(|m| &m.body == body));
        let face_power: f64 = heat.map_or(0.0, |h| h.face_sources.iter().filter(|s| &s.face.body == body).map(|s| s.power_w).sum());
        let part = crate::case_scene::part_name(body);
        let result = results.iter().find(|r| r.name == part || &r.name == body).cloned();
        let limit_k = material.and_then(|m| m.max_temperature_k).or_else(|| result.as_ref().and_then(|r| r.limit_k));
        let margin_k = match (limit_k, &result) {
            (Some(limit), Some(r)) => Some(limit - r.max_temperature_k),
            _ => None,
        };
        let status = match (margin_k, &result) {
            (Some(m), _) if m < 0.0 => ComponentStatus::Fail,
            (Some(_), Some(r)) if r.status == ComponentStatus::Warn => ComponentStatus::Warn,
            (Some(_), _) => ComponentStatus::Pass,
            (None, _) => ComponentStatus::NoLimit,
        };
        rows.push(TableRow {
            body: body.clone(),
            tags: info.map(|i| i.tags.iter().cloned().collect()).unwrap_or_default(),
            attributes: info.map(|i| i.attributes.clone()).unwrap_or_default(),
            in_study: in_study.contains(body),
            conductivity_w_per_m_k: material.map(|m| m.conductivity_w_per_m_k),
            power_w: material.map_or(0.0, |m| m.power_w) + face_power,
            limit_k,
            result,
            margin_k,
            status,
        });
    }
    let total_power_w = rows.iter().filter(|r| r.in_study).map(|r| r.power_w).sum();
    let worst = rows
        .iter()
        .filter(|r| r.in_study)
        .filter_map(|r| r.margin_k.map(|m| (m, &r.body)))
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, b)| b.clone());
    let stale = run.is_some_and(|r| r.study_revision != study.revision);
    Ok(ComponentTable { study: study.name.clone(), run: run.map(|r| r.id.clone()), stale, rows, total_power_w, worst })
}

/// One CSV cell: quoted where it must be, and a text that a spreadsheet would read as a formula
/// starts with an apostrophe.
/// A CSV cell: quoted where it must be, and never read as a formula.
pub(crate) fn cell(text: &str) -> String {
    let mut t = text.to_string();
    if t.starts_with(['=', '+', '-', '@', '\t', '\r']) && t.parse::<f64>().is_err() {
        t.insert(0, '\'');
    }
    if t.contains([',', '"', '\n', '\r']) { format!("\"{}\"", t.replace('"', "\"\"")) } else { t }
}

pub(crate) fn number(v: Option<f64>) -> String {
    v.filter(|v| v.is_finite()).map(|v| format!("{v}")).unwrap_or_default()
}

/// A verdict as CSV and reports write it.
pub(crate) fn status_name(status: ComponentStatus) -> &'static str {
    match status {
        ComponentStatus::Pass => "pass",
        ComponentStatus::Warn => "warn",
        ComponentStatus::Fail => "fail",
        ComponentStatus::NoLimit => "no_limit",
    }
}

/// The table as CSV: a header, one line per row, and a total line; attributes as `key=value` pairs.
pub fn to_csv(table: &ComponentTable) -> String {
    let mut out =
        String::from("body,in_study,tags,attributes,conductivity_w_per_m_k,power_w,limit_k,max_temperature_k,mean_temperature_k,margin_k,status\n");
    for r in &table.rows {
        let attributes: Vec<String> = r
            .attributes
            .iter()
            .map(|(k, v)| match v {
                AttrValue::Number(n) => format!("{k}={n}"),
                AttrValue::Text(t) => format!("{k}={t}"),
            })
            .collect();
        let status = status_name(r.status);
        let line = [
            cell(&r.body),
            r.in_study.to_string(),
            cell(&r.tags.join(" ")),
            cell(&attributes.join("; ")),
            number(r.conductivity_w_per_m_k),
            number(Some(r.power_w)),
            number(r.limit_k),
            number(r.result.as_ref().map(|x| x.max_temperature_k)),
            number(r.result.as_ref().map(|x| x.mean_temperature_k)),
            number(r.margin_k),
            status.to_string(),
        ];
        out.push_str(&line.join(","));
        out.push('\n');
    }
    out.push_str(&format!("total,,,,,{},,,,,\n", number(Some(table.total_power_w))));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{BodyHeat, FaceHeat, FaceRef, ModelHeat};
    use crate::parts::Retag;

    fn study_with_run(report: Option<&str>) -> Study {
        let mut project = crate::Project::default();
        let mut study = project.create("flow-channel", "Board").unwrap();
        let material = |body: &str, power: f64, limit: Option<f64>| BodyHeat {
            body: body.into(),
            conductivity_w_per_m_k: 200.0,
            power_w: power,
            density_kg_per_m3: None,
            specific_heat_j_per_kg_k: None,
            max_temperature_k: limit,
        };
        study.model = Some(ModelLink {
            bodies: vec!["U1".into(), "R 2".into()],
            heat: Some(ModelHeat {
                materials: vec![material("U1", 3.0, Some(358.0)), material("R 2", 0.5, Some(340.0))],
                face_sources: vec![FaceHeat { name: "pad".into(), face: FaceRef { body: "U1".into(), face: "top".into() }, power_w: 1.5 }],
                ..ModelHeat::default()
            }),
            ..ModelLink::default()
        });
        if let Some(report) = report {
            study.runs.push(RunEvidence {
                id: "run-1".into(),
                study_revision: study.revision,
                case_hash: String::new(),
                input_case: serde_json::Value::Null,
                manifest: serde_json::Value::Null,
                artifacts: vec![crate::Artifact { name: "thermal_report.json".into(), text: report.into() }],
            });
        }
        study
    }

    const REPORT: &str = r#"{"components":[
        {"name":"U1","power_w":4.5,"min_temperature_k":300,"mean_temperature_k":330,"max_temperature_k":350,"limit_k":null,"margin_k":null,"status":"no_limit"},
        {"name":"R_2","power_w":0.5,"min_temperature_k":300,"mean_temperature_k":335,"max_temperature_k":345,"limit_k":null,"margin_k":null,"status":"no_limit"}]}"#;

    #[test]
    fn rows_join_setup_and_results_and_find_the_worst_margin() {
        let study = study_with_run(Some(REPORT));
        let mut registry = Registry::default();
        registry.retag(&["U1".into()], &Retag { add_tags: vec!["power".into()], ..Retag::default() }).ok();
        let bodies = ["U1".to_string(), "R 2".to_string(), "Case".to_string()];
        let table = build(&study, &registry, &bodies, &bodies[..2], None).unwrap();
        assert_eq!(table.run.as_deref(), Some("run-1"));
        // 3 W in U1 plus 1.5 W on its face, and 0.5 W in R 2; the case is not in the study.
        assert_eq!(table.total_power_w, 5.0);
        assert_eq!(table.rows[0].power_w, 4.5);
        assert_eq!(table.rows[0].tags, vec!["power".to_string()]);
        assert_eq!(table.rows[0].margin_k, Some(8.0));
        // R 2 is reported as R_2 and runs 5 K over its 340 K limit.
        assert_eq!(table.rows[1].margin_k, Some(-5.0));
        assert_eq!(table.rows[1].status, ComponentStatus::Fail);
        assert!(!table.rows[2].in_study && table.rows[2].result.is_none());
        assert_eq!(table.worst.as_deref(), Some("R 2"));
        let csv = to_csv(&table);
        assert!(csv.lines().nth(1).is_some_and(|l| l.starts_with("U1,true,power,,200,4.5,358,350,330,8,pass")), "{csv}");
        assert!(csv.ends_with("total,,,,,5,,,,,\n"), "{csv}");
    }

    #[test]
    fn before_a_run_the_table_holds_the_setup_alone() {
        let study = study_with_run(None);
        let bodies = ["U1".to_string()];
        let table = build(&study, &Registry::default(), &bodies, &bodies, None).unwrap();
        assert_eq!((table.run.clone(), table.rows[0].margin_k, table.rows[0].status), (None, None, ComponentStatus::NoLimit));
        assert!(build(&study, &Registry::default(), &bodies, &bodies, Some("nope")).is_err());
    }

    #[test]
    fn csv_cells_cannot_become_formulas_or_break_lines() {
        assert_eq!(cell("=HYPERLINK(\"x\")"), "\"'=HYPERLINK(\"\"x\"\")\"");
        assert_eq!(cell("a,b"), "\"a,b\"");
        assert_eq!(cell("-3.5"), "-3.5");
        assert_eq!(cell("line\nbreak"), "\"line\nbreak\"");
    }
}
