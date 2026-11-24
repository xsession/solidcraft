//! Find: a stack of rules over a study's parts, on their setup (the part registry's tags, names
//! and attributes) and on a run's results (power, temperatures, margin, verdict), as in "the
//! parts tagged ic hotter than 75 °C". A part is found when it meets every rule.
//!
//! Portable: the engine supplies the component table (setup and results per part) and the
//! design's bodies as the registry sees them.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::component_table::{ComponentTable, TableRow};
use crate::parts::{Compare, PartFacts, Registry, Rule};
use crate::results::ComponentStatus;
use crate::{Error, Result};

/// A quantity of the component table a result rule compares.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Quantity {
    /// Heat released in the part and on its faces, W (setup; no run needed).
    PowerW,
    /// The part's temperature limit, °C (setup).
    LimitC,
    /// The run's highest temperature in the part, °C.
    MaxTemperatureC,
    /// The run's mean temperature in the part, °C.
    MeanTemperatureC,
    /// The limit less the run's highest temperature, K.
    MarginK,
}

/// One rule of a find.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(untagged)]
pub enum FindRule {
    /// A rule on the part's setup, as part selections use them.
    Setup(Rule),
    /// A rule on the study's numbers.
    Result(ResultRule),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "rule", rename_all = "snake_case", deny_unknown_fields)]
pub enum ResultRule {
    /// A quantity compared with a number; a part without that quantity does not meet it.
    Result { quantity: Quantity, compare: Compare, value: f64 },
    /// The run's verdict for the part against its limit.
    Verdict { status: ComponentStatus },
}

/// Most rules a find stacks.
pub const MAX_FIND_RULES: usize = 16;

fn quantity(row: &TableRow, q: Quantity) -> Option<f64> {
    let c = |k: f64| k - 273.15;
    match q {
        Quantity::PowerW => Some(row.power_w),
        Quantity::LimitC => row.limit_k.map(c),
        Quantity::MaxTemperatureC => row.result.as_ref().map(|r| c(r.max_temperature_k)),
        Quantity::MeanTemperatureC => row.result.as_ref().map(|r| c(r.mean_temperature_k)),
        Quantity::MarginK => row.margin_k,
    }
}

fn holds(have: f64, how: Compare, want: f64) -> bool {
    match how {
        Compare::Exists => true,
        Compare::Eq => (have - want).abs() <= 1e-9 * want.abs().max(1.0),
        Compare::Ne => (have - want).abs() > 1e-9 * want.abs().max(1.0),
        Compare::Lt => have < want,
        Compare::Le => have <= want,
        Compare::Gt => have > want,
        Compare::Ge => have >= want,
    }
}

/// Whether the find reads the built-in measured attributes (so the caller measures bodies).
pub fn needs_measure(rules: &[FindRule], registry: &Registry) -> bool {
    rules.iter().any(|r| matches!(r, FindRule::Setup(rule) if registry.rule_needs_measure(rule)))
}

/// Whether any rule reads a run's results.
pub fn needs_run(rules: &[FindRule]) -> bool {
    rules.iter().any(|r| {
        matches!(
            r,
            FindRule::Result(ResultRule::Verdict { .. })
                | FindRule::Result(ResultRule::Result { quantity: Quantity::MaxTemperatureC | Quantity::MeanTemperatureC | Quantity::MarginK, .. })
        )
    })
}

/// The rows of `table` that meet every rule, in its order. `bodies` are the design's bodies as
/// the registry sees them (measured where [`needs_measure`] says).
pub fn find(rules: &[FindRule], registry: &Registry, table: &ComponentTable, bodies: &[PartFacts]) -> Result<Vec<String>> {
    if rules.is_empty() || rules.len() > MAX_FIND_RULES {
        return Err(Error::Invalid(format!("a find has 1 to {MAX_FIND_RULES} rules")));
    }
    for rule in rules {
        if let FindRule::Result(ResultRule::Result { compare, value, .. }) = rule
            && *compare != Compare::Exists
            && !value.is_finite()
        {
            return Err(Error::Invalid("a result rule compares with a finite number".into()));
        }
    }
    let mut found = Vec::new();
    for row in &table.rows {
        let Some(facts) = bodies.iter().find(|b| b.name == row.body) else { continue };
        let mut all = true;
        for rule in rules {
            let meets = match rule {
                FindRule::Setup(rule) => registry.meets_rule(rule, facts, bodies)?,
                FindRule::Result(ResultRule::Result { quantity: q, compare, value }) => {
                    quantity(row, *q).is_some_and(|have| holds(have, *compare, *value))
                }
                FindRule::Result(ResultRule::Verdict { status }) => row.result.is_some() && row.status == *status,
            };
            if !meets {
                all = false;
                break;
            }
        }
        if all {
            found.push(row.body.clone());
        }
    }
    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{BodyHeat, ModelHeat, ModelLink};
    use crate::parts::Retag;

    const REPORT: &str = r#"{"components":[
        {"name":"U1","power_w":2,"min_temperature_k":300,"mean_temperature_k":330,"max_temperature_k":355,"limit_k":358,"margin_k":3,"status":"warn"},
        {"name":"U2","power_w":1,"min_temperature_k":300,"mean_temperature_k":310,"max_temperature_k":320,"limit_k":358,"margin_k":38,"status":"pass"}]}"#;

    fn setup() -> (Registry, ComponentTable, Vec<PartFacts>) {
        let mut project = crate::Project::default();
        let mut study = project.create("flow-channel", "Board").expect("study");
        let material = |body: &str, power: f64| BodyHeat {
            body: body.into(),
            conductivity_w_per_m_k: 1.0,
            power_w: power,
            density_kg_per_m3: None,
            specific_heat_j_per_kg_k: None,
            max_temperature_k: Some(358.0),
        };
        study.model = Some(ModelLink {
            heat: Some(ModelHeat { materials: vec![material("U1", 2.0), material("U2", 1.0)], ..ModelHeat::default() }),
            ..ModelLink::default()
        });
        study.runs.push(crate::RunEvidence {
            id: "run".into(),
            study_revision: study.revision,
            case_hash: String::new(),
            input_case: serde_json::Value::Null,
            manifest: serde_json::Value::Null,
            artifacts: vec![crate::Artifact { name: "thermal_report.json".into(), text: REPORT.into() }],
        });
        let mut registry = Registry::default();
        let tag = |t: &str| Retag { add_tags: vec![t.into()], remove_tags: vec![], attributes: Default::default() };
        registry.retag(&["U1".into(), "U2".into()], &tag("ic")).expect("tag");
        let bodies: Vec<String> = ["U1", "U2", "Board"].map(String::from).to_vec();
        let table = crate::component_table::build(&study, &registry, &bodies, &bodies, None).expect("table");
        let facts = bodies.iter().map(|b| PartFacts { name: b.clone(), ..PartFacts::default() }).collect();
        (registry, table, facts)
    }

    fn rules(json: serde_json::Value) -> Vec<FindRule> {
        serde_json::from_value(json).expect("rules")
    }

    #[test]
    fn setup_and_result_rules_stack() {
        let (registry, table, facts) = setup();
        // Tagged ic and hotter than 75 °C: only U1 (81.85 °C); U2 runs at 46.85 °C.
        let hot = rules(serde_json::json!([{"rule":"tag","tag":"ic"},{"rule":"result","quantity":"max_temperature_c","compare":"gt","value":75.0}]));
        assert_eq!(find(&hot, &registry, &table, &facts).expect("find"), ["U1"]);
        assert!(needs_run(&hot) && !needs_measure(&hot, &registry));
        let warned = rules(serde_json::json!([{"rule":"verdict","status":"warn"}]));
        assert_eq!(find(&warned, &registry, &table, &facts).expect("find"), ["U1"]);
        // Power is setup: the board, with none, is not over 0.5 W.
        let powered = rules(serde_json::json!([{"rule":"result","quantity":"power_w","compare":"ge","value":0.5}]));
        assert_eq!(find(&powered, &registry, &table, &facts).expect("find"), ["U1", "U2"]);
        assert!(!needs_run(&powered));
        let untagged = rules(serde_json::json!([{"rule":"not_tag","tag":"ic"}]));
        assert_eq!(find(&untagged, &registry, &table, &facts).expect("find"), ["Board"]);
    }

    #[test]
    fn hostile_finds_are_refused() {
        let (registry, table, facts) = setup();
        assert!(find(&[], &registry, &table, &facts).is_err());
        let nan = vec![FindRule::Result(ResultRule::Result { quantity: Quantity::MarginK, compare: Compare::Lt, value: f64::NAN })];
        assert!(find(&nan, &registry, &table, &facts).is_err());
        let many = vec![FindRule::Setup(Rule::Tag { tag: "ic".into() }); MAX_FIND_RULES + 1];
        assert!(find(&many, &registry, &table, &facts).is_err());
        let missing = rules(serde_json::json!([{"rule":"in_selection","name":"nope"}]));
        assert!(find(&missing, &registry, &table, &facts).is_err());
    }
}
