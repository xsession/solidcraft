//! Variant comparison: two runs' component tables side by side, part by part. The runs may be of
//! one study (two operating points) or of two studies of the model (two design variants, such as
//! a study and its clone or two points of a sweep). Each part's highest temperature in both, the
//! change between them, both margins and verdicts, and whether the verdict got better or worse;
//! the parts that changed most come first.
//!
//! Portable: the engine builds each run's component table and says which study and run it is.
use schemars::JsonSchema;
use serde::Serialize;

use crate::component_table::{ComponentTable, cell, number, status_name};
use crate::results::ComponentStatus;

/// A rise or fall smaller than this (K) counts as no change: a run's temperatures are not more
/// precise than that.
pub const SAME_K: f64 = 0.05;

/// One side of the comparison.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Side {
    pub study: String,
    pub run: String,
    /// The run was made from an earlier revision of its study.
    pub stale: bool,
    pub total_power_w: f64,
}

/// How a part's verdict moved from the first run to the second.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Change {
    /// Pass to warn or fail, warn to fail.
    Worse,
    /// Fail to warn or pass, warn to pass.
    Better,
    /// The same verdict, or one without a limit to judge by.
    Same,
    /// Only one of the runs reports the part.
    OneSide,
}

/// One part in both runs. Temperatures in °C, changes and margins in K.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct DeltaRow {
    pub body: String,
    pub first_c: Option<f64>,
    pub second_c: Option<f64>,
    /// The second run's highest temperature less the first's.
    pub delta_k: Option<f64>,
    pub first_margin_k: Option<f64>,
    pub second_margin_k: Option<f64>,
    pub first_status: Option<ComponentStatus>,
    pub second_status: Option<ComponentStatus>,
    pub change: Change,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Comparison {
    pub first: Side,
    pub second: Side,
    /// Parts in either run's study, the largest rise first, then the largest fall; parts only
    /// one run reports last.
    pub rows: Vec<DeltaRow>,
    /// Parts at least [`SAME_K`] warmer, and cooler, in the second run.
    pub warmer: usize,
    pub cooler: usize,
    /// Parts whose verdict got worse, and better.
    pub worse: usize,
    pub better: usize,
    /// The second run's total power less the first's, W.
    pub power_delta_w: f64,
}

fn celsius(k: f64) -> f64 {
    k - 273.15
}

/// How bad a verdict is, where it can be ranked.
fn severity(status: ComponentStatus) -> Option<u8> {
    match status {
        ComponentStatus::Pass => Some(0),
        ComponentStatus::Warn => Some(1),
        ComponentStatus::Fail => Some(2),
        ComponentStatus::NoLimit => None,
    }
}

/// `first` and `second` part by part: the parts of each table's study, matched by name.
pub fn compare(first: &ComponentTable, second: &ComponentTable) -> Comparison {
    let side = |table: &ComponentTable| Side {
        study: table.study.clone(),
        run: table.run.clone().unwrap_or_default(),
        stale: table.stale,
        total_power_w: table.total_power_w,
    };
    let reported = |table: &ComponentTable, body: &str| table.rows.iter().find(|r| r.in_study && r.body == body && r.result.is_some()).cloned();
    let mut names: Vec<&str> = Vec::new();
    for row in first.rows.iter().chain(&second.rows).filter(|r| r.in_study) {
        if !names.contains(&row.body.as_str()) {
            names.push(&row.body);
        }
    }
    let mut rows: Vec<DeltaRow> = names
        .into_iter()
        .map(|body| {
            let (a, b) = (reported(first, body), reported(second, body));
            let temperature = |row: &Option<crate::component_table::TableRow>| {
                row.as_ref().and_then(|r| r.result.as_ref()).map(|x| celsius(x.max_temperature_k)).filter(|t| t.is_finite())
            };
            let (first_c, second_c) = (temperature(&a), temperature(&b));
            let (first_status, second_status) = (a.as_ref().map(|r| r.status), b.as_ref().map(|r| r.status));
            let change = match (first_status, second_status) {
                (Some(x), Some(y)) => match (severity(x), severity(y)) {
                    (Some(x), Some(y)) if y > x => Change::Worse,
                    (Some(x), Some(y)) if y < x => Change::Better,
                    _ => Change::Same,
                },
                _ => Change::OneSide,
            };
            DeltaRow {
                body: body.to_string(),
                first_c,
                second_c,
                delta_k: first_c.zip(second_c).map(|(a, b)| b - a),
                first_margin_k: a.as_ref().and_then(|r| r.margin_k),
                second_margin_k: b.as_ref().and_then(|r| r.margin_k),
                first_status,
                second_status,
                change,
            }
        })
        .collect();
    // The largest rise first, down to the largest fall; parts one run leaves out at the end.
    rows.sort_by(|a, b| match (a.delta_k, b.delta_k) {
        (Some(x), Some(y)) => y.total_cmp(&x),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });
    let count = |keep: &dyn Fn(&DeltaRow) -> bool| rows.iter().filter(|r| keep(r)).count();
    Comparison {
        first: side(first),
        second: side(second),
        warmer: count(&|r| r.delta_k.is_some_and(|d| d >= SAME_K)),
        cooler: count(&|r| r.delta_k.is_some_and(|d| d <= -SAME_K)),
        worse: count(&|r| r.change == Change::Worse),
        better: count(&|r| r.change == Change::Better),
        power_delta_w: second.total_power_w - first.total_power_w,
        rows,
    }
}

/// The comparison as CSV, a line per part in its order.
pub fn to_csv(comparison: &Comparison) -> String {
    let mut out = String::from("body,first_c,second_c,delta_k,first_margin_k,second_margin_k,first_status,second_status,change\n");
    let status = |s: Option<ComponentStatus>| s.map(status_name).unwrap_or_default();
    // To a millionth of a degree: no noise from the kelvin to Celsius conversion.
    let number = |v: Option<f64>| number(v.map(|x| (x * 1e6).round() / 1e6));
    for row in &comparison.rows {
        let change = match row.change {
            Change::Worse => "worse",
            Change::Better => "better",
            Change::Same => "same",
            Change::OneSide => "one_side",
        };
        out.push_str(&format!(
            "{},{},{},{},{},{},{},{},{}\n",
            cell(&row.body),
            number(row.first_c),
            number(row.second_c),
            number(row.delta_k),
            number(row.first_margin_k),
            number(row.second_margin_k),
            status(row.first_status),
            status(row.second_status),
            change
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::component_table::TableRow;
    use crate::results::ComponentRow;

    fn row(body: &str, max_c: Option<f64>, limit_c: Option<f64>, status: ComponentStatus) -> TableRow {
        let limit_k = limit_c.map(|c| c + 273.15);
        let result = max_c.map(|c| ComponentRow {
            name: body.into(),
            power_w: 1.0,
            min_temperature_k: 300.0,
            mean_temperature_k: 310.0,
            max_temperature_k: c + 273.15,
            limit_k,
            margin_k: limit_k.map(|l| l - (c + 273.15)),
            status,
        });
        TableRow {
            body: body.into(),
            tags: vec![],
            attributes: Default::default(),
            in_study: true,
            conductivity_w_per_m_k: None,
            power_w: 1.0,
            limit_k,
            margin_k: result.as_ref().and_then(|r| r.margin_k),
            result,
            status,
        }
    }

    fn table(study: &str, run: &str, rows: Vec<TableRow>) -> ComponentTable {
        let total_power_w = rows.iter().filter(|r| r.in_study).map(|r| r.power_w).sum();
        ComponentTable { study: study.into(), run: Some(run.into()), stale: false, rows, total_power_w, worst: None }
    }

    #[test]
    fn two_runs_compare_part_by_part_with_the_largest_rise_first() {
        use ComponentStatus::*;
        let first = table(
            "Board",
            "fast",
            vec![
                row("U1", Some(74.6), Some(85.0), Pass),
                row("U2", Some(77.1), Some(85.0), Pass),
                row("Q1", Some(72.1), Some(80.0), Pass),
                row("J1", Some(53.3), Some(70.0), Pass),
                row("PCB1", Some(70.0), None, NoLimit),
            ],
        );
        let mut slow_rows = vec![
            row("U1", Some(83.0), Some(85.0), Warn),
            row("U2", Some(86.7), Some(85.0), Fail),
            row("Q1", Some(82.0), Some(80.0), Fail),
            row("J1", Some(53.32), Some(70.0), Pass),
            row("PCB1", Some(83.2), None, NoLimit),
            row("C1", None, None, NoLimit),
        ];
        slow_rows.push(TableRow { in_study: false, ..row("Duct", Some(56.0), None, NoLimit) });
        let second = table("Board", "slow", slow_rows);
        let c = compare(&first, &second);
        let order: Vec<&str> = c.rows.iter().map(|r| r.body.as_str()).collect();
        // C1 is in the study but neither run reports it; the duct is not in the study.
        assert_eq!(order, ["PCB1", "Q1", "U2", "U1", "J1", "C1"]);
        let q1 = &c.rows[1];
        assert!((q1.delta_k.unwrap() - 9.9).abs() < 1e-9);
        assert_eq!((q1.first_status, q1.second_status, q1.change), (Some(Pass), Some(Fail), Change::Worse));
        assert!((q1.second_margin_k.unwrap() + 2.0).abs() < 1e-9);
        assert_eq!(c.rows[0].change, Change::Same, "no limit, no verdict to compare");
        assert_eq!(c.rows[5].change, Change::OneSide);
        // J1 rose 0.02 K: no change worth the name.
        assert_eq!((c.warmer, c.cooler, c.worse, c.better), (4, 0, 3, 0));
        assert!((c.power_delta_w - 1.0).abs() < 1e-9, "C1's watt is the second table's only");
        // The other way round, the verdicts get better and the parts cooler.
        let back = compare(&second, &first);
        assert_eq!((back.warmer, back.cooler, back.worse, back.better), (0, 4, 0, 3));
        assert_eq!(back.rows[0].body, "J1", "the least fall first");
        let csv = to_csv(&c);
        assert!(csv.starts_with("body,first_c,second_c,delta_k,"));
        assert!(csv.lines().any(|l| l.starts_with("Q1,72.1,82,") && l.ends_with(",pass,fail,worse")), "{csv}");
        assert_eq!(csv.lines().count(), 7);
    }

    #[test]
    fn parts_of_two_studies_match_by_name() {
        use ComponentStatus::*;
        let first = table("Board", "a", vec![row("U1", Some(60.0), Some(85.0), Pass), row("R1", Some(40.0), None, NoLimit)]);
        let second = table("Board (2)", "b", vec![row("U1", Some(58.5), Some(85.0), Pass), row("U9", Some(70.0), Some(85.0), Pass)]);
        let c = compare(&first, &second);
        assert_eq!((c.first.study.as_str(), c.second.study.as_str()), ("Board", "Board (2)"));
        let order: Vec<(&str, Change)> = c.rows.iter().map(|r| (r.body.as_str(), r.change)).collect();
        assert_eq!(order, [("U1", Change::Same), ("R1", Change::OneSide), ("U9", Change::OneSide)]);
        assert!((c.rows[0].delta_k.unwrap() + 1.5).abs() < 1e-9);
        assert_eq!((c.warmer, c.cooler), (0, 1));
    }
}
