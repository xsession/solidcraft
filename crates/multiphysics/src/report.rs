//! Batch report: one post setup, the component table, applied to several runs of a study, so
//! design variants or operating points are read side by side. A summary per run (total power,
//! the hottest part, the least margin, how many parts pass, warn or fail) and every part's
//! numbers per run, as CSV and as a self-contained HTML page.
//!
//! Portable: the engine says which bodies the study holds and writes the files.
use schemars::JsonSchema;
use serde::Serialize;

use crate::component_table::{self, ComponentTable, cell, number, status_name};
use crate::parts::Registry;
use crate::results::ComponentStatus;
use crate::{Error, Result, Study};

/// Most runs one report reads.
pub const MAX_REPORT_RUNS: usize = 64;

/// The colours a report draws verdicts in, in its pictures: blue, yellow and vermilion, which
/// every eye tells apart (anyone may read a report), and grey, drawn see-through, for a part
/// without a limit or a result.
pub fn verdict_rgb(status: Option<ComponentStatus>) -> [u8; 3] {
    let [pass, warn, fail] = solidcraft_color::VERDICTS_FOR_EVERY_EYE;
    match status {
        Some(ComponentStatus::Pass) => pass,
        Some(ComponentStatus::Warn) => warn,
        Some(ComponentStatus::Fail) => fail,
        Some(ComponentStatus::NoLimit) | None => [150, 150, 150],
    }
}

/// A verdict in words and a mark after a cell's number, so the page does not lean on colour
/// alone; nothing for a part without a limit.
fn verdict_mark(status: ComponentStatus) -> &'static str {
    match status {
        ComponentStatus::Pass => " <span class=\"v\">✓ pass</span>",
        ComponentStatus::Warn => " <span class=\"v\">! warn</span>",
        ComponentStatus::Fail => " <span class=\"v\">✗ fail</span>",
        ComponentStatus::NoLimit => "",
    }
}

/// Each part's verdict in one run of the report (none where the run reports nothing on it).
pub fn verdicts(table: &ComponentTable) -> Vec<(String, Option<ComponentStatus>)> {
    table.rows.iter().map(|r| (r.body.clone(), r.result.as_ref().map(|_| r.status))).collect()
}

/// One run's line of the report.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct RunSummary {
    pub run: String,
    /// The run was made from an earlier revision of the study.
    pub stale: bool,
    pub total_power_w: f64,
    /// The part with the highest temperature, and that temperature (°C).
    pub hottest: Option<String>,
    pub hottest_c: Option<f64>,
    /// The part with the least margin to its limit, and that margin (K).
    pub worst: Option<String>,
    pub least_margin_k: Option<f64>,
    pub pass: usize,
    pub warn: usize,
    pub fail: usize,
    /// Parts of the study the run reports no temperature for.
    pub no_result: usize,
}

#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct BatchReport {
    pub study: String,
    pub runs: Vec<RunSummary>,
    /// Each run's component table, in the order of `runs`.
    #[serde(skip)]
    pub tables: Vec<ComponentTable>,
}

fn celsius(k: f64) -> f64 {
    k - 273.15
}

/// The report over `run_ids` (every run with a thermal report, oldest first, when empty), of
/// the study's parts among `in_study`.
pub fn build(study: &Study, registry: &Registry, in_study: &[String], run_ids: &[String]) -> Result<BatchReport> {
    let runs: Vec<String> = if run_ids.is_empty() {
        study.runs.iter().filter(|r| r.artifacts.iter().any(|a| a.name == "thermal_report.json")).map(|r| r.id.clone()).collect()
    } else {
        run_ids.to_vec()
    };
    if runs.is_empty() {
        return Err(Error::Invalid("the study has no run with a thermal report to report on".into()));
    }
    if runs.len() > MAX_REPORT_RUNS {
        return Err(Error::Invalid(format!("a report reads at most {MAX_REPORT_RUNS} runs")));
    }
    let mut summaries = Vec::with_capacity(runs.len());
    let mut tables = Vec::with_capacity(runs.len());
    for run in &runs {
        let table = component_table::build(study, registry, in_study, in_study, Some(run))?;
        let hottest =
            table.rows.iter().filter_map(|r| r.result.as_ref().map(|x| (r.body.clone(), x.max_temperature_k))).max_by(|a, b| a.1.total_cmp(&b.1));
        let count = |s: ComponentStatus| table.rows.iter().filter(|r| r.result.is_some() && r.status == s).count();
        summaries.push(RunSummary {
            run: run.clone(),
            stale: table.stale,
            total_power_w: table.total_power_w,
            hottest_c: hottest.as_ref().map(|h| celsius(h.1)),
            hottest: hottest.map(|h| h.0),
            least_margin_k: table.worst.as_ref().and_then(|w| table.rows.iter().find(|r| &r.body == w)).and_then(|r| r.margin_k),
            worst: table.worst.clone(),
            pass: count(ComponentStatus::Pass),
            warn: count(ComponentStatus::Warn),
            fail: count(ComponentStatus::Fail),
            no_result: table.rows.iter().filter(|r| r.result.is_none()).count(),
        });
        tables.push(table);
    }
    Ok(BatchReport { study: study.name.clone(), runs: summaries, tables })
}

/// Every part's numbers per run, one line each: run, part, power, limit, temperatures (°C),
/// margin (K) and verdict.
pub fn to_csv(report: &BatchReport) -> String {
    let mut out = String::from("run,body,power_w,limit_c,max_temperature_c,mean_temperature_c,margin_k,status\n");
    for (summary, table) in report.runs.iter().zip(&report.tables) {
        for r in &table.rows {
            let status = if r.result.is_some() { status_name(r.status) } else { "" };
            let line = [
                cell(&summary.run),
                cell(&r.body),
                number(Some(r.power_w)),
                number(r.limit_k.map(celsius)),
                number(r.result.as_ref().map(|x| celsius(x.max_temperature_k))),
                number(r.result.as_ref().map(|x| celsius(x.mean_temperature_k))),
                number(r.margin_k),
                status.to_string(),
            ];
            out.push_str(&line.join(","));
            out.push('\n');
        }
    }
    out
}

fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            c => out.push(c),
        }
    }
    out
}

fn fixed(v: Option<f64>, places: usize) -> String {
    v.filter(|v| v.is_finite()).map_or_else(|| "–".into(), |v| format!("{v:.places$}"))
}

/// A self-contained page: the runs' summary, then each part's highest temperature per run,
/// marked by its verdict, and each run's picture where `pictures` names one (a file beside the
/// page). No scripts and nothing fetched from elsewhere; every text is escaped.
pub fn to_html(report: &BatchReport, pictures: &[Option<String>]) -> String {
    let title = format!("{} — batch report", report.study);
    let mut out = String::new();
    out.push_str("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">");
    out.push_str(&format!("<title>{}</title>", escape(&title)));
    out.push_str(
        "<style>body{font-family:sans-serif;margin:2em}table{border-collapse:collapse;margin:1em 0}\
         th,td{border:1px solid #999;padding:.3em .6em;text-align:right}th:first-child,td:first-child{text-align:left}\
         .pass{background:#d6e6f7}.warn{background:#fbf6c8}.fail{background:#f8dcc8}.stale{color:#8a5a00}.v{font-size:.85em}</style></head><body>",
    );
    out.push_str(&format!("<h1>{}</h1>", escape(&title)));
    out.push_str("<h2>Runs</h2><table><tr><th>Run</th><th>Power (W)</th><th>Hottest part</th><th>Max (°C)</th><th>Least margin</th><th>Margin (K)</th><th>Pass</th><th>Warn</th><th>Fail</th><th>No result</th></tr>");
    for s in &report.runs {
        let run = if s.stale { format!("{} <span class=\"stale\">(older setup)</span>", escape(&s.run)) } else { escape(&s.run) };
        out.push_str(&format!(
            "<tr><td>{run}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            fixed(Some(s.total_power_w), 2),
            escape(s.hottest.as_deref().unwrap_or("–")),
            fixed(s.hottest_c, 1),
            escape(s.worst.as_deref().unwrap_or("–")),
            fixed(s.least_margin_k, 1),
            s.pass,
            s.warn,
            s.fail,
            s.no_result,
        ));
    }
    out.push_str("</table><h2>Highest temperature per part (°C)</h2><table><tr><th>Part</th>");
    for s in &report.runs {
        out.push_str(&format!("<th>{}</th>", escape(&s.run)));
    }
    out.push_str("</tr>");
    let parts: Vec<&str> = report.tables.first().map(|t| t.rows.iter().map(|r| r.body.as_str()).collect()).unwrap_or_default();
    for part in parts {
        out.push_str(&format!("<tr><td>{}</td>", escape(part)));
        for table in &report.tables {
            let row = table.rows.iter().find(|r| r.body == part);
            let result = row.and_then(|r| r.result.as_ref().map(|x| (celsius(x.max_temperature_k), r.status)));
            match result {
                Some((max, status)) => {
                    out.push_str(&format!("<td class=\"{}\">{}{}</td>", status_name(status), fixed(Some(max), 1), verdict_mark(status)))
                }
                None => out.push_str("<td>–</td>"),
            }
        }
        out.push_str("</tr>");
    }
    out.push_str("</table>");
    if pictures.iter().any(Option::is_some) {
        out.push_str("<h2>Parts by verdict</h2><p>Blue: pass, yellow: warn, vermilion: fail; see-through grey: no limit or no result.</p>");
        for (s, picture) in report.runs.iter().zip(pictures) {
            if let Some(file) = picture {
                out.push_str(&format!(
                    "<figure><img src=\"{}\" alt=\"{}\" width=\"640\"><figcaption>{}</figcaption></figure>",
                    escape(file),
                    escape(&s.run),
                    escape(&s.run)
                ));
            }
        }
    }
    out.push_str("</body></html>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{BodyHeat, ModelHeat, ModelLink};
    use crate::{Artifact, RunEvidence};

    fn report(u1: f64, r2: f64) -> String {
        format!(
            r#"{{"components":[
            {{"name":"U1","power_w":3,"min_temperature_k":300,"mean_temperature_k":320,"max_temperature_k":{u1},"limit_k":null,"margin_k":null,"status":"no_limit"}},
            {{"name":"R_2","power_w":0.5,"min_temperature_k":300,"mean_temperature_k":310,"max_temperature_k":{r2},"limit_k":null,"margin_k":null,"status":"no_limit"}}]}}"#
        )
    }

    fn study() -> Study {
        let mut project = crate::Project::default();
        let mut study = project.create("flow-channel", "Board <A&B>").unwrap();
        let material = |body: &str, power: f64, limit: f64| BodyHeat {
            body: body.into(),
            conductivity_w_per_m_k: 200.0,
            power_w: power,
            density_kg_per_m3: None,
            specific_heat_j_per_kg_k: None,
            max_temperature_k: Some(limit),
        };
        study.model = Some(ModelLink {
            bodies: vec!["U1".into(), "R 2".into()],
            heat: Some(ModelHeat { materials: vec![material("U1", 3.0, 358.0), material("R 2", 0.5, 340.0)], ..ModelHeat::default() }),
            ..ModelLink::default()
        });
        for (id, u1, r2) in [("slow fan", 350.0, 345.0), ("fast fan", 330.0, 320.0)] {
            study.runs.push(RunEvidence {
                id: id.into(),
                study_revision: study.revision,
                case_hash: String::new(),
                input_case: serde_json::Value::Null,
                manifest: serde_json::Value::Null,
                artifacts: vec![Artifact { name: "thermal_report.json".into(), text: report(u1, r2) }],
            });
        }
        study
    }

    #[test]
    fn a_report_reads_every_run_side_by_side() {
        let study = study();
        let parts = ["U1".to_string(), "R 2".to_string()];
        let report = build(&study, &Registry::default(), &parts, &[]).unwrap();
        let runs: Vec<&str> = report.runs.iter().map(|r| r.run.as_str()).collect();
        assert_eq!(runs, ["slow fan", "fast fan"], "oldest first");
        let slow = &report.runs[0];
        assert_eq!((slow.hottest.as_deref(), slow.worst.as_deref()), (Some("U1"), Some("R 2")));
        assert!((slow.hottest_c.unwrap() - 76.85).abs() < 1e-9);
        assert_eq!(slow.least_margin_k, Some(-5.0));
        assert_eq!((slow.pass, slow.warn, slow.fail, slow.no_result), (1, 0, 1, 0));
        assert_eq!((report.runs[1].pass, report.runs[1].fail), (2, 0), "the fast fan cools both");
        let csv = to_csv(&report);
        assert_eq!(csv.lines().count(), 5, "a header and two parts for each of two runs: {csv}");
        assert!(csv.contains("\nslow fan,R 2,0.5,"), "{csv}");
        assert!(csv.lines().nth(2).is_some_and(|l| l.ends_with(",-5,fail")), "{csv}");
        let html = to_html(&report, &[]);
        assert!(html.contains("<title>Board &lt;A&amp;B&gt; — batch report</title>"), "escaped: {html}");
        assert!(!html.contains("<A&B>") && !html.contains("<script"));
        assert!(html.contains("<td class=\"fail\">71.9 <span class=\"v\">✗ fail</span></td>"), "{html}");
        assert!(!html.contains("<img"), "no pictures");
        let pictured = to_html(&report, &[Some("run-1.png".into()), None]);
        assert!(pictured.contains("<img src=\"run-1.png\" alt=\"slow fan\""), "{pictured}");
        assert_eq!(pictured.matches("<img").count(), 1);
        let shown: Vec<_> = verdicts(&report.tables[0]).into_iter().map(|(_, v)| v).collect();
        assert_eq!(shown, [Some(ComponentStatus::Pass), Some(ComponentStatus::Fail)]);
        assert_ne!(verdict_rgb(Some(ComponentStatus::Fail)), verdict_rgb(Some(ComponentStatus::Pass)));
        // One run by name, and refusals.
        let one = build(&study, &Registry::default(), &parts, &["fast fan".to_string()]).unwrap();
        assert_eq!(one.runs.len(), 1);
        assert!(build(&study, &Registry::default(), &parts, &["nope".to_string()]).is_err());
        let many = vec!["fast fan".to_string(); MAX_REPORT_RUNS + 1];
        assert!(build(&study, &Registry::default(), &parts, &many).is_err());
        let mut bare = study.clone();
        bare.runs.clear();
        assert!(build(&bare, &Registry::default(), &parts, &[]).is_err());
    }
}
