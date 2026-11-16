//! ASCII SPICE RAW interchange (writer, reader, and transient/AC dataset builders). Ported from
//! upstream `cfd/circuit/raw` (MIT, cfd_solvers contributors).
use super::spice::{AcSolution, Complex, TransientPoint};
use crate::workflow::text::{format_g17, getline_split, strtod_prefix};
use crate::{require, Result, SolverError};

#[derive(Clone, Debug, PartialEq)]
pub struct RawTrace {
    pub name: String,
    pub kind: String,
    pub values: Vec<Complex>,
}

impl RawTrace {
    pub fn voltage(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            kind: "voltage".into(),
            values: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpiceRawDataset {
    pub title: String,
    pub plot_name: String,
    pub scale_name: String,
    pub scale_type: String,
    pub complex_values: bool,
    pub scale: Vec<f64>,
    pub traces: Vec<RawTrace>,
}

impl Default for SpiceRawDataset {
    fn default() -> Self {
        Self {
            title: "cfd_solvers".into(),
            plot_name: String::new(),
            scale_name: String::new(),
            scale_type: String::new(),
            complex_values: false,
            scale: Vec::new(),
            traces: Vec::new(),
        }
    }
}

/// The ASCII RAW text (`%.17g` values; complex values as `(re,im)`).
pub fn spice_raw_text(d: &SpiceRawDataset) -> Result<String> {
    require(!d.scale.is_empty(), "SPICE RAW dataset requires points")?;
    require(
        !d.scale_name.is_empty() && !d.scale_type.is_empty(),
        "SPICE RAW dataset requires scale metadata",
    )?;
    require(
        d.traces
            .iter()
            .all(|t| !t.name.is_empty() && t.values.len() == d.scale.len()),
        "SPICE RAW trace size/name mismatch",
    )?;
    let mut o = format!(
        "Title: {}\nPlotname: {}\nFlags: {}\nNo. Variables: {}\nNo. Points: {}\nVariables:\n\t0\t{}\t{}\n",
        d.title,
        d.plot_name,
        if d.complex_values { "complex" } else { "real" },
        d.traces.len() + 1,
        d.scale.len(),
        d.scale_name,
        d.scale_type
    );
    for (i, t) in d.traces.iter().enumerate() {
        o.push_str(&format!("\t{}\t{}\t{}\n", i + 1, t.name, t.kind));
    }
    o.push_str("Values:\n");
    for (point, &s) in d.scale.iter().enumerate() {
        if d.complex_values {
            o.push_str(&format!("{point}\t({},0)\n", format_g17(s)));
        } else {
            o.push_str(&format!("{point}\t{}\n", format_g17(s)));
        }
        for t in &d.traces {
            let v = t.values[point];
            if d.complex_values {
                o.push_str(&format!("\t({},{})\n", format_g17(v.re), format_g17(v.im)));
            } else {
                o.push_str(&format!("\t{}\n", format_g17(v.re)));
            }
        }
    }
    Ok(o)
}

pub fn write_spice_raw(output: &mut dyn std::io::Write, dataset: &SpiceRawDataset) -> Result<()> {
    let text = spice_raw_text(dataset)?;
    output
        .write_all(text.as_bytes())
        .map_err(|e| SolverError::Io(format!("cannot write SPICE RAW: {e}")))
}

fn trim(s: &str) -> &str {
    s.trim_matches([' ', '\t', '\r', '\n'])
}

fn after_colon(line: &str) -> Result<&str> {
    line.find(':')
        .map(|c| trim(&line[c + 1..]))
        .ok_or_else(|| SolverError::InvalidInput("malformed SPICE RAW header".into()))
}

fn leading_unsigned(text: &str) -> Result<usize> {
    let digits: String = text
        .trim_start()
        .chars()
        .take_while(char::is_ascii_digit)
        .collect();
    digits
        .parse()
        .map_err(|_| SolverError::InvalidInput("malformed SPICE RAW count".into()))
}

fn parse_complex(text: &str) -> Result<Complex> {
    let t = trim(text);
    require(!t.is_empty(), "empty SPICE RAW value")?;
    if !t.starts_with('(') {
        return Ok(Complex::new(strtod_prefix(t)?, 0.0));
    }
    require(t.ends_with(')'), "malformed SPICE RAW complex value")?;
    let comma = t
        .find(',')
        .ok_or_else(|| SolverError::InvalidInput("malformed SPICE RAW complex pair".into()))?;
    Ok(Complex::new(
        strtod_prefix(&t[1..comma])?,
        strtod_prefix(&t[comma + 1..t.len() - 1])?,
    ))
}

/// Read the ASCII RAW format written by [`write_spice_raw`] (and compatible simulators).
pub fn read_spice_raw(text: &str) -> Result<SpiceRawDataset> {
    let mut d = SpiceRawDataset::default();
    let (mut variables, mut points) = (0usize, 0usize);
    let (mut in_variables, mut found_values) = (false, false);
    let lines = getline_split(text, b'\n');
    let mut it = lines.into_iter();
    for line in it.by_ref() {
        let clean = trim(line);
        if clean.is_empty() {
            continue;
        }
        if clean.starts_with("Title:") {
            d.title = after_colon(clean)?.to_string();
        } else if clean.starts_with("Plotname:") {
            d.plot_name = after_colon(clean)?.to_string();
        } else if clean.starts_with("Flags:") {
            d.complex_values = after_colon(clean)?.contains("complex");
        } else if clean.starts_with("No. Variables:") {
            variables = leading_unsigned(after_colon(clean)?)?;
        } else if clean.starts_with("No. Points:") {
            points = leading_unsigned(after_colon(clean)?)?;
        } else if clean == "Variables:" {
            in_variables = true;
        } else if clean == "Values:" {
            found_values = true;
            break;
        } else if in_variables {
            let row: Vec<&str> = clean.split_whitespace().collect();
            let (Some(index), Some(name), Some(kind)) = (row.first(), row.get(1), row.get(2))
            else {
                return Err(SolverError::InvalidInput(
                    "malformed SPICE RAW variable row".into(),
                ));
            };
            let index: usize = index.parse().map_err(|_| {
                SolverError::InvalidInput("malformed SPICE RAW variable row".into())
            })?;
            if index == 0 {
                d.scale_name = name.to_string();
                d.scale_type = kind.to_string();
            } else {
                d.traces.push(RawTrace {
                    name: name.to_string(),
                    kind: kind.to_string(),
                    values: Vec::new(),
                });
            }
        }
    }
    require(
        found_values && variables != 0 && variables == d.traces.len() + 1 && points != 0,
        "incomplete SPICE RAW header",
    )?;
    let mut next_line = |what: &str| -> Result<String> {
        loop {
            let line = it
                .next()
                .ok_or_else(|| SolverError::InvalidInput(what.to_string()))?;
            let t = trim(line);
            if !t.is_empty() {
                return Ok(t.to_string());
            }
        }
    };
    for point in 0..points {
        let line = next_line("truncated SPICE RAW values")?;
        let mut parts = line.split_whitespace();
        let index = parts.next().and_then(|s| s.parse::<usize>().ok());
        let value = parts.next();
        let (Some(index), Some(value)) = (index, value) else {
            return Err(SolverError::InvalidInput(
                "SPICE RAW point index mismatch".into(),
            ));
        };
        require(index == point, "SPICE RAW point index mismatch")?;
        d.scale.push(parse_complex(value)?.re);
        for ti in 0..d.traces.len() {
            let line = next_line("truncated SPICE RAW trace values")?;
            let v = parse_complex(&line)?;
            d.traces[ti].values.push(v);
        }
    }
    Ok(d)
}

fn node_traces(node_names: &[String]) -> Vec<RawTrace> {
    node_names
        .iter()
        .skip(1)
        .map(|n| RawTrace::voltage(format!("v({n})")))
        .collect()
}

pub fn raw_from_transient(
    points: &[TransientPoint],
    node_names: &[String],
    title: &str,
) -> Result<SpiceRawDataset> {
    let first = points.first().ok_or_else(|| {
        SolverError::InvalidInput("cannot export empty transient RAW dataset".into())
    })?;
    require(
        first.node_voltage.len() == node_names.len(),
        "transient RAW node-name mismatch",
    )?;
    let mut out = SpiceRawDataset {
        title: title.into(),
        plot_name: "Transient Analysis".into(),
        scale_name: "time".into(),
        scale_type: "time".into(),
        traces: node_traces(node_names),
        ..Default::default()
    };
    for p in points {
        require(
            p.node_voltage.len() == node_names.len(),
            "transient RAW topology changed",
        )?;
        out.scale.push(p.time);
        for node in 1..node_names.len() {
            out.traces[node - 1]
                .values
                .push(Complex::new(p.node_voltage[node], 0.0));
        }
    }
    Ok(out)
}

pub fn raw_from_ac(
    points: &[AcSolution],
    node_names: &[String],
    title: &str,
) -> Result<SpiceRawDataset> {
    let first = points
        .first()
        .ok_or_else(|| SolverError::InvalidInput("cannot export empty AC RAW dataset".into()))?;
    require(
        first.node_voltage.len() == node_names.len(),
        "AC RAW node-name mismatch",
    )?;
    let mut out = SpiceRawDataset {
        title: title.into(),
        plot_name: "AC Analysis".into(),
        scale_name: "frequency".into(),
        scale_type: "frequency".into(),
        complex_values: true,
        traces: node_traces(node_names),
        ..Default::default()
    };
    for p in points {
        require(
            p.node_voltage.len() == node_names.len(),
            "AC RAW topology changed",
        )?;
        out.scale.push(p.frequency_hz);
        for node in 1..node_names.len() {
            out.traces[node - 1].values.push(p.node_voltage[node]);
        }
    }
    Ok(out)
}
