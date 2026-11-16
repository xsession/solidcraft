//! Unit-tagged case configuration, material database, result provenance and the validation-case
//! catalogue, ported from upstream `cfd/workflow/{schema,provenance,validation_catalog}`
//! (MIT, cfd_solvers contributors).
use super::text::{c_isspace, format_g17, strtod_prefix};
use crate::{require, Result, SolverError};
use std::{collections::BTreeMap, path::Path};

/// A finite value with a free-form unit string.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct UnitValue {
    pub value: f64,
    pub unit: String,
}

impl UnitValue {
    pub fn new(value: f64, unit: impl Into<String>) -> Self {
        Self {
            value,
            unit: unit.into(),
        }
    }
}

/// Validation rule for one case parameter. An empty `unit` accepts any unit.
#[derive(Clone, Debug, PartialEq)]
pub struct ParameterRule {
    pub key: String,
    pub unit: String,
    pub minimum: Option<f64>,
    pub maximum: Option<f64>,
    pub required: bool,
}

impl ParameterRule {
    /// Required rule without bounds (upstream default `required=true`).
    pub fn new(key: impl Into<String>, unit: impl Into<String>) -> Self {
        Self {
            key: key.into(),
            unit: unit.into(),
            minimum: None,
            maximum: None,
            required: true,
        }
    }
}

pub(crate) fn json_escape_basic(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for c in s.chars() {
        if c == '"' || c == '\\' {
            o.push('\\');
        }
        o.push(c);
    }
    o
}

/// Key → unit value map. Keys are kept sorted, which is the order upstream's `to_json` emits.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CaseConfig {
    values: BTreeMap<String, UnitValue>,
}

impl CaseConfig {
    pub fn set(&mut self, key: impl Into<String>, value: UnitValue) {
        self.values.insert(key.into(), value);
    }
    pub fn at(&self, key: &str) -> Result<&UnitValue> {
        self.values
            .get(key)
            .ok_or_else(|| SolverError::InvalidInput("case parameter not found".into()))
    }
    pub fn len(&self) -> usize {
        self.values.len()
    }
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
    pub fn validate(&self, rules: &[ParameterRule]) -> Result<()> {
        for r in rules {
            let Some(v) = self.values.get(&r.key) else {
                require(
                    !r.required,
                    format!("required case parameter missing: {}", r.key),
                )?;
                continue;
            };
            require(v.value.is_finite(), "non-finite case parameter")?;
            require(
                r.unit.is_empty() || v.unit == r.unit,
                format!("case parameter unit mismatch: {}", r.key),
            )?;
            if let Some(min) = r.minimum {
                require(v.value >= min, "case parameter below minimum")?;
            }
            if let Some(max) = r.maximum {
                require(v.value <= max, "case parameter above maximum")?;
            }
        }
        Ok(())
    }
    /// `{"parameters":{"key":{"value":v,"unit":"u"},...}}`, keys sorted, values `%.17g`.
    pub fn to_json(&self) -> String {
        let mut o = String::from("{\"parameters\":{");
        for (i, (k, v)) in self.values.iter().enumerate() {
            if i > 0 {
                o.push(',');
            }
            o.push_str(&format!(
                "\"{}\":{{\"value\":{},\"unit\":\"{}\"}}",
                json_escape_basic(k),
                format_g17(v.value),
                json_escape_basic(&v.unit)
            ));
        }
        o.push_str("}}");
        o
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Material {
    pub name: String,
    pub property: BTreeMap<String, UnitValue>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MaterialDatabase {
    materials: BTreeMap<String, Material>,
}

impl MaterialDatabase {
    pub fn add(&mut self, material: Material) -> Result<()> {
        require(!material.name.is_empty(), "material name is empty")?;
        require(
            material.property.values().all(|v| v.value.is_finite()),
            "non-finite material property",
        )?;
        self.materials.insert(material.name.clone(), material);
        Ok(())
    }
    pub fn at(&self, name: &str) -> Result<&Material> {
        self.materials
            .get(name)
            .ok_or_else(|| SolverError::InvalidInput("material not found".into()))
    }
    pub fn len(&self) -> usize {
        self.materials.len()
    }
    pub fn is_empty(&self) -> bool {
        self.materials.is_empty()
    }
    pub fn to_json(&self) -> String {
        let mut o = String::from("{\"materials\":[");
        for (i, (name, m)) in self.materials.iter().enumerate() {
            if i > 0 {
                o.push(',');
            }
            o.push_str(&format!(
                "{{\"name\":\"{}\",\"properties\":{{",
                json_escape_basic(name)
            ));
            for (j, (k, v)) in m.property.iter().enumerate() {
                if j > 0 {
                    o.push(',');
                }
                o.push_str(&format!(
                    "\"{}\":{{\"value\":{},\"unit\":\"{}\"}}",
                    json_escape_basic(k),
                    format_g17(v.value),
                    json_escape_basic(&v.unit)
                ));
            }
            o.push_str("}}");
        }
        o.push_str("]}");
        o
    }
}

fn skip_ws(b: &[u8], mut i: usize) -> usize {
    while i < b.len() && c_isspace(b[i]) {
        i += 1;
    }
    i
}

fn is_num_char(c: u8) -> bool {
    c.is_ascii_digit() || matches!(c, b'-' | b'+' | b'.' | b'e' | b'E')
}

/// Try the upstream item pattern
/// `"key"\s*:\s*\{\s*"value"\s*:\s*NUMBER\s*,\s*"unit"\s*:\s*"unit"\s*\}` at byte `i`.
/// Returns `(key, number text, unit, end)`.
fn match_json_item(s: &str, i: usize) -> Option<(&str, &str, &str, usize)> {
    let b = s.as_bytes();
    let expect = |j: usize, c: u8| (j < b.len() && b[j] == c).then_some(j + 1);
    let mut j = expect(i, b'"')?;
    let key_start = j;
    while j < b.len() && b[j] != b'"' {
        j += 1;
    }
    if j == key_start || j >= b.len() {
        return None;
    }
    let key = s.get(key_start..j)?;
    j = expect(skip_ws(b, j + 1), b':')?;
    j = expect(skip_ws(b, j), b'{')?;
    j = skip_ws(b, j);
    if !b.get(j..)?.starts_with(b"\"value\"") {
        return None;
    }
    j = expect(skip_ws(b, j + 7), b':')?;
    j = skip_ws(b, j);
    let num_start = j;
    while j < b.len() && is_num_char(b[j]) {
        j += 1;
    }
    if j == num_start {
        return None;
    }
    let num = s.get(num_start..j)?;
    j = expect(skip_ws(b, j), b',')?;
    j = skip_ws(b, j);
    if !b.get(j..)?.starts_with(b"\"unit\"") {
        return None;
    }
    j = expect(skip_ws(b, j + 6), b':')?;
    j = expect(skip_ws(b, j), b'"')?;
    let unit_start = j;
    while j < b.len() && b[j] != b'"' {
        j += 1;
    }
    if j >= b.len() {
        return None;
    }
    let unit = s.get(unit_start..j)?;
    j = expect(skip_ws(b, j + 1), b'}')?;
    Some((key, num, unit, j))
}

/// Scan `json` for every `"key":{"value":number,"unit":"..."}` item (upstream regex semantics:
/// leftmost, non-overlapping). Numbers are read like `std::stod` (longest valid prefix); a token
/// with no numeric prefix is an error.
pub fn case_config_from_json(json: &str) -> Result<CaseConfig> {
    let mut c = CaseConfig::default();
    let b = json.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'"' {
            if let Some((key, num, unit, end)) = match_json_item(json, i) {
                c.set(key, UnitValue::new(strtod_prefix(num)?, unit));
                i = end;
                continue;
            }
        }
        i += 1;
    }
    Ok(c)
}

/// Lines of the form `key: number unit` (optionally indented, trailing whitespace allowed, `#`
/// not allowed in the unit); every other line is ignored.
pub fn case_config_from_yaml(yaml: &str) -> Result<CaseConfig> {
    let mut c = CaseConfig::default();
    for line in super::text::getline_split(yaml, b'\n') {
        let b = line.as_bytes();
        let mut j = skip_ws(b, 0);
        let ks = j;
        while j < b.len() && (b[j].is_ascii_alphanumeric() || matches!(b[j], b'_' | b'.' | b'-')) {
            j += 1;
        }
        if j == ks {
            continue;
        }
        let key = &line[ks..j];
        j = skip_ws(b, j);
        if j >= b.len() || b[j] != b':' {
            continue;
        }
        j = skip_ws(b, j + 1);
        let ns = j;
        while j < b.len() && is_num_char(b[j]) {
            j += 1;
        }
        if j == ns {
            continue;
        }
        let num = &line[ns..j];
        let ws = j;
        j = skip_ws(b, j);
        if j == ws {
            continue;
        }
        let us = j;
        while j < b.len() && b[j] != b'#' && !c_isspace(b[j]) {
            j += 1;
        }
        if j == us || skip_ws(b, j) != b.len() {
            continue;
        }
        c.set(key, UnitValue::new(strtod_prefix(num)?, &line[us..j]));
    }
    Ok(c)
}

/// Solver identity, version, commit, UTC timestamp and free metadata (sorted by key).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ResultProvenance {
    pub solver: String,
    pub version: String,
    pub git_commit: String,
    pub timestamp_utc: String,
    pub metadata: BTreeMap<String, String>,
}

/// `YYYY-MM-DDTHH:MM:SSZ` for seconds since the Unix epoch (proleptic Gregorian, UTC).
pub fn format_utc_timestamp(unix_seconds: i64) -> String {
    let days = unix_seconds.div_euclid(86_400);
    let secs = unix_seconds.rem_euclid(86_400);
    // Howard Hinnant's civil_from_days.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        secs / 3600,
        (secs / 60) % 60,
        secs % 60
    )
}

/// Provenance stamped with the current UTC time (upstream default commit: `"unknown"`).
pub fn make_provenance(
    solver: impl Into<String>,
    version: impl Into<String>,
    git_commit: impl Into<String>,
) -> ResultProvenance {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| i64::try_from(d.as_secs()).unwrap_or(i64::MAX))
        .unwrap_or(0);
    ResultProvenance {
        solver: solver.into(),
        version: version.into(),
        git_commit: git_commit.into(),
        timestamp_utc: format_utc_timestamp(now),
        metadata: BTreeMap::new(),
    }
}

pub fn provenance_json(p: &ResultProvenance) -> String {
    let mut o = format!(
        "{{\"solver\":\"{}\",\"version\":\"{}\",\"git_commit\":\"{}\",\"timestamp_utc\":\"{}\",\"metadata\":{{",
        json_escape_basic(&p.solver),
        json_escape_basic(&p.version),
        json_escape_basic(&p.git_commit),
        p.timestamp_utc
    );
    for (i, (k, v)) in p.metadata.iter().enumerate() {
        if i > 0 {
            o.push(',');
        }
        o.push_str(&format!(
            "\"{}\":\"{}\"",
            json_escape_basic(k),
            json_escape_basic(v)
        ));
    }
    o.push_str("}}");
    o
}

pub fn write_provenance(path: &Path, p: &ResultProvenance) -> Result<()> {
    std::fs::write(path, provenance_json(p) + "\n")
        .map_err(|e| SolverError::Io(format!("cannot write provenance: {e}")))
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ValidationCase {
    pub name: &'static str,
    pub family: &'static str,
    pub command: &'static str,
    pub criterion: &'static str,
}

const fn vc(
    name: &'static str,
    family: &'static str,
    command: &'static str,
    criterion: &'static str,
) -> ValidationCase {
    ValidationCase {
        name,
        family,
        command,
        criterion,
    }
}

static VALIDATION_CATALOG: [ValidationCase; 6] = [
    vc(
        "taylor-green-fvm",
        "FVM",
        "cfd-solve fvm-taylor-green2d",
        "velocity refinement + divergence",
    ),
    vc(
        "poiseuille-lbm",
        "LBM",
        "cfd-solve lbm-poiseuille",
        "analytic parabolic profile",
    ),
    vc(
        "poisson-fem",
        "FEM",
        "cfd-solve fem-poisson3d",
        "manufactured solution",
    ),
    vc(
        "cpml-fdtd",
        "FDTD",
        "cfd-solve fdtd-cpml-tfsf",
        "absorbed residual energy",
    ),
    vc(
        "lens-optics",
        "Optics",
        "cfd-solve optics-lens",
        "paraxial/real-ray finite result",
    ),
    vc(
        "pnp-electrochem",
        "Electrochemistry",
        "cfd-solve electrochem-pnp-poly",
        "analytic charged-slab potential",
    ),
];

/// The six upstream reference validation cases.
pub fn validation_catalog() -> &'static [ValidationCase] {
    &VALIDATION_CATALOG
}
