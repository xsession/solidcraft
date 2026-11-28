//! Where a study of the model stands and what it breathes: the air around it, the dust in that
//! air, and the year of conditions it works through.
//!
//! - **Surroundings** (external studies): a breeze from one side or still air, open above and on
//!   the downwind side, a floor under the model. They set the six domain faces.
//! - **Dust**: sizes and their shares of the mass, how much is in the air, and where it is
//!   counted from: the air coming in at the domain's sides, or the air crossing into a box,
//!   such as the inside of an enclosure, so what is counted is what came in through its
//!   openings. The companion follows it through the settled flow (`dust_deposition_v1`).
//! - **Climate**: conditions over a year (a season's night and afternoon, say), each with the
//!   air's temperature and dew point, the breeze, the dust, and the hours of the year it stands
//!   for. One study per condition is run; [`summarize`] adds their results up over the year:
//!   the dust each part gathers, and the hours each part spends humid or wet.
//!
//! Portable: no geometry is read here; the engine resolves bodies and runs.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::{Error, Result};

/// Most dust sizes, and most climate conditions (one study each).
pub const MAX_DUST_SIZES: usize = 16;
pub const MAX_CONDITIONS: usize = 16;
/// Hours in a year.
pub const YEAR_HOURS: f64 = 8760.0;

fn bad(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

fn yes() -> bool {
    true
}

/// The side a breeze blows from.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub enum BreezeFrom {
    #[default]
    #[serde(rename = "-x")]
    MinusX,
    #[serde(rename = "+x")]
    PlusX,
    #[serde(rename = "-y")]
    MinusY,
    #[serde(rename = "+y")]
    PlusY,
}

impl BreezeFrom {
    /// The domain faces the breeze comes in at and goes out at, and the two it runs along.
    fn patches(self) -> (usize, usize, [usize; 2]) {
        match self {
            Self::MinusX => (0, 1, [2, 3]),
            Self::PlusX => (1, 0, [2, 3]),
            Self::MinusY => (2, 3, [0, 1]),
            Self::PlusY => (3, 2, [0, 1]),
        }
    }
}

/// The air around an external study, with Z up: a breeze (or still air), open above and
/// downwind, the sides along the breeze free to slip, and a floor under the model.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelSurroundings {
    /// 0 for still air: then the breeze's sides are open too.
    pub breeze_m_per_s: f64,
    #[serde(default)]
    pub breeze_from: BreezeFrom,
    /// A floor (a wall) under the model; otherwise the bottom is open as well.
    #[serde(default = "yes")]
    pub floor: bool,
}

impl ModelSurroundings {
    pub fn validate(&self) -> Result<()> {
        if !(self.breeze_m_per_s.is_finite() && (0.0..=100.0).contains(&self.breeze_m_per_s)) {
            return Err(bad("the breeze is 0 to 100 m/s"));
        }
        Ok(())
    }

    /// The air's starting velocity: the breeze, along its way.
    pub fn initial_velocity(&self) -> [f64; 3] {
        let b = self.breeze_m_per_s;
        match self.breeze_from {
            BreezeFrom::MinusX => [b, 0.0, 0.0],
            BreezeFrom::PlusX => [-b, 0.0, 0.0],
            BreezeFrom::MinusY => [0.0, b, 0.0],
            BreezeFrom::PlusY => [0.0, -b, 0.0],
        }
    }

    /// The role of each domain face (`case_scene::PATCHES` order: left, right, bottom, top,
    /// front, back; front and back are the Z faces, so front is the floor).
    pub fn roles(&self) -> [crate::case_scene::RoleInput; 6] {
        use crate::case_scene::RoleInput::{Inlet, Outlet, Symmetry, Wall};
        let open = Outlet { pressure: 0.0 };
        let mut roles = [open, open, open, open, open, open];
        let (upwind, _, along) = self.breeze_from.patches();
        if self.breeze_m_per_s > 0.0 {
            roles[upwind] = Inlet { speed: self.breeze_m_per_s };
            for side in along {
                roles[side] = Symmetry;
            }
        }
        if self.floor {
            roles[4] = Wall;
        }
        roles
    }
}

/// One size of dust and its share of the mass.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DustSize {
    pub diameter_um: f64,
    pub share: f64,
}

/// Where the dust is counted from.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DustEntry {
    /// The air coming in at the domain's sides.
    #[default]
    Air,
    /// The air crossing into a box of the design (mm): with the box's sides inside an
    /// enclosure's walls, the dust that comes in through its openings.
    Box { low_mm: [f64; 3], high_mm: [f64; 3] },
}

/// Dust in the air of a study.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelDust {
    /// All sizes together, mg/m³ of air.
    pub concentration_mg_per_m3: f64,
    #[serde(default = "default_sizes")]
    pub sizes: Vec<DustSize>,
    /// Of the particles themselves; mineral dust is about 2 600, organic dust about 1 000.
    #[serde(default = "default_density")]
    pub density_kg_per_m3: f64,
    /// Followed per size; more give smoother maps and take longer.
    #[serde(default = "default_particles")]
    pub particles_per_size: usize,
    #[serde(default)]
    pub entry: DustEntry,
    /// A particle still in the air after this long counts as airborne.
    #[serde(default = "default_max_time")]
    pub max_time_s: f64,
}

/// Coarse outdoor dust by mass: little below 2.5 µm, most of it from 5 to 20 µm.
pub fn default_sizes() -> Vec<DustSize> {
    [(1.0, 0.05), (2.5, 0.10), (5.0, 0.20), (10.0, 0.35), (20.0, 0.30)]
        .into_iter()
        .map(|(diameter_um, share)| DustSize { diameter_um, share })
        .collect()
}
fn default_density() -> f64 {
    2000.0
}
fn default_particles() -> usize {
    1000
}
fn default_max_time() -> f64 {
    3600.0
}

impl ModelDust {
    pub fn validate(&self) -> Result<()> {
        if !(self.concentration_mg_per_m3.is_finite() && (0.0..=1e6).contains(&self.concentration_mg_per_m3)) {
            return Err(bad("the dust concentration is 0 to 1e6 mg/m³"));
        }
        if self.sizes.is_empty() || self.sizes.len() > MAX_DUST_SIZES {
            return Err(bad(format!("give 1 to {MAX_DUST_SIZES} dust sizes")));
        }
        for size in &self.sizes {
            if !(size.diameter_um.is_finite() && (0.01..=1000.0).contains(&size.diameter_um)) {
                return Err(bad("dust diameters are 0.01 to 1000 µm"));
            }
            if !(size.share.is_finite() && size.share >= 0.0) {
                return Err(bad("dust shares are 0 or more"));
            }
        }
        if self.sizes.iter().map(|s| s.share).sum::<f64>() <= 0.0 {
            return Err(bad("the dust sizes' shares add up to nothing"));
        }
        if !(self.density_kg_per_m3.is_finite() && (1.0..=1e5).contains(&self.density_kg_per_m3)) {
            return Err(bad("the dust density is 1 to 1e5 kg/m³"));
        }
        if self.particles_per_size == 0 || self.particles_per_size.saturating_mul(self.sizes.len()) > 100_000 {
            return Err(bad("follow at least 1 particle per size and at most 100 000 in all"));
        }
        if !(self.max_time_s.is_finite() && self.max_time_s > 0.0 && self.max_time_s <= 1e7) {
            return Err(bad("the dust's time in the air is above 0 and at most 1e7 s"));
        }
        if let DustEntry::Box { low_mm, high_mm } = &self.entry
            && !(0..3).all(|a| low_mm[a].is_finite() && high_mm[a].is_finite() && low_mm[a] < high_mm[a])
        {
            return Err(bad("the dust box needs finite corners, the low one below the high one"));
        }
        Ok(())
    }

    /// The companion's `dust` block: `origin_mm` is the design point at the domain's origin,
    /// `gravity` the study's, `air_k` the air's temperature.
    pub fn case_block(&self, origin_mm: [f64; 3], gravity: [f64; 3], air_k: f64) -> Value {
        let total: f64 = self.sizes.iter().map(|s| s.share).sum();
        let kg_per_m3 = self.concentration_mg_per_m3 * 1e-6;
        let classes: Vec<Value> =
            self.sizes.iter().map(|s| json!({"diameter_m": s.diameter_um * 1e-6, "concentration_kg_per_m3": kg_per_m3 * s.share / total})).collect();
        let mut block = json!({
            "classes": classes,
            "gravity_m_per_s2": gravity,
            "density_kg_per_m3": self.density_kg_per_m3,
            "particles_per_class": self.particles_per_size,
            "max_time_s": self.max_time_s,
            "temperature_k": air_k,
        });
        if let DustEntry::Box { low_mm, high_mm } = &self.entry {
            block["release"] = json!("box");
            let low = [0, 1, 2].map(|a| (low_mm[a] - origin_mm[a]) / 1000.0);
            let high = [0, 1, 2].map(|a| (high_mm[a] - origin_mm[a]) / 1000.0);
            block["release_box"] = json!({"min_m": low, "max_m": high});
        }
        block
    }
}

/// One condition of the year: the air's temperature and dew point, and optionally a breeze and
/// dust different from the study's, for `hours_per_year` hours.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ClimateCondition {
    /// Letters, digits, spaces, `-` and `_`, as it names its study.
    pub name: String,
    pub hours_per_year: f64,
    pub ambient_c: f64,
    pub dew_point_c: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub breeze_m_per_s: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dust_mg_per_m3: Option<f64>,
}

impl ClimateCondition {
    pub fn validate(&self) -> Result<()> {
        let name = self.name.trim();
        if name.is_empty() || name.len() > 60 || !name.chars().all(|c| c.is_alphanumeric() || matches!(c, ' ' | '-' | '_')) {
            return Err(bad("a condition's name is 1 to 60 letters, digits, spaces, - and _"));
        }
        if !(self.hours_per_year.is_finite() && self.hours_per_year > 0.0 && self.hours_per_year <= YEAR_HOURS) {
            return Err(bad(format!("{}: its hours are above 0 and at most {YEAR_HOURS}", self.name)));
        }
        if !(self.ambient_c.is_finite() && (-90.0..=150.0).contains(&self.ambient_c)) {
            return Err(bad(format!("{}: the air is -90 to 150 °C", self.name)));
        }
        if !(self.dew_point_c.is_finite() && self.dew_point_c <= self.ambient_c + 1e-9 && self.dew_point_c >= -90.0) {
            return Err(bad(format!("{}: the dew point is -90 °C to the air's temperature", self.name)));
        }
        if self.breeze_m_per_s.is_some_and(|b| !(b.is_finite() && (0.0..=100.0).contains(&b))) {
            return Err(bad(format!("{}: the breeze is 0 to 100 m/s", self.name)));
        }
        if self.dust_mg_per_m3.is_some_and(|d| !(d.is_finite() && (0.0..=1e6).contains(&d))) {
            return Err(bad(format!("{}: the dust is 0 to 1e6 mg/m³", self.name)));
        }
        Ok(())
    }

    /// The air's relative humidity, %.
    pub fn relative_humidity_pct(&self) -> f64 {
        relative_humidity_pct(self.dew_point_c, self.ambient_c)
    }
}

/// The year of a study: its conditions, and the relative humidity at a surface above which dust
/// on it counts as humid (salts in dust take up water from about 60 % up).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ModelClimate {
    pub conditions: Vec<ClimateCondition>,
    #[serde(default = "default_humid")]
    pub humid_limit_pct: f64,
}

fn default_humid() -> f64 {
    60.0
}

impl ModelClimate {
    pub fn validate(&self) -> Result<()> {
        if self.conditions.is_empty() || self.conditions.len() > MAX_CONDITIONS {
            return Err(bad(format!("a climate has 1 to {MAX_CONDITIONS} conditions")));
        }
        let mut names = std::collections::BTreeSet::new();
        for c in &self.conditions {
            c.validate()?;
            if !names.insert(c.name.trim().to_lowercase()) {
                return Err(bad(format!("two conditions are called {}", c.name)));
            }
        }
        let hours: f64 = self.conditions.iter().map(|c| c.hours_per_year).sum();
        if hours > YEAR_HOURS * (1.0 + 1e-9) {
            return Err(bad(format!("the conditions add up to {hours:.0} hours; a year has {YEAR_HOURS}")));
        }
        if !(self.humid_limit_pct.is_finite() && (1.0..=100.0).contains(&self.humid_limit_pct)) {
            return Err(bad("the humid limit is 1 to 100 %"));
        }
        Ok(())
    }
}

/// Water's saturation vapour pressure over liquid water, hPa: the Magnus form with Alduchov and
/// Eskridge's coefficients (1996), within 0.4 % from -40 to 50 °C.
pub fn saturation_hpa(t_c: f64) -> f64 {
    6.1094 * (17.625 * t_c / (t_c + 243.04)).exp()
}

/// Relative humidity at a surface at `surface_c` in air whose dew point is `dew_point_c`, %:
/// the air's vapour pressure over the saturation pressure at the surface. Over 100 the surface
/// is below the dew point and water condenses on it.
pub fn relative_humidity_pct(dew_point_c: f64, surface_c: f64) -> f64 {
    100.0 * saturation_hpa(dew_point_c) / saturation_hpa(surface_c)
}

/// One condition's run, as the engine finds it.
pub struct ConditionRun<'a> {
    pub condition: &'a ClimateCondition,
    pub study: String,
    /// The run's id and its files (`dust_deposition.json`, `dust_deposits.csv`,
    /// `thermal_report.json`); none before it has run.
    pub run: Option<(String, &'a [crate::Artifact])>,
    /// The design point at the domain's origin, mm, for the deposits' places.
    pub origin_mm: [f64; 3],
}

/// What a condition contributed.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct ConditionSummary {
    pub name: String,
    pub study: String,
    pub run: Option<String>,
    pub hours_per_year: f64,
    pub ambient_c: f64,
    pub dew_point_c: f64,
    pub air_humidity_pct: f64,
    /// Dust carried in (entry `air`: into the domain; `box`: into the box), g a year.
    pub dust_in_g_per_year: f64,
    /// Of it, deposited anywhere, g a year.
    pub deposited_g_per_year: f64,
}

/// A part over the year.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct PartSummary {
    pub body: String,
    /// Dust that settles on it, g a year.
    pub dust_g_per_year: f64,
    /// Its share of all the dust deposited.
    pub dust_share: f64,
    /// Hours a year its coldest point is at or above the humid limit, and below the dew point.
    pub humid_hours: f64,
    pub condensing_hours: f64,
    /// Hours a year with dust on it at or above the humid limit, weighed by how much of its dust
    /// is: the hours its deposits spend humid.
    pub humid_dust_hours: f64,
    pub coldest_c: Option<f64>,
    pub warmest_c: Option<f64>,
}

/// The year added up.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct ClimateSummary {
    pub hours_covered: f64,
    pub humid_limit_pct: f64,
    pub conditions: Vec<ConditionSummary>,
    /// Most dust first; the domain's own faces (`left`, `front`: the floor, …) included.
    pub parts: Vec<PartSummary>,
    pub dust_in_g_per_year: f64,
    pub deposited_g_per_year: f64,
    /// Conditions with no run yet: the year above leaves them out.
    pub missing: Vec<String>,
    /// Each deposit over the year: x, y, z (design mm), the part, its diameter (µm), g a year
    /// and the hours a year it is humid, as CSV; written to a file, never sent in a reply.
    #[serde(skip)]
    pub deposits_csv: String,
}

fn artifact<'a>(files: &'a [crate::Artifact], name: &str) -> Option<&'a str> {
    files.iter().find(|a| a.name == name).map(|a| a.text.as_str())
}

/// Add the conditions' runs up over a year.
pub fn summarize(runs: &[ConditionRun<'_>], humid_limit_pct: f64) -> Result<ClimateSummary> {
    use std::collections::BTreeMap;
    #[derive(Default)]
    struct Part {
        dust: f64,
        humid: f64,
        condensing: f64,
        humid_dust: f64,
        coldest: Option<f64>,
        warmest: Option<f64>,
    }
    let mut parts: BTreeMap<String, Part> = BTreeMap::new();
    let mut conditions = Vec::new();
    let mut missing = Vec::new();
    let mut deposits = String::from("x_mm,y_mm,z_mm,part,diameter_um,g_per_year,humid_hours,condition\n");
    let (mut dust_in, mut deposited, mut covered) = (0.0, 0.0, 0.0);
    for r in runs {
        let c = r.condition;
        let mut summary = ConditionSummary {
            name: c.name.clone(),
            study: r.study.clone(),
            run: None,
            hours_per_year: c.hours_per_year,
            ambient_c: c.ambient_c,
            dew_point_c: c.dew_point_c,
            air_humidity_pct: c.relative_humidity_pct(),
            dust_in_g_per_year: 0.0,
            deposited_g_per_year: 0.0,
        };
        let Some((run_id, files)) = &r.run else {
            missing.push(c.name.clone());
            conditions.push(summary);
            continue;
        };
        summary.run = Some(run_id.clone());
        covered += c.hours_per_year;
        // kg/s over the condition's hours, in grams.
        let grams = |kg_per_s: f64| kg_per_s * c.hours_per_year * 3600.0 * 1000.0;
        if let Some(text) = artifact(files, "thermal_report.json") {
            let report: Value = serde_json::from_str(text).map_err(|e| Error::Json(e.to_string()))?;
            for comp in report["components"].as_array().into_iter().flatten() {
                let (Some(name), Some(min_k), Some(max_k)) =
                    (comp["name"].as_str(), comp["min_temperature_k"].as_f64(), comp["max_temperature_k"].as_f64())
                else {
                    continue;
                };
                let (cold, warm) = (min_k - 273.15, max_k - 273.15);
                let part = parts.entry(name.to_string()).or_default();
                part.coldest = Some(part.coldest.map_or(cold, |v| v.min(cold)));
                part.warmest = Some(part.warmest.map_or(warm, |v| v.max(warm)));
                let rh = relative_humidity_pct(c.dew_point_c, cold);
                if rh >= humid_limit_pct {
                    part.humid += c.hours_per_year;
                }
                if cold <= c.dew_point_c {
                    part.condensing += c.hours_per_year;
                }
            }
        }
        if let Some(text) = artifact(files, "dust_deposition.json") {
            let report: Value = serde_json::from_str(text).map_err(|e| Error::Json(e.to_string()))?;
            if report["amount_unit"].as_str() != Some("kg/s") {
                return Err(bad(format!("{}: its dust was not counted as a flow (kg/s)", c.name)));
            }
            summary.dust_in_g_per_year = grams(report["entering"].as_f64().unwrap_or(0.0));
            summary.deposited_g_per_year = grams(report["deposited"].as_f64().unwrap_or(0.0));
        }
        if let Some(text) = artifact(files, "dust_deposits.csv") {
            // class,diameter_m,x_m,y_m,z_m,target,amount,surface_temperature_k,time_s
            for line in text.lines().skip(1) {
                let f: Vec<&str> = line.split(',').collect();
                let (Some(d), Some(x), Some(y), Some(z), Some(target), Some(amount), Some(t)) =
                    (f.get(1), f.get(2), f.get(3), f.get(4), f.get(5), f.get(6), f.get(7))
                else {
                    continue;
                };
                let num = |s: &str| s.trim().parse::<f64>().ok().filter(|v| v.is_finite());
                let (Some(d), Some(x), Some(y), Some(z), Some(amount), Some(t)) = (num(d), num(x), num(y), num(z), num(amount), num(t)) else {
                    continue;
                };
                let g = grams(amount);
                let humid = relative_humidity_pct(c.dew_point_c, t - 273.15) >= humid_limit_pct;
                let part = parts.entry((*target).to_string()).or_default();
                part.dust += g;
                if humid {
                    part.humid_dust += g * c.hours_per_year;
                }
                let at = [x, y, z].iter().zip(r.origin_mm).map(|(v, o)| v * 1000.0 + o).collect::<Vec<_>>();
                deposits.push_str(&format!(
                    "{:.3},{:.3},{:.3},{target},{:.3},{:.6e},{},{}\n",
                    at[0],
                    at[1],
                    at[2],
                    d * 1e6,
                    g,
                    if humid { c.hours_per_year } else { 0.0 },
                    c.name
                ));
            }
        }
        dust_in += summary.dust_in_g_per_year;
        deposited += summary.deposited_g_per_year;
        conditions.push(summary);
    }
    let total: f64 = parts.values().map(|p| p.dust).sum();
    let mut parts: Vec<PartSummary> = parts
        .into_iter()
        .map(|(body, p)| PartSummary {
            dust_share: if total > 0.0 { p.dust / total } else { 0.0 },
            // Hours its deposits spend humid: weighed by mass.
            humid_dust_hours: if p.dust > 0.0 { p.humid_dust / p.dust } else { 0.0 },
            dust_g_per_year: p.dust,
            humid_hours: p.humid,
            condensing_hours: p.condensing,
            coldest_c: p.coldest,
            warmest_c: p.warmest,
            body,
        })
        .collect();
    parts.sort_by(|a, b| b.dust_g_per_year.total_cmp(&a.dust_g_per_year).then_with(|| a.body.cmp(&b.body)));
    Ok(ClimateSummary {
        hours_covered: covered,
        humid_limit_pct,
        conditions,
        parts,
        dust_in_g_per_year: dust_in,
        deposited_g_per_year: deposited,
        missing,
        deposits_csv: deposits,
    })
}

/// The summary's parts as CSV.
pub fn parts_csv(summary: &ClimateSummary) -> String {
    let mut out = String::from("part,dust_g_per_year,dust_share,humid_hours,condensing_hours,humid_dust_hours,coldest_c,warmest_c\n");
    let opt = |v: Option<f64>| v.map_or(String::new(), |v| format!("{v:.2}"));
    for p in &summary.parts {
        out.push_str(&format!(
            "{},{:.6e},{:.4},{:.0},{:.0},{:.0},{},{}\n",
            p.body,
            p.dust_g_per_year,
            p.dust_share,
            p.humid_hours,
            p.condensing_hours,
            p.humid_dust_hours,
            opt(p.coldest_c),
            opt(p.warmest_c)
        ));
    }
    out
}

/// A representative year in four seasons, each a night (6 h a day) and three more quarters of
/// the day, from a season's lowest and highest temperature, its dew point and its dust; a
/// starting point to edit, not a climate record.
pub fn four_seasons(seasons: &[(&str, f64, f64, f64, f64)]) -> Vec<ClimateCondition> {
    let mut out = Vec::new();
    let days = 365.0 / seasons.len().max(1) as f64;
    for &(season, low, high, dew, dust) in seasons {
        let mean = 0.5 * (low + high);
        let swing = 0.5 * (high - low);
        // The quarter days centred at 03:00, 09:00, 15:00 and 21:00, the warmest at 15:00.
        for (hour, label, breeze, activity) in
            [(3.0, "night", 0.05, 0.5), (9.0, "morning", 0.15, 1.2), (15.0, "afternoon", 0.25, 1.5), (21.0, "evening", 0.1, 0.8)]
        {
            let ambient = mean + swing * (std::f64::consts::PI * (hour - 15.0) / 12.0).cos();
            out.push(ClimateCondition {
                name: format!("{season} {label}"),
                hours_per_year: days * 6.0,
                ambient_c: (ambient * 10.0).round() / 10.0,
                dew_point_c: dew.min(ambient),
                breeze_m_per_s: Some(breeze),
                dust_mg_per_m3: Some(dust * activity),
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn humidity_follows_the_dew_point() {
        // At its own dew point a surface is saturated; 10 K warmer, Magnus gives about 52 %
        // at 10 °C dew point (6.1094·e^(17.625·10/253.04) / 6.1094·e^(17.625·20/263.04)).
        assert!((relative_humidity_pct(12.0, 12.0) - 100.0).abs() < 1e-9);
        let rh = relative_humidity_pct(10.0, 20.0);
        let expected = 100.0 * (17.625 * 10.0 / 253.04 - 17.625 * 20.0 / 263.04f64).exp();
        assert!((rh - expected).abs() < 1e-9 && (52.0..53.0).contains(&rh), "{rh}");
        // Saturation at 20 °C: 23.37 hPa in the tables (WMO); the Magnus form is good to 0.4 %.
        assert!((saturation_hpa(20.0) / 23.37 - 1.0).abs() < 4e-3, "{}", saturation_hpa(20.0));
    }

    #[test]
    fn surroundings_open_downwind_and_above_with_a_floor() {
        use crate::case_scene::RoleInput;
        let s = ModelSurroundings { breeze_m_per_s: 0.2, breeze_from: BreezeFrom::MinusX, floor: true };
        let roles = s.roles();
        assert!(matches!(roles[0], RoleInput::Inlet { speed } if speed == 0.2));
        assert!(matches!(roles[1], RoleInput::Outlet { .. }));
        assert!(matches!(roles[2], RoleInput::Symmetry) && matches!(roles[3], RoleInput::Symmetry));
        assert!(matches!(roles[4], RoleInput::Wall) && matches!(roles[5], RoleInput::Outlet { .. }));
        let still = ModelSurroundings { breeze_m_per_s: 0.0, breeze_from: BreezeFrom::PlusY, floor: false };
        assert!(still.roles().iter().all(|r| matches!(r, RoleInput::Outlet { .. })));
    }

    #[test]
    fn hostile_climates_and_dust_are_refused() {
        let c = |name: &str, hours: f64, ambient: f64, dew: f64| ClimateCondition {
            name: name.into(),
            hours_per_year: hours,
            ambient_c: ambient,
            dew_point_c: dew,
            breeze_m_per_s: None,
            dust_mg_per_m3: None,
        };
        let year = |conditions| ModelClimate { conditions, humid_limit_pct: 60.0 };
        assert!(year(vec![c("a", 4000.0, 20.0, 10.0), c("b", 4760.0, 5.0, 0.0)]).validate().is_ok());
        assert!(year(vec![c("a", 9000.0, 20.0, 10.0)]).validate().is_err());
        assert!(year(vec![c("a", 100.0, 20.0, 25.0)]).validate().is_err());
        assert!(year(vec![c("a", 100.0, f64::NAN, 0.0)]).validate().is_err());
        assert!(year(vec![c("a", 100.0, 20.0, 0.0), c("A", 100.0, 20.0, 0.0)]).validate().is_err());
        assert!(year(vec![c("x/y", 100.0, 20.0, 0.0)]).validate().is_err());
        assert!(year(vec![]).validate().is_err());
        let dust = |concentration: f64, sizes: Vec<DustSize>| ModelDust {
            concentration_mg_per_m3: concentration,
            sizes,
            density_kg_per_m3: 2000.0,
            particles_per_size: 100,
            entry: DustEntry::Air,
            max_time_s: 100.0,
        };
        assert!(dust(1.0, default_sizes()).validate().is_ok());
        assert!(dust(-1.0, default_sizes()).validate().is_err());
        assert!(dust(1.0, vec![]).validate().is_err());
        assert!(dust(1.0, vec![DustSize { diameter_um: 5.0, share: 0.0 }]).validate().is_err());
        assert!(dust(1.0, vec![DustSize { diameter_um: f64::INFINITY, share: 1.0 }]).validate().is_err());
    }

    #[test]
    fn dust_block_takes_the_shares_and_the_box_into_the_domain() {
        let dust = ModelDust {
            concentration_mg_per_m3: 2.0,
            sizes: vec![DustSize { diameter_um: 1.0, share: 1.0 }, DustSize { diameter_um: 10.0, share: 3.0 }],
            density_kg_per_m3: 2600.0,
            particles_per_size: 10,
            entry: DustEntry::Box { low_mm: [10.0, 20.0, 30.0], high_mm: [40.0, 50.0, 60.0] },
            max_time_s: 60.0,
        };
        let block = dust.case_block([-10.0, 0.0, 10.0], [0.0, 0.0, -9.81], 290.0);
        assert!((block["classes"][0]["concentration_kg_per_m3"].as_f64().unwrap() - 0.5e-6).abs() < 1e-18);
        assert!((block["classes"][1]["concentration_kg_per_m3"].as_f64().unwrap() - 1.5e-6).abs() < 1e-18);
        assert_eq!(block["release"], "box");
        assert_eq!(block["release_box"]["min_m"], json!([0.02, 0.02, 0.02]));
    }

    #[test]
    fn a_year_adds_each_condition_by_its_hours() {
        // Two conditions, each with one deposit on the board at 1e-12 kg/s: 3 000 and 5 000 h.
        let report =
            |entering: f64| json!({"schema": "cfd-solvers-rs/dust/1", "amount_unit": "kg/s", "entering": entering, "deposited": 1e-12}).to_string();
        let deposits = |t_k: f64| {
            format!("class,diameter_m,x_m,y_m,z_m,target,amount,surface_temperature_k,time_s\n0,1e-5,0.01,0.02,0.03,Board,1e-12,{t_k},5\n")
        };
        let thermal =
            |min_k: f64| json!({"components": [{"name": "Case", "min_temperature_k": min_k, "max_temperature_k": min_k + 1.0}]}).to_string();
        let files = |t: f64, entering: f64| {
            vec![
                crate::Artifact { name: "dust_deposition.json".into(), text: report(entering) },
                crate::Artifact { name: "dust_deposits.csv".into(), text: deposits(t) },
                crate::Artifact { name: "thermal_report.json".into(), text: thermal(t) },
            ]
        };
        let warm = ClimateCondition {
            name: "warm".into(),
            hours_per_year: 3000.0,
            ambient_c: 25.0,
            dew_point_c: 5.0,
            breeze_m_per_s: None,
            dust_mg_per_m3: None,
        };
        let damp = ClimateCondition { name: "damp".into(), hours_per_year: 5000.0, ambient_c: 5.0, dew_point_c: 4.0, ..warm.clone() };
        let (a, b) = (files(298.15, 2e-12), files(278.15, 4e-12));
        let runs = [
            ConditionRun { condition: &warm, study: "s1".into(), run: Some(("r1".into(), &a)), origin_mm: [0.0; 3] },
            ConditionRun { condition: &damp, study: "s2".into(), run: Some(("r2".into(), &b)), origin_mm: [0.0; 3] },
        ];
        let s = summarize(&runs, 60.0).unwrap();
        let g = |hours: f64| 1e-12 * hours * 3600.0 * 1000.0;
        assert!((s.dust_in_g_per_year - (g(3000.0) * 2.0 + g(5000.0) * 4.0)).abs() < 1e-15);
        let board = s.parts.iter().find(|p| p.body == "Board").unwrap();
        assert!((board.dust_g_per_year - g(8000.0)).abs() < 1e-15);
        // Humid only in the damp condition (93 %); dry in the warm one (27 %): by mass, the
        // board's dust is humid for 5/8 of its mass times 5 000 h.
        assert!((board.humid_dust_hours - 5000.0 * 5.0 / 8.0).abs() < 1e-6, "{}", board.humid_dust_hours);
        let case = s.parts.iter().find(|p| p.body == "Case").unwrap();
        assert_eq!((case.humid_hours, case.condensing_hours), (5000.0, 0.0));
        assert_eq!(s.deposits_csv.lines().count(), 3);
        assert!(s.missing.is_empty() && (s.hours_covered - 8000.0).abs() < 1e-9);
    }

    #[test]
    fn four_seasons_cover_the_year() {
        let year = four_seasons(&[("winter", -4.0, 3.0, -5.0, 0.1), ("summer", 14.0, 30.0, 13.0, 0.8)]);
        assert_eq!(year.len(), 8);
        let hours: f64 = year.iter().map(|c| c.hours_per_year).sum();
        assert!((hours - YEAR_HOURS).abs() < 1e-6);
        let afternoon = year.iter().find(|c| c.name == "summer afternoon").unwrap();
        assert!((afternoon.ambient_c - 30.0).abs() < 1e-9 && afternoon.dew_point_c == 13.0);
        assert!(ModelClimate { conditions: year, humid_limit_pct: 60.0 }.validate().is_ok());
    }
}
