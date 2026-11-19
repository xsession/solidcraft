//! Engineering goals for flow runs (`goals_v1`): named monitored quantities with a declared
//! convergence test and optional early stop.
//!
//! A goal is converged when its value changed by at most `relative_tolerance` of its scale
//! between the last two evaluation windows. For temperatures the scale is the rise above the
//! inflow temperature, so a 0.1% tolerance means 0.1% of the heating, not of 300 K. Goal
//! convergence describes these monitored quantities only. It is not mesh independence and not
//! physical validation.
use super::conjugate::EnergyState;
use super::RunManifest;
use crate::{
    fvm::{PolyMesh, PressureBoundaryType},
    require, Result,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, path::Path};

pub const MAX_GOALS: usize = 32;

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GoalKind {
    /// Outward volume flow through a boundary patch, m^3/s.
    PatchVolumeFlow,
    /// Area-weighted static pressure on a boundary patch, Pa (kinematic pressure times density).
    PatchMeanPressure,
    /// Conductive heat leaving through an exterior patch, W (thermal cases).
    PatchHeatOut,
    /// Maximum temperature of a solid, K (thermal cases).
    SolidMaxTemperature,
    /// Volume-mean temperature of a solid, K (thermal cases).
    SolidMeanTemperature,
    /// Maximum temperature of all fluid cells, K (thermal cases; `goals_v2`).
    GlobalMaxFluidTemperature,
    /// Volume-mean temperature of all fluid cells, K (thermal cases; `goals_v2`).
    GlobalMeanFluidTemperature,
    /// Maximum temperature of all solid cells, K (thermal cases; `goals_v2`).
    GlobalMaxSolidTemperature,
    /// Temperature of the cell nearest to `point`, K (thermal cases; `goals_v2`).
    PointTemperature,
    /// Static pressure of the fluid cell nearest to `point`, Pa (`goals_v2`).
    PointPressure,
    /// Arithmetic of other goals' current values (`goals_v2`): `+ - * /`, parentheses,
    /// numbers and goal names, e.g. `chip_max - inlet_temperature`.
    Equation,
}

impl GoalKind {
    fn thermal(self) -> bool {
        matches!(
            self,
            Self::PatchHeatOut
                | Self::SolidMaxTemperature
                | Self::SolidMeanTemperature
                | Self::GlobalMaxFluidTemperature
                | Self::GlobalMeanFluidTemperature
                | Self::GlobalMaxSolidTemperature
                | Self::PointTemperature
        )
    }

    /// Temperatures converge relative to their rise above the inflow temperature.
    fn temperature(self) -> bool {
        self.thermal() && self != Self::PatchHeatOut
    }

    fn unit(self) -> &'static str {
        match self {
            Self::PatchVolumeFlow => "m3/s",
            Self::PatchMeanPressure | Self::PointPressure => "Pa",
            Self::PatchHeatOut => "W",
            Self::Equation => "",
            _ => "K",
        }
    }
}

/// A parsed equation goal.
#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    Number(f64),
    Goal(String),
    Negate(Box<Expr>),
    Binary(char, Box<Expr>, Box<Expr>),
}

impl Expr {
    /// Parse `+ - * /`, unary minus, parentheses, numbers and goal names.
    pub fn parse(text: &str) -> Result<Self> {
        let tokens = tokenize(text)?;
        let mut at = 0;
        let expr = parse_sum(&tokens, &mut at, 0)?;
        require(at == tokens.len(), format!("unexpected text in equation: {text}"))?;
        Ok(expr)
    }

    /// The goal names the expression refers to.
    pub fn names(&self, out: &mut Vec<String>) {
        match self {
            Self::Number(_) => {}
            Self::Goal(name) => out.push(name.clone()),
            Self::Negate(inner) => inner.names(out),
            Self::Binary(_, a, b) => {
                a.names(out);
                b.names(out);
            }
        }
    }

    /// Evaluate with `value` giving each goal's current value (NaN propagates).
    pub fn eval(&self, value: &dyn Fn(&str) -> f64) -> f64 {
        match self {
            Self::Number(n) => *n,
            Self::Goal(name) => value(name),
            Self::Negate(inner) => -inner.eval(value),
            Self::Binary(op, a, b) => {
                let (a, b) = (a.eval(value), b.eval(value));
                match op {
                    '+' => a + b,
                    '-' => a - b,
                    '*' => a * b,
                    _ => a / b,
                }
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
enum Token {
    Number(f64),
    Name(String),
    Op(char),
}

fn tokenize(text: &str) -> Result<Vec<Token>> {
    require(text.len() <= 400, "equations are at most 400 characters")?;
    let chars: Vec<char> = text.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while let Some(&c) = chars.get(i) {
        if c.is_whitespace() {
            i += 1;
        } else if "+-*/()".contains(c) {
            out.push(Token::Op(c));
            i += 1;
        } else if c.is_ascii_digit() || c == '.' {
            let start = i;
            while chars.get(i).is_some_and(|d| d.is_ascii_digit() || *d == '.' || *d == 'e' || *d == 'E'
                || ((*d == '-' || *d == '+') && matches!(chars.get(i.wrapping_sub(1)), Some('e' | 'E'))))
            {
                i += 1;
            }
            let literal: String = chars[start..i].iter().collect();
            let value: f64 = literal
                .parse()
                .map_err(|_| crate::SolverError::InvalidInput(format!("bad number {literal} in equation")))?;
            out.push(Token::Number(value));
        } else if c.is_ascii_alphabetic() || c == '_' {
            let start = i;
            while chars.get(i).is_some_and(|d| d.is_ascii_alphanumeric() || *d == '_') {
                i += 1;
            }
            out.push(Token::Name(chars[start..i].iter().collect()));
        } else {
            return Err(crate::SolverError::InvalidInput(format!("unexpected {c:?} in equation")));
        }
    }
    require(!out.is_empty(), "empty equation")?;
    Ok(out)
}

fn parse_sum(tokens: &[Token], at: &mut usize, depth: usize) -> Result<Expr> {
    let mut left = parse_product(tokens, at, depth)?;
    while let Some(Token::Op(op @ ('+' | '-'))) = tokens.get(*at) {
        *at += 1;
        let right = parse_product(tokens, at, depth)?;
        left = Expr::Binary(*op, Box::new(left), Box::new(right));
    }
    Ok(left)
}

fn parse_product(tokens: &[Token], at: &mut usize, depth: usize) -> Result<Expr> {
    let mut left = parse_unary(tokens, at, depth)?;
    while let Some(Token::Op(op @ ('*' | '/'))) = tokens.get(*at) {
        *at += 1;
        let right = parse_unary(tokens, at, depth)?;
        left = Expr::Binary(*op, Box::new(left), Box::new(right));
    }
    Ok(left)
}

fn parse_unary(tokens: &[Token], at: &mut usize, depth: usize) -> Result<Expr> {
    require(depth < 32, "equation nests too deeply")?;
    match tokens.get(*at) {
        Some(Token::Op('-')) => {
            *at += 1;
            Ok(Expr::Negate(Box::new(parse_unary(tokens, at, depth + 1)?)))
        }
        Some(Token::Op('(')) => {
            *at += 1;
            let inner = parse_sum(tokens, at, depth + 1)?;
            require(tokens.get(*at) == Some(&Token::Op(')')), "missing ) in equation")?;
            *at += 1;
            Ok(inner)
        }
        Some(Token::Number(n)) => {
            *at += 1;
            Ok(Expr::Number(*n))
        }
        Some(Token::Name(name)) => {
            *at += 1;
            Ok(Expr::Goal(name.clone()))
        }
        _ => Err(crate::SolverError::InvalidInput("incomplete equation".into())),
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Goal {
    pub name: String,
    pub kind: GoalKind,
    /// Patch name, or solid name for solid goals; unused by global, point and equation goals.
    #[serde(default)]
    pub target: String,
    /// Location of a point goal, m.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub point: Option<[f64; 3]>,
    /// Expression of an equation goal.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expression: Option<String>,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GoalControl {
    /// Steps between evaluations that are compared.
    pub window_steps: usize,
    pub relative_tolerance: f64,
    /// Stop the run once every goal has converged (after `min_steps`).
    #[serde(default)]
    pub stop_when_converged: bool,
    #[serde(default)]
    pub min_steps: usize,
}

impl Default for GoalControl {
    fn default() -> Self {
        Self {
            window_steps: 50,
            relative_tolerance: 1e-3,
            stop_when_converged: false,
            min_steps: 0,
        }
    }
}

pub(super) struct GoalMonitor {
    goals: Vec<Goal>,
    /// Parsed expressions of equation goals, by goal index.
    equations: Vec<Option<Expr>>,
    control: GoalControl,
    inflow_k: f64,
    density: f64,
    pub steps: Vec<usize>,
    pub values: Vec<Vec<f64>>,
    pub converged: Vec<bool>,
    pub last_change: Vec<f64>,
    pub stopped_at: Option<usize>,
}

impl GoalMonitor {
    pub fn new(case: &super::CollocatedFlowCase, mesh: &PolyMesh) -> Result<Option<Self>> {
        if case.goals.is_empty() {
            return Ok(None);
        }
        let control = case.goal_control.unwrap_or_default();
        require(
            case.goals.len() <= MAX_GOALS,
            format!("at most {MAX_GOALS} goals"),
        )?;
        require(
            (1..=10_000).contains(&control.window_steps)
                && control.relative_tolerance.is_finite()
                && control.relative_tolerance > 0.0
                && control.relative_tolerance < 1.0,
            "goal window must be 1..10000 steps and tolerance in (0, 1)",
        )?;
        let solids = case.solid_names();
        let mut names = BTreeSet::new();
        for g in &case.goals {
            let n = g.name.as_bytes();
            require(
                !n.is_empty()
                    && n.len() <= 60
                    && n.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_')
                    && names.insert(g.name.as_str()),
                "goal names must be unique identifiers of letters, digits and underscores",
            )?;
            require(
                !g.kind.thermal() || case.thermal.is_some(),
                format!("goal {} needs a thermal block", g.name),
            )?;
            match g.kind {
                GoalKind::SolidMaxTemperature | GoalKind::SolidMeanTemperature => {
                    require(
                        solids.contains(&g.target.as_str()),
                        format!("goal {} names unknown solid {}", g.name, g.target),
                    )?;
                }
                GoalKind::PatchHeatOut => {
                    require(
                        ["left", "right", "bottom", "top", "front", "back"]
                            .contains(&g.target.as_str()),
                        format!("goal {} needs an exterior patch", g.name),
                    )?;
                }
                GoalKind::PatchVolumeFlow | GoalKind::PatchMeanPressure => {
                    mesh.patch_index(&g.target)?;
                }
                GoalKind::GlobalMaxFluidTemperature
                | GoalKind::GlobalMeanFluidTemperature
                | GoalKind::GlobalMaxSolidTemperature => {}
                GoalKind::PointTemperature | GoalKind::PointPressure => {
                    let lengths = [case.mesh.lengths.x, case.mesh.lengths.y, case.mesh.lengths.z];
                    require(
                        g.point.is_some_and(|p| {
                            p.iter().zip(lengths).all(|(v, l)| v.is_finite() && *v >= 0.0 && *v <= l)
                        }),
                        format!("point goal {} needs a point inside the domain", g.name),
                    )?;
                }
                GoalKind::Equation => {}
            }
        }
        let mut equations = Vec::with_capacity(case.goals.len());
        for g in &case.goals {
            if g.kind != GoalKind::Equation {
                equations.push(None);
                continue;
            }
            let text = g.expression.as_deref().unwrap_or_default();
            let expr = Expr::parse(text).map_err(|e| crate::SolverError::InvalidInput(format!("goal {}: {e}", g.name)))?;
            let mut refs = Vec::new();
            expr.names(&mut refs);
            for name in &refs {
                require(
                    case.goals.iter().any(|other| &other.name == name && other.kind != GoalKind::Equation),
                    format!("equation goal {} refers to {name}, which is not a non-equation goal", g.name),
                )?;
            }
            equations.push(Some(expr));
        }
        let n = case.goals.len();
        Ok(Some(Self {
            goals: case.goals.clone(),
            equations,
            control,
            inflow_k: case
                .thermal
                .as_ref()
                .map_or(0.0, |t| t.inflow_temperature_k),
            density: case.config.density,
            steps: vec![],
            values: vec![vec![]; n],
            converged: vec![false; n],
            last_change: vec![f64::NAN; n],
            stopped_at: None,
        }))
    }

    pub fn needs_thermal_refresh(&self, step: usize) -> bool {
        step % self.control.window_steps == 0 && self.goals.iter().any(|g| g.kind.thermal())
    }

    pub fn due(&self, step: usize) -> bool {
        step % self.control.window_steps == 0
    }

    /// Record goal values at an evaluation step; returns true when the run may stop.
    pub fn record(
        &mut self,
        step: usize,
        mesh: &PolyMesh,
        flux: &[f64],
        pressure: &[f64],
        pressure_fixed: &dyn Fn(usize) -> Option<f64>,
        energy: Option<&EnergyState>,
    ) -> Result<bool> {
        self.steps.push(step);
        for (gi, g) in self.goals.iter().enumerate() {
            let value = match g.kind {
                GoalKind::PatchVolumeFlow | GoalKind::PatchMeanPressure => {
                    let patch = mesh.patch_index(&g.target)?;
                    let (mut q, mut pa, mut area) = (0.0, 0.0, 0.0);
                    for (fi, f) in mesh
                        .faces()
                        .iter()
                        .enumerate()
                        .filter(|(_, f)| f.patch == Some(patch))
                    {
                        q += flux.get(fi).copied().unwrap_or(0.0);
                        let a = f.area.norm();
                        let p = pressure_fixed(patch)
                            .unwrap_or_else(|| pressure.get(f.owner).copied().unwrap_or(0.0));
                        pa += a * p;
                        area += a;
                    }
                    if g.kind == GoalKind::PatchVolumeFlow {
                        q
                    } else if area > 0.0 {
                        self.density * pa / area
                    } else {
                        0.0
                    }
                }
                GoalKind::PatchHeatOut => match energy {
                    Some(e) => e.patch_conduction(&g.target, flux, self.inflow_k)?,
                    None => f64::NAN,
                },
                GoalKind::SolidMaxTemperature | GoalKind::SolidMeanTemperature => {
                    match energy.and_then(|e| e.solid_temperature(&g.target, self.inflow_k)) {
                        Some((max, mean)) => {
                            if g.kind == GoalKind::SolidMaxTemperature {
                                max
                            } else {
                                mean
                            }
                        }
                        None => f64::NAN,
                    }
                }
                GoalKind::GlobalMaxFluidTemperature | GoalKind::GlobalMeanFluidTemperature => energy
                    .and_then(|e| e.region_temperature(false, self.inflow_k))
                    .map_or(f64::NAN, |(max, mean)| {
                        if g.kind == GoalKind::GlobalMaxFluidTemperature { max } else { mean }
                    }),
                GoalKind::GlobalMaxSolidTemperature => energy
                    .and_then(|e| e.region_temperature(true, self.inflow_k))
                    .map_or(f64::NAN, |(max, _)| max),
                GoalKind::PointTemperature => match (energy, g.point) {
                    (Some(e), Some(p)) => e.point_temperature(p, self.inflow_k).unwrap_or(f64::NAN),
                    _ => f64::NAN,
                },
                GoalKind::PointPressure => match g.point {
                    Some(p) => {
                        let at = crate::core::Vec3::new(p[0], p[1], p[2]);
                        (0..mesh.cell_count())
                            .min_by(|a, b| {
                                (mesh.cells()[*a].center - at)
                                    .norm()
                                    .total_cmp(&(mesh.cells()[*b].center - at).norm())
                            })
                            .and_then(|c| pressure.get(c))
                            .map_or(f64::NAN, |p| self.density * p)
                    }
                    None => f64::NAN,
                },
                // Evaluated below, once every other goal has its value for this step.
                GoalKind::Equation => f64::NAN,
            };
            if let Some(series) = self.values.get_mut(gi) {
                series.push(value);
            }
        }
        // Equation goals combine this step's values of the other goals.
        for gi in 0..self.goals.len() {
            let Some(Some(expr)) = self.equations.get(gi) else { continue };
            let current = |name: &str| {
                self.goals
                    .iter()
                    .position(|g| g.name == name)
                    .and_then(|i| self.values.get(i))
                    .and_then(|v| v.last())
                    .copied()
                    .unwrap_or(f64::NAN)
            };
            let value = expr.eval(&current);
            if let Some(slot) = self.values.get_mut(gi).and_then(|series| series.last_mut()) {
                *slot = value;
            }
        }
        for (gi, g) in self.goals.iter().enumerate() {
            if let Some(series) = self.values.get_mut(gi) {
                if let [.., previous, current] = series.as_slice() {
                    let base = if g.kind.temperature() {
                        (current - self.inflow_k).abs()
                    } else {
                        current.abs()
                    };
                    let change = (current - previous).abs() / base.max(f64::MIN_POSITIVE);
                    self.last_change[gi] = change;
                    self.converged[gi] =
                        change.is_finite() && change <= self.control.relative_tolerance;
                }
            }
        }
        let done = self.converged.iter().all(|c| *c);
        Ok(done && self.control.stop_when_converged && step >= self.control.min_steps)
    }

    pub fn write(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        let mut columns: Vec<(&str, Vec<f64>)> =
            vec![("step", self.steps.iter().map(|s| *s as f64).collect())];
        for (g, series) in self.goals.iter().zip(&self.values) {
            columns.push((g.name.as_str(), series.clone()));
        }
        let borrowed: Vec<(&str, &[f64])> =
            columns.iter().map(|(n, v)| (*n, v.as_slice())).collect();
        crate::io::write_csv(out.join("goals.csv"), &borrowed)?;
        m.files.push("goals.csv".into());
        let all = self.converged.iter().all(|c| *c);
        let goals: Vec<serde_json::Value> = self
            .goals
            .iter()
            .enumerate()
            .map(|(i, g)| {
                serde_json::json!({
                    "name": g.name, "kind": g.kind, "target": g.target,
                    "final_value": self.values.get(i).and_then(|v| v.last()).copied(),
                    "last_relative_change": self.last_change.get(i).copied().filter(|v| v.is_finite()),
                    "converged": self.converged.get(i).copied().unwrap_or(false),
                    "unit": g.kind.unit(),
                    "point": g.point,
                    "expression": g.expression
                })
            })
            .collect();
        let report = serde_json::json!({
            "schema": "cfd-solvers-rs/goals/1",
            "control": self.control,
            "all_converged": all,
            "stopped_early_at_step": self.stopped_at,
            "rule": "relative change between the last two evaluations; temperatures relative to their rise above inflow",
            "goals": goals,
            "meaning": "convergence of these monitored quantities only; not mesh independence or physical validation"
        });
        std::fs::write(out.join("goals.json"), serde_json::to_vec_pretty(&report)?)?;
        m.files.push("goals.json".into());
        m.metrics
            .insert("goals_converged".into(), if all { 1.0 } else { 0.0 });
        for (g, series) in self.goals.iter().zip(&self.values) {
            if let Some(v) = series.last().filter(|v| v.is_finite()) {
                m.metrics.insert(format!("goal_{}", g.name), *v);
            }
        }
        Ok(())
    }
}

pub(super) fn fixed_pressure_lookup(
    case: &super::CollocatedFlowCase,
    mesh: &PolyMesh,
) -> Vec<Option<f64>> {
    mesh.patches()
        .iter()
        .map(|p| {
            case.pressure_boundaries
                .iter()
                .find(|b| b.patch == p.name && b.kind == PressureBoundaryType::FixedValue)
                .map(|b| b.value)
        })
        .collect()
}
