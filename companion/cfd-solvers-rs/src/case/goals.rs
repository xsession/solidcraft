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
}

impl GoalKind {
    fn thermal(self) -> bool {
        matches!(
            self,
            Self::PatchHeatOut | Self::SolidMaxTemperature | Self::SolidMeanTemperature
        )
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Goal {
    pub name: String,
    pub kind: GoalKind,
    /// Patch name, or solid name for solid goals.
    pub target: String,
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
            }
        }
        let n = case.goals.len();
        Ok(Some(Self {
            goals: case.goals.clone(),
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
            };
            if let Some(series) = self.values.get_mut(gi) {
                series.push(value);
                if let [.., previous, current] = series.as_slice() {
                    let base = match g.kind {
                        GoalKind::SolidMaxTemperature | GoalKind::SolidMeanTemperature => {
                            (current - self.inflow_k).abs()
                        }
                        _ => current.abs(),
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
                    "unit": match g.kind {
                        GoalKind::PatchVolumeFlow => "m3/s",
                        GoalKind::PatchMeanPressure => "Pa",
                        GoalKind::PatchHeatOut => "W",
                        _ => "K",
                    }
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
