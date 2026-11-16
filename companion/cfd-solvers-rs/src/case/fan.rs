//! External fan pressure-flow curve on a planar inlet patch (`fan_curve_v1`).
//!
//! The fan supplies a uniform inward normal velocity on one exterior patch. Its volume
//! flow is not prescribed: after each step the system pressure rise between the fan patch
//! and an ambient reference patch is measured and the flow is relaxed towards the curve's
//! flow at that pressure. The curve is never extrapolated; an operating point outside the
//! supplied range is clamped and reported. This is a lumped external-fan model, not a
//! resolved rotor, swirl model or internal momentum-source fan.
use super::CollocatedFlowCase;
use crate::{
    core::Vec3,
    fvm::{PolyMesh, PressureBoundaryType},
    require, Result, SolverError,
};
use serde::{Deserialize, Serialize};

pub const MAX_CURVE_POINTS: usize = 64;

/// One manufacturer curve sample in SI units.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct FanCurvePoint {
    pub volume_flow_m3_per_s: f64,
    pub pressure_rise_pa: f64,
}

/// Whether the curve's pressure rise is static or total (static plus fan-face dynamic pressure).
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FanPressureConvention {
    Static,
    Total,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct FanInlet {
    /// Planar exterior patch through which the fan blows inward. The fan owns its velocity.
    pub patch: String,
    /// Patch with a fixed-value kinematic pressure representing the ambient reference.
    pub reference_patch: String,
    pub pressure_convention: FanPressureConvention,
    /// Strictly increasing flow with strictly decreasing pressure rise.
    pub curve: Vec<FanCurvePoint>,
    /// Starting flow; defaults to the middle of the curve's flow range.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial_volume_flow_m3_per_s: Option<f64>,
    /// Fraction of the curve correction applied after each step, in (0, 1].
    #[serde(default = "default_relaxation")]
    pub relaxation: f64,
    /// Steps between operating-point corrections, so the field can respond to each change.
    #[serde(default = "default_update_interval")]
    pub update_interval: usize,
    /// A correction waits until the per-step relative change of the measured system
    /// pressure is at most this value, so a lagging pressure field cannot drive the fan.
    #[serde(default = "default_settle_tolerance")]
    pub settle_tolerance: f64,
}

fn default_relaxation() -> f64 {
    0.7
}

fn default_update_interval() -> usize {
    25
}

fn default_settle_tolerance() -> f64 {
    2e-3
}

impl FanInlet {
    pub fn validate(&self) -> Result<()> {
        require(
            (2..=MAX_CURVE_POINTS).contains(&self.curve.len()),
            format!("a fan curve needs 2..{MAX_CURVE_POINTS} points"),
        )?;
        require(
            self.curve.iter().all(|p| {
                p.volume_flow_m3_per_s.is_finite()
                    && p.pressure_rise_pa.is_finite()
                    && p.volume_flow_m3_per_s >= 0.0
                    && p.pressure_rise_pa >= 0.0
            }),
            "fan curve points must be finite and non-negative",
        )?;
        require(
            self.curve.windows(2).all(|w| {
                w[1].volume_flow_m3_per_s > w[0].volume_flow_m3_per_s
                    && w[1].pressure_rise_pa < w[0].pressure_rise_pa
            }),
            "fan curve flow must strictly increase while pressure rise strictly decreases",
        )?;
        require(
            self.relaxation.is_finite() && self.relaxation > 0.0 && self.relaxation <= 1.0,
            "fan relaxation must be in (0, 1]",
        )?;
        require(
            (1..=10_000).contains(&self.update_interval),
            "fan update interval must be 1..10000 steps",
        )?;
        require(
            self.settle_tolerance.is_finite()
                && self.settle_tolerance > 0.0
                && self.settle_tolerance <= 1.0,
            "fan settle tolerance must be in (0, 1]",
        )?;
        require(
            self.patch != self.reference_patch,
            "fan and reference patches must differ",
        )?;
        if let Some(q) = self.initial_volume_flow_m3_per_s {
            let (lo, hi) = self.flow_range();
            require(
                q.is_finite() && q >= lo && q <= hi,
                "initial fan flow must lie within the curve's flow range",
            )?;
        }
        Ok(())
    }

    fn flow_range(&self) -> (f64, f64) {
        let lo = self.curve.first().map_or(0.0, |p| p.volume_flow_m3_per_s);
        let hi = self.curve.last().map_or(0.0, |p| p.volume_flow_m3_per_s);
        (lo, hi)
    }

    /// Curve pressure rise at a flow inside the curve range (clamped at the ends).
    pub fn pressure_at(&self, q: f64) -> f64 {
        let (Some(first), Some(last)) = (self.curve.first(), self.curve.last()) else {
            return 0.0;
        };
        if q <= first.volume_flow_m3_per_s {
            return first.pressure_rise_pa;
        }
        if q >= last.volume_flow_m3_per_s {
            return last.pressure_rise_pa;
        }
        for w in self.curve.windows(2) {
            let (a, b) = (w[0], w[1]);
            if q <= b.volume_flow_m3_per_s {
                let t = (q - a.volume_flow_m3_per_s)
                    / (b.volume_flow_m3_per_s - a.volume_flow_m3_per_s);
                return a.pressure_rise_pa + t * (b.pressure_rise_pa - a.pressure_rise_pa);
            }
        }
        last.pressure_rise_pa
    }

    /// Intersection of the curve with a resistance line `pressure = k * flow`.
    /// The curve decreases and the line increases, so the crossing is unique when it exists.
    pub fn intersect_resistance(&self, k: f64) -> f64 {
        let (Some(first), Some(last)) = (self.curve.first(), self.curve.last()) else {
            return 0.0;
        };
        if k == f64::INFINITY {
            return first.volume_flow_m3_per_s;
        }
        if !(k.is_finite() && k > 0.0) {
            // No measured resistance (or assisting pressure): the fan reaches free delivery.
            return last.volume_flow_m3_per_s;
        }
        if first.pressure_rise_pa <= k * first.volume_flow_m3_per_s {
            return first.volume_flow_m3_per_s;
        }
        for w in self.curve.windows(2) {
            let (a, b) = (w[0], w[1]);
            let fa = a.pressure_rise_pa - k * a.volume_flow_m3_per_s;
            let fb = b.pressure_rise_pa - k * b.volume_flow_m3_per_s;
            if fb <= 0.0 {
                let t = fa / (fa - fb);
                return a.volume_flow_m3_per_s
                    + t * (b.volume_flow_m3_per_s - a.volume_flow_m3_per_s);
            }
        }
        last.volume_flow_m3_per_s
    }

    /// Flow delivered at a required pressure rise and whether that pressure is on the curve.
    pub fn flow_at(&self, pressure: f64) -> (f64, bool) {
        let (Some(first), Some(last)) = (self.curve.first(), self.curve.last()) else {
            return (0.0, false);
        };
        if pressure > first.pressure_rise_pa {
            return (first.volume_flow_m3_per_s, false);
        }
        if pressure < last.pressure_rise_pa {
            return (last.volume_flow_m3_per_s, false);
        }
        for w in self.curve.windows(2) {
            let (a, b) = (w[0], w[1]);
            if pressure >= b.pressure_rise_pa {
                let t = (a.pressure_rise_pa - pressure) / (a.pressure_rise_pa - b.pressure_rise_pa);
                return (
                    a.volume_flow_m3_per_s + t * (b.volume_flow_m3_per_s - a.volume_flow_m3_per_s),
                    true,
                );
            }
        }
        (last.volume_flow_m3_per_s, true)
    }
}

/// One controller sample, recorded after the step whose pressure field it measured.
#[derive(Clone, Copy, Debug, Serialize, PartialEq)]
pub struct FanSample {
    pub step: usize,
    pub volume_flow_m3_per_s: f64,
    pub system_pressure_rise_pa: f64,
    pub curve_pressure_rise_pa: f64,
    pub within_curve: bool,
}

pub(super) struct FanController {
    fan: FanInlet,
    faces: Vec<usize>,
    area: f64,
    inward: Vec3,
    reference_pressure: f64,
    density: f64,
    flow: f64,
    next_update: usize,
    corrections: usize,
    pub history: Vec<FanSample>,
}

impl FanController {
    pub fn new(case: &CollocatedFlowCase, fan: &FanInlet, mesh: &PolyMesh) -> Result<Self> {
        fan.validate()?;
        require(
            !case
                .velocity_boundaries
                .iter()
                .any(|b| b.patch == fan.patch),
            "the fan owns its patch velocity; remove that patch from velocity_boundaries",
        )?;
        require(
            !case
                .pressure_boundaries
                .iter()
                .any(|b| b.patch == fan.patch && b.kind == PressureBoundaryType::FixedValue),
            "a fan patch cannot also have a fixed pressure",
        )?;
        let reference = case
            .pressure_boundaries
            .iter()
            .find(|b| b.patch == fan.reference_patch && b.kind == PressureBoundaryType::FixedValue)
            .ok_or_else(|| {
                SolverError::InvalidInput(
                    "the fan reference patch needs a fixed-value pressure boundary".into(),
                )
            })?;
        require(
            !case
                .obstacles
                .iter()
                .any(|o| o.name == fan.patch || o.name == fan.reference_patch),
            "fan and reference patches must be exterior patches",
        )?;
        let patch = mesh.patch_index(&fan.patch)?;
        let faces: Vec<usize> = mesh
            .faces()
            .iter()
            .enumerate()
            .filter(|(_, f)| f.patch == Some(patch))
            .map(|(i, _)| i)
            .collect();
        require(
            !faces.is_empty(),
            "the fan patch has no exposed fluid faces",
        )?;
        let mut sum = Vec3::ZERO;
        let mut area = 0.0;
        for &fi in &faces {
            let a = mesh.faces()[fi].area;
            sum += a;
            area += a.norm();
        }
        require(
            area > 0.0 && area.is_finite(),
            "fan patch area must be positive",
        )?;
        // A curved or folded patch would make one inward speed ambiguous.
        require(
            sum.norm() >= area * (1.0 - 1e-9),
            "the fan patch must be planar",
        )?;
        let inward = -(sum / sum.norm());
        let (lo, hi) = fan.flow_range();
        let flow = fan.initial_volume_flow_m3_per_s.unwrap_or(0.5 * (lo + hi));
        Ok(Self {
            fan: fan.clone(),
            faces,
            area,
            inward,
            reference_pressure: reference.value,
            density: case.config.density,
            flow,
            next_update: fan.update_interval,
            corrections: 0,
            history: vec![],
        })
    }

    pub fn patch(&self) -> &str {
        &self.fan.patch
    }

    pub fn velocity(&self) -> Vec3 {
        self.inward * (self.flow / self.area)
    }

    /// Pressure rise the fan must provide for the current flow and pressure field, in Pa.
    fn system_pressure(&self, mesh: &PolyMesh, kinematic_pressure: &[f64]) -> Result<f64> {
        let mut weighted = 0.0;
        for &fi in &self.faces {
            let f = mesh.faces()[fi];
            let p = kinematic_pressure
                .get(f.owner)
                .copied()
                .ok_or_else(|| SolverError::InvalidInput("fan pressure sample".into()))?;
            weighted += f.area.norm() * p;
        }
        let static_rise = self.density * (weighted / self.area - self.reference_pressure);
        let speed = self.flow / self.area;
        Ok(match self.fan.pressure_convention {
            FanPressureConvention::Static => static_rise,
            FanPressureConvention::Total => static_rise + 0.5 * self.density * speed * speed,
        })
    }

    /// Record the operating state for the latest field and, at each update interval,
    /// relax the flow towards the curve. Returns a new inlet velocity when it changes.
    pub fn update(
        &mut self,
        step: usize,
        mesh: &PolyMesh,
        kinematic_pressure: &[f64],
    ) -> Result<Option<Vec3>> {
        let system = self.system_pressure(mesh, kinematic_pressure)?;
        require(system.is_finite(), "nonfinite fan system pressure")?;
        let (direct, within) = self.fan.flow_at(system);
        // The measured point defines the current system resistance; its crossing with the
        // curve is a better-conditioned target than inverting the curve at a lagging pressure.
        let target = if self.flow > 0.0 {
            self.fan.intersect_resistance(system / self.flow)
        } else {
            direct
        };
        self.history.push(FanSample {
            step,
            volume_flow_m3_per_s: self.flow,
            system_pressure_rise_pa: system,
            curve_pressure_rise_pa: self.fan.pressure_at(self.flow),
            within_curve: within,
        });
        let previous = self
            .history
            .iter()
            .rev()
            .nth(1)
            .map(|s| s.system_pressure_rise_pa);
        let settled = previous.is_some_and(|p| {
            (system - p).abs() <= self.fan.settle_tolerance * system.abs().max(p.abs())
        });
        if step < self.next_update || !settled {
            return Ok(None);
        }
        self.next_update = step.saturating_add(self.fan.update_interval);
        self.corrections += 1;
        let (lo, hi) = self.fan.flow_range();
        self.flow = (self.flow + self.fan.relaxation * (target - self.flow)).clamp(lo, hi);
        Ok(Some(self.velocity()))
    }

    /// Final operating-point summary, measured on the last solved field.
    pub fn summary(&self) -> serde_json::Value {
        let last = self.history.last().copied();
        let (q, system, curve, within) = last.map_or((self.flow, 0.0, 0.0, false), |s| {
            (
                s.volume_flow_m3_per_s,
                s.system_pressure_rise_pa,
                s.curve_pressure_rise_pa,
                s.within_curve,
            )
        });
        let scale = curve.abs().max(system.abs());
        let mismatch = if scale > 0.0 {
            (curve - system).abs() / scale
        } else {
            0.0
        };
        serde_json::json!({
            "schema": "cfd-solvers-rs/fan-operating-point/1",
            "model": "external_fan_curve_uniform_normal_inlet",
            "patch": self.fan.patch,
            "reference_patch": self.fan.reference_patch,
            "pressure_convention": self.fan.pressure_convention,
            "curve": self.fan.curve,
            "relaxation": self.fan.relaxation,
            "update_interval": self.fan.update_interval,
            "settle_tolerance": self.fan.settle_tolerance,
            "corrections": self.corrections,
            "patch_area_m2": self.area,
            "inward_unit_normal": self.inward,
            "volume_flow_m3_per_s": q,
            "mean_normal_speed_m_per_s": q / self.area,
            "system_pressure_rise_pa": system,
            "curve_pressure_rise_pa": curve,
            "relative_operating_point_mismatch": mismatch,
            "within_curve_range": within,
            "samples": self.history.len(),
            "extrapolation": "refused_clamped_to_curve_range",
            "physical_validation": "not_established"
        })
    }
}

pub(super) fn write_report(
    controller: &FanController,
    out: &std::path::Path,
    m: &mut super::RunManifest,
) -> Result<()> {
    let summary = controller.summary();
    std::fs::write(
        out.join("fan_operating_point.json"),
        serde_json::to_vec_pretty(&summary)?,
    )?;
    m.files.push("fan_operating_point.json".into());
    let column = |f: fn(&FanSample) -> f64| controller.history.iter().map(f).collect::<Vec<f64>>();
    let step = column(|s| s.step as f64);
    let flow = column(|s| s.volume_flow_m3_per_s);
    let system = column(|s| s.system_pressure_rise_pa);
    let curve = column(|s| s.curve_pressure_rise_pa);
    crate::io::write_csv(
        out.join("fan_operating_point.csv"),
        &[
            ("step", &step),
            ("volume_flow_m3_per_s", &flow),
            ("system_pressure_rise_pa", &system),
            ("curve_pressure_rise_pa", &curve),
        ],
    )?;
    m.files.push("fan_operating_point.csv".into());
    for (key, name) in [
        ("volume_flow_m3_per_s", "fan_volume_flow_m3_per_s"),
        ("system_pressure_rise_pa", "fan_system_pressure_rise_pa"),
        ("curve_pressure_rise_pa", "fan_curve_pressure_rise_pa"),
        (
            "relative_operating_point_mismatch",
            "fan_relative_operating_point_mismatch",
        ),
    ] {
        if let Some(value) = summary.get(key).and_then(serde_json::Value::as_f64) {
            m.metrics.insert(name.into(), value);
        }
    }
    let within = summary
        .get("within_curve_range")
        .and_then(serde_json::Value::as_bool)
        == Some(true);
    m.metrics.insert(
        "fan_operating_point_within_curve".into(),
        if within { 1.0 } else { 0.0 },
    );
    Ok(())
}
