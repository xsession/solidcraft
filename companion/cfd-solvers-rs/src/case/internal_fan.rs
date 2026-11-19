//! Internal fans for `collocated_flow3d` (`internal_fans_v1`): a block of fluid cells that
//! raises the pressure of the air passing through it along one grid direction, following a
//! manufacturer pressure-flow curve (an actuator-block model).
//!
//! Each step the fan's volume flow is measured from the conservative face fluxes: the mean,
//! over the block's length, of the flux along n through its cross-section planes (the cell
//! velocities of a collocated grid can carry checkerboard noise; the fluxes cannot). The pressure rise is relaxed towards
//! the curve's rise at that flow, Δp ← Δp + r (Δp_curve(Q) − Δp), and applied as the body
//! force Δp / (ρ L) along n in the block's cells. A fluid particle crossing the block gains
//! Δp. The curve is never extrapolated: a flow outside its range takes the end value, and
//! the report says so. Swirl, blade loading and the fan's own losses are not modelled.
use super::{fan::FanCurvePoint, CollocatedFlowCase, RunManifest};
use crate::{core::Vec3, fvm::PolyMesh, require, Result, SolverError};
use serde::{Deserialize, Serialize};
use std::path::Path;

pub const MAX_INTERNAL_FANS: usize = 16;

/// The grid direction a fan blows along.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum FanDirection {
    #[serde(rename = "+x")]
    PlusX,
    #[serde(rename = "-x")]
    MinusX,
    #[serde(rename = "+y")]
    PlusY,
    #[serde(rename = "-y")]
    MinusY,
    #[serde(rename = "+z")]
    PlusZ,
    #[serde(rename = "-z")]
    MinusZ,
}

impl FanDirection {
    pub fn axis(self) -> usize {
        match self {
            Self::PlusX | Self::MinusX => 0,
            Self::PlusY | Self::MinusY => 1,
            Self::PlusZ | Self::MinusZ => 2,
        }
    }
    pub fn sign(self) -> f64 {
        match self {
            Self::PlusX | Self::PlusY | Self::PlusZ => 1.0,
            _ => -1.0,
        }
    }
    pub fn unit(self) -> Vec3 {
        let mut v = [0.0; 3];
        v[self.axis()] = self.sign();
        Vec3::new(v[0], v[1], v[2])
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct InternalFan {
    pub name: String,
    /// Half-open background-grid cell range of the fan block (all fluid).
    pub min: [usize; 3],
    pub max: [usize; 3],
    pub direction: FanDirection,
    /// Strictly increasing flow with strictly decreasing pressure rise (pascals).
    pub curve: Vec<FanCurvePoint>,
    /// Fraction of the curve correction applied each step, in (0, 1].
    #[serde(default = "default_relaxation")]
    pub relaxation: f64,
}

fn default_relaxation() -> f64 {
    0.2
}

/// Curve pressure rise at `q`, clamped to the end values; and whether `q` is on the curve.
pub fn curve_pressure(curve: &[FanCurvePoint], q: f64) -> (f64, bool) {
    let (Some(first), Some(last)) = (curve.first(), curve.last()) else {
        return (0.0, false);
    };
    if q < first.volume_flow_m3_per_s {
        return (first.pressure_rise_pa, false);
    }
    if q > last.volume_flow_m3_per_s {
        return (last.pressure_rise_pa, false);
    }
    for w in curve.windows(2) {
        let (a, b) = (w[0], w[1]);
        if q <= b.volume_flow_m3_per_s {
            let t =
                (q - a.volume_flow_m3_per_s) / (b.volume_flow_m3_per_s - a.volume_flow_m3_per_s);
            return (
                a.pressure_rise_pa + t * (b.pressure_rise_pa - a.pressure_rise_pa),
                true,
            );
        }
    }
    (last.pressure_rise_pa, true)
}

/// The flow at which the curve meets a system resistance line Δp = k Q (k ≥ 0).
pub fn operating_point(curve: &[FanCurvePoint], k: f64) -> f64 {
    let (Some(first), Some(last)) = (curve.first(), curve.last()) else {
        return 0.0;
    };
    if first.pressure_rise_pa <= k * first.volume_flow_m3_per_s {
        return first.volume_flow_m3_per_s;
    }
    for w in curve.windows(2) {
        let (a, b) = (w[0], w[1]);
        let fa = a.pressure_rise_pa - k * a.volume_flow_m3_per_s;
        let fb = b.pressure_rise_pa - k * b.volume_flow_m3_per_s;
        if fb <= 0.0 {
            return a.volume_flow_m3_per_s
                + fa / (fa - fb) * (b.volume_flow_m3_per_s - a.volume_flow_m3_per_s);
        }
    }
    last.volume_flow_m3_per_s
}

impl InternalFan {
    pub fn validate(&self, case: &CollocatedFlowCase) -> Result<()> {
        let d = case.mesh.dimensions;
        require(
            (0..3).all(|a| self.min[a] < self.max[a] && self.max[a] <= d[a]),
            format!(
                "internal fan {} must be a nonempty block of grid cells",
                self.name
            ),
        )?;
        require(
            (2..=super::fan::MAX_CURVE_POINTS).contains(&self.curve.len()),
            format!("internal fan {} needs 2..64 curve points", self.name),
        )?;
        require(
            self.curve.iter().all(|p| {
                p.volume_flow_m3_per_s.is_finite()
                    && p.pressure_rise_pa.is_finite()
                    && p.volume_flow_m3_per_s >= 0.0
                    && p.pressure_rise_pa >= 0.0
            }) && self.curve.windows(2).all(|w| {
                w[1].volume_flow_m3_per_s > w[0].volume_flow_m3_per_s && w[1].pressure_rise_pa < w[0].pressure_rise_pa
            }),
            format!(
                "internal fan {}: curve flow must strictly increase while the pressure rise strictly decreases, all finite and non-negative",
                self.name
            ),
        )?;
        require(
            self.relaxation.is_finite() && self.relaxation > 0.0 && self.relaxation <= 1.0,
            format!("internal fan {} relaxation must be in (0, 1]", self.name),
        )?;
        let overlaps = |min: [usize; 3], max: [usize; 3]| {
            (0..3).all(|a| self.min[a] < max[a] && min[a] < self.max[a])
        };
        require(
            !case.obstacles.iter().any(|o| overlaps(o.min, o.max)),
            format!(
                "internal fan {} overlaps a solid obstacle; its cells must be fluid",
                self.name
            ),
        )?;
        require(
            !case.porous_blocks.iter().any(|p| overlaps(p.min, p.max)),
            format!("internal fan {} overlaps a porous block", self.name),
        )?;
        Ok(())
    }
}

/// One controller sample.
#[derive(Clone, Copy, Debug, Serialize, PartialEq)]
pub struct InternalFanSample {
    pub step: usize,
    pub volume_flow_m3_per_s: f64,
    pub pressure_rise_pa: f64,
}

struct FanState {
    cells: Vec<usize>,
    /// Block length along the direction and cross-section area (metres, m²).
    length: f64,
    area: f64,
    pressure_rise: f64,
    flow: f64,
    within_curve: bool,
    history: Vec<InternalFanSample>,
}

pub(super) struct InternalFans {
    fans: Vec<InternalFan>,
    states: Vec<FanState>,
    density: f64,
}

impl InternalFans {
    pub fn new(case: &CollocatedFlowCase, mesh: &PolyMesh) -> Result<Option<Self>> {
        let fans = case.internal_fans.clone();
        if fans.is_empty() {
            return Ok(None);
        }
        require(
            fans.len() <= MAX_INTERNAL_FANS,
            "at most 16 internal fans are supported",
        )?;
        require(
            case.mesh.xy_shear == 0.0,
            "internal fans require zero xy_shear",
        )?;
        super::unique(fans.iter().map(|f| f.name.as_str()))?;
        for (i, fan) in fans.iter().enumerate() {
            fan.validate(case)?;
            for other in fans.iter().take(i) {
                require(
                    !(0..3).all(|a| fan.min[a] < other.max[a] && other.min[a] < fan.max[a]),
                    format!("internal fans {} and {} overlap", fan.name, other.name),
                )?;
            }
        }
        let mut cells = vec![Vec::new(); fans.len()];
        for (ci, cell) in mesh.cells().iter().enumerate() {
            let ijk = super::porous::background_index(case, cell.center).ok_or_else(|| {
                SolverError::InvalidInput("cell outside the background grid".into())
            })?;
            if let Some(fi) = fans
                .iter()
                .position(|f| (0..3).all(|a| f.min[a] <= ijk[a] && ijk[a] < f.max[a]))
            {
                if let Some(list) = cells.get_mut(fi) {
                    list.push(ci);
                }
            }
        }
        let extent = |a: usize, lo: usize, hi: usize| {
            let nodes = case.mesh.nodes(a);
            nodes.get(hi).copied().unwrap_or(0.0) - nodes.get(lo).copied().unwrap_or(0.0)
        };
        let mut states = Vec::new();
        for (fan, list) in fans.iter().zip(cells) {
            let expected: usize = (0..3).map(|a| fan.max[a] - fan.min[a]).product();
            require(
                list.len() == expected,
                format!(
                    "internal fan {} covers solid cells (a surface body); its cells must be fluid",
                    fan.name
                ),
            )?;
            let axis = fan.direction.axis();
            let length = extent(axis, fan.min[axis], fan.max[axis]);
            let area: f64 = (0..3)
                .filter(|a| *a != axis)
                .map(|a| extent(a, fan.min[a], fan.max[a]))
                .product();
            // Start half-way along the curve, as a fan spinning up would pass through it.
            let mid = 0.5
                * (fan.curve.first().map_or(0.0, |p| p.volume_flow_m3_per_s)
                    + fan.curve.last().map_or(0.0, |p| p.volume_flow_m3_per_s));
            let (pressure_rise, _) = curve_pressure(&fan.curve, mid);
            states.push(FanState {
                cells: list,
                length,
                area,
                pressure_rise,
                flow: 0.0,
                within_curve: false,
                history: Vec::new(),
            });
        }
        Ok(Some(Self {
            fans,
            states,
            density: case.config.density,
        }))
    }

    /// Add each fan's body force (per unit mass) to `acceleration`.
    pub fn add_to(&self, acceleration: &mut [Vec3]) {
        for (fan, state) in self.fans.iter().zip(&self.states) {
            let a = fan.direction.unit() * (state.pressure_rise / (self.density * state.length));
            for &c in &state.cells {
                if let Some(slot) = acceleration.get_mut(c) {
                    *slot = *slot + a;
                }
            }
        }
    }

    /// Measure each fan's flow on the solved face fluxes and relax its pressure rise to the
    /// curve. Each cell contributes the mean of the fluxes along n through its two faces
    /// across the axis, weighted by its length along the axis.
    pub fn update(&mut self, step: usize, mesh: &PolyMesh, face_flux: &[f64]) {
        for (fan, state) in self.fans.iter().zip(self.states.iter_mut()) {
            let n = fan.direction.unit();
            let mut weighted = 0.0;
            for &c in &state.cells {
                let (Some(cell), Some(faces)) = (mesh.cells().get(c), mesh.cell_faces().get(c))
                else {
                    continue;
                };
                let mut sum = 0.0;
                let mut across = 0.0;
                for &fi in faces {
                    let (Some(face), Some(phi)) = (mesh.faces().get(fi), face_flux.get(fi)) else {
                        continue;
                    };
                    let along = face.area.dot(n);
                    if along.abs() > 0.5 * face.area.norm() {
                        sum += 0.5 * phi * along.signum();
                        across = along.abs();
                    }
                }
                if across > 0.0 {
                    weighted += sum * cell.volume / across;
                }
            }
            state.flow = weighted / state.length;
            let (target, within) = curve_pressure(&fan.curve, state.flow);
            state.within_curve = within;
            state.pressure_rise += fan.relaxation * (target - state.pressure_rise);
            if step % 10 == 0 || state.history.is_empty() {
                state.history.push(InternalFanSample {
                    step,
                    volume_flow_m3_per_s: state.flow,
                    pressure_rise_pa: state.pressure_rise,
                });
            }
        }
    }

    pub fn report(&self) -> serde_json::Value {
        let fans: Vec<_> = self
            .fans
            .iter()
            .zip(&self.states)
            .map(|(fan, state)| {
                let (curve_rise, _) = curve_pressure(&fan.curve, state.flow);
                serde_json::json!({
                    "name": fan.name,
                    "direction": fan.direction,
                    "cells": state.cells.len(),
                    "length_m": state.length,
                    "area_m2": state.area,
                    "volume_flow_m3_per_s": state.flow,
                    "mean_speed_m_per_s": state.flow / state.area,
                    "pressure_rise_pa": state.pressure_rise,
                    "curve_pressure_rise_pa": curve_rise,
                    "within_curve_range": state.within_curve,
                    "curve": fan.curve,
                    "samples": state.history,
                })
            })
            .collect();
        serde_json::json!({
            "schema": "cfd-solvers-rs/internal-fans/1",
            "model": "actuator_block_pressure_rise_from_curve",
            "fans": fans,
            "extrapolation": "refused_clamped_to_curve_range",
            "physical_validation": "not_established"
        })
    }

    pub fn write_report(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        std::fs::write(
            out.join("internal_fans.json"),
            serde_json::to_vec_pretty(&self.report())?,
        )?;
        m.files.push("internal_fans.json".into());
        for (fan, state) in self.fans.iter().zip(&self.states) {
            m.metrics.insert(
                format!("internal_fan_{}_volume_flow_m3_per_s", fan.name),
                state.flow,
            );
            m.metrics.insert(
                format!("internal_fan_{}_pressure_rise_pa", fan.name),
                state.pressure_rise,
            );
        }
        Ok(())
    }
}
