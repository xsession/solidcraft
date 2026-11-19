//! Rotating regions for `collocated_flow3d` (`rotating_regions_v1`): the multiple reference
//! frame (MRF) model in the absolute-velocity formulation, as used for fans, impellers and
//! rotors in steady flow-analysis add-ins.
//!
//! In the cells of a region turning at Ω about an axis through `axis_point`:
//! - momentum gains the frame source −Ω×U (per unit mass, explicit, refreshed every step);
//! - momentum is convected by the flux relative to the frame, φ − (Ω×r)·S, on faces with both
//!   cells in the region and on the region's rotating walls;
//! - no-slip walls whose cell is in the region (obstacles, surface bodies and stationary
//!   exterior walls) move with the frame at Ω×r, so a rotor inside the region spins.
//!
//! Solid-body rotation of a closed box that lies wholly in a region is an exact steady
//! solution (U = Ω×r, p = ½Ω²r⊥²), which `tests/rotating.rs` checks.
use crate::{
    core::Vec3,
    fvm::{FaceProfiles, PolyMesh, VelocityBoundaryType},
    require, Result, SolverError,
};
use serde::{Deserialize, Serialize};

/// One rotating region: a block of background-grid cells turning about an axis.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RotatingRegion {
    pub name: String,
    /// Half-open background-grid cell range, as for obstacles.
    pub min: [usize; 3],
    pub max: [usize; 3],
    /// A point on the rotation axis (metres).
    pub axis_point: Vec3,
    /// The axis direction (normalised here; positive rotation by the right-hand rule).
    pub axis: Vec3,
    /// Speed of rotation; negative turns the other way.
    pub rpm: f64,
}

impl RotatingRegion {
    pub fn omega(&self) -> Vec3 {
        let n = self.axis.norm();
        self.axis * (self.rpm * std::f64::consts::TAU / 60.0 / n)
    }
    /// Frame velocity Ω×(x − axis point) at `x`.
    pub fn frame_velocity(&self, x: Vec3) -> Vec3 {
        self.omega().cross(x - self.axis_point)
    }
}

/// The regions resolved on the built mesh.
pub(super) struct RotatingFrames {
    regions: Vec<RotatingRegion>,
    /// Region of each mesh cell, if any.
    cell_region: Vec<Option<usize>>,
}

impl RotatingFrames {
    pub(super) fn new(case: &super::CollocatedFlowCase, mesh: &PolyMesh) -> Result<Self> {
        let regions = case.rotating_regions.clone();
        super::unique(regions.iter().map(|r| r.name.as_str()))?;
        let dims = case.mesh.dimensions;
        for r in &regions {
            require(
                (0..3).all(|a| r.min[a] < r.max[a] && r.max[a] <= dims[a]),
                format!("rotating region {} needs a non-empty cell range inside the grid", r.name),
            )?;
            require(r.axis.is_finite() && r.axis.norm() > 0.0, format!("rotating region {} needs a non-zero axis", r.name))?;
            require(r.axis_point.is_finite() && r.rpm.is_finite(), format!("rotating region {} needs a finite axis point and speed", r.name))?;
            // The frame source is explicit: keep the rotation per step small.
            let turn = r.omega().norm() * case.config.dt;
            require(
                turn <= 0.2,
                format!(
                    "rotating region {} turns {turn:.3} rad per time step; lower dt to at most {:.3e} s",
                    r.name,
                    0.2 / r.omega().norm().max(1e-300)
                ),
            )?;
        }
        for (i, a) in regions.iter().enumerate() {
            for b in &regions[i + 1..] {
                let overlap = (0..3).all(|k| a.min[k] < b.max[k] && b.min[k] < a.max[k]);
                require(!overlap, format!("rotating regions {} and {} overlap", a.name, b.name))?;
            }
        }
        let cell_region = mesh
            .cells()
            .iter()
            .map(|cell| {
                let index = [0, 1, 2].map(|a| case.mesh.cell_of(a, cell.center.component(a)));
                regions.iter().position(|r| (0..3).all(|a| index[a].is_some_and(|i| (r.min[a]..r.max[a]).contains(&i))))
            })
            .collect();
        Ok(Self { regions, cell_region })
    }

    pub(super) fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }

    /// Walls (fixed-value, zero-velocity patches, not transferred) that rotate: faces whose cell
    /// lies in a region. Returns the faces and their frame velocity, and fills `profiles`.
    pub(super) fn rotate_walls(&self, case: &super::CollocatedFlowCase, mesh: &PolyMesh, profiles: &mut FaceProfiles) -> Result<Vec<f64>> {
        let mut offset = vec![0.0; mesh.face_count()];
        if self.is_empty() {
            return Ok(offset);
        }
        if profiles.velocity.is_empty() {
            profiles.velocity = vec![None; mesh.face_count()];
        }
        let wall = |patch: &str| -> bool {
            case.solid_names().contains(&patch)
                || case.velocity_boundaries.iter().any(|b| b.patch == patch && b.kind == VelocityBoundaryType::FixedValue && b.value == Vec3::ZERO)
        };
        for (fi, face) in mesh.faces().iter().enumerate() {
            let Some(region) = self.cell_region.get(face.owner).copied().flatten() else { continue };
            let r = &self.regions[region];
            match (face.neighbour, face.patch) {
                (Some(n), _) => {
                    // Interior faces: relative convection only where both cells turn together.
                    if self.cell_region.get(n).copied().flatten() == Some(region) {
                        offset[fi] = r.frame_velocity(face.center).dot(face.area);
                    }
                }
                (None, Some(p)) => {
                    let name = mesh.patches()[p].name.as_str();
                    if !wall(name) {
                        continue;
                    }
                    require(
                        profiles.velocity.get(fi).copied().flatten().is_none(),
                        format!("patch {name} is both transferred and inside rotating region {}", r.name),
                    )?;
                    let v = r.frame_velocity(face.center);
                    profiles.velocity[fi] = Some(v);
                    offset[fi] = v.dot(face.area);
                }
                (None, None) => return Err(SolverError::InvalidInput("boundary face without patch".into())),
            }
        }
        Ok(offset)
    }

    /// Body acceleration: `base` plus the frame source −Ω×U in region cells.
    pub(super) fn acceleration(&self, base: &[Vec3], velocity: &[Vec3]) -> Vec<Vec3> {
        base.iter()
            .zip(velocity)
            .zip(&self.cell_region)
            .map(|((a, u), region)| match region {
                Some(i) => *a - self.regions[*i].omega().cross(*u),
                None => *a,
            })
            .collect()
    }

    /// Report: per region, the cells, the speed, and the fluid's mean angular velocity
    /// about the axis relative to the frame's (1 = turning with the frame).
    pub(super) fn report(&self, mesh: &PolyMesh, velocity: &[Vec3]) -> serde_json::Value {
        let rows: Vec<serde_json::Value> = self
            .regions
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let (mut cells, mut fluid, mut frame) = (0usize, 0.0, 0.0);
                let axis = r.axis * (1.0 / r.axis.norm());
                for (c, cell) in mesh.cells().iter().enumerate() {
                    if self.cell_region[c] != Some(i) {
                        continue;
                    }
                    cells += 1;
                    let d = cell.center - r.axis_point;
                    let radial = d - axis * d.dot(axis);
                    let tangent = axis.cross(radial);
                    fluid += velocity[c].dot(tangent) * cell.volume;
                    frame += r.frame_velocity(cell.center).dot(tangent) * cell.volume;
                }
                serde_json::json!({"name": r.name, "cells": cells, "rpm": r.rpm,
                    "omega_rad_per_s": r.omega().norm(),
                    "swirl_ratio": if frame.abs() > 0.0 { fluid / frame } else { 0.0 }})
            })
            .collect();
        serde_json::json!({"schema": "cfd-solvers-rs/rotating-regions/1", "model": "mrf_absolute_velocity", "regions": rows})
    }
}
