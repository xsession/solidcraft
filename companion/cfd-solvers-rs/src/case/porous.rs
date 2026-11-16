//! Grid-aligned isotropic porous resistance zones (`porous_blocks_v1`).
//!
//! Each zone keeps its cells in the fluid mesh and adds a Darcy-Forchheimer momentum
//! sink per unit mass `a = -(nu * D + 0.5 * C * |u|) u`, with viscous coefficient `D`
//! (1/m^2, the reciprocal permeability) and inertial coefficient `C` (1/m). The velocity
//! is the superficial (volume-averaged) velocity. Porosity effects on storage and
//! convection, anisotropic coefficients and thermal effects are not modelled. Coefficients
//! must come from a characterized vent, filter or grille; open-area fraction alone is not
//! sufficient to derive them.
use super::{CollocatedFlowCase, RunManifest};
use crate::{core::Vec3, fvm::PolyMesh, require, Result, SolverError};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, path::Path};

pub const MAX_POROUS_BLOCKS: usize = 32;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PorousBlock {
    pub name: String,
    /// Half-open background-grid cell indices, as for obstacles.
    pub min: [usize; 3],
    pub max: [usize; 3],
    /// Darcy coefficient D = 1/K in 1/m^2.
    pub viscous_resistance_per_m2: f64,
    /// Forchheimer coefficient C in 1/m.
    pub inertial_resistance_per_m: f64,
}

pub(super) struct PorousZones {
    /// (block, retained cell indices)
    cells: Vec<Vec<usize>>,
    viscosity: f64,
}

fn background_index(mesh: &CollocatedFlowCase, center: Vec3) -> Option<[usize; 3]> {
    let d = mesh.mesh.dimensions;
    let l = mesh.mesh.lengths;
    let axis = |x: f64, len: f64, n: usize| -> Option<usize> {
        let i = (x / len * n as f64).floor();
        (i >= 0.0 && i < n as f64).then_some(i as usize)
    };
    Some([
        axis(center.x, l.x, d[0])?,
        axis(center.y, l.y, d[1])?,
        axis(center.z, l.z, d[2])?,
    ])
}

impl PorousZones {
    pub fn new(case: &CollocatedFlowCase, mesh: &PolyMesh) -> Result<Self> {
        let blocks = &case.porous_blocks;
        require(
            blocks.len() <= MAX_POROUS_BLOCKS,
            "at most 32 porous blocks are supported",
        )?;
        require(
            case.mesh.xy_shear == 0.0,
            "porous blocks require zero xy_shear",
        )?;
        let mut names: BTreeSet<&str> = ["left", "right", "bottom", "top", "front", "back"]
            .into_iter()
            .collect();
        names.extend(case.obstacles.iter().map(|o| o.name.as_str()));
        for (bi, block) in blocks.iter().enumerate() {
            let name = block.name.as_bytes();
            require(
                !name.is_empty()
                    && name.len() <= 80
                    && name.first().is_some_and(u8::is_ascii_alphanumeric)
                    && name.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'_' || *c == b'-')
                    && names.insert(&block.name),
                "porous block names must be unique ASCII identifiers distinct from patches and obstacles",
            )?;
            require(
                (0..3).all(|a| block.min[a] < block.max[a] && block.max[a] <= case.mesh.dimensions[a]),
                "porous block bounds must be nonempty half-open grid cell indices within dimensions",
            )?;
            let (d, c) = (
                block.viscous_resistance_per_m2,
                block.inertial_resistance_per_m,
            );
            require(
                d.is_finite() && c.is_finite() && d >= 0.0 && c >= 0.0 && (d > 0.0 || c > 0.0),
                "porous coefficients must be finite, non-negative and not both zero",
            )?;
            let overlaps = |min: [usize; 3], max: [usize; 3]| {
                (0..3).all(|a| block.min[a] < max[a] && min[a] < block.max[a])
            };
            require(
                !case.obstacles.iter().any(|o| overlaps(o.min, o.max)),
                "porous blocks must not overlap solid obstacles",
            )?;
            require(
                !blocks.iter().take(bi).any(|o| overlaps(o.min, o.max)),
                "porous blocks must not overlap each other",
            )?;
        }
        let mut cells = vec![vec![]; blocks.len()];
        for (ci, cell) in mesh.cells().iter().enumerate() {
            let Some(ijk) = background_index(case, cell.center) else {
                return Err(SolverError::InvalidInput(
                    "cell outside the background grid".into(),
                ));
            };
            if let Some(bi) = blocks
                .iter()
                .position(|b| (0..3).all(|a| b.min[a] <= ijk[a] && ijk[a] < b.max[a]))
            {
                if let Some(list) = cells.get_mut(bi) {
                    list.push(ci);
                }
            }
        }
        for (block, list) in blocks.iter().zip(&cells) {
            let expected: usize = (0..3).map(|a| block.max[a] - block.min[a]).product();
            require(
                list.len() == expected,
                format!("porous block {} does not map to its grid cells", block.name),
            )?;
        }
        Ok(Self {
            cells,
            viscosity: case.config.kinematic_viscosity,
        })
    }

    /// Linearized drag coefficients for the current velocity (Picard iteration).
    pub fn drag(&self, case: &CollocatedFlowCase, velocity: &[Vec3]) -> Result<Vec<f64>> {
        let mut beta = vec![0.0; velocity.len()];
        for (block, list) in case.porous_blocks.iter().zip(&self.cells) {
            for &ci in list {
                let speed = velocity
                    .get(ci)
                    .map(|u| u.norm())
                    .ok_or_else(|| SolverError::InvalidInput("porous cell".into()))?;
                if let Some(b) = beta.get_mut(ci) {
                    *b = self.viscosity * block.viscous_resistance_per_m2
                        + 0.5 * block.inertial_resistance_per_m * speed;
                }
            }
        }
        Ok(beta)
    }

    pub fn write_report(
        &self,
        case: &CollocatedFlowCase,
        mesh: &PolyMesh,
        velocity: &[Vec3],
        out: &Path,
        m: &mut RunManifest,
    ) -> Result<()> {
        let beta = self.drag(case, velocity)?;
        let mut zones = vec![];
        let mut total_power = 0.0;
        for (block, list) in case.porous_blocks.iter().zip(&self.cells) {
            let mut volume = 0.0;
            let mut speed_volume = 0.0;
            let mut power = 0.0;
            for &ci in list {
                let v = mesh.cells()[ci].volume;
                let u = velocity[ci];
                volume += v;
                speed_volume += u.norm() * v;
                power += case.config.density * beta[ci] * u.dot(u) * v;
            }
            total_power += power;
            zones.push(serde_json::json!({
                "name": block.name, "min": block.min, "max": block.max, "cells": list.len(), "volume_m3": volume,
                "viscous_resistance_per_m2": block.viscous_resistance_per_m2,
                "inertial_resistance_per_m": block.inertial_resistance_per_m,
                "mean_superficial_speed_m_per_s": if volume > 0.0 { speed_volume / volume } else { 0.0 },
                "dissipated_power_W": power
            }));
        }
        let report = serde_json::json!({
            "schema": "cfd-solvers-rs/porous-zones/1",
            "model": "isotropic_darcy_forchheimer_superficial_velocity",
            "sink_per_unit_mass": "-(nu*D + 0.5*C*|u|) u",
            "zones": zones,
            "total_dissipated_power_W": total_power,
            "omitted": ["porosity_storage_and_convection", "anisotropy", "thermal_effects"],
            "physical_validation": "not_established"
        });
        std::fs::write(
            out.join("porous_zones.json"),
            serde_json::to_vec_pretty(&report)?,
        )?;
        m.files.push("porous_zones.json".into());
        m.metrics
            .insert("porous_dissipated_power_W".into(), total_power);
        Ok(())
    }
}
