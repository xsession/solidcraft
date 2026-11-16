//! Momentum-exchange force/torque on solid links and Darcy-Forchheimer drag, ported from upstream
//! `cfd/solvers/lbm/diagnostics` (MIT, cfd_solvers contributors).
use super::{D2Q9Solver, CX, CY};
use crate::{require, Result};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ForceTorque2D {
    pub force_x: f64,
    pub force_y: f64,
    pub torque_z: f64,
}

/// Sum of `2 f_q c_q` over fluid-to-solid links (no wrap at the domain edge); torque is about
/// `(origin_x, origin_y)` measured from the fluid-cell centre `(x + 0.5, y + 0.5)`.
pub fn momentum_exchange_force(solver: &D2Q9Solver, origin_x: f64, origin_y: f64) -> ForceTorque2D {
    let c = solver.config();
    let (nx, ny) = (c.nx, c.ny);
    let solid = solver.solid_mask();
    let f = solver.populations();
    let mut out = ForceTorque2D::default();
    for y in 0..ny {
        for x in 0..nx {
            let n = y * nx + x;
            if solid[n] {
                continue;
            }
            for q in 1..9 {
                let xx = x as i64 + i64::from(CX[q]);
                let yy = y as i64 + i64::from(CY[q]);
                if xx < 0 || yy < 0 || xx >= nx as i64 || yy >= ny as i64 {
                    continue;
                }
                if !solid[yy as usize * nx + xx as usize] {
                    continue;
                }
                let fq = f64::from(f[n][q]);
                let fx = 2.0 * fq * f64::from(CX[q]);
                let fy = 2.0 * fq * f64::from(CY[q]);
                out.force_x += fx;
                out.force_y += fy;
                let rx = x as f64 + 0.5 - origin_x;
                let ry = y as f64 + 0.5 - origin_y;
                out.torque_z += rx * fy - ry * fx;
            }
        }
    }
    out
}

/// Porous body acceleration `-(nu/K + c_F |u| / sqrt(K)) u`.
pub fn darcy_forchheimer_acceleration(
    ux: f64,
    uy: f64,
    kinematic_viscosity: f64,
    permeability: f64,
    forchheimer_coefficient: f64,
) -> Result<[f64; 2]> {
    require(
        kinematic_viscosity >= 0.0 && permeability > 0.0 && forchheimer_coefficient >= 0.0,
        "invalid porous drag controls",
    )?;
    let speed = ux.hypot(uy);
    let factor = -(kinematic_viscosity / permeability
        + forchheimer_coefficient * speed / permeability.sqrt());
    Ok([factor * ux, factor * uy])
}
