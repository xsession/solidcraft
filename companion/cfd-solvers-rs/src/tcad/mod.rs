//! Semiconductor TCAD: 1-D drift-diffusion devices, MOS capacitors and compact MOSFET/BJT/IGBT,
//! Gummel-Poon, MOS level 2/3, VDMOS and hysteretic-core baselines.
mod devices;
mod mos;
mod semiconductor;
pub use devices::*;
pub use mos::*;
pub use semiconductor::*;

use crate::Result;
use std::sync::Arc;

/// Static compact-device evaluation; currents are positive flowing out of each terminal and the
/// Jacobian is row-major `dI[row]/dV[col]` (upstream `cfd::circuit::StaticDeviceEvaluation`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StaticDeviceEvaluation {
    pub terminal_current: Vec<f64>,
    pub jacobian: Vec<f64>,
}

/// Upstream `cfd::circuit::StaticDeviceEvaluator`.
pub type StaticDeviceEvaluator =
    Arc<dyn Fn(&[f64]) -> Result<StaticDeviceEvaluation> + Send + Sync>;

/// Centered finite-difference Jacobian over `terminals` voltages; `step(v)` picks the step.
pub(crate) fn fd_evaluation(
    v: &[f64],
    terminals: usize,
    currents: &dyn Fn(&[f64]) -> Result<Vec<f64>>,
    step: &dyn Fn(f64) -> f64,
) -> Result<StaticDeviceEvaluation> {
    let current = currents(v)?;
    let mut jac = vec![0.0; terminals * terminals];
    for col in 0..terminals {
        let h = step(v[col]);
        let (mut vp, mut vm) = (v.to_vec(), v.to_vec());
        vp[col] += h;
        vm[col] -= h;
        let (ip, im) = (currents(&vp)?, currents(&vm)?);
        for row in 0..terminals {
            jac[row * terminals + col] = (ip[row] - im[row]) / (2.0 * h);
        }
    }
    Ok(StaticDeviceEvaluation {
        terminal_current: current,
        jacobian: jac,
    })
}
