//! Bounded case adapters for the domains ported after 0.6.0: electrical and electronics
//! (`electrical`), solid and multibody mechanics (`mechanics`) and particle, kinetic,
//! multiphase, phase-field and coupled physics (`physics`). Every adapter takes a small,
//! typed SI input, runs one ported method and writes CSV/JSON (and VTU where the method has
//! a spatial field). Inputs are deliberately narrow; the library API is broader.
#![allow(clippy::needless_range_loop)]
use super::RunManifest;
use crate::{io, Result};
use serde::Serialize;
use std::path::Path;

pub mod electrical;
pub mod mechanics;
pub mod physics;

pub(crate) fn csv(
    out: &Path,
    m: &mut RunManifest,
    file: &str,
    columns: &[(&str, &[f64])],
) -> Result<()> {
    io::write_csv(out.join(file), columns)?;
    m.files.push(file.into());
    Ok(())
}

pub(crate) fn json(
    out: &Path,
    m: &mut RunManifest,
    file: &str,
    value: &impl Serialize,
) -> Result<()> {
    std::fs::write(out.join(file), serde_json::to_vec_pretty(value)?)?;
    m.files.push(file.into());
    Ok(())
}

pub(crate) fn metric(m: &mut RunManifest, name: &str, value: f64) {
    if value.is_finite() {
        m.metrics.insert(name.into(), value);
    }
}

/// Linearly spaced points including both ends (`n >= 2`) or the single `start`.
pub(crate) fn linspace(start: f64, stop: f64, n: usize) -> Vec<f64> {
    if n < 2 {
        return vec![start];
    }
    (0..n)
        .map(|i| start + (stop - start) * i as f64 / (n - 1) as f64)
        .collect()
}

pub(crate) fn logspace(start: f64, stop: f64, n: usize) -> Vec<f64> {
    if n < 2 || start <= 0.0 || stop <= 0.0 {
        return vec![start];
    }
    let r = (stop / start).powf(1.0 / (n - 1) as f64);
    (0..n).map(|i| start * r.powi(i as i32)).collect()
}
