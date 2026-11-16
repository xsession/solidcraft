//! Admission checks for the explicitly discretized primitive airflow domain.
//! This does not classify CAD bodies or establish mesh connectivity.
use crate::{Error, MAX_CELLS, Result};
use serde::Deserialize;
use serde_json::Value;
use std::collections::BTreeSet;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GridObstacle {
    name: String,
    min: [u64; 3],
    max: [u64; 3],
}

pub(crate) fn validate(case: &Value) -> Result<()> {
    let Some(value) = case.get("obstacles") else { return Ok(()) };
    if case.get("solver").and_then(Value::as_str) != Some("collocated_flow3d") {
        return Err(Error::Invalid("grid obstacles require collocated_flow3d".into()));
    }
    let array = value.as_array().ok_or_else(|| Error::Invalid("obstacles must be an array".into()))?;
    if array.len() > 128 {
        return Err(Error::Limit("grid obstacles exceed 128".into()));
    }
    if array.is_empty() {
        return Ok(());
    }
    let mesh = case.get("mesh").ok_or_else(|| Error::Invalid("grid obstacles require an explicit mesh".into()))?;
    let dimensions: [u64; 3] = serde_json::from_value(mesh.get("dimensions").cloned().unwrap_or(Value::Null))
        .map_err(|_| Error::Invalid("obstacle mesh dimensions must be three positive integers".into()))?;
    let cells = dimensions.iter().try_fold(1u64, |total, size| {
        if *size == 0 {
            return Err(Error::Invalid("obstacle mesh dimensions must be positive".into()));
        }
        total.checked_mul(*size).ok_or_else(|| Error::Limit("obstacle mesh size overflow".into()))
    })?;
    if cells > MAX_CELLS {
        return Err(Error::Limit(format!("obstacle mesh exceeds {MAX_CELLS} cells")));
    }
    for axis in ["x", "y", "z"] {
        if mesh.get("lengths").and_then(|lengths| lengths.get(axis)).and_then(Value::as_f64).is_none_or(|v| !v.is_finite() || v <= 0.0) {
            return Err(Error::Invalid("obstacle mesh lengths must be finite and positive".into()));
        }
    }
    if mesh.get("xy_shear").is_some_and(|value| value.as_f64() != Some(0.0)) {
        return Err(Error::Invalid("grid obstacles require zero xy_shear".into()));
    }
    let mut obstacles: Vec<GridObstacle> = Vec::with_capacity(array.len());
    let mut names: BTreeSet<String> = ["left", "right", "bottom", "top", "front", "back"].into_iter().map(String::from).collect();
    for value in array {
        let obstacle: GridObstacle =
            serde_json::from_value(value.clone()).map_err(|error| Error::Invalid(format!("invalid grid obstacle: {error}")))?;
        if obstacle.name.is_empty()
            || obstacle.name.len() > 80
            || !obstacle.name.bytes().next().is_some_and(|byte| byte.is_ascii_alphanumeric())
            || !obstacle.name.bytes().all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-'))
        {
            return Err(Error::Invalid(
                "obstacle names must use 1–80 ASCII letters, digits, underscores or hyphens and start with a letter or digit".into(),
            ));
        }
        if !names.insert(obstacle.name.clone()) {
            return Err(Error::Invalid("duplicate or reserved obstacle name".into()));
        }
        if obstacle.min.iter().zip(obstacle.max).zip(dimensions).any(|((min, max), size)| *min >= max || max > size) {
            return Err(Error::Invalid("obstacle bounds must satisfy 0 <= min < max <= dimensions".into()));
        }
        if obstacles.iter().any(|other| {
            (0..3).all(|axis| {
                obstacle.min.get(axis).zip(other.max.get(axis)).is_some_and(|(min, max)| min < max)
                    && other.min.get(axis).zip(obstacle.max.get(axis)).is_some_and(|(min, max)| min < max)
            })
        }) {
            return Err(Error::Invalid("grid obstacles overlap".into()));
        }
        obstacles.push(obstacle);
    }
    Ok(())
}
