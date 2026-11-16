//! Offline self-verification suite (`cfd-solve verify`).
//!
//! Every check runs locally from embedded inputs, with no network access and no external
//! data. Each states its criterion next to the measured value. Passing establishes code
//! verification (discretization order, exact discrete conservation, determinism and
//! agreement with closed-form solutions) for these small cases only, not physical validation.
use super::{conjugate, Case, PolyMeshSpec};
use crate::{core::Vec3, Result, SolverError};
use serde_json::{json, Value};
use std::path::Path;

fn check(name: &str, criterion: &str, measured: Value, passed: bool) -> Value {
    json!({"name": name, "criterion": criterion, "measured": measured, "passed": passed})
}

fn wall(patch: &str, temperature_k: f64) -> conjugate::ThermalPatch {
    conjugate::ThermalPatch {
        patch: patch.into(),
        condition: conjugate::ThermalCondition::Temperature { temperature_k },
    }
}

fn energy_order() -> Result<Value> {
    let mut errors = vec![];
    for n in [8usize, 16, 32] {
        let spec = PolyMeshSpec {
            dimensions: [n, n, 1],
            lengths: Vec3::new(1.0, 1.0, 0.1),
            xy_shear: 0.0,
        };
        let (mesh, _) = spec.build()?;
        let pi = std::f64::consts::PI;
        let exact: Vec<f64> = mesh
            .cells()
            .iter()
            .map(|c| (pi * c.center.x).sin() * (pi * c.center.y).sin())
            .collect();
        let source: Vec<f64> = mesh
            .cells()
            .iter()
            .zip(&exact)
            .map(|(c, e)| 2.0 * pi * pi * e * c.volume)
            .collect();
        let patches: Vec<_> = ["left", "right", "bottom", "top"]
            .iter()
            .map(|p| wall(p, 300.0))
            .collect();
        let s = conjugate::solve_energy(
            spec,
            &vec![[1.0; 3]; n * n],
            &vec![0.0; mesh.face_count()],
            &source,
            &patches,
            300.0,
            1.0,
        )?;
        let volume: f64 = mesh.cells().iter().map(|c| c.volume).sum();
        let sum: f64 = mesh
            .cells()
            .iter()
            .zip(&s.temperature_k)
            .zip(&exact)
            .map(|((c, t), e)| (t - 300.0 - e).powi(2) * c.volume)
            .sum();
        errors.push((sum / volume).sqrt());
    }
    let orders: Vec<f64> = errors.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    let passed = orders.iter().all(|o| (o - 2.0).abs() < 0.1);
    Ok(check(
        "energy_operator_order",
        "observed L2 order 2 +/- 0.1 (8, 16, 32 cells)",
        json!({"l2_errors": errors, "orders": orders}),
        passed,
    ))
}

fn layered_wall() -> Result<Value> {
    let spec = PolyMeshSpec {
        dimensions: [6, 1, 1],
        lengths: Vec3::new(0.06, 0.01, 0.01),
        xy_shear: 0.0,
    };
    let (mesh, _) = spec.build()?;
    let k = [1.0, 1.0, 0.1, 0.1, 10.0, 10.0];
    let conductivity: Vec<[f64; 3]> = k.iter().map(|k| [*k; 3]).collect();
    let s = conjugate::solve_energy(
        spec,
        &conductivity,
        &vec![0.0; mesh.face_count()],
        &[0.0; 6],
        &[wall("left", 400.0), wall("right", 300.0)],
        300.0,
        1.0,
    )?;
    let exact = 100.0 / k.iter().map(|k| 0.01 / k).sum::<f64>() * 1e-4;
    let q = s
        .conduction_out_w
        .iter()
        .find(|(p, _)| p == "right")
        .map_or(f64::NAN, |(_, q)| *q);
    let relative = (q - exact).abs() / exact;
    Ok(check(
        "layered_wall_conduction",
        "heat flow within 1e-9 of the series-resistance solution",
        json!({"w": q, "exact_w": exact, "relative_error": relative}),
        relative < 1e-9,
    ))
}

fn run_case(value: Value, dir: &Path) -> Result<super::RunManifest> {
    let case: Case = serde_json::from_value(value)?;
    case.run(dir)
}

fn channel(extra: Value) -> Value {
    let mut base = json!({
        "solver": "collocated_flow3d", "method": "piso",
        "mesh": {"dimensions": [12, 6, 1], "lengths": {"x": 0.12, "y": 0.02, "z": 0.01}},
        "config": {"dt": 0.01, "density": 1.2, "kinematic_viscosity": 1.5e-5, "pressure_tolerance": 1e-12, "pressure_iterations": 4000},
        "output": {"steps": 300, "output_interval": 300},
        "velocity_boundaries": [
            {"patch": "right", "kind": "zero_gradient"}, {"patch": "bottom", "kind": "fixed_value"},
            {"patch": "top", "kind": "fixed_value"}, {"patch": "front", "kind": "slip"}, {"patch": "back", "kind": "slip"}],
        "pressure_boundaries": [{"patch": "right", "kind": "fixed_value", "value": 0.0}],
        "fan": {"patch": "left", "reference_patch": "right", "pressure_convention": "static",
            "curve": [{"volume_flow_m3_per_s": 0.0, "pressure_rise_pa": 0.02}, {"volume_flow_m3_per_s": 2e-5, "pressure_rise_pa": 0.0}]}
    });
    if let (Some(base), Some(extra)) = (base.as_object_mut(), extra.as_object()) {
        for (k, v) in extra {
            base.insert(k.clone(), v.clone());
        }
    }
    base
}

fn read_json(path: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&std::fs::read(path)?)?)
}

fn fan_consistency(dir: &Path) -> Result<Value> {
    run_case(channel(json!({})), &dir.join("fan"))?;
    let point = read_json(&dir.join("fan/fan_operating_point.json"))?;
    let q = point["volume_flow_m3_per_s"].as_f64().unwrap_or(f64::NAN);
    let p = point["system_pressure_rise_pa"]
        .as_f64()
        .unwrap_or(f64::NAN);
    let mut fixed = channel(json!({}));
    fixed["fan"]["initial_volume_flow_m3_per_s"] = json!(q);
    fixed["fan"]["update_interval"] = json!(10_000);
    run_case(fixed, &dir.join("fan_fixed"))?;
    let independent = read_json(&dir.join("fan_fixed/fan_operating_point.json"))?
        ["system_pressure_rise_pa"]
        .as_f64()
        .unwrap_or(f64::NAN);
    let mismatch = point["relative_operating_point_mismatch"]
        .as_f64()
        .unwrap_or(f64::NAN);
    let system = (independent - p).abs() / p.abs();
    Ok(check(
        "fan_operating_point",
        "curve mismatch < 1% and an independent fixed-flow run of the same system within 1%",
        json!({"volume_flow_m3_per_s": q, "pressure_pa": p, "curve_mismatch": mismatch, "independent_pressure_pa": independent, "system_difference": system}),
        mismatch < 0.01 && system < 0.01,
    ))
}

fn conjugate_balance(dir: &Path) -> Result<Value> {
    let case = channel(json!({
        "obstacles": [{"name": "chip", "min": [5, 0, 0], "max": [7, 2, 1]}],
        "thermal": {"fluid": {"specific_heat_j_per_kg_k": 1006.0, "conductivity_w_per_m_k": 0.026}, "inflow_temperature_k": 298.15,
            "solids": [{"obstacle": "chip", "conductivity_w_per_m_k": [150.0, 150.0, 150.0], "power_w": 0.01}],
            "boundaries": [{"patch": "bottom", "condition": {"kind": "convection", "h_w_per_m2_k": 20.0, "ambient_k": 298.15}}]}
    }));
    run_case(case, &dir.join("cht"))?;
    let r = read_json(&dir.join("cht/thermal_report.json"))?;
    let imbalance = r["relative_energy_imbalance"].as_f64().unwrap_or(f64::NAN);
    Ok(check(
        "conjugate_energy_balance",
        "relative global energy imbalance < 1e-9",
        json!({"relative_imbalance": imbalance}),
        imbalance < 1e-9,
    ))
}

fn porous_drop(dir: &Path) -> Result<Value> {
    let mut errors = vec![];
    for n in [20usize, 40] {
        let mut case = channel(json!({
            "mesh": {"dimensions": [n, 2, 2], "lengths": {"x": 0.2, "y": 0.02, "z": 0.02}},
            "output": {"steps": 200, "output_interval": 200},
            "porous_blocks": [{"name": "vent", "min": [2 * n / 5, 0, 0], "max": [3 * n / 5, 2, 2],
                "viscous_resistance_per_m2": 1.0e6, "inertial_resistance_per_m": 0.0}]
        }));
        case["config"]["dt"] = json!(0.02);
        case["velocity_boundaries"] = json!([
            {"patch": "right", "kind": "zero_gradient"}, {"patch": "bottom", "kind": "slip"},
            {"patch": "top", "kind": "slip"}, {"patch": "front", "kind": "slip"}, {"patch": "back", "kind": "slip"}]);
        case["fan"]["curve"] = json!([{"volume_flow_m3_per_s": 0.0, "pressure_rise_pa": 1.0}, {"volume_flow_m3_per_s": 1e-3, "pressure_rise_pa": 0.0}]);
        case["fan"]["initial_volume_flow_m3_per_s"] = json!(2e-5);
        case["fan"]["update_interval"] = json!(10_000);
        let folder = dir.join(format!("porous{n}"));
        run_case(case, &folder)?;
        let rise = read_json(&folder.join("fan_operating_point.json"))?["system_pressure_rise_pa"]
            .as_f64()
            .unwrap_or(f64::NAN);
        let exact = 1.2 * 1.5e-5 * 1.0e6 * 0.05 * 0.04;
        errors.push((rise - exact).abs() / exact);
    }
    let passed = errors
        .first()
        .zip(errors.get(1))
        .is_some_and(|(a, b)| *b < 0.03 && *b < 0.6 * a);
    Ok(check(
        "porous_darcy_plug",
        "Darcy plug drop error < 3% at 8 plug cells and shrinking under refinement",
        json!({"relative_errors": errors}),
        passed,
    ))
}

/// de Vahl Davis (1983) differentially heated square cavity, Ra = 1e3, Pr = 0.71.
/// The benchmark average Nusselt number is 1.118.
fn natural_convection(dir: &Path) -> Result<Value> {
    let (alpha, nu) = (0.01, 0.0071);
    let case = json!({
        "solver": "collocated_flow3d", "method": "piso",
        "mesh": {"dimensions": [16, 16, 1], "lengths": {"x": 1.0, "y": 1.0, "z": 0.1}},
        "config": {"dt": 0.04, "density": 1.0, "kinematic_viscosity": nu, "pressure_tolerance": 1e-12,
            "pressure_iterations": 4000, "pressure_preconditioner": "jacobi"},
        "output": {"steps": 2000, "output_interval": 2000},
        "velocity_boundaries": [
            {"patch": "left", "kind": "fixed_value"}, {"patch": "right", "kind": "fixed_value"},
            {"patch": "bottom", "kind": "fixed_value"}, {"patch": "top", "kind": "fixed_value"},
            {"patch": "front", "kind": "slip"}, {"patch": "back", "kind": "slip"}],
        "thermal": {"fluid": {"specific_heat_j_per_kg_k": 1.0, "conductivity_w_per_m_k": alpha}, "inflow_temperature_k": 300.0, "solids": [],
            "boundaries": [
                {"patch": "left", "condition": {"kind": "temperature", "temperature_k": 300.5}},
                {"patch": "right", "condition": {"kind": "temperature", "temperature_k": 299.5}}],
            "buoyancy": {"gravity_m_per_s2": [0.0, -1.0, 0.0], "expansion_per_k": 1.0e3 * nu * alpha, "reference_temperature_k": 300.0}}
    });
    run_case(case, &dir.join("cavity"))?;
    let r = read_json(&dir.join("cavity/thermal_report.json"))?;
    let mut hot = f64::NAN;
    for p in r["conducted_heat_out_w"].as_array().into_iter().flatten() {
        if p["patch"] == "left" {
            hot = -p["w"].as_f64().unwrap_or(f64::NAN);
        }
    }
    let nusselt = hot / (alpha * 1.0 * 0.1);
    let error = (nusselt - 1.118) / 1.118;
    Ok(check(
        "natural_convection_de_vahl_davis_ra1e3",
        "average hot-wall Nusselt number within 5% of the 1.118 benchmark on a 16x16 grid",
        json!({"nusselt": nusselt, "benchmark": 1.118, "relative_error": error}),
        error.abs() < 0.05,
    ))
}

fn determinism(dir: &Path) -> Result<Value> {
    let case = channel(json!({"output": {"steps": 40, "output_interval": 20}}));
    run_case(case.clone(), &dir.join("repeat_a"))?;
    run_case(case, &dir.join("repeat_b"))?;
    let mut compared = 0;
    let mut identical = true;
    for name in [
        "flow_000000.vtu",
        "flow_000020.vtu",
        "flow_000040.vtu",
        "fan_operating_point.csv",
    ] {
        let a = std::fs::read(dir.join("repeat_a").join(name))?;
        let b = std::fs::read(dir.join("repeat_b").join(name))?;
        identical &= a == b;
        compared += 1;
    }
    Ok(check(
        "bitwise_determinism",
        "two runs of the same input write byte-identical fields",
        json!({"files_compared": compared}),
        identical,
    ))
}

/// Run every check into `out` (which must not exist yet) and write `verification.json`.
pub fn run_suite(out: &Path) -> Result<Value> {
    if out.exists() {
        return Err(SolverError::InvalidInput(
            "verification output folder already exists; choose a new folder".into(),
        ));
    }
    std::fs::create_dir_all(out)?;
    let work = out.join("cases");
    let checks = vec![
        energy_order()?,
        layered_wall()?,
        fan_consistency(&work)?,
        porous_drop(&work)?,
        conjugate_balance(&work)?,
        natural_convection(&work)?,
        determinism(&work)?,
    ];
    let passed = checks.iter().all(|c| c["passed"] == json!(true));
    let report = json!({
        "schema": "cfd-solvers-rs/verification/1",
        "rust_port_version": env!("CARGO_PKG_VERSION"),
        "network": "not used",
        "scope": "code verification of small embedded cases; not physical validation",
        "passed": passed,
        "checks": checks
    });
    std::fs::write(
        out.join("verification.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    Ok(report)
}
