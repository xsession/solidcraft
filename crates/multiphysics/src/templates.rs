//! Case templates from the MIT cfd-solvers-rs 0.2.0 source release.
//! The companion performs full schema and physical validation at execution.
use crate::{Error, Result, Template, parse_json, profile_for_solver, validate_case};

const SOURCES: &[(&str, &str, &str)] = &[
    ("acoustic-kspace2d", "Acoustics · 2D k-space", include_str!("../data/templates/acoustic-kspace2d.json")),
    ("chemistry-batch", "Chemistry · isothermal batch", include_str!("../data/templates/chemistry-batch.json")),
    ("diffusion2d", "Diffusion · structured 2D", include_str!("../data/templates/diffusion2d.json")),
    ("flow-cavity-pimple", "Flow · transient cavity (PIMPLE)", include_str!("../data/templates/flow-cavity-pimple.json")),
    ("flow-cavity-simple", "Flow · steady cavity (SIMPLE)", include_str!("../data/templates/flow-cavity-simple.json")),
    ("flow-channel", "Flow · channel (PISO)", include_str!("../data/templates/flow-channel.json")),
    ("heat2d", "Heat · transient FEM 2D", include_str!("../data/templates/heat2d.json")),
    ("lbm-taylor-green", "Flow · D2Q9 Taylor–Green", include_str!("../data/templates/lbm-taylor-green.json")),
    ("maxwell1d", "EM · Yee Maxwell 1D", include_str!("../data/templates/maxwell1d.json")),
    ("optics-lens", "Optics · sequential lens", include_str!("../data/templates/optics-lens.json")),
    ("poisson1d", "Scalar · Poisson 1D", include_str!("../data/templates/poisson1d.json")),
    ("poisson2d", "Scalar · Poisson FEM 2D", include_str!("../data/templates/poisson2d.json")),
    ("poisson3d", "Scalar · Poisson FEM 3D", include_str!("../data/templates/poisson3d.json")),
    ("poly-scalar", "Transport · scalar FVM 3D", include_str!("../data/templates/poly-scalar.json")),
    ("rans-kepsilon", "Turbulence · k–epsilon transport", include_str!("../data/templates/rans-kepsilon.json")),
    ("rans-reynolds-stress", "Turbulence · Reynolds-stress transport", include_str!("../data/templates/rans-reynolds-stress.json")),
    ("rans-sa", "Turbulence · Spalart–Allmaras transport", include_str!("../data/templates/rans-sa.json")),
    ("rans-sst-des", "Turbulence · SST/DES transport", include_str!("../data/templates/rans-sst-des.json")),
    ("rf-tem", "RF · TEM two-port sweep", include_str!("../data/templates/rf-tem.json")),
    ("scalar-mixed", "Scalar · mixed-boundary FEM 2D", include_str!("../data/templates/scalar-mixed.json")),
    ("sod-weno", "Flow · Sod shock tube (WENO)", include_str!("../data/templates/sod-weno.json")),
    ("solid-heat", "Heat · solid conduction FVM 3D", include_str!("../data/templates/solid-heat.json")),
    ("thermal-radiation", "Heat · prescribed flow and radiation", include_str!("../data/templates/thermal-radiation.json")),
];

pub fn templates() -> Result<Vec<Template>> {
    SOURCES
        .iter()
        .map(|(id, name, text)| {
            let case = parse_json(text)?;
            validate_case(&case)?;
            let solver =
                case.get("solver").and_then(serde_json::Value::as_str).ok_or_else(|| Error::Invalid(format!("template {id} has no solver")))?;
            let profile = profile_for_solver(solver).ok_or_else(|| Error::Invalid(format!("template {id} has no supported profile")))?;
            Ok(Template { id: (*id).into(), name: (*name).into(), profile: profile.into(), case })
        })
        .collect()
}
