//! Case templates from the MIT cfd-solvers-rs source release and original examples.
//! The companion performs full schema and physical validation at execution.
use crate::{Error, Result, Template, parse_json, profile_for_solver, validate_case};

const SOURCES: &[(&str, &str, &str)] = &[
    ("acoustic-kspace2d", "Acoustics · 2D k-space", include_str!("../data/templates/acoustic-kspace2d.json")),
    ("bar-modes", "Structures · axial bar modes", include_str!("../data/templates/bar-modes.json")),
    ("battery-discharge", "Battery · SPMe cell 1C discharge and rest", include_str!("../data/templates/battery-discharge.json")),
    ("battery-pack", "Battery · four-cell series pack cycle", include_str!("../data/templates/battery-pack.json")),
    ("bioheat-hotspot", "Bioheat · SAR hot spot (Pennes)", include_str!("../data/templates/bioheat-hotspot.json")),
    ("cahn-hilliard", "Microstructure · Cahn–Hilliard spinodal decomposition", include_str!("../data/templates/cahn-hilliard.json")),
    ("chemistry-batch", "Chemistry · isothermal batch", include_str!("../data/templates/chemistry-batch.json")),
    ("cosserat-cantilever", "Rods · Cosserat cantilever under gravity", include_str!("../data/templates/cosserat-cantilever.json")),
    ("dem-settling", "Granular · DEM spheres settling in a box", include_str!("../data/templates/dem-settling.json")),
    ("dg-advection", "High order · DG advection", include_str!("../data/templates/dg-advection.json")),
    ("diffusion2d", "Diffusion · structured 2D", include_str!("../data/templates/diffusion2d.json")),
    ("edge-cavity", "EM · Nédélec FEM cavity port impedance", include_str!("../data/templates/edge-cavity.json")),
    ("elasticity-cantilever", "Structures · plate cantilever under gravity (FEM)", include_str!("../data/templates/elasticity-cantilever.json")),
    ("fdtd3d-cavity", "EM · 3D FDTD cavity pulse", include_str!("../data/templates/fdtd3d-cavity.json")),
    ("flow-cavity-pimple", "Flow · transient cavity (PIMPLE)", include_str!("../data/templates/flow-cavity-pimple.json")),
    ("flow-cavity-simple", "Flow · steady cavity (SIMPLE)", include_str!("../data/templates/flow-cavity-simple.json")),
    ("flow-channel", "Flow · channel (PISO)", include_str!("../data/templates/flow-channel.json")),
    (
        "flow-electronics-enclosure",
        "Flow · primitive electronics enclosure (isothermal)",
        include_str!("../data/templates/flow-electronics-enclosure.json"),
    ),
    (
        "flow-electronics-enclosure-cht",
        "Heat · fan-cooled enclosure with conducting, powered components",
        include_str!("../data/templates/flow-electronics-enclosure-cht.json"),
    ),
    (
        "flow-electronics-enclosure-fan",
        "Flow · fan-curve enclosure with exhaust filter (isothermal)",
        include_str!("../data/templates/flow-electronics-enclosure-fan.json"),
    ),
    (
        "flow-electronics-enclosure-natural",
        "Heat · sealed enclosure, natural convection",
        include_str!("../data/templates/flow-electronics-enclosure-natural.json"),
    ),
    ("galvanic-corrosion", "Electrochemistry · galvanic zinc/oxygen couple", include_str!("../data/templates/galvanic-corrosion.json")),
    ("heat2d", "Heat · transient FEM 2D", include_str!("../data/templates/heat2d.json")),
    ("j2-cyclic", "Plasticity · J2 cyclic uniaxial stress", include_str!("../data/templates/j2-cyclic.json")),
    ("joule-busbar", "Coupled · busbar Joule heating", include_str!("../data/templates/joule-busbar.json")),
    ("landau-damping", "Kinetic · Landau damping (Vlasov–Poisson)", include_str!("../data/templates/landau-damping.json")),
    ("lbm-taylor-green", "Flow · D2Q9 Taylor–Green", include_str!("../data/templates/lbm-taylor-green.json")),
    ("lennard-jones-md", "Atomistics · Lennard-Jones chain MD", include_str!("../data/templates/lennard-jones-md.json")),
    ("maxwell1d", "EM · Yee Maxwell 1D", include_str!("../data/templates/maxwell1d.json")),
    ("multiphase-bubble-column", "Multiphase · Euler–Euler bubble column", include_str!("../data/templates/multiphase-bubble-column.json")),
    ("multiphase-vof-stir", "Multiphase · stirred compressible VOF", include_str!("../data/templates/multiphase-vof-stir.json")),
    ("optics-lens", "Optics · sequential lens", include_str!("../data/templates/optics-lens.json")),
    ("pdn-ir-drop", "SI/PI · power-plane IR drop", include_str!("../data/templates/pdn-ir-drop.json")),
    ("phase-field-bar", "Fracture · phase-field notched bar", include_str!("../data/templates/phase-field-bar.json")),
    ("photonic-waveguide", "Photonics · silicon wire waveguide modes", include_str!("../data/templates/photonic-waveguide.json")),
    ("pic-two-stream", "Plasma · two-stream instability (PIC)", include_str!("../data/templates/pic-two-stream.json")),
    ("poisson1d", "Scalar · Poisson 1D", include_str!("../data/templates/poisson1d.json")),
    ("poisson2d", "Scalar · Poisson FEM 2D", include_str!("../data/templates/poisson2d.json")),
    ("poisson3d", "Scalar · Poisson FEM 3D", include_str!("../data/templates/poisson3d.json")),
    ("poly-scalar", "Transport · scalar FVM 3D", include_str!("../data/templates/poly-scalar.json")),
    ("power-network", "Power · five-bus flow, PV curve, dispatch and swing", include_str!("../data/templates/power-network.json")),
    ("rans-kepsilon", "Turbulence · k–epsilon transport", include_str!("../data/templates/rans-kepsilon.json")),
    ("rans-reynolds-stress", "Turbulence · Reynolds-stress transport", include_str!("../data/templates/rans-reynolds-stress.json")),
    ("rans-sa", "Turbulence · Spalart–Allmaras transport", include_str!("../data/templates/rans-sa.json")),
    ("rans-sst-des", "Turbulence · SST/DES transport", include_str!("../data/templates/rans-sst-des.json")),
    ("rf-tem", "RF · TEM two-port sweep", include_str!("../data/templates/rf-tem.json")),
    ("scalar-mixed", "Scalar · mixed-boundary FEM 2D", include_str!("../data/templates/scalar-mixed.json")),
    ("sod-weno", "Flow · Sod shock tube (WENO)", include_str!("../data/templates/sod-weno.json")),
    ("solid-heat", "Heat · solid conduction FVM 3D", include_str!("../data/templates/solid-heat.json")),
    ("spice-rc-filter", "Circuits · SPICE RC filter with diode clamp", include_str!("../data/templates/spice-rc-filter.json")),
    ("spray-evaporation", "Spray · evaporating injection in a cross-flow", include_str!("../data/templates/spray-evaporation.json")),
    ("tcad-pn-diode", "Semiconductors · PN diode IV and CV", include_str!("../data/templates/tcad-pn-diode.json")),
    ("thermal-radiation", "Heat · prescribed flow and radiation", include_str!("../data/templates/thermal-radiation.json")),
    ("thmc-column", "Porous · coupled THMC column", include_str!("../data/templates/thmc-column.json")),
    ("transport-depth-dose", "Radiation · slab depth dose (condensed history)", include_str!("../data/templates/transport-depth-dose.json")),
    ("wire-dipole", "Antennas · half-wave dipole impedance sweep", include_str!("../data/templates/wire-dipole.json")),
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
