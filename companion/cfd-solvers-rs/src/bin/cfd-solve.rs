use cfd_solvers::{
    case::{Case, SOLVERS},
    Result, SolverError, UPSTREAM_REVISION,
};
use std::path::PathBuf;
const HELP:&str="cfd-solve — native Rust reference solver port\n\n  cfd-solve list [--json]\n  cfd-solve run CASE.json --out DIRECTORY\n  cfd-solve example NAME --out DIRECTORY\n  cfd-solve coverage [--json]\n  cfd-solve verify --out DIRECTORY   (offline self-checks)\n  cfd-solve --version\n\nExamples: poisson1d, poisson2d, poisson3d, scalar-mixed, heat2d,\n          lbm-taylor-green, maxwell1d, sod-weno, diffusion2d,\n          rf-tem, optics-lens, chemistry-batch, acoustic-kspace2d,\n          poly-scalar, flow-channel, flow-electronics-enclosure,\n          flow-electronics-enclosure-fan, flow-electronics-enclosure-cht,\n          flow-electronics-enclosure-natural,\n          flow-cavity-simple, flow-cavity-pimple,\n          solid-heat, thermal-radiation, rans-sa, rans-kepsilon,\n          rans-sst-des, rans-reynolds-stress,\n          spice-rc-filter, tcad-pn-diode, battery-discharge, battery-pack,\n          wire-dipole, fdtd3d-cavity, edge-cavity, pdn-ir-drop,\n          galvanic-corrosion, elasticity-cantilever, bar-modes, j2-cyclic,\n          phase-field-bar, cosserat-cantilever, dem-settling, pic-two-stream,\n          transport-depth-dose, bioheat-hotspot, joule-busbar, cahn-hilliard,\n          dg-advection, landau-damping, lennard-jones-md, photonic-waveguide,\n          thmc-column, power-network, multiphase-bubble-column,\n          multiphase-vof-stir, spray-evaporation\n\nOutput includes manifest.json and solver-specific VTU/PVD, CSV, S2P or JSON.\nA new output directory is required for each completed run.\nThis release does not yet provide full upstream solver/API/backend parity.\n";
fn example(name: &str) -> Option<&'static str> {
    Some(match name {
        "poisson1d" => include_str!("../../examples/poisson1d.json"),
        "poisson2d" => include_str!("../../examples/poisson2d.json"),
        "poisson3d" => include_str!("../../examples/poisson3d.json"),
        "scalar-mixed" => include_str!("../../examples/scalar-mixed.json"),
        "heat2d" => include_str!("../../examples/heat2d.json"),
        "lbm-taylor-green" => include_str!("../../examples/lbm-taylor-green.json"),
        "maxwell1d" => include_str!("../../examples/maxwell1d.json"),
        "sod-weno" => include_str!("../../examples/sod-weno.json"),
        "diffusion2d" => include_str!("../../examples/diffusion2d.json"),
        "rf-tem" => include_str!("../../examples/rf-tem.json"),
        "optics-lens" => include_str!("../../examples/optics-lens.json"),
        "chemistry-batch" => include_str!("../../examples/chemistry-batch.json"),
        "acoustic-kspace2d" => include_str!("../../examples/acoustic-kspace2d.json"),
        "poly-scalar" => include_str!("../../examples/poly-scalar.json"),
        "flow-channel" => include_str!("../../examples/flow-channel.json"),
        "flow-electronics-enclosure" => {
            include_str!("../../examples/flow-electronics-enclosure.json")
        }
        "flow-electronics-enclosure-fan" => {
            include_str!("../../examples/flow-electronics-enclosure-fan.json")
        }
        "flow-electronics-enclosure-cht" => {
            include_str!("../../examples/flow-electronics-enclosure-cht.json")
        }
        "flow-electronics-enclosure-natural" => {
            include_str!("../../examples/flow-electronics-enclosure-natural.json")
        }
        "flow-cavity-simple" => include_str!("../../examples/flow-cavity-simple.json"),
        "flow-cavity-pimple" => include_str!("../../examples/flow-cavity-pimple.json"),
        "solid-heat" => include_str!("../../examples/solid-heat.json"),
        "thermal-radiation" => include_str!("../../examples/thermal-radiation.json"),
        "rans-sa" => include_str!("../../examples/rans-sa.json"),
        "rans-kepsilon" => include_str!("../../examples/rans-kepsilon.json"),
        "rans-sst-des" => include_str!("../../examples/rans-sst-des.json"),
        "rans-reynolds-stress" => include_str!("../../examples/rans-reynolds-stress.json"),
        "spice-rc-filter" => include_str!("../../examples/spice-rc-filter.json"),
        "tcad-pn-diode" => include_str!("../../examples/tcad-pn-diode.json"),
        "battery-discharge" => include_str!("../../examples/battery-discharge.json"),
        "battery-pack" => include_str!("../../examples/battery-pack.json"),
        "wire-dipole" => include_str!("../../examples/wire-dipole.json"),
        "fdtd3d-cavity" => include_str!("../../examples/fdtd3d-cavity.json"),
        "edge-cavity" => include_str!("../../examples/edge-cavity.json"),
        "pdn-ir-drop" => include_str!("../../examples/pdn-ir-drop.json"),
        "galvanic-corrosion" => include_str!("../../examples/galvanic-corrosion.json"),
        "elasticity-cantilever" => include_str!("../../examples/elasticity-cantilever.json"),
        "bar-modes" => include_str!("../../examples/bar-modes.json"),
        "j2-cyclic" => include_str!("../../examples/j2-cyclic.json"),
        "phase-field-bar" => include_str!("../../examples/phase-field-bar.json"),
        "cosserat-cantilever" => include_str!("../../examples/cosserat-cantilever.json"),
        "dem-settling" => include_str!("../../examples/dem-settling.json"),
        "pic-two-stream" => include_str!("../../examples/pic-two-stream.json"),
        "transport-depth-dose" => include_str!("../../examples/transport-depth-dose.json"),
        "bioheat-hotspot" => include_str!("../../examples/bioheat-hotspot.json"),
        "joule-busbar" => include_str!("../../examples/joule-busbar.json"),
        "cahn-hilliard" => include_str!("../../examples/cahn-hilliard.json"),
        "dg-advection" => include_str!("../../examples/dg-advection.json"),
        "landau-damping" => include_str!("../../examples/landau-damping.json"),
        "lennard-jones-md" => include_str!("../../examples/lennard-jones-md.json"),
        "photonic-waveguide" => include_str!("../../examples/photonic-waveguide.json"),
        "thmc-column" => include_str!("../../examples/thmc-column.json"),
        "power-network" => include_str!("../../examples/power-network.json"),
        "multiphase-bubble-column" => include_str!("../../examples/multiphase-bubble-column.json"),
        "multiphase-vof-stir" => include_str!("../../examples/multiphase-vof-stir.json"),
        "spray-evaporation" => include_str!("../../examples/spray-evaporation.json"),
        _ => return None,
    })
}
fn run() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let Some(command) = args.first() else {
        print!("{HELP}");
        return Ok(());
    };
    match command.as_str() {
        "--help" | "-h" => {
            print!("{HELP}");
        }
        "--version" | "-V" => println!(
            "cfd-solvers-rs {} (upstream {UPSTREAM_REVISION})",
            env!("CARGO_PKG_VERSION")
        ),
        "list" | "--list" => {
            if args.get(1).is_some_and(|s| s == "--json") {
                println!("{}", serde_json::to_string_pretty(&SOLVERS[..])?);
            } else {
                for solver in SOLVERS {
                    println!("{solver}");
                }
            }
        }
        "coverage" => {
            let ledger: serde_json::Value =
                serde_json::from_str(include_str!("../../parity/status.json"))?;
            if args.get(1).is_some_and(|s| s == "--json") {
                println!("{}", serde_json::to_string_pretty(&ledger)?);
            } else {
                println!("Complete clone: false\nUpstream public units: {}\nStatus counts: {}\nSee parity/status.json for every method family and limitation.",ledger["unit_count"],ledger["status_counts"]);
            }
        }
        "verify" => {
            if args.len() != 3 || args[1] != "--out" {
                return Err(SolverError::InvalidInput(
                    "usage: cfd-solve verify --out DIRECTORY".into(),
                ));
            }
            let report = cfd_solvers::case::verification::run_suite(&PathBuf::from(&args[2]))?;
            for check in report["checks"].as_array().into_iter().flatten() {
                println!(
                    "{} {}",
                    if check["passed"] == true {
                        "PASS"
                    } else {
                        "FAIL"
                    },
                    check["name"].as_str().unwrap_or("?")
                );
            }
            if report["passed"] != true {
                return Err(SolverError::NumericalFailure(
                    "verification failed; see verification.json".into(),
                ));
            }
        }
        "run" | "example" => {
            if args.len() != 4 || args[2] != "--out" {
                return Err(SolverError::InvalidInput(
                    "usage: cfd-solve run CASE.json --out DIRECTORY, or cfd-solve example NAME --out DIRECTORY".into(),
                ));
            }
            let input = if command == "run" {
                std::fs::read_to_string(&args[1])?
            } else {
                example(&args[1])
                    .ok_or_else(|| {
                        SolverError::Unsupported(format!(
                            "example {}; use --help for available examples",
                            args[1]
                        ))
                    })?
                    .to_string()
            };
            let case: Case = serde_json::from_str(&input)?;
            let folder = PathBuf::from(&args[3]);
            let manifest = case.run(&folder)?;
            println!("{}", serde_json::to_string_pretty(&manifest)?);
        }
        _ => {
            return Err(SolverError::Unsupported(format!(
                "command {command}; use --help"
            )))
        }
    }
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("cfd-solve: {error}");
        std::process::exit(2);
    }
}
