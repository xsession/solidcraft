//! Electrical and electronics case adapters: SPICE netlists, 1-D semiconductor diodes,
//! lithium-ion cells and series packs, thin-wire antennas, 3-D Yee FDTD, Nedelec cavity
//! FEM, PDN IR drop and galvanic corrosion.
use super::{csv, json, linspace, logspace, metric};
use crate::{
    battery::*,
    case::RunManifest,
    circuit::{Circuit, NewtonConfig, TransientMethod},
    core::Vec3,
    electrochemistry::*,
    em::*,
    fdtd::{Boundary3D, Maxwell3D, Maxwell3DConfig},
    fem::make_box_tet_mesh,
    io::{DataArray, VtuMesh},
    require,
    tcad::*,
    Result, SolverError,
};
use num_complex::Complex64 as C;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::path::Path;

const MAX_NETLIST_BYTES: usize = 64 * 1024;
const MAX_SWEEP_POINTS: usize = 2001;

fn points_ok(points: usize, what: &str) -> Result<()> {
    require(
        (1..=MAX_SWEEP_POINTS).contains(&points),
        format!("{what} needs 1..={MAX_SWEEP_POINTS} points"),
    )
}

// ---------------------------------------------------------------- SPICE

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpiceIntegration {
    BackwardEuler,
    #[default]
    Trapezoidal,
    Bdf2,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum SpiceAnalysis {
    Op,
    Ac {
        start_hz: f64,
        stop_hz: f64,
        points: usize,
    },
    Transient {
        dt_s: f64,
        steps: usize,
        #[serde(default)]
        method: SpiceIntegration,
    },
    DcSweep {
        source: String,
        start_v: f64,
        stop_v: f64,
        points: usize,
    },
}

/// A SPICE netlist with one or more analyses. Probes are node names; empty means every node.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpiceCircuitCase {
    pub netlist: String,
    pub analyses: Vec<SpiceAnalysis>,
    #[serde(default)]
    pub probes: Vec<String>,
}

impl SpiceCircuitCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.netlist.len() <= MAX_NETLIST_BYTES,
            "netlist exceeds 64 KiB",
        )?;
        require(
            !self.analyses.is_empty() && self.analyses.len() <= 8,
            "choose 1 to 8 analyses",
        )?;
        let mut circuit = Circuit::parse_spice(&self.netlist)?;
        let names: Vec<String> = circuit.node_names().to_vec();
        let probes: Vec<(String, usize)> = if self.probes.is_empty() {
            names
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, n)| (n.clone(), i))
                .collect()
        } else {
            self.probes
                .iter()
                .map(|p| Ok((p.clone(), circuit.find_node(p)?)))
                .collect::<Result<_>>()?
        };
        require(
            !probes.is_empty() && probes.len() <= 32,
            "choose 1 to 32 probe nodes",
        )?;
        let cfg = NewtonConfig::default();
        let mut linear_ok = true;
        for analysis in &self.analyses {
            match analysis {
                SpiceAnalysis::Op => {
                    let op = circuit.dc_operating_point_default_homotopy(&cfg)?;
                    linear_ok &= op.converged;
                    let voltages: serde_json::Map<String, serde_json::Value> = names
                        .iter()
                        .zip(&op.node_voltage)
                        .map(|(n, v)| (n.clone(), json!(v)))
                        .collect();
                    json(
                        out,
                        m,
                        "operating_point.json",
                        &json!({"node_voltage_V": voltages, "branch_current_A": op.branch_current, "iterations": op.iterations, "converged": op.converged}),
                    )?;
                    metric(m, "op_newton_iterations", op.iterations as f64);
                    for (name, index) in &probes {
                        metric(m, &format!("op_V_{name}"), op.node_voltage[*index]);
                    }
                }
                SpiceAnalysis::Ac {
                    start_hz,
                    stop_hz,
                    points,
                } => {
                    points_ok(*points, "AC sweep")?;
                    let sweep = circuit.ac_log_sweep(*start_hz, *stop_hz, *points, &cfg)?;
                    let f: Vec<f64> = sweep.iter().map(|s| s.frequency_hz).collect();
                    let mut cols: Vec<(String, Vec<f64>)> = vec![];
                    for (name, index) in &probes {
                        let v: Vec<C> = sweep.iter().map(|s| s.node_voltage[*index]).collect();
                        cols.push((
                            format!("{name}_dB"),
                            v.iter()
                                .map(|z| 20.0 * z.norm().max(1e-300).log10())
                                .collect(),
                        ));
                        cols.push((
                            format!("{name}_phase_deg"),
                            v.iter().map(|z| z.arg().to_degrees()).collect(),
                        ));
                    }
                    let mut table: Vec<(&str, &[f64])> = vec![("frequency_Hz", &f)];
                    table.extend(cols.iter().map(|(n, v)| (n.as_str(), v.as_slice())));
                    csv(out, m, "ac.csv", &table)?;
                }
                SpiceAnalysis::Transient {
                    dt_s,
                    steps,
                    method,
                } => {
                    require(
                        *steps >= 1 && *steps <= 20_000,
                        "transient needs 1..=20000 steps",
                    )?;
                    let method = match method {
                        SpiceIntegration::BackwardEuler => TransientMethod::BackwardEuler,
                        SpiceIntegration::Trapezoidal => TransientMethod::Trapezoidal,
                        SpiceIntegration::Bdf2 => TransientMethod::Bdf2,
                    };
                    let pts = circuit.transient(*dt_s, *steps, method, &cfg)?;
                    let t: Vec<f64> = pts.iter().map(|p| p.time).collect();
                    let cols: Vec<(String, Vec<f64>)> = probes
                        .iter()
                        .map(|(n, i)| {
                            (
                                format!("V_{n}"),
                                pts.iter().map(|p| p.node_voltage[*i]).collect(),
                            )
                        })
                        .collect();
                    let mut table: Vec<(&str, &[f64])> = vec![("time_s", &t)];
                    table.extend(cols.iter().map(|(n, v)| (n.as_str(), v.as_slice())));
                    csv(out, m, "transient.csv", &table)?;
                    metric(m, "transient_end_time_s", t.last().copied().unwrap_or(0.0));
                }
                SpiceAnalysis::DcSweep {
                    source,
                    start_v,
                    stop_v,
                    points,
                } => {
                    points_ok(*points, "DC sweep")?;
                    let sweep = circuit
                        .dc_sweep_voltage_source(source, *start_v, *stop_v, *points, &cfg)?;
                    linear_ok &= sweep.iter().all(|p| p.operating_point.converged);
                    let s: Vec<f64> = sweep.iter().map(|p| p.source_value).collect();
                    let cols: Vec<(String, Vec<f64>)> = probes
                        .iter()
                        .map(|(n, i)| {
                            (
                                format!("V_{n}"),
                                sweep
                                    .iter()
                                    .map(|p| p.operating_point.node_voltage[*i])
                                    .collect(),
                            )
                        })
                        .collect();
                    let label = format!("{source}_V");
                    let mut table: Vec<(&str, &[f64])> = vec![(label.as_str(), &s)];
                    table.extend(cols.iter().map(|(n, v)| (n.as_str(), v.as_slice())));
                    csv(out, m, "dc_sweep.csv", &table)?;
                }
            }
        }
        m.linear_solver_converged = Some(linear_ok);
        Ok(())
    }
}

// ---------------------------------------------------------------- TCAD

fn default_tcad_nodes() -> usize {
    101
}
fn default_area() -> f64 {
    1e-12
}
fn default_temperature() -> f64 {
    300.0
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VoltageSweep {
    pub start_v: f64,
    pub stop_v: f64,
    pub points: usize,
}

/// Silicon PN (or PIN, when `intrinsic_m` is positive) diode on a 1-D drift-diffusion mesh.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TcadDiodeCase {
    pub length_m: f64,
    #[serde(default = "default_tcad_nodes")]
    pub nodes: usize,
    #[serde(default = "default_area")]
    pub area_m2: f64,
    pub acceptor_m3: f64,
    pub donor_m3: f64,
    /// Junction at the centre; an intrinsic layer of this width is centred on it.
    #[serde(default)]
    pub intrinsic_m: f64,
    #[serde(default = "default_temperature")]
    pub temperature_k: f64,
    pub iv: VoltageSweep,
    #[serde(default)]
    pub cv: Option<VoltageSweep>,
}

impl TcadDiodeCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.nodes >= 11 && self.nodes <= 2001,
            "diode needs 11..=2001 nodes",
        )?;
        points_ok(self.iv.points, "IV sweep")?;
        let cfg = Semiconductor1DConfig {
            length_m: self.length_m,
            nodes: self.nodes,
            area_m2: self.area_m2,
            max_gummel_iterations: 1000,
            relative_tolerance: 2e-7,
            under_relaxation: 0.05,
            ..Default::default()
        };
        let mut d = SemiconductorDevice1D::new(cfg)?;
        if self.temperature_k != 300.0 {
            d.set_temperature(self.temperature_k)?;
        }
        let x = d.x().to_vec();
        let doping = if self.intrinsic_m > 0.0 {
            let mid = 0.5 * self.length_m;
            pin_diode_doping(
                &x,
                mid - 0.5 * self.intrinsic_m,
                mid + 0.5 * self.intrinsic_m,
                self.acceptor_m3,
                self.donor_m3,
                0.0,
            )?
        } else {
            pn_junction_doping(&x, 0.5 * self.length_m, self.acceptor_m3, self.donor_m3)?
        };
        d.set_net_doping(&doping)?;
        let eq = d.solve_equilibrium()?;
        let n = self.nodes;
        metric(
            m,
            "built_in_potential_V",
            d.potential()[n - 1] - d.potential()[0],
        );
        csv(
            out,
            m,
            "equilibrium.csv",
            &[
                ("x_m", &x),
                ("potential_V", d.potential()),
                ("electron_m3", d.electron_density()),
                ("hole_m3", d.hole_density()),
                ("net_doping_m3", &doping),
            ],
        )?;
        // The p side is the left contact (anode). Sweeping the right (cathode) contact by -V
        // applies +V to the anode; the reported terminal current is positive in forward bias.
        let iv = d.clone().sweep_right_contact(
            -self.iv.start_v,
            -self.iv.stop_v,
            self.iv.points.max(2),
        )?;
        let v: Vec<f64> = iv.iter().map(|p| -p.voltage_v).collect();
        let i: Vec<f64> = iv.iter().map(|p| p.current_a).collect();
        let ok: Vec<f64> = iv
            .iter()
            .map(|p| if p.converged { 1.0 } else { 0.0 })
            .collect();
        csv(
            out,
            m,
            "iv.csv",
            &[
                ("anode_voltage_V", &v),
                ("anode_current_A", &i),
                ("converged", &ok),
            ],
        )?;
        let mut converged = eq.converged && iv.iter().all(|p| p.converged);
        if let Some(cv) = &self.cv {
            points_ok(cv.points, "CV sweep")?;
            let c = d.sweep_cv(cv.start_v, cv.stop_v, cv.points.max(2), 1e-3)?;
            converged &= c.iter().all(|p| p.converged);
            let v: Vec<f64> = c.iter().map(|p| p.voltage_v).collect();
            let cap: Vec<f64> = c.iter().map(|p| p.capacitance_f).collect();
            csv(
                out,
                m,
                "cv.csv",
                &[("anode_voltage_V", &v), ("capacitance_F", &cap)],
            )?;
        }
        metric(m, "equilibrium_gummel_iterations", eq.iterations as f64);
        if let (Some(a), Some(b)) = (i.first(), i.last()) {
            metric(m, "anode_current_at_start_A", *a);
            metric(m, "anode_current_at_stop_A", *b);
        }
        m.linear_solver_converged = Some(converged);
        Ok(())
    }
}

// ---------------------------------------------------------------- battery

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BatteryModel {
    Spm,
    #[default]
    Spme,
    Dfn,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CurrentSegment {
    /// Positive discharges.
    pub current_a: f64,
    pub duration_s: f64,
}

fn default_cutoff_low() -> f64 {
    2.8
}
fn default_cutoff_high() -> f64 {
    4.25
}
fn default_ambient() -> f64 {
    298.15
}
fn default_dt() -> f64 {
    5.0
}
fn one_cell() -> usize {
    1
}

/// Graphite/NMC cell (SPM, SPMe or P2D DFN) or a balanced series pack of SPMe cells under a
/// piecewise-constant current protocol with voltage cut-offs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BatteryCellCase {
    #[serde(default)]
    pub model: BatteryModel,
    pub protocol: Vec<CurrentSegment>,
    #[serde(default = "default_dt")]
    pub dt_s: f64,
    #[serde(default = "default_cutoff_low")]
    pub cutoff_low_v: f64,
    #[serde(default = "default_cutoff_high")]
    pub cutoff_high_v: f64,
    #[serde(default = "default_ambient")]
    pub ambient_k: f64,
    #[serde(default = "one_cell")]
    pub cells_in_series: usize,
}

struct BatterySample {
    t: f64,
    i: f64,
    v: f64,
    ocv: f64,
    temp: f64,
    heat: f64,
}

impl BatteryCellCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            !self.protocol.is_empty() && self.protocol.len() <= 64,
            "battery protocol needs 1..=64 segments",
        )?;
        require(
            self.dt_s > 0.0 && self.dt_s.is_finite(),
            "battery time step must be positive",
        )?;
        require(
            self.cutoff_low_v < self.cutoff_high_v,
            "battery cut-offs must be ordered",
        )?;
        require(
            self.cells_in_series >= 1 && self.cells_in_series <= 64,
            "1..=64 cells in series",
        )?;
        let total: f64 = self.protocol.iter().map(|s| s.duration_s).sum();
        require(
            total.is_finite() && total > 0.0 && total / self.dt_s <= 20_000.0,
            "protocol exceeds 20000 steps",
        )?;
        let cell = default_graphite_nmc_config();
        let n = self.cells_in_series as f64;
        let mut samples = Vec::new();
        let mut stop = None;
        let mut charge_ah = 0.0;
        let mut energy_wh = 0.0;
        enum Engine {
            Spm(Box<SingleParticleModel>),
            Spme(Box<SingleParticleElectrolyteModel>),
            Dfn(Box<DoyleFullerNewmanModel>),
            Pack(Box<SeriesBatteryPack>),
        }
        let mut engine = if self.cells_in_series > 1 {
            // The pack's own protection sits outside the case cut-offs, which stop the run first.
            let control = PackControl {
                minimum_cell_voltage_v: self.cutoff_low_v - 0.5,
                maximum_cell_voltage_v: self.cutoff_high_v + 0.5,
                ..Default::default()
            };
            let mut pack =
                SeriesBatteryPack::new(&vec![cell.clone(); self.cells_in_series], control)?;
            pack.set_external_temperatures(&vec![self.ambient_k; self.cells_in_series])?;
            Engine::Pack(Box::new(pack))
        } else {
            match self.model {
                BatteryModel::Spm => {
                    let mut c = SingleParticleModel::new(cell)?;
                    c.set_external_temperature(self.ambient_k)?;
                    Engine::Spm(Box::new(c))
                }
                BatteryModel::Spme => {
                    let mut c = SingleParticleElectrolyteModel::new(cell)?;
                    c.set_external_temperature(self.ambient_k)?;
                    Engine::Spme(Box::new(c))
                }
                BatteryModel::Dfn => Engine::Dfn(Box::new(DoyleFullerNewmanModel::new(
                    default_graphite_nmc_dfn_config(),
                )?)),
            }
        };
        let cell_sample = |r: BatteryStepResult| BatterySample {
            t: r.time_s,
            i: r.current_a,
            v: r.terminal_voltage_v,
            ocv: r.open_circuit_voltage_v,
            temp: r.temperature_k,
            heat: r.heat_generation_w,
        };
        'protocol: for segment in &self.protocol {
            require(
                segment.current_a.is_finite() && segment.duration_s > 0.0,
                "invalid protocol segment",
            )?;
            let steps = (segment.duration_s / self.dt_s).ceil().max(1.0) as usize;
            let dt = segment.duration_s / steps as f64;
            for _ in 0..steps {
                let s = match &mut engine {
                    Engine::Spm(c) => cell_sample(c.step(segment.current_a, dt)?),
                    Engine::Spme(c) => cell_sample(c.step(segment.current_a, dt)?),
                    Engine::Dfn(c) => {
                        let r = c.step(segment.current_a, dt)?;
                        BatterySample {
                            t: r.time_s,
                            i: r.current_a,
                            v: r.terminal_voltage_v,
                            ocv: r.open_circuit_voltage_v,
                            temp: r.temperature_k,
                            heat: 0.0,
                        }
                    }
                    Engine::Pack(p) => {
                        let r = p.step(segment.current_a, dt)?;
                        let ocv = r.cells.iter().map(|c| c.open_circuit_voltage_v).sum();
                        let temp = r
                            .cells
                            .iter()
                            .map(|c| c.temperature_k)
                            .fold(f64::MIN, f64::max);
                        let heat = r.cells.iter().map(|c| c.heat_generation_w).sum();
                        BatterySample {
                            t: r.time_s,
                            i: segment.current_a,
                            v: r.terminal_voltage_v,
                            ocv,
                            temp,
                            heat,
                        }
                    }
                };
                charge_ah += s.i * dt / 3600.0;
                energy_wh += s.i * s.v * dt / 3600.0;
                let per_cell = s.v / n;
                samples.push(s);
                if per_cell <= self.cutoff_low_v {
                    stop = Some("low_voltage_cutoff");
                    break 'protocol;
                }
                if per_cell >= self.cutoff_high_v {
                    stop = Some("high_voltage_cutoff");
                    break 'protocol;
                }
            }
        }
        let col = |f: fn(&BatterySample) -> f64| samples.iter().map(f).collect::<Vec<f64>>();
        csv(
            out,
            m,
            "cycle.csv",
            &[
                ("time_s", &col(|s| s.t)),
                ("current_A", &col(|s| s.i)),
                ("terminal_voltage_V", &col(|s| s.v)),
                ("open_circuit_voltage_V", &col(|s| s.ocv)),
                ("temperature_K", &col(|s| s.temp)),
                ("heat_generation_W", &col(|s| s.heat)),
            ],
        )?;
        metric(m, "net_discharge_Ah", charge_ah);
        metric(m, "net_discharge_energy_Wh", energy_wh);
        metric(
            m,
            "minimum_voltage_V",
            samples.iter().map(|s| s.v).fold(f64::INFINITY, f64::min),
        );
        metric(m, "final_voltage_V", samples.last().map_or(0.0, |s| s.v));
        metric(
            m,
            "maximum_temperature_K",
            samples.iter().map(|s| s.temp).fold(f64::MIN, f64::max),
        );
        metric(
            m,
            "stopped_on_cutoff",
            if stop.is_some() { 1.0 } else { 0.0 },
        );
        json(
            out,
            m,
            "cycle_summary.json",
            &json!({"model": self.model, "cells_in_series": self.cells_in_series, "stop_reason": stop.unwrap_or("protocol_complete"), "samples": samples.len(), "dfn_thermal_note": "DFN runs keep their configured temperature; heat generation is reported as 0 for DFN."}),
        )?;
        Ok(())
    }
}

// ---------------------------------------------------------------- thin-wire antenna

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FrequencySweep {
    pub start_hz: f64,
    pub stop_hz: f64,
    pub points: usize,
    #[serde(default)]
    pub logarithmic: bool,
}

impl FrequencySweep {
    fn values(&self) -> Result<Vec<f64>> {
        points_ok(self.points, "frequency sweep")?;
        require(
            self.start_hz > 0.0 && self.stop_hz >= self.start_hz && self.stop_hz.is_finite(),
            "invalid frequency sweep",
        )?;
        Ok(if self.logarithmic {
            logspace(self.start_hz, self.stop_hz, self.points)
        } else {
            linspace(self.start_hz, self.stop_hz, self.points)
        })
    }
}

fn odd_segments() -> usize {
    31
}
fn fifty() -> f64 {
    50.0
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WireGround {
    /// Centre-fed dipole of total length `length_m`.
    #[default]
    FreeSpace,
    /// Base-fed monopole of height `length_m` on an infinite PEC plane (image theory).
    PecPlane,
}

/// Straight, centre-fed thin-wire dipole (or monopole over a PEC plane) by Harrington's
/// pulse-basis moment method with point matching, a reduced kernel and finite-difference
/// charge terms. The upstream point-matched Pocklington reference (`rf::thin_wire`) remains
/// a library API; it is not used here because its feed reactance is not reliable.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WireAntennaCase {
    pub length_m: f64,
    pub radius_m: f64,
    #[serde(default = "odd_segments")]
    pub segments: usize,
    #[serde(default)]
    pub ground: WireGround,
    pub sweep: FrequencySweep,
    #[serde(default = "fifty")]
    pub reference_impedance_ohm: f64,
}

const EPS0: f64 = 8.8541878128e-12;
const MU0: f64 = 1.25663706212e-6;

/// `(1/(4πΔ)) ∫_{z1}^{z2} exp(-jkR)/R dz'`, `R² = (z - z')² + a²`: the static part in closed
/// form and the smooth remainder by 8-point Gauss–Legendre on 4 panels.
fn psi(z: f64, z1: f64, z2: f64, a: f64, k: f64) -> C {
    const X: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const W: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    let delta = z2 - z1;
    let stat = ((z2 - z) / a).asinh() - ((z1 - z) / a).asinh();
    let mut smooth = C::new(0.0, 0.0);
    let panels = 4;
    let h = delta / panels as f64;
    for p in 0..panels {
        let (lo, hi) = (z1 + h * p as f64, z1 + h * (p + 1) as f64);
        let (mid, half) = (0.5 * (lo + hi), 0.5 * (hi - lo));
        for (x, w) in X.iter().zip(W) {
            for sgn in [-1.0, 1.0] {
                let zp = mid + sgn * x * half;
                let r = ((z - zp) * (z - zp) + a * a).sqrt();
                let term = if k * r < 1e-6 {
                    C::new(0.0, -k)
                } else {
                    (C::new(0.0, -k * r).exp() - 1.0) / r
                };
                smooth += term * (w * half);
            }
        }
    }
    (C::new(stat, 0.0) + smooth) / (4.0 * std::f64::consts::PI * delta)
}

/// Feed impedance and segment currents of a centre-fed straight dipole.
pub fn harrington_dipole(
    length_m: f64,
    radius_m: f64,
    segments: usize,
    frequency_hz: f64,
) -> Result<(C, Vec<f64>, Vec<C>)> {
    require(
        length_m > 0.0
            && radius_m > 0.0
            && radius_m < 0.05 * length_m
            && frequency_hz > 0.0
            && length_m.is_finite()
            && frequency_hz.is_finite(),
        "dipole needs positive length, frequency and a thin radius",
    )?;
    let n = segments;
    let d = length_m / n as f64;
    require(
        radius_m < 0.5 * d,
        "wire radius must be below half a segment length; use fewer segments",
    )?;
    let w = 2.0 * std::f64::consts::PI * frequency_hz;
    let k = w * (MU0 * EPS0).sqrt();
    let zc: Vec<f64> = (0..n)
        .map(|i| -0.5 * length_m + (i as f64 + 0.5) * d)
        .collect();
    let jwmu = C::new(0.0, w * MU0);
    let inv_jwe = C::new(0.0, -1.0 / (w * EPS0));
    let mut a = vec![C::new(0.0, 0.0); n * n];
    for mi in 0..n {
        for ni in 0..n {
            let zm = zc[mi];
            let zn = zc[ni];
            let p = |obs: f64, src: f64| psi(obs, src - 0.5 * d, src + 0.5 * d, radius_m, k);
            let vector = jwmu * d * d * p(zm, zn);
            let scalar = inv_jwe
                * (p(zm + 0.5 * d, zn + 0.5 * d)
                    - p(zm + 0.5 * d, zn - 0.5 * d)
                    - p(zm - 0.5 * d, zn + 0.5 * d)
                    + p(zm - 0.5 * d, zn - 0.5 * d));
            a[mi * n + ni] = vector + scalar;
        }
    }
    let feed = n / 2;
    let mut b = vec![C::new(0.0, 0.0); n];
    b[feed] = C::new(1.0, 0.0);
    // Dense complex Gaussian elimination with partial pivoting.
    for col in 0..n {
        let piv = (col..n)
            .max_by(|&i, &j| a[i * n + col].norm().total_cmp(&a[j * n + col].norm()))
            .unwrap_or(col);
        require(a[piv * n + col].norm() > 0.0, "singular moment matrix")?;
        if piv != col {
            for c in 0..n {
                a.swap(piv * n + c, col * n + c);
            }
            b.swap(piv, col);
        }
        let inv = 1.0 / a[col * n + col];
        for r in col + 1..n {
            let f = a[r * n + col] * inv;
            if f.norm() == 0.0 {
                continue;
            }
            for c in col..n {
                let v = a[col * n + c];
                a[r * n + c] -= f * v;
            }
            let v = b[col];
            b[r] -= f * v;
        }
    }
    let mut i = vec![C::new(0.0, 0.0); n];
    for r in (0..n).rev() {
        let mut sum = b[r];
        for c in r + 1..n {
            sum -= a[r * n + c] * i[c];
        }
        i[r] = sum / a[r * n + r];
    }
    require(
        i[feed].norm() > 0.0 && i.iter().all(|c| c.re.is_finite() && c.im.is_finite()),
        "moment solution is not finite",
    )?;
    Ok((C::new(1.0, 0.0) / i[feed], zc, i))
}

impl WireAntennaCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.segments >= 5 && self.segments <= 401 && self.segments % 2 == 1,
            "use an odd segment count from 5 to 401",
        )?;
        require(
            self.sweep.points <= 401,
            "wire sweep is limited to 401 points",
        )?;
        require(
            self.reference_impedance_ohm > 0.0,
            "reference impedance must be positive",
        )?;
        let (dipole_length, scale) = match self.ground {
            WireGround::FreeSpace => (self.length_m, 1.0),
            WireGround::PecPlane => (2.0 * self.length_m, 0.5),
        };
        let f = self.sweep.values()?;
        let z0 = self.reference_impedance_ohm;
        let (mut r, mut x, mut rl, mut vswr) = (vec![], vec![], vec![], vec![]);
        let mut centre = None;
        for (k, &freq) in f.iter().enumerate() {
            let (z, zc, current) =
                harrington_dipole(dipole_length, self.radius_m, self.segments, freq)?;
            let z = z * scale;
            let gamma = ((z - z0) / (z + z0)).norm();
            r.push(z.re);
            x.push(z.im);
            rl.push(-20.0 * gamma.max(1e-300).log10());
            vswr.push(if gamma < 1.0 - 1e-12 {
                (1.0 + gamma) / (1.0 - gamma)
            } else {
                1e12
            });
            if k == f.len() / 2 {
                centre = Some((zc, current));
            }
        }
        csv(
            out,
            m,
            "impedance.csv",
            &[
                ("frequency_Hz", &f),
                ("resistance_ohm", &r),
                ("reactance_ohm", &x),
                ("return_loss_dB", &rl),
                ("vswr", &vswr),
            ],
        )?;
        if let Some((zc, current)) = centre {
            let mag: Vec<f64> = current.iter().map(|i| i.norm()).collect();
            let ph: Vec<f64> = current.iter().map(|i| i.arg().to_degrees()).collect();
            csv(
                out,
                m,
                "current_distribution.csv",
                &[
                    ("z_m", &zc),
                    ("current_magnitude_A", &mag),
                    ("current_phase_deg", &ph),
                ],
            )?;
        }
        // Resonance: first sign change of reactance from negative to positive.
        if let Some(k) = (1..f.len()).find(|&k| x[k - 1] < 0.0 && x[k] >= 0.0) {
            let t = -x[k - 1] / (x[k] - x[k - 1]);
            metric(m, "resonant_frequency_Hz", f[k - 1] + t * (f[k] - f[k - 1]));
            metric(
                m,
                "resonant_resistance_ohm",
                r[k - 1] + t * (r[k] - r[k - 1]),
            );
        }
        if let Some((k, best)) = rl.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)) {
            metric(m, "best_return_loss_dB", *best);
            metric(m, "best_match_frequency_Hz", f[k]);
        }
        Ok(())
    }
}

// ---------------------------------------------------------------- 3-D FDTD

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FdtdWall {
    #[default]
    Pec,
    Pmc,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum FdtdExcitation {
    /// Initial Gaussian Ez blob in the centre (`width_fraction` of the smallest dimension).
    GaussianEz { amplitude: f64, width_fraction: f64 },
    /// Soft Ez Gaussian pulse at a cell: `amplitude · exp(-((t - delay)/width)²)`.
    PointPulse {
        cell: [usize; 3],
        amplitude: f64,
        delay_s: f64,
        width_s: f64,
    },
}

fn courant() -> f64 {
    0.95
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Maxwell3dCase {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    #[serde(default = "courant")]
    pub courant: f64,
    #[serde(default = "one")]
    pub epsilon_r: f64,
    #[serde(default)]
    pub boundary: FdtdWall,
    pub excitation: FdtdExcitation,
    pub steps: usize,
    #[serde(default)]
    pub probes: Vec<[usize; 3]>,
}

fn one() -> f64 {
    1.0
}

impl Maxwell3dCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.steps >= 1 && self.steps <= 20_000,
            "FDTD needs 1..=20000 steps",
        )?;
        require(self.probes.len() <= 16, "at most 16 probes")?;
        let mut s = Maxwell3D::new(Maxwell3DConfig {
            nx: self.nx,
            ny: self.ny,
            nz: self.nz,
            dx: self.dx,
            dy: self.dy,
            dz: self.dz,
            courant: self.courant,
            epsilon_r: self.epsilon_r,
            boundary: match self.boundary {
                FdtdWall::Pec => Boundary3D::Pec,
                FdtdWall::Pmc => Boundary3D::Pmc,
            },
            ..Default::default()
        })?;
        let index = |c: &[usize; 3]| -> Result<usize> {
            require(
                c[0] < self.nx && c[1] < self.ny && c[2] < self.nz,
                "FDTD cell index is outside the grid",
            )?;
            Ok((c[2] * self.ny + c[1]) * self.nx + c[0])
        };
        let probes: Vec<usize> = self.probes.iter().map(index).collect::<Result<_>>()?;
        if let FdtdExcitation::GaussianEz {
            amplitude,
            width_fraction,
        } = self.excitation
        {
            s.initialize_gaussian_ez(amplitude, width_fraction)?;
        }
        if let FdtdExcitation::PointPulse { cell, .. } = &self.excitation {
            index(cell)?;
        }
        let dt = s.dt();
        let mut time = vec![0.0];
        let mut energy = vec![s.energy()];
        let mut probe_series: Vec<Vec<f64>> = probes.iter().map(|&p| vec![s.ez()[p]]).collect();
        for n in 1..=self.steps {
            if let FdtdExcitation::PointPulse {
                cell,
                amplitude,
                delay_s,
                width_s,
            } = self.excitation
            {
                require(width_s > 0.0, "pulse width must be positive")?;
                let t = (n as f64 - 0.5) * dt;
                let a = ((t - delay_s) / width_s).powi(2);
                s.add_soft_ez_source(cell[0], cell[1], cell[2], amplitude * (-a).exp())?;
            }
            s.step(1);
            time.push(n as f64 * dt);
            energy.push(s.energy());
            for (series, &p) in probe_series.iter_mut().zip(&probes) {
                series.push(s.ez()[p]);
            }
        }
        require(
            energy.iter().all(|e| e.is_finite()),
            "FDTD field became non-finite",
        )?;
        csv(
            out,
            m,
            "energy.csv",
            &[("time_s", &time), ("field_energy_J", &energy)],
        )?;
        if !probes.is_empty() {
            let names: Vec<String> = (0..probes.len())
                .map(|i| format!("Ez_probe_{i}_V_per_m"))
                .collect();
            let mut table: Vec<(&str, &[f64])> = vec![("time_s", &time)];
            table.extend(
                names
                    .iter()
                    .zip(&probe_series)
                    .map(|(n, v)| (n.as_str(), v.as_slice())),
            );
            csv(out, m, "probes.csv", &table)?;
        }
        let mesh = VtuMesh::box_cells(
            [self.nx, self.ny, self.nz],
            Vec3::new(
                self.nx as f64 * self.dx,
                self.ny as f64 * self.dy,
                self.nz as f64 * self.dz,
            ),
            0.0,
        )?;
        let e: Vec<Vec3> = (0..s.cell_count())
            .map(|i| Vec3::new(s.ex()[i], s.ey()[i], s.ez()[i]))
            .collect();
        let h: Vec<Vec3> = (0..s.cell_count())
            .map(|i| Vec3::new(s.hx()[i], s.hy()[i], s.hz()[i]))
            .collect();
        crate::case::write_field(
            out,
            m,
            "fields.vtu",
            &mesh,
            vec![],
            vec![
                DataArray::vector("E_V_per_m", e),
                DataArray::vector("H_A_per_m", h),
            ],
        )?;
        metric(m, "time_step_s", dt);
        metric(m, "time_s", self.steps as f64 * dt);
        metric(m, "final_field_energy_J", s.energy());
        metric(
            m,
            "peak_field_energy_J",
            energy.iter().copied().fold(0.0, f64::max),
        );
        Ok(())
    }
}

// ---------------------------------------------------------------- Nedelec cavity

/// PEC box cavity on tetrahedra, driven by an impressed edge current between the mesh nodes
/// nearest to `port_a_m` and `port_b_m`; reports the port impedance `V/I` over frequency.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdgeCavityCase {
    pub width_m: f64,
    pub height_m: f64,
    pub depth_m: f64,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    #[serde(default = "one")]
    pub relative_permittivity: f64,
    #[serde(default)]
    pub conductivity_s_per_m: f64,
    pub port_a_m: [f64; 3],
    pub port_b_m: [f64; 3],
    pub sweep: FrequencySweep,
}

impl EdgeCavityCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.nx * self.ny * self.nz <= 512,
            "edge cavity is limited to 512 boxes (3072 tetrahedra)",
        )?;
        require(
            self.sweep.points <= 201,
            "edge cavity sweep is limited to 201 points",
        )?;
        let mesh = make_box_tet_mesh(
            self.nx,
            self.ny,
            self.nz,
            self.width_m,
            self.height_m,
            self.depth_m,
        )?;
        let nearest = |p: [f64; 3]| {
            mesh.nodes
                .iter()
                .enumerate()
                .map(|(i, n)| {
                    (
                        i,
                        (n.x - p[0]).powi(2) + (n.y - p[1]).powi(2) + (n.z - p[2]).powi(2),
                    )
                })
                .min_by(|a, b| a.1.total_cmp(&b.1))
                .map(|(i, _)| i)
        };
        let (a, b) = (nearest(self.port_a_m), nearest(self.port_b_m));
        let (Some(a), Some(b)) = (a, b) else {
            return Err(SolverError::InvalidInput("empty cavity mesh".into()));
        };
        require(a != b, "port ends snap to the same mesh node")?;
        let f = self.sweep.values()?;
        let port = [EdgeCurrentPort3D {
            node_a: a,
            node_b: b,
            current_a: C::new(1.0, 0.0),
        }];
        let (mut r, mut x, mut mag) = (vec![], vec![], vec![]);
        let mut last = None;
        let mut converged = true;
        for &freq in &f {
            let cfg = EdgeMaxwell3DConfig {
                frequency_hz: freq,
                relative_permittivity: self.relative_permittivity,
                conductivity_s_per_m: self.conductivity_s_per_m,
                physical_conductivity_term: true,
                ..Default::default()
            };
            let res = solve_driven_edge_maxwell_3d(&mesh, &cfg, None, &port)?;
            converged &= res.linear_result.converged;
            let key = [a.min(b), a.max(b)];
            let k = res.edges.iter().position(|e| *e == key).ok_or_else(|| {
                SolverError::InvalidInput("port nodes are not joined by a mesh edge".into())
            })?;
            require(!res.boundary_edge[k], "port edge lies on the PEC wall")?;
            let sign = if a < b { 1.0 } else { -1.0 };
            let z = -res.edge_voltage_v[k] * sign;
            r.push(z.re);
            x.push(z.im);
            mag.push(z.norm());
            last = Some((cfg, res));
        }
        csv(
            out,
            m,
            "port_impedance.csv",
            &[
                ("frequency_Hz", &f),
                ("resistance_ohm", &r),
                ("reactance_ohm", &x),
                ("magnitude_ohm", &mag),
            ],
        )?;
        if let Some((cfg, res)) = last {
            let e: Vec<Vec3> = res
                .electric_centroid_v_per_m
                .iter()
                .map(|v| Vec3::new(v[0].norm(), v[1].norm(), v[2].norm()))
                .collect();
            crate::case::write_field(
                out,
                m,
                "electric_field.vtu",
                &VtuMesh::from_tetrahedra(&mesh)?,
                vec![],
                vec![DataArray::vector("E_amplitude_V_per_m", e)],
            )?;
            let q = resonator_quality_3d(&mesh, &cfg, &res, 0.0, 1.0)?;
            metric(m, "electric_energy_J", q.electric_energy_j);
            metric(m, "magnetic_energy_J", q.magnetic_energy_j);
        }
        if let Some((k, peak)) = mag.iter().enumerate().max_by(|p, q| p.1.total_cmp(q.1)) {
            metric(m, "peak_impedance_ohm", *peak);
            metric(m, "peak_impedance_frequency_Hz", f[k]);
        }
        // A lossless cavity port's reactance has a pole at each parallel resonance, where it
        // jumps from + to −; 1/X crosses zero there, so interpolate that crossing.
        if let Some(k) = (1..x.len()).find(|&k| x[k - 1] > 0.0 && x[k] < 0.0) {
            let (y0, y1) = (1.0 / x[k - 1], 1.0 / x[k]);
            metric(
                m,
                "first_resonance_frequency_Hz",
                f[k - 1] + y0 / (y0 - y1) * (f[k] - f[k - 1]),
            );
        }
        metric(m, "tetrahedra", mesh.element_count() as f64);
        m.linear_solver_converged = Some(converged);
        Ok(())
    }
}

// ---------------------------------------------------------------- PDN IR drop

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PdnSource {
    pub node: [usize; 2],
    pub voltage_v: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PdnLoad {
    pub node: [usize; 2],
    pub current_a: f64,
}

/// Resistive power-delivery mesh: nodes on an `nx × ny` grid with VRM sources and sinks.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PdnIrDropCase {
    pub nx: usize,
    pub ny: usize,
    pub pitch_m: f64,
    pub segment_resistance_ohm: f64,
    pub sources: Vec<PdnSource>,
    #[serde(default)]
    pub loads: Vec<PdnLoad>,
    /// Spread evenly over every node in addition to `loads`.
    #[serde(default)]
    pub uniform_load_a: f64,
}

impl PdnIrDropCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            self.nx >= 2 && self.ny >= 2 && self.nx * self.ny <= 16_384,
            "PDN grid needs 2..=16384 nodes per side product",
        )?;
        require(
            !self.sources.is_empty() && self.sources.len() <= 64 && self.loads.len() <= 4096,
            "PDN needs 1..=64 sources and at most 4096 loads",
        )?;
        let id = |n: [usize; 2]| -> Result<usize> {
            require(
                n[0] < self.nx && n[1] < self.ny,
                "PDN node is outside the grid",
            )?;
            Ok(n[1] * self.nx + n[0])
        };
        let count = self.nx * self.ny;
        let mut load = vec![self.uniform_load_a / count as f64; count];
        for l in &self.loads {
            load[id(l.node)?] += l.current_a;
        }
        let sources = self
            .sources
            .iter()
            .map(|s| {
                Ok(PdnVoltageSource {
                    node: id(s.node)?,
                    voltage_v: s.voltage_v,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let cfg = PdnGridConfig {
            nx: self.nx,
            ny: self.ny,
            horizontal_resistance_ohm: self.segment_resistance_ohm,
            vertical_resistance_ohm: self.segment_resistance_ohm,
            max_iterations: 200_000,
            relative_tolerance: 1e-10,
        };
        let r = solve_pdn_ir_drop_grid(&cfg, &load, &sources)?;
        let vmax = self
            .sources
            .iter()
            .map(|s| s.voltage_v)
            .fold(f64::MIN, f64::max);
        let drop: Vec<f64> = r.voltage_v.iter().map(|v| vmax - v).collect();
        crate::case::write_field(
            out,
            m,
            "pdn_voltage.vtu",
            &VtuMesh::rectangle_cells(self.nx - 1, self.ny - 1, self.pitch_m, self.pitch_m)?,
            vec![
                DataArray::scalar("voltage_V", r.voltage_v.clone()),
                DataArray::scalar("ir_drop_V", drop),
                DataArray::scalar("load_A", load),
            ],
            vec![],
        )?;
        metric(m, "minimum_voltage_V", r.minimum_voltage_v);
        metric(m, "maximum_ir_drop_V", r.maximum_drop_v);
        metric(m, "gauss_seidel_iterations", r.iterations as f64);
        m.linear_solver_converged = Some(r.converged);
        Ok(())
    }
}

// ---------------------------------------------------------------- galvanic corrosion

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReactionSpec {
    pub name: String,
    pub equilibrium_potential_v: f64,
    pub exchange_current_density_a_m2: f64,
    #[serde(default = "one")]
    pub area_m2: f64,
    #[serde(default = "one")]
    pub electrons: f64,
    #[serde(default = "half")]
    pub alpha_anodic: f64,
    #[serde(default = "half")]
    pub alpha_cathodic: f64,
}

fn half() -> f64 {
    0.5
}
fn room() -> f64 {
    298.15
}

/// Butler–Volmer mixed (corrosion) potential of coupled electrode reactions and their
/// polarization curves.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GalvanicCorrosionCase {
    pub reactions: Vec<ReactionSpec>,
    #[serde(default = "room")]
    pub temperature_k: f64,
    pub polarization: VoltageSweep,
}

impl GalvanicCorrosionCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        require(
            !self.reactions.is_empty() && self.reactions.len() <= 16,
            "use 1..=16 reactions",
        )?;
        points_ok(self.polarization.points, "polarization sweep")?;
        let reactions: Vec<ElectrodeReaction> = self
            .reactions
            .iter()
            .map(|r| ElectrodeReaction {
                name: r.name.clone(),
                equilibrium_potential: r.equilibrium_potential_v,
                exchange_current_density: r.exchange_current_density_a_m2,
                area: r.area_m2,
                electrons: r.electrons,
                alpha_anodic: r.alpha_anodic,
                alpha_cathodic: r.alpha_cathodic,
            })
            .collect();
        let mp = solve_mixed_potential(&reactions, self.temperature_k, 0.0, 200, 1e-12)?;
        let e = linspace(
            self.polarization.start_v,
            self.polarization.stop_v,
            self.polarization.points,
        );
        let mut per: Vec<Vec<f64>> = vec![Vec::with_capacity(e.len()); reactions.len()];
        let mut total = Vec::with_capacity(e.len());
        for &v in &e {
            let mut sum = 0.0;
            for (k, r) in reactions.iter().enumerate() {
                let i = r.area
                    * butler_volmer_current_density(
                        r.exchange_current_density,
                        v - r.equilibrium_potential,
                        self.temperature_k,
                        r.electrons,
                        r.alpha_anodic,
                        r.alpha_cathodic,
                    )?;
                per[k].push(i);
                sum += i;
            }
            total.push(sum);
        }
        let names: Vec<String> = reactions.iter().map(|r| format!("{}_A", r.name)).collect();
        let mut table: Vec<(&str, &[f64])> = vec![("potential_V", &e), ("total_current_A", &total)];
        table.extend(
            names
                .iter()
                .zip(&per)
                .map(|(n, v)| (n.as_str(), v.as_slice())),
        );
        csv(out, m, "polarization.csv", &table)?;
        json(
            out,
            m,
            "mixed_potential.json",
            &json!({"corrosion_potential_V": mp.potential, "net_current_A": mp.total_current, "reaction_current_A": mp.reaction_current, "reaction_current_density_A_m2": mp.reaction_current_density, "iterations": mp.iterations, "converged": mp.converged}),
        )?;
        metric(m, "corrosion_potential_V", mp.potential);
        let anodic: f64 = mp.reaction_current.iter().filter(|i| **i > 0.0).sum();
        metric(m, "corrosion_current_A", anodic);
        m.linear_solver_converged = Some(mp.converged);
        Ok(())
    }
}

// ---------------------------------------------------------------- power network

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BusKind {
    Slack,
    #[default]
    Pq,
}

/// Per-unit bus. Specified powers are injections: loads are negative.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BusSpec {
    pub name: String,
    #[serde(default)]
    pub kind: BusKind,
    #[serde(default = "one")]
    pub voltage_pu: f64,
    #[serde(default)]
    pub angle_rad: f64,
    #[serde(default)]
    pub active_pu: f64,
    #[serde(default)]
    pub reactive_pu: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LineSpec {
    pub from: usize,
    pub to: usize,
    pub resistance_pu: f64,
    pub reactance_pu: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GeneratorSpec {
    pub minimum: f64,
    pub maximum: f64,
    pub linear_cost: f64,
    pub quadratic_cost: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DispatchSpec {
    pub generators: Vec<GeneratorSpec>,
    pub demand: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MachineSpec {
    pub inertia: f64,
    pub damping: f64,
    pub mechanical_power: f64,
    #[serde(default)]
    pub initial_angle_rad: f64,
    #[serde(default)]
    pub initial_speed: f64,
}

/// `None` (`null`) for an endpoint ties it to the infinite bus.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TieSpec {
    pub from: Option<usize>,
    pub to: Option<usize>,
    pub synchronizing_power: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SwingSpec {
    pub machines: Vec<MachineSpec>,
    pub ties: Vec<TieSpec>,
    pub dt_s: f64,
    pub steps: usize,
}

/// Per-unit AC network: Newton power flow at the base load, a load-scaling continuation (PV
/// curve) until Newton stops converging, and optional economic dispatch and multi-machine
/// swing transients.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PowerNetworkCase {
    pub buses: Vec<BusSpec>,
    pub lines: Vec<LineSpec>,
    #[serde(default = "one")]
    pub max_load_factor: f64,
    #[serde(default = "continuation_points")]
    pub continuation_points: usize,
    #[serde(default)]
    pub dispatch: Option<DispatchSpec>,
    #[serde(default)]
    pub swing: Option<SwingSpec>,
}

fn continuation_points() -> usize {
    21
}

impl PowerNetworkCase {
    pub(crate) fn run(&self, out: &Path, m: &mut RunManifest) -> Result<()> {
        use crate::power::*;
        require(
            self.buses.len() >= 2
                && self.buses.len() <= 200
                && !self.lines.is_empty()
                && self.lines.len() <= 1000,
            "use 2..=200 buses and 1..=1000 lines",
        )?;
        require(
            self.buses
                .iter()
                .filter(|b| b.kind == BusKind::Slack)
                .count()
                == 1,
            "the network needs exactly one slack bus",
        )?;
        require(
            (2..=MAX_SWEEP_POINTS).contains(&self.continuation_points)
                && self.max_load_factor >= 1.0,
            "use 2.. continuation points and a maximum load factor of at least 1",
        )?;
        let buses: Vec<PowerBus> = self
            .buses
            .iter()
            .map(|b| PowerBus {
                kind: if b.kind == BusKind::Slack {
                    PowerBusType::Slack
                } else {
                    PowerBusType::Pq
                },
                voltage: b.voltage_pu,
                angle: b.angle_rad,
                specified_active: b.active_pu,
                specified_reactive: b.reactive_pu,
            })
            .collect();
        let lines = self
            .lines
            .iter()
            .map(|l| PowerLine {
                from: l.from,
                to: l.to,
                resistance: l.resistance_pu,
                reactance: l.reactance_pu,
            })
            .collect();
        let network = MultiBusPowerFlow::new(buses, lines)?;
        let pq: Vec<usize> = (0..self.buses.len())
            .filter(|&i| self.buses[i].kind == BusKind::Pq)
            .collect();
        let unpack = |state: &[f64]| -> (Vec<f64>, Vec<f64>) {
            let mut v: Vec<f64> = self.buses.iter().map(|b| b.voltage_pu).collect();
            let mut a: Vec<f64> = self.buses.iter().map(|b| b.angle_rad).collect();
            for (k, &bus) in pq.iter().enumerate() {
                v[bus] = state[k];
                a[bus] = state[pq.len() + k];
            }
            (v, a)
        };
        // The finite-difference Newton can stall a little above its 1e-8 absolute target;
        // a residual below 1e-6 p.u. is accepted as a solved operating point.
        let solved = |r: &MultiBusPowerFlowResult| {
            r.solve.converged || (r.solve.residual_rms.is_finite() && r.solve.residual_rms < 1e-6)
        };
        let base = network.solve(1.0, None)?;
        require(
            base.state.len() == 2 * pq.len(),
            "unexpected power-flow state size",
        )?;
        let (v, a) = unpack(&base.state);
        let index: Vec<f64> = (0..self.buses.len()).map(|i| i as f64).collect();
        let deg: Vec<f64> = a.iter().map(|x| x.to_degrees()).collect();
        csv(
            out,
            m,
            "power_flow.csv",
            &[("bus", &index), ("voltage_pu", &v), ("angle_deg", &deg)],
        )?;
        json(
            out,
            m,
            "buses.json",
            &self
                .buses
                .iter()
                .map(|b| b.name.clone())
                .collect::<Vec<_>>(),
        )?;
        metric(m, "base_newton_iterations", base.solve.iterations as f64);
        metric(m, "base_residual_rms", base.solve.residual_rms);
        metric(
            m,
            "minimum_voltage_pu",
            v.iter().copied().fold(f64::MAX, f64::min),
        );
        let mut converged = solved(&base);
        metric(
            m,
            "base_newton_strictly_converged",
            if base.solve.converged { 1.0 } else { 0.0 },
        );
        // Step the load factor ourselves so the looser acceptance above applies to every point.
        let mut path = Vec::new();
        let mut previous: Option<Vec<f64>> = Some(base.state.clone());
        for i in 0..self.continuation_points {
            let factor = 1.0
                + (self.max_load_factor - 1.0) * i as f64 / (self.continuation_points - 1) as f64;
            let Ok(result) = network.solve(factor, previous.as_deref()) else {
                break;
            };
            let accepted = solved(&result);
            previous = Some(result.state.clone());
            path.push(ContinuationPoint {
                load_factor: factor,
                result,
            });
            if !accepted {
                break;
            }
        }
        let ok: Vec<&ContinuationPoint> = path.iter().take_while(|p| solved(&p.result)).collect();
        let lf: Vec<f64> = ok.iter().map(|p| p.load_factor).collect();
        let mut columns: Vec<(String, Vec<f64>)> = Vec::new();
        for &bus in &pq {
            columns.push((
                format!("V_{}_pu", self.buses[bus].name),
                ok.iter().map(|p| unpack(&p.result.state).0[bus]).collect(),
            ));
        }
        let mut table: Vec<(&str, &[f64])> = vec![("load_factor", &lf)];
        table.extend(columns.iter().map(|(n, v)| (n.as_str(), v.as_slice())));
        csv(out, m, "pv_curve.csv", &table)?;
        metric(
            m,
            "last_converged_load_factor",
            lf.last().copied().unwrap_or(0.0),
        );
        metric(
            m,
            "continuation_reached_maximum",
            if ok.len() == self.continuation_points {
                1.0
            } else {
                0.0
            },
        );
        if let Some(d) = &self.dispatch {
            require(
                !d.generators.is_empty() && d.generators.len() <= 100,
                "use 1..=100 generators",
            )?;
            let gens: Vec<QuadraticDispatchGenerator> = d
                .generators
                .iter()
                .map(|g| QuadraticDispatchGenerator {
                    minimum: g.minimum,
                    maximum: g.maximum,
                    linear_cost: g.linear_cost,
                    quadratic_cost: g.quadratic_cost,
                })
                .collect();
            let r = quadratic_dispatch(&gens, d.demand)?;
            json(
                out,
                m,
                "dispatch.json",
                &json!({"generation": r.generation, "served": r.served, "objective": r.objective, "feasible": r.feasible}),
            )?;
            metric(m, "dispatch_objective", r.objective);
            metric(m, "dispatch_served", r.served);
            converged &= r.feasible;
        }
        if let Some(sw) = &self.swing {
            require(
                !sw.machines.is_empty()
                    && sw.machines.len() <= 50
                    && sw.steps >= 1
                    && sw.steps <= 20_000,
                "use 1..=50 machines and 1..=20000 steps",
            )?;
            let machines: Vec<SwingBusParameters> = sw
                .machines
                .iter()
                .map(|g| SwingBusParameters {
                    inertia: g.inertia,
                    damping: g.damping,
                    mechanical_power: g.mechanical_power,
                    initial_angle: g.initial_angle_rad,
                    initial_speed: g.initial_speed,
                })
                .collect();
            let ties: Vec<SwingLine> = sw
                .ties
                .iter()
                .map(|t| SwingLine {
                    from: t.from,
                    to: t.to,
                    synchronizing_power: t.synchronizing_power,
                })
                .collect();
            let traj = integrate_swing_network(&machines, &ties, sw.dt_s, sw.steps)?;
            require(
                traj.states.iter().all(|s| s.len() == 2 * machines.len()),
                "unexpected swing state size",
            )?;
            let mut cols: Vec<(String, Vec<f64>)> = Vec::new();
            for i in 0..machines.len() {
                cols.push((
                    format!("angle_{i}_rad"),
                    traj.states.iter().map(|s| s[2 * i]).collect(),
                ));
                cols.push((
                    format!("speed_{i}"),
                    traj.states.iter().map(|s| s[2 * i + 1]).collect(),
                ));
            }
            let mut table: Vec<(&str, &[f64])> = vec![("time", &traj.times)];
            table.extend(cols.iter().map(|(n, v)| (n.as_str(), v.as_slice())));
            csv(out, m, "swing.csv", &table)?;
            if let Some(last) = traj.states.last() {
                metric(m, "swing_final_angle_0_rad", last[0]);
                metric(
                    m,
                    "swing_max_final_speed",
                    (0..machines.len())
                        .map(|i| last[2 * i + 1].abs())
                        .fold(0.0, f64::max),
                );
            }
            converged &= traj.converged;
        }
        m.linear_solver_converged = Some(converged);
        Ok(())
    }
}
