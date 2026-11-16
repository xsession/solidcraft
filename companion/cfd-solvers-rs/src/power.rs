//! AC two-bus and multi-bus PQ load flow with continuation, swing-network transients and
//! quadratic economic dispatch, ported from upstream `cfd/power/{ac_power_flow, multibus_power}`
//! (MIT, cfd_solvers contributors).
use crate::{
    core::{
        solve_matrix_free_newton, DaeTrajectory, HybridDaeGraph, MatrixFreeNewtonConfig,
        MatrixFreeNewtonResult, MatrixFreeNonlinearSystem,
    },
    require, Result,
};
use std::rc::Rc;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoBusPowerFlowConfig {
    pub line_resistance: f64,
    pub line_reactance: f64,
    pub slack_voltage: f64,
    pub specified_active_power: f64,
    pub specified_reactive_power: f64,
}

impl Default for TwoBusPowerFlowConfig {
    fn default() -> Self {
        Self {
            line_resistance: 0.01,
            line_reactance: 0.1,
            slack_voltage: 1.0,
            specified_active_power: -0.5,
            specified_reactive_power: -0.2,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TwoBusPowerFlowResult {
    pub voltage: f64,
    pub angle: f64,
    pub solve: MatrixFreeNewtonResult,
}

fn admittance(r: f64, x: f64) -> (f64, f64) {
    let d = r * r + x * x;
    (r / d, -x / d)
}

/// Polar PQ load flow for one load bus fed from a slack bus, with an analytic Jacobian action.
pub fn solve_two_bus_power_flow(
    config: &TwoBusPowerFlowConfig,
    initial: Option<[f64; 2]>,
) -> Result<TwoBusPowerFlowResult> {
    let c = *config;
    require(
        c.slack_voltage > 0.0
            && c.line_resistance >= 0.0
            && c.line_reactance.is_finite()
            && !(c.line_resistance == 0.0 && c.line_reactance == 0.0),
        "invalid two-bus power-flow configuration",
    )?;
    let mut state = initial.unwrap_or([c.slack_voltage, 0.0]);
    let (g, b) = admittance(c.line_resistance, c.line_reactance);
    let vs = c.slack_voltage;
    let residual = |x: &[f64], out: &mut [f64]| -> Result<()> {
        let (v, a) = (x[0], x[1]);
        let (cs, sn) = (a.cos(), a.sin());
        out[0] = v * v * g - v * vs * (g * cs + b * sn) - c.specified_active_power;
        out[1] = -v * v * b - v * vs * (g * sn - b * cs) - c.specified_reactive_power;
        Ok(())
    };
    let jacobian = |x: &[f64], d: &[f64], out: &mut [f64]| -> Result<()> {
        let (v, a) = (x[0], x[1]);
        let (cs, sn) = (a.cos(), a.sin());
        let dp_dv = 2.0 * v * g - vs * (g * cs + b * sn);
        let dp_da = vs * v * (g * sn - b * cs);
        let dq_dv = -2.0 * v * b - vs * (g * sn - b * cs);
        let dq_da = -vs * v * (g * cs + b * sn);
        out[0] = dp_dv * d[0] + dp_da * d[1];
        out[1] = dq_dv * d[0] + dq_da * d[1];
        Ok(())
    };
    let system = MatrixFreeNonlinearSystem {
        size: 2,
        residual: &residual,
        jacobian_action: Some(&jacobian),
        transpose_jacobian_action: None,
    };
    let nonlinear = MatrixFreeNewtonConfig {
        max_iterations: 30,
        linear_max_iterations: 12,
        linear_restart: 2,
        absolute_tolerance: 1e-12,
        relative_tolerance: 1e-11,
        linear_relative_tolerance: 1e-11,
        ..Default::default()
    };
    let solve = solve_matrix_free_newton(&mut state, &system, nonlinear)?;
    Ok(TwoBusPowerFlowResult {
        voltage: state[0],
        angle: state[1],
        solve,
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PowerBusType {
    Slack,
    Pq,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerBus {
    pub kind: PowerBusType,
    pub voltage: f64,
    pub angle: f64,
    pub specified_active: f64,
    pub specified_reactive: f64,
}

impl Default for PowerBus {
    fn default() -> Self {
        Self {
            kind: PowerBusType::Pq,
            voltage: 1.0,
            angle: 0.0,
            specified_active: 0.0,
            specified_reactive: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerLine {
    pub from: usize,
    pub to: usize,
    pub resistance: f64,
    pub reactance: f64,
}

impl Default for PowerLine {
    fn default() -> Self {
        Self {
            from: 0,
            to: 0,
            resistance: 0.01,
            reactance: 0.1,
        }
    }
}

#[derive(Clone, Debug)]
pub struct MultiBusPowerFlowResult {
    /// PQ voltages followed by PQ angles.
    pub state: Vec<f64>,
    pub pq_buses: usize,
    pub solve: MatrixFreeNewtonResult,
}

#[derive(Clone, Debug)]
pub struct ContinuationPoint {
    pub load_factor: f64,
    pub result: MultiBusPowerFlowResult,
}

/// Dense-admittance polar PQ network with one slack bus; the Jacobian action is the
/// finite-difference fallback, as upstream.
#[derive(Clone, Debug)]
pub struct MultiBusPowerFlow {
    buses: Vec<PowerBus>,
    conductance: Vec<f64>,
    susceptance: Vec<f64>,
    bus_to_pq: Vec<Option<usize>>,
    pq_buses: Vec<usize>,
}

impl MultiBusPowerFlow {
    pub fn new(buses: Vec<PowerBus>, lines: Vec<PowerLine>) -> Result<Self> {
        let n = buses.len();
        require(n >= 2, "multi-bus network needs two buses")?;
        require(
            buses.iter().all(|b| {
                b.voltage > 0.0
                    && b.voltage.is_finite()
                    && b.angle.is_finite()
                    && b.specified_active.is_finite()
                    && b.specified_reactive.is_finite()
            }),
            "invalid power bus",
        )?;
        let mut bus_to_pq = vec![None; n];
        let mut pq_buses = Vec::new();
        let mut slack = 0;
        for (i, b) in buses.iter().enumerate() {
            if b.kind == PowerBusType::Slack {
                slack += 1;
            } else {
                bus_to_pq[i] = Some(pq_buses.len());
                pq_buses.push(i);
            }
        }
        require(
            slack == 1 && !pq_buses.is_empty(),
            "network needs one slack and one PQ bus",
        )?;
        let (mut gm, mut bm) = (vec![0.0; n * n], vec![0.0; n * n]);
        for l in &lines {
            require(
                l.from < n
                    && l.to < n
                    && l.from != l.to
                    && l.resistance >= 0.0
                    && l.resistance.is_finite()
                    && l.reactance.is_finite()
                    && !(l.resistance == 0.0 && l.reactance == 0.0),
                "invalid power line",
            )?;
            let (g, b) = admittance(l.resistance, l.reactance);
            for (r, c, s) in [
                (l.from, l.from, 1.0),
                (l.to, l.to, 1.0),
                (l.from, l.to, -1.0),
                (l.to, l.from, -1.0),
            ] {
                gm[r * n + c] += s * g;
                bm[r * n + c] += s * b;
            }
        }
        Ok(Self {
            buses,
            conductance: gm,
            susceptance: bm,
            bus_to_pq,
            pq_buses,
        })
    }
    pub fn solve(
        &self,
        load_factor: f64,
        initial: Option<&[f64]>,
    ) -> Result<MultiBusPowerFlowResult> {
        require(
            load_factor >= 0.0 && load_factor.is_finite(),
            "invalid power load factor",
        )?;
        let npq = self.pq_buses.len();
        let unknowns = 2 * npq;
        let mut state = vec![1.0; unknowns];
        for (i, &bus) in self.pq_buses.iter().enumerate() {
            state[i] = self.buses[bus].voltage;
            state[npq + i] = self.buses[bus].angle;
        }
        if let Some(init) = initial {
            require(
                init.len() == unknowns,
                "multi-bus initial state size mismatch",
            )?;
            state.copy_from_slice(init);
        }
        let n = self.buses.len();
        let residual = |x: &[f64], out: &mut [f64]| -> Result<()> {
            for (i, &bus) in self.pq_buses.iter().enumerate() {
                let (vi, ai) = (x[i], x[npq + i]);
                let (mut p, mut q) = (0.0, 0.0);
                for j in 0..n {
                    let (vj, aj) = match self.bus_to_pq[j] {
                        None => (self.buses[j].voltage, self.buses[j].angle),
                        Some(k) => (x[k], x[npq + k]),
                    };
                    let d = ai - aj;
                    let (g, b) = (self.conductance[bus * n + j], self.susceptance[bus * n + j]);
                    p += vi * vj * (g * d.cos() + b * d.sin());
                    q += vi * vj * (g * d.sin() - b * d.cos());
                }
                out[2 * i] = p - load_factor * self.buses[bus].specified_active;
                out[2 * i + 1] = q - load_factor * self.buses[bus].specified_reactive;
            }
            Ok(())
        };
        let system = MatrixFreeNonlinearSystem {
            size: unknowns,
            residual: &residual,
            jacobian_action: None,
            transpose_jacobian_action: None,
        };
        let nonlinear = MatrixFreeNewtonConfig {
            max_iterations: 40,
            linear_max_iterations: 16usize.max(4 * unknowns),
            linear_restart: unknowns.min(8),
            absolute_tolerance: 1e-8,
            relative_tolerance: 1e-9,
            linear_relative_tolerance: 1e-7,
            finite_difference_step: 1e-7,
            ..Default::default()
        };
        let solve = solve_matrix_free_newton(&mut state, &system, nonlinear)?;
        Ok(MultiBusPowerFlowResult {
            state,
            pq_buses: npq,
            solve,
        })
    }
    /// Natural-parameter continuation warm-started from the previous point; stops after the first
    /// non-converged point.
    pub fn continuation(
        &self,
        first_factor: f64,
        last_factor: f64,
        points: usize,
    ) -> Result<Vec<ContinuationPoint>> {
        require(
            points >= 2
                && first_factor.is_finite()
                && last_factor.is_finite()
                && first_factor >= 0.0
                && last_factor >= 0.0,
            "invalid power continuation",
        )?;
        let mut path = Vec::with_capacity(points);
        let mut previous: Option<Vec<f64>> = None;
        for i in 0..points {
            let factor =
                first_factor + i as f64 / (points - 1) as f64 * (last_factor - first_factor);
            let result = self.solve(factor, previous.as_deref())?;
            previous = Some(result.state.clone());
            let converged = result.solve.converged;
            path.push(ContinuationPoint {
                load_factor: factor,
                result,
            });
            if !converged {
                break;
            }
        }
        Ok(path)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticDispatchGenerator {
    pub minimum: f64,
    pub maximum: f64,
    pub linear_cost: f64,
    pub quadratic_cost: f64,
}

impl Default for QuadraticDispatchGenerator {
    fn default() -> Self {
        Self {
            minimum: 0.0,
            maximum: 1.0,
            linear_cost: 1.0,
            quadratic_cost: 0.01,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct DispatchResult {
    pub generation: Vec<f64>,
    pub served: f64,
    pub objective: f64,
    pub feasible: bool,
}

/// Greedy marginal-cost loading from the minimum outputs: each pass fills the cheapest unit that
/// still has headroom (up to `generators.len()` passes), as upstream.
pub fn quadratic_dispatch(
    generators: &[QuadraticDispatchGenerator],
    demand: f64,
) -> Result<DispatchResult> {
    require(
        !generators.is_empty() && demand >= 0.0 && demand.is_finite(),
        "invalid quadratic dispatch problem",
    )?;
    require(
        generators.iter().all(|g| {
            g.minimum >= 0.0
                && g.maximum >= g.minimum
                && g.quadratic_cost >= 0.0
                && g.linear_cost.is_finite()
        }),
        "invalid dispatch generator",
    )?;
    let minimum_total: f64 = generators.iter().map(|g| g.minimum).sum();
    let maximum_total: f64 = generators.iter().map(|g| g.maximum).sum();
    if demand < minimum_total || demand > maximum_total {
        return Ok(DispatchResult {
            generation: vec![0.0; generators.len()],
            served: 0.0,
            objective: 0.0,
            feasible: false,
        });
    }
    let mut generation: Vec<f64> = generators.iter().map(|g| g.minimum).collect();
    let mut remaining = demand - minimum_total;
    let mut pass = 0;
    while pass < generators.len() && remaining > 0.0 {
        let mut best = 0;
        let mut best_cost = f64::INFINITY;
        for (i, g) in generators.iter().enumerate() {
            let marginal = g.linear_cost + 2.0 * g.quadratic_cost * generation[i];
            if generation[i] < g.maximum && marginal < best_cost {
                best = i;
                best_cost = marginal;
            }
        }
        let increment = remaining.min(generators[best].maximum - generation[best]);
        generation[best] += increment;
        remaining -= increment;
        pass += 1;
    }
    let served = generation.iter().sum();
    let objective = generators
        .iter()
        .zip(&generation)
        .map(|(g, p)| g.linear_cost * p + g.quadratic_cost * p * p)
        .sum();
    Ok(DispatchResult {
        generation,
        served,
        objective,
        feasible: remaining <= 1e-12,
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SwingBusParameters {
    pub inertia: f64,
    pub damping: f64,
    pub mechanical_power: f64,
    pub initial_angle: f64,
    pub initial_speed: f64,
}

impl Default for SwingBusParameters {
    fn default() -> Self {
        Self {
            inertia: 1.0,
            damping: 0.1,
            mechanical_power: 0.0,
            initial_angle: 0.0,
            initial_speed: 0.0,
        }
    }
}

/// `None` endpoint is the fixed reference (slack) angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SwingLine {
    pub from: Option<usize>,
    pub to: Option<usize>,
    pub synchronizing_power: f64,
}

/// Multi-machine swing equations `M w' + D w = P_m - sum P_s sin(d_i - d_j)` assembled through
/// the hybrid DAE graph and integrated with implicit Euler.
pub fn integrate_swing_network(
    buses: &[SwingBusParameters],
    lines: &[SwingLine],
    step: f64,
    steps: usize,
) -> Result<DaeTrajectory> {
    require(
        !buses.is_empty() && !lines.is_empty() && step > 0.0 && steps > 0 && step.is_finite(),
        "invalid swing-network integration",
    )?;
    require(
        buses.iter().all(|b| {
            b.inertia > 0.0
                && b.damping >= 0.0
                && b.inertia.is_finite()
                && b.damping.is_finite()
                && b.mechanical_power.is_finite()
                && b.initial_angle.is_finite()
                && b.initial_speed.is_finite()
        }),
        "invalid swing-bus parameters",
    )?;
    require(
        lines.iter().all(|l| {
            l.from.is_none_or(|i| i < buses.len())
                && l.to.is_none_or(|i| i < buses.len())
                && l.from != l.to
                && l.synchronizing_power > 0.0
                && l.synchronizing_power.is_finite()
        }),
        "invalid swing-network line",
    )?;
    let mut graph = HybridDaeGraph::new();
    let mut angles = Vec::with_capacity(buses.len());
    let mut speeds = Vec::with_capacity(buses.len());
    for (i, b) in buses.iter().enumerate() {
        angles.push(graph.add_variable(&format!("angle_{i}"), b.initial_angle, true)?);
        speeds.push(graph.add_variable(&format!("speed_{i}"), b.initial_speed, true)?);
    }
    let (angles, speeds) = (Rc::new(angles), Rc::new(speeds));
    let lines = Rc::new(lines.to_vec());
    for (i, bus) in buses.iter().copied().enumerate() {
        let (a, s, l) = (angles.clone(), speeds.clone(), lines.clone());
        graph.add_component(
            &format!("swing_{i}"),
            &[angles[i], speeds[i]],
            Box::new(move |_, state, derivative, out| {
                out[0] = derivative[a[i]] - state[s[i]];
                let mut electrical = 0.0;
                for line in l.iter() {
                    if line.from != Some(i) && line.to != Some(i) {
                        continue;
                    }
                    let other = if line.from == Some(i) {
                        line.to
                    } else {
                        line.from
                    };
                    let other_angle = other.map_or(0.0, |o| state[a[o]]);
                    electrical += line.synchronizing_power * (state[a[i]] - other_angle).sin();
                }
                out[1] = bus.inertia * derivative[s[i]] + bus.damping * state[s[i]]
                    - bus.mechanical_power
                    + electrical;
                Ok(())
            }),
        )?;
    }
    let nonlinear = MatrixFreeNewtonConfig {
        max_iterations: 30,
        linear_max_iterations: 32,
        linear_restart: 4,
        absolute_tolerance: 1e-6,
        relative_tolerance: 1e-7,
        linear_relative_tolerance: 1e-7,
        ..Default::default()
    };
    graph.integrate(0.0, step, steps, &[], nonlinear)
}
