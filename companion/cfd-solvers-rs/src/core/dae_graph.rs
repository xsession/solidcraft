//! Equation-based hybrid DAE graph with claimed residual rows and acausal connections, ported
//! from upstream `cfd/core/hybrid_dae_graph` (MIT, cfd_solvers contributors).
use super::{integrate_implicit_euler, DaeEvent, DaeSystem, DaeTrajectory, MatrixFreeNewtonConfig};
use crate::{require, Result};

#[derive(Clone, Debug, PartialEq)]
pub struct HybridDaeVariable {
    pub name: String,
    pub initial_value: f64,
    pub differential: bool,
}

/// `residual(t, state, derivative, local_rows)` for one component.
pub type HybridDaeResidual = Box<dyn Fn(f64, &[f64], &[f64], &mut [f64]) -> Result<()>>;

struct Component {
    #[allow(dead_code)]
    name: String,
    rows: Vec<usize>,
    residual: HybridDaeResidual,
}

/// Global variables; each component owns explicit residual rows, so every variable row must be
/// claimed exactly once before integration.
#[derive(Default)]
pub struct HybridDaeGraph {
    variables: Vec<HybridDaeVariable>,
    row_claimed: Vec<bool>,
    components: Vec<Component>,
}

impl HybridDaeGraph {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn add_variable(
        &mut self,
        name: &str,
        initial_value: f64,
        differential: bool,
    ) -> Result<usize> {
        require(
            !name.is_empty() && initial_value.is_finite(),
            "invalid hybrid DAE variable",
        )?;
        self.variables.push(HybridDaeVariable {
            name: name.into(),
            initial_value,
            differential,
        });
        self.row_claimed.push(false);
        Ok(self.variables.len() - 1)
    }
    pub fn add_component(
        &mut self,
        name: &str,
        residual_rows: &[usize],
        residual: HybridDaeResidual,
    ) -> Result<()> {
        require(
            !name.is_empty() && !residual_rows.is_empty(),
            "invalid hybrid DAE component",
        )?;
        let mut seen = std::collections::BTreeSet::new();
        require(
            residual_rows
                .iter()
                .all(|&r| r < self.variables.len() && !self.row_claimed[r] && seen.insert(r)),
            "hybrid DAE residual row is already claimed",
        )?;
        for &r in residual_rows {
            self.row_claimed[r] = true;
        }
        self.components.push(Component {
            name: name.into(),
            rows: residual_rows.to_vec(),
            residual,
        });
        Ok(())
    }
    /// Row `left - ratio * right - offset = 0`.
    pub fn add_connection(
        &mut self,
        name: &str,
        residual_row: usize,
        left_variable: usize,
        right_variable: usize,
        ratio: f64,
        offset: f64,
    ) -> Result<()> {
        require(
            ratio != 0.0
                && ratio.is_finite()
                && offset.is_finite()
                && left_variable < self.variables.len()
                && right_variable < self.variables.len(),
            "invalid hybrid DAE connection",
        )?;
        self.add_component(
            name,
            &[residual_row],
            Box::new(move |_, state, _, out| {
                out[0] = state[left_variable] - ratio * state[right_variable] - offset;
                Ok(())
            }),
        )
    }
    pub fn variables(&self) -> &[HybridDaeVariable] {
        &self.variables
    }
    pub fn size(&self) -> usize {
        self.variables.len()
    }
    pub fn initial_state(&self) -> Vec<f64> {
        self.variables.iter().map(|v| v.initial_value).collect()
    }
    /// Implicit-Euler integration of the assembled residual.
    pub fn integrate(
        &self,
        initial_time: f64,
        step: f64,
        steps: usize,
        events: &[DaeEvent<'_>],
        nonlinear: MatrixFreeNewtonConfig,
    ) -> Result<DaeTrajectory> {
        require(
            !self.variables.is_empty()
                && !self.components.is_empty()
                && steps > 0
                && step > 0.0
                && initial_time.is_finite()
                && step.is_finite(),
            "invalid hybrid DAE graph integration",
        )?;
        require(
            self.row_claimed.iter().all(|c| *c),
            "hybrid DAE graph has unclaimed residual rows",
        )?;
        let residual = |t: f64, state: &[f64], derivative: &[f64], out: &mut [f64]| -> Result<()> {
            out.fill(0.0);
            for c in &self.components {
                let mut local = vec![0.0; c.rows.len()];
                (c.residual)(t, state, derivative, &mut local)?;
                for (row, v) in c.rows.iter().zip(&local) {
                    out[*row] = *v;
                }
            }
            Ok(())
        };
        let system = DaeSystem {
            size: self.variables.len(),
            residual: &residual,
        };
        integrate_implicit_euler(
            &system,
            initial_time,
            step,
            &self.initial_state(),
            steps,
            events,
            nonlinear,
        )
    }
}
