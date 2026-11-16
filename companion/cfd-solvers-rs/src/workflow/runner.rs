use crate::{require, Result, SolverError};
use std::{collections::BTreeSet, path::Path};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SweepResult {
    pub parameter: f64,
    pub objective: f64,
}

pub fn parameter_sweep(values: &[f64], objective: impl Fn(f64) -> f64) -> Vec<SweepResult> {
    values
        .iter()
        .map(|&x| SweepResult {
            parameter: x,
            objective: objective(x),
        })
        .collect()
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct InverseResult {
    pub parameters: Vec<f64>,
    pub objective: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// Compass search: try `-step` then `+step` per coordinate, accept improvements, halve all steps
/// after a sweep without improvement; converged when the largest step is below `tolerance`.
/// Upstream defaults: 200 iterations, tolerance 1e-8.
pub fn inverse_coordinate_search(
    objective: impl Fn(&[f64]) -> f64,
    initial: Vec<f64>,
    step: Vec<f64>,
    max_iterations: usize,
    tolerance: f64,
) -> Result<InverseResult> {
    require(
        !initial.is_empty() && initial.len() == step.len() && max_iterations > 0 && tolerance > 0.0,
        "invalid inverse search",
    )?;
    let (mut x, mut step) = (initial, step);
    let mut best = objective(&x);
    let mut out = InverseResult {
        parameters: x.clone(),
        objective: best,
        iterations: 0,
        converged: false,
    };
    for it in 0..max_iterations {
        let mut improved = false;
        for j in 0..x.len() {
            for sign in [-1.0, 1.0] {
                let mut trial = x.clone();
                trial[j] += sign * step[j];
                let v = objective(&trial);
                if v < best {
                    x = trial;
                    best = v;
                    improved = true;
                }
            }
        }
        if !improved {
            for s in &mut step {
                *s *= 0.5;
            }
        }
        let max_step = step.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        out = InverseResult {
            parameters: x.clone(),
            objective: best,
            iterations: it + 1,
            converged: max_step < tolerance,
        };
        if out.converged {
            break;
        }
    }
    Ok(out)
}

pub type WorkflowAction = Box<dyn FnMut() -> Result<()>>;

struct Node {
    name: String,
    action: WorkflowAction,
    deps: Vec<usize>,
}

/// Ordered DAG of named actions; completed names are persisted (sorted, one per line) after every
/// action so an interrupted run resumes where it stopped.
#[derive(Default)]
pub struct RestartableWorkflow {
    nodes: Vec<Node>,
    completed: Vec<String>,
}

impl RestartableWorkflow {
    pub fn new() -> Self {
        Self::default()
    }
    /// Dependencies must refer to earlier nodes.
    pub fn add(
        &mut self,
        name: &str,
        action: WorkflowAction,
        dependencies: Vec<usize>,
    ) -> Result<usize> {
        require(!name.is_empty(), "invalid workflow node")?;
        require(
            dependencies.iter().all(|d| *d < self.nodes.len()),
            "workflow dependency must precede node",
        )?;
        self.nodes.push(Node {
            name: name.into(),
            action,
            deps: dependencies,
        });
        Ok(self.nodes.len() - 1)
    }
    pub fn run(&mut self, checkpoint: impl AsRef<Path>) -> Result<()> {
        let checkpoint = checkpoint.as_ref();
        let mut done: BTreeSet<String> = std::fs::read_to_string(checkpoint)
            .unwrap_or_default()
            .lines()
            .filter(|l| !l.is_empty())
            .map(String::from)
            .collect();
        self.completed = done.iter().cloned().collect();
        for i in 0..self.nodes.len() {
            if done.contains(&self.nodes[i].name) {
                continue;
            }
            if self.nodes[i]
                .deps
                .iter()
                .any(|d| !done.contains(&self.nodes[*d].name))
            {
                return Err(SolverError::NumericalFailure(
                    "workflow dependency incomplete".into(),
                ));
            }
            (self.nodes[i].action)()?;
            done.insert(self.nodes[i].name.clone());
            self.completed.push(self.nodes[i].name.clone());
            let text: String = done.iter().map(|n| format!("{n}\n")).collect();
            std::fs::write(checkpoint, text)?;
        }
        Ok(())
    }
    pub fn completed(&self) -> &[String] {
        &self.completed
    }
}
