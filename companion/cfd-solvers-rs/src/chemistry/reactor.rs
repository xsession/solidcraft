use super::ReactionNetwork;
use crate::{
    core::{newton_solve, CsrBuilder, NewtonConfig},
    positive, require, Result, SolverError,
};
#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ReactorConfig {
    pub initial_step: f64,
    pub minimum_step: f64,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
    pub max_steps: usize,
    pub newton_iterations: usize,
}
impl Default for ReactorConfig {
    fn default() -> Self {
        Self {
            initial_step: 1e-4,
            minimum_step: 1e-14,
            relative_tolerance: 1e-4,
            absolute_tolerance: 1e-10,
            max_steps: 100000,
            newton_iterations: 30,
        }
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct ReactorResult {
    pub time_advanced: f64,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
}
fn backward_euler(
    network: &ReactionNetwork,
    old: &[f64],
    temperature: f64,
    dt: f64,
    iterations: usize,
) -> Result<Option<Vec<f64>>> {
    let mut next = old.to_vec();
    let n = old.len();
    let residual = |x: &[f64], r: &mut [f64]| {
        if x.iter().any(|v| *v < 0.0 || !v.is_finite()) {
            r.fill(f64::INFINITY);
            return Ok(());
        }
        let source = network.source_terms(x, temperature)?;
        for i in 0..n {
            r[i] = x[i] - old[i] - dt * source[i];
        }
        Ok(())
    };
    let jacobian = |x: &[f64]| {
        let base = network.source_terms(x, temperature)?;
        let mut b = CsrBuilder::new(n, n);
        let mut perturbed = x.to_vec();
        for j in 0..n {
            let h = 1e-7 * x[j].abs().max(1.0);
            perturbed[j] = x[j] + h;
            let s = network.source_terms(&perturbed, temperature)?;
            for i in 0..n {
                b.add(
                    i,
                    j,
                    (if i == j { 1.0 } else { 0.0 }) - dt * (s[i] - base[i]) / h,
                )?;
            }
            perturbed[j] = x[j];
        }
        b.build(0.0)
    };
    let config = NewtonConfig {
        max_iterations: iterations,
        absolute_tolerance: 1e-13,
        relative_tolerance: 1e-11,
        linear_relative_tolerance: 1e-13,
        line_search_steps: 30,
        minimum_step: 1e-12,
        ..Default::default()
    };
    let result = newton_solve(&mut next, residual, jacobian, config)?;
    Ok(if result.converged { Some(next) } else { None })
}
/// On failure the caller's concentrations remain unchanged.
pub fn integrate_isothermal(
    network: &ReactionNetwork,
    c: &mut [f64],
    temperature: f64,
    duration: f64,
    config: ReactorConfig,
) -> Result<ReactorResult> {
    require(
        !c.is_empty()
            && c.len() == network.species().len()
            && duration.is_finite()
            && duration >= 0.0,
        "reactor state/duration",
    )?;
    positive(config.initial_step, "initial step")?;
    positive(config.minimum_step, "minimum step")?;
    positive(config.relative_tolerance, "relative tolerance")?;
    positive(config.absolute_tolerance, "absolute tolerance")?;
    require(
        config.minimum_step <= config.initial_step
            && config.max_steps > 0
            && config.newton_iterations > 0,
        "reactor step controls",
    )?;
    network.source_terms(c, temperature)?;
    let mut state = c.to_vec();
    let mut result = ReactorResult::default();
    let mut step = config.initial_step.min(duration);
    let mut attempt = 0;
    while result.time_advanced < duration {
        if attempt >= config.max_steps {
            return Err(SolverError::NumericalFailure(
                "reactor step budget exhausted".into(),
            ));
        }
        attempt += 1;
        step = step.min(duration - result.time_advanced);
        require(
            result.time_advanced + step > result.time_advanced,
            "reactor timestep underflow",
        )?;
        let mut candidate = None;
        let mut error: f64 = 0.0;
        if let Some(full) =
            backward_euler(network, &state, temperature, step, config.newton_iterations)?
        {
            if let Some(half) = backward_euler(
                network,
                &state,
                temperature,
                0.5 * step,
                config.newton_iterations,
            )? {
                if let Some(next) = backward_euler(
                    network,
                    &half,
                    temperature,
                    0.5 * step,
                    config.newton_iterations,
                )? {
                    for i in 0..state.len() {
                        error = error.max(
                            (full[i] - next[i]).abs()
                                / (config.absolute_tolerance
                                    + config.relative_tolerance
                                        * state[i].abs().max(next[i].abs())),
                        );
                    }
                    candidate = Some(next);
                }
            }
        }
        let accepted = if error <= 1.0 { candidate.take() } else { None };
        if let Some(accepted) = accepted {
            state = accepted;
            result.time_advanced += step;
            result.accepted_steps += 1;
            step *= (0.9 / error.max(1e-12).sqrt()).clamp(0.5, 2.0);
        } else {
            result.rejected_steps += 1;
            step *= if candidate.is_some() {
                (0.9 / error.sqrt()).clamp(0.1, 0.5)
            } else {
                0.5
            };
            if step < config.minimum_step {
                return Err(SolverError::NumericalFailure(
                    "reactor failed at minimum step".into(),
                ));
            }
        }
    }
    c.copy_from_slice(&state);
    Ok(result)
}
