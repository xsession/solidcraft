//! Distribution planning: weighted contiguous repartitioning and migration segments, ported
//! from upstream `cfd/distributed/repartition` (MIT, cfd_solvers contributors).
mod cartesian;
use crate::{require, Result};
pub use cartesian::*;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ContiguousRepartitionPlan {
    pub boundaries: Vec<usize>,
    pub loads: Vec<f64>,
    pub average_load: f64,
    pub imbalance: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MigrationSegment {
    pub source_part: usize,
    pub destination_part: usize,
    pub global_begin: usize,
    pub global_end: usize,
}

/// Prefix-sum cut points nearest to equal load (upstream default `minimum_cells_per_part = 1`).
pub fn plan_contiguous_repartition(
    weights: &[f64],
    parts: usize,
    minimum_cells_per_part: usize,
) -> Result<ContiguousRepartitionPlan> {
    require(
        !weights.is_empty()
            && parts > 0
            && minimum_cells_per_part > 0
            && parts
                .checked_mul(minimum_cells_per_part)
                .is_some_and(|m| weights.len() >= m),
        "invalid contiguous repartition dimensions",
    )?;
    let mut prefix = vec![0.0; weights.len() + 1];
    for (i, w) in weights.iter().enumerate() {
        require(*w >= 0.0 && w.is_finite(), "invalid repartition weight")?;
        prefix[i + 1] = prefix[i] + w;
    }
    let total = prefix[weights.len()];
    let mut plan = ContiguousRepartitionPlan {
        boundaries: vec![0; parts + 1],
        ..Default::default()
    };
    plan.boundaries[parts] = weights.len();
    for part in 1..parts {
        let minimum = plan.boundaries[part - 1] + minimum_cells_per_part;
        let maximum = weights.len() - (parts - part) * minimum_cells_per_part;
        if total.is_nan() || total <= 0.0 {
            let remaining = weights.len() - plan.boundaries[part - 1];
            plan.boundaries[part] = (plan.boundaries[part - 1] + remaining / (parts - part + 1))
                .clamp(minimum, maximum);
            continue;
        }
        let target = total * part as f64 / parts as f64;
        let mut selected = minimum + prefix[minimum..=maximum].partition_point(|v| *v < target);
        if selected > minimum
            && (prefix[selected - 1] - target).abs()
                <= (prefix.get(selected).copied().unwrap_or(f64::INFINITY) - target).abs()
        {
            selected -= 1;
        }
        plan.boundaries[part] = selected.clamp(minimum, maximum);
    }
    plan.loads = (0..parts)
        .map(|p| prefix[plan.boundaries[p + 1]] - prefix[plan.boundaries[p]])
        .collect();
    let maximum_load = plan.loads.iter().copied().fold(0.0, f64::max);
    plan.average_load = total / parts as f64;
    plan.imbalance = if plan.average_load > 0.0 {
        maximum_load / plan.average_load
    } else {
        1.0
    };
    Ok(plan)
}

fn validate_boundaries(b: &[usize]) -> Result<()> {
    require(
        b.len() >= 2 && b[0] == 0,
        "partition boundaries must begin at zero",
    )?;
    require(
        b.windows(2).all(|w| w[1] >= w[0]),
        "partition boundaries must be monotone",
    )
}

/// Ranges that change owner between two partitions of the same global range.
pub fn build_migration_segments(
    old_boundaries: &[usize],
    new_boundaries: &[usize],
) -> Result<Vec<MigrationSegment>> {
    validate_boundaries(old_boundaries)?;
    validate_boundaries(new_boundaries)?;
    require(
        old_boundaries.last() == new_boundaries.last(),
        "old and new partitions must cover the same global range",
    )?;
    let mut out = Vec::new();
    for s in 0..old_boundaries.len() - 1 {
        for d in 0..new_boundaries.len() - 1 {
            let begin = old_boundaries[s].max(new_boundaries[d]);
            let end = old_boundaries[s + 1].min(new_boundaries[d + 1]);
            if begin < end && s != d {
                out.push(MigrationSegment {
                    source_part: s,
                    destination_part: d,
                    global_begin: begin,
                    global_end: end,
                });
            }
        }
    }
    Ok(out)
}

/// Upstream defaults `trigger = 1.15`, `minimum_relative_improvement = 0.05`.
pub fn repartition_recommended(
    current_imbalance: f64,
    proposed_imbalance: f64,
    trigger_imbalance: f64,
    minimum_relative_improvement: f64,
) -> Result<bool> {
    require(
        current_imbalance >= 1.0
            && proposed_imbalance >= 1.0
            && trigger_imbalance >= 1.0
            && (0.0..1.0).contains(&minimum_relative_improvement),
        "invalid repartition decision controls",
    )?;
    if current_imbalance < trigger_imbalance || proposed_imbalance >= current_imbalance {
        return Ok(false);
    }
    Ok(
        (current_imbalance - proposed_imbalance) / current_imbalance
            >= minimum_relative_improvement,
    )
}
