//! 1-D continuum references: total-Lagrangian cubic-hardening bar (matrix-free Newton), bond-based
//! peridynamic bar and penalty impact.
use crate::{
    core::{
        solve_matrix_free_newton, MatrixFreeNewtonConfig, MatrixFreeNewtonResult,
        MatrixFreeNonlinearSystem,
    },
    require, Result,
};
use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearBarConfig {
    pub elements: usize,
    pub length: f64,
    pub area: f64,
    pub youngs_modulus: f64,
    pub cubic_modulus: f64,
}

impl Default for NonlinearBarConfig {
    fn default() -> Self {
        Self {
            elements: 8,
            length: 1.0,
            area: 1.0,
            youngs_modulus: 100.0,
            cubic_modulus: 100.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct NonlinearBarResult {
    pub displacement: Vec<f64>,
    pub solve: MatrixFreeNewtonResult,
}

/// Fixed-left bar with an end force; `initial` may be empty.
pub fn solve_nonlinear_bar(
    config: &NonlinearBarConfig,
    applied_force: f64,
    initial: &[f64],
) -> Result<NonlinearBarResult> {
    let c = *config;
    require(
        c.elements >= 2
            && c.length > 0.0
            && c.area > 0.0
            && c.youngs_modulus > 0.0
            && c.cubic_modulus >= 0.0
            && applied_force.is_finite(),
        "invalid nonlinear bar configuration",
    )?;
    let spacing = c.length / c.elements as f64;
    let mut unknown = vec![0.0; c.elements];
    if !initial.is_empty() {
        require(
            initial.len() == unknown.len(),
            "nonlinear bar initial state size mismatch",
        )?;
        unknown.copy_from_slice(initial);
    }
    let stress = |e: f64| c.youngs_modulus * e + c.cubic_modulus * e * e * e;
    let residual = move |state: &[f64], out: &mut [f64]| -> Result<()> {
        let mut nodes = vec![0.0; c.elements + 1];
        nodes[1..].copy_from_slice(state);
        let es: Vec<f64> = (0..c.elements)
            .map(|e| stress((nodes[e + 1] - nodes[e]) / spacing) * c.area)
            .collect();
        for node in 1..c.elements {
            out[node - 1] = es[node - 1] - es[node];
        }
        out[c.elements - 1] = es[c.elements - 1] - applied_force;
        Ok(())
    };
    let system = MatrixFreeNonlinearSystem {
        size: unknown.len(),
        residual: &residual,
        jacobian_action: None,
        transpose_jacobian_action: None,
    };
    let nl = MatrixFreeNewtonConfig {
        max_iterations: 30,
        linear_max_iterations: 20usize.max(c.elements * 4),
        linear_restart: 4usize.max(20usize.min(c.elements)),
        absolute_tolerance: 1e-11,
        relative_tolerance: 1e-10,
        linear_relative_tolerance: 1e-10,
        ..Default::default()
    };
    let solve = solve_matrix_free_newton(&mut unknown, &system, nl)?;
    Ok(NonlinearBarResult {
        displacement: unknown,
        solve,
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeridynamicBarConfig {
    pub points: usize,
    pub length: f64,
    pub horizon: f64,
    pub stiffness: f64,
    pub critical_stretch: f64,
}

impl Default for PeridynamicBarConfig {
    fn default() -> Self {
        Self {
            points: 32,
            length: 1.0,
            horizon: 0.1,
            stiffness: 100.0,
            critical_stretch: 0.08,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PeridynamicStepDiagnostics {
    pub broken_bonds: usize,
    pub kinetic_energy: f64,
    pub potential_energy: f64,
}

/// Bond-based peridynamic bar with irreversible bond breaking (unit masses, symplectic Euler).
#[derive(Clone, Debug, PartialEq)]
pub struct PeridynamicBar1D {
    config: PeridynamicBarConfig,
    positions: Vec<f64>,
    displacement: Vec<f64>,
    velocity: Vec<f64>,
    broken: BTreeSet<(usize, usize)>,
}

impl PeridynamicBar1D {
    pub fn new(config: PeridynamicBarConfig) -> Result<Self> {
        let c = &config;
        require(
            c.points >= 3
                && c.points <= 1 << 24
                && c.length > 0.0
                && c.horizon > 0.0
                && c.stiffness > 0.0
                && c.critical_stretch > 0.0,
            "invalid peridynamic bar configuration",
        )?;
        let spacing = c.length / (c.points - 1) as f64;
        Ok(Self {
            config,
            positions: (0..c.points).map(|i| spacing * i as f64).collect(),
            displacement: vec![0.0; c.points],
            velocity: vec![0.0; c.points],
            broken: BTreeSet::new(),
        })
    }
    pub fn set_displacement(&mut self, index: usize, value: f64) -> Result<()> {
        require(
            index < self.displacement.len() && value.is_finite(),
            "peridynamic displacement index",
        )?;
        self.displacement[index] = value;
        Ok(())
    }
    pub fn set_velocity(&mut self, index: usize, value: f64) -> Result<()> {
        require(
            index < self.velocity.len() && value.is_finite(),
            "peridynamic velocity index",
        )?;
        self.velocity[index] = value;
        Ok(())
    }
    pub fn displacement(&self) -> &[f64] {
        &self.displacement
    }
    pub fn velocity(&self) -> &[f64] {
        &self.velocity
    }
    /// `external_force` may be empty.
    pub fn step(&mut self, dt: f64, external_force: &[f64]) -> Result<PeridynamicStepDiagnostics> {
        let n = self.positions.len();
        require(
            dt > 0.0 && dt.is_finite() && (external_force.is_empty() || external_force.len() == n),
            "invalid peridynamic timestep",
        )?;
        let c = self.config;
        let mut force = vec![0.0; n];
        let mut potential = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                let r0 = self.positions[j] - self.positions[i];
                if r0 > c.horizon {
                    break;
                }
                if self.broken.contains(&(i, j)) {
                    continue;
                }
                let r = r0 + self.displacement[j] - self.displacement[i];
                if r.is_nan() || r <= 0.0 {
                    self.broken.insert((i, j));
                    continue;
                }
                let s = (r - r0) / r0;
                if s.abs() > c.critical_stretch {
                    self.broken.insert((i, j));
                    continue;
                }
                let f = c.stiffness * s;
                force[i] += f;
                force[j] -= f;
                potential += 0.5 * c.stiffness * s * s * r0;
            }
        }
        for i in 0..n {
            if !external_force.is_empty() {
                force[i] += external_force[i];
            }
            self.velocity[i] += dt * force[i];
            self.displacement[i] += dt * self.velocity[i];
        }
        let kinetic = self.velocity.iter().map(|v| 0.5 * v * v).sum();
        Ok(PeridynamicStepDiagnostics {
            broken_bonds: self.broken.len(),
            kinetic_energy: kinetic,
            potential_energy: potential,
        })
    }
    pub fn broken_bond_count(&self) -> usize {
        self.broken.len()
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ImpactState1D {
    pub left_position: f64,
    pub right_position: f64,
    pub left_velocity: f64,
    pub right_velocity: f64,
    pub contact_force: f64,
}

/// Two point masses with a one-sided penalty spring and approach-only damping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PenaltyImpact1D {
    left_mass: f64,
    right_mass: f64,
    stiffness: f64,
    damping: f64,
    state: ImpactState1D,
}

impl PenaltyImpact1D {
    pub fn new(
        left_mass: f64,
        right_mass: f64,
        initial_gap: f64,
        stiffness: f64,
        damping: f64,
    ) -> Result<Self> {
        require(
            left_mass > 0.0
                && right_mass > 0.0
                && initial_gap >= 0.0
                && stiffness > 0.0
                && damping >= 0.0,
            "invalid penalty impact configuration",
        )?;
        Ok(Self {
            left_mass,
            right_mass,
            stiffness,
            damping,
            state: ImpactState1D {
                right_position: initial_gap,
                ..Default::default()
            },
        })
    }
    pub fn set_velocities(&mut self, left: f64, right: f64) -> Result<()> {
        require(
            left.is_finite() && right.is_finite(),
            "invalid impact velocity",
        )?;
        self.state.left_velocity = left;
        self.state.right_velocity = right;
        Ok(())
    }
    pub fn state(&self) -> &ImpactState1D {
        &self.state
    }
    pub fn step(&mut self, dt: f64) -> Result<()> {
        require(dt > 0.0 && dt.is_finite(), "invalid impact timestep")?;
        let s = &mut self.state;
        let overlap = s.left_position - s.right_position;
        let rv = s.left_velocity - s.right_velocity;
        s.contact_force = if overlap > 0.0 {
            self.stiffness * overlap + self.damping * rv.max(0.0)
        } else {
            0.0
        };
        s.left_velocity -= dt * s.contact_force / self.left_mass;
        s.right_velocity += dt * s.contact_force / self.right_mass;
        s.left_position += dt * s.left_velocity;
        s.right_position += dt * s.right_velocity;
        Ok(())
    }
}
