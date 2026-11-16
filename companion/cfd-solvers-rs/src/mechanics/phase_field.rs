//! Quasi-static 1-D phase-field fracture (alternate minimization with a history field) and the
//! operator-based variant.
#![allow(clippy::needless_range_loop)]
use crate::{require, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseFieldFracture1DConfig {
    pub cells: usize,
    pub length: f64,
    pub area: f64,
    pub youngs_modulus: f64,
    pub fracture_energy: f64,
    pub length_scale: f64,
    pub residual_stiffness: f64,
    pub damage_iterations: usize,
    pub damage_tolerance: f64,
    /// `None` means no notch (upstream `size_t::max()`).
    pub notch_cell: Option<usize>,
    pub initial_notch_damage: f64,
}

impl Default for PhaseFieldFracture1DConfig {
    fn default() -> Self {
        Self {
            cells: 64,
            length: 1.0,
            area: 1.0,
            youngs_modulus: 100.0,
            fracture_energy: 0.01,
            length_scale: 0.05,
            residual_stiffness: 1e-8,
            damage_iterations: 100,
            damage_tolerance: 1e-10,
            notch_cell: None,
            initial_notch_damage: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PhaseFieldFracture1DDiagnostics {
    pub reaction_force: f64,
    pub elastic_energy: f64,
    pub fracture_energy: f64,
    pub total_energy: f64,
    pub maximum_damage: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PhaseFieldFracture1D {
    config: PhaseFieldFracture1DConfig,
    displacement: Vec<f64>,
    damage: Vec<f64>,
    history: Vec<f64>,
}

impl PhaseFieldFracture1D {
    pub fn new(config: PhaseFieldFracture1DConfig) -> Result<Self> {
        let c = &config;
        require(
            c.cells >= 2
                && c.cells <= 1 << 26
                && c.length > 0.0
                && c.area > 0.0
                && c.youngs_modulus > 0.0
                && c.fracture_energy > 0.0
                && c.length_scale > 0.0
                && c.residual_stiffness >= 0.0
                && c.residual_stiffness < 1.0
                && c.damage_iterations > 0
                && c.damage_tolerance > 0.0
                && c.initial_notch_damage >= 0.0
                && c.initial_notch_damage <= 1.0
                && c.notch_cell.is_none_or(|n| n < c.cells),
            "invalid phase-field fracture configuration",
        )?;
        let mut damage = vec![0.0; c.cells];
        if let Some(n) = c.notch_cell {
            damage[n] = c.initial_notch_damage;
        }
        Ok(Self {
            config,
            displacement: vec![0.0; c.cells + 1],
            damage,
            history: vec![0.0; c.cells],
        })
    }
    pub fn displacement(&self) -> &[f64] {
        &self.displacement
    }
    pub fn damage(&self) -> &[f64] {
        &self.damage
    }
    pub fn history(&self) -> &[f64] {
        &self.history
    }
    pub fn step(
        &mut self,
        prescribed_end_displacement: f64,
    ) -> Result<PhaseFieldFracture1DDiagnostics> {
        require(
            prescribed_end_displacement.is_finite(),
            "phase-field displacement must be finite",
        )?;
        let c = self.config;
        let dx = c.length / c.cells as f64;
        let strain = prescribed_end_displacement / c.length;
        let te = 0.5 * c.youngs_modulus * strain.max(0.0) * strain.max(0.0);
        self.history.iter_mut().for_each(|h| *h = h.max(te));
        let mut next = vec![0.0; self.damage.len()];
        let gc = c.fracture_energy * c.length_scale / (dx * dx);
        let n = self.damage.len();
        for _ in 0..c.damage_iterations {
            let mut max_change: f64 = 0.0;
            for cell in 0..n {
                let (hl, hr) = (cell > 0, cell + 1 < n);
                let sum = (if hl {
                    self.damage[cell - 1]
                } else {
                    self.damage[cell]
                }) + (if hr {
                    self.damage[cell + 1]
                } else {
                    self.damage[cell]
                });
                let count = hl as u8 as f64 + hr as u8 as f64;
                let diag =
                    2.0 * self.history[cell] + c.fracture_energy / c.length_scale + count * gc;
                let cand = (2.0 * self.history[cell] + gc * sum) / diag;
                next[cell] = self.damage[cell].max(cand).clamp(0.0, 1.0);
                max_change = max_change.max((next[cell] - self.damage[cell]).abs());
            }
            std::mem::swap(&mut self.damage, &mut next);
            if max_change <= c.damage_tolerance {
                break;
            }
        }
        for node in 0..self.displacement.len() {
            self.displacement[node] = prescribed_end_displacement * node as f64 / c.cells as f64;
        }
        Ok(self.diagnostics_at(strain, dx))
    }
    pub fn diagnostics(&self) -> PhaseFieldFracture1DDiagnostics {
        let dx = self.config.length / self.config.cells as f64;
        self.diagnostics_at(
            self.displacement[self.displacement.len() - 1] / self.config.length,
            dx,
        )
    }
    fn diagnostics_at(&self, strain: f64, dx: f64) -> PhaseFieldFracture1DDiagnostics {
        let c = &self.config;
        let n = self.damage.len();
        let (mut el, mut fr, mut avg, mut mx): (f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0);
        for cell in 0..n {
            let d = self.damage[cell];
            let g = (1.0 - d) * (1.0 - d) + c.residual_stiffness;
            let te = 0.5 * c.youngs_modulus * strain.max(0.0) * strain.max(0.0);
            el += g * te * dx * c.area;
            let l = if cell == 0 { d } else { self.damage[cell - 1] };
            let r = if cell + 1 == n {
                d
            } else {
                self.damage[cell + 1]
            };
            let grad = (r - l) / (2.0 * dx);
            fr += c.fracture_energy
                * (d * d / (2.0 * c.length_scale) + 0.5 * c.length_scale * grad * grad)
                * dx
                * c.area;
            avg += g * c.youngs_modulus * strain;
            mx = mx.max(d);
        }
        avg /= n as f64;
        PhaseFieldFracture1DDiagnostics {
            reaction_force: avg * c.area,
            elastic_energy: el,
            fracture_energy: fr,
            total_energy: el + fr,
            maximum_damage: mx,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledPhaseFieldFracture1DConfig {
    pub cells: usize,
    pub length: f64,
    pub area: f64,
    pub youngs_modulus: f64,
    pub fracture_energy: f64,
    pub length_scale: f64,
    pub residual_stiffness: f64,
}

impl Default for CoupledPhaseFieldFracture1DConfig {
    fn default() -> Self {
        Self {
            cells: 64,
            length: 1.0,
            area: 1.0,
            youngs_modulus: 100.0,
            fracture_energy: 0.01,
            length_scale: 0.05,
            residual_stiffness: 1e-8,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CoupledPhaseFieldFracture1DResult {
    pub residual_l2: f64,
    pub reaction_force: f64,
    pub elastic_energy: f64,
    pub fracture_energy: f64,
    pub total_energy: f64,
    pub maximum_damage: f64,
}

/// Matrix-free coupled AT2 displacement/damage operator; state layout `[u_0..u_n, d_0..d_(n-1)]`
/// with a frozen history field between loads.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledPhaseFieldFracture1D {
    config: CoupledPhaseFieldFracture1DConfig,
    history: Vec<f64>,
}

impl CoupledPhaseFieldFracture1D {
    pub fn new(config: CoupledPhaseFieldFracture1DConfig) -> Result<Self> {
        let c = &config;
        require(
            c.cells >= 2
                && c.cells <= 1 << 26
                && c.length > 0.0
                && c.area > 0.0
                && c.youngs_modulus > 0.0
                && c.fracture_energy > 0.0
                && c.length_scale > 0.0
                && c.residual_stiffness >= 0.0
                && c.residual_stiffness < 1.0
                && [
                    c.length,
                    c.area,
                    c.youngs_modulus,
                    c.fracture_energy,
                    c.length_scale,
                    c.residual_stiffness,
                ]
                .iter()
                .all(|v| v.is_finite()),
            "invalid coupled phase-field fracture configuration",
        )?;
        Ok(Self {
            config,
            history: vec![0.0; c.cells],
        })
    }
    pub fn displacement_dofs(&self) -> usize {
        self.config.cells + 1
    }
    pub fn damage_dofs(&self) -> usize {
        self.config.cells
    }
    pub fn size(&self) -> usize {
        self.displacement_dofs() + self.damage_dofs()
    }
    pub fn history(&self) -> &[f64] {
        &self.history
    }
    pub fn zero_state(&self) -> Vec<f64> {
        vec![0.0; self.size()]
    }
    fn dx(&self) -> f64 {
        self.config.length / self.config.cells as f64
    }
    fn degradation(&self, d: f64) -> f64 {
        (1.0 - d) * (1.0 - d) + self.config.residual_stiffness
    }
    fn validate_state(&self, s: &[f64]) -> Result<()> {
        require(s.len() == self.size(), "phase-field state size mismatch")?;
        require(
            s.iter().all(|v| v.is_finite()),
            "phase-field state must be finite",
        )?;
        require(
            s[self.displacement_dofs()..]
                .iter()
                .all(|d| (0.0..=1.0).contains(d)),
            "phase-field damage must be in [0,1]",
        )
    }
    fn validate_vector(&self, v: &[f64], name: &str) -> Result<()> {
        require(
            v.len() == self.size() && v.iter().all(|x| x.is_finite()),
            name,
        )
    }
    /// Takes the maximum of the stored and current tensile energy per cell.
    pub fn update_history(&mut self, state: &[f64]) -> Result<()> {
        self.validate_state(state)?;
        let dx = self.dx();
        for cell in 0..self.config.cells {
            let e = ((state[cell + 1] - state[cell]) / dx).max(0.0);
            self.history[cell] = self.history[cell].max(0.5 * self.config.youngs_modulus * e * e);
        }
        Ok(())
    }
    pub fn evaluate(
        &self,
        state: &[f64],
        prescribed_end_displacement: f64,
        residual: &mut [f64],
    ) -> Result<CoupledPhaseFieldFracture1DResult> {
        self.validate_state(state)?;
        self.validate_vector(residual, "phase-field residual")?;
        require(
            prescribed_end_displacement.is_finite(),
            "phase-field prescribed displacement must be finite",
        )?;
        residual.fill(0.0);
        let c = &self.config;
        let dx = self.dx();
        let off = self.displacement_dofs();
        let (mut el, mut fr, mut reaction, mut mx): (f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0);
        for cell in 0..c.cells {
            let strain = (state[cell + 1] - state[cell]) / dx;
            let d = state[off + cell];
            let g = self.degradation(d);
            let force = c.area * g * c.youngs_modulus * strain;
            residual[cell] -= force;
            residual[cell + 1] += force;
            if cell == 0 {
                reaction = force;
            }
            el += 0.5 * g * c.youngs_modulus * strain * strain * dx * c.area;
            let l = if cell == 0 { d } else { state[off + cell - 1] };
            let r = if cell + 1 == c.cells {
                d
            } else {
                state[off + cell + 1]
            };
            let lap = (r - 2.0 * d + l) / (dx * dx);
            residual[off + cell] = 2.0 * self.history[cell] * (1.0 - d)
                - c.fracture_energy * d / c.length_scale
                + c.fracture_energy * c.length_scale * lap;
            fr += c.fracture_energy
                * (d * d / (2.0 * c.length_scale)
                    + 0.5 * c.length_scale * (r - l) * (r - l) / (4.0 * dx * dx))
                * dx
                * c.area;
            mx = mx.max(d);
        }
        residual[0] = state[0];
        residual[off - 1] = state[off - 1] - prescribed_end_displacement;
        let norm = residual.iter().map(|v| v * v).sum::<f64>().sqrt();
        Ok(CoupledPhaseFieldFracture1DResult {
            residual_l2: norm,
            reaction_force: reaction,
            elastic_energy: el,
            fracture_energy: fr,
            total_energy: el + fr,
            maximum_damage: mx,
        })
    }
    pub fn tangent_action(
        &self,
        state: &[f64],
        direction: &[f64],
        prescribed_end_displacement: f64,
        output: &mut [f64],
    ) -> Result<()> {
        self.validate_state(state)?;
        self.validate_vector(direction, "phase-field tangent direction")?;
        self.validate_vector(output, "phase-field tangent output")?;
        require(
            prescribed_end_displacement.is_finite(),
            "phase-field prescribed displacement must be finite",
        )?;
        output.fill(0.0);
        let c = &self.config;
        let dx = self.dx();
        let off = self.displacement_dofs();
        for cell in 0..c.cells {
            let strain = (state[cell + 1] - state[cell]) / dx;
            let d = state[off + cell];
            let dd = direction[off + cell];
            let sd = (direction[cell + 1] - direction[cell]) / dx;
            let g = self.degradation(d);
            let gd = -2.0 * (1.0 - d) * dd;
            let fd = c.area * c.youngs_modulus * (g * sd + gd * strain);
            output[cell] -= fd;
            output[cell + 1] += fd;
            let l = if cell == 0 {
                dd
            } else {
                direction[off + cell - 1]
            };
            let r = if cell + 1 == c.cells {
                dd
            } else {
                direction[off + cell + 1]
            };
            let lap = (r - 2.0 * dd + l) / (dx * dx);
            output[off + cell] = -2.0 * self.history[cell] * dd
                - c.fracture_energy * dd / c.length_scale
                + c.fracture_energy * c.length_scale * lap;
        }
        output[0] = direction[0];
        output[off - 1] = direction[off - 1];
        Ok(())
    }
}
