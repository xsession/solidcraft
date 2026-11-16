//! Balanced series packs and a 3-D voxel conduction grid coupled to pack heat generation.
#![allow(clippy::needless_range_loop)]
use super::{invalid, BatteryStepResult, LithiumIonCellConfig, SingleParticleElectrolyteModel};
use crate::{core::matrix_free_conjugate_gradient, require, Result, SolverError};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PackControl {
    pub minimum_cell_voltage_v: f64,
    pub maximum_cell_voltage_v: f64,
    pub maximum_temperature_k: f64,
    pub balancing_current_a: f64,
    pub balancing_voltage_threshold_v: f64,
    pub balancing_minimum_voltage_v: f64,
}

impl Default for PackControl {
    fn default() -> Self {
        Self {
            minimum_cell_voltage_v: 2.5,
            maximum_cell_voltage_v: 4.3,
            maximum_temperature_k: 333.15,
            balancing_current_a: 0.0,
            balancing_voltage_threshold_v: 0.01,
            balancing_minimum_voltage_v: 3.5,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PackStepResult {
    pub cells: Vec<BatteryStepResult>,
    pub balancing_currents_a: Vec<f64>,
    pub time_s: f64,
    pub terminal_voltage_v: f64,
    pub terminal_power_w: f64,
    pub balancing_heat_w: f64,
    pub net_terminal_charge_ah: f64,
    pub terminal_energy_wh: f64,
    pub balancing_energy_wh: f64,
}

fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}

/// SPMe cells in series with passive bleed balancing and voltage/temperature cutoffs; every
/// step is transactional.
#[derive(Clone)]
pub struct SeriesBatteryPack {
    control: PackControl,
    cells: Vec<SingleParticleElectrolyteModel>,
    time_s: f64,
    charge_ah: f64,
    energy_wh: f64,
    balance_wh: f64,
}

impl SeriesBatteryPack {
    pub fn new(configs: &[LithiumIonCellConfig], control: PackControl) -> Result<Self> {
        let c = control;
        require(
            !configs.is_empty()
                && nonnegative(c.minimum_cell_voltage_v)
                && c.maximum_cell_voltage_v.is_finite()
                && c.maximum_cell_voltage_v > c.minimum_cell_voltage_v
                && c.maximum_temperature_k.is_finite()
                && c.maximum_temperature_k > 0.0
                && nonnegative(c.balancing_current_a)
                && nonnegative(c.balancing_voltage_threshold_v)
                && nonnegative(c.balancing_minimum_voltage_v),
            "invalid series-pack configuration",
        )?;
        let cells = configs
            .iter()
            .map(|cfg| SingleParticleElectrolyteModel::new(cfg.clone()))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            control,
            cells,
            time_s: 0.0,
            charge_ah: 0.0,
            energy_wh: 0.0,
            balance_wh: 0.0,
        })
    }
    fn check_limits(&self, s: &BatteryStepResult) -> Result<()> {
        let c = &self.control;
        if !s.terminal_voltage_v.is_finite()
            || !s.temperature_k.is_finite()
            || s.terminal_voltage_v < c.minimum_cell_voltage_v
            || s.terminal_voltage_v > c.maximum_cell_voltage_v
            || s.temperature_k <= 0.0
            || s.temperature_k > c.maximum_temperature_k
        {
            return Err(SolverError::NumericalFailure(
                "series-pack cell voltage or temperature cutoff".into(),
            ));
        }
        Ok(())
    }
    pub fn set_external_temperatures(&mut self, temperature_k: &[f64]) -> Result<()> {
        require(
            temperature_k.len() == self.cells.len(),
            "pack temperature count mismatch",
        )?;
        let mut trial = self.cells.clone();
        for (cell, t) in trial.iter_mut().zip(temperature_k) {
            cell.set_external_temperature(*t)?;
            self.check_limits(&cell.state(0.0)?)?;
        }
        self.cells = trial;
        Ok(())
    }
    pub fn reset(&mut self) {
        self.cells
            .iter_mut()
            .for_each(SingleParticleElectrolyteModel::reset);
        (self.time_s, self.charge_ah, self.energy_wh, self.balance_wh) = (0.0, 0.0, 0.0, 0.0);
    }
    pub fn step(&mut self, current: f64, dt: f64) -> Result<PackStepResult> {
        require(
            current.is_finite() && dt.is_finite() && dt > 0.0 && (self.time_s + dt).is_finite(),
            "invalid series-pack step",
        )?;
        let mut trial = self.cells.clone();
        let mut out = PackStepResult::default();
        let mut voltages = Vec::with_capacity(self.cells.len());
        let mut lowest = f64::INFINITY;
        for c in &self.cells {
            let s = c.state(0.0)?;
            self.check_limits(&s)?;
            voltages.push(s.terminal_voltage_v);
            lowest = lowest.min(s.terminal_voltage_v);
        }
        for (i, cell) in trial.iter_mut().enumerate() {
            let bleed = if voltages[i] >= self.control.balancing_minimum_voltage_v
                && voltages[i] - lowest > self.control.balancing_voltage_threshold_v
            {
                self.control.balancing_current_a
            } else {
                0.0
            };
            let cell_current = current + bleed;
            require(cell_current.is_finite(), "pack cell current overflow")?;
            self.check_limits(&cell.state(cell_current)?)?;
            let s = cell.step(cell_current, dt)?;
            self.check_limits(&s)?;
            out.terminal_voltage_v += s.terminal_voltage_v;
            out.balancing_heat_w += bleed * s.terminal_voltage_v;
            out.cells.push(s);
            out.balancing_currents_a.push(bleed);
        }
        out.time_s = self.time_s + dt;
        out.terminal_power_w = current * out.terminal_voltage_v;
        out.net_terminal_charge_ah = self.charge_ah + current * dt / 3600.0;
        out.terminal_energy_wh = self.energy_wh + out.terminal_power_w * dt / 3600.0;
        out.balancing_energy_wh = self.balance_wh + out.balancing_heat_w * dt / 3600.0;
        require(
            [
                out.terminal_voltage_v,
                out.balancing_heat_w,
                out.net_terminal_charge_ah,
                out.terminal_energy_wh,
                out.balancing_energy_wh,
            ]
            .iter()
            .all(|v| v.is_finite()),
            "series-pack accounting overflow",
        )?;
        self.cells = trial;
        self.time_s = out.time_s;
        self.charge_ah = out.net_terminal_charge_ah;
        self.energy_wh = out.terminal_energy_wh;
        self.balance_wh = out.balancing_energy_wh;
        Ok(out)
    }
    pub fn cells(&self) -> &[SingleParticleElectrolyteModel] {
        &self.cells
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalVoxel {
    pub volumetric_heat_capacity_j_per_m3_k: f64,
    pub conductivity_w_per_m_k: [f64; 3],
}

impl Default for ThermalVoxel {
    fn default() -> Self {
        Self {
            volumetric_heat_capacity_j_per_m3_k: 2e6,
            conductivity_w_per_m_k: [1.0; 3],
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalBoundary {
    pub convection_w_per_m2_k: f64,
    pub ambient_temperature_k: f64,
}

impl Default for ThermalBoundary {
    fn default() -> Self {
        Self {
            convection_w_per_m2_k: 0.0,
            ambient_temperature_k: 298.15,
        }
    }
}

/// Boundaries are ordered (-x, +x, -y, +y, -z, +z).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalGridConfig {
    pub cells: [usize; 3],
    pub spacing_m: [f64; 3],
    pub boundaries: [ThermalBoundary; 6],
    pub initial_temperature_k: f64,
    pub maximum_iterations: usize,
    pub relative_tolerance: f64,
}

impl Default for ThermalGridConfig {
    fn default() -> Self {
        Self {
            cells: [8; 3],
            spacing_m: [0.005; 3],
            boundaries: [ThermalBoundary::default(); 6],
            initial_temperature_k: 298.15,
            maximum_iterations: 2000,
            relative_tolerance: 1e-12,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermalStepResult {
    pub supplied_power_w: f64,
    pub outward_cooling_power_w: f64,
    pub stored_energy_change_j: f64,
    pub energy_balance_error_j: f64,
    pub maximum_temperature_k: f64,
    pub iterations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct ConductionFace {
    a: usize,
    b: usize,
    conductance: f64,
}

/// Implicit-Euler voxel conduction with harmonic face conductances and convective boundaries,
/// solved by matrix-free CG.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalGrid3D {
    faces: Vec<ConductionFace>,
    config: ThermalGridConfig,
    temperature: Vec<f64>,
    capacity: Vec<f64>,
    cooling: Vec<f64>,
    ambient_source: Vec<f64>,
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

impl ThermalGrid3D {
    /// One material is broadcast to every voxel (upstream default `[ThermalVoxel::default()]`).
    pub fn new(config: ThermalGridConfig, mut materials: Vec<ThermalVoxel>) -> Result<Self> {
        let c = config;
        let mut n: usize = 1;
        let mut volume = 1.0;
        for a in 0..3 {
            n = n
                .checked_mul(c.cells[a])
                .filter(|v| *v > 0 && *v <= 1 << 26)
                .ok_or_else(|| invalid("invalid thermal grid dimensions"))?;
            require(positive(c.spacing_m[a]), "invalid thermal grid dimensions")?;
            volume *= c.spacing_m[a];
        }
        require(
            positive(volume)
                && positive(c.initial_temperature_k)
                && c.maximum_iterations > 0
                && positive(c.relative_tolerance)
                && c.relative_tolerance < 1.0,
            "invalid thermal solver configuration",
        )?;
        require(
            c.boundaries
                .iter()
                .all(|b| nonnegative(b.convection_w_per_m2_k) && positive(b.ambient_temperature_k)),
            "invalid thermal boundary",
        )?;
        if materials.len() == 1 {
            materials = vec![materials[0]; n];
        }
        require(materials.len() == n, "thermal material count mismatch")?;
        let mut capacity = vec![0.0; n];
        for (cap, m) in capacity.iter_mut().zip(&materials) {
            *cap = volume * m.volumetric_heat_capacity_j_per_m3_k;
            require(positive(*cap), "invalid voxel heat capacity")?;
            require(
                m.conductivity_w_per_m_k.iter().all(|k| nonnegative(*k)),
                "invalid voxel conductivity",
            )?;
        }
        let stride = [1, c.cells[0], c.cells[0] * c.cells[1]];
        let (mut faces, mut cooling, mut ambient) = (Vec::new(), vec![0.0; n], vec![0.0; n]);
        for i in 0..n {
            for a in 0..3 {
                let coordinate = (i / stride[a]) % c.cells[a];
                let area = volume / c.spacing_m[a];
                let k = materials[i].conductivity_w_per_m_k[a];
                if coordinate + 1 < c.cells[a] {
                    let j = i + stride[a];
                    let kj = materials[j].conductivity_w_per_m_k[a];
                    let g = if k > 0.0 && kj > 0.0 {
                        area / (0.5 * c.spacing_m[a] / k + 0.5 * c.spacing_m[a] / kj)
                    } else {
                        0.0
                    };
                    require(nonnegative(g), "thermal face conductance overflow")?;
                    faces.push(ConductionFace {
                        a: i,
                        b: j,
                        conductance: g,
                    });
                }
                for side in 0..2 {
                    let on = if side == 0 {
                        coordinate == 0
                    } else {
                        coordinate + 1 == c.cells[a]
                    };
                    if on {
                        let b = c.boundaries[2 * a + side];
                        let h = b.convection_w_per_m2_k;
                        let g = if h > 0.0 && k > 0.0 {
                            area / (0.5 * c.spacing_m[a] / k + 1.0 / h)
                        } else {
                            0.0
                        };
                        cooling[i] += g;
                        ambient[i] += g * b.ambient_temperature_k;
                    }
                }
                require(
                    nonnegative(cooling[i]) && nonnegative(ambient[i]),
                    "thermal boundary conductance overflow",
                )?;
            }
        }
        Ok(Self {
            faces,
            config,
            temperature: vec![c.initial_temperature_k; n],
            capacity,
            cooling,
            ambient_source: ambient,
        })
    }
    pub fn reset(&mut self, temperature_k: &[f64]) -> Result<()> {
        require(
            temperature_k.len() == self.temperature.len()
                && temperature_k.iter().all(|t| positive(*t)),
            "invalid thermal reset field",
        )?;
        self.temperature = temperature_k.to_vec();
        Ok(())
    }
    fn conduction(&self, t: &[f64], out: &mut [f64]) {
        for i in 0..t.len() {
            out[i] = self.cooling[i] * t[i];
        }
        for f in &self.faces {
            let flux = f.conductance * (t[f.a] - t[f.b]);
            out[f.a] += flux;
            out[f.b] -= flux;
        }
    }
    /// Transactional implicit step with the energy-balance diagnostic.
    pub fn step(&mut self, dt: f64, voxel_power_w: &[f64]) -> Result<ThermalStepResult> {
        let n = self.temperature.len();
        require(
            positive(dt) && voxel_power_w.len() == n,
            "invalid thermal step",
        )?;
        let mut loss = vec![0.0; n];
        self.conduction(&self.temperature, &mut loss);
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            rhs[i] = dt * (voxel_power_w[i] + self.ambient_source[i] - loss[i]);
            require(
                voxel_power_w[i].is_finite() && rhs[i].is_finite(),
                "nonfinite thermal power",
            )?;
        }
        let mut delta = vec![0.0; n];
        let solved = matrix_free_conjugate_gradient(
            &rhs,
            &mut delta,
            |x, y| {
                self.conduction(x, y);
                for i in 0..n {
                    y[i] = self.capacity[i] * x[i] + dt * y[i];
                }
            },
            self.config.maximum_iterations,
            self.config.relative_tolerance,
        )?;
        require(solved.converged, "3D thermal CG did not converge")?;
        let mut next = self.temperature.clone();
        let mut out = ThermalStepResult {
            iterations: solved.iterations,
            ..Default::default()
        };
        for i in 0..n {
            next[i] += delta[i];
            require(positive(next[i]), "invalid thermal temperature")?;
            out.supplied_power_w += voxel_power_w[i];
            out.outward_cooling_power_w += self.cooling[i] * next[i] - self.ambient_source[i];
            out.stored_energy_change_j += self.capacity[i] * delta[i];
            out.maximum_temperature_k = out.maximum_temperature_k.max(next[i]);
        }
        out.energy_balance_error_j =
            out.stored_energy_change_j - dt * (out.supplied_power_w - out.outward_cooling_power_w);
        require(
            out.energy_balance_error_j.is_finite(),
            "thermal energy overflow",
        )?;
        self.temperature = next;
        Ok(out)
    }
    pub fn temperature(&self) -> &[f64] {
        &self.temperature
    }
    pub fn heat_capacity_j_per_k(&self) -> &[f64] {
        &self.capacity
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ElectrothermalPackResult {
    pub electrical: PackStepResult,
    pub thermal: ThermalStepResult,
}

/// Series pack whose cell heat drives a voxel grid; mapped voxel temperatures feed back as each
/// cell's external temperature.
#[derive(Clone)]
pub struct ElectrothermalBatteryPack {
    pack: SeriesBatteryPack,
    thermal: ThermalGrid3D,
    regions: Vec<Vec<usize>>,
    maximum_temperature_k: f64,
}

impl ElectrothermalBatteryPack {
    pub fn new(
        cells: &[LithiumIonCellConfig],
        control: PackControl,
        thermal: ThermalGrid3D,
        cell_voxels: Vec<Vec<usize>>,
    ) -> Result<Self> {
        let pack = SeriesBatteryPack::new(cells, control)?;
        require(
            cell_voxels.len() == cells.len(),
            "one thermal region required per cell",
        )?;
        let mut used = vec![false; thermal.temperature().len()];
        for (region, cell) in cell_voxels.iter().zip(cells) {
            require(!region.is_empty(), "empty cell thermal region")?;
            let mut capacity = 0.0;
            for &i in region {
                require(
                    i < used.len() && !used[i],
                    "overlapping or invalid cell thermal voxel",
                )?;
                used[i] = true;
                capacity += thermal.heat_capacity_j_per_k()[i];
            }
            let expected = cell.thermal.mass_kg * cell.thermal.heat_capacity_j_per_kg_k;
            require(
                positive(expected)
                    && positive(capacity)
                    && (capacity - expected).abs() <= 1e-8 * expected,
                "cell and mapped thermal heat capacities must agree",
            )?;
        }
        let hottest = thermal
            .temperature()
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        if hottest > control.maximum_temperature_k {
            return Err(SolverError::NumericalFailure(
                "initial pack thermal hotspot cutoff".into(),
            ));
        }
        let mut s = Self {
            pack,
            thermal,
            regions: cell_voxels,
            maximum_temperature_k: control.maximum_temperature_k,
        };
        let t = s.cell_temperatures(&s.thermal);
        s.pack.set_external_temperatures(&t)?;
        Ok(s)
    }
    fn cell_temperatures(&self, field: &ThermalGrid3D) -> Vec<f64> {
        self.regions
            .iter()
            .map(|region| {
                let (mut energy, mut capacity) = (0.0, 0.0);
                for &i in region {
                    let c = field.heat_capacity_j_per_k()[i];
                    energy += c * field.temperature()[i];
                    capacity += c;
                }
                energy / capacity
            })
            .collect()
    }
    pub fn step(&mut self, current: f64, dt: f64) -> Result<ElectrothermalPackResult> {
        let (mut trial_pack, mut trial_thermal) = (self.pack.clone(), self.thermal.clone());
        let mut electrical = trial_pack.step(current, dt)?;
        let mut q = vec![0.0; self.thermal.temperature().len()];
        for (c, region) in self.regions.iter().enumerate() {
            let power = electrical.cells[c].heat_generation_w / region.len() as f64;
            for &i in region {
                q[i] += power;
            }
        }
        let thermal = trial_thermal.step(dt, &q)?;
        if thermal.maximum_temperature_k > self.maximum_temperature_k {
            return Err(SolverError::NumericalFailure(
                "pack thermal hotspot cutoff".into(),
            ));
        }
        let temperatures = self.cell_temperatures(&trial_thermal);
        trial_pack.set_external_temperatures(&temperatures)?;
        for (cell, t) in electrical.cells.iter_mut().zip(&temperatures) {
            cell.temperature_k = *t;
        }
        self.pack = trial_pack;
        self.thermal = trial_thermal;
        Ok(ElectrothermalPackResult {
            electrical,
            thermal,
        })
    }
    pub fn electrical(&self) -> &SeriesBatteryPack {
        &self.pack
    }
    pub fn thermal(&self) -> &ThermalGrid3D {
        &self.thermal
    }
}
