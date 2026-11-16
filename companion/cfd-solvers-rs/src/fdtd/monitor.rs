//! Time probes, running DFT monitors, 1-D wave-port decomposition and S-parameters, ported from
//! upstream `cfd/solvers/fdtd/{monitor, port}` (MIT, cfd_solvers contributors).
use crate::{require, Complex64 as C, Result};
use std::f64::consts::PI;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TimeProbe {
    time: Vec<f64>,
    value: Vec<f64>,
}

impl TimeProbe {
    pub fn sample(&mut self, time: f64, value: f64) -> Result<()> {
        require(
            time.is_finite() && value.is_finite(),
            "non-finite FDTD probe sample",
        )?;
        require(
            self.time.last().is_none_or(|last| time >= *last),
            "FDTD probe time must be monotone",
        )?;
        self.time.push(time);
        self.value.push(value);
        Ok(())
    }
    pub fn clear(&mut self) {
        self.time.clear();
        self.value.clear();
    }
    pub fn time(&self) -> &[f64] {
        &self.time
    }
    pub fn value(&self) -> &[f64] {
        &self.value
    }
}

/// Trapezoidal running `integral v(t) exp(-i w t) dt` at one frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DftMonitor {
    frequency: f64,
    integral: C,
    previous_time: f64,
    previous_value: f64,
    have_previous: bool,
    samples: usize,
}

impl DftMonitor {
    pub fn new(frequency_hz: f64) -> Result<Self> {
        require(
            frequency_hz.is_finite() && frequency_hz >= 0.0,
            "invalid DFT frequency",
        )?;
        Ok(Self {
            frequency: frequency_hz,
            integral: C::new(0.0, 0.0),
            previous_time: 0.0,
            previous_value: 0.0,
            have_previous: false,
            samples: 0,
        })
    }
    pub fn sample(&mut self, time: f64, value: f64) -> Result<()> {
        require(
            time.is_finite() && value.is_finite(),
            "non-finite DFT sample",
        )?;
        let omega = 2.0 * PI * self.frequency;
        let phase = C::new(0.0, -omega * time).exp();
        if self.have_previous {
            let dt = time - self.previous_time;
            require(dt >= 0.0, "DFT sample time must be monotone")?;
            let prev_phase = C::new(0.0, -omega * self.previous_time).exp();
            self.integral += (prev_phase * self.previous_value + phase * value) * (0.5 * dt);
        }
        self.previous_time = time;
        self.previous_value = value;
        self.have_previous = true;
        self.samples += 1;
        Ok(())
    }
    pub fn clear(&mut self) {
        *self = Self {
            frequency: self.frequency,
            integral: C::new(0.0, 0.0),
            previous_time: 0.0,
            previous_value: 0.0,
            have_previous: false,
            samples: 0,
        };
    }
    pub fn integral(&self) -> C {
        self.integral
    }
    /// `2 |integral| / duration`: real-sinusoid amplitude estimate.
    pub fn amplitude(&self, duration: f64) -> Result<f64> {
        require(duration > 0.0, "DFT duration must be positive")?;
        Ok(2.0 * self.integral.norm() / duration)
    }
    pub fn samples(&self) -> usize {
        self.samples
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WaveAmplitudes1D {
    pub forward: f64,
    pub backward: f64,
}

/// TEM-like port for the `Maxwell1D` Ez/Hy convention: `+x` waves have `Hy = -Ez / Z`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WavePort1D {
    frequency_hz: f64,
    impedance_ohm: f64,
    direction: i32,
    forward: DftMonitor,
    backward: DftMonitor,
}

impl WavePort1D {
    /// `propagation_direction` is `+1` or `-1` (upstream default `+1`).
    pub fn new(frequency_hz: f64, impedance_ohm: f64, propagation_direction: i32) -> Result<Self> {
        require(
            frequency_hz.is_finite()
                && frequency_hz >= 0.0
                && impedance_ohm.is_finite()
                && impedance_ohm > 0.0
                && (propagation_direction == 1 || propagation_direction == -1),
            "invalid 1-D wave-port configuration",
        )?;
        Ok(Self {
            frequency_hz,
            impedance_ohm,
            direction: propagation_direction,
            forward: DftMonitor::new(frequency_hz)?,
            backward: DftMonitor::new(frequency_hz)?,
        })
    }
    pub fn instantaneous(&self, electric_z: f64, magnetic_y: f64) -> WaveAmplitudes1D {
        let signed_zh = f64::from(self.direction) * self.impedance_ohm * magnetic_y;
        WaveAmplitudes1D {
            forward: 0.5 * (electric_z - signed_zh),
            backward: 0.5 * (electric_z + signed_zh),
        }
    }
    pub fn sample(&mut self, time: f64, electric_z: f64, magnetic_y: f64) -> Result<()> {
        require(
            electric_z.is_finite() && magnetic_y.is_finite(),
            "non-finite wave-port sample",
        )?;
        let a = self.instantaneous(electric_z, magnetic_y);
        // Validate both before mutating either monitor.
        require(
            time.is_finite() && (!self.forward.have_previous || time >= self.forward.previous_time),
            "DFT sample time must be monotone",
        )?;
        self.forward.sample(time, a.forward)?;
        self.backward.sample(time, a.backward)
    }
    pub fn clear(&mut self) {
        self.forward.clear();
        self.backward.clear();
    }
    pub fn forward_spectrum(&self) -> C {
        self.forward.integral()
    }
    pub fn backward_spectrum(&self) -> C {
        self.backward.integral()
    }
    pub fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub fn impedance_ohm(&self) -> f64 {
        self.impedance_ohm
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SParameters1D {
    pub s11: C,
    pub s21: C,
}

/// `S11 = b_ref / a_ref`, `S21 = a_trans / a_ref`.
pub fn s_parameters(reference: &WavePort1D, transmitted: &WavePort1D) -> Result<SParameters1D> {
    let incident = reference.forward_spectrum();
    require(
        incident.norm() > 1e-30,
        "cannot compute S-parameters with zero incident spectrum",
    )?;
    Ok(SParameters1D {
        s11: reference.backward_spectrum() / incident,
        s21: transmitted.forward_spectrum() / incident,
    })
}
