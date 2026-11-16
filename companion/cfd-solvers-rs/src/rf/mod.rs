//! Two-port network algebra, transmission-line approximations and Touchstone S2P interchange.
mod antenna;
mod asymptotic;
mod nport;
mod peec;
mod surface_mom;
mod thin_wire;
use crate::{positive, require, Complex64 as C, Result, SolverError};
pub use antenna::*;
pub use asymptotic::*;
pub use nport::*;
pub use peec::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
pub use surface_mom::*;
pub use thin_wire::*;
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
pub struct Matrix2C {
    pub a11: C,
    pub a12: C,
    pub a21: C,
    pub a22: C,
}
impl Matrix2C {
    pub fn identity() -> Self {
        Self {
            a11: C::new(1.0, 0.0),
            a22: C::new(1.0, 0.0),
            ..Default::default()
        }
    }
    pub fn values(self) -> [C; 4] {
        [self.a11, self.a12, self.a21, self.a22]
    }
    pub fn inverse(self) -> Result<Self> {
        let d = self.a11 * self.a22 - self.a12 * self.a21;
        require(d.norm() >= 1e-30, "singular RF matrix")?;
        Ok(Self {
            a11: self.a22 / d,
            a12: -self.a12 / d,
            a21: -self.a21 / d,
            a22: self.a11 / d,
        })
    }
    pub fn multiply(self, b: Self) -> Self {
        Self {
            a11: self.a11 * b.a11 + self.a12 * b.a21,
            a12: self.a11 * b.a12 + self.a12 * b.a22,
            a21: self.a21 * b.a11 + self.a22 * b.a21,
            a22: self.a21 * b.a12 + self.a22 * b.a22,
        }
    }
    pub fn validate(self) -> Result<()> {
        require(
            self.values()
                .iter()
                .all(|z| z.re.is_finite() && z.im.is_finite()),
            "nonfinite RF matrix",
        )
    }
}
pub fn z_to_s(z: Matrix2C, z0: f64) -> Result<Matrix2C> {
    positive(z0, "reference impedance")?;
    z.validate()?;
    let n = Matrix2C {
        a11: z.a11 - z0,
        a22: z.a22 - z0,
        ..z
    };
    let d = Matrix2C {
        a11: z.a11 + z0,
        a22: z.a22 + z0,
        ..z
    };
    Ok(n.multiply(d.inverse()?))
}
pub fn s_to_z(s: Matrix2C, z0: f64) -> Result<Matrix2C> {
    positive(z0, "reference impedance")?;
    s.validate()?;
    let plus = Matrix2C {
        a11: 1.0 + s.a11,
        a22: 1.0 + s.a22,
        ..s
    };
    let minus = Matrix2C {
        a11: 1.0 - s.a11,
        a12: -s.a12,
        a21: -s.a21,
        a22: 1.0 - s.a22,
    };
    let r = plus.multiply(minus.inverse()?);
    Ok(Matrix2C {
        a11: z0 * r.a11,
        a12: z0 * r.a12,
        a21: z0 * r.a21,
        a22: z0 * r.a22,
    })
}
pub fn s_to_abcd(s: Matrix2C, z0: f64) -> Result<Matrix2C> {
    positive(z0, "reference impedance")?;
    s.validate()?;
    require(
        s.a21.norm() >= 1e-30,
        "ABCD conversion requires nonzero S21",
    )?;
    let d = 2.0 * s.a21;
    Ok(Matrix2C {
        a11: ((1.0 + s.a11) * (1.0 - s.a22) + s.a12 * s.a21) / d,
        a12: z0 * ((1.0 + s.a11) * (1.0 + s.a22) - s.a12 * s.a21) / d,
        a21: ((1.0 - s.a11) * (1.0 - s.a22) - s.a12 * s.a21) / (z0 * d),
        a22: ((1.0 - s.a11) * (1.0 + s.a22) + s.a12 * s.a21) / d,
    })
}
pub fn abcd_to_s(m: Matrix2C, z0: f64) -> Result<Matrix2C> {
    positive(z0, "reference impedance")?;
    m.validate()?;
    let d = m.a11 + m.a12 / z0 + m.a21 * z0 + m.a22;
    require(d.norm() >= 1e-30, "singular ABCD conversion")?;
    let det = m.a11 * m.a22 - m.a12 * m.a21;
    Ok(Matrix2C {
        a11: (m.a11 + m.a12 / z0 - m.a21 * z0 - m.a22) / d,
        a12: 2.0 * det / d,
        a21: 2.0 / d,
        a22: (-m.a11 + m.a12 / z0 - m.a21 * z0 + m.a22) / d,
    })
}
pub fn cascade_abcd(a: Matrix2C, b: Matrix2C) -> Matrix2C {
    a.multiply(b)
}
pub fn cascade_s(a: Matrix2C, b: Matrix2C, z0: f64) -> Result<Matrix2C> {
    abcd_to_s(s_to_abcd(a, z0)?.multiply(s_to_abcd(b, z0)?), z0)
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct TransmissionLine {
    pub characteristic_impedance: f64,
    pub propagation_constant: C,
    pub length: f64,
}
impl TransmissionLine {
    pub fn abcd(self) -> Result<Matrix2C> {
        positive(self.characteristic_impedance, "characteristic impedance")?;
        require(
            self.length.is_finite()
                && self.length >= 0.0
                && self.propagation_constant.re.is_finite()
                && self.propagation_constant.im.is_finite(),
            "invalid transmission line",
        )?;
        let gl = self.propagation_constant * self.length;
        let ch = gl.cosh();
        let sh = gl.sinh();
        let m = Matrix2C {
            a11: ch,
            a12: self.characteristic_impedance * sh,
            a21: sh / self.characteristic_impedance,
            a22: ch,
        };
        m.validate()?;
        Ok(m)
    }
    pub fn s_parameters(self, z0: f64) -> Result<Matrix2C> {
        abcd_to_s(self.abcd()?, z0)
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct TemTransmissionLine {
    pub characteristic_impedance: f64,
    pub phase_velocity: f64,
    pub attenuation_nepers_per_m: f64,
    pub length: f64,
}
impl TemTransmissionLine {
    pub fn at_frequency(self, frequency: f64) -> Result<TransmissionLine> {
        positive(frequency, "frequency")?;
        positive(self.phase_velocity, "phase velocity")?;
        positive(self.characteristic_impedance, "characteristic impedance")?;
        require(
            self.attenuation_nepers_per_m.is_finite()
                && self.attenuation_nepers_per_m >= 0.0
                && self.length.is_finite()
                && self.length >= 0.0,
            "invalid TEM line",
        )?;
        Ok(TransmissionLine {
            characteristic_impedance: self.characteristic_impedance,
            propagation_constant: C::new(
                self.attenuation_nepers_per_m,
                2.0 * PI * frequency / self.phase_velocity,
            ),
            length: self.length,
        })
    }
    pub fn s_parameters(self, frequency: f64, z0: f64) -> Result<Matrix2C> {
        self.at_frequency(frequency)?.s_parameters(z0)
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct PlanarLineQuasiStatic {
    pub characteristic_impedance: f64,
    pub effective_permittivity: f64,
    pub phase_velocity: f64,
    pub guided_wavelength: f64,
}
pub type MicrostripQuasiStatic = PlanarLineQuasiStatic;
pub fn microstrip_quasi_static(
    width: f64,
    height: f64,
    er: f64,
    frequency: f64,
) -> Result<MicrostripQuasiStatic> {
    positive(width, "width")?;
    positive(height, "height")?;
    positive(frequency, "frequency")?;
    require(
        er.is_finite() && er > 1.0,
        "microstrip permittivity must exceed one",
    )?;
    let u = width / height;
    let correction = if u < 1.0 {
        0.04 * (1.0 - u).powi(2)
    } else {
        0.0
    };
    let ee = 0.5 * (er + 1.0) + 0.5 * (er - 1.0) * (1.0 / (1.0 + 12.0 / u).sqrt() + correction);
    let z = if u <= 1.0 {
        60.0 / ee.sqrt() * (8.0 / u + 0.25 * u).ln()
    } else {
        120.0 * PI / (ee.sqrt() * (u + 1.393 + 0.667 * (u + 1.444).ln()))
    };
    let v = 299792458.0 / ee.sqrt();
    Ok(PlanarLineQuasiStatic {
        characteristic_impedance: z,
        effective_permittivity: ee,
        phase_velocity: v,
        guided_wavelength: v / frequency,
    })
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct CoaxQuasiStatic {
    pub characteristic_impedance: f64,
    pub phase_velocity: f64,
    pub propagation_constant_rad_per_m: f64,
    pub guided_wavelength: f64,
}
pub fn coax_quasi_static(
    inner: f64,
    outer: f64,
    er: f64,
    frequency: f64,
) -> Result<CoaxQuasiStatic> {
    positive(inner, "inner radius")?;
    require(
        outer.is_finite() && outer > inner,
        "outer radius must exceed inner radius",
    )?;
    positive(er, "permittivity")?;
    positive(frequency, "frequency")?;
    let v = 299792458.0 / er.sqrt();
    Ok(CoaxQuasiStatic {
        characteristic_impedance: 376.730313668 / (2.0 * PI * er.sqrt()) * (outer / inner).ln(),
        phase_velocity: v,
        propagation_constant_rad_per_m: 2.0 * PI * frequency / v,
        guided_wavelength: v / frequency,
    })
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct RectangularWaveguideTE10 {
    pub cutoff_frequency_hz: f64,
    pub propagation_constant_rad_per_m: f64,
    pub guide_wavelength: f64,
    pub wave_impedance: f64,
}
pub fn rectangular_waveguide_te10(
    a: f64,
    b: f64,
    er: f64,
    mur: f64,
    frequency: f64,
) -> Result<RectangularWaveguideTE10> {
    for (v, n) in [
        (a, "broad wall"),
        (b, "narrow wall"),
        (er, "permittivity"),
        (mur, "permeability"),
        (frequency, "frequency"),
    ] {
        positive(v, n)?;
    }
    let factor = (er * mur).sqrt();
    let cutoff = 299792458.0 / (2.0 * a * factor);
    require(frequency > cutoff, "TE10 frequency below cutoff")?;
    let root = (1.0 - (cutoff / frequency).powi(2)).sqrt();
    let beta = 2.0 * PI * frequency * factor / 299792458.0 * root;
    Ok(RectangularWaveguideTE10 {
        cutoff_frequency_hz: cutoff,
        propagation_constant_rad_per_m: beta,
        guide_wavelength: 2.0 * PI / beta,
        wave_impedance: 376.730313668 * (mur / er).sqrt() / root,
    })
}
pub fn stripline_quasi_static(
    width: f64,
    spacing: f64,
    er: f64,
    frequency: f64,
) -> Result<PlanarLineQuasiStatic> {
    for (v, n) in [
        (width, "width"),
        (spacing, "ground spacing"),
        (er, "permittivity"),
        (frequency, "frequency"),
    ] {
        positive(v, n)?;
    }
    let z = 30.0 * PI / er.sqrt() / (width / spacing + 0.441);
    let v = 299792458.0 / er.sqrt();
    Ok(PlanarLineQuasiStatic {
        characteristic_impedance: z,
        effective_permittivity: er,
        phase_velocity: v,
        guided_wavelength: v / frequency,
    })
}
fn elliptic_k(k: f64) -> Result<f64> {
    require((0.0..1.0).contains(&k), "elliptic modulus out of range")?;
    let (mut a, mut b) = (1.0, (1.0 - k * k).sqrt());
    for _ in 0..32 {
        if (a - b).abs() <= 1e-15 * a {
            break;
        }
        let next = 0.5 * (a + b);
        b = (a * b).sqrt();
        a = next;
    }
    Ok(PI / (2.0 * a))
}
pub fn coplanar_waveguide_quasi_static(
    width: f64,
    gap: f64,
    er: f64,
    frequency: f64,
) -> Result<PlanarLineQuasiStatic> {
    for (v, n) in [
        (width, "center width"),
        (gap, "gap"),
        (er, "permittivity"),
        (frequency, "frequency"),
    ] {
        positive(v, n)?;
    }
    let k = width / (width + 2.0 * gap);
    let ee = 0.5 * (er + 1.0);
    let z = 30.0 * PI / ee.sqrt() * elliptic_k((1.0 - k * k).sqrt())? / elliptic_k(k)?;
    let v = 299792458.0 / ee.sqrt();
    Ok(PlanarLineQuasiStatic {
        characteristic_impedance: z,
        effective_permittivity: ee,
        phase_velocity: v,
        guided_wavelength: v / frequency,
    })
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Touchstone2PortPoint {
    pub frequency_hz: f64,
    pub s: Matrix2C,
}
/// Touchstone v1 S2P, including wrapped records, comments, RI/MA/DB and all four frequency units.
pub fn read_touchstone_s2p(input: &str) -> Result<(Vec<Touchstone2PortPoint>, f64)> {
    let (mut scale, mut format, mut z0) = (1e9, "MA".to_string(), 50.0);
    let mut pending = vec![];
    let mut points = vec![];
    for line in input.lines() {
        let line = line.split('!').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if let Some(options) = line.strip_prefix('#') {
            require(pending.is_empty(), "Touchstone option within record")?;
            let words: Vec<_> = options.split_whitespace().map(str::to_uppercase).collect();
            require(words.len() >= 3, "incomplete Touchstone option line")?;
            scale = match words[0].as_str() {
                "HZ" => 1.0,
                "KHZ" => 1e3,
                "MHZ" => 1e6,
                "GHZ" => 1e9,
                _ => return Err(SolverError::Unsupported("Touchstone frequency unit".into())),
            };
            require(words[1] == "S", "Touchstone reader expects S parameters")?;
            format = words[2].clone();
            require(
                matches!(format.as_str(), "RI" | "MA" | "DB"),
                "unsupported Touchstone complex format",
            )?;
            if let Some(i) = words.iter().position(|s| s == "R") {
                z0 = words
                    .get(i + 1)
                    .ok_or_else(|| SolverError::InvalidInput("missing Touchstone R value".into()))?
                    .parse()
                    .map_err(|_| SolverError::InvalidInput("invalid Touchstone R".into()))?;
                positive(z0, "reference impedance")?;
            }
            continue;
        }
        for token in line.split_whitespace() {
            let v = token.parse::<f64>().map_err(|_| {
                SolverError::InvalidInput(format!("invalid Touchstone number {token}"))
            })?;
            require(v.is_finite(), "nonfinite Touchstone number")?;
            pending.push(v);
        }
        while pending.len() >= 9 {
            let decode = |a, b| match format.as_str() {
                "RI" => C::new(a, b),
                "MA" => C::from_polar(a, b * PI / 180.0),
                _ => C::from_polar(10.0_f64.powf(a / 20.0), b * PI / 180.0),
            };
            let s = Matrix2C {
                a11: decode(pending[1], pending[2]),
                a21: decode(pending[3], pending[4]),
                a12: decode(pending[5], pending[6]),
                a22: decode(pending[7], pending[8]),
            };
            s.validate()?;
            let frequency = pending[0] * scale;
            require(
                frequency.is_finite() && frequency >= 0.0,
                "invalid Touchstone frequency",
            )?;
            points.push(Touchstone2PortPoint {
                frequency_hz: frequency,
                s,
            });
            pending.drain(..9);
        }
    }
    require(pending.is_empty(), "truncated Touchstone record")?;
    Ok((points, z0))
}
pub fn write_touchstone_s2p(points: &[Touchstone2PortPoint], z0: f64) -> Result<String> {
    positive(z0, "reference impedance")?;
    let mut out = format!("# HZ S RI R {z0:.17e}\n");
    for p in points {
        require(
            p.frequency_hz.is_finite() && p.frequency_hz >= 0.0,
            "invalid Touchstone frequency",
        )?;
        p.s.validate()?;
        out += &format!(
            "{:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}\n",
            p.frequency_hz,
            p.s.a11.re,
            p.s.a11.im,
            p.s.a21.re,
            p.s.a21.im,
            p.s.a12.re,
            p.s.a12.im,
            p.s.a22.re,
            p.s.a22.im
        );
    }
    Ok(out)
}
