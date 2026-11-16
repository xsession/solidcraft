use crate::{positive, require, Result, SolverError};
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SellmeierMaterial {
    pub name: String,
    pub b: [f64; 3],
    pub c_um2: [f64; 3],
}
impl SellmeierMaterial {
    pub fn refractive_index_nm(&self, wavelength: f64) -> Result<f64> {
        positive(wavelength, "wavelength")?;
        let l2 = (1e-3 * wavelength).powi(2);
        let mut n2 = 1.0;
        for i in 0..3 {
            let d = l2 - self.c_um2[i];
            require(d.abs() >= 1e-14, "Sellmeier wavelength at resonance")?;
            n2 += self.b[i] * l2 / d;
        }
        positive(n2, "Sellmeier index squared")?;
        Ok(n2.sqrt())
    }
}
#[derive(Clone, Debug)]
pub struct MaterialCatalog {
    materials: Vec<SellmeierMaterial>,
}
impl Default for MaterialCatalog {
    fn default() -> Self {
        Self {
            materials: vec![
                SellmeierMaterial {
                    name: "N-BK7".into(),
                    b: [1.03961212, 0.231792344, 1.01046945],
                    c_um2: [0.00600069867, 0.0200179144, 103.560653],
                },
                SellmeierMaterial {
                    name: "FusedSilica".into(),
                    b: [0.6961663, 0.4079426, 0.8974794],
                    c_um2: [
                        0.0684043_f64.powi(2),
                        0.1162414_f64.powi(2),
                        9.896161_f64.powi(2),
                    ],
                },
            ],
        }
    }
}
impl MaterialCatalog {
    pub fn add(&mut self, m: SellmeierMaterial) -> Result<()> {
        require(
            !m.name.is_empty() && !self.materials.iter().any(|v| v.name == m.name),
            "empty or duplicate material name",
        )?;
        require(
            m.b.iter().chain(&m.c_um2).all(|v| v.is_finite()),
            "nonfinite Sellmeier coefficients",
        )?;
        self.materials.push(m);
        Ok(())
    }
    pub fn at(&self, name: &str) -> Result<&SellmeierMaterial> {
        self.materials
            .iter()
            .find(|m| m.name == name)
            .ok_or_else(|| SolverError::InvalidInput(format!("unknown optical material {name}")))
    }
    pub fn materials(&self) -> &[SellmeierMaterial] {
        &self.materials
    }
    pub fn nearest_index(&self, wavelength: f64, target: f64) -> Result<&SellmeierMaterial> {
        positive(target, "target index")?;
        let mut best = None;
        let mut error = f64::INFINITY;
        for m in &self.materials {
            let e = (m.refractive_index_nm(wavelength)? - target).abs();
            if e < error {
                error = e;
                best = Some(m);
            }
        }
        best.ok_or_else(|| SolverError::InvalidInput("empty material catalogue".into()))
    }
}
