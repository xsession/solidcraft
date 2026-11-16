use crate::{positive, require, Result};
#[derive(Clone, Copy, Debug)]
pub struct MaterialPoint1D {
    pub stress: f64,
    pub tangent: f64,
    pub plastic_strain: f64,
    pub equivalent_plastic_strain: f64,
}
pub trait UniaxialMaterialLaw {
    fn update(&mut self, total_strain: f64) -> Result<MaterialPoint1D>;
}
#[derive(Clone, Copy, Debug)]
pub struct LinearElastic1D {
    youngs: f64,
}
impl LinearElastic1D {
    pub fn new(youngs: f64) -> Result<Self> {
        positive(youngs, "Young modulus")?;
        Ok(Self { youngs })
    }
}
impl UniaxialMaterialLaw for LinearElastic1D {
    fn update(&mut self, strain: f64) -> Result<MaterialPoint1D> {
        require(strain.is_finite(), "nonfinite strain")?;
        Ok(MaterialPoint1D {
            stress: self.youngs * strain,
            tangent: self.youngs,
            plastic_strain: 0.0,
            equivalent_plastic_strain: 0.0,
        })
    }
}
#[derive(Clone, Copy, Debug)]
pub struct BilinearPlasticity1D {
    youngs: f64,
    yield_stress: f64,
    hardening: f64,
    plastic: f64,
    alpha: f64,
}
impl BilinearPlasticity1D {
    pub fn new(youngs: f64, yield_stress: f64, hardening: f64) -> Result<Self> {
        positive(youngs, "Young modulus")?;
        positive(yield_stress, "yield stress")?;
        require(
            hardening.is_finite() && hardening >= 0.0,
            "invalid hardening modulus",
        )?;
        Ok(Self {
            youngs,
            yield_stress,
            hardening,
            plastic: 0.0,
            alpha: 0.0,
        })
    }
}
impl UniaxialMaterialLaw for BilinearPlasticity1D {
    fn update(&mut self, strain: f64) -> Result<MaterialPoint1D> {
        require(strain.is_finite(), "nonfinite strain")?;
        let mut trial = self.youngs * (strain - self.plastic);
        let f = trial.abs() - (self.yield_stress + self.hardening * self.alpha);
        let mut tangent = self.youngs;
        if f > 0.0 {
            let dg = f / (self.youngs + self.hardening);
            let sign = if trial >= 0.0 { 1.0 } else { -1.0 };
            self.plastic += dg * sign;
            self.alpha += dg;
            trial -= self.youngs * dg * sign;
            tangent = self.youngs * self.hardening / (self.youngs + self.hardening);
        }
        Ok(MaterialPoint1D {
            stress: trial,
            tangent,
            plastic_strain: self.plastic,
            equivalent_plastic_strain: self.alpha,
        })
    }
}
#[derive(Clone, Copy, Debug)]
pub struct NeoHookeanResult {
    pub energy_density: f64,
    pub nominal_stress: [f64; 3],
}
pub fn compressible_neo_hookean(
    stretch: [f64; 3],
    shear_modulus: f64,
    lame_lambda: f64,
) -> Result<NeoHookeanResult> {
    positive(shear_modulus, "shear modulus")?;
    require(
        lame_lambda.is_finite()
            && lame_lambda >= 0.0
            && stretch.iter().all(|v| v.is_finite() && *v > 0.0),
        "invalid Neo-Hookean stretches/moduli",
    )?;
    let j = stretch.iter().product::<f64>();
    positive(j, "Jacobian")?;
    let log = j.ln();
    let i1 = stretch.iter().map(|v| v * v).sum::<f64>();
    let energy_density =
        0.5 * shear_modulus * (i1 - 3.0) - shear_modulus * log + 0.5 * lame_lambda * log * log;
    let nominal_stress = stretch.map(|v| shear_modulus * (v - 1.0 / v) + lame_lambda * log / v);
    require(
        energy_density.is_finite() && nominal_stress.iter().all(|v| v.is_finite()),
        "Neo-Hookean overflow",
    )?;
    Ok(NeoHookeanResult {
        energy_density,
        nominal_stress,
    })
}
pub fn penalty_contact_pressure(signed_gap: f64, penalty_stiffness: f64) -> Result<f64> {
    positive(penalty_stiffness, "penalty stiffness")?;
    require(signed_gap.is_finite(), "nonfinite gap")?;
    Ok(if signed_gap < 0.0 {
        -penalty_stiffness * signed_gap
    } else {
        0.0
    })
}
