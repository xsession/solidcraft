use crate::{positive, require, Complex64 as C, Result};
use std::f64::consts::PI;
#[derive(Clone, Copy, Debug)]
pub struct GaussianBeam {
    wavelength: f64,
    index: f64,
    q: C,
}
impl GaussianBeam {
    pub fn new(wavelength: f64, waist: f64, distance: f64, index: f64) -> Result<Self> {
        positive(wavelength, "wavelength")?;
        positive(waist, "waist radius")?;
        positive(index, "refractive index")?;
        require(distance.is_finite(), "nonfinite beam distance")?;
        let q = C::new(distance, PI * index * waist * waist / wavelength);
        positive(q.im, "Rayleigh range")?;
        Ok(Self {
            wavelength,
            index,
            q,
        })
    }
    pub fn q(self) -> C {
        self.q
    }
    pub fn refractive_index(self) -> f64 {
        self.index
    }
    pub fn transform(&mut self, a: f64, b: f64, c: f64, d: f64, output_index: f64) -> Result<()> {
        positive(output_index, "output index")?;
        require([a, b, c, d].iter().all(|v| v.is_finite()), "nonfinite ABCD")?;
        let det = a * d - b * c;
        let expected = self.index / output_index;
        require(
            det > 0.0 && (det - expected).abs() <= 1e-10 * expected,
            "ABCD determinant must equal input/output index ratio",
        )?;
        let q = (a * self.q + b) / (c * self.q + d);
        require(
            q.re.is_finite() && q.im.is_finite() && q.im > 0.0,
            "singular Gaussian beam transform",
        )?;
        self.q = q;
        self.index = output_index;
        Ok(())
    }
    pub fn propagate(&mut self, distance: f64) -> Result<()> {
        self.transform(1.0, distance, 0.0, 1.0, self.index)
    }
    pub fn thin_lens(&mut self, f: f64) -> Result<()> {
        require(f.is_finite() && f != 0.0, "invalid focal length")?;
        self.transform(1.0, 0.0, -1.0 / f, 1.0, self.index)
    }
    pub fn radius(self) -> f64 {
        (-self.wavelength / (PI * self.index * (1.0 / self.q).im)).sqrt()
    }
    pub fn curvature_radius(self) -> f64 {
        let inv = (1.0 / self.q).re;
        if inv == 0.0 {
            f64::INFINITY
        } else {
            1.0 / inv
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct JonesVector {
    pub x: C,
    pub y: C,
}
#[derive(Clone, Copy, Debug)]
pub struct StokesVector {
    pub i: f64,
    pub q: f64,
    pub u: f64,
    pub v: f64,
}
pub fn stokes_from_jones(j: JonesVector) -> StokesVector {
    let (ex, ey) = (j.x.norm_sqr(), j.y.norm_sqr());
    let c = j.x * j.y.conj();
    StokesVector {
        i: ex + ey,
        q: ex - ey,
        u: 2.0 * c.re,
        v: -2.0 * c.im,
    }
}
#[derive(Clone, Copy, Debug)]
pub struct FresnelCoefficients {
    pub cos_transmitted: f64,
    pub rs: f64,
    pub rp: f64,
    pub ts: f64,
    pub tp: f64,
    pub reflectance_s: f64,
    pub reflectance_p: f64,
    pub total_internal_reflection: bool,
}
pub fn fresnel_dielectric(n1: f64, n2: f64, cos_i: f64) -> Result<FresnelCoefficients> {
    positive(n1, "incident index")?;
    positive(n2, "transmitted index")?;
    require(cos_i.is_finite(), "nonfinite angle cosine")?;
    let ci = cos_i.clamp(0.0, 1.0);
    let st2 = (n1 / n2).powi(2) * (1.0 - ci * ci);
    if st2 > 1.0 {
        return Ok(FresnelCoefficients {
            cos_transmitted: 0.0,
            rs: 1.0,
            rp: 1.0,
            ts: 0.0,
            tp: 0.0,
            reflectance_s: 1.0,
            reflectance_p: 1.0,
            total_internal_reflection: true,
        });
    }
    let ct = (1.0 - st2).max(0.0).sqrt();
    let ds = n1 * ci + n2 * ct;
    let dp = n2 * ci + n1 * ct;
    require(ds > 0.0 && dp > 0.0, "singular grazing Fresnel interface")?;
    let rs = (n1 * ci - n2 * ct) / ds;
    let rp = (n2 * ci - n1 * ct) / dp;
    Ok(FresnelCoefficients {
        cos_transmitted: ct,
        rs,
        rp,
        ts: 2.0 * n1 * ci / ds,
        tp: 2.0 * n1 * ci / dp,
        reflectance_s: rs * rs,
        reflectance_p: rp * rp,
        total_internal_reflection: false,
    })
}
#[derive(Clone, Copy, Debug)]
pub struct ThinFilmLayer {
    pub refractive_index: f64,
    pub thickness_nm: f64,
}
pub fn thin_film_reflectance_normal(
    n0: f64,
    ns: f64,
    wavelength: f64,
    layers: &[ThinFilmLayer],
) -> Result<f64> {
    positive(n0, "incident index")?;
    positive(ns, "substrate index")?;
    positive(wavelength, "wavelength")?;
    let mut m = crate::rf::Matrix2C::identity();
    for l in layers {
        positive(l.refractive_index, "layer index")?;
        require(
            l.thickness_nm.is_finite() && l.thickness_nm >= 0.0,
            "invalid film thickness",
        )?;
        let delta = 2.0 * PI * l.refractive_index * l.thickness_nm / wavelength;
        let a = C::new(delta.cos(), 0.0);
        m = m.multiply(crate::rf::Matrix2C {
            a11: a,
            a12: C::new(0.0, delta.sin() / l.refractive_index),
            a21: C::new(0.0, l.refractive_index * delta.sin()),
            a22: a,
        });
    }
    let numerator = n0 * m.a11 + n0 * ns * m.a12 - m.a21 - ns * m.a22;
    let denominator = n0 * m.a11 + n0 * ns * m.a12 + m.a21 + ns * m.a22;
    let r = (numerator / denominator).norm_sqr();
    require(r.is_finite(), "thin-film transfer overflow")?;
    Ok(r)
}
