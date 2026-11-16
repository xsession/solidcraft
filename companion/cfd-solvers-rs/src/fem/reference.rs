use super::Point3;
use crate::{require, Result};
use serde::{Deserialize, Serialize};
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ElementType {
    Line2,
    Tri3,
    Quad4,
    Tet4,
    Hex8,
    Prism6,
    Pyramid5,
}
impl ElementType {
    pub fn dimension(self) -> usize {
        match self {
            Self::Line2 => 1,
            Self::Tri3 | Self::Quad4 => 2,
            _ => 3,
        }
    }
    pub fn node_count(self) -> usize {
        match self {
            Self::Line2 => 2,
            Self::Tri3 => 3,
            Self::Quad4 | Self::Tet4 => 4,
            Self::Hex8 => 8,
            Self::Prism6 => 6,
            Self::Pyramid5 => 5,
        }
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct ReferencePoint {
    pub xi: f64,
    pub eta: f64,
    pub zeta: f64,
}
impl ReferencePoint {
    pub fn new(xi: f64, eta: f64, zeta: f64) -> Self {
        Self { xi, eta, zeta }
    }
}
#[derive(Clone, Debug)]
pub struct ShapeEvaluation {
    pub value: Vec<f64>,
    pub gradient_reference: Vec<Point3>,
}
#[derive(Clone, Copy, Debug)]
pub struct QuadraturePoint {
    pub point: ReferencePoint,
    pub weight: f64,
}
#[derive(Clone, Debug)]
pub struct IsoparametricEvaluation {
    pub point: Point3,
    pub jacobian: [[f64; 3]; 3],
    pub measure: f64,
    pub gradient_physical: Vec<Point3>,
}
pub fn evaluate_shape(t: ElementType, p: ReferencePoint) -> ShapeEvaluation {
    let (x, y, z) = (p.xi, p.eta, p.zeta);
    let mut v = vec![0.0; t.node_count()];
    let mut g = vec![Point3::ZERO; t.node_count()];
    match t {
        ElementType::Line2 => {
            v = vec![0.5 * (1.0 - x), 0.5 * (1.0 + x)];
            g = vec![Point3::new(-0.5, 0.0, 0.0), Point3::new(0.5, 0.0, 0.0)];
        }
        ElementType::Tri3 => {
            v = vec![1.0 - x - y, x, y];
            g = vec![
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ];
        }
        ElementType::Tet4 => {
            v = vec![1.0 - x - y - z, x, y, z];
            g = vec![
                Point3::new(-1.0, -1.0, -1.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
            ];
        }
        ElementType::Quad4 | ElementType::Hex8 | ElementType::Pyramid5 => {
            let sx = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
            let sy = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0];
            let sz = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0];
            let count = if t == ElementType::Hex8 { 8 } else { 4 };
            for i in 0..count {
                let ax = 1.0 + sx[i] * x;
                let ay = 1.0 + sy[i] * y;
                if t == ElementType::Hex8 {
                    let az = 1.0 + sz[i] * z;
                    v[i] = 0.125 * ax * ay * az;
                    g[i] = Point3::new(
                        0.125 * sx[i] * ay * az,
                        0.125 * sy[i] * ax * az,
                        0.125 * sz[i] * ax * ay,
                    );
                } else {
                    let az = if t == ElementType::Pyramid5 {
                        1.0 - z
                    } else {
                        1.0
                    };
                    v[i] = 0.25 * ax * ay * az;
                    g[i] = Point3::new(
                        0.25 * sx[i] * ay * az,
                        0.25 * sy[i] * ax * az,
                        if t == ElementType::Pyramid5 {
                            -0.25 * ax * ay
                        } else {
                            0.0
                        },
                    );
                }
            }
            if t == ElementType::Pyramid5 {
                v[4] = z;
                g[4] = Point3::new(0.0, 0.0, 1.0);
            }
        }
        ElementType::Prism6 => {
            let tv = [1.0 - x - y, x, y];
            let tg = [
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ];
            for i in 0..3 {
                v[i] = tv[i] * 0.5 * (1.0 - z);
                v[i + 3] = tv[i] * 0.5 * (1.0 + z);
                g[i] = Point3::new(
                    tg[i].x * 0.5 * (1.0 - z),
                    tg[i].y * 0.5 * (1.0 - z),
                    -0.5 * tv[i],
                );
                g[i + 3] = Point3::new(
                    tg[i].x * 0.5 * (1.0 + z),
                    tg[i].y * 0.5 * (1.0 + z),
                    0.5 * tv[i],
                );
            }
        }
    }
    ShapeEvaluation {
        value: v,
        gradient_reference: g,
    }
}
pub fn gaussian_quadrature(t: ElementType, order: u32) -> Result<Vec<QuadraturePoint>> {
    require(
        order == 1 || order == 2,
        "quadrature supports orders 1 and 2",
    )?;
    let q = |x, y, z, w| QuadraturePoint {
        point: ReferencePoint::new(x, y, z),
        weight: w,
    };
    if order == 1 {
        return Ok(vec![match t {
            ElementType::Line2 => q(0.0, 0.0, 0.0, 2.0),
            ElementType::Tri3 => q(1.0 / 3.0, 1.0 / 3.0, 0.0, 0.5),
            ElementType::Quad4 => q(0.0, 0.0, 0.0, 4.0),
            ElementType::Tet4 => q(0.25, 0.25, 0.25, 1.0 / 6.0),
            ElementType::Hex8 => q(0.0, 0.0, 0.0, 8.0),
            ElementType::Prism6 => q(1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0),
            ElementType::Pyramid5 => q(0.0, 0.0, 0.25, 4.0 / 3.0),
        }]);
    }
    let g = 1.0 / 3.0_f64.sqrt();
    let mut out = vec![];
    match t {
        ElementType::Line2 => {
            out = vec![q(-g, 0.0, 0.0, 1.0), q(g, 0.0, 0.0, 1.0)];
        }
        ElementType::Tri3 => {
            out = vec![
                q(1.0 / 6.0, 1.0 / 6.0, 0.0, 1.0 / 6.0),
                q(2.0 / 3.0, 1.0 / 6.0, 0.0, 1.0 / 6.0),
                q(1.0 / 6.0, 2.0 / 3.0, 0.0, 1.0 / 6.0),
            ];
        }
        ElementType::Quad4 => {
            for y in [-g, g] {
                for x in [-g, g] {
                    out.push(q(x, y, 0.0, 1.0));
                }
            }
        }
        ElementType::Tet4 => {
            let (a, b, w) = (0.5854101966249685, 0.1381966011250105, 1.0 / 24.0);
            out = vec![q(b, b, b, w), q(a, b, b, w), q(b, a, b, w), q(b, b, a, w)];
        }
        ElementType::Hex8 => {
            for z in [-g, g] {
                for y in [-g, g] {
                    for x in [-g, g] {
                        out.push(q(x, y, z, 1.0));
                    }
                }
            }
        }
        ElementType::Prism6 => {
            for z in [-g, g] {
                for qt in gaussian_quadrature(ElementType::Tri3, 2)? {
                    out.push(q(qt.point.xi, qt.point.eta, z, qt.weight));
                }
            }
        }
        ElementType::Pyramid5 => {
            for z in [0.5 * (1.0 - g), 0.5 * (1.0 + g)] {
                let s = 1.0 - z;
                for y in [-g, g] {
                    for x in [-g, g] {
                        out.push(q(s * x, s * y, z, 0.5 * s * s));
                    }
                }
            }
        }
    }
    Ok(out)
}
pub fn evaluate_isoparametric(
    t: ElementType,
    nodes: &[Point3],
    p: ReferencePoint,
) -> Result<IsoparametricEvaluation> {
    require(
        nodes.len() == t.node_count() && nodes.iter().all(|x| x.is_finite()),
        "isoparametric node count/coordinates",
    )?;
    let s = evaluate_shape(t, p);
    let mut out = IsoparametricEvaluation {
        point: Point3::ZERO,
        jacobian: [[0.0; 3]; 3],
        measure: 0.0,
        gradient_physical: vec![Point3::ZERO; nodes.len()],
    };
    for i in 0..nodes.len() {
        out.point += nodes[i] * s.value[i];
        for a in 0..3 {
            for b in 0..3 {
                out.jacobian[a][b] += nodes[i].component(a) * s.gradient_reference[i].component(b);
            }
        }
    }
    let j = out.jacobian;
    if t.dimension() == 1 {
        let tangent = Point3::new(j[0][0], j[1][0], j[2][0]);
        out.measure = tangent.norm();
        require(out.measure > 1e-30, "singular line mapping")?;
        for i in 0..nodes.len() {
            out.gradient_physical[i] = tangent * (s.gradient_reference[i].x / out.measure.powi(2));
        }
    } else if t.dimension() == 2 {
        let d = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        require(d.abs() > 1e-30, "singular 2D mapping")?;
        out.measure = d.abs();
        for (i, g) in s.gradient_reference.iter().enumerate() {
            out.gradient_physical[i] = Point3::new(
                (j[1][1] * g.x - j[1][0] * g.y) / d,
                (-j[0][1] * g.x + j[0][0] * g.y) / d,
                0.0,
            );
        }
    } else {
        let a = Point3::new(j[0][0], j[1][0], j[2][0]);
        let b = Point3::new(j[0][1], j[1][1], j[2][1]);
        let c = Point3::new(j[0][2], j[1][2], j[2][2]);
        let d = a.dot(b.cross(c));
        require(d.abs() > 1e-30, "singular 3D mapping")?;
        out.measure = d.abs();
        let rx = b.cross(c) / d;
        let ry = c.cross(a) / d;
        let rz = a.cross(b) / d;
        for (i, g) in s.gradient_reference.iter().enumerate() {
            out.gradient_physical[i] = rx * g.x + ry * g.y + rz * g.z;
        }
    }
    require(
        out.point.is_finite()
            && out.measure.is_finite()
            && out.gradient_physical.iter().all(|g| g.is_finite()),
        "nonfinite FEM mapping",
    )?;
    Ok(out)
}
