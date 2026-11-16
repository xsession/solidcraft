use crate::{core::Vec3, positive, require, Result};
use serde::{Deserialize, Serialize};

pub type Point3 = Vec3;
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
pub struct Node2 {
    pub x: f64,
    pub y: f64,
}
impl Node2 {
    pub fn point3(self) -> Vec3 {
        Vec3::new(self.x, self.y, 0.0)
    }
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Tri3 {
    pub node: [usize; 3],
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct BoundaryEdge2 {
    pub node: [usize; 2],
    pub patch: i32,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Mesh2D {
    pub nodes: Vec<Node2>,
    pub triangles: Vec<Tri3>,
    pub boundary_edges: Vec<BoundaryEdge2>,
    pub boundary_node: Vec<bool>,
}
impl Mesh2D {
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    pub fn element_count(&self) -> usize {
        self.triangles.len()
    }
    pub fn validate(&self) -> Result<()> {
        require(
            !self.nodes.is_empty()
                && !self.triangles.is_empty()
                && self.boundary_node.len() == self.nodes.len(),
            "FEM mesh shape",
        )?;
        require(
            self.nodes
                .iter()
                .all(|p| p.x.is_finite() && p.y.is_finite()),
            "nonfinite mesh coordinate",
        )?;
        for e in &self.boundary_edges {
            require(
                e.node.iter().all(|&i| i < self.nodes.len()) && e.node[0] != e.node[1],
                "invalid boundary edge",
            )?;
        }
        for t in &self.triangles {
            require(
                t.node.iter().all(|&i| i < self.nodes.len()),
                "triangle index out of range",
            )?;
            let [a, b, c] = t.node.map(|i| self.nodes[i]);
            require(
                ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)).abs() > 1e-30,
                "degenerate triangle",
            )?;
        }
        Ok(())
    }
}
pub fn make_rectangle_tri_mesh(nx: usize, ny: usize, width: f64, height: f64) -> Result<Mesh2D> {
    require(nx > 0 && ny > 0, "rectangle subdivisions must be positive")?;
    positive(width, "width")?;
    positive(height, "height")?;
    let sx = nx
        .checked_add(1)
        .ok_or_else(|| crate::SolverError::InvalidInput("mesh size overflow".into()))?;
    let count = sx
        .checked_mul(
            ny.checked_add(1)
                .ok_or_else(|| crate::SolverError::InvalidInput("mesh size overflow".into()))?,
        )
        .ok_or_else(|| crate::SolverError::InvalidInput("mesh size overflow".into()))?;
    let mut mesh = Mesh2D {
        nodes: vec![Node2::default(); count],
        triangles: vec![],
        boundary_edges: vec![],
        boundary_node: vec![false; count],
    };
    let id = |i, j| j * sx + i;
    for j in 0..=ny {
        for i in 0..=nx {
            let n = id(i, j);
            mesh.nodes[n] = Node2 {
                x: width * i as f64 / nx as f64,
                y: height * j as f64 / ny as f64,
            };
            mesh.boundary_node[n] = i == 0 || i == nx || j == 0 || j == ny;
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            mesh.triangles.push(Tri3 {
                node: [id(i, j), id(i + 1, j), id(i + 1, j + 1)],
            });
            mesh.triangles.push(Tri3 {
                node: [id(i, j), id(i + 1, j + 1), id(i, j + 1)],
            });
        }
    }
    for j in 0..ny {
        mesh.boundary_edges.push(BoundaryEdge2 {
            node: [id(0, j + 1), id(0, j)],
            patch: 0,
        });
    }
    for j in 0..ny {
        mesh.boundary_edges.push(BoundaryEdge2 {
            node: [id(nx, j), id(nx, j + 1)],
            patch: 1,
        });
    }
    for i in 0..nx {
        mesh.boundary_edges.push(BoundaryEdge2 {
            node: [id(i, 0), id(i + 1, 0)],
            patch: 2,
        });
    }
    for i in 0..nx {
        mesh.boundary_edges.push(BoundaryEdge2 {
            node: [id(i + 1, ny), id(i, ny)],
            patch: 3,
        });
    }
    mesh.validate()?;
    Ok(mesh)
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Tet4 {
    pub node: [usize; 4],
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Mesh3D {
    pub nodes: Vec<Point3>,
    pub tetrahedra: Vec<Tet4>,
    pub boundary_node: Vec<bool>,
}
impl Mesh3D {
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    pub fn element_count(&self) -> usize {
        self.tetrahedra.len()
    }
    pub fn validate(&self) -> Result<()> {
        require(
            !self.nodes.is_empty()
                && !self.tetrahedra.is_empty()
                && self.boundary_node.len() == self.nodes.len(),
            "tetrahedral mesh shape",
        )?;
        require(
            self.nodes.iter().all(|p| p.is_finite()),
            "nonfinite mesh coordinate",
        )?;
        for t in &self.tetrahedra {
            require(
                t.node.iter().all(|&i| i < self.nodes.len()),
                "tetrahedron index out of range",
            )?;
            let [a, b, c, d] = t.node.map(|i| self.nodes[i]);
            require(
                (b - a).dot((c - a).cross(d - a)).abs() > 1e-30,
                "degenerate tetrahedron",
            )?;
        }
        Ok(())
    }
}
pub fn make_box_tet_mesh(
    nx: usize,
    ny: usize,
    nz: usize,
    width: f64,
    height: f64,
    depth: f64,
) -> Result<Mesh3D> {
    require(
        nx > 0 && ny > 0 && nz > 0,
        "box subdivisions must be positive",
    )?;
    positive(width, "width")?;
    positive(height, "height")?;
    positive(depth, "depth")?;
    let checked = || {
        nx.checked_add(1)?
            .checked_mul(ny.checked_add(1)?)?
            .checked_mul(nz.checked_add(1)?)
    };
    let count =
        checked().ok_or_else(|| crate::SolverError::InvalidInput("mesh size overflow".into()))?;
    let id = |i, j, k| (k * (ny + 1) + j) * (nx + 1) + i;
    let mut m = Mesh3D {
        nodes: vec![Vec3::ZERO; count],
        tetrahedra: vec![],
        boundary_node: vec![false; count],
    };
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                let n = id(i, j, k);
                m.nodes[n] = Vec3::new(
                    width * i as f64 / nx as f64,
                    height * j as f64 / ny as f64,
                    depth * k as f64 / nz as f64,
                );
                m.boundary_node[n] = i == 0 || i == nx || j == 0 || j == ny || k == 0 || k == nz;
            }
        }
    }
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let a = id(i, j, k);
                let z = id(i + 1, j + 1, k + 1);
                let ring = [
                    id(i + 1, j, k),
                    id(i + 1, j + 1, k),
                    id(i, j + 1, k),
                    id(i, j + 1, k + 1),
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                ];
                for r in 0..6 {
                    m.tetrahedra.push(Tet4 {
                        node: [a, ring[r], ring[(r + 1) % 6], z],
                    });
                }
            }
        }
    }
    m.validate()?;
    Ok(m)
}
