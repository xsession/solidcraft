//! One-way electro-thermal coupling on a shared Tri3 mesh:
//! DC conduction -> element current density -> Joule heat `|J|^2/sigma` -> transient heat FEM.
use crate::{
    core::LinearConfig,
    fem::{DcConduction2D, Heat2D, Heat2DConfig, Mesh2D, Node2, ScalarBoundaryCondition2D},
    require, Result, SolverError,
};

pub struct JouleHeatingCoupler2D {
    electrical: DcConduction2D,
    thermal: Heat2D,
    joule_heating_density: Vec<f64>,
    electrical_solved: bool,
    boundary_temperature: Option<Box<dyn Fn(Node2, f64) -> f64>>,
}

impl JouleHeatingCoupler2D {
    pub fn new(mesh: Mesh2D, thermal: Heat2DConfig, electrical: LinearConfig) -> Result<Self> {
        let n = mesh.element_count();
        Ok(Self {
            electrical: DcConduction2D::new(mesh.clone(), electrical)?,
            thermal: Heat2D::new(mesh, thermal)?,
            joule_heating_density: vec![0.0; n],
            electrical_solved: false,
            boundary_temperature: None,
        })
    }
    pub fn set_electrical_boundary(
        &mut self,
        patch: i32,
        condition: ScalarBoundaryCondition2D<'static>,
    ) {
        self.electrical.set_boundary(patch, condition);
        self.electrical_solved = false;
    }
    pub fn initialize_temperature(&mut self, temperature: impl Fn(Node2) -> f64) -> Result<()> {
        self.thermal.initialize(temperature)
    }
    pub fn set_thermal_dirichlet(
        &mut self,
        boundary_temperature: impl Fn(Node2, f64) -> f64 + 'static,
    ) {
        self.boundary_temperature = Some(Box::new(boundary_temperature));
    }
    pub fn solve_electrical(
        &mut self,
        conductivity: &dyn Fn(Node2) -> f64,
        current_source: Option<&dyn Fn(Node2) -> f64>,
    ) -> Result<()> {
        self.electrical_solved = false;
        self.electrical.solve(conductivity, current_source)?;
        let current = self.electrical.element_current_density(conductivity)?;
        let mesh = self.electrical.mesh();
        require(
            current.len() == mesh.element_count(),
            "electrical element field size mismatch",
        )?;
        self.joule_heating_density.clear();
        for (tri, j) in mesh.triangles.iter().zip(&current) {
            let [a, b, c] = tri.node.map(|i| mesh.nodes[i]);
            let sigma = conductivity(Node2 {
                x: (a.x + b.x + c.x) / 3.0,
                y: (a.y + b.y + c.y) / 3.0,
            });
            require(
                sigma.is_finite() && sigma > 0.0,
                "invalid conductivity in Joule-heating transfer",
            )?;
            self.joule_heating_density
                .push((j.x * j.x + j.y * j.y) / sigma);
        }
        self.electrical_solved = true;
        Ok(())
    }
    pub fn thermal_step(&mut self) -> Result<()> {
        if !self.electrical_solved
            || self.joule_heating_density.len() != self.thermal.mesh.element_count()
        {
            return Err(SolverError::InvalidInput(
                "Joule-heating field is not initialized".into(),
            ));
        }
        let boundary = self.boundary_temperature.as_ref().ok_or_else(|| {
            SolverError::InvalidInput("thermal Dirichlet boundary is not set".into())
        })?;
        self.thermal
            .step_element_source(&self.joule_heating_density, boundary)
    }
    pub fn thermal_run(&mut self, steps: usize) -> Result<()> {
        for _ in 0..steps {
            self.thermal_step()?;
        }
        Ok(())
    }
    pub fn electrical(&self) -> &DcConduction2D {
        &self.electrical
    }
    pub fn thermal(&self) -> &Heat2D {
        &self.thermal
    }
    pub fn joule_heating_density(&self) -> &[f64] {
        &self.joule_heating_density
    }
}
