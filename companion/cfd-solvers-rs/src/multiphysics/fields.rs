//! Named, unit-tagged field storage shared between coupled solvers.
use crate::{require, Result, SolverError};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldLocation {
    Cell,
    Face,
    Node,
    Edge,
    GridPoint,
    Surface,
    Global,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldTopology {
    Structured,
    Unstructured,
    Optical,
    MeshIndependent,
}

/// SI base-dimension exponents (kg, m, s, A, K, mol, cd), a scale to SI and a symbol.
#[derive(Clone, Debug, PartialEq)]
pub struct UnitSignature {
    pub exponent: [i32; 7],
    pub scale_to_si: f64,
    pub symbol: String,
}

impl Default for UnitSignature {
    fn default() -> Self {
        Self {
            exponent: [0; 7],
            scale_to_si: 1.0,
            symbol: "1".into(),
        }
    }
}

impl UnitSignature {
    pub fn dimensionally_compatible(&self, other: &Self) -> bool {
        self.exponent == other.exponent
    }
    pub fn dimensionless() -> Self {
        Self::default()
    }
    pub fn kelvin() -> Self {
        Self {
            exponent: [0, 0, 0, 0, 1, 0, 0],
            scale_to_si: 1.0,
            symbol: "K".into(),
        }
    }
    pub fn metre() -> Self {
        Self {
            exponent: [0, 1, 0, 0, 0, 0, 0],
            scale_to_si: 1.0,
            symbol: "m".into(),
        }
    }
    /// W/m^3 = kg m^-1 s^-3.
    pub fn watt_per_cubic_metre() -> Self {
        Self {
            exponent: [1, -1, -3, 0, 0, 0, 0],
            scale_to_si: 1.0,
            symbol: "W/m^3".into(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FieldMetadata {
    pub name: String,
    pub components: usize,
    pub entities: usize,
    pub location: FieldLocation,
    pub topology: FieldTopology,
    pub units: UnitSignature,
    pub producer: String,
}

impl FieldMetadata {
    pub fn new(name: &str, components: usize, entities: usize) -> Self {
        Self {
            name: name.into(),
            components,
            entities,
            location: FieldLocation::Cell,
            topology: FieldTopology::Unstructured,
            units: UnitSignature::default(),
            producer: String::new(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct FieldRegistry {
    fields: BTreeMap<String, (FieldMetadata, Vec<f64>)>,
}

impl FieldRegistry {
    pub fn add(&mut self, metadata: FieldMetadata, initial_value: f64) -> Result<()> {
        require(
            !metadata.name.is_empty() && metadata.components > 0 && metadata.entities > 0,
            "invalid field metadata",
        )?;
        require(
            metadata.units.scale_to_si.is_finite() && metadata.units.scale_to_si > 0.0,
            "field unit scale must be positive",
        )?;
        require(
            !self.fields.contains_key(&metadata.name),
            format!("field already exists: {}", metadata.name),
        )?;
        let size = metadata
            .components
            .checked_mul(metadata.entities)
            .ok_or_else(|| SolverError::InvalidInput("field size overflow".into()))?;
        self.fields
            .insert(metadata.name.clone(), (metadata, vec![initial_value; size]));
        Ok(())
    }
    pub fn contains(&self, name: &str) -> bool {
        self.fields.contains_key(name)
    }
    pub fn metadata(&self, name: &str) -> Result<&FieldMetadata> {
        self.fields
            .get(name)
            .map(|f| &f.0)
            .ok_or_else(|| SolverError::InvalidInput(format!("unknown field: {name}")))
    }
    pub fn data(&self, name: &str) -> Result<&[f64]> {
        self.fields
            .get(name)
            .map(|f| f.1.as_slice())
            .ok_or_else(|| SolverError::InvalidInput(format!("unknown field: {name}")))
    }
    pub fn data_mut(&mut self, name: &str) -> Result<&mut [f64]> {
        self.fields
            .get_mut(name)
            .map(|f| f.1.as_mut_slice())
            .ok_or_else(|| SolverError::InvalidInput(format!("unknown field: {name}")))
    }
    /// Field names in sorted order.
    pub fn names(&self) -> Vec<String> {
        self.fields.keys().cloned().collect()
    }
}
