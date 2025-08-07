pub mod compartment;
pub mod one_compartment;
pub mod two_compartment;
pub mod three_compartment;

pub use compartment::{CompartmentModel, ModelType, ModelParameters, ModelState};
pub use one_compartment::OneCompartmentModel;
pub use two_compartment::TwoCompartmentModel;
pub use three_compartment::ThreeCompartmentModel;

use thiserror::Error;

impl std::fmt::Display for ModelType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ModelType::OneCompartment => write!(f, "one-compartment"),
            ModelType::TwoCompartment => write!(f, "two-compartment"),
            ModelType::ThreeCompartment => write!(f, "three-compartment"),
        }
    }
}

#[derive(Error, Debug)]
pub enum ModelError {
    #[error("Invalid parameter value: {parameter} = {value}")]
    InvalidParameter { parameter: String, value: f64 },
    
    #[error("Model integration failed: {0}")]
    IntegrationFailed(String),
    
    #[error("Unsupported model type: {0}")]
    UnsupportedModel(String),
    
    #[error("Parameter bounds violation: {0}")]
    BoundsViolation(String),
}