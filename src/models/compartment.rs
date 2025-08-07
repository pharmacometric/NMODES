use super::{ModelError, OneCompartmentModel, TwoCompartmentModel, ThreeCompartmentModel};
use serde::{Deserialize, Serialize};
use nalgebra::{DVector, DMatrix};
use std::collections::HashMap;

// CORRECTED: Removed `Send` and `Sync` from derive macro
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ModelType {
    OneCompartment,
    TwoCompartment,
    ThreeCompartment,
}

#[derive(Debug, Clone)]
pub struct ModelParameters {
    pub fixed_effects: Vec<f64>,
    pub random_effects_variance: Vec<Vec<f64>>,
    pub residual_variance: f64,
    pub parameter_names: Vec<String>,
}

impl ModelParameters {
    pub fn new(n_params: usize, param_names: Vec<String>) -> Self {
        Self {
            fixed_effects: vec![0.0; n_params],
            random_effects_variance: {
                let mut matrix = vec![vec![0.0; n_params]; n_params];
                for i in 0..n_params {
                    matrix[i][i] = 1.0; // Identity matrix
                }
                matrix
            },
            residual_variance: 1.0,
            parameter_names: param_names,
        }
    }

    pub fn n_parameters(&self) -> usize {
        self.fixed_effects.len()
    }

    pub fn get_parameter(&self, name: &str) -> Option<f64> {
        self.parameter_names.iter()
            .position(|n| n == name)
            .map(|idx| self.fixed_effects[idx])
    }

    pub fn set_parameter(&mut self, name: &str, value: f64) -> Result<(), ModelError> {
        if let Some(idx) = self.parameter_names.iter().position(|n| n == name) {
            if value <= 0.0 {
                return Err(ModelError::BoundsViolation(
                    format!("Parameter {} must be positive, got {}", name, value)
                ));
            }
            self.fixed_effects[idx] = value;
            Ok(())
        } else {
            Err(ModelError::InvalidParameter {
                parameter: name.to_string(),
                value,
            })
        }
    }
    
    pub fn get_fixed_effects_vector(&self) -> DVector<f64> {
        DVector::from_vec(self.fixed_effects.clone())
    }
    
    pub fn get_random_effects_matrix(&self) -> DMatrix<f64> {
        let n = self.random_effects_variance.len();
        let mut matrix = DMatrix::<f64>::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                matrix[(i, j)] = self.random_effects_variance[i][j];
            }
        }
        matrix
    }
    
    pub fn set_fixed_effects(&mut self, effects: &DVector<f64>) {
        // Apply bounds checking - all PK parameters must be positive after exp transformation
        self.fixed_effects = effects.as_slice().iter()
            .map(|&x| x.max(-10.0)) // Prevent exp(x) from being too small (exp(-10) ≈ 4.5e-5)
            .collect();
    }
    
    pub fn set_random_effects_variance(&mut self, variance: &DMatrix<f64>) {
        let n = variance.nrows();
        self.random_effects_variance = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                self.random_effects_variance[i][j] = variance[(i, j)];
            }
        }
    }
}

// CORRECTED: Removed `Send` and `Sync` from derive macro
#[derive(Debug, Clone)]
pub struct ModelState {
    pub compartments: DVector<f64>,
    pub time: f64,
}

impl ModelState {
    pub fn new(n_compartments: usize) -> Self {
        Self {
            compartments: DVector::<f64>::zeros(n_compartments),
            time: 0.0,
        }
    }

    pub fn with_time(mut self, time: f64) -> Self {
        self.time = time;
        self
    }

    pub fn get_concentration(&self, compartment: usize) -> f64 {
        if compartment > 0 && compartment <= self.compartments.nrows() {
            self.compartments[compartment - 1]
        } else {
            0.0
        }
    }

    pub fn add_dose(&mut self, compartment: usize, amount: f64) {
        if compartment > 0 && compartment <= self.compartments.nrows() {
            self.compartments[compartment - 1] += amount;
        }
    }
}

pub trait CompartmentModelTrait {
    fn n_compartments(&self) -> usize;
    fn parameter_names(&self) -> Vec<String>;
    fn default_parameters(&self) -> ModelParameters;
    fn derivatives(&self, state: &ModelState, params: &ModelParameters) -> DVector<f64>;
    fn observation_function(&self, state: &ModelState, params: &ModelParameters, compartment: usize) -> f64;
    fn validate_parameters(&self, params: &ModelParameters) -> Result<(), ModelError>;
}

pub struct CompartmentModel {
    model_type: ModelType,
    inner: Box<dyn CompartmentModelTrait + Send + Sync>,
}

impl CompartmentModel {
    pub fn new(model_type: ModelType) -> Result<Self, ModelError> {
        let inner: Box<dyn CompartmentModelTrait + Send + Sync> = match model_type {
            ModelType::OneCompartment => Box::new(OneCompartmentModel::new()),
            ModelType::TwoCompartment => Box::new(TwoCompartmentModel::new()),
            ModelType::ThreeCompartment => Box::new(ThreeCompartmentModel::new()),
        };

        Ok(Self {
            model_type,
            inner,
        })
    }

    pub fn model_type(&self) -> &ModelType {
        &self.model_type
    }

    pub fn n_compartments(&self) -> usize {
        self.inner.n_compartments()
    }

    pub fn parameter_names(&self) -> Vec<String> {
        self.inner.parameter_names()
    }

    pub fn default_parameters(&self) -> ModelParameters {
        self.inner.default_parameters()
    }

    pub fn derivatives(&self, state: &ModelState, params: &ModelParameters) -> DVector<f64> {
        self.inner.derivatives(state, params)
    }

    pub fn observation_function(&self, state: &ModelState, params: &ModelParameters, compartment: usize) -> f64 {
        self.inner.observation_function(state, params, compartment)
    }

    pub fn validate_parameters(&self, params: &ModelParameters) -> Result<(), ModelError> {
        self.inner.validate_parameters(params)
    }
}

// Note: These unsafe impls are likely here because of the trait object `inner`.
// They are correct and should remain. The error was with the `derive` macro.
unsafe impl Send for CompartmentModel {}
unsafe impl Sync for CompartmentModel {}

impl std::fmt::Debug for CompartmentModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CompartmentModel")
            .field("model_type", &self.model_type)
            .field("n_compartments", &self.n_compartments())
            .finish()
    }
}