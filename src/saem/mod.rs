pub mod algorithm;
pub mod mcmc;

pub use algorithm::SaemEstimator;
pub use mcmc::{McmcSampler, McmcConfig};

use nalgebra::{DVector, DMatrix};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParameterStatistics {
    pub name: String,
    pub estimate: f64,
    pub rse_percent: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OmegaStatistics {
    pub parameter_i: String,
    pub parameter_j: String,
    pub estimate: f64,
    pub shrinkage_percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SaemResults {
    pub fixed_effects: Vec<f64>,
    pub random_effects_variance: Vec<Vec<f64>>,
    pub residual_variance: f64,
    pub log_likelihood_trajectory: Vec<f64>,
    pub parameter_trajectory: Vec<Vec<f64>>,
    pub final_log_likelihood: f64,
    pub objective_function_value: f64,
    pub converged: bool,
    pub n_iterations: usize,
    pub individual_parameters: HashMap<i32, Vec<f64>>,
    pub parameter_statistics: Vec<ParameterStatistics>,
    pub omega_statistics: Vec<OmegaStatistics>,
    pub parameter_names: Vec<String>,
}

impl SaemResults {
    pub fn new(n_params: usize, parameter_names: Vec<String>) -> Self {
        Self {
            fixed_effects: vec![0.0; n_params],
            random_effects_variance: vec![vec![0.0; n_params]; n_params],
            residual_variance: 1.0,
            log_likelihood_trajectory: Vec::new(),
            parameter_trajectory: Vec::new(),
            final_log_likelihood: f64::NEG_INFINITY,
            objective_function_value: f64::INFINITY,
            converged: false,
            n_iterations: 0,
            individual_parameters: HashMap::new(),
            parameter_statistics: Vec::new(),
            omega_statistics: Vec::new(),
            parameter_names,
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
        self.fixed_effects = effects.as_slice().to_vec();
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