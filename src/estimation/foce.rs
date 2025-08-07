use crate::data::{Dataset, Individual};
use crate::models::{CompartmentModel, ModelParameters, ModelState};
use crate::solver::{OdeSolver, OdeSystem, RungeKuttaSolver, SolverConfig};
use super::EstimationConfig;
use anyhow::{Context, Result};
use log::{info, debug, warn};
use nalgebra::{DVector, DMatrix};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FoceResults {
    pub fixed_effects: Vec<f64>,
    pub random_effects_variance: Vec<Vec<f64>>,
    pub residual_variance: f64,
    pub objective_function_value: f64,
    pub final_log_likelihood: f64,
    pub converged: bool,
    pub n_iterations: usize,
    pub individual_parameters: HashMap<i32, Vec<f64>>,
    pub parameter_names: Vec<String>,
    pub gradient_norm: f64,
    pub hessian_condition_number: f64,
    pub covariance_matrix: Vec<Vec<f64>>,
    pub standard_errors: Vec<f64>,
}

impl FoceResults {
    pub fn new(n_params: usize, parameter_names: Vec<String>) -> Self {
        Self {
            fixed_effects: vec![0.0; n_params],
            random_effects_variance: vec![vec![0.0; n_params]; n_params],
            residual_variance: 1.0,
            objective_function_value: f64::INFINITY,
            final_log_likelihood: f64::NEG_INFINITY,
            converged: false,
            n_iterations: 0,
            individual_parameters: HashMap::new(),
            parameter_names,
            gradient_norm: f64::INFINITY,
            hessian_condition_number: f64::INFINITY,
            covariance_matrix: vec![vec![0.0; n_params]; n_params],
            standard_errors: vec![0.0; n_params],
        }
    }
}

pub struct FoceEstimator {
    model: CompartmentModel,
    config: EstimationConfig,
    solver: Box<dyn OdeSolver + Send + Sync>,
}

impl FoceEstimator {
    pub fn new(model: CompartmentModel, config: EstimationConfig) -> Self {
        let solver = Box::new(RungeKuttaSolver::new());
        
        Self {
            model,
            config,
            solver,
        }
    }

    pub fn model(&self) -> &CompartmentModel {
        &self.model
    }

    pub fn fit(&mut self, dataset: &Dataset) -> Result<FoceResults> {
        info!("Starting FOCE estimation for {} individuals", dataset.n_individuals());
        
        let n_params = self.model.parameter_names().len();
        let parameter_names = self.model.parameter_names();
        let mut results = FoceResults::new(n_params, parameter_names);
        
        // Initialize parameters
        let mut current_params = self.model.default_parameters();
        let mut individual_params: HashMap<i32, Vec<f64>> = HashMap::new();
        
        // Initialize individual parameters to population means
        for (&id, _) in dataset.individuals() {
            individual_params.insert(id, current_params.fixed_effects.clone());
        }

        let mut previous_objective = f64::INFINITY;
        
        for iteration in 0..self.config.foce_max_iterations {
            debug!("FOCE iteration {}/{}", iteration + 1, self.config.foce_max_iterations);
            
            // E-step: Estimate individual parameters using first-order approximation
            self.estimate_individual_parameters(dataset, &current_params, &mut individual_params)?;
            
            // M-step: Update population parameters
            let objective = self.update_population_parameters(
                dataset,
                &individual_params,
                &mut current_params,
            )?;
            
            // Check convergence
            let objective_change = (previous_objective - objective).abs();
            let relative_change = objective_change / previous_objective.abs();
            
            if relative_change < self.config.foce_tolerance {
                info!("FOCE converged at iteration {} (relative change: {:.2e})", 
                      iteration + 1, relative_change);
                results.converged = true;
                break;
            }
            
            if iteration % 10 == 0 {
                info!("FOCE iteration {}: Objective = {:.3}, Change = {:.2e}", 
                      iteration + 1, objective, relative_change);
            }
            
            previous_objective = objective;
        }

        // Calculate final statistics
        let final_objective = self.calculate_objective_function(dataset, &individual_params, &current_params)?;
        
        // Estimate covariance matrix and standard errors
        let (covariance_matrix, standard_errors) = self.estimate_covariance_matrix(
            dataset, &individual_params, &current_params
        )?;

        // Populate results
        results.fixed_effects = current_params.fixed_effects;
        results.random_effects_variance = current_params.random_effects_variance;
        results.residual_variance = current_params.residual_variance;
        results.objective_function_value = final_objective;
        results.final_log_likelihood = -final_objective / 2.0;
        results.n_iterations = self.config.foce_max_iterations;
        results.individual_parameters = individual_params;
        results.covariance_matrix = covariance_matrix;
        results.standard_errors = standard_errors;

        info!("FOCE estimation completed. Objective function: {:.3}, Converged: {}", 
              results.objective_function_value, results.converged);

        Ok(results)
    }

    fn estimate_individual_parameters(
        &self,
        dataset: &Dataset,
        population_params: &ModelParameters,
        individual_params: &mut HashMap<i32, Vec<f64>>,
    ) -> Result<()> {
        for (&id, individual) in dataset.individuals() {
            let current_eta = individual_params.get(&id).unwrap().clone();
            
            // Newton-Raphson optimization for individual parameters
            let optimized_eta = self.optimize_individual_eta(
                individual,
                population_params,
                &current_eta,
            )?;
            
            individual_params.insert(id, optimized_eta);
        }
        
        Ok(())
    }

    fn optimize_individual_eta(
        &self,
        individual: &Individual,
        population_params: &ModelParameters,
        initial_eta: &[f64],
    ) -> Result<Vec<f64>> {
        let mut eta = initial_eta.to_vec();
        let max_inner_iterations = 20;
        
        for _iter in 0..max_inner_iterations {
            // Calculate gradient and Hessian of individual objective function
            let (gradient, hessian) = self.calculate_individual_derivatives(
                individual,
                population_params,
                &eta,
            )?;
            
            // Newton-Raphson step: eta_new = eta - H^(-1) * g
            let hessian_matrix = DMatrix::from_vec(eta.len(), eta.len(), hessian);
            let gradient_vector = DVector::from_vec(gradient);
            
            // Check if Hessian is positive definite (add regularization if needed)
            let regularized_hessian = self.regularize_hessian(&hessian_matrix);
            
            if let Some(chol) = regularized_hessian.cholesky() {
                let step = chol.solve(&gradient_vector);
                
                // Update eta with step size control
                let step_size = 1.0; // Could be adaptive
                for i in 0..eta.len() {
                    eta[i] -= step_size * step[i];
                    
                    // Apply bounds: keep individual deviations reasonable
                    eta[i] = eta[i].max(-5.0).min(5.0);
                }
                
                // Check convergence
                if gradient_vector.norm() < 1e-6 {
                    break;
                }
            } else {
                warn!("Hessian not positive definite for individual optimization");
                break;
            }
        }
        
        Ok(eta)
    }

    fn calculate_individual_derivatives(
        &self,
        individual: &Individual,
        population_params: &ModelParameters,
        eta: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let n_params = eta.len();
        let mut gradient = vec![0.0; n_params];
        let mut hessian = vec![0.0; n_params * n_params];
        
        // Calculate individual parameters: theta_i = theta + eta_i
        let mut individual_params = population_params.clone();
        for i in 0..n_params {
            individual_params.fixed_effects[i] = population_params.fixed_effects[i] + eta[i];
        }
        
        // Get predictions and residuals
        let predictions = self.predict_individual(individual, &individual_params)?;
        let mut residuals = Vec::new();
        
        for (obs, pred) in individual.observations().iter().zip(predictions.iter()) {
            residuals.push(obs.value - pred);
        }
        
        // Calculate derivatives using finite differences
        let h = 1e-6;
        
        for i in 0..n_params {
            // Forward difference for gradient
            let mut eta_plus = eta.to_vec();
            eta_plus[i] += h;
            
            let mut params_plus = population_params.clone();
            for j in 0..n_params {
                params_plus.fixed_effects[j] = population_params.fixed_effects[j] + eta_plus[j];
            }
            
            let predictions_plus = self.predict_individual(individual, &params_plus)?;
            
            // Gradient contribution from data likelihood
            let mut grad_data = 0.0;
            for (k, (obs, (pred, pred_plus))) in individual.observations().iter()
                .zip(predictions.iter().zip(predictions_plus.iter()))
                .enumerate()
            {
                let residual = obs.value - pred;
                let dpred_deta = (pred_plus - pred) / h;
                grad_data += residual * dpred_deta / population_params.residual_variance;
            }
            
            // Gradient contribution from prior (eta ~ N(0, Omega))
            let grad_prior = -eta[i] / population_params.random_effects_variance[i][i];
            
            gradient[i] = grad_data + grad_prior;
            
            // Diagonal Hessian approximation
            let mut hess_data = 0.0;
            for (pred, pred_plus) in predictions.iter().zip(predictions_plus.iter()) {
                let dpred_deta = (pred_plus - pred) / h;
                hess_data -= (dpred_deta * dpred_deta) / population_params.residual_variance;
            }
            
            let hess_prior = -1.0 / population_params.random_effects_variance[i][i];
            hessian[i * n_params + i] = hess_data + hess_prior;
        }
        
        Ok((gradient, hessian))
    }

    fn regularize_hessian(&self, hessian: &DMatrix<f64>) -> DMatrix<f64> {
        let mut regularized = hessian.clone();
        let regularization = 1e-6;
        
        // Add regularization to diagonal
        for i in 0..regularized.nrows() {
            regularized[(i, i)] += regularization;
        }
        
        regularized
    }

    fn update_population_parameters(
        &self,
        dataset: &Dataset,
        individual_params: &HashMap<i32, Vec<f64>>,
        current_params: &mut ModelParameters,
    ) -> Result<f64> {
        let n_individuals = individual_params.len() as f64;
        let n_params = current_params.n_parameters();
        
        // Update fixed effects (population means)
        let mut new_fixed_effects = vec![0.0; n_params];
        for params in individual_params.values() {
            for i in 0..n_params {
                new_fixed_effects[i] += params[i];
            }
        }
        for i in 0..n_params {
            new_fixed_effects[i] /= n_individuals;
            // Apply bounds to population parameters
            new_fixed_effects[i] = new_fixed_effects[i].max(-10.0);
        }
        current_params.fixed_effects = new_fixed_effects;
        
        // Update random effects variance (Omega matrix)
        let mut new_omega = vec![vec![0.0; n_params]; n_params];
        for params in individual_params.values() {
            for i in 0..n_params {
                for j in 0..n_params {
                    let eta_i = params[i] - current_params.fixed_effects[i];
                    let eta_j = params[j] - current_params.fixed_effects[j];
                    new_omega[i][j] += eta_i * eta_j;
                }
            }
        }
        for i in 0..n_params {
            for j in 0..n_params {
                new_omega[i][j] /= n_individuals;
            }
        }
        current_params.random_effects_variance = new_omega;
        
        // Update residual variance
        let mut residual_sum = 0.0;
        let mut total_observations = 0;
        
        for (&id, individual) in dataset.individuals() {
            if let Some(ind_params) = individual_params.get(&id) {
                let mut temp_params = current_params.clone();
                for i in 0..n_params {
                    temp_params.fixed_effects[i] = current_params.fixed_effects[i] + ind_params[i];
                }
                
                let predictions = self.predict_individual(individual, &temp_params)?;
                
                for (obs, pred) in individual.observations().iter().zip(predictions.iter()) {
                    let residual = (obs.value - pred).powi(2);
                    residual_sum += residual;
                    total_observations += 1;
                }
            }
        }
        
        if total_observations > 0 {
            current_params.residual_variance = residual_sum / total_observations as f64;
        }
        
        // Calculate objective function
        self.calculate_objective_function(dataset, individual_params, current_params)
    }

    fn calculate_objective_function(
        &self,
        dataset: &Dataset,
        individual_params: &HashMap<i32, Vec<f64>>,
        population_params: &ModelParameters,
    ) -> Result<f64> {
        let mut objective = 0.0;
        
        for (&id, individual) in dataset.individuals() {
            if let Some(eta) = individual_params.get(&id) {
                // Individual parameters: theta_i = theta + eta_i
                let mut ind_params = population_params.clone();
                for i in 0..eta.len() {
                    ind_params.fixed_effects[i] = population_params.fixed_effects[i] + eta[i];
                }
                
                // Data likelihood contribution
                let predictions = self.predict_individual(individual, &ind_params)?;
                for (obs, pred) in individual.observations().iter().zip(predictions.iter()) {
                    let residual = obs.value - pred;
                    objective += (residual * residual) / population_params.residual_variance;
                    objective += (2.0 * std::f64::consts::PI * population_params.residual_variance).ln();
                }
                
                // Prior likelihood contribution (eta ~ N(0, Omega))
                for i in 0..eta.len() {
                    objective += (eta[i] * eta[i]) / population_params.random_effects_variance[i][i];
                    objective += (2.0 * std::f64::consts::PI * population_params.random_effects_variance[i][i]).ln();
                }
            }
        }
        
        Ok(objective)
    }

    fn predict_individual(
        &self,
        individual: &Individual,
        params: &ModelParameters,
    ) -> Result<Vec<f64>> {
        let mut predictions = Vec::new();
        let solver_config = SolverConfig::default();
        
        let system = CompartmentSystem {
            model: &self.model,
            params,
        };
        
        let mut current_state = ModelState::new(self.model.n_compartments());
        let mut last_time = 0.0;
        
        // Apply dosing events
        for dose in individual.dosing_records() {
            if dose.time > last_time {
                let final_state = self.solver.solve_to_time(
                    &system,
                    last_time,
                    dose.time,
                    &current_state.compartments,
                    &solver_config,
                )?;
                current_state.compartments = final_state;
                current_state.time = dose.time;
            }
            
            current_state.add_dose(dose.compartment as usize, dose.amount);
            last_time = dose.time;
        }
        
        // Predict concentrations at observation times
        for obs in individual.observations() {
            if obs.time > last_time {
                let final_state = self.solver.solve_to_time(
                    &system,
                    last_time,
                    obs.time,
                    &current_state.compartments,
                    &solver_config,
                )?;
                current_state.compartments = final_state;
                current_state.time = obs.time;
                last_time = obs.time;
            }
            
            let concentration = self.model.observation_function(
                &current_state,
                params,
                obs.compartment as usize,
            );
            predictions.push(concentration);
        }
        
        Ok(predictions)
    }

    fn estimate_covariance_matrix(
        &self,
        dataset: &Dataset,
        individual_params: &HashMap<i32, Vec<f64>>,
        population_params: &ModelParameters,
    ) -> Result<(Vec<Vec<f64>>, Vec<f64>)> {
        let n_params = population_params.n_parameters();
        
        // Calculate Fisher Information Matrix using finite differences
        let mut fisher_matrix = vec![vec![0.0; n_params]; n_params];
        let h = 1e-6;
        
        for i in 0..n_params {
            for j in 0..n_params {
                // Calculate second derivatives
                let mut params_ij = population_params.clone();
                let mut params_i = population_params.clone();
                let mut params_j = population_params.clone();
                let mut params_base = population_params.clone();
                
                params_ij.fixed_effects[i] += h;
                params_ij.fixed_effects[j] += h;
                params_i.fixed_effects[i] += h;
                params_j.fixed_effects[j] += h;
                
                let obj_ij = self.calculate_objective_function(dataset, individual_params, &params_ij)?;
                let obj_i = self.calculate_objective_function(dataset, individual_params, &params_i)?;
                let obj_j = self.calculate_objective_function(dataset, individual_params, &params_j)?;
                let obj_base = self.calculate_objective_function(dataset, individual_params, &params_base)?;
                
                // Second derivative approximation
                let second_deriv = (obj_ij - obj_i - obj_j + obj_base) / (h * h);
                fisher_matrix[i][j] = second_deriv;
            }
        }
        
        // Invert Fisher matrix to get covariance matrix
        let fisher_dmatrix = DMatrix::from_vec(n_params, n_params, 
            fisher_matrix.iter().flatten().cloned().collect());
        
        let covariance_dmatrix = if let Some(inv) = fisher_dmatrix.clone().try_inverse() {
            inv
        } else {
            warn!("Fisher matrix not invertible, using regularized version");
            let regularized = &fisher_dmatrix + DMatrix::identity(n_params, n_params) * 1e-6;
            regularized.try_inverse().unwrap_or_else(|| DMatrix::identity(n_params, n_params))
        };
        
        // Convert back to Vec<Vec<f64>>
        let mut covariance_matrix = vec![vec![0.0; n_params]; n_params];
        let mut standard_errors = vec![0.0; n_params];
        
        for i in 0..n_params {
            for j in 0..n_params {
                covariance_matrix[i][j] = covariance_dmatrix[(i, j)];
            }
            standard_errors[i] = covariance_dmatrix[(i, i)].sqrt();
        }
        
        Ok((covariance_matrix, standard_errors))
    }
}

struct CompartmentSystem<'a> {
    model: &'a CompartmentModel,
    params: &'a ModelParameters,
}

impl<'a> OdeSystem for CompartmentSystem<'a> {
    fn derivatives(&self, t: f64, y: &DVector<f64>) -> DVector<f64> {
        let state = ModelState {
            compartments: y.clone(),
            time: t,
        };
        self.model.derivatives(&state, self.params)
    }

    fn dimension(&self) -> usize {
        self.model.n_compartments()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::{CompartmentModel, ModelType};
    use crate::data::*;
    use std::collections::HashMap;

    #[test]
    fn test_foce_estimator_creation() {
        let model = CompartmentModel::new(ModelType::OneCompartment).unwrap();
        let config = EstimationConfig::default()
            .with_method(super::super::EstimationMethod::Foce);
        
        let estimator = FoceEstimator::new(model, config);
        assert_eq!(estimator.model().n_compartments(), 1);
    }

    #[test]
    fn test_foce_results_creation() {
        let param_names = vec!["CL".to_string(), "V".to_string()];
        let results = FoceResults::new(2, param_names);
        
        assert_eq!(results.fixed_effects.len(), 2);
        assert_eq!(results.parameter_names.len(), 2);
        assert!(!results.converged);
    }
}