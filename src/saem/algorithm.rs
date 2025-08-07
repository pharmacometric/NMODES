use super::{SaemResults, McmcSampler, McmcConfig};
use super::{ParameterStatistics, OmegaStatistics};
use crate::data::Dataset;
use crate::models::{CompartmentModel, ModelParameters, ModelState};
use crate::estimation::EstimationConfig;
use crate::solver::{OdeSolver, OdeSystem, RungeKuttaSolver, SolverConfig};
use anyhow::{Context, Result};
use log::{info, debug, warn};
use std::collections::HashMap;
use nalgebra::DVector;

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

pub struct SaemEstimator {
    model: CompartmentModel,
    config: EstimationConfig,
    solver: Box<dyn OdeSolver + Send + Sync>,
}

impl SaemEstimator {
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

    // CORRECTED: Removed duplicate function definition
    pub fn fit(&mut self, dataset: &Dataset) -> Result<SaemResults> {
        info!("Starting SAEM estimation for {} individuals", dataset.n_individuals());
        
        let n_params = self.model.parameter_names().len();
        let parameter_names = self.model.parameter_names();
        let mut results = SaemResults::new(n_params, parameter_names.clone());
        
        let mut current_params = self.model.default_parameters();
        results.set_fixed_effects(&current_params.get_fixed_effects_vector());
        results.set_random_effects_variance(&current_params.get_random_effects_matrix());
        results.residual_variance = current_params.residual_variance;

        let mut individual_params: HashMap<i32, Vec<f64>> = HashMap::new();
        for (&id, _) in dataset.individuals() {
            individual_params.insert(id, current_params.fixed_effects.clone());
        }

        let mut sa_sum_theta = vec![0.0; n_params];
        let mut sa_sum_theta_sq = vec![vec![0.0; n_params]; n_params];
        let mut sa_sum_sigma = 0.0;

        for iteration in 0..self.config.n_iterations {
            debug!("SAEM iteration {}/{}", iteration + 1, self.config.n_iterations);
            
            let mut iteration_log_likelihood = 0.0;

            let gamma = if iteration < self.config.n_burnin {
                1.0
            } else {
                1.0 / ((iteration - self.config.n_burnin + 1) as f64).powf(0.7)
            };

            for (&id, individual) in dataset.individuals() {
                let mcmc_config = McmcConfig {
                    n_samples: self.config.mcmc_samples_per_iteration,
                    step_size: self.config.step_size,
                    target_acceptance: self.config.target_acceptance,
                    seed: self.config.seed.map(|s| s.wrapping_add(iteration as u64).wrapping_add(id as u64)),
                };

                let mut sampler = McmcSampler::new(
                    &self.model,
                    self.solver.as_ref(),
                    mcmc_config,
                );

                let (new_params, log_like) = sampler.sample_individual_parameters(
                    individual,
                    &current_params,
                    individual_params.get(&id).unwrap(),
                ).with_context(|| format!("MCMC sampling failed for individual {}", id))?;

                individual_params.insert(id, new_params);
                iteration_log_likelihood += log_like;
            }

            self.update_population_parameters(
                &individual_params,
                &mut current_params,
                &mut sa_sum_theta,
                &mut sa_sum_theta_sq,
                &mut sa_sum_sigma,
                gamma,
                dataset,
            );

            results.parameter_trajectory.push(current_params.fixed_effects.clone());
            results.log_likelihood_trajectory.push(iteration_log_likelihood);

            if iteration > self.config.n_burnin && iteration % 50 == 0 {
                if self.check_convergence(&results) {
                    info!("Convergence achieved at iteration {}", iteration);
                    results.converged = true;
                    break;
                }
            }

            if iteration % 100 == 0 {
                info!("Iteration {}: Log-likelihood = {:.3}", 
                      iteration, iteration_log_likelihood);
            }
        }

        results.set_fixed_effects(&current_params.get_fixed_effects_vector());
        results.set_random_effects_variance(&current_params.get_random_effects_matrix());
        results.residual_variance = current_params.residual_variance;
        results.final_log_likelihood = results.log_likelihood_trajectory.last().copied()
            .unwrap_or(f64::NEG_INFINITY);
        results.objective_function_value = -2.0 * results.final_log_likelihood;
        results.n_iterations = results.parameter_trajectory.len();
        results.individual_parameters = individual_params;

        // Calculate parameter statistics
        self.calculate_parameter_statistics(&mut results);
        self.calculate_omega_statistics(&mut results, dataset);

        info!("SAEM estimation completed. Final log-likelihood: {:.3}, Objective function: {:.3}", 
              results.final_log_likelihood, results.objective_function_value);

        Ok(results)
    }

    fn update_population_parameters(
        &self,
        individual_params: &HashMap<i32, Vec<f64>>,
        current_params: &mut ModelParameters,
        sa_sum_theta: &mut Vec<f64>,
        sa_sum_theta_sq: &mut Vec<Vec<f64>>,
        sa_sum_sigma: &mut f64,
        gamma: f64,
        dataset: &Dataset,
    ) {
        let n_individuals = individual_params.len() as f64;
        
        let mut mean_individual_params = vec![0.0; current_params.n_parameters()];
        for params in individual_params.values() {
            for (i, param) in params.iter().enumerate() {
                mean_individual_params[i] += param;
            }
        }
        for val in mean_individual_params.iter_mut() {
            *val /= n_individuals;
        }
        
        for i in 0..sa_sum_theta.len() {
            sa_sum_theta[i] = (1.0 - gamma) * sa_sum_theta[i] + gamma * mean_individual_params[i];
            // Apply bounds to prevent parameters from becoming too negative
            sa_sum_theta[i] = sa_sum_theta[i].max(-10.0);
        }
        current_params.fixed_effects = sa_sum_theta.clone();
        
        let mut sum_outer_products = vec![vec![0.0; current_params.n_parameters()]; current_params.n_parameters()];
        for params in individual_params.values() {
            for i in 0..params.len() {
                for j in 0..params.len() {
                    let centered_i = params[i] - current_params.fixed_effects[i];
                    let centered_j = params[j] - current_params.fixed_effects[j];
                    sum_outer_products[i][j] += centered_i * centered_j;
                }
            }
        }
        
        for i in 0..sum_outer_products.len() {
            for j in 0..sum_outer_products[i].len() {
                let mean_outer_product = sum_outer_products[i][j] / n_individuals;
                sa_sum_theta_sq[i][j] = (1.0 - gamma) * sa_sum_theta_sq[i][j] + gamma * mean_outer_product;
            }
        }
        current_params.random_effects_variance = sa_sum_theta_sq.clone();
        
        let mut residual_sum = 0.0;
        let mut total_observations = 0;
        
        for (&id, individual) in dataset.individuals() {
            if let Some(ind_params) = individual_params.get(&id) {
                let mut temp_params = current_params.clone();
                temp_params.fixed_effects = ind_params.clone();
                
                // CORRECTED: Handle potential errors from prediction
                let predicted = match self.predict_individual(individual, &temp_params) {
                    Ok(p) => p,
                    Err(e) => {
                        warn!("Could not predict for individual {}: {}. Skipping for residual variance update.", id, e);
                        continue;
                    }
                };
                
                for (obs, pred) in individual.observations().iter().zip(predicted.iter()) {
                    let residual = (obs.value - pred).powi(2);
                    residual_sum += residual;
                    total_observations += 1;
                }
            }
        }
        
        // CORRECTED: Add check to prevent division by zero
        if total_observations > 0 {
            let empirical_residual_var = residual_sum / total_observations as f64;
            *sa_sum_sigma = (1.0 - gamma) * (*sa_sum_sigma) + gamma * empirical_residual_var;
            current_params.residual_variance = *sa_sum_sigma;
        }
    }

    fn predict_individual(
        &self,
        individual: &crate::data::Individual,
        params: &ModelParameters,
    ) -> Result<Vec<f64>, anyhow::Error> {
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

    fn check_convergence(&self, results: &SaemResults) -> bool {
        let window_size = 50;
        if results.log_likelihood_trajectory.len() < window_size {
            return false;
        }

        let recent_values = &results.log_likelihood_trajectory[
            (results.log_likelihood_trajectory.len() - window_size)..
        ];
        
        let mean = recent_values.iter().sum::<f64>() / window_size as f64;
        
        // Avoid division by zero or a very small number
        if mean.abs() < 1e-9 {
            return false;
        }
        
        let variance = recent_values.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>() / (window_size - 1) as f64;
        
        let coefficient_of_variation = variance.sqrt() / mean.abs();
        
        coefficient_of_variation < self.config.convergence_tolerance
    }

    fn calculate_parameter_statistics(&self, results: &mut SaemResults) {
        let n_recent = 100.min(results.parameter_trajectory.len());
        if n_recent < 10 {
            return; // Not enough data for reliable statistics
        }

        let recent_params: Vec<&Vec<f64>> = results.parameter_trajectory
            .iter()
            .rev()
            .take(n_recent)
            .collect();

        for (param_idx, param_name) in results.parameter_names.iter().enumerate() {
            let param_values: Vec<f64> = recent_params
                .iter()
                .map(|params| params[param_idx])
                .collect();

            let mean = param_values.iter().sum::<f64>() / param_values.len() as f64;
            let variance = param_values.iter()
                .map(|&x| (x - mean).powi(2))
                .sum::<f64>() / (param_values.len() - 1) as f64;
            let std_error = variance.sqrt() / (param_values.len() as f64).sqrt();
            
            // Calculate %RSE (Relative Standard Error)
            let rse_percent = if mean.abs() > 1e-10 {
                (std_error / mean.abs()) * 100.0
            } else {
                0.0
            };

            results.parameter_statistics.push(ParameterStatistics {
                name: param_name.clone(),
                estimate: results.fixed_effects[param_idx],
                rse_percent,
            });
        }
    }

    fn calculate_omega_statistics(&self, results: &mut SaemResults, dataset: &Dataset) {
        let n_params = results.parameter_names.len();
        
        // Calculate empirical Bayes estimates (EBEs) for shrinkage calculation
        let mut empirical_variances = vec![0.0; n_params];
        let n_individuals = dataset.n_individuals() as f64;
        
        if n_individuals > 1.0 {
            for param_idx in 0..n_params {
                let individual_values: Vec<f64> = results.individual_parameters
                    .values()
                    .map(|params| params[param_idx])
                    .collect();
                
                let mean = individual_values.iter().sum::<f64>() / n_individuals;
                let empirical_var = individual_values.iter()
                    .map(|&x| (x - mean).powi(2))
                    .sum::<f64>() / (n_individuals - 1.0);
                
                empirical_variances[param_idx] = empirical_var;
            }
        }

        // Generate omega statistics
        for i in 0..n_params {
            for j in 0..n_params {
                let omega_estimate = results.random_effects_variance[i][j];
                
                let shrinkage_percent = if i == j && empirical_variances[i] > 1e-10 {
                    // Shrinkage = (1 - empirical_variance / omega) * 100%
                    let shrinkage = (1.0 - empirical_variances[i] / omega_estimate.abs()) * 100.0;
                    Some(shrinkage.max(0.0).min(100.0)) // Clamp between 0-100%
                } else {
                    None
                };

                results.omega_statistics.push(OmegaStatistics {
                    parameter_i: results.parameter_names[i].clone(),
                    parameter_j: results.parameter_names[j].clone(),
                    estimate: omega_estimate,
                    shrinkage_percent,
                });
            }
        }
    }
}