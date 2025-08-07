use crate::data::Individual;
use crate::models::{CompartmentModel, ModelParameters};
use crate::solver::{OdeSolver, SolverConfig};
use nalgebra::{DVector, DMatrix, Dynamic};
use rand::prelude::*;
use rand_distr::StandardNormal;
use rand::{SeedableRng, rngs::StdRng};

#[derive(Debug, Clone)]
pub struct McmcConfig {
    pub n_samples: usize,
    pub step_size: f64,
    pub target_acceptance: f64,
    pub seed: Option<u64>,
}

impl Default for McmcConfig {
    fn default() -> Self {
        Self {
            n_samples: 100,
            step_size: 0.1,
            target_acceptance: 0.44,
            seed: None,
        }
    }
}

pub struct McmcSampler<'a> {
    model: &'a CompartmentModel,
    solver: &'a dyn OdeSolver,
    config: McmcConfig,
    rng: StdRng,
}

impl<'a> McmcSampler<'a> {
    pub fn new(
        model: &'a CompartmentModel,
        solver: &'a dyn OdeSolver,
        config: McmcConfig,
    ) -> Self {
        let rng = if let Some(seed) = config.seed {
            StdRng::seed_from_u64(seed)
        } else {
            StdRng::from_entropy()
        };
        
        Self {
            model,
            solver,
            config,
            rng,
        }
    }

    pub fn sample_individual_parameters(
        &mut self,
        individual: &Individual,
        population_params: &ModelParameters,
        initial_params: &Vec<f64>,
    ) -> Result<(Vec<f64>, f64), anyhow::Error> {
        let mut current_params = initial_params.clone();
        let mut current_log_likelihood = self.log_likelihood(individual, &current_params, population_params)?;
        
        let mut n_accepted = 0;
        let n_params = current_params.len();
        
        for _ in 0..self.config.n_samples {
            // Propose new parameters
            let mut proposed_params = current_params.clone();
            
            for i in 0..n_params {
                let step: f64 = self.rng.sample(StandardNormal);
                proposed_params[i] += self.config.step_size * step;
                
                // Apply bounds: ensure exp(param) > 0 by keeping param > -10
                proposed_params[i] = proposed_params[i].max(-10.0);
            }
            
            // Calculate log-likelihood for proposed parameters
            let proposed_log_likelihood = self.log_likelihood(individual, &proposed_params, population_params)?;
            
            // Metropolis-Hastings acceptance
            let log_alpha = proposed_log_likelihood - current_log_likelihood;
            let alpha = log_alpha.exp().min(1.0);
            
            if self.rng.gen::<f64>() < alpha {
                current_params = proposed_params;
                current_log_likelihood = proposed_log_likelihood;
                n_accepted += 1;
            }
        }
        
        let _acceptance_rate = n_accepted as f64 / self.config.n_samples as f64;
        
        Ok((current_params, current_log_likelihood))
    }

    fn log_likelihood(
        &self,
        individual: &Individual,
        individual_params: &Vec<f64>,
        population_params: &ModelParameters,
    ) -> Result<f64, anyhow::Error> {
        // Log-likelihood = log p(y|θ) + log p(θ|μ,Ω)
        // where y are observations, θ are individual parameters, μ are population means, Ω is covariance
        
        let data_log_likelihood = self.data_log_likelihood(individual, individual_params)?;
        let prior_log_likelihood = self.prior_log_likelihood(individual_params, population_params);
        
        Ok(data_log_likelihood + prior_log_likelihood)
    }

    fn data_log_likelihood(
        &self,
        individual: &Individual,
        individual_params: &Vec<f64>,
    ) -> Result<f64, anyhow::Error> {
        let predictions = self.predict_concentrations(individual, individual_params)?;
        let mut log_likelihood = 0.0;
        
        for (obs, pred) in individual.observations().iter().zip(predictions.iter()) {
            if obs.value > 0.0 && *pred > 0.0 {
                // Log-normal error model
                let log_obs = obs.value.ln();
                let log_pred = pred.ln();
                let residual = log_obs - log_pred;
                
                // Assume proportional error model
                let sigma = 0.1; // This should come from population_params.residual_variance
                log_likelihood -= 0.5 * (residual / sigma).powi(2);
                log_likelihood -= 0.5 * (2.0 * std::f64::consts::PI * sigma.powi(2)).ln();
            }
        }
        
        Ok(log_likelihood)
    }

    fn prior_log_likelihood(
        &self,
        individual_params: &Vec<f64>,
        population_params: &ModelParameters,
    ) -> f64 {
        // Multivariate normal prior: θ ~ N(μ, Ω)
        let mut diff = vec![0.0; individual_params.len()];
        for i in 0..individual_params.len() {
            diff[i] = individual_params[i] - population_params.fixed_effects[i];
        }
        
        // Simplified calculation assuming diagonal covariance matrix
        let mut quadratic_form = 0.0;
        let mut det_omega = 1.0;
        for i in 0..diff.len() {
            let variance = population_params.random_effects_variance[i][i];
            quadratic_form += diff[i] * diff[i] / variance;
            det_omega *= variance;
        }
        
        -0.5 * quadratic_form - 0.5 * det_omega.ln() - 
        0.5 * (individual_params.len() as f64) * (2.0 * std::f64::consts::PI).ln()
    }

    fn predict_concentrations(
        &self,
        individual: &Individual,
        individual_params: &Vec<f64>,
    ) -> Result<Vec<f64>, anyhow::Error> {
        // Create temporary parameters for this individual
        let mut temp_params = self.model.default_parameters();
        temp_params.fixed_effects = individual_params.clone();
        
        let mut predictions = Vec::new();
        let solver_config = SolverConfig::default();
        
        // Simulate the PK profile
        let system = CompartmentSystem {
            model: &self.model,
            params: &temp_params,
        };
        
        let mut current_state = crate::models::ModelState::new(self.model.n_compartments());
        let mut last_time = 0.0;
        
        // Apply dosing events
        for dose in individual.dosing_records() {
            if dose.time > last_time {
                // Integrate from last_time to dose.time
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
            
            // Apply dose
            current_state.add_dose(dose.compartment as usize, dose.amount);
            last_time = dose.time;
        }
        
        // Predict concentrations at observation times
        for obs in individual.observations() {
            if obs.time > last_time {
                let (_, solutions) = self.solver.solve(
                    &system,
                    (last_time, obs.time),
                    &current_state.compartments,
                    &solver_config,
                )?;
                current_state.compartments = solutions.into_iter().last().unwrap_or(current_state.compartments.clone());
                current_state.time = obs.time;
                last_time = obs.time;
            }
            
            let concentration = self.model.observation_function(
                &current_state,
                &temp_params,
                obs.compartment as usize,
            );
            predictions.push(concentration);
        }
        
        Ok(predictions)
    }
}

struct CompartmentSystem<'a> {
    model: &'a CompartmentModel,
    params: &'a ModelParameters,
}

impl<'a> crate::solver::OdeSystem for CompartmentSystem<'a> {
    fn derivatives(&self, t: f64, y: &DVector<f64>) -> DVector<f64> {
        let state = crate::models::ModelState {
            compartments: y.clone(),
            time: t,
        };
        self.model.derivatives(&state, self.params)
    }

    fn dimension(&self) -> usize {
        self.model.n_compartments()
    }
}