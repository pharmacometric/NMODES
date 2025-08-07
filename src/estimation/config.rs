use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum EstimationMethod {
    Saem,
    Foce,
    FoceI, // FOCE with interaction
}

impl std::fmt::Display for EstimationMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EstimationMethod::Saem => write!(f, "SAEM"),
            EstimationMethod::Foce => write!(f, "FOCE"),
            EstimationMethod::FoceI => write!(f, "FOCE-I"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EstimationConfig {
    pub method: EstimationMethod,
    pub n_iterations: usize,
    pub n_burnin: usize,
    pub n_chains: usize,
    pub mcmc_samples_per_iteration: usize,
    pub step_size: f64,
    pub target_acceptance: f64,
    pub adaptation_interval: usize,
    pub convergence_tolerance: f64,
    pub max_retries: usize,
    pub seed: Option<u64>,
    // FOCE-specific parameters
    pub foce_max_iterations: usize,
    pub foce_tolerance: f64,
    pub foce_step_size: f64,
    pub foce_interaction: bool,
}

impl Default for EstimationConfig {
    fn default() -> Self {
        Self {
            method: EstimationMethod::Saem,
            n_iterations: 1000,
            n_burnin: 200,
            n_chains: 4,
            mcmc_samples_per_iteration: 10,
            step_size: 0.1,
            target_acceptance: 0.44,
            adaptation_interval: 50,
            convergence_tolerance: 0.001,
            max_retries: 3,
            seed: Some(12345), // Default seed for reproducibility
            foce_max_iterations: 100,
            foce_tolerance: 1e-6,
            foce_step_size: 1e-4,
            foce_interaction: false,
        }
    }
}

impl EstimationConfig {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_method(mut self, method: EstimationMethod) -> Self {
        self.method = method;
        self
    }

    pub fn with_iterations(mut self, n_iterations: usize) -> Self {
        self.n_iterations = n_iterations;
        self
    }

    pub fn with_burnin(mut self, n_burnin: usize) -> Self {
        self.n_burnin = n_burnin;
        self
    }

    pub fn with_chains(mut self, n_chains: usize) -> Self {
        self.n_chains = n_chains;
        self
    }

    pub fn with_step_size(mut self, step_size: f64) -> Self {
        self.step_size = step_size;
        self
    }

    pub fn with_seed(mut self, seed: Option<u64>) -> Self {
        self.seed = seed;
        self
    }

    pub fn with_foce_iterations(mut self, foce_max_iterations: usize) -> Self {
        self.foce_max_iterations = foce_max_iterations;
        self
    }

    pub fn with_foce_tolerance(mut self, foce_tolerance: f64) -> Self {
        self.foce_tolerance = foce_tolerance;
        self
    }

    pub fn with_foce_interaction(mut self, foce_interaction: bool) -> Self {
        self.foce_interaction = foce_interaction;
        self
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.n_iterations == 0 {
            return Err("Number of iterations must be positive".to_string());
        }
        
        if self.n_burnin >= self.n_iterations {
            return Err("Burn-in period must be less than total iterations".to_string());
        }
        
        if self.n_chains == 0 {
            return Err("Number of chains must be positive".to_string());
        }
        
        if self.step_size <= 0.0 {
            return Err("Step size must be positive".to_string());
        }
        
        if !(0.0..=1.0).contains(&self.target_acceptance) {
            return Err("Target acceptance rate must be between 0 and 1".to_string());
        }
        
        if self.foce_max_iterations == 0 {
            return Err("FOCE max iterations must be positive".to_string());
        }
        
        if self.foce_tolerance <= 0.0 {
            return Err("FOCE tolerance must be positive".to_string());
        }
        
        if self.foce_step_size <= 0.0 {
            return Err("FOCE step size must be positive".to_string());
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = EstimationConfig::default();
        assert!(config.validate().is_ok());
        assert_eq!(config.n_iterations, 1000);
        assert_eq!(config.n_burnin, 200);
    }

    #[test]
    fn test_config_validation() {
        let mut config = EstimationConfig::default();
        config.n_iterations = 0;
        assert!(config.validate().is_err());
        
        config.n_iterations = 100;
        config.n_burnin = 150;
        assert!(config.validate().is_err());
    }
}