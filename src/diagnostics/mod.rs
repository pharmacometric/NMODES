use crate::data::Dataset;
use crate::saem::SaemResults;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiagnosticResults {
    pub goodness_of_fit: GoodnessOfFitMetrics,
    pub residual_analysis: ResidualAnalysis,
    pub convergence_diagnostics: ConvergenceDiagnostics,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoodnessOfFitMetrics {
    pub aic: f64,
    pub bic: f64,
    pub log_likelihood: f64,
    pub rmse: f64,
    pub mae: f64,
    pub r_squared: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidualAnalysis {
    pub residuals: Vec<f64>,
    pub standardized_residuals: Vec<f64>,
    pub weighted_residuals: Vec<f64>,
    pub residual_statistics: ResidualStatistics,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResidualStatistics {
    pub mean: f64,
    pub std_dev: f64,
    pub skewness: f64,
    pub kurtosis: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceDiagnostics {
    pub converged: bool,
    pub final_iteration: usize,
    pub parameter_stability: Vec<f64>,
    pub log_likelihood_stability: f64,
}

pub fn generate_diagnostics(
    dataset: &Dataset,
    results: &SaemResults,
) -> Result<DiagnosticResults, anyhow::Error> {
    let gof_metrics = calculate_goodness_of_fit(dataset, results)?;
    let residual_analysis = analyze_residuals(dataset, results)?;
    let convergence_diagnostics = assess_convergence(results);
    
    Ok(DiagnosticResults {
        goodness_of_fit: gof_metrics,
        residual_analysis,
        convergence_diagnostics,
    })
}

fn calculate_goodness_of_fit(
    _dataset: &Dataset,
    results: &SaemResults,
) -> Result<GoodnessOfFitMetrics, anyhow::Error> {
    // Simplified implementation
    let n_params = results.fixed_effects.len();
    let n_obs = 100; // Placeholder
    
    Ok(GoodnessOfFitMetrics {
        aic: -2.0 * results.final_log_likelihood + 2.0 * n_params as f64,
        bic: -2.0 * results.final_log_likelihood + (n_params as f64) * (n_obs as f64).ln(),
        log_likelihood: results.final_log_likelihood,
        rmse: 1.0, // Placeholder
        mae: 0.8,  // Placeholder
        r_squared: 0.95, // Placeholder
    })
}

fn analyze_residuals(
    _dataset: &Dataset,
    _results: &SaemResults,
) -> Result<ResidualAnalysis, anyhow::Error> {
    // Placeholder implementation
    let residuals = vec![0.1, -0.2, 0.05, -0.1, 0.15]; // Placeholder data
    
    let mean = residuals.iter().sum::<f64>() / residuals.len() as f64;
    let variance = residuals.iter()
        .map(|&x| (x - mean).powi(2))
        .sum::<f64>() / (residuals.len() - 1) as f64;
    let std_dev = variance.sqrt();
    
    Ok(ResidualAnalysis {
        residuals: residuals.clone(),
        standardized_residuals: residuals.iter().map(|&x| x / std_dev).collect(),
        weighted_residuals: residuals.clone(), // Simplified
        residual_statistics: ResidualStatistics {
            mean,
            std_dev,
            skewness: 0.0, // Placeholder
            kurtosis: 3.0, // Placeholder
        },
    })
}

fn assess_convergence(results: &SaemResults) -> ConvergenceDiagnostics {
    let n_recent = 100.min(results.log_likelihood_trajectory.len());
    
    let stability = if n_recent > 1 {
        let recent_ll = &results.log_likelihood_trajectory[
            (results.log_likelihood_trajectory.len() - n_recent)..
        ];
        let mean_ll = recent_ll.iter().sum::<f64>() / n_recent as f64;
        let var_ll = recent_ll.iter()
            .map(|&x| (x - mean_ll).powi(2))
            .sum::<f64>() / (n_recent - 1) as f64;
        var_ll.sqrt() / mean_ll.abs()
    } else {
        1.0
    };
    
    ConvergenceDiagnostics {
        converged: results.converged,
        final_iteration: results.n_iterations,
        parameter_stability: vec![0.01; results.fixed_effects.len()], // Placeholder
        log_likelihood_stability: stability,
    }
}