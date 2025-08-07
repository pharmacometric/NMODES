use crate::saem::SaemResults;
use crate::diagnostics::DiagnosticResults;
use crate::data::Dataset;
use crate::models::{CompartmentModel, ModelParameters, ModelState};
use crate::solver::{RungeKuttaSolver, OdeSolver, SolverConfig, OdeSystem};
use nalgebra::DVector;
use std::path::Path;
use std::fs;
use log::info;

pub fn save_results(
    output_dir: &Path,
    results: &SaemResults,
    diagnostics: &DiagnosticResults,
    dataset: &Dataset,
    model: &CompartmentModel,
) -> Result<(), anyhow::Error> {
    info!("Saving results to {:?}", output_dir);
    
    // Ensure output directory exists
    fs::create_dir_all(output_dir)?;
    
    // Save parameter estimates
    save_parameter_estimates(output_dir, results)?;
    
    // Save diagnostics
    save_diagnostics(output_dir, diagnostics)?;
    
    // Save parameter trajectory
    save_parameter_trajectory(output_dir, results)?;
    
    // Save summary report
    save_summary_report(output_dir, results, diagnostics)?;
    
    // Save predictions CSV
    save_predictions_csv(output_dir, results, dataset, model)?;
    
    info!("All results saved successfully");
    Ok(())
}

fn save_parameter_estimates(
    output_dir: &Path,
    results: &SaemResults,
) -> Result<(), anyhow::Error> {
    let estimates_file = output_dir.join("parameter_estimates.json");
    let json_content = serde_json::to_string_pretty(results)?;
    fs::write(estimates_file, json_content)?;
    Ok(())
}

fn save_diagnostics(
    output_dir: &Path,
    diagnostics: &DiagnosticResults,
) -> Result<(), anyhow::Error> {
    let diagnostics_file = output_dir.join("diagnostics.json");
    let json_content = serde_json::to_string_pretty(diagnostics)?;
    fs::write(diagnostics_file, json_content)?;
    Ok(())
}

fn save_parameter_trajectory(
    output_dir: &Path,
    results: &SaemResults,
) -> Result<(), anyhow::Error> {
    let trajectory_file = output_dir.join("parameter_trajectory.csv");
    let mut wtr = csv::Writer::from_path(trajectory_file)?;
    
    // Write header
    let mut header = vec!["iteration".to_string()];
    for i in 0..results.fixed_effects.len() {
        header.push(format!("param_{}", i));
    }
    header.push("log_likelihood".to_string());
    wtr.write_record(&header)?;
    
    // Write data
    for (iter, (params, &log_like)) in results.parameter_trajectory.iter()
        .zip(results.log_likelihood_trajectory.iter())
        .enumerate()
    {
        let mut record = vec![iter.to_string()];
        for &param in params.iter() {
            record.push(param.to_string());
        }
        record.push(log_like.to_string());
        wtr.write_record(&record)?;
    }
    
    wtr.flush()?;
    Ok(())
}

fn save_summary_report(
    output_dir: &Path,
    results: &SaemResults,
    diagnostics: &DiagnosticResults,
) -> Result<(), anyhow::Error> {
    let report_file = output_dir.join("summary_report.txt");
    
    let mut report = String::new();
    report.push_str("NMODES SAEM Analysis Summary Report\n");
    report.push_str("=================================\n\n");
    
    report.push_str(&format!("Model Convergence: {}\n", results.converged));
    report.push_str(&format!("Total Iterations: {}\n", results.n_iterations));
    report.push_str(&format!("Final Log-Likelihood: {:.6}\n", results.final_log_likelihood));
    report.push_str(&format!("Objective Function Value: {:.6}\n", results.objective_function_value));
    report.push_str(&format!("Number of Individuals: {}\n", results.individual_parameters.len()));
    report.push_str(&format!("Number of Observations: {}\n", 
        results.individual_parameters.values().map(|_| 1).sum::<usize>())); // Simplified
    report.push_str(&format!("AIC: {:.6}\n", diagnostics.goodness_of_fit.aic));
    report.push_str(&format!("BIC: {:.6}\n", diagnostics.goodness_of_fit.bic));
    report.push_str(&format!("R-squared: {:.6}\n", diagnostics.goodness_of_fit.r_squared));
    report.push_str(&format!("RMSE: {:.6}\n", diagnostics.goodness_of_fit.rmse));
    
    report.push_str("\nFixed Effects Parameter Estimates:\n");
    report.push_str("----------------------------------\n");
    report.push_str(&format!("{:<10} {:<12} {:<10}\n", "Parameter", "Estimate", "%RSE"));
    report.push_str(&format!("{:<10} {:<12} {:<10}\n", "---------", "--------", "----"));
    for param_stat in &results.parameter_statistics {
        report.push_str(&format!("{:<10} {:<12.6} {:<10.2}\n", 
            param_stat.name, param_stat.estimate, param_stat.rse_percent));
    }
    
    report.push_str(&format!("\nResidual Error Variance: {:.6}\n", results.residual_variance));
    
    report.push_str("\nRandom Effects Variance (Omega):\n");
    report.push_str("-------------------------------\n");
    report.push_str(&format!("{:<15} {:<12} {:<12}\n", "Parameter", "Estimate", "Shrinkage%"));
    report.push_str(&format!("{:<15} {:<12} {:<12}\n", "---------", "--------", "----------"));
    for omega_stat in &results.omega_statistics {
        if omega_stat.parameter_i == omega_stat.parameter_j {
            let shrinkage_text = if let Some(shrinkage) = omega_stat.shrinkage_percent {
                format!("{:.1}", shrinkage)
            } else {
                "N/A".to_string()
            };
            report.push_str(&format!("{:<15} {:<12.6} {:<12}\n", 
                format!("{}({})", omega_stat.parameter_i, omega_stat.parameter_i),
                omega_stat.estimate, shrinkage_text));
        } else if omega_stat.estimate.abs() > 1e-10 {
            report.push_str(&format!("{:<15} {:<12.6} {:<12}\n", 
                format!("{}({})", omega_stat.parameter_i, omega_stat.parameter_j),
                omega_stat.estimate, "N/A"));
        }
    }
    
    fs::write(report_file, report)?;
    Ok(())
}

fn save_predictions_csv(
    output_dir: &Path,
    results: &SaemResults,
    dataset: &Dataset,
    model: &CompartmentModel,
) -> Result<(), anyhow::Error> {
    let predictions_file = output_dir.join("predictions.csv");
    let mut wtr = csv::Writer::from_path(predictions_file)?;
    
    // Write header
    wtr.write_record(&["ID", "TIME", "DV", "IPRED", "PRED"])?;
    
    let solver = RungeKuttaSolver::new();
    let solver_config = SolverConfig::default();
    
    // Calculate population predictions using population parameters
    let pop_params = model.default_parameters();
    let mut pop_params_final = pop_params.clone();
    pop_params_final.fixed_effects = results.fixed_effects.clone();
    
    for (&id, individual) in dataset.individuals() {
        // Get individual parameters
        let ind_params = results.individual_parameters.get(&id)
            .unwrap_or(&results.fixed_effects);
        
        // Calculate individual predictions (IPRED)
        let ipred = calculate_predictions(individual, ind_params, model, &solver, &solver_config)?;
        
        // Calculate population predictions (PRED) 
        let pred = calculate_predictions(individual, &results.fixed_effects, model, &solver, &solver_config)?;
        
        // Write data for each observation
        for (obs_idx, obs) in individual.observations().iter().enumerate() {
            let ipred_value = ipred.get(obs_idx).copied().unwrap_or(0.0);
            let pred_value = pred.get(obs_idx).copied().unwrap_or(0.0);
            
            wtr.write_record(&[
                id.to_string(),
                obs.time.to_string(),
                obs.value.to_string(),
                ipred_value.to_string(),
                pred_value.to_string(),
            ])?;
        }
    }
    
    wtr.flush()?;
    Ok(())
}

fn calculate_predictions(
    individual: &crate::data::Individual,
    params: &[f64],
    model: &CompartmentModel,
    solver: &dyn OdeSolver,
    solver_config: &SolverConfig,
) -> Result<Vec<f64>, anyhow::Error> {
    use crate::models::{ModelState, ModelParameters};
    use crate::solver::OdeSystem;
    
    // Create temporary parameters for this prediction
    let mut temp_params = model.default_parameters();
    temp_params.fixed_effects = params.to_vec();
    
    let system = CompartmentSystem {
        model,
        params: &temp_params,
    };
    
    let mut predictions = Vec::new();
    let mut current_state = ModelState::new(model.n_compartments());
    let mut last_time = 0.0;
    
    // Apply dosing events
    for dose in individual.dosing_records() {
        if dose.time > last_time {
            let final_state = solver.solve_to_time(
                &system,
                last_time,
                dose.time,
                &current_state.compartments,
                solver_config,
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
            let final_state = solver.solve_to_time(
                &system,
                last_time,
                obs.time,
                &current_state.compartments,
                solver_config,
            )?;
            current_state.compartments = final_state;
            current_state.time = obs.time;
            last_time = obs.time;
        }
        
        let concentration = model.observation_function(
            &current_state,
            &temp_params,
            obs.compartment as usize,
        );
        predictions.push(concentration);
    }
    
    Ok(predictions)
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