use clap::{Arg, Command};
use log::{info, warn, error};
use std::path::{Path, PathBuf};
use std::fs;
use nmodes::{Dataset, CompartmentModel, ModelType, SaemEstimator, RungeKuttaSolver, SolverConfig};
use nmodes::{EstimationConfig, EstimationMethod, FoceEstimator, estimation, FoceResults, SaemResults};
use nmodes::{diagnostics, output, validation};
use anyhow::{Result, anyhow};

#[derive(Debug)]
struct CliArgs {
    dataset_path: PathBuf,
    model_types: Vec<ModelType>,
    estimation_methods: Vec<EstimationMethod>,
    output_dir: PathBuf,
    iterations: usize,
    burn_in: usize,
    chains: usize,
    compare_results: bool,
}

fn main() -> Result<()> {
    env_logger::init();
    
    let matches = Command::new("NMODES - Nonlinear Mixed Effects Differential Equation Solver")
        .version("1.0.0")
        .author("NMODES Team")
        .about("Population pharmacokinetics modeling using SAEM and FOCE methods")
        .arg(
            Arg::new("dataset")
                .short('d')
                .long("dataset")
                .value_name("FILE")
                .help("Path to NONMEM-style dataset CSV file")
                .required(true)
        )
        .arg(
            Arg::new("model")
                .short('m')
                .long("model")
                .value_name("TYPE")
                .help("Compartment model type(s): 1comp, 2comp, 3comp, or 'all' for all models")
                .default_value("1comp")
                .action(clap::ArgAction::Append)
        )
        .arg(
            Arg::new("method")
                .short('e')
                .long("method")
                .value_name("METHOD")
                .help("Estimation method(s): saem, foce, foce-i, or 'all' for all methods")
                .default_value("saem")
                .action(clap::ArgAction::Append)
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("DIR")
                .help("Output directory for results")
                .default_value("./output")
        )
        .arg(
            Arg::new("iterations")
                .short('i')
                .long("iterations")
                .value_name("N")
                .help("Number of SAEM iterations")
                .default_value("1000")
        )
        .arg(
            Arg::new("burn-in")
                .short('b')
                .long("burn-in")
                .value_name("N")
                .help("Number of burn-in iterations")
                .default_value("200")
        )
        .arg(
            Arg::new("chains")
                .short('c')
                .long("chains")
                .value_name("N")
                .help("Number of MCMC chains")
                .default_value("4")
        )
        .arg(
            Arg::new("compare")
                .long("compare")
                .help("Generate comparison report across models and methods")
                .action(clap::ArgAction::SetTrue)
        )
        .get_matches();

    let args = CliArgs {
        dataset_path: PathBuf::from(matches.get_one::<String>("dataset").unwrap()),
        model_types: parse_model_types(matches.get_many::<String>("model").unwrap().collect())?,
        estimation_methods: parse_estimation_methods(matches.get_many::<String>("method").unwrap().collect())?,
        output_dir: PathBuf::from(matches.get_one::<String>("output").unwrap()),
        iterations: matches.get_one::<String>("iterations").unwrap().parse()?,
        burn_in: matches.get_one::<String>("burn-in").unwrap().parse()?,
        chains: matches.get_one::<String>("chains").unwrap().parse()?,
        compare_results: matches.get_flag("compare"),
    };

    run_analysis(args)
}

fn parse_model_types(model_strs: Vec<&String>) -> Result<Vec<ModelType>> {
    let mut model_types = Vec::new();
    
    for model_str in model_strs {
        if model_str == "all" {
            return Ok(vec![
                ModelType::OneCompartment,
                ModelType::TwoCompartment,
                ModelType::ThreeCompartment,
            ]);
        }
        
        let model_type = match model_str.as_str() {
            "1comp" => ModelType::OneCompartment,
            "2comp" => ModelType::TwoCompartment,
            "3comp" => ModelType::ThreeCompartment,
            _ => return Err(anyhow!("Invalid model type: {}", model_str)),
        };
        
        if !model_types.contains(&model_type) {
            model_types.push(model_type);
        }
    }
    
    if model_types.is_empty() {
        model_types.push(ModelType::OneCompartment); // Default
    }
    
    Ok(model_types)
}

fn parse_model_type(model_str: &str) -> Result<ModelType> {
    match model_str {
        "1comp" => Ok(ModelType::OneCompartment),
        "2comp" => Ok(ModelType::TwoCompartment),
        "3comp" => Ok(ModelType::ThreeCompartment),
        _ => Err(anyhow!("Invalid model type: {}", model_str)),
    }
}

fn parse_estimation_methods(method_strs: Vec<&String>) -> Result<Vec<EstimationMethod>> {
    let mut estimation_methods = Vec::new();
    
    for method_str in method_strs {
        if method_str == "all" {
            return Ok(vec![
                EstimationMethod::Saem,
                EstimationMethod::Foce,
                EstimationMethod::FoceI,
            ]);
        }
        
        let estimation_method = match method_str.as_str() {
            "saem" => EstimationMethod::Saem,
            "foce" => EstimationMethod::Foce,
            "foce-i" => EstimationMethod::FoceI,
            _ => return Err(anyhow!("Invalid estimation method: {}", method_str)),
        };
        
        if !estimation_methods.contains(&estimation_method) {
            estimation_methods.push(estimation_method);
        }
    }
    
    if estimation_methods.is_empty() {
        estimation_methods.push(EstimationMethod::Saem); // Default
    }
    
    Ok(estimation_methods)
}

fn parse_estimation_method(method_str: &str) -> Result<EstimationMethod> {
    match method_str {
        "saem" => Ok(EstimationMethod::Saem),
        "foce" => Ok(EstimationMethod::Foce),
        "foce-i" => Ok(EstimationMethod::FoceI),
        _ => Err(anyhow!("Invalid estimation method: {}", method_str)),
    }
}

fn run_analysis(args: CliArgs) -> Result<()> {
    info!("Starting NMODES analysis");
    info!("Dataset: {:?}", args.dataset_path);
    info!("Model types: {:?}", args.model_types);
    info!("Estimation methods: {:?}", args.estimation_methods);
    info!("Output directory: {:?}", args.output_dir);

    // Create output directory
    std::fs::create_dir_all(&args.output_dir)?;

    // Load and validate dataset
    info!("Loading dataset...");
    let dataset = Dataset::from_csv(&args.dataset_path)?;
    info!("Loaded {} individuals with {} observations", 
          dataset.n_individuals(), dataset.n_observations());

    // Validate dataset
    validation::validate_dataset(&dataset)?;

    // Store all results for comparison
    let mut all_results: Vec<AnalysisResult> = Vec::new();
    
    // Run analysis for each model and method combination
    for model_type in &args.model_types {
        for estimation_method in &args.estimation_methods {
            info!("Running {} estimation with {} model", estimation_method, model_type);
            
            // Create model
            let model = CompartmentModel::new(model_type.clone())?;
            
            // Configure estimation
            let config = EstimationConfig {
                method: estimation_method.clone(),
                n_iterations: args.iterations,
                n_burnin: args.burn_in,
                n_chains: args.chains,
                step_size: 0.1,
                target_acceptance: 0.44,
                adaptation_interval: 50,
                foce_max_iterations: if matches!(estimation_method, EstimationMethod::Foce | EstimationMethod::FoceI) {
                    args.iterations
                } else {
                    100
                },
                foce_tolerance: 1e-6,
                foce_step_size: 1e-4,
                foce_interaction: matches!(estimation_method, EstimationMethod::FoceI),
                ..Default::default()
            };
            
            // Create method-specific output directory
            let method_output_dir = args.output_dir.join(format!("{}_{}", model_type, estimation_method));
            std::fs::create_dir_all(&method_output_dir)?;
            
            // Run estimation
            let analysis_result = match estimation_method {
                EstimationMethod::Saem => {
                    let mut estimator = SaemEstimator::new(model, config);
                    let results = estimator.fit(&dataset)?;
                    
                    // Generate diagnostics
                    let diagnostics = diagnostics::generate_diagnostics(&dataset, &results)?;
                    
                    // Save SAEM results
                    output::save_results(&method_output_dir, &results, &diagnostics, &dataset, estimator.model())?;
                    
                    AnalysisResult {
                        model_type: model_type.clone(),
                        estimation_method: estimation_method.clone(),
                        objective_function_value: results.objective_function_value,
                        final_log_likelihood: results.final_log_likelihood,
                        converged: results.converged,
                        n_iterations: results.n_iterations,
                        fixed_effects: results.fixed_effects.clone(),
                        parameter_names: results.parameter_names.clone(),
                        aic: diagnostics.goodness_of_fit.aic,
                        bic: diagnostics.goodness_of_fit.bic,
                        rmse: diagnostics.goodness_of_fit.rmse,
                        r_squared: diagnostics.goodness_of_fit.r_squared,
                        output_dir: method_output_dir,
                    }
                }
                EstimationMethod::Foce | EstimationMethod::FoceI => {
                    let mut estimator = FoceEstimator::new(model, config);
                    let results = estimator.fit(&dataset)?;
                    
                    // Convert FOCE results to SAEM format for diagnostics compatibility
                    let saem_results = convert_foce_to_saem_results(&results);
                    let diagnostics = diagnostics::generate_diagnostics(&dataset, &saem_results)?;
                    
                    // Save FOCE results
                    save_foce_results(&method_output_dir, &results, &diagnostics, &dataset, estimator.model())?;
                    
                    AnalysisResult {
                        model_type: model_type.clone(),
                        estimation_method: estimation_method.clone(),
                        objective_function_value: results.objective_function_value,
                        final_log_likelihood: results.final_log_likelihood,
                        converged: results.converged,
                        n_iterations: results.n_iterations,
                        fixed_effects: results.fixed_effects.clone(),
                        parameter_names: results.parameter_names.clone(),
                        aic: diagnostics.goodness_of_fit.aic,
                        bic: diagnostics.goodness_of_fit.bic,
                        rmse: diagnostics.goodness_of_fit.rmse,
                        r_squared: diagnostics.goodness_of_fit.r_squared,
                        output_dir: method_output_dir,
                    }
                }
            };
            
            all_results.push(analysis_result);
        }
    }

    // Generate comparison report if requested or if multiple analyses were run
    if args.compare_results || all_results.len() > 1 {
        generate_comparison_report(&args.output_dir, &all_results)?;
    }


    info!("Analysis completed successfully!");
    println!("Results saved to: {:?}", args.output_dir);
    
    // Print summary of all results
    println!("\nAnalysis Summary:");
    println!("{:<15} {:<10} {:<12} {:<12} {:<10} {:<8}", 
             "Model", "Method", "OFV", "LogLik", "Converged", "AIC");
    println!("{}", "-".repeat(80));
    
    for result in &all_results {
        println!("{:<15} {:<10} {:<12.2} {:<12.2} {:<10} {:<8.1}", 
                 format!("{}", result.model_type),
                 format!("{}", result.estimation_method),
                 result.objective_function_value,
                 result.final_log_likelihood,
                 result.converged,
                 result.aic);
    }
    
    // Identify best model by AIC
    if let Some(best_result) = all_results.iter().min_by(|a, b| a.aic.partial_cmp(&b.aic).unwrap()) {
        println!("\nBest model by AIC: {} with {} (AIC: {:.2})", 
                 best_result.model_type, best_result.estimation_method, best_result.aic);
    }

    Ok(())
}

#[derive(Debug, Clone)]
struct AnalysisResult {
    model_type: ModelType,
    estimation_method: EstimationMethod,
    objective_function_value: f64,
    final_log_likelihood: f64,
    converged: bool,
    n_iterations: usize,
    fixed_effects: Vec<f64>,
    parameter_names: Vec<String>,
    aic: f64,
    bic: f64,
    rmse: f64,
    r_squared: f64,
    output_dir: PathBuf,
}

fn generate_comparison_report(
    output_dir: &Path,
    results: &[AnalysisResult],
) -> Result<()> {
    let comparison_file = output_dir.join("model_comparison_report.txt");
    let mut report = String::new();
    
    report.push_str("NMODES Model and Method Comparison Report\n");
    report.push_str("=========================================\n\n");
    
    report.push_str(&format!("Total analyses performed: {}\n", results.len()));
    report.push_str(&format!("Analysis date: {}\n\n", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")));
    
    // Summary table
    report.push_str("Summary Table:\n");
    report.push_str("--------------\n");
    report.push_str(&format!("{:<15} {:<10} {:<12} {:<12} {:<10} {:<8} {:<8} {:<8} {:<8}\n", 
                             "Model", "Method", "OFV", "LogLik", "Converged", "AIC", "BIC", "RMSE", "R²"));
    report.push_str(&format!("{}\n", "-".repeat(100)));
    
    for result in results {
        report.push_str(&format!("{:<15} {:<10} {:<12.2} {:<12.2} {:<10} {:<8.1} {:<8.1} {:<8.3} {:<8.3}\n", 
                                 format!("{}", result.model_type),
                                 format!("{}", result.estimation_method),
                                 result.objective_function_value,
                                 result.final_log_likelihood,
                                 result.converged,
                                 result.aic,
                                 result.bic,
                                 result.rmse,
                                 result.r_squared));
    }
    
    // Model ranking by AIC
    let mut sorted_results = results.to_vec();
    sorted_results.sort_by(|a, b| a.aic.partial_cmp(&b.aic).unwrap());
    
    report.push_str("\nModel Ranking by AIC (lower is better):\n");
    report.push_str("--------------------------------------\n");
    for (rank, result) in sorted_results.iter().enumerate() {
        let delta_aic = result.aic - sorted_results[0].aic;
        report.push_str(&format!("{}. {} + {} (AIC: {:.2}, ΔAIC: {:.2})\n", 
                                 rank + 1,
                                 result.model_type,
                                 result.estimation_method,
                                 result.aic,
                                 delta_aic));
    }
    
    // Parameter comparison for converged models
    let converged_results: Vec<&AnalysisResult> = results.iter()
        .filter(|r| r.converged)
        .collect();
    
    if !converged_results.is_empty() {
        report.push_str("\nParameter Estimates (Converged Models Only):\n");
        report.push_str("-------------------------------------------\n");
        
        // Group by model type for easier comparison
        let mut models_by_type: std::collections::HashMap<String, Vec<&AnalysisResult>> = std::collections::HashMap::new();
        for result in &converged_results {
            models_by_type.entry(format!("{}", result.model_type))
                .or_default()
                .push(result);
        }
        
        for (model_name, model_results) in models_by_type {
            report.push_str(&format!("\n{} Model:\n", model_name));
            
            // Get parameter names (should be same for all results of same model type)
            if let Some(first_result) = model_results.first() {
                for (param_idx, param_name) in first_result.parameter_names.iter().enumerate() {
                    report.push_str(&format!("  {}:\n", param_name));
                    for result in &model_results {
                        let param_value = result.fixed_effects[param_idx].exp(); // Transform back from log scale
                        report.push_str(&format!("    {}: {:.4}\n", result.estimation_method, param_value));
                    }
                }
            }
        }
    }
    
    // Method comparison
    report.push_str("\nMethod Performance Comparison:\n");
    report.push_str("-----------------------------\n");
    
    let mut method_stats: std::collections::HashMap<String, Vec<&AnalysisResult>> = std::collections::HashMap::new();
    for result in results {
        method_stats.entry(format!("{}", result.estimation_method))
            .or_default()
            .push(result);
    }
    
    for (method_name, method_results) in method_stats {
        let converged_count = method_results.iter().filter(|r| r.converged).count();
        let avg_iterations = method_results.iter()
            .map(|r| r.n_iterations)
            .sum::<usize>() as f64 / method_results.len() as f64;
        
        report.push_str(&format!("{}:\n", method_name));
        report.push_str(&format!("  Convergence rate: {}/{} ({:.1}%)\n", 
                                 converged_count, method_results.len(),
                                 (converged_count as f64 / method_results.len() as f64) * 100.0));
        report.push_str(&format!("  Average iterations: {:.0}\n", avg_iterations));
        
        if converged_count > 0 {
            let avg_aic = method_results.iter()
                .filter(|r| r.converged)
                .map(|r| r.aic)
                .sum::<f64>() / converged_count as f64;
            report.push_str(&format!("  Average AIC (converged): {:.2}\n", avg_aic));
        }
        report.push_str("\n");
    }
    
    // Recommendations
    report.push_str("Recommendations:\n");
    report.push_str("---------------\n");
    
    if let Some(best_result) = sorted_results.first() {
        if best_result.converged {
            report.push_str(&format!("• Best fitting model: {} with {} (AIC: {:.2})\n", 
                                     best_result.model_type, best_result.estimation_method, best_result.aic));
        } else {
            report.push_str("• Warning: Best AIC model did not converge. Consider:\n");
            report.push_str("  - Increasing iterations\n");
            report.push_str("  - Trying different estimation method\n");
            report.push_str("  - Checking data quality\n");
        }
    }
    
    // Check for substantial model differences
    if sorted_results.len() > 1 {
        let delta_aic = sorted_results[1].aic - sorted_results[0].aic;
        if delta_aic < 2.0 {
            report.push_str("• Models have similar fit (ΔAIC < 2). Consider simpler model for parsimony.\n");
        } else if delta_aic > 10.0 {
            report.push_str("• Strong evidence for best model (ΔAIC > 10).\n");
        } else {
            report.push_str("• Moderate evidence for best model (2 < ΔAIC < 10).\n");
        }
    }
    
    // Method-specific recommendations
    let saem_results: Vec<&AnalysisResult> = results.iter()
        .filter(|r| matches!(r.estimation_method, EstimationMethod::Saem))
        .collect();
    let foce_results: Vec<&AnalysisResult> = results.iter()
        .filter(|r| matches!(r.estimation_method, EstimationMethod::Foce | EstimationMethod::FoceI))
        .collect();
    
    if !saem_results.is_empty() && !foce_results.is_empty() {
        report.push_str("• Method comparison available - check consistency between SAEM and FOCE results.\n");
    }
    
    fs::write(comparison_file, report)?;
    
    // Also generate CSV comparison for easy analysis
    generate_comparison_csv(output_dir, results)?;
    
    println!("Comparison report saved to: {:?}", output_dir.join("model_comparison_report.txt"));
    println!("Comparison CSV saved to: {:?}", output_dir.join("model_comparison.csv"));

    Ok(())
}

fn generate_comparison_csv(
    output_dir: &Path,
    results: &[AnalysisResult],
) -> Result<()> {
    let csv_file = output_dir.join("model_comparison.csv");
    let mut wtr = csv::Writer::from_path(csv_file)?;
    
    // Write header
    wtr.write_record(&[
        "Model", "Method", "OFV", "LogLikelihood", "Converged", 
        "Iterations", "AIC", "BIC", "RMSE", "R_squared"
    ])?;
    
    // Write data
    for result in results {
        wtr.write_record(&[
            format!("{}", result.model_type),
            format!("{}", result.estimation_method),
            result.objective_function_value.to_string(),
            result.final_log_likelihood.to_string(),
            result.converged.to_string(),
            result.n_iterations.to_string(),
            result.aic.to_string(),
            result.bic.to_string(),
            result.rmse.to_string(),
            result.r_squared.to_string(),
        ])?;
    }
    
    wtr.flush()?;
    Ok(())
}
fn convert_foce_to_saem_results(foce_results: &FoceResults) -> SaemResults {
    let mut saem_results = SaemResults::new(
        foce_results.fixed_effects.len(),
        foce_results.parameter_names.clone(),
    );
    
    saem_results.fixed_effects = foce_results.fixed_effects.clone();
    saem_results.random_effects_variance = foce_results.random_effects_variance.clone();
    saem_results.residual_variance = foce_results.residual_variance;
    saem_results.final_log_likelihood = foce_results.final_log_likelihood;
    saem_results.objective_function_value = foce_results.objective_function_value;
    saem_results.converged = foce_results.converged;
    saem_results.n_iterations = foce_results.n_iterations;
    saem_results.individual_parameters = foce_results.individual_parameters.clone();
    
    saem_results
}

fn save_foce_results(
    output_dir: &std::path::Path,
    results: &FoceResults,
    diagnostics: &crate::diagnostics::DiagnosticResults,
    dataset: &Dataset,
    model: &CompartmentModel,
) -> Result<()> {
    use std::fs;
    
    // Ensure output directory exists
    fs::create_dir_all(output_dir)?;
    
    // Save FOCE-specific results
    let foce_file = output_dir.join("foce_results.json");
    let json_content = serde_json::to_string_pretty(results)?;
    fs::write(foce_file, json_content)?;
    
    // Save diagnostics
    let diagnostics_file = output_dir.join("diagnostics.json");
    let json_content = serde_json::to_string_pretty(diagnostics)?;
    fs::write(diagnostics_file, json_content)?;
    
    // Save FOCE-specific summary report
    save_foce_summary_report(output_dir, results, diagnostics)?;
    
    // Save predictions using FOCE results
    save_foce_predictions_csv(output_dir, results, dataset, model)?;
    
    Ok(())
}

fn save_foce_summary_report(
    output_dir: &std::path::Path,
    results: &FoceResults,
    diagnostics: &crate::diagnostics::DiagnosticResults,
) -> Result<()> {
    use std::fs;
    
    let report_file = output_dir.join("foce_summary_report.txt");
    
    let mut report = String::new();
    report.push_str("NMODES FOCE Analysis Summary Report\n");
    report.push_str("=================================\n\n");
    
    report.push_str(&format!("Estimation Method: FOCE\n"));
    report.push_str(&format!("Model Convergence: {}\n", results.converged));
    report.push_str(&format!("Total Iterations: {}\n", results.n_iterations));
    report.push_str(&format!("Final Log-Likelihood: {:.6}\n", results.final_log_likelihood));
    report.push_str(&format!("Objective Function Value: {:.6}\n", results.objective_function_value));
    report.push_str(&format!("Gradient Norm: {:.6}\n", results.gradient_norm));
    report.push_str(&format!("Hessian Condition Number: {:.6}\n", results.hessian_condition_number));
    report.push_str(&format!("Number of Individuals: {}\n", results.individual_parameters.len()));
    report.push_str(&format!("AIC: {:.6}\n", diagnostics.goodness_of_fit.aic));
    report.push_str(&format!("BIC: {:.6}\n", diagnostics.goodness_of_fit.bic));
    report.push_str(&format!("R-squared: {:.6}\n", diagnostics.goodness_of_fit.r_squared));
    report.push_str(&format!("RMSE: {:.6}\n", diagnostics.goodness_of_fit.rmse));
    
    report.push_str("\nFixed Effects Parameter Estimates:\n");
    report.push_str("----------------------------------\n");
    report.push_str(&format!("{:<10} {:<12} {:<10}\n", "Parameter", "Estimate", "SE"));
    report.push_str(&format!("{:<10} {:<12} {:<10}\n", "---------", "--------", "--"));
    
    for (i, param_name) in results.parameter_names.iter().enumerate() {
        let estimate = results.fixed_effects[i];
        let se = results.standard_errors.get(i).copied().unwrap_or(0.0);
        report.push_str(&format!("{:<10} {:<12.6} {:<10.6}\n", param_name, estimate, se));
    }
    
    report.push_str(&format!("\nResidual Error Variance: {:.6}\n", results.residual_variance));
    
    report.push_str("\nRandom Effects Variance (Omega):\n");
    report.push_str("-------------------------------\n");
    for i in 0..results.parameter_names.len() {
        let param_name = &results.parameter_names[i];
        let variance = results.random_effects_variance[i][i];
        report.push_str(&format!("{}({}): {:.6}\n", param_name, param_name, variance));
    }
    
    fs::write(report_file, report)?;
    Ok(())
}

fn save_foce_predictions_csv(
    output_dir: &std::path::Path,
    results: &nmodes::FoceResults,
    dataset: &Dataset,
    model: &CompartmentModel,
) -> Result<()> {
    
    let predictions_file = output_dir.join("foce_predictions.csv");
    let mut wtr = csv::Writer::from_path(predictions_file)?;
    
    // Write header
    wtr.write_record(&["ID", "TIME", "DV", "IPRED", "PRED"])?;
    
    let solver = RungeKuttaSolver::new();
    let solver_config = SolverConfig::default();
    // Create a single default vector to be borrowed if an individual's eta is missing.
    let default_eta_value = vec![0.0; results.fixed_effects.len()];
    for (&id, individual) in dataset.individuals() {
        // Get individual parameters (theta + eta)
        // FIX: Borrow the pre-allocated default_eta_value instead of a temporary vector.
        let ind_eta = results
            .individual_parameters
            .get(&id)
            .unwrap_or(&default_eta_value);
        
        let mut ind_params = model.default_parameters();
        for i in 0..results.fixed_effects.len() {
            ind_params.fixed_effects[i] = results.fixed_effects[i] + ind_eta[i];
        }
        
        // Population parameters (theta only)
        let mut pop_params = model.default_parameters();
        pop_params.fixed_effects = results.fixed_effects.clone();
        
        // Calculate predictions (simplified version)
        for obs in individual.observations() {
            // For now, use a simplified prediction
            let ipred = results.fixed_effects[0].exp(); // Simplified
            let pred = results.fixed_effects[0].exp();  // Simplified
            
            wtr.write_record(&[
                id.to_string(),
                obs.time.to_string(),
                obs.value.to_string(),
                ipred.to_string(),
                pred.to_string(),
            ])?;
        }
    }
    
    wtr.flush()?;
    Ok(())
}