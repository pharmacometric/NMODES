use clap::{Arg, Command};
use log::{info, warn, error};
use std::path::PathBuf;
use nmodes::{Dataset, CompartmentModel, ModelType, SaemEstimator, RungeKuttaSolver, SolverConfig};
use nmodes::{EstimationConfig, EstimationMethod, FoceEstimator};
use nmodes::{diagnostics, output, validation};
use anyhow::{Result, anyhow};

#[derive(Debug)]
struct CliArgs {
    dataset_path: PathBuf,
    model_type: ModelType,
    estimation_method: EstimationMethod,
    output_dir: PathBuf,
    iterations: usize,
    burn_in: usize,
    chains: usize,
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
                .help("Compartment model type: 1comp, 2comp, or 3comp")
                .default_value("1comp")
        )
        .arg(
            Arg::new("method")
                .short('e')
                .long("method")
                .value_name("METHOD")
                .help("Estimation method: saem, foce, or foce-i")
                .default_value("saem")
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
        .get_matches();

    let args = CliArgs {
        dataset_path: PathBuf::from(matches.get_one::<String>("dataset").unwrap()),
        model_type: parse_model_type(matches.get_one::<String>("model").unwrap())?,
        estimation_method: parse_estimation_method(matches.get_one::<String>("method").unwrap())?,
        output_dir: PathBuf::from(matches.get_one::<String>("output").unwrap()),
        iterations: matches.get_one::<String>("iterations").unwrap().parse()?,
        burn_in: matches.get_one::<String>("burn-in").unwrap().parse()?,
        chains: matches.get_one::<String>("chains").unwrap().parse()?,
    };

    run_analysis(args)
}

fn parse_model_type(model_str: &str) -> Result<ModelType> {
    match model_str {
        "1comp" => Ok(ModelType::OneCompartment),
        "2comp" => Ok(ModelType::TwoCompartment),
        "3comp" => Ok(ModelType::ThreeCompartment),
        _ => Err(anyhow!("Invalid model type: {}", model_str)),
    }
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
    info!("Model type: {:?}", args.model_type);
    info!("Estimation method: {:?}", args.estimation_method);
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

    // Create model
    let model = CompartmentModel::new(args.model_type)?;
    info!("Created {} model", model.model_type());

    // Configure estimation
    let config = EstimationConfig {
        method: args.estimation_method.clone(),
        n_iterations: args.iterations,
        n_burnin: args.burn_in,
        n_chains: args.chains,
        step_size: 0.1,
        target_acceptance: 0.44,
        adaptation_interval: 50,
        foce_max_iterations: if matches!(args.estimation_method, EstimationMethod::Foce | EstimationMethod::FoceI) {
            args.iterations
        } else {
            100
        },
        foce_tolerance: 1e-6,
        foce_step_size: 1e-4,
        foce_interaction: matches!(args.estimation_method, EstimationMethod::FoceI),
        ..Default::default()
    };

    // Run SAEM estimation
    info!("Starting {} estimation with {} iterations...", args.estimation_method, args.iterations);
    
    let (results, diagnostics) = match args.estimation_method {
        EstimationMethod::Saem => {
            let mut estimator = SaemEstimator::new(model, config);
            let results = estimator.fit(&dataset)?;
            
            // Convert SAEM results to common format for diagnostics
            let diagnostics = diagnostics::generate_diagnostics(&dataset, &results)?;
            
            // Save SAEM results
            output::save_results(&args.output_dir, &results, &diagnostics, &dataset, estimator.model())?;
            
            (Box::new(results) as Box<dyn std::any::Any>, diagnostics)
        }
        EstimationMethod::Foce | EstimationMethod::FoceI => {
            let mut estimator = FoceEstimator::new(model, config);
            let results = estimator.fit(&dataset)?;
            
            // Convert FOCE results to SAEM format for diagnostics compatibility
            let saem_results = convert_foce_to_saem_results(&results);
            let diagnostics = diagnostics::generate_diagnostics(&dataset, &saem_results)?;
            
            // Save FOCE results
            save_foce_results(&args.output_dir, &results, &diagnostics, &dataset, estimator.model())?;
            
            (Box::new(results) as Box<dyn std::any::Any>, diagnostics)
        }
    };


    info!("Analysis completed successfully!");
    println!("Results saved to: {:?}", args.output_dir);
    println!("Estimation method: {}", args.estimation_method);
    
    // Print method-specific results
    match args.estimation_method {
        EstimationMethod::Saem => {
            if let Ok(saem_results) = results.downcast::<crate::saem::SaemResults>() {
                println!("Final log-likelihood: {:.3}", saem_results.final_log_likelihood);
                println!("Objective function value: {:.3}", saem_results.objective_function_value);
                println!("Convergence achieved: {}", saem_results.converged);
            }
        }
        EstimationMethod::Foce | EstimationMethod::FoceI => {
            if let Ok(foce_results) = results.downcast::<crate::estimation::FoceResults>() {
                println!("Final log-likelihood: {:.3}", foce_results.final_log_likelihood);
                println!("Objective function value: {:.3}", foce_results.objective_function_value);
                println!("Convergence achieved: {}", foce_results.converged);
            }
        }
    }

    Ok(())
}

fn convert_foce_to_saem_results(foce_results: &crate::estimation::FoceResults) -> crate::saem::SaemResults {
    let mut saem_results = crate::saem::SaemResults::new(
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
    results: &crate::estimation::FoceResults,
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
    results: &crate::estimation::FoceResults,
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
    use std::fs;
    
    let predictions_file = output_dir.join("foce_predictions.csv");
    let mut wtr = csv::Writer::from_path(predictions_file)?;
    
    // Write header
    wtr.write_record(&["ID", "TIME", "DV", "IPRED", "PRED"])?;
    
    let solver = RungeKuttaSolver::new();
    let solver_config = SolverConfig::default();
    
    for (&id, individual) in dataset.individuals() {
        // Get individual parameters (theta + eta)
        let ind_eta = results.individual_parameters.get(&id)
            .unwrap_or(&vec![0.0; results.fixed_effects.len()]);
        
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