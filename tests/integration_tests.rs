use nmodes::data::Dataset;
use nmodes::models::{CompartmentModel, ModelType};
use nmodes::saem::SaemEstimator;
use nmodes::estimation::EstimationConfig;
use std::path::PathBuf;

#[test]
fn test_full_pipeline() {
    // Create test dataset
    let dataset_path = PathBuf::from("examples/example_dataset.csv");
    
    // This test requires the example dataset to exist
    if !dataset_path.exists() {
        println!("Skipping integration test - example dataset not found");
        return;
    }
    
    let dataset = Dataset::from_csv(&dataset_path).expect("Failed to load dataset");
    assert!(dataset.n_individuals() > 0);
    assert!(dataset.n_observations() > 0);
    
    // Create one-compartment model
    let model = CompartmentModel::new(ModelType::OneCompartment)
        .expect("Failed to create model");
    
    // Configure estimation (reduced iterations for testing)
    let config = EstimationConfig::default()
        .with_iterations(10)
        .with_burnin(2);
    
    // Run estimation
    let mut estimator = SaemEstimator::new(model, config);
    let results = estimator.fit(&dataset).expect("Estimation failed");
    
    // Basic checks
    assert_eq!(results.fixed_effects.len(), 2); // CL and V for 1-compartment
    assert!(results.n_iterations > 0);
    assert!(!results.log_likelihood_trajectory.is_empty());
}

#[test]
fn test_model_comparison() {
    // Test that different models can be created and have different characteristics
    let models = vec![
        ModelType::OneCompartment,
        ModelType::TwoCompartment,
        ModelType::ThreeCompartment,
    ];
    
    for model_type in models {
        let model = CompartmentModel::new(model_type.clone())
            .expect("Failed to create model");
        
        let expected_compartments = match model_type {
            ModelType::OneCompartment => 1,
            ModelType::TwoCompartment => 2,
            ModelType::ThreeCompartment => 3,
        };
        
        assert_eq!(model.n_compartments(), expected_compartments);
        
        let params = model.default_parameters();
        assert!(model.validate_parameters(&params).is_ok());
    }
}