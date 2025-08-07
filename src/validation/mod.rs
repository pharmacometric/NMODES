use crate::data::{Dataset, DataError};
use log::{info, warn};

pub fn validate_dataset(dataset: &Dataset) -> Result<(), DataError> {
    info!("Validating dataset with {} individuals", dataset.n_individuals());
    
    // Check minimum requirements
    if dataset.n_individuals() == 0 {
        return Err(DataError::InvalidFormat("No individuals in dataset".to_string()));
    }
    
    if dataset.n_observations() == 0 {
        return Err(DataError::InvalidFormat("No observations in dataset".to_string()));
    }
    
    let mut total_dose_events = 0;
    let mut individuals_with_doses = 0;
    let mut individuals_with_observations = 0;
    
    for (id, individual) in dataset.individuals() {
        // Validate individual data
        if individual.observations().is_empty() {
            warn!("Individual {} has no observations", id);
            continue;
        }
        individuals_with_observations += 1;
        
        if !individual.dosing_records().is_empty() {
            individuals_with_doses += 1;
            total_dose_events += individual.dosing_records().len();
        }
        
        // Check time ordering
        let obs_times = individual.observation_times();
        for i in 1..obs_times.len() {
            if obs_times[i] < obs_times[i-1] {
                return Err(DataError::InvalidTimeSequence(*id));
            }
        }
        
        // Check for reasonable concentration values
        for obs in individual.observations() {
            if obs.value < 0.0 {
                warn!("Individual {} has negative concentration at time {}", id, obs.time);
            }
            
            if obs.value > 1e6 {
                warn!("Individual {} has very high concentration ({}) at time {}", 
                      id, obs.value, obs.time);
            }
        }
    }
    
    info!("Dataset validation completed:");
    info!("  - {} individuals with observations", individuals_with_observations);
    info!("  - {} individuals with dosing records", individuals_with_doses);
    info!("  - {} total dose events", total_dose_events);
    
    if individuals_with_doses == 0 {
        warn!("No dosing information found in dataset");
    }
    
    Ok(())
}

pub fn validate_model_fit(
    predicted: &[f64],
    observed: &[f64],
) -> Result<(), String> {
    if predicted.len() != observed.len() {
        return Err("Predicted and observed vectors must have same length".to_string());
    }
    
    if predicted.is_empty() {
        return Err("No data to validate".to_string());
    }
    
    // Check for unreasonable predictions
    for (i, &pred) in predicted.iter().enumerate() {
        if !pred.is_finite() {
            return Err(format!("Non-finite prediction at index {}: {}", i, pred));
        }
        
        if pred < 0.0 {
            warn!("Negative prediction at index {}: {}", i, pred);
        }
    }
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::*;
    use std::collections::HashMap;

    #[test]
    fn test_validate_empty_dataset() {
        // Create a minimal test dataset using the public API
        let temp_file = std::env::temp_dir().join("empty_test.csv");
        std::fs::write(&temp_file, "ID,TIME,DV,AMT,EVID\n").unwrap();
        let dataset = Dataset::from_csv(&temp_file).unwrap_or_else(|_| {
            // If CSV parsing fails, we know it's empty, which is what we want to test
            panic!("Expected empty dataset error");
        });
        assert!(validate_dataset(&dataset).is_err());
        std::fs::remove_file(&temp_file).ok();
    }
}