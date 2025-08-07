use super::{Observation, DosingRecord};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Individual {
    pub id: i32,
    observations: Vec<Observation>,
    dosing_records: Vec<DosingRecord>,
    covariates: HashMap<String, f64>,
}

impl Individual {
    pub fn new(
        id: i32,
        observations: Vec<Observation>,
        dosing_records: Vec<DosingRecord>,
        covariates: HashMap<String, f64>,
    ) -> Self {
        Self {
            id,
            observations,
            dosing_records,
            covariates,
        }
    }

    pub fn observations(&self) -> &[Observation] {
        &self.observations
    }

    pub fn dosing_records(&self) -> &[DosingRecord] {
        &self.dosing_records
    }

    pub fn covariates(&self) -> &HashMap<String, f64> {
        &self.covariates
    }

    pub fn n_observations(&self) -> usize {
        self.observations.len()
    }

    pub fn observation_times(&self) -> Vec<f64> {
        self.observations.iter().map(|obs| obs.time).collect()
    }

    pub fn concentration_values(&self) -> Vec<f64> {
        self.observations.iter()
            .filter_map(|obs| {
                if obs.observation_type == super::ObservationType::Concentration {
                    Some(obs.value)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn get_covariate(&self, name: &str) -> Option<f64> {
        self.covariates.get(name).copied()
    }

    pub fn set_covariate(&mut self, name: String, value: f64) {
        self.covariates.insert(name, value);
    }

    pub fn first_dose_time(&self) -> Option<f64> {
        self.dosing_records.first().map(|dose| dose.time)
    }

    pub fn last_observation_time(&self) -> Option<f64> {
        self.observations.last().map(|obs| obs.time)
    }

    pub fn total_dose(&self) -> f64 {
        self.dosing_records.iter()
            .flat_map(|dose| dose.expand_multiple_doses())
            .map(|dose| dose.amount)
            .sum()
    }

    pub fn baseline_measurement(&self) -> Option<f64> {
        self.observations.first()
            .filter(|obs| obs.time == 0.0)
            .map(|obs| obs.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{ObservationType, DosingType};

    #[test]
    fn test_individual_creation() {
        let obs = vec![
            Observation::new(0.0, 0.0, 1, ObservationType::Concentration),
            Observation::new(1.0, 10.0, 1, ObservationType::Concentration),
        ];
        let doses = vec![
            DosingRecord::new(0.0, 100.0, 1, DosingType::Bolus),
        ];
        
        let individual = Individual::new(1, obs, doses, HashMap::new());
        assert_eq!(individual.id, 1);
        assert_eq!(individual.n_observations(), 2);
        assert_eq!(individual.total_dose(), 100.0);
    }
}