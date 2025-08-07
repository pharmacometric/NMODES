use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ObservationType {
    Concentration,
    Effect,
    Missing,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Observation {
    pub time: f64,
    pub value: f64,
    pub compartment: i32,
    pub observation_type: ObservationType,
}

impl Observation {
    pub fn new(time: f64, value: f64, compartment: i32, observation_type: ObservationType) -> Self {
        Self {
            time,
            value,
            compartment,
            observation_type,
        }
    }

    pub fn is_valid(&self) -> bool {
        self.time >= 0.0 && 
        self.value.is_finite() && 
        self.compartment > 0 &&
        !matches!(self.observation_type, ObservationType::Missing)
    }

    pub fn log_concentration(&self) -> Option<f64> {
        if self.value > 0.0 && matches!(self.observation_type, ObservationType::Concentration) {
            Some(self.value.ln())
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_observation_creation() {
        let obs = Observation::new(1.0, 10.0, 1, ObservationType::Concentration);
        assert_eq!(obs.time, 1.0);
        assert_eq!(obs.value, 10.0);
        assert!(obs.is_valid());
    }

    #[test]
    fn test_log_concentration() {
        let obs = Observation::new(1.0, 10.0, 1, ObservationType::Concentration);
        assert_eq!(obs.log_concentration(), Some(10.0_f64.ln()));
        
        let zero_obs = Observation::new(1.0, 0.0, 1, ObservationType::Concentration);
        assert_eq!(zero_obs.log_concentration(), None);
    }
}