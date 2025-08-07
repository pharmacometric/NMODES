use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum DosingType {
    Bolus,
    Infusion,
    Oral,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DosingRecord {
    pub time: f64,
    pub amount: f64,
    pub compartment: i32,
    pub rate: Option<f64>,
    pub dosing_type: DosingType,
    pub additional_doses: i32,
    pub interdose_interval: Option<f64>,
    pub steady_state: bool,
}

impl DosingRecord {
    pub fn new(
        time: f64,
        amount: f64,
        compartment: i32,
        dosing_type: DosingType,
    ) -> Self {
        Self {
            time,
            amount,
            compartment,
            rate: None,
            dosing_type,
            additional_doses: 0,
            interdose_interval: None,
            steady_state: false,
        }
    }

    pub fn is_valid(&self) -> bool {
        self.time >= 0.0 && 
        self.amount > 0.0 && 
        self.compartment > 0 &&
        self.additional_doses >= 0
    }

    pub fn infusion_duration(&self) -> Option<f64> {
        if matches!(self.dosing_type, DosingType::Infusion) {
            self.rate.map(|r| self.amount / r)
        } else {
            None
        }
    }

    pub fn expand_multiple_doses(&self) -> Vec<DosingRecord> {
        let mut doses = vec![self.clone()];
        
        if self.additional_doses > 0 {
            if let Some(ii) = self.interdose_interval {
                for i in 1..=self.additional_doses {
                    let mut dose = self.clone();
                    dose.time = self.time + (i as f64) * ii;
                    dose.additional_doses = 0; // Expanded doses are single
                    doses.push(dose);
                }
            }
        }
        
        doses
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dosing_record_creation() {
        let dose = DosingRecord::new(0.0, 100.0, 1, DosingType::Bolus);
        assert_eq!(dose.time, 0.0);
        assert_eq!(dose.amount, 100.0);
        assert!(dose.is_valid());
    }

    #[test]
    fn test_multiple_dose_expansion() {
        let mut dose = DosingRecord::new(0.0, 100.0, 1, DosingType::Oral);
        dose.additional_doses = 2;
        dose.interdose_interval = Some(12.0);
        
        let expanded = dose.expand_multiple_doses();
        assert_eq!(expanded.len(), 3);
        assert_eq!(expanded[0].time, 0.0);
        assert_eq!(expanded[1].time, 12.0);
        assert_eq!(expanded[2].time, 24.0);
    }
}