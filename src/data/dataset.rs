use super::{DataError, Individual, Observation, DosingRecord, ObservationType, DosingType};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

unsafe impl Send for Dataset {}
unsafe impl Sync for Dataset {}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NonmemRecord {
    #[serde(rename = "ID")]
    pub id: i32,
    #[serde(rename = "TIME")]
    pub time: f64,
    #[serde(rename = "DV")]
    pub dv: Option<f64>,
    #[serde(rename = "AMT")]
    pub amt: Option<f64>,
    #[serde(rename = "EVID")]
    pub evid: i32,
    #[serde(rename = "CMT")]
    pub cmt: Option<i32>,
    #[serde(rename = "RATE")]
    pub rate: Option<f64>,
    #[serde(rename = "II")]
    pub ii: Option<f64>,
    #[serde(rename = "ADDL")]
    pub addl: Option<i32>,
    #[serde(rename = "SS")]
    pub ss: Option<i32>,
}

#[derive(Debug, Clone)]
pub struct Dataset {
    individuals: HashMap<i32, Individual>,
    covariate_names: Vec<String>,
}

impl Dataset {
    pub fn from_csv<P: AsRef<Path>>(path: P) -> Result<Self, DataError> {
        let mut reader = csv::Reader::from_path(path)?;
        let headers = reader.headers()?.clone();
        
        // Validate required columns
        let required_cols = ["ID", "TIME", "DV", "AMT", "EVID"];
        for col in required_cols.iter() {
            if !headers.iter().any(|h| h == *col) {
                return Err(DataError::MissingColumn(col.to_string()));
            }
        }

        let mut individuals: HashMap<i32, Individual> = HashMap::new();
        let mut records_by_id: HashMap<i32, Vec<NonmemRecord>> = HashMap::new();

        // Parse all records
        for result in reader.deserialize() {
            let record: NonmemRecord = result?;
            
            // Validate basic constraints
            if record.time < 0.0 {
                return Err(DataError::NegativeTime(record.time));
            }
            
            if let Some(amt) = record.amt {
                if amt < 0.0 {
                    return Err(DataError::InvalidDose(amt));
                }
            }

            records_by_id.entry(record.id).or_default().push(record);
        }

        // Process records into individuals
        for (id, mut records) in records_by_id {
            // Sort by time
            records.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());
            
            let individual = Self::process_individual_records(id, records)?;
            individuals.insert(id, individual);
        }

        if individuals.is_empty() {
            return Err(DataError::InvalidFormat("No valid individuals found".to_string()));
        }

        // Extract covariate names (columns not in standard NONMEM set)
        let standard_cols = ["ID", "TIME", "DV", "AMT", "EVID", "CMT", "RATE", "II", "ADDL", "SS"];
        let covariate_names: Vec<String> = headers.iter()
            .filter(|h| !standard_cols.contains(h))
            .map(|h| h.to_string())
            .collect();

        Ok(Dataset {
            individuals,
            covariate_names,
        })
    }

    fn process_individual_records(
        id: i32, 
        records: Vec<NonmemRecord>
    ) -> Result<Individual, DataError> {
        let mut observations = Vec::new();
        let mut dosing_records = Vec::new();

        for record in records {
            match record.evid {
                0 => {
                    // Observation record
                    if let Some(dv) = record.dv {
                        let obs = Observation {
                            time: record.time,
                            value: dv,
                            compartment: record.cmt.unwrap_or(1),
                            observation_type: if dv > 0.0 { 
                                ObservationType::Concentration 
                            } else { 
                                ObservationType::Missing 
                            },
                        };
                        observations.push(obs);
                    }
                }
                1 => {
                    // Dosing record
                    if let Some(amt) = record.amt {
                        let dose = DosingRecord {
                            time: record.time,
                            amount: amt,
                            compartment: record.cmt.unwrap_or(1),
                            rate: record.rate,
                            dosing_type: if record.rate.is_some() && record.rate.unwrap() > 0.0 {
                                DosingType::Infusion
                            } else {
                                DosingType::Bolus
                            },
                            additional_doses: record.addl.unwrap_or(0),
                            interdose_interval: record.ii,
                            steady_state: record.ss.unwrap_or(0) == 1,
                        };
                        dosing_records.push(dose);
                    }
                }
                _ => {
                    // Other event types (reset, etc.)
                    continue;
                }
            }
        }

        if observations.is_empty() {
            return Err(DataError::NoObservations(id));
        }

        // Validate time sequence
        for i in 1..observations.len() {
            if observations[i].time < observations[i-1].time {
                return Err(DataError::InvalidTimeSequence(id));
            }
        }

        Ok(Individual::new(id, observations, dosing_records, HashMap::new()))
    }

    pub fn individuals(&self) -> &HashMap<i32, Individual> {
        &self.individuals
    }

    pub fn n_individuals(&self) -> usize {
        self.individuals.len()
    }

    pub fn n_observations(&self) -> usize {
        self.individuals.values()
            .map(|ind| ind.observations().len())
            .sum()
    }

    pub fn covariate_names(&self) -> &[String] {
        &self.covariate_names
    }

    pub fn get_individual(&self, id: i32) -> Option<&Individual> {
        self.individuals.get(&id)
    }

    pub fn get_all_times(&self) -> Vec<f64> {
    let mut times: Vec<f64> = self.individuals.values()
        .flat_map(|ind| ind.observations().iter().map(|obs| obs.time))
        .collect();
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    times.dedup_by(|a, b| (*a - *b).abs() < f64::EPSILON);
    times
}

    pub fn get_concentration_data(&self) -> Vec<(f64, f64)> {
        self.individuals.values()
            .flat_map(|ind| {
                ind.observations().iter()
                    .filter(|obs| matches!(obs.observation_type, ObservationType::Concentration))
                    .map(|obs| (obs.time, obs.value))
            })
            .collect()
    }
}