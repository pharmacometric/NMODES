pub mod dataset;
pub mod observation;
pub mod dosing;
pub mod individual;

pub use dataset::Dataset;
pub use observation::{Observation, ObservationType};
pub use dosing::{DosingRecord, DosingType};
pub use individual::Individual;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum DataError {
    #[error("CSV parsing error: {0}")]
    CsvError(#[from] csv::Error),
    
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    
    #[error("Invalid data format: {0}")]
    InvalidFormat(String),
    
    #[error("Missing required column: {0}")]
    MissingColumn(String),
    
    #[error("Invalid time sequence for individual {0}")]
    InvalidTimeSequence(i32),
    
    #[error("No observations found for individual {0}")]
    NoObservations(i32),
    
    #[error("Invalid dose amount: {0}")]
    InvalidDose(f64),
    
    #[error("Negative time value: {0}")]
    NegativeTime(f64),
}