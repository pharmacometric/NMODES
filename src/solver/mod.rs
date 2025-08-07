pub mod ode;
pub mod runge_kutta;

pub use ode::{OdeSolver, OdeSystem, SolverConfig};
pub use runge_kutta::RungeKuttaSolver;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum SolverError {
    #[error("Integration failed: {0}")]
    IntegrationFailed(String),
    
    #[error("Invalid time step: {0}")]
    InvalidTimeStep(f64),
    
    #[error("Maximum iterations exceeded")]
    MaxIterationsExceeded,
    
    #[error("Numerical instability detected")]
    NumericalInstability,
}