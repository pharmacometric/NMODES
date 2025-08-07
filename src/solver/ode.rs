use super::SolverError;
use nalgebra::DVector;

pub trait OdeSystem {
    fn derivatives(&self, t: f64, y: &DVector<f64>) -> DVector<f64>;
    fn dimension(&self) -> usize;
}

#[derive(Debug, Clone)]
pub struct SolverConfig {
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub max_step_size: f64,
    pub min_step_size: f64,
    pub max_iterations: usize,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            absolute_tolerance: 1e-8,
            relative_tolerance: 1e-6,
            max_step_size: 1.0,
            min_step_size: 1e-12,
            max_iterations: 10000,
        }
    }
}

pub trait OdeSolver {
    fn solve(
        &self,
        system: &dyn OdeSystem,
        t_span: (f64, f64),
        y0: &DVector<f64>,
        config: &SolverConfig,
    ) -> Result<(Vec<f64>, Vec<DVector<f64>>), SolverError>;

    fn solve_to_time(
        &self,
        system: &dyn OdeSystem,
        t_start: f64,
        t_end: f64,
        y0: &DVector<f64>,
        config: &SolverConfig,
    ) -> Result<DVector<f64>, SolverError> {
        let (_, solutions) = self.solve(system, (t_start, t_end), y0, config)?;
        Ok(solutions.into_iter().last().unwrap_or_else(|| y0.clone()))
    }
}
