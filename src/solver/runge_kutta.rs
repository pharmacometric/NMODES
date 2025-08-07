use super::{OdeSolver, OdeSystem, SolverConfig, SolverError};
use nalgebra::DVector;

pub struct RungeKuttaSolver;

impl RungeKuttaSolver {
    pub fn new() -> Self {
        Self
    }

    fn rk4_step(
        &self,
        system: &dyn OdeSystem,
        t: f64,
        y: &DVector<f64>,
        h: f64,
    ) -> DVector<f64> {
        let k1 = system.derivatives(t, y);
        let k2 = system.derivatives(t + h / 2.0, &(y + &k1 * (h / 2.0)));
        let k3 = system.derivatives(t + h / 2.0, &(y + &k2 * (h / 2.0)));
        let k4 = system.derivatives(t + h, &(y + &k3 * h));
        
        y + (&k1 + &k2 * 2.0 + &k3 * 2.0 + &k4) * (h / 6.0)
    }
}

impl OdeSolver for RungeKuttaSolver {
    fn solve(
        &self,
        system: &dyn OdeSystem,
        t_span: (f64, f64),
        y0: &DVector<f64>,
        config: &SolverConfig,
    ) -> Result<(Vec<f64>, Vec<DVector<f64>>), SolverError> {
        let dt = t_span.1 - t_span.0;
        if dt <= 0.0 {
            return Err(SolverError::InvalidTimeStep(dt));
        }

        // Determine number of steps
        let n_steps = ((dt / config.max_step_size).ceil() as usize).max(1);
        let step_size = dt / n_steps as f64;

        let mut times = Vec::with_capacity(n_steps + 1);
        let mut solutions = Vec::with_capacity(n_steps + 1);
        
        let mut t = t_span.0;
        let mut y = y0.clone();
        
        times.push(t);
        solutions.push(y.clone());
        
        for _ in 0..n_steps {
            y = self.rk4_step(system, t, &y, step_size);
            t += step_size;
            
            // Check for numerical issues
            if y.as_slice().iter().any(|&val| !val.is_finite()) {
                return Err(SolverError::NumericalInstability);
            }
            
            times.push(t);
            solutions.push(y.clone());
        }
        
        Ok((times, solutions))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct TestSystem;
    
    impl OdeSystem for TestSystem {
        fn derivatives(&self, _t: f64, y: &DVector<f64>) -> DVector<f64> {
            // Simple exponential decay: dy/dt = -y
            -y.clone()
        }
        
        fn dimension(&self) -> usize {
            1
        }
    }

    #[test]
    fn test_runge_kutta_solver() {
        let solver = RungeKuttaSolver::new();
        let system = TestSystem;
        let y0 = DVector::from_vec(vec![1.0]);
        let config = SolverConfig::default();
        
        let result = solver.solve(&system, (0.0, 1.0), &y0, &config);
        assert!(result.is_ok());
        
        let (times, solutions) = result.unwrap();
        assert!(!times.is_empty());
        assert_eq!(times.len(), solutions.len());
        
        // Final solution should be approximately e^(-1) ≈ 0.368
        let final_solution = solutions.last().unwrap()[0];
        assert!((final_solution - (-1.0_f64).exp()).abs() < 0.01);
    }
}