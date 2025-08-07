use super::compartment::{CompartmentModelTrait, ModelParameters, ModelState};
use super::ModelError;
use nalgebra::DVector;

pub struct OneCompartmentModel {
    // Model: dA/dt = -CL/V * A
    // Where A is amount in compartment, CL is clearance, V is volume
}

impl OneCompartmentModel {
    pub fn new() -> Self {
        Self {}
    }
}

impl CompartmentModelTrait for OneCompartmentModel {
    fn n_compartments(&self) -> usize {
        1
    }

    fn parameter_names(&self) -> Vec<String> {
        vec!["CL".to_string(), "V".to_string()]
    }

    fn default_parameters(&self) -> ModelParameters {
        let param_names = self.parameter_names();
        let mut params = ModelParameters::new(2, param_names);
        
        // Typical values for a one-compartment model
        params.fixed_effects[0] = 1.0_f64.ln(); // ln(CL) = ln(1.0 L/h) = 0.0
        params.fixed_effects[1] = 3.0_f64.ln(); // ln(V) = ln(20 L) ≈ 2.996
        
        // Inter-individual variability (diagonal omega matrix)
        params.random_effects_variance[0][0] = 0.09; // 30% CV for CL
        params.random_effects_variance[1][1] = 0.04; // 20% CV for V
        
        // Residual error (proportional)
        params.residual_variance = 0.01; // 10% CV
        
        params
    }

    fn derivatives(&self, state: &ModelState, params: &ModelParameters) -> DVector<f64> {
        let cl = params.fixed_effects[0].exp();
        let v = params.fixed_effects[1].exp();
        
        let ke = cl / v; // Elimination rate constant
        let mut derivatives = DVector::<f64>::zeros(1);
        
        // dA/dt = -ke * A
        derivatives[0] = -ke * state.compartments[0];
        
        derivatives
    }

    fn observation_function(&self, state: &ModelState, params: &ModelParameters, compartment: usize) -> f64 {
        if compartment != 1 {
            return 0.0;
        }
        
        let v = params.fixed_effects[1].exp();
        
        // Concentration = Amount / Volume
        state.compartments[0] / v
    }

    fn validate_parameters(&self, params: &ModelParameters) -> Result<(), ModelError> {
        if params.n_parameters() != 2 {
            return Err(ModelError::InvalidParameter {
                parameter: "n_parameters".to_string(),
                value: params.n_parameters() as f64,
            });
        }

        // Validate that CL and V are positive (after exp transformation)
        let cl = params.fixed_effects[0].exp();
        let v = params.fixed_effects[1].exp();

        if cl <= 0.0 {
            return Err(ModelError::InvalidParameter {
                parameter: "CL".to_string(),
                value: cl,
            });
        }

        if v <= 0.0 {
            return Err(ModelError::InvalidParameter {
                parameter: "V".to_string(),
                value: v,
            });
        }

        if params.residual_variance <= 0.0 {
            return Err(ModelError::InvalidParameter {
                parameter: "residual_variance".to_string(),
                value: params.residual_variance,
            });
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_one_compartment_model() {
        let model = OneCompartmentModel::new();
        assert_eq!(model.n_compartments(), 1);
        assert_eq!(model.parameter_names(), vec!["CL", "V"]);
        
        let params = model.default_parameters();
        assert!(model.validate_parameters(&params).is_ok());
    }

    #[test]
    fn test_derivatives() {
        let model = OneCompartmentModel::new();
        let params = model.default_parameters();
        let mut state = ModelState::new(1);
        state.compartments[0] = 100.0; // 100 mg
        
        let derivatives = model.derivatives(&state, &params);
        assert_eq!(derivatives.len(), 1);
        assert!(derivatives[0] < 0.0); // Should be decreasing
    }

    #[test]
    fn test_observation_function() {
        let model = OneCompartmentModel::new();
        let params = model.default_parameters();
        let mut state = ModelState::new(1);
        state.compartments[0] = 100.0; // 100 mg
        
        let conc = model.observation_function(&state, &params, 1);
        let v = params.fixed_effects[1].exp();
        assert!((conc - 100.0 / v).abs() < 1e-10);
    }
}