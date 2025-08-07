use super::compartment::{CompartmentModelTrait, ModelParameters, ModelState};
use super::ModelError;
use nalgebra::DVector;

pub struct TwoCompartmentModel {
    // Model: 
    // dA1/dt = -(CL/V1 + Q/V1) * A1 + Q/V2 * A2
    // dA2/dt = Q/V1 * A1 - Q/V2 * A2
    // Where A1, A2 are amounts in central and peripheral compartments
    // CL is clearance, V1 is central volume, V2 is peripheral volume, Q is intercompartmental clearance
}

impl TwoCompartmentModel {
    pub fn new() -> Self {
        Self {}
    }
}

impl CompartmentModelTrait for TwoCompartmentModel {
    fn n_compartments(&self) -> usize {
        2
    }

    fn parameter_names(&self) -> Vec<String> {
        vec![
            "CL".to_string(),
            "V1".to_string(),
            "Q".to_string(),
            "V2".to_string(),
        ]
    }

    fn default_parameters(&self) -> ModelParameters {
        let param_names = self.parameter_names();
        let mut params = ModelParameters::new(4, param_names);
        
        // Typical values for a two-compartment model
        params.fixed_effects[0] = 1.0_f64.ln();  // ln(CL) = ln(1.0 L/h) = 0.0
        params.fixed_effects[1] = 3.0_f64.ln();  // ln(V1) = ln(20 L) ≈ 2.996
        params.fixed_effects[2] = 0.5_f64.ln();  // ln(Q) = ln(0.5 L/h) ≈ -0.693
        params.fixed_effects[3] = 4.0_f64.ln();  // ln(V2) = ln(50 L) ≈ 3.912
        
        // Inter-individual variability
        for i in 0..4 {
            params.random_effects_variance[i][i] = 0.09; // 30% CV
        }
        
        // Residual error
        params.residual_variance = 0.01; // 10% CV
        
        params
    }

    fn derivatives(&self, state: &ModelState, params: &ModelParameters) -> DVector<f64> {
        let cl = params.fixed_effects[0].exp();
        let v1 = params.fixed_effects[1].exp();
        let q = params.fixed_effects[2].exp();
        let v2 = params.fixed_effects[3].exp();
        
        let a1 = state.compartments[0];
        let a2 = state.compartments[1];
        
        let mut derivatives = DVector::<f64>::zeros(2);
        
        // dA1/dt = -(CL/V1 + Q/V1) * A1 + Q/V2 * A2
        derivatives[0] = -(cl / v1 + q / v1) * a1 + (q / v2) * a2;
        
        // dA2/dt = Q/V1 * A1 - Q/V2 * A2
        derivatives[1] = (q / v1) * a1 - (q / v2) * a2;
        
        derivatives
    }

    fn observation_function(&self, state: &ModelState, params: &ModelParameters, compartment: usize) -> f64 {
        match compartment {
            1 => {
                // Central compartment concentration
                let v1 = params.fixed_effects[1].exp();
                state.compartments[0] / v1
            }
            2 => {
                // Peripheral compartment concentration
                let v2 = params.fixed_effects[3].exp();
                state.compartments[1] / v2
            }
            _ => 0.0,
        }
    }

    fn validate_parameters(&self, params: &ModelParameters) -> Result<(), ModelError> {
        if params.n_parameters() != 4 {
            return Err(ModelError::InvalidParameter {
                parameter: "n_parameters".to_string(),
                value: params.n_parameters() as f64,
            });
        }

        // Validate that all parameters are positive after exp transformation
        let param_values = vec![
            ("CL", params.fixed_effects[0].exp()),
            ("V1", params.fixed_effects[1].exp()),
            ("Q", params.fixed_effects[2].exp()),
            ("V2", params.fixed_effects[3].exp()),
        ];

        for (name, value) in param_values {
            if value <= 0.0 {
                return Err(ModelError::InvalidParameter {
                    parameter: name.to_string(),
                    value,
                });
            }
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
    fn test_two_compartment_model() {
        let model = TwoCompartmentModel::new();
        assert_eq!(model.n_compartments(), 2);
        assert_eq!(model.parameter_names().len(), 4);
        
        let params = model.default_parameters();
        assert!(model.validate_parameters(&params).is_ok());
    }

    #[test]
    fn test_two_compartment_derivatives() {
        let model = TwoCompartmentModel::new();
        let params = model.default_parameters();
        let mut state = ModelState::new(2);
        state.compartments[0] = 100.0; // Central compartment
        state.compartments[1] = 0.0;   // Peripheral compartment
        
        let derivatives = model.derivatives(&state, &params);
        assert_eq!(derivatives.len(), 2);
        assert!(derivatives[0] < 0.0); // Central compartment decreasing
        assert!(derivatives[1] > 0.0); // Peripheral compartment increasing
    }
}