use super::compartment::{CompartmentModelTrait, ModelParameters, ModelState};
use super::ModelError;
use nalgebra::DVector;

pub struct ThreeCompartmentModel {
    // Model:
    // dA1/dt = -(CL/V1 + Q2/V1 + Q3/V1) * A1 + Q2/V2 * A2 + Q3/V3 * A3
    // dA2/dt = Q2/V1 * A1 - Q2/V2 * A2
    // dA3/dt = Q3/V1 * A1 - Q3/V3 * A3
}

impl ThreeCompartmentModel {
    pub fn new() -> Self {
        Self {}
    }
}

impl CompartmentModelTrait for ThreeCompartmentModel {
    fn n_compartments(&self) -> usize {
        3
    }

    fn parameter_names(&self) -> Vec<String> {
        vec![
            "CL".to_string(),
            "V1".to_string(),
            "Q2".to_string(),
            "V2".to_string(),
            "Q3".to_string(),
            "V3".to_string(),
        ]
    }

    fn default_parameters(&self) -> ModelParameters {
        let param_names = self.parameter_names();
        let mut params = ModelParameters::new(6, param_names);
        
        // Typical values for a three-compartment model
        params.fixed_effects[0] = 1.0_f64.ln();   // ln(CL) = ln(1.0 L/h) = 0.0
        params.fixed_effects[1] = 3.0_f64.ln();   // ln(V1) = ln(20 L) ≈ 2.996
        params.fixed_effects[2] = 0.5_f64.ln();   // ln(Q2) = ln(0.5 L/h) ≈ -0.693
        params.fixed_effects[3] = 4.0_f64.ln();   // ln(V2) = ln(50 L) ≈ 3.912
        params.fixed_effects[4] = 0.2_f64.ln();   // ln(Q3) = ln(0.2 L/h) ≈ -1.609
        params.fixed_effects[5] = 5.0_f64.ln();   // ln(V3) = ln(150 L) ≈ 5.011
        
        // Inter-individual variability
        for i in 0..6 {
            params.random_effects_variance[i][i] = 0.09; // 30% CV
        }
        
        // Residual error
        params.residual_variance = 0.01; // 10% CV
        
        params
    }

    fn derivatives(&self, state: &ModelState, params: &ModelParameters) -> DVector<f64> {
        let cl = params.fixed_effects[0].exp();
        let v1 = params.fixed_effects[1].exp();
        let q2 = params.fixed_effects[2].exp();
        let v2 = params.fixed_effects[3].exp();
        let q3 = params.fixed_effects[4].exp();
        let v3 = params.fixed_effects[5].exp();
        
        let a1 = state.compartments[0];
        let a2 = state.compartments[1];
        let a3 = state.compartments[2];
        
        let mut derivatives = DVector::<f64>::zeros(3);
        
        // dA1/dt = -(CL/V1 + Q2/V1 + Q3/V1) * A1 + Q2/V2 * A2 + Q3/V3 * A3
        derivatives[0] = -(cl / v1 + q2 / v1 + q3 / v1) * a1 + (q2 / v2) * a2 + (q3 / v3) * a3;
        
        // dA2/dt = Q2/V1 * A1 - Q2/V2 * A2
        derivatives[1] = (q2 / v1) * a1 - (q2 / v2) * a2;
        
        // dA3/dt = Q3/V1 * A1 - Q3/V3 * A3
        derivatives[2] = (q3 / v1) * a1 - (q3 / v3) * a3;
        
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
                // First peripheral compartment concentration
                let v2 = params.fixed_effects[3].exp();
                state.compartments[1] / v2
            }
            3 => {
                // Second peripheral compartment concentration
                let v3 = params.fixed_effects[5].exp();
                state.compartments[2] / v3
            }
            _ => 0.0,
        }
    }

    fn validate_parameters(&self, params: &ModelParameters) -> Result<(), ModelError> {
        if params.n_parameters() != 6 {
            return Err(ModelError::InvalidParameter {
                parameter: "n_parameters".to_string(),
                value: params.n_parameters() as f64,
            });
        }

        // Validate that all parameters are positive after exp transformation
        let param_values = vec![
            ("CL", params.fixed_effects[0].exp()),
            ("V1", params.fixed_effects[1].exp()),
            ("Q2", params.fixed_effects[2].exp()),
            ("V2", params.fixed_effects[3].exp()),
            ("Q3", params.fixed_effects[4].exp()),
            ("V3", params.fixed_effects[5].exp()),
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
    fn test_three_compartment_model() {
        let model = ThreeCompartmentModel::new();
        assert_eq!(model.n_compartments(), 3);
        assert_eq!(model.parameter_names().len(), 6);
        
        let params = model.default_parameters();
        assert!(model.validate_parameters(&params).is_ok());
    }

    #[test]
    fn test_three_compartment_derivatives() {
        let model = ThreeCompartmentModel::new();
        let params = model.default_parameters();
        let mut state = ModelState::new(3);
        state.compartments[0] = 100.0; // Central compartment
        state.compartments[1] = 0.0;   // First peripheral
        state.compartments[2] = 0.0;   // Second peripheral
        
        let derivatives = model.derivatives(&state, &params);
        assert_eq!(derivatives.len(), 3);
        assert!(derivatives[0] < 0.0); // Central compartment decreasing
        assert!(derivatives[1] > 0.0); // First peripheral increasing
        assert!(derivatives[2] > 0.0); // Second peripheral increasing
    }
}