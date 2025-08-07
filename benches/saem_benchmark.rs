use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nmodes::models::{CompartmentModel, ModelType, ModelState};
use nmodes::solver::{RungeKuttaSolver, OdeSolver, SolverConfig};
use nalgebra::DVector;

struct BenchmarkSystem {
    model: CompartmentModel,
    params: nmodes::models::ModelParameters,
}

impl nmodes::solver::OdeSystem for BenchmarkSystem {
    fn derivatives(&self, t: f64, y: &DVector<f64>) -> DVector<f64> {
        let state = ModelState {
            compartments: y.clone(),
            time: t,
        };
        self.model.derivatives(&state, &self.params)
    }

    fn dimension(&self) -> usize {
        self.model.n_compartments()
    }
}

fn benchmark_ode_solving(c: &mut Criterion) {
    let model = CompartmentModel::new(ModelType::TwoCompartment).unwrap();
    let params = model.default_parameters();
    let system = BenchmarkSystem { model, params };
    
    let solver = RungeKuttaSolver::new();
    let config = SolverConfig::default();
    let y0 = DVector::from_vec(vec![100.0, 0.0]);
    
    c.bench_function("ode_solve_two_compartment", |b| {
        b.iter(|| {
            solver.solve(
                black_box(&system),
                black_box((0.0, 24.0)),
                black_box(&y0),
                black_box(&config),
            ).unwrap()
        })
    });
}

fn benchmark_model_derivatives(c: &mut Criterion) {
    let model = CompartmentModel::new(ModelType::ThreeCompartment).unwrap();
    let params = model.default_parameters();
    let state = ModelState {
        compartments: DVector::from_vec(vec![100.0, 50.0, 25.0]),
        time: 1.0,
    };
    
    c.bench_function("model_derivatives_three_compartment", |b| {
        b.iter(|| {
            model.derivatives(
                black_box(&state),
                black_box(&params),
            )
        })
    });
}

criterion_group!(benches, benchmark_ode_solving, benchmark_model_derivatives);
criterion_main!(benches);