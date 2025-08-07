# NMODES - Nonlinear Mixed Effects Differential Equation Solver

A comprehensive Rust implementation of population pharmacokinetics modeling using multiple estimation methods including Stochastic Approximation Expectation Maximization (SAEM) and First Order Conditional Estimation (FOCE). This library can fit NONMEM-style datasets to 1, 2, or 3 compartment pharmacokinetic models with robust parameter estimation and comprehensive diagnostics.

## Features

- **Multiple Compartment Models**: Support for 1, 2, and 3 compartment pharmacokinetic models
- **Multiple Estimation Methods**: 
  - **SAEM**: Stochastic Approximation Expectation Maximization for robust parameter estimation
  - **FOCE**: First Order Conditional Estimation for fast, deterministic estimation
  - **FOCE-I**: FOCE with interaction for improved accuracy with non-linear models
- **NONMEM Compatibility**: Reads standard NONMEM dataset formats and produces similar output
- **Adaptive ODE Solving**: High-performance numerical integration with error control
- **Comprehensive Diagnostics**: Goodness-of-fit metrics, residual analysis, and convergence assessment
- **Parameter Statistics**: %RSE (Relative Standard Error) and shrinkage calculations
- **Parallel Processing**: Multi-threaded MCMC sampling for improved performance
- **Production Ready**: Extensive error handling, validation, and logging

## Installation

### Prerequisites
- Rust 1.70 or later
- Git

### Build from Source

```bash
git clone <repository-url>
cd nmodes
cargo build --release
```

The compiled binary will be available at `./target/release/nmodes`.

## Quick Start

### Basic One-Compartment Model

```bash
# Fit a one-compartment model using SAEM (default)
./target/release/nmodes -d examples/example_dataset.csv -m 1comp -o results/

# Fit using FOCE method for faster estimation
./target/release/nmodes -d examples/example_dataset.csv -m 1comp -e foce -o results_foce/
```

### Multiple Model Comparison

```bash
# Compare all compartment models with SAEM
./target/release/nmodes -d examples/example_dataset.csv -m all -e saem -o model_comparison/

# Compare specific models with default method
./target/release/nmodes -d examples/example_dataset.csv -m 1comp -m 2comp -o two_model_comparison/

# Test all models with FOCE for quick assessment
./target/release/nmodes -d examples/example_dataset.csv -m all -e foce -o quick_model_screen/
```

### Multiple Method Comparison

```bash
# Compare all estimation methods on 2-compartment model
./target/release/nmodes -d examples/example_dataset.csv -m 2comp -e all -o method_comparison/

# Compare SAEM vs FOCE on specific model
./target/release/nmodes -d examples/example_dataset.csv -m 1comp -e saem -e foce -o saem_vs_foce/

# Test FOCE variants for method selection
./target/release/nmodes -d examples/example_dataset.csv -m 2comp -e foce -e foce-i -o foce_comparison/
```

### Comprehensive Analysis

```bash
# Full comparison: all models and all methods
./target/release/nmodes -d examples/example_dataset.csv -m all -e all --compare -o full_analysis/

# Production workflow: screen with FOCE, validate with SAEM
./target/release/nmodes -d examples/example_dataset.csv -m all -e foce -o screening/
./target/release/nmodes -d examples/example_dataset.csv -m 2comp -e saem -o final_model/
```

### Advanced Usage

```bash
# Single model with custom SAEM settings
./target/release/nmodes \
  --dataset data/my_study.csv \
  --model 2comp \
  --method saem \
  --iterations 2000 \
  --burn-in 400 \
  --chains 6 \
  --output results/two_comp_analysis/

# Multiple models with custom settings
./target/release/nmodes \
  --dataset data/my_study.csv \
  --model 1comp --model 2comp --model 3comp \
  --method saem \
  --iterations 2000 \
  --burn-in 400 \
  --compare \
  --output results/comprehensive_analysis/
```

## Command Line Interface

### Required Arguments

- `-d, --dataset <FILE>`: Path to NONMEM-style dataset CSV file

### Optional Arguments

- `-m, --model <TYPE>`: Compartment model type
  - `1comp`: One-compartment model (default)
  - `2comp`: Two-compartment model  
  - `3comp`: Three-compartment model
  - `all`: All compartment models (1comp, 2comp, 3comp)
  - **Multiple models**: Use multiple `-m` flags (e.g., `-m 1comp -m 2comp`)
- `-e, --method <METHOD>`: Estimation method
  - `saem`: Stochastic Approximation EM (default)
  - `foce`: First Order Conditional Estimation
  - `foce-i`: FOCE with interaction
  - `all`: All estimation methods (saem, foce, foce-i)
  - **Multiple methods**: Use multiple `-e` flags (e.g., `-e saem -e foce`)
- `-o, --output <DIR>`: Output directory for results (default: `./output`)
- `-i, --iterations <N>`: Number of SAEM iterations (default: 1000)
- `-b, --burn-in <N>`: Number of burn-in iterations (default: 200)
- `-c, --chains <N>`: Number of MCMC chains (default: 4)
- `--compare`: Force generation of comparison reports (automatic when multiple analyses run)

### Single Analysis Examples

```bash
# Quick analysis with defaults (SAEM)
./target/release/nmodes -d data.csv

# Fast FOCE analysis
./target/release/nmodes -d data.csv -e foce -i 100 -o foce_results/

# Production SAEM analysis with more iterations
./target/release/nmodes -d data.csv -m 2comp -e saem -i 3000 -b 600 -o saem_results/
```

### Multiple Model Examples

```bash
# Compare all compartment models with SAEM
./target/release/nmodes -d data.csv -m all -e saem -o model_comparison/

# Compare 1 and 2 compartment models with FOCE
./target/release/nmodes -d data.csv -m 1comp -m 2comp -e foce -o simple_models/

# Test all models with quick FOCE screening
./target/release/nmodes -d data.csv -m all -e foce -i 50 -o quick_screen/
```

### Multiple Method Examples

```bash
# Compare all methods on 2-compartment model
./target/release/nmodes -d data.csv -m 2comp -e all -o method_comparison/

# Compare SAEM vs FOCE methods
./target/release/nmodes -d data.csv -m 1comp -e saem -e foce -o saem_vs_foce/

# Test FOCE variants
./target/release/nmodes -d data.csv -m 2comp -e foce -e foce-i -o foce_variants/
```

### Comprehensive Analysis Examples

```bash
# Full model and method comparison
./target/release/nmodes -d data.csv -m all -e all --compare -o comprehensive/

# Production workflow: screen all models, then detailed analysis
./target/release/nmodes -d data.csv -m all -e foce -i 50 -o screening/
# Then run detailed analysis on best model from screening results
./target/release/nmodes -d data.csv -m 2comp -e saem -i 2000 -o final_analysis/
```

### Legacy Examples (still supported)

```bash

# FOCE-I for complex three-compartment model
./target/release/nmodes -d complex_data.csv -m 3comp -e foce-i -i 200 -o foce_i_results/

# Manual method comparison workflow (now automated with -e all)
./target/release/nmodes -d data.csv -m 2comp -e foce -o manual_foce/
./target/release/nmodes -d data.csv -m 2comp -e saem -o manual_saem/
```

## Output Structure

When running multiple models or methods, the output directory is organized as follows:

```
output_directory/
├── one-compartment_SAEM/          # Individual analysis results
│   ├── parameter_estimates.json
│   ├── predictions.csv
│   ├── diagnostics.json
│   └── summary_report.txt
├── one-compartment_FOCE/
│   ├── foce_results.json
│   ├── foce_predictions.csv
│   └── foce_summary_report.txt
├── two-compartment_SAEM/
├── two-compartment_FOCE/
├── model_comparison_report.txt    # Comprehensive comparison report
└── model_comparison.csv           # Machine-readable comparison data
```

### Comparison Report Contents

**Summary Table**: All analyses with key metrics (OFV, AIC, BIC, convergence status)

**Model Ranking**: Ordered by AIC with delta-AIC values for model selection

**Parameter Comparison**: Side-by-side parameter estimates for converged models

**Method Performance**: Convergence rates and computational efficiency by method

**Recommendations**: Automated suggestions based on statistical criteria

## Dataset Format

The program expects NONMEM-style CSV files with specific column names. All column names are case-sensitive.

### Required Columns

| Column | Description | Type | Example |
|--------|-------------|------|---------|
| `ID` | Individual identifier | Integer | 1, 2, 3, ... |
| `TIME` | Time of observation/dose (hours) | Float | 0.0, 0.5, 1.0, ... |
| `DV` | Dependent variable (concentration) | Float | 8.5, 7.2, 5.1, ... |
| `AMT` | Dose amount (mg) | Float | 100.0, 150.0, ... |
| `EVID` | Event ID | Integer | 0 (observation), 1 (dose) |

### Optional Columns

| Column | Description | Default |
|--------|-------------|---------|
| `CMT` | Compartment number | 1 |
| `RATE` | Infusion rate (mg/h) | 0 (bolus) |
| `II` | Interdose interval (h) | - |
| `ADDL` | Additional doses | 0 |
| `SS` | Steady state flag | 0 |

### Example Datasets

The repository includes three comprehensive example datasets with realistic demographics and PK profiles:

#### 1. One-Compartment Dataset (`examples/one_compartment_dataset.csv`)
- **Study Design**: 10 subjects, 100 mg oral doses every 4 weeks (5 doses total)
- **Sampling**: Weekly observations for 40 weeks (280 time points per subject)
- **PK Profile**: Classic one-compartment elimination pattern
- **Demographics**: Age 29-61, Weight 59.8-91.3 kg, Mixed gender/race, CRCL 65.2-125.8 mL/min

#### 2. Two-Compartment Dataset (`examples/two_compartment_dataset.csv`)
- **Study Design**: 10 subjects, 45-70 mg oral doses every 4 weeks (dose varies by subject)
- **Sampling**: Weekly observations for 40 weeks
- **PK Profile**: Two-compartment distribution with slower terminal elimination
- **Demographics**: Matched to one-compartment study for comparison

#### 3. Three-Compartment Dataset (`examples/three_compartment_dataset.csv`)
- **Study Design**: 10 subjects, 60 mg oral doses every 4 weeks
- **Sampling**: Weekly observations for 40 weeks
- **PK Profile**: Three-compartment distribution with complex elimination phases
- **Demographics**: Same population as other studies

**Demographic Variables:**
- `WEIGHT`: Body weight (kg) - Normal range 60-90 kg
- `AGE`: Age (years) - Range 29-61 years
- `SEX`: Gender (1=Male, 0=Female)
- `RACE`: Race (1=Caucasian, 2=African American, 3=Hispanic)
- `CRCL`: Creatinine clearance (mL/min) - Range 65-126 mL/min
- `BMI`: Body mass index (kg/m²) - Range 22.9-29.1 kg/m²

**Usage Examples:**
```bash
# Analyze one-compartment data with SAEM
./target/release/nmodes -d examples/one_compartment_dataset.csv -m 1comp -o results_1comp/

# Analyze two-compartment data with FOCE
./target/release/nmodes -d examples/two_compartment_dataset.csv -m 2comp -e foce -o results_2comp_foce/

# Analyze two-compartment data with SAEM
./target/release/nmodes -d examples/two_compartment_dataset.csv -m 2comp -o results_2comp/

# Analyze three-compartment data with FOCE-I
./target/release/nmodes -d examples/three_compartment_dataset.csv -m 3comp -e foce-i -o results_3comp/
```

### Dataset Format Example

Here's a sample of the dataset format:

```csv
ID,TIME,DV,AMT,EVID,CMT,RATE,II,ADDL,WEIGHT,AGE,SEX,RACE,CRCL,BMI
1,0.0,,100.0,1,1,0.0,672.0,4,75.2,45,1,1,95.3,24.8
1,1.0,8.45,,0,1,,,,,,,,,
1,7.0,6.12,,0,1,,,,,,,,,
1,14.0,4.23,,0,1,,,,,,,,,
2,0.0,,100.0,1,1,0.0,672.0,4,68.5,52,0,2,78.4,26.1
2,1.0,9.12,,0,1,,,,,,,,,
```

**Key Points:**
- **Dosing Records**: `EVID=1` with `AMT` specifying dose amount
- **Observation Records**: `EVID=0` with `DV` specifying concentration
- **Multiple Dosing**: `II` (interdose interval) and `ADDL` (additional doses) for repeated dosing
- **Demographics**: Additional columns for covariate analysis

### Data Validation Rules

- **Time values**: Must be non-negative and in ascending order within each individual
- **Dose amounts**: Must be positive for dosing events (EVID=1)
- **Concentrations**: Should be non-negative; negative values generate warnings
- **Individual IDs**: Must be unique integers
- **Event ordering**: Doses typically precede observations

## Pharmacokinetic Models

### One-Compartment Model

**Differential Equation:**
```
dA/dt = -(CL/V) × A
```

**Parameters:**
- `CL`: Clearance (L/h)
- `V`: Volume of distribution (L)

**Typical Values:**
- CL: 1.0 L/h
- V: 20 L

### Two-Compartment Model

**Differential Equations:**
```
dA1/dt = -(CL/V1 + Q/V1) × A1 + (Q/V2) × A2
dA2/dt = (Q/V1) × A1 - (Q/V2) × A2
```

**Parameters:**
- `CL`: Clearance from central compartment (L/h)
- `V1`: Central volume of distribution (L)
- `Q`: Intercompartmental clearance (L/h)
- `V2`: Peripheral volume of distribution (L)

**Typical Values:**
- CL: 1.0 L/h
- V1: 20 L
- Q: 0.5 L/h
- V2: 50 L

### Three-Compartment Model

**Differential Equations:**
```
dA1/dt = -(CL/V1 + Q2/V1 + Q3/V1) × A1 + (Q2/V2) × A2 + (Q3/V3) × A3
dA2/dt = (Q2/V1) × A1 - (Q2/V2) × A2
dA3/dt = (Q3/V1) × A1 - (Q3/V3) × A3
```

**Parameters:**
- `CL`: Clearance from central compartment (L/h)
- `V1`: Central volume of distribution (L)
- `Q2`: First intercompartmental clearance (L/h)
- `V2`: First peripheral volume (L)
- `Q3`: Second intercompartmental clearance (L/h)
- `V3`: Second peripheral volume (L)

## Estimation Methods

### SAEM (Stochastic Approximation Expectation Maximization)

**Advantages:**
- Robust convergence for complex models
- Handles non-linear mixed effects naturally
- Good performance with sparse data
- Provides uncertainty quantification

**Best For:**
- Complex pharmacokinetic models
- Small to medium datasets (< 500 individuals)
- When robustness is more important than speed
- Research and exploratory analysis

**Single Model Configuration:**
```bash
# Standard SAEM analysis
./target/release/nmodes -d data.csv -e saem -i 2000 -b 400 -c 6

# SAEM with example datasets
./target/release/nmodes -d examples/one_compartment_dataset.csv -m 1comp -e saem -i 1500 -b 300 -o saem_1comp/
./target/release/nmodes -d examples/two_compartment_dataset.csv -m 2comp -e saem -i 2000 -b 400 -o saem_2comp/
./target/release/nmodes -d examples/three_compartment_dataset.csv -m 3comp -e saem -i 2500 -b 500 -o saem_3comp/
```

**Multiple Model Configuration:**
```bash
# SAEM analysis across all models
./target/release/nmodes -d examples/example_dataset.csv -m all -e saem -i 1500 -b 300 -o saem_all_models/

# SAEM on specific models for comparison
./target/release/nmodes -d examples/example_dataset.csv -m 1comp -m 2comp -e saem -o saem_comparison/
```

### FOCE (First Order Conditional Estimation)

**Advantages:**
- Fast, deterministic estimation
- Industry standard method
- Excellent for large datasets
- Provides standard errors directly

**Best For:**
- Large datasets (> 500 individuals)
- Production pharmacokinetic analysis
- Regulatory submissions
- When speed is critical

**Single Model Configuration:**
```bash
# Fast FOCE analysis
./target/release/nmodes -d data.csv -e foce -i 100

# FOCE with example datasets
./target/release/nmodes -d examples/one_compartment_dataset.csv -m 1comp -e foce -i 50 -o foce_1comp/
./target/release/nmodes -d examples/two_compartment_dataset.csv -m 2comp -e foce -i 100 -o foce_2comp/
./target/release/nmodes -d examples/three_compartment_dataset.csv -m 3comp -e foce -i 150 -o foce_3comp/
```

**Multiple Model Configuration:**
```bash
# FOCE screening across all models
./target/release/nmodes -d examples/example_dataset.csv -m all -e foce -i 50 -o foce_screening/

# FOCE comparison on selected models
./target/release/nmodes -d examples/example_dataset.csv -m 1comp -m 3comp -e foce -o foce_selected/
```

### FOCE-I (FOCE with Interaction)

**Advantages:**
- More accurate than standard FOCE
- Better handling of non-linear models
- Accounts for parameter interactions
- Still faster than SAEM

**Best For:**
- Non-linear pharmacokinetic models
- Models with significant parameter correlations
- When FOCE shows bias but SAEM is too slow
- Complex covariate relationships

**Single Model Configuration:**
```bash
# FOCE with interaction
./target/release/nmodes -d data.csv -e foce-i -i 200

# FOCE-I with example datasets
./target/release/nmodes -d examples/one_compartment_dataset.csv -m 1comp -e foce-i -i 75 -o foce_i_1comp/
./target/release/nmodes -d examples/two_compartment_dataset.csv -m 2comp -e foce-i -i 125 -o foce_i_2comp/
./target/release/nmodes -d examples/three_compartment_dataset.csv -m 3comp -e foce-i -i 175 -o foce_i_3comp/
```

**Multiple Model Configuration:**
```bash
# FOCE-I analysis across all models
./target/release/nmodes -d examples/example_dataset.csv -m all -e foce-i -i 100 -o foce_i_all/

# FOCE-I on complex models only
./target/release/nmodes -d examples/example_dataset.csv -m 2comp -m 3comp -e foce-i -o foce_i_complex/
```

### Multiple Method Comparison

```bash
# Compare all methods on one model
./target/release/nmodes -d examples/example_dataset.csv -m 2comp -e all -o all_methods/

# Compare SAEM vs FOCE methods
./target/release/nmodes -d examples/example_dataset.csv -m 1comp -e saem -e foce -o saem_vs_foce/

# Compare FOCE variants
./target/release/nmodes -d examples/example_dataset.csv -m 2comp -e foce -e foce-i -o foce_variants/
```

### Method Comparison Guidelines

| Criterion | SAEM | FOCE | FOCE-I |
|-----------|------|------|--------|
| **Speed** | Slow | Fast | Medium |
| **Robustness** | High | Medium | High |
| **Large Datasets** | Poor | Excellent | Good |
| **Non-linear Models** | Excellent | Fair | Good |
| **Regulatory Acceptance** | Good | Excellent | Excellent |
| **Uncertainty Quantification** | Excellent | Good | Good |

**Recommended Workflows:**

#### 1. Model Selection Workflow
```bash
# Step 1: Quick screening of all models with FOCE
./target/release/nmodes -d data.csv -m all -e foce -i 50 -o screening/

# Step 2: Detailed analysis of top 2 models with SAEM
./target/release/nmodes -d data.csv -m 1comp -m 2comp -e saem -i 2000 -o detailed/

# Step 3: Final validation with FOCE-I on best model
./target/release/nmodes -d data.csv -m 2comp -e foce-i -i 150 -o final/
```

#### 2. Method Validation Workflow
```bash
# Step 1: Compare all methods on selected model
./target/release/nmodes -d data.csv -m 2comp -e all --compare -o method_comparison/

# Step 2: Extended SAEM run for final estimates
./target/release/nmodes -d data.csv -m 2comp -e saem -i 3000 -b 600 -o final_saem/
```

#### 3. Comprehensive Analysis Workflow
```bash
# Single command for complete analysis
./target/release/nmodes -d data.csv -m all -e all --compare -o comprehensive/

# Review comparison report for best model-method combination
cat comprehensive/model_comparison_report.txt
```

#### 4. Production Workflow
```bash
# Step 1: Automated model and method selection
./target/release/nmodes -d data.csv -m all -e foce -e foce-i -o selection/

# Step 2: Robust estimation with best model (based on selection results)
./target/release/nmodes -d data.csv -m 2comp -e saem -i 2000 -o production/

# Step 3: Cross-validation with alternative method
./target/release/nmodes -d data.csv -m 2comp -e foce-i -o validation/
```

### Complete Analysis Workflow Examples

#### Example 1: Automated Model Selection
```bash
# Screen all models with FOCE for speed
./target/release/nmodes -d examples/one_compartment_dataset.csv -m all -e foce -i 50 -o model_screen/

# Validate top models with SAEM
./target/release/nmodes -d examples/one_compartment_dataset.csv -m 1comp -m 2comp -e saem -i 1500 -b 300 -o validation/

# Final analysis with best model and method
./target/release/nmodes -d examples/one_compartment_dataset.csv -m 1comp -e saem -i 2000 -b 400 -o final/
```

#### Example 2: Method Comparison Study
```bash
# Single command to compare all methods
./target/release/nmodes -d examples/two_compartment_dataset.csv -m 2comp -e all --compare -o method_study/

# Review automated comparison report
cat method_study/model_comparison_report.txt

# Extract specific metrics for further analysis
grep "AIC" method_study/model_comparison.csv
```

#### Example 3: Comprehensive Analysis
```bash
# Complete analysis: all models and methods
./target/release/nmodes -d examples/three_compartment_dataset.csv -m all -e all --compare -o complete_analysis/

# Results automatically organized and compared
ls complete_analysis/
# Shows: one-compartment_SAEM/, one-compartment_FOCE/, one-compartment_FOCE-I/,
#        two-compartment_SAEM/, two-compartment_FOCE/, two-compartment_FOCE-I/,
#        three-compartment_SAEM/, three-compartment_FOCE/, three-compartment_FOCE-I/,
#        model_comparison_report.txt, model_comparison.csv

# View summary of all results
head -20 complete_analysis/model_comparison_report.txt
```

#### Example 4: Custom Comparison Workflow
```bash
# Custom selection: specific models and methods
./target/release/nmodes -d examples/example_dataset.csv \
  -m 1comp -m 2comp \
  -e saem -e foce-i \
  --iterations 1500 \
  --burn-in 300 \
  --compare \
  -o custom_comparison/

# Results in 4 analyses: 1comp+SAEM, 1comp+FOCE-I, 2comp+SAEM, 2comp+FOCE-I
# Plus automated comparison report

# Extract best model-method combination
grep "Best fitting model" custom_comparison/model_comparison_report.txt
```

### Console Output Examples

When running multiple analyses, the console shows progress and a summary table:

```
$ ./target/release/nmodes -d data.csv -m all -e saem -e foce -o comparison/

Starting NMODES analysis
Dataset: "data.csv"
Model types: [one-compartment, two-compartment, three-compartment]
Estimation methods: [SAEM, FOCE]
Output directory: "comparison/"

Running SAEM estimation with one-compartment model
Running FOCE estimation with one-compartment model
Running SAEM estimation with two-compartment model
Running FOCE estimation with two-compartment model
Running SAEM estimation with three-compartment model
Running FOCE estimation with three-compartment model

Analysis Summary:
Model           Method     OFV        LogLik     Converged  AIC    
--------------------------------------------------------------------------------
one-compartment SAEM       491.34     -245.67    true       495.3  
one-compartment FOCE       493.21     -246.61    true       497.2  
two-compartment SAEM       485.12     -242.56    true       493.1  
two-compartment FOCE       486.89     -243.45    true       494.9  
three-compartment SAEM     487.45     -243.73    true       499.5  
three-compartment FOCE     489.12     -244.56    false      501.1  

Best model by AIC: two-compartment with SAEM (AIC: 493.1)

Comparison report saved to: "comparison/model_comparison_report.txt"
Comparison CSV saved to: "comparison/model_comparison.csv"
```

## Output Files

The program generates comprehensive output files in the specified directory:

### SAEM Output Files

### 1. `parameter_estimates.json`
Complete parameter estimates with convergence information:
```json
{
  "fixed_effects": [0.693, 2.996],
  "random_effects_variance": [[0.09, 0.0], [0.0, 0.04]],
  "residual_variance": 0.01,
  "converged": true,
  "final_log_likelihood": -245.67,
  "parameter_statistics": [
    {
      "name": "CL",
      "estimate": 2.0,
      "rse_percent": 8.5
    }
  ]
}
```

### 2. `predictions.csv`
Original data with model predictions:
```csv
ID,TIME,DV,AMT,EVID,CMT,IPRED,PRED
1,0.0,,100.0,1,1,0.0,0.0
1,0.5,8.5,,0,1,8.2,8.0
1,1.0,7.2,,0,1,7.1,6.9
1,2.0,5.1,,0,1,5.0,4.8
```

**Column Descriptions:**
- `IPRED`: Individual predicted concentrations (using individual parameters)
- `PRED`: Population predicted concentrations (using population mean parameters)

### 3. `diagnostics.json`
Comprehensive model diagnostics:
```json
{
  "goodness_of_fit": {
    "aic": 495.34,
    "bic": 503.21,
    "log_likelihood": -245.67,
    "rmse": 0.85,
    "r_squared": 0.94
  },
  "convergence_diagnostics": {
    "converged": true,
    "parameter_stability": [0.008, 0.012]
  }
}
```

### 4. `parameter_trajectory.csv`
Parameter evolution during estimation:
```csv
iteration,CL,V,log_likelihood
0,1.000,20.000,-280.45
1,1.023,19.87,-275.23
...
999,2.001,18.95,-245.67
```

### FOCE Output Files

### 1. `foce_results.json`
Complete FOCE parameter estimates with covariance information:
```json
{
  "fixed_effects": [0.693, 2.996],
  "random_effects_variance": [[0.09, 0.0], [0.0, 0.04]],
  "residual_variance": 0.01,
  "converged": true,
  "objective_function_value": 491.34,
  "standard_errors": [0.058, 0.245],
  "covariance_matrix": [[0.0034, 0.0012], [0.0012, 0.0601]],
  "gradient_norm": 1.2e-7,
  "hessian_condition_number": 45.2
}
```

### 2. `foce_summary_report.txt`
FOCE-specific summary with standard errors:
```
PKPD FOCE Analysis Summary Report
=================================

Estimation Method: FOCE
Model Convergence: true
Total Iterations: 87
Objective Function Value: 491.340
Gradient Norm: 0.000000

Fixed Effects Parameter Estimates:
----------------------------------
Parameter  Estimate     SE        
---------  --------     --        
CL         2.001        0.058     
V          18.950       0.245     
```

### 3. `foce_predictions.csv`
FOCE predictions with individual and population estimates:
```csv
ID,TIME,DV,IPRED,PRED
1,0.0,,0.0,0.0
1,0.5,8.5,8.3,8.0
1,1.0,7.2,7.2,6.9
```

### 5. `summary_report.txt`
Human-readable summary with NONMEM-style formatting:
```
PKPD SAEM Analysis Summary Report
=================================

Model Convergence: true
Total Iterations: 1000
Final Log-Likelihood: -245.670
Objective Function Value: 491.340

Fixed Effects Parameter Estimates:
----------------------------------
Parameter  Estimate     %RSE      
---------  --------     ----      
CL         2.001        8.5       
V          18.950       12.3      

Random Effects Variance (Omega):
-------------------------------
Parameter       Estimate     Shrinkage%  
---------       --------     ----------  
CL(CL)          0.090        15.2        
V(V)            0.040        22.8        
```

## Programming Interface

### Basic Usage

```rust
use nmodes::*;

// Load dataset
let dataset = Dataset::from_csv("data.csv")?;

// Create two-compartment model
let model = CompartmentModel::new(ModelType::TwoCompartment)?;

// Configure SAEM estimation
let config = EstimationConfig::default()
    .with_method(EstimationMethod::Saem)
    .with_iterations(2000)
    .with_burnin(400)
    .with_chains(6);

// Run estimation
let mut estimator = SaemEstimator::new(model, config);
let results = estimator.fit(&dataset)?;

// Generate diagnostics
let diagnostics = diagnostics::generate_diagnostics(&dataset, &results)?;

// Save all results
output::save_results(
    Path::new("./results"),
    &results,
    &diagnostics,
    &dataset,
    estimator.model(),
)?;
```

### FOCE Estimation

```rust
use nmodes::*;

// Configure FOCE estimation
let config = EstimationConfig::default()
    .with_method(EstimationMethod::Foce)
    .with_foce_iterations(100)
    .with_foce_tolerance(1e-6);

// Run FOCE estimation
let mut estimator = FoceEstimator::new(model, config);
let results = estimator.fit(&dataset)?;

// FOCE provides standard errors directly
println!("Parameter estimates with standard errors:");
for (i, name) in results.parameter_names.iter().enumerate() {
    println!("{}: {:.3} ± {:.3}", 
             name, 
             results.fixed_effects[i], 
             results.standard_errors[i]);
}
```

### Complete Programming Examples

#### Example 1: Analyze Example Datasets
```rust
use nmodes::*;
use std::path::Path;

// Analyze one-compartment dataset
let dataset = Dataset::from_csv("examples/one_compartment_dataset.csv")?;

// Compare multiple models programmatically
let model_types = vec![
    ModelType::OneCompartment,
    ModelType::TwoCompartment,
    ModelType::ThreeCompartment,
];

let estimation_methods = vec![
    EstimationMethod::Foce,
    EstimationMethod::Saem,
];

let mut all_results = Vec::new();

for model_type in model_types {
    for estimation_method in &estimation_methods {
        let model = CompartmentModel::new(model_type.clone())?;
        
        let config = EstimationConfig::default()
            .with_method(estimation_method.clone())
            .with_iterations(match estimation_method {
                EstimationMethod::Saem => 1000,
                EstimationMethod::Foce => 50,
                EstimationMethod::FoceI => 75,
            });
        
        let (aic, converged) = match estimation_method {
            EstimationMethod::Saem => {
                let mut estimator = SaemEstimator::new(model, config);
                let results = estimator.fit(&dataset)?;
                let aic = -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64;
                (aic, results.converged)
            }
            EstimationMethod::Foce | EstimationMethod::FoceI => {
                let mut estimator = FoceEstimator::new(model, config);
                let results = estimator.fit(&dataset)?;
                let aic = -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64;
                (aic, results.converged)
            }
        };
        
        all_results.push((model_type.clone(), estimation_method.clone(), aic, converged));
    }
}

// Find best model-method combination
let best_result = all_results.iter()
    .filter(|(_, _, _, converged)| *converged)
    .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

if let Some((model_type, method, aic, _)) = best_result {
    println!("Best model: {:?} with {:?} (AIC: {:.2})", model_type, method, aic);
}
```

#### Example 2: Method Comparison Study

```rust
use nmodes::*;

// Load two-compartment dataset
let dataset = Dataset::from_csv("examples/two_compartment_dataset.csv")?;

// Compare all methods on two-compartment model
let methods = vec![
    EstimationMethod::Foce,
    EstimationMethod::FoceI,
    EstimationMethod::Saem,
];

let mut results_comparison = Vec::new();

for method in methods {
    let model = CompartmentModel::new(ModelType::TwoCompartment)?;
    let config = EstimationConfig::default()
        .with_method(method.clone())
        .with_iterations(match method {
            EstimationMethod::Saem => 2000,
            EstimationMethod::Foce => 100,
            EstimationMethod::FoceI => 125,
        })
        .with_burnin(if matches!(method, EstimationMethod::Saem) { 400 } else { 0 });
    
    let start_time = std::time::Instant::now();
    
    match method {
        EstimationMethod::Saem => {
            let mut estimator = SaemEstimator::new(model.clone(), config);
            let results = estimator.fit(&dataset)?;
            let duration = start_time.elapsed();
            
            results_comparison.push((
                method.to_string(),
                results.objective_function_value,
                results.converged,
                duration,
                results.fixed_effects.clone(),
            ));
        }
        EstimationMethod::Foce | EstimationMethod::FoceI => {
            let mut estimator = FoceEstimator::new(model.clone(), config);
            let results = estimator.fit(&dataset)?;
            let duration = start_time.elapsed();
            
            results_comparison.push((
                method.to_string(),
                results.objective_function_value,
                results.converged,
                duration,
                results.fixed_effects.clone(),
            ));
        }
    }
}

// Print comparison table
println!("\nMethod Comparison Results:");
println!("{:<10} {:<12} {:<10} {:<10} {:<15} {:<15}", 
         "Method", "OFV", "Converged", "Time(s)", "CL", "V1");
println!("{}", "-".repeat(80));

for (method, ofv, converged, duration, params) in results_comparison {
    println!("{:<10} {:<12.2} {:<10} {:<10.1} {:<15.3} {:<15.3}", 
             method, ofv, converged, duration.as_secs_f64(), 
             params[0].exp(), params[1].exp());
}
```

#### Example 3: Automated Model and Method Selection
```rust
use nmodes::*;
use std::path::Path;

fn run_automated_analysis(dataset_path: &str, output_base: &str) -> Result<()> {
    // Load and validate dataset
    let dataset = Dataset::from_csv(dataset_path)?;
    validation::validate_dataset(&dataset)?;
    
    println!("Dataset: {} individuals, {} observations", 
             dataset.n_individuals(), dataset.n_observations());
    
    // Test all combinations of models and methods
    let models = vec![
        (ModelType::OneCompartment, "1comp"),
        (ModelType::TwoCompartment, "2comp"),
        (ModelType::ThreeCompartment, "3comp"),
    ];
    
    let methods = vec![
        (EstimationMethod::Foce, "foce"),
        (EstimationMethod::FoceI, "foce-i"),
        (EstimationMethod::Saem, "saem"),
    ];
    
    let mut all_results = Vec::new();
    let mut best_aic = f64::INFINITY;
    let mut best_combination = None;
    
    for (model_type, model_name) in &models {
        for (estimation_method, method_name) in &methods {
            println!("Running {} with {} model", method_name, model_name);
            
            let model = CompartmentModel::new(model_type.clone())?;
            let config = EstimationConfig::default()
                .with_method(estimation_method.clone())
                .with_iterations(match estimation_method {
                    EstimationMethod::Saem => 1000,
                    EstimationMethod::Foce => 50,
                    EstimationMethod::FoceI => 75,
                });
            
            let (aic, converged, ofv) = match estimation_method {
                EstimationMethod::Saem => {
                    let mut estimator = SaemEstimator::new(model, config);
                    let results = estimator.fit(&dataset)?;
                    let aic = -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64;
                    (aic, results.converged, results.objective_function_value)
                }
                EstimationMethod::Foce | EstimationMethod::FoceI => {
                    let mut estimator = FoceEstimator::new(model, config);
                    let results = estimator.fit(&dataset)?;
                    let aic = -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64;
                    (aic, results.converged, results.objective_function_value)
                }
            };
            
            all_results.push((model_type.clone(), estimation_method.clone(), aic, converged, ofv));
            
            if aic < best_aic && converged {
                best_aic = aic;
                best_combination = Some((model_type.clone(), estimation_method.clone()));
            }
        }
    }
    
    // Print comparison table
    println!("\nComparison Results:");
    println!("{:<15} {:<10} {:<10} {:<10} {:<10}", "Model", "Method", "AIC", "OFV", "Converged");
    println!("{}", "-".repeat(70));
    
    for (model_type, method, aic, converged, ofv) in &all_results {
        println!("{:<15} {:<10} {:<10.1} {:<10.1} {:<10}", 
                 format!("{}", model_type), format!("{}", method), aic, ofv, converged);
    }
    
    // Run detailed analysis with best combination
    if let Some((best_model_type, best_method)) = best_combination {
        println!("\nRunning detailed analysis with best combination: {} + {}", 
                 best_model_type, best_method);
        
        let model = CompartmentModel::new(best_model_type.clone())?;
        let config = EstimationConfig::default()
            .with_method(best_method.clone())
            .with_iterations(match best_method {
                EstimationMethod::Saem => 2000,
                EstimationMethod::Foce => 200,
                EstimationMethod::FoceI => 250,
            })
            .with_burnin(if matches!(best_method, EstimationMethod::Saem) { 400 } else { 0 });
        
        match best_method {
            EstimationMethod::Saem => {
                let mut estimator = SaemEstimator::new(model, config);
                let results = estimator.fit(&dataset)?;
                
                // Generate comprehensive diagnostics
                let diagnostics = diagnostics::generate_diagnostics(&dataset, &results)?;
                
                // Save final results
                let final_output = Path::new(output_base).join("best_model_analysis");
                output::save_results(&final_output, &results, &diagnostics, &dataset, estimator.model())?;
                
                println!("Final SAEM analysis completed. Results saved to: {:?}", final_output);
            }
            EstimationMethod::Foce | EstimationMethod::FoceI => {
                let mut estimator = FoceEstimator::new(model, config);
                let results = estimator.fit(&dataset)?;
                
                println!("Final FOCE analysis completed. OFV: {:.2}", results.objective_function_value);
            }
        }
    } else {
        println!("No converged models found. Consider:");
        println!("- Increasing iterations");
        println!("- Checking data quality");
        println!("- Trying different initial values");
    }
    
    Ok(())
}

// Run automated analysis
run_automated_analysis("examples/two_compartment_dataset.csv", "automated_results")?;
```

#### Example 4: Production Pipeline with Multiple Models
```rust
use nmodes::*;

fn production_pipeline(dataset_path: &str) -> Result<()> {
    let dataset = Dataset::from_csv(dataset_path)?;
    
    // Phase 1: Quick screening with FOCE
    println!("Phase 1: Model screening with FOCE");
    let model_types = vec![
        ModelType::OneCompartment,
        ModelType::TwoCompartment,
        ModelType::ThreeCompartment,
    ];
    
    let mut screening_results = Vec::new();
    
    for model_type in &model_types {
        let model = CompartmentModel::new(model_type.clone())?;
        let config = EstimationConfig::default()
            .with_method(EstimationMethod::Foce)
            .with_foce_iterations(50);
        
        let mut estimator = FoceEstimator::new(model, config);
        let results = estimator.fit(&dataset)?;
        
        let aic = -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64;
        screening_results.push((model_type.clone(), aic, results.converged));
        
        println!("  {}: AIC = {:.1}, Converged = {}", model_type, aic, results.converged);
    }
    
    // Phase 2: Select top 2 models for detailed analysis
    screening_results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let top_models: Vec<_> = screening_results.iter()
        .filter(|(_, _, converged)| *converged)
        .take(2)
        .map(|(model_type, _, _)| model_type.clone())
        .collect();
    
    println!("\nPhase 2: Detailed analysis of top models with multiple methods");
    let detailed_methods = vec![EstimationMethod::Saem, EstimationMethod::FoceI];
    
    let mut detailed_results = Vec::new();
    
    for model_type in &top_models {
        for method in &detailed_methods {
            let model = CompartmentModel::new(model_type.clone())?;
            let config = EstimationConfig::default()
                .with_method(method.clone())
                .with_iterations(match method {
                    EstimationMethod::Saem => 1500,
                    EstimationMethod::FoceI => 150,
                    _ => 100,
                });
            
            let (aic, converged) = match method {
                EstimationMethod::Saem => {
                    let mut estimator = SaemEstimator::new(model, config);
                    let results = estimator.fit(&dataset)?;
                    let aic = -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64;
                    (aic, results.converged)
                }
                EstimationMethod::FoceI => {
                    let mut estimator = FoceEstimator::new(model, config);
                    let results = estimator.fit(&dataset)?;
                    let aic = -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64;
                    (aic, results.converged)
                }
                _ => unreachable!(),
            };
            
            detailed_results.push((model_type.clone(), method.clone(), aic, converged));
            println!("  {} + {}: AIC = {:.1}, Converged = {}", model_type, method, aic, converged);
        }
    }
    
    // Phase 3: Final model selection
    if let Some((best_model, best_method, best_aic, _)) = detailed_results.iter()
        .filter(|(_, _, _, converged)| *converged)
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
    {
        println!("\nFinal recommendation: {} with {} (AIC: {:.1})", 
                 best_model, best_method, best_aic);
    }
    
    Ok(())
}

// Run production pipeline
production_pipeline("examples/two_compartment_dataset.csv")?;
```

### Advanced Configuration

```rust
use pkpd_saem::estimation::EstimationConfig;

// SAEM configuration
let saem_config = EstimationConfig::default()
    .with_method(EstimationMethod::Saem)
    .with_iterations(5000)           // More iterations for complex models
    .with_burnin(1000)              // Longer burn-in period
    .with_chains(8)                 // More chains for better mixing
    .with_step_size(0.05)           // Smaller steps for better acceptance
    .with_seed(Some(42));           // Reproducible results

// FOCE configuration
let foce_config = EstimationConfig::default()
    .with_method(EstimationMethod::FoceI)
    .with_foce_iterations(200)      // FOCE iterations
    .with_foce_tolerance(1e-8)      // Tighter convergence
    .with_foce_interaction(true);   // Enable interaction terms

// Validate configuration
saem_config.validate()?;
foce_config.validate()?;
```

### Real-World Analysis Examples

#### Dose-Response Analysis
```rust
// Analyze two-compartment dataset with variable dosing
let dataset = Dataset::from_csv("examples/two_compartment_dataset.csv")?;

// Group individuals by dose level
let mut dose_groups = std::collections::HashMap::new();
for (&id, individual) in dataset.individuals() {
    let total_dose = individual.total_dose();
    dose_groups.entry(total_dose as i32).or_insert(Vec::new()).push(id);
}

println!("Dose groups:");
for (dose, ids) in dose_groups {
    println!("{}mg: {} subjects (IDs: {:?})", dose, ids.len(), ids);
}

// Fit model and examine dose-normalized parameters
let model = CompartmentModel::new(ModelType::TwoCompartment)?;
let config = EstimationConfig::default().with_method(EstimationMethod::Foce);
let mut estimator = FoceEstimator::new(model, config);
let results = estimator.fit(&dataset)?;

// Analyze individual parameters by dose group
for (dose, ids) in dose_groups {
    let mut cl_values = Vec::new();
    for id in ids {
        if let Some(params) = results.individual_parameters.get(&id) {
            cl_values.push(params[0].exp()); // CL parameter
        }
    }
    let mean_cl = cl_values.iter().sum::<f64>() / cl_values.len() as f64;
    println!("Dose {}mg: Mean CL = {:.3} L/h", dose, mean_cl);
}
```

#### Covariate Analysis
```rust
// Analyze relationship between demographics and PK parameters
let dataset = Dataset::from_csv("examples/one_compartment_dataset.csv")?;
let model = CompartmentModel::new(ModelType::OneCompartment)?;

// Fit model
let config = EstimationConfig::default().with_method(EstimationMethod::Saem);
let mut estimator = SaemEstimator::new(model, config);
let results = estimator.fit(&dataset)?;

// Analyze covariate relationships
println!("Individual Parameter Analysis:");
println!("{:<4} {:<8} {:<6} {:<8} {:<8} {:<8}", "ID", "Weight", "Age", "CL", "V", "CL/kg");

for (&id, individual) in dataset.individuals() {
    if let Some(params) = results.individual_parameters.get(&id) {
        let weight = individual.get_covariate("WEIGHT").unwrap_or(70.0);
        let age = individual.get_covariate("AGE").unwrap_or(40.0);
        let cl = params[0].exp();
        let v = params[1].exp();
        let cl_per_kg = cl / weight;
        
        println!("{:<4} {:<8.1} {:<6.0} {:<8.3} {:<8.1} {:<8.4}", 
                 id, weight, age, cl, v, cl_per_kg);
    }
}
```

#### Simulation and Validation
```rust
// Simulate data from known parameters for validation
use rand::prelude::*;

fn validate_estimation_method() -> Result<()> {
    // True population parameters
    let true_cl = 2.0;  // L/h
    let true_v = 20.0;  // L
    
    // Simulate dataset (simplified)
    let simulated_data = simulate_pk_data(true_cl, true_v, 50)?; // 50 subjects
    
    // Test both methods
    let methods = vec![
        ("FOCE", EstimationMethod::Foce),
        ("SAEM", EstimationMethod::Saem),
    ];
    
    for (name, method) in methods {
        let model = CompartmentModel::new(ModelType::OneCompartment)?;
        let config = EstimationConfig::default()
            .with_method(method)
            .with_iterations(if matches!(method, EstimationMethod::Saem) { 1000 } else { 100 });
        
        let (estimated_cl, estimated_v) = match method {
            EstimationMethod::Saem => {
                let mut estimator = SaemEstimator::new(model, config);
                let results = estimator.fit(&simulated_data)?;
                (results.fixed_effects[0].exp(), results.fixed_effects[1].exp())
            }
            EstimationMethod::Foce | EstimationMethod::FoceI => {
                let mut estimator = FoceEstimator::new(model, config);
                let results = estimator.fit(&simulated_data)?;
                (results.fixed_effects[0].exp(), results.fixed_effects[1].exp())
            }
        };
        
        let cl_bias = ((estimated_cl - true_cl) / true_cl) * 100.0;
        let v_bias = ((estimated_v - true_v) / true_v) * 100.0;
        
        println!("{}: CL bias = {:.1}%, V bias = {:.1}%", name, cl_bias, v_bias);
    }
    
    Ok(())
}
```

### Custom Model Parameters

```rust
// Get default parameters and modify
let mut params = model.default_parameters();

// Set specific parameter values
params.set_parameter("CL", 1.5)?;
params.set_parameter("V", 25.0)?;

// Validate parameters
model.validate_parameters(&params)?;
```

## Performance Characteristics

### Computational Complexity
- **SAEM Time Complexity**: O(n_individuals × n_iterations × n_mcmc_samples)
- **FOCE Time Complexity**: O(n_individuals × n_iterations × n_observations)
- **Memory Usage**: Linear in dataset size and number of parameters
- **Parallelization**: SAEM uses MCMC sampling across individuals; FOCE is sequential

### Benchmarks
Typical performance on modern hardware:

| Dataset Size | Model Type | SAEM Time | FOCE Time | Memory (MB) |
|--------------|------------|-----------|-----------|-------------|
| 50 individuals, 500 observations | 1-compartment | 2-5 min | 30-60 sec | 50-100 |
| 100 individuals, 1000 observations | 2-compartment | 8-15 min | 2-4 min | 100-200 |
| 200 individuals, 2000 observations | 3-compartment | 25-45 min | 5-10 min | 200-400 |
| 1000 individuals, 10000 observations | 2-compartment | 3-6 hours | 15-30 min | 500-1000 |

### Optimization Tips

#### For SAEM:
1. **Adequate Burn-in**: Use at least 20% of total iterations for burn-in
2. **Multiple Chains**: Use 4-8 chains for robust estimation
3. **Monitor Convergence**: Check parameter stability and log-likelihood trajectory
4. **Step Size Tuning**: Adjust for 40-50% acceptance rate

#### For FOCE:
1. **Start Simple**: Begin with FOCE before trying FOCE-I
2. **Iteration Monitoring**: FOCE typically converges in 50-200 iterations
3. **Gradient Checking**: Monitor gradient norm for convergence assessment
4. **Hessian Conditioning**: Check condition number for numerical stability

#### General:
1. **Method Selection**: Use FOCE for speed, SAEM for robustness
2. **Model Complexity**: Start with 1-compartment models
3. **Data Quality**: Ensure adequate sampling around Cmax and elimination phases
4. **Cross-Validation**: Compare results between methods for consistency

## Validation and Quality Control

### Dataset Validation
The program automatically validates:
- Required columns presence
- Time sequence ordering
- Positive dose amounts
- Reasonable concentration values
- Individual data completeness

### Model Validation
- Parameter bound checking
- Numerical stability monitoring
- Convergence assessment
- Residual analysis

### Convergence Diagnostics
- Parameter trajectory stability
- Log-likelihood convergence
- Between-chain variance (R-hat equivalent)
- Effective sample size estimation

## Error Handling

### Common Issues and Solutions

**Dataset Loading Errors:**
```
Error: Missing required column: TIME
Solution: Ensure your CSV has all required columns (ID, TIME, DV, AMT, EVID)
```

**SAEM Convergence Issues:**
```
Warning: SAEM did not converge after 1000 iterations
Solution: Increase iterations (-i 2000) or adjust burn-in period (-b 400)
```

**FOCE Convergence Issues:**
```
Warning: FOCE did not converge after 100 iterations
Solution: Increase FOCE iterations or check for model misspecification
```

**Numerical Instability:**
```
Error: Integration failed: Numerical instability detected
Solution: Check for extreme parameter values or data outliers
```

**Method-Specific Issues:**
```
FOCE Error: Hessian not positive definite
Solution: Try FOCE-I or switch to SAEM for more robust estimation

SAEM Error: Poor mixing in MCMC chains
Solution: Increase burn-in period or adjust step size
```

### Debugging Tips

1. **Enable Logging**: Set `RUST_LOG=debug` for detailed output
2. **Check Data**: Validate dataset format and values
3. **Method Selection**: Try FOCE first for quick assessment, then SAEM for final analysis
4. **Start Simple**: Begin with fewer iterations to test setup
5. **Monitor Progress**: Watch parameter trajectory (SAEM) or gradient norm (FOCE)
6. **Cross-Validation**: Compare results between methods

## Advanced Features

### Custom Covariate Models
```rust
// Add covariates to the dataset
let mut individual = dataset.get_individual_mut(1)?;
individual.set_covariate("WEIGHT".to_string(), 70.0);
individual.set_covariate("AGE".to_string(), 45.0);
```

### Model Comparison
```rust
// Fit multiple models with different methods
let models = vec![
    ModelType::OneCompartment,
    ModelType::TwoCompartment,
    ModelType::ThreeCompartment,
];

let methods = vec![
    EstimationMethod::Foce,
    EstimationMethod::Saem,
];

for model_type in models {
    for method in &methods {
        let model = CompartmentModel::new(model_type.clone())?;
        let config = EstimationConfig::default()
            .with_method(method.clone())
            .with_iterations(if matches!(method, EstimationMethod::Saem) { 1000 } else { 100 });
        
        let aic = match method {
            EstimationMethod::Saem => {
                let mut estimator = SaemEstimator::new(model, config);
                let results = estimator.fit(&dataset)?;
                -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64
            }
            EstimationMethod::Foce | EstimationMethod::FoceI => {
                let mut estimator = FoceEstimator::new(model, config);
                let results = estimator.fit(&dataset)?;
                -2.0 * results.final_log_likelihood + 2.0 * results.fixed_effects.len() as f64
            }
        };
        
        println!("Model: {:?}, Method: {}, AIC: {:.2}", model_type, method, aic);
    }
}
```

### Simulation Studies
```rust
// Simulate data from known parameters
let true_params = model.default_parameters();
let simulated_data = simulate_dataset(&model, &true_params, n_individuals, dose_schedule)?;

// Fit model to simulated data
let results = estimator.fit(&simulated_data)?;

// Compare estimated vs true parameters
for (i, &estimate) in results.fixed_effects.iter().enumerate() {
    let true_value = true_params.fixed_effects[i];
    let bias = ((estimate - true_value) / true_value) * 100.0;
    println!("Parameter {}: Bias = {:.1}%", i, bias);
}
```

## Statistical Methods

### SAEM (Stochastic Approximation Expectation Maximization)
The implementation follows the methodology described in:
- Kuhn & Lavielle (2005): Maximum likelihood estimation in nonlinear mixed effects models
- Delyon et al. (1999): Convergence of a stochastic approximation version of the EM algorithm

**Key Features:**
1. **Stochastic Approximation**: Robust convergence even with complex models
2. **Expectation Maximization**: Handles missing data and random effects naturally
3. **MCMC Sampling**: Flexible handling of non-linear mixed effects
4. **Adaptive Step Sizes**: Automatic tuning for optimal acceptance rates

### FOCE (First Order Conditional Estimation)
The implementation follows the methodology described in:
- Lindstrom & Bates (1990): Nonlinear mixed effects models for repeated measures data
- Beal & Sheiner (1982): Estimating population kinetics
- Wolfinger & Lin (1997): Two Taylor-series approximation methods

**Key Features:**
1. **First-Order Linearization**: Taylor expansion around conditional modes
2. **Newton-Raphson Optimization**: Fast convergence for individual parameters
3. **Fisher Information Matrix**: Direct calculation of standard errors
4. **Conditional Estimation**: Accounts for individual parameter uncertainty

**FOCE-I Enhancement:**
- Includes interaction terms between random effects and residual error
- More accurate for models with non-constant variance
- Better performance with non-linear observation functions

### Parameter Transformations
- **Log-normal Distribution**: Parameters are estimated on log-scale for positivity
- **Inter-individual Variability**: Modeled using multivariate normal distribution
- **Residual Error**: Proportional error model with log-normal observations

## Testing

### Unit Tests
```bash
# Run all unit tests
cargo test

# Run specific module tests
cargo test models::
cargo test saem::
```

### Integration Tests
```bash
# Run integration tests with example data
cargo test --test integration_tests

# Run with output
cargo test --test integration_tests -- --nocapture
```

### Benchmarks
```bash
# Run performance benchmarks
cargo bench

# Specific benchmark
cargo bench ode_solve
```

## Architecture

The codebase follows a modular design with clear separation of concerns:

```
src/
├── data/           # Dataset loading and validation
│   ├── dataset.rs  # Main dataset structure
│   ├── individual.rs # Individual subject data
│   ├── observation.rs # Observation records
│   └── dosing.rs   # Dosing event handling
├── models/         # Pharmacokinetic models
│   ├── compartment.rs # Base model traits
│   ├── one_compartment.rs
│   ├── two_compartment.rs
│   └── three_compartment.rs
├── solver/         # ODE solving
│   ├── ode.rs      # Solver traits
│   └── runge_kutta.rs # RK4 implementation
├── saem/           # SAEM algorithm
│   ├── algorithm.rs # Main SAEM implementation
│   └── mcmc.rs     # MCMC sampling
├── estimation/     # Estimation configuration
├── diagnostics/    # Model diagnostics
├── output/         # Result formatting
├── validation/     # Data validation
└── main.rs         # CLI interface
```

### Design Principles
- **Modularity**: Each module has a single responsibility
- **Type Safety**: Extensive use of Rust's type system for correctness
- **Error Handling**: Comprehensive error types with context
- **Performance**: Zero-cost abstractions and efficient algorithms
- **Testability**: Isolated components with clear interfaces

## Troubleshooting

### Common Performance Issues

**Slow Convergence:**
- Increase burn-in period: `-b 500`
- Reduce step size for better mixing
- Check for data outliers or model misspecification

**Memory Usage:**
- Large datasets: Process in batches
- Many parameters: Consider model simplification
- Long runs: Monitor memory usage with system tools

**Numerical Issues:**
- Check parameter bounds and initial values
- Verify data quality and outliers
- Consider alternative error models

### Getting Help

1. **Check Logs**: Enable debug logging with `RUST_LOG=debug`
2. **Validate Data**: Run with small subset first
3. **Compare Models**: Start with simpler models
4. **Review Literature**: Consult pharmacokinetic modeling references

## Contributing

### Development Setup
```bash
# Clone repository
git clone <repository-url>
cd nmodes

# Install development dependencies
cargo build

# Run tests
cargo test

# Check formatting
cargo fmt --check

# Run linter
cargo clippy
```

### Code Style
- Follow Rust standard formatting (`cargo fmt`)
- Use meaningful variable names
- Add comprehensive documentation
- Include unit tests for new features
- Maintain backwards compatibility

## References

### Primary Literature
- Kuhn, E., & Lavielle, M. (2005). Maximum likelihood estimation in nonlinear mixed effects models. *Computational Statistics & Data Analysis*, 49(4), 1020-1038.
- Delyon, B., Lavielle, M., & Moulines, E. (1999). Convergence of a stochastic approximation version of the EM algorithm. *Annals of Statistics*, 27(1), 94-128.
- Beal, S., Sheiner, L. B., Boeckmann, A., & Bauer, R. J. (Eds.). (2009). NONMEM User's Guides. Icon Development Solutions.

### Additional Resources
- Mentre, F., & Escolano, S. (2006). Prediction discrepancies for the evaluation of nonlinear mixed-effects models
- Savic, R. M., & Karlsson, M. O. (2009). Importance of shrinkage in empirical Bayes estimates for diagnostics
- Comets, E., Lavenu, A., & Lavielle, M. (2017). Parameter estimation in nonlinear mixed effect models using saemix

## License

MIT License - see LICENSE file for details.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{nmodes,
  title = {NMODES: Nonlinear Mixed Effects Differential Equation Solver},
  author = {NMODES Team},
  year = {2024},
  url = {https://github.com/your-repo/nmodes}
}
```