# Simulation Studies for condFGM

This repository contains the simulation pipeline for reproducing the results presented in our paper. The simulation studies are designed to evaluate the performance of the conditional Functional Graphical Model (condFGM) across different scenarios and compare it with competing algorithms.

## Repository Structure

```
Simulation_studies/
├── Function/                         # Helper functions for adjacency matrix simulation
├── Step 1 # Run the simualtion analysis in one setting with computational cost estimation and litterature comparison
├── Step 2 # Run simulations across different scenarius
├── Plot_simulations_resulrs.R        # Results visualization
└── Sumulation_full_pipeline.txt      # Terminal command pipeline
```

## Prerequisites

- R (version 4.0 or higher)
- SLURM job scheduler (for parallel execution)
- Required R packages (will be installed automatically by the scripts)

## Quick Start

For a complete pipeline execution, refer to the commands in `Sumulation_full_pipeline.txt`.

## Detailed Pipeline Instructions

### Configuration Setup

The main configuration file is `config_template.yaml`. Key parameters include:

#### Data Configuration
- **Line 17**: `data_save_folder` - Directory where simulation results will be saved
- **Line 41**: `name` - Prefix for saved result files
- **Lines 31-32**: Model specifications for the two groups

#### Algorithm Parameters
- **`rec_basis_type`**: `"fourier"` (reconstruction basis type)
- **`rec_basis_number`**: `15` (number of basis functions)
- **`M`**: `5` (number of functional principal components to retain)
- **`L`**: `100` (number of lambda values to test)
- **`K`**: `5` (number of folds in cross-validation optimization)
- **`thres_ctrl`**: `[0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]` (threshold values to test)

#### Simulation Parameters
- **Lines 17-19**: Specify different values of `p` (number of variables) and `n` (sample size)

### Step 1: Initial Simulation and Algorithm Comparison

#### 1.1 Run Main Simulation Pipeline
```bash
bash Sbatch_simulations_luncher_test.sh
```
This script automatically evaluates the pipeline performance across different combinations of `p` and `n` specified in the configuration file.

#### 1.2 Process Computational Time (Main Algorithm)
```bash
Rscript Process_logs_comp_time.R
```
This extracts and analyzes the computational time for each iteration of the main algorithm.

#### 1.3 Run Competing Algorithms
```bash
bash Sbatch_run_litt_comp_analysis.sh
```
Executes comparison with competing methods.

#### 1.4 Process Computational Time (Competing Algorithms)
```bash
Rscript Process_logs_comp_time_litt_comp.R
```
Analyzes computational time for FuDGE and other competing algorithms.

### Step 2: Detailed Analysis and Results Generation

For each configuration file you want to analyze:

#### 2.1 Generate Simulation Data
```bash
Rscript Data_simulator_tests.R "$CONFIG_FILE"
```

#### 2.2 Run Algorithm in Parallel
```bash
bash Sbatch_parallel_luncher_test.sh "$CONFIG_FILE"
```

#### 2.3 Verify Job Completion
Ensure all parallel jobs have completed successfully before proceeding.

#### 2.4 Run Multiple Iterations
To reproduce the manuscript results, repeat steps 2.1-2.2 for **10 iterations** by:
1. Modifying the iteration parameter (line 5 in the configuration file)
2. Running the pipeline for each iteration value (1-10)

#### 2.5 Check Results
```bash
Rscript Check_results_screening_procedure_tests.R "$CONFIG_FILE"
```
This validates the results for the specified configuration.

### Step 3: Generate Figures

After completing all simulations:

```bash
Rscript Plot_simulations_resulrs.R
```

This script generates all figures presented in the manuscript.

## Scenario Configuration

### Model Types
The simulation scenarios are defined by specifying models for the two groups in lines 31-32 of `config_template.yaml`. Available options are detailed in the README file located in the `Function/` folder.

### Parameter Customization

You can customize the following parameters based on your computational resources and research needs:

- **Sample sizes (`n`)**: Adjust based on your computational capacity
- **Number of variables (`p`)**: Modify to test different dimensionalities
- **Lambda grid (`L`)**: Increase for finer tuning, decrease for faster execution
- **Cross-validation folds (`K`)**: Standard is 5, but can be adjusted
- **Threshold values (`thres_ctrl`)**: Customize based on your specific requirements

## Output Structure

Results are saved in the directory specified by `data_save_folder` with the following structure:

```
results/
├── simulation_results_[name]_p[p]_n[n]_iter[i].rds
├── computational_times/
├── figures/
└── validation_reports/
```

## Computational Requirements

- **Storage**: Ensure sufficient disk space in the specified `data_save_folder`
- **Memory**: Requirements vary with `p` and `n` values
- **Time**: Full pipeline execution time depends on parameter choices and available cores

## Troubleshooting

1. **Job Failures**: Check SLURM logs in the output directory
2. **Memory Issues**: Reduce `p` or `n` values, or increase allocated memory
3. **Missing Dependencies**: Ensure all R packages are installed
4. **Configuration Errors**: Validate YAML syntax in configuration files

## Reproducibility Notes

- Set random seeds for reproducible results
- Use the same R version and package versions across runs
- Ensure consistent computational environment for timing comparisons

## Function Folder

The `Function/` folder contains utility functions for:
- Adjacency matrix simulation
- Data generation helpers
- Algorithm implementations

Refer to the README in the `Function/` folder for detailed descriptions of available simulation models and their parameters.

## Citation

If you use this simulation pipeline in your research, please cite our paper:

[Add your paper citation here]

## Support

For questions or issues with the simulation pipeline, please open an issue in this repository or contact [your contact information].

## License

[Add license information]

---

**Note**: This simulation pipeline is designed for high-performance computing environments with SLURM job scheduling. Modifications may be needed for other computing environments.
