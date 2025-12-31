# Simulation Studies for condFGM

This repository contains the simulation pipeline for reproducing the results presented in our paper. The simulation studies are designed to evaluate the performance of the conditional Functional Graphical Model (condFGM) across different scenarios and compare it with competing algorithms.

## Repository Structure

```
Simulation_studies/
├── Function/                         # Helper functions for adjacency matrix simulation
├── Step 1 # Run the simualtion analysis in one setting with computational cost estimation and litterature comparison
├── Step 2 # Run simulations across different scenarios
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

The main configuration file in eahc simulation folders are the `config_*.yaml` files. Key parameters include:

#### Data Configuration
- save_path: Directory where simulation results will be saved
- name_output: Prefix for saved result files
- model_g1 and model_g2: Model specifications for the two groups
- red_number: Define the number of edges affected by model_g1 and model_g2

#### Algorithm Parameters
- rec_basis_type: reconstruction basis type either "fourier" or "bsplines"
- rec_basis_number: number of basis functions
- M: number of functional principal components to retain
- L: number of lambda values to test
- K: number of folds in cross-validation optimization
- thres_ctrl: threshold values to test

### Step 1: Sigle setting Simulation and Algorithm Comparison

Run the simualtion analysis in one setting varying the number of nodes and samples with computational cost estimation and litterature comparison

#### 1.1 Run Main Simulation Pipeline
```bash
bash Sbatch_simulations_luncher_test.sh
```
This script automatically evaluates the pipeline performance across different combinations of `p` and `n` specified in the file at **lines 17-19** .

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

Fixing the number of nodes and sample, run simulations across different scenarios. For each configuration file you want to analyze:

```bash
$CONFIG_FILE="config_S1.yaml"
```

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
The simulation scenarios are defined by specifying models for the two groups in lines 31-32 of `config_*.yaml`. Available options are detailed in the README file located in the `Function/` folder.

## Output Structure

Results are saved in the directory specified by `save_path` with the following structure:

```
results/
└── p[p]_n[n]_n[n].rds
    ├── seed_[iteration]/
    ├── test_results_metices_[name_output].csv/
    └── litt_comp_test_results_metices_[name_output].csv/
```

## Function Folder

The `Function/` folder contains utility functions for:
- Adjacency matrix simulation
- Data generation helpers

## Citation

If you use this simulation pipeline in your research, please cite our paper:

---

**Note**: This simulation pipeline is designed for high-performance computing environments with SLURM job scheduling. Modifications may be needed for other computing environments.
