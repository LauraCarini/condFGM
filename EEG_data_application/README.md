# EEG Data Application for condFGM

This folder contains the real-world application of the conditional Functional Graphical Model (condFGM) method to electroencephalography (EEG) data for studying genetic predisposition to alcoholism.

## Dataset Description

The EEG data analyzed in this application was obtained from the repository associated with the work by Zhao, B., Wang, Y. S. & Kolar, M. (2022), "FuDGE: A method to estimate a functional differential graph in a high-dimensional setting", *Journal of Machine Learning Research* 23(82), 1–82.

**Data Source**: https://github.com/boxinz17/FuDGE/tree/master/EEG/

### Dataset Background

The dataset originates from a large study examining EEG correlates of genetic predisposition to alcoholism. It contains EEG recordings that capture neural activity differences between alcoholic and control subjects, making it an ideal test case for functional differential graph estimation methods.

**Key Characteristics:**
- **Electrodes**: 64 electrodes placed on the scalp
- **Sampling Rate**: 256 Hz
- **Duration**: 1 second per trial (256 time points)
- **Subjects**: Two groups - alcoholic and control
- **Study Focus**: Genetic predisposition to alcoholism
- **Original Source**: Neurodynamics Laboratory, State University of New York Health Center, Brooklyn

### Data Files

This application uses two versions of the EEG dataset:

#### 1. `alco_array.Rdata` (Raw Files)
- Contains the original, unprocessed EEG recordings
- Includes raw signal data from all 64 electrodes
- Suitable for custom preprocessing pipelines
- Format: R data array

#### 2. `alco_filtered_array.Rdata` (Preprocessed Data)
- Contains filtered and preprocessed EEG signals
- Ready for direct analysis with condFGM
- Preprocessing likely includes noise reduction and artifact removal
- Format: R data array

## Application Overview

This real-world application demonstrates the effectiveness of condFGM in:

1. **Functional Connectivity Analysis**: Identifying brain network differences between alcoholic and control subjects
2. **Differential Graph Estimation**: Detecting condition-specific functional connectivity patterns
3. **High-Dimensional Analysis**: Handling the complexity of multi-channel EEG data
4. **Temporal Dynamics**: Capturing functional relationships across time

## Scientific Context

### Alcoholism and Brain Function

Alcohol use disorders (AUDs) are associated with significant alterations in brain function, including:
- Impaired cognitive processing
- Reduced neural connectivity
- Disrupted brain rhythms
- Altered functional connectivity between brain regions

### EEG as a Biomarker

EEG signals provide valuable insights into:
- Real-time neural activity
- Brain network connectivity
- Functional differences between populations
- Potential biomarkers for neurological conditions

## Data Analysis Pipeline

### Prerequisites

```r
# Required R packages
library(condFGM)  # Main package for analysis
# Additional packages as required by your specific scripts
```

### Loading the Data

```r
# Load raw data
load("alco_array.Rdata")

# Or load preprocessed data
load("alco_filtered_array.Rdata")
```

### Typical Analysis Workflow

1. **Data Exploration**
   - Visualize EEG signals
   - Check data dimensions and structure
   - Explore group differences

2. **Preprocessing** (if using raw data)
   - Filter artifacts and noise
   - Standardize signals
   - Handle missing data

3. **condFGM Analysis**
   - Set up functional basis (e.g., Fourier, B-spline)
   - Configure model parameters
   - Estimate differential graphs
   - Validate results

4. **Interpretation**
   - Identify significant connectivity differences
   - Map results to brain regions
   - Interpret biological significance

## Expected Results

The condFGM analysis of this EEG dataset should reveal:

- **Differential Connectivity Patterns**: Brain regions showing different functional connectivity between alcoholic and control groups
- **Network Disruptions**: Areas where alcoholism affects normal brain network function
- **Temporal Dynamics**: How connectivity differences evolve over the 1-second recording period
- **Biomarker Identification**: Potential EEG-based markers for alcohol use disorders

## Computational Considerations

### Memory Requirements
- **Raw Data**: Moderate memory usage for 64 channels × 256 time points
- **Processing**: Memory requirements scale with the number of subjects and basis functions used

### Processing Time
- Analysis time depends on:
  - Number of subjects in the dataset
  - Choice of functional basis and parameters
  - Optimization settings in condFGM

## File Structure

```
EEG_data_application/
├── README.md                    # This file
├── alco_array.Rdata            # Raw EEG data
├── alco_filtered_array.Rdata   # Preprocessed EEG data
├── [analysis_script.R]         # Main analysis script
├── [preprocessing.R]           # Data preprocessing (if applicable)
├── [visualization.R]           # Results visualization
└── [results/]                  # Output directory
    ├── figures/                # Generated plots
    ├── differential_graphs/    # Estimated graph structures
    └── summary_statistics/     # Analysis summaries
```

## Citation and References

### Original Dataset
If you use this EEG dataset, please cite the original source:
- **FuDGE Paper**: Zhao, B., Wang, Y. S. & Kolar, M. (2022), 'FuDGE: A method to estimate a functional differential graph in a high-dimensional setting', *Journal of Machine Learning Research* 23(82), 1–82.
- **Data Repository**: https://github.com/boxinz17/FuDGE/tree/master/EEG/

### Original EEG Study
The EEG data originates from studies on genetic predisposition to alcoholism. Please also cite relevant papers from the Neurodynamics Laboratory at SUNY Health Center.

### condFGM Method
If you use the condFGM method for analysis, please cite:
[Your condFGM paper citation]

## Data Usage and Ethics

- This dataset is provided for research purposes
- Follow appropriate ethical guidelines when analyzing neurological data
- Ensure compliance with data sharing agreements from the original sources
- Respect privacy and confidentiality of study participants

## Troubleshooting

### Common Issues

1. **Loading Data**: Ensure R can locate the `.Rdata` files in your working directory
2. **Memory Errors**: Consider using data subsets for initial testing
3. **Package Dependencies**: Install all required packages before running analysis
4. **Data Format**: Verify the structure of loaded data matches expected format

### Support

For technical issues with:
- **condFGM Method**: Contact the method developers
- **Original Dataset**: Refer to the FuDGE repository or original papers
- **General Questions**: Check documentation and literature

## Reproducibility

To ensure reproducible results:
- Use consistent R versions and package versions
- Set random seeds where applicable
- Document parameter choices and preprocessing steps
- Save intermediate results for validation

## License and Distribution

Please respect the licensing terms of:
- The original EEG dataset
- The FuDGE repository
- The condFGM package

---

**Note**: This application demonstrates the practical utility of condFGM for neuroscience research and provides a benchmark for comparing functional differential graph estimation methods on real neurological data.
