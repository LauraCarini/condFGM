# File name: cFGGM_functions
## Function name: FGGReg_cov_estimation
#### Purpose
This function estimates a covariance network from the functional scores on a defined basis using sparse group lasso regression. 
#### Input
1. scores (Dataframe) – Functional score on a defined basis.
    - Rows represent subjects (nrow: n_samples).
    - Columns represent the scores for each of the functions considered (ncol: functions*n_basis) 



