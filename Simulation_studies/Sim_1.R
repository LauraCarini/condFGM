setwd("/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies")
rm(list=ls(all=TRUE))

library(parallel)
library(doParallel)
library(foreach)
#install.packages("rBayesianOptimization")
library(rBayesianOptimization)

data.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"
func.path <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Graph_estimation"
save.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"
runtime.path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1"

source(paste(func.path,"cFGGM_functions_two_groups.R", sep="/"))

scores <- read.csv(paste(data.path, "freq_scores_reord.csv", sep ="/"))[, -1]
n_basis <- 8
n_nodes <- ncol(scores)/n_basis
n_samples <- nrow(scores)
names <- rep(NA, ncol(scores))
for(l in 1:ncol(scores)){
  names[l] <- paste("f",ceiling(l/n_basis) ,".",l%%n_basis, sep ="")
}
colnames(scores) <- names
covariates <- data.frame(group= c(rep(1,24),rep(0,26)))
covariates$group <- as.factor(covariates$group, levels= c(0,1))
full_data <- cbind(covariates,scores)


numCores <- detectCores() - 1  # Use all but one core
cl <- makeCluster(numCores)
registerDoParallel(cl)

L = 10
K = 5
thres.ctrl = c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0)

G.mat <- Fast_FGGReg_diff_two_groups_SCV(scores, # functional score on a defined basis, nrow: subjects; ncol: functions*n_basis.
                                         n_basis = n_basis, #Number of bases considered 
                                         covariates  = covariates, #Additional covariates to regress on
                                         L = L, # How many penalization term in the Lasso to try
                                         K = K,
                                         thres.ctrl = thres.ctrl, # recognition threshold epsilon_n = thres.ctrl * lambda_n,
                                         verbose = FALSE,
                                         tol.abs =1e-4 ,
                                         tol.rel = 1e-4,
                                         eps = 1e-08,
                                         name_log= paste(save.path,"/logs/try2603", sep=""))

stopCluster(cl)
save(G.mat, file = paste(save.path,"res_2603.RData", spe="/"))


print( Sys.time() - start )
