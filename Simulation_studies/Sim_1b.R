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
covariates <- data.frame(group= c(rep(0,24),rep(1,26)))
covariates$group <- as.factor(covariates$group)
full_data <- cbind(covariates,scores)


numCores <- detectCores() - 1  # Use all but one core
cl <- makeCluster(15)
registerDoParallel(cl)


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
                                  name_log= "try1")

stopCluster(cl)
save(res, file = "Sim_1.RData")



print( Sys.time() - start )

L = 100 # How many penalization term in the Lasso to try
K = 5
thres.ctrl = c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0) # recognition threshold epsilon_n = thres.ctrl * lambda_n,
verbose = FALSE
tol.abs =1e-4 
tol.rel = 1e-4
eps = 1e-08
name_log = ""

len.t <- length(thres.ctrl)  
n <- nrow(scores)
M <- n_basis
Mp <- ncol(scores)
p <- ceiling(ncol(scores)/M)
if (is.null(covariates)) {
  C <- data.frame(Zeros = rep(0, n))
  iU <-data.frame(rep(1, n))
  colnames(iU) <- "(Intercept)"
  q <- 0
} else {
  numeric_columns <- sapply(covariates, is.numeric)
  covariates[, numeric_columns] <- scale(covariates[, numeric_columns])
  C <- model.matrix(~ ., covariates)
  q <- ncol(C) - 1
  iU <- C
}
#if (n < (q+1)*Mp) {
#  warning("The sample size is too small! Network estimate may be unreliable!")
#}

# G.mat is the optimal adjacency matrix to save
G.mat <- matrix(NA, p, (q+1)*p)

## DEFINITION OF THE DESIGN MATRIX

interM <- data.frame(Intercept = rep(1, n))
temp_groups <- c(0)

if (verbose) {
  cat("Comuputation of the design matrix \n ")
}
for(j in 1:(q + 1)){
  product <- scores * iU[, j]
  if (j != 1) {
    original_colnames <- colnames(product)
    iU_name <- colnames(iU)[j]
    new_colnames <- paste(iU_name, ":", original_colnames, sep = "")
    colnames(product) <- new_colnames
  }
  for (i in 1:p){temp_groups <- c(temp_groups, rep(i+(p*(j-1)),M))}
  interM <- cbind(interM, product)
}
interM <- interM[, -1]
temp_groups <- temp_groups[-1]

# Q: A cosa servono d.array e  d ? Quali sono le dimensioni corrette?
# Per il momento li ho sostuiti com dei vettori di uni della dimensione coretta (?)
d.array <- matrix(1, nrow=p, ncol=(p-1)*M*(q+1))
d.out <- list()
norm.adj <- rep(NA,p)
for(k in 1:p){
  d.out[[k]] <- d.array[k,]
  norm.adj[k] <- norm(d.out[[k]],"2")
}

P.def <- matrix(0, 2*(p-1)*M, M)
Q.def <- matrix(0.1, 2*(p-1)*M, M)
U.def <- matrix(0.01, 2*(p-1)*M, M)

j=1
cat(paste("Processing node ", j,"\n")) 
                                                                     
jth.range.y <- (j-1)*M+(1:M)
A.Y.out <- as.matrix(interM[, jth.range.y])
jth.range.x <- c(jth.range.y,(j+p-1)*M+(1:M))
A.X.out <- as.matrix(interM[, -jth.range.x])
groups <- temp_groups[-jth.range.x]
                                                                     
P.out <- P.def; Q.out <- Q.def; U.out <- U.def
                                                                     
lambda.max <- lambda.sup(A.X.out, A.Y.out)
lambdas <- exp(seq(log(lambda.max), log(1e-4), length.out = L))
                                                                     
thresholds <- c()
for(ind.t in 1:len.t){
  thresholds <- c(thresholds, lambdas * thres.ctrl[ind.t])}
thresholds <- unique(thresholds)[order(unique(thresholds))]
                                                                     
if(j==p){add.out=TRUE}else{add.out=FALSE}
M.out <- M
K.out <- K

start <- Sys.time()

opt_results <- BayesianOptimization(
  FUN = optimize_hyperparameters_wrapper,
  bounds = list(lambda = c(0,lambda.max),
                threshold = c(0,2)),
  init_points = 10,  # Number of random initial points
  n_iter = 30,       # Number of optimization iterations
  acq = "ucb", 
  kappa = 3, 
  verbose = verbose
)

print( Sys.time() - start )

best_lambda <- opt_results$Best_Par["lambda"]
best_threshold <- opt_results$Best_Par["threshold"]

N.hat.optimal <-  compute_Nj(A.X=A.X.out, A.Y=A.Y.out, M=M.out, d=d.out[[j]],
                             P=P.out, Q=Q.out, U=U.out, lambda =best_lambda, threshold=best_threshold, add=add.out)

N.hat.optimal$N.hat


