setwd("/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Graph_estimation")
rm(list=ls(all=TRUE))

library(glasso)
library(glassoFast)
library(sparseMatEst)
library(psych)
library("cglasso")
library("sparsegl")
library(parallel)
library(doParallel)
library(foreach)
library(stabs)

n <- 100  
p <- 5    
M <- 3
scores <- data.frame(matrix(rnorm(n * (p * M), mean = 0, sd = 1), nrow = n, ncol = p * M))
colnames(scores) <-c("f1.1", "f1.2", "f1.3", "f2.1", "f2.2", "f2.3", "f3.1", "f3.2", "f3.3",
              "f4.1", "f4.2", "f4.3", "f5.1", "f5.2", "f5.3")
covariates <- data.frame("pat_group" = c(rep(1,n/2), rep(2, n/2)))
covariates$pat_group <- as.factor(covariates$pat_group) 
scr = FALSE
verbose = TRUE

Mp <- ncol(scores)
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

iU
q

if (n < (q+1)*Mp) {
  warning("The sample size is too small! Network estimate may be unreliable!")
}

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


