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

# To be adapted to ours
CV_lambda <-
  function(x = NULL, # expression data, nrow: subjects; ncol: features.
           known_ppi = NULL, # previous known PPI
           covariates = NULL,       #covariates to regress on
           lambda = NULL,
           nlambda = 5,
           lambda.min.ratio = 0.1, # ratio for minial lambda in a grid search
           scr = TRUE, # optional screening to speed up
           gamma = NULL,
           rho = NULL,
           weight = 1.1,
           eta = 0,
           verbose = FALSE,
           eps = 1e-08,
           seed = NULL, # random seed for K fold sample split
           K = 5, # K fold to optimize lambda
           crit.cv = c("PBIC", "BIC", "loglik", "ploglik", "AIC",  "EBIC"),
           trace = c("progress", "print", "none"),
           ...) {
    n <- nrow(x)
    p <- ncol(x)
    

    # Initial checks on the covariates
    if (is.null(covariates)) {
      covariates = data.frame(Zeros = rep(0,n))
    }

    X = scale(x, center = TRUE, scale = TRUE)
    if(verbose){
      cat("Computing correlation matrix \n ")
    }
    S = cor(X)
    if(verbose){
      cat("Done with the matrix computation \n")
    }

    # Definition of the lambda sequence
    if (!is.null(lambda)){
      nlambda = length(lambda)}
    if (is.null(lambda)) {
      if (is.null(nlambda)) {
        nlambda = 15 }
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio = 0.1 }
      lambda.max = max(max(S - diag(p)), -min(S - diag(p)))
      lambda.min = lambda.min.ratio * lambda.max
      lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
      rm(lambda.max, lambda.min, lambda.min.ratio)
      gc()
    }
  
    
    #Match values
    crit.cv = match.arg(crit.cv)
    trace = match.arg(trace)
    lambda = sort(lambda)

    #Initialize the metrics
    CV_nLogLikelihood = CV_pNLogLikelihood = CV_Pbic.score = matrix(0, length(lambda), K)
    CV_ebic.score = CV_aic.score = CV_bic.score =matrix(0, length(lambda), K)
    #CV_density = CV_nLogLikelihood
    
    # set progress bar
    if (trace == "progress") {
      progress = txtProgressBar(max = K * length(lambda), style = 3)
    }
    # No need to create folds if K = 1
    if (K == 1) {
      # set sample size
      X.train = X
      empcov.train = S
      t = 1
      # loop over all tuning parameters
      for (loop.lambda in 1:length(lambda)) {
        # compute the penalized likelihood precision matrix
        # estimator
        fit = GGM_Estimation( x=X.train,
          known_ppi=known_ppi,
          covariates=covariates,
          scr=scr,
          gamma=gamma,
          lambda=lambda[loop.lambda],
          rho=rho,
          weight=weight,
          eta=eta,
          verbose=verbose,
          eps = eps
        )
        
        
        siginv = fit$invcov
        no.edge = (sum(abs(siginv) > eps)) / 2
        #CV_density[loop.lambda, t] = sparsePrec(fit$Adj)
        #Compute the observed negative loglikelihood
        CV_nLogLikelihood[loop.lambda, t] = tr(empcov.train %*% siginv) - 
                                  determinant(siginv, logarithm = T)$modulus
        
        #Compute the observed penalized negative loglikelihood
        CV_pNLogLikelihood[loop.lambda, t] = CV_nLogLikelihood[loop.lambda, t] + fit$lambda *
          sum(abs((1 - diag(p)) * siginv))
        
        #Compute the different metrices
        CV_bic.score[loop.lambda, t] = tr(empcov.train %*% siginv) -
          determinant(siginv, logarithm = T)$modulus + log(n) * no.edge / n
        CV_Pbic.score[loop.lambda, t] = CV_pNLogLikelihood[loop.lambda, t] + log(n) * no.edge / n
        CV_ebic.score[loop.lambda, t] = CV_bic.score[loop.lambda, t] + 4 * log(p) * no.edge / n
        CV_aic.score[loop.lambda, t] = CV_nLogLikelihood[loop.lambda, t] + 2 * no.edge / n
        
        # update progress bar
        if (trace == "progress") {
          setTxtProgressBar(progress, loop.lambda + (t - 1) *
                              length(lambda))
          
          # if not quiet, then print progress lambda
        } else if (trace == "print") {
          cat("\nFinished lambda = ", paste(loop.lambda, sep = ""))
        }
      }
      
      # if not quiet, then print progress kfold
      if (trace == "print") {
        cat("\nFinished fold ", paste(t, sep = ""))
      }
    }
    else {
      # designate folds and shuffle -- ensures randomized folds
      if (!is.null(seed)) {
        set.seed(seed)
      }
      ind = sample(n)
      
      # parse data into folds and perform CV
      for (t in 1:K) {
        # training set
        leave.out = ind[(1 + floor((t - 1) * n / K)):floor(t * n / K)]
        X.train = x[-leave.out, , drop = FALSE]
        covariates.train = covariates[-leave.out, , drop = FALSE]
        # # sample covariances
        X.train = scale(X.train, center = TRUE, scale = TRUE)
        empcov.train = cor(X.train)
        
        # loop over all tuning parameters
        for (loop.lambda in 1:length(lambda)) {
          # compute the penalized likelihood precision matrix
          # estimator
          fit = GGM_Estimation( x=X.train,
          known_ppi=known_ppi,
          covariates=covariates.train,
          scr=scr,
          gamma=gamma,
          lambda=lambda[loop.lambda],
          rho=rho,
          weight=weight,
          eta=eta,
          verbose=verbose,
          eps = eps
        )
        
          siginv = fit$invcov
          
          no.edge = (sum(abs(siginv) > eps)) / 2
          #CV_density[loop.lambda, t] = sparsePrec(fit$Adj)
          # compute the observed negative validation loglikelihood
          CV_nLogLikelihood[loop.lambda, t] = tr(empcov.train %*% siginv) - 
                                    determinant(siginv, logarithm = T)$modulus
          CV_pNLogLikelihood[loop.lambda, t] = CV_nLogLikelihood[loop.lambda, t] + 
                              fit$lambda *sum(abs((1 - diag(p)) * siginv))
          CV_bic.score[loop.lambda, t] = tr(empcov.train %*% siginv) - 
            determinant(siginv, logarithm = T)$modulus + log(n) * no.edge / n
          CV_Pbic.score[loop.lambda, t] = CV_pNLogLikelihood[loop.lambda, t] + 
                                          log(n) * no.edge / n
          CV_ebic.score[loop.lambda, t] = CV_bic.score[loop.lambda, t] + 4 *
                                          log(p) * no.edge / n
          CV_aic.score[loop.lambda, t] = CV_nLogLikelihood[loop.lambda, t] + 
                                          2 * no.edge / n
          # update progress bar
          if (trace == "progress") {
            setTxtProgressBar(progress, loop.lambda + (t - 1) * length(lambda))
            # if not quiet, then print progress lambda
          } else if (trace == "print") {
            cat("\nFinished lambda = ", paste(loop.lambda, sep = ""))
          }
        }
        # if not quiet, then print progress kfold
        if (trace == "print") {
          cat("\nFinished fold ", paste(t, sep = ""))
        }
      }
    }
    
    if (crit.cv == "loglik") {
      CV_errors = CV_nLogLikelihood
    }
    if (crit.cv == "ploglik") {
      CV_errors = CV_pNLogLikelihood
    }
    if (crit.cv == "PBIC") {
      CV_errors = CV_Pbic.score
    }
    if (crit.cv == "BIC") {
      CV_errors = CV_bic.score
    }
    if (crit.cv == "EBIC") {
      CV_errors = CV_ebic.score
    }
    if (crit.cv == "AIC") {
      CV_errors = CV_aic.score
    }
    
    allMetrics = cbind(
      lambda,
      #apply(CV_density, 1, mean),
      apply(CV_nLogLikelihood, 1, mean),
      apply(CV_pNLogLikelihood, 1, mean),
      apply(CV_bic.score, 1, mean),
      apply(CV_Pbic.score, 1, mean),
      apply(CV_ebic.score, 1, mean),
      apply(CV_aic.score, 1, mean)
    )
    colnames(allMetrics) <-
      c(
        "lambda",
        #"density",
        "nLogLikelihood",
        "pNLogLikelihood",
        "BIC",
        "PBIC",
        "EBIC",
        "AIC"
      )
    
    # determine optimal tuning parameters
    AVG = apply(CV_errors, 1, mean)
    STD = apply(CV_errors, 1, sd)
    if (K == 1) {
      bestLambda = lambda[which.min(AVG)]
      minerror = min(AVG)
    } else {
      minerror = min(AVG)
      minerror_sd = STD[which.min(AVG)]
      bestLambda = lambda[which.min(abs(AVG - (minerror + minerror_sd)))]
    }
    
    results = list(
      lambda = lambda,
      bestLambda = bestLambda,
      min.error = minerror,
      #CV_density = CV_density,
      avg.error = AVG,
      cv.error = CV_errors,
      CV_nLogLikelihood = CV_nLogLikelihood,
      CV_pNLogLikelihood = CV_pNLogLikelihood,
      CV_bic.score = CV_bic.score,
      CV_Pbic.score = CV_Pbic.score,
      CV_ebic.score = CV_ebic.score,
      CV_aic.score = CV_aic.score,
      allMetrics = allMetrics
    )
    class(results) = "CVahglasso"
    return(results)
    
  }



# To be adapted to ours
CV_lambda_ppi <-
  function(x, # expression data, nrow: subjects; ncol: features.
           known_ppi = NULL, # previous known PPI
           covariates = NULL,       #covariates to regress on
           lambda = NULL,
           nlambda = 15,
           lambda.min.ratio = 0.1, # ratio for minial lambda in a grid search
           scr = TRUE, # optional screening to speed up
           gamma = NULL,
           rho = NULL,
           weight = 0.1,
           eta = 0,
           verbose = FALSE,
           eps = 1e-08,
           seed = NULL, # random seed for K fold sample split
           K = 5, # K fold to optimize lambda
           crit.cv = c("PBIC", "BIC", "loglik", "ploglik", "AIC",  "EBIC"),
           trace = c("progress", "print", "none"),
           ...) {
    non_zero_columns <- colnames(known_ppi)[apply(known_ppi != 0, 2, any)]
    known_ppi_new <- known_ppi[non_zero_columns, non_zero_columns]
    x_new = x[,non_zero_columns]
    return(CV_lambda(x = x_new,
           known_ppi = known_ppi_new, 
           covariates = covariates,       
           lambda = lambda,
           nlambda = nlambda,
           lambda.min.ratio = lambda.min.ratio, 
           scr = scr, 
           gamma = gamma,
           rho = rho,
           weight = weight,
           eta = eta,
           verbose = verbose,
           eps = eps,
           seed = seed,
           K = K,
           crit.cv = crit.cv,
           trace = trace))
    
  }

# PARALLELIZED

GGReg_cov_estimation <- function(Z0, # expression data (residuals on the mean), nrow: subjects; ncol: features.
                                 known_ppi = NULL, # previously known PPI 
                                 covariates = NULL, # covariates to regress on
                                 scr = TRUE, # optional screening to speed up
                                 gamma = NULL, # Person correlation threshold for the screening
                                 lambda_prec = NULL, # Penalization term in the Lasso
                                 lambda_prec_type = "1se", # or "min"
                                 weight = 1.1, # Multiplicative factor for the penalization term when prior knowledge is not available
                                 asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
                                 verbose = FALSE,
                                 eps = 1e-08) {
  

  n <- nrow(Z0)
  p <- ncol(Z0)
  Ip <- diag(rep(1, p))

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

  if (n < (q+1)*p) {
    warning("The sample size is too small! Network estimate may be unreliable!")
  }
  
  if (!is.null(known_ppi)) {
    one <- known_ppi
    one[one > 0] <- 1
    rownames(one) <- colnames(Z0)
    penaltyfactor <- 1 - known_ppi
  } else {
    one <- matrix(0, p, p)
    rownames(one) <- colnames(Z0)
    penaltyfactor <- matrix(1, p, p)
  }

  if (scr) {
    if (verbose) {
      cat("Computing correlation matrix \n")
    }
    protdata_matrix <- as.matrix(Z0)
    Scor <- cor(protdata_matrix)
    if (verbose) {
      cat("Done computing correlation matrix \n")
    }
    scr_index <- matrix(1, p, p)
    if (is.null(gamma)) {
      gamma <- quantile(abs(Scor), probs = 0.1)
    }
    scr_index[abs(Scor) <= gamma] <- 0
    diag(scr_index) <- 0
    colnames(scr_index) <- colnames(Z0)
    rownames(scr_index) <- colnames(Z0)
  }

  interM <- data.frame(Intercept = rep(1, n))
  temp_groups <- c(0, rep(1, ncol(Z0)))
  temp_w_group <- c(0)


  # Register parallel backend for the first loop
  numCores <- detectCores()
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
    cat("Comuputation of the design matrix \n ")
  }
  # Parallelize the first loop - Need to be tuned 
  interM_results <- foreach(j = 1:(q + 1), .combine = 'cbind') %dopar% {
    product <- Z0 * iU[, j]
    if (j != 1) {
      original_colnames <- colnames(product)
      iU_name <- colnames(iU)[j]
      temp_groups <- rep(j, length(original_colnames))
      temp_w_group <- 1
      new_colnames <- paste(iU_name, ":", original_colnames, sep = "")
      colnames(product) <- new_colnames
    }
    list(product = product, temp_groups  = temp_groups , temp_w_group=temp_w_group)
  }

  stopCluster(cl)
  groups <- c()
  w_group <- c()
  for(j in 1:(q + 1)){
    groups <- c(groups, interM_results[,j]$temp_groups)
    w_group <- c(w_group, interM_results[,j]$temp_w_group)
    interM <- cbind(interM, interM_results[,j]$product)
    }

 # FROM HERE ON
  Delta <- matrix(0, p, ncol(interM))
  colnames(Delta) <- colnames(interM)
  rownames(Delta) <- colnames(Z0)
  Sigma_hat <- rep(1, p)
  names(Sigma_hat) <- colnames(Z0)
  No_sim_Delta_hat <- matrix(0, p, ncol(interM))
  colnames(No_sim_Delta_hat) <- colnames(interM)
  rownames(No_sim_Delta_hat) <- colnames(Z0)

  # Register parallel backend for the second loop
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
        cat("Comuputation of adjacency matrix \n ")
      }

  results <- foreach(i = 1:p, .combine = 'cbind', .packages = c('glmnet', 'sparsegl')) %dopar% {
    Y <- matrix(Z0[, i], ncol = 1)
    protname <- colnames(Z0)[i]
    pattern <- paste0(":", protname, "\\b")
    indexes_to_select <- !grepl(pattern, colnames(interM))
    indexes_to_select[colnames(interM) == protname] <- FALSE
    indexes_to_select_penalty <- rep(TRUE, ncol(Z0))
    indexes_to_select_penalty[colnames(Z0) == protname] <- FALSE
    if (scr) {
      zero_col_indices <- which(scr_index[i, ] == 0)
      zero_col_names <- colnames(scr_index)[zero_col_indices]
      for (j in zero_col_names) {
        pattern <- paste0(":", j, "\\b")
        indexes_to_select[grepl(pattern, colnames(interM))] <- FALSE
        indexes_to_select[colnames(interM) == j] <- FALSE
        indexes_to_select_penalty[colnames(Z0) == j] <- FALSE
      }
    }
    indexes_to_select_reg <- indexes_to_select
    indexes_to_select_reg[1] <- FALSE
    Xmat <- interM[, indexes_to_select_reg]
    Xgroup <- groups[indexes_to_select_reg]
    
    w_sparsity <- rep(weight, ncol(Xmat))
    w_sparsity[Xgroup == 1] <- penaltyfactor[i, indexes_to_select_penalty]
    w_group[which(w_group==1)] <- sqrt(sum(Xgroup == 1))

    
    if (is.null(lambda_prec)) {
      set.seed(2024)
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, intercept = FALSE)
      #mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup,
      #                  pf_sparse = w_sparsity, asparse = asparse, intercept = FALSE)
    } else {
      set.seed(2024)
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, lambda = lambda_prec, intercept = FALSE)
      #mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup,
      #                  pf_sparse = w_sparsity, asparse = asparse, lambda = lambda_prec, intercept = FALSE)
    }
    if (lambda_prec_type == "min"){
    delta <- coef(mod, s = "lambda.min")[, 1]
    }else{
      delta <- coef(mod, s = "lambda.1se")[, 1]
    }
    
    res.model <- Y - cbind(1, as.matrix(Xmat)) %*% delta
    df_j <- sum(delta != 0)
    sigma <- sum(res.model^2) / (n - df_j)
    Delta_row <- rep(0, ncol(interM))
    Delta_row[indexes_to_select] <- as.vector(delta)
    No_sim_Delta_hat_row <- rep(0, ncol(interM))
    No_sim_Delta_hat_row[indexes_to_select] <- -sigma * as.vector(delta)
    list(Delta_row = Delta_row, No_sim_Delta_hat_row = No_sim_Delta_hat_row, sigma = sigma)
  }

  stopCluster(cl)

  for (i in 1:p) {
    Delta[i, ] <- results[,i]$Delta_row
    No_sim_Delta_hat[i, ] <- results[,i]$No_sim_Delta_hat_row
    Sigma_hat[i] <- results[,i]$sigma
  }

  Dic_Delta_hat <- list()
  fact_groups <- as.factor(groups)
  if (!is.null(covariates)) {
    U_names <- colnames(C)[-1]
  } else {
    U_names <- c()
  }
  U_names <- c("Intercept", "Prot", U_names)

  for (k in 2:length(levels(fact_groups))) {
    id <- levels(fact_groups)[k]
    key <- U_names[k]
    value <- No_sim_Delta_hat[, fact_groups == id]
    Temp <- abs(value)
    Temp_transpose <- abs(t(value))
    max_matrix <- pmax(Temp, Temp_transpose)
    sim_value <- sign(value + t(value)) * max_matrix
    Dic_Delta_hat[[key]] <- sim_value
  }
  
  return(
    list(
      deltas = Delta,
      Sigma_hat = Sigma_hat,
      Dic_Delta_hat = Dic_Delta_hat
    )
  )
}

GGReg_adaptive_cov_estimation <- function(Z0, # expression data (residuals on the mean), nrow: subjects; ncol: features.
                                 known_ppi = NULL, # previously known PPI 
                                 covariates = NULL, # covariates to regress on
                                 scr = TRUE, # optional screening to speed up
                                 gamma = NULL, # Person correlation threshold for the screening
                                 lambda_prec = NULL, # Penalization term in the Lasso
                                 lambda_prec_type = "1se", # or "min"
                                 weight = 1.1, # Multiplicative factor for the penalization term when prior knowledge is not available
                                 asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
                                 verbose = FALSE,
                                 eps = 1e-08) {
  

  n <- nrow(Z0)
  p <- ncol(Z0)
  Ip <- diag(rep(1, p))

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

  if (n < (q+1)*p) {
    warning("The sample size is too small! Network estimate may be unreliable!")
  }
  
  if (!is.null(known_ppi)) {
    one <- known_ppi
    one[one > 0] <- 1
    rownames(one) <- colnames(Z0)
    penaltyfactor <- 1 - known_ppi
  } else {
    one <- matrix(0, p, p)
    rownames(one) <- colnames(Z0)
    penaltyfactor <- matrix(1, p, p)
  }

  if (scr) {
    if (verbose) {
      cat("Computing correlation matrix \n")
    }
    protdata_matrix <- as.matrix(Z0)
    Scor <- cor(protdata_matrix)
    if (verbose) {
      cat("Done computing correlation matrix \n")
    }
    scr_index <- matrix(1, p, p)
    if (is.null(gamma)) {
      gamma <- quantile(abs(Scor), probs = 0.1)
    }
    scr_index[abs(Scor) <= gamma] <- 0
    diag(scr_index) <- 0
    colnames(scr_index) <- colnames(Z0)
    rownames(scr_index) <- colnames(Z0)
  }

  interM <- data.frame(Intercept = rep(1, n))
  temp_groups <- c(0, rep(1, ncol(Z0)))
  temp_w_group <- c(0)


  # Register parallel backend for the first loop
  numCores <- detectCores()
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
    cat("Comuputation of the design matrix \n ")
  }
  # Parallelize the first loop - Need to be tuned 
  interM_results <- foreach(j = 1:(q + 1), .combine = 'cbind') %dopar% {
    product <- Z0 * iU[, j]
    if (j != 1) {
      original_colnames <- colnames(product)
      iU_name <- colnames(iU)[j]
      temp_groups <- rep(j, length(original_colnames))
      temp_w_group <- 1
      new_colnames <- paste(iU_name, ":", original_colnames, sep = "")
      colnames(product) <- new_colnames
    }
    list(product = product, temp_groups  = temp_groups , temp_w_group=temp_w_group)
  }

  stopCluster(cl)
  groups <- c()
  w_group <- c()
  for(j in 1:(q + 1)){
    groups <- c(groups, interM_results[,j]$temp_groups)
    w_group <- c(w_group, interM_results[,j]$temp_w_group)
    interM <- cbind(interM, interM_results[,j]$product)
    }

 # FROM HERE ON
  Delta <- matrix(0, p, ncol(interM))
  colnames(Delta) <- colnames(interM)
  rownames(Delta) <- colnames(Z0)
  Sigma_hat <- rep(1, p)
  names(Sigma_hat) <- colnames(Z0)
  No_sim_Delta_hat <- matrix(0, p, ncol(interM))
  colnames(No_sim_Delta_hat) <- colnames(interM)
  rownames(No_sim_Delta_hat) <- colnames(Z0)

  # Register parallel backend for the second loop
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
        cat("Comuputation of adjacency matrix \n ")
      }

  results <- foreach(i = 1:p, .combine = 'cbind', .packages = c('glmnet', 'sparsegl')) %dopar% {
    set.seed(2024)
    Y <- matrix(Z0[, i], ncol = 1)
    protname <- colnames(Z0)[i]
    pattern <- paste0(":", protname, "\\b")
    indexes_to_select <- !grepl(pattern, colnames(interM))
    indexes_to_select[colnames(interM) == protname] <- FALSE
    indexes_to_select_penalty <- rep(TRUE, ncol(Z0))
    indexes_to_select_penalty[colnames(Z0) == protname] <- FALSE
    if (scr) {
      zero_col_indices <- which(scr_index[i, ] == 0)
      zero_col_names <- colnames(scr_index)[zero_col_indices]
      for (j in zero_col_names) {
        pattern <- paste0(":", j, "\\b")
        indexes_to_select[grepl(pattern, colnames(interM))] <- FALSE
        indexes_to_select[colnames(interM) == j] <- FALSE
        indexes_to_select_penalty[colnames(Z0) == j] <- FALSE
      }
    }
    indexes_to_select_reg <- indexes_to_select
    indexes_to_select_reg[1] <- FALSE
    Xmat <- interM[, indexes_to_select_reg]
    Xgroup <- groups[indexes_to_select_reg]
    
    pc <- prcomp(Xmat, center=TRUE, scale.=TRUE )
    summ <- summary(pc)
    df <- summ$importance
    idx <- which(df[3,] >= 0.80)[1]
    scores <- pc$x[,1:idx]
    scores <- cbind(Y,scores)
    colnames(scores)[1] <- "Y"
    scores <- as.data.frame(scores)
    lm_model <- lm(Y~0+.,data=scores)
    w_sparsity <- 1/(abs(pc$rotation[,1:idx] %*% lm_model$coefficients)^(0.2))
    w_sparsity[Xgroup == 1] <- w_sparsity[Xgroup == 1]*penaltyfactor[i, indexes_to_select_penalty]
    for(k in unique(Xgroup)[-1]){
      w_group[k] <- sqrt(sum(Xgroup == k))/(norm((pc$rotation[,1:idx] %*% lm_model$coefficients)[Xgroup == k], type="2"))
      }
  
    
    if (is.null(lambda_prec)) {
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, intercept = FALSE)
    } else {
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, lambda = lambda_prec, intercept = FALSE)
    }
    if (lambda_prec_type == "min"){
    delta <- coef(mod, s = "lambda.min")[, 1]
    }else{
      delta <- coef(mod, s = "lambda.1se")[, 1]
    }
    
    res.model <- Y - cbind(1, as.matrix(Xmat)) %*% delta
    df_j <- sum(delta != 0)
    sigma <- sum(res.model^2) / (n - df_j)
    Delta_row <- rep(0, ncol(interM))
    Delta_row[indexes_to_select] <- as.vector(delta)
    No_sim_Delta_hat_row <- rep(0, ncol(interM))
    No_sim_Delta_hat_row[indexes_to_select] <- -sigma * as.vector(delta)
    list(Delta_row = Delta_row, No_sim_Delta_hat_row = No_sim_Delta_hat_row, sigma = sigma)
  }

  stopCluster(cl)

  for (i in 1:p) {
    Delta[i, ] <- results[,i]$Delta_row
    No_sim_Delta_hat[i, ] <- results[,i]$No_sim_Delta_hat_row
    Sigma_hat[i] <- results[,i]$sigma
  }

  Dic_Delta_hat <- list()
  fact_groups <- as.factor(groups)
  if (!is.null(covariates)) {
    U_names <- colnames(C)[-1]
  } else {
    U_names <- c()
  }
  U_names <- c("Intercept", "Prot", U_names)

  for (k in 2:length(levels(fact_groups))) {
    id <- levels(fact_groups)[k]
    key <- U_names[k]
    value <- No_sim_Delta_hat[, fact_groups == id]
    Temp <- abs(value)
    Temp_transpose <- abs(t(value))
    max_matrix <- pmax(Temp, Temp_transpose)
    sim_value <- sign(value + t(value)) * max_matrix
    Dic_Delta_hat[[key]] <- sim_value
  }
  
  return(
    list(
      deltas = Delta,
      Sigma_hat = Sigma_hat,
      Dic_Delta_hat = Dic_Delta_hat
    )
  )
}



#### TO DO: PREDICTION WITH THIS MODEL

# PARALLELIZED
GGReg_full_estimation <-
  function (x, # expression data, nrow: subjects; ncol: features.
            known_ppi = NULL, # previously known PPI 
            covariates = NULL,       #covariates to regress on
            scr = TRUE, # optional screening to speed up
            gamma = NULL, #Person correlation threshold for the screening
            lambda_mean = NULL, #Penalization term in the Lasso for the mean
            lambda_mean_type = "1se", # or "min"
            lambda_prec = NULL,
            lambda_prec_type = "1se", # or "min"
            weight = 1.1, #Multiplicative factor for the penalization term when prior knowledgw is not available
            asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
            verbose = FALSE,
            eps = 1e-08) {
    
        if (verbose) {
        cat("Estimating the mean \n")
        }
        res_mean_reg <- GGReg_mean_estimation(x=x, covariates = covariates, lambda_mean=lambda_mean, lambda_mean_type=lambda_mean_type,  verbose = verbose, eps = eps)
        if (verbose) {
        cat("Done estimating the mean \n")
        }
        # Z0 = scale(res_mean_reg$z, center = TRUE, scale = TRUE) try and remove the scaling of the outcome
        Z0 <-res_mean_reg$z
        if (verbose) {
        cat("Estimating the precision matrix \n")
        }
        res_cov_reg <- GGReg_cov_estimation (Z0=Z0, known_ppi = known_ppi, covariates = covariates,
                                     scr = scr, gamma = gamma, lambda_prec=lambda_prec,lambda_prec_type=lambda_prec_type,
                                     weight = weight, asparse = asparse, 
                                     verbose = verbose, eps = eps)

        if (verbose) {
        cat("Done estimating the precision matrix \n")
        }
        return(
         list(
            z = res_mean_reg$z,
            Cov_effect = res_mean_reg$Cov_effect,
            deltas = res_cov_reg$Delta,
            Sigma_hat = res_cov_reg$Sigma_hat,
            Dic_Delta_hat = res_cov_reg$Dic_Delta_hat
         )
        )

    }


GGReg_full_estimation_adapt <-
  function (x, # expression data, nrow: subjects; ncol: features.
            known_ppi = NULL, # previously known PPI 
            covariates = NULL,       #covariates to regress on
            scr = TRUE, # optional screening to speed up
            gamma = NULL, #Person correlation threshold for the screening
            lambda_mean = NULL, #Penalization term in the Lasso for the mean
            lambda_mean_type = "1se", # or "min"
            lambda_prec = NULL,
            lambda_prec_type = "1se", # or "min"
            weight = 1.1, #Multiplicative factor for the penalization term when prior knowledgw is not available
            asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
            verbose = FALSE,
            eps = 1e-08) {
    
        if (verbose) {
        cat("Estimating the mean \n")
        }
        res_mean_reg <- GGReg_mean_estimation(x=x, covariates = covariates, lambda_mean=lambda_mean, lambda_mean_type=lambda_mean_type,  verbose = verbose, eps = eps)
        if (verbose) {
        cat("Done estimating the mean \n")
        }
        # Z0 = scale(res_mean_reg$z, center = TRUE, scale = TRUE) try and remove the scaling of the outcome
        Z0 <-res_mean_reg$z
        if (verbose) {
        cat("Estimating the precision matrix \n")
        }
        res_cov_reg <- GGReg_adaptive_cov_estimation (Z0=Z0, known_ppi = known_ppi, covariates = covariates,
                                     scr = scr, gamma = gamma, lambda_prec=lambda_prec,lambda_prec_type=lambda_prec_type,
                                     weight = weight, asparse = asparse, 
                                     verbose = verbose, eps = eps)

        if (verbose) {
        cat("Done estimating the precision matrix \n")
        }
        return(
         list(
            z = res_mean_reg$z,
            Cov_effect = res_mean_reg$Cov_effect,
            deltas = res_cov_reg$Delta,
            Sigma_hat = res_cov_reg$Sigma_hat,
            Dic_Delta_hat = res_cov_reg$Dic_Delta_hat
         )
        )

    }

GGReg_full_estimation_stab_adapt <-
  function (x, # expression data, nrow: subjects; ncol: features.
            known_ppi = NULL, # previously known PPI 
            covariates = NULL,       #covariates to regress on
            scr = TRUE, # optional screening to speed up
            gamma = NULL, #Person correlation threshold for the screening
            lambda_mean = NULL, #Penalization term in the Lasso for the mean
            lambda_mean_type = "1se", # or "min"
            lambda_prec = NULL,
            lambda_prec_type = "1se", # or "min"
            weight = 1.1, #Multiplicative factor for the penalization term when prior knowledgw is not available
            asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
            verbose = FALSE,
            eps = 1e-08) {
    
        if (verbose) {
        cat("Estimating the mean \n")
        }
        res_mean_reg <- GGReg_stable_mean_estimation(x=x, covariates = covariates, lambda_mean=lambda_mean, lambda_mean_type=lambda_mean_type,  verbose = verbose, eps = eps)
        if (verbose) {
        cat("Done estimating the mean \n")
        }
        # Z0 = scale(res_mean_reg$z, center = TRUE, scale = TRUE) try and remove the scaling of the outcome
        Z0 <-res_mean_reg$z
        if (verbose) {
        cat("Estimating the precision matrix \n")
        }
        res_cov_reg <- GGReg_adaptive_cov_estimation (Z0=Z0, known_ppi = known_ppi, covariates = covariates,
                                     scr = scr, gamma = gamma, lambda_prec=lambda_prec,lambda_prec_type=lambda_prec_type,
                                     weight = weight, asparse = asparse, 
                                     verbose = verbose, eps = eps)

        if (verbose) {
        cat("Done estimating the precision matrix \n")
        }
        return(
         list(
            z = res_mean_reg$z,
            Cov_effect = res_mean_reg$Cov_effect,
            deltas = res_cov_reg$Delta,
            Sigma_hat = res_cov_reg$Sigma_hat,
            Dic_Delta_hat = res_cov_reg$Dic_Delta_hat
         )
        )

    }






GGReg_cov_estimation_diff_scr <- function(Z0, # expression data (residuals on the mean), nrow: subjects; ncol: features.
                                 known_ppi = NULL, # previously known PPI 
                                 covariates = NULL, # covariates to regress on
                                 scr = TRUE, # optional screening to speed up
                                 gamma = NULL, # Person correlation threshold for the screening
                                 lambda_prec = NULL, # Penalization term in the Lasso
                                 lambda_prec_type = "1se", # or "min"
                                 weight = 1.1, # Multiplicative factor for the penalization term when prior knowledge is not available
                                 asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
                                 verbose = FALSE,
                                 eps = 1e-08) {
  

  n <- nrow(Z0)
  p <- ncol(Z0)
  Ip <- diag(rep(1, p))

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

  if (n < (q+1)*p) {
    warning("The sample size is too small! Network estimate may be unreliable!")
  }
  
  if (!is.null(known_ppi)) {
    one <- known_ppi
    one[one > 0] <- 1
    rownames(one) <- colnames(Z0)
    penaltyfactor <- 1 - known_ppi
  } else {
    one <- matrix(0, p, p)
    rownames(one) <- colnames(Z0)
    penaltyfactor <- matrix(1, p, p)
  }

  if (scr) {
    if (verbose) {
      cat("Computing correlation matrix \n")
    }
    protdata_matrix <- as.matrix(Z0)
    Scor <- cor(protdata_matrix)
    if (verbose) {
      cat("Done computing correlation matrix \n")
    }
    scr_index <- matrix(1, p, p)
    if (is.null(gamma)) {
      gamma <- quantile(abs(Scor), probs = 0.1)
    }
    scr_index[abs(Scor) <= gamma] <- 0
    diag(scr_index) <- 0
    colnames(scr_index) <- colnames(Z0)
    rownames(scr_index) <- colnames(Z0)
  }

  interM <- data.frame(Intercept = rep(1, n))
  temp_groups <- c(0, rep(1, ncol(Z0)))
  temp_w_group <- c(0)


  # Register parallel backend for the first loop
  numCores <- detectCores()
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
    cat("Comuputation of the design matrix \n ")
  }
  # Parallelize the first loop - Need to be tuned 
  interM_results <- foreach(j = 1:(q + 1), .combine = 'cbind') %dopar% {
    product <- Z0 * iU[, j]
    if (j != 1) {
      original_colnames <- colnames(product)
      iU_name <- colnames(iU)[j]
      temp_groups <- rep(j, length(original_colnames))
      temp_w_group <- 1
      new_colnames <- paste(iU_name, ":", original_colnames, sep = "")
      colnames(product) <- new_colnames
    }
    list(product = product, temp_groups  = temp_groups , temp_w_group=temp_w_group)
  }

  stopCluster(cl)
  groups <- c()
  w_group <- c()
  for(j in 1:(q + 1)){
    groups <- c(groups, interM_results[,j]$temp_groups)
    w_group <- c(w_group, interM_results[,j]$temp_w_group)
    interM <- cbind(interM, interM_results[,j]$product)
    }

 # FROM HERE ON
  Delta <- matrix(0, p, ncol(interM))
  colnames(Delta) <- colnames(interM)
  rownames(Delta) <- colnames(Z0)
  Sigma_hat <- rep(1, p)
  names(Sigma_hat) <- colnames(Z0)
  No_sim_Delta_hat <- matrix(0, p, ncol(interM))
  colnames(No_sim_Delta_hat) <- colnames(interM)
  rownames(No_sim_Delta_hat) <- colnames(Z0)

  # Register parallel backend for the second loop
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
        cat("Comuputation of adjacency matrix \n ")
      }

  results <- foreach(i = 1:p, .combine = 'cbind', .packages = c('glmnet', 'sparsegl')) %dopar% {
    Y <- matrix(Z0[, i], ncol = 1)
    protname <- colnames(Z0)[i]
    pattern <- paste0(":", protname, "\\b")
    indexes_to_select <- !grepl(pattern, colnames(interM))
    indexes_to_select[colnames(interM) == protname] <- FALSE
    indexes_to_select_penalty <- rep(TRUE, ncol(Z0))
    indexes_to_select_penalty[colnames(Z0) == protname] <- FALSE
    if (scr) {
      zero_col_indices <- which(scr_index[i, ] == 0)
      zero_col_names <- colnames(scr_index)[zero_col_indices]
      for (j in zero_col_names) {
        indexes_to_select[colnames(interM) == j] <- FALSE
        indexes_to_select_penalty[colnames(Z0) == j] <- FALSE
      }
    }
    indexes_to_select_reg <- indexes_to_select
    indexes_to_select_reg[1] <- FALSE
    Xmat <- interM[, indexes_to_select_reg]
    Xgroup <- groups[indexes_to_select_reg]
    
    w_sparsity <- rep(weight, ncol(Xmat))
    w_sparsity[Xgroup == 1] <- penaltyfactor[i, indexes_to_select_penalty]
    w_group[which(w_group==1)] <- sqrt(sum(Xgroup == 1))

    
    if (is.null(lambda_prec)) {
      set.seed(2024)
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, intercept = FALSE)
      #mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup,
      #                  pf_sparse = w_sparsity, asparse = asparse, intercept = FALSE)
    } else {
      set.seed(2024)
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, lambda = lambda_prec, intercept = FALSE)
      #mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup,
      #                  pf_sparse = w_sparsity, asparse = asparse, lambda = lambda_prec, intercept = FALSE)
    }
    if (lambda_prec_type == "min"){
    delta <- coef(mod, s = "lambda.min")[, 1]
    }else{
      delta <- coef(mod, s = "lambda.1se")[, 1]
    }
    
    res.model <- Y - cbind(1, as.matrix(Xmat)) %*% delta
    df_j <- sum(delta != 0)
    sigma <- sum(res.model^2) / (n - df_j)
    Delta_row <- rep(0, ncol(interM))
    Delta_row[indexes_to_select] <- as.vector(delta)
    No_sim_Delta_hat_row <- rep(0, ncol(interM))
    No_sim_Delta_hat_row[indexes_to_select] <- -sigma * as.vector(delta)
    list(Delta_row = Delta_row, No_sim_Delta_hat_row = No_sim_Delta_hat_row, sigma = sigma)
  }

  stopCluster(cl)

  for (i in 1:p) {
    Delta[i, ] <- results[,i]$Delta_row
    No_sim_Delta_hat[i, ] <- results[,i]$No_sim_Delta_hat_row
    Sigma_hat[i] <- results[,i]$sigma
  }

  Dic_Delta_hat <- list()
  fact_groups <- as.factor(groups)
  if (!is.null(covariates)) {
    U_names <- colnames(C)[-1]
  } else {
    U_names <- c()
  }
  U_names <- c("Intercept", "Prot", U_names)

  for (k in 2:length(levels(fact_groups))) {
    id <- levels(fact_groups)[k]
    key <- U_names[k]
    value <- No_sim_Delta_hat[, fact_groups == id]
    Temp <- abs(value)
    Temp_transpose <- abs(t(value))
    max_matrix <- pmax(Temp, Temp_transpose)
    sim_value <- sign(value + t(value)) * max_matrix
    Dic_Delta_hat[[key]] <- sim_value
  }
  
  return(
    list(
      deltas = Delta,
      Sigma_hat = Sigma_hat,
      Dic_Delta_hat = Dic_Delta_hat
    )
  )
}

GGReg_adaptive_cov_estimation_diff_scr <- function(Z0, # expression data (residuals on the mean), nrow: subjects; ncol: features.
                                 known_ppi = NULL, # previously known PPI 
                                 covariates = NULL, # covariates to regress on
                                 scr = TRUE, # optional screening to speed up
                                 gamma = NULL, # Person correlation threshold for the screening
                                 lambda_prec = NULL, # Penalization term in the Lasso
                                 lambda_prec_type = "1se", # or "min"
                                 weight = 1.1, # Multiplicative factor for the penalization term when prior knowledge is not available
                                 asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
                                 verbose = FALSE,
                                 eps = 1e-08) {
  

  n <- nrow(Z0)
  p <- ncol(Z0)
  Ip <- diag(rep(1, p))

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

  if (n < (q+1)*p) {
    warning("The sample size is too small! Network estimate may be unreliable!")
  }
  
  if (!is.null(known_ppi)) {
    one <- known_ppi
    one[one > 0] <- 1
    rownames(one) <- colnames(Z0)
    penaltyfactor <- 1 - known_ppi
  } else {
    one <- matrix(0, p, p)
    rownames(one) <- colnames(Z0)
    penaltyfactor <- matrix(1, p, p)
  }

  if (scr) {
    if (verbose) {
      cat("Computing correlation matrix \n")
    }
    protdata_matrix <- as.matrix(Z0)
    Scor <- cor(protdata_matrix)
    if (verbose) {
      cat("Done computing correlation matrix \n")
    }
    scr_index <- matrix(1, p, p)
    if (is.null(gamma)) {
      gamma <- quantile(abs(Scor), probs = 0.1)
    }
    scr_index[abs(Scor) <= gamma] <- 0
    diag(scr_index) <- 0
    colnames(scr_index) <- colnames(Z0)
    rownames(scr_index) <- colnames(Z0)
  }

  interM <- data.frame(Intercept = rep(1, n))
  temp_groups <- c(0, rep(1, ncol(Z0)))
  temp_w_group <- c(0)


  # Register parallel backend for the first loop
  numCores <- detectCores()
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
    cat("Comuputation of the design matrix \n ")
  }
  # Parallelize the first loop - Need to be tuned 
  interM_results <- foreach(j = 1:(q + 1), .combine = 'cbind') %dopar% {
    product <- Z0 * iU[, j]
    if (j != 1) {
      original_colnames <- colnames(product)
      iU_name <- colnames(iU)[j]
      temp_groups <- rep(j, length(original_colnames))
      temp_w_group <- 1
      new_colnames <- paste(iU_name, ":", original_colnames, sep = "")
      colnames(product) <- new_colnames
    }
    list(product = product, temp_groups  = temp_groups , temp_w_group=temp_w_group)
  }

  stopCluster(cl)
  groups <- c()
  w_group <- c()
  for(j in 1:(q + 1)){
    groups <- c(groups, interM_results[,j]$temp_groups)
    w_group <- c(w_group, interM_results[,j]$temp_w_group)
    interM <- cbind(interM, interM_results[,j]$product)
    }

 # FROM HERE ON
  Delta <- matrix(0, p, ncol(interM))
  colnames(Delta) <- colnames(interM)
  rownames(Delta) <- colnames(Z0)
  Sigma_hat <- rep(1, p)
  names(Sigma_hat) <- colnames(Z0)
  No_sim_Delta_hat <- matrix(0, p, ncol(interM))
  colnames(No_sim_Delta_hat) <- colnames(interM)
  rownames(No_sim_Delta_hat) <- colnames(Z0)

  # Register parallel backend for the second loop
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  if (verbose) {
        cat("Comuputation of adjacency matrix \n ")
      }

  results <- foreach(i = 1:p, .combine = 'cbind', .packages = c('glmnet', 'sparsegl')) %dopar% {
    set.seed(2024)
    Y <- matrix(Z0[, i], ncol = 1)
    protname <- colnames(Z0)[i]
    pattern <- paste0(":", protname, "\\b")
    indexes_to_select <- !grepl(pattern, colnames(interM))
    indexes_to_select[colnames(interM) == protname] <- FALSE
    indexes_to_select_penalty <- rep(TRUE, ncol(Z0))
    indexes_to_select_penalty[colnames(Z0) == protname] <- FALSE
    if (scr) {
      zero_col_indices <- which(scr_index[i, ] == 0)
      zero_col_names <- colnames(scr_index)[zero_col_indices]
      for (j in zero_col_names) {
        indexes_to_select[colnames(interM) == j] <- FALSE
        indexes_to_select_penalty[colnames(Z0) == j] <- FALSE
      }
    }
    indexes_to_select_reg <- indexes_to_select
    indexes_to_select_reg[1] <- FALSE
    Xmat <- interM[, indexes_to_select_reg]
    Xgroup <- groups[indexes_to_select_reg]
    
    pc <- prcomp(Xmat, center=TRUE, scale.=TRUE )
    summ <- summary(pc)
    df <- summ$importance
    idx <- which(df[3,] >= 0.80)[1]
    scores <- pc$x[,1:idx]
    scores <- cbind(Y,scores)
    colnames(scores)[1] <- "Y"
    scores <- as.data.frame(scores)
    lm_model <- lm(Y~0+.,data=scores)
    w_sparsity <- 1/(abs(pc$rotation[,1:idx] %*% lm_model$coefficients)^(0.2))
    w_sparsity[Xgroup == 1] <- w_sparsity[Xgroup == 1]*penaltyfactor[i, indexes_to_select_penalty]
    for(k in unique(Xgroup)[-1]){
      w_group[k] <- sqrt(sum(Xgroup == k))/(norm((pc$rotation[,1:idx] %*% lm_model$coefficients)[Xgroup == k], type="2"))
      }
  
    
    if (is.null(lambda_prec)) {
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, intercept = FALSE)
    } else {
      mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, pf_group = w_group,
                        pf_sparse = w_sparsity, asparse = asparse, lambda = lambda_prec, intercept = FALSE)
    }
    if (lambda_prec_type == "min"){
    delta <- coef(mod, s = "lambda.min")[, 1]
    }else{
      delta <- coef(mod, s = "lambda.1se")[, 1]
    }
    
    res.model <- Y - cbind(1, as.matrix(Xmat)) %*% delta
    df_j <- sum(delta != 0)
    sigma <- sum(res.model^2) / (n - df_j)
    Delta_row <- rep(0, ncol(interM))
    Delta_row[indexes_to_select] <- as.vector(delta)
    No_sim_Delta_hat_row <- rep(0, ncol(interM))
    No_sim_Delta_hat_row[indexes_to_select] <- -sigma * as.vector(delta)
    list(Delta_row = Delta_row, No_sim_Delta_hat_row = No_sim_Delta_hat_row, sigma = sigma)
  }

  stopCluster(cl)

  for (i in 1:p) {
    Delta[i, ] <- results[,i]$Delta_row
    No_sim_Delta_hat[i, ] <- results[,i]$No_sim_Delta_hat_row
    Sigma_hat[i] <- results[,i]$sigma
  }

  Dic_Delta_hat <- list()
  fact_groups <- as.factor(groups)
  if (!is.null(covariates)) {
    U_names <- colnames(C)[-1]
  } else {
    U_names <- c()
  }
  U_names <- c("Intercept", "Prot", U_names)

  for (k in 2:length(levels(fact_groups))) {
    id <- levels(fact_groups)[k]
    key <- U_names[k]
    value <- No_sim_Delta_hat[, fact_groups == id]
    Temp <- abs(value)
    Temp_transpose <- abs(t(value))
    max_matrix <- pmax(Temp, Temp_transpose)
    sim_value <- sign(value + t(value)) * max_matrix
    Dic_Delta_hat[[key]] <- sim_value
  }
  
  return(
    list(
      deltas = Delta,
      Sigma_hat = Sigma_hat,
      Dic_Delta_hat = Dic_Delta_hat
    )
  )
}

GGReg_full_estimation_diff_scr <-
  function (x, # expression data, nrow: subjects; ncol: features.
            known_ppi = NULL, # previously known PPI 
            covariates = NULL,       #covariates to regress on
            scr = TRUE, # optional screening to speed up
            gamma = NULL, #Person correlation threshold for the screening
            lambda_mean = NULL, #Penalization term in the Lasso for the mean
            lambda_mean_type = "1se", # or "min"
            lambda_prec = NULL,
            lambda_prec_type = "1se", # or "min"
            weight = 1.1, #Multiplicative factor for the penalization term when prior knowledgw is not available
            asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
            verbose = FALSE,
            eps = 1e-08) {
    
        if (verbose) {
        cat("Estimating the mean \n")
        }
        res_mean_reg <- GGReg_mean_estimation(x=x, covariates = covariates, lambda_mean=lambda_mean, lambda_mean_type=lambda_mean_type,  verbose = verbose, eps = eps)
        if (verbose) {
        cat("Done estimating the mean \n")
        }
        # Z0 = scale(res_mean_reg$z, center = TRUE, scale = TRUE) try and remove the scaling of the outcome
        Z0 <-res_mean_reg$z
        if (verbose) {
        cat("Estimating the precision matrix \n")
        }
        res_cov_reg <- GGReg_cov_estimation_diff_scr (Z0=Z0, known_ppi = known_ppi, covariates = covariates,
                                     scr = scr, gamma = gamma, lambda_prec=lambda_prec,lambda_prec_type=lambda_prec_type,
                                     weight = weight, asparse = asparse, 
                                     verbose = verbose, eps = eps)

        if (verbose) {
        cat("Done estimating the precision matrix \n")
        }
        return(
         list(
            z = res_mean_reg$z,
            Cov_effect = res_mean_reg$Cov_effect,
            deltas = res_cov_reg$Delta,
            Sigma_hat = res_cov_reg$Sigma_hat,
            Dic_Delta_hat = res_cov_reg$Dic_Delta_hat
         )
        )

    }



GGReg_full_estimation_stab_adapt_diff_scr <-
  function (x, # expression data, nrow: subjects; ncol: features.
            known_ppi = NULL, # previously known PPI 
            covariates = NULL,       #covariates to regress on
            scr = TRUE, # optional screening to speed up
            gamma = NULL, #Person correlation threshold for the screening
            lambda_mean = NULL, #Penalization term in the Lasso for the mean
            lambda_mean_type = "1se", # or "min"
            lambda_prec = NULL,
            lambda_prec_type = "1se", # or "min"
            weight = 1.1, #Multiplicative factor for the penalization term when prior knowledgw is not available
            asparse = 0.75, # The relative weight to put on the `1-norm in sparse group lasso
            verbose = FALSE,
            eps = 1e-08) {
    
        if (verbose) {
        cat("Estimating the mean \n")
        }
        res_mean_reg <- GGReg_stable_mean_estimation(x=x, covariates = covariates, lambda_mean=lambda_mean, lambda_mean_type=lambda_mean_type,  verbose = verbose, eps = eps)
        if (verbose) {
        cat("Done estimating the mean \n")
        }
        # Z0 = scale(res_mean_reg$z, center = TRUE, scale = TRUE) try and remove the scaling of the outcome
        Z0 <-res_mean_reg$z
        if (verbose) {
        cat("Estimating the precision matrix \n")
        }
        res_cov_reg <- GGReg_adaptive_cov_estimation_diff_scr (Z0=Z0, known_ppi = known_ppi, covariates = covariates,
                                     scr = scr, gamma = gamma, lambda_prec=lambda_prec,lambda_prec_type=lambda_prec_type,
                                     weight = weight, asparse = asparse, 
                                     verbose = verbose, eps = eps)

        if (verbose) {
        cat("Done estimating the precision matrix \n")
        }
        return(
         list(
            z = res_mean_reg$z,
            Cov_effect = res_mean_reg$Cov_effect,
            deltas = res_cov_reg$Delta,
            Sigma_hat = res_cov_reg$Sigma_hat,
            Dic_Delta_hat = res_cov_reg$Dic_Delta_hat
         )
        )

    }