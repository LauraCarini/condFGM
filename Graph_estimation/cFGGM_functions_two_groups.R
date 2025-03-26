library(parallel)
library(doParallel)
library(foreach)

# To be adapted to ours
func.path.int <- "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Previous_litt/Boxin_Zhao_FGM_Neighboorhood/Functions"
source(paste(func.path.int,"ADMM.new.R", sep="/"))

#############################
### FUNCTIONS TO PERFROM TWO GROUPS COMPARISON
#############################

############################# SAME AS THEIRS lambda.sup
lambda.sup <- function(A.X, A.Y){
  n <- nrow(A.X)
  M <- ncol(A.Y)
  p <- ncol(A.X)/M
  candidates <- rep(0,p)
  for(j in 1:p){
    A.X.j <- A.X[,(j-1)*M + (1:M)]
    candidates[j] <- norm(t(A.X.j) %*% A.Y,"F") / n
  }
  return(max(candidates))
}

############################# UPDATE lambda.sup.global
lambda.sup.global.two.groups <- function(design.matrix, M, p){
  
  n <- nrow(design.matrix)
  
  Y <- matrix(NA, n, M)
  X <- matrix(NA, n, 2*(p-1)*M)
  l.max <- rep(0,p)
  for(j in 1:p){
    jth.range.y <- (j-1)*M+(1:M)
    jth.range.x <- c(jth.range.y,(j+p-1)*M+(1:M))
    for(i in 1:n){
      Y[i,] <- as.matrix(design.matrix[i, jth.range.y])
      X[i,] <- as.matrix(design.matrix[i, -jth.range.x])
    }
    l.max[j] <- lambda.sup(X, Y)
  }
  lambda.max <- max(l.max)
  return(lambda.max)
}


############################# UPDATE grp.soft.thres
grp.soft.thres.two.groups <- function(V.Grp, d, lambda, rho){
  # V.Grp: a group of p matrices (Vk), each dim M * M
  # d: vector, length pM.
  # lambda: penalty param of group lasso
  # rho: penalty param of augmented Lagrangian
  M <- nrow(V.Grp[[1]])
  n_blocks <- length(V.Grp)
  D.Grp <- list()
  P.Grp <- list() # output group of p matrices (Pk), each dim M * n
  
  for(k in 1:n_blocks){
    D.Grp[[k]] <- diag(d[(k-1)*M+(1:M)])
    left.plus <- matrix(0, M, M)
    norm.F <- norm(D.Grp[[k]] %*% V.Grp[[k]], "F")
    for(l in 1:M){
      left.plus[l,l] <- max(0, 1 - lambda / rho * D.Grp[[k]][l,l]^2 / norm.F )
    }
    
    P.Grp[[k]] <- left.plus %*% V.Grp[[k]]
  }
  
  return(P.Grp)
}

############################# UPDATE objective.function
objective.function.two.groups <- function(A.X, A.Y, d, B, lambda){
  # A.X: matrix, n * 2(p-1)M
  # A.Y: matrix, n * M
  # d: vector, length 2(p-1)M
  # B: matrix, 2(p-1)M * M
  # Q: matrix, 2(p-1)M * M
  # lambda: penalty parameter
  
  n <- nrow(A.X)
  M <- ncol(A.Y)
  n_blocks <- ncol(A.X) / M
  
  D <- diag(d)
  term.1 <- norm(A.X %*% D %*% B - A.Y, "F")^2 / (2*n)
  term.2 <- 0
  for(k in 1:n_blocks){
    B.k <- B[(k-1)*M + (1:M), ]
    Dkk <- diag(d[(k-1)*M + (1:M)])
    term.2 <- term.2 + lambda * norm(Dkk %*% B.k ,"F")
  }
  obj <- term.1 + term.2
  return(obj)
}

############################# UPDATE ADMM.grplasso -> miss a parantesis or one more

ADMM.grplasso.two.groups <- function(A.X, A.Y, d, lambda,
                                     rho.init=1, maxiter=2000, mu=10, tau=2,
                                     P.in, Q.in, U.in, tol.rel=1e-4, tol.abs=1e-4){
  
  # Input:
  #   A.X, n*2(p-1)M matrix.
  #   A.Y, n*M matrix
  #   d, a 2(p-1)M long vector
  # 1.1. Reading in the parameters
  n <- nrow(A.Y)
  M <- ncol(A.Y)
  n_blocks <- ncol(A.X) / M # corresponds to 2(p-1) in theory
  
  # 1.2. Recover matrices used in computation
  D <- diag(d)
  AX.D <- A.X %*% D
  DXXD.n <- (t(AX.D) %*% AX.D) / n
  DXY.n <- (t(AX.D) %*% A.Y)/n
  
  # 1.3. Setting initial values for P,Q,U
  P.old <- P.in
  Q.old <- Q.in
  U.old <- U.in
  
  # 1.4. Initialize output quantities
  prim.resid.vec <- numeric(0)
  dual.resid.vec <- numeric(0)
  rho.path <- rho.init
  tol.prim.path <- numeric(0)
  tol.dual.path <- numeric(0)
  obj.func.path <- numeric(0)
  
  # 1.5. Cached quantities -- initial
  # 1.5.1. Initial Vk's
  V <- Q.old - U.old # pM * M
  V.Grp <- list()
  for(k in 1:n_blocks){
    V.Grp[[k]] <- V[(k-1)*M + (1:M), ] # each Vk dim M * M
  }
  
  # 1.5.2. indicator of rho's change
  rho <- rho.init
  rho.changed <- TRUE
  ##############################
  ## 2. Main part: Iterations ##
  ##############################
  for(iter in 1:maxiter){
    
    #############################
    ## 7.1. Updating Variables ##
    #############################
    # 7.1.1. Update P with group soft thresholding
    P.new <- matrix(NA, nrow=n_blocks*M, ncol=M)
    P.Grp <- grp.soft.thres.two.groups(V.Grp, d, lambda, rho) # Group of p, each M * M
    for(k in 1:n_blocks){
      P.new[(k-1)*M + (1:M), ] <- P.Grp[[k]]
    }
    
    # 7.1.2. Update Q
    if(rho.changed){
      left.inv <- tryCatch(
        solve(DXXD.n + rho * diag(nrow(DXXD.n))),
        error = function(e) {
          stop("Matrix inversion error: Check dimensions of DXXD.n and rho term")
        }
      )
    }
    Q.new <- left.inv %*% (DXY.n + rho * P.new + rho * U.old)
    
    # 7.1.3. Update U
    U.new <- U.old + P.new - Q.new
    
    
    #####################################
    ## 7.2. Checking Stopping Criteria ##
    #####################################
    
    # 7.2.1. Calculate primal and dual residuals
    prim.resid <- norm( (P.new - Q.new) , "F")
    dual.resid <- norm( rho * (Q.new - Q.old) , "F")
    prim.resid.vec <- c(prim.resid.vec,prim.resid)
    dual.resid.vec <- c(dual.resid.vec,dual.resid)
    
    # 7.2.2. Calculate tolerance
    tol.prim <- tol.abs * sqrt(n_blocks*M*M) + tol.rel * max(norm(P.new, "F"), 
                                                             norm(Q.new, "F"))
    tol.dual <- tol.abs * sqrt(n_blocks*M*M) + tol.rel * norm(U.new, "F")
    # sqrt(n_blocksM^2) is because the Frob norms are in R^{n_blocksM * M}.

    # 7.2.3. Compare residuals with tolerance
    if(prim.resid < tol.prim && dual.resid < tol.dual){
      break
    }
    
    #####################################
    ## 7.3. Prepare for next iteration ##
    #####################################
    
    # 7.3.1. Update V
    V <- Q.new - U.new # pM * M
    V.Grp <- list()
    for(k in 1:n_blocks){
      V.Grp[[k]] <- V[((k-1)*M + 1):(k*M), ] # each Vk dim M * n
    }
    
    # 7.3.2. Update rho for next round
    if(prim.resid > mu * dual.resid){
      rho <- tau * rho
      rho.changed <- TRUE
    }else if(dual.resid > mu * prim.resid){
      rho <- rho / tau
      rho.changed <- TRUE
    }else{
      rho.changed <- FALSE
    }
    rho.path <- c(rho.path, rho)
    
    # 7.3.3. old <- new, end of one iterate
    P.old <- P.new
    Q.old <- Q.new
    U.old <- U.new
    
    obj <- tryCatch(
      objective.function.two.groups(A.X=A.X, A.Y=A.Y, d=d, B=P.new, lambda=lambda),
      error = function(e) {
        stop("Error in objective function calculation: Check objective.function.two.groups()")
      }
    )
    
    obj.func.path <- c(obj.func.path, obj)
  }
  
  #message("ADMM converged after ", iter, " iterations.")
  
  result <- list(P=P.new, Q=Q.new, U=U.new,
                 prim.res=prim.resid.vec, dual.res=dual.resid.vec,
                 rho.path=rho.path,
                 tol.prim.path=tol.prim.path, tol.dual.path=tol.dual.path,
                 obj.func.path=obj.func.path)
  return(result)}




# ############################# NOTATION ON THE OUTCOMES
# #The outcomes I'd like
# Delta_scores <- matrix(0, Mp, ncol(interM))
# colnames(Delta_scores) <- colnames(interM)
# rownames(Delta_scores) <- colnames(scores)
# No_sim_Delta_hat <- matrix(0, Mp, ncol(interM))
# colnames(No_sim_Delta_hat) <- colnames(interM)
# rownames(No_sim_Delta_hat) <- colnames(scores)
# 
# #The outcomes they have
# time.start.gX <- proc.time()
# 
# G.list.gX <- list()
# TPR.and.gX <- rep(0,L)
# FPR.and.gX <- rep(0,L)
# TPR.or.gX <- rep(0,L)
# FPR.or.gX <- rep(0,L)



FGGReg_diff_two_groups <- function(scores, # functional score on a defined basis, nrow: subjects; ncol: functions*n_basis.
                                   n_basis = 1, #Number of bases considered 
                                   covariates  = NULL, #Additional covariates to regress on
                                   lambda = 0 , # Penalization term in the Lasso 
                                   thres.ctrl = 0, # recognition threshold epsilon_n = thres.ctrl * lambda_n,
                                   verbose = FALSE,
                                   tol.abs =1e-4 ,
                                   tol.rel = 1e-4,
                                   eps = 1e-08) {
  
  time.start.gX <- proc.time()
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
  if (n < (q+1)*Mp) {
    warning("The sample size is too small! Network estimate may be unreliable!")
  }
  
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
  
  ## COMPUTE THE MULTIPLE REGRESSION
  
  # Q: A cosa servono d.array e  d ? Quali sono le dimensioni corrette?
  # Per il momento li ho sostuiti com dei vettori di uni della dimensione coretta (?)
  d.array <- matrix(1, nrow=p, ncol=(p-1)*M*(q+1))
  d_out <- list()
  norm.adj <- rep(NA,p)
  for(k in 1:p){
    d_out[[k]] <- d.array[k,]
    norm.adj[k] <- norm(d_out[[k]],"2")
  }
  
  # default warm start PQU
  P.def <- matrix(0, 2*(p-1)*M, M)
  Q.def <- matrix(0.1, 2*(p-1)*M, M)
  U.def <- matrix(0.01, 2*(p-1)*M, M)
  
  # Repeat the following loop for each node
  V.array <- foreach(j = 1:p, .combine="rbind", .packages="fda", .export=c('ADMM.grplasso.two.groups','lambda.sup',
                                                                           'lambda.sup.global.two.groups', 'grp.soft.thres.two.groups',
                                                                           'objective.function.two.groups')) %dopar% {
                                                                             message("Processing node ", j)
                                                                             jth.range.y <- (j-1)*M+(1:M)
                                                                             A.Y <- as.matrix(interM[, jth.range.y])
                                                                             jth.range.x <- c(jth.range.y,(j+p-1)*M+(1:M))
                                                                             A.X <- as.matrix(interM[, -jth.range.x])
                                                                             groups <- temp_groups[-jth.range.x]
                                                                             
                                                                             P <- P.def; Q <- Q.def; U <- U.def
                                                                             
                                                                             
                                                                             V.j <- matrix(NA, 1, 2*p)
                                                                             
                                                                             
                                                                             grp.lasso.result <- tryCatch({
                                                                                 ADMM.grplasso.two.groups(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
                                                                                                          lambda=lambda, rho.init=1,
                                                                                                          P.in=P, Q.in=Q, U.in=U,
                                                                                                          tol.abs = tol.abs, tol.rel = tol.rel,
                                                                                                          maxiter = 400)
                                                                               }, error = function(e) {
                                                                                 message("Error in ADMM.grplasso.two.groups: ", e$message)
                                                                                 return(NULL)
                                                                               })
                                                                               if (is.null(grp.lasso.result)) next
                                                                               P <- grp.lasso.result$P
                                                                               Q <- grp.lasso.result$Q
                                                                               U <- grp.lasso.result$U
                                                                               
                                                                               n_blocks <- nrow(P)/M
                                                                               # Process the estimated P into a neighborhood selection vector
                                                                               P.frob <- list()
                                                                               for(k in 1:n_blocks){
                                                                                 key <- unique(groups[(k-1)*M + (1:M)])
                                                                                 if(length(key) > 1){
                                                                                   message("Error in group definition for block ", k)
                                                                                   next
                                                                                 } else {
                                                                                   P.frob[[key]] <- norm(P[(k-1)*M + (1:M), ], "F")
                                                                                 }
                                                                               }
                                                                               
                                                                               threshold <- lambda * thres.ctrl
                                                                               
                                                                               # Neighbor recognization
                                                                               V.j <- rep(0, length(P.frob))
                                                                               for(n_block in 1:length(P.frob)){
                                                                                 if(!(is.null(P.frob[[n_block]]))){
                                                                                   if(P.frob[[n_block]] > threshold){
                                                                                     V.j[n_block] <- 1} 
                                                                                 }
                                                                               }
                                                                               if(j == p){V.j <- c(V.j,0)}
                                                                             
                                                                             V.j
                                                                           }
  
  time.end.gX <- proc.time()
  runtime.gX <- (time.end.gX - time.start.gX)[3]
  cat("Computational time of: ")
  cat(runtime.gX)
  return(V.array)
}
 
FGGReg_diff_two_groups_SCV <- function(scores, # functional score on a defined basis, nrow: subjects; ncol: functions*n_basis.
                                  n_basis = 1, #Number of bases considered 
                                  covariates  = NULL, #Additional covariates to regress on
                                  L = 100, # How many penalization term in the Lasso to try
                                  K = 5,
                                  thres.ctrl = c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0), # recognition threshold epsilon_n = thres.ctrl * lambda_n,
                                  verbose = FALSE,
                                  tol.abs =1e-4 ,
                                  tol.rel = 1e-4,
                                  eps = 1e-08) {
  
  
  time.start.gX <- proc.time()
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
  if (n < (q+1)*Mp) {
    warning("The sample size is too small! Network estimate may be unreliable!")
  }

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
  d_out <- list()
  norm.adj <- rep(NA,p)
  for(k in 1:p){
    d_out[[k]] <- d.array[k,]
    norm.adj[k] <- norm(d_out[[k]],"2")
  }
  
  
  lambda.max.gX <- lambda.sup.global.two.groups(interM,M, p)
  lambdas.gX <- seq(lambda.max.gX, 0, length.out=L)
  
  P.def <- matrix(0, 2*(p-1)*M, M)
  Q.def <- matrix(0.1, 2*(p-1)*M, M)
  U.def <- matrix(0.01, 2*(p-1)*M, M)
  
  results <- foreach(j = 1:p, .combine="rbind", .packages="fda", .export=c('ADMM.grplasso.two.groups','lambda.sup',
                                                                           'lambda.sup.global.two.groups', 'grp.soft.thres.two.groups',
                                                                           'objective.function.two.groups')) %dopar% {
    message("Processing node ", j)
    jth.range.y <- (j-1)*M+(1:M)
    A.Y <- as.matrix(interM[, jth.range.y])
    jth.range.x <- c(jth.range.y,(j+p-1)*M+(1:M))
    A.X <- as.matrix(interM[, -jth.range.x])
    groups <- temp_groups[-jth.range.x]
    
    P <- P.def; Q <- Q.def; U <- U.def
    
    SCV.mat <- matrix(NA, L, len.t)
    
    N.hat.j.dict <- list()
    
    for(l in 1:L){
      lambda <- lambdas.gX[l]
      # Step 1. Use the full dataset to estimate B.hat(lambda)
      cat(paste("Lambda ",l," of node ",j,"\n",sep=""))
      grp.lasso.result <- tryCatch({
        ADMM.grplasso.two.groups(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
                                 lambda=lambda, rho.init=1,
                                 P.in=P, Q.in=Q, U.in=U,
                                 tol.abs = tol.abs, tol.rel = tol.rel,
                                 maxiter = 400)
      }, error = function(e) {
        message("Error in ADMM.grplasso.two.groups: ", e$message)
        return(NULL)
      })
      if (is.null(grp.lasso.result)) next
      P <- grp.lasso.result$P
      Q <- grp.lasso.result$Q
      U <- grp.lasso.result$U
      

      n_blocks <- nrow(P)/M
      # Process the estimated P into a neighborhood selection vector
      P.frob <- list()
      for(k in 1:n_blocks){
        key <- unique(groups[(k-1)*M + (1:M)])
        if(length(key) > 1){
          message("Error in group definition for block ", k)
          next
        } else {
          P.frob[[key]] <- norm(P[(k-1)*M + (1:M), ], "F")
        }
      }
      N.hat.jl <- list()
      for(ind.t in 1:len.t){
        threshold <- lambda * thres.ctrl[ind.t]
        # Step 2. finding N.hat.j according to each threhold defined by a combination of lambda  and thres.ctrl
        N.hat.jlt <- rep(FALSE, length(P.frob))
        for(n_block in 1:length(P.frob)){
            if(!(is.null(P.frob[[n_block]]))){
              if(P.frob[[n_block]] > threshold){
              N.hat.jlt[n_block] <- TRUE} 
            }
        }
        if(j == p){N.hat.jlt <- c(N.hat.jlt,FALSE)}
        N.hat.jl[[ind.t]]<- N.hat.jlt
        card.N.hat.jlt <- sum(N.hat.jlt)
        slice.pM <- rep(N.hat.jlt, each=M)
        # Check that takes the right columns
        A.X.eff <- as.matrix(interM[, slice.pM]) # select only the columns of A.X with selected features
 
        SCV.single <- rep(NA, K)
        for(k in 1:K){
          if(card.N.hat.jlt > 0){  # Circumstance A. if |N.hat| >=1
            
            # Step A3: slicing out kth fold as CV set, the rest as training set
            fold.k.ind <- seq(floor((k-1)*n/K) + 1, floor(k*n/K)) # index set for the kth fold of SCV
            fold.k.size <- length(fold.k.ind) # cardinality of Ik under the paper's notation
            A.X.cv <- A.X.eff[fold.k.ind, ]
            A.X.train <- A.X.eff[-fold.k.ind, ]
            A.Y.cv <- A.Y[fold.k.ind, ]
            A.Y.train <- A.Y[-fold.k.ind, ]
            
            # Step A4: calculate B.tilde from a ridge regression
            B.tilde <- solve(t(A.X.train) %*% A.X.train + 0.1*diag(card.N.hat.jlt * M)) %*% t(A.X.train) %*% A.Y.train
            
            # Step A5: evaluate SCV for this single k
            residual <- A.Y.cv - A.X.cv %*% B.tilde # residual is n*M matrix
            emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
            SCV.single[k] <- norm(residual,"F")^2 + log(fold.k.size) * card.N.hat.jlt
            
          }else{ # Circumstance B. if cardinality of N.hat is 0
            # Step B3: slicing out kth fold as CV set, the rest as training set
            fold.k.ind <- seq(floor((k-1)*n/K) + 1, floor(k*n/K)) # index set for the kth fold of SCV
            fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
            A.Y.cv <- A.Y[fold.k.ind, ]
            A.Y.train <- A.Y[-fold.k.ind, ]
            
            # Step B5: evaluate SCV for this single k
            residual <- A.Y.cv # residual is n*M matrix
            emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
            SCV.single[k] <- norm(residual,"F")^2
          }
      }# end of the k'th fold
      # Step 6: averaging SCV over K folds for this specific lambda and epsilon
      SCV.mat[l, ind.t] <- mean(SCV.single)
      } # end of all thresholds (epsilon) 1:len.t
      N.hat.j.dict[[l]] <- N.hat.jl
    } # end of lambda 1:L
    scv.min <- which(SCV.mat == min(SCV.mat), arr.ind=T)
    index.optimal <- scv.min[dim(scv.min)[1], ] # there could be multiple. Take the last one
    l.optimal <- index.optimal[1]
    ind.t.optimal <- index.optimal[2]
    lambda.optimal <- lambdas.gX[l.optimal]
    t.optimal <- thres.ctrl[ind.t.optimal]

    N.hat.optimal <- N.hat.j.dict[[l.optimal]][[ind.t.optimal]]
    list(N.hat.optimal = N.hat.optimal, N.hat.all = N.hat.j.dict)
  }
  

  for (j in 1:p) {
    G.mat[j,] <- results[j,1][[1]]
  }
  G.all = list()
    for(l in 1:L){
      G.all[[l]] <- list()
      for(ind.t in 1:len.t){
        G.temp <- matrix(NA, p, (q+1)*p)
        for (j in 1:p) {
          G.temp[j,] <- results[j,2][[1]][[l]][[ind.t]]
        }
        G.all[[l]][[ind.t]] <-  G.temp
      }
      }


  time.end.gX <- proc.time()
  runtime.gX <- (time.end.gX - time.start.gX)[3]
  cat("Computational time of: ")
  cat(runtime.gX)
  return(list(
      G.mat = G.mat,
      G.all = G.all
    ))
}


compute_Nj <- function(A.X, A.Y, M, d, P, Q, U, lambda, threshold, add=FALSE){
  grp.lasso.result <- ADMM.grplasso.two.groups(A.X = A.X, A.Y = A.Y, d = d,
                             lambda=lambda, rho.init=1,
                             P.in=P, Q.in=Q, U.in=U,
                             tol.abs = 1e-4, tol.rel = 1e-4,
                             maxiter = 400)
  
  P <- grp.lasso.result$P
  Q <- grp.lasso.result$Q
  U <- grp.lasso.result$U

  P.out <<- grp.lasso.result$P
  Q.out <<- grp.lasso.result$Q
  U.out <<- grp.lasso.result$U
  
  n_blocks <- nrow(P)/M
  # Process the estimated P into a neighborhood selection vector
  P.frob <- list()
  for(k in 1:n_blocks){
    key <- unique(groups[(k-1)*M + (1:M)])
    if(length(key) > 1){
      cat(message("Error in group definition for block ", k))
      next
    } else {
      P.frob[[key]] <- norm(P[(k-1)*M + (1:M), ], "F")
    }
  }
  
  sel_columns = rep(FALSE, n_blocks)
  iter = 0
  N.hat.jlt <- rep(FALSE, length(P.frob))
  for(n_block in 1:length(P.frob)){
    if(!(is.null(P.frob[[n_block]]))){
      iter = iter + 1
      if(P.frob[[n_block]] > threshold){
        N.hat.jlt[n_block] <- TRUE
        sel_columns[iter] = TRUE } 
    }
  }
  if(add){N.hat.jlt <- c(N.hat.jlt,FALSE)}
  
  return(list(N.hat=N.hat.jlt,sel_columns=sel_columns))
}

compute_SCV <- function(N.hat.tot,K, A.X, A.Y){
  N.hat = N.hat.tot[[1]]
  sel_columns = N.hat.tot[[2]]
  card.N.hat <- sum(N.hat)
  slice.pM <- rep(sel_columns, each=M)
  # Check that takes the right columns
  A.X.eff <- as.matrix(A.X[, slice.pM]) # select only the columns of A.X with selected features
  
  SCV.single <- rep(NA, K)
  for(k in 1:K){
    if(card.N.hat > 0){  # Circumstance A. if |N.hat| >=1
      
      # Step A3: slicing out kth fold as CV set, the rest as training set
      fold.k.ind <- seq(floor((k-1)*n/K) + 1, floor(k*n/K)) # index set for the kth fold of SCV
      fold.k.size <- length(fold.k.ind) # cardinality of Ik under the paper's notation
      A.X.cv <- A.X.eff[fold.k.ind, ]
      A.X.train <- A.X.eff[-fold.k.ind, ]
      A.Y.cv <- A.Y[fold.k.ind, ]
      A.Y.train <- A.Y[-fold.k.ind, ]
      
      # Step A4: calculate B.tilde from a ridge regression
      B.tilde <- solve(t(A.X.train) %*% A.X.train + 0.1*diag(card.N.hat * M)) %*% t(A.X.train) %*% A.Y.train
      
      # Step A5: evaluate SCV for this single k
      residual <- A.Y.cv - A.X.cv %*% B.tilde # residual is n*M matrix
      emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
      SCV.single[k] <- norm(residual,"F")^2 + log(fold.k.size) * card.N.hat
      
    }else{ # Circumstance B. if cardinality of N.hat is 0
      # Step B3: slicing out kth fold as CV set, the rest as training set
      fold.k.ind <- seq(floor((k-1)*n/K) + 1, floor(k*n/K)) # index set for the kth fold of SCV
      fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
      A.Y.cv <- A.Y[fold.k.ind, ]
      A.Y.train <- A.Y[-fold.k.ind, ]
      
      # Step B5: evaluate SCV for this single k
      residual <- A.Y.cv # residual is n*M matrix
      emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
      SCV.single[k] <- norm(residual,"F")^2
    }
  }# end of the k'th fold
  # Step 6: averaging SCV over K folds for this specific lambda and epsilon
  return( mean(SCV.single) )
}

optimize_hyperparameters <- function(A.X, A.Y, M, d, P, Q, U, add, K, lambda, threshold) {
  N.hat.tot <- compute_Nj(A.X=A.X, A.Y=A.Y, M=M, d=d, P=P, Q=Q, U=U, lambda=lambda, threshold=threshold, add=add)
  # Calculate SCV (Split Cross-Validation) loss
  SCV.loss <- compute_SCV(N.hat.tot=N.hat.tot,K=K, A.X=A.X, A.Y=A.Y) # Define this function based on SCV computation
  
  return(list(Score = -SCV.loss, Pred = 0)) # Bayesian optimization maximize while we want to minimize, so negate SCV
}

optimize_hyperparameters_wrapper <- function(lambda, threshold) {
  res = optimize_hyperparameters(A.X=A.X.out, A.Y=A.Y.out, M=M.out, d=d.out[[j]], P=P.out, Q=Q.out,
                                   U=U.out, add=add.out, K=K.out, lambda, threshold)
  return(res)
}


Fast_FGGReg_diff_two_groups_SCV <- function(scores, # functional score on a defined basis, nrow: subjects; ncol: functions*n_basis.
                                            n_basis = 8, #Number of bases considered 
                                            covariates  = NULL, #Additional covariates to regress on
                                            L = 100, # How many penalization term in the Lasso to try
                                            K = 5,
                                            thresholds = c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0), 
                                            verbose = FALSE,
                                            tol.abs =1e-4 ,
                                            tol.rel = 1e-4,
                                            eps = 1e-08,
                                            name_log = "") {
  
  
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
  N <- foreach(j = 1:p, .combine="rbind", .packages="fda", 
  .export=c('ADMM.grplasso.two.groups','lambda.sup','grp.soft.thres.two.groups',
  'objective.function.two.groups', 'optimize_hyperparameters_wrapper', 'optimize_hyperparameters',
  'compute_Nj', 'compute_SCV')) %dopar% {

    time.start <- Sys.time()
    sink(paste(name_log,j, "log.txt",sep="_"), append=TRUE)  
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

    opt_results <- BayesianOptimization(
    FUN = optimize_hyperparameters_wrapper,
    bounds = list(lambda = c(0,lambda.max),
                  threshold = c(0,2)),
    init_points = 10,  # Number of random initial points
    n_iter = 30,       # Number of optimization iterations
    acq = "ucb",
    kappa = 3,         # Expected Improvement criterion
    verbose = verbose
    )

    best_lambda <- opt_results$Best_Par["lambda"]
    best_threshold <- opt_results$Best_Par["threshold"]
    
    N.hat.optimal <-  compute_Nj(A.X=A.X.out, A.Y=A.Y.out, M=M.out, d=d.out[[j]],
                           P=P.out, Q=Q.out, U=U.out, lambda =best_lambda, threshold=best_threshold, add=add.out)

    N.hat.optimal$N.hat
    cat("Computational time of: ")
    cat( Sys.time() - time.start )
    sink()
  }
  
  return(N)
}

############################# UPDATED prec.rec including F1 value
prec.rec <- function(G.true, G.mat, type=c("AND","OR")){
  # a function to calculate TP, FP, TN, FN
  # and return precision, recall, TPR and FPR.
  # Input:
  #   G.true, the true p*p adjacency matrix
  #   G.mat, the estimated p*p adjacency matrix
  #   type,
  #     AND: when two nodes both recognize each other as neighbors
  #     OR: when either of them recognize each other as neighbors
  # Output:
  #   prec, precision
  #   rec, recall
  #   TPR and FPR, as their names

  
  p <- nrow(G.true)
  TP <- 0; TN <- 0; FP <- 0; FN <- 0
  
  if(type=="AND"){
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(G.mat[i,j] == 1 & G.mat[j,i] == 1){
          if(G.true[i,j] == 1)
            TP <- TP + 1
          else
            FP <- FP + 1
        }else{
          if(G.true[i,j] == 1)
            FN <- FN + 1
          else
            TN <- TN + 1
        }
      }
    }
  }else if(type=="OR"){
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(G.mat[i,j] == 1 | G.mat[j,i] == 1){
          if(G.true[i,j] == 1)
            TP <- TP + 1
          else
            FP <- FP + 1
        }else{
          if(G.true[i,j] == 1)
            FN <- FN + 1
          else
            TN <- TN + 1
        }
      }
    }
  }
  prec <- TP / (TP + FP)
  if(TP+FP==0) prec <- 1
  
  rec <- TP / (TP + FN)
  if(TP+FN==0) rec <- 0
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  F1 <- 2*TP/(2*TP+FP+FN)
  return(list(prec=prec, rec=rec, TPR=TPR, FPR=FPR, F1=F1))
}


