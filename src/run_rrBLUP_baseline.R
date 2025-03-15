run_rrBLUP_baseline <- function(Y, sample_data, K, formula, 
                                return_all_BLUP_values = FALSE,
                                parallel = TRUE) 
{
  #' Run rrBLUP as a baseline comparison model
  #' 
  #' @param Y Training phenotype matrix
  #' @param sample_data Data frame with individual information
  #' @param K Kinship/relationship matrix
  #' @param formula Model formula for fixed effects
  #' @param return_all_BLUP_values Logical indicating whether to return all BLUP 
  #' values produced by mixed.solve.
  #' @param parallel Logical indicating whether to parallelize trait processing.
  #' Should be FALSE when called from parallel cross-validation functions.
  #' @return Matrix of rrBLUP predictions or list containing predictions and BLUP values
  #' ___________________________________________________________________________
  
  # Extract fixed effect formula parts
  fixed_formula <- as.formula(paste("~", as.character(formula)[2]))
  
  # Create design matrix
  X <- model.matrix(fixed_formula, sample_data)
  
  # Initialize predictions matrix
  rrBLUP_predictions <- matrix(NA, nrow(Y), ncol(Y), 
                               dimnames = dimnames(Y))
  
  # Initialize container for BLUP results
  BLUP_per_trait <- list()
  
  if (parallel) {
    # Use foreach with parallel backend
    library(foreach)
    library(doParallel)
    cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(cores)
    registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    
    results <- foreach(i = 1:ncol(Y), .packages = "rrBLUP") %dopar% {
      # Process trait code...
      X_i <- X
      if (ncol(X) > 1) {
        var_check <- apply(X[!is.na(Y[, i]), , drop = FALSE], 2, var)
        if (any(var_check == 0)) {
          X_i <- X[, var_check > 0, drop = FALSE]
        }
      }
      
      res <- mixed.solve(y = Y[, i], X = X_i, K = K)
      predictions <- c(X_i %*% res$beta) + res$u
      
      list(predictions = predictions, model_fit = res)
    }
  } else {
    # Use standard loop when called from other parallel functions
    results <- vector("list", ncol(Y))
    
    for (i in 1:ncol(Y)) {
      X_i <- X
      if (ncol(X) > 1) {
        var_check <- apply(X[!is.na(Y[, i]), , drop = FALSE], 2, var)
        if (any(var_check == 0)) {
          X_i <- X[, var_check > 0, drop = FALSE]
        }
      }
      
      res <- mixed.solve(y = Y[, i], X = X_i, K = K)
      predictions <- c(X_i %*% res$beta) + res$u
      
      results[[i]] <- list(predictions = predictions, model_fit = res)
    }
  }
  
  # Extract results from foreach/loop
  for (i in 1:ncol(Y)) {
    rrBLUP_predictions[, i] <- results[[i]]$predictions
    BLUP_per_trait[[colnames(Y)[i]]] <- results[[i]]$model_fit
  }
  
  if (return_all_BLUP_values) {
    return(list(predictions = rrBLUP_predictions, 
                BLUP_values = BLUP_per_trait))
  }
  else return(rrBLUP_predictions)
}