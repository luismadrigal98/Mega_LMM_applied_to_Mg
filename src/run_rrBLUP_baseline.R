run_rrBLUP_baseline <- function(Y, sample_data, K, formula, 
                                return_all_BLUP_values = FALSE) 
{
  #' Run rrBLUP as a baseline comparison model
  #' 
  #' @param Y Training phenotype matrix
  #' @param sample_data Data frame with individual information
  #' @param K Kinship/relationship matrix
  #' @param formula Model formula for fixed effects
  #' @param return_all_BLUP_values Logical indicating whether to return all BLUP 
  #' values produced by mixed.solve.
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
  
  # Run rrBLUP separately for each predictor
  results <- foreach(i = 1:ncol(Y), .packages = "rrBLUP") %dopar% {
    # Check if we need to drop some fixed effects due to zero variance
    X_i <- X
    if (ncol(X) > 1) {
      var_check <- apply(X[!is.na(Y[, i]), , drop = FALSE], 2, var)
      if (any(var_check == 0)) {
        X_i <- X[, var_check > 0, drop = FALSE]
      }
    }
    
    # Fit model
    res <- mixed.solve(y = Y[, i], X = X_i, K = K)
    
    # Make predictions
    predictions <- c(X_i %*% res$beta) + res$u
    
    # Return both predictions and model results
    list(predictions = predictions, model_fit = res)
  }
  
  # Extract results from foreach
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