run_rrBLUP_baseline <- function(Y_train, Y_test, sample_data, K, formula) 
{
  #' Run rrBLUP as a baseline comparison model
  #' 
  #' @param Y_train Training phenotype matrix
  #' @param Y_test Test phenotype matrix 
  #' @param sample_data Data frame with individual information
  #' @param K Kinship/relationship matrix
  #' @param formula Model formula for fixed effects
  #' @return Matrix of rrBLUP predictions
  #' ___________________________________________________________________________
  
  # Extract fixed effect formula parts
  fixed_formula <- as.formula(paste("~", as.character(formula)[2]))
  
  # Create design matrix
  X <- model.matrix(fixed_formula, sample_data)
  
  # Initialize predictions matrix
  rrBLUP_predictions <- matrix(NA, nrow(Y_train), ncol(Y_train), 
                               dimnames = dimnames(Y_train))
  
  # Run rrBLUP separately for each predictor
  for (i in 1:ncol(Y_train)) {
    # Check if we need to drop some fixed effects due to zero variance
    X_i <- X
    if (ncol(X) > 1) {
      var_check <- apply(X[!is.na(Y_train[, i]), , drop = FALSE], 2, var)
      if (any(var_check == 0)) {
        X_i <- X[, var_check > 0, drop = FALSE]
      }
    }
    
    # Fit model
    res <- mixed.solve(y = Y_train[, i], X = X_i, K = K)
    
    # Make predictions
    rrBLUP_predictions[, i] <- c(X_i %*% res$beta) + res$u
  }
  
  return(rrBLUP_predictions)
}