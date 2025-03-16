initialize_MegaLMM <- function(MegaLMM_state, 
                               retain_par = c('Lambda', 'F_h2', 'resid_h2', 'tot_Eta_prec', 'B1'),
                               retain_mean_par = c('Eta_mean'),
                               retain_function = list(
                                 U = 'U_F %*% Lambda + U_R + X1 %*% B1',
                                 G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
                                 R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])',
                                 h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'
                               ),
                               verbose = FALSE) 
{
  #' Initialize MegaLMM model and check memory requirements
  #' 
  #' @param MegaLMM_state MegaLMM state object
  #' @param retain_par Parameters to retain in the posterior
  #' @param retain_mean_par Mean parameters to retain in the posterior
  #' @param retain_function Functions to retain in the posterior
  #' @param verbose Print verbose information
  #' @return Initialized MegaLMM state object
  #' ___________________________________________________________________________
  
  # Initialize variables
  MegaLMM_state <- initialize_variables_MegaLMM(MegaLMM_state)
  
  # Estimate memory requirements
  mem_est <- estimate_memory_initialization_MegaLMM(MegaLMM_state)
  if (verbose) {
    cat("Estimated memory for initialization:", mem_est, "MB\n")
  }
  
  # Initialize model
  MegaLMM_state <- initialize_MegaLMM(MegaLMM_state, verbose = verbose)
  
  # Set up posterior storage for useful parameters and functions
  MegaLMM_state$Posterior$posteriorSample_params <- retain_par
  MegaLMM_state$Posterior$posteriorMean_params <- retain_mean_par
  MegaLMM_state$Posterior$posteriorFunctions <- retain_functions
  
  # Clear and initialize posterior database
  MegaLMM_state <- clear_Posterior(MegaLMM_state)
  
  return(MegaLMM_state)
}
