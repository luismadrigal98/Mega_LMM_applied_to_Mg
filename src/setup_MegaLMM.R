setup_MegaLMM <- function(Y, formula, sample_data, K, 
                          factors = 15, run_ID = "MegaLMM_run",
                          h2_divisions = 20, burn = 0, thin = 2,
                          scale_Y = T, h2_step_size = NULL,
                          save_current_state = T) 
{
  #' Setup MegaLMM model with optimized parameters
  #' 
  #' @param Y Training phenotype matrix
  #' @param formula Model formula
  #' @param sample_data Data frame with individual information
  #' @param K Kinship/relationship matrix
  #' @param factors Number of latent factors
  #' @param run_ID Identifier for the model run
  #' @param h2_divisions Number of divisions for h2 prior
  #' @param burn Burn-in iterations
  #' @param thin Thinning interval
  #' @param scale_Y Scale the phenotype matrix. Defaults to TRUE.
  #' @param h2_step_size Step size for h2 prior. Defaults to NULL.
  #' @param save_current_state Save the current state of the model. Defaults to TRUE.
  #' 
  #' @return Initialized MegaLMM state object
  #' ___________________________________________________________________________
  
  # Set control parameters
  run_parameters <- MegaLMM_control(
    h2_divisions = h2_divisions,
    burn = burn,  
    thin = thin,
    K = factors,
    scale_Y = scale_Y,
    h2_step_size = h2_step_size,
    save_current_state = save_current_state
  )
  
  # Extract first column name from K if it's a list
  relmat_var <- NULL
  if (is.list(K)) {
    relmat_var <- names(K)[1]
  } else {
    # Try to infer the variable from the formula
    relmat_var <- all.vars(formula[[3]])
    relmat_var <- relmat_var[grepl("\\|", formula[[3]])]
    if (length(relmat_var) > 0) {
      relmat_var <- gsub(".*\\|", "", relmat_var)
      relmat_var <- trimws(relmat_var)
    } else {
      stop("Could not determine random effect variable from formula")
    }
    K_list <- list()
    K_list[[relmat_var]] <- K
    K <- K_list
  }
  
  # Setup the model
  MegaLMM_state <- setup_model_MegaLMM(
    Y = Y,
    formula = formula,
    data = sample_data,
    relmat = K,
    run_parameters = run_parameters,
    run_ID = run_ID
  )
  
  return(MegaLMM_state)
}