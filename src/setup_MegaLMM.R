setup_MegaLMM <- function(Y, formula, sample_data, kinship_matrix, 
                          latent_factors = 15, run_ID = "MegaLMM_run",
                          h2_divisions = 20, burn = 0, thin = 2,
                          h2_step_size = NULL,
                          save_current_state = TRUE) 
{
  #' Setup MegaLMM model with optimized parameters
  #' 
  #' @param Y Training phenotype matrix
  #' @param formula Model formula
  #' @param sample_data Data frame with individual information
  #' @param kinship_matrix Kinship/relationship matrix
  #' @param latent_factors Number of latent factors
  #' @param run_ID Identifier for the model run
  #' @param h2_divisions Number of divisions for h2 prior
  #' @param burn Burn-in iterations
  #' @param thin Thinning interval
  #' @param h2_step_size Step size for h2 prior. Defaults to NULL.
  #' @param save_current_state Save the current state of the model. Defaults to TRUE.
  #' 
  #' @return Initialized MegaLMM state object
  #' ___________________________________________________________________________
  
  # Create MegaLMM control parameters first with very explicit naming
  control_params <- list(
    h2_divisions = h2_divisions,
    burn = burn,  
    thin = thin,
    K = latent_factors,  # K in MegaLMM_control refers to number of factors
    h2_step_size = h2_step_size,
    save_current_state = save_current_state
  )
  
  # Use do.call to avoid any evaluation errors
  run_parameters <- do.call(MegaLMM_control, control_params)
  
  # Process kinship matrix
  relmat_var <- NULL
  if (is.list(kinship_matrix)) {
    relmat_var <- names(kinship_matrix)[1]
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
    K_list[[relmat_var]] <- kinship_matrix
    kinship_matrix <- K_list
  }
  
  # Setup the model
  MegaLMM_state <- setup_model_MegaLMM(
    Y = Y,
    formula = formula,
    data = sample_data,
    relmat = kinship_matrix,
    run_parameters = run_parameters,
    run_ID = run_ID
  )
  
  return(MegaLMM_state)
}