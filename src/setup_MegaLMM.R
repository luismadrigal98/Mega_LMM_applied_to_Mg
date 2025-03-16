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
  
  # Create MegaLMM control parameters
  control_params <- list(
    h2_divisions = h2_divisions,
    burn = burn,  
    thin = thin,
    K = latent_factors,  # K in MegaLMM_control refers to number of factors
    h2_step_size = h2_step_size,
    save_current_state = save_current_state
  )
  
  run_parameters <- do.call(MegaLMM_control, control_params)
  
  # Process kinship matrix - check if it's already a list
  if (is.list(kinship_matrix) && !is.matrix(kinship_matrix)) {
    # Already a list, use as is
    relmat <- kinship_matrix
  } else {
    # Convert to list - need to extract random effect term
    # For a one-sided formula ~A+B+(1|G), formula[[1]] is ~, formula[[2]] is the right side
    formula_terms <- terms(formula)
    
    # Check if there are random effects terms
    random_terms <- attr(formula_terms, "term.labels")
    random_terms <- random_terms[grepl("\\|", random_terms)]
    
    if (length(random_terms) == 0) {
      stop("No random effect term (1|Group) found in formula")
    }
    
    # Extract the group variable (e.g., 'G' from '(1|G)')
    relmat_var <- sub(".*\\|", "", random_terms[1])
    relmat_var <- trimws(relmat_var)
    
    # Create the named list for relmat
    relmat <- list()
    relmat[[relmat_var]] <- kinship_matrix
  }
  
  # Setup the model
  MegaLMM_state <- setup_model_MegaLMM(
    Y = Y,
    formula = formula,
    data = sample_data,
    relmat = relmat,
    run_parameters = run_parameters,
    run_ID = run_ID
  )
  
  return(MegaLMM_state)
}