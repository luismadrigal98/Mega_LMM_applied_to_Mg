set_megalMM_priors <- function(MegaLMM_state, lambda_prior_type = "horseshoe") 
{
  #' Set priors for MegaLMM model
  #' 
  #' @param MegaLMM_state MegaLMM state object
  #' @param lambda_prior_type Type of prior for Lambda ('horseshoe', 'ARD', or 
  #' 'BayesC')
  #' @return MegaLMM state object with priors set
  #' ___________________________________________________________________________
  
  # Set Lambda prior based on type
  if (lambda_prior_type == "horseshoe") {
    Lambda_prior <- list(
      sampler = sample_Lambda_prec_horseshoe,
      Lambda_df = 3,
      delta_1 = list(shape = 2, rate = 1),
      delta_2 = list(shape = 3, rate = 1),
      delta_iterations_factor = 100
    )
  } else if (lambda_prior_type == "ARD") {
    Lambda_prior <- list(
      sampler = sample_Lambda_prec_ARD,
      prop_0 = 0.1,
      delta = list(shape = 3, scale = 1),
      delta_iterations_factor = 100
    )
  } else if (lambda_prior_type == "BayesC") {
    Lambda_prior <- list(
      sampler = sample_Lambda_prec_BayesC,
      prop_0 = 0.1,
      delta = list(shape = 3, scale = 1),
      delta_iterations_factor = 100
    )
  } else {
    stop("Invalid lambda_prior_type. Choose 'horseshoe', 'ARD', or 'BayesC'")
  }
  
  # Set remaining priors
  priors <- MegaLMM_priors(
    tot_Y_var = list(V = 0.5, nu = 5),
    tot_F_var = list(V = 18/20, nu = 20),
    h2_priors_resids_fun = function(h2s, n) 1,
    h2_priors_factors_fun = function(h2s, n) 1,
    Lambda_prior = Lambda_prior
  )
  
  # Apply priors to model
  MegaLMM_state <- set_priors_MegaLMM(MegaLMM_state, priors)
  
  return(MegaLMM_state)
}