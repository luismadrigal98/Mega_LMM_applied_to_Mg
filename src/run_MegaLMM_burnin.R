run_MegaLMM_burnin <- function(MegaLMM_state, burnin_rounds = 5, 
                               iterations_per_round = 100, 
                               drop_cor_threshold = 0.6,
                               make_plots = FALSE) 
{
  #' Run MegaLMM burnin phase with automatic factor reordering
  #' 
  #' @param MegaLMM_state MegaLMM state object
  #' @param burnin_rounds Number of burnin rounds
  #' @param iterations_per_round Iterations per burnin round
  #' @param drop_cor_threshold Threshold for dropping correlated factors
  #' @param make_plots Whether to create diagnostic plots
  #' @return MegaLMM state object after burnin
  #' ___________________________________________________________________________
  
  for (i in 1:burnin_rounds) {
    cat(sprintf('Burnin round %d of %d\n', i, burnin_rounds))
    
    # Reorder factors to help with MCMC mixing
    MegaLMM_state <- reorder_factors(MegaLMM_state, drop_cor_threshold = drop_cor_threshold)
    
    # Clear previous samples
    MegaLMM_state <- clear_Posterior(MegaLMM_state)
    
    # Run iterations
    MegaLMM_state <- sample_MegaLMM(MegaLMM_state, iterations_per_round)
    
    # Create diagnostic plots if requested
    if (make_plots) {
      traceplot_array(MegaLMM_state$Posterior$Lambda, 
                      name = file.path(MegaLMM_state$run_ID, sprintf('Lambda_burnin_%d.pdf', i)))
    }
    
    cat(sprintf('Completed %d burnin samples\n', MegaLMM_state$current_state$nrun))
  }
  
  # Final clear before sampling phase
  MegaLMM_state <- clear_Posterior(MegaLMM_state)
  
  return(MegaLMM_state)
}