run_MegaLMM_sampling <- function(MegaLMM_state, 
                                 sampling_rounds = 4, 
                                 iterations_per_round = 250) 
{
  #' Run MegaLMM sampling phase
  #' 
  #' @param MegaLMM_state MegaLMM state object
  #' @param sampling_rounds Number of sampling rounds
  #' @param iterations_per_round Iterations per sampling round
  #' @return MegaLMM state object after sampling
  #' ___________________________________________________________________________
  
  for (i in 1:sampling_rounds) {
    cat(sprintf('Sampling round %d of %d\n', i, sampling_rounds))
    
    # Run iterations
    MegaLMM_state <- sample_MegaLMM(MegaLMM_state, iterations_per_round)
    
    # Save posterior chunk to disk
    MegaLMM_state <- save_posterior_chunk(MegaLMM_state)
    
    cat(sprintf('Completed %d total samples\n', MegaLMM_state$current_state$nrun))
  }
  
  return(MegaLMM_state)
}