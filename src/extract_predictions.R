extract_predictions <- function(MegaLMM_state, prediction_type = "both") 
{
  #' Extract predictions from MegaLMM model
  #' 
  #' @param MegaLMM_state MegaLMM state object
  #' @param prediction_type Type of prediction to extract ('U', 'Eta_mean', or 'both')
  #' @return List containing predictions
  #' ___________________________________________________________________________
  
  result <- list()
  
  if (prediction_type %in% c("U", "both")) {
    U_samples <- load_posterior_param(MegaLMM_state, 'U')
    result$U_hat <- get_posterior_mean(U_samples)
  }
  
  if (prediction_type %in% c("Eta_mean", "both")) {
    result$Eta_mean <- load_posterior_param(MegaLMM_state, 'Eta_mean')
  }
  
  return(result)
}