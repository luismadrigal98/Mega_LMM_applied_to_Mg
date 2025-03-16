optimize_missing_data <- function(MegaLMM_state, max_groups = NULL) 
{
  #' Optimize missing data handling for MegaLMM
  #' 
  #' @param MegaLMM_state MegaLMM state object
  #' @param max_groups Maximum number of missing data groups to consider
  #' @return MegaLMM state object with optimized missing data map
  #' ___________________________________________________________________________
  
  if (is.null(max_groups)) {
    max_groups <- ncol(MegaLMM_state$Y) + 1
  }
  
  # Create missing data maps
  maps <- make_Missing_data_map(MegaLMM_state, 
                                max_NA_groups = max_groups, 
                                verbose = FALSE)
  
  # Find optimal map based on number of groups and NAs handled
  map_results <- maps$map_results
  
  # Select map based on diminishing returns principle
  na_diffs <- diff(map_results$n_NAs_dropped)
  cutoff_idx <- which(na_diffs/max(na_diffs) < 0.05)[1]
  if (is.na(cutoff_idx)) cutoff_idx <- nrow(map_results)
  optimal_idx <- min(cutoff_idx, nrow(map_results))
  
  # Apply selected map
  MegaLMM_state <- set_Missing_data_map(MegaLMM_state, 
                                        maps$Missing_data_map_list[[optimal_idx]])
  
  return(MegaLMM_state)
}