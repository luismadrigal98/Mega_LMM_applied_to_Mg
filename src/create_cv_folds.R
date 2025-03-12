create_cv_folds <- function(Y, sample_data, id_col, group_col = NULL, k = 5, 
                            seed = NULL) 
{
  #' Create k-fold cross-validation partitions
  #' 
  #' @param Y Phenotype matrix
  #' @param sample_data Data frame with individual information
  #' @param id_col Column name for individual IDs
  #' @param group_col Column name for grouping variable (e.g., Population)
  #' @param k Number of folds
  #' @param seed Random seed for reproducibility
  #' @return Matrix of fold assignments
  #' ___________________________________________________________________________
  
  if (!is.null(seed)) set.seed(seed)
  
  fold_ID_matrix <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y), dimnames = dimnames(Y))
  
  for (i in 1:ncol(fold_ID_matrix)) {
    # Get individuals with non-missing data for this environment
    observed_idx <- which(!is.na(Y[, i]))
    observed_lines <- sample_data[observed_idx, ]
    
    if (!is.null(group_col) && group_col %in% colnames(sample_data)) {
      # Use stratified sampling if group column is provided
      groups <- unique(observed_lines[[group_col]])
      for (grp in groups) {
        grp_idx <- which(observed_lines[[group_col]] == grp)
        n_lines <- length(grp_idx)
        if (n_lines > 0) {
          observed_lines$fold[grp_idx] <- sample(rep(1:k, length.out = n_lines))
        }
      }
    } else {
      # Regular random assignment
      n_lines <- nrow(observed_lines)
      observed_lines$fold <- sample(rep(1:k, length.out = n_lines))
    }
    
    # Assign fold IDs to the matrix
    fold_ID_matrix[match(observed_lines[[id_col]], rownames(fold_ID_matrix)), i] <- observed_lines$fold
  }
  
  return(fold_ID_matrix)
}