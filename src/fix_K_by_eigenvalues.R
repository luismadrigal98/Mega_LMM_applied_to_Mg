fix_K_by_eigenvalues <- function(K, max_iter = 10, min_eigenvalue = 1e-6, 
                                 verbose = FALSE, error_on_fail = TRUE) 
{
  #' Fix a Matrix to be Positive Semi-Definite by Eigenvalue Correction
  #'
  #' @description
  #' Ensures a matrix is positive semi-definite (PSD) by correcting eigenvalues below a
  #' specified threshold. This function is particularly useful for genetic relationship 
  #' matrices, variance-covariance matrices, or correlation matrices that must be PSD 
  #' for statistical analyses like GBLUP, REML, or mixed models.
  #'
  #' @param K A square, symmetric matrix to be made positive semi-definite. Typically
  #'   a relationship, kinship, or correlation matrix.
  #' @param max_iter An integer specifying the maximum number of iterations to attempt
  #'   while correcting the matrix. Default is 10.
  #' @param min_eigenvalue A numeric value specifying the minimum allowed eigenvalue.
  #'   Default is 1e-6, which balances numerical stability with minimal distortion
  #'   of the original matrix structure.
  #' @param verbose Logical; if TRUE, prints progress information during the correction
  #'   process. Default is FALSE.
  #' @param error_on_fail Logical; if TRUE, the function stops with an error if the matrix
  #'   couldn't be made PSD after max_iter iterations. If FALSE, it issues a warning
  #'   instead. Default is TRUE.
  #'
  #' @details
  #' The function performs eigendecomposition and replaces eigenvalues below the minimum
  #' threshold with that threshold value. The matrix is then reconstructed and symmetrized.
  #' If unsuccessful, the minimum threshold is automatically increased with each iteration.
  #' 
  #' The correction preserves the overall structure of the matrix while ensuring numerical
  #' properties required for statistical analyses. A common use case is ensuring genetic
  #' relationship matrices are PSD before variance component estimation.
  #'
  #' @return A list containing:
  #' \itemize{
  #'   \item \code{matrix}: The corrected matrix (or original if already PSD)
  #'   \item \code{success}: Logical indicating whether the correction was successful
  #'   \item \code{iterations}: Number of iterations performed
  #' }
  #'
  #' @examples
  #' # Create an example relationship matrix that's not PSD
  #' K <- matrix(c(1, 0.5, 0.7, 0.5, 1, 0.9, 0.7, 0.9, 1), 3, 3)
  #' K[1,3] <- K[3,1] <- 0.95  # Make it non-PSD
  #' 
  #' # Check if it's PSD
  #' # matrixcalc::is.positive.semi.definite(K)  # Would return FALSE
  #' 
  #' # Fix the matrix with default settings
  #' K_result <- fix_K_by_eigenvalues(K)
  #' K_fixed <- K_result$matrix
  #' 
  #' # Fix with more output
  #' K_result <- fix_K_by_eigenvalues(K, verbose = TRUE)
  #' 
  #' # Verify it's now PSD
  #' # matrixcalc::is.positive.semi.definite(K_result$matrix)  # Should return TRUE
  #'
  #' @export
  #' ___________________________________________________________________________
  
  # Save the dimnames
  dimnames_K <- dimnames(K)
  
  iter_count <- 0
  initial_min_eig <- min(eigen(K, symmetric = TRUE, only.values = TRUE)$values)
  
  while(iter_count < max_iter && !matrixcalc::is.positive.semi.definite(x = K))
  {
    iter_count <- iter_count + 1
    
    # Increase min_eigenvalue slightly each iteration if needed
    current_min <- min_eigenvalue * (1 + (iter_count - 1) * 0.5)
    
    if(verbose) {
      cat(sprintf("Iteration %d: Using min eigenvalue = %.8g\n", 
                  iter_count, current_min))
    }
    
    # Spectral decomposition
    eigen_K <- eigen(K, symmetric = TRUE)
    values <- eigen_K$values
    vectors <- eigen_K$vectors
    
    # Fix eigenvalues below threshold
    neg_indices <- values < current_min
    if(any(neg_indices)) {
      values[neg_indices] <- current_min
      K <- vectors %*% diag(values) %*% t(vectors)
      K <- (K + t(K))/2  # Ensure symmetry
    }
  }
  
  is_psd <- matrixcalc::is.positive.semi.definite(x = K)
  
  if(!is_psd) {
    msg <- sprintf("Failed to correct matrix after %d iterations", max_iter)
    if(error_on_fail) stop(msg) else warning(msg)
  }
  
  if(verbose) {
    final_min_eig <- min(eigen(K, symmetric = TRUE, only.values = TRUE)$values)
    cat(sprintf("Min eigenvalue changed from %.8g to %.8g\n", 
                initial_min_eig, final_min_eig))
  }
  
  ## Return the dimnames to K
  dimnames(K) <- dimnames_K
  
  return(list(
    matrix = K,
    success = is_psd,
    iterations = iter_count
  ))
}