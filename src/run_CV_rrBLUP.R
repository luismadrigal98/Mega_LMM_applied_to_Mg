run_CV_rrBLUP <- function(Y, fold_matrix, sample_data, K, formula)
{
  #' This function is designed to execute run_rrBLUP_baseline in parallel and using
  #' cross validation. This will produce a set of estimates per fold and provide
  #' a measure of accuracy that could be used to assess its performance.
  #' 
  #' @param Y Trait matrix, where rownames contain the Individual_ID
  #' @param fold_matrix A fold assignment matrix with same dimensions as Y
  #' @param sample_data Data frame with individual information needed for fixed effects
  #' @param K Kinship/relationship matrix
  #' @param formula Model formula for fixed effects
  #' 
  #' @return A list with the prediction per fold, the estimated accuracy per trait,
  #' and an overall accuracy across all folds and traits.
  #' ___________________________________________________________________________
  
  # Ensure fold_matrix has same dimensions as Y
  if(!identical(dim(Y), dim(fold_matrix)) || !identical(rownames(Y), rownames(fold_matrix))) {
    stop("fold_matrix must have same dimensions and rownames as Y")
  }
  
  # Retrieving the number of folds in the fold_matrix
  folds <- sort(unique(as.vector(fold_matrix)))
  n_folds <- length(folds)
  n_traits <- ncol(Y)
  
  # Prepare containers for results
  predictions <- Y * NA
  predictions_train <- Y * NA
  accuracies <- matrix(NA, n_folds, n_traits, 
                       dimnames = list(paste0("Fold_", folds), colnames(Y)))
  
  # Outer loop over folds
  results <- foreach(fold = folds, .combine = 'list', .multicombine = TRUE,
                     .packages = c("rrBLUP", "stats", "utils"),
                     .export = c("run_rrBLUP_baseline")) %dopar% {
                       
                       fold_predictions <- Y * NA
                       fold_predictions_train <- Y * NA
                       fold_accuracies <- numeric(n_traits)
                       names(fold_accuracies) <- colnames(Y)
                       
                       # Inner loop over traits
                       for(trait_idx in 1:n_traits) {
                         # Get trait name
                         trait <- colnames(Y)[trait_idx]
                         
                         # Create train/test split for this trait and fold
                         is_test <- fold_matrix[, trait_idx] == fold
                         is_train <- !is_test & !is.na(Y[, trait_idx])
                         
                         # Skip if no test samples for this fold
                         if(sum(is_test, na.rm = TRUE) == 0) next
                         
                         # Prepare training data
                         Y_train <- Y[is_train, , drop = FALSE]
                         sample_data_train <- sample_data[rownames(Y_train), , drop = FALSE]
                         
                         # Run rrBLUP on training data
                         predictions_all <- run_rrBLUP_baseline(Y_train, sample_data_train, 
                                                                K, formula,
                                                                parallel = FALSE)
                         
                         # FIXED: Store training predictions by row names to avoid dimension issues
                         if(is.list(predictions_all) && "predictions" %in% names(predictions_all)) {
                           fold_predictions_train[rownames(predictions_all$predictions), ] <- predictions_all$predictions
                         } else {
                           fold_predictions_train[rownames(predictions_all), ] <- predictions_all
                         }
                         
                         # Prepare test data
                         sample_data_test <- sample_data[is_test, , drop = FALSE]
                         
                         # Extract fixed effect formula parts
                         fixed_formula <- as.formula(paste("~", as.character(formula)[2]))
                         
                         # Create design matrices
                         X_train <- model.matrix(fixed_formula, sample_data_train)
                         X_test <- model.matrix(fixed_formula, sample_data_test)
                         
                         # FIXED: Ensure K has proper row/column names and dimensions
                         K_train <- K[rownames(Y_train), rownames(Y_train), drop = FALSE]
                         if(any(is.na(K_train))) warning("NA values in subsetted kinship matrix")
                         
                         # Get predictions for test samples
                         res <- mixed.solve(y = Y_train[, trait_idx], X = X_train, K = K_train)
                         
                         # FIXED: Handle potential mismatch in test predictions
                         matched_idx <- match(rownames(sample_data_test), rownames(Y_train))
                         test_u <- rep(0, nrow(sample_data_test))
                         valid_idx <- !is.na(matched_idx)
                         if(any(valid_idx)) {
                           test_u[valid_idx] <- res$u[matched_idx[valid_idx]]
                         }
                         test_preds <- c(X_test %*% res$beta) + test_u
                         
                         # Store test predictions
                         fold_predictions[is_test, trait_idx] <- test_preds
                         
                         # Calculate accuracy
                         observed <- Y[is_test, trait_idx]
                         predicted <- fold_predictions[is_test, trait_idx]
                         fold_accuracies[trait_idx] <- cor(observed, predicted, use = "pairwise.complete.obs")
                       }
                       
                       list(
                         predictions = fold_predictions,
                         predictions_train = fold_predictions_train,
                         accuracies = fold_accuracies
                       )
                     }
  
  # Combine results from all folds
  for(i in 1:n_folds) {
    fold_result <- results[[i]]
    
    # Update predictions - only fill in the test predictions from each fold
    is_pred <- !is.na(fold_result$predictions)
    if(any(is_pred)) predictions[is_pred] <- fold_result$predictions[is_pred]
    
    # Update training predictions
    is_train_pred <- !is.na(fold_result$predictions_train)
    if(any(is_train_pred)) predictions_train[is_train_pred] <- fold_result$predictions_train[is_train_pred]
    
    # Store accuracies
    accuracies[i,] <- fold_result$accuracies
  }
  
  # Calculate overall accuracies per trait
  overall_accuracies <- apply(accuracies, 2, mean, na.rm = TRUE)
  
  # Return results
  return(list(
    predictions = predictions,
    predictions_train = predictions_train,
    fold_accuracies = accuracies,
    trait_accuracies = overall_accuracies,
    overall_accuracy = mean(overall_accuracies, na.rm = TRUE)
  ))
}