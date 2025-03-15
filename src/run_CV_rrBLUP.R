run_CV_rrBLUP <- function(Y, fold_matrix, sample_data, K, formula) {
  #' Cross-validation for rrBLUP using run_rrBLUP_baseline
  #' 
  #' @param Y Trait matrix, where rownames contain the Individual_ID
  #' @param fold_matrix A fold assignment matrix with same dimensions as Y
  #' @param sample_data Data frame with individual information for fixed effects
  #' @param K Kinship/relationship matrix
  #' @param formula Model formula for fixed effects
  #' 
  #' @return A list with predictions, accuracies per trait, and overall accuracy
  #' ___________________________________________________________________________
  
  # Basic validation
  if(!identical(dim(Y), dim(fold_matrix)) || !identical(rownames(Y), rownames(fold_matrix))) {
    stop("fold_matrix must have same dimensions and rownames as Y")
  }
  
  # Setup
  folds <- sort(unique(as.vector(fold_matrix)))
  n_folds <- length(folds)
  n_traits <- ncol(Y)
  trait_names <- colnames(Y)
  all_IDs <- rownames(Y)
  
  # Results containers
  predictions <- Y * NA
  predictions_train <- Y * NA
  accuracies <- matrix(NA, n_folds, n_traits, 
                       dimnames = list(paste0("Fold_", folds), trait_names))
  
  # Process each fold sequentially
  for(f in 1:n_folds) {
    fold <- folds[f]
    cat(sprintf("\nProcessing fold %d/%d\n", f, n_folds))
    
    # For each trait, set up train/test datasets
    Y_train <- Y * NA    # Make a copy with NAs in test positions
    test_samples <- list()  # Track test samples by trait
    
    for(t in 1:n_traits) {
      # Get test samples for this trait/fold
      is_test <- fold_matrix[, t] == fold
      test_samples[[t]] <- rownames(Y)[is_test]
      
      # For this trait, set training values (non-test samples)
      Y_train[!is_test, t] <- Y[!is_test, t]
      
      cat(sprintf("  Trait %d (%s): %d test samples\n", 
                  t, trait_names[t], sum(is_test)))
    }
    
    # Important: No need to subset sample_data!
    # run_rrBLUP_baseline will automatically use only the rows needed
    # based on non-NA values in Y_train
    cat("  Running rrBLUP on training data...\n")
    
    tryCatch({
      # Run rrBLUP on the training data 
      # Note: We pass the COMPLETE sample_data and K matrix
      blup_results <- run_rrBLUP_baseline(
        Y = Y_train,           # Training data (NAs for test samples)
        sample_data = sample_data,  # Full sample data (used for all samples)
        K = K,                 # Full relationship matrix (for all samples)
        formula = formula,
        return_all_BLUP_values = TRUE,
        parallel = TRUE
      )
      
      # Store training predictions
      predictions_train[!is.na(Y_train)] <- blup_results$predictions[!is.na(Y_train)]
      
      # For each trait, extract test predictions and calculate accuracy
      for(t in 1:n_traits) {
        test_ids <- test_samples[[t]]
        if(length(test_ids) > 0) {
          # Store test predictions (already computed by run_rrBLUP_baseline)
          predictions[test_ids, t] <- blup_results$predictions[test_ids, t]
          
          # Calculate accuracy
          observed <- Y[test_ids, t]
          predicted <- predictions[test_ids, t]
          valid <- !is.na(observed) & !is.na(predicted)
          
          if(sum(valid) >= 5) {
            accuracies[f, t] <- cor(observed[valid], predicted[valid])
            cat(sprintf("  Trait %d (%s): Accuracy = %.4f (n=%d)\n", 
                        t, trait_names[t], accuracies[f, t], sum(valid)))
          } else {
            cat(sprintf("  Trait %d (%s): Not enough valid data points\n", 
                        t, trait_names[t]))
          }
        }
      }
      
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    })
  }
  
  # Calculate overall statistics
  trait_accuracies <- colMeans(accuracies, na.rm = TRUE)
  overall_accuracy <- mean(trait_accuracies, na.rm = TRUE)
  
  cat("\nOverall Results:\n")
  for(t in 1:n_traits) {
    cat(sprintf("  %s: %.4f\n", trait_names[t], trait_accuracies[t]))
  }
  cat(sprintf("  Average: %.4f\n", overall_accuracy))
  
  # Return results
  return(list(
    predictions = predictions,
    predictions_train = predictions_train,
    fold_accuracies = accuracies,
    trait_accuracies = trait_accuracies,
    overall_accuracy = overall_accuracy
  ))
}