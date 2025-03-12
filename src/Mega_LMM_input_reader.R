Mega_LMM_input_reader <- function(Y_path, K_path, C_path = NULL,
                                  sep = "\t")
{
  #' This function will read in all the data required for building a Mega_LMM
  #' object. Users must comply with the accepted input formats. We can make this 
  #' reader more flexible, but users can comply with format standards currently
  #' implemented easily.
  #' 
  #' @param Y_path A path to a file with the trait values. The first column contains the
  #' individuals or lines under analysis, and these are going to be transformed
  #' to the matrix representation by trimming out such metadata and including it
  #' as rownames. Notice that in the output, original data.frame and matrix re-
  #' presentations are retained.
  #' 
  #' @param K_path Path to relationship data. Usually is read in also as data frame. 
  #' Similar processing is done to get a matrix. No data.frame representation is 
  #' retained. User is reponsible of ensuring that the order of the entries here
  #' are the same that the order of entries in Y. In other words, Ind_ID order
  #' must be the same.
  #' 
  #' @param C_path NULL be default (no covariates are included). The user can especify
  #' a data source of covariates, that should in terms of individuals or lines
  #' analyzed the one provided as Y.
  #' 
  #' @param sep A character to specify the separator used in the input files.
  #' 
  #' @return A MegaLMM_in object.
  #' ___________________________________________________________________________
  
  # Create an instance of the container class
  Mega_object <- new("MegaLMM_in")
  
  # Read in the trait data
  Mega_object@Y_data <- read.table(file = Y_path,
                                   sep = sep)
  
  Ind_ID <- Mega_object@Y_data[, 1] ## First column contain the IDs
  
  if(all(grepl(pattern = "V", x = colnames(Mega_object@Y_data))))
  {
    # The file does not have header. We are going to assign the names.
    # First column refers to Ind_ID, the rest to traits
    
    col_m <- paste0("Trait_", 1:(length(colnames(Mega_object@Y_data))- 1)) # Colnames in matrix
    col_d <- c("Ind_ID", col_m)
    
    # Relabeling also the data.frame
    colnames(Mega_object@Y_data) <- col_d
  }
  
  else
  {
    # The file have header. We are going to assign the names.
    # First column refers to Ind_ID, the rest to traits
    
    col_d <- colnames(Mega_object@Y_data)
    col_m <- col_d[-1]
  }
  
  # Create the matrix version of Y, including the colnames and rownames
  Mega_object@Y <- as.matrix(Mega_object@Y_data[, -1])
  rownames(Mega_object@Y) <- Ind_ID
  colnames(Mega_object@Y) <- col_m
  
  # Read in the relationship data
  Mega_object@K <- read.table(file = K_path,
                              sep = sep) |>
    as.matrix()
  colnames(Mega_object@K) <- rownames(Mega_object@K) <- Ind_ID
  
  # Read in the covariates (if any)
  
  if (!is.null(C_path))
  {
    Mega_object@C_data <- read.table(file = C_path,
                                     sep = sep)
    
    if(all(grepl(pattern = "V", x = colnames(Mega_object@C_data))))
    {
      # The file does not have header. We are going to assign the names.
      # First column refers to Ind_ID, the rest to traits
      
      col_m <- paste0("Covariate_", 1:(length(colnames(Mega_object@C_data)) - 1)) # Colnames in matrix
      col_d <- c("Ind_ID", col_m)
      
      # Relabeling also the data.frame
      colnames(Mega_object@C_data) <- col_d
    }
    
    else
    {
      # The file have a header. We are going to assign the names.
      # First column refers to Ind_ID, the rest to traits
      
      col_d <- colnames(Mega_object@C_data)
      col_m <- col_d[-1]
    }
    
    # Create the matrix version of Y, including the colnames and rownames
    Mega_object@C <- as.matrix(Mega_object@C_data[, -1])
    rownames(Mega_object@C) <- Ind_ID
    colnames(Mega_object@C) <- col_m
  }
  
  return(Mega_object)
}