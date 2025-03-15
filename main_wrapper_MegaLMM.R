#' @title MegaLMM implementation in Mimulus guttatus data
#' 
#' @description
#' This wrapper will implement MegaLMM utilities for fitting multi-trait linear
#' mixed models.
#' 
#' @author Luis Javier Madrigal-Roca & John K. Kelly
#' 
#' @date 03/12/2025
#' _____________________________________________________________________________

## *****************************************************************************
## 1) Setting work directory ----
## _____________________________________________________________________________

setwd('./')

## *****************************************************************************
## 2) Setting up the environment ----
## _____________________________________________________________________________

# 2.1) Load the source files ----

# Check that the src directory is contained in current working directory

source_resources <- function(path)
{
  if(dir.exists(path = path))
  {
    lapply(list.files(path, full.names = T), FUN = source)
  }
  else
  {
    stop("Directory specified for sourcing resources does not exist.")
  }
}

source_resources("src")

# 2.2) Load the required libraries ----

required_libraries <- c("MegaLMM", "ggplot2", "tidyr", "rrBLUP", "foreach",
                        "parallel", "doParallel")

env <- set_environment(required_pckgs = required_libraries, 
                       parallel_backend = TRUE)

## *****************************************************************************
## 3) Importing the data ----
## _____________________________________________________________________________

#' Data that should be imported for this wrapper includes trait values, genetic
#' relationships, and covariates. User can omit covariates in this implementation

Data <- Mega_LMM_input_reader(Y_path = "data/Toy_version_Mg/Pheno6.txt",
                              K_path = "data/Toy_version_Mg/RelationshipMatrix.txt",
                              C_path = "data/Toy_version_Mg/mim.fe.txt",
                              sep = "\t")

## *****************************************************************************
## 4) Running BLUP using cross validation to assess model performance ----
## _____________________________________________________________________________

## 4.1) Creating the fold assignation per trait ----

#' By using different fold assignation per trait we ensure 

fold_matrix <- create_cv_folds(Y = Data@Y, sample_data = Data@C_data, 
                               k = 5, group_col = "Covariate_2", 
                               id_col = "Ind_ID")



## Cleaning-up environment ----
cleanup_parallel()