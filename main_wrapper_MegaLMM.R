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
                        "parallel", "doParallel", "matrixcalc")

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

## 3.1) Checking the data----

#' K must be positive semi-definite

if(!matrixcalc::is.positive.semi.definite(x = Data@K)) {
  Data@K <- fix_K_by_eigenvalues(Data@K, max_iter = 3, verbose = TRUE)$matrix
}

## *****************************************************************************
## 4) Running BLUP using cross validation to assess model performance ----
## _____________________________________________________________________________

## 4.1) Creating the fold assignation per trait ----

#' By using different fold assignation per trait we ensure robust statistical
#' individuals are used for training and testing in each traitrformance

fold_matrix <- create_cv_folds(Y = Data@Y, sample_data = Data@C_data, 
                               k = 5, group_col = "Covariate_2", 
                               id_col = "Ind_ID")

results_BLUP_CV <- run_CV_rrBLUP(
  Y = Data@Y, 
  fold_matrix = fold_matrix,
  sample_data = Data@C_data,
  K = Data@K,
  formula = ~ Covariate_1 + Covariate_2 + Covariate_3 + Covariate_3
)

# Create a directory if it doesn't exist
dir.create("./results/BLUP", recursive = TRUE, showWarnings = FALSE)

# Format results into a data frame
results_summary <- data.frame(
  Trait = names(results_BLUP_CV$trait_accuracies),
  Accuracy = results_BLUP_CV$trait_accuracies
)
# Add overall average
results_summary <- rbind(
  results_summary,
  data.frame(Trait = "Overall_Average", 
             Accuracy = results_BLUP_CV$overall_accuracy)
)

# Write as a single file
write.table(results_summary,
            file = "./results/BLUP/BLUP_CV_performance_all_traits.txt",
            sep = '\t', row.names = FALSE, quote = FALSE)

## Fitting a BLUP model using all the data ----

BLUP_result_all <- run_rrBLUP_baseline(
  Y = Data@Y, 
  sample_data = Data@C_data,
  K = Data@K,
  formula = ~ Covariate_1 + Covariate_2 + Covariate_3 + Covariate_3,
  return_all_BLUP_values = TRUE,
  parallel = TRUE
)

## Exporting the results
write.table(BLUP_result_all$predictions, 
            file = "./results/BLUP/BLUP_predicitons_all_traits.txt", 
            sep = '\t')

## *****************************************************************************
## 5) Setting up the MEGA_LMM run ----
## _____________________________________________________________________________

## Setting the MegaLMM state with all the data ----

MegaLMM_state <- setup_MegaLMM(
  Y = Data@Y,
  formula = ~ Covariate_1 + Covariate_2 + Covariate_3 + Covariate_4 + (1|Ind_ID),
  sample_data = Data@Design_dt,
  kinship_matrix = Data@K,  # Changed from K to kinship_matrix
  latent_factors = 7,       # Changed from factors to latent_factors
  run_ID = "MegaLMM_7_latent"
)

## Setting the priors for all the data ----
MegaLMM_state <- set_MegaLMM_priors(
  MegaLMM_state,
  lambda_prior_type = "horseshoe"
)

## Optimize the processing when missing data is present ----
if(any(is.na(Data@Y))){
  MegaLMM_state <- optimize_missing_data(MegaLMM_state)
}

## Initialize the MegaLMM instance ----

MegaLMM_state <- init_MegaLMM(MegaLMM_state, retain_par = c('Lambda', 'F_h2', 'resid_h2', 'tot_Eta_prec', 'B1'),
                                    retain_mean_par = c('Eta_mean'),
                                    retain_function = list(
                                      U = 'U_F %*% Lambda + U_R + X1 %*% B1',
                                      G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
                                      R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])',
                                      h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'
                                    ),
                                    verbose = TRUE)

## Run burnin ----

MegaLMM_state <- run_MegaLMM_burnin(
  MegaLMM_state,
  burnin_rounds = 3,
  iterations_per_round = 100,
  make_plots = TRUE
)

## Cleaning-up environment ----
cleanup_parallel()