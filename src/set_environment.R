set_environment <- function(required_pckgs,
                            automatic_download = FALSE,
                            personal_seed = as.numeric(Sys.time()),
                            parallel_backend = FALSE)
{
  #' This function will set up the working environment for performing all the
  #' analysis. It will load all the required packages, set the seed for
  #' reproducibility and set the parallel backend if required.
  #' 
  #' @param required_pckgs A character vector of the required packages.
  #' 
  #' @param automatic_download A logical value to set the automatic download of
  #' the required packages. Default set to FALSE.
  #' 
  #' @param personal_seed An integer to set the seed for reproducibility. Default
  #' set to system time.
  #' 
  #' @param parallel_backend A logical value to set the parallel backend. Default
  #' set to FALSE.
  #' 
  #' @return invisible
  #' ___________________________________________________________________________
  
  # Initialize return object
  env_objects <- list()
  
  ## Loading the required libraries ----
  
  message("Loading the required libraries")
  
  tryCatch(
    {
      for(pckg in required_pckgs)
      {
        
        if(!require(pckg, character.only = TRUE))
        {
          if (automatic_download == TRUE)
          {
            install.packages(pckg)
            library(pckg, character.only = TRUE)
          }
          else
          {
            message(paste0("The automatic download of the required package ", 
                           pckg,  " is disabled"))
            message("Install it manually and run the script again")
            stop()
          }
        }
        else
        {
          library(pckg, character.only = TRUE)
        }
      }
      
      message("The required libraries have been loaded")
    }, error = function(e) {
      message("Some packages cannot be installed through CRAN, trying with remotes")
      print(e)
    })
  
  ## Setting the seed ----
  set.seed(personal_seed)
  
  ## Setting the parallel backend ----
  if(parallel_backend == TRUE) {
    # Check dependencies first
    pkgs_needed <- c("future", "doFuture", "foreach", "doParallel", "parallelly")
    missing_pkgs <- pkgs_needed[!sapply(pkgs_needed, requireNamespace, quietly = TRUE)]
    if(length(missing_pkgs) > 0) {
      stop("Missing required parallel packages: ", paste(missing_pkgs, collapse=", "))
    }
    
    # Set up future
    options(doFuture.rng.onMisuse = "ignore")
    registerDoFuture()
    future::plan(future::multisession, workers = parallelly::availableCores() - 1)
    
    # Set up foreach
    cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    # Store cluster in return object
    env_objects$cluster <- cl
  }
  
  return(env_objects)
}