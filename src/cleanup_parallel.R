cleanup_parallel <- function() {
  if(!is.null(env$cluster)) {
    parallel::stopCluster(env$cluster)
    future::plan(future::sequential)
    message("Parallel backends have been shut down")
  }
}