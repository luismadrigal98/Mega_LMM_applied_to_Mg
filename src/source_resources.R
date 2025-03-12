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