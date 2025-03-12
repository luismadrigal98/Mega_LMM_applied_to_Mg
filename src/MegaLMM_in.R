#' This script contains the definition of the S4 class MegaLMM_in.
#' 
#' This will ease the manipulation of the inputs and the overall work with the 
#' package.
#' 
#' _____________________________________________________________________________

setClass(Class = "MegaLMM_in", slots = c(Y_data = "data.frame", 
                                         Y = "matrix",
                                         K = "matrix",
                                         C_data = "data.frame",
                                         C = "matrix"))