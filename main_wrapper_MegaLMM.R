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
## 1) Setting up the environment ----
## _____________________________________________________________________________

# 1.1) Load the required libraries ----

required_libraries <- c("MegaLMM", "ggplot2")

set_environment(required_pckgs = required_libraries)

## *****************************************************************************
## 2) Setting work directory ----
## _____________________________________________________________________________

setwd('./')

## *****************************************************************************
## 2) Importing the data ----
## _____________________________________________________________________________

#' Data that should be imported for this wrapper includes trait values, genetic
#' relationships, and covariates. User can omit covariates in this implementation

##