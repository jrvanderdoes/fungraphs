
#' Cumulative Intraday Returns
#'
#' Compute the cumulative intraday returns.
#'
#' @param dat Functional data as a matrix with the columns as individual
#'  observations.
#'
#' @return Functional data for the CIDRs
#' @export
#'
#' @references ice, G., Wirjanto, T., & Zhao, Y. (2019). Tests for conditional
#'  heteroscedasticity with functional data and goodness-of-fit tests for FGARCH
#'  models. IDEAS Working Paper Series from RePEc.
#'
#' @examples
#' compute_cidr(matrix(runif(5000),ncol=100))
compute_cidr <- function(dat){
  dat_cidr <- dat
  for(i in 1:nrow(dat_cidr)){
    dat_cidr[i,] <- 100*(log(dat[i,]) - log(dat[1,]))
  }
  dat_cidr
}
