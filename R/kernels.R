#' Bartlett Kernel Functions
#'
#' Kernel where \eqn{max(0,1-|x|), h\neq 0}. If
#'  \eqn{x=0/0} then the value \eqn{1} is given.
#'
#' @param x Numeric value(s) at which to evaluate kernel.
#'  It indicates current lag divided by window
#'
#' @name kernels
#'
#' @references Horvath, L., & Rice, G. (2024). Change point analysis for time
#'  series (1st ed. 2024.). Springer Nature Switzerland.
#'
#' @references L. Horvath, P. Kokoszka, G. Rice (2014) "Testing stationarity of functional
#'  time series", Journal of Econometrics, 179(1), 66-82.
#'
#' @references Politis, D. N. (2003). Adaptive bandwidth choice. Journal of Nonparametric
#'  Statistics, 15(4-5), 517-533.
#'
#' @references Politis, D. N. (2011). Higher-order accurate, positive semidefinite
#'  estimation of large-sample covariance and spectral density matrices.
#'  Econometric Theory, 27(4), 703-744.
#'
#' @return Values from given lag(s) in the kernel.
#'
#' @export
#'
#' @examples
#' bartlett_kernel(-20:20/15)
bartlett_kernel <- function(x) {
  pmax(0, ifelse(is.nan(1 - abs(x)),1,1 - abs(x)) )
}
