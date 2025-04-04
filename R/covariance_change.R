#' Covariance Change
#'
#' This method implements a method for detection of covariance changes in
#'  functional data.
#'
#' @param X Numeric data.frame of functional data observations--rows for
#'    evaluated values and columns indicating FD
#' @param kappa (Optional) Numeric used for weighting and such. Default is 1/4
#' @param len (Optional) Numeric for window/repetitions for covariance change.
#'    Default is 30.
#'
#' @return Location of change. NA for no change and numeric if there is a change
#' @export
#'
#' @references Change point analysis of covariance functions: a weighted
#'  cumulative sum approach, L. Horvath, G. Rice, Y. Zhao (2022) Journal of
#'  Multivariate Analysis.
cov_change <- function(X, kappa = 1 / 4, len = 30) {
  stat_d0 <- .weight_TNstat(X, kappa = kappa)
  cv_d0 <- .weight_criticalvalueMC(X, len = len, kappa = kappa)

  if (stat_d0[[1]] > cv_d0[2]) {
    return(stat_d0[[2]])
  }

  return(NA)
}

#' L2 Norm
#'
#' This (internal) function computes the L2 norm of the data.
#'
#' See use in tau_est
#'
#' @param vec Vector of numerics to compute L2 norm
#'
#' @return Numeric L2 norm value
#'
#' @noRd
.l2norm <- function(vec) {
  return(sqrt(sum(vec^2)))
}


#' Compute Zn Statistic for Functional Covariance Changes Under Null
#'
#' This (internal) function implements Theorem 2.1 / 2.2.
#'
#' @param xdm Data.frame of numerics. It is the finite realization of DEMEAN-ed
#'  functional time series data, where curves are stored in columns.
#' @param u Numeric in \eqn{(0, 1)}. That is, a fraction index over the interval
#' \eqn{[0, 1]}.
#'
#' @return Numeric Zn statistic
#'
#' @noRd
.ZNstat <- function(xdm, u) {
  grid_point <- nrow(xdm)
  N <- ncol(xdm)
  k <- floor(N * u)

  prek <- matrix(rowSums(apply(
    as.matrix(xdm[, 1:k]), 2,
    function(x) {
      x %o% x
    }
  )), grid_point, grid_point)
  fullk <- matrix(rowSums(apply(
    as.matrix(xdm), 2,
    function(x) {
      x %o% x
    }
  )), grid_point, grid_point)
  ZNu <- N^(-1 / 2) * (prek - (k / N) * fullk)

  ZNu
}


#' Compute Zn Statistic for Functional Covariance Changes Under Change
#'
#' @inheritParams .ZNstat
#'
#' @return Numeric Zn statistic
#'
#' @noRd
.ZNstat_cp <- function(xdm, u) {
  grid_point <- nrow(xdm)
  N <- ncol(xdm)
  k <- floor(N * u)

  prek <- matrix(rowSums(apply(as.matrix(xdm[, 1:k]), 2, function(x) {
    x %o% x
  })), grid_point, grid_point)
  fullk <- matrix(rowSums(apply(as.matrix(xdm), 2, function(x) {
    x %o% x
  })), grid_point, grid_point)
  ZNu <- (prek - (k / N) * fullk)

  ZNu
}


#' Compute the Weighted Tn Statistic
#'
#' The function computes the weighted Tn statistic, introduced after Theorem 2.3
#'
#' @param xf Data.frame of numerics. Finite realization of functional time
#'    series data, where curves are stored in columns.
#' @param kappa Numeric used for weighting and such
#'
#' @return List of two values: (stat) giving the weighted Tn statistic and
#'    (changepoint) numeric for change location
#'
#' @noRd
.weight_TNstat <- function(xf, kappa) {
  # .int_approx_tensor <- function(x) { # x is a 4-dimensional tensor
  #   dt <- length(dim(x))
  #   temp_n <- nrow(x)
  #   return(sum(x) / (temp_n^dt))
  # }

  grid_point <- nrow(xf)
  N <- ncol(xf)
  xdm <- apply(xf, 2, function(x, xmean) {
    x - xmean
  }, xmean = rowMeans(xf))
  uind <- seq(0, 1, length = N + 1)[2:(N + 1)]
  zn2 <- list()
  zn_cp <- c(rep(0, N))
  for (i in 1:(N - 1)) {
    zn2[[i]] <- (.ZNstat(xdm, uind[i]))^2 /
      ((uind[i] * (1 - uind[i]))^(2 * kappa)) ### kappa = 1/4

    zn_cp[i] <- (N / (i * (N - i)))^(kappa) *
      .int_approx_tensor((.ZNstat_cp(xdm, uind[i]))^2)
  }
  inm <- Reduce(`+`, zn2) / N
  stat <- (1 / grid_point)^2 * sum(inm)

  mcp <- max(zn_cp[(0.1 * N):(0.9 * N)])
  changepoint <- which(zn_cp == mcp)

  return(list(stat, changepoint))
}


## Critical values

#' Long-run covariance estimator using tensor operator
#'
#' @param dat An array with dimension (grid_point,grid_point,N)
#'
#' @return Numeric matrix for long-run covariance
#'
#' @noRd
long_run_covariance_4tensor <- function(dat) {
  grid_point <- dim(dat)[1]
  Tval <- dim(dat)[3]
  datmean <- apply(dat, c(1, 2), mean)
  center_dat <- sweep(dat, 1:2, datmean)

  .cov_l <- function(band, nval) {
    cov_sum <- .gamma_l(0, nval)

    for (ik in 1:(nval - 1)) {
      cov_sum <- cov_sum + .kweights(ik / band, kernel = "Bartlett") * (2 * .gamma_l(ik, nval)) # + .gamma_l(ik,nval))    ##aperm(.gamma_l(ik,nval),c(2,1,3,4)))
    }
    return(cov_sum)
  }

  .gamma_l <- function(lag, Tval) {
    gamma_lag_sum <- 0
    if (lag >= 0) {
      for (ij in 1:(Tval - lag)) {
        gamma_lag_sum <- gamma_lag_sum + center_dat[, , ij] %o% center_dat[, , (ij + lag)]
      }
    } else {
      for (ij in 1:(Tval + lag)) {
        gamma_lag_sum <- gamma_lag_sum + center_dat[, , (ij - lag)] %o% center_dat[, ij]
      }
    }
    return(gamma_lag_sum / (Tval - lag))
  }

  hat_h_opt <- Tval^(1 / 4)
  lr_covop <- .cov_l(band = hat_h_opt, nval = Tval)

  return(lr_covop)
}

#' Kernels for applying weights
#'
#' This (internal) function computes the kernels and related weights
#'
#' @param x <Add Information>
#' @param kernel String indicating kernel to use.
#' @param normalize (Optional) Boolean indicating if x should be normalized.
#'    Default is FALSE.
#'
#' @return Numeric for k weight
#'
#' @noRd
.kweights <- function(x,
                      kernel = c(
                        "Truncated", "Bartlett", "Parzen", "Tukey-Hanning",
                        "Quadratic Spectral"
                      ),
                      normalize = FALSE) {
  kernel <- match.arg(kernel)
  if (normalize) {
    ca <- switch(kernel,
                 Truncated = 2,
                 Bartlett = 2 / 3,
                 Parzen = 0.539285,
                 `Tukey-Hanning` = 3 / 4,
                 `Quadratic Spectral` = 1
    )
  } else {
    ca <- 1
  }
  switch(kernel,
         Truncated = {
           ifelse(ca * x > 1, 0, 1)
         },
         Bartlett = {
           ifelse(ca * x > 1, 0, 1 - abs(ca * x))
         },
         Parzen = {
           ifelse(ca * x > 1, 0, ifelse(ca * x < 0.5, 1 - 6 * (ca *
                                                                 x)^2 + 6 * abs(ca * x)^3, 2 * (1 - abs(ca * x))^3))
         },
         `Tukey-Hanning` = {
           ifelse(ca * x > 1, 0, (1 + cos(pi * ca * x)) / 2)
         },
         `Quadratic Spectral` = {
           y <- 6 * pi * x / 5
           ifelse(x < 1e-04, 1, 3 * (1 / y)^2 * (sin(y) / y - cos(y)))
         }
  )
}


#' Compute critical values for .weight_TNstat ( \eqn{T_N(\kappa)} )
#'
#' This (internal) function computes the critical values for weighted Tn statistic
#'
#' @inheritParams .weight_TNstat
#' @param len Numeric for window/repetitions for covariance change.
#'
#' @return Numeric critical values (0.9, 0.95, 0.99)
#'
#' @noRd
.weight_criticalvalueMC <- function(xf, len, kappa) {
  grid_point <- nrow(xf)
  N <- ncol(xf)

  ## cov weight function
  times <- 1:grid_point / grid_point
  wmat <- matrix(NA, grid_point - 2, grid_point - 2)
  for (i in 2:(grid_point - 1)) {
    for (j in 2:(grid_point - 1)) {
      wmat[i - 1, j - 1] <- (min(times[i], times[j]) - times[i] * times[j]) /
        ((times[i] * (1 - times[i]))^kappa * (times[j] * (1 - times[j]))^kappa)
    }
  }
  weig <- as.vector(svd(wmat / grid_point)$d)

  ## cov operators
  rref <- stats::runif(len, 0, 1)
  rref <- c(sort(rref), 1)
  rrefind <- round(rref * dim(xf)[1])
  rrefind[which(rrefind == 0)] <- 1
  xfMC <- xf[rrefind, ]

  xdm <- apply(xfMC, 2, function(x, xmean) {
    x - xmean
  }, xmean = rowMeans(xfMC))
  zi <- zm <- array(0, c((len + 1), (len + 1), N))
  for (i in 1:N) {
    zi[, , i] <- xdm[, i] %o% xdm[, i]
  }
  zimean <- apply(zi, c(1, 2), mean)
  for (i in 1:N) {
    zm[, , i] <- zi[, , i] - zimean
  }

  lrcov <- long_run_covariance_4tensor(zm)
  lrcov <- tensorA::as.tensor(round(lrcov / (len + 1)^2, 6) )
  eigvals <- tensorA::svd.tensor(lrcov, c(3, 4), by = "e")
  eigmat <- as.vector(eigvals$d)

  lim_sum <- 0
  for (ell in 1:length(eigmat)) {
    klim <- 0
    for (k in 1:length(weig)) {
      Nm <- stats::rnorm(2000, mean = 0, sd = 1)
      klim <- klim + eigmat[ell] * weig[k] * Nm^2
    }
    lim_sum <- lim_sum + klim
  }

  cv <- stats::quantile(lim_sum, probs = c(0.90, 0.95, 0.99))

  return(cv)
}


#' Approximate Integral of Tensor
#'
#' This (internal) function using a Riemann sum to approximate the integral
#'
#' @param x  4-dimensional tensor
#'
#' @return Approximate integral value
#'
#' @noRd
.int_approx_tensor <- function(x) {
  dt <- length(dim(x))
  temp_n <- nrow(x)

  return(sum(x) / (temp_n^dt))
}
