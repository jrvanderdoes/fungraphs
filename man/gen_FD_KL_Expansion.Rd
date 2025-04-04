% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmp_gen_fd.R
\name{gen_FD_KL_Expansion}
\alias{gen_FD_KL_Expansion}
\title{FD Generation - KL Expansion}
\usage{
gen_FD_KL_Expansion(
  ns,
  eigsList,
  basesList,
  meansList,
  distsArray,
  evals,
  kappasArray = c(0),
  burnin = 100,
  silent = F,
  dof = NULL,
  shape = NULL,
  ...
)
}
\arguments{
\item{ns}{Vector of Numerics. Each value in N is the number of observations
for a given group.}

\item{eigsList}{Vector of eigenvalues. Length 1 or m.}

\item{basesList}{A list of bases (eigenfunctions), length m.}

\item{meansList}{A vector of means, length 1 or N.}

\item{distsArray}{A vector of distributions, length 1 or m.}

\item{evals}{A vector of points indicating the points to evaluate the
functions on.}

\item{kappasArray}{Numeric [0,1] indicating strength of VAR(1) process.}

\item{burnin}{A numeric value indicating the number of burnin trials.}

\item{silent}{A Boolean that toggles running output.}

\item{dof}{Numeric for degrees of freedom with t-distribution.}

\item{shape}{Numeric for shape of gamma distribution.}

\item{...}{Additional parameters to send in. Unused.}
}
\value{
List with (1) data (N-by-m) and (2) previous errors.
}
\description{
Generation of FD using an autoregressive karhunen-loeve expansion. This can
include change points in any combination of the following:
\itemize{
\item Mean
\item Distribution
\item Eigenvalue(s)
\item Eigenvector(s)
}
In this sense, the function creates m 'groups' of discretely observed
functions with similar properties. See updated
version in the package fChange (\url{https://github.com/jrvanderdoes/fChange}).
}
\examples{
# result <- gen_FD_KL_Expansion(
#   N=c(100,50,25),
#   eigenvalues = list(rep(1,5),
#                      c(1/sqrt(1:5)),
#                      c(1/sqrt(1:5))),
#
#   parameters =list('bmotion'=list('N'=100, 'process'='bmotion', 'sd'=1),
#                    'bbridge'=list('N'=100, 'process'='bbridge', 'sd'=1),
#                    'kl'=list('process'='kl', 'N'=100,
#                              'distribution'='Normal',
#                              'eigenvalues'=1/1:4,
#                              'mean'=0, 'dependence'=0,
#                              'basis'=fda::create.bspline.basis(),
#                              'sd'=1),
#                    'ou'=list('N'=100, 'process'='ou', 'dependence'=0 ) ,
#                    'far1'=list('N'=100, 'process'='far1', 'dependence'=0,
#                                'sd'=1,'vary'=FALSE) )
#                  )
}
