% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_fd}
\alias{plot_fd}
\title{Plot functional data}
\usage{
plot_fd(
  data,
  curve_points = 1:nrow(data),
  plot_title = NULL,
  val_axis_title = "",
  res_axis_title = "",
  FD_axis_title = "",
  FDReps = as.Date(colnames(data), "\%Y_\%m_\%d"),
  eye = list(x = -0.5, y = -1.5, z = 0.5),
  aspectratio = list(x = 1, y = 0.75, z = 0.5),
  showticklabels = TRUE
)
}
\arguments{
\item{data}{Numeric data.frame with evaled points on rows and fd objects in columns.}

\item{curve_points}{Evaluation points.}

\item{plot_title, val_axis_title, res_axis_title, FD_axis_title}{Titles for use
in the plot.}

\item{FDReps}{Point to label individual curves.}

\item{eye}{Angle to view plot.}

\item{aspectratio}{Relative size of each dimension.}

\item{showticklabels}{Boolean if ticks should be shown on plot.}
}
\value{
A plotly object.
}
\description{
Function to display functional data. See vastly improved visualization in the
package fChange (\url{https://github.com/jrvanderdoes/fChange}).
}
