% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphs.R
\name{create_graph}
\alias{create_graph}
\title{Create Graph for Data}
\usage{
create_graph(
  data,
  K = 1,
  type = c("MDT", "MDP", "NNL"),
  errType = c("L2", "L1", "Lp")
)
}
\arguments{
\item{data}{Numeric data.frame with evaled points on rows and fd objects in columns}

\item{K}{Number of orthogonal trees}

\item{type}{Type of graph: "MDT", "MDP", or "NNL"}

\item{errType}{Type of errors: "L2", "L1", or "Lp"}
}
\value{
Graph for use in graph_segmentation
}
\description{
Create Graph for Data
}
