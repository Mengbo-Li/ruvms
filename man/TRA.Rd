% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TRA.R
\name{TRA}
\alias{TRA}
\title{Technical replicate agreement (TRA) plot}
\usage{
TRA(
  data,
  replicates = NULL,
  M = NULL,
  col = NULL,
  ylim = NULL,
  printLabel = TRUE,
  printLegend = FALSE,
  guides = NULL,
  title = NULL
)
}
\arguments{
\item{data}{The m-by-n data matrix with samples on the rows and measurements (
proteins etc.) on the columns. Missing data (NAs) are allowed.}

\item{replicates}{A vector indicating which samples are replicates of each
other. Obsolete when \code{M} is not \code{NULL}.}

\item{M}{The replication structure matrix.}

\item{col}{Color of the TRA boxes. Default is \code{NULL}.}

\item{ylim}{The y-axis limits of the output plot.}

\item{printLabel}{Logical. Whether to print the names of the technical
replicates.}

\item{printLegend}{Logical. Whether to print the color legend.}

\item{guides}{Numeric. Horizontal dashed lines for reference.}

\item{title}{Character. Title of the output plot.}
}
\value{
The TRA plot.
}
\description{
Technical replicate agreement (TRA) plots provide a visualisation diagnostic
that measures the closeness among technical replicates.
}
\details{
The assumption of TRA plots is that if all unwanted variations are
effectively removed, replications of a unique effective sample are
identical after normalisation, even if we do not declare them as replicates
of each other
}
\references{
Molania, R., Gagnon-Bartsch, J. A., Dobrovic, A., & Speed, T. P. (2019).
A new normalization for Nanostring nCounter gene expression data.
\emph{Nucleic acids research}, 47(12), 6073-6083.
}
