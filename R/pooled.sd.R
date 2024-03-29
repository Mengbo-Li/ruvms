#' Calculate pooled standard deviations for each protein with missing data
#'
#' The \code{pooled.sd} function calculates the pooled standard deviation in
#' each protein within technical replicates while accommodating possible missing
#' data. When there is not enough replicates to calculate pooled standard
#' deviations (number of non-missing samples is less than or equal to number of
#' groups for example), a simple sample standard deviation is used.
#' This step is needed when \code{standardise} is enabled in \code{ruvms}.
#'
#' @param Y The m-by-n data matrix with samples on the rows and measurements
#' (proteins etc. ) on the columns. Missing data (NAs) are allowed.
#'
#' @param M The replication structure. See \code{ruv::RUVIII} for details.
#'
#' @return A numeric vector of pooled standard deviation of each measurement.
#'
#' @export
pooled.sd <- function(Y, M) {
   if(!requireNamespace("Rfast", quietly = TRUE))
      stop("Rfast package required but is not installed (or can't be loaded)")
   group.df <- (t(M) %*% !is.na(Y)) - 1
   group.df <- ifelse(group.df < 0, 0, group.df)
   group.var <- do.call(cbind, lapply(1:ncol(M), function(i) Rfast::colVars(Y[M[, i] == 1, , drop = FALSE], na.rm = TRUE)))
   group.var[is.na(group.var)] <- 0
   var.j <- diag(group.var %*% group.df / colSums(group.df))
   if (any(is.na(var.j)))
      var.j[is.na(var.j)] <- Rfast::colVars(Y[, is.na(var.j), drop = FALSE], na.rm = TRUE)
   sqrt(var.j)
}
