#' Remove unwanted variation in mass-spectrometry data
#'
#' The `ruvms` function performs an `ruv::RUVIII` like normalisation on
#' mass-spectrometry data while allowing missing values (NAs) in input data.
#'
#' @param Y The data. An m-by-n matrix, where m is the number of observations
#' and n is the number of measurements.
#'
#' @param M The replication structure matrix.
#'
#' @param ctl A logical vector of length n indication the negative control
#' measurements.
#'
#' @param k1
#'
#' @param lambda
#'
#' @param eta Column-wise (e.g., protein-wise) covariates.
#'
#' @param include.intercept Logical.
#'
#' @param average Logical.
#'
#' @param return.info Logical. If \code{FALSE}, only the adjusted data matrix
#' is returned. If \code{TRUE}, additional information is returned (see below).
#'
#' @param input.check Logical. A basic sanity check on the inputs.
#'
#' @return
#' If \code{input.check = FASLE}, the normalised data matrix is
#' returned. Otherwise a list is returned which contains \tabular{llllll}{
#' \code{newY} \tab \tab \tab \tab \tab The normalised data matrix. \cr
#' \code{di} \tab \tab \tab \tab \tab Eigenvalues of the N matrix. \cr
#' \code{hi} \tab \tab \tab \tab \tab hi values in the generalised averging
#' operator. \cr
#' }
#'
#' @seealso
#' \code{\link[ruv]{RUV1}}, \code{\link[ruv]{RUVIII}}
#'
#' @examples
#' ## Not run
#'
#' @export
ruvms <- function (Y, M, ctl, k1 = 0, lambda = 1e-5, eta = NULL,
                   include.intercept = TRUE, average = FALSE,
                   return.info = FALSE, input.check = FALSE) {

   if (is.data.frame(Y)) Y <- data.matrix(Y)
   m <- nrow(Y)
   n <- ncol(Y)
   ctl <- ctl & colSums(is.na(Y)) < 1

   if (input.check) {
      if (m > n)
         warning("m is greater than n! This is not a problem itself,
              but may indicate that you need to transpose your data matrix.
              Please ensure that rows correspond to observations
              (e.g. sample ID) and columns correspond to features
              (e.g. proteins).")
      if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 0)
         stop("Y contains infinities. This is not supported.")
   }

   if (!is.null(eta)) {
      adjustments <- Y - RUV1(Y, eta, ctl, include.intercept = include.intercept)
      Y <- RUV1(Y, eta, ctl, include.intercept = include.intercept)
   }

   if (ncol(M) < m) {
      decomp <- eigen(Y[, ctl] %*% t(Y[, ctl]))
      di <- decomp$values
      U <- decomp$vectors
      k2 <- (1 - rankMatrix(Y[, ctl])/2/m) * rankMatrix(Y[, ctl])
      hi <- sapply(1:length(di), function(i)
         I(i > k1) * min(1/di[i] + lambda, 1/di[k2] + lambda))

      for (j in 1:ncol(Y)) {
         keep <- !is.na(Y[, j])
         Mj <- M[keep, , drop = FALSE]
         Mj <- Mj[, colSums(Mj) > 0, drop = FALSE]
         Uj <- U[keep, , drop = FALSE]
         Y[keep, j] <- Mj %*%
            solve(t(Mj) %*% Uj %*% diag(hi) %*% t(Uj) %*% Mj) %*%
            t(Mj) %*% Uj %*% diag(hi) %*% t(Uj) %*% Y[keep, j]
      }
   } else {stop("Not enough replications!")}

   if (!is.null(eta))
      Y <- Y + adjustments

   if (average)
      Y <- t(t(Y) %*% M / colSums(M))

   if (return.info) return(list(newY = Y, di = di, hi = hi)) else return(Y)
}
