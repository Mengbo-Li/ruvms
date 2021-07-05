#' Remove unwanted variation in mass-spectrometry data
#'
#' The `ruvms` function performs an `ruv::RUVIII` like normalisation on
#' mass-spectrometry data while allowing missing values (NAs) in input data.
#'
#' @param Y The data. An m-by-n matrix, where m is the number of observations
#' and n is the number of measurements.
#'
#' @param M The replication structure matrix. See \code{ruv::RUVIII} for
#' details.
#'
#' @param ctl A logical vector of length n indication the negative control
#' measurements (peptides, proteins or metabolites).
#'
#' @param k1 Tuning parameter. Default is 0.
#'
#' @param lambda Tuning parameter. Default is 1e-5, only needed when
#' \code{standardise} is \code{FALSE}.
#'
#' @param eta Column-wise (e.g., protein-wise) covariates. See \code{ruv::RUVI}
#' for more.
#'
#' @param include.intercept Logical. See \code{ruv::RUVI} for more.
#'
#' @param average Logical, whether to return averaged replicates after
#' normalisation.
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
#' \code{delta} \tab \tab \tab \tab \tab The delta tuning parameter. \cr
#' \code{lambda} \tab \tab \tab \tab \tab The lambda tuning parameter. \cr
#' }
#'
#' @seealso
#' \code{\link[ruv]{RUV1}}, \code{\link[ruv]{RUVIII}}
#'
#' @examples
#'
#' # See vignettes.
#'
#' @export
ruvms <- function (Y,
                   M,
                   ctl,
                   k1 = 0,
                   lambda = 1e-5,
                   standardise = TRUE,
                   delta = NULL,
                   eta = NULL,
                   include.intercept = TRUE,
                   average = FALSE,
                   return.info = FALSE,
                   inputcheck = FALSE) {

   if (is.data.frame(Y)) Y <- data.matrix(Y)
   m <- nrow(Y)
   n <- ncol(Y)
   ctl <- ctl & colSums(is.na(Y)) < 1

   if (inputcheck) {
      if (sum(ctl) < 2)
         stop("Too few complete negative control measurements. At least two are
           needed. ")
      if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 0)
         stop("Y contains infinities. This is not supported.")
      if (m > n)
         warning("m is greater than n! This is not a problem itself,
              but may indicate that you need to transpose your data matrix.
              Please ensure that rows correspond to observations
              (e.g. sample ID) and columns correspond to features
              (e.g. proteins).")
   }

   if (!is.null(eta)) {
      adjustments <- Y - RUV1(Y, eta, ctl, include.intercept = include.intercept)
      Y <- RUV1(Y, eta, ctl, include.intercept = include.intercept)
   }

   if (standardise) {
      mu.j <- colMeans(Y, na.rm = TRUE)
      Y <- sweep(Y, 2, mu.j, "-")
      sd.j <- pooled.sd(Y, M)
      Y <- sweep(Y, 2, sd.j, "/")
      Y[is.nan(Y)] <- 0
      lambda <- 0
   }

   if (ncol(M) < m) {
      decomp <- eigen(Y[, ctl] %*% t(Y[, ctl]))
      di <- decomp$values
      if (!standardise) delta <- 0
      if (standardise & min(di) > 0) delta <- 0
      if (standardise & min(di) <= 0 & is.null(delta))
         delta <- 0.1^floor(-log10(-min(di)))
      U <- decomp$vectors
      k2 <- (1 - Matrix::rankMatrix(Y[, ctl])/2/m) * rankMatrix(Y[, ctl])
      hi <- sapply(1:length(di), function(i)
         I(i > k1) * min(1/(di[i] + delta) + lambda,
                         1/(di[k2] + delta) + lambda))

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

   if (standardise) {
      Y <- sweep(Y, 2, sd.j, "*")
      Y <- sweep(Y, 2, mu.j, "+")
   }

   if (average) {
      Y <- do.call(rbind, lapply(1:ncol(M),
                                 function(i) colMeans(Y[M[, i] == 1, ,
                                                        drop = FALSE],
                                                      na.rm = TRUE)))
      rownames(Y) <- colnames(M)
   }

   if (return.info) return(list(newY = Y, di = di, hi = hi,
                                delta = delta, lambda = lambda)) else return(Y)
}
