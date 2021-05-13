ruvms <- function (Y,
                   M,
                   ctl,
                   k1 = 0,
                   lambda = 1e-5,
                   eta = NULL,
                   include.intercept = TRUE,
                   average = FALSE,
                   return.info = FALSE,
                   inputcheck = FALSE) {

   if (inputcheck) {
      if (m > n)
         warning("m is greater than n! This is not a problem itself, but may indicate that you need to transpose your data matrix. Please ensure that rows correspond to observations (e.g. sample ID) and columns correspond to features (e.g. proteins).")
      if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 0)
         stop("Y contains infinities. This is not supported.")
   }
   if (is.data.frame(Y)) { Y <- data.matrix(Y)  }

   m <- nrow(Y)
   n <- ncol(Y)
   ctl <- ctl & apply(Y, 2, noNA)

   if (!is.null(eta)) {
      adjustments <- Y - RUV1(Y, eta, ctl, include.intercept = include.intercept)
      Y <- RUV1(Y, eta, ctl, include.intercept = include.intercept)
   }

   if (ncol(M) < m) {

      decomp <- eigen(Y[, ctl] %*% t(Y[, ctl]))
      di <- decomp$values
      U <- decomp$vectors
      k2 <- (1 - rankMatrix(Y[, ctl])/2/m) * rankMatrix(Y[, ctl])
      hi <- sapply(1:length(di), function(i) I(i > k1) * min(1/di[i] + lambda, 1/di[k2] + lambda) )

      for (j in 1:ncol(Y)) {
         keep <- !is.na(Y[, j])
         Mj <- M[keep, , drop = FALSE]
         Mj <- Mj[, apply(Mj, 2, sum) > 0, drop = FALSE]
         Uj <- U[keep, , drop = FALSE]
         Y[keep, j] <- Mj %*% solve(t(Mj) %*% Uj %*% diag(hi) %*% t(Uj) %*% Mj) %*% t(Mj) %*% Uj %*% diag(hi) %*% t(Uj) %*% Y[keep, j]
      }

   } else { stop("Not enough replications!") }

   if (!is.null(eta)) { Y <- Y + adjustments }

   if (average) {
      subM <- apply(M, 2, function(xi) min(which(xi == 1)))
      Y <- Y[subM, ]
      rownames(Y) <- names(subM)
   }

   if (return.info) { return(list(newY = Y, di = di, hi = hi)) } else { return(Y) }

}


