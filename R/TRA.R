#' Technical replicate agreement (TRA) plot
#'
#' Technical replicate agreement (TRA) plots provide a visualisation diagnostic
#' that measures the closeness among technical replicates.
#'
#' @param data The m-by-n data matrix with samples on the rows and measurements (
#' proteins etc.) on the columns. Missing data (NAs) are allowed.
#'
#' @param replicates A vector indicating which samples are replicates of each
#' other. Obsolete when \code{M} is not \code{NULL}.
#'
#' @param M The replication structure matrix.
#'
#' @param col Color of the TRA boxes. Default is \code{NULL}.
#'
#' @param ylim The y-axis limits of the output plot.
#'
#' @param printLabel Logical. Whether to print the names of the technical
#' replicates.
#'
#' @param printLegend Logical. Whether to print the color legend.
#'
#' @param guides Numeric. Horizontal dashed lines for reference.
#'
#' @param title Character. Title of the output plot.
#'
#' @details
#' The assumption of TRA plots is that if all unwanted variations are
#' effectively removed, replications of a unique effective sample are
#' identical after normalisation, even if we do not declare them as replicates
#' of each other
#'
#' @return The TRA plot.
#'
#' @references
#' Molania, R., Gagnon-Bartsch, J. A., Dobrovic, A., & Speed, T. P. (2019).
#' A new normalization for Nanostring nCounter gene expression data.
#' \emph{Nucleic acids research}, 47(12), 6073-6083.
#'
#' @export
TRA <- function(data,
                replicates = NULL,
                M = NULL,
                col = NULL,
                ylim = NULL,
                printLabel = TRUE,
                printLegend = FALSE,
                guides = NULL,
                title = NULL) {
   if (is.null(replicates))
      if (is.null(M))
         stop("Need either M or replicates to indicate relications!") else {
            replicates <- matrix(colnames(M),
                                 nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
            replicates <- replicates[M == 1]
         }
   vars.in.replicates <- do.call(cbind, lapply(unique(replicates), function(ri)
      colVars(data[replicates == ri, ])))
   colnames(vars.in.replicates) <- paste("TR", unique(replicates), sep = "~")
   boxplot(vars.in.replicates,
           outline = FALSE, boxfill = col, xaxt = "none", ylim = ylim)
   if (printLabel)
      text(1:length(colnames(vars.in.replicates)),
           label = colnames(vars.in.replicates),
           par("usr")[3], srt = 60, xpd = TRUE, adj = 1.4, cex = 0.6)
   if (printLegend)
      legend("topright", legend = colnames(vars.in.replicates),
             fill = col, bty = "n", horiz = TRUE, cex = 0.4, x.intersp = 0.6)
   if (!is.null(guides))
      abline(h = guides, col = "darkgrey", lty = 2)
   if (!is.null(title))
      title(title)
}
