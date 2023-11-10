#' Relative log abundance (RLA) plot
#'
#' The \code{RLA} function is a wrapper function to help easily visualise
#' normalisation results on a data set. It is equivalent to the reltaive log
#' expression (RLE) plot under the gene expression contexts.
#'
#' @param data The m-by-n data matrix with samples on the rows and measurements (
#' proteins etc.) on the columns. Missing data (NAs) are allowed.
#'
#' @param smpID Sample IDs of the data. Default is \code{NULL}. This argument is
#' useful when one wishes to display the samples in an order different from the
#' sample orders of input data.
#'
#' @param repCol The replication color codes. Useful when one wishes to colour
#' the RLA boxes to indicate replications. Default is \code{NULL}.
#'
#' @param repLabel Color legend names.
#'
#' @param smpName Sample names. Sample names to display on the axis.
#' Default is \code{NULL}.
#'
#' @param ylim The y-axis labels. A numeric vector of length 2.
#'
#' @param guides Numeric. Horizontal dashed lines for reference.
#'
#' @param title Character. Title of the output plot.
#'
#' @param ... Extra parameters for text() for sample lable appearances. Eg,
#' adj = 1.
#'
#' @return The RLA plot.
#'
#' @references
#' Gandolfo, L. C., & Speed, T. P. (2018). RLE plots: Visualizing unwanted
#' variation in high dimensional data. \emph{PloS one}, 13(2), e0191629.
#'
#' @export
RLA <- function(data,
                smpID = NULL,
                repCol = NULL,
                repLabel = NULL,
                smpName = NULL,
                ylim = NULL,
                guides = NULL,
                title = NULL, ...) {
   y_med <- t(sweep(data, 2, Rfast::colMedians(data, na.rm = TRUE), "-"))
   if (is.null(smpID)) smpID <- rownames(data) else y_med <- y_med[, smpID]
   boxplot(y_med, outline = FALSE, boxfill = repCol, xaxt = "none", ylim = ylim)
   abline(h = 0, col = "darkgrey", lty = 2)
   if (!is.null(smpName))
      text(1:length(smpID), label = smpName,
           par("usr")[3], srt = 60, xpd = TRUE, cex = 0.6, ...)
   if (!is.null(guides))
      abline(h = guides, col = "darkgrey", lty = 2)
   if (!is.null(repLabel))
      legend("topright", legend = repLabel,
             fill = unique(repCol), bty = "n", horiz = TRUE,
             cex = 0.6, x.intersp = 0.6)
   if (!is.null(title))
      title(title)
}
