#' Median disagreement (MDA) plot
#'
#' For each measurement, calculate the shift in its median abundance level after
#' normalisation. The \code{MDA} plot is useful when choosing the tuning
#' parameters.
#'
#' @param raw Raw data before normalisation.
#'
#' @param normalised Normalised data.
#'
#' @param outline Logical. Whether to print the outliers.
#'
#' @param ylim The y-axis limits of the output plot.
#'
#' @param guides Numeric. Horizontal dashed lines for reference.
#'
#' @param title Character. Title of the plot.
#'
#' @return The MDA plot.
#'
#' @export
MDA <- function(raw,
                normalised,
                outline = FALSE,
                ylim = NULL,
                guides = NULL,
                title = NULL) {
   x <- colMedians(normalised, na.rm = TRUE) - colMedians(raw, na.rm = TRUE)
   boxplot(x, pch = 16, cex = 0.4, ylim = ylim,
           outline = outline, main = title)
   if (!is.null(guides))
      abline(h = guides, col = "darkgrey", lty = 2)
}
