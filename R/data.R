#' Log2 transformed LFQ intensity measurements of the embryonic stem cell (ESC)
#' proteomes post epiblast-like cell induction
#'
#' Processed log2 transformed LFQ intensities of the ESC data set with sample
#' information. 
#'
#' @docType data
#'
#' @format Matrix of log2 intensities and sample inforamation in a data frame.
#' \describe{
#' \item{esc}{Log2-transformed LFQ intensities. }
#' \item{smpinfo}{Sample information. }
#' }
#'
#' @usage
#' data(c("esc", "smpinfo"), package = "brainClass")
#'
#' @references
#' Yang, P., Humphrey, S. J., Cinghu, S., Pathania, R., Oldfield, A. J., Kumar,
#' D., ... & Jothi, R. (2019). Multi-omic profiling reveals dynamics of the
#' phased progression of pluripotency. \emph{Cell systems}, 8(5), 427-445.
#'
#' @examples
#' data(c("esc", "smpinfo"), package = "brainClass")
#' esc[1:5, 1:5]
#' smpinfo
"esc"
"smpinfo"