#' Proteomics data as raw LFQ intensities on the Sydney Heart Bank (SHB) cohort.
#'
#' LFQ intensities of the SHB proteomics data set with sample information.
#'
#' @docType data
#'
#' @format Matrix of raw LFQ intensities with sample information as a data
#' frame.
#' \describe{
#' \item{raw}{LFQ intensities. }
#' \item{smpinfo}{Sample information. }
#' }
#'
#' @usage
#' data("shb", package = "ruvms")
#'
#' @references
#' Li, M., Parker, B. L., Pearson, E., Hunter, B., Cao, J., Koay, Y. C., ... &
#' O’Sullivan, J. F. (2020). Core functional nodes and sex-specific pathways in
#' human ischaemic and dilated cardiomyopathy. \emph{Nature communications},
#' 11(1), 1-12.
#'
#' @examples
#' data("shb", package = "ruvms")
#' raw[1:5, 1:5]
#' head(smpinfo)
"raw"
"smpinfo"
