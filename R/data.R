#' Example data set derived from TCGA skin cutaneous melanoma (SKCM) data.
#'
#' A dataset containing the processed inputs used in the melanoma analysis within the CRSO publication.
#' @format A list with 3 items
#' \describe{
#'   \item{D}{Binary alteration matrix. Rows are candidate driver events, columns are samples.}
#'   \item{P}{Passenger probability matrix corresponding to D.}
#'   \item{cnv.dictionary}{Data frame containing copy number genes.}
#'   ...
#' }
#' @source Dataset derived from data generated by the TCGA Research Network: https://www.cancer.gov/tcga
"skcm.list"

## @source \url{http://www.diamondse.info/}
##' @name input_list_SKCM
