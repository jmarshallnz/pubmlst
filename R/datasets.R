#' MLST allelic profile data for Campylobacter jejuni and coli from PubMLST
#'
#' A dataset containing the sequence type number, 7 allele numbers at the housekeeping genes,
#' the clonal-complex number (if available) and whether the strain is like C.coli.
#'
#' @format A data frame with 8160 rows and 10 variables:
#' \describe{
#'   \item{ST}{Sequence type}
#'   \item{ASP}{aspA allele}
#'   \item{GLN}{glnA allele}
#'   \item{GLT}{gltA allele}
#'   \item{GLY}{glyA allele}
#'   \item{PGM}{pgm allele}
#'   \item{TKT}{tkt allele}
#'   \item{UNC}{uncA allele}
#'   \item{CC}{Clonal complex}
#'   \item{Coli}{TRUE if more isolates have been submitted to PubMLST as C.coli, NA if unknown.}
#' }
#' @source \url{http://pubmlst.org}
"pubmlst"