#' Compute BIC for MDS according to Lee (2001)
#'
#' \code{compute_BIC_MDS} is a generic function that computes the BIC for a specific number of dimensions of the
#'  MDS space according to Lee (2001). The function works on average or individual data, determined by the
#'  \code{\link{class}} of the first argument (i.e., `matrix` for average data or `list` for indidivdual data)
#'
#' @param x must be either a symmetric  *matrix* of averaged normalized distances or a *list* of individual symmetric normalized pair-wise distance matrices
#' @param min_dim minimum number of dimensions to test (default = 1)
#' @param max_dim maximum number of dimensions to test (default = 5)
#' @param s Assumed level of data precision. Has to be provided when d is an averaged matrix (default = 0.1). Otherwise s is computed from the individual level data.
#'
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @importFrom smacof mds
#'
#' @author David Izydorczyk
#'
#' @return data.frame with the number of dimensions `dim`, the stress `Stress`, the BIC for a specific number of dimensions `BIC`, the residual sum-of-squares `RSS` of the used solution, the variance of the distance values accounted for by the solution `VAF`,  the number of parameters `P` according to Lee (2001), and the used precision estiamte `s`
#'
#' @family compute_BIC_MDS
#' @examples
#'
#' \dontrun{
#'
#' compute_BIC_MDS(d,min_dim = 1, max_dim = 5)
#'
#' }
#'
#'
#' @export

compute_BIC_MDS <- function(x, ...) {
  UseMethod("compute_BIC_MDS",  x)
}
