#' Compute BIC for MDS according to Lee (2001)
#'
#' \code{BIC_MDS} computes the BIC for a specific number of dimensions of the MDS space according to Lee (2001)
#'
#' @param s       sample estimate of standard deviation
#' @param rss     residual sum-of-squares of current solution
#' @param n_dims  number of dimensions number
#' @param n_items number of items
#'
#'
#' @author David Izydorczyk
#'
#' @return BIC for a specific number of dimensions.
#' @examples
#'
#' \dontrun{
#'
#' BIC_MDS()
#'
#' }
#'
#'
#' @export
BIC_MDS <- function(s,rss,n_dims,n_items){

  BIC = 1/s^2 * rss + n_dims*n_items*log(n_items*(n_items-1)/2)

  return(BIC)
}
