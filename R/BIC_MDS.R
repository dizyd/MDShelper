#' Compute BIC for MDS according to Lee (2001)
#'
#' \code{BIC_MDS} computes the BIC for a specific number of dimensions of the
#' MDS space according to Lee (2001)
#'
#
#' @param s       Assumed level of data precision, on the same scale as the data (e.g., normalized 0-1 data).
#'  If individual level data is available, \code{\link{compute_BIC_MDS}} estimates s either as the mean
#'  standard error of the pooled cells (default; the precision of the averaged matrix) or as Lee's (2001)
#'  mean standard deviation of the pooled cells. If only aggregate level data is available, this value needs
#'  to be defined. Typical between-person standard deviations for normalized similarity data are
#'  about 0.10 to 0.20 (Lee, 2001; Gronau & Lee, 2020); the corresponding standard error is SD / sqrt(K)
#'  for K persons.
#' @param rss     Residual sum-of-squares of current solution, on the same scale as the data
#' @param n_dims  Number of dimensions number
#' @param n_items Number of items
#'
#'
#' @author David Izydorczyk
#'
#' @return BIC for a specific number of dimensions.
#' @examples
#'
#'
#'
#' BIC_MDS(s = 0.1, rss = 0.4, n_dims = 5, n_items = 16)
#'
#'
#'
#'
#' @export
BIC_MDS <- function(s, rss, n_dims, n_items){

  if (s <= 0) {
    stop("s should be positive")
  }

  if (n_dims < 1 || n_dims != round(n_dims)) {
    stop("n_dims must be a positive integer")
  }

  BIC <- (rss / s^2) + n_dims*n_items*log(n_items*(n_items-1)/2)

  return(BIC)

}
