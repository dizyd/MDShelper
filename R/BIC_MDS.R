#' Compute BIC for MDS according to Lee (2001)
#'
#' \code{BIC_MDS} computes the BIC for a specific number of dimensions of the
#' MDS space according to Lee (2001)
#'
#
#' @param s       Assumed level of data precision. If individual level data is available,
#'  this can is calculated as the average of the standard errors for each of the
#'  pooled cells in the final matrix. If only aggregate level data is available, this value needs to be defined.
#'  Sensible values from previous studies for normalized similarity matrices showed s values between approximately 0.05 and 0.15,
#'  which could be interpreted as ‘precise’ and ‘imprecise’ data sets.
#' @param rss     Residual sum-of-squares of current solution
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

