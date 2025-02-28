#' Compute BIC for MDS according to Lee (2001)
#'
#' Computes the BIC for a specific number of dimensions of the MDS space according to Lee (2001) based on
#' a symmetric  *matrix* of averaged normalized distances
#'
#' @param d must a symmetric *matrix* of averaged normalized distances
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
#' @examples
#'
#' \dontrun{
#'
#' compute_BIC_MDS(d,min_dim = 1, max_dim = 5, s = 0.1)
#'
#' }
#'
#'
#' @export
compute_BIC_MDS.matrix <- function(x, min_dim = 1, max_dim = 5, s = 0.1){


  if (!is.matrix(x) || nrow(x) != ncol(x)) {
    stop("data matrix must be square")
  }

  if (!all(x == t(x))) {
    stop("data matrix must be symmetric")
  }

  if (s <= 0) {
    stop("sigma should be positive")
  }

  if (min_dim > max_dim) {
    stop("min_dim should be smaller or equal to max_dim")
  }

  # rename
  d_mat <- x

  # init empty results data.frame
  res <- data.frame("dim"    = rep(0,max_dim),
                    "BIC"    = rep(0,max_dim),
                    "Stress" = rep(0,max_dim),
                    "RSS"    = rep(0,max_dim),
                    "VAF"    = rep(0,max_dim),
                    "P"      = rep(0,max_dim),
                    "s"      = rep(s,max_dim))


  n_items <- nrow(d_mat)

  for(dim in min_dim:max_dim){

    # Determine solution
    temp <- mds(d_mat,ndim=dim,type="ordinal")

    # get stress
    stress <- temp$stress

    # get residual sum of squares of MDS soluation with i dims
    rss    <-temp$rss

    # Calculate the variance of the distance matrix according to Lee (1999)
    dbar   <- (sum(d_mat) - sum(diag(d_mat))) / (n_items * (n_items - 1))
    temp   <- (d_mat - dbar)^2
    vard   <- 0.5 * (sum(temp) - sum(diag(temp)))
    VAF    <- 1-rss/vard

    # calculate BIC
    BIC    <- BIC_MDS(s = s, rss = rss, n_dims = dim, n_items = n_items)

    # safe
    res[dim,] <- c(dim,BIC,stress,rss,VAF,dim*n_items,s)

  }


  return(res)


}
