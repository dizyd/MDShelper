#' Compute BIC for MDS according to Lee (2001)
#'
#' Computes the BIC for a specific number of dimensions of the MDS space according to Lee (2001) based on
#' a symmetric  *matrix* of averaged normalized distances
#'
#' @param x must a symmetric *matrix* of averaged normalized distances
#' @param min_dim minimum number of dimensions to test (default = 1)
#' @param max_dim maximum number of dimensions to test (default = 5)
#' @param s Assumed level of data precision, on the same scale as x (default = 0.1). For an averaged matrix
#'  from K persons, use the mean between-person standard deviation divided by sqrt(K).
#'
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @importFrom smacof mds
#'
#' @author David Izydorczyk
#'
#' @return data.frame with the number of dimensions `dim`, the stress `Stress`, the BIC for a specific number of dimensions `BIC`, the residual sum-of-squares `RSS` of the used solution (on the scale of the data), the variance of the distance values accounted for by the solution `VAF`,  the number of parameters `P` according to Lee (2001), and the used precision estiamte `s`
#'
#' @family compute_BIC_MDS
#' @examples
#'
#' \dontrun{
#'
#' compute_BIC_MDS.matrix(d,min_dim = 1, max_dim = 5, s = 0.1)
#'
#' }
#'
#' @method compute_BIC_MDS matrix
#' @export
compute_BIC_MDS.matrix <- function(x, min_dim = 1, max_dim = 5, s = 0.1){


  if (!is.matrix(x) || nrow(x) != ncol(x)) {
    stop("data matrix must be square")
  }

  if (!isSymmetric(unname(x))) {
    stop("data matrix must be symmetric")
  }

  if (s <= 0) {
    stop("s should be positive")
  }

  if (min_dim > max_dim) {
    stop("min_dim should be smaller or equal to max_dim")
  }

  # rename
  d_mat   <- x
  n_items <- nrow(d_mat)
  dims    <- min_dim:max_dim

  # Calculate the variance of the distance matrix according to Lee (1999)
  dbar   <- (sum(d_mat) - sum(diag(d_mat))) / (n_items * (n_items - 1))
  temp   <- (d_mat - dbar)^2
  vard   <- 0.5 * (sum(temp) - sum(diag(temp)))

  # init empty results data.frame
  res <- data.frame("dim"    = dims,
                    "BIC"    = NA,
                    "Stress" = NA,
                    "RSS"    = NA,
                    "VAF"    = NA,
                    "P"      = dims * n_items,
                    "s"      = s)


  for(i in seq_along(dims)){

    # Determine solution
    fit <- mds(d_mat, ndim = dims[i], type = "ordinal")

    # get stress
    res$Stress[i] <- fit$stress

    # get residual sum of squares of MDS solution with i dims, on the scale
    # of the data (fit$rss is on smacof's internally normalized scale)
    rss           <- rss_data_scale(d_mat, fit$conf)$rss
    res$RSS[i]    <- rss

    # variance accounted for
    res$VAF[i]    <- 1 - rss / vard

    # calculate BIC
    res$BIC[i]    <- BIC_MDS(s = s, rss = rss, n_dims = dims[i], n_items = n_items)

  }


  return(res)


}
