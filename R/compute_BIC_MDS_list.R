#' Compute BIC for MDS according to Lee (2001)
#'
#' Computes the BIC for a specific number of dimensions of the MDS space according to Lee (2001) based on
#' a list of indidivdual normalized symmetric distance matrices
#'
#' @param x must a *list* of individual symmetric normalized pair-wise distance matrices
#' @param min_dim minimum number of dimensions to test (default = 1)
#' @param max_dim maximum number of dimensions to test (default = 5)
#' @param s_type how the data precision s is computed from the individual data: `"se"` (default) uses the
#'  mean standard error of the pooled cells, i.e. the precision of the averaged matrix; `"sd"` uses the
#'  mean standard deviation of the pooled cells, as in Lee (2001, Eq. 6)
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
#' compute_BIC_MDS.list(d,min_dim = 1, max_dim = 5)
#'
#' }
#'
#' @method compute_BIC_MDS list
#' @export
#'
compute_BIC_MDS.list <- function(x, min_dim = 1, max_dim = 5, s_type = c("se", "sd")){

  s_type <- match.arg(s_type)

  if (sapply(x,function(x){!is.matrix(x) || nrow(x) != ncol(x)}) |> any()) {
    stop("ALL data matrices must be square")
  }

  if (sapply(x,function(x){!isSymmetric(unname(x))}) |> any()) {
    stop("ALL data matrices must be symmetric")
  }

  if (min_dim > max_dim) {
    stop("min_dim should be smaller or equal to max_dim")
  }


  # One row per pair (lower triangular entries), one column per person
  temp <- sapply(x,function(x){
    x[lower.tri(x)]
  })

  # Between-person standard deviation of each pooled cell
  SDs <- apply(temp,1,function(x){
    sd(x,na.rm=TRUE)
  })

  # Determine s, either as the average of the standard errors ("se") or as
  # the average of the standard deviations ("sd", Lee, 2001, Eq. 6) of the
  # pooled cells
  if (s_type == "se") {
    s <- mean(SDs / sqrt(rowSums(!is.na(temp))))
  } else {
    s <- mean(SDs)
  }


  # Compute aggregate pair-wise distance matrix
  d_mat <- apply(simplify2array(x), 1:2, mean, na.rm = TRUE)


  # Compute BIC on the aggregate matrix with the estimated s
  res <- compute_BIC_MDS.matrix(d_mat, min_dim = min_dim, max_dim = max_dim, s = s)


  return(res)


}
