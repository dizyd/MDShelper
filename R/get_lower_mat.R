#' Make a lower traingular part of distance matrix
#'
#' \code{make_lower_mat} makes a lower triangular matrix out of a symmetric distance matrix
#'
#' @param d_mat symmetric matrix
#'
#' @author David Izydorczyk
#'
#' @return lower triangular matrix, with rest `NA`
#' @examples
#'
#' \dontrun{
#'
#' make_lower_mat(d_mat)
#'
#' }
#'
#'
#' @export
make_lower_mat <- function(d_mat){

  lower_mat <- d_mat

  diag(lower_mat) <- NA
  lower_mat[upper.tri(lower_mat)] <- NA


  return(lower_mat)

}


