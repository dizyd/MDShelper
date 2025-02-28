#' Fill distance matrix
#'
#' \code{fill_mat} transforms distance matrix with only lower triangle entries (rest `NA`) to full matrix with 0 on the diagonals
#'
#' @param lower_mat matrix with only lower triangle entries (rest `NA`)
#' @param diag_val  value for diagonal (default = 0)
#'
#' @author David Izydorczyk
#'
#' @return full symmetric distance matrix
#' @examples
#'
#' \dontrun{
#'
#' fill_mat(lower_dist_mat)
#'
#' }
#'
#'
#' @export
fill_mat <- function(lower_mat,diag_val = 0){

  diag(lower_mat) <- diag_val
  lower_mat[upper.tri(lower_mat)] <- t(lower_mat)[upper.tri(lower_mat)]

  return(lower_mat)

}
