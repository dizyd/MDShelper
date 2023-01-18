#' Compute pairwise distances between two vectors
#'
#' \code{compute_dists} is a wrapper around `pair_dist_cpp` to compute pairwise distances between two vectors.
#'
#' @param mat   matrix with coordinates of stimuli (i.e., the position in multidimensional space), with dimensions in the columns and stimuli in the rows
#' @param r     r parameter of minkowski space (default = 2)
#' @param diag  return diagonal matrix values (defualt = FALSE)
#' @param upper return upper matrix half (default = FALSE)
#'
#'
#' @author David Izydorczyk
#'
#' @return Pairwise distances. By defaul the diagonal and upper triangle of the matrix are `NA`
#' @examples
#' \dontrun{
#'
#' compute_dists(matrix(rnorm(10),nrow=5),r = 2)
#'
#' }
#'
#' @importFrom Rcpp sourceCpp
#'
#' @export
compute_dists <- function(mat,r=2,diag=F,upper=F){

  result <- pair_dist_cpp(mat,r)

  if(!diag){
    diag(result) <- NA
  }

  if(!upper){
    result[upper.tri(result)] <- NA
  }


  return(result)

}


