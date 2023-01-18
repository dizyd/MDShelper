#' Generate MDS data (i.e., true underlying dimensions and true pairwise distances)
#'
#' \code{gen_data_mds} generates a matrix with thetrue underlying dimensions and the corresponding true pairwise distances
#'
#' @param ndims number of latent dimensions/cues (default = 2)
#' @param n     number of itmes/rows (default = 16)
#' @param min   minimum possible cue value for the uniform distribution (default = -1)
#' @param max   maximum possible cue value for the uniform distribution (default = 1)
#' @param mat   matrix with coordinates of stimuli (i.e., the position in multidimensional space), with dimensions in the columns and stimuli in the rows
#' @param r     r parameter of minkowski space (default = 2)
#' @param diag  return diagonal matrix values (defualt = FALSE)
#' @param upper return upper matrix half (default = FALSE)
#' @param return_latent should the latent matrix of underlying cues and cue values be returned? (default = TRUE)
#'
#' @author David Izydorczyk
#'
#' @return matrix of pairwise-distances between objects (if return_latent = FALSE) or a list with the matrix of pairwise distances (`dist_mat`) and the underlying MDS space (`latent`)
#' @examples
#'
#' gen_data_mds()
#'
#' @importFrom Rcpp sourceCpp
#'
#' @export
gen_data_mds <- function(ndims=2, n=16, min=-1, max= 1, r=2, diag=F, upper=F, return_latent=T){

  # Create matrix with cue values for ndims cues
  latent_dims <- gen_lat_dim(ndims, n, min, max)

  # Make pair-wise distance matrix out of this
  dist_mat    <- compute_dists(latent_dims,r,diag,upper)


  if(return_latent){
    res_list <- list("latent" = latent_dims,"dist_mat" = dist_mat)
    return(res_list)
  } else {
    return(dist_mat)
  }

}


