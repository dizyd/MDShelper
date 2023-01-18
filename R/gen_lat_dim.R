#' Generate MDS Space
#'
#' @param ndims number of dimensions
#' @param n number of stimuli
#' @param min minimum value the coordinates can have (defualt -1)
#' @param max maximum value the coordinates can have (default 1)
#'
#' @importFrom stats runif
#'
#' @author David Izydorczyk
#'
#' @return matrix with `ndims` columns and `n` rows
#'
#' @examples
#' gen_lat_dim()
#' @export
gen_lat_dim <- function(ndims=2, n=16, min=-1, max= 1){

  # Create matrix with cue values for ndims cues
  latent_dims <- matrix(runif(ndims*n,min,max),nrow=n,ncol=ndims)

  # Return results
  return(latent_dims)

}
