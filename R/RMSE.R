#' Compute Root-Mean-Squared Error
#'
#' Function to compute the root-mean-squared-error (RMSE) between two numeric variables
#'
#' @param x numeric variable
#' @param y numeric variable
#'
#' @author David Izydorczyk
#'
#' @return RMSE
#'
#'
#' @examples
#' RMSE(rnorm(100),rnorm(100))
#'
#' @export
RMSE     <- function(x,y){

  RMSE <- sqrt(mean((x-y)^2))

  return(RMSE)

}
