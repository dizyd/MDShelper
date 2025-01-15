#' Leave-One-Out Cross-Validation for MDS
#'
#' \code{loo_CV_MDS} runs a leave-one-out (loo) cross-validation procedure to find the best fitting number of dimensions of the MDS space
#'
#' @param mat     pairwise distance matrix
#' @param max_dim maximum number of dimensions to test (default = 2)
#' @param verbose prints progress (default = FALSE)
#' @param parallel run the main computation in parallel
#' @param n_cores  if parallel == TRUE, how many cores should be used when registering the cluster
#'
#' @importFrom stats  dist
#' @importFrom stats  cor
#' @importFrom smacof mds
#' @importFrom purrr  map2_dbl
#' @importFrom purrr  map2
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom progress  progress_bar
#'
#' @author David Izydorczyk
#'
#' @return Matrix with number of dimensions `dim` and R-squared, RMSE and MSE between the true and predicted distances in each leave-one-out trial
#' @examples
#'
#' \dontrun{
#'
#' loo_CV_MDS(mat,max_dim = 5,verbose = TRUE)
#'
#' }
#'
#'
#' @export
loo_CV_MDS <- function(mat, max_dim = 2, verbose = FALSE, parallel = TRUE, n_cores = 4){

  dims_to_test <- 1:max_dim                  # define vector with dimensions to test
  res_MSE      <- vector("numeric",max_dim)  # create empty results vector
  res_RMSE     <- vector("numeric",max_dim)  # create empty results vector
  res_rsq      <- vector("numeric",max_dim)  # create empty results vector
  dims         <- dim(mat)                   # get dimensions of input matrix
  n_pairs      <- dims[1]*(dims[1]-1)/2      # get number of pairs

  # possible (lower triangular) matrix positions
  rows <- c()
  cols <- c()

  for(col in 1:(dims[1]-1)){
    for(row in (col+1):(dims[2])){

      rows <- c(rows,row)
      cols <- c(cols,col)

    }
  }

  positions <- cbind(rows,cols)

  # Create list of matrices where each entry (i.e., each pairwise distance)
  # is left out once
  # TODO: this code is based on the previous (non-loo) CV version, creating this
  #       list of matrices is probably unnecessary heavy for the RAM

  NA_mats       <- lapply(1:n_pairs,function(x){

    mat[positions[x,][1],positions[x,][2]] <- NA
    return(mat)

  })

  # Extract true pairwise distances
  TRUE_vals <-  mat[lower.tri(mat)] %>% as.vector()


  # if parallel == TRUE register a parallel backend with n_cores
  if(parallel){
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
  }


  # if verbose == TRUE init progressbar
  if(verbose){
    pb <- progress::progress_bar$new(
      format = "  Running [:bar] :percent in :elapsed",
      total = length(dims_to_test), clear = FALSE, width = 60)
  }


  for(dim in dims_to_test){


    # Calculate MDS solution
    MDS_solutions <- lapply(NA_mats,function(x,...){


      mds(x,ndim=dim,"ordinal") %>%
        .$conf %>%
        stats::dist(., method = "euclidean") %>%
        as.matrix()


    },dim,.parallel = parallel)


    # Extract predictions for missing distances based on MDS solution
    PRED_vals     <-  sapply(1:n_pairs,function(x){

      pred <- MDS_solutions[[x]][positions[x,][1],positions[x,][2]]

      return(pred)

    })

    # Calculate Metrics

    ## R squared
    res_rsq[dim]  <- cor(TRUE_vals,PRED_vals)^2

    ## RMSE
    res_RMSE[dim]  <- RMSE(PRED_vals, TRUE_vals)

    ## MSE
    res_MSE[dim]  <- res_RMSE[dim]^2

    # Print progress
    if(verbose){
      pb$tick()
    }

  }

  # stop cluste gain
  if(parallel){
    stopCluster(cl)
  }

  # Combine indices and return results
  return(cbind("dim"  = dims_to_test,
               "R_sq" = res_rsq,
               "RMSE" = res_RMSE,
               "MSE"  = res_MSE))
}
