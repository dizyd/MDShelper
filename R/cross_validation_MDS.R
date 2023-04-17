#' Cross-Validation for MDS
#'
#' \code{cross_validation_MDS} runs a cross-validation procedure to find the best fitting number of dimensions of the MDS space
#'
#' @param mat     pairwise distance matrix
#' @param reps    number of repetitions (default = 100)
#' @param max_dim maximum number of dimensions to test (default = 2)
#' @param NA_prob probability of turning an entry of `mat` to `NA` (default = 0.2)
#' @param verbose prints progress (default = FALSE)
#' @param metric  CV-metric to use. Either r (correlation) or RMSE (root-mean-squared-error)
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
#'
#' @author David Izydorczyk
#'
#' @return Matrix with number of dimensions `dim` and average RMSE `m_RMSE` over all `reps`.
#' @examples
#'
#' \dontrun{
#'
#' cross_validation_MDS(mat,max_dim = 5)
#'
#' }
#'
#'
#' @export
cross_validation_MDS <- function(mat, reps=100, max_dim=2, NA_prob=0.2, verbose = FALSE,
                                 metric = "r", parallel = TRUE, n_cores = 4){

  dims_to_test <- 1:max_dim                  # define vector with dimensions to test
  res_m        <- vector("numeric",max_dim)  # create empty results vector
  res_sd       <- vector("numeric",max_dim)  # create empty results vector
  res_avg_rsq  <- vector("numeric",max_dim)  # create empty results vector
  dims         <- dim(mat)                   # get dimensions of input matrix
  n_pairs      <- dims[1]*(dims[1]-1)/2      # get number of pairs

  # possible matrix positions
  rows <- c()
  cols <- c()

  for(col in 1:(dims[1]-1)){
    for(row in (col+1):(dims[2])){

      rows <- c(rows,row)
      cols <- c(cols,col)

    }
  }

  positions <- cbind(rows,cols)

  # Draw NA positions

  n_nas         <- round(n_pairs*NA_prob)

  NA_positions  <- replicate(reps,
                             sample(1:n_pairs,size = n_nas),
                             simplify = FALSE)

  # Make NAs
  NA_mats       <- lapply(NA_positions,function(x){

    for(i in x){

      mat[positions[i,][1],positions[i,][2]] <- NA

    }

    return(mat)


  })

  # True entries
  TRUE_vals     <- lapply(NA_positions,function(x){



    temp <- vector("numeric",n_nas)

    for(i in 1:n_nas){

      temp[i] <- mat[positions[x[i],][1],positions[x[i],][2]]

    }

    return(temp)


  })


  # if parallel == TRUE register a parallel backend with n_cores
  if(parallel){
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
  }


  for(dim in dims_to_test){

    MDS_solutions <- lapply(NA_mats,function(x,...){


      mds(x,ndim=dim,"ordinal") %>%
        .$conf %>%
        stats::dist(., method = "euclidean") %>%
        as.matrix()


    },dim,.parallel = parallel)


    PRED_vals     <- map2(NA_positions, MDS_solutions, function(pos,mat){

      temp <- vector("numeric",n_nas)

      for(i in 1:n_nas){

        temp[i] <- mat[positions[pos[i],][1],positions[pos[i],][2]]

      }

      return(temp)

    })

    if(metric == "RMSE"){

      temp_RMSE <- map2_dbl(PRED_vals, TRUE_vals, RMSE)
      temp_rsq  <- map2_dbl(PRED_vals, TRUE_vals, cor)

      res_m[dim]        <- mean(temp_RMSE)
      res_sd[dim]       <- sd(temp_RMSE)
      res_avg_rsq[dim]  <- mean(temp_rsq^2)

    } else if(metric == "r"){

      temp_r   <- map2_dbl(PRED_vals, TRUE_vals, cor)
      temp_rsq <- temp_r^2

      res_m[dim]        <- mean(temp_r)
      res_sd[dim]       <- sd(temp_r)
      res_avg_rsq[dim]  <- mean(temp_rsq)

    }


    if(verbose){
      print(paste0("Dim.: ", dim, "/", max_dim))
    }

  }

  # stop cluste gain
  if(parallel){
    stopCluster(cl)
  }

  return(cbind("dim"       = dims_to_test,
               "M"         = res_m,
               "SD"        = res_sd,
               "Avg. R_sq" = res_avg_rsq))
}
