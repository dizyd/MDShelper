#' Compute BIC for MDS according to Lee (2001)
#'
#' \code{BIC_MDS} computes the BIC for a specific number of dimensions of the MDS space according to Lee (2001)
#'
#' @param d_list  list of pair-wise distance matrices (upper + diag should be NA)
#' @param max_dim maximum number of dimensions to test (default = 5)
#'
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @importFrom smacof mds
#'
#' @author David Izydorczyk
#'
#' @return Matrix with the number of dimensions `dim`, the stress `stress`, the BIC for a specific number of dimensions `BIC`, the residual sum-of-squares `rss`, and the number of parameters `P` according to Lee (2001), and the estimated precision estiamte `s`
#' @examples
#'
#' \dontrun{
#'
#' compute_BIC_MDS(list_of_distances)
#'
#' }
#'
#'
#' @export
compute_BIC_MDS <- function(d_list,max_dim=5){


  check_NA <- sapply(d_list,function(x){

    d_T <- all(is.na(diag(x)))
    u_T <- all(is.na(x[upper.tri(x)]))

    all(d_T,u_T)

  }) %>% all()

  if(!check_NA){
    stop("Pair-wise distance matrix has non-NAs in upper triangle or diagonal")
  }

  # get number of items and persons
  n_items   <- sqrt(length(d_list[[1]]))

  lapply(d_list, function(x){
    # Normalize
    max_d <- x %>% as.vector() %>% na.omit() %>% max()
    min_d <- 0 # min(dist_mat)

    x %>% as.matrix() %>% apply(.,2,function(y){

      ifelse(is.na(y),NA,(y-min_d)/(max_d-min_d))

      # Check if scale goes from 1 - 7
      # 1 -> 0
      # 7 -> 1
      # 4 ->  0.5

    })
  })-> d_list_norm

  # Calculate s
  d_by_ID <- sapply(X = d_list_norm, as.vector, simplify=T) %>%
    na.omit() %>%
    as.data.frame()

  # Sum over row-wise SDs
  sum_sd_d <- sum(apply(d_by_ID,1,sd))

  # Compute s
  s <- 1/(n_items*(n_items-1)/2)*sum_sd_d
  #print(s) # same as mean(apply(d_by_ID,1,sd))

  # compute aggregate pair-wise distance matrix
  d_aggr <- apply(simplify2array(d_list_norm), 1:2, mean)


  # init empty results data.frame
  res <- data.frame("dim"    = rep(0,max_dim),
                    "stress" = rep(0,max_dim),
                    "BIC"    = rep(0,max_dim),
                    "rss"    = rep(0,max_dim),
                    "P"      = rep(0,max_dim),
                    "s"      = rep(s,max_dim))

  # for i to max_dims
  for(i in 1:max_dim){
    # save stress as well, since it does not cost much computation
    temp <- mds(d_aggr,ndim=i,type="ordinal")

    # get stress
    stress <- temp$stress

    # get residual sum of squares of MDS soluation with i dims
    rss    <-temp$rss

    # calculate BIC
    BIC    <- BIC_MDS(s,rss,i,n_items)

    # safe
    res[i,c("dim","stress","BIC","rss","P")] <- c(i,stress,BIC,rss,i*n_items)

  }



}
