#' Hold-Out Cross-Validation for MDS (Richie & Verheyen, 2020)
#'
#' \code{holdout_CV_MDS} runs the cross-validation procedure of Richie and Verheyen (2020)
#' to find the best fitting number of dimensions of the MDS space: in each of `reps`
#' random splits, a proportion `test_prop` of the pairs is held out, the MDS is fitted to
#' the remaining pairs, and the held-out dissimilarities are correlated with the
#' distances of the fitted configuration. The dimensionality with the highest median
#' correlation is chosen.
#'
#' @param mat       pairwise dissimilarity matrix
#' @param weightmat optional matrix of weights for each pair (e.g., inter-subject precision,
#'                  1 / variance of the ratings across subjects; default = NULL, i.e., equal weights)
#' @param max_dim   maximum number of dimensions to test (default = 2)
#' @param reps      number of random splits (default = 10)
#' @param test_prop proportion of pairs held out in each split (default = 0.2)
#' @param max_row_missing maximum proportion of a stimulus' pairs that may be held out;
#'                  splits that hold out more are redrawn (default = 0.75)
#' @param type      MDS type passed to \code{smacof::mds} (default = "ratio")
#' @param eps     convergence criterion passed to \code{smacof::mds} (default = 1e-8; smacof's default
#'                of 1e-6 stops too early in high dimensions when a cell is missing, see NEWS)
#' @param itmax   maximum number of iterations passed to \code{smacof::mds} (default = 1e5)
#' @param verbose   prints progress (default = FALSE)
#'
#' @importFrom stats  dist
#' @importFrom stats  cor
#' @importFrom stats  median
#' @importFrom smacof mds
#'
#' @author David Izydorczyk
#'
#' @references Richie, R., & Verheyen, S. (2020). Using cross-validation to determine
#'   dimensionality in multidimensional scaling. Proceedings of the 18th International
#'   Conference on Cognitive Modeling.
#'
#' @return Matrix with number of dimensions `dim`, the median (`Md_r`) and mean (`M_r`)
#'   correlation between the held-out dissimilarities and the predicted distances over all `reps`.
#' @examples
#'
#' \dontrun{
#'
#' holdout_CV_MDS(mat, max_dim = 5)
#'
#' }
#'
#'
#' @export
holdout_CV_MDS <- function(mat, weightmat = NULL, max_dim = 2, reps = 10, test_prop = 0.2,
                           max_row_missing = 0.75, type = "ratio", verbose = FALSE,
                           eps = 1e-8, itmax = 1e5){

  dims_to_test <- 1:max_dim                  # define vector with dimensions to test
  n            <- nrow(mat)                  # number of stimuli
  n_pairs      <- n*(n-1)/2                  # get number of pairs
  n_test       <- round(n_pairs*test_prop)   # number of held-out pairs per split

  if(is.null(weightmat)){
    weightmat <- matrix(1, n, n)
  }

  # possible (lower triangular) matrix positions
  positions <- which(lower.tri(mat), arr.ind = TRUE)

  # Draw the held-out pairs of each split; redraw a split if a stimulus
  # loses more than max_row_missing of its pairs
  draw_split <- function(){

    repeat{

      test       <- sample(1:n_pairs, size = n_test)
      n_missing  <- tabulate(c(positions[test, 1], positions[test, 2]), nbins = n)

      if(all(n_missing <= max_row_missing*(n-1))){
        return(test)
      }

    }

  }

  splits <- replicate(reps, draw_split(), simplify = FALSE)

  res_md <- matrix(NA, nrow = reps, ncol = max_dim)  # correlation per split and dimension

  for(r in 1:reps){

    test_pos <- positions[splits[[r]], , drop = FALSE]

    # Training matrix: held-out pairs get NA and weight 0 (smacof only sets
    # the weights of missing values to 0 if no weightmat is given)
    train   <- mat
    train_w <- weightmat
    train[test_pos]   <- NA
    train_w[test_pos] <- 0

    TRUE_vals <- mat[test_pos]

    for(dim in dims_to_test){

      fit       <- mds(train, ndim = dim, type = type, weightmat = train_w, eps = eps, itmax = itmax)
      PRED_vals <- as.matrix(dist(fit$conf))[test_pos]

      res_md[r, dim] <- cor(TRUE_vals, PRED_vals)

    }

    if(verbose){
      print(paste0("Split: ", r, "/", reps))
    }

  }

  # Combine indices and return results
  return(cbind("dim"  = dims_to_test,
               "Md_r" = apply(res_md, 2, median),
               "M_r"  = colMeans(res_md)))
}
