#' Residual sum of squares on the data scale
#'
#' Rescales the configuration distances by the least-squares factor b
#' (no intercept), so that they are on the scale of the input dissimilarities.
#' Cells that are NA are ignored.
#'
#' @param delta symmetric dissimilarity matrix (data scale)
#' @param conf  configuration (stimuli in rows, dimensions in columns)
#'
#' @importFrom stats dist
#' @keywords internal
#' @noRd
rss_data_scale <- function(delta, conf) {
  delta  <- as.matrix(delta)
  conf_d <- as.matrix(dist(conf))
  keep   <- lower.tri(delta) & !is.na(delta)
  b      <- sum(delta[keep] * conf_d[keep]) / sum(conf_d[keep]^2)
  list(rss = sum((delta[keep] - b * conf_d[keep])^2), b = b, conf_d = conf_d)
}
