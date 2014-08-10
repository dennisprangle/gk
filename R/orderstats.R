#' g&k order statistics
#'
#' Generates a subset of order statistics from independent g-and-k draws
#'
#' @export
#' @param n Total number of independent draws
#' @param orderstats Which order statistics to generate, in increasing order
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters.  Must be positive.
#' @param g Vector of g parameters.
#' @param k Vector of k parameters.  Must be greater than -0.5.
#' @param c Vector of c parameters.  Often fixed at 0.8 (see Rayner and MacGillivray) which is the default.
#' @param theta Vector or matrix of parameter values. If this is supplied all other parameter arguments are ignored. A vector is treated as a single row matrix.  The columns may correspond to either 1) (A,B,g,k) with c taken to equal 0.8 or 2) (A,B,c,g,k).
#' @details Uniform order statistics are generaed by the exponential spacings method (see Ripley for example) and then converted to g&k values.
#' @return A vector of order statistics equal in length to \code{orderstats}
#' @references Brian Ripley `Stochastic Simulation' Wiley (1987)
#' @examples
#' rgk.orderstats(100, c(25,50,75), 3, 1, 2, 0.5)
rgk.orderstats <- function(n, orderstats, A, B, g, k, c=0.8, theta=NULL) {
  ##Check that parameters are scalar.  Other checks are called later.
  if (is.matrix(theta)) {
    if (nrow(theta)>1) stop("Multiple parameter values not permitted")
    theta <- theta[1,]
  }
  if (is.null(theta)) {
    if (length(A) > 1 | length(B) > 1 | length(g) > 1 | length(k) > 1 | length(c) >1) stop("Multiple parameter values not permitted")
  }
  p <- length(orderstats)
  kk <- c(0,orderstats,n+1)
  w <- rgamma(p+1, kk[-1] - kk[1:(p+1)])
  u <- cumsum(w[1:p]) / sum(w)
  qgk(u, A, B, g, k, c, theta)
}

#' Uniform order statistics
#'
#' Generates a subset of order statistics from independent uniform draws
#'
#' @param n Total number of independent draws
#' @param orderstats Which order statistics to generate, in increasing order
#' @return A vector of order statistics
runif.orderstats <- function(n, orderstats) {
  p <- length(orderstats)
  if (orderstats[1]<1 || orderstats[p]>n) stop("orderstats must be between 1 and n")
  kk <- c(0,orderstats,n+1)
  diffs <- kk[-1] - kk[1:(p+1)]
  if (any(diffs<=0)) stop("orderstats must be strictly increasing")
  w <- rgamma(p+1, diffs)
  cumsum(w[1:p]) / sum(w)
}
