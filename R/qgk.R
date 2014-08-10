#' g&k distribution functions
#'
#' Draw from and calculate the quantile function for the g-and-k distribution
#'
#' @name g-and-k
#' @param n Number of draws to make.
#' @param p Vector of probabilities.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters.  Must be positive.
#' @param g Vector of g parameters.
#' @param k Vector of k parameters.  Must be greater than -0.5.
#' @param c Vector of c parameters.  Often fixed at 0.8 (see Rayner and MacGillivray) which is the default.
#' @param theta Vector or matrix of parameter values. If this is supplied all other parameter arguments are ignored. A vector is treated as a single row matrix.  The columns may correspond to either 1) (A,B,g,k) with c taken to equal 0.8 or 2) (A,B,c,g,k).
#' @details
#'   The g-and-k distribution is defined by its quantile function:
#'  \deqn{x(p) = A + B [1 + c h(z)][1 + z^2]^k z,}{x(p) = A + B [1 + c h(z)] [1 + z^2]^k z,}
#'  where z is the standard normal quantile of p and
#'  \deqn{h(z) = [1-\exp(-gz)]/[1+\exp(-gz)].}{h(z) = [1-exp(-gz)]/[1+exp(-gz)].}
#'  Parameter restrictions are B>0 and k>0.5. Typically c=0.8. For more
#'  background information see the references.
#' @return A vector of quantiles (\code{qgk}) or g-and-k draws (\code{rgk}).
#' @references
#' Haynes `Flexible distributions and statistical models in ranking and selection procedures, with applications' PhD Thesis QUT (1998)
#' Rayner and MacGillivray `Numerical maximum likelihood estimation for the g-and-k and generalized g-and-h distributions' Statistics and Computing, 12, 57-75 (2002)
#' @examples
#' p <- 1:9/10 ##Some probabilities
#' qgk(seq(0.1,0.9,0.1), A=3, B=1, g=2, k=0.5) ##g-and-k quantiles
#' rgk(5, A=3, B=1, g=2, k=0.5) ##g-and-k draws
NULL

#' @rdname g-and-k
#' @export
qgk <-function(p, A, B, g, k, c=0.8, theta=NULL){
  ##nb No need to check parameters here, done in z2gk
    if (any(p<0) || any(p>1)) stop("p values must be between zero and one")
  z <- qnorm(p)
  z2gk(z, A, B, g, k, c, theta)
}
