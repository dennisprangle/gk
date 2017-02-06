#' g-and-k distribution functions
#'
#' Density, distribution function, quantile function and random generation for the g-and-k distribution.
#'
#' @name g-and-k
#' @param n Number of draws to make.
#' @param p Vector of probabilities.
#' @param x Vector of quantiles.
#' @param q Vector of quantiles.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param k Vector of k parameters. Must be at least -0.5.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @param zscale If true the N(0,1) quantile of the cdf is returned.
#' @param log If true the log density is returned.
#' @details
#'   The g-and-k distribution is defined by its quantile function:
#'  \deqn{x(p) = A + B [1 + c \tanh(gz/2)] z(1 + z^2)^k,}{x(p) = A + B [1 + c tanh(gz/2)] z (1 + z^2)^k,}
#'  where z is the standard normal quantile of p.
#'  Parameter restrictions include \eqn{B>0} and \eqn{k \geq -0.5}{k >= -0.5}. Typically c=0.8. For more
#'  background information see the references.
#'
#' \code{rgk} and \code{qgk} use quick direct calculations. However \code{dgk} and \code{pgk} involve slower numerical inversion of the quantile function.
#'
#' Especially extreme values of the inputs will produce \code{pgk} output rounded to 0 or 1 (-Inf or Inf for \code{zscale=TRUE}).
#' The corresponding \code{dgk} output will be 0 or -Inf for \code{log=TRUE}.
#' @return \code{dgk} gives the density, \code{pgk} gives the distribution, \code{qgk} gives the quantile function, and \code{rgk} generates random deviates
#' @references
#' Haynes `Flexible distributions and statistical models in ranking and selection procedures, with applications' PhD Thesis QUT (1998)
#' Rayner and MacGillivray `Numerical maximum likelihood estimation for the g-and-k and generalized g-and-h distributions' Statistics and Computing, 12, 57-75 (2002)
#' @examples
#' p = 1:9/10 ##Some probabilities
#' x = qgk(seq(0.1,0.9,0.1), A=3, B=1, g=2, k=0.5) ##g-and-k quantiles
#' rgk(5, A=3, B=1, g=2, k=0.5) ##g-and-k draws
#' dgk(x, A=3, B=1, g=2, k=0.5) ##Densities of x under g-and-k
#' dgk(x, A=3, B=1, g=2, k=0.5, log=TRUE) ##Log densities of x under g-and-k
#' pgk(x, A=3, B=1, g=2, k=0.5) ##Distribution function of x under g-and-k
NULL

#' @rdname g-and-k
#' @export
qgk = function(p, A, B, g, k, c=0.8){
    z2gk(stats::qnorm(p), A, B, g, k, c)
}
