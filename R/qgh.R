#' g-and-h distribution functions
#'
#' Density, distribution function, quantile function and random generation for the generalised g-and-h distribution
#'
#' @name g-and-h
#' @param n Number of draws to make.
#' @param p Vector of probabilities.
#' @param x Vector of quantiles.
#' @param q Vector of quantiles.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param h Vector of h parameters. Must be non-negative.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @param zscale If true the N(0,1) quantile of the cdf is returned.
#' @param log If true the log density is returned.
#' @details
#'   The (generalised) g-and-h distribution is defined by its quantile function:
#'  \deqn{x(p) = A + B [1 + c \tanh(gz/2)] z \exp(hz^2/2),}{x(p) = A + B [1 + c tanh(gz/2)] z exp(hz^2/2),}
#'  where z is the standard normal quantile of p.
#'  Parameter restrictions include \eqn{B>0} and \eqn{h \geq 0}{h>=0}. Typically c=0.8. For more
#'  background information see the references.
#'
#' \code{rgh} and \code{qgh} use quick direct calculations. However \code{dgh} and \code{pgh} involve slower numerical inversion of the quantile function.
#'
#' Especially extreme values of the inputs will produce \code{pgh} output rounded to 0 or 1 (-Inf or Inf for \code{zscale=TRUE}).
#' The corresponding \code{dgh} output will be 0 or -Inf for \code{log=TRUE}.
#' @return \code{dgh} gives the density, \code{pgh} gives the distribution, \code{qgh} gives the quantile function, and \code{rgh} generates random deviates
#' @references
#' Haynes `Flexible distributions and statistical models in ranking and selection procedures, with applications' PhD Thesis QUT (1998)
#' Rayner and MacGillivray `Numerical maximum likelihood estimation for the g-and-k and generalized g-and-h distributions' Statistics and Computing, 12, 57-75 (2002)
#' @examples
#' p = 1:9/10 ##Some probabilities
#' x = qgh(seq(0.1,0.9,0.1), A=3, B=1, g=2, h=0.5) ##g-and-h quantiles
#' rgh(5, A=3, B=1, g=2, h=0.5) ##g-and-h draws
#' dgh(x, A=3, B=1, g=2, h=0.5) ##Densities of x under g-and-h
#' dgh(x, A=3, B=1, g=2, h=0.5, log=TRUE) ##Log densities of x under g-and-h
#' pgh(x, A=3, B=1, g=2, h=0.5) ##Distribution function of x under g-and-h
NULL

#' @rdname g-and-h
#' @export
qgh = function(p, A, B, g, h, c=0.8){
    z2gh(stats::qnorm(p), A, B, g, h, c)
}
