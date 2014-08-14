#' g-and-k maximum likelihood estimation
#'
#' Estimate the MLEs of A,B,g and k given independent identically distributed g-and-k values
#'
#' @export
#' @param x g-and-k values
#' @param theta0 Initial guess of (A,B,g,k)
#' @param details If TRUE the full \code{optim} output is returned
#' @param do.message If TRUE a message is shown that the numerical calculations may be time consuming
#' @param ... Further options for \code{optim}
#' @details
#' The maximum likelihood estimates of (A,B,g,k) are calculated taking c=0.8. This is done by numerical optimisation of \code{dgk} which itself must be estimated numerically. Therefore this function may be time-consuming and the results guaranteed to be accurate.
#' @return If \code{details==FALSE}, a vector of (A,B,g,k) MLEs. Otherwise, full output from \code{optim}.
#' @examples
#' x <- rgk(5, A=3, B=1, g=2, k=0.5)
#' gk.mle(x, c(5,5,5,5))
gk.mle <- function(x, theta0, details=FALSE, do.message=TRUE, ...) {
    if (do.message) cat("gk.mle must repeatedly estimate the g-and-k density function numerically and may be very slow to run\n")
    if (!is.vector(theta0) || length(theta0)!=4) stop("gk function called with wrong number of parameters")
    check.params(theta=theta0)
    if (!is.numeric(x) || !is.vector(x)) stop("gk.mle requires x to be a numeric vector")
    tomin <- function(theta){
        if(theta[2]<=0 || theta[4] <= 0.5) return(Inf)
        -sum(dgk(x, theta=c(theta[1:2],0.8,theta[3:4]), log=TRUE, do.message=FALSE))
    }
    y <- optim(theta0, tomin, method="Nelder-Mead", ...)
    if (details) return(y)
    if (y$convergence != 0) warning(paste("gk.mle optimisation failed, error", y$convergence, y$message))
    return(y$par)
}
