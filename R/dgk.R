#' @rdname g-and-k
#' @export
dgk <- function(x, A, B, g, k, c=0.8, theta=NULL, log=FALSE, do.message=TRUE){
    if (do.message) cat("dgk estimates the g-and-k density function numerically and may be slow to run\n")
    z <- gk2z(x, A, B, g, k, c, theta)
    params <- check.params(A,B,g,k,c,theta) ##Don't really need to check as done in gk2z.
    invJac <- z2gk.Jacobian(z, params)
    if (log) {
        return(dnorm(z, log = TRUE) - log(invJac))
    } else {
        return(dnorm(z)/invJac)
    }
}

#' g-and-k quantile Jacobian.
#'
#' Derivative of the g&k quantile function given z.
#'
#' @param z A standard normal value.
#' @param theta Vector (A,B,c,g,k).
#' @return The derivative of the g-and-k quantile function at z.
z2gk.Jacobian <- function(z, theta) {
  theta <- matrix(theta, ncol=5)
  colnames(theta) <- NULL
  temp <- (1 + theta[,3]*tanh(theta[,4]*z/2))*(1+(2*theta[,5]+1)*z*z) + (1+z*z)*z*theta[,4]*theta[,3]/(2*(cosh(theta[,4]*z/2)^2))
  theta[,2] * (1+z*z)^(theta[,5]-1) * temp
}
