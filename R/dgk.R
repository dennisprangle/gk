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
  theta <- as.matrix(theta, ncol=5)
  temp <- exp(-theta[,4]*z)
  infcases <- is.infinite(temp)
  tempA <- temp; tempB <- temp
  tempA[!infcases] <- (1-temp[!infcases])/(1+temp[!infcases])
  tempA[infcases] <- -1
  tempB[!infcases] <- temp[!infcases] / (1+temp[!infcases])^2
  tempB[infcases] <- 0
  temp <- (1/theta[,3] + tempA) * (1+(2*theta[,5]+1)*z*z) / (1+z*z) + 2*theta[,4]*z*tempB
  theta[,2]*theta[,3]*temp*(1+z*z)^theta[,5]
}
