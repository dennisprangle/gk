#' @rdname g-and-k
#' @export
dgk = function(x, A, B, g, k, c=0.8, log=FALSE){
    z = pgk(x, A, B, g, k, c, zscale=TRUE)
    if (log) {
        return(stats::dnorm(z, log=TRUE) - Qgk_log_deriv(z, A, B, g, k, c))
    } else {
        return(stats::dnorm(z) / Qgk_deriv(z, A, B, g, k, c))
    }
}

#' g-and-k Q derivative
#'
#' Derivative of the g-and-k Q function.
#'
#' @keywords internal
#' 
#' @param z A vector of normal quantiles.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param k Vector of k parameters. Must be at least -0.5.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @param getR When \code{TRUE} returns a faster partial calculation that has the same sign as the derivative (Used in checking parameter validity.)
#' @return The derivative of the g-and-k Q function at p.
Qgk_deriv = function(z, A, B, g, k, c=0.8, getR=FALSE) {
    ##Essentially this function calculates
    ##B*(1+z^2)^k*((1+c*tanh(g*z/2))*((1+(2*k+1)*z^2)/(1+z^2)) + g*z/(2*cosh(g*z/2)^2))
    ##But treats some edge cases carefully

    ##Recycle inputs to same length as output
    n = max(length(z), length(A), length(B), length(g), length(k), length(c))
    zeros = rep(0, n)
    z = z + zeros
    A = A + zeros
    B = B + zeros
    g = g + zeros
    k = k + zeros
    c = c + zeros

    ##Standard calculations
    z_squared = z^2
    if (!getR) {
        term1 = (1+z_squared)^k
    }
    term2 = 1+c*tanh(g*z/2)
    term3 = (1+(2*k+1)*z^2)/(1+z^2)
    term4 = c*g*z/(2*cosh(g*z/2)^2)

    ##Correct edge cases
    ##(likely to be rare so no need for particularly efficient code)
    gzero = (g==0)
    term2[gzero] = 1 ##Avoid possibility of 0*Inf
    term4[gzero] = 0
    zbig = is.infinite(z_squared)
    if (!getR) {
        term1[zbig] = abs(z[zbig])^(2*k)
    }
    term3[zbig] = 2*k+1
    term4[is.infinite(z)] = 0

    ##Return output
    if (getR) {
        return(term2*term3+term4)
    } else {
        return(B*term1*(term2*term3+term4))
    }
}

#' g-and-k Q log derivative
#'
#' Derivative of the g-and-k log(Q) function.
#'
#' @keywords internal
#' 
#' @param z A vector of normal quantiles.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param k Vector of k parameters. Must be greater than -0.5.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @return The derivative of the g-and-k Q function at p.
Qgk_log_deriv = function(z, A, B, g, k, c=0.8) {
    ##Essentially this function calculates
    ##log(B)+k*log(1+z^2)+log((1+c*tanh(g*z/2))*((1+(2*k+1)*z^2)/(1+z^2)) + g*z/(2*cosh(g*z/2)^2))
    ##But treats some edge cases carefully

    ##Recycle inputs to same length as output
    n = max(length(z), length(A), length(B), length(g), length(k), length(c))
    zeros = rep(0, n)
    z = z + zeros
    A = A + zeros
    B = B + zeros
    g = g + zeros
    k = k + zeros
    c = c + zeros

    ##Standard calculations
    z_squared = z^2
    term1 = k*log(1+z_squared)
    term2 = 1+c*tanh(g*z/2)
    term3 = (1+(2*k+1)*z^2)/(1+z^2)
    term4 = c*g*z/(2*cosh(g*z/2)^2)

    ##Correct edge cases
    ##(likely to be rare so no need for particularly efficient code)
    gzero = (g==0)
    term2[gzero] = 1 ##Avoid possibility of 0*Inf
    term4[gzero] = 0
    zbig = is.infinite(z_squared)
    term1[zbig] = 2*k*log(abs(z[zbig]))
    term1[k==0] = 0 ##Avoid possibility of 0*Inf
    term3[zbig] = 2*k+1
    term4[is.infinite(z)] = 0

    ##Return output
    return(log(B)+term1+log(term2*term3+term4))
}
