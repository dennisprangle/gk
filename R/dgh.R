#' @rdname g-and-h
#' @export
dgh = function(x, A, B, g, h, c=0.8, log=FALSE){
    z = pgh(x, A, B, g, h, c, zscale=TRUE)
    if (log) {
        return(stats::dnorm(z, log=TRUE) - Qgh_log_deriv(z, A, B, g, h, c))
    } else {
        return(stats::dnorm(z) / Qgh_deriv(z, A, B, g, h, c))
    }
}

#' g-and-h Q derivative
#'
#' Derivative of the g-and-h Q function.
#'
#' @param z A vector of normal quantiles.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param h Vector of h parameters. Must be positive.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @param getR When \code{TRUE} returns a faster partial calculation that has the same sign as the derivative (Used in checking parameter validity.)
#' @return The derivative of the g-and-h Q function at p.
Qgh_deriv = function(z, A, B, g, h, c=0.8, getR=FALSE) {
    ##Essentially this function calculates
    ##B*exp(h*z^2/2)*((1+c*tanh(g*z/2))*(1+h*z^2) + c*g*z/(2*cosh(g*z/2)^2))
    ##But treats some edge cases carefully

    ##Recycle inputs to same length as output
    n = max(length(z), length(A), length(B), length(g), length(h), length(c))
    zeros = rep(0, n)
    z = z + zeros
    A = A + zeros
    B = B + zeros
    g = g + zeros
    h = h + zeros
    c = c + zeros

    ##Standard calculations
    z_squared = z^2
    if (!getR) {
        term1 = exp(h*z_squared/2)
    }
    term2 = 1+c*tanh(g*z/2)
    term3 = 1+h*z_squared
    term4 = c*g*z/(2*cosh(g*z/2)^2)

    ##Correct edge cases
    ##(likely to be rare so no need for particularly efficient code)
    gzero = (g==0)
    term2[gzero] = 1 ##Avoid possibility of 0*Inf
    term4[gzero] = 0
    term4[is.infinite(z)] = 0
    
    ##Return output
    if (getR) {
        return(term2*term3+term4)
    } else {
        return(B*term1*(term2*term3+term4))
    }
}

#' g-and-h Q log derivative
#'
#' Derivative of the g-and-h log(Q) function.
#'
#' @param z A vector of normal quantiles.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param h Vector of k parameters. Must be positive.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @return The derivative of the g-and-h Q function at p.
Qgh_log_deriv = function(z, A, B, g, h, c=0.8) {
    ##Essentially this function calculates
    ##log(B)+(h*z^2)/2+log((1+c*tanh(g*z/2))*(1+h*z^2) + c*g*z/(2*cosh(g*z/2)^2))
    ##But treats some edge cases carefully

    ##Recycle inputs to same length as output
    n = max(length(z), length(A), length(B), length(g), length(h), length(c))
    zeros = rep(0, n)
    z = z + zeros
    A = A + zeros
    B = B + zeros
    g = g + zeros
    h = h + zeros
    c = c + zeros

    ##Standard calculations
    z_squared = z^2
    term1 = h*z_squared/2
    term2 = 1+c*tanh(g*z/2)
    term3 = 1+h*z_squared
    term4 = c*g*z/(2*cosh(g*z/2)^2)

    ##Correct edge cases
    ##(likely to be rare so no need for particularly efficient code)
    gzero = (g==0)
    term2[gzero] = 1 ##Avoid possibility of 0*Inf
    term4[gzero] = 0
    term4[is.infinite(z)] = 0
    
    ##Return output
    return(log(B)+term1+log(term2*term3+term4))
}
