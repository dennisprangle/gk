#' @rdname g-and-h
#' @export
dgh = function(x, A, B, g, h, c=0.8, log=FALSE, type=c("generalised", "tukey")){
    type = match.arg(type)
    z = pgh(x, A, B, g, h, c, zscale=TRUE, type=type)
    if (log) {
        base_log_density = stats::dnorm(z, log=TRUE)
        log_deriv = Qgh_log_deriv(z, A, B, g, h, c, type=type)
        out = base_log_density - log_deriv
        out[base_log_density == -Inf] = -Inf # To avoid -Inf + Inf
        return(out)
    } else {
        base_density = stats::dnorm(z)
        deriv = Qgh_deriv(z, A, B, g, h, c, type=type)
        out = base_density / deriv
        out[base_density == 0] = 0 # To avoid 0/0
        return(out)
    }
}

#' g-and-h Q derivative
#'
#' Derivative of the g-and-h Q function.
#'
#' @keywords internal
#' 
#' @param z A vector of normal quantiles.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param h Vector of h parameters. Must be positive.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @param getR When \code{TRUE} returns a faster partial calculation that has the same sign as the derivative (Used in checking parameter validity.)
#' @param type Can be "generalised" (default) or "tukey".
#' @return The derivative of the g-and-h Q function at p.
Qgh_deriv = function(z, A, B, g, h, c=0.8, getR=FALSE, type=c("generalised", "tukey")) {
    ##Essentially this function calculates
    ##B*exp(h*z^2/2)*((1+c*tanh(g*z/2))*(1+h*z^2) + c*g*z/(2*cosh(g*z/2)^2)) (generalised)
    ##B*exp(h*z^2/2)*(exp(g*z) + z*h*(exp(g*z)-1)/g) (Tukey)
    ##But treats some edge cases carefully

    type = match.arg(type)

    ##Recycle inputs to same length as output
    n = max(length(z), length(A), length(B), length(g), length(h), length(c))
    zeros = rep(0, n)
    z = z + zeros
    A = A + zeros
    B = B + zeros
    g = g + zeros
    h = h + zeros
    c = c + zeros

    gzero = (g==0)
    hzero = (h==0)
    zinf = is.infinite(z)

    ##Standard calculations
    z_squared = z^2
    if (!getR) {
        term1 = exp(h*z_squared/2)
        term1[hzero] = 1
    }
    if (type == "generalised") {
        term2 = 1+c*tanh(g*z/2)
        term3 = 1+h*z_squared
        term4 = c*g*z/(2*cosh(g*z/2)^2)

        ##Correct edge cases
        ##(likely to be rare so no need for particularly efficient code)
        term2[gzero] = 1 ##Avoid possibility of 0*Inf
        term4[gzero] = 0
        term4[zinf] = 0

        other_terms = term2*term3+term4
    } else {
        term2 = exp(g*z)
        term3 = expm1(g*z)/g
        term4 = z*h

        ##Correct edge cases
        term2[gzero] = 1  ##Avoid possibility of 0*Inf
        term3[gzero] = z
        term4[hzero] = 0

        other_terms = term2 + term3*term4
        other_terms[hzero] = term2
    }
    
    ##Return output
    if (getR) {
        return(other_terms)
    } else {
        return(B*term1*other_terms)
    }
}

#' g-and-h Q log derivative
#'
#' Derivative of the g-and-h log(Q) function.
#'
#' @keywords internal
#'
#' @param z A vector of normal quantiles.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param h Vector of k parameters. Must be positive.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @param type Can be "generalised" (default) or "tukey".
#' @return The derivative of the g-and-h Q function at p.
Qgh_log_deriv = function(z, A, B, g, h, c=0.8, type) {
    ##Essentially this function calculates
    ##log(B)+(h*z^2)/2+log((1+c*tanh(g*z/2))*(1+h*z^2) + c*g*z/(2*cosh(g*z/2)^2)) (generalised)
    ##log(B)+(h*z^2)/2+log(exp(g*z) + z*h*(exp(g*z)-1)/g) (Tukey)
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

    gzero = (g==0)
    hzero = (h==0)
    zinf = is.infinite(z)

    ##Standard calculations
    z_squared = z^2
    term1 = h*z_squared/2
    term1[hzero] = 0

    if (type=="generalised") {
        term2 = 1+c*tanh(g*z/2)
        term3 = 1+h*z_squared
        term4 = c*g*z/(2*cosh(g*z/2)^2)

        ##Correct edge cases
        ##(likely to be rare so no need for particularly efficient code)
        term2[gzero] = 1 ##Avoid possibility of 0*Inf
        term4[gzero] = 0
        term4[zinf] = 0

        other_terms = log(term2*term3+term4)
    } else {
        gz = g*z
        gz[gzero] = 0 # In case z infinite
        est1 = gz + log1p((1-exp(-gz))*z*h/g)
        est2 = log(exp(gz) + expm1(gz)*z*h/g)

        other_terms = ifelse(gz>1, est1, est2)

        ##Correct edge cases
        other_terms[gzero] = log1p(h*z_squared)
        other_terms[hzero] = gz
        other_terms[gzero & hzero] = 0
    }

    ##Return output
    return(log(B)+term1+other_terms)
}
