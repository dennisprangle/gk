#' Transform standard normal draws to g-and-h draws.
#'
#' @keywords internal
#' 
#' @param z Vector of N(0,1) draws.
#' @param A Vector of location parameters.
#' @param B Vector of scale parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param h Vector of h parameters.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @param type Can be "generalised" (default) or "tukey".
#' @return A vector of g-and-h values.
z2gh = function(z, A, B, g, h, c=0.8, type=c("generalised", "tukey")) {
    ##Essentially this function calculates
    ##x = A + B * (1+c*tanh(g*z/2)) * z*exp(0.5*h*z^2) (generalised) or
    ##x = A + B * ((exp(g*z) - 1)/g) * exp(0.5*h*z^2) (Tukey)
    ##but treats some edge cases carefully

    ##Recycle inputs to same length as output
    n = max(length(z), length(A), length(B), length(g), length(h), length(c))
    zeros = rep(0, n)
    z = z + zeros
    A = A + zeros
    B = B + zeros
    g = g + zeros
    h = h + zeros
    c = c + zeros

    if (type == "generalised") {
        ##Standard calculatations
        term1 = 1+c*tanh(g*z/2)
        term2 = z*exp(0.5*h*z^2)

        ##Correct edge cases
        ##(likely to be rare so no need for particularly efficient code)
        term1[g==0] = 1 ##Avoids possibility of 0*Inf
        term2[h==0] = z

        output = A + B*term1*term2
    } else if (type == "tukey") {
        ##Standard calculatations
        term1 = expm1(g*z)/g
        term2 = exp(0.5*h*z^2)

        ##Correct edge cases
        term1[g==0] = z

        output = A + B*term1*term2
    } else {
        stop("Only generalised and tukey type are implemented")
    }

    return(output)
}
