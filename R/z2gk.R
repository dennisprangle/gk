#' Transform standard normal draws to g-and-k draws.
#'
#' @keywords internal
#' 
#' @param z Vector of N(0,1) draws.
#' @param A Vector of location parameters.
#' @param B Vector of scale parameters. Must be positive.
#' @param g Vector of g parameters.
#' @param k Vector of k parameters. Must be at least -0.5.
#' @param c Vector of c parameters. Often fixed at 0.8 which is the default.
#' @return A vector of g-and-k values.
z2gk = function(z, A, B, g, k, c=0.8){
    ##Essentially this function calculates
    ##x = A + B * (1+c*tanh(g*z/2)) * z*(1+z^2)^k
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

    ##Standard calculatations
    z_squared = z^2
    term1 = (1+c*tanh(g*z/2))
    term2 = z*(1+z_squared)^k

    ##Correct edge cases
    ##(likely to be rare so no need for particularly efficient code)
    term1[g==0] = 1 ##Avoids possibility of 0*Inf
    zbig = which(is.infinite(z_squared))
    term2[zbig] = sign(z[zbig]) * abs(z[zbig])^(1 + 2*k[zbig]) ##Correct when abs(z) large or infinite

    ##Return output
    return(A + B*term1*term2)
}
