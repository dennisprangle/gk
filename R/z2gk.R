#' Transform standard normal draws to g-and-k.
#'
#' @param z A vector of standard normal draws.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters.  Must be positive.
#' @param g Vector of g parameters.
#' @param k Vector of k parameters.  Must be greater than -0.5.
#' @param c Vector of c parameters.  Often fixed at 0.8 (see Rayner and MacGillivray) which is the default.
#' @param theta Vector or matrix of parameter values. If this is supplied all other parameter arguments are ignored. A vector is treated as a single row matrix.  The columns may correspond to either 1) (A,B,g,k) with c taken to equal 0.8 or 2) (A,B,c,g,k).
#' @return A vector of g-and-k values.
z2gk <- function(z, A, B, g, k, c=0.8, theta=NULL){
  params <- check.params(A,B,g,k,c,theta)
  if (length(z) != length(params$A) & length(params$A) > 1) stop("Number of parameters supplied does not equal 1 or number of z values")
  temp <- exp(-params$g*z)
  infcases <- is.infinite(temp)
  temp[!infcases] <- (1-temp[!infcases])/(1+temp[!infcases])
  temp[infcases] <- -1 ##Otherwise we get NaNs
  temp <- params$A + params$B * (1+params$c*temp) * (1+z^2)^params$k * z
  if (params$k < 0) temp[is.infinite(z)] <- z ##Otherwise we get NaNs
  return(temp)
}

#' Transform one standard normal draw to g-and-k.
#'
#' @param z A standard normal draw.
#' @param theta Vector (A,B,c,g,k).
#' @return A g-and-k value.
z2gk.scalar <- function(z, theta) {
  temp <- exp(-theta[4]*z)
  infcases <- is.infinite(temp)
  temp[!infcases] <- (1-temp[!infcases])/(1+temp[!infcases])
  temp[infcases] <- -1 ##Otherwise we get NaNs
  temp <- theta[1] + theta[2] * (1+theta[3]*temp) * (1+z^2)^theta[5] * z
  return(temp)
}
