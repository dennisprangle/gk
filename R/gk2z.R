#' Transform g-and-k draws to standard normal by numerical inversion.
#'
#' @param x A vector of g-and-k draws.
#' @param A Vector of A (location) parameters.
#' @param B Vector of B (scale) parameters.  Must be positive.
#' @param g Vector of g parameters.
#' @param k Vector of k parameters.  Must be greater than -0.5.
#' @param c Vector of c parameters.  Often fixed at 0.8 (see Rayner and MacGillivray) which is the default.
#' @param theta Vector or matrix of parameter values. If this is supplied all other parameter arguments are ignored. A vector is treated as a single row matrix.  The columns may correspond to either 1) (A,B,g,k) with c taken to equal 0.8 or 2) (A,B,c,g,k).
#' @return A vector of standard normal values.
gk2z <-function(x, A, B, g, k, c=0.8, theta=NULL) {
  params <- check.params(A,B,g,k,c,theta)
  if (any(params$c<0)) stop("Case c<0 not supported")
  if (any(params$c>=1)) stop("Case c>=1 not supported")
  if (length(x) != length(params$A) & length(params$A) > 1) stop("Number of parameters supplied does not equal 1 or number of x values")
  if (nrow(params)==1) {
    gk2z_i <- function(i){
      gk2z.scalar(x[i], as.numeric(params[1,]))
    }
  } else {
    gk2z_i <- function(i){
      gk2z.scalar(x[i], as.numeric(params[i,]))
    }
  }
  sapply(seq(along.with=x), gk2z_i)
}

#' Transform one g-and-k draw to standard normal.
#'
#' @param x A g-and-k value
#' @param theta Vector (A,B,c,g,k).
#' @return A standard normal draw.
gk2z.scalar <- function(x, theta) {
  ##Convert to standard scale and location
  xx <- (x - theta[1]) / theta[2]
  if (xx == Inf) return(Inf)  
  xx <- as.numeric(xx)
  theta2 <- theta
  theta2[1:2] <- c(0,1)
  ##Find an interval to search for z
  if (xx == 0) return(0)
  if (theta[5] >= 0) {
    if (xx > 0) {
      inter <- c(0,xx)
    } else {
      inter <- c(xx/(1-theta[3]),0)
    }
  } else {
    if (xx > 0) {
      xx1 <- exp(-theta[4])
      if (xx1 < Inf) {
          xx1 <-  (1-xx1)/(1+xx1)
      } else {
          xx1 <- -1
      }
      xx1 <- (1 + theta[3]*xx1)*2^theta[5]
      if (xx <= xx1) {
        inter <- c(0,2*xx)
      } else {
        inter <- c(0, (2*xx)^(1/(1+2*theta[5])))
      }
    } else {
      xx1 <- exp(theta[4])
      if (xx1 < Inf) {
          xx1 <-  (1-xx1)/(1+xx1)
      } else {
          xx1 <- -1
      }
      xx1 <- -(1 + theta[3]*xx1)*2^theta[5]
      if (xx >= xx1) {
        inter <- c(2*xx / (1-theta[3]), 0)
      } else {
        inter <- c(-(-2*xx/(1-theta[3]))^(1/(1+2*theta[5])), 0)
      }
    }
  }        
  ##Numerically solve xx=qgk(z,theta2) for z
  uniroot(function(z){z2gk.scalar(z,theta2)-xx}, inter)$root
}
