#' Check parameters valid
check.params <-
function(A,B,g,k,c=0.8,theta){
  if(!is.null(theta)) {
    if(!is.matrix(theta)) {
      if(length(theta)==4) {
        A <- theta[1]
        B <- theta[2]
        g <- theta[3]
        k <- theta[4]
      } else if(length(theta)==5) {
        A <- theta[1]
        B <- theta[2]
        c <- theta[3]
        g <- theta[4]
        k <- theta[5]
      } else {
        stop("gk function called with wrong number of parameters")
      }
    } else {
      if(ncol(theta)==4) {
        A <- theta[,1]
        B <- theta[,2]
        c <- rep(0.8, nrow(theta))
        g <- theta[,3]
        k <- theta[,4]
      } else if(ncol(theta)==5) {
        A <- theta[,1]
        B <- theta[,2]
        c <- theta[,3]
        g <- theta[,4]
        k <- theta[,5]
      } else {
        stop("gk function called with wrong number of parameters")
      }
    }
  } else {
      if (length(c) == 1 & length(A)>1) {
          c <- rep(c, length(A))
      }
  }
  if (length(B) != length(A) | length(c) != length(A) | length(g) != length(A) | length(k) != length(A)) stop("gk function called with parameters vectors of different lengths")
  if (any(B<=0)) stop("gk functions require B>0")
  if (any(k<=-0.5)) stop("gk functions require k>-0.5")
  return(data.frame(A=A,B=B,c=c,g=g,k=k))
}

