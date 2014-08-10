#' @rdname g-and-k
#' @export
rgk <-function(n, A, B, g, k, c=0.8, theta=NULL){
  ##nb No need to check parameters here, done in z2gk
  z <- rnorm(n)
  z2gk(z, A, B, g, k, c, theta)
}
