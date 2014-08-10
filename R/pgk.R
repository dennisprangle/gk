#' @rdname g-and-k
#' @export
pgk <- function(q, A, B, g, k, c=0.8, theta=NULL, do.message=TRUE){
    if (do.message) cat("pgk estimates the g-and-k distribution function numerically and may be slow to run\n")
    ##nb No need to check parameters here, done in gk2z
    z <- gk2z(q, A, B, g, k, c, theta)
    pnorm(z)
}
