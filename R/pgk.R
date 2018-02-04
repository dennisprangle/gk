#' Distribution function for the g-and-k distribution
#'
#' @keywords internal
#' 
#' @param q Quantiles.
#' @param A A (location) parameter.
#' @param B B (scale) parameter. Must be positive.
#' @param g g parameter.
#' @param k k parameter. Must be at least -0.5.
#' @param c c parameter. Often fixed at 0.8 which is the default.
#' @param zscale When true returns the N(0,1) quantile of the cdf (needed by dgk).
#' @details
#'  This internal function performs the calculation assuming scalar inputs.
#'  The exported function is a vectorised wrapper of this.
#' @return The cumulative probability.
pgk_scalar = function(q, A, B, g, k, c=0.8, zscale=FALSE){
    toroot = function(p) {z2gk(p, A, B, g, k, c) - q}
    z = tryCatch(
        stats::uniroot(toroot, interval=c(-5,5), extendInt = "upX", check.conv=TRUE)$root,
        error=function(cond) {
            Inf*(q-A) ##If uniroot fails to converge assume root is at +/- Inf
        })
    if (zscale) {
        return(z)
    }
    return(stats::pnorm(z))
}

#' @rdname g-and-k
#' @export
pgk = function (q, A, B, g, k, c=0.8, zscale=FALSE) {
    mapply(pgk_scalar, q, A, B, g, k, c, zscale)
}
