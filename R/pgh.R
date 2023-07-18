#' Distribution function for the g-and-h distribution
#'
#' @keywords internal
#' 
#' @param q Quantiles.
#' @param A A (location) parameter.
#' @param B B (scale) parameter. Must be positive.
#' @param g g parameter.
#' @param h h parameter. Must be non-negative.
#' @param c c parameter. Often fixed at 0.8 which is the default.
#' @param zscale When true returns the N(0,1) quantile of the cdf (needed by dgh).
#' @param type Can be "generalised" (default) or "tukey".
#' @details
#'  This internal function performs the calculation assuming scalar inputs.
#'  The exported function is a vectorised wrapper of this.
#' @return The cumulative probability.
pgh_scalar = function(q, A, B, g, h, c=0.8, zscale=FALSE, type=c("generalised", "tukey")) {
    type = match.arg(type)
    toroot = function(p) {z2gh(p, A, B, g, h, c, type) - q}
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

#' @rdname g-and-h
#' @export
pgh = function (q, A, B, g, h, c=0.8, zscale=FALSE, type=c("generalised", "tukey")) {
    type = match.arg(type)
    mapply(pgh_scalar, q, A, B, g, h, c, MoreArgs=list(zscale=zscale, type=type))
}
