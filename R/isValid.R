#' Check validity of g-and-k or g-and-h parameters
#'
#' Check whether parameter choices produce a valid g-and-k or g-and-h distribution.
#'
#' @keywords internal
#' 
#' @param g A g parameter.
#' @param k_or_h A k or h parameter.
#' @param c A c parameter.
#' @param model Which model to check: "gk", "generalised_gh" or "tukey_gh".
#' For backwards compatibility, "gh" acts the same as "generalised_gh".
#' @param initial_z Vector of initial z values to use in optimisation.
#' @details
#'  This internal function performs the calculation using scalar parameter inputs.
#'  The exported function is a vectorised wrapper of this.
#' @return Logical vector denoting whether each parameter combination is valid
isValid_scalar = function(g, k_or_h, c=0.8, model=c("gk", "generalised_gh", "tukey_gh", "gh"), initial_z = seq(-1,1,0.2)) {
    model = match.arg(model)
    tomin = switch(
        model,
        "gk" = function(z) Qgk_deriv(z, 0, 1, g, k_or_h, c, getR=TRUE),
        "gh" = function(z) Qgh_deriv(z, 0, 1, g, k_or_h, c, getR=TRUE, type="generalised"),
        "generalised_gh" = function(z) Qgh_deriv(z, 0, 1, g, k_or_h, c, getR=TRUE, type="generalised"),
        "tukey_gh" = function(z) Qgh_deriv(z, 0, 1, g, k_or_h, c, getR=TRUE, type="tukey")
    )
    positive_minimum = sapply(initial_z, function(z0) (stats::optim(z0, tomin, method="BFGS")$value > 0)==TRUE)
    return(all(positive_minimum))
}

#' Check validity of g-and-k or g-and-h parameters
#'
#' Check whether parameter choices produce a valid g-and-k or g-and-h distribution.
#'
#' @param g Vector of g parameters.
#' @param k_or_h Vector of k or h parameters.
#' @param c Vector of c parameters.
#' @param model Which model to check: "gk", "generalised_gh" or "tukey_gh".
#' For backwards compatibility, "gh" acts the same as "generalised_gh".
#' @param initial_z Vector of initial z values to use in each optimisation.
#' @details
#' This function tests whether parameter choices provide a valid distribution.
#' Only g k and c parameters need be supplied as A and B>0 have no effect.
#' The function operates by numerically minimising the derivative of the quantile function,
#' and returning \code{TRUE} if the minimum is positive.
#' It is possible that a local minimum is found, so it is recommended to use multiple optimisation starting points, and to beware that false positive may still result!
#' @return Logical vector denoting whether each parameter combination is valid.
#' @examples
#' isValid(0:10, -0.5)
#' isValid(0:10, 0.5, c=0.9, model="generalised_gh")
#' isValid(0:10, 0.5, model="tukey_gh")
#' @export
isValid = function(g, k_or_h, c=0.8, model=c("gk", "generalised_gh", "tukey_gh", "gh"), initial_z = seq(-1,1,0.2)) {
    model = match.arg(model)
    mapply(isValid_scalar, g, k_or_h, c, MoreArgs=list(model=model, initial_z=initial_z))
}
