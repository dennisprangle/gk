#' Check validity of g-and-k or g-and-h parameters
#'
#' Check whether parameter choices produce a valid g-and-k or g-and-h distribution.
#'
#' @keywords internal
#' 
#' @param g A g parameter.
#' @param k_or_h A k or h parameter.
#' @param c A c parameter.
#' @param model Whether to check the g-and-k or g-and-h model.
#' @param initial_z Vector of initial z values to use in optimisation.
#' @details
#'  This internal function performs the calculation using scalar parameter inputs.
#'  The exported function is a vectorised wrapper of this.
#' @return Logical vector denoting whether each parameter combination is valid
isValid_scalar = function(g, k_or_h, c=0.8, model=c("gk","gh"), initial_z = seq(-1,1,0.2)) {
    if (model[1] == "gk") {
        tomin = function(z) Qgk_deriv(z, 0, 1, g, k_or_h, c, getR=TRUE)
    } else {
        tomin = function(z) Qgh_deriv(z, 0, 1, g, k_or_h, c, getR=TRUE)
    }
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
#' @param model Whether to check the g-and-k or g-and-h model.
#' @param initial_z Vector of initial z values to use in each optimisation.
#' @details
#' This function tests whether parameter choices provide a valid distribution.
#' Only g k and c parameters need be supplied as A and B>0 have no effect.
#' The function operates by numerically minimising the derivative of the quantile function,
#' and returning \code{TRUE} if the minimum is positive.
#' It is possible that a local minimum is found, so it is recommended to use multiple optimisation starting points, and to beware that false positive may still result!
#' @return Logical vector denoting whether each parameter combination is valid.
#' @references D. Prangle and K. Peel. gk: An R package for the g-and-k and g-and-h distributions, in preparation.
#' @examples
#' isValid(0:10, -0.5)
#' isValid(0:10, 0.5, c=0.9, model="gh")
#' @export
isValid = function(g, k_or_h, c=0.8, model=c("gk","gh"), initial_z = seq(-1,1,0.2)) {
    mapply(isValid_scalar, g, k_or_h, c, MoreArgs=list(model=model, initial_z=initial_z))
}
