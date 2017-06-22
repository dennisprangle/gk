#' Project into a region
#'
#' Project a vector elementwise into a constrained region
#'
#' @keywords internal
#'
#' @param x Vector
#' @param xmin Vector of lower bounds
#' @param xmax Vector of upper bounds
#'
#' @return Vector of closest values to x which satisfy the bounds
project = function(x, xmin, xmax) {
    pmin(pmax(x, xmin), xmax)
}

#' Calculate log sum exp safely
#'
#' Calculate log(exp(a)+exp(b)) avoiding numerical issues
#'
#' @keywords internal
#'
#' @param a Number
#' @param b Number
#'
#' @return Value of log(exp(a)+exp(b))
logSumExp = function(a, b) {
    A = max(a,b)
    B = min(a,b)
    A + log1p(exp(B-A))
}

#' Finite difference stochastic approximation inference
#'
#' Finite difference stochastic approximation (FDSA) inference for the g-and-k or g-and-h distribution
#'
#' @param x Vector of observations.
#' @param N number of iterations to perform.
#' @param model Whether to fit gk or gh model.
#' @param logB When true, the second parameter is log(B) rather than B.
#' @param theta0 Vector of initial value for 4 parameters.
#' @param batch_size Mini-batch size.
#' @param alpha Gain decay for step size.
#' @param gamma Gain decay for finite difference.
#' @param a0 Multiplicative step size tuning parameter (or vector of 4 values).
#' @param c0 Multiplicative finite difference step tuning parameter (or vector of 4 values).
#' @param A Additive step size tuning parameter.
#' @param theta_min Vector of minimum values for each parameter.
#' @param theta_max Vector of maximum values for each parameter.
#' @param silent When \code{FALSE} (the default) a progress bar and intermediate results plots are shown.
#' @param plotEvery How often to plot the results if \code{silent==FALSE}.
#' @details \code{fdsa} performs maximum likelihood inference for iid data from a g-and-k or g-and-h distribution, using simulataneous perturbation stochastic approximation. This should be faster than directly maximising the likelihood.
#' @return Matrix whose rows are FDSA states: the initial state \code{theta0} and N subsequent states.
#' The final row is the MLE estimate.
#' @references
#' D. Prangle gk: An R package for the g-and-k and generalised g-and-h distributions, 2017.
#' @examples
#' set.seed(1)
#' x = rgk(10, A=3, B=1, g=2, k=0.5) ##An unusually small dataset for fast execution of this example
#' out = fdsa(x, N=100, theta0=c(mean(x),sd(x),0,0), theta_min=c(-5,1E-5,-5,0), theta_max=c(5,5,5,5))
#' @export
fdsa = function(x, N, model=c("gk", "gh"), logB=FALSE, theta0, batch_size=100, alpha=1, gamma=0.49, a0=1, c0=NULL, A=100,
    theta_min=c(-Inf,ifelse(logB, -Inf, 1E-5),-Inf,0), theta_max=c(Inf,Inf,Inf,Inf),
    silent=FALSE, plotEvery=100) {
    if (!silent) { oldask = par(ask=FALSE) } ##Don't ask before progress plots
    cnames = c("A", ifelse(logB, "log B", "B"), "g", ifelse(model[1]=="gk", "k", "h"), "estimated log likelihood")
    if (model[1] == "gk") {
        if (logB) {
            get_log_densities = function(x, theta) dgk(batch, theta[1], exp(theta[2]), theta[3], theta[4], log=TRUE)
        } else {
            get_log_densities = function(x, theta) dgk(batch, theta[1], theta[2], theta[3], theta[4], log=TRUE)
        }
    } else {
        if (logB) {
            get_log_densities = function(x, theta) dgh(batch, theta[1], exp(theta[2]), theta[3], theta[4], log=TRUE)
        } else {
            get_log_densities = function(x, theta) dgh(batch, theta[1], theta[2], theta[3], theta[4], log=TRUE)
        }
    }
    nobs = length(x)
    batch_size = min(batch_size, nobs)
    nm_ratio = nobs/batch_size ##Multiplicative constant for later
    theta = theta0
    estimates = matrix(nrow=N+1, ncol=5)
    colnames(estimates) = cnames
    estimates[1,] = c(theta, NA)
    if (missing(c0)) {
        batch = x[sample(nobs, 100, replace=TRUE)] ##Allow replacement here just in case sample size small
        density_sample = get_log_densities(batch, theta0)
        c0 = stats::sd(density_sample) / sqrt(batch_size)
        c0 = pmin(c0, (theta_max - theta_min)/2)
    } else if (any(c0 > (theta_max-theta_min)/2)) {
        stop("c0 too large compared to parameter constraints")
    }
    if (!silent) { prog_bar = progress::progress_bar$new(total = N, format = "[:bar] :percent eta: :eta") }
    for (t in seq(0,N-1)) {
        at = a0*(t+1+A)^-alpha
        ct = c0*(t+1)^-gamma
        indices = sample(nobs, batch_size)
        batch = x[indices]
        gt = rep(0,4)        
        for (i in 1:4) {
            delta = rep(0,4)
            delta[i] = 1
            theta1 = project(theta+ct*delta, theta_min, theta_max)
            theta2 = project(theta-ct*delta, theta_min, theta_max)
            hatL1 = -nm_ratio*sum(get_log_densities(batch, theta1))
            hatL2 = -nm_ratio*sum(get_log_densities(batch, theta2))
            if (is.infinite(hatL1) || is.infinite(hatL2)) {
                stop(paste("Log likelihoods too small to calculate! Parameters probably became too extreme. Last values were", toString(theta), ". Try tighter theta_min or theta_max values."))
            }
            gt[i] = (hatL1 - hatL2)/(theta1[i]-theta2[i])
        }
        theta = project(theta-at*gt, theta_min, theta_max)
        estimates[t+2,1:4] = theta
        estimates[t+2,5] = log(2) - logSumExp(hatL1, hatL2) ##Log of mean likelihood estimate
        #cat("Iteration ", t+1, "\n Likelihood estimates ", hatL1, hatL2, "\n Estimated gradient ", gt, "\n Proposed step", at*gt, "\n New estimate ", theta, "\n\n")
        if (!silent && ((t+1) %% plotEvery == 0)) {
            graphics::par(mfrow=c(2,3))
            for (i in 1:5) {
                ylim = range(estimates[ceiling(t/10):(t+1),i])
                graphics::plot(estimates[,i], type='l', xlim=c(1,N), ylim=ylim, xlab="FDSA iteration", ylab=cnames[i])
            }
        }
        if (!silent) { prog_bar$tick() }
    }
    if (!silent) { par(ask=oldask) }
    return(estimates)
}
