#' Improper uniform log density
#'
#' Returns log density of an improper prior for the g-and-k or g-and-h distribution
#'
#' @param theta A vector of 4 parameters representing (A,B,g,k) or (A,B,g,h)
#' @return Value of an (unnormalised) log density
#' @details
#' \code{improper_uniform_log_density} takes a 4 parameter vector as input and returns a log density value.
#' The output corresponds to an improper uniform with constraints that the second and fourth parameters should be non-negative.
#' These ensure that the resulting parameters are valid to use in the g-and-k or g-and-h distribution is valid.
#' This function is supplied as a convenient default prior to use in the \code{mcmc} function.
#' @examples
#' improper_uniform_log_density(c(0,1,0,0)) ##Valid parameters - returns 0
#' improper_uniform_log_density(c(0,-1,0,0)) ##Invalid parameters - returns -Inf
#' @export
improper_uniform_log_density = function(theta) {
    if (theta[2]<0 || theta[4]<0) return(-Inf)
    return(0)
}

#' Markov chain Monte Carlo inference
#'
#' Markov chain Monte Carlo (MCMC) inference for the g-and-k or g-and-h distribution
#'
#' @param x Vector of observations.
#' @param N Number of MCMC steps to perform.
#' @param model Whether to fit g-and-k or g-and-h model.
#' @param logB When true, the second parameter is log(B) rather than B.
#' @param get_log_prior A function with one argument (corresponding to a vector of 4 parameters e.g. A,B,g,k) returning the log prior density. This should ensure the parameters are valid - i.e. return -Inf for invalid parameters - as the \code{mcmc} code will not check this.
#' @param theta0 Vector of initial value for 4 parameters.
#' @param Sigma0 MCMC proposal covariance matrix
#' @param t0 Tuning parameter (number of initial iterations without adaptation).
#' @param epsilon Tuning parameter (weight given to identity matrix in covariance calculation).
#' @param silent When \code{FALSE} (the default) a progress bar and intermediate results plots are shown.
#' @param plotEvery How often to plot the results if \code{silent==FALSE}.
#' @details \code{mcmc} performs Markov chain Monte Carlo inference for iid data from a g-and-k or g-and-h distribution, using the adaptive Metropolis algorithm of Haario et al (2001).
#' @return Matrix whose rows are MCMC states: the initial state \code{theta0} and N subsequent states.
#' @references
#' D. Prangle gk: An R package for the g-and-k and generalised g-and-h distributions, 2017.
#' H. Haario, E. Saksman, and J. Tamminen. An adaptive Metropolis algorithm. Bernoulli, 2001.
#' @examples
#' set.seed(1)
#' x = rgk(10, A=3, B=1, g=2, k=0.5) ##An unusually small dataset for fast execution of this example
#' out = mcmc(x, N=1000, theta0=c(mean(x),sd(x),0,0), Sigma0=0.1*diag(4))
#' @export
mcmc = function(x, N, model=c("gk", "gh"), logB=FALSE, get_log_prior=improper_uniform_log_density, theta0, Sigma0, t0=100, epsilon=1E-6, silent=FALSE, plotEvery=100) {
    if (!silent) { oldask = par(ask=FALSE) } ##Don't ask before progress plots
    output = matrix(nrow=N+1, ncol=4)
    colnames(output) = c("A", ifelse(logB, "log B", "B"), "g", ifelse(model[1]=="gk", "k", "h"))
    if (model[1] == "gk") {
        if (logB) {
            get_log_likelihood = function(theta) sum(dgk(x, theta[1], exp(theta[2]), theta[3], theta[4], log=TRUE))
        } else {
            get_log_likelihood = function(theta) sum(dgk(x, theta[1], theta[2], theta[3], theta[4], log=TRUE))
        }
    } else {
        if (logB) {
            get_log_likelihood = function(theta) sum(dgh(x, theta[1], exp(theta[2]), theta[3], theta[4], log=TRUE))
        } else {
            get_log_likelihood = function(theta) sum(dgh(x, theta[1], theta[2], theta[3], theta[4], log=TRUE))
        }
    }
    output[1,] = theta0
    theta = theta0
    Sigma = Sigma0
    C = chol(Sigma0)
    theta_bar = 0*theta0 ##Mean value of theta
    theta_mom2 = 0*Sigma ##2nd moment of theta
    if (!silent) { prog_bar = progress::progress_bar$new(total = N+1, format = "[:bar] :percent eta: :eta") }
    log_prior = get_log_prior(theta)
    log_likelihood = get_log_likelihood(theta)
    if (!silent) { prog_bar$tick() }
    for (i in 1:N) {
        theta_bar = theta_bar*(i-1)/i + theta/i
        theta_mom2 = theta_mom2*(i-1)/i + theta%*%t(theta)/i
        if (i > t0) {
            C = chol(sd_const*(theta_mom2 - theta_bar %*% t(theta_bar) + epsilon*diag(4)))
        }
        theta_prop = theta + C %*% stats::rnorm(4)
        log_prior_prop = get_log_prior(theta_prop)
        if (log_prior_prop > -Inf) { ##Skip following if proposal has zero prior
            log_likelihood_prop = get_log_likelihood(theta_prop)
            r = log_prior_prop + log_likelihood_prop - log_prior - log_likelihood
            if (stats::runif(1) < exp(r)) {
                theta = theta_prop
                log_prior = log_prior_prop
                log_likelihood = log_likelihood_prop
            }
        }
        output[i+1,] = theta
        if (!silent && ((i+1) %% plotEvery == 0)) {
            graphics::par(mfrow=c(2,2))
            for (j in 1:4) {
                ylim = range(output[ceiling(i/10):(i+1),j])
                graphics::plot(output[,j], type='s', xlim=c(1,N), ylim=ylim, xlab="MCMC iteration", ylab=colnames(output)[j])
            }
        }
        if (!silent) { prog_bar$tick() }
    }
    if (!silent) { par(ask=oldask) }
    output
}
