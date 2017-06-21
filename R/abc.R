#' Convert octiles to moment estimates
#'
#' Convert octiles to estimates of location, scale, skewness and kurtosis.
#'
#' @param octiles Vector of octiles.
#' @details Converts octiles to robust estimate of location, scale, skewness and kurtosis as used by Drovandi and Pettitt (2011).
#' @references C. C. Drovandi and A. N. Pettitt. Likelihood-free Bayesian estimation of multivariate quantile distributions. Computational Statistics & Data Analysis, 2011.
#' @return Vector of moment estimates.
momentEstimates = function(octiles) {
    momb = octiles[6]-octiles[2]
    c(octiles[4], momb, (octiles[6]+octiles[2]-2*octiles[4])/momb, (octiles[7]-octiles[5]+octiles[3] - octiles[1])/momb)
}

#' Approximate Bayesian computation inference
#'
#' Approximate Bayesian computation (ABC) inference for the g-and-k or g-and-h distribution.
#'
#' @param x Vector of observations.
#' @param N Number of iterations to perform.
#' @param model Whether to fit g-and-k or g-and-h model.
#' @param logB When true, the second parameter is log(B) rather than B.
#' @param rprior A function with single argument, n, which returns a matrix with n rows consisting of samples from the prior distribution for 4 parameters e.g. (A,B,g,k).
#' @param M Number of simulations to accept.
#' @param sumstats Which summary statistics to use.
#' @param silent When \code{FALSE} (the default) a progress bar is shown.
#' @details
#' This function performs approximate Bayesian inference for iid data from a g-and-k or g-and-h distribution, avoiding expensive density calculations.
#' The algorithm samples many parameter vectors from the prior and simulates corresponding data from the model.
#' The parameters are accepted or rejected based on how similar the simulations are to the observed data.
#' Similarity is measured using weighted Euclidean distance between summary vectors of the simulations and observations.
#' Several summaries can be used, including the complete order statistics or summaries based on octiles.
#' In the latter case only the corresponding order statistics are simulated, speeding up the method. 
#' @return Matrix whose rows are accepted parameter estimates plus a column giving the ABC distances.
#' @references
#' D. Prangle. gk: An R package for the g-and-k and generalised g-and-h distributions, 2017.
#' @examples
#' set.seed(1)
#' x = rgk(10, A=3, B=1, g=2, k=0.5) ##An unusually small dataset for fast execution of this example
#' rprior = function(n) { matrix(runif(4*n,0,10), ncol=4) }
#' abc(x, N=1E4, rprior=rprior, M=100)
#' @export
abc = function(x, N, model=c("gk", "gh"), logB=FALSE, rprior, M, sumstats=c("all order statistics", "octiles", "moment estimates"), silent=FALSE) {
    nobs = length(x)
    ##Define simStats: a function to simulate one set of summary statistics
    ##and sobs: the observed summary statistics
    if (sumstats[1] == "all order statistics") {
        sobs = sort(x)
        if (model[1] == "gk") {
            if (logB) {
                simStats = function(theta) sort(rgk(nobs, A=theta[1], B=exp(theta[2]), g=theta[3], k=theta[4]))
            } else {
                simStats = function(theta) sort(rgk(nobs, A=theta[1], B=theta[2], g=theta[3], k=theta[4]))
            }
        } else {
            if (logB) {
                simStats = function(theta) sort(rgh(nobs, A=theta[1], B=exp(theta[2]), g=theta[3], h=theta[4]))
            } else {
                simStats = function(theta) sort(rgh(nobs, A=theta[1], B=theta[2], g=theta[3], h=theta[4]))
            }
        }
    } else {
        indices = round(nobs*(1:7)/8)
        if (length(unique(indices)) < 7) stop("Too few observations to calculate distinct octiles")
        if (sumstats[1] == "octiles") {
            sobs = sort(x)[indices]
            if (model[1] == "gk") {
                if (logB) {
                    simStats = function(theta) qgk(orderstats(nobs, indices), A=theta[1], B=exp(theta[2]), g=theta[3], k=theta[4])
                } else {
                    simStats = function(theta) qgk(orderstats(nobs, indices), A=theta[1], B=theta[2], g=theta[3], k=theta[4])
                }
            } else {
                if (logB) {
                    simStats = function(theta) qgh(orderstats(nobs, indices), A=theta[1], B=exp(theta[2]), g=theta[3], h=theta[4])
                } else {
                    simStats = function(theta) qgh(orderstats(nobs, indices), A=theta[1], B=theta[2], g=theta[3], h=theta[4])
                }
            }
        } else {
            sobs = momentEstimates(sort(x)[indices])
            if (model[1] == "gk") {
                if (logB) {
                    simStats = function(theta) {
                        momentEstimates(qgk(orderstats(nobs, indices), A=theta[1], B=exp(theta[2]), g=theta[3], k=theta[4])) }
                } else {
                    simStats = function(theta) {
                        momentEstimates(qgk(orderstats(nobs, indices), A=theta[1], B=theta[2], g=theta[3], k=theta[4]))
                    }
                }
            } else {
                if (logB) {
                    simStats = function(theta) {
                        momentEstimates(qgh(orderstats(nobs, indices), A=theta[1], B=exp(theta[2]), g=theta[3], h=theta[4])) }
                } else {
                    simStats = function(theta) {
                        momentEstimates(qgh(orderstats(nobs, indices), A=theta[1], B=theta[2], g=theta[3], h=theta[4]))
                    }

                }
            }
        }
    }
    batch_size = 10^4
    if (N <= batch_size) { ##If N is small enough do ABC as a single batch
        samp = abc_batch(sobs, rprior(N), simStats, M)$samp
    } else { ##Otherwise split ABC into batches
        nbatches = ceiling(N / batch_size)
        last_batch_size = N %% batch_size
        if (last_batch_size == 0) { last_batch_size = batch_size }
        if (!silent) { prog_bar = progress::progress_bar$new(total = nbatches, format = "[:bar] :percent eta: :eta") }
        batch_out = abc_batch(sobs, rprior(batch_size), simStats, M)
        samp = batch_out$samp
        v = batch_out$v
        if (!silent) { prog_bar$tick() }
        next_batch_size = batch_size
        for (b in 2:nbatches) {
            if (b==nbatches) { next_batch_size = last_batch_size }
            next_samp = abc_batch(sobs, rprior(next_batch_size), simStats, M, v)$samp
            samp = rbind(samp, next_samp)
            toacc = order(samp[,5])[1:M]
            samp = samp[toacc,]
            if (!silent) { prog_bar$tick() }
        }
    }
    colnames(samp) = c("A", ifelse(logB, "log B", "B"), "g", ifelse(model[1]=="gk", "k", "h"), "distance")
    return(samp)
}

#' Single batch of ABC
#'
#' Interal function carrying out part of ABC inference calculations
#'
#' @param sobs Vector of observed summary statistics.
#' @param priorSims Matrix whose rows are vector of parameters drawn from the prior.
#' @param simStats Function mapping a vector of parameters to summary statistics.
#' @param M How many simulations to accept in this batch.
#' @param v Optional vector of estimated variances for each summary statistic
#' @details
#' Outputs a list containing:
#' 1) matrix whose rows are accepted parameter vectors, and a column of resulting distances
#' 2) value of v.
#' If v is not supplied then variances are estimated here and returned.
#' (The idea is that this is done for the first batch, and these values are reused from then on.)
abc_batch = function(sobs, priorSims, simStats, M, v=NULL) {
    summaries = apply(priorSims, 1, simStats)
    ##Calculate distances
    if (missing(v)) {
        v = apply(summaries, 1, stats::var)
    }
    d = apply(summaries, 2, function(s) { sum(((s-sobs)^2)/v) })
    ##Construct and return output
    toacc = order(d)[1:M]
    samp = cbind(priorSims[toacc,], d[toacc])
    list(samp=samp, v=v)
}
