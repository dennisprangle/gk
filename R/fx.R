#' Exchange rate example
#'
#' Performs a demo analysis of exchange rate data
#'
#' @param type What type of demo to perform. \code{standard} runs the full demo displaying plots. \code{quick} is a fast demo suitable for testing the package. \code{for paper} saves pdfs which can be used in the paper.
#' @details
#' \code{fx} performs the analysis of exchange rate data in the supporting reference (Prangle 2017).
#' @references
#' D. Prangle gk: An R package for the g-and-k and generalised g-and-h distributions, 2017.
#' @examples
#' \dontrun{
#' fx()
#' }
#' @import grDevices
#' @import graphics
#' @import stats
#' @export
fx = function(type=c("standard", "quick", "for paper")) {
    type = type[1]
    if (type == "standard") par(ask=TRUE)
    if (type == "quick") {
        Nscale = 1/100 ##Do this proportion of normal number of iterations in a quick run
    } else {
        Nscale = 1
    }
    
    ##Import data
    fx_rate = Ecdat::Garch$cd
    n = length(fx_rate)
    log_return = log(fx_rate[2:n] / fx_rate[1:(n-1)])

    ##Plot data as time series
    fx_date = lubridate::ymd(Ecdat::Garch$date[2:n])
    if (type == "for paper") pdf(file="fx_time_series.pdf", width=8, height=6)
    plot(fx_date, log_return, pch=16, xlab='Date', ylab='Log return')
    if (type == "for paper") dev.off()

    ##Preliminary ABC analysis
    if (type != "for paper") cat("ABC analysis\n")
    set.seed(1)
    rprior = function(i) { cbind(runif(i,-1,1), runif(i,0,1), runif(i,-5,5), runif(i,0,10)) }
    t0 = Sys.time()
    abc_out = abc(log_return, N=1E7*Nscale, rprior=rprior, M=200, sumstats='moment estimates', silent=(type == "for paper"))
    tabc = Sys.time() - t0
    cat("ABC took", tabc, attr(tabc, "units"), "\n")
    abc_est = apply(abc_out[,1:4], 2, mean)
    abc_out_tf = abc_out[,1:4]
    abc_out_tf[,2] = log(abc_out_tf[,2])
    abc_est_tf = colMeans(abc_out_tf)

    ##FDSA analysis
    if (type != "for paper") cat("FDSA analysis\n")
    a0_pilot = 2E-4
    fdsa_out_pilot = fdsa(log_return, N=1E4*Nscale, logB=TRUE, theta0=abc_est_tf, batch_size=100, a0=a0_pilot, silent=(type == "for paper"), plotEvery=100)

    a0 = c(1E-6, 1E-2, 1E-2, 1E-2)
    t0 = Sys.time()
    fdsa_out = fdsa(log_return, N=1E4*Nscale, logB=TRUE, theta0=abc_est_tf, batch_size=100, a0=a0, silent=(type == "for paper"), plotEvery=100)
    tfdsa = Sys.time() - t0
    cat("FDSA took", tfdsa, attr(tfdsa, "units"), "\n")
    fdsa_est_tf = fdsa_out[nrow(fdsa_out),1:4]
    fdsa_est = fdsa_est_tf
    fdsa_est[2] = exp(fdsa_est[2])

    if (type == "for paper") {
        pdf(file="fdsa_trace.pdf", width=8, height=6)
        par(mfrow=c(2,2), mar=c(5,4,1,1)+0.1)
        for (i in 1:4) {
            plot(fdsa_out_pilot[,i], type='s', xlab='Iteration', ylab=colnames(fdsa_out)[i],
                 ylim=range(c(fdsa_out_pilot[,i], fdsa_out[,i])))
            lines(fdsa_out[,i], type='s', col='red')
        }
        dev.off()
    }

    ##MCMC analysis
    if (type != "for paper") cat("MCMC analysis\n")
    if (type == "quick") {
        Sigma0 = var(fdsa_out[,1:4])
    } else {
        Sigma0 = var(fdsa_out[1E4 + (-1000:0),1:4])
    }

    log_prior = function(theta) {
        ##The prior is e^(log B) 1(k>0)
        ##So the log prior is log B for k>0 and -inf otherwise
        ##This is equivalent to a uniform prior in the original parameterisation
        if (theta[4]<0) return(-Inf)
        return(theta[2])
    }

    t0 = Sys.time()
    mcmc_out_tf = mcmc(log_return, N=1E4*Nscale, logB=TRUE, get_log_prior=log_prior, theta0=fdsa_est_tf, Sigma0=Sigma0, silent=(type == "for paper"))
    tmcmc = Sys.time() - t0
    cat("MCMC took", tmcmc, attr(tmcmc, "units"), "\n")

    if (type == "for paper") {
        pdf(file="mcmc_trace.pdf", width=8, height=6)
        par(mfrow=c(2,2), mar=c(5,4,1,1)+0.1)
        for (i in 1:4) plot(mcmc_out_tf[,i], type='s', xlab='Iteration', ylab=colnames(mcmc_out_tf)[i])
        dev.off()
    }

    if (type == "quick") {
        mcmc_out = mcmc_out_tf
    } else {
        mcmc_out = mcmc_out_tf[5000:10001,1:4]
    }
    mcmc_out[,2] = exp(mcmc_out[,2])
    colnames(mcmc_out)[2] = 'B'
    if (type != "for paper") {
        par(mfrow=c(2,2))
        for (i in 1:4) hist(mcmc_out[,i], xlab=colnames(mcmc_out)[i], main='')
    }
    mcmc_est = apply(mcmc_out, 2, mean)

    ##Plot posteriors
    if (type != "for paper") cat("Output plots\n")
    if (type == "for paper") pdf(file="posteriors.pdf", width=8, height=6)
    par(mfcol=c(2,4), mar=c(5,4,1,1)+0.1)
    for (i in 1:4) {
        hist(abc_out[,i], freq=FALSE, xlab=colnames(abc_out)[i], main='')
        points(x=fdsa_est[i], y=0, col='red', pch='X', cex=2)
        if (i == 4) mtext("ABC", side=4)
        hist(mcmc_out[,i], freq=FALSE, xlab=colnames(abc_out)[i], main='')
        points(x=fdsa_est[i], y=0, col='red', pch='X', cex=2)
        if (i == 4) mtext("MCMC", side=4)
    }
    if (type == "for paper") dev.off()

    ##Bivariate plots
    if (type != "for paper") {
        pairs(abc_out)
        pairs(mcmc_out)
    }

    ##Investigate model fit
    if (type == "for paper") {
        pdf("gk_fits.pdf", width=8, height=6)
        par(mfrow=c(1,2))
    } else {
        par(mfrow=c(1,1))
    }

    hist(log_return, freq=FALSE, breaks=60, main='', xlab='Log return')
    x = seq(min(log_return), max(log_return), length.out=1000)
    lines(x, dgk(x, abc_est[1], abc_est[2], abc_est[3], abc_est[4]), col='red', lwd=3)
    lines(x, dgk(x, mcmc_est[1], mcmc_est[2], mcmc_est[3], mcmc_est[4]), col='blue', lty=2, lwd=3)
    lines(x, dgk(x, fdsa_est[1], fdsa_est[2], fdsa_est[3], fdsa_est[4]), col='green', lty=3, lwd=3)
    norm_est = c(mean(log_return), var(log_return)*(n-1)/(n-2))
    lines(x, dnorm(x, norm_est[1], sqrt(norm_est[2])), col='black', lwd=3)
    legend('topright', legend=c('Normal', 'ABC', 'FDSA', 'MCMC'), lty=c(1,1,3,2), lwd=3, col=c('black','red','green','blue'), bty='n')

    sample_quantiles = sort(log_return)
    m = length(sample_quantiles)
    theoretical_quantiles_norm = qnorm((1:m)/(m+1), norm_est[1], sqrt(norm_est[2]))
    plot(sample_quantiles, theoretical_quantiles_norm, ylim=range(sample_quantiles), pch=1, cex=0.5,
         xlab='Sample quantile', ylab='Theoretical quantile')
    for (i in 1:30) {
        j = sample(nrow(mcmc_out), 1)
        theoretical_quantiles_mcmc = qgk((1:m)/(m+1), mcmc_out[j,1], mcmc_out[j,2], mcmc_out[j,3], mcmc_out[j,4])
        points(sample_quantiles, theoretical_quantiles_mcmc, col='blue', pch=4, cex=0.5)
    }
    theoretical_quantiles_abc = qgk((1:m)/(m+1), abc_est[1], abc_est[2], abc_est[3], abc_est[4])
    points(sample_quantiles, theoretical_quantiles_abc, col='red', pch=2, cex=0.5)
    theoretical_quantiles_fdsa = qgk((1:m)/(m+1), fdsa_est[1], fdsa_est[2], fdsa_est[3], fdsa_est[4])
    points(sample_quantiles, theoretical_quantiles_fdsa, col='green', pch=3, cex=0.5)
    points(sample_quantiles, theoretical_quantiles_norm, ylim=range(sample_quantiles), pch=1, cex=0.5)
    abline(a=0,b=1,lty=2)
    legend('topleft', legend=c('Normal', 'ABC', 'FDSA', 'MCMC'), pch=c(1,2,3,4), col=c('black','red','green','blue'), bty='n')
    
    if (type == "for paper") dev.off()
}
