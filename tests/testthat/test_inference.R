library(gk)
context("Test inference code executes")

test_that("Test abc", {
  set.seed(1)
  x = rgk(10, A=3, B=1, g=2, k=0.5)
  rprior = function(n) { matrix(runif(4*n,0,10), ncol=4) }
  for (model in c("gk", "generalised_gh", "tukey_gh", "gh")) {
    for (sumstats in c("all order statistics", "octiles", "moment estimates")) {
      abc(x, N=1E4, model=model, rprior=rprior, M=100, sumstats=sumstats)
      abc(x, N=1E4, model=model, logB=TRUE, rprior=rprior, M=100, sumstats=sumstats)
    }
  }
})

test_that("Test fdsa", {
  set.seed(1)
  x = rgk(10, A=3, B=1, g=2, k=0.5)
  for (model in c("gk", "generalised_gh", "tukey_gh", "gh")) {
    fdsa(x, N=100, model=model, theta0=c(mean(x),sd(x),0,0), theta_min=c(-5,1E-5,-5,1E-5), theta_max=c(5,5,5,5))
    fdsa(x, N=100, model=model, logB=TRUE, theta0=c(mean(x),sd(x),0,0), theta_min=c(-5,-5,-5,1E-5), theta_max=c(5,5,5,5))
  }
})