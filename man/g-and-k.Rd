% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgk.R, R/pgk.R, R/qgk.R, R/rgk.R
\name{dgk}
\alias{dgk}
\alias{pgk}
\alias{g-and-k}
\alias{qgk}
\alias{rgk}
\title{g-and-k distribution functions}
\usage{
dgk(x, A, B, g, k, c = 0.8, log = FALSE)

pgk(q, A, B, g, k, c = 0.8, zscale = FALSE)

qgk(p, A, B, g, k, c = 0.8)

rgk(n, A, B, g, k, c = 0.8)
}
\arguments{
\item{x}{Vector of quantiles.}

\item{A}{Vector of A (location) parameters.}

\item{B}{Vector of B (scale) parameters. Must be positive.}

\item{g}{Vector of g parameters.}

\item{k}{Vector of k parameters. Must be at least -0.5.}

\item{c}{Vector of c parameters. Often fixed at 0.8 which is the default.}

\item{log}{If true the log density is returned.}

\item{q}{Vector of quantiles.}

\item{zscale}{If true the N(0,1) quantile of the cdf is returned.}

\item{p}{Vector of probabilities.}

\item{n}{Number of draws to make.}
}
\value{
\code{dgk} gives the density, \code{pgk} gives the distribution, \code{qgk} gives the quantile function, and \code{rgk} generates random deviates
}
\description{
Density, distribution function, quantile function and random generation for the g-and-k distribution.
}
\details{
The g-and-k distribution is defined by its quantile function:
 \deqn{x(p) = A + B [1 + c \tanh(gz/2)] z(1 + z^2)^k,}{x(p) = A + B [1 + c tanh(gz/2)] z (1 + z^2)^k,}
 where z is the standard normal quantile of p.
 Parameter restrictions include \eqn{B>0} and \eqn{k \geq -0.5}{k >= -0.5}. Typically c=0.8. For more
 background information see the references.

\code{rgk} and \code{qgk} use quick direct calculations. However \code{dgk} and \code{pgk} involve slower numerical inversion of the quantile function.

Especially extreme values of the inputs will produce \code{pgk} output rounded to 0 or 1 (-Inf or Inf for \code{zscale=TRUE}).
The corresponding \code{dgk} output will be 0 or -Inf for \code{log=TRUE}.
}
\examples{
p = 1:9/10 ##Some probabilities
x = qgk(seq(0.1,0.9,0.1), A=3, B=1, g=2, k=0.5) ##g-and-k quantiles
rgk(5, A=3, B=1, g=2, k=0.5) ##g-and-k draws
dgk(x, A=3, B=1, g=2, k=0.5) ##Densities of x under g-and-k
dgk(x, A=3, B=1, g=2, k=0.5, log=TRUE) ##Log densities of x under g-and-k
pgk(x, A=3, B=1, g=2, k=0.5) ##Distribution function of x under g-and-k
}
\references{
Haynes `Flexible distributions and statistical models in ranking and selection procedures, with applications' PhD Thesis QUT (1998)
Rayner and MacGillivray `Numerical maximum likelihood estimation for the g-and-k and generalized g-and-h distributions' Statistics and Computing, 12, 57-75 (2002)
}
