% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isValid.R
\name{isValid_scalar}
\alias{isValid_scalar}
\title{Check validity of g-and-k or g-and-h parameters}
\usage{
isValid_scalar(
  g,
  k_or_h,
  c = 0.8,
  model = c("gk", "generalised_gh", "tukey_gh", "gh"),
  initial_z = seq(-1, 1, 0.2)
)
}
\arguments{
\item{g}{A g parameter.}

\item{k_or_h}{A k or h parameter.}

\item{c}{A c parameter.}

\item{model}{Which model to check: "gk", "generalised_gh" or "tukey_gh".
For backwards compatibility, "gh" acts the same as "generalised_gh".}

\item{initial_z}{Vector of initial z values to use in optimisation.}
}
\value{
Logical vector denoting whether each parameter combination is valid
}
\description{
Check whether parameter choices produce a valid g-and-k or g-and-h distribution.
}
\details{
This internal function performs the calculation using scalar parameter inputs.
 The exported function is a vectorised wrapper of this.
}
\keyword{internal}
