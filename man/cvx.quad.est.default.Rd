% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvx.quad.est.R
\name{cvx.quad.est.default}
\alias{cvx.quad.est.default}
\title{Nonparametric LSE for twice differentiable  function}
\usage{
\method{cvx.quad.est}{default}(x, y, L = NULL, run.optim = TRUE, ...)
}
\arguments{
\item{x}{\code{n x 1} numeric vector of observed covariate vector}

\item{y}{\code{n x 1} numeric vector of observed response vector}

\item{L}{A scalar for the quadratic parameter. If not specified, L is
data-adaptively selected via optimization or cross-validation}

\item{run.optim}{A Boolean parameter for selecting optimization or
cross-validation. If TRUE, it performs optimization. Default is TRUE.}

\item{...}{Additional control parameters}
}
\value{
An object with S3 class \code{cvx.quad.est}.
}
\description{
Nonparametric LSE for twice differentiable  function
}
\examples{
n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
res.est <- cvx.quad.est(x, y)
}
