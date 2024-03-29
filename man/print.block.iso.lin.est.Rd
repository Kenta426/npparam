% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/block.iso.lin.est.R
\name{print.block.iso.lin.est}
\alias{print.block.iso.lin.est}
\title{Inherit print function for \code{block.iso.lin.est}}
\usage{
\method{print}{block.iso.lin.est}(obj, ...)
}
\arguments{
\item{obj}{An object with S3 class \code{block.iso.lin.est}.}

\item{...}{Additional control parameter}
}
\description{
Inherit print function for \code{block.iso.lin.est}
}
\examples{
n1 <- 10; n2 <- 10
y <- matrix(rep(0, n1*n2), nrow = n1)
x <- matrix(rep(list(), n1*n2),nrow = n1, ncol =n2)
for (i in 1:n1){
  for (j in 1:n2){
    y[i,j] <- sin(10*norm(c(i/n1, j/n2), "2")) + rnorm(1, 0, 0.1)
    x[i,j] <- list(c(i/n1, j/n2))
  }
}
res.est <- block.iso.lin.est(x, y)
print(res.est)
}
