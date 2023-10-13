# npparam
Code repository for the paper "Semiparametric Shape-restricted Estimators for Nonparametric Regression". See https://arxiv.org/abs/2307.05732 for details.
To install this `R` package, first install the `devtools` package. Then type:

    library(devtools)
    devtools::install_github("Kenta426/npparam")
    library(npparam)

# Usage
All methods are developed by extending the standard S3 methods.

## Lipschitz function
```{R}
# sample data
n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
res.est <- iso.lin.est(x, y)
print(res.est) # display results
y.pred <- predict(res.est) # fitted value
y.pred <- predict(res.est, newdata=runif(n, -1, 1)) # predict new data
```

## Function with Lipschitz derivative 
```{R}
# sample data
n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
res.est <- cvx.quad.est(x, y)
print(res.est) # display results
y.pred <- predict(res.est) # fitted value
y.pred <- predict(res.est, newdata=runif(n, -1, 1)) # predict new data
```

## Additive function with Lipschitz components
```{R}
# sample data
n <- 100*2; x <- matrix(runif(n*2, -1, 1), ncol=2);
y <- abs(x[,1]) - abs(x[,2]) + rnorm(n, 0, 0.1)
res.est <- additive.iso.est(x, y)
print(res.est) # display results
y.pred <- predict(res.est) # fitted value
y.pred <- predict(res.est, newdata=matrix(runif(n*2, -1, 1), ncol=2)) # predict new data
```

## Single-index model with Lipschitz link function
```{R}
n <- 100; d <- 3; x <- matrix(runif(n*d, -1, 1), nrow=n);
beta <- c(1,0,0); y <- sin(4*x%*%beta) + rnorm(n, 0, 0.1)
res.est <- single.index.est(x, y)
y.pred <- predict(res.est) # fitted value
y.pred <- predict(res.est, newdata=matrix(runif(n*d, -1, 1))) # predict new data
```

