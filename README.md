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

## Lipschitz function over fixed-design 2d lattice based on block max-min algorithm (Deng and Zhang (2020))
```{R}
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
y.hat <- predict(res.est) # estimated values over the grid.
```
