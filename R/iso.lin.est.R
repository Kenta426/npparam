# Estimating L-Lipschitz function via isotonic regression + linear function.

#' Nonparametric LSE for L-Lipschitz function
#'
#' @param x \code{n x 1} numeric vector of the observed covariate vector.
#' @param y \code{n x 1} numeric vector of the observed response vector.
#' @param L A scalar for the linear parameter. If not specified, L is
#'  data-adaptively selected via numerical optimization.
#' @param L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 x log10(n), L0 x log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
#' @param run.cross.fit A Boolean parameter for selecting sample-splitting or
#'  cross-fitting If TRUE, it performs cross-fitting Default is TRUE.
#' @param rep When rep > 1, the sample-splitting is performed multiple times.
#'  This procedure can reduce the randomness due to the index of sample
#'  splitting.
#' @param ... Additional control parameters
#'
#' @return An object with S3 class \code{iso.lin.est}.
#' @export
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- iso.lin.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=runif(n, -1, 1)) # predict new data
iso.lin.est <- function(x, y, L=NULL, L0=100, run.cross.fit=TRUE, rep=1,...){
  UseMethod("iso.lin.est")
}

#' Helper function for fitting isotonic function
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector.
#' @param y \code{n x 1} numeric vector of observed response vector.
#' @param L A scalar for the linear parameter.
.fit.iso.reg <- function(x, y, L){
  # Use the min-max formula to find non-decreasing LSE.
  iso_final.incr <- isoreg(x, y - L*x)
  iso_model.incr <- stepfun(x=x, y=c(iso_final.incr$yf[1], iso_final.incr$yf))
  incr.risk <- mean((y - iso_model.incr(x) - L*x)^2)
  # Use the min-max formula to find non-increasing LSE.
  iso_final.decr <- isoreg(x, y + L*x)
  iso_model.decr <- stepfun(x = x, y=c(iso_final.incr$yf[1], iso_final.incr$yf))
  decr.risk <- mean((y - iso_model.decr(x) + L*x)^2)
  # pick estimator with smaller empirical risk.
  if(decr.risk < incr.risk){
    iso_model <- iso_model.decr; risk <- decr.risk
  }
  else{
    iso_model <- iso_model.incr; risk <- incr.risk
  }
  return(list(reg=iso_model, risk=risk))
}

#' Helper function for running an optimizer for the linear parameter
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector.
#' @param y \code{n x 1} numeric vector of observed response vector.
#' @param L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 x log10(n), L0 x log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
.run.optimizer.iso <- function(x, y, L0=100){
  # split data into train/validation
  n <- length(x); learn.n <- as.integer(n-1/2*(n))
  train.id <- sort(sample(seq_len(n), size=learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  sim.one <- function(l){
    reg.tr <- .fit.iso.reg(x.train, y.train, L=l)
    y.hat <- reg.tr$reg(x.val) + l*x.val
    mean((y.val-y.hat)^2)
  }
  # run optimizer to find the best L parameter
  opt.res <- optimise(sim.one, lower=-L0*log10(n), upper=L0*log10(n),
                      maximum = FALSE)
  # return L with the smallest validation error
  return(list(l.opt=opt.res$minimum, risk=opt.res$objective))
}

#' Helper function for running an optimizer
#' for the linear parameter with crossfitting
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector.
#' @param y \code{n x 1} numeric vector of observed response vector.
#' @param L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 x log10(n), L0 x log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
.run.optimizer.iso.cross.fit <- function(x, y, L0=100){
  # split data into train/validation
  n <- length(x); learn.n <- as.integer(n-1/2*(n))
  train.id <- sort(sample(seq_len(n), size=learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  # first compute one model
  sim.one <- function(l){
    reg.tr <- .fit.iso.reg(x.train, y.train, L=l)
    y.hat <- reg.tr$reg(x.val) + l*x.val
    mean((y.val-y.hat)^2)
  }
  # run optimizer to find the best L parameter
  opt.res.1 <- optimise(sim.one, lower=-L0*log10(n), upper=L0*log10(n),
                      maximum = FALSE)
  # then compute another model after swapping the role of train/val sets
  sim.one <- function(l){
    reg.tr <- .fit.iso.reg(x.val, y.val, L=l)
    y.hat <- reg.tr$reg(x.train) + l*x.train
    mean((y.train-y.hat)^2)
  }
  # run optimizer to find the best L parameter
  opt.res.2 <- optimise(sim.one, lower=-L0*log10(n), upper=L0*log10(n),
                      maximum = FALSE)
  # return two models
  return(list(l.opt=c(opt.res.1$minimum, opt.res.2$minimum),
              risk=c(opt.res.1$objective, opt.res.2$objective)))
}

#' Nonparametric LSE for Lipschitz function
#'
#' @param x \code{n x 1} numeric vector of the observed covariate vector.
#' @param y \code{n x 1} numeric vector of the observed response vector.
#' @param L A scalar for the linear parameter. If not specified, L is
#'  data-adaptively selected via numerical optimization.
#' @param L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 x log10(n), L0 x log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
#' @param run.cross.fit A Boolean parameter for selecting sample-splitting or
#'  cross-fitting If TRUE, it performs cross-fitting Default is TRUE.
#' @param rep When rep > 1, the sample-splitting is performed multiple times.
#'  This procedure can reduce the randomness due to the index of sample
#'  splitting.
#' @param ... Additional control parameters
#'
#' @return An object with S3 class \code{iso.lin.est}.
#' @export
#' @examples
#' n <- 100
#' x <- runif(n, -1, 1)
#' y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- iso.lin.est(x, y)
iso.lin.est.default <- function(x, y, L=NULL, L0=100, run.cross.fit=TRUE, rep=1,...){
  # TODO: check input format
  idx <- order(x); x <- x[idx]; y <- y[idx]
  # select parameter selection procedure
  if(run.cross.fit){
    param.selector <- .run.optimizer.iso.cross.fit; optimizer <- "Cross-fit"
  }else{
    param.selector <- .run.optimizer.iso; optimizer <- "Optimizer"
  }
  # run parameter selection if L is not specified
  if(is.null(L)){
    L <- param.selector(x, y, L0=L0)$l.opt
    if (run.cross.fit & (rep > 1)){
      # repeated sampling
      for (i in 1:(ceiling(rep)-1)){
        L <- c(L, param.selector(x, y)$l.opt)
      }
    }
  }
  # fit isotonic regression
  if(length(L) == 1){
    est.res <- .fit.iso.reg(x, y, L)
    iso_model <- est.res$reg; risk <- est.res$risk
    res <- list(x.values=x[order(idx)], y.values=y[order(idx)], l.value=L,
                g=iso_model, risk=risk, optimizer=optimizer)
  }
  else{
    iso_model <- c(); risk <- 0
    # for multiple Ls, create list of estimators for each L
    for (l in L){
      est.res <- .fit.iso.reg(x, y, l)
      iso_model <- c(iso_model, est.res$reg); risk <- est.res$risk/length(L)
    }
    res <- list(x.values=x[order(idx)], y.values=y[order(idx)], l.value=L,
                g=iso_model, risk=risk, optimizer=optimizer)
  }
  res$call <- match.call(); class(res) <- "iso.lin.est"
  return(res)
}

#' Inherit print function for \code{iso.lin.est}
#'
#' @param obj An object with S3 class \code{iso.lin.est}.
#' @param ... Additional control parameter
#'
#' @export
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- iso.lin.est(x, y)
#' print(res.est)
print.iso.lin.est <- function(obj, ...){
  cat("Call:\n")
  print(obj$call)
  cat("Empirical Risk:\n")
  print(obj$risk)
  cat("L:\n")
  print(obj$l.value)
  cat("Optimizer:\n")
  print(obj$optimizer)
}

#' Inherit predict function for \code{iso.lin.est}
#'
#' @param obj An object with S3 class \code{iso.lin.est}.
#' @param newdata \code{n x 1} numeric vector of new covariate vector. If not
#'  specified, use the covariate vector from training.
#' @return \code{n x 1} numeric vector of predicted values.
#'
#' @export
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- iso.lin.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=runif(n, -1, 1)) # predict new data
predict.iso.lin.est <- function(obj, newdata=NULL){
  # if newdata is not specified, return fitted y
  if(is.null(newdata)){
    newdata <- obj$x
  }
  if (length(obj$l.value) == 1){
    return(obj$g(newdata) + obj$l.value*newdata)
  }
  else{
    # for multiple estimators, prediction is an average over them
    y.hat <- rep(0, length(newdata)); Ls <- length(obj$l.value)
    for (i in 1:Ls){
      y.hat <- y.hat + (obj$g[[i]](newdata) + obj$l.value[[i]]*newdata)/Ls
    }
    return(y.hat)
  }
}
