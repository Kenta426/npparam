# Estimating L-Lipschitz function via isotonic regression + linear function.

#' Nonparametric LSE for Lipschitz function
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector
#' @param y \code{n x 1} numeric vector of observed response vector
#' @param L A scalar for the linear parameter. If not specified, L is
#'  data-adaptively selected via optimization or cross-validation
#' @param run.optim A Boolean parameter for selecting optimization or
#'  cross-validation. If TRUE, it performs optimization. Default is TRUE.
#' @param ... Additional control parameters
#' @return An object with S3 class \code{iso.lin.est}.
#'
#' @export
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- iso.lin.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=runif(n, -1, 1)) # predict new data
iso.lin.est <- function(x, y, L=NULL, run.optim=TRUE, ...){
  UseMethod("iso.lin.est")
}

#' Helper function for fitting isotonic function
.fit.iso.reg <- function(x, y, L){
  # Use the min-max formula to find non-decreasing LSE.
  iso_final.incr <- Iso::pava(y - L*x)
  iso_model.incr <- stepfun(x=x, y=c(iso_final.incr[1], iso_final.incr))
  incr.risk <- mean((y - iso_model.incr(x) - L*x)^2)
  # Use the min-max formula to find non-increasing LSE.
  iso_final.decr <- Iso::pava(y - L*x, decreasing = TRUE)
  iso_model.decr <- stepfun(x = x, y = c(iso_final.decr[1], iso_final.decr))
  decr.risk <- mean((y - iso_model.decr(x) - L*x)^2)
  # pick estimator with smaller empirical risk.
  if(decr.risk < incr.risk){
    iso_model <- iso_model.decr; risk <- decr.risk
  }
  else{
    iso_model <- iso_model.incr; risk <- incr.risk
  }
  return(list(reg=iso_model, risk=risk))
}

#' Helper function for performing CV for the linear parameter
.run.cross.validation.iso <- function(x, y, L0=0){
  # split data into train/validation
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size=learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  # create a list of L values to run CV
  L.list <- seq(-log10(n), log10(n), length.out=as.integer(n/20))
  sim.one <- function(l){
    reg.tr <- .fit.iso.reg(x.train, y.train, L=l)
    y.hat <- reg.tr$reg(x.val) + l*x.val
    mean((y.val-y.hat)^2)
  }
  val.risk <- sapply(L.list, sim.one); idx <- which.min(val.risk)
  # return L with the smallest validation error
  return(list(l.opt=L.list[idx], risk=val.risk[idx]))
}

#' Helper function for running optimizer for the linear parameter
.run.optimizer.iso <- function(x, y, L0=0){
  # split data into train/validation
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size=learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  sim.one <- function(l){
    reg.tr <- .fit.iso.reg(x.train, y.train, L=l)
    y.hat <- reg.tr$reg(x.val) + l*x.val
    mean((y.val-y.hat)^2)
  }
  # run optimizer to find the best L parameter
  opt.res <- optimise(sim.one, lower=-log10(n), upper=log10(n),
                      maximum = FALSE)
  # return L with the smallest validation error
  return(list(l.opt=opt.res$minimum, risk=opt.res$objective))
}

#' Nonparametric LSE for Lipschitz function
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector
#' @param y \code{n x 1} numeric vector of observed response vector
#' @param L A scalar for the linear parameter. If not specified, L is
#'  data-adaptively selected via optimization or cross-validation
#' @param run.optim A Boolean parameter for selecting optimization or
#'  cross-validation. If TRUE, it performs optimization. Default is TRUE.
#' @param ... Additional control parameters
#' @return An object with S3 class \code{iso.lin.est}.
#' @export
#' @examples
#' n <- 100
#' x <- runif(n, -1, 1)
#' y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- iso.lin.est(x, y)
iso.lin.est.default <- function(x, y, L=NULL, run.optim=TRUE, ...){
  # TODO: check input format
  idx <- order(x); x <- x[idx]; y <- y[idx]
  # select parameter selection procedure
  if(run.optim){
    param.selector <- .run.optimizer.iso; optimizer <- "optimise"
  }else{
    param.selector <- .run.cross.validation.iso; optimizer <- "CV"
  }
  # run parameter selection if L is not specified
  if(is.null(L)){
    L <- param.selector(x, y)$l.opt
  }
  # fit isotonic regression
  est.res <- .fit.iso.reg(x, y, L)
  iso_model <- est.res$reg; risk <- est.res$risk
  res <- list(x.values=x[order(idx)], y.values=y[order(idx)], l.value=L,
              g=iso_model, risk=risk, optimizer=optimizer)
  res$call <- match.call()
  class(res) <- "iso.lin.est"
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
  return(obj$g(newdata) + obj$l.value*newdata)
}
