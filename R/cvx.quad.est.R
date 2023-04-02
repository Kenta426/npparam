# Estimating twice differentiable function via convex regression +
# quadratic function.

#' Nonparametric LSE for twice differentiable function
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector
#' @param y \code{n x 1} numeric vector of observed response vector
#' @param L A scalar for the quadratic parameter. If not specified, L is
#'  data-adaptively selected via optimization or cross-validation
#' @param run.optim A Boolean parameter for selecting optimization or
#'  cross-validation. If TRUE, it performs optimization. Default is TRUE.
#' @param ...
#' @return An object with S3 class \code{cvx.quad.est}.
#' @export
#'
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- cvx.quad.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=runif(n, -1, 1)) # predict new data
cvx.quad.est <- function(x, y, L=NULL, run.optim=TRUE, ...){
  UseMethod("cvx.quad.est")
}

#' Helper function for fitting convex (or concave) function
.fit.cvx.reg <- function(x, y, L){
  # Find convex LSE.
  cvx_final <- simest::cvx.lse.reg(x, y - L*x^2)
  cvx.risk <- mean((y - cvx_final$fit.values - L*x^2)^2)
  # Find concave LSE.
  concave_final <- simest::cvx.lse.reg(x, -(y - L*x^2))
  concave.risk <- mean((y + concave_final$fit.values - L*x^2)^2)
  # pick estimator with smaller empirical risk.
  if(cvx.risk < concave.risk){
    cvx_model <- cvx_final; risk <- cvx.risk; convex <- TRUE
  }else{
    cvx_model <- concave_final; risk <- concave.risk; convex <- FALSE
  }
  return(list(reg=cvx_model, risk=risk, convex=convex))
}

#' Helper function for performing CV for the quadratic parameter
.run.cross.validation.cvx <- function(x, y, L0=0){
  # split data into train/validation
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  # create a list of L values to run CV
  L.list <- seq(L0, (L0+100)*log(n), length.out=n/2)
  sim.one <- function(l){
    reg.tr <- .fit.cvx.reg(x.train, y.train, L=l)
    if(reg.tr$convex){
      y.hat <- predict(reg.tr$reg, newdata=x.val) + l * x.val^2
    }else{
      y.hat <- -predict(reg.tr$reg, newdata=x.val) + l * x.val^2
    }
    mean((y.val-y.hat)^2)
  }
  val.risk <- sapply(L.list, sim.one); idx <- which.min(val.risk)
  # return L with the smallest validation error
  return(list(l.opt = L.list[idx], risk =val.risk[idx]))
}

#' Helper function for running optimizer for the quadratic parameter
.run.optimizer.cvx <- function(x, y, L0=0){
  # split data into train/validation
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  sim.one <- function(l){
    reg.tr <- .fit.cvx.reg(x.train, y.train, L=l)
    if(reg.tr$convex){
      y.hat <- predict(reg.tr$reg, newdata=x.val) + l * x.val^2
    }else{
      y.hat <- -predict(reg.tr$reg, newdata=x.val) + l * x.val^2
    }
    mean((y.val-y.hat)^2)
  }
  # run optimizer to find the best L parameter
  opt.res <- optimise(sim.one, lower = L0, upper = (100+L0)*log(n),
                      maximum = FALSE)
  # return L with the smallest validation error
  return(list(l.opt = opt.res$minimum, risk =opt.res$objective))
}


#' Nonparametric LSE for twice differentiable  function
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector
#' @param y \code{n x 1} numeric vector of observed response vector
#' @param L A scalar for the quadratic parameter. If not specified, L is
#'  data-adaptively selected via optimization or cross-validation
#' @param run.optim A Boolean parameter for selecting optimization or
#'  cross-validation. If TRUE, it performs optimization. Default is TRUE.
#' @param ... Additional control parameters
#' @return An object with S3 class \code{cvx.quad.est}.
#' @export
#'
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- cvx.quad.est(x, y)
cvx.quad.est.default <- function(x, y, L=NULL, run.optim=TRUE, ...){
  idx <- order(x); x <- x[idx]; y <- y[idx]
  # TODO: check input format
  if(run.optim){
    param.selector <- .run.optimizer.cvx; optimizer <- "optimise"
  }else{
    param.selector <- .run.cross.validation.cvx; optimizer <- "CV"
  }
  if(is.null(L)){
    L <- param.selector(x, y)$l.opt
    est.res <- .fit.cvx.reg(x, y, L)
    cvx_model <- est.res$reg; risk <- est.res$risk; convex=est.res$convex
  }else{
    est.res <- .fit.cvx.reg(x, y, L)
    cvx_model <- est.res$reg; risk <- est.res$risk; convex=est.res$convex
  }
  res <- list(x.values = x[order(idx)], y.values = y[order(idx)], l.value=L,
              g=cvx_model, risk=risk,
              convex=convex, optimizer = optimizer)
  res$call <- match.call()
  class(res) <- "cvx.quad.est"
  return(res)
}

#' Inherit print function for \code{cvx.quad.est}
#'
#' @param obj An object with S3 class \code{iso.lin.est}.
#' @param ... Additional control parameter
#' @export
#'
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- cvx.quad.est(x, y)
#' print(res.est)
print.cvx.quad.est <- function(obj, ...){
  cat("Call:\n")
  print(obj$call)
  cat("Empirical Risk:\n")
  print(obj$risk)
  cat("L:\n")
  print(obj$l.value)
  cat("Convex:\n")
  print(obj$convex)
  cat("Optimizer:\n")
  print(obj$optimizer)
}

#' Inherit predict function for \code{cvx.quad.est}
#'
#' @param obj An object with S3 class \code{cvx.quad.est}.
#' @param newdata \code{n x 1} numeric vector of new covariate vector. If not
#'  specified, use the covariate vector from training.
#' @return \code{n x 1} numeric vector of predicted values.
#'
#' @export
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- cvx.quad.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=runif(n, -1, 1)) # predict new data
predict.cvx.quad.est <- function(object, newdata=NULL){
  if(is.null(newdata)){
    newdata <- object$x
  }
  if(object$convex){
    y.hat <- predict(object$g, newdata=newdata) + object$l.value * newdata^2
  }else{
    y.hat <- -predict(object$g, newdata=newdata) + object$l.value * newdata^2
  }
  return(y.hat)
}
