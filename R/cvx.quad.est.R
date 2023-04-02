#' Title
#'
#' @param x
#' @param y
#' @param L
#' @param run.optim
#' @param ...
#'
#' @return
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
  cvx_final <- simest::cvx.lse.reg(x, y - L*x^2)
  cvx.risk <- mean((y - cvx_final$fit.values - L*x^2)^2)
  concave_final <- simest::cvx.lse.reg(x, -(y - L*x^2))
  concave.risk <- mean((y + concave_final$fit.values - L*x^2)^2)
  if(cvx.risk < concave.risk){
    cvx_model <- cvx_final; risk <- cvx.risk; convex <- TRUE
  }
  else{
    cvx_model <- concave_final; risk <- concave.risk; convex <- FALSE
  }
  res <- list(reg=cvx_model, risk=risk, convex=convex)
  return(res)
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
    }
    else{
      y.hat <- -predict(reg.tr$reg, newdata=x.val) + l * x.val^2
    }
    mean((y.val-y.hat)^2)
  }
  opt.res <- optimise(sim.one, lower = L0, upper = (100+L0)*log(n),
                      maximum = FALSE)
  list(l.opt = opt.res$minimum, risk =opt.res$objective)
}

#' Helper function for running optimizer for the quadratic parameter
.run.optimizer.cvx <- function(x, y, L0=0){
  n <- length(x); learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  sim.one <- function(l){
    reg.tr <- .fit.cvx.reg(x.train, y.train, L=l)
    if(reg.tr$convex){
      y.hat <- predict(reg.tr$reg, newdata=x.val) + l * x.val^2
    }
    else{
      y.hat <- -predict(reg.tr$reg, newdata=x.val) + l * x.val^2
    }
    mean((y.val-y.hat)^2)
  }
  opt.res <- optimise(sim.one, lower = L0, upper = (100+L0)*log(n),
                      maximum = FALSE)
  list(l.opt = opt.res$minimum, risk =opt.res$objective)
}


#' Title
#'
#' @param x
#' @param y
#' @param L
#' @param run.optim
#' @param ...
#'
#' @return
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
  }
  else{
    param.selector <- .run.cross.validation.cvx; optimizer <- "CV"
  }
  if(is.null(L)){
    L <- param.selector(x, y)$l.opt
    est.res <- .fit.cvx.reg(x, y, L)
    cvx_model <- est.res$reg; risk <- est.res$risk; convex=est.res$convex
  }
  else{
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

#' Title
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- cvx.quad.est(x, y)
#' print(res.est)
print.cvx.quad.est <- function(obj){
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

#' Title
#'
#' @param obj
#' @param newdata
#'
#' @return
#' @export
#'
#' @examples
#' n <- 100; x <- runif(n, -1, 1); y <- sin(x*4) + rnorm(n, 0, 0.1)
#' res.est <- cvx.quad.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=runif(n, -1, 1)) # predict new data
predict.cvx.quad.est <- function(object, newdata=NULL){
  if(is.null(newdata)){
    newdata <- object$x
  }
  idx <- order(newdata); newdata <- newdata[idx]
  if(object$convex){
    y.hat <- predict(object$g, newdata=newdata) + object$l.value * newdata^2
  }
  else{
    y.hat <- -predict(object$g, newdata=newdata) + object$l.value * newdata^2
  }
  return(y.hat[order(idx)])
}
