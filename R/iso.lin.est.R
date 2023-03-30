# Estimating L-Lipschitz function as isotonic regression + linear function

iso.lin.est <- function(x, y, L=NULL, run.optim=TRUE, ...){
  UseMethod("iso.lin.est")
}

.fit.iso.reg <- function(x, y, L){
  iso_final.incr <- Iso::pava(y - L*x)
  iso_model.incr <- stepfun(x = x, y = c(iso_final.incr[1], iso_final.incr))
  incr.risk <- mean((y - iso_model.incr(x) - L*x)^2)
  iso_final.decr <- Iso::pava(y - L*x, decreasing = TRUE)
  iso_model.decr <- stepfun(x = x, y = c(iso_final.decr[1], iso_final.decr))
  decr.risk <- mean((y - iso_model.decr(x) - L*x)^2)
  if(decr.risk < incr.risk){
    iso_model <- iso_model.decr; risk <- decr.risk
  }
  else{
    iso_model <- iso_model.incr; risk <- incr.risk
  }
  res <- list(reg=iso_model, risk=risk)
  return(res)
}

.run.cross.validation.iso <- function(x, y, val=1/4, L0=0){
  n <- length(x); learn.n <- n*(1-val)
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  L.list <- seq(L0, (L0+1)*log(n), length.out=n/2)
  sim.one <- function(l){
    reg.tr <- .fit.iso.reg(x.train, y.train, L=l)
    y.hat <- reg.tr$reg(x.val) + l*x.val
    mean((y.val-y.hat)^2)
  }
  val.risk <- sapply(L.list, sim.one); idx <- which.min(val.risk)
  list(l.opt = L.list[idx], risk =val.risk[idx])
}

.run.optimizer.iso <- function(x, y, val=1/4, L0=0){
  n <- length(x); learn.n <- n*(1-val)
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id]; x.val <- x[-train.id]
  y.train <- y[train.id]; y.val <- y[-train.id]
  sim.one <- function(l){
    reg.tr <- .fit.iso.reg(x.train, y.train, L=l)
    y.hat <- reg.tr$reg(x.val) + l*x.val
    mean((y.val-y.hat)^2)
  }
  opt.res <- optimise(sim.one, lower = L0, upper = (1+L0)*log(n),
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
iso.lin.est.default <- function(x, y, L=NULL, run.optim=TRUE, ...){
  idx <- order(x)
  x <- x[idx]; y <- y[idx]
  # TODO: check input format
  if(run.optim){
    param.selector <- .run.optimizer.iso; optimizer <- "optimise"
  }
  else{
    param.selector <- .run.cross.validation.iso; optimizer <- "CV"
  }
  if(is.null(L)){
    L <- param.selector(x, y)$l.opt
    est.res <- .fit.iso.reg(x, y, L)
    iso_model <- est.res$reg; risk <- est.res$risk
  }
  else{
    est.res <- .fit.iso.reg(x, y, L)
    iso_model <- est.res$reg; risk <- est.res$risk
  }
  res <- list(x.values = x, y.values = y, l.value=L, g=iso_model, risk=risk,
              optimizer = optimizer)
  res$call <- match.call()

  class(res) <- "iso.lin.est"
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
print.iso.lin.est <- function(obj){
  cat("Call:\n")
  print(obj$call)
  cat("Empirical Risk:\n")
  print(obj$risk)
  cat("L:\n")
  print(obj$l.value)
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
predict.iso.lin.est <- function(obj, newdata=NULL){
  if(is.null(newdata)){
    newdata <- obj$x
  }
  idx <- order(newdata); newdata <- newdata[idx]
  y.hat <- obj$g(newdata) + obj$l.value * newdata
  return(y.hat[order(idx)])
}

# # example
# n <- 5000
# x <- runif(n, -1,1)
# y <- abs(x) + rnorm(n,0,0.05)
# y <- sin(x*4) + rnorm(n,0,0.05)
#
# res <- iso.lin.est(x, y, run.optim = TRUE)
# res
# y.hat <- predict(res)
# plot(res$x.values, res$y.values, pch=".")
# lines(res$x.values, y.hat, col="red")

