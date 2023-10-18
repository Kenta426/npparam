# Estimating additive model where each component is L-Lipschitz

#' Nonparametric LSE for additive L-Lipschitz functions
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector
#' @param y \code{n x 1} numeric vector of observed response vector
#' @param ...
#' @return An object with S3 class \code{additive.iso.est}.
#' @export
#'
#' @examples
#' n <- 100*2; x <- matrix(runif(n*2, -1, 1), ncol=2);
#' y <- abs(x[,1]) - abs(x[,2]) + rnorm(n, 0, 0.1)
#' res.est <- additive.iso.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=matrix(runif(n*2, -1, 1), ncol=2))
additive.iso.est <- function(x, y, L=NULL, ...){
  UseMethod("additive.iso.est")
}

#' Helper function for predicting isotonic for given dimension.
.predict.iso.dim <- function(gs, x, skip.dim=NULL){
  fn <- gs$g; d <- gs$dim; n <- dim(x)[1]
  if (!is.null(skip.dim) && skip.dim == d){
    return(rep(0, n))
  }
  else{
    return(fn(x[, d]))
  }
}

#' Helper function for predicting isotonic functions for each dimension and
#' adding the predicted values.
.predict.additive.iso <- function(glist, x, skip.dim=NULL){
  d <- dim(x)[2]
  ys <- lapply(glist, function(l){.predict.iso.dim(l, x, skip.dim=skip.dim)})
  rowSums(matrix(unlist(ys), ncol = d))
}

#' Helper function for fitting bounded linear regression.
.update.linear <- function(x, y, glist, L.max=1, L.min=0){
  d <- dim(x)[2]
  g.hat <- .predict.additive.iso(glist, x)
  est.L <- bvls::bvls(x, y-g.hat, bu = rep(L.max, d), bl = rep(L.min, d))$x
  return(est.L)
}

#' Helper function for fitting isotonic function for given dimension.
.update.isotonic <- function(x, y, glist, dim, L){
  # get prediction from all isotonic functions but skip the current dimension
  g.minus.dim <- .predict.additive.iso(glist, x, skip.dim = dim)
  pseudo.y <- y - x%*%L - g.minus.dim
  # sort across given dimension
  gs <- glist[[dim]]
  x.idx <- gs$x.idx; y.idx <- gs$y.idx
  # fit non-decreasing function after sorting Y in the order of X
  iso.incr <- Iso::pava(pseudo.y[x.idx])
  iso.fn.incr <- stepfun(x = x[x.idx, dim], y = c(iso.incr[1], iso.incr))
  incr.risk <- mean((pseudo.y - iso.fn.incr(x[,dim]))^2)
  # fit non-increasing function after sorting Y in the order of X
  iso.decr <- Iso::pava(pseudo.y[x.idx], decreasing = TRUE)
  iso.fn.decr <- stepfun(x = x[x.idx, dim], y = c(iso.decr[1], iso.decr))
  decr.risk <- mean((pseudo.y - iso.fn.decr(x[,dim]))^2)
  # return functions with smaller risk
  if(decr.risk < incr.risk){
    glist[[dim]]$g <- iso.fn.decr
  }else{
    glist[[dim]]$g <- iso.fn.incr
  }
  return(glist)
}

#' Helper function for fitting isotonic function across dimension.
.update.additive.isotonic <- function(x, y, glist, L){
  for (i in 1:length(glist)) {
    glist <- .update.isotonic(x, y, glist, i, L)
  }
  return(glist)
}

#' Helper function for performing coordinate descent
.coordinate.descent <- function(x, y, L){
  n <- dim(x)[1]; d <- dim(x)[2]
  # prepare the initial functions
  g.list <- list()
  for (i in 1:d){
    g.list[[i]] <- list(g=function(t){sapply(t, function(u){0})},
                        dim=i, x.idx=order(x[,i]),
                        y.idx=order(order(x[,i])))
  }
  risk.diff <- Inf; risk.0 <- -1; itr.count <- 1; converged <- FALSE
  # main loop for coordinate descent
  # run loops until either (1) improvement of the risk between step is small or
  # (2) maximum iteration is reached.
  while (risk.diff > 1e-5 && itr.count < 1000){
    # update linear parameter with current list of isotonic functions in g.list
    # L <- .update.linear(x.train, y.train, g.list, L.max = L.max, L.min=L.min)
    # update list of isotonic functions with current L
    g.list <- .update.additive.isotonic(x, y, g.list, L)
    y.hat <- .predict.additive.iso(g.list, x) + x%*% L
    risk <- mean((y-y.hat)^2)
    # check if the improvement is small
    risk.diff <- abs(risk - risk.0)
    risk.0 <- risk; itr.count <- itr.count + 1
  }
  converged <- (itr.count < 1000)
  return(list(L=L, g.list=g.list, risk=risk, converged=converged,
              itr=itr.count, y.hat=y.hat))
}


.run.optimizer.additive <- function(x, y, L0=NULL){
  n <- dim(x)[1]; d <- dim(x)[2]
  L.min <- -10*log(n); L.max <- 10*log(n); d <- dim(x)[2]
  # split data into train/validation
  learn.n <- as.integer(n-sqrt(n))
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id, ]; x.val <- x[-train.id, ]
  y.train <- y[train.id]; y.val <- y[-train.id]
  sim.one <- function(l){
    reg.tr <- .coordinate.descent(x.train, y.train, L=l)
    y.hat <- .predict.additive.iso(reg.tr$g.list, x.val) + x.val%*% l
    mean((y.val-y.hat)^2)
  }
  if (is.null(L0)){
    opt.res <- optim(rep(0, d), sim.one, lower = L.min, upper = L.max, method = "L-BFGS-B")
    L0 <-opt.res$par
  }

  reg.tr <- .coordinate.descent(x.train, y.train, L=L0)
  y.hat <- .predict.additive.iso(reg.tr$g.list, x.val) + x.val%*% L0
  risk <- mean((y.val-y.hat)^2)
  return(list(L=L0, g.list=reg.tr$g.list,
              risk=risk, converged=opt.res$convergence,
              itr=0, y.hat=y.hat))
}
#' Nonparametric LSE for for additive L-Lipschitz functions
#'
#' @param x \code{n x 1} numeric vector of observed covariate vector
#' @param y \code{n x 1} numeric vector of observed response vector
#' @param ...
#' @return An object with S3 class \code{additive.iso.est}.
#' @export
#'
#' @examples
#' n <- 100*2; x <- matrix(runif(n*2, -1, 1), ncol=2);
#' y <- abs(x[,1]) - abs(x[,2]) + rnorm(n, 0, 0.1)
#' res.est <- additive.iso.est(x, y)
additive.iso.est.default <- function(x, y, L=NULL,...){

  # opt.res <- .coordinate.descent(x, y)
  opt.res <- .run.optimizer.additive(x, y, L=L)
  gs <- lapply(opt.res$g.list, function(l){list(g=l$g, dim=l$dim)})
  res <- list(gs = gs, L=opt.res$L, x.values=x, y.values=y,
              itr = opt.res$itr, converged=opt.res$converged,
              fitted=opt.res$y.hat, risk=opt.res$risk)
  res$call <- match.call()
  class(res) <- "additive.iso.est"
  return(res)
}

#' Inherit print function for \code{additive.iso.est}
#'
#' @param obj An object with S3 class \code{additive.iso.est}.
#' @param ... Additional control parameter
#' @export
#' @examples
#' n <- 100*2; x <- matrix(runif(n*2, -1, 1), ncol=2);
#' y <- abs(x[,1]) - abs(x[,2]) + rnorm(n, 0, 0.1)
#' res.est <- additive.iso.est(x, y)
#' print(res.est)
print.additive.iso.est <- function(obj, ...){
  cat("Call:\n")
  print(obj$call)
  cat("Empirical Risk:\n")
  print(obj$risk)
  cat("L:\n")
  print(obj$L)
  cat("# of iteration:\n")
  print(obj$itr)
  cat("Convergence status:\n")
  print(obj$converged)
}

#' Inherit predict function for \code{cvx.quad.est}
#'
#' @param obj n object with S3 class \code{cvx.quad.est}.
#' @param newdata \code{n x d} numeric vector of new covariate matrix If not
#'  specified, use the covariate vector from training.
#' @return \code{n x 1} numeric vector of predicted values.
#' @export
#' @examples
#' n <- 100*2; x <- matrix(runif(n*2, -1, 1), ncol=2);
#' y <- abs(x[,1]) - abs(x[,2]) + rnorm(n, 0, 0.1)
#' res.est <- additive.iso.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=matrix(runif(n*2, -1, 1), ncol=2))
predict.additive.iso.est <- function(obj, newdata=NULL){
  if(is.null(newdata)){
    newdata <- obj$x.values
  }
  # apply additive functions to each
  y.hat <- .predict.additive.iso(obj$gs, newdata) + newdata %*% obj$L
  return(y.hat)
}


# n <- 3000; data <- additive.debug(n, sigma=0)
# res.est <- additive.iso.est(data$x, data$y)
# x <- data$x; y <- data$y
# res.est
# plot(x, res.est$gs[[1]]$g(x)+res.est$L[1]*x)
# plot(x, res.est$gs[[2]]$g(x)+res.est$L[2]*x)
#
# test.x <- matrix(runif(1e4*2, 0, 1), ncol=2)
# idx <- order(test.x)
# test.x[idx]
# test.y <- data$fn(test.x)
#
# log10(mean((test.y-predict(res.est, test.x))^2))


# n <- 1000; data <- additive.debug(n, sigma=0)
# res.est <- additive.iso.est(data$x, data$y)
# x <- data$x; y <- data$y
# idx <- order(x[,1])
# res.est
# plot(x[,1], res.est$gs[[1]]$g(x[,1])+res.est$L[1]*x[,1])
# lines(x[idx,1], lipschitz.fn(1)$fn(x[idx,1]))
# idx <- order(x[,2])
# plot(x[,2], res.est$gs[[2]]$g(x[,2])+res.est$L[2]*x[,2])
# lines(x[idx,2], -lipschitz.fn(1)$fn(x[idx,2]))
#
# mean((test.y-predict(res.est, test.x))^2)
# plot(test.y, predict(res.est, test.x))

# plot(x, res.est$gs[[1]]$g(x)+res.est$L[1]*x)
# plot(x, res.est$gs[[2]]$g(x)+res.est$L[2]*x)
# plot(x, -lipschitz.fn(1)$fn(x))
