# Estimating additive model where each component is L-Lipscitz
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
additive.iso.est <- function(x, y, L=NULL, run.optim=TRUE, ...){
  UseMethod("additive.iso.est")
}

.predict.iso.dim <- function(gs, x, skip.dim=NULL){
  fn <- gs$g; d <- gs$dim; n <- dim(x)[1]
  if (!is.null(skip.dim) && skip.dim == d){
    return(rep(0, n))
  }
  else{
    return(fn(x[, d]))
  }
}

.predict.additive.iso <- function(glist, x, skip.dim=NULL){
  d <- dim(x)[2]
  ys <- lapply(glist, function(l){.predict.iso.dim(l, x, skip.dim=skip.dim)})
  rowSums(matrix(unlist(ys), ncol = d))
}

.update.linear <- function(x, y, glist, L.max=1, L.min=0){
  d <- dim(x)[2]
  g.hat <- .predict.additive.iso(glist, x)
  est.L <- bvls::bvls(x, y-g.hat, bu = rep(L.max, d), bl = rep(L.min, d))$x
  return(est.L)
}

.update.isotonic <- function(x, y, glist, dim, L){
  g.minus.dim <- .predict.additive.iso(glist, x, skip.dim = dim)
  pseudo.y <- y - x%*%L - g.minus.dim
  gs <- glist[[dim]]
  x.idx <- gs$x.idx; y.idx <- gs$y.idx
  iso.incr <- Iso::pava(pseudo.y[x.idx]) # fit isotonic after sorting Y in the order of X
  iso.fn.incr <- stepfun(x = x[x.idx, dim], y = c(iso.incr[1], iso.incr))
  incr.risk <- mean((y - x%*%L - g.minus.dim - iso.fn.incr(x[,dim]))^2)

  iso.decr <- Iso::pava(pseudo.y[x.idx], decreasing = TRUE) # fit isotonic after sorting Y in the order of X
  iso.fn.decr <- stepfun(x = x[x.idx, dim], y = c(iso.decr[1], iso.decr))
  decr.risk <- mean((y - x%*%L - g.minus.dim - iso.fn.decr(x[,dim]))^2)

  if(decr.risk < incr.risk){
    glist[[dim]]$g <- iso.fn.decr
  }else{
    glist[[dim]]$g <- iso.fn.incr
  }
  glist
}

.update.additive.isotonic <- function(x, y, glist, L){
  for (i in 1:length(glist)) {
    glist <- .update.isotonic(x, y, glist, i, L)
  }
  glist
}
.coordinate.descent <- function(x, y){
  n <- length(y)
  L.min <- 0; L.max <- (L.min + 1) * n; d <- dim(x)[2]
  g.list <- list()
  # identity mapping as base estimators
  # compute sort all at once
  for (i in 1:d){
    g.list[[i]] <- list(g=function(t){t},
                        dim=i, x.idx=order(x[,i]),
                        y.idx=order(order(x[,i])))
  }
  risk.diff <- Inf; risk.0 <- -1; itr <- 1
  while (risk.diff > 1e-05 || itr > 1000){
    L <- .update.linear(x, y, g.list, L.max = L.max, L.min=L.min)
    g.list <- .update.additive.isotonic(x, y, g.list, L)
    y.hat <- .predict.additive.iso(g.list, x) + x%*% L
    risk <- mean((y-y.hat)^2)
    risk.diff <- abs(risk - risk.0)
    risk.0 <- risk; itr <- itr + 1
  }
  converged <- (itr < 1000)
  return(list(L=L, g.list=g.list, risk=risk, converged=converged,
              itr=itr, y.hat=y.hat))
}
#' Title
#'
#' @param x
#' @param y
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
additive.iso.est.default <- function(x, y, ...){
  n <- dim(x)[1]; d <- dim(x)[2]
  opt.res <- .coordinate.descent(x, y)
  gs <- lapply(opt.res$g.list, function(l){list(g=l$g, dim=l$dim)})
  res <- list(gs = gs, L=opt.res$L, x.values=x, y.values=y,
              fitted=opt.res$y.hat, risk=opt.res$risk)
  res$call <- match.call()
  class(res) <- "additive.iso.est"
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
print.additive.iso.est <- function(obj){
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
predict.additive.iso.est <- function(obj, newdata=NULL){
  if(is.null(newdata)){
    newdata <- obj$x.values
  }
  # apply additive functions to each
  y.hat <- .predict.additive.iso(obj$gs, newdata) + newdata %*% obj$L
  return(y.hat)
}

# # example
# n <- 1000
# x <- matrix(runif(n*3, -1,1), ncol=3)
# y <- abs(x[,1]) -x[,2]^2 + rnorm(n, 0, 0.05)
# additive.iso <- additive.iso.est(x, y)
# y.hat <- predict.additive.iso.est(additive.iso)
# plot(y, y.hat)
# library(ggplot2)
# plot.df <- data.frame(x1=x[,1], x2=x[,2], y=y)
# ggplot(plot.df) + geom_point(aes(x=x1, y=x2, fill=y), size=10)
# y <- sin(x*4) + rnorm(n,0,0.05)
#
# res <- iso.lin.est(x, y, run.optim = TRUE)
# res
# y.hat <- predict(res)
# plot(res$x.values, res$y.values, pch=".")
# lines(res$x.values, y.hat, col="red")

