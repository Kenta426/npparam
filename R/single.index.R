# Estimating single-index model where the link function is L-lipschitz

#' Nonparametric LSE for single-index model with a Lipschitz link function
#'
#' @param x \code{n x d} numeric vector of the observed covariate matrix
#' @param y \code{n x 1} numeric vector of the observed response vector.
#' @param L0 L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 \times log10(n), L0 \times log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
#' @param ... Additional control parameters.
#'
#' @return An object with S3 class \code{single.index.est}.
#' @export
#'
#' @examples
#' n <- 100; d <- 3; x <- matrix(runif(n*d, -1, 1), nrow=n);
#' beta <- c(1,0,0); y <- sin(4*x%*%beta) + rnorm(n, 0, 0.1)
#' res.est <- single.index.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=matrix(runif(n*d, -1, 1))) # predict new data
single.index.est <- function(x, y, L0=100, ...){
  UseMethod("single.index.est")
}

#' Helper function for running an iterative optimizer
#'
#' @param x \code{n x d} numeric vector of the observed covariate matrix
#' @param y \code{n x 1} numeric vector of the observed response vector.
#' @param L0 L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 \times log10(n), L0 \times log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
.run.optimizer.single.index <- function(x, y, L0=100){
  n <- dim(x)[1]; d <- dim(x)[2]
  learn.n <- as.integer(n-n/2)
  train.id <- sort(sample(seq_len(n), size =learn.n))
  x.train <- x[train.id, ]; x.val <- x[-train.id, ]
  y.train <- y[train.id]; y.val <- y[-train.id]
  # initialize beta and link function
  beta <- runif(d); beta <- beta/norm(beta, "2")
  reg.est <- iso.lin.est(x.train%*%beta, y.train,
                     run.cross.fit=FALSE, L0=L0) # fit lipschitz link
  val.error <- Inf; val.diff <- Inf; itr.step <- 0
  while(val.diff > 1e-5 & itr.step < 30){
    sim.one <- function(beta.param){
      beta.param <- beta.param/norm(beta.param, "2")
      mean((y.val-predict(reg.est, newdata=x.val%*%beta.param))^2)
    }
    # optimize beta for a fixed link function
    res <- optim(c(beta), sim.one,method="L-BFGS-B", upper = 1, lower = -1)
    beta <- res$par/norm(res$par, "2")
    # optimize link function for a fixed beta parameter
    reg.est <- iso.lin.est(x.train%*%beta, y.train,
                          run.cross.fit=FALSE, L0=L0)
    # update stopping criteria
    val.error <- abs(val.error-res$value)
    val.error <- res$value; itr.step <- itr.step + 1
  }
  return(list(beta.opt=beta, risk=val.error))
}

#' Nonparametric LSE for for sing-index L-lipschitz functions
#'
#' @param x \code{n x d} numeric vector of the observed covariate matrix
#' @param y \code{n x 1} numeric vector of the observed response vector.
#' @param L0 L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 \times log10(n), L0 \times log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
#' @param ...
single.index.est.default <- function(x, y, L0=100, ...){
  opt.res <- .run.optimizer.single.index(x, y, L0=L0)
  beta <- opt.res$beta.opt
  est.res <- iso.lin.est(x%*%beta, y, run.cross.fit=F, L0=L0)
  x.vec <- x%*%beta; idx <- order(x.vec); x.vec <- x.vec[idx]; y <- y[idx]
  res <- list(x.values=x, y.values=y,
              g=est.res, beta=beta, risk=risk)
  res$call <- match.call()
  class(res) <- "single.index.est"
  return(res)
}

#' Inherit print function for \code{single.index.est}
#'
#' @param obj An object with S3 class \code{single.index.est}.
#' @param ... Additional control parameter
#'
#' @export
#' @examples
#' n <- 100; d <- 3; x <- matrix(runif(n*d, -1, 1), nrow=n);
#' beta <- c(1,0,0); y <- sin(4*x%*%beta) + rnorm(n, 0, 0.1)
#' res.est <- single.index.est(x, y)
#' print(res.est)
print.single.index.est <- function(obj, ...){
  cat("Call:\n")
  print(obj$call)
  cat("Empirical Risk:\n")
  print(obj$risk)
  cat("beta:\n")
  print(obj$beta)
}

#' Inherit predict function for \code{single.index.est}
#'
#' @param obj n object with S3 class \code{single.index.est}.
#' @param newdata \code{n x d} numeric vector of new covariate matrix If not
#'  specified, use the covariate vector from training.
#' @return \code{n x 1} numeric vector of predicted values.
#' @export
#' @examples
#' n <- 100; d <- 3; x <- matrix(runif(n*d, -1, 1), nrow=n);
#' beta <- c(1,0,0); y <- sin(4*x%*%beta) + rnorm(n, 0, 0.1)
#' res.est <- single.index.est(x, y)
#' y.pred <- predict(res.est) # fitted value
#' y.pred <- predict(res.est, newdata=matrix(runif(n*d, -1, 1))) # predict new data
predict.single.index.est <- function(obj, newdata=NULL){
  if(is.null(newdata)){
    newdata <- obj$x.values
  }
  newdata <- newdata %*% obj$beta
  y.hat <- predict.iso.lin.est(obj$g, newdata)
  return(y.hat)
}
