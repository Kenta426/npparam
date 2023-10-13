# Estimating L-Lipschitz function over 2d lattice.

#' Nonparametric LSE for L-Lipschitz function in L2 norm over the 2d lattice.
#' The implementation is based on the fixed design setting and currently does
#' not support random designs. The core implementation is based on block max-min
#' algorithm proposed by Deng and Zhang (2020), which runs O(n^3). Since
#' we have to run block max-min algorithm inside of the numerical optimizer,
#' this algorithm can take long time to complete.
#'
#' @param x \code{n1 x n2} numeric matrix of observed covariates. Each element
#' of this matrix is a list of vector. For instance, the (i,j)th element of the
#' matrix corresponds to x[i, j] = list(c(a, b)), which takes the value (a, b).
#' @param y \code{n1 x n2} numeric matrix of observed responses.
#' @param L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 x log10(n), L0 x log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
#' @param run.itr A Boolean parameter for selecting iterative optimization or
#' joint optimization. The iterative optimization gives an approximated minimum
#' but it tends to be faster. Default is TRUE.
#' @param ... Additional control parameters
#' @return An object with S3 class \code{block.iso.lin.est}.
#' @export
#' @examples
#' # creates sample
#' n1 <- 10; n2 <- 10
#' y <- matrix(rep(0, n1*n2), nrow = n1)
#' x <- matrix(rep(list(), n1*n2),nrow = n1, ncol =n2)
#' for (i in 1:n1){
#'   for (j in 1:n2){
#'     y[i,j] <- sin(10*norm(c(i/n1, j/n2), "2")) + rnorm(1, 0, 0.1)
#'     x[i,j] <- list(c(i/n1, j/n2))
#'   }
#' }
#' res.est <- block.iso.lin.est(x, y)
#' y.hat <- predict(res.est) # estimated values over the grid.
block.iso.lin.est <- function(x, y, L0=100, run.itr=TRUE,...){
  UseMethod("block.iso.lin.est")
}

#' Helper function for running block max min algorithm.
#' Initially implemented by Hang Deng who wrote the original block max min
#' paper (Hang and Zhang, 2022)
#'
#' @param y \code{n1 x n2} numeric matrix of observed responses
.fit.block.max.min <- function(y){
  # y_{ij} = f(x_{ij}) + eps_{ij}, where {x_{ij}} form a lattice
  ## 1. pre-processing
  n <- prod(dim(y)) # total number of points
  V <- which(array(rep(1,n), dim(y))==1, arr.ind = T) # vector of vertices
  A <- t(sapply(1:n, FUN = function(k){
    return(sapply(1:n, FUN =function(l){
      return(prod(V[k,] <= V[l,]))}))})) # adjacency matrix
  n1 <- dim(y)[1]
  E.coor <- which(A==1, arr.ind = T)
  # vector of edges, represented by (col1, col2) -> (col3, col4)
  E <- cbind((E.coor - 1) %% n1 + 1, (E.coor - 1) %/% n1 + 1)
  # vector of edges, represented by the indices of nonzero adjacency matrix entries
  E_index <- which(A==1)

  ## 2. calculate block average matrix mu, mu_[u,v] = Average of y_{ij} for {x_{ij} in [u,v]}
  mu = matrix(nrow = n, ncol = n)
  mu[E_index] = apply(E, MARGIN = 1, FUN = function(vec) {
    return(mean(y[vec[1] : vec[2], vec[3]:vec[4]]))
  })

  ## 3. calculate the block max-min estimator
  output = lapply(1:n, FUN = function(j){
    u.ind = which(A[,j]==1)
    v.ind = which(A[j,]==1)
    min_output = sapply(u.ind, function(x){
      return(min(mu[x, v.ind]))
    })
    return(max(min_output))
  })

  ## 4. return value
  return(array(unlist(output), dim = dim(y)))
}

#' Split two matrices with the same dimension into training/validation
#' in the shape of a "checker board"
#'
#' @param x \code{n1 x n2} numeric matrix of observed covariates. Each element
#' of this matrix is a list of vector. For instance, the (i,j)th element of the
#' matrix corresponds to x[i, j] = list(c(a, b)), which takes the value (a, b).
#' @param y \code{n1 x n2} numeric matrix of observed responses.
.split.matrix <- function(x, y){
  n <- prod(dim(y));   ncol <- dim(y)[2]
  V <- which(array(rep(1,n), dim(y))==1, arr.ind = T)
  train.set <- function(v){
    (v[1] %% 2) == (v[2] %% 2)
  }
  val.set <- function(v){
    (v[1] %% 2) != (v[2] %% 2)
  }
  V.train <- V[apply(V, 1, train.set),]
  V.val <- V[apply(V, 1, val.set),]

  x.train <- matrix(x[V.train], ncol=ncol)
  y.train <- matrix(y[V.train], ncol=ncol)
  x.val <- matrix(x[V.val], ncol=ncol)
  y.val <- matrix(y[V.val], ncol=ncol)
  list(x.train=x.train, y.train=y.train, x.val=x.val, y.val=y.val)
}

#' Helper function to perform dot product over list of items.
#'
#' @param beta 2d vector corresponding to the parametric components
.multiply.beta <- function(beta){
  function(M){unlist(M) %*% beta}
}

#' Helper function for running an optimizer for the linear parameter
#'
#' @param x \code{n1 x n2} numeric matrix of observed covariates. Each element
#' of this matrix is a list of vector. For instance, the (i,j)th element of the
#' matrix corresponds to x[i, j] = list(c(a, b)), which takes the value (a, b).
#' @param y \code{n1 x n2} numeric matrix of observed responses.
#' @param L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 x log10(n), L0 x log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
.run.optimizer.block.iso <- function(x, y, L0=100){
  # split matrix into a "checker board"
  data <- .split.matrix(x, y)
  x.train <- data$x.train; y.train <- data$y.train
  x.val <- data$x.val; y.val <- data$y.val; n <- prod(dim(x.train))
  # run optimization to minimize validation error
  sim.one <- function(beta){
    list.prod <- .multiply.beta(beta)
    y.hat <- y.train+apply(x.train, c(1, 2), list.prod)
    y.est <- .fit.block.max.min(y.hat)
    y.hat <- y.est - apply(x.val, c(1, 2), list.prod)
    mean((y.hat - y.val)^2)
  }
  opt.res <- optim(c(1,1), sim.one, method="L-BFGS-B",
                   lower = -L0*log10(n), upper = L0*log10(n))
  return(list(l.opt=opt.res$par, risk=opt.res$value))
}

#' Helper function for running an iterative optimizer for the linear parameter
#'
#' @param x \code{n1 x n2} numeric matrix of observed covariates. Each element
#' of this matrix is a list of vector. For instance, the (i,j)th element of the
#' matrix corresponds to x[i, j] = list(c(a, b)), which takes the value (a, b).
#' @param y \code{n1 x n2} numeric matrix of observed responses.
#' @param L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 x log10(n), L0 x log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
.run.iterative.optimizer.block.iso <- function(x, y, L0=100){
  # split matrix into a "checker board"
  data <- .split.matrix(x, y)
  x.train <- data$x.train; y.train <- data$y.train
  x.val <- data$x.val; y.val <- data$y.val
  n <- prod(dim(x.train)); beta0 <- c(0,0)
  list.prod <- .multiply.beta(beta0)
  # initialize estimator
  y.hat <- y.train+apply(x.train, c(1, 2), list.prod)
  y.est <- .fit.block.max.min(y.hat)
  # initialize stopping criteria
  val.error <- Inf; val.diff <- Inf; itr.step <- 0
  while(val.diff > 1e-5 & itr.step < 100){
    sim.one <- function(beta){
      list.prod <- .multiply.beta(beta)
      y.hat <- y.est - apply(x.val, c(1, 2), list.prod)
      mean((y.hat - y.val)^2)
    }
    # update parameter for fixed function
    opt.res <- optim(beta0, sim.one, method="L-BFGS-B",
                     lower = -L0*log10(n), upper = L0*log10(n))
    beta0 <- opt.res$par
    # update the nonparametric component
    list.prod <- .multiply.beta(beta0)
    y.hat <- y.train+apply(x.train, c(1, 2), list.prod)
    y.est <- .fit.block.max.min(y.hat)
    # update stopping criteria
    val.error <- abs(val.error-opt.res$value)
    val.error <- opt.res$value; itr.step <- itr.step + 1
  }
  return(list(l.opt=opt.res$par, risk=opt.res$value))
}

#' Estimating L-Lipschitz function over 2d lattice.
#'
#' @param x \code{n1 x n2} numeric matrix of observed covariates. Each element
#' of this matrix is a list of vector. For instance, the (i,j)th element of the
#' matrix corresponds to x[i, j] = list(c(a, b)), which takes the value (a, b).
#' @param y \code{n1 x n2} numeric matrix of observed responses.
#' @param L0 The upper bound of the interval of L parameters to search over
#'  the interval is in the form of [-L0 x log10(n), L0 x log10(n)].
#'  In practice, this should be a large enough constant so the set contains
#'  the true Lipschitz parameter L.
#' @param run.itr A Boolean parameter for selecting iterative optimization or
#'  joint optimization. The iterative optimization gives an approximated minimum
#'  but it tends to be faster. Default is TRUE.
#' @param ... Additional control parameters
#'
#' @return An object with S3 class \code{block.iso.lin.est}.
#' @export
#' @examples
#' # creates sample
#' n1 <- 10; n2 <- 10
#' y <- matrix(rep(0, n1*n2), nrow = n1)
#' x <- matrix(rep(list(), n1*n2),nrow = n1, ncol =n2)
#' for (i in 1:n1){
#'   for (j in 1:n2){
#'     y[i,j] <- sin(10*norm(c(i/n1, j/n2), "2")) + rnorm(1, 0, 0.1)
#'     x[i,j] <- list(c(i/n1, j/n2))
#'   }
#' }
#' res.est <- block.iso.lin.est(x, y)
#' y.hat <- predict(res.est) # estimated values over the grid.
block.iso.lin.est.default <- function(x, y, L0=100, run.itr=TRUE,...){
  if (!all(dim(x) == dim(y))){
    stop("Dimension of X and Y do not match")
  }
  if(dim(x)[1] %% 2 == 1){
    warning("The rows of X and Y are odd and the last rows will be removed")
    x <- x[-1, ]; y <- y[-1, ];
  }
  if(dim(x)[2] %% 2 == 1){
    warning("The cols of X and Y are odd and the last cols will be removed")
    x <- x[,-1]; y <- y[, -1];
  }

  # select parameter selection procedure
  if(run.itr){
    param.selector <- .run.iterative.optimizer.block.iso; optimizer <- "Iterative"
  }else{
    param.selector <- .run.optimizer.block.iso; optimizer <- "Optimizer"
  }
  opt.res <- param.selector(x, y, L0=L0)
  L <- opt.res$l.opt; risk <- opt.res$risk
  list.prod <- .multiply.beta(L)
  # get the residuals
  y.hat <- y + apply(x, c(1, 2), list.prod)
  # run the block max-min algorithm
  y.est <- .fit.block.max.min(y.hat)
  # subtract the parametric part
  y.hat <- y.est - apply(x, c(1, 2), list.prod)
  res <- list(x.values=x, y.values=y,
              y.hat=y.hat, l.values=L, risk=risk, optimizer=optimizer)
  res$call <- match.call()
  class(res) <- "block.iso.lin.est"
  return(res)
}

#' Inherit print function for \code{block.iso.lin.est}
#'
#' @param obj An object with S3 class \code{block.iso.lin.est}.
#' @param ... Additional control parameter
#'
#' @export
#' @examples
#' n1 <- 10; n2 <- 10
#' y <- matrix(rep(0, n1*n2), nrow = n1)
#' x <- matrix(rep(list(), n1*n2),nrow = n1, ncol =n2)
#' for (i in 1:n1){
#'   for (j in 1:n2){
#'     y[i,j] <- sin(10*norm(c(i/n1, j/n2), "2")) + rnorm(1, 0, 0.1)
#'     x[i,j] <- list(c(i/n1, j/n2))
#'   }
#' }
#' res.est <- block.iso.lin.est(x, y)
#' print(res.est)
print.block.iso.lin.est <- function(obj, ...){
  cat("Call:\n")
  print(obj$call)
  cat("Empirical Risk:\n")
  print(obj$risk)
  cat("L:\n")
  print(obj$l.value)
  cat("Optimizer:\n")
  print(obj$optimizer)
}

#' Inherit predict function for \code{block.iso.lin.est}
#'
#' @param obj An object with S3 class \code{block.iso.lin.est}.
#'
#' @return \code{n1 x n2} numeric matrix of observed responses.
#'
#' @export
#' @examples
#' n1 <- 10; n2 <- 10
#' y <- matrix(rep(0, n1*n2), nrow = n1)
#' x <- matrix(rep(list(), n1*n2),nrow = n1, ncol =n2)
#' for (i in 1:n1){
#'   for (j in 1:n2){
#'     y[i,j] <- sin(10*norm(c(i/n1, j/n2), "2")) + rnorm(1, 0, 0.1)
#'     x[i,j] <- list(c(i/n1, j/n2))
#'   }
#' }
#' res.est <- block.iso.lin.est(x, y)
#' y.pred <- predict(res.est) # fitted value
predict.block.iso.lin.est <- function(obj){
  obj$y.hat
}
