library(devtools)
setwd("../npparam/")
load_all()
library(npparam)
source("../sim-npparam/synthetic.R")
source("../sim-npparam/estimators.R")


simulate.one <- function(n){
  data <- additive.lipscitz(sample.size = n, sigma = 0.1)
  cat("Fitting additive isotonic for sample size", n, "\n")
  add.est <- additive.iso.est(data$x, data$y)
  cat("Fitting additive sieve for sample size", n, "\n")
  sieve.est <- .fit.additive.sieve(data$x, data$y)
  cat("Fitting GBM for sample size", n, "\n")
  gbm.est <- .fit.additive.gbm(data$x, data$y)
  gam.est <- .fit.gam.2d(data$x, data$y)
  if(n > 2000){
    cat("Skipping KRR for sample size", n, "\n")
    krr.est <- NULL
  }else{
    cat("Fitting KRR for sample size", n, "\n")
    krr.est <- .fit.additive.krr(data$x, data$y)
  }
  out <- list("add.est"=add.est,
              "sieve.est"=sieve.est,
              "krr.est"=krr.est,
              "gbm.est"=gbm.est,
              "gam.est"=gam.est)
  out
}

run.experiment <- function(seed, n){
  set.seed(seed)
  test.x <- matrix(runif(1e4*2, 0, 1), ncol=2)
  data <- additive.lipscitz(sample.size = 10, sigma = 0.1)
  test.y <- data$fn(test.x)
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  x <- c("n", "seed", "risk", "estimator", "L")
  colnames(df) <- x

  res <- replicate(10, simulate.one(n), simplify = FALSE)

  est.risk <- unlist(lapply(res, function(m){
    mean((predict(m$add.est, test.x)-test.y)^2)}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="additive iso",
                       L=0)
  df <- rbind(df, out.df)

  est.risk <- unlist(lapply(res, function(m){
    mean((m$sieve.est(test.x)-test.y)^2)}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="sieve",
                       L=0)
  df <- rbind(df, out.df)

  est.risk <- unlist(lapply(res, function(m){
    mean((m$gam.est(test.x)-test.y)^2)}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="gam",
                       L=0)
  df <- rbind(df, out.df)

  est.risk <- unlist(lapply(res, function(m){
    mean((m$gbm.est(test.x)-test.y)^2)}))
  out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="gbm",
                       L=0)
  df <- rbind(df, out.df)

  if(n <= 2000){
    est.risk <- unlist(lapply(res, function(m){
      mean((m$krr.est(test.x)-test.y)^2)}))
    out.df <- data.frame(n=n, seed=seed, risk=est.risk, estimator="krr",
                         L=0)
    df <- rbind(df, out.df)
  }

  return(df)
}

# call Rscript cross_validation.R 100 500
args <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(args[1]) # seed
cat("Base seed = ", seed, "\n")

sample.n <- as.numeric(args[2]) # seed
cat("Sample size = ", sample.n, "\n")
base.dir <- "../sim-npparam/results"
save.dir <- file.path(base.dir, "add-exp1-5d")

dir.create(save.dir)
out.df <- run.experiment(seed, sample.n)

filename <- paste0("seed", seed, "N", sample.n, ".RData")
save(out.df, file=file.path(save.dir, filename))
