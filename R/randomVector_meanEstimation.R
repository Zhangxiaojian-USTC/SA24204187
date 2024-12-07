#' @import nloptr
#' @import numDeriv
#' @import knitr
#' @import rmarkdown
#' @import Rcpp
#' @import RcppArmadillo
#' @import stats
#' @import xtable
#' @import MASS
#' @import boot
#' @import DAAG
#' @import bootstrap
#' @import stats4
#' @import purrr
#' @import Rglpk
#' @import microbenchmark
#' @import twosamples
#' @useDynLib SA24204187
NULL

#' @title The empirical mean using R
#' @description The empirical mean using R
#' @param matrix Design matrix over \eqn{R^{p \times n}}, \code{n} is the number of random vectors, \code{p} is the dimension of random vectors
#' @return The empirical mean of \code{n} column vectors of size \code{p} of the matrix
#' @examples
#' \dontrun{
#'     random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
#'     empirical.mean(random_matrix_normal)
#' }
#' @export
empirical.mean <- function(matrix) rowMeans(matrix)

#' @title The Catoni–Giulini mean estimator using R
#' @description The Catoni–Giulini mean estimator using R
#' @param matrix Design matrix over \eqn{R^{p \times n}}, \code{n} is the number of random vectors, \code{p} is the dimension of random vectors
#' @param alpha A small and positive Catoni–Giulini parameter
#' @return The Catoni–Giulini mean estimator of \code{n} column vectors of size \code{p} of the matrix
#' @examples
#' \dontrun{
#'     random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
#'     CG.mean(random_matrix_normal)
#' }
#' @export
CG.mean <- function(matrix,alpha=0.1){
  n <- ncol(matrix) 
  new_matrix <- as.matrix(matrix)
  for (i in 1:n) {
    l2_norm <- norm(as.matrix(matrix[,i]), type = "2")
    M <- min(1,1/(alpha*sqrt(l2_norm)))
    new_matrix[,i] <- matrix[,i]*M
  }
  output <- rowMeans(new_matrix); return(output)
}

#' @title The Catoni–Giulini mean estimator selected by MSE using R
#' @description The Catoni–Giulini mean estimator selected by MSE using R
#' @param matrix Design matrix over \eqn{R^{p \times n}}, \code{n} is the number of random vectors, \code{p} is the dimension of random vectors
#' @param delta_CE A small and positive parameter meaning a confident level like 0.05
#' @param level The level of precision of the selection
#' @return The Catoni–Giulini mean estimator of \code{n} column vectors of size \code{p} of the matrix selected by MSE
#' @examples
#' \dontrun{
#'     random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
#'     selectCG.mean(random_matrix_normal)
#' }
#' @export
selectCG.mean <- function(matrix,delta_CE = 0.05,level=10){
  n <- ncol(matrix); MSE <-  numeric(level)
  coff.MSE <- l2_norm <- numeric(n)
  for (i in 1:n) l2_norm[i] <- norm(as.matrix(matrix[,i]), type = "2")
  M <- sqrt(max(l2_norm))
  alpha <- sqrt(log(1/delta_CE)/n)/M
  times <- M/alpha
  eff <- 1:floor(times) / level
  for (i in 1:level) {
    for (j in 1:n) {
      coff.MSE[j] <- norm(as.matrix(matrix[,j]-CG.mean(matrix,eff[i]*alpha)), type = "2")
    }
    MSE[i] <- sum(coff.MSE)
  }
  MIN <- which.min(MSE)
  return(CG.mean(matrix,eff[MIN]*alpha))
}

#' @title The simple trimmed mean estimator using R
#' @description The simple trimmed mean estimator using R
#' @param matrix Design matrix over \eqn{R^{p \times n}}, \code{n} is the number of random vectors, \code{p} is the dimension of random vectors
#' @param epsilon Estimated proportion of outliers in a dataset
#' @return The simple trimmed mean estimator of \code{n} column vectors of size \code{p} of the matrix
#' @examples
#' \dontrun{
#'     random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
#'     trimmed.mean(random_matrix_normal)
#' }
#' @export
trimmed.mean <- function(matrix,epsilon=0.1){
  n <- ncol(matrix)
  p <- nrow(matrix)
  index <- sample(1:n,n,replace = FALSE)
  simulation.matrix <- matrix[,index[1:floor(n/2)]]
  test.matrix <- matrix[,index[(floor(n/2)+1):n]]
  dist <- numeric(floor(n/2))
  coff.dist <- numeric(n)
  for (i in 1:floor(n/2)) {
    for (j in 1:n) {
      coff.dist[j] <- norm(as.matrix(simulation.matrix[i]-matrix[,j]), type = "2")
    } 
    dist[i] <- sum(coff.dist)
  } 
  a=sort(dist); u=a[ceiling(epsilon*length(a))]; result=a[a<u]
  M <- max(result)
  dist.star <- numeric(n - floor(n/2))
  for (i in (floor(n/2)+1):n) {
    for (j in 1:n) {
      coff.dist[j] <- norm(as.matrix(test.matrix[i]-matrix[,j]), type = "2")
    } 
    dist.star[(i-floor(n/2))] <- sum(coff.dist)
  } 
  output <- test.matrix[,dist.star<M]
  output.median <- rowMeans(output)
  return(output.median)
} # epsilon = 32log(8/δ)/(3n), epsilon*n: deleted data number

huber_loss <- function(tau = 1.0, vector){
  l2_norm <- norm(as.matrix(vector), type = "2")
  if (l2_norm <= tau) {
    re <- l2_norm^2 / 2
  } else {
    re <- tau * l2_norm - tau^2 / 2
  } 
  return(re)
}

bisquare_loss <- function(tau = 1.0,vector){
  l2_norm <- norm(as.matrix(vector), type = "2")
  if (l2_norm <= tau) {
    re <- 1 - (1 - (l2_norm/tau)^2)^3
  } else re <- 1
  return(re)
}

wang_loss <- function(tau=1.0,vector){
  l2_norm <- norm(as.matrix(vector), type = "2")
  re <- 1 - exp(-l2_norm^2/tau); return(re)
}

#' @title The mean estimator with robust loss functions using R
#' @description The mean estimator with robust loss functions using R
#' @param matrix Design matrix over \eqn{R^{p \times n}}, \code{n} is the number of random vectors, \code{p} is the dimension of random vectors
#' @param robust.par the parameter of the selected loss function
#' @param method the special robust loss function name you choose, such as huber, bisquare and wang
#' @return The mean estimator with robust loss functions of \code{n} column vectors of size \code{p} of the matrix
#' @examples
#' \dontrun{
#'     random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
#'     robust.mean(random_matrix_normal)
#' }
#' @export
robust.mean <- function(matrix,robust.par=1.0,method = "huber"){
  n <- ncol(matrix)
  if (method == "huber") {
    huber.loss <- function(x){
      loss <- numeric(n)
      for (i in 1:n) {
        diff <- x-matrix[,i]
        loss[i] <- huber_loss(robust.par,diff)
      } 
      return(sum(loss))
    }
    gradient_function <- function(x) grad(huber.loss, x)
    res <- nloptr(x0 = rowMeans(matrix),eval_f = huber.loss,
                  eval_grad_f = gradient_function,
                  opts = list(algorithm = "NLOPT_LD_MMA",  # 选择优化算法
                              xtol_rel = 1.0e-8,           # 设定相对容差
                              maxeval = 1000             # 最大评估次数
                  ))
    return(res$solution)
  }
  if (method == "bisquare") {
    bisquare.loss <- function(x){
      loss <- numeric(n)
      for (i in 1:n) {
        diff <- x-matrix[,i]
        loss[i] <- bisquare_loss(robust.par,diff)
      } 
      return(sum(loss))
    }
    gradient_function <- function(x) grad(bisquare.loss, x)
    res <- nloptr(x0 = rowMeans(matrix),eval_f = bisquare.loss,
                  eval_grad_f = gradient_function,
                  opts = list(algorithm = "NLOPT_LD_MMA",  # 选择优化算法
                              xtol_rel = 1.0e-8,           # 设定相对容差
                              maxeval = 1000             # 最大评估次数
                  ))
    return(res$solution)
  }
  if (method == "wang") {
    wang.loss <- function(x){
      loss <- numeric(n)
      for (i in 1:n) {
        diff <- x-matrix[,i]
        loss[i] <- wang_loss(robust.par,diff)
      } 
      return(sum(loss))
    }
    gradient_function <- function(x) grad(wang.loss, x)
    res <- nloptr(x0 = rowMeans(matrix),eval_f = wang.loss,
                  eval_grad_f = gradient_function,
                  opts = list(algorithm = "NLOPT_LD_MMA",  # 选择优化算法
                              xtol_rel = 1.0e-8,           # 设定相对容差
                              maxeval = 1000             # 最大评估次数
                  ))
    return(res$solution)
  }
}

#' @title The median-of-mean estimator using R
#' @description The median-of-mean estimator using R
#' @param matrix Design matrix over \eqn{R^{p \times n}}, \code{n} is the number of random vectors, \code{p} is the dimension of random vectors
#' @param k the number of groups you want after the data is evenly grouped
#' @param method the special definition of "median" you choose, such as coordinate-wise, Euclidean-ball and geometric
#' @return The median-of-mean estimator of \code{n} column vectors of size \code{p} of the matrix
#' @examples
#' \dontrun{
#'     random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
#'     robust.mean(random_matrix_normal)
#' }
#' @export
MOM.mean <- function(matrix,k,method = "geometric"){
  n <- ncol(matrix)
  p <- nrow(matrix)
  if (k > floor(n/2)) {
    k = floor(n/2)
  } else if (k < 2) {
    k = 2
  } else {
    k <- k
  }
  index <- sample(1:n, size = n, replace = FALSE); remainder = n - k*floor(n/k)
  a <- vector("list", length = k)
  for (i in 1:k) {
    if (i <= remainder) {
      a[[i]] = matrix[,c(index[((i-1)*floor(n/k)+1):(i*floor(n/k))],(k*floor(n/k)+i))]
    } else {
      a[[i]] = matrix[,index[((i-1)*floor(n/k)+1):(i*floor(n/k))]]
    } 
  }
  b <- vector("list", length = k)
  output.mean <- matrix(0,nrow = p,ncol = k)
  for (i in 1:k) {
    b[[i]] <- rowMeans(a[[i]])
    output.mean[,i] <- b[[i]]
  } 
  output.median <- numeric(p)
  if (method == "coordinate-wise") {
    for (j in 1:p) output.median[j] <- median(output.mean[j,])
  } # k = ceiling(8/log(1/δ))
  else if (method == "Euclidean-ball") {
    ball.radius <- function(u){
      dist <- numeric(k)
      for (j in 1:k) {
        dist[j] <- norm(as.matrix(u-output.mean[,j]), type = "2")
      }
      radius <- median(dist)
      return(radius)
    }
    gradient_function <- function(x) {
      grad(ball.radius, x)
    }
    res <- nloptr(x0 = rowMeans(output.mean),eval_f = ball.radius,
                  eval_grad_f = gradient_function,
                  opts = list(algorithm = "NLOPT_LD_MMA",  # 选择优化算法
                              xtol_rel = 1.0e-8,           # 设定相对容差
                              maxeval = 1000             # 最大评估次数
                  ))
    output.median <- res$solution
  } # k = ceiling(8/log(1/δ))
  else if (method == "geometric") {
    ball.error <- function(u){
      dist <- numeric(k)
      for (j in 1:k) {
        dist[j] <- norm(as.matrix(u-output.mean[,j]), type = "2")
      }
      error <- sum(dist)
      return(error)
    }
    gradient_function <- function(x) {
      grad(ball.error, x)
    }
    res <- nloptr(x0 = rowMeans(output.mean),eval_f = ball.error,
                  eval_grad_f = gradient_function,
                  opts = list(algorithm = "NLOPT_LD_MMA",  # 选择优化算法
                              xtol_rel = 1.0e-8,           # 设定相对容差
                              maxeval = 1000             # 最大评估次数
                  ))
    output.median <- res$solution
  }# k = ceiling(8/log(1/δ))
  return(output.median)
}

#' @title A Gibbs sampler using R
#' @description A Gibbs sampler using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' rnR <- gibbsSamplingR(100,10)
#' par(mfrow=c(2,1));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
gibbsSamplingR <- function(N, thin) {
  mat <- matrix(nrow = N, ncol = 2)
  n = 10; a = 1; b = 1; x = n/2; y = 0.5
  for (i in 1:N) {
    for (j in 1:thin) {
      x = rbinom(1, n, y)
      y = rbeta(1, x + a,n - x + b)
    }
    mat[i, ] <- c(x, y)
  }
  mat
}
