#' Generate simulate data.
#'
#' @export
#' @param fun.list A list of linear and nonlinear function for link function.
#' @param n The number of sample in the generated data.
#' @param p The number of predictors in the generated data.
#' @param B.true The groundtruth of the coefficient matrix B.
#' @param snr The Signal Noise Ratio of data generation.
#' @param cor The corvariance of the generate data. Default value is 0.


data_generate <- function(fun.list, n, p, B.true, snr, cor=0){
  l <- length(fun.list)
  if(!cor){
    x <- matrix(rnorm(n*p, 0, 1), n, p)
  }else{
    mu <- rep(0, p)
    Sigma <- diag(p)
    Sigma <- cor^abs(row(Sigma)-col(Sigma))
    x <- mvrnorm(n, mu, Sigma)
  }
  y.true <- matrix(,n,l)
  y <- matrix(,n,l)

  for(i in 1:l){
    y.true[,i] <- fun.list[[i]](x%*%B.true[,i])
    sigma <- sqrt(var(y.true[,i]) / snr)
    y[,i] <- y.true[,i] + sigma*rnorm(n)
  }

  y.scaled <- scale(y)
  y.true.scaled <- scale(y.true, center=attributes(y.scaled)$'scaled:center', scale=attributes(y.scaled)$'scaled:scale')

  return(list(y.scaled=y.scaled, y=y, y.true.scaled=y.true.scaled, y.true=y.true, x=x))
}
