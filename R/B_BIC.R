#' Bayesian Information Criteria for parameter selection.
#' @export B_BIC
#' @param Y The index matrix
#' @param X The covariate matrix
#' @param B The initial coefficient matrix
#' @param tuning A coefficient matrix that need to be tuned.
#' @param linear Whether the link function is linear or not. The default value is False.
#' @return A list that contains information about the Log Mean Square Error and the Degree of Freedom.

B_BIC <- function(Y, X, B, tuning, linear=FALSE){
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(B)[1]
  SSE <- 0

  for(l in 1:q){
    if(linear==F){
      fit <- si.smooth(Y[,l], X, B[,l])
      h <- si.h(fit)
      Yl.est <- h(X%*%B[,l])
      SSE <- SSE + sum((Y[,l] - Yl.est$est)^2)
    }else{
      Yl.est <- X%*%B[,l]
      SSE <- SSE + sum((Y[,l] - Yl.est)^2)
    }
  }

  df.nai <- B.df.nai(X, B, tuning)
  logMSE <- log(SSE / (n*q))
  BIC.nai <- logMSE + df.nai * log(n*q) / (n*q)

  return(c(logMSE = logMSE, BIC.nai = BIC.nai, df.nai = df.nai))
}
