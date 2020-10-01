B.df <- function(X, B, tuning){
  decomp <- svd(B)


  if(tuning[2]==1){
    return(sum(decomp$u[,1] != 0))
  }
  else{
    B <- decomp$u[,1:tuning[2]] %*% diag(decomp$d[1:tuning[2]])
  }


  p <- dim(B)[1]
  q <- dim(B)[2]
  index <- apply(B==0, 1, sum)
  index <- (index != q)

  XI <- NULL
  delta <- diag(NULL)
  P <- diag(NULL)

  for(j in 1:p){
    if(index[j]){
      XI <- cbind(XI, X[,j])
      delta <- bdiag(delta, diag(q) / norm(B[j,],'2'))
      P <- bdiag(P, diag(q) - crossprod(t(B[j,])) / norm(B[j,],'2')^2)
    }
  }

  XTX <- kronecker(t(XI) %*% XI, diag(q))
  inv <- chol2inv(chol(XTX + tuning[1]*delta*P))
  df <- sum(diag(inv%*%XTX))
  return(df + q*tuning[2] - tuning[2]*(1+tuning[2])/2)
}

B.df.nai <- function(X, B, tuning){
  p <- min(dim(B)[1], rankMatrix(X)[1])
  q <- dim(B)[2]
  index <- apply(B==0, 1, sum)
  index <- (index != q)

  return((min(sum(index), p) + q - tuning[2])*tuning[2])
}
