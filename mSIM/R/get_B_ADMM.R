#' The main function of mSIM package
#'
#' @param Y The index matrix
#' @param X The covariate matrix
#' @param B The initial coefficient matrix
#' @param lambda The regualrize coefficient
#' @param rank The rank of B matrix
#' @param alpha The step in the optimization process
#' @param control1 A list that contains parameters for minimize loss function
#' @param control2 A list that contains parameters for regularize
#' @param select.method link funtion selection method
#' @param descend.method optimization method
#' @param plot whether or not plot the error during the iteration
#' @return A list that contains information such as the final B matrix, the training error.
#' @examples
#' test = get.B.ADMM(Y=Y, X=X, B=B.init, lambda=0.2, rank=3, alpha=1, control1 = list(max.iter=5e1, tol=1e-5), control2=list(ele.sparse=F, row.sparse=T, low.rank=T), select.method='linear', plot=F)


get_B_ADMM = function(Y, X, B = NULL, lambda, rank, alpha=1, control1=list(max.iter=1e2, tol=1e-2), control2=list(ele.sparse=F, row.sparse=T, low.rank=T), select.method='linear', descent.method='bfgs', plot=F){
  Y <- scale(Y)
  X <- scale(X)

  if(control2$ele.sparse & control2$row.sparse) stop('Only one sparse penalty is supported now')
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]

  if (is.null(B)){
    B <- get.B.ridge(Y, X)
  }

  if(select.method == 'linear'){
    ratio <- c()
    for(P in 1:p){
      temp <- 0
      for(Q in 1:q){
        temp <- temp + (cor(Y[,Q], X[,P]))^2
      }
      ratio <- c(ratio, temp)
    }
    fixed <- which.max(ratio)
  }
  if(select.method == 'nonlinear'){
    ratio <- c()
    for(P in 1:p){
      temp <- 0
      for(Q in 1:q){
        SST <- (n-1)*var(Y[,Q])
        nonlinear.fit <- pspline.gam(Y[,Q], X[,P])
        temp <- temp + sum((Y[,Q] - nonlinear.fit$fitted.values)^2)/SST
      }
      ratio <- c(ratio, temp)
    }
    fixed <- which.min(ratio)
  }

  C <- B
  A <- B
  D <- B
  W1 <- matrix(0, p, q)
  W2 <- matrix(0, p, q)
  W3 <- matrix(0, p, q)
  converge <- 0
  n.iter <- 1
  pri.err.save <- c()
  dual.err.save <- c()

  while(n.iter <= control1$max.iter){    ##2-step iteration
    temp <- B

    # update A
    for(j in 1:q){
      eta <- B[,j] - W1[,j] / alpha
      eta <- eta / norm(eta, '2')
      if(eta[fixed] != 0){
        eta <- sign(eta[fixed]) * eta
      }
      A[,j] <- eta
    }

    # update C
    newY <- B - W2 / alpha
    if(control2$ele.sparse){
      C <- sign(newY)*pmax(abs(newY)- lambda / alpha, 0)
      C[fixed,] <- newY[fixed,]
    }
    if(control2$row.sparse){
      for(k in 1:p){
        if(1){
          if(k != fixed){
            C[k,] <- max(0, 1 - lambda/alpha/norm(newY[k,], '2')) * newY[k,]
          }else{
            C[k,] <- newY[k,]
          }
        }
        if(0){
          C[k,] <- max(0, 1 - lambda/alpha/norm(newY[k,], '2')) * newY[k,]
          if(k == fixed & max(0, 1 - lambda/alpha/norm(newY[k,], '2')) == 0){
            C[k,] <- 1
          }
        }
      }
    }

    # update D
    if(control2$low.rank){
      D <- rank.norm(B - W3 / alpha, rank)
    }

    # update B
    for(j in 1:q){
      y <- Y[,j]
      x <- X
      fit0 <- si.smooth(y, x, C[,j])
      h0 <- si.h(fit0)
      if(descent.method == 'gd'){ # not very useful
        b <- grad.descent(y, x, B[,j], alpha, control2, h0, W1[,j], W2[,j], W3[,j], A[,j], C[,j], D[,j])
      }else if(descent.method == 'bfgs'){
        b <- bfgs.descent(y, x, B[,j], alpha, h0, W1[,j], W2[,j], W3[,j], A[,j], C[,j], D[,j])
      }
      B[,j] <- b$beta
    }

    # update W's
    W1 <- W1 + alpha*(A - B)
    if(control2$ele.sparse | control2$row.sparse){
      W2 <- W2 + alpha*(C - B)
    }
    if(control2$low.rank){
      W3 <- W3 + alpha*(D - B)
    }

    # error
    pri.err <- (sum((A - B)^2) + sum((C - B)^2)*(control2$ele.sparse | control2$row.sparse) + sum((D - B)^2)*control2$low.rank) / (p*q)
    dual.err <- (alpha * norm(temp-B, 'F'))^2 / (p*q)
    pri.err.save <- c(pri.err.save, pri.err)
    dual.err.save <- c(dual.err.save, dual.err)
    error.control <- (max(pri.err, dual.err) <= control1$tol)

    if(plot){
      par(mfrow <- c(2,1))
      plot(pri.err.save, type <- 'l')
      plot(dual.err.save, type <- 'l')
    }

    if(error.control){
      converge <- 1
      break
    }

    n.iter <- n.iter + 1
  }

  return(list(B.final <- col.norm(rank.norm(C, rank)), B.sparse <- C, B.proj <- A, B.lowrank <- D, converge <- converge, iteration <- n.iter, pri.err <- pri.err.save, dual.err <- dual.err.save, noPenIndex <- fixed))
}
