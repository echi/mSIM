#' @export data_generate

get.B.ridge = function(Y, X, lambda=NULL){
  # dim(B) = p*q; dim(Y) = n*q; dim(X) = n*p; length(lambda) = p
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  B <- matrix(,p,q)

  # ridge with CV
  for(j in 1:q){
    if(is.null(lambda)){
      temp <- cv.glmnet(x=X, y=Y[,j], alpha=0)
      B[,j] <- coef(temp, s=temp$lambda.min)[-1]
    }else{
      temp <- glmnet(x=X, y=Y[,j], alpha=0, lambda=lambda)
      B[,j] <- coef(temp, s=temp$lambda.min)[-1]
    }
  }

  B <- col.norm(B)

  return(B)
}


col.norm <- function(A){
  A <- apply(A, 2, function(x) x/norm(x, '2'))
  return(A)
}


gcv.wrapper <- function(Y, X, P) {

  Lambda <- exp(seq(-10, 10, length=21))
  K <- length(P)
  if (K == 1) {
    return(gcv.criteria(Y, X, crossprod(X), P[[1]], Lambda))
  }
  if (K > 1) {
    List = combn(length(Lambda),K-1)
    fit_list = list(length=ncol(List))

    B0 <- crossprod(X)

    for (i in 1:ncol(List)) {

      B <- B0
      seq <- List[,i]

      for (k in 1:(K-1)) {
        B <- B + Lambda[seq[k]]*P[[k]]
      }

      fit_list[[i]] <- gcv.criteria(Y, X, B, P[[K]], Lambda)
    }

    index <- which.min(sapply(fit_list, function(x) x$gcv))
    res <- fit_list[[index]]
    res$lambda <- c(Lambda[List[,index]], res$lambda)

    return(res)
  }
}

gcv.criteria <- function(Y, X, B, P, lambda) {

  EigenB <- eigen(B)
  EigenB$values <- sapply(EigenB$values, function(x) max(x, 1e-8))
  B_inv_half <- EigenB$vectors%*%diag(sqrt(1/EigenB$values))%*%t(EigenB$vectors)

  S <- B_inv_half%*%P%*%B_inv_half
  EigenS <- eigen(S)

  s <- EigenS$values
  s[s<1e-8] <- 0

  A <- X%*%B_inv_half%*%EigenS$vectors

  Y_square <- sum(Y^2)
  AtY <- t(A)%*%Y
  AtA <- t(A)%*%A
  AtA_diag <- diag(AtA)

  EigenA <- eigen(AtA,symmetric=TRUE)
  EigenA$values <- sapply(EigenA$values, function(x) max(x, 1e-8))


  AtA_half <- EigenA$vectors%*%diag(sqrt(EigenA$values))%*%t(EigenA$vectors)

  gcv <- function(x){
    d = 1/(1 + x*s)
    gcv = sum((AtA_half%*%(AtY*d))^2) - 2*sum(AtY^2*d) + Y_square
    gcv = gcv/(1-sum(d*AtA_diag)/length(Y))^2
    gcv
  }

  gcv_list <- sapply(lambda,gcv)
  index <- which.min(gcv_list)

  lambda.min <- lambda[index]
  trace <- sum(AtA_diag*(1/(1+lambda.min*s)))


  res <- list("lambda"=lambda.min,
             "trace" = trace,
             "gcv" = gcv_list[index]
  )
  return(res)
}

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

