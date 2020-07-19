#' @export data_generate

get.B.ridge = function(Y, X, lambda=NULL){        ##岭回归，用于B矩阵的初始化
  # dim(B) = p*q; dim(Y) = n*q; dim(X) = n*p; length(lambda) = p
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  B <- matrix(,p,q)

  # ridge with CV
  for(j in 1:q){
    if(is.null(lambda)){
      temp <- cv.glmnet(x=X, y=Y[,j], alpha=0)     ##广义线性回归
      B[,j] <- coef(temp, s=temp$lambda.min)[-1]
    }else{
      temp <- glmnet(x=X, y=Y[,j], alpha=0, lambda=lambda)     ##广义线性回归
      B[,j] <- coef(temp, s=temp$lambda.min)[-1]
    }
  }

  B <- col.norm(B)

  return(B)
}


col.norm <- function(A){
  A <- apply(A, 2, function(x) x/norm(x, '2'))    #apply函数，将一个函数作用到array/matrix的margin(行或者列)
  return(A)
}


gcv.wrapper <- function(Y, X, P) {   #交叉验证法的某一个部分，具体什么部分不知道？

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

gcv.criteria <- function(Y, X, B, P, lambda) {   ##交叉验证衡量标准

  EigenB <- eigen(B)
  EigenB$values <- sapply(EigenB$values, function(x) max(x, 1e-8))   ##作用与apply相同，将结果整理以向量，矩阵，列表 的形式输出
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

data_generate <- function(fun.list, n, p, B.true, snr, cor=0){   ##产生数据
  l <- length(fun.list)
  if(!cor){
    x <- matrix(rnorm(n*p, 0, 1), n, p)
  }else{
    mu <- rep(0, p)  ##重复p次
    Sigma <- diag(p)
    Sigma <- cor^abs(row(Sigma)-col(Sigma))
    x <- mvrnorm(n, mu, Sigma)    ##产生服从高斯分布的随机数
  }
  y.true <- matrix(,n,l)
  y <- matrix(,n,l)

  for(i in 1:l){
    y.true[,i] <- fun.list[[i]](x%*%B.true[,i])   ## %*%矩阵乘法
    sigma <- sqrt(var(y.true[,i]) / snr)
    y[,i] <- y.true[,i] + sigma*rnorm(n)    ## y是有误差的label， y.true是无误差label
  }

  y.scaled <- scale(y)
  y.true.scaled <- scale(y.true, center=attributes(y.scaled)$'scaled:center', scale=attributes(y.scaled)$'scaled:scale')   ##attributes添加属性关键字

  return(list(y.scaled=y.scaled, y=y, y.true.scaled=y.true.scaled, y.true=y.true, x=x))
}

