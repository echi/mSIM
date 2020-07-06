### Load
if("mgcv" %in% rownames(installed.packages()) == FALSE)
{install.packages("mgcv")
}
library("mgcv")

if("MASS" %in% rownames(installed.packages()) == FALSE)
{install.packages("MASS")
}
library("MASS")

if("gtools" %in% rownames(installed.packages()) == FALSE)
{install.packages("gtools")
}
library("gtools")

if("glmnet" %in% rownames(installed.packages()) == FALSE)
{install.packages("glmnet")
}
library("glmnet")

if("Matrix" %in% rownames(installed.packages()) == FALSE)
{install.packages("Matrix")
}
library("Matrix")

if("doParallel" %in% rownames(installed.packages()) == FALSE)
{install.packages("doParallel")
}
#library("doParallel")

if("parallel" %in% rownames(installed.packages()) == FALSE)
{install.packages("parallel")
}
library("parallel")

if("splines" %in% rownames(installed.packages()) == FALSE)
{install.packages("splines")
}
library("splines")


## data
get.data = function(fun.list, n, p, B.true, snr, cor=0){
  l = length(fun.list)
  if(!cor){
    x = matrix(rnorm(n*p, 0, 1), n, p)
  }else{
    mu = rep(0, p)
    Sigma = diag(p) 
    Sigma = cor^abs(row(Sigma)-col(Sigma))
    x = mvrnorm(n, mu, Sigma) 
  }
  y.true = matrix(,n,l)
  y = matrix(,n,l)
  
  for(i in 1:l){
    y.true[,i] = fun.list[[i]](x%*%B.true[,i]) 
    sigma = sqrt(var(y.true[,i]) / snr)
    y[,i] = y.true[,i] + sigma*rnorm(n)
  }
  
  y.scaled = scale(y)
  y.true.scaled=scale(y.true, center=attributes(y.scaled)$'scaled:center', scale=attributes(y.scaled)$'scaled:scale')
  
  return(list(y.scaled=y.scaled, y=y, y.true.scaled=y.true.scaled, y.true=y.true, x=x))
}

get.plot = function(y, y.true, x, alpha, beta){
  index.true = x%*%alpha
  index = x %*% beta
  plot(index.true, y, pch=16, cex=1, xlab='Index', ylab='y')
  
  od = order(index.true)
  lines(index.true[od], y.true[od], col='blue', lwd=2, lty=2)
  
  grid = seq(min(index),max(index),length=25)
  fit = si.smooth(y, x, as.vector(beta), k=10)
  h = si.h(fit)
  lines(grid, h(grid)$est, col='red', lwd=2)
}

get.plot.raw = function(y, x, beta, ylab, linear=F){
  index = x%*%beta
  indicator = (index < quantile(index, 0.95) & index > quantile(index, 0.05))
  plot(index[indicator], y[indicator], pch=16, cex=1, xlab='Index', ylab=ylab, ylim=c(-3,3))
  grid = seq(min(index),max(index),length=20)
  fit = si.smooth(y, x, as.vector(beta), k=5)
  h = si.h(fit)
  lines(grid, h(grid)$est, col='red', lwd=2)
  if(linear){
    fit.lm = lm(y~index-1)
    abline(fit.lm, col='blue', lwd=2)
  }
}

get.B.ridge = function(Y, X, lambda=NULL){
  # dim(B) = p*q; dim(Y) = n*q; dim(X) = n*p; length(lambda) = p
  n = dim(Y)[1]
  q = dim(Y)[2]
  p = dim(X)[2]
  B = matrix(,p,q)
  
  # ridge with CV
  for(j in 1:q){
    if(is.null(lambda)){
      temp = cv.glmnet(x=X, y=Y[,j], alpha=0)
      B[,j] = coef(temp, s=temp$lambda.min)[-1]
    }else{
      temp = glmnet(x=X, y=Y[,j], alpha=0, lambda=lambda)
      B[,j] = coef(temp, s=temp$lambda.min)[-1]
    }
  }
  
  B = col.norm(B)
  
  return(B)
}

col.norm = function(A){
  A = apply(A, 2, function(x) x/norm(x, '2'))
  return(A)
}

rank.norm = function(A, rank){
  decomp = svd(A)
  U = decomp$u
  U[abs(U)<1e-8] = 0
  V = decomp$v
  V[abs(V)<1e-8] = 0
  D = U[,1:rank] %*% diag(decomp$d[1:rank], rank) %*% t(V[,1:rank])
  return(D)
}


## P-spline
si.smooth = function(y, x, beta, k=7) {
  #beta = beta/sqrt(sum(beta^2))
  u = x%*%beta
  #fit = pspline.gam(y, u, k=k)
  fit = pspline.gam(y, u, k = k, knots.option='"equally-spaced"')
  return(fit)
}

si.h = function(fit) {
  
  stopifnot(class(fit)=="pspline.gam")
  
  par = fit$theta
  s.obj = fit$s.object
  knots = s.obj$knots
  
  coef = diff(c(0, as.numeric(par), 0))
  coef = 3*coef/diff(knots, lag=3)
  
  h = function(x) {
    est = as.vector(splineDesign(x = x, knots = knots, ord = 4, outer.ok=T)%*%par)
    der = as.vector(splineDesign(x = x, knots = knots, ord = 3, outer.ok=T)%*%coef)
    return(list(x = x, est = est, der = der))
  }
  
  return(h)
}

pspline.gam = function(y, x, bs="ps", k=7, knots.option="equally-spaced", lambda=NULL) {
  
  data = data.frame("x" = x)
  s.object = s(x, bs=bs, k = k, sp=lambda)
  knots = NULL
  if (knots.option == "quantile") knots = face::select.knots(x, knots = k-3, option="quantile")
  object = smooth.construct(s.object, data=data.frame(x=x), knots=list(x=knots))
  X = Matrix(object$X)
  S = object$S
  for (i in 1:length(S)) S[[i]] = Matrix(S[[i]])
  
  s.object = list(s.object=s.object,
                   knots = object$knots,
                   p.order = object$p.order,
                   knots.option = knots.option,
                   bs = bs)
  
  if(is.null(lambda)) {
    fit = gcv.wrapper(Y=y,X=object$X,P=object$S)
    lambda = fit$lambda
  } ##if
  
  temp = Matrix::tcrossprod(solve(Matrix::crossprod(X) + lambda*S[[1]]),X)
  theta = temp%*%y
  
  res = list(fitted.values = X%*%theta,theta=theta,
              lambda=lambda,
              s.object = s.object)
  class(res) = "pspline.gam"
  return(res)
}

gcv.wrapper = function(Y, X, P) {
  
  Lambda = exp(seq(-10, 10, length=21))
  K = length(P)
  if (K == 1) {
    return(gcv.criteria(Y, X, crossprod(X), P[[1]], Lambda))
  }
  if (K > 1) {
    List = combn(length(Lambda),K-1)
    fit_list = list(length=ncol(List))
    
    B0 = crossprod(X)
    
    for (i in 1:ncol(List)) {
      
      B = B0
      seq = List[,i]
      
      for (k in 1:(K-1)) {
        B = B + Lambda[seq[k]]*P[[k]]
      }
      
      fit_list[[i]] = gcv.criteria(Y, X, B, P[[K]], Lambda)
    }
    
    index = which.min(sapply(fit_list, function(x) x$gcv))
    res = fit_list[[index]]
    res$lambda = c(Lambda[List[,index]], res$lambda)
    
    return(res)
  }
}

gcv.criteria = function(Y, X, B, P, lambda) {
  
  EigenB = eigen(B)
  EigenB$values = sapply(EigenB$values, function(x) max(x, 1e-8))
  B_inv_half = EigenB$vectors%*%diag(sqrt(1/EigenB$values))%*%t(EigenB$vectors)
  
  S = B_inv_half%*%P%*%B_inv_half
  EigenS = eigen(S)
  
  s = EigenS$values
  s[s<1e-8] = 0
  
  A = X%*%B_inv_half%*%EigenS$vectors
  
  Y_square = sum(Y^2)
  AtY = t(A)%*%Y
  AtA = t(A)%*%A
  AtA_diag = diag(AtA)
  
  EigenA = eigen(AtA,symmetric=TRUE)
  EigenA$values = sapply(EigenA$values, function(x) max(x, 1e-8))
  
  
  AtA_half = EigenA$vectors%*%diag(sqrt(EigenA$values))%*%t(EigenA$vectors)
  
  gcv = function(x){
    d = 1/(1 + x*s)
    gcv = sum((AtA_half%*%(AtY*d))^2) - 2*sum(AtY^2*d) + Y_square
    gcv = gcv/(1-sum(d*AtA_diag)/length(Y))^2
    gcv
  }
  
  gcv_list = sapply(lambda,gcv)
  index = which.min(gcv_list)
  
  lambda.min = lambda[index]
  trace = sum(AtA_diag*(1/(1+lambda.min*s)))
  
  
  res = list("lambda"=lambda.min,
              "trace" = trace,
              "gcv" = gcv_list[index]
  )
  return(res)
}


## ADMM

get.B.ADMM = function(Y, X, B, lambda, rank, alpha=1, control1=list(max.iter=1e2, tol=1e-2), control2=list(ele.sparse=F, row.sparse=T, low.rank=T), select.method='linear', descent.method='bfgs', plot=F){
  if(control2$ele.sparse & control2$row.sparse) stop('Only one sparse penalty is supported now')
  n = dim(Y)[1]
  q = dim(Y)[2]
  p = dim(X)[2]
  
  if(select.method == 'linear'){
    ratio = c()
    for(P in 1:p){
      temp = 0
      for(Q in 1:q){
        temp = temp + (cor(Y[,Q], X[,P]))^2
      }
      ratio = c(ratio, temp)
    }
    fixed = which.max(ratio)
  }
  if(select.method == 'nonlinear'){
    ratio = c()
    for(P in 1:p){
      temp = 0
      for(Q in 1:q){
        SST = (n-1)*var(Y[,Q])
        nonlinear.fit = pspline.gam(Y[,Q], X[,P])
        temp = temp + sum((Y[,Q] - nonlinear.fit$fitted.values)^2)/SST
      }
      ratio = c(ratio, temp)
    }
    fixed = which.min(ratio)
  }
  
  C = B 
  A = B
  D = B
  W1 = matrix(0, p, q)
  W2 = matrix(0, p, q)
  W3 = matrix(0, p, q)
  converge = 0
  n.iter = 1
  pri.err.save = c()
  dual.err.save = c()
  
  while(n.iter <= control1$max.iter){
    temp = B
    
    # update A
    for(j in 1:q){
      eta = B[,j] - W1[,j] / alpha
      eta = eta / norm(eta, '2')
      if(eta[fixed] != 0){
        eta = sign(eta[fixed]) * eta
      }
      A[,j] = eta
    }
    
    # update C
    newY = B - W2 / alpha
    if(control2$ele.sparse){
      C = sign(newY)*pmax(abs(newY)- lambda / alpha, 0)
      C[fixed,] = newY[fixed,]
    }
    if(control2$row.sparse){
      for(k in 1:p){
        if(1){
          if(k != fixed){
            C[k,] = max(0, 1 - lambda/alpha/norm(newY[k,], '2')) * newY[k,]
          }else{
            C[k,] = newY[k,]
          }
        }
        if(0){
          C[k,] = max(0, 1 - lambda/alpha/norm(newY[k,], '2')) * newY[k,]
          if(k == fixed & max(0, 1 - lambda/alpha/norm(newY[k,], '2')) == 0){
            C[k,] = 1
          }
        }
      }
    }

    # update D
    if(control2$low.rank){
      D = rank.norm(B - W3 / alpha, rank)
    }
    
    # update B
    for(j in 1:q){
      y = Y[,j]
      x = X
      fit0 = si.smooth(y, x, C[,j])
      h0 = si.h(fit0)
      if(descent.method == 'gd'){ # not very useful
        b = grad.descent(y, x, B[,j], alpha, control2, h0, W1[,j], W2[,j], W3[,j], A[,j], C[,j], D[,j])
      }else if(descent.method == 'bfgs'){
        b = bfgs.descent(y, x, B[,j], alpha, h0, W1[,j], W2[,j], W3[,j], A[,j], C[,j], D[,j])
      }
      B[,j] = b$beta
    }
    
    # update W's
    W1 = W1 + alpha*(A - B)
    if(control2$ele.sparse | control2$row.sparse){
      W2 = W2 + alpha*(C - B)
    }
    if(control2$low.rank){
      W3 = W3 + alpha*(D - B)
    }

    # error
    pri.err = (sum((A - B)^2) + sum((C - B)^2)*(control2$ele.sparse | control2$row.sparse) + sum((D - B)^2)*control2$low.rank) / (p*q)
    dual.err = (alpha * norm(temp-B, 'F'))^2 / (p*q)
    pri.err.save = c(pri.err.save, pri.err)
    dual.err.save = c(dual.err.save, dual.err)
    error.control = (max(pri.err, dual.err) <= control1$tol)
    
    if(plot){
      par(mfrow=c(2,1))
      plot(pri.err.save, type='l')
      plot(dual.err.save, type='l')
    }
    
    if(error.control){
      converge=1
      break
    }
    
    n.iter = n.iter + 1
  }
  
  return(list(B.final=col.norm(rank.norm(C, rank)), B.sparse=C, B.proj=A, B.lowrank=D, converge=converge, iteration=n.iter, pri.err=pri.err.save, dual.err=dual.err.save, noPenIndex=fixed))
}

msim.pred = function(Y, X, B, Y.true, X.pred){
  Y.pred = matrix(, dim(X.pred)[1], dim(B)[2])
  for(j in 1:dim(B)[2]){
    beta = B[,j]
    fit = si.smooth(Y[,j], X, as.vector(beta))
    h = si.h(fit)
    Y.pred[,j] = h(X.pred%*%beta)$est
  }
  MSE = sum((Y.pred-Y.true)^2)/prod(dim(Y.true))
  return(list(Y.pred=Y.pred, MSE=MSE))
}

B.df = function(X, B, tuning){
  decomp = svd(B)
  
  if(1){
    if(tuning[2]==1){
      return(sum(decomp$u[,1] != 0))
    }
    else{
      B = decomp$u[,1:tuning[2]] %*% diag(decomp$d[1:tuning[2]])
    }
  }
  
  p = dim(B)[1]
  q = dim(B)[2]
  index = apply(B==0, 1, sum)
  index = (index != q)
  
  XI = NULL
  delta = diag(NULL)
  P = diag(NULL)
  
  for(j in 1:p){
    if(index[j]){
      XI = cbind(XI, X[,j])
      delta = bdiag(delta, diag(q) / norm(B[j,],'2'))
      P = bdiag(P, diag(q) - crossprod(t(B[j,])) / norm(B[j,],'2')^2)
    }
  }
  
  XTX = kronecker(t(XI) %*% XI, diag(q))
  inv = chol2inv(chol(XTX + tuning[1]*delta*P))
  df = sum(diag(inv%*%XTX))
  return(df + q*tuning[2] - tuning[2]*(1+tuning[2])/2)
}

B.df.nai = function(X, B, tuning){
  p = min(dim(B)[1], rankMatrix(X)[1])
  q = dim(B)[2]
  index = apply(B==0, 1, sum)
  index = (index != q)
  
  return((min(sum(index), p) + q - tuning[2])*tuning[2])
}

B.BIC = function(Y, X, B, tuning, linear=F){
  n = dim(Y)[1]
  q = dim(Y)[2]
  p = dim(B)[1]
  SSE = 0
  
  for(l in 1:q){
    if(linear==F){
      fit = si.smooth(Y[,l], X, B[,l])
      h = si.h(fit)
      Yl.est = h(X%*%B[,l])
      SSE = SSE + sum((Y[,l] - Yl.est$est)^2)
    }else{
      Yl.est = X%*%B[,l]
      SSE = SSE + sum((Y[,l] - Yl.est)^2)
    }
  }
  
  df = B.df(X, B, tuning)
  df.nai = B.df.nai(X, B, tuning)
  logSSE = log(SSE / (n*q))
  AIC = logSSE + df * 2 / (n*q)
  BIC = logSSE + df * log(n*q) / (n*q)
  BIC.nai = logSSE + df.nai * log(n*q) / (n*q)
  
  return(c(logSSE=logSSE, AIC=AIC, BIC=BIC, BIC.nai=BIC.nai, df=df, df.nai=df.nai))
}

B.BIC.new = function(Y, X, B, tuning, linear=F){
  n = dim(Y)[1]
  q = dim(Y)[2]
  p = dim(B)[1]
  SSE = 0
  
  for(l in 1:q){
    if(linear==F){
      fit = si.smooth(Y[,l], X, B[,l])
      h = si.h(fit)
      Yl.est = h(X%*%B[,l])
      SSE = SSE + sum((Y[,l] - Yl.est$est)^2)
    }else{
      Yl.est = X%*%B[,l]
      SSE = SSE + sum((Y[,l] - Yl.est)^2)
    }
  }
  
  #df = B.df(X, B, tuning)
  df.nai = B.df.nai(X, B, tuning)
  logMSE = log(SSE / (n*q))
  #AIC = logSSE + df * 2 / (n*q)
  #BIC = logSSE + df * log(n*q) / (n*q)
  BIC.nai = logMSE + df.nai * log(n*q) / (n*q)
  
  return(c(logMSE=logMSE, BIC.nai=BIC.nai, df.nai=df.nai))
}


bfgs.descent = function(y, x, beta, alpha, h0, omega1, omega2, omega3, eta, ksi, delta){
  n = dim(x)[1]
  p = dim(x)[2]
  
  grad = function(beta){
    h0.fit = h0(x%*%beta)
    f0 = h0.fit$est
    d0 = h0.fit$der
    return( t(x) %*% ((f0-y)*d0)/n - omega1 - omega2 - omega3 + alpha*(3*beta - eta - ksi - delta))
  }
  
  obj = function(beta){
    h0.fit = h0(x%*%beta)
    f0 = h0.fit$est
    #d0 = h0.fit$der
    return( norm(f0-y, '2')^2/(2*n) - crossprod(omega1+omega2+omega3, beta) + alpha*(norm(beta - eta, '2')^2 + norm(beta - ksi, '2')^2 + norm(beta - delta, '2')^2)/2)
  }
  
  output = optim(beta, obj, grad, method='BFGS')
  #output = lbfgs(obj, grad, beta, invisible=T, epsilon=1e-4)
  
  return(list(beta=output$par, converge=output$convergence, obj=output$value))
}

grad.descent = function(y, x, beta, alpha, control=list(max.iter=1e3, tol=1e-3, stepsize=1e-3), h0, omega1, omega2, eta, ksi){
  n = dim(x)[1]
  p = dim(x)[2]
  n.iter = 1
  converge = 0
  err.save = c()
  obj.save = c()
  
  while(n.iter <= control$max.iter){
    temp = beta
    grad = 0
    obj  = 0
    
    # calculate gradient
    for(i in 1:n){
      h0.fit = h0(x[i,]%*%beta)
      f0 = h0.fit$est
      d0 = h0.fit$der
      grad = grad + (f0 - y[i])*d0*x[i,]
      #obj = obj + (y[i] - f0)^2 
    }
    grad = grad / n - (omega1 + omega2) + alpha*(2*beta - eta - ksi)
    #obj  = obj / (2*n) + omega1 %*% (eta - b) + omega2 %*% (ksi - b) + alpha * (sum((eta - b)^2) + sum((ksi - b)^2)) / 2
    
    # update b
    beta = beta - control$stepsize*grad
    
    # error
    #err = norm(beta - temp, '2')
    err = max(abs(beta - temp))
    err.save = c(err.save, err)
    #obj.save = c(obj.save, obj)
    
    error.control = (err <= control$tol)
    if(error.control){
      converge=1
      break
    }
    
    n.iter = n.iter + 1
  }
  
  return(list(beta=beta, converge=converge, iteration=n.iter))#, err=err.save, obj=obj.save))
}

si.mgcv = function(alpha, y, x, opt=TRUE, k=10, fx=FALSE){
  ## Fit single index model using gam call, given theta (defines alpha). 
  ## Return ML if opt==TRUE and fitted gam with theta added otherwise.
  ## Suitable for calling from 'optim' to find optimal theta/alpha.
  kk = sqrt(sum(alpha^2))
  alpha = alpha/kk  ## so now ||alpha||=1
  a = x%*%alpha     ## argument of smooth
  b = gam(y~s(a,fx=fx,k=k),family=gaussian,method="ML") ## fit model
  if (opt) return(b$gcv.ubre) else {
    b$alpha = alpha  ## add alpha
    J = outer(alpha,-alpha/kk^2) ## compute Jacobian
    for (j in 1:length(alpha)) J[j+1,j] = J[j+1,j] + 1/kk
    b$J = J ## dalpha_i/dtheta_j 
    return(b)
  }
}

get.beta.naive = function(y, x, beta, alpha, k1=10, control1=list(max.iter=1e2, tol=1e-2), control2=list(max.iter=5e3, tol=1e-3, stepsize=1e-3), 
                           plot=F, select.method='NULL'){
  n = dim(x)[1]
  p = dim(x)[2]
  
  # ADMM initiate
  beta = as.vector(beta)
  #ksi = beta 
  eta = beta
  omega1 = rnorm(p, 0, 0.01)
  #omega2 = rnorm(p, 0, 0.01)
  
  n.iter = 1
  converge = 0
  pri.err.save = c()
  dual.err.save = c()
  
  while(n.iter <= control1$max.iter){
    temp = beta
    
    # update eta
    eta = beta - omega1 / alpha
    eta = eta / norm(eta, '2')
    eta = sign(eta[1]) * eta
    
    # update b
    # spline fit
    fit0 = si.smooth(y, x, eta, k=k1)
    h0 = si.h(fit0)
    b = grad.descent(y, x, beta, alpha, control2, h0, omega1, omega2=0, eta, ksi=0)
    beta = b$beta
    #print(round(beta, 3))
    
    # update omega's
    omega1 = omega1 + alpha*(eta - beta)
    
    # error
    pri.err = sum((eta - beta)^2) #+ sum((ksi - beta)^2))
    dual.err = (alpha * norm(temp-beta, '2'))^2
    pri.err.save = c(pri.err.save, pri.err)
    dual.err.save = c(dual.err.save, dual.err)
    error.control = (max(pri.err, dual.err) <= control1$tol)
    
    if(plot){
      par(mfrow=c(2,1))
      plot(pri.err.save, type='l')
      plot(dual.err.save, type='l')
    }
    
    if(error.control){
      converge=1
      break
    }
    
    n.iter = n.iter + 1
  }
  
  return(list(beta=beta/norm(beta,'2'), converge=converge, iteration=n.iter, pri.err=pri.err.save, dual.err=dual.err.save))
}

get.beta.ADMM = function(y, x, beta, lambda, alpha=1, control1=list(max.iter=1e2, tol=1e-2), select.method='linear', descent.method = 'bfgs', control2=list(max.iter=5e3, tol=1e-3, stepsize=1e-3), plot=F){
  n = dim(x)[1]
  p = dim(x)[2]
  
  # marginal screening
  if(select.method == 'linear'){
    ratio = c()
    for(i in 1:p){
      ratio = c(ratio, (cor(y, x[,i]))^2)
    }
    fixed = which.max(ratio)
  }
  
  if(select.method == 'nonlinear'){
    ratio = c()
    SST = (n-1)*var(y)
    for(i in 1:p){
      nonlinear.fit = pspline.gam(y, x[,i])
      ratio = c(ratio, sum((y - nonlinear.fit$fitted.values)^2)/SST)
    }
    fixed = which.min(ratio)
  }
  
  # ADMM initiate
  beta = as.vector(beta)
  ksi = beta 
  eta = beta
  omega1 = rep(0, p)
  omega2 = rep(0, p)
  
  n.iter = 1
  converge = 0
  pri.err.save = c()
  dual.err.save = c()
  
  while(n.iter <= control1$max.iter){
    temp = beta
    
    # update eta
    eta = beta - omega1 / alpha
    eta = eta / norm(eta, '2')
    eta = sign(eta[fixed]) * eta
    
    # update ksi
    newy = beta - omega2 / alpha
    ksi = sign(newy)*pmax(abs(newy) - lambda/alpha, 0)
    ksi[fixed] = newy[fixed]
    
    # update b
    # spline fit
    fit0 = si.smooth(y, x, ksi)
    h0 = si.h(fit0)
    if(descent.method == 'gd'){
      b = grad.descent(y, x, beta, alpha, control2, h0, omega1, omega2, eta, ksi)
    }else if(descent.method == 'bfgs'){
      b = bfgs.descent(y, x, beta, alpha, h0, omega1, omega2, eta, ksi)
    }
    
    #plot(b$beta, main='b')
    #plot(eta, main='eta')
    #plot(ksi, main='ksi')
    #plot(omega1)
    #plot(omega2)
    #plot(b$err, type='l')
    #invisible(readline(prompt="Press [enter] to continue"))
    #plot(b$err, type='l')
    #plot(b$obj, type='l')
    #Sys.sleep(1)
    beta = b$beta
    
    # update omega's
    omega1 = omega1 + alpha*(eta - beta)
    #plot(omega1, main='omega1')
    #Sys.sleep(1)
    omega2 = omega2 + alpha*(ksi - beta)
    #plot(omega2, main='omega2')
    #Sys.sleep(1)
    
    # error
    pri.err = (sum((eta - beta)^2) + sum((ksi - beta)^2))
    dual.err = (alpha * norm(temp-beta, '2'))^2
    pri.err.save = c(pri.err.save, pri.err)
    dual.err.save = c(dual.err.save, dual.err)
    error.control = (max(pri.err, dual.err) <= control1$tol)
    
    #if(n.iter %% 100 == 0)
    if(plot){
      par(mfrow=c(2,1))
      #cat('Lambda =', lambda, ';', 'Iteration =', n.iter, ';', '2-Norm =', norm(b, '2'), '\n')
      plot(pri.err.save, type='l')
      plot(dual.err.save, type='l')
    }
    
    if(error.control){
      converge=1
      break
    }
    
    n.iter = n.iter + 1
  }
  
  return(list(beta.sparse=ksi/norm(ksi, '2'), converge=converge, iteration=n.iter, noPenIndex=fixed, pri.err=pri.err.save, dual.err=dual.err.save))
  #return(list(beta=b, beta.sparse=ksi, beta.proj=eta, converge=converge, iteration=n.iter, noPenIndex=fixed))#pri.err=pri.err.save, dual.err=dual.err.save, noPenIndex=fixed))
}

beta.BIC = function(y, x, beta){
  n = length(y)
  fit = si.smooth(y, x, as.vector(beta))
  h = si.h(fit)
  h.new = h(x%*%beta)
  SSE = sum((y - h.new$est)^2)
  df = sum(beta != 0)
  logSSE = log(SSE / (n*q))
  AIC = logSSE + df * 2 / (n)
  BIC = logSSE + df * log(n) / n
  
  return(c(logSSE=logSSE, AIC=AIC, BIC=BIC, df=df))
}
