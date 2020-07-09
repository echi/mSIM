pspline.gam <- function(y, x, bs="ps", k=7, knots.option="equally-spaced", lambda=NULL) {   ##样条插值的具体算法

  data <- data.frame("x" = x)
  s.object <- s(x, bs <- bs, k <- k, sp <- lambda)
  knots = NULL
  if (knots.option == "quantile"){
    knots <- face::select.knots(x, knots <- k-3, option <- "quantile")
  }
  object <- smooth.construct(s.object, data <- data.frame(x <- x), knots <- list(x <- knots))
  X <- Matrix(object$X)
  S <- object$S
  for (i in 1:length(S)){
    S[[i]] <- Matrix(S[[i]])
  }

  s.object <- list(s.object <- s.object,
                  knots <- object$knots,
                  p.order <- object$p.order,
                  knots.option <- knots.option,
                  bs <- bs)

  if(is.null(lambda)) {
    fit <- gcv.wrapper(Y <- y,X <- object$X,P <- object$S)
    lambda <- fit$lambda
  } ##if

  temp <- Matrix::tcrossprod(solve(Matrix::crossprod(X) + lambda*S[[1]]),X)
  theta <- temp%*%y

  res <- list(fitted.values <- X%*%theta,theta <- theta,
             lambda <- lambda,
             s.object <- s.object)
  class(res) <- "pspline.gam"
  return(res)
}


rank.norm <- function(A, rank){
  decomp <- svd(A)
  U <- decomp$u
  U[abs(U)<1e-8] <- 0
  V <- decomp$v
  V[abs(V)<1e-8] <- 0
  D <- U[,1:rank] %*% diag(decomp$d[1:rank], rank) %*% t(V[,1:rank])
  return(D)
}

si.smooth <- function(y, x, beta, k=7) {    ##样条插值
  #beta = beta/sqrt(sum(beta^2))
  u <- x%*%beta
  #fit = pspline.gam(y, u, k=k)
  fit <- pspline.gam(y, u, k <- k, knots.option <- '"equally-spaced"')
  return(fit)
}

si.h <- function(fit) {

  stopifnot(class(fit)=="pspline.gam")

  par <- fit$theta
  s.obj <- fit$s.object
  knots <- s.obj$knots

  coef <- diff(c(0, as.numeric(par), 0))
  coef <- 3*coef/diff(knots, lag=3)

  h <- function(x) {
    est <- as.vector(splineDesign(x <- x, knots <- knots, ord <- 4, outer.ok <- T)%*%par)
    der = as.vector(splineDesign(x <- x, knots <- knots, ord <- 3, outer.ok <- T)%*%coef)
    return(list(x <- x, est <- est, der <- der))
  }

  return(h)
}

bfgs.descent <- function(y, x, beta, alpha, h0, omega1, omega2, omega3, eta, ksi, delta){   ##bfgs
  n <- dim(x)[1]
  p <- dim(x)[2]

  grad <- function(beta){
    h0.fit <- h0(x%*%beta)
    f0 <- h0.fit$est
    d0 <- h0.fit$der
    return( t(x) %*% ((f0-y)*d0)/n - omega1 - omega2 - omega3 + alpha*(3*beta - eta - ksi - delta))
  }

  obj <- function(beta){
    h0.fit <- h0(x%*%beta)
    f0 <- h0.fit$est
    #d0 = h0.fit$der
    return( norm(f0-y, '2')^2/(2*n) - crossprod(omega1+omega2+omega3, beta) + alpha*(norm(beta - eta, '2')^2 + norm(beta - ksi, '2')^2 + norm(beta - delta, '2')^2)/2)
  }

  output <- optim(beta, obj, grad, method='BFGS')
  #output = lbfgs(obj, grad, beta, invisible=T, epsilon=1e-4)

  return(list(beta <- output$par, converge <- output$convergence, obj <- output$value))
}

grad.descent <- function(y, x, beta, alpha, control = list(max.iter <- 1e3, tol <- 1e-3, stepsize <- 1e-3), h0, omega1, omega2, eta, ksi){  ##梯度下降
  n <- dim(x)[1]
  p <- dim(x)[2]
  n.iter <- 1
  converge <- 0
  err.save <- c()
  obj.save <- c()

  while(n.iter <= control$max.iter){
    temp <- beta
    grad <- 0
    obj <- 0

    # calculate gradient
    for(i in 1:n){
      h0.fit <- h0(x[i,]%*%beta)
      f0 <- h0.fit$est
      d0 <- h0.fit$der
      grad <- grad + (f0 - y[i])*d0*x[i,]
      #obj = obj + (y[i] - f0)^2
    }
    grad <- grad / n - (omega1 + omega2) + alpha*(2*beta - eta - ksi)
    #obj  = obj / (2*n) + omega1 %*% (eta - b) + omega2 %*% (ksi - b) + alpha * (sum((eta - b)^2) + sum((ksi - b)^2)) / 2

    # update b
    beta <- beta - control$stepsize*grad

    # error
    #err = norm(beta - temp, '2')
    err <- max(abs(beta - temp))
    err.save <- c(err.save, err)
    #obj.save = c(obj.save, obj)

    error.control <- (err <= control$tol)
    if(error.control){
      converge <- 1
      break
    }

    n.iter <- n.iter + 1
  }

  return(list(beta <- beta, converge <- converge, iteration <- n.iter))#, err=err.save, obj=obj.save))
}
