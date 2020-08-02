## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mSIM)

## -----------------------------------------------------------------------------
X_train = scale(X_train)
Y_train = scale(Y_train)
result = get_B_ADMM(Y = Y_train, X = X_train, lambda=0.03877, rank=3)

## -----------------------------------------------------------------------------
# set up tuning parameter
range = c(0.25, 0.75)
rank = c(3, 10)
tuning = list() # all combination of tuning parameters
len.tuning = 0
for(i in rank){
  for(j in range){
    len.tuning = len.tuning + 1
    tuning[[len.tuning]] = c(j, i)
  }
} 

# mSIM model.
temp = list()
for(i in 1:length(tuning)){
    temp[[i]] = get_B_ADMM(Y=Y_train, X=X_train, lambda=tuning[[i]][1],
rank=tuning[[i]][2], alpha=1, control1=list(max.iter=5e1, tol=1e-3), select.method='linear', plot=F)
}
# B for each tuning parameter.
B = list()
for(i in 1:length(tuning)){
    B[[i]] = temp[[i]]$B.sparse
}

BIC.msim = NULL
# BIC for each tuning parameter.
for(i in 1:length(tuning)){
    BIC.msim = rbind(BIC.msim, B_BIC(Y=Y_train, X=X_train, B=B[[i]], tuning=tuning[[i]], linear=F))
}

## -----------------------------------------------------------------------------
X_train = scale(X_train)
Y_train = scale(Y_train)
result = get_B_ADMM(Y = Y_train, X = X_train, lambda=0.03877, rank=3)
test_result = model_pred(Y = Y_train, X = X_train, B = result$B.sparse, Y.true = Y_test, X.pred = X_test)

