-----

<!-- README.md is generated from README.Rmd. Please edit that file -->

# mSIM

<!-- badges: start -->

<!-- badges: end -->

This package implements the Multivariate Single index Model(`mSIM`)
proposed by Feng et al. It can be formulated as following:

\(\mathbf{y}_{i}=\left\{f_{1}\left(\mathbf{x}_{i}^{\top} \mathbf{B}_{1}\right), \cdots, f_{q}\left(\mathbf{x}_{i}^{\top} \mathbf{B}_{\cdot q}\right)\right\}^{\top}+\boldsymbol{\epsilon}_{i}\)

The mSIM algorithm can both tackle potential nonlinearity in
multivariate response regression and challenges with high dimensional
data. For more details, please read Feng’s paper *Sparse Single Index
Models for Multivariate Responses*.

The package provides simulation data generation, mSIM model fitting,
Bayesian Information Criteria evaluation and prediction to solve
Multivariate Single Index Model problems.

## Installation

The current working version can be installed from Github:

``` r
library(devtools)
install_github("ecchi/mSIM")
```

## Example

We use a small dataset within our package to illustrate our package.

``` r
library(mSIM)
#> Warning: replacing previous import 'glmnet::na.replace' by 'gtools::na.replace'
#> when loading 'mSIM'
#> Warning: replacing previous import 'gtools::scat' by 'mgcv::scat' when loading
#> 'mSIM'
## basic example code

## fitting the mSIM model
X_train = scale(X_train)
Y_train = scale(Y_train)
result = get_B_ADMM(Y = Y_train, X = X_train, lambda=0.03877, rank=3)

## Predict index covariate matrix
test_result = model_pred(Y = Y_train, X = X_train, B = result$B.sparse, Y.true = Y_test, X.pred = X_test)
```

if you want to run multiple (lambda, rank) pair and use Bayesian
Information criteria to evaluate them, the following code is an example.

``` r
## Evaluating by using Bayesian Information criteria
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
```

## Authors

  - [Zirui Li](https://github.com/GryffindorLi)
    
    School of Electrical and Information Engineering, Tianjin University

  - \[Yuan Feng\]

  - [Luo Xiao](https://www4.stat.ncsu.edu/~xiao/)
    
    Department of Statistics, North Carolina State University

  - [Eric C. CHI](http://www.ericchi.com/)
    
    Department of Statistics, North Carolina State University

## References

Yuan Feng, Luo Xiao, and Eric C. Chi, Sparse Single Index Models for
Multivariate Responses, Journal of Computational and Graphical
Statistics, 2020.
