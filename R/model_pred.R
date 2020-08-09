#' The predict function for the mSIM package
#'
#' @export model_pred
#'
#' @param Y The training index matrix
#' @param X The training covariate matrix
#' @param B The coefficient matrix
#' @param Y.true The groundtruth for predicting index matrix
#' @param X.pred The covariate matrix for prediction
#' @return A list that contains the prediction result of Y and the MSE loss
#' @examples
#' pred = predict(Y = Y, X = X, real.MSIM$B.sparse, Y.true = Y.true, X.pred = X.pred)

model_pred = function(Y, X, B, Y.true, X.pred){
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
