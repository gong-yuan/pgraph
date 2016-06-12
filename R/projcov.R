#' Calculate the Projected covariance of two vectors
#'
#' \code{projcov} calculate the projected distance covariance of two vectors given
#' factors.
#'
#' @export
#' @importFrom energy dcov.test
#' @importFrom energy dcor
#' @importFrom glmnet cv.glmnet
#' @importFrom splines ns
#' @importFrom SAM samQL
#' @importFrom parcor ridge.cv
#' @param x first vector
#' @param y second vector
#' @param b factor matrix
#' @param method projection method. Default = 'linear'.
#' @return a list.
#' \item{test}{distance covariance test object}
#' \item{xeps}{residual of the projection of x}
#' \item{yeps}{residual of the projection of y}
#' @examples
#' K = 3
#' n = 100
#' b = matrix(rnorm(K*n),n,K)
#' bx = 1:3
#' by = c(1,2,2)
#' x = b%*%bx+rnorm(n)
#' y = b%*%by+rnorm(n)
#' fit1 = projcov(x, y, b)
#' fit2 = projcov(x, y, b, method = "sam")
projcov <- function(x, y, b, method = c("ridge","lasso","sam")){
  method = match.arg(method)
 # cor = match.arg(cor)
  if(method == 'ridge'){
    xfit = ridge.cv(b,x)
    xeps = x - b%*%xfit$coefficients-xfit$intercept
    yfit = ridge.cv(b,y)
    yeps = y - b%*%yfit$coefficients-yfit$intercept
}else if(method == 'lasso'){
xfit = cv.glmnet(b,x)
xeps = x - predict(xfit, b)
yfit = cv.glmnet(b,y)
yeps = y - predict(yfit, b)
} else if(method == 'sam'){
  xfit = cv.samQL(b,x)
  xeps = x - predict(xfit$sam.fit, b)$values[,xfit$lambda.min]
  yfit = cv.samQL(b,y)
  yeps = y - predict(yfit$sam.fit, b)$values[,yfit$lambda.min]
}

  test.pearson = abs(cor(xeps,yeps))
  test.dcov = dcor(xeps,yeps)
#  test.dcov = dcov.test(xeps,yeps)$estimate

return(list=list(test.pearson = test.pearson, test.dcov = test.dcov, xeps=xeps,yeps=yeps))
}
