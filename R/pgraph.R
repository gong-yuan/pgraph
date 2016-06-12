#' Compute the Projected Graph
#'
#' \code{pgraph} calculate the projected graph
#'
#' @export
# @importFrom parallel mcmapply
#' @param z n * p dimensional matrix
#' @param f n * q factor matrix. Default = 'NULL'.
#' @param method projection method. Default = 'linear'.
#' @param cond whether to create a conditional graph or unconditional graph.
#' Default = TRUE. If cond = FALSE, f must be provided.
#' @param trace whether to trace to estimation process.
#' @return a list.
#' \item{test}{distance covariance test object}
#' \item{xeps}{residual of the projection of x}
#' \item{yeps}{residual of the projection of y}
#' @examples
#'
#' p=30;
#' n=120;
#' tmp=runif(p-1,0.5,1)
#' s=c(0,cumsum(tmp));
#' s1=matrix(s,p,p)
#' cov.mat.true=exp(-abs(s1-t(s1)))
#' prec.mat.true=solve(cov.mat.true);
#' a=matrix(rnorm(p*n),n,p)
#' data.sa=a%*%chol(cov.mat.true);
#' fit = pgraph(data.sa)
#' fit2 = pgraph(data.sa, cor = 'dcov')
pgraph<-function(z, f = NULL, method = c("ridge","lasso","sam"), cond = TRUE , trace = FALSE){
 method = match.arg(method)
 n = nrow(z)
 p = ncol(z)
 statmat.pearson = matrix(0,p,p)
 statmat.dcov = matrix(0,p,p)
 for(i in 1:(p-1))
  for(j in (i+1):p){
    if(trace){
      cat('i=',i,'j=',j,'\n')
    }

x = z[,i]
y = z[,j]
if(cond) {
  b = cbind(z[,-c(i,j)],f)
} else {
  if(is.null(f)) {
    stop('f must be provided when cond is FALSE.')
  }
  b = f
}
fit = projcov(x,y,b,method)
statmat.pearson[i,j] = fit$test.pearson
statmat.dcov[i,j] = fit$test.dcov
  }
 statmat.pearson = statmat.pearson+t(statmat.pearson)
 statmat.dcov = statmat.dcov+t(statmat.dcov)
return(list(statmat.pearson = statmat.pearson,statmat.dcov = statmat.dcov))

}
