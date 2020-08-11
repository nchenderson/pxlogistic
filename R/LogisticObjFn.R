LogisticObjFn <- function(par, X, y) {
  
  ## Use Log-Sum-Exp trick for more stable computation of the log-likelihood
  x.beta <- as.numeric(X%*%par)
  ans <- sum(y*x.beta) - sum(log1p(exp(x.beta[x.beta < 0]))) - sum(x.beta[x.beta >= 0]) - sum(log1p(exp(-x.beta[x.beta >= 0])))
  return(ans)
}