LogisticObjFn <- function(par, X, y, n.trials=rep(1, length(y))) {
  
  ## Use Log-Sum-Exp trick for more stable computation of the log-likelihood
  x.beta <- as.numeric(X%*%par)
  lt.zero <- x.beta < 0
  ans <- sum(y*x.beta) - sum(n.trials[lt.zero]*log1p(exp(x.beta[lt.zero]))) - sum(n.trials[!lt.zero]*x.beta[!lt.zero]) - sum(n.trials[!lt.zero]*log1p(exp(-x.beta[!lt.zero])))
  return(ans)
}