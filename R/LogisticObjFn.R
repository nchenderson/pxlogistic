LogisticObjFn <- function(par, X, y, n.trials=rep(1, length(y)), weights=NULL) {
  
  ## Use Log-Sum-Exp trick for more stable computation of the log-likelihood
  x.beta <- as.numeric(X%*%par)
  if(is.null(weights)) {
    weights <- rep(1, length(y))
  }
  ans <- sum(y*weights*x.beta + n.trials*weights*plogis(-x.beta, log=TRUE) )
  return(ans)
}

L2PenLogisticObjFn <- function(par, X, y, n.trials=rep(1, length(y)), weights=NULL,
                               lambda=0) {
  
  ## Use Log-Sum-Exp trick for more stable computation of the log-likelihood
  x.beta <- as.numeric(X%*%par)
  if(is.null(weights)) {
    weights <- rep(1, length(y))
  }
  ans <- sum(y*weights*x.beta + n.trials*weights*plogis(-x.beta, log=TRUE) ) - (lambda/2)*sum(par*par)
  return(ans)
}

L1PenLogisticObjFn <- function(par, X, y, n.trials=rep(1, length(y)), weights=NULL,
                               lambda=0) {
  
  ## Use Log-Sum-Exp trick for more stable computation of the log-likelihood
  x.beta <- as.numeric(X%*%par)
  if(is.null(weights)) {
    weights <- rep(1, length(y))
  }
  ans <- sum(y*weights*x.beta + n.trials*weights*plogis(-x.beta, log=TRUE) ) - lambda*sum(abs(par))
  return(ans)
}


