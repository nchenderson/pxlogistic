gpxlogistic <- function(par, X, y, n.trials=rep(1, length(y)), lambda=NULL, control=list()) {
  
  control.default <- list(maxiter=2000, tol=1e-7, method="Newton")
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)
  
  maxiter <- control$maxiter
  tol <- control$tol
  method <- control$method 
  
  gam <- 4/max(eigen(crossprod(X))$values)  ## need to change this
  alpha.old <- 1
  beta.init <- par
  beta.old <- beta.init
  uvec <- y - n.trials/2
  iter <- 0
  alpha.old <- 1
  LogLikObserved <- function(alp, x.theta) {
    -sum(alp*y*x.theta - n.trials*log1p(exp(alp*x.theta)) )
  }
  PenLogLikObserved <- function(alp, x.theta, theta, lam) {
    -sum(alp*y*x.theta - n.trials*log1p(exp(alp*x.theta)) - lam*alp*alp*theta*theta/2)
  }
  ww <- rep(0, length(y))
  while(TRUE){
     phi <- as.numeric(X%*%beta.old)
     small.phi <- abs(phi) < 1e-4
     ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
     ww[!small.phi] <- as.numeric((n.trials[!small.phi]*tanh(phi[!small.phi]/2))/(2*phi[!small.phi])) 
    
     theta.new <- theta.old + gam/alpha.old*(xtk - crossprod(X*ww, phi))
     x.theta <- as.numeric(X%*%theta.new)
     alpha.new <- UpdateAlphaNewton(alpha.old, x.theta, y, n.trials, init=TRUE)
     beta.new <- theta.new*alpha.new
     
     iter <- iter + 1
     if(norm(beta.new-beta.old, "2") < tol | iter >= max.iter) break
     beta.old <- beta.new
     alpha.old <- alpha.new
     theta.old <- theta.new
  }
  return(list(beta=beta.new, iter=iter))
}


