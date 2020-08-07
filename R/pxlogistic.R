pxlogistic <- function(par, X, y, n.trials=rep(1, length(y)), lambda=NULL, control=list()) {
  
  control.default <- list(maxiter=2000, tol=1e-7, method="Newton")
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)
  
  maxiter <- control$maxiter
  tol <- control$tol
  method <- control$method 
  
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
  while(TRUE) {
      phi <- as.numeric(X%*%beta.old)
      ww <- as.numeric((n.trials/(2*phi))*tanh(phi/2)) ## Need to modify weight updates so no error is thrown when phi=0
      theta <- solve(crossprod(X, X*ww), crossprod(X, uvec))/alpha.old
      ## maybe add the MO in another function
      x.theta <- as.numeric(X%*%theta)
      if(method=="Newton") {
          alpha.new <- UpdateAlphaNewton(alpha.old, x.theta, y, n.trials, init=!iter)
      } else if(method=="BFGS") {
         it <- ifelse(iter==0, 1e4, 1)
         alpha.new <- optim(alpha.old, fn=LogLikObserved, method="BFGS", 
                            control=list(maxit=it), x.theta=x.theta)$par
      }
      ## why no uniroot or optimize as an option?
      beta.new <- alpha.new*theta
    
      iter <- iter + 1
      if(norm(beta.new-beta.old, "2") < tol | iter > maxiter) break
      beta.old <- beta.new
      alpha.old <- alpha.new
  }
  return(list(coef=beta.new, iter=iter))
}
