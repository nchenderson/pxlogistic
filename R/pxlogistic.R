pxlogistic <- function(par, X, y, n.trials=rep(1, length(y)), lambda=NULL, track.objective=FALSE,
                       control=list()) {
  
  control.default <- list(maxiter=2000, tol=1e-7, method="Newton")
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)
  
  maxiter <- control$maxiter
  tol <- control$tol
  method <- control$method 
  objfn.track <- NULL
  
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
  phi <- as.numeric(X%*%beta.old)
  if(track.objective) {
      objfn.track = rep(NA, maxiter + 1)
      objfn.track[1] <- LogLikObserved(1, phi)
  }
   
  
  while(TRUE) {
      ## Compute weight vector
      phat <- plogis(phi)
      small.phi <- abs(phi) < 1e-4
      ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
      ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
      
      theta <- solve(crossprod(X, X*ww), crossprod(X, uvec))/alpha.old
      x.theta <- as.numeric(X%*%theta)
      if(method=="Newton") {
          alpha.new <- UpdateAlphaNewton(alpha.old, x.theta, y, n.trials, init=TRUE)
      } else if(method=="BFGS") {
         it <- ifelse(iter==0, 1e4, 1)
         alpha.new <- optim(alpha.old, fn=LogLikObserved, method="BFGS", 
                            control=list(maxit=it), x.theta=x.theta)$par
      }
      ## why no uniroot or optimize as an option?
      beta.new <- as.numeric(alpha.new*theta)
    
      iter <- iter + 1
      if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter) break
      beta.old <- beta.new
      alpha.old <- alpha.new
      
      phi <- as.numeric(X%*%beta.old)
      if(track.objective) {
          objfn.track[iter+1] <- LogLikObserved(1, phi)
      }
  }
  if(track.objective) {
      objfn.track <- objfn.track[!is.na(objfn.track)]
  }
  return(list(coef=beta.new, iter=iter, objfn.track=objfn.track))
}


#UpdateAlphaNewton <- function(alpha.old, x.theta, y, n, init=FALSE){
#  if(init) alpha.old <- 0
  
#  alpha.new <- alpha.old
#  while(TRUE){
#    ept <- expit(alpha.old*x.theta)
#    dept <- ept*(1-ept)
#    f.st <- sum(x.theta*(y - n*ept))
#    if(abs(f.st) < 1e-8) break
#    f.nd <- sum(x.theta^2*n*dept)
#    alpha.new <- as.numeric(alpha.old + f.st/f.nd)
#    alpha.old <- alpha.new
#  }
#  ans <- alpha.new
#  return(ans)
#}
