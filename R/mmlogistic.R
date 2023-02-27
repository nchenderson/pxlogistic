pxmmlogistic <- function(par, X, y, n.trials=rep(1, length(y)), weights=NULL,
                       lambda=NULL, control=list()) {
  
  control.default <- list(maxiter=2000, tol=1e-7, method="Newton", objfn.track=TRUE)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)
  
  maxiter <- control$maxiter
  tol <- control$tol
  method <- control$method 
  track.objective <- control$objfn.track
  beta.init <- par
  
  if(is.null(weights)) {
    weights <- rep(1, length(y))
  }
  if(is.null(lambda)) {
      ans <- PXMM_Logistic(beta.init, X, y, n.trials, method, 
                            tol, maxiter, track.objective)
  } else {
      ans <- PXMM_Pen_Logistic(beta.init, X, y, n.trials, method, svec=weights,
                                tol, maxiter, track.objective, lambda)
  }
  return(list(coef=ans$coef, iter=ans$iter, objfn.track=ans$objfn.track))
}

PXMM_Logistic <- function(beta.init, X, y, n.trials, method, 
                             tol, maxiter, track.objective) {
  
  beta.old <- beta.init
  iter <- 0
  rho.old <- 1
  LogLikObserved <- function(rho, x.theta) {
    -sum(rho*y*x.theta - n.trials*log1p(exp(rho*x.theta)) )
  }
  ww <- rep(0, length(y))
  phi <- as.numeric(X%*%beta.old)
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 1)
    objfn.track[1] <- LogLikObserved(1, phi)
  }
  
  XtX.inv <- solve(crossprod(X))
  for(k in 1:maxiter) {
    ## Compute weight vector
    phat <- n.trials*plogis(phi) # This is expit(phi)
    small.phi <- abs(phi) < 1e-4
    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
    ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
    kappa.t <- max(ww)
    
    Xtr <- crossprod(X, (y - phat)/kappa.t)
    theta <- beta.old + XtX.inv%*%Xtr
    x.theta <- as.numeric(X%*%theta)
    if(method=="Newton") {
      rho.new <- UpdateRhoNewton(rho.old, x.theta, y, n.trials, init=TRUE,
                                 tol.newton = .01*tol)
    } else if(method=="BFGS") {
      it <- ifelse(iter==0, 1e4, 1)
      alpha.new <- optim(alpha.old, fn=LogLikObserved, method="BFGS", 
                         control=list(maxit=it), x.theta=x.theta)$par
    }
    ## why no uniroot or optimize as an option?
    beta.new <- as.numeric(rho.new*theta)
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter)  break
    
    beta.old <- beta.new
    rho.old <- rho.new
    
    phi <- as.numeric(X%*%beta.old)
    if(track.objective) {
      objfn.track[iter+1] <- LogLikObserved(1, phi)
    }
  }
  if(track.objective) {
    objfn.track <- objfn.track[!is.na(objfn.track)]
  } else {
    objfn.track <- FALSE
  }
  return(list(coef=beta.new, iter=iter, objfn.track=objfn.track))
}

PXMM_Pen_Logistic <- function(beta.init, X, y, n.trials, method, svec,
                              tol, maxiter, track.objective, lambda) {
  
  beta.old <- beta.init
  iter <- 0
  rho.old <- 1
  PenLogLikObserved <- function(rho, x.theta, theta, lam) {
    # plogis(x) = 1/(1 + exp(-x))
    -sum(rho*y*svec*x.theta + n.trials*svec*plogis(-rho*x.theta, log=TRUE) ) + sum(lam*rho*rho*theta*theta/2)
  }
  ww <- rep(0, length(y))
  phi <- as.numeric(X%*%beta.old)
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 1)
    objfn.track[1] <- PenLogLikObserved(1, phi, beta.old, lambda)
  }
  
  kappa.t <- max(n.trials/4)
  XtSX <- crossprod(X*sqrt(svec))
  diag(XtSX) <- diag(XtSX) + lambda/kappa.t
  XtX.inv.XtS <- solve(XtSX, t(X*svec))
  for(k in 1:maxiter) {
    ## Compute weight vector
    phat <- plogis(phi) # This is expit(phi)
    
    theta <- XtX.inv.XtS%*%(phi + (y - phat)/kappa.t) 
    x.theta <- as.numeric(X%*%theta)
    if(method=="Newton") {
      rho.new <- UpdateRhoNewtonPen(rho.old, x.theta, theta, y, n.trials, init=TRUE,
                                    tol.newton = .01*tol, svec, lambda)
    } else if(method=="BFGS") {
      it <- ifelse(iter==0, 1e4, 1)
      alpha.new <- optim(alpha.old, fn=LogLikObserved, method="BFGS", 
                         control=list(maxit=it), x.theta=x.theta)$par
    }
    ## why no uniroot or optimize as an option?
    beta.new <- as.numeric(rho.new*theta)
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter)  break
    
    beta.old <- beta.new
    rho.old <- rho.new
    
    phi <- as.numeric(X%*%beta.old)
    if(track.objective) {
      objfn.track[iter+1] <- PenLogLikObserved(1, phi, beta.old, lambda)
    }
  }
  if(track.objective) {
    objfn.track <- objfn.track[!is.na(objfn.track)]
  } else {
    objfn.track <- NULL
  }
  return(list(coef=beta.new, iter=iter, objfn.track=objfn.track))
}

