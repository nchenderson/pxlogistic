gpxlogistic <- function(par, X, y, n.trials=rep(1, length(y)), lambda=NULL, 
                        weights=NULL, control=list()) {
  
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
    ans <- GPX_ECME_Logistic(beta.init, X, y, n.trials, method, svec=weights,
                            tol, maxiter, track.objective)
  } else {
    ans <- GPX_ECME_Pen_Logistic(beta.init, X, y, n.trials, method, svec=weights,
                                tol, maxiter, track.objective, lambda)
  }
  return(list(coef=ans$coef, iter=ans$iter, objfn.track=ans$objfn.track))
}




GPX_ECME_Logistic <- function(beta.init, X, y, n.trials, method, svec,
                              tol, maxiter, track.objective) {
  beta.old <- beta.init
  iter <- 0
  rho.old <- 1
  LogLikObserved <- function(rho, x.theta) {
    # plogis(x) = 1/(1 + exp(-x))
    -sum(rho*y*svec*x.theta + n.trials*svec*plogis(-rho*x.theta, log=TRUE) )
  }
  phi <- as.numeric(X%*%beta.old)
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 1)
    objfn.track[1] <- LogLikObserved(phi)
  }
  kappa.t <- 4/(max(ntrials)*((norm(XX, "2")^2)))
  kXtu <- kappa.t*crossprod(X*svec, y - n.trials/2)
  ww <- rep(0, length(y))
  while(TRUE){
    phat <- plogis(phi) # This is expit(phi)
    small.phi <- abs(phi) < 1e-4
    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
    ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
    step <- kXtu - kappa.t*as.numeric(crossprod(X*sqrt(ww*svec))%*%beta.old)
    theta <- beta.old + step

    x.theta <- as.numeric(X%*%theta)
    if(method=="Newton") {
      rho.new <- UpdateRhoNewton(rho.old, x.theta, y, n.trials, init=TRUE,
                                 tol.newton = .01*tol)
    } else if(method=="BFGS") {
      it <- ifelse(iter==0, 1e4, 1)
      alpha.new <- optim(alpha.old, fn=LogLikObserved, method="BFGS", 
                         control=list(maxit=it), x.theta=x.theta)$par
    }
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
    objfn.track <- NULL
  }
  return(list(coef=as.numeric(beta.new), iter=iter, objfn.track=objfn.track))
}


GPX_ECME_Pen_Logistic <- function(beta.init, X, y, n.trials, method, svec,
                              tol, maxiter, track.objective, lambda) {
  beta.old <- beta.init
  iter <- 0
  rho.old <- 1
  PenLogLikObserved <- function(rho, x.theta, theta) {
    # plogis(x) = 1/(1 + exp(-x))
    -sum(rho*y*svec*x.theta + n.trials*svec*plogis(-rho*x.theta, log=TRUE) ) + sum(lambda*rho*rho*theta*theta/2)
  }
  phi <- as.numeric(X%*%beta.old)
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 1)
    objfn.track[1] <- PenLogLikObserved(1, phi, beta.old)
  }
  kappa.t <- 4/(max(n.trials)*((norm(XX, "2")^2)))
  kXtu <- kappa.t*crossprod(X*svec, y - n.trials/2)
  ww <- rep(0, length(y))
  while(TRUE){
    phat <- plogis(phi) # This is expit(phi)
    small.phi <- abs(phi) < 1e-4
    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
    ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
    step <- kXtu - kappa.t*as.numeric(crossprod(X*sqrt(ww*svec))%*%beta.old)
    theta <- beta.old + step
    
    x.theta <- as.numeric(X%*%theta)
    if(method=="Newton") {
      rho.new <- UpdateRhoNewtonPen(rho.old, x.theta, theta, y, n.trials, init=TRUE,
                                    tol.newton = .01*tol, svec, lambda)
    } else if(method=="BFGS") {
      it <- ifelse(iter==0, 1e4, 1)
      rho.new <- optim(1, fn=LogLikObserved, method="BFGS", 
                         control=list(maxit=it), x.theta=x.theta, theta=theta)$par
    }
    beta.new <- as.numeric(rho.new*theta)
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter)  break
    
    beta.old <- beta.new
    rho.old <- rho.new
    
    phi <- as.numeric(X%*%beta.old)
    if(track.objective) {
      objfn.track[iter+1] <- PenLogLikObserved(1, phi, beta.old)
    }
  }
  if(track.objective) {
    objfn.track <- objfn.track[!is.na(objfn.track)]
  } else {
    objfn.track <- NULL
  }
  return(list(coef=as.numeric(beta.new), iter=iter, objfn.track=objfn.track))
}
