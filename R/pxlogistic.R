pxlogistic <- function(par, X, y, n.trials=rep(1, length(y)), weights=NULL,
                       lambda=NULL, intermed=FALSE, aa.accelerate=FALSE, 
                       control=list()) {
  
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
  #beta.init <- par
  #beta.old <- beta.init
  #uvec <- y - n.trials/2
  #iter <- 0
  #alpha.old <- 1
  #LogLikObserved <- function(alp, x.theta) {
  #    -sum(alp*y*x.theta - n.trials*log1p(exp(alp*x.theta)) )
  #}
  #PenLogLikObserved <- function(alp, x.theta, theta, lam) {
  #    -sum(alp*y*x.theta - n.trials*log1p(exp(alp*x.theta)) - lam*alp*alp*theta*theta/2)
  #}
  #ww <- rep(0, length(y))
  #phi <- as.numeric(X%*%beta.old)
  #if(track.objective) {
  #    objfn.track = rep(NA, maxiter + 1)
  #    objfn.track[1] <- LogLikObserved(1, phi)
  #}
   
  #for(k in 1:maxiter) {
      ## Compute weight vector
   #   phat <- plogis(phi) # This is expit(phi)
  #    small.phi <- abs(phi) < 1e-4
  #    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
   #   ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
      
    #  theta <- solve(crossprod(X*sqrt(ww)), crossprod(X, uvec))/alpha.old
    #  x.theta <- as.numeric(X%*%theta)
    #  if(method=="Newton") {
    #      alpha.new <- UpdateAlphaNewton(alpha.old, x.theta, y, n.trials, init=TRUE,
    #                                     tol.newton = .01*tol)
    #  } else if(method=="BFGS") {
    #     it <- ifelse(iter==0, 1e4, 1)
    #     alpha.new <- optim(alpha.old, fn=LogLikObserved, method="BFGS", 
    #                        control=list(maxit=it), x.theta=x.theta)$par
     # }
      ## why no uniroot or optimize as an option?
    #  beta.new <- as.numeric(alpha.new*theta)
    
    #   iter <- iter + 1
    #  if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter)  break

     # beta.old <- beta.new
    #  alpha.old <- alpha.new
      
    #  phi <- as.numeric(X%*%beta.old)
    #  if(track.objective) {
    #      objfn.track[iter+1] <- LogLikObserved(1, phi)
    #  }
  #}
  #if(track.objective) {
  #    objfn.track <- objfn.track[!is.na(objfn.track)]
  #}
   if(is.null(lambda) & !aa.accelerate) {
       ans <- PX_ECME_Logistic(beta.init, X, y, n.trials, method, svec=weights, 
                                tol, maxiter, track.objective, intermed)
   } else if(!is.null(lambda) & !aa.accelerate) {
       ans <- PX_ECME_Pen_Logistic(beta.init, X, y, n.trials, method, svec=weights,
                                   tol, maxiter, track.objective, lambda)
   } else if(is.null(lambda) & aa.accelerate) {
       ans <- PX_AA_Logistic(beta.init, X, y, n.trials, method, svec=weights, 
                             tol, maxiter, track.objective, intermed)
   } else if(!is.null(lambda) & aa.accelerate) {
       ans <- PX_AA_Pen_Logistic(beta.init, X, y, n.trials, method, svec=weights, 
                             tol, maxiter, track.objective, intermed, lambda)
   }
   return(ans)
}

### Put order 1 AA Logistic regression here 

PX_ECME_Logistic <- function(beta.init, X, y, n.trials, method, svec,
                             tol, maxiter, track.objective, intermed) {
  # svec is the vector of weights
  beta.old <- beta.init
  Xtu <- crossprod(X, svec*(y - n.trials/2))
  iter <- 0
  rho.old <- 1
  LogLikObserved <- function(rho, x.theta) {
    # plogis(x) = 1/(1 + exp(-x))
     -sum(rho*y*svec*x.theta + n.trials*svec*plogis(-rho*x.theta, log=TRUE) )
  }
  ScoreObserved <- function(rho, x.theta) {
     -sum(svec*x.theta*(y - n.trials*plogis(rho*x.theta) ))
  }
  ww <- rep(0, length(y))
  phi <- as.numeric(X%*%beta.old)
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 1)
    objfn.track[1] <- LogLikObserved(1, phi)
  }
  BetaMat <- NULL
  if(intermed) {
    BetaMat <- matrix(0, nrow=maxiter + 1, ncol=length(beta.init))
    BetaMat[1,] <- beta.init
  }
  for(k in 1:maxiter) {
     ## Compute weight vector
     phat <- plogis(phi) # This is expit(phi)
     small.phi <- abs(phi) < 1e-4
     ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
     ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
     ww <- svec*ww
     
     theta <- solve(crossprod(X*sqrt(ww)), Xtu)
     x.theta <- as.numeric(X%*%theta)
     if(method=="Newton") {
        rho.new <- UpdateRhoNewton(1, x.theta, y, n.trials, init=TRUE,
                                   tol.newton = .01*tol, svec)
       # rho.low <- -6/(min(abs(x.theta)) + 0.01)
       # rho.high <- 6/(min(abs(x.theta)) + 0.01)
       # if(k==1) {
       #     rho.new <- 1
       # } else {
        #    rho.new <- uniroot(ScoreObserved, interval=c(rho.low, rho.high), x.theta=x.theta)$root
       # }
        #rho.new <- 1
     } else if(method=="BFGS") {
         it <- ifelse(iter==0, 1e4, 1)
         rho.new <- optim(1, fn=LogLikObserved, method="BFGS", 
                           control=list(maxit=it), x.theta=x.theta)$par
     } else if(method=="Brent") {
         rho.low <- -6/(min(abs(x.theta)) + 0.01)
         rho.high <- 6/(min(abs(x.theta)) + 0.01)
         rho.new <- uniroot(ScoreObserved, interval=c(rho.low, rho.high), x.theta=x.theta,
                            tol=1e-6)$root
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
     if(intermed) {
       BetaMat[iter+1,] <- beta.new
     }
  }
  if(track.objective) {
    objfn.track <- objfn.track[!is.na(objfn.track)]
  } else {
    objfn.track <- NULL
  }
  if(intermed) {
     BetaMat <- BetaMat[!is.na(objfn.track),]
  }
  return(list(coef=beta.new, iter=iter, objfn.track=objfn.track, intermed=BetaMat))
}

PX_ECME_Pen_Logistic <- function(beta.init, X, y, n.trials, method, svec,
                                 tol, maxiter, track.objective, lambda) {

  beta.old <- beta.init
  Xtu <- crossprod(X, svec*(y - n.trials/2))
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
  
  for(k in 1:maxiter) {
    ## Compute weight vector
    phat <- plogis(phi) # This is expit(phi)
    small.phi <- abs(phi) < 1e-4
    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
    ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
    ww <- svec*ww
    
    XtWX <- crossprod(X*sqrt(ww))
    diag(XtWX) <- diag(XtWX) + lambda
    theta <- solve(XtWX, Xtu)
    x.theta <- as.numeric(X%*%theta)
    if(method=="Newton") {
      rho.new <- UpdateRhoNewtonPen(rho.old, x.theta, theta, y, n.trials, init=TRUE,
                                   tol.newton = .01*tol, svec, lambda)
    } else if(method=="BFGS") {
      it <- ifelse(iter==0, 1e4, 1)
      rho.new <- optim(rho.old, fn=PenLogLikObserved, method="BFGS", 
                         control=list(maxit=it), x.theta=x.theta, theta=theta)$par
    }
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


PX_AA_Logistic <- function(beta.init, X, y, n.trials, method, svec,
                           tol, maxiter, track.objective, intermed) {
  # svec is the vector of weights
  beta.old <- beta.init
  Xtu <- crossprod(X, svec*(y - n.trials/2))
  iter <- 0
  rho.old <- 1
  LogLikObserved <- function(rho, x.theta) {
    # plogis(x) = 1/(1 + exp(-x))
    -sum(rho*y*svec*x.theta + n.trials*svec*plogis(-rho*x.theta, log=TRUE) )
  }
  ww <- rep(0, length(y))
  phi <- as.numeric(X%*%beta.old)
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 2)
    objfn.track[1] <- LogLikObserved(1, phi)
  }
  BetaMat <- NULL
  
  beta.lag <- beta.init  ## beta.lag is beta^(0)
  phat <- plogis(phi) # This is expit(phi)
  small.phi <- abs(phi) < 1e-4
  ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
  ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
  ww <- svec*ww
  
  theta.old <- solve(crossprod(X*sqrt(ww)), Xtu) ## theta.old is beta^(1), EM
  phi <- as.numeric(X%*%theta.old)
  loglik.old <- LogLikObserved(1, phi)
  if(intermed) {
    BetaMat <- matrix(0, nrow=maxiter + 2, ncol=length(beta.init))
    BetaMat[1,] <- beta.init
    BetaMat[2,] <- theta.old
  }
  if(track.objective) {
    objfn.track[2] <- loglik.old
  }
  for(k in 1:maxiter) {
    ## Compute weight vector
    phat <- plogis(phi) # This is expit(phi)
    small.phi <- abs(phi) < 1e-4
    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
    ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
    ww <- svec*ww
    
    theta.new <- solve(crossprod(X*sqrt(ww)), Xtu)
    
    rr <- theta.new - beta.old
    vv <- rr + beta.lag - theta.old
    gamma.new <- sum(rr*vv)/sum(vv*vv)
    ## just need to be sure theta.old is defined correctly.
    beta.prop <- as.numeric((1 - gamma.new)*theta.new) + as.numeric(gamma.new*theta.old)
    loglik.prop <- LogLikObserved(1, as.numeric(X%*%beta.prop))
    ## check monotonicity
    if(loglik.prop <= loglik.old) {
      beta.new <- beta.prop
      loglik.old <- loglik.prop
    } else {
      beta.new <- theta.new
      loglik.old <- LogLikObserved(1, as.numeric(X%*%beta.new))
    }
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter)  break
    
    beta.lag <- beta.old
    beta.old <- beta.new
    theta.old <- theta.new
    
    phi <- as.numeric(X%*%beta.old)
    if(track.objective) {
      objfn.track[iter+2] <- LogLikObserved(1, phi)
    }
    if(intermed) {
      BetaMat[iter+2,] <- beta.new
    }
  }
  if(track.objective) {
    objfn.track <- objfn.track[!is.na(objfn.track)]
  } else {
    objfn.track <- NULL
  }
  if(intermed) {
    BetaMat <- BetaMat[!is.na(objfn.track),]
  }
  return(list(coef=beta.new, iter=iter, objfn.track=objfn.track, intermed=BetaMat))
}



PX_AA_Pen_Logistic <- function(beta.init, X, y, n.trials, method, svec,
                           tol, maxiter, track.objective, intermed, lambda) {
  # svec is the vector of weights
  beta.old <- beta.init
  Xtu <- crossprod(X, svec*(y - n.trials/2))
  iter <- 0
  rho.old <- 1
  PenLogLikObserved <- function(rho, x.theta, theta) {
    # plogis(x) = 1/(1 + exp(-x))
    -sum(rho*y*svec*x.theta + n.trials*svec*plogis(-rho*x.theta, log=TRUE) ) + sum(lambda*rho*rho*theta*theta/2)
  }
  ww <- rep(0, length(y))
  phi <- as.numeric(X%*%beta.old)
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 2)
    objfn.track[1] <- PenLogLikObserved(1, phi, beta.old)
  }
  BetaMat <- NULL
  
  beta.lag <- beta.init  ## beta.lag is beta^(0)
  phat <- plogis(phi) # This is expit(phi)
  small.phi <- abs(phi) < 1e-4
  ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
  ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
  ww <- svec*ww
  
  XtWX <- crossprod(X*sqrt(ww))
  diag(XtWX) <- diag(XtWX) + lambda
  theta.old <- solve(XtWX, Xtu) ## theta.old is beta^(1), EM
  phi <- as.numeric(X%*%theta.old)
  loglik.old <- PenLogLikObserved(1, phi, theta.old)
  if(intermed) {
    BetaMat <- matrix(0, nrow=maxiter + 2, ncol=length(beta.init))
    BetaMat[1,] <- beta.init
    BetaMat[2,] <- theta.old
  }
  if(track.objective) {
    objfn.track[2] <- loglik.old
  }
  for(k in 1:maxiter) {
    ## Compute weight vector
    phat <- plogis(phi) # This is expit(phi)
    small.phi <- abs(phi) < 1e-4
    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
    ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
    ww <- svec*ww
    
    XtWX <- crossprod(X*sqrt(ww))
    diag(XtWX) <- diag(XtWX) + lambda
    theta.new <- solve(XtWX, Xtu) 
    
    rr <- theta.new - beta.old
    vv <- rr + beta.lag - theta.old
    gamma.new <- sum(rr*vv)/sum(vv*vv)
    ## just need to be sure theta.old is defined correctly.
    beta.prop <- as.numeric((1 - gamma.new)*theta.new) + as.numeric(gamma.new*theta.old)
    loglik.prop <- PenLogLikObserved(1, as.numeric(X%*%beta.prop), beta.prop)
    ## check monotonicity
    if(loglik.prop <= loglik.old) {
      beta.new <- beta.prop
      loglik.old <- loglik.prop
    } else {
      beta.new <- theta.new
      loglik.old <- PenLogLikObserved(1, as.numeric(X%*%beta.new), beta.new)
    }
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter)  break
    
    beta.lag <- beta.old
    beta.old <- beta.new
    theta.old <- theta.new
    
    phi <- as.numeric(X%*%beta.old)
    if(track.objective) {
      objfn.track[iter+2] <- PenLogLikObserved(1, phi, beta.old)
    }
    if(intermed) {
      BetaMat[iter+2,] <- beta.new
    }
  }
  if(track.objective) {
    objfn.track <- objfn.track[!is.na(objfn.track)]
  } else {
    objfn.track <- NULL
  }
  if(intermed) {
    BetaMat <- BetaMat[!is.na(objfn.track),]
  }
  return(list(coef=beta.new, iter=iter, objfn.track=objfn.track, intermed=BetaMat))
}

