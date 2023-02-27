emlogistic <- function(par, X, y, n.trials=rep(1, length(y)), weights=NULL, 
                       lambda=NULL, intermed=FALSE, control=list()) {
  
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
  svec <- weights

  beta.init <- par
  beta.old <- beta.init
  uvec <- y - n.trials/2
  iter <- 0
  phi <- as.numeric(X%*%beta.old)
  ww <- rep(0, length(y))
  Xtu <- crossprod(X, svec*(y - n.trials/2))
  
  LogLikObserved <- function(x.theta) {
    -sum(y*svec*x.theta + n.trials*svec*plogis(-x.theta, log=TRUE) )
  }
  PenLogLikObserved <- function(x.theta, beta.vals) {
    -sum(y*svec*x.theta + n.trials*svec*plogis(-x.theta, log=TRUE) ) + (lambda/2)*sum(beta.vals*beta.vals)
  }
  BetaMat <- NULL
  if(intermed) {
    BetaMat <- matrix(0, nrow=maxiter + 1, ncol=length(beta.init))
    BetaMat[1,] <- beta.init
  }
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 1)
    if(is.null(lambda)) {
       objfn.track[1] <- LogLikObserved(phi)
    } else {
       objfn.track[1] <- PenLogLikObserved(phi, beta.old)
    }
  }
  while(TRUE) {
    ## Compute weight vector
    phat <- plogis(phi) # This is expit(phi)
    small.phi <- abs(phi) < 1e-4
    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
    ww[!small.phi] <- as.numeric((n.trials[!small.phi]*(phat[!small.phi] - 1/2))/phi[!small.phi]) 
    ww <- svec*ww
    
    if(is.null(lambda)) {
        beta.new <- as.numeric(solve(crossprod(X*sqrt(ww)), Xtu))
    } else {
        XtWX <- crossprod(X*sqrt(ww))
        diag(XtWX) <- diag(XtWX) + lambda
        beta.new <- as.numeric(solve(XtWX, Xtu))
    }
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter) break
    beta.old <- beta.new
    phi <- as.numeric(X%*%beta.old)
    if(track.objective) {
      if(is.null(lambda)) {
         objfn.track[iter+1] <- LogLikObserved(phi)
      } else {
         objfn.track[iter+1] <- PenLogLikObserved(phi, beta.old)
      }
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
