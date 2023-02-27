pxgdlogisticlasso <- function(par, X, y, n.trials=rep(1, length(y)), weights=NULL, 
                              lambda=NULL, intermed=FALSE, control=list()) {
  ## This is like the "GEM" algorithm
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
  stplngth <- (2/norm(X, type="2"))^2
  
  beta.init <- par
  beta.old <- beta.init
  iter <- 0
  phi <- as.numeric(X%*%beta.old)
  ww <- rep(0, length(y))
  XSty <- crossprod(X, svec*y)
  LogLikObserved <- function(x.theta, bet) {
    sum(y*svec*x.theta + n.trials*svec*plogis(-x.theta, log=TRUE)) - lambda*sum(abs(bet))
  }
  ScoreObserved <- function(rho, x.theta, bet, lambda) {
    -sum(svec*x.theta*(y - n.trials*plogis(rho*x.theta) )) + lambda*rho*sum(abs(bet))
  }
  BetaMat <- NULL
  if(intermed) {
    BetaMat <- matrix(0, nrow=maxiter + 1, ncol=length(beta.init))
    BetaMat[1,] <- beta.init
  }
  if(track.objective) {
    objfn.track = rep(NA, maxiter + 1)
    objfn.track[1] <- LogLikObserved(phi, beta.init)
  }
  rho.old <- 1
  while(TRUE) {
    ## Compute weight vector
    phat <- plogis(phi) # This is expit(phi)
    
    #ww <- svec*ww
    ## How to incorporate the weights?
    theta <- SoftThresh(beta.old + stplngth*(XSty - crossprod(X, svec*phat)), lambda=lambda*stplngth)
    
    iter <- iter + 1
    
    x.theta <- as.numeric(X%*%theta)
    
    rho.pos <- UpdateRhoNewtonL1(1, x.theta, y, n.trials, l1norm=lambda*sum(abs(theta)), 
                                 tol.newton=0.01*tol, svec=svec)
    # rho.neg <- UpdateRhoNewton(1, x.theta, y, n.trials, cc=-lambda*sum(abs(theta)), tol.newton=0.01*tol)
    if(rho.pos > 0) {
      rho.new <- rho.pos
    } else if(rho.pos < 0) {
      rho.new <- 0
    }
    beta.new <- as.numeric(rho.new*theta)
    
    if(norm(beta.new-beta.old, "2") < tol | iter >= maxiter) break
    beta.old <- beta.new
    phi <- as.numeric(X%*%beta.old)
    rho.old <- rho.new
    
    if(track.objective) {
      objfn.track[iter+1] <- LogLikObserved(phi, beta.old)
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
