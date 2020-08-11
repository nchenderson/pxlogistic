emlogistic <- function(par, X, y, n.trials=rep(1, length(y)), lambda=NULL, control=list()) {
  
  control.default <- list(maxiter=2000, tol=1e-7)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)
  
  maxiter <- control$maxiter
  tol <- control$tol

  beta.init <- par
  beta.old <- beta.init
  uvec <- y - n.trials/2
  iter <- 0
  
  ww <- rep(0, length(y))
  while(TRUE) {
    ## Compute weight vector
    phi <- as.numeric(X%*%beta.old)
    small.phi <- abs(phi) < 1e-4
    ww[small.phi] <- n.trials[small.phi]*(1/4 - (phi[small.phi]^2)/48 + (phi[small.phi]^4)/480)
    ww[!small.phi] <- as.numeric((n.trials[!small.phi]*tanh(phi[!small.phi]/2))/(2*phi[!small.phi])) 
    
    beta.new <- as.numeric(solve(crossprod(X, X*ww), crossprod(X, uvec)))
  
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter > maxiter) break
    beta.old <- beta.new
  }
  return(list(coef=beta.new, iter=iter))
}
