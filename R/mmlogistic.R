mmlogistic <- function(par, X, y, n.trials=rep(1, length(y)), 
                       lambda=0, track.objective=FALSE,
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
  
  beta.old <- par
  
  xtk <- crossprod(X, y - n.trials/2)
  iter <- 0
  B.inv <- solve(crossprod(X)/4 + 2*lambda*diag(ncol(X)))
  while(TRUE) {
    phi <- as.numeric(X%*%beta.old)
    mu <- plogis(phi)
    beta.new <- beta.old + B.inv%*%(crossprod(X, y-mu) - 2*lambda*beta.old)
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter > maxiter) break
    beta.old <- beta.new
  }
  return(list(coef=beta.new, iter=iter, objfn.track=objfn.track))
}
