pxlogistic <- function(beta.init, X, y, n, omega=0, method="Newton", tol=1e-4, max.iter=500){
  beta.old <- theta <- beta.init
  kap <- y - n/2
  iter <- 0
  alpha.old <- 1
  obj <- function(alp) -sum(alp*y*x.theta-log1p(exp(alp*x.theta)))
  while(TRUE){
    phi <- as.numeric(X%*%beta.old)
    ww <- as.numeric((n/(2*phi))*tanh(phi/2))
    theta <- (1+omega)*solve(crossprod(X, X*ww), crossprod(X, kap))/alpha.old - omega*theta
    x.theta <- as.numeric(X%*%theta)
    if(method=="Newton") {
      alpha.new <- UpdateAlphaNewton(alpha.old, x.theta, y, init=!iter)
    } else if(method=="BFGS") {
      it <- ifelse(iter==0, 1e4, 1)
      alpha.new <- optim(alpha.old, obj, method="BFGS", control=list(maxit=it))$par
    }
    beta.new <- alpha.new*theta
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter > max.iter) break
    beta.old <- beta.new
    alpha.old <- alpha.new
  }
  return(list(beta=beta.new, iter=iter))
}