IRLS_Logistic <- function(par, X, y, n, tol=1e-4, max.iter=500){
  beta.old <- par
  iter <- 0
  linkinv <- binomial()$linkinv
  variance <- binomial()$variance
  while(TRUE){
    phi <- as.numeric(X%*%beta.old)
    mu <- linkinv(phi)
    v <- variance(mu)
    z <- phi + (y/n-mu)/v
    # beta.new <- solve(crossprod(X*v, X)) %*% crossprod(X, v*z)
    beta.new <- qr.solve(X*sqrt(v), sqrt(v)*z) ## this replicates glm
    
    iter <- iter + 1
    if(norm(beta.new-beta.old, "2") < tol | iter > max.iter) break
    beta.old <- beta.new 
  }
  return(list(beta=beta.new, iter=iter))
}