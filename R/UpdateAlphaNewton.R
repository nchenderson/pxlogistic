UpdateRhoNewton <- function(rho.old, x.theta, y, n.trials, init=FALSE,
                              tol.newton, svec) {
  if(init) rho.old <- 0

  x.theta2 <- x.theta*x.theta
  rho.new <- rho.old
  while(TRUE){
    ept <- plogis(rho.old*x.theta)
    dept <- ept*(1-ept)
    f.st <- sum(svec*x.theta*(y - n.trials*ept))
    #if(abs(f.st) < 1e-12) break
    if(abs(f.st) < tol.newton) break
    f.nd <- sum(svec*x.theta2*n.trials*dept)
    rho.new <- as.numeric(rho.old + f.st/f.nd)
    rho.old <- rho.new
  }
  return(rho.new)
}

UpdateRhoNewtonPen <- function(rho.old, x.theta, theta, y, n.trials, init=TRUE,
                            tol.newton, svec, lambda) {
  if(init) rho.old <- 1
  #print(rho.old)
  #print(x.theta)
  #print(theta)
  #print(y)
  #print(n.trials)
  #print(init)
  #print(tol.newton)
  #print(svec)
  #print(lambda)


  x.theta2 <- x.theta*x.theta
  theta2 <- sum(theta^2)
  rho.new <- rho.old
  count <- 0
  while(TRUE){
    ept <- plogis(rho.old*x.theta)
    dept <- ept*(1-ept)
    f.st <- sum(svec*x.theta*(y - n.trials*ept)) - rho.old*lambda*theta2
    #if(abs(f.st) < 1e-12) break
    if(abs(f.st) < tol.newton | count > 10) break
    f.nd <- sum(x.theta2*n.trials*svec*dept) + lambda*theta2
    rho.new <- as.numeric(rho.old + f.st/f.nd)
    rho.old <- rho.new
    count <- count + 1
  }
  return(rho.new)
}

UpdateRhoNewtonL1 <- function(rho.old, x.theta, y, n.trials,  l1norm, init=FALSE,
                            tol.newton, svec) {
  if(init) rho.old <- 0
  x.theta2 <- x.theta*x.theta
  rho.new <- rho.old
  while(TRUE){
    ept <- plogis(rho.old*x.theta)
    dept <- ept*(1-ept)
    f.st <- sum(svec*x.theta*(y - n.trials*ept)) - l1norm
    if(abs(f.st) < tol.newton) break
    f.nd <- sum(x.theta2*svec*n.trials*dept)
    rho.new <- as.numeric(rho.old + f.st/f.nd)
    rho.old <- rho.new
  }
  if(rho.new <= 0) {
     rho.new <- 1
  }
  return(rho.new)
}





